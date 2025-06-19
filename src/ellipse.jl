using GLMakie
using GLMakie: Gray, N0f8


function conic_to_axes(coeffs, normalise=true)
    # https://www.wikiwand.com/en/Ellipse#/General_ellipse
    # Extract  conic parameters
    if normalise && coeffs[1] != 0
        coeffs = coeffs / coeffs[1]
    end

    A, B, C, D, E, F = coeffs

    # Useful intermediates
    disc = B * B - 4 * A * C

    # Center
    x0 = (2 * C * D - B * E) / disc
    y0 = (2 * A * E - B * D) / disc

    # Axis lengths
    term1 = 2 * (A * E * E + C * D * D - B * D * E + disc * F)
    a = -sqrt(term1 * (A + C + sqrt((A - C)^2 + B^2))) / disc
    b = -sqrt(term1 * (A + C - sqrt((A - C)^2 + B^2))) / disc

    # if A < C
    #     if B == 0
    #         phi_b = 0
    #     else
    #         phi_b = atan(B / (A - C)) / 2
    #     end
    # else
    #     if B == 0
    #         phi_b = π / 2
    #     else
    #         phi_b = atan(B / (A - C)) / 2 - π / 2
    #     end
    # end
    phi_b = atan(-B, C - A) / 2

    return x0, y0, a, b, phi_b
end

function axes_to_conic(centre, axes, angle)
    x0, y0 = centre
    # a, b = sort(axes; rev=true)
    a, b = axes

    A = (a * sin(angle))^2 + (b * cos(angle))^2
    B = 2 * (b * b - a * a) * sin(angle) * cos(angle)
    C = (a * cos(angle))^2 + (b * sin(angle))^2
    D = -2 * A * x0 - B * y0
    E = -B * x0 - 2 * C * y0
    F = A * x0 * x0 + B * x0 * y0 + C * y0 * y0 - a * a * b * b
    return [A, B, C, D, E, F] ./ F
end

axes_to_conic(x0, y0, a, b, phi_b) = axes_to_conic((x0, y0), (a, b), phi_b)

"""
    fit_ellipse(positions::Vector{Point2-like}, weights = ones)

Document this function
"""
function fit_ellipse(positions, w=ones(length(positions)))
    length(positions) < 5 && return zeros(6) # infinite ellipse for less than 5 points
    b = -w
    zzz = point_conic_eq.(positions) .* w
    A = hcat(zzz...)'
    coeffs = A \ b
    return vcat(coeffs, [1])
end  # function fit_ellipse

# function fit_ellipse(positions)
#     length(positions) < 5 && return zeros(6) # infinite ellipse for less than 5 points
#     b = -ones(length(to_value(positions)))
#     zzz = point_conic_eq.(to_value(positions))
#     A = hcat(zzz...)'
#     coeffs = A \ b
#     return vcat(coeffs, [1])
# end

function point_conic_eq(x, y)
    return [x^2, x * y, y^2, x, y]
end

point_conic_eq(p::Point2) = point_conic_eq(p.data...)
point_conic_eq(p::CartesianIndex) = point_conic_eq(p.I...)


function getellipsepoints(cx, cy, rx, ry, θ)
    t = range(0, 2 * pi; length=100)
    ellipse_x_r = @. rx * cos(t)
    ellipse_y_r = @. ry * sin(t)
    R = [cos(θ) sin(θ); -sin(θ) cos(θ)]
    r_ellipse = [ellipse_x_r ellipse_y_r] * R
    x = @. cx + r_ellipse[:, 1]
    y = @. cy + r_ellipse[:, 2]
    return (x, y)
end

function _pos_to_elpoints(positions)
    if length(positions) > 4
        coeffs = fit_ellipse(to_value(positions))
        el = conic_to_axes(coeffs)
        ppp = getellipsepoints(el...)
        return elp = map(Point2f, zip(ppp...))
    else
        return Point2f[]
    end
end


function inside_ellipse(coeffs, position; shift=false)
    x, y = Tuple(position)
    if shift
        x -= 0.5
        y -= 0.5
    end

    return [x^2, x * y, y^2, x, y, 1]' * coeffs <= 0
end

function mask_ellipse(img, coeffs)
    mask = similar(img, Bool)
    for i in CartesianIndices(img)
        mask[i] = inside_ellipse(coeffs, i; shift=true)
    end
    return mask
end

"""
    draw_ellipse(img)

Display an image and let you draw an ellipse there.
    Controls :
    a+ click to add a point
    d+ click to delete a point
    mouse wheel or left mouse press and draw to zoom
    right mouse to pan
    Ctrl+click to reset zoom
    Esc or q to quit
"""
function draw_ellipse(img)
    ready = false
    displayhelp = Observable(true)
    state = Observable(:pass)

    helpmessage = """
    Controls :
    a to switch to "add a point" mode
    d to switch to "delete a point" mode
    p to switch to passive mode
    e to toggle ellipse/polygon mode
    h to show/hide this message
    mouse wheel or left mouse press and draw to zoom
    right mouse to pan
    Ctrl+click to reset zoom
    Esc or q to quit
    """

    isbig = size(img)[1] > 128 && size(img)[2] > 128

    fig, ax, implot = image(rotr90(img); axis=(aspect=DataAspect(),), interpolate=isbig)

    mode = :ellipse
    positions = Observable(Point2f[])
    current_el = lift(fit_ellipse, positions)
    el_points = map(_pos_to_elpoints, positions)
    overlay = lift(current_el) do current_el
        mask = mask_ellipse(rotr90(img), current_el)
        ap2mask(1 .- mask)
    end

    image!(overlay; colormap=(:Blues, 0.4))

    p = scatter!(ax, positions)
    c = lines!(ax, el_points)

    text!(
        0.1,
        0.9;
        text=helpmessage,
        space=:relative,
        color=:white,
        font=:bold,
        align=(:left, :top),
        visible=displayhelp,
        glowwidth=1,
        fontsize=24,
    )

    on(events(fig).mousebutton; priority=2) do event
        if event.button == Mouse.left && event.action == Mouse.press
            if state[] == :delete
                # Delete marker
                plt, i = pick(fig)
                if plt == p
                    deleteat!(positions[], i)
                    notify(positions)
                    return Consume(true)
                end
            elseif state[] == :add
                # Add marker
                push!(positions[], mouseposition(ax))
                notify(positions)
                return Consume(true)
            end
        end
        return Consume(false)
    end

    screen = display(fig)

    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat
            if event.key == Keyboard.q || event.key == Keyboard.escape
                println("q/esc is pressed, the data are saved")
                println("$(fit_ellipse(to_value(positions)))")
                ready = true
                close(screen)
                # return el .= current_el[]
            elseif event.key == Keyboard.h
                displayhelp[] = !displayhelp[]
            elseif event.key == Keyboard.a
                state[] = :add
            elseif event.key == Keyboard.d
                state[] = :delete
            elseif event.key == Keyboard.p
                state[] = :pass
            end

        end
    end
    wait(screen)
    return to_value(current_el), rotr90(mask_ellipse(rotr90(img), to_value(current_el)), 3)
end  # function draw_ellipse

"""
    draw_or_load_ellipse(elfile, img, apfile="", saveap=true)

    Loads an ellipse from jld2 file named `elfile`, and if not found, it will display the image `img` and let you draw an ellipse there.
    If saveap is true, it will save the aperture to `apfile` ("ap.png" by default) in the same directory as `elfile`.


Display an image and let you draw an ellipse there.
    Controls :
    a+ click to add a point
    d+ click to delete a point
    mouse wheel or left mouse press and draw to zoom
    right mouse to pan
    Ctrl+click to reset zoom
    Esc or q to quit
"""
function draw_or_load_ellipse(elfile, img, apfile=""; saveap=false)
    if apfile == ""
        apfile = joinpath(dirname(elfile), "ap.tif")
    end

    return try
        el = load(elfile, "el")
        @info "Aperture loaded"
        el
    catch
        @info "Aperture not found, please draw it"
        GLMakie.activate!()
        el, ap = draw_ellipse(img)
        saveap && (save(apfile, Gray{N0f8}.(ap));
        @info "$apfile is saved")
        jldsave(elfile; el)
        CairoMakie.activate!(; type="png")
        el
    end
end



# Updated draw_aperture with an additional convex hull mode.
function draw_aperture(img)
    ready = false
    displayhelp = Observable(true)
    state = Observable(:pass)

    helpmessage = """
    Controls :
    a to switch to "add a point" mode
    d to switch to "delete a point" mode
    p to switch to passive mode
    e to toggle through ellipse, polygon, and convex hull modes
    h to show/hide this message
    mouse wheel or left mouse press & drag: Zoom
    Right mouse: Pan
    Ctrl+click: Reset zoom
    Esc or q: Quit and save the data
    """

    rotated = rotr90(img)
    # sim = size(img)
    xrange, yrange = Base.axes(img)
    fig, ax, implot = image(rotated; axis=(aspect=DataAspect(),))

    # Observable mode can be :ellipse, :polygon or :hull.
    mode = Observable(:ellipse)
    positions = Observable(Point2f[])

    current_el = lift(fit_ellipse, positions)
    # Compute mask for user-drawn polygon.
    current_poly = lift(positions) do pos
        if length(pos) ≥ 3
            return draw_filled_polygon_on_ranges(xrange, yrange, reverse.(pos); fill_value=1)
        else
            return zeros(Float32, size(rotated)...)
        end
    end
    # Compute mask for the convex hull of the positions.
    current_hull = lift(positions) do pos
        if length(pos) ≥ 3
            hull_pts = convex_hull(reverse.(pos))
            return draw_filled_polygon_on_ranges(xrange, yrange, hull_pts; fill_value=1)
        else
            return zeros(Float32, size(rotated)...)
        end
    end

    # The overlay is built depending on the selected mode.
    overlay = lift(current_el, current_poly, current_hull, mode) do el, poly, hull, m

        if m == :ellipse && any(el .!= 0)
            mask = mask_ellipse(rotated, el)
            return ap2mask(1 .- mask)
        elseif m == :polygon
            return ap2mask(1 .- poly)
        elseif m == :hull
            return ap2mask(1 .- hull)
        else
            return zeros(eltype(rotated), size(rotated)...)
        end
    end

    image!(overlay; colormap=(:Blues, 0.4))

    p = scatter!(ax, positions)
    c = lines!(ax, map(_pos_to_elpoints, positions))

    text!(
        0.1,
        0.9;
        text=helpmessage,
        space=:relative,
        color=:white,
        font=:bold,
        align=(:left, :top),
        visible=displayhelp,
        glowwidth=1,
        fontsize=24,
    )

    on(events(fig).mousebutton; priority=2) do event
        if event.button == Mouse.left && event.action == Mouse.press
            if state[] == :delete
                plt, i = pick(fig)
                if plt == p
                    deleteat!(positions[], i)
                    notify(positions)
                    return Consume(true)
                end
            elseif state[] == :add
                push!(positions[], mouseposition(ax))
                notify(positions)
                return Consume(true)
            end
        end
        return Consume(false)
    end

    screen = display(fig)

    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat)
            if event.key == Keyboard.q || event.key == Keyboard.escape
                println("q/esc is pressed, the data are saved")
                println("$(fit_ellipse(to_value(positions)))")
                ready = true
                close(screen)
            elseif event.key == Keyboard.h
                displayhelp[] = !displayhelp[]
            elseif event.key == Keyboard.a
                state[] = :add
            elseif event.key == Keyboard.d
                state[] = :delete
            elseif event.key == Keyboard.p
                state[] = :pass
            elseif event.key == Keyboard.e
                # Cycle the mode: ellipse -> polygon -> hull -> ellipse -> ...
                mode[] = (
                    mode[] == :ellipse ? :polygon : (mode[] == :polygon ? :hull : :ellipse)
                )
            end
        end
    end

    wait(screen)
    # Return different values based on the mode.
    if mode[] == :ellipse
        return to_value(current_el), rotr90(mask_ellipse(rotated, to_value(current_el)), 3)
    elseif mode[] == :polygon
        return reverse.(to_value(positions)), rotr90(to_value(current_poly), 3)
    else  # mode == :hull
        let pts = reverse.(to_value(positions))
            hull_pts = length(pts) ≥ 3 ? convex_hull(pts) : pts
            return hull_pts, rotr90(to_value(current_hull), 3)
        end
    end
end  # function draw_aperture

####
# Approach through a special type and Makie recipes
#

struct Ellipse{T<:Real}
    A::T
    B::T
    C::T
    D::T
    E::T
    F::T
end

Ellipse(A, B, C, D, E, F) = Ellipse(promote(A, B, C, D, E, F)...)

Ellipse(coef::Vector) = Ellipse(coef...)

conic(x::Ellipse) = [x.A, x.B, x.C, x.D, x.E, x.F]

function centeraxesangle(el::Ellipse)
    cas = conic_to_axes(conic(el))
    return (center=cas[1:2], axes=cas[3:4], ϕb=cas[5])
end

mask_ellipse(img, el::Ellipse) = mask_ellipse(img, conic(el))


function scale(el::Ellipse, c)
    center, ax, angle = centeraxesangle(el)
    return Ellipse(axes_to_conic(center, ax .* c, angle))
end

scale_el(el::Vector, c) = conic(scale(Ellipse(el), c))


using CairoMakie
CairoMakie.convert_arguments(::Type{<:AbstractPlot}, x::Ellipse) =
    (map(Point2f, zip(getellipsepoints(conic_to_axes(conic(x))...)...)),)




"""
    point_inside_polygon(pos::Point2f, polygon::AbstractVector{<:Point2f})

Return `true` if the point `pos` is inside the polygon defined by the list
of vertices in `polygon` (using the ray-casting algorithm), and `false` otherwise.
"""
function point_inside_polygon(pos::Point2f, polygon::AbstractVector{<:Point2f})
    x, y = pos
    n = length(polygon)
    inside = false
    j = n
    for i in 1:n
        xi, yi = polygon[i].x, polygon[i].y
        xj, yj = polygon[j].x, polygon[j].y
        if (yi > y) != (yj > y) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end


"""
    draw_filled_polygon_on_ranges(x_range::AbstractRange, y_range::AbstractRange,
                                  polygon::AbstractVector{<:Point2f};
                                  fill_value = one(Float32))

Build a target array using the coordinate arrays defined by `x_range` (columns) and `y_range` (rows),
draw a filled polygon defined by `polygon` (whose vertices are in the same coordinate system) using
a scanline algorithm, and return the updated array.
"""
function draw_filled_polygon_on_ranges(
    x_range::AbstractRange,
    y_range::AbstractRange,
    polygon::AbstractVector{<:Point2f};
    fill_value=one(Float32),
)
    xvec = collect(x_range)
    yvec = collect(y_range)
    target = zeros(typeof(fill_value), length(yvec), length(xvec))
    draw_filled_polygon!(target, polygon, xvec, yvec; fill_value=fill_value)
    return target
end

"""
    draw_filled_polygon!(target::AbstractArray{T,2}, polygon::AbstractVector{<:Point2f};
                         fill_value::T = one(T))

Draw a filled polygon defined by the vertices in `polygon` onto the 2D array `target`.
The polygon coordinates are assumed to belong to the coordinate system represented by the pixel
indices of `target` (columns are 1:size(target,2) and rows are 1:size(target,1)). The algorithm
uses a scanline approach along vertical lines (at each x value) to fill in the pixels that fall within
the polygon, setting them to `fill_value`. The array is modified in-place.
"""
function draw_filled_polygon!(
    target::AbstractArray{T,2}, polygon::AbstractVector{<:Point2f}; fill_value::T=one(T)
) where {T}
    nrows, ncols = size(target)
    xvec = collect(1:ncols)
    yvec = collect(1:nrows)
    return draw_filled_polygon!(target, polygon, xvec, yvec; fill_value=fill_value)
end

"""
    draw_filled_polygon!(target::AbstractArray{T,2}, polygon::AbstractVector{<:Point2f},
                         xvec::AbstractVector, yvec::AbstractVector; fill_value::T = one(T))

Draw a filled polygon defined by the vertices in `polygon` onto the 2D array `target` (with dimensions
determined by `yvec` and `xvec`). The polygon coordinates are assumed to belong to the coordinate system
represented by `xvec` (columns) and `yvec` (rows). The algorithm uses a scanline approach along vertical
lines (at each x value) to fill in the pixels that fall within the polygon, setting them to `fill_value`.
The array is modified in-place.
"""
function draw_filled_polygon!(
    target::AbstractArray{T,2},
    polygon::AbstractVector{<:Point2f},
    xvec::AbstractVector,
    yvec::AbstractVector;
    fill_value::T=one(T),
) where {T}
    nrows, ncols = size(target)
    np = length(polygon)
    # For each column, use the corresponding x coordinate.
    for c in 1:ncols
        x = xvec[c]
        intersections = Float64[]
        j = np
        for i in 1:np
            x1, y1 = polygon[j].data[1], polygon[j].data[2]
            x2, y2 = polygon[i].data[1], polygon[i].data[2]
            # Check if the vertical line at x crosses the edge between (x1,y1) and (x2,y2)
            if (x1 ≤ x && x2 > x) || (x2 ≤ x && x1 > x)
                intersect_y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
                push!(intersections, intersect_y)
            end
            j = i
        end
        sort!(intersections)
        # Fill the pixels between pairs of intersections.
        k = 1
        while k < length(intersections)
            y_start, y_end = intersections[k], intersections[k + 1]
            if y_start > y_end
                y_start, y_end = y_end, y_start
            end
            lo = searchsortedfirst(yvec, y_start)
            hi = searchsortedlast(yvec, y_end)
            lo = max(lo, 1)
            hi = min(hi, nrows)
            for r in lo:hi
                target[r, c] = fill_value
            end
            k += 2
        end
    end
    return target
end

# Helper function: compute convex hull of points using Andrew’s monotone chain algorithm.
function convex_hull(pts::Vector{Point2f})
    if length(pts) < 3
        return pts
    end
    sorted_pts = sort(pts; by=p -> (p[1], p[2]))
    lower = Point2f[]
    for p in sorted_pts
        while length(lower) ≥ 2 &&
            (
                (lower[end][1] - lower[end - 1][1]) * (p[2] - lower[end - 1][2]) -
                (lower[end][2] - lower[end - 1][2]) * (p[1] - lower[end - 1][1])
            ) ≤ 0
            pop!(lower)
        end
        push!(lower, p)
    end
    upper = Point2f[]
    for p in reverse(sorted_pts)
        while length(upper) ≥ 2 &&
            (
                (upper[end][1] - upper[end - 1][1]) * (p[2] - upper[end - 1][2]) -
                (upper[end][2] - upper[end - 1][2]) * (p[1] - upper[end - 1][1])
            ) ≤ 0
            pop!(upper)
        end
        push!(upper, p)
    end
    pop!(lower)
    pop!(upper)
    return vcat(lower, upper)
end

"""
    transform_marker(p::Point2f/Vector, img) -> Point2f/Vector

Convert a marker coordinate `p`, specified in the original image space
(with cartesian coordinates as (row, col)), into the rotated array coordinate system.
Assuming the image is rotated by 90° counterclockwise, the new column becomes the original row,
and the new row is nrows - original column + 1.
This complies with Julia’s column-major order.
Adjust if your rotation convention (or display coordinate system) differs.
"""
function transform_marker(p, img)
    nrows, _ = size(img)   # nrows from original image.
    # Here we assume that p[1] is the original row and p[2] is the original column.
    # After rotr90 (90° counterclockwise), new row = nrows - original column + 1,
    # and new column = original row.
    return typeof(p)(nrows - p[2] + 1, p[1])
end

# When drawing filled polygons, we use ranges corresponding to the rotated image.
# (Rows: 1:nrows, Columns: 1:ncols)

# Refactored draw_aperture function using the improved transform_marker.
function draw_aperture2(img)
    ready = false
    displayhelp = Observable(true)
    state = Observable(:pass)

    helpmessage = """
    Controls :
      a     to switch to "add a point" mode
      d     to switch to "delete a point" mode
      p     to switch to passive mode
      e     to toggle through ellipse, polygon, and convex hull modes
      h     to show/hide this message
      Mouse wheel or left mouse press & drag: Zoom
      Right mouse: Pan
      Ctrl+click: Reset zoom
      Esc or q: Quit and save the data
    """

    # Rotate the image once.
    rotated = rotr90(img)
    # Note: rotated now has size (nrows_rot, ncols_rot) and is displayed as-is.
    fig, ax, implot = image(rotated; axis=(aspect=DataAspect(),))

    # mode will be one of :ellipse, :polygon, or :hull.
    mode = Observable(:ellipse)
    # Positions are stored in the coordinate system for "rotated" image: (row, col)
    positions = Observable(Point2f[])

    # When adding a marker, immediately convert the raw marker from the displayed coordinate system
    # (traditional image space) into the rotated (matrix index) space.
    on(events(fig).mousebutton; priority=2) do event
        if event.button == Mouse.left && event.action == Mouse.press
            if state[] == :delete
                plt, i = pick(fig)
                if plt == p
                    deleteat!(positions[], i)
                    notify(positions)
                    return Consume(true)
                end
            elseif state[] == :add
                raw_pos = mouseposition(ax)
                mpos = transform_marker(raw_pos, img)
                push!(positions[], mpos)
                notify(positions)
                return Consume(true)
            end
        end
        return Consume(false)
    end

    # Build the ellipse, polygon and convex hull masks based on the rotated image coordinates.
    current_el = lift(fit_ellipse, positions)
    current_poly = lift(positions) do pos
        if length(pos) ≥ 3
            nrows, ncols = size(rotated)
            return draw_filled_polygon_on_ranges(1:ncols, 1:nrows, pos; fill_value=1)
        else
            return zeros(Float32, size(rotated)...)
        end
    end
    current_hull = lift(positions) do pos
        if length(pos) ≥ 3
            hull_pts = convex_hull(pos)
            nrows, ncols = size(rotated)
            return draw_filled_polygon_on_ranges(1:ncols, 1:nrows, hull_pts; fill_value=1)
        else
            return zeros(Float32, size(rotated)...)
        end
    end

    # Build overlay according to the current mode.
    overlay = lift(current_el, current_poly, current_hull, mode) do el, poly, hull, m
        if m == :ellipse && any(el .!= 0)
            mask = mask_ellipse(rotated, el)
            return ap2mask(1 .- mask)
        elseif m == :polygon
            return ap2mask(1 .- poly)
        elseif m == :hull
            return ap2mask(1 .- hull)
        else
            return zeros(eltype(rotated), size(rotated)...)
        end
    end

    image!(overlay; colormap=(:Blues, 0.4))
    p = scatter!(ax, positions)
    c = lines!(ax, map(_pos_to_elpoints, positions))

    text!(
        0.1,
        0.9;
        text=helpmessage,
        space=:relative,
        color=:white,
        font=:bold,
        align=(:left, :top),
        visible=displayhelp,
        glowwidth=1,
        fontsize=24,
    )

    on(events(fig).keyboardbutton) do event
        if event.action in (Keyboard.press, Keyboard.repeat)
            if event.key == Keyboard.q || event.key == Keyboard.escape
                println("q/esc is pressed, the data are saved")
                println("$(fit_ellipse(to_value(positions)))")
                ready = true
                close(screen)
            elseif event.key == Keyboard.h
                displayhelp[] = !displayhelp[]
            elseif event.key == Keyboard.a
                state[] = :add
            elseif event.key == Keyboard.d
                state[] = :delete
            elseif event.key == Keyboard.p
                state[] = :pass
            elseif event.key == Keyboard.e
                # Cycle through modes.
                mode[] = (
                    mode[] == :ellipse ? :polygon : (mode[] == :polygon ? :hull : :ellipse)
                )
            end
        end
    end

    screen = display(fig)
    wait(screen)
    # Rotate the mask back to the original orientation (using rotr90 with appropriate count)
    if mode[] == :ellipse
        return to_value(current_el), rotr90(mask_ellipse(rotated, to_value(current_el)), 3)
    elseif mode[] == :polygon
        return to_value(positions), rotr90(to_value(current_poly), 3)
    else  # mode == :hull
        let pts = to_value(positions)
            hull_pts = length(pts) ≥ 3 ? convex_hull(pts) : pts
            return hull_pts, rotr90(to_value(current_hull), 3)
        end
    end
end  # function draw_aperture2
