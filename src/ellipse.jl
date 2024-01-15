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

    if A < C
        if B == 0
            phi_b = 0
        else
            phi_b = atan(B / (A - C)) / 2
        end
    else
        if B == 0
            phi_b = π / 2
        else
            phi_b = atan(B / (A - C)) / 2 - π / 2
        end
    end

    return x0, y0, a, b, phi_b
end

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


function inside_ellipse(coeffs, position)
    x, y = Tuple(position)
    return [x^2, x * y, y^2, x, y, 1]' * coeffs <= 0
end

function mask_ellipse(img, coeffs)
    mask = similar(img, Bool)
    for i in CartesianIndices(img)
        mask[i] = inside_ellipse(coeffs, i)
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
    a to swith to "add a point" mode
    d to swith to "delete a point" mode
    p to swith to passive mode
    h to show/hide this message
    mouse wheel or left mouse press and draw to zoom
    right mouse to pan
    Ctrl+click to reset zoom
    Esc or q to quit
    """

    fig, ax, implot = image(rotr90(img); axis=(aspect=DataAspect(),))

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
        fontssize=24,
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


export draw_ellipse

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

conic(x::Ellipse) = [x.A, x.B, x.C, x.D, x.E, x.F]

function centeraxesangle(el::Ellipse)
    cas = conic_to_axes(conic(el))
    return (center=cas[1:2], axes=cas[3:4], ϕb=cas[5])
end

mask_ellipse(img, el::Ellipse) = mask_ellipse(img, conic(el))

using CairoMakie
CairoMakie.convert_arguments(::Type{<:AbstractPlot}, x::Ellipse) =
    (map(Point2f, zip(getellipsepoints(conic_to_axes(conic(x))...)...)),)
