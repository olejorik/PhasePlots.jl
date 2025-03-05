phasemap = :cyclic_mygbm_30_95_c78_n256

# function showarray(arr, colormap=:viridis)
#     return heatmap(rotr90(arr); colormap=colormap, axis=(aspect=DataAspect(),))
# end

function showarray!(arr; colormap=:viridis, rot=1, args...)
    return heatmap!(rotr90(arr, rot); colormap=colormap, args...)
end

function showarray!(ax, arr; colormap=:viridis, rot=1, args...)
    return heatmap!(ax, rotr90(arr, rot); colormap=colormap, args...)
end

function showarray(arr; colormap=:viridis, rot=1, args...)
    return heatmap(
        rotr90(arr, rot);
        colormap=colormap,
        args...,
        axis=merge(get(args, :axis, (;)), (aspect=DataAspect(),)),
    )
end

function showphase(inarr; rot=1, fig=Figure(), picsize=512, cm=phasemap)
    # if max(size(rotr90(inarr))...) > picsize
    #     arr = imresize(inarr, picsize)
    # else
    arr = Array(inarr)
    # end

    ax = CairoMakie.Axis(fig[1, 1]; aspect=1)
    hm = heatmap!(ax, phwrap.(rotr90(arr, rot)); colormap=cm, colorrange=(-π, π))
    cb = Colorbar(fig[1, 2], hm; width=10, tellheight=true)
    return fig, ax, cb
end

function showphase!(ax, inarr; rot=1, picsize=512, cm=phasemap)
    arr = Array(inarr)
    hm = heatmap!(ax, phwrap.(rotr90(arr, rot)); colormap=cm, colorrange=(-π, π))
    return hm
end

function showphasetight(
    inarr, fig=Figure(); picsize=512, cm=phasemap, hidedec=true, kwarg...
)
    inarr = bboxview(Array(inarr))
    # if max(size(inarr)...) > picsize
    #     arr = imresize(inarr, picsize)
    # else
    arr = inarr
    # end

    if typeof(fig) == GridPosition
        pos = fig
    else
        pos = fig[1, 1]
    end
    ax = CairoMakie.Axis(pos; aspect=AxisAspect(1))
    hm = heatmap!(ax, phwrap.(rotr90(arr)); colormap=cm, colorrange=(-π, π), kwarg...)
    if hidedec
        hidedecorations!(ax; grid=false)
    end
    return fig, ax, hm
end

@recipe(PhasePlot, arr) do scene
    Attributes(; colormap=phasemap, colorrange=(-π, π), crop=true)
    # Theme(
    #         Axis = (
    #             aspect = 1,
    #             leftspinevisible = false,
    #             rightspinevisible = false,
    #             bottomspinevisible = false,
    #             topspinevisible = false,
    #             yticksvisible = false,
    #             xticksvisible = false,
    #             yticklabelsvisible = false,
    #             xticklabelsvisible = false,
    #             xautolimitmargin = (0, 0),
    #             yautolimitmargin = (0, 0),
    #         ),
    # )
end

function Makie.plot!(p::PhasePlot{<:Tuple{<:AbstractArray}})
    arr = p[:arr][]
    # arr =p[:arr][] |> rotr90
    if p[:crop][]
        @info "cropped array"
        arr = bboxview(p[:arr][])
    end
    # hm = heatmap!(p, rotr90(arr); colormap = p[:colormap],
    # colorrange=p[:colorrange][],
    # axis = (aspect =  1,)
    # )
    with_theme(phasetheme) do
        heatmap!(p, rotr90(arr))
    end

    # tightlimits!(p.plots.axis)
    # @info p.plots
    return p
end

phasetheme = Theme(;
    Axis=(
        aspect=1,
        leftspinevisible=false,
        rightspinevisible=false,
        bottomspinevisible=false,
        topspinevisible=false,
        yticksvisible=false,
        xticksvisible=false,
        yticklabelsvisible=false,
        xticklabelsvisible=false,
        xautolimitmargin=(0, 0),
        yautolimitmargin=(0, 0),
    ),
    colormap=phasemap,
)

# functions to show array of similar plots
#
#
"""
    plot_heatmaps_table(
    heatmaps_array;
    ncols::Int=0,
    width=150,
    height=150,
    colormap=:viridis,
    limits = (0,0),
    hidedecorations=false,
    rot=1, # rotate each array `rot` times 90° CCW
    kwargs...,
)

Show a vector of 2D arrays as a matrix of heatmaps with a common colorbar below.
"""
function plot_heatmaps_table(
    heatmaps_array;
    ncols::Int=0,
    width=150,
    height=150,
    colormap=:viridis,
    limits=(0, 0),
    hidedecorations=false,
    rot=1,
    titles="",
    title="",
    titlesize=20,
    kwargs...,
)

    if titles == ""
        titles = fill("", length(heatmaps_array))
    end

    if ncols == 0
        ## Calculate the number of columns based on the number of heatmaps
        ncols = ceil(Int, sqrt(length(heatmaps_array)))
    end
    ind(i) = divrem(i - 1 + ncols, ncols)

    _nrows = ceil(Int, length(heatmaps_array) / ncols)
    fig = Figure(; size=(width * ncols, height * _nrows))

    ## Define a common color range for all heatmaps
    if limits == (0, 0)
        min_val = minimum([minimum(filter(!isnan, hm)) for hm in heatmaps_array])
        max_val = maximum([maximum(filter(!isnan, hm)) for hm in heatmaps_array])
    else
        min_val, max_val = limits
    end

    ## Generate heatmaps
    for (i, hm) in enumerate(heatmaps_array)
        ax = Axis(fig[ind(i)...]; width=width, height=height, aspect=DataAspect())
        ax.title = titles[i]
        if hidedecorations
            hidedecorations!(ax)
        end
        heatmap!(
            ax, rotr90(hm, rot); colorrange=(min_val, max_val), colormap=colormap, kwargs...
        )
    end

    ## Add a common colorbar
    Colorbar(fig[end + 1, :]; limits=(min_val, max_val), colormap=colormap, vertical=false)

    if title != ""
        Label(fig[0, :], title; fontsize=titlesize)
    end


    resize_to_layout!(fig)
    return fig
end

# # to mark some details on the plot
circle_with_hole(r=0.9) = BezierPath([
    MoveTo(Point(1, 0)),
    EllipticalArc(Point(0, 0), 1, 1, 0, 0, 2pi),
    EllipticalArc(Point(0, 0), r, r, 0, 0, -2pi),
    ClosePath(),
])

export showphase, showphasetight, showarray, phasemap, showarray!, phaseplot, phaseplot!
export phasetheme
export plot_heatmaps_table, circle_with_hole
