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

export showphase, showphasetight, showarray, phasemap, showarray!, phaseplot, phaseplot!
export phasetheme
