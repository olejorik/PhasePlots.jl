phasemap = :cyclic_mygbm_30_95_c78_n256

# function showarray(arr, colormap=:viridis)
#     return heatmap(rotr90(arr); colormap=colormap, axis=(aspect=DataAspect(),))
# end

function showarray!(arr, colormap=:viridis)
    return heatmap!(rotr90(arr); colormap=colormap,)
end

function showarray(arr, colormap = :viridis; args...)
    return heatmap(rotr90(arr), colormap = colormap; axis=(
        aspect=DataAspect(),
    ), args...)
end

function showphase(inarr, fig=Figure(), picsize=512, cm=:cyclic_mygbm_30_95_c78_n256)
    if max(size(rotr90(inarr))...) > picsize
        arr = imresize(inarr, picsize)
    else
        arr = inarr
    end

    ax = CairoMakie.Axis(
        fig[1, 1];
        aspect=1,
        )
    hm = heatmap!(ax, phwrap.(rotr90(arr)); colormap=cm)
    cb = Colorbar(fig[1, 2], hm; width=10, tellheight=true)
    return fig, ax, cb
end

function showphasetight(
    inarr,
    fig=Figure();
    picsize=512,
    cm=:cyclic_mygbm_30_95_c78_n256,
    hidedec=true,
    kwarg...,
)
    inarr = bboxview(rotr90(inarr))
    if max(size(inarr)...) > picsize
        arr = imresize(inarr, picsize)
    else
        arr = inarr
    end
    rows, cols = size(arr)

    if typeof(fig) == GridPosition
        pos = fig
    else
        pos = fig[1, 1]
    end
    ax = CairoMakie.Axis(
        pos;
        aspect=AxisAspect(1),
        #    autolimitaspect = 1,
        #    xlabel = L"\sigma_x",
        #    ylabel = L"\sigma_y",
        #    xticks = ([0.5, cols / 2 + 0.5, cols + 0.5], [L"-\frac{1}{2s}", L"0", L"\frac{1}{2s}"]),
        #    yticks = ([0.5, rows / 2 + 0.5, rows + 0.5], [L"-\frac{1}{2s}", L"0", L"\frac{1}{2s}"])
    )
    # contourrange = map(x-> 0.5:1:x, size(arr))
    hm = heatmap!(ax, phwrap.(rotr90(arr)); colormap=cm, colorrange=(-π, π), kwarg...)
    # hm = heatmap!(ax, contourrange[2], contourrange[1], phwrap.(arr), colormap = cm, colorrange =(-π,π), kwarg...)
    if hidedec
        hidedecorations!(ax; grid=false)
    end
    # cb = Colorbar(fig[1,2], hm, width = 10, tellheight=true)
    # cb.ticks = (-π:π/2:π, ["-π","-π/2", "0","π/2","π"])
    return fig, ax, hm
end

export showphase, showphasetight, showarray, phasemap, showarray!
