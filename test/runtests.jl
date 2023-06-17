using PhasePlots
using Test
using PhaseUtils
using CairoMakie

@testset "PhasePlots.jl" begin
    # Write your tests here.
    s = (110,100)
    arr = zeros(s)
    y = range(-1.21,1.21,s[1])
    x = range(-1.1,1.1,s[2])
    
    arr[[x.^2 + y.^2 .<= 1 for y in y, x in x]] .= 1 
    showarray(arr)
    phaseplot(arr, axis = (aspect = DataAspect(),))

    mask = ap2mask(arr)
    phaseplot(arr .* mask,axis = (aspect = DataAspect(),))
    
    phase = copy(arr)
    phase = [-x^3 + 3x.^2 + y.^2 - 10y for y in y, x in x]
    phaseplot(phase .* mask,axis = (aspect = DataAspect(),))
    fig,ax, hm = phaseplot(phwrap(phase .* mask),axis = (aspect = DataAspect(),))
    ax.title = "Using `crop = true` option"
    hidedecorations!(ax, grid=false)
    fig |> display
    fig, ax, hm = phaseplot(phwrap(phase .* mask), crop = false,axis = (aspect = DataAspect(),))
    ax.title = "Using `crop = false` option"
    fig |> display
    fig, ax,hm = showphasetight(phase .* mask, hidedec = false)
    ax.title = "Using function `showphasetight`"
    ax.subtitle = "`hidedec` = false"
    fig |> display
    fig, ax,hm = showphasetight(phase .* mask, hidedec = true)
    ax.title = "Using function `showphasetight`"
    ax.subtitle = "`hidedec` = true"
    fig |> display
    
    

    fig, ax,hm = with_theme(phasetheme) do 
       heatmap(bboxview(phwrap(rotr90(phase .* mask)) )) 
    end
    ax.title = "Using phasetheme"
    fig |> display

    fig, ax,hm = with_theme(phasetheme) do 
        phaseplot(phwrap(phase .* mask)) 
     end
     ax.title = "Using phasetheme and `phaseplot`"
     fig |> display
    

end
