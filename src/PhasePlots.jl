module PhasePlots

using CairoMakie
using FFTViews
using PhaseUtils

include("fileutils.jl")
include("ellipse.jl")
include("phasedisplay.jl")



export dpng, save_figs_as_gif

end
