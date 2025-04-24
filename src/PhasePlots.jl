module PhasePlots

using CairoMakie
using FFTViews
using PhaseUtils
using JLD2

include("fileutils.jl")
include("ellipse.jl")
include("phasedisplay.jl")



export dpng, save_figs_as_gif, draw_ellipse, draw_or_load_ellipse

end
