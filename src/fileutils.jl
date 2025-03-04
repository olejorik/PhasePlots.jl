
"""
    dpng(x)

Displays `x` as PNG. Can be used when the default representation
(SVG or PDF) is too heavy.
"""
dpng(x) = display("image/png", x)

struct save2
    name::String
    subfolder::String
end

save2(name) = save2(name, "notebooks/phase_from_interferograms/images/")


# (s::save2)(fig) = (wsave("$(s.subfolder)$(s.name).pdf", fig); wsave("$(s.subfolder)$(s.name).png", fig);)




function save_figs_as_gif(figs, name, fps=2)
    if @isdefined Makie
        cb = Makie.current_backend()
    else
        error("No Makie found")
    end
    GLMakie.activate!()
    figanim = Figure(; size=figs[1].scene.camera.resolution[])
    record(figanim, name, figs; framerate=fps) do f
        figanim = f
        display(figanim)
    end
    return cb.activate!()
end
