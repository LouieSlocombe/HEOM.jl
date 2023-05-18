function os_display(p)
    """
    OS safe plot display, prevents issues with plotting when on a cluster
    and not on an X11 forwarding
    """
    if Sys.iswindows()
        display(p)
    end
end

function plot_wigner_heatmap(x, y, z; name="w_heatmap.pdf", title="")
    xlab = latexstring("Q \\: \\left[\\alpha_0 \\right]")
    ylab = latexstring("P \\: \\left[\\hbar / \\alpha_0 \\right]")
    z = permutedims(z, (2, 1))
    p = heatmap(
        x,
        y,
        z,
        seriescolor=:inferno,
        xlabel=xlab,
        ylabel=ylab,
        title=title,
        xlims=(minimum(x), maximum(x)),
        ylims=(minimum(y), maximum(y)),
        right_margin=10.0Plots.mm,
    )
    savefig(p, joinpath(plot_dump, name))
    os_display(p)
    return p
end

function animate_wigner_heatmap(
    q,
    p,
    sol;
    time=nothing,
    name="ani_wigner_heatmap.gif",
    title=""
)
    """
    Animate wigner function heatmap
    """
    println("Animating wigner heatmap...")
    if time === nothing
        nt = length(sol.t)
        z_vals = real.(sol.u)
    else
        nt = length(time)
        z_vals = real.(sol)
    end
    # Get the min max range to fix the plot
    z_max = maximum(z_vals)
    z_min = minimum(z_vals)

    xlab = latexstring("Q \\: \\left[\\alpha_0 \\right]")
    ylab = latexstring("P \\: \\left[\\hbar / \\alpha_0 \\right]")
    @time begin
        anim = @animate for i = 1:nt
            tmp = heatmap(
                q,
                p,
                z_vals[:, :, i],
                seriescolor=:inferno,
                xlabel=xlab,
                ylabel=ylab,
                title=title,
                xlims=(minimum(q), maximum(q)),
                ylims=(minimum(p), maximum(p)),
                zlims=(z_min, z_max),
                right_margin=10.0Plots.mm,
            )
        end
        p = gif(anim, joinpath(plot_dump, name), fps=60)
    end
    os_display(p)
    return p
end