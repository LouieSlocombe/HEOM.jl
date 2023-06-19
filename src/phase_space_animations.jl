function animate_wigner_heatmap(
    q,
    p,
    sol;
    name="ani_wigner_heatmap.gif",
    title="",
    xlab=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
    ylab=latexstring("P \\: \\left[\\hbar / \\alpha_0 \\right]")
)
    """
    Animate wigner function heatmap
    """
    println("Animating wigner heatmap...")
    nt = length(sol.t)
    z_vals = real.(sol)

    # Get the min max range to fix the plot
    z_max = maximum(z_vals)
    z_min = minimum(z_vals)


    @time begin
        anim = @animate for i = 1:nt
            heatmap(
                q,
                p,
                permutedims(z_vals[:, :, i], (2, 1)),
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
        fig = gif(anim, joinpath(plot_dump, name), fps=60)
    end
    os_display(fig)
    return fig
end

function animate_wigner_wq(
    q,
    p,
    sol;
    name="wq.gif",
    xlab=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
    ylab=latexstring("Wq = \\int W \\, dp"),
    yscale=:identity
)
    println("Animating Wq probability density...")
    time_sim = sol.t
    nt = length(time_sim)

    # Get the values
    Wq = [abs.(calc_wigner_wqp(q, p, sol[i])[1]) for i = 1:nt]

    # # Get the maximum value for plotting
    y_max = maximum([maximum(Wq[i]) for i = 2:nt])
    y_min = minimum([minimum(Wq[i]) for i = 2:nt])

    @time begin
        # Loop over time
        anim = @animate for i = 1:nt
            plot(
                q,
                Wq[i, :],
                ylims=(y_min, y_max),
                lw=2,
                xlabel=xlab,
                ylabel=ylab,
                legend=false,
                linecolor=:black,
                yscale=yscale,
            )
        end
        fig = gif(anim, joinpath(plot_dump, name), fps=60)
    end
    os_display(fig)
    return fig
end

function animate_wigner_wp(
    q,
    p,
    sol;
    name="wp.gif",
    xlab=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
    ylab=latexstring("Wp = \\int W \\, dq"),
    yscale=:identity
)
    println("Animating Wp probability density...")
    time_sim = sol.t
    nt = length(time_sim)

    # Get the values
    Wp = [abs.(calc_wigner_wqp(q, p, sol[i])[2]) for i = 1:nt]

    # Get the maximum value for plotting
    y_max = maximum([maximum(Wp[i]) for i = 2:nt])
    y_min = minimum([minimum(Wp[i]) for i = 2:nt])

    @time begin
        # Loop over time
        anim = @animate for i = 1:nt
            plot(
                q,
                Wp[i, :],
                ylims=(y_min, y_max),
                lw=2,
                xlabel=xlab,
                ylabel=ylab,
                legend=false,
                linecolor=:black,
                yscale=yscale,
            )
        end
        fig = gif(anim, joinpath(plot_dump, name), fps=60)
    end
    os_display(fig)
    return fig
end