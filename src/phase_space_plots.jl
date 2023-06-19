function os_display(p)
    """
    OS safe plot display, prevents issues with plotting when on a cluster
    and not on an X11 forwarding
    """
    if Sys.iswindows()
        display(p)
    end
end

function plot_wigner_heatmap(
    q,
    p,
    w;
    name="w_heatmap.pdf",
    title="",
    xlab=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
    ylab=latexstring("P \\: \\left[\\hbar / \\alpha_0 \\right]")
)
    # Fix the axes
    w = permutedims(w, (2, 1))
    # Plot the heatmap
    fig = heatmap(
        q,
        p,
        w,
        seriescolor=:inferno,
        xlabel=xlab,
        ylabel=ylab,
        title=title,
        xlims=(minimum(q), maximum(q)),
        ylims=(minimum(p), maximum(p)),
        right_margin=10.0Plots.mm,
    )
    # Save the figure
    savefig(fig, joinpath(plot_dump, name))
    os_display(fig)
    return fig
end

function animate_wigner_heatmap(
    q,
    p,
    sol;
    time=nothing,
    name="ani_wigner_heatmap.gif",
    title="",
    xlab=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
    ylab=latexstring("P \\: \\left[\\hbar / \\alpha_0 \\right]")
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

function plot_wigner_wq(
    q,
    p,
    W;
    qlab=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
    ylab=latexstring("\\int W \\, dp"),
    title="wq_expectation.pdf",
    yscale=:identity
)

    Wq, _ = calc_wigner_wqp(q, p, W)

    # Plot Q
    fig = plot_general(
        q,
        Wq,
        qlab,
        ylab,
    )
    fig = plot!(fig, yscale=yscale)
    os_display(fig)
    savefig(fig, joinpath(plot_dump, title))
    return fig
end

function plot_wigner_wp(
    q,
    p,
    W;
    plab=latexstring("P \\: \\left[\\hbar / \\alpha_0 \\right]"),
    ylab=latexstring("\\int W \\, dq"),
    title="wp_expectation.pdf",
    yscale=:identity
)

    _, Wp = calc_wigner_wqp(q, p, W)

    # Plot P
    fig = plot_general(
        p,
        Wp,
        plab,
        ylab,
    )
    fig = plot!(fig, yscale=yscale)
    os_display(fig)
    savefig(fig, joinpath(plot_dump, title))
end

function plot_wigner_wq_expectation(
    q,
    p,
    sol;
    ylab=latexstring("\\int \\int W q \\, dp"),
    title="wq_expectation_time.pdf",
    yscale=:identity,
    f_units="SI"
)
    time_sim = sol.t
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end

    # Calculate the expectation value
    q_ex, _ = [calc_wigner_wqp_expect(q, p, sol[i]) for i in sol]

    # Plot Q
    fig = plot_general(
        time_sim,
        q_ex,
        tlab,
        ylab,
    )
    fig = plot!(fig, yscale=yscale)
    os_display(fig)
    savefig(fig, joinpath(plot_dump, title))
    return fig
end

function plot_wigner_wp_expectation(
    q,
    p,
    sol;
    ylab=latexstring("\\int \\int W p \\, dp"),
    title="wp_expectation_time.pdf",
    yscale=:identity,
    f_units="SI"
)
    time_sim = sol.t
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end

    # Calculate the expectation value
    _, p_ex = [calc_wigner_wqp_expect(q, p, sol[i]) for i in sol]

    # Plot Q
    fig = plot_general(
        time_sim,
        p_ex,
        tlab,
        ylab,
    )
    fig = plot!(fig, yscale=yscale)
    os_display(fig)
    savefig(fig, joinpath(plot_dump, title))
    return fig
end

function plot_wigner_normalisation(    
    q,
    p,
    sol;
    ylab=latexstring("\\int \\int W p \\, dp"),
    title="wp_expectation_time.pdf",
    yscale=:identity,
    f_units="SI")
    time_sim = sol.t
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end
    norm = [int_2d(q, p, sol[i]) for i in sol]
    y_vals = [abs.(i - 1.0) for i in norm]
    replace!(y_vals, Inf => NaN)
    fig = plot_general(time_sim, y_vals, tlab, ylab)
    fig = plot!(fig, yscale=yscale)
    os_display(fig)
    savefig(fig, joinpath(plot_dump, title))
    return fig
end