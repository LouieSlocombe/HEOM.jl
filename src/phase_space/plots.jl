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
    ylab=latexstring("P \\: \\left[\\hbar / \\alpha_0 \\right]"),
    dir=nothing
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
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
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
    name="wq_expectation.pdf",
    yscale=:identity,
    dir=nothing
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
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_wigner_wp(
    q,
    p,
    W;
    plab=latexstring("P \\: \\left[\\hbar / \\alpha_0 \\right]"),
    ylab=latexstring("\\int W \\, dq"),
    name="wp_expectation.pdf",
    yscale=:identity,
    dir=nothing
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
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
end

function plot_wigner_wq_expectation(
    q,
    p,
    sol;
    ylab=latexstring("\\int \\int W q \\, dq dp"),
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
    q_ex, _ = [calc_wigner_wqp_expect(q, p, sol[i]) for i = 1:length(time_sim)]

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
    ylab=latexstring("\\int \\int W p \\, dq dp"),
    name="wp_expectation_time.pdf",
    yscale=:identity,
    f_units="SI",
    dir=nothing
)
    time_sim = sol.t
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end

    # Calculate the expectation value
    _, p_ex = [calc_wigner_wqp_expect(q, p, sol[i]) for i = 1:length(time_sim)]

    # Plot Q
    fig = plot_general(
        time_sim,
        p_ex,
        tlab,
        ylab,
    )
    fig = plot!(fig, yscale=yscale)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_wigner_normalisation(
    q,
    p,
    sol;
    ylab=latexstring("L_{2} \\; \\mathrm{Norm. \\: error}"),
    name="norm_error.pdf",
    yscale=:identity,
    f_units="SI",
    dir=nothing
)
    time_sim = sol.t
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end
    norm = [int_2d(q, p, sol[i]) for i = 1:length(time_sim)]
    y_vals = [abs.(i - 1.0) for i in norm]
    # Prevent log(0) errors
    replace!(y_vals, 0.0 => 1e-16)
    replace!(y_vals, Inf => NaN)
    fig = plot_general(time_sim, y_vals, tlab, ylab)
    fig = plot!(fig, yscale=yscale)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_wigner_purity(
    q,
    p,
    sol;
    name="wigner_purity.pdf",
    f_units="SI",
    yscale=:identity,
    ylab=latexstring("\\mathrm{Purity} \\: \\mathcal{P}"),
    dir=nothing
)
    # Get the time
    time_sim = sol.t
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end

    purity = [calc_wigner_purity(q, p, sol[i]) for i = 1:length(time_sim)]

    fig = plot_general(time_sim, purity, tlab, ylab)
    fig = plot!(fig, yscale=yscale)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_wigner_s2_entropy(
    q,
    p,
    sol;
    name="s2_entropy.pdf",
    f_units="SI",
    ylab=latexstring("\\mathrm{Entropy} \\, S_{2}"),
    yscale=:identity,
    dir=nothing
)
    # Get the time
    time_sim = sol.t
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end
    s2_entropy = [calc_wigner_s2_entropy(q, p, sol[i]) for i = 1:length(time_sim)]

    fig = plot_general(time_sim, s2_entropy, tlab, ylab)
    fig = plot!(fig, yscale=yscale)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_wigner_energy_expect(
    q,
    p,
    mass,
    v,
    sol;
    name="energy_expect.pdf",
    f_units="SI",
    yscale=:identity,
    ylab=latexstring(
        "\\mathrm{Energy \\; expectation}, \\; \\left< E(t)\\right>, \\; [E_\\mathrm{h}]",
    ),
    dir=nothing
)

    # Get the time
    time_sim = sol.t
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end

    # Get the energy expectation values
    e_expect = [calc_wigner_energy_expect(q, p, v, mass, sol[i]) for i = 1:length(time_sim)]

    # Plot them
    fig = plot_general(time_sim, e_expect, tlab, ylab)
    fig = plot!(fig, yscale=yscale)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_wigner_uncertainty_principle(
    q,
    p,
    sol;
    f_units="SI",
    ylab=latexstring("\\mathrm{Uncertainty}, \\, \\sigma_{Q}\\sigma_{P}"),
    title="uncertainty_principle.pdf",
    yscale=:identity,
    dir=nothing
)
    # Get the time
    time_sim = sol.t
    # Sort out the time saving
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end

    # Get the analytical variance
    uncert = [calc_wigner_uncertainty_principle(q, p, sol[i]) for i = 1:length(time_sim)]

    fig = plot_general(time_sim, uncert, tlab, ylab)
    fig = plot!(fig, yscale=yscale)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, title))
    else
        savefig(fig, joinpath(dir, title))
    end
    os_display(fig)
    return fig
end

function plot_wigner_step_occupation(
    q,
    p,
    q0,
    sol;
    f_units="SI",
    ylab=latexstring("\\int \\int W \\hat{h} \\, dq dp"),
    name="step_occupation.pdf",
    yscale=:identity,
    dir=nothing
)
    # Get the time
    time_sim = sol.t
    # Sort out the time saving
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end

    # Get the probability
    prob = [calc_wigner_probability(q, p, q0, sol[i]) for i = 1:length(time_sim)]

    fig = plot_general(time_sim, prob, tlab, ylab)
    fig = plot!(fig, yscale=yscale)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_k_qm(
    q,
    p,
    q0,
    sol;
    f_units="SI",
    ylab=latexstring("\\int \\int W \\, \\hat{h} \\, dq \\, dp"),
    name="k_qm.pdf",
    yscale=:identity,
    dir=nothing
)
    # Get the time
    time_sim = sol.t
    # Sort out the time saving
    if f_units == "SI"
        time_sim, prefix, _ = best_time_units(time_sim)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end

    # Get the probability flux
    k_qm = calc_wigner_k_qm(q, p, q0, sol)

    # Plot
    fig = plot_general(time_sim, k_qm, tlab, ylab)
    fig = plot!(fig, yscale=yscale)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end