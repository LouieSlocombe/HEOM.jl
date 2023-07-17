function plot_QSE_normalisation(
    q_vec,
    solu,
    time;
    name="QSE_deviation_norm.pdf",
    f_units="SI",
    dir=nothing
)
    # Track the deviation from normalisation conditon over time
    norm = QSE_norm_loop(q_vec, solu, time)

    if f_units == "SI"
        time, prefix, _ = best_time_units(time)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end
    ylab = latexstring("L_{2} \\; \\mathrm{Norm. \\: error}")

    # Loop over the time series
    y_vals = [abs.(i - 1.0) for i in norm]
    replace!(y_vals, Inf => NaN)
    fig = plot_general(time, y_vals, tlab, ylab)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_QSE_normalisation_log(
    q_vec,
    solu,
    time;
    name="QSE_deviation_norm_log.pdf",
    f_units="SI",
    dir=nothing
)
    # Track the deviation from normalisation conditon over time
    norm = QSE_norm_loop(q_vec, solu, time)

    if f_units == "SI"
        time, prefix, _ = best_time_units(time)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end
    ylab = latexstring("\\log \\: L_{2} \\; \\mathrm{Norm. \\: error}")

    # Loop over the time series
    y_vals = log10.([abs.(i - 1.0) for i in norm])
    replace!(y_vals, Inf => NaN)
    fig = plot_general(time, y_vals, tlab, ylab)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_QSE_P(q_vec, P; name="Pq.pdf", dir=nothing)
    fig = plot_general(
        q_vec,
        P,
        latexstring("Q \\: \\left[\\alpha_0 \\right]"),
        "P(q)",
    )
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_QSE_P_log(q_vec, P; name="Pq.pdf", max_y_min=-20, dir=nothing)
    y_vals = log10.(abs.(P))
    replace!(y_vals, Inf => NaN)
    # Get the maximum value for plotting
    y_min = minimum(y_vals)
    y_max = maximum(y_vals)
    # Prevent excessivly small values plotted
    if y_min < max_y_min
        y_min = max_y_min
    end

    fig = plot(
        q_vec,
        y_vals,
        ylims=(y_min, y_max),
        lw=3,
        xlabel=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
        ylabel="log P(q)",
        legend=false,
        linecolor=:black,
    )

    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_QSE_k_qm(time, k_qm; f_name="k_qm.pdf", dir=nothing)
    # Convert to better time array
    time, prefix, name = best_time_units(time)
    t_lab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    y_lab = latexstring("k_{\\mathrm{qm}}")
    fig = plot_general(time, k_qm, t_lab, y_lab)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_QSE_error_compare(
    q_vec,
    v,
    solu,
    time;
    name="QSE_error_compare.pdf",
    f_units="SI",
    dir=nothing,
    q0=nothing,
)
    # Track the deviation from normalisation conditon over time
    norm = QSE_norm_loop(q_vec, solu, time)

    # get the tunnelling probability
    if q0 == nothing
        p_p = calc_QSE_product_prob_loop(q_vec, solu, v, time)
    else
        p_p = calc_qse_k_qm(q_vec, q0, solu)
    end

    if f_units == "SI"
        time, prefix, _ = best_time_units(time)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end
    ylab = latexstring("\\log \\: L_{2} \\; \\mathrm{Norm. \\: error}")

    prob = log10.(abs.(p_p))
    error = log10.([abs.(i - 1.0) for i in norm])

    # Prevent crazy plotting
    replace!(error, Inf => NaN)
    replace!(prob, Inf => NaN)

    # Plot
    fig = plot(
        time,
        error,
        xlabel=tlab,
        ylabel=ylab,
        label="error",
        linecolor=:red,
        lw=3,
    )
    fig = plot!(time, prob, linecolor=:black, label="prob", lw=3)
    fig = plot!(legend=:bottomright)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_QSE_error_ratio(
    q_vec,
    v,
    solu,
    time;
    name="QSE_error_ratio.pdf",
    f_units="SI",
    dir=nothing
)
    # Track the deviation from normalisation conditon over time
    norm = QSE_norm_loop(q_vec, solu, time)

    # get the tunnelling probabilityF
    p_p = calc_QSE_product_prob_loop(q_vec, solu, v, time)

    if f_units == "SI"
        time, prefix, _ = best_time_units(time)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end
    ylab = "log error ratio"

    prob = abs.(p_p)
    error = [abs.(i - 1.0) for i in norm]
    ratio = log10.(prob ./ error)

    # Prevent crazy plotting
    replace!(ratio, Inf => NaN)

    # Plot
    fig = plot_general(time, ratio, tlab, ylab)

    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    os_display(fig)
    return fig
end

function plot_QSE_low_temp_terms(
    q_vec,
    v_vec,
    temperature,
    gamma,
    mass;
    f_approx=true,
    f_dutty=true,
    f_all=false,
    d_lim=nothing,
    D1_lim=nothing,
    D2_lim=nothing
)
    # Zero order terms
    d_0_val, D1_0_val, D2_0_val = QSE_low_temp_terms(
        0,
        q_vec,
        v_vec,
        temperature,
        gamma,
        mass;
        f_approx=f_approx,
        f_dutty=f_dutty
    )

    # First order terms
    d_1_val, D1_1_val, D2_1_val = QSE_low_temp_terms(
        1,
        q_vec,
        v_vec,
        temperature,
        gamma,
        mass;
        f_approx=f_approx,
        f_dutty=f_dutty
    )

    # Second order terms
    d_2_val, D1_2_val, D2_2_val = QSE_low_temp_terms(
        2,
        q_vec,
        v_vec,
        temperature,
        gamma,
        mass;
        f_approx=f_approx,
        f_dutty=f_dutty
    )

    x_lab = latexstring("Q \\: \\left[\\alpha_0 \\right]")

    # Make sure we can plot it ok
    replace!(d_0_val, Inf => NaN)
    replace!(D1_0_val, Inf => NaN)
    replace!(D2_0_val, Inf => NaN)
    replace!(d_1_val, Inf => NaN)
    replace!(D1_1_val, Inf => NaN)
    replace!(D2_1_val, Inf => NaN)
    replace!(d_2_val, Inf => NaN)
    replace!(D1_2_val, Inf => NaN)
    replace!(D2_2_val, Inf => NaN)

    if f_all
        os_display(plot_general(q_vec, d_0_val, x_lab, "d_0"))
        os_display(plot_general(q_vec, D1_0_val, x_lab, "D1_0"))
        os_display(plot_general(q_vec, D2_0_val, x_lab, "D2_0"))
        os_display(plot_general(q_vec, d_1_val, x_lab, "d_1"))
        os_display(plot_general(q_vec, D1_1_val, x_lab, "D1_1"))
        os_display(plot_general(q_vec, D2_1_val, x_lab, "D2_1"))
        os_display(plot_general(q_vec, d_2_val, x_lab, "d_2"))
        os_display(plot_general(q_vec, D1_2_val, x_lab, "D1_2"))
        os_display(plot_general(q_vec, D2_2_val, x_lab, "D2_2"))
    end

    # Set limits on the plotting ranges
    if d_lim == nothing
        d_lim = (-1.0 / (mass * gamma), 1.0 / (mass * gamma))
    end

    if D1_lim == nothing
        D1_lim = (minimum(D1_0_val), maximum(D1_0_val))
    end

    if D2_lim == nothing
        D2_lim = (-k_b * temperature, 3.0 * (k_b * temperature))
    end


    # Small d terms
    fig = plot(
        q_vec,
        d_0_val,
        xlabel=x_lab,
        ylabel="d terms",
        label="d_0",
        lw=3,
    )
    fig = plot!(fig, q_vec, d_1_val, label="d_1", lw=3)
    fig = plot!(fig, q_vec, d_2_val, label="d_2", lw=3, ylims=d_lim)
    os_display(fig)

    # Big D1 terms
    fig = plot(
        q_vec,
        D1_0_val,
        xlabel=x_lab,
        ylabel="D1 terms",
        label="D1_0",
        lw=3,
    )
    fig = plot!(fig, q_vec, D1_1_val, label="D1_1", lw=3)
    fig = plot!(fig, q_vec, D1_2_val, label="D1_2", lw=3, ylims=D1_lim)
    os_display(fig)

    # Big D2 terms
    fig = plot(
        q_vec,
        D2_0_val,
        xlabel=x_lab,
        ylabel="D2 terms",
        label="D2_0",
        lw=3,
    )
    fig = plot!(fig, q_vec, D2_1_val, label="D2_1", lw=3)
    fig = plot!(fig, q_vec, D2_2_val, label="D2_2", lw=3, ylims=D2_lim)
    D2_00 = QSE_d2_00(q_vec, v_vec, temperature, gamma, mass)
    fig = plot!(fig, q_vec, D2_00, label="D2_00", lw=3, ylims=D2_lim)
    os_display(fig)
end
