function convert_W_to_W_q(q_vec, p_vec, W)
    """
    Converts wigner distro to W_q by integrating over momentum
    """
    # Get matrix size
    N = length(q_vec)

    # Find W_q
    return [int_1d(p_vec, W[:, i]) for i = 1:N]
end

function QSE_norm(q_vec, P; k=5)
    return int_1d(q_vec, P; k=k)
end

function QSE_norm_loop(q_vec, solu, time; k=5)
    nt = length(time)
    norm = zeros(nt)
    for i = 1:nt
        norm[i] = QSE_norm(q_vec, solu[i][:]; k=k)
    end
    return norm
end

function calc_QSE_product_prob(q_vec, P, v; f_renorm=false, k=5)
    max_loc = find_double_well_barrier_loc(v)

    if f_renorm
        P = QSE_normalise(q_vec, P; k=k)
    end
    # Integrate for prob
    return int_1d(q_vec[max_loc:end], P[max_loc:end]; k=k)
end

function calc_QSE_product_prob_loop(
    q_vec,
    solu,
    v,
    time;
    f_renorm=false,
    k=5
)
    nt = length(time)
    Pp = zeros(nt)
    for i = 1:nt
        Pp[i] = calc_QSE_product_prob(
            q_vec,
            solu[i][:],
            v;
            f_renorm=f_renorm,
            k=k
        )
    end
    return Pp
end

function calc_QSE_react_prob(q_vec, P, v; f_renorm=false, k=5)
    max_loc = find_double_well_barrier_loc(v)

    if f_renorm
        P = QSE_normalise(q_vec, P; k=k)
    end
    # Integrate for prob
    return int_1d(q_vec[1:max_loc-1], P[1:max_loc-1]; k=k)
end

function calc_QSE_react_prob_loop(q_vec, solu, v, time; f_renorm=false, k=5)
    nt = length(time)
    Pr = zeros(nt)
    for i = 1:nt
        Pr[i] = calc_QSE_react_prob(
            q_vec,
            solu[i][:],
            v;
            f_renorm=f_renorm,
            k=k
        )
    end
    return Pr
end

function calc_QSE_reactive_correlation(q_vec, sol, v; f_renorm=false, k=5)
    # Get time
    time = sol.t
    # Get P(t)
    solu = sol.u

    dN_t = calc_QSE_product_prob_loop(
        q_vec,
        solu,
        v,
        time;
        f_renorm=f_renorm,
        k=k
    )

    # Take the first derivative
    return derivative_1d_interp(time, dN_t, 1; k=k)
end

function calc_QSE_reactive_correlation_eq(
    q_vec,
    v,
    sol,
    mass,
    temperature;
    k=5,
    f_print=true,
    f_low_t=false,
    f_renorm=true,
    f_denom=true,
    f_trunc=true
)
    """
    Calculate the reactive flux correlation function with
    thermodynamic equilibrium condition enforcement

    See papers
    2007 Craig J. Chem. Phys. 127, 144503 (2007);
    https://doi.org/10.1063/1.2772265

    See also
    Determines the forward rate equation using
    2009[Shi] Quantum rate dynamics for proton transfer reactions in
    condensed phase
    2011[Shi] Quantum rate dynamics for proton transfer reaction in
    a model system
    """

    # Get time
    time = sol.t
    solu = sol.u

    # Get the thermal state
    P_eq = QSE_thermal_distribution_wigner(
        q_vec,
        v,
        temperature
    )

    # Reactant and product partition functions
    Q_p = calc_QSE_product_prob(q_vec, P_eq, v; f_renorm=f_renorm, k=k)
    Q_r = 1.0 .- Q_p

    # Reactant and product populations
    p_p = calc_QSE_product_prob_loop(
        q_vec,
        sol,
        v,
        time;
        f_renorm=f_renorm,
        k=k
    )
    p_r = 1.0 .- p_p

    # Take the first derivative of the survival probability wrt. time
    dp_r_dt = derivative_1d_interp(time, p_r, 1; k=k)
    # Determine the thermal rate constants for the forward reactive processes
    # denom = @. (p_r - (Q_r / Q_p) * (1.0 - p_r))
    # rate = @. -1.0 * dp_r_dt / denom
    if f_denom
        denom = @. (p_r + (Q_r / Q_p) * (p_r - 1.0))
    else
        denom = 1.0
    end

    # Determine the rate
    rate = @. -1.0 * dp_r_dt / denom

    if f_print
        printout("Q_r = $(Q_r)")
        printout("Q_p = $(Q_p)")
        printout("denom = $(denom[end])")
    end

    # Prevent inf at zero
    rate[1] = rate[2]
    # Pack everything up
    p_rate = (time, p_r, log.(abs.(p_r)), dp_r_dt)
    return p_rate, abs.(rate)

end

function QSE_check_limits(temperature, gamma, q, v, mass)
    beta = 1.0 / (k_b * temperature)
    omega_0 = determine_best_local_spring(mass, v, q)
    printout(
        "γ / ω_0^2 ($(rnd(gamma / omega_0^2))) >> ħβ ($(rnd(h_bar * beta)))",
    )
    printout("γ^2 / ω_0^2 ($(rnd(gamma^2 / omega_0^2))) >> 1")
    printout("ħγβ ($(rnd(h_bar * gamma*beta))) >> 1")
    printout("QSE freq limit: γ/ω_0 ($(rnd(gamma/omega_0))) >>  1")
    printout("QSE low T limit: ħβω_0 ($(rnd(h_bar * beta *omega_0))) >> 1")
end

function q_correlation(q, v, temperature, mass, omega, gamma)
    """
    <q^2> q - correlation function

    2019[Ikeda] SI Low-Temperature Quantum Fokker−Planck and
    Smoluchowski Equations and Their Extension to Multistate Systems
    """
    beta = 1.0 / (k_b * temperature)
    omega_0 = determine_best_local_spring(mass, v, q)
    therm = @. 2.0 / (beta * h_bar * omega)
    nom = @. gamma * omega_0^2 * omega
    denom = @. (omega_0^2 - omega^2)^2 + gamma^2 * omega^2
    return @. 0.5 * therm / omega_0 * nom / denom
end

function harmonic_q2_anal(
    q,
    v,
    temperature,
    mass,
    gamma;
    k=5,
    n=Int(1e5),
    omega_min=0.001,
    omega_max=4.0
)
    """
    Analytical value of harmonic q spread

    2019[Ikeda] SI Low-Temperature Quantum Fokker−Planck and
    Smoluchowski Equations and Their Extension to Multistate Systems
    """
    omega_0 = determine_best_local_spring(mass, v, q)
    omega = linspace(omega_0 * omega_min, omega_0 * omega_max, n)
    q_corr = q_correlation(q, v, temperature, mass, omega, gamma)
    return 2.0 / pi * int_1d(omega, q_corr; k=k)
end

function QSE_thermal_distribution_harmonic(
    q,
    v,
    temperature,
    mass,
    gamma;
    k=5,
    n=Int(1e5),
    omega_min=0.001,
    omega_max=4.0
)
    """
    Analytical thermal solution of harmonic potential

    2019[Ikeda] SI Low-Temperature Quantum Fokker−Planck and
    Smoluchowski Equations and Their Extension to Multistate Systems
    """
    q2 = harmonic_q2_anal(
        q,
        v,
        temperature,
        mass,
        gamma;
        k=k,
        n=n,
        omega_min=omega_min,
        omega_max=omega_max
    )
    # Put the equation together
    P_ic = @. 1.0 / sqrt(2.0 * pi * q2) * exp(-q^2 / (2.0 * q2))
    # Normalise
    return QSE_normalise(q, P_ic; k=k)
end