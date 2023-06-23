function QSE_lambda(temperature, gamma, mass)
    beta = 1.0 / (k_b * temperature)
    x = (h_bar * beta * gamma) / (2.0 * pi)
    pre = h_bar / (pi * mass * gamma)
    # This is the approximate value
    lambda_apr = pre * log(x)
    # This is the exact value, this seems to be much larger than approx
    # lambda_apr = @. pre * (digamma(1 + x) + C_e)
    return lambda_apr
end

function QSE_big_lambda(temperature, mass)
    """
    small lambda but in the high temperature limit

    2007[Coffey] Semiclassical Klein-kramers and smoluchowski equations
    for the brownian motion of a particle in an external potential
    """
    beta = 1.0 / (k_b * temperature)
    return h_bar^2 * beta^2 / (24.0 * mass)
end

function QSE_lambda_0_exact(temperature, gamma, mass)
    """
    From 2010[Maier] Quantum Smoluchowski equation a systematic study
    See https://en.wikipedia.org/wiki/Digamma_function#Asymptotic_expansion
    psi(x) ~ ln(x) - 1/2x
    """
    # Calculate matsubara frequency
    nu = (2.0 * pi * k_b * temperature) / h_bar
    # Calculate prefactor
    pre = h_bar / (pi * mass * gamma)

    # This is the exact value
    return pre * (digamma(1.0 + gamma / nu) + C_e)
end

function QSE_lambda_1_exact(temperature, gamma, mass)
    """
    From 2010[Maier] Quantum Smoluchowski equation a systematic study
    https://en.wikipedia.org/wiki/Trigamma_function#Computation_and_approximation
    """
    beta = 1.0 / (k_b * temperature)
    nu = (2.0 * pi * k_b * temperature) / h_bar
    pre = 2.0 / (mass * beta * nu^2 * gamma^2)
    part1 = (2.0 * nu / gamma) * (digamma(1.0 + gamma / nu) + C_e)
    part2 = trigamma(1.0 + gamma / nu) - (pi^2 / 6.0)
    return pre * (part1 - part2)
end

function QSE_d_0(q_vec, v, temperature, gamma, mass)
    I = ones(length(v))
    return @. 1.0 / (mass * gamma) * I
end

function QSE_d_1(q_vec, v, temperature, gamma, mass; k=5, f_approx=true)
    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    if f_approx
        lambda = QSE_lambda(temperature, gamma, mass)
    else
        lambda = QSE_lambda_0_exact(temperature, gamma, mass)
    end

    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)
    return @. 1.0 / (mass * gamma) * (1.0 - beta^2 / 4.0 * dv2^2 * lambda^2)
end

function QSE_d_2(q_vec, v, temperature, gamma, mass; k=5, f_approx=true)
    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    if f_approx
        lambda = QSE_lambda(temperature, gamma, mass)
    else
        lambda = QSE_lambda_0_exact(temperature, gamma, mass)
    end

    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)

    return @. 1.0 / (mass * gamma) * (1.0 - beta^2 / 4.0 * dv2^2 * lambda^2)
end

function QSE_d1_0(q_vec, v, temperature, gamma, mass; k=5)
    # Get first derivative
    return derivative_1d_interp(q_vec, v, 1; k=k)
end

function QSE_d1_1(q_vec, v, temperature, gamma, mass; k=5, f_approx=true)
    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    if f_approx
        lambda = QSE_lambda(temperature, gamma, mass)
    else
        lambda = QSE_lambda_0_exact(temperature, gamma, mass)
    end

    # Get first derivative
    dv1 = derivative_1d_interp(q_vec, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)

    return @. dv1 * (1.0 + beta^2 / 4.0 * dv2^2 * lambda^2)
end

function QSE_d1_2(q_vec, v, temperature, gamma, mass; k=5, f_approx=true)
    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    if f_approx
        lambda = QSE_lambda(temperature, gamma, mass)
    else
        lambda = QSE_lambda_0_exact(temperature, gamma, mass)
    end

    # Get first derivative
    dv1 = derivative_1d_interp(q_vec, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)
    # Get third derivative
    dv3 = derivative_1d_interp(q_vec, v, 3; k=k)

    t1 = @. dv2^2 - 2.0 * dv1 * dv3
    return @. dv1 * (1.0 + beta^2 / 4.0 * t1 * lambda^2)
end

function QSE_d2_00(q_vec, v, temperature, gamma, mass)
    beta = 1.0 / (k_b * temperature)
    I = ones(length(v))
    return @. I / beta
end

function QSE_d2_0(
    q_vec,
    v,
    temperature,
    gamma,
    mass;
    k=5,
    f_approx=true,
    f_dutty=true
)
    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    if f_approx
        lambda = QSE_lambda(temperature, gamma, mass)
    else
        lambda = QSE_lambda_0_exact(temperature, gamma, mass)
    end

    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)

    # Approximate form
    if f_dutty
        rtn = @. 1.0 / beta * (1.0 + beta * dv2 * lambda)
    else
        rtn = @. 1.0 / (beta - beta^2 * dv2 * lambda)
    end

    return rtn
end

function QSE_d2_1(q_vec, v, temperature, gamma, mass; k=5, f_approx=true)
    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    if f_approx
        lambda = QSE_lambda(temperature, gamma, mass)
    else
        lambda = QSE_lambda_0_exact(temperature, gamma, mass)
    end

    # Get first derivative
    dv1 = derivative_1d_interp(q_vec, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)

    t1 = @. 1.0 / beta
    t2 = @. 1.0 - beta * dv2 * lambda
    t3 = @. 3.0 / 4.0 * beta^2 * dv2^2 * lambda^2

    return @. t1 / (t2 + t3)
end

function QSE_d2_2(q_vec, v, temperature, gamma, mass; k=5, f_approx=true)
    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    if f_approx
        lambda = QSE_lambda(temperature, gamma, mass)
    else
        lambda = QSE_lambda_0_exact(temperature, gamma, mass)
    end

    # Get first derivative
    dv1 = derivative_1d_interp(q_vec, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)
    # Get third derivative
    dv3 = derivative_1d_interp(q_vec, v, 3; k=k)

    t1 = @. 1.0 / beta
    t2 = @. 1.0 - beta * dv2 * lambda
    t3 = @. 3.0 / 4.0 * dv2^2 + dv1 * dv3
    t4 = @. beta^2 * t3 * lambda^2

    return @. t1 / (t2 + t4)
end

function QSE_low_temp_terms(
    order,
    q_vec,
    v,
    temperature,
    gamma,
    mass;
    k=5,
    f_approx=true,
    f_dutty=true
)
    if order == 0
        d = QSE_d_0(q_vec, v, temperature, gamma, mass)
        D1 = QSE_d1_0(q_vec, v, temperature, gamma, mass; k=k)
        D2 = QSE_d2_0(
            q_vec,
            v,
            temperature,
            gamma,
            mass;
            k=k,
            f_approx=f_approx,
            f_dutty=f_dutty
        )
    elseif order == 1
        d = QSE_d_1(
            q_vec,
            v,
            temperature,
            gamma,
            mass;
            k=k,
            f_approx=f_approx
        )
        D1 = QSE_d1_1(
            q_vec,
            v,
            temperature,
            gamma,
            mass;
            k=k,
            f_approx=f_approx
        )
        D2 = QSE_d2_1(
            q_vec,
            v,
            temperature,
            gamma,
            mass;
            k=k,
            f_approx=f_approx
        )
    elseif order == 2
        d = QSE_d_2(
            q_vec,
            v,
            temperature,
            gamma,
            mass;
            k=k,
            f_approx=f_approx
        )
        D1 = QSE_d1_2(
            q_vec,
            v,
            temperature,
            gamma,
            mass;
            k=k,
            f_approx=f_approx
        )
        D2 = QSE_d2_2(
            q_vec,
            v,
            temperature,
            gamma,
            mass;
            k=k,
            f_approx=f_approx
        )
    else
        d = QSE_d_0(q_vec, v, temperature, gamma, mass)
        D1 = QSE_d1_0(q_vec, v, temperature, gamma, mass; k=k)
        D2 = QSE_d2_00(q_vec, v, temperature, gamma, mass)
    end
    return d, D1, D2
end

function convert_W_to_W_q(q_vec, p_vec, W)
    """
    Converts wigner distro to W_q by integrating over momentum
    """
    # Get matrix size
    N = length(q_vec)

    # Find W_q
    return [int_1d(p_vec, W[:, i]) for i = 1:N]
end

function QSE_normalise(q_vec, P; k=5)
    # Normalise
    N = int_1d(q_vec, P; k=k)
    return P ./ N
end

function QSE_thermal_distribution_wigner(
    q_vec,
    v,
    temperature;
    k=5
)
    """
    2010[Maier] Quantum QSE equation a systematic study
    """
    # Get the thermal energy
    beta = @. 1.0 / (k_b * temperature)
    # Get the thermal distribution
    P_ic = @. exp(-beta * v)

    # Normalise
    return QSE_normalise(q_vec, P_ic; k=k)
end

function QSE_thermal_distribution_zeng(
    q_vec,
    v,
    temperature,
    gamma,
    mass;
    k=5,
    f_approx=true
)
    """
    2010[Zeng] Dynamical properties of an asymmetric
    bistable system with quantum fluctuations in the strong-friction limit
    Eq 3
    """
    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    if f_approx
        lambda = QSE_lambda(temperature, gamma, mass)
    else
        lambda = QSE_lambda_0_exact(temperature, gamma, mass)
    end

    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)

    # Eq 3 from paper
    P_ic = @. beta * (1.0 - lambda * beta * dv2) * exp(-beta * v)

    # Normalise
    return QSE_normalise(q_vec, P_ic; k=k)
end

function QSE_thermal_distribution_maier_0(
    q_vec,
    v,
    temperature,
    gamma,
    mass;
    k=5,
    f_approx=true
)
    """
    2010[Maier] Quantum Smoluchowski equation a systematic study
    Eq 47
    """
    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    if f_approx
        lambda = QSE_lambda(temperature, gamma, mass)
    else
        lambda = QSE_lambda_0_exact(temperature, gamma, mass)
    end

    # Get first derivative
    dv1 = derivative_1d_interp(q_vec, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)

    # Eq 47 from paper
    P_quantum =
        @. (1.0 - beta * dv2 * lambda) * exp(0.5 * beta^2 * dv1^2 * lambda)
    P_ic = @. P_quantum * exp(-beta * v)

    # Normalise
    return QSE_normalise(q_vec, P_ic; k=k)
end

function QSE_thermal_distribution_maier_1(
    q_vec,
    v,
    temperature,
    gamma,
    mass;
    k=5
)
    """
    2010[Maier] Quantum Smoluchowski equation a systematic study
    Eq 54
    """
    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    l_0 = QSE_lambda_0_exact(temperature, gamma, mass)
    l_1 = QSE_lambda_1_exact(temperature, gamma, mass)

    # Get first derivative
    dv1 = derivative_1d_interp(q_vec, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)

    # Eq 54 from paper
    term1 = @. 1.0 -
               beta *
               dv2 *
               (
                   l_0 + (3.0 * dv2) / (4.0 * mass) * l_1 -
                   (3.0 * beta * dv2) / 4.0 * l_0^2
               )
    term2 = @. exp(
        0.5 * beta^2 * dv1^2 * (l_0 + dv2 * l_1 / mass - beta * dv2 * l_0^2),
    )

    # Put everything together
    P_quantum = @. term1 * term2
    P_ic = @. P_quantum * exp(-beta * v)

    # Normalise
    return QSE_normalise(q_vec, P_ic; k=k)
end

function QSE_thermal_distribution_coffey_1(
    q_vec,
    v,
    temperature,
    gamma,
    mass;
    k=5,
    f_low_t=false,
    f_approx=true
)
    """
    2007[Coffey] Semiclassical Klein-kramers and smoluchowski equations for the
    brownian motion of a particle in an external potential
    Eq. 7
    """
    # Get beta
    beta = 1.0 / (k_b * temperature)

    if f_low_t
        # Get small lamda in the low temperature
        if f_approx
            # Use approximate form
            lambda = QSE_lambda(temperature, gamma, mass)
        else
            # Use exact
            lambda = QSE_lambda_0_exact(temperature, gamma, mass)
        end
    else
        # Get big lamda
        lambda = QSE_big_lambda(temperature, mass)
    end

    # Get first derivative
    dv1 = derivative_1d_interp(q_vec, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)

    P_class = @. exp(-beta * v)
    # Eq. 7
    P_quantum = @. 1.0 + lambda * (beta * dv1^2 - 2.0 * dv2)
    P_ic = @. P_quantum * P_class

    # Normalise
    return QSE_normalise(q_vec, P_ic; k=k)
end

function QSE_thermal_distribution_coffey_2(
    q_vec,
    v,
    temperature,
    gamma,
    mass;
    k=5,
    f_low_t=false,
    f_approx=true
)
    """
    2006[Coffey] Wigner function approach to the quantum Brownian
    motion of a particle in a potential
    Eq. 9
    """
    # Get beta
    beta = 1.0 / (k_b * temperature)

    if f_low_t
        # Get small lamda in the low temperature
        if f_approx
            # Use approximate form
            lambda = QSE_lambda(temperature, gamma, mass)
        else
            # Use exact
            lambda = QSE_lambda_0_exact(temperature, gamma, mass)
        end
    else
        # Get big lamda
        lambda = QSE_big_lambda(temperature, mass)
    end

    # Get first derivative
    dv1 = derivative_1d_interp(q_vec, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)
    # Get third derivative
    dv3 = derivative_1d_interp(q_vec, v, 3; k=k)
    # Get fourth derivative
    dv4 = derivative_1d_interp(q_vec, v, 4; k=k)

    P_class = @. exp(-beta * v)
    # Eq. 7
    term_1 = @. 1.0 + lambda * (beta * dv1^2 - 2.0 * dv2)


    s1 = @. 36.0 * dv2^2
    s2 = @. 48.0 * dv3 * dv1
    s3 = @. 44.0 * beta * dv2 * dv1^2
    s4 = @. 5.0 * beta^2 * dv1^4
    s5 = @. 24.0 * dv4 / beta

    term_2 = @. lambda^2 / 10.0 * (s1 + s2 - s3 + s4 - s5)

    P_quantum = @. term_1 + term_2
    P_ic = @. P_quantum * P_class

    # Normalise
    return QSE_normalise(q_vec, P_ic; k=k)
end

function QSE_thermal_distribution(q_vec, v, temperature, gamma, mass; k=5)
    """
    2010[Maier] Quantum QSE equation a systematic study
    """
    # Get beta
    beta = 1.0 / (k_b * temperature)
    # Get lambda
    lambda = QSE_lambda(temperature, gamma, mass)

    # Get first derivative
    dv1 = derivative_1d_interp(q_vec, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)

    P_class = @. exp(-beta * v)
    # P_quantum = @. (1.0 - beta *dv2* lambda) * exp(0.5*beta^2 * dv1^2* lambda)
    P_quantum = @. (1.0 - beta * dv2 * lambda) * beta

    P_ic = @. P_quantum * P_class

    # Normalise
    return QSE_normalise(q_vec, P_ic; k=k)
end

function QSE_thermal_distribution_reactant(
    q_vec,
    v,
    temperature,
    gamma,
    mass;
    k=5
)
    # Get the thermal distribution
    P_ic = QSE_thermal_distribution(q_vec, v, temperature, gamma, mass; k=k)
    # Make heaviside step function
    h = make_small_h(q_vec, v)
    # Apply the step function
    P_react = @. P_ic * (1.0 - h)
    # Normalise
    return QSE_normalise(q_vec, P_react; k=k)
end

function QSE_thermal_distribution_product(
    q_vec,
    v,
    temperature,
    gamma,
    mass;
    k=5
)
    # Get the thermal distribution
    P_ic = QSE_thermal_distribution(q_vec, v, temperature, gamma, mass; k=k)
    # Make heaviside step function
    h = make_small_h(q_vec, v)
    # Apply the step function
    P_product = @. P_ic * h
    # Normalise
    return QSE_normalise(q_vec, P_product; k=k)
end

function QSE_thermal_prepare_reactant(q_vec, v, P_init; k=5)
    """
    Prepares P in the reactant well
    """
    # Make heaviside step function
    h = make_small_h(q_vec, v)
    # Apply the step function
    P_react = @. P_init * (1.0 - h)
    # Ensure normalisation
    return QSE_normalise(q_vec, P_react; k=k)
end

function QSE_thermal_prepare_product(q_vec, v, P_init; k=5)
    """
    Prepares P in the product well
    """
    # Make heaviside step function
    h = make_small_h(q_vec, v)
    # Apply the step function
    P_product = @. P_init * h
    # Ensure normalisation
    return QSE_normalise(q_vec, P_product; k=k)
end

function QSE_thermal_distribution_freq_indep(q_vec, v, temperature, mass; k=5)
    """
    2007[Coffey] Semiclassical Klein-kramers and QSE equations for the brownian motion of a particle in an external potential
    """
    # Get beta
    beta = 1.0 / (k_b * temperature)
    # Get lambda
    Lambda = (h_bar^2 * beta^2) / (24.0 * mass)

    # Get first derivative
    dv1 = derivative_1d_interp(q_vec, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q_vec, v, 2; k=k)

    P_class = @. exp(-beta * v)
    P_quantum = @. (1.0 + Lambda * (beta * dv1^2 - 2.0 * dv2))
    P_ic = @. P_quantum * P_class

    # Normalise
    return QSE_normalise(q_vec, P_ic; k=k)
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

function prep_QSE(apx, q, v, mass, temperature, gamma)
    """
    Simple QSE
    """
    T = eltype(q)

    # Obtain the thermal energy
    beta = @. 1.0 / (k_b * temperature)

    # Calculate missing values
    dq = q[2] - q[1]
    N = length(q)

    # Calculate the terms
    t1 = 1.0 / (mass * gamma)
    t2 = 1.0 / (mass * gamma * beta)
    # Make derivatives
    d_dq = prepare_fd_dq(1, apx, dq, N)
    d_dq2 = prepare_fd_dq(2, apx, dq, N)

    dP_dq = zeros(N)
    d2P_dq2 = zeros(N)

    # Get the derivative of the potential
    dv1 = derivative_1d_interp(q, v, 1)
    dv2 = derivative_1d_interp(q, v, 2)

    p = (
        t1=t1,
        t2=t2,
        d_dq=d_dq,
        dP_dq=dP_dq,
        d_dq2=d_dq2,
        d2P_dq2=d2P_dq2,
        dv1=dv1,
        dv2=dv2,
    )
    return p
end

function QSE_eq!(du, u, p, t)
    """
    Simple QSE
    """
    # Calculate dP_dq
    mul!(p.dP_dq, p.d_dq, u)
    # Calculate d2P_dq2
    mul!(p.d2P_dq2, p.d_dq2, u)
    # Put the equation together
    @. du = p.t1 * (p.dv2 * u + p.dv1 * p.dP_dq) + p.t2 * p.d2P_dq2
end

function prep_QSE_LT_leading_o(
    apx,
    q,
    v,
    mass,
    temperature,
    gamma,
    order;
    k=5,
    f_approx=true
)
    T = eltype(q)

    # Calculate missing values
    dq = q[2] - q[1]
    N = length(q)

    # Get beta
    beta = 1.0 / (k_b * temperature)

    # Get lambda
    if f_approx
        lambda = QSE_lambda(temperature, gamma, mass)
    else
        lambda = QSE_lambda_0_exact(temperature, gamma, mass)
    end
    # Get first derivative
    dv1 = derivative_1d_interp(q, v, 1; k=k)
    # Get second derivative
    dv2 = derivative_1d_interp(q, v, 2; k=k)

    dv2_term = @. 1.0 / beta * (1.0 + beta * dv2 * lambda)

    outer = @. 1.0 / (mass * gamma)

    # Make derivatives
    d_dq = prepare_fd_dq(1, apx, dq, N)
    dP_dq = zeros(N)
    dv1_P = zeros(N)
    dPinner_dq = zeros(N)
    inner = zeros(N)

    p = (
        d_dq=d_dq,
        dP_dq=dP_dq,
        dv1_P=dv1_P,
        dv1=dv1,
        inner=inner,
        dv2_term=dv2_term,
        dPinner_dq=dPinner_dq,
        outer=outer,
    )
    return p
end

function QSE_LT_eq_leading_o(du, u, p, t)
    # Calculate derivative of P
    mul!(p.dP_dq, p.d_dq, u)

    # Calculate dV2 P
    @. p.dv1_P = p.dv1 * u

    # Calculate inner
    @. p.inner = p.dv1_P + p.dv2_term * p.dP_dq

    # Calculate derivative of inner
    mul!(p.dPinner_dq, p.d_dq, p.inner)

    # Put the equation together
    @. du = p.outer * p.dPinner_dq
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
