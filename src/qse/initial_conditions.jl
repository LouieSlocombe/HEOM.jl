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