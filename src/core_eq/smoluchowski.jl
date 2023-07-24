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

function prep_QSE_full(apx, q, v, mass, temperature, gamma)
    """
    Simple QSE
    """
    # Calculate missing values
    dq = q[2] - q[1]
    N = length(q)

    # Calculate the terms
    t1 = 1.0 / (mass * gamma)
    t2 = 1.0 / (k_b * temperature)

    # Make derivatives
    d_dq = prepare_fd_dq(1, apx, dq, N)

    dP_dq = zeros(N)
    inner = zeros(N)
    dPinner_dq = zeros(N)

    # Get the derivative of the potential
    dv1 = derivative_1d_interp(q, v, 1)

    p = (
        t1=t1,
        t2=t2,
        d_dq=d_dq,
        dP_dq=dP_dq,
        dv1=dv1,
        inner=inner,
        dPinner_dq=dPinner_dq,
    )
    return p
end

function QSE_eq_full!(du, u, p, t)
    """
    Simple QSE
    """
    # Calculate derivative of P
    mul!(p.dP_dq, p.d_dq, u)

    # Calculate inner
    @. p.inner = p.dv1 * u + p.t2 * p.dP_dq

    # Calculate derivative of inner
    mul!(p.dPinner_dq, p.d_dq, p.inner)

    # Put the equation together
    @. du = p.t1 * p.dPinner_dq
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