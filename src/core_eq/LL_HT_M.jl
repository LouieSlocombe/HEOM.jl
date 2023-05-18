function prep_LL_HT_M_fd(q_vec, p_vec, v_vec, mass, h_bar, gamma, beta; apx=2)
    n = length(q_vec)
    dq = q_vec[2] - q_vec[1]
    dp = p_vec[2] - p_vec[1]

    # Make the large grid matrices
    Q = repeat(q_vec, 1, n)
    P = repeat(reshape(p_vec, 1, :), n, 1)

    # Calculate the terms
    Pm = @. P / mass
    mh_bar = h_bar^2 / 24
    t_diss = gamma
    t_deco = gamma * mass / beta

    # Make the derivative matrices
    DqW = zeros(n, n)
    DpW = zeros(n, n)
    DppW = zeros(n, n)
    DqV = zeros(n, n)
    DqqqV = zeros(n, n)
    DqqqW = zeros(n, n)
    Wp = zeros(n, n)
    DpWp = zeros(n, n)

    # Make the large potential matrix
    V_vec = repeat(v_vec, 1, n)

    d_dq = prepare_fd_dq(1, apx, dq, n)
    d_dp = prepare_fd_dp(1, apx, dp, n)
    d_dpp = prepare_fd_dp(2, apx, dp, n)
    d_dqqq = prepare_fd_dq(3, apx, dq, n)
    d_dppp = prepare_fd_dp(3, apx, dp, n)

    # Take potential derivative
    dq_fd!(DqV, V_vec, d_dq)
    dp_fd!(DqqqV, V_vec, d_dqqq)

    DqqqV = @. DqqqV * mh_bar

    p = (
        Pm=Pm,
        t_diss=t_diss,
        t_deco=t_deco,
        d_dq=d_dq,
        d_dp=d_dp,
        d_dpp=d_dpp,
        d_dqqq=d_dqqq,
        d_dppp=d_dppp,
        Wp=Wp,
        DqW=DqW,
        DpW=DpW,
        DpWp=DpWp,
        DppW=DppW,
        DqV=DqV,
        DqqqV=DqqqV,
        DqqqW=DqqqW,
    )
    return p
end

function LL_HT_M_fd_trunc!(du, u, p, t)

    dq_fd!(p.DqW, u, p.d_dq)
    dp_fd!(p.DpW, u, p.d_dp)

    # Dissipation term
    @. p.Wp = p * u
    dp_fd!(p.DpWp, p.pW, p.d_dp)
    # Decoherence term
    dp_fd!(p.DppW, p.W, p.d_dpp)
    # Put it all together
    @. du = -p.Pm * p.DqW + p.DqV * p.DpW + p.t_diss * p.DpWp + p.t_deco * p.DppW
end

function LL_HT_M_fd!(du, u, p, t)
    dq_fd!(p.DqW, u, p.d_dq)
    dp_fd!(p.DpW, u, p.d_dp)

    # Higher order terms
    dq_fd!(p.DqqqW, u, p.d_dqqq)

    # Dissipation term
    @. p.Wp = p * u
    dp_fd!(p.DpWp, p.pW, p.d_dp)
    # Decoherence term
    dp_fd!(p.DppW, p.W, p.d_dpp)
    # Put it all together
    @. du = -p.Pm * p.DqW + p.DqV * p.DpW - p.DqqqV * p.DqqqW + p.t_diss * p.DpWp + p.t_deco * p.DppW
end
