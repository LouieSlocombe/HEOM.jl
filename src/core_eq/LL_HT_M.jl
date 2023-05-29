####################################################################################
# Simple LL_HT_M equation using finite difference method
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
    DpppW = zeros(n, n)
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
    dq_fd!(DqqqV, V_vec, d_dqqq)

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
        DpppW=DpppW,
        P=P,
    )
    return p
end

function LL_HT_M_fd_trunc!(du, u, p, t)
    # Calculate dW_dq
    dq_fd!(p.DqW, u, p.d_dq)
    # Calculate dW_dp
    dp_fd!(p.DpW, u, p.d_dp)

    # Dissipation term
    @. p.Wp = p.P * u
    dp_fd!(p.DpWp, p.Wp, p.d_dp)
    # Decoherence term
    dp_fd!(p.DppW, u, p.d_dpp)
    # Put it all together
    @. du = -p.Pm * p.DqW + p.DqV * p.DpW + p.t_diss * p.DpWp + p.t_deco * p.DppW
end

function LL_HT_M_fd!(du, u, p, t)
    # Calculate dW_dq
    dq_fd!(p.DqW, u, p.d_dq)
    # Calculate dW_dp
    dp_fd!(p.DpW, u, p.d_dp)
    # Higher order terms
    dq_fd!(p.DpppW, u, p.d_dppp)

    # Dissipation term
    @. p.Wp = p.P * u
    dp_fd!(p.DpWp, p.Wp, p.d_dp)
    # Decoherence term
    dp_fd!(p.DppW, u, p.d_dpp)
    # Put it all together
    @. du = -p.Pm * p.DqW + p.DqV * p.DpW - p.DqqqV * p.DpppW + p.t_diss * p.DpWp + p.t_deco * p.DppW
end

####################################################################################
# Simple LL_HT_M equation using FFT method
function prep_LL_HT_M_fft(q_vec, p_vec, v_vec, W0, mass, h_bar, gamma, beta; apx=2)
    n = length(q_vec)

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
    DpppW = zeros(n, n)
    Wp = zeros(n, n)
    DpWp = zeros(n, n)

    # Make the large potential matrix
    V_vec = repeat(v_vec, 1, n)

    d_dq_ft, d_dq_FD, d_dq_ik = prepare_fft_dq(q_vec, W0; order=1)
    # dq_fft!(p.DqW, u, p.d_dq_ft, p.d_dq_FD, p.d_dq_ik)

    d_dp_ft, d_dp_FD, d_dp_ik = prepare_fft_dp(p_vec, W0; order=1)
    # dp_fft!(p.DpW, u, p.d_dp_ft, p.d_dp_FD, p.d_dp_ik)

    d_dpp_ft, d_dpp_FD, d_dpp_ik = prepare_fft_dp(p_vec, W0; order=1)
    # dp_fft!(p.DppW, u, p.d_dpp_ft, p.d_dpp_FD, p.d_dpp_ik)

    d_dqqq_ft, d_dqqq_FD, d_dqqq_ik = prepare_fft_dq(q_vec, W0; order=3)
    # dq_fft!(p.DqqqW, u, p.d_dqqq_ft, p.d_dqqq_FD, p.d_dqqq_ik)

    d_dppp_ft, d_dppp_FD, d_dppp_ik = prepare_fft_dp(p_vec, W0; order=3)
    # dp_fft!(p.DpppW, u, p.d_dppp_ft, p.d_dppp_FD, p.d_dppp_ik)

    # Take potential derivative
    dq_fd!(DqV, V_vec, d_dq)
    dq_fd!(DqqqV, V_vec, d_dqqq)

    DqqqV = @. DqqqV * mh_bar

    p = (
        Pm=Pm,
        t_diss=t_diss,
        t_deco=t_deco,
        d_dq_ft=d_dq_ft,
        d_dq_FDq=d_dq_FD,
        d_dq_ik=d_dq_ik,
        d_dp_ft=d_dp_ft,
        d_dp_FDq=d_dp_FD,
        d_dp_ik=d_dp_ik,
        d_dpp_ft=d_dpp_ft, 
        d_dpp_FD=d_dpp_FD, 
        d_dpp_ik=d_dpp_ik,
        d_dqqq_ft=d_dqqq_ft,
        d_dqqq_FDq=d_dqqq_FD,
        d_dqqq_ik=d_dqqq_ik,
        d_dppp_ft=d_dppp_ft,
        d_dppp_FDq=d_dppp_FD,
        d_dppp_ik=d_dppp_ik,
        Wp=Wp,
        DqW=DqW,
        DpW=DpW,
        DpWp=DpWp,
        DppW=DppW,
        DqV=DqV,
        DqqqV=DqqqV,
        DpppW=DpppW,
        P=P,
    )
    return p
end

function LL_HT_M_fft_trunc!(du, u, p, t)
    # Calculate dW_dq
    dq_fft!(p.DqW, u, p.d_dq_ft, p.d_dq_FD, p.d_dq_ik)
    # Calculate dW_dp
    dp_fft!(p.DpW, u, p.d_dp_ft, p.d_dp_FD, p.d_dp_ik)

    # Dissipation term
    @. p.Wp = p.P * u
    dp_fft!(p.DpWp, p.Wp, p.d_dp_ft, p.d_dp_FD, p.d_dp_ik)
    # Decoherence term
    dp_fft!(p.DppW, u, p.d_dpp_ft, p.d_dpp_FD, p.d_dpp_ik)
    # Put it all together
    @. du = -p.Pm * p.DqW + p.DqV * p.DpW + p.t_diss * p.DpWp + p.t_deco * p.DppW
end

function LL_HT_M_fft!(du, u, p, t)
    # Calculate dW_dq
    dq_fft!(p.DqW, u, p.d_dq_ft, p.d_dq_FD, p.d_dq_ik)
    # Calculate dW_dp
    dp_fft!(p.DpW, u, p.d_dp_ft, p.d_dp_FD, p.d_dp_ik)
    # Calculate d3W_dp3
    dp_fft!(p.DpppW, u, p.d_dppp_ft, p.d_dppp_FD, p.d_dppp_ik)

    # Dissipation term
    @. p.Wp = p.P * u
    dp_fft!(p.DpWp, p.Wp, p.d_dp_ft, p.d_dp_FD, p.d_dp_ik)
    # Decoherence term
    dp_fft!(p.DppW, u, p.d_dpp_ft, p.d_dpp_FD, p.d_dpp_ik)
    # Put it all together
    @. du = -p.Pm * p.DqW + p.DqV * p.DpW - p.DqqqV * p.DpppW + p.t_diss * p.DpWp + p.t_deco * p.DppW
end