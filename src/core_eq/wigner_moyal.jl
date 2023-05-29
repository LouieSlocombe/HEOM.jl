####################################################################################
# Simple WM equation using finite difference method
function prep_wm_fd(q_vec, p_vec, v_vec, mass, h_bar; apx=2)
    n = length(q_vec)
    dq = q_vec[2] - q_vec[1]
    dp = p_vec[2] - p_vec[1]

    # Make the large grid matrices
    Q = repeat(q_vec, 1, n)
    P = repeat(reshape(p_vec, 1, :), n, 1)
    Pm = @. P / mass
    mh_bar = h_bar^2 / 24

    # Make the derivative matrices
    DqW = zeros(n, n)
    DpW = zeros(n, n)
    DqV = zeros(n, n)
    DqqqV = zeros(n, n)
    DqqqW = zeros(n, n)

    # Make the large potential matrix
    V_vec = repeat(v_vec, 1, n)

    # Prepare the derivative operators
    d_dq = prepare_fd_dq(1, apx, dq, n)
    d_dp = prepare_fd_dp(1, apx, dp, n)
    d_dqqq = prepare_fd_dq(3, apx, dq, n)
    d_dppp = prepare_fd_dp(3, apx, dp, n)

    # Take potential derivative
    dq_fd!(DqV, V_vec, d_dq)
    dq_fd!(DqqqV, V_vec, d_dqqq)

    DqqqV = @. DqqqV * mh_bar

    # Put everything together
    p = (
        Pm=Pm,
        d_dq=d_dq,
        d_dp=d_dp,
        d_dqqq=d_dqqq,
        d_dppp=d_dppp,
        DqW=DqW,
        DpW=DpW,
        DqV=DqV,
        DqqqV=DqqqV,
        DqqqW=DqqqW,
    )
    return p
end

function wm_fd_trunc!(du, u, p, t)
    # Calculate dW_dq
    dq_fd!(p.DqW, u, p.d_dq)
    # Calculate dW_dp
    dp_fd!(p.DpW, u, p.d_dp)

    # Main equation
    @. du = -p.Pm * p.DqW + p.DqV * p.DpW
end

function wm_fd!(du, u, p, t)
    # Calculate dW_dq
    dq_fd!(p.DqW, u, p.d_dq)
    # Calculate dW_dp
    dp_fd!(p.DpW, u, p.d_dp)
    # Calculate d3W_dp3
    dq_fd!(p.DqqqW, u, p.d_dqqq)

    # Main equation
    @. du = -p.Pm * p.DqW + p.DqV * p.DpW - p.DqqqV * p.DqqqW
end

####################################################################################
# Simple WM equation using FFT method
function prep_wm_fft(q_vec, p_vec, v_vec, W0, mass, h_bar)
    n = length(q_vec)
    dq = q_vec[2] - q_vec[1]
    dp = p_vec[2] - p_vec[1]

    # Make the large grid matrices
    Q = repeat(q_vec, 1, n)
    P = repeat(reshape(p_vec, 1, :), n, 1)
    Pm = @. P / mass
    mh_bar = h_bar^2 / 24

    # Make the derivative matrices
    DqW = zeros(n, n)
    DpW = zeros(n, n)
    DqV = zeros(n, n)
    DqqqV = zeros(n, n)
    DqqqW = zeros(n, n)

    # Make the large potential matrix
    V_vec = repeat(v_vec, 1, n)

    # Prepare the derivative operators
    d_dq_ft, d_dq_FDq, d_dq_ik = H.prepare_fft_dq(q_vec, W0; order=1)
    # H.dq_fft!(DqW, u, d_dq_ft, d_dq_FDq, d_dq_ik)

    d_dp_ft, d_dp_FDq, d_dp_ik = H.prepare_fft_dp(p_vec, W0; order=1)
    # H.dp_fft!(DpW, u, d_dp_ft, d_dp_FDq, d_dp_ik)

    d_dqqq_ft, d_dqqq_FDq, d_dqqq_ik = H.prepare_fft_dq(q_vec, W0; order=3)
    # H.dq_fft!(DqqqW, u, d_dqqq_ft, d_dqqq_FDq, d_dqqq_ik)

    d_dppp_ft, d_dppp_FDq, d_dppp_ik = H.prepare_fft_dp(p_vec, W0; order=3)
    # H.dp_fft!(DpppW, u, d_dppp_ft, d_dppp_FDq, d_dppp_ik)

    # Take potential derivative
    dq_fft!(DqV, V_vec, d_dq_ft, d_dq_FDq, d_dq_ik)
    dq_fd!(DqV, V_vec, d_dq)
    dq_fd!(DqqqV, V_vec, d_dqqq)

    DqqqV = @. DqqqV * mh_bar

    p = (
        Pm=Pm,
        d_dq=d_dq,
        d_dp=d_dp,
        d_dqqq=d_dqqq,
        d_dppp=d_dppp,
        DqW=DqW,
        DpW=DpW,
        DqV=DqV,
        DqqqV=DqqqV,
        DqqqW=DqqqW,
    )
    return p
end

function wm_fft_trunc!(du, u, p, t)
    # Calculate dW_dq
    dq_fd!(p.DqW, u, p.d_dq)
    # Calculate dW_dp
    dp_fd!(p.DpW, u, p.d_dp)

    # Main equation
    @. du = -p.Pm * p.DqW + p.DqV * p.DpW
end

function wm_fft!(du, u, p, t)
    # Calculate dW_dq
    dq_fd!(p.DqW, u, p.d_dq)
    # Calculate dW_dp
    dp_fd!(p.DpW, u, p.d_dp)
    # Calculate d3W_dp3
    dq_fd!(p.DqqqW, u, p.d_dqqq)

    # Main equation
    @. du = -p.Pm * p.DqW + p.DqV * p.DpW - p.DqqqV * p.DqqqW
end