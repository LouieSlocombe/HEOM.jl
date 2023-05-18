function prepare_fft_der(x, L; order=1)
    # https://docs.sciml.ai/SciMLOperators/stable/tutorials/fftw/
    n = length(x)

    # The frequency modes sampled by our finite grid
    k = rfftfreq(n, 2 * pi * n / L) |> Array

    # Plan the fft
    transform = plan_rfft(x)
    # Define our wrapper for the FFT object
    T = FunctionOperator((du, u, p, t) -> mul!(du, transform, u), x, im * k;
        isinplace=true,
        T=ComplexF64,
        op_adjoint=(du, u, p, t) -> ldiv!(du, transform, u),
        op_inverse=(du, u, p, t) -> ldiv!(du, transform, u),
        op_adjoint_inverse=(du, u, p, t) -> ldiv!(du, transform, u)
    )

    # Make the operator
    ik = @. (im * k)^order
    # Diagonalise the operator
    ik_diag = DiagonalOperator(ik)
    # Make the derivative operator
    Dx = T \ ik_diag * T
    # Cache the operator
    Dx = cache_operator(Dx, x)
    return Dx
end

function prepare_fft_der_2d(x, y, Lx, Ly; order=1, axis=1)
    # https://docs.sciml.ai/SciMLOperators/stable/tutorials/fftw/
    nx = size(x, 1)
    ny = size(y, 1)
    # Assert that the grid is square
    @assert nx == ny

    # The frequency modes sampled by our finite grid
    kx = rfftfreq(nx, 2 * pi * nx / Lx) |> Array
    ky = rfftfreq(ny, 2 * pi * ny / Ly) |> Array

    m = length(kx)

    # reshape the grid
    kx = repeat(kx, 1, nx)
    #kx = repeat(reshape(kx, 1, :), m, 1)

    # Plan the fft
    transform = plan_rfft(x, axis)
    # Define our wrapper for the FFT object
    T = FunctionOperator((du, u, p, t) -> mul!(du, transform, u), x, im * kx;
        isinplace=true,
        T=ComplexF64, # ComplexF64
        op_adjoint=(du, u, p, t) -> ldiv!(du, transform, u),
        op_inverse=(du, u, p, t) -> ldiv!(du, transform, u),
        op_adjoint_inverse=(du, u, p, t) -> ldiv!(du, transform, u)
    )

    # Make the operator
    ik = @. (im * kx)^order
    if axis == 2
        ik = transpose(ik)
    end

    # Diagonalise the operator
    ik_diag = DiagonalOperator(ik)
    # Make the derivative operator
    Dx = T \ ik_diag * T
    # Cache the operator
    Dx = cache_operator(Dx, x)
    return Dx
end

function prepare_fft_der_2d_y(x, y, Lx, Ly; order=1)
    # https://docs.sciml.ai/SciMLOperators/stable/tutorials/fftw/
    nx = size(x, 1)
    ny = size(y, 1)
    # Assert that the grid is square
    @assert nx == ny

    # The frequency modes sampled by our finite grid
    ky = rfftfreq(ny, 2 * pi * ny / Ly) |> Array


    # reshape the grid
    ky = repeat(ky, 1, nx)

    # Plan the fft
    transform = plan_rfft(y, 1:2)
    # Define our wrapper for the FFT object
    T = FunctionOperator((du, u, p, t) -> mul!(du, transform, u), y, im * ky;
        isinplace=true,
        T=ComplexF64,
        op_adjoint=(du, u, p, t) -> ldiv!(du, transform, u),
        op_inverse=(du, u, p, t) -> ldiv!(du, transform, u),
        op_adjoint_inverse=(du, u, p, t) -> ldiv!(du, transform, u)
    )

    # Make the operator
    ik = @. (im * ky)^order

    # Diagonalise the operator
    ik_diag = DiagonalOperator(ik')
    println("ik_diag =", size(ik_diag))
    println("T =", size(T))
    # Make the derivative operator
    Dy = T \ ik_diag * T
    #Dy = T \  T
    # Cache the operator
    Dy = cache_operator(Dy, y)
    return Dy
end


function prepare_fft_dq_bad(q_vec, eq_vec; order=1)
    # Preparing the derivative operator
    dq = q_vec[2] - q_vec[1]
    q_range = q_vec[end] - q_vec[1] + dq

    nq = length(q_vec)
    kq = rfftfreq(nq, 2 * pi * nq / q_range) |> Array
    m = length(kq)

    kq = Array{ComplexF64}(repeat(kq, 1, nq))
    ft = plan_rfft(eq_vec, 1)
    FDq = Array{ComplexF64}(zeros(m, nq))

    ik = @. (im * kq)^order # make the operator
    ik = DiagonalOperator(ik) # diagonalise the operator
    return ft, FDq, ik
end

function dq_fft!(du, u, ft, fdq, ik)
    # Taking the derivative
    mul!(fdq, ft, u) # Take FFT
    mul!(fdq, ik, fdq) # Multiply by the diagonalised operator
    ldiv!(du, ft, fdq) # Take inverse FFT
end


function prepare_fft_dp_bad(p_vec, eq_vec; order=1)
    # Preparing the derivative operator
    dp = p_vec[2] - p_vec[1]
    p_range = p_vec[end] - p_vec[1] + dp

    np = length(p_vec)
    kp = rfftfreq(np, 2 * pi * np / p_range) |> Array
    m = length(kp)

    kp = Array{ComplexF64}(repeat(kp, 1, np))'
    ft = plan_rfft(eq_vec, 2)
    FDp = Array{ComplexF64}(zeros(np, m))

    ik = @. (im * kp)^order # make the operator
    ik = DiagonalOperator(ik) # diagonalise the operator
    return ft, FDp, ik
end

function dp_fft!(du, u, ft, fdp, ik)
    # Taking the derivative
    mul!(fdp, ft, u) # Take FFT
    mul!(fdp, ik, fdp) # Multiply by the diagonalised operator
    ldiv!(du, ft, fdp) # Take inverse FFT
end

function prepare_fd_dq(order, apx, dq, n)
    A = sparse(Array(CenteredDifference{1}(order, apx, dq, n)))
    A = A[1:end, 1:end.!=1]
    A = A[1:end, 1:end.!=end]
    A = Array{eltype(dq)}(A)
    return A
end

function dq_fd!(du, u, d_dq)
    return mul!(du, d_dq, u)
end


function prepare_fd_dp(order, apx, dp, n)
    A = sparse(Array(CenteredDifference{1}(order, apx, dp, n)))
    A = A[1:end, 1:end.!=1]
    A = A[1:end, 1:end.!=end]
    A = Array{eltype(dp)}(A)
    return A'
end

function dp_fd!(du, u, d_dp)
    return mul!(du, u, d_dp)
end