# Independent Gaussian reference for the oscillator's finite exponential bath.
# This is the generalized Langevin construction used in test/analytic_benchmarks.jl:
# each positive real bath residue is an auxiliary Ornstein–Uhlenbeck force, with
# stationary initial variance real(cₖ) and retarded feedback -2imag(cₖ)q/ℏ.
# It retains the factorized initial slip and the bath counterterm from t = 0.
# Agreement tests the HEOM hierarchy and grid, not convergence of the bath poles.
using LinearAlgebra

function oscillator_gaussian_reference(times, mu0, covariance0, bath; mass, omega)
    @assert isfinite(mass) && mass > 0
    @assert isfinite(omega) && omega > 0
    @assert length(mu0) == 2 && all(isfinite, mu0)
    @assert size(covariance0) == (2, 2) && all(isfinite, covariance0)
    @assert issymmetric(covariance0) && isposdef(covariance0)
    @assert all(t -> isfinite(t) && t >= 0, times)
    # This simple stochastic embedding excludes negative residues and the
    # coupled correlation bases used near coincident Drude/Matsubara poles.
    @assert all(isone, bath.weights) && iszero(bath.mixing)
    @assert all(c -> real(c) > 0, bath.coefficients)

    modes = length(bath.rates)
    dim = modes + 2
    A, Q, initial_covariance = zeros(dim, dim), zeros(dim, dim), zeros(dim, dim)
    A[1, 2] = 1 / mass
    A[2, 1] = -mass * omega^2 - 2bath.counterterm
    Q[2, 2] = 2bath.diffusion
    initial_covariance[1:2, 1:2] = covariance0
    for k in 1:modes
        A[2, k+2] = 1
        A[k+2, 1] = -2imag(bath.coefficients[k]) / bath.hbar
        A[k+2, k+2] = -bath.rates[k]
        initial_covariance[k+2, k+2] = real(bath.coefficients[k])
        Q[k+2, k+2] = 2bath.rates[k] * real(bath.coefficients[k])
    end
    @assert all(z -> real(z) < 0, eigvals(A))

    # Solve AΣ∞ + Σ∞A' + Q = 0, then propagate only the decaying matrix exp(At).
    # Unlike a Van Loan block exponential, this does not construct exp(-A't),
    # whose rapidly growing bath modes lose accuracy on long trajectories.
    identity = Matrix{Float64}(I, dim, dim)
    stationary = reshape(-(kron(identity, A) + kron(A, identity)) \ vec(Q), dim, dim)
    stationary = (stationary + stationary') / 2
    @assert norm(A * stationary + stationary * A' + Q) < 1e-10 * max(1, norm(Q))
    initial_mean = [mu0; zeros(modes)]
    return map(collect(times)) do t
        propagator = exp(A * t)
        mean = propagator * initial_mean
        covariance =
            stationary + propagator * (initial_covariance - stationary) * propagator'
        reduced_covariance = Matrix(Symmetric(covariance[1:2, 1:2]))
        @assert all(isfinite, mean) && isposdef(reduced_covariance)
        (mean = mean[1:2], covariance = reduced_covariance)
    end
end
