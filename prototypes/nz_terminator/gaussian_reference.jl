# Exact Gaussian solution of a harmonic oscillator coupled to a specified finite
# exponential bath, including signed (negative) residues. Copied from
# test/cold_strong_benchmarks.jl so the prototype does not load the test suite.
using LinearAlgebra

function gaussian_wigner(q, p; mean, covariance)
    δ = [q, p] - collect(mean)
    return exp(-dot(δ, covariance \ δ) / 2) / (2π * sqrt(det(covariance)))
end

function cold_memory_matrices(bath; mass, omega)
    modes = length(bath.rates)
    dimension = modes + 2
    A, Q, initial = zeros(dimension, dimension),
    zeros(dimension, dimension),
    zeros(dimension, dimension)
    A[1, 2] = 1 / mass
    A[2, 1] = -mass * omega^2 - 2bath.counterterm
    Q[2, 2] = 2bath.diffusion
    damping = Matrix(Diagonal(bath.rates) + bath.mixing)
    A[2, 3:end] = bath.weights
    A[3:end, 1] = -2imag.(bath.coefficients) / bath.hbar
    A[3:end, 3:end] = -damping
    @assert all(iszero, bath.mixing) && all(isone, bath.weights)
    force_covariance = Matrix(Diagonal(real.(bath.coefficients)))
    initial[3:end, 3:end] = force_covariance
    Q[3:end, 3:end] = damping * force_covariance + force_covariance * damping'
    stationary = lyap(A, Q)
    return A, initial, stationary
end

function cold_memory_reference(t, mean, covariance, bath; mass, omega)
    A, initial, stationary = cold_memory_matrices(bath; mass, omega)
    initial[1:2, 1:2] = covariance
    F = exp(t * A)
    mean_t = F * [mean; zeros(length(bath.rates))]
    covariance_t = stationary + F * (initial - stationary) * F'
    return mean_t[1:2], Symmetric(covariance_t[1:2, 1:2])
end
