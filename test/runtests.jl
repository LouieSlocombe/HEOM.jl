using Aqua
using ADTypes
using HEOM
using LinearAlgebra
using LinearSolve
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqVerner
using SciMLBase: ODEProblem, remake
using SparseArrays
using Test

"""
Maximum absolute difference between two arrays.
"""
max_error(a, b) = maximum(abs, a - b)

"""
Normalised Gaussian Wigner function with the given mean and covariance.
"""
function gaussian_wigner(q, p; mean, covariance)
    δ = [q, p] - collect(mean)
    return exp(-dot(δ, covariance \ δ) / 2) / (2π * sqrt(det(covariance)))
end

"""
Evaluate the Wigner–Moyal right-hand side of `op` at `W` into a new matrix.
"""
function rhs(op, W)
    dW = similar(W)
    wigner_moyal!(dW, W, op, 0.0)
    return dW
end

"""
Bytes allocated by one evaluation of the right-hand side, after a warm-up call.
"""
function rhs_allocations(dW, W, op)
    wigner_moyal!(dW, W, op, 0.0)
    return @allocated wigner_moyal!(dW, W, op, 0.0)
end

@testset "HEOM" begin
    include("grid.jl")
    include("derivatives.jl")
    include("wigner_moyal_rhs.jl")
    include("harmonic_oscillator.jl")
    include("initial_states.jl")
    include("caldeira_leggett.jl")
    include("heom.jl")
    include("pade_bath.jl")
    include("brownian_bath.jl")
    include("composite_bath.jl")
    include("bath_diagnostics.jl")
    include("bath_degeneracy.jl")
    include("heom_scaling.jl")
    include("generalized_heom.jl")
    include("heom_solvers.jl")
    include("driven.jl")
    include("spectroscopy.jl")
    include("spectroscopy_benchmarks.jl")
    include("equilibrium.jl")
    include("analytic_benchmarks.jl")
    include("cold_strong_benchmarks.jl")
    include("hierarchy_stability.jl")
    include("equilibrium_harmonic.jl")
    include("bath_compression.jl")
    include("aaa_bath.jl")
    include("observables.jl")
    include("diagnostics.jl")
    include("populations.jl")
    include("tunnelling.jl")
    include("plotting.jl")
    include("animation.jl")

    @testset "Package quality" begin
        Aqua.test_all(HEOM)
    end
end
