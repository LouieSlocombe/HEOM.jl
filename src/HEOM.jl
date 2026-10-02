"""
    HEOM

Phase-space quantum dynamics. The package solves the Wigner–Moyal equation for a
particle of mass `m` in a potential `V(q)`,

    ∂W/∂t = -(p/m) ∂W/∂q + Σₛ cₛ V⁽²ˢ⁺¹⁾(q) ∂²ˢ⁺¹W/∂p²ˢ⁺¹,   cₛ = (-1)ˢ (ħ/2)²ˢ / (2s + 1)!,

by the method of lines. [`wigner_moyal_problem`](@ref) discretises phase space on a
[`PhaseSpaceGrid`](@ref) and returns an `ODEProblem` that any OrdinaryDiffEq solver can
integrate in time.
[`caldeira_leggett_problem`](@ref) adds Markovian friction and thermal momentum diffusion.
Analysis functions expose observables, grid diagnostics, trajectory summaries and populations.
Plotting recipes display Wigner functions, marginal densities and diagnostic trajectories
when Plots.jl is loaded. Animation helpers record Wigner and marginal trajectories.
"""
module HEOM

using FFTW: plan_brfft, plan_rfft, rfft, rfftfreq
using ForwardDiff: ForwardDiff
using LinearAlgebra: dot, kron, mul!
using RecipesBase: @recipe, @series, @userplot
using SciMLBase: AbstractODESolution, ODEProblem
using SparseArrays: SparseArrays, SparseMatrixCSC, sparse, spdiagm

export PhaseSpaceGrid, on_grid
export Spectral, FiniteDifference
export wigner_moyal_operator, wigner_moyal!, wigner_moyal_problem
export caldeira_leggett_operator, caldeira_leggett!, caldeira_leggett_problem
export phase_space_integral, expectation, purity, overlap, energy, wigner_negativity
export position_density, momentum_density, phase_space_mean, phase_space_covariance
export boundary_weight, spectral_tail, diagnostics
export probability, probability_current, probability_rate, expectation_rate
export harmonic_potential, coherent_wigner, fock_wigner, cat_wigner, harmonic_evolution
export wigneranimation, marginalanimation

include("grid.jl")
include("derivatives.jl")
include("wigner_moyal.jl")
include("caldeira_leggett.jl")
include("observables.jl")
include("diagnostics.jl")
include("populations.jl")
include("harmonic_oscillator.jl")
include("plotting.jl")
include("animation.jl")

end
