"""
    HEOM

Phase-space quantum dynamics. The package solves the Wigner–Moyal equation for a
particle of mass `m` in a potential `V(q,t)`,

    ∂W/∂t = -(p/m) ∂W/∂q + Σₛ cₛ V⁽²ˢ⁺¹⁾(q) ∂²ˢ⁺¹W/∂p²ˢ⁺¹,   cₛ = (-1)ˢ (ħ/2)²ˢ / (2s + 1)!,

by the method of lines. [`wigner_moyal_problem`](@ref) discretises phase space on a
[`PhaseSpaceGrid`](@ref) and returns an `ODEProblem` that any OrdinaryDiffEq solver can
integrate in time.
[`caldeira_leggett_problem`](@ref) adds Markovian friction and thermal momentum diffusion.
[`heom_problem`](@ref) propagates a hierarchy of auxiliary Wigner functions for a
Gaussian bath coupled linearly to position, including Drude–Lorentz and underdamped
Brownian thermal baths and combinations of independent components. Balanced model-order
reduction compresses a bath decomposition to fewer hierarchy modes.
Initial-state helpers transform wavefunctions and density kernels, and prepare
numerical energy eigenstates and isolated Gibbs states in arbitrary potentials.
Analysis functions expose observables, grid diagnostics, trajectory summaries and populations.
Population-relaxation fits extract forward and backward inter-well transfer rates.
Time-dependent potentials and separable dipole drives act on the full hierarchy.
Linear-response helpers propagate dipole perturbations and compute absorption spectra.
Plotting recipes display Wigner functions, marginal densities and diagnostic trajectories
when Plots.jl is loaded. Animation helpers record Wigner and marginal trajectories.
"""
module HEOM

using FFTW: fftfreq, ifft, plan_brfft, plan_rfft, rfft, rfftfreq
using ForwardDiff: ForwardDiff
using LinearAlgebra: Symmetric, dot, eigen, kron, mul!
using RecipesBase: @recipe, @series, @userplot
using SciMLBase: AbstractODESolution, ODEFunction, ODEProblem
using SciMLBase:
    CallbackSet, DiscreteCallback, ReturnCode, solve, successful_retcode, terminate!
using SparseArrays: SparseArrays, SparseMatrixCSC, sparse, spdiagm

export PhaseSpaceGrid, on_grid
export Spectral, FiniteDifference
export wigner_moyal_operator, wigner_moyal!, wigner_moyal_problem
export TimeDependentPotential, DrivenPotential
export caldeira_leggett_operator, caldeira_leggett!, caldeira_leggett_problem
export ExponentialBath, drude_lorentz_bath
export drude_lorentz_pade_bath
export brownian_oscillator_bath, combine_baths, bath_correlation, bath_spectrum
export harmonic_covariance, hankel_singular_values, compress_bath
export heom_operator, heom!, heom_problem, hierarchy_indices, physical_wigner
export hierarchy_size, rescale_hierarchy
export equilibrate, EquilibriumResult
export linear_response_problem, linear_response, LinearResponseResult, absorption_spectrum
export phase_space_integral, expectation, purity, overlap, energy, wigner_negativity
export position_density, momentum_density, phase_space_mean, phase_space_covariance
export boundary_weight, spectral_tail, diagnostics
export probability, probability_current, probability_rate, expectation_rate
export tunnelling_rates
export harmonic_potential, coherent_wigner, fock_wigner, cat_wigner, harmonic_evolution
export wavefunction_wigner, density_matrix_wigner
export eigenstates, eigenstate_wigner, thermal_wigner
export wigneranimation, marginalanimation

include("grid.jl")
include("derivatives.jl")
include("wigner_moyal.jl")
include("driven.jl")
include("caldeira_leggett.jl")
include("heom.jl")
include("pade_bath.jl")
include("brownian_bath.jl")
include("composite_bath.jl")
include("bath_compression.jl")
include("bath_diagnostics.jl")
include("heom_solvers.jl")
include("equilibrium.jl")
include("observables.jl")
include("diagnostics.jl")
include("populations.jl")
include("tunnelling.jl")
include("harmonic_oscillator.jl")
include("initial_states.jl")
include("stationary_states.jl")
include("spectroscopy.jl")
include("plotting.jl")
include("animation.jl")

end
