"""
    HEOM

Phase-space quantum dynamics. The package currently solves the Wigner–Moyal equation for a
particle of mass `m` in a potential `V(q)`,

    ∂W/∂t = -(p/m) ∂W/∂q + Σₛ cₛ V⁽²ˢ⁺¹⁾(q) ∂²ˢ⁺¹W/∂p²ˢ⁺¹,   cₛ = (-1)ˢ (ħ/2)²ˢ / (2s + 1)!,

by the method of lines. [`wigner_moyal_problem`](@ref) discretises phase space on a
[`PhaseSpaceGrid`](@ref) and returns an `ODEProblem` that any OrdinaryDiffEq solver can
integrate in time.
"""
module HEOM

using FFTW: plan_brfft, plan_rfft, rfftfreq
using ForwardDiff: ForwardDiff
using LinearAlgebra: kron, mul!
using SciMLBase: ODEProblem
using SparseArrays: SparseArrays, SparseMatrixCSC, sparse, spdiagm

export PhaseSpaceGrid, on_grid
export Spectral, FiniteDifference
export wigner_moyal_operator, wigner_moyal!, wigner_moyal_problem
export phase_space_integral, expectation, purity
export harmonic_potential, coherent_wigner, fock_wigner, cat_wigner, harmonic_evolution

include("grid.jl")
include("derivatives.jl")
include("wigner_moyal.jl")
include("observables.jl")
include("harmonic_oscillator.jl")

end
