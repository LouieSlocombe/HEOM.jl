module HEOM
using StaticArrays, LinearAlgebra, SparseArrays, BandedMatrices
using SciMLOperators, FFTW, Plots, Symbolics, ModelingToolkit
using BenchmarkTools, LaTeXStrings, Latexify
using MethodOfLines, OrdinaryDiffEq, DomainSets
using SciMLBase: AbstractDiffEqLinearOperator

if Sys.iswindows()
    const plot_dump =
        joinpath(homedir(), "OneDrive - University of Surrey/DUMP/")
else
    # Set to headless plotting
    ENV["GKSwstype"] = 100
    const plot_dump = pwd()
end

# Include unit system
include("unit_system.jl")

# Functions for grid generation
include("grid.jl")

# Replacement for diff eq operators to calculate the derivative of a vector/matrix
include("central_finite_differences.jl")

# Functions for plotting
include("phase_space_plots.jl")

# Functions for symbolic utilities
include("symbolic_utils.jl")

# Functions for derivatives
include("derivatives.jl")

# Functions for symbolic forms of the equations
include("core_eq/symbolic_eq.jl")

# Wigner Moyal equation
include("core_eq/wigner_moyal.jl")

# Linear-linear Coupling in the High-Temperature Markovian Limit
include("core_eq/LL_HT_M.jl")

# Linear-linear Coupling in the High-Temperature non-Markovian Limit
# include("core_eq/LL_HT_NM.jl")

# Linear-linear Coupling in the Low-Temperature non-Markovian Limit
# include("core_eq/LL_LT_NM.jl")

export HEOM

end
