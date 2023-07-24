module HEOM
using StaticArrays, LinearAlgebra, SparseArrays, BandedMatrices
using SciMLOperators, FFTW, Plots, Symbolics, ModelingToolkit, Dierckx, Interpolations
using BenchmarkTools, LaTeXStrings, Latexify
using MethodOfLines, OrdinaryDiffEq, DomainSets
using SciMLBase: AbstractDiffEqLinearOperator
using Unitful, UnitfulEquivalences, UnitfulAtomic

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

# Replacement for diff eq operators to calculate the derivative of a vector/matrix
include("central_finite_differences.jl")

# Misc functions
include("misc.jl")

# Functions for symbolic utilities
include("symbolic_utils.jl")

# Functions for derivatives
include("derivatives.jl")

# Functions for integration
include("integration.jl")

# Core equations
##################################################################################
# Functions for symbolic forms of the equations
include("core_eq/symbolic_eq.jl")

# Wigner Moyal equation
include("core_eq/wigner_moyal.jl")

# Quantum Smoluchowski equation
include("core_eq/smoluchowski.jl")

# Linear-linear Coupling in the High-Temperature Markovian Limit
include("core_eq/LL_HT_M.jl")

# Linear-linear Coupling in the High-Temperature non-Markovian Limit
# include("core_eq/LL_HT_NM.jl")

# Linear-linear Coupling in the Low-Temperature non-Markovian Limit
# include("core_eq/LL_LT_NM.jl")

# Wigner phase space
##################################################################################
# Initial conditions Wigner phase space
include("phase_space/initial_conditions.jl")

# Functions for phase space grid generation
include("phase_space/grid.jl")

# Functions for plotting in phase space
include("phase_space/plots.jl")

# Functions for animating phase space in time
include("phase_space/animations.jl")

# Functions for the wigner tools
include("phase_space/tools.jl")


# QSE
##################################################################################
# Initial conditions Quantum Smoluchowski equation
include("qse/initial_conditions.jl")

# Quantum Smoluchowski equation tools
include("qse/tools.jl")

# Quantum Smoluchowski equation plotting
include("qse/plot.jl")

# Quantum Smoluchowski equation animation
include("qse/animation.jl")

export HEOM

end
