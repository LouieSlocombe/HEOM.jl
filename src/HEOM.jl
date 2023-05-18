module HEOM
using StaticArrays, LinearAlgebra, SparseArrays, BandedMatrices
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

# Replacement for diff eq operators to calculate the derivative of a vector/matrix
include("central_finite_differences.jl")

# Functions for plotting
include("phase_space_plots.jl")

# Functions for symbolics
include("symbolic.jl")

# Functions for derivatives
include("derivatives.jl")

# Linear-linear Coupling in High-Temperature Markovian Limit
include("core_eq_LL_HT_M.jl")

# Linear-linear Coupling in High-Temperature non-Markovian Limit

export HEOM

end
