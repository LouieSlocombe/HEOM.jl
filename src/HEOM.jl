module HEOM
using StaticArrays, LinearAlgebra, SparseArrays, BandedMatrices
using SciMLBase: AbstractDiffEqLinearOperator

# Replacement for diff eq operators to calculate the derivative of a vector/matrix
include("central_finite_differences.jl")

include("phase_space_plots.jl")




export HEOM

end
