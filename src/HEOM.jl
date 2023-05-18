module HEOM
using StaticArrays, LinearAlgebra, SparseArrays, BandedMatrices
using SciMLBase: AbstractDiffEqLinearOperator

# Replacement for diff eq operators
include("central_finite_differences.jl")


export HEOM

end
