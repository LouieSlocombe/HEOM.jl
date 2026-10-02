using Aqua
using HEOM
using LinearAlgebra
using OrdinaryDiffEqVerner
using SciMLBase: remake
using SparseArrays
using Test

"""
Maximum absolute difference between two arrays.
"""
max_error(a, b) = maximum(abs, a - b)

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

    @testset "Package quality" begin
        Aqua.test_all(HEOM)
    end
end
