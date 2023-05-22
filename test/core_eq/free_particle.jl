using HEOM, Symbolics, ModelingToolkit, LinearAlgebra
using Test
const H = HEOM

# Set things up



####################################################################################
# Check the equation works for a simple case

println("Running...")
slow_prob = ODEProblem(wigner_moyal, W, (t0, t1), p)
@time sol =
    solve(slow_prob, alg, saveat = saveat, progress = true, reltol = tol, abstol = tol)
