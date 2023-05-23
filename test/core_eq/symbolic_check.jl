using HEOM, Symbolics, ModelingToolkit, LinearAlgebra, Latexify, LaTeXStrings
using Test
const H = HEOM

# Set things up

# Make the grid
n = 2^9
q_range = 12.0
p_range = 60.0
q_vec, p_vec, Q, P, dq, dp = H.create_basis_even(n, q_range, p_range)

# Define the equation and its derivatives
@parameters t q p γ m β ω ħ ζ σ q0 p0 L V
@variables W(..)
Dt = Differential(t)
Dq = Differential(q)
Dqq = Differential(q)^2
Dqqq = Differential(q)^3
Dp = Differential(p)
Dpp = Differential(p)^2
Dppp = Differential(p)^3

####################################################################################
# Check the the Wigner master equation
args = (p = p, m = m, Dq = Dq, Dp = Dp, ħ = ħ, V=V, Dqqq = Dqqq, Dppp = Dppp)
eq = generate_wm_eq(W, args; f_simple=false)
println("eq = ", eq)
println(latexify(eq))