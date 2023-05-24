using HEOM, Symbolics, ModelingToolkit, LinearAlgebra
using Test
const H = HEOM
# Set things up

# Make the grid
n = 2^9
q_range = 12.0
p_range = 60.0
q_vec, p_vec, Q, P, dq, dp = H.create_basis_even(n, q_range, p_range)

# Define the equation and its derivatives
@parameters t q p γ m β ω ħ ζ σ q0 p0 L
@variables W(..)
Dt = Differential(t)
Dq = Differential(q)
Dqq = Differential(q)^2
Dqqq = Differential(q)^3
Dp = Differential(p)
Dpp = Differential(p)^2
Dppp = Differential(p)^3

# Define the potential
v = 0.5 * m * ω^2 * q^2

# Parameters
params = [γ => gamma, m => mass, β => beta, ω => omega, ħ => h_bar, ζ => zeta, σ => sigma, q0 => q_0, p0 => p_0]

# Make symbolic intial gaussian
w_ic = 1 / (pi * ħ) * exp(-2 * σ * (q - q0)^2 - 1 / (2 * ħ^2 * σ) * (p - p0)^2)
println("w_ic = ", w_ic)
# Substitute in the parameters
w_ic = substitute(w_ic, Dict(ħ => h_bar, σ => sigma, q0 => q_0, p0 => p_0))
println("w_ic = ", w_ic)

####################################################################################
# Check the equation works for a simple case