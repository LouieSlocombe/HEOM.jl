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
@parameters t q p γ m β ω ħ ζ σ q0 p0 L
@variables W(..)
Dt = Differential(t)
Dq = Differential(q)
Dqq = Differential(q)^2
Dqqq = Differential(q)^3
Dp = Differential(p)
Dpp = Differential(p)^2
Dppp = Differential(p)^3

# Units
mass = 1836.0
h_bar = 1.0
si_k_b = 1.380649e-23 # Boltzmann constant
au_energy = 4.3597447222071e-18
k_b = si_k_b / au_energy #3.166811429e-6
temperature = 300.0
beta = 1.0 / (k_b * temperature)

# Harmonic potential
omega = 0.005
gamma = omega * 1.0
zeta = gamma

# Define the intial gaussian
sigma = 10.0
q_0 = 0.5
p_0 = 0.0

# Define the potential
V(q,t) = 0.5 * m * ω^2 * q^2

# Parameters
params = [γ => gamma, m => mass, β => beta, ω => omega, ħ => h_bar, ζ => zeta, σ => sigma, q0 => q_0, p0 => p_0]

# Make symbolic intial gaussian
w_ic = 1 / (pi * ħ) * exp(-2 * σ * (q - q0)^2 - 1 / (2 * ħ^2 * σ) * (p - p0)^2)
w_ic = exp(-q^2 - 0.05 * p^2)
println("w_ic = ", w_ic)
# Substitute in the parameters
# w_ic = substitute(w_ic, Dict(ħ => h_bar, σ => sigma, q0 => q_0, p0 => p_0))
# println("w_ic = ", w_ic)

####################################################################################
# Check the equation works for a simple case
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, V=V, Dqqq=Dqqq, Dppp=Dppp)
eq = H.generate_wm_eq(W, args; f_simple=false)
println("eq = ", eq)
println(latexify(eq))
println(H.latexify_nice(eq))
eq = substitute(eq, W(q,p,t) => w_ic)
println("eq = ", eq)
println(latexify(eq))
println(H.latexify_nice(eq))