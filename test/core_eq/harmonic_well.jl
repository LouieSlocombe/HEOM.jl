using HEOM, Symbolics, ModelingToolkit, LinearAlgebra, Latexify, LaTeXStrings, Plots
using Test, OrdinaryDiffEq
const H = HEOM
# Set things up
tol = 1e-2

# Solver settings
algo = Tsit5()

# Time grid
t0 = 0.0
t1 = 1.0 / H.aut_2_ps
nt = 100
saveat = range(t0, stop=t1, length=nt)

# n = 2^8
# q_range = 12.0
# p_range = 60.0

# Make the grid
n = 2^8
q_range = 6.0
p_range = 80.0
q_vec, p_vec, Q, P, dq, dp = H.create_basis_even(n, q_range, p_range)

# Define the equation and its derivatives
@parameters t q p γ m β ω ħ ζ σ q0 p0 L α V0 x0
@variables W(..)
Dt = Differential(t)
Dq = Differential(q)
Dqq = Differential(q)^2
Dqqq = Differential(q)^3
Dp = Differential(p)
Dpp = Differential(p)^2
Dppp = Differential(p)^3

# Units
mass = 58752.0
h_bar = 1.0
si_k_b = 1.380649e-23 # Boltzmann constant
au_energy = 4.3597447222071e-18
k_b = si_k_b / au_energy #3.166811429e-6
temperature = 300.0
beta = 1.0 / (k_b * temperature)

# Morse potential
x_0 = -1.1 #-4.7
v_0 = 0.0220
a = 2.5
gamma = 1.0 / 141341.0
zeta = gamma

# Define the intial gaussian
sigma = 10.0
q_0 = 0.0
p_0 = 0.0

# Parameters
params = [γ => gamma, m => mass, β => beta, ħ => h_bar, ζ => zeta, σ => sigma, q0 => q_0, p0 => p_0, α => a, x0 => x_0, V0 => v_0]

# Define the potential
V(q, t) = 0
v_vec = H.make_discretised_potential(V, q, t, q_vec, params)

# Make symbolic intial gaussian
w_ic = 1 / (pi * ħ) * exp(-2 * σ * (q - q0)^2 - 1 / (2 * ħ^2 * σ) * (p - p0)^2)
w_ic = substitute(w_ic, params)
W0 = H.make_discretised_2d(w_ic, [q, p], q_vec, p_vec, [])
# H.plot_wigner_heatmap(q_vec, p_vec, W0; title="W0")