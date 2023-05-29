using HEOM, Symbolics, ModelingToolkit, LinearAlgebra, Latexify, LaTeXStrings, Plots
using Test
const H = HEOM
# Set things up
tol = 1e-3

# Make the grid
n = 2^8
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

# Parameters
params = [γ => gamma, m => mass, β => beta, ω => omega, ħ => h_bar, ζ => zeta, σ => sigma, q0 => q_0, p0 => p_0]


# Define the potential
V(q, t) = 0.5 * m * ω^2 * q^2
v_vec = H.make_discretised_potential(V, q, t, q_vec, params)


# Make symbolic intial gaussian
w_ic = 1 / (pi * ħ) * exp(-2 * σ * (q - q0)^2 - 1 / (2 * ħ^2 * σ) * (p - p0)^2)
w_ic = exp(-q^2 - 0.05 * p^2)
w_ic = substitute(w_ic, params)
W0 = H.make_discretised_2d(w_ic, [q, p], q_vec, p_vec, [])

####################################################################################
# Check the generate_wm_eq and prep_wm_fd works
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, V=V, Dqqq=Dqqq, Dppp=Dppp)
eq = H.generate_wm_eq(W, args; f_simple=true)
# Insert the initial condition W0
eq = H.eq_inserter(eq, W(q, p, t), w_ic)

# Substitute in the parameters
eq = substitute(eq, params)

# Make the discretised equation
eq_vec = H.make_discretised_2d(eq, [q, p], q_vec, p_vec, [])

# prepare numerical eq
prep = H.prep_wm_fd(q_vec, p_vec, v_vec, mass, h_bar)
W_out = zeros(size(W0))
H.wm_fd_trunc!(W_out, W0, prep, 0.0)

# Check that they are the same
@test ≈(W_out, eq_vec; atol=tol)

####################################################################################
# Check the generate_wm_trunc_eq and prep_wm_fd works
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, V=V)
eq = H.generate_wm_trunc_eq(W, args; f_simple=true)
# Insert the initial condition W0
eq = H.eq_inserter(eq, W(q, p, t), w_ic)
# Substitute in the parameters
eq = substitute(eq, params)
# Make the discretised equation
eq_vec = H.make_discretised_2d(eq, [q, p], q_vec, p_vec, [])
# prepare numerical eq
W_out = zeros(size(W0))
H.wm_fd_trunc!(W_out, W0, prep, 0.0)

# Check that they are the same
@test ≈(W_out, eq_vec; atol=tol)

####################################################################################
# Check the generate_LL_HT_M_eq and prep_LL_HT_M_fd works
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, ζ=ζ, β=β, V=V, Dpp=Dpp, Dqqq=Dqqq, Dppp=Dppp)
eq = H.generate_LL_HT_M_eq(W, args; f_simple=true)
# Insert the initial condition W0
eq = H.eq_inserter(eq, W(q, p, t), w_ic)
# Substitute in the parameters
eq = substitute(eq, params)
# Make the discretised equation
eq_vec = H.make_discretised_2d(eq, [q, p], q_vec, p_vec, [])
# prepare numerical eq
prep = H.prep_LL_HT_M_fd(q_vec, p_vec, v_vec, mass, h_bar, gamma, beta)
W_out = zeros(size(W0))
H.LL_HT_M_fd!(W_out, W0, prep, 0.0)

# Check that they are the same
@test ≈(W_out, eq_vec; atol=tol)
println("Max difference: ", maximum(abs.(W_out .- eq_vec)))

# # plot the results
# H.plot_wigner_heatmap(q_vec, p_vec, W_out; title="numeric")
# H.plot_wigner_heatmap(q_vec, p_vec, eq_vec; title="symbolic")
# # plot the difference
# H.plot_wigner_heatmap(q_vec, p_vec, W_out .- eq_vec; title="difference")