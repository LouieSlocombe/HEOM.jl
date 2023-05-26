using HEOM, Symbolics, ModelingToolkit, LinearAlgebra, Latexify, LaTeXStrings, Plots
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

# Parameters
params = [γ => gamma, m => mass, β => beta, ω => omega, ħ => h_bar, ζ => zeta, σ => sigma, q0 => q_0, p0 => p_0]


# Define the potential
V(q,t) = 0.5 * m * ω^2 * q^2
# # V_tmp = substitute(V(q,t) , Dict(m => mass, ω => omega))
# V_tmp = substitute(V(q,t), params)
# v_vec = H.make_discretised_2d(V_tmp, [q,t], q_vec, [0.0], [])
# v_vec = v_vec[:,1]


v_vec = H.make_discretised_potential(V, q, t, q_vec, params)

fig = plot(q_vec, v_vec, xlabel="q", ylabel="V(q)", title="Potential", legend=false)
display(fig)

# V_tmp = 0.5 * m * ω^2 * q^2
# v_vec = H.make_discretised_1d(substitute(V_tmp, Dict(t=>0.0,m => mass, ω => omega)), [q], q_vec, [])





# Make symbolic intial gaussian
w_ic = 1 / (pi * ħ) * exp(-2 * σ * (q - q0)^2 - 1 / (2 * ħ^2 * σ) * (p - p0)^2)
w_ic = exp(-q^2 - 0.05 * p^2)
println("w_ic = ", w_ic)
w_ic = substitute(w_ic, params)
W0 = H.make_discretised_2d(w_ic, [q, p], q_vec, p_vec, [])

####################################################################################
# Check the equation works for generate_wm_eq
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, V=V, Dqqq=Dqqq, Dppp=Dppp)
eq = H.generate_wm_eq(W, args; f_simple=true)
eq = H.eq_inserter(eq, W(q,p,t), w_ic)

println("eq = ", eq)
println(H.latexify_nice(eq))
# Substitute in the parameters
eq = substitute(eq, params)
println("eq = ", eq)
println(H.latexify_nice(eq))

eq_vec = H.make_discretised_2d(eq, [q, p], q_vec, p_vec, [])
H.plot_wigner_heatmap(q_vec, p_vec, eq_vec, title="W0")

# prepare numerical eq
prep = H.prep_wm_fd(q_vec, p_vec, v_vec, mass, h_bar)
W_out = zeros(size(W0))
H.wm_fd_trunc!(W_out, W0, prep, 0.0)
H.plot_wigner_heatmap(q_vec, p_vec, W_out, title="W0")


H.plot_wigner_heatmap(q_vec, p_vec, eq_vec .- W_out, title="W0")