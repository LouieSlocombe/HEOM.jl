using HEOM, Symbolics, ModelingToolkit, LinearAlgebra, Latexify, LaTeXStrings, Plots
using Test
const H = HEOM
# Set things up
tol = 1e-2

# n = 2^8
# q_range = 12.0
# p_range = 60.0

# Make the grid
n = 2^8
q_range = 5.0
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
q_0 = 0.5
p_0 = 0.0

# Parameters
params = [γ => gamma, m => mass, β => beta, ħ => h_bar, ζ => zeta, σ => sigma, q0 => q_0, p0 => p_0, α => a, x0 => x_0, V0 => v_0]

# Define the potential
# V(q, t) = V0 * (exp(-2 * α * (q - x0)) - 2 * exp(-α * (q - x0)))
V(q, t) = V0 * (1.0 - exp(-α * (q - x0)))^2

v_vec = H.make_discretised_potential(V, q, t, q_vec, params)
v_vec = v_vec .- minimum(v_vec)
# # plot the potential
# fig = plot(q_vec, v_vec, ylims=(0.0, 1.0)) # xlims=(-1, q_vec[end]) (0.0, 0.025)
# fig = plot!(xlabel="q", ylabel="V(q)", label="V(q)", legend=:bottomright, title="Morse potential")
# display(fig)


# Make symbolic intial gaussian
w_ic = 1 / (pi * ħ) * exp(-2 * σ * (q - q0)^2 - 1 / (2 * ħ^2 * σ) * (p - p0)^2)
w_ic = substitute(w_ic, params)
W0 = H.make_discretised_2d(w_ic, [q, p], q_vec, p_vec, [])
# H.plot_wigner_heatmap(q_vec, p_vec, W0; title="W0")

####################################################################################
# Check the generate_wm_eq vs wm_fd and wm_fft
println("Preparing WM eq")
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, V=V, Dqqq=Dqqq, Dppp=Dppp)
eq = H.generate_wm_eq(W, args; f_simple=true)
# Insert the initial condition W0
eq = H.eq_inserter(eq, W(q, p, t), w_ic)
# Substitute in the parameters
eq = substitute(eq, params)
# Make the discretised equation
eq_vec = H.make_discretised_2d(eq, [q, p], q_vec, p_vec, [])

# Prepare numerical finite difference eq
prep = H.prep_wm_fd(q_vec, p_vec, v_vec, mass, h_bar)
W_out = zeros(size(W0))
H.wm_fd!(W_out, W0, prep, 0.0)
# Check that they are the same
println("FD max difference: ", maximum(abs.(W_out .- eq_vec)))
@test ≈(W_out, eq_vec; atol=tol)

# Prepare numerical FFT eq
prep = H.prep_wm_fft(q_vec, p_vec, v_vec, W0, mass, h_bar)
W_out = zeros(size(W0))
H.wm_fft!(W_out, W0, prep, 0.0)
# Check that they are the same
println("FFT max difference: ", maximum(abs.(W_out .- eq_vec)))
@test W_out ≈ eq_vec atol = tol broken = true


####################################################################################
# Check the generate_wm_trunc_eq vs wm_fd_trunc and wm_fft_trunc
println("Preparing WM trunc eq")
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, V=V)
eq = H.generate_wm_trunc_eq(W, args; f_simple=true)
# Insert the initial condition W0
eq = H.eq_inserter(eq, W(q, p, t), w_ic)
# Substitute in the parameters
eq = substitute(eq, params)
# Make the discretised equation
eq_vec = H.make_discretised_2d(eq, [q, p], q_vec, p_vec, [])

# Prepare numerical finite difference eq
prep = H.prep_wm_fd(q_vec, p_vec, v_vec, mass, h_bar)
W_out = zeros(size(W0))
H.wm_fd_trunc!(W_out, W0, prep, 0.0)
# Check that they are the same
println("FD max difference: ", maximum(abs.(W_out .- eq_vec)))
@test ≈(W_out, eq_vec; atol=tol)

# Prepare numerical FFT eq
prep = H.prep_wm_fft(q_vec, p_vec, v_vec, W0, mass, h_bar)
W_out = zeros(size(W0))
H.wm_fft_trunc!(W_out, W0, prep, 0.0)
# Check that they are the same
println("FFT max difference: ", maximum(abs.(W_out .- eq_vec)))
@test W_out ≈ eq_vec atol = tol broken = true


####################################################################################
# Check the generate_LL_HT_M_eq vs LL_HT_M_fd and LL_HT_M_fft
println("Preparing LL_HT_M eq")
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, ħ=ħ, ζ=ζ, β=β, V=V, Dpp=Dpp, Dqqq=Dqqq, Dppp=Dppp)
eq = H.generate_LL_HT_M_eq(W, args; f_simple=true)
# Insert the initial condition W0
eq = H.eq_inserter(eq, W(q, p, t), w_ic)
# Substitute in the parameters
eq = substitute(eq, params)
# Make the discretised equation
eq_vec = H.make_discretised_2d(eq, [q, p], q_vec, p_vec, [])

# Prepare numerical finite difference eq
prep = H.prep_LL_HT_M_fd(q_vec, p_vec, v_vec, mass, h_bar, gamma, beta)
W_out = zeros(size(W0))
H.LL_HT_M_fd!(W_out, W0, prep, 0.0)
# Check that they are the same
println("FD max difference: ", maximum(abs.(W_out .- eq_vec)))
@test ≈(W_out, eq_vec; atol=tol)

# Prepare numerical FFT eq
prep = H.prep_LL_HT_M_fft(q_vec, p_vec, v_vec, W0, mass, h_bar, gamma, beta)
W_out = zeros(size(W0))
H.LL_HT_M_fft!(W_out, W0, prep, 0.0)
# Check that they are the same
println("FFT Max difference: ", maximum(abs.(W_out .- eq_vec)))
@test W_out ≈ eq_vec atol = tol broken = true


####################################################################################
# Check the generate_LL_HT_M_trunc_eq vs LL_HT_M_fd_trunc and LL_HT_M_fft_trunc
println("Preparing LL_HT_M trunc eq")
args = (t=t, q=q, p=p, m=m, Dq=Dq, Dp=Dp, Dpp=Dpp, ħ=ħ, ζ=ζ, β=β, V=V)
eq = H.generate_LL_HT_M_trunc_eq(W, args; f_simple=true)
# Insert the initial condition W0
eq = H.eq_inserter(eq, W(q, p, t), w_ic)
# Substitute in the parameters
eq = substitute(eq, params)
# Make the discretised equation
eq_vec = H.make_discretised_2d(eq, [q, p], q_vec, p_vec, [])

# Prepare numerical finite difference eq
prep = H.prep_LL_HT_M_fd(q_vec, p_vec, v_vec, mass, h_bar, gamma, beta)
W_out = zeros(size(W0))
H.LL_HT_M_fd_trunc!(W_out, W0, prep, 0.0)
# Check that they are the same
println("FD max difference: ", maximum(abs.(W_out .- eq_vec)))
@test ≈(W_out, eq_vec; atol=tol)

# Prepare numerical FFT eq
prep = H.prep_LL_HT_M_fft(q_vec, p_vec, v_vec, W0, mass, h_bar, gamma, beta)
W_out = zeros(size(W0))
H.LL_HT_M_fft_trunc!(W_out, W0, prep, 0.0)
# Check that they are the same
println("FFT max difference: ", maximum(abs.(W_out .- eq_vec)))
@test W_out ≈ eq_vec atol = tol broken = true