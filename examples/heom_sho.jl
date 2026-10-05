# A harmonic oscillator coupled to a Drude–Lorentz bath in Wigner space.
# Requires HEOM and OrdinaryDiffEqVerner in the active Julia environment.
using HEOM, OrdinaryDiffEqVerner, LinearAlgebra

mass, omega, hbar = 1.0, 1.0, 1.0
reorganization, cutoff, kT = 0.12, 1.2, 0.8
grid = PhaseSpaceGrid((-8.0, 8.0), 48, (-8.0, 8.0), 48)
V = harmonic_potential(; mass, omega)
mu0 = [1.0, 0.0]
W0 = on_grid(
    (q, p) -> coherent_wigner(q, p; q0 = mu0[1], p0 = mu0[2], mass, omega, hbar),
    grid,
)

# The bath supplies the counterterm reorganization*q^2 internally. V is the
# physical harmonic potential. The default terminator supplies the diffusion
# contribution from the omitted Matsubara poles.
bath = drude_lorentz_bath(; reorganization, cutoff, kT, hbar, matsubara = 1)
op = heom_operator(grid; mass, potential = V, bath, depth = 5)
prob = heom_problem(W0, (0.0, 6.0), op)
sol = solve(prob, Vern7(); saveat = 0.1, abstol = 1e-9, reltol = 1e-9)

# Matrix W0 initializes the physical state; the higher ADOs start at zero.
# This factorized bare-bath initial condition includes the initial slip.
# For a harmonic well its exact centroid closes with one bath-force variable a:
# q' = p/m, p' = -(m*omega^2 + 2*reorganization)*q - a,
# a' = -cutoff*a - 2*reorganization*cutoff*q, with a(0) = 0.
A = [
    0.0 1/mass 0.0
    -(mass*omega^2+2reorganization) 0.0 -1.0
    -2reorganization*cutoff 0.0 -cutoff
]
initial_moments = [mu0; 0.0]
exact_means = [(exp(A*t)*initial_moments)[1:2] for t in sol.t]
means = [phase_space_mean(physical_wigner(U), grid) for U in sol.u]
mean_error =
    maximum(maximum(abs, mu - exact_mu) for (mu, exact_mu) in zip(means, exact_means))
d = diagnostics(sol; potential = V)
norm_drift = maximum(abs, d.norm .- first(d.norm))

# A second depth checks the effect of retaining more auxiliary states. This
# comparison concerns the whole Wigner function, not just its centroid.
coarse_op = heom_operator(grid; mass, potential = V, bath, depth = 3)
coarse_sol = solve(
    heom_problem(W0, prob.tspan, coarse_op),
    Vern7();
    saveat = sol.t,
    abstol = 1e-9,
    reltol = 1e-9,
)
depth_change = maximum(
    maximum(abs, physical_wigner(U) - physical_wigner(U_coarse)) for
    (U, U_coarse) in zip(sol.u, coarse_sol.u)
)

@assert sol.t[end] == 6.0
@assert coarse_sol.t[end] == 6.0
@assert mean_error < 1e-5
@assert norm_drift < 1e-8

println("Number of ADOs (including the physical state): ", length(hierarchy_indices(op)))
println("Maximum centroid error: ", mean_error)
println("Maximum norm drift: ", norm_drift)
println("Maximum Wigner-function change between depth 3 and 5: ", depth_change)
println("Final mean (q, p): ", means[end])
println("Final physical energy: ", d.energy[end])
println("Saved trajectory samples (time, mean q, mean p, norm):")
for i in 1:10:length(sol.t)
    println(round.((sol.t[i], means[i][1], means[i][2], d.norm[i]); digits = 6))
end

# The exact centroid alone does not test convergence of the full Wigner
# function. Repeat the depth comparison until the desired tolerance is reached,
# and vary matsubara, box size and grid resolution independently as well.
# Optional, with Plots installed:
# using Plots
# wignerplot(physical_wigner(sol), grid)
# diagnosticsplot(d; fields = (:mean_q, :mean_p, :energy))
