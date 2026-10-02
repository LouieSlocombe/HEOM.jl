# A displaced Gaussian relaxes in a harmonic well under Caldeira–Leggett damping.
# Requires HEOM and OrdinaryDiffEqVerner in the active Julia environment.
using HEOM, OrdinaryDiffEqVerner, LinearAlgebra

mass, omega, hbar = 1.0, 1.0, 0.2
friction, kT = 0.8, 2.0             # kT / (hbar * omega) = 10: high temperature
grid = PhaseSpaceGrid((-10.0, 10.0), 64, (-10.0, 10.0), 64)
V = harmonic_potential(; mass, omega)

mu0 = [2.0, 0.0]
sigma0 = [0.5 0.0; 0.0 0.5]        # a mixed Gaussian, with det(sigma0) > hbar^2 / 4
sigma_inf = [kT/(mass*omega^2) 0.0; 0.0 mass*kT]
A = [0.0 1/mass; -mass*omega^2 -friction]

function gaussian_on_grid(mu, sigma, grid)
    precision = inv(sigma)
    normalization = 1 / (2π * sqrt(det(sigma)))
    return on_grid(grid) do q, p
        dq, dp = q - mu[1], p - mu[2]
        exponent = precision[1, 1] * dq^2 + 2precision[1, 2] * dq * dp
        exponent += precision[2, 2] * dp^2
        normalization * exp(-exponent / 2)
    end
end

W0 = gaussian_on_grid(mu0, sigma0, grid)
prob =
    caldeira_leggett_problem(W0, (0.0, 25.0), grid; mass, potential = V, friction, kT, hbar)
sol = solve(prob, Vern9(); saveat = 0.25, abstol = 1e-10, reltol = 1e-10)

# The exact Ornstein–Uhlenbeck solution applies to under-, critically and
# overdamped oscillators without separate formulas for the damping regimes.
exact_means = [exp(A * t) * mu0 for t in sol.t]
exact_covariances = map(sol.t) do t
    propagator = exp(A * t)
    sigma_inf + propagator * (sigma0 - sigma_inf) * propagator'
end
mean_error = maximum(
    maximum(abs, phase_space_mean(W, grid) - mu) for (W, mu) in zip(sol.u, exact_means)
)
covariance_error = maximum(
    maximum(abs, phase_space_covariance(W, grid) - sigma) for
    (W, sigma) in zip(sol.u, exact_covariances)
)
state_error = maximum(
    maximum(abs, W - gaussian_on_grid(mu, sigma, grid)) for
    (W, mu, sigma) in zip(sol.u, exact_means, exact_covariances)
)
@assert sol.t[end] == 25.0
@assert max(mean_error, covariance_error, state_error) < 1e-6

# The centroid approaches the well bottom; finite-temperature energy approaches
# kT, with nonzero thermal width, rather than the quantum ground-state energy.
@assert norm(phase_space_mean(sol.u[end], grid)) < 2e-4
@assert abs(energy(sol.u[end], grid; mass, potential = V) - kT) < 1e-6

println("Maximum mean error: ", mean_error)
println("Maximum covariance error: ", covariance_error)
println("Maximum Wigner-function error: ", state_error)
println("Final mean (q, p): ", phase_space_mean(sol.u[end], grid))
println("Final energy: ", energy(sol.u[end], grid; mass, potential = V), " (kT = ", kT, ")")

# Optional, with Plots installed:
# using Plots
# diagnosticsplot(sol; potential = V, fields = (:mean_q, :mean_p, :energy))
# wignerplot(sol)
