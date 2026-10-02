# A coherent state displaced in position in a simple harmonic oscillator.
# Requires HEOM, OrdinaryDiffEqVerner and Plots in the active Julia environment.
# Optionally pass an output directory as the first command-line argument.
using HEOM, OrdinaryDiffEqVerner, Plots

mass, omega, hbar = 1.0, 1.0, 1.0
x0, p0 = 2.0, 0.0
period = 2π / omega
grid = PhaseSpaceGrid((-8.0, 8.0), 64, (-8.0, 8.0), 64)

# HEOM calls the position coordinate q; here q is x.
W0 = on_grid((q, p) -> coherent_wigner(q, p; q0 = x0, p0, mass, omega, hbar), grid)
V = harmonic_potential(; mass, omega)
prob = wigner_moyal_problem(W0, (0.0, period), grid; mass, potential = V, hbar)
sol = solve(
    prob,
    Vern9();
    saveat = range(0.0, period; length = 81),
    abstol = 1e-11,
    reltol = 1e-11,
)

# Check the numerical trajectory against the known coherent-state motion.
means = [phase_space_mean(W, grid) for W in sol.u]
x_error = maximum(abs(means[i][1] - x0 * cos(omega * sol.t[i])) for i in eachindex(sol.t))
p_error = maximum(
    abs(means[i][2] + mass * omega * x0 * sin(omega * sol.t[i])) for i in eachindex(sol.t)
)
return_error = maximum(abs, sol.u[end] - W0)
@assert sol.t[end] ≈ period
@assert max(x_error, p_error, return_error) < 1e-7

output_dir = isempty(ARGS) ? mktempdir(; prefix = "heom-sho-") : abspath(only(ARGS))
mkpath(output_dir)

# Omit the repeated endpoint for a seamless four-second loop at 20 frames/s.
indices = 1:(length(sol.u)-1)
wigner = wigneranimation(
    sol;
    indices,
    size = (640, 540),
    xlabel = "x",
    ylabel = "p",
    colorbar_title = "W(x, p)",
)
gif(wigner, joinpath(output_dir, "sho_wigner.gif"); fps = 20)

marginals = marginalanimation(
    sol;
    indices,
    size = (800, 360),
    xlabel = ["x" "p"],
    ylabel = ["P(x)" "P(p)"],
    linewidth = 2,
)
gif(marginals, joinpath(output_dir, "sho_marginals.gif"); fps = 20)

println("Maximum mean-position error: ", x_error)
println("Maximum mean-momentum error: ", p_error)
println("State error after one period: ", return_error)
println("Animations saved to ", output_dir)
