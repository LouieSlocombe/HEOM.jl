# A sinusoidally forced oscillator, its impulse response and absorption line.
# Requires HEOM and OrdinaryDiffEqVerner in the active Julia environment.
using HEOM, OrdinaryDiffEqVerner

mass, omega, hbar = 1.4, 1.2, 0.8
amplitude, frequency = 0.3, 0.75
grid = PhaseSpaceGrid((-7, 7), 64, (-7, 7), 64)
W0 = on_grid((q, p) -> coherent_wigner(q, p; mass, omega, hbar), grid)
V0 = harmonic_potential(; mass, omega)
field(t) = amplitude * cos(frequency * t)
potential = DrivenPotential(V0, field, identity)
# The general form is equivalent:
# potential = TimeDependentPotential((q, t) -> V0(q) - field(t)*q)

sol = solve(
    wigner_moyal_problem(W0, (0.0, 8.0), grid; mass, hbar, potential),
    Vern7();
    saveat = 0.05,
    abstol = 1e-10,
    reltol = 1e-10,
)
q_exact(t) =
    amplitude / (mass * (omega^2 - frequency^2)) * (cos(frequency * t) - cos(omega * t))
p_exact(t) =
    amplitude / (omega^2 - frequency^2) *
    (omega * sin(omega * t) - frequency * sin(frequency * t))
means = phase_space_mean.(sol.u, Ref(grid))
centroid_error = maximum(
    maximum(abs, mean - [q_exact(t), p_exact(t)]) for (t, mean) in zip(sol.t, means)
)
norm_drift = maximum(abs(phase_space_integral(W, grid) - 1) for W in sol.u)
@assert last(sol.t) ≈ 8.0
@assert centroid_error < 2e-7
@assert norm_drift < 1e-9

# Linear response uses the undriven generator and the equilibrium initial state.
op = wigner_moyal_operator(grid; mass, hbar, potential = V0)
response = linear_response(
    W0,
    (0.0, 30.0),
    op,
    Vern7();
    dipole = identity,
    saveat = 0.025,
    abstol = 1e-10,
    reltol = 1e-10,
)
response_error =
    maximum(abs, response.response - sin.(omega .* response.times) ./ (mass * omega))
@assert response_error < 2e-7

# Frequencies are angular frequencies; broadening is an inverse time.
broadening = 0.3
spectrum = absorption_spectrum(response; frequencies = 0.0:0.01:2.0, broadening)
analytic = @. 1 / (mass * (omega^2 - (spectrum.frequencies + im * broadening)^2))
@assert maximum(abs, spectrum.susceptibility - analytic) < 5e-4
peak = spectrum.frequencies[argmax(spectrum.intensity)]
println("Maximum forced-centroid error: ", centroid_error)
println("Maximum norm drift: ", norm_drift)
println("Maximum impulse-response error: ", response_error)
println(
    "Absorption peak: ",
    peak,
    " (broadened exact value ",
    sqrt(omega^2 + broadening^2),
    ")",
)
