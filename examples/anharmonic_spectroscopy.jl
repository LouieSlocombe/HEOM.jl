# Quartic oscillator: an anharmonic spectral shift and weak driven dynamics.
# Requires HEOM and OrdinaryDiffEqVerner in the active Julia environment.
# Quartic high-frequency modes can make explicit propagation comparatively slow.
using HEOM, OrdinaryDiffEqVerner, LinearAlgebra

mass, hbar = 1.0, 1.0
potential(q) = q^2 / 2 + 0.06q^4
grid = PhaseSpaceGrid((-6, 6), 56, (-6, 6), 48)
states = eigenstates(grid; mass, hbar, potential, nstates = 12)
ground = states.wavefunctions[:, 1]
W0 = wavefunction_wigner(ground, grid; hbar)
op = wigner_moyal_operator(grid; mass, hbar, potential)
dt = 0.04
response =
    linear_response(W0, (0.0, 24.0), op, Vern7(); saveat = dt, abstol = 1e-9, reltol = 1e-9)

# Independent spectral decomposition: R(t) = (2/ħ) sum_n |<n|q|0>|² sin(ω_n0*t).
gaps = (states.energies .- first(states.energies)) ./ hbar
dipoles = states.wavefunctions' * (grid.q .* ground) .* grid.dq
reference = [2 / hbar * sum(abs2.(dipoles) .* sin.(gaps .* t)) for t in response.times]
response_error = maximum(abs, response.response - reference)
@assert response_error < 3e-4
@assert abs(phase_space_integral(W0, grid) - 1) < 1e-12

broadening = 0.2
spectrum = absorption_spectrum(response; frequencies = 0.5:0.01:2.0, broadening)
reference_spectrum = absorption_spectrum(
    response.times,
    reference;
    frequencies = spectrum.frequencies,
    broadening,
)
peak = spectrum.frequencies[argmax(spectrum.intensity)]
reference_peak = reference_spectrum.frequencies[argmax(reference_spectrum.intensity)]
@assert abs(peak - reference_peak) <= 0.01
@assert abs(peak - gaps[2]) < 0.06
@assert gaps[2] > 1.1

# The induced dipole from a weak smooth field should equal the causal convolution
# ∫₀ᵗ R(t-s) E(s) ds. The even ground state has zero unperturbed mean position.
field(t) = 0.01 * exp(-((t - 1.0) / 0.35)^2)
driven = DrivenPotential(potential, field, identity)
sol = solve(
    wigner_moyal_problem(W0, (0.0, 4.0), grid; mass, hbar, potential = driven),
    Vern7();
    saveat = dt,
    abstol = 1e-9,
    reltol = 1e-9,
)
predicted = zeros(length(sol.t))
for i in 2:length(sol.t)
    integrand = [response.response[i-j+1] * field(sol.t[j]) for j in 1:i]
    predicted[i] = dt * (sum(integrand) - (first(integrand) + last(integrand)) / 2)
end
induced = [phase_space_mean(W, grid)[1] for W in sol.u]
drive_error = maximum(abs, induced - predicted)
@assert drive_error < 2e-5

println("First excitation gap / ħ: ", gaps[2], " (harmonic limit: 1)")
println(
    "Absorption peak: ",
    peak,
    " (finite-window eigenstate reference: ",
    reference_peak,
    ")",
)
println("Maximum Kubo response error: ", response_error)
println("Maximum weak-drive convolution error: ", drive_error)
# Converge position/momentum resolution and box size independently, then response
# duration and sampling. Broadening smooths the finite-time transform; it does not
# model a physical bath. For bath spectra use a converged full HEOM equilibrium.
