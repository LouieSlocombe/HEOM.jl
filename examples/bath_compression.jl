# Balanced compression of a cold Padé bath, with the checks a compressed bath needs.
# Run with HEOM and OrdinaryDiffEqVerner in the active Julia environment.
# The eight-mode reference hierarchy has 3003 members and takes several minutes.
using HEOM, OrdinaryDiffEqVerner, LinearAlgebra, Printf

# The cold, strongly coupled anharmonic oscillator of low_temperature_strong_coupling.jl.
mass, omega, hbar = 1.0, 1.0, 1.0
reorganization, cutoff, kT = 0.8, 0.5, 0.1
potential(q) = mass * omega^2 * q^2 / 2 + 0.08q^4
parameters = (; reorganization, cutoff, kT, hbar)

source = drude_lorentz_pade_bath(; parameters..., pade = 7)
σ = hankel_singular_values(source)
println("Hankel singular values of the 8-mode Padé bath:")
println("  ", join((@sprintf("%.2e", s) for s in σ), ", "))

candidates = (
    pade_3 = drude_lorentz_pade_bath(; parameters..., pade = 3),
    truncated_4 = compress_bath(source; modes = 4),
    residualized_4 = compress_bath(source; modes = 4, method = :residualize),
    truncated_3 = compress_bath(source; modes = 3),
    pade_7 = source,
)

# 1. Thermal spectrum S(ω) = 2ħJ(ω)/(1 - exp(-ħω/kT)), over the system's frequencies.
J(ω) = 2reorganization * cutoff * ω / (ω^2 + cutoff^2)
thermal(ω) =
    iszero(ω) ? 4reorganization * kT / cutoff : 2hbar * J(ω) / (-expm1(-hbar * ω / kT))
frequencies = range(-6, 6; length = 1201)
target = thermal.(frequencies)
spectral_error(bath) =
    maximum(abs, bath_spectrum.(Ref(bath), frequencies) - target) / maximum(target)

# 2. Harmonic equilibrium from the continuum fluctuation–dissipation theorem, with no
# bath decomposition: Σqq = (ħ/π)∫coth(ħΩ/2kT) Im χ(Ω) dΩ, and Σpp inserts (mΩ)².
# Gauss–Legendre nodes on Ω = tan(πu/2) cover the whole positive frequency axis.
function continuum_covariance(; nodes = 512)
    rule = eigen(SymTridiagonal(zeros(nodes), [j / sqrt(4j^2 - 1) for j in 1:(nodes-1)]))
    variances = zeros(2)
    for (u, weight) in zip((rule.values .+ 1) ./ 2, rule.vectors[1, :] .^ 2)
        Ω = tan(π * u / 2)
        χ = inv(mass * (omega^2 - Ω^2) - im * Ω * 2reorganization / (cutoff - im * Ω))
        integrand = hbar / π * imag(χ) / tanh(hbar * Ω / (2kT)) * (π / 2) * (1 + Ω^2)
        variances .+= weight .* [integrand, (mass * Ω)^2 * integrand]
    end
    return Diagonal(variances)
end
exact = continuum_covariance()
harmonic_error(bath) =
    maximum(abs, harmonic_covariance(bath; mass, omega) - exact) / maximum(exact)

println("\nBath checks (relative maximum errors):")
println("  bath            modes  bound*     spectrum   harmonic equilibrium")
for (name, bath) in pairs(candidates)
    modes = length(bath.rates)
    bound = 4sqrt(2) * sum(σ[(modes+1):end]; init = 0.0) / maximum(target)
    @printf(
        "  %-15s %5d  %s  %.2e   %.2e\n",
        name,
        modes,
        bath === source || name == :pade_3 ? "    —    " : @sprintf("%.2e", bound),
        spectral_error(bath),
        harmonic_error(bath)
    )
end
println("  * guaranteed bound on the spectral change from the 8-mode source")
# The source's own error is mostly its white-noise Padé remainder, which is why its
# equilibrium error is larger than its spectral error over this window.

# 3. Hierarchy dynamics. The quantity compared is the whole final Wigner function.
function cold_anharmonic_run(bath; depth = 6)
    grid = PhaseSpaceGrid((-7.0, 7.0), 40, (-7.0, 7.0), 40)
    initial = on_grid(
        (q, p) -> coherent_wigner(q, p; q0 = 0.7, p0 = -0.2, mass, omega, hbar),
        grid,
    )
    op = heom_operator(grid; mass, potential, bath, depth, scaled = true)
    seconds = @elapsed sol = solve(
        heom_problem(initial, (0.0, 0.8), op),
        Vern7();
        abstol = 1e-9,
        reltol = 1e-9,
        saveat = [0.8],
        save_start = false,
    )
    @assert last(sol.t) == 0.8
    state = copy(physical_wigner(sol))
    @assert abs(phase_space_integral(state, grid) - 1) < 1e-8
    return state, length(hierarchy_indices(op)), seconds
end

cold_anharmonic_run(candidates.truncated_3; depth = 1)  # compile before timing
reference, members, seconds = cold_anharmonic_run(source)
println("\nHierarchy dynamics at t = 0.8, depth 6, against the 8-mode source:")
@printf("  %-15s members %5d  %7.1f s\n", :pade_7, members, seconds)
states = Dict{Symbol,Matrix{Float64}}()
for name in (:pade_3, :truncated_4, :residualized_4, :truncated_3)
    states[name], auxiliaries, runtime = cold_anharmonic_run(candidates[name])
    @printf(
        "  %-15s members %5d  %7.1f s  max |ΔW| = %.2e\n",
        name,
        auxiliaries,
        runtime,
        maximum(abs, states[name] - reference)
    )
end
deeper, members, _ = cold_anharmonic_run(candidates.truncated_4; depth = 8)
depth_change = maximum(abs, deeper - states[:truncated_4])
@printf("  truncated_4 depth 6 → 8 (%d members): max |ΔW| = %.2e\n", members, depth_change)

compression_error(name) = maximum(abs, states[name] - reference)
# Equal mode counts: balanced modes are far more accurate than Padé poles.
@assert compression_error(:truncated_4) < compression_error(:pade_3) / 100
@assert compression_error(:residualized_4) < compression_error(:pade_3) / 100
@assert compression_error(:truncated_3) < compression_error(:pade_3) / 10
@assert spectral_error(candidates.truncated_4) < spectral_error(candidates.pade_3) / 10
@assert harmonic_error(candidates.truncated_4) < harmonic_error(candidates.pade_3) / 5
# Depth convergence must be repeated for the compressed bath.
@assert depth_change < 2e-6

# These checks establish this finite-time example. Use them again for other
# temperatures, potentials and times; converge the source decomposition first, since
# compression can only approach the source, not the physical bath.
