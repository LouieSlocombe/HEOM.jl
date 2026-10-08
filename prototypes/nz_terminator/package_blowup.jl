# Does the package's own hard-cutoff HEOM (default Spectral discretisation, plain
# heom_problem) blow up from a physical initial state? Time at which max|U| first
# exceeds 1e3 times its initial value, integrating to t = 40. Run with
#     julia --project=prototypes/nz_terminator prototypes/nz_terminator/package_blowup.jl
using HEOM, OrdinaryDiffEqVerner, Printf
include("gaussian_reference.jl")

function blowup_time(prob, limit)
    sol = solve(prob, Vern7(); abstol = 1e-9, reltol = 1e-9, save_everystep = false,
        unstable_check = (dt, u, p, t) -> !(maximum(abs, u) < limit))
    return HEOM.successful_retcode(sol) ? Inf : sol.t[end]
end

harmonic = harmonic_potential(; mass = 1.0, omega = 1.0)
for (name, bath) in (
    ("warm (λ=0.2, γ=1, kT=1), Padé 1",
        drude_lorentz_pade_bath(; reorganization = 0.2, cutoff = 1.0, kT = 1.0, pade = 1)),
    ("cold weak (λ=0.05, γ=0.5, kT=0.1), Padé 2",
        drude_lorentz_pade_bath(; reorganization = 0.05, cutoff = 0.5, kT = 0.1, pade = 2)),
    ("cold strong (λ=0.8, γ=0.5, kT=0.1), Padé 2",
        drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 2)),
), (L, n) in ((6.0, 48), (8.0, 64))
    grid = PhaseSpaceGrid((-L, L), n, (-L, L), n)
    W0 = on_grid((q, p) -> gaussian_wigner(q, p; mean = [1.0, 0.0], covariance = [0.5 0; 0 0.5]), grid)
    times = map((2, 4, 6)) do depth
        prob = heom_problem(W0, (0.0, 40.0), grid; mass = 1.0, potential = harmonic, bath, depth, scaled = true)
        blowup_time(prob, 1e3 * maximum(abs, W0))
    end
    @printf("Spectral %-44s box ±%.0f n=%d  blow-up time at depth 2/4/6: %s\n",
        name, L, n, join((@sprintf("%.1f", t) for t in times), " / "))
    flush(stdout)
end
