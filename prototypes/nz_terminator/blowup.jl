# Practical stability: time at which each closed hierarchy first exceeds 1e3 times its
# initial maximum, integrating a physical Gaussian state to t = 40. Complements the
# eigenvalue study in stability.jl on the production-sized 48×48 grid. Run with
#     julia --project=prototypes/nz_terminator prototypes/nz_terminator/blowup.jl
using OrdinaryDiffEqVerner, Printf
include("closure.jl")
include("gaussian_reference.jl")

function blowup_time(cl, W0; T = 40.0)
    prob = terminated_problem(W0, (0.0, T), cl)
    limit = 1e3 * maximum(abs, W0)
    sol = solve(prob, Vern7(); abstol = 1e-9, reltol = 1e-9, save_everystep = false,
        unstable_check = (dt, u, p, t) -> !(maximum(abs, u) < limit))
    return HEOM.successful_retcode(sol) ? Inf : sol.t[end]
end

grid = PhaseSpaceGrid((-6.0, 6.0), 48, (-6.0, 6.0), 48)
W0 = on_grid((q, p) -> gaussian_wigner(q, p; mean = [1.0, 0.0], covariance = [0.5 0; 0 0.5]), grid)
harmonic = harmonic_potential(; mass = 1.0, omega = 1.0)
for (name, bath) in (
    ("warm (λ=0.2, γ=1, kT=1), Padé 1",
        drude_lorentz_pade_bath(; reorganization = 0.2, cutoff = 1.0, kT = 1.0, pade = 1)),
    ("cold weak (λ=0.05, γ=0.5, kT=0.1), Padé 2",
        drude_lorentz_pade_bath(; reorganization = 0.05, cutoff = 0.5, kT = 0.1, pade = 2)),
    ("cold strong (λ=0.8, γ=0.5, kT=0.1), Padé 2",
        drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 2)),
)
    println("\n== harmonic, $name: time to exceed 1e3 × initial (Inf = stable to t=40) ==")
    @printf("%5s %s\n", "depth", join((@sprintf("%12s", k) for k in CLOSURES)))
    for depth in 0:4
        op = heom_operator(grid; mass = 1.0, potential = harmonic,
            discretization = FiniteDifference(4), moyal_terms = 1, bath, depth, scaled = true)
        times = [blowup_time(terminate(op, kind), W0) for kind in CLOSURES]
        @printf("%5d %s\n", depth, join((@sprintf("%12.2f", t) for t in times)))
        flush(stdout)
    end
end
