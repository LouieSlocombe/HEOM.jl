# Moment error of the hard cutoff against the exact Gaussian moments, against time.
# Exponential growth from t = 0 at the generator's rate means the unstable modes are
# seeded by the physical state itself, not by rounding error at the edge.
include("common.jl")
using HEOM: phase_space_mean, phase_space_covariance
mean0, cov0 = [1.0, 0.0], [0.5 0.0; 0.0 0.5]
for (name, bath, L, n) in (("warm", BATHS.warm, 8.0, 64), ("cold_strong", BATHS.cold_strong, 6.0, 32),
                           ("cold_strong", BATHS.cold_strong, 8.0, 64), ("cold_weak", BATHS.cold_weak, 8.0, 64))
    grid = box(L, n)
    W0 = displaced_gaussian(grid)
    op = heom_operator(grid; mass = 1.0, potential = HARMONIC, bath, depth = 2, scaled = true)
    times = collect(0.0:1.0:30.0)
    sol = solve(heom_problem(W0, (0.0, 30.0), op), Vern7(); abstol = 1e-11, reltol = 1e-11, saveat = times,
        unstable_check = (dt, u, p, t) -> !(maximum(abs, u) < 1e6))
    outer = abs.(grid.q) .> L / 2
    cells = map(zip(sol.t, sol.u)) do (t, U)
        m, C = cold_memory_reference(t, mean0, cov0, bath; mass = 1.0, omega = 1.0)
        W = physical_wigner(U)
        e = max(maximum(abs, phase_space_mean(W, grid) - m), maximum(abs, phase_space_covariance(W, grid) - C))
        edge = maximum(abs, U[outer, :, :])
        @sprintf("t=%2.0f %.0e/%.0e", t, e, edge)
    end
    println("== $name depth 2 ±$L/$n  (moment error / max|U| over |q|>L/2) ==")
    for chunk in Iterators.partition(cells, 8)
        println("  ", join(chunk, "  "))
    end
    flush(stdout)
end
