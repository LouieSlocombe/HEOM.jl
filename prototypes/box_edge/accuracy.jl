# Long-time accuracy of the hard cutoff and the auxiliary absorber, Spectral ±8/64.
# For a harmonic potential and depth ≥ 2 the root's mean and covariance are exact,
# so the Gaussian reference measures only discretisation and instability error.
#     julia --project=. accuracy.jl <bath> <depth>
include("remedies.jl")
using HEOM: phase_space_mean, phase_space_covariance

bname, depth = Symbol(ARGS[1]), parse(Int, ARGS[2])
bath = BATHS[bname]
grid = box(8.0, 64)
W0 = displaced_gaussian(grid)
mean0, cov0 = [1.0, 0.0], [0.5 0.0; 0.0 0.5]
op = heom_operator(grid; mass = 1.0, potential = HARMONIC, bath, depth, scaled = true)
times = collect(0.0:2.0:60.0)
moments = map(times) do t
    m, C = cold_memory_reference(t, mean0, cov0, bath; mass = 1.0, omega = 1.0)
    (m, Matrix(C))
end
limit = 1e3 * maximum(abs, W0)
function run(prob)
    solve(prob, Vern7(); abstol = 1e-10, reltol = 1e-10, saveat = times,
        unstable_check = (dt, u, p, t) -> !(maximum(abs, u) < limit))
end
hard = run(heom_problem(W0, (0.0, last(times)), op))
results = Dict{String,Any}("hard cutoff" => hard)
for (label, kw) in (("absorb q0=4 σ0=20", (absorb = 20.0, q0 = 4.0)),
                    ("absorb q0=5 σ0=20", (absorb = 20.0, q0 = 5.0)),
                    ("absorb q0=3 σ0=20", (absorb = 20.0, q0 = 3.0)))
    results[label] = run(remedied_problem(W0, (0.0, last(times)), remedied(op; kw...)))
end
println("== $bname depth $depth: moment error vs exact Gaussian (max over saved t in window) ==")
for (label, sol) in sort(collect(results); by = first)
    err(window) = maximum((max(maximum(abs, phase_space_mean(physical_wigner(U), grid) - moments[i][1]),
                               maximum(abs, phase_space_covariance(physical_wigner(U), grid) - moments[i][2]))
                           for (i, U) in enumerate(sol.u) if window[1] <= sol.t[i] <= window[2]); init = 0.0)
    status = HEOM.successful_retcode(sol) ? "bounded to t=60" : @sprintf("blew up at t=%.1f", sol.t[end])
    @printf("%-20s %-18s  t≤10 %.1e  t≤30 %.1e  t≤60 %.1e\n", label, status, err((0, 10)), err((0, 30)), err((0, 60)))
end
# Perturbation of the full Wigner function by the absorber, before the hard cutoff degrades.
for (label, sol) in sort(collect(results); by = first)
    label == "hard cutoff" && continue
    n = min(length(sol.u), length(hard.u), 6)   # t ≤ 10
    d = maximum(maximum(abs, physical_wigner(sol.u[i]) - physical_wigner(hard.u[i])) for i in 1:n) / maximum(abs, W0)
    @printf("%-20s max|W - W_hard|/max|W0| for t ≤ %.0f: %.1e\n", label, times[n], d)
end
