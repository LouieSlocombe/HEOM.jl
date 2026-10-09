# Long-time moment error of the scalar Markov closure, alone and with the absorber.
# The closure only adds ∂p terms, so depth-2 root moments stay exact for a harmonic
# potential: any error is instability (or absorber perturbation).
#     julia --project=. accuracy_closure.jl <bath> <depth>
include("remedies.jl")
include(joinpath(@__DIR__, "..", "nz_terminator", "closure.jl"))
member_generator(op::HEOM.WignerHEOM{<:HEOM.SpectralWignerMoyal}) = spzeros(1, 1)
using HEOM: phase_space_mean, phase_space_covariance

struct Combined{C,R}
    closure::C
    remedy::R
end
function combined_heom!(dU, U, c::Combined, t)
    terminated_heom!(dU, U, c.closure, t)
    for a in 2:size(U, 3)
        @views dU[:, :, a] .-= c.remedy.σ .* U[:, :, a]
    end
end

bname, depth = Symbol(ARGS[1]), parse(Int, ARGS[2])
bath = BATHS[bname]
grid = box(8.0, 64)
W0 = displaced_gaussian(grid)
mean0, cov0 = [1.0, 0.0], [0.5 0.0; 0.0 0.5]
op = heom_operator(grid; mass = 1.0, potential = HARMONIC, bath, depth, scaled = true)
times = collect(0.0:2.0:60.0)
U0 = cat(W0, zeros(size(W0)..., length(op.indices) - 1); dims = 3)
limit = 1e3 * maximum(abs, W0)
run(f, p) = solve(HEOM.ODEProblem(f, U0, (0.0, last(times)), p), Vern7(); abstol = 1e-10, reltol = 1e-10,
    saveat = times, unstable_check = (dt, u, p, t) -> !(maximum(abs, u) < limit))
cases = (
    "hard cutoff" => () -> run(heom!, op),
    "markov_full" => () -> run(terminated_heom!, terminate(op, :markov_full)),
    "markov_full + absorb q0=4" => () -> run(combined_heom!, Combined(terminate(op, :markov_full), remedied(op; absorb = 20.0, q0 = 4.0))),
    "markov_full + absorb q0=5" => () -> run(combined_heom!, Combined(terminate(op, :markov_full), remedied(op; absorb = 20.0, q0 = 5.0))),
)
println("== $bname depth $depth: root moment error vs exact Gaussian, max over saved t in window ==")
for (label, f) in cases
    sol = f()
    err(lo, hi) = maximum((begin
        m, C = cold_memory_reference(t, mean0, cov0, bath; mass = 1.0, omega = 1.0)
        W = physical_wigner(U)
        max(maximum(abs, phase_space_mean(W, grid) - m), maximum(abs, phase_space_covariance(W, grid) - C))
    end for (t, U) in zip(sol.t, sol.u) if lo <= t <= hi); init = 0.0)
    status = HEOM.successful_retcode(sol) ? "bounded to t=60" : @sprintf("blew up at t=%.1f", sol.t[end])
    @printf("%-28s %-18s t≤2 %.1e  t≤10 %.1e  t≤30 %.1e  t≤60 %.1e\n", label, status, err(0, 2), err(0, 10), err(0, 30), err(0, 60))
    flush(stdout)
end
