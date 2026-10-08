# Long-time relaxation of the root's second moments at low depth, by time
# integration from a displaced Gaussian. Tests Fay's equilibrium claim without
# relying on a null vector, which box-edge modes contaminate. Run with
#     julia --project=prototypes/nz_terminator prototypes/nz_terminator/relaxation.jl
using OrdinaryDiffEqVerner, Printf
include("closure.jl")
include("gaussian_reference.jl")

function variances(W, grid)
    w = W .* (grid.dq * grid.dp)
    Q, P = grid.q, grid.p'
    n = sum(w)
    mq, mp = sum(w .* Q) / n, sum(w .* P) / n
    return sum(w .* (Q .- mq) .^ 2) / n, sum(w .* (P .- mp) .^ 2) / n
end

grid = PhaseSpaceGrid((-6.0, 6.0), 48, (-6.0, 6.0), 48)
W0 = on_grid((q, p) -> gaussian_wigner(q, p; mean = [1.0, 0.0], covariance = [0.5 0; 0 0.5]), grid)
harmonic = harmonic_potential(; mass = 1.0, omega = 1.0)
times = [10.0, 20.0, 30.0, 40.0]
for (name, bath) in (
    ("warm (λ=0.2, γ=1, kT=1), Padé 1",
        drude_lorentz_pade_bath(; reorganization = 0.2, cutoff = 1.0, kT = 1.0, pade = 1)),
    ("cold weak (λ=0.05, γ=0.5, kT=0.1), Padé 2",
        drude_lorentz_pade_bath(; reorganization = 0.05, cutoff = 0.5, kT = 0.1, pade = 2)),
)
    _, _, Σ = cold_memory_matrices(bath; mass = 1.0, omega = 1.0)
    println("\n== harmonic, $name ==")
    @printf("exact finite-bath stationary ⟨δq²⟩ = %.4f  ⟨δp²⟩ = %.4f\n", Σ[1, 1], Σ[2, 2])
    println("closure      depth  ⟨δq²⟩,⟨δp²⟩ at t = ", join(Int.(times), ", "))
    for depth in 0:2, kind in CLOSURES
        op = heom_operator(grid; mass = 1.0, potential = harmonic,
            discretization = FiniteDifference(4), moyal_terms = 1, bath, depth, scaled = true)
        sol = solve(terminated_problem(W0, (0.0, last(times)), terminate(op, kind)), Vern7();
            abstol = 1e-9, reltol = 1e-9, saveat = times, save_start = false,
            unstable_check = (dt, u, p, t) -> !(maximum(abs, u) < 1e3))
        cells = map(sol.u) do U
            vq, vp = variances(physical_wigner(U), grid)
            @sprintf("%6.3f,%6.3f", vq, vp)
        end
        HEOM.successful_retcode(sol) || push!(cells, @sprintf("blew up t=%.1f", sol.t[end]))
        @printf("%-12s %5d  %s\n", kind, depth, join(cells, "  "))
        flush(stdout)
    end
end
