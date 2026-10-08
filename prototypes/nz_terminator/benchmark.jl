# Accuracy and cost of final-tier closures against hierarchy depth. Run with
#     julia --project=prototypes/nz_terminator prototypes/nz_terminator/benchmark.jl [cases...]
# where cases are any of: warm coldweak coldstrong quartic (default: all).
#
# Every closure is compared with a deep hard-cutoff hierarchy on the SAME grid and
# bath, which isolates depth-truncation error: the only error a terminator targets.
# For the harmonic cases an exact Gaussian solution of the same finite bath also
# gives the grid-discretisation floor. Moments are uninformative here: for a
# harmonic oscillator, second moments are already exact at depth two.
#
# The cold, strong hierarchies have box-edge modes that grow (see stability.jl).
# A deep reference integrated to t = 10 is destroyed by them, so those cases stop
# at t = 2.5, before the growth reaches the solution.
using OrdinaryDiffEqVerner, Printf, Statistics
include("closure.jl")
include("gaussian_reference.jl")

const TOL = (; abstol = 1e-10, reltol = 1e-10)

function rhs_seconds(cl, U)
    dU = similar(U)
    terminated_heom!(dU, U, cl, 0.0)
    return median([@elapsed(terminated_heom!(dU, U, cl, 0.0)) for _ in 1:15])
end

function propagate(cl, W0, saves)
    prob = terminated_problem(W0, (0.0, last(saves)), cl)
    seconds = @elapsed sol = solve(prob, Vern7(); saveat = saves, save_start = false,
        TOL..., unstable_check = (dt, u, p, t) -> !(maximum(abs, u) < 1e6))
    HEOM.successful_retcode(sol) || return nothing, seconds, sol.stats.nf
    return [copy(physical_wigner(U)) for U in sol.u], seconds, sol.stats.nf
end

function run_case(name; bath, potential, moyal_terms, depths, reference_depth, exact, saves)
    grid = PhaseSpaceGrid((-6.0, 6.0), 48, (-6.0, 6.0), 48)
    mean0, covariance0 = [1.0, 0.0], [0.5 0.0; 0.0 0.5]
    W0 = on_grid((q, p) -> gaussian_wigner(q, p; mean = mean0, covariance = covariance0), grid)
    common = (; mass = 1.0, potential, discretization = FiniteDifference(4), moyal_terms)
    build(depth) = heom_operator(grid; common..., bath, depth, scaled = true)

    println("\n== $name, t ≤ $(last(saves)) ==")
    println("rates    ", round.(bath.rates; sigdigits = 4))
    println("residues ", round.(bath.coefficients; sigdigits = 4), "  D = ", bath.diffusion)
    reference, ref_seconds, _ = propagate(terminate(build(reference_depth), :none), W0, saves)
    check, _, _ = propagate(terminate(build(reference_depth - 2), :none), W0, saves)
    scale = maximum(abs, W0)
    deviation(Ws) = maximum(maximum(abs, W - R) for (W, R) in zip(Ws, reference)) / scale
    @printf("reference depth %d (%.0f s); change from depth %d: %.1e\n",
        reference_depth, ref_seconds, reference_depth - 2, deviation(check))
    if exact
        floor = maximum(zip(saves, reference)) do (t, W)
            μ, Σ = cold_memory_reference(t, mean0, covariance0, bath; mass = 1.0, omega = 1.0)
            maximum(abs, W - on_grid((q, p) -> gaussian_wigner(q, p; mean = μ, covariance = Σ), grid))
        end / scale
        @printf("reference vs exact Gaussian (grid floor): %.1e\n", floor)
    end
    flush(stdout)

    @printf("%-12s %5s %7s %6s %9s %8s %8s %9s\n",
        "closure", "depth", "members", "solves", "rhs (ms)", "nf", "wall (s)", "error")
    rows = []
    for depth in depths
        op = build(depth)
        U = randn(size(grid)..., length(op.indices))
        for kind in CLOSURES
            cl = terminate(op, kind)
            Ws, seconds, nf = propagate(cl, W0, saves)
            err = Ws === nothing ? Inf : deviation(Ws)
            ms = 1e3 * rhs_seconds(cl, U)
            @printf("%-12s %5d %7d %6d %9.3f %8d %8.2f %9.2e\n",
                kind, depth, length(op.indices), solves_per_rhs(cl), ms, nf, seconds, err)
            flush(stdout)
            push!(rows, (; case = name, kind, depth, members = length(op.indices),
                solves = solves_per_rhs(cl), rhs_ms = ms, nf, seconds, error = err))
        end
    end
    return rows
end

harmonic = harmonic_potential(; mass = 1.0, omega = 1.0)
quartic(q) = q^2 / 2 + 0.08q^4
cold_strong = drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 2)
CASES = Dict(
    "warm" => () -> run_case("harmonic, warm (λ=0.2, γ=1, kT=1), Padé 1";
        bath = drude_lorentz_pade_bath(; reorganization = 0.2, cutoff = 1.0, kT = 1.0, pade = 1),
        potential = harmonic, moyal_terms = 1, depths = 0:6, reference_depth = 12,
        exact = true, saves = 0.5:0.5:10.0),
    "coldweak" => () -> run_case("harmonic, cold weak (λ=0.05, γ=0.5, kT=0.1), Padé 2";
        bath = drude_lorentz_pade_bath(; reorganization = 0.05, cutoff = 0.5, kT = 0.1, pade = 2),
        potential = harmonic, moyal_terms = 1, depths = 0:5, reference_depth = 10,
        exact = true, saves = 0.5:0.5:10.0),
    "coldstrong" => () -> run_case("harmonic, cold strong (λ=0.8, γ=0.5, kT=0.1), Padé 2";
        bath = cold_strong, potential = harmonic, moyal_terms = 1, depths = 0:6,
        reference_depth = 10, exact = true, saves = 0.25:0.25:2.5),
    "quartic" => () -> run_case("quartic q²/2+0.08q⁴, cold strong, Padé 2";
        bath = cold_strong, potential = quartic, moyal_terms = 2, depths = 0:6,
        reference_depth = 10, exact = false, saves = 0.25:0.25:2.5),
)
selected = isempty(ARGS) ? ["warm", "coldweak", "coldstrong", "quartic"] : ARGS
rows = reduce(vcat, [CASES[c]() for c in selected])
open(joinpath(@__DIR__, "benchmark_$(join(selected, "_")).csv"), "w") do io
    println(io, join(keys(first(rows)), ","))
    for r in rows
        println(io, join(("\"$(r.case)\"", Base.tail(values(r))...), ","))
    end
end
