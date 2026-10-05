# Cold, strongly coupled ANHARMONIC dynamics with separate convergence controls.
# Run with HEOM and OrdinaryDiffEqVerner in the active Julia environment.
# This example takes longer than heom_sho.jl: it solves seven full hierarchies.
using HEOM, OrdinaryDiffEqVerner, LinearAlgebra

# kT/(ħω)=0.1. The zero-frequency Drude friction is 2λ/(mγ)=3.2ω;
# neither a high-temperature nor a weak-coupling master equation is appropriate.
mass, omega, hbar = 1.0, 1.0, 1.0
reorganization, cutoff, kT = 0.8, 0.5, 0.1
potential(q) = mass * omega^2 * q^2 / 2 + 0.08q^4

function cold_anharmonic_run(; points = 40, extent = 7.0, pade = 2, depth = 6)
    grid = PhaseSpaceGrid((-extent, extent), points, (-extent, extent), points)
    bath = drude_lorentz_pade_bath(; reorganization, cutoff, kT, hbar, pade)
    initial = on_grid(
        (q, p) -> coherent_wigner(q, p; q0 = 0.7, p0 = -0.2, mass, omega, hbar),
        grid,
    )
    op = heom_operator(grid; mass, potential, bath, depth, scaled = true)
    sol = solve(
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
    @assert all(isfinite, state)
    println(
        "N=",
        pade,
        ", depth=",
        depth,
        ", grid=",
        points,
        "², extent=±",
        extent,
        ": members=",
        length(hierarchy_indices(op)),
        ", energy=",
        energy(state, grid; mass, potential),
        ", negativity=",
        wigner_negativity(state, grid),
    )
    return state
end

# Each change varies one numerical approximation only. The compared quantities
# are the whole Wigner distributions, not just moments (which can hide tier error).
coarse_depth = cold_anharmonic_run(; depth = 4)
base = cold_anharmonic_run()
fine_depth = cold_anharmonic_run(; depth = 8)
poles3 = cold_anharmonic_run(; pade = 3)
poles4 = cold_anharmonic_run(; pade = 4)

# Resolution changes at fixed box; every second fine-grid point matches the base.
fine_grid = cold_anharmonic_run(; points = 80)
# Box changes at fixed spacing 0.35; cropping aligns with the base coordinates.
wide_box = cold_anharmonic_run(; points = 48, extent = 8.4)

changes = (
    depth_4_to_6 = maximum(abs, coarse_depth - base),
    depth_6_to_8 = maximum(abs, base - fine_depth),
    pade_2_to_3 = maximum(abs, base - poles3),
    pade_3_to_4 = maximum(abs, poles3 - poles4),
    resolution_40_to_80 = maximum(abs, base - fine_grid[1:2:end, 1:2:end]),
    extent_7_to_8p4 = maximum(abs, base - wide_box[5:44, 5:44]),
)
println("Maximum pointwise changes in the final Wigner function:")
for (control, change) in pairs(changes)
    println("  ", control, ": ", change)
end

@assert changes.depth_6_to_8 < changes.depth_4_to_6
@assert changes.pade_3_to_4 < changes.pade_2_to_3
@assert changes.depth_6_to_8 < 2e-6
@assert changes.pade_3_to_4 < 2e-3
@assert changes.resolution_40_to_80 < 2e-4
@assert changes.extent_7_to_8p4 < 2e-4

# At these modest orders, the bath approximation dominates (~10⁻³ pointwise).
# For a tighter requested accuracy, increase pade further and repeat ALL controls
# around that refined calculation. These thresholds establish this finite-time
# example, not equilibrium or an error bound for different temperatures/potentials.
# Longer trajectories also require a longer convergence study and may benefit from
# a stiff solver with an exact Jacobian-vector product (see heom_problem).
# A 2D initial Wigner matrix prepares a factorized system/bath state. Relaxing
# the FULL hierarchy and restarting it prepares a correlated reduced equilibrium;
# retaining only its root and resetting the ADOs would reintroduce initial slip.
