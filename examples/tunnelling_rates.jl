# Fit transfer rates from the population dynamics of a bath-coupled double well.
# Requires HEOM and OrdinaryDiffEqVerner in the active Julia environment.
# Run: julia --project=/path/to/environment examples/tunnelling_rates.jl
# This modest calculation illustrates the workflow; it is not a converged benchmark.
using HEOM, OrdinaryDiffEqVerner, Printf

function double_well_trajectory(; points = 48, depth = 8, matsubara = 0)
    mass, hbar = 1.0, 0.5
    # Minima at ±1.5, barrier height 1, and local harmonic frequency sqrt(8/2.25).
    potential(q) = (q^2 / 2.25 - 1)^2
    grid = PhaseSpaceGrid((-5.0, 5.0), points, (-6.0, 6.0), points)
    initial = on_grid(
        (q, p) -> coherent_wigner(
            q,
            p;
            q0 = -1.5,
            p0 = 0.0,
            mass,
            omega = sqrt(8 / 2.25),
            hbar,
        ),
        grid,
    )
    bath =
        drude_lorentz_bath(; reorganization = 0.2, cutoff = 2.0, kT = 1.0, hbar, matsubara)
    op = heom_operator(grid; mass, potential, bath, depth, scaled = true, moyal_terms = 2)
    # For this quartic potential the Moyal series terminates at the hbar² term.
    # The operator includes the bath counterterm. Matrix initial data gives zero
    # ADOs and initial bath slip; the fitting windows exclude the preparation.
    sol = solve(
        heom_problem(initial, (0.0, 12.0), op),
        Vern7();
        saveat = 0.5,
        abstol = 1e-8,
        reltol = 1e-8,
    )
    @assert last(sol.t) == 12.0 "The solver did not reach the final time"
    @assert all(U -> all(isfinite, U), sol.u)
    return (; sol, grid, potential)
end

function report_transfer_rates(data)
    (; sol, grid, potential) = data
    # Reflection symmetry fixes this independently of the observed final sample,
    # even for the correlated equilibrium of the coordinate-coupled bath.
    equilibrium_product = 0.5
    println("Product: q > 0; equilibrium population = 0.5 by symmetry")
    println("Window           k_forward     k_backward    relaxation     log R²       RMS")
    for window in ((4.0, 9.0), (5.0, 10.0), (6.0, 12.0))
        rates = tunnelling_rates(sol; equilibrium_product, tspan = window)
        @printf(
            "[%4.1f, %4.1f]   %.6f      %.6f      %.6f      %.6f   %.3e\n",
            window...,
            rates.forward_rate,
            rates.backward_rate,
            rates.relaxation_rate,
            rates.r_squared,
            rates.rmse,
        )
    end
    d = diagnostics(sol; potential)
    @printf("Maximum norm error: %.3e\n", maximum(abs, d.norm .- 1))
    @printf("Maximum position boundary weight: %.3e\n", maximum(d.boundary_q))
    @printf("Maximum momentum boundary weight: %.3e\n", maximum(d.boundary_p))
    @printf("Maximum position spectral tail: %.3e\n", maximum(d.tail_q))
    @printf("Maximum momentum spectral tail: %.3e\n", maximum(d.tail_p))
    @printf(
        "Final product population: %.6f (not used as equilibrium)\n",
        probability(physical_wigner(sol), grid; q = (0.0, Inf)),
    )
    println("Rates have inverse simulation-time units. Inspect the window dependence.")
    println("Converge depth, bath poles, grid resolution and box size independently.")
    println("This fits total interwell transfer, including thermally activated transfer.")
    println("It does not establish a separate under-barrier tunnelling contribution.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    report_transfer_rates(double_well_trajectory())
end
