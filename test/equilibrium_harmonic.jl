@testset "Correlated harmonic preparation for a retained Drude bath" begin
    mass, omega, hbar, kT = 1.0, 1.0, 1.0, 1.0
    bath = drude_lorentz_bath(; reorganization = 0.1, cutoff = 2.0, kT, hbar, matsubara = 0)
    # The independent generalized-Langevin Lyapunov reference is exact for this
    # retained exponential plus Markovian tail, not for the continuum Drude bath.
    # cold_strong_benchmarks.jl separately checks bath convergence against FDT.
    _, _, joint_covariance = cold_memory_matrices(bath; mass, omega)
    covariance = joint_covariance[1:2, 1:2]
    grid = PhaseSpaceGrid((-8, 8), 40, (-8, 8), 40)
    bare_variance = hbar / (2tanh(hbar * omega / (2kT)))
    initial = on_grid(
        (q, p) ->
            gaussian_wigner(q, p; mean = zeros(2), covariance = bare_variance * I(2)),
        grid,
    )
    reference = on_grid((q, p) -> gaussian_wigner(q, p; mean = zeros(2), covariance), grid)
    op = heom_operator(
        grid;
        mass,
        potential = harmonic_potential(; mass, omega),
        bath,
        depth = 24,
        scaled = true,
    )
    prepared = equilibrate(
        initial,
        (0.0, 300.0),
        op,
        Vern7();
        stationarity_abstol = 1e-8,
        stationarity_reltol = 0,
        check_interval = 2.0,
        abstol = 1e-11,
        reltol = 1e-10,
    )
    @test prepared.converged
    @test prepared.status == :stationary
    @test 0 < prepared.restart.time < 300
    @test size(prepared.hierarchy, 3) == length(hierarchy_indices(op))
    @test all(prepared.residuals.scaled .<= 1)
    # Inspect every auxiliary independently of the result's convergence flag.
    derivative = similar(prepared.hierarchy)
    heom!(derivative, prepared.hierarchy, op, prepared.restart.time)
    for a in axes(derivative, 3)
        residual = maximum(abs, @view derivative[:, :, a])
        @test residual <= 1e-8
        @test prepared.residuals.absolute[a] ≈ residual
    end
    @test maximum(abs, @view prepared.hierarchy[:, :, 2:end]) > 0.01
    W = physical_wigner(prepared)
    @test phase_space_integral(W, grid) ≈ 1 atol = 2e-10
    @test phase_space_mean(W, grid) ≈ zeros(2) atol = 1e-7
    @test phase_space_covariance(W, grid) ≈ covariance atol = 5e-6
    @test max_error(W, reference) < 5e-7
    @test covariance[2, 2] - bare_variance > 0.02
    @test max_error(W, initial) > 1e-3

    # Reusing all auxiliaries preserves the equilibrium distribution and the bath
    # correlations. Keeping only its physical root reintroduces initial slip.
    duration = 0.5
    restart =
        heom_problem(prepared, (prepared.restart.time, prepared.restart.time + duration))
    continued =
        solve(restart, Vern7(); abstol = 1e-11, reltol = 1e-10, save_everystep = false)
    @test max_error(last(continued.u), prepared.hierarchy) < 1e-8
    @test phase_space_covariance(physical_wigner(continued), grid) ≈ covariance atol = 5e-6
    factorized = heom_problem(W, (0.0, duration), op)
    slipped =
        solve(factorized, Vern7(); abstol = 1e-11, reltol = 1e-10, save_everystep = false)
    @test max_error(physical_wigner(slipped), W) > 1e-3
end
