using SciMLBase: DiscreteCallback, ReturnCode, terminate!

@testset "Equilibrium stationarity covers the full hierarchy" begin
    grid = PhaseSpaceGrid((-2, 2), 8, (-2, 2), 8)
    bath = ExponentialBath([0.0], [1.0])
    op = heom_operator(grid; mass = 1, potential = q -> zero(q), bath, depth = 2)
    W0 = fill(1 / 8, size(grid))
    U0 = cat(W0, fill(0.2, size(grid)), fill(-0.1, size(grid)); dims = 3)
    original = copy(U0)
    atol, rtol = 1e-8, 1e-6

    # These constant auxiliaries obey exactly dU_n/dt = -n U_n. The physical
    # member is stationary from the start, but the hierarchy is not.
    result = equilibrate(
        U0,
        (2.0, 2.7),
        op,
        Vern7();
        stationarity_abstol = atol,
        stationarity_reltol = rtol,
        check_interval = 0.2,
        abstol = 1e-12,
        reltol = 1e-12,
        saveat = [2.1, 2.4],
    )
    @test result isa EquilibriumResult
    @test !result.converged
    @test result.status === :time_limit
    @test result.retcode === ReturnCode.Success
    @test result.restart.time == 2.7
    @test U0 == original
    @test result.hierarchy !== U0
    @test size(result.hierarchy) == size(U0)
    @test physical_wigner(result) == W0
    @test parent(physical_wigner(result)) === result.hierarchy
    @test phase_space_integral(physical_wigner(result), grid) == 2

    amplitudes = [1 / 8, 0.2exp(-0.7), 0.1exp(-1.4)]
    absolute = [0.0, amplitudes[2], 2amplitudes[3]]
    @test result.residuals.absolute ≈ absolute atol = 1e-11
    @test result.residuals.relative ≈ [0.0, 1.0, 2.0] atol = 1e-12
    @test result.residuals.scaled ≈ absolute ./ (atol .+ rtol .* amplitudes) rtol = 1e-9
    @test result.residuals.scaled[1] <= 1
    @test all(>(1), result.residuals.scaled[2:end])
    @test result.restart.operator === op
    @test result.restart.indices == hierarchy_indices(op)
    @test result.restart.scaled === false
    @test result.restart.depth == 2
    @test result.restart.jacobian === :matrixfree
    @test result.restart.tspan == (2.0, 2.7)
    @test result.restart.norm == result.restart.initial_norm == 2

    # Checking just the physical and first-tier members still misses an evolving
    # higher auxiliary. Here those first two derivatives vanish identically.
    high_tier = copy(U0)
    high_tier[:, :, 2] .= 0
    high_result = equilibrate(
        high_tier,
        (0.0, 0.1),
        op,
        Vern7();
        stationarity_abstol = atol,
        stationarity_reltol = rtol,
        abstol = 1e-12,
        reltol = 1e-12,
    )
    @test !high_result.converged
    @test high_result.status === :time_limit
    @test high_result.residuals.absolute[1:2] == [0.0, 0.0]
    @test high_result.residuals.relative[1:2] == [0.0, 0.0]
    @test high_result.residuals.scaled[3] > 1

    # Restarting preserves the auxiliaries, their convention and their ordering.
    restarted = heom_problem(result, (result.restart.time, 3.1))
    @test restarted.p === op
    @test restarted.u0 == result.hierarchy
    @test restarted.u0 !== result.hierarchy
    resumed = solve(restarted, Vern7(); abstol = 1e-12, reltol = 1e-12)
    uninterrupted =
        solve(heom_problem(U0, (2.0, 3.1), op), Vern7(); abstol = 1e-12, reltol = 1e-12)
    @test resumed.u[end] ≈ uninterrupted.u[end] atol = 1e-11
    restarted.u0[1, 1, 2] += 1
    @test restarted.u0[1, 1, 2] != result.hierarchy[1, 1, 2]
    result.restart.indices[1][1] = 99
    @test first(hierarchy_indices(op)) == [0]
end

@testset "Equilibrium detection and correlated initialization" begin
    grid = PhaseSpaceGrid((-2, 2), 8, (-2, 2), 8)
    bath = ExponentialBath([0.0], [1.0])
    W0 = fill(1 / 8, size(grid))
    for scaled in (false, true)
        op =
            heom_operator(grid; mass = 1, potential = q -> zero(q), bath, depth = 2, scaled)
        # Exact zero auxiliary amplitudes and residuals have relative residual 0.
        ready = equilibrate(W0, (3.0, 30.0), op, Vern7())
        @test ready.converged
        @test ready.status === :stationary
        @test ready.restart.time == 3.0
        @test ready.restart.scaled === scaled
        @test ready.residuals.absolute == zeros(3)
        @test ready.residuals.relative == zeros(3)
        @test ready.residuals.scaled == zeros(3)
        @test all(iszero, ready.hierarchy[:, :, 2:end])
        @test physical_wigner(ready) == W0
        physical_wigner(ready)[1, 1] += 1
        @test W0[1, 1] == 1 / 8

        U0 = cat(W0, fill(0.2, size(grid)), fill(-0.1, size(grid)); dims = 3)
        prepared = equilibrate(
            U0,
            (3.0, 30.0),
            op,
            Vern7();
            stationarity_abstol = 1e-5,
            stationarity_reltol = 0,
            check_interval = 0.5,
            abstol = 1e-12,
            reltol = 1e-12,
        )
        @test prepared.converged
        @test prepared.status === :stationary
        @test 3.0 < prepared.restart.time < 30.0
        @test all(<=(1), prepared.residuals.scaled)
        @test physical_wigner(prepared) == W0
        @test phase_space_integral(physical_wigner(prepared), grid) == 2
        elapsed = prepared.restart.time - 3.0
        @test prepared.hierarchy[:, :, 2] ≈ fill(0.2exp(-elapsed), size(grid)) atol = 1e-11
        @test prepared.hierarchy[:, :, 3] ≈ fill(-0.1exp(-2elapsed), size(grid)) atol =
            1e-11

        # The mandatory final check can certify convergence even if no scheduled
        # stationarity check falls inside this preparation interval.
        final_check = equilibrate(
            U0,
            (0.0, 10.0),
            op,
            Vern7();
            stationarity_abstol = 1e-5,
            stationarity_reltol = 0,
            check_interval = 20,
            abstol = 1e-12,
            reltol = 1e-12,
        )
        @test final_check.converged
        @test final_check.status === :stationary
        @test final_check.retcode === ReturnCode.Success
        @test final_check.restart.time == 10
    end

    # A factorized state can have a stationary physical root even while the bath
    # starts building correlations. It must not pass an initial root-only check.
    coupled = heom_operator(
        grid;
        mass = 1,
        potential = q -> zero(q),
        bath = ExponentialBath([0.4 - 0.2im], [1.0]),
        depth = 2,
    )
    building = equilibrate(
        W0,
        (0.0, 1e-9),
        coupled,
        Vern7();
        stationarity_abstol = 1e-7,
        stationarity_reltol = 0,
        abstol = 1e-12,
        reltol = 1e-12,
    )
    @test !building.converged
    @test building.status === :time_limit
    @test building.residuals.scaled[1] <= 1
    @test building.residuals.scaled[2] > 1
    @test maximum(abs, building.hierarchy[:, :, 2]) > 0

    assessed = equilibrate(W0, (2.0, 2.0), coupled, Vern7())
    @test !assessed.converged
    @test assessed.status === :time_limit
    @test assessed.restart.time == 2
    @test assessed.retcode === ReturnCode.Success
    @test assessed.residuals.absolute[1] == 0
    @test assessed.residuals.relative[2] === Inf
    @test assessed.residuals.absolute[2] > 0
    @test all(iszero, assessed.hierarchy[:, :, 2:end])

    finite_difference = heom_operator(
        grid;
        mass = 1,
        potential = q -> zero(q),
        bath,
        depth = 2,
        discretization = FiniteDifference(4),
        moyal_terms = 1,
    )
    sparse_result =
        equilibrate(W0, (0.0, 0.0), finite_difference, Vern7(); jacobian = :sparse)
    @test sparse_result.converged
    @test sparse_result.restart.jacobian === :sparse
    @test heom_problem(sparse_result, (0.0, 1.0)).f.jac_prototype !== nothing
    @test heom_problem(sparse_result, (0.0, 1.0); jacobian = :matrixfree).f.jac_prototype ===
          nothing

    # Finite input can still overflow the derivative evaluation. Assessment must
    # report unusable residuals rather than certify this hierarchy as stationary.
    huge = cat(W0, fill(floatmax(Float64), size(grid)); dims = 3)
    overflow_op = heom_operator(grid; mass = 1, potential = q -> zero(q), bath, depth = 1)
    overflow = equilibrate(huge, (0.0, 0.0), overflow_op, Vern7())
    @test !overflow.converged
    @test overflow.residuals.absolute[2] === Inf
    @test overflow.residuals.relative[2] === Inf
    @test overflow.residuals.scaled[2] === Inf
end

@testset "Equilibrium solver exits and validation" begin
    grid = PhaseSpaceGrid((-2, 2), 8, (-2, 2), 8)
    bath = ExponentialBath([0.0], [1.0])
    op = heom_operator(grid; mass = 1, potential = q -> zero(q), bath, depth = 1)
    W0 = fill(1 / 8, size(grid))
    U0 = cat(W0, fill(0.2, size(grid)); dims = 3)
    failed = equilibrate(
        U0,
        (0.0, 1.0),
        op,
        Vern7();
        stationarity_abstol = 1e-12,
        stationarity_reltol = 0,
        dt = 0.01,
        adaptive = false,
        maxiters = 1,
    )
    @test !failed.converged
    @test failed.status === :solver_failure
    @test failed.retcode === ReturnCode.MaxIters
    @test 0 <= failed.restart.time < 1
    @test failed.residuals.scaled[2] > 1
    @test failed.hierarchy[:, :, 2] ≈ fill(0.2exp(-failed.restart.time), size(grid)) atol =
        1e-11

    stop_early = DiscreteCallback(
        (u, t, integrator) -> t >= 0.2,
        terminate!;
        save_positions = (false, false),
    )
    stopped = equilibrate(
        U0,
        (0.0, 1.0),
        op,
        Vern7();
        stationarity_abstol = 1e-12,
        stationarity_reltol = 0,
        check_interval = 0.1,
        callback = stop_early,
        dt = 0.1,
        adaptive = false,
    )
    @test !stopped.converged
    @test stopped.status === :terminated
    @test stopped.retcode === ReturnCode.Terminated
    @test 0 < stopped.restart.time < 1
    @test stopped.residuals.scaled[2] > 1
    @test stopped.hierarchy[:, :, 2] ≈ fill(0.2exp(-stopped.restart.time), size(grid)) atol =
        1e-11

    for invalid in (-1.0, Inf, NaN)
        @test_throws ArgumentError equilibrate(
            W0,
            (0.0, 1.0),
            op,
            Vern7();
            stationarity_abstol = invalid,
        )
        @test_throws ArgumentError equilibrate(
            W0,
            (0.0, 1.0),
            op,
            Vern7();
            stationarity_reltol = invalid,
        )
    end
    @test_throws ArgumentError equilibrate(
        W0,
        (0.0, 1.0),
        op,
        Vern7();
        stationarity_abstol = 0,
    )
    for invalid in (0.0, -1.0, Inf, NaN)
        @test_throws ArgumentError equilibrate(
            W0,
            (0.0, 1.0),
            op,
            Vern7();
            check_interval = invalid,
        )
    end
    for invalid in ((1.0, 0.0), (0.0, Inf), (NaN, 1.0), (0.0,), (0.0, big"1e400"))
        @test_throws ArgumentError equilibrate(W0, invalid, op, Vern7())
    end
    for invalid in (zeros(size(grid)), -W0, fill(floatmax(Float64), size(grid)))
        @test_throws ArgumentError equilibrate(invalid, (0.0, 1.0), op, Vern7())
    end
    @test_throws ArgumentError equilibrate(complex.(W0), (0.0, 1.0), op, Vern7())
    @test_throws DimensionMismatch equilibrate(zeros(3, 4), (0.0, 1.0), op, Vern7())
    @test_throws DimensionMismatch equilibrate(zeros(8, 8, 3), (0.0, 1.0), op, Vern7())
    @test_throws ArgumentError equilibrate(W0, (0.0, 1.0), op, Vern7(); jacobian = :unknown)
    for options in ((; save_end = false), (; save_on = false), (; save_idxs = [1]))
        @test_throws ArgumentError equilibrate(W0, (0.0, 1.0), op, Vern7(); options...)
    end
    nonfinite = copy(U0)
    nonfinite[1, 1, 2] = NaN
    @test_throws ArgumentError equilibrate(nonfinite, (0.0, 1.0), op, Vern7())
end
