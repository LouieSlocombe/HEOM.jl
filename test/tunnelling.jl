using SciMLBase: DiscreteCallback, terminate!

@testset "Tunnelling rates from population relaxation" begin
    @testset "Independent reversible and irreversible master equations" begin
        times = [0.0, 0.13, 0.61, 1.4, 2.3, 4.8, 7.2]
        for (forward, backward) in ((0.21, 0.09), (0.05, 0.35), (0.3, 0.0), (0.0, 0.4))
            generator = [-forward backward; forward -backward]
            equilibrium_product = forward / (forward + backward)
            for initial_product in (0.1, 0.9)
                initial = [1 - initial_product, initial_product]
                populations = [(exp(t*generator)*initial)[2] for t in times]
                fit = tunnelling_rates(times, populations; equilibrium_product)
                @test fit.relaxation_rate ≈ forward + backward rtol = 1e-12
                @test fit.forward_rate ≈ forward atol = 1e-12
                @test fit.backward_rate ≈ backward atol = 1e-12
                @test fit.equilibrium_product == equilibrium_product
                @test fit.amplitude ≈ initial_product - equilibrium_product atol = 1e-12
                @test fit.r_squared ≈ 1 atol = 1e-14
                @test fit.rmse < 1e-13
                @test fit.times == times
                @test fit.population == populations
                @test fit.fitted_population ≈ populations atol = 1e-13
            end
        end

        populations = 0.6 .- 0.4 .* exp.(-0.3 .* times)
        shifted = tunnelling_rates(times .+ 17, populations; equilibrium_product = 0.6)
        rescaled = tunnelling_rates(1000times, populations; equilibrium_product = 0.6)
        @test shifted.relaxation_rate ≈ 0.3 rtol = 1e-12
        @test shifted.amplitude ≈ -0.4 atol = 1e-13
        @test rescaled.relaxation_rate ≈ 0.0003 rtol = 1e-12
        for unit_scale in (1e-300, 1e300)
            extreme = tunnelling_rates(
                unit_scale .* times,
                populations;
                equilibrium_product = 0.6,
            )
            @test extreme.relaxation_rate ≈ 0.3 / unit_scale rtol = 1e-12
            @test extreme.fitted_population ≈ populations atol = 1e-13
        end
        large_offset = tunnelling_rates(
            1e300 .+ 1e299 .* times,
            populations;
            equilibrium_product = 0.6,
        )
        @test large_offset.relaxation_rate ≈ 3e-300 rtol = 1e-12
    end

    @testset "Selected windows and fit diagnostics" begin
        times = collect(0:0.2:4)
        populations = 0.7 .- 0.5 .* exp.(-0.4 .* times)
        # A preparation transient before the kinetic regime must not enter the fit.
        populations[1:4] .= [0.10, 0.11, 0.15, 0.17]
        selected = 5:19
        fit = tunnelling_rates(
            times,
            populations;
            equilibrium_product = 0.7,
            tspan = (times[first(selected)], times[last(selected)]),
        )
        @test fit.times == times[selected]
        @test fit.population == populations[selected]
        @test fit.relaxation_rate ≈ 0.4 rtol = 1e-12
        @test fit.forward_rate ≈ 0.28 rtol = 1e-12
        @test fit.backward_rate ≈ 0.12 rtol = 1e-12
        @test fit.amplitude ≈ -0.5exp(-0.4times[first(selected)]) atol = 1e-13
        @test fit.fitted_population ≈ populations[selected] atol = 1e-13

        noisy_times = [0.0, 0.4, 1.1, 1.9, 3.0, 4.0]
        noisy = [0.1, 0.22, 0.28, 0.35, 0.44, 0.455]
        noisy_fit = tunnelling_rates(noisy_times, noisy; equilibrium_product = 0.5)
        @test 0 < noisy_fit.r_squared < 1
        @test noisy_fit.rmse > 0
        @test noisy_fit.rmse ≈
              sqrt(sum(abs2, noisy_fit.fitted_population - noisy) / length(noisy))
        logs = log.(abs.(noisy .- 0.5))
        fitted_logs = log.(abs.(noisy_fit.fitted_population .- 0.5))
        log_mean = sum(logs) / length(logs)
        @test noisy_fit.r_squared ≈
              1 - sum(abs2, logs - fitted_logs) / sum(abs2, logs .- log_mean)
        outside =
            tunnelling_rates(noisy_times, noisy; equilibrium_product = 0.5, tspan = (-1, 5))
        @test outside == noisy_fit
    end

    @testset "Invalid kinetic fits" begin
        times = [0.0, 1.0, 2.0, 3.0]
        populations = 0.6 .- 0.4 .* exp.(-0.3 .* times)
        @test_throws DimensionMismatch tunnelling_rates(
            times,
            populations[1:3];
            equilibrium_product = 0.6,
        )
        @test_throws ArgumentError tunnelling_rates(
            Float64[],
            Float64[];
            equilibrium_product = 0.6,
        )
        @test_throws ArgumentError tunnelling_rates(
            times[1:2],
            populations[1:2];
            equilibrium_product = 0.6,
        )
        for bad_times in
            ([0, 1, 1, 3], [3, 2, 1, 0], [0, 1, NaN, 3], [0, 1, Inf, 3], [0, 1, 2im, 3])
            @test_throws ArgumentError tunnelling_rates(
                bad_times,
                populations;
                equilibrium_product = 0.6,
            )
        end
        @test_throws ArgumentError tunnelling_rates(
            BigFloat[0, 1, big"1e400", big"2e400"],
            populations;
            equilibrium_product = 0.6,
        )
        @test_throws ArgumentError tunnelling_rates(
            [-1e308, -3e307, 3e307, 1e308],
            populations;
            equilibrium_product = 0.6,
        )
        for bad_populations in (
            [0.2, 0.3, NaN, 0.4],
            [0.2, 0.3, Inf, 0.4],
            [0.2, 0.3, 1im, 0.4],
            [-0.1, 0.2, 0.3, 0.4],
            [1.1, 0.9, 0.8, 0.7],
        )
            @test_throws ArgumentError tunnelling_rates(
                times,
                bad_populations;
                equilibrium_product = 0.6,
            )
        end
        @test_throws ArgumentError tunnelling_rates(
            times,
            BigFloat[0.2, 0.3, big"1e400", 0.4];
            equilibrium_product = 0.6,
        )
        for equilibrium_product in (-0.1, 1.1, NaN, Inf)
            @test_throws ArgumentError tunnelling_rates(
                times,
                populations;
                equilibrium_product,
            )
        end
        for tspan in ((1, 1), (2, 1), (0, Inf), (NaN, 2), (0,), (0, 1im), (4, 5), (0, 1))
            @test_throws ArgumentError tunnelling_rates(
                times,
                populations;
                equilibrium_product = 0.6,
                tspan,
            )
        end
        @test_throws ArgumentError tunnelling_rates(
            times,
            populations;
            equilibrium_product = 0.6,
            tspan = (big"1e400", big"2e400"),
        )
        @test_throws ArgumentError tunnelling_rates(
            times,
            populations;
            equilibrium_product = 0.6,
            min_deviation = 0,
        )
        @test_throws ArgumentError tunnelling_rates(
            [0, 1, 2],
            [1e308, 1e308, 1];
            equilibrium_product = 0,
            population_atol = 1e308,
        )
        for tolerance in (-1.0, NaN, Inf)
            @test_throws ArgumentError tunnelling_rates(
                times,
                populations;
                equilibrium_product = 0.6,
                min_deviation = tolerance,
            )
            @test_throws ArgumentError tunnelling_rates(
                times,
                populations;
                equilibrium_product = 0.6,
                population_atol = tolerance,
            )
        end
        # Equilibrium crossings, unresolved tails, growth and constant populations
        # do not establish a positive single-exponential relaxation rate.
        for bad_populations in (
            [0.2, 0.4, 0.6, 0.55],
            [0.2, 0.4, 0.49, 0.5],
            [0.2, 0.4, 0.49, 0.5 - 1e-10],
            [0.45, 0.4, 0.3, 0.2],
            fill(0.2, 4),
        )
            @test_throws ArgumentError tunnelling_rates(
                times,
                bad_populations;
                equilibrium_product = 0.5,
            )
        end
        # An allowed small overshoot is retained, never silently clipped to one.
        tolerated = 0.8 .+ (0.2 + 5e-9) .* exp.(-0.3 .* times)
        fit = tunnelling_rates(times, tolerated; equilibrium_product = 0.8)
        @test fit.population == tolerated
        @test fit.relaxation_rate ≈ 0.3 rtol = 1e-12
        @test_throws ArgumentError tunnelling_rates(
            times,
            tolerated;
            equilibrium_product = 0.8,
            population_atol = 1e-10,
        )
    end
end

@testset "Tunnelling rates from phase-space trajectories" begin
    # Smooth normalized densities whose half-space integrals are known analytically.
    grid = PhaseSpaceGrid((-π, π), 32, (-1, 1), 8)
    times = [0.0, 0.2, 0.7, 1.3, 2.0, 3.0]
    relaxation = 0.8
    equilibrium_coefficient, initial_coefficient = 0.25, -0.7
    density(coefficient) = on_grid((q, p) -> (1 + coefficient * sin(q)) / (4π), grid)
    coefficients =
        equilibrium_coefficient .+
        (initial_coefficient - equilibrium_coefficient) .* exp.(-relaxation .* times)
    states = density.(coefficients)
    wm = wigner_moyal_operator(grid; mass = 1, potential = q -> q^2 / 2)
    cl = caldeira_leggett_operator(
        grid;
        mass = 1,
        potential = q -> q^2 / 2,
        friction = 0.2,
        kT = 1,
    )
    bath = ExponentialBath([0.2 - 0.1im], [1.0])
    heom = heom_operator(grid; mass = 1, potential = q -> q^2 / 2, bath, depth = 1)
    # Auxiliary traces may be nonzero and must never enter a product population.
    hierarchies = [cat(W, fill(23.0, size(grid)); dims = 3) for W in states]

    @testset "Half-space orientation and physical hierarchy member" begin
        for (op, trajectory) in ((wm, states), (cl, states), (heom, hierarchies))
            for dividing_surface in (0.0, 0.31), product_side in (:right, :left)
                right_population(coefficient) =
                    (π - dividing_surface + coefficient * (1 + cos(dividing_surface))) /
                    (2π)
                expected = right_population.(coefficients)
                equilibrium_product = right_population(equilibrium_coefficient)
                if product_side === :left
                    expected = 1 .- expected
                    equilibrium_product = 1 - equilibrium_product
                end
                fit = tunnelling_rates(
                    trajectory,
                    op;
                    times,
                    dividing_surface,
                    product_side,
                    equilibrium_product,
                )
                @test fit.population ≈ expected atol = 1e-14
                @test fit.relaxation_rate ≈ relaxation rtol = 1e-12
                @test fit.forward_rate ≈ equilibrium_product * relaxation rtol = 1e-12
                @test fit.backward_rate ≈ (1 - equilibrium_product) * relaxation rtol =
                    1e-12
            end
        end
    end

    @testset "State validation and norm preservation" begin
        equilibrium_product = 0.5 + equilibrium_coefficient / π
        @test_throws DimensionMismatch tunnelling_rates(
            states,
            wm;
            times = times[1:3],
            equilibrium_product,
        )
        @test_throws ArgumentError tunnelling_rates(
            [2W for W in states],
            wm;
            times,
            equilibrium_product,
        )
        @test_throws DimensionMismatch tunnelling_rates(
            [ones(2, 2) for _ in times],
            wm;
            times,
            equilibrium_product,
        )
        @test_throws DimensionMismatch tunnelling_rates(
            [cat(W, W, W; dims = 3) for W in states],
            heom;
            times,
            equilibrium_product,
        )
        @test_throws ArgumentError tunnelling_rates(
            [reshape(W, size(W)..., 1) for W in states],
            wm;
            times,
            equilibrium_product,
        )
        for (op, trajectory) in ((wm, states), (heom, hierarchies))
            nonfinite = deepcopy(trajectory)
            nonfinite[2][1] = NaN
            @test_throws ArgumentError tunnelling_rates(
                nonfinite,
                op;
                times,
                equilibrium_product,
            )
            complex_states = [complex.(state) for state in trajectory]
            complex_states[2][1] += 1im
            @test_throws ArgumentError tunnelling_rates(
                complex_states,
                op;
                times,
                equilibrium_product,
            )
        end
        for product_side in (:both, :Left)
            @test_throws ArgumentError tunnelling_rates(
                states,
                wm;
                times,
                equilibrium_product,
                product_side,
            )
        end
        for dividing_surface in (-π - 0.1, -π, π, π + 0.1, Inf, NaN)
            @test_throws ArgumentError tunnelling_rates(
                states,
                wm;
                times,
                equilibrium_product,
                dividing_surface,
            )
        end
        for norm_atol in (-1.0, Inf, NaN)
            @test_throws ArgumentError tunnelling_rates(
                states,
                wm;
                times,
                equilibrium_product,
                norm_atol,
            )
        end
        driven =
            wigner_moyal_operator(grid; mass = 1, potential = (q, t) -> q^2 / 2 - t * q)
        @test_throws ArgumentError tunnelling_rates(
            states,
            driven;
            times,
            equilibrium_product,
        )

        scale = 1 + 2e-7
        scaled_states = [scale * W for W in states]
        fit = tunnelling_rates(
            scaled_states,
            wm;
            times,
            equilibrium_product = scale * equilibrium_product,
        )
        expected = scale .* (0.5 .+ coefficients ./ π)
        @test fit.population ≈ expected atol = 1e-14
        @test fit.relaxation_rate ≈ relaxation rtol = 1e-12
        @test_throws ArgumentError tunnelling_rates(
            scaled_states,
            wm;
            times,
            equilibrium_product,
            norm_atol = 1e-8,
        )
    end

    @testset "ODE solutions and incomplete integration" begin
        equilibrium_product = 0.5 + equilibrium_coefficient / π
        steady = density(equilibrium_coefficient)
        # Integrate a known linear relaxation independently of the fit, both as a
        # matrix and as a full hierarchy with deliberately large auxiliary traces.
        for (op, initial, target) in (
            (wm, first(states), steady),
            (heom, first(hierarchies), cat(steady, fill(-17.0, size(grid)); dims = 3)),
        )
            relaxation_rhs!(du, u, p, t) = (@. du = relaxation * (target - u))
            problem = ODEProblem(relaxation_rhs!, initial, (first(times), last(times)), op)
            solution =
                solve(problem, Vern7(); saveat = times, abstol = 1e-12, reltol = 1e-12)
            fit = tunnelling_rates(solution; equilibrium_product)
            @test fit.times == times
            @test fit.relaxation_rate ≈ relaxation rtol = 1e-9
            @test fit.population ≈ 0.5 .+ coefficients ./ π atol = 1e-10

            stop = DiscreteCallback((u, t, integrator) -> t >= 0.7, terminate!)
            partial = solve(problem, Vern7(); callback = stop, tstops = [0.7], saveat = 0.1)
            @test last(partial.t) < last(problem.tspan)
            @test_throws ArgumentError tunnelling_rates(partial; equilibrium_product)
            failed = solve(problem, Vern7(); maxiters = 1, dt = 0.01, adaptive = false)
            @test_throws ArgumentError tunnelling_rates(failed; equilibrium_product)
            unsaved_end =
                solve(problem, Vern7(); saveat = [0.0, 0.2, 0.7], save_end = false)
            @test_throws ArgumentError tunnelling_rates(unsaved_end; equilibrium_product)
        end

        scalar_problem = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 3.0))
        scalar_solution = solve(scalar_problem, Vern7(); saveat = times)
        @test_throws ArgumentError tunnelling_rates(
            scalar_solution;
            equilibrium_product = 0.5,
        )

        # Exercise a genuine HEOM solution as well. This short coherent transient
        # checks dispatch and population extraction, not a physical kinetic regime.
        physical_grid = PhaseSpaceGrid((-8, 8), 32, (-8, 8), 32)
        W0 = on_grid(
            (q, p) -> coherent_wigner(q, p; mass = 1, omega = 1, q0 = -0.5, p0 = 0.5),
            physical_grid,
        )
        heom_solution = solve(
            heom_problem(
                W0,
                (0, 0.15),
                physical_grid;
                mass = 1,
                potential = q -> q^2 / 2,
                bath,
                depth = 1,
            ),
            Vern7();
            saveat = 0.05,
            abstol = 1e-10,
            reltol = 1e-10,
        )
        physical_populations = [
            probability(physical_wigner(U), physical_grid; q = (0, Inf)) for
            U in heom_solution.u
        ]
        direct = tunnelling_rates(
            heom_solution.t,
            physical_populations;
            equilibrium_product = 0.5,
        )
        extracted = tunnelling_rates(heom_solution; equilibrium_product = 0.5)
        @test extracted.population ≈ direct.population atol = 1e-14
        @test extracted.forward_rate ≈ direct.forward_rate rtol = 1e-12
        @test extracted.backward_rate ≈ direct.backward_rate rtol = 1e-12
    end
end
