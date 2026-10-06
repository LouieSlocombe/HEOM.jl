@testset "Time-dependent generators" begin
    grid = PhaseSpaceGrid((-5.0, 5.0), 16, (-5.0, 5.0), 18)
    W = on_grid((q, p) -> exp(-q^2 - p^2) / π, grid)
    V0(q) = q^2 / 2 + 0.02q^4
    μ(q) = q + 0.03q^2
    E(t) = 0.2sin(1.3t)
    dE(t) = 0.26cos(1.3t)
    potential = DrivenPotential(V0, E, μ; field_derivative = dE)
    general = TimeDependentPotential((q, t) -> V0(q) - E(t) * μ(q))
    @test potential(0.3, 0.7) == general(0.3, 0.7)
    @test HEOM.is_time_dependent(potential)
    @test !HEOM.is_time_dependent(V0)
    # A static generator has zero partial time derivative even for a nonzero state.
    static = wigner_moyal_operator(grid; mass = 1.2, potential = V0)
    static_tgrad = fill(NaN, size(W))
    HEOM.wigner_tgrad!(static_tgrad, W, static, 0.7)
    @test all(iszero, static_tgrad)

    for (discretization, terms) in
        ((Spectral(), nothing), (Spectral(), 2), (FiniteDifference(4), 2))
        for V in (potential, general, (q, t) -> V0(q) - E(t) * μ(q))
            op = wigner_moyal_operator(
                grid;
                mass = 1.2,
                potential = V,
                discretization,
                moyal_terms = terms,
            )
            @test HEOM.is_time_dependent(op)
            # Solvers may revisit earlier times during rejected steps or interpolation.
            for t in (0.7, 1.1, 0.0, 0.7)
                frozen = wigner_moyal_operator(
                    grid;
                    mass = 1.2,
                    potential = q -> V0(q) - E(t) * μ(q),
                    discretization,
                    moyal_terms = terms,
                )
                actual, expected = similar(W), similar(W)
                wigner_moyal!(actual, W, op, t)
                wigner_moyal!(expected, W, frozen, t)
                @test actual ≈ expected atol = 2e-13
                dt, plus, minus = similar(W), similar(W), similar(W)
                HEOM.wigner_tgrad!(dt, W, op, t)
                wigner_moyal!(plus, W, op, t + 1e-5)
                wigner_moyal!(minus, W, op, t - 1e-5)
                @test dt ≈ (plus - minus) / 2e-5 atol = 2e-10
                if discretization isa FiniteDifference
                    @test_throws ArgumentError sparse(op)
                    @test sparse(op, t) * vec(W) ≈ vec(expected) atol = 2e-13
                end
            end
        end
    end

    @testset "Spatial symbols are cached for separable drives" begin
        calls = Ref(0)
        dipole(q) = (calls[] += 1; q)
        op = wigner_moyal_operator(
            grid;
            mass = 1.0,
            potential = DrivenPotential(V0, E, dipole),
        )
        count = calls[]
        wigner_moyal!(similar(W), W, op, 0.2)
        wigner_moyal!(similar(W), W, op, 0.4)
        @test calls[] == count
    end

    @testset "Analytic time derivatives and time domains" begin
        # Float64-only time inputs demonstrate no time AD when callbacks are supplied.
        field(t::Float64) = sin(t)
        value(q, t::Float64) = q^2 / 2 - field(t) * q
        for V in (
            DrivenPotential(q -> q^2 / 2, field, q -> q; field_derivative = cos),
            TimeDependentPotential(value; derivative = (q, t) -> -cos(t) * q),
        )
            op = wigner_moyal_operator(grid; mass = 1.0, potential = V, moyal_terms = 2)
            out = similar(W)
            HEOM.wigner_tgrad!(out, W, op, 0.3)
            @test all(isfinite, out)
            @test maximum(abs, out) > 0
        end
        op = wigner_moyal_operator(
            grid;
            mass = 1.0,
            potential = (q, t) -> q^2 / 2 - log(t) * q,
        )
        wigner_moyal!(similar(W), W, op, 2.0)
        bad = wigner_moyal_operator(
            grid;
            mass = 1.0,
            potential = DrivenPotential(V0, t -> NaN, μ),
        )
        @test_throws ArgumentError wigner_moyal!(similar(W), W, bad, 0.0)
        bad = wigner_moyal_operator(grid; mass = 1.0, potential = (q, t) -> (1 + im) * q)
        @test_throws ArgumentError wigner_moyal!(similar(W), W, bad, 0.0)
    end

    @testset "Caldeira–Leggett drive" begin
        for discretization in (Spectral(), FiniteDifference(4))
            op = caldeira_leggett_operator(
                grid;
                mass = 1.2,
                potential,
                friction = 0.15,
                kT = 1.0,
                discretization,
                moyal_terms = 2,
            )
            t = 0.43
            frozen = caldeira_leggett_operator(
                grid;
                mass = 1.2,
                potential = q -> potential(q, t),
                friction = 0.15,
                kT = 1.0,
                discretization,
                moyal_terms = 2,
            )
            actual, expected = similar(W), similar(W)
            caldeira_leggett!(actual, W, op, t)
            caldeira_leggett!(expected, W, frozen, t)
            @test actual ≈ expected atol = 2e-13
            if discretization isa FiniteDifference
                @test_throws ArgumentError sparse(op)
                @test sparse(op, t) * vec(W) ≈ vec(expected) atol = 2e-13
            end
        end
    end

    @testset "HEOM time derivatives and sparse Jacobians" begin
        bath = ExponentialBath([0.2 - 0.1im], [1.5]; counterterm = 0.07, diffusion = 0.02)
        for discretization in (Spectral(), FiniteDifference(4)), V in (potential, general)
            op = heom_operator(
                grid;
                mass = 1.2,
                potential = V,
                bath,
                depth = 1,
                discretization,
                moyal_terms = 2,
            )
            U = cat(W, 0.3W; dims = 3)
            t = 0.42
            frozen = heom_operator(
                grid;
                mass = 1.2,
                potential = q -> potential(q, t),
                bath,
                depth = 1,
                discretization,
                moyal_terms = 2,
            )
            actual, expected = similar(U), similar(U)
            heom!(actual, U, op, t)
            heom!(expected, U, frozen, t)
            @test actual ≈ expected atol = 2e-13
            prob = heom_problem(U, (0.0, 1.0), op)
            jv = similar(U)
            prob.f.jvp(jv, U, U, op, t)
            @test jv ≈ expected atol = 2e-13
            flat = similar(vec(U))
            prob.f.jvp(flat, vec(U), vec(U), op, t)
            @test flat ≈ vec(expected) atol = 2e-13
            plus, minus, dt = similar(U), similar(U), similar(U)
            prob.f.tgrad(dt, U, op, t)
            heom!(plus, U, op, t + 1e-5)
            heom!(minus, U, op, t - 1e-5)
            @test dt ≈ (plus - minus) / 2e-5 atol = 2e-10
            if discretization isa FiniteDifference
                @test_throws ArgumentError sparse(op)
                @test sparse(op, t) * vec(U) ≈ vec(expected) atol = 2e-13
                sparse_prob = heom_problem(U, (0.0, 1.0), op; jacobian = :sparse)
                J = copy(sparse_prob.f.jac_prototype)
                columns, rows = copy(J.colptr), copy(J.rowval)
                sparse_prob.f.jac(J, U, op, t)
                @test J * vec(U) ≈ vec(expected) atol = 2e-13
                @test J.colptr == columns
                @test J.rowval == rows
            end
        end
    end

    @testset "Driven stiff integration" begin
        small = PhaseSpaceGrid((-5.0, 5.0), 12, (-5.0, 5.0), 12)
        W0 = on_grid((q, p) -> exp(-q^2 - p^2) / π, small)
        bath = ExponentialBath([0.15 - 0.05im], [1.3])
        for V in (
            DrivenPotential(q -> q^2 / 2, t -> sin(t), q -> 0.03q^4),
            TimeDependentPotential((q, t) -> q^2 / 2 - 0.03sin(t) * q^4),
        )
            op = heom_operator(
                small;
                mass = 1.0,
                potential = V,
                bath,
                depth = 1,
                discretization = FiniteDifference(2),
                moyal_terms = 2,
            )
            problem = heom_problem(W0, (0.0, 0.3), op; jacobian = :sparse)
            implicit = solve(
                problem,
                Rodas5P();
                abstol = 1e-9,
                reltol = 1e-9,
                save_everystep = false,
            )
            explicit = solve(
                problem,
                Vern9();
                abstol = 1e-11,
                reltol = 1e-11,
                save_everystep = false,
            )
            @test implicit.u[end] ≈ explicit.u[end] atol = 2e-8
        end
        for V in (
            DrivenPotential(q -> q^2 / 2, t -> sin(t), q -> q),
            TimeDependentPotential((q, t) -> q^2 / 2 - sin(t) * q),
        )
            op = heom_operator(small; mass = 1.0, potential = V, bath, depth = 1)
            problem = heom_problem(W0, (0.0, 0.3), op)
            algorithm = Rodas5P(
                autodiff = AutoFiniteDiff(),
                linsolve = KrylovJL_GMRES(),
                concrete_jac = false,
            )
            implicit = solve(
                problem,
                algorithm;
                abstol = 1e-8,
                reltol = 1e-8,
                save_everystep = false,
            )
            explicit = solve(
                problem,
                Vern9();
                abstol = 1e-11,
                reltol = 1e-11,
                save_everystep = false,
            )
            @test implicit.u[end] ≈ explicit.u[end] atol = 2e-7
        end
    end
end
