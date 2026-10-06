@testset "Grid-health diagnostics" begin
    @testset "Boundary weight" begin
        grid = PhaseSpaceGrid((-8, 8), 64, (-8, 8), 80)
        uniform = ones(size(grid))
        @test boundary_weight(uniform, grid) == (q = 6 / 64, p = 8 / 80)
        @test boundary_weight(-uniform, grid; width = 0.001) == (q = 2 / 64, p = 2 / 80)
        for n in (2, 3, 7)
            small = PhaseSpaceGrid((-1, 1), n, (-1, 1), n)
            expected = 2 * (n ÷ 2) / n
            @test boundary_weight(ones(n, n), small; width = prevfloat(0.5)) ==
                  (q = expected, p = expected)
        end
        W = on_grid((q, p) -> coherent_wigner(q, p; mass = 1, omega = 1), grid)
        boundary = @inferred boundary_weight(W, grid)
        @test boundary.q < 1e-12
        @test boundary.p < 1e-12
        displaced =
            on_grid((q, p) -> coherent_wigner(q, p; q0 = 6, mass = 1, omega = 1), grid)
        @test boundary_weight(displaced, grid).q > 1e-2
        @test_throws ArgumentError boundary_weight(W, grid; width = 0)
        @test_throws ArgumentError boundary_weight(W, grid; width = 0.5)
        @test_throws DimensionMismatch boundary_weight(ones(1, 80), grid)
    end

    @testset "Spectral tail" begin
        grid = PhaseSpaceGrid((0, 2π), 45, (0, 2π), 64)
        f = [cos(2π * 2i / 45) + 0.3cos(2π * 20i / 45) for i in 0:44]
        g = [1 + 0.2(-1)^j for j in 0:63]
        W = f * transpose(g)
        tail = @inferred spectral_tail(W, grid)
        @test tail.q ≈ 0.3 / sqrt(1.09) atol = 1e-14
        @test tail.p ≈ 0.2 / sqrt(1.04) atol = 1e-14
        rounding_grid = PhaseSpaceGrid((0, 2π), 20, (0, 2π), 20)
        for k in (2, 3)
            mode = on_grid((q, p) -> cos(k * q), rounding_grid)
            @test spectral_tail(mode, rounding_grid; fraction = 0.8).q ≈ k - 2 atol = 1e-14
        end
        tails = map((32, 128)) do n
            cat_grid = PhaseSpaceGrid((-8, 8), n, (-8, 8), n)
            cat = on_grid(
                (q, p) -> cat_wigner(q, p; q0 = 2, mass = 1.5, omega = 0.8, hbar = 0.7),
                cat_grid,
            )
            spectral_tail(cat, cat_grid).p
        end
        @test tails[1] > 0.5
        @test tails[2] < 1e-10
        @test_throws ArgumentError spectral_tail(W, grid; fraction = 0)
        @test_throws ArgumentError spectral_tail(W, grid; fraction = 1)
        @test_throws DimensionMismatch spectral_tail(ones(1, 64), grid)
    end
end

@testset "Driven diagnostic times" begin
    grid = PhaseSpaceGrid((-8, 8), 64, (-8, 8), 64)
    mass, omega = 1.0, 1.0
    W = on_grid((q, p) -> coherent_wigner(q, p; mass, omega, q0 = 0.4), grid)
    V0 = harmonic_potential(; mass, omega)
    potential = DrivenPotential(V0, sin, identity)
    times = [0.0, 0.1, 0.2]
    @test_throws ArgumentError diagnostics(W, grid; mass, potential)
    @test_throws ArgumentError diagnostics([W, W], grid; mass, potential)
    @test_throws DimensionMismatch diagnostics([W, W], grid; mass, potential, times)
    @test_throws ArgumentError diagnostics(
        [W, W],
        grid;
        mass,
        potential,
        times = [0.0, Inf],
    )
    @test diagnostics(W, grid; mass, potential, t = 0.2).energy ≈
          energy(W, grid; mass, potential, t = 0.2) atol = 1e-12

    bath = ExponentialBath([0.0], [1.0])
    for prob in (
        wigner_moyal_problem(W, (0.0, 0.2), grid; mass, potential),
        heom_problem(W, (0.0, 0.2), grid; mass, potential, bath, depth = 1),
    )
        sol = solve(prob, Vern7(); saveat = times, abstol = 1e-11, reltol = 1e-11)
        values = diagnostics(sol; potential)
        states = HEOM.physical_states(sol.u, sol.prob.p)
        expected =
            [energy(state, grid; mass, potential, t) for (state, t) in zip(states, times)]
        @test values.t == times
        @test values.energy ≈ expected atol = 1e-12
        @test values.energy == diagnostics(states, grid; mass, potential, times).energy
    end
end

@testset "Diagnostic summaries" begin
    m, ω, ħ = 1.5, 0.8, 0.7
    q0, p0 = 1.2, -0.8
    grid = PhaseSpaceGrid((-8, 8), 128, (-8, 8), 128)
    V = harmonic_potential(; mass = m, omega = ω)
    W = on_grid(
        (q, p) -> coherent_wigner(q, p; q0, p0, mass = m, omega = ω, hbar = ħ),
        grid,
    )
    fields = (
        :norm,
        :mean_q,
        :mean_p,
        :var_q,
        :var_p,
        :cov_qp,
        :uncertainty,
        :robertson_schrodinger,
        :energy,
        :purity,
        :negativity,
        :boundary_q,
        :boundary_p,
        :tail_q,
        :tail_p,
    )
    d = @inferred diagnostics(W, grid; mass = m, potential = V, hbar = ħ)
    @test keys(d) == fields
    @test all(x -> x isa Float64, d)
    @test d.norm ≈ 1 atol = 1e-12
    @test d.mean_q ≈ q0 atol = 1e-12
    @test d.mean_p ≈ p0 atol = 1e-12
    @test d.var_q ≈ ħ / (2m * ω) atol = 1e-12
    @test d.var_p ≈ ħ * m * ω / 2 atol = 1e-12
    @test d.cov_qp ≈ 0 atol = 1e-12
    @test d.uncertainty ≈ ħ / 2 atol = 1e-12
    @test d.robertson_schrodinger ≈ ħ / 2 atol = 1e-12
    @test d.energy ≈ ħ * ω / 2 + p0^2 / (2m) + V(q0) atol = 1e-12
    @test d.purity ≈ 1 atol = 1e-12
    @test d.negativity == 0
    @test (q = d.boundary_q, p = d.boundary_p) == boundary_weight(W, grid)
    @test (q = d.tail_q, p = d.tail_p) == spectral_tail(W, grid)
    @test diagnostics(0.9W, grid; mass = m, potential = V, hbar = ħ).norm ≈ 0.9 atol = 1e-12
    @test_throws DimensionMismatch diagnostics(ones(1, 128), grid; mass = m, potential = V)

    T = 2π / ω
    prob = wigner_moyal_problem(W, (0.0, T / 2), grid; mass = m, potential = V, hbar = ħ)
    times = [0.0, T / 4, T / 2]
    sol = solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12, saveat = times)
    trajectory = @inferred diagnostics(sol; potential = V)
    @test keys(trajectory) == (:t, fields..., :autocorrelation)
    @test trajectory.t == sol.t
    @test trajectory.norm ≈ ones(3) atol = 1e-11
    @test trajectory.purity ≈ ones(3) atol = 1e-11
    @test trajectory.energy ≈ fill(d.energy, 3) atol = 1e-11
    @test trajectory.uncertainty ≈ fill(ħ / 2, 3) atol = 1e-11
    @test trajectory.robertson_schrodinger ≈ fill(ħ / 2, 3) atol = 1e-11
    @test trajectory.mean_q ≈ q0 .* cos.(ω .* times) + p0 / (m * ω) .* sin.(ω .* times) atol =
        1e-10
    α² = (m * ω * q0^2 + p0^2 / (m * ω)) / (2ħ)
    @test trajectory.autocorrelation ≈ exp.(-2α² .* (1 .- cos.(ω .* times))) atol = 1e-11
    vectors = @inferred diagnostics(sol.u, grid; mass = m, potential = V, hbar = ħ)
    @test keys(vectors) == (fields..., :autocorrelation)
    @test vectors.mean_q == trajectory.mean_q
    @test vectors.autocorrelation == trajectory.autocorrelation
    @test_throws ArgumentError diagnostics(Matrix{Float64}[], grid; mass = m, potential = V)

    other_prob = ODEProblem((du, u, p, t) -> fill!(du, 0.0), W, (0.0, 1.0))
    other_sol = solve(other_prob, Vern9(); save_everystep = false)
    @test_throws ArgumentError diagnostics(other_sol; potential = V)
end
