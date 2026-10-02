@testset "Observables" begin
    @testset "Correlated Gaussian moments and energy" begin
        grid = PhaseSpaceGrid((-8, 8), 64, (-8, 8), 80)
        μ = [0.5, -1.0]
        Σ = [0.36 0.2; 0.2 0.81]
        W = on_grid((q, p) -> gaussian_wigner(q, p; mean = μ, covariance = Σ), grid)
        m, ω, ħ = 1.5, 0.8, 0.7
        V = harmonic_potential(; mass = m, omega = ω)
        q_density = exp.(-(grid.q .- μ[1]) .^ 2 / (2Σ[1, 1])) / sqrt(2π * Σ[1, 1])
        p_density = exp.(-(grid.p .- μ[2]) .^ 2 / (2Σ[2, 2])) / sqrt(2π * Σ[2, 2])

        @test phase_space_integral(W, grid) ≈ 1 atol = 1e-12
        @test (@inferred position_density(W, grid)) ≈ q_density atol = 1e-12
        @test (@inferred momentum_density(W, grid)) ≈ p_density atol = 1e-12
        @test (@inferred phase_space_mean(W, grid)) ≈ μ atol = 1e-12
        @test (@inferred phase_space_covariance(W, grid)) ≈ Σ atol = 1e-12
        @test expectation((q, p) -> q, W, grid) ≈ μ[1] atol = 1e-12
        @test expectation((q, p) -> p, W, grid) ≈ μ[2] atol = 1e-12
        @test expectation((q, p) -> (q - μ[1])^2, W, grid) ≈ Σ[1, 1] atol = 1e-12
        @test purity(W, grid; hbar = ħ) ≈ ħ / (2sqrt(det(Σ))) atol = 1e-12
        exact_energy = (Σ[2, 2] + μ[2]^2) / (2m) + m * ω^2 * (Σ[1, 1] + μ[1]^2) / 2
        @test (@inferred energy(W, grid; mass = m, potential = V)) ≈ exact_energy atol =
            1e-12
        quartic(q) = 0.1q^4 - 0.5q^2 + 0.1q
        @test energy(W, grid; mass = m, potential = quartic) ≈
              expectation((q, p) -> p^2 / (2m) + quartic(q), W, grid) atol = 1e-12

        # Norm drift remains visible in every moment and overlap.
        c = 0.8
        scaled = c * W
        @test phase_space_integral(scaled, grid) ≈ c atol = 1e-12
        @test position_density(scaled, grid) ≈ c * q_density atol = 1e-12
        @test momentum_density(scaled, grid) ≈ c * p_density atol = 1e-12
        @test phase_space_mean(scaled, grid) ≈ c * μ atol = 1e-12
        @test phase_space_covariance(scaled, grid) ≈ c * (Σ + (1 - c)^2 * μ * μ') atol =
            1e-12
        @test expectation((q, p) -> q, scaled, grid) ≈ c * μ[1] atol = 1e-12
        @test energy(scaled, grid; mass = m, potential = V) ≈ c * exact_energy atol = 1e-12
        @test overlap(scaled, W, grid; hbar = ħ) ≈ c * purity(W, grid; hbar = ħ) atol =
            1e-12
        @test purity(scaled, grid; hbar = ħ) ≈ c^2 * purity(W, grid; hbar = ħ) atol = 1e-12
    end

    @testset "Harmonic states" begin
        grid = PhaseSpaceGrid((-10, 10), 128, (-10, 10), 128)
        m, ω, ħ = 1.5, 0.8, 0.7
        for n in 0:3
            W = on_grid((q, p) -> fock_wigner(n, q, p; mass = m, omega = ω, hbar = ħ), grid)
            @test phase_space_mean(W, grid) ≈ zeros(2) atol = 1e-12
            @test phase_space_covariance(W, grid) ≈
                  [ħ / (m * ω) * (n + 1 / 2) 0; 0 ħ * m * ω * (n + 1 / 2)] atol = 1e-12
        end
        cat =
            on_grid((q, p) -> cat_wigner(q, p; q0 = 2, mass = m, omega = ω, hbar = ħ), grid)
        @test phase_space_covariance(cat, grid) ≈
              [4.287463427657553 0; 0 0.4139473358268759] atol = 1e-12

        r, θ = 0.3, π / 6
        R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        S = Diagonal([1 / sqrt(m * ω), sqrt(m * ω)])
        Σ0 = (ħ / 2) * S * R * Diagonal([exp(-2r), exp(2r)]) * R' * S
        μ0 = [0.5, -0.3]
        W0(q, p) = gaussian_wigner(q, p; mean = μ0, covariance = Σ0)
        for t in (0.0, 0.4, 1.1)
            W = on_grid(harmonic_evolution(W0, t; mass = m, omega = ω), grid)
            s, c = sincos(ω * t)
            M = [c s / (m * ω); -m * ω * s c]
            Σt = phase_space_covariance(W, grid)
            @test phase_space_mean(W, grid) ≈ M * μ0 atol = 1e-12
            @test Σt ≈ M * Σ0 * M' atol = 1e-12
            @test det(Σt) ≈ ħ^2 / 4 atol = 1e-12
        end
    end

    @testset "Overlaps" begin
        grid = PhaseSpaceGrid((-10, 10), 128, (-10, 10), 128)
        m, ω, ħ = 1.5, 0.8, 0.7
        q1, p1, q2, p2 = 1.2, -0.8, -0.4, 0.6
        W1 = on_grid(
            (q, p) ->
                coherent_wigner(q, p; q0 = q1, p0 = p1, mass = m, omega = ω, hbar = ħ),
            grid,
        )
        W2 = on_grid(
            (q, p) ->
                coherent_wigner(q, p; q0 = q2, p0 = p2, mass = m, omega = ω, hbar = ħ),
            grid,
        )
        ΔQ, ΔP = (q1 - q2) * sqrt(m * ω / ħ), (p1 - p2) / sqrt(m * ω * ħ)
        @test overlap(W1, W2, grid; hbar = ħ) ≈ exp(-(ΔQ^2 + ΔP^2) / 2) atol = 1e-12
        @test overlap(W1, W2, grid; hbar = ħ) == overlap(W2, W1, grid; hbar = ħ)
        @test purity(W1, grid; hbar = ħ) == overlap(W1, W1, grid; hbar = ħ)
        α² = (m * ω * q1^2 / ħ + p1^2 / (m * ω * ħ)) / 2
        for n in 0:5
            fock =
                on_grid((q, p) -> fock_wigner(n, q, p; mass = m, omega = ω, hbar = ħ), grid)
            @test overlap(W1, fock, grid; hbar = ħ) ≈ exp(-α²) * α²^n / factorial(n) atol =
                1e-12
        end
    end

    @testset "Negativity" begin
        grid = PhaseSpaceGrid((-8, 8), 256, (-8, 8), 256)
        ground = on_grid((q, p) -> fock_wigner(0, q, p; mass = 1, omega = 1), grid)
        excited = on_grid((q, p) -> fock_wigner(1, q, p; mass = 1, omega = 1), grid)
        @test wigner_negativity(ground, grid) === 0.0
        @test wigner_negativity((ground + excited) / 2, grid) === 0.0
        @test wigner_negativity(excited, grid) ≈ 4exp(-1 / 2) - 2 atol = 1e-3
        @test wigner_negativity(0.8excited, grid) ≈ 0.8wigner_negativity(excited, grid) atol =
            1e-12
    end

    @testset "Size checks" begin
        grid = PhaseSpaceGrid((-8, 8), 64, (-8, 8), 80)
        W, wrong = ones(size(grid)), ones(1, length(grid.p))
        for observable in (
            phase_space_integral,
            purity,
            position_density,
            momentum_density,
            phase_space_mean,
            phase_space_covariance,
            wigner_negativity,
        )
            @test_throws DimensionMismatch observable(wrong, grid)
        end
        @test_throws DimensionMismatch expectation((q, p) -> q, wrong, grid)
        @test_throws DimensionMismatch energy(wrong, grid; mass = 1, potential = q -> q^2)
        @test_throws DimensionMismatch overlap(wrong, W, grid)
        @test_throws DimensionMismatch overlap(W, wrong, grid)
    end
end
