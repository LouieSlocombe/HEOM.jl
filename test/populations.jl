@testset "Populations and rates" begin
    @testset "Trigonometric interval weights" begin
        for n in (16, 15)
            lo, hi = -2.0, 3.0
            h = (hi - lo) / n
            x = [lo + (i - 1) * h for i in 1:n]
            a, b = -0.73, 1.26
            weights = HEOM.interval_weights(x, h, (a, b))
            @test sum(weights) ≈ b - a atol = 1e-14
            for k in 1:((n-1)÷2)
                κ = 2π * k / (hi - lo)
                f = cos.(κ .* (x .- 0.4))
                integral = (sin(κ * (b - 0.4)) - sin(κ * (a - 0.4))) / κ
                @test dot(weights, f) ≈ integral atol = 1e-14
            end
            if iseven(n)
                κ = π / h
                f = cos.(κ .* (x .- lo))
                integral = (sin(κ * (b - lo)) - sin(κ * (a - lo))) / κ
                @test dot(weights, f) ≈ integral atol = 1e-14
            end
            @test HEOM.interval_weights(x, h, (lo, hi)) == fill(h, n)
            @test HEOM.interval_weights(x, h, (-Inf, Inf)) == fill(h, n)
            @test HEOM.interval_weights(x, h, (-Inf, b)) ==
                  HEOM.interval_weights(x, h, (lo, b))
            @test HEOM.interval_weights(x, h, (a, Inf)) ==
                  HEOM.interval_weights(x, h, (a, hi))
            @test HEOM.interval_weights(x, h, (hi + 1, Inf)) == zeros(n)
            @test HEOM.interval_weights(x, h, (-Inf, lo - 1)) == zeros(n)
            @test HEOM.interval_weights(x, h, (a, a)) == zeros(n)
            @test_throws ArgumentError HEOM.interval_weights(x, h, (1, 0))
            @test_throws ArgumentError HEOM.interval_weights(x, h, (NaN, 0))
        end
    end

    @testset "Window probabilities" begin
        grid = PhaseSpaceGrid((-8, 8), 64, (-8, 8), 80)
        r = 0.4
        covariance = [0.64 0.4r; 0.4r 0.25]
        W = on_grid((q, p) -> gaussian_wigner(q, p; mean = [0.0, 0.0], covariance), grid)
        @test probability(W, grid) ≈ phase_space_integral(W, grid) atol = 1e-14
        @test probability(W, grid; q = (0, Inf)) ≈ 0.5 atol = 1e-12
        @test probability(W, grid; p = (0, Inf)) ≈ 0.5 atol = 1e-12
        @test probability(W, grid; q = (0, Inf), p = (0, Inf)) ≈ 1 / 4 + asin(r) / 2π atol =
            1e-12
        for dividing_surface in (-0.73, 0.0, 1.26)
            @test probability(W, grid; q = (-Inf, dividing_surface)) +
                  probability(W, grid; q = (dividing_surface, Inf)) ≈
                  phase_space_integral(W, grid) atol = 1e-14
            @test probability(W, grid; p = (-Inf, dividing_surface)) +
                  probability(W, grid; p = (dividing_surface, Inf)) ≈
                  phase_space_integral(W, grid) atol = 1e-14
        end
        # The signed test function exposes the error in grid-point indicator windows.
        signed = on_grid((q, p) -> q * exp(-q^2 - p^2) / sqrt(π), grid)
        a, b = -0.71, 1.13
        @test probability(signed, grid; q = (a, b)) ≈ (exp(-a^2) - exp(-b^2)) / 2 atol =
            1e-12
        @test probability(2W, grid; q = (0, Inf)) ≈ 1.0 atol = 1e-12
        @test_throws ArgumentError probability(W, grid; q = (1, 0))
        @test_throws ArgumentError probability(W, grid; p = (NaN, 0))
        @test_throws DimensionMismatch probability(ones(1, length(grid.p)), grid)
    end

    @testset "Coherent-state current and flux" begin
        mass, omega, hbar = 1.3, 0.8, 0.7
        q0, p0 = 0.4, -0.3
        grid = PhaseSpaceGrid((-8, 8), 128, (-8, 8), 128)
        W = on_grid((q, p) -> coherent_wigner(q, p; q0, p0, mass, omega, hbar), grid)
        current = probability_current(W, grid; mass)
        @test current ≈ (p0 / mass) .* position_density(W, grid) atol = 1e-12
        op = wigner_moyal_operator(
            grid;
            mass,
            hbar,
            potential = harmonic_potential(; mass, omega),
        )
        c = 0.17
        q_flux =
            (p0 / mass) *
            sqrt(mass * omega / (π * hbar)) *
            exp(-mass * omega * (c - q0)^2 / hbar)
        p_flux =
            -mass * omega^2 * q0 * exp(-(c - p0)^2 / (mass * omega * hbar)) /
            sqrt(π * mass * omega * hbar)
        @test probability_rate(W, op; q = (c, Inf)) ≈ q_flux atol = 1e-12
        @test probability_rate(W, op; p = (c, Inf)) ≈ p_flux atol = 1e-12
        @test probability_rate(W, op) ≈ 0 atol = 1e-12
        @test expectation_rate((q, p) -> 1, W, op) ≈ 0 atol = 1e-12
        @test_throws DimensionMismatch probability_current(
            ones(1, length(grid.p)),
            grid;
            mass,
        )
        @test_throws DimensionMismatch probability_rate(ones(1, length(grid.p)), op)
        @test_throws DimensionMismatch expectation_rate(
            (q, p) -> q,
            ones(1, length(grid.p)),
            op,
        )
    end

    @testset "Ehrenfest and semi-discrete rates" begin
        λ, μ, ε = 0.1, 0.5, 0.1
        mass, hbar = 1.3, 0.7
        V(q) = λ * q^4 - μ * q^2 + ε * q
        dV(q) = 4λ * q^3 - 2μ * q + ε
        H(q, p) = p^2 / (2mass) + V(q)
        grid = PhaseSpaceGrid((-7, 7), 64, (-7, 7), 64)
        W = on_grid((q, p) -> exp(-(q - 0.4)^2 - 0.8 * (p + 0.3)^2), grid)
        current = probability_current(W, grid; mass)
        for (discretization, moyal_terms) in
            ((Spectral(), nothing), (Spectral(), 2), (FiniteDifference(4), 2))
            op = wigner_moyal_operator(
                grid;
                mass,
                hbar,
                potential = V,
                discretization,
                moyal_terms,
            )
            @test expectation_rate((q, p) -> 1, W, op) ≈ 0 atol = 1e-12
            @test expectation_rate((q, p) -> q, W, op) ≈
                  expectation((q, p) -> p / mass, W, grid) atol = 1e-12
            @test expectation_rate((q, p) -> p, W, op) ≈
                  -expectation((q, p) -> dV(q), W, grid) atol = 1e-12
            @test expectation_rate(H, W, op) ≈ 0 atol = 1e-12
            c = grid.q[37]
            flux = probability_rate(W, op; q = (c, Inf))
            if discretization isa Spectral
                @test flux ≈ current[37] atol = 1e-13
            else
                # The finite-difference continuity equation uses D_q from the operator.
                Dq = HEOM.periodic_difference_matrix(1, 4, length(grid.q), grid.dq)
                @test position_density(rhs(op, W), grid) ≈ -Dq * current atol = 1e-13
                weights = HEOM.interval_weights(grid.q, grid.dq, (c, Inf))
                @test flux ≈ -dot(weights, Dq * current) atol = 1e-13
            end
            # A central finite difference of each linear observable agrees with its RHS
            # rate for both discretisations, including a joint quasi-probability window.
            dW = rhs(op, W)
            dt = 0.01
            window = (; q = (-0.71, 1.13), p = (-0.43, 0.89))
            finite_rate =
                (
                    probability(W + dt * dW, grid; window...) -
                    probability(W - dt * dW, grid; window...)
                ) / (2dt)
            @test probability_rate(W, op; window...) ≈ finite_rate atol = 1e-12
            finite_energy_rate =
                (expectation(H, W + dt * dW, grid) - expectation(H, W - dt * dW, grid)) /
                (2dt)
            @test expectation_rate(H, W, op) ≈ finite_energy_rate atol = 1e-12
            @test_throws DimensionMismatch probability_rate(ones(1, length(grid.p)), op)
            @test_throws DimensionMismatch expectation_rate(
                (q, p) -> q,
                ones(1, length(grid.p)),
                op,
            )
        end
    end
end

@testset "Driven instantaneous rates" begin
    grid = PhaseSpaceGrid((-8, 8), 64, (-8, 8), 64)
    mass, omega, q0, p0 = 1.3, 0.8, 0.4, -0.3
    W = on_grid((q, p) -> coherent_wigner(q, p; mass, omega, q0, p0), grid)
    V0 = harmonic_potential(; mass, omega)
    E(t) = 0.6sin(t)
    potential = DrivenPotential(V0, E, identity)
    bath = ExponentialBath([0.0], [1.0])
    t = 0.7
    c = 0.17
    force = -mass * omega^2 * q0 + E(t)
    momentum_flux = force * exp(-(c - p0)^2 / (mass * omega)) / sqrt(π * mass * omega)
    for (discretization, moyal_terms) in
        ((Spectral(), nothing), (Spectral(), 1), (FiniteDifference(4), 1))
        h = wigner_moyal_operator(grid; mass, potential, discretization, moyal_terms)
        cl = caldeira_leggett_operator(
            grid;
            mass,
            potential,
            discretization,
            moyal_terms,
            friction = 0,
            kT = 1,
        )
        heom = heom_operator(
            grid;
            mass,
            potential,
            discretization,
            moyal_terms,
            bath,
            depth = 1,
        )
        U = cat(W, zeros(size(grid)); dims = 3)
        for (op, state) in ((h, W), (cl, W), (heom, U))
            @test expectation_rate((q, p) -> p, state, op; t) ≈ force atol = 1e-11
            @test expectation_rate((q, p) -> p, state, op) ≈ -mass * omega^2 * q0 atol =
                1e-11
            @test expectation_rate((q, p) -> q, state, op; t) ≈ p0 / mass atol = 1e-11
            @test probability_rate(state, op; t) ≈ 0 atol = 1e-12
            if discretization isa Spectral
                @test probability_rate(state, op; t, p = (c, Inf)) ≈ momentum_flux atol =
                    1e-11
            end
        end
    end
end
