function caldeira_rhs(op, W)
    dW = similar(W)
    caldeira_leggett!(dW, W, op, 0.0)
    return dW
end

function caldeira_rhs_allocations(dW, W, op)
    caldeira_leggett!(dW, W, op, 0.0)
    return @allocated caldeira_leggett!(dW, W, op, 0.0)
end

@testset "Caldeira–Leggett right-hand side" begin
    m, ħ, γ, kT = 1.3, 0.7, 0.6, 0.9
    a, b, q0, p0 = 1.0, 0.8, 0.4, -0.3
    λ, μ, ε = 0.1, 0.5, 0.1
    V(q) = λ * q^4 - μ * q^2 + ε * q
    W0(q, p) = exp(-a * (q - q0)^2 - b * (p - p0)^2)
    function exact_rhs(q, p)
        u = p - p0
        kinetic = 2a * p / m * (q - q0)
        potential = -2b * u * (4λ * q^3 - 2μ * q + ε)
        quantum = -ħ^2 * λ * q * (12b^2 * u - 8b^3 * u^3)
        # γ ∂p(pW) + mγkT ∂p²W: γ is the full momentum damping rate.
        friction = γ * (1 - 2b * p * u)
        diffusion = m * γ * kT * (4b^2 * u^2 - 2b)
        return (kinetic + potential + quantum + friction + diffusion) * W0(q, p)
    end
    grid = PhaseSpaceGrid((-7, 7), 64, (-7, 7), 64)
    W, exact = on_grid(W0, grid), on_grid(exact_rhs, grid)
    for moyal_terms in (nothing, 2)
        op = caldeira_leggett_operator(
            grid;
            mass = m,
            potential = V,
            friction = γ,
            kT,
            hbar = ħ,
            moyal_terms,
        )
        @test max_error(caldeira_rhs(op, W), exact) <= 1e-11
    end

    for order in (2, 4, 6)
        errors = map((64, 128)) do n
            fine = PhaseSpaceGrid((-7, 7), n, (-7, 7), n)
            op = caldeira_leggett_operator(
                fine;
                mass = m,
                potential = V,
                friction = γ,
                kT,
                hbar = ħ,
                discretization = FiniteDifference(order),
                moyal_terms = 2,
            )
            max_error(caldeira_rhs(op, on_grid(W0, fine)), on_grid(exact_rhs, fine))
        end
        @test log2(errors[1] / errors[2]) ≈ order atol = 0.5
    end

    for (discretization, moyal_terms) in
        ((Spectral(), nothing), (Spectral(), 2), (FiniteDifference(4), 2))
        op = caldeira_leggett_operator(
            grid;
            mass = m,
            potential = V,
            friction = γ,
            kT,
            hbar = ħ,
            discretization,
            moyal_terms,
        )
        dW = caldeira_rhs(op, W)
        @test abs(sum(dW)) <= 1e-12 * sum(abs, dW)
        @test (@inferred caldeira_leggett!(dW, W, op, 0.0)) === nothing
        if discretization isa Spectral
            @test caldeira_rhs_allocations(dW, W, op) == 0
        end
        @test probability_rate(W, op) ≈ 0 atol = 1e-12
        @test expectation_rate((q, p) -> p^2, W, op) ≈ expectation((q, p) -> p^2, dW, grid) atol =
            1e-12

        isolated = caldeira_leggett_operator(
            grid;
            mass = m,
            potential = V,
            friction = 0,
            kT,
            hbar = ħ,
            discretization,
            moyal_terms,
        )
        wm = wigner_moyal_operator(
            grid;
            mass = m,
            potential = V,
            hbar = ħ,
            discretization,
            moyal_terms,
        )
        @test caldeira_rhs(isolated, W) ≈ rhs(wm, W) atol = 1e-12
        if discretization isa FiniteDifference
            L = sparse(op)
            @test L isa SparseMatrixCSC{Float64,Int}
            @test size(L) == (length(W), length(W))
            @test L * vec(W) ≈ vec(dW) atol = 1e-12
            @test maximum(abs, vec(sum(L; dims = 1))) < 1e-11
            L[1, 1] += 1
            @test caldeira_rhs(op, W) ≈ dW atol = 1e-12
        end
    end
end

@testset "Caldeira–Leggett construction and validation" begin
    grid = PhaseSpaceGrid((-5, 5), 32, (-5, 5), 32)
    V = harmonic_potential(; mass = 1, omega = 1)
    for value in (-1.0, 0.0, Inf, -Inf, NaN)
        @test_throws ArgumentError caldeira_leggett_operator(
            grid;
            mass = value,
            potential = V,
            friction = 1,
            kT = 1,
        )
        @test_throws ArgumentError caldeira_leggett_operator(
            grid;
            mass = 1,
            hbar = value,
            potential = V,
            friction = 1,
            kT = 1,
        )
    end
    for value in (-1.0, Inf, -Inf, NaN)
        @test_throws ArgumentError caldeira_leggett_operator(
            grid;
            mass = 1,
            potential = V,
            friction = value,
            kT = 1,
        )
        @test_throws ArgumentError caldeira_leggett_operator(
            grid;
            mass = 1,
            potential = V,
            friction = 1,
            kT = value,
        )
    end
    @test_throws ArgumentError caldeira_leggett_operator(
        grid;
        mass = 1e200,
        potential = V,
        friction = 1e200,
        kT = 1,
    )
    op = caldeira_leggett_operator(grid; mass = 1, potential = V, friction = 0.5, kT = 0)
    @test repr(op) ==
          "SpectralCaldeiraLeggett($(repr(grid)), mass = 1.0, " *
          "hbar = 1.0, friction = 0.5, kT = 0.0, moyal_terms = exact)"
    fd = caldeira_leggett_operator(
        grid;
        mass = 1,
        potential = V,
        friction = 0.5,
        kT = 2,
        discretization = FiniteDifference(4),
        moyal_terms = 1,
    )
    @test repr(fd) ==
          "FiniteDifferenceCaldeiraLeggett($(repr(grid)), order = 4, " *
          "mass = 1.0, hbar = 1.0, friction = 0.5, kT = 2.0, moyal_terms = 1)"
    W = on_grid((q, p) -> coherent_wigner(q, p; mass = 1, omega = 1), grid)
    prob = caldeira_leggett_problem(W, (0.0, 1.0), op)
    @test prob.u0 == W
    @test prob.p === op
    @test prob.tspan == (0.0, 1.0)
    @test_throws DimensionMismatch caldeira_leggett_problem(zeros(4, 4), (0.0, 1.0), op)
    @test_throws ArgumentError caldeira_leggett_operator(
        grid;
        mass = 1,
        potential = V,
        friction = 0.5,
        kT = 1,
        discretization = FiniteDifference(),
    )
    # With no diffusion, the dissipator alone is γ ∂p(pW).
    isolated = wigner_moyal_operator(grid; mass = 1, potential = V)
    expected = on_grid(
        (q, p) -> 0.5 * (1 - 2p^2) * coherent_wigner(q, p; mass = 1, omega = 1),
        grid,
    )
    @test max_error(caldeira_rhs(op, W) - rhs(isolated, W), expected) < 1e-8

    # The even-grid Nyquist mode has zero odd derivative, but its second derivative
    # must still damp. Subtracting the friction-only RHS isolates momentum diffusion.
    alternating = [(-1.0)^j for i in eachindex(grid.q), j in eachindex(grid.p)]
    warm = caldeira_leggett_operator(grid; mass = 1, potential = V, friction = 0.5, kT = 2)
    diffusion_only = caldeira_rhs(warm, alternating) - caldeira_rhs(op, alternating)
    @test diffusion_only ≈ -(π / grid.dp)^2 .* alternating atol = 1e-11
end

@testset "Damped harmonic oscillator: analytical Gaussian evolution" begin
    m, ω, ħ, kT = 1.3, 0.9, 0.4, 1.1
    V = harmonic_potential(; mass = m, omega = ω)
    H(q, p) = p^2 / (2m) + V(q)
    grid = PhaseSpaceGrid((-12, 12), 64, (-12, 12), 64)
    μ0 = [2.0, 0.6]
    Σ∞ = [kT / (m * ω^2) 0; 0 m * kT]
    Σ0 = 1.25Σ∞ + [0 0.1; 0.1 0]
    W0 = on_grid((q, p) -> gaussian_wigner(q, p; mean = μ0, covariance = Σ0), grid)
    equilibrium =
        on_grid((q, p) -> gaussian_wigner(q, p; mean = zeros(2), covariance = Σ∞), grid)

    # For the harmonic potential the full Wigner equation is an Ornstein–Uhlenbeck
    # equation. Its Gaussian solution has μ(t) = exp(At)μ₀ and
    # Σ(t) = Σ∞ + exp(At)(Σ₀ - Σ∞)exp(Aᵀt), valid in every damping regime.
    for γ in (1.2ω, 2ω, 3ω)
        underdamped = γ < 2ω
        times = underdamped ? [0.0, 0.5 / ω, 2 / ω, 6 / ω, 20 / ω] : [0.0, 0.5 / ω, 2 / ω]
        A = [0 1 / m; -m * ω^2 -γ]
        prob = caldeira_leggett_problem(
            W0,
            (0.0, last(times)),
            grid;
            mass = m,
            potential = V,
            friction = γ,
            kT,
            hbar = ħ,
        )
        @test maximum(abs, caldeira_rhs(prob.p, equilibrium)) < 1e-10
        @test expectation_rate((q, p) -> q, W0, prob.p) ≈ μ0[2] / m atol = 1e-8
        @test expectation_rate((q, p) -> p, W0, prob.p) ≈ -m * ω^2 * μ0[1] - γ * μ0[2] atol =
            1e-8
        @test expectation_rate(H, W0, prob.p) ≈ γ * (kT - (Σ0[2, 2] + μ0[2]^2) / m) atol =
            1e-8
        sol = solve(prob, Vern7(); abstol = 1e-10, reltol = 1e-10, saveat = times)
        @test sol.t ≈ times
        for (t, W) in zip(sol.t, sol.u)
            propagator = exp(A * t)
            μt = propagator * μ0
            Σt = Σ∞ + propagator * (Σ0 - Σ∞) * propagator'
            exact =
                on_grid((q, p) -> gaussian_wigner(q, p; mean = μt, covariance = Σt), grid)
            @test max_error(W, exact) < 1e-8
            @test phase_space_integral(W, grid) ≈ 1 atol = 1e-10
            @test phase_space_mean(W, grid) ≈ μt atol = 1e-7
            @test phase_space_covariance(W, grid) ≈ Σt atol = 1e-7
            exact_energy = (Σt[2, 2] + μt[2]^2) / (2m) + m * ω^2 * (Σt[1, 1] + μt[1]^2) / 2
            @test energy(W, grid; mass = m, potential = V) ≈ exact_energy atol = 1e-7
        end

        if underdamped
            trajectory = diagnostics(sol; potential = V)
            @test trajectory.t == times
            excess_energy = trajectory.energy .- kT
            @test all(diff(excess_energy) .< 0)
            @test abs(last(excess_energy)) < 1e-7
            @test norm(phase_space_mean(sol.u[end], grid)) < 5e-5
            # The center reaches the well bottom; a finite-temperature bath leaves
            # a thermal width and energy kT, rather than a delta function at q=p=0.
            @test phase_space_covariance(sol.u[end], grid) ≈ Σ∞ atol = 1e-7
            @test max_error(sol.u[end], equilibrium) < 5e-6
            @test last(trajectory.var_q) > 0.5
            @test last(trajectory.var_p) > 0.5
        end
    end
end
