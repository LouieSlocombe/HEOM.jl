@testset "Harmonic oscillator states" begin
    m, ω, ħ = 1.5, 0.8, 0.7
    H(q, p) = p^2 / (2m) + harmonic_potential(; mass = m, omega = ω)(q)
    grid = PhaseSpaceGrid((-10, 10), 128, (-10, 10), 128)

    for n in 0:3
        W = on_grid((q, p) -> fock_wigner(n, q, p; mass = m, omega = ω, hbar = ħ), grid)
        @test phase_space_integral(W, grid) ≈ 1 atol = 1e-10
        @test purity(W, grid; hbar = ħ) ≈ 1 atol = 1e-10
        @test expectation(H, W, grid) ≈ ħ * ω * (n + 1 / 2) rtol = 1e-10
    end
    @test fock_wigner(1, 0, 0; mass = m, omega = ω, hbar = ħ) ≈ -1 / (π * ħ)
    @test fock_wigner(0, 0.3, -0.2; mass = m, omega = ω, hbar = ħ) ≈
          coherent_wigner(0.3, -0.2; mass = m, omega = ω, hbar = ħ)
    @test_throws ArgumentError fock_wigner(-1, 0, 0; mass = m, omega = ω)

    q0, p0 = 1.2, -0.8
    W = on_grid(
        (q, p) -> coherent_wigner(q, p; q0, p0, mass = m, omega = ω, hbar = ħ),
        grid,
    )
    @test phase_space_integral(W, grid) ≈ 1 atol = 1e-10
    @test purity(W, grid; hbar = ħ) ≈ 1 atol = 1e-10
    @test expectation(H, W, grid) ≈ ħ * ω / 2 + H(q0, p0) rtol = 1e-10

    # The even cat |α⟩ + |-α⟩ has mean occupation α² tanh α², with α² = mωq₀²/2ħ.
    W = on_grid((q, p) -> cat_wigner(q, p; q0 = 2, mass = m, omega = ω, hbar = ħ), grid)
    α² = m * ω * 2^2 / (2ħ)
    @test phase_space_integral(W, grid) ≈ 1 atol = 1e-10
    @test purity(W, grid; hbar = ħ) ≈ 1 atol = 1e-10
    @test expectation(H, W, grid) ≈ ħ * ω * (α² * tanh(α²) + 1 / 2) rtol = 1e-10
    @test minimum(W) < 0

    W0(q, p) = coherent_wigner(q, p; q0, p0, mass = m, omega = ω, hbar = ħ)
    T = 2π / ω
    @test harmonic_evolution(W0, T; mass = m, omega = ω)(0.3, 0.4) ≈ W0(0.3, 0.4)
    @test harmonic_evolution(W0, T / 2; mass = m, omega = ω)(0.3, 0.4) ≈ W0(-0.3, -0.4)
end

@testset "Harmonic oscillator parameter validation and free limit" begin
    for bad in (-1.0, 0.0, Inf, -Inf, NaN)
        for state in (
            (m, ω, ħ) -> coherent_wigner(0, 0; mass = m, omega = ω, hbar = ħ),
            (m, ω, ħ) -> fock_wigner(1, 0, 0; mass = m, omega = ω, hbar = ħ),
            (m, ω, ħ) -> cat_wigner(0, 0; q0 = 1, mass = m, omega = ω, hbar = ħ),
        )
            @test_throws ArgumentError state(bad, 1, 1)
            @test_throws ArgumentError state(1, bad, 1)
            @test_throws ArgumentError state(1, 1, bad)
        end
        @test_throws ArgumentError harmonic_potential(; mass = bad, omega = 1)
        @test_throws ArgumentError harmonic_evolution(
            (q, p) -> q + p,
            1;
            mass = bad,
            omega = 1,
        )
    end
    for bad in (-1.0, Inf, -Inf, NaN)
        @test_throws ArgumentError harmonic_potential(; mass = 1, omega = bad)
        @test_throws ArgumentError harmonic_evolution(
            (q, p) -> q + p,
            1;
            mass = 1,
            omega = bad,
        )
    end
    @test_throws ArgumentError harmonic_evolution((q, p) -> q, Inf; mass = 1, omega = 1)
    @test harmonic_potential(; mass = 1.3, omega = 0)(2.0) == 0
    # The old sin(omega*t)/(m*omega) expression produced NaN at zero frequency.
    W0(q, p) = exp(-q^2 - p^2)
    for ω in (0.0, 1e-100)
        evolved = harmonic_evolution(W0, 0.7; mass = 1.3, omega = ω)
        @test evolved(0.4, 0.8) ≈ W0(0.4 - 0.8 * 0.7 / 1.3, 0.8)
    end
end

@testset "Harmonic oscillator dynamics: coherent and Fock states" begin
    T = 2π
    V = harmonic_potential(; mass = 1, omega = 1)
    H(q, p) = p^2 / 2 + V(q)
    W0(q, p) = coherent_wigner(q, p; q0 = 2, p0 = 1, mass = 1, omega = 1)
    grid = PhaseSpaceGrid((-8, 8), 64, (-8, 8), 64)
    prob = wigner_moyal_problem(on_grid(W0, grid), (0.0, T), grid; mass = 1, potential = V)
    @test prob.u0 == on_grid(W0, grid)
    @test prob.tspan == (0.0, T)

    times = [T / 4, T / 2, T]
    sol = solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12, saveat = times)
    @test sol.t ≈ times
    energy = expectation(H, prob.u0, grid)
    for (t, W) in zip(sol.t, sol.u)
        exact = on_grid(harmonic_evolution(W0, t; mass = 1, omega = 1), grid)
        @test max_error(W, exact) <= 1e-10
        # The centre follows the classical orbit through (2, 1).
        @test expectation((q, p) -> q, W, grid) ≈ 2cos(t) + sin(t) atol = 1e-10
        @test expectation((q, p) -> p, W, grid) ≈ cos(t) - 2sin(t) atol = 1e-10
        @test phase_space_integral(W, grid) ≈ 1 atol = 1e-12
        @test purity(W, grid) ≈ 1 atol = 1e-12
        @test expectation(H, W, grid) ≈ energy rtol = 1e-12
    end

    # The first excited state is stationary and keeps its negative value at the origin.
    fock = on_grid((q, p) -> fock_wigner(1, q, p; mass = 1, omega = 1), grid)
    stationary = remake(prob; u0 = fock, tspan = (0.0, T / 3))
    sol = solve(stationary, Vern9(); abstol = 1e-12, reltol = 1e-12, save_everystep = false)
    @test max_error(sol.u[end], fock) <= 1e-10
    origin = (findfirst(iszero, grid.q), findfirst(iszero, grid.p))
    @test sol.u[end][origin...] ≈ -1 / π atol = 1e-10
end

@testset "Harmonic oscillator dynamics: rotating cat state" begin
    m, ω, ħ = 1.5, 0.8, 0.7
    T = 2π / ω
    W0(q, p) = cat_wigner(q, p; q0 = 2, mass = m, omega = ω, hbar = ħ)
    # 128 points resolve the interference fringes, which oscillate as cos(2q₀p/ħ).
    grid = PhaseSpaceGrid((-8, 8), 128, (-8, 8), 128)
    op = wigner_moyal_operator(
        grid;
        mass = m,
        hbar = ħ,
        potential = harmonic_potential(; mass = m, omega = ω),
    )
    prob = wigner_moyal_problem(on_grid(W0, grid), (0.0, T / 2), op)
    sol = solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12, saveat = [T / 4, T / 2])
    quarter = on_grid(harmonic_evolution(W0, T / 4; mass = m, omega = ω), grid)
    @test minimum(quarter) < 0
    @test max_error(sol.u[1], quarter) <= 1e-10
    # An even cat is symmetric under the half-period rotation (q, p) → (-q, -p).
    @test max_error(sol.u[2], prob.u0) <= 1e-10
    @test purity(sol.u[2], grid; hbar = ħ) ≈ 1 atol = 1e-12
end

@testset "Harmonic oscillator dynamics: finite differences" begin
    T = 2π
    V = harmonic_potential(; mass = 1, omega = 1)
    W0(q, p) = coherent_wigner(q, p; q0 = 2, p0 = 1, mass = 1, omega = 1)
    errors = map((64, 128)) do n
        grid = PhaseSpaceGrid((-8, 8), n, (-8, 8), n)
        op = wigner_moyal_operator(
            grid;
            mass = 1,
            potential = V,
            discretization = FiniteDifference(4),
            moyal_terms = 1,
        )
        prob = wigner_moyal_problem(on_grid(W0, grid), (0.0, T / 4), op)
        sol =
            solve(prob, Vern7(); abstol = 1e-10, reltol = 1e-10, save_everystep = false)
        exact = on_grid(harmonic_evolution(W0, T / 4; mass = 1, omega = 1), grid)
        max_error(sol.u[end], exact)
    end
    @test 3.5 <= log2(errors[1] / errors[2]) <= 4.5
    @test errors[2] <= 2e-3
end
