# Right-hand-side checks against closed forms. The harmonic oscillator has V‴ = 0, so these
# are the tests of the quantum terms.

const m_test, ħ_test = 1.3, 0.7
const a_test, b_test, q0_test, p0_test = 1.0, 0.8, 0.4, -0.3

"""
Gaussian test function `exp(-a(q - q0)² - b(p - p0)²)` whose derivatives are known.
"""
gaussian(q, p) = exp(-a_test * (q - q0_test)^2 - b_test * (p - p0_test)^2)

"""
Kinetic term `-(p/m) ∂W/∂q` of the Wigner–Moyal equation for the Gaussian test function.
"""
gaussian_kinetic(q, p) = 2a_test * (p / m_test) * (q - q0_test) * gaussian(q, p)

@testset "Quartic double well: the Moyal series ends at ħ²" begin
    λ, μ, ε = 0.1, 0.5, 0.1
    V(q) = λ * q^4 - μ * q^2 + ε * q
    dV(q) = 4λ * q^3 - 2μ * q + ε
    # -(p/m)∂W/∂q + V′∂W/∂p - (ħ²/24)V‴∂³W/∂p³ with V‴ = 24λq and
    # ∂³exp(-bu²)/∂u³ = (12b²u - 8b³u³)exp(-bu²).
    function quartic_rhs(q, p; quantum = true)
        u, b = p - p0_test, b_test
        classical = gaussian_kinetic(q, p) - 2b * u * dV(q) * gaussian(q, p)
        correction = -ħ_test^2 * λ * q * (12b^2 * u - 8b^3 * u^3) * gaussian(q, p)
        return classical + quantum * correction
    end
    operator(grid; kwargs...) =
        wigner_moyal_operator(grid; mass = m_test, hbar = ħ_test, potential = V, kwargs...)

    grid = PhaseSpaceGrid((-7, 7), 64, (-7, 7), 64)
    W = on_grid(gaussian, grid)
    exact = on_grid(quartic_rhs, grid)
    classical = on_grid((q, p) -> quartic_rhs(q, p; quantum = false), grid)
    scale = maximum(abs, exact)

    for moyal_terms in (nothing, 2, 3, 4)
        @test max_error(rhs(operator(grid; moyal_terms), W), exact) <= 1e-12 * scale
    end
    @test max_error(rhs(operator(grid; moyal_terms = 1), W), classical) <= 1e-12 * scale
    @test max_error(classical, exact) > 0.1 * scale

    # Finite differences converge to the same right-hand side at their accuracy order.
    for order in (2, 4, 6)
        errors = map((64, 128)) do n
            fine = PhaseSpaceGrid((-7, 7), n, (-7, 7), n)
            op = operator(fine; discretization = FiniteDifference(order), moyal_terms = 2)
            max_error(rhs(op, on_grid(gaussian, fine)), on_grid(quartic_rhs, fine))
        end
        @test log2(errors[1] / errors[2]) ≈ order atol = 0.5
    end
end

@testset "Sine potential: every order of the Moyal series" begin
    A, k = 0.8, 1.3
    V(q) = A * sin(k * q)
    # The exact potential term shifts W in momentum by ±ħk/2.
    function sine_rhs(q, p)
        shifts = gaussian(q, p + ħ_test * k / 2) - gaussian(q, p - ħ_test * k / 2)
        return gaussian_kinetic(q, p) + A / ħ_test * cos(k * q) * shifts
    end
    operator(moyal_terms) = wigner_moyal_operator(
        grid;
        mass = m_test,
        hbar = ħ_test,
        potential = V,
        moyal_terms,
    )

    grid = PhaseSpaceGrid((-7, 7), 64, (-7, 7), 64)
    W = on_grid(gaussian, grid)
    exact = on_grid(sine_rhs, grid)
    scale = maximum(abs, exact)

    @test max_error(rhs(operator(nothing), W), exact) <= 1e-12 * scale
    # Each extra term of the truncated series brings it much closer to the exact operator.
    errors = [max_error(rhs(operator(moyal_terms), W), exact) for moyal_terms in 1:4]
    @test all(errors[2:end] .< errors[1:(end-1)] ./ 10)
end

@testset "Harmonic oscillator: the Moyal series ends at the classical term" begin
    grid = PhaseSpaceGrid((-8, 8), 64, (-8, 8), 64)
    V = harmonic_potential(; mass = 1.0, omega = 1.0)
    exact = wigner_moyal_operator(grid; mass = 1.0, potential = V)
    liouville = wigner_moyal_operator(grid; mass = 1.0, potential = V, moyal_terms = 1)
    coherent =
        on_grid((q, p) -> coherent_wigner(q, p; q0 = 2, p0 = 1, mass = 1, omega = 1), grid)
    @test max_error(rhs(exact, coherent), rhs(liouville, coherent)) <= 1e-12
    fock = on_grid((q, p) -> fock_wigner(1, q, p; mass = 1, omega = 1), grid)
    @test maximum(abs, rhs(exact, fock)) <= 1e-12
end

@testset "Conservation structure" begin
    V(q) = 0.1q^4 - 0.5q^2
    grid = PhaseSpaceGrid((-7, 7), 48, (-7, 7), 48)
    W = on_grid((q, p) -> gaussian(q, p) * (1 + 0.3 * sin(q * p)), grid)
    for (discretization, moyal_terms) in
        ((Spectral(), nothing), (Spectral(), 2), (FiniteDifference(4), 2))
        op = wigner_moyal_operator(
            grid;
            mass = m_test,
            hbar = ħ_test,
            potential = V,
            discretization,
            moyal_terms,
        )
        dW = rhs(op, W)
        # The operators are skew-symmetric, so the norm ΣW and the purity ΣW² are conserved.
        @test abs(sum(dW)) <= 1e-12 * sum(abs, dW)
        @test abs(sum(W .* dW)) <= 1e-12 * sum(abs, W .* dW)
    end

    fd = wigner_moyal_operator(
        grid;
        mass = m_test,
        hbar = ħ_test,
        potential = V,
        discretization = FiniteDifference(4),
        moyal_terms = 2,
    )
    L = sparse(fd)
    @test L isa SparseMatrixCSC{Float64,Int}
    @test size(L) == (48^2, 48^2)
    @test iszero(norm(L + transpose(L)))
    @test L * vec(W) ≈ vec(rhs(fd, W))
end

@testset "Spectral right-hand side is type stable and allocation free" begin
    grid = PhaseSpaceGrid((-6, 6), 32, (-6, 6), 48)
    op = wigner_moyal_operator(grid; mass = 1.0, potential = q -> q^4 / 4 - q^2 / 2)
    W = on_grid(gaussian, grid)
    dW = similar(W)
    @test (@inferred wigner_moyal!(dW, W, op, 0.0)) === nothing
    @test rhs_allocations(dW, W, op) == 0
end

@testset "Operator construction and validation" begin
    grid = PhaseSpaceGrid((-4, 4), 16, (-4, 4), 16)
    V = harmonic_potential(; mass = 1, omega = 1)
    op = wigner_moyal_operator(grid; mass = 2, hbar = 0.5, potential = V)
    @test repr(op) ==
          "SpectralWignerMoyal($(repr(grid)), mass = 2.0, hbar = 0.5, moyal_terms = exact)"
    fd = wigner_moyal_operator(
        grid;
        mass = 2,
        potential = V,
        discretization = FiniteDifference(6),
        moyal_terms = 3,
    )
    @test repr(fd) ==
          "FiniteDifferenceWignerMoyal($(repr(grid)), order = 6, mass = 2.0, " *
          "hbar = 1.0, moyal_terms = 3)"

    @test_throws ArgumentError wigner_moyal_operator(grid; mass = 0, potential = V)
    @test_throws ArgumentError wigner_moyal_operator(
        grid;
        mass = 1,
        hbar = -1,
        potential = V,
    )
    for moyal_terms in (0, 5)
        @test_throws ArgumentError wigner_moyal_operator(
            grid;
            mass = 1,
            potential = V,
            moyal_terms,
        )
    end
    # The exact operator is nonlocal in momentum, so finite differences need a truncation.
    @test_throws ArgumentError wigner_moyal_operator(
        grid;
        mass = 1,
        potential = V,
        discretization = FiniteDifference(),
    )
    # The potential overflows inside the box, and so do the shifted points of the exact
    # operator.
    overflowing(q) = exp(1000q)
    @test_throws ArgumentError wigner_moyal_operator(
        grid;
        mass = 1,
        potential = overflowing,
    )
    @test_throws ArgumentError wigner_moyal_operator(
        grid;
        mass = 1,
        potential = overflowing,
        moyal_terms = 1,
    )
    # A fourth-order first-derivative stencil needs five points.
    @test_throws ArgumentError wigner_moyal_operator(
        PhaseSpaceGrid((-1, 1), 4, (-1, 1), 4);
        mass = 1,
        potential = V,
        discretization = FiniteDifference(4),
        moyal_terms = 1,
    )
    @test_throws DimensionMismatch wigner_moyal_problem(zeros(8, 8), (0.0, 1.0), op)
end

@testset "Reject invalid Hamiltonians and malformed Wigner states" begin
    grid = PhaseSpaceGrid((-4, 4), 16, (-5, 5), 20)
    V(q) = q^2 / 2
    for invalid in (Inf, NaN, big"1e1000", big"1e-1000")
        @test_throws ArgumentError wigner_moyal_operator(
            grid;
            mass = invalid,
            potential = V,
        )
        # Even the classical truncation must reject an unrepresentable Planck constant.
        @test_throws ArgumentError wigner_moyal_operator(
            grid;
            mass = 1,
            potential = V,
            hbar = invalid,
            moyal_terms = 1,
        )
    end
    for (discretization, moyal_terms) in
        ((Spectral(), nothing), (Spectral(), 2), (FiniteDifference(4), 2))
        @test_throws ArgumentError wigner_moyal_operator(
            grid;
            mass = 1e-320,
            potential = V,
            discretization,
            moyal_terms,
        )
        @test_throws ArgumentError wigner_moyal_operator(
            PhaseSpaceGrid((-1e-310, 1e-310), 16, (-5, 5), 20);
            mass = 1,
            potential = V,
            discretization,
            moyal_terms,
        )
        for potential in (q -> im * q^2, q -> q^2 + im, q -> Inf, q -> NaN)
            @test_throws ArgumentError wigner_moyal_operator(
                grid;
                mass = 1,
                potential,
                discretization,
                moyal_terms,
            )
        end
        op = wigner_moyal_operator(
            grid;
            mass = 1,
            potential = V,
            discretization,
            moyal_terms,
        )
        W = zeros(size(grid))
        @test_throws ArgumentError wigner_moyal_problem(complex.(W), (0, 1), op)
        # Broadcasting a single column and reshaping the same number of entries both used
        # to bypass the phase-space axis contract in the direct RHS interfaces.
        for malformed in (zeros(16, 1), zeros(20, 16))
            @test_throws DimensionMismatch wigner_moyal!(similar(W), malformed, op, 0.0)
            @test_throws DimensionMismatch wigner_moyal!(malformed, W, op, 0.0)
        end
        for invalid in (Inf, NaN, big"1e1000")
            @test_throws ArgumentError wigner_moyal_problem(
                fill(invalid, size(grid)),
                (0, 1),
                op,
            )
        end
    end
end

@testset "Moyal skew symmetry on rectangular odd and even FFT grids" begin
    for (nq, np) in ((15, 18), (16, 19))
        grid = PhaseSpaceGrid((-3, 3), nq, (-4, 4), np)
        # Broad Fourier content exercises the last paired mode and the unpaired Nyquist
        # mode, which are almost absent from smooth Gaussian tests.
        A = [sin(0.7i * j) for i in 1:nq, j in 1:np]
        B = [cos(0.4i * j) for i in 1:nq, j in 1:np]
        op = wigner_moyal_operator(
            grid;
            mass = 1.3,
            hbar = 0.6,
            potential = q -> sin(q) + q^4,
        )
        LA, LB = rhs(op, A), rhs(op, B)
        @test abs(dot(A, LB) + dot(LA, B)) <= 1e-13 * norm(A) * norm(LB)
        @test abs(sum(LA)) <= 1e-13 * sum(abs, LA)
    end
end
