using SciMLBase: DiscreteCallback, terminate!

@testset "Dipole response seed and hierarchy correlations" begin
    grid = PhaseSpaceGrid((-6, 6), 48, (-6, 6), 48)
    W = on_grid((q, p) -> exp(-q^2-p^2) / π, grid)
    bath = ExponentialBath([0.2 - 0.1im], [1.0])
    for discretization in (Spectral(), FiniteDifference(4)), scaled in (false, true)
        op = heom_operator(
            grid;
            mass = 1.3,
            potential = q->q^2/2,
            bath,
            depth = 2,
            discretization,
            moyal_terms = 2,
            scaled,
        )
        U = cat(W, 0.2W, -0.3W; dims = 3)
        original = copy(U)
        μ(q) = q + q^3 / 7
        prob = linear_response_problem(U, (2, 3), op; dipole = μ)
        # Independent RHS subtraction cancels the kinetic term. A nonlinear dipole
        # checks that the commutator includes its third momentum derivative.
        dipole_op = wigner_moyal_operator(
            grid;
            mass = 1.3,
            potential = μ,
            discretization,
            moyal_terms = 2,
        )
        free_op = wigner_moyal_operator(
            grid;
            mass = 1.3,
            potential = q->zero(q),
            discretization,
            moyal_terms = 2,
        )
        expected = rhs(free_op, W) - rhs(dipole_op, W)
        @test prob.u0[:, :, 1] ≈ expected atol=1e-14
        @test prob.u0[:, :, 2] ≈ 0.2expected atol=1e-14
        @test prob.u0[:, :, 3] ≈ -0.3expected atol=1e-14
        @test abs(phase_space_integral(physical_wigner(prob.u0), grid)) < 1e-14
        @test U == original
        @test prob.u0 !== U
        factorized = linear_response_problem(W, (0, 1), op)
        @test iszero(factorized.u0[:, :, 2:end])
        if discretization isa Spectral
            @test factorized.u0[:, :, 1] ≈ on_grid((q, p)->2p*exp(-q^2-p^2)/π, grid) atol=1e-10
        else
            sparse_prob = linear_response_problem(U, (0, 1), op; jacobian = :sparse)
            @test sparse_prob.u0 ≈ linear_response_problem(U, (0, 1), op).u0
        end
    end
end

@testset "Response helpers and validation" begin
    grid = PhaseSpaceGrid((-5, 5), 16, (-5, 5), 16)
    W = on_grid((q, p) -> coherent_wigner(q, p; mass = 1, omega = 1), grid)
    bath = ExponentialBath(ComplexF64[], Float64[])
    op = heom_operator(grid; mass = 1, potential = q->q^2/2, bath, depth = 0)
    eq = equilibrate(W, (0, 0), op, Vern7(); stationarity_abstol = 0.1)
    @test eq.converged
    @test linear_response_problem(eq, (1, 2)).u0 ≈ linear_response_problem(W, (1, 2), op).u0
    response =
        linear_response(eq, (3, 3.1), Vern7(); saveat = 0.05, abstol = 1e-9, reltol = 1e-9)
    @test response isa LinearResponseResult
    @test first(response.times) == 0
    @test last(response.times) ≈ 0.1
    @test response.solution.prob.p === op
    @test absorption_spectrum(response; frequencies = [1.0]) ==
          absorption_spectrum(response.times, response.response; frequencies = [1.0])
    nonstationary = equilibrate(W .* reshape(1 .+ grid.q .^ 2, :, 1), (0, 0), op, Vern7())
    @test !nonstationary.converged
    @test_throws ArgumentError linear_response_problem(nonstationary, (0, 1))
    @test_throws ArgumentError linear_response(nonstationary, (0, 1), Vern7())
    for disc in (Spectral(), FiniteDifference(4))
        cl = caldeira_leggett_operator(
            grid;
            mass = 1,
            potential = q->q^2/2,
            friction = 0.2,
            kT = 1,
            discretization = disc,
            moyal_terms = 1,
        )
        sol = linear_response(W, (0, 0.1), cl, Vern7(); observable = q->2q, saveat = 0.05)
        @test all(isfinite, sol.response)
        @test sol.response[end] > 0
        @test_throws ArgumentError linear_response_problem(
            W,
            (0, 1),
            cl;
            jacobian = :sparse,
        )
    end
    wm = wigner_moyal_operator(grid; mass = 1, potential = q->q^2/2)
    @test_throws ArgumentError linear_response_problem(W, (0, 1), wm; jacobian = :sparse)
    @test_throws ArgumentError linear_response_problem(W, (0, 1), op; dipole = (q, t)->q*t)
    for interval in ((0, 0), (1, 0), (0, Inf), (0, 1im), (0,), (-1e308, 1e308))
        @test_throws ArgumentError linear_response_problem(W, interval, op)
    end
    driven = heom_operator(grid; mass = 1, potential = (q, t)->q^2/2-t*q, bath, depth = 0)
    @test_throws ArgumentError linear_response_problem(W, (0, 1), driven)
    for options in (
        (; save_start = false),
        (; save_end = false),
        (; save_on = false),
        (; save_idxs = 1),
    )
        @test_throws ArgumentError linear_response(W, (0, 1), op, Vern7(); options...)
    end
    @test_throws ArgumentError linear_response(W, (0, 1), op, Vern7(); observable = q->Inf)
    # Avoid silently accepting a partial, solver-terminated response.
    stop = DiscreteCallback((u, t, i)->t>0.01, terminate!)
    @test_throws ErrorException linear_response(W, (0, 1), op, Vern7(); callback = stop)
end

@testset "Causal absorption transform" begin
    # Analytic causal susceptibility of a damped unit-mass harmonic oscillator.
    times = collect(0:0.01:60)
    η, ω0 = 0.4, 1.3
    signal = sin.(ω0 .* times) ./ ω0
    frequencies = [0.0, 0.5, 1.3, 2.0]
    spectrum = absorption_spectrum(times, signal; frequencies, broadening = η)
    exact = 1 ./ (ω0^2 .+ (η .- im .* frequencies) .^ 2)
    @test maximum(abs, spectrum.susceptibility - exact) < 1e-5
    @test spectrum.intensity ≈ frequencies .* imag.(exact) atol=1e-8
    @test spectrum.intensity[1] == 0
    @test spectrum.intensity[3] > 0
    shifted = absorption_spectrum(times .+ 5, signal; frequencies, broadening = η)
    @test shifted.susceptibility ≈ spectrum.susceptibility atol=1e-14
    irregular = [0.0, 0.03, 0.2, 1.0]
    # A constant zero-frequency response integrates exactly even on irregular samples.
    @test absorption_spectrum(irregular, fill(2.0, 4); frequencies = [0.0]).susceptibility ==
          [2.0+0im]
    @test isempty(
        absorption_spectrum(irregular, ones(4); frequencies = Float64[]).intensity,
    )
    @test_throws DimensionMismatch absorption_spectrum([0, 1], [0]; frequencies)
    @test_throws ArgumentError absorption_spectrum([0], [0]; frequencies)
    for bad_times in ([0, 0], [1, 0], [0, Inf], [0, 1im], [-1e308, 1e308])
        @test_throws ArgumentError absorption_spectrum(bad_times, [0, 1]; frequencies)
    end
    for bad_signal in ([0, NaN], [0, 1im])
        @test_throws ArgumentError absorption_spectrum([0, 1], bad_signal; frequencies)
    end
    for bad_frequencies in ([-1.0], [Inf], [1im])
        @test_throws ArgumentError absorption_spectrum(
            [0, 1],
            [0, 1];
            frequencies = bad_frequencies,
        )
    end
    for broadening in (-1, Inf, NaN)
        @test_throws ArgumentError absorption_spectrum(
            [0, 1],
            [0, 1];
            frequencies,
            broadening,
        )
    end
    @test_throws ArgumentError absorption_spectrum(
        [0, 1],
        [1e308, 1e308];
        frequencies = [0.0],
    )
    @test_throws ArgumentError absorption_spectrum([0, 1e308], [0, 1]; frequencies = [10.0])
end
