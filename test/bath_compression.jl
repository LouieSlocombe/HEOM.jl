# Uses cold_memory_matrices, cold_memory_reference and drude_fdt_covariance from
# cold_strong_benchmarks.jl as independent harmonic references.

# The same correlation in another real basis, x' = Q*x for orthogonal Q.
function rotate_bath(bath, Q)
    Γ = Q * Matrix(Diagonal(bath.rates) + bath.mixing) * Q'
    rates = diag(Γ)
    return ExponentialBath(
        Q * bath.coefficients,
        rates;
        hbar = bath.hbar,
        diffusion = bath.diffusion,
        counterterm = bath.counterterm,
        weights = Q * bath.weights,
        mixing = Γ - Diagonal(rates),
    )
end

# A deterministic orthogonal matrix close to the identity.
rotation(n, angle) = exp(angle * [sign(j - i) / (i + j) for i in 1:n, j in 1:n])

# ∫₀^∞ C(t) dt of the retained exponentials, from the definition of the basis.
static_correlation(bath) =
    dot(bath.weights, (Diagonal(bath.rates) + Matrix(bath.mixing)) \ bath.coefficients)

# Frequencies spanning every rate of the decomposition, including ω = 0.
function spectral_samples(bath)
    largest = maximum(abs, eigvals(Diagonal(bath.rates) + Matrix(bath.mixing)))
    return sinh.(range(-1, 1; length = 1601) .* asinh(20largest))
end

@testset "Balanced bath compression" begin
    cold = (; reorganization = 0.8, cutoff = 0.5, kT = 0.1, hbar = 1.0)
    pade = drude_lorentz_pade_bath(; cold..., pade = 7)
    matsubara =
        drude_lorentz_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.3, matsubara = 40)
    brownian = brownian_oscillator_bath(;
        reorganization = 0.3,
        frequency = 1.5,
        damping = 0.4,
        kT = 0.1,
        matsubara = 6,
    )
    combined = combine_baths(pade, brownian)

    @testset "Hankel singular values of exact and redundant realizations" begin
        c, ν = 0.6 - 0.8im, 1.5
        # One exponential has P = |c|²/2ν and Q = 1/2ν, so σ = |c|/2ν.
        single = ExponentialBath([c], [ν]; diffusion = 0.2, counterterm = 0.3)
        @test hankel_singular_values(single) ≈ [abs(c) / (2ν)] rtol = 1e-14
        # Splitting one exponential into two equal halves adds an invisible state.
        split = ExponentialBath([c / 2, c / 2], [ν, ν]; diffusion = 0.2, counterterm = 0.3)
        σ = hankel_singular_values(split)
        @test σ[1] ≈ abs(c) / (2ν) rtol = 1e-14
        @test σ[2] < sqrt(eps()) * σ[1]
        for method in (:truncate, :residualize)
            reduced = compress_bath(split; modes = 1, method)
            @test length(reduced.rates) == 1
            @test reduced.rates ≈ [ν] rtol = 1e-14
            @test reduced.diffusion ≈ 0.2 atol = 1e-15
            @test reduced.counterterm ≈ 0.3 atol = 1e-15
            for t in (0.0, 0.4, 3.0)
                @test bath_correlation(reduced, t) ≈ bath_correlation(single, t) rtol =
                    1e-13
            end
        end
        @test hankel_singular_values(ExponentialBath(ComplexF64[], Float64[])) == Float64[]
        @test hankel_singular_values(ExponentialBath(zeros(2), [1.0, 2.0])) == zeros(2)
    end

    @testset "Invariance under a change of real basis" begin
        for bath in (pade, brownian)
            rotated = rotate_bath(bath, rotation(length(bath.rates), 0.4))
            @test !iszero(rotated.mixing)
            for t in (0.0, 0.3, 2.0)
                @test bath_correlation(rotated, t) ≈ bath_correlation(bath, t) rtol = 1e-12
            end
            σ, σ_rotated = hankel_singular_values(bath), hankel_singular_values(rotated)
            @test σ_rotated ≈ σ atol = 1e-12 * σ[1]
            for method in (:truncate, :residualize), modes in (2, 4)
                bath === brownian && method === :residualize && continue
                reduced = compress_bath(bath; modes, method)
                from_rotated = compress_bath(rotated; modes, method)
                @test from_rotated.diffusion ≈ reduced.diffusion rtol = 1e-10
                @test from_rotated.counterterm ≈ reduced.counterterm rtol = 1e-10
                for ω in (-2.0, 0.0, 0.7, 5.0)
                    @test bath_spectrum(from_rotated, ω) ≈ bath_spectrum(reduced, ω) rtol =
                        1e-9
                end
            end
        end
    end

    @testset "Frequency-uniform spectral bound" begin
        for (bath, methods) in (
            (pade, (:truncate, :residualize)),
            (matsubara, (:truncate, :residualize)),
            (brownian, (:truncate,)),
            (combined, (:truncate,)),
            # A Jordan block at an exact Drude–Matsubara coincidence.
            (
                drude_lorentz_bath(;
                    reorganization = 0.8,
                    cutoff = 2π * 0.1,
                    kT = 0.1,
                    matsubara = 3,
                ),
                (:truncate, :residualize),
            ),
        )
            σ = hankel_singular_values(bath)
            @test issorted(σ; rev = true)
            ω = spectral_samples(bath)
            source = bath_spectrum.(Ref(bath), ω)
            for modes in 1:count(>(1e-6σ[1]), σ), method in methods
                reduced = compress_bath(bath; modes, method)
                @test length(reduced.rates) == modes
                change = maximum(abs, bath_spectrum.(Ref(reduced), ω) - source)
                @test change <= 4sqrt(2) * sum(σ[(modes+1):end]) + 1e-12 * σ[1]
            end
        end
        # More retained modes approach the source monotonically for a Padé bath.
        ω = spectral_samples(pade)
        source = bath_spectrum.(Ref(pade), ω)
        changes = map(1:6) do modes
            maximum(abs, bath_spectrum.(Ref(compress_bath(pade; modes)), ω) - source)
        end
        @test issorted(changes; rev = true)
        @test changes[4] < 1e-3 * maximum(abs, source)
    end

    @testset "Truncation keeps diffusion and counterterm; residualization keeps statics" begin
        for bath in (pade, matsubara), modes in 1:4
            truncated = compress_bath(bath; modes)
            @test truncated.diffusion === bath.diffusion
            @test truncated.counterterm === bath.counterterm
            residualized = compress_bath(bath; modes, method = :residualize)
            @test residualized.diffusion > bath.diffusion
            # ∫C dt and its white-noise remainder: S(0) and the static force constant.
            @test bath_spectrum(residualized, 0) ≈ bath_spectrum(bath, 0) rtol = 1e-12
            static, reduced_static =
                static_correlation(bath), static_correlation(residualized)
            # The Drude counterterm cancels this force constant exactly.
            @test residualized.counterterm + imag(reduced_static) / bath.hbar ≈
                  bath.counterterm + imag(static) / bath.hbar atol = 1e-12
            @test residualized.diffusion + real(reduced_static) ≈
                  bath.diffusion + real(static) rtol = 1e-12
        end
        # Residualizing every mode of a terminated Drude bath leaves its exact Markovian
        # limit: noise 2λkT/γ and no net static potential.
        (; reorganization, cutoff, kT) = (; reorganization = 0.8, cutoff = 0.5, kT = 0.3)
        markovian = compress_bath(matsubara; modes = 0, method = :residualize)
        @test isempty(markovian.rates)
        @test markovian.diffusion ≈ 2reorganization * kT / cutoff rtol = 1e-12
        @test markovian.counterterm ≈ 0 atol = 1e-12
        dropped = compress_bath(matsubara; modes = 0)
        @test isempty(dropped.rates)
        @test dropped.diffusion === matsubara.diffusion
        @test dropped.counterterm === matsubara.counterterm
        # Fast negative Matsubara terms cannot become a nonnegative white noise.
        @test_throws "negative diffusion" compress_bath(
            brownian;
            modes = 3,
            method = :residualize,
        )
        # Nor can a fast dissipative mode without a counterterm become a static potential.
        uncompensated = ExponentialBath([0.5, 0.01 - 0.3im], [1.0, 50.0])
        @test_throws "negative counterterm" compress_bath(
            uncompensated;
            modes = 1,
            method = :residualize,
        )
        @test compress_bath(uncompensated; modes = 1).counterterm == 0
    end

    @testset "Real Schur form of the reduced rate matrix" begin
        for (bath, modes) in ((pade, 4), (brownian, 2), (brownian, 4), (combined, 5))
            reduced = compress_bath(bath; modes)
            Γ = Diagonal(reduced.rates) + Matrix(reduced.mixing)
            @test all(>(0), reduced.rates)
            @test all(iszero, diag(reduced.mixing))
            @test all(isfinite, reduced.weights)
            # Below the diagonal, only the coupling of a 2×2 complex-pair block.
            for j in 1:modes, i in (j+1):modes
                iszero(Γ[i, j]) && continue
                @test i == j + 1
                @test reduced.rates[i] ≈ reduced.rates[j] rtol = 1e-12
                @test Γ[i, j] * Γ[j, i] < 0
            end
        end
        # The underdamped oscillator is a complex pair, so it needs two real modes.
        @test !istriu(compress_bath(brownian; modes = 2).mixing)
        eigenvalues = eigvals(Matrix(Diagonal(pade.rates) + pade.mixing))
        @test all(isreal, eigenvalues)
    end

    @testset "Hierarchy roots do not depend on the reduced real basis" begin
        reduced = compress_bath(pade; modes = 3)
        rotated = rotate_bath(reduced, rotation(3, 0.15))
        @test all(>(0), rotated.rates)
        grid = PhaseSpaceGrid((-6, 6), 24, (-6, 6), 24)
        initial = on_grid(
            (q, p) ->
                coherent_wigner(q, p; q0 = 0.6, p0 = -0.3, mass = 1.0, omega = 1.0),
            grid,
        )
        potential(q) = q^2 / 2 + 0.08q^4
        for scaled in (false, true)
            roots = map((reduced, rotated)) do bath
                prob = heom_problem(
                    initial,
                    (0.0, 0.4),
                    grid;
                    mass = 1.0,
                    potential,
                    bath,
                    depth = 3,
                    scaled,
                )
                sol = solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12)
                copy(physical_wigner(sol))
            end
            @test max_error(roots...) < 1e-9
        end
    end

    @testset "Compressed cold hierarchy against exact Gaussian dynamics" begin
        mass, omega = 1.0, 1.0
        reduced = compress_bath(pade; modes = 3)
        grid = PhaseSpaceGrid((-8, 8), 48, (-8, 8), 48)
        mean0, covariance0 = [0.7, -0.2], [0.5 0.0; 0.0 0.5]
        initial = on_grid(
            (q, p) -> gaussian_wigner(q, p; mean = mean0, covariance = covariance0),
            grid,
        )
        prob = heom_problem(
            initial,
            (0.0, 0.8),
            grid;
            mass,
            potential = harmonic_potential(; mass, omega),
            bath = reduced,
            depth = 6,
            scaled = true,
        )
        sol = solve(prob, Vern7(); abstol = 1e-10, reltol = 1e-10, saveat = [0.8])
        W = physical_wigner(sol)
        own = cold_memory_reference(0.8, mean0, covariance0, reduced; mass, omega)
        source = cold_memory_reference(0.8, mean0, covariance0, pade; mass, omega)
        # The hierarchy is exact for the compressed bath it is given, to the depth and
        # grid errors of an uncompressed Padé bath with the same parameters.
        exact = on_grid(
            (q, p) -> gaussian_wigner(q, p; mean = own[1], covariance = own[2]),
            grid,
        )
        @test max_error(W, exact) < 1e-6
        @test phase_space_mean(W, grid) ≈ own[1] atol = 2e-6
        @test phase_space_covariance(W, grid) ≈ own[2] atol = 5e-6
        # ...and three balanced modes reproduce the eight-mode Padé source closely.
        @test own[1] ≈ source[1] atol = 1e-4
        @test own[2] ≈ source[2] atol = 1e-4
        three_poles = drude_lorentz_pade_bath(; cold..., pade = 2)
        few_poles = cold_memory_reference(0.8, mean0, covariance0, three_poles; mass, omega)
        @test maximum(abs, own[2] - source[2]) < maximum(abs, few_poles[2] - source[2]) / 10
    end

    @testset "Harmonic equilibrium covariance" begin
        mass, omega = 1.3, 0.9
        coincident = drude_lorentz_bath(;
            reorganization = 0.8,
            cutoff = 2π * 0.1,
            kT = 0.1,
            matsubara = 2,
        )
        crude = compress_bath(pade; modes = 1)
        for bath in (
            pade,
            matsubara,
            brownian,
            combined,
            coincident,
            crude,
            compress_bath(pade; modes = 3, method = :residualize),
            compress_bath(brownian; modes = 4),
        )
            covariance = harmonic_covariance(bath; mass, omega)
            _, _, joint = cold_memory_matrices(bath; mass, omega)
            @test covariance isa Symmetric
            @test covariance ≈ joint[1:2, 1:2] rtol = 1e-10
            bath === crude || @test det(covariance) >= bath.hbar^2 / 4
        end
        # A one-mode compression is a stable, bounded-error bath whose equilibrium
        # violates the uncertainty relation; the harmonic check exposes it.
        @test det(harmonic_covariance(crude; mass, omega)) < pade.hbar^2 / 4
        # Increasing Padé order converges to the continuum fluctuation–dissipation result.
        exact = drude_fdt_covariance(; mass, omega, cold...)
        errors = map((2, 4, 8)) do order
            bath = drude_lorentz_pade_bath(; cold..., pade = order)
            maximum(abs, harmonic_covariance(bath; mass, omega) - exact)
        end
        @test issorted(errors; rev = true)
        @test errors[3] < 1e-3 * maximum(exact)
        # A compressed bath inherits its source equilibrium, not merely its spectrum.
        for modes in (4, 5)
            reduced = compress_bath(pade; modes)
            @test harmonic_covariance(reduced; mass, omega) ≈
                  harmonic_covariance(pade; mass, omega) rtol = 1e-3
        end
        # The classical limit depends only on the static potential: kT/(mω²) and m*kT.
        hot = drude_lorentz_bath(;
            reorganization = 0.5,
            cutoff = 2.0,
            kT = 50.0,
            matsubara = 2,
        )
        classical = harmonic_covariance(hot; mass, omega)
        @test classical[1, 1] ≈ 50 / (mass * omega^2) rtol = 1e-3
        @test classical[2, 2] ≈ 50mass rtol = 1e-2
        @test abs(classical[1, 2]) < 1e-10
        # No friction: no unique stationary state.
        @test_throws ArgumentError harmonic_covariance(
            ExponentialBath(ComplexF64[], Float64[]);
            mass,
            omega,
        )
        @test_throws ArgumentError harmonic_covariance(
            ExponentialBath(ComplexF64[], Float64[]; diffusion = 0.2);
            mass,
            omega,
        )
        @test_throws ArgumentError harmonic_covariance(
            ExponentialBath([0.3 - 0.1im], [1.0]; weights = [0.0]);
            mass,
            omega,
        )
        @test_throws ArgumentError harmonic_covariance(pade; mass = 0, omega)
        @test_throws ArgumentError harmonic_covariance(pade; mass, omega = Inf)
    end

    @testset "Argument checks and copies" begin
        @test_throws ArgumentError compress_bath(pade; modes = 2, method = :balanced)
        @test_throws ArgumentError compress_bath(pade; modes = -1)
        @test_throws ArgumentError compress_bath(pade; modes = 9)
        # No coupling: every Hankel singular value is zero.
        uncoupled = ExponentialBath(zeros(2), [1.0, 2.0]; diffusion = 0.1)
        @test_throws ArgumentError compress_bath(uncoupled; modes = 1)
        @test compress_bath(uncoupled; modes = 0).diffusion == 0.1
        # Unresolved trailing values cannot define a retained subspace.
        σ = hankel_singular_values(brownian)
        @test σ[7] - σ[8] < sqrt(eps()) * σ[1]
        @test_throws ArgumentError compress_bath(brownian; modes = 7)
        empty = ExponentialBath(ComplexF64[], Float64[]; hbar = 0.7, counterterm = 0.2)
        @test compress_bath(empty; modes = 0).counterterm == 0.2
        @test compress_bath(empty; modes = 0).hbar == 0.7
        for method in (:truncate, :residualize)
            same = compress_bath(brownian; modes = 8, method)
            for field in fieldnames(ExponentialBath)
                @test getfield(same, field) == getfield(brownian, field)
            end
            same.coefficients[1] = 9
            same.mixing[1, 2] = 9
            @test brownian.coefficients[1] != 9
            @test brownian.mixing[1, 2] != 9
        end
    end
end
