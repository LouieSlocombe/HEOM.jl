# Uses drude_fdt_covariance and cold_memory_reference from cold_strong_benchmarks.jl,
# and static_correlation from bath_compression.jl.

# The exact thermal noise spectrum fitted by aaa_bath, at either sign of frequency.
thermal_spectrum(J, ω; kT, hbar = 1.0) =
    2hbar * sign(ω) * J(abs(ω)) / -expm1(-hbar * ω / kT)

# Continuum fluctuation–dissipation equilibrium of a harmonic oscillator, with no bath
# decomposition, as a sum over Matsubara frequencies νₙ = 2πn*kT/ħ. The bath enters
# through kernel(ν) = 2λ - χ(iν) = (2/π)∫J(ω)ν²/(ω(ω²+ν²))dω, the dynamic part of its
# restoring force with the counterterm λq² included (Grabert, Schramm and Ingold,
# Phys. Rep. 168, 115 (1988)). kernel(Inf) = 2λ gives the 1/n² tail beyond `terms`.
function matsubara_fdt_covariance(kernel; mass, omega, kT, hbar, terms = 100_000)
    q, p = 1 / (mass * omega^2), 1.0
    for n in 1:terms
        ν = 2π * n * kT / hbar
        stiffness = mass * omega^2 + kernel(ν)
        q += 2 / (stiffness + mass * ν^2)
        p += 2stiffness / (stiffness + mass * ν^2)
    end
    # Σ_{n>N} 2/(mνₙ²), from Euler–Maclaurin for Σ_{n>N} 1/n².
    tail = 2 / mass * (hbar / (2π * kT))^2 * (1 / terms - 1 / (2terms^2) + 1 / (6terms^3))
    stiffness = mass * omega^2 + kernel(Inf)
    return Diagonal([kT * (q + tail), mass * kT * (p + stiffness * tail)])
end

# Gauss–Legendre nodes and weights for ∫₀^∞ f(x)dx on x = tan(πu/2).
function half_line_rule(nodes)
    rule = eigen(SymTridiagonal(zeros(nodes), [j / sqrt(4j^2 - 1) for j in 1:(nodes-1)]))
    x = tan.(π .* (rule.values .+ 1) ./ 4)
    return x, rule.vectors[1, :] .^ 2 .* (π / 2) .* (1 .+ x .^ 2)
end

# λ + imag(∫₀^∞ C dt)/ħ: the static force constant left by the correlation and
# counterterm. It vanishes for an exactly balanced thermal bath.
static_imbalance(bath) = bath.counterterm + imag(static_correlation(bath)) / bath.hbar

@testset "AAA rational-fit baths" begin
    mass, omega, kT = 1.0, 1.0, 0.1
    frequencies = exp10.(range(-3, 3; length = 600))

    @testset "Continuum reference from Matsubara sums" begin
        (; reorganization, cutoff) = (; reorganization = 0.8, cutoff = 0.5)
        kernel(ν) = isinf(ν) ? 2reorganization : 2reorganization * ν / (cutoff + ν)
        exact = drude_fdt_covariance(; mass, omega, reorganization, cutoff, kT, hbar = 1.0)
        @test matsubara_fdt_covariance(kernel; mass, omega, kT, hbar = 1.0) ≈ exact rtol =
            1e-9
    end

    @testset "Cold Drude bath against a high-order Padé decomposition" begin
        (; reorganization, cutoff) = (; reorganization = 0.8, cutoff = 0.5)
        drude(ω) = 2reorganization * cutoff * ω / (ω^2 + cutoff^2)
        bath = aaa_bath(drude; kT, frequencies, reltol = 1e-8)
        # The fit criterion holds at every sample, at both signs of frequency.
        samples = [-reverse(frequencies); frequencies]
        target = thermal_spectrum.(drude, samples; kT)
        @test maximum(abs, bath_spectrum.(Ref(bath), samples) - target) <=
              1e-8 * maximum(target)
        # Between samples, over the frequencies a cold oscillator resolves.
        between = range(-6, 6; length = 1000)
        @test maximum(
            abs,
            bath_spectrum.(Ref(bath), between) - thermal_spectrum.(drude, between; kT),
        ) <= 1e-8 * maximum(target)
        # All poles of this spectrum are imaginary: only real modes, no oscillations.
        @test length(bath.rates) <= 18
        @test iszero(bath.mixing)
        @test all(isone, bath.weights)
        # Free poles recover the Drude pole and the slowest Matsubara frequencies.
        @test bath.rates[1] ≈ cutoff rtol = 1e-7
        @test bath.rates[2] ≈ 2π * kT rtol = 1e-5
        @test bath.counterterm ≈ reorganization rtol = 1e-12
        @test 0 <= bath.diffusion < 1e-3 * maximum(target)
        @test abs(static_imbalance(bath)) < 1e-10

        pade = drude_lorentz_pade_bath(; reorganization, cutoff, kT, pade = 30)
        times = range(0.5, 40; length = 80)
        @test maximum(
            abs,
            bath_correlation.(Ref(bath), times) - bath_correlation.(Ref(pade), times),
        ) < 5e-8 * abs(bath_correlation(pade, first(times)))
        exact = drude_fdt_covariance(; mass, omega, reorganization, cutoff, kT, hbar = 1.0)
        fitted = harmonic_covariance(bath; mass, omega)
        # About sixteen fitted modes reach the continuum equilibrium more closely than
        # 31 Padé modes.
        @test maximum(abs, fitted - exact) < 5e-8
        @test maximum(abs, fitted - exact) <
              maximum(abs, harmonic_covariance(pade; mass, omega) - exact) / 10

        # Tightening the tolerance adds modes and converges the equilibrium.
        fits = [aaa_bath(drude; kT, frequencies, reltol) for reltol in (1e-4, 1e-6, 1e-8)]
        @test issorted(length.(getfield.(fits, :rates)))
        errors =
            [maximum(abs, harmonic_covariance(fit; mass, omega) - exact) for fit in fits]
        @test issorted(errors; rev = true)
        @test errors[1] < 5e-6

        # The same convention holds for another ħ.
        hbar = 0.7
        scaled = aaa_bath(drude; kT, frequencies, hbar)
        exact = drude_fdt_covariance(; mass, omega, reorganization, cutoff, kT, hbar)
        @test harmonic_covariance(scaled; mass, omega) ≈ exact atol = 1e-7
        @test abs(static_imbalance(scaled)) < 1e-10
    end

    # A Drude background with an underdamped vibration, at kT ≪ ħω.
    background = (; reorganization = 0.3, cutoff = 0.5)
    vibration = (; reorganization = 0.2, frequency = 1.5, damping = 0.2)
    λd, γd = background.reorganization, background.cutoff
    λb, ωb, γb = vibration.reorganization, vibration.frequency, vibration.damping
    structured(ω) =
        2λd * γd * ω / (ω^2 + γd^2) + 2λb * γb * ωb^2 * ω / ((ωb^2 - ω^2)^2 + γb^2 * ω^2)
    # χ(iν) = 2λd*γd/(γd + ν) + 2λb*ωb²/(ωb² + ν² + γb*ν) for these two components.
    structured_kernel(ν) =
        isinf(ν) ? 2(λd + λb) :
        2λd * ν / (γd + ν) + 2λb * ν * (ν + γb) / (ωb^2 + ν^2 + γb * ν)
    structured_bath = aaa_bath(structured; kT, frequencies, reltol = 1e-6)
    structured_exact =
        matsubara_fdt_covariance(structured_kernel; mass, omega, kT, hbar = 1.0)

    @testset "Structured density at low temperature" begin
        bath = structured_bath
        # The vibration is one complex pole pair, so it costs two real modes.
        pair = findall(!iszero, bath.mixing)
        @test length(pair) == 2
        (i, j) = Tuple(pair[1])
        @test abs(i - j) == 1
        @test bath.rates[i] == bath.rates[j]
        @test bath.weights[max(i, j)] == 0
        @test bath.rates[i] ≈ γb / 2 rtol = 1e-6
        @test abs(bath.mixing[i, j]) ≈ sqrt(ωb^2 - γb^2 / 4) rtol = 1e-6
        @test bath.mixing[i, j] == -bath.mixing[j, i]
        @test length(bath.rates) == count(isone, bath.weights) + 1
        @test all(>(0), bath.rates)
        @test bath.counterterm ≈ λd + λb rtol = 1e-12
        @test abs(static_imbalance(bath)) < 1e-7

        samples = [-reverse(frequencies); frequencies]
        target = thermal_spectrum.(structured, samples; kT)
        @test maximum(abs, bath_spectrum.(Ref(bath), samples) - target) <=
              1e-6 * maximum(target)
        fitted = harmonic_covariance(bath; mass, omega)
        @test maximum(abs, fitted - structured_exact) < 1e-6
        @test det(fitted) >= 1 / 4
        # At strong coupling, the reduced equilibrium is not the isolated Gibbs state.
        bare_variance = 1 / (2tanh(1 / (2kT)))
        @test structured_exact[1, 1] < 0.9bare_variance
        @test structured_exact[2, 2] > 1.2bare_variance

        # A fit can be compressed; the harmonic check still holds to the discarded order.
        σ = hankel_singular_values(bath)
        reduced = compress_bath(bath; modes = 8)
        @test length(reduced.rates) == 8
        @test maximum(abs, harmonic_covariance(reduced; mass, omega) - structured_exact) <
              3e-4
        ω = range(-6, 6; length = 601)
        @test maximum(
            abs,
            bath_spectrum.(Ref(reduced), ω) - bath_spectrum.(Ref(bath), ω),
        ) <= 4sqrt(2) * sum(σ[9:end]) + 1e-12 * σ[1]
    end

    @testset "Hierarchy dynamics with a fitted oscillatory bath" begin
        bath = compress_bath(structured_bath; modes = 4)
        @test count(!iszero, bath.mixing) >= 2
        grid = PhaseSpaceGrid((-7, 7), 40, (-7, 7), 40)
        mean0, covariance0 = [0.7, -0.2], [0.5 0.0; 0.0 0.5]
        initial = on_grid(
            (q, p) -> gaussian_wigner(q, p; mean = mean0, covariance = covariance0),
            grid,
        )
        mean, covariance = cold_memory_reference(0.8, mean0, covariance0, bath; mass, omega)
        exact = on_grid((q, p) -> gaussian_wigner(q, p; mean, covariance), grid)
        errors = map((4, 6)) do depth
            prob = heom_problem(
                initial,
                (0.0, 0.8),
                grid;
                mass,
                potential = harmonic_potential(; mass, omega),
                bath,
                depth,
                scaled = true,
            )
            sol = solve(prob, Vern7(); abstol = 1e-10, reltol = 1e-10, saveat = [0.8])
            W = physical_wigner(sol)
            @test phase_space_integral(W, grid) ≈ 1 atol = 1e-9
            @test phase_space_mean(W, grid) ≈ mean atol = 1e-6
            @test phase_space_covariance(W, grid) ≈ covariance atol = 1e-6
            max_error(W, exact)
        end
        @test errors[2] < errors[1] / 10
        @test errors[2] < 1e-6
    end

    @testset "Sub-Ohmic density: low-frequency sampling controls thermalization" begin
        (; exponent, coupling, cutoff) = (; exponent = 0.5, coupling = 0.6, cutoff = 2.0)
        subohmic(ω) = coupling * ω^exponent * cutoff^(1 - exponent) * exp(-ω / cutoff)
        # λ = coupling*cutoff*Γ(s)/π, with Γ(1/2) = √π.
        reorganization = coupling * cutoff / sqrt(π)
        # With ω = cutoff*x², J(ω)/ω dω = 2*coupling*cutoff*exp(-x²)dx is smooth.
        x, w = half_line_rule(400)
        density = @. 2coupling * cutoff * exp(-x^2) * w
        kernel(ν) =
            isinf(ν) ? 2reorganization :
            (2 / π) * sum(@. density * ν^2 / ((cutoff * x^2)^2 + ν^2))
        @test kernel(Inf) ≈ (2 / π) * sum(density) rtol = 1e-12
        exact =
            matsubara_fdt_covariance(kernel; mass, omega, kT, hbar = 1.0, terms = 20_000)
        shifts, imbalances, errors = Float64[], Float64[], Float64[]
        for lowest in (1e-2, 1e-3, 1e-4)
            sampled = exp10.(range(log10(lowest), 2; length = 600))
            bath = aaa_bath(subohmic; kT, frequencies = sampled, reltol = 1e-4)
            # The quadrature integrates the ω^(s-1) singularity of J(ω)/ω.
            @test bath.counterterm ≈ reorganization rtol = 1e-10
            target = thermal_spectrum.(subohmic, sampled; kT)
            @test maximum(abs, bath_spectrum.(Ref(bath), sampled) - target) <=
                  1e-4 * maximum(target)
            covariance = harmonic_covariance(bath; mass, omega)
            push!(shifts, covariance[1, 1] - exact[1, 1])
            push!(imbalances, static_imbalance(bath))
            push!(errors, maximum(abs, covariance - exact))
        end
        # Every fit meets its spectral tolerance, but J(ω)/ω below the sampled range is
        # missing from the fitted correlation. The resulting static imbalance Δ stiffens
        # the classical (zero Matsubara frequency) part of ⟨q²⟩ = kT/(mω² + 2Δ) + …,
        # which dominates the equilibrium error until the spectral tolerance takes over.
        # Both shrink as the sampled range extends.
        @test all(>(0), imbalances)
        @test issorted(imbalances; rev = true)
        @test shifts[1:2] ≈ -2kT .* imbalances[1:2] ./ (mass * omega^2)^2 rtol = 0.1
        @test issorted(errors; rev = true)
        @test errors[1] > 1e-3
        @test errors[3] < errors[1] / 5
    end

    @testset "White-noise remainder and degenerate fits" begin
        # A coarse tolerance is met by white noise alone.
        hot(ω) = 2 * 0.1 * 1.0 * ω / (ω^2 + 1.0)
        white = aaa_bath(hot; kT = 5.0, frequencies, reltol = 0.9, max_modes = 0)
        @test isempty(white.rates)
        @test white.diffusion > 0
        @test white.counterterm ≈ 0.1 rtol = 1e-12
        # No coupling at the samples: no modes, but the counterterm remains.
        silent = aaa_bath(ω -> 0.0; kT = 0.1, frequencies)
        @test isempty(silent.rates)
        @test silent.counterterm == 0
        @test silent.diffusion == 0
        supplied = aaa_bath(ω -> 0.0; kT = 0.1, frequencies, reorganization = 0.3)
        @test supplied.counterterm == 0.3
        # Sample order and duplicates do not matter.
        shuffled = [reverse(frequencies); frequencies[1:7:end]]
        @test aaa_bath(hot; kT = 0.4, frequencies = shuffled).coefficients ==
              aaa_bath(hot; kT = 0.4, frequencies).coefficients
        # A supplied reorganization replaces the quadrature.
        @test aaa_bath(hot; kT = 0.4, frequencies, reorganization = 0.25).counterterm ==
              0.25
    end

    @testset "Argument checks" begin
        drude(ω) = 2 * 0.8 * 0.5 * ω / (ω^2 + 0.25)
        @test_throws ArgumentError aaa_bath(drude; kT = 0, frequencies)
        @test_throws ArgumentError aaa_bath(drude; kT = Inf, frequencies)
        @test_throws ArgumentError aaa_bath(drude; kT = 0.1, frequencies, hbar = 0)
        @test_throws ArgumentError aaa_bath(drude; kT = 0.1, frequencies, reltol = 0)
        @test_throws ArgumentError aaa_bath(drude; kT = 0.1, frequencies, reltol = NaN)
        @test_throws ArgumentError aaa_bath(drude; kT = 0.1, frequencies, max_modes = -1)
        @test_throws ArgumentError aaa_bath(drude; kT = 0.1, frequencies = Float64[])
        @test_throws ArgumentError aaa_bath(drude; kT = 0.1, frequencies = [0.0, 1.0])
        @test_throws ArgumentError aaa_bath(drude; kT = 0.1, frequencies = [1.0, Inf])
        @test_throws ArgumentError aaa_bath(
            drude;
            kT = 0.1,
            frequencies,
            reorganization = -1,
        )
        @test_throws ArgumentError aaa_bath(
            drude;
            kT = 0.1,
            frequencies,
            reorganization = NaN,
        )
        @test_throws "real, finite and nonnegative" aaa_bath(
            ω -> -drude(ω);
            kT = 0.1,
            frequencies,
        )
        @test_throws "real, finite and nonnegative" aaa_bath(
            ω -> drude(ω) * im;
            kT = 0.1,
            frequencies,
        )
        @test_throws "real, finite and nonnegative" aaa_bath(
            ω -> NaN;
            kT = 0.1,
            frequencies,
            reorganization = 0.8,
        )
        # J(ω)/ω ~ 1/ω is not integrable at zero frequency.
        flat(ω) = 1 / (1 + ω^2)
        @test_throws "did not converge" aaa_bath(flat; kT = 0.1, frequencies)
        supplied = aaa_bath(flat; kT = 0.1, frequencies, reltol = 1e-4, reorganization = 1)
        @test supplied.counterterm == 1
        # The thermal spectrum 2kT*J(ω)/ω overflows at a subnormal frequency.
        @test_throws "thermal spectrum must be finite" aaa_bath(
            ω -> 1.0;
            kT = 0.1,
            frequencies = [5e-324, 1.0],
            reorganization = 1,
        )
        @test_throws "did not reach reltol" aaa_bath(
            drude;
            kT = 0.1,
            frequencies,
            reltol = 1e-12,
            max_modes = 4,
        )
        # Too few samples for any rational fit.
        @test_throws "did not reach reltol" aaa_bath(
            drude;
            kT = 0.1,
            frequencies = [1.0],
            reorganization = 0.8,
        )
    end
end
