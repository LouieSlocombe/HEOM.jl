@testset "Public bath correlations and spectra" begin
    @testset "Independent modes, signed arguments and white noise" begin
        coefficients = [0.8 - 0.4im, -0.3 + 0.2im]
        rates = [0.7, 1.9]
        weights = [1.2, -0.4]
        bath = ExponentialBath(
            coefficients,
            rates;
            weights,
            diffusion = 0.13,
            counterterm = 0.9,
        )
        for t in (0.0, 0.1, 1.3, 10.0)
            expected = sum(weights .* coefficients .* exp.(-rates .* t))
            @test bath_correlation(bath, t) ≈ expected
            t > 0 && @test bath_correlation(bath, -t) ≈ conj(expected)
        end
        for ω in (-3.0, -0.2, 0.0, 0.2, 3.0)
            expected = 2real(sum(weights .* coefficients ./ (rates .- im * ω))) + 0.26
            @test bath_spectrum(bath, ω) ≈ expected
        end

        empty = ExponentialBath(ComplexF64[], Float64[])
        noise = ExponentialBath(ComplexF64[], Float64[]; diffusion = 0.3, counterterm = 2)
        for x in (-1.0, 0.0, 1.0)
            @test bath_correlation(empty, x) == 0.0 + 0.0im
            @test bath_spectrum(empty, x) === 0.0
            @test iszero(bath_correlation(noise, x))
            @test bath_spectrum(noise, x) == 0.6
        end
        for bad in (Inf, -Inf, NaN, big"1e1000")
            @test_throws ArgumentError bath_correlation(empty, bad)
            @test_throws ArgumentError bath_spectrum(empty, bad)
        end
        # Evaluating an unphysical supplied correlation must not hide its sign.
        negative = ExponentialBath([-1.0], [1.0])
        @test bath_spectrum(negative, 0) == -2
    end

    @testset "Polynomial and oscillatory correlation bases" begin
        c, γ = [0.5 - 0.2im, -0.7 + 0.1im], 1.2
        jordan = ExponentialBath(c, [γ, γ]; weights = [1, 0], mixing = [0 1; 0 0])
        for t in (0.0, 0.1, 0.7, 2.0)
            expected = (c[1] - t * c[2]) * exp(-γ * t)
            @test bath_correlation(jordan, t) ≈ expected rtol = 3e-15
            t > 0 && @test bath_correlation(jordan, -t) ≈ conj(expected) rtol = 3e-15
        end
        for ω in (-2.0, 0.0, 1.7)
            s = γ - im * ω
            @test bath_spectrum(jordan, ω) ≈ 2real(c[1] / s - c[2] / s^2)
        end

        a, b, damping, Ω = 0.7 + 0.1im, -0.3 - 0.2im, 0.3, 2.0
        # Uncoupled modes are interleaved with the oscillator to exercise the
        # reduced matrix exponential without assuming a contiguous block.
        oscillatory = ExponentialBath(
            [0.4, a, 0.8, b],
            [0.8, damping, 1.7, damping];
            weights = [0.7, 1, -0.2, 0],
            mixing = sparse([2, 4], [4, 2], [-Ω, Ω], 4, 4),
        )
        for t in (0.0, 0.2, 1.1, 3.0)
            expected =
                exp(-damping * t) * (a * cos(Ω * t) + b * sin(Ω * t)) + 0.28exp(-0.8t) -
                0.16exp(-1.7t)
            @test bath_correlation(oscillatory, t) ≈ expected
        end
        for ω in (-3.0, -0.2, 0.0, 0.5, 2.0)
            s = damping - im * ω
            expected = 2real(
                (a * s + b * Ω) / (s^2 + Ω^2) + 0.28 / (0.8 - im * ω) -
                0.16 / (1.7 - im * ω),
            )
            @test bath_spectrum(oscillatory, ω) ≈ expected
        end

        # General nonzero weights remain valid after a real basis change.
        rates, c = [0.7, 1.5], [0.4 - 0.2im, 0.8]
        diagonal = ExponentialBath(c, rates)
        change = [2.0 1.0; -1.0 1.0]
        generator = change * Diagonal(rates) / change
        transformed_rates = diag(generator)
        transformed = ExponentialBath(
            change * c,
            transformed_rates;
            weights = transpose(change) \ ones(2),
            mixing = generator - Diagonal(transformed_rates),
        )
        for t in (-1.1, 0.0, 0.2, 3.0)
            @test bath_correlation(transformed, t) ≈ bath_correlation(diagonal, t)
            @test bath_spectrum(transformed, t) ≈ bath_spectrum(diagonal, t)
        end
    end
end
