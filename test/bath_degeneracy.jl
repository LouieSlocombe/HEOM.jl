@testset "Coincident thermal poles" begin
    bath_matrix(b) = Matrix(Diagonal(b.rates) + b.mixing)
    correlation(b, t) = dot(b.weights, exp(-bath_matrix(b) * t) * b.coefficients)
    spectrum(b, ω) =
        2real(dot(b.weights, (bath_matrix(b) - im * ω * I) \ b.coefficients)) + 2b.diffusion

    @testset "Matsubara divided differences and exact double poles" begin
        λ, θ, ħ = 0.7, 0.03, 0.8
        # Independent arbitrary-precision evaluation of the original (singular
        # in Float64) exponential formula, without regularizing its cotangent.
        function reference(γ, K, t)
            setprecision(256) do
                coupling, rate, temperature, planck = BigFloat.((λ, γ, θ, ħ))
                c0 = coupling * planck * rate * (cot(planck * rate / (2temperature)) - im)
                result = c0 * exp(-rate * t)
                for k in 1:K
                    ν = 2big(π) * k * temperature / planck
                    c = 4coupling * rate * temperature * ν / (ν^2 - rate^2)
                    result += c * exp(-ν * t)
                end
                return ComplexF64(result)
            end
        end
        for pole in (1, 2, 5), displacement in (0.0, -1e-12, 1e-12, -1e-6, 1e-6)
            γ = 2π * pole * θ / ħ * (1 + displacement)
            bath = drude_lorentz_bath(;
                reorganization = λ,
                cutoff = γ,
                kT = θ,
                hbar = ħ,
                matsubara = pole + 2,
            )
            @test all(isfinite, bath.coefficients)
            @test bath.weights[pole+1] == 0
            @test bath.mixing[1, pole+1] == 1
            @test maximum(abs, bath.coefficients) < 10
            @test spectrum(bath, 0.0) ≈ 4λ * θ / γ rtol = 2e-12
            for t in (0.0, 0.25, 1.0, 5.0)
                @test correlation(bath, t) ≈ reference(γ, pole + 2, t) atol = 3e-14
            end
        end
        γ = 2π * θ / ħ
        exact = drude_lorentz_bath(;
            reorganization = λ,
            cutoff = γ,
            kT = θ,
            hbar = ħ,
            matsubara = 1,
            terminator = false,
        )
        @test exact.coefficients[1] ≈ λ * θ - im * λ * ħ * γ
        @test exact.coefficients[2] ≈ 2λ * γ * θ
        for t in (0.0, 0.2, 1.0, 4.0)
            @test correlation(exact, t) ≈
                  (λ * θ - im * λ * ħ * γ - 2λ * γ * θ * t) * exp(-γ * t)
        end
        @test_throws ArgumentError drude_lorentz_bath(;
            reorganization = λ,
            cutoff = γ,
            kT = θ,
            hbar = ħ,
            matsubara = 0,
            terminator = false,
        )
    end

    @testset "Padé double poles preserve the Bose spectrum" begin
        λ, θ, ħ = 0.6, 0.03, 0.7
        for N in (1, 4, 8), displacement in (0.0, -1e-12, 1e-12, -1e-5, 1e-5)
            ξ, η, R = HEOM.bose_pade(N)
            γ = first(ξ) * (θ / ħ) * (1 + displacement)
            bath = drude_lorentz_pade_bath(;
                reorganization = λ,
                cutoff = γ,
                kT = θ,
                hbar = ħ,
                pade = N,
            )
            @test all(isfinite, bath.coefficients)
            @test bath.weights[2] == 0
            @test bath.mixing[1, 2] == 1
            @test maximum(abs, bath.coefficients) < 10
            @test spectrum(bath, 0.0) ≈ 4λ * θ / γ rtol = 2e-12
            for ω in (0.02, 0.1, 0.5, 2.0)
                x = ħ * ω / θ
                coth_approx = 2 / x + 4sum(η .* x ./ (x^2 .+ ξ .^ 2)) + 2R * x
                spectral_density = 2λ * γ * ω / (ω^2 + γ^2)
                @test spectrum(bath, ω) ≈ ħ * spectral_density * (coth_approx + 1)
                @test spectrum(bath, -ω) ≈ ħ * spectral_density * (coth_approx - 1)
            end
        end
        # A physical Matsubara coincidence is now valid even at high Padé order,
        # where the first Padé pole equals the Matsubara pole to machine precision.
        γ = 2π * θ / ħ
        for N in (8, 16, 32)
            bath = drude_lorentz_pade_bath(;
                reorganization = λ,
                cutoff = γ,
                kT = θ,
                hbar = ħ,
                pade = N,
            )
            for ω in (0.03, 0.1, 0.3)
                target = 4ħ * λ * γ * ω / (ω^2 + γ^2) / (-expm1(-ħ * ω / θ))
                @test spectrum(bath, ω) ≈ target rtol = 1e-12
            end
        end
    end
end
