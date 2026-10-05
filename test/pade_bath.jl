@testset "Low-temperature Padé baths" begin
    @testset "Bose poles from an independent continued fraction" begin
        # Evaluate coth(x/2) by backwards continued-fraction recursion. This
        # neither computes an eigensystem nor uses the pole/residue formula.
        function continued_coth(x, N)
            denominator = 4N + 3.0
            for j in (2N):-1:1
                denominator = 2j + 1 + (x / 2)^2 / denominator
            end
            return 2 / x + (x / 2) / denominator
        end
        for N in (0, 1, 2, 4, 8, 32)
            ξ, η, R = HEOM.bose_pade(N)
            @test length(ξ) == length(η) == N
            @test issorted(ξ)
            @test all(>(0), ξ)
            @test all(>(0), η)
            for x in (0.01, 0.7, 4.0, 20.0, 100.0, 1000.0)
                approx = 2 / x + 4sum(η .* x ./ (x^2 .+ ξ .^ 2)) + 2R * x
                @test approx ≈ continued_coth(x, N) rtol = 3e-13
            end
        end
        ξ, η, R = HEOM.bose_pade(1)
        @test only(ξ) ≈ sqrt(42)
        @test only(η) ≈ 49 / 40
        @test R == 1 / 40
        ξ, η, R = HEOM.bose_pade(16)
        for x in range(0.1, 100.0; length = 40)
            approx = 2 / x + 4sum(η .* x ./ (x^2 .+ ξ .^ 2)) + 2R * x
            @test approx ≈ coth(x / 2) rtol = 2e-9
        end
    end

    @testset "Fluctuation–dissipation spectrum and quantum asymmetry" begin
        λ, γ, ħ = 0.8, 0.7, 0.9
        spectrum(b, ω) =
            2real(
                dot(
                    b.weights,
                    (Diagonal(b.rates) + b.mixing - im * ω * I) \ b.coefficients,
                ),
            ) + 2b.diffusion
        # The exact target uses only the spectral density and Bose factor.
        # It is independent of either exponential decomposition.
        J(ω) = 2λ * γ * ω / (ω^2 + γ^2)
        exact(ω, θ) = 2ħ * J(ω) / (-expm1(-ħ * ω / θ))
        for θ in (0.1, 0.01)
            errors = Float64[]
            for N in (4, 8, 16, 32)
                bath = drude_lorentz_pade_bath(;
                    reorganization = λ,
                    cutoff = γ,
                    kT = θ,
                    hbar = ħ,
                    pade = N,
                )
                push!(
                    errors,
                    maximum(
                        abs(spectrum(bath, ω) / exact(ω, θ) - 1) for
                        ω in range(0.01, 3.0; length = 101)
                    ),
                )
                @test bath.counterterm == λ
                @test bath.hbar == ħ
                @test bath.diffusion ≈ λ * γ * ħ^2 / (2θ * (N + 1) * (2N + 3))
                @test spectrum(bath, 0.0) ≈ 4λ * θ / γ rtol = 2e-11
                @test imag(bath.coefficients[1]) == -λ * ħ * γ
                for ω in (0.2, 1.0, 3.0)
                    @test spectrum(bath, ω) - spectrum(bath, -ω) ≈ 2ħ * J(ω)
                    @test spectrum(bath, -ω) >= -1e-12
                end
            end
            @test errors[2] < errors[1] / 5
            @test errors[3] < max(errors[2] / 100, 2e-12)
            @test errors[4] < 2e-11
        end

        # At βħγ = 63 there are negative thermal residues. Clipping them would
        # destroy both the thermal spectrum and the zero-frequency sum rule.
        bath = drude_lorentz_pade_bath(;
            reorganization = λ,
            cutoff = γ,
            kT = 0.01,
            hbar = ħ,
            pade = 16,
        )
        @test any(<(0), real.(bath.coefficients[2:end]))
        # An exact Matsubara collision need not be a collision in a Padé basis.
        bath =
            drude_lorentz_pade_bath(; reorganization = λ, cutoff = 2π, kT = 1.0, pade = 1)
        @test all(isfinite, bath.coefficients)
        @test all(isfinite, bath.rates)
    end

    @testset "Options and validation" begin
        args = (reorganization = 1.2, cutoff = 0.7, kT = 0.04, hbar = 1.3)
        b = drude_lorentz_pade_bath(; args...)
        omitted = drude_lorentz_pade_bath(; args..., terminator = false)
        @test length(b.rates) == 5
        @test omitted.diffusion == 0
        @test omitted.coefficients == b.coefficients
        @test omitted.rates == b.rates
        zero = drude_lorentz_pade_bath(; args..., reorganization = 0)
        @test isempty(zero.rates)
        @test isempty(zero.coefficients)
        @test zero.diffusion == zero.counterterm == 0
        @test zero.hbar == args.hbar
        classical = drude_lorentz_pade_bath(; args..., pade = 0)
        @test length(classical.rates) == 1
        @test real(classical.coefficients[1]) ≈
              2args.reorganization * args.kT - args.cutoff * classical.diffusion
        for bad in (-1.0, Inf, NaN)
            @test_throws ArgumentError drude_lorentz_pade_bath(;
                args...,
                reorganization = bad,
            )
        end
        for key in (:cutoff, :kT, :hbar), bad in (-1.0, 0.0, Inf, NaN)
            invalid = NamedTuple{(key,)}((bad,))
            @test_throws ArgumentError drude_lorentz_pade_bath(; args..., invalid...)
        end
        for bad in (-1, typemax(Int))
            @test_throws ArgumentError drude_lorentz_pade_bath(; args..., pade = bad)
        end
        @test_throws ArgumentError drude_lorentz_pade_bath(;
            reorganization = floatmax(Float64),
            cutoff = 2,
            kT = 1,
            pade = 1,
        )
    end
end
