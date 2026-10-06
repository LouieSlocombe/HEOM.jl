@testset "Underdamped Brownian oscillator baths" begin
    λ, ω₀, γ, θ, ħ = 0.35, 1.7, 0.6, 0.4, 0.8
    parameters = (; reorganization = λ, frequency = ω₀, damping = γ, kT = θ, hbar = ħ)
    J(ω) = 2λ * γ * ω₀^2 * ω / ((ω₀^2 - ω^2)^2 + γ^2 * ω^2)
    thermal_spectrum(ω) = iszero(ω) ? 4λ * γ * θ / ω₀^2 : 2ħ * J(ω) / (-expm1(-ħ * ω / θ))
    function reference_correlation(parameters, K, t)
        setprecision(256) do
            λb, ωb, γb, θb, ħb = BigFloat.((
                parameters.reorganization,
                parameters.frequency,
                parameters.damping,
                parameters.kT,
                parameters.hbar,
            ))
            αb = γb / 2
            Ωb = sqrt(ωb^2 - αb^2)
            thermal = coth(ħb * complex(Ωb, αb) / (2θb))
            value =
                λb * ħb * ωb^2 / Ωb *
                exp(-αb * t) *
                (real(thermal) * cos(Ωb * t) - (imag(thermal) + im) * sin(Ωb * t))
            for k in 1:K
                νb = 2big(π) * k * θb / ħb
                cb = -4λb * γb * ωb^2 * θb * νb / ((ωb^2 + νb^2)^2 - γb^2 * νb^2)
                value += cb * exp(-νb * t)
            end
            ComplexF64(value)
        end
    end

    @testset "Thermal spectrum and convergence" begin
        baths = [
            brownian_oscillator_bath(; parameters..., matsubara = K) for K in (0, 4, 16, 64)
        ]
        for (K, bath) in zip((0, 4, 16, 64), baths)
            @test length(bath.rates) == K + 2
            @test bath.counterterm == λ
            @test bath.diffusion == 0
            @test all(<(0), real.(bath.coefficients[3:end]))
            for ω in (0.03, 0.7, ω₀, 3.9)
                # The commutator, and therefore J, is independent of the
                # truncation of the real Matsubara contribution.
                @test bath_spectrum(bath, ω) - bath_spectrum(bath, -ω) ≈ 2ħ * J(ω) rtol =
                    2e-12 atol = 2e-14
            end
        end
        for ω in (-3.9, -ω₀, -0.1, 0.0, 0.1, ω₀, 3.9)
            errors = [bath_spectrum(b, ω) - thermal_spectrum(ω) for b in baths]
            # Omitted Brownian Matsubara residues are negative: truncation
            # approaches the physical thermal spectrum from above.
            @test all(>(0), errors)
            @test all(<(0), diff(errors))
            @test last(errors) < 3e-8
        end
        converged = brownian_oscillator_bath(; parameters..., matsubara = 128)
        @test bath_spectrum(converged, 0.0) ≈ 4λ * γ * θ / ω₀^2 atol = 4e-9
        for ω in (0.1, 0.7, ω₀, 2.5)
            @test bath_spectrum(converged, -ω) ≈
                  exp(-ħ * ω / θ) * bath_spectrum(converged, ω) atol = 4e-9
            @test bath_spectrum(converged, ω) ≈ thermal_spectrum(ω) atol = 4e-9
        end
    end

    @testset "Independent spectral integrals and time correlation" begin
        # Gauss-Legendre quadrature after ω=ω₀(1+x)/(1-x) maps the
        # semi-infinite spectral integral to [-1,1], with no new dependency.
        quadrature = eigen(SymTridiagonal(zeros(256), [k / sqrt(4k^2 - 1) for k in 1:255]))
        x = quadrature.values
        weights = 2 .* quadrature.vectors[1, :] .^ 2
        frequencies = @. ω₀ * (1 + x) / (1 - x)
        jacobian = @. 2ω₀ / (1 - x)^2
        reorganization_integral =
            sum(weights .* jacobian .* J.(frequencies) ./ frequencies) / π
        @test reorganization_integral ≈ λ rtol = 2e-13
        variance =
            ħ / π *
            sum(weights .* jacobian .* J.(frequencies) .* coth.(ħ .* frequencies ./ (2θ)))
        fine = brownian_oscillator_bath(; parameters..., matsubara = 128)
        coarse = brownian_oscillator_bath(; parameters..., matsubara = 8)
        # Here 2ω₀²>γ², so each omitted |cₖ| is bounded by
        # 4λγω₀²θ/νₖ³ and the integral test bounds the infinite tail.
        variance_tail_bound = 2λ * γ * ω₀^2 * θ / ((2π * θ / ħ)^3 * 128^2)
        @test 0 < real(bath_correlation(fine, 0)) - variance < variance_tail_bound
        @test abs(real(bath_correlation(fine, 0)) - variance) <
              abs(real(bath_correlation(coarse, 0)) - variance) / 100
        @test iszero(imag(bath_correlation(fine, 0)))
        α = γ / 2
        Ω = sqrt(ω₀^2 - α^2)
        for t in (0.0, 0.2, 0.9, 2.5)
            @test imag(bath_correlation(coarse, t)) ≈
                  -λ * ħ * ω₀^2 * exp(-α * t) * sin(Ω * t) / Ω atol = 2e-15
            @test bath_correlation(coarse, -t) ≈ conj(bath_correlation(coarse, t))
        end
        hot = brownian_oscillator_bath(; parameters..., kT = 1e5, matsubara = 0)
        for t in (0.0, 0.2, 0.9)
            classical = 2λ * 1e5 * exp(-α * t) * (cos(Ω * t) + α * sin(Ω * t) / Ω)
            @test real(bath_correlation(hot, t)) ≈ classical rtol = 3e-11
        end
    end

    @testset "Near-critical damping and argument validation" begin
        critical_parameters = (;
            reorganization = 0.2,
            frequency = 1.0,
            damping = prevfloat(2.0),
            kT = 0.7,
            hbar = 0.9,
        )
        nearcritical = brownian_oscillator_bath(; critical_parameters..., matsubara = 4)
        @test all(isfinite, nearcritical.coefficients)
        @test maximum(abs, nearcritical.coefficients) < 1
        for t in (0.0, 0.3, 1.2)
            reference = reference_correlation(critical_parameters, 4, t)
            @test bath_correlation(nearcritical, t) ≈ reference rtol = 2e-13 atol = 2e-15
        end
        isolated =
            brownian_oscillator_bath(; parameters..., reorganization = 0, matsubara = 4)
        @test isempty(isolated.rates)
        @test isolated.hbar == ħ
        @test isolated.counterterm == isolated.diffusion == 0
        @test bath_correlation(isolated, 0.3) == 0
        @test bath_spectrum(isolated, 0.3) == 0
        for keyword in (:frequency, :damping, :kT, :hbar), value in (0.0, -1.0, Inf, NaN)
            invalid = NamedTuple{(keyword,)}((value,))
            @test_throws ArgumentError brownian_oscillator_bath(; parameters..., invalid...)
        end
        for value in (-1.0, Inf, NaN)
            @test_throws ArgumentError brownian_oscillator_bath(;
                parameters...,
                reorganization = value,
            )
        end
        for damping in (2ω₀, 3ω₀)
            @test_throws ArgumentError brownian_oscillator_bath(; parameters..., damping)
        end
        @test_throws ArgumentError brownian_oscillator_bath(; parameters..., matsubara = -1)
        @test_throws ArgumentError brownian_oscillator_bath(;
            parameters...,
            matsubara = typemax(Int),
        )
        @test_throws ArgumentError brownian_oscillator_bath(;
            reorganization = floatmax(),
            frequency = 2,
            damping = 0.3,
            kT = 1,
        )
        @test_throws ArgumentError brownian_oscillator_bath(;
            parameters...,
            kT = nextfloat(0.0),
        )
    end

    @testset "Retained Matsubara pole near critical damping" begin
        coupling, temperature, planck = 0.2, 0.03, 0.8
        for pole in (1, 2), displacement in (-1e-12, 0.0, 1e-12), gap in (eps(), 1e-10)
            α = 2π * pole * temperature / planck * (1 + displacement)
            frequency = α * (1 + gap)
            damping = 2α
            K = pole + 2
            collision_parameters = (;
                reorganization = coupling,
                frequency,
                damping,
                kT = temperature,
                hbar = planck,
            )
            bath = brownian_oscillator_bath(; collision_parameters..., matsubara = K)
            @test all(isfinite, bath.coefficients)
            @test maximum(abs, bath.coefficients) < 1
            for t in (0.0, 0.5, 2.0)
                # Evaluate the singular residues independently at high precision;
                # their cancellation would lose almost every digit in Float64.
                reference = reference_correlation(collision_parameters, K, t)
                @test bath_correlation(bath, t) ≈ reference rtol = 3e-12 atol = 2e-15
            end
            for ω in (0.03, frequency, 0.8)
                density =
                    2coupling * damping * frequency^2 * ω /
                    ((frequency^2 - ω^2)^2 + damping^2 * ω^2)
                @test bath_spectrum(bath, ω) - bath_spectrum(bath, -ω) ≈ 2planck * density rtol =
                    3e-12 atol = 2e-15
            end
        end
        # Probe both sides of the regularized-basis threshold, including a
        # nonzero thermal-pole displacement in either direction.
        for distance in (0.05 * 0.999, 0.05 * 1.001), imaginary in (-0.03, 0.0, 0.03)
            realpart = sqrt(distance^2 - imaginary^2)
            α = 2π * temperature / planck + (2temperature / planck) * imaginary
            Ω = (2temperature / planck) * realpart
            boundary_parameters = (;
                reorganization = coupling,
                frequency = hypot(α, Ω),
                damping = 2α,
                kT = temperature,
                hbar = planck,
            )
            bath = brownian_oscillator_bath(; boundary_parameters..., matsubara = 3)
            @test iszero(bath.weights[3]) == (distance < 0.05)
            for t in (0.0, 0.5, 2.0)
                reference = reference_correlation(boundary_parameters, 3, t)
                @test bath_correlation(bath, t) ≈ reference rtol = 1e-11 atol = 2e-15
            end
        end
    end

    @testset "Brownian bath propagates with scaled and unscaled hierarchies" begin
        grid = PhaseSpaceGrid((-5, 5), 16, (-5, 5), 16)
        bath = brownian_oscillator_bath(; parameters..., matsubara = 1)
        options = (; mass = 1, potential = q -> q^2 / 2, bath, depth = 2)
        unscaled = heom_operator(grid; options...)
        scaled = heom_operator(grid; options..., scaled = true)
        U = zeros(size(grid)..., length(unscaled.indices))
        for a in axes(U, 3)
            U[:, :, a] .= on_grid((q, p) -> cos(a) * exp(-q^2 - p^2) * (1 + 0.1a * p), grid)
        end
        transformed = rescale_hierarchy(U, unscaled; scaled = true)
        du, ds = similar(U), similar(U)
        heom!(du, U, unscaled, 0.0)
        heom!(ds, transformed, scaled, 0.0)
        @test all(isfinite, du)
        @test ds ≈ rescale_hierarchy(du, unscaled; scaled = true) rtol = 2e-14
    end
end
