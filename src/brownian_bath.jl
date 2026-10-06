"""
    brownian_oscillator_bath(; reorganization, frequency, damping, kT,
                              matsubara = 0, hbar = 1.0)

Underdamped thermal [`ExponentialBath`](@ref) with spectral density

    J(ω) = 2λ*γ*ω₀²*ω / ((ω₀² - ω²)² + γ²*ω²),

where `λ = reorganization`, `ω₀ = frequency` and `γ = damping`. Frequencies are
angular frequencies; the oscillation decays at `γ/2`, with damped frequency
`Ω = sqrt(ω₀² - γ²/4)`. The correlation convention is the same as
[`drude_lorentz_bath`](@ref):
`C(t) = (ħ/π)∫₀∞ J(ω)[coth(ħω/(2kT))*cos(ωt) - im*sin(ωt)]dω`.
The counterterm is `λ*q²`; for dimensional position coupling, `λ` has units of
energy / position².

Two coupled real basis modes represent the damped oscillation. For
`z = ħ*(Ω + im*γ/2)/(2kT)` and `A = λ*ħ*ω₀²/Ω`, their correlation is

    Cosc(t) = A*exp(-γ*t/2) *
              (real(coth(z))*cos(Ω*t) - (imag(coth(z)) + im)*sin(Ω*t)).

The additional `matsubara = K` terms have

    νₖ = 2π*k*kT/ħ,
    cₖ = -4λ*γ*ω₀²*kT*νₖ / ((ω₀² + νₖ²)² - γ²*νₖ²),   k = 1,…,K.

These real coefficients are negative and are retained without clipping. The
omitted tail is dropped: it has negative integrated weight and cannot be replaced
by a nonnegative momentum-diffusion terminator. Increase `matsubara` to converge
[`bath_correlation`](@ref) and [`bath_spectrum`](@ref), especially at low temperature,
and converge hierarchy depth independently. A finite expansion need not obey
thermal detailed balance exactly.

Near critical damping and a retained Matsubara pole, a coupled three-mode basis
avoids cancellation between large residues without changing the correlation.

`λ` must be finite and nonnegative; `frequency`, `damping`, `kT` and `hbar` must
be finite and positive, with `0 < damping < 2frequency`. Critical/overdamped and
zero-temperature limits are not included. Zero coupling returns an empty bath;
otherwise the bath has `K + 2` modes and zero residual diffusion. Combine independent
components coupled to the same position with [`combine_baths`](@ref).

Reference: Lambert et al., Nature Communications **10**, 3721 (2019),
https://doi.org/10.1038/s41467-019-11656-1 (with coupling converted to `λ`).
"""
function brownian_oscillator_bath(;
    reorganization::Real,
    frequency::Real,
    damping::Real,
    kT::Real,
    matsubara::Integer = 0,
    hbar::Real = 1.0,
)
    λ, ω₀, γ = Float64(reorganization), Float64(frequency), Float64(damping)
    θ, ħ = Float64(kT), Float64(hbar)
    isfinite(λ) && λ >= 0 ||
        throw(ArgumentError("reorganization must be finite and nonnegative"))
    isfinite(ω₀) && ω₀ > 0 || throw(ArgumentError("frequency must be finite and positive"))
    isfinite(γ) && γ > 0 || throw(ArgumentError("damping must be finite and positive"))
    isfinite(θ) && θ > 0 || throw(ArgumentError("kT must be finite and positive"))
    isfinite(ħ) && ħ > 0 || throw(ArgumentError("hbar must be finite and positive"))
    α = γ / 2
    0 < α < ω₀ || throw(ArgumentError("require 0 < damping < 2frequency"))
    0 <= matsubara <= typemax(Int) - 2 ||
        throw(ArgumentError("matsubara must be a nonnegative representable order"))
    iszero(λ) && return ExponentialBath(ComplexF64[], Float64[]; hbar = ħ)

    # Factored frequency differences remain accurate close to critical damping.
    Ω = sqrt(ω₀ - α) * sqrt(ω₀) * sqrt(1 + α / ω₀)
    z = (ħ / θ / 2) * complex(Ω, α)
    isfinite(z) && real(z) > 0 ||
        throw(ArgumentError("the thermal oscillator frequency must be finite and positive"))
    thermal = coth(z)
    amplitude = λ * ħ * ω₀
    K = Int(matsubara)
    rates = [α; α; [2π * k * (θ / ħ) for k in 1:K]]
    coefficients = zeros(ComplexF64, K + 2)
    coefficients[1] = amplitude * (real(thermal) / (Ω / ω₀))
    coefficients[2] = complex(-amplitude * imag(thermal), -amplitude)
    weights = ones(K + 2)
    weights[2] = 0
    mixing = SparseArrays.spzeros(K + 2, K + 2)
    # The second basis function is ω₀*sin(Ω*t)/Ω. Unlike a pure rotation
    # basis, its initial coefficient stays finite as Ω tends to zero.
    mixing[1, 2] = -ω₀
    mixing[2, 1] = (ω₀ - α) * (1 + α / ω₀)

    # As Ω → 0 and α → νₖ, three poles coalesce. Transfer the cancelling
    # Matsubara mode into the oscillator block before evaluating its residue.
    pole = round(imag(z) / π)
    resonant = 0
    if 1 <= pole <= K
        k = Int(pole) + 2
        ν = rates[k]
        δ = ν - α
        ε = (ħ / θ / 2) * complex(Ω, -δ)
        if abs(ε) < 0.05
            resonant = k
            # The new third basis is
            # [exp(-νt) - exp(-αt)*(cos(Ωt) - δ*sin(Ωt)/Ω)] / (δ²+Ω²).
            # It tends to t²*exp(-αt)/2 at the triple pole. Transforming
            # (c₁+cₖ, c₂-δ*cₖ/ω₀, (δ²+Ω²)*cₖ) analytically avoids cancellation.
            # coth(ε) - 1/ε, evaluated through its analytic continuation.
            regular =
                ε * (
                    1 / 3 +
                    ε^2 * (-1 / 45 + ε^2 * (2 / 945 + ε^2 * (-1 / 4725 + ε^2 * 2 / 93555)))
                )
            scale = max(ν, ω₀)
            w = ω₀ / scale
            E = ((ν / scale) + (α / scale))^2 + (Ω / scale)^2
            coefficients[1] = 2λ * θ * w^2 / E + amplitude * (real(regular) / (Ω / ω₀))
            coefficients[2] = complex(
                -2λ * θ * w * (δ / scale) / E - amplitude * imag(regular),
                -amplitude,
            )
            coefficients[k] = -4λ * γ * θ * ν * w^2 / E
            weights[k] = 0
            mixing[2, k] = -inv(ω₀)
        end
    end
    for k in 3:length(rates)
        k == resonant && continue
        ν = rates[k]
        scale = max(ν, ω₀)
        v, a, o, w = ν / scale, α / scale, Ω / scale, ω₀ / scale
        # ((ν-α)²+Ω²)((ν+α)²+Ω²), scaled to avoid fourth powers
        # of dimensional frequencies and cancellation near critical damping.
        denominator = (hypot(v - a, o) * hypot(v + a, o))^2
        coefficients[k] = -4λ * θ * (γ / scale) * w^2 * v / denominator
    end
    return ExponentialBath(coefficients, rates; hbar = ħ, counterterm = λ, weights, mixing)
end
