"""
    bath_correlation(bath::ExponentialBath, t::Real)

Evaluate the regular force correlation `C(t) = ⟨B(t)B(0)⟩` of a bath. For `t ≥ 0`,
this is `transpose(weights) * exp(-Γ*t) * coefficients`, where
`Γ = Diagonal(rates) + mixing`. Negative times use `C(-t) = conj(C(t))`.
The time must be finite and representable as a `Float64`.

The white-noise contribution `2bath.diffusion*δ(t)` is excluded, including at `t = 0`;
it is included in [`bath_spectrum`](@ref). The counterterm does not contribute to
the force correlation. This evaluates the retained bath decomposition, including its
weights and mixing. An empty bath has zero regular correlation.

Use broadcasting for sampled times: `bath_correlation.(Ref(bath), times)`.
"""
function bath_correlation(bath::ExponentialBath, t::Real)
    time = Float64(t)
    isfinite(time) || throw(ArgumentError("time must be finite in Float64"))
    τ = abs(time)
    value = if iszero(τ)
        dot(bath.weights, bath.coefficients)
    elseif iszero(bath.mixing)
        correlation = 0.0 + 0.0im
        for k in eachindex(bath.rates)
            correlation += bath.weights[k] * bath.coefficients[k] * exp(-bath.rates[k] * τ)
        end
        correlation
    else
        mixed_bath_correlation(bath, τ)
    end
    return time < 0 ? conj(value) : value
end

function mixed_bath_correlation(bath, t)
    # Most structured baths couple only a few oscillator modes. Keep independent
    # thermal poles out of the dense matrix exponential, even for long tails.
    coupled = coupled_bath_modes(bath.mixing)
    indices = findall(coupled)
    generator = Matrix(bath.mixing[indices, indices])
    for (j, k) in enumerate(indices)
        generator[j, j] = bath.rates[k]
    end
    value = dot(bath.weights[indices], exp(-t * generator) * bath.coefficients[indices])
    for k in eachindex(bath.rates)
        coupled[k] && continue
        value += bath.weights[k] * bath.coefficients[k] * exp(-bath.rates[k] * t)
    end
    return value
end

"""
    bath_spectrum(bath::ExponentialBath, omega::Real)

Evaluate the two-sided, unsymmetrized force spectrum at angular frequency `omega`,
with Fourier convention `S(ω) = ∫₋∞∞ exp(iωt) C(t) dt`. Frequencies may be positive,
zero or negative, and must be finite and representable as a `Float64`.

The spectrum of the retained bath decomposition is
`S(ω) = 2real(transpose(weights) * ((Γ - iωI) \\ coefficients)) + 2bath.diffusion`,
where `Γ = Diagonal(rates) + mixing`. Thus the white-noise correlation
`2bath.diffusion*δ(t)` adds a frequency-independent `2bath.diffusion`, while the
counterterm does not contribute. No symmetrization or clipping is applied.

For an exact thermal bath with the spectral-density convention used by
[`drude_lorentz_bath`](@ref), `S(ω) = 2ħJ(ω)/(1 - exp(-ħω/kT))` for `ω > 0`, and
`S(-ω) = exp(-ħω/kT)*S(ω)`. Finite decompositions approximate this thermal relation.
Use broadcasting for sampled frequencies: `bath_spectrum.(Ref(bath), frequencies)`.
"""
function bath_spectrum(bath::ExponentialBath, omega::Real)
    frequency = Float64(omega)
    isfinite(frequency) || throw(ArgumentError("frequency must be finite in Float64"))
    integral = if iszero(bath.mixing)
        value = 0.0 + 0.0im
        for k in eachindex(bath.rates)
            value +=
                bath.weights[k] * bath.coefficients[k] / (bath.rates[k] - im * frequency)
        end
        value
    else
        generator = bath.mixing + spdiagm(0 => bath.rates .- im * frequency)
        dot(bath.weights, generator \ bath.coefficients)
    end
    return 2real(integral) + 2bath.diffusion
end
