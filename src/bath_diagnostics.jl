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

"""
    harmonic_covariance(bath::ExponentialBath; mass, omega)

Exact reduced equilibrium covariance `[⟨q²⟩ ⟨qp⟩; ⟨qp⟩ ⟨p²⟩]` (symmetrized) of the
harmonic oscillator `V(q) = mass*omega^2*q^2/2` coupled to the retained decomposition
of `bath`. Its correlation, diffusion and counterterm are included exactly, as in an
infinitely deep hierarchy. The covariance is that of the stationary Wigner function:
compare it with `phase_space_covariance(W, grid)` of an equilibrated root. No phase-space
grid or hierarchy is used.

The oscillator obeys a generalized quantum Langevin equation, where `weights'x` is the
bath force on `p`. The vector `x` decays with `Γ = Diagonal(rates) + mixing` and responds
to `q` with `-2imag(coefficients)/ħ`. Its noise has a signed covariance `S` with
`S*weights = real(coefficients)`, which reproduces the real correlation exactly. The
joint stationary covariance solves a Lyapunov equation.

This is the harmonic test of a bath decomposition proposed by Tokieda, Phys. Rev.
Research **7**, 043178 (2025). Compare the result for a compressed or fitted bath with
its source, and with an exact continuum fluctuation–dissipation result when one is
available. A finite decomposition need not produce a positive or uncertainty-respecting
covariance. Throws if the coupled linear system has no unique stationary state.
`mass` and `omega` must be finite and positive.
"""
function harmonic_covariance(bath::ExponentialBath; mass::Real, omega::Real)
    m, ω = Float64(mass), Float64(omega)
    isfinite(m) && m > 0 || throw(ArgumentError("mass must be finite and positive"))
    isfinite(ω) && ω > 0 || throw(ArgumentError("omega must be finite and positive"))
    n = length(bath.rates)
    A, Q = zeros(n + 2, n + 2), zeros(n + 2, n + 2)
    A[1, 2] = 1 / m
    A[2, 1] = -m * ω^2 - 2bath.counterterm
    Q[2, 2] = 2bath.diffusion
    if n > 0
        Γ, _, w = bath_realization(bath)
        residues = real.(bath.coefficients)
        # The symmetric solution of S*w = residues with the smallest Frobenius norm.
        norm_squared = dot(w, w)
        S = if iszero(norm_squared)
            zeros(n, n)
        else
            (residues * w' + w * residues') / norm_squared -
            dot(w, residues) * (w * w') / norm_squared^2
        end
        A[2, 3:end] = w
        A[3:end, 1] = -2imag.(bath.coefficients) / bath.hbar
        A[3:end, 3:end] = -Γ
        Q[3:end, 3:end] = Γ * S + S * Γ'
    end
    all(z -> real(z) < 0, eigvals(A)) ||
        throw(ArgumentError("the coupled oscillator and bath decomposition are not stable"))
    covariance = lyap(A, Q)
    return Symmetric(covariance[1:2, 1:2])
end
