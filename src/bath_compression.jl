using LinearAlgebra: Diagonal, diag, lyap, schur, svd

# For t ≥ 0, C(t) = wᵀexp(-Γt)c is the impulse response of a real linear system:
# state matrix -Γ, input columns real(c) and imag(c), and output wᵀ.
function bath_realization(bath)
    Γ = Matrix(bath.mixing)
    for k in eachindex(bath.rates)
        Γ[k, k] = bath.rates[k]
    end
    return Γ, [real.(bath.coefficients) imag.(bath.coefficients)], bath.weights
end

# Factor a positive-semidefinite Gramian as L*L'. Rounding can leave eigenvalues
# slightly negative, where a Cholesky factorization would fail; clip them.
function gramian_factor(X)
    F = eigen(Symmetric(X))
    return F.vectors * Diagonal(sqrt.(max.(F.values, 0)))
end

# Square-root balancing: the SVD of the product of Gramian factors gives the
# Hankel singular values and balancing projections without inverting a Gramian.
function balanced_factors(bath)
    Γ, inputs, outputs = bath_realization(bath)
    controllability = gramian_factor(lyap(-Γ, inputs * inputs'))
    observability = gramian_factor(lyap(-transpose(Γ), outputs * outputs'))
    return Γ, controllability, observability, svd(observability' * controllability)
end

"""
    hankel_singular_values(bath::ExponentialBath)

Hankel singular values `σ₁ ≥ σ₂ ≥ … ≥ 0` of the retained exponential correlation,
one per bath mode. For `t ≥ 0`, `C(t) = transpose(weights)*exp(-Γ*t)*coefficients`
with `Γ = Diagonal(rates) + mixing` is the impulse response of a real linear system
with state matrix `-Γ`, input columns `real(coefficients)` and `imag(coefficients)`,
and output `transpose(weights)`. `σᵢ` measures how strongly the `i`th balanced state
is both excited by, and visible in, that correlation. The diffusion and counterterm
are not included. An empty bath returns an empty vector.

[`compress_bath`](@ref) retains the `r` leading balanced states. Its spectrum then
differs from that of `bath` by at most `4√2 * sum(σ[r+1:end])` at every frequency.
Values below about `sqrt(eps())*σ₁` are not resolved in Float64.
"""
function hankel_singular_values(bath::ExponentialBath)
    isempty(bath.rates) && return Float64[]
    return balanced_factors(bath)[4].S
end

"""
    compress_bath(bath::ExponentialBath; modes, method = :truncate)

Reduce an [`ExponentialBath`](@ref) to `modes` real hierarchy modes by balanced
model-order reduction of its correlation. This cuts the hierarchy from
`binomial(length(bath.rates) + depth, depth)` members to `binomial(modes + depth, depth)`.
For example, a depth-8 hierarchy for 8 modes has 12870 members, and for 4 modes 495.
Every returned mode is real, so `modes` is also the hierarchy mode count; a damped
oscillatory pair uses two. Takahashi and Tanimura, J. Chem. Phys. **158**, 044115 (2023),
Appendix B, apply balanced truncation to Padé decompositions in HEOM in this way.

The correlation is the impulse response of the real system described in
[`hankel_singular_values`](@ref). The `modes` leading balanced states (Moore,
IEEE Trans. Autom. Control **26**, 17 (1981)) are retained. With `r = modes`, both
methods satisfy the frequency-uniform bound

    |bath_spectrum(bath, ω) - bath_spectrum(compressed, ω)| ≤ 4√2 * sum(σ[r+1:end]),

from the H∞ error bound `2∑σᵢ` (Glover, Int. J. Control **39**, 1115 (1984)). This
bounds the change from `bath`, not the error of `bath` itself.

  - `method = :truncate` (balanced truncation) discards the eliminated states. The
    diffusion and counterterm are copied unchanged, so the static balance between the
    counterterm and the retained imaginary correlation is only approximate.
  - `method = :residualize` (the singular-perturbation approximation of Liu and Anderson,
    Int. J. Control **50**, 1379 (1989)) replaces the eliminated states by their
    adiabatic limit. The compressed bath has exactly the same `∫₀^∞ C(t) dt`. The
    eliminated real part is added to `diffusion`, as a Markovian white-noise remainder
    like the terminators of [`drude_lorentz_bath`](@ref). The eliminated imaginary part
    `κ` is the static potential `(κ/ħ)*q²`, so it is added to `counterterm`. Thus the
    zero-frequency noise and total static force constant are unchanged, which can help
    when few modes are retained. Throws if the eliminated real part is negative beyond
    rounding, which would require negative diffusion. For example, this happens when
    negative Matsubara terms of [`brownian_oscillator_bath`](@ref) are eliminated.

The reduced rate matrix is returned in real Schur form. `rates` are the real parts of
its eigenvalues. A complex-conjugate pair is a 2×2 block with equal rates, and `mixing`
is otherwise strictly upper triangular. The root of a hierarchy truncated at total
depth does not depend on this choice of real basis.

A small spectral or correlation error does not establish accurate thermalization
(Tokieda, Phys. Rev. Research **7**, 043178 (2025)). Compare
[`bath_spectrum`](@ref) with the thermal target over the frequencies the system
resolves, and [`harmonic_covariance`](@ref) of both baths with an exact equilibrium
when one is available. Then repeat the hierarchy-depth convergence with the
compressed bath.

`modes` must lie between zero and `length(bath.rates)`; the full count returns an
unchanged copy. `modes = 0` keeps only the diffusion and counterterm. Throws if
`σ[modes]` and `σ[modes+1]` are not separated by more than `sqrt(eps())*σ₁`, since the
retained balanced subspace is then not determined.
"""
function compress_bath(bath::ExponentialBath; modes::Integer, method::Symbol = :truncate)
    method in (:residualize, :truncate) ||
        throw(ArgumentError("method must be :residualize or :truncate"))
    n = length(bath.rates)
    0 <= modes <= n || throw(ArgumentError("modes must be between 0 and $n"))
    r, ħ = Int(modes), bath.hbar
    diffusion, counterterm = bath.diffusion, bath.counterterm
    r == n && return ExponentialBath(
        bath.coefficients,
        bath.rates;
        hbar = ħ,
        diffusion,
        counterterm,
        bath.weights,
        bath.mixing,
    )
    Γ, controllability, observability, F = balanced_factors(bath)
    σ = F.S
    if r > 0
        σ[r] - σ[r+1] > sqrt(eps()) * σ[1] || throw(
            ArgumentError(
                "Hankel singular values $r and $(r + 1) are not separated beyond " *
                "numerical resolution; choose a different number of modes",
            ),
        )
    end
    c, w = bath.coefficients, bath.weights
    scale = Diagonal(inv.(sqrt.(σ[1:r])))
    right = controllability * F.V[:, 1:r] * scale
    left = scale * F.U[:, 1:r]' * observability'
    if method === :truncate
        reduced, coefficients, weights = left * Γ * right, left * c, right' * w
    else
        # Truncating the reciprocal system Γ⁻¹ in the same balanced coordinates and
        # inverting back is the singular-perturbation approximation. It needs only
        # the retained projections, never the ill-conditioned eliminated states.
        slow = Γ \ right
        static = left * (Γ \ c)
        reciprocal = left * slow
        reduced = inv(reciprocal)
        coefficients = reciprocal \ static
        observed = slow' * w
        weights = reduced' * observed
        total = dot(w, Γ \ c)
        retained = dot(observed, coefficients)
        remainder = total - retained
        rounding = 64eps() * (abs(total) + abs(retained))
        diffusion += real(remainder)
        counterterm += imag(remainder) / ħ
        -rounding <= diffusion < 0 && (diffusion = 0.0)
        -rounding / ħ <= counterterm < 0 && (counterterm = 0.0)
        diffusion >= 0 || throw(
            ArgumentError(
                "residualization requires negative diffusion $diffusion; the " *
                "eliminated real correlation has negative weight. Retain more " *
                "modes or use method = :truncate",
            ),
        )
        counterterm >= 0 || throw(
            ArgumentError(
                "residualization requires a negative counterterm $counterterm; " *
                "retain more modes or use method = :truncate",
            ),
        )
    end
    r == 0 &&
        return ExponentialBath(ComplexF64[], Float64[]; hbar = ħ, diffusion, counterterm)
    # A stable reduction has positive real Schur diagonals; ExponentialBath checks them.
    form = schur(reduced)
    rates = diag(form.T)
    mixing = copy(form.T)
    for k in eachindex(rates)
        mixing[k, k] = 0
    end
    return ExponentialBath(
        form.Z' * coefficients,
        rates;
        hbar = ħ,
        diffusion,
        counterterm,
        weights = form.Z' * weights,
        mixing = sparse(mixing),
    )
end
