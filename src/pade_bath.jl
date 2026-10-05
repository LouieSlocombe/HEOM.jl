using LinearAlgebra: SymTridiagonal, eigen

# The [N/N] Padé approximant of the Bose function is
# f(x) = 1/x + 1/2 + Σ 2η[k]x/(x² + ξ[k]²) + R*x.
# A (2N+1)-term truncation of its Stieltjes continued fraction gives the
# Jacobi matrix below. Paired eigenvalues yield ξ = 2/|eigenvalue|; the
# first-component spectral weights give η without products of pole gaps.
# See Hu et al., J. Chem. Phys. 134, 244106 (2011), and Ding et al.,
# J. Chem. Phys. 135, 164107 (2011), Eqs. (5)-(8).
function bose_pade(N::Int)
    R = inv(4.0 * (N + 1) * (2.0N + 3))
    N == 0 && return Float64[], Float64[], R
    offdiagonal = [inv(sqrt((2.0j + 1) * (2.0j + 3))) for j in 1:2N]
    decomposition = eigen(SymTridiagonal(zeros(2N + 1), offdiagonal))
    ξ = -2 ./ decomposition.values[1:N]
    η = [(decomposition.vectors[1, k] * ξ[k])^2 / 12 for k in 1:N]
    return ξ, η, R
end

"""
    drude_lorentz_pade_bath(; reorganization, cutoff, kT, pade = 4,
                            hbar = 1.0, terminator = true)

Thermal Drude–Lorentz [`ExponentialBath`](@ref), using the `[N/N]` Bose Padé
decomposition with `N = pade`. It has the same spectral density, units and
counterterm as [`drude_lorentz_bath`](@ref). Padé poles converge much faster than
Matsubara truncation over the frequency range relevant to low-temperature dynamics.
There is no weak-coupling or high-temperature approximation in the HEOM itself.

The Bose function is approximated consistently at every pole by

    f_N(x) = 1/x + 1/2 + Σₖ 2ηₖ*x/(x² + ξₖ²) + R_N*x,
    R_N = 1/[4(N+1)(2N+3)],   νₖ = ξₖ*kT/ħ.

The Drude coefficient uses `f_N(-im*ħ*cutoff/kT)`, not the exact cotangent with
only some of its cancelling poles retained. With `terminator = true`, the Padé
white-noise remainder adds momentum diffusion
`D = 2λ*γ*ħ²*R_N/kT`. This is part of the finite-order bath approximation and must
be converged by increasing `pade`; it does not terminate the hierarchy depth.
`terminator = false` drops that diffusion but keeps the same exponential terms.

Any finite positive `kT` is allowed. Negative real correlation coefficients at
low temperature are retained. Increase `pade`, hierarchy depth, grid extent and
grid resolution independently; a finite pole count cannot resolve all frequencies
at arbitrarily low temperature. `kT = 0` requires a different decomposition.
Coincident or nearly coincident Drude and Padé poles use a coupled correlation basis,
which represents the polynomial-exponential limit without divergent residues or
changing the physical bath. Zero coupling returns an empty bath.

References: Hu et al., J. Chem. Phys. **134**, 244106 (2011),
doi:10.1063/1.3602466; Ding et al., J. Chem. Phys. **135**, 164107 (2011),
doi:10.1063/1.3653479, https://arxiv.org/abs/1107.0249.
"""
function drude_lorentz_pade_bath(;
    reorganization::Real,
    cutoff::Real,
    kT::Real,
    pade::Integer = 4,
    hbar::Real = 1.0,
    terminator::Bool = true,
)
    λ, γ, θ, ħ = Float64(reorganization), Float64(cutoff), Float64(kT), Float64(hbar)
    isfinite(λ) && λ >= 0 ||
        throw(ArgumentError("reorganization must be finite and nonnegative"))
    isfinite(γ) && γ > 0 || throw(ArgumentError("cutoff must be finite and positive"))
    isfinite(θ) && θ > 0 || throw(ArgumentError("kT must be finite and positive"))
    isfinite(ħ) && ħ > 0 || throw(ArgumentError("hbar must be finite and positive"))
    0 <= pade <= (typemax(Int) - 1) ÷ 2 ||
        throw(ArgumentError("pade must be a nonnegative representable order"))
    iszero(λ) && return ExponentialBath(ComplexF64[], Float64[]; hbar = ħ)
    ξ, η, R = bose_pade(Int(pade))
    ν = [γ; ξ .* (θ / ħ)]
    c = zeros(ComplexF64, length(ν))
    weights = ones(length(ν))
    mixing = SparseArrays.spzeros(length(ν), length(ν))
    remainder = 2λ * γ * ħ * (ħ / θ) * R
    drude = 2λ * θ - γ * remainder
    for k in eachindex(η)
        v = ν[k+1]
        if abs(v - γ) <= 1e-3 * max(v, γ)
            # Change from two cancelling exponentials to their divided
            # difference. At equal rates exp(-Γt) contains -t*exp(-γt).
            drude += 4λ * θ * η[k] * (γ / (v + γ))
            c[k+1] = 4λ * γ * θ * η[k] * (v / (v + γ))
            weights[k+1] = 0
            mixing[1, k+1] = 1
        else
            # Factor the denominator to avoid squaring dimensional frequencies.
            ratio = γ / v
            c[k+1] = 4λ * θ * η[k] * ratio / ((1 - ratio) * (1 + ratio))
            drude -= real(c[k+1]) * ratio
        end
    end
    c[1] = ComplexF64(drude, -λ * ħ * γ)
    return ExponentialBath(
        c,
        ν;
        hbar = ħ,
        diffusion = terminator ? remainder : 0.0,
        counterterm = λ,
        weights,
        mixing,
    )
end
