using LinearAlgebra: eigvals!

# Multipliers of the hierarchy's first and second momentum derivatives on the
# nonnegative Fourier modes of p. Both operators are real, so the negative modes give
# complex-conjugate multipliers and need not be examined.
function momentum_multipliers(op::WignerHEOM)
    np = length(op.grid.p)
    if op.hamiltonian isa SpectralWignerMoyal
        return vec(op.momentum.first) .* np, vec(op.momentum.second) .* np
    end
    # A periodic difference matrix is circulant: its multipliers are the DFT of a row.
    function multipliers(D)
        columns, weights = SparseArrays.findnz(D[1, :])
        return [
            sum(w * cis(2π * (c - 1) * j / np) for (c, w) in zip(columns, weights)) for
            j in 0:(np÷2)
        ]
    end
    return multipliers(op.momentum.first), multipliers(op.momentum.second)
end

# Hierarchy block of the generator at a fixed position q and momentum wavenumber, where
# ∂p and ∂p² act as the multipliers `first` and `second`.
function frozen_hierarchy!(H, op::WignerHEOM, q, first, second)
    b = op.bath
    H .= op.transfer
    for a in eachindex(op.indices)
        H[a, a] += b.diffusion * second - op.damping[a]
        for k in eachindex(b.rates)
            below, above = op.lower[k, a], op.upper[k, a]
            below != 0 && (H[below, a] += op.lowering[k, a] * first)
            above != 0 && (
                H[above, a] +=
                    op.raising_derivative[k, a] * first + op.raising_coordinate[k, a] * q
            )
        end
    end
    return H
end

"""
    hierarchy_stability(op::WignerHEOM)

Frozen-coefficient stability indicator for the truncated hierarchy of `op`.

Every hierarchy coupling is local in position `q` and in the momentum wavenumber `κ`
conjugate to `p`: `∂p` multiplies by `iκ` (or the finite-difference symbol) and the
coordinate coupling multiplies by `q`. The Moyal potential term adds the same imaginary
multiplier to every member. Only the kinetic term `-(p/m)∂q` couples different `(q, κ)`.
Freezing `(q, κ)` therefore leaves one small `members × members` matrix per grid point.
This function returns the largest real part of its eigenvalues over the grid as a
named tuple `(rate, position, wavenumber, radius)`:

  - `rate`: the largest frozen growth rate over every grid position and nonnegative
    discrete wavenumber, in units of inverse time;
  - `position`, `wavenumber`: the `|q|` and `κ` at which it occurs;
  - `radius`: the smallest `|q|` at which any wavenumber is unstable, or `Inf`.

For a bath with nonnegative zero-frequency noise, the exact, untruncated hierarchy is
stable at every `(q, κ)`. A finite depth represents the bath faithfully only while
`|q|·κ` stays small compared with the depth and the bath rates. Beyond `radius` the
truncated hierarchy has growing modes, and the rate increases with `|q|`. On the periodic
box these modes concentrate near the edges `|q| ≈ L`. Rounding and the tails of the
physical state seed them, and they eventually overwhelm the state even when it never
approaches the edge. See the stability notes in [`heom_operator`](@ref).

`rate` is a heuristic, not a proven bound. In the package's tests and in the
`prototypes/box_edge` study, it was never below the spectral abscissa of the complete
generator, but it can lie far above it. The kinetic term carries modes through the
unstable region, so a box that extends only a little beyond `radius` can still be
stable. A rate of zero, to rounding, means that no local mode grows. A positive rate
means that the hard cutoff may diverge, as it did in most cases studied. The
divergence then develops over times of order `1/rate`. In the cases studied, the
root's second moments were wrong by 10⁻³ after 6 to 20 multiples of `1/rate`. Check
results beyond a few multiples against a different box and depth.

The cost is about `nq·np/4` dense eigenvalue problems, each of size
`length(hierarchy_indices(op))`. Drives and potentials do not affect the result.
"""
function hierarchy_stability(op::WignerHEOM)
    first, second = momentum_multipliers(op)
    # (q, κ) and (-q, -κ) are similar under a sign change of odd tiers, and -κ gives the
    # complex conjugate, so |q| and κ ≥ 0 cover every grid point.
    positions = sort!(unique(abs.(op.grid.q)))
    members = length(op.indices)
    H = zeros(ComplexF64, members, members)
    tolerance = sqrt(eps()) * max(1.0, maximum(op.damping))
    rate, position, wavenumber, radius = -Inf, NaN, NaN, Inf
    for q in positions, j in eachindex(first)
        growth = maximum(real, eigvals!(frozen_hierarchy!(H, op, q, first[j], second[j])))
        if growth > rate
            rate, position, wavenumber =
                growth, q, 2π * (j - 1) / (length(op.grid.p) * op.grid.dp)
        end
        growth > tolerance && (radius = min(radius, q))
    end
    # The physical root always has a zero eigenvalue at κ = 0.
    rate > tolerance || return (; rate = 0.0, position = NaN, wavenumber = NaN, radius)
    return (; rate, position, wavenumber, radius)
end
