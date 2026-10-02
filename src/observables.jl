"""
    phase_space_integral(W, grid::PhaseSpaceGrid)

Integral `∫∫ W(q, p) dq dp` of a function sampled on `grid`, by the rectangle rule. This
converges spectrally for smooth functions that decay at the edges of the periodic box. The
integral of a Wigner function is the norm of its state.
"""
function phase_space_integral(W::AbstractMatrix, grid::PhaseSpaceGrid)
    check_size(W, grid)
    return sum(W) * grid.dq * grid.dp
end

"""
    expectation(f, W, grid::PhaseSpaceGrid)

Phase-space average `∫∫ f(q, p) W(q, p) dq dp`. For a Wigner function `W` this is the
expectation value of the operator whose Weyl symbol is `f`. For example,
`expectation((q, p) -> p^2 / 2m + V(q), W, grid)` is the mean energy; see also
[`energy`](@ref). No normalisation is applied.
"""
function expectation(f, W::AbstractMatrix, grid::PhaseSpaceGrid)
    check_size(W, grid)
    return phase_space_integral(on_grid(f, grid) .* W, grid)
end

"""
    purity(W, grid::PhaseSpaceGrid; hbar = 1.0)

Purity `Tr ρ² = 2πħ ∫∫ W² dq dp` of the state with Wigner function `W`. It is `1` for a pure
state and less than `1` for a mixed state. The linear entropy is
`1 - purity(W, grid; hbar)`.
"""
function purity(W::AbstractMatrix, grid::PhaseSpaceGrid; hbar::Real = 1.0)
    return overlap(W, W, grid; hbar)
end

"""
    position_density(W, grid::PhaseSpaceGrid)

Position density `⟨q|ρ|q⟩ = ∫ W(q, p) dp` at each position point of `grid`, returned as a
vector. No normalisation is applied.
"""
function position_density(W::AbstractMatrix, grid::PhaseSpaceGrid)
    check_size(W, grid)
    return vec(sum(W; dims = 2)) * grid.dp
end

"""
    momentum_density(W, grid::PhaseSpaceGrid)

Momentum density `⟨p|ρ|p⟩ = ∫ W(q, p) dq` at each momentum point of `grid`, returned as a
vector. No normalisation is applied.
"""
function momentum_density(W::AbstractMatrix, grid::PhaseSpaceGrid)
    check_size(W, grid)
    return vec(sum(W; dims = 1)) * grid.dq
end

"""
    phase_space_mean(W, grid::PhaseSpaceGrid)

Return `[⟨q⟩, ⟨p⟩]`, where `⟨f⟩ = ∫∫ f W dq dp`. The moments are not divided by the norm
of `W`.
"""
function phase_space_mean(W::AbstractMatrix, grid::PhaseSpaceGrid)
    check_size(W, grid)
    q = dot(grid.q, position_density(W, grid)) * grid.dq
    p = dot(grid.p, momentum_density(W, grid)) * grid.dp
    return [q, p]
end

"""
    phase_space_covariance(W, grid::PhaseSpaceGrid)

Return the symmetrised covariance matrix `Σᵢⱼ = ⟨(xᵢ - ⟨xᵢ⟩)(xⱼ - ⟨xⱼ⟩)⟩`, with
`x = (q, p)`. The off-diagonal entry is the covariance of `(qp + pq)/2`. All moments are
integrals against `W`, without renormalisation.
"""
function phase_space_covariance(W::AbstractMatrix, grid::PhaseSpaceGrid)
    check_size(W, grid)
    q, p = phase_space_mean(W, grid)
    δq, δp = grid.q .- q, grid.p .- p
    var_q = dot(δq .^ 2, position_density(W, grid)) * grid.dq
    var_p = dot(δp .^ 2, momentum_density(W, grid)) * grid.dp
    cov_qp = dot(δq, W, δp) * grid.dq * grid.dp
    return [var_q cov_qp; cov_qp var_p]
end

"""
    energy(W, grid::PhaseSpaceGrid; mass, potential)

Mean energy `⟨H⟩ = ∫ p²/(2m) momentum_density(p) dp + ∫ V(q) position_density(q) dq`,
where `potential(q)` evaluates `V(q)`. The Weyl symbol of the Hamiltonian is exactly
`p²/(2m) + V(q)`. No normalisation is applied.
"""
function energy(W::AbstractMatrix, grid::PhaseSpaceGrid; mass::Real, potential)
    check_size(W, grid)
    kinetic = dot(grid.p .^ 2, momentum_density(W, grid)) * grid.dp / (2mass)
    potential_energy = dot(potential.(grid.q), position_density(W, grid)) * grid.dq
    return kinetic + potential_energy
end

"""
    overlap(W1, W2, grid::PhaseSpaceGrid; hbar = 1.0)

State overlap `Tr ρ₁ρ₂ = 2πħ ∫∫ W₁ W₂ dq dp`. For pure states this is
`|⟨ψ₁|ψ₂⟩|²`; when one state is an eigenstate, it gives that state's population. Neither
Wigner function is renormalised.
"""
function overlap(
    W1::AbstractMatrix,
    W2::AbstractMatrix,
    grid::PhaseSpaceGrid;
    hbar::Real = 1.0,
)
    check_size(W1, grid)
    check_size(W2, grid)
    return 2π * hbar * phase_space_integral(W1 .* W2, grid)
end

"""
    wigner_negativity(W, grid::PhaseSpaceGrid)

Negative volume `∫∫ (abs(W) - W) dq dp`. For a normalised Wigner function this is the
Kenfack–Życzkowski `δ`, and equals `-2N_W` in Cabrera's convention. It is exactly zero for
non-negative `W`. The kink on the zero contour makes this quadrature only second-order
accurate. No normalisation is applied.
"""
function wigner_negativity(W::AbstractMatrix, grid::PhaseSpaceGrid)
    check_size(W, grid)
    return sum(w -> abs(w) - w, W) * grid.dq * grid.dp
end
