"""
    phase_space_integral(W, grid::PhaseSpaceGrid)

Integral `∫∫ W(q, p) dq dp` of a function sampled on `grid`, by the rectangle rule. This
converges spectrally for smooth functions that decay at the edges of the periodic box. The
integral of a Wigner function is the norm of its state.
"""
function phase_space_integral(W::AbstractMatrix, grid::PhaseSpaceGrid)
    size(W) == size(grid) || throw(
        DimensionMismatch(
            "matrix has size $(size(W)), but the grid has size $(size(grid))",
        ),
    )
    return sum(W) * grid.dq * grid.dp
end

"""
    expectation(f, W, grid::PhaseSpaceGrid)

Phase-space average `∫∫ f(q, p) W(q, p) dq dp`. For a Wigner function `W` this is the
expectation value of the operator whose Weyl symbol is `f`. For example,
`expectation((q, p) -> p^2 / 2m + V(q), W, grid)` is the mean energy.
"""
function expectation(f, W::AbstractMatrix, grid::PhaseSpaceGrid)
    return phase_space_integral(on_grid(f, grid) .* W, grid)
end

"""
    purity(W, grid::PhaseSpaceGrid; hbar = 1.0)

Purity `Tr ρ² = 2πħ ∫∫ W² dq dp` of the state with Wigner function `W`. It is `1` for a pure
state and less than `1` for a mixed state.
"""
function purity(W::AbstractMatrix, grid::PhaseSpaceGrid; hbar::Real = 1.0)
    return 2π * hbar * phase_space_integral(W .^ 2, grid)
end
