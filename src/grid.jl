"""
    PhaseSpaceGrid(qlims, nq, plims, np)

Uniform periodic grid on the phase-space box `[qlims[1], qlims[2]) × [plims[1], plims[2])`
with `nq` position points and `np` momentum points. The right endpoints are excluded, so
`q[i] = qlims[1] + (i - 1) * dq` with `dq = (qlims[2] - qlims[1]) / nq`, and likewise for
`p`.

A function on the grid is an `nq × np` matrix `W` with `W[i, j] = W(q[i], p[j])`; see
[`on_grid`](@ref). Both discretisations treat the box as periodic, so a Wigner function must
decay to zero at its edges.

# Examples

```julia
grid = PhaseSpaceGrid((-8.0, 8.0), 64, (-8.0, 8.0), 64)
size(grid) # (64, 64)
grid.dq # 0.25
```
"""
struct PhaseSpaceGrid
    q::Vector{Float64}
    p::Vector{Float64}
    dq::Float64
    dp::Float64
    function PhaseSpaceGrid(
        qlims::Tuple{Real,Real},
        nq::Integer,
        plims::Tuple{Real,Real},
        np::Integer,
    )
        q, dq = periodic_points(qlims, nq)
        p, dp = periodic_points(plims, np)
        return new(q, p, dq, dp)
    end
end

# The `n` points and spacing of a periodic grid on `[lims[1], lims[2])`.
function periodic_points(lims::Tuple{Real,Real}, n::Integer)
    lo, hi = Float64.(lims)
    n >= 2 || throw(ArgumentError("a grid needs at least 2 points per axis, got $n"))
    isfinite(lo) && isfinite(hi) && hi > lo ||
        throw(ArgumentError("grid limits must be finite and increasing, got $lims"))
    h = (hi - lo) / n
    return [lo + (i - 1) * h for i in 1:n], h
end

Base.size(grid::PhaseSpaceGrid) = (length(grid.q), length(grid.p))

# Reject mismatched samples before reductions or broadcasting can hide their shape.
function check_size(W::AbstractMatrix, grid::PhaseSpaceGrid)
    size(W) == size(grid) || throw(
        DimensionMismatch(
            "matrix has size $(size(W)), but the grid has size $(size(grid))",
        ),
    )
    return nothing
end

function Base.show(io::IO, grid::PhaseSpaceGrid)
    nq, np = size(grid)
    qmax = first(grid.q) + nq * grid.dq
    pmax = first(grid.p) + np * grid.dp
    print(io, "PhaseSpaceGrid(q ∈ [", first(grid.q), ", ", qmax, ") × ", nq)
    print(io, ", p ∈ [", first(grid.p), ", ", pmax, ") × ", np, ")")
    return nothing
end

"""
    on_grid(f, grid::PhaseSpaceGrid)

Evaluate `f(q, p)` at every point of `grid` and return the `nq × np` matrix
`[f(q, p) for q in grid.q, p in grid.p]`.

# Examples

```julia
grid = PhaseSpaceGrid((-8.0, 8.0), 64, (-8.0, 8.0), 64)
W0 = on_grid((q, p) -> coherent_wigner(q, p; q0 = 2.0, mass = 1.0, omega = 1.0), grid)
```
"""
on_grid(f, grid::PhaseSpaceGrid) = [f(q, p) for q in grid.q, p in grid.p]
