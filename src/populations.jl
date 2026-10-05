# Integrals of the periodic trigonometric cardinal functions over a clipped interval.
# The even-n Nyquist cosine has half the multiplicity of a paired Fourier mode.
function interval_weights(x, h, limits::Tuple{Real,Real})
    limits[1] <= limits[2] ||
        throw(ArgumentError("window limits must be nondecreasing, got $limits"))
    n = length(x)
    lo, hi = first(x), first(x) + n * h
    a, b = clamp(limits[1], lo, hi), clamp(limits[2], lo, hi)
    a == lo && b == hi && return fill(h, n)
    weights = fill((b - a) / n, n)
    for k in 1:(n÷2)
        κ = 2π * k / (n * h)
        c = iseven(n) && k == n ÷ 2 ? 1 : 2
        for i in eachindex(x)
            weights[i] += c * (sin(κ * (b - x[i])) - sin(κ * (a - x[i]))) / (n * κ)
        end
    end
    return weights
end

"""
    probability(W, grid::PhaseSpaceGrid; q = (-Inf, Inf), p = (-Inf, Inf))

Window integral `∫q[1]^q[2] ∫p[1]^p[2] W(q′, p′) dp′ dq′`, evaluated exactly for the
periodic trigonometric interpolant of `W`. Limits are clipped to the grid's box and need
not be grid points.

A position-only or momentum-only window gives a genuine probability for a physical state,
including reactant or product populations across a dividing surface. A joint position and
momentum window gives only a quasi-probability and can be negative. Nothing is renormalised.
"""
function probability(
    W::AbstractMatrix,
    grid::PhaseSpaceGrid;
    q = (-Inf, Inf),
    p = (-Inf, Inf),
)
    check_size(W, grid)
    wq = interval_weights(grid.q, grid.dq, q)
    wp = interval_weights(grid.p, grid.dp, p)
    return dot(wq, W, wp)
end

"""
    probability_current(W, grid::PhaseSpaceGrid; mass)

Position probability current `j(qᵢ) = ∫ (p / mass) W(qᵢ, p) dp`, as a vector on `grid.q`.
The momentum integral uses the periodic rectangle rule.
"""
function probability_current(W::AbstractMatrix, grid::PhaseSpaceGrid; mass::Real)
    check_size(W, grid)
    return (W * grid.p) .* (grid.dp / mass)
end

# Use the same right-hand side as the integrator, including the selected discretisation.
function time_derivative(W::AbstractMatrix, op::AbstractPhaseSpaceOperator)
    check_size(W, op.grid)
    dW = similar(W, Float64)
    phase_space_rhs!(dW, W, op, 0.0)
    return dW
end

function time_derivative(U::AbstractArray{<:Real,3}, op::WignerHEOM)
    dU = similar(U, Float64)
    heom!(dU, U, op, 0.0)
    return physical_wigner(dU)
end

function phase_space_rhs!(dW, W, op::WignerHEOM, t)
    throw(ArgumentError("HEOM rates require the full 3D hierarchy, including auxiliaries"))
end

"""
    expectation_rate(f, W, op)

Instantaneous rate `d⟨f⟩/dt = ∫∫ f(q, p) ∂W/∂t dq dp` for a time-independent Weyl symbol
`f`, using the grid and right-hand side of `op`. This is exact for the semi-discrete
equations, since [`expectation`](@ref) is linear in `W`.
For a HEOM operator pass the full 3D hierarchy in place of `W`; the rate is
evaluated on its physical member and depends on the auxiliary members.
"""
function expectation_rate(f, W::AbstractMatrix, op::AbstractPhaseSpaceOperator)
    return expectation(f, time_derivative(W, op), op.grid)
end

function expectation_rate(f, U::AbstractArray{<:Real,3}, op::WignerHEOM)
    return expectation(f, time_derivative(U, op), op.grid)
end

"""
    probability_rate(W, op; q = (-Inf, Inf), p = (-Inf, Inf))

Instantaneous window population rate `dP/dt`, found by applying [`probability`](@ref) to
the right-hand side of `op`. This is exact for the semi-discrete equations because the
window integral is linear in `W`.
For HEOM, pass the full 3D hierarchy instead of only the physical matrix.

`probability_rate(W, op; q = (q‡, Inf))` is the reactive population flux into the product
region beyond `q‡`. For a resolved state that decays at the box edges, the spectral result
approximates [`probability_current`](@ref) at `q‡` to spectral accuracy. Finite differences
give the flux of their own semi-discrete evolution.
"""
function probability_rate(
    W::AbstractMatrix,
    op::AbstractPhaseSpaceOperator;
    q = (-Inf, Inf),
    p = (-Inf, Inf),
)
    return probability(time_derivative(W, op), op.grid; q, p)
end

function probability_rate(
    U::AbstractArray{<:Real,3},
    op::WignerHEOM;
    q = (-Inf, Inf),
    p = (-Inf, Inf),
)
    return probability(time_derivative(U, op), op.grid; q, p)
end
