"""
    boundary_weight(W, grid::PhaseSpaceGrid; width = 0.05)

Fraction of `sum(abs, W)` in the outer `k = max(1, round(Int, width*n))` points at
both ends of each axis, returned as `(q = ..., p = ...)`. The width must satisfy
`0 < width < 0.5`; boundary bands never overlap. Small values indicate that the Wigner
function has decayed at the edges of the periodic box. The ratios are undefined (`NaN`)
for an identically zero matrix.
"""
function boundary_weight(W::AbstractMatrix, grid::PhaseSpaceGrid; width::Real = 0.05)
    check_size(W, grid)
    0 < width < 0.5 || throw(ArgumentError("width must lie between 0 and 0.5, got $width"))
    nq, np = size(grid)
    # The cap also prevents overlap if floating-point multiplication rounds width*n up.
    kq = min(max(1, round(Int, width * nq)), nq ÷ 2)
    kp = min(max(1, round(Int, width * np)), np ÷ 2)
    total = sum(abs, W)
    q = (sum(abs, @view W[1:kq, :]) + sum(abs, @view W[(nq-kq+1):nq, :])) / total
    p = (sum(abs, @view W[:, 1:kp]) + sum(abs, @view W[:, (np-kp+1):np])) / total
    return (; q, p)
end

"""
    spectral_tail(W, grid::PhaseSpaceGrid; fraction = 1/3)

Relative Fourier `L²` norm in the highest-frequency fraction of each axis, returned as
`(q = ..., p = ...)`. For an axis of length `n`, the tail contains modes
`k > round(Int, (1 - fraction)*n/2)`. The ratio is `sqrt(sum_tail(abs2, F)/sum(abs2, F))`,
with the multiplicities of the real Fourier transform included. Small values indicate
that the grid resolves the Wigner function. Requires `0 < fraction < 1`; the ratios are
undefined (`NaN`) for an identically zero matrix.
"""
function spectral_tail(W::AbstractMatrix, grid::PhaseSpaceGrid; fraction::Real = 1 / 3)
    check_size(W, grid)
    0 < fraction < 1 ||
        throw(ArgumentError("fraction must lie between 0 and 1, got $fraction"))
    return (q = axis_spectral_tail(W, 1, fraction), p = axis_spectral_tail(W, 2, fraction))
end

function axis_spectral_tail(W, dim, fraction)
    n = size(W, dim)
    power = vec(sum(abs2, rfft(W, dim); dims = 3 - dim))
    cutoff = round(Int, (1 - fraction) * n / 2)
    total, tail = 0.0, 0.0
    for i in eachindex(power)
        k = i - 1
        weight = k == 0 || 2k == n ? 1 : 2
        value = weight * power[i]
        total += value
        if k > cutoff
            tail += value
        end
    end
    return sqrt(tail / total)
end

"""
    diagnostics(W, grid::PhaseSpaceGrid; mass, potential, hbar = 1.0)

State observables and grid-health indicators as a named tuple of `Float64` values:
`norm`, `mean_q`, `mean_p`, `var_q`, `var_p`, `cov_qp`, `uncertainty`,
`robertson_schrodinger`, `energy`, `purity`, `negativity`, `boundary_q`, `boundary_p`,
`tail_q`, `tail_p`.

`uncertainty = sqrt(var_q*var_p)` is σqσp and
`robertson_schrodinger = sqrt(var_q*var_p - cov_qp^2)` is √det Σ. Both are at least ħ/2
for a normalised physical state; √det Σ is invariant under harmonic evolution. No
observable is renormalised, so norm drift remains visible. Boundary and spectral
indicators use the defaults of [`boundary_weight`](@ref) and [`spectral_tail`](@ref).
"""
function diagnostics(
    W::AbstractMatrix,
    grid::PhaseSpaceGrid;
    mass::Real,
    potential,
    hbar::Real = 1.0,
)
    check_size(W, grid)
    mean = phase_space_mean(W, grid)
    covariance = phase_space_covariance(W, grid)
    var_q, var_p, cov_qp = covariance[1, 1], covariance[2, 2], covariance[1, 2]
    boundary = boundary_weight(W, grid)
    tail = spectral_tail(W, grid)
    return map(
        Float64,
        (
            norm = phase_space_integral(W, grid),
            mean_q = mean[1],
            mean_p = mean[2],
            var_q = var_q,
            var_p = var_p,
            cov_qp = cov_qp,
            uncertainty = sqrt(var_q * var_p),
            robertson_schrodinger = sqrt(var_q * var_p - cov_qp^2),
            energy = energy(W, grid; mass, potential),
            purity = purity(W, grid; hbar),
            negativity = wigner_negativity(W, grid),
            boundary_q = boundary.q,
            boundary_p = boundary.p,
            tail_q = tail.q,
            tail_p = tail.p,
        ),
    )
end

"""
    diagnostics(states::AbstractVector{<:AbstractMatrix}, grid; mass, potential, hbar = 1.0)

Evaluate [`diagnostics`](@ref) for each state and return a named tuple of vectors. The
additional `autocorrelation` column is the overlap of each state with the first state,
`Tr ρ(0)ρ(t)`, which is the survival probability for a pure initial state. The vector
of states must be nonempty.
"""
function diagnostics(
    states::AbstractVector{<:AbstractMatrix},
    grid::PhaseSpaceGrid;
    mass::Real,
    potential,
    hbar::Real = 1.0,
)
    isempty(states) && throw(ArgumentError("diagnostics needs at least one state"))
    rows = [diagnostics(W, grid; mass, potential, hbar) for W in states]
    autocorrelation = [overlap(first(states), W, grid; hbar) for W in states]
    return merge(columns(rows), (; autocorrelation))
end

function columns(rows::AbstractVector{<:NamedTuple{names}}) where {names}
    return NamedTuple{names}(map(n -> [getfield(r, n) for r in rows], names))
end

"""
    diagnostics(sol::AbstractODESolution; potential)

Evaluate [`diagnostics`](@ref) on the saved states of a Wigner–Moyal solution. The
result puts the saved times `t = sol.t` first, followed by the observable vectors and
`autocorrelation`. Grid, mass and ħ are read from the operator in `sol.prob.p`; the
potential must be supplied explicitly.
"""
function diagnostics(sol::AbstractODESolution; potential)
    op = sol.prob.p
    op isa AbstractWignerMoyal || throw(
        ArgumentError(
            "solution parameters must be a Wigner–Moyal operator; use " *
            "diagnostics(sol.u, grid; mass, potential, hbar) for other solutions",
        ),
    )
    values = diagnostics(sol.u, op.grid; mass = op.mass, potential, hbar = op.hbar)
    return merge((t = sol.t,), values)
end
