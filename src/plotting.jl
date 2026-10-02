# Recipes own their wrapper types, leaving generic matrix and ODE-solution plots alone.
# Plots.jl supplies the rendering methods when the user loads it.

function plotting_state(W::AbstractMatrix, grid::PhaseSpaceGrid)
    check_size(W, grid)
    eltype(W) <: Real || throw(ArgumentError("Wigner samples must be real"))
    all(isfinite, W) || throw(ArgumentError("Wigner samples must be finite"))
    return W, grid
end

function plotting_state(sol::AbstractODESolution, index::Integer = lastindex(sol.u))
    op = sol.prob.p
    op isa AbstractPhaseSpaceOperator ||
        throw(ArgumentError("solution parameters must be a phase-space operator"))
    return plotting_state(sol.u[index], op.grid)
end

"""
    wignerplot(W, grid; kwargs...)
    wignerplot(sol[, index]; kwargs...)
    wignerplot!([plot,] args...; kwargs...)

Plot a signed Wigner heatmap with position `q` horizontal and momentum `p` vertical.
Load `Plots` before calling this function. `W[i, j]` is plotted at `(grid.q[i], grid.p[j])`.
For a phase-space ODE solution, `index` selects a saved state and defaults to the last.

The default diverging colour scale is centred at zero, with limits `±maximum(abs, W)`
(or `(-1, 1)` for a zero matrix). Samples must be real and finite. Values are neither
clipped nor renormalised. All standard Plots attributes can be supplied, for example
`clims = (-0.3, 0.3)` to compare states using a shared scale. The `!` form adds to an
existing plot.
"""
@userplot WignerPlot

@recipe function f(wp::WignerPlot)
    W, grid = plotting_state(wp.args...)
    amplitude = maximum(abs, W)
    limit = iszero(amplitude) ? 1.0 : amplitude
    seriestype --> :heatmap
    xguide --> "q"
    yguide --> "p"
    colorbar_title --> "W(q, p)"
    colorbar --> true
    seriescolor --> :RdBu
    clims --> (-limit, limit)
    aspect_ratio --> 1
    legend --> false
    return grid.q, grid.p, permutedims(W)
end

"""
    marginalplot(W, grid; kwargs...)
    marginalplot(sol[, index]; kwargs...)
    marginalplot!([plot,] args...; kwargs...)

Plot position and momentum densities in separate panels, using [`position_density`](@ref)
and [`momentum_density`](@ref). Load `Plots` first. The solution form selects a saved
state, defaulting to the last. The densities retain the norm of `W`, including any
numerical drift, and are not clipped. Standard Plots attributes override the defaults.
The `!` form adds curves to the corresponding panels of an existing marginal plot.
"""
@userplot MarginalPlot

@recipe function f(mp::MarginalPlot)
    W, grid = plotting_state(mp.args...)
    layout --> (1, 2)
    size --> (800, 350)
    margin --> (5, :mm)
    legend --> false
    seriestype --> :path
    for (panel, x, density, coordinate) in (
        (1, grid.q, position_density(W, grid), "q"),
        (2, grid.p, momentum_density(W, grid), "p"),
    )
        @series begin
            subplot := panel
            xguide --> coordinate
            yguide --> "P($coordinate)"
            label --> "P($coordinate)"
            x, density
        end
    end
    return nothing
end

function plotting_diagnostics(d::NamedTuple; potential = nothing)
    haskey(d, :t) || throw(
        ArgumentError("diagnostics must include t; use diagnosticsplot(times, values)"),
    )
    return d
end

function plotting_diagnostics(times::AbstractVector, d::NamedTuple; potential = nothing)
    return merge(d, (t = times,))
end

function plotting_diagnostics(sol::AbstractODESolution; potential = nothing)
    isnothing(potential) &&
        throw(ArgumentError("supply potential to plot diagnostics of a solution"))
    return diagnostics(sol; potential)
end

# Validate only selected observables: diagnostic NaNs are intentionally left as gaps.
function diagnostic_fields(d, fields)
    d.t isa AbstractVector && !isempty(d.t) ||
        throw(ArgumentError("diagnostic times must be a nonempty vector"))
    selected = fields isa Symbol ? (fields,) : Tuple(fields)
    isempty(selected) && throw(ArgumentError("select at least one diagnostic field"))
    all(field -> field isa Symbol && field != :t && haskey(d, field), selected) ||
        throw(ArgumentError("fields must name diagnostic columns other than t"))
    allunique(selected) || throw(ArgumentError("diagnostic fields must be unique"))
    for field in selected
        values = getproperty(d, field)
        values isa AbstractVector ||
            throw(ArgumentError("diagnostic $field must be a vector"))
        length(values) == length(d.t) ||
            throw(DimensionMismatch("diagnostic $field and times must have equal length"))
    end
    return selected
end

"""
    diagnosticsplot(d; fields = (:norm, :energy, :purity, :negativity), kwargs...)
    diagnosticsplot(times, d; fields = ..., kwargs...)
    diagnosticsplot(sol; potential, fields = ..., kwargs...)
    diagnosticsplot!([plot,] args...; kwargs...)

Plot selected [`diagnostics`](@ref) columns against saved times, one panel per field.
Load `Plots` first. Pass the result of `diagnostics(sol; potential)`, explicit times and
`diagnostics(states, grid; ...)`, or a phase-space solution with its potential.

`fields` may be a symbol or a tuple/vector of distinct column names, such as
`(:mean_q, :mean_p)` or `(:boundary_q, :boundary_p, :tail_q, :tail_p)`. Fields are shown
in the requested order and labelled by name. Values retain their original units and
normalisation; zero values are shown on linear axes and `NaN`s are left as gaps.
Standard Plots attributes override defaults. The `!` form adds to existing panels;
use the same field selection and order as the original plot.
"""
@userplot DiagnosticsPlot

@recipe function f(
    dp::DiagnosticsPlot;
    fields = (:norm, :energy, :purity, :negativity),
    potential = nothing,
)
    d = plotting_diagnostics(dp.args...; potential)
    selected = diagnostic_fields(d, fields)
    n = length(selected)
    layout --> n
    size --> (800, 250 * cld(n, 2))
    margin --> (5, :mm)
    legend --> false
    seriestype --> :path
    xguide --> "t"
    for (panel, field) in enumerate(selected)
        @series begin
            subplot := panel
            yguide --> string(field)
            label --> string(field)
            d.t, getproperty(d, field)
        end
    end
    return nothing
end
