module HEOMPlotsExt

using HEOM:
    HEOM, PhaseSpaceGrid, marginalplot, momentum_density, position_density, wignerplot
using Plots: Plots
using SciMLBase: AbstractODESolution

# Validate the entire selection before creating the temporary frame directory.
function animation_data(states, grid; times = nothing, indices = eachindex(states))
    isempty(states) && throw(ArgumentError("animation states must be nonempty"))
    selected = collect(indices)
    isempty(selected) && throw(ArgumentError("select at least one animation frame"))
    all(i -> i isa Integer && !(i isa Bool), selected) ||
        throw(ArgumentError("animation indices must be integers"))
    for i in selected
        checkbounds(states, i)
        HEOM.plotting_state(states[i], grid)
    end
    if !isnothing(times)
        times isa AbstractVector || throw(ArgumentError("times must be a vector"))
        axes(times) == axes(states) ||
            throw(DimensionMismatch("times must have one entry per state"))
        all(t -> t isa Real && isfinite(t), times) ||
            throw(ArgumentError("animation times must be real and finite"))
    end
    return selected
end

function wigner_limits(states, selected)
    amplitude = maximum(i -> maximum(abs, states[i]), selected)
    limit = iszero(amplitude) ? 1.0 : amplitude
    return (-limit, limit)
end

function density_limits(density, states, grid, selected)
    lo, hi = 0.0, 0.0
    for i in selected
        lower, upper = extrema(density(states[i], grid))
        lo, hi = min(lo, lower), max(hi, upper)
    end
    iszero(lo) && iszero(hi) && return (-1.0, 1.0)
    padding = 0.05 * (hi - lo)
    return (lo - padding, hi + padding)
end

function marginal_limits(states, grid, selected)
    qlimits = density_limits(position_density, states, grid, selected)
    plimits = density_limits(momentum_density, states, grid, selected)
    return [qlimits plimits]
end

function record_animation(plotter, states, grid, selected, times; kwargs...)
    animation = Plots.Animation()
    for i in selected
        title = isnothing(times) ? "Frame $i" : "t = $(times[i])"
        figure = plotter(states[i], grid; title, kwargs...)
        Plots.frame(animation, figure)
    end
    return animation
end

function HEOM.wigneranimation(
    states::AbstractVector,
    grid::PhaseSpaceGrid;
    times = nothing,
    indices = eachindex(states),
    kwargs...,
)
    selected = animation_data(states, grid; times, indices)
    return record_animation(
        wignerplot,
        states,
        grid,
        selected,
        times;
        clims = wigner_limits(states, selected),
        xlims = (first(grid.q) - grid.dq / 2, last(grid.q) + grid.dq / 2),
        ylims = (first(grid.p) - grid.dp / 2, last(grid.p) + grid.dp / 2),
        kwargs...,
    )
end

function HEOM.marginalanimation(
    states::AbstractVector,
    grid::PhaseSpaceGrid;
    times = nothing,
    indices = eachindex(states),
    kwargs...,
)
    selected = animation_data(states, grid; times, indices)
    return record_animation(
        marginalplot,
        states,
        grid,
        selected,
        times;
        ylims = marginal_limits(states, grid, selected),
        kwargs...,
    )
end

function animation_operator(sol)
    op = sol.prob.p
    op isa HEOM.AbstractPhaseSpaceOperator ||
        throw(ArgumentError("solution parameters must be a phase-space operator"))
    return op
end

function HEOM.wigneranimation(sol::AbstractODESolution; kwargs...)
    op = animation_operator(sol)
    return HEOM.wigneranimation(sol.u, op.grid; times = sol.t, kwargs...)
end

function HEOM.marginalanimation(sol::AbstractODESolution; kwargs...)
    op = animation_operator(sol)
    return HEOM.marginalanimation(sol.u, op.grid; times = sol.t, kwargs...)
end

end
