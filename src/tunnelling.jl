"""
    tunnelling_rates(times, product_population; equilibrium_product, tspan = nothing,
                     min_deviation = 1e-8, population_atol = 1e-8)

Extract constant forward and backward inter-well transfer rates from population
relaxation. For normalised two-state kinetics,

    dP/dt = k_forward * (1 - P) - k_backward * P,
    P(t) = equilibrium_product + amplitude * exp(-relaxation_rate * (t - t₀)).

Fit `log(abs(P - equilibrium_product))` by unweighted linear least squares, with
`t₀` the first selected sample. Return a named tuple containing `relaxation_rate`,
`forward_rate = equilibrium_product * relaxation_rate`,
`backward_rate = (1 - equilibrium_product) * relaxation_rate`, `equilibrium_product`,
signed `amplitude`, log-space `r_squared`, population-space `rmse`, and the selected
`times`, `population`, and `fitted_population` vectors. Rates have inverse-time units.
Unevenly spaced samples are supported and each sample has equal fitting weight.

`equilibrium_product` is a fraction in `[0, 1]`, supplied independently (for example
from a converged coupled equilibrium); it is never inferred from the last sample.
The irreversible limits `0` and `1` are allowed. `tspan = (start, stop)` selects an
inclusive fitting window; `nothing` uses all samples. At least three selected
samples are required. Times must be finite, real and strictly increasing.
Populations must be finite real fractions in `[0, 1]` up to `population_atol`;
they are neither clipped nor renormalised.

Select a window after preparation transients but before deviations reach the
numerical noise floor. Every selected deviation must exceed the positive finite
`min_deviation` in magnitude and have the same sign. A nondecaying fit is rejected.
No samples are silently dropped. `r_squared` and `rmse` diagnose the fit, but do
not establish a kinetic regime: also compare rates across fitting windows and
converge the grid, hierarchy and bath expansion. An incorrect equilibrium value
can bias the rates even with a good fit.

These are total inter-well transfer rates. Population relaxation alone cannot
separate tunnelling from thermal activation. Coherent oscillations generally do
not admit constant rates; an energy splitting divided by `hbar` is an angular
frequency, not the rate fitted here. The model assumes an undriven, closed
reactant/product partition with conserved total population.
"""
function tunnelling_rates(
    times::AbstractVector,
    product_population::AbstractVector;
    equilibrium_product::Real,
    tspan = nothing,
    min_deviation::Real = 1e-8,
    population_atol::Real = 1e-8,
)
    length(times) == length(product_population) ||
        throw(DimensionMismatch("times and product_population must have equal lengths"))
    length(times) >= 3 || throw(ArgumentError("at least three samples are required"))
    all(t -> t isa Real && isfinite(t), times) ||
        throw(ArgumentError("sample times must be finite and real"))
    all(p -> p isa Real && isfinite(p), product_population) ||
        throw(ArgumentError("product populations must be finite and real"))
    t, population = Float64.(times), Float64.(product_population)
    all(isfinite, t) && all(isfinite, population) ||
        throw(ArgumentError("samples must be finite in Float64"))
    all(x -> isfinite(x) && x > 0, diff(t)) ||
        throw(ArgumentError("sample times must be strictly increasing"))
    peq, floor, atol = Float64.((equilibrium_product, min_deviation, population_atol))
    isfinite(peq) && 0 <= peq <= 1 ||
        throw(ArgumentError("equilibrium_product must be finite and in [0, 1]"))
    isfinite(floor) && floor > 0 ||
        throw(ArgumentError("min_deviation must be finite and positive"))
    isfinite(atol) && atol >= 0 ||
        throw(ArgumentError("population_atol must be finite and nonnegative"))
    all(p -> -atol <= p <= 1 + atol, population) || throw(
        ArgumentError("product populations must lie in [0, 1] within population_atol"),
    )

    if !isnothing(tspan)
        length(tspan) == 2 && all(x -> x isa Real && isfinite(x), tspan) ||
            throw(ArgumentError("tspan must contain two finite real times"))
        lo, hi = Float64.(tspan)
        isfinite(lo) && isfinite(hi) && lo < hi ||
            throw(ArgumentError("tspan must be finite and strictly increasing"))
        selected = findall(x -> lo <= x <= hi, t)
        t, population = t[selected], population[selected]
    end
    length(t) >= 3 ||
        throw(ArgumentError("the fitting window must contain at least three samples"))
    deviation = population .- peq
    all(x -> abs(x) > floor, deviation) || throw(
        ArgumentError("selected deviations reach min_deviation; choose an earlier window"),
    )
    all(x -> signbit(x) == signbit(first(deviation)), deviation) || throw(
        ArgumentError("selected populations cross equilibrium; choose a kinetic window"),
    )

    # Shift and scale time before regression, avoiding overflow of time variances
    # and of an amplitude extrapolated back to an arbitrary absolute time zero.
    duration = last(t) - first(t)
    isfinite(duration) && duration > 0 ||
        throw(ArgumentError("fitting duration must be finite and positive"))
    x = (t .- first(t)) ./ duration
    y = log.(abs.(deviation))
    xmean = sum(x) / length(x)
    ymean = first(y) + sum(y .- first(y)) / length(y)
    dx, dy = x .- xmean, y .- ymean
    total = sum(abs2, dy)
    total > 0 || throw(ArgumentError("constant populations do not identify a decay rate"))
    slope = dot(dx, dy) / sum(abs2, dx)
    rate = -slope / duration
    isfinite(rate) && rate > 0 ||
        throw(ArgumentError("the fitted population deviation must decay at a finite rate"))
    intercept = ymean - slope * xmean
    amplitude = sign(first(deviation)) * exp(intercept)
    fitted_log = intercept .+ slope .* x
    fitted = peq .+ sign(first(deviation)) .* exp.(fitted_log)
    r_squared = 1 - sum(abs2, y .- fitted_log) / total
    rmse = sqrt(sum(abs2, population .- fitted) / length(t))
    isfinite(amplitude) && all(isfinite, fitted) && isfinite(r_squared) && isfinite(rmse) ||
        throw(ArgumentError("population fit overflowed; choose another fitting window"))
    return (
        relaxation_rate = rate,
        forward_rate = peq * rate,
        backward_rate = (1 - peq) * rate,
        equilibrium_product = peq,
        amplitude = amplitude,
        r_squared = r_squared,
        rmse = rmse,
        times = t,
        population = population,
        fitted_population = fitted,
    )
end

"""
    tunnelling_rates(states, op; times, equilibrium_product, dividing_surface = 0.0,
                     product_side = :right, norm_atol = 1e-6, kwargs...)
    tunnelling_rates(sol::AbstractODESolution; equilibrium_product, kwargs...)

Extract position populations with [`probability`](@ref) and fit [`tunnelling_rates`](@ref).
The product occupies `q ≥ dividing_surface` for `product_side = :right`, or
`q ≤ dividing_surface` for `:left`; the complementary half-box is the reactant.
The finite dividing surface must lie strictly inside the grid's position box.
`equilibrium_product` must refer to this same product region.

Supply Wigner matrices for Wigner–Moyal or Caldeira–Leggett operators, or complete
3D hierarchies for HEOM. Only the physical hierarchy member is measured. Each
physical state must be finite and real, and have unit integral within the finite,
nonnegative `norm_atol`. No renormalisation is performed. All supplied states are
checked, including those outside the fitting window. Driven operators are rejected.

The solution overload obtains states, times and operator from the solution. The
integration must have succeeded and saved its requested endpoint. A final saved
sample is not assumed to be equilibrium. Remaining keywords select the fitting
window and tolerances as in the time/population method.
"""
function tunnelling_rates(
    states::AbstractVector{<:AbstractArray},
    op::AbstractPhaseSpaceOperator;
    times::AbstractVector,
    dividing_surface::Real = 0.0,
    product_side::Symbol = :right,
    norm_atol::Real = 1e-6,
    kwargs...,
)
    is_time_dependent(op) &&
        throw(ArgumentError("constant tunnelling rates require an undriven operator"))
    length(states) == length(times) ||
        throw(DimensionMismatch("times must contain one entry per state"))
    surface, atol = Float64.((dividing_surface, norm_atol))
    isfinite(surface) &&
    first(op.grid.q) < surface < first(op.grid.q) + length(op.grid.q) * op.grid.dq ||
        throw(ArgumentError("dividing_surface must be finite and strictly inside the box"))
    product_side in (:left, :right) ||
        throw(ArgumentError("product_side must be :left or :right"))
    isfinite(atol) && atol >= 0 ||
        throw(ArgumentError("norm_atol must be finite and nonnegative"))
    window = product_side === :right ? (surface, Inf) : (-Inf, surface)
    # The same weights apply to every saved state, including off-grid surfaces.
    weights = interval_weights(op.grid.q, op.grid.dq, window)
    populations = Float64[]
    for state in states
        op isa WignerHEOM && check_hierarchy_size(state, op)
        eltype(state) <: Real || throw(ArgumentError("Wigner samples must be real"))
        W = physical_state(state, op)
        W isa AbstractMatrix || throw(ArgumentError("physical states must be matrices"))
        check_size(W, op.grid)
        all(isfinite, W) || throw(ArgumentError("physical Wigner samples must be finite"))
        norm = phase_space_integral(W, op.grid)
        isfinite(norm) && abs(norm - 1) <= atol || throw(
            ArgumentError(
                "physical state integral must be one within norm_atol; got $norm",
            ),
        )
        push!(populations, dot(weights, position_density(W, op.grid)))
    end
    return tunnelling_rates(times, populations; kwargs...)
end

function tunnelling_rates(sol::AbstractODESolution; kwargs...)
    op = sol.prob.p
    op isa AbstractPhaseSpaceOperator ||
        throw(ArgumentError("solution parameters must be a phase-space operator"))
    successful_retcode(sol) && !isempty(sol.t) && last(sol.t) == last(sol.prob.tspan) ||
        throw(
            ArgumentError(
                "rate extraction requires a successful solution saved to its endpoint",
            ),
        )
    return tunnelling_rates(sol.u, op; times = sol.t, kwargs...)
end
