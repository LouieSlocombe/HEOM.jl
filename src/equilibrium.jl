"""
    EquilibriumResult

Result of [`equilibrate`](@ref). Fields:

  - `hierarchy`: the complete final 3D state, including correlated auxiliaries.
  - `residuals`: vectors `absolute`, `relative`, and `scaled`, in hierarchy order.
    For member `a`, these are `maximum(abs, dU[:, :, a])`, that value divided by
    `maximum(abs, U[:, :, a])`, and that value divided by the stationarity tolerance
    `stationarity_abstol + stationarity_reltol * maximum(abs, U[:, :, a])`.
    The relative residual of a zero member is zero if its RHS is zero, otherwise `Inf`.
  - `converged`: true only for a successful solver exit with every scaled residual ≤ 1
    and finite state with positive root integral.
  - `status`: `:stationary`, `:time_limit`, `:solver_failure`, or `:terminated` (another
    callback stopped integration before convergence).
  - `retcode`: the SciML solver return code (`ReturnCode.Success` if no solve was needed).
  - `restart`: time, requested `tspan`, `operator`, copied `indices`, `depth`, `scaled`,
    `jacobian`, final root `norm`, and `initial_norm`. The operator is retained by reference,
    including its grid, bath and work buffers; do not mutate it or use it concurrently.

Use `physical_wigner(result)` to view the root and `heom_problem(result, tspan)` to
restart with all correlations. Inspect `converged` before treating a result as equilibrium.
"""
struct EquilibriumResult{R,M}
    hierarchy::Array{Float64,3}
    residuals::R
    converged::Bool
    status::Symbol
    retcode::ReturnCode.T
    restart::M
end

function equilibrium_residuals!(dU, U, op, t, abstol, reltol)
    heom!(dU, U, op, t)
    absolute = Vector{Float64}(undef, size(U, 3))
    relative, scaled = similar(absolute), similar(absolute)
    for a in eachindex(absolute)
        amplitude = maximum(abs, @view U[:, :, a])
        residual = maximum(abs, @view dU[:, :, a])
        if !isfinite(amplitude) || !isfinite(residual)
            absolute[a] = relative[a] = scaled[a] = Inf
        else
            absolute[a] = residual
            relative[a] =
                iszero(amplitude) ? (iszero(residual) ? 0.0 : Inf) : residual / amplitude
            scaled[a] = residual / (abstol + reltol * amplitude)
        end
    end
    return (; absolute, relative, scaled)
end

function equilibrium_result(U, residuals, op, t, tspan, jacobian, initial_norm, retcode)
    integral = phase_space_integral(physical_wigner(U), op.grid)
    converged =
        successful_retcode(retcode) &&
        all(x -> x <= 1, residuals.scaled) &&
        all(isfinite, U) &&
        isfinite(integral) &&
        integral > 0
    status = if !successful_retcode(retcode)
        :solver_failure
    elseif converged
        :stationary
    elseif t >= last(tspan)
        :time_limit
    else
        :terminated
    end
    restart = (;
        time = t,
        tspan,
        operator = op,
        indices = hierarchy_indices(op),
        depth = op.depth,
        scaled = op.scaled,
        jacobian,
        norm = integral,
        initial_norm,
    )
    return EquilibriumResult(copy(U), residuals, converged, status, retcode, restart)
end

"""
    equilibrate(U0, tspan, op, alg; stationarity_abstol = 1e-8,
                stationarity_reltol = 1e-6, check_interval = 1.0,
                jacobian = :matrixfree, callback = nothing, kwargs...)

Prepare a correlated equilibrium by real-time relaxation under a [`heom_operator`](@ref).
`U0` is a Wigner matrix (zero initial auxiliaries) or a complete hierarchy in `op`'s
scaling convention. `alg` is a SciML ODE algorithm, supplied by the caller as for
`solve`. The finite, forward `tspan` bounds the preparation time. Equal endpoints
only assess the initial state. The input is copied, must have positive finite root
integral, and is never renormalised.

Evaluate the full HEOM RHS initially, after accepted steps separated by at least
`check_interval` time units, and at the final time. Stop when **every** hierarchy
member satisfies

    maximum(abs, dU[:, :, a]) ≤ stationarity_abstol +
                              stationarity_reltol * maximum(abs, U[:, :, a]).

`stationarity_abstol` must be positive and `stationarity_reltol` nonnegative; both
are finite. These are RHS tolerances (state/time and inverse time), separate from
the ODE solver's `abstol` and `reltol`. Residuals use the supplied auxiliary scaling
convention; changing that convention can change the stopping time. `check_interval`
must be positive and finite. A stationary root alone does not imply convergence.

Return an [`EquilibriumResult`](@ref), including on time exhaustion or solver failure.
No stationary state is guaranteed for an arbitrary potential, bath or truncation.
Stationarity tests the discretised, truncated equations; converge the bath expansion,
hierarchy depth and grid separately. The reduced equilibrium of a coupled harmonic
system generally differs from its isolated Gibbs state.

Extra keywords go to `solve`, including integration tolerances and `maxiters`. Saving
defaults to the final state only. `save_end = false`, `save_on = false`, and any
`save_idxs` other than `nothing` are rejected because a full final hierarchy is required. An optional
`callback` is combined with the stationarity check; it must preserve the autonomous
operator and its hierarchy convention. The operator owns shared work buffers, as in
[`heom_problem`](@ref).

# Examples

```julia
eq = equilibrate(W0, (0.0, 100.0), op, Vern7(); abstol = 1e-10, reltol = 1e-10)
eq.converged
maximum(eq.residuals.scaled)
prob = heom_problem(eq, (0.0, 10.0))  # keeps every auxiliary
```
"""
function equilibrate(
    U0::AbstractArray,
    tspan,
    op::WignerHEOM,
    alg;
    stationarity_abstol::Real = 1e-8,
    stationarity_reltol::Real = 1e-6,
    check_interval::Real = 1.0,
    jacobian::Symbol = :matrixfree,
    callback = nothing,
    kwargs...,
)
    atol, rtol, interval =
        Float64.((stationarity_abstol, stationarity_reltol, check_interval))
    isfinite(atol) && atol > 0 ||
        throw(ArgumentError("stationarity_abstol must be finite and positive"))
    isfinite(rtol) && rtol >= 0 ||
        throw(ArgumentError("stationarity_reltol must be finite and nonnegative"))
    isfinite(interval) && interval > 0 ||
        throw(ArgumentError("check_interval must be finite and positive"))
    length(tspan) == 2 && all(t -> t isa Real && isfinite(t), tspan) ||
        throw(ArgumentError("tspan must contain two finite real times"))
    times = Tuple(Float64.(tspan))
    all(isfinite, times) && times[2] >= times[1] ||
        throw(ArgumentError("tspan must be finite and forward in Float64"))
    get(kwargs, :save_end, true) === false &&
        throw(ArgumentError("equilibrate requires save_end = true"))
    get(kwargs, :save_on, true) === false &&
        throw(ArgumentError("equilibrate requires save_on = true"))
    get(kwargs, :save_idxs, nothing) === nothing ||
        throw(ArgumentError("equilibrate requires the complete hierarchy; omit save_idxs"))

    prob = heom_problem(U0, times, op; jacobian)
    initial_norm = phase_space_integral(physical_wigner(prob.u0), op.grid)
    isfinite(initial_norm) && initial_norm > 0 ||
        throw(ArgumentError("initial physical member must have positive finite integral"))
    dU = similar(prob.u0)
    residuals = equilibrium_residuals!(dU, prob.u0, op, times[1], atol, rtol)
    if all(x -> x <= 1, residuals.scaled) || times[1] == times[2]
        return equilibrium_result(
            prob.u0,
            residuals,
            op,
            times[1],
            times,
            jacobian,
            initial_norm,
            ReturnCode.Success,
        )
    end

    last_check = Ref(times[1])
    condition(u, t, integrator) = t - last_check[] >= interval
    function check_stationarity!(integrator)
        last_check[] = integrator.t
        current = equilibrium_residuals!(dU, integrator.u, op, integrator.t, atol, rtol)
        if all(x -> x <= 1, current.scaled)
            terminate!(integrator)
        end
        return nothing
    end
    stationarity =
        DiscreteCallback(condition, check_stationarity!; save_positions = (false, false))
    options = merge(
        (save_everystep = false, save_start = false, dense = false),
        (; kwargs...),
        (save_end = true, callback = CallbackSet(callback, stationarity)),
    )
    sol = solve(prob, alg; options...)
    U, t = last(sol.u), last(sol.t)
    residuals = equilibrium_residuals!(dU, U, op, t, atol, rtol)
    return equilibrium_result(
        U,
        residuals,
        op,
        t,
        times,
        jacobian,
        initial_norm,
        sol.retcode,
    )
end

"""
    physical_wigner(result::EquilibriumResult)

View of the physical member of an equilibrium preparation result.
"""
physical_wigner(result::EquilibriumResult) = physical_wigner(result.hierarchy)

"""
    heom_problem(result::EquilibriumResult, tspan; jacobian = result.restart.jacobian)

Restart a preparation result over `tspan` with its full hierarchy and retained
operator. The hierarchy is copied; the operator's work buffers are shared. This
also permits continuing a preparation that has not yet converged.
"""
function heom_problem(
    result::EquilibriumResult,
    tspan;
    jacobian::Symbol = result.restart.jacobian,
)
    return heom_problem(result.hierarchy, tspan, result.restart.operator; jacobian)
end
