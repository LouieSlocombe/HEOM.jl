"""
    LinearResponseResult

Impulse response from [`linear_response`](@ref). `times` contains delays from the
start of the requested interval, `response` contains the real causal response, and
`solution` retains the propagated perturbation (including every HEOM auxiliary).
These signed perturbations have zero trace; they are not normalised density states.
"""
struct LinearResponseResult{S}
    times::Vector{Float64}
    response::Vector{Float64}
    solution::S
end

response_hamiltonian(op::AbstractWignerMoyal) = op
response_hamiltonian(op::WignerHEOM) = op.hamiltonian
response_hamiltonian(op::AbstractCaldeiraLeggett) = op.hamiltonian

function response_problem(U, tspan, op::AbstractWignerMoyal, jacobian)
    jacobian === :matrixfree ||
        throw(ArgumentError("jacobian is only configurable for HEOM response problems"))
    return wigner_moyal_problem(U, tspan, op)
end

function response_problem(U, tspan, op::AbstractCaldeiraLeggett, jacobian)
    jacobian === :matrixfree ||
        throw(ArgumentError("jacobian is only configurable for HEOM response problems"))
    return caldeira_leggett_problem(U, tspan, op)
end

response_problem(U, tspan, op::WignerHEOM, jacobian) = heom_problem(U, tspan, op; jacobian)

function response_interval(tspan)
    length(tspan) == 2 && all(t -> t isa Real && isfinite(t), tspan) ||
        throw(ArgumentError("response tspan must have two finite real endpoints"))
    a, b = Float64.(tspan)
    isfinite(a) && isfinite(b) && isfinite(b - a) && b > a ||
        throw(ArgumentError("response tspan must be finite and strictly increasing"))
    return (a, b)
end

"""
    linear_response_problem(Ueq, tspan, op; dipole = identity, jacobian = :matrixfree)
    linear_response_problem(eq::EquilibriumResult, tspan; kwargs...)

Prepare the causal response to `V(q,t) = V₀(q) - E(t)*dipole(q)` under an
**undriven** Wigner–Moyal, Caldeira–Leggett, or HEOM operator. The initial
perturbation is `(i/ħ)[μ,ρeq]`, or `-L_μ W` in Wigner space, with no kinetic term.
For `μ(q)=q` it is `-∂p W`. The dipole uses the operator's spatial discretisation
and Moyal truncation; nonlinear dipoles retain the corresponding quantum terms.

`Ueq` is a matrix or, for HEOM, a complete 3D hierarchy. Every auxiliary is
perturbed, preserving its scaling and system–bath correlations. A HEOM matrix
initialises zero auxiliaries, which generally is not a coupled equilibrium.
The input is copied and never renormalised. The caller must establish stationarity
of array inputs; the `EquilibriumResult` overload requires `eq.converged`.

The returned `ODEProblem` propagates a signed, zero-trace perturbation, not a
physical state. Use [`linear_response`](@ref) to solve and measure it. `tspan`
must be finite and strictly increasing. Equilibrium spectroscopy requires a
static operator; driven operators are rejected. The operator's mutable work
buffers must not be shared by concurrent solves.
"""
function linear_response_problem(
    Ueq::AbstractArray,
    tspan,
    op::AbstractPhaseSpaceOperator;
    dipole = identity,
    jacobian::Symbol = :matrixfree,
)
    is_time_dependent(op) &&
        throw(ArgumentError("linear response requires an undriven operator"))
    interval = response_interval(tspan)
    prob = response_problem(Ueq, interval, op, jacobian)
    h = response_hamiltonian(op)
    discretization = h isa SpectralWignerMoyal ? Spectral() : FiniteDifference(h.order)
    kick = wigner_moyal_operator(
        op.grid;
        mass = op.mass,
        hbar = op.hbar,
        potential = dipole,
        discretization,
        moyal_terms = h.moyal_terms,
    )
    is_time_dependent(kick) &&
        throw(ArgumentError("dipole must be a time-independent function of position"))
    seed = similar(prob.u0)
    if op isa WignerHEOM
        for a in axes(seed, 3)
            potential_action!(@view(seed[:, :, a]), @view(prob.u0[:, :, a]), kick)
        end
    else
        potential_action!(seed, prob.u0, kick)
    end
    prob.u0 .= .-seed
    return prob
end

function linear_response_problem(eq::EquilibriumResult, tspan; kwargs...)
    eq.converged || throw(ArgumentError("linear response requires converged equilibrium"))
    return linear_response_problem(eq.hierarchy, tspan, eq.restart.operator; kwargs...)
end

"""
    linear_response(Ueq, tspan, op, alg; dipole = identity, observable = dipole,
                    jacobian = :matrixfree, kwargs...)
    linear_response(eq::EquilibriumResult, tspan, alg; kwargs...)

Solve [`linear_response_problem`](@ref) and return a [`LinearResponseResult`](@ref).
`observable(q)` is a finite real position observable; by default it is the dipole.
The signal is `R(t) = Tr[A exp(L₀t) (i/ħ)[μ,ρeq]]`, with delays measured from
`tspan[1]`. For a harmonic oscillator and `A=μ=q`,
`R(t)=sin(ω*t)/(mass*ω)`. To first order,
`δ⟨A(t)⟩ = ∫₀ᵗ R(t-s) E(s) ds`.

Extra keywords pass to `solve`; use `saveat` to resolve oscillations and a long
enough interval to resolve spectral lines. The initial sample and full states
must be saved (`save_start=true`, `save_on=true`, no `save_idxs`). Failed or
prematurely terminated integrations throw an error rather than returning a
partial spectrum. For HEOM, only the physical root is measured, but all
perturbed auxiliaries are propagated. See [`absorption_spectrum`](@ref).
"""
function linear_response(
    Ueq::AbstractArray,
    tspan,
    op::AbstractPhaseSpaceOperator,
    alg;
    dipole = identity,
    observable = dipole,
    jacobian::Symbol = :matrixfree,
    kwargs...,
)
    get(kwargs, :save_start, true) ||
        throw(ArgumentError("linear response requires save_start=true"))
    get(kwargs, :save_on, true) ||
        throw(ArgumentError("linear response requires save_on=true"))
    get(kwargs, :save_end, true) ||
        throw(ArgumentError("linear response requires save_end=true"))
    get(kwargs, :save_idxs, nothing) === nothing ||
        throw(ArgumentError("linear response requires full saved states"))
    values = Float64[potential_value(observable, q) for q in op.grid.q]
    all(isfinite, values) || throw(ArgumentError("observable must be finite in Float64"))
    prob = linear_response_problem(Ueq, tspan, op; dipole, jacobian)
    options = merge((; save_start = true, save_end = true), (; kwargs...))
    sol = solve(prob, alg; options...)
    successful_retcode(sol) && !isempty(sol.t) && last(sol.t) == last(prob.tspan) || throw(
        ErrorException(
            "linear-response integration did not reach its endpoint: $(sol.retcode)",
        ),
    )
    first(sol.t) == first(prob.tspan) ||
        throw(ArgumentError("linear response requires the initial saved sample"))
    times = Float64.(sol.t .- first(prob.tspan))
    response = [
        dot(values, position_density(physical_state(U, op), op.grid)) * op.grid.dq for
        U in sol.u
    ]
    all(isfinite, response) ||
        throw(ErrorException("linear response contains nonfinite values"))
    return LinearResponseResult(times, response, sol)
end

function linear_response(eq::EquilibriumResult, tspan, alg; kwargs...)
    eq.converged || throw(ArgumentError("linear response requires converged equilibrium"))
    return linear_response(eq.hierarchy, tspan, eq.restart.operator, alg; kwargs...)
end

"""
    absorption_spectrum(times, response; frequencies, broadening = 0.0)
    absorption_spectrum(result::LinearResponseResult; kwargs...)

Transform a real causal dipole response into the complex susceptibility

    χ(ω) = ∫₀ᵀ exp((iω - broadening)*t) R(t) dt

and relative absorption `intensity = ω*imag(χ)`. Return a named tuple with
`frequencies`, `susceptibility`, and `intensity`. Frequencies are **angular**
frequencies in inverse time, finite and nonnegative, in the caller's order.
No electromagnetic prefactor, refractive index, or arbitrary peak normalisation
is included. Negative numerical lobes are retained.

The trapezoid rule supports strictly increasing, nonuniform finite sample times;
`times[1]` is treated as the impulse time and subtracted. At least two samples
are required, including the impulse and end of the response interval. `response`
must be real and finite. `broadening ≥ 0` is an exponential damping rate (Lorentzian
half-width), not a frequency in cycles per time. It is added to physical damping.

Choose sample spacing to resolve the highest requested frequency, and converge
the time window and quadrature. Truncating a slowly decaying response causes
spectral ringing. A spectrum of equilibrium absorption requires a stationary
initial state and matching perturbing/measured dipoles in `linear_response`.
"""
function absorption_spectrum(
    times::AbstractVector,
    response::AbstractVector;
    frequencies::AbstractVector,
    broadening::Real = 0.0,
)
    length(times) == length(response) ||
        throw(DimensionMismatch("times and response must have equal lengths"))
    length(times) >= 2 || throw(ArgumentError("at least two response samples are required"))
    all(x -> x isa Real && isfinite(x), times) ||
        throw(ArgumentError("sample times must be finite and real"))
    all(x -> x isa Real && isfinite(x), response) ||
        throw(ArgumentError("response samples must be finite and real"))
    all(x -> x isa Real && isfinite(x) && x >= 0, frequencies) ||
        throw(ArgumentError("frequencies must be finite, real and nonnegative"))
    η = Float64(broadening)
    isfinite(η) && η >= 0 ||
        throw(ArgumentError("broadening must be finite and nonnegative"))
    t, R, ω = Float64.(times), Float64.(response), Float64.(frequencies)
    all(isfinite, t) && all(isfinite, R) && all(isfinite, ω) ||
        throw(ArgumentError("spectrum inputs must be finite in Float64"))
    t .-= first(t)
    steps = diff(t)
    all(x -> isfinite(x) && x > 0, steps) ||
        throw(ArgumentError("sample times must be strictly increasing"))
    weighted = R .* exp.(-η .* t)
    χ = zeros(ComplexF64, length(ω))
    for k in eachindex(ω)
        isfinite(ω[k] * last(t)) ||
            throw(ArgumentError("Fourier phase overflowed; rescale times or frequencies"))
        previous = weighted[1]
        for j in 2:length(t)
            current = weighted[j] * cis(ω[k] * t[j])
            χ[k] += (previous + current) * (steps[j-1] / 2)
            previous = current
        end
    end
    intensity = ω .* imag.(χ)
    all(isfinite, χ) && all(isfinite, intensity) ||
        throw(ArgumentError("spectrum overflowed; rescale times, frequencies or response"))
    return (; frequencies = ω, susceptibility = χ, intensity)
end

absorption_spectrum(result::LinearResponseResult; kwargs...) =
    absorption_spectrum(result.times, result.response; kwargs...)
