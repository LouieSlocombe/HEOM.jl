"""
    ExponentialBath(coefficients, rates; hbar = 1.0, diffusion = 0.0, counterterm = 0.0)

Gaussian bath coupled linearly to position, `H_SB = q B`, with force correlation
`C(t) = ⟨B(t)B(0)⟩ = sum(coefficients[k] * exp(-rates[k]*t))` for `t ≥ 0`.
Rates must be real, finite and positive; coefficients may be complex and must be finite.
The two vectors must have equal length; empty vectors describe an isolated system.
Physical consistency of a user-supplied correlation is the caller's responsibility.

`diffusion ≥ 0` adds a Markovian remainder `diffusion * ∂p²W` on every hierarchy member.
`counterterm ≥ 0` adds `counterterm*q²` to the supplied system potential.
`hbar` must be finite and positive. Inputs are copied. See [`drude_lorentz_bath`](@ref)
for a thermal bath with its counterterm and Matsubara remainder included.
"""
struct ExponentialBath
    coefficients::Vector{ComplexF64}
    rates::Vector{Float64}
    hbar::Float64
    diffusion::Float64
    counterterm::Float64
    function ExponentialBath(
        coefficients::AbstractVector,
        rates::AbstractVector;
        hbar::Real = 1.0,
        diffusion::Real = 0.0,
        counterterm::Real = 0.0,
    )
        length(coefficients) == length(rates) ||
            throw(DimensionMismatch("bath coefficients and rates must have equal length"))
        all(x -> x isa Real, rates) || throw(ArgumentError("bath rates must be real"))
        c, ν = ComplexF64.(coefficients), Float64.(rates)
        ħ, D, λ = Float64(hbar), Float64(diffusion), Float64(counterterm)
        all(isfinite, c) || throw(ArgumentError("bath coefficients must be finite"))
        all(x -> isfinite(x) && x > 0, ν) ||
            throw(ArgumentError("bath rates must be finite and positive"))
        isfinite(ħ) && ħ > 0 || throw(ArgumentError("hbar must be finite and positive"))
        isfinite(D) && D >= 0 ||
            throw(ArgumentError("diffusion must be finite and nonnegative"))
        isfinite(λ) && λ >= 0 ||
            throw(ArgumentError("counterterm must be finite and nonnegative"))
        return new(c, ν, ħ, D, λ)
    end
end

"""
    drude_lorentz_bath(; reorganization, cutoff, kT, matsubara = 0, hbar = 1.0,
                       terminator = true)

Thermal [`ExponentialBath`](@ref) with Drude spectral density
`J(ω) = 2λγω/(ω² + γ²)`, where `λ = reorganization` and `γ = cutoff`.
The convention is `C(t) = (ħ/π)∫₀∞ J(ω)[coth(ħω/(2kT))*cos(ωt) - i*sin(ωt)]dω`.
For dimensional position coupling, `λ` has units of energy / position².

The retained correlation coefficients are

    ν₀ = γ,                 c₀ = λħγ [cot(ħγ/(2kT)) - i],
    νₖ = 2πk*kT/ħ,          cₖ = 4λγ*kT*νₖ/(νₖ² - γ²),  k = 1,…,matsubara.

The system counterterm is `λq²`. With `terminator = true`, the omitted fast Matsubara
terms are approximated by momentum diffusion
`D = 2λ*kT/γ - sum(real(cₖ)/νₖ)`; this is a bath-correlation tail approximation,
not a closure of the hierarchy depth. Include enough Matsubara terms that the first
omitted rate exceeds `γ`. Increase both `matsubara` and hierarchy `depth` to check
convergence, especially at low temperature. `terminator = false` simply drops the tail.

`λ` must be finite and nonnegative; `γ`, `kT` and `ħ` finite and positive. Coincident
Drude and Matsubara poles are unsupported (their decomposition requires polynomial
exponentials). Zero coupling returns an empty bath. The zero-temperature limit requires
another correlation decomposition.
"""
function drude_lorentz_bath(;
    reorganization::Real,
    cutoff::Real,
    kT::Real,
    matsubara::Integer = 0,
    hbar::Real = 1.0,
    terminator::Bool = true,
)
    λ, γ, θ, ħ = Float64(reorganization), Float64(cutoff), Float64(kT), Float64(hbar)
    isfinite(λ) && λ >= 0 ||
        throw(ArgumentError("reorganization must be finite and nonnegative"))
    isfinite(γ) && γ > 0 || throw(ArgumentError("cutoff must be finite and positive"))
    isfinite(θ) && θ > 0 || throw(ArgumentError("kT must be finite and positive"))
    isfinite(ħ) && ħ > 0 || throw(ArgumentError("hbar must be finite and positive"))
    matsubara >= 0 || throw(ArgumentError("matsubara must be nonnegative"))
    iszero(λ) && return ExponentialBath(ComplexF64[], Float64[]; hbar = ħ)
    x = ħ * γ / (2θ)
    isfinite(x) && x > 0 ||
        throw(ArgumentError("ħ*cutoff/(2kT) must be finite and positive"))
    # Every positive integer x/π is a coincident pole, including poles outside the
    # retained sum. Separate exponentials are ill-conditioned near such a collision.
    pole = round(x / π)
    pole > 0 &&
        abs(x / π - pole) <= sqrt(eps(Float64)) * max(1, pole) &&
        throw(ArgumentError("Drude and Matsubara poles coincide; change cutoff or kT"))
    ν = [γ; [2π * k * θ / ħ for k in 1:matsubara]]
    # Factor out ν² to avoid squaring large dimensional frequencies.
    c = [
        ComplexF64(λ * ħ * γ * cot(x), -λ * ħ * γ);
        [ComplexF64(4λ * θ * (γ / v) / ((1 - γ / v) * (1 + γ / v))) for v in ν[2:end]]
    ]
    D = 0.0
    if terminator
        2π * (matsubara + 1) * θ / ħ > γ || throw(
            ArgumentError("increase matsubara so the first omitted rate exceeds cutoff"),
        )
        # Avoid cancellation between the classical area and cot(x) at high T.
        remainder = abs(x) < 1e-3 ? x / 3 + x^3 / 45 + 2x^5 / 945 : inv(x) - cot(x)
        D = λ * ħ * remainder - sum(real(c[k]) / ν[k] for k in 2:length(ν); init = 0.0)
        tolerance = 64eps(Float64) * max(abs(λ * ħ * remainder), floatmin(Float64))
        D >= -tolerance ||
            throw(ArgumentError("negative Matsubara tail; increase matsubara"))
        D = max(D, 0.0)
    end
    return ExponentialBath(c, ν; hbar = ħ, diffusion = D, counterterm = λ)
end

struct HEOMMomentum{D1,D2}
    first::D1
    second::D2
end

struct WignerHEOM{H<:AbstractWignerMoyal,D} <: AbstractPhaseSpaceOperator
    grid::PhaseSpaceGrid
    mass::Float64
    hbar::Float64
    bath::ExponentialBath
    depth::Int
    indices::Vector{Vector{Int}}
    upper::Matrix{Int}
    lower::Matrix{Int}
    damping::Vector{Float64}
    coordinate_coupling::Vector{Float64}
    hamiltonian::H
    momentum::D
    scratch::Matrix{Float64}
end

# Breadth-first traversal gives graded indices without enumerating a rectangular box.
function build_hierarchy(modes, depth)
    indices = [zeros(Int, modes)]
    lookup = Dict(Tuple(first(indices)) => 1)
    a = 1
    while a <= length(indices)
        n = indices[a]
        if sum(n) < depth
            for k in 1:modes
                next = copy(n)
                next[k] += 1
                key = Tuple(next)
                if !haskey(lookup, key)
                    push!(indices, next)
                    lookup[key] = length(indices)
                end
            end
        end
        a += 1
    end
    upper, lower = zeros(Int, modes, length(indices)), zeros(Int, modes, length(indices))
    for (a, n) in enumerate(indices), k in 1:modes
        next = copy(n)
        next[k] += 1
        upper[k, a] = get(lookup, Tuple(next), 0)
        next[k] -= 2
        lower[k, a] = get(lookup, Tuple(next), 0)
    end
    return indices, upper, lower
end

function heom_momentum(op::SpectralWignerMoyal)
    np = length(op.grid.p)
    κ = 2π .* collect(rfftfreq(np, 1 / op.grid.dp))
    return HEOMMomentum(
        reshape(im .* wavenumbers(np, op.grid.dp) ./ np, 1, :),
        reshape(-κ .^ 2 ./ np, 1, :),
    )
end

function heom_momentum(op::FiniteDifferenceWignerMoyal)
    np = length(op.grid.p)
    return HEOMMomentum(
        periodic_difference_matrix(1, op.order, np, op.grid.dp),
        periodic_difference_matrix(2, op.order, np, op.grid.dp),
    )
end

"""
    heom_operator(grid; mass, potential, bath, depth, discretization = Spectral(),
                  moyal_terms = nothing)

Wigner–Moyal hierarchical equations of motion for linear coordinate coupling to an
[`ExponentialBath`](@ref). The unscaled, real auxiliary Wigner functions obey

    dWₙ/dt = [L_WM(V + counterterm*q²) - Σₖ nₖνₖ + D ∂p²] Wₙ
             + Σₖ ∂p Wₙ₊ₑₖ
             + Σₖ nₖ [real(cₖ) ∂p + 2imag(cₖ)q/ħ] Wₙ₋ₑₖ.

This is the Wigner transform of the coordinate-coupled density-operator HEOM:
`[q,ρ]_W = iħ∂pW` and `{q,ρ}_W = 2qW`. The physical state is `W₀`.
`potential` is the physical system potential; the bath counterterm is added internally.
`hbar` comes from `bath`. Both Hamiltonian discretisations and their Moyal options are
supported (finite differences require integer `moyal_terms`).

All nonnegative multi-indices with `sum(n) ≤ depth` are retained, with the root first.
There are `binomial(length(bath.rates) + depth, depth)` members. Upward terms beyond
`depth` are set to zero (hard cutoff); `depth = 0` omits all explicit bath memory.
Converge the depth, bath decomposition, box and resolution independently. No rescaling
of auxiliaries or positivity correction is applied. The operator owns mutable work
buffers; give each concurrent solve its own operator. See [`heom_problem`](@ref).
"""
function heom_operator(
    grid::PhaseSpaceGrid;
    mass::Real,
    potential,
    bath::ExponentialBath,
    depth::Integer,
    discretization::Union{Spectral,FiniteDifference} = Spectral(),
    moyal_terms::Union{Nothing,Integer} = nothing,
)
    m = Float64(mass)
    isfinite(m) && m > 0 || throw(ArgumentError("mass must be finite and positive"))
    depth >= 0 || throw(ArgumentError("depth must be nonnegative"))
    # Copy and revalidate the bath, since its input vectors remain user-accessible.
    b = ExponentialBath(
        bath.coefficients,
        bath.rates;
        hbar = bath.hbar,
        diffusion = bath.diffusion,
        counterterm = bath.counterterm,
    )
    V(q) = potential(q) + b.counterterm * q^2
    h = wigner_moyal_operator(
        grid;
        mass = m,
        potential = V,
        hbar = b.hbar,
        discretization,
        moyal_terms,
    )
    indices, upper, lower = build_hierarchy(length(b.rates), Int(depth))
    damping = [sum(n .* b.rates) for n in indices]
    all(isfinite, damping) || throw(ArgumentError("hierarchy damping rates must be finite"))
    coordinate_coupling = 2 .* (imag.(b.coefficients) ./ b.hbar)
    all(isfinite, coordinate_coupling) &&
    all(isfinite, depth .* coordinate_coupling) &&
    all(isfinite, depth .* real.(b.coefficients)) ||
        throw(ArgumentError("hierarchy coupling factors must be finite"))
    return WignerHEOM(
        grid,
        m,
        b.hbar,
        b,
        Int(depth),
        indices,
        upper,
        lower,
        damping,
        coordinate_coupling,
        h,
        heom_momentum(h),
        zeros(size(grid)),
    )
end

"""
    hierarchy_indices(op)

Copy of the multi-indices in the order used by the third axis of a HEOM state.
The first is the all-zero physical member; subsequent entries are ordered by tier.
"""
hierarchy_indices(op::WignerHEOM) = deepcopy(op.indices)

function heom_derivative!(out, W, h::SpectralWignerMoyal, symbol)
    h.W .= W
    mul!(h.Wp, h.forward_p, h.W)
    h.Wp .*= symbol
    mul!(out, h.backward_p, h.Wp)
    return nothing
end

function heom_derivative!(out, W, h::FiniteDifferenceWignerMoyal, matrix)
    mul!(out, W, transpose(matrix))
    return nothing
end

function check_hierarchy_size(U, op)
    expected = (size(op.grid)..., length(op.indices))
    size(U) == expected ||
        throw(DimensionMismatch("hierarchy must have size $expected, got $(size(U))"))
    return nothing
end

"""
    heom!(dU, U, op, t)

Evaluate the Wigner HEOM in place. `U` and `dU` have shape `(nq, np, number_of_ADOs)`
and must not alias. The time is unused. Only the root is a normalised Wigner function;
auxiliaries encode bath correlations and may have nonzero integrals.
"""
function heom!(dU, U, op::WignerHEOM, t)
    check_hierarchy_size(U, op)
    check_hierarchy_size(dU, op)
    h, b = op.hamiltonian, op.bath
    for a in eachindex(op.indices)
        W, dW = @view(U[:, :, a]), @view(dU[:, :, a])
        wigner_moyal!(dW, W, h, t)
        @. dW -= op.damping[a] * W
        if !iszero(b.diffusion)
            heom_derivative!(op.scratch, W, h, op.momentum.second)
            @. dW += b.diffusion * op.scratch
        end
    end
    # Differentiate each source once, then scatter its upward/downward contributions.
    for a in eachindex(op.indices)
        W = @view U[:, :, a]
        heom_derivative!(op.scratch, W, h, op.momentum.first)
        for k in eachindex(b.rates)
            below, above = op.lower[k, a], op.upper[k, a]
            if below != 0
                dW = @view dU[:, :, below]
                dW .+= op.scratch
            end
            if above != 0
                dW = @view dU[:, :, above]
                n = op.indices[above][k]
                cr, ci = real(b.coefficients[k]), op.coordinate_coupling[k]
                @. dW += n * (cr * op.scratch + ci * op.grid.q * W)
            end
        end
    end
    return nothing
end

"""
    heom_problem(W0, tspan, grid; mass, potential, bath, depth, kwargs...)
    heom_problem(W0, tspan, op)

Construct an in-place SciML `ODEProblem` for [`heom_operator`](@ref). An `nq × np`
Wigner matrix initialises the physical root and sets every auxiliary to zero: a
factorized initial system and bare bath state. With a counterterm this includes the
physical initial-slip transient. An entire 3D hierarchy may instead be supplied to
restart a calculation or initialise a correlated state. The input is copied to a
real, finite `Float64` array, without renormalisation.

Use `physical_wigner(sol)` to extract the final root. Explicit solvers such as `Vern7()`
work; high Matsubara rates and large depth can make the hierarchy stiff. FFTW buffers
cannot accept dual numbers and no sparse Jacobian is provided.
"""
function heom_problem(U0::AbstractArray, tspan, grid::PhaseSpaceGrid; kwargs...)
    return heom_problem(U0, tspan, heom_operator(grid; kwargs...))
end

function heom_problem(U0::AbstractArray, tspan, op::WignerHEOM)
    eltype(U0) <: Real || throw(ArgumentError("initial hierarchy must be real"))
    all(isfinite, U0) || throw(ArgumentError("initial hierarchy must be finite"))
    if ndims(U0) == 2
        check_size(U0, op.grid)
        U = zeros(size(op.grid)..., length(op.indices))
        U[:, :, 1] .= U0
    else
        check_hierarchy_size(U0, op)
        U = Array{Float64,3}(U0)
    end
    all(isfinite, U) || throw(ArgumentError("initial hierarchy must be finite in Float64"))
    return ODEProblem(heom!, U, tspan, op)
end

"""
    physical_wigner(U::AbstractArray{<:Real,3})
    physical_wigner(sol[, index])

View of the physical (first) hierarchy member, usable with all matrix-and-grid
observables and plotting functions. For a HEOM solution, select a saved state with
`index` (the final state by default). Mutating the view changes the underlying state.
"""
function physical_wigner(U::AbstractArray{<:Real,3})
    size(U, 3) > 0 || throw(ArgumentError("hierarchy needs a physical member"))
    return @view U[:, :, 1]
end

function physical_wigner(sol::AbstractODESolution, index::Integer = lastindex(sol.u))
    sol.prob.p isa WignerHEOM ||
        throw(ArgumentError("solution must come from heom_problem"))
    return physical_wigner(sol.u[index])
end

# The existing solution diagnostics/plots consume physical states only.
physical_states(states, ::AbstractPhaseSpaceOperator) = states
physical_states(states, ::WignerHEOM) = physical_wigner.(states)
physical_state(W, ::AbstractPhaseSpaceOperator) = W
physical_state(U, ::WignerHEOM) = physical_wigner(U)

function Base.show(io::IO, op::WignerHEOM)
    print(
        io,
        "WignerHEOM(",
        op.grid,
        ", depth = ",
        op.depth,
        ", modes = ",
        length(op.bath.rates),
        ", members = ",
        length(op.indices),
        ")",
    )
    return nothing
end
