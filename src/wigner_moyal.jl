# Largest supported truncation of the Moyal series. Nested forward-mode differentiation of
# V⁽²ᴺ⁻¹⁾ carries 2^(2N - 1) partials, and four terms are already exact for polynomial
# potentials up to degree eight.
const MAX_MOYAL_TERMS = 4

# Common interface for Hamiltonian and dissipative phase-space evolution.
abstract type AbstractPhaseSpaceOperator end

"""
    AbstractWignerMoyal

Supertype of the semi-discrete Wigner–Moyal operators returned by
[`wigner_moyal_operator`](@ref).
"""
abstract type AbstractWignerMoyal <: AbstractPhaseSpaceOperator end

# Pseudo-spectral operator. The FFTs are planned on the stored buffers, and only those
# buffers are ever passed to them. The backward transforms are unnormalised, so the symbols
# carry the 1/n factors. The arrays live on the host or on a device (`R` real matrices,
# `CV` complex vectors, `C` complex matrices) in one floating-point type.
struct SpectralWignerMoyal{
    R<:AbstractMatrix,
    CV<:AbstractVector,
    C<:AbstractMatrix,
    FQ,
    BQ,
    FP,
    BP,
    D,
} <: AbstractWignerMoyal
    grid::PhaseSpaceGrid
    mass::Float64
    hbar::Float64
    moyal_terms::Union{Nothing,Int}
    velocity::R                          # 1 × np row of -p/m
    kinetic_symbol::CV                   # iκ_q / nq
    potential_symbol::C                  # Moyal symbol M(q, κ_p) / np
    W::R
    Wq::C
    Wp::C
    Tq::R
    Tp::R
    forward_q::FQ
    backward_q::BQ
    forward_p::FP
    backward_p::BP
    drive::D
end

# Finite-difference operator, stored as a sparse matrix acting on vec(W). The coefficients
# cₛ V⁽²ˢ⁺¹⁾(qᵢ) of the momentum derivatives (`length(q) × moyal_terms`) are kept so that
# the same operator can be rebuilt as matrix-free stencils on a device.
struct FiniteDifferenceWignerMoyal{D} <: AbstractWignerMoyal
    grid::PhaseSpaceGrid
    mass::Float64
    hbar::Float64
    moyal_terms::Int
    order::Int
    matrix::SparseMatrixCSC{Float64,Int}
    coefficients::Matrix{Float64}
    drive::D
end

"""
    wigner_moyal_operator(grid; mass, potential, hbar = 1.0, discretization = Spectral(),
                          moyal_terms = nothing, device = Array, float_type = Float64)

Semi-discretise the Wigner–Moyal right-hand side on `grid` for a particle of `mass` in the
`potential` `V(q)`, and return the operator used by [`wigner_moyal!`](@ref).
Use `V(q,t)`, [`TimeDependentPotential`](@ref), or [`DrivenPotential`](@ref) for a
time-dependent potential. Separable drives precompute their spatial coefficients.
`mass` and `hbar` must be finite and positive in Float64, and `V` must return finite real
values at every evaluation point.

With `moyal_terms = nothing` (the default) the whole Moyal series is kept. In momentum
Fourier space, where `∂/∂p` becomes `iκ`, the potential term is then exactly

    i [V(q + ħκ/2) - V(q - ħκ/2)] / ħ.

This needs no derivatives of `V`, but evaluates it up to `πħ/(2dp)` beyond the box.
`moyal_terms = N`, with `1 ≤ N ≤ 4`, keeps only the terms `s = 0, …, N - 1` instead.
`N = 1` is the classical Liouville equation. `N = 2` adds `-(ħ²/24) V‴ ∂³W/∂p³` and is
exact for potentials up to quartic. ForwardDiff takes the derivatives of `V`, so the
potential must accept dual numbers.

`discretization` is [`Spectral`](@ref)`()` or [`FiniteDifference`](@ref)`(order)`; the
latter needs an integer `moyal_terms`. Strongly anharmonic potentials make the system stiff
because the symbol grows like `ħ² V‴ κ³`.

`device` is the array type that holds the operator and its states (`Array`, or for example
`CuArray` once a GPU package and KernelAbstractions are loaded) and `float_type` is
`Float64` or `Float32`. The symbols are evaluated on the host in `Float64` and rounded once,
so `Float32` loses no accuracy to the setup. A `FiniteDifference` operator that leaves the
host `Float64` case is applied as matrix-free stencils instead of a sparse matrix, which
needs KernelAbstractions (see [`stencil_operator`](@ref)); `sparse(op)` is then unavailable.

The operator holds work buffers, so give each parallel task its own operator.
"""
function wigner_moyal_operator(
    grid::PhaseSpaceGrid;
    mass::Real,
    potential,
    hbar::Real = 1.0,
    discretization::Union{Spectral,FiniteDifference} = Spectral(),
    moyal_terms::Union{Nothing,Integer} = nothing,
    device = Array,
    float_type = Float64,
)
    mass, hbar = Float64(mass), Float64(hbar)
    isfinite(mass) && mass > 0 ||
        throw(ArgumentError("mass must be finite and positive, got $mass"))
    isfinite(hbar) && hbar > 0 ||
        throw(ArgumentError("hbar must be finite and positive, got $hbar"))
    moyal_terms === nothing ||
        1 <= moyal_terms <= MAX_MOYAL_TERMS ||
        throw(
            ArgumentError(
                "moyal_terms must be between 1 and $MAX_MOYAL_TERMS, or nothing for the " *
                "exact operator, got $moyal_terms",
            ),
        )
    op = build_operator(
        discretization,
        grid,
        Float64(mass),
        Float64(hbar),
        normalize_potential(potential),
        moyal_terms,
    )
    return move_operator(op, check_device(device), check_float_type(float_type))
end

function build_operator(::Spectral, grid, mass, hbar, potential, moyal_terms)
    nq, np = size(grid)
    W = zeros(nq, np)
    Wq = zeros(ComplexF64, nq ÷ 2 + 1, np)
    Wp = zeros(ComplexF64, nq, np ÷ 2 + 1)
    kinetic_symbol = im .* wavenumbers(nq, grid.dq) ./ nq
    velocity = reshape(-grid.p ./ mass, 1, np)
    all(isfinite, velocity) && all(isfinite, kinetic_symbol) ||
        throw(ArgumentError("kinetic coefficients must be finite; rescale mass or grid"))
    κp = wavenumbers(np, grid.dp)
    potential_symbol = moyal_symbol(potential, grid.q, κp, hbar, moyal_terms) ./ np
    check_finite(potential_symbol)
    return SpectralWignerMoyal(
        grid,
        mass,
        hbar,
        moyal_terms,
        velocity,
        kinetic_symbol,
        potential_symbol,
        W,
        Wq,
        Wp,
        zeros(nq, np),
        zeros(nq, np),
        plan_rfft(W, 1),
        plan_brfft(Wq, nq, 1),
        plan_rfft(W, 2),
        plan_brfft(Wp, np, 2),
        nothing,
    )
end

function build_operator(fd::FiniteDifference, grid, mass, hbar, potential, moyal_terms)
    moyal_terms === nothing && throw(
        ArgumentError(
            "FiniteDifference needs a truncated Moyal series, for example " *
            "moyal_terms = 2, because the exact operator is nonlocal in momentum",
        ),
    )
    nq, np = size(grid)
    derivatives = moyal_derivatives(potential, grid.q, moyal_terms)
    # vec(A * W * transpose(B)) == kron(B, A) * vec(W): A acts along q and B along p.
    coefficients = moyal_coefficients(derivatives, hbar)
    L = kron(spdiagm(-grid.p ./ mass), periodic_difference_matrix(1, fd.order, nq, grid.dq))
    for s in 0:(moyal_terms-1)
        Dp = periodic_difference_matrix(2s + 1, fd.order, np, grid.dp)
        L += kron(Dp, spdiagm(coefficients[:, s+1]))
    end
    all(isfinite, L.nzval) || throw(
        ArgumentError(
            "finite-difference coefficients must be finite; rescale parameters or grid",
        ),
    )
    return FiniteDifferenceWignerMoyal(
        grid,
        mass,
        hbar,
        moyal_terms,
        fd.order,
        L,
        coefficients,
        nothing,
    )
end

# Operators are built on the host in Float64; `move_operator` converts one to the requested
# array type and precision, with fresh buffers and FFT plans on the target arrays.
move_operator(op::SpectralWignerMoyal, ::Type{Array}, ::Type{Float64}) = op
move_operator(op::FiniteDifferenceWignerMoyal, ::Type{Array}, ::Type{Float64}) = op

function move_operator(op::SpectralWignerMoyal, device, ::Type{T}) where {T<:AbstractFloat}
    nq, np = size(op.grid)
    velocity = to_device(device, T, op.velocity)
    kinetic_symbol = to_device(device, T, op.kinetic_symbol)
    potential_symbol = to_device(device, T, op.potential_symbol)
    W = fill!(similar(velocity, T, nq, np), 0)
    Wq = fill!(similar(potential_symbol, Complex{T}, nq ÷ 2 + 1, np), 0)
    Wp = fill!(similar(potential_symbol, Complex{T}, nq, np ÷ 2 + 1), 0)
    return SpectralWignerMoyal(
        op.grid,
        op.mass,
        op.hbar,
        op.moyal_terms,
        velocity,
        kinetic_symbol,
        potential_symbol,
        W,
        Wq,
        Wp,
        fill!(similar(W), 0),
        fill!(similar(W), 0),
        plan_rfft(W, 1),
        plan_brfft(Wq, nq, 1),
        plan_rfft(W, 2),
        plan_brfft(Wp, np, 2),
        move_drive(op.drive, device, T),
    )
end

element_type(op::SpectralWignerMoyal) = eltype(op.velocity)
element_type(op::FiniteDifferenceWignerMoyal) = Float64
state_prototype(op::SpectralWignerMoyal) = op.velocity
state_prototype(op::FiniteDifferenceWignerMoyal) = Matrix{Float64}(undef, 0, 0)

# Fourier symbol of the Moyal potential term: M[i, j] multiplies the mode exp(iκⱼp) of
# W(qᵢ, ⋅).
function moyal_symbol(V, q, κ, hbar, ::Nothing)
    return ComplexF64[
        im * (potential_value(V, x + hbar * k / 2) - potential_value(V, x - hbar * k / 2)) /
        hbar for x in q, k in κ
    ]
end

# A complex potential is not a Hermitian Hamiltonian; the real FFT would silently discard
# part of its dynamics. Check V itself too: differentiating an infinite constant can give 0.
function potential_value(V, q)
    value = V(q)
    value isa Real && isfinite(value) ||
        throw(ArgumentError("potential must return finite real values, got $value at q=$q"))
    return value
end

# Truncated series Σₛ cₛ V⁽²ˢ⁺¹⁾(q) (iκ)²ˢ⁺¹, the Fourier form of Σₛ cₛ V⁽²ˢ⁺¹⁾ ∂²ˢ⁺¹/∂p²ˢ⁺¹.
function moyal_symbol(V, q, κ, hbar, moyal_terms::Integer)
    derivatives = moyal_derivatives(V, q, moyal_terms)
    symbol = zeros(ComplexF64, length(q), length(κ))
    for s in 0:(moyal_terms-1)
        c = moyal_coefficient(s, hbar)
        symbol .+= c .* derivatives[:, s+1] .* transpose((im .* κ) .^ (2s + 1))
    end
    return symbol
end

# Coefficient cₛ = (-1)ˢ (ħ/2)²ˢ / (2s + 1)! of V⁽²ˢ⁺¹⁾ ∂²ˢ⁺¹W/∂p²ˢ⁺¹ in the Moyal series.
moyal_coefficient(s::Integer, hbar::Real) = (-1)^s * (hbar / 2)^(2s) / factorial(2s + 1)

# The coefficients cₛ V⁽²ˢ⁺¹⁾(qᵢ) of the odd momentum derivatives, as a `length(q) × terms`
# matrix, from the derivatives of the potential.
function moyal_coefficients(derivatives, hbar)
    return derivatives .*
           transpose([moyal_coefficient(s, hbar) for s in axes(derivatives, 2) .- 1])
end

# Odd derivatives V⁽²ˢ⁺¹⁾(qᵢ) for s = 0, …, terms - 1, as a `length(q) × terms` matrix.
function moyal_derivatives(V, q, terms)
    foreach(x -> potential_value(V, x), q)
    derivatives = [Float64(potential_derivative(V, x, 2s + 1)) for x in q, s in 0:(terms-1)]
    check_finite(derivatives)
    return derivatives
end

# k-th derivative of V at q by nested forward-mode automatic differentiation.
function potential_derivative(V, q::Real, k::Integer)
    k == 0 && return V(q)
    return ForwardDiff.derivative(x -> potential_derivative(V, x, k - 1), q)
end

function check_finite(values)
    return all(isfinite, values) || throw(
        ArgumentError(
            "the Moyal potential term is not finite on this grid; the exact operator " *
            "also evaluates the potential up to πħ/(2dp) beyond the box",
        ),
    )
end

"""
    wigner_moyal!(dW, W, op, t)

Evaluate the semi-discrete Wigner–Moyal right-hand side `dW = ∂W/∂t` in place, for the
Wigner function `W` on the grid of the operator `op` from [`wigner_moyal_operator`](@ref).
The potential is evaluated at `t` when it depends on time. The signature is that of an
in-place `ODEProblem` function whose
parameter is `op`.
"""
function wigner_moyal!(dW, W, op::AbstractWignerMoyal, t)
    update_generator!(op, t)
    return apply_wigner_moyal!(dW, W, op)
end

function apply_wigner_moyal!(dW, W, op::SpectralWignerMoyal)
    check_size(W, op.grid)
    check_size(dW, op.grid)
    op.W .= W
    mul!(op.Wq, op.forward_q, op.W)
    op.Wq .*= op.kinetic_symbol
    mul!(op.Tq, op.backward_q, op.Wq)
    mul!(op.Wp, op.forward_p, op.W)
    op.Wp .*= op.potential_symbol
    mul!(op.Tp, op.backward_p, op.Wp)
    @. dW = op.velocity * op.Tq + op.Tp
    return nothing
end

function apply_wigner_moyal!(dW, W, op::FiniteDifferenceWignerMoyal)
    check_size(W, op.grid)
    check_size(dW, op.grid)
    mul!(vec(dW), op.matrix, vec(W))
    return nothing
end

phase_space_rhs!(dW, W, op::AbstractWignerMoyal, t) = wigner_moyal!(dW, W, op, t)

"""
    wigner_moyal_problem(W0, tspan, grid::PhaseSpaceGrid; mass, potential, kwargs...)
    wigner_moyal_problem(W0, tspan, op)

Method-of-lines `ODEProblem` for the Wigner–Moyal equation, starting from the Wigner
function `W0` (an `nq × np` matrix on the grid) over `tspan`.

The first form builds the operator with [`wigner_moyal_operator`](@ref), passing on the
keywords. The second reuses an existing operator `op`, which becomes the problem's
parameter.

The spectral right-hand side cannot take dual numbers, and the finite-difference one has no
sparse Jacobian set up yet. Use an explicit solver, for example `Tsit5()` or `Vern9()`.

# Examples

```julia
using HEOM, OrdinaryDiffEqVerner

grid = PhaseSpaceGrid((-8.0, 8.0), 64, (-8.0, 8.0), 64)
W0 = on_grid((q, p) -> coherent_wigner(q, p; q0 = 2.0, mass = 1.0, omega = 1.0), grid)
V = harmonic_potential(; mass = 1.0, omega = 1.0)
prob = wigner_moyal_problem(W0, (0.0, 2π), grid; mass = 1.0, potential = V)
sol = solve(prob, Vern9(); abstol = 1e-10, reltol = 1e-10)
W = sol.u[end] # back to W0 after one period
```
"""
function wigner_moyal_problem(W0::AbstractMatrix, tspan, grid::PhaseSpaceGrid; kwargs...)
    return wigner_moyal_problem(W0, tspan, wigner_moyal_operator(grid; kwargs...))
end

function wigner_moyal_problem(W0::AbstractMatrix, tspan, op::AbstractWignerMoyal)
    eltype(W0) <: Real || throw(ArgumentError("initial Wigner function must be real"))
    size(W0) == size(op.grid) || throw(
        DimensionMismatch(
            "initial Wigner function has size $(size(W0)), but the grid has size " *
            "$(size(op.grid))",
        ),
    )
    T = element_type(op)
    initial = state_like(state_prototype(op), T, W0)
    all(isfinite, initial) || throw(ArgumentError("initial Wigner function must be finite"))
    return ODEProblem(wigner_moyal!, initial, convert_tspan(T, tspan), op)
end

"""
    sparse(op::FiniteDifferenceWignerMoyal)
    sparse(op::FiniteDifferenceWignerMoyal, t)

Return a copy of the finite-difference Wigner–Moyal operator as a sparse matrix `L` acting on
`vec(W)`, so that `vec(dW) == L * vec(W)`. `L` is exactly skew-symmetric, so the discrete
norm `sum(W)` and purity `sum(abs2, W)` are conserved.
A driven operator requires the explicit time argument.
"""
function SparseArrays.sparse(op::FiniteDifferenceWignerMoyal)
    is_time_dependent(op) &&
        throw(ArgumentError("use sparse(op, t) for a driven generator"))
    return copy(op.matrix)
end

function SparseArrays.sparse(op::FiniteDifferenceWignerMoyal, t::Real)
    update_generator!(op, t)
    return copy(op.matrix)
end

function Base.show(io::IO, op::SpectralWignerMoyal)
    print(io, "SpectralWignerMoyal(", op.grid, ", mass = ", op.mass, ", hbar = ", op.hbar)
    print(io, ", moyal_terms = ", something(op.moyal_terms, "exact"), ")")
    return nothing
end

function Base.show(io::IO, op::FiniteDifferenceWignerMoyal)
    print(io, "FiniteDifferenceWignerMoyal(", op.grid, ", order = ", op.order)
    print(io, ", mass = ", op.mass, ", hbar = ", op.hbar)
    print(io, ", moyal_terms = ", op.moyal_terms, ")")
    return nothing
end
