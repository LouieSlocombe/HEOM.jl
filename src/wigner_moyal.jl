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
# carry the 1/n factors.
struct SpectralWignerMoyal{FQ,BQ,FP,BP} <: AbstractWignerMoyal
    grid::PhaseSpaceGrid
    mass::Float64
    hbar::Float64
    moyal_terms::Union{Nothing,Int}
    velocity::Matrix{Float64}            # 1 × np row of -p/m
    kinetic_symbol::Vector{ComplexF64}   # iκ_q / nq
    potential_symbol::Matrix{ComplexF64} # Moyal symbol M(q, κ_p) / np
    W::Matrix{Float64}
    Wq::Matrix{ComplexF64}
    Wp::Matrix{ComplexF64}
    Tq::Matrix{Float64}
    Tp::Matrix{Float64}
    forward_q::FQ
    backward_q::BQ
    forward_p::FP
    backward_p::BP
end

# Finite-difference operator, stored as a sparse matrix acting on vec(W).
struct FiniteDifferenceWignerMoyal <: AbstractWignerMoyal
    grid::PhaseSpaceGrid
    mass::Float64
    hbar::Float64
    moyal_terms::Int
    order::Int
    matrix::SparseMatrixCSC{Float64,Int}
end

"""
    wigner_moyal_operator(grid; mass, potential, hbar = 1.0, discretization = Spectral(),
                          moyal_terms = nothing)

Semi-discretise the Wigner–Moyal right-hand side on `grid` for a particle of `mass` in the
`potential` `V(q)`, and return the operator used by [`wigner_moyal!`](@ref).

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

The operator holds work buffers, so give each parallel task its own operator.
"""
function wigner_moyal_operator(
    grid::PhaseSpaceGrid;
    mass::Real,
    potential,
    hbar::Real = 1.0,
    discretization::Union{Spectral,FiniteDifference} = Spectral(),
    moyal_terms::Union{Nothing,Integer} = nothing,
)
    mass > 0 || throw(ArgumentError("mass must be positive, got $mass"))
    hbar > 0 || throw(ArgumentError("hbar must be positive, got $hbar"))
    moyal_terms === nothing ||
        1 <= moyal_terms <= MAX_MOYAL_TERMS ||
        throw(
            ArgumentError(
                "moyal_terms must be between 1 and $MAX_MOYAL_TERMS, or nothing for the " *
                "exact operator, got $moyal_terms",
            ),
        )
    return build_operator(
        discretization,
        grid,
        Float64(mass),
        Float64(hbar),
        potential,
        moyal_terms,
    )
end

function build_operator(::Spectral, grid, mass, hbar, potential, moyal_terms)
    nq, np = size(grid)
    W = zeros(nq, np)
    Wq = zeros(ComplexF64, nq ÷ 2 + 1, np)
    Wp = zeros(ComplexF64, nq, np ÷ 2 + 1)
    kinetic_symbol = im .* wavenumbers(nq, grid.dq) ./ nq
    κp = wavenumbers(np, grid.dp)
    potential_symbol = moyal_symbol(potential, grid.q, κp, hbar, moyal_terms) ./ np
    check_finite(potential_symbol)
    return SpectralWignerMoyal(
        grid,
        mass,
        hbar,
        moyal_terms,
        reshape(-grid.p ./ mass, 1, np),
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
    L = kron(spdiagm(-grid.p ./ mass), periodic_difference_matrix(1, fd.order, nq, grid.dq))
    for s in 0:(moyal_terms-1)
        Dp = periodic_difference_matrix(2s + 1, fd.order, np, grid.dp)
        L += kron(Dp, spdiagm(moyal_coefficient(s, hbar) .* derivatives[:, s+1]))
    end
    return FiniteDifferenceWignerMoyal(grid, mass, hbar, moyal_terms, fd.order, L)
end

# Fourier symbol of the Moyal potential term: M[i, j] multiplies the mode exp(iκⱼp) of
# W(qᵢ, ⋅).
function moyal_symbol(V, q, κ, hbar, ::Nothing)
    return [im * (V(x + hbar * k / 2) - V(x - hbar * k / 2)) / hbar for x in q, k in κ]
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

# Odd derivatives V⁽²ˢ⁺¹⁾(qᵢ) for s = 0, …, terms - 1, as a `length(q) × terms` matrix.
function moyal_derivatives(V, q, terms)
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
The time `t` is unused. The signature is that of an in-place `ODEProblem` function whose
parameter is `op`.
"""
function wigner_moyal!(dW, W, op::SpectralWignerMoyal, t)
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

function wigner_moyal!(dW, W, op::FiniteDifferenceWignerMoyal, t)
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
    size(W0) == size(op.grid) || throw(
        DimensionMismatch(
            "initial Wigner function has size $(size(W0)), but the grid has size " *
            "$(size(op.grid))",
        ),
    )
    return ODEProblem(wigner_moyal!, Matrix{Float64}(W0), tspan, op)
end

"""
    sparse(op::FiniteDifferenceWignerMoyal)

Return a copy of the finite-difference Wigner–Moyal operator as a sparse matrix `L` acting on
`vec(W)`, so that `vec(dW) == L * vec(W)`. `L` is exactly skew-symmetric, so the discrete
norm `sum(W)` and purity `sum(abs2, W)` are conserved.
"""
SparseArrays.sparse(op::FiniteDifferenceWignerMoyal) = copy(op.matrix)

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
