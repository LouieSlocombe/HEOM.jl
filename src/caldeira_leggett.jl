abstract type AbstractCaldeiraLeggett <: AbstractPhaseSpaceOperator end

struct SpectralCaldeiraLeggett{H<:SpectralWignerMoyal} <: AbstractCaldeiraLeggett
    grid::PhaseSpaceGrid
    mass::Float64
    hbar::Float64
    friction::Float64
    kT::Float64
    hamiltonian::H
    momentum::Matrix{Float64}
    friction_symbol::Matrix{ComplexF64}
    diffusion_symbol::Matrix{Float64}
    bath::Matrix{ComplexF64}
end

struct FiniteDifferenceCaldeiraLeggett <: AbstractCaldeiraLeggett
    grid::PhaseSpaceGrid
    mass::Float64
    hbar::Float64
    friction::Float64
    kT::Float64
    moyal_terms::Int
    order::Int
    matrix::SparseMatrixCSC{Float64,Int}
end

"""
    caldeira_leggett_operator(grid; mass, potential, friction, kT, hbar = 1.0,
                             discretization = Spectral(), moyal_terms = nothing)

Semi-discretise the high-temperature, Markovian Caldeira–Leggett equation,

    ∂W/∂t = L_WM W + friction * ∂(pW)/∂p + mass * friction * kT * ∂²W/∂p²,

where `L_WM` is the Hamiltonian operator from [`wigner_moyal_operator`](@ref).
`friction` is the momentum damping rate: `d⟨p⟩/dt = -⟨V′⟩ - friction * ⟨p⟩`.
Some conventions call this rate `2γ`. `kT` is the thermal energy `k_B T`, in the same
units as the Hamiltonian. Both must be finite and nonnegative; zero friction recovers
Wigner–Moyal evolution. `mass` and `hbar` must be finite and positive.

This is the high-temperature Ohmic-bath approximation, not a zero-temperature relaxation
model. For a harmonic well it approaches the classical thermal Gaussian with mean energy
`kT`, centred at the well minimum. Setting `kT = 0` gives formal friction-only evolution,
which need not preserve a physical quantum state. No extra Lindblad diffusion is added.

`discretization` and `moyal_terms` have the same meaning as in
[`wigner_moyal_operator`](@ref). Friction is discretised in conservative flux form so that
the discrete integral of `W` is conserved. The periodic box must contain the evolving
distribution with negligible boundary weight. Diffusion can make fine grids stiff.
The spectral operator holds work buffers; use a separate operator for each parallel task.
"""
function caldeira_leggett_operator(
    grid::PhaseSpaceGrid;
    mass::Real,
    potential,
    friction::Real,
    kT::Real,
    hbar::Real = 1.0,
    discretization::Union{Spectral,FiniteDifference} = Spectral(),
    moyal_terms::Union{Nothing,Integer} = nothing,
)
    m, ħ, γ, θ = Float64(mass), Float64(hbar), Float64(friction), Float64(kT)
    isfinite(m) && m > 0 || throw(ArgumentError("mass must be finite and positive"))
    isfinite(ħ) && ħ > 0 || throw(ArgumentError("hbar must be finite and positive"))
    isfinite(γ) && γ >= 0 || throw(ArgumentError("friction must be finite and nonnegative"))
    isfinite(θ) && θ >= 0 || throw(ArgumentError("kT must be finite and nonnegative"))
    diffusion = m * γ * θ
    isfinite(diffusion) ||
        throw(ArgumentError("momentum diffusion coefficient must be finite"))
    hamiltonian = wigner_moyal_operator(
        grid;
        mass = m,
        potential,
        hbar = ħ,
        discretization,
        moyal_terms,
    )
    return build_caldeira_leggett(hamiltonian, γ, θ, diffusion)
end

function build_caldeira_leggett(op::SpectralWignerMoyal, friction, kT, diffusion)
    np = length(op.grid.p)
    # Unlike odd derivatives, the second derivative retains the even-n Nyquist mode.
    κ = 2π .* collect(rfftfreq(np, 1 / op.grid.dp))
    return SpectralCaldeiraLeggett(
        op.grid,
        op.mass,
        op.hbar,
        friction,
        kT,
        op,
        reshape(copy(op.grid.p), 1, np),
        reshape(im .* friction .* wavenumbers(np, op.grid.dp) ./ np, 1, :),
        reshape(-diffusion .* κ .^ 2 ./ np, 1, :),
        similar(op.Wp),
    )
end

function build_caldeira_leggett(op::FiniteDifferenceWignerMoyal, friction, kT, diffusion)
    nq, np = size(op.grid)
    Dp = periodic_difference_matrix(1, op.order, np, op.grid.dp)
    Dpp = periodic_difference_matrix(2, op.order, np, op.grid.dp)
    # Dp * diag(p) differentiates the flux pW; a discrete product rule would lose norm.
    bath = friction * Dp * spdiagm(op.grid.p) + diffusion * Dpp
    matrix = op.matrix + kron(bath, spdiagm(ones(nq)))
    return FiniteDifferenceCaldeiraLeggett(
        op.grid,
        op.mass,
        op.hbar,
        friction,
        kT,
        op.moyal_terms,
        op.order,
        matrix,
    )
end

"""
    caldeira_leggett!(dW, W, op, t)

Evaluate the Caldeira–Leggett right-hand side in place using an operator returned by
[`caldeira_leggett_operator`](@ref). The time `t` is unused. Both discretisations conserve
the discrete integral of `W`; energy and purity generally change through the bath.
"""
function caldeira_leggett!(dW, W, op::SpectralCaldeiraLeggett, t)
    h = op.hamiltonian
    wigner_moyal!(dW, W, h, t)
    # Reuse the Hamiltonian's planned FFT buffers. Accumulate both bath terms before
    # transforming back, keeping every transform on the buffers its plan was made for.
    h.W .= W
    mul!(h.Wp, h.forward_p, h.W)
    @. op.bath = op.diffusion_symbol * h.Wp
    @. h.W = op.momentum * W
    mul!(h.Wp, h.forward_p, h.W)
    @. h.Wp = op.friction_symbol * h.Wp + op.bath
    mul!(h.Tp, h.backward_p, h.Wp)
    dW .+= h.Tp
    return nothing
end

function caldeira_leggett!(dW, W, op::FiniteDifferenceCaldeiraLeggett, t)
    check_size(W, op.grid)
    check_size(dW, op.grid)
    mul!(vec(dW), op.matrix, vec(W))
    return nothing
end

phase_space_rhs!(dW, W, op::AbstractCaldeiraLeggett, t) = caldeira_leggett!(dW, W, op, t)

"""
    caldeira_leggett_problem(W0, tspan, grid::PhaseSpaceGrid; mass, potential, friction, kT,
                            kwargs...)
    caldeira_leggett_problem(W0, tspan, op)

Method-of-lines `ODEProblem` for Caldeira–Leggett evolution of the `nq × np` initial
Wigner matrix `W0`. The grid form passes keywords to [`caldeira_leggett_operator`](@ref);
the operator form reuses an existing operator. The problem owns a real, finite `Float64`
copy of `W0`.

As with [`wigner_moyal_problem`](@ref), use an explicit solver such as `Vern9()`.
The spectral buffers do not accept dual numbers and no sparse Jacobian is supplied for
finite differences. Strong diffusion may require small time steps.
"""
function caldeira_leggett_problem(
    W0::AbstractMatrix,
    tspan,
    grid::PhaseSpaceGrid;
    kwargs...,
)
    return caldeira_leggett_problem(W0, tspan, caldeira_leggett_operator(grid; kwargs...))
end

function caldeira_leggett_problem(W0::AbstractMatrix, tspan, op::AbstractCaldeiraLeggett)
    check_size(W0, op.grid)
    eltype(W0) <: Real || throw(ArgumentError("initial Wigner function must be real"))
    W = Matrix{Float64}(W0)
    all(isfinite, W) ||
        throw(ArgumentError("initial Wigner function must be finite in Float64"))
    return ODEProblem(caldeira_leggett!, W, tspan, op)
end

"""
    sparse(op::FiniteDifferenceCaldeiraLeggett)

Return a copy of the finite-difference Caldeira–Leggett matrix acting on `vec(W)`.
Its columns sum to zero, conserving the discrete norm. It is generally not skew-symmetric.
"""
SparseArrays.sparse(op::FiniteDifferenceCaldeiraLeggett) = copy(op.matrix)

function Base.show(io::IO, op::SpectralCaldeiraLeggett)
    print(
        io,
        "SpectralCaldeiraLeggett(",
        op.grid,
        ", mass = ",
        op.mass,
        ", hbar = ",
        op.hbar,
    )
    print(io, ", friction = ", op.friction, ", kT = ", op.kT)
    print(io, ", moyal_terms = ", something(op.hamiltonian.moyal_terms, "exact"), ")")
    return nothing
end

function Base.show(io::IO, op::FiniteDifferenceCaldeiraLeggett)
    print(io, "FiniteDifferenceCaldeiraLeggett(", op.grid, ", order = ", op.order)
    print(io, ", mass = ", op.mass, ", hbar = ", op.hbar)
    print(io, ", friction = ", op.friction, ", kT = ", op.kT)
    print(io, ", moyal_terms = ", op.moyal_terms, ")")
    return nothing
end
