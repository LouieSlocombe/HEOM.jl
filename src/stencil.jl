# Finite differences without a sparse matrix, for devices and for other precisions.
#
# The sparse operator of `FiniteDifference` is a sum of Kronecker products of circulant
# bands. Applied to the grid it is a sum of one-dimensional stencils: the q-derivative
# weighted by -p/m, and the odd p-derivatives weighted by the per-position Moyal
# coefficients. A `StencilTable` holds the bands, `stencil_moyal!` and its relatives
# (defined by the KernelAbstractions extension) apply them in one pass, and the operators
# below mirror `FiniteDifferenceWignerMoyal` and `FiniteDifferenceCaldeiraLeggett`.

# Several circulant stencils on one axis: stencil `s` occupies entries
# `pointers[s]:pointers[s+1]-1`, and applies `(Df)[k] = Σ weights[p] f[mod1(k + offsets[p], n)]`.
struct StencilTable{I<:AbstractVector,V<:AbstractVector}
    pointers::I
    offsets::I
    weights::V
end

# Stencils for the given derivative orders on an axis of spacing `h`. The grid is known to
# have enough points, since the sparse operator was built from the same weights.
function stencil_table(device, ::Type{T}, derivatives, order, h) where {T}
    pointers, offsets, weights = Int[1], Int[], Float64[]
    for derivative in derivatives
        band = central_difference_weights(derivative, order)
        for (offset, weight) in zip((-(length(band)÷2)):(length(band)÷2), band)
            iszero(weight) && continue
            push!(offsets, offset)
            push!(weights, Float64(weight) / h^derivative)
        end
        push!(pointers, length(offsets) + 1)
    end
    return StencilTable(
        to_device(device, Int32, pointers),
        to_device(device, Int32, offsets),
        to_device(device, T, weights),
    )
end

# Matrix-free finite-difference Wigner–Moyal operator: the stencils of `kinetic` (along q)
# and `potential` (the odd p-derivatives) weighted by `velocity` and `coefficients`.
struct StencilWignerMoyal{V<:AbstractVector,R<:AbstractMatrix,K<:StencilTable,D} <:
       AbstractWignerMoyal
    grid::PhaseSpaceGrid
    mass::Float64
    hbar::Float64
    moyal_terms::Int
    order::Int
    velocity::V                # -p / m
    coefficients::R            # length(q) × moyal_terms, refreshed by a drive
    kinetic::K
    potential::K
    drive::D
end

struct StencilCaldeiraLeggett{
    H<:StencilWignerMoyal,
    V<:AbstractVector,
    K<:StencilTable,
    T<:AbstractFloat,
} <: AbstractCaldeiraLeggett
    grid::PhaseSpaceGrid
    mass::Float64
    hbar::Float64
    friction::Float64
    kT::Float64
    hamiltonian::H
    momentum::V                # p
    first::K
    second::K
    friction_rate::T
    diffusion::T
end

"""
    stencil_operator(op, device, T)

Matrix-free copy of a `FiniteDifference` operator on the array type `device` in precision
`T`. This is what `move_operator` returns for a finite-difference operator that leaves the
host or the `Float64` precision.
"""
function stencil_operator(
    op::FiniteDifferenceWignerMoyal,
    device,
    ::Type{T},
) where {T<:AbstractFloat}
    require_backend(Base.get_extension(@__MODULE__, :HEOMKernelAbstractionsExt))
    nq, np = size(op.grid)
    return StencilWignerMoyal(
        op.grid,
        op.mass,
        op.hbar,
        op.moyal_terms,
        op.order,
        to_device(device, T, -op.grid.p ./ op.mass),
        to_device(device, T, op.coefficients),
        stencil_table(device, T, (1,), op.order, op.grid.dq),
        stencil_table(device, T, 1:2:(2op.moyal_terms-1), op.order, op.grid.dp),
        stencil_drive(op, op.drive, device, T),
    )
end

stencil_drive(op, ::Nothing, device, ::Type{T}) where {T} = nothing
stencil_drive(op, drive::GeneralMoyalDrive, device, ::Type{T}) where {T} = drive
function stencil_drive(op, drive::SeparableMoyalDrive, device, ::Type{T}) where {T}
    coupling = stencil_coefficients(op, drive.potential.dipole)
    return SeparableMoyalDrive(
        drive.potential,
        to_device(device, T, op.coefficients),
        to_device(device, T, coupling),
    )
end

# Coefficients of the potential `V` for the momentum derivatives of this operator.
function stencil_coefficients(op, V)
    return moyal_coefficients(moyal_derivatives(V, op.grid.q, op.moyal_terms), op.hbar)
end

function stencil_operator(
    op::FiniteDifferenceCaldeiraLeggett,
    device,
    ::Type{T},
) where {T<:AbstractFloat}
    hamiltonian = stencil_operator(op.hamiltonian, device, T)
    diffusion = op.mass * op.friction * op.kT
    return StencilCaldeiraLeggett(
        op.grid,
        op.mass,
        op.hbar,
        op.friction,
        op.kT,
        hamiltonian,
        to_device(device, T, op.grid.p),
        stencil_table(device, T, (1,), op.order, op.grid.dp),
        stencil_table(device, T, (2,), op.order, op.grid.dp),
        T(op.friction),
        T(diffusion),
    )
end

# Finite differences keep the sparse matrices on the host in `Float64`.
move_operator(op::FiniteDifferenceWignerMoyal, device, ::Type{T}) where {T} =
    stencil_operator(op, device, T)
move_operator(op::FiniteDifferenceCaldeiraLeggett, device, ::Type{T}) where {T} =
    stencil_operator(op, device, T)

element_type(op::StencilWignerMoyal) = eltype(op.velocity)
state_prototype(op::StencilWignerMoyal) = op.velocity

# `stencil_moyal!`, `stencil_derivative!`, `stencil_bath!` and `hierarchy_stencil!` are
# implemented by the KernelAbstractions extension.

"""
    stencil_moyal!(dW, W, op, coefficients, kinetic)

Kernel (defined by the KernelAbstractions extension): `dW = L W` for the stencil operator
`op`, with the given per-position `coefficients` of the potential terms, and the kinetic
term included if `kinetic`. `W` and `dW` are a matrix or a hierarchy of members.
"""
function stencil_moyal! end

"""
    stencil_derivative!(out, W, table)

Kernel (defined by the KernelAbstractions extension): apply the first stencil of `table`
along the momentum axis of every member of `W`.
"""
function stencil_derivative! end

"""
    stencil_bath!(dW, W, op)

Kernel (defined by the KernelAbstractions extension): add the Caldeira–Leggett friction
and diffusion terms of `op` to `dW`.
"""
function stencil_bath! end

"""
    hierarchy_stencil!(dU, U, op, batch)

Kernel (defined by the KernelAbstractions extension): add the hierarchy terms of a
stencil operator to `dU`, which already holds the Hamiltonian part.
"""
function hierarchy_stencil! end

function apply_wigner_moyal!(dW, W, op::StencilWignerMoyal)
    check_size(W, op.grid)
    check_size(dW, op.grid)
    stencil_moyal!(dW, W, op, op.coefficients, true)
    return nothing
end

function update_generator!(op::StencilWignerMoyal, drive::SeparableMoyalDrive, t)
    E = checked_drive_value(drive.potential.field(user_time(t)), eltype(op.coefficients))
    @. op.coefficients = drive.baseline - E * drive.coupling
    return nothing
end

function update_generator!(op::StencilWignerMoyal, drive::GeneralMoyalDrive, t)
    coefficients = stencil_coefficients(op, potential_at(drive.potential, user_time(t)))
    copyto!(op.coefficients, array_like(op.coefficients, coefficients))
    return nothing
end

function potential_action!(out, W, op::StencilWignerMoyal, t = 0.0)
    update_generator!(op, t)
    check_size(out, op.grid)
    check_size(W, op.grid)
    stencil_moyal!(out, W, op, op.coefficients, false)
    return nothing
end

# Coefficients and scale of ∂H/∂t acting on the potential stencils.
function stencil_tgrad(op, drive::SeparableMoyalDrive, t)
    factor = -field_time_derivative(drive.potential, t, eltype(op.coefficients))
    return drive.coupling, factor
end
function stencil_tgrad(op, drive::GeneralMoyalDrive, t)
    coefficients =
        stencil_coefficients(op, q -> potential_time_derivative(drive.potential, q, t))
    return array_like(op.coefficients, coefficients), one(eltype(op.coefficients))
end

function wigner_tgrad!(
    out,
    W,
    op::StencilWignerMoyal,
    drive::Union{SeparableMoyalDrive,GeneralMoyalDrive},
    t,
)
    coefficients, factor = stencil_tgrad(op, drive, t)
    stencil_moyal!(out, W, op, coefficients, false)
    out .*= factor
    return nothing
end

function heom_derivative!(out, W, h::StencilWignerMoyal, table)
    stencil_derivative!(out, W, table)
    return nothing
end

function caldeira_leggett!(dW, W, op::StencilCaldeiraLeggett, t)
    wigner_moyal!(dW, W, op.hamiltonian, t)
    stencil_bath!(dW, W, op)
    return nothing
end

element_type(op::StencilCaldeiraLeggett) = element_type(op.hamiltonian)
state_prototype(op::StencilCaldeiraLeggett) = state_prototype(op.hamiltonian)

function no_sparse_matrix()
    throw(
        ArgumentError(
            "a sparse matrix needs the host operator; build it with the default " *
            "device = Array and float_type = Float64",
        ),
    )
end
SparseArrays.sparse(::StencilWignerMoyal) = no_sparse_matrix()
SparseArrays.sparse(::StencilWignerMoyal, ::Real) = no_sparse_matrix()
SparseArrays.sparse(::StencilCaldeiraLeggett) = no_sparse_matrix()
SparseArrays.sparse(::StencilCaldeiraLeggett, ::Real) = no_sparse_matrix()

# Momentum stencils of the hierarchy for the finite-difference operator on a device.
function move_momentum(h::StencilWignerMoyal, op::WignerHEOM, device, ::Type{T}) where {T}
    return HEOMMomentum(
        stencil_table(device, T, (1,), h.order, op.grid.dp),
        stencil_table(device, T, (2,), h.order, op.grid.dp),
    )
end

# The hierarchy terms of a stencil operator in one kernel, on top of the Hamiltonian part.
struct BatchedStencilHierarchy{T<:AbstractFloat,V<:AbstractVector,G<:GatherTable}
    diffusion::T
    q::V
    damping::V
    derivative_table::G
    coordinate_table::G
    transfer_table::G
end

function build_batched(
    op::WignerHEOM,
    h::StencilWignerMoyal,
    device,
    ::Type{T},
) where {T<:AbstractFloat}
    derivative, coordinate, transfer = gather_tables(op)
    return BatchedStencilHierarchy(
        T(op.bath.diffusion),
        to_device(device, T, op.grid.q),
        to_device(device, T, op.damping),
        gather_table(device, T, derivative),
        gather_table(device, T, coordinate),
        gather_table(device, T, transfer),
    )
end

function heom_rhs!(batch::BatchedStencilHierarchy, dU, U, op::WignerHEOM)
    stencil_moyal!(dU, U, op.hamiltonian, op.hamiltonian.coefficients, true)
    hierarchy_stencil!(dU, U, op, batch)
    return nothing
end

function heom_tgrad_rhs!(::BatchedStencilHierarchy, out, U, op::WignerHEOM, t)
    h = op.hamiltonian
    coefficients, factor = stencil_tgrad(h, h.drive, t)
    stencil_moyal!(out, U, h, coefficients, false)
    out .*= factor
    return nothing
end

function response_seed!(::BatchedStencilHierarchy, seed, U, kick)
    stencil_moyal!(seed, U, kick, kick.coefficients, false)
    return nothing
end

function Base.show(io::IO, op::StencilWignerMoyal)
    print(io, "StencilWignerMoyal(", op.grid, ", order = ", op.order)
    print(io, ", mass = ", op.mass, ", hbar = ", op.hbar)
    print(io, ", moyal_terms = ", op.moyal_terms, ")")
    return nothing
end

function Base.show(io::IO, op::StencilCaldeiraLeggett)
    print(io, "StencilCaldeiraLeggett(", op.grid, ", order = ", op.hamiltonian.order)
    print(io, ", mass = ", op.mass, ", hbar = ", op.hbar)
    print(io, ", friction = ", op.friction, ", kT = ", op.kT)
    print(io, ", moyal_terms = ", op.hamiltonian.moyal_terms, ")")
    return nothing
end
