# Batched evaluation of the Wigner HEOM right-hand side for the whole hierarchy at once.
#
# The per-member loop of `heom_rhs!(::Nothing, …)` is the reference. It launches a handful
# of small kernels per member, which leaves a GPU idle. Here the members mix only through
# real matrices acting on the third axis, and those commute with the momentum transform, so
# the Moyal potential term, ∂p and ∂p² share one forward and one backward p-transform, and
# the whole hierarchy needs four batched FFTs per evaluation. With `A`, `B`, `T` the gather
# tables built by `gather_tables`:
#
#   Ûp  = rfft_p(U)
#   Yp  = (M + D s2) Ûp + s1 A Ûp                      (kernel `hierarchy_mix!`)
#   dU  = vel Tq + brfft_p(Yp) - γ U + q (B U) + T U   (kernel `hierarchy_combine!`)
#
# where Tq is the q-derivative term. The p axis is moved to the front with a transposed copy
# so every FFT runs over the contiguous first dimension of a batch.

# Compressed rows of a sparse N × N matrix: output member `a` gathers the inputs
# `colidx[p]` with weights `values[p]` for `p in rowptr[a]:rowptr[a+1]-1`.
struct GatherTable{I<:AbstractVector,V<:AbstractVector}
    rowptr::I
    colidx::I
    values::V
end

struct BatchedHierarchy{
    T<:AbstractFloat,
    V<:AbstractVector,
    CV<:AbstractVector,
    CM<:AbstractMatrix,
    K<:AbstractArray,
    G<:GatherTable,
    R3<:AbstractArray,
    C3<:AbstractArray,
    FQ,
    BQ,
    FP,
    BP,
}
    diffusion::T
    potential_symbol::CM     # Moyal symbol in (κ, q) order, refreshed for driven potentials
    first_symbol::CV         # ∂p multiplier, length np ÷ 2 + 1
    second_symbol::V         # ∂p² multiplier
    kinetic_symbol::K        # iκ_q / nq, shaped (nq ÷ 2 + 1, 1, 1)
    velocity::V              # -p / m
    q::V
    damping::V               # Σ nₖνₖ per member
    derivative_table::G      # lowering and raising_derivative
    coordinate_table::G      # raising_coordinate
    transfer_table::G        # bath mixing
    qwork::R3                # copy of the state for the q-FFT, then the q-derivative
    pwork::R3                # transposed state (np, nq, N), then the transposed p-terms
    qhat::C3
    phat::C3
    mixed::C3
    forward_q::FQ
    backward_q::BQ
    forward_p::FP
    backward_p::BP
end

"""
    hierarchy_mix!(batch)

Kernel of the batched right-hand side (defined by the KernelAbstractions extension):
combine the transformed members `batch.phat` into `batch.mixed`.
"""
function hierarchy_mix! end

"""
    hierarchy_combine!(dU, U, batch)

Kernel of the batched right-hand side (defined by the KernelAbstractions extension):
assemble `dU` from the transformed derivative terms in `batch` and the members of `U`.
"""
function hierarchy_combine! end

function require_backend(extension)
    extension === nothing && throw(
        ArgumentError(
            "batched hierarchies need KernelAbstractions.jl: run `using KernelAbstractions` " *
            "(GPU packages such as CUDA.jl load it)",
        ),
    )
    return nothing
end

function gather_tables(op::WignerHEOM)
    members = length(op.indices)
    rows, columns, values = Int[], Int[], Float64[]
    coordinate_rows, coordinate_columns, coordinate_values = Int[], Int[], Float64[]
    for a in 1:members, k in eachindex(op.bath.rates)
        below, above = op.lower[k, a], op.upper[k, a]
        if below != 0
            push!(rows, below)
            push!(columns, a)
            push!(values, op.lowering[k, a])
        end
        if above != 0
            push!(rows, above)
            push!(columns, a)
            push!(values, op.raising_derivative[k, a])
            push!(coordinate_rows, above)
            push!(coordinate_columns, a)
            push!(coordinate_values, op.raising_coordinate[k, a])
        end
    end
    derivative = sparse(rows, columns, values, members, members)
    coordinate =
        sparse(coordinate_rows, coordinate_columns, coordinate_values, members, members)
    return derivative, coordinate, op.transfer
end

# Row-compressed copy of `M` (the column-compressed form of its transpose).
function gather_table(device, ::Type{T}, M::SparseMatrixCSC) where {T}
    rows = copy(transpose(M))
    index(v) = to_device(device, Int32, v)
    return GatherTable(
        index(rows.colptr),
        index(rows.rowval),
        to_device(device, T, rows.nzval),
    )
end

function build_batched(
    op::WignerHEOM,
    h::SpectralWignerMoyal,
    device,
    ::Type{T},
) where {T<:AbstractFloat}
    require_backend(Base.get_extension(@__MODULE__, :HEOMKernelAbstractionsExt))
    host = op.hamiltonian
    nq, np = size(op.grid)
    members = length(op.indices)
    derivative, coordinate, transfer = gather_tables(op)
    to(x) = to_device(device, T, x)
    symbol = to(permutedims(host.potential_symbol))
    kinetic = reshape(to(host.kinetic_symbol), :, 1, 1)
    real_3d(dims...) = fill!(similar(h.velocity, T, dims...), 0)
    complex_3d(dims...) = fill!(similar(symbol, Complex{T}, dims...), 0)
    qwork = real_3d(nq, np, members)
    pwork = real_3d(np, nq, members)
    qhat = complex_3d(nq ÷ 2 + 1, np, members)
    phat = complex_3d(np ÷ 2 + 1, nq, members)
    return BatchedHierarchy(
        T(op.bath.diffusion),
        symbol,
        to(vec(op.momentum.first)),
        to(vec(op.momentum.second)),
        kinetic,
        to(vec(host.velocity)),
        to(op.grid.q),
        to(op.damping),
        gather_table(device, T, derivative),
        gather_table(device, T, coordinate),
        gather_table(device, T, transfer),
        qwork,
        pwork,
        qhat,
        phat,
        similar(phat),
        plan_rfft(qwork, 1),
        plan_brfft(qhat, nq, 1),
        plan_rfft(pwork, 1),
        plan_brfft(phat, np, 1),
    )
end

# Moyal symbol in the batch's (κ, q) layout, for a drive that changed it.
function refresh_symbol!(batch::BatchedHierarchy, h::SpectralWignerMoyal)
    permutedims!(batch.potential_symbol, h.potential_symbol, (2, 1))
    return nothing
end

function heom_rhs!(batch::BatchedHierarchy, dU, U, op::WignerHEOM)
    h = op.hamiltonian
    h.drive === nothing || refresh_symbol!(batch, h)
    copyto!(batch.qwork, U)
    mul!(batch.qhat, batch.forward_q, batch.qwork)
    batch.qhat .*= batch.kinetic_symbol
    mul!(batch.qwork, batch.backward_q, batch.qhat)
    permutedims!(batch.pwork, U, (2, 1, 3))
    mul!(batch.phat, batch.forward_p, batch.pwork)
    hierarchy_mix!(batch)
    mul!(batch.pwork, batch.backward_p, batch.mixed)
    hierarchy_combine!(dU, U, batch)
    return nothing
end

# Apply a Moyal-type symbol in (κ, q) order to every member: the time derivative of a drive.
function batched_symbol_action!(out, U, batch::BatchedHierarchy, symbol)
    permutedims!(batch.pwork, U, (2, 1, 3))
    mul!(batch.phat, batch.forward_p, batch.pwork)
    batch.phat .*= reshape(symbol, size(symbol, 1), size(symbol, 2), 1)
    mul!(batch.pwork, batch.backward_p, batch.phat)
    permutedims!(out, batch.pwork, (2, 1, 3))
    return nothing
end

# Symbol and scalar factor of ∂H/∂t acting on the members, in the operator's precision.
function tgrad_symbol(h::SpectralWignerMoyal, drive::SeparableMoyalDrive, t)
    factor = -field_time_derivative(drive.potential, t, real(eltype(h.potential_symbol)))
    return permutedims(drive.coupling), factor
end
function tgrad_symbol(h::SpectralWignerMoyal, drive::GeneralMoyalDrive, t)
    symbol =
        moyal_symbol(
            q -> potential_time_derivative(drive.potential, q, t),
            h.grid.q,
            wavenumbers(length(h.grid.p), h.grid.dp),
            h.hbar,
            h.moyal_terms,
        ) ./ length(h.grid.p)
    check_finite(symbol)
    return permutedims(array_like(h.potential_symbol, symbol)),
    one(real(eltype(h.potential_symbol)))
end

function heom_tgrad_rhs!(batch::BatchedHierarchy, out, U, op::WignerHEOM, t)
    symbol, factor = tgrad_symbol(op.hamiltonian, op.hamiltonian.drive, t)
    batched_symbol_action!(out, U, batch, symbol)
    out .*= factor
    return nothing
end
