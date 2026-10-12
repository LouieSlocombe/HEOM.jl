module HEOMKernelAbstractionsExt

using HEOM: HEOM, BatchedHierarchy, BatchedStencilHierarchy, StencilCaldeiraLeggett
using HEOM: StencilWignerMoyal
using KernelAbstractions: CPU, @Const, @index, @kernel, get_backend, synchronize

# Sum of weighted members: Σₚ values[p] X[…, colidx[p]] over row `a` of a gather table. The
# member is the last index of the 3D array `X`, so one method serves both layouts.
@inline function gather(rowptr, colidx, values, X, i, j, a)
    total = zero(eltype(X))
    lo, hi = rowptr[a], rowptr[a+1] - one(eltype(rowptr))
    for p in lo:hi
        total += values[p] * X[i, j, colidx[p]]
    end
    return total
end

# The `@kernel` shells below are excluded from line coverage: the macro rewrites their line
# information, so Julia attributes no executed lines to them. The `@inline` point functions
# they call hold the code and are covered.

# Momentum-space terms of one point (κ, q, member): the Moyal potential and ∂p² terms of
# the member itself, and ∂p of the members it is coupled to.
@inline function mix_point!(mixed, phat, symbol, s1, s2, D, rowptr, colidx, values, κ, i, a)
    own = (symbol[κ, i] + D * s2[κ]) * phat[κ, i, a]
    coupled = s1[κ] * gather(rowptr, colidx, values, phat, κ, i, a)
    mixed[κ, i, a] = own + coupled
    return nothing
end

# COV_EXCL_START
@kernel function mix_kernel!(
    mixed,
    @Const(phat),
    @Const(symbol),
    @Const(s1),
    @Const(s2),
    D,
    @Const(rowptr),
    @Const(colidx),
    @Const(values),
)
    κ, i, a = @index(Global, NTuple)
    mix_point!(mixed, phat, symbol, s1, s2, D, rowptr, colidx, values, κ, i, a)
end
# COV_EXCL_STOP

# Real-space assembly of one point (q, p, member).
@inline function combine_point!(
    dU,
    U,
    qwork,
    pwork,
    velocity,
    q,
    damping,
    brow,
    bcol,
    bval,
    trow,
    tcol,
    tval,
    i,
    j,
    a,
)
    coordinate = gather(brow, bcol, bval, U, i, j, a)
    transfer = gather(trow, tcol, tval, U, i, j, a)
    dU[i, j, a] =
        velocity[j] * qwork[i, j, a] + pwork[j, i, a] - damping[a] * U[i, j, a] +
        q[i] * coordinate +
        transfer
    return nothing
end

# COV_EXCL_START
@kernel function combine_kernel!(
    dU,
    @Const(U),
    @Const(qwork),
    @Const(pwork),
    @Const(velocity),
    @Const(q),
    @Const(damping),
    @Const(brow),
    @Const(bcol),
    @Const(bval),
    @Const(trow),
    @Const(tcol),
    @Const(tval),
)
    i, j, a = @index(Global, NTuple)
    combine_point!(
        dU,
        U,
        qwork,
        pwork,
        velocity,
        q,
        damping,
        brow,
        bcol,
        bval,
        trow,
        tcol,
        tval,
        i,
        j,
        a,
    )
end
# COV_EXCL_STOP

# Finite-difference stencils. Arrays are a single matrix, possibly a view into a hierarchy,
# or a hierarchy of members along the third axis; `entry` and `store!` hide the difference.
@inline entry(W::AbstractArray{<:Any,2}, i, j, a) = W[i, j]
@inline entry(W::AbstractArray{<:Any,3}, i, j, a) = W[i, j, a]
@inline function store!(W::AbstractArray{<:Any,2}, value, i, j, a)
    W[i, j] = value
    return nothing
end
@inline function store!(W::AbstractArray{<:Any,3}, value, i, j, a)
    W[i, j, a] = value
    return nothing
end

# Σₚ weights[p] W[i, j + offsets[p], a] over the entries `first:last` of a table, with the
# index wrapped periodically into the `n` points of the momentum axis.
@inline function along_p(W, weights, offsets, first, last, n, i, j, a)
    total = zero(eltype(W))
    for p in first:last
        total += weights[p] * entry(W, i, mod1(j + offsets[p], n), a)
    end
    return total
end

@inline function along_q(W, weights, offsets, first, last, n, i, j, a)
    total = zero(eltype(W))
    for p in first:last
        total += weights[p] * entry(W, mod1(i + offsets[p], n), j, a)
    end
    return total
end

# The Wigner–Moyal operator at one point: the kinetic stencil weighted by -p/m, and each
# odd momentum derivative weighted by its coefficient at this position.
@inline function moyal_point!(
    dW,
    W,
    velocity,
    coefficients,
    kptr,
    koff,
    kw,
    pptr,
    poff,
    pw,
    kinetic,
    i,
    j,
    a,
)
    one_index = one(eltype(kptr))
    total = zero(eltype(dW))
    if kinetic
        total +=
            velocity[j] *
            along_q(W, kw, koff, kptr[1], kptr[2] - one_index, size(W, 1), i, j, a)
    end
    for s in 1:size(coefficients, 2)
        derivative =
            along_p(W, pw, poff, pptr[s], pptr[s+1] - one_index, size(W, 2), i, j, a)
        total += coefficients[i, s] * derivative
    end
    store!(dW, total, i, j, a)
    return nothing
end

# COV_EXCL_START
@kernel function moyal_kernel!(
    dW,
    @Const(W),
    @Const(velocity),
    @Const(coefficients),
    @Const(kptr),
    @Const(koff),
    @Const(kw),
    @Const(pptr),
    @Const(poff),
    @Const(pw),
    kinetic,
)
    i, j, a = @index(Global, NTuple)
    moyal_point!(
        dW,
        W,
        velocity,
        coefficients,
        kptr,
        koff,
        kw,
        pptr,
        poff,
        pw,
        kinetic,
        i,
        j,
        a,
    )
end
# COV_EXCL_STOP

@inline function derivative_point!(out, W, ptr, off, w, i, j, a)
    derivative = along_p(W, w, off, ptr[1], ptr[2] - one(eltype(ptr)), size(W, 2), i, j, a)
    store!(out, derivative, i, j, a)
    return nothing
end

# COV_EXCL_START
@kernel function derivative_kernel!(out, @Const(W), @Const(ptr), @Const(off), @Const(w))
    i, j, a = @index(Global, NTuple)
    derivative_point!(out, W, ptr, off, w, i, j, a)
end
# COV_EXCL_STOP

# Friction `γ ∂p(pW)` and diffusion `D ∂p²W` of Caldeira–Leggett, added to `dW`.
@inline function bath_point!(
    dW,
    W,
    momentum,
    friction,
    diffusion,
    ptr1,
    off1,
    w1,
    ptr2,
    off2,
    w2,
    i,
    j,
    a,
)
    np = size(W, 2)
    flux = zero(eltype(dW))
    for p in ptr1[1]:(ptr1[2]-one(eltype(ptr1)))
        k = mod1(j + off1[p], np)
        flux += w1[p] * momentum[k] * entry(W, i, k, a)
    end
    spread = along_p(W, w2, off2, ptr2[1], ptr2[2] - one(eltype(ptr2)), np, i, j, a)
    store!(dW, entry(dW, i, j, a) + friction * flux + diffusion * spread, i, j, a)
    return nothing
end

# COV_EXCL_START
@kernel function bath_kernel!(
    dW,
    @Const(W),
    @Const(momentum),
    friction,
    diffusion,
    @Const(ptr1),
    @Const(off1),
    @Const(w1),
    @Const(ptr2),
    @Const(off2),
    @Const(w2),
)
    i, j, a = @index(Global, NTuple)
    bath_point!(
        dW,
        W,
        momentum,
        friction,
        diffusion,
        ptr1,
        off1,
        w1,
        ptr2,
        off2,
        w2,
        i,
        j,
        a,
    )
end
# COV_EXCL_STOP

# The hierarchy terms of one point, added to the Hamiltonian part already in `dU`.
@inline function hierarchy_point!(
    dU,
    U,
    damping,
    q,
    D,
    ptr1,
    off1,
    w1,
    ptr2,
    off2,
    w2,
    arow,
    acol,
    aval,
    brow,
    bcol,
    bval,
    trow,
    tcol,
    tval,
    i,
    j,
    a,
)
    np = size(U, 2)
    spread = along_p(U, w2, off2, ptr2[1], ptr2[2] - one(eltype(ptr2)), np, i, j, a)
    coupled = zero(eltype(U))
    lo, hi = arow[a], arow[a+1] - one(eltype(arow))
    for p in lo:hi
        derivative =
            along_p(U, w1, off1, ptr1[1], ptr1[2] - one(eltype(ptr1)), np, i, j, acol[p])
        coupled += aval[p] * derivative
    end
    coordinate = gather(brow, bcol, bval, U, i, j, a)
    transfer = gather(trow, tcol, tval, U, i, j, a)
    dU[i, j, a] +=
        -damping[a] * U[i, j, a] + D * spread + coupled + q[i] * coordinate + transfer
    return nothing
end

# COV_EXCL_START
@kernel function hierarchy_stencil_kernel!(
    dU,
    @Const(U),
    @Const(damping),
    @Const(q),
    D,
    @Const(ptr1),
    @Const(off1),
    @Const(w1),
    @Const(ptr2),
    @Const(off2),
    @Const(w2),
    @Const(arow),
    @Const(acol),
    @Const(aval),
    @Const(brow),
    @Const(bcol),
    @Const(bval),
    @Const(trow),
    @Const(tcol),
    @Const(tval),
)
    i, j, a = @index(Global, NTuple)
    hierarchy_point!(
        dU,
        U,
        damping,
        q,
        D,
        ptr1,
        off1,
        w1,
        ptr2,
        off2,
        w2,
        arow,
        acol,
        aval,
        brow,
        bcol,
        bval,
        trow,
        tcol,
        tval,
        i,
        j,
        a,
    )
end
# COV_EXCL_STOP

# Views and reshapes are launched on the backend of the array underneath.
backend_of(x) = parent(x) === x ? get_backend(x) : backend_of(parent(x))

ndrange_of(W) = (size(W, 1), size(W, 2), size(W, 3))

function HEOM.stencil_moyal!(dW, W, op::StencilWignerMoyal, coefficients, kinetic::Bool)
    backend = backend_of(dW)
    k, p = op.kinetic, op.potential
    moyal_kernel!(backend)(
        dW,
        W,
        op.velocity,
        coefficients,
        k.pointers,
        k.offsets,
        k.weights,
        p.pointers,
        p.offsets,
        p.weights,
        kinetic;
        ndrange = ndrange_of(W),
    )
    backend isa CPU && synchronize(backend)
    return nothing
end

function HEOM.stencil_derivative!(out, W, table)
    backend = backend_of(out)
    derivative_kernel!(backend)(
        out,
        W,
        table.pointers,
        table.offsets,
        table.weights;
        ndrange = ndrange_of(W),
    )
    backend isa CPU && synchronize(backend)
    return nothing
end

function HEOM.stencil_bath!(dW, W, op::StencilCaldeiraLeggett)
    backend = backend_of(dW)
    first, second = op.first, op.second
    bath_kernel!(backend)(
        dW,
        W,
        op.momentum,
        op.friction_rate,
        op.diffusion,
        first.pointers,
        first.offsets,
        first.weights,
        second.pointers,
        second.offsets,
        second.weights;
        ndrange = ndrange_of(W),
    )
    backend isa CPU && synchronize(backend)
    return nothing
end

function HEOM.hierarchy_stencil!(dU, U, op, batch::BatchedStencilHierarchy)
    backend = get_backend(dU)
    first, second = op.momentum.first, op.momentum.second
    a, b, t = batch.derivative_table, batch.coordinate_table, batch.transfer_table
    hierarchy_stencil_kernel!(backend)(
        dU,
        U,
        batch.damping,
        batch.q,
        batch.diffusion,
        first.pointers,
        first.offsets,
        first.weights,
        second.pointers,
        second.offsets,
        second.weights,
        a.rowptr,
        a.colidx,
        a.values,
        b.rowptr,
        b.colidx,
        b.values,
        t.rowptr,
        t.colidx,
        t.values;
        ndrange = size(dU),
    )
    backend isa CPU && synchronize(backend)
    return nothing
end

# Kernel launches are asynchronous on a GPU, ordered with the FFTs on its stream. The CPU
# backend needs no ordering, but is synchronised to be safe.
function HEOM.hierarchy_mix!(batch::BatchedHierarchy)
    backend = get_backend(batch.mixed)
    table = batch.derivative_table
    mix_kernel!(backend)(
        batch.mixed,
        batch.phat,
        batch.potential_symbol,
        batch.first_symbol,
        batch.second_symbol,
        batch.diffusion,
        table.rowptr,
        table.colidx,
        table.values;
        ndrange = size(batch.mixed),
    )
    backend isa CPU && synchronize(backend)
    return nothing
end

function HEOM.hierarchy_combine!(dU, U, batch::BatchedHierarchy)
    backend = get_backend(dU)
    coordinate, transfer = batch.coordinate_table, batch.transfer_table
    combine_kernel!(backend)(
        dU,
        U,
        batch.qwork,
        batch.pwork,
        batch.velocity,
        batch.q,
        batch.damping,
        coordinate.rowptr,
        coordinate.colidx,
        coordinate.values,
        transfer.rowptr,
        transfer.colidx,
        transfer.values;
        ndrange = size(dU),
    )
    backend isa CPU && synchronize(backend)
    return nothing
end

end
