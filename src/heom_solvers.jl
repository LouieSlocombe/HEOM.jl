"""
    sparse(op::WignerHEOM)

Assemble the complete, constant HEOM generator for a `FiniteDifference` operator.
It acts on `vec(U)` and includes the selected auxiliary scaling, diffusion and
hierarchy couplings. Storage grows with the number of grid points and auxiliaries;
use the matrix-free Jacobian-vector product from [`heom_problem`](@ref) for large
hierarchies or spectral discretisation.
"""
function SparseArrays.sparse(op::WignerHEOM)
    h = op.hamiltonian
    h isa FiniteDifferenceWignerMoyal || throw(
        ArgumentError("a sparse HEOM generator requires FiniteDifference discretisation"),
    )
    nq, np = size(op.grid)
    members = length(op.indices)
    Iq, Ip = spdiagm(ones(nq)), spdiagm(ones(np))
    Iphase, Iado = spdiagm(ones(nq * np)), spdiagm(ones(members))
    derivative = kron(op.momentum.first, Iq)
    coordinate = kron(Ip, spdiagm(op.grid.q))
    diagonal = h.matrix + op.bath.diffusion * kron(op.momentum.second, Iq)
    L = kron(Iado, diagonal) + kron(op.transfer - spdiagm(op.damping), Iphase)
    for k in eachindex(op.bath.rates)
        below = findall(!iszero, @view op.lower[k, :])
        above = findall(!iszero, @view op.upper[k, :])
        down = sparse(op.lower[k, below], below, op.lowering[k, below], members, members)
        up_derivative = sparse(
            op.upper[k, above],
            above,
            op.raising_derivative[k, above],
            members,
            members,
        )
        up_coordinate = sparse(
            op.upper[k, above],
            above,
            op.raising_coordinate[k, above],
            members,
            members,
        )
        L += kron(down + up_derivative, derivative) + kron(up_coordinate, coordinate)
    end
    all(isfinite, L.nzval) || throw(ArgumentError("HEOM generator must be finite"))
    return L
end

# The system is linear and autonomous. These exact derivatives avoid differentiating
# FFTW buffers, and work with both flattened Krylov vectors and the 3D ODE state.
function heom_jvp!(Jv, v, u, op, t)
    dimensions = (size(op.grid)..., length(op.indices))
    heom!(reshape(Jv, dimensions), reshape(v, dimensions), op, t)
    return nothing
end

function heom_tgrad!(dT, u, op, t)
    fill!(dT, 0)
    return nothing
end

function heom_ode_function(op, jacobian)
    if jacobian === :matrixfree
        return ODEFunction(heom!; jvp = heom_jvp!, tgrad = heom_tgrad!)
    elseif jacobian === :sparse
        L = sparse(op)
        function jac!(J, u, p, t)
            p === op || throw(
                ArgumentError(
                    "rebuild heom_problem when changing a cached sparse operator",
                ),
            )
            copyto!(J, L)
            return nothing
        end
        return ODEFunction(
            heom!;
            jac = jac!,
            jac_prototype = L,
            jvp = heom_jvp!,
            tgrad = heom_tgrad!,
        )
    end
    throw(ArgumentError("jacobian must be :matrixfree or :sparse"))
end
