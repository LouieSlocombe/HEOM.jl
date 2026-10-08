# Final-tier closures for the Wigner HEOM (prototype, not part of the package).
#
# A hard cutoff sets every member at tier depth+1 to zero. The closures here instead
# eliminate that tier adiabatically. In the unscaled convention a child m at tier
# depth+1 obeys
#
#     dXₘ/dt = (L_c - γₘ) Xₘ + Σⱼ mⱼ Φⱼ W_{m-eⱼ},   Φⱼ = real(cⱼ)∂p + 2imag(cⱼ)q/ħ,
#
# where L_c is the Wigner–Moyal generator plus the bath diffusion and γₘ = Σₖ mₖνₖ.
# Setting dXₘ/dt = 0 gives Xₘ = R(γₘ) Σⱼ mⱼ Φⱼ W_{m-eⱼ}, which returns to each
# retained parent m-eₖ as wₖ ∂p Xₘ. The variants differ in R and in which parents
# feed each child:
#
#   :none         Xₘ = 0 (the package's hard cutoff)
#   :markov_diag  R(z) = 1/z, each parent sees only its own contribution
#   :markov_full  R(z) = 1/z, every parent of the child contributes
#   :nz_diag      R(z) = (z - L_c)⁻¹, own contribution only. This is the diagonal
#                 Nakajima–Zwanzig terminator of Fay, JCP 157, 054108 (2022), Eq. (14)
#   :nz_full      R(z) = (z - L_c)⁻¹, every parent. This is the complete second-order
#                 Markovian kernel K of Fay's Eq. (10), equivalent to the
#                 time-derivative truncation.
#
# A diagonal closure is the full one with each child split into one pseudo-child per
# parent, so both share the code below. Mixing (coupled bath modes) is not supported.

using HEOM, LinearAlgebra, SparseArrays

const CLOSURES = (:none, :markov_diag, :markov_full, :nz_diag, :nz_full)

struct Child
    z::Float64
    parents::Vector{Int}           # hierarchy positions of the retained parents
    modes::Vector{Int}             # parent = child - e_mode
    source_factors::Vector{Float64} # mₖ * S_parent / S_child
    target_factors::Vector{Float64} # wₖ * S_child / S_parent
end

struct TerminatedHEOM{F}
    op::HEOM.WignerHEOM
    kind::Symbol
    children::Vector{Child}
    resolvents::Dict{Float64,F}
    final::Vector{Int}             # final-tier hierarchy positions
    slot::Dict{Int,Int}            # position -> index into derivatives
    cr::Vector{Float64}
    ci::Vector{Float64}
    derivatives::Array{Float64,3}  # ∂p W for every final-tier member
    source::Matrix{Float64}
    X::Matrix{Float64}
    DX::Matrix{Float64}
end

uses_resolvent(kind) = kind in (:nz_diag, :nz_full)
is_diagonal(kind) = kind in (:markov_diag, :nz_diag)

# Generator of a single member without its hierarchy damping: L_WM + D ∂p².
function member_generator(op)
    h = op.hamiltonian
    h isa HEOM.FiniteDifferenceWignerMoyal ||
        throw(ArgumentError("the prototype resolvent needs FiniteDifference"))
    nq = length(op.grid.q)
    return h.matrix + op.bath.diffusion * kron(op.momentum.second, spdiagm(ones(nq)))
end

function shifted_factorization(Lc, z)
    return lu(spdiagm(fill(z, size(Lc, 1))) - Lc)
end

function terminate(op::HEOM.WignerHEOM, kind::Symbol)
    kind in CLOSURES || throw(ArgumentError("unknown closure $kind"))
    b = op.bath
    all(iszero, b.mixing.nzval) || throw(ArgumentError("mixing is not supported"))
    K = length(b.rates)
    scales = HEOM.coefficient_scale.(b.coefficients)
    final = [a for (a, n) in enumerate(op.indices) if sum(n) == op.depth]
    lookup = Dict(Tuple(n) => a for (a, n) in enumerate(op.indices))
    groups = Dict{Vector{Int},Vector{Tuple{Int,Int}}}()
    if kind !== :none
        for a in final, k in 1:K
            m = copy(op.indices[a])
            m[k] += 1
            push!(get!(groups, m, Tuple{Int,Int}[]), (a, k))
        end
    end
    children = Child[]
    for (m, links) in groups
        z = sum(m .* b.rates)
        # log S of the child in the scaled convention, consistently from any parent.
        a1, k1 = first(links)
        log_child = op.log_scales[a1] + log(m[k1]) / 2 + log(scales[k1])
        ratio(a) = op.scaled ? exp(op.log_scales[a] - log_child) : 1.0
        pieces = is_diagonal(kind) ? [[link] for link in links] : [links]
        for piece in pieces
            parents, modes = first.(piece), last.(piece)
            push!(
                children,
                Child(
                    z,
                    parents,
                    modes,
                    [m[k] * ratio(a) for (a, k) in piece],
                    [b.weights[k] / ratio(a) for (a, k) in piece],
                ),
            )
        end
    end
    sort!(children; by = c -> (c.z, c.parents))
    Lc = member_generator(op)
    F = typeof(shifted_factorization(Lc, 1.0))
    resolvents = Dict{Float64,F}()
    if uses_resolvent(kind)
        for z in unique(c.z for c in children)
            resolvents[z] = shifted_factorization(Lc, z)
        end
    end
    nq, np = size(op.grid)
    return TerminatedHEOM{F}(
        op,
        kind,
        children,
        resolvents,
        final,
        Dict(a => i for (i, a) in enumerate(final)),
        real.(b.coefficients),
        2 .* imag.(b.coefficients) ./ b.hbar,
        zeros(nq, np, length(final)),
        zeros(nq, np),
        zeros(nq, np),
        zeros(nq, np),
    )
end

function apply_resolvent!(X, cl::TerminatedHEOM, z, S)
    if uses_resolvent(cl.kind)
        ldiv!(vec(X), cl.resolvents[z], vec(S))
    else
        @. X = S / z
    end
    return X
end

function terminated_heom!(dU, U, cl::TerminatedHEOM, t)
    op = cl.op
    heom!(dU, U, op, t)
    cl.kind === :none && return nothing
    h, q, Dp = op.hamiltonian, op.grid.q, op.momentum.first
    for (i, a) in enumerate(cl.final)
        HEOM.heom_derivative!(@view(cl.derivatives[:, :, i]), @view(U[:, :, a]), h, Dp)
    end
    for child in cl.children
        S = cl.source
        fill!(S, 0)
        for (a, j, f) in zip(child.parents, child.modes, child.source_factors)
            W, DW = @view(U[:, :, a]), @view(cl.derivatives[:, :, cl.slot[a]])
            cr, ci = f * cl.cr[j], f * cl.ci[j]
            @. S += cr * DW + ci * q * W
        end
        apply_resolvent!(cl.X, cl, child.z, S)
        HEOM.heom_derivative!(cl.DX, cl.X, h, Dp)
        for (a, g) in zip(child.parents, child.target_factors)
            dW = @view dU[:, :, a]
            @. dW += g * cl.DX
        end
    end
    return nothing
end

solves_per_rhs(cl::TerminatedHEOM) = uses_resolvent(cl.kind) ? length(cl.children) : 0

function terminated_problem(U0, tspan, cl::TerminatedHEOM)
    U = zeros(size(cl.op.grid)..., length(cl.op.indices))
    if ndims(U0) == 2
        U[:, :, 1] .= U0
    else
        U .= U0
    end
    return HEOM.ODEProblem(terminated_heom!, U, tspan, cl)
end

# Dense effective generator acting on vec(U), for small eigenvalue studies.
function dense_generator(cl::TerminatedHEOM)
    op = cl.op
    L = Matrix(sparse(op))
    cl.kind === :none && return L
    nq, np = size(op.grid)
    N = nq * np
    Dp = kron(op.momentum.first, spdiagm(ones(nq)))
    Q = kron(spdiagm(ones(np)), spdiagm(op.grid.q))
    Lc = member_generator(op)
    block(a) = ((a-1)*N+1):(a*N)
    for child in cl.children
        R =
            uses_resolvent(cl.kind) ? inv(Matrix(spdiagm(fill(child.z, N)) - Lc)) :
            Matrix(I / child.z, N, N)
        DR = Matrix(Dp * R)
        for (a, j, f) in zip(child.parents, child.modes, child.source_factors)
            Φ = Matrix(f * (cl.cr[j] * Dp + cl.ci[j] * Q))
            M = DR * Φ
            for (target, g) in zip(child.parents, child.target_factors)
                L[block(target), block(a)] .+= g .* M
            end
        end
    end
    return L
end
