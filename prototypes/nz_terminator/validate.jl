# Consistency checks for closure.jl. Run with
#     julia --project=prototypes/nz_terminator prototypes/nz_terminator/validate.jl
#
# 1. :nz_full at depth d is exact adiabatic elimination of tier d+1: inserting the
#    eliminated children into a depth-(d+1) hierarchy reproduces the closed tier-d
#    derivative, and makes every child stationary.
# 2. The dense effective generator matches the matrix-free right-hand side.
# 3. With a single bath mode every child has one parent, so diagonal == full.
using Printf, Random
include("closure.jl")

Random.seed!(1)
grid = PhaseSpaceGrid((-6.0, 6.0), 20, (-6.0, 6.0), 22)
V(q) = q^2 / 2 + 0.05q^4
common = (; mass = 1.0, potential = V, discretization = FiniteDifference(4), moyal_terms = 2)
bath = drude_lorentz_pade_bath(; reorganization = 0.6, cutoff = 0.7, kT = 0.15, pade = 2)
println("bath rates ", bath.rates, "\nresidues ", bath.coefficients)

for scaled in (false, true), depth in (0, 1, 3)
    op = heom_operator(grid; common..., bath, depth, scaled)
    big_op = heom_operator(grid; common..., bath, depth = depth + 1, scaled)
    cl = terminate(op, :nz_full)
    U = randn(size(grid)..., length(op.indices))
    dU = similar(U)
    terminated_heom!(dU, U, cl, 0.0)

    # Rebuild the children explicitly and embed them in the deeper hierarchy.
    position = Dict(Tuple(n) => a for (a, n) in enumerate(big_op.indices))
    Ubig = zeros(size(grid)..., length(big_op.indices))
    for (a, n) in enumerate(op.indices)
        Ubig[:, :, position[Tuple(n)]] .= U[:, :, a]
    end
    for child in cl.children
        cl_source = zeros(size(grid))
        for (a, j, f) in zip(child.parents, child.modes, child.source_factors)
            W = U[:, :, a]
            DW = W * transpose(op.momentum.first)
            cl_source .+= f .* (cl.cr[j] .* DW .+ cl.ci[j] .* grid.q .* W)
        end
        X = similar(cl_source)
        apply_resolvent!(X, cl, child.z, cl_source)
        a, k = first(child.parents), first(child.modes)
        m = copy(op.indices[a])
        m[k] += 1
        Ubig[:, :, position[Tuple(m)]] .= X
    end
    dUbig = similar(Ubig)
    heom!(dUbig, Ubig, big_op, 0.0)
    retained = maximum(
        maximum(abs, dU[:, :, a] - dUbig[:, :, position[Tuple(n)]]) for
        (a, n) in enumerate(op.indices)
    )
    children = maximum(
        maximum(abs, dUbig[:, :, b]) for
        (b, n) in enumerate(big_op.indices) if sum(n) == depth + 1
    )
    scale = maximum(abs, dU)
    @printf(
        "scaled=%-5s depth=%d  retained mismatch %.2e  child residual %.2e  (|dU| %.2e)\n",
        scaled,
        depth,
        retained,
        children,
        scale
    )
    @assert retained < 1e-9 * scale && children < 1e-9 * scale
end

op = heom_operator(grid; common..., bath, depth = 2, scaled = true)
U = randn(size(grid)..., length(op.indices))
for kind in CLOSURES
    cl = terminate(op, kind)
    dU = similar(U)
    terminated_heom!(dU, U, cl, 0.0)
    mismatch = maximum(abs, dense_generator(cl) * vec(U) - vec(dU))
    @printf("%-12s dense vs matrix-free %.2e  children %d\n", kind, mismatch, length(cl.children))
    @assert mismatch < 1e-9 * maximum(abs, dU)
end

single = ExponentialBath(bath.coefficients[1:1], bath.rates[1:1]; hbar = 1.0)
op = heom_operator(grid; common..., bath = single, depth = 3)
U = randn(size(grid)..., length(op.indices))
for (diag, full) in ((:markov_diag, :markov_full), (:nz_diag, :nz_full))
    a, b = similar(U), similar(U)
    terminated_heom!(a, U, terminate(op, diag), 0.0)
    terminated_heom!(b, U, terminate(op, full), 0.0)
    @printf("one mode: %s == %s to %.1e\n", diag, full, maximum(abs, a - b))
    @assert a ≈ b
end
println("all checks passed")
