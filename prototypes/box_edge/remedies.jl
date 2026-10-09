# Candidate remedies as wrappers around the package operator (prototype).
#
#   absorb: dWₙ/dt -= σ(q) Wₙ for every auxiliary n ≠ 0 (never the root, so the trace of
#           the physical state is untouched). σ rises smoothly from 0 at |q| = q₀ to σ₀
#           at the box edge.
#   taper:  the upward coordinate coupling uses q̃(q) instead of q. q̃ = q for |q| ≤ q₀
#           and falls smoothly to 0 at the box edge, which also removes the jump of q
#           across the periodic boundary.
include("common.jl")

ramp(x) = x <= 0 ? 0.0 : x >= 1 ? 1.0 : x^3 * (10 - 15x + 6x^2)   # C² smoothstep

struct Remedied{O}
    op::O
    σ::Vector{Float64}       # per-q auxiliary damping
    δq::Vector{Float64}      # q̃ - q
end

edge_distance(q, q0, L) = (abs(q) - q0) / (L - q0)

function remedied(op; absorb = 0.0, taper = false, q0 = 0.5maximum(abs, op.grid.q))
    q = op.grid.q
    L = maximum(abs, q)
    σ = [absorb * ramp(edge_distance(x, q0, L)) for x in q]
    δq = taper ? [x * (1 - ramp(edge_distance(x, q0, L))) - x for x in q] : zeros(length(q))
    return Remedied(op, σ, δq)
end

function remedied_heom!(dU, U, r::Remedied, t)
    op = r.op
    heom!(dU, U, op, t)
    any(!iszero, r.σ) && for a in 2:length(op.indices)
        @views dU[:, :, a] .-= r.σ .* U[:, :, a]
    end
    any(!iszero, r.δq) && for a in eachindex(op.indices), k in eachindex(op.bath.rates)
        above = op.upper[k, a]
        above == 0 && continue
        ci = op.raising_coordinate[k, a]
        @views dU[:, :, above] .+= ci .* r.δq .* U[:, :, a]
    end
    return nothing
end

# Same for local-symbol analysis.
remedy_profiles(r::Remedied) = (
    taper = x -> x + r.δq[argmin(abs.(r.op.grid.q .- x))],
    absorb = x -> r.σ[argmin(abs.(r.op.grid.q .- x))])

function dense_generator(r::Remedied)
    op = r.op
    L = dense_generator(op)
    N = prod(size(op.grid))
    σ = repeat(r.σ, length(op.grid.p))
    δq = repeat(r.δq, length(op.grid.p))
    block(a) = ((a-1)*N+1):(a*N)
    for a in 2:length(op.indices)
        L[block(a), block(a)] .-= Diagonal(σ)
    end
    for a in eachindex(op.indices), k in eachindex(op.bath.rates)
        above = op.upper[k, a]
        above == 0 && continue
        L[block(above), block(a)] .+= op.raising_coordinate[k, a] .* Diagonal(δq)
    end
    return L
end

function growth_rate(r::Remedied; kwargs...)
    return growth_rate(r.op; rhs! = (dU, U, _, t) -> remedied_heom!(dU, U, r, t), kwargs...)
end

remedied_problem(W0, tspan, r::Remedied) = HEOM.ODEProblem(
    remedied_heom!, cat(W0, zeros(size(W0)..., length(r.op.indices) - 1); dims = 3), tspan, r)
