# Frozen-point hierarchies in density-matrix variables (ħ = 1). At fixed x = q + y/2 and
# x' = q - y/2 every bath coupling is a scalar; the kinetic energy is the only term
# that couples different (x, x'). For real rates νₖ and coefficients cₖ:
#
#   :package   K modes,  Ẇₙ = -n·ν Wₙ - iyΣWₙ₊ₑ - iΣnₖ(cₖx - cₖ*x')Wₙ₋ₑ        (heom.jl)
#   :standard  2K modes, left/right split, Gatto et al. Eq. (19)
#   :gatto     2K modes, Gaussian + Bogoliubov transformed, Gatto et al. Eq. (60) with
#              γᵢ = 0. The sine sectors (N, P) are then invariant and start at zero, so
#              they are dropped. Physical state: Σ w_q R_q over all-even occupations.
#
# Every form also carries the Markovian remainder -D y² on each member.
using LinearAlgebra, Printf

function simplex(modes, depth)
    out = [zeros(Int, modes)]
    for total in 1:depth
        stack = [(Int[], total)]
        level = Vector{Int}[]
        function rec(prefix, left)
            if length(prefix) == modes - 1
                push!(level, [prefix; left])
                return
            end
            for v in left:-1:0
                rec([prefix; v], left - v)
            end
        end
        rec(Int[], total)
        append!(out, level)
    end
    return out, Dict(n => i for (i, n) in enumerate(out))
end

function frozen_matrix(form, c, ν, D, x, xp, depth)
    K = length(ν)
    y, q = x - xp, (x + xp) / 2
    modes = form === :package ? K : 2K
    idx, lookup = simplex(modes, depth)
    H = zeros(ComplexF64, length(idx), length(idx))
    at(n) = all(>=(0), n) ? get(lookup, n, 0) : 0
    shift(n, pairs...) = (m = copy(n); for (i, d) in pairs; m[i] += d; end; m)
    for (r, n) in enumerate(idx)
        H[r, r] -= D * y^2
        if form === :package
            for k in 1:K
                H[r, r] -= n[k] * ν[k]
                (s = at(shift(n, k => 1))) != 0 && (H[r, s] += -im * y)
                (s = at(shift(n, k => -1))) != 0 && (H[r, s] += -im * n[k] * (c[k] * x - conj(c[k]) * xp))
            end
        elseif form === :standard
            for k in 1:K
                u, v = k, K + k
                H[r, r] -= (n[u] + n[v]) * ν[k]
                for j in (u, v)
                    (s = at(shift(n, j => 1))) != 0 && (H[r, s] += -im * y)
                end
                (s = at(shift(n, u => -1))) != 0 && (H[r, s] += -im * n[u] * c[k] * x)
                (s = at(shift(n, v => -1))) != 0 && (H[r, s] += im * n[v] * conj(c[k]) * xp)
            end
        elseif form === :gatto
            for k in 1:K
                M, O = k, K + k
                H[r, r] -= (n[M] + n[O]) * ν[k]
                for j in (M, O)
                    (s = at(shift(n, j => 2))) != 0 && (H[r, s] += ν[k])
                    (s = at(shift(n, j => 1))) != 0 && (H[r, s] += -im * sqrt(2) * y)
                end
                # -(i/√2)[ηx(M R₋ - R₊) - η*x'(O R₋ - R₊)]
                (s = at(shift(n, M => -1))) != 0 && (H[r, s] += -im / sqrt(2) * c[k] * x * n[M])
                (s = at(shift(n, M => 1))) != 0 && (H[r, s] += im / sqrt(2) * c[k] * x)
                (s = at(shift(n, O => -1))) != 0 && (H[r, s] += im / sqrt(2) * conj(c[k]) * xp * n[O])
                (s = at(shift(n, O => 1))) != 0 && (H[r, s] += -im / sqrt(2) * conj(c[k]) * xp)
            end
        end
    end
    weights = zeros(length(idx))
    for (r, n) in enumerate(idx)
        if form === :gatto
            all(iseven, n) && (weights[r] = prod(1 / (2.0^(m ÷ 2) * factorial(m ÷ 2)) for m in n))
        else
            weights[r] = r == 1 ? 1.0 : 0.0
        end
    end
    return H, weights
end

# Exact frozen-path influence functional: ρ(t) = exp(-y Σ(cₖx - cₖ*x')Aₖ(t) - D y² t),
# Aₖ(t) = t/νₖ - (1 - e^{-νₖt})/νₖ².
function exact_frozen(c, ν, D, x, xp, t)
    y = x - xp
    A = @. t / ν - (1 - exp(-ν * t)) / ν^2
    return exp(-y * sum(@. (c * x - conj(c) * xp) * A) - D * y^2 * t)
end

function frozen_root(form, c, ν, D, x, xp, depth, t)
    H, w = frozen_matrix(form, c, ν, D, x, xp, depth)
    ψ = zeros(ComplexF64, size(H, 1)); ψ[1] = 1
    return dot(w, exp(t * H) * ψ)
end

frozen_abscissa(form, c, ν, D, x, xp, depth) =
    maximum(real, eigvals(frozen_matrix(form, c, ν, D, x, xp, depth)[1]))

# Package form closed with the free scalar Markovian closure using every parent
# (:markov_full in prototypes/nz_terminator/closure.jl), in unscaled frozen variables:
# a tier-(N+1) child m is Xₘ = Σⱼ mⱼ Φⱼ W_{m-eⱼ}/(m·ν) and returns -iy Xₘ to each parent.
function markov_frozen_matrix(c, ν, D, x, xp, depth)
    H, w = frozen_matrix(:package, c, ν, D, x, xp, depth)
    K = length(ν)
    idx, lookup = simplex(K, depth)
    y = x - xp
    Φ = [-im * (c[j] * x - conj(c[j]) * xp) for j in 1:K]
    children = Set{Vector{Int}}()
    for n in idx
        sum(n) == depth || continue
        for k in 1:K
            m = copy(n); m[k] += 1; push!(children, m)
        end
    end
    for m in children
        z = sum(m .* ν)
        parents = [(j, lookup[(p = copy(m); p[j] -= 1; p)]) for j in 1:K if m[j] > 0]
        for (_, target) in parents, (j, source) in parents
            H[target, source] += -im * y * m[j] * Φ[j] / z
        end
    end
    return H, w
end
