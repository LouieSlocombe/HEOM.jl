# Shared tools for the box-edge stability study (prototype, not part of the package).
#
# Every hierarchy coupling in heom! is local in the mixed representation (q, κ), where κ
# is conjugate to p: ∂p becomes a multiplier s₁(κ) and q stays a multiplier. The Moyal
# potential term is a purely imaginary multiplier there, the same on every member, so it
# moves no real parts. Only the kinetic term -(p/m)∂q couples different (q, κ). Freezing
# (q, κ) therefore leaves a small members×members matrix H(q, κ): its eigenvalues are
# the "local symbol" of the hierarchy, exact in the limit m → ∞.
using HEOM, LinearAlgebra, SparseArrays, Printf, Random
using OrdinaryDiffEqVerner
include(joinpath(@__DIR__, "..", "nz_terminator", "gaussian_reference.jl"))

const HARMONIC = harmonic_potential(; mass = 1.0, omega = 1.0)

const BATHS = (
    warm = drude_lorentz_pade_bath(; reorganization = 0.2, cutoff = 1.0, kT = 1.0, pade = 1),
    cold_weak = drude_lorentz_pade_bath(; reorganization = 0.05, cutoff = 0.5, kT = 0.1, pade = 2),
    cold_strong = drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 2),
    cold_strong1 = drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 1),
)

box(L, n) = PhaseSpaceGrid((-L, L), n, (-L, L), n)

# Multipliers of the first and second momentum derivative on the nonnegative
# wavenumbers; the operators are real, so the negative ones are complex conjugates.
function momentum_symbols(op)
    np = length(op.grid.p)
    m = op.momentum
    if op.hamiltonian isa HEOM.SpectralWignerMoyal
        return vec(m.first) .* np, vec(m.second) .* np
    end
    j = 0:(np ÷ 2)
    phase = [cis(2π * (c - 1) * k / np) for k in j, c in 1:np]
    return phase * Vector(m.first[1, :]), phase * Vector(m.second[1, :])
end

# Frozen-coefficient hierarchy matrix at one (q, κ). `extra_damping` is added to every
# auxiliary (not the root); `coordinate` replaces q in the upward coupling.
function local_matrix(op, coordinate, s1, s2; extra_damping = 0.0)
    b = op.bath
    M = length(op.indices)
    H = Matrix{ComplexF64}(op.transfer)
    for a in 1:M
        H[a, a] += -op.damping[a] + b.diffusion * s2 - (a > 1 ? extra_damping : 0.0)
        for k in eachindex(b.rates)
            below, above = op.lower[k, a], op.upper[k, a]
            below != 0 && (H[below, a] += op.lowering[k, a] * s1)
            above != 0 &&
                (H[above, a] += op.raising_derivative[k, a] * s1 +
                                op.raising_coordinate[k, a] * coordinate)
        end
    end
    return H
end

"""
Largest real part of the frozen hierarchy matrix over every grid q and κ ≥ 0.
Returns (rate, q, κ index). `taper(q)` and `absorb(q)` model remedies.
"""
function local_abscissa(op; taper = identity, absorb = q -> 0.0, qs = op.grid.q)
    s1, s2 = momentum_symbols(op)
    best = (-Inf, NaN, 0)
    for q in qs, j in eachindex(s1)
        H = local_matrix(op, taper(q), s1[j], s2[j]; extra_damping = absorb(q))
        r = maximum(real, eigvals!(H))
        r > best[1] && (best = (r, q, j))
    end
    return best
end

# Dense generator acting on vec(U), for either discretisation.
function dense_generator(op)
    if op.hamiltonian isa HEOM.FiniteDifferenceWignerMoyal
        return Matrix(sparse(op))
    end
    shape = (size(op.grid)..., length(op.indices))
    N = prod(shape)
    L = zeros(N, N)
    U, dU = zeros(shape), zeros(shape)
    for j in 1:N
        U[j] = 1
        heom!(dU, U, op, 0.0)
        L[:, j] .= vec(dU)
        U[j] = 0
    end
    return L
end

"""
Leading eigenvalues of a dense generator, with the fraction of each eigenvector's
squared weight in the outer quarter of q and of p (|x| > L/2).
"""
function leading_modes(L, op; count = 3)
    F = eigen(L)
    order = sortperm(real.(F.values); rev = true)[1:count]
    nq, np = size(op.grid)
    Lq = maximum(abs, op.grid.q)
    Lp = maximum(abs, op.grid.p)
    outer_q = abs.(op.grid.q) .> Lq / 2
    outer_p = abs.(op.grid.p) .> Lp / 2
    return map(order) do i
        v = reshape(abs2.(F.vectors[:, i]), nq, np, :)
        w = dropdims(sum(v; dims = 3); dims = 3)
        total = sum(w)
        (λ = F.values[i], outer_q = sum(w[outer_q, :]) / total,
            outer_p = sum(w[:, outer_p]) / total,
            root = sum(v[:, :, 1]) / sum(v))
    end
end

"""
Asymptotic growth rate of the propagator, measured by integrating a random hierarchy
and renormalising at every interval of length `dt`. Returns the mean log-growth over the
last half of the run (an estimate of the spectral abscissa that includes any
non-normal transient only through its effect on that window).
"""
function growth_rate(op; T = 60.0, dt = 2.0, seed = 1, rhs! = heom!)
    rng = MersenneTwister(seed)
    U = randn(rng, size(op.grid)..., length(op.indices))
    U ./= norm(U)
    rates = Float64[]
    prob = HEOM.ODEProblem(rhs!, U, (0.0, dt), op)
    t = 0.0
    while t < T - 1e-9
        sol = solve(remake(prob; u0 = U), Vern7(); abstol = 1e-10, reltol = 1e-8,
            save_everystep = false, save_start = false)
        U = sol.u[end]
        g = norm(U)
        push!(rates, log(g) / dt)
        U ./= g
        t += dt
    end
    return sum(rates[(end÷2+1):end]) / length(rates[(end÷2+1):end]), U
end

"""
Time at which max|U| first exceeds `factor` times its initial value, or Inf.
"""
function blowup_time(prob; factor = 1e3, abstol = 1e-9, reltol = 1e-9)
    limit = factor * maximum(abs, prob.u0)
    sol = solve(prob, Vern7(); abstol, reltol, save_everystep = false,
        unstable_check = (dt, u, p, t) -> !(maximum(abs, u) < limit))
    return HEOM.successful_retcode(sol) ? Inf : sol.t[end]
end

displaced_gaussian(grid) = on_grid(
    (q, p) -> gaussian_wigner(q, p; mean = [1.0, 0.0], covariance = [0.5 0; 0 0.5]), grid)

function variances(W, grid)
    w = W .* (grid.dq * grid.dp)
    Q, P = grid.q, grid.p'
    n = sum(w)
    mq, mp = sum(w .* Q) / n, sum(w .* P) / n
    return sum(w .* (Q .- mq) .^ 2) / n, sum(w .* (P .- mp) .^ 2) / n
end
