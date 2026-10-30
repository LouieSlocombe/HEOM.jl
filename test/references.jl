# Independent cold, strongly coupled quantum Brownian-motion references.
# Hänggi and Ingold, Chaos 15, 026105 (2005), Eqs. (5), (9), (11)-(14):
# https://arxiv.org/pdf/quant-ph/0412052
# Fleming, Roura and Hu, Ann. Phys. 326, 1207 (2011):
# https://arxiv.org/pdf/1004.1603

# A linear generalized Langevin representation of a specified exponential bath.
# The force covariance is SIGNED: individual quantum correlation coefficients may
# be negative. The auxiliary matrices are algebraic covariance representations,
# not covariances of independently sampleable classical noises. Their physical
# 2×2 marginal is positive for a sufficiently accurate physical bath expansion.
function cold_memory_matrices(bath; mass, omega)
    modes = length(bath.rates)
    dimension = modes + 2
    A, Q, initial = zeros(dimension, dimension),
    zeros(dimension, dimension),
    zeros(dimension, dimension)
    A[1, 2] = 1 / mass
    A[2, 1] = -mass * omega^2 - 2bath.counterterm
    Q[2, 2] = 2bath.diffusion
    damping = Matrix(Diagonal(bath.rates) + bath.mixing)
    A[2, 3:end] = bath.weights
    A[3:end, 1] = -2imag.(bath.coefficients) / bath.hbar
    A[3:end, 3:end] = -damping
    residues = real.(bath.coefficients)
    if all(iszero, bath.mixing) && all(isone, bath.weights)
        force_covariance = Matrix(Diagonal(residues))
    else
        # For a confluent (polynomial-exponential) or weighted basis, any symmetric S with
        # S*w=real(c) generates exactly w' exp(-Γt)c. This explicit solution
        # makes the reference independent of a sampled stochastic bath and
        # remains finite at a Drude/thermal-pole collision.
        w = bath.weights
        norm_squared = dot(w, w)
        force_covariance =
            (residues * w' + w * residues') / norm_squared -
            dot(w, residues) * (w * w') / norm_squared^2
    end
    initial[3:end, 3:end] = force_covariance
    Q[3:end, 3:end] = damping * force_covariance + force_covariance * damping'
    # AΣ + ΣA' + Q = 0. This avoids exponentially growing Van Loan blocks when
    # low-temperature decompositions include large rates or long trajectories.
    stationary = lyap(A, Q)
    return A, initial, stationary
end

function cold_memory_reference(t, mean, covariance, bath; mass, omega)
    A, initial, stationary = cold_memory_matrices(bath; mass, omega)
    initial[1:2, 1:2] = covariance
    F = exp(t * A)
    mean_t = F * [mean; zeros(length(bath.rates))]
    covariance_t = stationary + F * (initial - stationary) * F'
    return mean_t[1:2], Symmetric(covariance_t[1:2, 1:2])
end

# Direct continuum FDT integral, with NO Matsubara/Padé or hierarchy expansion.
# χ(Ω) = [m(ω²-Ω²) - iΩ 2λ/(γ-iΩ)]⁻¹ for the counterterm convention used here.
# Σqq = ħ/π ∫ coth(ħΩ/2kT) Imχ dΩ; Σpp inserts m²Ω².
# Gauss-Legendre quadrature on Ω=tan(πu/2) covers the whole positive axis.
function drude_fdt_covariance(;
    mass,
    omega,
    reorganization,
    cutoff,
    kT,
    hbar,
    quadrature = 512,
)
    N = quadrature
    rule = eigen(SymTridiagonal(zeros(N), [j / sqrt(4j^2 - 1) for j in 1:(N-1)]))
    variances = zeros(2)
    for (u, weight) in zip((rule.values .+ 1) ./ 2, rule.vectors[1, :] .^ 2)
        frequency = tan(π * u / 2)
        susceptibility = inv(
            mass * (omega^2 - frequency^2) -
            im * frequency * 2reorganization / (cutoff - im * frequency),
        )
        integrand =
            hbar / π * imag(susceptibility) / tanh(hbar * frequency / (2kT)) *
            (π / 2) *
            (1 + frequency^2)
        variances .+= weight .* [integrand, (mass * frequency)^2 * integrand]
    end
    return Diagonal(variances)
end

# ∫₀^∞ C(t) dt of the retained exponentials, from the definition of the basis.
static_correlation(bath) =
    dot(bath.weights, (Diagonal(bath.rates) + Matrix(bath.mixing)) \ bath.coefficients)
