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

@testset "Cold, strong Drude bath: coupled equilibrium from continuum FDT" begin
    physical = (;
        mass = 1.0,
        omega = 1.0,
        hbar = 1.0,
        reorganization = 0.8,
        cutoff = 0.5,
        kT = 0.1,
    )
    # kT/(ħω)=0.1, with zero-frequency friction 2λ/(mγ)=3.2ω.
    exact = drude_fdt_covariance(; physical...)
    @test exact ≈ drude_fdt_covariance(; physical..., quadrature = 256) atol = 1e-10
    @test det(exact) >= physical.hbar^2 / 4
    # At strong coupling the reduced equilibrium is not the isolated Gibbs state.
    bare_variance = physical.hbar / (2tanh(physical.hbar / (2physical.kT)))
    @test exact[1, 1] < 0.9bare_variance
    @test exact[2, 2] > 1.2bare_variance
    parameters = (;
        reorganization = physical.reorganization,
        cutoff = physical.cutoff,
        kT = physical.kT,
        hbar = physical.hbar,
    )
    grid = PhaseSpaceGrid((-5, 5), 40, (-6, 6), 48)
    reference =
        on_grid((q, p) -> gaussian_wigner(q, p; mean = zeros(2), covariance = exact), grid)
    for (constructor, orders) in
        ((drude_lorentz_bath, (2, 8, 32)), (drude_lorentz_pade_bath, (2, 4, 8)))
        errors = map(orders) do order
            bath =
                constructor === drude_lorentz_bath ?
                constructor(; parameters..., matsubara = order) :
                constructor(; parameters..., pade = order)
            _, _, covariance =
                cold_memory_matrices(bath; mass = physical.mass, omega = physical.omega)
            equilibrium = on_grid(
                (q, p) -> gaussian_wigner(
                    q,
                    p;
                    mean = zeros(2),
                    covariance = covariance[1:2, 1:2],
                ),
                grid,
            )
            max_error(equilibrium, reference)
        end
        @test errors[2] < errors[1] / 2
        @test errors[3] < errors[2] / 2
        @test errors[3] < 1e-4
    end
end

@testset "Cold, strong HEOM: full transient with signed quantum noise" begin
    mass, omega, hbar = 1.0, 1.0, 1.0
    # Negative Drude and Matsubara residues, then an exactly coincident pole.
    for (cutoff, matsubara, negative_index) in ((0.5, 1, 1), (0.8, 2, 2), (2π * 0.1, 2, 0))
        bath = drude_lorentz_bath(; reorganization = 0.8, cutoff, kT = 0.1, hbar, matsubara)
        if negative_index > 0
            @test real(bath.coefficients[negative_index]) < 0
        else
            @test !all(iszero, bath.mixing)
        end
        grid = PhaseSpaceGrid((-8, 8), 56, (-8, 8), 56)
        mean0, covariance0 = [0.7, -0.2], [0.5 0.0; 0.0 0.5]
        initial = on_grid(
            (q, p) -> gaussian_wigner(q, p; mean = mean0, covariance = covariance0),
            grid,
        )
        times = [0.4, 0.8]
        reference = map(times) do t
            mean, covariance =
                cold_memory_reference(t, mean0, covariance0, bath; mass, omega)
            @test det(covariance) >= hbar^2 / 4
            on_grid((q, p) -> gaussian_wigner(q, p; mean, covariance), grid)
        end
        errors = map((2, 4, 6)) do depth
            prob = heom_problem(
                initial,
                (0.0, last(times)),
                grid;
                mass,
                potential = harmonic_potential(; mass, omega),
                bath,
                depth,
                scaled = true,
            )
            sol = solve(prob, Vern7(); abstol = 1e-10, reltol = 1e-10, saveat = times)
            @test last(sol.t) == last(times)
            for (t, U) in zip(sol.t, sol.u)
                W = physical_wigner(U)
                mean, covariance =
                    cold_memory_reference(t, mean0, covariance0, bath; mass, omega)
                @test phase_space_integral(W, grid) ≈ 1 atol = 2e-9
                @test phase_space_mean(W, grid) ≈ mean atol = 2e-7
                @test phase_space_covariance(W, grid) ≈ covariance atol = 2e-6
            end
            maximum(
                max_error(physical_wigner(U), exact) for
                (U, exact) in zip(sol.u, reference)
            )
        end
        @test errors[2] < errors[1] / 20
        @test errors[3] < errors[2] / 5
        @test errors[3] < 2e-7
    end
end
