# Initial-state preparation for tunnelling and a molecular potential.
# Requires only HEOM. Run: julia --project=. examples/initial_states.jl
using HEOM, Printf

function double_well_preparation()
    mass, hbar = 1.0, 1.0
    potential(q) = (q^2 - 4)^2 / 8
    grid = PhaseSpaceGrid((-6.0, 6.0), 192, (-8.0, 8.0), 128)
    states = eigenstates(grid; mass, potential, hbar, nstates = 2)
    psi0, psi1 = states.wavefunctions[:, 1], states.wavefunctions[:, 2]

    # Eigenvector phases are arbitrary. Align the transition matrix element
    # <0|q|1> to be positive so (psi0 - psi1)/sqrt(2) occupies the left well.
    transition = sum(conj.(psi0) .* grid.q .* psi1) * grid.dq
    psi1 = psi1 * (conj(transition) / abs(transition))
    splitting = states.energies[2] - states.energies[1]
    transfer_time = π * hbar / splitting
    @printf("Double-well energies: %.8f, %.8f\n", states.energies...)
    @printf("Left-to-right transfer time: %.6f\n", transfer_time)
    println("Time, probability in the right well:")

    # Evolve the two stationary components by their relative phase. This is an
    # exact solution in the numerical two-state subspace and needs no ODE solver.
    W_left = wavefunction_wigner((psi0 - psi1) / sqrt(2), grid; hbar)
    for time in range(0.0, transfer_time; length = 5)
        psi = (psi0 - exp(-im * splitting * time / hbar) * psi1) / sqrt(2)
        W = wavefunction_wigner(psi, grid; hbar)
        @assert abs(phase_space_integral(W, grid) - 1) < 1e-10
        @printf("  %.6f, %.6f\n", time, probability(W, grid; q = (0.0, 6.0)))
    end
    @assert probability(W_left, grid; q = (0.0, 6.0)) < 0.02

    # Kernel samples, not basis coefficients: multiply their diagonal sum by dq
    # to obtain the trace. A classical mixture has no superposition interference.
    rho = 0.7 * (psi0 * psi0') + 0.3 * (psi1 * psi1')
    W_mixture = density_matrix_wigner(rho, grid; hbar)
    mixture_reference =
        0.7 * wavefunction_wigner(psi0, grid; hbar) +
        0.3 * wavefunction_wigner(psi1, grid; hbar)
    @assert maximum(abs, W_mixture - mixture_reference) < 1e-10

    # This is the isolated double well's Gibbs state. It can seed a factorised
    # HEOM preparation; use equilibrate to prepare correlated bath auxiliaries.
    W_thermal = thermal_wigner(grid; mass, potential, hbar, kT = 0.3)
    @printf("Thermal-state norm: %.12f\n", phase_space_integral(W_thermal, grid))
    @printf("Thermal-state purity: %.8f\n", purity(W_thermal, grid; hbar))
    return (; grid, potential, W_left, W_mixture, W_thermal)
end

function morse_preparation()
    mass, hbar = 1.0, 1.0
    dissociation, range_parameter = 8.0, 0.4
    potential(q) = dissociation * (1 - exp(-range_parameter * q))^2
    grid = PhaseSpaceGrid((-4.0, 12.0), 256, (-8.0, 8.0), 128)
    W_ground = eigenstate_wigner(0, grid; mass, potential, hbar)
    ground_energy = eigenstates(grid; mass, potential, hbar, nstates = 1).energies[1]
    exact_energy =
        hbar * range_parameter * sqrt(2dissociation / mass) / 2 -
        (hbar * range_parameter)^2 / (8mass)
    @printf("Morse ground energy: %.8f (analytic %.8f)\n", ground_energy, exact_energy)
    @assert abs(ground_energy - exact_energy) < 1e-6
    @assert abs(phase_space_integral(W_ground, grid) - 1) < 1e-10

    # Callable input can evaluate shifted positions directly. The sampled
    # equivalent uses local interpolation and converges as dq is reduced.
    width, q0, p0 = 0.8, 0.8, 0.6
    psi(q) = exp(-(q - q0)^2 / (2width^2) + im * p0 * q / hbar) / sqrt(sqrt(π) * width)
    W_packet = wavefunction_wigner(psi, grid; hbar)
    W_density = density_matrix_wigner((q, qp) -> psi(q) * conj(psi(qp)), grid; hbar)
    @assert maximum(abs, W_packet - W_density) < 1e-10
    W_samples = wavefunction_wigner(psi.(grid.q), grid; hbar)
    @printf(
        "Sampled/callable packet difference: %.3e\n",
        maximum(abs, W_samples - W_packet)
    )

    # Check convergence in both boxes and spacings before using quantitative
    # results. Morse has a continuum above dissociation: a thermal_wigner here
    # would be a finite-box Gibbs state, not infinite-line thermal equilibrium.
    return (; grid, potential, W_ground, W_packet)
end

double_well = double_well_preparation()
morse = morse_preparation()
