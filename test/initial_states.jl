@testset "Wavefunction and density-matrix Wigner transforms" begin
    mass, omega, hbar = 1.4, 0.9, 0.7
    q0, p0 = 0.65, 1.15
    a = mass * omega / hbar
    psi(q) = (a / π)^(1 / 4) * exp(-a * (q - q0)^2 / 2 + im * p0 * q / hbar)
    rho(x, y) = psi(x) * conj(psi(y))

    # Neither axis needs to be centred on zero, and the momentum count can be odd.
    for np in (95, 96)
        grid = PhaseSpaceGrid((-7.3, 8.7), 128, (-6.35, 7.85), np)
        reference =
            on_grid((q, p) -> coherent_wigner(q, p; q0, p0, mass, omega, hbar), grid)
        W = wavefunction_wigner(psi, grid; hbar)
        @test size(W) == size(grid)
        @test eltype(W) <: Real
        @test max_error(W, reference) < 1e-10
        @test phase_space_integral(W, grid) ≈ 1 atol = 1e-10
        @test purity(W, grid; hbar) ≈ 1 atol = 1e-10
        @test phase_space_mean(W, grid) ≈ [q0, p0] atol = 1e-10
        @test position_density(W, grid) ≈ abs2.(psi.(grid.q)) atol = 1e-12
        @test density_matrix_wigner(rho, grid; hbar) ≈ W atol = 1e-13
        @test wavefunction_wigner(q -> cis(0.43) * psi(q), grid; hbar) ≈ W atol = 1e-13

        # Transforming is linear in the density kernel and does not renormalize it.
        @test wavefunction_wigner(q -> 2psi(q), grid; hbar) ≈ 4W atol = 1e-12
        @test density_matrix_wigner((x, y) -> 0.4rho(x, y), grid; hbar) ≈ 0.4W atol = 1e-12
    end

    grid = PhaseSpaceGrid((-8, 8), 128, (-8, 8), 128)
    ground(q) = (a / π)^(1 / 4) * exp(-a * q^2 / 2)
    excited(q) = sqrt(2a) * q * ground(q)
    W1 = wavefunction_wigner(excited, grid; hbar)
    reference1 = on_grid((q, p) -> fock_wigner(1, q, p; mass, omega, hbar), grid)
    @test max_error(W1, reference1) < 1e-10
    @test W1[65, 65] ≈ -1 / (π * hbar) atol = 1e-10

    weight, trace = 0.3, 1.7
    mixture(x, y) =
        trace * (weight * ground(x) * ground(y) + (1 - weight) * excited(x) * excited(y))
    mixed = density_matrix_wigner(mixture, grid; hbar)
    @test phase_space_integral(mixed, grid) ≈ trace atol = 1e-10
    @test purity(mixed, grid; hbar) ≈ trace^2 * (weight^2 + (1 - weight)^2) atol = 1e-10
    reference0 = on_grid((q, p) -> fock_wigner(0, q, p; mass, omega, hbar), grid)
    @test mixed ≈ trace * (weight * reference0 + (1 - weight) * reference1) atol = 1e-10

    # Sampled inputs converge under position refinement; complex conjugation and
    # the position-density quadrature convention must match the callable API.
    errors = map((64, 128)) do nq
        sampled_grid = PhaseSpaceGrid((-7.3, 8.7), nq, (-6.35, 7.85), 96)
        samples = psi.(sampled_grid.q)
        W = wavefunction_wigner(samples, sampled_grid; hbar)
        exact = wavefunction_wigner(psi, sampled_grid; hbar)
        @test density_matrix_wigner(samples * samples', sampled_grid; hbar) ≈ W atol =
            1e-12
        @test position_density(W, sampled_grid) ≈ abs2.(samples) atol = 1e-12
        @test phase_space_integral(W, sampled_grid) ≈ sum(abs2, samples) * sampled_grid.dq
        max_error(W, exact)
    end
    @test errors[2] < errors[1] / 8
    @test errors[2] < 3e-4

    # With fewer than four position samples the interpolant lowers its degree.
    # Zero extension still preserves the grid values and their quadrature norm.
    for nq in (2, 3)
        small = PhaseSpaceGrid((-1, 2), nq, (-4, 4), 16)
        samples = ComplexF64.(1:nq) .+ im
        W = wavefunction_wigner(samples, small)
        @test position_density(W, small) ≈ abs2.(samples) atol = 1e-13
        @test density_matrix_wigner(samples * samples', small) ≈ W atol = 1e-13
    end
    @test all(iszero, wavefunction_wigner(zeros(128), grid))
    @test all(iszero, density_matrix_wigner(zeros(128, 128), grid))
end

@testset "Numerical eigenstate preparation" begin
    mass, omega, hbar = 1.3, 0.85, 0.7
    grid = PhaseSpaceGrid((-9, 9), 128, (-8, 8), 128)
    potential = harmonic_potential(; mass, omega)
    states = eigenstates(grid; mass, potential, hbar, nstates = 4)
    @test size(states.wavefunctions) == (length(grid.q), 4)
    @test states.energies ≈ hbar * omega .* (collect(0:3) .+ 1 / 2) atol = 1e-10
    @test states.wavefunctions' * states.wavefunctions * grid.dq ≈ I(4) atol = 1e-12

    for n in 0:2
        W = eigenstate_wigner(n, grid; mass, potential, hbar)
        reference = on_grid((q, p) -> fock_wigner(n, q, p; mass, omega, hbar), grid)
        @test max_error(W, reference) < 5e-4
        @test phase_space_integral(W, grid) ≈ 1 atol = 1e-12
        @test purity(W, grid; hbar) ≈ 1 atol = 2e-3
        @test W ≈ wavefunction_wigner(states.wavefunctions[:, n+1], grid; hbar) atol = 1e-12
    end

    # A periodic free particle checks the finite-box spectrum, including the
    # unpaired even-grid Nyquist mode, which has nonzero kinetic energy.
    for nq in (15, 16)
        free_grid = PhaseSpaceGrid((-2.3, 4.7), nq, (-5, 5), 24)
        free = eigenstates(free_grid; mass, potential = q -> 0.0, hbar)
        modes = (-fld(nq, 2)):(cld(nq, 2)-1)
        reference = sort(hbar^2 .* (2π .* modes ./ 7) .^ 2 ./ (2mass))
        @test free.energies ≈ reference atol = 1e-11
        @test free.wavefunctions' * free.wavefunctions * free_grid.dq ≈ I(nq) atol = 1e-12
    end

    # Quartic wells exercise numerical preparation without a harmonic ansatz.
    double_well(q) = (q^2 - 2.25)^2 / 8
    molecular_grid = PhaseSpaceGrid((-7, 7), 128, (-7, 7), 128)
    molecular =
        eigenstates(molecular_grid; mass, potential = double_well, hbar, nstates = 4)
    @test issorted(molecular.energies)
    @test all(diff(molecular.energies) .> 0)
    @test molecular.wavefunctions' * molecular.wavefunctions * molecular_grid.dq ≈ I(4) atol =
        1e-12
    # An independent sixth-order derivative verifies the Schrödinger residual.
    D2 = HEOM.periodic_difference_matrix(2, 6, length(molecular_grid.q), molecular_grid.dq)
    for n in 1:4
        psi = molecular.wavefunctions[:, n]
        residual =
            -hbar^2 / (2mass) * D2 * psi + double_well.(molecular_grid.q) .* psi -
            molecular.energies[n] * psi
        @test norm(residual) * sqrt(molecular_grid.dq) < 3e-5
    end
    W = eigenstate_wigner(0, molecular_grid; mass, potential = double_well, hbar)
    @test phase_space_integral(W, molecular_grid) ≈ 1 atol = 1e-12
    op = wigner_moyal_operator(molecular_grid; mass, potential = double_well, hbar)
    @test maximum(abs, rhs(op, W)) < 5e-3
end

@testset "Numerical thermal-state preparation" begin
    mass, omega, hbar = 1.3, 0.8, 0.7
    potential = harmonic_potential(; mass, omega)
    grid = PhaseSpaceGrid((-10, 10), 192, (-10, 10), 128)
    H(q, p) = p^2 / (2mass) + potential(q)
    for kT in (0.12, 1.7)
        W = thermal_wigner(grid; mass, potential, hbar, kT)
        a = tanh(hbar * omega / (2kT))
        reference =
            on_grid((q, p) -> a / (π * hbar) * exp(-2a * H(q, p) / (hbar * omega)), grid)
        @test max_error(W, reference) < 5e-4
        @test phase_space_integral(W, grid) ≈ 1 atol = 1e-11
        @test purity(W, grid; hbar) ≈ a atol = 2e-3
        @test expectation(H, W, grid) ≈ hbar * omega / (2a) atol = 1e-2
    end

    ground = eigenstate_wigner(0, grid; mass, potential, hbar)
    cold = thermal_wigner(grid; mass, potential, hbar, kT = 1e-4)
    @test cold ≈ ground atol = 1e-12
    @test thermal_wigner(grid; mass, potential, hbar, kT = 3.0, nstates = 1) ≈ ground atol =
        1e-12

    # Truncated Gibbs weights are normalized over the requested states, and
    # subtracting E₀ makes them independent of a large potential-energy offset.
    kT = 0.6
    states = eigenstates(grid; mass, potential, hbar, nstates = 4)
    weights = exp.(-(states.energies .- states.energies[1]) / kT)
    weights ./= sum(weights)
    reference = sum(
        weights[n] * wavefunction_wigner(states.wavefunctions[:, n], grid; hbar) for
        n in 1:4
    )
    W = thermal_wigner(grid; mass, potential, hbar, kT, nstates = 4)
    @test W ≈ reference atol = 1e-12
    shifted = thermal_wigner(
        grid;
        mass,
        potential = q -> potential(q) - 1e4,
        hbar,
        kT,
        nstates = 4,
    )
    @test shifted ≈ W atol = 1e-10
end

@testset "Initial-state input validation" begin
    grid = PhaseSpaceGrid((-6, 6), 32, (-6, 6), 32)
    psi(q) = exp(-q^2 / 2)
    rho(x, y) = psi(x) * psi(y)
    potential(q) = q^2 / 2
    for bad in (0.0, -1.0, Inf, -Inf, NaN)
        @test_throws ArgumentError wavefunction_wigner(psi, grid; hbar = bad)
        @test_throws ArgumentError density_matrix_wigner(rho, grid; hbar = bad)
        @test_throws ArgumentError eigenstates(grid; mass = bad, potential)
        @test_throws ArgumentError eigenstates(grid; mass = 1, potential, hbar = bad)
        @test_throws ArgumentError thermal_wigner(grid; mass = 1, potential, kT = bad)
    end
    @test_throws DimensionMismatch wavefunction_wigner(ones(31), grid)
    @test_throws DimensionMismatch density_matrix_wigner(ones(32, 31), grid)
    @test_throws ArgumentError wavefunction_wigner(fill(NaN, 32), grid)
    @test_throws ArgumentError wavefunction_wigner(q -> Inf, grid)
    @test_throws ArgumentError wavefunction_wigner(q -> "invalid", grid)
    @test_throws ArgumentError wavefunction_wigner(q -> big"1e400", grid)
    @test_throws ArgumentError wavefunction_wigner(psi, grid; hbar = big"1e-1000")
    @test_throws ArgumentError density_matrix_wigner(fill(NaN, 32, 32), grid)
    @test_throws ArgumentError density_matrix_wigner((x, y) -> NaN, grid)
    @test_throws ArgumentError density_matrix_wigner((x, y) -> im * rho(x, y), grid)
    nonhermitian = Matrix{ComplexF64}(I, 32, 32)
    nonhermitian[1, 2] = im
    @test_throws ArgumentError density_matrix_wigner(nonhermitian, grid)
    tiny_momentum = PhaseSpaceGrid((-1, 1), 4, (-1e-310, 1e-310), 4)
    @test_throws ArgumentError wavefunction_wigner(psi, tiny_momentum)
    @test_throws ArgumentError density_matrix_wigner((x, y) -> 1e308, grid)
    for badpotential in (q -> NaN, q -> Inf, q -> 1im)
        @test_throws ArgumentError eigenstates(grid; mass = 1, potential = badpotential)
    end
    @test_throws ArgumentError eigenstates(grid; mass = 1e-320, potential)
    @test_throws ArgumentError eigenstates(grid; mass = big"1e-1000", potential)
    @test_throws ArgumentError eigenstates(grid; mass = 1, potential, hbar = big"1e-1000")
    @test_throws ArgumentError thermal_wigner(grid; mass = 1, potential, kT = big"1e-1000")
    for nstates in (-1, 0, 33)
        @test_throws ArgumentError eigenstates(grid; mass = 1, potential, nstates)
        @test_throws ArgumentError thermal_wigner(
            grid;
            mass = 1,
            potential,
            kT = 1,
            nstates,
        )
    end
    for n in (-1, 32)
        @test_throws ArgumentError eigenstate_wigner(n, grid; mass = 1, potential)
    end
end
