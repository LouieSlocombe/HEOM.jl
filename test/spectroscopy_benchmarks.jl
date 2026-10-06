# Physical references for driven dynamics and spectroscopy. The harmonic flow and
# the eigenstate Kubo sum are independent of the phase-space propagation.
@testset "Driven harmonic oscillator follows the exact forced Gaussian" begin
    mass, omega, hbar = 1.4, 1.2, 0.8
    amplitude, frequency = 0.3, 0.75
    q0, p0 = 0.3, -0.2
    grid = PhaseSpaceGrid((-7, 7), 56, (-7, 7), 56)
    V0 = harmonic_potential(; mass, omega)
    field(t) = amplitude * cos(frequency * t)
    W0 = on_grid((q, p) -> coherent_wigner(q, p; q0, p0, mass, omega, hbar), grid)
    times = 0.0:0.3:2.7
    for potential in (
        DrivenPotential(V0, field, identity),
        TimeDependentPotential((q, t) -> V0(q) - field(t) * q),
    )
        prob = wigner_moyal_problem(W0, (0.0, last(times)), grid; mass, hbar, potential)
        sol = solve(prob, Vern7(); abstol = 1e-10, reltol = 1e-10, saveat = times)
        @test last(sol.t) ≈ last(times)
        for (t, W) in zip(sol.t, sol.u)
            # m*q'' + m*omega^2*q = amplitude*cos(frequency*t).
            qcentre =
                q0 * cos(omega * t) +
                p0 / (mass * omega) * sin(omega * t) +
                amplitude / (mass * (omega^2 - frequency^2)) *
                (cos(frequency * t) - cos(omega * t))
            pcentre =
                p0 * cos(omega * t) - mass * omega * q0 * sin(omega * t) +
                amplitude / (omega^2 - frequency^2) *
                (omega * sin(omega * t) - frequency * sin(frequency * t))
            exact = on_grid(grid) do q, p
                coherent_wigner(q, p; q0 = qcentre, p0 = pcentre, mass, omega, hbar)
            end
            @test maximum(abs, W - exact) < 2e-7
            @test phase_space_mean(W, grid) ≈ [qcentre, pcentre] atol = 2e-8
            @test phase_space_integral(W, grid) ≈ 1 atol = 2e-10
        end
    end
end

@testset "Harmonic response and absorption have the correct sign and scale" begin
    mass, omega, hbar = 1.3, 1.1, 0.8
    grid = PhaseSpaceGrid((-6, 6), 48, (-6, 6), 48)
    V = harmonic_potential(; mass, omega)
    W0 = on_grid((q, p) -> coherent_wigner(q, p; mass, omega, hbar), grid)
    op = wigner_moyal_operator(grid; mass, hbar, potential = V)
    result = linear_response(
        W0,
        (0.0, 30.0),
        op,
        Vern7();
        saveat = 0.025,
        abstol = 1e-10,
        reltol = 1e-10,
    )
    exact = sin.(omega .* result.times) ./ (mass * omega)
    @test maximum(abs, result.response - exact) < 2e-7
    @test maximum(abs(phase_space_integral(W, grid)) for W in result.solution.u) < 1e-10

    broadening = 0.3
    frequencies = 0.0:0.01:2.0
    spectrum = absorption_spectrum(result; frequencies, broadening)
    # Infinite-time retarded susceptibility with exponential damping exp(-eta*t).
    analytic = @. 1 / (mass * (omega^2 - (frequencies + im * broadening)^2))
    @test maximum(abs, spectrum.susceptibility - analytic) < 5e-4
    @test maximum(abs, spectrum.intensity - frequencies .* imag.(analytic)) < 5e-4
    peak = spectrum.frequencies[argmax(spectrum.intensity)]
    @test abs(peak - sqrt(omega^2 + broadening^2)) < 0.015
    @test iszero(first(spectrum.intensity))
end

@testset "Quartic response converges to the eigenstate Kubo sum" begin
    mass, hbar = 1.0, 1.0
    potential(q) = q^2 / 2 + 0.06q^4
    errors = Float64[]
    for points in (40, 56)
        grid = PhaseSpaceGrid((-6, 6), points, (-6, 6), 48)
        states = eigenstates(grid; mass, hbar, potential, nstates = 12)
        ground = states.wavefunctions[:, 1]
        W0 = wavefunction_wigner(ground, grid; hbar)
        @test phase_space_integral(W0, grid) ≈ 1 atol = 1e-12
        op = wigner_moyal_operator(grid; mass, hbar, potential)
        result = linear_response(
            W0,
            (0.0, 8.0),
            op,
            Vern7();
            saveat = 0.05,
            abstol = 1e-10,
            reltol = 1e-10,
        )
        gaps = (states.energies .- first(states.energies)) ./ hbar
        dipoles = states.wavefunctions' * (grid.q .* ground) .* grid.dq
        reference =
            [2 / hbar * sum(abs2.(dipoles) .* sin.(gaps .* t)) for t in result.times]
        push!(errors, maximum(abs, result.response - reference))
        @test gaps[2] > 1.1 # Resolved anharmonic shift from the harmonic frequency 1.
        @test maximum(abs(phase_space_integral(W, grid)) for W in result.solution.u) < 1e-10
    end
    @test errors[2] < errors[1] / 2
    @test errors[2] < 2e-4
end
