# Independent references are defined in references.jl.

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
    errors = map((2, 8, 32)) do order
        bath = drude_lorentz_bath(; parameters..., matsubara = order)
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

@testset "Cold, strong HEOM: full transient with signed quantum noise" begin
    mass, omega, hbar = 1.0, 1.0, 1.0
    # A negative Drude residue, then an exactly coincident pole.
    for (cutoff, matsubara, negative_index) in ((0.5, 1, 1), (2π * 0.1, 2, 0))
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
