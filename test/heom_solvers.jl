@testset "Exact HEOM sparse generator and Jacobian products" begin
    grid = PhaseSpaceGrid((-4, 4), 8, (-5, 5), 10)
    baths = (
        ExponentialBath(ComplexF64[], Float64[]),
        ExponentialBath(
            [0.4 - 0.3im, -0.06],
            [0.7, 9.0];
            diffusion = 0.08,
            counterterm = 0.2,
        ),
        ExponentialBath(
            [0.4 - 0.3im, -0.06],
            [0.7, 0.7];
            weights = [1.0, 0.0],
            mixing = [0 1; 0 0],
            diffusion = 0.08,
            counterterm = 0.2,
        ),
    )
    for bath in baths, scaled in (false, true)
        op = heom_operator(
            grid;
            mass = 1.2,
            potential = q -> q^4 / 10 + q^2 / 2,
            bath,
            depth = 2,
            scaled,
            discretization = FiniteDifference(4),
            moyal_terms = 2,
        )
        dimensions = (size(grid)..., length(hierarchy_indices(op)))
        U = reshape(sin.(0.3 .* (1:prod(dimensions))), dimensions)
        dU = hierarchy_rhs(op, U)
        L = sparse(op)
        @test L * vec(U) ≈ vec(dU) rtol = 1e-13
        prob = heom_problem(U, (0, 1), op; jacobian = :sparse)
        @test prob.f.jac_prototype == L
        J = zero(L)
        prob.f.jac(J, U, op, 0.0)
        @test J == L
        @test_throws ArgumentError prob.f.jac(J, U, deepcopy(op), 0.0)
        product = zeros(length(U))
        prob.f.jvp(product, vec(U), U, op, 0.0)
        @test product ≈ vec(dU)
        prob.f.tgrad(product, U, op, 0.0)
        @test iszero(product)
        @test_throws ArgumentError heom_problem(U, (0, 1), op; jacobian = :unknown)
    end
    op = heom_operator(grid; mass = 1, potential = q -> q^2 / 2, bath = baths[2], depth = 1)
    @test_throws ArgumentError sparse(op)
    @test_throws ArgumentError heom_problem(
        zeros(size(grid)),
        (0, 1),
        op;
        jacobian = :sparse,
    )
    invalid = heom_operator(
        grid;
        mass = 1,
        potential = q -> q^2 / 2,
        bath = baths[2],
        depth = 1,
        discretization = FiniteDifference(4),
        moyal_terms = 1,
    )
    invalid.damping[end] = Inf
    @test_throws ArgumentError sparse(invalid)
end

@testset "Stiff HEOM integration with exact derivatives" begin
    grid = PhaseSpaceGrid((-5, 5), 16, (-5, 5), 16)
    bath = ExponentialBath(
        [0.2 - 0.1im, 0.04],
        [1.0, 40.0];
        diffusion = 0.1,
        counterterm = 0.1,
    )
    W0 = on_grid((q, p) -> coherent_wigner(q, p; mass = 1, omega = 1, q0 = 0.4), grid)
    for (discretization, jacobian) in
        ((Spectral(), :matrixfree), (FiniteDifference(4), :sparse))
        prob = heom_problem(
            W0,
            (0.0, 0.2),
            grid;
            mass = 1,
            potential = q -> q^2 / 2,
            bath,
            depth = 2,
            scaled = true,
            discretization,
            moyal_terms = 1,
            jacobian,
        )
        # Krylov uses the supplied exact JVP, not finite differences of the FFT.
        alg =
            jacobian === :sparse ? Rodas5P() :
            Rodas5P(
                autodiff = AutoFiniteDiff(),
                linsolve = KrylovJL_GMRES(),
                concrete_jac = false,
            )
        sol = solve(prob, alg; abstol = 1e-9, reltol = 1e-9, save_everystep = false)
        reference =
            solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12, save_everystep = false)
        @test last(sol.t) == last(prob.tspan)
        @test max_error(sol.u[end], reference.u[end]) < 1e-9
        @test phase_space_integral(physical_wigner(sol), grid) ≈
              phase_space_integral(W0, grid) atol = 1e-12
    end
end
