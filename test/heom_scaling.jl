@testset "Scaled auxiliaries preserve the complete correlated hierarchy" begin
    grid = PhaseSpaceGrid((-6, 6), 20, (-5, 5), 24)
    bath = ExponentialBath(
        [2.4 - 1.3im, -0.7 + 0.2im, 0.0im],
        [0.8, 2.1, 4.3];
        hbar = 0.7,
        diffusion = 0.11,
        counterterm = 0.4,
    )
    for discretization in (Spectral(), FiniteDifference(4))
        options = (;
            mass = 1.3,
            potential = q -> q^2 / 2 + 0.01q^4,
            bath,
            depth = 4,
            discretization,
            moyal_terms = 2,
        )
        unscaled = heom_operator(grid; options...)
        scaled = heom_operator(grid; options..., scaled = true)
        @test !unscaled.scaled
        @test scaled.scaled
        @test hierarchy_indices(scaled) == hierarchy_indices(unscaled)
        U = zeros(size(grid)..., length(unscaled.indices))
        for a in axes(U, 3)
            U[:, :, a] .= on_grid(
                (q, p) -> sin(a) * exp(-q^2 - p^2) * (1 + 0.1a * q - 0.05a * p),
                grid,
            )
        end
        # Independent closed-form normalization, including the zero-residue mode.
        factors = [
            sqrt(
                prod(
                    factorial(n[k]) *
                    (iszero(bath.coefficients[k]) ? 1 : abs(bath.coefficients[k]))^n[k]
                    for k in eachindex(n)
                ),
            ) for n in unscaled.indices
        ]
        expected = U ./ reshape(factors, 1, 1, :)
        transformed = HEOM.rescale_hierarchy(U, unscaled; scaled = true)
        @test transformed ≈ expected rtol = 3e-15
        @test physical_wigner(transformed) == physical_wigner(U)
        @test HEOM.rescale_hierarchy(transformed, scaled; scaled = false) ≈ U rtol = 3e-15
        duplicate = HEOM.rescale_hierarchy(U, unscaled; scaled = false)
        @test duplicate == U
        @test duplicate !== U
        du, ds = similar(U), similar(U)
        heom!(du, U, unscaled, 0.0)
        heom!(ds, transformed, scaled, 0.0)
        @test ds ≈ du ./ reshape(factors, 1, 1, :) rtol = 2e-14
        prob = heom_problem(transformed, (0.0, 0.1), scaled)
        @test prob.u0 == transformed
        @test prob.u0 !== transformed
    end
end

@testset "Deep scaling and representable conversion" begin
    grid = PhaseSpaceGrid((-4, 4), 8, (-4, 4), 8)
    # Both factorial(180) and abs(c) overflow Float64, but all scaled local
    # transition factors and the hierarchy remain representable.
    bath = ExponentialBath([complex(1.3e308, -1.3e308)], [0.5])
    op = heom_operator(
        grid;
        mass = 1,
        potential = q -> zero(q),
        bath,
        depth = 180,
        scaled = true,
    )
    @test all(isfinite, op.lowering)
    @test all(isfinite, op.raising_derivative)
    @test all(isfinite, op.raising_coordinate)
    @test op.log_scales[end] > log(floatmax(Float64))
    U = zeros(size(grid)..., length(op.indices))
    U[:, :, 1] .= 0.01
    dU = similar(U)
    heom!(dU, U, op, 0.0)
    @test all(isfinite, dU)
    @test HEOM.rescale_hierarchy(U, op; scaled = false) == U
    # The scale itself can overflow even when the converted value does not.
    U[1, 1, 3] = 1e-308
    converted = HEOM.rescale_hierarchy(U, op; scaled = false)
    @test converted[1, 1, 3] ≈ 2.6 rtol = 3e-13
    U[1, 1, end] = 1
    @test_throws ArgumentError HEOM.rescale_hierarchy(U, op; scaled = false)

    tiny = ExponentialBath([1e-308], [0.5])
    tiny_op =
        heom_operator(grid; mass = 1, potential = q -> zero(q), bath = tiny, depth = 3)
    tiny_state = zeros(size(grid)..., length(tiny_op.indices))
    tiny_state[1, 1, 4] = 1e-300
    tiny_scaled = HEOM.rescale_hierarchy(tiny_state, tiny_op; scaled = true)
    @test tiny_scaled[1, 1, 4] ≈ 1e162 / sqrt(6) rtol = 3e-13
    tiny_state[1, 1, 4] = 1
    @test_throws ArgumentError HEOM.rescale_hierarchy(tiny_state, tiny_op; scaled = true)
    @test_throws DimensionMismatch HEOM.rescale_hierarchy(zeros(8, 8), op; scaled = true)
    @test_throws ArgumentError HEOM.rescale_hierarchy(complex.(U), op; scaled = true)
    @test_throws ArgumentError HEOM.rescale_hierarchy(fill(Inf, size(U)), op; scaled = true)
    @test_throws ArgumentError HEOM.rescale_hierarchy(
        fill(big"1e400", size(U)),
        op;
        scaled = true,
    )
end

@testset "Hierarchy resource estimates and allocation guard" begin
    @test HEOM.hierarchy_size(0, 1000) == 1
    @test HEOM.hierarchy_size(1000, 0) == 1
    @test HEOM.hierarchy_size(7, 12) == 50388
    @test HEOM.hierarchy_size(40, 40) == binomial(big(80), 40)
    @test HEOM.hierarchy_size(40, 40) > typemax(Int)
    @test_throws ArgumentError HEOM.hierarchy_size(-1, 2)
    @test_throws ArgumentError HEOM.hierarchy_size(1, -2)
    grid = PhaseSpaceGrid((-4, 4), 8, (-4, 4), 8)
    bath = ExponentialBath([1, 2], [1, 2])
    options = (; mass = 1, potential = q -> q^2 / 2, bath)
    @test length(heom_operator(grid; options..., depth = 3, max_ados = 10).indices) == 10
    @test_throws ArgumentError heom_operator(grid; options..., depth = 3, max_ados = 9)
    @test_throws ArgumentError heom_operator(grid; options..., depth = 3, max_ados = 0)
    @test_throws ArgumentError heom_operator(
        grid;
        options...,
        depth = big(typemax(Int)) + 1,
    )
    @test_throws ArgumentError heom_operator(grid; options..., depth = typemax(Int))
    # State byte counts are checked before any hierarchy enumeration or FFT plans.
    @test_throws ArgumentError heom_operator(grid; options..., depth = 1000000000)
end

@testset "Cold strong-coupling evolution is invariant under scaling" begin
    grid = PhaseSpaceGrid((-8, 8), 32, (-8, 8), 32)
    bath =
        drude_lorentz_bath(; reorganization = 1.5, cutoff = 0.7, kT = 0.15, matsubara = 3)
    options = (; mass = 1, potential = q -> q^2 / 2, bath, depth = 4)
    W0 = on_grid(
        (q, p) -> coherent_wigner(q, p; q0 = 0.3, p0 = -0.2, mass = 1, omega = 1),
        grid,
    )
    solutions = map((false, true)) do scaled
        prob = heom_problem(W0, (0.0, 0.25), grid; options..., scaled)
        solve(prob, Vern7(); abstol = 2e-11, reltol = 2e-10, save_everystep = false)
    end
    original, transformed = solutions
    @test physical_wigner(original) ≈ physical_wigner(transformed) rtol = 1e-8
    @test original.u[end] ≈
          HEOM.rescale_hierarchy(transformed.u[end], transformed.prob.p; scaled = false) rtol =
        3e-8
end
