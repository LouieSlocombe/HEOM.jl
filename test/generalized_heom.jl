@testset "Generalized real correlation basis" begin
    c = [0.5 - 0.2im, -0.7 + 0.1im]
    rates = [1.2, 1.2]
    weights = [1.0, 0.0]
    mixing = [0.0 1.0; 0.0 0.0]
    bath = ExponentialBath(c, rates; weights, mixing)
    @test bath.weights == weights
    @test bath.mixing == mixing
    # A Jordan block is a polynomial exponential, including at exact degeneracy.
    for t in (0.0, 0.1, 0.7, 2.0)
        correlation =
            transpose(bath.weights) * exp(-Matrix(Diagonal(rates) + bath.mixing) * t) * c
        @test correlation ≈ (c[1] - t * c[2]) * exp(-rates[1] * t) rtol = 3e-15
    end
    weights[1] = 3
    mixing[1, 2] = 5
    @test bath.weights == [1, 0]
    @test bath.mixing == [0 1; 0 0]
    independent = ExponentialBath(c, rates)
    @test independent.weights == ones(2)
    @test iszero(independent.mixing)
    copied = ExponentialBath(c, rates; mixing = bath.mixing)
    copied.mixing[1, 2] = 9
    @test bath.mixing[1, 2] == 1
    large = ExponentialBath(ones(10000), ones(10000))
    @test Base.summarysize(large.mixing) < 1000000
    # A real rotation block permits damped oscillatory correlations without
    # introducing complex Wigner auxiliaries.
    oscillatory =
        ExponentialBath([1.0, 0.0], [0.3, 0.3]; weights = [1, 0], mixing = [0 -2; 2 0])
    for t in (0.0, 0.2, 1.1)
        correlation =
            transpose(oscillatory.weights) *
            exp(-Matrix(Diagonal(oscillatory.rates) + oscillatory.mixing) * t) *
            oscillatory.coefficients
        @test correlation ≈ exp(-0.3t) * cos(2t)
    end
    # Adding a long independent tail must not turn oscillator stability
    # validation into a dense eigensystem of all Matsubara modes.
    modes = 2002
    large_coefficients, large_rates = ones(modes), ones(modes)
    large_mixing = sparse([1, 2], [2, 1], [-2.0, 2.0], modes, modes)
    ExponentialBath(large_coefficients, large_rates; mixing = large_mixing)
    allocated =
        @allocated ExponentialBath(large_coefficients, large_rates; mixing = large_mixing)
    @test allocated < 2_000_000
    @test_throws DimensionMismatch ExponentialBath(c, rates; weights = [1])
    @test_throws DimensionMismatch ExponentialBath(c, rates; mixing = zeros(2, 1))
    @test_throws ArgumentError ExponentialBath(c, rates; weights = [1 + im, 0])
    @test_throws ArgumentError ExponentialBath(c, rates; mixing = [0 1im; 0 0])
    @test_throws ArgumentError ExponentialBath(c, rates; weights = [Inf, 0])
    @test_throws ArgumentError ExponentialBath(c, rates; mixing = [0 NaN; 0 0])
    @test_throws ArgumentError ExponentialBath(c, rates; mixing = [1 0; 0 0])
    @test_throws ArgumentError ExponentialBath(c, rates; mixing = [0 2; 2 0])
end

@testset "Changing the bath basis preserves physical propagation" begin
    grid = PhaseSpaceGrid((-6, 6), 24, (-6, 6), 24)
    rates = [0.7, 1.5]
    coefficients = [0.4 - 0.2im, 0.8]
    diagonal = ExponentialBath(coefficients, rates)
    # T = [1 1; 0 (ν-γ)] changes diagonal exponentials to a divided-difference
    # basis. The transformed generator has mixing[1,2]=1 and weights=[1,0].
    generalized = ExponentialBath(
        [sum(coefficients), (rates[2] - rates[1]) * coefficients[2]],
        rates;
        weights = [1, 0],
        mixing = [0 1; 0 0],
    )
    W0 = on_grid((q, p) -> coherent_wigner(q, p; mass = 1, omega = 1, q0 = 0.3), grid)
    for scaled in (false, true)
        final_states = map((diagonal, generalized)) do bath
            prob = heom_problem(
                W0,
                (0.0, 0.4),
                grid;
                mass = 1,
                potential = q -> q^2 / 2,
                bath,
                depth = 3,
                scaled,
            )
            sol = solve(
                prob,
                Vern7();
                reltol = 1e-10,
                abstol = 1e-12,
                save_everystep = false,
            )
            copy(physical_wigner(sol))
        end
        @test final_states[1] ≈ final_states[2] rtol = 3e-10
    end
end

@testset "Analytical same-tier generalized hierarchy" begin
    grid = PhaseSpaceGrid((-4, 4), 8, (-4, 4), 8)
    # Constant phase-space members isolate exact same-tier damping and transfer.
    bath = ExponentialBath(
        [0.4, -0.8],
        [1.3, 0.7];
        weights = [1.0, 0.2],
        mixing = [0 -0.4; 0.3 0],
    )
    op = heom_operator(grid; mass = 1, potential = q -> zero(q), bath, depth = 3)
    values = [sin(0.7a) for a in eachindex(op.indices)]
    U = repeat(reshape(values, 1, 1, :), size(grid)...)
    dU = similar(U)
    heom!(dU, U, op, 0.0)
    lookup = Dict(Tuple(n) => a for (a, n) in enumerate(op.indices))
    for (a, n) in enumerate(op.indices)
        expected = -dot(n, bath.rates) * values[a]
        for k in eachindex(n), j in eachindex(n)
            iszero(n[k]) && continue
            shifted = copy(n)
            shifted[k] -= 1
            shifted[j] += 1
            expected -= n[k] * bath.mixing[k, j] * values[lookup[Tuple(shifted)]]
        end
        @test all(x -> isapprox(x, expected; atol = 3e-15), @view(dU[:, :, a]))
    end
    @test iszero(op.transfer[1, :])
    @test iszero(op.transfer[:, 1])
    bath.weights[1] = 12
    bath.mixing[1, 2] = 17
    @test op.bath.weights == [1.0, 0.2]
    @test op.bath.mixing == [0 -0.4; 0.3 0]
end

@testset "Generalized hierarchy scaling is an exact similarity" begin
    grid = PhaseSpaceGrid((-5, 5), 16, (-5, 5), 16)
    for mixing in ([0.0 1.0; 0.0 0.0], [0.0 -0.8; 0.4 0.0])
        bath =
            ExponentialBath([0.5 - 0.2im, -0.8], [1.1, 1.1]; weights = [1.0, -0.3], mixing)
        options = (; mass = 1, potential = q -> q^2 / 2, bath, depth = 4)
        unscaled = heom_operator(grid; options...)
        scaled = heom_operator(grid; options..., scaled = true)
        U = zeros(size(grid)..., length(unscaled.indices))
        for a in axes(U, 3)
            U[:, :, a] .= on_grid((q, p) -> cos(a) * exp(-q^2 - p^2) * (1 + 0.1a * p), grid)
        end
        transformed = rescale_hierarchy(U, unscaled; scaled = true)
        du, ds = similar(U), similar(U)
        heom!(du, U, unscaled, 0.0)
        heom!(ds, transformed, scaled, 0.0)
        @test ds ≈ rescale_hierarchy(du, unscaled; scaled = true) rtol = 2e-14
        factors = exp.(unscaled.log_scales)
        @test Matrix(scaled.transfer) ≈
              Diagonal(1 ./ factors) * Matrix(unscaled.transfer) * Diagonal(factors)
        @test scaled.bath.weights == bath.weights
        @test scaled.bath.mixing == bath.mixing
    end
    # Intermediate scale ratios can overflow even though the full transfer is finite.
    bath = ExponentialBath(
        [5e-324, complex(1.3e308, -1.3e308)],
        [1.0, 1.0];
        mixing = [0 1e-315; 0 0],
    )
    op = heom_operator(
        grid;
        mass = 1,
        potential = q -> zero(q),
        bath,
        depth = 3,
        scaled = true,
    )
    a = findfirst(==([1, 0]), op.indices)
    b = findfirst(==([0, 1]), op.indices)
    expected =
        -BigFloat(bath.mixing[1, 2]) * sqrt(
            abs(Complex{BigFloat}(bath.coefficients[2])) /
            BigFloat(real(bath.coefficients[1])),
        )
    @test op.transfer[a, b] ≈ expected rtol = 3e-13
end
