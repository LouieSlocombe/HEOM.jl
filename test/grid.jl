@testset "Phase-space grid" begin
    grid = PhaseSpaceGrid((-8, 8), 64, (-6.0, 4.0), 40)
    @test size(grid) == (64, 40)
    @test grid.q isa Vector{Float64}
    @test grid.dq == 0.25
    @test grid.dp == 0.25
    @test first(grid.q) == -8 && last(grid.q) == 8 - grid.dq
    @test first(grid.p) == -6 && last(grid.p) == 4 - grid.dp
    @test diff(grid.q) ≈ fill(grid.dq, 63)
    @test repr(grid) == "PhaseSpaceGrid(q ∈ [-8.0, 8.0) × 64, p ∈ [-6.0, 4.0) × 40)"

    values = @inferred on_grid((q, p) -> q + 10p, grid)
    @test values isa Matrix{Float64}
    @test size(values) == size(grid)
    @test values[3, 5] == grid.q[3] + 10 * grid.p[5]

    @test_throws ArgumentError PhaseSpaceGrid((-1, 1), 1, (-1, 1), 8)
    @test_throws ArgumentError PhaseSpaceGrid((-1, 1), 8, (-1, 1), 1)
    @test_throws ArgumentError PhaseSpaceGrid((1, -1), 8, (-1, 1), 8)
    @test_throws ArgumentError PhaseSpaceGrid((-1, 1), 8, (-1, Inf), 8)
end

@testset "Observables" begin
    grid = PhaseSpaceGrid((-8, 8), 64, (-8, 8), 64)
    σq, σp, q0, p0 = 0.6, 0.9, 0.5, -1.0
    W = on_grid(
        (q, p) -> exp(-(q - q0)^2 / (2σq^2) - (p - p0)^2 / (2σp^2)) / (2π * σq * σp),
        grid,
    )
    @test phase_space_integral(W, grid) ≈ 1 atol = 1e-12
    @test expectation((q, p) -> q, W, grid) ≈ q0 atol = 1e-12
    @test expectation((q, p) -> p, W, grid) ≈ p0 atol = 1e-12
    @test expectation((q, p) -> (q - q0)^2, W, grid) ≈ σq^2 atol = 1e-12
    # A Gaussian with uncorrelated widths σq and σp has purity ħ/(2σqσp).
    @test purity(W, grid; hbar = 0.7) ≈ 0.7 / (2σq * σp) atol = 1e-12
    @test_throws DimensionMismatch phase_space_integral(zeros(3, 3), grid)
end
