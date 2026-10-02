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
