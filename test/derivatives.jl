@testset "Finite-difference stencils" begin
    @test HEOM.central_difference_weights(1, 2) == [-1 // 2, 0, 1 // 2]
    @test HEOM.central_difference_weights(1, 4) == [1 // 12, -2 // 3, 0, 2 // 3, -1 // 12]
    @test HEOM.central_difference_weights(3, 2) == [-1 // 2, 1, 0, -1, 1 // 2]
    @test HEOM.central_difference_weights(2, 2) == [1, -2, 1]

    for order in (2, 4, 6), derivative in (1, 3)
        D = HEOM.periodic_difference_matrix(derivative, order, 32, 0.1)
        @test D isa SparseMatrixCSC{Float64,Int}
        # Odd central differences are exactly antisymmetric.
        @test iszero(norm(D + transpose(D)))

        # On a periodic function the error falls as hᵒʳᵈᵉʳ.
        errors = map((32, 64)) do n
            h = 2π / n
            x = h .* (0:(n-1))
            f = @. sin(x) + cos(2x) / 2
            exact = derivative == 1 ? (@. cos(x) - sin(2x)) : (@. -cos(x) + 4sin(2x))
            max_error(HEOM.periodic_difference_matrix(derivative, order, n, h) * f, exact)
        end
        @test log2(errors[1] / errors[2]) ≈ order atol = 0.5
    end

    # A seven-point stencil does not fit on six points.
    @test_throws ArgumentError HEOM.periodic_difference_matrix(3, 4, 6, 0.1)
    @test FiniteDifference().order == 4
    @test FiniteDifference(6).order == 6
    @test_throws ArgumentError FiniteDifference(3)
    @test_throws ArgumentError FiniteDifference(0)
end

@testset "Spectral wavenumbers" begin
    # The unpaired Nyquist wavenumber of an even grid is dropped.
    @test HEOM.wavenumbers(4, 0.5) ≈ [0, π, 0]
    @test HEOM.wavenumbers(5, 1.0) ≈ [0, 0.4π, 0.8π]
end
