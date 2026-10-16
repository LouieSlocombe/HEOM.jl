@testset "Composite baths" begin
    @testset "Generalized components retain their independent bases" begin
        jordan = ExponentialBath(
            [0.4 - 0.2im, 0.3],
            [0.8, 0.8];
            weights = [1.0, 0.0],
            mixing = [0.0 1.0; 0.0 0.0],
            hbar = 0.7,
            diffusion = 0.04,
            counterterm = 0.2,
        )
        rotation = ExponentialBath(
            [0.6, -0.1im],
            [0.3, 0.3];
            weights = [1.0, 0.0],
            mixing = [0.0 -1.5; 1.5 0.0],
            hbar = 0.7,
            diffusion = 0.06,
            counterterm = 0.3,
        )
        exponential = ExponentialBath([0.2 - 0.1im], [1.4]; hbar = 0.7)
        components = [jordan, rotation, exponential]
        combined = combine_baths(components...)
        from_vector = combine_baths(components)
        @test combined.coefficients == vcat((b.coefficients for b in components)...)
        @test combined.rates == vcat((b.rates for b in components)...)
        @test combined.weights == vcat((b.weights for b in components)...)
        @test combined.mixing isa SparseMatrixCSC
        @test combined.mixing == blockdiag((b.mixing for b in components)...)
        @test nnz(combined.mixing) == sum(nnz(b.mixing) for b in components)
        @test combined.hbar == 0.7
        @test combined.diffusion ≈ 0.1
        @test combined.counterterm == 0.5
        for field in fieldnames(ExponentialBath)
            @test getfield(from_vector, field) == getfield(combined, field)
        end
        for t in (0.0, 0.2, 1.0, 4.0)
            @test bath_correlation(combined, t) ≈
                  sum(bath_correlation(b, t) for b in components) atol = 5e-15
        end
        for ω in (-2.0, -0.1, 0.0, 0.1, 2.0)
            @test bath_spectrum(combined, ω) ≈ sum(bath_spectrum(b, ω) for b in components) atol =
                5e-15
        end

        # Editing one component after composition cannot change the result.
        saved = deepcopy(combined)
        jordan.coefficients[1] = 9
        jordan.rates[1] = 9
        jordan.weights[1] = 9
        jordan.mixing[1, 2] = 9
        for field in fieldnames(ExponentialBath)
            @test getfield(combined, field) == getfield(saved, field)
        end
    end

    @testset "Single components and empty mode sets are copied" begin
        source = ExponentialBath(
            [1.0, 0.2im],
            [0.5, 0.5];
            mixing = [0.0 1.0; 0.0 0.0],
            weights = [1.0, 0.0],
        )
        for combined in (combine_baths(source), combine_baths([source]))
            for field in (:coefficients, :rates, :weights, :mixing)
                @test getfield(combined, field) == getfield(source, field)
                @test getfield(combined, field) !== getfield(source, field)
            end
            @test combined.mixing.nzval !== source.mixing.nzval
            @test combined.mixing.rowval !== source.mixing.rowval
            @test combined.mixing.colptr !== source.mixing.colptr
        end

        empty = ExponentialBath(ComplexF64[], Float64[]; diffusion = 0.1, counterterm = 0.2)
        for combined in (combine_baths(empty), combine_baths([empty, empty]))
            @test isempty(combined.coefficients)
            @test isempty(combined.rates)
            @test isempty(combined.weights)
            @test size(combined.mixing) == (0, 0)
            @test combined.coefficients !== empty.coefficients
            @test combined.rates !== empty.rates
            @test combined.weights !== empty.weights
            @test combined.mixing !== empty.mixing
        end
        combined = combine_baths(empty, source, empty)
        @test combined.coefficients == source.coefficients
        @test combined.mixing == source.mixing
        @test combined.diffusion == 0.2
        @test combined.counterterm == 0.4
    end

    @testset "Invalid combinations" begin
        @test_throws ArgumentError combine_baths()
        @test_throws ArgumentError combine_baths(ExponentialBath[])
        one = ExponentialBath([1.0], [1.0])
        different = ExponentialBath([1.0], [1.0]; hbar = nextfloat(1.0))
        @test_throws ArgumentError combine_baths(one, different)
    end
end
