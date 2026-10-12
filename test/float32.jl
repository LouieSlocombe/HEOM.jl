@testset "Float32 precision and device keywords" begin
    grid = PhaseSpaceGrid((-6.0, 6.0), 24, (-6.0, 6.0), 20)
    V = harmonic_potential(; mass = 1.2, omega = 0.9)
    W0 = on_grid(
        (q, p) -> coherent_wigner(q, p; q0 = 1.0, p0 = 0.5, mass = 1.2, omega = 0.9),
        grid,
    )
    bath = drude_lorentz_bath(;
        reorganization = 0.12,
        cutoff = 1.2,
        kT = 0.8,
        hbar = 1.0,
        matsubara = 2,
        terminator = true,
    )
    relative(a, b) = maximum(abs, Float64.(a) - b) / maximum(abs, b)

    @testset "Keywords are validated" begin
        build(; kwargs...) =
            wigner_moyal_operator(grid; mass = 1.2, potential = V, kwargs...)
        @test_throws ArgumentError build(float_type = Float16)
        @test_throws ArgumentError build(float_type = Int)
        @test_throws ArgumentError build(float_type = 1)
        @test_throws ArgumentError build(device = Int)
        @test_throws ArgumentError build(device = 3)
        # Finite differences become matrix-free stencils away from Float64 on the host.
        fd = (; discretization = FiniteDifference(2), moyal_terms = 2)
        stencils = build(; fd..., float_type = Float32)
        @test stencils isa HEOM.StencilWignerMoyal
        @test HEOM.element_type(stencils) === Float32
        @test build(; fd...) isa HEOM.FiniteDifferenceWignerMoyal
        # The default reproduces the host Float64 operator, a subtype of Array is the same.
        op = build()
        @test HEOM.move_operator(op, Array, Float64) === op
        @test build(device = Array{Float32}) isa typeof(op)
        @test HEOM.element_type(op) === Float64
        @test HEOM.element_type(build(float_type = Float32)) === Float32
    end

    @testset "Wigner–Moyal and Caldeira–Leggett in Float32" begin
        exact = wigner_moyal_operator(grid; mass = 1.2, potential = V)
        single =
            wigner_moyal_operator(grid; mass = 1.2, potential = V, float_type = Float32)
        @test eltype(single.W) === Float32
        @test eltype(single.potential_symbol) === ComplexF32
        W32 = Float32.(W0)
        dW = similar(W32)
        wigner_moyal!(dW, W32, single, 0.0f0)
        @test eltype(dW) === Float32
        @test relative(dW, rhs(exact, W0)) < 1e-5

        prob = wigner_moyal_problem(W0, (0.0, 0.1), single)
        @test eltype(prob.u0) === Float32
        @test prob.tspan === (0.0f0, 0.1f0)
        sol = solve(prob, Vern7(); abstol = 1.0f-6, reltol = 1.0f-6)
        reference = solve(
            wigner_moyal_problem(W0, (0.0, 0.1), exact),
            Vern7();
            abstol = 1e-10,
            reltol = 1e-10,
        )
        @test eltype(sol.u[end]) === Float32
        @test relative(sol.u[end], reference.u[end]) < 1e-4

        damped = caldeira_leggett_operator(
            grid;
            mass = 1.2,
            potential = V,
            friction = 0.1,
            kT = 0.5,
        )
        damped32 = caldeira_leggett_operator(
            grid;
            mass = 1.2,
            potential = V,
            friction = 0.1,
            kT = 0.5,
            float_type = Float32,
        )
        d64, d32 = similar(W0), similar(W32)
        caldeira_leggett!(d64, W0, damped, 0.0)
        caldeira_leggett!(d32, W32, damped32, 0.0f0)
        @test relative(d32, d64) < 1e-5
        problem = caldeira_leggett_problem(W0, (0.0, 0.1), damped32)
        @test eltype(problem.u0) === Float32
        @test HEOM.element_type(damped32) === Float32
    end

    @testset "HEOM in Float32" begin
        for scaled in (false, true)
            op = heom_operator(grid; mass = 1.2, potential = V, bath, depth = 3, scaled)
            op32 = heom_operator(
                grid;
                mass = 1.2,
                potential = V,
                bath,
                depth = 3,
                scaled,
                float_type = Float32,
            )
            @test HEOM.element_type(op32) === Float32
            @test eltype(op32.scratch) === Float32
            @test eltype(op32.qcol) === Float32
            U = randn(size(grid)..., length(op.indices))
            expected = hierarchy_rhs(op, U)
            dU = similar(Float32.(U))
            heom!(dU, Float32.(U), op32, 0.0f0)
            @test eltype(dU) === Float32
            @test relative(dU, expected) < 1e-5
        end

        op = heom_operator(grid; mass = 1.2, potential = V, bath, depth = 3, scaled = true)
        op32 = heom_operator(
            grid;
            mass = 1.2,
            potential = V,
            bath,
            depth = 3,
            scaled = true,
            float_type = Float32,
        )
        prob = heom_problem(W0, (0.0, 0.5), op)
        prob32 = heom_problem(W0, (0.0, 0.5), op32)
        @test eltype(prob32.u0) === Float32
        @test prob32.tspan === (0.0f0, 0.5f0)
        full = heom_problem(prob.u0, (0.0, 0.5), op32)
        @test eltype(full.u0) === Float32 && full.u0 == Float32.(prob.u0)
        sol = solve(prob, Vern7(); saveat = 0.25, abstol = 1e-10, reltol = 1e-10)
        sol32 = solve(prob32, Vern7(); saveat = 0.25f0, abstol = 1.0f-6, reltol = 1.0f-6)
        @test relative(sol32.u[end], sol.u[end]) < 1e-4
        # Analysis functions return Float64 results from a Float32 state.
        W = physical_wigner(sol32)
        reference = physical_wigner(sol)
        @test phase_space_integral(W, grid) isa Float64
        integral, integral_reference = phase_space_integral(W, grid), 1.0
        integral_reference = phase_space_integral(reference, grid)
        @test isapprox(integral, integral_reference; atol = 1e-4)
        energy32 = energy(W, grid; mass = 1.2, potential = V)
        energy64 = energy(reference, grid; mass = 1.2, potential = V)
        @test isapprox(energy32, energy64; rtol = 1e-3)
        measured = diagnostics(sol32; potential = V)
        expected = diagnostics(sol; potential = V)
        @test all(
            maximum(abs, a - b) < 2e-3 for (a, b) in zip(values(measured), values(expected))
        )

        # A stationarity check in the state's own time type must not report a time-out as
        # an early termination, even when the end time is not representable.
        result = equilibrate(
            W0,
            (0.0, 0.1),
            op32,
            Vern7();
            stationarity_abstol = 1e-12,
            abstol = 1.0f-6,
            reltol = 1.0f-6,
        )
        @test eltype(result.hierarchy) === Float32
        @test result.status === :time_limit
        eq = equilibrate(
            W0,
            (0.0, 0.2),
            op,
            Vern7();
            stationarity_abstol = 1e-3,
            abstol = 1e-9,
            reltol = 1e-9,
        )
        eq32 = equilibrate(
            W0,
            (0.0, 0.2),
            op32,
            Vern7();
            stationarity_abstol = 1e-3,
            abstol = 1.0f-6,
            reltol = 1.0f-6,
        )
        @test relative(eq32.hierarchy, eq.hierarchy) < 1e-3
        response = linear_response(
            eq.hierarchy,
            (0.0, 0.5),
            op,
            Vern7();
            saveat = 0.25,
            abstol = 1e-9,
            reltol = 1e-9,
        )
        response32 = linear_response(
            eq32.hierarchy,
            (0.0, 0.5),
            op32,
            Vern7();
            saveat = 0.25f0,
            abstol = 1.0f-6,
            reltol = 1.0f-6,
        )
        @test maximum(abs, response32.response - response.response) < 1e-3
        @test eltype(HEOM.heom_problem(eq32, (0.0, 0.1)).u0) === Float32

        # Rescaling keeps the array type and checks the range of the narrower type.
        converted = rescale_hierarchy(prob32.u0, op32; scaled = false)
        @test eltype(converted) === Float32
        @test relative(converted, rescale_hierarchy(prob.u0, op; scaled = false)) < 1e-6
        @test rescale_hierarchy(prob32.u0, op32; scaled = true) == prob32.u0
    end

    @testset "Drives in Float32" begin
        field(t) = 0.3 * sin(1.3t)
        for potential in (
            DrivenPotential(V, field, q -> q),
            TimeDependentPotential((q, t) -> 0.5q^2 + 0.2 * sin(t) * q),
        )
            op = heom_operator(grid; mass = 1.2, potential, bath, depth = 2, scaled = true)
            op32 = heom_operator(
                grid;
                mass = 1.2,
                potential,
                bath,
                depth = 2,
                scaled = true,
                float_type = Float32,
            )
            U = randn(size(grid)..., length(op.indices))
            for t in (0.0, 0.7)
                dU = similar(U)
                heom!(dU, U, op, t)
                dU32 = similar(Float32.(U))
                heom!(dU32, Float32.(U), op32, Float32(t))
                @test relative(dU32, dU) < 1e-5
                expected, measured = similar(U), similar(dU32)
                HEOM.heom_tgrad!(expected, U, op, t)
                HEOM.heom_tgrad!(measured, Float32.(U), op32, Float32(t))
                @test relative(measured, expected) < 1e-5
            end
        end
        @test HEOM.checked_drive_value(0.1, Float32) === 0.1f0
        @test_throws ArgumentError HEOM.checked_drive_value(1e300, Float32)
        @test HEOM.user_time(0.1f0) === Float64(0.1f0)
        @test HEOM.user_time(0.25) === 0.25
    end

    @testset "Coefficients outside the range of Float32" begin
        huge = ExponentialBath([1e37], [1.0])
        build(T) = heom_operator(
            grid;
            mass = 1,
            potential = V,
            bath = huge,
            depth = 40,
            float_type = T,
        )
        @test build(Float64) isa HEOM.WignerHEOM
        @test_throws ArgumentError build(Float32)
        # Beyond the range of Float32 only after the conversion of a Float64 hierarchy.
        wide = heom_operator(
            PhaseSpaceGrid((-1.0, 1.0), 4, (-1.0, 1.0), 4);
            mass = 1,
            potential = q -> zero(q),
            bath = ExponentialBath([1e6], [1.0]),
            depth = 20,
            scaled = true,
            float_type = Float32,
        )
        @test_throws ArgumentError rescale_hierarchy(
            ones(Float32, 4, 4, length(wide.indices)),
            wide;
            scaled = false,
        )
    end

    @testset "Host conversion helpers" begin
        W = rand(3, 3)
        @test HEOM.to_host(W) === W
        @test HEOM.to_host(view(W, :, 1)) isa SubArray
        @test HEOM.to_host(Float32.(W)) isa Matrix{Float64}
        @test HEOM.to_host(Float32.(W)) == Float64.(Float32.(W))
        # Complex arrays are copied, not narrowed, so that real-only checks still fire.
        @test HEOM.to_host(complex.(W)) == complex.(W)
        @test_throws ArgumentError HEOM.plotting_state(
            complex.(W),
            PhaseSpaceGrid((0, 1), 3, (0, 1), 3),
        )
        @test HEOM.convert_tspan(Float64, (0, 1)) === (0, 1)
        @test HEOM.convert_tspan(Float32, (0, 1.5)) === (0.0f0, 1.5f0)
        @test HEOM.physical_save_idxs(
            heom_operator(grid; mass = 1, potential = V, bath, depth = 1),
        ) == 1:prod(size(grid))
    end
end
