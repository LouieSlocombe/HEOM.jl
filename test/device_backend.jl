# The batched hierarchy engine and the device code paths, on the KernelAbstractions CPU
# backend (`Array` with `batched = true`) and on JLArrays, the reference GPU array type.
# The shim turns scalar indexing and host/device mixing into errors.
include("jlarray_shim.jl")

@testset "Batched hierarchy and device arrays" begin
    @test !isnothing(Base.get_extension(HEOM, :HEOMKernelAbstractionsExt))
    grid = PhaseSpaceGrid((-6.0, 6.0), 24, (-6.0, 6.0), 20)
    V = harmonic_potential(; mass = 1.2, omega = 0.9)
    W0 = on_grid(
        (q, p) -> coherent_wigner(q, p; q0 = 1.0, p0 = 0.5, mass = 1.2, omega = 0.9),
        grid,
    )
    relative(a, b) = maximum(abs, Array(a) - b) / maximum(abs, b)
    function device_state(device, T, x)
        return HEOM.to_device(device, T, x)
    end

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
    baths = Dict(
        "Drude with terminator" => drude_lorentz_bath(;
            reorganization = 0.12,
            cutoff = 1.2,
            kT = 0.8,
            hbar = 1.0,
            matsubara = 2,
            terminator = true,
        ),
        "Padé" => drude_lorentz_pade_bath(;
            reorganization = 0.3,
            cutoff = 1.0,
            kT = 0.5,
            pade = 2,
        ),
        "Brownian oscillator" => brownian_oscillator_bath(;
            reorganization = 0.35,
            frequency = 1.7,
            damping = 0.6,
            kT = 0.4,
            hbar = 0.8,
            matsubara = 2,
        ),
        "Jordan and rotation blocks" => combine_baths(
            jordan,
            rotation,
            ExponentialBath([0.2 - 0.1im], [1.4]; hbar = 0.7),
        ),
    )

    @testset "Batched RHS equals the per-member loop" begin
        for (name, bath) in baths, scaled in (false, true), device in (Array, JLArray)
            kwargs = (; mass = 1.2, potential = V, bath, depth = 2, scaled)
            loop = heom_operator(grid; kwargs...)
            batched = heom_operator(grid; kwargs..., device, batched = true)
            @test batched.batch isa HEOM.BatchedHierarchy
            U = randn(size(grid)..., length(loop.indices))
            expected = hierarchy_rhs(loop, U)
            dU = similar(device_state(device, Float64, U))
            heom!(dU, device_state(device, Float64, U), batched, 0.0)
            @test dU isa device
            @test relative(dU, expected) < 1e-13
            # Float32 on the device, rounded once from the Float64 setup.
            batched32 = heom_operator(grid; kwargs..., device, float_type = Float32)
            U32 = device_state(device, Float32, U)
            dU32 = similar(U32)
            heom!(dU32, U32, batched32, 0.0f0)
            @test eltype(dU32) === Float32
            @test relative(dU32, expected) < 1e-5
        end
        # The loop path also runs on a device array, one member at a time.
        loop =
            heom_operator(grid; mass = 1.2, potential = V, bath = baths["Padé"], depth = 2)
        slow = heom_operator(
            grid;
            mass = 1.2,
            potential = V,
            bath = baths["Padé"],
            depth = 2,
            device = JLArray,
            batched = false,
        )
        @test isnothing(slow.batch)
        U = randn(size(grid)..., length(loop.indices))
        dU = similar(JLArray(U))
        heom!(dU, JLArray(U), slow, 0.0)
        @test relative(dU, hierarchy_rhs(loop, U)) < 1e-13
        # Size checks still come first.
        @test_throws DimensionMismatch heom!(dU, JLArray(U[:, :, 2:end]), slow, 0.0)
    end

    @testset "Driven hierarchies on a device" begin
        field(t) = 0.3 * sin(1.3t)
        for potential in (
                DrivenPotential(V, field, q -> q),
                TimeDependentPotential((q, t) -> 0.5q^2 + 0.2 * sin(t) * q),
            ),
            device in (Array, JLArray),
            T in (Float64, Float32)

            kwargs =
                (; mass = 1.2, potential, bath = baths["Padé"], depth = 2, scaled = true)
            loop = heom_operator(grid; kwargs...)
            batched = heom_operator(grid; kwargs..., device, float_type = T, batched = true)
            U = randn(size(grid)..., length(loop.indices))
            tolerance = T === Float64 ? 1e-13 : 1e-5
            # The generator is refreshed at every call, also when time runs backwards.
            for t in (0.0, 0.7, 0.2)
                expected = similar(U)
                heom!(expected, U, loop, t)
                Ud = device_state(device, T, U)
                dU = similar(Ud)
                heom!(dU, Ud, batched, T(t))
                @test relative(dU, expected) < tolerance
                tgrad, tgrad_expected = similar(Ud), similar(U)
                HEOM.heom_tgrad!(tgrad, Ud, batched, T(t))
                HEOM.heom_tgrad!(tgrad_expected, U, loop, t)
                @test relative(tgrad, tgrad_expected) < tolerance
            end
        end
        # A static hierarchy has no time derivative.
        static = heom_operator(
            grid;
            mass = 1.2,
            potential = V,
            bath = baths["Padé"],
            depth = 1,
            device = JLArray,
        )
        U = JLArray(randn(size(grid)..., length(static.indices)))
        dT = similar(U)
        HEOM.heom_tgrad!(dT, U, static, 0.0)
        @test all(iszero, Array(dT))
    end

    @testset "Single-slice operators on a device" begin
        exact = wigner_moyal_operator(grid; mass = 1.2, potential = V)
        device_op = wigner_moyal_operator(
            grid;
            mass = 1.2,
            potential = V,
            device = JLArray,
            float_type = Float32,
        )
        W = JLArray(Float32.(W0))
        dW = similar(W)
        wigner_moyal!(dW, W, device_op, 0.0f0)
        @test relative(dW, rhs(exact, W0)) < 1e-5
        damped = caldeira_leggett_operator(
            grid;
            mass = 1.2,
            potential = V,
            friction = 0.1,
            kT = 0.5,
        )
        damped_device = caldeira_leggett_operator(
            grid;
            mass = 1.2,
            potential = V,
            friction = 0.1,
            kT = 0.5,
            device = JLArray,
        )
        expected = similar(W0)
        caldeira_leggett!(expected, W0, damped, 0.0)
        dW64 = similar(JLArray(W0))
        caldeira_leggett!(dW64, JLArray(W0), damped_device, 0.0)
        @test relative(dW64, expected) < 1e-13
        prob = wigner_moyal_problem(W0, (0.0, 0.2), device_op)
        @test prob.u0 isa JLArray{Float32}
        @test prob.tspan === (0.0f0, 0.2f0)
        @test caldeira_leggett_problem(W0, (0.0, 0.2), damped_device).u0 isa
              JLArray{Float64}
        # A driven potential acts through its dipole symbol on device arrays too.
        driven = wigner_moyal_operator(
            grid;
            mass = 1.2,
            potential = DrivenPotential(V, t -> 0.2t, q -> q^2),
            device = JLArray,
        )
        host_driven = wigner_moyal_operator(
            grid;
            mass = 1.2,
            potential = DrivenPotential(V, t -> 0.2t, q -> q^2),
        )
        out, out_device = similar(W0), similar(JLArray(W0))
        HEOM.potential_action!(out, W0, host_driven, 0.4)
        HEOM.potential_action!(out_device, JLArray(W0), driven, 0.4)
        @test relative(out_device, out) < 1e-13
    end

    @testset "Problems, solvers and analysis on a device" begin
        bath = baths["Drude with terminator"]
        kwargs = (; mass = 1.2, potential = V, bath, depth = 3, scaled = true)
        host = heom_operator(grid; kwargs...)
        for device in (Array, JLArray)
            op = heom_operator(grid; kwargs..., device, batched = true)
            prob_host = heom_problem(W0, (0.0, 0.5), host)
            prob = heom_problem(W0, (0.0, 0.5), op)
            @test prob.u0 isa device{Float64,3}
            @test_throws ArgumentError heom_problem(W0, (0.0, 0.5), op; jacobian = :sparse)
            options = (; saveat = 0.25, abstol = 1e-10, reltol = 1e-10)
            sol_host = solve(prob_host, Vern7(); options...)
            sol = solve(prob, Vern7(); options...)
            @test relative(sol.u[end], sol_host.u[end]) < 1e-9
            @test physical_wigner(sol) isa AbstractMatrix
            @test relative(physical_wigner(sol), physical_wigner(sol_host)) < 1e-9

            # Saving only the physical member keeps device memory small.
            reduced =
                solve(prob, Vern7(); options..., save_idxs = HEOM.physical_save_idxs(op))
            @test size(reduced.u[end]) == (prod(size(grid)),)
            @test relative(physical_wigner(reduced), physical_wigner(sol_host)) < 1e-9
            @test HEOM.physical_state(reduced.u[end], op) isa Matrix{Float64}
            measured, expected =
                diagnostics(reduced; potential = V), diagnostics(sol_host; potential = V)
            @test all(
                maximum(abs, a - b) < 1e-8 for
                (a, b) in zip(values(measured), values(expected))
            )

            # Observables take device states and return host Float64 results.
            W = physical_wigner(sol)
            Wh = physical_wigner(sol_host)
            @test phase_space_integral(W, grid) ≈ phase_space_integral(Wh, grid)
            @test expectation((q, p) -> q^2 + p, W, grid) ≈
                  expectation((q, p) -> q^2 + p, Wh, grid)
            @test purity(W, grid) ≈ purity(Wh, grid)
            @test overlap(W, W0 .+ 0.0, grid) ≈ overlap(Wh, W0, grid)
            @test position_density(W, grid) isa Vector{Float64}
            @test position_density(W, grid) ≈ position_density(Wh, grid)
            @test momentum_density(W, grid) ≈ momentum_density(Wh, grid)
            @test phase_space_mean(W, grid) ≈ phase_space_mean(Wh, grid)
            @test phase_space_covariance(W, grid) ≈ phase_space_covariance(Wh, grid)
            @test energy(W, grid; mass = 1.2, potential = V) ≈
                  energy(Wh, grid; mass = 1.2, potential = V)
            @test wigner_negativity(W, grid) ≈ wigner_negativity(Wh, grid)
            @test all(values(boundary_weight(W, grid)) .≈ values(boundary_weight(Wh, grid)))
            @test all(values(spectral_tail(W, grid)) .≈ values(spectral_tail(Wh, grid)))
            @test probability(W, grid; q = (0.0, Inf)) ≈
                  probability(Wh, grid; q = (0.0, Inf))
            @test probability_current(W, grid; mass = 1.2) ≈
                  probability_current(Wh, grid; mass = 1.2)
            @test eltype(HEOM.plotting_state(W, grid)[1]) === Float64
            @test !(HEOM.plotting_state(W, grid)[1] isa JLArray)
            @test hierarchy_stability(op).rate ≈ hierarchy_stability(host).rate

            # Rescaling returns an array of the input's type.
            unscaled = rescale_hierarchy(sol.u[end], op; scaled = false)
            @test unscaled isa typeof(sol.u[end])
            @test relative(
                unscaled,
                rescale_hierarchy(sol_host.u[end], host; scaled = false),
            ) < 1e-12
            @test Array(rescale_hierarchy(sol.u[end], op; scaled = true)) ==
                  Array(sol.u[end])
        end
    end

    @testset "Equilibrium and linear response on a device" begin
        bath = baths["Drude with terminator"]
        kwargs = (; mass = 1.2, potential = V, bath, depth = 2, scaled = true)
        host = heom_operator(grid; kwargs...)
        options = (; abstol = 1e-10, reltol = 1e-10)
        equilibrium_host = equilibrate(
            W0,
            (0.0, 0.5),
            host,
            Vern7();
            stationarity_abstol = 1e-3,
            options...,
        )
        for device in (Array, JLArray)
            op = heom_operator(grid; kwargs..., device, batched = true)
            equilibrium = equilibrate(
                W0,
                (0.0, 0.5),
                op,
                Vern7();
                stationarity_abstol = 1e-3,
                options...,
            )
            @test equilibrium.hierarchy isa device{Float64,3}
            @test equilibrium.status === equilibrium_host.status
            @test relative(equilibrium.hierarchy, equilibrium_host.hierarchy) < 1e-9
            @test maximum(
                abs,
                equilibrium.residuals.absolute - equilibrium_host.residuals.absolute,
            ) < 1e-9
            @test physical_wigner(equilibrium) isa AbstractMatrix
            @test heom_problem(equilibrium, (0.0, 0.1)).u0 isa device{Float64,3}

            response_options = (; saveat = 0.25, options...)
            expected = linear_response(
                equilibrium_host.hierarchy,
                (0.0, 0.5),
                host,
                Vern7();
                response_options...,
            )
            response = linear_response(
                equilibrium.hierarchy,
                (0.0, 0.5),
                op,
                Vern7();
                response_options...,
            )
            @test maximum(abs, response.response - expected.response) < 1e-9
            @test response.times == expected.times
            # Only the physical member is saved from a device hierarchy.
            expected_size = device === Array ? ndims(expected.solution.u[1]) : 1
            @test ndims(response.solution.u[1]) == expected_size
            # Single-slice operators measure the response too.
            single = wigner_moyal_operator(grid; mass = 1.2, potential = V, device)
            single_host = wigner_moyal_operator(grid; mass = 1.2, potential = V)
            @test maximum(
                abs,
                linear_response(W0, (0.0, 0.5), single, Vern7(); response_options...).response -
                linear_response(W0, (0.0, 0.5), single_host, Vern7(); response_options...).response,
            ) < 1e-9
        end
    end

    @testset "Finite differences as stencils" begin
        fd = (; discretization = FiniteDifference(4), moyal_terms = 2)
        quartic(q) = V(q) + 0.05q^4
        for (device, T) in ((JLArray, Float64), (JLArray, Float32), (Array, Float32))
            tolerance = T === Float64 ? 1e-13 : 1e-5
            host = wigner_moyal_operator(grid; mass = 1.2, potential = quartic, fd...)
            stencil = wigner_moyal_operator(
                grid;
                mass = 1.2,
                potential = quartic,
                fd...,
                device,
                float_type = T,
            )
            @test stencil isa HEOM.StencilWignerMoyal
            @test occursin("StencilWignerMoyal", sprint(show, stencil))
            W = device_state(device, T, W0)
            dW = similar(W)
            wigner_moyal!(dW, W, stencil, T(0))
            @test relative(dW, rhs(host, W0)) < tolerance
            # The same operator acts on a view into a hierarchy.
            U = device_state(device, T, cat(W0, 2 .* W0; dims = 3))
            slice = similar(U, size(grid)...)
            wigner_moyal!(slice, view(U, :, :, 2), stencil, T(0))
            @test relative(slice, 2 .* rhs(host, W0)) < tolerance

            damped = caldeira_leggett_operator(
                grid;
                mass = 1.2,
                potential = quartic,
                friction = 0.1,
                kT = 0.5,
                fd...,
            )
            damped_stencil = caldeira_leggett_operator(
                grid;
                mass = 1.2,
                potential = quartic,
                friction = 0.1,
                kT = 0.5,
                fd...,
                device,
                float_type = T,
            )
            @test occursin("StencilCaldeiraLeggett", sprint(show, damped_stencil))
            expected = similar(W0)
            caldeira_leggett!(expected, W0, damped, 0.0)
            caldeira_leggett!(dW, W, damped_stencil, T(0))
            @test relative(dW, expected) < tolerance
            prob = caldeira_leggett_problem(W0, (0.0, 0.1), damped_stencil)
            @test eltype(prob.u0) === T
            @test_throws ArgumentError sparse(stencil)
            @test_throws ArgumentError sparse(stencil, 0.0)
            @test_throws ArgumentError sparse(damped_stencil)
            @test_throws ArgumentError sparse(damped_stencil, 0.0)
        end

        field(t) = 0.3 * sin(1.3t)
        potentials = (
            V,
            DrivenPotential(V, field, q -> q),
            TimeDependentPotential((q, t) -> 0.5q^2 + 0.2 * sin(t) * q),
        )
        for potential in potentials
            kwargs = (;
                mass = 1.2,
                potential,
                bath = baths["Drude with terminator"],
                depth = 2,
                scaled = true,
                fd...,
            )
            loop = heom_operator(grid; kwargs...)
            U = randn(size(grid)..., length(loop.indices))
            for (device, T, batched) in (
                (JLArray, Float64, true),
                (JLArray, Float64, false),
                (Array, Float32, false),
                (Array, Float64, true),
            )
                op = heom_operator(grid; kwargs..., device, float_type = T, batched)
                @test op.hamiltonian isa HEOM.StencilWignerMoyal
                @test (op.batch isa HEOM.BatchedStencilHierarchy) == batched
                tolerance = T === Float64 ? 1e-13 : 1e-5
                Ud = device_state(device, T, U)
                for t in (0.0, 0.7)
                    expected = similar(U)
                    heom!(expected, U, loop, t)
                    dU = similar(Ud)
                    heom!(dU, Ud, op, T(t))
                    @test relative(dU, expected) < tolerance
                    tgrad, tgrad_expected = similar(Ud), similar(U)
                    HEOM.heom_tgrad!(tgrad, Ud, op, T(t))
                    HEOM.heom_tgrad!(tgrad_expected, U, loop, t)
                    HEOM.is_time_dependent(loop) &&
                        @test relative(tgrad, tgrad_expected) < tolerance
                end
            end
        end

        # The linear response and the stability indicator use the same operators.
        kwargs = (;
            mass = 1.2,
            potential = V,
            bath = baths["Drude with terminator"],
            depth = 2,
            scaled = true,
            fd...,
        )
        host = heom_operator(grid; kwargs...)
        device_op = heom_operator(grid; kwargs..., device = JLArray)
        options = (; saveat = 0.25, abstol = 1e-10, reltol = 1e-10)
        response = linear_response(
            heom_problem(W0, (0.0, 0.5), host).u0,
            (0.0, 0.5),
            host,
            Vern7();
            options...,
        )
        measured = linear_response(
            heom_problem(W0, (0.0, 0.5), device_op).u0,
            (0.0, 0.5),
            device_op,
            Vern7();
            options...,
        )
        @test maximum(abs, measured.response - response.response) < 1e-9
        @test hierarchy_stability(device_op).rate ≈ hierarchy_stability(host).rate
        @test_throws ArgumentError heom_problem(
            W0,
            (0.0, 0.1),
            device_op;
            jacobian = :sparse,
        )
        @test_throws ArgumentError sparse(device_op)
        # A potential action on a device, with a driven dipole.
        driven = DrivenPotential(V, t -> 0.2t, q -> q^2)
        stencil = wigner_moyal_operator(
            grid;
            mass = 1.2,
            potential = driven,
            fd...,
            device = JLArray,
        )
        reference = wigner_moyal_operator(grid; mass = 1.2, potential = driven, fd...)
        out, out_device = similar(W0), similar(JLArray(W0))
        HEOM.potential_action!(out, W0, reference, 0.4)
        HEOM.potential_action!(out_device, JLArray(W0), stencil, 0.4)
        @test relative(out_device, out) < 1e-13
    end

    @testset "Errors" begin
        @test_throws ArgumentError HEOM.require_backend(nothing)
        @test isnothing(
            HEOM.require_backend(Base.get_extension(HEOM, :HEOMKernelAbstractionsExt)),
        )
    end
end
