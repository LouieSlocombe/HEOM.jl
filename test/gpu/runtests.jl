# Equivalence of the CUDA device paths with the CPU reference. CI has no GPU, so this suite
# is run by hand on a machine with one:
#
#     julia --project=test/gpu -e 'using Pkg; Pkg.instantiate()'
#     julia --project=test/gpu test/gpu/runtests.jl
#
# The same cases run on a CPU in `test/device_backend.jl` with JLArrays; here they run on a
# real `CuArray` with scalar indexing disabled. Float64 must agree to rounding error, and
# Float32 to about 1e-5.
using ADTypes
using CUDA
using HEOM
using KernelAbstractions
using LinearSolve
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqVerner
using Test

CUDA.functional() || error("no functional CUDA device")
CUDA.allowscalar(false)

relative(a, b) = maximum(abs, Array(a) - b) / maximum(abs, b)
on_gpu(T, x) = CuArray(T.(x))

grid = PhaseSpaceGrid((-6.0, 6.0), 24, (-6.0, 6.0), 20)
V = harmonic_potential(; mass = 1.2, omega = 0.9)
W0 = on_grid(
    (q, p) -> coherent_wigner(q, p; q0 = 1.0, p0 = 0.5, mass = 1.2, omega = 0.9),
    grid,
)
drude = drude_lorentz_bath(;
    reorganization = 0.12,
    cutoff = 1.2,
    kT = 0.8,
    hbar = 1.0,
    matsubara = 2,
    terminator = true,
)
pade = drude_lorentz_pade_bath(; reorganization = 0.3, cutoff = 1.0, kT = 0.5, pade = 2)
brownian = brownian_oscillator_bath(;
    reorganization = 0.35,
    frequency = 1.7,
    damping = 0.6,
    kT = 0.4,
    hbar = 0.8,
    matsubara = 2,
)
baths = Dict("Drude" => drude, "Padé" => pade, "Brownian" => brownian)
field(t) = 0.3 * sin(1.3t)
fd = (; discretization = FiniteDifference(4), moyal_terms = 2)

@testset "CUDA" begin
    @testset "Single-slice operators" begin
        exact = wigner_moyal_operator(grid; mass = 1.2, potential = V)
        damped = caldeira_leggett_operator(
            grid;
            mass = 1.2,
            potential = V,
            friction = 0.1,
            kT = 0.5,
        )
        for T in (Float64, Float32)
            tolerance = T === Float64 ? 1e-13 : 1e-5
            W = on_gpu(T, W0)
            dW = similar(W)
            wigner_moyal!(
                dW,
                W,
                wigner_moyal_operator(
                    grid;
                    mass = 1.2,
                    potential = V,
                    device = CuArray,
                    float_type = T,
                ),
                T(0),
            )
            expected = similar(W0)
            wigner_moyal!(expected, W0, exact, 0.0)
            @test relative(dW, expected) < tolerance
            caldeira_leggett!(
                dW,
                W,
                caldeira_leggett_operator(
                    grid;
                    mass = 1.2,
                    potential = V,
                    friction = 0.1,
                    kT = 0.5,
                    device = CuArray,
                    float_type = T,
                ),
                T(0),
            )
            caldeira_leggett!(expected, W0, damped, 0.0)
            @test relative(dW, expected) < tolerance
        end
    end

    @testset "Batched RHS: $name, scaled = $scaled" for (name, bath) in baths,
        scaled in (false, true)

        kwargs = (; mass = 1.2, potential = V, bath, depth = 3, scaled)
        loop = heom_operator(grid; kwargs...)
        U = randn(size(grid)..., length(loop.indices))
        expected = similar(U)
        heom!(expected, U, loop, 0.0)
        for T in (Float64, Float32)
            op = heom_operator(grid; kwargs..., device = CuArray, float_type = T)
            @test op.batch isa HEOM.BatchedHierarchy
            Ud = on_gpu(T, U)
            dU = similar(Ud)
            heom!(dU, Ud, op, T(0))
            @test dU isa CuArray{T,3}
            @test relative(dU, expected) < (T === Float64 ? 1e-13 : 1e-5)
        end
    end

    @testset "The per-member loop on a device" begin
        kwargs = (; mass = 1.2, potential = V, bath = pade, depth = 2, scaled = true)
        loop = heom_operator(grid; kwargs...)
        slow = heom_operator(grid; kwargs..., device = CuArray, batched = false)
        U = randn(size(grid)..., length(loop.indices))
        expected, dU = similar(U), similar(CuArray(U))
        heom!(expected, U, loop, 0.0)
        heom!(dU, CuArray(U), slow, 0.0)
        @test relative(dU, expected) < 1e-13
    end

    @testset "Drives: $(nameof(typeof(potential)))" for potential in (
        DrivenPotential(V, field, q -> q),
        TimeDependentPotential((q, t) -> 0.5q^2 + 0.2 * sin(t) * q),
    )
        kwargs = (; mass = 1.2, potential, bath = pade, depth = 2, scaled = true)
        loop = heom_operator(grid; kwargs...)
        U = randn(size(grid)..., length(loop.indices))
        for T in (Float64, Float32)
            op = heom_operator(grid; kwargs..., device = CuArray, float_type = T)
            tolerance = T === Float64 ? 1e-13 : 1e-5
            for t in (0.0, 0.7, 0.2)
                expected, tgrad_expected = similar(U), similar(U)
                heom!(expected, U, loop, t)
                HEOM.heom_tgrad!(tgrad_expected, U, loop, t)
                Ud = on_gpu(T, U)
                dU, tgrad = similar(Ud), similar(Ud)
                heom!(dU, Ud, op, T(t))
                HEOM.heom_tgrad!(tgrad, Ud, op, T(t))
                @test relative(dU, expected) < tolerance
                @test relative(tgrad, tgrad_expected) < tolerance
            end
        end
    end

    @testset "Finite differences as stencils: $(nameof(typeof(potential)))" for potential in
                                                                                (
        V,
        DrivenPotential(V, field, q -> q),
    )
        kwargs = (; mass = 1.2, potential, bath = drude, depth = 2, scaled = true, fd...)
        loop = heom_operator(grid; kwargs...)
        U = randn(size(grid)..., length(loop.indices))
        for T in (Float64, Float32), batched in (true, false)
            op = heom_operator(grid; kwargs..., device = CuArray, float_type = T, batched)
            expected = similar(U)
            heom!(expected, U, loop, 0.7)
            Ud = on_gpu(T, U)
            dU = similar(Ud)
            heom!(dU, Ud, op, T(0.7))
            @test relative(dU, expected) < (T === Float64 ? 1e-13 : 1e-5)
        end
    end

    @testset "Solvers, analysis and workflows" begin
        kwargs = (; mass = 1.2, potential = V, bath = drude, depth = 3, scaled = true)
        host = heom_operator(grid; kwargs...)
        op = heom_operator(grid; kwargs..., device = CuArray)
        options = (; saveat = 0.25, abstol = 1e-10, reltol = 1e-10)
        reference = solve(heom_problem(W0, (0.0, 0.5), host), Vern7(); options...)
        sol = solve(heom_problem(W0, (0.0, 0.5), op), Vern7(); options...)
        @test sol.u[end] isa CuArray{Float64,3}
        @test relative(sol.u[end], reference.u[end]) < 1e-9
        @test physical_wigner(sol) isa AbstractMatrix
        measured, expected =
            diagnostics(sol; potential = V), diagnostics(reference; potential = V)
        @test all(
            maximum(abs, a - b) < 1e-8 for (a, b) in zip(values(measured), values(expected))
        )

        # A stiff solver with a matrix-free Krylov linear solve, entirely on the device.
        stiff = Rodas5P(
            autodiff = AutoFiniteDiff(),
            linsolve = KrylovJL_GMRES(),
            concrete_jac = false,
        )
        coarse = (; abstol = 1e-8, reltol = 1e-8)
        stiff_host = solve(heom_problem(W0, (0.0, 0.5), host), stiff; coarse...)
        stiff_device = solve(heom_problem(W0, (0.0, 0.5), op), stiff; coarse...)
        @test stiff_device.retcode == stiff_host.retcode
        @test relative(stiff_device.u[end], stiff_host.u[end]) < 1e-6

        # Only the physical member is saved, so device memory stays small.
        reduced = solve(
            heom_problem(W0, (0.0, 0.5), op),
            Vern7();
            options...,
            save_idxs = physical_save_idxs(op),
        )
        @test relative(physical_wigner(reduced), physical_wigner(reference)) < 1e-9

        equilibrium_host = equilibrate(
            W0,
            (0.0, 0.5),
            host,
            Vern7();
            stationarity_abstol = 1e-3,
            abstol = 1e-10,
            reltol = 1e-10,
        )
        equilibrium = equilibrate(
            W0,
            (0.0, 0.5),
            op,
            Vern7();
            stationarity_abstol = 1e-3,
            abstol = 1e-10,
            reltol = 1e-10,
        )
        @test equilibrium.hierarchy isa CuArray{Float64,3}
        @test equilibrium.status === equilibrium_host.status
        @test relative(equilibrium.hierarchy, equilibrium_host.hierarchy) < 1e-9
        response_options = (; saveat = 0.25, abstol = 1e-10, reltol = 1e-10)
        response_host = linear_response(
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
        @test maximum(abs, response.response - response_host.response) < 1e-9
        @test hierarchy_stability(op).rate ≈ hierarchy_stability(host).rate
    end

    @testset "Float32 propagation" begin
        kwargs = (; mass = 1.2, potential = V, bath = drude, depth = 3, scaled = true)
        host = heom_operator(grid; kwargs...)
        op = heom_operator(grid; kwargs..., device = CuArray, float_type = Float32)
        reference = solve(
            heom_problem(W0, (0.0, 0.5), host),
            Vern7();
            abstol = 1e-10,
            reltol = 1e-10,
        )
        sol = solve(
            heom_problem(W0, (0.0, 0.5), op),
            Vern7();
            abstol = 1.0f-5,
            reltol = 1.0f-5,
        )
        @test eltype(sol.u[end]) === Float32
        @test relative(sol.u[end], reference.u[end]) < 1e-3
    end
end
