# The hard depth cutoff is unstable on a finite box: the upward coordinate coupling
# grows with |q|, and the truncated hierarchy has growing modes near the box edge.
# See prototypes/box_edge for the full study.

# Dense generator of either discretisation, built column by column from heom!.
function dense_heom_generator(op)
    shape = (size(op.grid)..., length(hierarchy_indices(op)))
    L = zeros(prod(shape), prod(shape))
    U, dU = zeros(shape), zeros(shape)
    for j in eachindex(U)
        U[j] = 1
        heom!(dU, U, op, 0.0)
        L[:, j] .= vec(dU)
        U[j] = 0
    end
    return L
end

@testset "Hierarchy stability indicator" begin
    cold = drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 1)
    warm = drude_lorentz_pade_bath(; reorganization = 0.2, cutoff = 1.0, kT = 1.0, pade = 1)
    grid = PhaseSpaceGrid((-6, 6), 12, (-6, 6), 12)
    harmonic = harmonic_potential(; mass = 1.0, omega = 1.0)
    fd = (; discretization = FiniteDifference(4), moyal_terms = 1)
    operator(bath, depth; kwargs...) =
        heom_operator(grid; mass = 1.0, potential = harmonic, bath, depth, kwargs...)

    # Without auxiliaries, or without a bath, nothing can grow.
    for op in (operator(cold, 0), operator(ExponentialBath(ComplexF64[], Float64[]), 3))
        s = hierarchy_stability(op)
        @test s.rate == 0
        @test isnan(s.position) && isnan(s.wavenumber)
        @test s.radius == Inf
    end

    # Here the frozen rate exceeds the spectral abscissa of the complete generator, which
    # is positive: the hard cutoff has growing modes, concentrated at large |q|.
    for (bath, depth, options) in
        ((cold, 1, fd), (warm, 2, fd), (cold, 2, (;)), (warm, 2, (;)))
        op = operator(bath, depth; options..., scaled = true)
        s = hierarchy_stability(op)
        F = eigen(options == fd ? Matrix(sparse(op)) : dense_heom_generator(op))
        leading = argmax(real.(F.values))
        abscissa = real(F.values[leading])
        @test 0 < abscissa <= s.rate < 3abscissa
        @test s.position == 6
        @test 0 < s.radius < 6
        weight = sum(reshape(abs2.(F.vectors[:, leading]), size(grid)..., :); dims = (2, 3))
        @test sum(weight[abs.(grid.q) .> 3]) > 0.75sum(weight)
        # Amplitude scaling is a diagonal similarity transformation.
        @test hierarchy_stability(operator(bath, depth; options..., scaled = false)).rate ≈
              s.rate rtol = 1e-10
    end

    # Growth increases with depth and with the half-width of the box.
    rates = [hierarchy_stability(operator(cold, depth; fd...)).rate for depth in 1:4]
    @test issorted(rates) && rates[1] > 0
    widths = map((4.0, 6.0, 8.0)) do L
        wide = PhaseSpaceGrid((-L, L), round(Int, 2L), (-6, 6), 12)
        op = heom_operator(wide; mass = 1.0, potential = harmonic, bath = cold, depth = 2)
        hierarchy_stability(op).rate
    end
    @test issorted(widths) && widths[1] > 0
end

@testset "Cold, strong HEOM: long-time box-edge growth (known defect)" begin
    mass, omega = 1.0, 1.0
    # The bath of the cold, strong transient benchmarks, which stop at t = 0.8.
    bath = drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 2)
    grid = PhaseSpaceGrid((-6, 6), 32, (-6, 6), 32)
    mean0, covariance0 = [1.0, 0.0], [0.5 0.0; 0.0 0.5]
    W0 = on_grid(
        (q, p) -> gaussian_wigner(q, p; mean = mean0, covariance = covariance0),
        grid,
    )
    op = heom_operator(
        grid;
        mass,
        potential = harmonic_potential(; mass, omega),
        bath,
        depth = 2,
        scaled = true,
    )
    s = hierarchy_stability(op)
    @test s.rate > 1
    @test s.radius < 1.5
    limit = 1e3 * maximum(abs, W0)
    sol = solve(
        heom_problem(W0, (0.0, 30.0), op),
        Vern7();
        abstol = 1e-9,
        reltol = 1e-9,
        saveat = 0:0.5:30,
        unstable_check = (dt, u, p, t) -> !(maximum(abs, u) < limit),
    )
    # Second moments are exact at depth 2 for a harmonic potential, so their error
    # measures only the discretisation and the growing modes.
    errors = map(zip(sol.t, sol.u)) do (t, U)
        mean, covariance = cold_memory_reference(t, mean0, covariance0, bath; mass, omega)
        W = physical_wigner(U)
        return max(
            maximum(abs, phase_space_mean(W, grid) - mean),
            maximum(abs, phase_space_covariance(W, grid) - covariance),
        )
    end
    @test first(errors) < 1e-9
    # Known defect: the growing modes, seeded by the tails of the physical state,
    # corrupt the moments within a few time units and exceed 10³ times the initial
    # amplitude near the box edge before t = 30. A fix should make these pass.
    @test_broken maximum(errors) < 1e-6
    @test_broken last(sol.t) == 30
    if last(sol.t) < 30
        U = last(sol.u)
        @test abs(grid.q[argmax(abs.(U))[1]]) > 3
        @test 10 < last(sol.t) * s.rate < 60
    end
end
