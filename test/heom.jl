"""
Evaluate the full Wigner hierarchy into a fresh state array.
"""
function hierarchy_rhs(op, U)
    dU = similar(U)
    heom!(dU, U, op, 0.0)
    return dU
end

@testset "Exponential bath and Drude–Lorentz decomposition" begin
    coefficients = [0.4 - 0.2im, -0.1 + 0.03im]
    rates = [0.6, 1.7]
    bath = ExponentialBath(
        coefficients,
        rates;
        hbar = 0.7,
        diffusion = 0.05,
        counterterm = 0.2,
    )
    @test bath.coefficients == coefficients
    @test bath.rates == rates
    @test bath.hbar == 0.7
    @test bath.diffusion == 0.05
    @test bath.counterterm == 0.2
    coefficients[1] = 9
    rates[1] = 9
    @test bath.coefficients[1] == 0.4 - 0.2im
    @test bath.rates[1] == 0.6
    @test ExponentialBath(Float64[], Float64[]).coefficients == ComplexF64[]
    @test ExponentialBath([1], [2]).coefficients == ComplexF64[1]

    λ, γ, kT, ħ = 0.3, 0.8, 1.2, 0.7
    for terms in (0, 1, 4)
        drude = drude_lorentz_bath(;
            reorganization = λ,
            cutoff = γ,
            kT,
            hbar = ħ,
            matsubara = terms,
        )
        ν = [2π * k * kT / ħ for k in 1:terms]
        expected = ComplexF64[λ*ħ*γ*(cot(ħ*γ/(2kT))-im)]
        append!(expected, [4λ * γ * kT * v / (v^2 - γ^2) for v in ν])
        @test drude.coefficients ≈ expected
        @test drude.rates ≈ [γ; ν]
        @test drude.counterterm == λ
        @test drude.hbar == ħ
        @test drude.diffusion > 0
        # The delta-correlated remainder restores the exact integrated real C(t).
        @test sum(real.(drude.coefficients) ./ drude.rates) + drude.diffusion ≈ 2λ * kT / γ
        bare = drude_lorentz_bath(;
            reorganization = λ,
            cutoff = γ,
            kT,
            hbar = ħ,
            matsubara = terms,
            terminator = false,
        )
        @test bare.coefficients == drude.coefficients
        @test bare.diffusion == 0
    end
    cold = drude_lorentz_bath(; reorganization = λ, cutoff = γ, kT, matsubara = 0)
    refined = drude_lorentz_bath(; reorganization = λ, cutoff = γ, kT, matsubara = 4)
    @test refined.diffusion < cold.diffusion
    absent = drude_lorentz_bath(; reorganization = 0, cutoff = γ, kT)
    @test all(iszero, absent.coefficients)
    @test iszero(absent.diffusion)
    @test iszero(absent.counterterm)

    @test_throws DimensionMismatch ExponentialBath([1, 2], [1])
    @test_throws ArgumentError ExponentialBath([1], [1 + im])
    for invalid in (Inf, -Inf, NaN)
        @test_throws ArgumentError ExponentialBath([invalid], [1])
        @test_throws ArgumentError ExponentialBath([complex(0, invalid)], [1])
    end
    for invalid in (0, -1, Inf, -Inf, NaN)
        @test_throws ArgumentError ExponentialBath([1], [invalid])
        @test_throws ArgumentError ExponentialBath([1], [1]; hbar = invalid)
        @test_throws ArgumentError drude_lorentz_bath(;
            reorganization = λ,
            cutoff = invalid,
            kT,
        )
        @test_throws ArgumentError drude_lorentz_bath(;
            reorganization = λ,
            cutoff = γ,
            kT = invalid,
        )
        @test_throws ArgumentError drude_lorentz_bath(;
            reorganization = λ,
            cutoff = γ,
            kT,
            hbar = invalid,
        )
    end
    for invalid in (-1, Inf, -Inf, NaN)
        @test_throws ArgumentError ExponentialBath([1], [1]; diffusion = invalid)
        @test_throws ArgumentError ExponentialBath([1], [1]; counterterm = invalid)
        @test_throws ArgumentError drude_lorentz_bath(;
            reorganization = invalid,
            cutoff = γ,
            kT,
        )
    end
    @test_throws ArgumentError drude_lorentz_bath(;
        reorganization = λ,
        cutoff = γ,
        kT,
        matsubara = -1,
    )
    # Coincident Drude and Matsubara poles need a different exponential expansion.
    @test_throws ArgumentError drude_lorentz_bath(;
        reorganization = λ,
        cutoff = 2π * kT / ħ,
        kT,
        hbar = ħ,
        matsubara = 1,
    )
end

@testset "Drude tail edge cases" begin
    λ, γ, kT, ħ = 0.3, 0.8, 1.2, 0.7
    @test_throws ArgumentError drude_lorentz_bath(;
        reorganization = λ,
        cutoff = 1e300,
        kT,
        hbar = 1e300,
    )
    @test_throws ArgumentError drude_lorentz_bath(; reorganization = λ, cutoff = 8, kT = 1)
    high_cutoff =
        drude_lorentz_bath(; reorganization = λ, cutoff = 8, kT = 1, terminator = false)
    @test iszero(high_cutoff.diffusion)
    # The classical high-temperature limit uses a stable series for the small tail.
    hot = drude_lorentz_bath(; reorganization = λ, cutoff = γ, kT = 1e8, hbar = ħ)
    @test hot.diffusion ≈ λ * ħ^2 * γ / (6e8) rtol = 1e-14
    # Large Matsubara rates remain finite even when their squares overflow Float64.
    very_hot = drude_lorentz_bath(;
        reorganization = λ,
        cutoff = γ,
        kT = 1e160,
        hbar = ħ,
        matsubara = 1,
    )
    @test real(very_hot.coefficients[2]) ≈ 2λ * γ * ħ / π
    @test very_hot.diffusion > 0
    @test_throws ArgumentError drude_lorentz_bath(;
        reorganization = λ,
        cutoff = 1e-300,
        kT,
        hbar = 1e-300,
    )
end

@testset "Hierarchy indexing and problem construction" begin
    grid = PhaseSpaceGrid((-7, 7), 32, (-6, 6), 40)
    V(q) = q^2 / 2
    for modes in 0:3, depth in 0:3
        bath = ExponentialBath(fill(0.2 - 0.1im, modes), collect(1:modes))
        op = heom_operator(grid; mass = 1, potential = V, bath, depth)
        indices = hierarchy_indices(op)
        @test length(indices) == binomial(modes + depth, depth)
        @test first(indices) == zeros(Int, modes)
        @test length(unique(Tuple.(indices))) == length(indices)
        @test all(n -> length(n) == modes && all(>=(0), n) && sum(n) <= depth, indices)
        @test size(op.upper) == (modes, length(indices))
        @test size(op.lower) == (modes, length(indices))
        for (i, n) in enumerate(indices), k in 1:modes
            up, down = op.upper[k, i], op.lower[k, i]
            if sum(n) == depth
                @test iszero(up)
            else
                expected = copy(n)
                expected[k] += 1
                @test indices[up] == expected
                @test op.lower[k, up] == i
            end
            if iszero(n[k])
                @test iszero(down)
            else
                expected = copy(n)
                expected[k] -= 1
                @test indices[down] == expected
                @test op.upper[k, down] == i
            end
        end
        if modes > 0
            indices[1][1] = 100
            @test first(hierarchy_indices(op)) == zeros(Int, modes)
        end
    end

    bath = ExponentialBath([0.3 - 0.1im], [0.7]; hbar = 0.6)
    op = heom_operator(grid; mass = 1.3, potential = V, bath, depth = 2)
    W0 = on_grid((q, p) -> exp(-q^2 - p^2), grid)
    prob = heom_problem(W0, (0.0, 1.0), op)
    @test prob.p === op
    @test prob.u0 isa Array{Float64,3}
    @test size(prob.u0) == (size(grid)..., 3)
    @test physical_wigner(prob.u0) == W0
    @test all(iszero, prob.u0[:, :, 2:end])
    @test prob.tspan == (0.0, 1.0)
    @test parent(physical_wigner(prob.u0)) === prob.u0
    prob.u0[1, 1, 1] += 1
    @test prob.u0[1, 1, 1] != W0[1, 1]
    U0 = cat(W0, 0.2W0, -0.1W0; dims = 3)
    full = heom_problem(U0, (0.0, 1.0), op)
    @test full.u0 == U0
    @test full.u0 !== U0
    built = heom_problem(W0, (0.0, 1.0), grid; mass = 1.3, potential = V, bath, depth = 2)
    @test physical_wigner(built.u0) == W0
    built_full =
        heom_problem(U0, (0.0, 1.0), grid; mass = 1.3, potential = V, bath, depth = 2)
    @test built_full.u0 == U0
    @test occursin("HEOM", repr(op))
    bath.coefficients[1] = 99
    bath.rates[1] = 99
    @test op.bath.coefficients == [0.3 - 0.1im]
    @test op.bath.rates == [0.7]
    @test op.coordinate_coupling ≈ [-1 / 3]

    @test_throws ArgumentError heom_operator(
        grid;
        mass = 1,
        potential = V,
        bath,
        depth = -1,
    )
    for invalid in (0, -1, Inf, -Inf, NaN)
        @test_throws ArgumentError heom_operator(
            grid;
            mass = invalid,
            potential = V,
            bath,
            depth = 1,
        )
    end
    @test_throws ArgumentError heom_operator(
        grid;
        mass = 1,
        potential = V,
        bath,
        depth = 1,
        discretization = FiniteDifference(),
    )
    @test_throws DimensionMismatch heom_problem(zeros(2, 2), (0.0, 1.0), op)
    @test_throws DimensionMismatch heom_problem(zeros(32, 40, 2), (0.0, 1.0), op)
    @test_throws DimensionMismatch heom_problem(zeros(31, 40, 3), (0.0, 1.0), op)
    @test_throws DimensionMismatch heom!(zeros(32, 40, 2), U0, op, 0.0)
    @test_throws DimensionMismatch heom!(similar(U0), zeros(32, 40, 2), op, 0.0)
    @test_throws ArgumentError physical_wigner(zeros(32, 40, 0))
    @test_throws ArgumentError heom_problem(complex.(U0), (0.0, 1.0), op)
    @test_throws DimensionMismatch heom_problem(zeros(32, 40, 3, 1), (0.0, 1.0), op)
    @test_throws ArgumentError heom_problem(fill(big"1e400", size(W0)), (0.0, 1.0), op)
    rapid = ExponentialBath([1], [1e308])
    @test_throws ArgumentError heom_operator(
        grid;
        mass = 1,
        potential = V,
        bath = rapid,
        depth = 2,
    )
    tiny_hbar = ExponentialBath([1im], [1]; hbar = 1e-309)
    @test_throws ArgumentError heom_operator(
        grid;
        mass = 1,
        potential = q -> zero(q),
        bath = tiny_hbar,
        depth = 1,
    )
    huge_real = ExponentialBath([1e308], [1])
    @test_throws ArgumentError heom_operator(
        grid;
        mass = 1,
        potential = q -> zero(q),
        bath = huge_real,
        depth = 2,
    )
    huge_imaginary = ExponentialBath([4e307im], [1])
    @test_throws ArgumentError heom_operator(
        grid;
        mass = 1,
        potential = q -> zero(q),
        bath = huge_imaginary,
        depth = 3,
    )
    for invalid in (NaN, Inf, -Inf)
        bad_W = copy(W0)
        bad_W[1] = invalid
        @test_throws ArgumentError heom_problem(bad_W, (0.0, 1.0), op)
        bad_U = copy(U0)
        bad_U[1, 1, 2] = invalid
        @test_throws ArgumentError heom_problem(bad_U, (0.0, 1.0), op)
    end
end

@testset "Wigner HEOM: analytical Gaussian hierarchy right-hand side" begin
    m, ħ, a, b, q0, p0 = 1.3, 0.7, 1.0, 0.8, 0.4, -0.3
    λ, μ, ε = 0.1, 0.5, 0.1
    V(q) = λ * q^4 - μ * q^2 + ε * q
    bath = ExponentialBath(
        [0.4 - 0.2im, -0.1 + 0.03im],
        [0.6, 1.7];
        hbar = ħ,
        diffusion = 0.07,
        counterterm = 0.12,
    )
    G(q, p) = exp(-a * (q - q0)^2 - b * (p - p0)^2)

    function gaussian_hierarchy(grid; kwargs...)
        op = heom_operator(grid; mass = m, potential = V, bath, depth = 2, kwargs...)
        indices = hierarchy_indices(op)
        weights = [(-1.0)^i / (i + 1) for i in eachindex(indices)]
        U = cat([weights[i] * on_grid(G, grid) for i in eachindex(indices)]...; dims = 3)
        exact = similar(U)
        for (i, n) in enumerate(indices)
            exact[:, :, i] = on_grid(grid) do q, p
                u = p - p0
                dp = -2b * u
                dpp = 4b^2 * u^2 - 2b
                dppp = 12b^2 * u - 8b^3 * u^3
                effective_force = 4λ * q^3 - 2μ * q + ε + 2bath.counterterm * q
                liouvillian =
                    2a * p / m * (q - q0) + effective_force * dp - ħ^2 * λ * q * dppp
                value =
                    weights[i] * (liouvillian - dot(n, bath.rates) + bath.diffusion * dpp)
                for k in eachindex(n)
                    upper = copy(n)
                    upper[k] += 1
                    j = findfirst(==(upper), indices)
                    if j !== nothing
                        value += weights[j] * dp
                    end
                    if n[k] > 0
                        lower = copy(n)
                        lower[k] -= 1
                        j = findfirst(==(lower), indices)
                        c = bath.coefficients[k]
                        value += n[k] * weights[j] * (real(c) * dp + 2imag(c) * q / ħ)
                    end
                end
                return value * G(q, p)
            end
        end
        return op, U, exact
    end

    grid = PhaseSpaceGrid((-7, 7), 64, (-7, 7), 72)
    for moyal_terms in (nothing, 2)
        op, U, exact = gaussian_hierarchy(grid; moyal_terms)
        dU = hierarchy_rhs(op, U)
        @test max_error(dU, exact) < 2e-11
        @test abs(sum(physical_wigner(dU))) < 1e-12 * sum(abs, physical_wigner(dU))
        @test (@inferred heom!(dU, U, op, 0.0)) === nothing
        @test probability_rate(U, op) ≈ 0 atol = 1e-12
        @test probability_rate(U, op; q = (0, Inf)) ≈
              probability(physical_wigner(dU), grid; q = (0, Inf)) atol = 1e-12
        @test expectation_rate((q, p) -> p^2, U, op) ≈
              expectation((q, p) -> p^2, physical_wigner(dU), grid) atol = 1e-12
        @test_throws ArgumentError probability_rate(physical_wigner(U), op)
        @test_throws ArgumentError expectation_rate((q, p) -> p, physical_wigner(U), op)
        @test U ==
              cat([(-1.0)^i / (i + 1) * on_grid(G, grid) for i in axes(U, 3)]...; dims = 3)
    end

    for order in (2, 4, 6)
        errors = map((64, 128)) do n
            fine = PhaseSpaceGrid((-7, 7), n, (-7, 7), n)
            op, U, exact = gaussian_hierarchy(
                fine;
                discretization = FiniteDifference(order),
                moyal_terms = 2,
            )
            max_error(hierarchy_rhs(op, U), exact)
        end
        @test log2(errors[1] / errors[2]) ≈ order atol = 0.5
    end

    # A zero bath leaves the isolated Wigner–Moyal dynamics and never populates ADOs.
    for coefficients in (ComplexF64[], ComplexF64[0, 0])
        empty_bath = ExponentialBath(coefficients, ones(length(coefficients)); hbar = ħ)
        for (discretization, moyal_terms) in
            ((Spectral(), nothing), (FiniteDifference(4), 2))
            isolated = wigner_moyal_operator(
                grid;
                mass = m,
                potential = V,
                hbar = ħ,
                discretization,
                moyal_terms,
            )
            prob = heom_problem(
                on_grid(G, grid),
                (0.0, 1.0),
                grid;
                mass = m,
                potential = V,
                bath = empty_bath,
                depth = 2,
                discretization,
                moyal_terms,
            )
            dU = hierarchy_rhs(prob.p, prob.u0)
            @test physical_wigner(dU) ≈ rhs(isolated, physical_wigner(prob.u0)) atol = 1e-12
            @test all(iszero, dU[:, :, 2:end])
        end
    end

    # Momentum diffusion must damp the Nyquist mode even though its first derivative is zero.
    nyquist_grid = PhaseSpaceGrid((-4, 4), 16, (-4, 4), 24)
    diffusion = 0.17
    diffusion_bath = ExponentialBath(Float64[], Float64[]; diffusion)
    diffusion_op = heom_operator(
        nyquist_grid;
        mass = 1,
        potential = q -> zero(q),
        bath = diffusion_bath,
        depth = 0,
    )
    alternating = [
        (-1.0)^j for
        i in eachindex(nyquist_grid.q), j in eachindex(nyquist_grid.p), k in 1:1
    ]
    @test hierarchy_rhs(diffusion_op, alternating) ≈
          -diffusion * (π / nyquist_grid.dp)^2 .* alternating atol = 1e-11
end

@testset "Coupled harmonic HEOM: independent closed moment evolution" begin
    m, ω, ħ = 1.3, 0.9, 0.7
    c, γ, D, λ = 0.35 - 0.12im, 1.1, 0.03, 0.15
    bath = ExponentialBath([c], [γ]; hbar = ħ, diffusion = D, counterterm = λ)
    V = harmonic_potential(; mass = m, omega = ω)
    grid = PhaseSpaceGrid((-9, 9), 64, (-9, 9), 64)
    μ0, Σ0 = [0.8, -0.4], [0.8 0.08; 0.08 0.9]
    W0 = on_grid((q, p) -> gaussian_wigner(q, p; mean = μ0, covariance = Σ0), grid)
    times = [0.0, 0.3, 0.7, 1.0]
    K, α = m * ω^2 + 2λ, 2imag(c) / ħ
    # The physical first moments close with the first ADO's integral.
    Amean = [0 1/m 0; -K 0 -1; α 0 -γ]
    # Raw second moments [<q²>, <qp>, <p²>, ∫qW₁, ∫pW₁, ∫W₂, 1].
    # Integration by parts derives these equations without discretising phase space.
    Asecond = zeros(7, 7)
    Asecond[1, 2] = 2 / m
    Asecond[2, [1, 3, 4]] = [-K, 1 / m, -1]
    Asecond[3, [2, 5, 7]] = [-2K, -2, 2D]
    Asecond[4, [1, 4, 5]] = [α, -γ, 1 / m]
    Asecond[5, [2, 4, 5, 6, 7]] = [α, -K, -γ, -1, -real(c)]
    Asecond[6, [4, 6]] = [2α, -2γ]
    mean_initial = [μ0; 0]
    second_initial =
        [Σ0[1, 1] + μ0[1]^2, Σ0[1, 2] + prod(μ0), Σ0[2, 2] + μ0[2]^2, 0, 0, 0, 1]

    prob =
        heom_problem(W0, (0.0, last(times)), grid; mass = m, potential = V, bath, depth = 2)
    sol = solve(prob, Vern7(); abstol = 1e-10, reltol = 1e-10, saveat = times)
    @test sol.t ≈ times
    @test physical_wigner(sol) == physical_wigner(last(sol.u))
    @test physical_wigner(sol, 1) == W0
    @test HEOM.plotting_state(sol) == (physical_wigner(sol), grid)
    @test HEOM.plotting_state(sol, 1) == (W0, grid)
    for (t, U) in zip(sol.t, sol.u)
        W = physical_wigner(U)
        μ = (exp(Amean*t)*mean_initial)[1:2]
        second = exp(Asecond * t) * second_initial
        Σ = [second[1] second[2]; second[2] second[3]] - μ * μ'
        @test phase_space_integral(W, grid) ≈ 1 atol = 1e-10
        @test phase_space_mean(W, grid) ≈ μ atol = 2e-8
        @test phase_space_covariance(W, grid) ≈ Σ atol = 2e-8
        @test phase_space_integral(U[:, :, 2], grid) ≈ (exp(Amean*t)*mean_initial)[3] atol =
            2e-8
        @test phase_space_integral(U[:, :, 3], grid) ≈ second[6] atol = 2e-8
    end
    # The hierarchy does transfer population to auxiliary states and alters motion.
    @test norm(sol.u[end][:, :, 2]) > 0.01
    isolated_mean = exp([0 1/m; -m*ω^2 0] * last(times)) * μ0
    @test norm(phase_space_mean(physical_wigner(sol), grid) - isolated_mean) > 0.01
    trajectory = diagnostics(sol; potential = V)
    @test trajectory.t == times
    @test trajectory.norm ≈ ones(length(times)) atol = 1e-10
    @test trajectory.mean_q[end] ≈ phase_space_mean(physical_wigner(sol), grid)[1] atol =
        1e-12
    @test trajectory.energy[end] ≈
          energy(physical_wigner(sol), grid; mass = m, potential = V) atol = 1e-12
    isolated =
        wigner_moyal_problem(W0, (0.0, 0.01), grid; mass = m, potential = V, hbar = ħ)
    isolated_sol = solve(isolated, Vern7())
    @test_throws ArgumentError physical_wigner(isolated_sol)
end
