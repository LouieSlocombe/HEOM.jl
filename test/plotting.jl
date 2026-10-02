ENV["GKSwstype"] = "100"
using Plots
using RecipesBase

function plotting_recipe(kind, args...; attributes...)
    return RecipesBase.apply_recipe(Dict{Symbol,Any}(attributes), kind(args))
end

@testset "Plotting recipes" begin
    # Different axis lengths, spacings and values expose transposition mistakes that
    # square grids and rotationally symmetric Wigner functions would conceal.
    grid = PhaseSpaceGrid((-2, 4), 4, (-5, 10), 5)
    W = on_grid((q, p) -> 2q - 3p, grid)

    @testset "Signed Wigner heatmap" begin
        series = only(plotting_recipe(HEOM.WignerPlot, W, grid))
        q, p, z = series.args
        @test q == grid.q
        @test p == grid.p
        @test size(z) == reverse(size(grid))
        @test z[4, 2] == W[2, 4]
        @test z == permutedims(W)
        @test minimum(z) < 0 < maximum(z)
        @test series.plotattributes[:seriestype] == :heatmap
        @test series.plotattributes[:xguide] == "q"
        @test series.plotattributes[:yguide] == "p"
        @test series.plotattributes[:colorbar_title] == "W(q, p)"
        @test series.plotattributes[:colorbar] === true
        @test series.plotattributes[:seriescolor] == :RdBu
        @test series.plotattributes[:clims] == (-maximum(abs, W), maximum(abs, W))

        overridden = only(
            plotting_recipe(
                HEOM.WignerPlot,
                W,
                grid;
                clims = (-2, 3),
                xguide = "position",
                yguide = "momentum",
                seriescolor = :viridis,
                colorbar = false,
                title = "Signed state",
            ),
        )
        @test overridden.plotattributes[:clims] == (-2, 3)
        @test overridden.plotattributes[:xguide] == "position"
        @test overridden.plotattributes[:yguide] == "momentum"
        @test overridden.plotattributes[:seriescolor] == :viridis
        @test overridden.plotattributes[:colorbar] === false
        @test overridden.plotattributes[:title] == "Signed state"
        zero_series = only(plotting_recipe(HEOM.WignerPlot, zero(W), grid))
        @test zero_series.plotattributes[:clims] == (-1, 1)
        @test all(iszero, last(zero_series.args))
    end

    @testset "Marginals retain the grid quadrature" begin
        series = plotting_recipe(HEOM.MarginalPlot, W, grid)
        @test length(series) == 2
        @test series[1].args == (grid.q, position_density(W, grid))
        @test series[2].args == (grid.p, momentum_density(W, grid))
        @test series[1].plotattributes[:xguide] == "q"
        @test series[2].plotattributes[:xguide] == "p"
        @test series[1].plotattributes[:subplot] != series[2].plotattributes[:subplot]
        @test sum(last(series[1].args)) * grid.dq ≈ phase_space_integral(W, grid)
        @test sum(last(series[2].args)) * grid.dp ≈ phase_space_integral(W, grid)
        scaled = plotting_recipe(HEOM.MarginalPlot, 0.4W, grid)
        @test last(scaled[1].args) ≈ 0.4last(series[1].args)
        @test last(scaled[2].args) ≈ 0.4last(series[2].args)
        overridden = plotting_recipe(HEOM.MarginalPlot, W, grid; linewidth = 3)
        @test all(s -> s.plotattributes[:linewidth] == 3, overridden)
    end

    @testset "State validation" begin
        for kind in (HEOM.WignerPlot, HEOM.MarginalPlot)
            @test_throws DimensionMismatch plotting_recipe(kind, ones(2, 5), grid)
            @test_throws ArgumentError plotting_recipe(kind, complex.(W), grid)
            @test_throws ArgumentError plotting_recipe(kind, fill(NaN, size(grid)), grid)
            @test_throws ArgumentError plotting_recipe(kind, fill(Inf, size(grid)), grid)
            @test_throws MethodError plotting_recipe(kind, W)
            @test_throws MethodError plotting_recipe(kind, W, grid, 1)
        end
    end

    @testset "Saved solution states" begin
        V = harmonic_potential(; mass = 1, omega = 1)
        prob = wigner_moyal_problem(W, (0.0, 0.01), grid; mass = 1, potential = V)
        sol = solve(prob, Vern9(); saveat = [0.0, 0.01], abstol = 1e-10, reltol = 1e-10)
        for kind in (HEOM.WignerPlot, HEOM.MarginalPlot)
            final = plotting_recipe(kind, sol)
            expected_final = plotting_recipe(kind, sol.u[end], grid)
            initial = plotting_recipe(kind, sol, 1)
            expected_initial = plotting_recipe(kind, sol.u[1], grid)
            @test map(s -> s.args, final) == map(s -> s.args, expected_final)
            @test map(s -> s.args, initial) == map(s -> s.args, expected_initial)
            @test_throws BoundsError plotting_recipe(kind, sol, 0)
            @test_throws BoundsError plotting_recipe(kind, sol, length(sol.u) + 1)
            @test_throws MethodError plotting_recipe(kind, sol, 1.5)
        end
        other_prob = ODEProblem((du, u, p, t) -> fill!(du, 0.0), W, (0.0, 0.01))
        other_sol = solve(other_prob, Vern9(); save_everystep = false)
        @test_throws ArgumentError plotting_recipe(HEOM.WignerPlot, other_sol)
        @test_throws ArgumentError plotting_recipe(HEOM.MarginalPlot, other_sol)
        @test_throws ArgumentError plotting_recipe(
            HEOM.DiagnosticsPlot,
            other_sol;
            potential = V,
        )
    end

    @testset "Caldeira–Leggett solution integration" begin
        V = harmonic_potential(; mass = 1, omega = 1)
        damped_grid = PhaseSpaceGrid((-6, 6), 24, (-7, 7), 28)
        W0 = on_grid(
            (q, p) -> coherent_wigner(q, p; q0 = 0.6, mass = 1, omega = 1),
            damped_grid,
        )
        prob = caldeira_leggett_problem(
            W0,
            (0.0, 0.01),
            damped_grid;
            mass = 1,
            potential = V,
            friction = 0.5,
            kT = 2,
        )
        sol = solve(prob, Vern9(); save_everystep = false)
        for kind in (HEOM.WignerPlot, HEOM.MarginalPlot)
            actual = plotting_recipe(kind, sol)
            expected = plotting_recipe(kind, sol.u[end], damped_grid)
            @test map(s -> s.args, actual) == map(s -> s.args, expected)
        end
        actual = plotting_recipe(HEOM.DiagnosticsPlot, sol; potential = V)
        expected = plotting_recipe(HEOM.DiagnosticsPlot, diagnostics(sol; potential = V))
        @test map(s -> s.args, actual) == map(s -> s.args, expected)
        @test Base.get_extension(HEOM, :HEOMPlotsExt).animation_operator(sol) === prob.p
    end

    @testset "Diagnostic trajectories" begin
        times = [0.0, 0.2, 0.7]
        values = (
            norm = [0.8, 0.79, 0.77],
            energy = [2.0, 1.9, 1.8],
            purity = [0.64, 0.62, 0.59],
            negativity = [0.0, 0.01, 0.02],
            boundary_q = [0.0, NaN, 1e-5],
        )
        d = merge((; t = times), values)
        fields = (:norm, :energy, :purity, :negativity)
        series = plotting_recipe(HEOM.DiagnosticsPlot, d)
        @test length(series) == length(fields)
        for (i, field) in enumerate(fields)
            @test series[i].args == (times, values[field])
            @test series[i].plotattributes[:subplot] == i
            @test series[i].plotattributes[:yguide] == string(field)
        end
        explicit_times = plotting_recipe(HEOM.DiagnosticsPlot, times, values)
        @test map(s -> s.args, explicit_times) == map(s -> s.args, series)
        chosen = plotting_recipe(HEOM.DiagnosticsPlot, d; fields = [:purity, :norm])
        @test map(s -> last(s.args), chosen) == [values.purity, values.norm]
        gaps = only(plotting_recipe(HEOM.DiagnosticsPlot, d; fields = :boundary_q))
        @test isequal(last(gaps.args), values.boundary_q)
        overridden = only(
            plotting_recipe(
                HEOM.DiagnosticsPlot,
                d;
                fields = :norm,
                xguide = "time / fs",
                linewidth = 4,
            ),
        )
        @test overridden.plotattributes[:xguide] == "time / fs"
        @test overridden.plotattributes[:linewidth] == 4

        @test_throws ArgumentError plotting_recipe(HEOM.DiagnosticsPlot, values)
        @test_throws ArgumentError plotting_recipe(HEOM.DiagnosticsPlot, Float64[], values)
        @test_throws ArgumentError plotting_recipe(HEOM.DiagnosticsPlot, d; fields = ())
        @test_throws ArgumentError plotting_recipe(
            HEOM.DiagnosticsPlot,
            d;
            fields = (:norm, :norm),
        )
        @test_throws ArgumentError plotting_recipe(HEOM.DiagnosticsPlot, d; fields = :t)
        @test_throws ArgumentError plotting_recipe(
            HEOM.DiagnosticsPlot,
            d;
            fields = :absent,
        )
        @test_throws ArgumentError plotting_recipe(HEOM.DiagnosticsPlot, d; fields = (1,))
        @test_throws ArgumentError plotting_recipe(
            HEOM.DiagnosticsPlot,
            times,
            (norm = 1.0,);
            fields = :norm,
        )
        @test_throws DimensionMismatch plotting_recipe(
            HEOM.DiagnosticsPlot,
            times,
            (norm = [1.0, 2.0],);
            fields = :norm,
        )

        physical_grid = PhaseSpaceGrid((-6, 6), 24, (-7, 7), 28)
        W0 = on_grid(
            (q, p) -> coherent_wigner(q, p; q0 = 0.6, p0 = -0.2, mass = 1, omega = 1),
            physical_grid,
        )
        V = harmonic_potential(; mass = 1, omega = 1)
        prob = wigner_moyal_problem(W0, (0.0, 0.02), physical_grid; mass = 1, potential = V)
        sol = solve(prob, Vern9(); saveat = [0.0, 0.02], abstol = 1e-10, reltol = 1e-10)
        from_solution = plotting_recipe(HEOM.DiagnosticsPlot, sol; potential = V)
        expected = plotting_recipe(HEOM.DiagnosticsPlot, diagnostics(sol; potential = V))
        @test map(s -> s.args, from_solution) == map(s -> s.args, expected)
        @test_throws ArgumentError plotting_recipe(HEOM.DiagnosticsPlot, sol)
    end

    @testset "Plots integration and rendering" begin
        heatmap = wignerplot(W, grid; xlabel = "position", ylabel = "momentum")
        @test heatmap isa Plots.Plot
        @test heatmap.series_list[1][:x] == grid.q
        @test heatmap.series_list[1][:y] == grid.p
        @test heatmap.series_list[1][:z].surf == permutedims(W)
        @test heatmap.subplots[1][:xaxis][:guide] == "position"
        @test heatmap.subplots[1][:yaxis][:guide] == "momentum"
        @test heatmap.subplots[1][:colorbar] != :none
        hidden_colorbar = wignerplot(W, grid; colorbar = false)
        @test hidden_colorbar.subplots[1][:colorbar] == :none
        marginals = marginalplot(
            W,
            grid;
            layout = (2, 1),
            xlabel = "coordinate",
            ylabel = "density",
        )
        @test length(marginals.subplots) == 2
        @test size(marginals.layout) == (2, 1)
        @test marginals.series_list[1][:y] == position_density(W, grid)
        @test marginals.series_list[2][:y] == momentum_density(W, grid)
        @test all(sp -> sp[:xaxis][:guide] == "coordinate", marginals.subplots)
        @test all(sp -> sp[:yaxis][:guide] == "density", marginals.subplots)
        d = (t = [0.0, 0.5], norm = [0.9, 0.85], purity = [0.8, 0.7])
        trajectory = diagnosticsplot(d; fields = (:norm, :purity), xlabel = "time / fs")
        @test length(trajectory.subplots) == 2
        @test trajectory.series_list[1][:x] == d.t
        @test trajectory.series_list[1][:y] == d.norm
        @test trajectory.series_list[2][:y] == d.purity
        three = diagnosticsplot(
            merge(d, (energy = [2.0, 1.9],));
            fields = (:norm, :purity, :energy),
        )
        @test length(three.subplots) == 3
        @test all(sp -> sp[:xaxis][:guide] == "time / fs", trajectory.subplots)

        # Both explicit-target and current-plot mutators must append the same data.
        for (figure, mutate, args, kwargs, count) in (
            (heatmap, wignerplot!, (W, grid), (;), 1),
            (marginals, marginalplot!, (W, grid), (;), 2),
            (trajectory, diagnosticsplot!, (d,), (; fields = (:norm, :purity)), 2),
        )
            @test mutate(figure, args...; kwargs...) === figure
            @test length(figure.series_list) == 2count
            @test mutate(args...; kwargs...) === figure
            @test length(figure.series_list) == 3count
        end

        mktempdir() do directory
            for (name, figure) in
                (("wigner", heatmap), ("marginals", marginals), ("diagnostics", trajectory))
                path = joinpath(directory, "$name.png")
                savefig(figure, path)
                @test read(path)[1:8] ==
                      UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
            end
        end
    end
end
