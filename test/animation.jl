@testset "Animation helpers" begin
    extension = Base.get_extension(HEOM, :HEOMPlotsExt)
    @test extension !== nothing

    # Unequal axes and signed, asymmetric states expose transposition errors and
    # distinguish each marginal's shared scale from the final frame's scale.
    grid = PhaseSpaceGrid((-2, 4), 4, (-5, 10), 5)
    W = on_grid((q, p) -> 2q - 3p, grid)
    states = [W, -2W, 100W]
    original_states = deepcopy(states)
    times = [0.0, 0.2, 3.0]
    selected = [2, 1]

    @testset "Selection and shared scales" begin
        @test extension.animation_data(states, grid) == [1, 2, 3]
        @test extension.animation_data(states, grid; indices = 3:-1:1) == [3, 2, 1]
        @test extension.animation_data(states, grid; indices = [2, 1, 2]) == [2, 1, 2]
        @test extension.animation_data(states, grid; times, indices = selected) == selected
        @test extension.wigner_limits(states, selected) == (-50.0, 50.0)
        @test extension.wigner_limits([zero(W)], [1]) == (-1.0, 1.0)

        limits = extension.marginal_limits(states, grid, selected)
        @test collect(limits[1]) ≈ [-120.75, 225.75]
        @test collect(limits[2]) ≈ [-207.6, 267.6]
        zeros_limits = extension.marginal_limits([zero(W)], grid, [1])
        @test all(==((-1.0, 1.0)), zeros_limits)

        # Positive densities must still include zero, with a little space below it.
        positive_limits = extension.marginal_limits([ones(size(grid))], grid, [1])
        @test collect(positive_limits[1]) ≈ [-0.75, 15.75]
        @test collect(positive_limits[2]) ≈ [-0.3, 6.3]

        # Unselected states need not be plotted or included in validation/scaling.
        partial = [W, fill(NaN, size(grid))]
        @test extension.animation_data(partial, grid; indices = [1]) == [1]
    end

    @testset "Input validation" begin
        for animate in (wigneranimation, marginalanimation)
            @test_throws ArgumentError animate(Matrix{Float64}[], grid)
            @test_throws ArgumentError animate(states, grid; indices = Int[])
            for indices in ([true], [1.5], ["1"])
                @test_throws ArgumentError animate(states, grid; indices)
            end
            for indices in ([0], [4])
                @test_throws BoundsError animate(states, grid; indices)
            end
            @test_throws DimensionMismatch animate([ones(2, 5)], grid)
            @test_throws ArgumentError animate([complex.(W)], grid)
            @test_throws ArgumentError animate([fill(NaN, size(grid))], grid)
            @test_throws ArgumentError animate([fill(Inf, size(grid))], grid)
            @test_throws DimensionMismatch animate(states, grid; times = [0.0])
            for bad_times in
                ((0.0, 0.2, 3.0), [0.0, NaN, 3.0], [0.0, Inf, 3.0], complex.(times))
                @test_throws ArgumentError animate(states, grid; times = bad_times)
            end
            # Time validation covers all saved states, including omitted frames.
            @test_throws ArgumentError animate(
                states,
                grid;
                times = [0.0, 0.2, NaN],
                indices = selected,
            )
        end
    end

    animations = Plots.Animation[]
    try
        @testset "Frames, data and plotting attributes" begin
            heatmap =
                wigneranimation(states, grid; times, indices = selected, size = (240, 200))
            push!(animations, heatmap)
            @test heatmap isa Plots.Animation
            @test length(heatmap.frames) == length(selected)
            figure = Plots.current()
            @test figure.series_list[1][:x] == grid.q
            @test figure.series_list[1][:y] == grid.p
            @test figure.series_list[1][:z].surf == permutedims(W)
            @test figure.subplots[1][:clims] == (-50.0, 50.0)
            # Grid coordinates are cell centres; include the complete outer cells.
            @test figure.subplots[1][:xaxis][:lims] == (-2.75, 3.25)
            @test figure.subplots[1][:yaxis][:lims] == (-6.5, 8.5)
            @test figure.subplots[1][:title] == "t = 0.0"

            marginals =
                marginalanimation(states, grid; indices = selected, size = (400, 180))
            push!(animations, marginals)
            @test marginals isa Plots.Animation
            @test length(marginals.frames) == length(selected)
            figure = Plots.current()
            @test figure.series_list[1][:x] == grid.q
            @test figure.series_list[2][:x] == grid.p
            @test figure.series_list[1][:y] == position_density(W, grid)
            @test figure.series_list[2][:y] == momentum_density(W, grid)
            @test collect(figure.subplots[1][:yaxis][:lims]) ≈ [-120.75, 225.75]
            @test collect(figure.subplots[2][:yaxis][:lims]) ≈ [-207.6, 267.6]
            @test all(sp -> sp[:title] == "Frame 1", figure.subplots)

            overridden_heatmap = wigneranimation(
                states,
                grid;
                indices = [1],
                times,
                clims = (-3, 4),
                title = "Custom Wigner title",
                xlims = (-10, 10),
                ylims = (-12, 12),
                size = (240, 200),
            )
            push!(animations, overridden_heatmap)
            figure = Plots.current()
            @test figure.subplots[1][:clims] == (-3, 4)
            @test figure.subplots[1][:title] == "Custom Wigner title"
            @test figure.subplots[1][:xaxis][:lims] == (-10, 10)
            @test figure.subplots[1][:yaxis][:lims] == (-12, 12)

            overridden_marginals = marginalanimation(
                states,
                grid;
                indices = [1],
                title = "Custom marginal title",
                ylims = (-9, 11),
                size = (400, 180),
            )
            push!(animations, overridden_marginals)
            figure = Plots.current()
            @test all(sp -> sp[:title] == "Custom marginal title", figure.subplots)
            @test all(sp -> sp[:yaxis][:lims] == (-9, 11), figure.subplots)
            @test states == original_states

            png_signature = UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
            for animation in (heatmap, marginals), filename in animation.frames
                @test read(joinpath(animation.dir, filename))[1:8] == png_signature
            end

            mktempdir() do directory
                path = joinpath(directory, "wigner.gif")
                gif(heatmap, path; fps = 2, show_msg = false)
                @test String(read(path)[1:6]) in ("GIF87a", "GIF89a")
                @test filesize(path) > 100
            end
        end

        @testset "Saved solution states and times" begin
            V = harmonic_potential(; mass = 1, omega = 1)
            prob = wigner_moyal_problem(W, (0.0, 0.01), grid; mass = 1, potential = V)
            sol = solve(prob, Vern9(); saveat = [0.0, 0.01], abstol = 1e-10, reltol = 1e-10)
            for animate in (wigneranimation, marginalanimation)
                animation = animate(sol; indices = [2, 1], size = (400, 200))
                push!(animations, animation)
                @test length(animation.frames) == 2
                figure = Plots.current()
                @test all(sp -> sp[:title] == "t = $(sol.t[1])", figure.subplots)
                if animate === wigneranimation
                    @test figure.series_list[1][:z].surf == permutedims(sol.u[1])
                else
                    @test figure.series_list[1][:y] == position_density(sol.u[1], grid)
                    @test figure.series_list[2][:y] == momentum_density(sol.u[1], grid)
                end
            end

            bath = ExponentialBath([0.2 - 0.1im], [1.0])
            hierarchy =
                heom_problem(W, (0.0, 0.01), grid; mass = 1, potential = V, bath, depth = 1)
            # Large auxiliaries must not influence physical plot limits or frames.
            hierarchy.u0[:, :, 2] .= 100W
            hsol = solve(hierarchy, Vern9(); save_everystep = false)
            for animate in (wigneranimation, marginalanimation)
                animation = animate(hsol; indices = [1], size = (400, 200))
                push!(animations, animation)
                @test length(animation.frames) == 1
                figure = Plots.current()
                if animate === wigneranimation
                    @test figure.series_list[1][:z].surf == permutedims(W)
                    @test figure.subplots[1][:clims] == (-25.0, 25.0)
                else
                    @test figure.series_list[1][:y] == position_density(W, grid)
                    @test figure.series_list[2][:y] == momentum_density(W, grid)
                end
            end

            other_prob = ODEProblem((du, u, p, t) -> fill!(du, 0.0), W, (0.0, 0.01))
            other_sol = solve(other_prob, Vern9(); save_everystep = false)
            @test_throws ArgumentError wigneranimation(other_sol)
            @test_throws ArgumentError marginalanimation(other_sol)
        end
    finally
        for animation in animations
            rm(animation.dir; recursive = true, force = true)
        end
    end
end
