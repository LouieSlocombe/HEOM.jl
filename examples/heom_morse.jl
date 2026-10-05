# Finite-time vibrational relaxation in a Morse well with a Drude–Lorentz bath.
# Requires HEOM, OrdinaryDiffEqVerner and Plots in the active Julia environment.
# Run: julia --project=/path/to/environment examples/heom_morse.jl [output_dir]
# On headless machines, set GKSwstype=100 before starting Julia.
using HEOM, OrdinaryDiffEqVerner, Plots, Printf

function morse_trajectory(; depth = 8, matsubara = 1, nq = 64, np = 48)
    mass, hbar = 1.0, 0.5
    dissociation, range_parameter = 8.0, 0.25
    potential(q) = dissociation * (1 - exp(-range_parameter * q))^2
    # The well minimum is q=0, V=0; V approaches dissociation as q→+∞.
    # Its curvature defines the initial Gaussian width, not a harmonic evolution.
    omega = range_parameter * sqrt(2dissociation / mass)
    grid = PhaseSpaceGrid((-4.0, 8.0), nq, (-5.0, 5.0), np)
    W0 = on_grid(
        (q, p) -> coherent_wigner(q, p; q0 = 1.5, p0 = 0.0, mass, omega, hbar),
        grid,
    )
    # This displaced Gaussian is not a Morse eigenstate. Its mean energy is
    # well below dissociation, but that does not make every component bound.
    # Weak coupling leaves anharmonic quantum structure visible over two periods.
    bath =
        drude_lorentz_bath(; reorganization = 0.03, cutoff = 1.2, kT = 0.5, hbar, matsubara)
    # The bath adds its counterterm internally. Matrix W0 initializes zero ADOs,
    # a factorized bare-bath preparation that includes the initial slip.
    # Unlike a harmonic well, Morse requires more than one Moyal term. `nothing`
    # retains the full spectral symbol, with no Taylor truncation of the potential.
    op = heom_operator(
        grid;
        mass,
        potential,
        bath,
        depth,
        scaled = true,
        moyal_terms = nothing,
    )
    println(
        "Solving Morse HEOM: depth=",
        depth,
        ", Matsubara=",
        matsubara,
        ", ADOs=",
        length(hierarchy_indices(op)),
    )
    sol = solve(
        heom_problem(W0, (0.0, 12.0), op),
        Vern7();
        saveat = 0.1,
        abstol = 1e-9,
        reltol = 1e-9,
    )
    @assert last(sol.t) == 12.0 "The solver did not reach the final time"
    @assert all(U -> all(isfinite, U), sol.u)
    return (; sol, grid, potential, mass, hbar, dissociation, range_parameter)
end

function check_morse_trajectory(data, coarse)
    (; sol, potential) = data
    d = diagnostics(sol; potential)
    @assert sol.t ≈ coarse.sol.t
    depth_change = maximum(
        maximum(abs, physical_wigner(U) - physical_wigner(Ucoarse)) for
        (U, Ucoarse) in zip(sol.u, coarse.sol.u)
    )
    checks = (;
        depth_6_to_8 = depth_change,
        norm_drift = maximum(abs, d.norm .- first(d.norm)),
        boundary_weight = max(maximum(d.boundary_q), maximum(d.boundary_p)),
        spectral_tail = max(maximum(d.tail_q), maximum(d.tail_p)),
    )
    for (name, value) in pairs(checks)
        println(name, " = ", value)
    end
    @assert checks.depth_6_to_8 < 1e-4
    @assert checks.norm_drift < 1e-8
    @assert checks.boundary_weight < 1e-4
    @assert checks.spectral_tail < 1e-3
    # These checks support this finite-time example, not arbitrary parameters.
    # Vary matsubara, nq/np and box size separately for quantitative convergence.
    # The exact symbol evaluates V beyond the grid by π*hbar/(2dp): a finer p grid
    # or a box extending further left makes the exponential wall much stiffer.
    # An unconfined Morse well has a continuum, so do not interpret this short
    # trajectory as a globally normalizable thermal equilibrium or a rate study.
    return d, checks
end

function render_morse_trajectory(data, d, checks, output_dir)
    (; sol, grid, potential, dissociation, range_parameter) = data
    mkpath(output_dir)
    qview = (-2.5, 4.5)
    W0 = physical_wigner(sol, 1)
    color_limit = maximum(maximum(abs, physical_wigner(U)) for U in sol.u)
    density_limit =
        1.1 * maximum(maximum(position_density(physical_wigner(U), grid)) for U in sol.u)
    animation = Animation()
    for i in eachindex(sol.t)
        W = physical_wigner(sol, i)
        phase = wignerplot(
            W,
            grid;
            xlims = qview,
            ylims = (-3.5, 3.5),
            clims = (-color_limit, color_limit),
            title = "Wigner distribution",
            xlabel = "Position q",
            ylabel = "Momentum p",
        )
        plot!(phase, d.mean_q[1:i], d.mean_p[1:i]; color = :black, lw = 1, label = "")
        scatter!(phase, [d.mean_q[i]], [d.mean_p[i]]; color = :orange, ms = 4, label = "")
        density = plot(
            grid.q,
            position_density(W0, grid);
            color = :gray,
            linestyle = :dash,
            lw = 2,
            label = "Initial",
            xlabel = "Position q",
            ylabel = "P(q)",
            title = "Position probability",
            xlims = qview,
            ylims = (0, density_limit),
        )
        plot!(
            density,
            grid.q,
            position_density(W, grid);
            lw = 2.5,
            color = "#176b9c",
            fillrange = 0,
            fillalpha = 0.12,
            label = "Current",
        )
        q = range(qview...; length = 200)
        well = plot(
            q,
            potential.(q);
            color = "#7353a6",
            lw = 2.5,
            label = "Morse V(q)",
            title = "Asymmetric well",
            xlabel = "Position q",
            ylabel = "Energy",
            xlims = qview,
            ylims = (0, 3.5),
        )
        hline!(well, [d.energy[i]]; color = :gray, linestyle = :dash, label = "Mean energy")
        scatter!(
            well,
            [d.mean_q[i]],
            [potential(d.mean_q[i])];
            color = :orange,
            ms = 5,
            label = "V(mean q)",
        )
        fig = plot(
            phase,
            density,
            well;
            layout = (1, 3),
            size = (1320, 460),
            plot_title = @sprintf(
                "HEOM | Morse oscillator | D = %.0f, a = %.2f | t = %.1f",
                dissociation,
                range_parameter,
                sol.t[i]
            ),
            margin = 5Plots.mm,
            titlefontsize = 11,
            guidefontsize = 10,
            tickfontsize = 9,
            legendfontsize = 8,
            background_color = "#fcfcfd",
        )
        frame(animation, fig)
        if i in (1, 31, length(sol.t))
            savefig(fig, joinpath(output_dir, "morse_frame_$(i).png"))
        end
    end
    gif(animation, joinpath(output_dir, "heom_morse.gif"); fps = 15)
    mp4(animation, joinpath(output_dir, "heom_morse.mp4"); fps = 15)
    savefig(
        diagnosticsplot(d; fields = (:mean_q, :mean_p, :energy, :purity)),
        joinpath(output_dir, "morse_diagnostics.png"),
    )
    open(joinpath(output_dir, "observables.csv"), "w") do io
        fields = (:t, :mean_q, :mean_p, :energy, :purity, :norm, :negativity)
        println(io, join(string.(fields), ","))
        for i in eachindex(sol.t)
            println(io, join((getproperty(d, f)[i] for f in fields), ","))
        end
    end
    open(joinpath(output_dir, "validation.txt"), "w") do io
        println(io, "Morse: D=8, a=0.25, mass=1, hbar=0.5; initial q=1.5, p=0")
        println(io, "Drude bath: lambda=0.03, cutoff=1.2, kT=0.5, Matsubara=1")
        println(io, "Depth=8 (45 ADOs); full Moyal symbol; 64 x 48 grid")
        println(io, "Grid: q in [-4, 8), p in [-5, 5); time in [0, 12]")
        for (name, value) in pairs(checks)
            println(io, name, " = ", value)
        end
        println(io, "Initial / final energy: ", first(d.energy), " / ", last(d.energy))
        println(
            io,
            "Depth comparison and grid diagnostics do not establish bath-pole convergence.",
        )
    end
    println("Morse animations and diagnostics saved to ", output_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_dir = isempty(ARGS) ? mktempdir(; prefix = "heom-morse-") : abspath(only(ARGS))
    coarse = morse_trajectory(; depth = 6)
    data = morse_trajectory()
    d, checks = check_morse_trajectory(data, coarse)
    render_morse_trajectory(data, d, checks, output_dir)
end
