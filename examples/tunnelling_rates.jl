# Low-temperature tunnelling in a bath-coupled double well.
# Requires HEOM, OrdinaryDiffEqVerner, LinearSolve, SciMLBase and Plots.
# Run: julia --project=/path/to/environment examples/tunnelling_rates.jl [output_dir]
# On headless machines, set GKSwstype=100 before starting Julia.
# The stationary preparation uses about 2 GB; the trajectory takes a few minutes.
# Refine the independent numerical
# controls below before interpreting this finite-time illustration quantitatively.
using HEOM, OrdinaryDiffEqVerner, Plots, Printf
using LinearAlgebra, LinearSolve
using SciMLBase: successful_retcode
include("double_well_equilibrium.jl")

function double_well_preparation(grid; mass, hbar, potential)
    states = eigenstates(grid; mass, potential, hbar, nstates = 2)
    psi0, psi1 = states.wavefunctions[:, 1], states.wavefunctions[:, 2]
    # Align <0|q|1> so the minus superposition is localised in the left well.
    transition = sum(psi0 .* grid.q .* psi1) * grid.dq
    psi1 = sign(transition) * psi1
    initial = wavefunction_wigner((psi0 - psi1) / sqrt(2), grid; hbar)
    @assert states.energies[2] < potential(0) "Both prepared levels must be below the barrier"
    @assert abs(phase_space_integral(initial, grid) - 1) < 1e-10
    return (; initial, energies = states.energies)
end

function double_well_trajectory(;
    points = 49,
    depth = 4,
    pade = 2,
    qextent = 5.0,
    pextent = 6.0,
    stop = 40.0,
    saveat = 0.5,
    kT = 0.1,
    reorganization = 0.01,
    cutoff = 1.0,
    preparation = :biased_equilibrium,
    bias = 0.002,
    abstol = 1e-8,
    reltol = 1e-8,
)
    mass, hbar = 1.0, 0.75
    # Minima at ±1.5, barrier height 1; kT/(hbar*omega_well) = 0.0707.
    # The old warm example also used a different hbar and bath coupling:
    # these examples are not a controlled temperature-only comparison.
    potential(q) = (q^2 / 2.25 - 1)^2
    grid = PhaseSpaceGrid((-qextent, qextent), points, (-pextent, pextent), points)
    energies = eigenstates(grid; mass, potential, hbar, nstates = 2).energies
    # Retain quantum thermal poles and vary pade independently of depth and grid.
    bath = drude_lorentz_pade_bath(; reorganization, cutoff, kT, hbar, pade)
    op = heom_operator(grid; mass, potential, bath, depth, scaled = true, moyal_terms = 2)
    # The quartic Moyal series terminates at hbar²; HEOM adds the counterterm.
    equilibrium = nothing
    if preparation === :biased_equilibrium
        isfinite(bias) && bias > 0 ||
            throw(ArgumentError("bias must be finite and positive"))
        # A positive tilt favours the left well. Prepare the entire correlated
        # hierarchy, then quench only the potential back to the symmetric one.
        tilted(q) = potential(q) + bias * q
        preparation_op = heom_operator(
            grid;
            mass,
            potential = tilted,
            bath,
            depth,
            scaled = true,
            moyal_terms = 2,
        )
        seed = thermal_wigner(grid; mass, potential = tilted, kT, hbar)
        println("Preparing correlated equilibrium with tilt +", bias, " q")
        flush(stdout)
        equilibrium = double_well_equilibrium(seed, preparation_op)
        initial = equilibrium.hierarchy
        @assert hierarchy_indices(preparation_op) == hierarchy_indices(op)
        depletion = 0.5 - probability(physical_wigner(initial), grid; q = (0.0, Inf))
        depletion > 1e-6 || throw(
            ArgumentError(
                "bias gives no resolvable right-well depletion (must exceed 1e-6)",
            ),
        )
        @printf(
            "Preparation: P(right) = %.9f; maximum full RHS = %.3e\n",
            probability(physical_wigner(initial), grid; q = (0.0, Inf)),
            maximum(equilibrium.residuals.absolute)
        )
    elseif preparation === :localized_doublet
        initial = double_well_preparation(grid; mass, hbar, potential).initial
    else
        throw(
            ArgumentError("preparation must be :biased_equilibrium or :localized_doublet"),
        )
    end
    println("Solving cold double well: ", op)
    flush(stdout)
    sol = solve(
        heom_problem(initial, (0.0, stop), op),
        Vern7();
        saveat,
        abstol,
        reltol,
        dense = false,
    )
    @assert successful_retcode(sol) && last(sol.t) == stop "Integration did not complete"
    @assert all(U -> all(isfinite, U), sol.u)
    return (;
        sol,
        times = sol.t,
        states = [copy(physical_wigner(U)) for U in sol.u],
        grid,
        potential,
        mass,
        hbar,
        kT,
        reorganization,
        cutoff,
        depth,
        pade,
        energies,
        preparation,
        bias,
        equilibrium,
    )
end

function report_transfer_rates(data; io = stdout, deviation_floor = 1e-6)
    (; times, states, grid, potential, mass, hbar, kT, energies) = data
    population = [probability(W, grid; q = (0.0, Inf)) for W in states]
    isfinite(deviation_floor) && deviation_floor > 0 ||
        throw(ArgumentError("deviation_floor must be finite and positive"))
    abs(first(population) - 0.5) > deviation_floor ||
        throw(ArgumentError("initial depletion is below the reporting noise floor"))
    d = diagnostics(states, grid; mass, hbar, potential, times)
    splitting = energies[2] - energies[1]
    transfer_time = π * hbar / splitting
    # Exact for the isolated localised doublet. For the thermal preparation this
    # is only a frequency reference: it omits higher levels and bath correlations.
    isolated = 0.5 .+ (first(population) - 0.5) .* cos.(splitting .* times ./ hbar)
    peak = argmax(population)
    @printf(io, "kT = %.3f; barrier = 1; hbar = %.3f\n", kT, hbar)
    @printf(io, "kT/(hbar*omega_well) = %.6f\n", kT / (hbar * sqrt(8 / 2.25)))
    @printf(io, "Sub-barrier doublet energies: %.9f, %.9f\n", energies...)
    @printf(
        io,
        "Isolated splitting/hbar: %.9f (angular frequency, not a rate)\n",
        splitting / hbar
    )
    @printf(io, "Isolated first transfer time: %.6f\n", transfer_time)
    @printf(
        io,
        "Right-well population: initial %.6f; peak %.6f at t=%.2f; final %.6f\n",
        first(population),
        population[peak],
        times[peak],
        last(population)
    )
    near_equilibrium = data.preparation === :biased_equilibrium
    if near_equilibrium
        @printf(
            io,
            "Initial state: correlated equilibrium of V(q) + %.6g q; tilt removed at t=0.\n",
            data.bias
        )
        @printf(
            io,
            "Initial right-well depletion: %.6f percentage points\n",
            100 * (0.5 - first(population))
        )
        @printf(
            io,
            "Preparation maximum full-hierarchy RHS: %.3e\n",
            maximum(data.equilibrium.residuals.absolute)
        )
        println(
            io,
            "The isolated-doublet curve is a frequency reference, not an exact thermal trajectory.",
        )
    else
        println(
            io,
            "Initial system state: pure, left-localised doublet; kT describes the bath.",
        )
    end
    # Symmetry fixes equilibrium at one half. Do not manufacture a rate by
    # fitting a fraction of a coherent oscillation which crosses equilibrium.
    deviations = filter(x -> abs(x) > deviation_floor, population .- 0.5)
    crossings = count(i -> deviations[i-1] * deviations[i] < 0, 2:length(deviations))
    fit = nothing
    if minimum(deviations) < 0 < maximum(deviations)
        println(
            io,
            "Resolved equilibrium crossings: ",
            crossings,
            " (population floor ",
            deviation_floor,
            ")",
        )
        println(
            io,
            "No single constant rate describes this full recorded trajectory: it crosses equilibrium.",
        )
        println(
            io,
            "Inspect any longer-time tail separately before assigning a kinetic rate.",
        )
        near_equilibrium && println(
            io,
            "A small initial displacement does not by itself establish a kinetic regime.",
        )
    else
        # For other parameters this is only a candidate fit. Even high R²
        # requires stability across several post-transient fitting windows.
        window = (last(times) / 2, last(times))
        try
            fit = tunnelling_rates(
                times,
                population;
                equilibrium_product = 0.5,
                tspan = window,
            )
            @printf(
                io,
                "Candidate late-window forward/backward rate: %.6g / %.6g; log R² %.6f\n",
                fit.forward_rate,
                fit.backward_rate,
                fit.r_squared
            )
            println(
                io,
                "Check window dependence before interpreting this candidate as a kinetic rate.",
            )
        catch error
            error isa ArgumentError || rethrow()
            println(io, "No usable late-window rate fit: ", sprint(showerror, error))
        end
    end
    @printf(io, "Maximum norm error: %.3e\n", maximum(abs, d.norm .- 1))
    @printf(
        io,
        "Maximum boundary weights (q, p): %.3e, %.3e\n",
        maximum(d.boundary_q),
        maximum(d.boundary_p)
    )
    @printf(
        io,
        "Maximum spectral tails (q, p): %.3e, %.3e\n",
        maximum(d.tail_q),
        maximum(d.tail_p)
    )
    @printf(
        io,
        "Minimum marginal densities (q, p): %.3e, %.3e\n",
        minimum(minimum(position_density(W, grid)) for W in states),
        minimum(minimum(momentum_density(W, grid)) for W in states)
    )
    println(
        io,
        "Vary depth, Padé order, grid spacings, box sizes and ODE tolerances independently.",
    )
    println(
        io,
        "Bath-coupled populations measure total transfer; low temperature alone does not partition its mechanisms.",
    )
    normalized_response = (population .- 0.5) ./ (first(population) - 0.5)
    return (;
        population,
        isolated,
        normalized_response,
        peak,
        transfer_time,
        fit,
        diagnostics = d,
    )
end

function render_double_well(data, result, output_dir)
    (; times, states, grid, potential, energies, kT, hbar) = data
    (; population, isolated, peak) = result
    near_equilibrium = data.preparation === :biased_equilibrium
    mkpath(output_dir)
    qview = (-3.0, 3.0)
    colors = ("#2563eb", "#d97706", "#0d9488")
    qcurve = range(qview...; length = 601)
    well = plot(
        qcurve,
        potential.(qcurve);
        color = "#334155",
        lw = 2.5,
        xlims = qview,
        ylims = (0.0, 1.55),
        xlabel = "Position q",
        ylabel = "Energy",
        title = "Two levels below the barrier",
        label = "V(q)",
    )
    for i in eachindex(energies)
        hline!(
            well,
            [energies[i]];
            color = colors[i],
            ls = :dash,
            lw = 1.5,
            label = "E$(i-1) = $(round(energies[i]; digits=4))",
        )
    end
    hline!(well, [1.0]; color = :gray, ls = :dot, label = "Barrier = 1")
    if near_equilibrium
        plot!(
            well,
            qcurve,
            potential.(qcurve) .+ data.bias .* qcurve;
            color = colors[3],
            ls = :dot,
            lw = 2,
            label = "Preparation: V(q) + $(data.bias) q",
        )
    end

    populations = plot(
        times,
        population;
        color = colors[1],
        lw = 2.5,
        xlabel = "Time",
        ylabel = "Right-well probability",
        ylims = near_equilibrium ? :auto : (0, 1),
        title = near_equilibrium ? "Relaxation after a small right-well depletion" :
                "Damped tunnelling oscillation",
        label = "HEOM, cold bath",
        legend = :topright,
    )
    plot!(
        populations,
        times,
        isolated;
        color = "#94a3b8",
        ls = :dash,
        lw = 2,
        label = near_equilibrium ? "Doublet frequency reference" : "Isolated doublet",
    )
    hline!(populations, [0.5]; color = :gray, ls = :dot, label = "Equilibrium")
    scatter!(
        populations,
        [times[peak]],
        [population[peak]];
        color = colors[2],
        ms = 4,
        label = "Peak $(round(population[peak]; digits=3))",
    )

    density = plot(;
        xlabel = "Position q",
        ylabel = "Position density",
        xlims = qview,
        title = near_equilibrium ? "Nearly thermal position density" :
                "Transfer from the left well",
        legend = :topright,
    )
    for (j, i) in enumerate((1, peak, length(times)))
        plot!(
            density,
            grid.q,
            position_density(states[i], grid);
            color = colors[j],
            lw = 2,
            label = "t = $(round(times[i]; digits=1))",
        )
    end
    if near_equilibrium
        phase = plot(
            times,
            result.normalized_response;
            color = colors[1],
            lw = 2.5,
            xlabel = "Time",
            ylabel = "[P(right,t) - 0.5] / [P(right,0) - 0.5]",
            title = "A positive rate would give monotone decay",
            label = "Normalized population response",
            legend = :topright,
        )
        hline!(phase, [0.0]; color = :gray, ls = :dot, label = "Equilibrium")
    else
        Wpeak = states[peak]
        limit = maximum(abs, Wpeak)
        phase = wignerplot(
            Wpeak,
            grid;
            xlims = qview,
            ylims = (-3.0, 3.0),
            clims = (-limit, limit),
            title = "Wigner function at peak transfer",
            colorbar = true,
        )
    end
    figure = plot(
        well,
        populations,
        density,
        phase;
        layout = (2, 2),
        size = (1300, 850),
        dpi = 150,
        margin = 5Plots.mm,
        background_color = :white,
        plot_title = near_equilibrium ?
                     "Near-equilibrium double well | kT=$(kT), tilt=$(data.bias)" :
                     "Low-temperature double well | kT=$(kT), hbar=$(hbar)",
        guidefontsize = 10,
        tickfontsize = 9,
        titlefontsize = 12,
        legendfontsize = 8,
    )
    savefig(figure, joinpath(output_dir, "cold_double_well.png"))
    open(joinpath(output_dir, "populations.tsv"), "w") do io
        println(
            io,
            "time\tproduct_population\tdoublet_frequency_reference\tnormalized_response\tnorm\tenergy",
        )
        for i in eachindex(times)
            println(
                io,
                join(
                    (
                        times[i],
                        population[i],
                        isolated[i],
                        result.normalized_response[i],
                        result.diagnostics.norm[i],
                        result.diagnostics.energy[i],
                    ),
                    '\t',
                ),
            )
        end
    end
    open(joinpath(output_dir, "results.txt"), "w") do io
        println(io, "Grid: ", grid)
        println(
            io,
            "Bath: lambda=$(data.reorganization), cutoff=$(data.cutoff), Padé=$(data.pade); depth=$(data.depth)",
        )
        report_transfer_rates(data; io)
    end
    println("Figure, population data and results saved to ", output_dir)
    return figure
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_dir =
        isempty(ARGS) ? mktempdir(; prefix = "heom-cold-double-well-") : abspath(only(ARGS))
    data = double_well_trajectory()
    result = report_transfer_rates(data)
    render_double_well(data, result, output_dir)
end
