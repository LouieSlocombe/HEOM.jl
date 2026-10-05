# A displaced quantum oscillator relaxing in a Drude–Lorentz bath.
# Requires HEOM, OrdinaryDiffEqVerner and Plots in the active environment.
# Run: julia --project=/path/to/environment examples/animated_heom_sho.jl [output_dir]
# For headless machines, set GKSwstype=100 before starting Julia.
using HEOM, OrdinaryDiffEqVerner, Plots, LinearAlgebra, Printf

include("heom_gaussian_reference.jl")

function simulate_heom_oscillator(; depth = 10)
    mass, omega, hbar = 1.0, 1.0, 1.0
    reorganization, cutoff, kT = 0.2, 1.2, 0.8
    grid = PhaseSpaceGrid((-8.0, 8.0), 64, (-8.0, 8.0), 64)
    potential = harmonic_potential(; mass, omega)
    mu0 = [2.0, 0.0]
    covariance0 = [hbar/(2mass*omega) 0.0; 0.0 hbar*mass*omega/2]
    W0 = on_grid(
        (q, p) -> coherent_wigner(q, p; q0 = mu0[1], p0 = mu0[2], mass, omega, hbar),
        grid,
    )
    bath = drude_lorentz_bath(; reorganization, cutoff, kT, hbar, matsubara = 1)
    # The bath supplies the counterterm internally. A matrix initial state sets
    # higher ADOs to zero, retaining the factorized-bath initial slip.
    # One Moyal term is exact for a quadratic potential; scaling is numerical only.
    op = heom_operator(grid; mass, potential, bath, depth, scaled = true, moyal_terms = 1)
    println("Solving HEOM with ", length(hierarchy_indices(op)), " ADOs...")
    sol = solve(
        heom_problem(W0, (0.0, 20.0), op),
        Vern7();
        saveat = 0.1,
        abstol = 1e-9,
        reltol = 1e-9,
    )
    @assert sol.t[end] == 20.0 "The solver did not reach the final time"
    states = [copy(physical_wigner(U)) for U in sol.u]
    d = diagnostics(sol; potential)

    # Check the full Wigner distribution, rather than only its low moments,
    # against the exact Gaussian solution of this finite exponential bath.
    reference = oscillator_gaussian_reference(sol.t, mu0, covariance0, bath; mass, omega)
    mean_error, state_error = 0.0, 0.0
    for (i, exact) in enumerate(reference)
        precision = inv(exact.covariance)
        amplitude = 1 / (2π * sqrt(det(exact.covariance)))
        Wexact = on_grid(grid) do q, p
            dq, dp = q - exact.mean[1], p - exact.mean[2]
            amplitude * exp(
                -(
                    precision[1, 1] * dq^2 +
                    2precision[1, 2] * dq * dp +
                    precision[2, 2] * dp^2
                ) / 2,
            )
        end
        mean_error = max(mean_error, maximum(abs, [d.mean_q[i], d.mean_p[i]] - exact.mean))
        state_error = max(state_error, maximum(abs, states[i] - Wexact))
    end
    norm_drift = maximum(abs, d.norm .- first(d.norm))
    boundary = max(maximum(d.boundary_q), maximum(d.boundary_p))
    println("Maximum centroid error: ", mean_error)
    println("Maximum full Wigner error: ", state_error)
    println("Maximum norm drift: ", norm_drift)
    println("Maximum boundary weight: ", boundary)
    @assert mean_error < 1e-6
    @assert state_error < 1e-5
    @assert norm_drift < 1e-8
    @assert boundary < 1e-7
    # This tests hierarchy/grid accuracy for the retained bath expansion.
    # Matsubara convergence is a separate approximation, as in examples/heom_sho.jl.
    return (;
        grid,
        states,
        d,
        mass,
        omega,
        hbar,
        reorganization,
        cutoff,
        kT,
        depth,
        nados = length(hierarchy_indices(op)),
        checks = (; mean_error, state_error, norm_drift, boundary),
    )
end

function oscillator_frame(data, i)
    (; grid, states, d, mass, omega) = data
    blue, orange, purple = "#176b9c", "#d46a25", "#7353a6"
    W = states[i]
    qmean, pmean = d.mean_q[i], d.mean_p[i]
    phase = wignerplot(
        W,
        grid;
        xlims = (-4.5, 4.5),
        ylims = (-4.5, 4.5),
        clims = (-1 / π, 1 / π),
        color = :RdBu,
        title = "Wigner distribution",
        xlabel = "Position q",
        ylabel = "Momentum p",
        colorbar_title = "W(q, p)",
    )
    plot!(phase, d.mean_q[1:i], d.mean_p[1:i]; color = "#333333", lw = 1.5, label = "")
    scatter!(phase, [qmean], [pmean]; color = orange, ms = 5, label = "")

    density = plot(
        grid.q,
        position_density(first(states), grid);
        color = "#a3acb7",
        linestyle = :dash,
        lw = 2,
        label = "Initial",
        title = "Position probability",
        xlabel = "Position q",
        ylabel = "P(q)",
        xlims = (-4.5, 4.5),
        ylims = (0, 0.65),
        legend = :topleft,
    )
    plot!(
        density,
        grid.q,
        position_density(W, grid);
        color = blue,
        lw = 3,
        fillrange = 0,
        fillalpha = 0.12,
        label = "Current",
    )
    vline!(density, [qmean]; color = orange, lw = 1.5, label = "Mean q")

    motion = plot(
        d.t,
        d.mean_q;
        color = blue,
        alpha = 0.18,
        lw = 2,
        label = "",
        title = "Damped oscillation",
        xlabel = "Time (1 / omega)",
        ylabel = "Mean q, p",
        xlims = (0, last(d.t)),
        ylims = (-2.5, 2.5),
        legend = :topright,
    )
    plot!(motion, d.t, d.mean_p; color = orange, alpha = 0.18, lw = 2, label = "")
    plot!(motion, d.t[1:i], d.mean_q[1:i]; color = blue, lw = 2.5, label = "Mean q")
    plot!(motion, d.t[1:i], d.mean_p[1:i]; color = orange, lw = 2.5, label = "Mean p")
    vline!(motion, [d.t[i]]; color = "#777777", lw = 1, label = "")

    coherent_energy = (d.mean_p .^ 2) ./ (2mass) + (mass * omega^2 / 2) .* d.mean_q .^ 2
    energies = plot(
        d.t,
        d.energy;
        color = purple,
        alpha = 0.18,
        lw = 2,
        label = "",
        title = "Oscillator energy",
        xlabel = "Time (1 / omega)",
        ylabel = "Energy (hbar omega)",
        xlims = (0, last(d.t)),
        ylims = (0, 3.1),
        legend = :topright,
    )
    plot!(energies, d.t[1:i], d.energy[1:i]; color = purple, lw = 2.5, label = "Oscillator")
    plot!(
        energies,
        d.t[1:i],
        coherent_energy[1:i];
        color = blue,
        lw = 2,
        label = "Centroid motion",
    )
    vline!(energies, [d.t[i]]; color = "#777777", lw = 1, label = "")

    return plot(
        phase,
        density,
        motion,
        energies;
        layout = (2, 2),
        size = (1040, 800),
        plot_title = @sprintf(
            "HEOM | Quantum harmonic oscillator in a thermal bath | t = %.1f",
            d.t[i]
        ),
        plot_titlefontsize = 15,
        titlefontsize = 12,
        guidefontsize = 10,
        tickfontsize = 9,
        legendfontsize = 9,
        background_color = "#fcfcfd",
        foreground_color = "#243446",
        margin = 5Plots.mm,
    )
end

function render_heom_oscillator(data, output_dir)
    mkpath(output_dir)
    anim = Animation()
    for i in eachindex(data.states)
        fig = oscillator_frame(data, i)
        frame(anim, fig)
        if i in (1, 51, length(data.states))
            savefig(fig, joinpath(output_dir, "heom_oscillator_frame_$(i).png"))
        end
        i % 50 == 0 && println("Rendered ", i, " / ", length(data.states), " frames")
    end
    # The final state is held for one second, then the animation restarts.
    for _ in 1:20
        frame(anim, oscillator_frame(data, length(data.states)))
    end
    gif(anim, joinpath(output_dir, "heom_oscillator.gif"); fps = 20)
    mp4(anim, joinpath(output_dir, "heom_oscillator.mp4"); fps = 20)
    open(joinpath(output_dir, "observables.csv"), "w") do io
        println(io, "time,mean_q,mean_p,energy,purity,norm")
        for i in eachindex(data.d.t)
            println(
                io,
                join(
                    (
                        getproperty(data.d, f)[i] for
                        f in (:t, :mean_q, :mean_p, :energy, :purity, :norm)
                    ),
                    ",",
                ),
            )
        end
    end
    open(joinpath(output_dir, "validation.txt"), "w") do io
        println(io, "m = omega = hbar = 1; q0 = 2; p0 = 0")
        println(
            io,
            "Drude bath: lambda = ",
            data.reorganization,
            "; cutoff = ",
            data.cutoff,
            "; kT = ",
            data.kT,
        )
        println(
            io,
            "Matsubara poles = 1; scaled hierarchy depth = ",
            data.depth,
            "; ADOs = ",
            data.nados,
            "; grid = 64 x 64 on [-8, 8)^2",
        )
        for (name, value) in pairs(data.checks)
            println(io, name, " = ", value)
        end
        println(
            io,
            "Initial / final energy = ",
            first(data.d.energy),
            " / ",
            last(data.d.energy),
        )
        println(io, "The Gaussian reference uses the same finite bath expansion.")
    end
    println("Animations and observables saved to ", output_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_dir =
        isempty(ARGS) ? mktempdir(; prefix = "animated-heom-") : abspath(only(ARGS))
    render_heom_oscillator(simulate_heom_oscillator(), output_dir)
end
