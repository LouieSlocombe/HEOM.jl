function animate_QSE_P(q_vec, solu, time; name="Pq.gif", dir=nothing)
    # Get the time loop
    nt = length(time)

    # Get the maximum value for plotting
    y_min = minimum([minimum(i) for i in solu[:][:]])
    y_max = maximum([maximum(i) for i in solu[:][:]])
    @time begin
        # Loop over time
        anim = @animate for i = 1:nt
            plot(
                q_vec,
                solu[i][:],
                ylims=(y_min, y_max),
                lw=3,
                xlabel=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
                ylabel="P(q)",
                legend=false,
                linecolor=:black,
            )
        end
        if dir === nothing
            fig = gif(anim, joinpath(plot_dump, name), fps=60)
        else
            fig = gif(anim, joinpath(dir, name), fps=60)
        end
    end
    os_display(p)
    return p
end

function animate_QSE_P_log(
    q_vec,
    solu,
    time;
    name="Pq_log.gif",
    max_y_min=-20,
    dir=nothing
)
    # Get the time loop
    nt = length(time)
    y_vals = [log10.(abs.(solu[i][:])) for i = 1:nt]
    replace!(y_vals, Inf => NaN)
    # Get the maximum value for plotting
    y_min = minimum([minimum(i) for i in y_vals[:][:]])
    y_max = maximum([maximum(i) for i in y_vals[:][:]])
    # Prevent excessivly small values plotted
    if y_min < max_y_min
        y_min = max_y_min
    end
    @time begin
        # Loop over time
        anim = @animate for i = 1:nt
            plot(
                q_vec,
                y_vals[i][:],
                ylims=(y_min, y_max),
                lw=3,
                xlabel=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
                ylabel="log P(q)",
                legend=false,
                linecolor=:black,
            )
        end
        if dir === nothing
            fig = gif(anim, joinpath(plot_dump, name), fps=60)
        else
            fig = gif(anim, joinpath(dir, name), fps=60)
        end
    end
    os_display(p)
    return p
end