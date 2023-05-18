function plot_wigner_heatmap(x, y, z; name="w_heatmap.pdf", title="")
    xlab = latexstring("Q \\: \\left[\\alpha_0 \\right]")
    ylab = latexstring("P \\: \\left[\\hbar / \\alpha_0 \\right]")
    z = permutedims(z, (2, 1))
    p = heatmap(
        x,
        y,
        z,
        seriescolor=:inferno,
        xlabel=xlab,
        ylabel=ylab,
        title=title,
        xlims=(minimum(x), maximum(x)),
        ylims=(minimum(y), maximum(y)),
        right_margin=10.0Plots.mm,
    )
    savefig(p, joinpath(plot_dump, name))
    display(p)
    return p
end