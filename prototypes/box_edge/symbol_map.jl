# Where in (q, κ) is the frozen hierarchy unstable? Production grids.
include("common.jl")
function growth_map(op)
    s1, s2 = momentum_symbols(op)
    q = op.grid.q
    G = [maximum(real, eigvals!(local_matrix(op, q[i], s1[j], s2[j]))) for i in eachindex(q), j in eachindex(s1)]
    return G
end
for bname in (:cold_strong, :warm, :cold_weak), depth in (2, 4, 6)
    bath = BATHS[bname]
    grid = box(8.0, 64)
    op = heom_operator(grid; mass = 1.0, potential = HARMONIC, bath, depth, scaled = true)
    s1, _ = momentum_symbols(op)
    κ = imag.(s1)
    G = growth_map(op)
    println("\n== $bname depth $depth, Spectral ±8/64 (κmax = $(round(maximum(κ), digits=2))), max local rate $(round(maximum(G), digits=3)) ==")
    # Smallest |q| at which any κ is unstable, and the unstable κ range there.
    unstable_q = [maximum(G[i, :]) > 1e-10 for i in axes(G, 1)]
    qc = minimum(abs.(grid.q[unstable_q]); init = Inf)
    println("smallest unstable |q| = $qc")
    for qv in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, -8.0)
        i = argmin(abs.(grid.q .- qv))
        row = G[i, :]
        js = findall(>(1e-10), row)
        rng = isempty(js) ? "none" : @sprintf("κ ∈ [%.2f, %.2f]", κ[first(js)], κ[last(js)])
        @printf("  q=%5.2f  max rate %7.3f at κ=%5.2f   unstable %s\n", grid.q[i], maximum(row), κ[argmax(row)], rng)
    end
    flush(stdout)
end
