# Frozen-symbol bound for the absorber and taper, production grid ±8/64 Spectral.
include("remedies.jl")
for bname in (:cold_strong, :warm, :cold_weak), depth in (2, 4, 6)
    op = heom_operator(box(8.0, 64); mass = 1.0, potential = HARMONIC, bath = BATHS[bname],
        depth, scaled = true)
    println("\n== $bname depth $depth: hard cutoff bound $(round(local_abscissa(op)[1], digits = 3)) ==")
    for q0 in (2.0, 4.0, 6.0)
        cells = String[]
        for σ0 in (2.0, 5.0, 10.0, 20.0, 50.0)
            p = remedy_profiles(remedied(op; absorb = σ0, q0))
            r = local_abscissa(op; absorb = p.absorb)
            push!(cells, @sprintf("σ₀=%-4g %6.3f@q=%5.2f", σ0, r[1], r[2]))
        end
        p = remedy_profiles(remedied(op; taper = true, q0))
        r = local_abscissa(op; taper = p.taper)
        @printf("q₀=%.0f absorb: %s | taper %6.3f@q=%5.2f\n", q0, join(cells, "  "), r[1], r[2])
    end
    flush(stdout)
end
