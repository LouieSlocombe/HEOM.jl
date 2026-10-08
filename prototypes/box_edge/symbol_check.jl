# Does the frozen (q, κ) symbol predict the dense spectral abscissa?
include("common.jl")
for (bname, bath) in pairs((cold_strong1 = BATHS.cold_strong1, warm = BATHS.warm))
    for disc in (FiniteDifference(4), Spectral()), L in (6.0, 10.0), depth in (1, 2, 3)
        grid = box(L, 20)
        moyal = disc isa Spectral ? nothing : 1
        op = heom_operator(grid; mass = 1.0, potential = HARMONIC, bath, depth,
            discretization = disc, moyal_terms = moyal, scaled = true)
        loc = local_abscissa(op)
        modes = leading_modes(dense_generator(op), op; count = 1)
        m = modes[1]
        @printf("%-13s %-22s L=%4.1f d=%d  local %7.3f (q=%5.2f, j=%2d)  dense %7.3f%+7.3fi  outer q %.2f  outer p %.2f\n",
            bname, string(disc), L, depth, loc[1], loc[2], loc[3], real(m.λ), imag(m.λ), m.outer_q, m.outer_p)
        flush(stdout)
    end
end
