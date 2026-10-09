include("common.jl")
for disc in (Spectral(), FiniteDifference(4)), (name, bath) in pairs(merge(BATHS, (brownian = brownian_oscillator_bath(; reorganization = 0.2, frequency = 1.0, damping = 0.5, kT = 1.0),))), depth in (0, 2, 3)
    op = heom_operator(box(6.0, 24); mass = 1.0, potential = HARMONIC, bath, depth, discretization = disc,
        moyal_terms = disc isa Spectral ? nothing : 1, scaled = true)
    s = hierarchy_stability(op)
    loc = local_abscissa(op)
    @printf("%-20s %-13s d=%d  package %.6f (|q|=%.2f κ=%.2f radius=%.2f)  prototype %.6f (q=%.2f)\n",
        string(disc), name, depth, s.rate, s.position, s.wavenumber, s.radius, loc[1], loc[2])
end
