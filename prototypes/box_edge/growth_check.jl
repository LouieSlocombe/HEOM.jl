# Validate the power-iteration growth rate against dense eigenvalues.
include("remedies.jl")
bath = BATHS.cold_strong1
for disc in (FiniteDifference(4), Spectral()), depth in (1, 3)
    grid = box(6.0, 20)
    op = heom_operator(grid; mass = 1.0, potential = HARMONIC, bath, depth,
        discretization = disc, moyal_terms = disc isa Spectral ? nothing : 1, scaled = true)
    λ = leading_modes(dense_generator(op), op; count = 1)[1].λ
    g, _ = growth_rate(op; T = 40.0)
    @printf("%-20s depth %d  dense %.4f   power iteration %.4f\n", string(disc), depth, real(λ), g)
end
