# Timing and blow-up of candidate long-time regression tests.
include("common.jl")
for (L, n) in ((6.0, 48), (6.0, 32), (5.0, 32)), (name, bath) in (("cold strong Padé 2", BATHS.cold_strong), ("cold strong Matsubara 1", drude_lorentz_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, matsubara = 1)))
    grid = box(L, n)
    W0 = displaced_gaussian(grid)
    prob = heom_problem(W0, (0.0, 30.0), grid; mass = 1.0, potential = HARMONIC, bath, depth = 2, scaled = true)
    t = @elapsed tb = blowup_time(prob)
    @printf("%-24s ±%.0f/%d depth 2: blow-up at t=%.2f  [%.1fs]\n", name, L, n, tb, t)
end
