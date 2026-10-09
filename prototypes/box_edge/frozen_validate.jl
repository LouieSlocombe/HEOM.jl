# Check the three frozen forms against the exact influence functional.
include("frozen_forms.jl")
using HEOM
for (name, b) in (("warm Padé 1", drude_lorentz_pade_bath(; reorganization = 0.2, cutoff = 1.0, kT = 1.0, pade = 1)),
                  ("cold strong Padé 1", drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 1)))
    c, ν, D = b.coefficients, b.rates, b.diffusion
    println("== $name ==")
    for (x, xp) in ((0.3, -0.2), (1.0, 0.5), (2.0, -1.0))
        ref = exact_frozen(c, ν, D, x, xp, 3.0)
        cells = [@sprintf("%s d=%d %.1e", f, d, abs(frozen_root(f, c, ν, D, x, xp, d, 3.0) - ref) / abs(ref))
                 for f in (:package, :standard, :gatto), d in (2, 6, 10)]
        @printf("x=%4.1f x'=%4.1f |ρ|=%.3e  rel err: %s\n", x, xp, abs(ref), join(cells, "  "))
    end
end
