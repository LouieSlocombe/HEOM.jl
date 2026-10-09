# Frozen-point abscissa over the ±8 box (κ ≤ 12) for each hierarchy form.
include("frozen_forms.jl")
using HEOM
const FORMS = (:package, :markov_full, :standard, :gatto)
function form_matrix(form, c, ν, D, x, xp, d)
    form === :markov_full ? markov_frozen_matrix(c, ν, D, x, xp, d)[1] : frozen_matrix(form, c, ν, D, x, xp, d)[1]
end
function scan(form, b, depth; qs = -8:0.5:8, κs = 0:0.4:12)
    c, ν, D = b.coefficients, b.rates, b.diffusion
    best, onset = (-Inf, NaN, NaN), Inf
    for q in qs, κ in κs
        x, xp = q - κ / 2, q + κ / 2
        r = maximum(real, eigvals(form_matrix(form, c, ν, D, x, xp, depth)))
        r > best[1] && (best = (r, q, κ))
        r > 1e-8 && (onset = min(onset, abs(q)))
    end
    return best, onset
end
for (name, b) in (("warm Padé 1", drude_lorentz_pade_bath(; reorganization = 0.2, cutoff = 1.0, kT = 1.0, pade = 1)),
                  ("cold weak Padé 2", drude_lorentz_pade_bath(; reorganization = 0.05, cutoff = 0.5, kT = 0.1, pade = 2)),
                  ("cold strong Padé 2", drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 2)))
    println("\n== $name: max frozen Re λ over |q| ≤ 8, κ ≤ 12 (smallest unstable |q|) ==")
    for depth in (2, 4, 6)
        cells = map(FORMS) do form
            form === :gatto && depth == 6 && occursin("Padé 2", name) && return @sprintf("%-12s %s", form, "skipped (C(12,6)=924)")
            (r, q, κ), onset = scan(form, b, depth)
            @sprintf("%-12s %7.3f at q=%5.1f κ=%4.1f (onset |q|=%4.1f)", form, r, q, κ, onset)
        end
        println("depth $depth\n  ", join(cells, "\n  "))
        flush(stdout)
    end
end
