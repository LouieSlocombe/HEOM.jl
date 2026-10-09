# Rows of the characterisation lost when the first sweep was interrupted.
#     julia --project=. characterise2.jl <sweep>
include("common.jl")
function report(label, op; T = 30.0)
    t = @elapsed g, _ = growth_rate(op; T)
    s = hierarchy_stability(op)
    @printf("%-44s members %4d  growth %8.4f  indicator %8.4f (radius %.2f)  [%5.0fs]\n",
        label, length(op.indices), g, s.rate, s.radius, t)
    flush(stdout)
end
make(grid, bath, depth) = heom_operator(grid; mass = 1.0, potential = HARMONIC, bath, depth, scaled = true)
sweep = Symbol(ARGS[1])
println("== sweep: $sweep (harmonic m=ω=ħ=1; Spectral, scaled, ±8/64 unless stated) ==")
if sweep === :box
    for L in (4.0, 6.0, 8.0, 10.0, 12.0)
        report("warm depth 4 L=$L n=$(Int(8L))", make(box(L, Int(8L)), BATHS.warm, 4))
    end
elseif sweep === :coupling
    for λ in (0.05, 0.2, 0.8)
        report("γ=0.5 kT=1 λ=$λ Padé 2 depth 4", make(box(8.0, 64), drude_lorentz_pade_bath(; reorganization = λ, cutoff = 0.5, kT = 1.0, pade = 2), 4))
    end
elseif sweep === :baths
    report("cold strong Padé 2, no terminator, depth 4", make(box(8.0, 64), drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 2, terminator = false), 4))
    report("Brownian λ=0.2 Ω=1 γ=0.5 kT=1, depth 4", make(box(8.0, 64), brownian_oscillator_bath(; reorganization = 0.2, frequency = 1.0, damping = 0.5, kT = 1.0), 4))
    report("cold weak depth 6", make(box(8.0, 64), BATHS.cold_weak, 6))
end
