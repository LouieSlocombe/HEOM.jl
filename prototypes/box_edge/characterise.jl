# Growth rate of the hard-cutoff generator (power iteration on the matrix-free
# operator) and the frozen-symbol bound, swept along one axis at a time.
#     julia --project=. characterise.jl <sweep>
include("common.jl")

function report(label, op; T = 30.0)
    t = @elapsed g, _ = growth_rate(op; T)
    loc = local_abscissa(op)
    @printf("%-44s members %4d  growth %8.4f  local bound %8.4f (q=%6.2f)  [%5.0fs]\n",
        label, length(op.indices), g, loc[1], loc[2], t)
    flush(stdout)
end

make(grid, bath, depth; disc = Spectral(), scaled = true) = heom_operator(grid; mass = 1.0,
    potential = HARMONIC, bath, depth, discretization = disc,
    moyal_terms = disc isa Spectral ? nothing : 1, scaled)

sweep = Symbol(ARGS[1])
println("== sweep: $sweep (harmonic m=ω=ħ=1; Spectral, scaled, ±8/64 unless stated) ==")
if sweep === :depth
    for bname in (:cold_strong, :warm, :cold_weak), depth in 0:2:8
        report("$bname depth $depth", make(box(8.0, 64), BATHS[bname], depth))
    end
elseif sweep === :box
    # fixed spacing dq = dp = 0.25
    for bname in (:cold_strong, :warm), L in (3.0, 4.0, 6.0, 8.0, 10.0, 12.0)
        report("$bname depth 4 L=$L n=$(Int(8L))", make(box(L, Int(8L)), BATHS[bname], 4))
    end
elseif sweep === :spacing
    for bname in (:cold_strong, :warm), n in (32, 48, 64, 96, 128)
        report("$bname depth 4 L=8 n=$n (dp=$(16/n))", make(box(8.0, n), BATHS[bname], 4))
    end
elseif sweep === :discretisation
    for bname in (:cold_strong, :warm), disc in (Spectral(), FiniteDifference(2), FiniteDifference(4), FiniteDifference(8))
        report("$bname depth 4 $(disc)", make(box(8.0, 64), BATHS[bname], 4; disc))
    end
elseif sweep === :coupling
    for kT in (0.1, 1.0), λ in (0.05, 0.1, 0.2, 0.4, 0.8, 1.6)
        bath = drude_lorentz_pade_bath(; reorganization = λ, cutoff = 0.5, kT, pade = 2)
        report("γ=0.5 kT=$kT λ=$λ depth 4", make(box(8.0, 64), bath, 4))
    end
elseif sweep === :scaling
    for bname in (:cold_strong, :warm), scaled in (true, false), depth in (2, 6)
        report("$bname depth $depth scaled=$scaled", make(box(8.0, 64), BATHS[bname], depth; scaled))
    end
elseif sweep === :baths
    # Other decompositions of the cold, strong Drude bath and an underdamped Brownian bath.
    for (name, bath) in (
        ("Padé 1", BATHS.cold_strong1),
        ("Padé 4", drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 4)),
        ("Matsubara 2", drude_lorentz_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, matsubara = 2)),
        ("Padé 2, no terminator", drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 2, terminator = false)),
        ("Brownian λ=0.2 Ω=1 γ=0.5 kT=1", brownian_oscillator_bath(; reorganization = 0.2, frequency = 1.0, damping = 0.5, kT = 1.0)),
    )
        report("cold strong $name depth 4", make(box(8.0, 64), bath, 4))
    end
end
