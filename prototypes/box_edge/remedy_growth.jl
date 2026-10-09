# Growth rate of the real (matrix-free) generator with each remedy, Spectral ±8/64.
#     julia --project=. remedy_growth.jl <bath> <depth>
include("remedies.jl")
include(joinpath(@__DIR__, "..", "nz_terminator", "closure.jl"))
# The scalar Markov closure needs no resolvent, so allow the spectral operator.
member_generator(op::HEOM.WignerHEOM{<:HEOM.SpectralWignerMoyal}) = spzeros(1, 1)

struct Combined{C,R}
    closure::C
    remedy::R
end
function combined_heom!(dU, U, c::Combined, t)
    terminated_heom!(dU, U, c.closure, t)
    σ = c.remedy.σ
    for a in 2:size(U, 3)
        @views dU[:, :, a] .-= σ .* U[:, :, a]
    end
end

bname, depth = Symbol(ARGS[1]), parse(Int, ARGS[2])
op = heom_operator(box(8.0, 64); mass = 1.0, potential = HARMONIC, bath = BATHS[bname], depth, scaled = true)
cases = (
    "hard cutoff" => (op, heom!),
    "absorb q0=4 σ0=20" => (r = remedied(op; absorb = 20.0, q0 = 4.0); (r, remedied_heom!)),
    "absorb q0=3 σ0=20" => (r = remedied(op; absorb = 20.0, q0 = 3.0); (r, remedied_heom!)),
    "taper q0=4" => (r = remedied(op; taper = true, q0 = 4.0); (r, remedied_heom!)),
    "markov_full" => (terminate(op, :markov_full), terminated_heom!),
    "markov_full + absorb q0=4" => (Combined(terminate(op, :markov_full), remedied(op; absorb = 20.0, q0 = 4.0)), combined_heom!),
)
for (label, (p, f)) in cases
    t = @elapsed g, _ = growth_rate(op; T = 30.0, rhs! = (dU, U, _, t) -> f(dU, U, p, t))
    @printf("%-12s depth %d  %-28s growth %8.4f  [%4.0fs]\n", bname, depth, label, g, t)
    flush(stdout)
end
