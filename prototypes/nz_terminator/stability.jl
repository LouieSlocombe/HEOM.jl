# Spectral abscissa (largest real eigenvalue) of the closed hierarchy generator, by
# dense diagonalisation on a small grid. Run with
#     julia --project=prototypes/nz_terminator prototypes/nz_terminator/stability.jl
#
# A physical closed hierarchy has one zero eigenvalue (the steady state) and all
# others in the left half-plane. A positive abscissa is a long-time instability.
# This does not measure transient (non-normal) amplification.
using Printf
include("closure.jl")

function abscissa(cl)
    λ = eigvals(dense_generator(cl))
    sort!(λ; by = real, rev = true)
    return λ[1:3]
end

function run_case(name; bath, potential, moyal_terms, depths)
    grid = PhaseSpaceGrid((-6.0, 6.0), 20, (-6.0, 6.0), 20)
    common = (; mass = 1.0, potential, discretization = FiniteDifference(4), moyal_terms)
    println("\n== $name ==")
    println("rates    ", round.(bath.rates; sigdigits = 4))
    println("residues ", round.(bath.coefficients; sigdigits = 4), "  D = ", bath.diffusion)
    @printf("%-12s %5s %7s %12s %12s %12s\n", "closure", "depth", "size", "Re λ₁", "Re λ₂", "Re λ₃")
    for depth in depths
        op = heom_operator(grid; common..., bath, depth, scaled = true)
        for kind in CLOSURES
            λ = abscissa(terminate(op, kind))
            @printf("%-12s %5d %7d %12.3e %12.3e %12.3e\n",
                kind, depth, prod(size(grid)) * length(op.indices), real.(λ)...)
        end
        flush(stdout)
    end
end

harmonic = harmonic_potential(; mass = 1.0, omega = 1.0)
run_case("harmonic, cold strong (λ=0.8, γ=0.5, kT=0.1), Padé 1";
    bath = drude_lorentz_pade_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, pade = 1),
    potential = harmonic, moyal_terms = 1, depths = 0:5)
run_case("harmonic, very strong (λ=2, γ=0.5, kT=0.1), Padé 1";
    bath = drude_lorentz_pade_bath(; reorganization = 2.0, cutoff = 0.5, kT = 0.1, pade = 1),
    potential = harmonic, moyal_terms = 1, depths = 0:5)
run_case("harmonic, cold strong, Matsubara 1 (negative Drude residue)";
    bath = drude_lorentz_bath(; reorganization = 0.8, cutoff = 0.5, kT = 0.1, matsubara = 1),
    potential = harmonic, moyal_terms = 1, depths = 0:5)
run_case("quartic q²/2+0.08q⁴, very strong, Padé 1";
    bath = drude_lorentz_pade_bath(; reorganization = 2.0, cutoff = 0.5, kT = 0.1, pade = 1),
    potential = q -> q^2 / 2 + 0.08q^4, moyal_terms = 2, depths = 0:5)
