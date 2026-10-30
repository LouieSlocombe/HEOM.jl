# Select files with Pkg.test(test_args = ["grid", "heom_solvers"]).
# With no file names, run everything, including the package quality checks.
const TEST_FILES = [
    "grid",
    "derivatives",
    "wigner_moyal_rhs",
    "harmonic_oscillator",
    "initial_states",
    "caldeira_leggett",
    "heom",
    "pade_bath",
    "brownian_bath",
    "composite_bath",
    "bath_diagnostics",
    "bath_degeneracy",
    "heom_scaling",
    "generalized_heom",
    "heom_solvers",
    "driven",
    "spectroscopy",
    "spectroscopy_benchmarks",
    "equilibrium",
    "analytic_benchmarks",
    "cold_strong_benchmarks",
    "hierarchy_stability",
    "equilibrium_harmonic",
    "bath_compression",
    "aaa_bath",
    "observables",
    "diagnostics",
    "populations",
    "tunnelling",
    "plotting",
    "animation",
    "float32",
    "device_backend",
    "quality",
]
const TEST_TIMINGS = "--timings" in ARGS
requested_tests = unique(replace(arg, r"\.jl$" => "") for arg in ARGS if arg != "--timings")
unknown_tests = setdiff(requested_tests, TEST_FILES)
isempty(unknown_tests) || throw(
    ArgumentError(
        "Unknown test selection: $(join(unknown_tests, ", ")). Choose from: $(join(TEST_FILES, ", ")).",
    ),
)
selected_tests = isempty(requested_tests) ? TEST_FILES : requested_tests

using ADTypes
using HEOM
using JLArrays
using KernelAbstractions
using LinearAlgebra
using LinearSolve
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqVerner
using SciMLBase: ODEProblem
using SparseArrays
using Test

include("helpers.jl")
include("references.jl")

@testset "HEOM" begin
    for name in selected_tests
        elapsed = @elapsed @testset "$name" begin
            include("$name.jl")
        end
        if TEST_TIMINGS
            println("Test file ", rpad(name, 26), " ", round(elapsed; digits = 3), " s")
            flush(stdout)
        end
    end
end
