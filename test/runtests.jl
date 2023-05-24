using HEOM
using Test, SafeTestsets

begin
    @time @safetestset "Derivatives" begin
        include("derivatives.jl")
    end
    @time @safetestset "Symbolic eq check" begin
        include("core_eq/symbolic_check.jl")
    end

    @time @safetestset "Simple eq check" begin
        include("core_eq/simple_check.jl")
    end
    @time @safetestset "Free Particle" begin
        include("core_eq/free_particle.jl")
    end
    # @time @safetestset "Harmonic well" begin include("core_eq/harmonic_well.jl") end
    # Morse potential
    # Double well case

end