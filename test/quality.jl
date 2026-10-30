ENV["GKSwstype"] = "100"
using Aqua
using Plots

@testset "Package quality" begin
    Aqua.test_all(HEOM)
end
