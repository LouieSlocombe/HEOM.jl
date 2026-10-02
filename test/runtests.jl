using Aqua
using HEOM
using Test

"""
Run the command-line entry point with separate output and error streams.
"""
function capture_main(args::Vector{String})
    output = IOBuffer()
    errors = IOBuffer()
    status = HEOM.main(args; io = output, err = errors)
    return (status = status, output = String(take!(output)), errors = String(take!(errors)))
end

"""
Run the installed Julia executable without a shell or user startup file.
"""
function capture_cli(args::String...)
    project = pkgdir(HEOM)
    script = joinpath(project, "bin", "heom.jl")
    command = `$(Base.julia_cmd()) --startup-file=no --project=$project $script $args`
    output = IOBuffer()
    errors = IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout = output, stderr = errors))
    return (
        status = process.exitcode,
        output = String(take!(output)),
        errors = String(take!(errors)),
    )
end

@testset "HEOM" begin
    @testset "Greeting" begin
        output = IOBuffer()
        @test print_hello(; io = output) === nothing
        @test String(take!(output)) == "Hello, World!\n"

        for name in ("Julia", "Ada Lovelace", "世界", "")
            @test print_hello(name; io = output) === nothing
            @test String(take!(output)) == "Hello, $(name)!\n"
        end
    end

    @testset "Evenly spaced samples" begin
        samples = @inferred line()
        @test samples isa Vector{Float64}
        @test length(samples) == 100
        @test first(samples) == 0.0
        @test last(samples) == 1.0
        @test all(isapprox.(diff(samples), 1 / 99))

        @test (@inferred line(-1, 1; num = 5)) == [-1.0, -0.5, 0.0, 0.5, 1.0]
        @test line(1, -1; num = 5) == [1.0, 0.5, 0.0, -0.5, -1.0]
        @test line(2, 2; num = 4) == fill(2.0, 4)
        @test line(1 // 4, 3 // 4; num = 3) == [0.25, 0.5, 0.75]
        @test line(2; num = 3) == [2.0, 1.5, 1.0]
        @test line(-2, 8; num = 2) == [-2.0, 8.0]

        empty_samples = @inferred line(2, 8; num = 0)
        @test empty_samples isa Vector{Float64}
        @test isempty(empty_samples)
        @test (@inferred line(2, 8; num = 1)) == [2.0]
        @test_throws ArgumentError line(; num = -1)
        @test_throws ArgumentError line(1, -1; num = -10)
    end

    @testset "Command-line argument handling" begin
        for (args, name) in (
            (String[], "World"),
            (["Julia"], "Julia"),
            (["Ada Lovelace"], "Ada Lovelace"),
            (["世界"], "世界"),
            ([""], ""),
            (["--"], "World"),
            (["--", "Julia"], "Julia"),
            (["--", "--help"], "--help"),
            (["--", "-h"], "-h"),
        )
            result = capture_main(args)
            @test result.status === 0
            @test result.output == "Hello, $(name)!\n"
            @test isempty(result.errors)
        end

        for option in ("-h", "--help")
            result = capture_main([option])
            @test result.status === 0
            @test occursin("usage", lowercase(result.output))
            @test occursin("Print a friendly greeting.", result.output)
            @test isempty(result.errors)
        end

        for args in (
            ["--unknown"],
            ["-x"],
            ["Alice", "Bob"],
            ["--help", "Alice"],
            ["Alice", "--help"],
            ["--", "Alice", "Bob"],
        )
            result = capture_main(args)
            @test result.status === 2
            @test isempty(result.output)
            @test !isempty(result.errors)
        end
    end

    @testset "Executable script" begin
        greeting = capture_cli("Ada Lovelace")
        @test greeting.status == 0
        @test greeting.output == "Hello, Ada Lovelace!\n"
        @test isempty(greeting.errors)

        help = capture_cli("--help")
        @test help.status == 0
        @test occursin("Print a friendly greeting.", help.output)
        @test isempty(help.errors)

        invalid = capture_cli("--unknown")
        @test invalid.status == 2
        @test isempty(invalid.output)
        @test !isempty(invalid.errors)
    end

    @testset "Package quality" begin
        Aqua.test_all(HEOM)
    end
end
