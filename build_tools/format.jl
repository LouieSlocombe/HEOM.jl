using JuliaFormatter

function main(args)
    if !(isempty(args) || args == ["--fix"] || args == ["--check"])
        println(
            stderr,
            "Usage: julia --project=build_tools build_tools/format.jl [--check|--fix]",
        )
        return 2
    end

    check = args == ["--check"]
    root = dirname(@__DIR__)
    paths = [joinpath(root, dir) for dir in ("src", "test", "bin", "build_tools")]
    already_formatted =
        format(paths; overwrite = !check, verbose = true, throw_on_error = true)
    if check && !already_formatted
        println(stderr, "Formatting differs. Run build_tools/format.jl --fix to apply it.")
        return 1
    end
    return 0
end

exit(main(ARGS))
