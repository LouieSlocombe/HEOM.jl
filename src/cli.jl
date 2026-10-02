const CLI_HELP = """
Usage: heom [name]

Print a friendly greeting.

Arguments:
  name        Name to greet (default: World)

Options:
  -h, --help  Show this help message
  --          Treat the following argument as a name
"""

"""
    main(args::AbstractVector{<:AbstractString}=ARGS; io::IO=stdout, err::IO=stderr)

Run the greeting command. Return `0` on success or help and `2` for invalid
arguments. Output goes to `io`; usage errors go to `err`. This function does
not exit the Julia process.
"""
function main(
    args::AbstractVector{<:AbstractString} = ARGS;
    io::IO = stdout,
    err::IO = stderr,
)
    if length(args) == 1 && first(args) in ("-h", "--help")
        print(io, CLI_HELP)
        return 0
    end

    positional = args
    if !isempty(args) && first(args) == "--"
        positional = args[2:end]
    elseif any(arg -> startswith(arg, "-"), args)
        println(err, "Error: unknown option. Use --help for usage or -- before a name.")
        return 2
    end

    if length(positional) > 1
        println(err, "Error: expected at most one name. Use --help for usage.")
        return 2
    end

    name = isempty(positional) ? "World" : only(positional)
    print_hello(name; io = io)
    return 0
end
