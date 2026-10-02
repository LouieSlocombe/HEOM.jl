"""
    print_hello(name::AbstractString="World"; io::IO=stdout)

Print a friendly greeting to `io` and return `nothing`.

# Examples

```julia
print_hello(\"Ada\") # Hello, Ada!
```
"""
function print_hello(name::AbstractString = "World"; io::IO = stdout)
    println(io, "Hello, ", name, "!")
    return nothing
end

"""
    line(start::Real=0.0, stop::Real=1.0; num::Integer=100)

Return `num` evenly spaced `Float64` values over the closed interval from
`start` to `stop`. Both endpoints are included when `num >= 2`. With `num=0`
the result is empty; with `num=1` it contains only `start`. A negative `num`
throws an `ArgumentError`.

# Examples

```julia
line(-1.0, 1.0; num = 3) # [-1.0, 0.0, 1.0]
```
"""
function line(start::Real = 0.0, stop::Real = 1.0; num::Integer = 100)
    num >= 0 || throw(ArgumentError("num must be nonnegative"))
    num == 0 && return Float64[]
    num == 1 && return [Float64(start)]
    return collect(range(Float64(start), Float64(stop); length = num))
end
