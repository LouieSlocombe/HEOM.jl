# Element-type and array-type helpers. Operators are always built on the host in Float64 and
# then converted once by `move_operator`, so the default path stays untouched and symbols
# are rounded a single time into the requested precision.

check_float_type(::Type{T}) where {T} =
    T === Float32 || T === Float64 ? T :
    throw(ArgumentError("float_type must be Float32 or Float64, got $T"))
check_float_type(T) = throw(ArgumentError("float_type must be Float32 or Float64, got $T"))

check_device(::Type{<:Array}) = Array
check_device(device::Type{<:AbstractArray}) = device
check_device(device) = throw(
    ArgumentError("device must be an array type such as Array or CuArray, got $device"),
)

# Host copy with the requested floating-point type (shares `A` when it already matches).
float_like(::Type{T}, A::AbstractArray{<:Real}) where {T<:AbstractFloat} =
    convert(Array{T,ndims(A)}, A)
float_like(::Type{T}, A::AbstractArray{<:Complex}) where {T<:AbstractFloat} =
    convert(Array{Complex{T},ndims(A)}, A)

"""
    to_device(device, T, A)

Copy the host array `A` to the array type `device` (`Array`, `CuArray`, …) with element type
`T` (`Complex{T}` for complex input).
"""
to_device(device, ::Type{T}, A::AbstractArray) where {T<:AbstractFloat} =
    device === Array ? float_like(T, A) : device(float_like(T, A))
to_device(device, ::Type{Int32}, v::AbstractVector{<:Integer}) =
    device === Array ? convert(Vector{Int32}, v) : device(convert(Vector{Int32}, v))

"""
    to_host(x)

Host `Float64` array holding the real array `x`. A host `Array{Float64}`, or a view of one,
is returned as it is; any other real array (a device array, another element type) is
copied. Arrays that are not real are copied without a change of type, so that callers still
reject them.
"""
to_host(x::Array{Float64}) = x
to_host(x::SubArray{Float64,<:Any,<:Array{Float64}}) = x
to_host(x::AbstractArray{<:Real}) = convert(Array{Float64,ndims(x)}, Array(x))
to_host(x::AbstractArray) = Array(x)

"""
    state_like(prototype, T, host)

New array of the type of `prototype` (a device array or an `Array`) and element type `T`
holding the host data `host`.
"""
function state_like(prototype, ::Type{T}, host::AbstractArray) where {T<:AbstractFloat}
    out = similar(prototype, T, size(host)...)
    copyto!(out, float_like(T, host))
    return out
end

# `host` as an array of the type and precision of `like` (shared when it already matches).
array_like(::Array{ComplexF64}, host::Array{ComplexF64}) = host
function array_like(like, host)
    return copyto!(similar(like, size(host)...), float_like(real(eltype(like)), host))
end

# The array type (`Array`, `CuArray`, …) without its parameters.
array_wrapper(x) = Base.typename(typeof(x)).wrapper

# The `T` time span avoids promoting device kernels through Float64 step sizes.
convert_tspan(::Type{Float64}, tspan) = tspan
convert_tspan(::Type{T}, tspan) where {T<:AbstractFloat} = map(T, Tuple(tspan))

"""
    physical_save_idxs(op)

Linear indices of the physical (first) hierarchy member of a HEOM state, for the
`save_idxs` keyword of `solve`. Saving only these keeps device memory small; the saved
states are vectors that [`physical_wigner`](@ref) and the analysis functions reshape.
"""
physical_save_idxs(op) = 1:prod(size(op.grid))

# Convert a state computed on the host back to the array type and precision of `like`.
function convert_like(like::Array{Float64}, host)
    return host
end
function convert_like(like, host)
    T = eltype(like)
    converted = float_like(T, host)
    # Overflow or underflow in the narrower type would silently change the hierarchy.
    all(isfinite, converted) &&
    all(i -> iszero(host[i]) || !iszero(converted[i]), eachindex(host)) ||
        throw(ArgumentError("converted values are outside the range of $T"))
    return state_like(like, T, converted)
end
