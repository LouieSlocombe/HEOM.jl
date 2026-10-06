"""
    combine_baths(bath::ExponentialBath, baths::ExponentialBath...)
    combine_baths(baths::AbstractVector{<:ExponentialBath})

Combine independent bath components coupled to the same position into one
[`ExponentialBath`](@ref). The resolved correlation and spectrum are the sums
of the corresponding component diagnostics; the spectrum includes the summed
Markovian remainders. Diffusion and potential counterterms are added. This does
not describe baths coupled to different system operators or correlated bath
components.

All components must have exactly the same `hbar`, and at least one component
must be supplied. Components with no exponential modes are allowed. Coefficients,
rates, and weights retain their component ordering; mixing matrices form a sparse
block diagonal matrix. All arrays are copied, including for a single component.
"""
combine_baths(bath::ExponentialBath, baths::ExponentialBath...) =
    _combine_baths((bath, baths...))

combine_baths(baths::AbstractVector{<:ExponentialBath}) = _combine_baths(baths)

combine_baths() = throw(ArgumentError("at least one bath component is required"))

function _combine_baths(baths)
    isempty(baths) && throw(ArgumentError("at least one bath component is required"))
    ħ = first(baths).hbar
    all(bath -> bath.hbar == ħ, baths) ||
        throw(ArgumentError("bath components must have exactly the same hbar"))
    return ExponentialBath(
        vcat((bath.coefficients for bath in baths)...),
        vcat((bath.rates for bath in baths)...);
        hbar = ħ,
        diffusion = sum(bath.diffusion for bath in baths),
        counterterm = sum(bath.counterterm for bath in baths),
        weights = vcat((bath.weights for bath in baths)...),
        mixing = SparseArrays.blockdiag((bath.mixing for bath in baths)...),
    )
end
