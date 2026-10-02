"""
    Spectral()

Fourier pseudo-spectral discretisation of the phase-space derivatives on a periodic
[`PhaseSpaceGrid`](@ref). This is the default. It converges spectrally for smooth Wigner
functions and is the only discretisation that supports the exact Moyal operator.
"""
struct Spectral end

"""
    FiniteDifference(order::Integer = 4)

Central finite-difference discretisation on a periodic [`PhaseSpaceGrid`](@ref), with the
given even accuracy `order` for every derivative. The semi-discrete Wigner–Moyal operator is
then a sparse matrix, available as `sparse(op)`. It needs a truncated Moyal series
(`moyal_terms`), because the exact operator is nonlocal in momentum.
"""
struct FiniteDifference
    order::Int
    function FiniteDifference(order::Integer = 4)
        order >= 2 && iseven(order) || throw(
            ArgumentError(
                "finite-difference order must be even and at least 2, got $order",
            ),
        )
        return new(order)
    end
end

# Exact weights for offsets -r:r of the central difference approximating the
# `derivative`-th derivative to accuracy `order` on a unit grid, r = (derivative + order - 1) ÷ 2.
# They solve the moment equations Σₖ wₖ kʲ = derivative! δⱼ for j = 0, …, 2r in rationals.
function central_difference_weights(derivative::Integer, order::Integer)
    r = (derivative + order - 1) ÷ 2
    moments = [big(k)^j // 1 for j in 0:(2r), k in (-r):r]
    targets = [j == derivative ? factorial(big(j)) // 1 : big(0) // 1 for j in 0:(2r)]
    return moments \ targets
end

# Sparse `n × n` central-difference matrix for the `derivative`-th derivative on a periodic
# grid of spacing `h`.
function periodic_difference_matrix(
    derivative::Integer,
    order::Integer,
    n::Integer,
    h::Real,
)
    weights = central_difference_weights(derivative, order)
    r = length(weights) ÷ 2
    n >= 2r + 1 || throw(
        ArgumentError(
            "order $order differences of derivative $derivative need at least " *
            "$(2r + 1) points per axis, got $n",
        ),
    )
    rows, cols, values = Int[], Int[], Float64[]
    for (offset, weight) in zip((-r):r, weights)
        iszero(weight) && continue
        append!(rows, 1:n)
        append!(cols, mod1.((1:n) .+ offset, n))
        append!(values, fill(Float64(weight) / h^derivative, n))
    end
    return sparse(rows, cols, values, n, n)
end

# Angular wavenumbers of `rfft` along an axis of `n` points with spacing `h`. For even `n`
# the unpaired Nyquist wavenumber is set to zero, so that the odd operators built from it map
# real data to real data.
function wavenumbers(n::Integer, h::Real)
    κ = 2π .* collect(rfftfreq(n, 1 / h))
    iseven(n) && (κ[end] = 0)
    return κ
end
