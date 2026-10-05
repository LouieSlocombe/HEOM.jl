"""
    wavefunction_wigner(psi, grid::PhaseSpaceGrid; hbar = 1.0)

Transform a position wavefunction into a real `nq × np` Wigner matrix on `grid`.
`psi` can be a callable `psi(q)` or a vector containing its values at `grid.q`.
The convention is

    B(q, θ) = ψ(q - ħθ/2) conj(ψ(q + ħθ/2)),
    W(q, p) = (1/2π) ∫ B(q, θ) exp(i p θ) dθ.

The discrete transform samples `θ = 2π k / (np * grid.dp)` and applies an inverse
FFT with the physical scaling `1 / grid.dp`. A momentum-origin phase evaluates
the result at `grid.p`, including grids that are not centred at zero. For even
`np`, the two opposite Nyquist samples are averaged after applying that phase;
this is the real, symmetric trapezoidal convention at the transform endpoints.

Sampled wavefunctions use local cubic Lagrange interpolation (degree `nq - 1`
when `nq < 4`) and are zero outside `[first(grid.q), last(grid.q)]`. They must be
resolved on the position grid and negligible at its edges. Refining the grid
reduces the interpolation error; no periodic copies of the wavefunction are
introduced. Callable wavefunctions are instead evaluated directly at all shifted
positions, including positions outside the grid. The momentum interval must
resolve the state and its tails to avoid momentum aliasing.

The state is **not** normalised: `position_density(W, grid)` equals
`abs2.(psi.(grid.q))` (or `abs2.(psi)` for samples), to roundoff. Negative Wigner
values are retained. Wavefunction values must be finite, and `hbar` finite and
positive. The output can be passed directly to the phase-space problem builders.
"""
function wavefunction_wigner(psi, grid::PhaseSpaceGrid; hbar::Real = 1.0)
    evaluate = _wavefunction_evaluator(psi, grid)
    kernel = (x, y) -> evaluate(x) * conj(evaluate(y))
    return _kernel_wigner(kernel, grid, hbar)
end

"""
    density_matrix_wigner(rho, grid::PhaseSpaceGrid; hbar = 1.0)

Transform a position-space density matrix into a real `nq × np` Wigner matrix.
`rho` is either a callable kernel `rho(x, y) = ⟨x|ρ|y⟩` or an `nq × nq` matrix of
that kernel's values at `grid.q`. Thus a pure state has
`rho[i, j] = psi[i] * conj(psi[j])`, and a sampled kernel's trace is
`grid.dq * sum(real(rho[i, i]) for i in eachindex(grid.q))`. Matrix entries are
kernel samples, not coefficients in an orthonormal discrete basis.

The transform uses `B(q, θ) = rho(q - ħθ/2, q + ħθ/2)` with the same Fourier
scaling, momentum-origin phase and even-grid Nyquist convention as
[`wavefunction_wigner`](@ref). Sampled matrices use tensor-product local cubic
Lagrange interpolation (degree `nq - 1` for `nq < 4`) and zero outside the sampled
position interval in either argument. Resolve the kernel and make it negligible
at the position boundaries. A callable kernel must be defined at the shifted
positions, which can lie outside that interval.

All values must be finite. Sampled matrices are checked for Hermiticity over the
whole matrix; callable kernels are checked at every pair used by the transform.
Hermiticity errors larger than `1e-12` times the largest sampled real or imaginary
component raise
an `ArgumentError`; smaller roundoff errors are symmetrised before transforming.
Positive semidefiniteness is not checked. The transform preserves the supplied
trace and negative Wigner values without normalising or clipping them.
`hbar` must be finite and positive.
"""
function density_matrix_wigner(rho, grid::PhaseSpaceGrid; hbar::Real = 1.0)
    return _kernel_wigner(_density_kernel_evaluator(rho, grid), grid, hbar)
end

# Convert at the evaluation boundary so invalid complex values cannot be hidden
# by taking real parts or by the Hermitian completion used for the FFT.
function _finite_state_value(value)
    value isa Number && isfinite(value) ||
        throw(ArgumentError("state values must be finite numbers"))
    converted = ComplexF64(value)
    isfinite(converted) ||
        throw(ArgumentError("state values must remain finite in ComplexF64"))
    return converted
end

_wavefunction_evaluator(psi, ::PhaseSpaceGrid) = x -> _finite_state_value(psi(x))

function _wavefunction_evaluator(psi::AbstractVector, grid::PhaseSpaceGrid)
    length(psi) == length(grid.q) ||
        throw(DimensionMismatch("wavefunction samples must have length $(length(grid.q))"))
    values = _finite_state_value.(psi)
    return x -> _interpolate_wavefunction(values, grid, x)
end

_density_kernel_evaluator(rho, ::PhaseSpaceGrid) = rho

function _density_kernel_evaluator(rho::AbstractMatrix, grid::PhaseSpaceGrid)
    nq = length(grid.q)
    size(rho) == (nq, nq) ||
        throw(DimensionMismatch("density kernel samples must have size ($nq, $nq)"))
    values = _finite_state_value.(rho)
    scale = maximum(_state_maxabs, values)
    maximum(_state_maxabs, values .- values') <= 1e-12 * scale ||
        throw(ArgumentError("density kernel must be Hermitian"))
    return (x, y) -> _interpolate_density_kernel(values, grid, x, y)
end

_state_maxabs(value) = max(abs(real(value)), abs(imag(value)))

# Halve the factor that cannot underflow first, avoiding an intermediate overflow
# when ħ*θ is too large but the requested half-displacement remains representable.
_wigner_shift(ħ, theta) = ħ >= 1 ? (ħ / 2) * theta : ħ * (theta / 2)

# A fixed-size tuple avoids allocating a tiny weights vector at every shifted
# sample. Smaller grids use their available two or three interpolation nodes.
function _state_interpolation_stencil(grid::PhaseSpaceGrid, x)
    if x < first(grid.q) || x > last(grid.q)
        return 0, (0.0, 0.0, 0.0, 0.0)
    end
    n = length(grid.q)
    m = min(4, n)
    u = (x - first(grid.q)) / grid.dq
    firstindex = clamp(floor(Int, u) + 1 - fld(m - 1, 2), 1, n - m + 1)
    t = u - (firstindex - 1)
    weights = ntuple(4) do k
        k > m && return 0.0
        weight = 1.0
        for l in 1:m
            if l != k
                weight *= (t - (l - 1)) / (k - l)
            end
        end
        return weight
    end
    return firstindex, weights
end

function _interpolate_wavefunction(values, grid, x)
    firstindex, weights = _state_interpolation_stencil(grid, x)
    firstindex == 0 && return 0.0 + 0.0im
    value = 0.0 + 0.0im
    for k in 1:min(4, length(grid.q))
        value += weights[k] * values[firstindex+k-1]
    end
    return value
end

function _interpolate_density_kernel(values, grid, x, y)
    firstx, wx = _state_interpolation_stencil(grid, x)
    firsty, wy = _state_interpolation_stencil(grid, y)
    (firstx == 0 || firsty == 0) && return 0.0 + 0.0im
    value = 0.0 + 0.0im
    for l in 1:min(4, length(grid.q)), k in 1:min(4, length(grid.q))
        value += wx[k] * wy[l] * values[firstx+k-1, firsty+l-1]
    end
    return value
end

function _kernel_wigner(kernel, grid::PhaseSpaceGrid, hbar::Real)
    ħ = Float64(hbar)
    isfinite(ħ) && ħ > 0 ||
        throw(ArgumentError("hbar must be finite and positive in Float64"))
    nq, np = size(grid)
    dtheta = (2π / np) / grid.dp
    max_theta = (np ÷ 2) * dtheta
    max_shift = _wigner_shift(ħ, max_theta)
    isfinite(dtheta) &&
    dtheta > 0 &&
    isfinite(max_theta) &&
    isfinite(max_shift) &&
    isfinite(first(grid.q) - max_shift) &&
    isfinite(last(grid.q) + max_shift) &&
    isfinite(first(grid.p) * max_theta) ||
        throw(ArgumentError("the grid and hbar must give finite transform coordinates"))
    transformed = Matrix{ComplexF64}(undef, nq, np)
    scale = 0.0
    hermitian_error = 0.0
    for k in 0:(np÷2)
        theta = k * dtheta
        shift = _wigner_shift(ħ, theta)
        phase = cis(first(grid.p) * theta)
        for i in eachindex(grid.q)
            q = grid.q[i]
            forward = _finite_state_value(kernel(q - shift, q + shift))
            backward = _finite_state_value(kernel(q + shift, q - shift))
            scale = max(scale, _state_maxabs(forward), _state_maxabs(backward))
            hermitian_error = max(hermitian_error, _state_maxabs(forward - conj(backward)))
            value = _finite_state_value((forward / 2 + conj(backward) / 2) * phase)
            if k == 0 || (iseven(np) && k == np ÷ 2)
                # Zero frequency and the shared ±Nyquist bin must be real. At
                # Nyquist this averages the independently checked endpoint pair.
                transformed[i, k+1] = real(value)
            else
                transformed[i, k+1] = value
                transformed[i, np-k+1] = conj(value)
            end
        end
    end
    hermitian_error <= 1e-12 * scale ||
        throw(ArgumentError("density kernel must be Hermitian at the transform samples"))
    W = real.(ifft(transformed, 2)) ./ grid.dp
    all(isfinite, W) || throw(ArgumentError("Wigner transform must remain finite"))
    return W
end
