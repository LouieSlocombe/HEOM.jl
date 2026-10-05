using LinearAlgebra: eigvals, istril, istriu

"""
    ExponentialBath(coefficients, rates; hbar = 1.0, diffusion = 0.0, counterterm = 0.0,
                    weights = ones(length(rates)), mixing = spzeros(length(rates), length(rates)))

Gaussian bath coupled linearly to position, `H_SB = q B`, with force correlation
`C(t) = ⟨B(t)B(0)⟩ = sum(coefficients[k] * exp(-rates[k]*t))` for `t ≥ 0`.
Rates must be real, finite and positive; coefficients may be complex and must be finite.
The two vectors must have equal length; empty vectors describe an isolated system.
Physical consistency of a user-supplied correlation is the caller's responsibility.

The optional real basis weights and mixing matrix generalise the correlation to
`C(t) = transpose(weights) * exp(-(Diagonal(rates) + mixing)*t) * coefficients`.
Mixing must have a zero diagonal, and the complete rate matrix must have eigenvalues
with positive real part. This supports polynomial exponentials at coincident poles,
and damped oscillatory correlations through real two-mode blocks. Default unit
weights and zero mixing give the independent exponential sum above. All supplied
arrays are copied; weights and mixing must be real and finite. Mixing is stored
sparsely, so a large Matsubara decomposition does not allocate a dense rate matrix.
The associated hierarchy follows the general basis construction of Ikeda and Scholes,
J. Chem. Phys. 152, 204101 (2020), DOI: 10.1063/5.0007327.

`diffusion ≥ 0` adds a Markovian remainder `diffusion * ∂p²W` on every hierarchy member.
`counterterm ≥ 0` adds `counterterm*q²` to the supplied system potential.
`hbar` must be finite and positive. Inputs are copied. See [`drude_lorentz_bath`](@ref)
for a thermal bath with its counterterm and Matsubara remainder included.
"""
struct ExponentialBath
    coefficients::Vector{ComplexF64}
    rates::Vector{Float64}
    weights::Vector{Float64}
    mixing::SparseMatrixCSC{Float64,Int}
    hbar::Float64
    diffusion::Float64
    counterterm::Float64
    function ExponentialBath(
        coefficients::AbstractVector,
        rates::AbstractVector;
        hbar::Real = 1.0,
        diffusion::Real = 0.0,
        counterterm::Real = 0.0,
        weights::AbstractVector = ones(length(rates)),
        mixing::AbstractMatrix = SparseArrays.spzeros(length(rates), length(rates)),
    )
        length(coefficients) == length(rates) ||
            throw(DimensionMismatch("bath coefficients and rates must have equal length"))
        all(x -> x isa Real, rates) || throw(ArgumentError("bath rates must be real"))
        length(weights) == length(rates) ||
            throw(DimensionMismatch("bath weights and rates must have equal length"))
        size(mixing) == (length(rates), length(rates)) || throw(
            DimensionMismatch("bath mixing must be a square matrix matching the rates"),
        )
        all(x -> x isa Real, weights) || throw(ArgumentError("bath weights must be real"))
        mixing_values = mixing isa SparseMatrixCSC ? mixing.nzval : mixing
        all(x -> x isa Real, mixing_values) ||
            throw(ArgumentError("bath mixing must be real"))
        c, ν = ComplexF64.(coefficients), Float64.(rates)
        w, M = Float64.(weights), copy(SparseMatrixCSC{Float64,Int}(sparse(mixing)))
        ħ, D, λ = Float64(hbar), Float64(diffusion), Float64(counterterm)
        all(isfinite, c) || throw(ArgumentError("bath coefficients must be finite"))
        all(x -> isfinite(x) && x > 0, ν) ||
            throw(ArgumentError("bath rates must be finite and positive"))
        all(isfinite, w) || throw(ArgumentError("bath weights must be finite"))
        all(isfinite, M.nzval) || throw(ArgumentError("bath mixing must be finite"))
        all(i -> iszero(M[i, i]), eachindex(ν)) || throw(
            ArgumentError("bath mixing must have a zero diagonal; use rates for damping"),
        )
        if !istriu(M) && !istril(M)
            rate_matrix = Matrix(M)
            for i in eachindex(ν)
                rate_matrix[i, i] = ν[i]
            end
            all(z -> isfinite(z) && real(z) > 0, eigvals(rate_matrix)) || throw(
                ArgumentError(
                    "the complete bath rate matrix must have positive-real eigenvalues",
                ),
            )
        end
        isfinite(ħ) && ħ > 0 || throw(ArgumentError("hbar must be finite and positive"))
        isfinite(D) && D >= 0 ||
            throw(ArgumentError("diffusion must be finite and nonnegative"))
        isfinite(λ) && λ >= 0 ||
            throw(ArgumentError("counterterm must be finite and nonnegative"))
        return new(c, ν, w, M, ħ, D, λ)
    end
end

# Evaluate Σ_{k=K+1}∞ 2x/(π²k²-x²) directly instead of subtracting retained
# poles from 1/x-cot(x). The latter loses accuracy near a Drude–Matsubara
# collision, even when that pole is retained and the omitted tail is smooth.
function matsubara_tail(x, K)
    a = x / π
    k = Float64(K) + 1
    # Accurate argument reduction also protects k-a when the *first omitted*
    # pole is close. Direct subtraction of x/π would lose its small separation.
    offset = rem2pi(2x, RoundNearest) / (2π)
    nearest = round(a - offset)
    left, right = (k - nearest) - offset, (k + nearest) + offset
    tail = 0.0
    # Move the Euler–Maclaurin endpoint away from the nearest pole. At most
    # sixteen positive terms are needed because the caller requires K+1 > a.
    while left < 16
        tail += (2a / right) / left
        left += 1
        right += 1
    end
    # Apply Euler–Maclaurin to 1/(k-a)-1/(k+a). Its integral is this logarithm;
    # log1p/expm1 keep both it and all derivative differences accurate as a → 0.
    logratio = log1p(2a / left)
    tail += logratio + (a / right) / left
    inverse = inv(left)
    for (n, coefficient) in
        enumerate((1 / 12, -1 / 120, 1 / 252, -1 / 240, 1 / 132, -691 / 32760))
        tail += coefficient * inverse^(2n) * (-expm1(-2n * logratio))
    end
    return tail / π
end

# Analytic continuation of cot(x)-1/x through x=0. Evaluating the two terms
# separately loses the finite part of a colliding Drude/Matsubara residue.
function regularized_cot(x)
    return -x * (1 / 3 + x^2 * (1 / 45 + x^2 * (2 / 945 + x^2 / 4725)))
end

"""
    drude_lorentz_bath(; reorganization, cutoff, kT, matsubara = 0, hbar = 1.0,
                       terminator = true)

Thermal [`ExponentialBath`](@ref) with Drude spectral density
`J(ω) = 2λγω/(ω² + γ²)`, where `λ = reorganization` and `γ = cutoff`.
The convention is `C(t) = (ħ/π)∫₀∞ J(ω)[coth(ħω/(2kT))*cos(ωt) - i*sin(ωt)]dω`.
For dimensional position coupling, `λ` has units of energy / position².

The retained correlation coefficients are

    ν₀ = γ,                 c₀ = λħγ [cot(ħγ/(2kT)) - i],
    νₖ = 2πk*kT/ħ,          cₖ = 4λγ*kT*νₖ/(νₖ² - γ²),  k = 1,…,matsubara.

The system counterterm is `λq²`. With `terminator = true`, the omitted fast Matsubara
terms are approximated by momentum diffusion
`D = 2λ*kT/γ - weights' * (Γ \\ real(coefficients))`, where
`Γ = Diagonal(rates) + mixing`; away from collisions this reduces to
`2λ*kT/γ - sum(real(cₖ)/νₖ)`. This is a bath-correlation tail approximation,
not a closure of the hierarchy depth. Include enough Matsubara terms that the first
omitted rate exceeds `γ` and is fast compared with the relevant system frequencies.
Increase both `matsubara` and hierarchy `depth` to check
convergence, especially at low temperature. `terminator = false` simply drops the tail.

`λ` must be finite and nonnegative; `γ`, `kT` and `ħ` finite and positive. Coincident
or nearly coincident Drude and Matsubara poles use a coupled correlation basis that
includes their finite polynomial-exponential limit. The colliding Matsubara pole must
be retained. Zero coupling returns an empty bath. The zero-temperature limit requires
another correlation decomposition.
"""
function drude_lorentz_bath(;
    reorganization::Real,
    cutoff::Real,
    kT::Real,
    matsubara::Integer = 0,
    hbar::Real = 1.0,
    terminator::Bool = true,
)
    λ, γ, θ, ħ = Float64(reorganization), Float64(cutoff), Float64(kT), Float64(hbar)
    isfinite(λ) && λ >= 0 ||
        throw(ArgumentError("reorganization must be finite and nonnegative"))
    isfinite(γ) && γ > 0 || throw(ArgumentError("cutoff must be finite and positive"))
    isfinite(θ) && θ > 0 || throw(ArgumentError("kT must be finite and positive"))
    isfinite(ħ) && ħ > 0 || throw(ArgumentError("hbar must be finite and positive"))
    matsubara >= 0 || throw(ArgumentError("matsubara must be nonnegative"))
    iszero(λ) && return ExponentialBath(ComplexF64[], Float64[]; hbar = ħ)
    x = ħ * γ / (2θ)
    isfinite(x) && x > 0 ||
        throw(ArgumentError("ħ*cutoff/(2kT) must be finite and positive"))
    # Use a divided-difference basis near a repeated pole, before either large
    # residue is evaluated. The regularized cotangent uses accurate reduction.
    pole = round(x / π)
    offset = rem2pi(2x, RoundNearest) / 2
    nearby = pole > 0 && abs(offset) < 1e-3
    nearby &&
        pole > matsubara &&
        abs(offset) <= 32eps(Float64) * max(1, x) &&
        throw(
            ArgumentError(
                "increase matsubara to retain the coincident Drude/Matsubara pole",
            ),
        )
    collision = nearby && pole <= matsubara
    resonant = collision ? Int(pole) + 1 : 0
    ν = [γ; [2π * k * θ / ħ for k in 1:matsubara]]
    c = zeros(ComplexF64, length(ν))
    weights = ones(length(ν))
    mixing = SparseArrays.spzeros(length(ν), length(ν))
    if collision
        # c₀+cₖ = λħγ[cot(δ)-1/δ+1/(x+kπ)], δ=x-kπ.
        # At coincidence this is λ*kT-iλħγ, while (νₖ-γ)cₖ=2λγ*kT.
        c[1] = ComplexF64(
            λ * ħ * γ * (regularized_cot(offset) + inv(x + pole * π)),
            -λ * ħ * γ,
        )
        v = ν[resonant]
        c[resonant] = 4λ * γ * θ * (v / (v + γ))
        weights[resonant] = 0
        mixing[1, resonant] = 1
    else
        c[1] = ComplexF64(λ * ħ * γ * cot(x), -λ * ħ * γ)
    end
    for k in 2:length(ν)
        k == resonant && continue
        ratio = γ / ν[k]
        c[k] = 4λ * θ * ratio / ((1 - ratio) * (1 + ratio))
    end
    D = 0.0
    if terminator
        2π * (matsubara + 1) * θ / ħ > γ || throw(
            ArgumentError("increase matsubara so the first omitted rate exceeds cutoff"),
        )
        D = λ * ħ * matsubara_tail(x, matsubara)
    end
    return ExponentialBath(c, ν; hbar = ħ, diffusion = D, counterterm = λ, weights, mixing)
end

struct HEOMMomentum{D1,D2}
    first::D1
    second::D2
end

struct WignerHEOM{H<:AbstractWignerMoyal,D} <: AbstractPhaseSpaceOperator
    grid::PhaseSpaceGrid
    mass::Float64
    hbar::Float64
    bath::ExponentialBath
    depth::Int
    scaled::Bool
    indices::Vector{Vector{Int}}
    upper::Matrix{Int}
    lower::Matrix{Int}
    damping::Vector{Float64}
    coordinate_coupling::Vector{Float64}
    lowering::Matrix{Float64}
    raising_derivative::Matrix{Float64}
    raising_coordinate::Matrix{Float64}
    log_scales::Vector{Float64}
    transfer::SparseMatrixCSC{Float64,Int}
    hamiltonian::H
    momentum::D
    scratch::Matrix{Float64}
end

"""
    hierarchy_size(modes, depth)

Exact number of auxiliary members, including the physical root, for `modes` bath
exponentials and total occupation at most `depth`. Returns a `BigInt` so resource
estimates do not overflow before a hierarchy is built. One Float64 state requires
`8 * prod(size(grid)) * hierarchy_size(modes, depth)` bytes; integration needs several
state-sized work arrays, and saved trajectories need additional storage.
"""
function hierarchy_size(modes::Integer, depth::Integer)
    modes >= 0 || throw(ArgumentError("modes must be nonnegative"))
    depth >= 0 || throw(ArgumentError("depth must be nonnegative"))
    return binomial(big(modes) + big(depth), big(depth))
end

# Breadth-first traversal gives graded indices without enumerating a rectangular box.
function build_hierarchy(modes, depth)
    indices = [zeros(Int, modes)]
    lookup = Dict(Tuple(first(indices)) => 1)
    a = 1
    while a <= length(indices)
        n = indices[a]
        if sum(n) < depth
            for k in 1:modes
                next = copy(n)
                next[k] += 1
                key = Tuple(next)
                if !haskey(lookup, key)
                    push!(indices, next)
                    lookup[key] = length(indices)
                end
            end
        end
        a += 1
    end
    upper, lower = zeros(Int, modes, length(indices)), zeros(Int, modes, length(indices))
    for (a, n) in enumerate(indices), k in 1:modes
        next = copy(n)
        next[k] += 1
        upper[k, a] = get(lookup, Tuple(next), 0)
        next[k] -= 2
        lower[k, a] = get(lookup, Tuple(next), 0)
    end
    return indices, upper, lower
end

# sqrt(abs(c)) without overflowing abs(c) for finite complex components.
function coefficient_scale(c)
    largest = max(abs(real(c)), abs(imag(c)))
    iszero(largest) && return 1.0
    return sqrt(largest) * sqrt(hypot(real(c) / largest, imag(c) / largest))
end

function hierarchy_couplings(bath, indices, upper, lower, scaled)
    scales = coefficient_scale.(bath.coefficients)
    denominators = scaled ? scales : ones(length(scales))
    cr = real.(bath.coefficients) ./ denominators
    ci = 2 .* ((imag.(bath.coefficients) ./ denominators) ./ bath.hbar)
    lowering = zeros(size(lower))
    raising_derivative, raising_coordinate = zeros(size(upper)), zeros(size(upper))
    log_scales = zeros(length(indices))
    for (a, n) in enumerate(indices)
        # Breadth-first order guarantees that the selected parent was already visited.
        k = findfirst(!iszero, n)
        if k !== nothing
            log_scales[a] = log_scales[lower[k, a]] + log(n[k]) / 2 + log(scales[k])
        end
        for k in eachindex(scales)
            if lower[k, a] != 0
                lowering[k, a] = bath.weights[k] * (scaled ? sqrt(n[k]) * scales[k] : 1.0)
            end
            if upper[k, a] != 0
                factor = scaled ? sqrt(n[k] + 1) : n[k] + 1
                raising_derivative[k, a] = factor * cr[k]
                raising_coordinate[k, a] = factor * ci[k]
            end
        end
    end
    all(isfinite, ci) &&
    all(isfinite, lowering) &&
    all(isfinite, raising_derivative) &&
    all(isfinite, raising_coordinate) &&
    all(isfinite, log_scales) ||
        throw(ArgumentError("hierarchy coupling factors must be finite"))
    return ci, lowering, raising_derivative, raising_coordinate, log_scales
end

function hierarchy_transfer(bath, indices, upper, lower, scaled)
    members = length(indices)
    rows, columns, values = Int[], Int[], Float64[]
    all(iszero, bath.mixing.nzval) && return sparse(rows, columns, values, members, members)
    scales = coefficient_scale.(bath.coefficients)
    # G[k,j] contributes -n_k G[k,j] U_(n-e_k+e_j) to the target n.
    # A source loses an occupation in j and gains one in k at the same total tier.
    for (source, n) in enumerate(indices), j in eachindex(bath.rates)
        parent = lower[j, source]
        iszero(parent) && continue
        for pointer in SparseArrays.nzrange(bath.mixing, j)
            k = bath.mixing.rowval[pointer]
            mixing = bath.mixing.nzval[pointer]
            iszero(mixing) && continue
            target = upper[k, parent]
            factor = if scaled
                candidate = -mixing * sqrt(n[k] + 1) * sqrt(n[j]) * (scales[j] / scales[k])
                if isfinite(candidate) && !iszero(candidate)
                    candidate
                else
                    magnitude =
                        log(abs(mixing)) +
                        (log(n[k] + 1) + log(n[j])) / 2 +
                        log(scales[j]) - log(scales[k])
                    copysign(exp(magnitude), -mixing)
                end
            else
                -(n[k] + 1) * mixing
            end
            isfinite(factor) ||
                throw(ArgumentError("hierarchy mixing factors must be finite"))
            push!(rows, target)
            push!(columns, source)
            push!(values, factor)
        end
    end
    return sparse(rows, columns, values, members, members)
end

function heom_momentum(op::SpectralWignerMoyal)
    np = length(op.grid.p)
    κ = 2π .* collect(rfftfreq(np, 1 / op.grid.dp))
    return HEOMMomentum(
        reshape(im .* wavenumbers(np, op.grid.dp) ./ np, 1, :),
        reshape(-κ .^ 2 ./ np, 1, :),
    )
end

function heom_momentum(op::FiniteDifferenceWignerMoyal)
    np = length(op.grid.p)
    return HEOMMomentum(
        periodic_difference_matrix(1, op.order, np, op.grid.dp),
        periodic_difference_matrix(2, op.order, np, op.grid.dp),
    )
end

"""
    heom_operator(grid; mass, potential, bath, depth, discretization = Spectral(),
                  moyal_terms = nothing, scaled = false, max_ados = nothing)

Wigner–Moyal hierarchical equations of motion for linear coordinate coupling to an
[`ExponentialBath`](@ref). The unscaled, real auxiliary Wigner functions obey

    dWₙ/dt = [L_WM(V + counterterm*q²) - Σₖ nₖνₖ + D ∂p²] Wₙ
             + Σₖ weights[k] ∂p Wₙ₊ₑₖ
             + Σₖ nₖ [real(cₖ) ∂p + 2imag(cₖ)q/ħ] Wₙ₋ₑₖ
             - Σₖⱼ nₖ mixing[k,j] Wₙ₋ₑₖ₊ₑⱼ.

This is the Wigner transform of the coordinate-coupled density-operator HEOM:
`[q,ρ]_W = iħ∂pW` and `{q,ρ}_W = 2qW`. The physical state is `W₀`.
`potential` is the physical system potential; the bath counterterm is added internally.
`hbar` comes from `bath`. Both Hamiltonian discretisations and their Moyal options are
supported (finite differences require integer `moyal_terms`).

All nonnegative multi-indices with `sum(n) ≤ depth` are retained, with the root first.
There are `binomial(length(bath.rates) + depth, depth)` members. Upward terms beyond
`depth` are set to zero (hard cutoff); `depth = 0` omits all explicit bath memory.
Converge the depth, bath decomposition, box and resolution independently.

`scaled = true` propagates `W̃ₙ = Wₙ / sqrt(∏ₖ nₖ! aₖ^nₖ)`, where `aₖ = abs(cₖ)`
for nonzero coefficients and `aₖ = 1` for zero coefficients. The physical root is
unchanged. This is the amplitude/factorial scaling of Shi et al., J. Chem. Phys. 130,
084105 (2009), adapted to Wigner functions. It improves numerical conditioning for
deep hierarchies and strong coupling; it does not change the retained equations or
guarantee the accuracy or stability of a finite truncation. No factorials or powers
are formed during construction. Use [`rescale_hierarchy`](@ref) to convert an existing
3D initial state or restart between conventions.

`max_ados` optionally limits the number of members before allocating the hierarchy.
Use [`hierarchy_size`](@ref) to estimate storage beforehand. No positivity correction
is applied. The operator owns mutable work buffers; give each concurrent solve its
own operator. See [`heom_problem`](@ref).
"""
function heom_operator(
    grid::PhaseSpaceGrid;
    mass::Real,
    potential,
    bath::ExponentialBath,
    depth::Integer,
    discretization::Union{Spectral,FiniteDifference} = Spectral(),
    moyal_terms::Union{Nothing,Integer} = nothing,
    scaled::Bool = false,
    max_ados::Union{Nothing,Integer} = nothing,
)
    m = Float64(mass)
    isfinite(m) && m > 0 || throw(ArgumentError("mass must be finite and positive"))
    depth >= 0 || throw(ArgumentError("depth must be nonnegative"))
    depth <= typemax(Int) ||
        throw(ArgumentError("depth exceeds the supported integer range"))
    members = hierarchy_size(length(bath.rates), depth)
    members <= typemax(Int) ||
        throw(ArgumentError("hierarchy size exceeds the integer range"))
    prod(big.(size(grid))) * members <= typemax(Int) ÷ sizeof(Float64) ||
        throw(ArgumentError("one hierarchy state exceeds the addressable array size"))
    if max_ados !== nothing
        max_ados > 0 || throw(ArgumentError("max_ados must be positive"))
        members <= max_ados || throw(
            ArgumentError(
                "hierarchy requires $members members, exceeding max_ados=$max_ados",
            ),
        )
    end
    # Copy and revalidate the bath, since its input vectors remain user-accessible.
    b = ExponentialBath(
        bath.coefficients,
        bath.rates;
        hbar = bath.hbar,
        diffusion = bath.diffusion,
        counterterm = bath.counterterm,
        weights = bath.weights,
        mixing = bath.mixing,
    )
    V(q) = potential(q) + b.counterterm * q^2
    h = wigner_moyal_operator(
        grid;
        mass = m,
        potential = V,
        hbar = b.hbar,
        discretization,
        moyal_terms,
    )
    indices, upper, lower = build_hierarchy(length(b.rates), Int(depth))
    damping = [sum(n .* b.rates) for n in indices]
    all(isfinite, damping) || throw(ArgumentError("hierarchy damping rates must be finite"))
    coordinate_coupling, lowering, raising_derivative, raising_coordinate, log_scales =
        hierarchy_couplings(b, indices, upper, lower, scaled)
    transfer = hierarchy_transfer(b, indices, upper, lower, scaled)
    return WignerHEOM(
        grid,
        m,
        b.hbar,
        b,
        Int(depth),
        scaled,
        indices,
        upper,
        lower,
        damping,
        coordinate_coupling,
        lowering,
        raising_derivative,
        raising_coordinate,
        log_scales,
        transfer,
        h,
        heom_momentum(h),
        zeros(size(grid)),
    )
end

"""
    hierarchy_indices(op)

Copy of the multi-indices in the order used by the third axis of a HEOM state.
The first is the all-zero physical member; subsequent entries are ordered by tier.
"""
hierarchy_indices(op::WignerHEOM) = deepcopy(op.indices)

"""
    rescale_hierarchy(U, op; scaled)

Copy a real 3D hierarchy in `op`'s convention into the requested scaling convention.
The grid, bath and index ordering are unchanged, and the physical root is copied
exactly. For example, convert an unscaled restart `U` using
`rescale_hierarchy(U, unscaled_op; scaled = true)` before passing it to a matching
scaled operator. To convert back, supply that scaled operator with `scaled = false`.

Conversion uses logarithmic scales and never forms factorials. Throws `ArgumentError`
if a nonzero converted value would overflow or underflow Float64; propagation in the
scaled convention itself does not require representable unscaled auxiliaries.
"""
function rescale_hierarchy(U::AbstractArray, op::WignerHEOM; scaled::Bool)
    check_hierarchy_size(U, op)
    eltype(U) <: Real || throw(ArgumentError("hierarchy must be real"))
    all(isfinite, U) || throw(ArgumentError("hierarchy must be finite"))
    out = Array{Float64,3}(U)
    all(isfinite, out) || throw(ArgumentError("hierarchy must be finite in Float64"))
    scaled == op.scaled && return out
    for a in 2:length(op.indices)
        exponent = (scaled ? -1 : 1) * op.log_scales[a]
        factor = exp(exponent)
        W = @view out[:, :, a]
        for i in eachindex(W)
            x = W[i]
            iszero(x) && continue
            y =
                isfinite(factor) && !iszero(factor) ? x * factor :
                copysign(exp(log(abs(x)) + exponent), x)
            isfinite(y) && !iszero(y) ||
                throw(ArgumentError("converted auxiliary $a is outside the Float64 range"))
            W[i] = y
        end
    end
    return out
end

function heom_derivative!(out, W, h::SpectralWignerMoyal, symbol)
    h.W .= W
    mul!(h.Wp, h.forward_p, h.W)
    h.Wp .*= symbol
    mul!(out, h.backward_p, h.Wp)
    return nothing
end

function heom_derivative!(out, W, h::FiniteDifferenceWignerMoyal, matrix)
    mul!(out, W, transpose(matrix))
    return nothing
end

function check_hierarchy_size(U, op)
    expected = (size(op.grid)..., length(op.indices))
    size(U) == expected ||
        throw(DimensionMismatch("hierarchy must have size $expected, got $(size(U))"))
    return nothing
end

"""
    heom!(dU, U, op, t)

Evaluate the Wigner HEOM in place. `U` and `dU` have shape `(nq, np, number_of_ADOs)`
and must not alias. The time is unused. Only the root is a normalised Wigner function;
auxiliaries encode bath correlations and may have nonzero integrals.
"""
function heom!(dU, U, op::WignerHEOM, t)
    check_hierarchy_size(U, op)
    check_hierarchy_size(dU, op)
    h, b = op.hamiltonian, op.bath
    for a in eachindex(op.indices)
        W, dW = @view(U[:, :, a]), @view(dU[:, :, a])
        wigner_moyal!(dW, W, h, t)
        @. dW -= op.damping[a] * W
        if !iszero(b.diffusion)
            heom_derivative!(op.scratch, W, h, op.momentum.second)
            @. dW += b.diffusion * op.scratch
        end
    end
    # Differentiate each source once, then scatter its upward/downward contributions.
    for a in eachindex(op.indices)
        W = @view U[:, :, a]
        for pointer in SparseArrays.nzrange(op.transfer, a)
            target = op.transfer.rowval[pointer]
            factor = op.transfer.nzval[pointer]
            dW = @view dU[:, :, target]
            @. dW += factor * W
        end
        heom_derivative!(op.scratch, W, h, op.momentum.first)
        for k in eachindex(b.rates)
            below, above = op.lower[k, a], op.upper[k, a]
            if below != 0
                dW = @view dU[:, :, below]
                factor = op.lowering[k, a]
                @. dW += factor * op.scratch
            end
            if above != 0
                dW = @view dU[:, :, above]
                cr, ci = op.raising_derivative[k, a], op.raising_coordinate[k, a]
                @. dW += cr * op.scratch + ci * op.grid.q * W
            end
        end
    end
    return nothing
end

"""
    heom_problem(W0, tspan, grid; mass, potential, bath, depth, jacobian = :matrixfree, kwargs...)
    heom_problem(W0, tspan, op; jacobian = :matrixfree)

Construct an in-place SciML `ODEProblem` for [`heom_operator`](@ref). An `nq × np`
Wigner matrix initialises the physical root and sets every auxiliary to zero: a
factorized initial system and bare bath state. With a counterterm this includes the
physical initial-slip transient. An entire 3D hierarchy may instead be supplied to
restart a calculation or initialise a correlated state. The input is copied to a
real, finite `Float64` array, without renormalisation.

Use `physical_wigner(sol)` to extract the final root. Explicit solvers such as `Vern7()`
work; high bath rates and large depth can make the hierarchy stiff. The default
`jacobian = :matrixfree` supplies the exact Jacobian-vector product and time derivative,
including for spectral operators. For example, `Rodas5P(autodiff = AutoFiniteDiff(), linsolve = KrylovJL_GMRES(), concrete_jac = false)` integrates without constructing a
dense Jacobian or differentiating FFTW buffers. These solver types are supplied by
OrdinaryDiffEqRosenbrock, ADTypes and LinearSolve, respectively.

`jacobian = :sparse` assembles the constant finite-difference generator and provides
an exact sparse Jacobian; it requires `FiniteDifference` and can use `Rodas5P()`.
It stores a matrix in addition to the hierarchy. Rebuild the problem rather than
replacing its operator parameter when using this cached sparse Jacobian.
"""
function heom_problem(
    U0::AbstractArray,
    tspan,
    grid::PhaseSpaceGrid;
    jacobian::Symbol = :matrixfree,
    kwargs...,
)
    return heom_problem(U0, tspan, heom_operator(grid; kwargs...); jacobian)
end

function heom_problem(
    U0::AbstractArray,
    tspan,
    op::WignerHEOM;
    jacobian::Symbol = :matrixfree,
)
    eltype(U0) <: Real || throw(ArgumentError("initial hierarchy must be real"))
    all(isfinite, U0) || throw(ArgumentError("initial hierarchy must be finite"))
    if ndims(U0) == 2
        check_size(U0, op.grid)
        U = zeros(size(op.grid)..., length(op.indices))
        U[:, :, 1] .= U0
    else
        check_hierarchy_size(U0, op)
        U = Array{Float64,3}(U0)
    end
    all(isfinite, U) || throw(ArgumentError("initial hierarchy must be finite in Float64"))
    return ODEProblem(heom_ode_function(op, jacobian), U, tspan, op)
end

"""
    physical_wigner(U::AbstractArray{<:Real,3})
    physical_wigner(sol[, index])

View of the physical (first) hierarchy member, usable with all matrix-and-grid
observables and plotting functions. For a HEOM solution, select a saved state with
`index` (the final state by default). Mutating the view changes the underlying state.
"""
function physical_wigner(U::AbstractArray{<:Real,3})
    size(U, 3) > 0 || throw(ArgumentError("hierarchy needs a physical member"))
    return @view U[:, :, 1]
end

function physical_wigner(sol::AbstractODESolution, index::Integer = lastindex(sol.u))
    sol.prob.p isa WignerHEOM ||
        throw(ArgumentError("solution must come from heom_problem"))
    return physical_wigner(sol.u[index])
end

# The existing solution diagnostics/plots consume physical states only.
physical_states(states, ::AbstractPhaseSpaceOperator) = states
physical_states(states, ::WignerHEOM) = physical_wigner.(states)
physical_state(W, ::AbstractPhaseSpaceOperator) = W
physical_state(U, ::WignerHEOM) = physical_wigner(U)

function Base.show(io::IO, op::WignerHEOM)
    print(
        io,
        "WignerHEOM(",
        op.grid,
        ", depth = ",
        op.depth,
        ", modes = ",
        length(op.bath.rates),
        ", members = ",
        length(op.indices),
        ", scaled = ",
        op.scaled,
        ")",
    )
    return nothing
end
