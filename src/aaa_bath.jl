using LinearAlgebra: diagm, eigvals, norm, svd

# A sampled spectral density must be a physical, finite J(ω) ≥ 0 for ω > 0.
function spectral_density_value(spectral_density, ω)
    value = spectral_density(ω)
    value isa Real && isfinite(value) && value >= 0 || throw(
        ArgumentError(
            "the spectral density must be real, finite and nonnegative; " *
            "J($ω) = $value",
        ),
    )
    return Float64(value)
end

# (1/π)∫₀^∞ J(ω)/ω dω by exp-sinh quadrature. With ω = scale*exp((π/2)sinh(t)),
# J(ω)/ω dω = (π/2)J(ω)cosh(t) dt decays double exponentially at both ends whenever J
# behaves algebraically at zero and infinity. The window |t| ≤ 4.5 spans scale*10^±30.
function reorganization_integral(spectral_density, scale)
    term(t) =
        spectral_density_value(spectral_density, scale * exp(π / 2 * sinh(t))) * cosh(t)
    edges = max(term(-4.5), term(4.5))
    h = 0.5
    total = h * sum(term, -4.5:h:4.5)
    for _ in 1:12
        h /= 2
        refined = total / 2 + h * sum(term, (-4.5+h):(2h):4.5)
        converged = abs(refined - total) <= 1e-10 * refined && edges <= 1e-12 * refined
        total = refined
        converged && return total / 2
    end
    throw(
        ArgumentError(
            "(1/π)∫J(ω)/ω dω did not converge over frequencies $(scale)*10^±30; " *
            "check that J(ω)/ω is integrable or pass reorganization explicitly",
        ),
    )
end

# Zeros of the barycentric denominator Σⱼ wⱼ/(u - zⱼ) are the finite eigenvalues of an
# arrowhead pencil (Nakatsukasa, Sète and Trefethen, SIAM J. Sci. Comput. 40, A1494
# (2018)). The pencil loses relative accuracy for poles much smaller than the largest
# support point, so Newton steps on the well-conditioned barycentric sum refine them.
function barycentric_poles(support, weights)
    m = length(support)
    pencil = diagm(0 => [0; support])
    pencil[1, 2:end] = weights
    pencil[2:end, 1] .= 1
    poles = ComplexF64.(filter(isfinite, eigvals(pencil, diagm(0 => [0; ones(m)]))))
    for k in eachindex(poles)
        for _ in 1:10
            r = inv.(poles[k] .- support)
            step = -sum(weights .* r) / sum(weights .* r .^ 2)
            poles[k] -= step
            abs(step) <= 4eps() * abs(poles[k]) && break
        end
    end
    return poles
end

# The fitted spectral parts are rational in u = ω², so a pole u gives the rate
# z = sqrt(-u) with positive real part. A negative real u is one real mode. A conjugate
# pair u, conj(u) gives exp(-z*t) and exp(-conj(z)*t), represented by one real two-mode
# block; the pole below the real axis is therefore skipped. Positive real u lies on the
# real frequency axis, has no decaying exponential, and is discarded.
function thermal_modes(poles)
    modes = Tuple{Float64,Float64}[]
    for u in poles
        isfinite(u) && imag(u) >= 0 || continue
        z = sqrt(-u)
        real(z) > 0 && push!(modes, (real(z), abs(imag(z))))
    end
    return sort!(modes)
end

function scaled_least_squares(A, b)
    scales = [norm(column) for column in eachcol(A)]
    return ((A ./ transpose(scales)) \ b) ./ scales
end

# Least-squares coefficients of fixed decaying modes. A real mode contributes
# c*exp(-a*t) for t ≥ 0. A block with rates (a, a), weights (1, 0) and mixing
# [0 -b; b 0] contributes exp(-a*t)*(c₁*cos(b*t) + c₂*sin(b*t)), as in
# brownian_oscillator_bath. The spectrum 2real(response*c) has even part
# 2real(response)*real(c) and odd part -2imag(response)*imag(c), fitted separately.
# A constant even part is a white-noise remainder 2D, kept only when D ≥ 0.
function fit_thermal_modes(modes, ω, even, odd)
    n = sum(mode -> iszero(mode[2]) ? 1 : 2, modes; init = 0)
    rates, weights = zeros(n), ones(n)
    mixing = SparseArrays.spzeros(n, n)
    response = zeros(ComplexF64, length(ω), n)
    k = 0
    for (a, b) in modes
        α = a .- im .* ω
        if iszero(b)
            k += 1
            rates[k] = a
            response[:, k] = inv.(α)
        else
            rates[k+1], rates[k+2], weights[k+2] = a, a, 0
            mixing[k+1, k+2], mixing[k+2, k+1] = -b, b
            response[:, k+1] = α ./ (α .^ 2 .+ b^2)
            response[:, k+2] = b ./ (α .^ 2 .+ b^2)
            k += 2
        end
    end
    symmetric, antisymmetric = 2real.(response), -2imag.(response)
    even_fit = scaled_least_squares([symmetric ones(length(ω))], even)
    if even_fit[end] < 0
        even_fit = [scaled_least_squares(symmetric, even); 0]
    end
    odd_fit = scaled_least_squares(antisymmetric, ω .* odd)
    residual = maximum(
        abs.(symmetric * even_fit[1:n] .+ even_fit[end] - even) +
        abs.(antisymmetric * odd_fit - ω .* odd),
    )
    coefficients = complex.(even_fit[1:n], odd_fit)
    return (; coefficients, rates, weights, mixing, diffusion = even_fit[end] / 2, residual)
end

"""
    aaa_bath(spectral_density; kT, frequencies, hbar = 1.0, reltol = 1e-8,
             max_modes = 40, reorganization = nothing)

Thermal [`ExponentialBath`](@ref) for a general spectral density `J(ω)`, obtained by
fitting its noise spectrum with the AAA rational algorithm, as in free-pole HEOM
(Xu, Yan, Shi, Ankerhold and Stockburger, Phys. Rev. Lett. **129**, 230601 (2022),
doi:10.1103/PhysRevLett.129.230601). The fitted rates are not tied to the poles of
`J` or to Matsubara frequencies, so a fit typically needs far fewer modes than a
Matsubara expansion, especially at low temperature.

`spectral_density(ω)` is called only at positive angular frequencies and must return a
finite real `J(ω) ≥ 0`; `J(-ω) = -J(ω)` is implied. The correlation convention is that
of [`drude_lorentz_bath`](@ref),
`C(t) = (ħ/π)∫₀∞ J(ω)[coth(ħω/(2kT))*cos(ωt) - im*sin(ωt)]dω`, whose spectrum is
`S(ω) = 2ħJ(ω)/(1 - exp(-ħω/kT))` as in [`bath_spectrum`](@ref).

`frequencies` are positive sample frequencies; duplicates are ignored. The fit matches
`S(ω)` and `S(-ω)` at each of them. The even part `ħJ(ω)coth(ħω/(2kT))` and the odd
part `ħJ(ω)` divided by `ω` are both even functions, so AAA fits them jointly as
rational functions of `ω²` with shared poles (set-valued AAA; Lietaert et al., IMA J.
Numer. Anal. **42**, 1087 (2022)). Each pole `u` gives a rate `sqrt(-u)` with positive
real part, so complex rates occur in conjugate pairs; poles on the positive real `u` axis
do not decay and are discarded. The coefficients of the retained rates are refitted by
linear least squares on the samples, with a constant even part as nonnegative
white-noise `diffusion`. The degree increases until

    |bath_spectrum(bath, ±ω) - S(±ω)| ≤ reltol * maximum(S)

at every sample, and throws if that needs more than `max_modes` modes. Sample densely
where `S` varies: on the scale `kT/ħ` near zero, across features of `J`, and until `J`
has decayed. Logarithmic spacing, such as `exp10.(range(-3, 2; length = 500))` for
`J` on unit scales, is usually suitable. The fit is not controlled outside the sampled
range. Very low sample frequencies produce correspondingly slow modes, and singular
low-frequency behaviour, as for sub-Ohmic `J`, needs many poles; do not sample below
the inverse of the longest relevant time.

Every mode is real. A real rate is one hierarchy mode, and a conjugate pair is a real
two-mode block with equal `rates`, `weights` `(1, 0)` and antisymmetric `mixing`, as in
[`brownian_oscillator_bath`](@ref). Hence `length(bath.rates)` is the hierarchy mode
count, with two modes for each damped oscillation. The fit can be followed by
[`compress_bath`](@ref).

The counterterm is the reorganization `λ = (1/π)∫₀∞ J(ω)/ω dω`, consistent with the
other constructors. It is computed by exp-sinh quadrature over the frequencies
`sqrt(minimum(frequencies)*maximum(frequencies)) * 10^±30`, unless `reorganization`
supplies it. The quadrature throws if it cannot integrate `J(ω)/ω` there, as for a
strongly sub-Ohmic `J` or a very narrow resonance; then pass the exact `λ`. The
imaginary part of the fitted correlation balances this counterterm only to the accuracy
of the fit of `J`, and misses `J(ω)/ω` below the lowest sample.

A small spectral residual does not establish accurate thermalization (Tokieda, Phys.
Rev. Research **7**, 043178 (2025)). Compare [`harmonic_covariance`](@ref) with an
exact continuum equilibrium where possible, and converge results in `reltol`, the
sampled range and the hierarchy depth.

`kT` and `hbar` must be finite and positive, `reltol` finite and positive, and
`max_modes` nonnegative. A spectral density that vanishes at every sample returns a
bath with no modes and the counterterm `λ`.
"""
function aaa_bath(
    spectral_density;
    kT::Real,
    frequencies::AbstractVector{<:Real},
    hbar::Real = 1.0,
    reltol::Real = 1e-8,
    max_modes::Integer = 40,
    reorganization::Union{Nothing,Real} = nothing,
)
    θ, ħ, tolerance = Float64(kT), Float64(hbar), Float64(reltol)
    isfinite(θ) && θ > 0 || throw(ArgumentError("kT must be finite and positive"))
    isfinite(ħ) && ħ > 0 || throw(ArgumentError("hbar must be finite and positive"))
    isfinite(tolerance) && tolerance > 0 ||
        throw(ArgumentError("reltol must be finite and positive"))
    max_modes >= 0 || throw(ArgumentError("max_modes must be nonnegative"))
    ω = sort!(unique(Float64.(frequencies)))
    !isempty(ω) && all(x -> isfinite(x) && x > 0, ω) ||
        throw(ArgumentError("frequencies must be nonempty, finite and positive"))
    J = [spectral_density_value(spectral_density, x) for x in ω]
    λ = if isnothing(reorganization)
        reorganization_integral(spectral_density, sqrt(first(ω)) * sqrt(last(ω)))
    else
        Float64(reorganization)
    end
    isfinite(λ) && λ >= 0 ||
        throw(ArgumentError("reorganization must be finite and nonnegative"))
    # S(±ω) = even ± ω*odd, with both parts even functions of ω.
    even = @. ħ * J / tanh(ħ * ω / (2θ))
    odd = @. ħ * J / ω
    all(isfinite, even) && all(isfinite, odd) ||
        throw(ArgumentError("the thermal spectrum must be finite at every frequency"))
    scale = maximum(even + ω .* odd)
    iszero(scale) &&
        return ExponentialBath(ComplexF64[], Float64[]; hbar = ħ, counterterm = λ)

    # Set-valued AAA in u = ω². Rows of the stacked Loewner matrix measure both parts in
    # units of S, matching the stopping test. At most n ÷ 2 support points keep the
    # Loewner system overdetermined.
    u, n = ω .^ 2, length(ω)
    support = Int[]
    even_approximation, odd_approximation = fill(sum(even) / n, n), fill(sum(odd) / n, n)
    best, best_modes = Inf, 0
    for _ in 0:min(max_modes, n÷2-1)
        deviation = @. abs(even - even_approximation) + ω * abs(odd - odd_approximation)
        deviation[support] .= 0
        push!(support, argmax(deviation))
        rest = setdiff(1:n, support)
        cauchy = inv.(u[rest] .- transpose(u[support]))
        loewner = [
            (even[rest] .- transpose(even[support])) .* cauchy
            ω[rest] .* (odd[rest] .- transpose(odd[support])) .* cauchy
        ]
        w = svd(loewner).V[:, end]
        denominator = cauchy * w
        even_approximation[rest] = (cauchy * (w .* even[support])) ./ denominator
        odd_approximation[rest] = (cauchy * (w .* odd[support])) ./ denominator
        modes = thermal_modes(barycentric_poles(u[support], w))
        (; coefficients, rates, weights, mixing, diffusion, residual) =
            fit_thermal_modes(modes, ω, even, odd)
        residual <= tolerance * scale && return ExponentialBath(
            coefficients,
            rates;
            hbar = ħ,
            diffusion,
            counterterm = λ,
            weights,
            mixing,
        )
        residual < best && ((best, best_modes) = (residual, length(rates)))
    end
    throw(
        ArgumentError(
            "the AAA fit did not reach reltol = $tolerance within max_modes = " *
            "$max_modes modes and $n samples; the best relative residual was " *
            "$(best / scale) with $best_modes modes. Increase max_modes or reltol, " *
            "or adjust the sampled frequencies",
        ),
    )
end
