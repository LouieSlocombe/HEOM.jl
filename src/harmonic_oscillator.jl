"""
    harmonic_potential(; mass, omega)

Return the potential `q -> mass * omega^2 * q^2 / 2` of a harmonic oscillator with angular
frequency `omega`. `mass` must be finite and positive; `omega` finite and nonnegative.
Setting `omega = 0` gives the free-particle potential.
"""
function harmonic_potential(; mass::Real, omega::Real)
    check_oscillator_parameters(mass, omega; allow_zero_frequency = true)
    return q -> mass * omega^2 * q^2 / 2
end

"""
    coherent_wigner(q, p; q0 = 0.0, p0 = 0.0, mass, omega, hbar = 1.0)

Wigner function of the oscillator's coherent state centred at `(q0, p0)`,

    W(q, p) = exp(-mω(q - q₀)²/ħ - (p - p₀)²/(mωħ)) / (πħ).

It is the ground state displaced in phase space, and the oscillator carries it rigidly along
the classical orbit through `(q0, p0)`. `mass`, `omega` and `hbar` must be finite and positive.
"""
function coherent_wigner(
    q::Real,
    p::Real;
    q0::Real = 0.0,
    p0::Real = 0.0,
    mass::Real,
    omega::Real,
    hbar::Real = 1.0,
)
    Q, P = oscillator_coordinates(q - q0, p - p0, mass, omega, hbar)
    return exp(-Q^2 - P^2) / (π * hbar)
end

"""
    fock_wigner(n, q, p; mass, omega, hbar = 1.0)

Wigner function of the oscillator's `n`-th energy eigenstate,

    Wₙ(q, p) = (-1)ⁿ exp(-2H/ħω) Lₙ(4H/ħω) / (πħ),   H = p²/2m + mω²q²/2,

where `Lₙ` is a Laguerre polynomial. It is stationary, and negative at the origin for odd
`n`. `mass`, `omega` and `hbar` must be finite and positive.
"""
function fock_wigner(
    n::Integer,
    q::Real,
    p::Real;
    mass::Real,
    omega::Real,
    hbar::Real = 1.0,
)
    n >= 0 || throw(ArgumentError("Fock state index must be non-negative, got $n"))
    Q, P = oscillator_coordinates(q, p, mass, omega, hbar)
    r2 = Q^2 + P^2 # 2H/ħω
    return (-1)^n * exp(-r2) * laguerre(n, 2r2) / (π * hbar)
end

"""
    cat_wigner(q, p; q0, mass, omega, hbar = 1.0)

Wigner function of the normalised even cat state `|α⟩ + |-α⟩`, the superposition of the
coherent states at `(±q0, 0)`:

    W = [G(Q - Q₀, P) + G(Q + Q₀, P) + 2G(Q, P) cos(2PQ₀)] / (2πħ(1 + exp(-Q₀²))).

Here `G(x, y) = exp(-x² - y²)`, `Q = q√(mω/ħ)`, `P = p/√(mωħ)` and `Q₀ = q0√(mω/ħ)`. The
interference term makes `W` negative between the two peaks. `mass`, `omega` and `hbar`
must be finite and positive.
"""
function cat_wigner(q::Real, p::Real; q0::Real, mass::Real, omega::Real, hbar::Real = 1.0)
    Q, P = oscillator_coordinates(q, p, mass, omega, hbar)
    Q0 = first(oscillator_coordinates(q0, 0, mass, omega, hbar))
    peaks = exp(-(Q - Q0)^2 - P^2) + exp(-(Q + Q0)^2 - P^2)
    fringes = 2 * exp(-Q^2 - P^2) * cos(2 * P * Q0)
    return (peaks + fringes) / (2π * hbar * (1 + exp(-Q0^2)))
end

"""
    harmonic_evolution(W0, t; mass, omega)

Exact solution at time `t` of the Wigner–Moyal equation for the harmonic oscillator,
starting from the function `W0(q, p)`. The oscillator's Moyal series stops at the classical
term, so `W` is carried along the classical flow:

    W(q, p, t) = W0(q cos ωt - p sin ωt/(mω), p cos ωt + mωq sin ωt).

Returns the function `(q, p) -> W(q, p, t)`. `mass` must be finite and positive, `omega`
finite and nonnegative, and `t` finite. At `omega = 0`, the continuous free-particle limit
is `W0(q - p*t/mass, p)`.
"""
function harmonic_evolution(W0, t::Real; mass::Real, omega::Real)
    check_oscillator_parameters(mass, omega; allow_zero_frequency = true)
    isfinite(t) || throw(ArgumentError("evolution time must be finite"))
    angle = omega * t
    s, c = sincos(angle)
    # sin(angle)/angle has a continuous value of one at zero. Use the same sine as
    # the momentum map so the two maps retain a consistent phase even at large times.
    transport = t * (iszero(angle) ? one(s) : s / angle) / mass
    return (q, p) -> W0(q * c - p * transport, p * c + mass * omega * q * s)
end

# Dimensionless oscillator coordinates Q = q√(mω/ħ) and P = p/√(mωħ), in which the
# oscillator's flow is a rotation and the ground state is exp(-Q² - P²)/(πħ).
function oscillator_coordinates(q, p, mass, omega, hbar)
    check_oscillator_parameters(mass, omega)
    isfinite(hbar) && hbar > 0 || throw(ArgumentError("hbar must be finite and positive"))
    return q * sqrt(mass * omega / hbar), p / sqrt(mass * omega * hbar)
end

function check_oscillator_parameters(mass, omega; allow_zero_frequency = false)
    isfinite(mass) && mass > 0 || throw(ArgumentError("mass must be finite and positive"))
    if allow_zero_frequency
        isfinite(omega) && omega >= 0 ||
            throw(ArgumentError("omega must be finite and nonnegative"))
    else
        isfinite(omega) && omega > 0 ||
            throw(ArgumentError("omega must be finite and positive"))
    end
    return nothing
end

# Laguerre polynomial Lₙ(x) from the three-term recurrence.
function laguerre(n::Integer, x::Real)
    n == 0 && return one(x)
    previous, current = one(x), 1 - x
    for k in 1:(n-1)
        previous, current = current, ((2k + 1 - x) * current - k * previous) / (k + 1)
    end
    return current
end
