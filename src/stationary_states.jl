"""
    eigenstates(grid::PhaseSpaceGrid; mass, potential, hbar = 1.0,
                nstates = length(grid.q))

Compute the lowest `nstates` eigenstates of `H = P²/(2mass) + potential(Q)` on the
position grid. Return `(; energies, wavefunctions)`, with ascending energies and one
wavefunction per column. The columns obey `wavefunctions' * wavefunctions * grid.dq ≈ I`.
Eigenvector signs, and the basis within a degenerate eigenspace, are arbitrary.

The kinetic energy uses Fourier collocation with periodic boundary conditions on
`[first(grid.q), first(grid.q) + length(grid.q)*grid.dq)`. The potential is sampled
only at `grid.q` and need not accept automatic differentiation. `mass` and `hbar`
must be finite and positive; the potential must return finite real values.
`nstates` is between one and the number of position points. Momentum grid limits
and spacing do not affect this eigensolve.

This dense eigensolve uses O(nq²) storage and O(nq³) work. For localized states,
enlarge the position box until the wavefunctions decay at its edges and converge
the position spacing. Continuum states and unconfined potentials describe this
finite periodic box, not normalizable eigenstates on the whole real line.
Use [`wavefunction_wigner`](@ref) to transform columns or coherent superpositions;
[`eigenstate_wigner`](@ref) directly prepares one eigenstate on the phase-space grid.
"""
function eigenstates(
    grid::PhaseSpaceGrid;
    mass::Real,
    potential,
    hbar::Real = 1.0,
    nstates::Integer = length(grid.q),
)
    mass, hbar = Float64(mass), Float64(hbar)
    isfinite(mass) && mass > 0 || throw(ArgumentError("mass must be finite and positive"))
    isfinite(hbar) && hbar > 0 || throw(ArgumentError("hbar must be finite and positive"))
    nq = length(grid.q)
    1 <= nstates <= nq ||
        throw(ArgumentError("nstates must be between 1 and $nq, got $nstates"))

    # F⁻¹ diag(ħ²k²/2m) F is a real symmetric circulant matrix. Keep the even-grid
    # Nyquist kinetic energy: unlike an odd derivative, k² has no sign ambiguity.
    k = 2π .* fftfreq(nq, 1 / grid.dq)
    kinetic = real.(ifft((hbar .* k) .^ 2 ./ (2mass)))
    H = [kinetic[mod(i-j, nq)+1] for i in 1:nq, j in 1:nq]
    for i in 1:nq
        H[i, i] += potential_value(potential, grid.q[i])
    end
    all(isfinite, H) || throw(
        ArgumentError("Hamiltonian entries must be finite; rescale mass, hbar or grid"),
    )
    states = eigen(Symmetric(H), 1:Int(nstates))
    return (; energies = states.values, wavefunctions = states.vectors ./ sqrt(grid.dq))
end

"""
    eigenstate_wigner(n, grid::PhaseSpaceGrid; mass, potential, hbar = 1.0)

Wigner function of the `n`-th numerical energy eigenstate in an arbitrary real
potential. The index is zero-based, as in [`fock_wigner`](@ref): `n = 0` is the
ground state. It must be smaller than `length(grid.q)`.

Calls [`eigenstates`](@ref) followed by [`wavefunction_wigner`](@ref), returning a
real matrix of size `size(grid)` with unit phase-space integral up to roundoff.
The eigensolve uses a finite periodic position box; the transform interpolates
the samples with zero extension. Converge both grid spacings and box sizes and
check that localized states decay at the boundaries.
"""
function eigenstate_wigner(
    n::Integer,
    grid::PhaseSpaceGrid;
    mass::Real,
    potential,
    hbar::Real = 1.0,
)
    0 <= n < length(grid.q) ||
        throw(ArgumentError("eigenstate index must be between 0 and $(length(grid.q) - 1)"))
    states = eigenstates(grid; mass, potential, hbar, nstates = n + 1)
    return wavefunction_wigner(view(states.wavefunctions, :, n + 1), grid; hbar)
end

"""
    thermal_wigner(grid::PhaseSpaceGrid; mass, potential, kT, hbar = 1.0,
                   nstates = length(grid.q))

Prepare the isolated-system Gibbs state `ρ = exp(-H/kT) / Tr(exp(-H/kT))` for an
arbitrary real potential, and return its real Wigner matrix on `grid`. `kT` is
Boltzmann's constant times temperature and must be finite and positive. For a
pure ground state use [`eigenstate_wigner`](@ref) with index zero.

Diagonalizes the finite periodic position-grid Hamiltonian with [`eigenstates`](@ref).
Boltzmann weights are evaluated relative to the lowest energy to avoid underflow
at low temperature and to remove dependence on the potential's energy zero.
By default all position-grid eigenstates are included. Supplying fewer `nstates`
normalizes the Gibbs state within that subspace: increase it until omitted levels
have negligible thermal weight. Grid and finite-box errors require independent
convergence checks, especially at high temperature.

The position density kernel is passed to [`density_matrix_wigner`](@ref), which
preserves its unit quadrature trace. For an unconfined potential (including Morse
above dissociation), this is a finite-box thermal state; an infinite-volume Gibbs
state need not exist. This isolated Gibbs preparation does not include system–bath
correlations. Use [`equilibrate`](@ref) to prepare a correlated HEOM equilibrium.
"""
function thermal_wigner(
    grid::PhaseSpaceGrid;
    mass::Real,
    potential,
    kT::Real,
    hbar::Real = 1.0,
    nstates::Integer = length(grid.q),
)
    kT = Float64(kT)
    isfinite(kT) && kT > 0 || throw(ArgumentError("kT must be finite and positive"))
    states = eigenstates(grid; mass, potential, hbar, nstates)
    weights = exp.(-(states.energies .- first(states.energies)) ./ kT)
    weights ./= sum(weights)
    ψ = states.wavefunctions
    ρ = (ψ .* transpose(weights)) * transpose(ψ)
    return density_matrix_wigner(ρ, grid; hbar)
end
