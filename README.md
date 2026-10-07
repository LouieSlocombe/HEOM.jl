# HEOM.jl

Phase-space quantum dynamics in Julia, including hierarchical equations of motion
(HEOM) for a particle coupled linearly to a Gaussian bath. The package solves the
Wigner–Moyal equation, its bath hierarchy, and the Caldeira–Leggett model with
friction and thermal diffusion. The Wigner–Moyal equation is the
exact quantum evolution of the Wigner function $W(q, p, t)$ of a particle of mass
$m$ in a static or time-dependent potential $V(q,t)$:

```math
\frac{\partial W}{\partial t} = -\frac{p}{m}\frac{\partial W}{\partial q}
+ \sum_{s=0}^{\infty} \frac{(-1)^s}{(2s+1)!}\left(\frac{\hbar}{2}\right)^{2s}
\frac{\partial^{2s+1} V}{\partial q^{2s+1}}\frac{\partial^{2s+1} W}{\partial p^{2s+1}}
```

It uses the method of lines. Phase space is discretised on a periodic grid, and the
result is a SciML `ODEProblem`, so any solver from the DifferentialEquations.jl
ecosystem integrates it in time.

## Requirements

- Julia 1.10 or newer in the 1.x series. Install it with
  [Juliaup](https://julialang.org/downloads/).
- An ODE solver package in your environment, for example `OrdinaryDiffEqVerner`,
  `OrdinaryDiffEqTsit5`, or the full `OrdinaryDiffEq`.
- Optional: `Plots` for Wigner heatmaps, marginal densities, diagnostic plots and animations.
- Optional: [pre-commit](https://pre-commit.com/#installation) for Git hooks.

CI tests the minimum supported Julia version and the latest stable Julia on
Linux, macOS, and Windows.

## Quick start

From the repository root, instantiate and precompile the package:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

A coherent state in a harmonic well, solved over one period:

```julia
using HEOM, OrdinaryDiffEqVerner

mass, omega = 1.0, 1.0
grid = PhaseSpaceGrid((-8.0, 8.0), 64, (-8.0, 8.0), 64)
W0 = on_grid((q, p) -> coherent_wigner(q, p; q0 = 2.0, p0 = 1.0, mass, omega), grid)
V = harmonic_potential(; mass, omega)

prob = wigner_moyal_problem(W0, (0.0, 2π / omega), grid; mass, potential = V)
sol = solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12)

W = sol.u[end]                           # back to W0 after one period
phase_space_integral(W, grid)            # norm, 1
expectation((q, p) -> q, W, grid)        # ⟨q⟩, 2
purity(W, grid)                          # Tr ρ², 1 for a pure state
d = diagnostics(sol; potential = V)      # observables and grid health at every saved time
```

A Wigner function on the grid is an `nq × np` matrix with `W[i, j] = W(q[i], p[j])`.
The solution's states are `sol.u` at times `sol.t`.

To work on this package from another Julia project, use
`Pkg.develop(path="/path/to/HEOM.jl")` in that project's environment. Every public
function has a docstring available through Julia's help mode, for example
`?wigner_moyal_operator`.

## Driven dynamics and spectroscopy

Wrap a general potential as `TimeDependentPotential((q, t) -> V(q, t))`, or use
`DrivenPotential(V0, field, dipole)` for `V(q,t) = V0(q) - field(t)*dipole(q)`.
Both work with the Wigner–Moyal, HEOM and Caldeira–Leggett problem builders and
with spectral or finite-difference discretisation. A plain callable `q -> V(q)`
continues to describe a static potential; plain two-argument callables are also
accepted for driven dynamics. For separable drives the spatial
operators are cached and only the field amplitude changes during propagation.

```julia
using HEOM, OrdinaryDiffEqVerner

mass, omega, hbar = 1.0, 1.0, 1.0
grid = PhaseSpaceGrid((-7, 7), 64, (-7, 7), 64)
W0 = on_grid((q, p) -> coherent_wigner(q, p; mass, omega, hbar), grid)
V0 = harmonic_potential(; mass, omega)
field(t) = 0.3 * cos(0.75t)
potential = DrivenPotential(V0, field, q -> q)
# Equivalent general form:
# potential = TimeDependentPotential((q, t) -> V0(q) - field(t)*q)
prob = wigner_moyal_problem(W0, (0.0, 8.0), grid; mass, hbar, potential)
sol = solve(prob, Vern7(); saveat = 0.05, abstol = 1e-10, reltol = 1e-10)
d = diagnostics(sol; potential)          # instantaneous driven-system energy
```

The drive acts on every HEOM auxiliary, while the bath and its counterterm remain
fixed. The generator is evaluated at the solver's current time. HEOM also supplies
time-aware Jacobian products and partial time derivatives for stiff methods.
Time derivatives are obtained with ForwardDiff unless supplied explicitly as
`TimeDependentPotential(V; derivative = (q,t) -> dVdt(q,t))` or
`DrivenPotential(V0, field, dipole; field_derivative = t -> dEdt(t))`.
For discontinuous or sharply varying fields, provide solver `tstops` at the
switches and resolve the pulse with appropriate integration tolerances and steps.
At one saved time use `energy(W, grid; mass, potential, t)` or
`diagnostics(W, grid; mass, hbar, potential, t)`; rates likewise accept `t`.
Finite-difference operators provide the instantaneous matrix as `sparse(op, t)`;
`sparse(op)` requires a static generator. HEOM's `jacobian = :sparse` also updates
its matrix at the current solver time for driven problems.

Linear response starts from an equilibrium state and an **undriven** operator:

```julia
op = wigner_moyal_operator(grid; mass, hbar, potential = V0)
response = linear_response(W0, (0.0, 30.0), op, Vern7();
    dipole = q -> q, saveat = 0.025, abstol = 1e-10, reltol = 1e-10)
spectrum = absorption_spectrum(response;
    frequencies = 0.0:0.01:2.0, broadening = 0.3)
response.times                         # elapsed times since tspan[1]
response.response                      # causal impulse response R(t)
spectrum.susceptibility                 # complex retarded susceptibility
spectrum.intensity                     # omega * imag(susceptibility)
```

For the perturbation `-E(t)*mu(q)`, the helper propagates
`delta_rho(0) = (im/hbar)*[mu, rho_eq]` and evaluates
`R(t) = Tr[observable*delta_rho(t)]`. `observable` defaults to `dipole`; both are
functions of position. The physical root determines the signal. The response
state is signed and has zero trace; it must not be normalised as a density.
`linear_response_problem(initial, tspan, op; dipole)` exposes the same seeded
problem for custom solve workflows. The helper accepts isolated Wigner matrices
or full HEOM hierarchies. For a coupled bath, prepare a converged correlated
equilibrium with `equilibrate`, then pass its full hierarchy so the commutator
also acts on every correlated auxiliary, or call
`linear_response(eq, tspan, alg; ...)` directly on the converged result.
Applying the excitation to the full hierarchy preserves system–bath correlations;
see [Tanimura's HEOM review, Appendix A](https://arxiv.org/html/2006.05501#A1).
A factorized initial state describes the response of that preparation and need
not yield an equilibrium spectrum.

The spectrum convention is
`chi(omega) = integral exp(im*omega*t - broadening*t)*R(t) dt`, evaluated by
trapezoidal quadrature over the saved delays. Frequencies are angular frequencies;
`broadening >= 0` has units of inverse time. Returned intensity is proportional
to absorption; experimental electromagnetic prefactors are not included.
`absorption_spectrum(times, response; frequencies, broadening)` also accepts
sampled response data. Increase the response duration, refine time sampling and
converge the phase-space grid independently. Exponential broadening controls
finite-window ringing and does not represent a physical bath by itself.

Start with [the driven harmonic benchmark](examples/driven_harmonic.jl), which
checks the forced centroid, `R(t) = sin(omega*t)/(mass*omega)` and the analytic
broadened susceptibility. Then run
[the anharmonic spectroscopy example](examples/anharmonic_spectroscopy.jl), which
compares a quartic oscillator's response against an eigenstate Kubo sum, resolves
its shifted absorption peak and compares a weak pulse with response convolution.
Both run in an environment containing HEOM and OrdinaryDiffEqVerner:

```sh
julia --project=/path/to/environment examples/driven_harmonic.jl
julia --project=/path/to/environment examples/anharmonic_spectroscopy.jl
```

## General initial states

Prepare a Wigner function from a wavefunction or a position-space density kernel:

```julia
grid = PhaseSpaceGrid((-8.0, 8.0), 192, (-8.0, 8.0), 128)
hbar = 1.0
psi(q) = π^(-1/4) * exp(-q^2 / 2 + im * q / hbar)
W = wavefunction_wigner(psi, grid; hbar)
W_samples = wavefunction_wigner(psi.(grid.q), grid; hbar)

rho(q, qp) = psi(q) * conj(psi(qp))
W_rho = density_matrix_wigner(rho, grid; hbar)
samples = psi.(grid.q)
W_matrix = density_matrix_wigner(samples * samples', grid; hbar)
```

The convention is

```math
W(q,p) = \frac{1}{2\pi\hbar}\int e^{ipy/\hbar}
\rho(q-y/2,q+y/2)\,dy,\qquad
\rho(q,q')=\psi(q)\psi(q')^*.
```

A sampled density matrix contains kernel values `rho[i, j] = rho(q[i], q[j])`;
its trace is `sum(diag(rho))*grid.dq`. Wavefunction samples have squared norm
`sum(abs2, psi)*grid.dq`. The transforms preserve these quantities without
normalising the input, and return real `nq × np` arrays that retain Wigner negativity.
Sampled inputs use local cubic interpolation (lower order on very small grids)
and are zero outside the sampled position interval. Resolve their spatial
structure and make their boundary values negligible. Callable inputs are evaluated
at shifted positions up to approximately `π*hbar/(2*grid.dp)` beyond each `q`.
Both forms need convergence checks in the position and momentum boxes and spacings.

Numerical eigenstates and Gibbs states accept any finite real potential on `grid.q`:

```julia
mass = 1.0
V(q) = (q^2 - 4)^2 / 8                 # quartic double well
states = eigenstates(grid; mass, potential = V, hbar, nstates = 2)
states.energies                        # ascending energies
psi0 = states.wavefunctions[:, 1]      # quadrature-normalised columns
psi1 = states.wavefunctions[:, 2]
W_superposition = wavefunction_wigner((psi0 + psi1) / sqrt(2), grid; hbar)
W_ground = eigenstate_wigner(0, grid; mass, potential = V, hbar)
W_thermal = thermal_wigner(grid; mass, potential = V, kT = 0.3, hbar)
```

`eigenstate_wigner` uses zero-based state indices, as `fock_wigner` does.
`eigenstates` diagonalises a dense Fourier spectral Hamiltonian with periodic
position boundaries; it is intended for moderate position grids and localised
states that decay before those boundaries. Eigenvector phases are arbitrary,
so set relative phases explicitly when choosing a particular superposition.
`thermal_wigner` prepares the normalised, isolated-system finite-box Gibbs state
for positive `kT = k_B T`. It includes the full discrete spectrum by default;
`nstates` truncates the Gibbs sum and must be increased to check convergence.
For an unconfined potential such as Morse, the continuum is discretised by the box:
its finite-box thermal state does not establish a normalisable Gibbs state on the
infinite line. For correlated system–bath equilibrium, use [`equilibrate`](#correlated-equilibrium-preparation).

The [initial-state example](examples/initial_states.jl) prepares a tunnelling
superposition in a double well, compares wavefunction and density-matrix transforms,
and prepares a Morse ground state. It needs only HEOM:

```sh
julia --project=. examples/initial_states.jl
```

## Wigner-space HEOM

HEOM evolves the physical Wigner function together with auxiliary density
operators (ADOs), represented on the same phase-space grid. Independent exponential
bath correlations are `C(t) = sum(c[k] * exp(-rates[k] * t))` for `t ≥ 0`.
The supported coupling is linear in the one-dimensional coordinate. The generalized
real basis also supports damped oscillations and coincident poles through weights
and mixing: `C(t) = transpose(weights)*exp(-Γ*t)*coefficients`, with
`Γ = Diagonal(rates) + mixing`. For independent exponentials with unit weights,
the unscaled auxiliary Wigner functions obey

```math
\partial_t W_{\mathbf n} =
\left(\mathcal L_{\mathrm{WM}}[V+\Lambda q^2]
-\sum_k n_k\nu_k+D\partial_p^2\right)W_{\mathbf n}
+\sum_k\partial_p W_{\mathbf n+\mathbf e_k}
+\sum_k n_k\left(\operatorname{Re}c_k\partial_p
+\frac{2\operatorname{Im}c_k}{\hbar}q\right)W_{\mathbf n-\mathbf e_k}.
```

Here `W₀` is physical and all ADOs are real arrays; complex bath coefficients
enter through their real and imaginary parts. The hierarchy retains every
multi-index with `sum(n) ≤ depth` and sets upward neighbours beyond that depth
to zero. This is a hard cutoff at the final tier.

The Drude–Lorentz constructor uses the spectral-density convention
`J(ω) = 2λγω / (ω² + γ²)` and
`C(t) = (hbar/π) ∫₀∞ J(ω)[coth(hbar*ω/(2kT))*cos(ω*t) - im*sin(ω*t)] dω`.
Its pole has rate `γ` and coefficient
`c₀ = λ*hbar*γ*(cot(hbar*γ/(2kT)) - im)`. The `matsubara = K` additional poles
have `νₖ = 2π*k*kT/hbar` and `cₖ = 4λγ*kT*νₖ / (νₖ² - γ²)`.

```julia
bath = drude_lorentz_bath(;
    reorganization = 0.12, cutoff = 1.2, kT = 0.8, matsubara = 1, hbar = 1.0)
op = heom_operator(grid; mass, potential = V, bath, depth = 5)
prob = heom_problem(W0, (0.0, 6.0), op)
sol = solve(prob, Vern7(); saveat = 0.1, abstol = 1e-9, reltol = 1e-9)

W = physical_wigner(sol)              # view of the final physical Wigner function
W_initial = physical_wigner(sol, 1)   # first saved physical state
indices = hierarchy_indices(op)      # ADO multi-indices, root first
d = diagnostics(sol; potential = V)   # diagnostics of the physical states
```

`sol.u[i]` is an `nq × np × nados` array, with the physical state in its first
slice. Passing a matrix to `heom_problem` fills all higher ADOs with zero: this
is a factorized system/bare-bath initial condition. A full three-dimensional
initial state can instead supply nonzero ADOs or restart a saved trajectory.
`physical_wigner(U)` returns a view, so copy it if it must be modified independently.

For the Drude bath, `reorganization = λ` sets the counterterm `Λ = λ`:
the operator adds `λ*q²` internally to the supplied physical potential `V`.
For dimensional position coupling, `λ` has units of energy divided by position
squared. Do not add this term to `V` yourself. Zero initial ADOs produce an initial
slip as the bath adjusts to the system. The default `terminator = true` approximates
omitted Matsubara poles by momentum diffusion,
`D = 2λ*kT/γ - sum(real(c[k])/rates[k])` for independent exponentials;
in a mixed basis the sum becomes `transpose(weights)*(Γ\real(coefficients))`,
where `Γ = Diagonal(rates) + mixing`. `terminator = false` sets this diffusion
to zero. This correction is separate from the final-tier hard cutoff.

For another exponential correlation expansion, use
`ExponentialBath(coefficients, rates; hbar = 1, diffusion = 0, counterterm = 0)`.
Here `counterterm` is the coefficient of `q²`, and `diffusion` is the coefficient
of `∂p²`. Constructing a bath does not establish that an arbitrary correlation
expansion describes a physical thermal environment.

Converge results independently in `depth`, `matsubara`, phase-space box and grid
resolution. The Drude constructor requires positive finite temperature and
handles coincident Drude and retained Matsubara poles with a finite mixed basis.
With the terminator enabled, include
enough poles that the first omitted Matsubara rate exceeds `cutoff` and is fast
compared with the system frequencies of interest.
Low temperatures generally require more Matsubara poles; their fast decay rates
can restrict explicit time steps.
The number of ADOs is `binomial(depth + K + 1, K + 1)` for a nonzero Drude bath with
`K` Matsubara poles. Each ADO uses the same spectral or finite-difference Wigner–Moyal
operator, including quantum terms for anharmonic potentials.

See [the HEOM harmonic oscillator example](examples/heom_sho.jl), which checks
the centroid against an exact bath-memory equation, compares two hierarchy depths,
and prints the norm drift. Run it with `julia --project=/path/to/environment
examples/heom_sho.jl` in an environment containing HEOM and OrdinaryDiffEqVerner.
The physical model is discussed in
[Tanimura's Wigner-space HEOM derivation](https://arxiv.org/pdf/1502.04077).
The equation above uses coordinate-coupled, unscaled ADOs; their definitions and
initial conditions differ from the transformed momentum form in that paper.

### Structured baths and correlation diagnostics

`brownian_oscillator_bath` adds a thermal vibrational resonance with spectral density
`J(ω) = 2λ*γ*ω₀²*ω / ((ω₀² - ω²)² + γ²*ω²)`, using the same correlation convention
and counterterm `λ*q²` as the Drude constructor. Set `frequency = ω₀` and
`damping = γ`, with `0 < γ < 2ω₀`; the oscillation decays at `γ/2` and its damped
angular frequency is `sqrt(ω₀² - γ²/4)`. Temperature and `hbar` must be positive.
The bath uses two coupled real modes plus `matsubara` thermal modes. A finite
mixed basis also handles near-coincident oscillator and retained thermal poles.

```julia
vibration = brownian_oscillator_bath(;
    reorganization = 0.08, frequency = 2.0, damping = 0.3,
    kT = 0.8, matsubara = 3)
background = drude_lorentz_bath(;
    reorganization = 0.12, cutoff = 1.2, kT = 0.8, matsubara = 1)
bath = combine_baths(background, vibration)

times = range(0, 10; length = 201)
frequencies = range(-5, 5; length = 301)
correlation = bath_correlation.(Ref(bath), times)
spectrum = bath_spectrum.(Ref(bath), frequencies)
members = hierarchy_size(length(bath.rates), 4)
# Pass bath to heom_operator or heom_problem as usual.
```

`combine_baths(baths...)` or `combine_baths([bath1, bath2, ...])` sums independent
components coupled to the same position. It preserves their bases, adds diffusion
and counterterms, copies all arrays, and requires matching `hbar`. Components at
different temperatures can be combined, but the sum generally has no single
thermal temperature. Correlated components and different coupling operators are
outside this constructor's scope.

Brownian Matsubara coefficients are negative. The constructor retains their signs
and drops the omitted tail, with zero residual diffusion; that tail cannot be
represented by a positive diffusion terminator. Increase `matsubara` to converge
the correlation and spectrum, particularly at low temperature, as well as `depth`
for dynamics. The decomposition follows the
[underdamped Brownian correlation](https://doi.org/10.1038/s41467-019-11656-1),
with coupling expressed as the reorganization coefficient `λ`.

`bath_correlation(bath, t)` evaluates the retained regular correlation, including
weights and mixing; negative times use `C(-t) = conj(C(t))`. It excludes the
white-noise term `2bath.diffusion*δ(t)`, even at zero time.
`bath_spectrum(bath, ω)` evaluates the unsymmetrized, two-sided Fourier transform
`S(ω) = ∫ exp(im*ω*t)*C(t) dt`, including the constant `2bath.diffusion` contribution.
It uses angular frequencies and does not clip negative values. For an exact thermal
bath, `S(ω) = 2hbar*J(ω)/(1-exp(-hbar*ω/kT))` for positive `ω`, and
`S(-ω) = exp(-hbar*ω/kT)*S(ω)`. These diagnostics expose errors from a finite
decomposition instead of replacing it with the ideal thermal spectrum.

### Correlated equilibrium preparation

`equilibrate` relaxes a matrix or full hierarchy under an existing HEOM operator
and tests stationarity of **every** auxiliary, including the physical root:

```julia
eq = equilibrate(W0, (0.0, 100.0), op, Vern7();
    stationarity_abstol = 1e-8, stationarity_reltol = 1e-6,
    check_interval = 1.0, abstol = 1e-10, reltol = 1e-10)

eq.converged                         # check before treating the state as equilibrium
eq.status                            # :stationary, :time_limit, :solver_failure, :terminated
eq.residuals.absolute                # maximum |dU/dt| for each hierarchy member
maximum(eq.residuals.scaled)          # ≤ 1 for full-hierarchy stationarity
W_eq = physical_wigner(eq)            # view of the final root
U_eq = eq.hierarchy                   # includes all correlated auxiliaries

restart = heom_problem(eq, (0.0, 10.0))
sol = solve(restart, Vern7(); abstol = 1e-10, reltol = 1e-10)
```

For each member `a`, the stopping test is
`maximum(abs, dU[:,:,a]) ≤ stationarity_abstol + stationarity_reltol*maximum(abs, U[:,:,a])`.
These RHS tolerances are separate from integration tolerances and apply in the
operator's auxiliary scaling convention. Checks run initially, after accepted steps
at least `check_interval` apart, and at the final time. The time span is a finite
preparation budget; reaching its end does not imply convergence. Failed solves also
return their final hierarchy, residuals and solver `retcode` with `converged = false`.
The input must have positive finite root integral and is copied without normalisation.

`eq.restart` records the final time, requested time span, root norm and initial norm,
hierarchy indices, depth, scaling, Jacobian mode, and the operator (including its
bath and grid). The operator is retained by reference and shares its work buffers
with restarts. For a longer preparation, pass `eq.hierarchy` to `equilibrate` again
with that same operator. Saving only `W_eq` would discard the equilibrium correlations.

Stationarity is a property of the chosen finite hierarchy and grid. Check convergence
in depth, bath expansion, box and resolution separately, and inspect boundary and
spectral diagnostics of the physical state. A coupled harmonic oscillator generally
relaxes to a reduced equilibrium different from the isolated oscillator's Gibbs state.

### Low temperature and strong coupling

Use `drude_lorentz_pade_bath(...; pade=N)` for a compact `[N/N]` Bose Padé expansion
and `heom_operator(...; scaled=true)` for factorial/amplitude-scaled auxiliaries.
Both negative bath residues and repeated poles are supported. This path retains the
full Wigner–Moyal potential and works for anharmonic systems. Increase `N` and
`depth` independently; scaling changes conditioning, not the retained physics.

```julia
bath = drude_lorentz_pade_bath(;
    reorganization = 0.8, cutoff = 0.5, kT = 0.1, hbar = 1, pade = 4)
members = hierarchy_size(length(bath.rates), 6)
op = heom_operator(grid; mass, potential = V, bath, depth = 6,
                   scaled = true, max_ados = 10_000)
prob = heom_problem(W0, (0.0, 0.8), op)
```

`heom_problem` supplies an exact Jacobian-vector product for stiff Krylov solvers.
With `OrdinaryDiffEqRosenbrock`, `LinearSolve` and `ADTypes` installed:

```julia
using OrdinaryDiffEqRosenbrock, LinearSolve, ADTypes
alg = Rodas5P(autodiff = AutoFiniteDiff(),
             linsolve = KrylovJL_GMRES(), concrete_jac = false)
sol = solve(prob, alg; abstol = 1e-9, reltol = 1e-9)
```

The supplied product is exact; no finite differences of the FFT are used.
Finite-difference operators also support `sparse(op)` and
`heom_problem(W0, tspan, op; jacobian=:sparse)` with `Rodas5P()`.
Use `rescale_hierarchy(U, old_op; scaled=true)` when converting an unscaled
correlated state or restart. Passing a matrix still sets higher auxiliaries to zero.

The [cold, strongly coupled anharmonic example](examples/low_temperature_strong_coupling.jl)
separately varies depth, Padé count, grid spacing and box size. The
[support and validation notes](references/low_temperature_strong_coupling.md)
give the equations, exact quantum Brownian-motion comparisons and convergence limits.

### Compressing bath decompositions

The hierarchy size is `binomial(modes + depth, depth)`, so the bath mode count
dominates the cost of a converged low-temperature calculation. At depth 8, eight
modes need 12870 auxiliaries and four need 495. `compress_bath` reduces any
`ExponentialBath` by balanced model-order reduction of `C(t) = wᵀexp(-Γt)c`. For a
Padé bath this is the truncated Padé decomposition of
[Takahashi and Tanimura, J. Chem. Phys. 158, 044115 (2023), Appendix B](https://doi.org/10.1063/5.0135725):

```julia
source = drude_lorentz_pade_bath(;
    reorganization = 0.8, cutoff = 0.5, kT = 0.1, hbar = 1, pade = 7)  # 8 modes
σ = hankel_singular_values(source)        # decide how many modes to keep
bath = compress_bath(source; modes = 4)   # 4 real hierarchy modes
bound = 4sqrt(2) * sum(σ[5:end])          # max over ω of the spectral change
covariance = harmonic_covariance(bath; mass = 1, omega = 1)
```

Every returned mode is real, so `modes` is the hierarchy mode count; a damped
oscillatory pair uses two. The reduced rate matrix is in real Schur form, with
upper-triangular `mixing`. The default `method = :truncate` keeps the diffusion and
counterterm unchanged. `method = :residualize` instead preserves `∫₀^∞ C(t) dt`. It
moves the eliminated real part into `diffusion` and the eliminated static potential
into `counterterm`, and refuses to create negative diffusion, as from Brownian
Matsubara terms. The spectral bound measures the change from `source`, not the error
of `source` itself.

A small spectral or correlation residual does not establish accurate thermalization
([Tokieda, Phys. Rev. Research 7, 043178 (2025)](https://doi.org/10.1103/bv19-dtb1)).
Check a compressed bath three ways:

1. Compare `bath_spectrum` with the thermal target over the frequencies that the
   system resolves.
2. Compare `harmonic_covariance`, the exact equilibrium of a harmonic oscillator
   coupled to the decomposition, for the source, the compressed bath and a continuum
   result.
3. Repeat the hierarchy-depth convergence with the compressed bath.

The [bath compression example](examples/bath_compression.jl) performs all three
checks on the cold anharmonic oscillator above. Four balanced modes from the
eight-mode source change the final Wigner function by about `2e-6`, while four Padé
modes differ from the source by `1.4e-3`.
For spectra without a thermal constructor, a rational fit such as AAA
([Xu et al., Phys. Rev. Lett. 129, 230601 (2022)](https://doi.org/10.1103/PhysRevLett.129.230601))
can supply the source; count the real modes of any complex poles it returns.

## Caldeira–Leggett damping

The high-temperature, Markovian Caldeira–Leggett model couples the particle to a
thermal bath:

```math
\frac{\partial W}{\partial t} = \mathcal{L}_{\mathrm{WM}} W
+ \gamma\frac{\partial(pW)}{\partial p}
+ m\gamma k_B T\frac{\partial^2 W}{\partial p^2}.
```

Here `friction = γ` is the full momentum damping rate, so the mean momentum obeys
`d⟨p⟩/dt = −⟨V′(q)⟩ − γ⟨p⟩`. `kT` is the thermal energy `k_B T`, in the same units
as the potential. Both must be finite and nonnegative. The quantum Hamiltonian
term uses the same `hbar`, `discretization` and `moyal_terms` options as
`wigner_moyal_problem`. This convention matches the thermal diffusion coefficient
`Dpp = mass * friction * kT` in
[García-Palacios and Zueco, Eq. (4)](https://arxiv.org/pdf/cond-mat/0407454), with
the mixed diffusion coefficient set to zero.

```julia
prob = caldeira_leggett_problem(W0, (0.0, 25.0), grid;
    mass, potential = V, friction = 0.8, kT = 2.0)
sol = solve(prob, Vern9(); abstol = 1e-10, reltol = 1e-10, saveat = 0.25)
d = diagnostics(sol; potential = V)
```

For a harmonic well, the centroid decays to its minimum. At finite temperature,
the state retains a thermal width: its equilibrium variances are
`⟨q²⟩ = kT / (mass * omega²)` and `⟨p²⟩ = mass * kT`, with mean energy `kT`.
This model is a high-temperature approximation (`kT ≫ hbar * omega` for the
oscillator); it does not describe cooling into the quantum ground state.

A Gaussian initial state remains Gaussian. Its exact mean and covariance are

```math
\mu(t)=e^{At}\mu_0,\qquad
\Sigma(t)=\Sigma_\infty+e^{At}(\Sigma_0-\Sigma_\infty)e^{A^\mathsf{T}t},
\quad
A=\begin{pmatrix}0&1/m\\-m\omega^2&-\gamma\end{pmatrix},\quad
\Sigma_\infty=\begin{pmatrix}k_BT/(m\omega^2)&0\\0&mk_BT\end{pmatrix}.
```

See [the damped harmonic oscillator example](examples/damped_sho.jl) for a
comparison of the full Wigner function and its moments against this solution.
The existing observables, diagnostics, plots and rate functions also accept
Caldeira–Leggett operators and solutions.

## API

| Function | Purpose |
|---|---|
| `PhaseSpaceGrid(qlims, nq, plims, np)` | Uniform periodic grid on `[qmin, qmax) × [pmin, pmax)` |
| `on_grid(f, grid)` | Sample `f(q, p)` on the grid |
| `wigner_moyal_problem(W0, tspan, grid; mass, potential, ...)` | `ODEProblem` for the Wigner–Moyal equation |
| `wigner_moyal_operator(grid; mass, potential, hbar, discretization, moyal_terms)` | Reusable semi-discrete operator; `wigner_moyal_problem(W0, tspan, op)` takes it |
| `wigner_moyal!(dW, W, op, t)` | In-place right-hand side |
| `TimeDependentPotential(V; derivative)`, `DrivenPotential(V0, field, dipole; field_derivative)` | General `V(q,t)` or separable `V0(q) - field(t)*dipole(q)` |
| `linear_response(initial, tspan, op, alg; dipole, observable, ...)` | Propagate the impulse response with an undriven generator |
| `linear_response_problem(initial, tspan, op; dipole)` | Response problem seeded by the dipole commutator |
| `absorption_spectrum(response; frequencies, broadening)` | Retarded susceptibility and absorption intensity from saved response data |
| `ExponentialBath(coefficients, rates; hbar, diffusion, counterterm)` | Gaussian bath described by an exponential correlation expansion |
| `drude_lorentz_pade_bath(; reorganization, cutoff, kT, pade, ...)` | Compact quantum Drude bath for low temperatures |
| `hierarchy_size(modes, depth)`, `rescale_hierarchy(U, op; scaled)` | Estimate hierarchy size and convert auxiliary scaling |
| `drude_lorentz_bath(; reorganization, cutoff, kT, matsubara, hbar, terminator)` | Drude–Lorentz bath with Matsubara poles and optional residual diffusion |
| `brownian_oscillator_bath(; reorganization, frequency, damping, kT, matsubara, hbar)` | Underdamped Brownian resonance with retained thermal poles |
| `combine_baths(baths...)`, `combine_baths(baths)` | Sum independent bath components coupled to the same position |
| `bath_correlation(bath, t)`, `bath_spectrum(bath, omega)` | Retained force correlation and unsymmetrized noise spectrum |
| `compress_bath(bath; modes, method)`, `hankel_singular_values(bath)` | Balanced reduction of an exponential bath to fewer real hierarchy modes |
| `harmonic_covariance(bath; mass, omega)` | Exact harmonic-oscillator equilibrium covariance for a bath decomposition |
| `heom_problem(W0, tspan, grid; mass, potential, bath, depth, ...)` | `ODEProblem` for the Wigner-space hierarchy |
| `heom_operator(grid; mass, potential, bath, depth, ...)` | Reusable hierarchy operator; `heom_problem(W0, tspan, op)` takes it |
| `heom!(dU, U, op, t)` | In-place hierarchy right-hand side |
| `hierarchy_indices(op)` | Multi-indices corresponding to the third state-array dimension |
| `physical_wigner(U)`, `physical_wigner(sol[, index])` | View of the physical Wigner function in a hierarchy state or solution |
| `caldeira_leggett_problem(W0, tspan, grid; mass, potential, friction, kT, ...)` | `ODEProblem` with Caldeira–Leggett friction and thermal diffusion |
| `caldeira_leggett_operator(grid; mass, potential, friction, kT, ...)` | Reusable operator; `caldeira_leggett_problem(W0, tspan, op)` takes it |
| `caldeira_leggett!(dW, W, op, t)` | In-place dissipative right-hand side |
| `phase_space_integral`, `expectation`, `energy` | Norm, Weyl-symbol averages and mean energy |
| `position_density`, `momentum_density` | Marginal densities, integrating over the other axis |
| `phase_space_mean`, `phase_space_covariance` | Mean position and momentum, and symmetrised covariance |
| `purity`, `overlap`, `wigner_negativity` | `Tr ρ²`, `Tr ρ₁ρ₂` and the integral of `abs(W) − W` |
| `boundary_weight`, `spectral_tail` | Fractions of absolute weight near each boundary and Fourier amplitude in high modes |
| `diagnostics(W, grid; mass, potential, hbar)`, `diagnostics(states, grid; ...)`, `diagnostics(sol; potential)` | State or trajectory summaries, with autocorrelation for trajectories |
| `wignerplot(W, grid)`, `wignerplot(sol[, index])` | Signed Wigner heatmap, with position horizontal and momentum vertical |
| `marginalplot(W, grid)`, `marginalplot(sol[, index])` | Position and momentum densities in separate panels |
| `diagnosticsplot(d)`, `diagnosticsplot(times, d)`, `diagnosticsplot(sol; potential)` | Diagnostic time series, selectable with `fields` |
| `wigneranimation(sol)`, `wigneranimation(states, grid; times)` | Wigner animation with a fixed signed colour scale |
| `marginalanimation(sol)`, `marginalanimation(states, grid; times)` | Marginal density animation with fixed scales in both panels |
| `probability(W, grid; q, p)`, `probability_current(W, grid; mass)` | Window populations and position probability current |
| `expectation_rate(f, W, op)`, `probability_rate(W, op; q, p)` | Instantaneous rates from the semi-discrete equation |
| `tunnelling_rates(sol; equilibrium_product, tspan)` | Forward and backward transfer constants from population relaxation |
| `wavefunction_wigner(psi, grid; hbar)`, `density_matrix_wigner(rho, grid; hbar)` | Wigner transforms of callable or sampled position-space states |
| `eigenstates(grid; mass, potential, hbar, nstates)` | Numerical energies and quadrature-normalised wavefunction columns |
| `eigenstate_wigner(n, grid; mass, potential, hbar)` | Numerical eigenstate Wigner function, with zero-based `n` |
| `thermal_wigner(grid; mass, potential, kT, hbar, nstates)` | Normalised isolated-system Gibbs state on the finite position box |
| `harmonic_potential`, `coherent_wigner`, `fock_wigner`, `cat_wigner` | Harmonic oscillator and its standard states |
| `harmonic_evolution(W0, t; mass, omega)` | Exact harmonic-oscillator solution, for testing |

## Analysis

Nothing is renormalised: norm drift remains visible in every observable. The
marginals, means, central covariance and energy use the same grid quadrature as
`expectation`. `overlap(W1, W2, grid; hbar)` gives `Tr ρ₁ρ₂`, which is the squared
state overlap for pure states and an eigenstate population when one state is an
energy eigenstate. `wigner_negativity` integrates `abs(W) − W`; it vanishes for a
nonnegative Wigner function. Its quadrature is only second-order accurate at the
zero contours, so resolve these before interpreting small changes.

`diagnostics` collects norm, means, variances, covariance, energy, purity and
negativity alongside two uncertainty measures: `uncertainty` is σqσp, and
`robertson_schrodinger` is √det Σ. Both are at least ħ/2 for a normalised physical
state, with √det Σ invariant under harmonic evolution. For a vector of states or
an ODE solution, each field is a vector and `autocorrelation` records overlap with
the first saved state. For a pure initial state this is the survival probability.
The solution method adds `t` and reads the grid, mass and ħ from the operator.

Keep `boundary_q`, `boundary_p`, `tail_q` and `tail_p` small. Boundary weights are
fractions of the integral of `abs(W)` in the outer strips (5% of the points at each
end by default); enlarge the box when they grow. Spectral tails measure the
relative Fourier amplitude in the highest modes (the upper third by default);
refine the grid when they grow. These indicators depend on the state and the
chosen box, so check that observables converge as the box and grid are enlarged.

`probability(W, grid; q = (a, b))` gives a genuine position-window probability,
and `p = (a, b)` gives a momentum-window probability. Specifying both gives a
phase-space quasi-probability, which can be negative. Window limits need not be
grid points: the integral uses the grid's trigonometric interpolant. Omitting an
axis integrates its whole box, and limits outside the box are clipped.

The rate functions apply the same linear observables to `∂W/∂t`, so they are exact
for the semi-discrete equations. For a state that decays at the box boundary,
`probability_rate(W, op; q = (q_surface, Inf))` gives the reactive flux into the
product region. For a resolved state, spectral derivatives give agreement with the
interpolated probability current at the dividing surface to spectral accuracy;
finite differences approximate that continuum relation as the grid is refined.

For example, analyse transfer across `q = 0` in a tilted double well:

```julia
using HEOM, OrdinaryDiffEqVerner

mass, hbar = 1.3, 0.7
V(q) = 0.1q^4 - 0.5q^2 + 0.1q
grid = PhaseSpaceGrid((-7.0, 7.0), 64, (-7.0, 7.0), 64)
W0 = on_grid(
    (q, p) -> coherent_wigner(q, p; q0 = -1.0, p0 = 0.5, mass, omega = 1.0, hbar),
    grid,
)
op = wigner_moyal_operator(grid; mass, potential = V, hbar)
prob = wigner_moyal_problem(W0, (0.0, 0.2), op)
sol = solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12, saveat = 0.001)

d = diagnostics(sol; potential = V)
product_population = [probability(W, grid; q = (0.0, Inf)) for W in sol.u]
product_flux = [probability_rate(W, op; q = (0.0, Inf)) for W in sol.u]
survival = d.autocorrelation
norm_drift = maximum(abs, d.norm .- first(d.norm))
energy_drift = maximum(abs, d.energy .- first(d.energy))
```

### Tunnelling rates from population relaxation

`tunnelling_rates` estimates constant forward and backward transfer rates when
reactant and product populations obey a two-state kinetic equation:

```math
\dot P_P = k_f(1-P_P)-k_bP_P,\qquad
P_P(t)=P_P^{\mathrm{eq}}+A e^{-\lambda(t-t_0)},\qquad
\lambda=k_f+k_b.
```

The helper fits `log(abs(P_P - equilibrium_product))` by unweighted linear least
squares in a chosen time window. It returns `relaxation_rate = λ`,
`forward_rate = equilibrium_product * λ` and
`backward_rate = (1 - equilibrium_product) * λ`, in inverse time units. This
two-state description is consistent with the population rate equation in
[Lindoy, Mandal and Reichman, Methods, Eqs. 16–18](https://www.nature.com/articles/s41467-023-38368-x).

For a saved HEOM trajectory in a symmetric, undriven double well, with reactant
at `q < 0` and product at `q > 0`, equilibrium symmetry supplies
`equilibrium_product = 0.5`:

```julia
rates = tunnelling_rates(sol;
    equilibrium_product = 0.5, dividing_surface = 0.0,
    product_side = :right, tspan = (10.0, 20.0))
rates.forward_rate                      # reactant → product
rates.backward_rate                     # product → reactant
rates.r_squared                         # goodness of the log-population fit
rates.rmse                              # population-space residual RMS
rates.times, rates.population, rates.fitted_population
```

The helper integrates the physical HEOM root over the product half of the box;
`product_side = :left` selects the opposite side. It accepts isolated Wigner and
Caldeira–Leggett trajectories as well. For saved states use
`tunnelling_rates(states, op; times, equilibrium_product, ...)`. State norms must
be within `norm_atol = 1e-6` of one; the helper does not renormalise them and
rejects time-dependent operators. Already integrated populations can be fitted
with `tunnelling_rates(times, product_population; equilibrium_product, tspan)`.

Supply the equilibrium population independently. For an asymmetric, bath-coupled
system, use the product probability of a converged correlated equilibrium from
`equilibrate`, with the same operator and dividing surface. An isolated Gibbs
state need not have the coupled equilibrium population. The final sample of a
short trajectory is not an equilibrium estimate.

Choose `tspan` after intrawell relaxation and initial bath slip, but before the
population difference reaches numerical noise. The inclusive window must contain
at least three samples with deviations of one sign, all larger than
`min_deviation = 1e-8`, and a positive fitted decay rate. Samples are not silently
removed. With `tspan = nothing`, all supplied samples are fitted. The returned
`amplitude` is the signed fitted deviation at the first retained time. Compare
the fitted populations with the trajectory, vary both window endpoints, and
converge the grid, box, bath expansion and hierarchy depth. A high `r_squared`
alone does not establish a physical rate.

These are total interwell transfer rates. Population relaxation does not separate
under-barrier tunnelling from thermally activated passage over the barrier.
Coherent tunnelling oscillations generally do not admit constant two-state rates;
use populations and fluxes to describe that regime. See the
[runnable double-well example](examples/tunnelling_rates.jl) for a low-temperature
HEOM calculation with `kT = 0.1`, barrier height 1, and `hbar = 0.75`. Its default
preparation is close to thermal equilibrium: first solve the stationary HEOM for
the weakly tilted potential `V(q) + bias*q`, with `bias = 0.002`, then remove the
tilt at `t = 0`. The positive tilt slightly reduces the right-well population
(to approximately 0.48939). Every auxiliary is retained on release, preserving
system–bath correlations. The example uses a Padé bath expansion (`pade = 2`,
`depth = 4`) and the same weak coupling throughout preparation and propagation.
It saves a figure, population data and a text report:

```sh
julia --project=/path/to/environment examples/tunnelling_rates.jl [output_dir]
```

Set `GKSwstype=100` on headless machines. Without an output path, files go to a
temporary directory printed by the script. The default run takes a few minutes
and its stationary preconditioner needs about 2 GB of memory. The companion
[`double_well_equilibrium.jl`](examples/double_well_equilibrium.jl) solves the
trace-constrained stationary equations and checks the full HEOM residual before
release. This accelerator is intended for this small example. It uses odd grid
sizes (`points = 49`) to avoid conserved Nyquist modes of even spectral grids;
its dense block factors grow rapidly with grid size and hierarchy depth.

Repeat with `double_well_trajectory(; bias=0.001)` and compare
`(P_right(t) - 0.5)/(P_right(0) - 0.5)` to check the small-perturbation regime.
Near-equilibrium preparation does not guarantee a constant rate: coherent modes
can still oscillate through equilibrium. With the default parameters, the right
population peaks near 0.50843 at `t = 20.5` and crosses 0.5 twice through `t = 40`.
Halving the bias changes the normalized response by at most approximately
`5.1e-6` on this interval, while preserving those oscillations. A single exponential
`P_right(t) - 0.5 = A*exp(-lambda*t)` would give equal directional rates
`k_left_to_right = k_right_to_left = lambda/2`. Require stable fits across
post-transient time windows. The plotted isolated-doublet cosine is a frequency
reference, not an exact trajectory for the coupled thermal preparation.

After including the file, vary `depth=5`, `pade=3`, odd `points=65`, `qextent`,
`pextent`, and solver tolerances separately to check convergence. Small boundary
weights and Fourier tails alone do not establish hierarchy or bath-expansion
convergence. To recover the strongly displaced pure preparation, use
`double_well_trajectory(; preparation=:localized_doublet)`; that option starts
with zero auxiliaries and has a bath transient. This cold setup also changes
Planck's constant and coupling from the earlier warm example; it is not a
controlled comparison varying temperature alone.

## Plotting

HEOM provides lightweight RecipesBase recipes; `Plots` is optional and does not
load when you use HEOM for computation alone. Install plotting and solver packages
in a separate consumer environment, for example:

```julia
using Pkg
Pkg.activate("heom-examples")
Pkg.develop(path = "/path/to/HEOM.jl")
Pkg.add(["Plots", "OrdinaryDiffEqVerner", "LinearSolve", "SciMLBase"])
```

After running the quick-start simulation in that environment:

```julia
using Plots

fig = wignerplot(sol)                     # final saved state
wignerplot(sol, 1; title = "Initial state")
wignerplot(sol.u[end], grid)              # equivalent matrix-and-grid form
marginalplot(sol)                         # position and momentum density panels
diagnosticsplot(sol; potential = V)       # norm, energy, purity and negativity
diagnosticsplot(d; fields = (:norm, :energy))
diagnosticsplot(sol.t, diagnostics(sol.u, grid; mass, potential = V);
    fields = :purity)
savefig(fig, "wigner.png")
```

`wignerplot` preserves negative values and uses a diverging `:RdBu` colour map
with symmetric colour limits centred on zero. `marginalplot` integrates over the
other axis using the grid quadrature. Neither plot renormalises the data. Both
accept an integer saved-state index for a solution and default to its final state.
`diagnosticsplot` accepts a symbol, tuple or vector of symbols in `fields`; its
default is `(:norm, :energy, :purity, :negativity)`.

All three accept standard Plots attributes such as `size`, `title` and `linewidth`,
and have `!` forms for adding to an existing plot. Repeat custom attributes on `!`
calls to keep them in place, and use the same diagnostic field order when adding
curves to a diagnostic plot. Set `clims` explicitly to compare Wigner heatmaps on
the same colour scale.

## Animations

With `Plots` loaded, `wigneranimation` and `marginalanimation` record saved states
as standard `Plots.Animation` objects. Export them with Plots' `gif` or `mp4`:

```julia
using Plots

# Uniform saved times give constant-speed playback of the physical trajectory.
sol = solve(prob, Vern9(); abstol = 1e-12, reltol = 1e-12,
    saveat = range(prob.tspan...; length = 101))
animation = wigneranimation(sol; size = (600, 500))
gif(animation, "wigner.gif"; fps = 20)
mp4(animation, "wigner.mp4"; fps = 20)

marginals = marginalanimation(sol; indices = 1:2:length(sol.u))
gif(marginals, "marginals.gif"; fps = 10)

# Explicit trajectories also work; omit times to label frames by state index.
wigneranimation(sol.u, grid; times = sol.t, clims = (-0.3, 0.3))
```

By default, Wigner animations use one symmetric colour scale across all selected
states. Marginal animations fix each panel's density limits across those states,
including zero and negative densities. Neither helper clips or renormalises data.
Both add time labels and accept ordinary Plots attributes to override defaults,
including `title`, `clims` for Wigner heatmaps, or `ylims` for marginal densities.

`indices` selects frames in playback order and supports subsets, reversal and
repetition. All selected states are validated before rendering. Each selected
state becomes one frame with equal playback duration; the helpers do not
interpolate irregular saved times. `fps` belongs to `gif`/`mp4`, not the animation
helper. Frame PNGs are stored in Plots' temporary directory. The animation
extension loads only when both HEOM and Plots are loaded.

For a complete runnable example, see [the displaced harmonic oscillator](examples/displaced_sho.jl).
It starts a coherent state at `x = 2`, `p = 0`, checks its motion over one period,
and writes Wigner and marginal GIFs. Run it in the consumer environment above;
an optional command-line argument selects the output directory.

For bath-coupled dynamics, [the animated HEOM oscillator](examples/animated_heom_sho.jl)
starts the same displaced coherent state in a Drude–Lorentz thermal bath and shows
its Wigner distribution, position density, damped centroid motion and oscillator
energy over 20 time units. It checks the full distribution against the Gaussian
solution of the retained bath expansion and writes GIF/MP4 animations, snapshots,
observables and validation results. Run it from the repository in the same
consumer environment, with an optional output directory:

```sh
julia --project=/path/to/heom-examples examples/animated_heom_sho.jl [output_dir]
```

[The Morse oscillator example](examples/heom_morse.jl) follows a displaced quantum
wavepacket in an asymmetric Morse well coupled to a Drude–Lorentz bath. It retains
the full spectral Moyal operator, compares hierarchy depths 6 and 8, and checks
normalization, boundary weight and Fourier tails. It writes GIF/MP4 animations of
the Wigner distribution, position density and potential, along with diagnostics,
CSV observables and validation results:

```sh
julia --project=/path/to/heom-examples examples/heom_morse.jl [output_dir]
```

## Numerical method

**Discretisation.** `discretization = Spectral()`, the default, uses Fourier
pseudo-spectral derivatives. It converges spectrally for smooth Wigner functions,
and its static right-hand side allocates nothing. `discretization = FiniteDifference(order)`
uses central differences of even accuracy `order` (default 4). The Wigner–Moyal and
Caldeira–Leggett finite-difference operators are available as `sparse(op)`; the
Wigner–Moyal part is exactly skew-symmetric, while Caldeira–Leggett damping adds
dissipative terms. HEOM applies the hierarchy couplings directly. Both discretisations
treat the box as periodic, so make the grid large enough that `W` decays to zero at
its edges.

**Moyal series.** `moyal_terms = nothing`, the default, keeps every order. In
momentum Fourier space, where `∂/∂p → iκ`, the potential term is applied exactly as
`i[V(q + ħκ/2) − V(q − ħκ/2)]/ħ`. This needs no derivatives of `V`, but evaluates it
up to `πħ/(2dp)` beyond the box. Only `Spectral()` supports this. `moyal_terms = N`
keeps only the first `N` terms (1 ≤ N ≤ 4), with derivatives of `V` from ForwardDiff:

- `N = 1` is classical Liouville dynamics.
- `N = 2` adds `−(ħ²/24) V‴ ∂³W/∂p³`, which is exact for potentials up to quartic.

**Units.** `hbar` defaults to `1`, as in atomic units. Use any consistent unit system.

**Solvers.** Explicit methods such as `Tsit5()`, `Vern7()` or `Vern9()` work for all
three evolution models. HEOM additionally provides exact matrix-free or sparse
Jacobians for stiff solvers, as described above. The standalone Wigner and CL
problem constructors do not provide these Jacobians. The spectral right-hand side
runs on FFTW and does not accept dual numbers. Strongly anharmonic potentials make the problem
stiff, because the potential symbol grows like `ħ² V‴ κ³`. An operator holds work
buffers, so give each parallel task (for example in an `EnsembleProblem`) its own.
Thermal diffusion also restricts explicit time steps as momentum resolution or
`mass * friction * kT` increases.

**Validation.** For the harmonic oscillator, $V''' = 0$, the Moyal series stops at
the classical term, and every Wigner function rotates rigidly in phase space. The
tests compare against these analytic solutions:

- coherent states;
- the stationary Fock state $|1\rangle$, which has negative values;
- a cat state with interference fringes.

They also check that norm, energy and purity are conserved. The quantum terms are
checked against closed forms instead:

- For a quartic double well, the series terminates at ħ², so the exact operator and
  the `N = 2` truncation must agree.
- For a sine potential, every order contributes, and the exact operator is a pair
  of momentum shifts.
- Correlated Gaussian moments, harmonic-state overlaps and Fock-state negativity
  check the analysis formulas against closed forms.
- Fourier modes test boundary and spectral diagnostics, and exact window integrals.
- Harmonic trajectories test the summary fields and survival probability; population
  fluxes, Ehrenfest relations and energy rates test both discretisations.
- Damped harmonic Gaussians test Caldeira–Leggett evolution against the exact
  time-dependent mean, covariance and full state, including relaxation to thermal
  equilibrium at the bottom of the well.
- Free and uniformly accelerated wavepackets test exact spreading and force signs.
- An analytical sextic eigenstate, transformed directly from its wavefunction,
  tests stationarity through the ħ⁴ Moyal term.
- A non-Markovian Gaussian Langevin solution tests the full Drude HEOM distribution
  as depth increases. Accurate first and second moments alone do not establish
  convergence of the full state.
- A damped Fock state tests the decay of non-Gaussian structure, and the Drude
  fluctuation–dissipation spectrum independently tests Matsubara convergence.

The [Wigner–HEOM audit](references/wigner_heom_audit.md) records the literature
conventions, repaired defects, analytical references, measured errors, and remaining
convergence limits.

## Development

Install the separate formatting and coverage tools once:

```bash
julia --project=build_tools -e 'using Pkg; Pkg.instantiate()'
```

Run the local checks from the repository root:

```bash
julia --project=build_tools build_tools/format.jl --check
julia --project=. -e 'using Pkg; Pkg.test()'
```

`Pkg.test()` installs the test-only dependencies, including `OrdinaryDiffEqVerner`
and `Plots`. It runs the numerical and plotting tests, selected type-inference and
allocation checks, and Aqua's package quality checks. Aqua checks issues such as
method ambiguities, undefined exports, stale dependencies, and missing compatibility
bounds.

Apply formatting changes with:

```bash
julia --project=build_tools build_tools/format.jl --fix
```

Run the coverage gate used by CI:

```bash
julia --project=build_tools -e 'using Coverage; foreach(clean_folder, ("src", "ext"))'
julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'
julia --project=build_tools build_tools/coverage.jl
```

The gate requires 100% coverage of executable lines in `src/` and `ext/` and writes
`lcov.info`. Julia reports line coverage, so this is not the Python template's
branch-coverage metric. Clearing previous coverage first prevents stale runs
from hiding missing tests. Coverage reports do not require a hosted service or
secret token.

Install the optional Git hooks once, then run them on all files:

```bash
pre-commit install
pre-commit run --all-files
```

The hooks check file hygiene and Julia formatting. They require `julia` on your
`PATH` and the tooling environment instantiated as above. See
[build_tools/README.md](build_tools/README.md) for environment maintenance.

CI also tests loading this package from a fresh consumer environment. Julia
packages are distributed as source; there is no wheel-building step.

## Project layout

```text
.
├── .github/workflows/ci.yml       # formatting, tests, coverage, consumer smoke test
├── build_tools/                   # separate development tools and scripts
├── ext/HEOMPlotsExt.jl             # optional Plots animation implementation
├── src/
│   ├── HEOM.jl                    # package module and public exports
│   ├── grid.jl                    # periodic phase-space grid
│   ├── derivatives.jl             # spectral and finite-difference discretisations
│   ├── wigner_moyal.jl            # Wigner–Moyal operators and ODEProblem
│   ├── caldeira_leggett.jl        # thermal friction and diffusion operators
│   ├── heom.jl                    # exponential baths and Wigner-space hierarchy
│   ├── brownian_bath.jl           # thermal underdamped Brownian oscillators
│   ├── composite_bath.jl          # independent bath combinations
│   ├── bath_compression.jl        # balanced reduction of exponential baths
│   ├── bath_diagnostics.jl        # correlations, spectra and harmonic equilibrium
│   ├── observables.jl             # marginals, moments, energy, overlaps, negativity
│   ├── diagnostics.jl             # grid health and state/trajectory summaries
│   ├── populations.jl             # window populations, currents and rates
│   ├── tunnelling.jl              # two-state rates from population relaxation
│   ├── plotting.jl                # optional Plots interface through RecipesBase
│   ├── animation.jl               # public animation API and documentation
│   ├── initial_states.jl          # wavefunction and density-matrix Wigner transforms
│   ├── stationary_states.jl       # numerical eigenstates and finite-box Gibbs states
│   └── harmonic_oscillator.jl     # analytic harmonic-oscillator states and evolution
├── test/                          # numerical and package quality tests
├── examples/initial_states.jl      # double-well tunnelling and Morse state preparation
├── examples/heom_sho.jl            # bath-coupled oscillator and exact centroid check
├── .JuliaFormatter.toml           # formatting rules
└── Project.toml                   # package metadata, compatibility, test dependencies
```

## Development starting point

The scaffold comes from `LouieSlocombe/template_julia` at commit
`98142bfc232b5efc2d19d7362a36c7d6a9b05b3e`, adapted to the `HEOM` package name
and existing package UUID.

Keep the quality checks passing as functionality is added. Add runtime
dependencies with `Pkg.add` and give them explicit compatibility bounds that
still resolve on the oldest supported Julia. Update the supported Julia versions
in `[compat]` and the CI matrix together.

Manifests are ignored because this is a reusable library: its tests
resolve dependencies within the declared compatibility bounds. If you turn it
into an application that needs a fixed environment, commit the relevant
`Manifest.toml` and update `.gitignore` accordingly.

## License

Released under the [MIT License](LICENSE).
