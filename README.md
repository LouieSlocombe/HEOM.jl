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

## Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start), with the package [conventions](#conventions)
- [Choosing an evolution model](#choosing-an-evolution-model)
- [Driven dynamics and spectroscopy](#driven-dynamics-and-spectroscopy)
- [General initial states](#general-initial-states)
- [Wigner-space HEOM](#wigner-space-heom): baths, equilibrium, low temperature,
  stability, compression and spectral-density fits
- [Caldeira–Leggett damping](#caldeiraleggett-damping)
- [Examples](#examples)
- [API](#api)
- [Analysis](#analysis), including tunnelling rates
- [Plotting](#plotting) and [Animations](#animations)
- [Numerical method](#numerical-method)
- [Scope and limitations](#scope-and-limitations)
- [Troubleshooting](#troubleshooting)
- [References and background notes](#references-and-background-notes)
- [Development](#development), with contributing notes and exploratory prototypes
- [Project layout](#project-layout)

## Requirements

- Julia 1.10 or newer in the 1.x series. Install it with
  [Juliaup](https://julialang.org/downloads/).
- An ODE solver package in your environment, for example `OrdinaryDiffEqVerner`,
  `OrdinaryDiffEqTsit5`, or the full `OrdinaryDiffEq`.
- Optional: `Plots` for Wigner heatmaps, marginal densities, diagnostic plots and animations.
- Optional: [pre-commit](https://pre-commit.com/#installation) for Git hooks.

CI runs the test suite on Julia 1.13 on Linux only. The declared lower bound of
Julia 1.10 in `Project.toml` is not exercised in CI.

## Installation

HEOM.jl is not registered. Add it to a project environment directly from GitHub,
together with an ODE solver package:

```julia
using Pkg
Pkg.add(url = "https://github.com/LouieSlocombe/HEOM.jl")
Pkg.add("OrdinaryDiffEqVerner")
```

To track local changes, or to work on the package itself from another project, use
`Pkg.develop(path = "/path/to/HEOM.jl")` in that project's environment instead. From
a clone, instantiate and precompile the package's own environment with:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

That environment holds only the runtime dependencies. Solver and plotting packages
belong in the consumer environment that runs your scripts; [Plotting](#plotting)
shows a complete example environment.

## Quick start

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

Every public function has a docstring available through Julia's help mode, for
example `?wigner_moyal_operator`.

### Conventions

- **Arrays.** A Wigner function is an `nq × np` matrix with `W[i, j] = W(q[i], p[j])`:
  position along the first dimension, momentum along the second. HEOM states are
  `nq × np × nados` arrays with the physical Wigner function in the first slice.
- **Grid.** `PhaseSpaceGrid(qlims, nq, plims, np)` is uniform and periodic on the
  half-open box `[qmin, qmax) × [pmin, pmax)`. Its fields are `q`, `p`, `dq` and `dp`,
  and `size(grid)` is `(nq, np)`. Every integral uses the uniform grid quadrature
  `sum(W) * dq * dp`.
- **Normalisation.** A physical state has `phase_space_integral(W, grid) == 1`.
  Nothing is renormalised during propagation or analysis, so drift stays visible.
- **Units.** `hbar` defaults to `1`. Mass, energy, time and `kT = k_B T` must share one
  consistent unit system; the package performs no unit conversion.
- **Potentials.** A one-argument callable `q -> V(q)` is static. A two-argument callable
  `(q, t) -> V(q, t)`, a `TimeDependentPotential` or a `DrivenPotential` is
  time-dependent and is evaluated at the solver's current time.
- **Baths.** Bath correlations are `C(t) = ⟨B(t)B(0)⟩` for the coupling `H_SB = q B`.
  The Drude–Lorentz spectral density is `J(ω) = 2λγω/(ω² + γ²)`. Thermal constructors
  add their counterterm `λq²` internally; do not add it to `V` yourself.
- **Solutions.** Every problem is an ordinary SciML `ODEProblem`. `sol.u[i]` is the state
  at `sol.t[i]`, and `physical_wigner(sol, i)` views a hierarchy's physical slice.

## Choosing an evolution model

The package provides three semi-discrete generators. They share the grid, potentials,
discretisations, observables, diagnostics and plotting.

| Model | Problem builder | Environment | State | Main cost |
|---|---|---|---|---|
| Wigner–Moyal | `wigner_moyal_problem` | None. Exact closed-system quantum dynamics with the full Moyal series, or a truncation of it | `nq × np` | Grid size; two FFT pairs per evaluation |
| Caldeira–Leggett | `caldeira_leggett_problem` | Markovian Ohmic bath at high temperature, `kT ≫ hbar*omega`, through friction and momentum diffusion | `nq × np` | Grid size; diffusion stiffens fine momentum grids |
| HEOM | `heom_problem` | Non-Markovian Gaussian bath with an exponential correlation expansion, at any positive temperature and coupling strength | `nq × np × nados` | `binomial(modes + depth, depth)` auxiliaries, each a full Wigner function |

Start with the Wigner–Moyal problem to converge the grid and box for the isolated
system. Use Caldeira–Leggett as a quick damped reference when the bath is hot and
memoryless. Use HEOM when bath memory, low temperature, strong coupling or system–bath
correlations matter, and converge its depth and bath expansion separately. Nothing in
the HEOM path is restricted to harmonic systems: every auxiliary uses the same
Wigner–Moyal operator as the isolated problem, including the exact Moyal series.

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

Estimate the size and memory of a hierarchy before building it:

```julia
members = hierarchy_size(length(bath.rates), depth)   # BigInt, root included
bytes_per_state = 8 * prod(size(grid)) * members      # one Float64 hierarchy state
op = heom_operator(grid; mass, potential = V, bath, depth, max_ados = 20_000)
```

An explicit solver holds several state-sized work arrays, and `saveat` multiplies the
storage by the number of saved times. `max_ados` makes `heom_operator` throw instead of
allocating a hierarchy larger than intended.

See [the HEOM harmonic oscillator example](examples/heom_sho.jl), which checks
the centroid against an exact bath-memory equation, compares two hierarchy depths,
and prints the norm drift. Run it with `julia --project=/path/to/environment
examples/heom_sho.jl` in an environment containing HEOM and OrdinaryDiffEqVerner.
The physical model is discussed in
[Tanimura's Wigner-space HEOM derivation](https://arxiv.org/pdf/1502.04077).
The equation above uses coordinate-coupled, unscaled ADOs; their definitions and
initial conditions differ from the transformed momentum form in that paper.

### Bath constructors at a glance

Every constructor returns an `ExponentialBath`. The hierarchy mode count is
`length(bath.rates)`; a damped oscillation occupies two real modes.

| Constructor | Describes | Modes | When to use |
|---|---|---|---|
| `ExponentialBath(c, ν; ...)` | Any `C(t) = Σ cₖ exp(-νₖ t)`, or a general real basis through `weights` and `mixing` | As supplied | Custom or externally fitted decompositions; thermal consistency is not checked |
| `drude_lorentz_bath` | Drude–Lorentz `J(ω)` with `K` Matsubara poles and an optional diffusion terminator | `1 + K` | Warm baths, where few Matsubara poles are needed; the cheapest thermal constructor |
| `drude_lorentz_pade_bath` | Drude–Lorentz `J(ω)` with an `[N/N]` Bose Padé expansion | `1 + N` | Cold baths; far fewer modes than Matsubara at the same accuracy |
| `brownian_oscillator_bath` | Underdamped vibrational resonance with `K` thermal poles | `2 + K` | Structured environments with a discrete mode; its Matsubara coefficients are negative |
| `aaa_bath` | Any `J(ω)`, through an AAA rational fit of the thermal noise spectrum | Chosen by the fit | Spectral densities without a closed-form decomposition, or very cold baths |
| `combine_baths` | Independent components coupled to the same position | Sum of the components | A background plus a resonance, or several fitted bands |
| `compress_bath` | Balanced reduction of any of the above | Chosen | Reduce the mode count before building a deep hierarchy |

Inspect any bath with `bath_correlation`, `bath_spectrum` and `harmonic_covariance`
before using it, and estimate the hierarchy with `hierarchy_size`.

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

Keep such propagations short. In this cold, strong regime the hard cutoff's growing
box-edge modes corrupt the moments within a few time units and diverge by t ≈ 6–13,
depending on depth. See [Long propagations and the box-edge
instability](#long-propagations-and-the-box-edge-instability) for the mechanism, the
`hierarchy_stability` indicator and how to stop a diverging run.

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

### Long propagations and the box-edge instability

The hard depth cutoff is not stable over long times on a finite box. The upward
coupling `2imag(cₖ)q/hbar` grows linearly with `|q|`, so a fixed depth represents the
bath faithfully only where `|q|` times the momentum wavenumber is small compared with
the depth and the bath rates. Beyond that radius the truncated hierarchy has growing
modes. On a periodic box they concentrate near the edges `|q| ≈ L`, and the tails of
the physical state seed them even when it never approaches the edge. The growth rate
rises with coupling strength, depth and box half-width, affects both discretisations,
and depends only weakly on momentum resolution. Amplitude scaling is a diagonal
similarity transformation and leaves it unchanged.

Estimate the growth before a long run, and stop a diverging solve with `unstable_check`:

```julia
s = hierarchy_stability(op)
s.rate                      # largest frozen-coefficient growth rate; 0 means no local growth
s.radius                    # smallest |q| at which any wavenumber is unstable, or Inf
s.position, s.wavenumber    # where the largest rate occurs

limit = 1e3 * maximum(abs, W0)
sol = solve(heom_problem(W0, (0.0, 30.0), op), Vern7();
    abstol = 1e-9, reltol = 1e-9,
    unstable_check = (dt, u, p, t) -> !(maximum(abs, u) < limit))
sol.retcode                 # ReturnCode.Unstable once the limit is reached
```

`rate` is a heuristic, not a bound. In the package tests and in the box-edge study it
was never below the spectral abscissa of the complete generator, but it can lie well
above it, because the kinetic term carries modes out of a narrow unstable band before
they grow. A positive rate means the hard cutoff may diverge over times of order
`1/rate`; in the cases studied the root's second moments were wrong by `1e-3` after 6 to
20 multiples of `1/rate`. Warm, weak baths on moderate boxes often stay stable for a
whole run, and not oversizing the box helps. Confirm any long-time result with a
different box and depth, and watch the physical state's boundary weights.

For scale, the cold Padé bath above (`λ = 0.8`, `γ = 0.5`, `kT = 0.1`, `pade = 2`) with a
unit harmonic oscillator on a ±8, 64-point grid has second moments wrong by `1e-2` at
`t = 5` at depth 2, and exceeds `10³` times its initial amplitude at `t ≈ 12.5`, `7.6`
and `5.9` for depths 2, 4 and 6. A warm bath (`λ = 0.2`, `γ = kT = 1`, `pade = 1`)
reaches an error of `1e-3` near `t = 10` and diverges near `t = 40`. The test suite
records the cold case as a known defect with `@test_broken`. No remedy evaluated so far
repairs the cold, strong regime; see [Exploratory prototypes](#exploratory-prototypes).

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
For spectral densities without a thermal constructor, `aaa_bath` below fits the source.

### Fitting general spectral densities

`aaa_bath` builds a thermal bath for any spectral density `J(ω)` by fitting the noise
spectrum `S(ω) = 2hbar*J(ω)/(1 - exp(-hbar*ω/kT))` with the AAA rational algorithm, as
in free-pole HEOM
([Xu et al., Phys. Rev. Lett. 129, 230601 (2022)](https://doi.org/10.1103/PhysRevLett.129.230601)).
The fit chooses its own decay rates instead of the poles of `J` and the Matsubara
frequencies, so a cold bath needs far fewer modes:

```julia
drude(ω) = 2 * 0.3 * 0.5 * ω / (ω^2 + 0.5^2)
vibration(ω) = 2 * 0.2 * 0.2 * 1.5^2 * ω / ((1.5^2 - ω^2)^2 + 0.2^2 * ω^2)
J(ω) = drude(ω) + vibration(ω)

bath = aaa_bath(J; kT = 0.1, hbar = 1, reltol = 1e-6,
                frequencies = exp10.(range(-3, 3; length = 600)))
modes = length(bath.rates)                # 12 real hierarchy modes
covariance = harmonic_covariance(bath; mass = 1, omega = 1)
reduced = compress_bath(bath; modes = 8)  # optional balanced reduction
```

`J` is called only at the positive sample `frequencies`, and `J(-ω) = -J(ω)` is
implied. The correlation convention is that of the Drude constructor, and the
counterterm is the reorganization `λ = (1/π)∫₀^∞ J(ω)/ω dω`. It is computed by
quadrature unless `reorganization = λ` is passed, which is needed when `J(ω)/ω` is
too singular or too sharply peaked to integrate numerically. The fit raises its
degree until `|bath_spectrum(bath, ±ω) - S(±ω)| ≤ reltol*maximum(S)` at every sample,
and throws beyond `max_modes = 40` modes. Sample the thermal scale `kT/hbar` near zero,
every feature of `J`, and its decay; the fit is not controlled outside the sampled
range.

The even part of `S` and its odd part divided by `ω` are fitted together as rational
functions of `ω²`. All rates therefore decay, and complex rates come in conjugate
pairs. A real rate is one hierarchy mode. A damped oscillation is a real two-mode
block, like `brownian_oscillator_bath`, so the vibration above uses two of the twelve
modes and `length(bath.rates)` is always the hierarchy mode count. The coefficients
are refitted by least squares on the samples. A nonnegative constant remainder of the
even spectrum becomes white-noise `diffusion`.

For the cold Drude bath of the previous sections (`λ = 0.8`, `γ = 0.5`, `kT = 0.1`),
eight fitted modes at `reltol = 1e-4` give a harmonic equilibrium within a relative
`1.1e-6` of the continuum result, compared with `3.7e-4` for the eight-mode Padé bath.
Sixteen modes at `reltol = 1e-8` reach `1.1e-8`. For very few modes, fit the band the
system resolves rather than compressing a much wider fit: balanced truncation weighs
all frequencies equally.

As for compression, a small spectral residual does not establish thermalization.
`J(ω)/ω` below the lowest sample is missing from the fit, while the counterterm
contains it. The resulting static imbalance, `bath.counterterm + imag(∫₀^∞ C dt)/hbar`
with `∫₀^∞ C dt = transpose(weights) * ((Diagonal(rates) + mixing) \ coefficients)`,
stiffens the equilibrium. For a sub-Ohmic `J ∝ sqrt(ω)exp(-ω/2)` at `kT = 0.1`,
fits sampled down to `ω = 1e-2` and `1e-4` both meet `reltol = 1e-4`. Their harmonic
equilibria differ from the continuum result by a relative `3.8e-3` and `3.5e-4`. Check
`harmonic_covariance`, extend the sampled range, and repeat the hierarchy-depth
convergence with the fitted bath.

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

## Examples

Each script in [`examples/`](examples) runs on its own and prints or writes its own
checks. Run one from the repository root with
`julia --project=/path/to/environment examples/<script>.jl`, in an environment that
contains HEOM and the packages listed below. Scripts that write figures accept an
optional output directory as their first argument and need `GKSwstype=100` on
headless machines.

| Script | What it demonstrates | Also needs |
|---|---|---|
| [`displaced_sho.jl`](examples/displaced_sho.jl) | A coherent state over one harmonic period; writes Wigner and marginal GIFs | OrdinaryDiffEqVerner, Plots |
| [`damped_sho.jl`](examples/damped_sho.jl) | Caldeira–Leggett relaxation against the exact Gaussian mean and covariance | OrdinaryDiffEqVerner |
| [`driven_harmonic.jl`](examples/driven_harmonic.jl) | Forced centroid, impulse response and the analytic broadened susceptibility | OrdinaryDiffEqVerner |
| [`anharmonic_spectroscopy.jl`](examples/anharmonic_spectroscopy.jl) | Quartic-oscillator response against an eigenstate Kubo sum; a weak pulse against response convolution | OrdinaryDiffEqVerner |
| [`initial_states.jl`](examples/initial_states.jl) | Double-well tunnelling superposition, wavefunction and density-matrix transforms, Morse ground state | nothing |
| [`heom_sho.jl`](examples/heom_sho.jl) | Drude–Lorentz HEOM oscillator against an exact bath-memory centroid; two depths; norm drift | OrdinaryDiffEqVerner |
| [`animated_heom_sho.jl`](examples/animated_heom_sho.jl) | Displaced oscillator relaxing in a Drude–Lorentz bath, checked against a Gaussian reference; GIF and MP4 output | OrdinaryDiffEqVerner, Plots |
| [`heom_morse.jl`](examples/heom_morse.jl) | Morse wavepacket relaxation at depths 6 and 8, with animations and grid diagnostics | OrdinaryDiffEqVerner, Plots |
| [`low_temperature_strong_coupling.jl`](examples/low_temperature_strong_coupling.jl) | Cold, strongly coupled anharmonic dynamics with separate depth, Padé, spacing and box sweeps; solves seven hierarchies | OrdinaryDiffEqVerner |
| [`bath_compression.jl`](examples/bath_compression.jl) | Balanced compression of an eight-mode Padé bath, with spectral, covariance and depth checks; several minutes | OrdinaryDiffEqVerner |
| [`tunnelling_rates.jl`](examples/tunnelling_rates.jl) | Low-temperature double-well tunnelling rates from a near-equilibrium HEOM preparation; about 2 GB and a few minutes | OrdinaryDiffEqVerner, LinearSolve, SciMLBase, Plots |

Two files are helpers included by the scripts above rather than examples in their own
right. [`heom_gaussian_reference.jl`](examples/heom_gaussian_reference.jl) is the
generalized Langevin Gaussian solution for a finite exponential bath, the same
construction the tests use, and
[`double_well_equilibrium.jl`](examples/double_well_equilibrium.jl) is the
trace-constrained stationary solver used by the tunnelling example.

## API

| Function | Purpose |
|---|---|
| `PhaseSpaceGrid(qlims, nq, plims, np)` | Uniform periodic grid on `[qmin, qmax) × [pmin, pmax)` |
| `on_grid(f, grid)` | Sample `f(q, p)` on the grid |
| `Spectral()`, `FiniteDifference(order = 4)` | Discretisation options for every operator; finite differences need an integer `moyal_terms` |
| `wigner_moyal_problem(W0, tspan, grid; mass, potential, ...)` | `ODEProblem` for the Wigner–Moyal equation |
| `wigner_moyal_operator(grid; mass, potential, hbar, discretization, moyal_terms)` | Reusable semi-discrete operator; `wigner_moyal_problem(W0, tspan, op)` takes it |
| `wigner_moyal!(dW, W, op, t)` | In-place right-hand side |
| `TimeDependentPotential(V; derivative)`, `DrivenPotential(V0, field, dipole; field_derivative)` | General `V(q,t)` or separable `V0(q) - field(t)*dipole(q)` |
| `linear_response(initial, tspan, op, alg; dipole, observable, ...)` | Propagate the impulse response with an undriven generator; returns a `LinearResponseResult` with `times`, `response` and `solution` |
| `linear_response_problem(initial, tspan, op; dipole)` | Response problem seeded by the dipole commutator |
| `absorption_spectrum(response; frequencies, broadening)` | Retarded susceptibility and absorption intensity from saved response data |
| `ExponentialBath(coefficients, rates; hbar, diffusion, counterterm)` | Gaussian bath described by an exponential correlation expansion |
| `drude_lorentz_pade_bath(; reorganization, cutoff, kT, pade, ...)` | Compact quantum Drude bath for low temperatures |
| `hierarchy_size(modes, depth)`, `rescale_hierarchy(U, op; scaled)` | Estimate hierarchy size and convert auxiliary scaling |
| `hierarchy_stability(op)` | Frozen-coefficient growth rate and stability radius of the truncated hierarchy |
| `drude_lorentz_bath(; reorganization, cutoff, kT, matsubara, hbar, terminator)` | Drude–Lorentz bath with Matsubara poles and optional residual diffusion |
| `brownian_oscillator_bath(; reorganization, frequency, damping, kT, matsubara, hbar)` | Underdamped Brownian resonance with retained thermal poles |
| `combine_baths(baths...)`, `combine_baths(baths)` | Sum independent bath components coupled to the same position |
| `bath_correlation(bath, t)`, `bath_spectrum(bath, omega)` | Retained force correlation and unsymmetrized noise spectrum |
| `aaa_bath(J; kT, frequencies, hbar, reltol, max_modes, reorganization)` | Thermal bath for a general spectral density, from an AAA fit of its noise spectrum |
| `compress_bath(bath; modes, method)`, `hankel_singular_values(bath)` | Balanced reduction of an exponential bath to fewer real hierarchy modes |
| `harmonic_covariance(bath; mass, omega)` | Exact harmonic-oscillator equilibrium covariance for a bath decomposition |
| `heom_problem(W0, tspan, grid; mass, potential, bath, depth, ...)` | `ODEProblem` for the Wigner-space hierarchy |
| `heom_operator(grid; mass, potential, bath, depth, ...)` | Reusable hierarchy operator; `heom_problem(W0, tspan, op)` takes it |
| `heom!(dU, U, op, t)` | In-place hierarchy right-hand side |
| `sparse(op)`, `sparse(op, t)` | Full sparse generator of a finite-difference Wigner–Moyal, Caldeira–Leggett or HEOM operator |
| `hierarchy_indices(op)` | Multi-indices corresponding to the third state-array dimension |
| `equilibrate(U0, tspan, op, alg; stationarity_abstol, stationarity_reltol, check_interval, ...)` | Relax to a correlated equilibrium and test stationarity of every member; returns an `EquilibriumResult` |
| `heom_problem(eq, tspan)`, `linear_response(eq, tspan, alg; ...)`, `physical_wigner(eq)` | Restart, excite or view a converged `EquilibriumResult` with all its correlations |
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

| Field | Meaning |
|---|---|
| `norm` | `∫W dq dp`, which should stay at 1 |
| `mean_q`, `mean_p` | Centroid |
| `var_q`, `var_p`, `cov_qp` | Central second moments, with the symmetrised covariance |
| `uncertainty` | σqσp, at least ħ/2 for a normalised physical state |
| `robertson_schrodinger` | √det Σ, at least ħ/2 and invariant under harmonic evolution |
| `energy` | `⟨p²/2m + V⟩`, instantaneous for a driven potential |
| `purity` | `Tr ρ²`, 1 for a pure state |
| `negativity` | `∫(abs(W) − W) dq dp`, zero for a nonnegative Wigner function |
| `boundary_q`, `boundary_p` | Fraction of `∫abs(W)` in the outer 5% of points at both ends of the axis |
| `tail_q`, `tail_p` | Relative Fourier amplitude in the upper third of modes along the axis |
| `autocorrelation` | Trajectories only: `overlap` with the first saved state |
| `t` | Solutions only: the saved times |

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

## Scope and limitations

- **One dimension.** The phase space is a single position `q` and momentum `p`. There
  is no multi-particle or multi-dimensional grid.
- **Linear position coupling.** Baths couple through `H_SB = q B`. Several independent
  components can be summed, but all couple to the same operator; nonlinear or momentum
  coupling is not supported.
- **Periodic box.** Both discretisations treat the box as periodic. States must decay to
  negligible weight at the edges, and there are no absorbing boundaries; watch
  `boundary_q` and `boundary_p` and enlarge the box when they grow.
- **Hard depth cutoff.** The hierarchy sets tiers beyond `depth` to zero. On a finite
  box this truncation has growing modes near the edges, so long HEOM propagations need
  the checks in [Long propagations and the box-edge
  instability](#long-propagations-and-the-box-edge-instability).
- **Positive temperature.** The thermal bath constructors require `kT > 0`. Exactly zero
  temperature needs a different correlation decomposition, which is not provided.
- **High-temperature Caldeira–Leggett.** The friction-and-diffusion model is a Markovian,
  high-temperature approximation and does not cool into the quantum ground state.
- **No positivity correction.** Nothing renormalises the state or enforces positivity
  of the reduced density operator; inaccuracies show up in norm, purity and negativity.
- **Finite expansions.** Bath decompositions, hierarchy depth, grid, box and ODE
  tolerances are separate approximations and must be converged independently. A small
  spectral residual of a fitted or compressed bath does not establish thermalization.
- **Differentiation.** The spectral right-hand side runs on FFTW and does not accept dual
  numbers, so the ODE right-hand side cannot be differentiated by ForwardDiff. With
  `moyal_terms = N`, the potential itself must accept dual numbers.
- **Dense eigenstates.** `eigenstates` and `thermal_wigner` diagonalise a dense `nq × nq`
  Hamiltonian; they are meant for moderate position grids and localised states.
- **Operators are stateful.** Every operator owns FFT or scratch buffers. Construct a
  separate operator for each concurrent solve, for example in an `EnsembleProblem`, and
  do not mutate an operator retained by an `EquilibriumResult`.

## Troubleshooting

- **Norm, energy or purity drift grows in time.** The state is reaching the box edge or
  is under-resolved. Inspect `boundary_q`, `boundary_p`, `tail_q` and `tail_p` from
  `diagnostics`. Enlarge the box when boundary weights grow and refine the grid when
  spectral tails grow, then tighten `abstol` and `reltol`.
- **A HEOM solve stops with `ReturnCode.Unstable` or produces huge values.** The
  hard-cutoff box-edge modes have grown. Check `hierarchy_stability(op)`, shorten the
  time span, avoid oversizing the box, compare a different depth, and locate the
  hierarchy's largest values; see [Long propagations and the box-edge
  instability](#long-propagations-and-the-box-edge-instability).
- **An explicit solver takes very small steps.** The problem is stiff: fast Matsubara
  rates, momentum diffusion on a fine grid, or a strongly anharmonic potential whose
  Moyal symbol grows like `ħ²V‴κ³`. Switch to a Padé or AAA bath with fewer, slower
  modes, or use the stiff path with `Rodas5P` and the Krylov or sparse Jacobian from
  [Low temperature and strong coupling](#low-temperature-and-strong-coupling).
- **`heom_operator` refuses a hierarchy because of `max_ados`.** The requested depth
  and bath produce more members than the limit. Lower `depth`, reduce the mode count
  with `compress_bath` or `aaa_bath`, or raise `max_ados` after checking
  `hierarchy_size`.
- **A finite-difference operator rejects `moyal_terms = nothing`.** The exact Moyal
  operator is nonlocal in momentum and exists only for `Spectral()`. Pass an integer
  `moyal_terms` with `FiniteDifference`, or use the spectral discretisation.
- **`sparse(op)` throws for a driven problem.** The generator changes in time. Use
  `sparse(op, t)` for the instantaneous matrix, or `jacobian = :sparse` in
  `heom_problem`, which updates the matrix itself.
- **`equilibrate` returns `converged = false`.** Read `status`. `:time_limit` means
  the budget ran out, so pass `eq.hierarchy` back in with a longer span.
  `:solver_failure` carries the solver `retcode`, and `:terminated` means another
  callback stopped the run. Loosen the stationarity tolerances only after checking
  `residuals` member by member.
- **`wigneranimation` or `marginalanimation` is undefined.** The animation extension
  loads only when both HEOM and Plots are loaded in the same session.
- **Plotting fails on a headless machine.** Set the environment variable
  `GKSwstype=100` before starting Julia.
- **The potential errors on dual numbers.** `moyal_terms = N` differentiates `V` with
  ForwardDiff, so `V` must be generic in its argument type. Use `moyal_terms = nothing`
  with `Spectral()` to avoid derivatives entirely.

## References and background notes

The [`references/`](references) directory holds implementation notes written for this
package. They record conventions, derivations, measured errors and open problems, and
they separate what each paper states from what the implementation derives.

- [`wigner_heom_audit.md`](references/wigner_heom_audit.md): the audit of the
  Wigner–Moyal, hierarchy and Caldeira–Leggett conventions, the repaired defects and the
  analytical validations behind the test suite.
- [`low_temperature_strong_coupling.md`](references/low_temperature_strong_coupling.md):
  scaled auxiliaries, Padé baths, repeated poles, stiff integration and the cold,
  strong-coupling benchmarks.
- [`cabrera_2015_wigner_implementation.md`](references/cabrera_2015_wigner_implementation.md):
  notes on Cabrera, Bondar, Jacobs and Rabitz, *Phys. Rev. A* **92**, 042122 (2015), the
  spectral Wigner propagation scheme and Caldeira–Leggett terms that the Hamiltonian
  operator follows.
- [`tanimura_2020_heom_implementation_notes.md`](references/tanimura_2020_heom_implementation_notes.md):
  notes on Tanimura, *J. Chem. Phys.* **153**, 020901 (2020), the HEOM review that
  supplies the density-operator hierarchy, its Wigner form and the Brownian-oscillator
  acceptance tests.
- [`stabilized_heom_implementation.md`](references/stabilized_heom_implementation.md):
  notes on Gatto, Rudge, Hou, Rabani and Thoss, arXiv:2609.35484 (2026), a transformed
  bosonic hierarchy with improved truncation stability. It is not implemented; the
  box-edge study found its frozen symbol worse than the package form in the cold,
  strong case.

Methods used by the package, with the links cited in the sections above:

- Wigner-space HEOM and its harmonic benchmarks: [Tanimura, *J. Chem. Phys.* **142**, 144110 (2015)](https://arxiv.org/pdf/1502.04077).
- Factorial and amplitude scaling of auxiliaries: Shi et al., *J. Chem. Phys.* **130**, 084105 (2009), [doi:10.1063/1.3077918](https://doi.org/10.1063/1.3077918).
- Bose Padé decomposition: Hu et al., *J. Chem. Phys.* **134**, 244106 (2011), [doi:10.1063/1.3602466](https://doi.org/10.1063/1.3602466), and [Ding et al., *J. Chem. Phys.* **135**, 164107 (2011)](https://arxiv.org/pdf/1107.0249).
- General real correlation bases with weights and mixing: Ikeda and Scholes, *J. Chem. Phys.* **152**, 204101 (2020), [doi:10.1063/5.0007327](https://doi.org/10.1063/5.0007327).
- Underdamped Brownian correlation decomposition: [doi:10.1038/s41467-019-11656-1](https://doi.org/10.1038/s41467-019-11656-1).
- Truncated Padé decomposition and balanced compression: [Takahashi and Tanimura, *J. Chem. Phys.* **158**, 044115 (2023)](https://doi.org/10.1063/5.0135725); on the limits of spectral residuals, [Tokieda, *Phys. Rev. Research* **7**, 043178 (2025)](https://doi.org/10.1103/bv19-dtb1).
- Free-pole HEOM with AAA fits: [Xu et al., *Phys. Rev. Lett.* **129**, 230601 (2022)](https://doi.org/10.1103/PhysRevLett.129.230601).
- Caldeira–Leggett coefficients: [García-Palacios and Zueco](https://arxiv.org/pdf/cond-mat/0407454). Two-state rate equations: [Lindoy, Mandal and Reichman, *Nat. Commun.* (2023)](https://www.nature.com/articles/s41467-023-38368-x).
- Nakajima–Zwanzig terminators, evaluated and not adopted: [Fay, *J. Chem. Phys.* **157**, 054108 (2022)](https://arxiv.org/abs/2205.09270).

There is no registered release or archived DOI yet. When citing the package, give the
repository URL and the commit used, together with the method papers above.

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

### Contributing

- Add a test for every new public function or convention, and keep the coverage gate
  at 100% of executable lines in `src/` and `ext/`.
- Give every public function a docstring and a row in the [API](#api) table, and
  update [Project layout](#project-layout) when adding a source file.
- Validate numerics against a closed form or an independent reference where one exists;
  the analytic benchmarks in `test/` show the pattern.
- Run the formatter before committing; the pre-commit hook enforces it.

### Exploratory prototypes

[`prototypes/`](prototypes) holds studies that informed the package but are not part of
it. They are not tested, not covered and not loaded by `using HEOM`. Each has its own
environment that develops the package from the repository root, so run a script with
`julia --project=prototypes/<study> prototypes/<study>/<script>.jl`.

- [`nz_terminator`](prototypes/nz_terminator/README.md) asked whether the
  Nakajima–Zwanzig terminator of Fay (2022) or a scalar Markov closure should replace
  the hard depth cutoff. Neither was adopted. The terminator gains at most 2.4× at a
  fixed depth, costs 9–64× per evaluation and destabilises the hierarchy, while one
  extra depth level gains 2–18×. The only useful variant was the scalar closure fed by
  all parents, and the study also uncovered the box-edge instability.
- [`box_edge`](prototypes/box_edge/README.md) characterised that instability across
  baths, depths, boxes, discretisations and scalings, derived the frozen-coefficient
  symbol that became `hierarchy_stability`, and tested absorbers, tapers, closures and
  the transformed hierarchy of Gatto et al. None repaired the cold, strong regime; warm
  and weak baths can be stabilised by the scalar closure or by not oversizing the box.

## Project layout

```text
.
├── .github/
│   ├── dependabot.yml             # monthly GitHub Actions updates
│   └── workflows/ci.yml           # tests, consumer smoke test, formatting, coverage
├── build_tools/                   # separate JuliaFormatter and Coverage environment
├── examples/                      # runnable scripts, listed under Examples
├── ext/HEOMPlotsExt.jl            # Plots-only animation implementation
├── prototypes/                    # exploratory studies, not part of the package
│   ├── box_edge/                  # box-edge instability characterisation and remedies
│   └── nz_terminator/             # final-tier closure benchmarks
├── references/                    # implementation notes and audits
├── src/
│   ├── HEOM.jl                    # package module and public exports
│   ├── grid.jl                    # periodic phase-space grid
│   ├── derivatives.jl             # spectral and finite-difference discretisations
│   ├── wigner_moyal.jl            # Wigner–Moyal operators and ODEProblem
│   ├── driven.jl                  # time-dependent and separable driven potentials
│   ├── caldeira_leggett.jl        # thermal friction and diffusion operators
│   ├── heom.jl                    # exponential baths, Drude–Lorentz bath, hierarchy operator
│   ├── hierarchy_stability.jl     # frozen-coefficient stability indicator
│   ├── pade_bath.jl               # Bose Padé Drude–Lorentz bath
│   ├── brownian_bath.jl           # thermal underdamped Brownian oscillators
│   ├── composite_bath.jl          # independent bath combinations
│   ├── bath_compression.jl        # balanced reduction of exponential baths
│   ├── aaa_bath.jl                # AAA rational fits of general spectral densities
│   ├── bath_diagnostics.jl        # correlations, spectra and harmonic equilibrium
│   ├── heom_solvers.jl            # sparse generator, Jacobian products, ODEFunction
│   ├── equilibrium.jl             # equilibrate and EquilibriumResult
│   ├── observables.jl             # marginals, moments, energy, overlaps, negativity
│   ├── diagnostics.jl             # grid health and state/trajectory summaries
│   ├── populations.jl             # window populations, currents and rates
│   ├── tunnelling.jl              # two-state rates from population relaxation
│   ├── harmonic_oscillator.jl     # analytic harmonic-oscillator states and evolution
│   ├── initial_states.jl          # wavefunction and density-matrix Wigner transforms
│   ├── stationary_states.jl       # numerical eigenstates and finite-box Gibbs states
│   ├── spectroscopy.jl            # linear response and absorption spectra
│   ├── plotting.jl                # optional Plots interface through RecipesBase
│   └── animation.jl               # public animation API and documentation
├── test/                          # one file per source area, plus Aqua
├── .JuliaFormatter.toml           # formatting rules
├── .pre-commit-config.yaml        # file hygiene and formatting hooks
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
