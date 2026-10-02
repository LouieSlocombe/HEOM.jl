# HEOM.jl

Phase-space quantum dynamics in Julia, developed towards hierarchical equations of
motion (HEOM). The package currently solves the Wigner–Moyal equation. This is the
exact quantum evolution of the Wigner function $W(q, p, t)$ of a particle of mass
$m$ in a potential $V(q)$:

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
```

A Wigner function on the grid is an `nq × np` matrix with `W[i, j] = W(q[i], p[j])`.
The solution's states are `sol.u` at times `sol.t`.

To work on this package from another Julia project, use
`Pkg.develop(path="/path/to/HEOM.jl")` in that project's environment. Every public
function has a docstring available through Julia's help mode, for example
`?wigner_moyal_operator`.

## API

| Function | Purpose |
|---|---|
| `PhaseSpaceGrid(qlims, nq, plims, np)` | Uniform periodic grid on `[qmin, qmax) × [pmin, pmax)` |
| `on_grid(f, grid)` | Sample `f(q, p)` on the grid |
| `wigner_moyal_problem(W0, tspan, grid; mass, potential, ...)` | `ODEProblem` for the Wigner–Moyal equation |
| `wigner_moyal_operator(grid; mass, potential, hbar, discretization, moyal_terms)` | Reusable semi-discrete operator; `wigner_moyal_problem(W0, tspan, op)` takes it |
| `wigner_moyal!(dW, W, op, t)` | In-place right-hand side |
| `phase_space_integral`, `expectation`, `purity` | Norm, phase-space averages and `Tr ρ²` |
| `harmonic_potential`, `coherent_wigner`, `fock_wigner`, `cat_wigner` | Harmonic oscillator and its standard states |
| `harmonic_evolution(W0, t; mass, omega)` | Exact harmonic-oscillator solution, for testing |

## Numerical method

**Discretisation.** `discretization = Spectral()`, the default, uses Fourier
pseudo-spectral derivatives. It converges spectrally for smooth Wigner functions,
and its right-hand side allocates nothing. `discretization = FiniteDifference(order)`
uses central differences of even accuracy `order` (default 4). The whole operator is
a sparse, exactly skew-symmetric matrix, available as `sparse(op)`. Both treat the
box as periodic, so make the grid large enough that `W` decays to zero at its edges.

**Moyal series.** `moyal_terms = nothing`, the default, keeps every order. In
momentum Fourier space, where `∂/∂p → iκ`, the potential term is applied exactly as
`i[V(q + ħκ/2) − V(q − ħκ/2)]/ħ`. This needs no derivatives of `V`, but evaluates it
up to `πħ/(2dp)` beyond the box. Only `Spectral()` supports this. `moyal_terms = N`
keeps only the first `N` terms (1 ≤ N ≤ 4), with derivatives of `V` from ForwardDiff:

- `N = 1` is classical Liouville dynamics.
- `N = 2` adds `−(ħ²/24) V‴ ∂³W/∂p³`, which is exact for potentials up to quartic.

**Units.** `hbar` defaults to `1`, as in atomic units. Use any consistent unit system.

**Solvers.** Use an explicit method such as `Tsit5()`, `Vern7()` or `Vern9()`. The
spectral right-hand side runs on FFTW and does not accept dual numbers, and no
sparse Jacobian is wired up yet. Strongly anharmonic potentials make the problem
stiff, because the potential symbol grows like `ħ² V‴ κ³`. An operator holds work
buffers, so give each parallel task (for example in an `EnsembleProblem`) its own.

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

`Pkg.test()` installs the test-only dependencies, including `OrdinaryDiffEqVerner`.
It runs the numerical tests, selected type-inference and allocation checks, and Aqua's
package quality checks. Aqua checks issues such as method ambiguities, undefined
exports, stale dependencies, and missing compatibility bounds.

Apply formatting changes with:

```bash
julia --project=build_tools build_tools/format.jl --fix
```

Run the coverage gate used by CI:

```bash
julia --project=build_tools -e 'using Coverage; clean_folder("src")'
julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'
julia --project=build_tools build_tools/coverage.jl
```

The gate requires 100% coverage of executable lines in `src/` and writes
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
├── src/
│   ├── HEOM.jl                    # package module and public exports
│   ├── grid.jl                    # periodic phase-space grid
│   ├── derivatives.jl             # spectral and finite-difference discretisations
│   ├── wigner_moyal.jl            # Wigner–Moyal operators and ODEProblem
│   ├── observables.jl             # norm, expectation values, purity
│   └── harmonic_oscillator.jl     # analytic harmonic-oscillator states and evolution
├── test/                          # numerical and package quality tests
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
