# Box-edge instability of the hard-cutoff Wigner HEOM

Exploratory study, not part of the package. It characterises the growing modes of
`heom_operator`'s hard depth cutoff, which [`../nz_terminator`](../nz_terminator/README.md)
first found, and evaluates candidate remedies. Run any script with
`julia --project=prototypes/box_edge prototypes/box_edge/<script> [args]`.

## Conclusions

1. **Mechanism.** Every hierarchy coupling is local in position `q` and in the momentum
   wavenumber `κ` (the variable conjugate to `p`). The Moyal potential term is a purely
   imaginary multiplier there, the same on every member. Only the kinetic term `-(p/m)∂q`
   couples different `(q, κ)`. At a frozen `(q, κ)`, the exact hierarchy is a single
   damped mode, stable for every `q`. Its stationary auxiliaries form a coherent state with
   amplitude about `|Φ|/ν ∝ |q|κ/ν`. A fixed-depth truncation reproduces it only while
   that amplitude is small. Beyond a radius the truncated chain has eigenvalues with
   positive real part, and their size grows roughly linearly with `|q|`. On the periodic
   box this region always reaches the edge `|q| = L`, where the growth is fastest.
2. **The frozen symbol predicts the instability.** Its largest real part (now
   `hierarchy_stability(op)` in the package) was never below the dense or power-iteration
   spectral abscissa. It was 1.05–2× the abscissa on small grids and up to 4× for warm or
   weak baths on production grids. It can even be positive when the full generator is
   stable: if the box ends only a little beyond the radius, the kinetic term carries
   modes out of the unstable band before they grow.
3. **Seeding.** The growing modes are seeded by the tails of the physical state, not
   only by rounding at the edge. Moment errors grow at the generator's rate from the
   start: in the cold, strong case, 9e-7 at t = 1, 6e-4 at t = 3 and 1e-2 at t = 5. The
   existing cold, strong transient tests (t ≤ 0.8) sit just inside the safe window.
4. **Every regime tested is affected.** This covers warm, cold weak and cold strong
   Drude baths with Padé or Matsubara poles, with and without the terminator, and an
   underdamped Brownian bath with mode mixing. Both discretisations and both scalings
   are affected. The growth rate increases with depth, coupling strength and box
   half-width. It is nearly independent of momentum resolution, and identical with and
   without amplitude scaling, which is a diagonal similarity transformation. The Padé
   terminator's diffusion `D∂p²` damps the highest wavenumbers; without it the cold,
   strong rate nearly doubles (9.3 instead of 4.9 at depth 4).
5. **No remedy tested fixes the cold, strong regime.**
   - That hierarchy is locally unstable from `|q| ≈ 1`, inside the physical state, at
     every depth tested.
   - Edge remedies (absorber, taper) only cap the growth there.
   - The transformed hierarchy of Gatto et al. has a *worse* frozen symbol than the
     package form, 3× in the cold, strong case.
   - `:markov_full` combined with an absorber makes depth 2 stable, but the moments then
     settle at about 1e-2 error. Depth 4 is still unstable.
6. **Warm and weak baths can be stabilised.**
   - The scalar Markov closure `:markov_full` (issue #37) is the best option tested. It
     makes depths 2–4 stable in the warm and cold, weak cases. At depth 6 it leaves
     small rates (warm 0.30, cold weak 0.07), which an absorber removes. The long-time moment error to t = 60 drops to 2e-4 (warm,
     depth 2), 7e-6 (warm, depth 4) and 2e-11 (cold weak, depth 2), where the hard
     cutoff diverges or loses 10⁻³ accuracy.
   - An auxiliary absorber also stabilises these cases. In warm baths, however, it
     perturbs the moments by 1e-4 to 1e-2, because the auxiliaries carry `|q|ⁿ`-weighted
     tails of the state into the layer.
   - Not oversizing the box also helps. At fixed spacing, warm depth 4 is stable for
     `L ≤ 6`, with growth rates 1.19, 2.65 and 3.98 at L = 8, 10 and 12. A weak, warm
     bath (λ = 0.05, kT = 1) is stable even on ±8.

Package changes made from these findings: the `hierarchy_stability` diagnostic, the
stability notes in the `heom_operator`, `heom_problem` and `equilibrate` docstrings and
the README, and `test/hierarchy_stability.jl`. That file checks the indicator against
dense spectra. It also holds a long-time cold, strong regression marked `@test_broken`,
which a fix will turn into an unexpected pass. Default behaviour is unchanged.

## The frozen-coefficient symbol

At fixed `q` and `κ ≥ 0`, `∂p → s₁(κ)` (`iκ` spectrally, the finite-difference symbol
otherwise) and `∂p² → s₂(κ)`. Each grid point then gives a `members × members` matrix:

    H(q, κ) = transfer − diag(damping) + D s₂ + Σₖ [lowering·s₁ + raising_derivative·s₁ + raising_coordinate·q].

`common.jl` (`local_abscissa`) and the package's `hierarchy_stability` both compute it.
They agree to every printed digit (`check_package_indicator.jl`). `(q, κ)` and `(−q, −κ)`
are similar under a sign change on the odd tiers, and `−κ` gives the complex conjugate,
so `|q|` and `κ ≥ 0` cover the grid.

Against the dense generator on 20×20 grids (`symbol_check.log`), the frozen symbol
is the larger of the two in all 24 cases. The ratio of dense to frozen is 0.53–0.95,
and the frozen maximum always sits at `|q| = L`. Power iteration with the
matrix-free operator matches the dense abscissa to 1% (`growth_check.log`). It is
used for every production-grid rate below.

Smallest unstable `|q|` (the radius) on ±8/64 at depth 2 / 4 / 6 (`symbol_map.log`):

| bath | radius | largest frozen rate |
|---|---|---|
| cold strong (λ=0.8, γ=0.5, kT=0.1, Padé 2) | 0.75 / 1.0 / 1.25 | 3.4 / 6.1 / 8.2 |
| warm (λ=0.2, γ=kT=1, Padé 1) | 2.5 / 3.75 / 4.5 | 2.0 / 2.9 / 3.1 |
| cold weak (λ=0.05, γ=0.5, kT=0.1, Padé 2) | 2.0 / 3.25 / 4.25 | 0.8 / 1.2 / 1.3 |

The unstable wavenumbers are moderate (`κ ≈ 0.4–10`). The bath diffusion `D∂p²` damps
the highest ones, so the band overlaps the physical state's own momentum content.

## Characterisation

Growth rate of the complete generator by power iteration (`characterise.jl`,
`characterise2.jl`). Unless stated: harmonic `m = ω = ħ = 1`, Spectral, scaled, ±8 with
64 points.

| sweep | values | growth rate |
|---|---|---|
| depth 0 / 2 / 4 / 6, cold strong | | stable / 2.55 / 4.85 / 6.70 |
| depth 2 / 4 / 6, warm | | 0.55 / 1.19 / 1.34 |
| depth 2 / 4 / 6, cold weak | | 0.20 / 0.50 / 0.59 |
| box half-width at dq = dp = 0.25, cold strong depth 4 | L = 3 / 4 / 6 / 8 / 10 | 1.09 / 1.99 / 3.51 / 4.85 / 6.12 |
| box half-width at dq = dp = 0.25, warm depth 4 | L = 4 / 6 / 8 / 10 / 12 | stable / stable / 1.19 / 2.65 / 3.98 |
| momentum points at L = 8, cold strong depth 4 | n = 32 / 48 / 64 / 96 | 5.22 / 4.88 / 4.85 / 4.66 |
| discretisation, cold strong depth 4 | Spectral / FD2 / FD4 / FD8 | 4.85 / 5.11 / 5.42 / 5.46 |
| discretisation, warm depth 4 | Spectral / FD2 / FD4 / FD8 | 1.19 / 2.18 / 2.25 / 2.16 |
| coupling λ at γ = 0.5, kT = 0.1, Padé 2, depth 4 | 0.05 / 0.1 / 0.2 / 0.4 / 0.8 / 1.6 | 0.50 / 0.82 / 1.55 / 2.97 / 4.85 / 7.35 |
| coupling λ at γ = 0.5, kT = 1, Padé 2, depth 4 | 0.05 / 0.2 / 0.8 | stable / 0.60 / 3.36 |
| scaled / unscaled, cold strong depth 2 and 6 | | 2.55 / 2.52 and 6.70 / 6.71 |
| decomposition, cold strong depth 4 | Padé 1 / 2 / 4, Matsubara 2 | 4.03 / 4.85 / 5.36 / 4.26 |
| Padé 2 without terminator, cold strong depth 4 | | 9.29 |
| Brownian (λ=0.2, Ω=1, γ=0.5, kT=1), depth 4 | | 1.66 |

Blow-up times of a plain `heom_problem` (10³ × the initial maximum) are in
[`../nz_terminator/package_blowup.log`](../nz_terminator/package_blowup.log). For
cold strong on ±8/64 they are t = 12.5 / 7.6 / 5.9 at depth 2 / 4 / 6.

`seeding.log` tracks the depth-2 moment error against the exact Gaussian moments.
These moments are exact at depth 2 for a harmonic potential, so any error comes from
the growing modes. The error grows from t ≈ 1 (cold strong) or t ≈ 3 (warm) at the
generator's rate. Extrapolated back to t = 0, the unstable component starts near 1e-7
(cold strong) or 1e-9 (warm). Both are far above rounding, so the physical state itself
feeds the unstable modes.

## Remedies

Implementations are in `remedies.jl`, plus `:markov_full` from
`../nz_terminator/closure.jl`:

- **absorb**: `dWₙ/dt −= σ(q)Wₙ` on every auxiliary, never the root. σ rises with a C²
  smoothstep from 0 at `|q| = q₀` to σ₀ = 20 at the box edge.
- **taper**: the upward coordinate coupling uses `q̃(q)`. It equals `q` for `|q| ≤ q₀`
  and falls smoothly to 0 at the edge.
- **markov_full**: the scalar Markov closure of the final tier using every parent.
  At a frozen point and depth 1 it is exact, because the stationary auxiliaries form a
  coherent state.

Frozen-symbol bound (`remedy_symbol.log`): with σ₀ ≥ 5–10, the absorber removes all
growth beyond q₀. What remains is the hard cutoff's own rate at `|q| ≈ q₀`. An edge
remedy can therefore succeed only if `q₀` lies inside the radius above.

Growth rate of the complete generator, ±8/64 Spectral (`remedy_growth.jl`, `rg_*.log`):

| case | hard | absorb q₀=4 | absorb q₀=3 | taper q₀=4 | markov_full | markov_full + absorb q₀=4 |
|---|---|---|---|---|---|---|
| cold strong d2 | 2.55 | 0.81 | 0.56 | 0.87 | 0.52 | stable |
| cold strong d4 | 4.85 | 1.40 | 1.00 | 1.52 | 3.25 | 0.65 |
| warm d2 | 0.55 | stable | 0.002 | 0.005 | stable | stable |
| warm d4 | 1.19 | stable | stable | stable | stable | stable |
| warm d6 | 1.34 | stable | stable | stable | 0.30 | stable |
| cold weak d2 | 0.20 | stable | stable | stable | stable | stable |
| cold weak d4 | 0.50 | stable | stable | stable | stable | — |
| cold weak d6 | 0.59 | stable | stable | stable | 0.07 | stable |

"Stable" means the power iteration converged to the zero steady-state eigenvalue,
with a rate between −0.08 and 0. The entries 0.002 and 0.005 are not distinguishable
from it in a T = 30 run.

Long-time root-moment error against the exact Gaussian, to t = 60 (`accuracy.jl`,
`accuracy_closure.jl`). The hard cutoff and `markov_full` leave the depth-2 moments
exact in exact arithmetic, so their error is instability. The absorber's error also
includes its perturbation of the auxiliaries.

| case | method | t ≤ 10 | t ≤ 30 | t ≤ 60 |
|---|---|---|---|---|
| warm d2 | hard cutoff | 1.7e-3 | 5.0 | diverges at t = 40 |
| | markov_full | 7.3e-5 | 1.5e-4 | 2.0e-4 |
| | absorb q₀ = 5 / 4 / 3 | 2.7e-3 / 7.0e-3 / 2.5e-2 | 8.3e-3 / 1.9e-2 / 5.5e-2 | 1.2e-2 / 2.1e-2 / 5.6e-2 |
| warm d4 | hard cutoff | 1.7e-4 | 0.27 | diverges at t = 44 |
| | markov_full | 1.2e-6 | 2.6e-6 | 7.1e-6 |
| | absorb q₀ = 5 / 4 | 4.1e-4 / 1.1e-3 | 6.0e-4 / 2.2e-3 | 7.8e-4 / 2.2e-3 |
| cold weak d2 | hard cutoff | 8.9e-12 | 2.1e-10 | 2.2e-9 |
| | markov_full | 1.7e-11 | 1.7e-11 | 1.9e-11 |
| | absorb q₀ = 5 / 4 | 1.4e-8 / 1.7e-6 | 1.2e-7 / 3.5e-6 | 1.2e-7 / 3.5e-6 |
| cold weak d4 | hard cutoff / absorb q₀ = 5 | 2.0e-12 / 1.8e-9 | 2.1e-12 / 2.3e-9 | 1.0e-10 / 2.3e-9 |
| cold strong d2 | hard cutoff | 2.0 | diverges at t = 12 | |
| | markov_full | 1.2e-3 | 6.6e-3 | diverges at t = 52 |
| | markov_full + absorb q₀ = 5 | 1.1e-3 | 9.7e-3 | 1.7e-2 |

### Transformed hierarchy of Gatto et al.

`frozen_forms.jl` builds the frozen-point hierarchy in density-matrix variables
`x = q + y/2`, `x' = q − y/2` (ħ = 1) in three forms:

- the package's K-mode form;
- the left/right 2K-mode standard form, Eq. (19);
- the Gaussian and Bogoliubov transformed form, Eq. (60) of
  [`../../references/stabilized_heom_implementation.md`](../../references/stabilized_heom_implementation.md).

For real rates the sine sectors are invariant and start at zero, so 2K transformed
modes remain. All three converge to the exact frozen-path influence functional
(`frozen_validate.log`), but the transformed truncation converges much more slowly.
For the cold, strong Padé 1 bath at x = 1, x' = 0.5, depth 10 gives a relative error
of 1.3e-12 in the package form and 1.5e-4 in the transformed form. At x = 2, x' = −1
the errors are 2.4e-4 and 8.5.

Largest frozen real part over |q| ≤ 8, κ ≤ 12 (`gatto_symbol.log`):

| bath, depth | package (= standard) | markov_full | Gatto transformed |
|---|---|---|---|
| warm, 2 / 4 / 6 | 2.04 / 2.87 / 3.14 | 0.13 / 1.48 / 2.02 | 2.90 / 3.96 / 4.61 |
| cold weak, 2 / 4 / 6 | 0.81 / 1.17 / 1.30 | stable / 0.50 / 0.75 | 0.91 / 1.32 / — |
| cold strong, 2 / 4 | 3.43 / 6.09 | 1.25 / 4.34 | 10.4 / 16.3 |

The transformation balances the creation-like coupling in the auxiliary space, which
Gatto et al. show suppresses non-normal transient amplification for bounded system
operators. Here the coupling grows without bound with `|q|` and `y`. The transformed
truncation has the larger frozen growth rate and converges more slowly. It would also
double the number of modes, or increase the member count about 11× for K = 3 at depth 6.

### Not tested

A position-dependent displacement of the generating function (the Drude-only
transformation of [`../../references/wigner_heom_audit.md`](../../references/wigner_heom_audit.md),
Tanimura's Wigner form) removes the `q` coupling. It replaces it with an upward
coupling `(2λ/m)p`, which is unbounded at the momentum edge. On a square box with
`m = 1`, it is `1/γ` times as strong there as `2imag(c)q/ħ` is at the q edge (2× for
γ = 0.5). The frozen-symbol argument
predicts the same instability at the p edge, so it was not pursued.

## Recommendations

1. Keep the hard cutoff as the default, with the new diagnostic and warnings.
2. Land `:markov_full` (#37) as an option. It is cheap (1.2–1.6× a right-hand side)
   and more accurate per tier. It is also the most effective stabiliser found for warm
   and weak baths.
3. If an absorber is added, keep it opt-in, and put `q₀` well outside the state and
   its auxiliaries (`q₀ ≥ 5` on ±8 here). Report the perturbation it introduces.
4. Cold, strong coupling over long times needs a representation in which the depth
   required does not grow with `|q|`. None of the remedies tested provides one.

## Files

| File | Purpose |
|---|---|
| `common.jl` | Frozen symbol, dense generator, power-iteration growth rate, blow-up timer |
| `remedies.jl` | Absorber and taper wrappers (matrix-free, dense, frozen) |
| `symbol_check.jl`, `growth_check.jl` | Frozen symbol and power iteration against dense eigenvalues |
| `symbol_map.jl` | Unstable region in `(q, κ)` |
| `characterise.jl`, `characterise2.jl` | Growth-rate sweeps (`char_*.log`, `char2_*.log`) |
| `seeding.jl` | Moment error against time: what seeds the modes |
| `remedy_symbol.jl`, `remedy_growth.jl` | Remedies: frozen bound and full generator |
| `accuracy.jl`, `accuracy_closure.jl` | Long-time moment accuracy with remedies |
| `frozen_forms.jl`, `frozen_validate.jl`, `gatto_symbol.jl` | Package, standard and transformed frozen hierarchies |
| `check_package_indicator.jl` | `hierarchy_stability` against the prototype |
| `regression_candidate.jl` | Cost and blow-up time of candidate regression tests |
