# Final-tier closures for the Wigner HEOM: prototype results

Exploratory code, not part of the package. It asks whether the Nakajima–Zwanzig (NZ)
terminator of [Fay, JCP 157, 054108 (2022)](https://arxiv.org/abs/2205.09270), Eq. (14),
should replace the hard depth cutoff in `heom_operator`.

## Conclusions

1. **Fay's terminator is not worth adding.** At a fixed depth it reduces the error by
   at most 2.4×, and in the cold, strong cases often not at all. One extra depth level
   reduces it by 2–18× instead. Each evaluation costs 9–64× a hard-cutoff evaluation,
   because every eliminated child needs a sparse LU solve. The diagonal form published
   as Eq. (14) is worse than the hard cutoff at every depth ≥ 1 in the cold cases,
   by up to four orders of magnitude for the quartic potential.
2. **It destabilises the hierarchy.** It always raises the largest real eigenvalue
   slightly. In time integration it blows up where the hard cutoff stays bounded:
   - Warm bath, depths 1–4: blows up at t ≈ 32–38; the hard cutoff is bounded to t = 40.
   - Cold, strong bath at depth 0: blows up at t ≈ 16; the hard cutoff is stable.
   - Depth 0 in both warm and cold, weak baths: negative variances by t = 10.

   This is the minimal-hierarchy regime where Fay reports his main gains.
3. **The free scalar closure with every parent (`:markov_full`) is the only useful variant.**
   - Accuracy: 1.2–3.5× lower error at depth ≥ 2 for about 1.2–1.6× the cost of a
     right-hand-side evaluation. That is less than one extra tier gains at weak
     coupling (4–18×). In the strong-coupling cases, where a tier gains only 2–3×, it
     is comparable to adding one. It is worse than the hard cutoff at depths 0–1 in
     cold transients.
   - Equilibrium: in the warm case at depth 1 it relaxes to the correct moments,
     ⟨δq²⟩ = 1.079 against an exact 1.082. The hard cutoff gives about 0.76.
   - Stability: it removes the depth-1 instability exactly, and it is the only closure
     that delays the growth of the box-edge modes rather than accelerating it.
4. **Separate finding about the package itself:** the hard cutoff has growing
   modes localised at the box edge in every regime tested, with both discretisations.
   In the cold, strong regime used by `test/cold_strong_benchmarks.jl`, a plain
   `heom_problem` with the default `Spectral()` discretisation blows up at
   t = 12.5 / 7.6 / 5.9 for depth 2 / 4 / 6 (box ±8, 64 points). The existing tests
   stop at t = 0.8.

## Files

| File | Purpose |
|---|---|
| `closure.jl` | Five final-tier closures, a matrix-free right-hand side, and a dense generator |
| `validate.jl` | Consistency checks (all pass) |
| `benchmark.jl`, `summarize.jl` | Accuracy and cost against depth; condenses the logs into tables |
| `stability.jl` | Largest real eigenvalue of each closed generator (20×20 grid) |
| `blowup.jl`, `package_blowup.jl` | Time to exceed 1000× the initial maximum, prototype and package |
| `relaxation.jl` | Long-time second moments at depths 0–2 |
| `gaussian_reference.jl` | Exact Gaussian solution for a finite bath (copied from the tests) |

Run any script with `julia --project=prototypes/nz_terminator prototypes/nz_terminator/<script>`.

## Closures

A hard cutoff sets tier `d+1` to zero. Each closure instead eliminates tier `d+1`
adiabatically. A child `m` satisfies `Xₘ = R(γₘ) Σⱼ mⱼ Φⱼ W_{m−eⱼ}` and returns
`wₖ ∂p Xₘ` to each parent `m−eₖ`, where `Φⱼ = real(cⱼ)∂p + 2imag(cⱼ)q/ħ`.

| Kind | `R(z)` | Parents feeding each child |
|---|---|---|
| `:none` | 0 | – (the package today) |
| `:markov_diag` | `1/z` | own parent only |
| `:markov_full` | `1/z` | all parents |
| `:nz_diag` | `(z − L_WM − D∂p²)⁻¹` | own parent only (Fay Eq. 14) |
| `:nz_full` | same | all parents (Fay's full kernel, Eq. 10, the time-derivative truncation) |

`validate.jl` confirms the following:

- `:nz_full` equals exact elimination of tier `d+1`. Embedding the solved children
  in a depth-`d+1` hierarchy reproduces the closed derivative, and each child is
  stationary to 1e-12 relative, in both scaled and unscaled conventions.
- The dense generator matches the matrix-free one.
- The diagonal and full variants coincide for a single-mode bath.

Two implementation notes:

- With a matrix-free right-hand side, the full closure needs *fewer* solves than the
  diagonal one: one per child, rather than one per (parent, mode) pair. Fay's cost
  argument for the diagonal form applies to precomputed dense blocks, not here.
- The bath diffusion `D∂p²` is part of each child's own equation, so it belongs
  inside the resolvent. A depth closure does not double-count it. Only Fay's
  separate tail correction (Eqs. 12 and 24) would replace it.

## Benchmark setup

- Oscillator: m = ω = ħ = 1, starting from a Gaussian with mean (1, 0) and covariance 0.5·I.
- Grid: `FiniteDifference(4)` with 48×48 points on [−6, 6]², scaled hierarchy, Vern7 at tolerance 1e-10.
- Error: the maximum over saved times of `max|W − W_ref| / max|W₀|`.
- Reference: a deep hard cutoff on the same grid. This isolates depth error, the only
  error a terminator targets. For a harmonic oscillator second moments are already
  exact at depth 2, so the full Wigner function is compared instead.
- Gain: given in parentheses, relative to the hard cutoff at the same depth (values
  below 1× mean the closure is worse).

### Warm (λ = 0.2, γ = 1, kT = 1), Padé 1, t ≤ 10
Reference: depth 12, converged to 2e-7.

| depth | none | markov_full | nz_diag | nz_full | none at depth+1 |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.53 | 0.38 (1.4×) | 0.22 (2.4×) | 0.22 (2.4×) | 0.091 |
| 1 | 0.091 | 0.05 (1.8×) | 0.052 (1.7×) | 0.051 (1.8×) | 0.016 |
| 2 | 0.016 | 0.013 (1.2×) | 0.0082 (1.9×) | 0.0079 (2.0×) | 0.0039 |
| 3 | 0.0039 | 0.0026 (1.5×) | 0.0019 (2.0×) | 0.0017 (2.2×) | 0.001 |
| 4 | 0.001 | 0.00051 (2.0×) | 0.00048 (2.1×) | 0.00045 (2.2×) | 0.00022 |
| 5 | 0.00022 | 0.00014 (1.6×) | 0.00011 (2.1×) | 0.0001 (2.2×) | 6.0e-5 |
| 6 | 6.0e-5 | 2.8e-5 (2.1×) | 2.9e-5 (2.1×) | 2.7e-5 (2.2×) | — |

Cost per right-hand-side evaluation, relative to none: markov_full 1.2×, nz_diag 12×, nz_full 9×.

### Cold, weak (λ = 0.05, γ = 0.5, kT = 0.1), Padé 2, t ≤ 10 (Fay's regime)
Reference: depth 10, converged to 5e-10.

| depth | none | markov_full | nz_diag | nz_full | none at depth+1 |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.18 | 0.38 (0.46×) | 0.13 (1.3×) | 0.13 (1.3×) | 0.01 |
| 1 | 0.01 | 0.013 (0.8×) | 0.012 (0.89×) | 0.0067 (1.5×) | 0.00081 |
| 2 | 0.00081 | 0.00032 (2.6×) | 0.0016 (0.51×) | 0.00047 (1.7×) | 4.8e-5 |
| 3 | 4.8e-5 | 2.0e-5 (2.4×) | 9.1e-5 (0.52×) | 2.7e-5 (1.8×) | 3.5e-6 |
| 4 | 3.5e-6 | 1.0e-6 (3.4×) | 6.5e-6 (0.54×) | 1.6e-6 (2.1×) | 2.4e-7 |
| 5 | 2.4e-7 | 6.8e-8 (3.5×) | 4.3e-7 (0.54×) | 1.1e-7 (2.1×) | — |

Cost relative to none: markov_full 1.6×, nz_diag 31×, nz_full 19×.

### Cold, strong (λ = 0.8, γ = 0.5, kT = 0.1), Padé 2, t ≤ 2.5
Reference: depth 10, converged to 2e-4. The deep hierarchy blows up before t = 10.

| depth | none | markov_full | nz_diag | nz_full | none at depth+1 |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.21 | 0.67 (0.31×) | 0.28 (0.75×) | 0.28 (0.75×) | 0.09 |
| 1 | 0.09 | 0.16 (0.57×) | 0.17 (0.54×) | 0.088 (1.0×) | 0.043 |
| 2 | 0.043 | 0.023 (1.9×) | 0.12 (0.36×) | 0.037 (1.1×) | 0.017 |
| 3 | 0.017 | 0.011 (1.5×) | 0.051 (0.33×) | 0.012 (1.4×) | 0.0068 |
| 4 | 0.0068 | 0.0034 (2.0×) | 0.057 (0.12×) | 0.0055 (1.2×) | 0.0029 |
| 5 | 0.0029 | 0.00093 (3.2×) | 0.18 (0.017×) | 0.0022 (1.3×) | 0.0011 |
| 6 | 0.0011 | 0.00042 (2.7×) | 0.32 (0.0036×) | 0.00086 (1.3×) | — |

Cost relative to none: markov_full 1.3×, nz_diag 38×, nz_full 17×.

### Quartic q²/2 + 0.08q⁴, cold, strong, Padé 2, t ≤ 2.5
Reference: depth 10, converged only to 1.2e-3, so rows at depth ≥ 4 are at the
reference's precision.

| depth | none | markov_full | nz_diag | nz_full | none at depth+1 |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.17 | 0.7 (0.25×) | 0.21 (0.81×) | 0.21 (0.81×) | 0.066 |
| 1 | 0.066 | 0.14 (0.47×) | 0.13 (0.51×) | 0.069 (0.95×) | 0.028 |
| 2 | 0.028 | 0.019 (1.5×) | 0.094 (0.29×) | 0.028 (0.99×) | 0.0079 |
| 3 | 0.0079 | 0.0066 (1.2×) | 0.59 (0.013×) | 0.0091 (0.87×) | 0.0028 |

Cost relative to none: markov_full 1.3×, nz_diag 64×, nz_full 40×.

The `markov_diag` closure is omitted from these tables. It is no better than
`markov_full` in the warm case, worse in the cold weak case, and blows up in the
cold strong cases.

### Resolvent cost

For L_WM + D∂p² with a quartic potential and 2 Moyal terms:

| Grid, FD order | 32², 4 | 48², 4 | 64², 4 | 128², 4 | 48², 8 | 128², 8 |
|---|---:|---:|---:|---:|---:|---:|
| One LU solve / one L_WM application | 15× | 27× | 35× | 42× | 73× | 100× |

Each eliminated child costs as much as roughly 15–100 extra hierarchy members, and
there is one child per tier-`d+1` member. The spectral discretisation would need
iterative solves instead, which this prototype does not attempt.

## Stability

### Eigenvalues
`stability.jl` reports the largest real eigenvalue of each closed generator on a 20×20
grid. For the cold, strong Padé-1 bath:

| depth | none | markov_full | nz_full |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0.47 |
| 1 | 0.89 | 0 | 1.19 |
| 2 | 1.69 | 0.39 | 1.98 |
| 3 | 2.41 | 1.19 | 2.70 |
| 5 | 3.65 | 2.57 | 3.90 |

The very strong, Matsubara and quartic cases show the same ordering. In every case
`nz_*` is above `none`, and `markov_full` is lowest from depth 1.

Where the unstable modes live:

- 84–100% of each unstable eigenvector's weight sits in the outer quarter of the box,
  and each mode oscillates fast.
- Growth increases with depth and with box half-width L (depth 2: 1.7 at L = 6,
  2.7 at L = 10).
- Every finite bath here has a positive noise spectrum, so the cause is not
  anti-diffusion in the bath. The upward coupling `2imag(c)q/ħ` grows with |q|.
- The default spectral discretisation gives the same picture. Even the warm and
  cold, weak baths have positive growth rates there (0.2–1.8, depending on box and depth).

A plausible but untested reason why the resolvent makes things worse: over its memory
time 1/ν, the free flow carries edge regions across the periodic boundary, where `q`
changes sign.

### Time to exceed 1000× the initial maximum
`blowup.jl`, FD 48², harmonic, integrated to t = 40. "Inf" means it stayed bounded.

| Bath | depth | none | markov_full | nz_diag | nz_full |
|---|---:|---:|---:|---:|---:|
| warm | 1–4 | Inf | Inf | 32–37 | 33–38 |
| cold weak | 0–4 | Inf | Inf | Inf | Inf |
| cold strong | 0 | Inf | Inf | 15.8 | 15.8 |
| cold strong | 1 | 23.1 | Inf | 7.6 | 12.6 |
| cold strong | 2 | 12.5 | Inf | 6.0 | 10.7 |
| cold strong | 3 | 9.5 | 16.0 | 4.4 | 8.7 |
| cold strong | 4 | 7.9 | 11.4 | 3.5 | 7.4 |

`markov_diag` blows up at t ≈ 1.4–3.9 in the cold, strong case.

The package's own default `Spectral()` hard cutoff (`package_blowup.jl`) gives the
following blow-up times, for depth 2 / 4 / 6:

- Cold, strong: 14.0 / 8.2 / 6.5 on ±6 with 48 points, and 12.5 / 7.6 / 5.9 on ±8 with 64 points.
- Warm, ±8 with 64 points: bounded / bounded / 36.6.
- Every other combination tested stays bounded to t = 40.

### Relaxation (`relaxation.jl`)
Variances ⟨δq²⟩, ⟨δp²⟩ at t = 10, 20, 30, 40.

Warm bath; the exact values are (1.082, 1.111):

| depth | none | markov_full | nz_full |
|---:|---|---|---|
| 0 | 0.51 → 0.78 (no friction, oscillates) | heats toward a uniform box | negative variances by t = 10 |
| 1 | 0.76, 0.79 at t = 10–30, then corrupted | **1.079, 1.107**, stable | 1.04 at t = 20, then diverges |
| 2 | correct until t = 20, then −4356 by t = 40 | **1.081, 1.110**, stable | correct until t = 20, then blows up |

In the cold, weak case, depths 1–2 agree across closures to within the slow
relaxation still under way by t = 40. Depth-0 NZ again gives negative variances.

## Limitations

- Only Drude baths (Padé or Matsubara) without coupled modes, so no underdamped
  Brownian baths. All results use one grid family, finite differences for every
  closure, and a harmonic or one quartic potential.
- Timings were measured with other jobs running, so cost ratios are approximate. The
  ordering between closures is robust.
- I first tried to compute stationary states as dense null vectors. The growing
  box-edge modes contaminate them, so I abandoned that approach in favour of time
  integration.
- Fay's separate tail correction (Eqs. 12 and 24), which would replace the `D∂p²`
  terminator, was not implemented. It would need resolvent solves on every
  hierarchy member, not only at the final tier.
