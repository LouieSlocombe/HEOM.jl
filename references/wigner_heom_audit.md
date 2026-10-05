# Wigner–Moyal HEOM audit

Audit date: 2026-10-05. Scope: Hamiltonian Wigner evolution, real-FFT and central-difference discretisations, exponential/Drude baths, hierarchy couplings and cutoff, Caldeira–Leggett (CL) evolution, and the analytical states used for validation. Tests run with Julia 1.13.1. The results below are implementation checks, not reproductions of published figures.

Final verification: **1,726/1,726 assertions passed**, **656/656 executable source
lines covered**, JuliaFormatter check passed, and `git diff --check` passed.

## Literature and convention checks

The core Hamiltonian and hierarchy signs were consistent; they did not require reversal or additional factors of ħ.

| Item | Convention checked | Finding |
|---|---|---|
| Wigner generator | `−(p/m)∂q + V′∂p − ħ²V‴∂p³/24 + ħ⁴V⁽⁵⁾∂p⁵/1920 + …` | Correct, including mass and non-unit ħ |
| Exact potential symbol | `i[V(q+ħκ/2)−V(q−ħκ/2)]/ħ`, with `∂p → iκ` | Correct Fourier sign and normalisation |
| Periodic derivatives | Odd derivative has zero unpaired Nyquist mode; second derivative retains it | Correct; unitary generator is skew-symmetric |
| CL damping | `friction*∂p(pW) + mass*friction*kT*∂p²W` | Conservative flux and thermal coefficient correct |
| Hierarchy | Upward `∂p`; downward `nₖ[Re(cₖ)∂p + 2Im(cₖ)q/ħ]` | Correct for the package's direct density-HEOM auxiliaries |
| Drude counterterm | Add `λq²` once to the supplied physical potential | Correct; zero initial auxiliaries imply a factorized bare-bath preparation and initial slip |
| Truncations | Missing upper tiers set to zero; omitted Matsubara poles replaced by diffusion | Distinct approximations; each needs convergence checks |

The Wigner symbol and CL coefficient comparison use [Cabrera et al., *Physical Review A* 92, 042122 (2015)](https://arxiv.org/pdf/1212.3406), especially the Moyal and open-system equations. Their damping parameter is half the package's `friction`. This conversion accounts for their `2γ` drift and `2mγkT` diffusion.

[Tanimura, *J. Chem. Phys.* 153, 020901 (2020)](https://arxiv.org/pdf/2006.05501) supplies the density-HEOM framework. Applying `[q,ρ] → iħ∂pW` and `{q,ρ} → 2qW` gives the package's real auxiliary equations. Some Wigner-space presentations transform the auxiliaries so that their downward operator contains momentum/friction directly; individual auxiliary coefficients must not be compared across these representations without that transformation.

[Tanimura, *J. Chem. Phys.* 142, 144110 (2015)](https://arxiv.org/pdf/1502.04077), Sections III–V, provides Wigner HEOM and harmonic Brownian-oscillator benchmarks. Its correlated reduced equilibrium generally differs from the isolated oscillator Gibbs state. Neither the CL thermal Gaussian nor an isolated Gibbs state is therefore a universal quantum HEOM equilibrium target.

For example, the Drude-only auxiliary transformation derived in this audit is
`W̃ₙ = Σⱼ₌₀ⁿ binomial(n,j)(2λq)ⁿ⁻ʲ Wⱼ`. It leaves the root unchanged, absorbs the
explicit counterterm, and changes the downward operator to
`n[(2λ/m)p + Re(c₀)∂p]`. Thus the equivalent friction convention is
`ζ=2λ/(mγ)`. Zero raw auxiliaries transform to `W̃ₙ(0)=(2λq)ⁿW₀`, which explains
why setting both representations' auxiliaries to zero describes different preparations.

## Defects repaired

1. **Drude remainder cancellation.** The old expression subtracted retained Matsubara residues from `1/x−cot(x)`, where `x=ħγ/(2kT)`. Near a retained pole, large cancelling contributions corrupted a small, smooth remainder. Tests reproduce order-one relative error with many retained poles. The implementation now evaluates the positive omitted sum

   ```math
   D=\lambda\hbar\sum_{k=K+1}^{\infty}\frac{2x}{\pi^2k^2-x^2}
   ```

   using a short direct sum and Euler–Maclaurin tail. For `λ=0.3, ħ=kT=1,
   γ=2π(1+1e-6), K=10000`, the old result was `9.34821e-6`; the corrected value is
   `1.9097657434490765e-5`. Accurate argument reduction protects the distance to a
   nearby omitted pole. Independent 512-bit cot-identity regressions cover small
   `x`, nearby poles, half-integer arguments, and `K=10000`. Coincident poles remain
   unsupported.

2. **Invalid Hamiltonians could be accepted.** Complex potentials are incompatible with this real Wigner commutator implementation; constant imaginary offsets could even disappear silently. Infinite constant potentials could disappear under differentiation. Potential evaluations now require finite real values. Mass and ħ must remain finite and positive after Float64 conversion, and constructed kinetic/finite-difference coefficients must be finite.

3. **State geometry could be silently changed.** Direct RHS calls could broadcast a single column or reinterpret a differently shaped matrix with the same number of entries. Wigner and CL RHS calls now enforce both input and output grid shapes. Their problem constructors reject non-real/nonfinite initial states, including Float64 overflow.

4. **Unrepresentable grids and analytical parameters.** Finite endpoints could still produce infinite/zero spacing or duplicate Float64 points. These grids are rejected. Harmonic-state parameters are validated, and `harmonic_evolution(...; omega=0)` now returns the continuous free-streaming limit instead of NaNs.

## Independent analytical comparisons

The new [analytical benchmark tests](../test/analytic_benchmarks.jl) supplement the existing harmonic coherent/Fock/cat, quartic, sinusoidal, CL Gaussian and hierarchy-moment tests. All quoted errors are maximum absolute grid errors from this audit, not accuracy guarantees for other parameters.

| Reference | Measured result | What it tests |
|---|---:|---|
| Free Gaussian spreading | `9.80e-13` | Streaming sign, mass and correlated covariance |
| Constant-force Gaussian | `8.92e-13` | Force sign and acceleration |
| Exact sextic ground state | Stationary RHS residual `1.91e-13` | Non-Gaussian stationary state and ħ⁴ correction |
| Full Drude Gaussian, depth 2 | `8.87e-4` | Full distribution beyond closed low moments |
| Same bath, depth 4 | `1.33e-5` | Hierarchy convergence |
| Same bath, depth 6 | `1.77e-7` | Hierarchy convergence |

For the sextic case, choose `ψ(q)=N exp(−a q⁴−b q²)` and derive the potential directly from the Schrödinger equation:

```math
V(q)=\frac{\hbar^2}{2m}
\left[16a^2q^6+16abq^4+(4b^2-12a)q^2\right],\qquad E=\frac{\hbar^2b}{m}.
```

This is a member of the exactly solvable family in [Lévai and Ishkhanyan, *Modern Physics Letters A* 34, 1950134 (2019), Eqs. (2)–(3)](https://arxiv.org/pdf/1904.09488). The test uses `m=1.3, ħ=0.8, a=0.04, b=0.45` and direct cosine quadrature of `ψ(q+y/2)ψ(q−y/2)`, independently of the solver's FFT and automatic derivatives. Classical-only and ħ²-only residuals are respectively `5.15e-2` and `7.84e-3`; the test is sensitive to the missing quantum terms.

For the Drude case, an independent linear Langevin representation generates each positive exponential noise component. Matrix exponentials give the mean and covariance, hence the full Gaussian density. This implements the Gaussian propagation structure described by [Fleming, Roura and Hu, *Annals of Physics* 326, 1207 (2011), Eqs. (35)–(36)](https://arxiv.org/pdf/1004.1603). Parameters are `m=1.3, ω=0.85, ħ=0.9, λ=0.22, γ=1.1, kT=0.75, K=1`; comparisons extend to `t=1.2`. The reference includes initial slip and is exact for the chosen exponential bath plus white-noise remainder. It does not establish infinite-Matsubara convergence. First and second moments agree already at depth two, despite its visibly larger full-state error.

The [CL tests](../test/caldeira_leggett.jl) also propagate an initially negative `n=1` Fock Wigner function using an independent Ornstein–Uhlenbeck characteristic-function solution. They compare the entire distribution, covariance, purity and central negativity. [Roy and Venugopalan, *Exact Solutions of the Caldeira–Leggett Master Equation*](https://arxiv.org/pdf/quant-ph/9910004), Eq. (20), gives the relevant factorization into transported initial data and Gaussian noise. Their `γ` equals `friction/2`, and their symbol `D` is four times the package's momentum diffusion coefficient.

Separately, [bath tests](../test/heom.jl) reconstruct the symmetric spectrum and compare it with `ħJ(ω)coth(ħω/(2kT))`. The antisymmetric part reproduces `ħJ(ω)`. These checks include negative low-temperature residues, and verify decreasing finite-frequency error for `K=2,8,32`. They test the physical bath definition independently of the pole formulas.

An additional classical Drude equilibrium check uses the correlated hierarchy
`Wₙ=(-2λq)ⁿWβ`, with the normalized Boltzmann root. Its RHS vanishes in every
nonterminal tier; the nonzero top-tier residual explicitly exposes the hard cutoff.

## Remaining numerical limits

- Converge hierarchy depth, Matsubara count, box extent, mesh resolution and ODE tolerances independently. A trace or moment check alone cannot certify the full distribution.
- The first omitted Matsubara rate exceeding the cutoff makes the omitted tail positive; the white-noise approximation also needs that rate to be fast relative to the system dynamics. The README and bath docstring now make this distinction explicit.
- Hard hierarchy truncation and unscaled auxiliaries remain the implemented method. Very deep/strong-coupling cases can require scaling or a more advanced closure; this audit does not claim convergence throughout that regime.
- The grid is periodic. Edge weight and spectral tails must remain small; no absorbing boundary, positivity projection or automatic renormalisation is applied.
- General low-temperature quantum equilibrium, response functions from correlated equilibrium, and literal published-figure reproduction remain additional benchmarks. The new finite-bath Gaussian test should not be presented as those validations.

Reproduce the complete numerical/package checks with `julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'`, then run `julia --project=build_tools build_tools/coverage.jl` and `julia --project=build_tools build_tools/format.jl --check`. Install the build-tool environment as described in the README if needed.
