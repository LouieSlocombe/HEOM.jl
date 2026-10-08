# Low-temperature and strong-coupling support

This extends the [Wigner–HEOM audit](wigner_heom_audit.md). The implementation supports
finite positive temperatures and strong linear coupling to a Gaussian bath for harmonic
and anharmonic system potentials. It does not use a weak-coupling, rotating-wave or
high-temperature approximation to the system dynamics. Finite bath expansions,
hierarchy depth, periodic grids and ODE tolerances remain numerical approximations.
Exactly zero temperature is a separate correlation-decomposition problem; neither
thermal constructor accepts `kT=0`.

## Numerical capabilities

**Scaled auxiliaries.** `heom_operator(...; scaled=true)` propagates

```math
\widetilde W_{\mathbf n}=W_{\mathbf n}/\sqrt{\prod_k n_k! |c_k|^{n_k}}.
```

This is a diagonal similarity transformation of the finite hierarchy. Raising/lowering
factors are computed locally, and factorials and large powers are never formed.
Zero coefficients use a unit scale. This follows the conditioning strategy of
[Shi et al., JCP 130, 084105 (2009)](https://doi.org/10.1063/1.3077918); a detailed
scaled formulation is also given by [Xu, Xu and Yan](https://arxiv.org/pdf/0905.1432).
It preserves the physical root and accepts negative or complex coefficients.
`rescale_hierarchy(U, old_op; scaled=true)` converts a correlated initial hierarchy
or restart. This changes the scaling convention only, not the bath expansion or basis.
The default remains `scaled=false` for compatibility.

**Efficient thermal expansion.** `drude_lorentz_pade_bath(...; pade=N)` implements
the `[N/N]` Bose Padé expansion, including its consistent Drude residue and
white-noise remainder. Jacobi-matrix eigenvalues and spectral weights give the poles
and residues without unstable products of pole differences. The coefficient and
remainder formulas follow [Ding et al., JCP 135, 164107 (2011), Eqs. (5)–(8)](https://arxiv.org/pdf/1107.0249).
The Matsubara constructor remains available as an independent bath approximation.
Increase `pade` or `matsubara` independently of hierarchy `depth`.

**Repeated poles and general real bases.** `ExponentialBath` accepts `weights` and
sparse `mixing`, describing

```math
C(t)=w^T e^{-\Gamma t}c,\qquad\Gamma=\operatorname{diag}(\nu)+M.
```

Mixing adds the unscaled same-tier term
`−Σₖⱼ nₖ Mₖⱼ Wₙ₋ₑₖ₊ₑⱼ` and replaces the upward derivative by
`Σₖ wₖ ∂p Wₙ₊ₑₖ`. Diagonal damping and the downward real/imaginary coefficients
retain their previous conventions. Scaled couplings apply the corresponding
similarity factors. This is a finite real basis realization of the generalized
correlation framework in [Ikeda and Scholes, JCP 152, 204101 (2020)](https://arxiv.org/pdf/2003.06134).
Real rotation blocks also represent damped oscillatory correlations.

Both thermal constructors automatically combine near-colliding poles. At the exact
Drude–Matsubara coincidence, the combined contribution is

```math
(\lambda kT-i\lambda\hbar\gamma-2\lambda\gamma kT\,t)e^{-\gamma t}.
```

A Jordan block represents this expression directly, without shifting cutoff or
temperature and without divergent individual residues. Matsubara truncation must
retain the coincident pole. Mixing storage is sparse: a 10,000-pole construction
requires roughly 1 MB, rather than a dense 10,000-square matrix. Rate matrices must
have positive-real eigenvalues; physical consistency of an arbitrary user-supplied
correlation is still the caller's responsibility.

**Stiff integration.** Every `heom_problem` supplies exact Jacobian-vector products
and a zero time derivative. Spectral problems can use a Krylov stiff solver without
automatic differentiation through FFTW:

```julia
using HEOM, OrdinaryDiffEqRosenbrock, LinearSolve, ADTypes
alg = Rodas5P(autodiff=AutoFiniteDiff(),
             linsolve=KrylovJL_GMRES(), concrete_jac=false)
sol = solve(prob, alg; abstol=1e-9, reltol=1e-9)
```

Despite the fallback differentiation setting, the supplied JVP is exact. For
finite differences, `sparse(op)` assembles the complete constant generator, including
scaling and same-tier mixing. `heom_problem(W0,tspan,op; jacobian=:sparse)` supplies
it to implicit solvers such as `Rodas5P()`. Sparse factorization trades memory for
speed; the matrix-free path avoids storing the full Jacobian. Solver packages are
optional consumer dependencies and test dependencies, not runtime requirements for
loading HEOM. See the [SciML Rosenbrock interface](https://docs.sciml.ai/OrdinaryDiffEq/stable/semiimplicit/Rosenbrock/).

`hierarchy_size(modes, depth)` estimates the exact member count before allocating.
`max_ados` optionally rejects a calculation exceeding a user-selected budget.
One Float64 hierarchy array costs `8*nq*np*members` bytes; solver stages require
additional arrays. Scaling and Padé improve feasibility but do not eliminate
combinatorial hierarchy growth.

## Independent physical validation

The [cold-system tests](../test/cold_strong_benchmarks.jl) use `m=ω=ħ=1`, `kT=0.1`
and `λ=0.8`. With `γ=0.5`, the zero-frequency damping rate is `3.2ω`, so these are
low-temperature, strongly coupled tests.

The equilibrium reference integrates the continuum fluctuation–dissipation relation,
without a Padé, Matsubara or HEOM expansion:

```math
\chi(\Omega)=\left[m(\omega^2-\Omega^2)
  -\frac{2i\lambda\Omega}{\gamma-i\Omega}\right]^{-1},\qquad
\langle q^2\rangle=\frac{\hbar}{\pi}\int_0^\infty
\coth\left(\frac{\hbar\Omega}{2kT}\right)\operatorname{Im}\chi(\Omega)\,d\Omega.
```

Inserting `m²Ω²` gives the momentum variance. These are the coupled-system equilibrium
relations discussed by [Hänggi and Ingold, Chaos 15, 026105 (2005)](https://arxiv.org/pdf/quant-ph/0412052).
The resulting variances are `0.4013479887` and `0.7669323102`, rather than the
bare-Gibbs value `0.5000454` for both. Exact stationary solutions of the finite-bath
linear reference converge toward these values as bath order increases. This checks
the decomposition and coupled equilibrium; it is not a claim that every finite-depth
grid hierarchy has already relaxed to that equilibrium.

Actual HEOM trajectories are compared with independent Gaussian memory solutions
for negative Drude residues, negative Matsubara residues and exact repeated poles.
The reference uses matrix exponentials and a signed covariance representation;
negative quantum components are not treated as independently sampleable classical
noises. For the generalized basis it constructs symmetric `S` satisfying
`S*w=real(c)` and uses `Γ*S+S*Γ'` for the signed noise matrix. The scalar force
correlation is then exactly `w'*exp(-Γ*t)*real(c)`. This implements the Gaussian
propagation structure of [Fleming, Roura and Hu](https://arxiv.org/pdf/1004.1603).
Depths 2, 4 and 6 give convergent full Wigner distributions, with final error below
`2e-7` in all three cases; the checks also cover mean and covariance.

Additional tests verify scaled/unscaled equivalence, an independent bath-basis
transformation, correlation functions against arbitrary precision at repeated poles,
low-temperature quantum spectra, and stiff solves against tightly integrated
explicit references. Harmonic validation is supplemented by the full quartic
evolution in the example below.

## Anharmonic convergence example

Run `examples/low_temperature_strong_coupling.jl` in an environment containing HEOM
and OrdinaryDiffEqVerner. It uses `V(q)=q²/2+0.08q⁴`, the cold/strong parameters above,
and a displaced pure initial state evolved to `t=0.8`. Seven complete hierarchies
vary each numerical control independently. Measured final-state pointwise changes:

| Control | Maximum change |
|---|---:|
| Depth 4 → 6 | `1.20e-6` |
| Depth 6 → 8 | `9.39e-8` |
| Padé 2 → 3 | `3.19e-3` |
| Padé 3 → 4 | `9.30e-4` |
| Grid 40² → 80² at fixed box | `1.51e-4` |
| Box ±7 → ±8.4 at fixed spacing | `1.30e-4` |

The bath expansion dominates this modest-cost example. Tighter accuracy requires
more poles and repeating **all** controls around the refined calculation; these
numbers do not certify arbitrary potentials, longer times or lower temperatures.
Also check ODE tolerances, boundary weight, Fourier tails and the desired observables.
For correlated equilibrium or response calculations, retain and restart the entire
relaxed hierarchy; resetting auxiliaries to zero changes the preparation.

Hard depth truncation remains a convergence-controlled approximation. No finite
depth is promised stable for every parameter set, and neither a small moment error
nor norm conservation establishes convergence of the full distribution. The new
features provide the numerical paths needed to perform these calculations and
their convergence studies without a high-temperature or weak-coupling restriction.

## Balanced bath compression

`compress_bath(bath; modes)` treats the retained correlation as the impulse response
of a real linear system,

```math
C(t)=w^T e^{-\Gamma t}c,\qquad
\dot x=-\Gamma x+\begin{pmatrix}\operatorname{Re}c&\operatorname{Im}c\end{pmatrix}u,
\qquad y=w^Tx,
```

and retains its `modes` leading balanced states. Balanced truncation of a Padé
decomposition is the truncated Padé decomposition of
[Takahashi and Tanimura, JCP 158, 044115 (2023), Appendix B](https://doi.org/10.1063/5.0135725).
With Hankel singular values `σᵢ`, the H∞ bound `2Σ_{i>r}σᵢ` of
Glover, Int. J. Control 39, 1115 (1984) gives

```math
\max_\omega |S(\omega)-S_r(\omega)|\le 4\sqrt2\sum_{i>r}\sigma_i,
```

because `S(ω)` is twice the real part of a fixed linear combination of the two
transfer-function columns. Gramian square roots and the SVD avoid inverting either
Gramian. `method=:residualize` is the singular-perturbation approximation. It is
computed by truncating the reciprocal system `Γ⁻¹` in the same coordinates, so it needs
only the retained projections. The eliminated `∫₀^∞ C dt` becomes momentum diffusion
(real part) and a `q²` potential (imaginary part over `ħ`). Thus `S(0)` and the net
static force constant are exact.

Truncation at total hierarchy depth is invariant under any linear change of real
bath basis. In generating-function form, tier `N` is the space of degree-`N`
polynomials in the mode variables, and a linear substitution preserves degree.
The reduced rate matrix is therefore returned in real Schur form, without
changing the root at any depth. The tests confirm this with rotated bases.

`harmonic_covariance(bath; mass, omega)` solves the stationary Lyapunov equation of
the generalized Langevin system above with signed noise `S`, `S*w=real(c)`. It is the
exact infinite-depth equilibrium of a harmonic oscillator coupled to that
decomposition, the test advocated by
[Tokieda, Phys. Rev. Research 7, 043178 (2025)](https://doi.org/10.1103/bv19-dtb1).
The source of the measurements below is the 8-mode `[7/7]` Padé bath with the parameters above.
Its Hankel singular values are `0.429, 0.118, 1.62e-2, 2.02e-3, 1.66e-4, 7.91e-6,
1.86e-7, 1.53e-9`. The spectral error is relative to the exact thermal maximum on
`|ω| ≤ 6`. The equilibrium error is relative to the continuum variances above.

| Bath | Modes | Spectral error | Equilibrium error | `max ΔW`, depth 6 | Members | Time |
|---|---:|---:|---:|---:|---:|---:|
| Padé 7 (source) | 8 | `2.33e-5` | `3.73e-4` | — | 3003 | 207 s |
| Padé 3 | 4 | `1.67e-2` | `6.52e-3` | `1.43e-3` | 210 | 13 s |
| Truncated | 4 | `4.58e-4` | `5.59e-4` | `1.97e-6` | 210 | 14 s |
| Residualized | 4 | `3.54e-4` | `7.72e-4` | `1.87e-6` | 210 | 13 s |
| Truncated | 3 | `5.69e-3` | `8.67e-3` | `2.02e-5` | 84 | 5 s |

`max ΔW` is the final-state pointwise difference from the source hierarchy for the
quartic example at `t=0.8`. The compressed 4-mode bath changes by `4.39e-7` from depth
6 to 8. Its depth convergence resembles that of an uncompressed Padé bath of equal
size. Times are single-thread measurements and scale with member count. The source's
equilibrium error exceeds its spectral error because its white-noise Padé remainder is
not resolved by the window. A one-mode compression is stable and has a bounded
spectral change, yet its harmonic covariance violates `det Σ ≥ ħ²/4`. This is why
compression should be accepted only after spectral, harmonic-equilibrium and depth
checks. Residualization cannot be applied when the eliminated real part is negative,
as for fast Brownian Matsubara terms. Use truncation there.

## AAA rational fits

`aaa_bath(J; kT, frequencies)` fits the thermal noise spectrum of a general spectral
density, following the free-pole HEOM of
[Xu, Yan, Shi, Ankerhold and Stockburger, PRL 129, 230601 (2022)](https://doi.org/10.1103/PhysRevLett.129.230601).
For `ω > 0`,

```math
S(\pm\omega)=E(\omega)\pm\omega\,O(\omega),\qquad
E=\hbar J\coth\frac{\hbar\omega}{2kT},\qquad O=\frac{\hbar J}{\omega}.
```

`E` and `O` are even, so both are fitted as rational functions of `u = ω²` with a
common barycentric denominator (set-valued AAA). The two Loewner blocks are weighted
as `E` and `ωO`, and greedy support selection maximizes `|ΔE| + ω|ΔO|`, which equals
`max(|ΔS(ω)|, |ΔS(-ω)|)`. A pole `u` gives the rate `z = sqrt(-u)` with `Re z > 0`, so
complex rates come in conjugate pairs. Fitting `S` directly in `ω` would give lower
half-plane poles without this pairing, so a damped oscillation would generally cost
four modes instead of two. Poles on the positive real `u` axis lie on the real frequency
axis and are discarded.

The finite eigenvalues of the arrowhead pencil lose relative accuracy for poles much
smaller than the largest support point. With samples spanning twelve decades in `u`,
the slowest cold Drude poles are only accurate to `5e-3`. Newton steps on the
barycentric denominator, which is evaluated accurately near each pole, restore
partial-fraction agreement with the barycentric fit to about `1e-12`. For fixed rates, the
real part of each mode response is even in `ω` and the imaginary part odd. Least
squares therefore fits `real(c)` to `E` and `imag(c)` to `ωO` separately, with a
constant even term as white-noise diffusion when it is nonnegative. The degree rises
until this refitted bath, not merely the barycentric interpolant, meets `reltol`.
The counterterm `(1/π)∫J/ω dω` uses exp-sinh quadrature over `10^±30` times the
geometric mean sample frequency.

For the cold Drude parameters above, with 600 logarithmic samples on `[1e-3, 1e3]`
and the same error measures as the compression table:

| Bath | Modes | Spectral error | Equilibrium error |
|---|---:|---:|---:|
| Padé 7 | 8 | `2.33e-5` | `3.73e-4` |
| Padé 30 | 31 | `1.58e-15` | `1.79e-6` |
| AAA, `reltol = 1e-2` | 4 | `1.15e-3` | `6.48e-4` |
| AAA, `reltol = 1e-3` | 5 | `3.35e-4` | `5.99e-4` |
| AAA, `reltol = 1e-4` | 8 | `1.63e-5` | `1.07e-6` |
| AAA, `reltol = 1e-6` | 11 | `6.79e-7` | `1.50e-7` |
| AAA, `reltol = 1e-8` | 16 | `3.98e-10` | `1.06e-8` |
| AAA `1e-8`, truncated | 4 | `4.48e-3` | `8.11e-3` |

The fitted rates include the Drude pole and the first Matsubara frequency, to
`4e-10` and `8e-7` relative error. The static imbalance between the counterterm and
the fitted `imag ∫C dt` is `1e-12` at `reltol = 1e-8`. At equal mode count, the fit
reaches the continuum equilibrium more than 300 times more closely than Padé. Direct
four- and five-mode fits are comparable to the balanced Padé reductions above. A
four-mode truncation of the 16-mode fit is worse, because its Hankel singular values
also weigh frequencies up to `10³` that the oscillator never resolves.

For a Drude background with an underdamped vibration (`λ = 0.3, γ = 0.5` and
`λ = 0.2, ω₀ = 1.5, γ = 0.2`) at `kT = 0.1`, `reltol = 1e-6` gives 12 modes.
The vibration is a single two-mode block with rate `γ/2` and frequency
`sqrt(ω₀² - γ²/4)`, both to `2e-8`. The harmonic equilibrium agrees with the
continuum result to a relative `3.3e-7`. That result is a Matsubara sum with
closed-form friction kernels, which matches the real-axis Drude integral to `3e-11`.
Compressing the fit to eight modes keeps the equilibrium within `1.4e-4`. A depth-6
hierarchy with a four-mode compression reproduces its exact Gaussian dynamics to
`2e-7`.

For sub-Ohmic `J = 0.6√(2ω) exp(-ω/2)` at `kT = 0.1`, the quadrature reproduces
`λ = 1.2/√π` to rounding despite the `ω^{-1/2}` singularity of `J/ω`. Fits sampled
from `ω = 10^{-2}, 10^{-3}, 10^{-4}` all meet `reltol = 1e-4`, with 14, 13 and 15 modes.
Their static imbalances are `1.4e-2`, `5.0e-3` and `1.6e-3`, from `J/ω` below the
lowest sample. The equilibrium errors are `3.8e-3`, `1.3e-3` and `3.5e-4` relative to
the continuum. For the two narrower ranges, the `⟨q²⟩` error is the classical shift
`2kTΔ/(mω²)²` to within `4%`. Their thermalization is limited by the sampled range,
not by the spectral tolerance, and the static imbalance diagnoses it. From `10^{-4}`,
the larger error is in `⟨p²⟩` and comes from the tolerance: with 1000 samples,
`reltol = 1e-5` uses 24 modes and reduces it from `2.6e-4` to `3e-6`, leaving the
classical `⟨q²⟩` shift to within `1%`. Samples spanning many more decades in `u` make
the slowest poles ill-conditioned in Float64; such fits fail `reltol` and throw
instead of returning an inaccurate bath.

## Verification

The complete package suite passes **4,315 tests**, with **1,651/1,651 executable
source lines covered (100%)**. The standalone seven-run quartic example also passes its
assertions and produces the convergence changes reported above. JuliaFormatter and
`git diff --check` pass. Coverage records exercised code; the independent physical
references and convergence studies establish the numerical checks described here.
