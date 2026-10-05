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

## Verification

The complete package suite passes **2,542 tests**, with **854/854 executable source
lines covered (100%)**. The standalone seven-run quartic example also passes its
assertions and produces the convergence changes reported above. JuliaFormatter and
`git diff --check` pass. Coverage records exercised code; the independent physical
references and convergence studies establish the numerical checks described here.
