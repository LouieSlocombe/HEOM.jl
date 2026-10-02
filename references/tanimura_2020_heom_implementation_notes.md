# Tanimura (2020): Hierarchical equations of motion

## Source, purpose, and scope

**Paper:** Yoshitaka Tanimura, *Numerically "exact" approach to open quantum dynamics: The hierarchical equations of motion (HEOM)*, **The Journal of Chemical Physics 153**, 020901 (2020). DOI: `10.1063/5.0011599`.

**Source file:** `020901_1_5.0011599.pdf`.

**Purpose:** Repository context for designing and beginning an HEOM implementation, especially finite-temperature bosonic baths and continuous reaction coordinates. This is a perspective/review, not a new solver release or a complete specification of every variant it surveys. The paper states that it created or analyzed no new data. [T20, abstract; Data Availability]

References such as **[T20, Sec. III B, Eq. (5)]** refer to this paper. Printed article page 1 is PDF page 2; the first PDF page is a publisher cover. The article runs to printed page 22. Reference numbers such as **[T20, Ref. 44]** identify papers in its bibliography, not additional sources consulted for this summary.

This document separates three kinds of material:

- **Paper:** What the uploaded article explicitly presents or reports.
- **Derived:** Algebraic consequences and convention checks added here to make an implementation unambiguous; these are not additional equations quoted from the article.
- **Design:** Proposed software structure, discretizations, and tests; the paper does not prescribe them.

**Important:** Several printed expressions have internal sign or normalization inconsistencies. Section 13 records them, including the original expressions. The proposed reference kernel below fixes its own conventions explicitly; it is **not** an author-approved erratum. A small kernel-level check is documented in Section 10, but the paper's full benchmark suite has not been reproduced here.

## 1. Main message and paper structure

**Paper.** HEOM replaces a nonlocal reduced-dynamics problem with a coupled hierarchy of time-local differential equations. The physical reduced density operator is the hierarchy's zeroth member. Auxiliary density operators (ADOs) retain information about bath memory and system-bath correlations, called **bathentanglement** in the paper. The physical member includes contributions from all orders of the interaction; higher tiers are not simply successive perturbative corrections to an uncoupled root. [T20, Secs. I and III B; Fig. 2]

For the specified system-bath model, a converged calculation can describe strong coupling, non-Markovian dynamics, finite-temperature quantum fluctuations, correlated equilibrium, and driven dynamics without a weak-coupling or rotating-wave approximation. "Numerically exact" means controllable numerical accuracy, supported by comparison with analytical results, not that an arbitrary finite hierarchy is exact. [T20, Secs. I, III, and IV]

The paper's organization is useful for implementation:

| Source section | Main content | Repository relevance |
| --- | --- | --- |
| II | Hamiltonian, harmonic baths, spectral distributions | Model and units contract |
| III A-B | Quantum noise and density-operator HEOM | Bath decomposition and RHS |
| III C-E | Wigner, classical, and multistate equations; truncation | Coordinate-space extensions |
| III F-G | Correlated states, positivity, numerical techniques | Initialization and performance |
| IV | Four Brownian-oscillator tests | Scientific acceptance criteria |
| V | Alternative HEOM formulations | Extension boundaries |
| VI-VII | Applications and future directions | What the method can address |
| Appendix | Nonlinear response functions | Pulse/propagate/readout workflow |

**Design.** Start with a finite-dimensional, deterministic, single-bosonic-bath Drude implementation. Validate that kernel before adding coordinate grids, multiple baths, optimized decompositions, or nonlinear spectra. A single-bath prototype is not an implementation of every HEOM variant in the review.

## 2. Physical model and required inputs

### 2.1 Hamiltonian

**Paper.** The total Hamiltonian is

$$
H_{\mathrm{tot}}(t)=H_A(t)+\sum_a\sum_j
\left[
\frac{p_{aj}^2}{2m_{aj}}+
\frac{m_{aj}\omega_{aj}^2x_{aj}^2}{2}
-\alpha_{aj}V_a x_{aj}
\right].
$$

Here, $H_A$ is the system Hamiltonian; $V_a$ is its coupling operator to bath $a$; and the bath coordinates are harmonic oscillators. Different baths are assumed independent. [T20, Eqs. (1)-(2), Sec. II]

The system may be a discrete-state Hamiltonian, a coordinate Hamiltonian

$$
H_A(t)=\frac{p^2}{2m}+U(q;t),
$$

or a matrix of coupled electronic potential-energy surfaces. The system potential need not be harmonic. Likewise, a nonlinear **system-side** coupling such as $V(q)=q+cq^2$ does not itself make a harmonic, linearly coupled bath non-Gaussian. The review discusses both this coupling and $V(q)=e^{-cq}$. [T20, Secs. II, III C, and III E]

The standard construction relies on Gaussian bath statistics: higher bath cumulants vanish, so two-point bath correlations determine the influence. An anharmonic bath, nonlinear dependence on bath coordinates, or correlated interacting electron bath requires additional treatment; it is not covered merely by choosing a complicated $V(q)$. [T20, Secs. II and V F]

### 2.2 Spectral distribution and temperature

**Paper.** The spectral distribution function uses the normalization

$$
J_a(\omega)=\sum_j
\frac{\hbar\alpha_{aj}^2}{2m_{aj}\omega_{aj}}
\delta(\omega-\omega_{aj}),
\qquad
\beta=\frac{1}{k_B T}.
$$

The continuum bath represents an effectively unlimited thermal reservoir. [T20, Eq. (3), Sec. II]

**Design.** Make the following explicit model inputs: $H_A(t)$, each $V_a$, bath temperature, spectral parameters and convention, initial preparation, observables, and counterterm policy. Keep the system representation separate from the bath decomposition. A potential curve alone is not a complete open-system model.

Use either a coherent physical unit system or fully nondimensionalized quantities. Every occurrence of $H/\hbar$, $\gamma$, and $\nu_k$ must have units of inverse time. The paper uses angular frequencies; do not mix them directly with energies or spectroscopic wavenumbers.

### 2.3 Counterterm: a model choice, not a numerical stabilizer

**Paper.** For Brownian coordinate models the paper introduces

$$
H_{\mathrm{ct}}=\sum_{a,j}
\frac{\alpha_{aj}^2 V_a^2}{2m_{aj}\omega_{aj}^2}.
$$

It completes the oscillator squares and preserves translational symmetry for a free coordinate. The paper explicitly distinguishes its counterterm-containing QHFPE from a straightforward Wigner transformation of the regular HEOM. [T20, Sec. II; Sec. III C, continuing onto printed p. 7]

**Derived.** Under Eq. (3)'s normalization, write

$$
H_{\mathrm{ct}}=\sum_a\kappa_a V_a^2,
\qquad
\kappa_a=\frac{1}{\hbar}\int_0^\infty
\frac{J_a(\omega)}{\omega}\,d\omega.
$$

For the Drude spectrum below, $\kappa=\eta\gamma/2$.

**Design.** In the matrix reference kernel use

$$
H_{\mathrm{sys}}(t)=H_A(t)+H_{\mathrm{ct}}
$$

when modeling the completed-square Hamiltonian. Otherwise use the specified bare $H_A$. Record whether the supplied potential already contains this contribution. Do not add it twice or compare counterterm and no-counterterm models as though they were identical.

## 3. Drude bath and quantum thermal memory

### 3.1 Expressions presented in the paper

**Paper.** The main worked bath is

$$
J(\omega)=\frac{\hbar\eta}{\pi}
\frac{\gamma^2\omega}{\gamma^2+\omega^2}.
$$

The parameter $\gamma$ sets the inverse bath-memory time; $\eta$ sets the coupling strength. For coordinate coupling the paper also uses $\zeta=\eta/m$. [T20, Eq. (4); Eq. (9)]

Writing $\Omega=\sum_j\alpha_jx_j$, the bath response and symmetric correlation are defined as

$$
L_1(t)=\frac{i}{\hbar}\langle[\Omega(t),\Omega(0)]\rangle_B,
\qquad
L_2(t)=\frac12\langle\{\Omega(t),\Omega(0)\}\rangle_B.
$$

The Drude expansion for the symmetric correlation is

$$
L_2(t)=c_0e^{-\gamma t}+\sum_{k=1}^{\infty}c_ke^{-\nu_k t},
\qquad t>0,
$$

with

$$
\nu_k=\frac{2\pi k}{\beta\hbar},
\qquad
c_0=\frac{\hbar\eta\gamma^2}{2}
\cot\left(\frac{\beta\hbar\gamma}{2}\right),
$$

$$
c_k=\frac{2\eta\gamma^2\nu_k}{\beta(\nu_k^2-\gamma^2)},
\qquad k\ge1.
$$

The last expression is algebraically identical to the paper's version with a minus sign and denominator $\gamma^2-\nu_k^2$. The paper's separate printed response coefficient $\bar c_0$ needs the convention check in Section 13. [T20, Secs. II and III, printed pp. 3-4]

### 3.2 Why the Matsubara terms matter

**Paper.** There are two sources of memory: the mechanical bath decay $e^{-\gamma t}$, and quantum thermal terms $e^{-\nu_k t}$. At low temperature $\nu_1$ decreases, so thermal memory can remain important even for a mechanically fast bath. Figure 1, PDF page 5, shows negative portions of the symmetric correlation at low temperature. These are not negative probabilities and must not be clipped. [T20, Sec. III A; Fig. 1]

Taking a fast-bath limit does not by itself justify a low-temperature Markovian approximation. The paper obtains its high-temperature Markovian QFPE only while maintaining $\beta\hbar\gamma\ll1$, and warns that using it outside that regime can violate positivity. [T20, Secs. III A and III D, Eq. (13)]

### 3.3 Convention-explicit complex correlation for the reference kernel

**Derived from the oscillator Hamiltonian and Eq. (3).** Define

$$
C_B(t)=\langle\Omega(t)\Omega(0)\rangle_B
=\int_0^\infty J(\omega)
\left[
\coth\left(\frac{\beta\hbar\omega}{2}\right)\cos\omega t
-i\sin\omega t
\right]d\omega.
$$

This fixes all normalization factors:

$$
C_B(t)=L_2(t)-\frac{i\hbar}{2}L_1(t).
$$

For the Drude spectrum, define real positive decay rates $r_k$ and complex coefficients $d_k$ by

$$
r_0=\gamma,\qquad r_k=\nu_k\quad(k\ge1),
$$

$$
d_0=\frac{\hbar\eta\gamma^2}{2}
\left[\cot\left(\frac{\beta\hbar\gamma}{2}\right)-i\right],
\qquad d_k=c_k\quad(k\ge1).
$$

Then

$$
C_B(t)\approx\sum_{k=0}^{K}d_ke^{-r_k t}.
$$

These definitions imply $L_1(t)=\eta\gamma^2e^{-\gamma t}$ for $t>0$. They are the **derived implementation convention**, not a transcription of the conflicting printed $\bar c_0$ formulas.

**Design safeguards.** Require $\beta>0$, $\gamma>0$, $\eta\ge0$, and integer $K\ge0$. When $\gamma$ approaches a Matsubara frequency, individual terms become singular while their combined limit requires cancellation. Reject such parameters in an initial coefficient builder with a clear error; later implement the combined limiting expression. Do not silently shift $\gamma$ or cap the coefficients.

**Derived endpoint caveat.** Eq. (4) has a $1/\omega$ high-frequency tail. Its ideal quantum equal-time force variance is not finite: the Matsubara coefficients decay as $1/k$. Compare correlations at $t>0$ or compare integrated kernels, rather than requiring convergence of a finite value of $C_B(0)$.

## 4. Density-operator hierarchy

### 4.1 Indexing and the equation as printed

**Paper, with index relabeling.** Combine the paper's indices $(n,j_1,\ldots,j_K)$ into

$$
\mathbf n=(n_0,n_1,\ldots,n_K),
\qquad
|\mathbf n|=\sum_{k=0}^{K}n_k,
\qquad
\Gamma_{\mathbf n}=\sum_{k=0}^{K}n_kr_k.
$$

The physical reduced state is $\rho_A=\rho_{\mathbf0}$. The other members are ADOs, not independently normalized physical density matrices. Figure 2, PDF page 6, illustrates how they retain interaction histories. [T20, Sec. III B; Fig. 2]

Equation (5) has the structure

$$
\dot\rho_{\mathbf n}
=-\left(i\mathcal L_A+\Gamma_{\mathbf n}+\Xi\right)\rho_{\mathbf n}
+\sum_{k=0}^{K}\Phi\rho_{\mathbf n+\mathbf e_k}
+\sum_{k=0}^{K}n_k\Theta_k\rho_{\mathbf n-\mathbf e_k},
$$

where $\mathcal L_A X=[H_A,X]/\hbar$ and

$$
V^\times X=VX-XV,\qquad V^\circ X=VX+XV.
$$

**As printed**, its operators are

$$
\Phi=-\frac{i}{\hbar}V^\times,
\qquad
\Theta_0=-\frac{i}{\hbar}(c_0V^\times-i\bar c_0V^\circ),
\qquad
\Theta_k=-\frac{ic_k}{\hbar}V^\times,
$$

$$
\Xi=-\left(\sum_{k>K}\frac{c_k}{\nu_k}\right)
\frac{V^\times V^\times}{\hbar^2}.
$$

The printed $\bar c_0$ normalization and the sign of $\Xi$ should **not** be copied into code without resolving the issues in Section 13.

### 4.2 Proposed reference equation

**Derived convention for implementation.** With $C_B(t)$, $d_k$, and $r_k$ defined in Section 3.3, use

$$
\boxed{
\begin{aligned}
\dot\rho_{\mathbf n}={}&
-\frac{i}{\hbar}[H_{\mathrm{sys}}(t),\rho_{\mathbf n}]
-\Gamma_{\mathbf n}\rho_{\mathbf n}\\
&-\frac{i}{\hbar}\sum_{k=0}^{K}
[V,\rho_{\mathbf n+\mathbf e_k}]\\
&-\frac{i}{\hbar}\sum_{k=0}^{K}n_k
\left(d_kV\rho_{\mathbf n-\mathbf e_k}
-d_k^*\rho_{\mathbf n-\mathbf e_k}V\right)
+\mathcal D_{\mathrm{tail}}\rho_{\mathbf n}.
\end{aligned}}
$$

The sign of the paper's $-V\Omega$ interaction can be absorbed into the definition of the bath force and ADOs; the two-point correlation is unchanged. The above choice makes the upward coupling $-i[V,\cdot]/\hbar$. The choice must be held fixed if bath observables are later reconstructed from ADOs.

The downward term is equivalently

$$
-\frac{i n_k}{\hbar}
\left[\operatorname{Re}(d_k)V^\times
+i\operatorname{Im}(d_k)V^\circ\right]
\rho_{\mathbf n-\mathbf e_k}.
$$

Only the Drude mode has a nonzero imaginary coefficient here. That term carries the bath response; dropping it would not describe the same dissipative model.

**Scope restriction:** This simple use of $d_k^*$ assumes real decay rates. A decomposition with complex-conjugate decay pairs needs consistent pairing of the expansions of $C_B$ and $C_B^*$. Do not generalize this kernel by merely allowing complex numbers in the rate array.

### 4.3 Optional fast-Matsubara residual

**Paper.** The article compensates omitted fast terms with a local renormalization operator. This is separate from cutting off the hierarchy depth. [T20, Secs. III B and III D]

**Derived convention.** If the omitted Matsubara terms are fast relative to the dynamics of interest, their local contribution in the reference equation is

$$
\mathcal D_{\mathrm{tail}}X
=-\frac{\Delta_K}{\hbar^2}[V,[V,X]],
\qquad
\Delta_K=\sum_{k=K+1}^{\infty}\frac{c_k}{\nu_k}.
$$

Using the Drude coefficients,

$$
\Delta_K=
\frac{\eta}{\beta}-\frac{c_0}{\gamma}
-\sum_{k=1}^{K}\frac{c_k}{\nu_k}.
$$

This follows from the partial-fraction expansion of the cotangent. It avoids summing an infinite residual explicitly, but numerical cancellation is possible near poles or when the residual is very small.

For a positive residual, this term damps coherences in the $V$ eigenbasis. Its sign can also be obtained by eliminating a fast ADO from the upward and downward couplings. This is opposite to inserting the printed $\Xi$ into the printed $-\Xi$ term without modification.

**Design.** Provide two explicit modes: `tail = none`, converging the explicit Matsubara expansion, and `tail = fast_matsubara`, using the above approximation only after checking its timescale separation. The residual is not a substitute for sufficient $K$, and its sign should not be silently clipped.

## 5. Minimal implementation architecture

Everything in this section is **Design**, based on the hierarchy structure rather than a supplied software implementation.

### 5.1 State and connectivity

Let $D$ be the system Hilbert-space dimension, $M=K+1$ the number of exponentials, and $L$ the chosen maximum total tier. Keep all nonnegative multi-indices satisfying $|\mathbf n|\le L$.

The number of retained ADOs is

$$
N_{\mathrm{ADO}}=\binom{M+L}{L}.
$$

For example, $M=6,L=6$ gives 924 ADOs; $M=20,L=6$ gives 230,230. A `complex128` state needs approximately $16N_{\mathrm{ADO}}D^2$ bytes, before integrator work arrays and saved states. This combinatorial count is derived from the proposed total-tier cutoff; it is not a performance measurement from the paper.

Store:

```text
rho[A, i, j]        complex ADO matrices; root A = 0
occupation[A, k]   integer multi-indices
up[A, k]           index of n + e_k, or an explicit absent marker
down[A, k]         index of n - e_k, or an explicit absent marker
damping[A]         sum_k occupation[A, k] * rate[k]
rate[k]            real positive decay rates
coefficient[k]     complex d_k values
```

Enumerate weak compositions tier by tier rather than constructing and filtering the full Cartesian product. The retained index set must be downward closed. Precompute connectivity, but not time-dependent Hamiltonians.

### 5.2 RHS pseudocode

The following is language-independent pseudocode for Section 4.2, **not executable Python**. The absent upper neighbors implement a simple hard tier cutoff; they do not implement the paper's hybrid terminator.

```text
function rhs(t, rho):
    H = H_system(t)                         # includes counterterm if chosen
    drho = zero_array_like(rho)

    for A in all_ADOs:
        X = rho[A]
        drho[A] = -(i / hbar) * (H @ X - X @ H)
        drho[A] -= damping[A] * X

        if tail_is_enabled:
            C = V @ X - X @ V
            drho[A] -= (Delta_K / hbar**2) * (V @ C - C @ V)

        for k in all_bath_modes:
            if up[A, k] exists:
                Y = rho[up[A, k]]
                drho[A] -= (i / hbar) * (V @ Y - Y @ V)

            if occupation[A, k] > 0:
                Y = rho[down[A, k]]
                d = coefficient[k]
                drho[A] -= (i / hbar) * occupation[A, k] * (
                    d * (V @ Y) - conjugate(d) * (Y @ V)
                )

    return drho
```

All RHS entries must read the same input state. Do not update ADOs sequentially in place. An absent marker such as `-1` must be checked before array indexing.

For independent baths, concatenate the $(a,k)$ modes, use the matching $V_a$ for each mode, and sum their residual operators. There is one joint hierarchy, including mixed-bath multi-indices, not separate independently propagated system density matrices. This is a derived extension of the independent-bath model in Sec. II.

### 5.3 Suggested module boundaries

```text
models/       H_system(t), coupling operators, counterterm policy, units
baths/        spectral conventions, decomposition coefficients, residuals
hierarchy/    index enumeration, neighbors, optional rescaling
rhs/          matrix implementation; later coordinate or Wigner backend
propagation/  integration, equilibrium preparation, restart checkpoints
observables/  root-state measurements and hierarchy-wide pulse operations
tests/        algebraic checks, convergence sweeps, paper benchmarks
```

An initial explicit Runge-Kutta implementation is straightforward, but large $\Gamma_{\mathbf n}$ values can make the equations stiff. Reassess the integrator as $K$ and $L$ increase. Evaluate $H(t)$ at the integrator's internal stage times. Integration tolerance is only one error source; it does not control bath-expansion or hierarchy errors.

Save observables frequently and the full hierarchy at restart or response-branch points. Saving every ADO at every output time can dominate memory and storage.

## 6. Preparation, equilibrium, and observables

### 6.1 Factorized preparation

**Design, consistent with the factorized preparation discussed in the paper.** For a system state initially uncorrelated with a thermal bath, start with

$$
\rho_{\mathbf0}(0)=\rho_A(0),\qquad
\rho_{\mathbf n\ne\mathbf0}(0)=0.
$$

Require unit trace, Hermiticity, and positive semidefiniteness of the initial root. Do not normalize higher ADOs individually. This preparation includes the transient creation of system-bath correlations; it is not generally coupled thermal equilibrium. [T20, Secs. III D, III F, and IV]

### 6.2 Correlated equilibrium

**Paper.** Propagate the undriven HEOM until **all** hierarchy elements become stationary, then retain that entire hierarchy as the correlated initial state. The higher ADOs encode correlations that are needed for subsequent dynamics and nonlinear response. [T20, Sec. III F]

**Derived normalization made explicit.** The target root is

$$
\rho_A^{\mathrm{eq}}=
\frac{\operatorname{Tr}_B e^{-\beta H_{\mathrm{tot}}}}
{\operatorname{Tr}_{A+B}e^{-\beta H_{\mathrm{tot}}}},
$$

not generally $e^{-\beta H_A}/\operatorname{Tr}_A e^{-\beta H_A}$. The source displays the partial Boltzmann trace without its normalization denominator.

**Design.** Check an appropriately scaled full-hierarchy RHS residual and changes in observables over a time window. Converged root populations alone do not establish a converged correlated state. Save the full hierarchy and its model/decomposition metadata together.

**Derived qualification.** Relaxation from an arbitrary state gives a unique thermal state only when the model actually equilibrates across the relevant sectors. Pure dephasing with conserved populations, for example, does not redistribute arbitrary initial populations into Gibbs weights. For periodic driving, a limit cycle rather than a time-independent equilibrium may be the relevant long-time state.

**Paper.** Imaginary-time HEOM provide another route to equilibrium and thermodynamic quantities. Their equations are not supplied in full here; the review points to Refs. 43-44. A driven steady state still requires real-time treatment. [T20, Sec. V G]

### 6.3 Readout

For a system operator $O$,

$$
\langle O\rangle(t)=\operatorname{Tr}_A[O\rho_{\mathbf0}(t)].
$$

Compute populations, coherences, coordinate moments, and system energy from this root. Bath observables and strong-coupling heat currents generally need additional ADO expressions. The paper warns that heat definitions based only on system-energy changes can give misleading thermodynamic conclusions; it does not provide a complete heat-current implementation in this review. [T20, Secs. III B and VI D]

## 7. Coordinate-space and Wigner implementations

### 7.1 Direct position-basis route

**Paper.** The system Hamiltonian may explicitly contain a continuous coordinate and an arbitrary potential. [T20, Sec. II]

**Derived.** For $H_A=p^2/(2m)+U(q;t)$, the Hamiltonian part of each ADO equation in the position basis is

$$
\left.\partial_t\rho_{\mathbf n}(q,q')\right|_H
=\frac{i\hbar}{2m}(\partial_q^2-\partial_{q'}^2)\rho_{\mathbf n}(q,q')
-\frac{i}{\hbar}[U(q;t)-U(q';t)]\rho_{\mathbf n}(q,q').
$$

For coordinate-diagonal coupling,

$$
[V,X](q,q')=[V(q)-V(q')]X(q,q'),
$$

$$
\{V,X\}(q,q')=[V(q)+V(q')]X(q,q').
$$

Add the counterterm to the potential differences when that is the selected model.

**Design.** A first grid implementation can reuse the matrix kernel. On an orthonormalized uniform interior grid with zero wavefunction boundary values, a second-order finite-difference kinetic operator has

$$
T_{ii}=\frac{\hbar^2}{m\Delta q^2},
\qquad
T_{i,i\pm1}=-\frac{\hbar^2}{2m\Delta q^2},
$$

with $U_{ij}=U(q_i)\delta_{ij}$ and $V_{ij}=V(q_i)\delta_{ij}$. These discretization choices are not specified by the paper.

In this convention $\operatorname{Tr}\rho=1$ and $P(q_i)\approx\rho_{ii}/\Delta q$. A raw coordinate-kernel representation instead carries quadrature weights in its trace. Do not mix these conventions.

This route accepts a tabulated reaction-coordinate potential without deriving a new bath hierarchy. Supply an effective mass, coupling function, bath parameters, preparation, and boundaries as well. Increase the grid extent and resolution independently; verify that the dynamics do not reflect from artificial boundaries over the observation window.

### 7.2 Wigner representation

**Paper.** Equation (6) defines

$$
W_{\mathbf n}(p,q;t)=\frac{1}{2\pi\hbar}\int_{-\infty}^{\infty}
 e^{ipx/\hbar}\rho_{\mathbf n}(q-x/2,q+x/2;t)\,dx.
$$

The physical Wigner function is real but may be negative. Its position marginal is $P(q)=\int W_{\mathbf0}(p,q)\,dp$. Wigner-space evolution is useful for boundary conditions and quantum/classical comparisons. [T20, Sec. III C]

The Hamiltonian generator is [T20, Eq. (7)]

$$
\mathcal Q W=-\frac{p}{m}\partial_q W
-\frac{1}{\hbar}\int\frac{dp'}{2\pi\hbar}
U_W(p-p',q)W(p',q),
$$

$$
U_W(p,q)=2\int_0^\infty
\sin(px/\hbar)[U(q+x/2)-U(q-x/2)]\,dx.
$$

**Derived sign check.** For smooth potentials for which the expansion is appropriate,

$$
\mathcal QW=-\frac{p}{m}\partial_q W+
\sum_{\ell=0}^{\infty}
\frac{(-1)^\ell(\hbar/2)^{2\ell}}{(2\ell+1)!}
U^{(2\ell+1)}(q)\partial_p^{2\ell+1}W.
$$

Consequently, the classical Hamiltonian limit is

$$
\mathcal Q_{\mathrm{cl}}W=-\frac{p}{m}\partial_qW+U'(q)\partial_pW.
$$

For a harmonic potential, higher Moyal terms vanish. Thermal quantum effects in the bath do **not** vanish for that reason. The printed classical Liouvillian following Eq. (9) has a sign issue noted in Section 13.

**Paper.** The QHFPE, Eq. (8), have the same hierarchy connectivity as Eq. (5), but act on phase-space functions and incorporate the counterterm. Nonlinear couplings and multistate Wigner functions are discussed in Secs. III C and III E. The multistate representation stores a matrix $W_{jk}(p,q)$ and uses the matrix potential $U_{jk}(q;t)$. [T20, Eqs. (8), (14)-(16)]

**Design limitation.** The finite-temperature operators printed below Eq. (8) have unresolved coefficient inconsistencies. Do not construct a production QHFPE solver by copying those coefficients from this review alone. Implementing the position-basis matrix route first avoids relying on that ambiguous operator listing. Exact reproduction of the author's QHFPE convention should consult the derivations cited in Refs. 44 and 89-90.

### 7.3 High-temperature and classical reference limits

**Paper.** For linear coordinate coupling and the stated high-temperature approximation, Eq. (9) is

$$
\partial_tW^{(n)}=
-(\mathcal L_{\mathrm{QM}}+n\gamma)W^{(n)}
+\partial_pW^{(n+1)}
+n\gamma\zeta\left(p+\frac{m}{\beta}\partial_p\right)W^{(n-1)}.
$$

Replacing the Hamiltonian generator by its classical limit gives the classical hierarchical Fokker-Planck equations. Taking the fast-bath limit while preserving the high-temperature condition gives [T20, Eq. (13)]

$$
\partial_tW=\mathcal QW+
\zeta\partial_p\left(pW+\frac{m}{\beta}\partial_pW\right).
$$

This is a useful limiting test, not a low-temperature replacement for the finite-temperature hierarchy. For $V=q$, the derived positive fast-Matsubara residual of Section 4.3 transforms into $+\Delta_K\partial_p^2W$, providing another sign check.

## 8. Correlation functions and nonlinear spectra

### 8.1 Hierarchy-wide operations

**Paper.** The Appendix gives a concrete algorithm: equilibrate, apply an interaction operator to **every** ADO, propagate for an inter-pulse interval, apply the next interaction, propagate again, then read the observable from the root. The paper's Raman-IR-IR example first applies $\Pi^\times$, then $\mu^\times$, and finally measures $\mu$. [T20, Appendix, printed pp. 18-19]

Resetting higher ADOs after a pulse discards system-bath correlations. Figures 7 and 15 illustrate the resulting loss of memory and spectral structure in approximate factorized treatments. A root-only regression calculation is therefore not interchangeable with the full-hierarchy procedure. [T20, Sec. III F; Appendix]

### 8.2 A consistent response convention

**Derived convention.** Let $\boldsymbol\rho$ denote the whole hierarchy and $\mathcal G(t)$ its undriven propagator. Define an explicit drive Hamiltonian

$$
H_{\mathrm{drive}}(t)=-f(t)\mu.
$$

The corresponding kick on every ADO is

$$
(\mathcal K_\mu\boldsymbol\rho)_{\mathbf n}
=\frac{i}{\hbar}[\mu,\rho_{\mathbf n}].
$$

Then, with the final subscript $\mathbf0$ selecting the root,

$$
R^{(1)}(t)=\operatorname{Tr}\left\{
\mu[\mathcal G(t)\mathcal K_\mu\boldsymbol\rho^{\mathrm{eq}}]_{\mathbf0}
\right\},
$$

$$
R_{\mathrm{TTR}}^{(2)}(t_2,t_1)=\operatorname{Tr}\left\{
\Pi[\mathcal G(t_2)\mathcal K_\mu\mathcal G(t_1)
\mathcal K_\mu\boldsymbol\rho^{\mathrm{eq}}]_{\mathbf0}
\right\},
$$

$$
R^{(3)}(t_3,t_2,t_1)=\operatorname{Tr}\left\{
\mu[\mathcal G(t_3)\mathcal K_\mu\mathcal G(t_2)\mathcal K_\mu
\mathcal G(t_1)\mathcal K_\mu\boldsymbol\rho^{\mathrm{eq}}]_{\mathbf0}
\right\}.
$$

The prefactors are thus $i/\hbar$, $-1/\hbar^2$, and $-i/\hbar^3$. The second- and third-order factors agree with the Appendix. Its printed first-order forms are inconsistent with one another and with Sec. IV; see Section 13.

For a symmetric coordinate autocorrelation, initialize the response hierarchy with

$$
X_{\mathbf n}(0)=\frac12(q\rho_{\mathbf n}^{\mathrm{eq}}
+\rho_{\mathbf n}^{\mathrm{eq}}q),
$$

propagate it, and measure $\operatorname{Tr}[qX_{\mathbf0}(t)]$. This is a derived application of the same operator/propagator construction.

**Design.** A response hierarchy is a linear-response object, not necessarily a normalized physical state. Do not renormalize it or project it onto positive density matrices.

### 8.3 Spectral assembly

**Paper.** Double Fourier transformation of the appropriate time intervals produces 2D spectra. For third-order spectra, transform $t_1$ and $t_3$ while holding the waiting time $t_2$ fixed. Specific Liouville pathways can be selected by left or right operator multiplication instead of the full commutator. [T20, Appendix; Figs. 14-15]

**Design.** Record the Fourier sign, angular-frequency conversion, selected pathways, time windows, and any apodization. These numerical choices are not fully specified by the review. Store the full hierarchy at branching times so each subsequent interval starts from the correct correlated state. A transform of the full commutator response is not automatically every experimentally phase-matched spectrum.

## 9. Truncation, convergence, and performance

### 9.1 Separate the approximation controls

**Paper.** The formal hierarchy is infinite. Sec. III D discusses time-nonlocal asymptotic closures, time-local closures for fast thermal terms, and a hybrid QHFPE scheme. The cutoff must be checked by increasing the retained hierarchy. Weak coupling can reduce the importance of high tiers; that argument alone is insufficient at strong coupling. [T20, Sec. III D]

**Design.** Keep separate controls for the number of correlation terms $K$, total tier depth $L$, system basis/grid, integration error, and preparation/observation times. Do not use one variable called `N` for all of them.

The source suggests $K\gg\omega_0/\nu_1$, where $\omega_0$ is a characteristic system frequency. It is a scale-separation guide, not an error guarantee for arbitrary models. The hybrid scheme also defines

$$
K_\gamma=
\begin{cases}
\operatorname{int}(K\nu_1/\gamma),&\nu_1>\gamma,\\
K,&\nu_1\le\gamma,
\end{cases}
$$

and gives boundary equations (11)-(12). These equations belong to the paper's QHFPE convention; the hard total-tier cutoff in Section 5 does not reproduce that scheme. [T20, Sec. III D]

A useful convergence record reports changes in the actual target observable when each control is tightened, including a joint $K,L$ check. Establish the tolerances before interpreting small population changes, rates, or spectral peaks.

### 9.2 Improvements discussed by the paper

**Paper.** Sec. III G discusses Pade and Fano decompositions, hierarchy rescaling, optimized hierarchical bases, tensor networks, GPU implementations, distributed-memory MPI, and low-storage Runge-Kutta integration. These are the review's 2020 numerical landscape, not a current ranking of software or integrators. [T20, Sec. III G]

**Design.** Optimize only after a transparent reference kernel passes tests. Batch matrix products and avoid assembling a dense global matrix over all $N_{\mathrm{ADO}}D^2$ unknowns. A matrix-free RHS naturally follows the sparse hierarchy graph.

A possible rescaling is

$$
\widetilde\rho_{\mathbf n}=\rho_{\mathbf n}/s_{\mathbf n},
\qquad
s_{\mathbf n}=\sqrt{\prod_k n_k!\,|d_k|^{n_k}}.
$$

This displayed choice is an implementation derivation, not a formula printed in the review. Upward edges gain $s_{\mathbf n+\mathbf e_k}/s_{\mathbf n}$ and downward edges gain $s_{\mathbf n-\mathbf e_k}/s_{\mathbf n}$. Handle zero coefficients explicitly. Rescaling alone does not change the number of ADOs or establish convergence; filtering adds another approximation to test.

## 10. Validation plan

### 10.1 Kernel and representation checks

The following are **Design/Derived** checks, distinct from reproducing the paper's figures.

| Check | Required behavior |
| --- | --- |
| Zero coupling | Root agrees with unitary evolution under the selected $H_A(t)$; zero-initialized ADOs stay zero |
| Trace | Root trace remains one for ordinary state propagation |
| Hermiticity | Root Hermiticity error decreases with numerical error controls |
| Positivity | Root minimum eigenvalue approaches a nonnegative value within converged numerical tolerance |
| Pure dephasing | Constant populations and the analytic coherence below |
| Counterterm | The implemented Hamiltonian matches the declared completed-square or bare model |
| Coordinate limits | Correct harmonic motion, boundary convergence, and phase-space normalization |
| Correlated restart | Restarting a saved full hierarchy reproduces uninterrupted propagation |
| Pulse operation | Direct weak-field evolution agrees with the chosen response sign convention |

Do not repair failed tests by clipping the physical state, forcing all ADO traces to zero, or deleting negative Wigner values. Positivity of the root is a different requirement from positivity of a phase-space quasiprobability. The paper also notes that its numerical positivity evidence is broader than the simple cases for which an analytical verification was available. [T20, Secs. III C and III F]

### 10.2 Analytic pure-dephasing check for the reference kernel

**Derived.** If $H_{\mathrm{sys}}$ and $V$ commute, let their simultaneous eigenvalues be $E_i$ and $v_i$. For factorized preparation and the specified exponential correlation define

$$
g(t)=\frac{1}{\hbar^2}\sum_{k=0}^{K}
\frac{d_k}{r_k^2}\left(e^{-r_kt}+r_kt-1\right).
$$

For the local residual option add $\Delta_Kt/\hbar^2$. The coherence is

$$
\rho_{ij}(t)=\rho_{ij}(0)
\exp\left[
-\frac{i(E_i-E_j)t}{\hbar}
-(v_i-v_j)^2\operatorname{Re}g(t)
-i(v_i^2-v_j^2)\operatorname{Im}g(t)
\right].
$$

The energies here must include any selected counterterm. A coupling with unequal squared eigenvalues, such as $V=\operatorname{diag}(0,1)$, tests the imaginary bath coefficient as well as dephasing. Using only eigenvalues $\pm v$ would miss that phase check.

**Preparation-time smoke test performed for this guide.** An independently coded matrix RHS using the displayed reference equation was compared with this expression for $H=\operatorname{diag}(0,0.4)$, $V=\operatorname{diag}(0,1)$, $\rho_{ij}(0)=1/2$, $\hbar=1.3$, $\beta=1.5$, $\eta=0.2$, $\gamma=0.7$, $K=2$, no counterterm, and $0\le t\le3$. Increasing the depth from $L=4$ to $L=6$ reduced the maximum sampled coherence error from approximately $1.23\times10^{-8}$ to $9.94\times10^{-13}$. The same check with the optional residual gave approximately $1.22\times10^{-8}$ and $9.85\times10^{-13}$. Trace and Hermiticity errors were zero at the recorded precision.

This checks indexing, factors of $\hbar$, the imaginary coefficient, and the residual sign for the chosen finite expansion. It does **not** establish convergence to the full Drude bath, thermalization accuracy, arbitrary-coupling validity at finite depth, or reproduction of the paper's benchmarks.

### 10.3 The paper's four scientific acceptance tests

**Paper.** Sec. IV uses a Brownian oscillator,

$$
H_A=\frac{p^2}{2m}+\frac12m\omega_0^2q^2,
$$

and compares QHFPE, analytical solutions, and TCL Redfield calculations. Figure 8, PDF page 12, provides four complementary tests:

| Test | Observable | What the paper uses it to test |
| --- | --- | --- |
| (a) | Steady-state distribution | Correlated thermal equilibrium |
| (b) | $C_q(t)=\langle\{q(t),q\}\rangle/2$ | Fluctuations |
| (c) | $R^{(1)}(t)=i\langle[q(t),q]\rangle/\hbar$ | Dissipation |
| (d) | $R_{\mathrm{TTR}}^{(2)}(t_2,t_1)=-\langle[[q^2(t_1+t_2),q(t_1)],q]\rangle/\hbar^2$ | Dynamical bathentanglement |

The caption's base parameters are $\omega_0=1$, $\gamma=1$, $\zeta=1$, and $\beta\hbar=3$, with stated changes including $\beta\hbar=1$ for the autocorrelation and $\zeta=3$ for the linear-response test. The harmonic linear response is temperature-independent. The figure shows agreement of QHFPE with the analytical results and missing structure in the approximate nonlinear response. [T20, Sec. IV; Fig. 8]

The article does not give the full analytical benchmark formulas, grid sizes, hierarchy depths, or integration tolerances needed to regenerate all panels. It points to the Brownian-oscillator literature and especially Ref. 44. Treat these as required follow-on validation sources, not as benchmark data already supplied by this PDF.

The author specifically argues that FMO exciton-transfer calculations are useful for scalability but are less discriminating than these tests for quantum thermal noise and general open-system accuracy. [T20, Sec. VI A]

## 11. Other HEOM families covered by the review

These are **Paper** summaries and extension boundaries, not implemented features of the reference kernel.

**Other spectral distributions.** Brownian, Lorentzian, super-Ohmic, Drude-Lorentz, and combined distributions can be treated using suitable decaying/oscillatory correlation decompositions. Fourier and quadrature-based extensions also appear. The review warns that a finite nondecaying basis may fail to approach long-time thermal equilibrium. A fit adequate over a short time interval is not automatically a thermalization model. [T20, Secs. III F and V A]

**Stochastic and wavefunction hierarchies.** Stochastic HEOM can trade memory for trajectory sampling. HOPS, stochastic Schrodinger hierarchies, and HSEOM reduce the per-state representation from a density matrix to wavefunctions, but introduce sampling or contour/basis costs. The paper discusses difficulties at low temperature and long times for the formulations it surveys; those comments should not be read as universal statements about later implementations. [T20, Secs. V B-C]

**Grand-canonical electron baths.** Fermionic HEOM introduce chemical potentials, electronic creation/annihilation operators, and fermionic thermal correlations. Equation (22) uses ordered auxiliary labels and parity signs. This is not the bosonic integer-occupation hierarchy with a changed spectral function. Electrode hybridization models, charge transport, and combinations with electronic-structure methods are discussed. [T20, Sec. V D, Eqs. (17)-(22)]

**Other models and non-Gaussian extensions.** The review covers Holstein and deformation-potential models, rotationally symmetric coupling, phenomenological extensions, and higher-cumulant corrections. Bath harmonicity/Gaussianity is an actual assumption of the standard solver. Preserving rotational symmetry or representing a non-Gaussian bath needs the corresponding specialized formulation. [T20, Secs. V E-F]

**Imaginary-time HEOM.** These give equilibrium/partition-function information and can support free-energy and entropy calculations. The review does not supply a turnkey imaginary-time RHS or a conversion procedure to every real-time ADO convention. [T20, Sec. V G]

## 12. Applications and what remains to be supplied

**Paper.** Applications include proton, electron, and exciton transfer; coupled electronic/nuclear dynamics; nonlinear vibrational and electronic spectra; quantum ratchets and resonant tunneling; molecular motors; quantum information; and heat transport/engines. [T20, Sec. VI]

Figure 5, PDF page 8, illustrates a proton-transfer model and a relationship between off-diagonal vibrational spectral peaks and a transition rate in that model. Figure 6 shows a driven molecular motor, and Figure 11 compares a nonlinear Brownian-oscillator description with molecular-dynamics-derived vibrational spectra. These are examples from cited studies, not complete new parameterizations provided by this review. [T20, Figs. 5, 6, and 11]

For a reaction-coordinate repository, the immediate implementation contract is therefore

```text
inputs:
    coordinate representation and effective mass
    potential U(q, t), or coupled surfaces U_jk(q, t)
    coupling V(q), or the corresponding operator matrix
    spectral distribution, temperature, and counterterm convention
    initial state or equilibrium-preparation procedure
    measured observable, dividing surface, or spectroscopic operators
    numerical controls and convergence targets

outputs:
    root reduced density matrix or Wigner function
    populations, moments, and selected response functions
    full-hierarchy checkpoints
    numerical diagnostics and convergence comparisons
```

**Design boundary.** A product-state population can be defined with a projector, for example $P_R(t)=\operatorname{Tr}(\Pi_R\rho_{\mathbf0})$, once a dividing surface is specified. A rate-extraction prescription, realistic potential, effective mass, bath fit, and molecular dipole/polarizability model remain application inputs. The review does not supply a universal reaction rate formula or an atomistic-to-HEOM parameter-fitting pipeline.

## 13. Source-equation audit: do not propagate these ambiguities silently

These checks were made against the rendered equations, not just extracted PDF text. They identify internal consistency issues in the uploaded article; resolving the author's intended notation definitively requires the cited derivations or an erratum.

### 13.1 Response coefficient and anticommutator normalization

On printed p. 3, the source defines $L_1=i\langle[\Omega(t),\Omega]\rangle/\hbar$. On printed p. 4, it writes

$$
L_1(t)=\bar c_0e^{-\gamma|t|},\qquad \bar c_0=\hbar\eta\gamma^2.
$$

But evaluation from its oscillator Hamiltonian and Eqs. (3)-(4) gives $L_1(t)=\eta\gamma^2e^{-\gamma t}$ for $t>0$. The complex force-correlation coefficient has imaginary part $-\hbar\eta\gamma^2/2$, whereas inserting the printed $\bar c_0$ in the printed $\Theta_0$ gives a different anticommutator strength. These are not interchangeable definitions. Section 3.3 deliberately starts from $C_B(t)$ to avoid that ambiguity.

### 13.2 Residual sign in the density hierarchy

Equation (5) contains $-\Xi\rho$, while the following definition is $\Xi=-\Delta_K V^\times V^\times/\hbar^2$. Together these give $+\Delta_K[V,[V,\rho]]/\hbar^2$. For positive $\Delta_K$ this amplifies rather than damps coherences in the $V$ eigenbasis. Eliminating a fast real-correlation mode gives the negative double commutator in Section 4.3. The reference kernel uses that explicitly derived sign.

### 13.3 Finite-temperature QHFPE coefficients

Below Eq. (8), the source prints

$$
\bar\Theta_0=\eta\left[\frac{p}{m}
+c_0\cot\left(\frac{\beta\hbar\gamma}{2}\right)\partial_p\right],
\qquad
\Xi'=-\frac{\eta}{\beta}\left(\sum_{k>K}c_k\right)\partial_p^2.
$$

The previously defined $c_0$ already contains a cotangent, and these formulas do not transparently reduce to Eq. (9) under those definitions. A different expression for $\Xi'$ appears after Eq. (12). No clear intervening coefficient redefinition is given. This guide therefore preserves the structural QHFPE summary but does not pretend that this operator listing is an unambiguous finite-temperature implementation specification.

### 13.4 Classical drift sign

Following Eq. (9), the printed $\mathcal L_{\mathrm{cl}}$ has $+(\partial_qU)\partial_p$, while evolution is written with $-\mathcal L_{\mathrm{cl}}$. The Hamiltonian limit derived from Eq. (7) instead gives $-p\partial_q/m+U'\partial_p$. A harmonic-oscillator trajectory or the Wigner expansion detects this discrepancy. Section 7.2 displays the derived generator explicitly.

### 13.5 Linear-response prefactors

Sec. IV defines linear response with $+i/\hbar$ multiplying the commutator. The Appendix first writes a Heisenberg commutator without $i$ and then a Liouville form with $-i/\hbar$. These cannot all be the same convention with the stated $O^\times X=[O,X]$. Section 8 defines the external-field Hamiltonian and derives the kick from it. This also fixes signs for later response orders.

### 13.6 Conventions that should not be treated as numerical guarantees

The unnormalized partial Boltzmann trace in Sec. III F needs its normalization for a density matrix. The hierarchy scale arguments in Secs. III B-D do not provide a universal cutoff or error estimate. Finally, the review's statements about numerical exactness and positivity apply to properly represented and converged physical dynamics, not to arbitrary negative/complex coefficients, arbitrary ADO initial conditions, or an unconverged hard cutoff.

## 14. Implementation milestones

**Design.** A practical sequence is:

1. **Reference matrix kernel.** Fix units and counterterms; implement Drude coefficients, total-tier indexing, the RHS, factorized initialization, diagnostics, and the dephasing test. Keep the equation convention next to the code.
2. **Scientific validation.** Add equilibration, full-hierarchy checkpoints, independent $K/L$/basis/time-step studies, and the four Brownian-oscillator tests from Sec. IV using their cited analytical sources.
3. **Target application.** Add a reaction-coordinate grid or multiple states, then hierarchy-wide response operations. Optimize decomposition, rescaling, storage, and parallel execution only while preserving the reference results.

A run record should contain the source/convention identifier, all model parameters, units, counterterm policy, exact correlation coefficients, hierarchy rule, solver settings, initial preparation, boundary conditions, and the relevant convergence comparisons. Without those, "HEOM" alone is not enough to reproduce a calculation.

## 15. Bibliography and follow-on implementation references

### Primary paper

```bibtex
@article{Tanimura2020HEOM,
  author  = {Tanimura, Yoshitaka},
  title   = {Numerically "exact" approach to open quantum dynamics:
             The hierarchical equations of motion (HEOM)},
  journal = {The Journal of Chemical Physics},
  volume  = {153},
  pages   = {020901},
  year    = {2020},
  doi     = {10.1063/5.0011599}
}
```

### References identified in the paper, not independently reviewed here

| T20 reference | Bibliographic locator as given in T20 | Why it is a follow-on source |
| --- | --- | --- |
| 42 | Y. Tanimura, J. Phys. Soc. Jpn. **75**, 082001 (2006) | Broader HEOM background and conventions |
| 43 | Y. Tanimura, J. Chem. Phys. **141**, 044114 (2014) | Correlated preparation and imaginary-time treatment |
| 44 | Y. Tanimura, J. Chem. Phys. **142**, 144110 (2015) | QHFPE conventions, equilibrium, and the main benchmark source |
| 86 | Y. Tanimura and R. Kubo, J. Phys. Soc. Jpn. **58**, 101 (1989) | Original hierarchy formulation |
| 88 | A. Ishizaki and Y. Tanimura, J. Phys. Soc. Jpn. **74**, 3131 (2005) | Low-temperature corrections and local truncation |
| 89-90 | Y. Tanimura and P. G. Wolynes, Phys. Rev. A **43**, 4131 (1991); J. Chem. Phys. **96**, 8485 (1992) | Coordinate-space hierarchy and reaction applications |
| 115 | J. Zhang, R. Borrelli, and Y. Tanimura, J. Chem. Phys. **152**, 214114 (2020) | Exponential-linear coupling and proton-transfer application |
| 125 | T. Ikeda and Y. Tanimura, J. Chem. Theory Comput. **15**, 2517 (2019) | Low-temperature multistate equations |
| 157-159 | J. Chem. Phys. **133**, 101106 (2010); **133**, 114112 (2010); **134**, 244106 (2011) | Pade decomposition methods |
| 162 | Q. Shi et al., J. Chem. Phys. **130**, 084105 (2009) | Hierarchy rescaling |
| 185 | Y.-a. Yan, Chin. J. Chem. Phys. **30**, 277 (2017) | Low-storage integration |
| 196 | J. Jin, X. Zheng, and Y. Yan, J. Chem. Phys. **128**, 234703 (2008) | Fermionic HEOM |

**Implementation takeaway:** Preserve the full hierarchy as the dynamical state; represent quantum bath memory rather than only mechanical damping; specify all coefficient and counterterm conventions; and establish convergence with equilibrium, fluctuation, dissipation, and nonlinear-response tests before interpreting application results.
