# Stabilized bosonic HEOM: paper summary and implementation specification

## Source and scope

**[P]** Salvatore Gatto, Samuel L. Rudge, Bokang Hou, Eran Rabani, and Michael Thoss, *Improving the Stability of the Hierarchical Equations of Motion for Open Quantum Systems with Strong Coupling to Structured Bosonic Baths*, arXiv:2609.35484v1, 28 September 2026. Source file: `2609.35484v1(1).pdf`, 17 pages.

This document summarizes the supplied version for use as research-code context. References such as **[P, Eq. (60), p. 9]** refer to that PDF, not to equation numbers in this document. Sections marked **Implementation note**, **Derived**, or **Proposed** contain algebraic translations or engineering choices, not additional results reported by the authors. No outside literature, software documentation, or unpublished data were used as sources.

**Implementation target:** a matrix-free solver for the transformed hierarchy in Eq. (60), with factorized initial conditions and physical density-matrix reconstruction using Eqs. (65)-(66). Start with a user-provided or independently generated complex bath-pole expansion; implement spectral-density fitting as a separate component.

**Critical contract:** the propagated transformed zeroth-tier ADO is generally **not** the physical reduced density matrix. Physical observables require a weighted sum over retained ADOs whose **individual occupations are all even**.

## 1. Main contribution and applicability

The paper addresses numerical instabilities caused by finite-tier truncation of bosonic hierarchical equations of motion (HEOM). Its analysis distinguishes spectral instability from non-normal transient amplification: even a generator with no eigenvalues in the right half-plane can strongly amplify errors over a finite time. Increasing the hierarchy depth does not necessarily remove this problem. [P, Secs. I and III A, pp. 1-2 and 4-7.]

The proposed construction is an auxiliary-space change of representation:

1. Rewrite the hierarchy in a bosonic Fock space with separate cosine/sine and left/right sectors.
2. Apply a non-unitary Gaussian similarity transformation and a Bogoliubov-like rotation.
3. Truncate and propagate in the transformed representation, then reconstruct physical observables.

The exact, untruncated dynamics is unchanged in the paper's formulation. The finite truncated generator changes, substantially reducing the instabilities in the demonstrated calculations. This is not a weak-coupling, Markovian, or modified-bath approximation. The demonstrated improvement is also not a theorem guaranteeing stability or positivity for every finite truncation. [P, Secs. III B-C and V, pp. 7-9 and 12.]

The assumptions are a Gaussian environment, coupling linear in bath operators, and an initially factorized system-bath state with a thermal bath. The system Hamiltonian and system coupling operator need not be two-dimensional. Multiple independent Gaussian coupling channels are explicitly allowed, although the displayed algorithm and benchmarks use one channel. The paper does not provide a correlated-initial-state implementation or a general non-Gaussian-bath extension. [P, Sec. II, pp. 2-3; Sec. III C, p. 9.]

## 2. Physical model and bath representation

### 2.1 Hamiltonian and initial state

The model is

$$
H=H_S+H_B+H_{SB},\qquad
H_B=\sum_k\omega_k b_k^\dagger b_k,
$$

$$
H_{SB}=O_S\otimes O_B,\qquad
O_B=\sum_k\lambda_k(b_k+b_k^\dagger).
$$

The bath is initially thermal and uncorrelated with the system:

$$
\rho(0)=\rho_S(0)\otimes\rho_B,\qquad
\rho_B=\frac{e^{-\beta H_B}}{\operatorname{Tr}_B e^{-\beta H_B}},\qquad
\beta=(k_BT)^{-1}.
$$

The numerical examples use

$$
H_S=\epsilon\sigma_z+\Delta\sigma_x,\qquad O_S=\sigma_z.
$$

There is **no factor of one-half** in this definition of the spin Hamiltonian. The text assigns energies $-\epsilon$ and $+\epsilon$ to $|0\rangle$ and $|1\rangle$, respectively. [P, Eqs. (1)-(7), pp. 2-3.]

**Implementation note:** in the ordered basis $(|0\rangle,|1\rangle)$, that labeling corresponds to

$$
H_S=\begin{pmatrix}-\epsilon&\Delta\\\Delta&\epsilon\end{pmatrix},\qquad
O_S=\begin{pmatrix}-1&0\\0&1\end{pmatrix}.
$$

Make the state ordering explicit instead of assuming a library's Pauli-matrix convention. Do not introduce an additional counterterm or another factor in $\Delta$ while claiming to implement the displayed model; none appears in Eqs. (1)-(5).

### 2.2 Spectral-density and correlation conventions

For positive frequencies,

$$
J(\omega)=\pi\sum_k|\lambda_k|^2\delta(\omega-\omega_k),
$$

$$
C(t)=\frac{1}{\pi}\int_0^\infty d\omega\,J(\omega)
\left[\coth\left(\frac{\beta\omega}{2}\right)\cos(\omega t)
-i\sin(\omega t)\right].
$$

For a two-sided frequency representation, use the **odd continuation**
$J(-\omega)=-J(\omega)$, giving

$$
S_\beta(\omega)=J(\omega)
\left[\coth\left(\frac{\beta\omega}{2}\right)+1\right],
\qquad
C(t)=\frac{1}{2\pi}\int_{-\infty}^{\infty}
S_\beta(\omega)e^{-i\omega t}\,d\omega.
$$

These factors of $\pi$, the Fourier sign, and the negative imaginary part of $C(t)$ must remain consistent throughout fitting and propagation. [P, Eqs. (8)-(11) and (15), pp. 3-4.]

The Brownian-oscillator benchmark uses

$$
J_{\mathrm{Br}}(\omega)=
\frac{2\Lambda\gamma\Omega^2\omega}
{(\omega^2-\Omega^2)^2+(\gamma\omega)^2},
\qquad
\Lambda=\frac{1}{\pi}\int_0^\infty\frac{J(\omega)}{\omega}\,d\omega.
$$

Here $\gamma$ is the physical Brownian damping parameter, $\Omega$ its natural frequency, and $\Lambda$ the reorganization energy. Do not confuse this scalar $\gamma$ with the fitted complex exponents $\gamma_\ell$ below. [P, Eqs. (12)-(13), p. 3.]

The structured spectral density $J_{\mathrm{Str}}$ is imported from the paper's Ref. [10]. Figure 1, p. 3, shows a broad low-frequency band and a narrow higher-frequency feature. An explicit analytic formula or tabulated numerical dataset for this spectrum is **not supplied in the PDF**.

### 2.3 Exponential decomposition: the solver's bath input

HEOM consumes a finite approximation

$$
C(t)\approx\widetilde C(t)=\sum_{\ell=1}^{K}
\eta_\ell e^{-\gamma_\ell t},\qquad
\gamma_\ell=\gamma_{\ell,r}+i\gamma_{\ell,i},
$$

where $K=N_\ell^{\max}$ in the paper. Both $\eta_\ell$ and $\gamma_\ell$ may be complex. [P, Eq. (14), p. 3.]

| Bath | Decomposition used in the paper | Reported detail |
|---|---|---|
| Brownian oscillator | Approximate $S_\beta(\omega)$ by rational functions using AAA; obtain exponential terms by residues. | Convergence generally reached at $K=3$. |
| Structured spectrum | Obtain a highly accurate reference $C(t)$ from a larger AAA expansion; compress the time signal with ESPRIT over $[0,t_{\max}]$. | The final pole count, pole table, sampling grid, and fit tolerances are not reported in the supplied PDF. |

[P, Sec. II, Eq. (16), p. 4.]

**Derived implementation convention:** with the paper's Fourier convention, a simple pole $z_\ell$ of a rational approximation to $S_\beta$ in the lower half-plane contributes, for $t>0$,

$$
\gamma_\ell=i z_\ell,\qquad
\eta_\ell=-i\operatorname{Res}_{z=z_\ell}S_\beta(z).
$$

This follows by closing the contour clockwise below the real axis. It assumes that the rational representation permits the contour argument; verify the resulting $\widetilde C(t)$ against an independently evaluated reference rather than relying on signs by inspection. This residue formula is an implementation derivation from Eq. (11), not an explicit algorithm given in the paper.

**Implementation notes:** validate decaying fitted components through $\operatorname{Re}\gamma_\ell>0$ for these damped-bath benchmarks, retain the sign of $\operatorname{Im}\gamma_\ell$, and compare both real and imaginary correlation components over the entire propagation window. A finite-window ESPRIT fit must not be assumed accurate beyond its fitted window. Since $C(t)$ is generally complex, do not enforce conjugate pairing merely to make it real. Right-action conjugation is already built into the hierarchy.

For a fixed spectral shape and temperature, scaling $\Lambda$ scales $J$, $C$, and the residues linearly; an existing pole set can therefore be amplitude-rescaled before checking its fit error. This is a derived consequence of Eqs. (9) and (12), not a reported fitting shortcut.

### 2.4 Units and scaling

The paper sets $e=\hbar=1$, while figures display energies in meV/eV and times in fs. [P, end of Sec. I, p. 2; figure captions.]

**Proposed implementation convention:** fix an energy scale $E_*$ and use dimensionless inputs throughout the hierarchy. For dimensionless $O_S$,

$$
\tau=\frac{E_*t_{\mathrm{phys}}}{\hbar},\qquad
\widehat H_S=H_S/E_*,\qquad
\widehat\eta_\ell=\eta_{\ell,\mathrm{energy}^2}/E_*^2.
$$

If exponents are stored as energies in the paper's $\hbar=1$ convention, use
$\widehat\gamma_\ell=\gamma_{\ell,\mathrm{energy}}/E_*$; if stored as physical rates, use $\widehat\gamma_\ell=\hbar\gamma_{\ell,\mathrm{rate}}/E_*$. The equations below then use dimensionless quantities, with hats suppressed.

Record the energy scale and ADO scaling in every run. The Gaussian transformation mixes auxiliary orders, so a change of numerical scale should not be treated as an innocuous finite-hierarchy rescaling without checking convergence. Do not mix a meV Hamiltonian with a femtosecond integrator while omitting the conversion by $\hbar$. The paper does not specify all internal unit/scaling choices of its implementation.

## 3. Standard HEOM: reference formulation

The original hierarchy has indices $j=(\ell,s)$, where $s=0$ denotes left action and $s=1$ right action. Define

$$
\mathcal L_S X=[H_S,X],\qquad
\mathcal L_O X=[O_S,X],\qquad
O_S^L X=O_SX,\qquad O_S^R X=XO_S.
$$

For the two sectors, the damping exponents are $\gamma_{(\ell,0)}=\gamma_\ell$ and $\gamma_{(\ell,1)}=\gamma_\ell^*$. The downward superoperator is

$$
\mathcal C_{(\ell,s)}X=
\eta_\ell O_SX\,\delta_{s0}-\eta_\ell^*XO_S\,\delta_{s1}.
$$

[P, Eqs. (19)-(21), p. 4.]

**Derived occupation-number transcription for a baseline implementation:** let $Z_{\mathbf u,\mathbf v}$ contain $u_\ell$ left and $v_\ell$ right occurrences of each pole. Then

$$
\begin{aligned}
\dot Z_{\mathbf u,\mathbf v}={}&-i[H_S,Z_{\mathbf u,\mathbf v}]
-\sum_\ell(u_\ell\gamma_\ell+v_\ell\gamma_\ell^*)Z_{\mathbf u,\mathbf v}\\
&-i\sum_\ell[O_S,
 Z_{\mathbf u+\mathbf e_\ell,\mathbf v}
 +Z_{\mathbf u,\mathbf v+\mathbf e_\ell}]\\
&-i\sum_\ell\left(
 u_\ell\eta_\ell O_SZ_{\mathbf u-\mathbf e_\ell,\mathbf v}
 -v_\ell\eta_\ell^*Z_{\mathbf u,\mathbf v-\mathbf e_\ell}O_S
\right).
\end{aligned}
$$

For factorized initial conditions, $Z_{\mathbf0,\mathbf0}(0)=\rho_S(0)$ and all other $Z$ vanish; here, unlike the transformed hierarchy, the physical density is simply $Z_{\mathbf0,\mathbf0}(t)$. The paper's standard formulation therefore has **$2K$ occupation modes**, not the $4K$ modes of the stabilized construction. This baseline is useful for weak-coupling and short-time cross-checks; do not silently substitute another HEOM convention. [P, Eqs. (17)-(21), p. 4; Eq. (61), p. 9.]

## 4. Stability mechanism and transformation

### 4.1 Why increasing the tier can fail

The paper isolates the minimal auxiliary-space generator

$$
\mathcal L=-\gamma a^\dagger a+\eta a^\dagger,\qquad \gamma>0.
$$

After projection onto occupations $0,\ldots,N$, its eigenvalues are $-\gamma n$, yet the propagator can have a large Euclidean norm. The paper derives the estimate

$$
\sup_{t\ge0}\|e^{t\mathcal L_N}\|_2
\gtrsim\frac{(|\eta|/\gamma)^N}{\sqrt{N!}}.
$$

The interpretation is competition between auxiliary damping and an unbalanced creation-like term. Strong couplings and long bath memory are particularly problematic. This estimate is not a claim of monotonic growth at every $N$, nor a universal stability threshold for the full HEOM. The full truncated hierarchy can also have genuine long-time spectral instabilities. [P, Eqs. (33)-(40), pp. 6-7; Appendix B, p. 14.]

For the minimal model, the Gaussian transformation

$$
e^{-S}=\exp[-(a+a^\dagger)^2/4]
$$

turns $\eta a^\dagger$ into $\eta(a^\dagger-a)/2$. Figure 2, p. 7, compares propagator norms for $N=50$, $\gamma=10^{-4}$, and $\eta=2\times10^{-4}(1+i)$, showing substantially smaller amplification after transformation. The plotted norm is a norm of an auxiliary representation, not a physical density-matrix observable. [P, Eqs. (41)-(46), p. 7.]

### 4.2 Full Gaussian and Bogoliubov transformations

Splitting each exponential into damped cosine and sine contributions produces four auxiliary modes per pole:

| Sector | Original ladder | Rotated ladder | Occupation |
|---|---|---|---|
| Left, cosine | $a_\ell$ | $A_\ell$ | $M_\ell$ |
| Left, sine | $b_\ell$ | $B_\ell$ | $N_\ell$ |
| Right, cosine | $c_\ell$ | $C_\ell$ | $O_\ell$ |
| Right, sine | $d_\ell$ | $D_\ell$ | $P_\ell$ |

The scalar occupation $O_\ell$ is unrelated to the system operator $O_S$. The transformed extended state is

$$
|\widetilde\rho_{\mathrm{ext}}\rangle\rangle
=\mathcal N_0e^{-S}|\rho_{\mathrm{ext}}\rangle\rangle,\qquad
\mathcal N_0=2^K,
$$

$$
e^{-S}=\prod_{\ell=1}^K\prod_{\alpha\in\{a,b,c,d\}}
\exp\left[-\frac{(\alpha_\ell+\alpha_\ell^\dagger)^2}{4}\right].
$$

For example,

$$
A_\ell=\frac{3a_\ell+a_\ell^\dagger}{2\sqrt2},\qquad
A_\ell^\dagger=\frac{3a_\ell^\dagger+a_\ell}{2\sqrt2},
$$

with analogous definitions for $B,C,D$. These preserve the canonical bosonic commutators. [P, Eqs. (22)-(30), pp. 5-6; Eqs. (47)-(57), p. 8; Appendix D, p. 15.]

**Implementation consequence:** derive the transformed infinite-space stencil first and then truncate it. Exponentiating a matrix $S$ inside an already truncated original hierarchy is not the operation analyzed in the paper; projection and the transformation do not commute. The explicit stencil below avoids constructing $e^{\pm S}$ entirely. [P, discussion following Eq. (46), p. 7; Sec. III B, p. 8.]

The coordinate variables in Appendices A, C, and D describe the **auxiliary bath-memory space**. They are not automatically a physical particle's position or a molecular reaction coordinate. The system basis remains a separate choice for $H_S$ and $O_S$.

## 5. Stabilized equations to implement

### 5.1 State convention, indexing, and truncation

Let

$$
\mathbf q=(M_1,N_1,O_1,P_1,\ldots,M_K,N_K,O_K,P_K),
\qquad R_{\mathbf q}=\widetilde\rho^{\mathbf M|\mathbf N}_{\mathbf O|\mathbf P}.
$$

This interleaved storage order is a proposed software convention; the four named occupation families are those of the paper. Write

$$
\mathbf q!\equiv\prod_{\ell=1}^K M_\ell!N_\ell!O_\ell!P_\ell!,
\qquad
|\widetilde\rho_{\mathrm{ext}}\rangle\rangle
=\sum_{\mathbf q}\frac{R_{\mathbf q}}{\sqrt{\mathbf q!}}
\otimes|\mathbf q\rangle.
$$

Thus $R_{\mathbf q}$ is a **factorial-unscaled ADO**, not the normalized Fock coefficient $F_{\mathbf q}=R_{\mathbf q}/\sqrt{\mathbf q!}$. The integer occupation coefficients below correspond to $R$. [P, Eq. (58), p. 9.]

Use the paper's **L-truncation**:

$$
\mathcal Q=\{\mathbf q\in\mathbb N_0^{4K}:|\mathbf q|_1\le n_{\max}\}.
$$

Define $R_{\mathbf q}=0$ whenever any occupation is negative or $|\mathbf q|_1>n_{\max}$. This is a total-tier cutoff, not a separate per-mode cutoff. No additional terminator is specified. [P, Eq. (67), p. 9.]

Let $\mathbf e_{M\ell},\mathbf e_{N\ell},\mathbf e_{O\ell},\mathbf e_{P\ell}$ increment the corresponding single occupations. All multi-index shifts below are simultaneous.

### 5.2 Complete component equation

The following is Eq. (60) rewritten with explicit matrix multiplication and additive occupation shifts. It also follows from the compact generator in Eq. (D12). Identity superoperators are suppressed.

$$
\begin{aligned}
\dot R_{\mathbf q}={}&-i[H_S,R_{\mathbf q}]
-\sum_\ell\gamma_{\ell,r}(M_\ell+N_\ell+O_\ell+P_\ell)R_{\mathbf q}\\[2pt]
&+\sum_\ell\gamma_{\ell,r}\sum_{X\in\{M,N,O,P\}}
R_{\mathbf q+2\mathbf e_{X\ell}}\\[2pt]
&-i\sqrt2\sum_\ell\left[O_S,
\sum_{X\in\{M,N,O,P\}}R_{\mathbf q+\mathbf e_{X\ell}}\right]\\[2pt]
&-\frac{i}{\sqrt2}\sum_\ell\Big[
\eta_\ell O_S\big(M_\ell R_{\mathbf q-\mathbf e_{M\ell}}
-R_{\mathbf q+\mathbf e_{M\ell}}\big)\\
&\hspace{90pt}-\eta_\ell^*
\big(O_\ell R_{\mathbf q-\mathbf e_{O\ell}}
-R_{\mathbf q+\mathbf e_{O\ell}}\big)O_S\Big]\\[2pt]
&-i\sum_\ell\gamma_{\ell,i}\Big[
M_\ell R_{\mathbf q-\mathbf e_{M\ell}+\mathbf e_{N\ell}}
+N_\ell R_{\mathbf q+\mathbf e_{M\ell}-\mathbf e_{N\ell}}
-2R_{\mathbf q+\mathbf e_{M\ell}+\mathbf e_{N\ell}}\Big]\\[2pt]
&+i\sum_\ell\gamma_{\ell,i}\Big[
O_\ell R_{\mathbf q-\mathbf e_{O\ell}+\mathbf e_{P\ell}}
+P_\ell R_{\mathbf q+\mathbf e_{O\ell}-\mathbf e_{P\ell}}
-2R_{\mathbf q+\mathbf e_{O\ell}+\mathbf e_{P\ell}}\Big].
\end{aligned}
$$

[P, Eq. (60), p. 9; Eq. (D12), p. 15. See the transcription cautions in Section 10 below.]

The equation needs same-tier, tier-minus-one, tier-plus-one, and tier-plus-two neighbors. In particular, retain both kinds of $+2$ shifts: same-mode shifts and mixed cosine/sine shifts. The right-coupling term multiplies $O_S$ **on the right**. Left cosine/sine mixing carries $-i\gamma_{\ell,i}$, while right mixing carries $+i\gamma_{\ell,i}$.

**Derived coefficient check:** for this ADO convention, projecting $A_\ell^\dagger$ produces $M_\ell R_{\mathbf q-\mathbf e_{M\ell}}$, projecting $A_\ell$ produces $R_{\mathbf q+\mathbf e_{M\ell}}$, and projecting $A_\ell^2$ produces $R_{\mathbf q+2\mathbf e_{M\ell}}$. Do not insert square roots of occupations into this stencil.

### 5.3 Independent operator form for testing

With all auxiliary and system tensor products understood, Eq. (D12) is

$$
\begin{aligned}
\widetilde{\mathcal L}={}&-i\mathcal L_S
-\sum_\ell\gamma_{\ell,r}\sum_{X\in\{A,B,C,D\}}
(X_\ell^\dagger X_\ell-X_\ell^2)\\
&-i\sqrt2\sum_\ell\mathcal L_O(A_\ell+B_\ell+C_\ell+D_\ell)\\
&-\frac{i}{\sqrt2}\sum_\ell\left[
\eta_\ell O_S^L(A_\ell^\dagger-A_\ell)
-\eta_\ell^*O_S^R(C_\ell^\dagger-C_\ell)\right]\\
&-i\sum_\ell\gamma_{\ell,i}
(A_\ell^\dagger B_\ell-2A_\ell B_\ell+A_\ell B_\ell^\dagger)\\
&+i\sum_\ell\gamma_{\ell,i}
(C_\ell^\dagger D_\ell-2C_\ell D_\ell+C_\ell D_\ell^\dagger).
\end{aligned}
$$

**Derived implementation caution:** an operator assembled with standard square-root ladder matrices acts on $F$, not on $R$. If $D_{\mathbf q\mathbf q}=\sqrt{\mathbf q!}$, then $R=DF$ and $L_R=DL_FD^{-1}$.

Also perform infinite-space algebra **before** projecting. For example, normal-order $A_\ell B_\ell^\dagger$ as $B_\ell^\dagger A_\ell$ before multiplying cutoff ladder matrices. The naive product $(PAP)(PB^\dagger P)$ can discard a temporary excursion above the cutoff and incorrectly remove a valid same-tier transfer. The direct component stencil does not have this ambiguity.

## 6. Initialization and physical reconstruction

For the factorized initial state,

$$
R_{\mathbf0}(0)=\rho_S(0),\qquad
R_{\mathbf q\ne\mathbf0}(0)=0.
$$

Although the Gaussian transformation alone mixes original Fock occupations, the combined Gaussian/Bogoliubov representation makes this initialization sparse. Do not initialize nonzero tiers by applying an additional squeezed-state formula. [P, Eqs. (61)-(62), Sec. III C, p. 9.]

At each requested output time, reconstruct

$$
\rho_S(t)=\sum_{\mathbf q\in\mathcal Q}w_{\mathbf q}R_{\mathbf q}(t),
$$

$$
w_{\mathbf q}=\begin{cases}
\displaystyle\prod_{j=1}^{4K}\frac{(q_j-1)!!}{q_j!},
&\text{every }q_j\text{ is even},\\
0,&\text{otherwise}.
\end{cases}
$$

Use $(-1)!!=1$, so $w_{\mathbf0}=1$. These weights already include the normalization needed by the inverse transformation; do not multiply the reconstruction by another $2^{-K}$. [P, Eqs. (63)-(66), p. 9; Appendix E, Eqs. (E1)-(E4), pp. 15-16.]

**Derived numerically convenient form:** for $q_j=2k_j$,

$$
\frac{(q_j-1)!!}{q_j!}=\frac{1}{2^{k_j}k_j!}.
$$

Precompute weights using this formula, a recurrence $w_{2k+2}=w_{2k}/(2k+2)$ for one occupation, or log-factorials. Example checks are

$$
w_{\mathbf0}=1,\qquad
w_{2\mathbf e_i}=\tfrac12,\qquad
w_{4\mathbf e_i}=\tfrac18,\qquad
w_{2\mathbf e_i+2\mathbf e_j}=\tfrac14\quad(i\ne j).
$$

A state with two odd occupations has even **total tier** but weight zero. Odd-occupation ADOs must still be propagated because they mediate dynamics. If storing normalized coefficients $F$ instead, reconstruction weights become $w_{\mathbf q}\sqrt{\mathbf q!}$.

For any system observable $A_S$, evaluate $\langle A_S\rangle=\operatorname{Tr}(A_S\rho_S)$. The benchmark quantity is $P_1(t)=\rho_{S,11}(t)$. Apply physical diagnostics to reconstructed $\rho_S$, not to the transformed root or individual ADOs.

## 7. Proposed implementation structure

This section is an engineering translation of Sections 5-6, not a description of the authors' source code.

### 7.1 Minimal interfaces and storage

Keep the following components separate:

```text
BathExpansion:
    eta[K]                  complex correlation amplitudes
    gamma[K]                complex decay exponents
    energy_scale            documented nondimensionalization
    fit_window, fit_error   provenance, when available

HierarchyTopology:
    occupations[num_ados, 4*K]
    tuple_to_index
    tier[num_ados]
    reconstruction_weights[num_ados]
    neighbor metadata       optional; trade memory against lookup cost

State:
    R[num_ados, d, d]        complex matrices in the raw-ADO convention

Core operations:
    build_topology(K, n_max)
    rhs(t, R, H_S, O_S, bath, topology)
    reconstruct_density(R, topology)
    diagnose(rho_S, R)
```

Validate that $H_S$, $O_S$, and $\rho_S(0)$ have matching square shapes, that the initial density is Hermitian and normalized, that pole arrays have equal lengths and finite values, and that $n_{\max}$ is a nonnegative integer. Keep complex arithmetic throughout; do not discard small imaginary parts before measuring numerical errors.

Generate weak compositions of each total tier into $4K$ occupations. This enumerates the retained simplex directly rather than constructing and filtering a $(n_{\max}+1)^{4K}$ Cartesian box. Include the root once, order tuples deterministically, and map absent neighbors to a read-only zero matrix.

For a simultaneous shift such as $\mathbf q-\mathbf e_M+\mathbf e_N$, apply **all** changes before validating the target. Do not sequentially test an intermediate increment against the tier boundary.

The number of retained ADOs, derived from this truncation, is

$$
N_{\mathrm{ADO}}=\binom{4K+n_{\max}}{n_{\max}}.
$$

For $K=3$ and $n_{\max}=14$, the unreduced hierarchy contains **9,657,700 ADOs**. A single array of $2\times2$ complex128 matrices occupies approximately **589.5 MiB**, before Runge-Kutta stages, derivatives, topology, or output. This is a combinatorial estimate, not a memory measurement from the paper; pruning inactive modes or exploiting additional structure can change it. Full neighbor tables may themselves be expensive. Do not use a dense global Liouvillian for these production sizes.

The standard Eq. (19) representation has only $2K$ occupation modes. Consequently, equal numerical $n_{\max}$ does not imply equal storage, work, or equivalent finite approximations. The paper explicitly warns that its standard and transformed truncations are not term-by-term identical. [P, discussion above Sec. IV, p. 10.]

### 7.2 Reference RHS pseudocode

This is language-independent pseudocode. `at(q, shifts)` returns the stored ADO at the simultaneously shifted tuple, or a zero matrix if the final tuple is outside the retained set. `@` means matrix multiplication; `comm(A,B) = A @ B - B @ A`.

```text
function rhs(R):
    dR = zeros_like(R)
    for row, q in retained_occupations:
        dR[row] = -i * comm(H_S, R[row])

        for ell in 0 .. K-1:
            m, n, o, p = 4*ell + (0, 1, 2, 3)
            gr = real(gamma[ell])
            gi = imag(gamma[ell])
            eta_l = eta[ell]

            dR[row] -= gr * (q[m]+q[n]+q[o]+q[p]) * R[row]

            for j in (m, n, o, p):
                dR[row] += gr * at(q, {j: +2})
                dR[row] -= i*sqrt(2) * comm(O_S, at(q, {j: +1}))

            left = q[m]*at(q, {m: -1}) - at(q, {m: +1})
            right = q[o]*at(q, {o: -1}) - at(q, {o: +1})
            dR[row] -= i/sqrt(2) * (
                eta_l * (O_S @ left)
                - conjugate(eta_l) * (right @ O_S)
            )

            dR[row] -= i*gi * (
                q[m]*at(q, {m: -1, n: +1})
                + q[n]*at(q, {m: +1, n: -1})
                - 2*at(q, {m: +1, n: +1})
            )
            dR[row] += i*gi * (
                q[o]*at(q, {o: -1, p: +1})
                + q[p]*at(q, {o: +1, p: -1})
                - 2*at(q, {o: +1, p: +1})
            )
    return dR
```

A direct reference implementation may use tuple lookups. Optimize only after validating it: precompute selected neighbors, vectorize system-matrix actions, parallelize independent target-ADO updates, or use a structured/tensor representation. These optimizations are not prescribed or benchmarked in the paper.

### 7.3 Propagation and output

The authors use an **adaptive-timestep fourth-order Runge-Kutta** method. The PDF does not specify the full adaptive tableau/controller, tolerances, or timestep limits. [P, Sec. III C, step 2, p. 9.]

Expose an integrator-independent complex-valued RHS. First establish accuracy with conservative integration settings, then converge the hierarchy and bath fit independently. Store physical outputs and selected diagnostics rather than every ADO at every output time. Checkpoints intended for restart must contain the full hierarchy state and its convention metadata, not only $\rho_S$.

Do not enforce unit trace on $R_{\mathbf0}$, impose Hermiticity on every ADO, or silently clip negative populations. Such operations are not part of the paper's algorithm and can hide reconstruction, truncation, or propagation errors.

## 8. Reported numerical results and reproduction targets

### 8.1 Common spin-boson setup

All population benchmarks are unbiased, with

$$
\epsilon=0,\quad\Delta=10\ \mathrm{meV},\quad
k_BT=25.85\ \mathrm{meV},\quad
\rho_S(0)=|1\rangle\langle1|.
$$

The paper expects equal stationary populations, $P_0=P_1=0.5$, for these cases. This is not an equilibrium prescription for arbitrary $H_S$, bias, or coupling operators. [P, Fig. 3 caption and Eq. (68), p. 10.]

| Target | Bath/coupling parameters | Tier information | Displayed interval |
|---|---|---|---|
| Fig. 3, p. 10 | Brownian: $\gamma=4.9$ meV, $\Omega=13.89$ meV; $\Lambda=20,40,60$ meV | Converged transformed tier varies with coupling; largest is $n_{\max}=14$ for $\Lambda=60$ meV. | 0-2000 fs |
| Fig. 4, p. 10 | Same Brownian bath; $\Lambda=60$ meV | $n_{\max}=10,12,14$ | 0-2000 fs |
| Fig. 5, p. 11 | Structured bath; $\Lambda=60,80,100$ meV | $n_{\max}=8$ for all three couplings | 0-3000 fs |
| Fig. 6, p. 12 | Structured bath; $\Lambda=100$ meV | $n_{\max}=4,6,8$ | 0-3000 fs |

The structured bath inherits the common spin/temperature parameters, not the Brownian spectral function. The original spectrum must be obtained independently before claiming to reproduce Figures 5-6.

### 8.2 What the figures establish

For the Brownian bath, Figures 3-4 show transformed populations relaxing smoothly toward $0.5$, while the standard hierarchy eventually leaves the physical population interval and diverges. Stronger coupling moves the standard-hierarchy divergence earlier. In the tier sweep, increasing the standard tier does not systematically stabilize it; the transformed results instead show a convergence pattern. [P, Sec. IV A, pp. 10-11.]

For the structured bath, Figures 5-6 show stable transformed dynamics over the plotted interval at all displayed couplings. The standard hierarchy again becomes unstable, with timing sensitive to both coupling and tier. The transformed curves in Figure 6 nearly coincide for tiers 4, 6, and 8. The paper attributes the more demanding dynamics to combined broad, short-memory and narrow, long-memory spectral contributions. [P, Sec. IV B, pp. 11-12.]

These are stability/convergence demonstrations for the reported model and windows. The PDF does not provide wall-clock speedups, a general positivity guarantee, raw trajectories for numerical regression, or proof of stability for arbitrary coupling and propagation time.

## 9. Proposed validation and convergence plan

These tests are implementation recommendations derived from the equations; they are not an additional benchmark suite supplied by the authors.

### 9.1 Algebra and topology tests

Check that tuple counts equal $\binom{4K+n_{\max}}{n_{\max}}$, indices are unique, and absent/negative neighbors contribute zero. Test every shift type, including mixed same-tier transfers **at the truncation boundary**.

For small hierarchies, assemble Eq. (D12) independently in the normalized Fock basis and verify $L_R=DL_FD^{-1}$ against the component RHS on random complex matrices. Use complex residues, both signs of $\operatorname{Im}\gamma_\ell$, and a non-diagonal Hermitian $O_S$ so left/right multiplication mistakes cannot cancel accidentally. Apply the normal-ordering caution from Section 5.3.

Test reconstruction weights against the examples in Section 6, including a tuple with even total tier but two odd entries. Check that the initialized state reconstructs to exactly the specified $\rho_S(0)$.

### 9.2 Dynamical checks

With all $\eta_\ell=0$ and non-root ADOs initially zero, verify

$$
\rho_S(t)=e^{-iH_St}\rho_S(0)e^{+iH_St}
$$

and that non-root ADOs remain zero. For the paper's unbiased two-state initial condition, this gives $P_1(t)=\cos^2(\Delta t)$ in natural units. This is a derived isolated-system check.

For any factorized initial state with the stated zero-mean thermal bath, verify the initial derivative

$$
\left.\frac{d\rho_S}{dt}\right|_{0}=-i[H_S,\rho_S(0)].
$$

Next compare independently converged standard and transformed hierarchies at weak coupling and short times, where the standard calculation is still well behaved. At a common finite tier, exact equality is not the correct requirement because their truncations differ.

### 9.3 Diagnostics and independent convergence axes

For reconstructed $\rho_S$, record trace error, Hermiticity error, populations, and the smallest eigenvalue of its Hermitian part. Also record ADO norms by tier, the highest-tier contribution to reconstruction, failed integration steps, and non-finite values. Individual ADOs are not physical density matrices.

Refine the bath fit/pole count, hierarchy tier, and integration tolerances separately. A stable trajectory is not automatically accurate, and a small bath-fit residual does not prove hierarchy convergence. Stability improvement does not remove the need to converge the system basis in applications larger than the spin benchmark.

Do not substitute post-step normalization, Hermiticity projection, or population clipping for an error diagnosis. Acceptance thresholds must be chosen for the scientific observable and recorded; the paper does not provide a universal tolerance.

### 9.4 Checks performed while preparing this transcription

An independent local algebra check compared the component RHS with a normally ordered Eq. (D12) matrix acting on factorial-normalized coefficients. Random complex inputs were tested at $(K,n_{\max})=(1,2),(1,4),(2,3)$; relative discrepancies were below $3\times10^{-16}$. A zero-coupling propagation check agreed with direct unitary evolution to about $2.3\times10^{-16}$ in Frobenius norm.

These checks validate the transcription and normalization convention, not the physical bath fitting or the paper's convergence claims. **The paper's Brownian and structured benchmark curves were not reproduced in preparing this document.**

## 10. Source ambiguities and missing reproducibility details

### 10.1 Displayed-equation cautions

The following are explicit transcription decisions, not silent corrections to the source.

**Factorial normalization:** the projection notation in Eq. (59), and similarly in some earlier projection displays, omits the factorial multiplier implied by the expansion in Eq. (58). This document defines $R$ by Eq. (58), for which

$$
R_{\mathbf q}=\sqrt{\mathbf q!}\,
\langle\mathbf q|\widetilde\rho_{\mathrm{ext}}\rangle\rangle.
$$

That convention produces the integer coefficients printed in Eq. (60). A bare projection instead defines $F$, requiring square-root coefficients and different reconstruction weights. Do not combine the two conventions.

**Trailing symbols in Eq. (60):** the printed cosine/sine-mixing lines contain dangling plus signs before their closing brackets. No further terms are specified there. This document transcribes the three mixing terms per sector that are explicitly present in Eq. (D12); it does not invent an omitted term.

**Auxiliary coordinate representation:** the Fock/coordinate material in Appendices A and D uses different Hermite conventions for the original and rotated modes. A component-ADO implementation does not need to implement these integrals. Deriving a separate coordinate-space solver would require additional care and is outside this specification.

### 10.2 What the PDF does not supply

The paper supplies enough information to implement the transformed hierarchy, its initialization, and its reconstruction. It does **not** supply everything needed for a byte-for-byte or curve-by-curve reproduction. Missing items include the numerical structured spectral density, fitted pole/residue tables, the structured fit's final pole count, AAA/ESPRIT sampling and tolerances, the precise fitting window, full integration settings, all internal ADO scaling choices, and machine-readable reference trajectories.

The analytic Brownian spectrum and its physical parameters are available, so an independently fitted and converged Brownian implementation is feasible. Treat $K=3$ as the reported experience for that case, not as a guarantee for every temperature, fit method, error target, or simulation window. For the structured case, retain an explicit external-data dependency rather than inventing a substitute spectrum and labeling it a reproduction.

## 11. Minimal implementation sequence

**First milestone:** implement topology, the raw-ADO RHS, factorized initialization, and even-occupation reconstruction for a small supplied pole set. Pass normalization, shift, left/right-action, and isolated-system tests before adding a fitting algorithm.

**Second milestone:** implement the standard Eq. (19) baseline, build a validated Brownian correlation expansion, and compare converged weak-coupling dynamics. Then reproduce the qualitative stability and tier trends of Figures 3-4 with documented numerical choices.

**Third milestone:** acquire the original structured spectrum, validate AAA-to-ESPRIT compression over the intended window, and address Figures 5-6. Optimize memory, parallelism, or tensor structure only after preserving the reference implementation as a regression oracle.

**Completion criterion:** validated physical-density reconstruction and demonstrated convergence of the target observables with respect to bath representation, tier, and integration error. Merely obtaining a bounded transformed root ADO does not establish a working physical simulation.

## 12. Citation and source navigation

```bibtex
@misc{gatto2026stabilizedheom,
  title = {Improving the Stability of the Hierarchical Equations of Motion
           for Open Quantum Systems with Strong Coupling to Structured
           Bosonic Baths},
  author = {Gatto, Salvatore and Rudge, Samuel L. and Hou, Bokang
            and Rabani, Eran and Thoss, Michael},
  year = {2026},
  eprint = {2609.35484},
  archivePrefix = {arXiv},
  primaryClass = {quant-ph},
  note = {Version 1, 28 September 2026}
}
```

| Implementation concern | Location in [P] |
|---|---|
| Model, spectral-density convention, exponential bath expansion | Sec. II, Eqs. (1)-(16), pp. 2-4 |
| Standard HEOM | Sec. III, Eqs. (17)-(21), p. 4 |
| Auxiliary-space construction and instability mechanism | Sec. III A, Eqs. (22)-(46), pp. 5-7; Appendix B, p. 14 |
| Gaussian transformation and rotated ladders | Sec. III B, Eqs. (47)-(58), pp. 8-9 |
| Production component equation | Eq. (60), p. 9 |
| Initial state, reconstruction, total-tier cutoff | Sec. III C, Eqs. (61)-(67), p. 9 |
| Independent compact generator | Appendix D, Eq. (D12), p. 15 |
| Reconstruction derivation and parity condition | Appendix E, Eqs. (E1)-(E4), pp. 15-16 |
| Numerical reproduction targets | Sec. IV, Figs. 3-6, pp. 10-12 |

The paper cites Ref. [64] for the frequency-domain decomposition framework, Refs. [49]-[51] for AAA/ESPRIT-related methods, and Ref. [10] for the structured spectral density. Those works are referenced here only as dependencies identified by the supplied paper, not as independently reviewed sources.
