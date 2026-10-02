# Cabrera et al. (2015): FFT propagation of the Wigner function

## Citation and purpose

**Paper [P]:** Renan Cabrera, Denys I. Bondar, Kurt Jacobs, and Herschel A. Rabitz, "Efficient method to generate time evolution of the Wigner function for open quantum systems," *Physical Review A* **92**, 042122 (2015). DOI: `10.1103/PhysRevA.92.042122`.

**Source:** The supplied ten-page published article. References below use its equation numbers and PDF page numbers. The Python supplementary material mentioned in Ref. [70] was not supplied or inspected.

**Purpose:** Repository context sufficient to implement a one-dimensional spectral Wigner solver, add the paper's decoherence and Caldeira-Leggett terms, and plan its two-particle extension.

**Provenance convention:** Statements labelled **Paper** summarize [P]. Statements labelled **Implementation derivation** or **Implementation choice** are explicit translations, checks, or proposed engineering decisions, not additional claims made by the authors. The code is newly written for this note, not a transcription of their supplementary code. No outside literature is used.

**Important:** The printed article contains several internally inconsistent formulas and incompletely specified benchmark inputs. They are flagged where relevant, rather than silently repaired. In particular, do not copy Eqs. (58), (59), (61), or the Morse benchmark directly into code without reading the associated notes.

## 1. Main contribution and supported scope

**Paper.** The method combines Hilbert phase space with spectral split-operator propagation. Rather than discretizing the infinite derivative expansion of the Moyal equation, it switches between Fourier-related representations in which individual pieces of the generator become pointwise multipliers. A potential is evaluated at displaced coordinates, so its quantum contribution does not require explicit high-order derivatives. The per-step cost is $O(N\log N)$ for $N$ stored phase-space samples. [P, abstract; Secs. II A and III; Eqs. (22)-(29), (50)-(60).]

The worked models are a particle with

$$
\hat H=\frac{\hat P^2}{2m}+V(\hat X),
$$

position-basis decoherence, and the time-independent Caldeira-Leggett approximation. The same formalism describes classical Liouville, Koopman-von Neumann, and diffusion-only Fokker-Planck dynamics. The paper also demonstrates a coupled two-particle system with a bath acting directly on only one particle. [P, Secs. II-V.]

The computational advantage is efficient propagation of a **full grid**, not removal of dimensional scaling: for $d$ spatial coordinates, $W$ has $2d$ phase-space coordinates. The potential-difference representation retains the quantum generator; the numerical method still has finite-domain, finite-grid, and time-splitting errors.

## 2. Formalism: representations and conventions

### 2.1 Coordinates and units

Use $\hat X,\hat P$ for physical quantum operators and $x,p,\lambda,\theta$ for Hilbert-phase-space coordinates. This distinguishes symbols which the paper differentiates partly through typography.

| Symbol | Meaning | Units or normalization |
|---|---|---|
| $x,p$ | Position and momentum coordinates | Length; momentum |
| $\lambda$ | Fourier variable conjugate to $x$ | Inverse length |
| $\theta$ | Fourier variable conjugate to $p$ | Inverse momentum |
| $\hbar\theta$ | Position separation in a density-matrix element | Length |
| $W(x,p)$ | Wigner function | $\int W\,dx\,dp=1$ |
| $B(x,\theta)$ | Double-configuration-space / Blokhintsev representation | Defined below |
| $D$ | Momentum diffusion coefficient | Momentum squared per time |
| $\gamma$ | Paper's Caldeira-Leggett damping parameter | Inverse time; momentum damping rate is $2\gamma$ |

The equations below retain $\hbar$. The numerical examples in the paper use atomic units. [P, Eqs. (4)-(9), (15)-(21), (39)-(41); Sec. IV.]

### 2.2 Wigner and Blokhintsev representations

**Paper, Eqs. (5)-(9):**

$$
B(x,\theta)=
\left\langle x-\frac{\hbar\theta}{2}\middle|\rho\middle|
 x+\frac{\hbar\theta}{2}\right\rangle,
$$

$$
B(x,\theta)=\int_{-\infty}^{\infty}W(x,p)e^{-ip\theta}\,dp,
\qquad
W(x,p)=\frac{1}{2\pi}\int_{-\infty}^{\infty}
B(x,\theta)e^{+ip\theta}\,d\theta.
$$

There is **no additional $1/\hbar$ in this Fourier pair**: the density-matrix displacement already contains $\hbar\theta$. If instead the integration variable is the physical separation $s=\hbar\theta$, the inverse prefactor is $1/(2\pi\hbar)$.

For a Hermitian density operator, $W$ is real and $B(x,-\theta)=B(x,\theta)^*$. Negative $W$ values are allowed and must not be clipped. Useful direct consequences are

$$
\int W(x,p)\,dp=B(x,0)=\langle x|\rho|x\rangle,
\qquad
\int W(x,p)\,dx=\langle p|\rho|p\rangle.
$$

For a normalized pure state, initialization is

$$
B(x,\theta)=\psi\!\left(x-\frac{\hbar\theta}{2}\right)
\psi^*\!\left(x+\frac{\hbar\theta}{2}\right).
$$

The signs, conjugation, and transform direction must be kept together. A wavefunction with phase $e^{+ip_0x/\hbar}$ must produce a Wigner distribution centered at **positive** $p_0$. [P, Eqs. (5), (7), (8); initialization discussed on p. 5.]

### 2.3 Other representations and the operator algebra

**Paper.** The additional representations are

$$
Z(\lambda,p)=\int e^{-i\lambda x}W(x,p)\,dx,
\qquad
A(\lambda,\theta)=\int e^{-i\lambda x}B(x,\theta)\,dx.
$$

The diagram in Eq. (29), on PDF p. 3, connects $W$, $B$, $Z$, and $A$ by partial Fourier transforms. In implementation, $W\leftrightarrow Z$ transforms the position axis, while $W\leftrightarrow B$ transforms the momentum axis. [P, Eqs. (27)-(29).]

The extended algebra includes

$$
[\hat x,\hat p]=0,\quad [\hat x,\hat\lambda]=i,
\quad [\hat p,\hat\theta]=i,\quad [\hat\lambda,\hat\theta]=0.
$$

In the $x,p$ representation,

$$
\hat\lambda=-i\partial_x,\qquad \hat\theta=-i\partial_p.
$$

The physical and mirror variables are represented by

$$
\hat X_L=\hat x-\frac{\hbar}{2}\hat\theta,
\quad \hat P_L=\hat p+\frac{\hbar}{2}\hat\lambda,
\qquad
\hat X_R=\hat x+\frac{\hbar}{2}\hat\theta,
\quad \hat P_R=\hat p-\frac{\hbar}{2}\hat\lambda.
$$

These represent left and right density-operator action, leading to

$$
i\hbar\,\partial_t|\rho\rangle=
[H(\hat X_L,\hat P_L)-H(\hat X_R,\hat P_R)]|\rho\rangle.
$$

[P, Eqs. (15)-(22).] The implementation does not need to construct matrices for these operators.

### 2.4 Governing equation for a separable Hamiltonian

Define the **scalar** potential difference in $x,\theta$ space:

$$
\Delta V(x,\theta)=V^- -V^+
=V\!\left(x-\frac{\hbar\theta}{2}\right)
-V\!\left(x+\frac{\hbar\theta}{2}\right).
$$

Then

$$
\partial_t B=-\frac{i}{m}\partial_x\partial_\theta B
-\frac{i}{\hbar}\Delta V B.
$$

Equivalently, the abstract generator is

$$
\partial_t|\rho\rangle=
-i\left(\frac{\hat p\hat\lambda}{m}
+\frac{V^- -V^+}{\hbar}\right)|\rho\rangle.
$$

The kinetic term is diagonal in $\lambda,p$ space; the potential term is diagonal in $x,\theta$ space. This is the core implementation mechanism. [P, Eqs. (22)-(25), (50).]

## 3. Open-system generators

### 3.1 General Lindblad structure

**Paper.** A dissipative contribution has the form

$$
\mathcal L[\rho]=A\rho A^\dagger
-\tfrac12 A^\dagger A\rho-\tfrac12\rho A^\dagger A.
$$

Left and right actions can be expressed using the physical and mirror operators above. [P, Eqs. (30)-(35).]

**Implementation limitation.** This formal mapping does not make every arbitrary jump operator a pointwise multiplier. Noncommuting products require the actual left/right operator ordering to be retained. The explicit numerical recipes developed in the paper are for the simpler generators below; a general Lindblad engine is a separate implementation task.

### 3.2 Position-basis decoherence / momentum diffusion

**Implementation convention anchored to the paper's evolution equations.** Define the actual contribution to $\dot\rho$ as

$$
\mathcal L_D\rho=-\frac{D}{\hbar^2}[\hat X,[\hat X,\rho]].
$$

It is equivalent to a position jump operator $A=\sqrt{2D}\,\hat X/\hbar$. In the useful representations,

$$
\left.\partial_tB\right|_D=-D\theta^2B,
\qquad
\left.\partial_tW\right|_D=D\partial_p^2W.
$$

Hence the full $B$ equation is

$$
\partial_tB=
\left[-\frac{i}{m}\partial_x\partial_\theta
-\frac{i}{\hbar}\Delta V-D\theta^2\right]B.
$$

The diffusion multiplier over a duration $h$ is exactly $e^{-hD\theta^2}$. Potential evolution and this diffusion term commute because both are multipliers in $x,\theta$. [P, Eqs. (38), (60).]

**Printed-notation warning.** Eqs. (36)-(37) put $i\hbar\mathcal L$ on their left-hand sides while displaying the real diffusion-generator expression on the right. Literal division by $i\hbar$ would not yield the damping in Eq. (38) or Eq. (60). This note uses the time-derivative convention above, fixed by those latter equations. Do not insert an extra factor of $i$ or $\hbar$ into the diffusion exponential.

### 3.3 Caldeira-Leggett damping and diffusion

Using the same actual-time-derivative convention, the model implemented by Eqs. (40)-(41) is

$$
\partial_t\rho=-\frac{i}{\hbar}[\hat H,\rho]
-\frac{i\gamma}{\hbar}[\hat X,\{\hat P,\rho\}]
-\frac{2m\gamma k_BT}{\hbar^2}[\hat X,[\hat X,\rho]].
$$

Thus

$$
D=2m\gamma k_BT,
$$

and the additional terms are

$$
\left.\partial_tW\right|_{\rm bath}
=2\gamma\partial_p(pW)+D\partial_p^2W,
$$

$$
\left.\partial_tB\right|_{\rm bath}
=-2\gamma\theta\partial_\theta B-D\theta^2B.
$$

[P, Eqs. (39)-(41), (63)-(65).] Eq. (39) has a similar $i\hbar\hat{\mathcal D}$ notation issue; the convention here follows Eqs. (40)-(41) and the numerical update.

**Paper limitation.** The authors explicitly note that this Caldeira-Leggett approximation is not Lindblad and does not generally guarantee density-matrix positivity. They use it as an approximate damping/thermalization model. Stability of a numerical calculation does not remove this physical limitation. [P, Sec. II B, p. 4.]

**Implementation consequence.** Do not confuse $\gamma$ with the decay rate of mean momentum: with the force omitted, $d\langle p\rangle/dt=-2\gamma\langle p\rangle$. Also, diffusion alone is not thermalization; it broadens momentum without the opposing friction term.

## 4. Classical dynamics and regularization

**Paper.** Taking $\hbar\rightarrow0$ gives

$$
\frac{\Delta V}{\hbar}\longrightarrow-\theta V'(x),
$$

and therefore the Liouville/Koopman-von Neumann transport equation

$$
\partial_t f=-\frac{p}{m}\partial_xf+V'(x)\partial_pf.
$$

For a classical probability density, $f=\rho_c\geq0$ and $\int\rho_c\,dx\,dp=1$. For a Koopman-von Neumann amplitude, $f=\Psi$ and the physical density is $\rho_c=|\Psi|^2$, normalized by $\int|\Psi|^2\,dx\,dp=1$. Without diffusion, both obey the same first-order transport equation. [P, Eqs. (42)-(44).]

Adding the paper's classical diffusion term gives

$$
\partial_t\rho_c=-\frac{p}{m}\partial_x\rho_c
+V'(x)\partial_p\rho_c+D\partial_p^2\rho_c.
$$

This is the Fokker-Planck model in Eq. (49). Its classical potential-plus-diffusion multiplier is

$$
Q^{\rm cl}_h(x,\theta)
=\exp\!\left[+ih\theta V'(x)-hD\theta^2\right].
$$

**Printed-equation warning.** Eq. (61) prints a factor $\exp[-i\,dt\,V'(x)-\delta D\theta^2]$, without the $\theta$ in the force phase. Both the missing factor and the sign conflict with the classical limit of Eq. (60), given the paper's Fourier convention. The expression above is an explicit derivation from Eqs. (42)-(43), (49), and (60), not a literal transcription of Eq. (61).

**Paper.** Classical evolution develops fine phase-space structure, termed velocity filamentation, that eventually becomes unresolved on a fixed grid. The authors suppress unresolved structure with a weak factor $e^{-\delta D\theta^2}$. Their Fig. 3 shows that this reduces spurious drift in the negative volume of the signed Koopman-von Neumann amplitude. [P, Sec. III A, Eq. (61); Sec. IV, Fig. 3.]

**Implementation choices.** Keep numerical filtering separate from physical diffusion. A filter $e^{-\epsilon\theta^2}$ has per-step strength $\epsilon$; a diffusion *rate* $D_{\rm reg}$ requires $\epsilon=D_{\rm reg}dt$. The printed Eq. (61) contains no $dt$ multiplying $\delta D$, so do not silently identify that symbol with a timestep-independent rate.

Do not conflate a signed amplitude with a probability density. A Wigner function with negative regions can be rescaled into an initial real $\Psi$, but then its associated classical density is $|\Psi|^2$, not $W$. Likewise, after adding diffusion to $\Psi$, squaring the evolved amplitude does not generally produce the solution obtained by applying the same diffusion equation to $\rho_c$; the equivalence of the two transport equations relies on their first-order derivative structure.

## 5. Spectral split-operator propagation

### 5.1 Elementary substeps

For duration $h$, define two operations on a state stored in $x,p$ space:

$$
\mathsf K_hW=\mathcal F^{-1}_{x\to\lambda}
\left[e^{-ihp\lambda/m}\,\mathcal F_{x\to\lambda}W\right],
$$

$$
\mathsf Q_hW=\mathcal F^{-1}_{p\to\theta}
\left[e^{-ih\Delta V/\hbar-hD\theta^2}
\,\mathcal F_{p\to\theta}W\right].
$$

Here each forward transform has a negative exponential and its inverse a positive exponential. The inverse-transform notation means returning from the indicated Fourier variable to the original coordinate.

**Paper's first-order step:**

$$
W_{n+1}=\mathsf K_{dt}\mathsf Q_{dt}W_n.
$$

Apply the rightmost operation first: potential/diffusion, then kinetic transport. Eq. (53) gives this for unitary dynamics and Eq. (60) adds diffusion. The local error is $O(dt^2)$; the accumulated error over a fixed duration is generally $O(dt)$.

**Paper's second-order unitary step, Eq. (52):**

$$
W_{n+1}=\mathsf K_{dt/2}\mathsf Q_{dt}\mathsf K_{dt/2}W_n.
$$

The local error is $O(dt^3)$ and fixed-duration global error is generally $O(dt^2)$. Including diffusion inside $\mathsf Q$ uses the same second-order splitting, since the potential and diffusion terms commute. That combined second-order assembly is an implementation derivation from Eqs. (52) and (60).

### 5.2 Discrete grid contract

**Implementation choice:** Use complex arrays with shape `(n_p, n_x)`. Axis 0 is momentum or $\theta$; axis 1 is position or $\lambda$. Use an even number of points in each direction, as in the paper.

The centered physical grids are

$$
x_j=-L_x+j\,dx,\quad dx=\frac{2L_x}{N_x},\qquad
p_k=-L_p+k\,dp,\quad dp=\frac{2L_p}{N_p}.
$$

The positive endpoint is excluded. Store coordinate and state arrays in FFT order, with zero first and negative coordinates in the second half. Use `ifftshift` once to convert centered inputs to this order, and `fftshift` for display. Avoid shifting arrays inside each substep. For even sizes these shifts coincide; using their directional names makes the convention explicit. [P, Eqs. (54)-(57), p. 5.]

**Implementation derivation: angular Fourier grids**

```python
lam = 2.0 * np.pi * np.fft.fftfreq(n_x, d=dx)
theta = 2.0 * np.pi * np.fft.fftfreq(n_p, d=dp)
```

Consequently,

$$
d\lambda=\frac{2\pi}{N_xdx}=\frac{\pi}{L_x},
\qquad
d\theta=\frac{2\pi}{N_pdp}=\frac{\pi}{L_p}.
$$

**Printed-grid warning.** Eqs. (58)-(59) instead state $d\theta=2\pi/L_p$ and $d\lambda=2\pi/L_x$. With the half-width definitions in Eqs. (54)-(56), these are twice the spacings implied by the discrete Fourier pair. The code below uses the `fftfreq` construction, not those two printed expressions. This is an identified internal inconsistency, not a verified author erratum.

### 5.3 Physical transform normalization

For FFT-ordered arrays, the continuum-normalized discrete transforms are

```python
B = dp * np.fft.fft(W, axis=0)
W = np.fft.ifft(B, axis=0) / dp
Z = dx * np.fft.fft(W, axis=1)
W = np.fft.ifft(Z, axis=1) / dx
```

The inverse $B\to W$ factor follows from $N_p\,d\theta/(2\pi)=1/dp$.

Inside a split step, `fft` and `ifft` scaling cancels around a multiplier, so the explicit `dx` or `dp` factors are unnecessary. Those unscaled intermediate arrays must not be treated as physical $B$ or $Z$ when inspecting values or initializing states. In particular, applying an inverse FFT to an analytically sampled $B$ requires division by `dp`.

### 5.4 Potential domain and boundaries

**Implementation consequences.** The solver needs a callable $V$ at $x\pm\hbar\theta/2$, not only values at the physical position grid. The largest displacement is approximately $\hbar\pi/(2dp)$. If a potential is tabulated, interpolation and out-of-range behavior must be explicitly specified; periodic extrapolation of the physical potential is not implied by the paper.

FFTs make the numerical state representation periodic on the chosen windows. For a nonperiodic physical problem, choose windows large enough that wraparound does not contaminate the occupied region. Check position and momentum boundaries throughout propagation, not only at initialization. Refining resolution and enlarging a window are distinct convergence operations.

Large shifted coordinates can also cause potential overflow or cancellation in $V^- -V^+$. Check finite values, and derive a stable analytic difference for a chosen potential when necessary. Nyquist components and unresolved tails can produce non-negligible imaginary residues; retaining complex arrays makes that error visible rather than hiding it with `.real`.

## 6. Friction update and complete Caldeira-Leggett step

**Paper, Eqs. (62)-(65).** Separate friction from the already treated potential/diffusion part. With

$$
\mathsf C W=2\gamma\partial_p(pW),
$$

the paper approximates $e^{h\mathsf C}$ by

$$
W^{(1)}=W+\frac{h}{2}\mathsf C W,
\qquad
W_{\rm new}=W+h\mathsf C W^{(1)}.
$$

For this linear generator, the result is $(1+h\mathsf C+h^2\mathsf C^2/2)W$, with local error $O(h^3)$. This is an explicit second-order Runge-Kutta/Taylor update, not an exact exponential or an implicit method.

Evaluate the derivative spectrally:

$$
\mathsf C W=2\gamma\,\operatorname{IFFT}_p
\left[i\theta\,\operatorname{FFT}_p(pW)\right].
$$

Multiply by $p$ **before** differentiating. Replacing $\partial_p(pW)$ by $p\partial_pW$ drops the $W$ term and breaks the intended damping generator.

**Implementation choice: overall second-order assembly.** Let $\widetilde{\mathsf F}_h$ denote the friction update above and $\mathsf S_{dt}=\mathsf K_{dt/2}\mathsf Q_{dt}\mathsf K_{dt/2}$. Use

$$
W_{n+1}=\widetilde{\mathsf F}_{dt/2}
\mathsf S_{dt}
\widetilde{\mathsf F}_{dt/2}W_n.
$$

This places second-order friction updates around a second-order Hamiltonian/diffusion step. The paper supplies the component formulas but does not specify this complete composition. Simply placing a full friction step after a full Hamiltonian/diffusion step would generally make their mutual splitting first order.

**Stability qualification.** The paper reports good stability for its examples, including its two-particle timestep. It does not establish an unconditional stability or positivity guarantee for this explicit polynomial update. Require timestep convergence, bounded numerical errors, and adequately resolved boundaries for each new problem.

## 7. Minimal implementation derived for this note

The following code implements the conventions above for a static real potential. It includes unitary, diffusion, classical-force, and Caldeira-Leggett modes. It uses complex128 arrays and out-of-place transforms for clarity; it is not the memory-optimized implementation used in the paper's four-dimensional demonstration.

`D` and `gamma` are accepted independently. For the Caldeira-Leggett model, the caller must impose `D = 2 * mass * gamma * kBT` in consistent units. Setting `Vprime` selects classical transport; it is the derivative of the potential, **not** the signed force.

```python
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import numpy as np

Potential = Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True)
class Grid1D:
    # W has shape (n_p, n_x); coordinate arrays broadcast to that shape.
    x: np.ndarray
    p: np.ndarray
    lam: np.ndarray
    theta: np.ndarray
    dx: float
    dp: float

    @classmethod
    def build(cls, n_x: int, n_p: int, L_x: float, L_p: float) -> Grid1D:
        for n in (n_x, n_p):
            if not isinstance(n, (int, np.integer)) or n < 4 or n % 2:
                raise ValueError("Grid sizes must be even integers >= 4.")
        if not np.all(np.isfinite([L_x, L_p])) or min(L_x, L_p) <= 0:
            raise ValueError("Half-widths must be finite and positive.")
        dx, dp = 2.0 * L_x / n_x, 2.0 * L_p / n_p
        x = np.fft.ifftshift((np.arange(n_x) - n_x // 2) * dx)[None, :]
        p = np.fft.ifftshift((np.arange(n_p) - n_p // 2) * dp)[:, None]
        lam = (2.0 * np.pi * np.fft.fftfreq(n_x, d=dx))[None, :]
        theta = (2.0 * np.pi * np.fft.fftfreq(n_p, d=dp))[:, None]
        return cls(x, p, lam, theta, dx, dp)

    @property
    def shape(self) -> tuple[int, int]:
        return self.p.size, self.x.size


def checked_state(W: np.ndarray, g: Grid1D) -> np.ndarray:
    W = np.asarray(W, dtype=np.complex128)
    if W.shape != g.shape or not np.all(np.isfinite(W)):
        raise ValueError("W must be finite and have shape (n_p, n_x).")
    return W


def from_wavefunction(psi: Potential, g: Grid1D, hbar: float = 1.0) -> np.ndarray:
    """Sample Eq. (5), then apply Eq. (8). psi must already be normalized."""
    if not np.isfinite(hbar) or hbar <= 0:
        raise ValueError("hbar must be finite and positive.")
    B = psi(g.x - 0.5 * hbar * g.theta) * np.conj(
        psi(g.x + 0.5 * hbar * g.theta)
    )
    B = checked_state(B, g)
    # This conversion needs physical scaling; internal split steps do not.
    return np.fft.ifft(B, axis=0) / g.dp


class SplitStep1D:
    """Static real potential, fixed step; optional diffusion and CL friction.

    Set Vprime to a callable to use classical force transport instead of
    the quantum potential difference. Rates D and gamma are independent;
    enforce D = 2*mass*gamma*kBT outside this class for the CL model.
    """

    def __init__(
        self, g: Grid1D, V: Potential, *, mass: float, dt: float,
        hbar: float = 1.0, D: float = 0.0, gamma: float = 0.0,
        Vprime: Potential | None = None,
    ) -> None:
        values = [mass, dt, hbar, D, gamma]
        if (not np.all(np.isfinite(values)) or min(mass, dt, hbar) <= 0
                or min(D, gamma) < 0):
            raise ValueError("Use positive mass/dt/hbar and nonnegative D/gamma.")
        self.g, self.dt, self.gamma = g, float(dt), float(gamma)
        self.K_half = np.exp(-0.5j * dt * g.p * g.lam / mass)
        if Vprime is None:
            # V must be defined at these shifted positions, not just at g.x.
            delta_V = (np.asarray(V(g.x - 0.5 * hbar * g.theta))
                       - np.asarray(V(g.x + 0.5 * hbar * g.theta)))
            omega = delta_V / hbar
        else:
            # Classical limit: (V_minus - V_plus)/hbar -> -theta*V'(x).
            omega = -g.theta * np.asarray(Vprime(g.x))
        omega = np.broadcast_to(omega, g.shape)
        if np.iscomplexobj(omega) or not np.all(np.isfinite(omega)):
            raise ValueError("Potential/derivative must give finite real values.")
        self.Q_full = np.exp(-1j * dt * omega - dt * D * g.theta**2)

    def _kinetic(self, W: np.ndarray) -> np.ndarray:
        return np.fft.ifft(self.K_half * np.fft.fft(W, axis=1), axis=1)

    def friction_rhs(self, W: np.ndarray) -> np.ndarray:
        # C W = 2*gamma*d_p(p*W), with multiplication BEFORE differentiation.
        return 2.0 * self.gamma * np.fft.ifft(
            1j * self.g.theta * np.fft.fft(self.g.p * W, axis=0), axis=0
        )

    def friction_rk2(self, W: np.ndarray, h: float) -> np.ndarray:
        # Eqs. (62)-(65): explicit second-order polynomial, not an exact flow.
        k1 = self.friction_rhs(W)
        return W + h * self.friction_rhs(W + 0.5 * h * k1)

    def step(self, W: np.ndarray) -> np.ndarray:
        W = checked_state(W, self.g)
        if self.gamma:
            W = self.friction_rk2(W, 0.5 * self.dt)
        W = self._kinetic(W)
        W = np.fft.ifft(self.Q_full * np.fft.fft(W, axis=0), axis=0)
        W = self._kinetic(W)
        if self.gamma:
            W = self.friction_rk2(W, 0.5 * self.dt)
        # Do not clip negativity, discard imaginary parts, or renormalize here.
        return W
```

A basic initialization and unitary smoke test is:

```python
g = Grid1D.build(n_x=128, n_p=128, L_x=10.0, L_p=8.0)
psi = lambda x: np.pi**(-0.25) * np.exp(-0.5 * (x - 0.6)**2 + 1.1j * x)
W = from_wavefunction(psi, g)  # hbar = 1
solver = SplitStep1D(g, lambda x: 0.5 * x**2, mass=1.0, dt=0.01)
for _ in range(100):
    W = solver.step(W)
assert abs(W.sum() * g.dx * g.dp - 1.0) < 1e-10
assert np.max(np.abs(W.imag)) < 1e-10
```

Mixed states can be initialized by summing pure-state Wigner functions with nonnegative weights summing to one, or by sampling the supplied density matrix in the definition of $B$. This follows from the linearity of Eqs. (5)-(8).

## 8. Diagnostics and acceptance tests

### 8.1 Quantities to record

Record the following with a declared quadrature measure $dx\,dp$:

$$
\mathcal N=\sum_{k,j}W_{kj}\,dx\,dp,
\qquad
\langle H\rangle=\sum_{k,j}
\left[\frac{p_k^2}{2m}+V(x_j)\right]W_{kj}\,dx\,dp.
$$

The paper's signed negative volume is

$$
N_W=\int_{W<0}W\,dx\,dp\leq0.
$$

If the repository also uses a nonnegative convention, give it a different name, for example $\nu=-N_W=(\int|W|\,dx\,dp-1)/2$ for a normalized real state. Do not change conventions between plots. [P, Eq. (48).]

Single-particle purity is

$$
\mathcal P=2\pi\hbar\int W(x,p)^2\,dx\,dp.
$$

This is Eq. (70), with $\hbar$ retained. Apply the formula to the physical real Wigner function; a substantial imaginary part is an error to diagnose, not another physical component. For $d$ spatial coordinates, the corresponding full-state prefactor is $(2\pi\hbar)^d$.

**Implementation choices.** Also record imaginary residuals, the fraction of $\int|W|$ near each boundary, high-frequency spectral weight, and marginal densities. Do not automatically renormalize every step: that can hide transport or boundary errors. Wigner negativity is not the same as failure of density-matrix positivity. For small validation systems, positivity can be checked by reconstructing and testing a density matrix, with the interpolation between coordinate systems explicitly controlled.

### 8.2 Analytical unit tests derived from the equations

| Test | Expected result | Main error detected |
|---|---|---|
| Fourier round trip | `ifft(fft(W))` returns `W`; physical $B$ also returns $W$ with the stated scaling | Transform direction, axis, or normalization |
| Gaussian wavefunction | Correct $x_0,p_0$, trace one, and purity one, including $\hbar\neq1$ | Initial-state conjugation, sign, or missing factors |
| Free particle | $W(x,p,t)=W_0(x-pt/m,p)$ | Kinetic phase sign and mass scaling |
| Constant force, $V=-Fx$ | $p$ advances by $Ft$ and $x$ by $p_0t/m+Ft^2/(2m)$ | Potential-difference sign and $\hbar$ convention |
| Quadratic potential | Quantum and classical-force multipliers agree | Incorrect classical limit or shifted potential |
| Diffusion only | $d\langle p^2\rangle/dt=2D$, mean momentum unchanged | Diffusion coefficient or sign |
| Friction only | $W(x,p,t)=e^{2\gamma t}W_0(x,e^{2\gamma t}p)$ | Missing product-rule term or factor of two |
| Free friction plus diffusion | Mean and variance follow the formulas below | Bath composition and overall timestep order |

For the last test, with no force and $\gamma>0$,

$$
\langle p\rangle_t=\langle p\rangle_0e^{-2\gamma t},
$$

$$
\operatorname{Var}(p)_t=
\operatorname{Var}(p)_0e^{-4\gamma t}
+\frac{D}{2\gamma}(1-e^{-4\gamma t}).
$$

Under the Caldeira-Leggett coefficient relation, the limiting momentum variance is $mk_BT$. This is a moment test, not proof that the approximate quantum master equation has the exact quantum Gibbs state as its stationary state.

For closed-system propagation, trace and the discrete squared norm should be conserved to numerical precision by the FFT/phase factors. Energy need not be exactly conserved at finite timestep. Validate harmonic motion against its exact phase-space rotation and verify approximately fourfold error reduction when a second-order method's timestep is halved, before spatial error dominates. Then converge grid resolution and domain sizes independently.

### 8.3 Local checks performed on the embedded code

**New implementation checks, not paper benchmark reproductions.** The embedded code was tested on a $128\times128$ grid with $L_x=10$, $L_p=8$, using a displaced Gaussian

$$
W_0=\pi^{-1}\exp[-(x-0.6)^2-(p-1.1)^2].
$$

Wavefunction initialization was also tested with $\hbar=0.7$. Initialization, free drift, constant force, quadratic quantum/classical agreement, and isolated diffusion tests agreed with their analytical results to approximately $10^{-14}$ or better.

For $V(x)=\tfrac12(1.3)^2x^2$, $m=\hbar=1$, and final time $t=1$, the absolute phase-space $L^2$ errors were:

| Timestep | Error against exact harmonic rotation |
|---|---:|
| 0.1 | $1.67692\times10^{-3}$ |
| 0.05 | $4.18475\times10^{-4}$ |
| 0.025 | $1.04572\times10^{-4}$ |

The error ratios were approximately 4.007 and 4.002. A friction-only test also showed approximately fourfold convergence. A free-particle friction-plus-diffusion test with $D=0.15$, $\gamma=0.2$, $dt=0.002$, and $t=0.4$ gave mean-momentum and variance errors below $3\times10^{-8}$. These checks do not establish arbitrary-grid stability, density-matrix positivity, or fidelity to the article's under-specified examples.

## 9. Single-particle demonstration: reported inputs and findings

### 9.1 Reported setup

**Paper, Sec. IV, PDF pp. 6-7.** A heavy particle in a Morse potential represents diatomic vibrational dynamics. The initial state is described as a displaced first-excited Morse state. The authors compare unitary evolution, position decoherence, Caldeira-Leggett dynamics, regularized Koopman-von Neumann evolution, and classical Fokker-Planck evolution.

| Quantity | Value reported in the article |
|---|---|
| Mass | $m=58\,752$ a.u. |
| Potential depth parameter | $V_0=0.6$ eV, reported as $0.0220$ a.u. |
| Range parameter | $a=2.5$ a.u. |
| Equilibrium-coordinate parameter | $r_e=-4.7$ a.u. |
| Initial-state displacement parameter | $x_0=4.3$ a.u., as printed |
| Displayed final time | $t=40\,400$ a.u. |
| Quantum decoherence coefficient | $D=2.70\times10^{-3}$ a.u. |
| Caldeira-Leggett bath | $T=300$ K; $\gamma^{-1}=41\,341$ a.u., stated as 1 ps |
| Quantum simulation grid | $512\times1024$ |
| Regularized Koopman-von Neumann grid | $768\times6144$ |
| Printed regularization parameter | $\delta D=1.5\times10^{-6}$ a.u. |
| Fokker-Planck grid | $512\times1024$ |
| Fokker-Planck coefficient in Fig. 2 caption | $D=2.61\times10^{-3}$ a.u. |

The dimension sizes are reported as in the captions; the table does not infer an undocumented axis assignment or simulation window from cropped plots. [P, Figs. 1-2; Eqs. (66)-(67).]

### 9.2 What the figures demonstrate

**Paper.** Fig. 1 on PDF p. 6 shows persistent negative structure under unitary evolution, but suppression of negativity in the displayed open-system results. Diffusion spreads the state, while adding damping changes its energy distribution. Fig. 2 compares the classical calculations; the decohering quantum distribution resembles the diffusive Fokker-Planck result rather than the fine filamentary Koopman-von Neumann amplitude. Fig. 3 on PDF p. 7 illustrates spurious growth in the magnitude of the signed negative volume without classical regularization. [P, Sec. IV; Figs. 1-3.]

The interpretation is that decoherence introduces backaction diffusion. The appropriate classical comparison therefore retains that diffusion, unless it is negligible on the action scale of interest. The figures illustrate this behavior for the reported examples; positivity of a Wigner function alone is not an implementation test of exact equality with classical evolution. [P, Sec. IV, pp. 6-7.]

### 9.3 Reproduction blockers: do not silently resolve them

**Potential sign.** Printed Eq. (66) is

$$
V(x)=V_0\{e^{-2a(x-r_e)}-2e^{+a(x-r_e)}\}.
$$

As a direct check, its derivative at $r_e$ is $-4aV_0$, so $r_e$ is not a minimum of the printed function. A plausible alternative consistent with the stated Morse-well interpretation is

$$
V_{\rm candidate}(x)=V_0\{e^{-2a(x-r_e)}-2e^{-a(x-r_e)}\}.
$$

This latter expression is an **inferred repair**, not the formula printed in the supplied article or a verified author correction. It must be documented as a modeling choice unless resolved from the authors' implementation.

**Initial state.** Eq. (67) contains $z$ without defining it in the displayed surrounding text; its exponent also retains a state-index symbol $n$. The stated first-excited state suggests $n=1$, but the absent $z$ definition still prevents direct evaluation. The positive printed $x_0$ and the negative-coordinate location in Fig. 1 also require an explicit displacement convention. Do not substitute a guessed analytic wavefunction and call it an exact reproduction.

**Diffusion coefficients.** Fig. 1 reports $2.70\times10^{-3}$, while Fig. 2 reports $2.61\times10^{-3}$, even though the prose describes using the same coefficient. Preserve both reported values until the discrepancy is resolved.

**Other missing details.** The article does not provide complete one-particle timestep, domain extents, and initialization/normalization procedures for every classical comparison. In particular, its statement that the same initial state is propagated does not fully specify a numerical conversion between a signed Wigner function, a normalized Koopman-von Neumann amplitude, and a positive classical density. The regularization-rate versus per-step-filter distinction also remains relevant.

**Implementation choice:** Validate the method using the analytical cases in Section 8 before attempting this benchmark. A numerical eigenstate of an explicitly chosen Morse potential is a useful substitute experiment, but not automatically the state of Eq. (67).

## 10. Two-particle demonstration and extension

### 10.1 Reported system

**Paper, Sec. V, PDF p. 8.** The state is $W_2(x,p_x,y,p_y)$ and the potential is

$$
V(x,y)=\frac12(x^2+y^2)+\frac1{10}(x^4+y^4+xy).
$$

Only the $x$ particle is directly coupled to the Caldeira-Leggett bath, with

$$
D_x=0.04,\qquad \gamma_x=1/12.5
$$

in atomic units. The $y$ particle is affected through its coupling to the first particle. The initial wavefunction is described as the antisymmetric, "fermionic-like" entangled state

$$
\psi_F(x,y)=\frac{1}{\sqrt2}
[\psi_1(x)\psi_2(y)-\psi_1(y)\psi_2(x)],
$$

with Gaussians centered at $+1$ and $-1$. The reported grid is $128\times192\times128\times192$, and the calculation uses single-precision arithmetic. The authors report approximately 4.7 GB for the stored state, require two state copies for the friction update, and report stable evolution at $dt=0.01$ a.u. The displayed reduced-state time is $t=5$ a.u. [P, Eqs. (68)-(72); Figs. 4-5.]

Fig. 4 on PDF p. 7 shows initially equal reduced states that later differ. The directly bath-coupled particle's reduced Wigner function is smoother and less negative. Fig. 5 on PDF p. 8 shows a larger purity decrease for that particle and a smaller decrease for its partner. These are reduced-state results, not a claim that reduced purity alone isolates decoherence from changes in interparticle entanglement.

### 10.2 Implementing the higher-dimensional substeps

**Implementation derivation.** For

$$
H=\frac{p_x^2}{2m_x}+\frac{p_y^2}{2m_y}+V(x,y),
$$

use the kinetic multiplier

$$
K_h=\exp\left[-ih\left(\frac{p_x\lambda_x}{m_x}
+\frac{p_y\lambda_y}{m_y}\right)\right]
$$

and the potential difference

$$
\Delta V=V\!\left(x-\frac{\hbar\theta_x}{2},y-\frac{\hbar\theta_y}{2}\right)
-V\!\left(x+\frac{\hbar\theta_x}{2},y+\frac{\hbar\theta_y}{2}\right).
$$

Both coordinates are shifted together, including the coupling term. For the paper's bath placement,

$$
Q_h=\exp[-ih\Delta V/\hbar-hD_x\theta_x^2],
\qquad
\mathsf C_xW_2=2\gamma_x\partial_{p_x}(p_xW_2).
$$

There is no direct $D_y$ or $\gamma_y$ term in that example. Potential evolution requires momentum-axis transforms; kinetic evolution requires position-axis transforms. The friction algorithm is unchanged apart from its selected momentum axis.

**Implementation choice:** One explicit storage contract is `(N_px, N_x, N_py, N_y)`. Position FFTs then act on axes `(1, 3)` and momentum FFTs on axes `(0, 2)`. This is a proposed contract, not an axis order specified by the reported four grid sizes.

The reduced states are

$$
W_x(x,p_x)=\int W_2\,dy\,dp_y,
\qquad
W_y(y,p_y)=\int W_2\,dx\,dp_x.
$$

For the storage contract above:

```python
W_x = W_2.sum(axis=(2, 3)) * dy * dp_y
W_y = W_2.sum(axis=(0, 1)) * dx * dp_x
```

Compute each reduced purity using $2\pi\hbar\int W_{x/y}^2$, and full-state purity using $(2\pi\hbar)^2\int W_2^2$. A pure entangled full state can have mixed reduced states. [P, Eqs. (69)-(70).]

### 10.3 Initialization and memory caveats

**Under-specified inputs.** Gaussian widths, full grid extents, and an explicit mass assignment for the two-particle example are not given in Sec. V. Unit masses may be a reasonable chosen convention, but must not be presented as an explicitly supplied parameter.

**Normalization derivation.** The factor $1/\sqrt2$ in Eq. (72) normalizes a determinant of orthonormal orbitals. For individually normalized Gaussians with overlap $S=\langle\psi_1|\psi_2\rangle$, the normalization factor is instead

$$
\frac{1}{\sqrt{2(1-|S|^2)}}.
$$

Separated Gaussians are not automatically orthogonal. Record their widths and either their overlap-normalized state or an explicit orthonormalization. Do not silently assume zero overlap.

**Memory arithmetic, derived here.** The reported dimensions contain 603,979,776 samples. A real float32 array alone takes about 2.42 decimal GB; a complex64 array about 4.83 decimal GB; a complex128 array about 9.66 decimal GB. Account separately for state copies, FFT workspace, cached multipliers, and temporary products. The article's approximately 4.7 GB report is not the total memory requirement of a straightforward multi-buffer implementation.

The simple one-dimensional code above is not suitable for scaling directly to that grid unchanged. Check backend dtype preservation and workspace allocation, use broadcast coordinate arrays rather than full coordinate meshes, and control which full-size phase arrays are cached. These are engineering recommendations derived from the array sizes, not measured performance claims.

## 11. Suggested repository implementation plan

**Implementation choices.** Keep the first version narrow: static separable kinetic energy, a callable real potential, one phase-space degree of freedom, and the paper's two bath terms. Establish correctness before introducing additional master equations or the full four-dimensional grid.

A useful module split is:

```text
src/<package>/
    grids.py          # Axis contracts, FFT-order coordinates, quadrature measures.
    states.py         # Wavefunction/density-matrix initialization and transforms.
    propagators.py    # Kinetic, potential/diffusion, and friction substeps.
    diagnostics.py    # Trace, marginals, energy, purity, negativity, boundaries.
tests/
    test_transforms.py
    test_analytic_dynamics.py
    test_timestep_convergence.py
examples/
    harmonic_validation.py
    morse_comparison.py
    coupled_particles.py
docs/papers/
    cabrera_2015_wigner_implementation.md
```

Store grid sizes and half-widths, axis order, units, $\hbar$, masses, potential parameters, $dt$, final time, $D$, $\gamma$, precision, initial-state definition, and any filter strength with every result. For a thermal bath also store $k_BT$ and verify the relation to $D$; do not substitute a temperature expressed in kelvin directly for an energy.

The practical progression is: validate transforms and Gaussian initialization; implement and converge closed-system Strang splitting; add diffusion; add and converge friction; add the classical-force option with explicit amplitude/density semantics; then test a small coupled two-particle grid before attempting the paper's dimensions. Treat literal figure reproduction as a separate milestone requiring resolution of the benchmark ambiguities.

## 12. Boundaries of what this paper supplies

The article supplies the Hilbert-phase-space formulation, explicit FFT-based propagation ingredients, a practical second-order friction approximation, and illustrative quantum/classical comparisons. It does **not** supply a complete unambiguous parameter set for every demonstration, a generally positivity-preserving Caldeira-Leggett discretization, an unconditional friction stability proof, or a complete ready-made efficient implementation for arbitrary Lindblad operators.

The implementation above does not add non-Markovian memory, derive a bath from molecular data, fit a potential, support a coordinate-dependent mass, or implement time-dependent Hamiltonians. Those would be separately specified extensions, not capabilities established by this summary.

## 13. Source map and bibliography

| Topic | Location in [P] |
|---|---|
| Density matrix, $B$, and $W$ | Sec. II A, Eqs. (3)-(14), PDF pp. 1-2 |
| Hilbert phase-space operators and representation diagram | Sec. II A, Eqs. (15)-(29), PDF pp. 2-3 |
| Lindblad mapping, decoherence, Caldeira-Leggett model | Sec. II B, Eqs. (30)-(41), PDF pp. 3-4 |
| Classical amplitudes, densities, cumulative functions, diffusion | Sec. II C, Eqs. (42)-(49), PDF p. 4 |
| Split operators, discrete grids, open-system steps | Sec. III, Eqs. (50)-(65), PDF pp. 4-6 |
| Morse example and quantum/classical comparisons | Sec. IV, Eqs. (66)-(67), Figs. 1-3, PDF pp. 6-7 |
| Coupled particles, reduced states, and purity | Sec. V, Eqs. (68)-(72), Figs. 4-5, PDF pp. 7-8 |
| Pointer to supplementary Python programs | Ref. [70], PDF p. 9; not inspected for this note |

```bibtex
@article{Cabrera2015Wigner,
  author  = {Cabrera, Renan and Bondar, Denys I. and Jacobs, Kurt
             and Rabitz, Herschel A.},
  title   = {Efficient method to generate time evolution of the Wigner
             function for open quantum systems},
  journal = {Physical Review A},
  volume  = {92},
  pages   = {042122},
  year    = {2015},
  doi     = {10.1103/PhysRevA.92.042122}
}
```
