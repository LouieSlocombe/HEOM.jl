# HEOM.jl

Hierarchical equations of motion (HEOM) is a theoretical framework used in open quantum systems to describe the dynamics of quantum systems interacting with their environment. In open quantum systems, the system of interest interacts with its surrounding environment, causing decoherence and dissipation effects. HEOM provides a systematic approach to modelling the time evolution of the reduced density matrix of the system hierarchically. The HEOM approach decomposes the reduced density matrix into a hierarchy of auxiliary density matrices, each representing a different level of correlation between the system and the environment. The equations of motion for these auxiliary density matrices form a set of coupled differential equations that describe the system's time evolution.

The hierarchy is typically truncated at a certain level to make the calculations computationally feasible. Therefore, the level of truncation determines the accuracy of the results obtained using the HEOM approach.

By solving the hierarchical equations of motion, one can obtain the dynamics of the reduced density matrix, which provides information about the system's evolution, including decoherence, relaxation, and other environmental dissipative effects. This information is crucial for understanding and predicting the behaviour of open quantum systems and their interactions with the surrounding environment.

In this code repo, we work in phase-space!

# Linear-linear Coupling in High-Temperature Markovian Limit
$$
\frac{\partial}{\partial t} W^{(0)}=  -\hat{\mathcal{L}}_{\mathrm{QM}} W^{(0)}
+\zeta \frac{\partial}{\partial p}\left(p+\frac{m}{\beta} \frac{\partial}{\partial p}\right) W^{(0)}
$$

# Linear-linear Coupling in High-Temperature non-Markovian Limit
$$
\frac{\partial}{\partial t} W^{(n)}=
-\left(\hat{\mathcal{L}}_{\mathrm{QM}}+n \gamma\right) W^{(n)}\\
+\frac{\partial}{\partial p} W^{(n+1)}\\
+n \gamma \zeta\left(p+\frac{m}{\beta} \frac{\partial}{\partial p}\right) W^{(n-1)}
$$
