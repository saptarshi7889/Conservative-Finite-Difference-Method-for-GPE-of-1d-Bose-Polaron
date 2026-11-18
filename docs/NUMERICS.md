# Equations and numerical method

[Quick start](../README.md) · [Validation and limitations](../VALIDATION.md)

## Equation and units

The coefficient structure implements the dimensionless equations

```math
i\partial_t\phi = \left[-\frac{1+m/M}{2}\partial_x^2
 + |\phi|^2 + \widetilde G\delta(x) + i V(t)\partial_x\right]\phi,
```

```math
V(t)=V_0+\frac{i(m/M)}{\sqrt\gamma}\int\phi^*\partial_x\phi\,dx.
```

Here $x$ is measured in healing lengths, $t$ in inverse chemical-potential units,
$\phi$ in square roots of the background density, and $V$ in sound-speed units.
These conventions correspond to the coefficients of the supplied code.
The kinetic factor $1+m/M$ selects Eq. (5), rather than the alternative Eq. (7)
discussed in the paper. The discrete velocity routine uses symmetrization over both adjacent time slices,
as detailed below.


## Numerical algorithm

This section summarizes *GPE_simulation.pdf*, Chapter 1, printed pages 2–4,
Eqs. (1.3)–(1.7), before Section 1.1, and compares those notes with
[`src/solver.cpp`](../src/solver.cpp) and [`src/numerics.hpp`](../src/numerics.hpp).
The PDF is not bundled with this repository. Its Eqs. (1.6)/(1.7) describe the
updates referenced as Eqs. (4.6)/(4.7) in the inherited source comments.
Tildes on dimensionless quantities are omitted here; write $\mu=m/M$,
$h=\Delta x$, $\tau=\Delta t$, and $N=L/h$.

### Spatial and temporal discretization

The field is stored on $N$ periodic points, with indices interpreted modulo $N$.
Define the unscaled second difference and the velocity-selected first difference:

```math
(L_1\Phi)_j=\Phi_{j+1}-2\Phi_j+\Phi_{j-1},\qquad
(L_2(V)\Phi)_j=
\begin{cases}
\Phi_{j+1}-\Phi_j,&\mathrm{Re}(V)\geq0,\\
\Phi_j-\Phi_{j-1},&\mathrm{Re}(V)<0.
\end{cases}
```

For the ordered field vector $\Phi=(\Phi_0,\Phi_1,\ldots,\Phi_{N-1})^{\mathsf T}$,
the corresponding matrices are

```math
L_1=
\begin{pmatrix}
-2 & 1 & 0 & \cdots & 1 \\
1 & -2 & 1 & \ddots & 0 \\
0 & 1 & \ddots & \ddots & \vdots \\
\vdots & \ddots & \ddots & -2 & 1 \\
1 & 0 & \cdots & 1 & -2
\end{pmatrix}_{N\times N}.
```

For $\mathrm{Re}(V)\geq0$, the forward-difference matrix is

```math
L_2^{+}=
\begin{pmatrix}
-1 & 1 & 0 & \cdots & 0 \\
0 & -1 & 1 & \ddots & \vdots \\
\vdots & \ddots & \ddots & \ddots & 0 \\
0 & \cdots & 0 & -1 & 1 \\
1 & 0 & \cdots & 0 & -1
\end{pmatrix}_{N\times N}.
```

For $\mathrm{Re}(V)<0$, the backward-difference matrix is

```math
L_2^{-}=
\begin{pmatrix}
1 & 0 & \cdots & 0 & -1 \\
-1 & 1 & 0 & \ddots & 0 \\
0 & -1 & \ddots & \ddots & \vdots \\
\vdots & \ddots & \ddots & 1 & 0 \\
0 & \cdots & 0 & -1 & 1
\end{pmatrix}_{N\times N},\qquad
L_2(V)=
\begin{cases}
L_2^{+},&\mathrm{Re}(V)\geq0,\\
L_2^{-},&\mathrm{Re}(V)<0.
\end{cases}
```

The corner entries implement periodic boundaries. These are the unscaled
matrices shown in the notes: the derivative operators are $L_1/h^2$ and
$L_2(V)/h$. Their dimension is $N\times N$, where $N=L/h$ is the number of
grid points, rather than the physical system length $L$.

The forward stencil for positive velocity follows the sign of the drift in
$\partial_t\phi=V\partial_x\phi+\cdots$, whose transport speed is $-V$.
The delta potential is represented by $D_j=(\widetilde G/h)\delta_{j,j_*}$,
where $j_*=\lfloor N/2\rfloor-1$ is the inherited implementation location
(see the grid-offset note above).

The kinetic and drift terms are averaged between the old and new time slices.
The nonlinear term uses the average of the two densities multiplied by the
average field. With $W=V^n$ held fixed during the step, the implicit update is

```math
i\frac{\Phi^{n+1}-\Phi^n}{\tau}
=-\frac{1+\mu}{4h^2}L_1(\Phi^{n+1}+\Phi^n)
+\frac{iW}{2h}L_2(W)(\Phi^{n+1}+\Phi^n)
+\frac12\left[D+\frac{|\Phi^{n+1}|^2+|\Phi^n|^2}{2}\right]
\odot(\Phi^{n+1}+\Phi^n).
```

Here $\odot$ and the density operations act componentwise. This is a nonlinear
implicit scheme; the one-sided spatial drift stencil is first order in $h$.

### Picard iteration and cyclic solve

For each time step, set $\Psi_0=\Phi^n$ and form

```math
 a=\frac{i\tau(1+\mu)}{4h^2},\qquad b=\frac{\tau W}{2h},\qquad
 A_1=I-aL_1-bL_2(W),\qquad A_2=I+aL_1+bL_2(W).
```

Each Picard iteration solves the linear system

```math
A_1\Psi_{k+1}=A_2\Phi^n
-\frac{i\tau}{2}\left[D+\frac{|\Psi_k|^2+|\Phi^n|^2}{2}\right]
\odot(\Psi_k+\Phi^n).
```

`Temp` stores $A_2\Phi^n$; `Phi_k`, `Phi_tilde`, and `Phi_hat` store
$\Psi_k$, the right-hand side, and $\Psi_{k+1}$, respectively. Only vectors
and scalar stencil coefficients are stored, rather than dense matrices.
Periodicity makes $A_1$ cyclic tridiagonal. `Thomas_tridiagonal_solver_PBC`
uses two ordinary Thomas solves and a rank-one correction, giving $O(N)$
work and storage per Picard iteration.

Accept $\Phi^{n+1}=\Psi_{k+1}$ when
$\|\Psi_{k+1}-\Psi_k\|_2\leq\varepsilon$. This is an absolute, unweighted
Euclidean norm, not a relative error or an equation-residual test. The notes
specify $\varepsilon=10^{-9}$; the packaged CLI defaults to `--tolerance 1e-8`
and at most `--iterations 12` iterations. It stops with an error if convergence
is not reached within that limit.

### Coupled velocity update

After accepting the field, define $S_j=\Phi_j^n+\Phi_j^{n+1}$. The routine
`Velocity_calculation` evaluates the notes' symmetrized expression as

```math
V^{n+1}=V_0+\frac{i\mu}{8\sqrt\gamma}
\sum_{j=0}^{N-1}\left[
S_j^*(S_{j+1}-S_j)-S_j(S_{j+1}^*-S_j^*)\right].
```

The factors of $h$ from the spatial integral and first difference cancel.
The initial fresh-run field is $\Phi_j^0=1$, with $V^0=V_0$. The updated
velocity selects the stencil and coefficients for the **next** time step;
it is not iterated together with the field inside the current Picard loop.

### Conservation and scope

The notes base the nondrift discretization on the conservative finite-difference
scheme of V. A. Trofimov and N. V. Peskov, “Comparison of finite-difference schemes
for the Gross-Pitaevskii equation,” *Mathematical Modelling and Analysis*
**14**(1), 109–126 (2009), their reference [1]. They describe norm conservation
without drift and negligible norm decay with drift. These statements do not
establish exact conservation for the implemented one-sided drift scheme or a
finite Picard tolerance. Monitor $h\sum_j|\Phi_j|^2$ and check spatial, temporal,
and iteration convergence for the parameters of interest. The measured drift
and limits of the existing checks are documented in [VALIDATION.md](../VALIDATION.md).

