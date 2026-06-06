# MARS Pump Flow Solver — CFD Engineer's Guide

A working reference for the incompressible Navier–Stokes pump-flow capability in
MARS: what it solves, how it discretizes, why it makes the choices it does, how to
read its output, and where it is fragile. Written for someone with a CFD/numerics
background who needs to defend the method and compare it against an external code.

> **New to this?** For a from-scratch explanation (what a pump simulation even is,
> mesh → method → boundary conditions → running → reading output, assuming no prior
> CFD), read [pump_cfd_tutorial.md](pump_cfd_tutorial.md) first. For making flow
> videos, see [scripts/render_pump_flow.py](../scripts/render_pump_flow.py).

Conventions used below:
- **[verified]** — implemented and checked (unit test, invariant, or measurement).
- **[weak]** — known limitation / open issue.
- **[inferred]** — plausible but not directly confirmed; treat as hypothesis.

---

## 1. Problem class and governing equations

Incompressible, isothermal, Newtonian, single-phase flow. Constant density $\rho$
and kinematic viscosity $\nu$.

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu\,\nabla^2 \mathbf{u}$$

$$\nabla\cdot\mathbf{u} = 0$$

The defining difficulty is the **constraint** $\nabla\cdot\mathbf{u} = 0$. There is no evolution
equation for pressure; pressure is the Lagrange multiplier that enforces
incompressibility. Every incompressible scheme is, at heart, a strategy for
coupling `p` and `u` so the constraint holds. All design choices below follow from
that.

**Flow regime.** `Re = U L / nu`. For the pump runs `Re = 500` is an *imposed test
value* — `nu` is back-computed as `U L / Re` after the bounding-box length `L` is
known. Re = 500 is transitional: inertial enough for separation and recirculation,
laminar enough that **no turbulence model is used**. The runs are DNS-like
(DNS = Direct Numerical Simulation: resolving all flow scales with no turbulence
model) but under-resolved; we do not claim resolved turbulence.

**Phase scope.** This is **flow-only**: the diaphragm and check valves are *frozen*
(rigid). The real pump is a fluid–structure interaction (FSI) problem (flexing
diaphragm + moving valves); that is deferred. Flow-only is the right first
validation step and is what the developers' legacy code can also produce, enabling
a direct comparison.

---

## 2. Spatial discretization — CVFEM (vertex-centered, median-dual)

MARS uses **CVFEM (Control-Volume Finite Element Method)** on **tetrahedra**. This is the
edge-based / median-dual finite-volume family (the Nalu-Wind lineage), *not*
Galerkin FEM and *not* cell-centered FV.

- **Unknowns at vertices** (P1-like): one `(u, v, w, p)` per node.
- **Median-dual control volumes**: each vertex owns a polyhedral CV built from
  edge-midpoints, face-centroids, and cell-centroids. The CV faces are
  **sub-control-surfaces (SCS)** — one quad facet per tet edge.
- **Per-edge SCS area vectors $\mathbf{A}_f$**: precomputed once
  (`precomputeTetAreaVectorsGpu`), oriented L→R with a consistent sign convention
  ($\mathbf{A}\cdot(\mathbf{x}_R - \mathbf{x}_L) > 0$). A tet has 6 edges → 6 SCS facets.

All discrete operators are **edge/face-based scatters**: each flux $(\cdot)\cdot\mathbf{A}_f$ is
accumulated to the two endpoint CVs with opposite sign. This gives **discrete local
conservation** by construction (telescoping).

**Discrete operators:**
- **Divergence $D$**: $(\nabla\cdot\mathbf{u})_i \approx \frac{1}{V_i}\sum_f (\mathbf{u}_f\cdot\mathbf{A}_f)$, with
  $\mathbf{u}_f = \tfrac{1}{2}(\mathbf{u}_L + \mathbf{u}_R)$.
- **Gradient $G$** — built two ways (this distinction matters, see §4):
  - *Green–Gauss / SCS*: $(\nabla p)_i \approx \frac{1}{V_i}\sum_f p_f\,\mathbf{A}_f$.
  - *$D^T$ (discrete transpose of divergence)*: same-sign scatter
    $\tfrac{1}{2}(p_L - p_R)\,\mathbf{A}_f$, integrating to $-V\,\nabla p$. **This is the gradient used
    in the projection** (`useLegacyGradient = false`).
- **Lumped mass $M$**: $M_i = V_i$ (CV volume), diagonal.

**Verified invariants** [verified]: $D(\text{const}) = 0$, $G(\text{linear}) = \text{const}$, and the
area-vector closure $\sum_f \mathbf{A}_f\cdot(\mathbf{x}_R - \mathbf{x}_L) = 3V$ per tet. Lumped mass is
rank-consistent: $\sum M$ over owned nodes is bit-identical on 1 and 4 ranks and
equals the domain volume.

---

## 3. Pressure–velocity coupling — Chorin projection (fractional step)

Segregated projection (not monolithic, not SIMPLE/PISO iteration). Per timestep:

**(a) Predictor** — momentum without the constraint:
$$\mathbf{u}^* = \mathbf{u}^n - \Delta t\,(\mathbf{u}\cdot\nabla)\mathbf{u} - \frac{\Delta t}{\rho}\nabla p^n \quad (+\text{ viscous, handled in b})$$
$\mathbf{u}^*$ is generally not divergence-free.

**(b) Implicit diffusion** — the viscous term is stiff ($\nu\,\Delta t / h^2$ blows up on
fine meshes), so it is solved implicitly:
$$\left(\frac{M}{\Delta t} + \nu K\right)\mathbf{u}^{**} = \text{rhs}(\mathbf{u}^*)$$
$K$ = assembled CVFEM stiffness (Laplacian). Solved by CG (Conjugate Gradient;
`cg_uvw` in the output).
The operator diagonal is rank-consistent [verified] (measured identical on 1 vs 4
ranks).

**(c) Pressure Poisson** — find $\phi$ making the corrected velocity divergence-free:
$$A\,\phi = \frac{\rho}{\Delta t}\,D\,\mathbf{u}^{**}$$

**(d) Corrector**:
$$\mathbf{u}^{n+1} = \mathbf{u}^{**} - \frac{\Delta t}{\rho}\,G\phi, \qquad p^{n+1} = p^n + \phi$$

**Time integration — BDF2 / EXT2** (BDF2 = 2nd-order Backward Differentiation
Formula for the time derivative; EXT2 = 2nd-order explicit extrapolation of the
advection term):
- Time derivative: $\dfrac{3\mathbf{u}^{n+1} - 4\mathbf{u}^n + \mathbf{u}^{n-1}}{2\,\Delta t}$ (2nd order).
- Advection extrapolated explicitly: $N_{\text{ext}} = 2N^n - N^{n-1}$ (Adams–Bashforth-like).
- **Consequence:** the explicit advection has a real **advective CFL limit**
  (CFL = Courant–Friedrichs–Lewy: the stability bound $u\,\Delta t / h \lesssim O(1)$).
  `--bdf1` switches to BDF1/Chorin, whose explicit term is more
  stable. A developing jet can cross the CFL limit at fixed `dt`; remedies are a
  smaller `dt`, `--bdf1`, or the skew scheme (no strict limit). [weak — CFL is real;
  see also §8 for a multi-rank issue that is *not* CFL.]

---

## 4. The pressure operator — matrix-free DDT (`D M^-1 D^T`)

The most important design choice.

**Why not the Galerkin Laplacian `K`?** The natural pressure operator is the
assembled stiffness `K`. On tets this is **inconsistent with the corrector** and
blows up. Mechanism:
- The corrector applies the SCS gradient $G = M^{-1}D^T$.
- For the projection to zero divergence you need $D(\text{corrector}) = 0$, i.e. the
  Poisson operator must be exactly $D M^{-1} D^T$ — the same $D$ and $G$ you actually
  use.
- $K$ is not $D M^{-1} D^T$. On tets the SCS gradient and the Galerkin gradient are
  nearly **orthogonal** (measured $\cos \approx -0.21$ on tets vs $-0.73$ on hex). Solving
  $K\phi = \text{rhs}$ then correcting with $M^{-1}D^T\phi$ leaves a large residual
  divergence — on tets it *amplifies* divergence each step → blow-up.

**The fix — DDT.** Use the literal algebraic composition $A = D M^{-1} D^T$. Then
$$D\,\mathbf{u}^{n+1} = D\,\mathbf{u}^{**} - \frac{\Delta t}{\rho}\,D M^{-1} D^T \phi = 0$$
holds **by construction** when $\phi$ solves $A\phi = \frac{\rho}{\Delta t}\,D\,\mathbf{u}^{**}$. The projection
is algebraically exact; $\text{div}_{\max} \sim \text{roundoff}$ in the ideal (single-rank) case.

**Matrix-free.** $A$ is never assembled. Applying $A$ inside CG is three kernels:
$D^T$ (gradient scatter) → $M^{-1}$ (divide by lumped mass) → $D$ (divergence
scatter), each with the halo exchanges needed for cross-rank consistency. This is
**halo-complete and verified correct** [verified]. Matrix-free is required for the
billion-element target — assembling $D M^{-1} D^T$ as CSR would exhaust VRAM.

**Solver.** Jacobi-preconditioned CG (`cg_p`). The Jacobi diagonal is computed
matrix-free: $\text{diag}_i = \sum_f 0.25\,|\mathbf{A}_f|^2\,(1/M_L + 1/M_R)$. [weak] The diagonal is a
weak preconditioner on stretched boundary-layer tets (diag spread of several orders
of magnitude → 300–600 CG iterations). AMG (Hypre BoomerAMG) is the obvious lever
but fights the matrix-free design and was reverted; it remains a future option.

**Null space.** For all-wall (pure Neumann) or periodic cases, `A` has a constant-
mode null space, removed by subtracting the mean of the RHS and of `phi`. For the
pump (Dirichlet pressure at the outlet) it is pinned — no null space.

---

## 5. Advection schemes — the developer-comparison crux

The convective term $(\mathbf{u}\cdot\nabla)\mathbf{u}$ discretization dominates accuracy and stability.
Face mass flux $\dot{m}_f = \rho\,\mathbf{u}_f\cdot\mathbf{A}_f$; schemes differ in the convected value at
the SCS:

1. **1st-order upwind** — $q_f = q_{\text{upwind}}$. Robust but $O(h)$ numerical diffusion
   $\sim |\mathbf{u}|\,h$; on a coarse mesh it dominates physical $\nu$ (effective Re much lower
   than the true Re). The developers' verdict — "crap unless you have shocks" — is
   correct.

2. **Skew-symmetric (Verstappen–Veldman)** — central, with the discrete identity
   $\sum_i q_i\,(dq/dt)_i = 0$ for any velocity field: **discretely kinetic-energy
   conserving**, zero artificial dissipation. The DNS-quality choice (NekRS-style).
   No strict advective CFL (injects no energy). [verified] Default; stable to steady
   state on multiple pump meshes single-rank.

3. **Barth–Jespersen (2nd-order limited upwind)** — linear reconstruction
   $q_f = q_{\text{up}} + \psi\,\nabla q_{\text{up}}\cdot(\mathbf{x}_f - \mathbf{x}_{\text{up}})$ with a slope limiter $\psi \in [0,1]$
   that drops toward 1st-order near extrema to prevent over/undershoots (monotone /
   TVD-like; TVD = Total Variation Diminishing). **This is the developers' scheme**, implemented for apples-to-apples
   comparison. Uses a Green–Gauss nodal gradient plus per-node neighbor min/max for
   the limiter.

**Why it matters for comparison.** Running the *same* scheme (Barth–Jespersen)
means result differences are physics/setup, not "different numerics." Skew vs BJ
also brackets dissipation: skew = zero numerical diffusion, BJ = a little near
gradients. If they disagree, the gap measures how dissipation-sensitive this flow
is.

[inferred / open] BJ accuracy on stretched boundary-layer tets depends on the
gradient reconstruction. MARS uses Green–Gauss, which is only ~1st-order and biased
on skewed tets; least-squares is the more accurate alternative. **Ask the developers
which they use** — it can be a real source of difference.

---

## 6. Boundary conditions

| Surface                        | Mathematics                          | Implementation |
|--------------------------------|--------------------------------------|----------------|
| Inlet (suction opening)        | $\mathbf{u} = U\,\hat{\mathbf{n}}_{\text{in}}$ (Dirichlet velocity), $\partial p/\partial n = 0$ | velocity nodes pinned; $U = 0.5$, $\hat{\mathbf{n}}_{\text{in}}$ = area-weighted **inward face normal** from the side-set triangles, **globally reduced over all ranks** so the direction is identical on every rank |
| Outlet (pressure opening)      | "do-nothing": $p = 0$ (Dirichlet pressure), $\mathbf{u}$ free | pressure pinned on the face; velocity evolves naturally |
| Walls + frozen diaphragm       | $\mathbf{u} = 0$ (no-slip)                    | every side-set **not** named inlet/outlet → no-slip bucket |

**BC assignment rule.** The driver takes two side-set names (`--inlet-ss=`,
`--outlet-ss=`); **every other side-set becomes a no-slip wall automatically**. So
the frozen-diaphragm opening, being neither inlet nor outlet, is walled.

**Well-posedness.** Velocity-Dirichlet at inlet + walls and pressure-Dirichlet at
the outlet is the standard inflow/outflow setup. A mass-conserving outlet variant
(velocity-Dirichlet outflow + single pressure pin) over-constrains and was seen to
blow up; the "do-nothing" outlet is the stable choice.

**Non-dimensionalization.** $L$ = bbox diagonal, $U$ = inlet speed, $\nu = U L / Re$.
Reported: $u_{\max}/U$, $\;u_{\text{rms}}/U = \dfrac{\left(\int |\mathbf{u}|^2\,dV\right)^{1/2}}{\sqrt{V}\,U}$, $\;\text{div}\,L/U$.
$Re = 500$ is an **input** (chosen test value); the rest are **outputs**.

**Two honest gaps for the developer conversation:**
1. `Re = 500` is a test value, **not** the pump's operating Re — ask what Re to match.
2. "Diaphragm-chamber opening" (our mesh side-set name) vs "diaphragm outlet" (their
   term): same word, **unconfirmed same surface** — ask them to confirm.

---

## 7. Reading a run line

```
# under-developed (early transient):
Step 200   ft=0.23  u_rms=3.26e-4  u_max=1.50  uMax/U=3.00  d(u_rms)=2.6e-6  div*L/U=381  cg_p=624  cg_uvw=364
# developed / converged (4-GPU, Re=100, skew, do-nothing outlet):
Step 3000  ft=3.48  u_rms=3.02e-1  u_max=6.68  uMax/U=13.4  d(u_rms)=9.8e-6  div*L/U=2.66 cg_p=554  cg_uvw=63
```
The contrast is the lesson: at `ft=0.23` the flow has barely entered (`div*L/U=381`,
a big *transient*); by `ft=3.5` it is developed and converged (`div*L/U=2.7`,
`d(u_rms)≈1e-5`). Always judge health on a developed run, not the first few steps.

- `ft` — flow-throughs = $t / (L/U)$. Convergence needs several (the transient must
  flush). `ft = 0.23` means **not developed yet**; a comparable small-mesh run
  reached `ft ~ 2`.
- `u_rms`, `uMax/U` — volume-RMS and peak velocity over inlet speed. `uMax/U = 3`
  means 3x acceleration through the constriction.
- `d(u_rms)` — steady-state residual. `2.6e-6` is flat — but possibly *early*-flat
  (at `ft = 0.23`), not fully-flushed-flat.
- `div*L/U` — non-dimensional divergence = **mass-conservation error**; should be
  ~0. `381` is large (a small-mesh run reached ~27). See §8.
- `cg_p`, `cg_uvw` — CG iterations for the pressure and velocity solves. High `cg_p`
  reflects the weak Jacobi preconditioner plus the consistency problem.

**Diagnostic lines** (when enabled) split a step:
`ENTRY` ($\mathbf{u}^n$) → `PRED` ($\mathbf{u}^*$) → `DIFF` ($\mathbf{u}^{**}$) → `POIS` ($\phi$) → `CORR` ($\mathbf{u}^{n+1}$), each
with its kinetic energy and divergence, so the stage that injects error is visible.

---

## 8. Where the method is fragile (and where the real work is)

**(a) Discrete consistency / checkerboard.** [weak] Equal-order P1–P1
velocity–pressure on an unstructured collocated layout is inf-sup (LBB) unstable →
admits a spurious checkerboard pressure mode. DDT projects exactly onto the discrete
divergence-free space. The *nodal* `div*L/U` print is a conservative upper bound, not
roundoff: it reads a few units even on a converged run (a developed 4-GPU run settles
to `div*L/U≈2.7`; early/under-developed runs read much higher, e.g. ~380 at `ft=0.2`
— that is transient, not the steady value). Whether the residual grows on *very*
large meshes is still open; the standard cures are:
- **Rhie–Chow momentum interpolation** — *tried; the compact $|\mathbf{A}|^2/(\mathbf{A}\cdot d\mathbf{x})$ form is
  geometrically unsafe on skewed tets ($\mathbf{A}$ not aligned with the edge → blows up).
  Abandoned for tets.*
- **PSPG / Bochev–Dohrmann pressure stabilization** (PSPG = Pressure-Stabilizing
  Petrov–Galerkin) — the FEM-canonical, geometry-
  robust fix. *Implemented for hex only; a tet version is the real remaining work.*

This is the sharpest topic to probe the developers on: how do they keep divergence
low on large pump meshes?

**(b) Explicit-advection CFL.** [weak] BJ/upwind + EXT2 are conditionally stable; a
developing jet can cross the limit. Remedies: smaller `dt`, `--bdf1`, or skew. An
adaptive-CFL `dt` controller was attempted but it breaks BDF2's constant-`dt`
coefficients and was shelved.

**(c) Multi-rank correctness — RESOLVED (2026-06-06).** [verified] Multi-rank now
matches single-rank: a 4-GPU run was verified bit-identical to single-rank at the
predictor (`div(u*)rms = 1.3171936`, `KE* = 5.6565327e-9` on both 1 and 4 ranks),
and a developed 4-GPU run (Re=100, `--skew`, do-nothing outlet, 3000 steps) reaches
`ft≈3.5`, `uMax/U≈13`, `d(u_rms)≈1e-5` (converged), `div*L/U≈2.7`, with `cg_p`
converging to ~1e-9 each step.

History worth keeping (the bug was real, then fixed, then *regressed*, then fixed
again — and the regression had nothing to do with the partition):
- The original multi-rank inlet-resolve bug (Exodus vs SFC index space) was fixed on
  2026-06-02 and verified 1-rank ≡ 4-rank (commits `d1b324e`/`8bb0402`).
- It later appeared to break again. After a long mis-diagnosis (chasing ownership /
  orphan-node theories on the inlet BC), the `MARS_PUMP_BC_DEBUG` per-rank audit
  proved the **inlet BC was never the problem** (every inlet node resolved and was
  owned by exactly one rank; e.g. 14+5 = 19 = the global inlet count). The real
  regression was in the **shared `mars_ns_solver.hpp`**: the TGV multi-rank-periodic
  work (≈55 commits) altered the DDT pressure path and made `cg_p` diverge (`-2`) for
  the pump. Fix: the pump now includes a **forked, verified-good solver**
  `mars_ns_pump_solver.hpp` pinned to `8bb0402` (commit `7ea9978`); TGV keeps the
  current shared header. Reconciling the two back into one solver is deferred.
- **Lesson for readers:** rank-0-local diagnostic prints (`BC: 0 inlet nodes`,
  `0 masked DOFs on rank 0`) are *not* global truth — a rank simply may hold no
  inlet nodes. Trust the rank-count-invariant sums (`[bc-count check]`,
  `MARS_PUMP_BC_DEBUG` global owned-of-resolved) instead.

**(d) Under-resolution.** [weak] No turbulence model + transitional Re + tet mesh →
a converged laminar/transitional solution, not resolved physics. No mesh-convergence
study yet.

---

## 9. Current validated status

- **Single- AND multi-rank: correct, stable, converges to steady state.** Verified
  1-rank ≡ 4-rank at the predictor; a developed 4-GPU run reaches `ft≈3.5`,
  `uMax/U≈13`, `div*L/U≈2.7`, `d(u_rms)≈1e-5` (converged). Physical inlet→outlet jet.
- The developers' Barth–Jespersen scheme runs, enabling a direct comparison.
- GPU-native, matrix-free, cornerstone-octree partitioned (built for scale); the
  multi-rank path is now trusted.
- The pump runs from a **forked solver** (`mars_ns_pump_solver.hpp`, pinned to the
  verified-good `8bb0402`) to insulate it from in-flight changes to the shared
  `mars_ns_solver.hpp`. Reconciling the two is deferred.
- Open: tet pressure stabilization for divergence consistency on very large meshes
  (the *nodal* `div*L/U` diagnostic stays a few units even when the projection is
  exact — a conservative upper bound, not "not converged").

---

## 10. Mental model (one paragraph)

> MARS solves incompressible Navier–Stokes by **median-dual CVFEM on tets**,
> vertex-centered P1, via a **Chorin projection** whose pressure Poisson operator is
> the **matrix-free $D M^{-1} D^T$** — chosen because the naive Galerkin Laplacian is
> inconsistent with the SCS gradient on tets and blows up. Time stepping is
> **BDF2/EXT2**; advection is selectable (**skew** = KE-conserving/DNS, **Barth–
> Jespersen** = the developers' 2nd-order limited scheme). It is GPU-native,
> matrix-free, and cornerstone-octree-partitioned for billion-element scale. It runs
> **stable to a steady state on multiple pump meshes single-rank**, with a physical
> jet, **single- and multi-rank** (1-rank ≡ 4-rank verified). The remaining open
> issue is **discrete divergence consistency** on very large meshes (a tet-robust
> pressure stabilization — PSPG; Rhie–Chow failed on tets).

---

## 11. Sharpest questions for the developers

1. **Pressure stabilization** — how do you suppress the checkerboard / keep
   $\nabla\cdot\mathbf{u}$ small on large pump meshes (Rhie–Chow, PSPG, a div-consistent operator)?
2. **Gradient reconstruction** for Barth–Jespersen — Green–Gauss or least-squares?
3. **Operating Re** for your pump case (so we match the regime).
4. **Outlet condition** — fixed pressure, convective ($\partial \mathbf{u}/\partial t + U_c\,\partial \mathbf{u}/\partial n = 0$), or
   velocity-outflow?
5. **Steady or unsteady** at this Re — and if unsteady, do you time-average?
6. **Peak/inlet velocity ratio** you observe.
7. **Confirm** the diaphragm side-set identity (same-surface check).

---

## Appendix: glossary

- **CVFEM** — Control-Volume Finite Element Method; vertex-centered finite volume on
  an FE mesh.
- **SCS** — sub-control-surface; a facet of the median-dual control-volume boundary.
- **DDT** — the $D M^{-1} D^T$ matrix-free pressure operator ($D$ = divergence, $M$ = lumped
  mass, $D^T$ = discrete gradient).
- **EXT2** — 2nd-order explicit (Adams–Bashforth-like) extrapolation of the advection
  term, paired with BDF2.
- **ft** — flow-through time, $t / (L/U)$.
- **inf-sup / LBB** — the stability condition equal-order velocity–pressure spaces
  violate, producing the checkerboard mode that pressure stabilization fixes.
