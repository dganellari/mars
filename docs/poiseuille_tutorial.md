# Poiseuille Channel Flow with MARS — A Validation Tutorial

This tutorial explains the `mars_poiseuille_flow` example: what Poiseuille flow
is, why it is the standard first validation case for any incompressible CFD
code, how to build and run it, how to read every line of its output, and the
one numerical lesson it teaches that generalizes to every inlet-driven flow in
MARS (including the pump). It assumes no prior CFD knowledge at the start and
builds the numerical machinery as it goes.

---

## Part 0 — What is Poiseuille flow, and why validate with it?

Push water through a long straight channel between two parallel walls. The
fluid sticks to the walls (the *no-slip* condition: velocity is exactly zero
at a wall), so the layers near the walls are slow and the middle is fast.
Far enough downstream, the velocity settles into a profile that no longer
changes — and that profile is exactly a **parabola**:

```
u(y) = U_max * (1 - (2y/H)^2)
```

where `H` is the channel height, `y` is measured from the centerline, and
`U_max` is the speed in the middle. This is Hagen–Poiseuille flow (the plane-
channel variant), one of the very few solutions of the Navier–Stokes
equations that can be written in closed form.

Two facts about the parabola matter for validation:

1. **`U_max = 1.5 × U_mean`.** If you push fluid in uniformly at speed `U`,
   mass conservation forces the developed centerline speed to be exactly
   `1.5 U`. No fitting, no tuning — the number is fixed by the math.
2. **It is reached after a known distance.** A uniform inflow needs roughly
   the *entrance length* `Le ≈ 0.05 * Re * H` to relax into the parabola.
   At Re = 100 and H = 1 that is `Le ≈ 5`, so in our 10-unit channel the
   profile is fully developed well before the outlet.

Because the exact answer is known, this case can fail loudly. A CFD code that
produces a parabola with the right `U_max` at the right distance has
demonstrated, in one run: the advection and diffusion operators, the pressure
projection, the inlet/outlet/wall boundary conditions, and global mass
conservation. That is why every CFD code's regression suite starts here, and
why we compare against a reference solution (`report.html`, produced by FLUYA,
the STK-based CVFEM code MARS mirrors) with a hard tolerance: **RMS error
below 6.0e-3**.

---

## Part 1 — The case set-up

### The mesh

`poiseuille_hex_14k_elem.e` (Exodus format): 14,751 hexahedra, 30,000 nodes.

```
x: -0.5 .. 10.5   streamwise   (~150 node planes)
y:  0.0 .. 1.0    wall-normal  (~100 node planes)  -> channel height H = 1
z:  0.0 .. 0.066  ONE element thick                (2 node planes)
```

The z-direction is a single element: this is a **quasi-2D** mesh. Plane
Poiseuille flow is two-dimensional, and the thin z-extrusion exists only
because the solver works in 3D. This has one important consequence (Part 3).

### The physics parameters

```
density          rho = 1
kinematic visc.  nu  = 0.01
inflow speed     U   = 1
Reynolds number  Re  = U*H/nu = 100   (laminar -- the parabola is stable)
```

### The boundary conditions (the FLUYA reference configuration)

| Boundary       | Condition                                            |
|----------------|------------------------------------------------------|
| inlet (x=xmin) | velocity Dirichlet u=(1,0,0) — uniform plug inflow   |
| outlet (x=xmax)| pressure Dirichlet p=0, velocity **free**            |
| walls (y=0, 1) | no-slip u=(0,0,0)                                    |
| z faces        | symmetry (free) — NOT walls                          |

Two of these are easy to get wrong:

- **The outlet is a pressure condition, not a velocity condition.** The
  outlet velocity must stay free so the parabola can leave the domain
  undisturbed; prescribing it over-constrains the exit and the interior flow
  dies. Pinning p=0 on the whole outlet plane anchors the pressure level and
  lets mass balance itself through the pressure field.
- **The thin z faces are symmetry planes, not walls.** On a one-element-thick
  mesh every single node touches a z face. If z faces are treated as no-slip
  walls (the default 3D channel marking), every node in the mesh becomes a
  Dirichlet node and the entire velocity field is frozen — nothing can ever
  flow. The driver rebuilds the boundary masks so only the y walls are
  no-slip.

---

## Part 2 — Equations, exact solution, discretization, solvers (for starters)

### 2.1 The equations

Incompressible flow obeys the Navier–Stokes equations:

```
∂u/∂t + (u·∇)u = −(1/ρ)∇p + ν∇²u      (momentum)
∇·u = 0                                (mass: fluid neither piles up nor vanishes)
```

In words, a fluid parcel accelerates because neighboring fluid carries it
along (*advection*, the `(u·∇)u` term), because pressure differences push it
(`∇p`), and because viscous friction drags it (`ν∇²u`). The second equation
is the hard one: there is no equation "for" pressure — pressure is whatever
it must be so that the velocity stays divergence-free. Every incompressible
solver is, at heart, a strategy for finding that pressure.

### 2.2 Why the exact solution is a parabola (3-line derivation)

Far from the inlet the flow is *steady* (`∂u/∂t = 0`) and *fully developed*
(`u = (u(y), 0, 0)`, nothing changes with x). Then advection vanishes
(`u·∂u/∂x = 0`) and the x-momentum equation collapses to a balance between
the pressure push and the wall drag:

```
0 = −dp/dx + μ d²u/dy²
```

The left term cannot depend on y, the right cannot depend on x, so both equal
a constant: `dp/dx = −G`. Integrate twice and apply no-slip (`u = 0` at both
walls):

```
u(y) = G/(2μ) · y(h−y)        — the parabola
```

Everything else follows by integration: flow rate `Q = Gh³/(12μ)`, mean
velocity `U_mean = Gh²/(12μ)`, centerline `U_max = Gh²/(8μ) = 1.5·U_mean`.
This is why the case validates so sharply — every number is pinned by the
input parameters alone.

### 2.3 The discretization (CVFEM in one breath)

The channel is tiled into hexahedral **elements** with **nodes** at the
corners; velocity *and* pressure live at the nodes (equal-order,
vertex-centered — one set of points for everything). Around every node sits
a small **control volume**, and the physics is enforced in integral form:
whatever flows into a control volume must flow out. The "FEM" in CVFEM is
how the fluxes through the control-volume faces are evaluated: with finite-
element shape functions interpolated inside each hex, on the sub-control
surfaces (SCS) that the element contributes.

From this one idea come the discrete operators the solver uses:

| Operator | Discrete meaning |
|---|---|
| divergence `D` | net flux `Σ u·A` out of a node's control volume (the mass check) |
| gradient `G` | how p varies across the control-volume faces (the pressure push) |
| Laplacian `K` | viscous exchange between neighboring control volumes (the drag) |

Two consequences matter for this example. First, the SCS faces are
*interior* faces — boundary opening faces need the explicit source of
Part 3, or inflow is invisible. Second, equal-order velocity/pressure is
not naturally stable (the LBB issue): it works here on a clean mesh, but
harder cases (the pump) need pressure stabilization.

### 2.4 Marching in time: the projection method

Each time step splits the physics into manageable pieces (Chorin
projection, BDF2 time accuracy):

1. **Predict**: move velocity by advection + the old pressure gradient
   (explicit — cheap, but ignores incompressibility).
2. **Diffuse**: apply viscosity implicitly — one linear solve per velocity
   component (implicit so large time steps stay stable).
3. **Project**: the predicted field violates `∇·u = 0`; solve a Poisson
   equation for a pressure correction `phi` whose gradient removes exactly
   that violation.
4. **Correct**: subtract the gradient, update `p += phi`, re-impose boundary
   values.

The projection (step 3) is the heart and the cost: it is a global problem —
a flux imbalance at the inlet must be felt instantly at the outlet — which
is why its linear system is the hard one.

### 2.5 The linear solvers

Each implicit step is a sparse linear system `A x = b` with one row per
node. At 30k nodes (or 10⁹ on Alps) you never factor A — you iterate:

- **CG (conjugate gradient)** only needs matrix-vector products, which are
  perfectly GPU-shaped. Each iteration improves the answer; you stop when
  the residual `|Ax−b|/|b|` drops below `--tol`.
- A **preconditioner** is a cheap approximate inverse applied each iteration
  to speed convergence. Here: **Jacobi** (divide by the diagonal) — the
  simplest possible one.
- The three **velocity** solves are easy (mass-dominated matrices, a dozen
  iterations). The **pressure Poisson** solve is hard: its conditioning
  worsens with mesh size and cell-aspect-ratio, and on this thin-z channel
  Jacobi-PCG needs ~3600 iterations per step. That is the entire runtime.
  Stronger preconditioners (multigrid/AMG) exist, but reject this particular
  matrix-free operator — a known open item, shared with the wing/pump work.
- The pressure solve here is **matrix-free**: the operator `D M⁻¹ Dᵀ` is
  applied as three scatters per CG iteration instead of being assembled —
  cheaper in memory, and exactly consistent with the divergence/gradient
  used elsewhere in the step.

Deeper material on the CVFEM operators lives in
[CVFEM-Kernels.md](CVFEM-Kernels.md) and [FEM-Assembly.md](FEM-Assembly.md).

---

## Part 3 — The lesson: the opening-flux source

### The bug this example exposed

The discrete divergence in MARS is assembled from **interior** sub-control-
surface faces only: every internal face between two nodes contributes
`+flux` to one node and `-flux` to the other, so by **summation by parts** the
interior contributions cancel in pairs and only the boundary survives — except
the *opening faces* (inlet, outlet) at the domain boundary are never integrated.

The consequence is invisible until you run an inlet-driven flow: the
prescribed inlet velocity **never enters the divergence bookkeeping**, so the
pressure solve never "sees" any fluid entering. It builds no pressure
gradient down the channel, the inlet value sits painted on the first node
plane, and nothing ever flows. Symptoms (all observed before the fix):

- `|u|` frozen at the inlet-sliver value, forever — no downstream propagation
  even after many flow-through times;
- `div_max` pinned at a constant (the un-cancellable one-sided boundary
  residual);
- with a natural outlet, the pressure CG *stalls* at a residual floor and
  `div_max` scales like `1/dt` — the fingerprint of a right-hand side the
  operator cannot represent (smaller dt makes it *worse*, ruling out a
  time-step problem).

An analogy: a warehouse inventory system that logs every pallet moved between
aisles but has no scanner at the receiving dock. Trucks unload all day, the
system shows zero incoming stock, so it never schedules anything to ship out.

This is the same root cause that left the pump's passage dead — the fix below
was developed for the pump and is validated here against the analytic answer.

### The fix: a balanced opening-flux source

After the interior divergence scatter (and its reverse-halo exchange) and
before the pressure RHS is built, add the missing opening fluxes directly to
the divergence accumulator:

```
divAccNode[i] += U_prescribed . areaVec_outward[i]     (inlet and outlet nodes)
```

with two non-negotiable details, both learned the hard way:

1. **Use prescribed velocities, never the solved field.** At step 0 the
   solved outlet velocity is zero, so a solved-field source is one-sided —
   the system becomes inconsistent and blows up immediately.
2. **Balance the source to machine zero.** The inlet and outlet
   contributions must cancel *exactly* — to the last floating-point bit. The
   RHS multiplies the divergence by `rho/dt` (here 100), so even a tiny
   numerical imbalance is amplified every step and explodes the solve. The
   recipe (OpenFOAM's `adjustPhi` idea): measure both discrete fluxes, then
   rescale the outlet by `oScale = -Qin/Qout` so the sum is zero by
   construction.

In code: `addOpeningFluxSourceKernel` + the `oScale` block in
`runPressureSolveStep` of **`mars_ns_channel_solver.hpp`** — a fork of the
shared `mars_ns_solver.hpp` (same pattern as the pump's fork), so the shared
solver used by cavity/channel/TGV is untouched. The source is gated behind
`NSStepper::useOpeningFluxSource`, default off; the driver enables it.

The per-node outward area vectors are built by the driver: for each opening
plane node, area = (Voronoi interval in y) × (dz/2), which sums exactly to
the face area `H*dz`. For domain-aligned planes the direction is just `-x`
(inlet) and `+x` (outlet).

---

## Part 4 — Build and run

Build (the target links only against the `mars` library):

```bash
make mars_poiseuille_flow -j
```

Run (single rank; the area lumping assumes the opening planes are rank-local):

```bash
MARS_NODEHALO_V2=1 srun --account=<acct> --time=00:30:00 \
  --nodes=1 --ntasks-per-node=1 \
  ./examples/distributed/unstructured/mars_poiseuille_flow \
  --mesh=/path/to/poiseuille_hex_14k_elem.e \
  --uinf=1.0 --nu=0.01 --dt=0.01 --tol=1e-6 --max-iter=4000 --num-steps=1200 \
  --vtu-output=poiseuille
```

Notes on the numbers:
- `--num-steps=1200` (t = 12) is comfortably past convergence (~t = 10).
- `--max-iter=4000`: the matrix-free DDT pressure solve with Jacobi
  preconditioning needs ~3600 CG iterations per step on this mesh. That is
  the price of the un-preconditioned thin-channel Poisson operator (AMG
  rejects it); it is slow but completely stable.
- `--tol=1e-6` for the pressure solve is sufficient; 1e-10 is unreachable
  for Jacobi-PCG here and would FAIL every step.
- Budget ~25 minutes wall; a 15-minute limit dies around step 1250.

Useful flags:

| Flag | Meaning |
|------|---------|
| `--drive=inlet` (default) | FLUYA reference config: velocity inlet + free pressure outlet |
| `--drive=bodyforce`       | constant streamwise force instead (study mode; needs periodic x to be meaningful) |
| `--no-opening-flux-source`| disable the Part-3 source — reproduces the dead channel (A/B baseline) |
| `--no-seed-interior`      | start from rest instead of seeded plug flow (slower development) |
| `--profile-x=X`           | move the validation plane (default 90% down the channel) |
| `--cross-axis=y\|z`       | which axis the parabola varies over (default y) |
| `--check`                 | regression mode: exit 1 unless RMS and flux ratios pass (Part 6) |
| `--rms-tol=VAL`           | RMS pass threshold for `--check` (default 6e-3, the reference tol) |
| `--flux-tol=VAL`          | flux-ratio tolerance for `--check` (default 0.10 = ±10%) |

---

## Part 5 — Reading the output

With `--vtu-output=PREFIX` the run writes three ParaView timelines:

| Files | Field | Use |
|---|---|---|
| `PREFIX_u.pvd` + steps | u — streamwise velocity component | **the one to visualize** (clean, validated) |
| `PREFIX_umag.pvd` | umag — velocity *magnitude* √(u²+v²+w²) | avoid: contaminated by the v/w artifact (Part 6) |
| `PREFIX_p.pvd` | p — pressure | avoid: carries the accumulated startup ramp (Part 6) |

Per-step line:

```
Step  500: t=5.0000 |u|=8.167e-01 div_max=2.947e+00 cg_iter_p=3630
```

- `|u|` — mass-weighted L2 norm of streamwise velocity over the domain. For
  this mesh (volume 0.6) the developed parabola gives **|u| → 0.849**;
  watching it rise from the seeded 0.77 and plateau there *is* watching the
  parabola form. If it is frozen near 0.045 the channel is dead (Part 3).
- `div_max` — max nodal divergence after projection. It should fall steadily
  (7.4 → 1.6 over the run). The nonzero floor lives at boundary nodes where
  the one-sided stencil cannot vanish; the interior is much cleaner.
- `cg_iter_p` — pressure CG iterations. ~3600 is normal here. `FAIL` means
  it hit `--max-iter` or a breakdown — see Part 6.

Final validation block:

```
  probe plane x = 9.0000 ...
  U_max analytic = 1.5000  (= 1.5*Uinf)
  RMS error      = 5.781047e-03  (normalized: 3.854031e-03)
```

The RMS compares solved `u(y)` on the probe plane against the analytic
parabola. **Pass criterion: RMS < 6.0e-3** (the reference report's
tolerance). The probe sits at 90% of the channel — past the entrance length,
upstream of any outlet influence.

Interior-flux probe:

```
  Q(inlet) = 6.0000e-02
  Q(25%) = 5.9529e-02  ratio=0.992
  Q(50%) = 5.9529e-02  ratio=0.992
  Q(75%) = 6.0377e-02  ratio=1.006
```

This is the *honest* through-flow diagnostic: the volumetric flux through
interior cross-sections, computed from the **solved** velocity. It cannot be
faked by boundary values (a flux measured at a Dirichlet plane just reports
the BC back at you — a lesson from the pump debugging). Ratios near 1.0 at
25/50/75% mean real, mass-conserving flow through the whole channel.

### Expected results (the validated run)

| Quantity | Expected |
|----------|----------|
| `|u|` plateau | 0.847 (analytic 0.849) |
| RMS vs parabola | 5.4e-3 raw, 3.6e-3 normalized (tol 6e-3) |
| flux ratios 25/50/75% | 0.99 – 1.01 |
| velocity-fit G | 0.121 (exact 0.12, ~101%) |
| `div_max` at convergence | ~1.5–1.8 (frozen boundary-bookkeeping residual, not interior error) |
| wall time, 1500 steps, 1 GH200 | ~25 min |

---

## Part 6 — Regression test, profile plot, and animation

**Regression test.** With `--check` the driver grades itself at the end and
sets the exit code:

```
VALIDATION PASS: RMS=5.781e-03 < 6e-3, flux ratios 0.992/0.992/1.006 within 1 +/- 0.1
```

The case is registered with ctest as `marsPoiseuilleValidation` (labels
`validation;gpu;long`, 40-minute timeout) whenever the mesh sits at the repo
root. Run it with `ctest -L validation`. This is the canary for any change to
the projection, the boundary conditions, or the opening-flux source: if one
of them regresses, the parabola degrades and the test fails loudly.

**Profile plot.** At the end of every run the driver writes
`<prefix>_profile.csv` — solved u(y) sampled at ONE fixed x station (the node
plane nearest the probe x; the RMS uses a small slab, but a figure must show a
single station). Turn it into the same figure the reference report shows:

```bash
python3 scripts/plot_poiseuille_profile.py poiseuille_profile.csv
```

(plain matplotlib, no VTK). The exact curve is drawn from the closed-form
Wikipedia formula — `u(y) = G/(2μ)·y(h−y)` with `G = 8μU_max/h²` stated in the
legend — not fitted to the data, and the script warns if the driver's analytic
column ever disagrees with it.

**The exact solution, four ways.** The Wikipedia "Plane Poiseuille flow"
section defines the solution by four related quantities; the run checks all of
them independently:

| Wikipedia quantity | Checked by |
|---|---|
| `u(y) = G/(2μ)·y(h−y)` | profile RMS at fixed x + the figure above |
| `Q = Gh³/(12μ)` | interior-flux probe (Q ≈ Q_in at 25/50/75%) |
| `G = −dp/dx` | **velocity-fit G**: quadratic fit of the core profile, `G = −μ·u″` (validated: 100.8% of exact at convergence) |
| `U_max = Gh²/(8μ)` | identical to the `U_max = 1.5·U_mean` target (same parabola, anchored by the inlet velocity instead of G) |

The G line in the run log reads:

```
G from velocity (core parabola fit, 48 nodes): 1.2096e-01   exact 1.2000e-01   (100.8% of exact)
```

A match here means the velocity profile, the flow rate, *and* the momentum
balance (which sets G) independently agree with the exact solution — stronger
than the profile RMS alone. Note the fit tracks the actual state: on a
half-developed run it honestly reports ~96%, reaching ~100% only at
convergence.

**Why G is measured from the velocity, not the solved pressure (a known
artifact).** Incremental Chorin updates `p += phi` every step and nothing ever
removes what accumulates, so the violent startup transient leaves a giant
frozen smooth ramp in `p` (~1e5 × the physical pressure; the driver prints the
raw `dp/dx` under an ARTIFACT label for tracking). Worse, the mismatch between
the predictor's SCS gradient and the corrector's `D^T` gradient leaks this
giant field into the weakly-anchored `v`/`w` components — `umag` in the VTU
output reaches O(100) garbage while `u` underneath is the clean validated
parabola. Consequences: **visualize the `u` field, never `umag` or `p`** (the
render script defaults to `u` for this reason), and trust only velocity-derived
diagnostics. The clean structural fix is a boundary-complete divergence /
stabilized formulation — active work on the pump side.

**Animation.** `--vtu-output=PREFIX` writes a `PREFIX_umag.pvd` timeline;
render the channel developing — uniform plug at the inlet bending into the
red-core/blue-wall parabola — with:

```bash
pvbatch scripts/render_poiseuille.py --pvd PREFIX     # -> PREFIX_movie.mp4
```

Use `--vtu-every=10` on the run for a smooth 121-frame movie (default 50
gives 25 frames).

---

## Part 7 — Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `|u|` frozen ~0.045, `div_max` constant | opening-flux source off, or all nodes Dirichlet (z faces marked as walls) | default config already handles both; check the `Drive:` and `Opening-flux source:` banner lines |
| `cg_iter_p=FAIL`, residual stalls ~1e-2..1e-5 | tolerance unreachable for Jacobi-PCG, or `--max-iter` too low | `--tol=1e-6 --max-iter=4000` |
| `cg_iter_p=FAIL` at iteration 1 (`pAp<=0`) | Jacobi diagonal clip too small for extreme cell-size ratios | `MARS_DDT_JACOBI_CLIP_FRAC=1e-2` (env var, no rebuild) |
| blows up within a few steps after flow develops | unbalanced opening source (should not happen with `oScale`), or — on harder meshes — the equal-order checkerboard instability | this clean hex mesh does not checkerboard through t=12; for meshes that do, pressure stabilization (PSPG/VMS) is required — active work on the pump side |
| run dies near step 1250 with no validation block | srun wall-time limit | `--time=00:30:00` |
| `umag`/`p` look like garbage in ParaView while the run PASSes | known incremental-pressure artifact leaking into v/w (Part 6) | render the `u` field (script default); trust velocity diagnostics only |
| rerun reproduces identical numbers to many digits after a code change | stale binary — build dir not rebuilt | rebuild `mars_poiseuille_flow`; bit-identical output after an intended change is the tell |

---

## Part 8 — Where to go next

- [CVFEM-Kernels.md](CVFEM-Kernels.md) and [FEM-Assembly.md](FEM-Assembly.md) —
  the discrete CVFEM operators and how they assemble into the projection method.
- `examples/distributed/unstructured/mars_poiseuille_flow.cu` — the driver:
  BC rebuild, area lumping, validation, flux probe. Deliberately small.
- `backend/distributed/unstructured/fem/mars_ns_channel_solver.hpp` — the
  channel fork; its only delta vs the shared solver is the opening-flux
  source (search "Balanced opening-flux source").
- `report.html` — the FLUYA reference report this case is validated against,
  including the reference solver's input file.
