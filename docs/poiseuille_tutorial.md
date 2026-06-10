# Poiseuille Channel Flow with MARS — A Validation Tutorial

This tutorial explains the `mars_poiseuille_flow` example: what Poiseuille flow
is, why it is the standard first validation case for any incompressible CFD
code, how to build and run it, how to read every line of its output, and the
one numerical lesson it teaches that generalizes to every inlet-driven flow in
MARS (including the pump). It assumes no prior CFD knowledge at the start; the
deeper companion is `docs/pump_cfd_tutorial.md`, which builds the full
numerical machinery from scratch.

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

## Part 2 — How the solver works (the short version)

Each time step of the Chorin projection method does:

1. **Predict**: move the velocity by advection + old pressure gradient
   (explicit, BDF2/EXT2).
2. **Diffuse**: implicit viscous solve, one CG solve per velocity component.
3. **Project**: compute how much the predicted field violates `div u = 0`,
   solve a pressure Poisson equation, and correct the velocity so mass is
   conserved.
4. **Correct**: update velocity and pressure; re-impose boundary values.

The pressure projection is the heart of incompressible CFD: pressure is
whatever it must be to keep the velocity divergence-free. The full derivation
is in `docs/pump_cfd_tutorial.md`. Here the only step that needs detail is
the projection — because it is where this example's big lesson lives.

---

## Part 3 — The lesson: the opening-flux source

### The bug this example exposed

The discrete divergence in MARS is assembled from **interior** sub-control-
surface faces only: every internal face between two nodes contributes
`+flux` to one node and `-flux` to the other, and the sum telescopes — except
at the domain boundary, where the *opening faces* (inlet, outlet) are never
integrated.

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

---

## Part 5 — Reading the output

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
| `|u|` plateau | 0.846 (analytic 0.849) |
| RMS vs parabola | 5.8e-3 raw, 3.9e-3 normalized (tol 6e-3) |
| flux ratios 25/50/75% | 0.99 – 1.01 |
| `div_max` at t=12 | ~1.6, still falling |
| wall time, 1200 steps, 1 GH200 | ~20 min |

---

## Part 6 — Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `|u|` frozen ~0.045, `div_max` constant | opening-flux source off, or all nodes Dirichlet (z faces marked as walls) | default config already handles both; check the `Drive:` and `Opening-flux source:` banner lines |
| `cg_iter_p=FAIL`, residual stalls ~1e-2..1e-5 | tolerance unreachable for Jacobi-PCG, or `--max-iter` too low | `--tol=1e-6 --max-iter=4000` |
| `cg_iter_p=FAIL` at iteration 1 (`pAp<=0`) | Jacobi diagonal clip too small for extreme cell-size ratios | `MARS_DDT_JACOBI_CLIP_FRAC=1e-2` (env var, no rebuild) |
| blows up within a few steps after flow develops | unbalanced opening source (should not happen with `oScale`), or — on harder meshes — the equal-order checkerboard instability | this clean hex mesh does not checkerboard through t=12; for meshes that do, pressure stabilization (PSPG/VMS) is required — active work on the pump side |
| run dies near step 1250 with no validation block | srun wall-time limit | `--time=00:30:00` |

---

## Part 7 — Where to go next

- `docs/pump_cfd_tutorial.md` — the full from-scratch CFD background: CVFEM,
  the discrete operators, the projection algebra, stabilization.
- `examples/distributed/unstructured/mars_poiseuille_flow.cu` — the driver:
  BC rebuild, area lumping, validation, flux probe. Deliberately small.
- `backend/distributed/unstructured/fem/mars_ns_channel_solver.hpp` — the
  channel fork; its only delta vs the shared solver is the opening-flux
  source (search "Balanced opening-flux source").
- `report.html` — the FLUYA reference report this case is validated against,
  including the reference solver's input file.
