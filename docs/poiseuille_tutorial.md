# Poiseuille Channel Flow with MARS — A Validation Tutorial

> The 1500-step validation passes on 1, 2 and 4 GPUs; see the
> [results and recipe](../tests/reference/poiseuille/planar_validation.md).

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

`tests/data/poiseuille/poiseuille_hex_14k_elem.e` (Exodus format, shipped with MARS):
14,751 hexahedra, 30,000 nodes. Reading it needs a MARS build with netCDF.

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
  walls, every node in the mesh becomes a Dirichlet node and the entire
  velocity field is frozen — nothing can ever flow. The solver marks only the
  y faces as walls; on the z faces it keeps w = 0 and leaves u and v free.

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
| face flux `F` | the volume flow `u·A` through one sub-control face, stored per face |
| divergence `D_F` | net flux `Σ F` out of a node's control volume (the mass check) |
| gradient `G` | how p varies across the control-volume faces (the pressure push) |
| Laplacian `K` | exchange between neighboring control volumes through the faces (viscous drag, and the pressure equation) |

Two consequences matter for this example. First, the SCS faces are
*interior* faces — boundary opening faces need the explicit source of
Part 3, or inflow is invisible. Second, equal-order velocity/pressure is
not naturally stable. The flux through a face computed from the plain
average of its two nodes, `A·(u_L + u_R)/2`, cannot see a pressure that
alternates from node to node (a checkerboard). The solver therefore adds a
Rhie–Chow term to every face flux:

```
F = A·(u_L + u_R)/2 − h [ (∇p·A)_f − A·(G p_L + G p_R)/2 ],   h = dt_eff/ρ
```

The bracket is the difference between the pressure gradient on the face
(from the shape functions of the hex) and the average of the two nodal
gradients. It is tiny for a smooth pressure and large for a checkerboard,
so it couples every pressure node to its neighbours.

### 2.4 Marching in time: the projection method

Each time step splits the physics into manageable pieces (Chorin
projection, BDF2 time accuracy):

1. **Predict**: move velocity by advection + the old pressure gradient
   (explicit — cheap, but ignores incompressibility).
2. **Diffuse**: apply viscosity implicitly — one linear solve per velocity
   component (implicit so large time steps stay stable).
3. **Project**: the face fluxes `F**` of the predicted field violate
   continuity; solve the Poisson equation `K φ = −(ρ/dt_eff)(D_F F** + openings)`
   for a pressure correction `φ`.
4. **Correct**: `F = F** − h (∇φ·A)` on every face, `u = u** − h G φ` at the
   free nodes, `p += φ`.

Step 4 applies to the fluxes the same face gradients from which `K` is
built, so `D_F F + openings = 0` holds to the solver tolerance on any number
of GPUs. The nodal velocity follows the fluxes up to the stabilization term.

The projection (step 3) is the heart and the cost: it is a global problem —
a flux imbalance at the inlet must be felt instantly at the outlet — which
is why its linear system is the hard one.

### 2.5 The linear solvers

Each implicit step is a sparse linear system `A x = b` with one row per
node. At 30k nodes (or 10⁹ on Alps) you never factor A — you iterate:

- **CG (conjugate gradient)** only needs matrix-vector products, which are
  perfectly GPU-shaped. You stop when the residual `|Ax−b|/|b|` drops below
  `--tol`.
- A **preconditioner** is a cheap approximate inverse applied each iteration.
  The simplest, Jacobi (divide by the diagonal), is not enough for the
  pressure: its iteration count grows with the mesh and the number of GPUs
  (thousands of iterations per step on this channel).
- **Algebraic multigrid (BoomerAMG from Hypre)** solves the error on a
  hierarchy of coarser problems it builds from the matrix itself. Its
  iteration count stays flat as the mesh and the GPU count grow. Here it takes
  about 18 iterations per step for the pressure and 10–13 per velocity component.

Both matrices are constant in time, so the solver assembles them once and
builds the multigrid hierarchies once:

| System | Matrix |
|---|---|
| velocity (each of u, v) | `M/dt_eff + nu K`, fixed velocity rows and columns removed |
| pressure | `K`, the rows of outlet nodes (p = 0) removed |

Each GPU assembles the matrix of its own elements, and Hypre adds the
contributions of the nodes that several GPUs share (the same assembly the
Taylor–Green tutorial explains for periodic points). At setup the solver
compares the assembled `K` with the matrix-free flux correction and stops if
they differ.

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

### The fix: the opening fluxes

After the interior divergence scatter (and its reverse-halo exchange), add the
missing flux of every inlet and outlet node:

```
div[i] += A_in[i] · u_prescribed[i] + A_out[i] · u[i]
```

`A_in` and `A_out` are the node's outward face areas on the inlet and outlet
planes (negative at the inlet, whose outward normal is −x). The inlet uses the
prescribed velocity; the outlet uses the computed one, so the outlet flux is
part of the operator and needs no rescaling. The skew-symmetric advection gets
the matching term: an opening face with flux `m` contributes `−m q / 2` to its
node.

The areas are built on the GPU from the element faces on each opening plane:
every face node gets the area of its sub-quad (node, edge midpoints, face
centre), a quarter of the face for these rectangles. Each rank adds the faces
of its own elements and the reverse halo completes the shared nodes, so the
areas sum exactly to `H*dz` on any number of ranks. In code:
`nsOpeningAreaKernel` and `nsOpeningFlux` in
`backend/distributed/unstructured/fem/mars_navier_stokes.hpp`.

---

## Part 4 — Build and run

The example needs a CUDA build with `MARS_ENABLE_HYPRE=ON`:

```bash
cmake --build . --target mars_poiseuille_flow --parallel 32
```

Run the validation case on one GPU (1500 steps; the time loop takes about 24 s on a GH200):

```bash
srun --account=<acct> --time=00:30:00 --nodes=1 --ntasks-per-node=1 --export=ALL \
  ./examples/distributed/unstructured/mars_poiseuille_flow \
  --mesh=/path/to/mars/tests/data/poiseuille/poiseuille_hex_14k_elem.e \
  --uinf=1 --nu=0.01 --dt=0.01 --num-steps=1500 --report-every=100 --check \
  --vtu-output=poiseuille --vtu-every=50
```

The same command runs on any number of GPUs (`--ntasks-per-node=4`). For
scaling runs, generate the channel on every rank instead of reading a mesh:
`--cells=NX,NY` meshes `[0,10] × [0,1] × [0,0.06]` with one cell in z.

| Flag | Meaning |
|------|---------|
| `--mesh=FILE` / `--cells=NX,NY` | read a mesh, or generate the channel |
| `--y-grading=S` | cluster the generated y spacing at the walls, S in [0,1) |
| `--uinf --rho --nu --dt --num-steps` | inflow velocity, density, viscosity, time step, steps |
| `--bdf1` | first-order time stepping (default BDF2) |
| `--tol --max-iter` | relative AMG-PCG tolerance (default 1e-10) and iteration cap |
| `--report-every=N` | progress line interval |
| `--vtu-output=PREFIX --vtu-every=N` | ParaView frames |
| `--check` | release gate: exit 1 unless every check of Part 6 passes |
| `--profile-x --profile-xtol` | move the validation plane (default 90% down the channel) |
| `--comparison-output=PREFIX` | u, v, w, p at full precision for comparisons with other codes |

---

## Part 5 — Reading the output

Setup prints the operator check:

```
Pressure operator: assembled vs matrix-free, max |difference| / max |Kx| = 1.145854e-15
```

Then one line per `--report-every` steps:

```
Step   1500  t=15.0000  |u|_M=8.4091374298e-01  continuity=6.224e-12  amg(u,v,p)=10/13/18
```

- `|u|_M` — mass-weighted L2 norm of the streamwise velocity. For this mesh
  (volume 0.6) the developed parabola gives about 0.84; watching it rise from
  the seeded uniform flow and level off *is* watching the parabola form.
- `continuity` — max of `|D_F F + openings| / M` over the nodes the projection
  constrains: the discrete mass balance of every control volume, at the solver
  tolerance.
- `amg(u,v,p)` — AMG-PCG iterations of the last step. A solve that misses its
  tolerance stops the run.

At the end:

```
[timing] ranks=1 nodes=30000 ... ms/step: total=15.903 predictor=0.054 viscous=7.212 pressure=8.594 ...
Poiseuille validation
  profile RMS at x=9.0000 +/- 0.2000: 4.551373e-04 (U_max=1.500000e+00, 1200 nodes)
  flux Q(x)/Q(inlet) at 25/50/75%: 1.0000 / 1.0000 / 1.0000
  -dp/dx from p: 1.2178e-01   from the u profile: 1.1992e-01   exact: 1.2000e-01
```

- The profile RMS compares the computed `u(y)` on the probe slab with the
  analytic parabola. **Pass: RMS < 6.0e-3** (the reference report's
  tolerance). The probe sits at 90% of the channel — past the entrance length,
  upstream of any outlet influence.
- The flux ratios are the volumetric flux through interior cross-sections,
  computed from the solved velocity, over the inlet flux. They cannot be faked
  by boundary values; 1.0000 at 25/50/75% means the flow carries all the mass
  through the channel.
- `-dp/dx` is measured twice: from the computed pressure between 60% and 90% of
  the channel, and from a parabola fit of the velocity core (`G = −μ u″`). Both
  should match the exact `12 ρ ν U / H²`.

With `--vtu-output=PREFIX` the run writes `PREFIX.pvd` with the point fields
`u`, `v` and `p`.

---

## Part 6 — Regression test, profile plot, and animation

**Regression test.** With `--check` the example grades itself at the end:

```
VALIDATION PASS: RMS=4.551e-04 < 6.000e-03, flux PASS, steady=2.627e-08 PASS, continuity*H/U=8.050e-15 balance=8.797e-08 PASS, projection PASS
```

It requires the profile RMS and the flux ratios above, plus:

- **continuity**: the RMS of `D_F F + openings` per unit volume, times H/U, at
  most 1e-6, and the inlet and outlet fluxes balancing to 1e-6 (the outlet flux
  uses the nodal velocity, so the balance shows the stabilization term, about 1e-7);
- **steadiness**: the velocity change over the last 20 steps, divided by U,
  at most 1e-6;
- **projection**: at steps 1, 2, 3 and the last, the corrected fluxes
  satisfy `D_F F = D_F F** + (dt_eff/ρ) K φ` to 1e-7 (printed as
  `[channel-projection]` lines).

The case is registered with ctest as `marsPoiseuilleValidation` (labels
`validation;gpu;long`) when you configure with
`-DMARS_ENABLE_VALIDATION_TESTS=ON`. Run it with `ctest -L validation` inside
a GPU allocation. It is the canary for any change to the projection, the
boundary conditions or the opening fluxes.

**Profile plot.** Every run writes `PREFIX_profile.csv` (default
`poiseuille_profile.csv`): the computed u(y) on the one node plane nearest the
probe. Plot it against the exact curve:

```bash
python3 scripts/plot_poiseuille_profile.py poiseuille_profile.csv
```

**The exact solution, four ways.** The Wikipedia "Plane Poiseuille flow"
section defines the solution by four related quantities; the run checks all of
them independently:

| Wikipedia quantity | Checked by |
|---|---|
| `u(y) = G/(2μ)·y(h−y)` | profile RMS at fixed x + the figure above |
| `Q = Gh³/(12μ)` | flux ratios at 25/50/75% |
| `G = −dp/dx` | the two `-dp/dx` values |
| `U_max = Gh²/(8μ)` | the `U_max = 1.5·U_mean` target (same parabola, anchored by the inlet) |

**Animation.** Render the channel developing — uniform plug at the inlet
bending into the parabola:

```bash
pvbatch scripts/render_poiseuille.py --pvd PREFIX     # -> PREFIX_movie.mp4
```

Use `--vtu-every=10` for a smooth movie.

---

## Part 7 — Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `NavierStokes: planar flow needs one layer of elements between two z planes` | the example solves planar flow (u, v, p) on a one-element-thick mesh | use a mesh with every node on one of two z planes |
| `the assembled pressure operator is wrong` at setup | assembly and matrix-free operator disagree | a real bug: report it with the rank count and mesh |
| `the ... solve did not converge` | `--max-iter` too low, or a broken input (NaN) | check the last `amg(u,v,p)` line; raise `--max-iter` |
| identical numbers after a code change | stale binary | rebuild `mars_poiseuille_flow` |
| more GPUs are slower on the tutorial mesh | 30k nodes are too few to share | scale with `--cells=NX,NY` |

---

## Part 8 — Where to go next

- `examples/distributed/unstructured/mars_poiseuille_flow.cu` — the example:
  mesh and domain, solver, time loop, output, in that order.
- `backend/distributed/unstructured/fem/mars_navier_stokes.hpp` — the solver,
  shared with the Taylor–Green vortex and the lid-driven cavity: boundary
  conditions, DOF numbering for Hypre, assembly, the four stages of a step, and
  the projection check.
- [periodic_tgv_tutorial.md](periodic_tgv_tutorial.md) — the same solver on a
  periodic box, and how the unknowns are shared between GPUs.
- `examples/distributed/unstructured/mars_poiseuille_validation.hpp` — the
  validation against the exact solution and the `--check` gate.
- [CVFEM-Kernels.md](CVFEM-Kernels.md) and [FEM-Assembly.md](FEM-Assembly.md) —
  the discrete CVFEM operators.
- `report.html` — the FLUYA reference report this case is validated against,
  including the reference solver's input file.
