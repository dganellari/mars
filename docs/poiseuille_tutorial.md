# Poiseuille Channel Flow

This tutorial runs the `mars_poiseuille_flow` example: incompressible flow through a plane
channel, checked against the exact solution. It explains the case, the numerical method of the
MARS Navier–Stokes solver, how to run the example, and how to read every line it prints. No
prior CFD knowledge is assumed.

> The 1500-step validation passes on 1, 2 and 4 GPUs; see the
> [validation record](https://github.com/dganellari/mars/blob/master/tests/reference/poiseuille/planar_validation.md).

---

## 1. What Poiseuille flow is, and why it is a validation case

Push fluid through a long channel between two parallel walls. The fluid sticks to the walls
(the *no-slip* condition: the velocity is zero at a wall), so the layers near the walls are slow
and the middle is fast. Far enough downstream, the velocity profile stops changing, and it is a
**parabola**:

```
u(y) = U_max (1 - (2y/H)^2)
```

`H` is the channel height, `y` is measured from the centre line, and `U_max` is the speed in the
middle. This is plane Poiseuille flow, one of the few exact solutions of the Navier–Stokes
equations.

Two facts matter for the validation:

1. **`U_max = 1.5 U`** for a uniform inflow at speed `U`. Mass conservation fixes this number.
2. **The parabola forms after the entrance length** `Le ≈ 0.05 Re H`. At `Re = 100` and `H = 1`
   that is about 5, so in the 10-long channel of this example the flow is developed well before
   the outlet.

Because the answer is known, the case can fail clearly. A correct parabola at the right place
tests, in one run, advection and diffusion, the pressure projection, the inlet, outlet and wall
conditions, and global mass conservation.

---

## 2. The case

### The mesh

`tests/data/poiseuille/poiseuille_hex_14k_elem.e` (Exodus II, shipped with MARS; reading it needs a
MARS build with netCDF): 14,751 hexahedra and 30,000 nodes.

```
x: 0 .. 10     streamwise     (150 node planes)
y: 0 .. 1      wall-normal    (100 node planes)  -> channel height H = 1
z: 0 .. 0.06   one element    (2 node planes)
```

The mesh is one element thick in z: plane Poiseuille flow is two-dimensional, and the thin
extrusion exists only because the solver works in 3D. The solver runs in *planar* mode on such a
mesh: it solves for `u`, `v` and `p` and keeps `w = 0`.

### The physical parameters

```
density              rho = 1
kinematic viscosity  nu  = 0.01
inflow speed         U   = 1
Reynolds number      Re  = U H / nu = 100   (laminar)
```

### The boundary conditions

| Boundary | Condition |
|---|---|
| inlet (x = xmin) | velocity fixed, u = (1, 0, 0): uniform inflow |
| outlet (x = xmax) | pressure fixed, p = 0; velocity free |
| walls (y = ymin, ymax) | no-slip, u = (0, 0, 0) |
| z faces | nothing: planar mode keeps w = 0 |

The example finds these boundaries from the node coordinates, so the mesh needs no named side
sets. Two of them are easy to get wrong:

- **The outlet is a pressure condition.** The outlet velocity must stay free so that the
  parabola can leave the domain. Fixing it over-constrains the exit. `p = 0` on the outlet plane
  fixes the pressure level.
- **The z faces are not walls.** On a one-element-thick mesh every node lies on a z face. If the
  z faces were no-slip walls, every node would have a fixed velocity and nothing could flow.

---

## 3. The method

### 3.1 The equations

```
du/dt + (u . grad) u = -(1/rho) grad p + nu laplace(u)    (momentum)
div u = 0                                                  (mass)
```

A fluid parcel accelerates because the flow carries it (advection, `(u . grad) u`), because
pressure differences push it (`grad p`), and because viscous friction slows it
(`nu laplace(u)`). There is no equation "for" the pressure: the pressure is whatever keeps the
velocity divergence-free. Every incompressible solver is a strategy to find that pressure.

### 3.2 Why the exact solution is a parabola

Far from the inlet the flow is steady and does not change with x: `u = (u(y), 0, 0)`. Advection
then vanishes, and the x-momentum equation is a balance between the pressure push and the
viscous drag (`mu = rho nu`):

```
0 = -dp/dx + mu d2u/dy2
```

The first term cannot depend on y and the second cannot depend on x, so both are a constant:
`dp/dx = -G`. Integrate twice and apply `u = 0` at both walls (`y` now measured from the lower
wall):

```
u(y) = G / (2 mu) * y (H - y)
```

Then the flow rate is `Q = G H^3 / (12 mu)`, the mean velocity `U = G H^2 / (12 mu)`, and the
centre-line speed `U_max = G H^2 / (8 mu) = 1.5 U`. For this case `G = 12 mu U / H^2 = 0.12`.

### 3.3 The discretization (CVFEM)

The channel is split into hexahedral **elements** with **nodes** at their corners. Velocity and
pressure live at the same nodes (equal order). Around every node sits a small **control volume**,
and the equations are enforced in integral form: what flows into a control volume flows out.
The fluxes through the control-volume faces are evaluated with finite-element shape functions
inside each hex; each hex contributes 12 such sub-control faces.

| Operator | Meaning |
|---|---|
| face flux `F` | volume flow `u . A` through one sub-control face, stored per face |
| divergence `D_F` | net flux out of each node's control volume (the mass balance) |
| gradient `G` | nodal pressure gradient, `G p = -M^-1 D^T p` (`M`: lumped volume) |
| Laplacian `K` | exchange between neighbouring control volumes through the faces |

Equal-order velocity and pressure need stabilization. A face flux computed from the plain average
of its two nodes, `A . (u_L + u_R) / 2`, cannot see a pressure that alternates from node to node
(a checkerboard). The solver therefore adds a Rhie–Chow term to every face flux:

```
F = A . (u_L + u_R)/2 - h [ (grad p . A)_f - A . (G p_L + G p_R)/2 ],    h = dt_eff / rho
```

The bracket is the difference between the pressure gradient on the face (from the shape
functions) and the average of the two nodal gradients. It is tiny for a smooth pressure and large
for a checkerboard, so it couples every pressure node to its neighbours.

### 3.4 One time step: the projection method

Each step is a BDF2 incremental pressure correction (BDF1 on the first step), with
`dt_eff = dt` for BDF1 and `2 dt / 3` for BDF2:

1. **Predictor**: move the velocity with the advection term and the old pressure gradient.
2. **Viscous step**: apply viscosity implicitly, one linear solve per velocity component.
3. **Projection**: the fluxes `F**` of the predicted velocity violate mass balance. Solve
   `K phi = -(rho / dt_eff) (D_F F** + openings)` for a pressure correction `phi`.
4. **Corrector**: `F = F** - h (grad phi . A)` on every face, `u = u** - h G phi` at the nodes
   with free velocity, and `p = p + phi`.

Step 4 applies to the fluxes the same face gradients from which `K` is built, so
`D_F F + openings = 0` holds to the solver tolerance on any number of GPUs. The nodal velocity
follows the fluxes up to the stabilization term.

### 3.5 The linear solvers

Each implicit step is a sparse linear system `A x = b` with one row per node. The solver never
factors `A`; it iterates. Both systems use preconditioned conjugate gradients (PCG) with Hypre's
algebraic multigrid (BoomerAMG) as preconditioner, and stop when the relative residual is below
`--tol` (default 1e-10). A simple Jacobi preconditioner is not enough for the pressure: on this
channel an earlier Jacobi-CG version needed thousands of iterations per step. With BoomerAMG the
last step of the validation run takes 18 pressure iterations and 10 to 13 per velocity
component, on 1, 2 and 4 GPUs.

Both matrices are constant in time, so the solver assembles them once and builds the multigrid
hierarchies once:

| System | Matrix |
|---|---|
| velocity (each of u, v) | `M / dt_eff + nu K`, rows and columns of fixed velocities removed |
| pressure | `K`, rows of the outlet nodes (p = 0) removed |

Each GPU assembles the matrix of its own elements. Nodes shared by several GPUs are one unknown,
and their rows are summed on the GPU that owns the node (the
[Taylor–Green tutorial](periodic_tgv_tutorial.md) explains this for periodic points too). At setup
the solver compares the assembled `K` with the matrix-free flux correction and stops if they
differ.

### 3.6 The opening fluxes

The discrete divergence `D_F` sums the fluxes of the **interior** sub-control faces. Each interior
face adds `+F` to one node and `-F` to the other, so interior fluxes cancel in pairs. The faces on
the inlet and outlet planes are not interior faces. Without an extra term, the prescribed inflow
never enters the mass balance: the pressure solve sees no fluid entering, builds no pressure
gradient along the channel, and nothing flows.

The solver therefore adds the flux through the openings to every inlet and outlet node, after
the interior divergence:

```
div[i] += A_in[i] . u_prescribed[i] + A_out[i] . u[i]
```

`A_in` and `A_out` are the node's outward face areas on the inlet and outlet planes (the inlet
normal points to -x). The inlet uses the prescribed velocity; the outlet uses the computed one,
so the outlet flux is part of the operator. The skew-symmetric advection gets the matching term:
an opening face with flux `m` adds `-m q / 2` to its node.

The areas are built on the GPU from the element faces on each opening plane. Every face node gets
the area of its quarter of the face (node, edge midpoints and face centre). Each rank adds the
faces of its own elements, and the reverse halo completes the shared nodes, so the areas sum to
`H dz` on any number of ranks. In code: `nsOpeningAreaKernel` and `nsOpeningFlux` in
`backend/distributed/unstructured/fem/mars_navier_stokes.hpp`.

---

## 4. Build and run

The example needs a CUDA build with Hypre (`-DMARS_ENABLE_HYPRE=ON`) and, to read the Exodus
mesh, netCDF:

```bash
cmake -B build -DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_FEM_EXAMPLES=ON -DMARS_ENABLE_HYPRE=ON \
      -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build build --target mars_poiseuille_flow -j
```

Run the validation case on one GPU, from the repository root:

```bash
mpirun -np 1 ./build/examples/distributed/unstructured/mars_poiseuille_flow \
  --mesh=tests/data/poiseuille/poiseuille_hex_14k_elem.e \
  --uinf=1 --nu=0.01 --dt=0.01 --num-steps=1500 --report-every=100 --check \
  --vtu-output=poiseuille --vtu-every=50
```

The same command runs on any number of GPUs (`-np 4`); each rank uses GPU `rank % deviceCount`.
With Slurm, start it with `srun` (for example `srun --ntasks=4 --export=ALL ...`). On one GH200
a step takes about 16 ms, so the 1500 steps take about 24 s. On this 30,000-node mesh more GPUs
are slower, not faster. For scaling runs, generate the channel on every rank instead of reading
a mesh: `--cells=NX,NY` meshes `[0,10] x [0,1] x [0,0.06]` with one cell in z.

| Option | Meaning |
|---|---|
| `--mesh=FILE` or `--cells=NX,NY` | read a mesh, or generate the channel |
| `--y-grading=S` | cluster the generated y spacing at the walls, S in [0, 1) |
| `--uinf --rho --nu --dt --num-steps` | inflow velocity, density, viscosity, time step, number of steps |
| `--bdf1` | first-order time stepping (default BDF2) |
| `--tol --max-iter` | relative AMG-PCG tolerance (default 1e-10) and iteration cap (default 1000) |
| `--report-every=N` | progress line interval (default 1) |
| `--vtu-output=PREFIX --vtu-every=N` | ParaView frames (default every 50 steps) |
| `--check` | validation gate: exit 1 unless every check of section 6 passes |
| `--profile-x --profile-xtol` | move the validation plane (default 90% down the channel, half width 2% of the length) |
| `--comparison-output=PREFIX` | u, v, w, p at full precision, for comparisons with other codes |

`--help` prints all options, including the tolerances of `--check`.

---

## 5. Reading the output

At the start, among the setup lines:

```
Pressure operator: assembled vs matrix-free, max |difference| / max |Kx| = ...
Pressure matrix symmetry: max |a_ij - a_ji| / max |a_ij| = ...
Poiseuille channel: 30000 nodes on 1 ranks, Re = U H / nu = 100
```

The first line is the setup check of section 3.5. It must be at round-off level (1.1e-15 to
1.3e-15 in the validation record); above 1e-10 the run stops. The second line checks that the
matrices are symmetric, which PCG assumes. On box-shaped hexahedra it is at round-off level; on
distorted ones it is not, and the solver prints a warning.

Then one line every `--report-every` steps:

```
Step   1500  t=15.0000  |u|_M=...  continuity=...  amg(u,v,p)=10/13/18
```

- `|u|_M` is the volume-weighted L2 norm of the streamwise velocity, `sqrt(sum M u^2)`. It rises
  from the uniform start and levels off when the flow is developed. For a parabola over the
  whole channel (volume 0.6) it would be `sqrt(1.2 * 0.6) ≈ 0.85`; the entrance region keeps it
  a little lower.
- `continuity` is the largest `|D_F F + openings| / M` over the nodes the projection enforces:
  the mass balance of every control volume, at the solver tolerance.
- `amg(u,v,p)` are the AMG-PCG iterations of the last step. A solve that misses its tolerance
  stops the run.

At the end:

```
[timing] ranks=1 nodes=30000 nodes/rank=30000 steps=1499 ms/step: total=... predictor=... viscous=... pressure=... corrector=... | ...
Profile CSV: poiseuille_profile.csv (plane x=..., 200 nodes)

Poiseuille validation
  profile RMS at x=9.0000 +/- 0.2000: 4.551373e-04 (U_max=1.500000e+00, 1200 nodes)
  flux Q(x)/Q(inlet) at 25/50/75%: 1.0000 / 1.0000 / 1.0000
  -dp/dx from p: 1.2178e-01   from the u profile: 1.1992e-01   exact: 1.2000e-01
```

- `[timing]` gives the time per step of each stage, for the slowest rank, without the first
  step. The validation run measured `total` = 15.9 ms on one GH200.
- The profile RMS compares the computed `u(y)` with the exact parabola on the nodes within
  0.2 of `x = 9` (90% of the channel: past the entrance length, upstream of the outlet).
- The flux ratios are the volume flux through cross-sections at 25, 50 and 75% of the length,
  computed from the solved velocity, over the inlet flux. 1.0000 means the flow carries all the
  mass through the channel.
- `-dp/dx` is measured twice: from the computed pressure between 60% and 90% of the length, and
  from a parabola fit of the velocity. The exact value is `12 rho nu U / H^2 = 0.12`. These two
  are diagnostics, not part of the gate.

With `--vtu-output=PREFIX` the run writes `PREFIX.pvd` with the point fields `u`, `v` and `p`.

---

## 6. The validation gate, the profile plot and an animation

**Validation gate.** With `--check` the example grades itself at the end:

```
VALIDATION PASS: RMS=4.551e-04 < 6.000e-03, flux PASS, steady=2.627e-08 PASS, continuity*H/U=8.050e-15 balance=8.797e-08 PASS, projection PASS
```

It requires:

- **profile**: RMS error below 6e-3;
- **flux**: the three flux ratios within 10% of 1;
- **continuity**: the volume-weighted RMS of `D_F F + openings`, times `H / U`, at most 1e-6,
  and the inlet and outlet fluxes balancing to 1e-6 (the outlet flux uses the nodal velocity, so
  the balance shows the stabilization term, about 1e-7);
- **steadiness**: the velocity change over the last 20 steps, divided by `U`, at most 1e-6;
- **projection**: at steps 1, 2, 3 and the last, the corrected fluxes satisfy
  `D_F F = D_F F** + (dt_eff / rho) K phi` to 1e-7 (printed as `[channel-projection]` lines).

The case is registered with ctest as `marsPoiseuilleValidation` (labels `validation;gpu;long`)
when you configure with `-DMARS_ENABLE_VALIDATION_TESTS=ON`. Run it with `ctest -L validation`
on a node with a GPU. The release tests (`ctest -L release`) run a shorter generated channel
(`--cells=200,40`, 10 steps) on 1 and N ranks.

**Profile plot.** Every run writes `PREFIX_profile.csv` (`poiseuille_profile.csv` without
`--vtu-output`): the computed and exact `u(y)` on the node plane nearest the probe. Plot it with

```bash
python3 scripts/plot_poiseuille_profile.py poiseuille_profile.csv
```

**The exact solution, four ways.** The run checks four related forms of the exact solution:

| Quantity | Checked by |
|---|---|
| `u(y) = G / (2 mu) y (H - y)` | profile RMS and the plot |
| `Q = G H^3 / (12 mu)` | the flux ratios |
| `G = -dp/dx` | the two `-dp/dx` values |
| `U_max = G H^2 / (8 mu) = 1.5 U` | the exact profile the RMS uses |

**Animation.** Render the channel developing, from the uniform inflow to the parabola, with
ParaView (needs `ffmpeg` for the movie):

```bash
pvbatch scripts/render_poiseuille.py --pvd poiseuille     # writes poiseuille_movie.mp4
```

Use `--vtu-every=10` for a smooth movie.

---

## 7. Troubleshooting

| Message or symptom | Cause | Fix |
|---|---|---|
| `NavierStokes: planar flow needs one layer of elements between two z planes` | the example solves planar flow on a one-element-thick mesh | use a mesh with every node on one of two z planes |
| `NavierStokes: the assembled pressure operator is wrong` | assembly and matrix-free operator disagree | a bug: report it with the rank count and the mesh |
| `NavierStokes: the ... solve did not converge` | `--max-iter` too low, or broken input (NaN) | check the last `amg(u,v,p)` line; raise `--max-iter` |
| more GPUs are slower on this mesh | 30,000 nodes are too few to share | scale with `--cells=NX,NY` |

---

## 8. Where to go next

- `examples/distributed/unstructured/mars_poiseuille_flow.cu`: the example, in the order mesh
  and domain, solver, time loop, result.
- `examples/distributed/unstructured/mars_poiseuille_validation.hpp`: the comparison with the
  exact solution and the `--check` gate.
- `backend/distributed/unstructured/fem/mars_navier_stokes.hpp`: the solver, shared with the
  Taylor–Green vortex and the lid-driven cavity.
- [Taylor–Green vortex](periodic_tgv_tutorial.md): the same solver on a periodic box, and how
  the unknowns are shared between GPUs.
- [FEM Assembly](FEM-Assembly.md) and [CVFEM Kernels](CVFEM-Kernels.md): the assembly layer.
