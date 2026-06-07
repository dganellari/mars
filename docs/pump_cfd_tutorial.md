# Simulating Flow Through a Pump with MARS — A From-Scratch Tutorial

This tutorial explains, from the ground up, what the MARS pump simulation does and
how it works. It assumes **no prior CFD knowledge** at the start and builds toward
the real numerical methods, so it doubles as an introduction to computational fluid
dynamics (CFD) for internal pump flow. Later sections add enough depth (the
discrete operators, the projection algebra, stabilization, accuracy orders) to serve
as a working reference once you know the basics.

---

## Part 0 — What are we even computing?

A pump moves liquid: fluid enters through an **inlet**, gets pushed through the
pump's internal passages, and leaves through an **outlet**. We want to predict, at
every point inside the pump, **how fast the fluid moves and in which direction**
(the *velocity field* `u`), and **what the pressure is** (the *pressure field* `p`).

"Computational fluid dynamics" means: instead of solving these by hand (impossible
for a real shape), we chop the pump's interior into millions of tiny pieces, write
the laws of fluid motion on each piece, and let a computer solve the resulting giant
system of equations. MARS does this **on GPUs**, which is what makes millions–
billions of pieces tractable.

The laws of motion are the **incompressible Navier–Stokes equations**:

```
∂u/∂t + (u·∇)u = −(1/ρ)∇p + ν∇²u        (momentum: how velocity changes)
∇·u = 0                                  (incompressibility: mass is conserved)
```

In words:
- **Momentum**: a fluid parcel accelerates due to being pushed by neighboring fluid
  (the `(u·∇)u` *advection* term — fluid carrying itself along), pressure
  differences (`∇p`), and viscous friction (`ν∇²u`, `ν` = how "thick"/sticky the
  fluid is).
- **Incompressibility**: `∇·u = 0` says the fluid neither piles up nor vanishes
  anywhere — **whatever volume flows into a region must flow out**. Liquids are
  (nearly) incompressible, so this holds.

**The whole difficulty of incompressible CFD is that second equation.** There is no
separate equation that tells you the pressure. Pressure is whatever it has to be to
keep the velocity divergence-free. Every method below is, at heart, a trick for
finding that pressure. Hold onto this idea — it explains every design choice.

---

## Part 1 — The mesh: chopping the pump into tetrahedra

### Nodes and elements

We can't store the fluid at every one of infinitely many points, so we pick a finite
set of points — **nodes** (also called vertices) — and store `u` and `p` there. The
nodes are connected into small building-block volumes — **elements** — that tile the
entire interior of the pump.

MARS uses **tetrahedral elements** ("tets"): each is a little pyramid with **4
corner nodes** and 4 triangular faces. A pump's interior — curved passages, an
impeller, a volute — is geometrically messy, and **tetrahedra can fill any shape**
automatically (mesh-generation software does this). That flexibility is why tets are
the standard choice for complicated real-world geometry. (Cube-shaped "hexahedral"
elements are more accurate per element, but meshing a twisty pump with hexes by hand
is extremely hard.)

Key points to internalize:
- **Nodes carry the solution** (`u, v, w, p` — three velocity components and
  pressure — one set per node).
- **Elements carry no values of their own**; an element is just a record of *which 4
  nodes* are its corners. Many tets share each node.
- The counts differ: a mesh has, say, a few hundred thousand nodes and roughly six
  times as many tets.

### The mesh file (Exodus) and side-sets

The mesh comes from a file in **Exodus** format (`.exo` / `.e`) — a standard
finite-element mesh format. MARS reads it directly. The file contains:
- **Node coordinates** — the (x, y, z) of every node (stored as three separate
  arrays X, Y, Z, which is friendly to GPU memory access).
- **Connectivity** — for each tet, the 4 node indices that are its corners.
- **Side-sets** — *named* groups of element faces that mark physical surfaces. This
  is how the mesh author labels "this patch of the boundary is the inlet," "this is
  the outlet," "this is a wall." Side-sets are the bridge from CAD geometry to
  boundary conditions: without them the solver wouldn't know which surface is the
  inlet.

You tell the pump driver the side-set names on the command line
(`--inlet-ss=NAME`, `--outlet-ss=NAME`); **every other side-set is automatically
treated as a solid wall**.

### Why physical scale matters (non-dimensionalization)

A real pump can be tiny (sub-millimeter passages). If we report raw numbers
(velocities, divergence) in physical units, they come out as awkward tiny or huge
values that are hard to interpret. So MARS divides everything by natural scales — a
velocity scale `U` (the inlet speed) and a length scale `L` (the size of the pump,
its bounding-box diagonal) — to get **dimensionless, order-1 numbers** you can read
at a glance. The most important one is the **Reynolds number**:

```
Re = U·L/ν
```

`Re` measures inertia vs. viscosity — *the* knob that sets the flow regime (smooth
and laminar at low `Re`, swirling and turbulent at high `Re`). Rather than guess a
viscosity `ν`, you specify `--Re=` and MARS back-computes `ν = U·L/Re`. (More
dimensionless diagnostics in Part 5.)

---

## Part 2 — Splitting the mesh across GPUs (partitioning)

To run on many GPUs, the mesh is divided so each GPU (MPI rank) owns a chunk. MARS
uses a **space-filling curve (SFC)** — specifically a **Hilbert curve** — from the
cornerstone-octree library.

A space-filling curve is a single 1-D thread that snakes through all of 3-D space
visiting nearby points consecutively. Give every node/element an integer "key" =
its position along the curve, sort by that key, and you get a 1-D ordering of the
whole mesh in which **spatially-close things are close in the list**. Then:
- **Load balance**: cut the sorted list into equal contiguous chunks → every GPU
  gets about the same number of elements, no matter how irregular the geometry.
- **Locality**: each chunk is a compact spatial blob, so neighbors mostly live on
  the same GPU and inter-GPU communication is small.

Each GPU has **owned** elements (its chunk) plus a thin layer of **halo / ghost**
elements — read-only copies of neighbors' elements right across the partition
boundary, needed so a GPU can compute fluxes at the edge of its region. A node sitting
on a partition boundary is **owned by exactly one GPU** and seen as a ghost by the
others; this "owned by exactly one" rule is what keeps global sums correct.

Everything — coordinates, connectivity, the solution — lives in **GPU memory** and
stays there; there is no copying back to the CPU each step. That is the core of
MARS being "GPU-native," and it's what makes large meshes fast.

> *Subtlety (you can skip on first read):* after partitioning, MARS renumbers nodes
> by their Hilbert key, so the original Exodus file-order node numbers are no longer
> valid. To find which local node an inlet-surface point is, MARS recomputes the
> point's Hilbert key from its (x,y,z) and looks it up. Coordinates are the durable
> identity that survives partitioning; file indices are not.

---

## Part 3 — The numerical method

### 3.1 CVFEM: where the equations actually live

MARS uses the **Control-Volume Finite Element Method (CVFEM)** — the same family as
Sandia's Nalu-Wind. The idea:

- Put the unknowns at the **nodes** (vertices).
- Around each node, build a little **control volume** (its "median dual" cell) by
  connecting edge midpoints and element/face centroids. Every node owns a small
  polyhedral cell; these cells tile the whole domain with no gaps.
- The boundaries between neighboring nodes' cells are **sub-control-surfaces (SCS)**.
- The discrete equations are **flux balances**: for each node, sum what flows in and
  out across its SCS facets. Because each facet is shared by exactly two cells with
  opposite sign, **mass and momentum are conserved locally by construction** (what
  leaves one cell enters its neighbor exactly).

This is a hybrid: it has the **local conservation** of finite volume *and* the
**geometric flexibility** of finite elements on unstructured tets — which is exactly
what you want for a messy pump geometry.

Each SCS facet has a precomputed **area vector** `A_f` (direction = facet normal,
length = facet area). All the operators below are "scatter the flux `(·)·A_f` to the
two endpoint nodes" — simple, conservative, GPU-friendly.

### 3.2 The three discrete operators

- **Divergence `D`** — for each node, `(∇·u)_i ≈ (1/V_i) Σ_f (u_f · A_f)`, with the
  face value `u_f = ½(u_L + u_R)`. This is the discrete "is mass conserved here?"
  measure.
- **Gradient `Dᵀ`** — the *exact transpose* of `D`: scatter `½(p_L − p_R)·A_f` to each
  facet's two endpoints. Used to turn a pressure field into a force on the velocity.
- **Lumped mass `M`** — a diagonal matrix whose entry `M_i` is just the volume of
  node `i`'s control cell.

Why "scatter to two endpoints with opposite sign" matters: it makes the discrete
operators **telescope**, so they satisfy the same identities as the continuous ones.
MARS checks these as unit invariants: `D(constant) = 0` (a uniform flow has no
divergence), `gradient(linear) = constant`, and the area-vector closure
`Σ_f A_f·(x_R − x_L) = 3V` per tetrahedron (a discrete divergence theorem). The lumped
mass is also rank-consistent — `Σ M` over owned nodes is bit-identical on 1 and N GPUs
and equals the domain volume. These invariants are how you trust the geometry before
trusting the flow.

These three operators combine into the pressure operator (Part 3.4), which is the
heart of the method.

### 3.3 Marching in time: BDF2 + projection

We step the solution forward in small time increments `dt`. Each step has to (a)
advance momentum and (b) find the pressure that keeps `∇·u = 0`. MARS uses a
**fractional-step (Chorin) projection** — solve momentum *ignoring* the constraint,
then *project* the result back onto the divergence-free state. The per-step recipe:

1. **Predictor** — compute a provisional velocity `u*` from the momentum equation
   using the advection term and the *old* pressure. `u*` is **not**
   divergence-free (it doesn't yet know the new pressure).
2. **Implicit diffusion** — the viscous term is "stiff" (it would force a
   microscopic `dt` if done explicitly), so solve it implicitly: three linear
   systems (one each for `u, v, w`) of the form `(M/dt + νK) u** = …`, solved with
   Conjugate Gradient (CG). This is the `cg_uvw` count in the output.
3. **Pressure solve** — find the pressure increment `φ` that will remove the
   divergence of `u**` (Part 3.4). Solved with CG; this is `cg_p`.
4. **Corrector** — subtract the pressure gradient: `u^{n+1} = u** − (dt/ρ)·∇φ`, and
   update `p`. Now `∇·u^{n+1} = 0`.

**Why projection works (the one big idea):** any vector field can be split into a
divergence-free part plus a pure gradient (the *Helmholtz decomposition*). `u**` has
some leftover gradient (divergent) part. The pressure solve finds exactly the
potential `φ` whose gradient *is* that divergent part; subtracting `∇φ` in the
corrector removes it, leaving the divergence-free field. **Pressure is the
projector.**

**Time accuracy — BDF2:** the time derivative uses a 2nd-order backward formula
(BDF2) over the last three time levels, and the advection is extrapolated to 2nd
order (EXT2). This is the same well-tested scheme NekRS and deal.II use. (A
first-order option, `--bdf1`, exists but is less accurate; plain first-order
time-stepping was found unstable for through-flow, which is *why* BDF2 is the
default.)

### 3.4 The pressure operator `D M⁻¹ Dᵀ` — the key choice

This is the single most important detail of the pump solver.

For the projection to actually zero the divergence, the pressure operator and the
corrector's gradient have to be **consistent** — the gradient you correct with must
be the exact transpose of the divergence you're trying to kill. The "obvious"
pressure operator (a standard finite-element Laplacian `K`) is **not** consistent
with the median-dual gradient on tetrahedra. On tets the two are nearly at right
angles to each other — the angle between the finite-element gradient and the
median-dual gradient was *measured* at `cos ≈ −0.21` on tets (vs `−0.73` on hexes) —
so correcting with one while solving the other **leaves the divergence almost
untouched, and on tets it amplifies it each step — the run blows up**. This is a
*consistency* failure, not a coding bug: the geometry can be perfectly correct and
the projection still fails because the two operators don't match.

MARS instead uses the **literal composition `A = D M⁻¹ Dᵀ`** as the pressure
operator. Then the algebra closes exactly:

```
∇·u^{n+1} = D u** − (dt/ρ) · D M⁻¹ Dᵀ φ
          = D u** − (dt/ρ) · (ρ/dt) D u**      (because φ solves A φ = (ρ/dt) D u**)
          = 0
```

Because the *same* `D`, `M⁻¹`, `Dᵀ` kernels build both the pressure operator and the
corrector gradient, the divergence cancels **algebraically**. This operator is never
stored as a matrix — applying it inside CG is just three GPU passes (`Dᵀ` then `M⁻¹`
then `D`), which is essential for huge meshes (a stored matrix would exhaust GPU
memory). The CG uses a **Jacobi preconditioner** built matrix-free from the same area
vectors (each facet contributes `+0.25·|A_f|²·(1/M_L + 1/M_R)` to the diagonal,
strictly positive so the operator stays well-posed). On stretched boundary-layer
tetrahedra the diagonal spans several orders of magnitude, so Jacobi is a weak
preconditioner and the pressure CG can take a few hundred iterations — that is normal
for this operator, not a sign of trouble (a stronger multigrid preconditioner is a
possible future improvement).

> **Aside — the checkerboard mode (for the curious).** Putting velocity *and*
> pressure at the same nodes (equal-order, "collocated") technically violates the
> inf-sup / LBB stability condition, which admits a spurious oscillating "checkerboard"
> pressure pattern. The `D M⁻¹ Dᵀ` projection suppresses it well in practice, but on
> very large meshes a dedicated *pressure stabilization* may be needed. The textbook
> cure (Rhie–Chow momentum interpolation) turns out to be geometrically unsafe on
> skewed tetrahedra; the more robust route is a PSPG / Bochev–Dohrmann pressure
> stabilization. This is an active area, not a solved one — worth knowing if you push
> to very fine meshes.

### 3.5 Advection schemes (how fluid carries itself)

The `(u·∇)u` term — fluid transporting its own momentum — is the trickiest to
discretize. MARS offers three, selectable on the command line:

- **Skew-symmetric (`--skew`, the default).** Discretizes advection in a form that
  **conserves discrete kinetic energy exactly** — advection can only *move* energy
  around, never create it. That makes it stable with no artificial smearing. This is
  the right partner for BDF2 and the recommended scheme for the pump (a through-flow
  has no walls to dissipate spurious energy, so an energy-conserving scheme is
  essential — non-conserving schemes blow up).
- **First-order upwind (`--upwind`).** Takes the value from the upstream node. Robust
  but adds heavy *numerical diffusion* that smears smooth flow — only 1st-order
  accurate. Fine near shocks, poor for smooth internal flow, and here it actually
  blew up (no wall dissipation to absorb its error).
- **Barth–Jespersen (`--bj`).** A 2nd-order *limited* reconstruction:
  `q_face = q_up + φ·(∇q·r)` where the limiter `φ ∈ [0,1]` keeps the reconstruction
  from overshooting the local neighbor values (no spurious oscillations) — a
  monotone, TVD-like 2nd-order scheme widely used in production finite-volume CFD.

### 3.6 Orders of accuracy (summary)

| Aspect | Scheme | Order |
|---|---|---|
| Space | CVFEM median-dual (linear tets) | ~2nd order |
| Time | BDF2 + EXT2 | 2nd order |
| Advection — skew | energy-conserving | 2nd order |
| Advection — upwind | upwind | 1st order |
| Advection — Barth–Jespersen | limited reconstruction | 2nd order |

Default pump config (`--skew` + BDF2 + CVFEM) is **2nd-order in space and time**.

---

## Part 4 — Boundary conditions (telling the solver about the openings and walls)

An internal pump flow is sealed except for the inlet and outlet; everything else is
solid wall. MARS maps the named side-sets onto:

- **Inlet — prescribed velocity (Dirichlet).** Fluid enters at speed `--inlet-velocity`
  along the **inward normal** of the inlet face. MARS computes that normal direction
  from the inlet surface triangles (area-weighted), summed across all GPUs so every
  rank agrees on the direction.
- **Outlet — two options.**
  - **Do-nothing / pressure outlet (`--outlet=do-nothing`, default):** the velocity
    is left free and the outlet face is held at `p = 0`. Fluid leaves naturally at
    whatever rate the solution wants. This is the **stable** choice.
  - **Mass-conserving velocity outlet (`--outlet=mass-conserving`):** force the
    outflow speed to exactly balance the inflow (`U_out·A_out = U_in·A_in`). Tidier
    on paper but over-constrains the problem and can go unstable; use the do-nothing
    outlet unless you have a reason not to.
- **Walls — no-slip (`u = 0`).** Every non-inlet/outlet side-set. Fluid sticks to
  solid surfaces.

**Why you must "pin" the pressure.** With velocity prescribed almost everywhere, the
pressure equation only ever sees pressure *differences* — so the absolute pressure
is undefined (you can add any constant). Mathematically the pressure system is
**singular** (it has a constant "null mode"). You fix this by anchoring one
reference: the `p = 0` outlet face (do-nothing) or a single pinned pressure node
(mass-conserving). Skip the pin and the pressure solve fails. (This is a classic
gotcha for newcomers — "my pressure solver won't converge" on an all-velocity-BC
problem almost always means a missing pressure reference.)

**Mass balance.** Incompressible means in = out. The do-nothing outlet lets the
right amount leave on its own (the pressure reference makes it consistent); the
mass-conserving outlet enforces it exactly. Either way the pressure projection gets a
divergence it can actually remove.

---

## Part 5 — Running the pump and reading the output

### A typical command

```bash
# 4 GPUs via srun (Alps); use mpirun -np 4 elsewhere
srun -N1 -n4 ./mars_pump \
  --mesh=PUMP.exo \
  --inlet-ss=<INLET_SS> --outlet-ss=<OUTLET_SS> \
  --inlet-velocity=0.5 \
  --Re=100 \
  --dt=1e-3 \
  --num-steps=3000 \
  --skew \
  --outlet=do-nothing \
  --vtu-output=pump_out --vtu-every=20
```

Common flags: `--mesh=` (Exodus file), `--inlet-ss=` / `--outlet-ss=` (required side-set
names), `--inlet-velocity=` (inflow speed `U`), `--Re=` (sets `ν`), `--dt=` (time
step), `--num-steps=`, the advection scheme (`--skew` / `--bj` / `--upwind`),
`--outlet=` (`do-nothing` or `mass-conserving`), and `--vtu-output=PREFIX`
`--vtu-every=N` for visualization output. (`--vtu-output` is a file *prefix*, not a
folder — it writes `PREFIX.pvd` and `PREFIX_step####.pvtu` in the current directory.)

### The setup banner — sanity checks before stepping

On startup, rank 0 prints invariants that must hold; they catch a broken setup
before you waste a long run:
- `[owned-node check]` — total owned DOFs across ranks must equal the global node
  count (else boundary nodes are double-counted).
- `[mass-sum check]` — summed control-volume mass = the domain volume, identical on
  1 vs N ranks.
- `[bc-count check]` — number of boundary DOFs, identical on 1 vs N ranks (else BC
  nodes were dropped or doubled).
- `[ss-resolve] found G … resolved R` — how many inlet/outlet nodes were matched.
- It **aborts early** with a clear message if the inlet or outlet matched zero nodes
  (otherwise you'd get a singular pressure system and a cryptic crash).

### The per-step line

```
Step 3000  ft=3.48  u_rms=3.02e-1  u_max=6.68  uMax/U=13.4  d(u_rms)=9.8e-6  div*L/U=2.66  cg_p=554  pres_res=9.8e-9  cg_uvw=63
```

- `ft` — **flow-throughs** elapsed (`= t / (L/U)`): how many times fluid has crossed
  the pump. Developed flow needs several; `ft=3.5` is well-developed.
- `u_rms`, `u_max`, `uMax/U` — volume-RMS speed, peak interior speed, and peak over
  inlet speed. `uMax/U = 13.4` means the jet accelerates to 13× the inlet speed
  through a constriction — a strong, real jet.
- `div*L/U` — **dimensionless divergence** = how badly mass conservation is violated;
  should be small. (Note: this *nodal* divergence reads a few units even when the
  projection is doing its job; it's a conservative upper-bound diagnostic.)
- `d(u_rms)` — steady-state residual (how much the flow changed this step). Shrinking
  toward zero means converging to steady state; `9.8e-6` is essentially converged.
- `cg_p` — pressure-CG iterations. A real number that converges (here `pres_res ≈
  1e-9`) is healthy. **`cg_p = -2` means the pressure solve FAILED** — the run is
  blowing up.
- `cg_uvw` — velocity-diffusion CG iterations.

**Healthy run:** `cg_p` converging (not `-2`), `div*L/U` bounded and dropping,
`uMax/U` building as the jet develops, `d(u_rms)` shrinking. **Blowing up:** `cg_p =
-2`, kinetic energy / `u_max` exploding to huge numbers (often from `--upwind` or
`--outlet=mass-conserving`, or too large a `dt`).

### What a converged pump run looks like

A healthy 4-GPU developed run (Re=100, skew, do-nothing outlet, 3000 steps) reaches
`ft≈3.5`, `uMax/U≈13` (strong jet), `d(u_rms)≈1e-5` (converged), `div*L/U≈2.7`, with
`cg_p` converging to `~1e-9` every step. That's a steady, physical inlet→outlet jet.

---

## Part 6 — Visualizing the flow (ParaView)

With `--vtu-output=`, the pump writes per-step `.pvtu` files plus a `.pvd` time
collection that ParaView opens as an animation. Each frame carries `u, v, w, p` and a
`velocity` vector at every node.

There's a render script, [scripts/render_pump_flow.py](../scripts/render_pump_flow.py)
(headless ParaView / `pvbatch`), with several ways to show the flow:
- **Threshold isovolume** (default, fast/robust): shows the fast-moving jet region,
  which grows as the flow develops.
- **Streamlines** (`--streamlines`): traces instantaneous flow paths from the inlet.
- **Pathlines / particle tracks** (`--pathlines`): releases particles at the inlet
  and lets the flow carry them — looks like real fluid moving through.

Practical notes learned the hard way:
- Run the heavier streamline/pathline modes **on a GPU compute node** (under `srun`),
  not the login node — the integrator needs the GPU and memory.
- **Render the developed flow** (a long run, several flow-throughs). On an
  under-developed run the jet has only penetrated partway, so streamlines look short
  and "stop" — that's physics (the flow genuinely hasn't filled the domain yet), not
  a bug.
- If the cluster has no `ffmpeg`, render PNG frames and encode the video on your
  laptop:
  `ffmpeg -framerate 30 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p out.mp4`.

---

## One-paragraph mental model

> MARS simulates pump flow by chopping the interior into **tetrahedra**, putting the
> velocity and pressure at the **nodes**, and writing **conservation balances** on a
> control volume around each node (**CVFEM**). It marches forward in time with a
> **projection method**: guess the velocity from momentum, then solve a **pressure**
> equation — built as the consistent **`D M⁻¹ Dᵀ`** operator so the correction
> *exactly* removes any mass-conservation error — and subtract the pressure gradient.
> Time stepping is **2nd-order (BDF2)**; advection is **energy-conserving skew-
> symmetric** by default. The mesh is split across GPUs by a **Hilbert space-filling
> curve**, everything lives on the GPU, and a healthy run develops a steady,
> physical inlet→outlet jet over a few flow-throughs.

---

## Where to go next

- Mesh reading & partitioning internals: [Mesh-Reading-and-Partitioning.md](Mesh-Reading-and-Partitioning.md),
  [SFC-Mapping.md](SFC-Mapping.md), [Multi-Rank-Support.md](Multi-Rank-Support.md).
- The related periodic Taylor–Green tutorial (same solver family, a canonical CFD
  test case): [periodic_tgv_tutorial.md](periodic_tgv_tutorial.md).
