# Periodic Boundary Conditions in the MARS Navier–Stokes Solver

### A beginner's tour of the Taylor–Green vortex, and how one periodic point stays one unknown on any number of GPUs

This tutorial is for a developer who is new to FEM/CVFEM and to domain
decomposition, but who needs to understand how MARS handles **periodic boundary
conditions** in its GPU incompressible Navier–Stokes solver. It explains one idea,
the **reduced periodic space**, shows how the solver's vectors and matrices are
built on it, and tells the story of the bug this design removes on several ranks.

The solver is the same one the Poiseuille channel and the lid-driven cavity use;
only its boundary description differs. You should be able to read this in about
20 minutes.

| File | What lives there |
|------|------------------|
| `backend/distributed/unstructured/fem/mars_periodic_bc.hpp` | Pairing: which node slot is the periodic image of which, and the cross-rank pair tables |
| `backend/distributed/unstructured/fem/mars_dof_space.hpp` | `DofSpace`: the unknowns of a nodal field, `prolong` (P), `restrict` (Pᵀ), dot products |
| `backend/distributed/unstructured/fem/mars_navier_stokes.hpp` | `NavierStokes`: the operators, the matrices and the time step |
| `examples/distributed/unstructured/mars_tgv.cu` | The Taylor–Green example: mesh, periodic pairs, solver, time loop, output |

---

## 1. What is the Taylor–Green vortex, and what does "periodic" mean?

The **Taylor–Green vortex (TGV)** is the standard test for an incompressible
Navier–Stokes solver. You start from a smooth, swirling velocity field in a cube
and let viscosity drain its energy. It is the "hello world" of turbulence codes.

The initial field is set in the example (`tgvInitialConditionKernel`), with one
wavelength per box length, `k = 2π / L`:

```cpp
u[i] = V0 * sin(X) * cos(Y) * cos(Z);
v[i] = -V0 * cos(X) * sin(Y) * cos(Z);
w[i] = T(0);
p[i] = rho * V0 * V0 / T(16) * (cos(T(2) * X) + cos(T(2) * Y)) * (cos(T(2) * Z) + T(2));
```

**"Periodic" means the box wraps around.** The face at `x = lo` is not a wall: it
*is* the same surface as the face at `x = hi`. Flow that leaves on the right comes
back on the left, as in Pac-Man. The same holds in y and z. Glue all three pairs
of opposite faces together and the cube becomes a 3-D torus with no boundary.

So there are no boundary conditions at all. One consequence: the pressure is only
defined up to a constant. The solver fixes it by setting `p = 0` at one unknown.

---

## 2. Two storage slots, one physical point

The cube is filled with hexahedral elements; their corners are the **nodes**. A
node stores one value of every field.

Look at the periodic faces. A node on the `x = lo` face and the node on the
`x = hi` face with the same `(y, z)` are **the same physical point**, because the
box wraps. The mesh still stores them as **two separate nodes**. That is the single
most important fact of this tutorial:

> **Two storage slots, one physical point.**

MARS names the two slots:

- the node on a **min** face is the **master**,
- the node on the matching **max** face is the **slave**.

`buildPeriodicMap` marks every slave with the faces it lies on (`d_periodicMask`).
A node on an edge or a corner of the box is a slave in two or three directions at
once; its **final master** is the node with each of those coordinates moved to `lo`.

Every node is identified by its SFC key, the space-filling-curve code of its
position. So the master's key follows from the slave's: decode the slave key into
its integer coordinates, set each marked axis to the `lo` face, and encode again
(`periodicMasterKey`). No search and no coordinate matching is needed, and the master
does not have to be on the slave's rank.

```
     master (x = lo)                         slave (x = hi)
     key = K with ix = 0                     key = K
         o  ------------ same (y, z) ------------  o
                 one physical point, two slots
```

---

## 3. The reduced periodic space: P and Pᵀ

A periodic point is **one unknown** of every field, however many slots store it.
The unknown lives in one place: the master slot, on the rank that owns the master.
We call that slot the point's **DOF** (degree of freedom). An interior node is its
own DOF.

Two maps connect the slots and the DOFs (`DofSpace`):

- **prolong (P)**: copy the DOF value into every slot of the point, including
  copies on other ranks. After `prolong`, every slot of a periodic point holds the
  same value.
- **restrict (Pᵀ)**: add what every slot of the point received into the DOF slot.

They are exact transposes of each other: copying one value to several places, and
summing several places into one.

```
    DOF value  --prolong (P)-->  master slot, slave slots, ghost copies  (all equal)
    DOF sum    <--restrict (Pᵀ)--  per-slot contributions of an element scatter
```

Now the rule that makes periodicity correct, the same rule MFEM and deal.II use:

> **Every field is kept in range(P), and every operator is Pᵀ A P.**

Here `A` is the plain element operator that knows nothing about periodicity: a loop
over elements that reads its corner nodes and adds a contribution to each of them.
All periodic logic lives in P and Pᵀ. Without a periodic map the same `DofSpace`
handles an ordinary mesh: the only copies are then ghosts on other ranks.

Sums over the unknowns (norms, kinetic energy) count each periodic point once: they
add only over DOF slots, the owned nodes that are not slaves (`DofSpace::isDof`).

---

## 4. The operators and the time step

`NavierStokes` keeps velocity and pressure at the nodes. Its operators come from
the 12 sub-control faces of each hex (area vector `A_f` from node L to node R) and
the sub-control volumes (lumped mass `M`):

| Symbol | Operator |
|--------|----------|
| `F` | the flux through each sub-control face, stored per face |
| `D_F` | face flux divergence: `F` leaves node L and enters node R |
| `G` | nodal gradient `G p = −M⁻¹ Dᵀ p`, where `Dᵀ` adds `(p_L − p_R) / 2 · A_f` to both ends of a face |
| `K` | CVFEM Laplacian from the compact face gradient `(∇p · A)_f` |
| `N(F)` | skew-symmetric advection by the face fluxes |

The face flux is **stabilized** (Rhie–Chow):

```
F = A · (u_L + u_R) / 2  −  h [ (∇p · A)_f − A · (G p_L + G p_R) / 2 ],   h = dtEff / ρ
```

With equal-order nodes, the plain average `A · (u_L + u_R) / 2` cannot see a
pressure that alternates from node to node (a checkerboard). The bracket is the
difference between the compact face gradient and the average of the nodal
gradients: tiny for a smooth pressure, large for a checkerboard. It couples every
pressure node to its neighbours.

One time step is the BDF2 incremental pressure correction (BDF1 on the first step),
with `dtEff = dt` for BDF1 and `2 dt / 3` for BDF2:

1. **Predictor**: `u*` from the BDF history, the extrapolated advection by `F` and the
   old pressure gradient.
2. **Viscous step**: `(M / dtEff + ν K) u** = (M / dtEff) u*`, one solve per component.
3. **Projection**: `F**` is the stabilized flux of `u**` and `p^n`; solve
   `K φ = −(ρ / dtEff) D_F F**`.
4. **Corrector**: `F = F** − h (∇φ · A)`, `u = u** − h G φ`, `p = p^n + φ`.

Apply `D_F` to the corrected flux:

```
D_F F = D_F F** − h D_F (∇φ · A) = D_F F** + h K φ = D_F F** − D_F F** = 0
```

The face fluxes are divergence-free **to the solver tolerance** on every rank
count: the right-hand side, the matrix `K` and the flux correction use the same
face gradients. In the output this is the `div=` column, at 10⁻¹⁵.

**Every vector operation goes through P.** An element loop scatters into the local
slots, `restrict` sums the slots into the DOFs, a node loop updates the DOFs, and
`prolong` copies the result back to every slot. The new velocity is single-valued at
every periodic point by construction.

**Every matrix goes through P too.** `K` and `M / dtEff + ν K` are constant in time,
so they are assembled once and solved with Hypre PCG + BoomerAMG. Each rank
assembles its own elements into a matrix over its local slots, and Hypre forms the
matrix over the DOFs as the product `Pᵀ A_local P`, with `P` stored as a Hypre matrix
(one entry per slot: the global id of its DOF). This is how MFEM assembles. Rank
boundaries and periodic seams are then the same thing: slots that share a DOF.
BoomerAMG keeps the iteration count flat as the mesh and the number of GPUs grow
(about 20 pressure iterations per step from 16³ to 64³).

**Skew advection conserves kinetic energy.** Node L receives `−F q_R / 2` and node R
`+F q_L / 2`, so the kinetic energy production `Σ q_i (N q)_i` cancels face by face.
Only viscosity removes energy.

---

## 5. What changes on several GPUs

With several ranks, each rank owns some elements and some nodes and also holds
**ghost** copies of nodes owned elsewhere. The cstone node halo moves values between
an owner and its ghosts:

- `exchangeNodeHalo`: owner → ghosts (a copy),
- `reverseExchangeNodeHaloAdd`: ghosts → owner (a sum).

The example makes the cornerstone box periodic:

```cpp
amr.initialize(opt.mesh, rank, numRanks, /*periodicAxesMask=*/7, opt.boxLo, opt.boxHi);
```

Node coordinates stay real, and every rank also receives the elements on the other
side of each periodic face. So a rank that owns a slave also holds its master,
usually as a ghost.

A ghost and its owned node share a key; a slave and its master do not. Both cases
follow the same rule: every slot that is not a DOF (a ghost, a slave, or a ghost of a
slave) knows the key of its DOF, its own key or its master's key, and the rank that
holds the DOF is the SFC owner of that key (SFC node ownership: the rank whose range
of the space-filling curve contains it). Every rank can evaluate that without
communication.

At setup `DofSpace` groups its non-DOF slots by that rank and sends each rank the
list of keys it needs, once. The owner answers with the slots that hold those DOFs.
From then on P and Pᵀ are each one exchange, straight between a slot and the rank of
its DOF:

```cpp
// prolong: every copy takes the value of its DOF
dofLocalKernel<<<...>>>(localCopy, localDof, n, fields, /*add=*/false);   // DOF on this rank
exchange(fields, /*reverse=*/false);                                    // DOF on another rank

// restrict: every copy adds into its DOF
exchange(fields, /*reverse=*/true);
dofLocalKernel<<<...>>>(localCopy, localDof, n, fields, /*add=*/true);
```

A ghost copy of a slave on a third rank receives its value directly from the rank of
the master, not through the slave's rank, so there is no order to get right. Several
fields share one exchange: the three velocity components and the pressure travel in
one message per neighbour rank.

The setup stops if a requested key is not a DOF on the rank that owns it: the ranks
disagree on ownership, or a slave has no matching master (box bounds, `faceEps`). The
matrices need nothing extra: the global id of each slot's DOF is itself a field,
prolonged once at setup.

Nothing else in the solver knows about ranks or periodicity. The time step is the
same code on 1 and on N ranks, which is why the results agree to roundoff.

---

## 6. Walking through `mars_tgv.cu`

`runTgv()` follows the steps of the method:

1. **Mesh and domain.** `AmrManager` reads the hex mesh and builds the distributed
   domain with a periodic cornerstone box.
2. **Periodic DOFs.** `buildPeriodicMap` pairs every max-face node with its master.
3. **Solver.** `NavierStokes` with no boundary conditions (`FreeNodes`), no openings
   and the periodic map. It builds the geometry, the matrices over the periodic
   unknowns and their AMG hierarchies, and checks the assembled pressure matrix
   against the matrix-free one. It prints the number of periodic DOFs, which must be
   the same on every rank count (4096 on a 16³ mesh).
4. **Time loop.** `solver->step()` once per step.
5. **Output.** Kinetic energy, continuity and AMG iterations every `--report-every`
   steps, and VTU frames with `--vtu-output` (velocity, pressure, vorticity).

The validation output (`EnergyReport`) is kept apart from these steps. At low
Reynolds number the initial field is an eigenfunction of the Laplacian with
eigenvalue `−3k²`, so the kinetic energy decays like the Stokes solution,
`KE(t) ≈ KE(0) · exp(−6 ν k² t)`. The report prints `KE / KE_Stokes`; on a resolved
mesh at low Re it stays within a few 10⁻³ of 1. At Re = 1600 it does not: the vortex
breaks down, and the ratio only shows that transition.

`mainTgv` destroys the domain, the periodic map and the solver before
`MPI_Finalize`: they hold MPI communicators and Hypre objects.

Mesh adaptation (`--adapt-every`) is not supported by the solver yet. Refinement
leaves hanging nodes on faces between a refined and an unrefined element, and those
need constraints (a hanging node takes the average of its coarse neighbours) that
are not implemented.

---

## 7. The bug story: what went wrong before

Before this design, periodic TGV was correct on 1 rank and blew up on several. It
is worth understanding why, because the mistake is easy to make.

The old solver (`NSStepper` in `mars_ns_solver.hpp`) collapsed the **pressure** to
one unknown per periodic point, including across ranks, but kept the **velocity** as
two independent slots, master and slave. Its corrector subtracted a per-slot
pressure gradient from each slot separately, so the new velocity had
`u[slave] ≠ u[master]`: two values for one physical point.

The next step's advection and divergence read node slots, so the code then copied
`u[slave] := u[master]`. That copy threw away exactly the part of the correction
that made the per-slot field divergence-free. The divergence the projection had
removed came back every step. Three more defects sat on the multi-rank path: a
pressure right-hand side written (not added) from both slots of a pair, so the
result depended on which GPU write landed last; assembled velocity rows that lost
the slave side of every cross-rank pair; and a halo exchange that ran before a copy
it should follow. On every rank count, the BDF2 pressure right-hand side used
`ρ / dt` while the corrector used `2 dt / 3`, so each step removed only two thirds
of the divergence.

None of these can occur with the reduced space. There is no per-slot velocity to
copy, no right-hand side written per slot, no matrix row assembled from a partial
set of elements (every rank assembles only its own elements, and Pᵀ A P sums the
rest), and P and Pᵀ are the only two places that know the order of the halo and the
periodic copy.

---

## 8. Running and checking

Build the example (`make mars_tgv`, needs Hypre) and generate a 16³ unit cube:

```bash
python3 scripts/generate_hex_cube.py --nx 16 --ny 16 --nz 16 --output hex16
```

Low-Reynolds-number check on the unit box, on 1 and 4 ranks:

```bash
mpirun -np 1 ./examples/distributed/unstructured/mars_tgv --mesh=hex16 --box-lo=0 --box-hi=1 \
    --nu=0.05 --dt=1e-4 --num-steps=300 --report-every=50
mpirun -np 4 ./examples/distributed/unstructured/mars_tgv --mesh=hex16 --box-lo=0 --box-hi=1 \
    --nu=0.05 --dt=1e-4 --num-steps=300 --report-every=50
```

What to check:

- `periodic DOFs=4096` on every rank count (one unknown per periodic point, none counted twice).
- `KE` at each report step is the same on every rank count.
- `div` (the continuity of the face fluxes, `max |D_F F| / M`) stays at roundoff, a few 10⁻¹⁵.
- `KE / KE_Stokes` follows the reference, measured on GH200 with 1, 2 and 4 GPUs
  (identical to all printed digits on each):

| step | t | KE | KE / KE_Stokes |
|------|---|----|----------------|
| 0 | 0.000 | 1.2500000000e-01 | 1.00000000 |
| 100 | 0.010 | 1.1120541324e-01 | 1.00150170 |
| 200 | 0.020 | 9.8927086261e-02 | 1.00294417 |
| 300 | 0.030 | 8.7995777160e-02 | 1.00429002 |

The ratio grows slowly above 1 because the discrete Laplacian of a 16-point
wavelength is about 1.3 % weaker than the exact one. The release tests
`marsReleaseTgv_np1` and `_npN` run 100 steps and require `KE / KE_Stokes` in
`[1.0013, 1.0017]`.

---

## 9. Takeaways

1. **A periodic point is one unknown for every field.** Pressure and velocity
   alike. Collapsing only the pressure leaves two velocity values for one point, and
   no copy afterwards can make both divergence-free.
2. **Put all periodic logic in one map, P.** Vectors go through `prolong` and
   `restrict`; matrices are `Pᵀ A P`. Consistency between the right-hand side, the
   operator and the corrector then follows by construction instead of by care.
3. **Let every copy talk to its DOF directly.** A slave's master and a node's owner
   both follow from keys, so P and Pᵀ are one exchange each, with no chain of copies
   whose order could go wrong.
4. **Stabilize equal-order pressure.** Without the face-flux stabilization the
   pressure matrix of a periodic box has checkerboard null modes, and multigrid
   breaks down on them.
5. **Check invariants with one number.** The DOF count must not depend on the rank
   count, and the KE history must agree across rank counts to roundoff.

---

### Quick reference

| Role | Symbol | File |
|------|--------|------|
| Max-face marks | `buildPeriodicMap`, `d_periodicMask` | `mars_periodic_bc.hpp` |
| Master key | `periodicMasterKey` | `mars_periodic_bc.hpp` |
| P, Pᵀ | `DofSpace::prolong`, `DofSpace::restrict` | `mars_dof_space.hpp` |
| Slot to DOF lists (setup) | `DofSpace::build` | `mars_dof_space.hpp` |
| Solver | `NavierStokes` | `mars_navier_stokes.hpp` |
| Stabilized face flux | `nsFaceFluxKernel` | `mars_navier_stokes.hpp` |
| Matrices over the DOFs | `assembleReduced`, `hypreGalerkin` | `mars_navier_stokes.hpp`, `mars_hypre_amg_pcg_solver.hpp` |
| Example | `runTgv` | `mars_tgv.cu` |
