# Periodic Boundary Conditions in the MARS Navier–Stokes Solver

### A beginner's tour of the Taylor–Green vortex, and how one periodic point stays one unknown on any number of GPUs

This tutorial is for a developer who is new to FEM/CVFEM and to domain
decomposition, but who needs to understand how MARS handles **periodic boundary
conditions** in its GPU incompressible Navier–Stokes solver. It explains one idea,
the **reduced periodic space**, shows how every operator is built on it, and tells
the story of the bug that this design removes on several ranks.

You should be able to read it in about 20 minutes. Code is quoted from the repo
with `file:line` references so you can jump straight to the source.

| File | What lives there |
|------|------------------|
| `backend/distributed/unstructured/fem/mars_periodic_bc.hpp` | Pairing: which node slot is the periodic image of which, and the cross-rank pair tables |
| `backend/distributed/unstructured/fem/mars_periodic_space.hpp` | The reduced periodic space: `prolong` (P), `restrict` (Pᵀ), dot products |
| `backend/distributed/unstructured/fem/mars_periodic_ns.hpp` | `PeriodicNavierStokes`: the operators and the time step |
| `examples/distributed/unstructured/mars_tgv.cu` | The Taylor–Green driver: mesh, periodic DOFs, solver, time loop, output |
| `tests/periodic/check_periodic_space.py` | Host model of the same algebra on 1–4 simulated ranks (no GPU needed) |

---

## 1. What is the Taylor–Green vortex, and what does "periodic" mean?

The **Taylor–Green vortex (TGV)** is the standard test for an incompressible
Navier–Stokes solver. You start from a smooth, swirling velocity field in a cube
and let viscosity drain its energy. It is the "hello world" of turbulence codes.

The initial field is set in the driver (`mars_tgv.cu:146`), with one wavelength per
box length, `k = 2π / L`:

```cpp
u[i] = V0 * sin(X) * cos(Y) * cos(Z);
v[i] = -V0 * cos(X) * sin(Y) * cos(Z);
w[i] = T(0);
p[i] = rho * V0 * V0 / T(16) * (cos(T(2) * X) + cos(T(2) * Y)) * (cos(T(2) * Z) + T(2));
```

**"Periodic" means the box wraps around.** The face at `x = lo` is not a wall: it
*is* the same surface as the face at `x = hi`. Flow that leaves on the right comes
back on the left, exactly as in Pac-Man. The same holds in y and z. Glue all three
pairs of opposite faces together and the cube becomes a 3-D torus with no boundary.

So there are no boundary conditions at all. One consequence: the pressure is only
defined up to a constant. MARS fixes the constant by giving the pressure a zero mean.

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

`buildPeriodicMap` (`mars_periodic_bc.hpp:925`) finds the pairs by coordinates and
stores them in `d_periodicPartner` (`mars_periodic_bc.hpp:109`): for every slave
slot, the slot of its master; `-1` for every other slot. A node on an edge or a
corner of the box is a slave in two or three directions at once;
`flattenPartnerChainKernel` (`mars_periodic_bc.hpp:274`) follows the chain so that
the partner is always the **final master**, the node with every coordinate at `lo`.

```
     master (x = lo)                         slave (x = hi)
     partner = -1                            partner = master slot
         o  ------------ same (y, z) ------------  o
                 one physical point, two slots
```

---

## 3. The reduced periodic space: P and Pᵀ

A periodic point is **one unknown** of every field, however many slots store it.
The unknown lives in one place: the master slot, on the rank that owns the master.
We call that slot the point's **DOF** (degree of freedom). An interior node is its
own DOF.

Two maps connect the slots and the DOFs (`mars_periodic_space.hpp:168`):

- **prolong (P)**: copy the DOF value into every slot of the point, including
  copies on other ranks. After `prolong`, every slot of a periodic point holds the
  same value. We say the field is *in range(P)*.
- **restrict (Pᵀ)**: add what every slot of the point received into the DOF slot.

They are exact transposes of each other: copying one value to several places, and
summing several places into one.

```
    DOF value  --prolong (P)-->  master slot, slave slots, ghost copies  (all equal)
    DOF sum    <--restrict (Pᵀ)--  per-slot contributions of an element scatter
```

Now the rule that makes periodicity correct, the same rule MFEM and deal.II use:

> **Every field is kept in range(P), and every operator is `restrict(A(prolong(x)))` = Pᵀ A P.**

Here `A` is the plain element operator that knows nothing about periodicity: a loop
over elements that reads its corner nodes and adds a contribution to each of them.
All periodic logic lives in P and Pᵀ.

Inner products count each periodic point once: `dot` (`mars_periodic_space.hpp:219`)
sums only over DOF slots, which are the owned nodes that are not slaves
(`periodicDofMaskKernel`, `mars_periodic_space.hpp:50`).

---

## 4. The operators, and why the projection closes exactly

`PeriodicNavierStokes` (`mars_periodic_ns.hpp:398`) builds five operators, all as Pᵀ A P:

| Symbol | Operator | Element kernel |
|--------|----------|----------------|
| `M` | lumped mass: each corner gets 1/8 of the element volume | `pnsLumpedMassKernel` (`mars_periodic_ns.hpp:66`) |
| `D` | divergence: the flux `A_f · (u_L + u_R) / 2` through each sub-control face leaves node L and enters node R | `pnsDivergenceKernel` (`:92`) |
| `G` | gradient, **the transpose of D**: `(p_L − p_R) / 2 · A_f` added to both ends | `pnsDivergenceTransposeKernel` (`:115`) |
| `K` | CVFEM viscous stiffness, stored per element | `pnsStiffnessKernel` (`:181`) |
| `N(u)` | explicit advection, skew-symmetric (default) or upwind | `pnsAdvectionKernel` (`:144`) |

Because `D` and `G` are built from the same faces and the same P, `G = Dᵀ` holds
for the reduced operators too. `M⁻¹ G p` approximates `−∇p`.

One time step (`step()`, `mars_periodic_ns.hpp:456`) is the BDF2 incremental pressure
correction (BDF1 on the first step), with `dtEff = dt` for BDF1 and `2 dt / 3` for BDF2:

1. **Predictor** (`:464`): `u*` from the BDF extrapolation, the extrapolated advection
   and the old pressure gradient.
2. **Viscous step** (`:477`): `(M / dtEff + ν K) u** = (M / dtEff) u*`, one CG per component.
3. **Projection** (`:492`): solve `A φ = −(ρ / dtEff) D u**` with `A = D M⁻¹ G`.
4. **Corrector** (`:509`): `u^{n+1} = u** + (dtEff / ρ) M⁻¹ G φ`, `p^{n+1} = p^n + φ`.

Apply `D` to the corrector and use the projection equation:

```
D u^{n+1} = D u** + (dtEff/ρ) D M⁻¹ G φ = D u** + (dtEff/ρ) A φ = D u** − D u** = 0
```

This holds **exactly** (to the CG tolerance) because the right-hand side, the
operator `A` and the corrector use the same `D`, the same `G = Dᵀ` and the same `M`,
all built with the same P. And because the corrector computes one value per DOF and
then prolongs it, the new velocity is single-valued at every periodic point without
any extra copy.

The pressure operator, applied inside the CG (`applyProjection`, `mars_periodic_ns.hpp:634`):

```cpp
space_.prolong(x);          // P
gradientOf(x);              // restrict(G x), times M^-1
space_.prolong(gx_);        // P again: the gradient is a field too
space_.prolong(gy_);
space_.prolong(gz_);
divergenceOf(gx_, gy_, gz_, out);   // restrict(D g)
```

Two more consequences of the rule:

- **The operators are symmetric.** `D M⁻¹ G` and `M / dtEff + ν K` stay symmetric
  in the reduced space, so plain Jacobi-preconditioned CG (`pcg`, `:648`) works on
  every rank count.
- **Skew advection conserves kinetic energy.** The skew-symmetric flux produces
  kinetic energy only in proportion to `D u`. Once the projection makes `D u = 0`
  per DOF, advection neither creates nor destroys energy, and only viscosity
  dissipates it. This is why the default is `--skew=1`.

---

## 5. What changes on several GPUs

With several ranks, each rank owns some elements and some nodes and also holds
**ghost** copies of nodes owned elsewhere. The cstone node halo moves values between
an owner and its ghosts:

- `exchangeNodeHalo`: owner → ghosts (a copy),
- `reverseExchangeNodeHaloAdd`: ghosts → owner (a sum).

The driver makes the cornerstone box periodic (`mars_tgv.cu:346`):

```cpp
amr.initialize(opt.mesh, rank, numRanks, /*periodicAxesMask=*/7, opt.boxLo, opt.boxHi);
```

Node coordinates stay real, and every rank also receives the elements on the other
side of each periodic face. So a rank that owns a slave also holds its master,
usually as a ghost.

The halo alone cannot connect a slave to its master: they have **different SFC
keys**, so for cstone they are two unrelated nodes. When the slave is owned on rank A
and the master on rank B, `PeriodicMap` keeps a small pair table
(`buildCrossRankPeriodicMap`, `mars_periodic_bc.hpp:588`) and two exchanges use it:
`crossRankPeriodicBroadcast` (master value to the slave's rank, `:1369`) and
`crossRankPeriodicPairSum` (slave contribution to the master's rank, `:1027`).

With these, P and Pᵀ are three calls each (`mars_periodic_space.hpp:197` and `:209`):

```cpp
void prolong(Vector& v) const
{
    periodicBroadcastSameRankKernel<RealType><<<...>>>(partner, ownership, n_, v.data());  // slave <- master, same rank
    crossRankPeriodicBroadcast<KeyType, RealType>(map_, v);                              // slave <- master, other rank
    domain_.exchangeNodeHalo(v);                                                         // ghosts <- owners
}

void restrict(Vector& acc) const
{
    domain_.reverseExchangeNodeHaloAdd(acc);                                             // owners += ghosts
    periodicPairSumKernel<RealType><<<...>>>(partner, ownership, n_, acc.data());        // master += slave, same rank
    crossRankPeriodicPairSum<KeyType, RealType>(map_, acc, /*broadcastBack=*/false);     // master += slave, other rank
}
```

**The order matters, and it is forced by what each step reads.** In `prolong`,
the slave copies are updated first and the halo runs last, so that ghost copies of a
slave on a third rank also receive the new value. In `restrict`, the halo runs first,
so that every owned slave holds its complete sum before it is added into its master.

The constructor checks the two facts these maps rely on and stops if either fails
(`checkPeriodicPairing`, `mars_periodic_space.hpp:135`): every owned slave reaches a
final master in its local halo, and every slave whose master is owned elsewhere is in
the cross-rank table.

Nothing else in the solver knows about ranks. The time step is the same code on 1
and on N ranks, which is why the results agree to roundoff.

---

## 6. Walking through `mars_tgv.cu`

`main()` (`mars_tgv.cu:322`) follows the five steps of the method:

1. **Mesh and domain** (`:337`). `AmrManager` reads the hex mesh and builds the
   distributed domain with a periodic cornerstone box.
2. **Periodic DOFs** (`:349`). `buildPeriodicMap` pairs every max-face node with its
   master.
3. **Operators** (`:355`). `PeriodicNavierStokes` builds the reduced periodic space,
   the element geometry (sub-control face area vectors, element viscous matrices),
   the lumped mass and the Jacobi diagonals. It prints the number of periodic DOFs,
   which must be the same on every rank count (4096 on a 16³ mesh).
4. **Time loop** (`:375`). `solver->step()` once per step.
5. **Output** (`:393`). Kinetic energy, divergence and CG iterations every
   `--report-every` steps, and VTU frames with `--vtu-output`.

The validation output (`EnergyReport`, `mars_tgv.cu:181`) is kept apart from these
steps. At low Reynolds number the initial field is an eigenfunction of the
Laplacian with eigenvalue `−3k²`, so the kinetic energy decays like the Stokes
solution, `KE(t) ≈ KE(0) · exp(−6 ν k² t)`. The report prints `KE / KE_Stokes`; on a
resolved mesh at low Re it stays within a few 10⁻³ of 1. At Re = 1600 it does not:
the vortex breaks down, and the ratio only shows that transition.

Mesh adaptation (`--adapt-every`) is optional and separate (`adaptMesh`,
`mars_tgv.cu:281`): it refines where the velocity is largest, rebuilds the periodic
pairing and the solver on the new mesh, and restarts the BDF history.

---

## 7. The bug story: what went wrong before

Before this design, periodic TGV was correct on 1 rank and blew up on several. It
is worth understanding why, because the mistake is easy to make.

### The pressure was collapsed, the velocity was not

The old solver (`NSStepper` in `mars_ns_solver.hpp`) did collapse the **pressure**
to one unknown per periodic point, including across ranks (the "owner-migration"
DOF collapse, `mars_ns_solver.hpp:2749-2771`). It kept the **velocity** as two
independent slots, master and slave. Its multi-rank pressure operator divided the
gradient by the mass slot by slot, without first summing it into the DOF and
copying it back, and the corrector subtracted that per-slot gradient from each
slot separately. The result was a
velocity with `u[slave] ≠ u[master]`: two values for one physical point.

The next step's advection and divergence read node slots, so the code then copied
`u[slave] := u[master]` (`mars_ns_solver.hpp:8932-8951` and `:9168-9184`). That copy
throws away exactly the part of the correction that made the per-slot field
divergence-free. The divergence the projection had removed came back every step.

### Four concrete defects

The host model (`tests/periodic/check_periodic_space.py`, part 1) transcribes the old
multi-rank kernels and measures each defect on an 8³ mesh with 2 and 4 simulated ranks:

| Defect | Where | Host-model measurement |
|--------|-------|------------------------|
| Per-slot correction, then `u[S] := u[M]` | reduced-path corrector | `|D u^{n+1}| / |D u**|` = **0.22** even with an exact solve (should be 0) |
| The pressure RHS is written, not added, from both slots of a same-rank pair; the slave slot is 0 after the fold, and which write survives is undefined on a GPU | `buildPressureRhsKernel`, `mars_ns_solver.hpp:1989` | if the slave wins: 0.54–0.64 |
| Assembled velocity rows lose the slave side of every cross-rank pair: the slave's row maps to a ghost DOF and is dropped, and on the master's rank the slave is a ghost | DOF collapse + `mars_cvfem_hex_kernel_tensor.hpp:255` | 47–53 % row error at those points, 10⁻¹⁶ elsewhere |
| The fully folded operator (the single-rank path) looked non-symmetric on several ranks only because the halo ran **before** the master → slave copy of the gradient | `mars_ns_solver.hpp:5409-5412` then `:5484-5520` | symmetry error 1.2–1.3·10⁻², and 10⁻¹⁶ with the order swapped |

On Alps the old projection probe measured `|D u^{n+1}| / |D u**| ≈ 1.0` on 4 ranks,
consistent with the first two rows together.

A fifth, smaller issue affected every rank count: with BDF2 the old pressure
right-hand side used `ρ / dt` while the corrector used `2 dt / 3`
(`mars_ns_solver.hpp:7725` vs `:8455-8456`), so each step removed only two thirds of
the divergence. The kinetic energy was barely affected (the same to 7 digits on the
reference case), but `div` was 10⁻⁶ instead of 10⁻⁹.

### What the design changes

None of the five is fixed by a local patch in `PeriodicNavierStokes`; they cannot
occur. There is no per-slot velocity to copy, no right-hand side assembled per slot,
no assembled matrix whose rows can miss elements, and P and Pᵀ are the only two
places that know the order of the halo and the periodic copy.

---

## 8. Running and checking

Build the driver (`make mars_tgv`) and generate a 16³ unit cube:

```bash
python3 scripts/generate_hex_cube.py --nx 16 --ny 16 --nz 16 --output hex16
```

Low-Reynolds-number check on the unit box, on 1, 2 and 4 ranks:

```bash
mpirun -np 1 ./examples/distributed/unstructured/mars_tgv --mesh=hex16 --box-lo=0 --box-hi=1 \
    --nu=0.05 --dt=1e-4 --num-steps=300 --report-every=50
mpirun -np 4 ./examples/distributed/unstructured/mars_tgv --mesh=hex16 --box-lo=0 --box-hi=1 \
    --nu=0.05 --dt=1e-4 --num-steps=300 --report-every=50
```

What to check:

- `periodic DOFs=4096` on every rank count (one unknown per periodic point, none counted twice).
- `KE` at each report step is the same on 1, 2 and 4 ranks to about 10⁻¹⁰ relative or better.
- `div` stays near 10⁻⁹ (CG tolerance level), on every rank count.
- `KE / KE_Stokes` follows the host model reference, printed by
  `python3 tests/periodic/check_periodic_space.py --reference --n 16 --steps 300`:

| step | t | KE | KE / KE_Stokes |
|------|---|----|----------------|
| 0 | 0.000 | 1.2500000000e-01 | 1.00000000 |
| 100 | 0.010 | 1.1120567302e-01 | 1.00150404 |
| 200 | 0.020 | 9.8927523095e-02 | 1.00294860 |
| 300 | 0.030 | 8.7996328224e-02 | 1.00429631 |

The ratio grows slowly above 1 because the discrete Laplacian of a 16-point
wavelength is about 1.3 % weaker than the exact one. With `--skew=0` (upwind) the
ratio drops to 0.965 at step 300: upwind adds numerical dissipation of the same size
as the physical viscosity on this mesh.

The host checks need only python3 and numpy:

```bash
python3 tests/periodic/check_periodic_space.py
```

Part 2 checks the reduced space on 1–4 simulated ranks: the operators agree across
rank counts to 10⁻¹⁶, both are symmetric, the projection leaves `|D u| / |D u**|` at
10⁻¹⁶, and skew advection produces no kinetic energy. Part 3 runs six TGV steps on
1, 2 and 4 ranks and compares the velocity field node by node.

---

## 9. Takeaways

1. **A periodic point is one unknown for every field.** Pressure and velocity
   alike. Collapsing only the pressure leaves two velocity values for one point, and
   no copy afterwards can make both divergence-free.
2. **Put all periodic logic in two maps, P and Pᵀ.** Every operator is then
   `restrict(A(prolong(x)))`, and consistency between the right-hand side, the
   operator and the corrector follows by construction instead of by care.
3. **Order matters in P and Pᵀ.** Copy to slaves before the halo; sum over the halo
   before folding slaves into masters.
4. **Check invariants with one number.** The DOF count must not depend on the rank
   count, and the KE history must agree across rank counts to roundoff.
5. **Model the data layout on the host.** A small numpy model of per-rank slots,
   ghosts and pair tables found every defect above in seconds, without a GPU.

---

### Quick reference

| Role | Symbol | Location |
|------|--------|----------|
| Partner table | `d_periodicPartner` | `mars_periodic_bc.hpp:109` |
| Pairing | `buildPeriodicMap` | `mars_periodic_bc.hpp:925` |
| Cross-rank pair table | `buildCrossRankPeriodicMap` | `mars_periodic_bc.hpp:588` |
| P (prolong) | `PeriodicSpace::prolong` | `mars_periodic_space.hpp:197` |
| Pᵀ (restrict) | `PeriodicSpace::restrict` | `mars_periodic_space.hpp:209` |
| Setup checks | `checkPeriodicPairing` | `mars_periodic_space.hpp:135` |
| Solver | `PeriodicNavierStokes` | `mars_periodic_ns.hpp:398` |
| Time step | `PeriodicNavierStokes::step` | `mars_periodic_ns.hpp:456` |
| Pressure operator `D M⁻¹ G` | `applyProjection` | `mars_periodic_ns.hpp:634` |
| Driver | `main` | `mars_tgv.cu:322` |
| Host checks | `check_periodic_space.py` | `tests/periodic/` |
