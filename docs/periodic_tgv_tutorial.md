# Periodic Boundary Conditions in the MARS Navier–Stokes Solver

### A beginner's tour of the Taylor–Green Vortex case, and how a multi-rank periodic bug was found and fixed

This tutorial is for a developer who is new to FEM/CVFEM and to domain
decomposition, but who needs to understand how MARS handles **periodic boundary
conditions** in its GPU incompressible Navier–Stokes (NS) solver, and the story
of a real multi-rank bug that took a while to corner.

You should be able to read it in about 20 minutes. Code snippets are quoted from
the repo with `file:line` references so you can jump straight to the source.

The relevant files:

| File | What lives there |
|------|------------------|
| `backend/distributed/unstructured/fem/mars_periodic_bc.hpp` | All periodic-pair logic: the partner table, the P / Pᵀ kernels, the cross-rank MPI exchange |
| `backend/distributed/unstructured/fem/mars_ns_solver.hpp` | The NS time-stepper: lumped mass, the pressure operator, the velocity solve, the corrector |
| `backend/distributed/unstructured/solvers/mars_cg_solver.hpp` | The Conjugate-Gradient (CG) linear solver and its callback hooks |
| `examples/distributed/unstructured/mars_tgv.cu` | The TGV driver: initial condition + time loop |

---

## 1. What is the Taylor–Green Vortex, and what does "periodic" mean physically?

The **Taylor–Green Vortex (TGV)** is the standard validation case for an
incompressible Navier–Stokes solver. You start with a known, smooth, swirling
velocity field on a cube and let viscosity slowly drain its energy. The kinetic
energy decay over time is well studied, so you can check your solver against a
reference curve. It is the "hello world" of turbulence-resolving solvers.

In MARS the cube is `[0, 2π]³` and the initial field is set analytically in the
driver (`mars_tgv.cu:64-88`):

```cpp
d_u[i] =  V0 * sx * cy * cz;
d_v[i] = -V0 * cx * sy * cz;
d_w[i] = RealType(0);
d_p[i] = (rho * V0 * V0 / RealType(16)) *
         (cos(RealType(2) * x) + cos(RealType(2) * y)) *
         (cos(RealType(2) * z) + RealType(2));
```

**"Periodic" means the box wraps around.** The face at `x=0` is not a wall — it
*is* the same physical surface as the face at `x=1` (here, `x=2π`). A vortex that
swirls off the right edge re-enters from the left edge, unchanged. The same holds
for y and z. There are no walls anywhere; all six faces are periodic.

Two everyday analogies:

- **Pac-Man.** Walk off the right edge of the screen, reappear on the left. The
  two edges are the same place.
- **Rolling the cube into a torus.** Glue `x=0` to `x=1`, then `y=0` to `y=1`,
  then `z=0` to `z=1`. The flat cube becomes a 3-D torus with no boundary at all.

```
   A 2-D slice of the periodic box (the "seam" is where x=0 meets x=1):

        y=1  +----------------------+
             |                      |
             |     flow that        |
   x=0  ===> |     leaves here  ----+---> ... re-enters at x=0
   (seam)    |                      |    (same physical face!)
             |                      |
        y=0  +----------------------+
            x=0                    x=1
            \_____ these two edges are ONE place _____/
```

Because there are no walls, there are **no Dirichlet boundary conditions** on
velocity. The driver sets this up explicitly (`mars_tgv.cu:383`):

```cpp
s.bcKind = NSStepper<KeyType, RealType>::BCKind::Periodic;
```

One consequence: the pressure is only defined up to an additive constant (a
"null space"). MARS removes it by subtracting the global mean each step instead
of pinning a node — see `removeMean` in the driver (`mars_tgv.cu:446`).

---

## 2. The mesh's role: what is a "node", and what is the periodic seam?

The cube is filled with hexahedral (brick) elements. The corners of those bricks
are **nodes**. A node is just a point in space where we store one unknown per
field (one `u`, one `v`, one `w`, one `p`). The solver's job is to find the
field values at every node.

Now look at the seam. A node sitting on the `x=0` face and a node sitting on the
`x=1` face, at the same `(y, z)`, are — physically — **the same point**, because
the box wraps. But in the mesh data they are stored as **two separate nodes with
two separate storage slots**. This is the single most important idea in the
whole tutorial:

> **Two storage slots, one physical unknown.**

MARS gives these two slots names:

- The node on the **min** face (`x=0`) is the **master**.
- The node on the **max** face (`x=1`) is the **slave**.

A device array, `d_periodicPartner`, records the pairing. For every slave it
stores the index of its master; for everything else it stores `-1`
(`mars_periodic_bc.hpp:62-63`):

```cpp
// d_periodicPartner[slave_node] = master_node_index (within local node array,
// including ghosts). For non-slaves: -1.
```

How are pairs found? Each node carries a **space-filling-curve (SFC) key** from
the cornerstone-octree library that MARS is built on. You can mostly ignore SFC
keys for this tutorial — just know they are a unique integer label per node, and
that the master and slave of a periodic pair have **distinct** SFC keys (the box
is periodic but the coordinates are *not* shifted, so `x=0` and `x=1` stay
genuinely different points to the octree). The pairing itself is done by spatial
search: a slave on `x=1` looks for the master on `x=0` with the same `(y, z)`
(`mars_periodic_bc.hpp:147-281`, the `classify` / `packMasterKeys` /
`matchSlaves` kernels).

```
   Master/slave pairing across the x-seam:

        master (x=0)                       slave (x=1)
        d_periodicPartner = -1             d_periodicPartner = master_idx
            o                                   o
            |  same (y,z)  ===  one physical point  ===  |
            o                                   o
```

Corner nodes are special: a node at `(1,1,1)` is a slave in all three directions
at once, so several slaves can chain to a single master. The
`flattenPartnerChainKernel` (`mars_periodic_bc.hpp:270-281`) rewrites those
chains so one lookup reaches the ultimate master.

---

## 3. DOF collapse: one unknown per periodic pair

If master and slave are one physical unknown, the linear system should have **one
equation for them, not two**. This is the same idea MFEM and deal.II use:
*eliminate the slave*. MARS does it through a `node_to_dof` map — a slave node's
"DOF" (degree of freedom = a row/column index in the linear system) is set to its
**master's** DOF. The reduced system has size = number of interior + master DOFs.
From the header (`mars_periodic_bc.hpp:11-14`):

```
//   2) During DOF numbering, slave nodes do NOT get their own DOF; their
//      "DOF" is the master's DOF. The resulting linear system size = number
//      of owned interior+master DOFs ...
```

```
   Two nodes, one DOF:

      node array:   [ ... master_node ...  slave_node ... ]
                          |                     |
      node_to_dof:        v                     v
                        DOF 42  <------------- DOF 42   (same number!)
                          |
      solver vector:    x[42]   (the single shared unknown)
```

### The prolongation P and its transpose Pᵀ

To move between the **node** picture (two slots) and the **DOF** picture (one
unknown) we need two operators that are exact transposes of each other.

**P (prolongation): master → slave copy.** Take the one DOF value and write it
into *both* the master slot and the slave slot, so the per-node kernels see a
consistent field across the seam. This is `periodicBroadcastKernel`
(`mars_periodic_bc.hpp:311-320`):

```cpp
template<typename RealType>
__global__ void periodicBroadcastKernel(const int* d_partner, size_t numNodes,
                                        RealType* d_field)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    int master = d_partner[i];
    if (master < 0) return;
    d_field[i] = d_field[master];   // slave := master
}
```

**Pᵀ (restriction): slave → master sum.** Going the other way, contributions that
landed on the slave slot must be *added* onto the master, and the slave slot
zeroed. This is the transpose of "copy master to slave". MARS has two flavors:

`periodicPairSumKernel` — the gated, single-rank-safe sum used after assembly
(`mars_periodic_bc.hpp:296-309`):

```cpp
__global__ void periodicPairSumKernel(const int*     d_partner,
                                      const uint8_t* d_ownership,
                                      size_t         numNodes,
                                      RealType*      d_field)
{
    ...
    int master = d_partner[i];
    if (master < 0) return;
    if (d_ownership[master] != 1) return;   // only when master is owned here
    atomicAdd(&d_field[master], d_field[i]);
    d_field[i] = RealType(0);
}
```

`periodicFoldToMasterKernel` — the *ungated* sum-only fold used inside the
matrix-free reduced operator, because there the master may be a ghost
(`mars_periodic_bc.hpp:330-340`):

```cpp
// Exact TRANSPOSE of periodicBroadcastKernel: for every slave (partner>=0),
// field[partner] += field[slave]; field[slave] = 0.
__global__ void periodicFoldToMasterKernel(const int* d_partner, size_t numNodes,
                                           RealType* d_field)
{
    ...
    int master = d_partner[i];
    if (master < 0) return;
    atomicAdd(&d_field[master], d_field[i]);
    d_field[i] = RealType(0);
}
```

**The operator the solver actually inverts is `Pᵀ A P`.** You take the reduced
DOF vector, prolong it to nodes with P, apply the bare element operator A, then
restrict back with Pᵀ. Because P and Pᵀ are exact transposes, `Pᵀ A P` stays
symmetric — which CG requires. This is exactly how `applyDDTReduced` is built
(`mars_ns_solver.hpp:5784`):

```
//   p_dof  --P-->  phiNode  --A(bare)-->  ApNode  --P^T(sum-only)-->  Ap_dof
```

---

## 4. The CVFEM pipeline, stage by stage — and what periodicity needs at each

MARS solves incompressible NS with **Chorin projection**. Conceptually each time
step has four stages:

1. **Predictor** — advance velocity ignoring the
   incompressibility constraint, giving a provisional velocity `u**`.
2. **Implicit velocity diffusion solve** — solve a Helmholtz-like system per
   velocity component (this is the **velocity CG**) to handle viscosity
   implicitly.
3. **Pressure Poisson solve** — find a pressure correction `φ` that, when its
   gradient is subtracted, makes the velocity divergence-free.
4. **Corrector** — subtract `(dt/ρ) ∇φ` from `u**` to get the final
   divergence-free `u^{n+1}`.

The discrete identity the projection enforces is:

```
   D ( u** − (dt/ρ) M⁻¹ Dᵀ φ ) = 0
```

where:

- **D** is the discrete **divergence** operator (turns a velocity field into a
  per-node divergence).
- **Dᵀ** is its transpose, the discrete **gradient** operator (turns pressure
  into a force).
- **M** is the **lumped mass** matrix: a diagonal matrix whose entry at a node is
  the control volume around that node. `M⁻¹` just divides by that volume.

Rearranging, `φ` solves the pressure-Poisson system with operator
`A = D M⁻¹ Dᵀ`. This is the "DDT" operator (the default for periodic TGV; see
`mars_tgv.cu:245`).

### Why every stage must agree on the *same* discrete D

This is the subtle part. The RHS of the pressure solve uses D (the divergence of
`u**`). The operator `A = D M⁻¹ Dᵀ` uses D and Dᵀ. The corrector subtracts
`(dt/ρ) Dᵀ φ`. **If any of these uses a slightly different discrete D or Dᵀ, the
projection does not close** — the corrected velocity is not actually
divergence-free, and the error feeds back step after step. The driver is careful
to use the literal `Dᵀ` everywhere (`mars_tgv.cu:417`):

```cpp
s.useLegacyGradient = false;
```

with the comment:

```
// This makes the predictor's grad(p^n) consistent with the
// DDT projection operator D M^{-1} D^T, so div(u^{n+1})=0 holds algebraically
// exactly.
```

### Lumped mass M and the periodic seam

`M` is built per node by scattering each element's volume to its corner nodes,
then summing. For periodic pairs, the master and slave together bound the same
physical control volume, so **the seam node must carry the combined pair mass**.

Here is the tricky bit that caused the first bug (see §6). After assembly,
`maybePeriodicSum` folds the slave's mass onto the master *and zeroes the slave
slot*. That zero is poison: the gradient-normalize kernel divides by mass, and a
zero mass makes it **zero out the gradient at the slave node**
(`normalizeGradientPerNodeKernel`, `mars_ns_solver.hpp:1645-1646`):

```cpp
RealType m = lumpedMassNode[i];
if (m == RealType(0)) { gx[i] = gy[i] = gz[i] = RealType(0); return; }
```

So MARS explicitly **mirrors** the combined mass back onto the slave after the
DOF gather (`mars_ns_solver.hpp:3041-3049`):

```cpp
if (s.periodicMap)
{
    mars::fem::periodicBroadcastKernel<<<nBlocks, s.blockSize>>>(
        s.periodicMap->d_periodicPartner.data(), s.nodeCount, s.d_massNode.data());
    cudaDeviceSynchronize();
    // Re-sync ghost slots: a cross-rank master-ghost slave just got the
    // master value; keep ghost copies consistent for the assembler.
    s.domain.exchangeNodeHalo(s.d_massNode);
}
```

The comment above it (`mars_ns_solver.hpp:3023-3040`) spells out the failure mode
in full: a zeroed seam mass makes the operator's internal D drop the seam flux,
so only about *half* the divergence gets removed.

### The pressure operator as Pᵀ A P, matrix-free

The pressure operator is never assembled as a matrix. It is applied **matrix-free**
inside the CG, recomputed every iteration, as the `Pᵀ A P` chain shown in §3
(`applyDDTReduced`, `mars_ns_solver.hpp:5807-5871`). The key correctness rule,
repeated throughout the file: the RHS, the operator, and the corrector must all
use the **identical** discrete D, or the projection identity `D u^{n+1} = 0` will
not hold.

---

## 5. Single rank vs multiple ranks — the escalation

Everything above works on **one GPU (one rank)** without any MPI. On one rank a
periodic pair is *same-rank*: the master node is owned locally, and `node_to_dof`
simply maps the slave to the master's DOF. The collapse is trivial — both node
slots literally share one solver unknown.

```
   SINGLE RANK (master and slave both local):

        master_node ----\
                          >-- node_to_dof --> DOF 42   (one unknown, done)
        slave_node  ----/
```

Now **split the cube across N GPUs** (domain decomposition). Each rank owns a
chunk of the mesh. The master of a periodic pair can now live on a **different
rank** than its slave. This is the **cross-rank seam**. Enforcing "slave equals
master" now requires the two ranks to *talk* — that means MPI.

```
   MULTI-RANK (the cross-rank seam):

      rank A (owns x=1 face)           rank B (owns x=0 face)
      -----------------------          -----------------------
        slave_node (owned)               master_node (owned)
        master_ghost  <----- cstone halo delivers the image ----- (master)
              |                                |
        node_to_dof = its OWN dof        node_to_dof = its own dof
        (NOT collapsed!)                 (the real shared unknown)
```

The communication primitive MARS leans on is the **cstone node halo**:

- `exchangeNodeHalo` — owner → ghost. Copies each owned node's value into the
  ghost copies that other ranks hold.
- `reverseExchangeNodeHaloAdd` — ghost → owner, summing. Adds each ghost's
  accumulated contribution back into the owner.

The design goal is to **minimize bespoke MPI** and use the halo for as much as
possible. The driver turns the periodic box on at the octree level so cstone
itself delivers opposite-face nodes as halo ghosts (`mars_tgv.cu:352-354`):

```cpp
amr.initialize(meshFile, rank, numRanks,
               /*periodicAxesMask=*/ 7,
               RealType(boxLo), RealType(boxHi));
```

But there is a catch that the rest of this tutorial is about: **the halo cannot
bridge a periodic pair by itself.** The master and slave have *distinct SFC keys*
(no coordinate shift), so cstone does not consider them the same shared node. The
halo will deliver the master's *image* as a ghost, but it will not automatically
sum the slave's contribution onto the master. That gap is what the cross-rank
periodic code (and the bug below) is all about. From the header
(`mars_periodic_bc.hpp:66-70`):

```
// reverseExchangeNodeHaloAdd cannot
// bridge this because slave and master have DISTINCT SFC keys by design
// (no coord shift); they are NOT the same cstone-shared node.
```

---

## 6. The bug story

This is the heart of the tutorial. It is a detective story, and the moral is:
**stop guessing, measure.**

### The symptom

On **1 rank**, periodic TGV ran beautifully: the divergence of the final
velocity was at roundoff, around `1e-14`. The projection was closing perfectly.

On **4 ranks**, the same case **blew up**. The pressure magnitude `|p|` exploded
from a healthy `0.13` to `7447` within a handful of steps. Same physics, same
mesh, more ranks — diverging.

### The discipline: build on-device probes

Rather than theorize about which stage was wrong, the fix was to add **on-device
probes** that measure the two things that *must* be true if the projection works.
They are gated behind `MARS_PROJ_PROBE` so they cost nothing when off
(`mars_ns_solver.hpp:7031`):

```cpp
static const bool g_projProbe = (std::getenv("MARS_PROJ_PROBE") != nullptr);
```

- **P1** checks the *linear solve*: `||A·φ − b|| / ||b||`. If this is tiny, CG
  genuinely solved `A φ = b`, so any remaining failure is *not* in the solve
  (`mars_ns_solver.hpp:8010-8019`):

  ```cpp
  applyDDTReduced<KeyType, RealType, ElementTag>(s, phi_dof, Aphi, pn, apn, g1, g2, g3);
  axpyDofKernel<RealType><<<dofBlocks, s.blockSize>>>(
      rdif.data(), Aphi.data(), b.data(), RealType(-1), n);
  ...
  std::cout << "  [PROJ-P1] |b|=" << nb << " |A*phi|=" << nA
            << " |A*phi-b|/|b|=" << (nr / (nb > 0 ? nb : RealType(1)))
  ```

- **P3** checks the *projection*: `||D u^{n+1}|| / ||D u**||`. This is the ratio
  of divergence *after* the corrector to divergence *before*. It should be ≈ 0
  (projection removed the divergence). A ratio of ≈ 0.5 means only half the
  divergence was removed (`mars_ns_solver.hpp:8508-8510`):

  ```cpp
  std::cout << "  [PROJ-P2] |Du**|=" << dIn << " |Du^{n+1}|=" << dOut
            << "  [PROJ-P3] ratio=" << (dOut / (dIn > 0 ? dIn : RealType(1)))
  ```

### Fix 1 — the seam mass mirror

The first measurement pointed at the lumped mass. The slave's mass had been
zeroed (§4), so `normalizeGradientPerNodeKernel` zeroed the slave's gradient,
so the operator's internal D dropped the seam flux. The result: P3 ≈ 0.5 — the
projection removed only **half** the divergence at the seam, and pressure ran
away.

The fix was the **seam mass mirror** shown in §4 (`mars_ns_solver.hpp:3041-3049`):
mirror the combined pair mass back onto the slave so every mass consumer — the
operator's M⁻¹, the corrector's M⁻¹, the RHS normalize, the DDT diagonal, and the
probe — sees the same non-zero seam mass. This killed the catastrophic blowup
(`|p|` no longer went to thousands).

### The pivotal redirect — "we were debugging the wrong stage"

Fixing the mass helped, but multi-rank still drifted. The breakthrough was a
direct **single-rank vs 4-rank baseline**: run the identical step on both and
compare the numbers entering each stage.

The result was damning: the divergence of `u**` *entering* the pressure solve was
already **0.0096 on 4 ranks vs 1e-5 on 1 rank**, and the velocity diffusion
solve's iteration count `cg_uvw` was **19–20 on 4 ranks vs 2–3 on 1 rank**. The
TGV handoff notes (`TGV_HANDOFF.md`) record exactly this kind of `cg_uvw`
breakdown and the eventual blow-up "around step 30".

The conclusion: **the velocity diffusion solve was already broken *before* the
pressure stage ever ran.** All the earlier effort had been aimed at the pressure
operator — the wrong stage.

### Root cause — the half-strength seam rows (the symcheck)

A second on-device probe, `MARS_PERIODIC_XR_SYMCHECK`
(`mars_ns_solver.hpp:3906-3911`), inspects the **assembled velocity matrix**
across a cross-rank seam. It compares `A[slave_row, master_ghost_col]` on the
slave's rank against `A[master_row, slave_ghost_col]` on the master's rank. If
the matrix were symmetric across the seam, these would match.

They did not. The symcheck found `A_local = 0` and `A_peer ≈ 0` at the seam: the
master row did **not** carry the slave-side stiffness, and vice versa. The
assembler had emitted **two independent half-strength rows** — one per side of
the seam — that *never coupled across ranks*. The reason is the same SFC-key
issue from §5: the assembler could not route the periodic-image column entries
across ranks, so each rank built only its local half of the equation. The comment
records it (`mars_periodic_bc.hpp:1066-1070`):

```
// the matrix has half-strength rows on the seam (assembler couldn't
// route the periodic-image column entries cross-rank), so the resulting
// Ap slot is half what the merged equation requires.
```

### The fix — Strategy B: mirror the operator matrix-free, every CG iteration

The clean fix does **not** try to assemble the missing cross-rank entries into
the matrix. Instead it keeps both seam rows half-strength and **reconstructs the
merged operator on the fly, inside CG**, using one cross-rank primitive:
`crossRankPeriodicPairSumDof` (`mars_periodic_bc.hpp:1074-1186`). It does a
two-leg MPI exchange — slave → master sum, then master → slave broadcast — so
after every matvec the two half-rows add up to the full equation on both sides.

CG exposes a hook for exactly this: `setSpmvPostCallback`, called right after each
sparse matvec (`mars_cg_solver.hpp:239-246`):

```cpp
// Ap = A * p
spmv(A, p, Ap);

// Post-SpMV cross-rank coupling enforcement (periodic-pair Ap sync,
// etc). Without this, cross-rank periodic slave+master rows hold
// independent residual values and pAp can go negative.
if (spmvPostCallback_) {
    spmvPostCallback_(Ap);
}
```

The velocity solve wires it up (`mars_ns_solver.hpp:5146-5151`):

```cpp
solver.setSpmvPostCallback(
    [pm, n2d](cstone::DeviceVector<RealType>& Ap)
    {
        mars::fem::crossRankPeriodicPairSumDof<KeyType, RealType>(*pm, n2d, Ap);
    });
```

There is one more subtlety. After the pair-sum makes `slave == master` on two
different ranks, a normal dot product would count that **one** physical unknown
**twice** (once per rank). So CG also takes a **dot mask** that skips the slave
copy (`mars_cg_solver.hpp:367-371`):

```cpp
// Owned entries with mask==1 are skipped in dot(). Used for cross-rank
// periodic pairs: after the post-SpMV pair-sum makes slave==master on two
// different ranks, dot() would count that one physical unknown twice. Mask
// out the slave copy so each merged DOF contributes once. nullptr = no mask.
void setDotMask(const uint8_t* mask) { dotMask_ = mask; }
```

wired as (`mars_ns_solver.hpp:5138-5139`):

```cpp
solver.setDotMask(s.d_isPeriodicXRSlaveDof.empty()
                  ? nullptr : s.d_isPeriodicXRSlaveDof.data());
```

The masked dot itself is a small free function (`mars_cg_solver.hpp:26-35`):

```cpp
template<typename RealType>
RealType maskedOwnedDot(const RealType* x, const RealType* y,
                        const uint8_t* mask, int n)
{
    return thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<int>(0),
        thrust::counting_iterator<int>(n),
        [x, y, mask] __device__ (int i) -> RealType {
            return mask[i] ? RealType(0) : x[i] * y[i];
        }, RealType(0), thrust::plus<RealType>());
}
```

### The parity completion — four seam quantities, not two

Strategy B made the operator and the dot product seam-consistent — two of the
quantities CG touches. A longer **110-step** run revealed it was *still* not
done: the run was stable for about **30 steps, then blew up**. The residual was
**compounding**: a small per-step inconsistency that quietly grew until it
dominated.

The diagnosis (`mars_ns_solver.hpp:5154-5164`): the proven *pressure* path makes
**four** seam quantities consistent — the operator/matvec, the dot product, the
**Jacobi preconditioner diagonal**, and the **RHS**. Strategy B had only fixed
the first two. The Jacobi diagonal at a seam DOF was wrong (the slave's mass was
zeroed before the DOF gather, so its `M/dt` diagonal term was 0), which biased
the CG search direction a little every iteration — exactly the kind of error that
compounds.

The cure was to merge the **diagonal** too, using the **same primitive**, fed
through a single-shot override (`mars_ns_solver.hpp:5165-5174`):

```cpp
typename NSStepper<...>::Matrix::DeviceVectorReal diagVelMerged;
if (s.bcKind == ...::BCKind::Periodic
    && s.periodicMap != nullptr
    && !s.periodicMap->cross_.peers_.empty())
{
    diagVelMerged = A.getDiagonal();
    mars::fem::crossRankPeriodicPairSumDof<KeyType, RealType>(
        *s.periodicMap, s.d_node_to_dof.data(), diagVelMerged);
    solver.setPreconditionerDiagOverride(diagVelMerged);
}
```

CG uses the override in place of the matrix diagonal (`mars_cg_solver.hpp:67`,
`79`):

```cpp
void setPreconditionerDiagOverride(const Vector& diag) { diagOverride_ = &diag; }
...
Vector diag = diagOverride_ ? *diagOverride_ : A.getDiagonal();
```

The RHS gets the same pair-sum treatment. With all **four** CG-visible
quantities — operator, dot, diagonal, RHS — made seam-consistent through the one
primitive, the multi-rank run reached **single-rank parity** and the energy decay
matched.

---

## 7. Takeaways for beginners

If you remember five things from this tutorial, make it these:

1. **Periodicity = "two storage slots, one physical unknown."** Every operator
   that reads those slots — mass, divergence, gradient, the matrix rows, the
   preconditioner, the dot product — must agree about the seam, or the math
   silently breaks.

2. **Single-rank can hide a multi-rank bug.** The same code that was perfect on
   1 rank was catastrophically wrong on 4. Always cross-check ranks; a
   single-rank-vs-N-rank baseline is one of the most powerful diagnostics you
   have.

3. **Measure with on-device probes, don't theorize.** `MARS_PROJ_PROBE` (P1 for
   the solve, P3 for the projection) and `MARS_PERIODIC_XR_SYMCHECK` (for matrix
   symmetry) turned a vague "it blows up" into precise numbers that pointed at the
   exact broken stage. The pivotal moment was realizing the velocity solve, not
   the pressure solve, was the culprit.

4. **For a distributed assembled operator, every CG-visible quantity must be made
   seam-consistent — not just the matvec.** The operator, the preconditioner
   diagonal, the RHS, and the dot product all have to agree. Fixing only two of
   the four bought ~30 stable steps before the residual compounded into a blow-up.

5. **Reuse one proven primitive everywhere.** The whole cross-rank fix rides on a
   single building block, `crossRankPeriodicPairSumDof`, applied to the operator,
   the diagonal, and the RHS. One well-tested primitive used in four places beats
   four bespoke patches — it is easier to reason about and far easier to get
   right.

---

### Quick reference: the kernels and hooks

| Role | Symbol | Location |
|------|--------|----------|
| Partner table | `d_periodicPartner` | `mars_periodic_bc.hpp:62` |
| P (master → slave copy) | `periodicBroadcastKernel` | `mars_periodic_bc.hpp:311` |
| Pᵀ (slave → master sum, gated) | `periodicPairSumKernel` | `mars_periodic_bc.hpp:296` |
| Pᵀ (slave → master sum, fold) | `periodicFoldToMasterKernel` | `mars_periodic_bc.hpp:330` |
| Cross-rank merge (per-DOF) | `crossRankPeriodicPairSumDof` | `mars_periodic_bc.hpp:1074` |
| Cross-rank broadcast (per-DOF) | `crossRankPeriodicBroadcastDof` | `mars_periodic_bc.hpp:1193` |
| Reduced pressure operator Pᵀ A P | `applyDDTReduced` | `mars_ns_solver.hpp:5807` |
| Seam mass mirror | (lumped mass build) | `mars_ns_solver.hpp:3041` |
| Velocity solve + CG wiring | `solveOneComponent` | `mars_ns_solver.hpp:5030` |
| CG post-SpMV hook | `setSpmvPostCallback` | `mars_cg_solver.hpp:363` |
| CG dot mask | `setDotMask` / `maskedOwnedDot` | `mars_cg_solver.hpp:371` / `26` |
| CG diagonal override | `setPreconditionerDiagOverride` | `mars_cg_solver.hpp:67` |
| Projection probes | `MARS_PROJ_PROBE` (P1/P3) | `mars_ns_solver.hpp:7031`, `8010`, `8508` |
| Matrix symmetry probe | `MARS_PERIODIC_XR_SYMCHECK` | `mars_ns_solver.hpp:3906` |
