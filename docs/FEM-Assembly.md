# FEM Assembly — From Mesh to Linear System

MARS gives you the pipeline from a partitioned GPU mesh to an **assembled sparse
linear system** — a degree-of-freedom (DOF) map, a CSR sparsity pattern, an
assembled CSR matrix, and a right-hand-side vector, all in device memory. **Solving
that system is your choice**: plug in your own Krylov solver, an external library
(Hypre is supported optionally), or one of the simple bundled solvers. The boundary
MARS owns and guarantees is *the assembled CSR + RHS*; everything downstream is
swappable.

This page is the conceptual map of that assembly layer — what each piece is, why it
exists, and how the data flows between them. It is written to be read top-to-bottom
the first time and used as a reference afterward.

---

## The data flow at a glance

```
ElementDomain                     GPU-native mesh: connectivity, node coords,
   │                              per-node ownership (owned / ghost)
   │   node → DOF map
   ▼
node-to-DOF numbering             owned DOFs first [0, numDofs), ghosts after
   │                              (buildDofMappingGpu, two device scans)
   │   sparsity pattern
   ▼
CSR structure                     rowPtr[numDofs+1], colIndices[nnz], diagPtr[numDofs]
   │                              (CvfemSparsityBuilder: 7-NNZ graph or 27-NNZ full)
   │   assemble
   ▼
assembled CSR matrix + RHS        values[nnz] filled by element-loop kernels that
   │                              scatter local element matrices (atomicAdd)
   ▼
 ───────────────── HANDOFF ─────────────────
   │   your choice
   ▼
solver                            CG / BiCGSTAB / GMRES (bundled), Hypre (optional),
                                  or bring your own — consumes CSR + RHS
```

Every stage lives in `backend/distributed/unstructured/fem/` (namespace
`mars::fem`). The stages pass **plain device arrays** (raw pointers / CSR triplets),
so the boundaries are explicit and language-agnostic.

---

## 1. Degrees of freedom: the node → DOF map

A continuous (H1) Lagrange discretization puts **one DOF per mesh node**. Before you
can build a matrix you need a numbering: which equation row does each node own, and
which nodes are ghosts (owned by another rank but needed locally)?

The numbering is built **on the GPU**, with `buildDofMappingGpu`
(`mars_cvfem_utils.hpp`). Two `exclusive_scan`s over the owned/ghost masks produce a
contiguous local numbering entirely on the device — owned nodes → `[0, numOwned)`,
ghosts → `[numOwned, numNodes)` — with no host round-trip beyond reading back the two
final counts. Every CVFEM and Navier–Stokes driver uses this path
(`mars_cvfem_graph.cu`, `mars_amr_*`, the NS solvers). It returns `numOwned`, the
number of equation rows this rank owns.

The lightweight `H1FESpace` class (`mars_h1_fe_space.hpp`) ties this to the
assemblers — it counts the owned DOFs and exposes the reference element. It is the
H1 space wrapper, not a separate map owner.

> **Ownership convention.** The mesh ownership map has three states — `0` ghost, `1`
> owned, `2` shared (owned and on a halo boundary). For assembly, **owned and shared
> both count as owned** (they get an equation row); only pure ghost (`0`) does not.
> This `(== 1 || == 2)` test appears throughout the kernels.

---

## 2. Sparsity: the shape of the matrix

A sparse matrix needs its **structure** (which entries are nonzero) before its
values. `CvfemSparsityBuilder` (`mars_sparsity_builder.hpp`) builds the CSR
structure from element connectivity, on the GPU, via a Thrust pipeline:

```
edge-list kernel  →  append diagonals  →  remove invalid  →  sort by (row,col)
                  →  unique (dedupe)    →  histogram rows   →  exclusive_scan (rowPtr)
                  →  binary-search each row for diagPtr
```

The **I/O contract** is uniform: in = element connectivity columns + `numElements`
+ the node→DOF map + `numDofs`; out = `rowPtr[numDofs+1]`, `colIndices[nnz]`
(pass `nullptr` on a first pass to just *count* `nnz`), and an optional
`diagPtr[numDofs]`; the return value is `nnz`. You typically call it twice — once to
count, once to fill — which is the standard two-pass CSR idiom.

### Two patterns: 7 NNZ/row (graph) vs 27 NNZ/row (full)

- **Graph sparsity (`buildGraphSparsity`)** records only **edge-adjacent**
  couplings — the sub-control-surface (SCS) edges of the control-volume dual (12 per
  hex, 6 per tet). An interior hex node ends up coupling to itself + its 6 axis
  face-neighbors ≈ **7 nonzeros per row**. This matches the STK / Nalu-Wind layout.
- **Full sparsity (`buildFullSparsity`)** records **every** corner-pair within each
  element (8×8 = 64 per hex, 4×4 = 16 per tet) ≈ **27 nonzeros per row** for hexes.

Why two? **Graph + diagonal lumping is what scales** — ~4× fewer
nonzeros means ~4× less memory and bandwidth, which is what makes billion-element
runs feasible on a GPU cluster. **Full** is for consumers that need the complete
element-local coupling (e.g. a symmetric Poisson stiffness solved with CG).

There is also a CUB-backed fast path (`buildGraphSparsityCub`, opt-in via
`MARS_SPARSITY_CUB=1`) that replaces the Thrust sort/unique/scan with
`cub::DeviceRadixSort` / `DeviceSelect` / `DeviceScan` for a ~2–3× speedup on the
sparsity phase at hundreds of millions of nonzeros. The Thrust path is the default
until CUB is validated bit-for-bit.

---

## 3. The matrix: device CSR

The assembled matrix is **compressed sparse row (CSR)** in device memory. MARS has
two CSR representations, used at different layers:

- **`CSRMatrix`** (`mars_cvfem_hex_kernel.hpp`) — a lightweight device struct of raw
  pointers (`rowPtr, colInd, values, diagPtr, numRows, nnz, numOwnedRows`) passed *by
  pointer* into the CVFEM kernels. It carries `__device__ addValue()` and an O(1)
  `diagPtr`-based diagonal accessor (so a kernel can write the diagonal without
  searching the row). The `numOwnedRows` field lets kernels skip scatters into
  ghost-only rows. This is what the assembly kernels fill directly.
- **`SparseMatrix`** (`mars_sparse_matrix.hpp`) — an owning, allocatable device CSR
  container around the same `rowOffsets / colIndices / values` arrays, with
  convenience operations: `allocate()`, `zero()` (clear values, keep the pattern — so
  you can re-assemble each timestep without rebuilding the structure), `getDiagonal()`,
  and raw-pointer accessors for kernels.

---

## 4. Assembly: scatter element matrices into the CSR

Assembly is an **element loop**: one GPU thread (or team) per element computes that
element's local matrix and **scatters** its entries into the global CSR with
`atomicAdd`. Because the control-volume operators telescope (each SCS flux is added
to two nodes with opposite sign), this scatter is locally conservative by
construction.

MARS provides GPU assemblers per element type:

- **`CvfemHexAssembler` / `CvfemTetAssembler`** — the control-volume finite-element
  assembly. Two entry points: `assembleGraphLump(...)` (graph sparsity + diagonal
  lumping — the scalable path) and `assembleFull(...)` (full element-local matrix).
  The tet assembler's `assembleFullPerip` variant pre-resolves the CSR positions
  once per element to avoid a per-(i,j) search in the hot loop — the recommended
  production tet path. The hex assembler additionally dispatches a dozen kernel
  variants (tensor-product, shared-memory, coloring, tensor-core) selected by a
  `CvfemKernelVariant` enum — pure performance variants of the same math, all running
  on the device.

The Navier–Stokes solvers build on these same control-volume operators (the discrete
divergence `D`, gradient `Dᵀ`, and lumped mass `M`) — see the
[Poiseuille channel-flow tutorial](poiseuille_tutorial.md) for how they compose into a projection
method.

A typical re-assembly loop (e.g. per timestep) is: `matrix.zero()` → call the
assembler → the values are refreshed against the unchanged sparsity pattern.

---

## 5. The handoff: your solver

After assembly you hold the **CSR matrix + RHS vector in device memory** — and that
is the contract. From here you use whatever solver you like: a bundled
CG / BiCGSTAB / GMRES, an external library (Hypre BoomerAMG is supported optionally
via `MARS_ENABLE_HYPRE`), or your own. Nothing in the assembly layer hard-requires a
particular solver.

Boundary conditions are applied by marking the essential (Dirichlet) DOFs during
assembly so the system handed off already reflects them.

For multi-rank solves, a distributed sparse mat-vec reconciles ghost contributions
after each product; MARS exposes the ghost scatter-add
(`scatterAddGhostData`, GPU-aware MPI) for exactly this. That is the solver↔matrix
interaction point, not part of assembly proper.

---

## A minimal end-to-end example

The shortest real path from mesh to assembled-and-handed-off system is the CVFEM
graph example,
[`examples/distributed/unstructured/mars_cvfem_graph.cu`](../examples/distributed/unstructured/mars_cvfem_graph.cu).
Its skeleton:

```cpp
// 1. mesh (GPU-native, partitioned)
ElementDomain<HexTag, double, uint64_t, cstone::GpuTag> domain(meshFile, rank, nRanks);

// 2. node -> DOF (GPU-native local numbering)
int numOwned = buildDofMappingGpu(domain.getNodeOwnershipMap(), d_nodeToDof);

// 3. sparsity (two passes: count, then fill) -> CSR structure
int nnz = CvfemSparsityBuilder<uint64_t>::buildGraphSparsity(
              conn..., numElements, d_nodeToDof, numDofs,
              d_rowPtr, /*colInd=*/nullptr, d_diagPtr);   // count
CvfemSparsityBuilder<uint64_t>::buildGraphSparsity(
              conn..., numElements, d_nodeToDof, numDofs,
              d_rowPtr, d_colInd, d_diagPtr);             // fill

// 4. wrap as a device CSR and assemble (fills values + RHS)
CSRMatrix<double> A{d_rowPtr, d_colInd, d_values, d_diagPtr, numDofs, nnz, numOwned};
CvfemHexAssembler<HexTag, double>::assembleGraphLump(/* fields, coords, */ &A, d_rhs, ...);

// 5. HANDOFF: A (CSR) + d_rhs -> your solver of choice
```

That is the whole MARS responsibility: after step 4 you hold an assembled CSR system
in device memory, ready for whatever solver you bring.

---

## What MARS owns vs. what you bring

| Stage | MARS provides | Notes |
|-------|---------------|-------|
| Mesh + partition | ✓ | GPU-native, SFC-partitioned (see [Mesh Reading and Partitioning](Mesh-Reading-and-Partitioning.md)) |
| Node → DOF map | ✓ | GPU-native (`buildDofMappingGpu`, two device scans) |
| Sparsity (CSR structure) | ✓ | 7-NNZ graph or 27-NNZ full; GPU (Thrust or CUB) |
| Assembled CSR matrix + RHS | ✓ | CVFEM assembly kernels (graph-lumped or full) |
| **Linear solver** | optional | **Bring your own**, or use bundled CG/GMRES, or Hypre |
| **Preconditioner** | optional | Jacobi (trivial from the diagonal), or external (Hypre AMG) |

The whole pipeline above — DOF numbering, sparsity, and assembly — runs **on the
GPU**: the mesh is built on-device and stays there to be numbered, patterned, and
assembled. That is what lets the assembly scale to very large meshes on a GPU
cluster.

---

## See also

- [Mesh Reading and Partitioning](Mesh-Reading-and-Partitioning.md) — how the mesh
  upstream of assembly is built and split across GPUs.
- [SFC Mapping](SFC-Mapping.md) — the space-filling-curve numbering that DOF locality
  is built on.
- [Multi-Rank Support](Multi-Rank-Support.md) — the distributed/MPI side, including
  ghost exchange.
- [Poiseuille channel-flow tutorial](poiseuille_tutorial.md) — a worked CFD application built on this
  assembly layer (the CVFEM Navier–Stokes solver).
