# FEM Assembly: From Mesh to Linear System

MARS takes a partitioned GPU mesh to an **assembled sparse linear system**: a degree-of-freedom
(DOF) map, a CSR sparsity pattern, an assembled CSR matrix and a right-hand-side vector, all in
device memory. **Solving that system is your choice**: your own Krylov solver, an external
library (Hypre is supported as an option), or one of the bundled solvers.

This page is the map of that assembly layer: what each piece is, why it exists, and how the data
flows between them. [Quickstart](Quickstart.md) runs the whole pipeline on a generated cube.

---

## The data flow at a glance

```
ElementDomain                     GPU-native mesh: connectivity, node coordinates,
   │                              node ownership (owned / ghost)
   │   node → DOF map
   ▼
node-to-DOF numbering             owned DOFs first [0, numOwned), ghosts after
   │                              (buildDofMappingGpu, two device scans)
   │   sparsity pattern
   ▼
CSR structure                     rowPtr, colInd, diagPtr
   │                              (CvfemSparsityBuilder: 7-nonzero graph or 27-nonzero full)
   │   assemble
   ▼
assembled CSR matrix + RHS        values filled by element-loop kernels that scatter
   │                              the element contributions (atomicAdd)
   ▼
 ───────────────── HANDOFF ─────────────────
   │   your choice
   ▼
solver                            CG / BiCGSTAB / GMRES (bundled), Hypre (optional),
                                  or your own; it consumes CSR + RHS
```

Every stage lives in `backend/distributed/unstructured/fem/` (namespace `mars::fem`). The stages
pass plain device arrays (raw pointers, CSR arrays), so the boundaries are explicit.

---

## 1. Degrees of freedom: the node → DOF map

A continuous linear (P1/Q1) discretization puts **one DOF per mesh node**. Before you can build
a matrix you need a numbering: which equation row does each node own, and which nodes are ghosts
(owned by another rank but needed locally)?

`buildDofMappingGpu` (`mars_cvfem_utils.hpp`) builds the numbering **on the GPU**. Two
`exclusive_scan`s over the owned and ghost masks give a contiguous local numbering on the device:
owned nodes → `[0, numOwned)`, ghosts → `[numOwned, numNodes)`. It returns `numOwned`, the number
of equation rows this rank owns. The CVFEM drivers `mars_cvfem_graph`, `mars_cvfem_graph_tet`,
`mars_cvfem_poisson` and `mars_ex1_poisson` use it. The Navier–Stokes solver numbers its unknowns
with `DofSpace` (`mars_dof_space.hpp`) instead, which also handles periodic points (see the
[Taylor–Green tutorial](periodic_tgv_tutorial.md)).

> **Ownership.** `ElementDomain::getNodeOwnershipMap()` holds one byte per local node: `1` if
> this rank owns the node, `0` for a ghost. On several ranks a node belongs, by default, to the
> rank whose space-filling-curve range contains it, so every node is owned by exactly one rank.
> Only owned nodes get an equation row.

---

## 2. Sparsity: the shape of the matrix

A sparse matrix needs its **structure** (which entries are nonzero) before its values.
`CvfemSparsityBuilder` (`mars_sparsity_builder.hpp`) builds the CSR structure from the element
connectivity, on the GPU, with a Thrust pipeline:

```
edge-list kernel  →  append diagonals  →  remove invalid  →  sort by (row, col)
                  →  unique              →  count per row     →  exclusive_scan (rowPtr)
                  →  binary search of each row for diagPtr
```

In: the element connectivity columns, `numElements`, the node → DOF map and the number of rows.
Out: `rowPtr[numRows+1]`, `colInd[nnz]` and, optionally, `diagPtr[numRows]`; the return value is
`nnz`. Pass `nullptr` for `colInd` to only count `nnz`. You call it twice, once to count and once
to fill: the usual two-pass CSR pattern.

### Two patterns: 7 nonzeros per row (graph) or 27 (full)

- **Graph sparsity (`buildGraphSparsity`)** keeps only **edge-adjacent** couplings: the edges of
  the control-volume dual (12 per hex, 6 per tet). An interior hex node couples to itself and its
  6 axis neighbours: **7 nonzeros per row**. This is the layout of STK / Nalu-Wind.
- **Full sparsity (`buildFullSparsity`)** keeps **every** corner pair of each element (8 × 8 per
  hex, 4 × 4 per tet): **27 nonzeros per row** for an interior hex node.

The graph pattern has about 4 times fewer nonzeros, so it needs about 4 times less memory and
bandwidth in every matrix-vector product. The full pattern keeps the complete element coupling,
for example for a symmetric Poisson matrix solved with CG.

`MARS_SPARSITY_CUB=1` selects a CUB implementation of the graph pattern
(`buildGraphSparsityCub`: radix sort on packed (row, col) keys). The Thrust path stays the
default until the CUB path is validated bit for bit.

---

## 3. The matrix: device CSR

The assembled matrix is **compressed sparse row (CSR)** in device memory. Two representations are
used at different layers:

- **`CSRMatrix`** (`mars_cvfem_hex_kernel.hpp`): a small struct of raw device pointers
  (`rowPtr, colInd, values, diagPtr, numRows, nnz, numOwnedRows`) that the CVFEM kernels receive
  by pointer. `diagPtr` gives the position of each diagonal entry without a search.
  `numOwnedRows` lets kernels skip rows of ghost DOFs. The assembly kernels fill this struct.
- **`SparseMatrix`** (`mars_sparse_matrix.hpp`): an owning device CSR container around the same
  arrays, with `allocate()`, `zero()` (clear the values and keep the pattern, so you can
  re-assemble every time step without rebuilding the structure) and `getDiagonal()`.

---

## 4. Assembly: scatter element contributions into the CSR

Assembly is an **element loop**: one GPU thread (or a team of threads) per element computes the
element's contributions and **scatters** them into the global CSR with `atomicAdd`. The
control-volume operators add each sub-control-surface flux to two nodes with opposite signs, so
interior contributions cancel in pairs and the scatter is locally conservative.

The GPU assemblers per element type:

- **`CvfemHexAssembler` / `CvfemTetAssembler`** (`mars_cvfem_assembler.hpp`,
  `mars_cvfem_tet_assembler.hpp`). Two entry points: `assembleGraphLump(...)` (graph sparsity;
  contributions outside the pattern are lumped onto the diagonal) and `assembleFull(...)` (the
  full element matrix). The tet assembler also has `assembleFullPerip`, which looks up the CSR
  positions once per element instead of once per entry. The hex assembler selects one of twelve
  kernel variants with the `CvfemKernelVariant` enum; they compute the same matrix and differ only
  in performance. See [CVFEM Kernels](CVFEM-Kernels.md).

The Navier–Stokes solver builds on the same control-volume operators (face fluxes, their
divergence, the nodal gradient, the CVFEM Laplacian and the lumped mass `M`); see the
[Poiseuille tutorial](poiseuille_tutorial.md) for how they form a projection method.

A typical re-assembly loop (for example per time step) is: zero the values, call the assembler,
and the values are refreshed against the unchanged sparsity pattern.

---

## 5. The handoff: your solver

After assembly you hold the **CSR matrix and the right-hand side in device memory**. From here
you can use a bundled solver (`backend/distributed/unstructured/solvers/`: CG, BiCGSTAB, GMRES),
Hypre BoomerAMG (with `-DMARS_ENABLE_HYPRE=ON`), or your own.

The drivers apply Dirichlet conditions to the assembled system: `mars_cvfem_poisson.cu` turns
each owned boundary row into an identity row with a zero right-hand side.

On several ranks, a matrix-vector product needs the current values of the ghost DOFs. The
domain's node halo provides them: `exchangeNodeHalo` copies owner values to the ghosts, and
`reverseExchangeNodeHaloAdd` adds ghost contributions to the owners. The bundled CG takes the
exchange as a callback, and the number of owned rows for its global dot products:

```cpp
solver.setOwnedSize(numOwned);   // dot products sum owned rows only, then MPI_Allreduce
solver.setHaloExchangeCallback([&domain, dofMap](cstone::DeviceVector<double>& p) {
    domain.exchangeNodeHalo(p, dofMap);
});
```

---

## A minimal end-to-end example

The shortest real path from mesh to assembled system is the CVFEM graph example,
`examples/distributed/unstructured/mars_cvfem_graph.cu`. Its skeleton:

```cpp
using KeyType   = uint64_t;
using Assembler = CvfemHexAssembler<KeyType, double>;

// 1. mesh (GPU-native, partitioned)
ElementDomain<HexTag, double, KeyType, cstone::execution::Gpu> domain(meshFile, rank, numRanks);
const auto& d_ownership = domain.getNodeOwnershipMap();
size_t nodeCount = domain.getNodeCount();

// 2. node -> DOF (GPU-native local numbering): owned rows first, then ghosts
cstone::DeviceVector<int> d_nodeToDof(nodeCount);
int numOwned = buildDofMappingGpu<KeyType>(d_ownership.data(), d_nodeToDof.data(), nodeCount);

// 3. sparsity, two passes: count nnz, then fill (rows for owned and ghost DOFs)
const auto& c = domain.getElementToNodeConnectivity();   // 8 device columns for a hex
int nnz = CvfemSparsityBuilder<KeyType>::buildGraphSparsity(
    /* 8 connectivity columns */ ..., elementCount, d_nodeToDof.data(), int(nodeCount),
    d_rowPtr.data(), nullptr, nullptr);
CvfemSparsityBuilder<KeyType>::buildGraphSparsity(
    /* same */ ..., d_rowPtr.data(), d_colInd.data(), d_diagPtr.data());

// 4. assemble the values and the right-hand side
Assembler::assembleGraphLump(/* connectivity, coordinates, fields, DOF map, ownership, */ d_matrix, d_rhs.data(), config);

// 5. HANDOFF: CSR + d_rhs -> your solver
```

---

## What MARS provides and what you bring

| Stage | MARS provides | Notes |
|-------|---------------|-------|
| Mesh and partition | yes | GPU-native, SFC-partitioned |
| Node → DOF map | yes | `buildDofMappingGpu`, two device scans |
| Sparsity (CSR structure) | yes | 7-nonzero graph or 27-nonzero full, on the GPU |
| Assembled CSR matrix + RHS | yes | CVFEM assembly kernels (graph-lumped or full) |
| Linear solver | optional | your own, the bundled CG/BiCGSTAB/GMRES, or Hypre |
| Preconditioner | optional | Jacobi (from the diagonal), or Hypre BoomerAMG |

The DOF numbering, the sparsity and the assembly all run **on the GPU**: the mesh is built on the
device and stays there.

---

## See also

- [Quickstart](Quickstart.md): build and run `mars_cvfem_graph`.
- [CVFEM Kernels](CVFEM-Kernels.md): the hex assembly kernel variants.
- [Poiseuille tutorial](poiseuille_tutorial.md): a CFD solver built on these operators.
