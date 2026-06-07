# CVFEM Assembly Kernels — GPU Optimization

The element-loop assembly (see [FEM Assembly](FEM-Assembly.md)) is the bandwidth-
and-atomics-bound heart of a CVFEM solver. MARS ships **many GPU kernel variants that
compute the same hex assembly** but trade differently against the GPU's real
bottlenecks. They are selected at runtime through a `CvfemKernelVariant` enum
(`mars_cvfem_assembler.hpp`), so you can pick the strategy that fits your hardware
without changing the physics.

This page explains *what each kernel does differently and why* — organized by the
bottleneck it attacks — so you can choose one deliberately. Every variant computes
the identical operator (a 12-sub-control-surface flux balance per hex: an advection
flux and a diffusion flux scattered into the global CSR); they differ only in *how*
they get there.

> **A naming caveat to know up front.** The variant called `Tensor` is **not** a
> tensor-product / sum-factorization kernel — it is a straightforward scalar
> thread-per-element kernel that assembles the full 8×8 local matrix. "Tensor" here
> means "assemble the full local tensor," and it is the staging ground for the
> tensor-*core* kernels. The only kernels that actually factor the operator into a
> matrix multiply are the WMMA/WGMMA ones. (For trilinear p=1 hexes the classical
> sum-factorization FLOP argument barely applies anyway.)

---

## The bottlenecks (the "why" axis)

Every kernel below is a different point on these trade-offs:

- **Atomic scatter contention.** Neighboring elements share nodes, so threads
  `atomicAdd` into the same CSR entries. Both the contention *and* the per-entry CSR
  position search cost. → addressed by coloring (remove atomics) or warp-cooperative
  reduction before scattering.
- **Register pressure vs occupancy.** The full local matrix is `double lhs[64]` —
  ~64 registers that often spill and crush occupancy. Eliminating it is the single
  biggest occupancy lever. → addressed by the "low-register" and "Perip" kernels.
- **Node gather traffic.** Each thread gathers 9 nodal fields at 8 *scattered* node
  indices. → addressed by AoS packing or a block-level shared-memory node cache.
- **CSR position lookup in the hot loop.** Finding `(row, col)` in CSR per entry. →
  addressed by sorted-column binary search, or pre-resolving all positions once.
- **Matrix size (memory/bandwidth).** Hex CVFEM's natural stencil is 27 columns/row;
  the graph stencil keeps only 7. → a *solver-side* memory tradeoff (see last
  section), not an assembly-speed one.

> The performance figures quoted in the kernel headers (occupancy %, register
> counts, L2-traffic ratios) are the authors' profiling notes and design estimates,
> not independently verified benchmarks. Treat them as intent, and measure on your
> own hardware.

---

## 1. Baseline and incremental upgrades

- **`Original`** — the reference thread-per-element kernel. Computes area vectors on
  the fly, carries the full `lhs[64]`, and scatters with a **linear scan** of each
  CSR row. Correct and readable; the slowest. Use it as the ground truth.
- **`Optimized`** — same structure, but uses the fact that the shifted shape function
  has only two nonzero entries per sub-control-surface (so interpolation is a 2-term
  average, not an 8-term loop), unrolls the loops for instruction-level parallelism,
  and replaces the linear scatter scan with a **binary search** on sorted columns. A
  cheap, safe upgrade over the baseline. Still register-heavy.

## 2. Lowering register pressure (occupancy)

- **`Shmem`** — eliminates `lhs[64]`. It keeps only `rhs[8]` and `diag[8]`
  accumulators and **scatters each off-diagonal contribution immediately** into
  pre-resolved CSR positions; only the diagonal is accumulated locally. It also reads
  **precomputed** area vectors instead of recomputing them. The header reports a large
  register drop and occupancy gain. (Despite the name it uses no shared memory — it is
  a low-register / direct-scatter kernel; the name is historical.) The tradeoff is
  **many more atomics** — it wins when atomics are cheap (modern hardware atomics) and
  occupancy was the limiter.
- **`Team`** — a **warp-per-element** mapping: 32 threads cooperate, lanes process the
  12 sub-control-surfaces in parallel and reduce the local matrix in **shared memory**
  before scattering, so per-thread state collapses and there are fewer *global*
  atomics. The tradeoff is that lanes 12–31 idle during the surface loop (low lane
  utilization) and shared-memory capacity caps elements per block.

## 3. Reducing node-gather traffic

- **`TensorAoS`** — numerically identical to `Tensor`, but reads the 9 nodal fields
  from a **packed, aligned `NodeData` struct** instead of 9 separate arrays. The node
  read is a scatter/gather (8 random node indices per thread); AoS packing turns ~9
  scattered memory sectors per node into ~3, for a large cut in gather traffic. Costs
  a one-time packing pass and a second copy of node data in memory; the scatter side
  is unchanged.
- **`SmemCache`** — the most elaborate memory kernel. A block of (Morton-ordered)
  elements touches far fewer *unique* nodes than node-references; this kernel builds a
  **block-level deduplication cache in dynamic shared memory**, loads each unique
  node's fields and CSR row metadata into shared memory **once**, then serves every
  thread's reads from shared memory (~5 cycles) instead of L2 (~30 cycles). It needs a
  large shared-memory carveout (>48 KB), and it is the most fragile production-track
  kernel: it relies on spatial (Morton) ordering keeping the unique-node count per
  block under a fixed cap.

## 4. Avoiding the race a different way: coloring

- **`TensorColored`** — byte-for-byte the `Tensor` kernel **with every `atomicAdd`
  replaced by `+=`**, launched **once per color**. A host-side greedy graph coloring
  guarantees that within one color no two elements share a node, so the scatters are
  provably race-free and atomics are **eliminated entirely**. This is the big win on
  older GPUs where FP64 atomics are slow. The cost is ~8–12 sequential kernel launches
  (less parallelism per launch) plus a host coloring precompute. On hardware with fast
  FP64 atomics the benefit shrinks — which is why the atomic-based `Tensor` is the
  claimed current best there, not `Colored`.

## 5. Pre-resolving CSR positions: the "Perip" trick

- **`TensorPerip`** (and `TensorPeripLb2`) — instead of accumulating `lhs[64]` and
  scattering at the end, it **pre-looks-up all 64 CSR positions once** (binary search,
  before the surface loop) and **scatters each contribution immediately** to the known
  position. The register trade is explicit: drop the 64-double local matrix, add a
  64-int position array — a net register reduction — and the expensive search leaves
  the hot loop entirely. The `Lb2` variant additionally forces two blocks per SM; with
  the local matrix gone, forcing that occupancy costs almost nothing. Same tradeoff as
  `Shmem` (more, smaller atomics). `SmemCache` is essentially "Perip + a shared-memory
  node/row cache."

## 6. Tensor cores (experimental)

The diffusion block of the local matrix can be written as a matrix multiply
`Sᵀ · D`, where `S` is the element-independent surface scatter (±1) and `D` is the
element-specific diffusion-coefficient matrix. Two kernels offload that 8×8 onto the
**FP64 tensor cores**:

- **`WmmaTensor`** — one warp per element, using the `wmma` FP64 mma API (SM80+,
  i.e. A100/H100/GH200), with a non-tensor-core fallback below SM80.
- **`WgmmaTensor`** — Hopper-only (SM90a), using asynchronous `wgmma` via inline PTX,
  packing several elements' rows into one wide matmul so the FP64 tensor pipe is better
  fed, overlapping with the scalar advection work.

These are **experimental / capability-demonstration** paths: the 8×8 FP64 tiles are
small (the tensor cores are under-fed for a trilinear element), and the WGMMA path
hand-rolls architecture-pinned PTX. They show MARS can drive the tensor cores from a
CVFEM assembly; they are not the recommended production path on current hardware.

---

## Tetrahedra: a separate, smaller family

Tets (4 nodes, 6 sub-control-surfaces) have only three kernels — `Full`, `FullPerip`,
and `Graph` — and **none** of the tensor / shared-memory / coloring hex variants. Two
things to know:

- **`FullPerip`** is the recommended production tet kernel: it applies the same
  pre-resolved-positions trick as the hex Perip (pre-resolve all 16 CSR positions,
  diagonal via the direct pointer, `__ldg` on read-only loads), and is bit-identical
  to `Full` at double precision.
- The tet **diffusion** operator is the classical constant-gradient linear-tet
  stiffness `K_ij = V·γ·(∇N_i·∇N_j)` (exact one-point quadrature), while the
  **advection** uses the per-surface edge form. So the tet `Graph` kernel and the hex
  `Graph` kernel are **not** the same operator — do not treat them interchangeably.

---

## Graph (7 NNZ/row) vs full (27 NNZ/row): a solver tradeoff

This is independent of which kernel you pick. The hex CVFEM stencil naturally couples
27 columns per row; the **graph** path keeps only the 7-point edge-connected stencil
and **lumps** any element contribution whose column is not in the pattern onto the
diagonal:

- **Full (27 NNZ/row)** — the complete element-local coupling. Use for symmetric /
  Poisson-type systems solved with CG that need every coupling.
- **Graph + lumped (7 NNZ/row)** — ~4× fewer nonzeros → ~4× smaller matrix in memory
  and ~4× less bandwidth in every solver mat-vec. At the billion-element scale, where
  the matrix and the SpMV dominate cost, this is the path that matters. The price is
  a *lumped* (more diagonally-dominant, less accurate) operator.

So "graph vs full" is an **accuracy/memory** choice made at the sparsity-and-assemble
level; the kernel variants above are **how fast** you assemble whichever you chose.

---

## Choosing a kernel

| If you want… | Use |
|---|---|
| Correctness reference | `Original` (hex), `Full` (tet) |
| A solid default on modern GPUs | `Tensor` (hex), `FullPerip` (tet) |
| Maximum occupancy / large meshes | `TensorPerip` / `TensorPeripLb2` |
| Lower node-gather traffic | `TensorAoS` |
| Squeeze L2 node traffic further | `SmemCache` (fragile; needs spatial ordering) |
| No FP64 atomics (older GPUs) | `TensorColored` |
| Smallest matrix / extreme scale | the **graph-lumped** assembly path |
| Drive the tensor cores (experimental) | `WmmaTensor` / `WgmmaTensor` |

Start from the default (`Tensor` / `FullPerip`) and the graph-lumped sparsity, then
profile — the right kernel depends on your GPU's atomic throughput, your mesh's
spatial locality, and whether you are occupancy- or bandwidth-bound.

---

## See also

- [FEM Assembly](FEM-Assembly.md) — the pipeline these kernels sit inside
  (DOF map → sparsity → assemble → handoff).
- [GPU Acceleration](GPU-Acceleration.md) — the broader device-execution model.
- [Pump Flow tutorial](pump_cfd_tutorial.md) — the CVFEM operators in a full solver.
