# Matrix-free HO-CVFEM — scaling achievement (2026-06-22)

GPU-native high-order matrix-free CVFEM (Knaus Alg 2) on Alps GH200, cstone branch.

## Operator (validated)
- High-order sum-factorized hex operator, **p=1–7**, bit-exact host+GPU; gates A–D pass.
- **~p+1 convergence**, flat throughput (~4–8 GDOF/s, p=7≈p=1), memory crossover 336× at p=7.

## Where matrix-free wins (measured p=1, GH200, ncells=256, single GPU)
- `ho_compare` SpMV throughput / bytes-per-DOF:
  - assembled **CVFEM graph (7-nnz)** — production: **24.7 GDOF/s, 88 B/DOF**
  - assembled CVFEM **full (27-nnz CSR)**: 8.74 GDOF/s, 325 B/DOF
  - matrix-free **PA** (store `d_G`): 5.32 GDOF/s, 278 B/DOF
  - matrix-free **MF** (recompute, no stored metric): 0.78 GDOF/s, 185 B/DOF
- **At p=1 assembled CVFEM wins — faster *and* leaner. Do not use matrix-free for p=1
  production.** Matrix-free is a **high-order tool (p≥2–3)**: the stored matrix explodes
  ~(2p+1)³ (**8,748 B/DOF at p=4, 40,500 at p=7**) while matrix-free holds ~100–120 B/DOF
  → **~340× leaner at p=7**, throughput flat across order (sum-factorization, O(p⁴) apply).
- The assembled baseline here is **CVFEM** (shape functions, 7-nnz graph), *not* EBVC —
  EBVC is a separate edge-based scheme (polynomials replace shape functions) per `mars.pdf`.
- This is the **Knaus / Nalu-Wind split**: store the matrix at low order, go matrix-free
  (partial assembly) for high order, on CPU and GPU alike. MARS's one extension is the
  **--MF** path (recompute the metric, nothing stored); Nalu/Knaus always store it.

## Roofline — fair comparison (bytes-moved, not storage)
- Both kernels run at **~71–72% of HBM peak** (comparably tuned). The earlier "matrix-free
  41% of peak" was a **storage-anchoring error**; against the *true moved traffic*
  (~544 B/DOF at p=1) the apply is ~72%.
- The winner is set by **bytes moved** — structural, tuning-independent. **DMMA (FP64 tensor
  cores) is irrelevant to throughput**: the apply is bandwidth-bound at every p≤7 (arithmetic
  intensity 1.2→4, never near the 8.5 FLOP/byte CUDA ridge).
- On paper matrix-free reads *less* than even the 7-nnz CVFEM at p=1 (geometry ~44 B/DOF
  cached < 88), so its roofline *ceiling* is higher — but reaching it is a research gamble
  (irregular gather/scatter vs a solved cuSPARSE SpMV). **The deck does not depend on a p=1 win.**

## Distributed scaling
- **Validated matvec headline (the only validated matvec scaling number): 40 billion DOF / 64
  GH200**, bit-exact `A·1 = 1e-18` (p=2/3/4), apply ~98% weak-scaling. Everything larger below is
  either *setup-scaled* (640B build) or a *projection* (1T) — never a >40B matvec/`A·1` result.
- **Apply weak-scales 95%** to 64M DOF / 32 GPU (procedural cube, held per-GPU).
- **Communication self-resolves** with per-GPU size (surface/volume, V^-1/3):
  comm fraction **51.9 → 41.9 → 33.5 → 28.6 → 19.2 %** at 2 / 8 / 16 / 32 / 64 M elems/GPU.
- Bit-exact comm/compute overlap (forward halo hidden behind interior apply).

## Distributed high-order operator (multi-rank, validated 2026-06-22)
- Distributed HO DOF numbering + ownership + halo: `HODofHandler::buildDistributed`
  (corner = cstone P1 ownership, interior = element-local, **edge/face = min-rank-among-holders**
  via a peer key exchange `resolveHoDofOwnership`) + `HoHalo` (receiver-driven forward/reverse,
  with a **device gather/scatter + GPUDirect** path for scaling).
- Gates on Alps GH200, 4 ranks, p=2/3/4: global DOF count == analytic (no orphans, no
  double-ownership); **A·1 = HO-Laplacian·const = 1.7e-18 / 2.5e-17 / 2.8e-17** (full matvec:
  forward halo → owned-element apply → reverse-add); device halo **bit-exact vs host**.
- Weak scaling p=4, ~4.2M DOF/GPU held: **apply 98%** (1→8 GPU; 5.9 GDOF/s/GPU = the single-GPU
  rate, so going distributed costs the operator nothing); full matvec 78% (blocking halo).
  **Comm self-resolves**: at 8 GPU, 4M→11M DOF/GPU cuts comm 22→15%, full-matvec 78→87%.
  90.5M DOF at p=4 on 8 GPUs, clean.
- **Per-GPU ceiling (measured, single GH200, p=4)**: apply dead-flat at ~6 GDOF/s from 135M to
  **~620M DOF/GPU clean** (cube216); 721M OOMs (per-point metric `d_G` ~7.2 KB/elem is the wall).
- **Why it matters**: high-order is the trillion-DOF path that *fits* Alps. At the measured ~620M
  DOF/GPU, **10¹² DOF projects to ~1,700 GPUs (~430 nodes)** — under half the allocation, using the
  general per-point metric (works for curved/pump meshes). The affine-cube metric (512× less `d_G`)
  would push to ~1–2B DOF/GPU → ~500–1000 GPUs. (vs ~9,300 GPUs for the p=1 element trillion.)
  **This is a projection from the measured per-GPU ceiling, not a run.**
- Tests: `mars_ho_dist_dof_test.cu` (numbering+halo), `mars_ho_dist_apply_test.cu` (operator +
  device-halo timing). Figure: `docs/figures/fig_ho_scaling.png`.

## Halo fix (64M was crashing)
- Root cause (via `MARS_HALO_DEBUG` instrument): **`MPI_ERR_TRUNCATE`** — a per-rank
  ownership over-claim in `NodeHaloTopology` (reverse exchange), *not* GPUDirect.
- Fix: **receiver-driven send-list symmetrization** in `buildFromCstoneHalos`
  (`A.send[B] ≡ B.recv[A]` by construction, build-time only, ownership untouched).
- Verified: 64M runs clean, gates A·1=0 / A·linear=0 at **1e-18** over 498M rows,
  comm 19.2%, pump path no-regression.

## >4.3B-element blocker — diagnosed, fixed, setup-validated at 640B (2026-06-23 overnight)
- **Diagnosis (hypothesis rejected first)**: the crash was *not* a `uint32`/`LocalIndex` overflow.
  It was cstone's replicated global-octree count `MPI_Allreduce` exceeding **2 GiB**: ~805M leaf
  nodes × 4 B = **3.22 GB** in one collective, overflowing Cray-MPICH's chunked-collective 32-bit
  **byte** count (each individual count is a fine uint32; the message size is the wall).
- **Fix**: auto-scale the global `bucketSize` (env `MARS_GLOBAL_BUCKETSIZE`) — a coarser global
  tree means fewer leaf nodes, so the count message stays under 2 GiB. Free for the apply (the
  matvec is entirely local; only the global partition tree gets coarser).
- **Setup-validated at 640B DOF / 1024 GH200**: the overnight run executed the full distributed
  **setup** — domain decomposition + GPU DOF numbering (**numDof = 628,857,939 per rank**) +
  ownership + halo — with **no Allreduce failure** (global bucketSize = 256, `extractDistDof` 4.2 s,
  `buildDistributedGpu` 7.858 s). **This is a SETUP/BUILD milestone, not a matvec.** The run then
  timed out in the host-side `A·1` cross-check before the apply (since fixed by a host-gate skip).
- **Trillion operator BUILT — 1.005T DOF / 1,600 GH200, ~16 s** (numbering + ownership resolve +
  halo, all on device, no crash). **This is the BUILD, not the matvec.** The full matvec is queued
  on **448 nodes (1,792 GPU, 561M DOF/GPU — fits the ~620M `d_G` ceiling)**; it validates with the
  device `A·1 ~ 1e-18` gate.

## Per-GPU ceilings (8-rank, Alps GH200, 95 GiB)
| mode | per-GPU ceiling |
|---|---|
| assembly (CSR) | ~70M elems |
| matrix-free apply | ~77M elems |
| **domain decomposition (DD build: SFC partition + ownership + halo)** | **~108M elems** |

DD is leaner (skips the ~15 GB CSR + ~30 GB area-vector/operator buffers); the build
is sync-dominated (~0.55 GiB/M-elem).

## Trillion-element domain decomposition — reachable on Alps
- **Validated 108.7M elems/rank clean** (`mars_cvfem_scale --build-only`, balanced 0.4%,
  18 GiB headroom, 48.5 s build).
- 10¹² elements ÷ ~108M/GPU ≈ **9,300 GPUs** — fits full Alps (~10.7k GH200).
- Run: `mars_cvfem_scale --ncells=10000 --build-only` → **10¹²-element FEM domain
  decomposition** (SFC partition + node ownership + halo topology).
- Supporting near-trillion FEM numbers on full Alps: assembly ~0.7T, apply ~0.8T.

## Figures
`docs/figures/*.png` — memory crossover, throughput, convergence, weak-scaling,
comm self-resolve (now incl. 64M), and `fig_trillion` (path-to-trillion ceilings).

## Next lever
Ghost-halo balancing: SFC-interior ranks hold ~2× the recv-ghosts of SFC-endpoint
ranks, which sets the per-GPU ceiling. Balancing it raises the ceiling and cuts the
GPU count for 10¹². (Deeper cstone domain-sync memory slim is the further lever.)
