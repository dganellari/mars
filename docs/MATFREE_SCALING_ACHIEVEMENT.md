# Matrix-free HO-CVFEM — scaling achievement (2026-06-22)

GPU-native high-order matrix-free CVFEM (Knaus Alg 2) on Alps GH200, cstone branch.

## Operator (validated)
- High-order sum-factorized hex operator, **p=1–7**, bit-exact host+GPU; gates A–D pass.
- **~p+1 convergence**, flat throughput (~4–8 GDOF/s, p=7≈p=1), memory crossover 336× at p=7.

## Distributed scaling
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
  **647M DOF/GPU clean** (cube216); 721M OOMs (per-point metric `d_G` ~7.2 KB/elem is the wall).
- **Why it matters**: high-order is the trillion-DOF path that *fits* Alps. At the measured 647M
  DOF/GPU, **10¹² DOF = ~1,550 GPUs (~390 nodes)** — under half the allocation, using the general
  per-point metric (works for curved/pump meshes). The affine-cube metric (512× less `d_G`) would
  push to ~1–2B DOF/GPU → ~500–1000 GPUs. (vs ~9,300 GPUs for the p=1 element trillion.)
- Tests: `mars_ho_dist_dof_test.cu` (numbering+halo), `mars_ho_dist_apply_test.cu` (operator +
  device-halo timing). Figure: `docs/figures/fig_ho_scaling.png`.

## Halo fix (64M was crashing)
- Root cause (via `MARS_HALO_DEBUG` instrument): **`MPI_ERR_TRUNCATE`** — a per-rank
  ownership over-claim in `NodeHaloTopology` (reverse exchange), *not* GPUDirect.
- Fix: **receiver-driven send-list symmetrization** in `buildFromCstoneHalos`
  (`A.send[B] ≡ B.recv[A]` by construction, build-time only, ownership untouched).
- Verified: 64M runs clean, gates A·1=0 / A·linear=0 at **1e-18** over 498M rows,
  comm 19.2%, pump path no-regression.

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
