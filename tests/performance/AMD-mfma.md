# AMD P7 MFMA operator

The scalar Laplacian applies A = sum_ab D_a^T G_ab D_b with the existing GLL
quadrature, shared DOF gather and additive scatter. G contains the weighted
geometric metric in per-point mode; affine mode multiplies the constant metric
by the three quadrature weights. Boundary elimination remains the caller's job.

On HIP the public per-point and affine launchers use the FP64 MFMA plane kernel
(`mars_ho_laplacian_mfma_direct.hpp`) for P7/double with 64 threads per block.
On CUDA SM80+ they use the same schedule on `mma.sync` (see H100 below). AMD
targets without FP64 MFMA (other than gfx90a and gfx942) and older NVIDIA GPUs
take the staged kernel.

## Plane schedule

Two elements share a wavefront. Each lane owns one (j, k) position of both
elements, which is where the in-plane MFMA results land, and keeps u and the
output along i in registers. Per i-plane, the j and k contractions and their
transposes run on MFMA through LDS copies of the current plane; the i
contraction, its transpose and the metric are lane-local. LDS holds three 8x8
planes per element: 3072 bytes per block.

The previous kernel kept u and the three gradients of both elements in LDS,
32768 bytes per block. Profiling showed why that was slow: only about two
wavefronts fit on a CU (LDS limit) and they waited 81-87% of the time.

## Evidence

Beverin, ROCm 7.2.3, MI300A, public cube P7/E32: 32768 elements, 11390625 DOFs,
one GPU, 30 applies per run, three runs, timings include output zeroing.
Run time in ms per apply:

| Kernel | Per-point | Affine | LDS/block | VGPR + AGPR | Waves/CU |
| --- | ---: | ---: | ---: | ---: | ---: |
| Staged | 3.083-3.086 | 2.482-2.537 | 35840 | 48 + 8 | |
| Previous MFMA | 1.301-1.306 | 0.946-0.966 | 32768 | 48 + 216 | 1.9 |
| Plane MFMA | 0.591-0.595 | 0.346-0.356 | 3072 | 108 + 28 | 11-14 |

Both MFMA kernels match the staged per-point result to within 6.2e-16
(relative). These are operator timings, not a measured full-solver speedup.

## H100

The same schedule runs on NVIDIA SM80+ with `mma.sync` m8n8k4
(`mars_ho_laplacian_dmma_plane.hpp`): one warp per element, each lane owning its
two C-fragment positions of the plane, 1536 bytes of shared memory per warp.
Daint, CUDA 12.9, GH200, same case, 4 warps per block (2 and 8 within 2%):

| Kernel | Per-point | Affine | Shared/block | Registers | Occupancy |
| --- | ---: | ---: | ---: | ---: | ---: |
| Staged WMMA | 0.828-0.834 | 0.688-0.693 | 35840 (2 warps) | 48-56 | 18.75% |
| Plane | 0.426-0.430 | 0.213-0.214 | 6144 (4 warps) | 176-182 | 12.5% |

On H100 the plane kernel is register-limited, not faster by occupancy: its
per-point apply reaches 62% of DRAM throughput against 30% for the staged one.
With the same schedule on both GPUs, MI300A is 1.37 (per-point) and 1.58 (affine)
times slower than H100 on this case, and 1.28 and 1.60 in the round below. Bank
conflicts do not explain it: padding the shared-memory rows from 8 to 10 doubles
removes all of them (`SQ_LDS_BANK_CONFLICT` 1.57e7 to 0 per dispatch), but the
per-point apply takes 0.558-0.563 ms against 0.569-0.572 ms and the affine one
0.340-0.367 ms against 0.329-0.337 ms (3 runs each).

## Where the MI300A time goes

A timing-only ablation (not kept in the tree) replaced the scatter of the plane
kernels on both GPUs: same case, 3 runs each, ms per apply, run time (2026-10-31).
The plain `+=` on interior points matches the staged result (5.2e-16 on H100,
6.2e-16 on MI300A); store and none give wrong results.

| Scatter | H100 per-point | MI300A per-point | H100 affine | MI300A affine |
| --- | ---: | ---: | ---: | ---: |
| atomic (shipped) | 0.430-0.431 | 0.549-0.552 | 0.213-0.214 | 0.341-0.343 |
| plain `+=` on interior points | 0.479-0.482 | 0.564-0.570 | 0.265-0.271 | 0.338-0.341 |
| plain store | 0.400-0.402 | 0.470-0.478 | 0.178-0.183 | 0.209-0.212 |
| none (all work, no writes) | 0.336-0.338 | 0.388-0.392 | 0.150-0.155 | 0.173-0.176 |

- Without the atomics, MI300A is 1.17-1.18 times slower than H100 with plain stores
  and 1.15-1.16 times without writes, on both paths.
- The atomics take 0.077 ms (per-point) and 0.132 ms (affine) more than plain stores
  on MI300A, and 0.029 and 0.034 ms on H100. That is 40% of the per-point gap and
  76% of the affine one. The work without writes is 44% and 18% of it.
- On MI300A every atomic goes past the L2 (`TCC_EA0_ATOMIC` = `TCC_ATOMIC` = 6.50e6
  per apply), as in the graph assembly (`AMD-assembly.md`). On H100 they are L2
  reductions, and the L2 serves 89% of their 1.14e7 sectors (ncu `lts`).
- No other element shares an element's interior points (216 of its 512), so they need
  no atomic. Adding them with a plain `+=` is slower on H100 (12% per-point, 25%
  affine). On MI300A it is 3% slower per-point and the same affine: the atomic requests
  fall only 17% (to 5.42e6), and the reads past the L2 rise 6% (8.56e6 to 9.03e6).

Open: an atomic-free scatter, as the row map is for the graph assembly. On MI300A its
saving is bounded by the gap between the atomic scatter and plain stores (0.08 ms
per-point, 0.13 ms affine), before its own cost. These runs do not explain the
remaining 1.15-1.18 times.

## Local checks

```bash
python3 scripts/test_amd_mfma_direct.py
python3 scripts/test_cuda_dmma_plane.py
python3 scripts/test_ho_dmma_dispatch.py
```

The first two run the actual kernels through host lane emulators (64 lanes with
the measured gfx942 MFMA layout; 32 lanes with the m8n8k4 layout verified on
GH200), under address and undefined-behaviour sanitizers, against an independent
scalar reference (up to 9 elements, both metric modes, constant fields). The third
compiles the public wrappers against recording stubs and checks the HIP/CUDA
selection with and without FP64 matrix instructions, arguments, stream and error
forwarding. Neither replaces compiling and running with the GPU toolchain.
