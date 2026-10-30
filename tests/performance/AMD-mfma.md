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
times slower than H100 on this case.

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
