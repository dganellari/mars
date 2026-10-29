# AMD P7 MFMA operator

The scalar Laplacian applies A = sum_ab D_a^T G_ab D_b with the existing GLL
quadrature, shared DOF gather and additive scatter. G contains the weighted
geometric metric in per-point mode; affine mode multiplies the constant metric
by the three quadrature weights. Boundary elimination remains the caller's job.

On HIP the public per-point and affine launchers use the FP64 MFMA plane kernel
(`mars_ho_laplacian_mfma_direct.hpp`) for P7/double with 64 threads per block.
The CUDA launchers keep their staged WMMA kernels. AMD targets without FP64 MFMA
(other than gfx90a and gfx942) take the staged kernel.

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

For comparison, the H100 staged WMMA kernel takes 0.831-0.833 and 0.688-0.693
ms on the same case (Daint, CUDA 12.9). Both MFMA kernels match the staged
per-point result to within 6.2e-16 (relative). These are operator timings, not
a measured full-solver speedup.

## Local checks

```bash
python3 scripts/test_amd_mfma_direct.py
python3 scripts/test_ho_dmma_dispatch.py
```

The first runs the actual kernel through a host 64-lane emulator with the
measured gfx942 MFMA lane layout, under address and undefined-behaviour
sanitizers, against an independent scalar reference (1, 2, 3 and 5 elements,
both metric modes, constant fields). The second compiles the public wrappers
against recording stubs and checks HIP/CUDA selection, arguments, stream and
error forwarding. Neither replaces compiling and running with the GPU toolchain.
