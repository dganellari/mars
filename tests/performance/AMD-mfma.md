# AMD P7 MFMA integration

The scalar Laplacian applies A = sum_ab D_a^T G_ab D_b with the existing GLL
quadrature, shared DOF gather and additive scatter. G contains the weighted
geometric metric in per-point mode; affine mode multiplies the constant metric
by the three quadrature weights. Boundary elimination remains the caller's job.
No pressure, flow, time integration or halo algorithm changes here.

On HIP the public per-point and affine launchers use the direct FP64 MFMA kernel
for P7/double with 64 threads per block. The CUDA launchers retain their staged
WMMA kernels. Unsupported AMD instruction targets trap rather than produce an
empty result. Supported device bodies remain gfx90a and gfx942.

Two elements share a wavefront. MFMA inputs and results go directly between
registers and the input/gradient arrays; intermediate A/B/C tiles and the shared
copy of D are removed. The working storage is 32768 instead of 35840 bytes/block.
The odd final element has a dummy partner that does not scatter. The caller must
zero output when required; launches remain asynchronous on the supplied stream.

## Evidence

Beverin job 687389 tested the direct experimental entry point with ROCm 7.2.3.
Public cube P7/E32: 32768 elements, 11390625 DOFs, one GPU, 30 applies per run,
three runs. Timings include output zeroing, not a solver iteration.

| Metric | Staged median | Direct median | Ratio |
| --- | ---: | ---: | ---: |
| Per-point | 3.086 ms | 1.314 ms | 2.35 |
| Affine | 2.446 ms | 0.931 ms | 2.63 |

GPU parity, nullspace, symmetry and constant-metric linear-field cancellation
passed for E1/2/3. Host lane emulation also covers varying SPD metrics and shared
DOFs under address/undefined-behavior sanitizers. These are bounded fixtures,
not proof for all geometries or a measured full-solver speedup.

The integration adds public-launcher tests for both metric modes, additive
output, non-default stream ordering, a varying metric and empty HIP partitions.
These public-launcher checks passed in jobs 687397 and 687942. The existing test
and benchmark programs are the only callers of this public P7 API in this tree;
no full flow solver is switched to it.

Job 687397 measured staged/direct medians of 3.085/1.302 ms (per-point) and
2.452/0.930 ms (affine), on the same public P7/E32 configuration: 2.37 and 2.64
times faster respectively. HIP integration passed the full build, release tests,
other tests and installed-consumer checks. CUDA job 5104760 also passed the full
build, DMMA, 48 release tests, 40 other tests (12 skipped), kernel/weak-scaling checks
and installed consumers on one/four ranks. Tensor assembly is described in
AMD-tensor-early.md; rejected experiments are excluded from this release.

## Local checks

```bash
python3 scripts/test_amd_mfma_direct.py
python3 scripts/test_ho_dmma_dispatch.py
```

The first runs actual device indexing through a host 64-lane emulator. The
second compiles the actual public wrappers against recording stubs and verifies
HIP/CUDA selection, arguments, stream and error forwarding. Neither substitutes
for compilation and execution with the GPU toolchain.

## Cluster regression

HIP jobs 687397 and 687942 validate the selected implementation. CUDA job 5104760
validates the cleaned release source; source hashes were checked against the local
tree. The production kernels are byte-identical to the HIP-tested kernels, and
preprocessing checks confirm unchanged dispatch after experimental code removal.
The CUDA implementation is unchanged, but the shared test sources were rebuilt.
All launches used one node; no multi-node or full-flow performance claim follows.
