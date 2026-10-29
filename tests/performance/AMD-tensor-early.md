# AMD tensor assembly with earlier row scatter

Status: integrated into the HIP Tensor dispatch; HIP and CUDA regression checks passed.
CUDA retains the previous tensor kernel. The previous row-per-thread
experiment passed correctness but was 2.56 times slower in job 687417. This candidate
keeps one thread per element, twelve face evaluations and one scatter per matrix entry.
The P7 MFMA operator (AMD-mfma.md) is a separate optimization.

## Discrete algebra and scheduling

The scalar steady CVFEM equation is div(rho*u*phi - gamma*grad(phi)) = source.
Here mdot is mass/time through the oriented internal surface; gamma is the diffusion
coefficient, phi the nodal scalar and beta the dimensionless upwind reconstruction
factor. Source/boundary/time terms remain outside this kernel. Each internal flux
enters the neighboring control volumes with opposite signs.

Faces are visited in order 0,3,8,1,9,2,10,11,4,7,5,6. Row r is scattered after its
last incident face, at stages 2,4,6,7,9,10,11,11 respectively. A compile-time assertion
checks every face occurs once and every row is released at its last contribution.
At most five rows have begun but not finished, instead of retaining all eight rows
until the end. The compiler determines the actual register lifetimes.

For J[physical][reference], grad_x(N) = J^-T grad_reference(N), hence

    grad_x(N) dot area = grad_reference(N) dot (J^-1 area).

The kernel computes J^-1 area once per face instead of materializing all eight
three-component physical gradients. It uses the same Jacobian, inverse, shape tables,
area and diffusion sign. This changes floating-point association, not the discrete
operator. Upwind beta/gradient and scatter ownership/DOFs are loaded when needed;
this avoids retaining their arrays but can increase input reads. There are no new
barriers, shuffles or shared-memory allocations.

Owned-row bounds, negative DOFs, periodic DOF-equality diagonals, sparse CSR lookup,
additive output and the supplied stream retain the baseline semantics. Face reordering
also changes floating-point summation order; bitwise equality is not required.

## Evidence and run

`python3 scripts/test_tensor_early.py` passes the shared 36 exact-source host
cases under ASan/UBSan: affine/warped adjacent hexes, variable fields, positive/negative/
zero flux, periodic aliases, excluded DOFs, ownership/ghost row bounds, dense/sparse CSR,
additive repeats, empty input and block sizes 128/256/512. Dense all-owned cases check
internal RHS balance and matrix column sums. Separate syntax checks cover HIP/CUDA
stream selection. These are not GPU execution results.

The local LLVM19 gfx942 stub compilation with 64-bit connectivity gives baseline/
candidate scratch 872/792 bytes per lane, VGPR spills 317/266, both 128 VGPR and
occupancy 4. This modest local result is not a ROCm7.2 or runtime measurement; the
previous experiment demonstrated that lower scratch can still mean slower assembly.

Job 687425 passed all 36 entrywise GPU cases and the public hex16 norm checks on
1/4 ranks with block sizes 128/256/512. On public hex192, one MI300A rank on nid002672,
block size 256 and 20 assembly iterations, three alternating repetitions gave:

| | Trial 1 (ms) | Trial 2 (ms) | Trial 3 (ms) | Median (ms) |
|---|---:|---:|---:|---:|
| Previous tensor | 30.16892 | 30.04540 | 29.87233 | 30.04540 |
| Early scatter | 25.77398 | 25.65225 | 25.77512 | 25.77398 |

This is 1.166 times faster assembly (14.2 percent less time), excluding output
zeroing, not a full-solver speedup. ROCm compiler diagnostics give scratch storage
1912 to 1084 bytes/lane and VGPR spills 770 to 494; both use 128 VGPRs and report
occupancy 4. Spills remain.

Both public Tensor assembler entry points select early scatter on HIP. CUDA's
preprocessed dispatch remains unchanged. The kernel itself is byte-identical to
the tested candidate. MARS_ENABLE_TENSOR_BASELINE_CHECK builds the previous kernel
as mars_cvfem_graph_tensor_baseline for regression comparisons; normal production
selection needs no option. The early fixture builds when HIP and tests are enabled.
The rejected row-per-thread experiment is excluded from the release source.

Job 687942 validated production HIP dispatch: full build, 36 entrywise GPU cases,
public hex16 comparisons on 1/4 ranks and all three block sizes, MFMA parity,
48 release tests, 40 other tests (12 skipped), and installed consumers including
the GPU example on four ranks. Release cleanup removes experimental dispatch and
makes the same fixture standalone; preprocessing checks confirm the same HIP/host
fixture and production dispatcher, and both optimized kernels are byte-identical.

CUDA job 5104760 passed the full build, DMMA checks, all 48 release tests, 40 other
tests (12 skipped), assembly-kernel and weak-scaling checks, and installed consumers
including the GPU example on four ranks. Its source manifest matches the cleaned
release implementation. These results validate the tested integration; they do not
establish full-solver performance or parity on arbitrary geometries.
