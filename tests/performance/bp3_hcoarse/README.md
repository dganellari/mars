# Geometric coarse levels for BP3

The normal `mars_marsir_bp` driver now includes the validated geometric
hierarchy below p=1. Native integration passed jobs 5281137/5281138 on base
35085f9e. The algorithm passed jobs 5278195/5278196; both results are recorded below.
Compatibility with the newer operator interface passed job 5281607.
The accepted Hypre consolidation adapter and its 50,000-row cap are unchanged.

## Numerical contract

The operator is scalar diffusion, `-div(k grad u) = f`, with homogeneous
Dirichlet boundaries on the procedural cube. The continuous Qp space uses GLL
nodes and p+2 Gauss quadrature. The fine operator, geometry, right-hand side
and accuracy checks are unchanged; only the preconditioner changes.

Each geometric level groups aligned 2x2x2 local hexes. Trilinear interpolation
P maps eight parent corner values to 64 child corner copies. Restriction is
its raw transpose, with no DSS or multiplicity weighting before restriction.
The cached parent element matrix is `P^T diag(K_child) P`. The first matrices
come from the existing generated p=1 action on eight unit columns. Each level
uses two symmetric Chebyshev pre/post steps with a diagonal/Gershgorin bound;
the bottom uses the existing P1Amg settings. Symmetry and positive definiteness
of the complete preconditioner require a symmetric positive bottom action.
The native checks exercise this assumption for the tested Hypre recipe.

All added matrices, vectors and transfers remain on the GPU. Solves reuse the
existing halo exchange and add no global reductions or host vector copies.
Setup uses collective spectral bounds; the inherited P1Amg assembly still
passes remote setup entries through the host. Automatic selection uses a
collective hardware eligibility check once before setup.

## Selection

`--h-levels -1` is the default. It uses the following measured choices only
with the validated SM90, CUDA/Hypre consolidation build and AMG recipe:

| Degree | Ranks | Global cells per axis (`--ne`) | Added levels |
|---|---:|---:|---:|
| 3 | 8 | 128 | 2 |
| 3 | 64 | 128 or 256 | 2 |
| 7 | 64 | 128 | 1 |
| 7 | 8 | 64, with `--solves >= 5` | 1 |

Automatic selection also requires unblocked procedural `bp3`, the default
0.05 deformation, unit stretch, discrete RHS, two pre/post Jacobi-Chebyshev
steps, the MARS interval, no power-step override, three AMG cycles, default
AMG coarsening, adaptive consolidation, GPU-aware MPI and Hypre's own sparse
kernels. Other workloads keep the existing p=1 AMG path. This is an explicit
calibration table for the benchmark driver, not a universal scaling model or
a silent extrapolation to unmeasured meshes, ranks or hardware. The printed
selection receipt identifies the actual depth.

`--h-levels 0` retains the previous path, including for historical MFEM
comparisons. Explicit depths 1, 2 or 3 are available for experiments on aligned
unblocked procedural BP3 at p=3/7 with symmetric smoothing. Each rank must have
at least two cells per axis at the last level. Unsupported requests fail before
solver setup. The h hierarchy changes iteration counts and must be reported
separately from the earlier equal-algorithm comparison.

The p=7 eight-rank exception accounts for its cold cost: the worst measured
setup penalty requires four later solves to recover. One-GPU and small p=7
cases failed the original 5% gain gate and retain the previous default.

## Integration gate

The frozen job uses the exact previous candidate executable from preparation
5278195 as the reference, the same image, Hypre library and 61 PTX files.
Preparation builds the normal driver and another P1Amg consumer with the
consolidation adapter ON/AUTO and OFF, the adapter fixture, and an instrumented
copy of the identical normal source. Only that test target accepts `--h-check`
and runs direct GPU algebra validation; production runs no test instrumentation.

Preparation runs CTest, 36 solver invocations on one/four ranks and six rejected
commands. It checks p=1 fallbacks and explicit depths 1/2/3 against the reference,
including six direct GPU hierarchy checks. A sealed pass enables the dependent
job: 54 interleaved invocations on 1/8/64 GPUs, covering every automatic choice
and the important fallbacks. It compares the chosen depth only; it does not
repeat the search over candidate depths. Both sides use three solves, or five
for the p=7 eight-rank amortization case. The metric is the final solve's maximum
host wall time over ranks, from three alternating rounds. Iterations/hierarchies
must match and median time must stay within 3% in either direction. True residual
must be <=1e-10 and relative FCG measure <=1e-13. Cold costs and all failed samples
remain in the record; passing this gate does not claim additional speedup.

After transfer and checksum verification, the user runs on Daint:

```bash
bash /capstor/scratch/cscs/gandanie/bp3-hcoarse-resume-kit-20261106/hcoarse/submit.sh
```

The recovery script submits one node to finish preparation (thirty-minute limit),
then sixteen nodes for integration comparison (ninety-minute limit), only after
preparation passes. It reuses the six verified binaries, CTest and thirty completed
solves from job 5279173; only six solves and six rejected commands remain. Results are under `bp3-hcoarse-integrated/<job-id>/` in the same scratch.
The uenv is explicit; submission from an existing uenv is handled by
`--uenv-passthrough=ignore`. The old completed kit is left untouched.

Local checks use the actual CLI source, target registration, parser mutations,
preparation seals, complete mocked dispatch including retained failures, and
spooled submission scripts. CUDA host/device syntax is checked with and without
Hypre/consolidation and with test instrumentation. These are local checks, not
native runtime validation. The transfer oracle compiles the production transfer
header; the earlier full host oracle used an exact bottom solve and therefore
is not evidence about GPU Hypre.

## Integration preparation 5279173

All six native builds, CTest and thirty solver invocations passed their checks.
The runner then stopped because it compared total AMG cycle counts between
normal and instrumented runs. The four direct validation applications each use
three bottom AMG cycles: the instrumented run correctly counted 156 instead of
144. The level, rows, nonzeros, FCG iterations and residual histories match.
This was a runner defect, not a failing numerical check.

The correction keeps the complete cycle-count equality for normal/reference
comparisons. For the instrumented preparation only, it requires exactly twelve
additional cycles when consolidation is active, or zero when inactive. It also
requires the direct validation result and identical hierarchy and iterations.
Raw counters are preserved. Negative tests reject any other delta or a changed
hierarchy. All thirty completed solves were independently reparsed; worst true
residual is 7.683e-12 and worst direct check is 2.65e-16. Six downloaded binaries,
257 raw artifact hashes, 66 process receipts and 51 rank audits were verified.

The recovery kit verifies the original source snapshot, binaries and completed
receipts before copying them into a new result directory. Numerical sources,
PTX, dependencies and test instrumentation must match exactly. It makes no edits
to the failed job or its raw evidence and does not recompile. The recovery completed as job 5281137 and its timing comparison as 5281138,
with the results below. The original failure record is retained.

## Completed normal-build integration

Recovery 5281137 and comparison 5281138 pass. All 36 preparation results and
54 timed invocations were independently reparsed, including six direct GPU
checks and six rejected commands. The recovery reused thirty earlier solves
and all six binaries; it did not rebuild them. Verified 337 sealed preparation
files, six downloaded executable hashes, 1,782 process receipts, 1,767 rank
audits and 3,888 raw artifacts. Exact commands, summaries and all nine verdicts
reproduce. Rank placement is identical within each timing case.

The largest true residual is 8.28e-11, relative FCG measure 9.122e-14 and direct
hierarchy error 2.65e-16. Iterations and hierarchies match the validated candidate.
The largest median timing difference is 0.34%, within the 3% equivalence gate.
No sample was removed. Whole warm solve medians in milliseconds:

| Case | Candidate reference | Normal executable |
|---|---:|---:|
| p3, 1 GPU, fallback | 163.780 | 163.822 |
| p7, 1 GPU, fallback | 241.088 | 241.602 |
| p3 weak, 8 GPUs, h2 | 209.504 | 209.702 |
| p7 short, 8 GPUs, fallback | 312.145 | 311.741 |
| p7 repeated, 8 GPUs, h1 | 284.204 | 284.238 |
| p3 weak, 64 GPUs, h2 | 232.490 | 232.465 |
| p7 weak, 64 GPUs, h1 | 321.701 | 321.590 |
| p3 strong, 64 GPUs, h2 | 83.296 | 83.579 |
| p7 small, 64 GPUs, fallback | 101.738 | 101.521 |

These compare the same hierarchy and operator in two executables, establishing
integration equivalence. They do not measure additional speedup over stock AMG.
The separate algorithm experiment below remains the source of that claim.
The later constant-memory operator/loader change bf087b34 was not in this frozen
build; combined-source compatibility is a separate check.

## Completed GH200 results

Preparation 5278195 and comparison 5278196 pass all 162 solver invocations,
30 direct GPU hierarchy checks at 1/4/8/64 ranks, CTest and six invalid-CLI
checks. The largest true residual is `8.28e-11`, relative FCG measure `9.122e-14`,
and direct hierarchy error `2.87e-16`. The three downloaded executable hashes,
332 preparation evidence files, 4,416 process receipts and 4,416 rank audits
were independently checked. All summaries, command lines and timing verdicts
were reproduced from 9,367 downloaded files. Rank placement stayed the same
within each comparison case. Builds report no warnings or unused options.

Whole warm solve medians in milliseconds, three interleaved rounds:

| Case | GPUs | Accepted solver | Private h=0 | h=1 | h=2 | h=3 |
|---|---:|---:|---:|---:|---:|---:|
| Weak p=3, ne=64 | 1 | 163.730 | 162.968 | 173.595 | 172.915 | 172.035 |
| Weak p=7, ne=32 | 1 | 239.038 | 241.079 | 243.806 | 242.459 | 241.067 |
| Weak p=3, ne=128 | 8 | 221.944 | 222.000 | 235.959 | **209.147** | 214.443 |
| Weak p=7, ne=64 | 8 | 312.028 | 312.076 | **284.318** | 290.005 | 292.154 |
| Weak p=3, ne=256 | 64 | 260.989 | 260.208 | 259.468 | **229.676** | 232.368 |
| Weak p=7, ne=128 | 64 | 347.625 | 347.836 | **318.028** | 320.277 | 324.440 |
| Strong p=3, ne=128 | 64 | 94.747 | 94.712 | 86.092 | **83.199** | 86.400 |
| Strong p=7, ne=64 | 64 | **101.244** | 102.242 | 98.522 | 103.591 | 107.388 |

Bold candidates pass the predefined median gain and separated-range gates;
the accepted solver is retained in the final row because h=1's gain is below
the 5% threshold. None of the single-GPU candidates passes. Several slower
multi-GPU alternatives also pass the gate, but the highlighted candidate is
the fastest accepted warm result for that case. The private control matches
accepted iteration counts and hierarchies, with at most 0.99% median difference.
Outliers remain in the raw samples; none was removed.

At 64 GPUs, h=2 improves p=3 weak solves by 12.0% and strong solves by 12.2%
against the accepted executable. The weak case uses 14 FCG iterations instead
of 16; its complete coarse call remains about 6.1 ms. Thus this gain mainly
comes from fewer iterations, not a faster coarse call. For weak p=7, h=1
improves solve time by 8.5%, with 15 instead of 16 iterations and the coarse
call falling from 4.069 to 3.492 ms. These coarse timings include the entire
h hierarchy and bottom AMG; they are not pure communication or AMG timings.

Using the accepted one-GPU solver and these selected 64-GPU candidates, whole
warm-solve weak-scaling efficiency is approximately 71% for p=3 and 75% for
p=7, versus 63% and 69% for the accepted solver in this batch. Work is fixed
in elements per GPU. These efficiencies include the iteration-count changes;
they do not claim ideal scaling or a constant per-iteration cost.

Cold cost matters. For weak p=7 on eight GPUs, h=1's median setup plus first
solve is 1,010 ms against 935 ms for h=0, while each later solve saves about
28 ms. The median difference amortizes after three additional warm solves;
using the least favorable observed pairing needs four. On the 64-GPU weak
cases, the selected candidates also reduce median setup plus first solve.
All memory snapshots retain more headroom than their h=0 controls; none is a
peak-memory measurement. A single solve or a different local problem size
therefore needs its own choice rather than a blanket new default.

These are MARS-versus-MARS comparisons with a changed preconditioner. They do
not update the historical equal-algorithm MFEM comparison. Normal-build integration
and the current operator compatibility checks pass. These results do not establish
performance outside the measured degrees, sizes, rank counts,
geometry and environment.

## Compatibility with the current operator interface

While the integration jobs ran, bf087b34 changed the generated operator's launch
arguments and moved the one-dimensional matrices to module constants. Its loader
and twelve BP PTX files must travel together. The public C++ `apply` signature is
unchanged, and the coarse code still supplies fixed matrices on one stream, but
the older frozen binary/old PTX test does not execute this new interface.

The h implementation was committed locally and rebased over those three incoming
commits without changing their files. Local CUDA syntax passes with the newer
loader. A separate two-node check builds only the normal and instrumented drivers
and runs 24 solves on 1/4/8 GPUs, including eight direct hierarchy checks, p=1 and
small-problem fallbacks, and automatic versus explicit depth selection. It uses
the current loader and its matching twelve PTX files with the same Hypre/image.
It does not repeat the completed 64-GPU equivalence benchmark or claim new scaling
numbers. Current-source performance must remain distinct from the older measurements.

After transfer, the user runs:

```bash
bash /capstor/scratch/cscs/gandanie/bp3-hcoarse-compat-kit-20261106/hcoarse/submit-compat.sh
```

One two-node job has a two-hour limit and writes under
`/capstor/scratch/cscs/gandanie/bp3-hcoarse-compat/<job-id>/`. Native combined-source
compatibility passed job 5281607. The historical exporter rejects the new loader
with old PTX; `export_compat.py` freezes the matching current pair.

All 24 solves and eight direct hierarchy checks pass on 1/4/8 GH200 GPUs.
The largest true residual is `2.098e-11` (limit `1e-10`), relative FCG measure
`9.835e-14` (limit `1e-13`), and direct hierarchy error `2.649236e-16`.
The two downloaded binary hashes, sealed artifacts, exact commands, process
receipts and rank audits were independently verified, and every summary was
reproduced from raw output. Automatic and explicit selections agree; normal
and instrumented runs agree after accounting for the four diagnostic applies.
The builds have no compiler warnings or unused CMake options. This closes the
combined-source correctness gate; it does not remeasure 64-GPU performance.
