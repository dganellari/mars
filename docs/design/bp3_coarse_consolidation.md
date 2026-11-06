# BP3 coarse consolidation

The BP3 driver can move a small distributed AMG level onto rank zero, continue
AMG there, and return the correction to its original owners. This removes
communication on levels with too few rows per GPU. It does not eliminate the
distributed work above that level.

The driver solves the same scalar diffusion problem and uses the same continuous
Qp spaces with GLL nodes and `p+2` Gauss points per axis, Dirichlet elimination,
p hierarchy, smoothing and FCG.
This is a steady operator; there is no time discretization. If `G` gathers the
uniquely owned rows at the selected level, the root matrix is `G A G^T` and the
returned correction is `G^T B_root G r`. The child correction starts from zero
on every application. Matrix entries and vectors stay on GPUs; only counts and
partition metadata are exchanged on the host. Persistent GPU-aware MPI requests
transfer vectors after their producing stream completes.

`MARS_ENABLE_BP3_COARSE_CONSOLIDATION` defaults to `AUTO`. A compatible Linux
CUDA build with MPI, Hypre, unstructured, FEM examples and MARSIR enables the
adapter. Unsupported configurations retain ordinary Hypre. `ON` requires support
and fails configuration if it is unavailable; `OFF` excludes the adapter.
Existing CMake caches retain their previous value: select `AUTO` explicitly when
upgrading a build previously configured with `OFF`.

In a supported build, `mars_marsir_bp --mg --amg` automatically uses the accepted
adaptive policy with the default AMG settings, GPU-aware MPI and Hypre's own
SpMV/SpGEMM. Changed hierarchy settings or vendor sparse kernels retain ordinary
AMG unless consolidation is explicitly requested. `--amg-consolidate 0` disables
it. A build without the adapter rejects a nonzero runtime request.

`-1` selects `min(50000, max(5000, 1000*ranks))` global rows. A positive value
from 10 through 50000 selects a fixed cutoff; `0` disables consolidation. One
rank uses ordinary AMG. If coarsening stalls above the cutoff, the distributed
coarse relaxation remains in use. The experimental two-stage `-2` mode is not
included.

Support is deliberately limited to the validated Linux CUDA build of shared
Hypre 2.33.0, real FP64, with 32-bit local indices. The adapter uses two private
Hypre callbacks, so CMake checks that the linked library exposes their dynamic
relocations. Compile-time layout checks and runtime AMG-recipe checks reject
unsupported combinations. In particular, the tested scalar recipe uses coarsen
type 8, interpolation 6, no aggressive levels, l1-Jacobi smoothing and four
coarsest sweeps. This is not a general Hypre extension or HIP support claim.

Only `mars_marsir_bp` links the adapter. Other applications using `P1Amg` keep
their ordinary Hypre calls. The driver reports host wall time after completing
each solve, taking the maximum over ranks. Setup and warm solve times remain
separate.

`tests/bp3_coarse/fixture.cu` checks gathered matrix action, repeated correction,
symmetry, positivity, reset/lifetime handling and manufactured solves, including
empty-rank partitions. `tests/bp3_coarse/integration/` builds the actual MARS
targets with the option OFF and AUTO, checks callback isolation and links, and
compares the integrated executable with the accepted CUDA executable at 1, 8
and 64 ranks. Integration is accepted only after native numerical checks and
no more than a 3% median warm-solve regression on those matched cases.

The earlier explicit-ON build passed this gate on GH200 with CUDA 12.9 and shared Hypre
2.33.0. Preparation 5269058 passed the CTest, 15 fixture invocations, 30 solver
checks on one and four GPUs, and three rejected CLI modes. Comparison 5269059
passed all 48 solver invocations and all eight timing comparisons on 1, 8 and
64 GPUs. Iterations and coarse hierarchies match the accepted executable; the
largest true residual is `8.008e-11` against the `1e-10` limit.

For 64 GPUs, median warm solve times in milliseconds were:

| Case | Accepted executable | Normal MARS build |
|---|---:|---:|
| Weak scaling, p=3, ne=256 | 258.655 | 260.335 |
| Weak scaling, p=7, ne=128 | 347.011 | 346.518 |
| Strong scaling, p=3, ne=128 | 97.084 | 97.385 |
| Strong scaling, p=7, ne=64 | 103.107 | 104.006 |

These are three alternating rounds, reporting the third solve in each run and
the maximum host wall time over ranks. Both executables use consolidation;
this verifies integration equivalence, not a new speedup against stock AMG.
The largest median slowdown across all eight cases is 0.87%. Setup is excluded.

Automatic selection passed its native gate on the same GH200 environment.
Preparation 5269652 built seven normal/private binaries, passed the production
CTest, 20 fixture invocations, 48 solver checks and three CLI rejection cases.
Comparison 5269653 passed all 72 solver invocations. Across all 120 solver
checks, the largest true residual was `8.017e-11` and the largest relative FCG
stopping measure was `8.773e-14`. Automatic selection matched accepted iteration
counts and hierarchies in all eight timing cases; its largest median regression
was 0.45%, within the 3% gate. This confirms selection equivalence, not a new
algorithmic speedup.

The same comparison rejected raising the root cutoff to 500000. The following
64-GPU medians compare the accepted adaptive cutoff and the larger fixed cutoff
on the same private executable, in three alternating rounds. Times are the
maximum host wall time over ranks for the third solve; setup is excluded.

| Case | Adaptive cutoff | Larger cutoff | Change |
|---|---:|---:|---:|
| Weak scaling, p=3, ne=256 | 262.843 ms | 283.463 ms | +7.8% |
| Weak scaling, p=7, ne=128 | 354.841 ms | 373.379 ms | +5.2% |
| Strong scaling, p=3, ne=128 | 96.298 ms | 117.391 ms | +21.9% |
| Strong scaling, p=7, ne=64 | 104.738 ms | 105.527 ms | +0.8% |

The first three cases moved an earlier level onto one GPU and remained correct,
with unchanged iteration counts, but every larger-root sample was slower than
every adaptive-control sample. The last case selected the same level under
both limits. Setup plus first solve was cheaper in median with the larger root;
that does not meet the repeated warm-solve acceptance gate. No larger cutoff
is promoted. Root memory snapshots were taken after setup and are not peak
allocation measurements.

`tests/performance/bp3_coarse_policy/` retains the reproducible comparison.
The production adaptive cutoff remains capped at 50000. Distributed work above
the selected level, DSS and global solver reductions remain; this result does
not establish perfect weak or strong scaling.
