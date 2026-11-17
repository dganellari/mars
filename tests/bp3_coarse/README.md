# BP3 coarse consolidation checks

The production adapter is in
`backend/distributed/unstructured/solvers/mars_cellwise_amg_consolidate.cu`.
The fixture includes that implementation so it can check its gathered matrix and
lifetime state directly. It is not linked into the MARS library.

With `MARS_ENABLE_BP3_COARSE_CONSOLIDATION=ON` and tests enabled, build
`mars_bp3_coarse_test`. The `marsBP3CoarseSingle` CTest launches its single-rank
case through CMake's `MPIEXEC_EXECUTABLE` and `MPIEXEC_PREFLAGS`. Run CTest from
the batch process, so its MPI test can start its own scheduler step.
On an allocated CUDA node with GPU-aware MPI, also run the executable on four
ranks with `balanced`, `empty-root` and `empty-tail`, followed by the row count
511. Set `BP3_CONSOLIDATE_ROWS` to 0, 100 and 50000 for each partition. A balanced
32767-row case with cutoffs 0, 5000 and 50000 checks the larger coarse limit.
These 15 cases perform 60 manufactured solves, with action, symmetry,
positivity, repeated-application and lifetime checks. Their true residual must
be below `1e-10`.

## Whole-bottom consolidation

`--amg-consolidate root` gathers the complete assembled bottom matrix onto
rank 0 before AMG setup. Each correction gathers the RHS once, runs the
configured number of AMG cycles on that GPU, and scatters the result once.
The p and h transfers and the assembled matrix are unchanged; the AMG hierarchy
can differ from distributed setup. This is an explicit best-configuration
option, not the matched MFEM configuration or the automatic default.

The mode refuses matrices above 300000 rows or 30000000 entries before
allocating the root matrix. These bounds do not include AMG hierarchy storage.
On one rank it uses ordinary AMG. It has been checked on CUDA with shared
Hypre 2.33, including empty owners; this does not establish HIP support or
arbitrary-scale weak scaling.

Set `BP3_CONSOLIDATE_ROWS=-2` for the fixture's whole-bottom mode. Run balanced
511-row cases on 1, 4 and 8 ranks, empty-root and empty-tail cases on 4 ranks,
and a balanced 32767-row case on 4 ranks. Besides the checks above, the fixture
requires three cycles per root exchange and checks changed vector bindings.
Keep the disabled and nested-consolidation controls as well.

The site validation under `integration/` freezes a code-only source snapshot,
the pinned Cornerstone/GoogleTest sources and the exact accepted PTX. It builds
the normal top-level MARS project twice, with consolidation OFF and ON. Both
builds also compile the other cell-wise solver, and verify that it exports no
adapter callbacks. Fresh builds use separate directories.

The preparation job checks runtime Hypre linkage, callback scope, the CTest,
all 15 fixtures, 30 BP3 smoke invocations (degrees 1, 3 and 7; one and four
ranks), and three CLI rejection cases. Stock and consolidated runs compare
against the corresponding mode of the accepted executable. A dependent job
runs 48 alternating accepted/integrated invocations: degrees 3 and 7, weak
scaling on 1/8/64 GPUs, plus strong scaling on 64 GPUs. It requires matching
iterations and coarse hierarchies, true residual at most `1e-10`, and at most
a 3% median warm host-solve regression. Every timing is retained.

The exact accepted control is the user-run preparation 5253439, built with
CUDA 12.9 and shared Hypre 2.33.0. The site paths and hashes are in
`integration/common.py`. This gate checks integration equivalence; it does not
establish a new algorithmic speedup, cold-start gain or support for another
Hypre version. Optional audits go to node-local memory and are copied out only
after the solver exits.

Run local orchestration checks with
`python3 -B tests/bp3_coarse/integration/test_integration.py`. They exercise the
actual extracted C++ CLI, CMake rejection guards, recorded native output parsing,
audit exports and relocated Slurm wrappers. They do not compile or run CUDA.
The historical result files used by the output-reader check are local evidence
under the preceding architecture worktree.

After the immutable bundle has been transferred to Daint, the user submits both
jobs with:

```bash
bash /capstor/scratch/cscs/gandanie/bp3-coarse-integration-retry-20261105/kit/submit.sh
```

The wrapper selects the same CUDA image for the batch process and each MPI
step. A failed preparation prevents the 16-node comparison from starting.

For the reviewed native build 5266939, pass its result directory as the sole
argument to `submit.sh`. This copies the five hash-checked executables into new
results without recompiling or modifying the old build. Source checks permit
only the CTest registration change and validation-script changes. The resumed
CTest uses the same registration file as the normal build. All numerical and
timing gates still run; no earlier numerical pass is assumed.
