# BP3 coarse-level profile on GH200

This experiment identifies which BoomerAMG levels account for the p=1 path's
increase from one to eight GPUs. It keeps the numerical solver, AMG settings,
matrix and installed Hypre library unchanged. It is a diagnostic, not a proposed
speedup or evidence that geometric h-coarsening will win.

`hypre_ranges.cpp` replaces Hypre 2.33's disabled profiling functions at dynamic
link time. It emits NVTX ranges from the existing `AMGCycle`, `AMG Level-N`,
relaxation and transfer markers. It adds no GPU wait, MPI operation or numerical
work. Its per-process audit rejects missing or unbalanced hooks. The build
preflight rejects static Hypre, unsupported versions, internal calls that cannot
be interposed, and native Hypre profiling builds. No dependency rebuild or MARS
source patch is required; the job rebuilds the existing MARS target once.

Hypre's native NVTX implementation in 2.33 calls `cudaStreamSynchronize(0)` in
`hypre_GpuProfilingPopRange`. Enabling it would change the scheduling under study.
The replacement emits annotations only. This relies on the installed library
retaining externally bound profiling calls; the preflight and small smoke trace
must establish that on Daint before the large runs proceed.

## Run

Use the existing CUDA 12.9 uenv with Nsight Systems. The MARS checkout must include
`339e6552` (the true-residual acceptance gate). The job rebuilds
`mars_marsir_bp` in `/capstor/scratch/cscs/gandanie/git/mars/mlir`; let any other
build there finish first. It does not use `hip1/mars`. GPU selection remains in
the procedural BP driver's existing initialization.

Laptop, from the MARS repository root:

```bash
rsync -av tests/performance/bp3_coarse_profile/ \
  daint-alps3:/capstor/scratch/cscs/gandanie/bp3-coarse-profile-kit/
```

Daint, inside that uenv (the passthrough option keeps the active uenv in the job):

```bash
sbatch --uenv-passthrough=use -A csstaff -N2 -n8 --ntasks-per-node=4 --export=ALL \
  --exclude=nid005323,nid005617,nid005669,nid005679,nid005742,nid005820,nid005854 \
  --output=/capstor/scratch/cscs/gandanie/bp3_coarse_%j.log \
  /capstor/scratch/cscs/gandanie/bp3-coarse-profile-kit/job.sh pilot
```

The script has a one-hour limit. This is an allocation limit, not a measured
runtime estimate. It first runs a small stock/trace pair, then the following:

| GPUs | p=3 global ne | p=7 global ne |
|---:|---:|---:|
| 1 | 64 | 32 |
| 8 | 128 | 64 |

Each case runs three alternating stock/probe pairs without the profiler, then
one trace on every rank. Each invocation uses `--solves 3`: two warm-ups and one
reported solve. The trace captures only the last solve. Settings are explicit:
deformation 0.05, discrete RHS, MARS Chebyshev interval, (2,2) smoothing, three
AMG cycles, FCG tolerance 1e-13 and true residual <=1e-10. The tighter FCG
tolerance avoids accepting historical results just above the true-residual
threshold. Failed accuracy, changed iteration counts or changed printed history
stop the experiment. A trace is never used as an unprofiled timing result.

Output is `/capstor/scratch/cscs/gandanie/bp3-coarse-profile/JOBID/`:

- `job.exit`, `job.log`, `overview.txt`: completion and compact comparison.
- `build/manifest.json`: exact binary, Hypre library, config and hook hashes;
  source revision/status and PTX hashes are recorded alongside it.
- Each case's `solver.log`, `summary.json`, per-rank audit, Nsight report and CSV
  exports. Zero markers fail; missing MPI data is allowed only on a single rank
  with an explicit no-MPI-events report.

If this identifies the relevant levels and probe overhead is acceptable, the
same experiment can confirm the diagnosis at 64 ranks, with p3/ne256 and
p7/ne128. This separate job is not part of the pilot:

```bash
sbatch --uenv-passthrough=use -A csstaff -N16 -n64 --ntasks-per-node=4 --export=ALL \
  --exclude=nid005323,nid005617,nid005669,nid005679,nid005742,nid005820,nid005854 \
  --output=/capstor/scratch/cscs/gandanie/bp3_coarse_%j.log \
  /capstor/scratch/cscs/gandanie/bp3-coarse-profile-kit/job.sh scale
```

`MARS_BUILD`, `BP3_PROFILE_KIT`, and `BP3_PROFILE_RESULTS` can override the paths.
If automatic header discovery fails, `NVTX_INCLUDE` must name the directory
containing `nvtx3/nvToolsExt.h`. Do not substitute another Hypre build to make the
preflight pass; record the failure and adapt the probe to the actual library.

## Interpretation and next implementation

`overview.txt` reports the three unprofiled stock/probe times and the median
paired overhead ratio. For each AMG level it reports median and maximum across
ranks of host range duration, correlated kernel duration and time in MPI waits,
normalized by captured AMG cycle count. A p=1 adapter call contains three of
these cycles. The JSON retains per-rank values and projected GPU spans.

CPU, GPU and MPI measurements overlap. Do not add them. GPU projection spans
include gaps and can overlap between levels; kernel totals exclude copy time.
The analyzer counts level roots, excluding duplicate child summaries. MPI wait
time can include remote rank imbalance; it is not proof of network latency.
Hypre labels interpolation with the departing coarse level, although the work
writes the next finer vector. Inspect that operation before assigning its cost
to tiny-level work. The trace omits the adapter's outer DSS from AMG ranges;
the driver coarse timer includes them.

- If added cost is in the large first AMG levels, with appreciable sparse GPU
  work, test one geometric h-level below p=1 with a Galerkin coarse matrix.
- If it concentrates in small levels and MPI/launch waits, investigate late
  rank reduction or a different small-level solve before adding h smoothing.
- If ranks differ sharply, inspect ownership and per-rank timelines first.
- If probe overhead is substantial or attribution is incomplete, refine the
  measurement before deciding on an algorithm. No gain is inferred from markers.

These are decision rules, not findings from GPU execution.

## Local validation

```bash
python3 -B tests/performance/bp3_coarse_profile/test_profile.py
bash -n tests/performance/bp3_coarse_profile/job.sh tests/performance/bp3_coarse_profile/rank.sh
```

Five parser/launcher tests pass, covering nested-range double counting, duplicate MPI
send/receive rows, single-rank absent MPI data, missing hooks, unbalanced hooks,
failed convergence, nonfinite residuals, residuals above the target, and retaining
Nsight's own preload when attaching the probe. The C++
probe compiled with NVIDIA NVTX 3.1.1 headers and the real Hypre 2.33 CPU library.
`host_fixture.cpp` produced bit-identical 512-entry vectors in stock/probe builds;
nine cycles visited four levels, with 747 balanced pushes/pops and no audit errors.
That check used static linking to substitute the hook on macOS. It does not
validate Linux preload behavior, CUDA execution or Nsight export schemas on Daint.
Those remain explicit gates in the job. No cluster run is claimed.

## Source references

- [Hypre 2.33 cycle](https://github.com/hypre-space/hypre/blob/v2.33.0/src/parcsr_ls/par_cycle.c)
  and [profiling hooks](https://github.com/hypre-space/hypre/blob/v2.33.0/src/utilities/nvtx.c).
- [Nsight Systems analysis guide](https://docs.nvidia.com/nsight-systems/AnalysisGuide/index.html):
  NVTX GPU attribution uses the enclosing launch API range, not temporal
  containment of GPU execution. Parent and child reports overlap.

Suggested commit message for the MARS `cstone` diagnostic files, after review:
“Add a BP3 profiling experiment that attributes AMG work and MPI waits to each
level. Check solver agreement and probe overhead before using the traces to
choose a coarse-solve change.”
