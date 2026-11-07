# BP3 paper comparison with the current operators

This freezes a merged MARS revision, its matching generated operators, MFEM 4.8
(`a01719101027383954b69af1777dc828bf795d62`), the same CUDA image and Hypre 2.33
library. Fresh preparation builds both drivers in an isolated result directory; the current
resume reuses the verified binaries and library described below. Existing paper
builds and prior benchmark results are untouched. Native results are pending.

## Comparisons

Every case measures four explicit configurations:

- `matched`: MARS, h0, consolidation disabled, MFEM Chebyshev interval.
- `mfem`: MFEM CG, with the MARS Hypre settings, GPU-aware MPI and its own SpMV.
- `best_h0`: MARS, its own Chebyshev interval and adaptive coarse consolidation.
- `best_h`: the same, plus h2 at p3 or h1 at p7, explicitly at every rank count.

The first pair matches the preconditioner settings; MARS uses FCG and MFEM CG,
so do not describe their finite-precision algorithms as identical. The h variants
change the preconditioner. `best_h` names a candidate configuration, not an assertion
that it wins every case. Retain its slower single-GPU and cold-start samples. The
h0 control isolates the additional h hierarchy's effect. The automatic depth table
is not used, and no unmeasured depth is described as a general selection policy.

Each invocation runs five zero-initialized solves. Report the fifth maximum-rank
host wall time, three cyclically ordered rounds, with all three samples and their
minimum/median/maximum. Report setup and setup-plus-first-solve separately. True
residual must be at most 1e-10 and the relative stopping measure at most 1e-13, with
convergence before 2000 iterations. All failures remain in the report; after a
failed configuration its remaining rounds are skipped, not replaced. There is no
speedup acceptance threshold in this descriptive comparison.

Exact constant-element-work weak scaling uses 1/8/64 ranks: p3 ne64/128/256 and
p7 ne32/64/128, 16,777,216 element-local DoFs per GPU. Strong scaling uses
1/4/16/64/256 ranks at p3 ne128 and p7 ne64 (134,217,728 element-local DoFs globally).
Unique DoF counts differ by degree and are recorded separately. Report whole-solve
and per-iteration scaling separately, with the iteration counts. Never divide
historical MFEM times by current MARS times to claim a measured speedup.

The 256-rank job first checks all four configurations at ne64, p3/p7, plus two
instrumented direct hierarchy checks. It proceeds to strong-scaling timing only
if all pass. This is not a 256-GPU weak-scaling claim. The cube-only driver cannot
exactly preserve the existing per-GPU work at 4/16/256 ranks; rectangular weak
scaling would require a separate numerical change and validation in both drivers.

## Jobs

After the kit is transferred, on Daint:

```bash
bash /capstor/scratch/cscs/gandanie/bp3-paper-kit-20261107-v3/paper/submit.sh
```

This submits three jobs:

1. Two nodes, one-hour limit: resume preparation 5286400, reuse the MARS
   binaries and 24 completed checks, reuse the installed MFEM library, build only
   the MFEM driver, then run the 18 remaining baseline checks on 1/4/8 GPUs.
2. Sixteen nodes, four-hour limit, after preparation passes: 168 invocations
   covering the weak and strong cases above through 64 GPUs.
3. Sixty-four nodes, 90-minute limit, after preparation passes and the preceding
   job ends: ten pilot invocations, then 24 strong-scaling invocations at 256 GPUs.

The third job independently validates preparation and the pilot. A failed sample
in the standard comparison remains a failure but does not prevent this separate
pilot. No build occurs in the large allocations. Logs are `bp3_paper_*_<job>.log`;
results are `/capstor/scratch/cscs/gandanie/bp3-paper/<job>/`. Pass the printed job IDs
when returning the results. The user submits and manages every cluster job.

The kit manifest seals source, PTX, scripts and public dependency inputs. Build
records seal executables and runtime linkage; rank wrappers record executable
hashes, exit status, host placement and adapter audits. The parsers check numerical
receipts and timing fields. A final analysis must also verify the raw hashes and
commands, and distinguish pending, failed, and accepted measurements.

After retrieving the results, `analyze.py` checks the evidence seal and writes the
comparison tables. It refuses mixed manifests or duplicate cases. Efficiencies
use each configuration's own one-GPU baseline; explicitly label a curve assembled
from the fastest configuration at each scale if the paper also shows one. This
prevents a slower h-enabled one-GPU baseline from disguising a scaling cost.

## Four-rank MFEM failure handling

The old MFEM batch reported a shutdown assertion on four ranks. The exact raw
signature has not been verified here, so a printed PASS followed by a nonzero exit
is never promoted to success. Only the two four-rank MFEM preparation checks are
optional: any failure is retained in `optional-baseline-failures.json`, with raw
logs and process receipts, and has no accepted timing. One- and eight-rank MFEM
preparation checks must finish cleanly, as must every MARS check. This lets the
independent scaling work proceed without pretending that a crash was validated.
The standard timed runs and the 256-rank pilot still require clean exits.

The kits use separate immutable directories; they do not alter a previously
submitted job. The v3 command supersedes the earlier submission commands. Numerical source, stopping tolerances and timing acceptance rules are unchanged.
The MFEM consumer build is corrected below.

Each MARS paper-run invocation prints and stores `generated_sha256`: SHA-256 of the compact,
sorted JSON map from relative filenames under `generated/` to their SHA-256 values.
MFEM records null/unused because it does not load these operators. The manifest
continues to seal every individual PTX file.

## Resume after preparation 5286400

The current-operator MARS builds, 24 solves and eight direct hierarchy checks pass.
Their two executable hashes, 280 sealed artifacts, 97 process receipts and 97 rank
audits were independently verified. The largest true residual is 2.098e-11 and
the direct hierarchy error is 2.649236e-16. The MFEM library built and installed,
but its driver's Makefile invocation omitted the CUDA include directory and failed
on `cusparse.h`; no MFEM solve or scaling comparison ran.

The driver now uses an installed-package CMake consumer linked to CUDA targets,
which supply the missing public headers and runtime libraries. The source and
compiler remain the same; this is a build-recipe fix, not a numerical change.
A local compiled regression checks include propagation and reproduces the error
when the CUDA targets are removed. GPU compilation of the corrected driver is
still pending.

The v3 kit includes a receipt sealing the previous binaries, installation and raw
checks. Resume verifies it and unchanged production source, dependencies and PTX,
then copies the completed work into its new result directory. It does not rebuild
MARS or MFEM's library or rerun the completed 24 solves. The failed logs remain
unchanged. Only a passing resumed preparation enables the dependent scaling jobs.
