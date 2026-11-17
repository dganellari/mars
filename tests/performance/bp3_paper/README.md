# Final BP3 paper comparison

This freezes MARS `ea1c0ef923246f7489ece4b847a11e43f78ebe4c`, all its tracked
generated operators, MFEM 4.8
(`a01719101027383954b69af1777dc828bf795d62`), the same CUDA image and Hypre 2.33
library. Preparation builds MARS and MFEM from scratch in a new result directory.
It also builds `mars_marsir_cellwise_solve` for the separate CVFEM campaign.
No historical executable, build directory or result is reused. Native results
for this final revision are pending.

The source checkout must be clean and at that exact revision. `export_paper.py
--repo SOURCE --reference V3_BUNDLE --target NEW_BUNDLE` verifies the sealed v3
reference before reusing its public dependency sources. Every generated file
comes from the clean source checkout, including CVFEM and DSS support kernels.
The manifest records production source and generated revisions separately from
the reviewed runner files, which may contain uncommitted harness changes.

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
bash /capstor/scratch/cscs/gandanie/bp3-paper-kit-20261109-final/paper/submit.sh
```

This submits three jobs:

1. Two nodes, five-hour limit: fresh MARS and MFEM builds, 24 hierarchy solves
   including eight direct checks, then 18 baseline checks on 1/4/8 GPUs.
   The MARS build step permits 125 minutes and MFEM 140 minutes; the remaining
   allocation time is for the runtime probes and small validation solves.
2. Sixteen nodes, four-hour limit, after preparation passes: 168 invocations
   covering the weak and strong cases above through 64 GPUs.
3. Sixty-four nodes, 90-minute limit, after preparation passes and the preceding
   job ends: ten pilot invocations, then 24 strong-scaling invocations at 256 GPUs.

The third job independently validates preparation and the pilot. A failed sample
in the standard comparison remains a failure but does not prevent this separate
pilot. No build occurs in the large allocations. Logs are `bp3_paper_*_<job>.log`;
results are `/capstor/scratch/cscs/gandanie/bp3-paper/<job>/`. Pass the printed job IDs
when returning the results. The user submits and manages every cluster job.

The submission checks the bundle and pinned image's presence before submitting
anything. The kit manifest seals source, PTX, scripts and public dependency inputs. Build
records seal executables and runtime linkage; rank wrappers record executable
hashes, exit status, host placement and adapter audits. The parsers check numerical
receipts and timing fields. A final analysis must also verify the raw hashes and
commands, and distinguish pending, failed, and accepted measurements.

Before builds or measurements, one probe on every allocated node records all four
GPU UUIDs/models, the driver, compiler versions, image path and Hypre hashes.
Preparation and both comparison allocations must have the same driver/toolchain
signature. A mismatch stops that phase before timing. `runtime.json` and the raw
per-node receipts remain in the evidence seal; these are not profiler runs.

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
submitted job. This final command supersedes v3 for new measurements. The frozen
operators are newer; workload definitions, tolerances and acceptance rules stay
the same. Keep v3's accepted results and failed samples under their old manifest.

Each MARS paper-run invocation prints and stores `generated_sha256`: SHA-256 of the compact,
sorted JSON map from relative filenames under `generated/` to their SHA-256 values.
MFEM records null/unused because it does not load these operators. The manifest
continues to seal every individual PTX file.

## Shared build for CVFEM and BP1

`paper-build.json` seals the normal BP driver, diagnostic BP driver, CVFEM driver
and MFEM driver. The normal MARS executables are under
`build/on/examples/distributed/unstructured/` in the preparation result directory;
generated kernels are under `generated/` in the bundle. No later job should
reconfigure or rebuild that directory.

The separate CVFEM/BP1 runner must call `driver.verify_prepared(prepared, bundle)`
before using these inputs and retain the manifest, binary and PTX hashes. It can
call `runtime.capture(bundle, results, node_count, prepared)` to require the same
driver/toolchain, with a new, empty result directory. Keep its outcomes separate:
building CVFEM does not validate a CVFEM solve, and the BP3 jobs do not exercise
`hl_full_p7_sm90.ptx`. CVFEM solve checks, wing scaling and BP1 Table 7 remain in
their own campaign.
