# Automatic BP3 consolidation and larger-root comparison

The accepted adaptive policy uses at most 50,000 global coarse rows. This kit
checks its automatic selection, then tests gathering an earlier AMG level
directly onto one GPU with a private fixed limit of 500,000 rows. It does not
repeat the staged 64-to-8-to-1 experiment. The production adaptive limit is
unchanged.

The measured remaining cost includes distributed AMG levels above the selected
root, p=1 DSS, and returning the correction. A larger root removes distributed
levels but increases work on one GPU; the comparison measures that tradeoff. BP3's PDE, spaces,
quadrature, fine operator, p hierarchy and stopping tests stay unchanged. At
the selected level the root matrix is exactly the owned-row permutation
`G A G^T`; the preconditioner may change because AMG is rebuilt there. Candidate
iteration counts are therefore allowed to change. Every solve must meet the
same accuracy gates.

`export_policy.py` freezes the current code, pinned public dependencies, accepted
PTX and private candidate sources. `policy.py` reuses the already validated
integration launcher, MPI descriptor handling, audit export and log parser.
The CTest coordinator runs in the batch process and launches a fresh MPI step.
No external Hypre library is modified. The image and Hypre library/header
hashes remain those of the accepted integration.

Preparation builds normal OFF/AUTO targets, the other P1Amg consumer and the
private pair of solver/fixture targets. It runs the production CTest, 20 fixture
invocations, 48 solver checks at p=1/3/7 on one/four GPUs, and three invalid CLI
cases. The large fixtures include empty root and empty tail partitions and
must actually select more than 50,000 rows. It checks automatic selection
against explicit `-1`, and private `-1` against the production default. Only a
fully sealed passing preparation enables the timing job.

The comparison makes 72 solver invocations. Eight weak/strong cases at 1/8/64
GPUs compare the accepted executable with the normal executable with no
consolidation flag. The four 64-GPU cases additionally compare private `-1`
and private `500000`, on the same binary. Three alternating rounds report the
third solve's maximum host time over ranks. All timings, first solves, setup
costs and failures are retained. Root memory is a device-memory snapshot after
AMG setup, not a peak allocation measurement.

Acceptance is fixed before the run:

- Every solver invocation must pass true residual `<=1e-10` and FCG relative
  stopping measure `<=1e-13`, plus all-rank lifecycle and ownership checks.
- Automatic selection must preserve iterations and hierarchy, with median
  warm-solve regression no greater than 3% against the accepted executable.
- The private `-1` control must preserve iterations/hierarchy and stay within
  3% of the normal automatic build.
- When the larger limit selects an earlier level, it needs at least 5% median
  whole-solve improvement over its private control, with all three samples
  faster than every private-control and normal-default sample. An unchanged
  level permits at most 3% regression.
- The script never promotes the larger limit. Promotion also requires review
  of setup plus first-solve cost and memory headroom. In particular, a cold-cost
  increase must be justified by the measured number of solves needed to amortize
  it. A timing rejection does not invalidate a passing automatic-selection gate.

Local CPU/MPI fixtures and orchestration tests do not establish GPU performance.
After transferring and verifying the frozen bundle, the user runs on Daint:

```bash
bash /capstor/scratch/cscs/gandanie/bp3-coarse-policy-kit/policy/submit.sh
```

This submits one node for preparation (three-hour limit), followed only on
success by sixteen nodes for the comparison (90-minute limit). Results go to
`/capstor/scratch/cscs/gandanie/bp3-coarse-policy/<job-id>/`.

## Completed GH200 comparison

Preparation 5269652 and comparison 5269653 passed all numerical checks. Automatic
selection passed all eight equivalence gates, with at most 0.45% median regression
against the accepted executable. The private control also passed equivalence.

The larger root was rejected: in the three cases selecting an earlier level,
median warm solves were 7.8%, 5.2% and 21.9% slower, with separated ranges and
unchanged iteration counts. The strong p=7 case selected the same level under
both limits and stayed within the 3% regression gate. Setup plus first solve
was cheaper in median, but the measured objective was repeated warm solves.
The production cutoff remains unchanged. Details and the timing table are in
`docs/design/bp3_coarse_consolidation.md`; a rerun is not needed to accept AUTO
or reject this larger-cutoff proposal.
