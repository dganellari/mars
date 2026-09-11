# Opt-in outlet and Hypre phase timings

Author: GPT/Codex. Date: 2026-09-12.
Status: implemented and executed in a 200-step public one-rank GPU run; local
profiler tests also pass on 1/2/4 real MPI ranks with CUDA synchronization stubbed.
Multi-rank GPU profiling remains unverified.

`MARS_OUTLET_PROFILE=1` enables one `[outlet-profile]` line per physical outlet
correction step. This measures the repeated setup identified in the
[cost review](gpt_outlet_setup_cost_review_2026-09-12.md). No matrix, solver
tolerance, pressure update, or boundary condition changes.

## Measurement boundaries

| Field prefix | Timed work |
|---|---|
| context | Build/publish the frozen VMS context |
| assembly | Restore K, assemble compact derivatives, anchor checks, publish Apre |
| hypre_prepare | Hypre initialization, old-state cleanup, DOF map copy, matrix/vector creation and data preparation |
| hypre_setup | Preconditioner/Krylov creation and configuration, AMG/Krylov setup, existing setup barrier |
| hypre_solve | Existing solve barrier and Hypre GMRES/FlexGMRES solve |
| hypre_finish | Solution extraction, residual/finite-value checks and existing reporting |
| preconditioner | Entire outer-FGMRES preconditioner call, including scaling, all Hypre phases, destruction, scatter/halo and output copy |
| correction | Entire pressure-correction routine, including validation and final diagnostics; excludes momentum predictor/diffusion and the timing report itself |

Every prefix has `_ms`, `_calls_min`, and `_calls_max` fields. Durations are
accumulated locally across all corrections in a step, then reduced with MPI_MAX.
Call counts report rank minima and maxima. Normally all four Hypre phase counts
equal the preconditioner count; assembly/context/correction each occur once.
The profiler owns its step counter: the solver's BDF startup flag saturates at 1
and cannot label physical steps.

Preconditioner and correction times include their child phases. Do not add them
to the children. Different phase maxima can come from different ranks, so their
sum is not a rank-consistent breakdown. Small gaps for timing boundaries and
setup-error reporting are included in the parent, not assigned to Hypre phases.

Only enabled timings perform checked GPU synchronizations at phase boundaries.
The three reporting reductions occur once after the correction loop; the new
timers add no inner-loop MPI collectives. Existing Hypre barriers remain. Profiling
serializes GPU work and adds overhead: use it to locate cost, not as an uninstrumented
performance benchmark. Reports contain durations and call counts, no mesh metadata.

Rank 0's flag is broadcast once per stepper when the average-pressure correction
is first entered. This prevents different environments from selecting different
collective paths. With profiling disabled, this single broadcast remains; there
are no new timing synchronizations, per-step reporting reductions or output.
The legacy outlet path does not enter this initialization. Other Hypre callers
retain a null profiling pointer.

## Validation completed locally

Compiled the production profiler header using C++20 and `-Wall -Wextra -Werror`.
Executed the test below on one, two and four real MPI ranks: all passed. It checks
deliberately conflicting rank-local environment flags, no synchronization/output
when disabled, rank-zero-only reporting, count minima/maxima, and step reset.
Only CUDA synchronization is stubbed; these checks do not compile the surrounding
CUDA kernels or validate Hypre on a GPU.

```bash
mpicxx -std=c++20 -Wall -Wextra -Werror -Itests/outlet_profile_stubs -I. \
  tests/outlet_profile_host_check.cpp -o /tmp/mars-outlet-profile-check
mpirun -np 1 /tmp/mars-outlet-profile-check
mpirun -np 2 /tmp/mars-outlet-profile-check
mpirun -np 4 /tmp/mars-outlet-profile-check
```

The stub include directory is test-only and is not registered with any production
target. The production changes introduce no new device lambdas or kernels.

## Executed public GPU result

GPT retrieved the complete `channel-C-profile.log` by authorized rsync and
inspected it locally. The header confirms the public fixture, one rank, beta=1,
rho=1000, nu=1e-4, dt=2e-6 and the expected relaxation settings. There are 200
consecutive profile records, no solve-trace records and no failure/nonfinite
markers. Every phase has matching rank-min/max call counts, and all four Hypre
phase counts equal the preconditioner count. There are 3569 preconditioner calls
over the complete run. The header does not embed an executable revision.

Means over steps 3-200, excluding BDF startup:

| Phase | ms/step | Share of correction time |
|---|---:|---:|
| VMS context | 0.043 | 0.016% |
| Outlet matrix assembly | 0.101 | 0.038% |
| Hypre preparation | 16.404 | 6.124% |
| Hypre setup | 48.072 | 17.947% |
| Hypre solve | 195.315 | 72.917% |
| Hypre finish | 1.695 | 0.633% |
| Whole preconditioner, including the Hypre phases | 262.871 | 98.137% |
| Whole correction | 267.860 | 100% |

This is a single-rank run, so phase shares have no cross-rank-max ambiguity.
The whole-preconditioner and whole-correction rows are inclusive totals.
Steps 101-200 give essentially the same shares: solve 72.92%, preparation/setup
24.09%, assembly 0.036%. Step 1 is exceptional: correction 759.369 ms, including
394.525 ms of setup. It should not represent steady per-step cost.

The run finishes in 55628.7 ms (278.1 ms/step). Final full continuity RMS is
3.86e-13 /s, maximum 2.43e-12 /s and net boundary imbalance 4e-15 volume/s.
The final maximum speed is 1.233 and all three stabilized cuts print 0.5000.
These agree with the previous fixed-trace results at reported precision; they
do not constitute a field-level or multi-rank GPU comparison. Profiling is on
and solve tracing is off, so do not interpret its wall time as a controlled
speedup against the earlier trace-enabled baseline.

```
remote: /capstor/scratch/cscs/gandanie/git/mars/daint-gpu/channel-C-profile.log
local:  /private/tmp/mars-outlet-audit-20260911/channel-C-profile.log
SHA256: 9765b1fe190223e63cfe5908fe4e497db98762967a5522d15f1b0c7ec733c2d5
```

The measured ranking settles the earlier hypotheses. Matrix assembly is not a
meaningful performance target here. Persistent Hypre preparation/setup has a
useful but limited opportunity: those phases account for about 24% of correction
time, and preparation includes per-RHS work that cannot all disappear. This is
an opportunity estimate, not a promised speedup. The dominant cost is the nested
Hypre solve, which must also be addressed for a large improvement.

Recommended next scope: High effort, one agent, a prepared Hypre/AMG context with
explicit invalidation and a separately selectable direct AMG-cycle preconditioner
experiment for the existing outer FGMRES. First verify reusable solves preserve
the current equations and residual gates; then measure whether avoiding the
inner GMRES reduces total work without excessive outer iterations. Keep the
true-J solve and full continuity/flux acceptance unchanged. Neither a single
AMG cycle nor persistence has an established speedup yet; no cache or algorithm
change is part of this profiling commit.

## Reproduce the public GPU run

From the existing Daint build directory, use the fixed-trace public case. This
keeps the stricter original tolerance and disables solve tracing so per-solve
printing does not dominate the profile. Its wall time is therefore not directly
comparable to the earlier trace-enabled runs.

```bash
set -o pipefail
env -u MARS_SOLVE_TRACE MARS_OUTLET_PROFILE=1 \
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_pump \
  --mesh=../tests/data/public_outlet_channel/outlet_channel.exo --inlet-ss=inlet --outlet-ss=outlet \
  --solver=hypre --vms-stab --rc-implicit --outlet=do-nothing --outlet-beta=1 \
  --rho=1000 --nu=1e-4 --inlet-velocity=0.5 --dt=2e-6 --num-steps=200 --source-ramp-steps=100 \
  --relax-u=0.3 --relax-mass=1 --tol=1e-8 \
  --outlet-max-corrections=200 --outlet-rtol=1e-9 --outlet-div-tol=1e-8 --outlet-flux-tol=1e-10 \
  2>&1 | tee channel-C-profile.log
```

Confirm successful completion and the existing full continuity/flux reports,
then inspect the phase sums and call counts. Treat startup separately from later
BDF2 steps. Actual speedup from a future cache requires a matched run with profiling
and tracing both disabled; this commit implements no cache.
