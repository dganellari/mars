# Opt-in outlet and Hypre phase timings

Author: GPT/Codex. Date: 2026-09-12.
Status: implemented; profiler header compiled and tested with real local MPI and
a CUDA synchronization stub. Full CUDA/Hypre compilation and GPU timing pending.

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

## Daint public run after pulling and rebuilding

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
