# Prepared Hypre solves and direct AMG-cycle experiment

Author: GPT/Codex. Date: 2026-09-12.
Status: implementation; 46 local host FGMRES checks pass. User-executed Daint
GPU lifecycle gates pass on one and four ranks. Public flow and timing evidence
is recorded below; CUDA/Hypre execution was on Daint, not on the laptop.

The measured public profile spends about 24% of correction time on Hypre
preparation/setup and 73% on inner solves. This change exposes two opt-in paths
to test those costs separately. It implements no cross-step cache or SIMPLE
solver and does not change the beta=0.05 stability finding.

| `MARS_OUTLET_PRECOND` | Behavior |
|---|---|
| unset or `legacy` | Existing fresh Hypre GMRES solve for each outer preconditioner call |
| `reuse` | Same inner GMRES target and tolerances; reuse matrix, vectors and solver/AMG setup within a physical step |
| `amg-cycle` | Reuse BoomerAMG setup and apply one cycle directly for each outer preconditioner call |

Rank 0 selects the mode and broadcasts it once. Unknown values abort collectively.
The new paths use plain inner GMRES configuration; an inner FlexGMRES override
is rejected. Direct cycles require BoomerAMG, so a Jacobi override is rejected
for that mode. The existing outer FGMRES is unchanged.

## Numerical contract

The physical-pressure Chorin/BDF2 equations, velocity mask, Tet4 control-volume
quadrature, opening flux samples, and frozen trace/gradient remain unchanged.
With h=dt_eff/rho and S_ii=1/sqrt(V_i), the existing outer solve is

```
(S J) phi = -S R/h.
```

The legacy/reuse preconditioner approximately solves `Apre*z=S^-1*v` using
inner GMRES. The new cycle applies BoomerAMG once to the same RHS on the same
Apre. BoomerAMG retains the existing coefficient settings, zero tolerance and
one maximum iteration. Each application starts from zero, as before. A successful
cycle is an approximate preconditioner action, not a converged pressure solve.

FGMRES stores the actual preconditioned directions and verifies the residual
against the true J. The correction loop still checks every owned continuity row
and the stabilized boundary flux before accepting a step. None of those gates
were relaxed. Direct cycles reject nonfinite output or a globally zero action
on a nonzero RHS. They do not apply the inner-GMRES residual or amplitude-ratio
acceptance tests to a deliberately approximate action.

The reusable GMRES path retains the existing residual/null-solution acceptance
and `MARS_HYPRE_ACCEPT_RES` fallback. It bypasses the common wrapper's intermediate
node scatter: the preconditioner consumes owned DOF values, and the subsequent
true-J application calls `set_pressure` and publishes the needed node halo itself.
This avoids an unused scatter without changing the matrix action.

## Ownership and invalidation

`OutletKrylovOps` owns the prepared solver. It is created after the physical
step's matrix refresh and destroyed when that step's correction routine exits,
before MPI finalization. This bounds the cache to an interval where Apre, its
coefficient context and partition are frozen. The next physical step rebuilds;
BDF startup, changed timestep/coefficients and any between-step mesh changes
cannot inherit a stale hierarchy.

The wrapper also compares matrix/storage identity, row/column/NNZ dimensions,
global row/column ranges and device DOF-map storage. One rank's mismatch triggers
a rebuild on all ranks. A changed underlying ParVector after reinitialization
also triggers a collective rebuild. In-place matrix/map edits require collective
`invalidate_setup()`; pointer identity cannot detect them. The GPU gate exercises
that contract. Reuse does not cross a solver/preconditioner configuration change.

Persistent storage covers row indices, local/global mapping and cast buffers;
RHS and initial-guess values are refreshed on device each call. The Hypre IJ
vector lifecycle follows its documented initialize/set/assemble sequence.
See the [official IJ interface](https://hypre.readthedocs.io/en/latest/ch-ij.html)
and [BoomerAMG interface](https://hypre.readthedocs.io/en/latest/api-sol-parcsr.html).
No new matrix/field device-to-host copy is introduced in production; the existing
setup summaries remain. Metadata, API errors and finite-value decisions use
global scalar reductions. The reuse checks add collectives per application;
their scaling cost must be measured rather than assumed negligible.

Shared wrapper factoring moves vector updates and solve/extraction into common
methods. Legacy callers still rebuild their setup. Three existing matrix-building
kernel launches now skip zero owned rows; their surrounding collectives do not.
The public GPU gate performs small fixture uploads/downloads for verification.

## Validation and expected profiling

Completed locally: C++20 build with strict warnings and all 46 production FGMRES
controller checks, including varying preconditioners, nonfinite actions, failed
preconditioners and true-residual verification. Scoped diff checks pass. This does
not validate the modified Hypre wrapper or GPU execution.

New target `mars_outlet_hypre_gate` is built when Hypre is enabled. It uses a
procedural nonsymmetric diagonally dominant 512-row cyclic matrix, with genuine
off-rank couplings on multiple ranks and an analytic solution. It checks:

- fresh and reused GMRES solutions against each other and the analytic solution;
- new RHS data, including a sign/scale change for a cached AMG action;
- repeatability of an AMG action from zero and a useful fresh action on this fixture;
- one setup across repeated applications, then a second after explicit matrix invalidation;
- collective rebuilding when only rank 0 changes its map allocation.

Fresh AMG hierarchies need not be identical, so the gate does not demand bitwise
agreement between independently built cycles. The user supplied PASS transcripts
for this gate on one and four ranks. A separate two-rank run was not supplied.

With profiling enabled, legacy has one setup per preconditioner call. Reuse and
cycle modes should have one setup per physical step that actually calls the
preconditioner. Preparation, solve and finish still occur per call. In cycle mode,
`cg_p` counts AMG cycles; the mode banner identifies those units. Cached teardown
is included in total correction time, outside individual preconditioner calls.

## Daint commands after pulling and rebuilding

Run the GPU lifecycle gate before the flow comparison:

```bash
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_outlet_hypre_gate
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=4 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_outlet_hypre_gate
```

Then run the existing public fixed-trace case with reuse first:

```bash
set -o pipefail
env -u MARS_SOLVE_TRACE MARS_OUTLET_PROFILE=1 MARS_OUTLET_PRECOND=reuse \
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_pump \
  --mesh=../tests/data/public_outlet_channel/outlet_channel.exo --inlet-ss=inlet --outlet-ss=outlet \
  --solver=hypre --vms-stab --rc-implicit --outlet=do-nothing --outlet-beta=1 \
  --rho=1000 --nu=1e-4 --inlet-velocity=0.5 --dt=2e-6 --num-steps=200 --source-ramp-steps=100 \
  --relax-u=0.3 --relax-mass=1 --tol=1e-8 \
  --outlet-max-corrections=200 --outlet-rtol=1e-9 --outlet-div-tol=1e-8 --outlet-flux-tol=1e-10 \
  2>&1 | tee channel-C-reuse.log
```

For the separate cycle experiment, change only `MARS_OUTLET_PRECOND=amg-cycle`
and the output name to `channel-C-amg-cycle.log`. Keep beta=1 and all tolerances.
Confirm full continuity and boundary balance, and compare total work rather than
only inner iteration counts. A matched timing claim requires subsequent runs
with both profiling and solve tracing disabled.

## Daint validation and measured comparison

Recorded by GPT/Codex on 2026-09-12 from user-executed public-channel runs.
Both lifecycle gates passed, including collective invalidation. Both one-rank
flow modes completed 200 steps. The complete profiled logs were pulled with
authorized rsync and inspected locally:

| Log | SHA256 |
|---|---|
| `channel-C-reuse.log` | `0e0fa14bff07e02bb8958949cbec32aae616875fd35e399d07d85acb06185fd4` |
| `channel-C-amg-cycle.log` | `95ef2c8b3771f550c4deac701bd06da279b5114962bf3a941dc487050904b5b7` |

Both modes report one setup per physical step. Mean correction times over steps
3–200 are 214.10 ms (reuse), 29.38 ms (cycle), and 267.86 ms in the previously
recorded legacy profile. The cycle run uses 5,466 preconditioner applications
versus 3,571 for reuse: cheap approximate actions win despite more outer work.

The subsequent matched one-rank commands disable both `MARS_OUTLET_PROFILE` and
`MARS_SOLVE_TRACE`. The supplied complete transcripts give:

| Mode | 200-step runtime | ms/step | Final continuity RMS /s | Final signed boundary flux |
|---|---:|---:|---:|---:|
| legacy | 57,935.4 ms | 289.7 | 3.86e-13 | 4.11e-15 |
| amg-cycle | 8,000.4 ms | 40.0 | 5.33e-13 | 2.38e-14 |

This pair measures 7.24x lower total runtime (86.2% reduction) on this fixture.
It is not a repeated-run statistical benchmark. Both final reports have
`u_rms=0.6615`, `u_max=1.233`, and stabilized cut fluxes 0.500/0.500/0.500.
This is agreement at printed precision, not a field-norm comparison.

The user's four-rank AMG transcript covers steps 21–200 and reaches step 200
with the same printed velocities and cuts, RMS 6.10e-13 /s and signed boundary
flux -6.64e-14. Every supplied profile row has one setup on every rank. Mean
correction time over those steps is 85.98 ms versus 29.35 ms for one rank over
the same interval. This tiny fixture does not demonstrate useful strong scaling.
The four-rank excerpt does not contain the final runtime/completion footer.

These results support the opt-in cycle mode for further beta=1 work. They do
not validate beta=0.05, larger meshes, long-time physics, or a SIMPLE solver.
