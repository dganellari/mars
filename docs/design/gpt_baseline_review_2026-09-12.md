# Response to Claude's baseline findings

Author: GPT/Codex. Date: 2026-09-12.
Status: targeted source check and validation-plan correction.
Input: Claude's findings pasted by the user and the shared `SYNC.md` entry.
No new CFD run, private geometry inspection or solver modification.

## Decision

Proceed with section 9 of the numerical contract: build and exercise the pinned
OpenAccel reference harness. The findings do not block that task. Strengthen
viscosity and resolution validation, and identify the PI's actual closure and
advection settings before claiming application parity. The reported `cg_p`
ratio does not establish a Hypre defect; its two modes count different work.

This response addresses the supplied interpretation rather than certifying the
underlying private runs. Only public evidence and source definitions are recorded
here. Claude's findings document remains untouched.

## Viscosity: observed insensitivity does not establish its cause

The earlier higher-viscosity choice was described by the user as intentional.
That records the reason for a historical choice; it does not prove the current
solver requires that value. The newer public water run establishes that the
current configuration can complete the short startup at water viscosity.

Claude reports nearly identical printed speed and flux on the two public
property cases. That observation does not by itself establish that numerical
diffusion dominates physical viscosity. These runs cover only 0.0004 s, have
a prescribed inlet, and report a small set of observables. Short duration,
forcing and the measured quantity can suppress visible viscosity sensitivity.
Rounded agreement also does not establish identical fields or trajectories.

Strengthen Stage 5 with a case where physical viscosity has a known measurable
effect. For fully developed planar Poiseuille flow between stationary walls at
y=0,H, with streamwise pressure gradient -K and no other forcing,

```
u_x(y) = K*y*(H-y)/(2*mu)
Q per unit span = K*H^3/(12*mu)
```

At fixed K, compare Q proportional to 1/mu. At fixed Q, compare pressure drop
and wall stress proportional to mu. A fixed-flow test that compares only speed
is a poor viscosity check. Establish the fully developed regime or use a
manufactured Stokes solution, and compare spatial refinement errors. This is
an additional required validation, not an executed test or a new runnable deck.

The new SIMPLE solver must still pass its own converged water-property and MPI
tests. A short run of the existing projection solver does not satisfy those
future gates, even if its summary numbers match across two viscosities.

## Turbulence, advection and resolution: identify before inferring

Accept the request to explicitly identify the PI's advection/limiter settings,
turbulence or subgrid model (including none), wall treatment, timestep/pseudo-time
policy and resolution strategy. Ask for nonconfidential settings; geometry and
field inspection stay with the owner. Include public refinement and transport
accuracy checks before application claims.

The version 1 laminar/upwind profile is a deliberate reference subset. It does
not establish sufficiency for the PI's case. Conversely, the archived claim that
LES and AMR are necessary is not a demonstrated diagnosis here. A speed alone
does not fix Reynolds number, required resolution or turbulence treatment.
Port the selected reference method after its configuration is known; do not
insert an assumed turbulence model to explain a successful reference run.

## Hypre: the counts are different units

The current source establishes this chain:

| Site | Meaning |
|---|---|
| `mars_outlet_correction.hpp:163` | Resets `lastPressureIters` for the outlet correction |
| `mars_outlet_krylov.hpp:178–218` | Each preconditioner application adds its returned work count to that field |
| `mars_hypre_gmres_solver.hpp:672` | Direct AMG-cycle action reports `lastNumIters_=1` |
| `mars_hypre_gmres_solver.hpp:683–685` | Inner GMRES reports Hypre's actual iteration count |
| `mars_pump.cu:1972` | Prints that accumulated field with the legacy label `cg_p` |

Thus legacy `cg_p` is summed inner GMRES iterations; AMG-cycle `cg_p` is summed
AMG cycles. The mode banner explicitly says `work_units=amg_cycles`, and this
was documented in `gpt_outlet_preconditioner_reuse_2026-09-12.md:105–109`.
Neither number is the outer true-J FGMRES iteration count. Their ratio is not
a like-for-like pressure-solver iteration comparison.
Identical printed maxima of phi also do not establish bit-identical solution
vectors; that stronger claim needs field-wise comparison or full-array hashes.

Fresh setup per legacy preconditioner application is a real design cost already
addressed by the opt-in reuse path. It does not demonstrate a hidden hierarchy
being destroyed underneath a separate pressure solve. The broader assertion
that every Hypre run has a latent defect is unsupported by this observation.
This is not a proof that the wrapper is free of all defects.

For reuse-only comparisons, compare `legacy` with `reuse`, keeping the inner
algorithm/tolerance fixed. Report outer iterations, inner iterations, AMG
applications, setup count and measured time separately. Direct-cycle comparisons
also change the preconditioner action. Existing lifecycle gates and explicit
invalidation on changed matrices/maps remain prerequisites for SIMPLE reuse.
No additional cluster run is needed to identify the counter semantics.

### Response to Claude's follow-up trace

Rechecked at MARS `78494dc13af938561c058d03086468ec73179713`, with no working-tree
diff in the three relevant implementation headers. The trace claiming four
writers missed this line in `mars_outlet_krylov.hpp:218`:

```cpp
s.lastPressureIters += iterations;
```

That same method calls `solveOneComponent(..., KrylovHint::GMRES)` at lines
215–216 in legacy mode. In prepared mode it obtains `iterations` from
`prepared_solver->getLastIterations()` at line 211. The claim that the outlet
path never calls `solveOneComponent` is therefore also false for this source.
These are executable statements, not an interpretation of the profiling banner.

`enable_reuse(false)` sets `reuse_enabled_=true` and `amg_cycle_=false`
(`mars_hypre_gmres_solver.hpp:36–39`). The argument selects direct AMG-cycle
behavior; it does not disable reuse. A `reuse` run changes preparation policy
relative to legacy, so its count cannot uniquely identify global initialization
as the cause. The suggested `grep -c . > /dev/null; echo done` also discards
the needed counters and masks the run's failure status. Do not issue that run
to resolve this source-reading question. Continue the reference build.

## Relative mesh size is not a resolution test

A peak-cell length comparable to a global median does not prove adequate or
inadequate resolution. Resolution is relative to local passage width, shear-layer
thickness, curvature and other physical scales, and needs refinement evidence.
A scalar length also does not exclude every element-quality defect. Keep the
owner's peak-location/quality observations separate from a proven mechanism.
No new private-case diagnosis is adopted from these summaries.

## Claude handoff

The local source already exists at `/Users/gandanie/scratch/santis/OpenAccel`,
clean at `0d69041ba1afda63e9e4328d9e0d9834bba37756` when checked for this response.
Claude initialized its solver submodule during this review; it now checks out
the required `e351ba5eeaf3537dcc53d0aba09a8347f0a44cd0`. The temporary clone
`/private/tmp/mars-openaccel-liblin-stage0` already contains the pinned
`e351ba5eeaf3537dcc53d0aba09a8347f0a44cd0` Git object; it was cloned without a
worktree checkout. Verify these temporary objects still exist before using them.
Use an isolated reference build/worktree for instrumentation and preserve the pin.

Recommended next effort remains Medium, one Claude implementation agent.
Escalate only a named build/dependency obstacle if needed. No production SIMPLE
kernel work is requested before the reference outputs can be reviewed.

Both agents share this workspace: Codex can already read an uncommitted
`SYNC.md` update. A push is needed for another checkout, not for local visibility.
`SYNC.md` is currently excluded through `.git/info/exclude`; preserve that policy
unless the user changes it. Put public detailed handoffs in tracked documents
and update the local sync pointer. Neither agent should stage the other's files.
