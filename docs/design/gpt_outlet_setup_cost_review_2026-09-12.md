# Outlet tolerance experiment and repeated Hypre setup

Author: GPT/Codex. Date: 2026-09-12.
Status: complete public-log comparison and targeted source inspection; no timings
of individual phases and no production changes.

Follow-up: [opt-in profiling](gpt_outlet_profiling_2026-09-12.md) now implements
the measurement points below. The complete public GPU result is now recorded
there: after startup, solve is 72.9% of correction time, preparation/setup 24.1%,
and matrix assembly 0.038%.

The looser-tolerance run is 5.25% faster, but it does not remove 45% of the total
Krylov work. It needs an additional pressure correction on 74 of its 200 steps.
Full-run counts contradict the conclusion drawn from its final-step excerpt.

## Executed log comparison

Both complete logs were retrieved with authorized rsync from the user's Daint
build directory and parsed locally. Their headers show the same public channel,
one rank, rho=1000, nu=1e-4, dt=2e-6, beta=1, skew advection, relax_u=0.3,
relax_mass=1 and 200 requested steps. The user identifies the second run as
`--tol=1e-6`; neither header embeds the full command or executable revision.

| Whole-run quantity | Fixed-trace baseline | Looser tolerance |
|---|---:|---:|
| Runtime, ms | 56366.0 | 53406.2 |
| Physical steps | 200 | 200 |
| Pressure correction solves | 200 | 274 |
| FGMRES iterations / Hypre preconditioner calls | 3571 | 4009 |
| Inner Hypre iterations | 74280 | 71108 |
| Momentum solve calls | 600 | 600 |
| Momentum iterations | 997 | 798 |
| Final full continuity RMS, 1/s | 3.86e-13 | 9.69e-11 |
| Final full continuity maximum, 1/s | 2.43e-12 | 8.12e-10 |
| Final signed boundary imbalance, volume/s | 4.22e-15 | 3.79e-12 |

The parser grouped solve-trace entries before each accepted correction. Every
first correction has exactly three momentum calls plus the reported number of
FGMRES preconditioner calls; every second correction has only the latter.
The baseline has 200 first corrections and no second corrections; the looser
run has 200 first corrections and 74 second corrections. No failed-solve or
nonfinite marker appears in either complete log.

The inner Hypre iteration reduction is 4.27%, while its solve/setup call count
increases 12.27%. Iterations do not all have identical cost. These totals alone
cannot apportion wall time between assembly, setup, solves, synchronization and
diagnostics. Both runs enable solve tracing. Final printed velocity summaries
agree; this is not proof of unchanged field accuracy.

Local evidence directory: `/private/tmp/mars-outlet-audit-20260911/`.
Remote directory: `$SCRATCH/git/mars/daint-gpu/`.

```
channel-C-fixedtrace.log
SHA256 174ea09c5dc82a5e23ac5c2e9cdc6c010e245e8db59efad507ed5f327f43c5d8
channel-C-tol1e6.log
SHA256 d3051bca872de3f580a0c605753cd734451f8db5a0746c460fd13781a891fcdf
```

## Where setup actually occurs

The source path is:

1. `mars_outlet_correction.hpp` builds the step context and calls
   `refresh_outlet_pressure_operator` once before its correction loop.
2. `mars_ns_pump_solver.hpp:3800` restores bare K, adds compact derivatives,
   checks the anchor collectively, reapplies BC rows and copies values to Apre.
   This function does not construct Hypre or set up AMG.
3. Every outer FGMRES preconditioner application calls `solveOneComponent`
   through `mars_outlet_krylov.hpp:174`.
4. The GMRES branch of `solveOneComponent` constructs a local HypreGMRESSolver.
5. Its device-map `solve` overload calls `destroy()` at
   `mars_hypre_gmres_solver.hpp:77`, rebuilds the Hypre matrix/vectors, and calls
   `HYPRE_ParCSRGMRESSetup` (or FlexGMRESSetup) at lines 486-488 before solving.
   The default BoomerAMG preconditioner is registered with its setup callback.

Thus the current default pressure preconditioner executes setup 3571 times in
the baseline, not just 200 times. Caching only `refresh_outlet_pressure_operator`
would leave those setups untouched. Keeping the solver object alive alone also
does not suffice: the current `solve` API destroys its state on each call.

## Safe scope for subsequent work

The constancy argument is useful under the current supported configuration:
after BDF startup, fixed dt and fixed momentum matrices give fixed nodal tau and
compact derivative values. The old pressure gradient and trace still change;
their per-step context refresh cannot be skipped wholesale. Atomic reassembly
may differ at roundoff even when mathematical entries are unchanged.

A reusable prepared Hypre solve needs explicit invalidation for matrix values,
sparsity/partition, boundary masks, coefficient mode and solver/preconditioner
settings. BDF1-to-BDF2 is a mandatory invalidation. A pointer to an unchanged
allocation is not evidence that its contents are unchanged. Keep the ownership
and invalidation rules collective across ranks.

Recommended next step: Medium effort, one agent, gated phase timings and call
counts separating outlet assembly, Hypre matrix/vector preparation, Krylov/AMG
setup, inner solve, and total correction. Synchronize GPU phase boundaries only
when profiling is enabled and report rank-max elapsed times outside inner loops.
The setup path is demonstrably repeated; its share of wall time remains unmeasured.
Use those measurements to scope a prepared-solve API, then compare unchanged
continuity/flux gates and one-/multi-rank results before claiming a speedup.
