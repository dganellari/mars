# Baseline measurements that bear on the port plan

Author: Claude. Date: 2026-09-12.
Status: measurements from user-executed Alps runs plus archived project records.
Three items correct or extend claims in the GPU port plan and BC contract.
No confidential geometry or field data inspected; pump numbers are printed
scalars from the driver's own diagnostics.

## 1. `nu=1e-4` was not a load-bearing workaround

Port plan line 23 and BC contract lines 30-32 state that `nu=1e-4` was an
intentional workaround because water viscosity could not previously be run.
Three independent measurements contradict this:

| | `nu=1e-4` | `nu=1e-6` |
|---|---|---|
| public channel, step 200 | `u_rms=6.615e-01 u_max=1.233`, cuts `0.500`, RMS `5.33e-13` | identical at printed precision, RMS `5.36e-13` |
| pump, step 30 raw cuts | `-2.937e-05 / -3.239e-05 / -3.110e-05` | identical |
| pump, step 70 | `u_rms=1.519e-01 u_max=8.463` | `u_rms=1.519e-01 u_max=8.469` (0.07%) |

`mars_pump.cu:133` declares `nu = 1.0e-6` as the driver default, commented
"water kinematic viscosity". Archived finding, 2026-08-27: "nu 1e-4 vs 1e-5
agree to 0.5% on EVERY percentile, so the scheme's numerical dissipation
dominates nu ... and nu=1e-6 is not worth running."

A workaround that is load-bearing fails when removed. This one changes nothing,
on the real geometry, at the target viscosity. Consequence for the plan: Stage 5's
exit evidence "converged public cases at water properties" is satisfied trivially
and has no discriminating power. It needs a case where viscosity is dynamically
active.

Related: the BC contract reports the water GPU transcript (38.5 ms/step, RMS
5.36e-13, max speed 1.233, cuts 0.500) without noting the `nu=1e-4` clean run
gave 40.0 ms/step, 5.33e-13, 1.233, 0.500. They are the same result.

## 2. No sub-grid or turbulence item exists in the plan

The archived blocker for real water on the pump is not viscosity:

> the wall is NOT pressure, NOT geometry/side-sets -- it is the HIGH-Re ADVECTION
> INSTABILITY: no sub-grid dissipation, so at the Re needed for real flow it piles
> energy at the grid scale and blows up. 430 m/s needs LES/turbulence model +
> near-passage AMR.

Neither document has a sub-grid, turbulence or resolution item. Stage 0 item 3
covers advection, diffusion, linearization and under-relaxation, but not what the
reference does that lets it run at Re~2e5. If OpenAccel succeeds on the PI's case
it is doing something specific there -- a turbulence model, a more dissipative
advection scheme, or a much finer mesh. Stages 0-6 could reach term-by-term parity
on a public channel and the pump still blow up.

Suggested: add it to the Stage 0 contract's required identifications.

## 3. RETRACTED: the cg_p gap is a work-unit artifact, not a defect

This section previously reported that `MARS_OUTLET_PRECOND=legacy` costs the pressure
Poisson solve ~12x its iterations (`cg_p` 288-465 vs 20-32) and suspected a latent
defect in the Hypre wrapper. That is wrong. GPT/Codex identified the cause:

`mars_outlet_krylov.hpp:218` does `s.lastPressureIters += iterations;` inside the
outlet preconditioner apply, and the same function calls `solveOneComponent` in the
legacy branch. So `cg_p` is the pressure Poisson count PLUS accumulated outlet
preconditioner work -- reported as inner GMRES iterations in legacy mode and as AMG
cycles in amg-cycle mode. Different units; the ratio measures nothing.

My supporting analysis was wrong because two greps were mis-scoped: the pattern
`lastPressureIters *=` cannot match `+=`, and I searched the correction and flux
files but not the krylov file where both the accumulation and the legacy call live.
Null grep results were treated as proof of absence.

No Hypre wrapper defect is claimed. The amg-cycle timing improvement (278.1 -> 38.5
ms/step) is a separate, still-valid wall-clock measurement.

## 4. Pump velocity peak: not a mesh defect, probably under-resolved

300-step pump probe, `MARS_UMAX_PROFILE=1`, water properties, 4 ranks:

- `max/p999` flat at ~2.1 across 150 steps; `peak-node=interior`; `h_peak/h_med=1.00`.
- Percentile growth exponents in step number: p50 1.73, p90 1.56, p99 1.39,
  p999 1.35, max 1.37, u_rms 1.23, Q_out 1.00 (exact). Lower percentiles grow
  fastest -- the distribution narrows; `p99/p50` falls 42.0 to 12.2.
- Outlet clean throughout: mass balance -0.063% to -0.003%, RC `0.000%`,
  `corrections=1`, `rms~1e-10`.

At >1 rank the percentiles are `MPI_MAX`-reduced so their levels are upper bounds;
`max`, `peak-node` and `h_peak` are exact (owner-broadcast, `mars_ns_pump_solver.hpp`
:12286-12294). `h_med` is an upper bound, so `h_peak/h_med=1.00` is a LOWER bound.
Per that function's own comment, a ratio at or above 1 means the peak sits where the
mesh is coarse and "is probably not resolved" -- consistent with item 2, and not a
sliver artifact.

Open: whether the gain `u_rms/Q_out` plateaus. It ran 901 to 1720 through step 150,
increments decelerating but non-monotone over uneven intervals. A one-rank long run
is in progress; `scripts/pump_steady_trend.py` extracts the trend.
