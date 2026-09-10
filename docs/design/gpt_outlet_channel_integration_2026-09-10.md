# Public channel: full outlet-stepper integration gate

Author: GPT/Codex. Date: 2026-09-10.
Status: implemented; mesh read-back and log-checker fault tests pass locally.
Compilation of the new driver check and integrated CUDA/MPI runs remain pending.

Update: the user compiled and ran the first public-channel step on Daint. Momentum
converged, but the pressure line search rejected all backtracks before any channel
step passed. The [line-search diagnosis and fix](gpt_outlet_line_search_fix_2026-09-10.md)
records the confirmed acceptance-rule defect and the exact SFC volume reference.
The corrected integrated run is pending; this is not a passing flow result.

## Equations and bounded purpose

The existing driver solves incompressible Navier–Stokes with physical pressure:
`rho*(du/dt + u.grad(u)) = -grad(p) + rho*nu*Laplacian(u)`, `div(u)=0`.
The test sets `rho=1`, `nu=0.1`, and thus dynamic viscosity `mu=0.1` in consistent
SI units. The public channel is `[0,4] x [0,1] x [0,1]` metres. Inlet speed ramps
to `0.1 m/s` over four steps; the Reynolds number using channel width is 1.
No bbox-based `--Re` override is used.

Velocity starts from rest with the driver's existing boundary initialization.
The inlet prescribes axial velocity. The four side walls use no-slip; inlet/outlet
nodes retain the driver's existing precedence at wall intersections. The outlet
velocity is free and its frozen pressure trace has mean zero and beta 0.05.
Consequently this is a transient rectangular-duct integration test with the driver's
existing corner treatment, not the parabolic-inlet manufactured channel in the spec.
There is no claimed analytic velocity or pressure solution for these eight steps.

The production nodal Tet4 CVFEM/SCS discretization, lumped mass, per-vertex opening
samples with triangle area/3, nodal VMS coefficients, and boundary-aware gradients
are reused unchanged. `dt=0.01 s`; the first step must use BDF1 (`dtEff=dt`), followed
by BDF2 (`dtEff=2*dt/3`). The full Hypre correction path is used, including its
measured damping and failure conditions. No solver terms or tolerance defaults change.

## Checks and limits

The opt-in `--outlet-channel-check` captures velocity before each physical step and
checks the stored BDF history afterward. History equality is exact because this is
a copy invariant, not a floating-point reduction. It forward-exchanges copies of
the solved fields using the actual domain halo and requires that no ghost value
changes. These checks use persistent device scratch and scalar host reductions.

Every step reassembles the actual full continuity residual with the frozen context
and real reverse halo. Require RMS and max <= `1e-7 /s`, boundary imbalance <=
`1e-8 m^3/s`, and agreement between the residual sum and boundary reporter within
`1e-10 m^3/s`. The prescribed inlet flux is independently `-U_in*1 m^2` at the
current ramp. Require positive outlet flow, total lumped volume matching the
production SFC-encoded/decoded box (4.0000005722048915 for the current 64-bit keys), scalar outlet
area 1, and the prescribed frozen trace mean. The checker prints full-precision
velocity/pressure means and RMS values for partition comparisons.

The log tool requires all steps and an explicit completion marker, rejects failures
and nonfinite data, and checks BDF startup, source flux, and residual tolerances.
Physical integral comparisons use `atol=1e-8`, `rtol=5e-6`; these tolerate solver
and reduction differences and are not a pointwise field comparison. Empty-opening
rank coverage is measured, not assumed from rank count. `--require-empty` fails if
none of the runs exercises it; additional partition coverage is then required.

The new scalar records (including empty-opening rank count) are explicitly gated.
Use the new flag only with this public fixture. No VTU output is enabled or needed.
Passing does not establish spatial accuracy, temporal order, fresh-trace convergence,
or performance. Those remain later gates.

## Daint run sequence

Pull `cstone` and rebuild `mars_pump` using the existing CUDA/Hypre/netCDF build.
No new target or Python package is required. Run from `daint-gpu/`.
Stop at the first failure; do not increase damping or relax tolerances to hide it.
The shell `pipefail` setting preserves `srun` failure when output is saved with `tee`.

```bash
set -o pipefail
MARS_SOLVE_TRACE=1 MARS_NS_DEBUG_STEPS=1 srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_pump \
  --mesh=../tests/data/public_outlet_channel/outlet_channel.exo --inlet-ss=inlet --outlet-ss=outlet \
  --solver=hypre --vms-stab --rc-implicit --outlet=do-nothing --outlet-beta=0.05 --outlet-channel-check \
  --rho=1 --nu=0.1 --inlet-velocity=0.1 --dt=0.01 --num-steps=8 --source-ramp-steps=4 \
  --relax-u=0.3 --relax-mass=1 --tol=1e-11 --max-iter=2000 \
  --outlet-max-corrections=200 --outlet-rtol=1e-9 --outlet-div-tol=1e-8 --outlet-flux-tol=1e-10 \
  2>&1 | tee outlet-channel-1.log
```

After the one-rank run passes:

```bash
MARS_NS_DEBUG_STEPS=1 srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=2 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_pump \
  --mesh=../tests/data/public_outlet_channel/outlet_channel.exo --inlet-ss=inlet --outlet-ss=outlet \
  --solver=hypre --vms-stab --rc-implicit --outlet=do-nothing --outlet-beta=0.05 --outlet-channel-check \
  --rho=1 --nu=0.1 --inlet-velocity=0.1 --dt=0.01 --num-steps=8 --source-ramp-steps=4 \
  --relax-u=0.3 --relax-mass=1 --tol=1e-11 --max-iter=2000 \
  --outlet-max-corrections=200 --outlet-rtol=1e-9 --outlet-div-tol=1e-8 --outlet-flux-tol=1e-10 \
  2>&1 | tee outlet-channel-2.log

MARS_NS_DEBUG_STEPS=1 srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=4 --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_pump \
  --mesh=../tests/data/public_outlet_channel/outlet_channel.exo --inlet-ss=inlet --outlet-ss=outlet \
  --solver=hypre --vms-stab --rc-implicit --outlet=do-nothing --outlet-beta=0.05 --outlet-channel-check \
  --rho=1 --nu=0.1 --inlet-velocity=0.1 --dt=0.01 --num-steps=8 --source-ramp-steps=4 \
  --relax-u=0.3 --relax-mass=1 --tol=1e-11 --max-iter=2000 \
  --outlet-max-corrections=200 --outlet-rtol=1e-9 --outlet-div-tol=1e-8 --outlet-flux-tol=1e-10 \
  2>&1 | tee outlet-channel-4.log

python3 ../scripts/outlet_channel_check.py --require-empty \
  outlet-channel-1.log outlet-channel-2.log outlet-channel-4.log
```

All three solver runs must print `PASS: public outlet channel integration steps=8`.
The log comparison must pass too; if its only failure is missing empty-rank coverage,
record that limitation and prepare a partition that exercises it before claiming it.
