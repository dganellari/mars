# Planar Poiseuille release validation

`mars_poiseuille_flow` runs the planar CVFEM projection solver of
`backend/distributed/unstructured/fem/mars_channel_flow.hpp` on the public mesh
`tests/data/poiseuille/poiseuille_hex_14k_elem.e`. The solver accepts axis-aligned
Hex8 cells with one element through z and runs on any number of ranks. Both linear
systems are assembled once and solved with Hypre PCG + BoomerAMG. This is separate
from the Tet4 SIMPLE solver.

## Result, 2026-09-30

1500 steps (t = 15), rho = 1, nu = 0.01, U = 1, dt = 0.01, relative solver
tolerance 1e-10, on one Daint node. All three runs print `VALIDATION PASS`:

| Quantity | 1 GPU | 2 GPUs | 4 GPUs | Limit |
|----------|-------|--------|--------|-------|
| Profile RMS error (m/s) | 4.553029e-4 | 4.553029e-4 | 4.553029e-4 | < 6e-3 |
| Through-flow ratios at 25/50/75% | 1.0000 x3 | 1.0000 x3 | 1.0000 x3 | 0.9 to 1.1 |
| Continuity RMS * H/U | 6.0e-13 | 2.3e-13 | 1.5e-13 | <= 1e-6 |
| Relative net boundary flux | 4.6e-16 | 1.2e-16 | 4.6e-16 | <= 1e-6 |
| Final 20-step velocity change / U | 3.0e-10 | 3.0e-10 | 3.0e-10 | <= 1e-6 |
| Projection identity, step 1500 | 2.09e-13 | 2.10e-13 | 2.07e-13 | <= 1e-7 |
| AMG iterations u/v/p, step 1500 | 15/14/15 | 14/14/16 | 15/14/16 | |
| Time per step (ms) | 16.3 | 65.3 | 87.6 | |

The pressure gradient from the solved pressure is 0.12178 Pa/m and from the
velocity profile 0.11992 Pa/m, against 0.12 Pa/m exactly; these are diagnostics,
not acceptance tests. At setup the assembled pressure operator matches the
matrix-free D Q M^-1 D^T to 1.4e-15 (relative, max norm). On this 30k-node mesh
more GPUs are slower; scaling runs use the generated channel (`--cells`).

### Equivalence with the previous solver

The rewrite replaced a 10.9k-line solver fork. Both ran 200 steps with the same
AMG solvers and tolerance 1e-10 on 1, 2 and 4 GPUs, and the fields were compared
node by node with `compare_fields.py`:

| | max \|du\| / U | max \|dp\| off the inlet plane | max \|dp\| on the inlet plane |
|---|---|---|---|
| old vs new, 1 / 2 / 4 GPUs | 1.0e-9 / 3.6e-10 / 1.6e-10 | 6.5e-10 / 2.0e-10 / 1.5e-10 Pa | 1.0e-5 / 5.2e-6 / 6.2e-6 Pa |
| old, 1 vs 4 GPUs | 9.8e-10 | 7.9e-10 Pa | 9.7e-6 Pa |

The inlet-plane differences sit where the pressure next to the walls is large
(3.5e3 Pa at step 200) and equal the old solver's own difference between rank counts.

### Earlier record

The same checks passed with the previous solver and Jacobi-preconditioned CG
(tolerance 1e-6) on 1, 2 and 4 ranks on 2026-09-27: profile RMS 4.562586e-4 m/s,
continuity * H/U about 7e-9, 20-step change 3.1e-7, 4202 pressure iterations in the
final step, 25 minutes to 2 hours per run.

The inlet pressure next to the walls keeps growing with time on every rank count.
It does not enter the acceptance metrics.

## Discrete equations

Pressure and its increment are in Pa; density is in kg/m^3 and `nu` is kinematic
viscosity in m^2/s. Velocity is prescribed at the inlet, zero on the y walls, free
at the p = 0 outlet, and w = 0 everywhere. No-slip takes precedence at inlet-wall
corners.

Let D be the unnormalised sub-control-face divergence, M the lumped nodal volume,
Q the free-velocity mask and P the active pressure-row selection. The pressure
action is

```
A = P D Q M^-1 D^T P^T
A phi = -(rho/dt_eff) P R(u**)
u_new = u** + (dt_eff/rho) Q M^-1 D^T P^T phi
p_new = p_old + phi
```

`dt_eff` is dt in BDF1 and 2dt/3 in BDF2. The predictor uses the same adjoint
gradient. R adds the opening fluxes: prescribed nodal inlet velocities and the
computed outlet velocity, without rescaling. The pressure mask removes outlet rows
and inlet-corner rows, whose columns of A are zero. Diffusion moves the fixed inlet
values to the right-hand side (lift) before eliminating their columns. Skew
advection includes the matching opening momentum flux.

## Daint release run

From the build directory of a CUDA build with `MARS_ENABLE_HYPRE=ON`:

```bash
cmake --build . --target mars_poiseuille_flow --parallel 32
R=$(mktemp -d /capstor/scratch/cscs/gandanie/poiseuille-release-XXXXXX)
git rev-parse HEAD > $R/mars-revision.txt
for np in 1 2 4; do
  srun --account=csstaff --time=00:30:00 --nodes=1 --ntasks-per-node=$np --export=ALL --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_poiseuille_flow \
    --mesh=../tests/data/poiseuille/poiseuille_hex_14k_elem.e \
    --uinf=1 --rho=1 --nu=0.01 --dt=0.01 --num-steps=1500 --report-every=100 --check \
    --vtu-output=$R/np$np/flow --vtu-every=500 --comparison-output=$R/np$np/fields --comparison-every=500 \
    > $R/np$np.log 2>&1
  echo "np$np exit $?"
done
```

Compare two runs:

```bash
python3 tests/reference/poiseuille/compare_fields.py <run1>/fields_step1500.pvtu <run2>/fields_step1500.pvtu
```

`VALIDATION PASS` requires all of:

- The downstream nodal profile RMS below .006 m/s.
- The 25/50/75% through-flow ratios within 10% of the inlet flow.
- Volume-weighted continuity RMS times H/U at most 1e-6.
- Absolute net boundary flux divided by the inlet flux at most 1e-6.
- Volume-weighted velocity change over the final 20 steps, divided by U, at most 1e-6.
- The projection check at steps 1, 2, 3 and 1500: D u = D u** + h A phi to 1e-7,
  flux balance to 1e-10, zero continuity at the inlet-wall corners.

A failed linear solve stops the run with exit code 1. The equivalent CTest,
`marsPoiseuilleValidation`, is enabled with `MARS_ENABLE_VALIDATION_TESTS=ON`; run
it only inside a GPU allocation.
