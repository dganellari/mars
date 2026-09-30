# Planar Poiseuille release validation

`mars_poiseuille_flow` runs the CVFEM Navier-Stokes solver of
`backend/distributed/unstructured/fem/mars_navier_stokes.hpp` (the same solver as
the Taylor-Green vortex and the lid-driven cavity) in planar mode on the public mesh
`tests/data/poiseuille/poiseuille_hex_14k_elem.e`. Planar mode needs Hex8 cells with
one element through z and runs on any number of ranks. Both linear systems are
assembled once and solved with Hypre PCG + BoomerAMG. This is separate from the
Tet4 SIMPLE solver.

## Result, 2026-09-30

1500 steps (t = 15), rho = 1, nu = 0.01, U = 1, dt = 0.01, relative solver
tolerance 1e-10, on one Daint node. All three runs print `VALIDATION PASS`:

| Quantity | 1 GPU | 2 GPUs | 4 GPUs | Limit |
|----------|-------|--------|--------|-------|
| Profile RMS error (m/s) | 4.551373e-4 | 4.551373e-4 | 4.551375e-4 | < 6e-3 |
| Through-flow ratios at 25/50/75% | 1.0000 x3 | 1.0000 x3 | 1.0000 x3 | 0.9 to 1.1 |
| Continuity RMS * H/U | 8.0e-15 | 1.7e-14 | 1.1e-14 | <= 1e-6 |
| Relative net boundary flux | 8.8e-8 | 8.8e-8 | 8.8e-8 | <= 1e-6 |
| Final 20-step velocity change / U | 2.6e-8 | 2.6e-8 | 2.6e-8 | <= 1e-6 |
| Projection identity, step 1500 | 3.2e-12 | 3.3e-12 | 3.3e-12 | <= 1e-7 |
| AMG iterations u/v/p, step 1500 | 10/13/18 | 10/13/18 | 11/13/18 | |
| Time per step (ms) | 15.9 | 63.8 | 82.5 | |

The pressure gradient from the solved pressure is 0.12178 Pa/m and from the
velocity profile 0.11992 Pa/m, against 0.12 Pa/m exactly; these are diagnostics,
not acceptance tests. At setup the assembled pressure matrix matches the
matrix-free flux correction to 1.1e-15 to 1.3e-15 (relative, max norm). On this
30k-node mesh more GPUs are slower; scaling runs use the generated channel
(`--cells`).

The net boundary flux is not at roundoff because the outlet flux uses the nodal
velocity, which follows the corrected face fluxes only up to the Rhie-Chow term.
The face fluxes themselves balance to roundoff (the continuity row).

### Earlier records

The exact projection (pressure matrix D Q M^-1 D^T, no stabilization) passed the
same checks on 1, 2 and 4 GPUs earlier on 2026-09-30, with profile RMS
4.553029e-4 m/s, a net boundary flux of 5e-16 and 15 to 16 pressure iterations. It was
replaced because its pressure matrix has checkerboard null modes on closed and
periodic boxes, where multigrid breaks down. The stabilized solver changes the
profile RMS by 1.7e-7 m/s.

With the older solver and Jacobi-preconditioned CG (tolerance 1e-6) on 2026-09-27:
profile RMS 4.562586e-4 m/s, continuity * H/U about 7e-9, 20-step change 3.1e-7,
4202 pressure iterations in the final step, 25 minutes to 2 hours per run.

## Discrete equations

Pressure and its increment are in Pa; density is in kg/m^3 and `nu` is kinematic
viscosity in m^2/s. Velocity is prescribed at the inlet, zero on the y walls, free
at the p = 0 outlet, and w = 0 everywhere. No-slip takes precedence at inlet-wall
corners.

Each hex has 12 sub-control faces f with area vector A_f from node L to node R.
D_F is the face flux divergence, M the lumped nodal volume, K the CVFEM Laplacian
built from the compact face gradients, G p = -M^-1 D^T p the nodal gradient, Q the
free-velocity mask, and h = dt_eff/rho. The face flux is stabilized (Rhie-Chow):

```
F = A_f . (u_L + u_R)/2 - h [ (grad p . A)_f - A_f . (G p_L + G p_R)/2 ]
K phi = -(rho/dt_eff) R(F**)
F_new = F** - h (grad phi . A)_f
u_new = u** - h Q G phi
p_new = p_old + phi
```

`dt_eff` is dt in BDF1 and 2dt/3 in BDF2. F** is the stabilized flux of u** and
p_old. The predictor uses Q G p_old. R adds the opening fluxes to D_F: prescribed
nodal inlet velocities and the computed outlet velocity, without rescaling. The
outlet rows of K (p = 0) are removed. Diffusion moves the fixed inlet values to
the right-hand side (lift) before eliminating their columns. Skew advection uses
the face fluxes F and includes the matching opening momentum flux.

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
- The projection check at steps 1, 2, 3 and 1500: D_F F = D_F F** + h K phi to
  1e-7, and the boundary flux balance consistent with it to 1e-10.

A failed linear solve stops the run with exit code 1. The equivalent CTest,
`marsPoiseuilleValidation`, is enabled with `MARS_ENABLE_VALIDATION_TESTS=ON`; run
it only inside a GPU allocation.
