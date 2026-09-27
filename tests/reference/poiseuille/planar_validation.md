# Planar Poiseuille release validation

The repaired path is `mars_poiseuille_flow --planar-ddt`. It uses the native
MARS reader and GPU Hex8 projection solver on the public mesh in
`tests/data/poiseuille/poiseuille_hex_14k_elem.e`. It is restricted to one rank,
rectangular cells and one element through z. This is separate from the Tet4
SIMPLE solver. General 3D symmetry and distributed opening-area reconstruction
are not implemented by this change.

The three-step GPU startup gate passed previously. A full steady-profile pass
on the committed repair is still required before declaring this regression fixed.

## Discrete equations

Pressure and its increment are in Pa; density is in kg/m^3 and `nu` is kinematic
viscosity in m^2/s. The validation uses rho=1, nu=.01, U=1, dt=.01. Velocity is
prescribed at the inlet, zero on the y walls, free at the p=0 outlet, and w=0
on both z planes. No-slip takes precedence at inlet-wall corners.

Let D be the unnormalised SCS divergence, M the lumped nodal volume, Q the
free-velocity mask and P the active pressure-row selection. The pressure action is

```
A = P D Q M^-1 D^T P^T
A phi = -(rho/dt_eff) P R(u**)
u_new = u** + (dt_eff/rho) Q M^-1 D^T P^T phi
p_new = p_old + phi
```

`dt_eff` is dt in BDF1 and 2dt/3 in BDF2. The predictor uses the same adjoint
gradient. R includes actual nodal inlet targets and solved outlet flux, without
outlet rescaling. The pressure mask removes outlet rows and inlet-corner rows
with no free velocity columns. Both the pressure operator and correction apply Q.
Diffusion restores the nonzero inlet lift before velocity column elimination.
Skew advection includes the matching opening momentum flux.

These changes are in `mars_ns_channel_solver.hpp` and
`mars_channel_projection.hpp`; the separate pump and SIMPLE solvers are unchanged.

## Local arithmetic checks

From the repository root:

```bash
python3 scripts/cvfem_gradient_host_check.py
python3 scripts/channel_bdf_scaling_check.py
python3 scripts/channel_projection_host_check.py
python3 scripts/vtu_precision_host_check.py
```

These compile production arithmetic with host launch stubs. They check gradients,
the constrained projection identity, symmetry/positivity on small invented grids,
opening flux cancellation, nonzero inlet lift, skew-advection energy, BDF startup,
DOF permutation and failure of invalid convergence metrics. The VTU check verifies
Float64 field/coordinate/time round-trips using invented data. They do not execute CUDA.
Temporary build files stay under the repository's `.local-worktrees/poiseuille-host/`.

## Daint release run

Use the existing MARS CUDA environment and configured `mars-v010-check/build-hypre`
directory. No OpenAccel build or environment change is needed. The public mesh is
tracked, so there is no Python preparation step. Pull and build this target:

```bash
git pull --ff-only &&
cmake --build . --target mars_poiseuille_flow --parallel 4
```

Then run interactively:

```bash
(
set -euo pipefail
poiseuille_run=$(mktemp -d /capstor/scratch/cscs/gandanie/poiseuille-release-XXXXXX)
git rev-parse HEAD > "$poiseuille_run/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_poiseuille_flow \
  ../tests/data/poiseuille/poiseuille_hex_14k_elem.e > "$poiseuille_run/sha256.txt"
printf 'Results: %s\n' "$poiseuille_run"

MARS_NS_DEBUG_STEPS=3 MARS_DDT_CG_PRINT_EVERY=1000 \
srun --account=csstaff --time=04:00:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_poiseuille_flow \
  --mesh=../tests/data/poiseuille/poiseuille_hex_14k_elem.e \
  --planar-ddt --uinf=1 --rho=1 --nu=0.01 --dt=0.01 \
  --tol=1e-6 --max-iter=8000 --num-steps=1500 --check \
  --vtu-output="$poiseuille_run/flow" --vtu-every=500 \
  --comparison-output="$poiseuille_run/fields" --comparison-every=50 \
  2>&1 | tee "$poiseuille_run/run.log"
)
```

The existing binding helper is read from home; generated files all go to capstor.
The wall-time request is a budget, not a measured runtime for the repair.

`VALIDATION PASS` requires all of:

- The unchanged downstream nodal profile RMS below .006 m/s.
- The unchanged 25/50/75% through-flow ratios within 10% of actual inlet flow.
- Volume-weighted active continuity RMS times H/U at most 1e-6.
- Absolute net boundary flux divided by actual inlet flux at most 1e-6.
- Volume-weighted velocity change over the final 20 steps, divided by U, at most 1e-6.

Failed momentum or pressure solves abort before their output reaches the next stage.
Startup and final projection checks also verify zero normal velocity, compatible
corner rows and the boundary-flux sum identity. The solved pressure-gradient report
is a separate diagnostic; an empty moving-core sample prints NaN, not zero.

`--check` returns nonzero if the fixed run has not met all criteria. Preserve its
log and Float64 fields when it fails. No release claim follows from exit zero
without the validation line. The equivalent long CTest is enabled with
`MARS_ENABLE_VALIDATION_TESTS=ON`; invoke it only inside a GPU allocation.
