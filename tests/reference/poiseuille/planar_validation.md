# Planar Poiseuille release validation

The repaired path is `mars_poiseuille_flow --planar-ddt`. It uses the native
MARS reader and GPU Hex8 projection solver on the public mesh in
`tests/data/poiseuille/poiseuille_hex_14k_elem.e`. It is restricted to rectangular
cells and one element through z, and runs on any number of ranks (validated on
1, 2 and 4; see [Multi-rank result](#multi-rank-result)). This is separate from
the Tet4 SIMPLE solver. General 3D symmetry is not implemented.

The full 1500-step Daint run passed on 2026-09-27, as reported in the user's
terminal transcript for `poiseuille-release-nWKnbQ`. This closes the single-rank
planar regression gate. The saved revision and binary hashes were not retrieved;
the transcript follows the release recipe supplied after repair `cdffd654`.

## Recorded CUDA result

At t=15, the run printed `VALIDATION PASS` with the unchanged acceptance limits:

| Quantity | Reported value | Acceptance limit |
|----------|----------------|------------------|
| Profile RMS error | 4.562586e-4 m/s | < 6e-3 m/s |
| Profile RMS / analytic peak speed | 3.041724e-4 | Diagnostic |
| Through-flow ratios at 25/50/75% | 1.000 / 1.000 / 1.000 (printed) | 0.9 to 1.1 |
| Active continuity RMS * H/U | 7.3849e-9 | <= 1e-6 |
| Relative net boundary flux | 2.4104e-10 | <= 1e-6 |
| Final 20-step velocity change / U | 3.10352119e-7 | <= 1e-6 |

The final projection gate also passed. The velocity-fit pressure gradient was
0.11992 Pa/m against 0.12 Pa/m analytically; the solved-pressure gradient was
0.12178 Pa/m (about 1.48% high). These gradient reports are diagnostics, not
additional acceptance tests. Pressure CG took 4202 iterations in the final step;
this result establishes correctness for the stated case, not solver efficiency.

## Multi-rank result

Commit 4a2b5a62 builds the opening and cut-plane areas from each rank's owned
element faces on the GPU (reverse halo for shared nodes) and counts only owned
nodes in every validation sum. The same 1500-step case ran on 1, 2 and 4 GPUs of
one Daint node on 2026-09-27, from `/capstor/scratch/cscs/gandanie/poiseuille-mr`:

| Quantity | 1 rank | 2 ranks | 4 ranks |
|----------|--------|---------|---------|
| Profile RMS error | 4.562586e-4 | 4.562586e-4 | 4.562586e-4 |
| Through-flow ratios | 1.000 x3 | 1.000 x3 | 1.000 x3 |
| Continuity * H/U | 7.4569e-9 | 6.6959e-9 | 7.1049e-9 |
| Relative net boundary flux | 2.3331e-10 | 1.9911e-10 | 1.9827e-10 |
| Final 20-step change / U | 3.10352183e-7 | 3.10352133e-7 | 3.10352144e-7 |
| Wall time (sacct elapsed) | 0:25:39 | 1:09:19 | 2:01:46 |

All three print `VALIDATION PASS`. Startup is rank-invariant: 30000 owned nodes,
volume 0.6, 796 velocity-Dirichlet nodes, inlet and outlet areas 0.06, and the
same projection identity at steps 1 to 3. The 1-rank rerun reproduces the
single-rank result above to all printed digits.

Final fields, compared node by node with `compare_fields.py` (pieces merged by
coordinates; ghost copies must be identical):

| Pair | max \|du\|/U | max \|dp\| off the inlet plane | max \|dp\| on the inlet plane |
|------|---------|------------------|-----------------|
| 2 vs 1 rank | 2.0e-9 | 2.3e-7 Pa | 1.7e-5 Pa |
| 4 vs 1 rank | 1.9e-9 | 2.3e-7 Pa | 1.7e-5 Pa |
| 4 vs 2 ranks | 1.8e-9 | 4.5e-9 Pa | 1.0e-6 Pa |

At steps 500 and 1000 the partition-seam nodes differed by at most 2e-8 Pa, so
the gap is not at the seams. It fits the pressure preconditioner, although no run
has isolated it: one rank uses the diagonal of the assembled matrix, several ranks
the matrix-free diagonal (`hexMultirankDiag` in `mars_ns_channel_solver.hpp`), and
CG stops at a relative residual of 1e-6. Two and four ranks share the
preconditioner and agree more closely. The inlet pressure next to the walls keeps growing with time on every
rank count (2.1e4 Pa at step 1500); the interior pressure is 0 to 1.96 Pa. This
does not enter the acceptance metrics.

Compare two runs:

```bash
python3 tests/reference/poiseuille/compare_fields.py \
  <run1>/fields_step1500.pvtu <run2>/fields_step1500.pvtu
```

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
