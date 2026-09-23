# Ordered SIMPLE update replay

This slice adds device-callable pressure, velocity, mass-flux and outlet-state
updates for the pinned OpenAccel reference. It does not expose a new flow solver
or change the existing projection solvers. The public fixture is the same
425-node, 1536-Tet4 channel used for the interior, node and boundary captures.

## Equations and state order

The target is steady, fixed-frame, incompressible flow:
`div(rho*u)=0` and `div(rho*u*u)=-grad(p)+div(mu*(grad(u)+grad(u)^T))`.
Pressure is physical pressure in Pa, density is kg/m³, dynamic viscosity is Pa s,
velocity is m/s, and stored flux is kg/s. Nodal momentum influence `d=V/a_hat`
is m³ s/kg and already includes momentum relaxation. There is no physical
`timestep/rho` factor in these updates.

At OpenAccel `0d69041ba1afda63e9e4328d9e0d9834bba37756`, the relevant sequence in
`src/equation/flow/segregatedFlowEquations.cpp::solve` is:

1. Momentum preparation, assembly and solve update the predicted velocity.
2. Pressure preparation refreshes the outlet trace from current pressure and
   the previous reversal flags. The trace stays fixed through this correction.
3. Solve pressure correction. The generic corrector writes
   `p_new=p_old+alpha_p*effective_increment`. The pressure-specific override
   separately stores the full raw linear-solve increment `p'`, publishes ghosts
   and reconstructs `grad(p')` with gradient relaxation one. The effective and
   raw increments are captured separately; acceleration is not implemented here.
4. Recompute and relax stored interior/inlet/outlet mass flux. This uses predicted
   velocity, **updated pressure and the old reconstructed pressure gradient**.
5. Update outlet reversal flags/flux, then assemble mass divergence.
6. Correct velocity componentwise: `u_new=u_pred-d_selected*grad(p')`.
   `d_selected` is `du` for SIMPLE and `duTilde` for SIMPLEC. There is no second
   `alpha_u`, and no `alpha_p` on this velocity increment.
7. Refresh the pressure gradient for the next momentum solve; refresh velocity
   scale/gradients/blending. These reconstructions remain reference inputs.

The capture records these phases at the actual calls. It is restricted to one
rank, steady 3D, fixed mesh/frame, incompressible flow and one pressure
subiteration. The public runner requires every expected node/sample in both
iterations. It does not certify a MARS driver implements that ordering yet.

## Discrete update contract

| Stage | Reference site | MARS arithmetic |
|---|---|---|
| 0 pressure | `equation.h::correctField_` | old plus relaxed effective increment |
| 1 raw increment | `pressureCorrectionEquation.h::correctField_` | store the unrelaxed solve value |
| 2 velocity | `segregatedFlowEquations.cpp::solve` | subtract selected influence times full increment gradient |
| 3 interior flux | `flowModel.cpp::updateMassFlowRateInterior_` | recompute RC flux and relax against stored flux |
| 4 inlet flux | `updateMassFlowRateBoundaryFieldInletSpecifiedVelocity_` | prescribed side velocity flux, then mass relaxation |
| 5 outlet flux | `updateMassFlowRateBoundaryFieldOutletSpecifiedPressure_` | mixed pressure/RC flux; reversed samples become zero |
| 6 reversal | outlet pressure branch of `updateFlowReversalFlag_` | face-level wall transition and sample clipping |
| 7 trace moments | constant/time-table `updatePressureBoundarySideFieldAverageStaticPressure_` | open sample pressure times scalar sample area |
| 8 trace application | same routine | nearest-node fluctuation about the patch mean |
| 9 patch mean | same routine | divide total pressure moment by total open area |

For an interior sample with oriented area vector `A`, the fresh mass flux is

```
m_fresh = rho * [u - d*(g_compact-g_reconstructed)
                  + d*(F_original-F_reconstructed)] dot A
m_new   = alpha_m*m_fresh + (1-alpha_m)*m_old
```

Products with `d` are componentwise. Body forces have pressure-gradient units.
The interior and outlet kernels preserve the respective reference operation
order. All interpolation, compact/reconstructed gradients, force averages and
sample density are frozen inputs to this gate. In particular, interior density
upwinding in the flux **update** uses the old stored-flux sign; the pressure
assembly uses a newly reconstructed flux sign. Do not blindly reuse its flux
routine when building interpolation for this update.

Both interior and outlet flux updates read `du`, even when SIMPLEC selects
`duTilde` elsewhere. The velocity update receives the selected coefficient.
These are distinct choices in the pinned code.

For an inlet, fresh flux is `rho*U_prescribed dot A`; solved nodal velocity is
not substituted. For an already reversed outlet sample, stored flux becomes
zero immediately, including its relaxed history.

The pressure outlet reversal rule first sums the **unclipped** three sample
fluxes, then clips each sample to nonnegative outflow. If the face was open and
its net flux is below `-SMALL`, all three flags become one and all three fluxes
become zero. A closed face reopens only when the face-average velocity points
outward and the face-average trace does not exceed face-average pressure.
Otherwise its flux remains zero. The `ignoreFlagUpdate` option preserves flags
but still executes the native flux clipping/zeroing. `SMALL` is FP64 epsilon
at this pin. Nonzero normals and uniform Tri3 face flags are caller preconditions.

The trace mean includes only samples not flagged as reversed, weighted by
`|A_s|`. The mean uses the reference pressure shape weights, while application
uses the nearest-node value:

```
mean = sum_open(p_interpolated,s * |A_s|) / sum_open(|A_s|)
t_s  = p_prescribed + (1-beta)*(p_nearest,s - mean)
```

A reversed sample retains its previous trace. If the whole patch is closed,
the pinned division has zero area; capture/replay rejects it. No fallback has
been invented. Interpolation of the resulting side trace to boundary nodes and
normal/corner boundary constraints are not part of this gate.

## Validation

Local checks at implementation time:

- 135 independent algebra checks, including pressure/velocity relaxation,
  compact-gradient sensitivity, balanced forces, scalar area weights,
  mixed-sign sample fluxes and wall creation/reopening.
- 18,000 comparisons against arithmetic extracted from the actual pinned
  functions: zero observed differences. These are synthetic states, not a CFD run.
- The actual inserted capture expressions compile with mock mesh plumbing.
  Their writer, strict parser and production host replay pass 276 comparisons
  in two synthetic ordered sequences. This is not an STK compilation.
- Transport tests reject incorrect phases, trace order, modified flux/flag
  history, missing samples, duplicate records, incorrect pins and nonfinite data.
  Runner tests reject missing update exports and preserve limited success claims.

Run the local checks with:

```bash
c++ -std=c++20 -Wall -Wextra -Werror \
  -Ibackend/distributed/unstructured/fem/segregated \
  tests/reference/openaccel/update_replay/algebra_check.cpp -o /tmp/mars-update-algebra
/tmp/mars-update-algebra
python3 tests/reference/openaccel/update_replay/capture_check.py \
  --source ../OpenAccel --output /tmp/mars-update-capture-check
PYTHONPATH=scripts python3 -m unittest \
  scripts/test_openaccel_update_check.py scripts/test_openaccel_public_run.py
```

Use a fresh output directory. Actual instrumented STK compilation/capture and
CUDA execution remain pending. The public channel may not exercise reversal
transitions; synthetic branch checks must not be reported as actual-flow
coverage. Global geometry/reconstruction, CSR assembly, constraints/gauge,
MPI halos and a complete SIMPLE iteration remain separate gates.

## Interactive Daint commands after review and publication

From `mars/mlir`, with the working **OpenAccel Cray-MPICH environment**. The
helper creates an instrumented reference worktree using the existing dependency
cache; it does not install a new environment. Execute the guarded group so a
failed build cannot fall through to the launch.

```bash
update_work=$(mktemp -d /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-updates-XXXXXX)
printf 'Update capture directory: %s\n' "$update_work"
(
set -euo pipefail
python3 ../scripts/prepare_openaccel_reference.py \
  --source /capstor/scratch/cscs/gandanie/git/OpenAccel \
  --public-case /capstor/scratch/cscs/gandanie/git/OpenAccel/mars-reference-inputs-20260920-v2 \
  --output "$update_work/bundle" --include-updates
python3 "$update_work/bundle/build_reference.py" \
  --source /capstor/scratch/cscs/gandanie/git/OpenAccel \
  --baseline-build /capstor/scratch/cscs/gandanie/git/OpenAccel/prgenv \
  --destination "$update_work/source" --jobs 4 2>&1 | tee "$update_work/build.log"
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --export=ALL --kill-on-bad-exit=1 \
  env OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
  python3 "$update_work/bundle/run_openaccel_public.py" \
  --executable "$update_work/source/build/openaccel-3D.exe" \
  --output "$update_work/run" --capture-interior
)
```

Require `update_capture_completed`. Node and boundary block recapture are not
needed; interior capture remains the existing harness's public-case check.
Retain the printed `update_work` path when switching shells/environments.

Then restore the configured **MARS CUDA environment**, still in `mars/mlir`:

```bash
(
set -euo pipefail
cmake -S .. -B .
cmake --build . --target mars_segregated_update_replay mars_segregated_update_algebra_check -j4
python3 ../scripts/openaccel_update_check.py "$update_work/run/exports/updates" \
  --pack "$update_work/update-inputs.txt"
./examples/distributed/unstructured/mars_segregated_update_algebra_check
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_update_replay \
  "$update_work/update-inputs.txt" 2>&1 | tee "$update_work/update-replay.log"
git rev-parse HEAD > "$update_work/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_update_replay \
  ../backend/distributed/unstructured/fem/segregated/mars_segregated_update.hpp \
  "$update_work/update-inputs.txt" > "$update_work/replay-sha256.txt"
)
```

Only inputs go to the GPU. Expected outputs stay on the host. The replay's
computed-output copy back is a test-harness operation. Each scalar is checked
with `abs(actual-expected)/max(1,abs(expected)) <= 1e-12`; unused output padding
is checked too. Passing this gate validates the captured updates, not nonlinear
convergence or a complete assembled MARS solve.
