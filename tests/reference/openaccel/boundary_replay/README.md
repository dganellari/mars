# Frozen boundary assembly replay

This slice implements the pinned OpenAccel incompressible Tet4 inlet, outlet and
wall blocks. Momentum wall assembly uses a native Tri3 block. It adds a device
operator and a standalone replay target; existing flow solvers are unchanged.

The steady equations are `div(rho*u) = 0` and
`div(rho*u*u) = -grad(p) + div(mu*(grad(u)+grad(u)^T))`. Pressure is physical
pressure (Pa), `rho` is kg/m³, `mu` is dynamic viscosity (Pa s), sample areas
point outward (m²), and momentum influence components `d` have units m³ s/kg.
The unknowns in these local blocks are velocity and pressure increments.
There is no physical timestep or additional `rho/dt` in the boundary terms.

The capture hook runs immediately before the six native boundary `applyCoeff_`
calls. It records the workspaces **after boundary side values have been
substituted**, the actual LHS/RHS, connected global node IDs and face sample
maps. The pressure/velocity interpolation weights and gradients are supplied
reference inputs. On an unshifted linear triangle the nearest-node weight is
11/18 and the other two are 7/36. No face-mean substitution is made.
Moving meshes, rotating frames, compressibility, outlet body forces and harmonic
outlet-gradient blending are rejected by the capture hook.

| Stage | Captured reference routine suffix | Block |
|---|---|---|
| 0 pressure inlet | `InletSpecifiedVelocity_` | 4 × 4 |
| 1 pressure outlet | `OutletSpecifiedPressure_` | 4 × 4 |
| 2 pressure wall | `WallNoSlip_` | 4 × 4 |
| 3 momentum inlet | `InletSpecifiedVelocity_` | 12 × 12 |
| 4 momentum outlet | `OutletSpecifiedPressure_` | 12 × 12 |
| 5 momentum wall | `WallNoSlip_` | 9 × 9 |

Both assemblers live under `src/assemble/flow/segregatedFlow/` in the pinned
reference. In MARS, each stage calls `boundary_block` in
`fem/segregated/mars_segregated_boundary.hpp` (under the unstructured backend).

## Discrete terms

Let `A_s` be an outward sample area, `N_sf` its face weights and `r` its nearest
node. Inlet and wall continuity use `m_s = rho_s*Ubc_s·A_s`, scatter `-m_s` to
the nearest RHS row, and add no LHS. For the stationary public wall `Ubc=0`.
The pressure inlet's reversal flag suppresses that sample.

The outlet uses distinct interpolated LHS and RHS influence vectors:

```
g_s = 0.5*(g_nearest + g_opposite)
m_s = sum_j rho_s * (u_s,j - d_RHS,s,j*(grad(p_mix)_j-g_s,j)) * A_s,j
H[r,c] += -sum_j rho_s*d_LHS,s,j*grad(N_c)_j*A_s,j * bc_multiplier[c]
b[r]   -= m_s
```

`p_mix` contains frozen face trace values and the volume opposite-node value.
The mask is zero on the three face columns and one on the opposite column.
This is column elimination, not replacement of the continuity row. A unit-tet
check gives total flux changes `+0.6` for an opposite pressure change of `0.2`
and `-0.6` for the same common face-trace change, with `rho=1,d=2`.
Constant volume pressure has nonzero outlet action, but that alone does not
prove the assembled system is nonsingular.

Momentum inlet advection contributes `-m_old*Ubc` to the RHS. Its full viscous
traction is `-mu*(grad(u)+grad(u)^T)*A_s`, with face columns masked in the LHS
and the prescribed side velocities retained in the residual. This momentum
inlet routine does not branch on reversal flags.

Momentum outlet advection contributes `m_old` to the nearest diagonal and
`-m_old*u_nearest` to the RHS. Its viscous traction is projected by
`P = I - n*n^T`; the cross-component stress terms remain. A reversed outlet
sample contributes zero to both pressure and momentum blocks.

The wall contributes `c_s*P*N_sf` to the velocity LHS and
`-c_s*P*(sum_f N_sf*u_f-Ubc_s)` to the RHS. `c_s` is the captured wall
coefficient, in kg/s. Its construction from wall distance is outside this gate.
This tangential block alone is not the complete wall boundary condition.

## Interactive Daint commands, after review and publication

Use the working **OpenAccel Cray-MPICH environment**, from `mars/mlir`.
No dependency installation or MPI-family change is needed. The public case is
the same 425-node/1536-tet channel as the previous interior/node captures.

```bash
set -o pipefail
boundary_work=$(mktemp -d /capstor/scratch/cscs/gandanie/git/OpenAccel-reference-boundary-XXXXXX)
printf 'Boundary run directory: %s\n' "$boundary_work"
python3 ../scripts/prepare_openaccel_reference.py \
  --source /capstor/scratch/cscs/gandanie/git/OpenAccel \
  --public-case /capstor/scratch/cscs/gandanie/git/OpenAccel/mars-reference-inputs-20260920-v2 \
  --output "$boundary_work/bundle" --include-boundary
python3 "$boundary_work/bundle/build_reference.py" \
  --source /capstor/scratch/cscs/gandanie/git/OpenAccel \
  --baseline-build /capstor/scratch/cscs/gandanie/git/OpenAccel/prgenv \
  --destination "$boundary_work/source" --jobs 4 2>&1 | tee "$boundary_work/build.log"
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --cpus-per-task=1 --cpu-bind=cores --export=ALL --kill-on-bad-exit=1 \
  env OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
  python3 "$boundary_work/bundle/run_openaccel_public.py" \
  --executable "$boundary_work/source/build/openaccel-3D.exe" \
  --output "$boundary_work/run" --capture-interior
```

Execute each command after its predecessor succeeds. The reference helper uses
the existing baseline cache to create/build a separate instrumented worktree.
It must report `boundary_capture_completed`. The capture includes interior
records because the established reference harness checks the public mesh and
two-iteration execution through those records. Node recapture is not required.

Then use the configured **MARS CUDA environment**, still from `mars/mlir`.
If switching shells, set `boundary_work` to the exact directory printed above.

```bash
set -o pipefail
cmake -S .. -B .
cmake --build . --target mars_segregated_boundary_replay mars_segregated_boundary_algebra_check -j4
python3 ../scripts/openaccel_boundary_check.py "$boundary_work/run/exports/boundary" \
  --pack "$boundary_work/boundary-inputs.txt"
./examples/distributed/unstructured/mars_segregated_boundary_algebra_check
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_boundary_replay \
  "$boundary_work/boundary-inputs.txt" 2>&1 | tee "$boundary_work/boundary-replay.log"
printf '%s\n' "${PIPESTATUS[0]}" > "$boundary_work/boundary-replay.exit"
git rev-parse HEAD > "$boundary_work/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_boundary_replay \
  ../backend/distributed/unstructured/fem/segregated/mars_segregated_boundary.hpp \
  ../examples/distributed/unstructured/mars_segregated_boundary_replay.cu \
  "$boundary_work/boundary-inputs.txt" > "$boundary_work/replay-sha256.txt"
```

This is a one-GPU local-block replay, not distributed assembly. Only inputs go
to the GPU; expected outputs stay on the host. The test harness copies computed
blocks back for comparison. Each LHS and RHS uses scaled error
`abs(actual-expected)/max(1,max(abs(expected))) <= 1e-12`; zero padding is also
checked. Nonfinite values, bad maps, truncated data and missing stages fail.
The Python checker verifies source pins, unique owned faces, complete rank/call
files, and pressure/momentum face identity and outward areas by global node ID.
The replay reports reversed-sample counts per stage; zero counts mean the
captured CUDA run has not exercised that branch.

## Local evidence and remaining scope

727 independent host algebra checks pass: affine pressure, compact derivative,
distinct influence vectors, nonuniform interpolation, reversed samples,
momentum finite-difference derivatives and tangential stress.
92,400 comparisons pass against arithmetic loops extracted from the pinned
reference; maximum scaled difference is `2.08e-15`. The source hashes are checked
before extraction. This exercises 1,200 synthetic blocks, not an STK run.

The exact capture expressions compile with mock mesh plumbing. Their writer,
parser and production host replay pass 1,872 scalar comparisons across all six
stages and two calls. Separate hand-worked transport fixtures pass 936 checks.
Transport rejection tests and public runner failure-path tests are also included.

To reproduce the source-expression checks locally (with an MPI C++ compiler):

```bash
python3 tests/reference/openaccel/boundary_replay/extract_oracle.py \
  --source ../OpenAccel --output /tmp/mars_boundary_reference.hpp
c++ -std=c++20 -Wall -Wextra -Werror -Wno-unused-variable \
  -Ibackend/distributed/unstructured/fem/segregated -I/tmp \
  tests/reference/openaccel/boundary_replay/native_check.cpp -o /tmp/mars-boundary-native
/tmp/mars-boundary-native
python3 tests/reference/openaccel/boundary_replay/capture_check.py \
  --source ../OpenAccel --output /tmp/mars-boundary-capture-check
```

The capture test requires a fresh output directory. Actual instrumented STK
compilation/capture and CUDA execution remain pending. Reversal flags are frozen
inputs: this gate does not validate their update rule or transition history.
Also pending: trace updates, wall-distance coefficients, normal inlet direction,
corner boundary ordering, geometry/reconstruction, global constraints/gauges,
CSR/halo assembly, a complete SIMPLE iteration, convergence and real-case parity.
