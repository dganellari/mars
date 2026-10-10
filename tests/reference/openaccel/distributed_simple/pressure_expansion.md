# Retained pressure corrections in SIMPLE

`--pressure-expansion 1` is an opt-in extension of the distributed steady,
laminar Tet4 SIMPLE runner. It preserves a pressure value as two FP64
components, `high + low`, including the part lost when a small correction is
added to a larger value. The default remains disabled.

The user-run `simple-pressure-expansion-zj3Ndd` frozen GPU replay passed the
original pressure target with this representation. Its rounded one-component
candidate and reference-library replay still failed. That result validates one
frozen equation, not nonlinear pump convergence or OpenAccel field parity.

## Equations and acceptance

The incompressible equations remain `rho (u dot grad) u = -grad p + mu laplace u`
and `div u = 0`, with physical pressure in Pa and dynamic viscosity `mu`.
The pressure increment solves the existing assembled CVFEM equation
`A phi = b`. Geometry, Tet4/shifted Tri3 quadrature, boundary elimination,
momentum influence coefficients and pressure reference are unchanged.

The normal pressure solve runs first. If it fails acceptance, recovery solves
`A delta = b - A x` with the existing matrix and AMG hierarchy, using a
compensated defect. There are at most four correction solves, each bounded by
the original iteration cap. Their inner relative target is 0.1. Original
stopping controls are restored on exit.

The final candidate is accepted only when **both** the Hypre CSR check and the
independent original MARS owned-row check pass
`||b - A high - A low|| <= max(atol, rtol ||b||)`.
Products of the two components are accumulated separately. Evaluation and norm
rounding are bounded; RHS lower bounds make the target conservative. The shared
norm enclosure includes an additive squared-norm underflow bound under gradual
FP64 underflow. Nonfinite or uncertifiable values fail. No tolerance is relaxed
and a backend convergence flag alone cannot accept a candidate.

Pressure is updated as `p += alpha_p phi` in the paired representation. The
shifted incremental pressure gradient gathers complete owned node stars, using
`0.5 (p_right - p_left) SCS_area / volume` at both edge endpoints; its shifted
Tri3 boundary numerator is zero. The pair is retained through outlet pressure
moments and traces, compact/reconstructed pressure-gradient cancellation,
momentum RHS assembly, and mass-flux relaxation. Velocity and mass flux remain
FP64, rounded after the pressure action. Outlet reversal uses the sign after
paired pressure differences have been summed. Pressure-change diagnostics
include the low component.

The existing SIMPLE order is preserved: new pressure, predicted velocity and
old pressure gradient update fluxes; reversal follows; then the pressure
increment gradient corrects velocity. The pseudo-time diagonal and all three
relaxation parameters retain their existing meanings.

This changes arithmetic and adds linear recovery beyond OpenAccel's recorded
path. It must not be labelled identical solver behavior or matched iteration
history. Ordinary snapshot comparison rejects expanded runs. CLI combinations
with first-step audits, snapshots, failure capture or legacy refinement are
rejected. Explicit pressure rtol/atol are required. A supported
`--pressure-solver-profile` may be used with expansion.

## Device and output scope

High/low fields, node incidence and defect scratch stay on device. Ownership and
complete-star checks remain in force. Both components are exchanged together
through existing sparse peer halos. Hypre is initialized before GPU-aware MPI
is enabled. Recovery reuses the matrix and AMG setup; its temporary vectors are
allocated per recovery invocation and reused across correction rounds.

Outlet averaging reduces four host control scalars with paired addition; it
does not copy a pressure field to the host. Acceptance also reduces scalar
norms. No device-wide synchronization is added to the iteration; the Hypre
compute stream is completed before copying an accepted pair device-to-device.
This implementation has no measured performance or multi-node scaling claim.

Normal CSV output contains the normalized high pressure component only. It
cannot reconstruct the retained pair or certify the expanded linear residual,
and is not an expanded restart format. Detailed private fields are not required
for the public gates below.

## Validation and GPU command

Local tests cover cancellation, a nonzero low pressure gradient, outlet moments
and reversal, velocity correction, momentum forcing, pressure-change metrics,
and owned-row/ghost residual rejection. Three-step public-channel runs on
1/2/4 CPU ranks match the ordinary implementation for upwind and shifted
high-resolution configurations. Real CPU Hypre replay regressions exercise the
shared recovery implementation, including a rational rounding oracle.

The first user-run GPU check on 2026-10-10 passed the one-rank matrix/recovery
gate: 127 passed, zero failed, eight skipped. The following operator gate
aborted with a CUDA launch failure reported during cleanup. Source inspection
found that its single-thread fixture passed local arrays to assembly atomics.
The matrix, RHS and boundary-balance destinations now use device allocations;
test-only synchronization reports kernel errors before cleanup. Host operator
checks pass on 1/2/4 ranks. The corrected GPU operator gate and remaining
multi-rank/flow gates still need to run. Full-flow convergence and pump field
parity remain unverified.

A subsequent GPU flow comparison reached the snapshots and failed one
`blend@1` entry at 1.35e-8 against the 1e-8 field threshold. The test had used
pressure stopping rtol 1e-12 for the reference and 1e-10 for the expanded run.
Both now explicitly use rtol 1e-12 and atol zero; the field threshold is
unchanged. This removes a solver-accuracy difference, but its role in that
mismatch still needs GPU confirmation. With `--write-reference`, the expansion
flag selects the matched pressure target while retaining ordinary arithmetic.
The target is recorded in the reference header, so old references are rejected.

Run this block in the **MARS terminal only**, in one tmux pane.
It builds the production executable and runs public synthetic tests on 1/2/4
GPUs. No OpenAccel environment, Python modules or private geometry are needed.
The matrix gate must produce a nonzero retained low component and verify that
the paired residual passes while the rounded residual fails.

```bash
(
set -euo pipefail
repo=/capstor/scratch/cscs/gandanie/git/mars-v010-check
cd "$repo"
test "$(git branch --show-current)" = cstone
git fetch origin cstone
git merge --ff-only origin/cstone
cd build-hypre
cmake -S .. -B . \
  -DCMAKE_PROJECT_mars_INCLUDE="$repo/tests/reference/openaccel/distributed_matrix/inject.cmake"
cmake --build . --parallel 4 --target mars_segregated_simple \
  mars_distributed_matrix_cuda_gate mars_simple_pressure_expansion_cuda_gate \
  mars_distributed_simple_cuda_gate

run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-expansion-gate-XXXXXX)
printf 'Public results: %s\n' "$run"
git rev-parse HEAD > "$run/mars-revision.txt"
sha256sum ./mars_distributed_matrix_cuda_gate ./mars_simple_pressure_expansion_cuda_gate \
  ./mars_distributed_simple_cuda_gate ./examples/distributed/unstructured/mars_segregated_simple \
  > "$run/executable.sha256"
export MARS_HYPRE_SPMV_VENDOR=0 MARS_HYPRE_FLEXGMRES=0 MARS_HYPRE_ABSTOL=0
export MARS_HYPRE_VERBOSE=0 MARS_HYPRE_RESIDUAL_AUDIT=0
unset MARS_HYPRE_MINITER MARS_HYPRE_MAXX_RATIO MARS_HYPRE_NULLX_RATIO CUDA_LAUNCH_BLOCKING
unset MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  "$HOME/affinity/bind_numa.sh" ./mars_distributed_simple_cuda_gate \
  --mesh 4x2x2 --iterations 3 --high-resolution 1 --velocity-interpolation linear-linear \
  --pressure-expansion 1 \
  --write-reference "$run/reference.bin" 2>&1 | tee "$run/reference.log"
for np in 1 2 4; do
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    "$HOME/affinity/bind_numa.sh" bash -c '
      set -e
      "$1/mars_distributed_matrix_cuda_gate"
      "$1/mars_simple_pressure_expansion_cuda_gate"
      "$1/mars_distributed_simple_cuda_gate" --mesh 4x2x2 --iterations 3 \
        --high-resolution 1 --velocity-interpolation linear-linear \
        --pressure-expansion 1 --tolerance 1e-8 --reference "$2/reference.bin"
    ' _ "$PWD" "$run" 2>&1 | tee "$run/np$np.log"
done
printf 'Public gates passed. Results: %s\n' "$run"
)
```
