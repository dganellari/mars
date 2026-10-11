# First-iteration assembly and solve audit

The fresh startup pair agrees at state zero but differs after iteration one.
That result does not distinguish different equations from different linear-solver
accuracy. This opt-in audit captures one iteration, preserves numerical controls
and solver settings, and compares each stage by source node identity.

All matrices, fields, logs, identifiers and detailed errors are **private**.
Only `comparison-public.json` may be shared. The numerical solver stays on the
GPU; the explicit diagnostic option copies arrays to the host for file output.
Ordinary runs perform none of these copies. Audit timings are not benchmarks.

## Discrete quantities and interpretation

For steady, incompressible, laminar Tet4 CVFEM, the momentum system solves a
velocity increment: `A_u delta_u = b_u`, `u* = u0 + delta_u`. The pressure system
solves `A_p phi = b_p`, followed by `p1 = p0 + alpha_p phi` and
`u1 = u* - d grad(phi)`. Pressure and phi are physical pressure in Pa; velocity
is m/s. `d` has units m³ s/kg. These are the actual segregated updates, not a
replacement pressure Laplacian or a different boundary treatment.

At public OpenAccel revision `0d69041ba1afda63e9e4328d9e0d9834bba37756`:

- `linearSystem.h::writeLinearSystem_` writes scalar CSR and RHS before the
  backend solve, after optional normalization/diagonal row scaling.
- `meshNodeGraph.cpp` publishes full-graph solver node numbers in `aux`.
  The supported one-zone case uses that full graph. Missing numbering or a
  numbering that is not a complete permutation rejects the audit.
- `pressureCorrectionEquation.h::correctField_` stores unrelaxed
  `pressure_correction` and its gradient. `segregatedFlowEquations.cpp` then
  subtracts `du * pressure_correction_gradient` from velocity. Thus the reference
  predictor is reconstructed from the final velocity plus that correction.
  This reconstruction relies on the inspected supported control path; it is not
  a separately instrumented reference predictor.

The compatible PETSc-enabled reference executable is reused. Preparation enables
its existing `write_system` and field output, and stops after one iteration.
It does not replace a solver family, change tolerances or disable scaling.
MARS's `--first-step-audit 1` records matrices/RHS and intermediate fields through
the existing runner observer. It requires one iteration and distributed state
snapshots. Per-rank binary files contain the graph, source IDs, owned nodes,
momentum system/increment/predictor/influence, pressure system/increment and
correction gradient. No global gather is added to MARS. Ghost rows can be partial
or poisoned; only owned rows participate in the audit.

The comparator verifies launch/input/output hashes, initial fields, node coverage,
update identities and matrix dimensions. It canonicalizes both matrices by
source-node/component indices. Each row and its RHS are divided by that row's
maximum absolute coefficient before comparison, so positive row scaling does
not produce a false assembly mismatch. Matrix tolerance is `1e-10`; RHS errors
use the corresponding unknown scale as well. Stored zero entries do not matter.
Field errors use `U`, `rho U²`, `rho U²/L` and `L/(rho U)` as appropriate, with a
`1e-5` threshold. No pressure offset is removed.

Both saved linear solutions are checked against their own captured equations.
Public flags report whether `||b-Ax||₂/||b||₂ <= 1e-8` (null if RHS is zero) and
whether the maximum row backward error is at most `1e-8`. The latter denominator
is `|b_i| + sum_j |A_ij x_j| + max_j |A_ij| * unknown_scale`.
These are common diagnostic thresholds, not the backends' configured stopping
tests, and small residuals do not bound field error without conditioning evidence.
Reference residuals use the reference's possibly scaled system.

The saved local MARS increments are also checked before replacing ghost values
with globally owned values. `mars_referenced_copies_equal_owners` requires exact
agreement for every entry referenced by an owned row; binary halo exchange should
not round these values. A disagreement separates a publication/mapping problem
from a residual computed using a different copy of the solution.

`pressure_solve_checks` compares recomputed residuals against the saved controls.
Explicit MARS pressure options use `max(atol, rtol*||b||)`. Without them, its
independent acceptance check uses `1e-13 + 1e-10*||b||`, separately from the Hypre
stopping request (`rtol=1e-12`, `MARS_HYPRE_ABSTOL` or zero). The report distinguishes
these limits from the common `1e-8` diagnostic. The reference's declared limit
uses its own saved, possibly row-scaled RHS; its backend stopping norm and
convergence reason are not verified. In particular PETSc may use a preconditioned
norm, and the inspected OpenAccel wrapper does not check `KSPGetConvergedReason`.
CPU recomputation can also differ from GPU reductions by roundoff. These flags
do not by themselves prove a backend defect or an acceptable field error.

`pressure_common_system_checks` additionally evaluates both increments on the
same reference matrix/RHS and original reference target, and checks the reverse
application on MARS rows. It distinguishes unscaled matrices from positive
row-scaled equivalence. These are FP64 residual diagnostics without an evaluation
error bound; the field verdict and solver acceptance are unchanged. See the
[saved-result recipe](pressure_profile.md#recheck-the-saved-original-pressure-result)
for interpretation and a comparison that launches neither solver.

`first_differing_stage` gives the first captured mismatch in execution order.
A pressure-stage difference can follow a differing momentum predictor; it is not
automatically a pressure-kernel bug. A reconstructed predictor mismatch may also
require checking reference output semantics. No root cause or pump parity is
claimed until this experiment is run.

## OpenAccel terminal

Use the working OpenAccel runtime and Python environment. This creates a new
pair from the already successful startup pair and leaves that evidence intact.
No OpenAccel rebuild is needed for the inspected built-in output paths.

```bash
(
set -euo pipefail
root=/capstor/scratch/cscs/gandanie/git
scratch=/capstor/scratch/cscs/gandanie
git -C "$root/mars-v010-check" pull --ff-only
python3 -c 'import numpy, netCDF4, yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
previous="$scratch/simple-startup-y0Be1K/pair"
pair=$(mktemp -d "$scratch/simple-first-step-XXXXXX")/pair
probe="$root/mars-v010-check/scripts/simple_startup_probe.py"
exe="$root/OpenAccel/prgenv/openaccel-stk.exe"
python3 "$probe" prepare --first-step-audit \
  --case "$previous/case.json" --reference-dir "$previous/reference" --output-dir "$pair"
printf '%s\n' "$pair" > "$scratch/simple-first-step-current.txt"
printf 'Private pair: %s\n' "$pair"
python3 "$probe" run --pair "$pair" --solver openaccel --ranks 4 --executable "$exe" -- \
  srun --account=csstaff --time=00:15:00 --nodes=1 --ntasks-per-node=4 \
  --cpus-per-task=1 --cpu-bind=cores --export=ALL --kill-on-bad-exit=1
)
```

## MARS terminal, after reference capture

Use the MARS CUDA/Hypre runtime and Python environment. Build the driver and the
public output test; run the four-rank output test before the private capture.
The private solver normally exits 2 at the single-iteration limit. Its launcher
uses `--kill-on-bad-exit=0` so peers can finish writing files. No raw log is printed.

```bash
(
set -euo pipefail
root=/capstor/scratch/cscs/gandanie/git
scratch=/capstor/scratch/cscs/gandanie
cd "$root/mars-v010-check/build-hypre"
git pull --ff-only
python3 -c 'import numpy, netCDF4, yaml'
cmake -S .. -B . -DMARS_ENABLE_SEGREGATED=ON
cmake --build . --parallel 4 --target mars_segregated_simple mars_simple_output_profile_cuda_gate
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
unset CUDA_LAUNCH_BLOCKING MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
export MARS_HYPRE_SPMV_VENDOR=0 MARS_HYPRE_FLEXGMRES=0 MARS_HYPRE_ABSTOL=0
export MARS_HYPRE_VERBOSE=0 MARS_HYPRE_RESIDUAL_AUDIT=0 MARS_SIMPLE_PRESSURE_AUDIT=1
gate=$(mktemp -d "$scratch/simple-first-step-io-XXXXXX")
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=4 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  "$HOME/affinity/bind_numa.sh" ./mars_simple_output_profile_cuda_gate "$gate"
pair=$(cat "$scratch/simple-first-step-current.txt")
test -f "$pair/reference/launch.json"
probe=../scripts/simple_startup_probe.py
python3 "$probe" run --pair "$pair" --solver mars --ranks 4 \
  --executable "$PWD/examples/distributed/unstructured/mars_segregated_simple" -- \
  srun --account=csstaff --time=00:15:00 --nodes=1 --ntasks-per-node=4 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=0 \
  "$HOME/affinity/bind_numa.sh"
status=0
python3 "$probe" compare --pair "$pair" --output "$pair/comparison-public.json" || status=$?
cat "$pair/comparison-public.json"
printf 'Share only: %s\n' "$pair/comparison-public.json"
exit "$status"
)
```

## Recheck an existing capture without launching either solver

The first MPQGzD audit matches both assembled systems and the momentum predictor,
then differs at the pressure increment. Both pressure residuals exceed the common
diagnostic threshold. This localizes the discrepancy but does not identify its
cause. Use this saved-file check to distinguish declared targets and local versus
owned pressure values. All input and launch hashes are checked again; earlier
captures and reports are preserved. No C++ build or GPU allocation is needed.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
git -C "$repo" pull --ff-only
python3 -c 'import numpy, netCDF4, yaml'
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
report=$(mktemp -d "$scratch/simple-pressure-audit-XXXXXX")
status=0
python3 "$repo/scripts/simple_startup_probe.py" compare \
  --pair "$scratch/simple-first-step-MPQGzD/pair" \
  --detail-dir "$report/private" --output "$report/public.json" || status=$?
cat "$report/public.json"
printf 'Share only: %s\n' "$report/public.json"
exit "$status"
)
```

## Local validation

The saved-file target check now reports both pressure solutions inside their
declared tolerance, with matching MARS owner/ghost copies. The next controlled
experiment is the [pressure-only accuracy comparison](pressure_accuracy.md).
It uses fresh captures and tighter pressure targets; baseline evidence is retained.

Synthetic paired evidence covers reordered global/local nodes, overlapping ghost
copies, poisoned ghost rows, positive row scaling, 32/64-bit reference indices,
independently wrong momentum/pressure matrices and RHS, bad solutions, missing or
truncated data, nonfinite owned coefficients, duplicate ownership, missing fields,
bad solver numbering, inconsistent updates and changed launch artifacts.
Further tests separate local/owned residuals with deliberately stale ghosts,
loose declared targets from the common threshold, and absolute tolerance floors;
reanalysis preserves earlier reports and still rejects tampered captures.
The C++ output gate exercises the actual binary writer and overwrite rejection
on one and four host MPI ranks with address/undefined-behavior sanitizers.
CUDA compilation and the actual paired capture remain Alps validation.

```bash
PYTHONPATH=scripts python3 -m unittest test_simple_first_step_audit test_simple_startup_probe test_simple_snapshot_compare test_prepare_simple_deck
```
