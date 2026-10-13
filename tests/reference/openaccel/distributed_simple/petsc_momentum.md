# PETSc momentum with Hypre pressure

The saved reference declares PETSc `ksp_type: fgmres`, `pc_type: bjacobi`
for momentum, with no other PETSc options. SIMPLE can now use that backend
through `MARS_ENABLE_PETSC=ON` and `--momentum-backend petsc`. Hypre remains
the default, and pressure keeps its separately verified Hypre profile.

This adds a linear backend to the existing momentum equation `A du = b`.
The nodal Tet4 CVFEM velocity increment `du` has units m/s; matrix coefficients
have units kg/s and `b` has units N. Assembly, boundary terms, shifted quadrature,
pseudo-time, relaxation and the subsequent pressure correction are unchanged.
The original rtol, atol and maximum iterations follow OpenAccel's
`src/solver/PETSc/ContextPETSc.h` at reference revision
`0d69041ba1afda63e9e4328d9e0d9834bba37756`: context tolerances override options,
and its generic minimum-iteration setting is not forwarded to PETSc.

Only this two-option, unscaled profile is supported. Preparation rejects other
PETSc options rather than interpreting them as Hypre options. PETSc supplies
the omitted defaults: FGMRES with right preconditioning and an unpreconditioned
residual norm; one block per rank, with local PREONLY and ILU(0).
These choices are checked after setup. The actual restart, divergence tolerance,
version and local solver settings are recorded privately for every selected step
and rank. PETSc option files and ambient options do not override this profile.

GPU matrices and vectors use AIJCUSPARSE and CUDA storage; the local ILU uses
cuSPARSE. Owned scalar rows are converted to COO indices on device, and numeric
values are refreshed without copying MARS matrices or vectors to host. PETSc
3.24's COO setup internally stages index arrays on host once per graph build;
this API setup cost is retained, not claimed to be GPU-only. Library
boundaries synchronize their work before memory is consumed by another backend.
The existing independent MARS check still exchanges ghosts and reduces only
owned rows. Acceptance requires both a positive PETSc convergence reason and
the explicit original-system residual to satisfy `max(atol, rtol ||b||)`.
The unchanged MARS check must also pass. Recovery is disabled in this experiment.

This aligns the solver family and configured algorithm, not the entire numerical
trajectory. CPU/GPU factorization implementations, PETSc versions, row ordering
and MPI partitions can differ. In particular, block Jacobi depends on the rows
owned by each rank. Matching rank counts does not prove matching blocks. Neither
solver settings nor completion establishes pump field parity or convergence.

## Local validation and limits

The real PETSc 3.24.0 CPU build passes synthetic tests on 1/2/4 ranks: coupled
three-component rows, uneven ownership, ghost completion, three matrix/RHS scales,
numeric refresh with and without caching, zero RHS, and rejection of capped or
inconsistent solves. A two-step public channel additionally checks the mixed
SIMPLE runner against its independent host solver. Ambient PETSc options are
cleared in the same way as the production driver. The Python preparation/comparison tests check library
identity, preservation of targets, output coverage, settings tampering and
unsupported options. These do not validate CUDA execution.

The GPU gate below uses the same public synthetic test with device buffers.
It checks actual CUDA vector, matrix and cuSPARSE factor backends. All three
rank counts must pass before the private one-step run. The comparison reuses
the original OpenAccel first-step capture; no OpenAccel run is needed here.

## Run once in the MARS terminal

Keep the MARS uenv and working Python environment active. Do not use the
OpenAccel terminal's CPU PETSc or MPI libraries. The helper builds pinned PETSc
3.24.0 from its SHA-256-verified official archive with the current MPI wrappers,
CUDA and downloaded BLAS/LAPACK. It requires compiler and download access but
does not launch a cluster job. The new build directory is retained for reuse.
The CUDA dependency build is not validated by the local CPU build.

Run the following block in just one pane. Public gate logs may be shared;
private solver logs, profiles and fields remain local. The existing capture
checks permit adding PETSc libraries but reject any replacement of the recorded
MARS MPI, Hypre or other dependencies before launching the private case.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
umask 077
mkdir -p "$scratch/tmp"
export TMPDIR="$scratch/tmp" PYTHONDONTWRITEBYTECODE=1
python3 -c 'import numpy, netCDF4, yaml; print("Python dependencies OK")'
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
run=$(mktemp -d "$scratch/simple-petsc-momentum-XXXXXX")
printf 'Results: %s\n' "$run"
bash "$repo/scripts/build_simple_petsc.sh" "$run/dependency"
petsc_root="$run/dependency/petsc-3.24.0"
cmake -S "$repo" -B "$repo/build-hypre" \
  -DMARS_ENABLE_CUDA=ON -DMARS_ENABLE_MPI=ON -DMARS_ENABLE_TESTS=ON \
  -DMARS_ENABLE_UNSTRUCTURED=ON -DMARS_ENABLE_HYPRE=ON \
  -DMARS_ENABLE_FEM_EXAMPLES=ON -DMARS_ENABLE_SEGREGATED=ON \
  -DCMAKE_CUDA_ARCHITECTURES=90 \
  -DMARS_ENABLE_PETSC=ON -DMARS_PETSC_ROOT="$petsc_root" \
  -DMARS_PETSC_ARCH=arch-mars-cuda \
  -DCMAKE_PROJECT_mars_INCLUDE="$repo/tests/reference/openaccel/distributed_matrix/inject.cmake"
cmake --build "$repo/build-hypre" --parallel 4 --target \
  mars_segregated_simple mars_simple_petsc_momentum_gate
cd "$repo/build-hypre"
mkdir "$run/public-gates"
git rev-parse HEAD > "$run/public-gates/mars-revision.txt"
sha256sum ./mars_simple_petsc_momentum_gate \
  ./examples/distributed/unstructured/mars_segregated_simple > "$run/public-gates/executables.sha256"
for np in 1 2 4; do
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    "$HOME/affinity/bind_numa.sh" ./mars_simple_petsc_momentum_gate \
    2>&1 | tee "$run/public-gates/petsc-$np.log"
done
prepare_status=0
python3 "$repo/scripts/simple_pressure_profile.py" prepare \
  --baseline-pair "$scratch/simple-first-step-MPQGzD/pair" \
  --defaults-probe "$scratch/simple-hypre-defaults-PrkNXY/probe" \
  --executable "$repo/build-hypre/examples/distributed/unstructured/mars_segregated_simple" \
  --output-dir "$run/pair" --both-linear-solvers \
  --output "$run/preparation-public.json" || prepare_status=$?
cat "$run/preparation-public.json"
if (( prepare_status != 0 )); then exit "$prepare_status"; fi
run_status=0
python3 "$repo/scripts/simple_pressure_profile.py" run \
  --pair "$run/pair" --output "$run/launch-public.json" || run_status=$?
cat "$run/launch-public.json"
if (( run_status != 0 )); then exit "$run_status"; fi
python3 "$repo/scripts/simple_pressure_profile.py" compare \
  --pair "$run/pair" --detail-dir "$run/private-comparison" \
  --output "$run/comparison-public.json"
cat "$run/comparison-public.json"
printf 'Share only: %s\n' "$run/comparison-public.json"
)
```

`simple_pressure_profile.py run` executes the saved launcher once. Do not wrap
it in another `srun`. A successful report should verify the recorded momentum
and pressure settings and preserve both original targets. Its first differing
stage and field checks determine whether this alignment changed agreement;
success is not presumed.

For local CPU development, configure the standalone `distributed_matrix`
test directory with `MARS_ENABLE_PETSC=ON`, `MARS_PETSC_ROOT` and (for an
uninstalled build) `MARS_PETSC_ARCH`. Build `mars_simple_petsc_momentum_gate`
and run `ctest -R marsPetscMomentumHost`. Production SIMPLE requires CUDA PETSc.

## Preparation fails after the GPU gate passes

`private_preflight_or_capture` alone does not identify the failed check. Retry
preparation with `--output` to write a public JSON report. It preserves all input
and library checks and exports only fixed error categories and file/line locations
in the public Python source. Exception messages, private paths, values and frame
locals are not exported. Missing reports are errors; existing reports are never
overwritten. A failed attempt may leave a partial private pair, so use a new
directory instead of editing or deleting the previous evidence.

This step needs the MARS environment and existing executable, but no rebuild,
GPU allocation or repeat of passing GPU gates. Kokkos's `OMP_PROC_BIND` warning
in a completed gate does not identify a Python preparation failure. Keep the
recorded solver environment intact while locating the failed check.

```bash
(
set -euo pipefail
scratch=/capstor/scratch/cscs/gandanie
repo="$scratch/git/mars-v010-check"
test "$(git -C "$repo" branch --show-current)" = cstone
git -C "$repo" fetch origin cstone
git -C "$repo" merge --ff-only refs/remotes/origin/cstone
python3 -c 'import numpy, netCDF4, yaml'
umask 077
run=$(mktemp -d "$scratch/simple-petsc-preparation-XXXXXX")
status=0
python3 "$repo/scripts/simple_pressure_profile.py" prepare \
  --baseline-pair "$scratch/simple-first-step-MPQGzD/pair" \
  --defaults-probe "$scratch/simple-hypre-defaults-PrkNXY/probe" \
  --executable "$repo/build-hypre/examples/distributed/unstructured/mars_segregated_simple" \
  --output-dir "$run/pair" --both-linear-solvers \
  --output "$run/public.json" || status=$?
cat "$run/public.json"
printf 'Share only: %s\n' "$run/public.json"
exit "$status"
)
```

`preparation_complete` verifies preparation only. It does not claim that the flow
solver ran or that fields agree. If preparation succeeds, use the new pair with
the `run` and `compare` commands above.
