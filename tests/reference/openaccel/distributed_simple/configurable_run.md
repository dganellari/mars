# Configurable distributed SIMPLE

`mars_segregated_simple` now uses the native distributed path on one or more
ranks. It reads Exodus with the MARS C++ reader, restores exact coordinates
through ElementDomain's SFC keys, and uses complete owned element stars.
Setup geometry, assembly, Hypre solves, halo fields and numerical reductions
stay on the GPU. File I/O, launch/control APIs and scalar reports use the CPU.
The old `mars_segregated_simple_mpi` gate target remains a compatibility entry
point to the same source. Prepared-text inputs remain in the reference gates;
the production executable accepts native Exodus only.

The equations are steady, incompressible, laminar Navier-Stokes with upwind
advection by default (or opt-in [high-resolution advection](high_resolution.md)),
physical pressure in Pa and dynamic viscosity in Pa s. `--mu` sets
dynamic viscosity; kinematic viscosity is mu/rho. For water, rho=1000 and
mu=.001 give nu=1e-6. This specifies material properties, not a turbulence
model or a convergence guarantee at high Reynolds number.

On inlet faces, velocity is `-U*A/|A|` for outward area vector A. Walls are
no-slip. `--outlet-pressure` is the area-mean pressure target in the existing
outlet trace update, including its close/reopen treatment. Mesh coordinates
must be in metres. `--reference-length` changes residual normalization only.
Pseudo-time and relaxation are steady iteration controls, not physical time.
The existing defaults reproduce the public channel configuration.

Current scope is one 3D Tet4 element block, one inlet set, one outlet set and
one or more explicitly selected wall sets. Every exterior face must be tagged
exactly once; unknown, duplicate or missing selections are rejected. Rank zero
reads the file; element partitions and coordinate/tag requests use bounded GPU
buffers. Root retains the source mesh during setup; this is not parallel I/O.
Empty owned-row ranks and periodic/multi-block configurations remain unsupported here.

Local validation on 2026-09-27 used a strict C++20 build and 62 CPU/MPI tests.
The first pass found two oblique four-rank fixture failures: the test's slab
partition used rotated coordinates and left a rank empty. Preserving the logical
partition fixed that test setup; all eight configured-case tests then passed,
including reversal on 1/2/4 ranks. The other 54 tests passed unchanged. The
configured single-rank CPU reference converged at iteration 2050. Another 1644
SIMPLE algebra checks and 12 field-comparator tests passed; the 133 control and
inlet-normal checks also passed AddressSanitizer and UndefinedBehaviorSanitizer.
The subsequent Daint run at `simple-configured-SQbfsO` passed the default
configuration on 1/2/4 GPUs (1277 iterations and field parity), then failed in
the oblique single-rank pressure solve. After the repair below, user-reported
2/4-rank oblique runs converged at 2050 iterations and matched the one-rank
field file within 2.006e-15 in velocity/U and 1.277e-13 in pressure/(rho U^2).
The supplied excerpt did not include the one-rank convergence line or binary
hash. The subsequent default four-rank regression converged at 1277 and matched
its saved baseline within 5.551e-13 in scaled pressure. These are reported CUDA
results for those fixtures, not a general material/geometry validation.

## Small-system AMG repair

The failed run printed only iteration 0 because reports were 100 iterations
apart. Its RHS maximum, .0051740964707256282, matches pressure correction at
iteration 11: a local replay from the native Exodus reader gave
.0051740964707257947. This identifies the likely stage without relying on
the wrapper's former, unsupported null-mode diagnosis.

Using CPU Hypre 3.1.0 and the actual frozen 81-row pressure matrix reproduces
`Solve=1`, three iterations, relative residual 1 and a zero solution. Setup
produces one AMG level with no l1 row norms. Hypre's single-level cycle uses
the requested l1-Jacobi smoother, while setup prepared the default coarse
type 9. The missing l1 data makes that preconditioner action zero. Setting
the coarse relaxation explicitly to 18 prepares the row norms and restores
the solve. Removing the GMRES minimum iteration count alone does not fix it.
The relevant upstream paths are
[single-level cycling](https://github.com/hypre-space/hypre/blob/v3.1.0/src/parcsr_ls/par_cycle.c),
[row-norm setup](https://github.com/hypre-space/hypre/blob/v3.1.0/src/parcsr_ls/par_amg_setup.c)
and [l1-Jacobi](https://github.com/hypre-space/hypre/blob/v3.1.0/src/parcsr_ls/par_relax.c).
The installed Daint Hypre version has not been retrieved.

Both SIMPLE runners now explicitly select coarse l1-Jacobi (18), which also
avoids the host Gaussian coarse solve. Other users of `HypreGMRESSolver`
retain their existing defaults. This changes only the preconditioner in
the momentum and pressure solves, not their matrices, RHS, boundary terms,
relaxation or acceptance tolerances. The independent true-residual check
remains mandatory. Failures now identify the stage and SIMPLE iteration.

Local frozen-system checks with the repaired settings passed for momentum
and pressure at iterations 1, 11 and 101; the largest true relative residual
was 3.96e-13. Nine configured CPU/MPI regression tests passed. These tests
do not validate device execution. The existing CUDA matrix gate also tests
81-node scalar and three-component systems with the new coarse setting.
Use the native oblique run below as the end-to-end regression; run it first,
then repeat the default-channel comparisons because AMG has changed.

From the configured MARS CUDA build, after pulling and building as below:

```bash
(
set -euo pipefail
simple_run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-amg-fixed-XXXXXX)
printf 'Results: %s\n' "$simple_run"
git rev-parse HEAD > "$simple_run/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_simple > "$simple_run/sha256.txt"
for np in 1 2 4; do
  srun --account=csstaff --time=00:20:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
    --mesh ../tests/data/public_simple_oblique/channel.exo \
    --inlet-ss feed --outlet-ss exit --wall-ss casing,cover \
    --rho 2 --mu .4 --inlet-velocity .2 --outlet-pressure .03 \
    --pseudo-dt .005 --relax-u .4 --relax-p .2 --relax-mass .8 --outlet-beta .1 \
    --reference-length 1 --iterations 3000 --report-every 100 \
    --output-prefix "$simple_run/oblique-$np" \
    2>&1 | tee "$simple_run/oblique-$np.log"
  python3 ../tests/reference/openaccel/distributed_simple/compare_fields.py \
    "$simple_run/oblique-1-fields.csv" "$simple_run/oblique-$np-fields.csv" \
    --rho 2 --inlet-velocity .2
done
)
```

## Daint build and default regression

From the existing `mars-v010-check/build-hypre` MARS CUDA environment. The normal
executable does not require a CMake injection or an OpenAccel rebuild:

```bash
git pull --ff-only &&
cmake -S .. -B . &&
cmake --build . --target mars_segregated_simple --parallel 4
```

First reproduce the recorded public channel on 1/2/4 ranks with unchanged
settings. All generated artifacts go under capstor:

```bash
(
set -euo pipefail
unset MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
simple_run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-configured-XXXXXX)
printf 'Results: %s\n' "$simple_run"
git rev-parse HEAD > "$simple_run/mars-revision.txt"
sha256sum ./examples/distributed/unstructured/mars_segregated_simple \
  ../tests/data/public_outlet_channel/outlet_channel.exo \
  ../tests/data/public_simple_oblique/channel.exo > "$simple_run/sha256.txt"
for np in 1 2 4; do
  srun --account=csstaff --time=00:20:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
    --mesh ../tests/data/public_outlet_channel/outlet_channel.exo \
    --output-prefix "$simple_run/default-$np" --iterations 2000 --report-every 100 \
    2>&1 | tee "$simple_run/default-$np.log"
  python3 ../tests/reference/openaccel/distributed_simple/compare_fields.py \
    /capstor/scratch/cscs/gandanie/git/mars/mlir/simple-native-bMEZdp/channel-fields.csv \
    "$simple_run/default-$np-fields.csv"
done

for np in 1 2 4; do
  srun --account=csstaff --time=00:20:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
    --mesh ../tests/data/public_simple_oblique/channel.exo \
    --inlet-ss feed --outlet-ss exit --wall-ss casing,cover \
    --rho 2 --mu .4 --inlet-velocity .2 --outlet-pressure .03 \
    --pseudo-dt .005 --relax-u .4 --relax-p .2 --relax-mass .8 --outlet-beta .1 \
    --reference-length 1 --iterations 3000 --report-every 100 \
    --output-prefix "$simple_run/oblique-$np" \
    2>&1 | tee "$simple_run/oblique-$np.log"
  python3 ../tests/reference/openaccel/distributed_simple/compare_fields.py \
    "$simple_run/oblique-1-fields.csv" "$simple_run/oblique-$np-fields.csv" \
    --rho 2 --inlet-velocity .2
done
)
```

Require `CONVERGED` and field comparison `PASS`, including absolute pressure;
exit zero alone does not establish parity. Tolerances remain 1e-6. The new
configuration does not yet have CUDA execution evidence; the saved default
channel baseline does. A local CPU/MPI oracle tests nondefault controls and
outlet reversal independently of these end-to-end GPU runs.
