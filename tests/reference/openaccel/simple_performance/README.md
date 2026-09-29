# SIMPLE performance validation

This work preserves the existing steady, laminar Tet4 CVFEM equations, quadrature,
boundary closure, physical pressure, relaxation and high-resolution limiter.
Momentum and pressure matrices still use current coefficients on every iteration.
Changes concern storage, communication scheduling, diagnostics, file distribution
and GPU SpMV backend selection.

## Changes and limits

- Hypre's native GPU SpMV is the default; `MARS_HYPRE_SPMV_VENDOR=1` opts back
  into its vendor path. This avoids the observed residual mismatch without
  changing either linear acceptance check. The exact vendor-path cause remains
  unresolved; see the [backend audit](HYPRE_HOST_TEST.md#gpu-spmv-backend-probe).
- `--linear-cache 1` (default) retains the Hypre wrapper's CSR packing, outer
  IJ/ParCSR objects, vectors, GMRES and AMG handles. Every solve overwrites all
  numerical entries, including zeros, and reruns setup. It does **not** freeze
  the pressure matrix or lag the AMG hierarchy. Hypre 3.1's GPU IJ refresh still
  reconstructs its internal diagonal/off-diagonal CSR arrays. Set the option to
  `0` for the former fresh-wrapper behavior. See [the lifecycle gate](HYPRE_HOST_TEST.md).
- `--halo-overlap 1` (default) packs owned fields, posts sparse GPU-buffer MPI,
  launches momentum contributions whose complete parent elements are owned,
  then waits, unpacks and assembles contributions depending on ghosts. Storage
  and MPI requests persist. Default-stream ordering protects all dependencies.
  Whether useful overlap occurs depends on the MPI transport and available
  interior work; the option alone establishes no speedup. Set `0` for the
  synchronous assembly order.
- Diagnostic SUM and MAX reductions are posted together using device buffers.
  Dependent outlet/solver reductions and coordinated error checks remain. This
  does not remove Hypre's communication or establish GPUDirect transport.
- Native C++ Exodus input is read on rank zero. Initial element partitions and
  exact-coordinate/side-tag requests move through bounded GPU buffers, rather
  than replicating the entire source mesh on every GPU. Root retains O(global
  mesh) storage during setup and services requests sequentially. This is not
  parallel file I/O or an unbounded-scale setup algorithm.
- `--field-output distributed` writes each rank's owned rows plus a root
  `PREFIX-fields.json` manifest. The comparator accepts that manifest or the
  existing CSV. `gathered` remains the compatible default; `none` skips fields.
  The native reader remains C++; Python below only creates public test inputs
  or compares saved output.

All production mesh/field computation and bulk MPI payloads remain on the GPU.
Host work consists of file I/O, API orchestration, peer/count metadata and scalar
reports. CPU code in the gates supplies independent reference calculations.
No new turbulence model or arbitrary-pump accuracy claim follows from this work.

## Timing contract

Every run prints rank-maximum setup, iteration-loop and field-output wall times.
Iteration time includes initial diagnostics, all iterations and metric-file I/O;
the printed per-iteration quotient includes startup/warmup. These three maxima
can come from different ranks and are not additive global elapsed time.

`--profile 1 --profile-warmup 10` additionally prints seven phase totals, Hypre
preparation/packing/setup/solve/finish, and application halo pack/readiness/wait/
unpack totals. It also records processor names and GPU models. Iterations 0–10
are excluded from phase/halo totals; `samples` gives the denominator. CUDA event
intervals include stream idle gaps and MPI waits, **not SM busy time**. Profiling
adds event collection fences; repeat wall-time runs with profiling disabled.

Hypre subphases are wall API times without added synchronization; packing is a
subset of preparation. Halo counts exclude startup metadata and Hypre's own
messages. Halo times are rank maxima, bytes and rounds are rank sums. Graph,
numeric-update and setup counters are lifetime rank maxima, including warmup.
A stable graph should report one wrapper graph build per equation with caching,
and one numerical refresh/setup per solve. That does not count Hypre's internal
GPU structural allocations.

## Local evidence

On 2026-09-28, a strict C++20 `-Wall -Wextra -Werror` host build passed all 95
CPU/MPI gates, including synchronous/overlapped comparisons, 1/2/4 ranks,
subcommunicators, invalid halo lifetimes, exact source routing and distributed
output. The field comparator's 13 cases passed. The production Hypre wrapper,
with CUDA primitives emulated on the host and real sequential Hypre 3.1, passed
ASan/UBSan lifecycle/numeric-refresh checks. Synthetic meshes were independently
read back for positive volumes, exterior-face coverage and normal orientation.

Those host checks do not validate CUDA. Subsequent user-reported Daint results
in `simple-release-nI3XXr` passed optimized high-resolution channel convergence
and saved-baseline field parity on 1/2/4 ranks, all at 1318 iterations. The
[upwind duct study](../simple_duct/DAINT_RESULTS.md) subsequently passed an
independent recheck of saved GPU results with cache, overlap and distributed
output enabled. This establishes correctness for those cases, not speedup,
kernel occupancy or multi-node scaling. The fine duct's saved loop times
increase from one to four ranks; see the result's timing scope.

The commands below reproduce the optimization checks. Existing public
high-resolution agreement with OpenAccel is their numerical baseline;
OpenAccel needs no rebuild or rerun for them.

## Daint: public correctness first

Use the configured MARS CUDA/Hypre build at
`/capstor/scratch/cscs/gandanie/git/mars-v010-check/build-hypre`, with its working
CUDA-aware MPI environment. Pull and rebuild these targets:

```bash
git pull --ff-only
cmake -S .. -B .
cmake --build . --parallel 4 --target mars_segregated_simple \
  mars_distributed_simple_cuda_gate mars_distributed_halo_cuda_gate \
  mars_simple_output_profile_cuda_gate
```

The following uses only the public channel. It compares a fresh-wrapper control
and optimized 1/2/4-rank native runs against the saved high-resolution baseline.
The shell retains `perf_run` for the next block. With `set -e`, a failed run stops
the sequence. This is an interactive sequence of explicit `srun` calls.

```bash
set -euo pipefail
perf_run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-performance-XXXXXX)
channel_mesh=/capstor/scratch/cscs/gandanie/git/OpenAccel-simple-converged-20260925-151406/channel.exo
baseline=/capstor/scratch/cscs/gandanie/simple-highres-yQnE5j/np1/channel-fields.csv
compare=../tests/reference/openaccel/distributed_simple/compare_fields.py
git rev-parse HEAD > "$perf_run/mars-revision.txt"
git -C _deps/cornerstone_fetch-src rev-parse HEAD > "$perf_run/cornerstone-revision.txt"
sha256sum "$channel_mesh" ./examples/distributed/unstructured/mars_segregated_simple > "$perf_run/sha256.txt"
printf 'Results: %s\n' "$perf_run"

mkdir "$perf_run/control"
srun --account=csstaff --time=00:20:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
  --mesh "$channel_mesh" --output-prefix "$perf_run/control/channel" \
  --advection high-resolution --iterations 2000 --report-every 100 \
  --linear-cache 0 --halo-overlap 0 --field-output distributed \
  2>&1 | tee "$perf_run/control/run.log"
python3 "$compare" "$baseline" "$perf_run/control/channel-fields.json"

for np in 1 2 4; do
  mkdir "$perf_run/np$np"
  srun --account=csstaff --time=00:20:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
    --mesh "$channel_mesh" --output-prefix "$perf_run/np$np/channel" \
    --advection high-resolution --iterations 2000 --report-every 100 \
    --linear-cache 1 --halo-overlap 1 --profile 1 --profile-warmup 10 \
    --field-output distributed 2>&1 | tee "$perf_run/np$np/run.log"
  python3 "$compare" "$baseline" "$perf_run/np$np/channel-fields.json"
done
```

Require convergence, field parity at the existing 1e-6 scaled tolerance, no
linear/continuity failure, one cached graph build per equation and refreshed
numerics/setups on every solve. The 425-node case is a correctness test, not a
scaling benchmark. If CMake uses a Cornerstone source override, record that
checkout's revision instead of `_deps/cornerstone_fetch-src`.

The next short gate uses water (`rho=1000`, dynamic `mu=.001`), high-resolution
advection and a synthetic reversing outlet field. It checks four iterations
against an independently assembled one-rank reference, not water-flow convergence.
It never opens a pump file. SIMPLE uses the explicit true linear residual,
not the legacy solution/RHS magnitude ratio: the first pressure increment is
about `3.07e9 Pa` on this abrupt-start fixture because of its small pseudo-time
step. That magnitude alone neither invalidates the linear solve nor establishes
physical accuracy. See [the residual regression](HYPRE_HOST_TEST.md#simple-convergence-with-physical-units).

```bash
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./mars_distributed_simple_cuda_gate \
  --write-reference "$perf_run/water.bin" --water 1 --high-resolution 1 \
  --backflow -1 --iterations 4 --mesh 8x2x2 2>&1 | tee "$perf_run/water-reference.log"
for np in 1 2 4; do
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_distributed_simple_cuda_gate \
    --reference "$perf_run/water.bin" --water 1 --high-resolution 1 \
    --backflow -1 --iterations 4 --mesh 8x2x2 --builder 1 \
    2>&1 | tee "$perf_run/water-$np.log"
done
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=4 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./mars_distributed_halo_cuda_gate
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=4 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./mars_simple_output_profile_cuda_gate "$perf_run/output-gate"
python3 "$compare" "$perf_run/output-gate/field-fields.csv" "$perf_run/output-gate/field-fields.json"
```

## Daint: measurement after correctness passes

Generate a public 98,304-element mesh with NumPy/netCDF4 in the active Python:

```bash
scale_run=$(mktemp -d /capstor/scratch/cscs/gandanie/simple-scaling-XXXXXX)
python3 ../tests/reference/openaccel/simple_performance/generate_channel.py \
  --nx 64 --ny 16 --nz 16 --output "$scale_run/channel.exo"
```

For each fixed mesh, measure 1/2/4 GPUs on one node, then 8 GPUs over two nodes
and a 4-GPU two-node control. Repeat at least three times in the same allocation
where practical. Measure cache and overlap separately using configurations
`0/0`, `1/0`, `1/1`; all use identical numerical controls and iteration counts.
Below is one 80-iteration sample. Exit 2 with the explicit iteration-limit line
is expected; any other failure is rejected. Retain fields for rank comparisons.

```bash
mkdir "$scale_run/np8"
set +e
srun --account=csstaff --time=00:20:00 --nodes=2 --ntasks-per-node=4 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./examples/distributed/unstructured/mars_segregated_simple \
  --mesh "$scale_run/channel.exo" --output-prefix "$scale_run/np8/channel" \
  --advection high-resolution --iterations 80 --report-every 20 \
  --linear-cache 1 --halo-overlap 1 --profile 1 --profile-warmup 10 \
  --field-output distributed 2>&1 | tee "$scale_run/np8/run.log"
scale_status=${PIPESTATUS[0]}
set -e
test "$scale_status" -eq 0 || {
  test "$scale_status" -eq 2 && grep -q '^NOT CONVERGED: iteration limit iterations=80 ' "$scale_run/np8/run.log"
}
```

Use fresh output prefixes for each sample. Compare equal-iteration fields with
`compare_fields.py`. Report rank-maximum loop seconds/iteration, phase totals,
Krylov iterations and application bytes/wait, alongside revision, node names,
mesh size and timing scope. Repeat representative samples with `--profile 0`.
For weak scaling, use `(nx,ny,nz)` of `(64,16,16)`, `(64,32,16)`, `(64,32,32)`,
`(128,32,32)` on 1/2/4/8 GPUs: approximately constant elements per GPU, fixed
physical geometry. Refinement changes conditioning, so iteration counts must
accompany timings. Larger repeats are warranted only when these measurements
show the per-rank work is too small.

Kernel occupancy/atomics, AMG hierarchy policies and residual collective counts
remain measurement targets. None are declared optimized by the local gates.
