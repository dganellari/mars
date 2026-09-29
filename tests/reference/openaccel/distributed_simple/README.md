# Distributed SIMPLE validation

The distributed runner retains the single-rank steady, laminar SIMPLE kernels.
Upwind remains the default; the opt-in [high-resolution path](high_resolution.md)
adds limited velocity reconstruction. It assembles complete owned rows, solves through the owned-row Hypre
adapter, publishes ghost fields and reduces diagnostics over unique owners.
The separate [velocity-interpolation option](velocity_interpolation.md) selects
standard or shifted field sampling; its GPU/reference checks are still pending.
Controlled-partition CUDA gates passed on 1/2/4 ranks at revision `3f7ce1e4`,
including reversal and split communicators. Subsequent user-reported Daint runs
also passed native Exodus/ElementDomain convergence and field parity on 1/2/4
ranks: 1277 iterations on the fixed public channel, with maximum scaled velocity
and pressure differences below 6e-13. Later results cover the
[configured oblique fixture](configurable_run.md), the high-resolution public
channel, and the [upwind duct refinement study](../simple_duct/DAINT_RESULTS.md).
These results used standard velocity interpolation and validate the listed cases,
not the shifted option, arbitrary meshes or multi-node scaling.

The distributed path is now the normal `mars_segregated_simple` executable.
See [configurable controls and interactive runs](configurable_run.md). The shared
ElementDomain halo code is unchanged. The compatibility target below still works.

## Ownership and communication contract

Each rank supplies every element touching its owned nodes, all relevant boundary
faces, and a unique owner for each node, element and boundary face. Local node
order may differ between ranks; solver node IDs are contiguous within each rank.
Boundary ownership need not match element ownership, but the owner must hold the
face. The partition builder uses the owner of the smallest-key face node.

The setup checks peer counts using sparse synchronous sends and a nonblocking
barrier; no rank-sized all-to-all metadata is built. It rejects malformed offsets,
duplicate peers, duplicate receive slots and asymmetric lists collectively.
The partition builders then verify that each ghost has one owner, send slots are
owned, and peer slots carry exactly matching integer node keys. Solver IDs use
integer messages too. Finally, reverse-added own-element incidence must equal
held incidence at owned nodes. This last check assumes unique element ownership
and a valid local mesh; counts alone are not a general identity proof.

The steady iteration publishes four fused messages per peer (14 doubles per ghost):

1. reconstructed pressure gradient (3);
2. momentum increment and influence coefficient (3 + 3);
3. pressure increment (1);
4. corrected velocity and pressure (3 + 1).

High-resolution mode bundles velocity gradients (9) and limiter coefficients (3)
into the first round, for 26 doubles per ghost and the same four rounds. Only
owners compute complete limiter bounds; partial ghost bounds are never read. Dual
volumes and boundary factors are complete at owned nodes; ghost partial values
are not used. No reverse-add assembly occurs during the iteration.

Flux history and reversal state are recomputed on each holder from published
inputs. The gate compares every copy by identity. Global continuity cancellation
is an additional necessary check; it cannot detect all compensating copy errors.

Hypre policies, row numbering, exchanges and reductions use the runner's supplied
communicator. Expected list/assembly validation errors are reduced before throwing.
Invalid field layouts, CUDA exchange errors and MPI exchange errors abort the
communicator, before invalid buffers enter communication. Callers must abort on
other unexpected rank-local exceptions (including allocation errors). The gate
only recovers the deliberately tested, globally determined all-outlets-closed case.

## Executed host checks

On 2026-09-27, AppleClang 21 and OpenMPI 5.0.5 passed a strict build with
`-Wall -Wextra -Wpedantic -Wshadow -Werror` and all 43 initial CTest cases. A later index-overflow guard passed the focused
12-case exchange suite, including the added overflow test (44 distinct cases in
total). The independent matrix-adapter suite passed all three 1/2/4-rank tests.

The SIMPLE tests compare complete matrix rows, RHS, velocity, absolute pressure,
influence coefficients, flux copies, outlet traces, reversal flags and diagnostics
against the unchanged single-rank host runner. Coverage includes 1/2/4 ranks,
outlet closure/reopening, convergence on a smaller generated channel, missing
stars, stale ghosts, duplicate faces and independent two-rank subcommunicators.

Additional exchange tests exercise ordinary host-list uploads, fused forward
publication, reverse addition, an isolated rank, integer keys above 2^53, swapped
halo slots and malformed/asymmetric lists. One-rank field errors must terminate
the job promptly; a timeout is a failure. The CSV comparator rejects nonfinite
values, empty files, duplicate IDs, invalid tolerances and malformed records.

```bash
cmake -S tests/reference/openaccel/distributed_simple -B build-dsimple \
  -DCMAKE_CXX_FLAGS="-Wall -Wextra -Wpedantic -Wshadow -Werror"
cmake --build build-dsimple --parallel 4
ctest --test-dir build-dsimple --output-on-failure
```

These are CPU MPI results. The new device partition kernels also have a host
oracle checked against the independent partition builder on 1/2/4 ranks, including
rejection of missing stars, swapped identities and missing boundary tags. The
native Exodus reader is tested with packed/split coordinates and eight malformed
fixtures on 1/2 ranks under ASan/UBSan. Host execution does not validate CUDA.

## Short Daint CUDA gates

From the configured MARS CUDA/Hypre build directory (`mars/mlir`), in its existing
environment. This adds test targets through one cache entry, without replacing
other build options. It requires CUDA-aware MPI. No OpenAccel rebuild is needed.
Remove the injection later with `cmake -S .. -B . -UCMAKE_PROJECT_mars_INCLUDE`.

```bash
set -euo pipefail
git pull --ff-only
cmake -S .. -B . \
  -DCMAKE_PROJECT_mars_INCLUDE="$PWD/../tests/reference/openaccel/distributed_matrix/inject.cmake"
cmake --build . --parallel 4 --target mars_distributed_matrix_cuda_gate \
  mars_distributed_halo_cuda_gate mars_distributed_simple_cuda_gate mars_segregated_simple_mpi
simple_gate_run=$(mktemp -d "$PWD/simple-mpi-gates-XXXXXX")
printf 'Results: %s\n' "$simple_gate_run"
git rev-parse HEAD > "$simple_gate_run/mars-revision.txt"

for np in 1 2 4; do
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_distributed_halo_cuda_gate \
    2>&1 | tee "$simple_gate_run/halo-$np.log"
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_distributed_matrix_cuda_gate \
    2>&1 | tee "$simple_gate_run/matrix-$np.log"
done

srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=1 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./mars_distributed_simple_cuda_gate \
  --write-reference "$simple_gate_run/reference.bin" \
  2>&1 | tee "$simple_gate_run/reference.log"
for np in 1 2 4; do
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_distributed_simple_cuda_gate \
    --reference "$simple_gate_run/reference.bin" --builder 1 \
    2>&1 | tee "$simple_gate_run/simple-$np.log"
done
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=4 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./mars_distributed_simple_cuda_gate \
  --reference "$simple_gate_run/reference.bin" --builder 1 --split 1 \
  2>&1 | tee "$simple_gate_run/subcommunicators.log"
```

The reference above is the existing MARS single-rank runner, which already matched
OpenAccel. Two native distributed iterations must match its full state. Before
long runs, also exercise the sheared closure/reopening case: write a separate
one-rank reference with `--backflow -1 --iterations 12`, then use those same options
for the distributed comparisons. Use `--backflow 1 --iterations 8` for the expected
all-closed rejection. These saved references are tied to the gate executable.

## ElementDomain public driver

`mars_segregated_simple_mpi` reads native Exodus through the MARS C++ reader.
No Python mesh preprocessing or OpenAccel rerun is needed. Rank zero reads the
file arrays; the GPU converts indices, validates side sets and constructs input
for the existing device-data ElementDomain constructor. Cornerstone redistributes
elements and completes stars using SFC node ownership. The vote fallback is rejected.

The lazy local node-key map is built **before** deriving local array sizes. The
constructor's input node count can exceed the nodes held after redistribution;
using that stale count caused the earlier multi-rank download to overrun buffers.
`mars_simple_domain_view_host_test` exercises shrinking, growing and empty local
maps and malformed arrays. The native driver no longer downloads those arrays.

SFC-key matching restores exact file coordinates on the GPU. Face matching,
side-set tags, owned entity lists and solver IDs are built on the GPU. Integer
keys and IDs travel directly through CUDA-aware MPI. Device coverage counters
check unique ownership of every public source node and boundary face. Missing
stars or inconsistent halo identities are errors, not reasons to enlarge halos
blindly. The shared ElementDomain halo implementation is unchanged.

During iteration, fields, outlet moments, residual reductions and convergence
calculations stay on the GPU; MPI receives device buffers. Scratch is reused.
The CPU controls CUDA/MPI/Hypre APIs, handles peer/count metadata and reads small
error/convergence reports for logging. At the end, field rows are packed, gathered
and sorted on the GPU, then downloaded once for CSV file output.

Current scope: one 3D Tet4 block with one selected inlet, one outlet and one or
more wall side sets, configurable physical controls and MPI_COMM_WORLD. Initial source arrays are broadcast
to every rank for coordinate/tag matching and released before iteration. This is
not yet scalable distributed file ingestion, and no performance improvement is
claimed without GPU measurements. Empty owned-row ranks remain rejected by the
Hypre adapter. This establishes no pump or unrestricted mesh capability.

From the existing CUDA/Hypre build directory, build the new targets without
replacing its other options:

```bash
git pull --ff-only
cmake -S .. -B . \
  -DCMAKE_PROJECT_mars_INCLUDE="$PWD/../tests/reference/openaccel/distributed_matrix/inject.cmake"
cmake --build . --parallel 4 --target mars_segregated_simple_mpi \
  mars_simple_mesh_cuda_gate mars_simple_exodus_cuda_gate
```

Run the small topology/input gates and native startup before a long solve:

```bash
(
set -euo pipefail
unset MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
native_run=$(mktemp -d "$PWD/simple-device-XXXXXX")
mesh=/capstor/scratch/cscs/gandanie/git/mars/mlir/simple-native-bMEZdp/channel.exo
printf 'Results: %s\n' "$native_run"
git rev-parse HEAD > "$native_run/mars-revision.txt"
for np in 1 2 4; do
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_simple_mesh_cuda_gate \
    2>&1 | tee "$native_run/mesh-$np.log"
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_simple_exodus_cuda_gate "$native_run/fixtures-$np" \
    2>&1 | tee "$native_run/exodus-$np.log"
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_segregated_simple_mpi \
    --mesh "$mesh" --mesh-format exodus --setup-only 1 \
    --output-prefix "$native_run/setup-$np" \
    2>&1 | tee "$native_run/setup-$np.log"
done
)
```

A setup PASS builds the full runner but runs no iteration and writes no field CSV.
After the short gates pass, run the same public case to convergence and compare
the absolute pressure and velocity against the saved baseline:

```bash
(
set -euo pipefail
unset MARS_OWNERSHIP MARS_HALO_FACTOR MARS_NODEHALO_ALLOW_INCONSISTENT
baseline=/capstor/scratch/cscs/gandanie/git/mars/mlir/simple-native-bMEZdp
simple_mpi_run=$(mktemp -d "$PWD/simple-mpi-XXXXXX")
printf 'Results: %s\n' "$simple_mpi_run"
for np in 1 2 4; do
  srun --account=csstaff --time=00:30:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_segregated_simple_mpi \
    --mesh "$baseline/channel.exo" --mesh-format exodus \
    --output-prefix "$simple_mpi_run/channel-$np" --iterations 2000 --report-every 100 \
    --residual-tol 1e-6 --mass-tol 1e-6 --change-tol 1e-6 \
    2>&1 | tee "$simple_mpi_run/run-$np.log"
  python3 ../tests/reference/openaccel/distributed_simple/compare_fields.py \
    "$baseline/channel-fields.csv" "$simple_mpi_run/channel-$np-fields.csv" --tol 1e-6
done
)
```

Require convergence and field agreement; an identical iteration count is not
required. The comparator uses U=0.1 and rho U^2=0.01, removes no pressure mean,
and needs only the Python standard library.
