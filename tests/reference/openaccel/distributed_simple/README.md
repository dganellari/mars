# Distributed SIMPLE validation

The distributed runner retains the single-rank steady, laminar, upwind SIMPLE
kernels. It assembles complete owned rows, solves through the owned-row Hypre
adapter, publishes ghost fields and reduces diagnostics over unique owners.
This integration is a public validation path. CUDA execution, distributed
ElementDomain convergence and arbitrary-mesh ingestion remain unverified.
The existing single-rank driver and shared ElementDomain halo code are unchanged.

## Ownership and communication contract

Each rank supplies every element touching its owned nodes, all relevant boundary
faces, and a unique owner for each node, element and boundary face. Local node
order may differ between ranks; solver node IDs are contiguous within each rank.
Boundary ownership need not match element ownership, but the owner must hold the
face. The partition builder uses the owner of the smallest-key face node.

The setup checks peer counts using sparse synchronous sends and a nonblocking
barrier; no rank-sized all-to-all metadata is built. It rejects malformed offsets,
duplicate peers, duplicate receive slots and asymmetric lists collectively.
`simple_partition` then verifies that each ghost has one owner, send slots are
owned, and peer slots carry exactly matching integer node keys. Solver IDs use
integer messages too. Finally, reverse-added own-element incidence must equal
held incidence at owned nodes. This last check assumes unique element ownership
and a valid local mesh; counts alone are not a general identity proof.

The steady iteration publishes four fused messages per peer (14 doubles per ghost):

1. reconstructed pressure gradient (3);
2. momentum increment and influence coefficient (3 + 3);
3. pressure increment (1);
4. corrected velocity and pressure (3 + 1).

Velocity gradients are not published because this profile has zero reconstructed
advection blend. Enable their exchange before adding a nonzero blend. Dual
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

These are CPU MPI results. Neither nvcc compilation nor device execution has been
validated for this integration. Earlier compile-only checks do not cover the repairs.

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

`mars_segregated_simple_mpi` uses cornerstone's element range, node ownership and
NodeHaloTopology. Its setup prototype downloads domain state, restores exact
coordinates and boundary tags from the replicated public input, builds ownership
on the host, and uploads the runtime arrays. The iteration uses device arrays;
this is not yet fully device-native distributed ingestion.

Startup verifies ownership of every public node and face by identity. The builder
validates key correspondence and element incidence without changing shared halo
code or choosing a new ownership scheme. An incomplete star is a failure to fix
in the mesh/halo layer, not a reason to skip the check or tune a halo factor blindly.

After the short CUDA gates pass, use the saved mesh-only public input. No reference
exports or OpenAccel invocation are used inside this solve:

```bash
simple_mpi_run=$(mktemp -d "$PWD/simple-mpi-XXXXXX")
printf 'Results: %s\n' "$simple_mpi_run"
for np in 1 2 4; do
  srun --account=csstaff --time=00:30:00 --nodes=1 --ntasks-per-node="$np" \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_segregated_simple_mpi \
    --mesh "$PWD/simple-channel-4VjvlT/channel.txt" \
    --output-prefix "$simple_mpi_run/channel-$np" --iterations 2000 --report-every 100 \
    --residual-tol 1e-6 --mass-tol 1e-6 --change-tol 1e-6 \
    2>&1 | tee "$simple_mpi_run/run-$np.log"
done
for np in 1 2 4; do
  python3 ../tests/reference/openaccel/distributed_simple/compare_fields.py \
    "$PWD/simple-native-bMEZdp/channel-fields.csv" \
    "$simple_mpi_run/channel-$np-fields.csv" --tol 1e-6
done
```

Require convergence and field agreement with the saved native single-rank baseline;
an identical iteration count is not required. The comparator uses U=0.1 and
rho U^2=0.01, removes no pressure mean, and needs only the Python standard library.

Production side-set ingestion, GPU-built partition metadata and integration into
the single driver remain separate milestones. This channel gate establishes no
pump, turbulence, arbitrary-mesh or scaling claim.
