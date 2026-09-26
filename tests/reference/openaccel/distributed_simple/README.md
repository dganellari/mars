# Distributed SIMPLE on explicit ownership

`backend/distributed/unstructured/fem/segregated/mars_segregated_simple_distributed.hpp` runs
`SimpleRunner`'s steady SIMPLE iteration on one rank's partition. The kernels are unchanged;
only assembly scope, exchanges and reductions differ. Linear solves go through the owned-row
adapter (`../distributed_matrix`), and ghost refreshes through the fused
`mars_segregated_halo_exchange.hpp`.

Ownership is an explicit input. Native ElementDomain ingestion and element-halo completion
provide it, and are not duplicated here. This establishes the distributed algorithm and its
parity on generated meshes. It does not establish GPU execution, the public-mesh driver, or
pump readiness.

## Contract

Inputs per rank (`SimpleOwnership` plus a `SimpleInput`-shaped local mesh):
- **Complete stars:** every element that touches an owned node, plus every boundary face of the
  elements this rank holds. Element node order is the source order on every rank, so sample
  and face orientation agree across ranks.
- **Nodes and their owners:**
  - `owned_nodes` in solver order;
  - `solver_node` for every local node (ghosts carry their owner's id; wider id types are
    range-checked, not narrowed);
  - halo lists covering every non-owned local node, each received from that node's actual owner.
- **Owned elements and faces:** `owned_elements` and `owned_faces` are unique across ranks; a
  boundary face belongs to its element's owner.

Assembly scope and reductions:

| Quantity | Computed from | Ghost entries |
|---|---|---|
| dual volume, boundary factor | all held elements/faces | partial; never read at ghosts |
| grad p | complete star | published every iteration (3) |
| grad u | complete star | not published: only multiplied by `velocity_blend = 0` in this profile. Publish it (9) before enabling a blend. |
| momentum/pressure rows, d | owned rows only (adapter drops ghost rows) | d published (3) |
| stored interior/boundary flux, reversal flags, outlet trace | recomputed on every rank that holds the element or face, from identical published inputs | copies must agree; the gate checks every copy |
| mass divergence | complete star | partial; never read at ghosts |
| diagnostics, outlet area moments | owned nodes, faces and element samples | one `Allreduce` each |

There are no reverse additions anywhere, so nothing can be counted twice. The global identity
`sum(continuity) = Q_in + Q_out` holds only if every copy of a flux history is the same. The
existing cancellation check therefore detects any divergence between copies.

Per iteration:
- **Ghost exchange:** 4 fused rounds, 14 doubles per ghost:
  1. grad p;
  2. [du, d];
  3. phi;
  4. [u, p].
- **Small collectives:** 7, each one int or up to 3 doubles. These are the momentum-assembly
  flag, 2 update gates, 2 true residuals, the outlet moment and the flux-update flag. Evaluating
  diagnostics adds 2 more. Decisions that are already identical on every rank (the solver's
  return, the residual verdict, program order) are checked locally, without a collective.
- **Nothing at setup:** the mutation run below showed that a setup exchange of volume and
  boundary factor was unnecessary, so it was removed.

Every error is collective:
- assembly and flux-update flags;
- nonfinite owned values;
- Hypre failure;
- the true residual;
- all outlet faces closed: the global area moment is zero on every rank at once.

## Gates (host, executed)

**Fixture:** a generated Tet4 channel. 16×4×4 hexes give the public channel's sizes: 425 nodes,
1536 tets, and 32/32/512 inlet/outlet/wall faces. Its one-rank reference also converges at
iteration 1277. It is still a generated mesh, not the public fixture.
- **Partition:** uneven x slabs, plus a y split at 0.7 on 4 ranks, with jagged seams.
- **Node ownership:** the lowest rank among the elements touching the node.
- **Local order:** nodes, elements and faces are shuffled per rank.
- **Poison:** ghost-only arrays that must never be read are NaN-poisoned. The velocity gradient
  gets a huge finite value instead, because it is multiplied by the zero blend.

**Reference:** the unchanged one-rank `SimpleRunner` with its dense host solver. The distributed
side uses a test-only gathered dense LU.

**Compared by global node/element/face and component, entry by entry:**
- the momentum and pressure block rows, plus the scalar rows the adapter hands to the solver;
- RHS;
- u, p, d, divergence;
- every held copy of the interior and boundary flux, trace and reversal flags;
- all diagnostics.

`ctest --test-dir build-dsimple`: 22/22 passed.

| Case | 1 rank | 2 ranks | 4 ranks |
|---|---|---|---|
| 2 iterations, 16×4×4, every entry | pass, worst 8e-13 | pass | pass |
| sheared backflow, 12 iterations: 12 faces close and gradually reopen; on 4 ranks the closed faces span 2 ranks | pass | pass | pass |
| uniform backflow: all outlet faces close; both sides throw the same error at the same iteration, all ranks together | pass | pass | pass |
| convergence 8×2×2: same iteration (2742); final u, p within 3e-15; all 2742 diagnostics within 9e-14 | pass | pass | pass |
| convergence 16×4×4 (run once, outside ctest): iteration 1277 on both sides; final u, p within 2e-15; all diagnostics within 1.3e-13 | reference | – | pass |
| fault: one stale ghost after a publish | refuses (no ghosts) | detected | detected |
| fault: halo element missing from a star | refuses (no halo) | detected | detected |
| fault: an owned face counted twice | detected | detected | detected |

Mutation check (4 ranks, 2-iteration and backflow cases): each deliberate change to the runner
must fail the gate.

| Change | Result |
|---|---|
| grad p not published | detected |
| d not published | detected |
| final [u, p] not published | detected |
| pressure increment applied to owned nodes only | detected |
| outlet moment over held faces instead of owned faces | detected |
| outlet moment not all-reduced | detected: ranks without outlet faces throw alone, the others hang in the pressure solve, and the timeout kills the run (exit 124) |
| node diagnostics over all local nodes | detected |
| setup volume/factor publish removed | not detected; the round was unnecessary and has been removed |

## Not executed: CUDA/Hypre gate

`simple_gate.cu` builds the same gate with `MARS_REPLAY_CUDA`. Both sides then solve with
Hypre, and the default tolerance is 1e-8. Only the compile check below has been done:

- **Toolchain:** clang 18 in CUDA mode, CUDA 12.9 headers and `ptxas`, cornerstone `4eb195a`,
  Hypre `master` headers.
- **Result:** host and sm_90 device compile, with no warnings from these files.

The CUDA reference must be its own one-rank job, because `SimpleRunner`'s Hypre solve uses
`MPI_COMM_WORLD`. The target is added by `../distributed_matrix/inject.cmake`. From the MARS
CUDA/Hypre build directory, after checking out the branch:

```bash
(
set -euo pipefail
gates=$(cd .. && pwd)/tests/reference/openaccel/distributed_matrix
cmake -S .. -B . -DCMAKE_PROJECT_mars_INCLUDE="$gates/inject.cmake"
cmake --build . --target mars_distributed_simple_cuda_gate -j4
run=$(mktemp -d "$PWD/distributed-simple-XXXXXX"); printf 'Results: %s\n' "$run"
git -C .. rev-parse HEAD > "$run/mars-revision.txt"
status=0
launch() { local n=$1; shift
  srun --account=csstaff --time=00:10:00 --nodes=1 --ntasks-per-node=$n \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_distributed_simple_cuda_gate "$@"; }
launch 1 --write-reference "$run/ref-2it.bin" 2>&1 | tee "$run/ref-2it.log" || status=1
launch 1 --write-reference "$run/ref-shear.bin" --backflow -1 --iterations 12 2>&1 | tee "$run/ref-shear.log" || status=1
for n in 1 2 4; do
  launch $n --reference "$run/ref-2it.bin" 2>&1 | tee "$run/2it-$n.log" || status=1
  launch $n --reference "$run/ref-shear.bin" --backflow -1 --iterations 12 2>&1 | tee "$run/shear-$n.log" || status=1
done
printf 'Results: %s (status %s)\n' "$run" "$status"
)
```

A pass is `PASS:` in every log.

## What remains for production

1. **Build `SimpleOwnership` from ElementDomain.** This belongs to the ingestion and halo work:
   - owned node list and solver ids: `Exscan`, plus one ghost exchange of ids;
   - owned elements: the local range `[startIndex, endIndex)`;
   - owned faces: faces of owned elements;
   - halo lists: `getNodeHaloTopology()`.

   The runner's local node order must then be ElementDomain's order.
2. **Resolve the NodeHaloTopology send-list issue first.** An owner that misses a requested key,
   or answers for a key it doesn't own, gives shifted or stale ghosts.
   `FieldExchange` rejects count mismatches, but not a wrong owner.
3. **Switch the driver.** `mars_segregated_simple` constructs `DistributedSimpleRunner` for
   more than one rank. `SimpleRunner` stays as the one-rank reference.
