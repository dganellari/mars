# Owned-row distributed matrix adapter

This implements the matrix and solve plumbing for a 1/2/4-rank SIMPLE solve:
`backend/distributed/unstructured/fem/segregated/mars_segregated_distributed_matrix.hpp`.
It converts a rank's assembled block CSR into the scalar CSR, RHS and column map that
`HypreGMRESSolver`'s device-map `solve` overload takes. It also unpacks the solution and
evaluates the true distributed residual. It does **not** establish distributed SIMPLE,
nonlinear convergence or pump readiness. Node ownership, element-star completion and native
ingestion belong to other work; this code takes their results as explicit inputs.

Base revision: `eb6a605fd3b61b90bb80692d23ba98cc1009e266` (`cstone`). Shared CMake, the
SIMPLE runtime and drivers, domain/halo, native-input and Poiseuille files are unchanged.
The runtime integration is a separate, unapplied diff.

## Contract

Owned rows must already contain their complete assembled contributions, i.e. every element
and boundary face that touches an owned node. The adapter keeps owned rows only, drops ghost
rows, and never repairs a missing element or coupling. A missing coupling is caught by
the reference and residual gates, not by the adapter.

| Identity | Type | Meaning |
|---|---|---|
| local node `v` | `int` | Index into this rank's runtime arrays (ElementDomain SFC-local order or runner order). Owned and ghost nodes may be interleaved. |
| owned list `owned[k]` | device `int[owned_count]` | Owned local nodes in solver order. Hypre row `C*(first+k)+c` is `(owned[k], c)`. |
| solver node | `GlobalId` (`HYPRE_BigInt`) | Contiguous global id. Rank r owns `[first_r, first_r+owned_r)`, where `first_r` is the exclusive prefix sum of owned counts in rank order. The adapter recomputes `first_r` (`MPI_Exscan`) and checks it. |
| `solver_node[v]` | device `GlobalId[nodes]` | Solver node of every local node, ghosts included. `-1` is allowed only for nodes that no owned row references. |
| solver DOF | `GlobalId` | `C*solver_node + c`, node-major, matching `setPointBlock(C)`. |
| local column | `int` | `C*v + c`. `A.colIndices` holds these, and the separate `solver_dof_map()` translates them. The wrapper uses these as two different arrays. |
| source id | never an input | Exodus/public-file id, for I/O and comparisons only. SFC keys are not solver ids either. |

The real wrapper call, with `System = OwnedRowSystem<C, Solver::Matrix, HYPRE_BigInt>` and
`Solver = mars::fem::HypreGMRESSolver<double,int,cstone::GpuTag>`:

```cpp
System s(MPI_COMM_WORLD, graph.view<C>(nullptr,nullptr), d_owned, owned_count,
         d_solver_node, nodes);                                        // once per mesh/partition
s.update(assembled_view, b.data(), b.size(), local_assembly_failed);   // every solve
cudaMemset(x.data(), 0, s.rows()*sizeof(double));
bool ok = solve_owned(solver, s, b, x);
//   == solver.solve(s.matrix(), b, x, C*first, C*(first+owned), 0, C*total, s.solver_dof_map())
s.unpack(x.data(), x.size(), increment, C*nodes);                      // owned entries only
domain.exchangeNodeHaloBlock(increment_vector, C);                     // ghosts from owners; local order = ElementDomain's, RealType double
auto r = s.residual(halo_complete(increment, C*nodes), b.data(), {1e-13, 1e-10});
ensure(ok && r.passed, ...);
```

The residual check is `sqrt(sum r^2) <= abs + rel*sqrt(sum b^2)`, summed over owned scalar
rows of all ranks before the norm is formed, with the same defaults as `SimpleRunner`.
`halo_complete` is the caller's explicit statement that ghost entries were refreshed after
the last change of owned entries.

Every check ends in a collective, so all ranks throw at the same call and none is left
stranded:
- **Build:** capacity, overflow, owned range, duplicates, row contiguity, ghost ids in range
  and outside the owner's own range, columns without ids, missing diagonal, and the empty-rank
  policy.
- **`update`:** a changed graph, RHS capacity, nonfinite owned values or RHS, and the caller's
  own assembly error.
- **`unpack` and `residual`:** capacity, and a residual requested before the first update.

## Wrapper findings (wrapper not modified)

1. **Unmapped columns are dropped silently.** `buildParCsr` discards any column whose map
   entry is `-1`. The adapter therefore rejects a referenced column without a solver id at
   build time.
2. **A NaN matrix can hang other ranks.** When the matrix contains NaN, `buildParCsr` returns
   before the collective `HYPRE_IJMatrixAssemble`, which can strand the other ranks. The
   adapter's `update` rejects nonfinite owned values collectively before the solve.
3. **Row ranges are `int`.** `IndexType` is `int` for the row range, even though the map
   holds `HYPRE_BigInt`. Rows beyond `min(INT_MAX, max HYPRE_BigInt)` are rejected, not
   narrowed.
4. **Zero-row ranks are unproven.** The wrapper has no explicit empty-range handling. Its
   `m>0` guards and Hypre's empty IJ ranges suggest it works, but this has never been run. The
   default `EmptyRanks::reject` rejects such a rank collectively. `--probe-empty-rank`
   (below) tests the real solve. Switch the default only after that probe passes.

## GPU design

All structure, maps and scratch are built once on the device and reused. Beyond the CSR
itself, the adapter keeps a 4-byte gather index per scalar entry, an 8-byte solver DOF per
local scalar column, the owned list, and a 16 KB reduction scratch buffer.

| Operation | Device work | Host transfers | Collectives |
|---|---|---|---|
| build (once) | 3 validation kernels, 1 scan, 3 fill kernels | 1 int (graph size), 1 long long (nnz), 1 int (status) | `Exscan`, 1 `Allreduce` (total), 2x `Allreduce(BOR)` |
| `update` (per solve) | 1 fused gather: values, RHS pack, finiteness | 1 int (status) | 1 `Allreduce(BOR)`, which can carry the caller's assembly error |
| `unpack` | 1 scatter kernel, asynchronous | none | none (errors deferred to `residual`) |
| `residual` | `W` lanes per row (4 for C=1, 8 for C=3), fixed-order two-stage reduction (reproducible) | 16 bytes | 1 `Allreduce` of 3 doubles |

Compared with the current one-rank solve:
- Offsets and columns are written once, not on every solve.
- The residual no longer needs the `2*rows` squares array or two separate reductions.
- The one extra device pass is the value gather, which moves 20 bytes per entry.

Two more changes would remove the gather:
- For C=1 with owned-first local numbering, the owned block rows already are the scalar CSR.
- For C=3, assembly would need to write scalar-row order directly.

Both are runtime changes and are out of scope here.

## Gates

`gate_common.hpp` runs the same code on host buffers (`host_gate.cpp`, CPU MPI) and on device
buffers (`cuda_gate.cu`, CUDA-aware MPI plus Hypre).

**Fixture:** a generated 6x5x4 node grid with 27-point coupling, nonsymmetric and diagonally
dominant.
- **Blocks:** full CxC blocks, cross-component entries included.
- **Entries:** functions of nontrivial source ids (`1000003+17*perm`).
- **Ownership:** uneven, with irregular seams (e.g. 20/27/30/43 nodes on 4 ranks).
- **Local order:** shuffled; the solver order is a different shuffle.
- **Ghosts:** a complete first ghost layer plus unreferenced second-layer ghosts.
- **Poison:** ghost rows and ghost RHS are NaN, so any read of them fails.

What the gates check:
- **Structure:** every CSR entry, RHS entry and map entry, compared by global row and column
  identity against an independent oracle built from the definition. This includes component
  order.
- **Reuse:** a second update keeps the structure bitwise and overwrites every value.
- **Unpack:** owned entries are written in interleaved local order; ghosts are untouched.
- **Residual:**
  - A known solution passes after a real ghost exchange.
  - Stale ghosts (exchange skipped) fail. The adapter's r² matches the independent oracle r².
  - A single stale referenced ghost fails.
  - With a zero RHS, the absolute tolerance accepts x=0 and rejects x≠0.
- **Omitted coupling:** a dropped ghost-column coupling fails both the oracle and the residual.
- **Empty work:** a rank with no owned rows, and a rank with no local nodes at all. Under
  `reject`, all ranks throw together; under `allow`, results are exact.
- **Collective rejection:** 14 faults injected on one rank only, plus a residual requested
  before the first update. The ctest timeout turns a stranded rank into a failure.
- **CUDA only:** the real Hypre GMRES device-map solve for C=1 and C=3 (`setPointBlock(3)`),
  run twice. The second round reuses the structure with new values, RHS and solution. Checks:
  true residual after exchange, max `|x - x*|` including ghosts `<= 1e-8`, and that the
  unexchanged residual is rejected.

### Executed here

Toolchain: g++ 13.3.0, Open MPI 4.1.6, CMake 3.28.3. No CUDA, Hypre or GPU was available.

| Run | Result |
|---|---|
| standalone CMake + ctest (1/2/4 ranks, `-Wall -Wextra -Wpedantic -Wshadow -Werror`) | 3/3 passed |
| host gate, 1 rank | 42 passed, 0 failed, 8 skipped (cases that need ghosts or 2 ranks) |
| host gate, 2 ranks / 4 ranks | 52 passed, 0 failed, 0 skipped each |
| same under ASan + UBSan, 1/2/4 ranks | same counts, no sanitizer reports |
| `inject.cmake` against a stand-in project named `mars` | deferred include adds the targets; the CUDA branch is entered when `MARS_ENABLE_CUDA`/`HYPRE` are set |

Mutation check on 2 ranks: each of 9 deliberate bugs in a scratch copy of the adapter makes
the gate fail.

| Injected bug | How it was caught |
|---|---|
| swapped component order | oracle |
| values accumulated instead of overwritten | second-update gate |
| local ids used as solver ids | oracle and map check |
| last residual row skipped | residual vs oracle |
| unpack ignoring the owned list | unpack gate |
| contiguity check removed | collective-rejection gate |
| missing-id check removed | collective-rejection gate |
| ghost RHS packed | NaN poison triggers the nonfinite rejection |
| rejection on the faulting rank only | deadlock, killed by the timeout |

These are host results. The device kernels (atomics, the subwarp residual and the thrust scan)
run only in the CUDA gate.

### Compile-only CUDA check (not nvcc, not run)

Clang 18.1.3 in CUDA mode was used, with CUDA 12.9 headers, libdevice and `ptxas` from
NVIDIA's pip wheels. Headers: cornerstone at the pinned `4eb195a`, and Hypre `master`
(`5b58261`, configured without CUDA, 32-bit `HYPRE_BigInt`).

| Translation unit | Host | Device (through `ptxas`) |
|---|---|---|
| adapter + shared gates on device buffers (C=1/3, 64- and 32-bit ids) | compiles, no warnings | sm_80 and sm_90, no warnings |
| `cuda_gate.cu` with the real `HypreGMRESSolver`, `domain.hpp`, `mars.hpp` | compiles | sm_90 |
| `mars_segregated_simple.cu`, before and after the proposed integration diff | compiles | sm_90 |

Clang needed one change, applied only to a scratch copy of `domain.hpp`: its out-of-line
`getConnectivity<I>(size_t)` definition lacks the `MARS_HOST_DEVICE` of its declaration, which
nvcc accepts and clang rejects. nvcc was not available, so nvcc-specific diagnostics remain
possible, and nothing was linked or run on a GPU.

### Not executed: CUDA/Hypre gate on Daint

This code has never been compiled with nvcc or run with Hypre. Apply the patch to the MARS
checkout (`git -C .. am <patch>`). Then, from the configured MARS CUDA/Hypre build directory
(`mars/mlir`), run the block below. The configure step adds one cache entry,
`CMAKE_PROJECT_mars_INCLUDE`. Remove it later with
`cmake -S .. -B . -UCMAKE_PROJECT_mars_INCLUDE`. The same `inject.cmake` also adds the
distributed SIMPLE gates (`../distributed_simple/README.md`), so one configure covers both.

```bash
(
set -euo pipefail
gates=$(cd .. && pwd)/tests/reference/openaccel/distributed_matrix
cmake -S .. -B . -DCMAKE_PROJECT_mars_INCLUDE="$gates/inject.cmake"
cmake --build . --target mars_distributed_matrix_cuda_gate mars_distributed_matrix_host_gate -j4
run=$(mktemp -d "$PWD/distributed-matrix-XXXXXX"); printf 'Results: %s\n' "$run"
git -C .. rev-parse HEAD > "$run/mars-revision.txt"
sha256sum ./mars_distributed_matrix_cuda_gate ./mars_distributed_matrix_host_gate > "$run/sha256.txt"
status=0
for n in 1 2 4; do
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=$n \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_distributed_matrix_cuda_gate 2>&1 | tee "$run/cuda-$n.log" || status=1
  srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=$n \
    --export=ALL --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_distributed_matrix_host_gate 2>&1 | tee "$run/host-$n.log" || status=1
done
# Zero-row rank through the real Hypre path (policy allow). A hang ends at the time limit.
for n in 2 4; do
  srun --account=csstaff --time=00:03:00 --nodes=1 --ntasks-per-node=$n \
    --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
    ~/affinity/bind_numa.sh ./mars_distributed_matrix_cuda_gate --probe-empty-rank 2>&1 | tee "$run/empty-$n.log" || status=1
done
# Information only: timed update/exchange/residual on a 55k-node fixture.
srun --account=csstaff --time=00:05:00 --nodes=1 --ntasks-per-node=4 \
  --export=ALL,MPICH_GPU_SUPPORT_ENABLED=1 --kill-on-bad-exit=1 \
  ~/affinity/bind_numa.sh ./mars_distributed_matrix_cuda_gate --bench 50 2>&1 | tee "$run/bench-4.log" || status=1
printf 'Results: %s (status %s)\n' "$run" "$status"
)
```

A pass means every log ends in `PASS:`. If an empty-rank probe fails or hangs, keep
`EmptyRanks::reject`: the wrapper limitation is then confirmed.

## Minimum later changes to SimpleRunner

Most of steps 2–4 below are now implemented, outside the active runtime, as
`DistributedSimpleRunner` (`../distributed_simple/README.md`). On the host it matches the
one-rank `SimpleRunner` entry by entry on 1, 2 and 4 ranks. What remains is building its
ownership input from ElementDomain. The list below is the original minimum for adapting
`SimpleRunner` itself.

1. **Adopt the adapter (one rank).** Apply `proposed-simple-runtime-integration.diff` (kept
   outside this patch). `LinearSystem` then uses the adapter with the one-rank identity maps.
   Hypre input is identical for one rank; only the residual reduction order changes.
2. **Supply ownership inputs.**
   - Pass the owned list and the `solver_node` map (built from ElementDomain ownership, an
     `Exscan` and one ghost exchange).
   - Size `rows`, `b` and `x` by owned rows.
   - Call `domain.exchangeNodeHaloBlock(increment, C)` between `unpack` and `residual`.
3. **Collective error flags.** Pass assembly error flags into `update(..., failed)`, and make
   `SimpleRunner::check` and the all-outlets-closed test collective.
4. **Global reductions.** Restrict diagnostics, trace moments and flux sums to uniquely owned
   nodes, faces and elements, then all-reduce them.
