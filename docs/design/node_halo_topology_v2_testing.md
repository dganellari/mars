# NodeHaloTopology v2 — Cluster Testing Guide

The v2 device-resident NodeHaloTopology has been implemented but **not
compiled or run anywhere** (developer machine has no CUDA/MPI). Use this
guide on Alps/Santis to drive it through Gates 1–7.

## State on disk

- `backend/distributed/unstructured/domain.hpp`
  - `NodeHaloTopology` struct gained `buildFromCstoneHalos()` and
    `validateAgainstHost()` declarations (around L439–457)
  - `exchangeNodeHalo()` hot-path rewritten:
    - removed T(0) sentinel — now index-driven (L807 area)
    - removed hot-path resize of `sendBuf_/recvBuf_` (L789–791 area)
    - epoch-counter MPI tags (L848 area) instead of hardcoded 777
- `backend/distributed/unstructured/domain.cu`
  - new free function `buildNodeHaloTopologyHostPath()` wraps the
    original O(global) host body unchanged
  - `NodeHaloTopology` ctor now selects path via env vars; pre-sizes
    `sendBuf_/recvBuf_`
  - new `emitHaloNodesKernel` + `buildFromCstoneHalos` + `validateAgainstHost`
    implementations (~300 lines after the ctor)
- `.attic/nodehalotopo_v1_broken.diff` — earlier abandoned device-resident
  attempt, archived for reference

## Build

Apply cstone patches (skip if already applied):
```bash
bash scripts/patch_cornerstone.sh
```

Then standard CUDA build:
```bash
cd build
cmake .. \
  -DMARS_ENABLE_KOKKOS=OFF \
  -DMARS_ENABLE_CUDA=ON \
  -DMARS_ENABLE_TESTS=ON \
  -DMARS_ENABLE_UNSTRUCTURED=ON \
  -DMARS_ENABLE_FEM_EXAMPLES=ON \
  -DCMAKE_CUDA_ARCHITECTURES=90
make -j mars_cvfem_graph
```

## Runtime selection

Three modes via env vars:

| Env | Behavior |
|-----|----------|
| (unset) | Original host O(global) path. Default. Production. |
| `MARS_NODEHALO_VALIDATE=1` | Build BOTH paths, diff per-peer node-id lists, abort on mismatch. Use for Gate 1. |
| `MARS_NODEHALO_V2=1` | v2 device-resident path only. Use after Gate 1 passes. |

If `VALIDATE` is set, `V2` is implied (validation builds v2 into the live
struct after the diff passes).

## Gate 1 — CSR construction matches host

```bash
MARS_NODEHALO_VALIDATE=1 mpirun -np 4 ./bin/examples/mars_cvfem_graph \
    --mesh /path/to/cube16.mesh
```

**Pass criterion**: stdout shows
```
[Rank N] NodeHaloTopo v2 == host (Gate 1 pass)
```
on every rank, and the program continues past topology construction
without aborting.

**Failure mode 1**: per-peer node-count mismatch.
Means v2 is missing or adding node IDs at some peer. Most likely cause:
ownership filter mismatch with host's tiebreaker. Recovery path: add a
device-side tiebreaker before §4.3 of the design.

**Failure mode 2**: peer-list mismatch.
Means the v2 derivation of peers from `incomingHaloIndices`/
`outgoingHaloIndices` disagrees with the host's `unordered_map` traversal.
Investigate whether some rank has only outgoing or only incoming halo
(my filter requires non-empty in either direction; that's correct).

**Try first on cube16, 4 ranks.** Tiny mesh, cheap to iterate.
Then cube64, then cube128.

## Gate 2 — pack output matches host

This isn't a separate run — Gates 1 and 4 together exercise pack/unpack
on real data. After Gate 1 passes, run a single-iteration smoke test:

```bash
MARS_NODEHALO_V2=1 mpirun -np 4 ./bin/examples/mars_cvfem_graph \
    --mesh /path/to/cube16.mesh --max-iter 1 --tol 1.0
```

Compare `||u||` after iteration 1 with the same flag set vs unset.
Should match to ~1e-12 (double).

## Gate 3 — single MPI round-trip

If Gate 2 deviates, isolate by setting array values to a known pattern
(e.g., `arr[i] = i * 3.14`) before exchangeNodeHalo, then dump the
post-exchange ghost slot values from each rank and compare. The
hot-path code is the same in both modes; if Gate 2 fails it's likely
the v2 CSRs themselves, not the MPI plumbing.

## Gate 4 — one CG iter matches host

```bash
MARS_NODEHALO_V2=1 mpirun -np 4 ./bin/examples/mars_cvfem_graph \
    --mesh cube64.mesh --max-iter 1
```
vs unset. `||u||` after 1 iter must match.

## Gate 5 — full Poisson solve on cube16/64/128

```bash
for sz in 16 64 128; do
  for mode in "" "MARS_NODEHALO_V2=1"; do
    eval $mode mpirun -np 4 ./bin/examples/mars_cvfem_graph \
        --mesh cube${sz}.mesh > out_${sz}_${mode:-host}.log
  done
done
diff out_16_host.log out_16_MARS_NODEHALO_V2=1.log
```
`||u||` and `Max` must match.

## Gate 6 — tet mesh

```bash
MARS_NODEHALO_V2=1 mpirun -np 4 ./bin/examples/mars_ex1_poisson \
    --mesh /path/to/tet_mesh
```
Compare against unset. NPC=4 codepath through emitHaloNodesKernel.

## Gate 7 — scale test

Only after Gates 1–6 pass cleanly:

```bash
for n in 4 16 64; do
  MARS_NODEHALO_V2=1 srun -N$((n/4)) -n$n \
      ./bin/examples/mars_cvfem_graph --mesh cube512.mesh
done
```

Look for the line `Rank 0: NodeHaloTopo[v2] N peers, M send nodes, K recv nodes`
in stdout, with the **wallclock for that phase reported by the existing
PhaseTimer instrumentation** — should be sub-second, not 70 s.

## If validation reveals ownership-tiebreaker is needed

Most likely outcome on real meshes: the host-path tiebreaker via
`MPI_Allgatherv` of all owned keys + `unordered_map` resolves
corner-shared duplicates that v2 doesn't. Two paths forward:

**Option A (simpler)**: keep `HaloData::buildNodeOwnership()` Steps 1+2
as-is, then add a **point-to-point peer-bounded tiebreaker** before
v2's CSR build:
- Each rank exchanges (sfc_key, ownership) with peers from
  `findPeersMac` only
- Apply lowest-rank-wins
- Then run v2 builder on the resolved ownership

This is O(boundary), still device-resident with CUDA-aware MPI.

**Option B (cheaper but less general)**: assume cstone's element halo
already covers all corner-shared neighbors. On cube meshes with
sufficient halo width this holds. Run with `bucketSize >= 2*NPC` and
verify Gate 5 numerically. Document the assumption.

Decide after seeing the actual mismatch on cube256 with realistic rank
counts.

## Rolling back

The host path is unchanged and remains the default. To roll back v2
entirely:
```bash
git checkout backend/distributed/unstructured/domain.cu
git checkout backend/distributed/unstructured/domain.hpp
```
The `.attic/nodehalotopo_v1_broken.diff` is a **separate** abandoned
attempt and should NOT be reapplied — it's there for reference only.
