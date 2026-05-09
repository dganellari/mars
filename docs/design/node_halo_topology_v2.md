# NodeHaloTopology v2 — Device-Resident, Gordon-Bell-Suitable Design

**Status:** Plan, not yet implemented. v1 (broken) archived at
`.attic/nodehalotopo_v1_broken.diff`. Working tree restored to commit
`1d90501` (host-driven version that produces correct results but is the
known scaling bottleneck at cube512^3).

**Authoring context:** Cstone branch, Alps/Santis H100s, targeting 10^9
hex elements with CVFEM graph assembly + AMR. The host-driven version
costs ~70 s in NodeHaloTopology construction at cube512^3 regardless of
rank count — that's the wall we have to break.

## 1. Goals

| Goal | Why |
|------|-----|
| **Correctness gate first**: must reproduce host version's `||u||` exactly on cube16/64 before any optimization | v1 failed because correctness and performance changed simultaneously |
| **Device-resident**: no host arrays of SFC keys, ownership claims, peer key tables | Source of the O(70 s) bottleneck — host hash-builds and unordered_map lookups |
| **Compose cstone primitives, don't reinvent** | The deep read showed cstone already exposes everything we need; the failed v1 reimplemented `gatherRanges` and CUDA-aware MPI from scratch |
| **Index-driven pack/unpack** (no NaN sentinel, no `T(0)` fallback) | v1 bug at [domain.hpp:794](../../backend/distributed/unstructured/domain.hpp) silently corrupted ghost zeros |
| **Pre-sized buffers** | v1 race at [domain.hpp:776–780](../../backend/distributed/unstructured/domain.hpp) on `mutable sendBuf_/recvBuf_` resize |
| **Peer-bounded work**: O(boundary) not O(global) | At 1024 ranks, O(global) Allgatherv is the single biggest scaling killer |
| **CUDA-aware MPI by default**, host staging fallback under `CSTONE_HAVE_GPU_AWARE_MPI=OFF` | Cstone already gives this for free if we use its primitives |

## 2. Non-Goals

- Not changing the public `exchangeNodeHalo()` signature. Existing CG /
  CVFEM call sites must continue to work without modification.
- Not changing `HaloData::buildNodeOwnership()` Steps 1–2. Those are
  correct and not the bottleneck.
- Not introducing Hypre / new dependencies. Pure cstone + thrust.

## 3. Architecture (Honest Reuse Tiering)

This isn't "compose cstone primitives." Be specific about what's reused
and what's hand-rolled. Three tiers:

**Tier A — real reuse, types compile without adapters:**
- `mpiSendGpuDirect` / `mpiRecvGpuDirect`
  ([cstone/primitives/mpi_cuda.cuh:31,62](../../cornerstone-octree/include/cstone/primitives/mpi_cuda.cuh)) —
  saves the CUDA-aware MPI + host-staging fallback code path. Direct call.
- `incomingHaloIndices()` / `outgoingHaloIndices()`
  ([cstone/halos/halos.hpp:111-112](../../cornerstone-octree/include/cstone/halos/halos.hpp))
  — bookkeeping accessors we patched in. Direct read.

**Tier B — pattern reuse, structure copied but rewritten:**
- `haloExchangeGpu` ([cstone/halos/exchange_halos_gpu.cuh:34](../../cornerstone-octree/include/cstone/halos/exchange_halos_gpu.cuh))
  — its post-recv / post-send / waitall structure is the template for
  §5's hot path. We don't call it; the buffer layout is element-halo
  shaped, ours is node-halo shaped.

**Tier C — no reuse, hand-rolled:**
- §4 CSR construction (per-peer node-id sets from element-halo ranges)
- §5 pack / unpack kernels — small, ~5 lines each. `gatherRanges` shape
  doesn't match (ours is sparse node IDs, not contiguous ranges); forcing
  it would be degenerate. Verified by mechanical signature check below.

### Signature & call-site evidence (the mechanical check)

Recorded once here so future-us can verify by re-opening the cstone files.

`gatherRanges` ([cstone/halos/gather_halos_gpu.h:22](../../cornerstone-octree/include/cstone/halos/gather_halos_gpu.h)):
```cpp
template<class T, class IndexType>
void gatherRanges(const IndexType* rangeScan, const IndexType* rangeOffsets,
                  int numRanges, const T* src, T* buffer, size_t bufferSize);
```
Inputs are *contiguous range descriptors*. Our send list, after the §4
sort+unique, is a *dense list of arbitrary node IDs*. Wrong shape; do not
reuse. Hand-rolled `dst[i] = src[ids[i]]` instead.

`mpiSendGpuDirect` ([cstone/primitives/mpi_cuda.cuh:31](../../cornerstone-octree/include/cstone/primitives/mpi_cuda.cuh)):
```cpp
template<class T>
auto mpiSendGpuDirect(T* data, size_t count, int rank, int tag,
                     std::vector<MPI_Request>& requests, ...);
```
Our v2 call site:
```cpp
cstone::mpiSendGpuDirect(d_sendBuf_.data() + h_sendOffsets_[i],
                         sendCount, peer, tag, reqs);
```
Types compile. Direct reuse.

### Composition diagram

```
            HOST                           DEVICE
   +----------------+
   | peers (vec<int>)|---------------- already in cstone Domain
   | incomingHalos[p]| (RecvList, single range per peer)
   | outgoingHalos[p]| (SendList, multiple ranges per peer)
   +----------------+
            |
            v
   per-peer construction (one-time, in NodeHaloTopology ctor)
   ----------------------------------------------------------
   For each peer p, derive on device:
     d_sendNodeIds[p]   - dense local node IDs to send to p (CSR)
     d_recvNodeIds[p]   - dense local node IDs to overwrite from p (CSR)

   Pre-allocate:
     d_sendBuf  size = sendOffsets_.back() * sizeof(RealType)
     d_recvBuf  size = recvOffsets_.back() * sizeof(RealType)

            |
            v
   per-call exchangeNodeHalo<T>(arr, nodeToDof) (hot path, no alloc)
   ---------------------------------------------------------------
   1. Pack:    one kernel reads arr[nodeToDof[d_sendNodeIds[i]]] -> d_sendBuf[i]
               (or use cstone gatherRanges if signature fits)
   2. MPI:     mpiSendGpuDirect / mpiRecvGpuDirect on device pointers
   3. Unpack:  one kernel writes d_recvBuf[i] -> arr[nodeToDof[d_recvNodeIds[i]]]
```

Key insight: with **dense node-id CSRs computed once**, the hot path is a
plain gather/scatter + MPI on contiguous device buffers. There is no
sentinel; the unpack writes only at indices that the recv list explicitly
contains. Indices that aren't in the recv list aren't touched.

## 4. The Construction Algorithm (One-Time)

This runs once when the NodeHaloTopology is built (or rebuilt after AMR).
**It is the only non-trivial part of v2.** The hot path is trivial.

### Inputs (all device-resident)

- `peers : vector<int>` (host ok, tiny)
- `incomingHaloIndices_[p]` : `IndexPair<LocalIndex>` — single contiguous
  element range received from peer p ([halos.hpp:111](../../cornerstone-octree/include/cstone/halos/halos.hpp))
- `outgoingHaloIndices_[p]` : `SendManifest` — possibly multiple element
  ranges sent to peer p ([halos.hpp:112](../../cornerstone-octree/include/cstone/halos/halos.hpp))
- `d_conn_local_ids_` : tuple of `NPC` device vectors with **dense local
  node IDs** per element (already built by AdjacencyData)
- `d_nodeOwnership_` : per-node ownership flag from
  `HaloData::buildNodeOwnership()` (already built)
- `startIdx, endIdx` : owned element range

### Algorithm sketch (per peer p)

#### 4.1 Build candidate `recvNodeIds` from incoming element halo

```
incoming element range for p = [S, E) = incomingHaloIndices_[p]
candidate ghost nodes touched by p =
    flatten d_conn_local_ids_[*][S:E]      (NPC * (E-S) entries)
```

Filter to nodes **this rank does not own** (`d_nodeOwnership_[n] == 0`).
Sort + unique. That's `d_recvNodeIds[p]`.

Why this works: cstone already guarantees that for every ghost node N
on this rank, the element halo from N's owner contains at least one
element touching N. So scanning the per-peer incoming element range and
filtering on ownership=0 enumerates exactly the ghost nodes coming from p.

**Kernel 1**: per element in [S, E), emit (peer=p, local_node_id) tuples
for the NPC corner nodes whose ownership==0. Stream-compact + sort
unique per-peer.

#### 4.2 Build `sendNodeIds` from outgoing element halo

```
outgoing element ranges for p = outgoingHaloIndices_[p].ranges()
candidate owned nodes that p needs =
    flatten d_conn_local_ids_[*] over those ranges
```

Filter to nodes **this rank owns** (`d_nodeOwnership_[n] == 1`). Sort +
unique. That's `d_sendNodeIds[p]`.

**Kernel 2**: per element in each outgoing range, emit (peer=p,
local_node_id) for ownership==1. Stream-compact + sort unique per-peer.

#### 4.3 Build CSRs

After per-peer dedup, build flat CSRs:

```
d_sendOffsets : DeviceVector<int>, size numPeers+1
d_sendNodeIds : DeviceVector<LocalIndex>, size sendOffsets.back()
d_recvOffsets : DeviceVector<int>, size numPeers+1
d_recvNodeIds : DeviceVector<LocalIndex>, size recvOffsets.back()
```

Each peer's slice is `d_sendNodeIds[d_sendOffsets[p] : d_sendOffsets[p+1]]`.

#### 4.4 Pre-size buffers

```
d_sendBuf.resize(sendOffsets.back());
d_recvBuf.resize(recvOffsets.back());
```

Done once. **Never resized on hot path.**

### Correctness invariant we want to verify

For each owned node N that's also touched by an incoming halo element
from peer p, that node must NOT also appear in `recvNodeIds[p]` (we
already own it; the value is authoritative locally). The ownership filter
in 4.1 enforces this by construction.

### What about the corner-shared-node bug?

The host version has a known issue ([domain.cu:1002–1008](../../backend/distributed/unstructured/domain.cu)) where ~9 corner-shared
nodes on a 4-rank cube remain doubly-owned because cstone's element halo
doesn't include opposite-rank elements that touch corner-only nodes.

In v2, that issue **does not arise from this construction** — but it's
not fixed either. The `d_nodeOwnership_` we consume is what
`buildNodeOwnership()` Steps 1+2 produced, with that same bug. If a node
is doubly-owned, both ranks include it in their `sendNodeIds` to each
other but neither has it in their `recvNodeIds`, so values remain
self-consistent. The Phase 5 tiebreaker exchange in the host version
exists to fix this; we'll port it as a separate, optional pre-step in
4.0 if validation shows it matters.

For the cube16/64 validation gate, this won't matter (mesh is small
enough we can compare `||u||` exactly against the host version, which
has the same bug). For Gordon-Bell scale we'll need to evaluate.

## 5. The Hot Path: exchangeNodeHalo()

```cpp
template<class VectorType>
void exchangeNodeHalo(VectorType& arr, const int* nodeToDof = nullptr) const
{
    // 1. Pack: one kernel
    packNodeHaloKernel<<<...>>>(
        rawPtr(arr), nodeToDof,
        rawPtr(d_sendNodeIds_), rawPtr(d_sendOffsets_),
        rawPtr(d_sendBuf_),
        numPeers_);

    // 2. MPI: post all recvs first, then all sends, then waitall
    std::vector<MPI_Request> reqs;
    reqs.reserve(2 * numPeers_);
    for (int i = 0; i < numPeers_; ++i) {
        int p = peers_[i];
        size_t recvCount = h_recvOffsets_[i+1] - h_recvOffsets_[i];
        if (recvCount == 0) continue;
        cstone::mpiRecvGpuDirect(
            d_recvBuf_.data() + h_recvOffsets_[i], recvCount,
            p, tagBase_ + epoch_, reqs);
    }
    for (int i = 0; i < numPeers_; ++i) {
        int p = peers_[i];
        size_t sendCount = h_sendOffsets_[i+1] - h_sendOffsets_[i];
        if (sendCount == 0) continue;
        cstone::mpiSendGpuDirect(
            d_sendBuf_.data() + h_sendOffsets_[i], sendCount,
            p, tagBase_ + epoch_, reqs);
    }
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);

    // 3. Unpack: one kernel — index-driven, NO sentinel
    unpackNodeHaloKernel<<<...>>>(
        rawPtr(d_recvBuf_),
        rawPtr(d_recvNodeIds_), rawPtr(d_recvOffsets_),
        rawPtr(arr), nodeToDof,
        numPeers_);

    ++epoch_;  // unique MPI tag per exchange (cstone convention)
}
```

Notes:
- `epoch_` is `mutable`; same convention as cstone's `Halos::haloEpoch_`.
  An MPI collective is required between successive `exchangeNodeHalo()`
  calls if no other collective happens. CG already has Allreduce per
  iter, so this is satisfied.
- `nodeToDof` indirection is preserved exactly as the host version handles
  it.
- We start with **hand-rolled pack/unpack kernels** rather than
  `gatherRanges` because our send is a single index list per peer (not
  multiple disjoint ranges). If profiling later shows `gatherRanges` is
  faster (it's coalesced and tuned), we switch — but the simple kernel
  is correct first.

## 6. Validation Plan (Step-by-Step Gate)

**Rule:** each step must pass before the next is touched. Don't bundle
multiple changes between checksums — that's exactly what hid the v1 bugs.

### Gate 1 — CSR construction matches host
Build `d_sendNodeIds_`, `d_recvNodeIds_`, offsets on device. Then on
cube16, 4 ranks:
1. Build host version (existing code).
2. Build v2 device version, copy back to host.
3. **Per peer**, sort both lists, compare:
   - per-peer node count: must match exactly
   - per-peer sorted node-id sequence: must match exactly

If they don't match, the construction in §4 is wrong; fix that before
touching anything else.

### Gate 2 — pack output matches host
With CSRs validated, run pack on a known input (e.g. `arr[i] = i*3.14`).
Copy `d_sendBuf_` back, compare to host pack. Bit-exact for float, exact
or ULP-1 for double.

### Gate 3 — single MPI exchange round-trip
Pack on rank A, send to B, B unpacks into a zeroed array. Then B re-packs
the ghost slots and sends back to A. `arr` on A at the ghost-mirror slots
should match the original send.

### Gate 4 — one CG iteration matches host
Run one CG iteration with v2 `exchangeNodeHalo()`. Compare `||p||` and
`||Ap||` after the iteration. Must match host within solver tolerance.

### Gate 5 — full Poisson solve matches host
Cube16 + cube64 + cube128. `||u||` and `Max` must match host's working
output to ~1e-12 (double precision Poisson on small mesh is essentially
deterministic).

### Gate 6 — tet mesh validation
Same as Gate 5 but on a tet mesh. This catches NPC=4 vs NPC=8 bugs
before scale-testing exposes them.

### Gate 7 — scale test
Only after Gates 1–6 pass: cube256, cube512 on 4/16/64 ranks. Measure:
- NodeHaloTopology construction time (should be <1 s, not 70 s)
- per-CG-iter halo exchange time (should be ~1 ms latency-bound)
- weak-scaling efficiency

## 7. Where Each Piece Lives in the Code

Files to touch:

| File | What changes |
|------|--------------|
| [backend/distributed/unstructured/domain.hpp](../../backend/distributed/unstructured/domain.hpp) | Replace `NodeHaloTopology` body. Public surface unchanged. Remove `mutable sendBuf_/recvBuf_` resize on hot path; pre-size in ctor. |
| [backend/distributed/unstructured/domain.cu](../../backend/distributed/unstructured/domain.cu) | New construction kernels (`buildSendNodeIdsKernel`, `buildRecvNodeIdsKernel`). New hot-path kernels (`packNodeHaloKernel`, `unpackNodeHaloKernel`). |

Files we **do not** touch:

- `HaloData::buildNodeOwnership()` Steps 1+2 — keep as-is, correct.
- `AdjacencyData` — `d_conn_local_ids_` is already what we need.
- Public `exchangeNodeHalo<T>(...)` signature — internal body changes only.
- Cstone — we use existing accessors and primitives, no patches needed
  beyond the ones we already applied.

## 8. Pitfalls to Avoid (Lessons from v1 + Audit Findings)

1. **Do not change behavior and performance in the same commit.** v1
   merged "make device-resident" with "remove host bookkeeping" with
   "change pack semantics" all at once. When NaN appeared we couldn't
   bisect.
2. **Do not use a value-based sentinel** (NaN, T(0), INT_MAX in float
   bits). Always index-driven: pack and unpack iterate over explicit ID
   lists.
3. **Do not resize device buffers on the hot path.** Mutable shared state
   + concurrent calls = race. Pre-size in ctor.
4. **Do not assume `RecvList` and `SendList` are symmetric.** RecvList
   per peer is one contiguous range; SendList per peer is multiple
   disjoint ranges (`SendManifest` with `offsets_[]` and `scan_[]`).
   Iterate accordingly.
5. **Do not skip Gate 1.** If the CSRs are wrong, every later gate
   fails for a misleading reason.
6. **Do not call `gatherRanges` for our pack step.** Mechanically wrong
   shape — see §3 "signature & call-site evidence." Hand-roll instead.

### Audit finding — pre-existing bugs in `1d90501`

Reading the committed working version closely revealed the T(0) sentinel
([domain.hpp:794](../../backend/distributed/unstructured/domain.hpp))
and `mutable resize` race ([domain.hpp:776–780](../../backend/distributed/unstructured/domain.hpp))
**already exist in `1d90501`**, not just the broken v1. They're inert
under current usage because:
- T(0) is only hit when `nodeToDof[n] < 0` (node is not a DOF),
  which happens cleanly for pure-Dirichlet H1 problems.
- The race needs concurrent `exchangeNodeHalo` calls; CG is
  single-threaded on host, MPI Allreduce serializes between halo
  exchanges.

v2 §5 fixes both regardless: index-driven unpack means no sentinel; ctor
pre-sizing means no hot-path resize. Also: hardcoded `tag=777`
([domain.hpp:813,825](../../backend/distributed/unstructured/domain.hpp))
gets replaced with an epoch counter (cstone convention) so multi-field
exchange in one timestep doesn't tag-collide.

## 9. Open Questions

- **Per-peer kernel launch vs single fused kernel?** Two design choices:
  (a) one launch per peer, simpler but has launch overhead at high peer
  count; (b) single launch over flat CSRs with peer index from
  binary-search of offsets. Start with (a); if profiling shows launch
  overhead matters, switch to (b).
- **CUDA-aware MPI fallback path.** When `CSTONE_HAVE_GPU_AWARE_MPI=OFF`,
  cstone's `mpiSendGpuDirect` automatically host-stages. Verify Alps's
  cray-mpich actually has it on; if not, we eat the staging cost.
- **AMR rebuild cost.** Every adaptMesh re-runs NodeHaloTopology
  construction. The new construction is O(local boundary) so should be
  fast, but verify we're not doing host syncs we missed.
- **Phase 5 tiebreaker port.** Host version has a tiebreaker exchange
  ([domain.cu:1374–1476](../../backend/distributed/unstructured/domain.cu)) for the corner-shared bug. Decide after Gate 5
  whether to port it. If `||u||` matches host on the cube cases, the bug
  doesn't bite our validation; we'd port it only when scale-testing
  reveals it.

## 10. Out-of-Scope Improvements (Future)

These are tempting but explicitly **not** in v2. Land v2 working first.

- Async overlap of halo exchange with local SpMV.
- Persistent MPI requests (`MPI_Send_init` / `MPI_Start`).
- Combining multiple field exchanges into one MPI message.
- Replacing thrust sort+unique with a hash-table-based dedup on device.
