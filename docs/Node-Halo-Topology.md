# Node Halo Topology

`NodeHaloTopology` is a per-rank table of "which local nodes do I send to which peer, and which local nodes do I receive from which peer." It enables direct CUDA-aware MPI exchange of node-keyed DOF arrays without going through the cstone element halo.

For the semantics of node ownership and the broader halo layering, see [Halo Management](Halo-Management.md).

## Why a separate node-keyed topology?

`cstone::Domain` is **element-keyed**: it partitions and load-balances elements across ranks via SFC, and its halo exchange transports per-element data. Distributed FEM, by contrast, computes on **nodes** (DOFs). Three issues arise when you try to use only the element halo for per-CG-iteration ghost-DOF exchange:

1. **Wire-data redundancy.** A hex8 has 8 corners. Sending per-element corner arrays via cstone means 8 floats per halo element. But most boundary nodes are shared by 4–8 boundary elements, so the receiver gets the same value 4–8 times. For our cube16/4-rank measurements: ~13K floats per CG iter via cstone vs ~4500 via per-node — about 3× less data on the wire.
2. **Wrong granularity.** cstone halo is keyed by element ownership; FEM needs node ownership semantics. Naive packs require a NaN-sentinel ownership-gating + ghost-zeroing dance to avoid corruption.
3. **Halo-coverage gaps.** cstone's halo width doesn't always cover all elements touching corner-shared nodes (4-way partition junctions in particular). The `Allgatherv((sfc_key, owner))` topology build also resolves ownership duplicates that cstone halo couldn't see — `globalDofs` 4922 → 4913 (matches single-rank), `||u||` matches single-rank to 4 digits.

This is the same dual decomposition pattern used by other distributed FEM frameworks (PETSc DMPlex, Trilinos Tpetra, MFEM): an element/cell distribution for assembly + load balance, plus a separate vertex/DOF distribution for solver communication.

## Data layout

```cpp
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
struct NodeHaloTopology {
    std::vector<int>      peers_;          // peer ranks we send/receive with
    std::vector<int>      sendOffsets_;    // CSR, host (size = peers_.size() + 1)
    std::vector<int>      recvOffsets_;    // CSR, host
    DeviceVector<int>     sendNodeIds_;    // flat owned-node local IDs
    DeviceVector<int>     recvNodeIds_;    // flat ghost-node local IDs

    mutable DeviceVector<RealType> sendBuf_;  // staging
    mutable DeviceVector<RealType> recvBuf_;

    NodeHaloTopology(const ElementDomain&);
};
```

Per-peer index `p`:
- Send to `peers_[p]`: `sendNodeIds_[sendOffsets_[p] : sendOffsets_[p+1]]`.
- Receive from `peers_[p]`: `recvNodeIds_[recvOffsets_[p] : recvOffsets_[p+1]]`.

The CSR offsets and peer list live on the host (used during MPI calls); the node-ID arrays live on the device (used during pack/unpack thrust kernels).

## Construction (one-shot, at first `ensureHalo()`)

`NodeHaloTopology` is built lazily inside `ensureHalo()`, immediately after `HaloData`. The build runs three MPI collectives:

### Step 1: gather ownership table — `MPI_Allgatherv`

Each rank publishes its OWNED `(sfc_key, owner_rank)` pairs:

```cpp
// each rank's owned-key list
std::vector<KeyType> myOwnedKeys;        // node SFC keys this rank claims to own
std::vector<int>     myOwnedLocalIds;    // matching local node IDs

MPI_Allgather(&myOwnedCount, 1, MPI_INT, ownedCounts.data(), ...);
MPI_Allgatherv(myOwnedKeys.data(), myOwnedCount, mpiKeyType,
               globalKeys.data(), ownedCounts.data(), ownedDispls.data(),
               mpiKeyType, MPI_COMM_WORLD);
```

After Allgatherv, every rank has the full `(sfc_key, owner_rank)` table for all globally-owned nodes.

### Step 2: lowest-rank-wins tiebreaker

If multiple ranks claim the same SFC key (which happens for corner-shared nodes that cstone halo didn't expose to all touching ranks), the lowest-rank wins. Other ranks yield those nodes' ownership locally:

```cpp
std::unordered_map<KeyType, int> keyOwner;
for each (rank, key) in globalKeys:
    keyOwner[key] = min(keyOwner[key], rank);

// fix this rank's d_nodeOwnership_:
for each owned node n on this rank:
    if keyOwner[h_sfc[n]] != myrank:
        h_owned[n] = 0;  // yield
```

The fixed ownership is written back to `halo_->d_nodeOwnership_` so subsequent reads (DOF mapping, sparsity) see the corrected state. **This must happen before any caller reads ownership** — `ensureHalo()` builds `NodeHaloTopology` before returning to ensure the corrected map is what callers see.

### Step 3: resolve ghost requests — `MPI_Alltoallv`

Each rank groups its ghost SFC keys by inferred owner (looked up in `keyOwner`), sends those keys to the owners. Owners receive the requested keys, look them up in their own SFC map, and remember the requester as a peer:

```cpp
// per-rank: keys I want from each owner
std::vector<std::vector<KeyType>> reqKeysPerPeer(numRanks);

// alltoall to learn each peer's request sizes
MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, ...);

// alltoallv to ship the actual keys
MPI_Alltoallv(reqKeysFlat.data(), sendCounts.data(), sendDispls.data(), mpiKeyType,
              recvReqKeys.data(),  recvCounts.data(), recvDispls.data(), mpiKeyType,
              MPI_COMM_WORLD);

// owner side resolves received keys to local node IDs
for each received key from rank p:
    sendIdsPerPeer[p].push_back(myKeyToLocal[key]);
```

After this, `recvNodeIds_` (my ghosts grouped by owner) and `sendNodeIds_` (my owned nodes that peers want) are fully built.

## Exchange path: `exchangeNodeHalo`

The hot path. Runs once per CG iteration (and post-AMR for solution-transfer ghost fill):

```cpp
template<class VectorType>
void ElementDomain::exchangeNodeHalo(VectorType& nodeArray,
                                      const int* nodeToDof = nullptr) const;
```

Five steps:

1. **Zero ghost slots.** Nodes never reached by the receive list (rare, but possible on weird partitions) keep `0` instead of stale memory.
2. **Pack.** `thrust::for_each` over `sendNodeIds_` gathers `nodeArray[nodeToDof[localId]]` into `sendBuf_`.
3. **Post receives + sends.** Per-peer non-blocking `MPI_Irecv` and `MPI_Isend` on **device pointers** (CUDA-aware MPI; no host round-trip).
4. **`MPI_Waitall`.**
5. **Scatter.** `thrust::for_each` over `recvNodeIds_` writes received values back into `nodeArray[nodeToDof[localId]]`.

```cpp
// (simplified) pack
thrust::for_each(thrust::device, counting_iterator(0), counting_iterator(sendTotal),
    [arr, snds, sbuf, nodeToDof] __device__ (size_t i) {
        int n = snds[i];
        int idx = (nodeToDof != nullptr) ? nodeToDof[n] : n;
        sbuf[i] = (idx >= 0) ? arr[idx] : T(0);
    });

// post Irecv/Isend on device pointers
for (size_t p = 0; p < topo.peers_.size(); ++p) {
    int peer = topo.peers_[p];
    int rcnt = topo.recvOffsets_[p+1] - topo.recvOffsets_[p];
    if (rcnt > 0)
        MPI_Irecv(rbuf + topo.recvOffsets_[p], rcnt, mpiType, peer, 777, ..., &r);
    // ... matching MPI_Isend ...
}
MPI_Waitall(...);

// scatter
thrust::for_each(thrust::device, counting_iterator(0), counting_iterator(recvTotal),
    [arr, rnds, rbuf, nodeToDof] __device__ (size_t i) {
        int n = rnds[i];
        int idx = (nodeToDof != nullptr) ? nodeToDof[n] : n;
        if (idx >= 0) arr[idx] = rbuf[i];
    });
```

## Distributed CG integration

Wire `exchangeNodeHalo` into the CG halo callback. The solver also needs `setOwnedSize(numOwnedDofs)` to do `MPI_Allreduce` on dot products:

```cpp
ConjugateGradientSolver<RealType, int, cstone::GpuTag> solver(maxIter, tolerance);
solver.setOwnedSize(numOwnedDofs);
solver.setHaloExchangeCallback(
    [&domain, dofMap = d_node_to_dof.data()](cstone::DeviceVector<RealType>& p) {
        domain.exchangeNodeHalo(p, dofMap);
    });
solver.solve(A, b, x);
```

`d_node_to_dof` here is the `nodeCount`-sized array produced by `buildDofMappingGpu` (owned → `[0, numOwnedDofs)`, ghost → `[numOwnedDofs, numTotalDofs)`).

## Verified results (cube16, 4 ranks)

| Metric | Value |
|---|---|
| `globalDofs` | 4913 (matches 1-rank exact) |
| `||u||` | 1.588786 (1-rank: 1.588588, 0.012% off) |
| Topology size | ~3 peers/rank, ~1000 send/recv nodes per peer |
| Wire data per CG iter | ~4500 floats (vs ~13000 via per-element pack) |

Rank 3 yielded 9 duplicate-owned nodes via the Allgatherv tiebreaker — the same 9 corner-shared nodes that previously made `globalDofs` come out as 4922.

## Future direction: dual cstone Domain

Instead of hand-rolling Allgatherv + Alltoallv to discover node ownership and halos, maintain two `cstone::Domain` instances — one keyed by elements (current), one keyed by nodes — sharing the same SFC encoding. cstone's halo discovery on the node-domain would replace `NodeHaloTopology`'s build with cstone-native machinery. Both domains stay aligned because both are built from the same physical mesh via the same Hilbert curve.

This is a future cleanup direction. The current dual-decomposition (cstone for elements + custom `NodeHaloTopology` for nodes) is correct and fast enough for shipping; the cstone-native alternative is a maintainability win, not a correctness/perf one. Tracked for post-Gordon-Bell.

## See also

- [Halo Management](Halo-Management.md) — overall halo layering
- [Multi-Rank Support](Multi-Rank-Support.md) — broader MPI integration
- [AMR Module](AMR-Module.md) — what triggers topology rebuilds
- [Solution Transfer](Solution-Transfer.md) — uses `exchangeNodeHalo` to fill ghost slots after refinement
