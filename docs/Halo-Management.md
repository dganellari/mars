# Halo Management

Halo management in MARS handles ghost-element / ghost-DOF data exchange for distributed unstructured FEM. Three layers cooperate:

| Layer | Owns | Purpose |
|---|---|---|
| `HaloData` | node ownership flags | classifies each local node as owned (`1`) or ghost (`0`) |
| `cstone::Domain` | element halos | element-keyed exchange for connectivity, gradients, per-element scalars |
| `NodeHaloTopology` | per-peer node send/recv lists | per-CG-iteration node-DOF exchange via direct CUDA-aware MPI |

Each layer is built lazily: `HaloData` and `NodeHaloTopology` are constructed on the first call to `getNodeOwnershipMap()` (which triggers `ensureHalo()`). The cstone element halo is built during `ElementDomain::sync()`.

## Why three layers?

`cstone::Domain` is element-keyed: it partitions and load-balances elements across ranks. Its halo machinery transports data attached to elements (e.g., per-element error indicators, AMR refinement marks).

Distributed FEM computes on **nodes** (DOFs). A single node is shared by multiple elements; its ownership rank is determined by which rank owns the surrounding elements. To exchange per-node data using only the element halo, you'd send `NodesPerElement` values per element (8 floats per hex8) — most of which are redundant since boundary nodes are shared by 4–8 boundary elements.

`NodeHaloTopology` builds a separate per-node halo on top of cstone's element decomposition. The per-CG-iteration cost drops by 3–8× compared to packing per-element corner arrays, and ghost-DOF semantics become explicit instead of being inferred from element-corner ownership.

For motivation in detail, see [Node Halo Topology](Node-Halo-Topology.md).

## Node ownership: `HaloData`

`HaloData` lives at `backend/distributed/unstructured/domain.hpp` and is built by `buildNodeOwnership()` in `domain.cu`:

```cpp
template<typename ElementTag, typename RealType, typename KeyType, typename AcceleratorTag>
struct HaloData {
    DeviceVector<uint8_t> d_nodeOwnership_;        // size=nodeCount, 0=ghost, 1=owned
    DeviceVector<KeyType> d_haloElementIndices_;   // halo element indices for fast iteration

    HaloData(const ElementDomain&);
    void buildNodeOwnership(const ElementDomain&);
    void buildHaloElementIndices(const ElementDomain&);
};
```

Ownership is determined in two passes:

1. **Mark all locally-owned-element nodes as owned**: every node appearing in this rank's owned element range `[startIndex, endIndex)` gets `1`.
2. **Yield boundary nodes to lower-SFC-key ranks**: any node also appearing in halo elements at indices `[0, startIndex)` (which are owned by lower-rank peers) gets reset to `0`.

Combined with the lowest-rank-wins tiebreaker performed by `NodeHaloTopology` during its construction, every globally-shared node ends up owned by exactly one rank.

### Access

```cpp
const auto& d_nodeOwnership = domain.getNodeOwnershipMap();  // triggers ensureHalo()
const auto& d_haloElems     = domain.getHaloElementIndices();
```

`getNodeOwnershipMap()` is `const` and returns a `DeviceVector<uint8_t>&` of size `domain.getNodeCount()`. Lazy build cost is bounded by the cstone halo width and an `Allgatherv` of size `globalOwnedNodes` (see [Node Halo Topology](Node-Halo-Topology.md)).

## Element-halo exchange (cstone)

cstone provides element-keyed exchange via `Domain::exchangeHalos`:

```cpp
domain.exchangeHalos(std::tie(d_perElementArray), sendBuf, recvBuf);
```

Used in MARS for:

- **Connectivity sync** during `ElementDomain::sync()` (post-redistribute, propagates element-to-node connectivity to halo elements).
- **AMR refinement marks** during `AmrManager::adaptMeshMultiField` (each rank's owned elements broadcast their `uint8_t`/`int` mark to halo peers).

Each array passed in must have size `bufDesc_.size` (= `domain.getElementCount()`, owned + halos).

This element-halo path is **not** used for per-CG-iteration ghost-DOF exchange in the solve loop — that's handled by `NodeHaloTopology` for the reasons above.

## Per-node halo: `exchangeNodeHalo`

The hot path for distributed solvers. Direct CUDA-aware MPI on packed device buffers:

```cpp
template<class VectorType>
void ElementDomain::exchangeNodeHalo(VectorType& nodeArray,
                                      const int* nodeToDof = nullptr) const;
```

- **`nodeArray`**: `DeviceVector<RealType>` of size `nodeCount` (or any DOF-indexed size when `nodeToDof` is provided).
- **`nodeToDof`**: optional remap from local node ID → DOF index. Pass `nullptr` to index `nodeArray` directly by local node ID.

Owned-node slots must hold authoritative local values on entry; ghost-node slots are overwritten with their owner's values on return.

### Typical use in CG

```cpp
// Distributed CG with halo exchange before each SpMV
solver.setOwnedSize(numOwnedDofs);  // enables MPI_Allreduce on dot products
solver.setHaloExchangeCallback(
    [&domain, dofMap = d_node_to_dof.data()](cstone::DeviceVector<RealType>& p) {
        domain.exchangeNodeHalo(p, dofMap);
    });
solver.solve(A, b, x);
```

Each call:

1. Pack via `thrust::for_each` over `sendNodeIds_`.
2. CUDA-aware `MPI_Isend`/`MPI_Irecv` on device pointers (no host round-trip).
3. `MPI_Waitall`.
4. Scatter via `thrust::for_each` over `recvNodeIds_`.

Ghost slots in `nodeArray` are zeroed at the start of every call; nodes not reached by the receive list keep `0` instead of stale memory. This was a real correctness fix — without zeroing, uninitialized memory in unreachable ghost slots produced non-deterministic NaN/inf downstream.

For full topology details and the multi-step build via `MPI_Allgatherv` + `MPI_Alltoallv`, see [Node Halo Topology](Node-Halo-Topology.md).

## Halo for AMR

After AMR refinement and `rebuildDomainFromDevice`, both `HaloData` and `NodeHaloTopology` are rebuilt lazily on the first `getNodeOwnershipMap()` call against the new mesh. The AMR manager triggers this implicitly by calling `amrOctree_.initialize(*domain_)` post-rebuild.

For solution-field carry-over across AMR levels, see [Solution Transfer](Solution-Transfer.md).

## Performance characteristics

| Operation | Frequency | Cost |
|---|---|---|
| `HaloData` build | once per mesh | thrust passes over local elements |
| `NodeHaloTopology` build | once per mesh | `MPI_Allgatherv` (~globalOwnedNodes) + `MPI_Alltoallv` (~boundary nodes) |
| `exchangeNodeHalo` per CG iter | many | one round-trip Isend/Irecv on `~boundary node count` floats per peer |
| cstone element halo | rare (sync, AMR) | per-element data, batched |

For Gordon Bell-scale runs, the dominant cost is `exchangeNodeHalo`. The Allgatherv/Alltoallv at topology build is amortized across all CG iterations of the level.

## See also

- [Node Halo Topology](Node-Halo-Topology.md) — design rationale + topology construction
- [Multi-Rank Support](Multi-Rank-Support.md) — broader MPI integration
- [AMR Module](AMR-Module.md) — refinement pipeline that drives halo rebuilds
- [Solution Transfer](Solution-Transfer.md) — field interpolation across AMR levels
