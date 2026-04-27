# AMR Module

The AMR (Adaptive Mesh Refinement) module in `backend/distributed/unstructured/amr/` provides GPU-native, multi-rank h-refinement for unstructured hex and tet meshes built on top of `ElementDomain`. The pipeline is octree-aligned: each rank emits children for its locally-owned, marked elements, and `cstone` redistributes them by SFC during the rebuild.

## What's in the module

| File | Class / kernels | Purpose |
|---|---|---|
| `mars_amr.hpp` | umbrella include | exposes the public API |
| `mars_amr_octree.hpp` | `AmrOctree` | refinement-level tracking, Doerfler marking, 1-irregular enforcement |
| `mars_amr_octree_native.hpp` | `OctreeNativeAmr` | drives marking from the cstone leaf decomposition |
| `mars_amr_octree_refine.hpp` | `OctreeAlignedRefine` | hex8 1→8 child emission with new edge/face/body nodes; solution interpolation kernel |
| `mars_amr_hex_refine.hpp` | `HexRefiner` | legacy hex refiner (replaced by `OctreeAlignedRefine` for the multi-rank path) |
| `mars_amr_tet_refine.hpp` | `TetRefiner` | tet4 1→8 red refinement (Bey's algorithm) |
| `mars_amr_error_indicator.hpp` | `HexErrorIndicator`, `TetErrorIndicator`, `ErrorIndicator` | per-element gradient-jump error proxy + Doerfler marking helpers |
| `mars_amr_solution_transfer.hpp` | `SolutionTransfer` | parentage-based interpolation (legacy hex/tet path) |
| `mars_amr_manager.hpp` | `AmrManager` | orchestrates the full cycle: mark → refine → rebuild → transfer |

## Refinement cycle (one level)

```
   d_errorPerElement (size = numElements)
            │
            ▼
   ┌─────────────────────────────────────┐
   │  marking (Doerfler or OctreeNative) │  AmrOctree / OctreeNativeAmr
   │  → d_marks (size = numElements)     │
   └─────────────────────────────────────┘
            │
            ▼
   ┌─────────────────────────────────────┐
   │  octree-aligned refinement          │  OctreeAlignedRefine::refineLocal
   │  → 8 children per marked parent     │
   │    + 19 new node coords / parent    │
   │    (no dedup; cstone handles it)    │
   └─────────────────────────────────────┘
            │
            ▼
   ┌─────────────────────────────────────┐
   │  solution transfer (per field)      │  OctreeAlignedRefine::transferSolution
   │  → trilinear interp at the same     │
   │    19 positions emitted by refine   │
   └─────────────────────────────────────┘
            │
            ▼
   ┌─────────────────────────────────────┐
   │  rebuildDomainFromDevice            │  ElementDomain device-data ctor
   │  → cstone re-derives box, sorts by  │
   │    SFC, redistributes, dedupes,     │
   │    rebuilds halos                   │
   └─────────────────────────────────────┘
            │
            ▼
   ┌─────────────────────────────────────┐
   │  post-sync key re-encode + lookup   │
   │  → match transferred values onto    │
   │    new local node IDs by SFC key    │
   │  → exchangeNodeHalo to fill ghosts  │
   └─────────────────────────────────────┘
            │
            ▼
   d_nodeSolution (warm start for next level's CG)
```

## Public API

### `AmrManager` (the orchestrator)

```cpp
template<typename ElementTag, typename KeyType, typename RealType>
class AmrManager {
public:
    struct Config {
        int    maxLevels        = 3;
        RealType refineFraction = 0.3;     // Doerfler fraction
        RealType errorTolerance = 1e-6;    // stop refining below this error
        int    blockSize        = 256;
        int    bucketSize       = 64;
        MarkingStrategy strategy = MarkingStrategy::OctreeNative;
    };

    AmrManager(Config cfg = {});

    void initialize(const std::string& meshFile, int rank, int numRanks);

    // Single-field overload
    AmrStats adaptMesh(const RealType* d_errorPerElement,
                       const RealType* d_nodeSolution,
                       cstone::DeviceVector<RealType>& d_newNodeSolution);

    // Multi-field overload — takes any number of node-keyed fields
    // (e.g. pressure, velocity components)
    AmrStats adaptMeshMultiField(const RealType* d_errorPerElement,
                                  const std::vector<const RealType*>& oldFields,
                                  const std::vector<DevVector*>& newFields);

    bool shouldContinue(const AmrStats& stats) const;
    int  currentLevel() const;
    ElementDomain<ElementTag, RealType, KeyType, cstone::GpuTag>& domain();
    static void printStats(const AmrStats& stats, int rank);
};
```

`AmrStats` reports per-level: elements before/after (local + global), nodes before/after, elements refined, error norm, total time.

`shouldContinue(stats)` returns false when any of: `currentLevel >= maxLevels`, `errorNorm <= errorTolerance`, or `elementsRefined == 0`.

### Marking strategies

```cpp
enum class MarkingStrategy {
    Doerfler,        // mark elements with err > refineFraction * maxError
    OctreeNative,    // mark via cstone leaf-error reduction + computeNodeOpsGpu
};
```

`Doerfler` runs `AmrOctree::markForRefinement` then enforces 1-irregular at mesh level (with an extra halo exchange of `int`-promoted marks for cross-rank consistency).

`OctreeNative` reduces per-leaf errors via `MPI_Allreduce(MAX)` on the global cstone tree, then calls cstone's `computeNodeOpsGpu` for split decisions. Per-leaf split conformity is the octree's invariant — no extra mark exchange needed.

### `OctreeAlignedRefine`

```cpp
template<typename KeyType, typename RealType>
struct OctreeAlignedRefine {
    struct Result {
        DeviceVector<KeyType>  d_conn0..d_conn7;   // hex connectivity
        DeviceVector<RealType> d_x, d_y, d_z;       // node coords
        size_t numNodes, numElements;
    };

    static Result refineLocal(/* parent conn, coords, marks, sizes */);

    static void transferSolution(/* parent conn, oldSolution, marks, sizes */,
                                 DeviceVector<RealType>& newSolution);
};
```

`refineLocal` emits 8 children per marked element with 19 new local nodes per parent (12 edge midpoints + 6 face centers + 1 body center). No dedup — `cstone::Domain::sync` handles it during the rebuild.

`transferSolution` runs trilinear interpolation at the same 19 positions used by the geometry kernel. See [Solution Transfer](Solution-Transfer.md).

### `ErrorIndicator` + `HexErrorIndicator`

```cpp
template<typename KeyType, typename RealType>
class HexErrorIndicator {
public:
    static cstone::DeviceVector<RealType> computeError(/* hex conn, sol, coords, n2d, count */);
};

template<typename KeyType, typename RealType>
class ErrorIndicator {
public:
    // sums over ALL local elements (incl. halos) — double-counts at boundaries
    static RealType globalErrorNorm(DeviceVector<RealType>& d_error);

    // sums only over [startIdx, endIdx) — owned elements only, correct for distributed
    static RealType globalErrorNormOwned(DeviceVector<RealType>& d_error,
                                          size_t startIdx, size_t endIdx);
};
```

Use `globalErrorNormOwned` whenever `d_error` includes halo elements, which is always the case when computed via per-element kernels over the full local element array. `globalErrorNorm` is kept for cases where `d_error` is already filtered.

## Working example: `mars_amr_cvfem_graph.cu`

Distributed Poisson with adaptive refinement and warm-start CG. Located at `examples/distributed/unstructured/mars_amr_cvfem_graph.cu`.

```cpp
// Setup
amr.initialize(meshFile, rank, numRanks);

cstone::DeviceVector<RealType> d_nodeSolution;
cstone::DeviceVector<RealType> d_errorPerElement;

for (int level = 0; level <= amrLevels; ++level) {
    // Solve the current level. d_nodeSolution arrives non-empty from the
    // previous level (transferred + halo-filled), so CG warm-starts from it.
    solvePoissonGraph(amr.domain(), sourceTerm, kernelVariant, ...,
                      d_nodeSolution, d_errorPerElement);

    if (level >= amrLevels) break;

    // Refine + transfer
    cstone::DeviceVector<RealType> d_newSolution;
    auto stats = amr.adaptMesh(d_errorPerElement.data(),
                                d_nodeSolution.data(), d_newSolution);
    AmrManager<HexTag, KeyType, RealType>::printStats(stats, rank);

    d_nodeSolution = std::move(d_newSolution);
    if (!amr.shouldContinue(stats)) break;
}
```

The solver itself (`solvePoissonGraph`) uses the standard distributed CG path:

```cpp
ConjugateGradientSolver<RealType, int, cstone::GpuTag> solver(maxIter, tolerance);
solver.setOwnedSize(numOwnedDofs);
solver.setHaloExchangeCallback(
    [&domain, dofMap = d_node_to_dof.data()](cstone::DeviceVector<RealType>& p) {
        domain.exchangeNodeHalo(p, dofMap);
    });
solver.solve(A, b, x);

// Halo-exchange the per-node solution before the error indicator,
// so per-element gradient kernels see correct ghost-corner values.
domain.exchangeNodeHalo(d_nodeSolution);
```

See [Halo Management](Halo-Management.md) for the halo layering, [Node Halo Topology](Node-Halo-Topology.md) for `exchangeNodeHalo` internals, [Solution Transfer](Solution-Transfer.md) for the warm-start mechanism.

## Verified results

**cube16 mesh (16³ hex), 4 ranks, `--amr-levels=2`:**

| Level | Elements (global) | `||u||` | `||err||` | Solve time |
|---|---|---|---|---|
| L0 | 4096 | 1.588786 | 0.0753 | 384 ms (cold) |
| L1 | 32768 | 4.513941 | 0.0796 | 70 ms (warm-started) |

L1 is ~3× faster than a cold-start L1 (~209 ms) thanks to the trilinear-interpolated initial guess for CG. See [Solution Transfer](Solution-Transfer.md).

**Single-rank baseline:**
- L0: `||u||=1.588588`, matches 4-rank to 0.012%.
- `globalDofs` matches exactly (4913) thanks to the lowest-rank-wins ownership tiebreaker in `NodeHaloTopology`.

## Multi-rank correctness invariants

Maintained by the AMR pipeline:

1. **Each globally-shared node has exactly one owner.** Enforced by Allgatherv tiebreaker in `NodeHaloTopology` build, triggered automatically by `ensureHalo()` after every refinement.
2. **Ghost-DOF values are sourced from each node's TRUE owner**, not from intermediate ranks holding stale ghost data. Enforced by the per-node send/recv lists.
3. **Per-rank emitted children with coincident coords across ranks dedupe via cstone SFC.** Both ranks compute the same trilinear value (same parent corners), so dedup preserves correctness.
4. **Solution transfer preserves the converged value on coincident pre/post-sync nodes.** Enforced by post-sync key re-encode (re-encode pre-sync coords with cstone's new box, look up by SFC key).

## Known limitations

- **No coarsening (de-refinement) yet.** AMR is monotonic: marked elements split, others stay. Adding coarsening is straightforward but unimplemented.
- **`||err||` discrepancy.** 4-rank `||err||` runs ~2.17× higher than single-rank for the same mesh. Traced to ~9 corner-shared nodes whose RHS gets contributions only from one rank's owned elements (cstone halo-width gap). Doesn't block AMR; tracked separately.
- **Tet path partially tested.** `TetRefiner` and `TetErrorIndicator` exist; multi-rank tet AMR end-to-end is on the validation list.

## See also

- [Halo Management](Halo-Management.md)
- [Node Halo Topology](Node-Halo-Topology.md)
- [Solution Transfer](Solution-Transfer.md)
- [ElementDomain Overview](ElementDomain-Overview.md)
- [Multi-Rank Support](Multi-Rank-Support.md)
