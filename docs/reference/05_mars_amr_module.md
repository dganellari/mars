# MARS AMR Module — Deep Reference

`backend/distributed/unstructured/amr/`. For ElementDomain/cstone integration
see [02_mars_overview.md](02_mars_overview.md).

## 1. Directory Contents

| File | Purpose |
|------|---------|
| `mars_amr.hpp` | Top-level include |
| `mars_amr_manager.hpp` | `AmrManager<ElementTag, KeyType, RealType>` orchestrator |
| `mars_amr_octree.hpp` | Doerfler marking + 1-irregular enforcement |
| `mars_amr_octree_native.hpp` | cstone-driven AMR (unified DD + AMR) — recommended |
| `mars_amr_octree_refine.hpp` | `OctreeAlignedRefine` — emits 8 children per marked parent (no dedup) |
| `mars_amr_hex_refine.hpp` | `HexRefiner` — GPU-native isotropic hex with edge/face dedup |
| `mars_amr_tet_refine.hpp` | `TetRefiner` — Bey's red 1→8 |
| `mars_amr_solution_transfer.hpp` | P0 trilinear transfer by parentage |
| `mars_amr_error_indicator.hpp` | Hex/Tet gradient-based error + Doerfler marking |

## 2. AmrManager — `mars_amr_manager.hpp:190–649`

### Config (struct)
```cpp
struct Config {
    int maxLevels             = 5;
    RealType refineFraction   = 0.3;
    RealType coarsenFraction  = 0.03;
    RealType errorTolerance   = 1e-6;
    int bucketSize            = 64;
    int blockSize             = 256;
    int irregularityIters     = 3;
    MarkingStrategy strategy  = MarkingStrategy::OctreeNative;  // or Doerfler
};
```

### Constructor
```cpp
AmrManager(const Config& config = Config{}) : config_(config), currentLevel_(0) {
    typename AmrOctree<KeyType, RealType>::Config oct;
    oct.refineFraction  = config_.refineFraction;
    oct.coarsenFraction = config_.coarsenFraction;
    oct.maxLevel        = config_.maxLevels;
    amrOctree_.config() = oct;
}
```

### Members
```cpp
Config config_;
int currentLevel_;
int rank_, numRanks_;
InitTimings initTimings_;
DomainPtr domain_;            // ElementDomain<...>
AmrOctree<KeyType, RealType> amrOctree_;
OctreeNativeAmr<KeyType, RealType> octreeNativeAmr_;
```

### initialize() — L248–306
Two overloads:
- `initialize(meshFile, rank, numRanks)` (L248–287) — load mesh, build cstone domain, cache coords, init octrees
- `initialize(domain, rank, numRanks)` (L289–306) — accept pre-built domain

Lap timing (L252–260):
```cpp
auto lap = [&t0]() -> float {
    cudaDeviceSynchronize();
    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = std::chrono::high_resolution_clock::now();
    float ms = std::chrono::duration<float, std::milli>(t1 - t0).count();
    t0 = t1;
    return ms;
};
```

### adaptMesh() — L312–570

```cpp
AmrStats adaptMesh(const RealType* d_errorPerElement,
                   const RealType* d_nodeSolution,
                   DeviceVector& d_newNodeSolution);
```

**Six phases** (each timed with cudaSync + MPI_Barrier):

#### Phase 1: Marking (L356–411, `markTimeMs`)

**OctreeNative path** (L357–367):
```cpp
octreeNativeAmr_.markFromOctree(d_errorPerElement, *domain_);
// no halo exchange, no 1-irregular (octree guarantees conformity)
```

**Doerfler path** (L369–411):
```cpp
amrOctree_.markForRefinement();   // per-element uint8_t marks
exchangeMarks();                  // promote to int, exchangeHalos, convert back
amrOctree_.enforce1Irregular();   // iterative propagation
exchangeMarks();
```

#### Phase 2: Refinement (L449–456, `refineTimeMs`)
```cpp
auto refined = OctreeAlignedRefine::refineLocal(
    coords, conn, d_marks, startIdx, endIdx);
// 8 children per marked, 19 new nodes per child for hex
// numElements_after = 8*numMarked + numUnmarked
// numNodes_after = oldNumNodes + 19*numMarked
```

#### Phase 3: Solution Transfer (L465–486, `transferTimeMs`)
Loop over fields:
```cpp
transferred[f] = SolutionTransfer<...>::transferByParentage(
    refined, oldFields[f], ...);
```
**Pre-sync coord save** (L481–485) — to re-encode SFC keys after cstone redistributes.

#### Phase 4: Rebuild (L488–490, `rebuildTimeMs`)
```cpp
rebuildDomainFromDevice(refined);    // L607–637
amrOctree_.reset(*domain_);
```

#### Phase 5: Reorder/SFC lookup (L492–551, `reorderTimeMs`)
1. Encode pre-sync coords in NEW box → `d_preKeys` (L496–499)
2. `thrust::sort_by_key((preKey, value))` per field (L504–514)
3. For each new local node, binary-search its key in pre-sync sorted table
   (L534–550)

#### Phase 6: Halo Fill (L553–559, `haloFillTimeMs`)
```cpp
for (auto& field : newFields) domain_->exchangeNodeHalo(*field);
```
Fills ghosts not reached by binary search (owner on different rank).

### AmrStats — L22–43

```cpp
struct AmrStats {
    size_t elementsBeforeLocal, elementsAfterLocal;
    size_t elementsBeforeGlobal, elementsAfterGlobal;     // via MPI_Allreduce
    size_t nodesBeforeGlobal, nodesAfterGlobal;
    size_t elementsRefined;
    int level;
    double errorNorm;                                      // globalErrorNorm
    float totalTimeMs, markTimeMs, refineTimeMs,
          transferTimeMs, rebuildTimeMs, reorderTimeMs, haloFillTimeMs;
};
```

### InitTimings — L238–244
```cpp
struct InitTimings {
    float domainSyncTimeMs;   // file read + bbox + cstone sync
    float haloTopoTimeMs;     // HaloData + NodeHaloTopology
    float adjacencyTimeMs;
    float coordCacheTimeMs;
    float octreeTimeMs;
    float totalMs;
};
```

### shouldContinue()
```cpp
bool shouldContinue(const AmrStats& stats) const {
    return currentLevel_ < config_.maxLevels &&
           stats.errorNorm > config_.errorTolerance &&
           stats.elementsRefined > 0;
}
```

## 3. Refinement Strategies

### Doerfler Marking — `mars_amr_octree.hpp:168–223`

```cpp
RealType maxError = thrust::reduce(d_error, +, thrust::maximum());
RealType refineThresh  = refineFraction  * maxError;
RealType coarsenThresh = coarsenFraction * maxError;
// per-element kernel: mark = (error > refineThresh && level < maxLevel) ? 1 :
//                            (error < coarsenThresh && level > 0) ? 2 : 0
```

**1-irregular enforcement** (L189–223): iterate
```cpp
for (int it = 0; it < irregularityIters; ++it) {
    propagateMarksKernel<<<...>>>();   // node-adjacency neighbor check
    if (thrust::equal(marks_old, marks_new)) break;
}
```

### OctreeNative Marking — `mars_amr_octree_native.hpp:156–295`

1. **Scatter errors to leaves** (L210–230): per element, atomically max
   error into leaf indexed by representative SFC key. atomicCAS-emulated
   `atomicMax<double>` (L70–80).
2. **MPI_Allreduce per-leaf errors** (L237–250): copy d_leafError H→ →
   `MPI_Allreduce(MAX)` → copy back. (Could be GPU-direct.)
3. **Convert to weighted counts** (L265–269):
   ```cuda
   leafCount_weighted = count * (1 + error/maxError * bucketSize);
   ```
4. **cstone rebalance decision** (L272–276):
   ```cpp
   cstone::computeNodeOpsGpu(treeLeaves, numLeaves, d_weightedCounts,
                             bucketSize, d_nodeOps);
   ```
   Output `nodeOps[]`: 0/1/8/64/... per leaf
5. **Map leaf decisions to element marks** (L278–291): per local element,
   look up its leaf, check `nodeOps`. ≥8 → mark=1, ==0 → mark=2, else 0

**No mark exchange** — each rank independently applies leaf ops to its
elements (cstone determinism).

## 4. Hex Refinement (1 → 8)

### Node numbering — `mars_amr_hex_refine.hpp:26–41`

```
        7--------6
       /|       /|       Edges:    0-1, 0-2, 0-3, 0-4 (4×)
      / |      / |                 1-2, 1-5         (2×)
     4--------5  |                 2-3, 2-6         (2×)
     |  3-----|--2                 3-7              (1×)
     | /      | /                 4-5, 4-7         (2×)
     |/       |/                  5-6, 6-7         (2×) → 12 total
     0--------1
```

**Per marked element, 19 new nodes** (8 corners reused):
- 12 edge midpoints
- 6 face centers (bot, top, front, back, left, right)
- 1 body center

### HexRefiner::refine() — L411–600

1. **Count children** (L427–464):
   `childCounts[e] = (marks[e] > 0) ? 8 : 1`
   Exclusive scan → `elemPrefix`. Count marked → `markedPrefix`.

2. **Edge keys + dedup** (L497–508):
   - `emitEdgeKeysKernel`: 12 edge keys per marked
   - `EdgeKey::encode(lo, hi) = (hi << 32) | lo`
   - `thrust::sort + unique` → `numUniqueEdges` (typically 6–7× compression)

3. **Face keys + dedup** (L510–532):
   - 6 face keys per marked, sorted 4 nodes packed into Lo/Hi pair
   - `thrust::sort` on `zip(Lo, Hi)` + unique → `numUniqueFaces`
   (10–12× compression typical)

4. **Allocate + compute coords** (L534–566):
   - `numNewNodes = oldNumNodes + numUniqueEdges + numUniqueFaces + numMarked`
   - Edge: midpoint of two corner coords
   - Face: 0.25× sum of 4 corners
   - Body: 0.125× sum of 8 corners

5. **Build child connectivity** (L568–583): `buildChildConnectivityKernel`
   per marked element, 8 children each with connectivity into the 27-point
   subdivision (corners + edge mids + face centers + body)

### Result struct
```cpp
struct Result {
    DeviceVector d_conn0..7, d_x, d_y, d_z;
    DeviceVector d_edgeKeys, d_faceKeysLo, d_faceKeysHi, d_markedPrefix, d_marks;
    size_t numNodes, numElements, numUniqueEdges, numUniqueFaces, numMarked, oldNumNodes;
    KeyType edgeBaseNode() { return oldNumNodes; }
    KeyType faceBaseNode() { return oldNumNodes + numUniqueEdges; }
    KeyType bodyBaseNode() { return oldNumNodes + numUniqueEdges + numUniqueFaces; }
};
```

### OctreeAlignedRefine::refineLocal() — `mars_amr_octree_refine.hpp:240–306`

GPU-native, no edge/face dedup. Emits 19 nodes per child directly:
- Output node IDs: `oldNumNodes + 19*markedPrefix[e] + {0..18}`
- cstone sync handles dedup via SFC ordering (coincident coords → same key → dedup)
- **Pro**: no GPU dedup overhead
- **Con**: 19×numMarked nodes pre-dedup; cstone reduces post-sync

## 5. Tet Refinement (Bey's Red 1 → 8)

`mars_amr_tet_refine.hpp:95–171`

- 4 **corner tets** keep one original vertex + 3 adjacent edge midpoints
- 4 **octahedron tets** triangulate 6 edge midpoints with consistent
  diagonal (m01-m23) for quality preservation

`TetRefiner::refine()` (L198–322) is hex-pattern but only edges (no
face/body). `numNewNodes = oldNumNodes + numUniqueEdges`.

## 6. Solution Transfer

`mars_amr_solution_transfer.hpp`

### Hex `transferByParentage()` — L354–395

Layout post-refinement:
```
[0, oldNumNodes)                                              old nodes
[oldNumNodes, +numUniqueEdges)                                edge mids
[+numUniqueEdges, +numUniqueFaces)                            face centers
[+numUniqueFaces, +numMarked)                                 body centers
```

Two-phase:
1. `copyOldNodeSolutionsKernel` (L247–255): copy `[0, oldNumNodes)`
2. `transferByParentageKernel` (L107–196): for each marked element,
   compute new node values from corners:
   ```cuda
   newSol[findEdge(a,b)]   = 0.5  * (u[a] + u[b])
   newSol[findFace(a,b,c,d)] = 0.25 * (u[a]+u[b]+u[c]+u[d])
   newSol[bodyNode]        = 0.125 * sum_8(u)
   ```
   `findEdge/findFace` = binary search in sorted key arrays.

Multi-field: loop calling per field.

### Tet `transferByParentageTet()` — L397–420
6 edge midpoints only. No face/body.

## 7. Error Indicators

`mars_amr_error_indicator.hpp`

### Hex (least-squares, L19–102)
Per element:
- Build 3×3 Gram `A = Σ(dx_i dx_i^T)`, RHS `b = Σ(dx_i du_i)`
- Cramer's rule for `grad = A^{-1} b`
- Det singularity guard: `det threshold = h^4 * 1e-12`
- `error[e] = h^2 * ||grad||`, h = ³√volume
- isfinite check

### Tet (exact, L105–161)
- `J = [dx1 dx2 dx3]^T` (edges from node 0)
- `grad = J^{-T} [du1 du2 du3]^T` via Cramer
- `error[e] = h^2 * ||grad||`, h = (6V)^{1/3}, V = det(J)/6

### Doerfler marking — L242–261
```cpp
RealType maxError = thrust::reduce(d_error, max);
RealType thresh = refineFraction * maxError;
markElementsKernel: mark[e] = (error[e] > thresh) ? 1 : 0
```

### globalErrorNorm — L263–301
```cpp
RealType localSum = thrust::transform_reduce(
    d_error, d_error+n,
    [] __device__(RealType e) { return e*e; },
    RealType(0), thrust::plus<RealType>());
MPI_Allreduce(&localSum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
return std::sqrt(globalSum);
```

`globalErrorNormOwned()` — sums only `[startIdx, endIdx)` to avoid
double-counting halos.

## 8. Rebuild Flow — `rebuildDomainFromDevice` (L607–637)

```cpp
template<typename ResultType>
void rebuildDomainFromDevice(ResultType& refined) {
    // 1. Wrap refined coords + conn
    auto coords = std::make_tuple(std::move(refined.d_x), ..., std::move(refined.d_z));
    auto conn   = std::make_tuple(std::move(refined.d_conn0), ..., std::move(refined.d_conn7));

    // 2. ElementDomain device-data ctor:
    //    cstone ingests coords + conn → SFC keys → SFC sort → MPI redistribute → element halos
    domain_ = std::make_unique<Domain>(
        std::move(coords), std::move(conn), rank_, numRanks_, config_.bucketSize);

    // 3. Lazy init in the order the solver wants
    (void)domain_->getNodeOwnershipMap();
    cudaDeviceSynchronize();
    (void)domain_->getElementToNodeConnectivity();
    domain_->cacheNodeCoordinates();
}
```

**No host round-trip**: coords + conn move directly device → cstone.

**Node dedup**: multiple ranks may emit same (x,y,z). SFC ordering +
distribution naturally dedup.

## 9. MPI Patterns

| Phase | Pattern | Frequency |
|-------|---------|-----------|
| Doerfler mark exchange | Alltoallv (cstone halo) | 2× per level (before + after 1-irregular) |
| 1-irregular iter | Alltoallv per iter | up to `irregularityIters` |
| OctreeNative leaf errors | Allreduce(MAX) | 1× per level |
| Rebuild (cstone sync) | Allgather + Alltoallv | 1× per level |
| Halo fill | Alltoallv | 1× per field |

**Big win**: OctreeNative eliminates Doerfler's 2–3× Alltoallv on marks.

## 10. Octree Integration with cstone

cstone's octree **is** the AMR data structure. Leaf split/merge directly
controls element refinement.

`markFromOctree` accesses cstone's global tree (L195–202):
```cpp
const auto& globalTree = domain.getDomain().globalTree();
```

Element keys = min of 8 corner-node SFC keys (representative node
convention). After refineLocal, new child coords emitted; cstone re-derives
SFC keys with new bbox during sync. Children's keys fall into parent's
octree range by construction → SFC sort respects hierarchy.

**Unified DD + AMR**: domain decomposition and refinement use the same
tree → automatic load balancing, no asynchrony.

## 11. Gordon-Bell Scaling Gotchas

| Concern | Fix / Recommendation |
|---------|----------------------|
| Doerfler 2× Alltoallv on marks | Use OctreeNative (zero mark exchange) |
| O(N) host loops in 1-irregular | OctreeNative skips this entirely |
| Reorder binary search per node | typically 5–10% of adapt time; linear in `newNumNodes` |
| Solution transfer per-field serialization | batch transfers if multifield |
| Hex dedup `thrust::sort/unique` overhead | use `OctreeAlignedRefine` (cstone deduplicates via SFC) |
| Element halo rebuild via Alltoallv | unavoidable; ~30–40% of adapt time |
| Memory blowup during transfer | limit refineFraction to <10% if memory-constrained |
| Error indicator compute (gradient LSQ) | ~2 ms / 10^6 elem on H100; coarsen if hot path |

## 12. Canonical AMR Loop — `mars_amr_poisson.cu:407–450`

```cpp
using KeyType  = uint64_t;
using RealType = double;

AmrManager<HexTag, KeyType, RealType>::Config cfg;
cfg.maxLevels = 3; cfg.refineFraction = 0.3; cfg.bucketSize = 64;
AmrManager<HexTag, KeyType, RealType> amr(cfg);
amr.initialize(meshFile, rank, numRanks);

cstone::DeviceVector<RealType> d_nodeSolution;
AmrStats stats; stats.elementsRefined = 1;

for (int level = 0; level <= amrLevels; ++level) {
    // 1. SOLVE
    solvePoisson<KeyType, RealType>(amr.domain(), sourceTerm, kernelVariant,
                                    blockSize, maxIter, tolerance, rank,
                                    d_nodeSolution, d_nodeToDof, numOwnedDofs);
    if (level >= amrLevels) break;

    // 2. ERROR INDICATOR
    auto d_error = HexErrorIndicator<KeyType, RealType>::computeError(
        std::get<0>(d_conn).data(), ..., std::get<7>(d_conn).data(),
        d_nodeSolution.data(), d_x.data(), d_y.data(), d_z.data(),
        d_nodeToDof.data(), numElements, blockSize);

    // 3. ADAPT
    cstone::DeviceVector<RealType> d_newSolution;
    stats = amr.adaptMesh(d_error.data(), d_nodeSolution.data(), d_newSolution);
    AmrManager<HexTag, KeyType, RealType>::printStats(stats, rank);
    d_nodeSolution = std::move(d_newSolution);

    // 4. CONVERGENCE
    if (!amr.shouldContinue(stats)) break;
}
```

## 13. Roadmap to 10^9 hex

| Concern | Recommendation |
|---------|----------------|
| Marking | OctreeNative (no mark exchange) |
| Element type | HexTag (19 new nodes per child = smooth interp) |
| Error indicator | coarse approximation (jumps or residual) for hot loops |
| Transfer | parentage-based trilinear |
| Refine fraction | <10% to bound memory peak |
| MPI scaling | OctreeNative + cstone halo rebuild ~100–200 ms / level at 1000 ranks |

## 14. File Index
```
amr/
  mars_amr.hpp                              top-level include
  mars_amr_manager.hpp                      L22–649    AmrManager + Stats + InitTimings
  mars_amr_octree.hpp                       L168–223   Doerfler + 1-irregular
  mars_amr_octree_native.hpp                L70–295    OctreeNative (cstone-driven)
  mars_amr_octree_refine.hpp                L240–306   OctreeAlignedRefine::refineLocal
  mars_amr_hex_refine.hpp                   L26–600    HexRefiner (with edge/face dedup)
  mars_amr_tet_refine.hpp                   L95–322    TetRefiner (Bey's red)
  mars_amr_solution_transfer.hpp            L107–420   transferByParentage(Hex/Tet)
  mars_amr_error_indicator.hpp              L19–301    Hex/Tet gradient + Doerfler + globalErrorNorm
```
