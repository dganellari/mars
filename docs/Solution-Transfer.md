# Solution Transfer

Solution transfer carries node-keyed fields (e.g. the CG iterate, pressure, velocity) from a coarse AMR level onto its refined successor. Without it, every level starts CG from zero and pays the full convergence cost. With it, level L+1's CG warm-starts from the trilinearly-interpolated level-L solution and converges in a fraction of the iterations.

For the broader AMR pipeline, see [AMR Module](AMR-Module.md).

## Where it lives

Two pieces:

| File | Component |
|---|---|
| `backend/distributed/unstructured/amr/mars_amr_octree_refine.hpp` | `interpolateChildSolutionKernel`, `OctreeAlignedRefine::transferSolution` |
| `backend/distributed/unstructured/amr/mars_amr_manager.hpp` | wiring inside `AmrManager::adaptMeshMultiField` |

The legacy `SolutionTransfer` class in `mars_amr_solution_transfer.hpp` predates the octree-aligned refiner; it's kept for the older hex/tet path but is not used by the distributed pipeline.

## What gets interpolated

`OctreeAlignedRefine::refineLocal` emits **19 new nodes per refined hex parent**:

- 12 edge midpoints (one per parent edge)
- 6 face centers (one per parent face)
- 1 body center (parent's centroid)

The geometry kernel (`emitChildHexesKernel`) writes these positions into `Result.d_x/y/z` at indices `[oldNumNodes + 19*markedPrefix[e], ...)`. The matching solution kernel (`interpolateChildSolutionKernel`) writes solution values at the **same indices** using trilinear interpolation:

```cpp
// edge midpoints: average of 2 corner values
newSol[base + 0]  = 0.5  * (u[0] + u[1]);  // e01
newSol[base + 1]  = 0.5  * (u[1] + u[2]);  // e12
// ... 10 more edges ...

// face centers: average of 4 corner values
newSol[base + 12] = 0.25 * (u[0] + u[1] + u[2] + u[3]);  // fBot
// ... 5 more faces ...

// body center: average of 8 corner values
newSol[base + 18] = 0.125 * (u[0] + u[1] + ... + u[7]);
```

These averages ARE the trilinear shape-function values at edge mid / face center / body center positions for a hex8 with corner-stored DOFs.

## API

```cpp
template<typename KeyType, typename RealType>
struct OctreeAlignedRefine {
    static void transferSolution(
        const KeyType* conn0..conn7,
        const RealType* oldSolution,
        const uint8_t* marks,
        size_t numElements, size_t numNodes,
        cstone::DeviceVector<RealType>& newSolution,  // out, resized to numNodes + 19*numMarked
        int blockSize = 256);
};
```

Output layout (matches `refineLocal`'s coord output):

- `newSolution[0 : numNodes)` — copied from `oldSolution`.
- `newSolution[numNodes : numNodes + 19*numMarked)` — interpolated values at the new nodes.

Multi-field is handled by the manager — call `transferSolution` per field, one per loop iteration:

```cpp
std::vector<DeviceVector<RealType>> transferred(oldFields.size());
for (size_t f = 0; f < oldFields.size(); ++f) {
    OctreeAlignedRefine<KeyType, RealType>::transferSolution(
        /* connectivity slices */, oldFields[f], d_marks.data() + startIdx,
        localCount, numNodes, transferred[f]);
}
```

## The post-sync reordering problem

After `transferSolution` runs, `transferred[f]` is sized for the **pre-cstone-sync** node layout (`numNodes + 19*numMarked`). Then `rebuildDomainFromDevice` is called, which:

1. Hands coords + connectivity to a fresh `cstone::Domain`.
2. cstone re-derives the bounding box from refined coords.
3. cstone sorts by SFC, redistributes nodes by rank, dedupes coincident-coord emissions.
4. The new domain has its own local node ID layout with size `domain_->getNodeCount()`.

The transferred solution arrays don't ride along with cstone's particle set — they're separate device vectors. So after rebuild, `transferred[f]` has the right values but in the wrong order for the new domain.

## Solving it: post-sync key re-encode + binary-search lookup

Implemented in `AmrManager::adaptMeshMultiField`. Five steps:

```cpp
// 1. Save pre-sync coords before cstone consumes them
size_t preNumNodes = refined.numNodes;
DeviceVector<RealType> d_preX(preNumNodes), d_preY(preNumNodes), d_preZ(preNumNodes);
cudaMemcpy(d_preX.data(), refined.d_x.data(), preNumNodes * sizeof(RealType), D2D);
// ... same for Y, Z ...

// 2. Run rebuildDomainFromDevice (cstone moves refined.d_x/y/z into itself,
//    re-derives box, redistributes, dedupes)
rebuildDomainFromDevice(refined);

// 3. Re-encode pre-sync keys with the NEW box that cstone just derived.
//    Same physical coords -> same keys cstone computed during sync.
const cstone::Box<RealType>& newBox = domain_->getBoundingBox();
DeviceVector<KeyType> d_preKeys(preNumNodes);
generateSfcKeys<KeyType, RealType>(d_preX.data(), d_preY.data(), d_preZ.data(),
                                    d_preKeys.data(), preNumNodes, newBox);

// 4. Sort (preKey, value) pairs by key for binary search post-sync
DeviceVector<KeyType> sortedKeys = d_preKeys;
DeviceVector<RealType> sortedVals = transferred[f];
thrust::sort_by_key(thrust::device, sortedKeys, sortedKeys + preNumNodes, sortedVals);

// 5. Per new local node: binary-search its SFC key in sortedKeys, pull the
//    matching value from sortedVals. Owner nodes find their own emission;
//    ghost nodes find a peer's emission (boundary nodes are emitted by all
//    ranks touching them).
const auto& d_newKeys = domain_->getLocalToGlobalSfcMap();
size_t newNumNodes = domain_->getNodeCount();
newField->resize(newNumNodes);

thrust::for_each(thrust::device,
    counting_iterator(0), counting_iterator(newNumNodes),
    [sortedKeys, sortedVals, newKeys, outVals, nSorted] __device__ (size_t i) {
        KeyType target = newKeys[i];
        size_t lo = 0, hi = nSorted;
        while (lo < hi) {  // lower_bound
            size_t mid = (lo + hi) >> 1;
            if (sortedKeys[mid] < target) lo = mid + 1;
            else hi = mid;
        }
        outVals[i] = (lo < nSorted && sortedKeys[lo] == target)
                      ? sortedVals[lo]
                      : RealType(0);
    });

// 6. Fill any node not reached by lookup via the node halo
domain_->exchangeNodeHalo(*newField);
```

Multiple ranks may emit the same SFC key for a shared boundary node — both compute the same trilinear value (same parent corners), so dedup preserves correctness.

## Why option B (re-encode after rebuild) and not option A (pin box across rebuilds)?

Two ways to keep SFC keys consistent across `rebuildDomainFromDevice`:

| | Option A: pin box | Option B (current): re-encode |
|---|---|---|
| Cstone autonomy | constrained — caller dictates box | preserved — cstone derives box freely |
| Cross-level coupling | yes — every level uses the parent's box | none — each level fully independent |
| Moving meshes / contracting domains | breaks (box never adapts) | works (box adapts) |
| MPI cost saved | 6 × `MPI_Allreduce` of 1 RealType | – |
| Extra cost | – | one `generateSfcKeys` on small array (sub-ms) |

We chose **Option B** because:

- Cstone retains full ownership of the box. Pinning would create cross-level coupling that's brittle if geometry ever changes.
- Cstone's box-derivation is 6 batched 1-element Allreduces, rare (per AMR level). Saving them isn't load-bearing.
- The decode-then-re-encode is tiny (sub-ms at our sizes).

## Warm-start in CG

The warm-start needs the example/solver to actually use `d_nodeSolution` as the initial guess. In `mars_amr_cvfem_graph.cu`:

```cpp
cstone::DeviceVector<RealType> b(numOwnedDofs), x(numTotalDofs);
cudaMemcpy(b.data(), d_rhs.data(), ...);
cudaMemset(x.data(), 0, ...);

// If d_nodeSolution arrives non-empty (from the previous level's transfer),
// scatter its values into x's DOF layout. Boundary DOFs get overwritten
// by BC enforcement.
if (d_nodeSolution.size() == nodeCount) {
    thrust::for_each(thrust::device, counting_iterator(0), counting_iterator(nodeCount),
        [nodeToDof = d_node_to_dof.data(),
         nodeSol   = d_nodeSolution.data(),
         xOut      = x.data()] __device__ (size_t i) {
            int dof = nodeToDof[i];
            if (dof >= 0) xOut[dof] = nodeSol[i];
        });
}
```

CG converges to the same fixed point regardless of initial guess — the win is iteration count, not final value.

## Verified speedup

**cube16, 4 ranks, `--amr-levels=2`:**

| Level | CG initial guess | Solve time |
|---|---|---|
| L0 | zero (cold start) | 384 ms |
| L1 (cold-start, no transfer) | zero | ~209 ms |
| L1 (warm-start, transferred) | trilinear-interpolated L0 | **70 ms (~3× faster)** |

L1's `||u||` is bit-identical between cold and warm starts (CG converges to the same answer either way; only iterations differ).

## Limitations

- **Hex8 only.** `interpolateChildSolutionKernel` and `transferSolution` are written for the hex8 element corner topology. Tet equivalents would use barycentric averages on tet edges (existing `SolutionTransfer::transferByParentageTet` handles this for the legacy refiner; not yet wired for `OctreeAlignedRefine`).
- **No high-order or non-trilinear interpolation.** The trilinear assumption is exact for piecewise-trilinear hex8 FE spaces but not for higher-order (Q2, Q3, …).
- **Each field requires its own sort.** `sort_by_key` permutes the value array, so multi-field calls re-sort the keys per field. At our sizes the cost is negligible; for many-field cases (e.g. RANS with k-epsilon + 6 stress components) consider sorting a permutation array once.

## See also

- [AMR Module](AMR-Module.md)
- [Halo Management](Halo-Management.md)
- [Node Halo Topology](Node-Halo-Topology.md) — used to fill ghost slots after lookup
