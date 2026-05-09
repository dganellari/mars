# Cornerstone Focus Tree, Global Assignment & Tree Subsystems — Deep Reference

Algorithmic depth on cstone internals load-bearing for MARS. For the API
surface map, see [01_cstone_overview.md](01_cstone_overview.md).

## 1. Hilbert vs Morton Encoding

`include/cstone/sfc/hilbert.hpp`, `include/cstone/sfc/morton.hpp`

Both map 3D coords to 1D SFC keys via bit interleaving. **21 levels for
uint64_t**, **10 levels for uint32_t** (3 bits per level = octree subdivision).

### Morton (magic-number method)
`iMorton<KeyType>(ix, iy, iz)`:
- Expand each coord to 30/63 bits by inserting 2 zero bits between each bit
  via multiply-and-mask (`expandBits`)
- Interleave: `key = xx*4 + yy*2 + zz` (Z-order curve)
- Decode: `compactBits` inverts via shift-XOR parallel prefix

### Hilbert
`iHilbert<KeyType>(px, py, pz)`:
- Per level coarse-to-fine, extract 3-bit octant
- Apply 3D rotation matrix lookup `mortonToHilbert[8] = {0,1,3,2,7,6,4,5}`
- Apply in-place rotations/reflections (L71–82) to maintain continuity

Default is Hilbert because it preserves locality better (fewer jumps
between adjacent cells).

**Invariant**: `keyEnd - keyStart` is always a power of 8. `treeLevel(range)
= log_8(range)`.

Decoding inverts level-by-level (`decodeHilbert` L130–172).

## 2. Cornerstone Tree (csarray)

`include/cstone/tree/csarray.hpp`

**Three invariants** (L19–33):
1. Vector contains keys 0 and max (2^30 or 2^61)
2. Sorted ascending
3. Difference between consecutive keys is a power of 8

Length-N array → N−1 leaf nodes. Leaf `i` spans `[tree[i], tree[i+1])`.

Data: `leaves_: std::vector<KeyType>` of length numNodes+1. **No internal
nodes stored** — only leaves.

### Rebalance — `csarray.hpp:375–389`
1. For each node, compute op: merge (0), keep (1), or split (8/64/512/4096)
2. Exclusive scan `nodeOps` → output indices
3. For each input node, `processNode()` emits 1, 8, or N new leaves
4. Last element always copy of original max

### Search
- `findNodeAbove(tree, n, key)` — lower_bound, first node ≥ key
- `findNodeBelow(tree, n, key)` — upper_bound − 1, last node ≤ key

## 3. Full Octree

`include/cstone/tree/octree.hpp`

Stores both **internal and leaf nodes** in one linear array (no recursion).

### Members of `OctreeData` (L337–355)
- `prefixes` — SFC keys with **Warren-Salmon placeholder bits** encoding
  level (L341, L89). `encodePlaceholderBit(key, 3*level)` packs level
  into high bits.
- `childOffsets` — first-child index (0 for leaves)
- `parents` — parent index, one per group of 8 siblings
- `levelRange` — first node index per level (level-order access)
- `leafToInternal` / `internalToLeaf` — between cornerstone (SFC-sorted)
  and internal (level-order) layouts

### Construction `buildOctreeCpu` — L168–196
1. `createUnsortedLayoutCpu` — emit leaves and internals; for each leaf
   pair, compute common prefix at multiples of 3 bits, emit internal at
   that level (L87–100)
2. `sort_by_key` reorder by prefix (level-order then SFC)
3. `linkTreeCpu` — for each internal, find first child via binary search
   in prefixes one level down; set `childOffsets[i]`, record parent
   (L115–148)

### Navigation
- `child(node, octant)` = `childOffsets[node] + octant`
- `parent(node)` = `parents[(node-1)/8]`
- `isLeaf(node)` = `childOffsets[node] == 0`

## 4. GlobalAssignment

`include/cstone/domain/assignment.hpp`

### Members (L54–300)
- `myRank_, numRanks_, bucketSize_, box_`
- `leaves_` — global tree leaves replicated across all ranks
- `nodeCounts_` — particles per global leaf
- `tree_` — full octree on host; `d_csTree_, d_nodeCounts_` on GPU
- `assignment_: SfcAssignment<KeyType>` — rank → SFC-key boundary
- `exchanges_: SendRanges` — initial particle distribution

### `assign()` — L93–146
1. Compute SFC keys for assigned particles, sort
2. Update global tree via `updateOctreeGlobal` (uses `MPI_Allreduce` MAX
   per node — `tree/update_mpi.hpp:33–52`)
3. Recompute rank boundaries via `makeSfcAssignment`
4. Return required buffer size for next exchange

### Accessors
- `treeLeaves()` L219 — GPU/CPU view
- `nodeCounts()` L226
- `octreeHost()` L243 — always CPU-resident
- `assignment()` L257
- `postExchangeStart(bufDesc)` L262 — first locally-assigned particle index

## 5. SfcAssignment

`include/cstone/domain/domaindecomp.hpp:58–127`

### Members
- `rankBoundaries_` — SFC key per rank boundary (numRanks+1)
- `counts_` — particles per rank
- `numNodesPerRank_` — global leaves per rank
- `treeOffsets_` — cumulative scan of `numNodesPerRank_`

### `findRank(key)` — L96–100
```cpp
auto it = std::upper_bound(begin(rankBoundaries_), end(rankBoundaries_), key);
return int(it - begin(rankBoundaries_)) - 1;
```
O(log numRanks).

### `treeOffsets()` L83
Maps rank `r` → leaf index range `[treeOffsets[r], treeOffsets[r+1])`.

## 6. makeSfcAssignment — Load-Balanced Bucket Split

`domaindecomp.hpp:112–127`

1. `uniformBins(counts, bins, binCounts)` (L34–55):
   - Inclusive scan `counts` → cumulative distribution
   - Bin target = total / numRanks
   - Per rank boundary, binary search in scan for leaf closest to target
   - Residuals: last bin gets fewer particles
2. Gather SFC keys at bin boundaries
3. `numNodesPerRank[i] = treeOffsets[i+1] - treeOffsets[i]`

**Invariant**: each rank gets approx equal counts **respecting tree node
boundaries** (no splitting).

## 7. FocusedOctree Lifecycle

`include/cstone/focus/octree_focus_mpi.hpp`

### Constructor (L55–79)
Single-node tree, octree on host or GPU.

### `converge()` — L599–617
```cpp
while (converged != numRanks_) {
    updateMinMac(assignment, invThetaEff, false);
    converged = updateTree(peers, assignment, globalTreeLeaves, box, scratch);
    updateCounts(particleKeys, globalTreeLeaves, globalCounts, scratch);
    updateGeoCenters();
    MPI_Allreduce(...converged...);
}
```
Per iteration:
- Recompute MAC acceptance radii from geometric centers
- Rebalance based on counts and MACs
- Update leaf counts
- Update geometric centers

### `updateTree()` — L93–193
- Transfer leaves from other ranks via `focusTransfer`
- Enforce mandatory keys (rank/peer boundaries)
- `CombinedUpdate<KeyType>::updateFocus` (`octree_focus.hpp:65–118`):
  rebalance loop with MAC refinement
- Translate global assignment to focus-tree node indices via
  `translateAssignment` (`domaindecomp.hpp:197–220`)
- Sync treelets (peer-owned portions) via `syncTreelets`

### `discoverHalos()` — L530–586
For each assigned focus leaf in `[firstNode, lastNode)`:
1. Compute bounding box from particle positions+h × searchExt × 2
2. Walk octree from root; at each node: if its geometric box overlaps
   bbox AND node is unassigned leaf → mark `macsAcc_[nodeIdx] = 1`

### `computeLayout()` — L588–595
For each leaf: if assigned or halo → include cumulative range, else zero.
Returns nonzero error if any halo not matched to a peer.

### `updateCounts()` — L210–269
1. Local leaf counts via `computeNodeCounts`
2. Global gather for halo counts via `rangeCount`
3. Upsweep (sum-of-children) for internal nodes
4. Exchange counts with peers (MPI p2p)
5. Second upsweep with peer data

## 8. exchangeRequestKeys — Halo Negotiation

`include/cstone/domain/exchange_keys.hpp:44–99`

1. Each rank extracts SFC key ranges of halo nodes via
   `extractMarkedElements` (`layout.hpp:111–141`):
   - Scan layout array; where `layout[i] != layout[i+1]`, add pair
     `(treeLeaves[i], treeLeaves[i+1])`
2. Send to each peer via `mpiSendAsync` (tag `haloRequestKeys`)
3. Receive from peers in any order (`MPI_ANY_SOURCE`)
4. For each `(lowerKey, upperKey)`:
   - Binary search in own tree → `lowerIdx, upperIdx`
   - Add range `[layout[lowerIdx], layout[upperIdx])` to peer's `SendList`
5. `MPI_Waitall` on sends

**Invariant**: peer A's outgoing halos to B = B's incoming halos from A.

## 9. findPeersMac — Dual-Tree Traversal

`include/cstone/traversal/peers.hpp:47–104`

Inputs: `myRank`, `assignment`, `domainTree`, `box`, `invThetaEff` (1/θ + margin).

Algorithm:
1. Span assigned SFC range with intermediate keys (L85–87)
2. Parallel loop over spanning nodes:
   - `locateNode(key, domainTree)` → containing internal node
   - `dualTraversal(...)` recurses on (focus_node, global_node) pairs:
     - `crossFocusPairs(a, b)`: a intersects focus range, b lies outside,
       MAC fails (geometric boxes too close)
     - `p2p(a, b)`: extract rank of b from assignment, mark as peer
3. Return sorted, deduplicated peer list

**MAC test** (L73): `!minMacMutualInt(aBox, bBox, ellipse, pbc)` where
`ellipse = Vec3<T>{box.ilx(), box.ily(), box.ilz()} * box.maxExtent() *
invThetaEff * roundOff` (L58).

## 10. discoverHalos Kernel Logic

`include/cstone/focus/octree_focus_mpi.hpp:530–586`

Per assigned focus leaf:
1. **Bounding box** (L560–578):
   - `inclusive_scan(leafCountsAcc_)` → cumulative counts
   - `computeBoundingBox(x, y, z, h, layout[i], layout[i+1], 2*searchExt, ...)`
2. **Tree traversal** via `findHalos()` CPU or `findHalosGpu`:
   - Start at root
   - Per node: if geom box overlaps bbox AND node is halo (not in
     assigned range) → set `macsAcc_[nodeIdx] = 1`; else if internal →
     recurse

If `accumulate==true`, OR with existing flags.

## 11. Layout Array

`include/cstone/domain/layout.hpp:143–182`

For focus tree of N leaves, `layout[0..N]` is cumulative particle count:
- `layout[i]` = buffer index where leaf i's particles start
- `layout[i+1] - layout[i]` = particles in leaf i (assigned + halos)
- `layout[N]` = total particles present

Per leaf i: assigned → contiguous start; halo → particles present too;
neither → zero-width.

`computeNodeLayout` (L152–182):
```cpp
for (i = 0; i < numLeaves; ++i) {
    bool present = (assignedRange covers i) || flags[leafToInternal[i]];
    layout[i] = present ? leafCounts[i] : 0;
}
exclusive_scan(layout);
```

## 12. Full sync() Data Flow

`include/cstone/domain/domain.hpp:164–206`

1. `distribute()` (`assignment.hpp`): SFC encode → sort → MPI exchange to
   assigned ranks. Returns `keyView` (SFC-sorted assigned keys), `exchangeStart`
2. **Reorder to SFC post-exchange** (L182): `gatherArrays` reorders x,y,z,h
3. `findPeersMac` (L185)
4. `focusTree.converge()` (L189) — first call only
5. `focusTree.updateMinMac()` (L193)
6. `focusTree.updateTree()` (L194)
7. `focusTree.updateCounts()` (L195)
8. `focusTree.discoverHalos()` (L198)
9. `focusTree.computeLayout()` (L200)
10. `halos_.exchangeRequests()` (L201)
11. `updateLayout()` (L203) — reorder to final positions
12. `setupHalos()` (L204) — exchange particle data

Output: x,y,z,h,keys SFC-sorted; assigned `[startIndex, endIndex)`; halos
beyond. Properties at assigned filled; halo property values undefined
until next `exchangeHalos()`.

## 13. Reusable Primitives for MARS NodeHaloTopology v2

| Primitive | File:Line | Use |
|-----------|-----------|-----|
| `findPeersMac` | peers.hpp:47 | Peer list (already exposed via `domain.peers()` if accessor exists) |
| `incomingHaloIndices()` | halos.hpp:111 | per-peer single contiguous element range |
| `outgoingHaloIndices()` | halos.hpp:112 | per-peer multiple ranges (SendManifest) |
| `gatherRanges` | halos/gather_halos_gpu.h:22 | pack scattered device ranges |
| `mpiSendGpuDirect` | primitives/mpi_cuda.cuh:31 | GPU-aware send (or staging) |
| `mpiRecvGpuDirect` | primitives/mpi_cuda.cuh:62 | GPU-aware recv (or staging) |
| `haloexchange` / `haloExchangeGpu` | exchange_halos.hpp / exchange_halos_gpu.cuh | full pack/MPI/unpack pipeline |
| `exchangeRequestKeys` | exchange_keys.hpp:44 | negotiate halo ranges via MPI |
| `computeNodeLayout` | layout.hpp:151 | cumulative offsets from halo flags |
| `extractMarkedElements` | layout.hpp:111 | extract SFC key pairs for flagged elements |

**No primitive for** "send arbitrary keyed data to peers and merge by key"
— MARS would need to implement that for node-level halos. Composition idea
for v2:

```
1. Have peer list (findPeersMac) and focus tree
2. Per iteration:
   a. discoverHalos → mark macsAcc_
   b. exchangeRequests → incomingHaloIndices, outgoingHaloIndices
   c. For each peer:
      - outgoingHaloIndices[peer].scan() → element-buffer offsets
      - element_id = leaves[offsets[range]]
      - Launch GPU kernel: gather scattered element data → contiguous peer buffer
   d. mpiSendGpuDirect (with CUDA-aware MPI)
   e. mpiRecvGpuDirect for incoming
   f. Apply received halos to device particle/node buffers
```

For node halos specifically: derive node sets from element halo ranges,
then build per-peer (sendNodeIds_, recvNodeIds_) CSRs once at sync, reuse
each iteration with a thin pack/MPI/unpack wrapper.

## 14. Gordon-Bell Limitations

| Issue | Where | Impact |
|-------|-------|--------|
| `MPI_Allgather` of tree | update_mpi.hpp:43, assignment.hpp:117 | full tree replicated; scales to ~1000 ranks; latency dominates beyond |
| `MPI_Allreduce` in converge | octree_focus_mpi.hpp:615 | global barrier per iter |
| Host-side rebalancing | csarray.hpp:375 | CPU exclusive_scan + processNode loop; no GPU path for cornerstone rebalance |
| `computeNodeCounts` per rank | csarray.hpp:182 | O(numLeaves) binary searches in particle keys; OMP-parallel but host |
| `extractMarkedElements` | layout.hpp:111 | linear scan, no SIMD/GPU |
| Global octree replicated | assignment.hpp:62 | O(numRanks × bucketSize) memory per rank |
| SPH bbox loop in discoverHalos | octree_focus_mpi.hpp:530 | kernel per pair, batching limited |

Mitigations:
- GPU acceleration for tree rebalance (some functions in
  `cstone/tree/octree_gpu.h` already)
- Distributed global tree (only own subtree + peer portions)
- Async MPI collectives (defer Allreduce to compute phases)

## 15. Algorithmic Citations

- **Cornerstone**: Keller et al., SIAM J. Sci. Comput. 2016
- **Hilbert**: Miki & Umemura, New Astron. 2016 (cited in cstone)
- **Warren-Salmon placeholder bits**: Warren & Salmon SC93 — internal
  tree linking (L115–148)
- **Dual-tree traversal**: standard for FMM/SPH; cstone uses it for peer
  discovery to avoid O(N²) pair evaluations

## 16. Data Structure Diagram

```
Global Tree (all ranks have copy)
├─ leaves[]              cornerstone array (SFC-sorted, invariants 1-3)
├─ nodeCounts[]          particles per leaf
└─ octree                Warren-Salmon prefixes, internal+leaf

Focus Tree (per-rank, distributed)
├─ leaves[]              cornerstone subset (assigned + mandatory)
├─ octree                full internal+leaf
├─ leafCountsAcc[]       local + peer particle counts
├─ centersAcc[] / geoCentersAcc[]    com / geometric centers
├─ macsAcc[]             MAC pass/fail per node (also halo flags)
└─ assignment_           peer rank → leaf index ranges

Halo Bookkeeping (Halos<>)
├─ incomingHaloIndices[rank] = {start, end}    single contiguous range
├─ outgoingHaloIndices[rank] = SendManifest    multiple ranges
└─ layout[0..numLeaves]                         cumulative particle offsets

Per-iteration exchange:
1. discoverHalos → set macsAcc_
2. exchangeRequests → MPI p2p negotiation
3. haloexchange → send/recv data (CPU or GPU-direct)
```
