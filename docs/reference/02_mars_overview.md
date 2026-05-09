# MARS Unstructured Backend Reference

This is the API-surface map of `backend/distributed/unstructured/`. For
deep dives see:
- [03_mars_cvfem_kernels.md](03_mars_cvfem_kernels.md) — CVFEM kernels
- [05_mars_amr_module.md](05_mars_amr_module.md) — AMR
- [06_mars_solvers.md](06_mars_solvers.md) — solvers + sparse matrix

## 1. Architecture

| Subdir | Purpose |
|--------|---------|
| `.` | `domain.hpp` / `domain.cu` — `ElementDomain<Tag,Real,Key,Acc>` wrapping cstone |
| `fem/` | CVFEM kernels (12 variants), tensor ops, coloring, sparsity, H1 FE space |
| `solvers/` | CG, BiCGSTAB, GMRES, Hypre PCG, AMG preconditioner |
| `amr/` | Octree-driven hex/tet refinement, transfer, error indicators |
| `utils/` | Mesh I/O (binary, MFEM, Exodus), reference elements, BC handling |

cstone integration: `domain.hpp:490`
`using DomainType = cstone::Domain<KeyType, RealType, AcceleratorTag>`.
`incomingHaloIndices()` / `outgoingHaloIndices()` accessed at `domain.cu:1244–1245`.

## 2. ElementDomain — `domain.hpp:485–1091`

### Member layout
| Member | Type | Lazy? | Purpose |
|--------|------|-------|---------|
| `d_elemSfcCodes_` | `DeviceVector<KeyType>` | no | element SFC codes |
| `d_conn_keys_` | `Tuple(NPC × DV<KeyType>)` | no | element→node connectivity in **SFC keys** (NPC=4 tet, 8 hex) |
| `d_props_` | `DeviceVector<RealType>` | no | per-node characteristic size `h` |
| `d_boundary_` | `DeviceVector<uint8_t>` | no | per-node boundary flag |
| `box_` | `cstone::Box<RealType>` | no | bounding box |
| `d_elemToNodeMap_` | `DeviceVector<KeyType>` | no | representative node (min SFC) per elem |
| `d_localToGlobalSfcMap_` | `DeviceVector<KeyType>` | yes | local-id → global SFC (sorted) |
| `d_localToGlobalNodeMap_` | `DeviceVector<KeyType>` | yes | local-id → global node idx |
| `adjacency_` | `unique_ptr<AdjacencyData>` | yes | local-id connectivity + node→elem CSR |
| `halo_` | `unique_ptr<HaloData>` | yes | node ownership + halo elem indices |
| `nodeHaloTopo_` | `unique_ptr<NodeHaloTopology>` | yes | per-peer node send/recv lists |
| `coordCache_` | `unique_ptr<CoordinateCache>` | yes | decoded (x,y,z) from SFC |
| `originalCoords_` | `unique_ptr<OriginalCoordinates>` | yes | exact mesh-file coords |

### Element index model
SFC keys are globally unique scalars. cstone partitions SFC space:
```
[0, startIndex())              low-rank halos
[startIndex(), endIndex())     owned elements (this rank assembles)
[endIndex(), getElementCount())  high-rank halos
```
Local node IDs `[0, nodeCount)` are dense and contiguous per rank, built
during `sync()`. `d_localToGlobalSfcMap_` is sorted ascending → binary
searchable on device. Built by sort+unique on element connectivity SFC keys
(`domain.cu:2346–2380`).

### Constructor flow
File constructor (`domain.hpp:513`):
1. Read mesh (binary `.float32/.double`, MFEM `.mesh`, Exodus `.exo`)
2. MPI_Allreduce global bbox (`domain.hpp:1266–1276`)
3. Move to GPU; call device constructor

Device constructor (`domain.hpp:540`):
1. Move device coords + connectivity in
2. Call `sync()`

`sync()`:
1. Compute element SFC keys from centroids
2. cstone::Domain::sync → SFC sort + MPI redistribute + element halo
3. Extract unique node SFC keys
4. Sort, assign dense local IDs
5. Allocate per-rank node arrays, scatter coords/boundary by local ID

### Public accessors (selected)
```cpp
size_t getNodeCount() const;           // not lazy
size_t getElementCount() const;        // not lazy
auto getConnectivity<I>() const;       // d_conn_keys_[I]
auto getElementToNodeConnectivity();   // LAZY: tuple of local-id arrays
auto getLocalToGlobalSfcMap();         // LAZY
auto getNodeOwnershipMap();            // LAZY
auto getNodeX/Y/Z();                   // LAZY
auto getNodeToElementOffsets/List();   // LAZY (CSR)
pair<size_t,size_t> localElementRange();   // [startIdx, endIdx)
vector<pair> haloElementRanges();          // [0,startIdx) ∪ [endIdx,N)
bool isLocalElement(idx);
```

### Template instantiations (8)
`{TetTag, HexTag} × {float, double} × {unsigned, uint64_t} × {cstone::GpuTag}`

## 3. AdjacencyData — CSR connectivity

Triggered via `buildAdjacency()` or lazily via `ensureAdjacency()`
(`domain.cu:2345, 2380`).

### Build flow
1. **buildNodeToElementMap()** — for each element, iterate corner local IDs,
   accumulate per-node element list, dedup via `thrust::sort_by_key + unique`
2. **createElementToNodeLocalIdMap()** — for each element corner, look up
   SFC key → binary search `d_localToGlobalSfcMap_` → store local id

### Storage
| Field | Type | Size |
|-------|------|------|
| `d_nodeToElementOffsets_` | `DV<KeyType>` | nodeCount + 1 |
| `d_nodeToElementList_` | `DV<KeyType>` | sum of element degrees |
| `d_conn_local_ids_` | `Tuple(NPC × DV<KeyType>)` | elementCount × NPC |

Consumed by FEM assembly (gather), coloring (conflict detection), AMR error
indicators.

## 4. HaloData / NodeHaloTopology

### HaloData::buildNodeOwnership() — `domain.cu:939–1011`

**Step 1** (`domain.cu:979–987`): Mark nodes in `[startIdx, endIdx)` as
owned (=1) via `markElementNodesOwnership<KeyType, NPC>` kernel.

**Step 2** (`domain.cu:989–1000`): For elements in `[0, startIdx)` (halos
from lower-ranked partitions), mark their corner nodes as not-owned (=0).
Implements "lowest-SFC-key rank wins" tiebreaker.

**Known limitation** (`domain.cu:1002–1008`): ~9 corner-shared nodes on
4-rank cube remain doubly-owned because cstone's halo width misses
opposite-rank elements that touch corner-only nodes. NodeHaloTopology Phase 5
fixes this with a tiebreaker exchange.

`buildHaloElementIndices()` caches halo elem indices for fast iteration.

### NodeHaloTopology — `domain.cu:1214–1725` (~500 lines)

Replaces O(global) `MPI_Allgatherv` with O(boundary) using cstone's
`incomingHaloIndices()` / `outgoingHaloIndices()`.

**Six phases**:
1. **Per-element source-rank array** (L1254–1278): `d_elemOwnerRank[e]` =
   peer rank or INT_MAX, filled from `incomingHaloIndices` ranges
2. **`atomicMin` lowest peer touching each node** (L1305–1329): walk halo
   element corners, atomicMin into `d_lowestPeerTouching[nodeId]`
3. **Recompute ownership** (L1331–1341): node owned ⇔ locally-touched AND
   `rank ≤ lowestPeerTouching`. Fixes corner-shared bug.
4. **Boundary mask** (L1343–1372): mark nodes appearing in any halo elem
5. **Tiebreaker exchange via CUDA-aware MPI** (L1374–1476): publish
   boundary node SFC keys + ownership claims among peers; merge claims;
   apply lowest-rank tiebreaker. Point-to-point `Isend/Irecv` on device
   pointers, 3-tag pattern (count, keys, claims)
6. **Build send/recv lists** (L1500+): per-peer
   `recvNodeIds_[p]` (ghost nodes owned by p) and `sendNodeIds_[p]`
   (owned nodes sent to p). Stored as CSR.

### exchangeNodeHalo() — `domain.hpp:765`
```cpp
template<class VectorType>
void exchangeNodeHalo(VectorType& nodeArray, const int* nodeToDof = nullptr) const;
```
1. **Pack** (L786–797): gather `arr[nodeToDof?nodeToDof[n]:n]` into `sendBuf_`
2. **MPI** (L800–829): post all `MPI_Irecv` on `recvBuf_` device ptrs, all
   `MPI_Isend` on `sendBuf_` device ptrs, `MPI_Waitall`
3. **Scatter** (L831–844): write `recvBuf_` back into `arr` at ghost nodes

**Invariant**: For every ghost node N, cstone's element halo must include
≥1 element owned by N's owner rank. Holds for typical FEM partitions with
halo width ≥ full-element-layer.

## 5. CVFEM Assembly — `fem/`

See [03_mars_cvfem_kernels.md](03_mars_cvfem_kernels.md) for kernel-level
detail. Files:
- `mars_cvfem_assembler.hpp` — dispatch enum + launcher
- `mars_cvfem_hex_kernel*.hpp` — 12 variants
- `mars_cvfem_coloring.hpp` — graph coloring
- `mars_sparse_matrix.hpp` + `mars_sparsity_builder.hpp`
- `mars_h1_fe_space.hpp`
- `mars_cvfem_node_data.hpp` — AoS struct (72 B aligned)

Halo exchange in solve: `exchangeNodeHalo(Ap)` before each SpMV in CG loop.

## 6. Solvers — `solvers/`

See [06_mars_solvers.md](06_mars_solvers.md). Available:
- CG (single-rank, distributed via `setHaloExchangeCallback`)
- BiCGSTAB (non-symmetric)
- GMRES(restart)
- Hypre PCG + BoomerAMG

CG uses cuSPARSE SpMV + cuBLAS dot + `MPI_Allreduce` (when
`ownedSize_ > 0`). **Watch**: `CUSPARSE_INDEX_32I` limits to ~2.1B
nonzeros (fix when scaling beyond).

## 7. AMR — `amr/`

See [05_mars_amr_module.md](05_mars_amr_module.md).
- `mars_amr_manager.hpp` — orchestration, InitTimings, AmrStats
- `mars_amr_octree.hpp` — Doerfler marking + 1-irregular
- `mars_amr_octree_native.hpp` — cstone-driven (recommended at scale)
- `mars_amr_hex_refine.hpp` — 1→8 hex with 19 new nodes
- `mars_amr_tet_refine.hpp` — Bey's red 1→8 tet
- `mars_amr_solution_transfer.hpp` — P0 trilinear by parentage
- `mars_amr_error_indicator.hpp` — gradient-based + Doerfler

## 8. Examples

| File | Element | Status |
|------|---------|--------|
| `mars_ex1_poisson.cu` | tet | working |
| `mars_cvfem_full.cu` | hex | working (Hypre AMG) |
| `mars_cvfem_graph.cu` | hex | working |
| `mars_amr_cvfem_graph.cu` | hex | working |
| `mars_amr_poisson.cu` | hex/tet | working |
| `mars_cvfem_poisson.cu` | hex | **⚠ segfault** — TODO #7 |
| `mars_ex_beam_tet.cu` | tet | working |

`PhaseTimer` pattern (RAII): cudaDeviceSynchronize on dtor + duration print.

## 9. Helpers

### extractConnPtrs — `domain.cu:932–936`
```cpp
template <typename KeyType, typename Tuple, std::size_t... Is>
void extractConnPtrs(const Tuple& t, const KeyType* ptrs[], std::index_sequence<Is...>) {
    ((ptrs[Is] = thrust::raw_pointer_cast(std::get<Is>(t).data())), ...);
}
```
Used at `domain.cu:973`.

### cudaCheckError
Wraps every kernel launch in domain.cu (lines 995, 1007, 1028, 1045, 1059,
1077, 1088, 1326, 1370, 1405, ...).

### Mesh I/O
- `read_mesh_binary` — `.float32/.double` SoA with separate connectivity
- `read_mfem_mesh` — MFEM `.mesh` ASCII
- `read_exodus_mesh` — Exodus II (if `MARS_HAVE_NETCDF`)

## 10. Connectivity Conventions
| Space | Type | Range | Sparse? |
|-------|------|-------|---------|
| Global SFC key | KeyType | [0, 2^32) or [0, 2^64) | yes |
| Dense local ID | KeyType | [0, nodeCount) | no |

Conversion via `mapSfcToLocalIdKernel` (`domain.cu:829–848`):
```cpp
KeyType local_id = cub::LowerBound(sorted_sfc, num_nodes, sfc_key);
assert(sorted_sfc[local_id] == sfc_key);
```

Naming: `d_*` device, `h_*` host, `*Sfc / *_keys` global SFC, `*_local_ids`
dense local, `nodeToDof` indirection map.

## 11. Identified Bugs (current state)

### Bug 1: NaN sentinel changed to T(0) — `domain.hpp:794`
```cpp
sbuf[i] = (idx >= 0) ? arr[idx] : T(0);
```
Comments at L754–757 state pack should use NaN for nodes this rank doesn't
own. Changing to `T(0)` corrupts ghost values that should remain owned-rank
data — zero is indistinguishable from real zero.
**Action**: revert to NaN sentinel with NaN-aware unpack, OR redesign so
unpack is index-driven (only writes nodes the recv list says belong to peer).

### Bug 2: Race on mutable buffer resize — `domain.hpp:776–780`
```cpp
mutable DeviceVector<RealType> sendBuf_, recvBuf_;
...
if (topo.sendBuf_.size() != sendTotal) topo.sendBuf_.resize(sendTotal);
if (topo.recvBuf_.size() != recvTotal) topo.recvBuf_.resize(recvTotal);
```
Concurrent calls to `exchangeNodeHalo()` race the resize. No mutex.
**Action**: pre-size during NodeHaloTopology construction OR add mutex.

### Bug 3: Broken device-resident NodeHaloTopology rewrite (current state)
Mentioned in compact summary. Current code at `domain.cu:1214–1725` is the
host-driven version. An earlier device-resident attempt was abandoned due
to NaN/non-deterministic results.

### Bug 4: Corner-shared double ownership (minor)
~9 nodes on 4-rank cube remain doubly-owned. NodeHaloTopology Phase 5
tiebreaker handles this in working version.

### Bug 5: Lazy field initialization order
`getNodeX/Y/Z` checks `if (originalCoords_) return ... else
ensureCoordinateCache()`. If `box_` is mutated post-construction, cached
coords go stale. Rare but possible.

### Bug 6: mars_cvfem_poisson.cu segfault — TODO #7
Likely uninitialized adjacency access or coordinate caching path. Needs
investigation; deferred.

## 12. Build System

CMake target: `mars_unstructured` (STATIC). Sources: `domain.cu` + `domain_cuda_impl.cpp`.

Options:
- `MARS_ENABLE_CUDA` — gates nvcc compilation
- `MARS_ENABLE_HIP` — AMD path; mutually exclusive with CUDA
- `MARS_ENABLE_UNSTRUCTURED` — this backend
- `MARS_ENABLE_FEM_EXAMPLES` — examples
- `MARS_ENABLE_TESTS` — test harness
- `MARS_ENABLE_MPI` — required for distributed
- `MARS_ENABLE_DISTRIBUTED_BACKEND` — distributed (Kokkos SFC) backend

cstone fetched via `FetchContent`; available at
`_deps/cornerstone_fetch-src/include/cstone/`. Patches applied via sed
(see [01_cstone_overview.md](01_cstone_overview.md) §11).

CUDA flags: `--diag-suppress=20281`, `CUDA_SEPARABLE_COMPILATION ON`,
`-Xnvlink=--suppress-stack-size-warning`. Architectures `80;90` typical.

ctest: `ctest -R marsReadMeshBinary` and similar filters.

## 13. Scaling Notes (10^9 elements target)

Per-rank at 10^9 / 1024 ranks = ~10^6 elements:
- Connectivity: 8 × 8B = 64 B/elem → 64 MB
- Coordinates: ~192 B/elem → 192 MB
- CSR adjacency: ~192 MB
- **Total ~500 MB/rank** (fits H100 80 GB easily)

MPI:
- Halo exchange: O(√(elem/rank)) boundary nodes → ~10k floats/iter, 40 KB
- Dominant: assembly (~4.6 ms for 10^6 hex with tensor kernel on GH200)

cstone partitioning: SFC sort gives spatial locality → bounded peer count
and halo size → linear scaling.
