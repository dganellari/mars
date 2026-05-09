# Cornerstone-Octree Reference for MARS

Header-only + CUDA library at `/Users/gandanie/scratch/santis/mars/cornerstone-octree/`.
This document is the API-surface map. For algorithmic depth (Hilbert, focus
tree, sync data flow), see [04_cstone_focus_tree.md](04_cstone_focus_tree.md).

## 1. Repo Layout

| Directory | Purpose |
|-----------|---------|
| `cstone/cuda/` | GPU abstractions: `DeviceVector` (PIMPL), errorcheck, thrust wrappers |
| `cstone/domain/` | `Domain` class, `GlobalAssignment` SFC splitting, buffer layouts, halo book-keeping |
| `cstone/fields/` | Particle field utilities |
| `cstone/focus/` | Focused octree on GPU: construction, convergence, MAC updates |
| `cstone/halos/` | Halo exchange: `Halos`, CPU `haloexchange`, GPU `haloExchangeGpu` |
| `cstone/primitives/` | sort, gather/scatter, MPI wrappers, math, CUB integrations |
| `cstone/sfc/` | Hilbert/Morton encoding (32/64-bit), `Box`, key ops |
| `cstone/traversal/` | Dual-tree MAC, `findPeersMac`, collision/overlap detection |
| `cstone/tree/` | `Octree`, `OctreeView`, B-tree leaves, `TreeNodeIndex`, tree construction |
| `cstone/util/` | Allocators, tuple metaprogramming, strong types, `array` container |

No `src/`; everything is in headers / `.cuh`.

## 2. Core Types

### Index types — `tree/definitions.h`
| Type | Where | Notes |
|------|-------|-------|
| `TreeNodeIndex` | definitions.h:25 | `int`, octree node index |
| `LocalIndex` | definitions.h:27 | `unsigned`, max ~4B particles per rank |
| `maxTreeLevel<unsigned>` | definitions.h:51 | 10 levels (32-bit) |
| `maxTreeLevel<unsigned long long>` | definitions.h:54 | 21 levels (64-bit) |

### SFC keys — `sfc/sfc.hpp:27–112`
```cpp
template<class IntegerType> using MortonKey  = StrongType<IntegerType, struct MortonKeyTag>;
template<class IntegerType> using HilbertKey = StrongType<IntegerType, struct HilbertKeyTag>;
template<class IntegerType> using SfcKind    = HilbertKey<IntegerType>;  // default
```
Instantiations: `HilbertKey<unsigned>`, `HilbertKey<unsigned long long>`,
`MortonKey<...>`. Encode/decode at `sfc.hpp:142–194`.

### IndexPair / IndexRanges — `domain/index_ranges.hpp:29–107`
```cpp
template<class T> struct IndexPair : public std::tuple<T,T> {
    T start() const; T end() const; T count() const;
};
using TreeIndexPair = IndexPair<TreeNodeIndex>;

class IndexRanges<Index> {           // CSR-like: multiple (start,end) ranges
    void addRange(IndexType lo, IndexType hi);
    IndexType rangeStart(size_t i) const;
    IndexType rangeEnd(size_t i) const;
    std::size_t count(size_t i) const;
    std::size_t totalCount() const;
};
using SendManifest = IndexRanges<LocalIndex>;
using RecvList     = std::vector<IndexPair<LocalIndex>>;       // 1 contiguous range per rank
using SendList     = std::vector<SendManifest>;                // multiple ranges per rank
```

**Critical asymmetry**: `incomingHaloIndices_[peer]` is a **single
contiguous range** (RecvList element). `outgoingHaloIndices_[peer]` is a
SendManifest, possibly **multiple disjoint ranges**.

### BufferDescription — `domain/buffer_description.hpp:29–45`
```cpp
struct BufferDescription { LocalIndex start; LocalIndex end; LocalIndex size; };
```
Layout: `|-- halos --|-- assigned --|-- halos --|`
                    `0         start           end         size`

### Box — `sfc/box.hpp:94–149`
Float bounding box with per-axis `BoundaryType { open, periodic, fixed }`
(`sfc.hpp:80–85`). Accessors: `xmin/xmax/.../lx/ly/lz/ilx/ily/ilz`.

### SfcAssignment — `domain/domaindecomp.hpp:57–110`
```cpp
template<class KeyType> class SfcAssignment {
    KeyType operator[](int rank) const;
    int findRank(KeyType key) const;            // binary search
    LocalIndex totalCount(int rank) const;
    int numRanks() const;
    std::span<const TreeNodeIndex> numNodesPerRank() const;
    std::span<const TreeNodeIndex> treeOffsetsConst() const;
};
```

### DeviceVector — `cuda/device_vector.h:32–72`
PIMPL wrapper hiding thrust headers. Explicit instantiations:
- arithmetic: `char, uint8_t, int, unsigned, uint64_t, float, double`
- arrays: `util::array<int,2/3>`, `util::array<unsigned,1/2>`,
  `util::array<float/double,3/4>`

`rawPtr(deviceVector) -> T*`. **`double` is instantiated**, but for
gatherRanges and exchange paths, support depends on the kernel's own
instantiations (see Section 4).

## 3. Domain Class — `domain/domain.hpp:38–641`

### Constructor — `domain.hpp:63–81`
```cpp
Domain(int rank, int nRanks, unsigned bucketSize, unsigned bucketSizeFocus,
       float theta, const Box<T>& box = Box<T>{0, 1});
```
Asserts `bucketSize >= bucketSizeFocus`.

### Members
- `BufferDescription prevBufDesc_, bufDesc_;`
- `FocusedOctree<KeyType, T, Accelerator> focusTree_;`
- `AccVector<LocalIndex> layoutAcc_;` (GPU/CPU layout)
- `std::vector<LocalIndex> layout_;` (CPU mirror)
- `GlobalAssignment<KeyType, T, Accelerator> global_;`
- `Halos<KeyType, Accelerator> halos_{myRank_};`

### sync() — `domain.hpp:164–206`
```cpp
template<class KeyVec, class VectorX, class VectorH, class... Vectors1, class... Vectors2>
void sync(KeyVec& particleKeys, VectorX& x, VectorX& y, VectorX& z, VectorH& h,
          std::tuple<Vectors1&...> particleProperties,
          std::tuple<Vectors2&...> scratchBuffers);
```
**Preconditions** (domain.hpp:95–109): array sizes equal `localNParticles_`
from previous call (any size on first call); ≥3 scratch buffers.

**Postconditions**: arrays resized to `nParticles + nHalos`; `x,y,z,h` at
`[startIndex():endIndex()]` filled and **SFC-sorted**; halo property values
**undefined** until `exchangeHalos(property)`.

**Step sequence** (domain.hpp:150–162, 174–206):
1. Compute global bbox
2. Build global octree from assigned-range particles
3. Compute per-node max interaction radius
4. Assign leaves to ranks via SFC split
5. Discover halos
6. Compute layout
7. Resize arrays to assigned + halo
8. MPI exchange coords + properties
9. Sort to SFC order
10. Halo-exchange remaining properties

### syncGrav() — `domain.hpp:208–282`
Adds mass `m`; iterates focus-tree refinement using expansion centers
(center-of-mass) until convergence MPI_Allreduce reports 0.

### Public accessors — `domain.hpp:333–390`
```cpp
const auto& incomingHaloIndices() const;   // forward to halos_
const auto& outgoingHaloIndices() const;
LocalIndex startIndex() const;
LocalIndex endIndex() const;
LocalIndex nParticles() const;
LocalIndex nParticlesWithHalos() const;
OctreeView<const KeyType> globalTree() const;
const FocusedOctree<KeyType, T, Accelerator>& focusTree() const;
TreeNodeIndex startCell() const;
TreeNodeIndex endCell() const;
std::span<const LocalIndex> layout() const;
const Box<T>& box() const;
KeyType assignmentStart() const;
```

### exchangeHalos() / reapplySync() — `domain.hpp:291–331`
Repeat previous halo pattern for new arrays; uses `mutable haloEpoch_` for
unique MPI tags.

## 4. Halo Subsystem

### Halos<KeyType, Accelerator> — `halos/halos.hpp:38–125`
```cpp
int exchangeRequests(std::span<const KeyType> leaves,
                     std::span<const TreeIndexPair> assignment,
                     std::span<const int> peers,
                     std::span<const LocalIndex> layout);

template<class Scratch1, class Scratch2, class... Vectors>
void exchangeHalos(std::tuple<Vectors&...> arrays,
                   Scratch1& sendBuffer, Scratch2& receiveBuffer) const;

const RecvList& incomingHaloIndices() const;   // per-peer single range
const SendList& outgoingHaloIndices() const;   // per-peer multiple ranges

mutable int haloEpoch_{0};   // <-- requires MPI collective between calls
```

### CPU path — `halos/exchange_halos.hpp:26–93`
```cpp
template<class... Arrays>
void haloexchange(int epoch, const RecvList&, const SendList&, Arrays...);
```
Pack to uint64_t-aligned host buffers, `MPI_Isend` per peer, `MPI_Recv` with
`MPI_ANY_SOURCE`. Tag = `P2pTags::haloExchange + epoch`.

### GPU path — `halos/exchange_halos_gpu.cuh:34–120`
```cpp
template<class DevVec1, class DevVec2, class... Arrays>
void haloExchangeGpu(int epoch, const RecvList&, const SendList&,
                     DevVec1& sendScratch, DevVec2& recvScratch, Arrays... arrays);
```
- `TransferType = uint64_t` for alignment + 32-bit count workaround
- `gatherRanges(d_rangeScan, d_rangeOffsets, numRanges, src, buffer, bufferSize)`
  packs scattered ranges into contiguous device buffer (`halos/gather_halos_gpu.h:22`)
- `mpiSendGpuDirect` / `mpiRecvGpuDirect` from `primitives/mpi_cuda.cuh:31–74`
  branch on `CSTONE_HAVE_GPU_AWARE_MPI` at compile time

### Type instantiation gaps
`gatherRanges` is explicitly instantiated for: `float, double, int,
unsigned, uint64_t, char, uint8_t`, and small `util::array<...,3>` variants.
**Custom struct types require an explicit kernel instantiation in a new .cu**
or byte-level packing.

### haloEpoch_ pitfall
`mutable int haloEpoch_{0}` is **not thread-safe**. Multiple threads calling
exchangeHalos concurrently corrupt the counter. Also: between two
`exchangeHalos()` calls the user **must execute an MPI collective** (e.g.
`MPI_Barrier`), or tag mismatches deadlock.

## 5. SFC Machinery — `sfc/sfc.hpp`

```cpp
KeyType sfc3D<KeyType, T>(T x, T y, T z, const Box<T>& box);   // L142–178
util::tuple<unsigned,unsigned,unsigned> decodeSfc<KeyType>(KeyType k);  // L181–194
void computeSfcKeys<T,KeyType>(const T* x, const T* y, const T* z,
                               KeyType* keys, size_t n, const Box<T>& box);  // L268
```
`computeSfcKeys` skips entries equal to `removeKey<KeyType>::value`.

## 6. Tree / Octree
- `TreeNodeIndex` = `int` index into octree node arrays
- Cornerstone octree = SFC-key array marking leaf starts (3 invariants in `csarray.hpp:19–33`)
- Full octree (`tree/octree.hpp`) stores internal+leaf with Warren-Salmon
  placeholder bits encoding level
- See [04_cstone_focus_tree.md](04_cstone_focus_tree.md) §2–3 for full detail

## 7. MPI Utilities — `primitives/mpi_wrappers.hpp` + `mpi_cuda.cuh`
- `MpiType<T>` mapping (L25–94) for double, float, char, uchar, short, ushort,
  int, unsigned, long, ulong, ulonglong
- `mpiSendAsync<T>(data, count, rank, tag, requests)` — `MPI_Isend`
- `mpiRecvAsync<T>(...)` — `MPI_Irecv`
- `mpiAllgatherv` (L202–226) — handles struct/array types via scaled counts
- GPU variants in `mpi_cuda.cuh`; conditional on
  `CSTONE_HAVE_GPU_AWARE_MPI` (L24–28)

## 8. Reusable Primitives Table (for MARS)

| Function | File:Line | Purpose |
|----------|-----------|---------|
| `sortByKey<useGpu>` | primitives_acc.hpp:101 | Sort + permutation, GPU/CPU dual path |
| `gatherRanges` | halos/gather_halos_gpu.h:22 | Pack disjoint device ranges to contiguous buffer |
| `gatherAcc<useGpu>` | primitives_acc.hpp:58 | Gather: `dst[i] = src[ord[i]]` |
| `scatterAcc<useGpu>` | primitives_acc.hpp:65 | Scatter: `dst[ord[i]] = src[i]` |
| `sequence<useGpu>` | primitives_acc.hpp:92 | Iota fill |
| `sfc3D` / `decodeSfc` / `computeSfcKeys` | sfc/sfc.hpp | SFC encode/decode |
| `findPeersMac<T,KeyType>` | traversal/peers.hpp:46 | Dual-tree MAC peer discovery |
| `Domain::sync` | domain/domain.hpp:164 | Full sync |
| `Domain::exchangeHalos` | domain/domain.hpp:327 | Repeat halo for new arrays |
| `halos_.exchangeRequests` | halos/halos.hpp:60 | Negotiate halo index ranges |
| `haloExchangeGpu` | halos/exchange_halos_gpu.cuh:34 | GPU-direct MPI halo exchange |
| `mpiSendGpuDirect` | primitives/mpi_cuda.cuh:31 | GPU-aware send (or staging) |
| `mpiRecvGpuDirect` | primitives/mpi_cuda.cuh:62 | GPU-aware recv (or staging) |
| `DeviceVector<T>` | cuda/device_vector.h:32 | PIMPL device buffer |

## 9. Canonical Usage Example — `test/integration_mpi/exchange_halos_gpu.cpp:70–150`
```cpp
RecvList incomingHalos(numRanks);
SendList outgoingHalos(numRanks);
if (rank == 0) {
    outgoingHalos[1].addRange(0, 1);
    outgoingHalos[1].addRange(1, 3);
    incomingHalos[1] = {3, 10};
}
// ... mirror on rank 1 ...

DeviceVector<double> d_x = x;
DeviceVector<char> sendBuffer(7 * 24);
DeviceVector<char> receiveBuffer(7 * 24);

haloExchangeGpu(0, incomingHalos, outgoingHalos, sendBuffer, receiveBuffer,
                rawPtr(d_x), rawPtr(d_y), rawPtr(d_z));
```

## 10. Limitations / Gotchas
1. `haloEpoch_` mutable, not thread-safe; needs MPI collective between calls
2. Halo discovery is **particle-level**, not node-level. MARS builds node
   halos on top via cstone's element halo bookkeeping.
3. `DeviceVector` copy = deep copy. Prefer swap.
4. `gatherRanges` instantiation gaps for custom types
5. Particle array size invariant across syncs (use `removeKey` sentinel)
6. SFC functions assume global rectangular `Box`; coords outside clipped
7. MPI collective required between `exchangeHalos()` calls

## 11. Build / Fetch
| Flag | Effect |
|------|--------|
| `CSTONE_WITH_CUDA` (or USE_CUDA) | Enable CUDA paths |
| `CSTONE_WITH_GPU_AWARE_MPI` | Direct device-pointer MPI vs host staging |
| `CSTONE_WITH_HIP` | AMD path; mutually exclusive with CUDA |

MARS pulls cstone via `FetchContent_Declare`. Output at
`_deps/cornerstone_fetch-src/include/cstone/`. **Patches we apply** (sed-based,
preserving everything else): added accessors in `halos.hpp:111–112` and
forwarders in `domain.hpp` for `incomingHaloIndices()` /
`outgoingHaloIndices()`. Headers only — no link library.

C++ standards: host C++20, CUDA C++17, HIP C++20.

## 12. Key File Index
```
include/cstone/
  domain/domain.hpp                 L38–641     Domain orchestration
  domain/assignment.hpp             L36–210     GlobalAssignment
  domain/index_ranges.hpp           L28–107     IndexPair, SendList, RecvList
  domain/buffer_description.hpp     L29–45      BufferDescription
  domain/domaindecomp.hpp           L57–281     SfcAssignment, makeSfcAssignment
  halos/halos.hpp                   L38–125     Halos class (+ our patch L111–112)
  halos/exchange_halos.hpp          L26–93      CPU haloexchange
  halos/exchange_halos_gpu.cuh      L34–120     GPU-direct haloExchangeGpu
  halos/gather_halos_gpu.h          L21–28      gatherRanges signature
  sfc/sfc.hpp                       L27–277     sfc3D, decodeSfc, key types
  sfc/box.hpp                       L94–149     Box<T>
  tree/definitions.h                L19–116     index types, removeKey
  tree/octree.hpp                   L38–196     full octree
  primitives/gather.hpp             L29–148     CPU gather/scatter, SfcSorter
  primitives/primitives_acc.hpp     L27–109     unified GPU/CPU wrappers
  primitives/mpi_wrappers.hpp       L24–226     MPI type maps, Isend/Irecv
  primitives/mpi_cuda.cuh           L15–141     GPU-aware MPI
  cuda/device_vector.h              L26–102     PIMPL DeviceVector
  traversal/peers.hpp               L21–161     findPeersMac
test/integration_mpi/
  exchange_halos_gpu.cpp            L1–150      canonical usage
```
