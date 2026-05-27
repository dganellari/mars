#pragma once
// Mesh-level periodic vertex collapse (MFEM / Nek5000 style).
//
// Identifies periodic image vertices on opposite faces of a box-shaped domain,
// unions them into one logical vertex with a union-find, and rewrites the
// element-to-node connectivity to use the merged vertex IDs. The merged
// mesh is then handed to cstone for SFC partitioning -- because the
// connectivity already has periodic-paired vertices identified, cstone's
// existing shared-vertex / halo machinery handles cross-rank periodic
// communication with no special periodic logic at the solver layer.
//
// References:
//   * MFEM Mesh::CreatePeriodicVertexMapping (mesh.cpp:6048-6217) -- the
//     canonical implementation. Uses a KD-tree; we use a flat sweep + sort
//     since boundary vertex counts are O(N^(2/3)) and small.
//   * Nek5000 setvert3d (core/navier8.f) -- same pattern with glo_num
//     collisions.
//
// API:
//   applyPeriodicCollapse(h_x, h_y, h_z, h_conn, box_lo, box_hi, eps)
//       In-place modifies h_conn so that every periodic-image-pair vertex
//       maps to its primary; does NOT compact the coordinate arrays
//       (cstone handles unused vertices during sync). Returns the number
//       of merged pairs (for reporting).

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <tuple>
#include <vector>
#include <mpi.h>

namespace mars { namespace fem {

// Union-find with path compression, in a plain std::vector for portability.
namespace detail {
struct UnionFind {
    std::vector<int> parent;
    explicit UnionFind(int n) : parent(n) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }
    int find(int i) {
        while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
        return i;
    }
    // Make i a replica of j. Caller ensures j is a primary.
    void unite(int i, int j) {
        int ri = find(i), rj = find(j);
        if (ri == rj) return;
        // Smaller-index wins so we get a deterministic primary.
        if (ri < rj) parent[rj] = ri; else parent[ri] = rj;
    }
};
} // namespace detail

// Apply periodic vertex collapse to a hex8 mesh. Modifies h_conn in place;
// h_x / h_y / h_z are read only (we don't compact them -- cstone removes
// unused nodes during its own dedup pass).
//
// boxLo / boxHi define the periodic box (inclusive). Vertices on the
// "max" face along each axis are merged with their image on the "min" face.
//
// epsRel is the matching tolerance as a fraction of the box width (1e-6
// is a reasonable default; mesh nodes are point-exact in binary format).
//
// Returns the number of vertex IDs that were demoted to replicas of a
// smaller-index primary. Equals the number of distinct periodic merges.
template<typename RealType, typename HostConnT>
inline std::size_t applyPeriodicCollapse(const std::vector<RealType>& h_x,
                                         const std::vector<RealType>& h_y,
                                         const std::vector<RealType>& h_z,
                                         HostConnT& h_conn,
                                         RealType boxLo, RealType boxHi,
                                         bool periodicX = true,
                                         bool periodicY = true,
                                         bool periodicZ = true,
                                         RealType epsRel = RealType(1e-6))
{
    const std::size_t numNodes = h_x.size();
    if (numNodes == 0) return 0;

    const RealType L = boxHi - boxLo;
    const RealType eps = epsRel * L;

    detail::UnionFind uf{int(numNodes)};

    // Helper: for one axis (axis = 0,1,2), gather (perpKey, nodeIndex) for
    // all nodes on the min-face and the max-face, sort each list by perpKey,
    // then sweep-match. Quantize perp coords to int64 buckets so floating
    // noise doesn't break the match.
    auto quantize = [&](RealType a, RealType b) -> uint64_t {
        // Two perpendicular coords, each in [boxLo, boxHi]. Map to a 32-bit
        // integer and concatenate into 64. Resolution = L / 2^31 -- well
        // below mesh node spacing for any reasonable mesh.
        constexpr uint64_t MASK = 0x7FFFFFFFu;
        uint64_t ai = uint64_t(((a - boxLo) / L) * RealType(double(MASK))) & MASK;
        uint64_t bi = uint64_t(((b - boxLo) / L) * RealType(double(MASK))) & MASK;
        return (ai << 32) | bi;
    };

    auto matchAxis = [&](int axis) {
        std::vector<std::pair<uint64_t, int>> minFace, maxFace;
        minFace.reserve(numNodes / 32);
        maxFace.reserve(numNodes / 32);

        for (std::size_t i = 0; i < numNodes; ++i) {
            RealType ci = (axis == 0) ? h_x[i] : (axis == 1) ? h_y[i] : h_z[i];
            // Use the OTHER two coords as the perpendicular key.
            RealType pa = (axis == 0) ? h_y[i] : (axis == 1) ? h_x[i] : h_x[i];
            RealType pb = (axis == 0) ? h_z[i] : (axis == 1) ? h_z[i] : h_y[i];
            if (std::fabs(ci - boxLo) < eps) {
                minFace.emplace_back(quantize(pa, pb), int(i));
            } else if (std::fabs(ci - boxHi) < eps) {
                maxFace.emplace_back(quantize(pa, pb), int(i));
            }
        }

        std::sort(minFace.begin(), minFace.end());
        std::sort(maxFace.begin(), maxFace.end());

        // Linear merge: each max-face vertex looks up its perp key in the
        // sorted min-face list. Equal-key entries are periodic images.
        std::size_t merged = 0;
        std::size_t i = 0, j = 0;
        while (i < maxFace.size() && j < minFace.size()) {
            if      (maxFace[i].first < minFace[j].first) ++i;
            else if (maxFace[i].first > minFace[j].first) ++j;
            else {
                // Found pair. union: max-face node becomes replica of min-face.
                uf.unite(maxFace[i].second, minFace[j].second);
                ++merged;
                ++i; ++j;
            }
        }
        return merged;
    };

    std::size_t total = 0;
    if (periodicX) total += matchAxis(0);
    if (periodicY) total += matchAxis(1);
    if (periodicZ) total += matchAxis(2);

    // Rewrite connectivity. h_conn is a tuple of 8 std::vector<KeyType>
    // (one per corner of a hex8). After this loop, every element node ID
    // points at its periodic primary.
    constexpr int NCORNERS = std::tuple_size<HostConnT>::value;
    static_assert(NCORNERS == 8 || NCORNERS == 4,
                  "applyPeriodicCollapse: only hex8 / tet4 supported");
    const std::size_t numElems = std::get<0>(h_conn).size();

    auto remapColumn = [&](auto& col) {
        for (std::size_t e = 0; e < numElems; ++e) {
            col[e] = static_cast<typename std::remove_reference<decltype(col[e])>::type>(
                uf.find(int(col[e])));
        }
    };

    // For each corner column, rewrite IDs.
    if constexpr (NCORNERS == 8) {
        remapColumn(std::get<0>(h_conn));
        remapColumn(std::get<1>(h_conn));
        remapColumn(std::get<2>(h_conn));
        remapColumn(std::get<3>(h_conn));
        remapColumn(std::get<4>(h_conn));
        remapColumn(std::get<5>(h_conn));
        remapColumn(std::get<6>(h_conn));
        remapColumn(std::get<7>(h_conn));
    } else if constexpr (NCORNERS == 4) {
        remapColumn(std::get<0>(h_conn));
        remapColumn(std::get<1>(h_conn));
        remapColumn(std::get<2>(h_conn));
        remapColumn(std::get<3>(h_conn));
    }

    return total;
}

// =============================================================================
// MFEM-style mesh-level periodic identification for MARS / cstone.
//
// **Strategy chosen after tracing the MARS data flow**:
//   - cstone dedups vertices by SFC key, which is a pure coord->key
//     quantization (cstone::sfc3D is NOT periodic-aware).
//   - Two vertices at the same physical position get the same SFC key.
//   - cstone's host-fallback peer discovery (MARS_NODEHALO_V2=0) detects
//     cross-rank shared vertices by MPI_Allgatherv of owned-node SFC keys.
//   - Therefore: collapse periodic pairs by rewriting the SLAVE'S
//     COORDINATES to match the master's. Then cstone naturally identifies
//     the two physical locations as one vertex (locally via thrust::unique
//     on SFC keys, cross-rank via the Allgatherv key-match in
//     buildNodeHaloTopologyHostPath at domain.cu:1174-1198).
//
// This means **no MPI is needed at collapse time** -- each rank's slaves
// just shift their own (x,y,z) by the periodicity vector. The cross-rank
// pairing happens later, automatically, in cstone's existing halo build.
//
// Inputs:
//   h_x, h_y, h_z   -- per-rank local node coordinates (modified in-place)
//   boxLo, boxHi    -- periodic domain extent (same on every rank)
//   periodicX/Y/Z   -- which axes are periodic
//   epsRel          -- face-detection tolerance, fraction of box width
//
// Returns: number of local vertices whose coords were shifted (one per
// max-face boundary node per active axis).
//
// h_conn is NOT modified. Element connectivity still references the same
// local indices; what changes is that two distinct local indices now
// resolve to the same SFC key and get merged by cstone's sync().
// =============================================================================
template<typename RealType>
inline std::size_t collapsePeriodicCoordinatesInPlace(
    std::vector<RealType>& h_x,
    std::vector<RealType>& h_y,
    std::vector<RealType>& h_z,
    RealType boxLo, RealType boxHi,
    bool periodicX = true,
    bool periodicY = true,
    bool periodicZ = true,
    RealType epsRel = RealType(1e-6))
{
    const std::size_t numNodes = h_x.size();
    if (numNodes == 0) return 0;

    const RealType L   = boxHi - boxLo;
    const RealType eps = epsRel * L;

    std::size_t shifted = 0;

    // For each axis, find local nodes on the max-face (x=boxHi) and shift
    // their corresponding coord to boxLo. Since the master at boxLo has
    // identical perpendicular coords, the resulting vertex has identical
    // SFC key as the master's image (whether the master lives locally or
    // on another rank).
    if (periodicX) {
        for (std::size_t i = 0; i < numNodes; ++i) {
            if (std::fabs(h_x[i] - boxHi) < eps) {
                h_x[i] = boxLo;
                ++shifted;
            }
        }
    }
    if (periodicY) {
        for (std::size_t i = 0; i < numNodes; ++i) {
            if (std::fabs(h_y[i] - boxHi) < eps) {
                h_y[i] = boxLo;
                ++shifted;
            }
        }
    }
    if (periodicZ) {
        for (std::size_t i = 0; i < numNodes; ++i) {
            if (std::fabs(h_z[i] - boxHi) < eps) {
                h_z[i] = boxLo;
                ++shifted;
            }
        }
    }

    return shifted;
}

// =============================================================================
// MPI-parallel periodic vertex collapse.
//
// Each rank holds a per-rank slice of the mesh: local node coordinates
// h_x/h_y/h_z and a connectivity h_conn where each entry is a GLOBAL node
// index (the binary format stores global IDs). After this call, h_conn
// has periodic-image global IDs rewritten to their primary -- consistently
// across all ranks -- so handing the result to cstone's ElementDomain
// gives every periodic pair the same logical vertex.
//
// Algorithm:
//   1. Each rank scans its local nodes, gathers (perp_key, my_global_id)
//      pairs for boundary nodes on the min-face along each periodic axis.
//   2. MPI_Allgatherv combines per-rank min-face pair lists into a sorted
//      global master table.
//   3. Each rank scans its local max-face nodes, binary-searches the
//      global master table, and records (own_global_id -> master_global_id)
//      mappings.
//   4. MPI_Allgatherv combines the mappings; path-compress so corner nodes
//      that are simultaneously max-faces on 2 or 3 axes resolve to the
//      single lowest-index primary.
//   5. Each rank rewrites its h_conn through the global mapping.
//
// `localToGlobal` is the per-rank table mapping local node index -> global
// node ID (the file-format ID). Pass an empty vector if h_conn already
// holds local indices (single-rank case).
//
// Returns the total number of distinct global IDs that were merged into
// a smaller primary (across all ranks).
// =============================================================================
template<typename RealType, typename HostConnT>
inline std::size_t applyPeriodicCollapseMpi(const std::vector<RealType>& h_x,
                                            const std::vector<RealType>& h_y,
                                            const std::vector<RealType>& h_z,
                                            const std::vector<std::int64_t>& localToGlobal,
                                            HostConnT& h_conn,
                                            RealType boxLo, RealType boxHi,
                                            bool periodicX,
                                            bool periodicY,
                                            bool periodicZ,
                                            MPI_Comm comm,
                                            RealType epsRel = RealType(1e-6))
{
    const std::size_t numLocalNodes = h_x.size();
    const RealType L   = boxHi - boxLo;
    const RealType eps = epsRel * L;

    int rank = 0, numRanks = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &numRanks);

    // Quantize two perpendicular coords into one 64-bit key so floating
    // noise doesn't break matches across ranks.
    auto quantize = [&](RealType a, RealType b) -> std::uint64_t {
        constexpr std::uint64_t MASK = 0x7FFFFFFFu;
        std::uint64_t ai = std::uint64_t(((a - boxLo) / L) * RealType(double(MASK))) & MASK;
        std::uint64_t bi = std::uint64_t(((b - boxLo) / L) * RealType(double(MASK))) & MASK;
        return (ai << 32) | bi;
    };

    // For each axis, build a global master table via Allgatherv of
    // (perp_key, master_global_id) pairs from each rank's min-face nodes,
    // then each rank rewrites its max-face nodes' global IDs.
    //
    // Net change in h_conn after all three axes: every max-face vertex's
    // global ID is replaced by its periodic primary's global ID. Path
    // compression for corners (vertex on 2 or 3 max-faces) is handled by
    // running the axes in sequence and using the *latest* h_conn IDs each
    // pass, with the global master table re-gathered each axis (cheap;
    // boundary count is O(N^(2/3))).

    auto matchAxisGlobal = [&](int axis) -> std::size_t {
        // --- 1. Collect min-face masters on this rank ---
        struct Pair { std::uint64_t key; std::int64_t gid; };
        std::vector<Pair> localMasters;
        std::vector<Pair> localSlaves;   // (perp_key, this rank's own global_id for the max-face vertex)
        localMasters.reserve(numLocalNodes / 32);
        localSlaves.reserve(numLocalNodes / 32);

        for (std::size_t i = 0; i < numLocalNodes; ++i) {
            RealType ci = (axis == 0) ? h_x[i] : (axis == 1) ? h_y[i] : h_z[i];
            RealType pa = (axis == 0) ? h_y[i] : (axis == 1) ? h_x[i] : h_x[i];
            RealType pb = (axis == 0) ? h_z[i] : (axis == 1) ? h_z[i] : h_y[i];
            std::int64_t gid = localToGlobal.empty() ? std::int64_t(i)
                                                     : localToGlobal[i];
            if (std::fabs(ci - boxLo) < eps) {
                localMasters.push_back({quantize(pa, pb), gid});
            } else if (std::fabs(ci - boxHi) < eps) {
                localSlaves.push_back({quantize(pa, pb), gid});
            }
        }

        // --- 2. Allgatherv of (perp_key, master_global_id) ---
        int nLocal = int(localMasters.size());
        std::vector<int> recvCounts(numRanks), recvDispls(numRanks + 1, 0);
        MPI_Allgather(&nLocal, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, comm);
        for (int r = 0; r < numRanks; ++r)
            recvDispls[r + 1] = recvDispls[r] + recvCounts[r];
        int globalCount = recvDispls[numRanks];

        // Pack Pair as two parallel arrays for portable MPI types.
        std::vector<std::uint64_t> localKeys(nLocal);
        std::vector<std::int64_t>  localGids(nLocal);
        for (int i = 0; i < nLocal; ++i) {
            localKeys[i] = localMasters[i].key;
            localGids[i] = localMasters[i].gid;
        }

        std::vector<std::uint64_t> globalKeys(globalCount);
        std::vector<std::int64_t>  globalGids(globalCount);
        MPI_Allgatherv(localKeys.data(), nLocal, MPI_UINT64_T,
                       globalKeys.data(), recvCounts.data(), recvDispls.data(),
                       MPI_UINT64_T, comm);
        MPI_Allgatherv(localGids.data(), nLocal, MPI_INT64_T,
                       globalGids.data(), recvCounts.data(), recvDispls.data(),
                       MPI_INT64_T, comm);

        // --- 3. Sort the global master table by perp_key ---
        std::vector<int> order(globalCount);
        for (int i = 0; i < globalCount; ++i) order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](int a, int b) { return globalKeys[a] < globalKeys[b]; });
        std::vector<std::uint64_t> sortedKeys(globalCount);
        std::vector<std::int64_t>  sortedGids(globalCount);
        for (int i = 0; i < globalCount; ++i) {
            sortedKeys[i] = globalKeys[order[i]];
            sortedGids[i] = globalGids[order[i]];
        }

        // --- 4. For each local slave, binary-search and record mapping ---
        // mappings[i] = (own_gid, master_gid)
        std::vector<std::pair<std::int64_t, std::int64_t>> mappings;
        mappings.reserve(localSlaves.size());
        for (const auto& s : localSlaves) {
            auto it = std::lower_bound(sortedKeys.begin(), sortedKeys.end(), s.key);
            if (it != sortedKeys.end() && *it == s.key) {
                std::int64_t masterGid = sortedGids[it - sortedKeys.begin()];
                // Make the smaller-index gid the primary, like the union-find
                // version does locally. This keeps the result independent
                // of which rank owned which side of the pair.
                std::int64_t primary = std::min(s.gid, masterGid);
                std::int64_t replica = std::max(s.gid, masterGid);
                if (primary != replica)
                    mappings.emplace_back(replica, primary);
            }
        }

        // --- 5. Allgatherv the mappings so every rank knows the full
        //         (replica_gid -> primary_gid) table for this axis ---
        int nMaps = int(mappings.size());
        std::vector<int> mapCounts(numRanks), mapDispls(numRanks + 1, 0);
        MPI_Allgather(&nMaps, 1, MPI_INT, mapCounts.data(), 1, MPI_INT, comm);
        for (int r = 0; r < numRanks; ++r)
            mapDispls[r + 1] = mapDispls[r] + mapCounts[r];
        int globalMaps = mapDispls[numRanks];

        std::vector<std::int64_t> myRep(nMaps), myPri(nMaps);
        for (int i = 0; i < nMaps; ++i) {
            myRep[i] = mappings[i].first;
            myPri[i] = mappings[i].second;
        }
        std::vector<std::int64_t> gRep(globalMaps), gPri(globalMaps);
        MPI_Allgatherv(myRep.data(), nMaps, MPI_INT64_T,
                       gRep.data(), mapCounts.data(), mapDispls.data(),
                       MPI_INT64_T, comm);
        MPI_Allgatherv(myPri.data(), nMaps, MPI_INT64_T,
                       gPri.data(), mapCounts.data(), mapDispls.data(),
                       MPI_INT64_T, comm);

        // --- 6. Path compression: build a global gid -> primary gid map.
        // Use std::sort + linear sweep instead of std::map for speed.
        std::vector<std::pair<std::int64_t, std::int64_t>> remap(globalMaps);
        for (int i = 0; i < globalMaps; ++i) remap[i] = {gRep[i], gPri[i]};
        std::sort(remap.begin(), remap.end());
        // Repeatedly resolve chains: replica -> primary -> primary's primary -> ...
        bool changed = true;
        int pcIters = 0;
        while (changed && pcIters < 8) {
            changed = false;
            for (auto& m : remap) {
                auto it = std::lower_bound(remap.begin(), remap.end(),
                                           std::make_pair(m.second, std::int64_t(-1)));
                if (it != remap.end() && it->first == m.second) {
                    if (it->second != m.second) {
                        m.second = it->second;
                        changed = true;
                    }
                }
            }
            ++pcIters;
        }

        // --- 7. Apply remap to local h_conn entries (which hold global IDs) ---
        const std::size_t numElems = std::get<0>(h_conn).size();
        constexpr int NCORNERS = std::tuple_size<HostConnT>::value;
        auto remapCol = [&](auto& col) {
            for (std::size_t e = 0; e < numElems; ++e) {
                auto it = std::lower_bound(remap.begin(), remap.end(),
                                           std::make_pair(std::int64_t(col[e]),
                                                          std::int64_t(-1)));
                if (it != remap.end() && it->first == std::int64_t(col[e])) {
                    col[e] = static_cast<typename std::remove_reference<decltype(col[e])>::type>(
                        it->second);
                }
            }
        };
        if constexpr (NCORNERS == 8) {
            remapCol(std::get<0>(h_conn));
            remapCol(std::get<1>(h_conn));
            remapCol(std::get<2>(h_conn));
            remapCol(std::get<3>(h_conn));
            remapCol(std::get<4>(h_conn));
            remapCol(std::get<5>(h_conn));
            remapCol(std::get<6>(h_conn));
            remapCol(std::get<7>(h_conn));
        }

        return std::size_t(globalMaps);
    };

    std::size_t total = 0;
    if (periodicX) total += matchAxisGlobal(0);
    if (periodicY) total += matchAxisGlobal(1);
    if (periodicZ) total += matchAxisGlobal(2);

    return total;
}

}} // namespace mars::fem
