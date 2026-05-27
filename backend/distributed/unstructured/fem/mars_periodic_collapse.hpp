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

}} // namespace mars::fem
