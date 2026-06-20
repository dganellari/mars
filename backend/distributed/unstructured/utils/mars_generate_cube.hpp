#pragma once

// Procedural per-rank structured-cube generation, drop-in for
// readMeshWithElementPartitioning (mars_read_mesh_binary.hpp): each rank computes
// ONLY its element slice of an Ncells^3 hex cube on [0,1]^3 -- no mesh file. This
// is the prerequisite for trillion-DOF runs, where an O(10 TB) mesh file is
// infeasible. The output format is identical to the file reader:
//   - per-rank coords for the nodes this rank's elements reference (local-indexed)
//   - 8 connectivity arrays remapped to those local node indices
// cstone's sync() then merges nodes shared across ranks by SFC key from coords,
// exactly as for a file mesh. Element->rank assignment is a flat contiguous split
// (rank r owns elements [r*per, (r+1)*per)); cstone redistributes by SFC anyway.
//
// Node dedup is host sort+unique here (simple, correct). For very large per-rank
// slices (hundreds of millions of elements) swap in the GPU dedup the reader uses
// (dedupNodesGpu) -- noted, not yet wired.

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

namespace mars {

// Returns (nodeCount, elementCount, x, y, z, localConn[8]) for this rank's slice
// of an Ncells^3 hex cube. coords are local-node-indexed; localConn holds local
// node indices in hex corner order matching generate_hex_cube.py:
//   0:(i,j,k) 1:(i+1,j,k) 2:(i+1,j+1,k) 3:(i,j+1,k)
//   4:(i,j,k+1) 5:(i+1,j,k+1) 6:(i+1,j+1,k+1) 7:(i,j+1,k+1)
template<typename RealType, typename KeyType>
inline auto generateCubeElementPartition(size_t Ncells, int rank, int numRanks)
{
    const size_t Np1 = Ncells + 1;
    const size_t totalElems = Ncells * Ncells * Ncells;
    const size_t per   = totalElems / static_cast<size_t>(numRanks);
    const size_t start = static_cast<size_t>(rank) * per;
    const size_t end   = (rank == numRanks - 1) ? totalElems : start + per;
    const size_t elementCount = end - start;

    auto gid = [Np1](size_t cx, size_t cy, size_t cz) -> KeyType {
        return static_cast<KeyType>((cx * Np1 + cy) * Np1 + cz);
    };

    std::array<std::vector<KeyType>, 8> gconn;
    for (auto& v : gconn) v.resize(elementCount);
    std::vector<KeyType> allNodes;
    allNodes.reserve(elementCount * 8);

    for (size_t e = start; e < end; ++e)
    {
        const size_t ex = e / (Ncells * Ncells);
        const size_t r  = e % (Ncells * Ncells);
        const size_t ey = r / Ncells;
        const size_t ez = r % Ncells;
        const size_t li = e - start;
        const KeyType c[8] = {
            gid(ex, ey, ez),     gid(ex + 1, ey, ez),     gid(ex + 1, ey + 1, ez),     gid(ex, ey + 1, ez),
            gid(ex, ey, ez + 1), gid(ex + 1, ey, ez + 1), gid(ex + 1, ey + 1, ez + 1), gid(ex, ey + 1, ez + 1)};
        for (int k = 0; k < 8; ++k) { gconn[k][li] = c[k]; allNodes.push_back(c[k]); }
    }

    // Unique global node ids this rank references.
    std::sort(allNodes.begin(), allNodes.end());
    allNodes.erase(std::unique(allNodes.begin(), allNodes.end()), allNodes.end());
    const size_t nodeCount = allNodes.size();

    // Coords of the needed nodes (decode global id -> (cx,cy,cz) -> physical).
    const RealType h = RealType(1) / static_cast<RealType>(Ncells);
    std::vector<RealType> x(nodeCount), y(nodeCount), z(nodeCount);
    for (size_t i = 0; i < nodeCount; ++i)
    {
        const KeyType g = allNodes[i];
        const size_t cz = static_cast<size_t>(g) % Np1;
        const size_t cy = (static_cast<size_t>(g) / Np1) % Np1;
        const size_t cx = static_cast<size_t>(g) / (Np1 * Np1);
        x[i] = static_cast<RealType>(cx) * h;
        y[i] = static_cast<RealType>(cy) * h;
        z[i] = static_cast<RealType>(cz) * h;
    }

    // Remap connectivity to local node indices (binary search into allNodes).
    std::array<std::vector<KeyType>, 8> lconn;
    for (int k = 0; k < 8; ++k)
    {
        lconn[k].resize(elementCount);
        for (size_t j = 0; j < elementCount; ++j)
        {
            auto it = std::lower_bound(allNodes.begin(), allNodes.end(), gconn[k][j]);
            lconn[k][j] = static_cast<KeyType>(std::distance(allNodes.begin(), it));
        }
    }

    return std::make_tuple(nodeCount, elementCount, std::move(x), std::move(y), std::move(z), std::move(lconn));
}

} // namespace mars
