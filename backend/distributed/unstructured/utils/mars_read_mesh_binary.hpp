#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <filesystem>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <mpi.h>
#include "mars.hpp"

#ifdef __CUDACC__
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/copy.h>
#endif


namespace mars
{

template<int N, typename KeyType = unsigned>
auto createNVectors(size_t size)
{
    if constexpr (N == 3)
    {
        return std::tuple<std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>>(
            std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size));
    }
    else if constexpr (N == 4)
    {
        return std::tuple<std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>>(
            std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size));
    }
    else if constexpr (N == 8)
    {
        return std::tuple<std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>,
                          std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>, std::vector<KeyType>>(
            std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size),
            std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size), std::vector<KeyType>(size));
    }
}

// GPU-accelerated sort+unique for the node-dedup step. Only enabled when this
// header is compiled by nvcc (CUDA TU); plain C++ TUs (e.g. host-only unit
// tests) fall back to serial std::sort. The GPU path saves several seconds per
// rank at 100M+ keys.
template<typename KeyType>
inline std::vector<KeyType> dedupNodesGpu(std::vector<KeyType>& allNodes)
{
#ifdef __CUDACC__
    constexpr size_t cpuThreshold = 1ull << 20;  // ~1M keys: H2D round-trip not worth it
    if (allNodes.size() >= cpuThreshold)
    {
        thrust::device_vector<KeyType> d_keys(allNodes.size());
        thrust::copy(allNodes.begin(), allNodes.end(), d_keys.begin());
        thrust::sort(d_keys.begin(), d_keys.end());
        auto newEnd = thrust::unique(d_keys.begin(), d_keys.end());
        size_t uniqueCount = newEnd - d_keys.begin();

        std::vector<KeyType> result(uniqueCount);
        thrust::copy(d_keys.begin(), d_keys.begin() + uniqueCount, result.begin());
        return result;
    }
#endif
    std::sort(allNodes.begin(), allNodes.end());
    auto last = std::unique(allNodes.begin(), allNodes.end());
    allNodes.erase(last, allNodes.end());
    return std::move(allNodes);
}

// Header info the reader needs across ranks. Rank 0 detects from disk; all
// ranks receive via a single MPI_Bcast to avoid metadata-server pile-ups.
struct MeshFileHeader
{
    int      useInt64;        // bool packed as int for portability
    int      coordIsFloat32;
    uint64_t totalElements;
};

inline MeshFileHeader detectHeaderRank0(const std::string& meshDir)
{
    MeshFileHeader h{0, 0, 0};

    {
        std::ifstream f(meshDir + "/x.float32", std::ios::binary);
        if (f) { h.coordIsFloat32 = 1; }
        else
        {
            std::ifstream g(meshDir + "/x.double", std::ios::binary);
            if (!g) throw std::runtime_error("Failed to find coordinate files: neither x.float32 nor x.double in " + meshDir);
        }
    }

    std::string indexExt;
    {
        std::ifstream f(meshDir + "/i0.int64", std::ios::binary);
        if (f) { h.useInt64 = 1; indexExt = "int64"; }
        else
        {
            std::ifstream g(meshDir + "/i0.int32", std::ios::binary);
            if (!g) throw std::runtime_error("Failed to find index files: neither i0.int64 nor i0.int32 in " + meshDir);
            indexExt = "int32";
        }
    }

    size_t indexSize = h.useInt64 ? sizeof(int64_t) : sizeof(int32_t);
    std::ifstream f(meshDir + "/i0." + indexExt, std::ios::binary);
    f.seekg(0, std::ios::end);
    h.totalElements = static_cast<uint64_t>(f.tellg()) / indexSize;
    return h;
}

inline MeshFileHeader broadcastHeader(const std::string& meshDir, int rank)
{
    MeshFileHeader h{0, 0, 0};
    if (rank == 0) h = detectHeaderRank0(meshDir);
    MPI_Bcast(&h, sizeof(MeshFileHeader), MPI_BYTE, 0, MPI_COMM_WORLD);
    return h;
}

// Per-rank coord read. Each rank reads the contiguous [minNode, maxNode] range
// from the coord file (or chunks of it if too sparse), then scatters into
// neededNodes. This is the original path; a previous experiment with rank-0
// MPI_Scatterv added 12 minutes of pack-loop overhead at 256 ranks because
// the per-rank node ranges are dense slabs (not sparse), so each rank reading
// its own slab is already efficient.
template<typename KeyType, typename RealType>
inline void readCoordsPerRank(const std::string& path,
                              bool coordIsFloat32,
                              const std::vector<KeyType>& neededNodes,
                              std::vector<RealType>& out)
{
    const size_t nodeCount = neededNodes.size();
    out.resize(nodeCount);

    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("Failed to open coord file: " + path);

    size_t coordDiskSize = coordIsFloat32 ? sizeof(float) : sizeof(double);

    KeyType minNode = neededNodes.front();
    KeyType maxNode = neededNodes.back();
    size_t rangeCount = static_cast<size_t>(maxNode - minNode) + 1;

    double density = static_cast<double>(nodeCount) / static_cast<double>(rangeCount);
    // Buffer cap: 4 GB. For a 1024^3 mesh on 256 ranks, per-rank range is
    // ~17 MB so dense path triggers easily. Larger meshes may need adjustment.
    size_t maxBufBytes = 4ULL * 1024 * 1024 * 1024;
    size_t rangeBufBytes = rangeCount * coordDiskSize;

    if (density > 0.25 && rangeBufBytes <= maxBufBytes)
    {
        file.seekg(static_cast<std::streamoff>(static_cast<size_t>(minNode) * coordDiskSize));
        if (coordIsFloat32)
        {
            std::vector<float> buf(rangeCount);
            file.read(reinterpret_cast<char*>(buf.data()), rangeCount * sizeof(float));
            for (size_t j = 0; j < nodeCount; j++)
                out[j] = static_cast<RealType>(buf[static_cast<size_t>(neededNodes[j] - minNode)]);
        }
        else
        {
            std::vector<double> buf(rangeCount);
            file.read(reinterpret_cast<char*>(buf.data()), rangeCount * sizeof(double));
            for (size_t j = 0; j < nodeCount; j++)
                out[j] = static_cast<RealType>(buf[static_cast<size_t>(neededNodes[j] - minNode)]);
        }
    }
    else
    {
        // Sparse case: chunked reads sized to fit the buffer cap.
        size_t maxChunkNodes = maxBufBytes / coordDiskSize;
        size_t j = 0;
        while (j < nodeCount)
        {
            size_t chunkStart = j;
            KeyType chunkMinNode = neededNodes[j];
            size_t chunkEnd = j + 1;
            while (chunkEnd < nodeCount)
            {
                size_t chunkRange = static_cast<size_t>(neededNodes[chunkEnd] - chunkMinNode) + 1;
                if (chunkRange > maxChunkNodes) break;
                chunkEnd++;
            }
            size_t chunkRange = static_cast<size_t>(neededNodes[chunkEnd - 1] - chunkMinNode) + 1;

            file.seekg(static_cast<std::streamoff>(static_cast<size_t>(chunkMinNode) * coordDiskSize));

            if (coordIsFloat32)
            {
                std::vector<float> buf(chunkRange);
                file.read(reinterpret_cast<char*>(buf.data()), chunkRange * sizeof(float));
                for (size_t k = chunkStart; k < chunkEnd; k++)
                    out[k] = static_cast<RealType>(buf[static_cast<size_t>(neededNodes[k] - chunkMinNode)]);
            }
            else
            {
                std::vector<double> buf(chunkRange);
                file.read(reinterpret_cast<char*>(buf.data()), chunkRange * sizeof(double));
                for (size_t k = chunkStart; k < chunkEnd; k++)
                    out[k] = static_cast<RealType>(buf[static_cast<size_t>(neededNodes[k] - chunkMinNode)]);
            }
            j = chunkEnd;
        }
    }
}

/**
 * Read mesh using element-based partitioning to ensure all elements have their nodes.
 *
 * Optimizations for large-scale meshes (1B+ elements):
 *   - Rank 0 detects file format and broadcasts the header (avoids P metadata stats).
 *   - Each rank reads its connectivity slice independently (8 small sequential reads).
 *   - Node dedup runs on GPU via thrust::sort + thrust::unique (was serial std::sort).
 *   - Each rank reads its own contiguous coord-file slab (dense path) or chunks (sparse).
 */
template<int N, typename RealType = float, typename KeyType = unsigned>
inline auto readMeshWithElementPartitioning(const std::string& meshDir, int rank, int numRanks)
{
    // STEP 0: Detect file formats once on rank 0, broadcast to everyone.
    MeshFileHeader header = broadcastHeader(meshDir, rank);
    bool useInt64        = header.useInt64 != 0;
    bool coordIsFloat32  = header.coordIsFloat32 != 0;
    size_t indexSize     = useInt64 ? sizeof(int64_t) : sizeof(int32_t);
    std::string indexExt = useInt64 ? "int64" : "int32";
    std::string coordExt = coordIsFloat32 ? "float32" : "double";
    size_t total_elements = header.totalElements;

    if (rank == 0)
    {
        std::cout << "Binary mesh: " << total_elements << " elements, index=" << indexExt
                  << ", coords=" << coordExt << std::endl;
    }

    // Per-rank element slice
    size_t elemPerRank  = total_elements / static_cast<size_t>(numRanks);
    size_t elemStartIdx = static_cast<size_t>(rank) * elemPerRank;
    size_t elemEndIdx   = (rank == numRanks - 1) ? total_elements : elemStartIdx + elemPerRank;
    size_t elementCount = elemEndIdx - elemStartIdx;

    // STEP 1: Read connectivity (per-rank slice from each i*.intXX file).
    auto rawConnectivity = createNVectors<N, KeyType>(elementCount);

    int connFileIdx = 0;
    auto read_connectivity = [&](auto& vec)
    {
        std::ifstream idx_file(meshDir + "/i" + std::to_string(connFileIdx++) + "." + indexExt, std::ios::binary);
        if (!idx_file) { throw std::runtime_error("Failed to open index file i" + std::to_string(connFileIdx - 1) + "." + indexExt); }

        idx_file.seekg(static_cast<std::streamoff>(elemStartIdx * indexSize));

        if (useInt64)
        {
            std::vector<int64_t> tempVec(elementCount);
            idx_file.read(reinterpret_cast<char*>(tempVec.data()), elementCount * sizeof(int64_t));
            std::transform(tempVec.begin(), tempVec.end(), vec.begin(),
                           [](int64_t val) { return static_cast<KeyType>(val); });
        }
        else
        {
            std::vector<int32_t> tempVec(elementCount);
            idx_file.read(reinterpret_cast<char*>(tempVec.data()), elementCount * sizeof(int32_t));
            std::transform(tempVec.begin(), tempVec.end(), vec.begin(),
                           [](int32_t val) { return static_cast<KeyType>(val); });
        }
    };

    std::apply([&](auto&... vecs) { (read_connectivity(vecs), ...); }, rawConnectivity);

    // STEP 2: Identify unique nodes via GPU sort+unique. Serial std::sort on
    // 134M keys per rank cost ~3.5 s; GPU does it in ~50 ms.
    std::vector<KeyType> allNodes;
    allNodes.reserve(elementCount * N);
    auto collect_nodes = [&](const auto& vec)
    {
        allNodes.insert(allNodes.end(), vec.begin(), vec.end());
    };
    std::apply([&](const auto&... vecs) { (collect_nodes(vecs), ...); }, rawConnectivity);

    std::vector<KeyType> neededNodes = dedupNodesGpu<KeyType>(allNodes);
    size_t nodeCount = neededNodes.size();

    // STEP 3: Read coordinates per rank from contiguous file slabs.
    std::vector<RealType> x_data, y_data, z_data;
    readCoordsPerRank<KeyType, RealType>(meshDir + "/x." + coordExt, coordIsFloat32, neededNodes, x_data);
    readCoordsPerRank<KeyType, RealType>(meshDir + "/y." + coordExt, coordIsFloat32, neededNodes, y_data);
    readCoordsPerRank<KeyType, RealType>(meshDir + "/z." + coordExt, coordIsFloat32, neededNodes, z_data);

    // STEP 4: Adjust connectivity to local indices via binary search.
    auto localConnectivity = createNVectors<N, KeyType>(elementCount);

    auto adjust_indices = [&](const auto& global_vec, auto& local_vec)
    {
        for (size_t j = 0; j < elementCount; j++)
        {
            auto it = std::lower_bound(neededNodes.begin(), neededNodes.end(), global_vec[j]);
            local_vec[j] = static_cast<KeyType>(std::distance(neededNodes.begin(), it));
        }
    };

    std::apply(
        [&](const auto&... global_vecs)
        {
            std::apply([&](auto&... local_vecs) { (adjust_indices(global_vecs, local_vecs), ...); }, localConnectivity);
        },
        rawConnectivity);

    return std::make_tuple(nodeCount, elementCount, std::move(x_data), std::move(y_data), std::move(z_data),
                           std::move(localConnectivity), neededNodes);
}

/**
 * Adapter function that provides backwards compatibility with the old API
 * but uses the new element-based partitioning internally.
 */
template<typename RealType = float, typename KeyType = unsigned>
inline auto readMeshCoordinatesBinary(const std::string& meshDir, int rank, int numRanks)
{
    // Try to read all data with element-based partitioning
    auto [nodeCount, elementCount, x_data, y_data, z_data, _, localToGlobal] =
        readMeshWithElementPartitioning<4, RealType, KeyType>(meshDir, rank, numRanks);

    // Return in the old format for compatibility
    return std::make_tuple(nodeCount, 0, // nodeStartIdx is now always 0 in the new approach
                           std::move(x_data), std::move(y_data), std::move(z_data));
}

/**
 * Adapter function for connectivity data that provides backwards compatibility
 * with the old API but uses the new element-based partitioning internally.
 */
template<int N, typename RealType = float, typename KeyType = unsigned>
inline auto readMeshConnectivityBinaryTuple(const std::string& meshDir, int nodeStartIdx, int rank, int numRanks)
{
    // Try to read all data with element-based partitioning - reuse existing data
    auto [nodeCount, elementCount, x_data, y_data, z_data, connectivity, localToGlobal] =
        readMeshWithElementPartitioning<N, RealType, KeyType>(meshDir, rank, numRanks);

    // Return in the old format for compatibility - nodeStartIdx is ignored
    return std::make_tuple(elementCount, std::move(connectivity));
}

/**
 * Adapter for vector-based connectivity API
 */
template<typename RealType = float, typename KeyType = unsigned>
inline auto
readMeshConnectivityBinary(const std::string& meshDir, int nodesPerElement, size_t nodeStartIdx, int rank, int numRanks)
{
    // Call the appropriate tuple version based on nodesPerElement
    if (nodesPerElement == 4)
    {
        auto [elemCount, conn] = readMeshConnectivityBinaryTuple<4, RealType, KeyType>(meshDir, nodeStartIdx, rank, numRanks);
        std::vector<std::vector<KeyType>> result(4);

        // Convert from tuple to vector-of-vectors format
        result[0] = std::get<0>(conn);
        result[1] = std::get<1>(conn);
        result[2] = std::get<2>(conn);
        result[3] = std::get<3>(conn);

        return std::make_tuple(elemCount, result);
    }
    else if (nodesPerElement == 3)
    {
        // Same implementation for triangles...
        auto [elemCount, conn] = readMeshConnectivityBinaryTuple<3, RealType, KeyType>(meshDir, nodeStartIdx, rank, numRanks);
        std::vector<std::vector<KeyType>> result(3);

        result[0] = std::get<0>(conn);
        result[1] = std::get<1>(conn);
        result[2] = std::get<2>(conn);

        return std::make_tuple(elemCount, result);
    }
    else if (nodesPerElement == 8)
    {
        // Same implementation for hexahedra...
        auto [elemCount, conn] = readMeshConnectivityBinaryTuple<8, RealType, KeyType>(meshDir, nodeStartIdx, rank, numRanks);
        std::vector<std::vector<KeyType>> result(8);

        result[0] = std::get<0>(conn);
        result[1] = std::get<1>(conn);
        result[2] = std::get<2>(conn);
        result[3] = std::get<3>(conn);
        result[4] = std::get<4>(conn);
        result[5] = std::get<5>(conn);
        result[6] = std::get<6>(conn);
        result[7] = std::get<7>(conn);

        return std::make_tuple(elemCount, result);
    }
    else { throw std::runtime_error("Unsupported number of nodes per element: " + std::to_string(nodesPerElement)); }
}

} // namespace mars
