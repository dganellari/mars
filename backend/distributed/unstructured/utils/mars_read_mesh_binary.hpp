#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <filesystem>
#include <algorithm>
#include "mars.hpp"


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

/**
 * Read mesh using element-based partitioning to ensure all elements have their nodes
 * @tparam N Number of nodes per element
 * @tparam RealType Type for coordinate data (float or double)
 * @param meshDir Directory containing mesh files
 * @param rank Current MPI rank
 * @param numRanks Total MPI ranks
 * @return Tuple containing: node count, element count, coordinates (x,y,z), and element connectivity
 */
template<int N, typename RealType = float, typename KeyType = unsigned>
inline auto readMeshWithElementPartitioning(const std::string& meshDir, int rank, int numRanks)
{
    // Detect coordinate file format on disk (float32 or double)
    // The file may not match RealType — e.g., RealType=double but files are x.float32
    bool coordIsFloat32 = false;
    std::string coordExt;
    {
        std::ifstream test_coord(meshDir + "/x.float32", std::ios::binary);
        if (test_coord)
        {
            coordExt = "float32";
            coordIsFloat32 = true;
        }
        else
        {
            std::ifstream test_coord2(meshDir + "/x.double", std::ios::binary);
            if (test_coord2)
            {
                coordExt = "double";
                coordIsFloat32 = false;
            }
            else
            {
                throw std::runtime_error("Failed to find coordinate files: neither x.float32 nor x.double in " + meshDir);
            }
        }
    }
    size_t coordDiskSize = coordIsFloat32 ? sizeof(float) : sizeof(double);

    // STEP 1: Detect connectivity file format (int64 vs int32) and read element count
    bool useInt64 = false;
    size_t indexSize = sizeof(int32_t);
    std::string indexExt = "int32";

    // Try int64 first, fall back to int32
    {
        std::ifstream test_file(meshDir + "/i0.int64", std::ios::binary);
        if (test_file)
        {
            useInt64 = true;
            indexSize = sizeof(int64_t);
            indexExt = "int64";
        }
        else
        {
            test_file.open(meshDir + "/i0.int32", std::ios::binary);
            if (!test_file) { throw std::runtime_error("Failed to open index file: neither i0.int64 nor i0.int32 found in " + meshDir); }
        }
    }

    // Determine total element count from first index file
    size_t total_elements;
    {
        std::ifstream test_file(meshDir + "/i0." + indexExt, std::ios::binary);
        test_file.seekg(0, std::ios::end);
        total_elements = static_cast<size_t>(test_file.tellg()) / indexSize;
        test_file.close();
    }

    if (rank == 0)
    {
        std::cout << "Binary mesh: " << total_elements << " elements, index=" << indexExt << ", coords=" << coordExt << std::endl;
    }

    // Calculate this rank's element portion
    size_t elemPerRank  = total_elements / numRanks;
    size_t elemStartIdx = static_cast<size_t>(rank) * elemPerRank;
    size_t elemEndIdx   = (rank == numRanks - 1) ? total_elements : elemStartIdx + elemPerRank;
    size_t elementCount = elemEndIdx - elemStartIdx;

    // Read element connectivity
    auto rawConnectivity = createNVectors<N, KeyType>(elementCount);

    int i = 0;
    auto read_connectivity = [&](auto& vec)
    {
        std::ifstream idx_file(meshDir + "/i" + std::to_string(i++) + "." + indexExt, std::ios::binary);
        if (!idx_file) { throw std::runtime_error("Failed to open index file i" + std::to_string(i - 1) + "." + indexExt); }

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

    // STEP 2: Identify unique nodes needed by this rank using sorted vector (no hash tables)
    std::vector<KeyType> allNodes;
    allNodes.reserve(elementCount * N);

    auto collect_nodes = [&](const auto& vec)
    {
        allNodes.insert(allNodes.end(), vec.begin(), vec.end());
    };

    std::apply([&](const auto&... vecs) { (collect_nodes(vecs), ...); }, rawConnectivity);

    std::sort(allNodes.begin(), allNodes.end());
    auto last = std::unique(allNodes.begin(), allNodes.end());
    allNodes.erase(last, allNodes.end());

    std::vector<KeyType> neededNodes = std::move(allNodes);
    size_t nodeCount = neededNodes.size();

    // STEP 3: Read coordinates using contiguous range read
    std::vector<RealType> x_data(nodeCount), y_data(nodeCount), z_data(nodeCount);

    std::ifstream x_file(meshDir + "/x." + coordExt, std::ios::binary);
    std::ifstream y_file(meshDir + "/y." + coordExt, std::ios::binary);
    std::ifstream z_file(meshDir + "/z." + coordExt, std::ios::binary);

    if (!x_file || !y_file || !z_file) { throw std::runtime_error("Failed to open coordinate files in " + meshDir); }

    // Read coordinates: use contiguous range if dense, chunked reads if sparse
    KeyType minNode = neededNodes.front();
    KeyType maxNode = neededNodes.back();
    size_t rangeCount = static_cast<size_t>(maxNode - minNode) + 1;

    // Density ratio: if needed nodes are >25% of the range, read contiguous block
    // Otherwise, read in chunks to avoid huge temporary buffers
    double density = static_cast<double>(nodeCount) / static_cast<double>(rangeCount);
    // Also cap the buffer at 512 MB to avoid OOM
    size_t maxBufBytes = 512ULL * 1024 * 1024;
    size_t rangeBufBytes = rangeCount * coordDiskSize;

    if (density > 0.25 && rangeBufBytes <= maxBufBytes)
    {
        // Dense case: read contiguous range and scatter
        auto readAndScatter = [&](std::ifstream& file, std::vector<RealType>& out)
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
        };

        readAndScatter(x_file, x_data);
        readAndScatter(y_file, y_data);
        readAndScatter(z_file, z_data);
    }
    else
    {
        // Sparse case: read in chunks along the sorted neededNodes list
        // Each chunk reads a contiguous sub-range of the coordinate file
        size_t maxChunkNodes = maxBufBytes / coordDiskSize;

        auto readChunked = [&](std::ifstream& file, std::vector<RealType>& out)
        {
            size_t j = 0;
            while (j < nodeCount)
            {
                // Find chunk: start from neededNodes[j], include as many consecutive-ish nodes
                // as fit in the buffer
                size_t chunkStart = j;
                KeyType chunkMinNode = neededNodes[j];
                // Extend chunk until the range exceeds buffer or we run out of nodes
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
        };

        readChunked(x_file, x_data);
        readChunked(y_file, y_data);
        readChunked(z_file, z_data);
    }

    // STEP 4: Adjust connectivity to use local indices (binary search instead of hash map)
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
