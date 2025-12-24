#pragma once

#include <fstream>
#include <vector>
#include <stdexcept>
#include <cstdint>

namespace mars {
namespace fem {

// Simple CRS graph reader for STK-compatible binary graph files
// Format (binary):
//   - uint64_t: number of nodes (owned)
//   - uint64_t: number of indices (NNZ)
//   - uint64_t[nNodes+1]: row pointers
//   - uint64_t[NNZ]: column indices
class GraphReader {
public:
    using IndexType = uint64_t;

    std::vector<IndexType> rowPtr;
    std::vector<IndexType> colInd;
    IndexType numNodes;
    IndexType numIndices;

    GraphReader() : numNodes(0), numIndices(0) {}

    // Read from STK-format binary graph file
    bool readBinary(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Cannot open graph file: " + filename);
        }

        // STK format has header: [version?, numNodes, numNodes_dup?, numIndices]
        // Skip version/header (first uint64)
        IndexType header;
        file.read(reinterpret_cast<char*>(&header), sizeof(IndexType));
        if (!file) {
            throw std::runtime_error("Failed to read header from graph file");
        }

        // Read number of nodes
        file.read(reinterpret_cast<char*>(&numNodes), sizeof(IndexType));
        if (!file) {
            throw std::runtime_error("Failed to read numNodes from graph file");
        }

        // Read what appears to be duplicate numNodes or global offset (skip it)
        IndexType skip_value;
        file.read(reinterpret_cast<char*>(&skip_value), sizeof(IndexType));
        if (!file) {
            throw std::runtime_error("Graph file format error reading skip value");
        }
        // Note: This value might be global_offset or other metadata, not validating it

        // Read number of indices
        file.read(reinterpret_cast<char*>(&numIndices), sizeof(IndexType));
        if (!file) {
            throw std::runtime_error("Failed to read numIndices from graph file");
        }

        // Read row pointers (numNodes + 1 entries)
        rowPtr.resize(numNodes + 1);
        file.read(reinterpret_cast<char*>(rowPtr.data()),
                  (numNodes + 1) * sizeof(IndexType));
        if (!file) {
            throw std::runtime_error("Failed to read rowPtr from graph file");
        }

        // Verify last entry matches numIndices
        if (rowPtr[numNodes] != numIndices) {
            throw std::runtime_error("Graph file inconsistent: rowPtr[numNodes]=" +
                                   std::to_string(rowPtr[numNodes]) +
                                   " != numIndices=" + std::to_string(numIndices));
        }

        // Read column indices
        colInd.resize(numIndices);
        file.read(reinterpret_cast<char*>(colInd.data()),
                  numIndices * sizeof(IndexType));
        if (!file) {
            throw std::runtime_error("Failed to read colInd from graph file");
        }

        file.close();
        return true;
    }

    // Get average NNZ per row
    double avgNnzPerRow() const {
        return numNodes > 0 ? static_cast<double>(numIndices) / numNodes : 0.0;
    }

    // Check if entry (row, col) exists in sparsity pattern
    bool hasEntry(IndexType row, IndexType col) const {
        if (row >= numNodes) return false;

        IndexType start = rowPtr[row];
        IndexType end = rowPtr[row + 1];

        // Binary search in sorted column indices
        for (IndexType i = start; i < end; ++i) {
            if (colInd[i] == col) return true;
            if (colInd[i] > col) return false;  // Assumes sorted
        }
        return false;
    }

    // Print statistics
    void printStats() const {
        std::cout << "Graph statistics:" << std::endl;
        std::cout << "  Nodes: " << numNodes << std::endl;
        std::cout << "  Indices (NNZ): " << numIndices << std::endl;
        std::cout << "  Average NNZ per row: " << avgNnzPerRow() << std::endl;
    }
};

} // namespace fem
} // namespace mars
