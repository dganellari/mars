#pragma once

#include <adios2.h>
#include <algorithm>
#include "mars.hpp"
#include "mars_globals.hpp"

namespace {
    // Helper function to open mesh file
    bool openMeshFile(const std::string& filename, adios2::ADIOS& adios, adios2::IO& io, adios2::Engine& engine) {
        try {
            io = adios.DeclareIO("MeshIO");
            engine = io.Open(filename, adios2::Mode::ReadRandomAccess);

            // Check if the engine is valid/open
            if (!engine) {
                return false;
            }
            return true;
        } catch (std::exception& e) {
            std::cerr << "Error opening ADIOS2 file: " << e.what() << std::endl;
            return false;
        }
    }

    // Helper function to read mesh metadata
    void readMeshMetadata(adios2::IO& io, int& numPoints, int& numCells, bool printInfo = false, int rank = 0) {
        auto attrNumPoints = io.InquireAttribute<int>("NumPoints");
        if (attrNumPoints) {
            numPoints = attrNumPoints.Data()[0];
        }

        auto attrNumCells = io.InquireAttribute<int>("NumCells");
        if (attrNumCells) {
            numCells = attrNumCells.Data()[0];
        }

        // Print metadata (conditionally)
        if (printInfo && rank == 0) {
            std::cout << "Mesh metadata: " << numPoints << " points, " << numCells << " cells\n";
        }
    }

    // Helper function to inquire mesh variables
    bool inquireMeshVariables(adios2::IO& io,
                              adios2::Variable<double>& varNodes,
                              adios2::Variable<int>& varTets,
                              adios2::Variable<uint8_t>& varTypes,
                              adios2::Variable<int>& varOffsets) {
        varNodes = io.InquireVariable<double>("nodes");
        varTets = io.InquireVariable<int>("tets");
        varTypes = io.InquireVariable<uint8_t>("types");
        varOffsets = io.InquireVariable<int>("offsets");

        return (varNodes && varTets);
    }

    // Helper function to calculate partitioning
    void calculatePartitioning(const adios2::Variable<double>& varNodes,
                               const adios2::Variable<int>& varTets,
                               int rank,
                               int size,
                               std::size_t& startNode,
                               std::size_t& countNode,
                               std::size_t& startTet,
                               std::size_t& countTet) {
        const std::size_t totalNodes = varNodes.Shape()[0];
        const std::size_t nodesPerProcess = totalNodes / size;
        const std::size_t remainderNodes = totalNodes % size;
        startNode = rank * nodesPerProcess + std::min<std::size_t>(rank, remainderNodes);
        countNode = nodesPerProcess + (rank < remainderNodes ? 1 : 0);

        const std::size_t totalTets = varTets.Shape()[0];
        const std::size_t tetsPerProcess = totalTets / size;
        const std::size_t remainderTets = totalTets % size;
        startTet = rank * tetsPerProcess + std::min<std::size_t>(rank, remainderTets);
        countTet = tetsPerProcess + (rank < remainderTets ? 1 : 0);
    }
}  // namespace

// Utility function to read mesh data - returns data
void readMeshData(const std::string& filename,
                  int rank,
                  int size,
                  std::vector<double>& nodes,
                  std::vector<int>& tets,
                  std::size_t& localNodeCount,
                  std::size_t& localTetCount) {
    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO io;
    adios2::Engine engine;

    if (!openMeshFile(filename, adios, io, engine)) {
        localNodeCount = 0;
        localTetCount = 0;
        return;
    }

    // Read metadata
    int numPoints = -1, numCells = -1;
    readMeshMetadata(io, numPoints, numCells, true, rank);

    // Inquire variables
    adios2::Variable<double> varNodes;
    adios2::Variable<int> varTets;
    adios2::Variable<uint8_t> varTypes;
    adios2::Variable<int> varOffsets;

    if (!inquireMeshVariables(io, varNodes, varTets, varTypes, varOffsets)) {
        std::cerr << "Error: Required variables not found in file\n";
        engine.Close();
        localNodeCount = 0;
        localTetCount = 0;
        return;
    }

    // Calculate partitioning
    std::size_t startNode, startTet;
    calculatePartitioning(varNodes, varTets, rank, size, startNode, localNodeCount, startTet, localTetCount);

    // Set selection and prepare vectors
    varNodes.SetSelection({{startNode, 0}, {localNodeCount, 3}});
    varTets.SetSelection({{startTet, 0}, {localTetCount, 4}});

    nodes.resize(localNodeCount * 3);
    tets.resize(localTetCount * 4);

    // Read data
    engine.Get(varNodes, nodes.data());
    engine.Get(varTets, tets.data());
    engine.PerformGets();
    engine.Close();
}

// Original function - prints data
void readMeshInParallel(const std::string& filename, int rank, int size) {
    // Use the data-returning version to read the mesh
    std::vector<double> nodes;
    std::vector<int> tets;
    std::size_t countNode, countTet;

    // Read the mesh data
    readMeshData(filename, rank, size, nodes, tets, countNode, countTet);

    // Early return if reading failed
    if (countNode == 0 || countTet == 0) {
        return;
    }

    // Read additional data (types and offsets) if needed
    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO io;
    adios2::Engine engine;

    if (openMeshFile(filename, adios, io, engine)) {
        std::vector<uint8_t> types;
        std::vector<int> offsets;

        auto varTypes = io.InquireVariable<uint8_t>("types");
        auto varOffsets = io.InquireVariable<int>("offsets");

        int numCells = -1;
        auto attrNumCells = io.InquireAttribute<int>("NumCells");
        if (attrNumCells) {
            numCells = attrNumCells.Data()[0];
        }

        if (varTypes) {
            types.resize(numCells);
            engine.Get(varTypes, types.data());
        }

        if (varOffsets) {
            offsets.resize(numCells);
            engine.Get(varOffsets, offsets.data());
        }

        engine.PerformGets();
        engine.Close();

        // Print the read data
        std::cout << "Process " << rank << " read " << countNode << " nodes and " << countTet << " tetrahedra.\n";

        // Print nodes
        std::cout << "Nodes:\n";
        for (std::size_t i = 0; i < std::min<std::size_t>(countNode, 8); ++i) {
            std::cout << i << ": (" << nodes[i * 3] << ", " << nodes[i * 3 + 1] << ", " << nodes[i * 3 + 2] << ")\n";
        }

        // Print tetrahedra
        std::cout << "Tetrahedra:\n";
        for (std::size_t i = 0; i < std::min<std::size_t>(countTet, 6); ++i) {
            std::cout << i << ": " << tets[i * 4] << " " << tets[i * 4 + 1] << " " << tets[i * 4 + 2] << " "
                      << tets[i * 4 + 3] << "\n";
        }

        // Print types and offsets if available
        if (!types.empty()) {
            std::cout << "Cell Types: ";
            for (std::size_t i = 0; i < std::min<std::size_t>(types.size(), 6); ++i) {
                std::cout << static_cast<int>(types[i]) << " ";
            }
            std::cout << "\n";
        }

        if (!offsets.empty()) {
            std::cout << "Offsets: ";
            for (std::size_t i = 0; i < std::min<std::size_t>(offsets.size(), 6); ++i) {
                std::cout << offsets[i] << " ";
            }
            std::cout << "\n";
        }
    }
}
