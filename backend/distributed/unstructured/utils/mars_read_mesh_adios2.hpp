#include <adios2.h>
#include <algorithm>  // Include this header for std::min
#include "mars.hpp"
#include "mars_globals.hpp"

void readMeshInParallel(const std::string& filename, int rank, int size) {
    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO io = adios.DeclareIO("MeshIO");

    // Open the file
    adios2::Engine engine = io.Open(filename, adios2::Mode::ReadRandomAccess);

    // Read metadata attributes
    int numPoints = -1;
    int numCells = -1;

    auto attrNumPoints = io.InquireAttribute<int>("NumPoints");
    if (attrNumPoints) {
        numPoints = attrNumPoints.Data()[0];
    }

    auto attrNumCells = io.InquireAttribute<int>("NumCells");
    if (attrNumCells) {
        numCells = attrNumCells.Data()[0];
    }

    // Print the metadata
    std::cout << "Mesh metadata: " << numPoints << " points, " << numCells << " cells\n";

    // Inquire variables - using exact names from writer
    auto varNodes = io.InquireVariable<double>("nodes");
    auto varTets = io.InquireVariable<int>("tets");
    auto varTypes = io.InquireVariable<uint8_t>("types");
    auto varOffsets = io.InquireVariable<int>("offsets");

    if (!varNodes || !varTets) {
        std::cerr << "Error: Required variables not found in file\n";
        engine.Close();
        return;
    }

    // Set selection for reading - exactly match the structure in the writer
    const std::size_t totalNodes = varNodes.Shape()[0];
    const std::size_t nodesPerProcess = totalNodes / size;
    const std::size_t remainderNodes = totalNodes % size;
    const std::size_t startNode = rank * nodesPerProcess + std::min<std::size_t>(rank, remainderNodes);
    const std::size_t countNode = nodesPerProcess + (rank < remainderNodes ? 1 : 0);

    const std::size_t totalTets = varTets.Shape()[0];
    const std::size_t tetsPerProcess = totalTets / size;
    const std::size_t remainderTets = totalTets % size;
    const std::size_t startTet = rank * tetsPerProcess + std::min<std::size_t>(rank, remainderTets);
    const std::size_t countTet = tetsPerProcess + (rank < remainderTets ? 1 : 0);

    varNodes.SetSelection({{startNode, 0}, {countNode, 3}});
    varTets.SetSelection({{startTet, 0}, {countTet, 4}});

    // Read data
    std::vector<double> nodes(countNode * 3);
    std::vector<int> tets(countTet * 4);
    std::vector<uint8_t> types;
    std::vector<int> offsets;

    engine.Get(varNodes, nodes.data());
    engine.Get(varTets, tets.data());

    // Also read types and offsets if available
    if (varTypes) {
        types.resize(numCells);
        engine.Get(varTypes, types.data());
    }

    if (varOffsets) {
        offsets.resize(numCells);
        engine.Get(varOffsets, offsets.data());
    }

    // Perform the reads
    engine.PerformGets();

    // Close the engine
    engine.Close();

    // Print the read data
    std::cout << "Process " << rank << " read " << countNode << " nodes and " << countTet << " tetrahedra.\n";

    // Print nodes
    std::cout << "Nodes:\n";
    for (std::size_t i = 0; i < std::min<std::size_t>(countNode, 8); ++i) {  // Limit output to first 8 nodes
        std::cout << i << ": (" << nodes[i * 3] << ", " << nodes[i * 3 + 1] << ", " << nodes[i * 3 + 2] << ")\n";
    }

    // Print tetrahedra
    std::cout << "Tetrahedra:\n";
    for (std::size_t i = 0; i < std::min<std::size_t>(countTet, 6); ++i) {  // Limit output to first 6 tets
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
