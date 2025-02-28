#include <adios2.h>
#include "mars.hpp"
#include "mars_globals.hpp"
#include <algorithm> // Include this header for std::min

void readMeshInParallel(const std::string& filename, int rank, int size) {
    adios2::ADIOS adios(MPI_COMM_WORLD);
    adios2::IO io = adios.DeclareIO("MeshIO");

    // Open the file
    adios2::Engine engine = io.Open(filename, adios2::Mode::Read);

    // Inquire variables
    auto varNodes = io.InquireVariable<double>("nodes");
    auto varTets = io.InquireVariable<int>("tets");

    // Set selection for reading
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

    engine.Get(varNodes, nodes.data());
    engine.Get(varTets, tets.data());

    // Perform the reads
    engine.PerformGets();

    // Close the engine
    engine.Close();

    // Print the read data for each rank
    std::cout << "Rank " << rank << " read " << countNode << " nodes and " << countTet << " tets.\n";
    std::cout << "Nodes:\n";
    for (std::size_t i = 0; i < countNode; ++i) {
        std::cout << nodes[i * 3] << " " << nodes[i * 3 + 1] << " " << nodes[i * 3 + 2] << "\n";
    }
    std::cout << "Tets:\n";
    for (std::size_t i = 0; i < countTet; ++i) {
        std::cout << tets[i * 4] << " " << tets[i * 4 + 1] << " " << tets[i * 4 + 2] << " " << tets[i * 4 + 3] << "\n";
    }
}

