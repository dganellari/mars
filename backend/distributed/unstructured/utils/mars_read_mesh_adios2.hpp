#include <adios2.h>

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
    const std::size_t startNode = rank * nodesPerProcess + std::min(rank, remainderNodes);
    const std::size_t countNode = nodesPerProcess + (rank < remainderNodes ? 1 : 0);

    const std::size_t totalTets = varTets.Shape()[0];
    const std::size_t tetsPerProcess = totalTets / size;
    const std::size_t remainderTets = totalTets % size;
    const std::size_t startTet = rank * tetsPerProcess + std::min(rank, remainderTets);
    const std::size_t countTet = tetsPerProcess + (rank < remainderTets ? 1 : 0);

    varNodes.SetSelection({{startNode, 0}, {countNode, 3}});
    varTets.SetSelection({{startTet, 0}, {countTet, 4}});

    // Read data
    std::vector<std::array<double, 3>> nodes(countNode);
    std::vector<std::array<int, 4>> tets(countTet);

    engine.Get(varNodes, nodes.data());
    engine.Get(varTets, tets.data());

    // Perform the reads
    engine.PerformGets();

    // Close the engine
    engine.Close();

    // Print the data read by each process
    std::cout << "Process " << rank << " read nodes:\n";
    for (const auto& node : nodes) {
        std::cout << node[0] << " " << node[1] << " " << node[2] << "\n";
    }

    std::cout << "Process " << rank << " read tetrahedra:\n";
    for (const auto& tet : tets) {
        std::cout << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << "\n";
    }
}

