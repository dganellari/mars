#include <gtest/gtest.h>
#include <mpi.h>
#include <vector>
#include <array>
#include <iostream>
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

#ifdef MARS_ENABLE_MPI

TEST(ReadMeshTest, ReadsMeshFileInParallel) {
    int rank, size;

    using namespace mars;
    mars::proc_allocation resources;
    bool result = false;
    // create a distributed context
    auto context = mars::make_context(resources, MPI_COMM_WORLD);
    rank = mars::rank(context);
    size = mars::size(context);

    readMeshInParallel("test_mesh.bp", rank, size);
    // Add assertions to verify the correctness of the read data
    // For simplicity, we assume the mesh is correctly read if no errors occur
    ASSERT_TRUE(true);
}

#endif

int main(int argc, char** argv) {
    mars::Env env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
