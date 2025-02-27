#include <gtest/gtest.h>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <adios2.h>

void generateTestMesh() {
    // 8 vertices of unit cube
    std::vector<std::array<double,3>> nodes = {
        {0.0, 0.0, 0.0}, // 0
        {1.0, 0.0, 0.0}, // 1
        {1.0, 1.0, 0.0}, // 2
        {0.0, 1.0, 0.0}, // 3
        {0.0, 0.0, 1.0}, // 4
        {1.0, 0.0, 1.0}, // 5
        {1.0, 1.0, 1.0}, // 6
        {0.0, 1.0, 1.0}  // 7
    };

    // 5 tetrahedra sharing vertex 4
    std::vector<std::array<int,4>> tets = {
        {4, 0, 1, 2},
        {4, 2, 3, 0},
        {4, 5, 6, 2},
        {4, 6, 7, 2},
        {4, 2, 7, 3}
    };

    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("MeshIO");

    // Define variables
    auto varNodes = io.DefineVariable<double>("nodes", {nodes.size(), 3}, {0, 0}, {nodes.size(), 3});
    auto varTets = io.DefineVariable<int>("tets", {tets.size(), 4}, {0, 0}, {tets.size(), 4});

    // Create engine
    adios2::Engine engine = io.Open("test_mesh.bp", adios2::Mode::Write);

    // Write data
    engine.Put(varNodes, nodes.data());
    engine.Put(varTets, tets.data());

    // Close the engine
    engine.Close();
}

TEST(GenerateMeshTest, GeneratesMeshFile) {
    generateTestMesh();

    std::ifstream infile("test_mesh.bp");
    ASSERT_TRUE(infile.good());
    infile.close();
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
