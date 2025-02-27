#include <gtest/gtest.h>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <adios2.h>

void generateTestMesh() {
    try {
        // Initialize ADIOS2
        adios2::ADIOS adios;

        // Define the cube's 8 vertices (x, y, z coordinates)
        std::vector<double> vertices = {
            0.0, 0.0, 0.0,  // V0: (0, 0, 0)
            1.0, 0.0, 0.0,  // V1: (1, 0, 0)
            1.0, 1.0, 0.0,  // V2: (1, 1, 0)
            0.0, 1.0, 0.0,  // V3: (0, 1, 0)
            0.0, 0.0, 1.0,  // V4: (0, 0, 1)
            1.0, 0.0, 1.0,  // V5: (1, 0, 1)
            1.0, 1.0, 1.0,  // V6: (1, 1, 1)
            0.0, 1.0, 1.0   // V7: (0, 1, 1)
        };

        // Define the five tetrahedra (each with 4 vertex indices)
        std::vector<int> tetrahedra = {
            4, 0, 1, 2,  // T1
            4, 2, 3, 0,  // T2
            4, 5, 6, 2,  // T3
            4, 6, 7, 2,  // T4
            4, 2, 7, 3   // T5
        };

        // Declare an ADIOS2 IO object
        adios2::IO io = adios.DeclareIO("MeshIO");

        // Open a file for writing in ADIOS2 format
        adios2::Engine writer = io.Open("tetrahedral_mesh.bp", adios2::Mode::Write);

        // Define ADIOS2 variables for vertices and tetrahedra
        adios2::Variable<double> varVertices =
            io.DefineVariable<double>("vertices", {8, 3}, {0, 0}, {8, 3}  // Shape: 8 vertices, 3 coordinates each
            );
        adios2::Variable<int> varTetrahedra =
            io.DefineVariable<int>("tetrahedra", {5, 4}, {0, 0}, {5, 4}  // Shape: 5 tetrahedra, 4 indices each
            );

        // Write the data to the file
        writer.Put(varVertices, vertices.data());
        writer.Put(varTetrahedra, tetrahedra.data());

        // Close the writer
        writer.Close();

        std::cout << "Tetrahedral mesh written to 'tetrahedral_mesh.bp' successfully.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
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
