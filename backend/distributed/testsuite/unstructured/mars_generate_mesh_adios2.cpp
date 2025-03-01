#include <gtest/gtest.h>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <adios2.h>
#include <string>

void generateTestVTKMesh() {
    try {
        // Define the cube's 8 vertices (x, y, z coordinates)
        std::vector<std::array<double, 3>> points = {
            {0.0, 0.0, 0.0},  // V0: (0, 0, 0) - bottom front left
            {1.0, 0.0, 0.0},  // V1: (1, 0, 0) - bottom front right
            {1.0, 1.0, 0.0},  // V2: (1, 1, 0) - bottom back right
            {0.0, 1.0, 0.0},  // V3: (0, 1, 0) - bottom back left
            {0.0, 0.0, 1.0},  // V4: (0, 0, 1) - top front left
            {1.0, 0.0, 1.0},  // V5: (1, 0, 1) - top front right
            {1.0, 1.0, 1.0},  // V6: (1, 1, 1) - top back right
            {0.0, 1.0, 1.0}   // V7: (0, 1, 1) - top back left
        };

        // Define six tetrahedra to completely fill the cube
        std::vector<std::array<int, 4>> tetrahedra = {
            {0, 1, 3, 7},  // T1: bottom front left, bottom front right, bottom back left, top back left
            {1, 2, 3, 6},  // T2: bottom front right, bottom back right, bottom back left, top back right
            {1, 5, 6, 7},  // T3: bottom front right, top front right, top back right, top back left
            {0, 1, 4, 7},  // T4: bottom front left, bottom front right, top front left, top back left
            {1, 4, 5, 7},  // T5: bottom front right, top front left, top front right, top back left
            {1, 3, 6, 7}   // T6: bottom front right, bottom back left, top back right, top back left
        };

        // Write standard VTK format which ParaView will definitely read correctly
        std::ofstream vtkFile("tetrahedral_mesh.vtu");
        
        // VTK XML UnstructuredGrid header
        vtkFile << "<?xml version=\"1.0\"?>\n";
        vtkFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        vtkFile << "  <UnstructuredGrid>\n";
        vtkFile << "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << tetrahedra.size() << "\">\n";
        
        // Write points
        vtkFile << "      <Points>\n";
        vtkFile << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        vtkFile << "          ";
        for (const auto& point : points) {
            vtkFile << point[0] << " " << point[1] << " " << point[2] << " ";
        }
        vtkFile << "\n";
        vtkFile << "        </DataArray>\n";
        vtkFile << "      </Points>\n";
        
        // Write cells
        vtkFile << "      <Cells>\n";
        vtkFile << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        vtkFile << "          ";
        for (const auto& tet : tetrahedra) {
            vtkFile << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << " ";
        }
        vtkFile << "\n";
        vtkFile << "        </DataArray>\n";
        
        // Write offsets
        vtkFile << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        vtkFile << "          ";
        for (int i = 1; i <= tetrahedra.size(); i++) {
            vtkFile << i * 4 << " ";
        }
        vtkFile << "\n";
        vtkFile << "        </DataArray>\n";
        
        // Write cell types (10 = VTK_TETRA)
        vtkFile << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        vtkFile << "          ";
        for (size_t i = 0; i < tetrahedra.size(); i++) {
            vtkFile << "10 ";
        }
        vtkFile << "\n";
        vtkFile << "        </DataArray>\n";
        vtkFile << "      </Cells>\n";
        
        // Add cell IDs for better visualization - MOVED OUTSIDE CELLS TAG
        vtkFile << "      <CellData>\n";
        vtkFile << "        <DataArray type=\"Int32\" Name=\"CellID\" format=\"ascii\">\n";
        vtkFile << "          ";
        for (size_t i = 0; i < tetrahedra.size(); i++) {
            vtkFile << i << " ";
        }
        vtkFile << "\n";
        vtkFile << "        </DataArray>\n";
        vtkFile << "      </CellData>\n";
        
        // End the file
        vtkFile << "    </Piece>\n";
        vtkFile << "  </UnstructuredGrid>\n";
        vtkFile << "</VTKFile>\n";
        
        vtkFile.close();
        std::cout << "Tetrahedral mesh written to 'tetrahedral_mesh.vtu' successfully.\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error generating mesh: " << e.what() << std::endl;
    }
}

TEST(GenerateVTKMesh, GeneratesVTKMeshFile) {
    generateTestVTKMesh();

    std::ifstream infile("tetrahedral_mesh.vtu");
    ASSERT_TRUE(infile.good());
    infile.close();
}

void generateTestMesh() {
    try {
        // Initialize ADIOS2
        adios2::ADIOS adios;

        // Define the cube's 8 vertices (x, y, z coordinates)
        std::vector<double> vertices = {
            0.0, 0.0, 0.0,  // V0
            1.0, 0.0, 0.0,  // V1
            1.0, 1.0, 0.0,  // V2
            0.0, 1.0, 0.0,  // V3
            0.0, 0.0, 1.0,  // V4
            1.0, 0.0, 1.0,  // V5
            1.0, 1.0, 1.0,  // V6
            0.0, 1.0, 1.0   // V7
        };

        // Define tetrahedra connectivity (which vertices make up each tetrahedron)
        std::vector<int> connectivity = {
            0, 1, 3, 7,  // T1
            1, 2, 3, 6,  // T2
            1, 5, 6, 7,  // T3
            0, 1, 4, 7,  // T4
            1, 4, 5, 7,  // T5
            1, 3, 6, 7   // T6
        };

        // Define cell types (10 = VTK_TETRA)
        std::vector<uint8_t> types(6, 10);

        // Define offsets (where each cell ends in connectivity array)
        std::vector<int> offsets = {4, 8, 12, 16, 20, 24};

        // Declare an IO object - make it clear this is for UnstructuredGrid
        adios2::IO io = adios.DeclareIO("UnstructuredGrid");

        // Open a file for writing
        adios2::Engine writer = io.Open("tetrahedral_mesh.bp", adios2::Mode::Write);

        // Store mesh type to indicate this is not image data
        io.DefineAttribute<std::string>("MeshType", "UnstructuredGrid");
        io.DefineAttribute<int>("NumPoints", 8);
        io.DefineAttribute<int>("NumCells", 6);
        io.DefineAttribute<std::string>("Description", "Tetrahedral mesh of a cube");

        // Define ADIOS2 variables with 2D shape for nodes/connectivity, 1D for types/offsets
        adios2::Variable<double> varVertices =
            io.DefineVariable<double>("nodes", {8, 3}, {0, 0}, {8, 3});

        adios2::Variable<int> varConnectivity =
            io.DefineVariable<int>("tets", {6, 4}, {0, 0}, {6, 4});

        adios2::Variable<uint8_t> varTypes =
            io.DefineVariable<uint8_t>("types", {6}, {0}, {6});

        adios2::Variable<int> varOffsets =
            io.DefineVariable<int>("offsets", {6}, {0}, {6});

        // Write the data
        writer.Put(varVertices, vertices.data());
        writer.Put(varConnectivity, connectivity.data());
        writer.Put(varTypes, types.data());
        writer.Put(varOffsets, offsets.data());

        // Add VTK schema with more explicit type information
        std::string vtkSchema = R"(
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="8" NumberOfCells="6">
      <Points>
        <DataArray Name="nodes" NumberOfComponents="3" type="Float64"/>
      </Points>
      <Cells>
        <DataArray Name="tets" type="Int32"/>
        <DataArray Name="types" type="UInt8"/>
        <DataArray Name="offsets" type="Int32"/>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)";
        io.DefineAttribute<std::string>("vtk.xml", vtkSchema);

        // Close the writer
        writer.Close();

        std::cout << "Tetrahedral mesh written to 'tetrahedral_mesh.bp' successfully.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error generating mesh: " << e.what() << std::endl;
    }
}

TEST(GenerateMesh, GeneratesMeshFile) {
    generateTestMesh();

    std::ifstream infile("tetrahedral_mesh.bp");
    ASSERT_TRUE(infile.good());
    infile.close();
}

void generateComplexMesh() {
    try {
        // Initialize ADIOS2
        adios2::ADIOS adios;

        // Define vertices for a complex shape (prism with additional features)
        // Base layer (z=0)
        std::vector<double> vertices = {
            0.0, 0.0, 0.0,   // V0: base layer, corner
            1.0, 0.0, 0.0,   // V1: base layer, corner
            1.0, 1.0, 0.0,   // V2: base layer, corner
            0.0, 1.0, 0.0,   // V3: base layer, corner
            // Middle layer (z=0.5)
            0.0, 0.0, 0.5,   // V4: middle layer, corner
            1.0, 0.0, 0.5,   // V5: middle layer, corner
            1.0, 1.0, 0.5,   // V6: middle layer, corner
            0.0, 1.0, 0.5,   // V7: middle layer, corner
            0.5, 0.5, 0.5,   // V8: middle layer, center
            // Top layer (z=1.0)
            0.0, 0.0, 1.0,   // V9: top layer, corner
            1.0, 0.0, 1.0,   // V10: top layer, corner
            1.0, 1.0, 1.0,   // V11: top layer, corner
            0.0, 1.0, 1.0,   // V12: top layer, corner
            0.5, 0.5, 1.5    // V13: apex point
        };

        // Define tetrahedra connectivity (20 tetrahedra)
        std::vector<int> connectivity = {
            // Lower section (between z=0 and z=0.5)
            0, 1, 3, 8,      // T1: bottom front left section
            1, 2, 3, 8,      // T2: bottom front right section
            0, 1, 4, 8,      // T3: connecting base to middle
            1, 2, 5, 8,      // T4: connecting base to middle
            2, 3, 6, 8,      // T5: connecting base to middle
            3, 0, 7, 8,      // T6: connecting base to middle
            4, 5, 7, 8,      // T7: middle section
            5, 6, 7, 8,      // T8: middle section
            // Upper section (between z=0.5 and z=1.0)
            4, 5, 8, 9,      // T9: connecting middle to top
            5, 6, 8, 10,     // T10: connecting middle to top
            6, 7, 8, 11,     // T11: connecting middle to top
            7, 4, 8, 12,     // T12: connecting middle to top
            9, 10, 12, 8,    // T13: upper section
            10, 11, 12, 8,   // T14: upper section
            // Pyramid on top
            9, 10, 12, 13,   // T15: pyramid section
            10, 11, 12, 13,  // T16: pyramid section
            9, 10, 13, 8,    // T17: pyramid section
            10, 11, 13, 8,   // T18: pyramid section
            11, 12, 13, 8,   // T19: pyramid section
            12, 9, 13, 8     // T20: pyramid section
        };

        // Define cell types (10 = VTK_TETRA)
        std::vector<uint8_t> types(20, 10);

        // Define offsets
        std::vector<int> offsets;
        for (int i = 1; i <= 20; i++) {
            offsets.push_back(i * 4);
        }

        // Declare an IO object
        adios2::IO io = adios.DeclareIO("UnstructuredGrid");

        // Open a file for writing
        adios2::Engine writer = io.Open("complex_mesh.bp", adios2::Mode::Write);

        // Store mesh type to indicate this is not image data
        io.DefineAttribute<std::string>("MeshType", "UnstructuredGrid");
        io.DefineAttribute<int>("NumPoints", static_cast<int>(vertices.size()/3));
        io.DefineAttribute<int>("NumCells", static_cast<int>(connectivity.size()/4));
        io.DefineAttribute<std::string>("Description", "Complex unstructured mesh with ~20 tetrahedral elements");

        // Define ADIOS2 variables
        size_t numVertices = vertices.size() / 3;
        size_t numCells = connectivity.size() / 4;
        
        adios2::Variable<double> varVertices =
            io.DefineVariable<double>("nodes", {numVertices, 3}, {0, 0}, {numVertices, 3});

        adios2::Variable<int> varConnectivity =
            io.DefineVariable<int>("tets", {numCells, 4}, {0, 0}, {numCells, 4});

        adios2::Variable<uint8_t> varTypes =
            io.DefineVariable<uint8_t>("types", {numCells}, {0}, {numCells});

        adios2::Variable<int> varOffsets =
            io.DefineVariable<int>("offsets", {numCells}, {0}, {numCells});

        // Write the data
        writer.Put(varVertices, vertices.data());
        writer.Put(varConnectivity, connectivity.data());
        writer.Put(varTypes, types.data());
        writer.Put(varOffsets, offsets.data());

        // Add VTK schema with more explicit type information
        std::string vtkSchema = R"(
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints=")" + std::to_string(numVertices) + R"(" NumberOfCells=")" + std::to_string(numCells) + R"(">
      <Points>
        <DataArray Name="nodes" NumberOfComponents="3" type="Float64"/>
      </Points>
      <Cells>
        <DataArray Name="tets" type="Int32"/>
        <DataArray Name="types" type="UInt8"/>
        <DataArray Name="offsets" type="Int32"/>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)";
        io.DefineAttribute<std::string>("vtk.xml", vtkSchema);

        // Close the writer
        writer.Close();

        std::cout << "Complex mesh written to 'complex_mesh.bp' successfully.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error generating complex mesh: " << e.what() << std::endl;
    }
}

TEST(GenerateMesh, GeneratesComplexMeshFile) {
    generateComplexMesh();

    std::ifstream infile("complex_mesh.bp");
    ASSERT_TRUE(infile.good());
    infile.close();
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
