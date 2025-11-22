#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <cstdint>

namespace mars
{
namespace fem
{

// VTK writer for unstructured tetrahedral meshes
class VTKWriter
{
public:
    VTKWriter() = default;

    // Write unstructured mesh and solution to VTU file
    template<typename DomainType, typename VectorType>
    void writeVTU(const std::string& filename,
                  const DomainType& domain,
                  const VectorType& solution,
                  const std::string& solutionName = "solution") const
    {
        // Get mesh data
        auto nodes = domain.getNodeCoordinates();
        auto elements = domain.getElementConnectivity();

        size_t numNodes = nodes.size() / 3;  // x,y,z per node
        size_t numElements = elements.size() / 4;  // 4 nodes per tet

        // Open file
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        // VTK header
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "<UnstructuredGrid>\n";
        file << "<Piece NumberOfPoints=\"" << numNodes << "\" NumberOfCells=\"" << numElements << "\">\n";

        // Points
        file << "<Points>\n";
        file << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t i = 0; i < numNodes; ++i) {
            file << nodes[3*i] << " " << nodes[3*i+1] << " " << nodes[3*i+2] << "\n";
        }
        file << "</DataArray>\n";
        file << "</Points>\n";

        // Cells
        file << "<Cells>\n";
        file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (size_t i = 0; i < numElements; ++i) {
            file << elements[4*i] << " " << elements[4*i+1] << " " 
                 << elements[4*i+2] << " " << elements[4*i+3] << "\n";
        }
        file << "</DataArray>\n";

        file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (size_t i = 1; i <= numElements; ++i) {
            file << (4 * i) << "\n";
        }
        file << "</DataArray>\n";

        file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (size_t i = 0; i < numElements; ++i) {
            file << "10\n";  // VTK_TETRA = 10
        }
        file << "</DataArray>\n";
        file << "</Cells>\n";

        // Point data
        file << "<PointData Scalars=\"" << solutionName << "\">\n";
        file << "<DataArray type=\"Float32\" Name=\"" << solutionName << "\" format=\"ascii\">\n";
        for (size_t i = 0; i < numNodes; ++i) {
            file << solution[i] << "\n";
        }
        file << "</DataArray>\n";
        file << "</PointData>\n";

        // Close
        file << "</Piece>\n";
        file << "</UnstructuredGrid>\n";
        file << "</VTKFile>\n";

        file.close();
    }
};

} // namespace fem
} // namespace mars
