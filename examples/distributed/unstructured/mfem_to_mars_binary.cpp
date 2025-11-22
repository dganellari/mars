/**
 * Convert MFEM mesh to MARS binary format with optional refinement
 * 
 * Usage: mfem_to_mars_binary input.mesh output_dir [refine_levels]
 * 
 * This tool:
 * 1. Reads MFEM mesh format
 * 2. Optionally refines uniformly
 * 3. Writes in MARS binary format (x.float32, y.float32, z.float32, i0.int32, etc.)
 */

#include "mars_read_mfem_mesh.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstring>

// Write vectors to binary files
template<typename T>
void writeBinaryFile(const std::string& filename, const std::vector<T>& data) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Failed to open file for writing: " + filename);
    }
    file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(T));
    file.close();
    std::cout << "Wrote " << data.size() << " values to " << filename << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <input.mesh> <output_dir> [refine_levels]\n";
        std::cerr << "\nConverts MFEM mesh to MARS binary format\n";
        std::cerr << "  input.mesh    - MFEM mesh file\n";
        std::cerr << "  output_dir    - Output directory for binary files\n";
        std::cerr << "  refine_levels - Optional: number of uniform refinement levels\n";
        std::cerr << "\nExample: " << argv[0] << " beam-tet.mesh mesh_binary 3\n";
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputDir = argv[2];
    int refLevels = (argc == 4) ? std::atoi(argv[3]) : 0;

    std::cout << "Converting MFEM mesh to MARS binary format\n";
    std::cout << "  Input:  " << inputFile << "\n";
    std::cout << "  Output: " << outputDir << "\n";
    std::cout << "  Refine: " << refLevels << " levels\n\n";

    try {
        // Read MFEM mesh (single rank reads all)
        // Note: For serial conversion, use rank=0, numRanks=1 to get full mesh
        auto [nodeCount, elementCount, x_data, y_data, z_data, connectivity, local_to_global] =
            mars::readMFEMMeshWithElementPartitioning<4, float, unsigned>(inputFile, 0, 1);

        std::cout << "Initial mesh: " << elementCount << " elements, " << nodeCount << " nodes\n";

        // Apply refinement if requested
        if (refLevels < 0) {
            // Auto-calculate refinement to reach ~50k elements
            refLevels = (int)std::floor(std::log(50000.0 / elementCount) / std::log(2.0) / 3.0);
            std::cout << "Auto-refinement: " << refLevels << " levels to reach ~50k elements\n";
        }

        for (int l = 0; l < refLevels; l++) {
            mars::uniformRefineMFEMMesh<4, float, unsigned>(x_data, y_data, z_data, connectivity, nodeCount, elementCount);
            std::cout << "  After refinement " << (l+1) << ": " 
                      << elementCount << " elements, " << nodeCount << " nodes\n";
        }

        std::cout << "\nFinal mesh: " << elementCount << " elements, " << nodeCount << " nodes\n";

        // Create output directory (simple approach - assumes it exists or mkdir was called)
        std::cout << "\nWriting binary files to " << outputDir << "/\n";

        // Write coordinates
        writeBinaryFile(outputDir + "/x.float32", x_data);
        writeBinaryFile(outputDir + "/y.float32", y_data);
        writeBinaryFile(outputDir + "/z.float32", z_data);

        // Write connectivity
        // Convert unsigned to int32 for MARS format
        auto write_connectivity_index = [&](int idx, const std::vector<unsigned>& data) {
            std::vector<int32_t> int32_data(data.size());
            for (size_t i = 0; i < data.size(); i++) {
                int32_data[i] = static_cast<int32_t>(data[i]);
            }
            writeBinaryFile(outputDir + "/i" + std::to_string(idx) + ".int32", int32_data);
        };

        write_connectivity_index(0, std::get<0>(connectivity));
        write_connectivity_index(1, std::get<1>(connectivity));
        write_connectivity_index(2, std::get<2>(connectivity));
        write_connectivity_index(3, std::get<3>(connectivity));

        std::cout << "\nConversion complete!\n";
        std::cout << "To use with MARS, run:\n";
        std::cout << "  mpirun -np 1 mars_ex1_poisson --mesh " << outputDir << "\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}
