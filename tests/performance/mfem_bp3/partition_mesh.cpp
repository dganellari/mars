#include "mfem.hpp"
#include "exodus_hex_mesh.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sys/resource.h>

// Offline partitioning: mesh topology is held once, never by every solver rank.
static int partition(int argc, char** argv)
{
    if (argc != 4)
    {
        std::cerr << "usage: mfem_partition input.e output-prefix parts\n";
        return 1;
    }
    const int parts = std::stoi(argv[3]);
    if (parts < 1) return 1;
    for (int rank = 0; rank < parts; ++rank)
        if (std::ifstream(mfem::MakeParFilename(argv[2], rank)))
            throw std::runtime_error("Output partition already exists");
    if (std::ifstream(std::string(argv[2]) + ".parts"))
        throw std::runtime_error("Output metadata already exists");
    const auto start = std::chrono::steady_clock::now();
    auto input       = bp3_mesh::read(argv[1]);
    mfem::Mesh& mesh = *input;
    if (mesh.Dimension() != 3 || mesh.SpaceDimension() != 3 || mesh.GetNodes() || mesh.GetNE() < parts ||
        mesh.CheckElementOrientation(false))
        return 1;
    if (!mesh.GetNBE())
    {
        std::cerr << "Missing exterior boundary\n";
        return 1;
    }
    std::vector<int> partition(mesh.GetNE());
    {
        mfem::Array<int> position;
        mesh.GetHilbertElementOrdering(position);
        // MFEM returns the position of each original element in Hilbert order.
        std::vector<bool> seen(mesh.GetNE(), false);
        for (int e = 0; e < mesh.GetNE(); ++e)
        {
            const int pos = position[e];
            if (pos < 0 || pos >= mesh.GetNE() || seen[pos])
                throw std::runtime_error("Invalid MFEM Hilbert permutation");
            seen[pos] = true;
            partition[e] = static_cast<int>(int64_t(pos) * parts / mesh.GetNE());
        }
    }
    mfem::MeshPartitioner partitioner(mesh, parts, partition.data());
    mfem::MeshPart part;
    for (int rank = 0; rank < parts; ++rank)
    {
        partitioner.ExtractPart(rank, part);
        std::ofstream out(mfem::MakeParFilename(argv[2], rank));
        out << std::setprecision(17);
        part.Print(out);
        out.close();
        if (!out) return 1;
    }
    std::ofstream metadata(std::string(argv[2]) + ".parts");
    metadata << parts << ' ' << mesh.GetNE() << '\n';
    metadata.close();
    if (!metadata) return 1;
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    double peak = usage.ru_maxrss;
#ifndef __APPLE__
    peak *= 1024;
#endif
    std::cout << "MFEM HILBERT PARTITION PASS elements " << mesh.GetNE() << " vertices " << mesh.GetNV() << " boundary_faces "
              << mesh.GetNBE() << " parts " << parts << " peak_host_bytes " << std::fixed << std::setprecision(0)
              << peak << " elapsed_seconds " << std::setprecision(6)
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << '\n';
    return 0;
}

int main(int argc, char** argv)
{
    try
    {
        return partition(argc, argv);
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
