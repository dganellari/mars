#pragma once
#include "graded_mesh.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>

namespace amr_test
{
struct PublicTetMesh
{
    GradedMesh<4> mesh;
    std::map<int, std::vector<std::array<uint64_t, 3>>> boundary;
};
inline PublicTetMesh read_public_tet(const std::string& path)
{
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot read public MFEM fixture");
    std::ostringstream clean;
    std::string line;
    while (std::getline(file, line))
    {
        auto at = line.find('#');
        clean << line.substr(0, at) << '\n';
    }
    std::istringstream in(clean.str());
    std::string word;
    PublicTetMesh result;
    bool dimension = false, elements = false, boundary = false, vertices = false;
    while (in >> word)
    {
        if (word == "dimension")
        {
            int d;
            in >> d;
            if (d != 3) throw std::runtime_error("Expected 3D fixture");
            dimension = true;
        }
        else if (word == "elements")
        {
            size_t n;
            in >> n;
            for (size_t i = 0; i < n; ++i)
            {
                int tag, type;
                in >> tag >> type;
                if (type != 4 || tag != 1) throw std::runtime_error("Expected single-block Tet4 fixture");
                for (auto& c : result.mesh.conn)
                {
                    uint64_t v;
                    in >> v;
                    c.push_back(v);
                }
            }
            elements = true;
        }
        else if (word == "boundary")
        {
            size_t n;
            in >> n;
            for (size_t i = 0; i < n; ++i)
            {
                int tag, type;
                std::array<uint64_t, 3> v;
                in >> tag >> type >> v[0] >> v[1] >> v[2];
                if (type != 2) throw std::runtime_error("Expected triangle boundary");
                result.boundary[tag].push_back(v);
            }
            boundary = true;
        }
        else if (word == "vertices")
        {
            size_t n;
            int d;
            in >> n >> d;
            if (d != 3) throw std::runtime_error("Expected 3D vertices");
            for (size_t i = 0; i < n; ++i)
                for (auto& x : result.mesh.xyz)
                {
                    double v;
                    in >> v;
                    x.push_back(v);
                }
            vertices = true;
        }
        else if (word == "nodes")
            throw std::runtime_error("Curved geometry is outside this Tet4 fixture");
    }
    if (!dimension || !elements || !boundary || !vertices) throw std::runtime_error("Incomplete MFEM fixture");
    const auto& xyz = result.mesh.xyz;
    auto& conn      = result.mesh.conn;
    for (size_t e = 0; e < conn[0].size(); ++e)
    {
        double a[3], b[3], c[3];
        for (int j = 0; j < 4; ++j)
            if (conn[j][e] >= xyz[0].size()) throw std::runtime_error("Invalid fixture corner");
        for (int d = 0; d < 3; ++d)
        {
            a[d] = xyz[d][conn[1][e]] - xyz[d][conn[0][e]];
            b[d] = xyz[d][conn[2][e]] - xyz[d][conn[0][e]];
            c[d] = xyz[d][conn[3][e]] - xyz[d][conn[0][e]];
        }
        double det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
                     a[2] * (b[0] * c[1] - b[1] * c[0]);
        if (!std::isfinite(det) || det == 0) throw std::runtime_error("Degenerate fixture tet");
        if (det < 0) std::swap(conn[1][e], conn[2][e]);
    }
    return result;
}
} // namespace amr_test
