#pragma once
#include <array>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cmath>

namespace amr_test
{
template<int N>
struct GradedMesh
{
    std::array<std::vector<double>, 3> xyz;
    std::array<std::vector<uint64_t>, N> conn;
};

// A conforming tensor mesh with a 90:1 spacing change, sheared away from the coordinate axes.
template<int N>
GradedMesh<N> graded_mesh()
{
    GradedMesh<N> mesh;
    constexpr double grid[5] = {0, .005, .02, .1, 1};
    for (int z = 0; z < 5; ++z)
        for (int y = 0; y < 5; ++y)
            for (int x = 0; x < 5; ++x)
            {
                double a = grid[x], b = y * .25, c = z * .25;
                mesh.xyz[0].push_back(a + .2 * b);
                mesh.xyz[1].push_back(b + .1 * c);
                mesh.xyz[2].push_back(c + .05 * a);
            }
    auto id                 = [](int x, int y, int z) { return uint64_t(x + 5 * (y + 5 * z)); };
    constexpr int tet[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    for (int z = 0; z < 4; ++z)
        for (int y = 0; y < 4; ++y)
            for (int x = 0; x < 4; ++x)
            {
                uint64_t h[8] = {id(x, y, z),     id(x + 1, y, z),     id(x + 1, y + 1, z),     id(x, y + 1, z),
                                 id(x, y, z + 1), id(x + 1, y, z + 1), id(x + 1, y + 1, z + 1), id(x, y + 1, z + 1)};
                if constexpr (N == 8)
                    for (int c = 0; c < 8; ++c)
                        mesh.conn[c].push_back(h[c]);
                else
                    for (const auto& t : tet)
                        for (int c = 0; c < 4; ++c)
                            mesh.conn[c].push_back(h[t[c]]);
            }
    return mesh;
}

// Seven tets triangulate each root tetrahedron. The second interior vertex is close to AB but not on it.
// The old geometric detector incorrectly constrains it to AB in the neighbouring tetrahedron ABPD.
inline GradedMesh<4> thin_mesh()
{
    GradedMesh<4> mesh;
    const double p[3] = {.5, .01, .01};
    std::array<std::array<double, 3>, 6> points{{{0, 0, 0},
                                                 {1, 0, 0},
                                                 {0, 1, 0},
                                                 {0, 0, 1},
                                                 {p[0], p[1], p[2]},
                                                 {.995 * p[0], .995 * p[1] + .005 / 3, .995 * p[2] + .005 / 3}}};
    const int cells[7][4] = {{4, 1, 2, 3}, {0, 1, 4, 3}, {0, 1, 2, 4}, {5, 4, 2, 3},
                             {0, 5, 2, 3}, {0, 4, 5, 3}, {0, 4, 2, 5}};
    for (int block = 0; block < 64; ++block)
    {
        uint64_t base = mesh.xyz[0].size();
        for (auto point : points)
            for (int d = 0; d < 3; ++d)
                mesh.xyz[d].push_back(point[d] + 2 * ((block >> (2 * d)) & 3));
        for (auto cell : cells)
        {
            int v[4] = {cell[0], cell[1], cell[2], cell[3]};
            double a[3], b[3], c[3];
            for (int d = 0; d < 3; ++d)
            {
                a[d] = points[v[1]][d] - points[v[0]][d];
                b[d] = points[v[2]][d] - points[v[0]][d];
                c[d] = points[v[3]][d] - points[v[0]][d];
            }
            double det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
                         a[2] * (b[0] * c[1] - b[1] * c[0]);
            if (det < 0) std::swap(v[1], v[2]);
            for (int j = 0; j < 4; ++j)
                mesh.conn[j].push_back(base + v[j]);
        }
    }
    return mesh;
}
} // namespace amr_test
