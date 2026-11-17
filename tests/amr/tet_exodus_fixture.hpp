#pragma once

#include <algorithm>
#include <array>
#include <vector>

namespace amr_test
{
struct ExodusFixture
{
    std::vector<std::array<double, 3>> points;
    std::vector<std::array<int, 4>> cells;
};

inline ExodusFixture exodus_fixture()
{
    constexpr int n = 4, nodes = 125, cells = 384;
    ExodusFixture mesh;
    mesh.points.resize(nodes);
    mesh.cells.resize(cells);
    auto node = [](int x, int y, int z) { return (37 * (x + 5 * (y + 5 * z)) + 11) % nodes; };
    for (int z = 0; z <= n; ++z)
        for (int y = 0; y <= n; ++y)
            for (int x = 0; x <= n; ++x)
                mesh.points[node(x, y, z)] = {double(x) / n, double(y) / n, double(z) / n};

    int element = 0;
    for (int z = 0; z < n; ++z)
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x)
            {
                // The six monotone paths share face diagonals between neighbouring cubes.
                std::array<int, 3> axes{0, 1, 2};
                do
                {
                    std::array<int, 3> corner{x, y, z};
                    std::array<int, 4> tet{node(x, y, z)};
                    for (int j = 0; j < 3; ++j)
                    {
                        ++corner[axes[j]];
                        tet[j + 1] = node(corner[0], corner[1], corner[2]);
                    }
                    int inversions = (axes[0] > axes[1]) + (axes[0] > axes[2]) + (axes[1] > axes[2]);
                    if (inversions % 2) std::swap(tet[0], tet[1]);
                    // Even permutations preserve orientation and exercise every Exodus side ordinal.
                    std::rotate(tet.begin(), tet.begin() + element % 3, tet.begin() + 3);
                    mesh.cells[(73 * element + 19) % cells] = tet;
                    ++element;
                } while (std::next_permutation(axes.begin(), axes.end()));
            }
    return mesh;
}
} // namespace amr_test
