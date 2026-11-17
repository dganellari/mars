#include "core/mars_tet_red.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace
{
using Point = std::array<double, 3>;
using Tet = std::array<unsigned, 4>;
using Face = std::array<unsigned, 3>;

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

double signed_volume(const Tet& tet, const std::vector<Point>& points)
{
    Point a{}, b{}, c{};
    for (int d = 0; d < 3; ++d)
    {
        a[d] = points[tet[1]][d] - points[tet[0]][d];
        b[d] = points[tet[2]][d] - points[tet[0]][d];
        c[d] = points[tet[3]][d] - points[tet[0]][d];
    }
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
            a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
}

void check_red_template()
{
    std::vector<Point> points{{0.2, 0.3, 0.4}, {2.1, 0.1, 0.3}, {0.4, 1.6, 0.1}, {0.1, 0.2, 1.9}};
    Tet parent{0, 1, 2, 3};
    unsigned permutations = 0;
    do
    {
        unsigned nodes[10] = {parent[0], parent[1], parent[2], parent[3], 4, 5, 6, 7, 8, 9};
        constexpr int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
        for (const auto& edge : edges)
        {
            Point midpoint;
            for (int d = 0; d < 3; ++d)
                midpoint[d] = (points[parent[edge[0]]][d] + points[parent[edge[1]]][d]) / 2;
            points.push_back(midpoint);
        }
        double volume = signed_volume(parent, points), sum = 0;
        std::map<Face, int> faces;
        for (int c = 0; c < 8; ++c)
        {
            Tet child;
            mars::amr::red_child_connectivity(c, nodes, child.data());
            double child_volume = signed_volume(child, points);
            require(std::abs(child_volume / volume - 0.125) < 1.e-13, "Red child orientation or volume");
            sum += child_volume;
            for (int opposite = 0; opposite < 4; ++opposite)
            {
                Face face;
                for (int i = 0, j = 0; i < 4; ++i)
                    if (i != opposite) face[j++] = child[i];
                std::sort(face.begin(), face.end());
                ++faces[face];
            }
        }
        int boundary = 0, interior = 0;
        for (const auto& [face, count] : faces)
        {
            boundary += count == 1;
            interior += count == 2;
            require(count == 1 || count == 2, "Red child face multiplicity");
        }
        require(boundary == 16 && interior == 8, "Red child faces do not match");
        require(std::abs(sum / volume - 1) < 1.e-13, "Red children do not preserve volume");
        points.resize(4);
        ++permutations;
    } while (std::next_permutation(parent.begin(), parent.end()));
    require(permutations == 24, "Missing parent permutations");
}
} // namespace

int main()
{
    try
    {
        check_red_template();
        std::cout << "PASS: 24 parent permutations, 192 child orientations, volumes and face counts\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
