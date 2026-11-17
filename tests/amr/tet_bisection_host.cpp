#include "core/mars_simplex_bisection.hpp"
#include "core/mars_tet_red.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace
{
using Point = std::array<double, 3>;
using Tet   = std::array<unsigned, 4>;
using Edge  = std::array<unsigned, 2>;
using Face  = std::array<unsigned, 3>;

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
            a[2] * (b[0] * c[1] - b[1] * c[0])) /
           6;
}

double affine(const Point& point) { return 2 + 3 * point[0] - 5 * point[1] + 7 * point[2]; }

void check_longest_edge_ties()
{
    using Key                         = uint64_t;
    const std::array<Key, 2> edges[6] = {{19, 7}, {7, 11}, {11, 19}, {3, 19}, {7, 3}, {11, 3}};
    std::array<int, 6> order{0, 1, 2, 3, 4, 5};
    do
    {
        mars::bisection::LongestEdge<Key, double> choice;
        for (int i : order)
            choice.consider(1, edges[i][0], edges[i][1], i);
        require(choice.edge == mars::bisection::Edge<Key>(3, 7), "Longest-edge ties depend on traversal order");
        choice.consider(2, 19, 11, 6);
        require(choice.index == 6, "Global-key tie rule overrode a longer edge");
    } while (std::next_permutation(order.begin(), order.end()));
}

void check_template()
{
    std::vector<Point> points{{0.2, 0.3, 0.4}, {2.1, 0.1, 0.3}, {0.4, 1.6, 0.1}, {0.1, 0.2, 1.9}};
    Tet tet{0, 1, 2, 3};
    int cases = 0;
    do
    {
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b)
            {
                Point midpoint;
                for (int d = 0; d < 3; ++d)
                    midpoint[d] = (points[tet[a]][d] + points[tet[b]][d]) / 2;
                points.push_back(midpoint);
                Tet left, right;
                mars::bisection::child_connectivity<4>(tet, b, 4u, left);
                mars::bisection::child_connectivity<4>(tet, a, 4u, right);
                double volume = signed_volume(tet, points);
                require(std::abs(signed_volume(left, points) / volume - 0.5) < 1e-13, "left child orientation/volume");
                require(std::abs(signed_volume(right, points) / volume - 0.5) < 1e-13,
                        "right child orientation/volume");
                require(std::abs(affine(midpoint) - (affine(points[tet[a]]) + affine(points[tet[b]])) / 2) < 1e-13,
                        "affine midpoint transfer");
                points.pop_back();
                ++cases;
            }
    } while (std::next_permutation(tet.begin(), tet.end()));
    require(cases == 144, "all vertex permutations and edges tested");
}

void check_red_template()
{
    std::vector<Point> points{{0.2, 0.3, 0.4}, {2.1, 0.1, 0.3}, {0.4, 1.6, 0.1}, {0.1, 0.2, 1.9}};
    Tet parent{0, 1, 2, 3};
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
                for (int i = 0, j = 0; i < 4; ++i) if (i != opposite) face[j++] = child[i];
                std::sort(face.begin(), face.end());
                ++faces[face];
            }
        }
        int boundary = 0, interior = 0;
        for (const auto& [face, count] : faces)
        {
            boundary += count == 1;
            interior += count == 2;
            require(count == 1 || count == 2, "Red children overlap");
        }
        require(boundary == 16 && interior == 8, "Red child faces do not match");
        require(std::abs(sum / volume - 1) < 1.e-13, "Red children do not preserve volume");
        points.resize(4);
    } while (std::next_permutation(parent.begin(), parent.end()));
}

struct Mesh
{
    std::vector<Point> points{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    std::vector<Tet> tets{{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    std::vector<double> values;

    Mesh()
    {
        for (const auto& point : points)
            values.push_back(affine(point));
        for (auto& tet : tets)
            if (signed_volume(tet, points) < 0) std::swap(tet[0], tet[1]);
    }

    std::vector<Edge> edges() const
    {
        std::vector<Edge> result;
        for (const auto& tet : tets)
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b)
                    result.push_back({std::min(tet[a], tet[b]), std::max(tet[a], tet[b])});
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        return result;
    }

    void split(Edge edge)
    {
        std::vector<Tet> children;
        Point midpoint;
        for (int d = 0; d < 3; ++d)
            midpoint[d] = (points[edge[0]][d] + points[edge[1]][d]) / 2;
        unsigned id     = points.size();
        size_t incident = 0;
        for (const auto& tet : tets)
        {
            auto a = std::find(tet.begin(), tet.end(), edge[0]);
            auto b = std::find(tet.begin(), tet.end(), edge[1]);
            if (a == tet.end() || b == tet.end())
                children.push_back(tet);
            else
            {
                Tet left, right;
                mars::bisection::child_connectivity<4>(tet, int(b - tet.begin()), id, left);
                mars::bisection::child_connectivity<4>(tet, int(a - tet.begin()), id, right);
                children.push_back(left);
                children.push_back(right);
                ++incident;
            }
        }
        require(incident > 0, "selected edge has an incident element");
        require(children.size() == tets.size() + incident, "whole edge star is bisected");
        points.push_back(midpoint);
        values.push_back((values[edge[0]] + values[edge[1]]) / 2);
        tets = std::move(children);
    }

    void check() const
    {
        std::map<Face, unsigned> faces;
        double volume = 0;
        for (const auto& tet : tets)
        {
            double local = signed_volume(tet, points);
            require(local > 0, "positive child orientation");
            volume += local;
            for (int omit = 0; omit < 4; ++omit)
            {
                Face face;
                int j = 0;
                for (int i = 0; i < 4; ++i)
                    if (i != omit) face[j++] = tet[i];
                std::sort(face.begin(), face.end());
                ++faces[face];
            }
        }
        require(std::abs(volume - 1) < 1e-12, "total cube volume");
        for (const auto& [face, count] : faces)
        {
            bool boundary = false;
            for (int d = 0; d < 3; ++d)
                for (double side : {0., 1.})
                    boundary = boundary ||
                               (points[face[0]][d] == side && points[face[1]][d] == side && points[face[2]][d] == side);
            require(count == (boundary ? 1u : 2u), "matching complete interior faces, including edge stars");
        }
        for (size_t i = 0; i < points.size(); ++i)
            require(std::abs(values[i] - affine(points[i])) < 1e-13, "repeated affine parentage transfer");
    }
};

void check_meshes()
{
    Mesh base;
    base.check();
    const auto initial_edges = base.edges();
    for (const auto& edge : initial_edges)
    {
        Mesh mesh;
        mesh.split(edge);
        mesh.check();
    }
    for (bool reverse : {false, true})
    {
        Mesh mesh;
        auto edges = initial_edges;
        if (reverse) std::reverse(edges.begin(), edges.end());
        for (const auto& edge : edges)
        {
            mesh.split(edge);
            mesh.check();
        }
        for (int level = 0; level < 3; ++level)
        {
            auto current_edges = mesh.edges();
            for (size_t i = 0; i < current_edges.size(); i += 7)
                mesh.split(current_edges[i]);
            mesh.check();
        }
    }
}
} // namespace

int main()
{
    try
    {
        check_longest_edge_ties();
        check_template();
        check_red_template();
        check_meshes();
        std::cout
            << "PASS: 144 oriented edge templates, complete edge stars, matching faces, volume and affine transfer\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
