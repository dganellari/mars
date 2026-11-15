#pragma once
#include "public_tet_mesh.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_tet_red_boundary.hpp"
#include "core/mars_tet_red.hpp"
#include <fstream>
#include <iomanip>
#include <set>

namespace amr_test
{
template<class Execution>
void check_red_boundary(MPI_Comm comm, bool staged, const std::string& export_prefix = {})
{
    using namespace mars::amr;
    using Key  = uint64_t;
    using Face = LeppFaceKey<Key>;
    using Edge = std::pair<Key, Key>;
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    auto input = read_public_tet(MARS_PUBLIC_TET_MESH);
    auto xyz   = input.mesh.xyz;
    std::vector<std::array<Key, 4>> cells;
    for (size_t e = 0; e < input.mesh.conn[0].size(); ++e)
        cells.push_back({input.mesh.conn[0][e], input.mesh.conn[1][e], input.mesh.conn[2][e], input.mesh.conn[3][e]});
    std::map<Face, int> boundary;
    for (const auto& [tag, faces] : input.boundary)
        for (auto f : faces)
            boundary.emplace(Face(f[0], f[1], f[2]), tag);
    std::map<Key, Edge> ancestry;
    auto upload = []<class T>(const std::vector<T>& h)
    {
        typename Execution::template Vector<T> d(h.size());
        Execution::upload(h.data(), h.size(), Execution::data(d));
        return d;
    };
    auto owner = [=] HOST_DEVICE_FUN(Key k) { return int(k % ranks); };
    for (int level = 0; level < 3; ++level)
    {
        std::vector<uint8_t> marks(cells.size());
        std::set<uint64_t> edges;
        for (size_t e = 0; e < cells.size(); ++e)
        {
            marks[e] = level == 1 || e % (level ? 5 : 3) == 0;
            if (marks[e])
                for (int a = 0; a < 4; ++a)
                    for (int b = a + 1; b < 4; ++b)
                        edges.insert((std::max(cells[e][a], cells[e][b]) << 32) | std::min(cells[e][a], cells[e][b]));
        }
        const size_t old_nodes = xyz[0].size();
        std::vector<uint64_t> sorted_edges(edges.begin(), edges.end());
        RedBoundaryMidpoint<Key> midpoint{sorted_edges.data(), sorted_edges.size(), Key(old_nodes)};
        for (auto edge : sorted_edges)
        {
            Key a = edge & UINT32_MAX, b = edge >> 32, node = xyz[0].size();
            ancestry[node] = {a, b};
            for (int d = 0; d < 3; ++d)
                xyz[d].push_back(.5 * (xyz[d][a] + xyz[d][b]));
        }
        std::array<std::vector<Key>, 4> local_conn;
        std::vector<uint8_t> local_marks;
        std::vector<LeppTaggedFace> labels;
        for (size_t e = 0; e < cells.size(); ++e)
            if (int(e % ranks) == rank)
            {
                size_t local = local_marks.size();
                local_marks.push_back(marks[e]);
                for (int v = 0; v < 4; ++v)
                    local_conn[v].push_back(cells[e][v]);
                for (int f = 0; f < 4; ++f)
                {
                    Key v[3];
                    for (int c = 0, j = 0; c < 4; ++c)
                        if (c != f) v[j++] = cells[e][c];
                    auto label = boundary.find(Face(v[0], v[1], v[2]));
                    if (label != boundary.end()) labels.push_back({local, f, label->second});
                }
            }
        auto d_labels = upload(labels);
        auto d_marks  = upload(local_marks);
        auto d_edges  = upload(sorted_edges);
        std::array<typename Execution::template Vector<Key>, 4> d_conn;
        std::array<const Key*, 4> pointers;
        for (int v = 0; v < 4; ++v)
        {
            d_conn[v]   = upload(local_conn[v]);
            pointers[v] = Execution::data(d_conn[v]);
        }
        auto saved =
            save_red_boundary<Execution, Key>(d_labels, pointers, Execution::data(d_marks),
                                              {Execution::data(d_edges), sorted_edges.size(), Key(old_nodes)}, comm);
        // All ranks must reject a missing midpoint, including ranks with no labelled face.
        bool rejected = false;
        try
        {
            save_red_boundary<Execution, Key>(d_labels, pointers, Execution::data(d_marks),
                                              {Execution::data(d_edges), 0, Key(old_nodes)}, comm);
        }
        catch (const std::runtime_error&)
        {
            rejected = true;
        }
        require_collectively(rejected, "Missing boundary midpoint was accepted", comm);
        std::vector<std::array<Key, 4>> children;
        std::map<Face, int> next_boundary;
        for (size_t e = 0; e < cells.size(); ++e)
        {
            auto parent   = cells[e];
            Key nodes[10] = {parent[0], parent[1], parent[2], parent[3]};
            if (marks[e])
                for (int a = 0, j = 4; a < 4; ++a)
                    for (int b = a + 1; b < 4; ++b)
                        nodes[j++] = midpoint(parent[a], parent[b]);
            for (int child = 0; child < (marks[e] ? 8 : 1); ++child)
            {
                std::array<Key, 4> cell = parent;
                if (marks[e]) red_child_connectivity(child, nodes, cell.data());
                children.push_back(cell);
                // Independent incidence oracle: all corners of a child face descend from one parent face.
                for (int f = 0; f < 4; ++f)
                {
                    Key v[3];
                    for (int c = 0, j = 0; c < 4; ++c)
                        if (c != f) v[j++] = cell[c];
                    for (int pf = 0; pf < 4; ++pf)
                    {
                        Key pv[3];
                        for (int c = 0, j = 0; c < 4; ++c)
                            if (c != pf) pv[j++] = parent[c];
                        auto label = boundary.find(Face(pv[0], pv[1], pv[2]));
                        if (label == boundary.end()) continue;
                        std::set<Key> corners(pv, pv + 3);
                        bool on = true;
                        for (Key node : v)
                            if (node < old_nodes)
                                on &= corners.count(node) != 0;
                            else
                            {
                                auto [a, b] = ancestry.at(node);
                                on &= corners.count(a) && corners.count(b);
                            }
                        if (on) next_boundary.emplace(Face(v[0], v[1], v[2]), label->second);
                    }
                }
            }
        }
        // Simulate rebuild identity merging; a neighbouring coarse cell may emit an existing midpoint.
        std::map<std::array<double, 3>, Key> identity;
        std::vector<Key> keys(xyz[0].size());
        for (size_t i = 0; i < keys.size(); ++i)
        {
            auto [it, inserted] = identity.emplace(std::array<double, 3>{xyz[0][i], xyz[1][i], xyz[2][i]}, Key(i));
            keys[i]             = it->second;
            if (!inserted) ancestry.erase(i);
        }
        for (auto& cell : children)
            for (auto& node : cell)
                node = keys[node];
        boundary.clear();
        for (auto [face, tag] : next_boundary)
            boundary.emplace(Face(keys[face.keys[0]], keys[face.keys[1]], keys[face.keys[2]]), tag);
        auto d_keys  = upload(keys);
        auto records = Execution::data(saved);
        auto k       = Execution::data(d_keys);
        Execution::for_each(saved.size(),
                            [=] HOST_DEVICE_FUN(size_t i)
                            {
                                auto f          = records[i].face;
                                records[i].face = Face(k[f.keys[0]], k[f.keys[1]], k[f.keys[2]]);
                            });
        cells = std::move(children);
        std::vector<LeppTaggedFace> expected;
        for (auto& c : local_conn)
            c.clear();
        for (size_t e = 0; e < cells.size(); ++e)
            if (int((e + level + 1) % ranks) == rank)
            {
                size_t local = local_conn[0].size();
                for (int v = 0; v < 4; ++v)
                    local_conn[v].push_back(cells[e][v]);
                for (int f = 0; f < 4; ++f)
                {
                    Key v[3];
                    for (int c = 0, j = 0; c < 4; ++c)
                        if (c != f) v[j++] = cells[e][c];
                    auto label = boundary.find(Face(v[0], v[1], v[2]));
                    if (label != boundary.end()) expected.push_back({local, f, label->second});
                }
            }
        for (int v = 0; v < 4; ++v)
        {
            d_conn[v]   = upload(local_conn[v]);
            pointers[v] = Execution::data(d_conn[v]);
        }
        auto matched = match_lepp_faces<Execution, Key>(saved, Execution::data(d_keys), pointers, local_conn[0].size(),
                                                        owner, comm, staged);
        std::vector<LeppTaggedFace> got(matched.size());
        Execution::download(Execution::data(matched), got.size(), got.data());
        auto canonical = [](const auto& list)
        {
            std::set<std::tuple<size_t, int, int>> out;
            for (auto f : list)
                out.emplace(f.element, f.face, f.tag);
            return out;
        };
        require_collectively(got.size() == expected.size() && canonical(got) == canonical(expected),
                             "Red side-set inheritance/repartition mismatch", comm);
        if (rank == 0 && !export_prefix.empty())
        {
            std::ofstream out(export_prefix + std::to_string(level) + ".json");
            out << std::setprecision(17) << "{\"points\":[";
            for (size_t i = 0; i < keys.size(); ++i)
            {
                if (i) out << ',';
                out << '[' << xyz[0][i] << ',' << xyz[1][i] << ',' << xyz[2][i] << ']';
            }
            out << "],\"cells\":[";
            for (size_t e = 0; e < cells.size(); ++e)
            {
                if (e) out << ',';
                auto c = cells[e];
                out << '[' << c[0] << ',' << c[1] << ',' << c[2] << ',' << c[3] << ']';
            }
            out << "],\"boundary\":[";
            bool comma = false;
            for (size_t e = 0; e < cells.size(); ++e)
                for (int f = 0; f < 4; ++f)
                {
                    Key v[3];
                    for (int c = 0, j = 0; c < 4; ++c)
                        if (c != f) v[j++] = cells[e][c];
                    if (boundary.count(Face(v[0], v[1], v[2])))
                    {
                        if (comma) out << ',';
                        comma = true;
                        out << '[' << e << ',' << f << ']';
                    }
                }
            out << "],\"ancestry\":[";
            comma = false;
            for (auto [node, edge] : ancestry)
            {
                if (comma) out << ',';
                comma = true;
                out << '[' << node << ',' << edge.first << ',' << edge.second << ']';
            }
            out << "]}\n";
        }
        if (rank == 0)
            std::cout << "red boundary level=" << level << " cells=" << cells.size() << " labels=" << boundary.size()
                      << " staged=" << staged << " PASS\n";
    }
}
} // namespace amr_test
