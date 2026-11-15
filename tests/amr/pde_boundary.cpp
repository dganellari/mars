#include "public_tet_mesh.hpp"
#include "examples/distributed/unstructured/mars_amr_galerkin_boundary.hpp"
#include <fstream>
#include <iomanip>
#include <set>
#include <iostream>
using Key  = uint64_t;
using Exec = mars::amr::LeppHostExecution;
struct NodeKey
{
    Key operator()(double x, double y, double z) const
    {
        auto fixed = [](double a) { return Key((a + 8) * 65536 + .5); };
        return fixed(x) | (fixed(y) << 21) | (fixed(z) << 42);
    }
};
struct Monomial
{
    double u(double x, double, double) const { return x * x; }
};
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try
    {
        auto input = amr_test::read_public_tet(MARS_PUBLIC_TET_MESH);
        auto xyz   = input.mesh.xyz;
        std::vector<Key> keys;
        for (size_t i = 0; i < xyz[0].size(); ++i)
            keys.push_back(NodeKey{}(xyz[0][i], xyz[1][i], xyz[2][i]));
        std::map<std::array<Key, 3>, int> labelled;
        for (const auto& [tag, faces] : input.boundary)
            for (auto face : faces)
            {
                std::sort(face.begin(), face.end());
                labelled[face] = tag;
            }
        std::array<std::vector<Key>, 4> conn;
        std::vector<mars::amr::LeppTaggedFace> labels;
        for (size_t e = 0; e < input.mesh.conn[0].size(); ++e)
            if (int(e % ranks) == rank)
            {
                size_t local = conn[0].size();
                for (int v = 0; v < 4; ++v)
                    conn[v].push_back(input.mesh.conn[v][e]);
                for (int f = 0; f < 4; ++f)
                {
                    std::array<Key, 3> face;
                    for (int v = 0, j = 0; v < 4; ++v)
                        if (v != f) face[j++] = conn[v][local];
                    std::sort(face.begin(), face.end());
                    if (labelled.count(face)) labels.push_back({local, f, labelled.at(face)});
                }
            }
        auto owner = [ranks](Key k) { return int(k % ranks); }; // Reference routing, unrelated to the device partition.
        for (int level = 0; level <= 4; ++level)
        {
            std::array<const Key*, 4> c{conn[0].data(), conn[1].data(), conn[2].data(), conn[3].data()};
            auto faces = mars::amr_pde::check_boundary<Exec, Key>(keys.data(), c, conn[0].size(), labels, owner,
                                                                  MPI_COMM_WORLD, false);
            for (int failure = 0; failure < 2; ++failure)
            {
                auto broken = labels;
                if (failure == 0 && rank == 0 && !broken.empty()) broken.pop_back();
                if (failure == 1)
                {
                    std::set<std::pair<size_t, int>> tagged;
                    for (auto f : labels)
                        tagged.insert({f.element, f.face});
                    for (size_t e = 0; e < conn[0].size(); ++e)
                        for (int f = 0; f < 4; ++f)
                            if (!tagged.count({e, f}))
                            {
                                broken.push_back({e, f, 999});
                                e = conn[0].size();
                                break;
                            }
                }
                bool rejected = false;
                try
                {
                    mars::amr_pde::check_boundary<Exec, Key>(keys.data(), c, conn[0].size(), broken, owner,
                                                             MPI_COMM_WORLD, false);
                }
                catch (const std::runtime_error&)
                {
                    rejected = true;
                }
                if (!rejected) throw std::runtime_error("Missing or interior boundary label accepted");
            }
            if (ranks == 1 && argc == 2)
            {
                std::ofstream out(std::string(argv[1]) + std::to_string(level) + ".json");
                out << std::setprecision(17);
                out << "{\"points\":[";
                for (size_t i = 0; i < keys.size(); ++i)
                {
                    if (i) out << ",";
                    out << "[" << xyz[0][i] << "," << xyz[1][i] << "," << xyz[2][i] << "]";
                }
                out << "],\"cells\":[";
                for (size_t e = 0; e < conn[0].size(); ++e)
                {
                    if (e) out << ",";
                    out << "[" << conn[0][e] << "," << conn[1][e] << "," << conn[2][e] << "," << conn[3][e] << "]";
                }
                out << "],\"boundary\":[";
                for (size_t i = 0; i < labels.size(); ++i)
                {
                    if (i) out << ",";
                    out << "[" << labels[i].element << "," << labels[i].face << "]";
                }
                out << "]}\n";
            }
            if (rank == 0) std::cout << "PDE boundary level=" << level << " faces=" << faces << " PASS\n";
            if (level == 4) break;
            std::vector<uint8_t> marks(conn[0].size(), 1);
            mars::amr::TetLeppView<Key, double> view{
                {c[0], c[1], c[2], c[3]}, keys.data(), xyz[0].data(), xyz[1].data(), xyz[2].data(), marks.data(),
                conn[0].size(),           keys.size()};
            auto next = mars::amr::refine_tet_lepp<Exec>(view, NodeKey{}, owner, MPI_COMM_WORLD, false, labels);
            conn      = std::move(next.conn);
            keys      = std::move(next.keys);
            labels    = std::move(next.tagged_faces);
            xyz       = {std::move(next.x), std::move(next.y), std::move(next.z)};
        }
        double x[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, u[4] = {};
        double integral = mars::amr_pde::tet_error(x, u, Monomial{});
        if (std::abs(integral - 1. / 210) > 1e-16) throw std::runtime_error("Quartic integration failed");
        MPI_Finalize();
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
