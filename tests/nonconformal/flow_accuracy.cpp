#include "examples/distributed/unstructured/mars_nonconformal_flow_runtime.hpp"
#include "channel_mesh.hpp"
#include <tuple>

// Shared and matched cases use identical cells; only the interface identification changes.
static Mesh channel_mesh(int n, const std::string& mode)
{
    Mesh m = mesh(n, 1., n, mode == "nonmatching" ? 2 * n : n);
    for (auto& x : m.x)
        x.x += 1;
    if (mode == "shared")
    {
        std::map<std::tuple<double, double, double>, int> ids;
        std::vector<int> map(m.x.size()), side;
        std::vector<Vec3<double>> coordinates;
        for (size_t i = 0; i < m.x.size(); ++i)
        {
            auto x              = m.x[i];
            auto [it, inserted] = ids.emplace(std::make_tuple(x.x, x.y, x.z), int(coordinates.size()));
            if (inserted)
            {
                coordinates.push_back(x);
                side.push_back(m.side[i]);
            }
            map[i] = it->second;
        }
        for (auto& id : m.cells)
            id = map[id];
        m.x    = std::move(coordinates);
        m.side = std::move(side);
        m.tags.clear();
        m.dofs.resize(m.x.size());
        std::iota(m.dofs.begin(), m.dofs.end(), 0);
    }
    for (int e = 0; e < int(m.elements.size()); ++e)
        for (int opp = 0; opp < 4; ++opp)
        {
            bool planes[6]{true, true, true, true, true, true};
            for (int k = 0; k < 4; ++k)
                if (k != opp)
                {
                    const auto& x = m.x[m.cells[4 * e + k]];
                    planes[0]     = planes[0] && x.x == 0;
                    planes[1]     = planes[1] && x.x == 2;
                    planes[2]     = planes[2] && x.y == 0;
                    planes[3]     = planes[3] && x.y == 1;
                    planes[4]     = planes[4] && x.z == 0;
                    planes[5]     = planes[5] && x.z == 1;
                }
            for (int j = 0; j < 6; ++j)
                if (planes[j]) m.outer.push_back({e, opp, j < 2 ? j : 2});
        }
    return m;
}
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankId);
    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
    int status = 1;
    try
    {
        if (argc != 4) throw std::runtime_error("expected resolution shared|matched|nonmatching output-prefix");
        size_t used{};
        int n            = std::stoi(argv[1], &used);
        std::string mode = argv[2], prefix = argv[3];
        if (used != std::string(argv[1]).size() || n < 2 || n > 32 ||
            (mode != "shared" && mode != "matched" && mode != "nonmatching"))
            throw std::runtime_error("invalid accuracy case");
        auto m = channel_mesh(n, mode);
        ncflow::Channel flow(m, false, .1, .1);
        bool ok = flow.solve(m, 3000, 1e-9);
        if (ok)
        {
            flow.write(m, prefix);
            if (rankId == 0)
            {
                std::ofstream out(prefix + ".cells.csv");
                out << "element,n0,n1,n2,n3\n";
                for (size_t e = 0; e < m.elements.size(); ++e)
                    out << e << ',' << m.cells[4 * e] + 1 << ',' << m.cells[4 * e + 1] + 1 << ','
                        << m.cells[4 * e + 2] + 1 << ',' << m.cells[4 * e + 3] + 1 << '\n';
                out.close();
                if (!out) throw std::runtime_error("accuracy connectivity output failed");
                std::printf("PASS accuracy solve: n=%d mode=%s ranks=%d\n", n, mode.c_str(), rankCount);
            }
        }
        status = ok ? 0 : 1;
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "rank %d: %s\n", rankId, e.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
    return status;
}
