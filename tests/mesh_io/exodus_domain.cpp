#include "domain.hpp"
#include <netcdf.h>
#include <filesystem>
#include <set>
#include <cstdlib>

namespace
{
void require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}
void nc(int result)
{
    if (result != NC_NOERR) throw std::runtime_error(nc_strerror(result));
}
template<int N>
constexpr auto boundary_faces()
{
    if constexpr (N == 4)
        return std::array<std::array<int, 3>, 4>{{{0, 1, 3}, {1, 2, 3}, {0, 3, 2}, {0, 2, 1}}};
    else
        return std::array<std::array<int, 4>, 6>{
            {{0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}, {0, 3, 2, 1}, {4, 5, 6, 7}}};
}
const std::array<std::string, 8> names = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax", "empty", "single_face"};
template<int N>
struct Fixture
{
    static constexpr auto face_nodes = boundary_faces<N>();
    std::array<std::array<double, 3>, 125> points;
    std::vector<std::array<int, N>> cells = std::vector<std::array<int, N>>(N == 8 ? 64 : 384);
    std::array<std::vector<int>, 8> elements, sides;
    Fixture()
    {
        auto node = [](int x, int y, int z) { return (37 * (x + 5 * (y + 5 * z)) + 11) % 125; };
        for (int z = 0; z <= 4; ++z)
            for (int y = 0; y <= 4; ++y)
                for (int x = 0; x <= 4; ++x)
                    points[node(x, y, z)] = {x / 4., y / 4., z / 4.};
        for (int z = 0; z < 4; ++z)
            for (int y = 0; y < 4; ++y)
                for (int x = 0; x < 4; ++x)
                {
                    std::array<int, 8> hex = {
                        node(x, y, z),     node(x + 1, y, z),     node(x + 1, y + 1, z),     node(x, y + 1, z),
                        node(x, y, z + 1), node(x + 1, y, z + 1), node(x + 1, y + 1, z + 1), node(x, y + 1, z + 1)};
                    int cube = (17 * (x + 4 * (y + 4 * z)) + 11) % 64;
                    if constexpr (N == 8)
                        cells[cube] = hex;
                    else
                    {
                        constexpr int split[6][4] = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
                                                     {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
                        for (int t = 0; t < 6; ++t)
                            for (int c = 0; c < 4; ++c)
                                cells[6 * cube + t][c] = hex[split[t][c]];
                    }
                }
        for (size_t e = 0; e < cells.size(); ++e)
            for (size_t s = 0; s < face_nodes.size(); ++s)
                for (int tag = 0; tag < 6; ++tag)
                {
                    bool on = true;
                    for (int c : face_nodes[s])
                        on &= points[cells[e][c]][tag / 2] == tag % 2;
                    if (on)
                    {
                        elements[tag].push_back(e + 1);
                        sides[tag].push_back(s + 1);
                    }
                }
        elements[7] = {elements[0].front()};
        sides[7]    = {sides[0].front()};
    }
    void write(const std::string& path, bool multi) const
    {
        int file;
        nc(nc_create(path.c_str(), NC_CLOBBER, &file));
        auto dim = [&](const std::string& name, size_t size)
        {
            int id;
            nc(nc_def_dim(file, name.c_str(), size, &id));
            return id;
        };
        auto var = [&](const std::string& name, nc_type type, std::initializer_list<int> dims)
        {
            int id;
            nc(nc_def_var(file, name.c_str(), type, int(dims.size()), dims.begin(), &id));
            return id;
        };
        int nn = dim("num_nodes", 125);
        dim("num_elem", cells.size());
        int nb     = dim("num_el_blk", multi ? 3 : 1);
        int status = var("eb_status", NC_INT, {nb});
        std::vector<int> conn;
        for (int block = 0; block < (multi ? 3 : 1); ++block)
        {
            if (block == 1) continue;
            std::string suffix = std::to_string(block + 1);
            int ne             = dim("num_el_in_blk" + suffix, multi ? cells.size() / 2 : cells.size()),
                ncorners       = dim("num_nod_per_el" + suffix, N);
            conn.push_back(var("connect" + suffix, NC_INT, {ne, ncorners}));
        }
        int coords[3];
        const char* axes[] = {"coordx", "coordy", "coordz"};
        for (int a = 0; a < 3; ++a)
            coords[a] = var(axes[a], NC_DOUBLE, {nn});
        int ns = dim("num_side_sets", 8), len = dim("len_name", 16);
        int sn = var("ss_names", NC_CHAR, {ns, len}), si = var("ss_prop1", NC_INT, {ns}),
            st = var("ss_status", NC_INT, {ns});
        int elem[8], side[8];
        for (int t = 0; t < 8; ++t)
        {
            if (t == 6) continue;
            std::string suffix = std::to_string(t + 1);
            int n              = dim("num_side_ss" + suffix, elements[t].size());
            elem[t]            = var("elem_ss" + suffix, NC_INT, {n});
            side[t]            = var("side_ss" + suffix, NC_INT, {n});
        }
        nc(nc_enddef(file));
        const int active[] = {1, 0, 1};
        nc(nc_put_var_int(file, status, active));
        for (size_t b = 0; b < conn.size(); ++b)
        {
            std::vector<int> values;
            for (size_t e = b * (cells.size() / 2); e < (multi ? (b + 1) * (cells.size() / 2) : cells.size()); ++e)
                for (int n : cells[e])
                    values.push_back(n + 1);
            nc(nc_put_var_int(file, conn[b], values.data()));
        }
        for (int a = 0; a < 3; ++a)
        {
            std::array<double, 125> values;
            for (size_t i = 0; i < points.size(); ++i)
                values[i] = points[i][a];
            nc(nc_put_var_double(file, coords[a], values.data()));
        }
        std::array<char, 128> labels{};
        int ids[8], ss_active[8];
        for (int t = 0; t < 8; ++t)
        {
            std::copy(names[t].begin(), names[t].end(), labels.begin() + t * 16);
            ids[t]       = 100 + t * 13;
            ss_active[t] = t != 6;
        }
        nc(nc_put_var_text(file, sn, labels.data()));
        nc(nc_put_var_int(file, si, ids));
        nc(nc_put_var_int(file, st, ss_active));
        for (int t = 0; t < 8; ++t)
        {
            if (t == 6) continue;
            nc(nc_put_var_int(file, elem[t], elements[t].data()));
            nc(nc_put_var_int(file, side[t], sides[t].data()));
        }
        nc(nc_close(file));
    }
};

template<int N>
void check(const std::string& path, const Fixture<N>& fixture, bool multi, bool block_identity, int rank, int ranks)
{
    using Key           = uint64_t;
    using Tag           = std::conditional_t<N == 8, mars::HexTag, mars::TetTag>;
    using Domain        = mars::ElementDomain<Tag, double, Key, cstone::execution::Gpu>;
    const bool separate = multi && block_identity;
    if (!rank)
        std::cout << "DOMAIN_CASE corners=" << N << " multiblock=" << multi << " separate=" << separate << std::endl;
    // Hex8 block identities use the SFC-coordinate path; stored coordinates cannot co-move block tags yet.
    Domain domain(path, rank, ranks, !(N == 8 && separate), 64, 8, 0, 0., 0., true);
    const auto& d_keys = domain.getLocalToGlobalSfcMap();
    std::vector<Key> keys(d_keys.size());
    thrust::copy(thrust::device_pointer_cast(d_keys.data()), thrust::device_pointer_cast(d_keys.data() + d_keys.size()),
                 keys.begin());
    require(domain.numBlocks() == (separate ? 3 : 1), "unexpected node identity mode");
    std::vector<int> blocks(keys.size(), 0);
    if (separate)
    {
        const auto& d = domain.getLocalToGlobalBlock();
        require(d.size() == keys.size(), "bad block map");
        thrust::copy(thrust::device_pointer_cast(d.data()), thrust::device_pointer_cast(d.data() + d.size()),
                     blocks.begin());
    }
    if (!rank) require(domain.hasSideSet("single_face"), "rank-local getter lost a boundary");
    MPI_Barrier(MPI_COMM_WORLD);
    require(domain.sideSetNames().size() == 8, "domain lost side-set names");
    for (int tag = 0; tag < 8; ++tag)
    {
        std::set<std::pair<Key, int>> global;
        if (tag != 6)
            for (size_t f = 0; f < fixture.elements[tag].size(); ++f)
            {
                int e = fixture.elements[tag][f] - 1, side = fixture.sides[tag][f] - 1;
                for (int c : Fixture<N>::face_nodes[side])
                {
                    auto p  = fixture.points[fixture.cells[e][c]];
                    Key key = cstone::sfc3D<cstone::HilbertKey<Key>>(p[0], p[1], p[2], domain.getBoundingBox()).value();
                    global.emplace(key, separate && size_t(e) >= fixture.cells.size() / 2 ? 2 : 0);
                }
            }
        std::set<int> expected;
        for (size_t i = 0; i < keys.size(); ++i)
            if (global.count({keys[i], blocks[i]})) expected.insert(int(i));
        const auto& d = domain.sideSetNodes(names[tag]);
        std::vector<int> actual(d.size());
        if (!actual.empty())
            thrust::copy(thrust::device_pointer_cast(d.data()), thrust::device_pointer_cast(d.data() + d.size()),
                         actual.begin());
        const auto& dk = domain.sideSetNodeKeys(names[tag]);
        std::vector<Key> actual_keys(dk.size());
        if (!actual_keys.empty())
            thrust::copy(thrust::device_pointer_cast(dk.data()), thrust::device_pointer_cast(dk.data() + dk.size()),
                         actual_keys.begin());
        require(std::set<int>(actual.begin(), actual.end()) == expected && actual.size() == expected.size(),
                "domain boundary nodes differ from geometric oracle");
        require(actual_keys.size() == actual.size(), "boundary arrays differ");
        for (size_t i = 0; i < actual.size(); ++i)
            require(actual_keys[i] == keys.at(actual[i]), "boundary key/local alignment changed");
        if constexpr (N == 4)
        {
            std::vector<Key> expected_faces;
            for (size_t f = 0; f < fixture.elements[tag].size(); ++f)
            {
                int e = fixture.elements[tag][f] - 1, side = fixture.sides[tag][f] - 1;
                int block = separate && size_t(e) >= fixture.cells.size() / 2 ? 2 : 0;
                std::array<Key, 3> face;
                bool complete = true;
                for (int c = 0; c < 3; ++c)
                {
                    auto p  = fixture.points[fixture.cells[e][Fixture<N>::face_nodes[side][c]]];
                    face[c] = cstone::sfc3D<cstone::HilbertKey<Key>>(p[0], p[1], p[2], domain.getBoundingBox()).value();
                    bool found = false;
                    for (size_t i = 0; i < keys.size(); ++i)
                        found |= keys[i] == face[c] && blocks[i] == block;
                    complete &= found;
                }
                if (complete) expected_faces.insert(expected_faces.end(), face.begin(), face.end());
            }
            require(domain.has_side_set_face_keys(names[tag]), "Tet4 boundary face ancestry missing");
            const auto& df = domain.side_set_face_keys(names[tag]);
            std::vector<Key> actual_faces(df.size());
            if (!actual_faces.empty())
                thrust::copy(thrust::device_pointer_cast(df.data()), thrust::device_pointer_cast(df.data() + df.size()),
                             actual_faces.begin());
            require(actual_faces == expected_faces, "Tet4 boundary face winding/order/coverage changed");
        }
        int present = int(!expected.empty()), total = 0;
        MPI_Allreduce(&present, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        require(tag == 6 ? total == 0 : total > 0, "global boundary coverage lost");
    }
    cudaCheckCall(cudaDeviceSynchronize());
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try
    {
        require(argc == 3, "expected public fixture directory and conformal|blocks mode");
        const std::string mode = argv[2];
        require(mode == "conformal" || mode == "blocks", "unknown node identity mode");
        const bool block_identity = mode == "blocks";
        // The domain caches this option on its first use, so modes run in separate processes.
        require((block_identity ? setenv("MARS_BLOCK_NODE_IDENTITY", "1", 1) : unsetenv("MARS_BLOCK_NODE_IDENTITY")) ==
                    0,
                "cannot select node identity mode");
        MPI_Comm local_comm;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &local_comm);
        int local_rank, devices;
        MPI_Comm_rank(local_comm, &local_rank);
        cudaCheckCall(cudaGetDeviceCount(&devices));
        require(local_rank < devices, "one GPU per rank required");
        cudaCheckCall(cudaSetDevice(local_rank));
        MPI_Comm_free(&local_comm);
        std::filesystem::create_directories(argv[1]);
        auto run = [&]<int N>()
        {
            Fixture<N> fixture;
            for (bool multi : {false, true})
            {
                std::string path =
                    std::string(argv[1]) + "/n" + std::to_string(N) + (multi ? "-blocks.e" : "-single.e");
                if (!rank) fixture.write(path, multi);
                MPI_Barrier(MPI_COMM_WORLD);
                check(path, fixture, multi, block_identity, rank, ranks);
            }
        };
        run.template operator()<8>();
        run.template operator()<4>();
        if (!rank)
            std::cout
                << "PASS: Hex8/Tet4 ElementDomain boundary labels after redistribution, including ghosts and block "
                   "identity mode="
                << mode << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "rank " << rank << ": " << e.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
