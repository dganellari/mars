#include "tet_checks.hpp"
#include "tet_exodus_fixture.hpp"
#include "backend/distributed/unstructured/utils/mars_read_exodus_mesh.hpp"

#ifndef MARS_HAVE_NETCDF
#error "The Exodus gate requires the actual netCDF reader"
#endif

namespace amr_test
{
void check_import(const char* path, int rank, int ranks)
{
    auto [nodes, count, x, y, z, connectivity, file_ids, boundary] =
        mars::readExodusMeshWithElementPartitioning<4, RealType, KeyType>(path, rank, ranks);
    const auto fixture = exodus_fixture();
    require(nodes == file_ids.size() && count > 0, "Invalid imported partition size");
    std::vector<Tet> local(count);
    const std::array<const KeyType*, 4> conn = {std::get<0>(connectivity).data(), std::get<1>(connectivity).data(),
                                                std::get<2>(connectivity).data(), std::get<3>(connectivity).data()};
    for (size_t e = 0; e < count; ++e)
        for (int c = 0; c < 4; ++c)
        {
            const auto i = conn[c][e];
            require(i < nodes && file_ids[i] < fixture.points.size(), "Imported connectivity is out of range");
            local[e].points[c] = {x[i], y[i], z[i]};
            require(local[e].points[c] == fixture.points[file_ids[i]],
                    "File node IDs no longer identify their coordinates");
        }
    const auto global = gather(local);
    require(global.size() == fixture.cells.size(), "Imported element count changed");
    check_mesh(global);
    std::set<std::array<Point, 4>> actual, expected;
    for (const auto& tet : global)
    {
        auto corners = tet.points;
        std::sort(corners.begin(), corners.end());
        actual.insert(corners);
    }
    for (const auto& cell : fixture.cells)
    {
        std::array<Point, 4> corners;
        for (int c = 0; c < 4; ++c)
            corners[c] = fixture.points[cell[c]];
        std::sort(corners.begin(), corners.end());
        expected.insert(corners);
    }
    require(actual == expected, "Imported connectivity differs from the fixture");
    auto [blocks, number] = mars::readExodusElementBlocks(path, rank, ranks);
    require(number == 1 && blocks.size() == count &&
                std::all_of(blocks.begin(), blocks.end(), [](auto block) { return block == 0; }),
            "Nonsequential Exodus block ID was mistaken for a block index");

    auto sets = mars::readExodusSideSetsTet4(path, rank);
    require(sets.trianglesByName.size() == 7, "Imported boundary names changed");
    for (int tag = 0; tag < 7; ++tag)
    {
        auto name          = "boundary" + std::to_string(tag);
        const auto& ids    = sets.trianglesByName.at(name);
        const auto& points = sets.triangleCoordsByName.at(name);
        require(ids.size() == points.size() && ids.size() % 3 == 0, "Misaligned Exodus face arrays");
        std::set<Face> faces, wanted;
        std::set<uint64_t> face_nodes;
        for (size_t i = 0; i < ids.size(); i += 3)
        {
            Face face;
            for (int c = 0; c < 3; ++c)
            {
                require(ids[i + c] < fixture.points.size() && points[i + c] == fixture.points[ids[i + c]],
                        "Imported triangle coordinates disagree with file node IDs");
                face[c] = points[i + c];
                face_nodes.insert(ids[i + c]);
            }
            std::sort(face.begin(), face.end());
            require(faces.insert(face).second, "Repeated face inside an imported side set");
        }
        for (const auto& tet : global)
            for (int opposite = 0; opposite < 4; ++opposite)
            {
                auto face = tet_face(tet, opposite);
                if (boundary_plane(face) == (tag == 6 ? 0 : tag)) wanted.insert(face);
            }
        require(faces == wanted && faces.size() == 32, "Imported side-set labels do not cover the expected plane");
        const auto& node_ids    = sets.nodesByName.at(name);
        const auto& node_points = sets.nodeCoordsByName.at(name);
        require(std::set<uint64_t>(node_ids.begin(), node_ids.end()) == face_nodes && node_ids.size() == 25 &&
                    node_ids.size() == node_points.size(),
                "Imported boundary node list changed");
        for (size_t i = 0; i < node_ids.size(); ++i)
            require(node_points[i] == fixture.points[node_ids[i]], "Imported boundary node coordinates changed");
    }
}
} // namespace amr_test

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    char host[MPI_MAX_PROCESSOR_NAME];
    int host_length = 0;
    MPI_Get_processor_name(host, &host_length);
    std::cout << "CPU rank=" << rank << " host=" << std::string(host, host_length) << std::endl;
    try
    {
        amr_test::require(argc == 2, "Supply the public Exodus fixture path");
        amr_test::check_import(argv[1], rank, ranks);
        if (rank == 0) std::cout << "Tet4 Exodus reader gate PASS on " << ranks << " ranks\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << "rank " << rank << ": " << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
}
