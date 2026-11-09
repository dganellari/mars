#include "tet_reference.hpp"
#include "tet_gpu_checks.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_tet_bisection.hpp"

namespace amr_test
{
void check_coarsening_ignored(Domain& domain, int rank)
{
    auto marks   = upload(std::vector<uint8_t>(domain.getElementCount(), 2));
    auto refined = mars::amr::TetBisectionRefiner<KeyType, RealType>::refine(domain, marks.data());
    require(refined.numMarked == 0 && refined.numUniqueEdges == 0 && refined.d_edgeKeys.size() == 0,
            "Coarsening marks requested bisection edges");
    require(refined.numElements == domain.localElementCount() && refined.numNodes == domain.getNodeCount(),
            "Coarsening marks changed mesh counts");
    const auto& conn                    = domain.getElementToNodeConnectivity();
    std::array<const KeyType*, 4> input = {std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                                           std::get<3>(conn).data()};
    for (int c = 0; c < 4; ++c)
        require(read_device(input[c] + domain.startIndex(), domain.localElementCount()) ==
                    read_device(refined.d_conn[c].data(), refined.numElements),
                "Coarsening marks changed connectivity");
    require(read_device(domain.getNodeX().data(), domain.getNodeCount()) ==
                    read_device(refined.d_x.data(), refined.numNodes) &&
                read_device(domain.getNodeY().data(), domain.getNodeCount()) ==
                    read_device(refined.d_y.data(), refined.numNodes) &&
                read_device(domain.getNodeZ().data(), domain.getNodeCount()) ==
                    read_device(refined.d_z.data(), refined.numNodes),
            "Coarsening marks changed coordinates");
    if (rank == 0) std::cout << "Coarsening marks leave bisection mesh unchanged PASS\n";
}

void check_single_request(Domain& domain, int rank, int ranks, bool staged)
{
    auto local      = snapshot(domain, rank);
    auto global     = gather(local);
    RealType volume = check_mesh(global);
    std::map<Edge, std::set<int>> edge_ranks;
    std::map<Edge, size_t> edge_stars;
    for (const auto& tet : global)
        for (int a = 0; a < 4; ++a)
            for (int b = a + 1; b < 4; ++b)
            {
                Edge edge = {tet.keys[a], tet.keys[b]};
                std::sort(edge.begin(), edge.end());
                edge_ranks[edge].insert(tet.rank);
                ++edge_stars[edge];
            }
    auto chosen = global.end();
    for (auto it = global.begin(); it != global.end(); ++it)
        if (ranks == 1 || edge_ranks[longest_edge(*it)].size() > 1)
        {
            chosen = it;
            break;
        }
    require(chosen != global.end(), "Fixture has no longest-edge star crossing an MPI partition");
    auto identity = chosen->keys;
    std::sort(identity.begin(), identity.end());
    std::vector<uint8_t> marks(domain.getElementCount(), 0);
    for (size_t e = 0; e < local.size(); ++e)
    {
        auto keys = local[e].keys;
        std::sort(keys.begin(), keys.end());
        if (keys == identity) marks[domain.startIndex() + e] = 1;
    }
    auto device_marks = upload(marks);
    auto fields       = initial_fields(domain);
    auto refined      = mars::amr::TetBisectionRefiner<KeyType, RealType>::refine(domain, device_marks.data(), staged);
    auto output       = snapshot(
        {refined.d_conn[0].data(), refined.d_conn[1].data(), refined.d_conn[2].data(), refined.d_conn[3].data()}, 0,
        refined.numElements, refined.d_x.data(), refined.d_y.data(), refined.d_z.data(), refined.numNodes, nullptr,
        rank);
    auto all_output = gather(output);
    require(all_output.size() > global.size(), "LEPP did not refine the selected tetrahedron");
    // Run the same shared core on one CPU rank; canonical connectivity must agree.
    std::map<KeyType, Point> unique;
    for (const auto& tet : global)
        for (int c = 0; c < 4; ++c)
            unique.emplace(tet.keys[c], tet.points[c]);
    std::map<KeyType, KeyType> local_ids;
    std::vector<KeyType> keys;
    std::vector<RealType> x, y, z;
    for (const auto& [key, point] : unique)
    {
        local_ids[key] = keys.size();
        keys.push_back(key);
        x.push_back(point[0]);
        y.push_back(point[1]);
        z.push_back(point[2]);
    }
    std::array<std::vector<KeyType>, 4> conn;
    std::vector<uint8_t> reference_marks;
    for (const auto& tet : global)
    {
        for (int c = 0; c < 4; ++c)
            conn[c].push_back(local_ids.at(tet.keys[c]));
        auto sorted = tet.keys;
        std::sort(sorted.begin(), sorted.end());
        reference_marks.push_back(sorted == identity);
    }
    mars::amr::TetLeppView<KeyType, RealType> view{{conn[0].data(), conn[1].data(), conn[2].data(), conn[3].data()},
                                                   keys.data(),
                                                   x.data(),
                                                   y.data(),
                                                   z.data(),
                                                   reference_marks.data(),
                                                   global.size(),
                                                   keys.size()};
    auto reference = refine_reference(view, domain.getBoundingBox());
    std::vector<std::array<Point, 4>> expected, actual;
    for (size_t e = 0; e < reference.conn[0].size(); ++e)
    {
        std::array<Point, 4> points;
        for (int c = 0; c < 4; ++c)
        {
            auto i    = reference.conn[c][e];
            points[c] = {reference.x[i], reference.y[i], reference.z[i]};
        }
        std::sort(points.begin(), points.end());
        expected.push_back(points);
    }
    for (auto tet : all_output)
    {
        std::sort(tet.points.begin(), tet.points.end());
        actual.push_back(tet.points);
    }
    std::sort(expected.begin(), expected.end());
    std::sort(actual.begin(), actual.end());
    require(actual == expected, "Distributed device LEPP differs from the single-rank CPU core");
    require(std::abs(check_mesh(all_output) - volume) < 1.e-12, "Bisection did not preserve volume");
    std::array<Values, 2> transferred;
    for (int f = 0; f < 2; ++f)
        transferred[f] = mars::amr::transfer_lepp_field<mars::amr::LeppDeviceExecution>(
            fields[f].data(), refined.oldNumNodes, refined.d_edgeKeys.data(), refined.numUniqueEdges,
            refined.parent_birth_offsets, MPI_COMM_WORLD);
    check_fields(refined.d_x.data(), refined.d_y.data(), refined.d_z.data(), refined.numNodes, transferred, 1.e-12);
    if (rank == 0)
        std::cout << "LEPP closure with CPU parity: initial star=" << edge_stars[longest_edge(*chosen)]
                  << ", participating ranks=" << edge_ranks[longest_edge(*chosen)].size() << " PASS\n";
}

} // namespace amr_test
