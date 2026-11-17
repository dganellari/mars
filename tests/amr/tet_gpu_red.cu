#include "tet_gpu_checks.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_manager.hpp"

namespace amr_test
{
using Manager = mars::amr::AmrManager<mars::TetTag, KeyType, RealType>;

void check_default_red(Domain& domain, int rank, int ranks)
{
    require(Manager::Config{}.tet_refinement == mars::amr::TetRefinement::Red, "Default tet policy changed");
    using Refiner = mars::amr::OctreeAlignedRefine<KeyType, RealType, mars::TetTag>;
    size_t start = domain.startIndex(), count = domain.localElementCount();
    const auto& conn = domain.getElementToNodeConnectivity();
    auto marks       = upload(std::vector<uint8_t>(count, 1));
    auto result      = Refiner::refineLocal(std::get<0>(conn).data() + start, std::get<1>(conn).data() + start,
                                            std::get<2>(conn).data() + start, std::get<3>(conn).data() + start,
                                            domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(),
                                            marks.data(), count, domain.getNodeCount());
    require(!result.transfer_from_parent_edges && result.numElements == 8 * count, "Default red refinement changed");
    auto fields = initial_fields(domain);
    // Unrefined ranks have no marks array even when other ranks split elements.
    for (int f = 0; f < 2; ++f)
    {
        auto unchanged = mars::amr::SolutionTransfer<KeyType, RealType>::transferByParentageTet(
            nullptr, nullptr, nullptr, nullptr, nullptr, count, fields[f].data(), domain.getNodeCount(),
            nullptr, 0, domain.getNodeCount());
        require(read_device(unchanged.data(), unchanged.size()) == read_device(fields[f].data(), fields[f].size()),
                "Unrefined-rank field copy changed values");
    }
    std::array<Values, 2> transferred;
    for (int f = 0; f < 2; ++f)
        Refiner::transferSolution(std::get<0>(conn).data() + start, std::get<1>(conn).data() + start,
                                  std::get<2>(conn).data() + start, std::get<3>(conn).data() + start, fields[f].data(),
                                  result, count, transferred[f]);
    check_fields(result.d_x.data(), result.d_y.data(), result.d_z.data(), result.numNodes, transferred, 1.e-12);
    auto children = snapshot(
        {result.d_conn[0].data(), result.d_conn[1].data(), result.d_conn[2].data(), result.d_conn[3].data()}, 0,
        result.numElements, result.d_x.data(), result.d_y.data(), result.d_z.data(), result.numNodes, nullptr, rank);
    check_mesh(gather(children));
    Domain::DeviceCoordsTuple coordinates{std::move(result.d_x), std::move(result.d_y), std::move(result.d_z)};
    Domain::DeviceConnectivityTuple connectivity{std::move(result.d_conn[0]), std::move(result.d_conn[1]),
                                                 std::move(result.d_conn[2]), std::move(result.d_conn[3])};
    Domain rebuilt(std::move(coordinates), std::move(connectivity), rank, ranks, 8);
    rebuilt.cacheNodeCoordinates();
    (void)rebuilt.getNodeOwnershipMap();
    check_mesh(gather(snapshot(rebuilt, rank)));
    if (rank == 0) std::cout << "Default red refinement, parentage and device rebuild PASS\n";
}

} // namespace amr_test
