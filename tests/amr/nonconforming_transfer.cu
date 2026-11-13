#include "tet_gpu_checks.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_manager.hpp"
#include "backend/distributed/unstructured/fem/mars_dof_space.hpp"
#include "graded_mesh.hpp"

namespace amr_test
{
void require(bool valid, const char* message)
{
    if (!valid) throw std::runtime_error(message);
}
}

namespace
{
using namespace amr_test;
using Fields = std::array<Values, 5>;
struct NodeField
{
    KeyType key;
    std::array<double, 5> values;
    int owner;
};

template<class DomainType>
std::vector<NodeField> owned_fields(DomainType& domain, const Fields& fields, int rank)
{
    auto keys = read_device(domain.getLocalToGlobalSfcMap().data(), domain.getNodeCount());
    auto owned = read_device(domain.getNodeOwnershipMap().data(), keys.size());
    std::array<std::vector<double>, 5> values;
    for (int f = 0; f < 5; ++f) values[f] = read_device(fields[f].data(), keys.size());
    std::vector<NodeField> local;
    for (size_t i = 0; i < keys.size(); ++i)
        if (owned[i])
            local.push_back({keys[i], {values[0][i], values[1][i], values[2][i], values[3][i], values[4][i]}, rank});
    auto all = gather(local);
    std::sort(all.begin(), all.end(), [](const auto& a, const auto& b) { return a.key < b.key; });
    for (size_t i = 1; i < all.size(); ++i) require(all[i - 1].key != all[i].key, "Duplicate field owner");
    return all;
}

template<class DomainType>
std::unique_ptr<DomainType> empty_input_cube(int rank, int ranks)
{
    constexpr int N = DomainType::NodesPerElement;
    std::array<std::vector<KeyType>, N> columns;
    std::array<std::vector<double>, 3> coords;
    if (rank == 0)
    {
        for (int z = 0; z < 5; ++z)
            for (int y = 0; y < 5; ++y)
                for (int x = 0; x < 5; ++x)
                {
                    coords[0].push_back(x * .25); coords[1].push_back(y * .25); coords[2].push_back(z * .25);
                }
        auto id = [](int x, int y, int z) { return KeyType(x + 5 * (y + 5 * z)); };
        constexpr int tet[6][4] = {{0,1,2,6}, {0,2,3,6}, {0,3,7,6}, {0,7,4,6}, {0,4,5,6}, {0,5,1,6}};
        for (int z = 0; z < 4; ++z)
            for (int y = 0; y < 4; ++y)
                for (int x = 0; x < 4; ++x)
                {
                    KeyType h[8] = {id(x,y,z), id(x+1,y,z), id(x+1,y+1,z), id(x,y+1,z),
                                    id(x,y,z+1), id(x+1,y,z+1), id(x+1,y+1,z+1), id(x,y+1,z+1)};
                    if constexpr (N == 8)
                        for (int c = 0; c < N; ++c) columns[c].push_back(h[c]);
                    else
                        for (const auto& t : tet)
                            for (int c = 0; c < N; ++c) columns[c].push_back(h[t[c]]);
                }
    }
    typename DomainType::DeviceCoordsTuple xyz{upload(coords[0]), upload(coords[1]), upload(coords[2])};
    typename DomainType::DeviceConnectivityTuple conn;
    [&]<size_t... I>(std::index_sequence<I...>) { ((std::get<I>(conn) = upload(columns[I])), ...); }
    (std::make_index_sequence<N>{});
    std::printf("[amr-transfer-input] rank=%d elements=%zu nodes=%zu\n", rank, columns[0].size(), coords[0].size());
    std::fflush(stdout);
    return std::make_unique<DomainType>(std::move(xyz), std::move(conn), rank, ranks, 8);
}

template<class DomainType>
std::unique_ptr<DomainType> general_input(const std::string& fixture, int rank, int ranks)
{
    constexpr int N = DomainType::NodesPerElement;
    auto mesh       = graded_mesh<N>();
    if constexpr (N == 4)
        if (fixture == "thin") mesh = thin_mesh();
    typename DomainType::DeviceCoordsTuple xyz{upload(mesh.xyz[0]), upload(mesh.xyz[1]), upload(mesh.xyz[2])};
    typename DomainType::DeviceConnectivityTuple conn;
    std::array<std::vector<KeyType>, N> local;
    for (size_t e = 0; e < mesh.conn[0].size(); ++e)
        if (int(e % ranks) == rank)
            for (int c = 0; c < N; ++c)
                local[c].push_back(mesh.conn[c][e]);
    [&]<size_t... I>(std::index_sequence<I...>)
    { ((std::get<I>(conn) = upload(local[I])), ...); }(std::make_index_sequence<N>{});
    return std::make_unique<DomainType>(std::move(xyz), std::move(conn), rank, ranks, 8);
}

template<class DomainType>
void check_constraints(DomainType& domain, double tolerance, bool initial, bool staged)
{
    mars::fem::HangingNodes<KeyType, double> hanging(domain, staged);
    if (initial) require(hanging.globalCount() == 0, "Conforming input has false hanging constraints");
    mars::fem::DofSpace<KeyType, double, DomainType> space(domain, nullptr, &hanging);
    const size_t n = domain.getNodeCount();
    auto keys      = read_device(domain.getLocalToGlobalSfcMap().data(), n);
    auto mask = read_device(space.isDof(), n), current = read_device(space.isSynced(), n);
    auto x = read_device(domain.getNodeX().data(), n), y = read_device(domain.getNodeY().data(), n),
         z = read_device(domain.getNodeZ().data(), n);
    std::vector<double> u(n, 0), v(n, 0), affine(n, 0), constant(n, 0);
    for (size_t i = 0; i < n; ++i)
    {
        if (mask[i])
        {
            u[i]        = double(keys[i] % 1237) / 1237.;
            affine[i]   = 1 + 2 * x[i] - 3 * y[i] + 4 * z[i];
            constant[i] = 1;
        }
        if (current[i]) v[i] = double((keys[i] + 17) % 997) / 997.;
    }
    auto du = upload(u), dv = upload(v), da = upload(affine), dc = upload(constant);
    space.prolong({&du, &da, &dc});
    space.restrict(dv);
    auto pu = read_device(du.data(), n), ptv = read_device(dv.data(), n);
    affine     = read_device(da.data(), n);
    constant   = read_device(dc.data(), n);
    double lhs = 0, rhs = 0, error = 0, constant_error = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (current[i])
        {
            lhs += pu[i] * v[i];
            error          = std::max(error, std::abs(affine[i] - (1 + 2 * x[i] - 3 * y[i] + 4 * z[i])));
            constant_error = std::max(constant_error, std::abs(constant[i] - 1));
        }
        if (mask[i]) rhs += u[i] * ptv[i];
    }
    double totals[2] = {lhs, rhs};
    MPI_Allreduce(MPI_IN_PLACE, totals, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &constant_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    require(std::abs(totals[0] - totals[1]) <= 1.e-11 * std::max(1., std::abs(totals[0])),
            "Constraint P/P^T adjoint test failed");
    require(error <= 10 * tolerance && constant_error < 1.e-13, "Constraint affine or constant reproduction failed");
    if (domain.rank() == 0)
        std::printf(
            "[amr-constraints] hanging=%lld affine=%.17g affine_bound=%.17g adjoint=%.17g adjoint_bound=%.17g PASS\n",
            hanging.globalCount(), error, 10 * tolerance, std::abs(totals[0] - totals[1]),
            1.e-11 * std::max(1., std::abs(totals[0])));
}

template<class ElementTag>
void check_manager(const std::string& mesh,
                   bool empty_input,
                   bool octree,
                   bool staged,
                   int rank,
                   int ranks,
                   const std::string& fixture)
{
    using Manager = mars::amr::AmrManager<ElementTag, KeyType, double>;
    typename Manager::Config config;
    config.strategy = octree ? mars::amr::MarkingStrategy::OctreeNative : mars::amr::MarkingStrategy::Doerfler;
    config.maxLevels       = fixture.empty() ? 64 : 3;
    config.bucketSize = 8;
    config.coarsenFraction = 0;
    config.refineFraction = .5;
    config.stage_refinement_messages = staged;
    Manager manager(config);
    if (!fixture.empty())
        manager.initialize(general_input<typename Manager::Domain>(fixture, rank, ranks), rank, ranks);
    else if (empty_input)
        manager.initialize(empty_input_cube<typename Manager::Domain>(rank, ranks), rank, ranks);
    else manager.initialize(mesh, rank, ranks);
    if (rank == 0) { std::puts("[amr-transfer-domain] ready"); std::fflush(stdout); }
    if (!fixture.empty())
    {
        auto& domain      = manager.domain();
        const auto before = domain.localElementCount();
        Values zero(domain.getElementCount(), 0.);
        auto stats = manager.adaptMeshMultiField(zero.data(), {}, {});
        require(stats.elementsRefined == 0 && before == manager.domain().localElementCount(),
                "Zero error refined a graded mesh");
        check_constraints(manager.domain(), 1.e-6, true, staged);
    }
    Fields fields;
    auto& first = manager.domain();
    auto x = read_device(first.getNodeX().data(), first.getNodeCount());
    auto y = read_device(first.getNodeY().data(), x.size());
    auto z = read_device(first.getNodeZ().data(), x.size());
    auto keys = read_device(first.getLocalToGlobalSfcMap().data(), x.size());
    std::array<std::vector<double>, 5> values{std::vector<double>(x.size(), 2.5), x, y, z,
                                            std::vector<double>(x.size())};
    for (size_t i = 0; i < keys.size(); ++i) values[4][i] = double(keys[i] % 997) / 1024;
    for (int f = 0; f < 5; ++f) fields[f] = upload(values[f]);
    double affine_bound = 1.e-12;
    size_t moved_total = 0;
    for (int round = 0; round < 3; ++round)
    {
        auto& domain = manager.domain();
        const auto before = owned_fields(domain, fields, rank);
        const size_t nodes = domain.getNodeCount(), elements = domain.getElementCount();
        const auto* owned = domain.getNodeOwnershipMap().data();
        for (auto& field : fields)
        {
            auto* data = field.data();
            thrust::for_each_n(thrust::device, thrust::counting_iterator<size_t>(0), nodes,
                              [=] __device__(size_t i) { if (!owned[i]) data[i] = -1.e8; });
        }
        Values errors(elements, 0.);
        const auto& conn = domain.getElementToNodeConnectivity();
        const auto* corner = std::get<0>(conn).data();
        const auto* coord = round == 0 ? domain.getNodeX().data() : round == 1 ? domain.getNodeY().data()
                                                                                 : domain.getNodeZ().data();
        auto* error = errors.data();
        size_t begin = domain.startIndex();
        const bool continuous = !fixture.empty();
        thrust::for_each_n(thrust::device, thrust::counting_iterator<size_t>(begin), domain.localElementCount(),
                           [=] __device__(size_t e)
                           {
                               double distance = (coord[corner[e]] - .35) / .4;
                               error[e] =
                                   continuous ? std::exp(-distance * distance) : (coord[corner[e]] < .4 ? 1. : 0.);
                           });
        Fields next;
        auto stats = manager.adaptMeshMultiField(errors.data(),
            {fields[0].data(), fields[1].data(), fields[2].data(), fields[3].data(), fields[4].data()},
            {&next[0], &next[1], &next[2], &next[3], &next[4]});
        require(stats.elementsRefined > 0, "Transfer fixture did not refine");
        auto& updated = manager.domain();
        const auto after = owned_fields(updated, next, rank);
        size_t moved = 0;
        for (const auto& node : before)
        {
            auto found = std::lower_bound(after.begin(), after.end(), node.key,
                [](const NodeField& a, KeyType key) { return a.key < key; });
            require(found != after.end() && found->key == node.key && found->values == node.values,
                    "Adaptation changed an existing node field");
            moved += found->owner != node.owner;
        }
        moved_total += moved;
        auto current = read_device(updated.getNodeOwnershipMap().data(), updated.getNodeCount());
        const auto& receive = updated.getNodeHaloTopology().recvNodeIds_;
        for (int i : read_device(receive.data(), receive.size())) current.at(i) = 1;
        x = read_device(updated.getNodeX().data(), current.size());
        y = read_device(updated.getNodeY().data(), current.size());
        z = read_device(updated.getNodeZ().data(), current.size());
        const auto& box = updated.getBoundingBox();
        auto bits = box.getBoxDimBits(cstone::maxTreeLevel<KeyType>{});
        double cell = std::max({box.lx()/double((uint64_t(1)<<bits[0])-1),
                               box.ly()/double((uint64_t(1)<<bits[1])-1),
                               box.lz()/double((uint64_t(1)<<bits[2])-1)});
        // Convex interpolation cannot amplify the old error; snapping adds one key-cell interval.
        affine_bound += cell;
        double constant_error = 0, affine_error = 0;
        for (int f = 0; f < 4; ++f)
        {
            auto actual = read_device(next[f].data(), current.size());
            for (size_t i = 0; i < current.size(); ++i)
                if (current[i])
                {
                    require(std::isfinite(actual[i]), "Nonfinite transferred field");
                    double expected = f == 0 ? 2.5 : f == 1 ? x[i] : f == 2 ? y[i] : z[i];
                    double& maximum = f == 0 ? constant_error : affine_error;
                    maximum = std::max(maximum, std::abs(actual[i] - expected));
                }
        }
        MPI_Allreduce(MPI_IN_PLACE, &constant_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &affine_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        require(constant_error <= 1.e-13 && affine_error <= affine_bound, "Immediate post-adapt field check failed");
        if (!fixture.empty())
        {
            check_constraints(updated, affine_bound, false, staged);
            auto generations = read_device(updated.refinement_generations().data(), updated.getElementCount());
            require(
                std::all_of(generations.begin(), generations.end(), [&](auto g) { return g <= KeyType(round + 1); }),
                "Generation did not survive rebuild and halo exchange");
        }
        fields = std::move(next);
        if (rank == 0)
            std::printf("[amr-transfer] round=%d old_nodes=%zu moved=%zu constant=%.17g affine=%.17g bound=%.17g PASS\n",
                        round, before.size(), moved, constant_error, affine_error, affine_bound);
    }
    if (!fixture.empty())
    {
        bool stopped = false;
        for (int round = 0; round < 4; ++round)
        {
            auto& domain = manager.domain();
            Values uniform(domain.getElementCount(), 1.);
            auto stats          = manager.adaptMeshMultiField(uniform.data(), {}, {});
            const auto& current = manager.domain();
            auto generations    = read_device(current.refinement_generations().data(), current.getElementCount());
            require(std::all_of(generations.begin(), generations.end(), [](auto g) { return g <= 3; }),
                    "Generation cap exceeded");
            if (stats.elementsRefined == 0)
            {
                stopped = true;
                break;
            }
        }
        require(stopped, "Uniform requests did not stop at the generation limit");
        const auto& capped = manager.domain();
        auto generations   = read_device(capped.refinement_generations().data(), capped.getElementCount());
        require(std::all_of(generations.begin(), generations.end(), [](auto g) { return g == 3; }),
                "Uniform refinement stopped below the generation cap");
        uint64_t cells = capped.localElementCount();
        MPI_Allreduce(MPI_IN_PLACE, &cells, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
        const uint64_t roots = fixture == "thin" ? 448 : ElementTag::NodesPerElement == 4 ? 384 : 64;
        require(cells == roots * 512, "Uniform generation cap has the wrong element count");
        if (rank == 0) std::printf("[amr-generation] maximum=3 elements=%llu PASS\n", (unsigned long long)cells);
        check_constraints(manager.domain(), 1.e-5, false, staged);
    }
    require(ranks == 1 || empty_input || !fixture.empty() || moved_total > 0,
            "Transfer fixture did not change node owners");
    if (rank == 0) std::puts("[amr-transfer] PASS");
}
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try
    {
        MPI_Comm local;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &local);
        int local_rank, devices;
        MPI_Comm_rank(local, &local_rank); cudaCheckCall(cudaGetDeviceCount(&devices));
        amr_test::require(devices == 1 || local_rank < devices, "No device for local rank");
        int device = devices == 1 ? 0 : local_rank;
        cudaCheckCall(cudaSetDevice(device)); MPI_Comm_free(&local);
        char host[MPI_MAX_PROCESSOR_NAME], pci[32]; int length;
        MPI_Get_processor_name(host, &length);
#if defined(__HIPCC__)
        cudaCheckCall(hipDeviceGetPCIBusId(pci, sizeof(pci), device));
#else
        cudaCheckCall(cudaDeviceGetPCIBusId(pci, sizeof(pci), device));
#endif
        std::printf("[amr-transfer-device] rank=%d host=%s pci=%s\n", rank, host, pci); std::fflush(stdout);
        std::string mesh, fixture, element = "tet";
        bool empty = false, octree = false, staged = false;
        for (int i = 1; i < argc; ++i)
        {
            std::string a = argv[i];
            if (a.rfind("--mesh=", 0) == 0) mesh = a.substr(7);
            else if (a.rfind("--element=", 0) == 0) element = a.substr(10);
            else if (a.rfind("--fixture=", 0) == 0)
                fixture = a.substr(10);
            else if (a == "--empty-input") empty = true;
            else if (a == "--octree") octree = true;
            else if (a == "--staged") staged = true;
            else throw std::runtime_error("Unknown transfer-test argument");
        }
        amr_test::require((empty || !mesh.empty() || fixture == "graded" || (fixture == "thin" && element == "tet")) &&
                              (element == "tet" || element == "hex"),
                          "Invalid transfer-test arguments");
        if (element == "tet")
            check_manager<mars::TetTag>(mesh, empty, octree, staged, rank, ranks, fixture);
        else
            check_manager<mars::HexTag>(mesh, empty, octree, staged, rank, ranks, fixture);
    }
    catch (const std::exception& e)
    {
        std::fprintf(stderr, "[amr-transfer] FAIL: %s\n", e.what()); MPI_Abort(MPI_COMM_WORLD, 82);
    }
    MPI_Finalize();
}
