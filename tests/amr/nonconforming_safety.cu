#include "tet_gpu_checks.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_octree.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_solution_transfer.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_octree_refine.hpp"
#include "field_migration_checks.hpp"
#include "nonconforming_topology_checks.hpp"

namespace
{
void require_all(bool valid, const char* message)
{
    int bad = !valid;
    MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (bad)
    {
        std::fprintf(stderr, "[amr-safety] FAIL: %s\n", message);
        MPI_Abort(MPI_COMM_WORLD, 80);
    }
}

template<class KeyType, class RealType>
void check_node_keys()
{
    using namespace amr_test;
    const cstone::Box<RealType> box(0, 1);
    const std::vector<RealType> x{0, RealType(0.25), RealType(0.5), RealType(0.75)},
        y{RealType(0.75), RealType(0.5), RealType(0.25), 0}, z(4, RealType(0.125));
    auto d_x = upload(x), d_y = upload(y), d_z = upload(z);
    std::vector<KeyType> expected;
    for (size_t i = 0; i < x.size(); ++i)
        expected.push_back(cstone::sfc3D<cstone::SfcKind<KeyType>>(x[i], y[i], z[i], box));
    for (bool mixed : {false, true})
    {
        auto poison = std::vector<KeyType>(x.size(), cstone::removeKey<KeyType>::value);
        if (mixed) poison[1] = 0;
        auto keys = upload(poison);
        mars::generateSfcKeys<KeyType, RealType>(d_x.data(), d_y.data(), d_z.data(), keys.data(), keys.size(), box);
        require_all(read_device(keys.data(), keys.size()) == expected,
                    "Fresh node keys retained a removal marker");
    }
    mars::generateSfcKeys<KeyType, RealType>(nullptr, nullptr, nullptr, nullptr, 0, box);
}

template<class Real>
void check_marking(int rank, int ranks)
{
    using namespace amr_test;
    std::vector<KeyType> leaves(9), keys(8);
    std::vector<Real> errors{0, Real(.1), Real(.2999), Real(.3), Real(.3001), Real(.5), Real(.9), 1};
    for (int level : {1, int(cstone::maxTreeLevel<KeyType>{})})
    {
        for (int i = 0; i <= 8; ++i)
            leaves[i] = KeyType(i) * cstone::nodeRange<KeyType>(level);
        for (int i = 0; i < 8; ++i)
            keys[i] = leaves[i];
        size_t first = size_t(8 * rank / ranks), last = size_t(8 * (rank + 1) / ranks);
        Real maximum = 0;
        for (size_t i = first; i < last; ++i)
            maximum = std::max(maximum, errors[i]);
        MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, std::is_same_v<Real, double> ? MPI_DOUBLE : MPI_FLOAT, MPI_MAX,
                      MPI_COMM_WORLD);
        auto dl = upload(leaves), dk = upload(keys);
        auto de = upload(errors);
        for (unsigned bucket : {1u, 64u, 1024u})
            for (Real fraction : {Real(0), Real(.3), Real(.5), Real(1)})
                for (int threads : {1, 64, 256})
                {
                    cstone::DeviceVector<unsigned> counts(8);
                    cstone::DeviceVector<cstone::TreeNodeIndex> ops(9);
                    mars::amr::weight_owned_leaves<<<8, threads>>>(dl.data(), 8, dk.data(), first, last, de.data(),
                                                                   maximum, bucket, fraction, counts.data());
                    cudaCheckError();
                    cstone::computeNodeOpsGpu(cstone::execution::gpuDefaultStream, dl.data(), 8, counts.data(), bucket,
                                              ops.data());
                    cstone::DeviceVector<uint8_t> marks(8, 0);
                    mars::amr::mark_owned_leaves<<<8, threads>>>(dl.data(), 8, dk.data(), first, last, ops.data(),
                                                                 counts.data(), bucket, marks.data());
                    cudaCheckError();
                    auto count = read_device(counts.data(), 8);
                    auto mark  = read_device(marks.data(), 8);
                    bool valid = true;
                    for (size_t i = 0; i < 8; ++i)
                    {
                        bool selected = i >= first && i < last && errors[i] > fraction * maximum;
                        valid &= count[i] == bucket + unsigned(selected) && mark[i] == uint8_t(selected);
                    }
                    require_all(valid, "Native non-binary octree threshold differs from refineFraction");
                }
    }
}

void check_empty_refiners(int rank)
{
    using namespace amr_test;
    using Tet = mars::amr::OctreeAlignedRefine<KeyType, RealType, mars::TetTag>;
    using Hex = mars::amr::OctreeAlignedRefine<KeyType, RealType, mars::HexTag>;
    const std::vector<RealType> x{0, 1, 1, 0, 0, 1, 1, 0}, y{0, 0, 1, 1, 0, 0, 1, 1},
        z{0, 0, 0, 0, 1, 1, 1, 1};
    for (int mode = 0; mode < 3; ++mode)
    {
        size_t nodes = mode == 0 ? 0 : 8;
        size_t count = mode == 2 && rank == 0 ? 1 : 0;
        auto dx = upload(nodes ? x : std::vector<RealType>{});
        auto dy = upload(nodes ? y : std::vector<RealType>{});
        auto dz = upload(nodes ? z : std::vector<RealType>{});
        std::array<cstone::DeviceVector<KeyType>, 8> conn;
        for (size_t i = 0; i < conn.size(); ++i)
            conn[i] = upload(count ? std::vector<KeyType>{i} : std::vector<KeyType>{});
        auto marks = upload(std::vector<uint8_t>(count, 1));
        std::vector<RealType> values(nodes);
        for (size_t i = 0; i < nodes; ++i) values[i] = 1 + 2*x[i] - 3*y[i] + 4*z[i];
        auto field = upload(values);
        auto tet = Tet::refineLocal(conn[0].data(), conn[1].data(), conn[3].data(), conn[4].data(),
                                   dx.data(), dy.data(), dz.data(), marks.data(), count, nodes);
        auto hex = Hex::refineLocal(conn[0].data(), conn[1].data(), conn[2].data(), conn[3].data(),
                                   conn[4].data(), conn[5].data(), conn[6].data(), conn[7].data(),
                                   dx.data(), dy.data(), dz.data(), marks.data(), count, nodes);
        cstone::DeviceVector<RealType> tet_field, hex_field;
        Tet::transferSolution(conn[0].data(), conn[1].data(), conn[3].data(), conn[4].data(),
                              field.data(), tet, count, tet_field);
        Hex::transferSolution(conn[0].data(), conn[1].data(), conn[2].data(), conn[3].data(),
                              conn[4].data(), conn[5].data(), conn[6].data(), conn[7].data(),
                              field.data(), marks.data(), count, nodes, hex_field);
        auto check = [&](const auto& refined, const auto& transferred, size_t added)
        {
            bool valid = refined.numElements == 8 * count && refined.numMarked == count &&
                         refined.oldNumNodes == nodes && refined.numNodes == nodes + added * count &&
                         transferred.size() == refined.numNodes;
            auto rx = read_device(refined.d_x.data(), refined.numNodes);
            auto ry = read_device(refined.d_y.data(), refined.numNodes);
            auto rz = read_device(refined.d_z.data(), refined.numNodes);
            auto actual = read_device(transferred.data(), transferred.size());
            for (size_t i = 0; i < actual.size(); ++i)
                valid &= actual[i] == 1 + 2*rx[i] - 3*ry[i] + 4*rz[i];
            require_all(valid, "Empty-rank refinement or affine transfer changed data");
        };
        check(tet, tet_field, 6);
        check(hex, hex_field, 19);
    }
    mars::amr::AmrOctree<KeyType, RealType> octree;
    require_all(octree.markForRefinement(nullptr, 0).empty(), "Empty-rank marking produced data");
    cudaCheckCall(cudaDeviceSynchronize());
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks, local_rank, devices;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    MPI_Comm local;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &local);
    MPI_Comm_rank(local, &local_rank);
    cudaCheckCall(cudaGetDeviceCount(&devices));
    require_all(devices == 1 || local_rank < devices, "No device for local rank");
    cudaCheckCall(cudaSetDevice(devices == 1 ? 0 : local_rank));
    MPI_Comm_free(&local);
    char host[MPI_MAX_PROCESSOR_NAME], pci[32];
    int host_length, device;
    MPI_Get_processor_name(host, &host_length);
    cudaCheckCall(cudaGetDevice(&device));
#if defined(__HIPCC__)
    cudaCheckCall(hipDeviceGetPCIBusId(pci, sizeof(pci), device));
#else
    cudaCheckCall(cudaDeviceGetPCIBusId(pci, sizeof(pci), device));
#endif
    std::printf("[amr-safety-device] rank=%d host=%s pci=%s\n", rank, host, pci);
    std::fflush(stdout);
    {
        check_node_keys<uint32_t, float>();
        check_node_keys<uint32_t, double>();
        check_node_keys<uint64_t, float>();
        check_node_keys<uint64_t, double>();
        using namespace amr_test;
        using Transfer = mars::amr::SolutionTransfer<KeyType, RealType>;
        auto values = upload(std::vector<RealType>{1, 2, 4, 8});
        auto unchanged = Transfer::transferByParentageTet(
            nullptr, nullptr, nullptr, nullptr, nullptr, 2, values.data(), 4, nullptr, 0, 4);
        require_all(read_device(unchanged.data(), unchanged.size()) == std::vector<RealType>({1, 2, 4, 8}),
                    "Unrefined-rank field copy changed values");
        auto empty = Transfer::transferByParentageTet(
            nullptr, nullptr, nullptr, nullptr, nullptr, 0, nullptr, 0, nullptr, 0, 0);
        require_all(empty.empty(), "Empty transfer produced nodes");

        auto c0 = upload(std::vector<KeyType>{0, 0}), c1 = upload(std::vector<KeyType>{1, 1});
        auto c2 = upload(std::vector<KeyType>{2, 2}), c3 = upload(std::vector<KeyType>{3, 3});
        auto marks = upload(std::vector<uint8_t>{1, 1});
        auto edges = upload(std::vector<uint64_t>{1ull << 32, 2ull << 32, (2ull << 32) | 1,
                                                 3ull << 32, (3ull << 32) | 1, (3ull << 32) | 2});
        auto refined = Transfer::transferByParentageTet(
            c0.data(), c1.data(), c2.data(), c3.data(), marks.data(), 2, values.data(), 4, edges.data(), 6, 10);
        require_all(read_device(refined.data(), refined.size()) ==
                        std::vector<RealType>({1, 2, 4, 8, 1.5, 2.5, 3, 4.5, 5, 6}),
                    "Shared-edge interpolation changed values");
        for (bool staged : {false, true})
            amr_test::check_nonconforming_topology<mars::amr::LeppDeviceExecution>(MPI_COMM_WORLD, staged);
        check_marking<float>(rank, ranks);
        check_marking<double>(rank, ranks);
        check_empty_refiners(rank);
        for (bool staged : {false, true})
            amr_test::check_field_sources<mars::amr::LeppDeviceExecution>(MPI_COMM_WORLD, staged);
        cudaCheckCall(cudaDeviceSynchronize());
        std::printf("[amr-safety] rank=%d PASS\n", rank);
    }
    MPI_Finalize();
}
