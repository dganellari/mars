#include "cstone/domain/domain.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_tree_layout.hpp"
#include <iostream>
#include <vector>

template<class T>
std::vector<T> download(const T* device, size_t size)
{
    std::vector<T> host(size);
    cstone::memcpyD2HAsync(cstone::execution::gpuDefaultStream, device, size, host.data());
    cstone::syncGpu(cstone::execution::gpuDefaultStream);
    return host;
}

template<class KeyType, class RealType>
void check_tree(int rank, int ranks)
{
    // The pinned GPU domain supports 32-bit keys with float coordinates; error precision is independent.
    using CoordType = std::conditional_t<sizeof(KeyType) == 4, float, RealType>;
    using Vec = cstone::DeviceVector<CoordType>;
    using Keys = cstone::DeviceVector<KeyType>;
    const auto exec = cstone::execution::gpuDefaultStream;
    for (int rebuild = 0; rebuild < 3; ++rebuild)
    {
        const int n = 8 + rebuild * 4;
        std::vector<CoordType> hx, hy, hz, hh;
        for (int i = 0; i < n * n * n; ++i)
        {
            if (i % ranks != rank) continue;
            const int j = i / 3;
            hx.push_back((j % n + CoordType(.5)) / n);
            hy.push_back((j / n % n + CoordType(.5)) / n);
            hz.push_back((j / (n * n) + CoordType(.5)) / n);
            hh.push_back(CoordType(.08));
        }
        Vec x, y, z, h, s1, s2, s3;
        x = hx; y = hy; z = hz; h = hh;
        Keys keys(x.size(), KeyType(0));
        cstone::Domain<KeyType, CoordType, cstone::execution::Gpu> domain(exec, rank, ranks, 32, 8, .5,
                                                                       MPI_COMM_WORLD, cstone::Box<CoordType>(0, 1));
        for (int sync = 0; sync < 2; ++sync)
        {
            domain.sync(keys, x, y, z, h, std::tie(), std::tie(s1, s2, s3));
            const auto focus = domain.octreeProperties();
            const auto global = domain.globalTree();
            const auto leaf_keys = download(focus.leaves, focus.numLeafNodes + 1);
            const auto global_keys = download(global.leaves, global.numLeafNodes + 1);
            const auto sorted = download(keys.data(), keys.size());
            const size_t first = domain.startIndex(), last = domain.endIndex();
            int invalid = 0;
            cstone::DeviceVector<int> levels(keys.size(), -1);
            mars::amr::assign_leaf_levels<<<focus.numLeafNodes, 256>>>(focus, levels.data());
            checkGpuErrors(cudaGetLastError());
            const auto actual_levels = download(levels.data(), levels.size());
            for (size_t e = 0; e < sorted.size(); ++e)
            {
                const auto leaf = std::upper_bound(leaf_keys.begin(), leaf_keys.end(), sorted[e]) - leaf_keys.begin() - 1;
                if (leaf < 0 || leaf >= focus.numLeafNodes) { ++invalid; continue; }
                int depth = 0;
                const auto range = leaf_keys[leaf + 1] - leaf_keys[leaf];
                for (auto span = cstone::nodeRange<KeyType>(0); span > range; span >>= 3) ++depth;
                invalid += actual_levels[e] != depth;
            }
            std::vector<RealType> errors(keys.size(), RealType(1000)), local(global.numLeafNodes, 0), full;
            RealType maximum = 0;
            for (size_t e = first; e < last; ++e)
            {
                errors[e] = RealType((sorted[e] / 997) % 101) / 100;
                auto leaf = std::upper_bound(global_keys.begin(), global_keys.end(), sorted[e]) - global_keys.begin() - 1;
                local[leaf] = std::max(local[leaf], errors[e]);
                maximum = std::max(maximum, errors[e]);
            }
            auto type = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
            full.resize(local.size());
            MPI_Allreduce(local.data(), full.data(), int(local.size()), type, MPI_MAX, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, type, MPI_MAX, MPI_COMM_WORLD);
            cstone::DeviceVector<RealType> d_error; d_error = errors;
            for (int threads : {1, 64, 256})
            {
                constexpr unsigned bucket = 64;
                cstone::DeviceVector<unsigned> counts(global.numLeafNodes);
                cstone::DeviceVector<cstone::TreeNodeIndex> ops(global.numLeafNodes + 1);
                mars::amr::weight_owned_leaves<<<global.numLeafNodes, threads>>>(global.leaves, global.numLeafNodes,
                    keys.data(), first, last, d_error.data(), maximum, bucket, counts.data());
                checkGpuErrors(cudaGetLastError());
                auto actual_counts = download(counts.data(), counts.size());
                std::vector<unsigned> reference(global.numLeafNodes);
                for (int leaf = 0; leaf < global.numLeafNodes; ++leaf)
                {
                    auto weighted = [=](RealType error) {
                        const RealType scale = maximum > 0 ? RealType(1) + error / maximum * bucket : RealType(1);
                        return std::max(bucket, unsigned(bucket * scale));
                    };
                    invalid += actual_counts[leaf] != weighted(local[leaf]);
                    reference[leaf] = weighted(full[leaf]);
                }
                cstone::computeNodeOpsGpu(exec, global.leaves, global.numLeafNodes, counts.data(), bucket, ops.data());
                cstone::DeviceVector<uint8_t> marks(keys.size(), 0);
                mars::amr::mark_owned_leaves<<<global.numLeafNodes, threads>>>(global.leaves, global.numLeafNodes,
                    keys.data(), first, last, ops.data(), marks.data());
                checkGpuErrors(cudaGetLastError());
                const auto actual = download(marks.data(), marks.size());
                for (size_t e = 0; e < sorted.size(); ++e)
                {
                    uint8_t expected = 0;
                    if (e >= first && e < last)
                    {
                        auto leaf = std::upper_bound(global_keys.begin(), global_keys.end(), sorted[e]) - global_keys.begin() - 1;
                        const auto op = cstone::calculateNodeOp(global_keys.data(), leaf, reference.data(), bucket);
                        expected = op >= 8 ? 1 : op == 0 ? 2 : 0;
                    }
                    invalid += actual[e] != expected;
                }
            }
            MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (invalid) throw std::runtime_error("Device tree layout or owned-leaf marking differs from reference");
            if (rank == 0) std::cout << "tree-gpu key=" << sizeof(KeyType)*8 << " coord=" << sizeof(CoordType)*8 << " real=" << sizeof(RealType)*8
                << " rebuild=" << rebuild << " sync=" << sync << " ranks=" << ranks << " PASS\n";
            hx = download(x.data(), x.size());
            for (size_t e = first; e < last; ++e) hx[e] *= CoordType(.7);
            x = hx;
        }
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks, local_rank, devices;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    MPI_Comm local;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &local);
    MPI_Comm_rank(local, &local_rank);
    checkGpuErrors(cudaGetDeviceCount(&devices));
    if (devices == 0 || (devices != 1 && local_rank >= devices)) MPI_Abort(MPI_COMM_WORLD, 2);
    const int device = devices == 1 ? 0 : local_rank;
    checkGpuErrors(cudaSetDevice(device));
    char host[MPI_MAX_PROCESSOR_NAME], pci[32]; int length;
    MPI_Get_processor_name(host, &length);
#if defined(__HIPCC__)
    checkGpuErrors(hipDeviceGetPCIBusId(pci, sizeof(pci), device));
#else
    checkGpuErrors(cudaDeviceGetPCIBusId(pci, sizeof(pci), device));
#endif
    std::cout << "[tree-device] rank=" << rank << " host=" << host << " pci=" << pci << std::endl;
    MPI_Comm_free(&local);
    try
    {
        check_tree<uint32_t, float>(rank, ranks);
        check_tree<uint32_t, double>(rank, ranks);
        check_tree<uint64_t, float>(rank, ranks);
        check_tree<uint64_t, double>(rank, ranks);
    }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; MPI_Abort(MPI_COMM_WORLD, 1); }
    MPI_Finalize();
}
