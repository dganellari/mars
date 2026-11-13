#include "mars_amr_tet_audit.hpp"
#include "mars_amr_keyed_position.hpp"
#include "../utils/mars_domain_utils.hpp"
#include "../utils/mars_gpu_reduce.hpp"
#include <thrust/device_ptr.h>
#include <thrust/iterator/counting_iterator.h>
#include <mpi.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <vector>

namespace mars::amr
{
template<class KeyType, class RealType>
__global__ void audit_tets_kernel(const KeyType* a,
                                  const KeyType* b,
                                  const KeyType* c,
                                  const KeyType* d,
                                  size_t begin,
                                  size_t count,
                                  KeyedPosition<KeyType, RealType> position,
                                  const RealType* x,
                                  const RealType* y,
                                  const RealType* z,
                                  size_t num_nodes,
                                  bool keyed_geometry,
                                  const KeyType* node_keys,
                                  TetAuditRecord<KeyType>* records)
{
    size_t i = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (i >= count) return;
    KeyType keys[4] = {a[begin + i], b[begin + i], c[begin + i], d[begin + i]};
    double p[4][3];
    for (int j = 0; j < 4; ++j)
    {
        RealType px, py, pz;
        if (x)
        {
            size_t node = size_t(keys[j]);
            if (node >= num_nodes)
            {
                records[i]         = {};
                records[i].invalid = 1;
                return;
            }
            px = x[node];
            py = y[node];
            pz = z[node];
            if (!::isfinite(px) || !::isfinite(py) || !::isfinite(pz))
            {
                records[i]         = {};
                records[i].invalid = 1;
                return;
            }
            keys[j] = node_keys ? node_keys[node] : position.encode(px, py, pz);
        }
        if (!x || keyed_geometry) position(keys[j], px, py, pz);
        p[j][0] = px;
        p[j][1] = py;
        p[j][2] = pz;
    }
    records[i] = audit_record(keys, p);
    records[i].roundoff_scale *=
        double(std::numeric_limits<RealType>::epsilon()) / std::numeric_limits<double>::epsilon();
}

template<class KeyType, class RealType>
TetAuditSnapshot<KeyType> audit_tets(const char* stage,
                                     std::array<const KeyType*, 4> conn,
                                     size_t begin,
                                     size_t end,
                                     const cstone::Box<RealType>& box,
                                     const RealType* x,
                                     const RealType* y,
                                     const RealType* z,
                                     size_t num_nodes,
                                     bool keyed_geometry,
                                     const KeyType* node_keys)
{
    TetAuditSnapshot<KeyType> snapshot;
    snapshot.stage     = stage;
    const size_t count = end - begin;
    snapshot.records.resize(count * sizeof(TetAuditRecord<KeyType>));
    auto* records = reinterpret_cast<TetAuditRecord<KeyType>*>(snapshot.records.data());
    if (count)
    {
        audit_tets_kernel<<<unsigned((count + 255) / 256), 256>>>(conn[0], conn[1], conn[2], conn[3], begin, count,
                                                                  KeyedPosition<KeyType, RealType>{box}, x, y, z,
                                                                  num_nodes, keyed_geometry, node_keys, records);
        cudaCheckError();
    }
    // Attribute asynchronous errors to this stage, before entering its collectives.
    cudaCheckCall(cudaDeviceSynchronize());
    TetAuditTotals local;
    if (count)
        local =
            mars::transformReduce(thrust::device_pointer_cast(records), thrust::device_pointer_cast(records) + count,
                                  TetAuditValue<KeyType>{}, TetAuditTotals{}, TetAuditFold{});
    uint64_t values[5] = {local.count, local.hash_a, local.hash_b, local.invalid, local.inverted};
    MPI_Allreduce(MPI_IN_PLACE, values, 5, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    double volumes[2] = {local.volume, local.absolute_volume};
    MPI_Allreduce(MPI_IN_PLACE, volumes, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &local.min_abs_volume, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    snapshot.totals = {values[0], values[1],  values[2],  values[3],
                       values[4], volumes[0], volumes[1], local.min_abs_volume};
    int rank        = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        std::printf("[amr-mesh-audit] stage=%s count=%llu hash_a=%016llx hash_b=%016llx invalid=%llu "
                    "inverted=%llu volume=%.17g absolute_volume=%.17g min_abs_volume=%.17g\n",
                    stage, (unsigned long long)values[0], (unsigned long long)values[1], (unsigned long long)values[2],
                    (unsigned long long)values[3], (unsigned long long)values[4], volumes[0], volumes[1],
                    local.min_abs_volume);
        std::fflush(stdout);
    }
    return snapshot;
}

template<class KeyType>
void audit_dump(const TetAuditSnapshot<KeyType>& snapshot, const char* side)
{
    const char* directory = std::getenv("MARS_AMR_AUDIT_DIR");
    if (!directory || !*directory) return;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    std::string filename =
        std::string(directory) + "/" + side + "-" + snapshot.stage + "-rank" + std::to_string(rank) + ".csv";
    std::ofstream out(filename);
    if (!out)
    {
        std::fprintf(stderr, "[amr-mesh-audit] cannot write %s\n", filename.c_str());
        return;
    }
    // Full records cross to the host only on failure, and only with the explicit dump directory.
    std::vector<TetAuditRecord<KeyType>> records(snapshot.records.size() / sizeof(TetAuditRecord<KeyType>));
    if (!records.empty())
        cudaCheckCall(
            cudaMemcpy(records.data(), snapshot.records.data(), snapshot.records.size(), cudaMemcpyDeviceToHost));
    out << "key0,key1,key2,key3,volume,invalid\n" << std::setprecision(17);
    for (const auto& r : records)
        out << r.keys[0] << ',' << r.keys[1] << ',' << r.keys[2] << ',' << r.keys[3] << ',' << r.volume << ','
            << r.invalid << '\n';
    out.close();
    if (!out) std::fprintf(stderr, "[amr-mesh-audit] incomplete dump %s\n", filename.c_str());
}

template<class KeyType>
void audit_require(bool condition,
                   const char* reason,
                   const TetAuditSnapshot<KeyType>& current,
                   const TetAuditSnapshot<KeyType>* previous)
{
    int failed = !condition;
    MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (!failed) return;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::fprintf(stderr, "[amr-mesh-audit] FAIL stage=%s reason=%s\n", current.stage.c_str(), reason);
    if (previous) audit_dump(*previous, "before");
    audit_dump(current, "after");
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 77);
    std::abort();
}

template<class KeyType>
void audit_compare(const TetAuditSnapshot<KeyType>& before, const TetAuditSnapshot<KeyType>& after)
{
    const auto& a    = before.totals;
    const auto& b    = after.totals;
    const double tol = audit_volume_tolerance(a, b);
    bool same        = audit_same_cells(a, b) && a.inverted == b.inverted && std::abs(a.volume - b.volume) <= tol &&
                std::abs(a.absolute_volume - b.absolute_volume) <= tol;
    audit_require(same, "owned cell set or canonical geometry changed", after, &before);
}

template<class KeyType>
void audit_red_children(const TetAuditSnapshot<KeyType>& parents,
                        const TetAuditSnapshot<KeyType>& children,
                        const uint8_t* marks,
                        const uint64_t* marked_prefix)
{
    size_t count       = parents.records.size() / sizeof(TetAuditRecord<KeyType>);
    size_t child_count = children.records.size() / sizeof(TetAuditRecord<KeyType>);
    auto* p            = reinterpret_cast<const TetAuditRecord<KeyType>*>(parents.records.data());
    auto* c            = reinterpret_cast<const TetAuditRecord<KeyType>*>(children.records.data());
    uint64_t bad       = mars::transformReduce(
        thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(count),
        [=] __host__ __device__(size_t i) -> uint64_t
        {
            size_t first = i + (marked_prefix ? 7 * marked_prefix[i] : 0);
            return !audit_red_children_valid(p[i], c, child_count, first, marks[i] > 0);
        },
        uint64_t(0), thrust::plus<uint64_t>());
    MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    audit_require(bad == 0, "child volume or orientation differs from its parent", children, &parents);
}

#define MARS_INSTANTIATE_TET_AUDIT(Key, Real)                                                                          \
    template TetAuditSnapshot<Key> audit_tets<Key, Real>(const char*, std::array<const Key*, 4>, size_t, size_t,       \
                                                         const cstone::Box<Real>&, const Real*, const Real*,           \
                                                         const Real*, size_t, bool, const Key*);
MARS_INSTANTIATE_TET_AUDIT(uint32_t, float)
MARS_INSTANTIATE_TET_AUDIT(uint32_t, double)
MARS_INSTANTIATE_TET_AUDIT(uint64_t, float)
MARS_INSTANTIATE_TET_AUDIT(uint64_t, double)
#undef MARS_INSTANTIATE_TET_AUDIT

template void
audit_require<uint32_t>(bool, const char*, const TetAuditSnapshot<uint32_t>&, const TetAuditSnapshot<uint32_t>*);
template void
audit_require<uint64_t>(bool, const char*, const TetAuditSnapshot<uint64_t>&, const TetAuditSnapshot<uint64_t>*);
template void audit_compare<uint32_t>(const TetAuditSnapshot<uint32_t>&, const TetAuditSnapshot<uint32_t>&);
template void audit_compare<uint64_t>(const TetAuditSnapshot<uint64_t>&, const TetAuditSnapshot<uint64_t>&);
template void audit_red_children<uint32_t>(const TetAuditSnapshot<uint32_t>&,
                                           const TetAuditSnapshot<uint32_t>&,
                                           const uint8_t*,
                                           const uint64_t*);
template void audit_red_children<uint64_t>(const TetAuditSnapshot<uint64_t>&,
                                           const TetAuditSnapshot<uint64_t>&,
                                           const uint8_t*,
                                           const uint64_t*);
} // namespace mars::amr
