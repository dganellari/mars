#pragma once

// The unknowns of a nodal field on a distributed mesh.
//
// Every rank stores one value per node slot: its owned nodes and ghost copies
// of nodes owned by other ranks. On a periodic mesh a point on the box
// boundary also has several slots on the same rank or on different ranks: the
// master on a min face and its slaves (the same point) on the max faces. The
// unknown (DOF) of a point lives in ONE slot: the owned slot that is not a
// periodic slave. Two maps connect slots and DOFs:
//
//   prolong  (P)   copy each DOF value into every slot of its point
//   restrict (P^T) add the per-slot contributions of an element scatter into
//                  the DOF slot (the other slots end up holding nothing useful)
//
// A discrete operator is then restrict(A_local(prolong(x))) = P^T A P, with
// A_local the element scatter over this rank's own elements. Without a
// periodic map the only copies are ghosts, and P^T / P are the cstone reverse
// and forward node halo exchanges.

#include "backend/distributed/unstructured/fem/mars_periodic_bc.hpp"

#include <thrust/count.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>

#include <algorithm>
#include <cstdio>
#include <type_traits>
#include <utility>
#include <vector>
#include <mpi.h>

namespace mars
{
namespace fem
{

template<typename RealType>
inline MPI_Datatype mpiDatatype()
{
    return std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
}

// A slot is a DOF when it is owned and is not a periodic slave. partner == nullptr: no periodic pairs.
// remote marks slaves paired only through the cross-rank tables (their partner is -1).
template<typename IndexType>
__global__ void dofMaskKernel(const IndexType* partner, const uint8_t* remote, const uint8_t* ownership, size_t n,
                              uint8_t* isDof)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    bool slave = (partner != nullptr && partner[i] >= 0) || (remote != nullptr && remote[i] != 0);
    isDof[i]   = (ownership[i] == 1 && !slave) ? 1 : 0;
}

template<typename IndexType>
__global__ void dofFlagSlotsKernel(const IndexType* ids, size_t count, uint8_t* flag)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k < count) flag[ids[k]] = 1;
}

// Two dot products over DOF slots in one pass: block partial sums, reduced in
// a second single-block pass so the result does not depend on scheduling.
// b == nullptr sums a; c == nullptr skips the second product.
template<typename RealType, int BlockSize>
__global__ void dofDot2PartialKernel(const RealType* a, const RealType* b, const RealType* c, const RealType* d,
                                     const uint8_t* isDof, size_t n, RealType* partial)
{
    __shared__ RealType s0[BlockSize];
    __shared__ RealType s1[BlockSize];
    RealType ab = 0, cd = 0;
    for (size_t i = size_t(blockIdx.x) * BlockSize + threadIdx.x; i < n; i += size_t(gridDim.x) * BlockSize)
    {
        if (!isDof[i]) continue;
        ab += b ? a[i] * b[i] : a[i];
        if (c) cd += c[i] * d[i];
    }
    s0[threadIdx.x] = ab;
    s1[threadIdx.x] = cd;
    __syncthreads();
    for (unsigned stride = BlockSize / 2; stride > 0; stride /= 2)
    {
        if (threadIdx.x < stride)
        {
            s0[threadIdx.x] += s0[threadIdx.x + stride];
            s1[threadIdx.x] += s1[threadIdx.x + stride];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
    {
        partial[2 * blockIdx.x]     = s0[0];
        partial[2 * blockIdx.x + 1] = s1[0];
    }
}

template<typename RealType, int BlockSize>
__global__ void dofDot2FinalKernel(RealType* partial, int numPartials)
{
    __shared__ RealType s0[BlockSize];
    __shared__ RealType s1[BlockSize];
    RealType ab = 0, cd = 0;
    for (int i = threadIdx.x; i < numPartials; i += BlockSize)
    {
        ab += partial[2 * i];
        cd += partial[2 * i + 1];
    }
    s0[threadIdx.x] = ab;
    s1[threadIdx.x] = cd;
    __syncthreads();
    for (unsigned stride = BlockSize / 2; stride > 0; stride /= 2)
    {
        if (threadIdx.x < stride)
        {
            s0[threadIdx.x] += s0[threadIdx.x + stride];
            s1[threadIdx.x] += s1[threadIdx.x + stride];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
    {
        partial[0] = s0[0];
        partial[1] = s1[0];
    }
}

template<typename RealType>
__global__ void dofShiftKernel(RealType* v, RealType shift, size_t n)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) v[i] -= shift;
}

// Setup-time guarantee prolong and restrict rely on: every owned slave is paired
// exactly once, either with a final master (a node on no max face) owned on this
// rank, or through the cross-rank table. A violation means the mesh or the halo
// cannot represent a periodic point, so stop instead of running with a split point.
template<typename KeyType, typename RealType, typename DomainT>
void checkPeriodicPairing(const DomainT& domain, const PeriodicMap<KeyType, RealType>& map)
{
    const size_t n        = domain.getNodeCount();
    const auto& table     = map.cross_.d_sendOwnedSlaveIds_;
    const int* partner    = map.d_periodicPartner.data();
    const uint8_t* mask   = map.d_periodicMask.data();
    const uint8_t* own    = domain.getNodeOwnershipMap().data();
    cstone::DeviceVector<uint8_t> inTable(n);
    thrust::fill(thrust::device, inTable.data(), inTable.data() + n, uint8_t(0));
    if (table.size() > 0)
        dofFlagSlotsKernel<int><<<int((table.size() + 255) / 256), 256>>>(table.data(), table.size(), inTable.data());
    cudaCheckError();
    const uint8_t* tab = inTable.data();

    auto first = thrust::counting_iterator<size_t>(0);
    long long bad[2];
    bad[0] = thrust::count_if(thrust::device, first, first + n, [partner, mask, own, tab] __device__(size_t i) -> bool {
        if (own[i] != 1 || mask[i] == 0) return tab[i] != 0;
        int m         = partner[i];
        bool local    = m >= 0 && mask[m] == 0 && own[m] == 1;
        bool viaTable = tab[i] != 0 && (m < 0 || mask[m] == 0);
        return local == viaTable;
    });
    long long tableSize = (long long)table.size();
    long long flagged   = thrust::count(thrust::device, tab, tab + n, uint8_t(1));
    bad[1]              = tableSize != flagged ? 1 : 0;
    long long badGlobal[2];
    MPI_Allreduce(bad, badGlobal, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (badGlobal[0] > 0 || badGlobal[1] > 0)
    {
        if (domain.rank() == 0)
            std::fprintf(stderr,
                         "DofSpace: %lld owned slaves are not paired exactly once with a final master, %lld ranks "
                         "list a slave twice in the cross-rank table. Check the periodic box bounds and faceEps.\n",
                         badGlobal[0], badGlobal[1]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

template<typename KeyType, typename RealType, typename DomainT>
class DofSpace
{
public:
    using Vector = cstone::DeviceVector<RealType>;
    using Map    = PeriodicMap<KeyType, RealType>;
    static constexpr int ReduceBlock = 256;
    static constexpr int ReduceGrid  = 256;

    // periodic == nullptr: no periodic pairs.
    DofSpace(const DomainT& domain, const Map* periodic, int blockSize = 256)
        : domain_(domain)
        , map_(periodic)
        , n_(domain.getNodeCount())
        , blockSize_(blockSize)
        , isDof_(n_)
        , partial_(2 * ReduceGrid)
        , hostPartial_(2)
    {
        const auto& own = domain_.getNodeOwnershipMap();
        const uint8_t* remote = map_ && map_->d_remoteSlave.size() == n_ ? map_->d_remoteSlave.data() : nullptr;
        dofMaskKernel<int><<<grid(), blockSize_>>>(partner(), remote, own.data(), n_, isDof_.data());
        cudaCheckError();
        if (map_) checkPeriodicPairing<KeyType, RealType>(domain_, *map_);

        long long local = thrust::count(thrust::device, isDof_.data(), isDof_.data() + n_, uint8_t(1));
        MPI_Allreduce(&local, &numDofs_, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    }

    // P: every slot of a point takes the DOF value. Owned slaves first (locally,
    // or from the master's owner rank), then the halo refreshes all ghosts.
    void prolong(Vector& v) const
    {
        if (map_)
        {
            periodicBroadcastSameRankKernel<RealType><<<grid(), blockSize_>>>(
                partner(), domain_.getNodeOwnershipMap().data(), n_, v.data());
            cudaCheckError();
            crossRankPeriodicBroadcast<KeyType, RealType>(*map_, v);
        }
        domain_.exchangeNodeHalo(v);
    }

    // P^T: the reverse halo completes every owned slot, then every owned slave
    // adds into its master (locally, or on the master's owner rank). Must be
    // applied to a scatter over this rank's own elements.
    void restrict(Vector& acc) const
    {
        domain_.reverseExchangeNodeHaloAdd(acc);
        if (!map_) return;
        periodicPairSumKernel<RealType><<<grid(), blockSize_>>>(partner(), domain_.getNodeOwnershipMap().data(), n_,
                                                               acc.data());
        cudaCheckError();
        crossRankPeriodicPairSum<KeyType, RealType>(*map_, acc, /*broadcastBack=*/false);
    }

    // Inner products over DOFs: each point counts once globally.
    RealType dot(const Vector& a, const Vector& b) const { return reduce2(a.data(), b.data(), a.data(), b.data()).first; }
    RealType sum(const Vector& a) const { return reduce2(a.data(), nullptr, nullptr, nullptr).first; }

    // Shift v by its mean over the DOFs. Applied to every slot, so a prolonged field stays prolonged.
    void removeMean(Vector& v) const
    {
        RealType mean = sum(v) / RealType(numDofs_);
        dofShiftKernel<RealType><<<grid(), blockSize_>>>(v.data(), mean, n_);
        cudaCheckError();
    }

    long long numDofs() const { return numDofs_; }
    const uint8_t* isDof() const { return isDof_.data(); }
    const int* partner() const { return map_ ? map_->d_periodicPartner.data() : nullptr; }
    bool periodic() const { return map_ != nullptr; }
    size_t numSlots() const { return n_; }

private:
    int grid() const { return std::max(1, int((n_ + blockSize_ - 1) / blockSize_)); }

    std::pair<RealType, RealType> reduce2(const RealType* a, const RealType* b, const RealType* c,
                                          const RealType* d) const
    {
        dofDot2PartialKernel<RealType, ReduceBlock><<<ReduceGrid, ReduceBlock>>>(a, b, c, d, isDof_.data(), n_,
                                                                                 partial_.data());
        dofDot2FinalKernel<RealType, ReduceBlock><<<1, ReduceBlock>>>(partial_.data(), ReduceGrid);
        cudaCheckError();
        RealType local[2];
        cudaMemcpy(local, partial_.data(), 2 * sizeof(RealType), cudaMemcpyDeviceToHost);
        MPI_Allreduce(local, hostPartial_.data(), 2, mpiDatatype<RealType>(), MPI_SUM, MPI_COMM_WORLD);
        return {hostPartial_[0], hostPartial_[1]};
    }

    const DomainT& domain_;
    const Map* map_;
    size_t n_;
    int blockSize_;
    cstone::DeviceVector<uint8_t> isDof_;
    mutable Vector partial_;
    mutable std::vector<RealType> hostPartial_;
    long long numDofs_ = 0;
};

} // namespace fem
} // namespace mars
