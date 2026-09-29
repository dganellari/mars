#pragma once

// The reduced periodic space: one unknown per periodic point, for every field.
//
// A mesh of a periodic box keeps separate node slots on opposite faces: the
// master on a min face and its slaves (the same physical point) on the max
// faces. On several ranks a point also has ghost copies. The unknown of a
// point lives in ONE slot, the master slot on the master's owner rank; that
// slot is its DOF. Two maps connect the slots and the DOFs:
//
//   prolong  (P)   copy the DOF value into every slot of the point
//   restrict (P^T) add the per-slot contributions of an element scatter into
//                  the DOF slot (the other slots end up holding nothing useful)
//
// Every discrete operator is then restrict(A_elem(prolong(x))) = P^T A P, with
// A_elem the plain element scatter over owned elements. P and P^T are exact
// transposes, so a symmetric element operator stays symmetric, and the
// gradient and the divergence remain transposes of each other. This is what
// makes the projection D u = 0 hold exactly on any number of ranks.
//
// Communication is the cstone node halo (owner <-> ghost) plus the cross-rank
// pair tables of PeriodicMap (a slave owned on one rank whose master is owned
// on another rank; cstone cannot link them because their SFC keys differ).

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

// A slot is a DOF when it is owned and is not a periodic slave.
__global__ void periodicDofMaskKernel(const int* partner, const uint8_t* ownership, size_t n, uint8_t* isDof)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    isDof[i] = (ownership[i] == 1 && partner[i] < 0) ? 1 : 0;
}

// Two dot products over DOF slots in one pass: block partial sums, reduced in
// a second single-block pass so the result does not depend on scheduling.
// b == nullptr sums a; c == nullptr skips the second product.
template<typename RealType, int BlockSize>
__global__ void periodicDot2PartialKernel(const RealType* a, const RealType* b,
                                          const RealType* c, const RealType* d,
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
__global__ void periodicDot2FinalKernel(RealType* partial, int numPartials)
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
__global__ void periodicShiftKernel(RealType* v, RealType shift, size_t n)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    v[i] -= shift;
}

// Setup-time guarantees prolong and restrict rely on. A violation means the
// mesh or the halo cannot represent a periodic point, so stop instead of
// running with a silently split point.
template<typename KeyType, typename RealType, typename DomainT>
void checkPeriodicPairing(const DomainT& domain, const PeriodicMap<KeyType, RealType>& map)
{
    const int* partner  = map.d_periodicPartner.data();
    const uint8_t* mask = map.d_periodicMask.data();
    const uint8_t* own  = domain.getNodeOwnershipMap().data();
    auto first          = thrust::counting_iterator<size_t>(0);
    auto last           = thrust::counting_iterator<size_t>(domain.getNodeCount());

    // every owned slave must reach its final master, a node on no max face
    long long unpaired = thrust::count_if(thrust::device, first, last, [partner, mask, own] __device__(size_t i) {
        if (own[i] != 1 || mask[i] == 0) return false;
        int m = partner[i];
        return m < 0 || mask[m] != 0;
    });
    // every owned slave with a remote master must be in the cross-rank table
    long long remote = thrust::count_if(thrust::device, first, last, [partner, own] __device__(size_t i) {
        return own[i] == 1 && partner[i] >= 0 && own[partner[i]] != 1;
    });
    long long bad[2] = {unpaired, remote != (long long)map.cross_.d_sendOwnedSlaveIds_.size() ? 1LL : 0LL};
    long long badGlobal[2];
    MPI_Allreduce(bad, badGlobal, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (badGlobal[0] > 0 || badGlobal[1] > 0)
    {
        if (domain.rank() == 0)
            std::fprintf(stderr,
                         "PeriodicSpace: %lld owned slaves have no master in the local halo, %lld ranks have an "
                         "incomplete cross-rank pair table. Check the periodic box bounds and faceEps.\n",
                         badGlobal[0], badGlobal[1]);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

template<typename KeyType, typename RealType, typename DomainT>
class PeriodicSpace
{
public:
    using Vector = cstone::DeviceVector<RealType>;
    static constexpr int ReduceBlock = 256;
    static constexpr int ReduceGrid  = 256;

    PeriodicSpace(const DomainT& domain, const PeriodicMap<KeyType, RealType>& map, int blockSize = 256)
        : domain_(domain)
        , map_(map)
        , n_(domain.getNodeCount())
        , blockSize_(blockSize)
        , isDof_(n_)
        , partial_(2 * ReduceGrid)
        , hostPartial_(2)
    {
        const auto& own = domain_.getNodeOwnershipMap();
        periodicDofMaskKernel<<<grid(), blockSize_>>>(map_.d_periodicPartner.data(), own.data(), n_,
                                                                isDof_.data());
        cudaCheckError();
        checkPeriodicPairing<KeyType, RealType>(domain_, map_);

        long long local = thrust::count(thrust::device, isDof_.data(), isDof_.data() + n_, uint8_t(1));
        MPI_Allreduce(&local, &numDofs_, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    }

    // P: every slot of a periodic point takes the DOF value. Owned slaves first
    // (locally, or from the master's owner rank), then the halo refreshes all
    // ghosts, including ghost copies of slaves.
    void prolong(Vector& v) const
    {
        periodicBroadcastSameRankKernel<RealType><<<grid(), blockSize_>>>(
            map_.d_periodicPartner.data(), domain_.getNodeOwnershipMap().data(), n_, v.data());
        cudaCheckError();
        crossRankPeriodicBroadcast<KeyType, RealType>(map_, v);
        domain_.exchangeNodeHalo(v);
    }

    // P^T: the reverse halo completes every owned slot, then every owned slave
    // adds into its master (locally, or on the master's owner rank). Must be
    // applied to a scatter over OWNED elements only.
    void restrict(Vector& acc) const
    {
        domain_.reverseExchangeNodeHaloAdd(acc);
        periodicPairSumKernel<RealType><<<grid(), blockSize_>>>(
            map_.d_periodicPartner.data(), domain_.getNodeOwnershipMap().data(), n_, acc.data());
        cudaCheckError();
        crossRankPeriodicPairSum<KeyType, RealType>(map_, acc, /*broadcastBack=*/false);
    }

    // Inner products over DOFs: each periodic point is counted once globally.
    RealType dot(const Vector& a, const Vector& b) const { return dot2(a, b, a, b).first; }

    // (a, b) and (c, d) with one reduction pass and one Allreduce.
    std::pair<RealType, RealType> dot2(const Vector& a, const Vector& b, const Vector& c, const Vector& d) const
    {
        return reduce2(a.data(), b.data(), c.data(), d.data());
    }

    RealType sum(const Vector& a) const { return reduce2(a.data(), nullptr, nullptr, nullptr).first; }

    // Pressure is defined up to a constant: fix it by a zero mean over DOFs.
    // The shift is applied to every slot, so a prolonged field stays prolonged.
    void removeMean(Vector& v) const
    {
        RealType mean = sum(v) / RealType(numDofs_);
        periodicShiftKernel<RealType><<<grid(), blockSize_>>>(v.data(), mean, n_);
        cudaCheckError();
    }

    long long numDofs() const { return numDofs_; }
    const uint8_t* isDof() const { return isDof_.data(); }
    size_t numSlots() const { return n_; }
    int grid() const { return std::max(1, int((n_ + blockSize_ - 1) / blockSize_)); }

private:
    std::pair<RealType, RealType> reduce2(const RealType* a, const RealType* b, const RealType* c,
                                          const RealType* d) const
    {
        periodicDot2PartialKernel<RealType, ReduceBlock><<<ReduceGrid, ReduceBlock>>>(a, b, c, d, isDof_.data(), n_,
                                                                                      partial_.data());
        periodicDot2FinalKernel<RealType, ReduceBlock><<<1, ReduceBlock>>>(partial_.data(), ReduceGrid);
        cudaCheckError();
        RealType local[2];
        cudaMemcpy(local, partial_.data(), 2 * sizeof(RealType), cudaMemcpyDeviceToHost);
        MPI_Allreduce(local, hostPartial_.data(), 2, mpiDatatype<RealType>(), MPI_SUM, MPI_COMM_WORLD);
        return {hostPartial_[0], hostPartial_[1]};
    }

    const DomainT& domain_;
    const PeriodicMap<KeyType, RealType>& map_;
    size_t n_;
    int blockSize_;
    cstone::DeviceVector<uint8_t> isDof_;
    mutable Vector partial_;
    mutable std::vector<RealType> hostPartial_;
    long long numDofs_ = 0;
};

} // namespace fem
} // namespace mars
