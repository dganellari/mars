#pragma once

// The unknowns of a nodal field on a distributed mesh.
//
// Every rank stores one value per node slot: its owned nodes and ghost copies
// of nodes owned by other ranks. On a periodic box a point on a max face has
// more slots, on any ranks: the slaves of its master on the min faces. The
// unknown (DOF) of a point lives in ONE slot, the owned slot of its master.
// Two maps connect slots and DOFs:
//
//   prolong  (P)   copy each DOF value into every slot of its point
//   restrict (P^T) add the per-slot contributions of an element scatter into
//                  the DOF slot (the other slots end up holding nothing useful)
//
// A discrete operator is then restrict(A_local(prolong(x))) = P^T A P, with
// A_local the element scatter over this rank's own elements.
//
// Every slot finds its DOF from keys alone. The DOF key is the slot's SFC key,
// with each max-face axis moved to the min face on a periodic box
// (periodicMasterKey), and the DOF lives on the SFC owner of that key. So P and
// P^T are one exchange each, straight between a slot and the owner of its DOF,
// however a point is split over ranks and periodic faces.

#include "backend/distributed/unstructured/fem/mars_periodic_bc.hpp"

#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <iterator>
#include <type_traits>
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

// Up to Max node fields that share one exchange.
template<typename T>
struct DofFields
{
    static constexpr int Max = 8;
    T* f[Max];
    int count = 0;
};

// The DOF key of every slot (its own key, or its master's on a periodic max face),
// and whether the slot holds the DOF: owned and on no max face. mask == nullptr: no periodic box.
template<typename KeyType, typename RealType>
__global__ void dofKeyKernel(const KeyType* keys, const uint8_t* mask, const uint8_t* own, size_t n,
                             cstone::Box<RealType> box, unsigned minX, unsigned minY, unsigned minZ, KeyType* dofKey,
                             uint8_t* isDof)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    uint8_t m = mask ? mask[i] : 0;
    dofKey[i] = m ? periodicMasterKey<KeyType>(keys[i], m, box, minX, minY, minZ) : keys[i];
    isDof[i]  = (own[i] == 1 && m == 0) ? 1 : 0;
}

template<typename KeyType>
struct SingleRankOwner
{
    HOST_DEVICE_FUN int operator()(KeyType) const { return 0; }
};

// The rank holding the DOF of each slot, -1 for the DOF slots themselves.
template<typename KeyType, typename Owner>
__global__ void dofOwnerKernel(const KeyType* dofKey, const uint8_t* isDof, size_t n, Owner owner, int* rank)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    rank[i] = isDof[i] ? -1 : owner(dofKey[i]);
}

struct DofIsCopy
{
    HOST_DEVICE_FUN bool operator()(int rank) const { return rank >= 0; }
};

// The slot of the DOF with key query[k] on this rank, -1 if this rank does not hold it as a DOF.
template<typename KeyType>
__global__ void dofSlotKernel(const KeyType* keys, const uint8_t* isDof, size_t n, const KeyType* query,
                              size_t count, int* slot)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= count) return;
    const KeyType* pos = thrust::lower_bound(thrust::seq, keys, keys + n, query[k]);
    size_t j           = pos - keys;
    slot[k]            = (j < n && keys[j] == query[k] && isDof[j]) ? int(j) : -1;
}

// buf[i * count + c] = field c at slot ids[i]
template<typename T>
__global__ void dofPackKernel(const int* ids, size_t n, DofFields<T> f, T* buf)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    for (int c = 0; c < f.count; ++c)
        buf[i * f.count + c] = f.f[c][ids[i]];
}

// field c at slot ids[i] = (or +=) buf[i * count + c]; a DOF can receive from several copies, so += is atomic
template<typename T>
__global__ void dofUnpackKernel(const int* ids, size_t n, DofFields<T> f, const T* buf, bool add)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    for (int c = 0; c < f.count; ++c)
    {
        if (add) atomicAdd(&f.f[c][ids[i]], buf[i * f.count + c]);
        else f.f[c][ids[i]] = buf[i * f.count + c];
    }
}

// Copies whose DOF is on this rank: copy <- DOF, or DOF += copy
template<typename T>
__global__ void dofLocalKernel(const int* copy, const int* dof, size_t n, DofFields<T> f, bool add)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= n) return;
    for (int c = 0; c < f.count; ++c)
    {
        if (add) atomicAdd(&f.f[c][dof[k]], f.f[c][copy[k]]);
        else f.f[c][copy[k]] = f.f[c][dof[k]];
    }
}

struct DofNonNegative
{
    HOST_DEVICE_FUN bool operator()(int v) const { return v >= 0; }
};

// Length of the CSR row of each listed slot.
template<typename IndexType>
__global__ void dofRowLengthKernel(const IndexType* slots, size_t count, const IndexType* rowPtr, IndexType* length)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k < count) length[k] = rowPtr[slots[k] + 1] - rowPtr[slots[k]];
}

// The entries of each listed slot row from offset[k] on: the global id of the column slot's DOF
// and the value; with outRow, also the row the entries go to (target[k]).
template<typename T>
__global__ void dofRowEntriesKernel(const int* slots, const int* target, const long long* offset, size_t count,
                                    const int* rowPtr, const int* cols, const T* values, const long long* slotGid,
                                    int* outRow, long long* outCol, T* outVal)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= count) return;
    int s       = slots[k];
    long long o = offset[k];
    for (int j = rowPtr[s]; j < rowPtr[s + 1]; ++j, ++o)
    {
        if (outRow) outRow[o] = target[k];
        outCol[o] = slotGid[cols[j]];
        outVal[o] = values[j];
    }
}

// The row of each received entry: those of request k belong to DOF slot dofSlot[k].
template<typename IndexType>
__global__ void dofReceivedRowsKernel(const IndexType* dofSlot, const IndexType* dofRow, const IndexType* length,
                                      const long long* offset, size_t count, IndexType* outRow)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= count) return;
    int r = dofRow[dofSlot[k]];
    for (int j = 0; j < length[k]; ++j)
        outRow[offset[k] + j] = r;
}

template<typename KeyType, typename RealType, typename DomainT>
class DofSpace
{
public:
    using Vector = cstone::DeviceVector<RealType>;
    using Map    = PeriodicMap<KeyType, RealType>;
    using Fields = DofFields<RealType>;

    // periodic == nullptr: no periodic box. On several ranks the domain must use SFC node ownership.
    DofSpace(const DomainT& domain, const Map* periodic, int blockSize = 256)
        : n_(domain.getNodeCount())
        , blockSize_(blockSize)
        , isDof_(n_)
    {
        MPI_Comm_dup(MPI_COMM_WORLD, &comm_);
        build(domain, periodic);
    }
    DofSpace(const DofSpace&)            = delete;
    DofSpace& operator=(const DofSpace&) = delete;
    ~DofSpace() { MPI_Comm_free(&comm_); }

    // P: every copy takes the value of its DOF. The fields share one exchange.
    void prolong(Vector* const* fields, int count) const
    {
        Fields f = pointers(fields, count);
        local(f, false);
        exchange(f, false);
    }

    // P^T: every copy adds into its DOF. Must be applied to a scatter over this rank's own elements.
    void restrict(Vector* const* fields, int count) const
    {
        Fields f = pointers(fields, count);
        exchange(f, true);
        local(f, true);
    }

    void prolong(std::initializer_list<Vector*> fields) const { prolong(fields.begin(), int(fields.size())); }
    void restrict(std::initializer_list<Vector*> fields) const { restrict(fields.begin(), int(fields.size())); }
    void prolong(Vector& v) const { prolong({&v}); }
    void restrict(Vector& acc) const { restrict({&acc}); }

    long long numDofs() const { return numDofs_; }
    const uint8_t* isDof() const { return isDof_.data(); }

    // P^T A P for a matrix A over this rank's slots: the row of every copy is added into the
    // row of its DOF, on whatever rank owns it, in one exchange (row lengths, then entries),
    // the reverse direction of prolong. rowPtr/cols/values: A as a CSR over the local slots.
    // slotGid: the global id of each slot's DOF. dofRow: the local row of each DOF slot, -1
    // elsewhere. Out: the entries of this rank's DOF rows (local row, global column, value),
    // unsorted and with duplicates.
    template<typename T>
    void restrictMatrix(const int* rowPtr, const int* cols, const T* values, const long long* slotGid,
                        const int* dofRow, thrust::device_vector<int>& outRow, thrust::device_vector<long long>& outCol,
                        thrust::device_vector<T>& outVal) const
    {
        // Rows that stay: the DOF slots, then the copies of DOFs on this rank.
        const size_t numLocal = localCopy_.size();
        thrust::device_vector<int> stay(n_ + numLocal), target(n_ + numLocal);
        const size_t numDofs  = thrust::copy_if(thrust::device, thrust::counting_iterator<int>(0),
                                                thrust::counting_iterator<int>(int(n_)), dofRow, stay.begin(),
                                                DofNonNegative{}) -
                               stay.begin();
        thrust::gather(thrust::device, stay.begin(), stay.begin() + numDofs, dofRow, target.begin());
        if (numLocal > 0)
        {
            thrust::copy(thrust::device, localCopy_.data(), localCopy_.data() + numLocal, stay.begin() + numDofs);
            thrust::gather(thrust::device, localDof_.data(), localDof_.data() + numLocal, dofRow,
                           target.begin() + numDofs);
        }
        const size_t numStay = numDofs + numLocal;

        // Rows that leave: the copies of DOFs on other ranks, grouped by peer.
        const size_t numSend = recvCopy_.size(), numRecv = sendDof_.size();
        thrust::device_vector<int> stayLength(numStay), sendLength(numSend), recvLength(numRecv);
        if (numStay > 0)
            dofRowLengthKernel<int><<<int((numStay + 255) / 256), 256>>>(thrust::raw_pointer_cast(stay.data()), numStay,
                                                                   rowPtr, thrust::raw_pointer_cast(stayLength.data()));
        if (numSend > 0)
            dofRowLengthKernel<int><<<int((numSend + 255) / 256), 256>>>(recvCopy_.data(), numSend, rowPtr,
                                                                   thrust::raw_pointer_cast(sendLength.data()));
        cudaCheckError();
        thrust::device_vector<long long> stayOffset(numStay), sendOffset(numSend);
        thrust::exclusive_scan(thrust::device, stayLength.begin(), stayLength.end(), stayOffset.begin(), 0LL);
        thrust::exclusive_scan(thrust::device, sendLength.begin(), sendLength.end(), sendOffset.begin(), 0LL);
        const long long stayTotal = thrust::reduce(thrust::device, stayLength.begin(), stayLength.end(), 0LL);
        const long long sendTotal = thrust::reduce(thrust::device, sendLength.begin(), sendLength.end(), 0LL);

        // Row lengths to the owners, so each side knows the entry counts per peer.
        cudaDeviceSynchronize();
        std::vector<MPI_Request> requests;
        requests.reserve(4 * peers_.size());
        for (size_t p = 0; p < peers_.size(); ++p)
        {
            int count = sendOffsets_[p + 1] - sendOffsets_[p];
            if (count == 0) continue;
            requests.emplace_back();
            MPI_Irecv(thrust::raw_pointer_cast(recvLength.data()) + sendOffsets_[p], count, MPI_INT, peers_[p], 0x4e60,
                      comm_, &requests.back());
        }
        for (size_t p = 0; p < peers_.size(); ++p)
        {
            int count = recvOffsets_[p + 1] - recvOffsets_[p];
            if (count == 0) continue;
            requests.emplace_back();
            MPI_Isend(thrust::raw_pointer_cast(sendLength.data()) + recvOffsets_[p], count, MPI_INT, peers_[p], 0x4e60,
                      comm_, &requests.back());
        }
        MPI_Waitall(int(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        requests.clear();
        thrust::device_vector<long long> recvOffset(numRecv);
        thrust::exclusive_scan(thrust::device, recvLength.begin(), recvLength.end(), recvOffset.begin(), 0LL);
        const long long recvTotal = thrust::reduce(thrust::device, recvLength.begin(), recvLength.end(), 0LL);

        // Entries per peer, on the host: the row lengths summed over each peer's rows.
        std::vector<int> hostSend(numSend), hostRecv(numRecv);
        thrust::copy(sendLength.begin(), sendLength.end(), hostSend.begin());
        thrust::copy(recvLength.begin(), recvLength.end(), hostRecv.begin());
        auto entriesOf = [](const std::vector<int>& length, int begin, int end) {
            long long sum = 0;
            for (int k = begin; k < end; ++k)
                sum += length[k];
            return sum;
        };

        outRow.resize(stayTotal + recvTotal);
        outCol.resize(stayTotal + recvTotal);
        outVal.resize(stayTotal + recvTotal);
        thrust::device_vector<long long> sendCol(sendTotal);
        thrust::device_vector<T> sendVal(sendTotal);
        if (numStay > 0)
            dofRowEntriesKernel<T><<<int((numStay + 255) / 256), 256>>>(
                thrust::raw_pointer_cast(stay.data()), thrust::raw_pointer_cast(target.data()),
                thrust::raw_pointer_cast(stayOffset.data()), numStay, rowPtr, cols, values, slotGid,
                thrust::raw_pointer_cast(outRow.data()), thrust::raw_pointer_cast(outCol.data()),
                thrust::raw_pointer_cast(outVal.data()));
        if (numSend > 0)
            dofRowEntriesKernel<T><<<int((numSend + 255) / 256), 256>>>(
                recvCopy_.data(), nullptr, thrust::raw_pointer_cast(sendOffset.data()), numSend, rowPtr, cols, values,
                slotGid, nullptr, thrust::raw_pointer_cast(sendCol.data()), thrust::raw_pointer_cast(sendVal.data()));
        cudaCheckError();
        cudaDeviceSynchronize();

        // The entries themselves: received straight into the output, after this rank's own.
        long long* recvColPtr = thrust::raw_pointer_cast(outCol.data()) + stayTotal;
        T* recvValPtr         = thrust::raw_pointer_cast(outVal.data()) + stayTotal;
        long long recvAt = 0, sendAt = 0;
        for (size_t p = 0; p < peers_.size(); ++p)
        {
            long long count = entriesOf(hostRecv, sendOffsets_[p], sendOffsets_[p + 1]);
            if (count > 0)
            {
                requests.emplace_back();
                MPI_Irecv(recvColPtr + recvAt, int(count), MPI_LONG_LONG, peers_[p], 0x4e61, comm_, &requests.back());
                requests.emplace_back();
                MPI_Irecv(recvValPtr + recvAt, int(count), mpiDatatype<T>(), peers_[p], 0x4e62, comm_,
                          &requests.back());
            }
            recvAt += count;
        }
        for (size_t p = 0; p < peers_.size(); ++p)
        {
            long long count = entriesOf(hostSend, recvOffsets_[p], recvOffsets_[p + 1]);
            if (count > 0)
            {
                requests.emplace_back();
                MPI_Isend(thrust::raw_pointer_cast(sendCol.data()) + sendAt, int(count), MPI_LONG_LONG, peers_[p],
                          0x4e61, comm_, &requests.back());
                requests.emplace_back();
                MPI_Isend(thrust::raw_pointer_cast(sendVal.data()) + sendAt, int(count), mpiDatatype<T>(), peers_[p],
                          0x4e62, comm_, &requests.back());
            }
            sendAt += count;
        }
        MPI_Waitall(int(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        if (numRecv > 0)
            dofReceivedRowsKernel<int><<<int((numRecv + 255) / 256), 256>>>(
                sendDof_.data(), dofRow, thrust::raw_pointer_cast(recvLength.data()),
                thrust::raw_pointer_cast(recvOffset.data()), numRecv, thrust::raw_pointer_cast(outRow.data()) + stayTotal);
        cudaCheckError();
    }

    // The exchanges of prolong and restrict on this rank since the last reset: how many,
    // their time from pack to the end of the MPI wait, and the MPI part of it. With
    // MARS_EXCHANGE_PROFILE set, a barrier before the MPI part splits it: waitMs is the wait for the
    // last rank to arrive (skew from the work before), mpiMs the transfer after it.
    struct ExchangeStats
    {
        long count    = 0;
        double ms     = 0;
        double mpiMs  = 0;
        double waitMs = 0;
    };
    const ExchangeStats& exchangeStats() const { return stats_; }
    void resetExchangeStats() { stats_ = {}; }
    bool profiling() const { return profile_; }

    // This rank's peers and the values it sends per field in a forward exchange, split by whether
    // the peer runs on the same node. Collective.
    struct ExchangeShape
    {
        int peers = 0, peersOnNode = 0;
        long long sendValues = 0, sendValuesOffNode = 0;
    };
    ExchangeShape exchangeShape() const
    {
        int rank = 0, size = 1;
        MPI_Comm_rank(comm_, &rank);
        MPI_Comm_size(comm_, &size);
        MPI_Comm node;
        MPI_Comm_split_type(comm_, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node);
        int leader = rank;
        MPI_Allreduce(MPI_IN_PLACE, &leader, 1, MPI_INT, MPI_MIN, node);
        MPI_Comm_free(&node);
        std::vector<int> nodeOf(size);
        MPI_Allgather(&leader, 1, MPI_INT, nodeOf.data(), 1, MPI_INT, comm_);
        ExchangeShape s;
        for (size_t p = 0; p < peers_.size(); ++p)
        {
            const long long sent = sendOffsets_[p + 1] - sendOffsets_[p];
            const bool onNode    = nodeOf[peers_[p]] == nodeOf[rank];
            s.peers += 1;
            s.peersOnNode += onNode ? 1 : 0;
            s.sendValues += sent;
            if (!onNode) s.sendValuesOffNode += sent;
        }
        return s;
    }

    // reps exchanges of these fields back to back, each after a barrier, so no skew from earlier work
    // is in them: the cost of the exchange pattern alone. Per exchange, this rank's time in ms from
    // the barrier to the end of the unpack. Collective; the exchange statistics are left unchanged.
    std::vector<double> benchmark(Vector* const* fields, int count, int reps, bool reverse) const
    {
        const ExchangeStats saved = stats_;
        Fields f                  = pointers(fields, count);
        std::vector<double> t(reps);
        for (int r = 0; r < reps; ++r)
        {
            MPI_Barrier(comm_);
            const double start = MPI_Wtime();
            exchange(f, reverse);
            cudaDeviceSynchronize();
            t[r] = 1e3 * (MPI_Wtime() - start);
        }
        stats_ = saved;
        return t;
    }

private:
    int grid() const { return std::max(1, int((n_ + blockSize_ - 1) / blockSize_)); }

    static Fields pointers(Vector* const* fields, int count)
    {
        if (count > Fields::Max)
        {
            std::fprintf(stderr, "DofSpace: %d fields in one exchange, at most %d\n", count, Fields::Max);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        Fields f;
        f.count = count;
        for (int c = 0; c < count; ++c)
            f.f[c] = fields[c]->data();
        return f;
    }

    // Classifies the slots, then asks the owner of every remote DOF key which slot holds it
    // (one sparse exchange): peers_, the DOF slots sent to each peer and the copies received.
    void build(const DomainT& domain, const Map* periodic)
    {
        const int numRanks  = domain.numRanks();
        const int rank      = domain.rank();
        const KeyType* keys = domain.getLocalToGlobalSfcMap().data();
        const uint8_t* own  = domain.getNodeOwnershipMap().data();
        const auto& box     = domain.getBoundingBox();
        if (numRanks > 1 && !domain.sfcOwnership())
        {
            if (rank == 0)
                std::fprintf(stderr, "DofSpace: several ranks need SFC node ownership; unset MARS_OWNERSHIP\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        const uint8_t* mask = nullptr;
        unsigned minX = 0, minY = 0, minZ = 0;
        if (periodic)
        {
            if (periodic->d_periodicMask.size() != n_)
            {
                std::fprintf(stderr, "DofSpace: the periodic map belongs to another mesh\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            mask = periodic->d_periodicMask.data();
            auto [x0, y0, z0] = cstone::decodeSfc(
                cstone::sfc3D<cstone::SfcKind<KeyType>>(periodic->xmin, periodic->ymin, periodic->zmin, box),
                box.getBoxDimBits(cstone::maxTreeLevel<KeyType>{}));
            minX = x0;
            minY = y0;
            minZ = z0;
        }

        cstone::DeviceVector<KeyType> dofKey(n_);
        thrust::device_vector<int> owner(n_);
        int* ownerPtr = thrust::raw_pointer_cast(owner.data());
        if (n_ > 0)
        {
            dofKeyKernel<KeyType, RealType><<<grid(), blockSize_>>>(keys, mask, own, n_, box, minX, minY, minZ,
                                                                     dofKey.data(), isDof_.data());
            if (numRanks > 1)
                dofOwnerKernel<<<grid(), blockSize_>>>(dofKey.data(), isDof_.data(), n_, domain.sfcNodeOwner(),
                                                       ownerPtr);
            else
                dofOwnerKernel<<<grid(), blockSize_>>>(dofKey.data(), isDof_.data(), n_, SingleRankOwner<KeyType>{},
                                                       ownerPtr);
        }
        cudaCheckError();
        long long dofs = thrust::count(thrust::device, isDof_.data(), isDof_.data() + n_, uint8_t(1));
        MPI_Allreduce(&dofs, &numDofs_, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        // The copies grouped by the rank of their DOF; within a rank, slot ids ascend.
        size_t numCopies = thrust::count_if(thrust::device, owner.begin(), owner.end(), DofIsCopy{});
        thrust::device_vector<int> copySlot(numCopies), copyOwner(numCopies);
        thrust::device_vector<KeyType> copyKey(numCopies);
        thrust::copy_if(thrust::device, thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(int(n_)),
                        owner.begin(), copySlot.begin(), DofIsCopy{});
        thrust::gather(thrust::device, copySlot.begin(), copySlot.end(), owner.begin(), copyOwner.begin());
        thrust::stable_sort_by_key(thrust::device, copyOwner.begin(), copyOwner.end(), copySlot.begin());
        thrust::gather(thrust::device, copySlot.begin(), copySlot.end(), thrust::device_pointer_cast(dofKey.data()),
                       copyKey.begin());

        thrust::device_vector<int> d_owners(numCopies), d_counts(numCopies);
        auto ends        = thrust::reduce_by_key(thrust::device, copyOwner.begin(), copyOwner.end(),
                                                 thrust::constant_iterator<int>(1), d_owners.begin(), d_counts.begin());
        size_t numOwners = ends.first - d_owners.begin();
        std::vector<int> owners(numOwners), counts(numOwners);
        thrust::copy(d_owners.begin(), d_owners.begin() + numOwners, owners.begin());
        thrust::copy(d_counts.begin(), d_counts.begin() + numOwners, counts.begin());

        // This rank's segment is the local copies; the others go to their owners in rank order.
        size_t localBegin = 0, localCount = 0, offset = 0;
        std::vector<int> dests, destCounts;
        for (size_t o = 0; o < numOwners; ++o)
        {
            if (owners[o] == rank)
            {
                localBegin = offset;
                localCount = counts[o];
            }
            else
            {
                dests.push_back(owners[o]);
                destCounts.push_back(counts[o]);
            }
            offset += counts[o];
        }

        const int* copySlotPtr   = thrust::raw_pointer_cast(copySlot.data());
        const KeyType* copyKeyPtr = thrust::raw_pointer_cast(copyKey.data());
        localCopy_.resize(localCount);
        localDof_.resize(localCount);
        if (localCount > 0)
        {
            cudaMemcpy(localCopy_.data(), copySlotPtr + localBegin, localCount * sizeof(int), cudaMemcpyDeviceToDevice);
            dofSlotKernel<KeyType><<<int((localCount + 255) / 256), 256>>>(keys, isDof_.data(), n_,
                                                                          copyKeyPtr + localBegin, localCount,
                                                                          localDof_.data());
        }

        size_t numRemote = numCopies - localCount, tail = numCopies - localBegin - localCount;
        recvCopy_.resize(numRemote);
        thrust::device_vector<KeyType> remoteKey(numRemote);
        KeyType* remoteKeyPtr = thrust::raw_pointer_cast(remoteKey.data());
        if (localBegin > 0)
        {
            cudaMemcpy(recvCopy_.data(), copySlotPtr, localBegin * sizeof(int), cudaMemcpyDeviceToDevice);
            cudaMemcpy(remoteKeyPtr, copyKeyPtr, localBegin * sizeof(KeyType), cudaMemcpyDeviceToDevice);
        }
        if (tail > 0)
        {
            cudaMemcpy(recvCopy_.data() + localBegin, copySlotPtr + localBegin + localCount, tail * sizeof(int),
                       cudaMemcpyDeviceToDevice);
            cudaMemcpy(remoteKeyPtr + localBegin, copyKeyPtr + localBegin + localCount, tail * sizeof(KeyType),
                       cudaMemcpyDeviceToDevice);
        }

        cstone::DeviceVector<KeyType> requested;
        constexpr int tagDofRequests = 0x4d43;
        auto requests = sparseExchange<KeyType>(dests, destCounts, remoteKeyPtr, requested, tagDofRequests, comm_);
        size_t numRequested = requested.size();
        thrust::device_vector<int> requestedDof(numRequested);
        int* requestedDofPtr = thrust::raw_pointer_cast(requestedDof.data());
        if (numRequested > 0)
            dofSlotKernel<KeyType><<<int((numRequested + 255) / 256), 256>>>(keys, isDof_.data(), n_, requested.data(),
                                                                            numRequested, requestedDofPtr);
        cudaCheckError();

        // A key this rank does not hold as a DOF means the ranks disagree on ownership, or a
        // periodic slave has no matching master (box bounds, faceEps).
        long long missing = thrust::count(thrust::device, localDof_.data(), localDof_.data() + localCount, -1) +
                            thrust::count(thrust::device, requestedDof.begin(), requestedDof.end(), -1);
        MPI_Allreduce(MPI_IN_PLACE, &missing, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (missing > 0)
        {
            if (rank == 0)
                std::fprintf(stderr, "DofSpace: %lld node copies have no DOF on the rank that owns its key\n", missing);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Peers: ranks we receive DOFs from (dests) and ranks that asked us (requests), both ascending.
        std::vector<int> requesters;
        for (const auto& r : requests)
            requesters.push_back(r.source);
        peers_.clear();
        std::set_union(dests.begin(), dests.end(), requesters.begin(), requesters.end(), std::back_inserter(peers_));
        sendOffsets_.assign(1, 0);
        recvOffsets_.assign(1, 0);
        sendDof_.resize(numRequested);
        size_t destIdx = 0, requestIdx = 0;
        for (int peer : peers_)
        {
            int sendCount = 0, recvCount = 0;
            if (requestIdx < requests.size() && requests[requestIdx].source == peer)
            {
                const auto& r = requests[requestIdx++];
                sendCount     = r.count;
                cudaMemcpy(sendDof_.data() + sendOffsets_.back(), requestedDofPtr + r.offset, sendCount * sizeof(int),
                           cudaMemcpyDeviceToDevice);
            }
            if (destIdx < dests.size() && dests[destIdx] == peer) recvCount = destCounts[destIdx++];
            sendOffsets_.push_back(sendOffsets_.back() + sendCount);
            recvOffsets_.push_back(recvOffsets_.back() + recvCount);
        }
    }

    void local(Fields f, bool add) const
    {
        size_t n = localCopy_.size();
        if (n == 0 || f.count == 0) return;
        dofLocalKernel<RealType><<<int((n + 255) / 256), 256>>>(localCopy_.data(), localDof_.data(), n, f, add);
        cudaCheckError();
    }

    // Forward: DOF values to the copies on other ranks. Reverse: copies added into their DOFs.
    void exchange(Fields f, bool reverse) const
    {
        if (peers_.empty() || f.count == 0) return;
        const int k               = f.count;
        const auto& packIds       = reverse ? recvCopy_ : sendDof_;
        const auto& packOffsets   = reverse ? recvOffsets_ : sendOffsets_;
        const auto& unpackIds     = reverse ? sendDof_ : recvCopy_;
        const auto& unpackOffsets = reverse ? sendOffsets_ : recvOffsets_;
        size_t packTotal = size_t(packOffsets.back()), unpackTotal = size_t(unpackOffsets.back());
        if (sendBuf_.size() < packTotal * k) sendBuf_.resize(packTotal * k);
        if (recvBuf_.size() < unpackTotal * k) recvBuf_.resize(unpackTotal * k);

        // The sync also finishes the previous unpack before its receive buffer is reused.
        double start = MPI_Wtime();
        if (packTotal > 0)
            dofPackKernel<RealType><<<int((packTotal + 255) / 256), 256>>>(packIds.data(), packTotal, f,
                                                                           sendBuf_.data());
        cudaDeviceSynchronize();
        double packed = MPI_Wtime();
        if (profile_) MPI_Barrier(comm_);
        double synced = MPI_Wtime();

        auto type     = mpiDatatype<RealType>();
        const int tag = reverse ? 0x4e51 : 0x4e50;
        std::vector<MPI_Request> requests;
        requests.reserve(2 * peers_.size());
        for (size_t p = 0; p < peers_.size(); ++p)
        {
            int count = unpackOffsets[p + 1] - unpackOffsets[p];
            if (count == 0) continue;
            requests.emplace_back();
            MPI_Irecv(recvBuf_.data() + size_t(unpackOffsets[p]) * k, count * k, type, peers_[p], tag, comm_,
                      &requests.back());
        }
        for (size_t p = 0; p < peers_.size(); ++p)
        {
            int count = packOffsets[p + 1] - packOffsets[p];
            if (count == 0) continue;
            requests.emplace_back();
            MPI_Isend(sendBuf_.data() + size_t(packOffsets[p]) * k, count * k, type, peers_[p], tag, comm_,
                      &requests.back());
        }
        MPI_Waitall(int(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
        double done = MPI_Wtime();
        stats_.count += 1;
        stats_.ms += 1e3 * (done - start);
        stats_.waitMs += 1e3 * (synced - packed);
        stats_.mpiMs += 1e3 * (done - synced);

        if (unpackTotal > 0)
            dofUnpackKernel<RealType><<<int((unpackTotal + 255) / 256), 256>>>(unpackIds.data(), unpackTotal, f,
                                                                               recvBuf_.data(), reverse);
        cudaCheckError();
    }

    size_t n_;
    int blockSize_;
    MPI_Comm comm_ = MPI_COMM_NULL;
    cstone::DeviceVector<uint8_t> isDof_;
    long long numDofs_ = 0;
    cstone::DeviceVector<int> localCopy_, localDof_; // copies whose DOF is on this rank
    std::vector<int> peers_, sendOffsets_, recvOffsets_;
    cstone::DeviceVector<int> sendDof_, recvCopy_;   // per peer: DOF slots sent, copies received
    mutable Vector sendBuf_, recvBuf_;
    mutable ExchangeStats stats_;
    bool profile_ = std::getenv("MARS_EXCHANGE_PROFILE") != nullptr;
};

} // namespace fem
} // namespace mars
