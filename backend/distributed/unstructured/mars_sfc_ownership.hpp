#pragma once

// SFC node ownership for multi-rank meshes.
//
// A node belongs to the rank whose cornerstone SFC range contains the node's position. Every rank can evaluate
// this from the replicated assignment (numRanks + 1 keys), so ownership needs no communication and every node has
// exactly one owner. The owner may hold no element of its own that touches the node, and cornerstone's
// distance-based halo search does not guarantee that it receives all of them either. So after the element sync each
// rank sends every local element to the owners of its corners that do not already get it as a halo, and those
// owners add it as a regular cornerstone halo (Domain::addHalos). Then each owner holds the complete element star of
// its nodes and can assemble their rows.
//
// Host and device code: the same functions run in the GPU kernels and in the host MPI test.

#include <mpi.h>

#include <algorithm>
#include <span>
#include <utility>
#include <vector>

#include "cstone/cuda/annotation.hpp"
#include "cstone/domain/index_ranges.hpp"
#include "cstone/primitives/mpi_wrappers.hpp"
#include "cstone/primitives/stl.hpp"
#include "cstone/sfc/box.hpp"
#include "cstone/sfc/sfc.hpp"

namespace mars
{

template<class KeyType, class RealType>
struct SfcNodeOwner
{
    cstone::Box<RealType> meshBox; // box of the node keys
    cstone::Box<RealType> sfcBox;  // cornerstone box of the element decomposition
    const KeyType* rankBounds;     // cornerstone assignment, numRanks + 1 keys
    int numRanks;

    HOST_DEVICE_FUN int operator()(KeyType nodeKey) const
    {
        // Mixed-dimension keys: each axis has its own number of bits, set by the box they were encoded with
        const auto bits   = meshBox.getBoxDimBits(cstone::maxTreeLevel<KeyType>{});
        auto [ix, iy, iz] = cstone::decodeSfc(cstone::sfcKey(nodeKey), bits);

        // Cornerstone fits its box to the element representative corners, so nodes can lie outside it
        RealType x = meshBox.xmin() + ix * (RealType(1) / ((1u << bits[0]) - 1)) * (meshBox.xmax() - meshBox.xmin());
        RealType y = meshBox.ymin() + iy * (RealType(1) / ((1u << bits[1]) - 1)) * (meshBox.ymax() - meshBox.ymin());
        RealType z = meshBox.zmin() + iz * (RealType(1) / ((1u << bits[2]) - 1)) * (meshBox.zmax() - meshBox.zmin());
        x          = ::stl::min(::stl::max(x, sfcBox.xmin()), sfcBox.xmax());
        y          = ::stl::min(::stl::max(y, sfcBox.ymin()), sfcBox.ymax());
        z          = ::stl::min(::stl::max(z, sfcBox.zmin()), sfcBox.zmax());

        KeyType key = cstone::sfc3D<cstone::SfcKind<KeyType>>(x, y, z, sfcBox);
        return int(::stl::upper_bound(rankBounds, rankBounds + numRanks + 1, key) - rankBounds) - 1;
    }
};

// Local element index ranges that cornerstone already sends to each rank as halos, CSR over all ranks
struct HaloSendRanges
{
    const int* rankOffsets;
    const cstone::LocalIndex* start; // sorted within each rank
    const cstone::LocalIndex* end;

    HOST_DEVICE_FUN bool contains(int rank, cstone::LocalIndex i) const
    {
        int first = rankOffsets[rank];
        int lo = first, hi = rankOffsets[rank + 1];
        while (lo < hi)
        {
            int mid = (lo + hi) / 2;
            if (start[mid] <= i) { lo = mid + 1; }
            else { hi = mid; }
        }
        return lo > first && i < end[lo - 1];
    }
};

struct HaloSendRangesHost
{
    std::vector<int> rankOffsets;
    std::vector<cstone::LocalIndex> start, end;
};

inline HaloSendRangesHost flattenHaloSendRanges(const cstone::SendList& outgoing, int numRanks)
{
    HaloSendRangesHost r;
    r.rankOffsets.assign(numRanks + 1, 0);
    for (int rank = 0; rank < numRanks; ++rank)
    {
        std::vector<std::pair<cstone::LocalIndex, cstone::LocalIndex>> ranges;
        if (size_t(rank) < outgoing.size())
        {
            const auto& m = outgoing[rank];
            for (size_t i = 0; i < m.nRanges(); ++i)
            {
                if (m.count(i) > 0) { ranges.emplace_back(m.rangeStart(i), m.rangeEnd(i)); }
            }
        }
        std::sort(ranges.begin(), ranges.end());
        for (auto [a, b] : ranges)
        {
            r.start.push_back(a);
            r.end.push_back(b);
        }
        r.rankOffsets[rank + 1] = int(r.start.size());
    }
    return r;
}

// Ranks other than myRank that own a corner of local element `elem` and do not receive it as a cornerstone halo.
// Returns their number; dests needs room for NumCorners entries.
template<int NumCorners, class KeyType, class RealType>
HOST_DEVICE_FUN int missingStarRanks(const KeyType* corners,
                                     cstone::LocalIndex elem,
                                     const SfcNodeOwner<KeyType, RealType>& owner,
                                     const HaloSendRanges& sent,
                                     int myRank,
                                     int* dests)
{
    int count = 0;
    for (int c = 0; c < NumCorners; ++c)
    {
        int rank  = owner(corners[c]);
        bool skip = rank == myRank || sent.contains(rank, elem);
        for (int j = 0; j < count; ++j)
        {
            skip = skip || dests[j] == rank;
        }
        if (!skip) { dests[count++] = rank; }
    }
    return count;
}

// One message received by sparseExchange: `count` items at `offset` in the receive buffer
struct ReceivedMessage
{
    int source;
    size_t offset;
    int count;
};

// Sends messages to a sparse set of ranks without knowing which ranks send to us (NBX: Hoefler, Siebert,
// Lumsdaine, PPoPP 2010). A synchronous send completes only once it is matched, so when the non-blocking barrier
// completes, every message addressed to this rank has arrived; there is no collective over all ranks otherwise.
// sendBuf holds the messages back to back in the order of dests. With CUDA-aware MPI, sendBuf and recv may be
// device memory. Messages are returned sorted by source; the buffer keeps arrival order. Call sites need distinct
// tags unless a collective separates them.
template<class T, class Vector>
std::vector<ReceivedMessage> sparseExchange(std::span<const int> dests,
                                            std::span<const int> counts,
                                            const T* sendBuf,
                                            Vector& recv,
                                            int tag,
                                            MPI_Comm comm)
{
    std::vector<MPI_Request> sendRequests(dests.size());
    size_t sendOffset = 0;
    for (size_t i = 0; i < dests.size(); ++i)
    {
        MPI_Issend(sendBuf + sendOffset, counts[i], ::MpiType<T>{}, dests[i], tag, comm, &sendRequests[i]);
        sendOffset += counts[i];
    }

    std::vector<ReceivedMessage> arrived;
    recv.resize(0);
    MPI_Request barrier = MPI_REQUEST_NULL;
    bool barrierActive  = false;
    int done            = 0;
    while (!done)
    {
        int flag = 0;
        MPI_Status status;
        MPI_Iprobe(MPI_ANY_SOURCE, tag, comm, &flag, &status);
        if (flag)
        {
            int count = 0;
            MPI_Get_count(&status, ::MpiType<T>{}, &count);
            size_t offset = recv.size();
            recv.resize(offset + count);
            MPI_Recv(recv.data() + offset, count, ::MpiType<T>{}, status.MPI_SOURCE, tag, comm,
                     MPI_STATUS_IGNORE);
            arrived.push_back({status.MPI_SOURCE, offset, count});
        }

        if (barrierActive) { MPI_Test(&barrier, &done, MPI_STATUS_IGNORE); }
        else
        {
            int allSent = 0;
            MPI_Testall(int(sendRequests.size()), sendRequests.data(), &allSent, MPI_STATUSES_IGNORE);
            if (allSent)
            {
                MPI_Ibarrier(comm, &barrier);
                barrierActive = true;
            }
        }
    }

    std::sort(arrived.begin(), arrived.end(), [](const auto& a, const auto& b) { return a.source < b.source; });
    return arrived;
}

} // namespace mars
