#pragma once

// Hanging nodes of a non-conforming mesh, and the regular nodes each one takes its value from.
//
// Refining an element but not its neighbour leaves nodes in the middle of the neighbour's edges and, on
// hexahedra, of its faces. A continuous P1/Q1 field is linear along that edge and bilinear on that face, so its
// value at such a node is the mean of the edge ends (face corners): the node is not an unknown. When an edge end is
// itself hanging on a coarser element, its own constraint is substituted, down to regular nodes.
//
// Detection uses active split-edge/face records emitted by the refiner and matched to current coarse
// entities after redistribution. Key equality identifies nodes; parent incidence establishes the constraint.
// Each rank sends these records to the SFC node owner, then queries the nodes used by its owned elements.

#include "../domain.hpp"
#include "../amr/mars_amr_sparse_plan.hpp"
#include "../mars_sfc_ownership.hpp"
#include "../utils/mars_gpu_reduce.hpp"

#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/gather.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

#include <cstdio>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <mpi.h>

namespace mars
{
namespace fem
{

// The regular nodes behind one hanging node. A face centre whose four corners hang on edges of an element two
// levels coarser has 8; 16 leaves room for one more level.
constexpr int MaxHangingLinks = 16;

// A constraint as it travels between ranks: the node key, how many nodes it hangs on (2 or 4; 0: a regular node)
// and their keys in ascending order, unused ones ~0.
constexpr int HangingRecordSize = 6;

template<typename KeyType, int N>
struct HangingCorners
{
    const KeyType* c[N];
};

// Of records sorted by node key: keep[r] for the first of each key, conflict[r] if it differs from that first one
template<typename KeyType>
__global__ void hangingFirstOfKeyKernel(const KeyType* records, size_t n, uint8_t* keep, uint8_t* conflict)
{
    const size_t r = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (r >= n) return;
    const KeyType* cur = records + size_t(HangingRecordSize) * r;
    size_t first       = r;
    while (first > 0 && records[size_t(HangingRecordSize) * (first - 1)] == cur[0])
        --first;
    keep[r]                = first == r;
    const KeyType* ref     = records + size_t(HangingRecordSize) * first;
    bool differs           = false;
    for (int v = 1; v < HangingRecordSize; ++v)
        differs |= cur[v] != ref[v];
    conflict[r] = differs;
}

template<typename KeyType, typename Owner>
__global__ void hangingOwnerKernel(const KeyType* keys, size_t stride, size_t n, Owner owner, int* rank)
{
    const size_t r = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (r < n) rank[r] = owner(keys[stride * r]);
}

// The record of each query key from a table of records sorted by node key; a regular node if absent
template<typename KeyType>
__global__ void hangingLookupKernel(const KeyType* query, size_t numQuery, const KeyType* tableKeys,
                                    const KeyType* tableRecords, size_t tableSize, KeyType* answer)
{
    const size_t q = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (q >= numQuery) return;
    const KeyType* pos = thrust::lower_bound(thrust::seq, tableKeys, tableKeys + tableSize, query[q]);
    const size_t j     = pos - tableKeys;
    KeyType* out       = answer + size_t(HangingRecordSize) * q;
    if (j < tableSize && tableKeys[j] == query[q])
    {
        for (int v = 0; v < HangingRecordSize; ++v)
            out[v] = tableRecords[size_t(HangingRecordSize) * j + v];
    }
    else
    {
        out[0] = query[q];
        out[1] = 0;
        for (int v = 2; v < HangingRecordSize; ++v)
            out[v] = ~KeyType(0);
    }
}

template<typename KeyType>
struct HangingKeyOf
{
    const KeyType* records;
    __host__ __device__ KeyType operator()(size_t r) const { return records[size_t(HangingRecordSize) * r]; }
};

struct HangingIsSet
{
    __host__ __device__ bool operator()(uint8_t v) const { return v != 0; }
};

// One flattening round over the hanging slots: every link to a node that hangs itself is replaced by that node's
// links, scaled by its weight; equal keys merge. status: records of every key the links name, sorted by key.
// overflow: more than MaxHangingLinks regular nodes behind one node.
template<typename KeyType, typename RealType>
__global__ void hangingFlattenKernel(size_t numHanging, int* numLinks, KeyType* linkKey, RealType* linkWeight,
                                     const KeyType* statusKeys, const KeyType* statusRecords, size_t statusSize,
                                     int* changed, int* overflow)
{
    const size_t h = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (h >= numHanging) return;
    KeyType keys[MaxHangingLinks];
    RealType weights[MaxHangingLinks];
    int count    = 0;
    bool expand  = false;
    bool tooMany = false;
    auto add     = [&](KeyType k, RealType w)
    {
        for (int i = 0; i < count; ++i)
            if (keys[i] == k)
            {
                weights[i] += w;
                return;
            }
        if (count == MaxHangingLinks)
        {
            tooMany = true;
            return;
        }
        keys[count]      = k;
        weights[count++] = w;
    };
    KeyType* myKey    = linkKey + size_t(MaxHangingLinks) * h;
    RealType* myW     = linkWeight + size_t(MaxHangingLinks) * h;
    for (int l = 0; l < numLinks[h]; ++l)
    {
        const KeyType* pos = thrust::lower_bound(thrust::seq, statusKeys, statusKeys + statusSize, myKey[l]);
        const size_t j     = pos - statusKeys;
        const KeyType* rec = (j < statusSize && statusKeys[j] == myKey[l]) ? statusRecords + HangingRecordSize * j : nullptr;
        const int cnt      = rec ? int(rec[1]) : 0;
        if (cnt == 0)
        {
            add(myKey[l], myW[l]);
            continue;
        }
        expand = true;
        for (int v = 0; v < cnt; ++v)
            add(rec[2 + v], myW[l] / RealType(cnt));
    }
    for (int i = 0; i < count; ++i)
    {
        myKey[i] = keys[i];
        myW[i]   = weights[i];
    }
    numLinks[h] = count;
    if (expand) atomicExch(changed, 1);
    if (tooMany) atomicExch(overflow, 1);
}

template<typename KeyType>
struct HangingMarkCorner
{
    const KeyType* corner;
    uint8_t* used;
    __host__ __device__ void operator()(size_t e) const { used[corner[e]] = 1; }
};

struct HangingIsOwned
{
    __host__ __device__ uint8_t operator()(uint8_t ownership) const { return ownership == 1; }
};

template<typename KeyType>
struct HangingRecordHangs
{
    const KeyType* records;
    __host__ __device__ uint8_t operator()(size_t r) const { return records[size_t(HangingRecordSize) * r + 1] != 0; }
};

template<typename KeyType>
struct HangingIsKey
{
    __host__ __device__ bool operator()(KeyType k) const { return k != ~KeyType(0); }
};

// out[i] = items[order[i]], width keys each (scatter == false); out[order[i]] = items[i] (scatter == true)
template<typename KeyType>
struct HangingPermute
{
    const KeyType* items;
    const int* order;
    KeyType* out;
    int width;
    bool scatter;
    __host__ __device__ void operator()(size_t i) const
    {
        const size_t from = scatter ? i : size_t(order[i]);
        const size_t to   = scatter ? size_t(order[i]) : i;
        for (int v = 0; v < width; ++v)
            out[size_t(width) * to + v] = items[size_t(width) * from + v];
    }
};

template<typename KeyType, typename RealType>
struct HangingFirstLinks
{
    const int* which;
    const KeyType* records;
    int* numLinks;
    KeyType* keys;
    RealType* weights;
    __host__ __device__ void operator()(size_t h) const
    {
        const KeyType* r = records + size_t(HangingRecordSize) * which[h];
        const int cnt    = int(r[1]);
        numLinks[h]      = cnt;
        for (int v = 0; v < cnt; ++v)
        {
            keys[size_t(MaxHangingLinks) * h + v]    = r[2 + v];
            weights[size_t(MaxHangingLinks) * h + v] = RealType(1) / RealType(cnt);
        }
    }
};

namespace hanging_detail
{
template<class KeyType>
struct Packet
{
    KeyType keys[HangingRecordSize];
};
template<class KeyType>
struct Pack
{
    const KeyType* records;
    Packet<KeyType>* packets;
    __device__ void operator()(size_t i) const
    {
        for (int j = 0; j < HangingRecordSize; ++j)
            packets[i].keys[j] = records[HangingRecordSize * i + j];
    }
};
template<class KeyType>
struct Unpack
{
    const Packet<KeyType>* packets;
    KeyType* records;
    __device__ void operator()(size_t i) const
    {
        for (int j = 0; j < HangingRecordSize; ++j)
            records[HangingRecordSize * i + j] = packets[i].keys[j];
    }
};
template<class KeyType>
struct Questions
{
    const KeyType* keys;
    Packet<KeyType>* packets;
    __device__ void operator()(size_t i) const
    {
        packets[i]         = {};
        packets[i].keys[0] = keys[i];
    }
};
template<class KeyType>
struct QuestionKeys
{
    const Packet<KeyType>* packets;
    KeyType* keys;
    __device__ void operator()(size_t i) const { keys[i] = packets[i].keys[0]; }
};
} // namespace hanging_detail

// The hanging nodes among a rank's slots and the regular nodes behind each, by key, with weights summing to one.
// Built for the slots this rank's own elements use and the slots it owns: what a DofSpace classifies.
template<typename KeyType, typename RealType>
class HangingNodes
{
public:
    template<typename ElementTag>
    explicit HangingNodes(const ElementDomain<ElementTag, RealType, KeyType, cstone::execution::Gpu>& domain,
                          bool staged_messages = false)
        : staged_messages_(staged_messages)
    {
        MPI_Comm_dup(MPI_COMM_WORLD, &comm_);
        build(domain);
    }
    HangingNodes(const HangingNodes&)            = delete;
    HangingNodes& operator=(const HangingNodes&) = delete;
    ~HangingNodes() { MPI_Comm_free(&comm_); }

    // Hanging slots that own elements use or this rank owns, ascending. Hanging slot h has numLinks()[h] regular
    // nodes at linkKeys()/linkWeights() + MaxHangingLinks * h.
    size_t size() const { return slots_.size(); }
    const int* slots() const { return thrust::raw_pointer_cast(slots_.data()); }
    const int* numLinks() const { return thrust::raw_pointer_cast(numLinks_.data()); }
    const KeyType* linkKeys() const { return thrust::raw_pointer_cast(linkKeys_.data()); }
    const RealType* linkWeights() const { return thrust::raw_pointer_cast(linkWeights_.data()); }
    // Hanging nodes in the whole mesh, each counted once
    long long globalCount() const { return globalCount_; }
    // Slot count of the domain they were built for
    size_t slotCount() const { return n_; }

private:
    using Execution = amr::LeppDeviceExecution;
    using Packet    = hanging_detail::Packet<KeyType>;

    static Execution::Vector<Packet> pack(const KeyType* records, size_t n)
    {
        Execution::Vector<Packet> packets(n);
        auto* out = Execution::data(packets);
        Execution::for_each(n, hanging_detail::Pack<KeyType>{records, out});
        return packets;
    }

    static void unpack(const Execution::Vector<Packet>& packets, thrust::device_vector<KeyType>& records)
    {
        records.resize(HangingRecordSize * packets.size());
        auto* out      = thrust::raw_pointer_cast(records.data());
        const auto* in = Execution::data(packets);
        Execution::for_each(packets.size(), hanging_detail::Unpack<KeyType>{in, out});
    }

    template<typename DomainT>
    void query(const thrust::device_vector<KeyType>& keys, thrust::device_vector<KeyType>& answers,
               const DomainT& domain) const
    {
        const size_t n = keys.size();
        thrust::device_vector<int> owners(n);
        ownerOf(keys, 1, owners, domain);
        cstone::DeviceVector<int> destinations(n);
        Execution::copy(thrust::raw_pointer_cast(owners.data()), n, destinations.data());
        Execution::Vector<Packet> questions(n);
        auto* q         = Execution::data(questions);
        const auto* key = thrust::raw_pointer_cast(keys.data());
        Execution::for_each(n, hanging_detail::Questions<KeyType>{key, q});
        amr::AmrSparsePlan<Execution, Packet> plan(destinations, comm_, staged_messages_);
        const auto& incoming = plan.forward(questions);
        const auto* in       = Execution::data(incoming);
        cstone::DeviceVector<KeyType> asked(incoming.size());
        auto* ask = asked.data();
        Execution::for_each(incoming.size(), hanging_detail::QuestionKeys<KeyType>{in, ask});
        thrust::device_vector<KeyType> found(HangingRecordSize * incoming.size());
        lookup(ask, incoming.size(), thrust::raw_pointer_cast(found.data()));
        auto replies = pack(thrust::raw_pointer_cast(found.data()), incoming.size());
        unpack(plan.reverse(replies), answers);
    }

    template<typename ElementTag>
    void build(const ElementDomain<ElementTag, RealType, KeyType, cstone::execution::Gpu>& domain)
    {
        constexpr int N     = ElementTag::NodesPerElement;
        n_                 = domain.getNodeCount();
        rank_              = domain.rank();
        numRanks_          = domain.numRanks();
        const KeyType* key = domain.getLocalToGlobalSfcMap().data();
        int missing_history = !domain.has_refinement_history();
        MPI_Allreduce(MPI_IN_PLACE, &missing_history, 1, MPI_INT, MPI_MAX, comm_);
        if (missing_history)
        {
            MPI_Comm_free(&comm_);
            throw std::runtime_error("HangingNodes requires AMR refinement history from an initially conforming mesh");
        }
        if (numRanks_ > 1 && !domain.sfcOwnership())
        {
            if (rank_ == 0)
                std::fprintf(stderr, "HangingNodes: several ranks need SFC node ownership; unset MARS_OWNERSHIP\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        const auto& conn = domain.getElementToNodeConnectivity();
        HangingCorners<KeyType, N> c;
        [&]<std::size_t... I>(std::index_sequence<I...>) { ((c.c[I] = std::get<I>(conn).data()), ...); }
        (std::make_index_sequence<N>{});
        const auto& ancestry = domain.refinement_constraints();
        thrust::device_vector<KeyType> found(ancestry.size());
        if (!ancestry.empty())
            thrust::copy(thrust::device, thrust::device_pointer_cast(ancestry.data()),
                         thrust::device_pointer_cast(ancestry.data() + ancestry.size()), found.begin());
        long long conflicts = uniqueRecords(found);

        // 2. Each record goes to the SFC owner of its node, which keeps the union: the registry.
        thrust::device_vector<int> owner(numRecords(found));
        ownerOf(found, HangingRecordSize, owner, domain);
        sendToOwners(found, owner, registry_);
        conflicts += uniqueRecords(registry_);
        MPI_Allreduce(MPI_IN_PLACE, &conflicts, 1, MPI_LONG_LONG, MPI_SUM, comm_);
        if (conflicts > 0)
        {
            if (rank_ == 0)
                std::fprintf(stderr,
                             "HangingNodes: %lld nodes hang on different edges or faces of different elements;"
                             " the mesh is not 2:1 balanced\n",
                             conflicts);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        registryKeys_.resize(numRecords(registry_));
        thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0),
                          thrust::counting_iterator<size_t>(registryKeys_.size()), registryKeys_.begin(),
                          HangingKeyOf<KeyType>{thrust::raw_pointer_cast(registry_.data())});
        globalCount_ = (long long)registryKeys_.size();
        MPI_Allreduce(MPI_IN_PLACE, &globalCount_, 1, MPI_LONG_LONG, MPI_SUM, comm_);

        // 3. The slots own elements use and the owned ones (an owned node can have halo elements only), and their
        // records from the owners.
        thrust::device_vector<uint8_t> used(n_, 0);
        const uint8_t* owned = domain.getNodeOwnershipMap().data();
        thrust::transform(thrust::device, owned, owned + n_, used.begin(), HangingIsOwned{});
        const size_t first = domain.startIndex(), own = domain.localElementCount();
        for (int i = 0; i < N && own > 0; ++i)
            thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(first),
                             thrust::counting_iterator<size_t>(first + own),
                             HangingMarkCorner<KeyType>{c.c[i], thrust::raw_pointer_cast(used.data())});
        thrust::device_vector<int> usedSlots(n_);
        const size_t numUsed =
            thrust::copy_if(thrust::device, thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(int(n_)),
                            used.begin(), usedSlots.begin(), HangingIsSet{}) -
            usedSlots.begin();
        usedSlots.resize(numUsed);
        thrust::device_vector<KeyType> usedKeys(numUsed), usedStatus;
        thrust::gather(thrust::device, usedSlots.begin(), usedSlots.end(), thrust::device_pointer_cast(key),
                       usedKeys.begin());
        query(usedKeys, usedStatus, domain);

        // 4. The hanging used slots, with the nodes they hang on as first links.
        thrust::device_vector<uint8_t> isHanging(numUsed);
        thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(numUsed),
                          isHanging.begin(), HangingRecordHangs<KeyType>{thrust::raw_pointer_cast(usedStatus.data())});
        const size_t numHanging = mars::count(isHanging.begin(), isHanging.end(), uint8_t(1));
        thrust::device_vector<int> which(numHanging);
        thrust::copy_if(thrust::device, thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(int(numUsed)),
                        isHanging.begin(), which.begin(), HangingIsSet{});
        slots_.resize(numHanging);
        thrust::gather(thrust::device, which.begin(), which.end(), usedSlots.begin(), slots_.begin());
        numLinks_.resize(numHanging);
        linkKeys_.assign(size_t(MaxHangingLinks) * numHanging, ~KeyType(0));
        linkWeights_.assign(size_t(MaxHangingLinks) * numHanging, RealType(0));
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(numHanging),
                         HangingFirstLinks<KeyType, RealType>{
                             thrust::raw_pointer_cast(which.data()), thrust::raw_pointer_cast(usedStatus.data()),
                             thrust::raw_pointer_cast(numLinks_.data()), thrust::raw_pointer_cast(linkKeys_.data()),
                             thrust::raw_pointer_cast(linkWeights_.data())});

        // 5. A link to a node that hangs itself is replaced by that node's links, until all are regular.
        thrust::device_vector<KeyType> names, status;
        thrust::device_vector<int> flags(2);
        for (int round = 0;; ++round)
        {
            names.resize(linkKeys_.size());
            const size_t numNames = thrust::copy_if(thrust::device, linkKeys_.begin(), linkKeys_.end(), names.begin(),
                                                    HangingIsKey<KeyType>{}) -
                                    names.begin();
            names.resize(numNames);
            thrust::sort(thrust::device, names.begin(), names.end());
            names.resize(thrust::unique(thrust::device, names.begin(), names.end()) - names.begin());
            query(names, status, domain);

            thrust::fill(flags.begin(), flags.end(), 0);
            if (numHanging > 0)
                hangingFlattenKernel<KeyType, RealType><<<grid(numHanging), 256>>>(
                    numHanging, thrust::raw_pointer_cast(numLinks_.data()), thrust::raw_pointer_cast(linkKeys_.data()),
                    thrust::raw_pointer_cast(linkWeights_.data()), thrust::raw_pointer_cast(names.data()),
                    thrust::raw_pointer_cast(status.data()), names.size(), thrust::raw_pointer_cast(flags.data()),
                    thrust::raw_pointer_cast(flags.data()) + 1);
            cudaCheckError();
            int f[2];
            thrust::copy(flags.begin(), flags.end(), f);
            MPI_Allreduce(MPI_IN_PLACE, f, 2, MPI_INT, MPI_MAX, comm_);
            if (f[1])
            {
                if (rank_ == 0)
                    std::fprintf(stderr, "HangingNodes: more than %d regular nodes behind a hanging node\n",
                                 MaxHangingLinks);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            if (!f[0]) break;
            if (round == 8)
            {
                if (rank_ == 0) std::fprintf(stderr, "HangingNodes: hanging-node chains deeper than 8 levels\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    static size_t numRecords(const thrust::device_vector<KeyType>& records) { return records.size() / HangingRecordSize; }
    static int grid(size_t n) { return int((n + 255) / 256); }

    static void permute(const thrust::device_vector<KeyType>& items, const thrust::device_vector<int>& order, int width,
                        thrust::device_vector<KeyType>& out, bool scatter)
    {
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(order.size()),
                         HangingPermute<KeyType>{thrust::raw_pointer_cast(items.data()),
                                                 thrust::raw_pointer_cast(order.data()),
                                                 thrust::raw_pointer_cast(out.data()), width, scatter});
    }

    void lookup(const KeyType* keys, size_t n, KeyType* out) const
    {
        if (n == 0) return;
        hangingLookupKernel<KeyType><<<grid(n), 256>>>(keys, n, thrust::raw_pointer_cast(registryKeys_.data()),
                                                       thrust::raw_pointer_cast(registry_.data()), registryKeys_.size(),
                                                       out);
        cudaCheckError();
    }

    // Sorts records by node key and keeps one per key; returns the number of keys whose records differ.
    long long uniqueRecords(thrust::device_vector<KeyType>& records) const
    {
        const size_t n = numRecords(records);
        if (n == 0) return 0;
        thrust::device_vector<KeyType> keys(n);
        thrust::device_vector<int> order(n);
        thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                          keys.begin(), HangingKeyOf<KeyType>{thrust::raw_pointer_cast(records.data())});
        thrust::sequence(thrust::device, order.begin(), order.end());
        thrust::stable_sort_by_key(thrust::device, keys.begin(), keys.end(), order.begin());
        thrust::device_vector<KeyType> sorted(records.size());
        permute(records, order, HangingRecordSize, sorted, false);
        thrust::device_vector<uint8_t> keep(n), conflict(n);
        hangingFirstOfKeyKernel<KeyType><<<grid(n), 256>>>(thrust::raw_pointer_cast(sorted.data()), n,
                                                           thrust::raw_pointer_cast(keep.data()),
                                                           thrust::raw_pointer_cast(conflict.data()));
        cudaCheckError();
        const long long conflicts = mars::count(conflict.begin(), conflict.end(), uint8_t(1));
        thrust::device_vector<int> kept(n);
        const size_t numKept =
            thrust::copy_if(thrust::device, thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(int(n)),
                            keep.begin(), kept.begin(), HangingIsSet{}) -
            kept.begin();
        kept.resize(numKept);
        records.resize(size_t(HangingRecordSize) * numKept);
        permute(sorted, kept, HangingRecordSize, records, false);
        return conflicts;
    }

    template<typename DomainT>
    void ownerOf(const thrust::device_vector<KeyType>& keys, size_t stride, thrust::device_vector<int>& owner,
                 const DomainT& domain) const
    {
        const size_t n = owner.size();
        if (n == 0) return;
        if (numRanks_ > 1)
            hangingOwnerKernel<<<grid(n), 256>>>(thrust::raw_pointer_cast(keys.data()), stride, n,
                                                 domain.sfcNodeOwner(), thrust::raw_pointer_cast(owner.data()));
        else
            thrust::fill(owner.begin(), owner.end(), 0);
        cudaCheckError();
    }

    // Routes and payloads remain on device; only the sparse peer list reaches MPI on the host.
    void sendToOwners(const thrust::device_vector<KeyType>& records, thrust::device_vector<int>& owner,
                      thrust::device_vector<KeyType>& out) const
    {
        cstone::DeviceVector<int> destinations(owner.size());
        Execution::copy(thrust::raw_pointer_cast(owner.data()), owner.size(), destinations.data());
        auto packets = pack(thrust::raw_pointer_cast(records.data()), owner.size());
        amr::AmrSparsePlan<Execution, Packet> plan(destinations, comm_, staged_messages_);
        unpack(plan.forward(packets), out);
    }

    size_t n_ = 0;
    int rank_ = 0, numRanks_ = 1;
    MPI_Comm comm_         = MPI_COMM_NULL;
    long long globalCount_ = 0;
    thrust::device_vector<KeyType> registry_, registryKeys_; // records of the hanging nodes this rank owns
    thrust::device_vector<int> slots_, numLinks_;
    thrust::device_vector<KeyType> linkKeys_;
    thrust::device_vector<RealType> linkWeights_;
    bool staged_messages_ = false;
};

} // namespace fem
} // namespace mars
