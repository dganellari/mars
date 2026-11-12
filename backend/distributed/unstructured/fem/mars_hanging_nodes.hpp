#pragma once

// Hanging nodes of a non-conforming mesh, and the regular nodes each one takes its value from.
//
// Refining an element but not its neighbour leaves nodes in the middle of the neighbour's edges and, on
// hexahedra, of its faces. A continuous P1/Q1 field is linear along that edge and bilinear on that face, so its
// value at such a node is the mean of the edge ends (face corners): the node is not an unknown. When an edge end is
// itself hanging on a coarser element, its own constraint is substituted, down to regular nodes.
//
// Detection: each element E compares the corners of the elements around its corners with the midpoints of its
// edges and the centres of its faces. Positions match within a tenth of the edge, not by key: every AMR rebuild
// encodes the coordinates against a new box. Each rank sends what it found to the SFC owner of the node, then asks
// the owners about every node its own elements use and every node a constraint refers to, so all ranks agree
// whatever their halos hold.

#include "../domain.hpp"
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

template<typename ElementTag>
struct HangingGeometry;

// Corner numbering of the hexahedral refiner (mars_amr_octree_refine.hpp)
template<>
struct HangingGeometry<HexTag>
{
    static constexpr int numNodes = 8, numEdges = 12, numFaces = 6;
    __host__ __device__ static void edge(int k, int& a, int& b)
    {
        constexpr int e[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                  {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
        a = e[k][0];
        b = e[k][1];
    }
    __host__ __device__ static void face(int k, int (&c)[4])
    {
        constexpr int f[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 2, 6, 7}, {0, 3, 7, 4}, {1, 2, 6, 5}};
        for (int i = 0; i < 4; ++i)
            c[i] = f[k][i];
    }
};

// Red refinement of tetrahedra only puts nodes on edges.
template<>
struct HangingGeometry<TetTag>
{
    static constexpr int numNodes = 4, numEdges = 6, numFaces = 0;
    __host__ __device__ static void edge(int k, int& a, int& b)
    {
        constexpr int e[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
        a = e[k][0];
        b = e[k][1];
    }
    __host__ __device__ static void face(int, int (&)[4]) {}
};

template<typename KeyType, int N>
struct HangingCorners
{
    const KeyType* c[N];
};

// The hanging nodes on the edges and faces of each element E, from the corners of the elements around E's
// corners. Count pass (records == nullptr): count[e]. Write pass: from offset[e], records of
// {node slot, 2 or 4, master slots, -1 for unused}.
template<typename ElementTag, typename KeyType, typename RealType>
__global__ void hangingFindKernel(HangingCorners<KeyType, HangingGeometry<ElementTag>::numNodes> conn,
                                  size_t numElements, const RealType* x, const RealType* y, const RealType* z,
                                  const KeyType* n2eOffsets, const KeyType* n2eList, const unsigned* offset,
                                  unsigned* count, int* records)
{
    using G                = HangingGeometry<ElementTag>;
    constexpr int N        = G::numNodes;
    constexpr int maxFound = G::numEdges + G::numFaces; // one node per edge and face in a 2:1 mesh
    const size_t e         = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    int n[N];
    RealType X[N][3], lo[3], hi[3], shortest2 = 0;
    for (int i = 0; i < N; ++i)
    {
        n[i]    = int(conn.c[i][e]);
        X[i][0] = x[n[i]];
        X[i][1] = y[n[i]];
        X[i][2] = z[n[i]];
    }
    for (int t = 0; t < 3; ++t)
    {
        lo[t] = hi[t] = X[0][t];
        for (int i = 1; i < N; ++i)
        {
            lo[t] = X[i][t] < lo[t] ? X[i][t] : lo[t];
            hi[t] = X[i][t] > hi[t] ? X[i][t] : hi[t];
        }
    }
    for (int q = 0; q < G::numEdges; ++q)
    {
        int a, b;
        G::edge(q, a, b);
        RealType l2 = 0;
        for (int t = 0; t < 3; ++t)
            l2 += (X[b][t] - X[a][t]) * (X[b][t] - X[a][t]);
        shortest2 = (q == 0 || l2 < shortest2) ? l2 : shortest2;
    }
    // A node on an edge or face lies in E's box; the slack covers the tolerance of the tests below.
    const RealType slack = RealType(0.1) * sqrt(shortest2);

    int found[maxFound];
    int numFound = 0;
    for (int i = 0; i < N && numFound < maxFound; ++i)
    {
        for (KeyType j = n2eOffsets[n[i]]; j < n2eOffsets[n[i] + 1] && numFound < maxFound; ++j)
        {
            const size_t f = size_t(n2eList[j]);
            if (f == e) continue;
            for (int k = 0; k < N && numFound < maxFound; ++k)
            {
                const int m = int(conn.c[k][f]);
                bool skip   = false;
                for (int v = 0; v < N; ++v)
                    skip |= m == n[v];
                for (int v = 0; v < numFound; ++v)
                    skip |= m == found[v];
                if (skip) continue;
                const RealType xm[3] = {x[m], y[m], z[m]};
                bool inside          = true;
                for (int t = 0; t < 3; ++t)
                    inside &= xm[t] >= lo[t] - slack && xm[t] <= hi[t] + slack;
                if (!inside) continue;

                int masters[4]  = {-1, -1, -1, -1};
                int numMasters  = 0;
                for (int q = 0; q < G::numEdges && numMasters == 0; ++q)
                {
                    int a, b;
                    G::edge(q, a, b);
                    RealType d2 = 0, l2 = 0;
                    for (int t = 0; t < 3; ++t)
                    {
                        const RealType mid = RealType(0.5) * (X[a][t] + X[b][t]);
                        d2 += (xm[t] - mid) * (xm[t] - mid);
                        l2 += (X[b][t] - X[a][t]) * (X[b][t] - X[a][t]);
                    }
                    if (d2 < RealType(0.01) * l2)
                    {
                        masters[0] = n[a];
                        masters[1] = n[b];
                        numMasters = 2;
                    }
                }
                for (int q = 0; q < G::numFaces && numMasters == 0; ++q)
                {
                    int c[4];
                    G::face(q, c);
                    RealType d2 = 0, diag2 = 0;
                    for (int t = 0; t < 3; ++t)
                    {
                        const RealType centre = RealType(0.25) * (X[c[0]][t] + X[c[1]][t] + X[c[2]][t] + X[c[3]][t]);
                        d2 += (xm[t] - centre) * (xm[t] - centre);
                        diag2 += (X[c[2]][t] - X[c[0]][t]) * (X[c[2]][t] - X[c[0]][t]);
                    }
                    // a tenth of the half diagonal
                    if (d2 < RealType(0.0025) * diag2)
                    {
                        for (int v = 0; v < 4; ++v)
                            masters[v] = n[c[v]];
                        numMasters = 4;
                    }
                }
                if (numMasters == 0) continue;
                if (records)
                {
                    int* r = records + size_t(HangingRecordSize) * (offset[e] + numFound);
                    r[0]   = m;
                    r[1]   = numMasters;
                    for (int v = 0; v < 4; ++v)
                        r[2 + v] = masters[v];
                }
                found[numFound++] = m;
            }
        }
    }
    if (!records) count[e] = unsigned(numFound);
}

// Slot records to key records, the master keys ascending, so equal constraints compare equal
template<typename KeyType>
__global__ void hangingKeyRecordsKernel(const int* slotRecords, size_t numRecords, const KeyType* keys,
                                        KeyType* keyRecords, KeyType* nodeKeys)
{
    const size_t r = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (r >= numRecords) return;
    const int* in = slotRecords + size_t(HangingRecordSize) * r;
    KeyType* out  = keyRecords + size_t(HangingRecordSize) * r;
    const int cnt = in[1];
    KeyType m[4];
    for (int v = 0; v < 4; ++v)
        m[v] = v < cnt ? keys[in[2 + v]] : ~KeyType(0);
    for (int a = 1; a < 4; ++a) // insertion sort of four
        for (int b = a; b > 0 && m[b] < m[b - 1]; --b)
        {
            KeyType t = m[b];
            m[b]      = m[b - 1];
            m[b - 1]  = t;
        }
    out[0] = keys[in[0]];
    out[1] = KeyType(cnt);
    for (int v = 0; v < 4; ++v)
        out[2 + v] = m[v];
    nodeKeys[r] = out[0];
}

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

// The hanging nodes among a rank's slots and the regular nodes behind each, by key, with weights summing to one.
// Built for the slots this rank's own elements use and the slots it owns: what a DofSpace classifies.
template<typename KeyType, typename RealType>
class HangingNodes
{
public:
    template<typename ElementTag>
    explicit HangingNodes(const ElementDomain<ElementTag, RealType, KeyType, cstone::execution::Gpu>& domain)
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
    // Collective. The record of each key from its owner: the nodes it hangs on, or count 0 for a regular node.
    template<typename DomainT>
    void query(const thrust::device_vector<KeyType>& keys, thrust::device_vector<KeyType>& answers,
               const DomainT& domain) const
    {
        constexpr int W = HangingRecordSize;
        const size_t n  = keys.size();
        thrust::device_vector<int> owner(n), order(n);
        ownerOf(keys, 1, owner, domain);
        thrust::sequence(thrust::device, order.begin(), order.end());
        thrust::stable_sort_by_key(thrust::device, owner.begin(), owner.end(), order.begin());
        thrust::device_vector<KeyType> sortedKeys(n);
        permute(keys, order, 1, sortedKeys, false);

        // runs of equal owner in the sorted keys
        std::vector<int> h_owner(n);
        thrust::copy(owner.begin(), owner.end(), h_owner.begin());
        std::vector<int> dests, counts;
        std::vector<size_t> begins;
        size_t selfBegin = 0, selfCount = 0;
        for (size_t i = 0; i < n;)
        {
            size_t j = i;
            while (j < n && h_owner[j] == h_owner[i])
                ++j;
            if (h_owner[i] == rank_)
            {
                selfBegin = i;
                selfCount = j - i;
            }
            else
            {
                dests.push_back(h_owner[i]);
                counts.push_back(int(j - i));
                begins.push_back(i);
            }
            i = j;
        }

        // questions to the owners
        thrust::device_vector<KeyType> sendKeys(n - selfCount);
        thrust::copy(sortedKeys.begin(), sortedKeys.begin() + selfBegin, sendKeys.begin());
        thrust::copy(sortedKeys.begin() + selfBegin + selfCount, sortedKeys.end(), sendKeys.begin() + selfBegin);
        cudaCheckCall(cudaDeviceSynchronize()); // GPU-aware MPI reads sendKeys
        cstone::DeviceVector<KeyType> asked;
        auto questions = sparseExchange<KeyType>(dests, counts, thrust::raw_pointer_cast(sendKeys.data()), asked,
                                                 tagQueries, comm_);

        // answers, one message per asker in source order
        const size_t numAsked = asked.size();
        thrust::device_vector<KeyType> lookedUp(size_t(W) * numAsked), replyBuf(size_t(W) * numAsked);
        lookup(asked.data(), numAsked, thrust::raw_pointer_cast(lookedUp.data()));
        std::vector<int> replyDests, replyCounts;
        size_t at = 0;
        for (const auto& q : questions)
        {
            thrust::copy(lookedUp.begin() + size_t(W) * q.offset, lookedUp.begin() + size_t(W) * (q.offset + q.count),
                         replyBuf.begin() + at);
            at += size_t(W) * q.count;
            replyDests.push_back(q.source);
            replyCounts.push_back(W * q.count);
        }
        cudaCheckCall(cudaDeviceSynchronize());
        cstone::DeviceVector<KeyType> got;
        auto replies = sparseExchange<KeyType>(replyDests, replyCounts, thrust::raw_pointer_cast(replyBuf.data()), got,
                                               tagReplies, comm_);
        if (replies.size() != dests.size())
        {
            std::fprintf(stderr, "HangingNodes: rank %d asked %zu ranks and got %zu replies\n", rank_, dests.size(),
                         replies.size());
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // in the sorted-by-owner order: own keys answered here, the others from the replies (sources ascend, as dests)
        thrust::device_vector<KeyType> sortedAnswers(size_t(W) * n);
        if (selfCount > 0)
            lookup(thrust::raw_pointer_cast(sortedKeys.data()) + selfBegin, selfCount,
                   thrust::raw_pointer_cast(sortedAnswers.data()) + size_t(W) * selfBegin);
        for (size_t d = 0; d < dests.size(); ++d)
        {
            auto src = thrust::device_pointer_cast(got.data()) + replies[d].offset;
            thrust::copy(src, src + size_t(W) * counts[d], sortedAnswers.begin() + size_t(W) * begins[d]);
        }
        answers.resize(size_t(W) * n);
        permute(sortedAnswers, order, W, answers, true);
    }

    template<typename ElementTag>
    void build(const ElementDomain<ElementTag, RealType, KeyType, cstone::execution::Gpu>& domain)
    {
        constexpr int N    = HangingGeometry<ElementTag>::numNodes;
        n_                 = domain.getNodeCount();
        rank_              = domain.rank();
        numRanks_          = domain.numRanks();
        const KeyType* key = domain.getLocalToGlobalSfcMap().data();
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
        const size_t numElements = domain.getElementCount();
        const auto& n2eOffsets   = domain.getNodeToElementOffsets();
        const auto& n2eList      = domain.getNodeToElementList();
        const RealType *x = domain.getNodeX().data(), *y = domain.getNodeY().data(), *z = domain.getNodeZ().data();

        // 1. Every element this rank holds, own and halo, finds the hanging nodes on its edges and faces.
        thrust::device_vector<unsigned> count(numElements + 1, 0u), offset(numElements + 1, 0u);
        if (numElements > 0)
            hangingFindKernel<ElementTag, KeyType, RealType>
                <<<grid(numElements), 256>>>(c, numElements, x, y, z, n2eOffsets.data(), n2eList.data(), nullptr,
                                             thrust::raw_pointer_cast(count.data()), nullptr);
        cudaCheckError();
        thrust::exclusive_scan(thrust::device, count.begin(), count.end(), offset.begin(), 0u);
        unsigned numFound = offset.back();
        thrust::device_vector<int> slotRecords(size_t(HangingRecordSize) * numFound);
        if (numFound > 0)
            hangingFindKernel<ElementTag, KeyType, RealType><<<grid(numElements), 256>>>(
                c, numElements, x, y, z, n2eOffsets.data(), n2eList.data(), thrust::raw_pointer_cast(offset.data()),
                nullptr, thrust::raw_pointer_cast(slotRecords.data()));
        cudaCheckError();
        thrust::device_vector<KeyType> found(size_t(HangingRecordSize) * numFound), foundKeys(numFound);
        if (numFound > 0)
            hangingKeyRecordsKernel<KeyType><<<grid(numFound), 256>>>(
                thrust::raw_pointer_cast(slotRecords.data()), numFound, key, thrust::raw_pointer_cast(found.data()),
                thrust::raw_pointer_cast(foundKeys.data()));
        cudaCheckError();
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

    // Records to the owners of their nodes; out: this rank's own records, then the received ones.
    void sendToOwners(const thrust::device_vector<KeyType>& records, thrust::device_vector<int>& owner,
                      thrust::device_vector<KeyType>& out) const
    {
        constexpr int W = HangingRecordSize;
        const size_t n  = owner.size();
        thrust::device_vector<int> order(n);
        thrust::sequence(thrust::device, order.begin(), order.end());
        thrust::stable_sort_by_key(thrust::device, owner.begin(), owner.end(), order.begin());
        thrust::device_vector<KeyType> sorted(records.size());
        permute(records, order, W, sorted, false);
        std::vector<int> h_owner(n), dests, counts;
        thrust::copy(owner.begin(), owner.end(), h_owner.begin());
        size_t selfBegin = 0, selfCount = 0;
        for (size_t i = 0; i < n;)
        {
            size_t j = i;
            while (j < n && h_owner[j] == h_owner[i])
                ++j;
            if (h_owner[i] == rank_)
            {
                selfBegin = i;
                selfCount = j - i;
            }
            else
            {
                dests.push_back(h_owner[i]);
                counts.push_back(int(W * (j - i)));
            }
            i = j;
        }
        thrust::device_vector<KeyType> sendBuf(size_t(W) * (n - selfCount));
        auto s = sorted.begin();
        thrust::copy(s, s + W * selfBegin, sendBuf.begin());
        thrust::copy(s + W * (selfBegin + selfCount), sorted.end(), sendBuf.begin() + W * selfBegin);
        cudaCheckCall(cudaDeviceSynchronize()); // GPU-aware MPI reads sendBuf
        cstone::DeviceVector<KeyType> recv;
        sparseExchange<KeyType>(dests, counts, thrust::raw_pointer_cast(sendBuf.data()), recv, tagRecords, comm_);
        out.resize(size_t(W) * selfCount + recv.size());
        thrust::copy(s + W * selfBegin, s + W * (selfBegin + selfCount), out.begin());
        if (recv.size() > 0)
            thrust::copy(thrust::device_pointer_cast(recv.data()), thrust::device_pointer_cast(recv.data()) + recv.size(),
                         out.begin() + size_t(W) * selfCount);
    }

    size_t n_ = 0;
    int rank_ = 0, numRanks_ = 1;
    MPI_Comm comm_         = MPI_COMM_NULL;
    long long globalCount_ = 0;
    thrust::device_vector<KeyType> registry_, registryKeys_; // records of the hanging nodes this rank owns
    thrust::device_vector<int> slots_, numLinks_;
    thrust::device_vector<KeyType> linkKeys_;
    thrust::device_vector<RealType> linkWeights_;
    static constexpr int tagRecords = 0x4e70, tagQueries = 0x4e71, tagReplies = 0x4e72;
};

} // namespace fem
} // namespace mars
