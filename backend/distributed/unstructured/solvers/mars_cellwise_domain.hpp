#pragma once
// Cell-wise tables (mars_cellwise_topology.hpp) from a hexahedral ElementDomain: this
// rank's elements (the domain's local range, in its order) and, after them, the ghost
// layer, every element of another rank that shares a vertex with one of them.
//
// The ghost layer is completed through the node owners (SfcNodeOwner), so it does not
// depend on how far the cornerstone halo search reaches: every rank sends each of its
// elements to the owners of the element's vertices, so each owner holds the complete
// star of every node it owns, and each owner then sends every star to the ranks that
// have an element in it. A ghost needs only its id, owner rank and corner keys: the
// operator never runs on it, and the DSS reads its shared values from the exchange.

#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/mars_sfc_ownership.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"

#include <thrust/adjacent_difference.h>
#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/find.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/remove.h>
#include <thrust/sort.h>
#include <thrust/tuple.h>
#include <thrust/unique.h>

#include <memory>
#include <span>
#include <utility>
#include <vector>

namespace mars {
namespace cellwise {

namespace domain_detail {

using topo_detail::blocks;
using topo_detail::raw;

constexpr int kRecord = 10;   // element id, owner rank, 8 corner keys

template <class KeyType>
struct Corners {
    const KeyType* key[8];
};

template <class KeyType, class Domain, std::size_t... I>
Corners<KeyType> corner_keys(const Domain& d, std::index_sequence<I...>)
{
    return {{thrust::raw_pointer_cast(d.template indices<I>().data())...}};
}

template <class KeyType>
__global__ void local_records_kernel(Corners<KeyType> corners, const KeyType* __restrict__ codes, long long first,
                                     long long n, int rank, KeyType* __restrict__ rec)
{
    const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    KeyType* r = rec + i * kRecord;
    r[0] = codes[first + i];
    r[1] = (KeyType)rank;
    for (int c = 0; c < 8; ++c) r[2 + c] = corners.key[c][first + i];
}

// The ranks other than `self` that own a corner of record i, each once; -1 fills.
template <class KeyType, class RealType>
__global__ void owner_targets_kernel(const KeyType* __restrict__ rec, long long n,
                                     SfcNodeOwner<KeyType, RealType> owner, int self, int* __restrict__ dest)
{
    const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int seen[8];
    for (int c = 0; c < 8; ++c) {
        const int o = owner(rec[i * kRecord + 2 + c]);
        bool skip = o == self;
        for (int j = 0; j < c; ++j) skip = skip || seen[j] == o;
        seen[c] = o;
        dest[i * 8 + c] = skip ? -1 : o;
    }
}

// (vertex key, record) for every corner that `self` owns; key ~0 elsewhere.
template <class KeyType, class RealType>
__global__ void owned_corners_kernel(const KeyType* __restrict__ rec, long long n,
                                     SfcNodeOwner<KeyType, RealType> owner, int self, KeyType* __restrict__ key,
                                     long long* __restrict__ which)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n * 8) return;
    const KeyType k = rec[(t / 8) * kRecord + 2 + t % 8];
    key[t] = owner(k) == self ? k : ~KeyType(0);
    which[t] = t / 8;
}

// For each star of an owned vertex: every member to every rank (other than `self`)
// that owns another member. COUNT counts per star; FILL writes from the offsets.
template <bool FILL, class KeyType>
__global__ void star_targets_kernel(const KeyType* __restrict__ rec, const long long* __restrict__ which,
                                    const long long* __restrict__ off, long long stars, int self,
                                    long long* __restrict__ at, int* __restrict__ dest, long long* __restrict__ elem)
{
    const long long s = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= stars) return;
    long long o = FILL ? at[s] : 0;
    for (long long j = off[s]; j < off[s + 1]; ++j) {
        const int r = (int)rec[which[j] * kRecord + 1];
        bool first = r != self;   // the first member owned by r
        for (long long i = off[s]; i < j && first; ++i) first = (int)rec[which[i] * kRecord + 1] != r;
        if (!first) continue;
        for (long long m = off[s]; m < off[s + 1]; ++m)
            if ((int)rec[which[m] * kRecord + 1] != r) {
                if (FILL) {
                    dest[o] = r;
                    elem[o] = which[m];
                }
                ++o;
            }
    }
    if (!FILL) at[s] = o;
}

template <class KeyType>
__global__ void gather_records_kernel(const KeyType* __restrict__ rec, const long long* __restrict__ which,
                                      long long n, KeyType* __restrict__ out)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n * kRecord) out[t] = rec[which[t / kRecord] * kRecord + t % kRecord];
}

// 1 when a record shares a vertex with this rank's elements (sorted, unique keys).
template <class KeyType>
__global__ void touches_kernel(const KeyType* __restrict__ rec, long long n, const KeyType* __restrict__ local,
                               long long nlocal, int* __restrict__ flag)
{
    const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int hit = 0;
    for (int c = 0; c < 8 && !hit; ++c) {
        const KeyType k = rec[i * kRecord + 2 + c];
        long long lo = 0, hi = nlocal;
        while (lo < hi) {
            const long long mid = (lo + hi) / 2;
            if (local[mid] < k) lo = mid + 1;
            else hi = mid;
        }
        hit = lo < nlocal && local[lo] == k;
    }
    flag[i] = hit;
}

// Element arrays of the tables from the records: corner keys, global id, owner.
template <class KeyType>
__global__ void unpack_records_kernel(const KeyType* __restrict__ rec, long long n, unsigned long long* const* key,
                                      unsigned long long* __restrict__ gid, int* __restrict__ owner)
{
    const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    gid[i] = (unsigned long long)rec[i * kRecord];
    owner[i] = (int)rec[i * kRecord + 1];
    for (int c = 0; c < 8; ++c) key[c][i] = (unsigned long long)rec[i * kRecord + 2 + c];
}

// Sends the records rec[which[k]] to dest[k] (both sorted by dest) and returns what arrives.
template <class KeyType>
thrust::device_vector<KeyType> send_records(const thrust::device_vector<int>& dest,
                                            const thrust::device_vector<long long>& which, const KeyType* rec,
                                            int tag, MPI_Comm comm)
{
    const long long n = (long long)dest.size();
    thrust::device_vector<KeyType> buf(n * kRecord);
    if (n > 0) gather_records_kernel<<<blocks(n * kRecord), kThreads>>>(rec, raw(which), n, raw(buf));
    MARS_CELLWISE_CK(cudaGetLastError());
    thrust::device_vector<int> ranks(n), counts(n);
    const auto end = thrust::reduce_by_key(thrust::device, dest.begin(), dest.end(), thrust::constant_iterator<int>(1),
                                           ranks.begin(), counts.begin());
    const long long k = end.first - ranks.begin();
    std::vector<int> h_ranks(k), h_counts(k);
    thrust::copy(ranks.begin(), ranks.begin() + k, h_ranks.begin());
    thrust::copy(counts.begin(), counts.begin() + k, h_counts.begin());
    for (int& c : h_counts) c *= kRecord;
    MARS_CELLWISE_CK(cudaDeviceSynchronize());   // MPI reads buf next
    cstone::DeviceVector<KeyType> recv;
    sparseExchange<KeyType>(std::span<const int>(h_ranks), std::span<const int>(h_counts), raw(buf), recv, tag, comm);
    thrust::device_vector<KeyType> out(recv.size());
    if (recv.size() > 0)
        thrust::copy(thrust::device, recv.data(), recv.data() + recv.size(), out.begin());
    return out;
}

// Sorts (dest, which) by dest, then which, and drops repeated pairs.
inline void sort_unique_pairs(thrust::device_vector<int>& dest, thrust::device_vector<long long>& which)
{
    thrust::stable_sort_by_key(thrust::device, which.begin(), which.end(), dest.begin());
    thrust::stable_sort_by_key(thrust::device, dest.begin(), dest.end(), which.begin());
    const auto pairs = thrust::make_zip_iterator(thrust::make_tuple(dest.begin(), which.begin()));
    const long long n = thrust::unique(thrust::device, pairs, pairs + dest.size()) - pairs;
    dest.resize(n);
    which.resize(n);
}

}  // namespace domain_detail

// The tables of one rank, its exchange (null on one rank), and where its elements sit
// in the domain: local element i is domain element first + i.
struct DomainTables {
    UnstructuredTopology topo;
    std::unique_ptr<UnstructuredHalo> halo;
    long long first = 0;
};

template <class RealType, class KeyType>
DomainTables build_tables(const ElementDomain<HexTag, RealType, KeyType, cstone::execution::Gpu>& domain,
                          MPI_Comm comm)
{
    using namespace domain_detail;
    int rank = 0, size = 1;
    MARS_CELLWISE_MPI(MPI_Comm_rank(comm, &rank));
    MARS_CELLWISE_MPI(MPI_Comm_size(comm, &size));
    DomainTables out;
    out.first = (long long)domain.startIndex();
    const long long L = (long long)domain.localElementCount();

    thrust::device_vector<KeyType> local(L * kRecord);
    if (L > 0)
        local_records_kernel<<<blocks(L), kThreads>>>(corner_keys<KeyType>(domain, std::make_index_sequence<8>{}),
                                                      thrust::raw_pointer_cast(domain.getElementSfcCodes().data()),
                                                      out.first, L, rank, raw(local));
    MARS_CELLWISE_CK(cudaGetLastError());

    thrust::device_vector<KeyType> ghosts;
    if (size > 1) {
        if (!domain.sfcOwnership()) {   // the node owners are the meeting points of step 1
            if (rank == 0) fprintf(stderr, "cell-wise tables need SFC node ownership (MARS_OWNERSHIP=vote is set)\n");
            MPI_Abort(comm, 1);
        }
        const auto owner = domain.sfcNodeOwner();
        constexpr int kTagStars = 0x3340, kTagReturn = 0x3341;   // clear of the halo (0x3300) and MARS tags
        // 1. Every element to the owners of its vertices.
        thrust::device_vector<int> dest(L * 8);
        thrust::device_vector<long long> which(L * 8);
        if (L > 0) owner_targets_kernel<<<blocks(L), kThreads>>>(raw(local), L, owner, rank, raw(dest));
        MARS_CELLWISE_CK(cudaGetLastError());
        thrust::transform(thrust::device, thrust::counting_iterator<long long>(0),
                          thrust::counting_iterator<long long>(L * 8), which.begin(),
                          [] __host__ __device__(long long t) { return t / 8; });
        {
            const auto pairs = thrust::make_zip_iterator(thrust::make_tuple(dest.begin(), which.begin()));
            const long long n =
                thrust::remove_if(thrust::device, pairs, pairs + dest.size(),
                                  [] __host__ __device__(const thrust::tuple<int, long long>& p) {
                                      return thrust::get<0>(p) < 0;
                                  }) -
                pairs;
            dest.resize(n);
            which.resize(n);
        }
        sort_unique_pairs(dest, which);
        const thrust::device_vector<KeyType> pushed = send_records(dest, which, raw(local), kTagStars, comm);
        const long long A = (long long)pushed.size() / kRecord;

        // 2. The stars of the owned vertices, from this rank's elements and the pushed ones.
        thrust::device_vector<KeyType> pool(local);
        pool.insert(pool.end(), pushed.begin(), pushed.end());
        const long long P = L + A;
        thrust::device_vector<KeyType> vkey(P * 8);
        thrust::device_vector<long long> vrec(P * 8);
        if (P > 0)
            owned_corners_kernel<<<blocks(P * 8), kThreads>>>(raw(pool), P, owner, rank, raw(vkey), raw(vrec));
        MARS_CELLWISE_CK(cudaGetLastError());
        thrust::stable_sort_by_key(thrust::device, vkey.begin(), vkey.end(), vrec.begin());
        const long long owned =
            thrust::find(thrust::device, vkey.begin(), vkey.end(), ~KeyType(0)) - vkey.begin();
        thrust::device_vector<KeyType> star_key(owned);
        thrust::device_vector<long long> star_size(owned);
        const auto star_end = thrust::reduce_by_key(thrust::device, vkey.begin(), vkey.begin() + owned,
                                                    thrust::constant_iterator<long long>(1), star_key.begin(),
                                                    star_size.begin());
        const long long stars = star_end.first - star_key.begin();
        thrust::device_vector<long long> off(stars + 1, 0);
        thrust::inclusive_scan(thrust::device, star_size.begin(), star_size.begin() + stars, off.begin() + 1);

        // 3. Every star to the ranks with an element in it.
        thrust::device_vector<long long> at(stars);
        if (stars > 0)
            star_targets_kernel<false><<<blocks(stars), kThreads>>>(raw(pool), raw(vrec), raw(off), stars, rank,
                                                                    raw(at), nullptr, nullptr);
        MARS_CELLWISE_CK(cudaGetLastError());
        const long long n = stars > 0 ? thrust::reduce(thrust::device, at.begin(), at.end(), 0LL) : 0;
        thrust::exclusive_scan(thrust::device, at.begin(), at.end(), at.begin());
        thrust::device_vector<int> sdest(n);
        thrust::device_vector<long long> selem(n);
        if (stars > 0)
            star_targets_kernel<true><<<blocks(stars), kThreads>>>(raw(pool), raw(vrec), raw(off), stars, rank,
                                                                   raw(at), raw(sdest), raw(selem));
        MARS_CELLWISE_CK(cudaGetLastError());
        sort_unique_pairs(sdest, selem);
        const thrust::device_vector<KeyType> returned = send_records(sdest, selem, raw(pool), kTagReturn, comm);

        // 4. Ghosts: the pushed and returned elements, once each, that share a vertex
        // with this rank's elements.
        thrust::device_vector<KeyType> cand(pushed);
        cand.insert(cand.end(), returned.begin(), returned.end());
        const long long C = (long long)cand.size() / kRecord;
        thrust::device_vector<KeyType> id(C);
        thrust::device_vector<long long> order(C);
        const KeyType* cr = raw(cand);
        thrust::transform(thrust::device, thrust::counting_iterator<long long>(0),
                          thrust::counting_iterator<long long>(C), id.begin(),
                          [cr] __host__ __device__(long long i) { return cr[i * kRecord]; });
        thrust::sequence(thrust::device, order.begin(), order.end());
        thrust::stable_sort_by_key(thrust::device, id.begin(), id.end(), order.begin());
        const long long U = thrust::unique_by_key(thrust::device, id.begin(), id.end(), order.begin()).first - id.begin();
        order.resize(U);
        thrust::device_vector<KeyType> unique_rec(U * kRecord);
        if (U > 0) gather_records_kernel<<<blocks(U * kRecord), kThreads>>>(raw(cand), raw(order), U, raw(unique_rec));
        MARS_CELLWISE_CK(cudaGetLastError());

        thrust::device_vector<KeyType> lkeys(L * 8);
        const KeyType* lr = raw(local);
        thrust::transform(thrust::device, thrust::counting_iterator<long long>(0),
                          thrust::counting_iterator<long long>(L * 8), lkeys.begin(),
                          [lr] __host__ __device__(long long t) { return lr[(t / 8) * kRecord + 2 + t % 8]; });
        thrust::sort(thrust::device, lkeys.begin(), lkeys.end());
        const long long nl = thrust::unique(thrust::device, lkeys.begin(), lkeys.end()) - lkeys.begin();
        thrust::device_vector<int> touch(U, 0);
        if (U > 0) touches_kernel<<<blocks(U), kThreads>>>(raw(unique_rec), U, raw(lkeys), nl, raw(touch));
        MARS_CELLWISE_CK(cudaGetLastError());
        thrust::device_vector<long long> keep(U);
        const long long G =
            thrust::copy_if(thrust::device, thrust::counting_iterator<long long>(0),
                            thrust::counting_iterator<long long>(U), touch.begin(), keep.begin(),
                            [] __host__ __device__(int t) { return t != 0; }) - keep.begin();
        keep.resize(G);
        ghosts.resize(G * kRecord);
        if (G > 0) gather_records_kernel<<<blocks(G * kRecord), kThreads>>>(raw(unique_rec), raw(keep), G, raw(ghosts));
        MARS_CELLWISE_CK(cudaGetLastError());
    }

    // The tables over this rank's elements and the ghosts.
    thrust::device_vector<KeyType> all(local);
    all.insert(all.end(), ghosts.begin(), ghosts.end());
    const long long E = (long long)all.size() / kRecord;
    thrust::device_vector<unsigned long long> key[8], gid(E);
    thrust::device_vector<int> lid[8], owner(E);
    std::vector<unsigned long long*> kp(8);
    for (int c = 0; c < 8; ++c) {
        key[c].resize(E);
        kp[c] = raw(key[c]);
    }
    thrust::device_vector<unsigned long long*> d_kp(kp.begin(), kp.end());
    if (E > 0) unpack_records_kernel<<<blocks(E), kThreads>>>(raw(all), E, raw(d_kp), raw(gid), raw(owner));
    MARS_CELLWISE_CK(cudaGetLastError());

    // The exchange orders shared nodes by element id, so ids must not repeat here.
    {
        thrust::device_vector<unsigned long long> sorted(gid);
        thrust::sort(thrust::device, sorted.begin(), sorted.end());
        int repeated = E > 1 && thrust::unique(thrust::device, sorted.begin(), sorted.end()) - sorted.begin() != E;
        MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &repeated, 1, MPI_INT, MPI_MAX, comm));
        if (repeated) {
            if (rank == 0) fprintf(stderr, "cell-wise tables: two elements share an SFC code on one rank\n");
            MPI_Abort(comm, 1);
        }
    }

    // Local corner ids: ranks of the corner keys among this rank's distinct keys.
    thrust::device_vector<unsigned long long> distinct(E * 8);
    for (int c = 0; c < 8; ++c) thrust::copy(key[c].begin(), key[c].end(), distinct.begin() + c * E);
    thrust::sort(thrust::device, distinct.begin(), distinct.end());
    const long long nd = thrust::unique(thrust::device, distinct.begin(), distinct.end()) - distinct.begin();
    std::vector<int*> lp(8);
    for (int c = 0; c < 8; ++c) {
        lid[c].resize(E);
        thrust::lower_bound(thrust::device, distinct.begin(), distinct.begin() + nd, key[c].begin(), key[c].end(),
                            lid[c].begin());
        lp[c] = raw(lid[c]);
    }
    const std::vector<const unsigned long long*> ckp(kp.begin(), kp.end());
    const std::vector<const int*> clp(lp.begin(), lp.end());
    const thrust::device_vector<const unsigned long long*> d_ckp(ckp.begin(), ckp.end());
    const thrust::device_vector<const int*> d_clp(clp.begin(), clp.end());
    out.topo = build_topology(thrust::raw_pointer_cast(d_ckp.data()), thrust::raw_pointer_cast(d_clp.data()), E, L);
    if (size > 1) out.halo.reset(new UnstructuredHalo(out.topo, raw(gid), raw(owner), comm));
    return out;
}

}  // namespace cellwise
}  // namespace mars
