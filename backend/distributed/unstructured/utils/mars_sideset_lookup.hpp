#pragma once

#include "mars_read_exodus_mesh.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#include <thrust/transform.h>
#include <thrust/execution_policy.h>
#define MARS_SIDESET_HD __host__ __device__
#else
#define MARS_SIDESET_HD
#endif

namespace mars
{
template<class KeyType>
struct SideSetHit
{
    uint64_t id;
    int block;
    int local;
    KeyType key;
};

// File blocks distinguish coincident interface nodes. Search the current map,
// including ghosts; the initial file partition no longer identifies local nodes.
template<class KeyType, class Encode>
struct SideSetLookup
{
    const KeyType* keys;
    const int* blocks;
    size_t count;
    Encode encode;

    MARS_SIDESET_HD SideSetHit<KeyType> operator()(const ExodusSideSetNode& node) const
    {
        const KeyType key = encode(node);
        size_t lo = 0, hi = count;
        while (lo < hi)
        {
            const size_t mid = lo + (hi - lo) / 2;
            if (keys[mid] < key)
                lo = mid + 1;
            else
                hi = mid;
        }
        if (blocks)
            while (lo < count && keys[lo] == key && blocks[lo] != node.block)
                ++lo;
        const int local = lo < count && keys[lo] == key ? int(lo) : -1;
        // Material tags do not split nodes unless the domain keeps block identities.
        return {node.id, blocks ? node.block : 0, local, key};
    }
};

template<class KeyType>
struct SideSetHitValid
{
    MARS_SIDESET_HD bool operator()(const SideSetHit<KeyType>& hit) const { return hit.local >= 0; }
};

template<class Encode, class Prefix>
struct SideSetNodePrefix
{
    Encode encode;
    Prefix prefix;
    MARS_SIDESET_HD uint64_t operator()(const ExodusSideSetNode& node) const { return prefix(encode(node)); }
};

template<class Encode, class Prefix>
class SideSetBatchPrefixes
{
    SideSetNodePrefix<Encode, Prefix> encode_;
#if defined(__CUDACC__) || defined(__HIPCC__)
    thrust::device_vector<ExodusSideSetNode> d_nodes_;
    thrust::device_vector<uint64_t> d_prefixes_;
#endif
public:
    SideSetBatchPrefixes(Encode encode, Prefix prefix)
        : encode_{encode, prefix}
    {
    }

    std::vector<uint64_t> operator()(const ExodusSideSetBatch& batch)
    {
        std::vector<uint64_t> result(batch.nodes.size());
#if defined(__CUDACC__) || defined(__HIPCC__)
        // Match the domain's GPU key encoding even at floating-point cell boundaries.
        d_nodes_.assign(batch.nodes.begin(), batch.nodes.end());
        d_prefixes_.resize(batch.nodes.size());
        thrust::transform(thrust::device, d_nodes_.begin(), d_nodes_.end(), d_prefixes_.begin(), encode_);
        // Routing needs only these prefixes on the host, not the domain's node-key array.
        thrust::copy(d_prefixes_.begin(), d_prefixes_.end(), result.begin());
#else
        std::transform(batch.nodes.begin(), batch.nodes.end(), result.begin(), encode_);
#endif
        return result;
    }
};

template<class KeyType, class Encode>
class SideSetBatchLookup
{
    using Hit = SideSetHit<KeyType>;
    SideSetLookup<KeyType, Encode> lookup_;
#if defined(__CUDACC__) || defined(__HIPCC__)
    thrust::device_vector<ExodusSideSetNode> d_nodes_;
    thrust::device_vector<Hit> d_hits_, d_selected_;
#endif

public:
    explicit SideSetBatchLookup(SideSetLookup<KeyType, Encode> lookup)
        : lookup_(lookup)
    {
    }

    // Face consumers keep misses to retain the face-index alignment. Node-only
    // consumers download hits alone. Scratch is reused across side sets/batches.
    std::vector<Hit> resolve(const ExodusSideSetBatch& batch, bool keep_misses = false)
    {
#if defined(__CUDACC__) || defined(__HIPCC__)
        d_nodes_.resize(batch.nodes.size());
        d_hits_.resize(batch.nodes.size());
        thrust::copy(batch.nodes.begin(), batch.nodes.end(), d_nodes_.begin());
        thrust::transform(thrust::device, d_nodes_.begin(), d_nodes_.end(), d_hits_.begin(), lookup_);
        if (keep_misses)
        {
            std::vector<Hit> result(d_hits_.size());
            thrust::copy(d_hits_.begin(), d_hits_.end(), result.begin());
            return result;
        }
        d_selected_.resize(batch.nodes.size());
        auto end = thrust::copy_if(thrust::device, d_hits_.begin(), d_hits_.end(), d_selected_.begin(),
                                   SideSetHitValid<KeyType>{});
        std::vector<Hit> result(end - d_selected_.begin());
        thrust::copy(d_selected_.begin(), end, result.begin());
        return result;
#else
        std::vector<Hit> result;
        for (const auto& node : batch.nodes)
        {
            const auto hit = lookup_(node);
            if (keep_misses || hit.local >= 0) result.push_back(hit);
        }
        return result;
#endif
    }
};
} // namespace mars

#undef MARS_SIDESET_HD
