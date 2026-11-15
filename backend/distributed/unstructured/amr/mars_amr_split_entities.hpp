#pragma once

#include "mars_amr_sparse_plan.hpp"
#include <array>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_ENTITY_HD __host__ __device__
#else
#define MARS_ENTITY_HD
#endif

namespace mars::amr
{
// Sorted parent keys identify an edge or a quadrilateral face independently of rank and orientation.
template<class KeyType>
struct RefinementEntity
{
    KeyType corners[4] = {~KeyType(0), ~KeyType(0), ~KeyType(0), ~KeyType(0)};
    MARS_ENTITY_HD int size() const { return corners[2] == ~KeyType(0) ? 2 : 4; }
    MARS_ENTITY_HD void sort()
    {
        for (int a = 1; a < 4; ++a)
            for (int b = a; b && corners[b] < corners[b - 1]; --b)
            {
                auto tmp       = corners[b];
                corners[b]     = corners[b - 1];
                corners[b - 1] = tmp;
            }
    }
    MARS_ENTITY_HD bool operator==(const RefinementEntity& rhs) const
    {
        for (int i = 0; i < 4; ++i)
            if (corners[i] != rhs.corners[i]) return false;
        return true;
    }
    MARS_ENTITY_HD bool operator<(const RefinementEntity& rhs) const
    {
        for (int i = 0; i < 4; ++i)
        {
            if (corners[i] < rhs.corners[i]) return true;
            if (corners[i] > rhs.corners[i]) return false;
        }
        return false;
    }
};

template<class KeyType>
struct SplitEntity
{
    RefinementEntity<KeyType> parent;
    KeyType node = ~KeyType(0);
};

namespace split_detail
{
template<class KeyType>
struct IdentityPacket
{
    SplitEntity<KeyType> split;
    int kind;
};
template<class KeyType>
struct MatchPacket
{
    SplitEntity<KeyType> split;
    uint64_t slot;
};
} // namespace split_detail

// Same ordering as the red refiners. No physical length enters incidence or balance.
template<int N>
struct RefinementEntities
{
    static_assert(N == 4 || N == 8);
    static constexpr int edges = N == 4 ? 6 : 12;
    static constexpr int count = N == 4 ? 6 : 18;
    MARS_ENTITY_HD static int vertices(int entity, int (&v)[4])
    {
        constexpr int tet[6][2]  = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
        constexpr int hex[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                    {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
        constexpr int face[6][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 2, 6, 7}, {0, 3, 7, 4}, {1, 2, 6, 5}};
        if (entity < edges)
        {
            v[0] = N == 4 ? tet[entity][0] : hex[entity][0];
            v[1] = N == 4 ? tet[entity][1] : hex[entity][1];
            return 2;
        }
        for (int i = 0; i < 4; ++i)
            v[i] = face[entity - edges][i];
        return 4;
    }
};

template<class Execution, class KeyType>
typename Execution::template Vector<SplitEntity<KeyType>>
compact_split_entities(const typename Execution::template Vector<SplitEntity<KeyType>>& input)
{
    using U64 = typename Execution::template Vector<uint64_t>;
    size_t n  = input.size();
    U64 count(n + 1), offset(n + 1);
    auto in = Execution::data(input);
    auto c  = Execution::data(count);
    auto o  = Execution::data(offset);
    Execution::for_each(n + 1, [=] MARS_ENTITY_HD(size_t i) { c[i] = i < n && in[i].node != ~KeyType(0); });
    Execution::scan(c, n + 1, o);
    typename Execution::template Vector<SplitEntity<KeyType>> result(Execution::read(o + n));
    auto out = Execution::data(result);
    Execution::for_each(n,
                        [=] MARS_ENTITY_HD(size_t i)
                        {
                            if (c[i]) out[o[i]] = in[i];
                        });
    Execution::synchronize();
    return result;
}

// A new midpoint may reuse an old key only when that old node already has this exact parent entity.
// This catches quantization collisions with unrelated vertices before a rebuild can merge their keys.
template<class Execution, class KeyType, class Owner>
void validate_split_identity(const typename Execution::template Vector<SplitEntity<KeyType>>& sources,
                             size_t previous_count,
                             const KeyType* old_keys,
                             const uint8_t* owned,
                             size_t nodes,
                             Owner owner,
                             MPI_Comm comm,
                             bool staged,
                             const KeyType* body_keys = nullptr,
                             size_t bodies            = 0)
{
    using Packet   = split_detail::IdentityPacket<KeyType>;
    const size_t n = sources.size() + nodes + bodies;
    typename Execution::template Vector<Packet> packets(n);
    typename Execution::template Vector<int> destinations(n);
    auto p = Execution::data(packets);
    auto d = Execution::data(destinations);
    auto s = Execution::data(sources);
    int rank;
    MPI_Comm_rank(comm, &rank);
    const size_t ns = sources.size();
    Execution::for_each(n,
                        [=] MARS_ENTITY_HD(size_t i)
                        {
                            if (i < ns)
                            {
                                p[i] = {s[i], i < previous_count ? 1 : 2};
                                d[i] = owner(s[i].node);
                            }
                            else if (i < ns + nodes)
                            {
                                const size_t j = i - ns;
                                p[i]           = {SplitEntity<KeyType>{{}, old_keys[j]}, owned[j] ? 0 : -1};
                                d[i]           = rank;
                            }
                            else
                            {
                                p[i] = {SplitEntity<KeyType>{{}, body_keys[i - ns - nodes]}, 3};
                                d[i] = owner(p[i].split.node);
                            }
                        });
    AmrSparsePlan<Execution, Packet> plan(destinations, comm, staged);
    auto received   = plan.forward(packets);
    auto r          = Execution::data(received);
    const size_t nr = received.size();
    Execution::sort(r, nr, [] MARS_ENTITY_HD(const Packet& a, const Packet& b)
                    { return a.split.node < b.split.node || (a.split.node == b.split.node && a.kind < b.kind); });
    typename Execution::template Vector<uint64_t> invalid(nr);
    auto bad = Execution::data(invalid);
    Execution::for_each(nr,
                        [=] MARS_ENTITY_HD(size_t i)
                        {
                            bad[i] = 0;
                            if (i && r[i - 1].split.node == r[i].split.node) return;
                            size_t end = i + 1;
                            while (end < nr && r[end].split.node == r[i].split.node)
                                ++end;
                            bool exists       = false;
                            size_t body_count = 0, split_count = 0;
                            for (size_t j = i; j < end; ++j)
                            {
                                exists |= r[j].kind == 0;
                                body_count += r[j].kind == 3;
                                split_count += r[j].kind == 2;
                            }
                            bad[i] |= body_count && (exists || body_count > 1 || split_count);
                            for (size_t j = i; j < end; ++j)
                                if (r[j].kind == 2)
                                {
                                    bool known = false;
                                    for (size_t k = i; k < end; ++k)
                                    {
                                        if (r[k].kind == 1 && r[k].split.parent == r[j].split.parent) known = true;
                                        if (r[k].kind == 2 && !(r[k].split.parent == r[j].split.parent)) bad[i] = 1;
                                    }
                                    if (exists && !known) bad[i] = 1;
                                }
                        });
    require_collectively(Execution::sum(bad, nr) == 0,
                         "AMR midpoint key collides with an unrelated node or split entity", comm);
}

// Query only current coarse entities. Dropping unmatched sources bounds retained ancestry by the live mesh.
// Local traffic stays on device; shared entities rendezvous at the SFC owner of their smallest corner.
template<class Execution, class KeyType, class NodeOwner>
typename Execution::template Vector<SplitEntity<KeyType>>
match_split_entities(const typename Execution::template Vector<SplitEntity<KeyType>>& sources,
                     const typename Execution::template Vector<RefinementEntity<KeyType>>& queries,
                     NodeOwner owner,
                     MPI_Comm comm,
                     bool staged)
{
    using Packet    = split_detail::MatchPacket<KeyType>;
    using Packets   = typename Execution::template Vector<Packet>;
    const size_t nq = queries.size(), n = nq + sources.size();
    Packets packets(n);
    typename Execution::template Vector<int> destinations(n);
    auto p = Execution::data(packets);
    auto d = Execution::data(destinations);
    auto q = Execution::data(queries);
    auto s = Execution::data(sources);
    Execution::for_each(n,
                        [=] MARS_ENTITY_HD(size_t i)
                        {
                            p[i].split = i < nq ? SplitEntity<KeyType>{q[i], ~KeyType(0)} : s[i - nq];
                            d[i]       = owner(p[i].split.parent.corners[0]);
                        });
    AmrSparsePlan<Execution, Packet> plan(destinations, comm, staged);
    auto received   = plan.forward(packets);
    auto r          = Execution::data(received);
    const size_t nr = received.size();
    Execution::for_each(nr, [=] MARS_ENTITY_HD(size_t i) { r[i].slot = i; });
    Execution::sort(r, nr,
                    [] MARS_ENTITY_HD(const Packet& a, const Packet& b)
                    {
                        return a.split.parent < b.split.parent ||
                               (a.split.parent == b.split.parent && a.split.node < b.split.node);
                    });
    Packets answers(nr);
    typename Execution::template Vector<uint64_t> invalid(nr);
    auto a   = Execution::data(answers);
    auto bad = Execution::data(invalid);
    Execution::for_each(nr,
                        [=] MARS_ENTITY_HD(size_t i)
                        {
                            bad[i] = 0;
                            if (i && r[i - 1].split.parent == r[i].split.parent) return;
                            size_t end = i + 1;
                            while (end < nr && r[end].split.parent == r[i].split.parent)
                                ++end;
                            const KeyType node = r[i].split.node;
                            for (size_t j = i; j < end; ++j)
                            {
                                bad[i] |= r[j].split.node != ~KeyType(0) && r[j].split.node != node;
                                a[r[j].slot]            = r[j];
                                a[r[j].slot].split.node = node;
                            }
                        });
    require_collectively(Execution::sum(bad, nr) == 0, "Conflicting split nodes for one AMR entity", comm);
    const auto& returned = plan.reverse(answers);
    auto back            = Execution::data(returned);
    typename Execution::template Vector<SplitEntity<KeyType>> matched(nq);
    auto m = Execution::data(matched);
    Execution::for_each(nq, [=] MARS_ENTITY_HD(size_t i) { m[i] = back[i].split; });
    return compact_split_entities<Execution, KeyType>(matched);
}
} // namespace mars::amr
#undef MARS_ENTITY_HD
