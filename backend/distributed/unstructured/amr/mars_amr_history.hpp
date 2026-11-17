#pragma once
#include "mars_amr_split_entities.hpp"
#include "mars_amr_generation_balance.hpp"
#include "mars_amr_tet_refine.hpp"
#include "../fem/mars_hanging_nodes.hpp"

namespace mars::amr::history_detail
{
// A named operation avoids nvcc's capture-count failure for the Tet/Hex compile-time branches.
template<class KeyType, int N>
struct SaveSplitRecords
{
    static constexpr int Q = RefinementEntities<N>::count;
    using Split            = SplitEntity<KeyType>;
    uint64_t* bad;
    size_t prior;
    const KeyType* previous;
    const uint8_t* mark;
    Split* out;
    const KeyType* old_keys;
    std::array<const KeyType*, N> conn;
    size_t first, old_nodes;
    const uint64_t* prefix;
    const uint64_t* edges;
    size_t ne;
    const KeyType* new_keys;

    __device__ void operator()(size_t i) const
    {
        bad[i] = 0;
        Split record;
        if (i < prior)
        {
            record.node = previous[6 * i];
            for (int v = 0; v < 4; ++v)
                record.parent.corners[v] = previous[6 * i + 2 + v];
        }
        else
        {
            const size_t e = (i - prior) / Q;
            const int q    = int((i - prior) % Q);
            if (!mark[e])
            {
                out[i] = record;
                return;
            }
            int v[4];
            const int nc = RefinementEntities<N>::vertices(q, v);
            for (int j = 0; j < nc; ++j)
                record.parent.corners[j] = old_keys[conn[v[j]][first + e]];
            record.parent.sort();
            size_t node = 0;
            if constexpr (N == 8)
                node = old_nodes + 19 * prefix[e] + q;
            else
            {
                uint64_t edge = EdgeKey::encode(conn[v[0]][first + e], conn[v[1]][first + e]);
                size_t lo = 0, hi = ne;
                while (lo < hi)
                {
                    size_t mid = lo + (hi - lo) / 2;
                    if (edges[mid] < edge)
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                if (lo == ne || edges[lo] != edge)
                {
                    bad[i] = 1;
                    out[i] = record;
                    return;
                }
                node = old_nodes + lo;
            }
            record.node = new_keys[node];
            for (int j = 0; j < nc; ++j)
                bad[i] |= record.node == record.parent.corners[j];
        }
        out[i] = record;
    }
};

// Keep device-lambda orchestration outside private methods: nvcc needs an addressable enclosing function.
template<class KeyType, class RealType, class ElementTag, class Domain>
std::unique_ptr<GenerationBalance<LeppDeviceExecution, KeyType>> make_generation_balance(const Domain& domain,
                                                                                         bool staged)
{
    using Execution = LeppDeviceExecution;
    using Balance   = GenerationBalance<Execution, KeyType>;
    fem::HangingNodes<KeyType, RealType> hanging(domain, staged);
    constexpr int N = ElementTag::NodesPerElement;
    std::array<const KeyType*, N> conn;
    const auto& columns = domain.getElementToNodeConnectivity();
    [&]<size_t... I>(std::index_sequence<I...>)
    { ((conn[I] = std::get<I>(columns).data()), ...); }(std::make_index_sequence<N>{});
    const size_t first = domain.startIndex(), count = domain.localElementCount();
    const auto* keys    = domain.getLocalToGlobalSfcMap().data();
    const auto* slots   = hanging.slots();
    const auto* links   = hanging.numLinks();
    const auto* masters = hanging.linkKeys();
    const size_t nh     = hanging.size();
    cstone::DeviceVector<uint64_t> counts(count + 1), offsets(count + 1);
    auto* c       = counts.data();
    auto* o       = offsets.data();
    auto linksFor = [=] __device__(KeyType slot) -> size_t
    {
        size_t lo = 0, hi = nh;
        while (lo < hi)
        {
            size_t mid = lo + (hi - lo) / 2;
            if (slots[mid] < slot)
                lo = mid + 1;
            else
                hi = mid;
        }
        return lo < nh && slots[lo] == slot ? lo : nh;
    };
    Execution::for_each(count + 1,
                        [=] __device__(size_t i)
                        {
                            c[i] = 0;
                            if (i == count) return;
                            c[i] = N;
                            for (int v = 0; v < N; ++v)
                            {
                                size_t h = linksFor(conn[v][first + i]);
                                if (h < nh) c[i] += links[h];
                            }
                        });
    Execution::scan(c, count + 1, o);
    Execution::Vector<GenerationIncident<KeyType>> incidents(Execution::read(o + count));
    auto* out = Execution::data(incidents);
    Execution::for_each(count,
                        [=] __device__(size_t i)
                        {
                            const size_t e = first + i;
                            size_t at      = o[i];
                            for (int v = 0; v < N; ++v)
                            {
                                KeyType slot = conn[v][e];
                                out[at++]    = {keys[slot], e};
                                size_t h     = linksFor(slot);
                                if (h < nh)
                                    for (int j = 0; j < links[h]; ++j)
                                        out[at++] = {masters[fem::MaxHangingLinks * h + j], e};
                            }
                        });
    auto owner = domain.sfcNodeOwner();
    return std::make_unique<Balance>(
        std::move(incidents), AmrNodeOwner<KeyType, decltype(owner)>{owner, domain.numRanks()}, MPI_COMM_WORLD, staged);
}

template<class KeyType, class ElementTag, class Domain, class Result>
LeppDeviceExecution::Vector<SplitEntity<KeyType>> save_split_entities(const Domain& domain,
                                                                      const Result& refined,
                                                                      const cstone::DeviceVector<uint8_t>& marks,
                                                                      const KeyType* new_keys,
                                                                      cstone::DeviceVector<KeyType>& child_generations,
                                                                      bool staged)
{
    using Execution = LeppDeviceExecution;
    using Split     = SplitEntity<KeyType>;
    constexpr int N = ElementTag::NodesPerElement, Q = RefinementEntities<N>::count;
    const size_t first = domain.startIndex(), count = domain.localElementCount();
    const auto* mark = marks.data() + first;
    cstone::DeviceVector<uint64_t> counts(count + 1), prefixes(count + 1);
    auto* c      = counts.data();
    auto* prefix = prefixes.data();
    Execution::for_each(count + 1, [=] __device__(size_t i) { c[i] = i < count && mark[i] == 1; });
    Execution::scan(c, count + 1, prefix);
    child_generations.resize(refined.numElements);
    auto* generation           = child_generations.data();
    const auto* old_generation = domain.refinement_generations().data();
    Execution::for_each(count,
                        [=] __device__(size_t i)
                        {
                            size_t begin = i + 7 * prefix[i];
                            for (int j = 0; j < (mark[i] == 1 ? 8 : 1); ++j)
                                generation[begin + j] = old_generation[first + i] + (mark[i] == 1);
                        });

    const auto& old    = domain.refinement_constraints();
    const size_t prior = old.size() / 6;
    Execution::Vector<Split> sources(prior + count * Q);
    auto* out            = Execution::data(sources);
    const auto* previous = old.data();
    std::array<const KeyType*, N> conn;
    const auto& columns = domain.getElementToNodeConnectivity();
    [&]<size_t... I>(std::index_sequence<I...>)
    { ((conn[I] = std::get<I>(columns).data()), ...); }(std::make_index_sequence<N>{});
    const auto* old_keys   = domain.getLocalToGlobalSfcMap().data();
    const auto* edges      = refined.d_edgeKeys.data();
    const size_t ne        = refined.d_edgeKeys.size();
    const size_t old_nodes = domain.getNodeCount();
    cstone::DeviceVector<uint64_t> invalid(sources.size());
    auto* bad = invalid.data();
    Execution::for_each(sources.size(), SaveSplitRecords<KeyType, N>{bad, prior, previous, mark, out, old_keys, conn,
                                                                     first, old_nodes, prefix, edges, ne, new_keys});
    require_collectively(Execution::sum(bad, sources.size()) == 0,
                         "AMR split entity collapsed or lost its emitted node", MPI_COMM_WORLD);
    auto compact = compact_split_entities<Execution, KeyType>(sources);
    cstone::DeviceVector<KeyType> body_keys;
    if constexpr (N == 8)
    {
        body_keys.resize(Execution::read(prefix + count));
        auto* bodies = body_keys.data();
        Execution::for_each(count,
                            [=] __device__(size_t i)
                            {
                                if (mark[i]) bodies[prefix[i]] = new_keys[old_nodes + 19 * prefix[i] + 18];
                            });
    }
    auto owner = domain.sfcNodeOwner();
    validate_split_identity<Execution, KeyType>(compact, prior, old_keys, domain.getNodeOwnershipMap().data(),
                                                old_nodes,
                                                AmrNodeOwner<KeyType, decltype(owner)>{owner, domain.numRanks()},
                                                MPI_COMM_WORLD, staged, body_keys.data(), body_keys.size());
    return compact;
}

template<class KeyType, class ElementTag, class Domain>
void restore_split_entities(Domain& domain,
                            const LeppDeviceExecution::Vector<SplitEntity<KeyType>>& sources,
                            bool staged)
{
    using Execution = LeppDeviceExecution;
    using Entity    = RefinementEntity<KeyType>;
    constexpr int N = ElementTag::NodesPerElement, Q = RefinementEntities<N>::count;
    const size_t first = domain.startIndex(), count = domain.localElementCount();
    std::array<const KeyType*, N> conn;
    const auto& columns = domain.getElementToNodeConnectivity();
    [&]<size_t... I>(std::index_sequence<I...>)
    { ((conn[I] = std::get<I>(columns).data()), ...); }(std::make_index_sequence<N>{});
    const auto* keys = domain.getLocalToGlobalSfcMap().data();
    Execution::Vector<Entity> queries(count * Q);
    auto* query = Execution::data(queries);
    Execution::for_each(queries.size(),
                        [=] __device__(size_t i)
                        {
                            Entity entity;
                            int v[4];
                            const int nc = RefinementEntities<N>::vertices(int(i % Q), v);
                            for (int j = 0; j < nc; ++j)
                                entity.corners[j] = keys[conn[v[j]][first + i / Q]];
                            entity.sort();
                            query[i] = entity;
                        });
    auto owner  = domain.sfcNodeOwner();
    auto active = match_split_entities<Execution, KeyType>(
        sources, queries, AmrNodeOwner<KeyType, decltype(owner)>{owner, domain.numRanks()}, MPI_COMM_WORLD, staged);
    cstone::DeviceVector<KeyType> records(6 * active.size());
    auto* out      = records.data();
    const auto* in = Execution::data(active);
    Execution::for_each(active.size(),
                        [=] __device__(size_t i)
                        {
                            out[6 * i]     = in[i].node;
                            out[6 * i + 1] = in[i].parent.size();
                            for (int j = 0; j < 4; ++j)
                                out[6 * i + 2 + j] = in[i].parent.corners[j];
                        });
    domain.set_refinement_constraints(std::move(records));
}

} // namespace mars::amr::history_detail
