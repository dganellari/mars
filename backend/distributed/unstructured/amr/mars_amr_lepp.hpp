#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <span>
#include "core/mars_simplex_bisection.hpp"
#include "mars_amr_sparse_plan.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_LEPP_HD __host__ __device__
#else
#define MARS_LEPP_HD
#endif

namespace mars::amr
{
struct LeppTaggedFace
{
    uint64_t element;
    int face, tag;
};

template<class KeyType, class RealType>
struct TetLeppView
{
    const KeyType* conn[4];
    const KeyType* keys;
    const RealType* x;
    const RealType* y;
    const RealType* z;
    const uint8_t* marks;
    size_t elements, nodes;
    const KeyType* generations = nullptr;
};

template<class Execution, class KeyType, class RealType>
struct TetLeppResult
{
    template<class T>
    using Vector = typename Execution::template Vector<T>;
    std::array<Vector<KeyType>, 4> conn;
    Vector<KeyType> keys, generations;
    Vector<RealType> x, y, z;
    Vector<uint64_t> parent_edges;
    Vector<LeppTaggedFace> tagged_faces;
    std::vector<size_t> birth_offsets;
    size_t old_nodes = 0, refined_parents = 0;
};

namespace lepp_detail
{
MARS_LEPP_HD inline void edge_vertices(int edge, int& a, int& b)
{
    a = 0;
    b = 1;
    for (int i = 0, k = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j, ++k)
            if (k == edge)
            {
                a = i;
                b = j;
            }
}

template<class KeyType, class RealType>
struct Incidence
{
    bisection::Edge<KeyType> edge;
    uint64_t element;
    int slot;
};

template<class KeyType>
struct EdgePacket
{
    bisection::Edge<KeyType> edge;
    int requested = 0, terminal = 0, shared = 0;
};

template<class Execution, class KeyType>
class EdgeAgreement
{
    using Packet = EdgePacket<KeyType>;
    template<class T>
    using Vector = typename Execution::template Vector<T>;

public:
    EdgeAgreement(const Vector<int>& destinations, MPI_Comm comm, bool staged)
        : plan_(destinations, comm, staged)
        , order_(plan_.incoming_size())
        , answers_(plan_.incoming_size())
    {
    }

    const Vector<Packet>& exchange(const Vector<Packet>& packets)
    {
        auto received   = Execution::data(plan_.forward(packets));
        auto order      = Execution::data(order_);
        auto answer     = Execution::data(answers_);
        size_t incoming = order_.size();
        if (initial_)
        {
            Execution::for_each(incoming, [=] MARS_LEPP_HD(size_t i) { order[i] = i; });
            Execution::sort(order, incoming,
                            [=] MARS_LEPP_HD(uint64_t a, uint64_t b) { return received[a].edge < received[b].edge; });
            initial_ = false;
        }
        Execution::for_each(incoming,
                            [=] MARS_LEPP_HD(size_t i)
                            {
                                if (i && received[order[i - 1]].edge == received[order[i]].edge) return;
                                size_t end      = i + 1;
                                Packet combined = received[order[i]];
                                while (end < incoming && received[order[end]].edge == combined.edge)
                                {
                                    combined.terminal &= received[order[end]].terminal;
                                    combined.requested |= received[order[end]].requested;
                                    ++end;
                                }
                                // One registration per incident rank makes this an exact local-star test.
                                combined.shared = end - i > 1;
                                for (size_t j = i; j < end; ++j)
                                    answer[order[j]] = combined;
                            });
        return plan_.reverse(answers_);
    }

private:
    AmrSparsePlan<Execution, Packet> plan_;
    Vector<uint64_t> order_;
    Vector<Packet> answers_;
    bool initial_ = true;
};

template<class KeyType, class RealType>
struct IncidentSelection
{
    const Incidence<KeyType, RealType>* incidents;
    const int* selected;
    const uint64_t* cell_edges;
    const EdgePacket<KeyType>* packets;
    size_t begin;

    MARS_LEPP_HD bisection::Edge<KeyType> operator()(size_t i) const
    {
        size_t element = incidents[begin + i].element;
        return packets[cell_edges[6 * element + selected[element]]].edge;
    }
};

template<class KeyType, class RealType>
struct NodeRecord
{
    KeyType key;
    bisection::Edge<KeyType> parent;
    RealType x, y, z;
};

// Each key is stored at one coordinator. Only new births are exchanged after
// initialization; different parent edges must never acquire one node identity.
template<class Execution, class KeyType, class RealType, class NodeOwner>
class NodeRegistry
{
    using Record = NodeRecord<KeyType, RealType>;
    template<class T>
    using Vector = typename Execution::template Vector<T>;

public:
    NodeRegistry(NodeOwner owner, MPI_Comm comm, bool staged)
        : comm_(comm)
        , staged_(staged)
        , owner_(owner)
    {
    }

    void insert(const Vector<Record>& births)
    {
        Vector<int> destinations(births.size());
        auto records = Execution::data(births);
        auto dest    = Execution::data(destinations);
        auto owner   = owner_;
        Execution::for_each(births.size(), [=] MARS_LEPP_HD(size_t i) { dest[i] = owner(records[i].key); });
        AmrSparsePlan<Execution, Record> plan(destinations, comm_, staged_);
        const auto& incoming = plan.forward(births);
        size_t old           = records_.size();
        require_collectively(incoming.size() <= SIZE_MAX - old - 1, "LEPP registry exceeds index capacity", comm_);
        size_t n = old + incoming.size();
        records_.resize(n);
        records = Execution::data(records_);
        if (incoming.size())
            Execution::copy(Execution::data(incoming), incoming.size(), Execution::data(records_) + old);
        auto mutable_records = Execution::data(records_);
        Execution::sort(mutable_records, n,
                        [] MARS_LEPP_HD(const Record& a, const Record& b) { return a.key < b.key; });
        Vector<uint64_t> invalid(n), flags(n + 1), offsets(n + 1);
        auto bad    = Execution::data(invalid);
        auto flag   = Execution::data(flags);
        auto offset = Execution::data(offsets);
        Execution::for_each(n + 1,
                            [=] MARS_LEPP_HD(size_t i)
                            {
                                bool first = i < n && (i == 0 || records[i - 1].key != records[i].key);
                                flag[i]    = first;
                                if (i < n)
                                {
                                    bad[i] = i > 0 && !first &&
                                             (records[i - 1].parent != records[i].parent ||
                                              records[i - 1].x != records[i].x || records[i - 1].y != records[i].y ||
                                              records[i - 1].z != records[i].z);
                                }
                            });
        require_collectively(Execution::sum(bad, n) == 0,
                             "LEPP node identity collision or inconsistent shared coordinates", comm_);
        Execution::scan(flag, n + 1, offset);
        Vector<Record> unique(Execution::read(offset + n));
        auto target = Execution::data(unique);
        Execution::for_each(n,
                            [=] MARS_LEPP_HD(size_t i)
                            {
                                if (flag[i]) target[offset[i]] = records[i];
                            });
        Execution::synchronize();
        records_.swap(unique);
    }

private:
    MPI_Comm comm_;
    bool staged_;
    NodeOwner owner_;
    Vector<Record> records_;
};
} // namespace lepp_detail

// Re-encoding in a new box must not merge distinct coordinates before rebuild.
template<class Execution, class KeyType, class RealType, class NodeOwner>
void validate_lepp_node_keys(const KeyType* keys,
                             const RealType* x,
                             const RealType* y,
                             const RealType* z,
                             size_t nodes,
                             NodeOwner node_owner,
                             MPI_Comm comm,
                             bool staged)
{
    using Record = lepp_detail::NodeRecord<KeyType, RealType>;
    typename Execution::template Vector<Record> records(nodes);
    auto output = Execution::data(records);
    Execution::for_each(nodes, [=] MARS_LEPP_HD(size_t i)
                        { output[i] = {keys[i], bisection::Edge<KeyType>(keys[i], keys[i]), x[i], y[i], z[i]}; });
    AmrCommunicator private_comm(comm);
    lepp_detail::NodeRegistry<Execution, KeyType, RealType, NodeOwner> registry(node_owner, private_comm, staged);
    registry.insert(records);
}

// Local legacy LEPP decisions and distributed edge agreement share one
// implementation on CPU, CUDA and HIP. Cornerstone owns the subsequent rebuild.
template<class Execution, class KeyType, class RealType, class NodeKey, class NodeOwner>
TetLeppResult<Execution, KeyType, RealType>
refine_tet_lepp(TetLeppView<KeyType, RealType> input,
                NodeKey node_key,
                NodeOwner node_owner,
                MPI_Comm parent_comm,
                bool host_staged                                                       = false,
                const typename Execution::template Vector<LeppTaggedFace>& input_faces = {},
                KeyType maximum_generation = std::numeric_limits<KeyType>::max())
{
    using namespace lepp_detail;
    using Edge         = bisection::Edge<KeyType>;
    using Incident     = Incidence<KeyType, RealType>;
    using Packet       = EdgePacket<KeyType>;
    using Record       = NodeRecord<KeyType, RealType>;
    using U64Vector    = typename Execution::template Vector<uint64_t>;
    using KeyVector    = typename Execution::template Vector<KeyType>;
    using PacketVector = typename Execution::template Vector<Packet>;
    using FaceVector   = typename Execution::template Vector<LeppTaggedFace>;
    AmrCommunicator comm(parent_comm);
    require_collectively(input.nodes <= UINT32_MAX && input.elements <= SIZE_MAX / (6 * sizeof(Incident)),
                         "LEPP local mesh exceeds index capacity", comm);
    bool pointers_valid = !input.nodes || (input.keys && input.x && input.y && input.z);
    pointers_valid &=
        !input.elements || (input.marks && input.conn[0] && input.conn[1] && input.conn[2] && input.conn[3]);
    require_collectively(pointers_valid, "Missing LEPP input array", comm);

    TetLeppResult<Execution, KeyType, RealType> result;
    result.generations.resize(input.elements);
    auto initial_generation = Execution::data(result.generations);
    Execution::for_each(input.elements, [=] MARS_LEPP_HD(size_t e)
                        { initial_generation[e] = input.generations ? input.generations[e] : 0; });
    result.old_nodes = input.nodes;
    result.birth_offsets.push_back(input.nodes);
    result.keys.resize(input.nodes);
    result.x.resize(input.nodes);
    result.y.resize(input.nodes);
    result.z.resize(input.nodes);
    Execution::copy(input.keys, input.nodes, Execution::data(result.keys));
    Execution::copy(input.x, input.nodes, Execution::data(result.x));
    Execution::copy(input.y, input.nodes, Execution::data(result.y));
    Execution::copy(input.z, input.nodes, Execution::data(result.z));
    for (int c = 0; c < 4; ++c)
    {
        result.conn[c].resize(input.elements);
        Execution::copy(input.conn[c], input.elements, Execution::data(result.conn[c]));
    }
    result.tagged_faces = input_faces;
    U64Vector needed(input.elements), original(input.elements), invalid(input.elements + input_faces.size());
    auto need  = Execution::data(needed);
    auto old   = Execution::data(original);
    auto bad   = Execution::data(invalid);
    auto faces = Execution::data(input_faces);
    Execution::for_each(input.elements,
                        [=] MARS_LEPP_HD(size_t e)
                        {
                            need[e] = input.marks[e] == 1;
                            old[e]  = 1;
                            bad[e]  = 0;
                            for (int c = 0; c < 4; ++c)
                                bad[e] |= input.conn[c][e] >= input.nodes;
                        });
    Execution::for_each(input_faces.size(),
                        [=] MARS_LEPP_HD(size_t i)
                        {
                            bad[input.elements + i] = faces[i].element >= input.elements || faces[i].face < 0 ||
                                                      faces[i].face >= 4 || faces[i].tag < 0;
                        });
    require_collectively(Execution::sum(bad, invalid.size()) == 0, "Invalid LEPP connectivity or tagged face", comm);

    invalid.resize(input.nodes + input.elements);
    bad                   = Execution::data(invalid);
    KeyVector sorted_keys = result.keys;
    auto sorted           = Execution::data(sorted_keys);
    Execution::sort(sorted, input.nodes, [] MARS_LEPP_HD(KeyType a, KeyType b) { return a < b; });
    Execution::for_each(input.nodes,
                        [=] MARS_LEPP_HD(size_t i)
                        {
                            bad[i] = (i && sorted[i - 1] == sorted[i]) || !std::isfinite(input.x[i]) ||
                                     !std::isfinite(input.y[i]) || !std::isfinite(input.z[i]);
                        });
    Execution::for_each(input.elements,
                        [=] MARS_LEPP_HD(size_t e)
                        {
                            RealType a[3], b[3], c[3];
                            const RealType* coordinates[3] = {input.x, input.y, input.z};
                            for (int d = 0; d < 3; ++d)
                            {
                                a[d] = coordinates[d][input.conn[1][e]] - coordinates[d][input.conn[0][e]];
                                b[d] = coordinates[d][input.conn[2][e]] - coordinates[d][input.conn[0][e]];
                                c[d] = coordinates[d][input.conn[3][e]] - coordinates[d][input.conn[0][e]];
                            }
                            RealType det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
                                           a[2] * (b[0] * c[1] - b[1] * c[0]);
                            bad[input.nodes + e] = !std::isfinite(det) || det == RealType(0);
                        });
    require_collectively(Execution::sum(bad, invalid.size()) == 0,
                         "LEPP requires unique local node keys and finite nondegenerate tetrahedra", comm);

    NodeRegistry<Execution, KeyType, RealType, NodeOwner> registry(node_owner, comm, host_staged);
    typename Execution::template Vector<Record> births(input.nodes);
    auto birth = Execution::data(births);
    Execution::for_each(
        input.nodes, [=] MARS_LEPP_HD(size_t i)
        { birth[i] = {input.keys[i], Edge(input.keys[i], input.keys[i]), input.x[i], input.y[i], input.z[i]}; });
    registry.insert(births);

    typename Execution::template Vector<Incident> incidents;
    typename Execution::template Vector<int> longest, destinations;
    U64Vector flags, offsets, edge_offsets, cell_edges, child_counts, child_offsets, edge_counts, midpoint_offsets;
    PacketVector packets, answers;
    auto global_sum = [&](uint64_t local)
    {
        uint64_t global;
        MPI_Allreduce(&local, &global, 1, MPI_UINT64_T, MPI_SUM, comm);
        return global;
    };
    while (global_sum(Execution::sum(Execution::data(needed), needed.size())))
    {
        size_t cells = needed.size(), incidence_count = 6 * cells;
        TetLeppView<KeyType, RealType> mesh{{Execution::data(result.conn[0]), Execution::data(result.conn[1]),
                                             Execution::data(result.conn[2]), Execution::data(result.conn[3])},
                                            Execution::data(result.keys),
                                            Execution::data(result.x),
                                            Execution::data(result.y),
                                            Execution::data(result.z),
                                            nullptr,
                                            cells,
                                            result.keys.size()};
        incidents.resize(incidence_count);
        longest.resize(cells);
        invalid.resize(cells);
        auto incident = Execution::data(incidents);
        auto selected = Execution::data(longest);
        bad           = Execution::data(invalid);
        Execution::for_each(
            cells,
            [=] MARS_LEPP_HD(size_t e)
            {
                bisection::LongestEdge<KeyType, RealType> choice;
                KeyType nodes[4];
                for (int c = 0; c < 4; ++c)
                    nodes[c] = mesh.conn[c][e];
                for (int a = 0, k = 0; a < 4; ++a)
                    for (int b = a + 1; b < 4; ++b, ++k)
                    {
                        RealType dx = mesh.x[nodes[a]] - mesh.x[nodes[b]], dy = mesh.y[nodes[a]] - mesh.y[nodes[b]],
                                 dz = mesh.z[nodes[a]] - mesh.z[nodes[b]];
                        choice.consider(dx * dx + dy * dy + dz * dz, mesh.keys[nodes[a]], mesh.keys[nodes[b]], k);
                        incident[6 * e + k] = {Edge(mesh.keys[nodes[a]], mesh.keys[nodes[b]]), e, k};
                    }
                selected[e] = choice.index;
                bad[e]      = choice.index < 0 || !(choice.squared_length > 0) ||
                         choice.squared_length > std::numeric_limits<RealType>::max();
            });
        require_collectively(Execution::sum(bad, cells) == 0, "LEPP requires finite nonzero edge lengths", comm);
        Execution::sort(incident, incidence_count, [] MARS_LEPP_HD(const Incident& a, const Incident& b)
                        { return a.edge < b.edge || (a.edge == b.edge && a.element < b.element); });
        flags.resize(incidence_count + 1);
        offsets.resize(incidence_count + 1);
        auto flag   = Execution::data(flags);
        auto offset = Execution::data(offsets);
        Execution::for_each(incidence_count + 1, [=] MARS_LEPP_HD(size_t i)
                            { flag[i] = i < incidence_count && (i == 0 || incident[i - 1].edge != incident[i].edge); });
        Execution::scan(flag, incidence_count + 1, offset);
        size_t edges = Execution::read(offset + incidence_count);
        edge_offsets.resize(edges + 1);
        cell_edges.resize(incidence_count);
        packets.resize(edges);
        destinations.resize(edges);
        auto edge_offset = Execution::data(edge_offsets);
        auto cell_edge   = Execution::data(cell_edges);
        auto packet      = Execution::data(packets);
        auto destination = Execution::data(destinations);
        Execution::for_each(incidence_count,
                            [=] MARS_LEPP_HD(size_t i)
                            {
                                size_t group                                          = offset[i] + flag[i] - 1;
                                cell_edge[6 * incident[i].element + incident[i].slot] = group;
                                if (flag[i])
                                {
                                    edge_offset[group] = i;
                                    packet[group].edge = incident[i].edge;
                                }
                            });
        Execution::for_each(size_t(1), [=] MARS_LEPP_HD(size_t) { edge_offset[edges] = incidence_count; });
        need = Execution::data(needed);
        Execution::for_each(edges,
                            [=] MARS_LEPP_HD(size_t i)
                            {
                                size_t begin = edge_offset[i], end = edge_offset[i + 1];
                                IncidentSelection<KeyType, RealType> selection{incident, selected, cell_edge, packet,
                                                                               begin};
                                packet[i].terminal  = bisection::terminal_edge(packet[i].edge, end - begin, selection);
                                packet[i].requested = 0;
                                for (size_t j = begin; j < end; ++j)
                                    if (need[incident[j].element] && selected[incident[j].element] == incident[j].slot)
                                        packet[i].requested = 1;
                                destination[i] = node_owner(packet[i].edge.first);
                            });
        // Discover complete edge stars once per topology round. Interior edges
        // then stay local even when their SFC coordinator is on another rank.
        {
            EdgeAgreement<Execution, KeyType> registration(destinations, comm, host_staged);
            answers = registration.exchange(packets);
        }
        auto agreed = Execution::data(answers);
        U64Vector shared_counts(edges + 1), shared_offsets(edges + 1), changed(edges);
        auto shared_count  = Execution::data(shared_counts);
        auto shared_offset = Execution::data(shared_offsets);
        auto change        = Execution::data(changed);
        Execution::for_each(edges + 1, [=] MARS_LEPP_HD(size_t i) { shared_count[i] = i < edges && agreed[i].shared; });
        Execution::scan(shared_count, edges + 1, shared_offset);
        size_t shared_edges = Execution::read(shared_offset + edges);
        PacketVector shared_packets(shared_edges);
        typename Execution::template Vector<int> shared_destinations(shared_edges);
        auto shared_packet = Execution::data(shared_packets);
        auto shared_dest   = Execution::data(shared_destinations);
        Execution::for_each(edges,
                            [=] MARS_LEPP_HD(size_t i)
                            {
                                if (shared_count[i]) shared_dest[shared_offset[i]] = destination[i];
                            });
        EdgeAgreement<Execution, KeyType> shared_plan(shared_destinations, comm, host_staged);
        for (;;)
        {
            Execution::for_each(edges,
                                [=] MARS_LEPP_HD(size_t i)
                                {
                                    packet[i] = agreed[i];
                                    for (size_t j = edge_offset[i]; j < edge_offset[i + 1]; ++j)
                                    {
                                        size_t e = incident[j].element;
                                        if (selected[e] != incident[j].slot) continue;
                                        for (int k = 0; k < 6; ++k)
                                            packet[i].requested |= agreed[cell_edge[6 * e + k]].requested;
                                    }
                                    change[i] = packet[i].requested != agreed[i].requested;
                                });
            if (global_sum(Execution::sum(change, edges)) == 0) break;
            Execution::for_each(edges,
                                [=] MARS_LEPP_HD(size_t i)
                                {
                                    if (shared_count[i]) shared_packet[shared_offset[i]] = packet[i];
                                });
            auto shared_answer = Execution::data(shared_plan.exchange(shared_packets));
            Execution::for_each(edges, [=] MARS_LEPP_HD(size_t i)
                                { agreed[i] = shared_count[i] ? shared_answer[shared_offset[i]] : packet[i]; });
        }

        edge_counts.resize(edges + 1);
        midpoint_offsets.resize(edges + 1);
        auto edge_count = Execution::data(edge_counts);
        auto mid_offset = Execution::data(midpoint_offsets);
        Execution::for_each(edges + 1, [=] MARS_LEPP_HD(size_t i)
                            { edge_count[i] = i < edges && packet[i].requested && packet[i].terminal; });
        Execution::scan(edge_count, edges + 1, mid_offset);
        size_t new_nodes = Execution::read(mid_offset + edges), node_base = result.keys.size();
        require_collectively(new_nodes <= UINT32_MAX - node_base, "LEPP exceeds 32-bit local parent-edge indices",
                             comm);
        require_collectively(global_sum(new_nodes) != 0, "Distributed LEPP found no terminal edge", comm);
        result.keys.resize(node_base + new_nodes);
        result.x.resize(node_base + new_nodes);
        result.y.resize(node_base + new_nodes);
        result.z.resize(node_base + new_nodes);
        auto keys          = Execution::data(result.keys);
        auto x             = Execution::data(result.x);
        auto y             = Execution::data(result.y);
        auto z             = Execution::data(result.z);
        size_t prior_edges = result.parent_edges.size();
        result.parent_edges.resize(prior_edges + new_nodes);
        auto parents = Execution::data(result.parent_edges);
        births.resize(new_nodes);
        birth = Execution::data(births);
        Execution::for_each(edges,
                            [=] MARS_LEPP_HD(size_t i)
                            {
                                if (!edge_count[i]) return;
                                const auto first_incident = incident[edge_offset[i]];
                                int a, b;
                                edge_vertices(first_incident.slot, a, b);
                                uint64_t na = mesh.conn[a][first_incident.element],
                                         nb = mesh.conn[b][first_incident.element];
                                // Keep midpoint arithmetic identical across local node orderings.
                                if (keys[nb] < keys[na])
                                {
                                    auto tmp = na;
                                    na       = nb;
                                    nb       = tmp;
                                }
                                size_t index                         = node_base + mid_offset[i];
                                x[index]                             = bisection::midpoint_value(x[na], x[nb]);
                                y[index]                             = bisection::midpoint_value(y[na], y[nb]);
                                z[index]                             = bisection::midpoint_value(z[na], z[nb]);
                                keys[index]                          = node_key(x[index], y[index], z[index]);
                                parents[prior_edges + mid_offset[i]] = na | (nb << 32);
                                birth[mid_offset[i]] = {keys[index], packet[i].edge, x[index], y[index], z[index]};
                            });
        registry.insert(births);
        result.birth_offsets.push_back(node_base + new_nodes);

        child_counts.resize(cells + 1);
        child_offsets.resize(cells + 1);
        invalid.resize(cells);
        auto count        = Execution::data(child_counts);
        auto child_offset = Execution::data(child_offsets);
        old               = Execution::data(original);
        bad               = Execution::data(invalid);
        Execution::for_each(cells + 1,
                            [=] MARS_LEPP_HD(size_t e)
                            {
                                if (e == cells)
                                {
                                    count[e] = 0;
                                    return;
                                }
                                bool split = edge_count[cell_edge[6 * e + selected[e]]];
                                count[e]   = 1 + split;
                                bad[e]     = split && old[e];
                            });
        result.refined_parents += Execution::sum(bad, cells);
        auto generation = Execution::data(result.generations);
        Execution::for_each(
            cells, [=] MARS_LEPP_HD(size_t e)
            { bad[e] = generation[e] > maximum_generation || (count[e] == 2 && generation[e] >= maximum_generation); });
        require_collectively(Execution::sum(bad, cells) == 0,
                             "LEPP conformity closure exceeds the maximum refinement generation", comm);
        Execution::scan(count, cells + 1, child_offset);
        size_t new_cells = Execution::read(child_offset + cells);
        require_collectively(new_cells <= SIZE_MAX / (6 * sizeof(Incident)), "LEPP child mesh exceeds index capacity",
                             comm);
        std::array<KeyVector, 4> children;
        for (auto& column : children)
            column.resize(new_cells);
        KeyType* child[4] = {Execution::data(children[0]), Execution::data(children[1]), Execution::data(children[2]),
                             Execution::data(children[3])};
        KeyVector child_generations(new_cells);
        auto next_generation = Execution::data(child_generations);
        U64Vector next_needed(new_cells), next_original(new_cells);
        auto next_need = Execution::data(next_needed);
        auto next_old  = Execution::data(next_original);
        Execution::for_each(cells,
                            [=] MARS_LEPP_HD(size_t e)
                            {
                                KeyType nodes[4];
                                for (int c = 0; c < 4; ++c)
                                    nodes[c] = mesh.conn[c][e];
                                size_t out = child_offset[e];
                                for (size_t j = 0; j < count[e]; ++j)
                                    next_generation[out + j] = generation[e] + (count[e] == 2);
                                if (count[e] == 1)
                                {
                                    for (int c = 0; c < 4; ++c)
                                        child[c][out] = nodes[c];
                                    next_need[out] = need[e];
                                    next_old[out]  = old[e];
                                }
                                else
                                {
                                    int a, b;
                                    edge_vertices(selected[e], a, b);
                                    KeyType midpoint = KeyType(node_base + mid_offset[cell_edge[6 * e + selected[e]]]);
                                    KeyType left[4], right[4];
                                    bisection::child_connectivity<4>(nodes, b, midpoint, left);
                                    bisection::child_connectivity<4>(nodes, a, midpoint, right);
                                    for (int c = 0; c < 4; ++c)
                                    {
                                        child[c][out]     = left[c];
                                        child[c][out + 1] = right[c];
                                    }
                                    next_need[out] = next_need[out + 1] = next_old[out] = next_old[out + 1] = 0;
                                }
                            });
        size_t face_count = result.tagged_faces.size();
        flags.resize(face_count + 1);
        offsets.resize(face_count + 1);
        flag   = Execution::data(flags);
        offset = Execution::data(offsets);
        faces  = Execution::data(result.tagged_faces);
        Execution::for_each(
            face_count + 1,
            [=] MARS_LEPP_HD(size_t i)
            {
                if (i == face_count)
                {
                    flag[i] = 0;
                    return;
                }
                int a, b;
                edge_vertices(selected[faces[i].element], a, b);
                flag[i] = count[faces[i].element] == 2 && bisection::face_children(faces[i].face, a, b) == 3 ? 2 : 1;
            });
        Execution::scan(flag, face_count + 1, offset);
        FaceVector child_faces(Execution::read(offset + face_count));
        auto next_face = Execution::data(child_faces);
        Execution::for_each(face_count,
                            [=] MARS_LEPP_HD(size_t i)
                            {
                                auto face  = faces[i];
                                size_t out = offset[i], e = face.element;
                                if (count[e] == 1)
                                {
                                    next_face[out] = {child_offset[e], face.face, face.tag};
                                    return;
                                }
                                int a, b;
                                edge_vertices(selected[e], a, b);
                                unsigned mask = bisection::face_children(face.face, a, b);
                                for (int c = 0; c < 2; ++c)
                                    if (mask & (1u << c)) next_face[out++] = {child_offset[e] + c, face.face, face.tag};
                            });
        Execution::synchronize();
        for (int c = 0; c < 4; ++c)
            result.conn[c].swap(children[c]);
        result.generations.swap(child_generations);
        result.tagged_faces.swap(child_faces);
        needed.swap(next_needed);
        original.swap(next_original);
    }
    return result;
}

template<class Execution, class RealType>
typename Execution::template Vector<RealType> transfer_lepp_field(const RealType* old_values,
                                                                  size_t old_nodes,
                                                                  const uint64_t* parent_edges,
                                                                  size_t edge_count,
                                                                  std::span<const size_t> birth_offsets,
                                                                  MPI_Comm comm)
{
    bool valid = edge_count <= SIZE_MAX - old_nodes && (!old_nodes || old_values) && (!edge_count || parent_edges) &&
                 !birth_offsets.empty() && birth_offsets.front() == old_nodes &&
                 birth_offsets.back() == old_nodes + edge_count;
    for (size_t i = 1; i < birth_offsets.size(); ++i)
        valid &= birth_offsets[i - 1] <= birth_offsets[i];
    uint64_t batches = birth_offsets.size(), min_batches, max_batches;
    MPI_Allreduce(&batches, &min_batches, 1, MPI_UINT64_T, MPI_MIN, comm);
    MPI_Allreduce(&batches, &max_batches, 1, MPI_UINT64_T, MPI_MAX, comm);
    valid &= min_batches == max_batches;
    require_collectively(valid, "Invalid LEPP field parentage", comm);
    typename Execution::template Vector<RealType> result(old_nodes + edge_count);
    auto values = Execution::data(result);
    Execution::copy(old_values, old_nodes, values);
    typename Execution::template Vector<uint64_t> invalid;
    for (size_t batch = 1; batch < birth_offsets.size(); ++batch)
    {
        size_t begin = birth_offsets[batch - 1], count = birth_offsets[batch] - begin;
        invalid.resize(count);
        auto bad = Execution::data(invalid);
        Execution::for_each(count,
                            [=] MARS_LEPP_HD(size_t i)
                            {
                                uint64_t edge = parent_edges[begin - old_nodes + i];
                                size_t a = edge & UINT64_C(0xffffffff), b = edge >> 32;
                                bad[i] = a >= begin || b >= begin;
                                if (!bad[i]) values[begin + i] = bisection::midpoint_value(values[a], values[b]);
                            });
        require_collectively(Execution::sum(bad, count) == 0, "LEPP midpoint has a forward parent dependency", comm);
    }
    Execution::synchronize();
    return result;
}
} // namespace mars::amr

#undef MARS_LEPP_HD
