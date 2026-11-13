#pragma once
#include "mars_amr_sparse_plan.hpp"
#include <memory>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_BALANCE_HD __host__ __device__
#else
#define MARS_BALANCE_HD
#endif

namespace mars::amr
{
template<class KeyType>
struct GenerationIncident
{
    KeyType node;
    uint64_t element;
};

namespace generation_detail
{
template<class KeyType>
struct Packet
{
    KeyType node, generation;
    uint64_t slot;
};
template<class KeyType, class Owner>
struct BuildRoutes
{
    const GenerationIncident<KeyType>* in;
    int* destinations;
    Packet<KeyType>* packets;
    Owner owner;
    MARS_BALANCE_HD void operator()(size_t i) const
    {
        destinations[i] = owner(in[i].node);
        packets[i]      = {in[i].node, 0, i};
    }
};
struct Sequence
{
    uint64_t* order;
    MARS_BALANCE_HD void operator()(size_t i) const { order[i] = i; }
};
template<class KeyType>
struct NodeLess
{
    const Packet<KeyType>* packets;
    MARS_BALANCE_HD bool operator()(uint64_t a, uint64_t b) const { return packets[a].node < packets[b].node; }
};
} // namespace generation_detail

// Each owned cell registers its vertices and the regular masters of its hanging vertices.
// The plan survives all closure sweeps; it does not depend on distance-search halos.
template<class Execution, class KeyType>
class GenerationBalance
{
    template<class T>
    using Vector = typename Execution::template Vector<T>;
    using Packet = generation_detail::Packet<KeyType>;
    Vector<GenerationIncident<KeyType>> incidents_;
    Vector<Packet> messages_, answers_;
    Vector<uint64_t> order_, changed_, invalid_;
    std::unique_ptr<AmrSparsePlan<Execution, Packet>> plan_;
    MPI_Comm comm_;

public:
    template<class Owner>
    GenerationBalance(Vector<GenerationIncident<KeyType>> incidents, Owner owner, MPI_Comm comm, bool staged)
        : incidents_(std::move(incidents))
        , messages_(incidents_.size())
        , changed_(incidents_.size())
        , invalid_(incidents_.size())
        , comm_(comm)
    {
        Vector<int> destinations(incidents_.size());
        auto in = Execution::data(incidents_);
        auto d  = Execution::data(destinations);
        auto p  = Execution::data(messages_);
        Execution::for_each(incidents_.size(), generation_detail::BuildRoutes<KeyType, Owner>{in, d, p, owner});
        plan_                = std::make_unique<AmrSparsePlan<Execution, Packet>>(destinations, comm, staged);
        const auto& received = plan_->forward(messages_);
        const size_t n       = received.size();
        answers_.resize(n);
        order_.resize(n);
        auto r     = Execution::data(received);
        auto order = Execution::data(order_);
        Execution::for_each(n, generation_detail::Sequence{order});
        Execution::sort(order, n, generation_detail::NodeLess<KeyType>{r});
    }

    bool sweep(const KeyType* generations, uint8_t* marks, size_t elements, KeyType limit)
    {
        auto in = Execution::data(incidents_);
        auto p  = Execution::data(messages_);
        Execution::for_each(incidents_.size(), [=] MARS_BALANCE_HD(size_t i)
                            { p[i].generation = generations[in[i].element] + (marks[in[i].element] == 1); });
        const auto& received = plan_->forward(messages_);
        auto r               = Execution::data(received);
        auto order           = Execution::data(order_);
        auto a               = Execution::data(answers_);
        const size_t n       = received.size();
        Execution::for_each(n,
                            [=] MARS_BALANCE_HD(size_t i)
                            {
                                if (i && r[order[i - 1]].node == r[order[i]].node) return;
                                size_t end      = i + 1;
                                KeyType maximum = r[order[i]].generation;
                                while (end < n && r[order[end]].node == r[order[i]].node)
                                {
                                    if (r[order[end]].generation > maximum) maximum = r[order[end]].generation;
                                    ++end;
                                }
                                for (size_t j = i; j < end; ++j)
                                {
                                    a[order[j]]            = r[order[j]];
                                    a[order[j]].generation = maximum;
                                }
                            });
        const auto& returned = plan_->reverse(answers_);
        auto back            = Execution::data(returned);
        // Incidents are emitted in element order. One thread writes each cell's mark.
        auto c          = Execution::data(changed_);
        auto bad        = Execution::data(invalid_);
        const size_t ni = incidents_.size();
        Execution::for_each(ni,
                            [=] MARS_BALANCE_HD(size_t i)
                            {
                                c[i] = bad[i] = 0;
                                if (i && in[i - 1].element == in[i].element) return;
                                const size_t e = in[i].element;
                                if (e >= elements)
                                {
                                    bad[i] = 1;
                                    return;
                                }
                                KeyType maximum = 0;
                                for (size_t j = i; j < ni && in[j].element == e; ++j)
                                    if (back[j].generation > maximum) maximum = back[j].generation;
                                const bool refine = maximum > generations[e] + 1;
                                bad[i]            = maximum > limit || generations[e] > limit ||
                                         (refine && generations[e] >= limit) || maximum > generations[e] + 2;
                                if (refine && marks[e] != 1)
                                {
                                    marks[e] = 1;
                                    c[i]     = 1;
                                }
                            });
        require_collectively(Execution::sum(bad, ni) == 0, "AMR generation limit or initial 2:1 balance violated",
                             comm_);
        uint64_t total = Execution::sum(c, ni), global = 0;
        MPI_Allreduce(&total, &global, 1, MPI_UINT64_T, MPI_SUM, comm_);
        return global != 0;
    }
};
} // namespace mars::amr
#undef MARS_BALANCE_HD
