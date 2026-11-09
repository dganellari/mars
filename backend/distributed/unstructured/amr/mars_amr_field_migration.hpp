#pragma once

#include <memory>
#include "mars_amr_sparse_plan.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_AMR_FIELD_HD __host__ __device__
#else
#define MARS_AMR_FIELD_HD
#endif

namespace mars::amr
{
template<class KeyType, class RealType>
struct AmrFieldRecord
{
    KeyType key;
    RealType value;
};

namespace amr_detail
{
struct MarkFieldSources
{
    uint64_t* flags;
    const uint8_t* valid;
    size_t count;
    MARS_AMR_FIELD_HD void operator()(size_t i) const { flags[i] = i < count && (!valid || valid[i]); }
};
template<class KeyType, class RealType, class Owner>
struct PackFieldSources
{
    const uint64_t* flags;
    const uint64_t* offsets;
    const KeyType* keys;
    uint64_t* indices;
    AmrFieldRecord<KeyType, RealType>* packets;
    int* destinations;
    Owner owner;
    MARS_AMR_FIELD_HD void operator()(size_t i) const
    {
        if (!flags[i]) return;
        size_t j        = offsets[i];
        indices[j]      = i;
        packets[j]      = {keys[i], RealType(0)};
        destinations[j] = owner(keys[i]);
    }
};
struct InitializeFieldOrder
{
    uint64_t* order;
    MARS_AMR_FIELD_HD void operator()(size_t i) const { order[i] = i; }
};
template<class KeyType, class RealType>
struct FieldKeyLess
{
    const AmrFieldRecord<KeyType, RealType>* records;
    MARS_AMR_FIELD_HD bool operator()(uint64_t a, uint64_t b) const { return records[a].key < records[b].key; }
};
template<class KeyType, class RealType>
struct FindFieldSources
{
    const AmrFieldRecord<KeyType, RealType>* records;
    const uint64_t* order;
    const KeyType* keys;
    const uint8_t* owned;
    uint64_t* lookup;
    uint64_t* invalid;
    size_t count;
    MARS_AMR_FIELD_HD void operator()(size_t i) const
    {
        lookup[i]  = count;
        invalid[i] = 0;
        if (!owned[i]) return;
        size_t lo = 0, hi = count;
        while (lo < hi)
        {
            size_t mid = lo + (hi - lo) / 2;
            if (records[order[mid]].key < keys[i])
                lo = mid + 1;
            else
                hi = mid;
        }
        if (lo == count || records[order[lo]].key != keys[i])
            invalid[i] = 1;
        else
            lookup[i] = order[lo];
    }
};
} // namespace amr_detail

// The node layout changes once per rebuild, while every field uses the same
// sources, destinations and owned-node lookup. Preserve that plan across fields.
template<class Execution, class KeyType, class RealType>
class AmrFieldMigration
{
    template<class T>
    using Vector = typename Execution::template Vector<T>;
    using Record = AmrFieldRecord<KeyType, RealType>;

public:
    template<class Owner>
    AmrFieldMigration(const KeyType* pre_keys,
                      size_t pre_count,
                      const uint8_t* source_valid,
                      const KeyType* new_keys,
                      const uint8_t* owned,
                      size_t new_count,
                      Owner owner,
                      MPI_Comm comm,
                      bool staged = false)
        : comm_(comm)
        , input_count_(pre_count)
        , lookup_(new_count)
    {
        require_collectively(pre_count < SIZE_MAX && (!pre_count || pre_keys) && (!new_count || (new_keys && owned)),
                             "Invalid AMR migration arrays", comm_);
        Vector<uint64_t> flags(pre_count + 1), offsets(pre_count + 1);
        auto flag   = Execution::data(flags);
        auto offset = Execution::data(offsets);
        Execution::for_each(pre_count + 1, amr_detail::MarkFieldSources{flag, source_valid, pre_count});
        Execution::scan(flag, pre_count + 1, offset);
        size_t count = Execution::read(offset + pre_count);
        source_indices_.resize(count);
        packets_.resize(count);
        Vector<int> destinations(count);
        auto index       = Execution::data(source_indices_);
        auto packets     = Execution::data(packets_);
        auto destination = Execution::data(destinations);
        Execution::for_each(pre_count, amr_detail::PackFieldSources<KeyType, RealType, Owner>{
                                           flag, offset, pre_keys, index, packets, destination, owner});
        plan_                = std::make_unique<AmrSparsePlan<Execution, Record>>(destinations, comm_, staged);
        const auto& received = plan_->forward(packets_);
        auto records         = Execution::data(received);
        size_t n             = received.size();
        received_order_.resize(n);
        auto order = Execution::data(received_order_);
        Execution::for_each(n, amr_detail::InitializeFieldOrder{order});
        Execution::sort(order, n, amr_detail::FieldKeyLess<KeyType, RealType>{records});
        auto lookup = Execution::data(lookup_);
        invalid_.resize(new_count);
        auto invalid = Execution::data(invalid_);
        Execution::for_each(new_count, amr_detail::FindFieldSources<KeyType, RealType>{records, order, new_keys, owned,
                                                                                       lookup, invalid, n});
        require_collectively(Execution::sum(invalid, new_count) == 0, "AMR migration is missing an owned node", comm_);
    }

    void apply(const RealType* values, size_t count, Vector<RealType>& field)
    {
        require_collectively(count == input_count_ && (!count || values),
                             "AMR field size does not match its migration plan", comm_);
        auto packets = Execution::data(packets_);
        auto indices = Execution::data(source_indices_);
        Execution::for_each(packets_.size(),
                            [=] MARS_AMR_FIELD_HD(size_t i) { packets[i].value = values[indices[i]]; });
        const auto& received = plan_->forward(packets_);
        auto records         = Execution::data(received);
        auto order           = Execution::data(received_order_);
        size_t n             = received.size();
        invalid_.resize(n);
        auto invalid = Execution::data(invalid_);
        Execution::for_each(n,
                            [=] MARS_AMR_FIELD_HD(size_t i)
                            {
                                invalid[i] = i > 0 && records[order[i - 1]].key == records[order[i]].key &&
                                             records[order[i - 1]].value != records[order[i]].value;
                            });
        require_collectively(Execution::sum(invalid, n) == 0, "AMR authoritative field copies disagree", comm_);
        field.resize(lookup_.size());
        auto output = Execution::data(field);
        auto lookup = Execution::data(lookup_);
        Execution::for_each(field.size(), [=] MARS_AMR_FIELD_HD(size_t i)
                            { output[i] = lookup[i] == n ? RealType(0) : records[lookup[i]].value; });
        Execution::synchronize();
    }

private:
    AmrCommunicator comm_;
    size_t input_count_;
    Vector<uint64_t> source_indices_, received_order_, lookup_, invalid_;
    Vector<Record> packets_;
    std::unique_ptr<AmrSparsePlan<Execution, Record>> plan_;
};

#if defined(__CUDACC__) || defined(__HIPCC__)
template<class Domain, class KeyVector, class RealVector>
void migrate_amr_field(Domain& domain,
                       const KeyVector& pre_keys,
                       const RealVector& pre_values,
                       RealVector& field,
                       const uint8_t* source_valid = nullptr)
{
    using KeyType  = typename KeyVector::value_type;
    using RealType = typename RealVector::value_type;
    int ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    require_collectively(ranks == 1 || domain.sfcOwnership(), "AMR migration requires SFC ownership", MPI_COMM_WORLD);
    auto owner = domain.sfcNodeOwner();
    AmrFieldMigration<LeppDeviceExecution, KeyType, RealType> migration(
        pre_keys.data(), pre_keys.size(), source_valid, domain.getLocalToGlobalSfcMap().data(),
        domain.getNodeOwnershipMap().data(), domain.getNodeCount(),
        AmrNodeOwner<KeyType, decltype(owner)>{owner, ranks}, MPI_COMM_WORLD);
    migration.apply(pre_values.data(), pre_values.size(), field);
    domain.exchangeNodeHalo(field);
}
#endif
} // namespace mars::amr

#undef MARS_AMR_FIELD_HD
