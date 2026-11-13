#pragma once

#include "backend/distributed/unstructured/amr/mars_amr_field_migration.hpp"

namespace amr_test
{
struct FieldTestOwner
{
    HOST_DEVICE_FUN int operator()(uint64_t) const { return 0; }
};
template<class Execution, class T>
auto field_upload(const std::vector<T>& input)
{
    typename Execution::template Vector<T> result(input.size());
    Execution::upload(input.data(), input.size(), Execution::data(result));
    return result;
}

template<class Execution>
void check_field_sources(MPI_Comm comm, bool staged)
{
    using KeyType = uint64_t;
    using Migration = mars::amr::AmrFieldMigration<Execution, KeyType, double>;
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    // The last rank is the only old owner. Other ranks have no old nodes,
    // and the new owner is rank zero, beyond any assumed local lookup.
    std::vector<KeyType> keys;
    std::vector<double> values;
    if (rank == ranks - 1)
    {
        keys = {10, 20};
        values = {3, 7};
    }
    const size_t old_count = keys.size();
    // A new interpolant at an existing hanging node must lose to its old owner.
    keys.insert(keys.end(), {20, 30, 30});
    values.insert(values.end(), {-99, 11 + double(rank), -77});
    auto pre_keys = field_upload<Execution>(keys);
    auto valid = field_upload<Execution>(std::vector<uint8_t>(keys.size(), 1));
    auto new_keys = field_upload<Execution>(std::vector<KeyType>{30, 10, 20});
    auto owned = field_upload<Execution>(std::vector<uint8_t>(3, rank == 0));
    auto owner = FieldTestOwner{};
    Migration plan(Execution::data(pre_keys), keys.size(), Execution::data(valid), Execution::data(new_keys),
                   Execution::data(owned), 3, owner, comm, staged, old_count);
    for (double scale : {1., -2., 0.})
    {
        auto scaled = values;
        for (auto& v : scaled) v *= scale;
        auto input = field_upload<Execution>(scaled);
        typename Execution::template Vector<double> output;
        plan.apply(Execution::data(input), input.size(), output);
        std::vector<double> actual(3);
        Execution::download(Execution::data(output), actual.size(), actual.data());
        const std::vector<double> expected = rank == 0 ? std::vector<double>{11 * scale, 3 * scale, 7 * scale}
                                                       : std::vector<double>(3, 0);
        mars::amr::require_collectively(actual == expected, "Field authority or reusable plan failed", comm);
    }

    Migration empty(nullptr, 0, nullptr, nullptr, nullptr, 0, owner, comm, staged, 0);
    typename Execution::template Vector<double> output;
    empty.apply(nullptr, 0, output);
    mars::amr::require_collectively(output.empty(), "Empty field migration produced data", comm);

    auto single_source = field_upload<Execution>(rank == ranks - 1 ? std::vector<KeyType>{40}
                                                                  : std::vector<KeyType>{});
    auto single_value = field_upload<Execution>(std::vector<double>(single_source.size(), 17));
    auto single_key = field_upload<Execution>(std::vector<KeyType>{40});
    auto single_owned = field_upload<Execution>(std::vector<uint8_t>{uint8_t(rank == 0)});
    Migration empty_senders(Execution::data(single_source), single_source.size(), nullptr,
                            Execution::data(single_key), Execution::data(single_owned), 1,
                            owner, comm, staged, single_source.size());
    empty_senders.apply(Execution::data(single_value), single_value.size(), output);
    double value;
    Execution::download(Execution::data(output), 1, &value);
    mars::amr::require_collectively(value == (rank == 0 ? 17 : 0), "Empty sender or receiver failed", comm);

    auto ghost_key = field_upload<Execution>(std::vector<KeyType>{40});
    auto ghost_valid = field_upload<Execution>(std::vector<uint8_t>{uint8_t(rank == ranks - 1)});
    auto ghost_value = field_upload<Execution>(std::vector<double>{rank == ranks - 1 ? 17. : -1.e8});
    Migration ghosts(Execution::data(ghost_key), 1, Execution::data(ghost_valid),
                     Execution::data(single_key), Execution::data(single_owned), 1, owner, comm, staged, 1);
    ghosts.apply(Execution::data(ghost_value), 1, output);
    Execution::download(Execution::data(output), 1, &value);
    mars::amr::require_collectively(value == (rank == 0 ? 17 : 0), "Stale ghost became a field source", comm);

    // Missing sources and duplicate old owners must fail collectively.
    for (bool duplicate : {false, true})
    {
        auto sources = field_upload<Execution>(rank == 0 && duplicate ? std::vector<KeyType>{10, 10}
                                                                      : std::vector<KeyType>{});
        auto target = field_upload<Execution>(std::vector<KeyType>{10});
        auto target_owned = field_upload<Execution>(std::vector<uint8_t>{uint8_t(rank == 0)});
        bool rejected = false;
        try
        {
            Migration invalid(Execution::data(sources), sources.size(), nullptr, Execution::data(target),
                              Execution::data(target_owned), 1, owner, comm, staged, sources.size());
        }
        catch (const std::runtime_error& e)
        {
            rejected = std::string(e.what()) == (duplicate ? "AMR node has multiple old owners"
                                                         : "AMR migration is missing an owned node");
        }
        mars::amr::require_collectively(rejected, "Invalid field sources were not rejected", comm);
    }

    // LEPP retains its stricter duplicate-value contract.
    auto strict_keys = field_upload<Execution>(std::vector<KeyType>{10, 10});
    auto target = field_upload<Execution>(std::vector<KeyType>{10});
    auto target_owned = field_upload<Execution>(std::vector<uint8_t>{uint8_t(rank == 0)});
    Migration strict(Execution::data(strict_keys), 2, nullptr, Execution::data(target),
                     Execution::data(target_owned), 1, owner, comm, staged);
    auto equal = field_upload<Execution>(std::vector<double>{5, 5});
    strict.apply(Execution::data(equal), 2, output);
    auto unequal = field_upload<Execution>(std::vector<double>{5, 6});
    bool rejected = false;
    try { strict.apply(Execution::data(unequal), 2, output); }
    catch (const std::runtime_error& e)
    {
        rejected = std::string(e.what()) == "AMR authoritative field copies disagree";
    }
    mars::amr::require_collectively(rejected, "Strict field-source check changed", comm);
}
} // namespace amr_test
