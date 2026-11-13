#pragma once
#include "backend/distributed/unstructured/amr/mars_amr_split_entities.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_generation_balance.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_tree_layout.hpp"
#include <cmath>
#include <iostream>

namespace amr_test
{
inline void topology_require(bool valid, const char* message, MPI_Comm comm)
{
    mars::amr::require_collectively(valid, message, comm);
}

template<class Execution, class T>
auto topology_upload(const std::vector<T>& host)
{
    typename Execution::template Vector<T> result(host.size());
    Execution::upload(host.data(), host.size(), Execution::data(result));
    return result;
}

template<class Execution>
void check_nonconforming_topology(MPI_Comm comm, bool staged)
{
    using namespace mars::amr;
    using Key    = uint64_t;
    using Entity = RefinementEntity<Key>;
    using Split  = SplitEntity<Key>;
    int rank, ranks;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
    // Contiguous synthetic SFC ranges, including empty ranks in the small fixture.
    auto owner  = [=] HOST_DEVICE_FUN(Key key) { return int((key / 100) % ranks); };
    auto entity = [](Key a, Key b, Key c = ~Key(0), Key d = ~Key(0))
    {
        Entity result{{a, b, c, d}};
        result.sort();
        return result;
    };
    const Entity edge = entity(20, 10), face = entity(101, 110, 120, 130), unrelated = entity(11, 21);
    for (int round = 0; round < 3; ++round)
    {
        std::vector<Split> sources;
        std::vector<Entity> queries;
        if (rank == round % ranks) sources = {{edge, 15}, {face, 115}, {edge, 15}};
        if (rank == (round + 1) % ranks) queries = {edge, face, unrelated};
        auto actual = match_split_entities<Execution, Key>(topology_upload<Execution>(sources),
                                                           topology_upload<Execution>(queries), owner, comm, staged);
        std::vector<Split> got(actual.size());
        Execution::download(Execution::data(actual), got.size(), got.data());
        bool valid = rank != (round + 1) % ranks ? got.empty() : got.size() == 2;
        if (got.size() == 2)
            valid &= got[0].parent == edge && got[0].node == 15 && got[1].parent == face && got[1].node == 115;
        topology_require(valid, "Split-entity incidence or repartition changed constraints", comm);
        auto removed = match_split_entities<Execution, Key>(
            actual, topology_upload<Execution>(std::vector<Entity>{unrelated}), owner, comm, staged);
        topology_require(removed.empty(), "Retained an obsolete coarse entity", comm);
    }
    // A nearby unrelated vertex passed the former 0.1-edge test. It has no split-parent record.
    const double distance2 = .001 * .001, old_threshold = .01;
    topology_require(distance2 < old_threshold, "Thin-mesh negative control no longer triggers old position test",
                     comm);
    auto false_match = match_split_entities<Execution, Key>(
        topology_upload<Execution>(std::vector<Split>{{unrelated, 15}}),
        topology_upload<Execution>(std::vector<Entity>{edge}), owner, comm, staged);
    topology_require(false_match.empty(), "Key coincidence invented parent incidence", comm);

    // Coincident keys are permitted only for an already recorded parent. All ranks must reject a false reuse.
    std::vector<Key> old_keys;
    std::vector<uint8_t> old_owned;
    if (rank == owner(15))
    {
        old_keys  = {15};
        old_owned = {1};
    }
    auto dk  = topology_upload<Execution>(old_keys);
    auto do_ = topology_upload<Execution>(old_owned);
    std::vector<Split> identity_sources;
    size_t prior = 0;
    if (rank == 0)
    {
        identity_sources = {{edge, 15}, {edge, 15}};
        prior            = 1;
    }
    validate_split_identity<Execution, Key>(topology_upload<Execution>(identity_sources), prior, Execution::data(dk),
                                            Execution::data(do_), old_keys.size(), owner, comm, staged);
    if (rank == 0) identity_sources = {{unrelated, 15}};
    bool rejected = false;
    try
    {
        validate_split_identity<Execution, Key>(topology_upload<Execution>(identity_sources), 0, Execution::data(dk),
                                                Execution::data(do_), old_keys.size(), owner, comm, staged);
    }
    catch (const std::runtime_error&)
    {
        rejected = true;
    }
    topology_require(rejected, "Unrelated old node reused as a midpoint", comm);

    std::vector<Split> none;
    auto body = topology_upload<Execution>(rank == 0 ? std::vector<Key>{15} : std::vector<Key>{});
    rejected  = false;
    try
    {
        validate_split_identity<Execution, Key>(topology_upload<Execution>(none), 0, Execution::data(dk),
                                                Execution::data(do_), old_keys.size(), owner, comm, staged,
                                                Execution::data(body), body.size());
    }
    catch (const std::runtime_error&)
    {
        rejected = true;
    }
    topology_require(rejected, "Hex body centre reused an old vertex key", comm);
    auto conflicts =
        topology_upload<Execution>(rank == 0 ? std::vector<Split>{{edge, 15}, {edge, 16}} : std::vector<Split>{});
    rejected = false;
    try
    {
        match_split_entities<Execution, Key>(conflicts, topology_upload<Execution>(std::vector<Entity>{edge}), owner,
                                             comm, staged);
    }
    catch (const std::runtime_error&)
    {
        rejected = true;
    }
    topology_require(rejected, "One parent acquired two midpoint identities", comm);

    // Three cells meet through two vertices; closure must cross rank boundaries twice.
    std::vector<GenerationIncident<Key>> incidents;
    std::vector<Key> generation;
    std::vector<uint8_t> marks;
    for (int cell = 0; cell < 3; ++cell)
        if (cell % ranks == rank)
        {
            size_t e = generation.size();
            generation.push_back(cell);
            marks.push_back(cell == 2);
            if (cell < 2) incidents.push_back({10, e});
            if (cell > 0) incidents.push_back({110, e});
        }
    auto d_gen   = topology_upload<Execution>(generation);
    auto d_marks = topology_upload<Execution>(marks);
    GenerationBalance<Execution, Key> balance(topology_upload<Execution>(incidents), owner, comm, staged);
    int sweeps = 0;
    while (balance.sweep(Execution::data(d_gen), Execution::data(d_marks), generation.size(), 3))
        topology_require(++sweeps < 5, "Generation closure did not settle", comm);
    Execution::download(Execution::data(d_marks), marks.size(), marks.data());
    topology_require(std::all_of(marks.begin(), marks.end(), [](auto m) { return m == 1; }),
                     "Missing generation closure", comm);
    // Initially graded cells all have generation zero; physical sizes are deliberately absent from the plan.
    std::fill(generation.begin(), generation.end(), 0);
    std::fill(marks.begin(), marks.end(), 0);
    Execution::upload(generation.data(), generation.size(), Execution::data(d_gen));
    Execution::upload(marks.data(), marks.size(), Execution::data(d_marks));
    topology_require(!balance.sweep(Execution::data(d_gen), Execution::data(d_marks), generation.size(), 3),
                     "Zero marks refined an initially graded mesh", comm);
    Execution::download(Execution::data(d_marks), marks.size(), marks.data());
    topology_require(std::all_of(marks.begin(), marks.end(), [](auto m) { return m == 0; }), "Zero marks changed",
                     comm);
    for (unsigned bucket : {1u, 64u, 1024u})
        for (double fraction : {0., .3, .5, 1.})
            for (double error : {0., .0003, .1, .2999, .3, .3001, .5, .9, 1.})
                topology_require(weighted_leaf_count(error, 1., bucket, fraction) ==
                                     bucket + unsigned(error > fraction),
                                 "Octree threshold depends on bucket size or ignores refineFraction", comm);
    if (rank == 0) std::cout << "[amr-topology] ranks=" << ranks << " staged=" << staged << " PASS\n";
}
} // namespace amr_test
