#include "backend/distributed/unstructured/amr/mars_amr_tree_layout.hpp"
#include "cstone/domain/domain.hpp"

#include <array>
#include <iostream>
#include <vector>

template<class KeyType, class RealType>
void check_marks(const cstone::Domain<KeyType, double>& domain, const std::vector<KeyType>& keys)
{
    const auto global = domain.globalTree();
    const size_t first = domain.startIndex(), last = domain.endIndex();
    std::vector<RealType> errors(keys.size(), 0), local(global.numLeafNodes, 0), full(global.numLeafNodes);
    // Deliberately unrelated to key ordering: maxima must cover each entire occupied leaf.
    for (size_t e = first; e < last; ++e) errors[e] = RealType((keys[e] / 997) % 101) / 100;
    RealType maximum = 0;
    std::vector<int> owners(global.numLeafNodes, 0);
    for (int leaf = 0; leaf < global.numLeafNodes; ++leaf)
    {
        auto range = mars::amr::owned_leaf_range(global.leaves, leaf, keys.data(), first, last);
        for (size_t e = range[0]; e < range[1]; ++e) local[leaf] = std::max(local[leaf], errors[e]);
        owners[leaf] = range[0] != range[1];
        maximum = std::max(maximum, local[leaf]);
    }
    auto type = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(local.data(), full.data(), global.numLeafNodes, type, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, type, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, owners.data(), global.numLeafNodes, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    int errors_found = 0;
    for (unsigned bucket : {1u, 64u, 1024u})
        for (bool zero : {false, true})
        {
            std::vector<unsigned> reference(global.numLeafNodes), candidate(global.numLeafNodes);
            for (int leaf = 0; leaf < global.numLeafNodes; ++leaf)
            {
                reference[leaf] = bucket + unsigned(!zero && maximum > 0 && full[leaf] > RealType(.3) * maximum);
                candidate[leaf] = mars::amr::weighted_leaf_count(zero ? RealType(0) : local[leaf],
                                                                 zero ? RealType(0) : maximum, bucket);
                errors_found += owners[leaf] > 1;
            }
            for (int leaf = 0; leaf < global.numLeafNodes; ++leaf)
            {
                auto range = mars::amr::owned_leaf_range(global.leaves, leaf, keys.data(), first, last);
                if (range[0] == range[1]) continue;
                const int old_op = cstone::calculateNodeOp(global.leaves, leaf, reference.data(), bucket);
                const int new_op = cstone::calculateNodeOp(global.leaves, leaf, candidate.data(), bucket);
                errors_found += old_op != new_op;
            }
        }
    MPI_Allreduce(MPI_IN_PLACE, &errors_found, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (errors_found) throw std::runtime_error("Owned-leaf decisions differ from replicated global errors");
}

// Cornerstone must keep a unit-range leaf even when its selected elements still need refinement.
template<class KeyType>
void check_terminal_leaf()
{
    std::array<KeyType, 9> leaves;
    for (size_t i = 0; i < leaves.size(); ++i) leaves[i] = KeyType(i);
    for (unsigned bucket : {1u, 64u, 1024u})
        for (bool selected : {false, true})
        {
            std::array<unsigned, 8> counts;
            counts.fill(bucket + unsigned(selected));
            for (int i = 0; i < 8; ++i)
            {
                const auto op = cstone::calculateNodeOp(leaves.data(), i, counts.data(), bucket);
                if (op != 1 || mars::amr::mesh_leaf_mark(KeyType(1), op, counts[i], bucket) != unsigned(selected))
                    throw std::runtime_error("SFC resolution incorrectly limits mesh refinement");
                // The old op>=8 predicate misses every selected finest-resolution leaf.
                if (selected && op >= 8) throw std::runtime_error("Terminal-leaf negative control did not trigger");
            }
        }
}

// Check the actual post-sync layout against an independent per-key search.
template<class KeyType>
void check_tree(int rank, int ranks)
{
    check_terminal_leaf<KeyType>();
    cstone::Box<double> box(0, 1);
    for (int round = 0; round < 3; ++round)
    {
        std::vector<double> x, y, z, h;
        const int n = 8 + 4 * round;
        for (int i = 0; i < n * n * n; ++i)
        {
            if (i % ranks != rank) continue;
            // Repeated positions exercise several elements in one SFC bucket.
            const int j = i / 3;
            x.push_back((j % n + .5) / n);
            y.push_back((j / n % n + .5) / n);
            z.push_back((j / (n * n) + .5) / n);
            h.push_back(.08);
        }
        cstone::Domain<KeyType, double> domain(cstone::execution::cpu, rank, ranks, 32, 8, .5,
                                               MPI_COMM_WORLD, box);
        std::vector<KeyType> keys(x.size(), 0);
        std::vector<double> s1, s2, s3;
        std::vector<KeyType> generations(x.size()), identities(x.size()), s4, s5;
        for (size_t i = 0; i < x.size(); ++i)
        {
            identities[i]  = KeyType(i * ranks + rank);
            generations[i] = identities[i] % 7;
        }
        for (int sync = 0; sync < 2; ++sync)
        {
            domain.sync(keys, x, y, z, h, std::tie(generations, identities), std::tie(s1, s2, s3, s4, s5));
            domain.exchangeHalos(std::tie(generations, identities), s4, s5);
            int generation_errors = 0;
            for (size_t i = 0; i < keys.size(); ++i)
                generation_errors += generations[i] != identities[i] % 7;
            MPI_Allreduce(MPI_IN_PLACE, &generation_errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (generation_errors)
                throw std::runtime_error("Generation property changed through Cornerstone sync/halos");
            check_marks<KeyType, float>(domain, keys);
            check_marks<KeyType, double>(domain, keys);
            const auto tree = domain.octreeProperties();
            std::vector<int> levels(keys.size(), -1);
            for (int leaf = 0; leaf < tree.numLeafNodes; ++leaf)
                for (unsigned lane = 0; lane < 7; ++lane)
                    mars::amr::assign_leaf_level(tree, leaf, lane, 7, levels.data());
            int errors = tree.layout[tree.numLeafNodes] != keys.size();
            int empty = 0;
            for (int leaf = 0; leaf < tree.numLeafNodes; ++leaf)
                empty += tree.layout[leaf] == tree.layout[leaf + 1];
            for (size_t e = 0; e < keys.size(); ++e)
            {
                auto it = std::upper_bound(tree.leaves, tree.leaves + tree.numLeafNodes + 1, keys[e]);
                int leaf = int(it - tree.leaves) - 1;
                if (leaf < 0 || leaf >= tree.numLeafNodes) { ++errors; continue; }
                auto range = tree.leaves[leaf + 1] - tree.leaves[leaf];
                int expected = 0;
                for (auto full = cstone::nodeRange<KeyType>(0); full > range; full >>= 3) ++expected;
                errors += levels[e] != expected;
                errors += e < tree.layout[leaf] || e >= tree.layout[leaf + 1];
            }
            MPI_Allreduce(MPI_IN_PLACE, &errors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &empty, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (errors || !empty) throw std::runtime_error("Leaf layout/key-search parity or empty-leaf coverage failed");
            if (rank == 0) std::cout << "tree-layout bits=" << 8 * sizeof(KeyType) << " rebuild=" << round
                                     << " sync=" << sync << " ranks=" << ranks << " PASS\n";
            // Move owned records before the second sync to force new keys and a refreshed tree view.
            for (size_t e = domain.startIndex(); e < domain.endIndex(); ++e) x[e] *= .7;
        }
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try { check_tree<uint32_t>(rank, ranks); check_tree<uint64_t>(rank, ranks); }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; MPI_Abort(MPI_COMM_WORLD, 1); }
    MPI_Finalize();
}
