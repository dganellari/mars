#pragma once

#include "mars_read_exodus_mesh.hpp"
#include "../mars_sfc_ownership.hpp"
#include <bit>
#include <exception>
#include <set>

namespace mars
{
// Prefixes only choose candidate ranks. The existing GPU lookup still decides
// exact node and block membership, including ghosts and coincident nodes.
template<class KeyType>
struct SideSetPrefix
{
    unsigned shift;
    HOST_DEVICE_FUN KeyType operator()(KeyType key) const { return key >> shift; }
};

template<class KeyType>
SideSetPrefix<KeyType> side_set_prefix(int ranks)
{
    unsigned bits  = 6;
    size_t buckets = 64;
    while (buckets < size_t(ranks) * 64)
    {
        buckets *= 2;
        ++bits;
    }
    return {unsigned(std::numeric_limits<KeyType>::digits) -
            std::min(bits, unsigned(std::numeric_limits<KeyType>::digits))};
}

// The directory is distributed by prefix, not replicated across all ranks.
// Only compact prefix metadata leaves the GPU. Face records originate in the
// host-only netCDF API and travel directly to candidate ranks before GPU lookup.
class SideSetDistribution
{
    using Word     = uint64_t;
    using Messages = std::map<int, std::vector<Word>>;
    MPI_Comm comm_ = MPI_COMM_NULL;
    int ranks_     = 1;
    std::map<Word, std::vector<int>> holders_;

    static void check_mpi(int error)
    {
        if (error != MPI_SUCCESS) throw std::runtime_error("Side-set MPI exchange failed");
    }

    std::vector<ReceivedMessage> exchange(const Messages& messages, std::vector<Word>& received, int tag) const
    {
        std::vector<int> peers, counts;
        std::vector<Word> send;
        for (const auto& [peer, values] : messages)
        {
            if (values.empty()) continue;
            if (peer < 0 || peer >= ranks_ || values.size() > size_t(std::numeric_limits<int>::max()))
                throw std::overflow_error("Side-set message exceeds MPI count or rank limits");
            peers.push_back(peer);
            counts.push_back(int(values.size()));
            send.insert(send.end(), values.begin(), values.end());
        }
        return sparseExchange<Word>(peers, counts, send.data(), received, tag, comm_);
    }

public:
    SideSetDistribution(const std::vector<Word>& prefixes, MPI_Comm comm)
    {
        check_mpi(MPI_Comm_dup(comm, &comm_));
        check_mpi(MPI_Comm_size(comm_, &ranks_));
        Messages outgoing;
        for (Word prefix : prefixes)
            outgoing[int(prefix % ranks_)].push_back(prefix);
        std::vector<Word> received;
        for (const auto& message : exchange(outgoing, received, 0))
            for (size_t i = message.offset; i < message.offset + size_t(message.count); ++i)
                holders_[received[i]].push_back(message.source);
        for (auto& [prefix, ranks] : holders_)
        {
            std::sort(ranks.begin(), ranks.end());
            ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
        }
    }

    SideSetDistribution(const SideSetDistribution&)            = delete;
    SideSetDistribution& operator=(const SideSetDistribution&) = delete;
    ~SideSetDistribution()
    {
        // Abort before collective cleanup if one rank failed while peers were exchanging.
        if (std::uncaught_exceptions() && ranks_ > 1) MPI_Abort(comm_, 1);
        if (comm_ != MPI_COMM_NULL) MPI_Comm_free(&comm_);
    }

    template<int Corners>
    ExodusSideSetBatch route(const ExodusSideSetBatch& batch, std::span<const Word> prefixes) const
    {
        static_assert(Corners == 3 || Corners == 4);
        constexpr size_t words = 1 + Corners * 6;
        if (batch.face_nodes.size() != batch.face_ids.size() * Corners || prefixes.size() != batch.nodes.size())
            throw std::invalid_argument("Side-set face IDs and corners differ");
        auto node = [](const Word* record, int corner)
        {
            const Word* data = record + 1 + 6 * corner;
            if (data[1] > Word(std::numeric_limits<int>::max()))
                throw std::overflow_error("Side-set block index exceeds int capacity");
            return ExodusSideSetNode{data[0], int(data[1]), std::bit_cast<double>(data[2]),
                                     std::bit_cast<double>(data[3]), std::bit_cast<double>(data[4])};
        };
        Messages owners;
        for (size_t f = 0; f < batch.face_ids.size(); ++f)
        {
            std::array<Word, words> record;
            record[0] = batch.face_ids[f];
            std::set<int> destinations;
            for (int c = 0; c < Corners; ++c)
            {
                const auto& n = batch.nodes.at(batch.face_nodes[f * Corners + c]);
                Word* data    = record.data() + 1 + 6 * c;
                data[0]       = n.id;
                data[1]       = Word(n.block);
                data[2]       = std::bit_cast<Word>(n.x);
                data[3]       = std::bit_cast<Word>(n.y);
                data[4]       = std::bit_cast<Word>(n.z);
                data[5]       = prefixes[batch.face_nodes[f * Corners + c]];
                destinations.insert(int(data[5] % ranks_));
            }
            for (int destination : destinations)
            {
                auto& values = owners[destination];
                values.insert(values.end(), record.begin(), record.end());
            }
        }
        std::vector<Word> requests;
        exchange(owners, requests, 1);
        if (requests.size() % words) throw std::runtime_error("Incomplete routed side-set face");
        Messages copies;
        for (size_t i = 0; i < requests.size(); i += words)
        {
            std::set<int> destinations;
            for (int c = 0; c < Corners; ++c)
            {
                auto found = holders_.find(requests[i + 1 + 6 * c + 5]);
                if (found != holders_.end()) destinations.insert(found->second.begin(), found->second.end());
            }
            for (int destination : destinations)
            {
                auto& values = copies[destination];
                values.insert(values.end(), requests.begin() + i, requests.begin() + i + words);
            }
        }
        std::vector<Word> received;
        exchange(copies, received, 2);
        if (received.size() % words) throw std::runtime_error("Incomplete delivered side-set face");
        // Different corner prefixes can send the same face here. Its file entry
        // identifies that duplicate, while distinct repeated entries remain intact.
        std::map<Word, size_t> order;
        for (size_t i = 0; i < received.size(); i += words)
            order.emplace(received[i], i);
        ExodusSideSetBatch result;
        for (auto [id, offset] : order)
        {
            if (result.nodes.size() > size_t(std::numeric_limits<unsigned>::max()) - Corners)
                throw std::overflow_error("Routed side-set batch exceeds index capacity");
            result.face_ids.push_back(id);
            for (int c = 0; c < Corners; ++c)
            {
                result.face_nodes.push_back(unsigned(result.nodes.size()));
                result.nodes.push_back(node(received.data() + offset, c));
            }
        }
        return result;
    }
};
} // namespace mars
