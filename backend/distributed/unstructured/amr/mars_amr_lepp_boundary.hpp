#pragma once

#include "mars_amr_lepp.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_LEPP_FACE_HD __host__ __device__
#else
#define MARS_LEPP_FACE_HD
#endif

namespace mars::amr
{
template<class KeyType>
struct LeppFaceKey
{
    KeyType keys[3];
    MARS_LEPP_FACE_HD LeppFaceKey() = default;
    MARS_LEPP_FACE_HD LeppFaceKey(KeyType a, KeyType b, KeyType c)
        : keys{a, b, c}
    {
        for (int i = 1; i < 3; ++i)
            for (int j = i; j > 0 && keys[j] < keys[j - 1]; --j)
            {
                auto tmp    = keys[j];
                keys[j]     = keys[j - 1];
                keys[j - 1] = tmp;
            }
    }
    MARS_LEPP_FACE_HD bool operator==(const LeppFaceKey& other) const
    {
        return keys[0] == other.keys[0] && keys[1] == other.keys[1] && keys[2] == other.keys[2];
    }
    MARS_LEPP_FACE_HD bool operator<(const LeppFaceKey& other) const
    {
        for (int i = 0; i < 3; ++i)
        {
            if (keys[i] < other.keys[i]) return true;
            if (other.keys[i] < keys[i]) return false;
        }
        return false;
    }
};

template<class KeyType>
struct LeppFaceRecord
{
    LeppFaceKey<KeyType> face;
    int tag;
};

namespace lepp_detail
{
template<class KeyType>
struct FacePacket
{
    LeppFaceKey<KeyType> face;
    uint64_t query;
    int rank, tag;
};
struct FaceReply
{
    uint64_t query;
    int rank, tag;
};
} // namespace lepp_detail

template<class Execution, class KeyType, class NodeOwner>
typename Execution::template Vector<LeppTaggedFace>
match_lepp_faces(const typename Execution::template Vector<LeppFaceRecord<KeyType>>& sources,
                 const KeyType* keys,
                 std::array<const KeyType*, 4> conn,
                 size_t elements,
                 NodeOwner node_owner,
                 MPI_Comm parent_comm,
                 bool staged = false)
{
    using Face   = LeppFaceKey<KeyType>;
    using Packet = lepp_detail::FacePacket<KeyType>;
    using Reply  = lepp_detail::FaceReply;
    using U64    = typename Execution::template Vector<uint64_t>;
    AmrCommunicator comm(parent_comm);
    int rank;
    MPI_Comm_rank(comm, &rank);
    require_collectively(elements <= (SIZE_MAX - sources.size()) / 4, "LEPP face count exceeds index capacity", comm);
    size_t queries = 4 * elements, count = queries + sources.size();
    typename Execution::template Vector<Packet> packets(count);
    typename Execution::template Vector<int> destinations(count);
    auto packet = Execution::data(packets);
    auto dest   = Execution::data(destinations);
    auto source = Execution::data(sources);
    U64 invalid(count);
    auto bad = Execution::data(invalid);
    Execution::for_each(count,
                        [=] MARS_LEPP_FACE_HD(size_t i)
                        {
                            Packet item{};
                            item.rank  = rank;
                            item.query = i;
                            item.tag   = -1;
                            if (i < queries)
                            {
                                KeyType face_keys[3];
                                for (int c = 0, j = 0; c < 4; ++c)
                                    if (c != int(i % 4)) face_keys[j++] = keys[conn[c][i / 4]];
                                item.face = Face(face_keys[0], face_keys[1], face_keys[2]);
                            }
                            else
                            {
                                item.face = source[i - queries].face;
                                item.tag  = source[i - queries].tag;
                            }
                            bad[i] = (i >= queries && item.tag < 0) || item.face.keys[0] == item.face.keys[1] ||
                                     item.face.keys[1] == item.face.keys[2];
                            packet[i] = item;
                            dest[i]   = node_owner(item.face.keys[0]);
                        });
    require_collectively(Execution::sum(bad, count) == 0, "Invalid LEPP face key or label", comm);
    AmrSparsePlan<Execution, Packet> plan(destinations, comm, staged);
    auto received = plan.forward(packets);
    auto records  = Execution::data(received);
    size_t n      = received.size();
    Execution::sort(records, n, [] MARS_LEPP_FACE_HD(const Packet& a, const Packet& b)
                    { return a.face < b.face || (a.face == b.face && a.tag < b.tag); });
    U64 counts(n + 1), offsets(n + 1), group_end(n);
    auto counts_ptr  = Execution::data(counts);
    auto offsets_ptr = Execution::data(offsets);
    auto end_ptr     = Execution::data(group_end);
    invalid.resize(n);
    bad = Execution::data(invalid);
    Execution::for_each(n + 1,
                        [=] MARS_LEPP_FACE_HD(size_t i)
                        {
                            counts_ptr[i] = 0;
                            if (i == n) return;
                            bad[i] = 0;
                            if (i && records[i - 1].face == records[i].face) return;
                            size_t end = i + 1;
                            while (end < n && records[end].face == records[i].face)
                                ++end;
                            size_t tags = i;
                            while (tags < end && records[tags].tag < 0)
                                ++tags;
                            bad[i] = tags == i; // Every carried label must still identify a mesh face.
                            // A group leader writes only its group's queries.
                            for (size_t j = i; j < tags; ++j)
                            {
                                end_ptr[j] = end;
                            }
                        });
    // Separate initialization from group writes: device threads need no ordering.
    Execution::for_each(n,
                        [=] MARS_LEPP_FACE_HD(size_t i)
                        {
                            if (records[i].tag >= 0) return;
                            size_t lo = i + 1, hi = end_ptr[i];
                            while (lo < hi)
                            {
                                size_t mid = lo + (hi - lo) / 2;
                                if (records[mid].tag < 0)
                                    lo = mid + 1;
                                else
                                    hi = mid;
                            }
                            uint64_t unique = 0;
                            for (size_t j = lo; j < end_ptr[i]; ++j)
                                if (j == lo || records[j - 1].tag != records[j].tag) ++unique;
                            counts_ptr[i] = unique;
                        });
    require_collectively(Execution::sum(bad, n) == 0, "A carried LEPP face disappeared during rebuild", comm);
    Execution::scan(counts_ptr, n + 1, offsets_ptr);
    size_t replies_count = Execution::read(offsets_ptr + n);
    typename Execution::template Vector<Reply> replies(replies_count);
    destinations.resize(replies_count);
    dest       = Execution::data(destinations);
    auto reply = Execution::data(replies);
    Execution::for_each(n,
                        [=] MARS_LEPP_FACE_HD(size_t i)
                        {
                            if (!counts_ptr[i]) return;
                            size_t out = offsets_ptr[i];
                            for (size_t j = i + 1; j < end_ptr[i]; ++j)
                            {
                                if (records[j].tag < 0 || records[j - 1].tag == records[j].tag) continue;
                                reply[out]  = {records[i].query, records[i].rank, records[j].tag};
                                dest[out++] = records[i].rank;
                            }
                        });
    AmrSparsePlan<Execution, Reply> return_plan(destinations, comm, staged);
    const auto& returned = return_plan.forward(replies);
    typename Execution::template Vector<LeppTaggedFace> result(returned.size());
    auto output   = Execution::data(result);
    auto incoming = Execution::data(returned);
    Execution::for_each(result.size(), [=] MARS_LEPP_FACE_HD(size_t i)
                        { output[i] = {incoming[i].query / 4, int(incoming[i].query % 4), incoming[i].tag}; });
    Execution::synchronize();
    return result;
}
} // namespace mars::amr

#undef MARS_LEPP_FACE_HD
