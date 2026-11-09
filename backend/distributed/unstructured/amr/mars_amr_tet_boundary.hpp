#pragma once

#include "mars_amr_lepp_boundary.hpp"

namespace mars::amr
{
template<class KeyType, class Domain>
auto lepp_node_owner(const Domain& domain)
{
    auto owner = domain.sfcNodeOwner();
    return AmrNodeOwner<KeyType, decltype(owner)>{owner, owner.numRanks};
}

inline void check_lepp_side_names(const std::vector<std::string>& names, MPI_Comm comm)
{
    std::string joined;
    bool valid = names.size() <= size_t(INT_MAX);
    for (const auto& name : names)
    {
        valid &= name.find('\0') == std::string::npos;
        joined += name;
        joined += '\0';
    }
    require_collectively(valid && joined.size() <= size_t(INT_MAX), "Invalid LEPP side-set names", comm);
    int rank;
    MPI_Comm_rank(comm, &rank);
    int bytes = int(joined.size());
    MPI_Bcast(&bytes, 1, MPI_INT, 0, comm);
    std::string reference = rank == 0 ? joined : std::string(bytes, '\0');
    MPI_Bcast(reference.data(), bytes, MPI_CHAR, 0, comm);
    require_collectively(reference == joined, "LEPP side-set names differ between ranks", comm);
}

template<class KeyType, class Domain>
LeppDeviceExecution::Vector<LeppTaggedFace> collect_lepp_boundary(const Domain& domain, bool staged)
{
    using Execution = LeppDeviceExecution;
    using Record    = LeppFaceRecord<KeyType>;
    auto names      = domain.sideSetNames();
    check_lepp_side_names(names, MPI_COMM_WORLD);
    if (names.empty()) return {};
    bool valid   = true;
    size_t count = 0;
    for (const auto& name : names)
    {
        valid &= domain.has_side_set_face_keys(name) && domain.side_set_face_keys(name).size() % 3 == 0;
        count += domain.side_set_face_keys(name).size() / 3;
    }
    require_collectively(valid, "Conforming AMR requires face labels; node-only side sets have no face ancestry",
                         MPI_COMM_WORLD);
    Execution::Vector<Record> sources(count);
    auto records  = Execution::data(sources);
    size_t offset = 0;
    for (size_t tag = 0; tag < names.size(); ++tag)
    {
        const auto& faces = domain.side_set_face_keys(names[tag]);
        auto keys         = faces.data();
        size_t n          = faces.size() / 3;
        Execution::for_each(
            n, [=] __host__ __device__(size_t i)
            { records[offset + i] = {LeppFaceKey<KeyType>(keys[3 * i], keys[3 * i + 1], keys[3 * i + 2]), int(tag)}; });
        offset += n;
    }
    const auto& conn = domain.getElementToNodeConnectivity();
    size_t start     = domain.startIndex();
    std::array<const KeyType*, 4> pointers{std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                                           std::get<3>(conn).data()};
    if (start)
        for (auto& pointer : pointers)
            pointer += start;
    return match_lepp_faces<Execution, KeyType>(sources, domain.getLocalToGlobalSfcMap().data(), pointers,
                                                domain.localElementCount(), lepp_node_owner<KeyType>(domain),
                                                MPI_COMM_WORLD, staged);
}

// Save local node indices before connectivity moves into cornerstone. They
// index the saved pre-rebuild coordinates, whose keys use the rebuilt box.
template<class KeyType, class Result>
LeppDeviceExecution::Vector<LeppFaceRecord<KeyType>> save_lepp_boundary(const Result& refined)
{
    using Execution = LeppDeviceExecution;
    Execution::Vector<LeppFaceRecord<KeyType>> saved(refined.d_tagged_faces.size());
    auto faces  = Execution::data(refined.d_tagged_faces);
    auto output = Execution::data(saved);
    std::array<const KeyType*, 4> conn{refined.d_conn[0].data(), refined.d_conn[1].data(), refined.d_conn[2].data(),
                                       refined.d_conn[3].data()};
    Execution::for_each(saved.size(),
                        [=] __host__ __device__(size_t i)
                        {
                            KeyType nodes[3];
                            for (int c = 0, j = 0; c < 4; ++c)
                                if (c != faces[i].face) nodes[j++] = conn[c][faces[i].element];
                            output[i] = {LeppFaceKey<KeyType>(nodes[0], nodes[1], nodes[2]), faces[i].tag};
                        });
    Execution::synchronize();
    return saved;
}

template<class KeyType, class Domain>
void restore_lepp_boundary(Domain& domain,
                           const std::vector<std::string>& names,
                           LeppDeviceExecution::Vector<LeppFaceRecord<KeyType>>& saved,
                           const KeyType* pre_keys,
                           bool staged)
{
    using Execution = LeppDeviceExecution;
    if (names.empty()) return;
    auto records = Execution::data(saved);
    Execution::for_each(saved.size(),
                        [=] __host__ __device__(size_t i)
                        {
                            auto face       = records[i].face;
                            records[i].face = LeppFaceKey<KeyType>(pre_keys[face.keys[0]], pre_keys[face.keys[1]],
                                                                   pre_keys[face.keys[2]]);
                        });
    const auto& conn = domain.getElementToNodeConnectivity();
    std::array<const KeyType*, 4> pointers{std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                                           std::get<3>(conn).data()};
    auto matched = match_lepp_faces<Execution, KeyType>(saved, domain.getLocalToGlobalSfcMap().data(), pointers,
                                                        domain.getElementCount(), lepp_node_owner<KeyType>(domain),
                                                        MPI_COMM_WORLD, staged);
    auto faces   = Execution::data(matched);
    auto keys    = domain.getLocalToGlobalSfcMap().data();
    size_t n     = matched.size();
    Execution::Vector<uint64_t> counts(n + 1), offsets(n + 1);
    auto count  = Execution::data(counts);
    auto offset = Execution::data(offsets);
    for (size_t tag = 0; tag < names.size(); ++tag)
    {
        Execution::for_each(n + 1, [=] __host__ __device__(size_t i) { count[i] = i < n && faces[i].tag == int(tag); });
        Execution::scan(count, n + 1, offset);
        Execution::Vector<KeyType> face_keys(3 * Execution::read(offset + n));
        auto output = Execution::data(face_keys);
        Execution::for_each(n,
                            [=] __host__ __device__(size_t i)
                            {
                                if (!count[i]) return;
                                for (int c = 0, j = 0; c < 4; ++c)
                                    if (c != faces[i].face)
                                        output[3 * offset[i] + j++] = keys[pointers[c][faces[i].element]];
                            });
        domain.set_side_set_from_face_keys(names[tag], std::move(face_keys));
    }
    Execution::synchronize();
}
} // namespace mars::amr
