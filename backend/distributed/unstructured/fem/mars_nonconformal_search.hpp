#pragma once

#include "mars_coarse_search.hpp"
#include "mars_nonconformal_donor_exchange.hpp"
#include "mars_nonconformal_scalar_graph.hpp"
#include "mars_nonconformal_overlap.hpp"
#include <map>
#include <stdexcept>
#if defined(__CUDACC__) || defined(__HIPCC__)
#include "mars_nonconformal_scalar_device.hpp"
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#endif

namespace mars::fem
{
namespace nonconformal_search_detail
{
// Metadata uses the existing host GhostRegistry API. Coordinates and coefficients
// never leave device storage on the device path.
template<class T, bool Device>
struct Storage
{
    using type = std::vector<T>;
};
#if defined(__CUDACC__) || defined(__HIPCC__)
template<class T>
struct Storage<T, true>
{
    using type = thrust::device_vector<T>;
};
#endif
template<class T, bool Device>
using Array = typename Storage<T, Device>::type;
template<class T>
T* pointer(std::vector<T>& a)
{
    return a.data();
}
template<class T>
const T* pointer(const std::vector<T>& a)
{
    return a.data();
}
#if defined(__CUDACC__) || defined(__HIPCC__)
template<class T>
T* pointer(thrust::device_vector<T>& a)
{
    return thrust::raw_pointer_cast(a.data());
}
template<class T>
const T* pointer(const thrust::device_vector<T>& a)
{
    return thrust::raw_pointer_cast(a.data());
}
#endif
template<bool Device, class Function>
void each(int count, Function f)
{
#if defined(__CUDACC__) || defined(__HIPCC__)
    if constexpr (Device)
    {
        thrust::for_each(thrust::device, thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(count), f);
    }
    else
#endif
    {
        for (int i = 0; i < count; ++i)
            f(i);
    }
}
template<class T>
struct Boxes
{
    const InterfaceTetFace<T>* faces;
    const int* owners;
    const uint64_t* ids;
    int rank;
    T tolerance;
    Aabb<T>* query;
    Aabb<T>* range;
    int* valid;
    MARS_FS_FN void operator()(int i) const
    {
        const auto& f = faces[i];
        nonconformal_detail::Geometry<T> g;
        valid[i] = f.id == ids[i] && (f.side == 0 || f.side == 1) &&
                   nonconformal_detail::geometry(f, g) == InterfaceStatus::valid;
        auto box = Aabb<T>::empty();
        if (valid[i])
            for (int j = 0; j < 3; ++j)
            {
                const auto& x = f.x[g.corner[j]];
                T v[]{x.x, x.y, x.z};
                for (int d = 0; d < 3; ++d)
                {
                    if (v[d] < box.lo[d]) box.lo[d] = v[d];
                    if (v[d] > box.hi[d]) box.hi[d] = v[d];
                }
            }
        for (int d = 0; d < 3; ++d)
        {
            box.lo[d] -= tolerance;
            box.hi[d] += tolerance;
        }
        query[i] = box;
        range[i] = owners[i] == rank ? box : Aabb<T>::empty();
    }
};
template<class T>
MARS_FS_FN bool faceArea(const InterfaceTetFace<T>& f, T& area)
{
    nonconformal_detail::Geometry<T> g;
    if (nonconformal_detail::geometry(f, g) != InterfaceStatus::valid) return false;
    auto n =
        nonconformal_detail::cross(sub(f.x[g.corner[1]], f.x[g.corner[0]]), sub(f.x[g.corner[2]], f.x[g.corner[0]]));
    area = std::sqrt(dot(n, n)) / T(2);
    return nonconformal_detail::finite(area) && area > T(0);
}
template<class T>
MARS_FS_FN bool intersection(const InterfaceTetFace<T>& a, const InterfaceTetFace<T>& b, T tolerance, T& area)
{
    using namespace nonconformal_detail;
    using namespace nonconformal_overlap_detail;
    Geometry<T> ga, gb;
    if (geometry(a, ga) != InterfaceStatus::valid || geometry(b, gb) != InterfaceStatus::valid) return false;
    auto origin = a.x[ga.corner[0]], u = sub(a.x[ga.corner[1]], origin);
    u      = scale(u, T(1) / std::sqrt(dot(u, u)));
    auto v = cross(ga.normal, u);
    Point2<T> pa[3], pb[3], center;
    for (int j = 0; j < 3; ++j)
    {
        auto da = sub(a.x[ga.corner[j]], origin), db = sub(b.x[gb.corner[j]], origin);
        if (std::abs(dot(db, ga.normal)) > tolerance) return false;
        pa[j] = {dot(da, u), dot(da, v)};
        pb[j] = {dot(db, u), dot(db, v)};
    }
    if (cross2(pa[0], pa[1], pa[2]) < 0)
    {
        auto t = pa[1];
        pa[1]  = pa[2];
        pa[2]  = t;
    }
    if (cross2(pb[0], pb[1], pb[2]) < 0)
    {
        auto t = pb[1];
        pb[1]  = pb[2];
        pb[2]  = t;
    }
    return overlap(pa, pb, area, center);
}
// Positive same-side overlap is checked separately: area sum alone could let
// a duplicated patch compensate for a hole elsewhere on the face.
template<class T>
struct Coverage
{
    const InterfaceTetFace<T>* faces;
    const int* offsets;
    const int* donors;
    T distanceTolerance, relativeAreaTolerance;
    int* valid;
    MARS_FS_FN void operator()(int i) const
    {
        const auto& a = faces[i];
        T area = 0, covered = 0;
        bool ok = faceArea(a, area);
        for (int k = offsets[i]; k < offsets[i + 1] && ok; ++k)
        {
            const auto& b = faces[donors[k]];
            if (a.id == b.id) continue;
            T part = 0;
            ok     = intersection(a, b, distanceTolerance, part);
            if (!ok) break;
            if (a.side == b.side)
            {
                if (part > relativeAreaTolerance * area) ok = false;
            }
            else
            {
                nonconformal_detail::Geometry<T> ga, gb;
                nonconformal_detail::geometry(a, ga);
                nonconformal_detail::geometry(b, gb);
                if (part > relativeAreaTolerance * area && T(1) + dot(ga.normal, gb.normal) > T(1e-12)) ok = false;
                covered += part;
            }
        }
        valid[i] =
            ok && nonconformal_detail::finite(covered) && std::abs(covered - area) <= relativeAreaTolerance * area;
    }
};
template<class T>
struct PairRows
{
    const InterfaceTetFace<T>* faces;
    const InterfaceFacePair* pairs;
    const int64_t* rows;
    int rowCount;
    InterfacePairingOptions<T> options;
    InterfaceStencil<T>* out;
    int* valid;
    InterfacePatch<T>* patches = nullptr;
    MARS_FS_FN void operator()(int i) const
    {
        auto pair   = pairs[i];
        auto* block = out + size_t(i) * nonconformal_pair_capacity;
        for (int j = 0; j < nonconformal_pair_capacity; ++j)
            block[j] = {};
        const auto& query = faces[pair.first];
        const auto& donor = faces[pair.second];
        valid[i]          = 1;
        if (query.id == donor.id || query.side == donor.side) return;
        bool first  = query.id < donor.id;
        int count   = 0;
        auto status = nonconformalConservativePair(
            first ? query : donor, first ? donor : query, options, block, nonconformal_pair_capacity, count,
            patches ? patches + size_t(i) * (nonconformal_pair_capacity / 2) : nullptr);
        valid[i] = status == InterfaceStatus::valid;
        if (!valid[i]) return;
        // Every query face contributes only its own rows. The reverse search
        // pair contributes the opposite rows, using the same canonical order.
        for (int j = 0; j < count; ++j)
            if ((j % 2 == 0) != first || nonconformal_detail::columnSlot(rows, 0, rowCount, block[j].row) < 0)
                block[j].status = InterfaceStatus::invalid_input;
    }
};
template<class T>
struct PatchRow
{
    InterfacePatch<T> patch;
    int side = 0;
};
template<class T>
struct ExpandPatches
{
    const InterfacePatch<T>* input;
    PatchRow<T>* output;
    MARS_FS_FN void operator()(int i) const { output[i] = {input[i / 2], i % 2}; }
};
struct KeepStencil
{
    template<class T>
    MARS_FS_FN bool operator()(const InterfaceStencil<T>& s) const
    {
        return s.status == InterfaceStatus::valid;
    }
};
} // namespace nonconformal_search_detail

// Face sets must include every incident interface face of the owned rows and
// one uniquely owned source for each face. Peers are a complete reciprocal
// candidate set from mesh/partition discovery; no global rank table is built.
// This first implementation accepts stationary coplanar P1 triangular interfaces.
template<class T, bool Device>
class NonconformalInterfacePlan
{
    using Face = InterfaceTetFace<T>;
    template<class V>
    using Array = nonconformal_search_detail::Array<V, Device>;

public:
    explicit NonconformalInterfacePlan(MPI_Comm comm)
        : exchange_(comm)
    {
        MPI_Comm_dup(comm, &comm_);
    }
    NonconformalInterfacePlan(const NonconformalInterfacePlan&)            = delete;
    NonconformalInterfacePlan& operator=(const NonconformalInterfacePlan&) = delete;
    ~NonconformalInterfacePlan()
    {
        int done = 0;
        MPI_Finalized(&done);
        if (!done) MPI_Comm_free(&comm_);
    }
    const Array<InterfaceStencil<T>>& stencils() const { return stencils_; }
    const Array<nonconformal_search_detail::PatchRow<T>>& patches() const { return patches_; }
    // Geometry is copied device-to-device on the default path. Only search
    // identities are staged for the existing host GhostRegistry setup API.
    bool build(const Array<Face>& input,
               const std::vector<uint64_t>& ids,
               const std::vector<int>& owners,
               const std::vector<int64_t>& ownedRows,
               const std::vector<int>& peers,
               T distanceTolerance,
               T relativeAreaTolerance = T(1e-10),
               bool flowPatches        = false)
    {
        try
        {
            return buildImpl(input, ids, owners, ownedRows, peers, distanceTolerance, relativeAreaTolerance,
                             flowPatches);
        }
        catch (...)
        {
            MPI_Abort(comm_, 1);
            throw;
        }
    }

private:
    // Invalid inputs return false collectively. Backend/allocation failures abort
    // the communicator so other ranks cannot remain in a pending exchange.
    bool buildImpl(const Array<Face>& input,
                   const std::vector<uint64_t>& ids,
                   const std::vector<int>& owners,
                   const std::vector<int64_t>& ownedRows,
                   const std::vector<int>& peers,
                   T distanceTolerance,
                   T relativeAreaTolerance,
                   bool flowPatches)
    {
        using namespace nonconformal_search_detail;
        stencils_.clear();
        patches_.clear();
        int rank;
        MPI_Comm_rank(comm_, &rank);
        bool ok = input.size() == ids.size() && ids.size() == owners.size() &&
                  ids.size() <= size_t(INT_MAX / sizeof(Face)) && ownedRows.size() <= size_t(INT_MAX) &&
                  std::isfinite(distanceTolerance) && distanceTolerance >= 0 && std::isfinite(relativeAreaTolerance) &&
                  relativeAreaTolerance > 0 && relativeAreaTolerance < T(.01) &&
                  std::is_sorted(ownedRows.begin(), ownedRows.end()) &&
                  std::adjacent_find(ownedRows.begin(), ownedRows.end()) == ownedRows.end();
        if (!ownedRows.empty()) ok = ok && ownedRows.front() >= 0;
        if (!allValid(ok) || !exchange_.buildFromHost(ids, owners, peers)) return false;
        faces_ = input;
        if (!forward()) return false;
        int n = int(ids.size());
        Array<Aabb<T>> query(n), range(n);
        Array<int> valid(n), owner(owners.begin(), owners.end());
        Array<uint64_t> identity(ids.begin(), ids.end());
        each<Device>(n, Boxes<T>{pointer(faces_), pointer(owner), pointer(identity), rank, distanceTolerance,
                                 pointer(query), pointer(range), pointer(valid)});
        if (!flags(valid)) return false;
        std::vector<SearchPair> found;
#if defined(__CUDACC__) || defined(__HIPCC__)
        if constexpr (Device)
        {
            Array<uint64_t> did(ids.begin(), ids.end());
            Array<SearchPair> result;
            coarse_.search(pointer(query), pointer(did), n, pointer(range), pointer(did), n, peers, comm_, result,
                           false);
            found.resize(result.size());
            thrust::copy(result.begin(), result.end(), found.begin());
        }
        else
#endif
            coarseSearchHost(query, ids, range, ids, peers, comm_, found, false);
        std::vector<uint64_t> allIds = ids;
        std::vector<int> allOwners   = owners;
        std::map<uint64_t, int> slot;
        for (int i = 0; i < n; ++i)
            slot[ids[i]] = i;
        for (const auto& pair : found)
        {
            auto q              = slot.find(pair.domain.id);
            ok                  = ok && pair.domain.proc == rank && q != slot.end() && q->second < n;
            auto [it, inserted] = slot.emplace(pair.range.id, int(allIds.size()));
            if (inserted)
            {
                allIds.push_back(pair.range.id);
                allOwners.push_back(pair.range.proc);
            }
            else
                ok = ok && allOwners[it->second] == pair.range.proc;
        }
        ok = ok && found.size() <= size_t(INT_MAX / nonconformal_pair_capacity);
        if (!allValid(ok) || !exchange_.buildFromHost(allIds, allOwners, peers)) return false;
        faces_.resize(allIds.size());
        if (!forward()) return false;
        std::vector<InterfaceFacePair> hp;
        std::vector<int> offsets(n + 1, 0), donors;
        std::vector<std::vector<int>> byQuery(n);
        for (const auto& p : found)
            byQuery[slot.at(p.domain.id)].push_back(slot.at(p.range.id));
        for (int i = 0; i < n; ++i)
        {
            auto& d = byQuery[i];
            std::sort(d.begin(), d.end());
            d.erase(std::unique(d.begin(), d.end()), d.end());
            for (int j : d)
            {
                donors.push_back(j);
                hp.push_back({i, j});
            }
            offsets[i + 1] = int(donors.size());
        }
        Array<int> off(offsets.begin(), offsets.end()), don(donors.begin(), donors.end());
        each<Device>(n, Coverage<T>{pointer(faces_), pointer(off), pointer(don), distanceTolerance,
                                    relativeAreaTolerance, pointer(valid)});
        if (!flags(valid)) return false;
        Array<InterfaceFacePair> pairs(hp.begin(), hp.end());
        Array<int64_t> rows(ownedRows.begin(), ownedRows.end());
        Array<InterfaceStencil<T>> candidate(hp.size() * nonconformal_pair_capacity);
        Array<InterfacePatch<T>> geometry(flowPatches ? candidate.size() / 2 : 0);
        valid.resize(hp.size());
        each<Device>(int(hp.size()), PairRows<T>{pointer(faces_),
                                                 pointer(pairs),
                                                 pointer(rows),
                                                 int(rows.size()),
                                                 {distanceTolerance, T(1e-12)},
                                                 pointer(candidate),
                                                 pointer(valid),
                                                 flowPatches ? pointer(geometry) : nullptr});
        if (!flags(valid)) return false;
        Array<PatchRow<T>> patchRows(flowPatches ? candidate.size() : 0);
        if (flowPatches)
        {
            each<Device>(int(candidate.size()), ExpandPatches<T>{pointer(geometry), pointer(patchRows)});
            patches_.resize(candidate.size());
        }
        stencils_.resize(candidate.size());
#if defined(__CUDACC__) || defined(__HIPCC__)
        if constexpr (Device)
        {
            auto end =
                thrust::copy_if(thrust::device, candidate.begin(), candidate.end(), stencils_.begin(), KeepStencil{});
            stencils_.resize(end - stencils_.begin());
            if (flowPatches)
            {
                auto last = thrust::copy_if(thrust::device, patchRows.begin(), patchRows.end(), candidate.begin(),
                                            patches_.begin(), KeepStencil{});
                patches_.resize(last - patches_.begin());
            }
        }
        else
#endif
        {
            auto end = std::copy_if(candidate.begin(), candidate.end(), stencils_.begin(), KeepStencil{});
            stencils_.erase(end, stencils_.end());
            if (flowPatches)
            {
                size_t kept = 0;
                for (size_t i = 0; i < candidate.size(); ++i)
                    if (KeepStencil{}(candidate[i])) patches_[kept++] = patchRows[i];
                patches_.resize(kept);
            }
        }
#if defined(__CUDACC__) || defined(__HIPCC__)
        if constexpr (Device)
        {
            bool complete = cudaGetLastError() == cudaSuccess;
            complete      = (cudaStreamSynchronize(0) == cudaSuccess) && complete;
            if (!allValid(complete))
            {
                stencils_.clear();
                return false;
            }
        }
#endif
        return true;
    }
    bool allValid(bool ok) const
    {
        int a = ok, b = 0;
        MPI_Allreduce(&a, &b, 1, MPI_INT, MPI_MIN, comm_);
        return b != 0;
    }
    bool forward()
    {
#if defined(__CUDACC__) || defined(__HIPCC__)
        if constexpr (Device)
            return exchange_.forward(nonconformal_search_detail::pointer(faces_), faces_.size());
        else
#endif
            return exchange_.forwardHost(faces_);
    }
    bool flags(const Array<int>& valid) const
    {
        bool ok = true;
#if defined(__CUDACC__) || defined(__HIPCC__)
        if constexpr (Device)
        {
            ok = cudaGetLastError() == cudaSuccess;
            ok = (cudaStreamSynchronize(0) == cudaSuccess) && ok;
            if (ok) ok = thrust::count(valid.begin(), valid.end(), 0) == 0;
        }
        else
#endif
            ok = std::find(valid.begin(), valid.end(), 0) == valid.end();
        return allValid(ok);
    }
    MPI_Comm comm_;
    NonconformalDonorExchange<T> exchange_;
    Array<Face> faces_;
    Array<InterfaceStencil<T>> stencils_;
    Array<nonconformal_search_detail::PatchRow<T>> patches_;
#if defined(__CUDACC__) || defined(__HIPCC__)
    CoarseSearch<T> coarse_;
#endif
};
template<class T>
using NonconformalInterfaceHost = NonconformalInterfacePlan<T, false>;
#if defined(__CUDACC__) || defined(__HIPCC__)
template<class T>
using NonconformalInterface = NonconformalInterfacePlan<T, true>;
#endif
} // namespace mars::fem
