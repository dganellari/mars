#pragma once

#include "mars_amr_lepp.hpp"
#include "mars_amr_tet_boundary.hpp"
#include "mars_amr_octree_refine.hpp"

namespace mars::amr
{
template<class KeyType, class RealType>
struct TetNodeKey
{
    cstone::Box<RealType> box;
    __host__ __device__ KeyType operator()(RealType x, RealType y, RealType z) const
    {
        return cstone::sfc3D<cstone::SfcKind<KeyType>>(x, y, z, box);
    }
};

template<class KeyType, class RealType>
struct TetBisectionRefiner
{
    using Domain = ElementDomain<TetTag, RealType, KeyType, cstone::execution::Gpu>;
    using Result = typename OctreeAlignedRefine<KeyType, RealType, TetTag>::Result;

    static Result refine(const Domain& domain, const uint8_t* marks, bool host_staged = false)
    {
        int ranks;
        MPI_Comm_size(MPI_COMM_WORLD, &ranks);
        require_collectively(
            domain.numBlocks() == 1 && domain.getPeriodicAxesMask() == 0 && !domain.hasBoundaryInfo() &&
                (ranks == 1 || domain.sfcOwnership()),
            "Conforming LEPP needs nonperiodic, single-block tetrahedra with SFC ownership and explicit face tags.",
            MPI_COMM_WORLD);
        auto faces       = collect_lepp_boundary<KeyType>(domain, host_staged);
        const auto& conn = domain.getElementToNodeConnectivity();
        size_t start = domain.startIndex(), count = domain.localElementCount();
        TetLeppView<KeyType, RealType> view{
            {std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(), std::get<3>(conn).data()},
            domain.getLocalToGlobalSfcMap().data(),
            domain.getNodeX().data(),
            domain.getNodeY().data(),
            domain.getNodeZ().data(),
            marks ? marks + start : nullptr,
            count,
            domain.getNodeCount()};
        if (start)
            for (auto& pointer : view.conn)
                pointer += start;
        auto refined =
            refine_tet_lepp<LeppDeviceExecution>(view, TetNodeKey<KeyType, RealType>{domain.getBoundingBox()},
                                                 lepp_node_owner<KeyType>(domain), MPI_COMM_WORLD, host_staged, faces);
        Result result;
        result.transfer_from_parent_edges = true;
        result.oldNumNodes                = refined.old_nodes;
        result.numNodes                   = refined.keys.size();
        result.numElements                = refined.conn[0].size();
        result.numMarked                  = refined.refined_parents;
        result.numUniqueEdges             = refined.parent_edges.size();
        result.d_x                        = std::move(refined.x);
        result.d_y                        = std::move(refined.y);
        result.d_z                        = std::move(refined.z);
        for (int c = 0; c < 4; ++c)
            result.d_conn[c] = std::move(refined.conn[c]);
        result.d_edgeKeys           = std::move(refined.parent_edges);
        result.parent_birth_offsets = std::move(refined.birth_offsets);
        result.d_tagged_faces       = std::move(refined.tagged_faces);
        return result;
    }
};
} // namespace mars::amr
