#ifndef MESHTESTTOOLS_H
#define MESHTESTTOOLS_H

#include "asmBase.h"
#include <stk_io/FillMesh.hpp>
#include <stk_mesh/base/ForEachEntity.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>

using MeshToolsExecSpace =
    stk::mesh::NgpMesh::MeshExecSpace; // default exec space. Priority: Cuda >
                                       // OMP > Serial

template <typename scalar>
void set_bucket_id_field_on_host(const stk::mesh::BulkData& bulkData,
                                 stk::mesh::Field<scalar>& stkField)
{
    std::cout << "set bucket field on host for entity rank "
              << stkField.entity_rank() << " ..." << std::flush;

    const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
    const stk::mesh::BucketVector& buckets = bulkData.get_buckets(
        stkField.entity_rank(), metaData.locally_owned_part());

    stkField.sync_to_host();

    for (const stk::mesh::Bucket* bptr : buckets)
    {
        unsigned bucket_id = bptr->bucket_id();
        for (stk::mesh::Entity elem : *bptr)
        {
            *stk::mesh::field_data(stkField, elem) = bucket_id;
        }
    }

    stkField.modify_on_host();

    std::cout << " DONE" << std::endl;
}

template <typename scalar>
void set_field_on_host(const stk::mesh::BulkData& bulkData,
                       const std::string& fieldName,
                       const stk::mesh::EntityRank& entityRank,
                       const scalar& fieldVal,
                       const int verbose = 0)
{
    if (verbose > 0)
    {
        std::cout << "set " << fieldName << " field on host ..." << std::flush;
    }
    const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();

    const stk::mesh::Field<scalar>& stkField =
        *metaData.get_field<scalar>(entityRank, fieldName);

    stk::mesh::Selector select(stkField);
    stkField.sync_to_host();

    stk::mesh::for_each_entity_run(
        bulkData,
        stkField.entity_rank(),
        select,
        [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity entity)
    { *stk::mesh::field_data(stkField, entity) = fieldVal; });

    stkField.modify_on_host();

    if (verbose > 0)
    {
        std::cout << " DONE" << std::endl;
    }
}

template <typename scalar>
void fill_phi_gradPhi(const stk::mesh::BulkData& bulkData,
                      const int verbose = 0)
{
    if (verbose > 0)
    {
        std::cout << "set phi field on host ..." << std::flush;
    }

    const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();

    const stk::mesh::Field<scalar>& phiField =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "phi");
    const stk::mesh::Field<scalar>& gradPhiField =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "phiGrad");
    const stk::mesh::Field<scalar>& nodeCoords =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "coordinates");

    stk::mesh::Selector select(phiField);
    phiField.sync_to_host();

    // fill phi
    stk::mesh::for_each_entity_run(
        bulkData,
        phiField.entity_rank(),
        select,
        [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity entity)
    {
        scalar x = stk::mesh::field_data(nodeCoords, entity)[0];
        scalar y = stk::mesh::field_data(nodeCoords, entity)[1];
        scalar z = stk::mesh::field_data(nodeCoords, entity)[2];
        scalar fieldVal = sin(x) + 3 * cos(y) + 4 * sin(5 * x * y * z);
        *stk::mesh::field_data(phiField, entity) = fieldVal;
    });

    if (verbose > 0)
    {
        std::cout << " DONE" << std::endl;
    }

    stk::mesh::Selector select2(gradPhiField);
    if (verbose > 0)
    {
        std::cout << "set gradPhi field on host ..." << std::flush;
    }
    gradPhiField.sync_to_host();

    // fill gradPhi
    stk::mesh::for_each_entity_run(
        bulkData,
        gradPhiField.entity_rank(),
        select2,
        [&](const stk::mesh::BulkData& bulk, stk::mesh::Entity entity)
    {
        scalar x = stk::mesh::field_data(nodeCoords, entity)[0];
        scalar y = stk::mesh::field_data(nodeCoords, entity)[1];
        scalar z = stk::mesh::field_data(nodeCoords, entity)[2];

        scalar ddx = cos(x) + 20 * y * z * cos(5 * x * y * z);
        scalar ddy = 20 * x * z * cos(5 * x * y * z) - 3 * sin(y);
        scalar ddz = 20 * x * y * cos(5 * x * y * z);
        stk::mesh::field_data(gradPhiField, entity)[0] = ddx;
        stk::mesh::field_data(gradPhiField, entity)[1] = ddy;
        stk::mesh::field_data(gradPhiField, entity)[2] = ddz;
    });

    phiField.modify_on_host();
    gradPhiField.modify_on_host();

    if (verbose > 0)
    {
        std::cout << " DONE" << std::endl;
    }
}

template <typename scalar>
void set_field_on_device(const stk::mesh::NgpMesh& ngpMesh,
                         stk::mesh::NgpField<scalar>& ngpField,
                         scalar fieldVal)
{
    std::cout << "set field on device ..." << std::flush;
    stk::mesh::Selector select(*ngpField.get_field_base());

    ngpField.sync_to_device();

    stk::mesh::for_each_entity_run(
        ngpMesh,
        ngpField.get_rank(),
        select,
        KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entityIndex) {
            ngpField(entityIndex, 0) = fieldVal;
        });

    ngpField.modify_on_device();
    std::cout << " DONE" << std::endl;
}

template <typename scalar>
void set_bucketOrdinal_on_device(const stk::mesh::NgpMesh& ngpMesh,
                                 stk::mesh::NgpField<scalar>& ngpField)
{
    std::cout << "set bucket field on device for entity rank "
              << ngpField.get_rank() << std::flush;
    stk::mesh::Selector select(*ngpField.get_field_base());

    ngpField.sync_to_device();

    stk::mesh::for_each_entity_run(
        ngpMesh,
        ngpField.get_rank(),
        select,
        KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& entityIndex) {
            ngpField(entityIndex, 0) = entityIndex.bucket_ord;
        });

    ngpField.modify_on_device();
    std::cout << " DONE" << std::endl;
}

template <typename scalar>
void add_scalar_field(const std::string& field_name,
                      stk::mesh::MetaData& meta,
                      stk::mesh::EntityRank rank)
{
    stk::mesh::Field<scalar>& field = meta.declare_field<scalar>(rank, field_name);

    // initialize field without value (nullptr instead of value)
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), nullptr);
}

template <typename scalar, unsigned N>
void add_scalar_N_field(const std::string& field_name,
                        stk::mesh::MetaData& meta,
                        stk::mesh::EntityRank rank)
{
    constexpr unsigned lengthPerEntity = N;

    stk::mesh::Field<scalar>& field = meta.declare_field<scalar>(rank, field_name);

    // initialize field without value (nullptr instead of value)
    stk::mesh::put_field_on_mesh(
        field, meta.universal_part(), lengthPerEntity, nullptr);
}

template <typename scalar, unsigned M, unsigned N>
void add_scalar_MxN_field(const std::string& field_name,
                          stk::mesh::MetaData& meta,
                          stk::mesh::EntityRank rank)
{
    stk::mesh::Field<scalar>& field = meta.declare_field<scalar>(rank, field_name);

    // initialize field without value (nullptr instead of value)
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), M, N, nullptr);
}

template <typename scalar>
stk::mesh::Field<scalar>& add_vector_field(const std::string& field_name,
                                      stk::mesh::MetaData& meta,
                                      stk::mesh::EntityRank rank)
{
    constexpr unsigned vectorFieldLengthPerEntity = 3;

    stk::mesh::Field<scalar>& field = meta.declare_field<scalar>(rank, field_name);

    // initialize field without value (nullptr instead of value)
    stk::mesh::put_field_on_mesh(
        field, meta.universal_part(), vectorFieldLengthPerEntity, nullptr);
    stk::io::set_field_output_type(field, stk::io::FieldOutputType::VECTOR_3D);

    return field;
}

void read_mesh_from_file(const std::string& filename,
                         const stk::ParallelMachine& pm,
                         std::unique_ptr<stk::mesh::BulkData>& bulkData,
                         const std::unique_ptr<stk::io::StkMeshIoBroker>& stkIo,
                         const int verbose = 0)
{
    if (verbose > 0)
    {
        std::cout << "Reading Mesh ... " << std::flush;
    }

    stk::mesh::MeshBuilder builder(pm);
    builder.set_aura_option(stk::mesh::BulkData::AUTO_AURA);
    bulkData = builder.create();

    // some field options
    bulkData->mesh_meta_data().use_simple_fields();
    bulkData->mesh_meta_data().enable_late_fields();

    // read mesh from file
    stk::io::fill_mesh_preexisting(
        *stkIo, filename, *bulkData, stk::io::READ_MESH);

    if (verbose > 0)
    {
        std::cout << "    Creating Edges ... " << std::flush;
    }
    // // just for areaV
    stk::mesh::create_edges(*bulkData.get());

    if (verbose > 0)
    {
        std::cout << "DONE!" << std::endl;
    }

    unsigned numNodes =
        stk::mesh::count_entities(*bulkData,
                                  stk::topology::NODE_RANK,
                                  bulkData->mesh_meta_data().universal_part());
    unsigned numEdges =
        stk::mesh::count_entities(*bulkData,
                                  stk::topology::EDGE_RANK,
                                  bulkData->mesh_meta_data().universal_part());
    unsigned numFaces =
        stk::mesh::count_entities(*bulkData,
                                  stk::topology::FACE_RANK,
                                  bulkData->mesh_meta_data().universal_part());
    unsigned numElems =
        stk::mesh::count_entities(*bulkData,
                                  stk::topology::ELEMENT_RANK,
                                  bulkData->mesh_meta_data().universal_part());
    std::cout << "Mesh Size (Nodes):    " << numNodes << std::endl;
    std::cout << "Mesh Size (Edges):    " << numEdges << std::endl;
    std::cout << "Mesh Size (Faces):    " << numFaces << std::endl;
    std::cout << "Mesh Size (Elements): " << numElems << std::endl;
}

template <typename scalar, typename label, size_t BLOCKSIZE, size_t SPATIAL_DIM>
void calculateEdgeAreaVectors(const stk::mesh::BulkData& bulkData)
{
    using Context = AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::Context;
    using ScratchViewIndex =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::
            ScratchViewIndex;
    using ScratchViewEntity =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::
            ScratchViewEntity;
    using ScratchView2D =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchView2D;
    using ScratchView3D =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchView3D;
    using TeamHandleType =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::TeamHandleType;

    const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();

    const stk::mesh::Field<scalar>& areaVEdgeField =
        *metaData.get_field<scalar>(stk::topology::EDGE_RANK, "areaVEdge");
    const stk::mesh::Field<scalar>& nodeCoords =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "coordinates");

    const stk::mesh::NgpMesh& ngpMesh =
        stk::mesh::get_updated_ngp_mesh(bulkData);
    stk::mesh::NgpField<double>& ngpCoordinatesRef =
        stk::mesh::get_updated_ngp_field<double>(nodeCoords);
    stk::mesh::NgpField<double>& ngpAreaVEdgeRef =
        stk::mesh::get_updated_ngp_field<double>(areaVEdgeField);

    const stk::mesh::Selector selAllElements = metaData.universal_part();

    stk::NgpVector<unsigned> bucketIds =
        ngpMesh.get_bucket_ids(stk::topology::ELEMENT_RANK, selAllElements);
    unsigned numBuckets = bucketIds.size();

    std::set<stk::topology> all_topos;
    for (unsigned i = 0; i < numBuckets; ++i)
    {
        const stk::mesh::NgpMesh::BucketType& myelementBucket =
            ngpMesh.get_bucket(stk::topology::ELEM_RANK, i);
        all_topos.insert(myelementBucket.topology());
    }

    accel::MasterElement* meSCS_ptrs[stk::topology::NUM_TOPOLOGIES] = {nullptr};
    for (const stk::topology& topo : all_topos)
    {
        label topo_index = topo.value();
        meSCS_ptrs[topo_index] =
            accel::MasterElementRepo::get_surface_master_element_on_dev(topo);
    }

    auto teamPolicy =
        stk::ngp::TeamPolicy<MeshToolsExecSpace>(numBuckets, Kokkos::AUTO);

    static constexpr int nodesPerElementMax = 8;
    static constexpr int numScsIpMax = 12;

    // entity array sizes
    const label connectedNodesSize = nodesPerElementMax;

    const label totalEntitySize = connectedNodesSize;

    const label totalIndexSize = 2 * connectedNodesSize;

    const size_t entityScratchSize =
        ScratchViewEntity::shmem_size(totalEntitySize);

    // DoubleType array sizes (for master element data)
    const size_t coordinatesScratchSize =
        ScratchView2D::shmem_size(nodesPerElementMax, SPATIAL_DIM);
    const size_t scsAreaVScratchSize =
        ScratchView2D::shmem_size(numScsIpMax, SPATIAL_DIM);

    // scratch size required for each thread
    const size_t totalThreadScratchSize =
        (entityScratchSize + coordinatesScratchSize + scsAreaVScratchSize);

    // scratch size required for each team
    const size_t totalTeamScratchSize = 0;

    const int SCRATCH_SPACE_LEVEL = 1;

    Kokkos::parallel_for(
        "areaVloop",
        teamPolicy.set_scratch_size(SCRATCH_SPACE_LEVEL,
                                    Kokkos::PerTeam(totalTeamScratchSize),
                                    Kokkos::PerThread(totalThreadScratchSize)),
        KOKKOS_LAMBDA( // TODO: class lambda here for numNodes_
            const TeamHandleType& teamMember) {
            const int bucketIndex =
                bucketIds.get<MeshToolsExecSpace>(teamMember.league_rank());
            const stk::mesh::NgpMesh::BucketType& elementBucket =
                ngpMesh.get_bucket(stk::topology::ELEM_RANK, bucketIndex);
            unsigned nElementsPerBucket = elementBucket.size();

            ScratchViewEntity connectedNodes(
                teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                connectedNodesSize);

            // geometry related to populate
            ScratchView2D ws_coordinates(
                teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                nodesPerElementMax,
                SPATIAL_DIM);
            ScratchView2D ws_scs_areav(
                teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                numScsIpMax,
                SPATIAL_DIM);

            label topology_index = elementBucket.topology();
            accel::MasterElement* meSCS = meSCS_ptrs[topology_index];

            // extract master element specifics
            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->num_integration_points();

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(teamMember, 0u, nElementsPerBucket),
                [&](const int& iElement)
            {
                // get elem
                stk::mesh::Entity elem = elementBucket[iElement];
                stk::mesh::FastMeshIndex elemFastIndex =
                    ngpMesh.fast_mesh_index(elem);

                stk::mesh::NgpMesh::ConnectedNodes nodeRels =
                    ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemFastIndex);
                label numNodes = nodeRels.size();

                // sanity check on num nodes
                STK_NGP_ThrowAssert(numNodes == nodesPerElement);

                for (label ni = 0; ni < numNodes; ++ni)
                {
                    stk::mesh::Entity node = nodeRels[ni];

                    // set connected nodes
                    connectedNodes[ni] = node;

                    stk::mesh::FastMeshIndex nodeFastIndex =
                        ngpMesh.fast_mesh_index(node);

                    for (unsigned i = 0; i < SPATIAL_DIM; ++i)
                    {
                        ws_coordinates(ni, i) =
                            ngpCoordinatesRef(nodeFastIndex, i);
                    }
                }

                meSCS->determinant(ws_coordinates, ws_scs_areav);

                const int* lrscv = meSCS->adjacentNodes();
                const int* scsIpEdgeMap = meSCS->scsIpEdgeOrd();

                const auto edges = ngpMesh.get_edges(
                    stk::topology::ELEM_RANK, ngpMesh.fast_mesh_index(elem));
                for (int ip = 0; ip < numScsIp; ++ip)
                {
                    // Edge for this integration point
                    const int nedge = scsIpEdgeMap[ip];
                    // Index of "left" node in the element relations
                    const int iLn = lrscv[2 * ip];

                    // Nodes connected to this edge
                    const auto edgeID = ngpMesh.fast_mesh_index(edges[nedge]);
                    const auto edge_nodes =
                        ngpMesh.get_nodes(stk::topology::EDGE_RANK, edgeID);

                    // Left node comparison
                    const auto lnElemId = connectedNodes[iLn];
                    const auto lnEdgeId = edge_nodes[0];

                    const double sign = (lnElemId == lnEdgeId) ? 1.0 : -1.0;

                    for (unsigned d = 0; d < SPATIAL_DIM; ++d)
                    {
                        Kokkos::atomic_add(&ngpAreaVEdgeRef(edgeID, d),
                                           ws_scs_areav(ip, d)._data.get() *
                                               sign);
                    }
                }
            });
        });
    areaVEdgeField.modify_on_device();
    areaVEdgeField.sync_to_host();
}

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          bool includeAdv = true,
          bool isShifted = true>
void preCalcAreaV(const stk::mesh::BulkData& bulkData)
{
    using ScratchView2D =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchView2D;
    using ScratchView3D =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::ScratchView3D;
    using TeamHandleType =
        typename AsmBase<scalar, label, BLOCKSIZE, SPATIAL_DIM>::TeamHandleType;

    // init ngp data
    const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
    const stk::mesh::Field<scalar>& areaVField =
        *metaData.get_field<scalar>(stk::topology::ELEMENT_RANK, "areaV");
    const stk::mesh::Field<scalar>& nodeCoords =
        *metaData.get_field<scalar>(stk::topology::NODE_RANK, "coordinates");

    const stk::mesh::NgpMesh& ngpMesh =
        stk::mesh::get_updated_ngp_mesh(bulkData);
    stk::mesh::NgpField<scalar>& ngpAreaVSTKFieldRef =
        stk::mesh::get_updated_ngp_field<scalar>(areaVField);

    stk::mesh::NgpField<double>& ngpCoordinatesRef =
        stk::mesh::get_updated_ngp_field<double>(nodeCoords);

    const stk::mesh::Selector selAllElements = metaData.universal_part();

    stk::NgpVector<unsigned> bucketIds =
        ngpMesh.get_bucket_ids(stk::topology::ELEMENT_RANK, selAllElements);
    unsigned numBuckets = bucketIds.size();

    // initialize master element pointer array
    std::set<stk::topology> all_topos;
    for (unsigned i = 0; i < numBuckets; ++i)
    {
        const stk::mesh::NgpMesh::BucketType& myelementBucket =
            ngpMesh.get_bucket(stk::topology::ELEM_RANK, i);
        all_topos.insert(myelementBucket.topology());
    }

    accel::MasterElement* meSCS_ptrs[stk::topology::NUM_TOPOLOGIES] = {nullptr};
    for (const stk::topology& topo : all_topos)
    {
        label topo_index = topo.value();
        meSCS_ptrs[topo_index] =
            accel::MasterElementRepo::get_surface_master_element_on_dev(topo);
    }

    auto teamPolicy =
        stk::ngp::TeamPolicy<MeshToolsExecSpace>(numBuckets, Kokkos::AUTO);

    static constexpr int nodesPerElementMax =
        8; // accel::AlgtraitsHex8::nodesPerElement_; //8;
    static constexpr int numScsIpMax =
        12; // accel::AlgtraitsHex8::numScsIp_; //12;

    // DoubleType array sizes (for master element data)
    const size_t coordinatesScratchSize =
        ScratchView2D::shmem_size(nodesPerElementMax, SPATIAL_DIM);
    const size_t scsAreaVScratchSize =
        ScratchView2D::shmem_size(numScsIpMax, SPATIAL_DIM);

    // scratch size required for each thread
    const size_t totalThreadScratchSize =
        (coordinatesScratchSize + scsAreaVScratchSize);

    // scratch size required for each team
    const size_t totalTeamScratchSize = 0;

    const int SCRATCH_SPACE_LEVEL = 1;

    Kokkos::parallel_for(
        "team_assembly_loop",
        teamPolicy.set_scratch_size(SCRATCH_SPACE_LEVEL,
                                    Kokkos::PerTeam(totalTeamScratchSize),
                                    Kokkos::PerThread(totalThreadScratchSize)),
        KOKKOS_LAMBDA(const TeamHandleType& teamMember) {
            const int bucketIndex =
                bucketIds.get<MeshToolsExecSpace>(teamMember.league_rank());
            const stk::mesh::NgpMesh::BucketType& elementBucket =
                ngpMesh.get_bucket(stk::topology::ELEM_RANK, bucketIndex);
            unsigned nElementsPerBucket = elementBucket.size();

            // geometry related to populate
            ScratchView2D ws_coordinates(
                teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                nodesPerElementMax,
                SPATIAL_DIM);
            ScratchView2D ws_scs_areav(
                teamMember.thread_scratch(SCRATCH_SPACE_LEVEL),
                numScsIpMax,
                SPATIAL_DIM);

            label topology_index = elementBucket.topology();
            accel::MasterElement* meSCS = meSCS_ptrs[topology_index];

            // extract master element specifics
            const label nodesPerElement = meSCS->nodesPerElement_;
            const label numScsIp = meSCS->num_integration_points();
            const label* lrscv = meSCS->adjacentNodes();

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange(teamMember, 0u, nElementsPerBucket),
                [&](const int& iElement)
            {
                // get elem
                stk::mesh::Entity elem = elementBucket[iElement];
                stk::mesh::FastMeshIndex elemFastIndex =
                    ngpMesh.fast_mesh_index(elem);

                stk::mesh::NgpMesh::ConnectedNodes nodeRels =
                    ngpMesh.get_nodes(stk::topology::ELEM_RANK, elemFastIndex);
                label numNodes = nodeRels.size(); // bulkData.num_nodes(elem);

                // sanity check on num nodes
                STK_NGP_ThrowAssert(numNodes == nodesPerElement);

                for (label ni = 0; ni < numNodes; ++ni)
                {
                    stk::mesh::Entity node = nodeRels[ni];

                    stk::mesh::FastMeshIndex nodeFastIndex =
                        ngpMesh.fast_mesh_index(node);

                    // gather vectors
                    for (unsigned i = 0; i < SPATIAL_DIM; ++i)
                    {
                        ws_coordinates(ni, i) =
                            ngpCoordinatesRef(nodeFastIndex, i);
                    }
                }

                meSCS->determinant(ws_coordinates, ws_scs_areav);
                for (label ip = 0; ip < numScsIp; ++ip)
                {
                    for (unsigned j = 0; j < SPATIAL_DIM; ++j)
                    {
                        const scalar axj = ws_scs_areav(ip, j)._data.get();
                        const label idx = ip * SPATIAL_DIM + j;
                        ngpAreaVSTKFieldRef(elemFastIndex, idx) = axj;
                    }
                }
            });
        });

    ngpAreaVSTKFieldRef.modify_on_device();
    ngpAreaVSTKFieldRef.sync_to_host();
}

#endif // MESHTESTTOOLS_H
