#ifndef ASM_BASE_H
#define ASM_BASE_H

#include "linearSolverContext.h"

// #include "Kokkos_Core.hpp"

// #include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetNgpField.hpp>
#include <stk_mesh/base/GetNgpMesh.hpp>

#include <MasterElement.h>
#include <MasterElementRepo.h>

#include <map>
#include <string>

template <typename scalar,
          typename label,
          size_t BLOCKSIZE,
          size_t SPATIAL_DIM,
          bool includeAdv = true,
          bool isShifted = true>
class AsmBase
{
public:
    using Context = ::linearSolver::Context<BLOCKSIZE>;
    using Matrix = typename Context::Matrix;
    using Vector = typename Context::Vector;

    using DataType = scalar;
    using ExecSpace = stk::mesh::NgpMesh::MeshExecSpace;
    using ScratchViewScalar =
        Kokkos::View<scalar*, ExecSpace::scratch_memory_space>;
    using ScratchViewScalar2 =
        Kokkos::View<scalar**, ExecSpace::scratch_memory_space>;
    using ScratchViewScalar3 =
        Kokkos::View<scalar***, ExecSpace::scratch_memory_space>;
    using ScratchViewLabel =
        Kokkos::View<label*, ExecSpace::scratch_memory_space>;
    using ScratchViewEntity =
        Kokkos::View<stk::mesh::Entity*, ExecSpace::scratch_memory_space>;
    using ScratchViewIndex =
        Kokkos::View<typename Context::Index*, ExecSpace::scratch_memory_space>;
    using ScratchView2D =
        accel::SharedMemView<accel::DoubleType**,
                             ExecSpace::scratch_memory_space>; // for master
                                                               // element data
    using ScratchView3D =
        accel::SharedMemView<accel::DoubleType***,
                             ExecSpace::scratch_memory_space>; // for master
                                                               // element data
    using TeamHandleType = stk::ngp::TeamPolicy<ExecSpace>::member_type;

    AsmBase(const stk::mesh::BulkData& bulk)
    {
        numNodes_ =
            stk::mesh::count_entities(bulk,
                                      stk::topology::NODE_RANK,
                                      bulk.mesh_meta_data().universal_part());
    }

    virtual ~AsmBase() = default;

    virtual std::string name() const = 0;

    virtual void assemble(Matrix& A,
                          Vector& b,
                          const stk::mesh::MetaData& metaData,
                          const stk::mesh::BulkData& bulkData) const = 0;

    virtual void getScratchSizes(size_t& perTeam, size_t& perThread) const = 0;

protected:
    unsigned numNodes_;
};

#endif // ASM_BASE_H
