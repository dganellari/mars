#pragma once

// Main FEM header - includes all components

#include "mars_reference_element.hpp"
#include "mars_sparse_matrix.hpp"
#include "mars_h1_fe_space.hpp"
#include "mars_stiffness_assembler.hpp"
#include "mars_mass_assembler.hpp"
#include "mars_boundary_conditions.hpp"
#include "mars_dof_elimination.hpp"
#include "mars_unstructured_dof_handler.hpp"
#include "../solvers/mars_cg_solver.hpp"
#include "../solvers/mars_bicgstab_solver.hpp"
#include "../solvers/mars_gmres_solver.hpp"
#include "mars_debug_utils.hpp"
// #include "mars_vtk_writer.hpp"

namespace mars {
namespace fem {

// Convenience aliases for common configurations
template<typename RealType = float, typename KeyType = uint64_t>
using TetFESpace = H1FESpace<TetTag, RealType, KeyType, cstone::GpuTag>;

template<typename RealType = float, typename KeyType = uint64_t>
using TetStiffnessAssembler = StiffnessAssembler<TetTag, RealType, KeyType, cstone::GpuTag>;

template<typename RealType = float, typename KeyType = uint64_t>
using TetMassAssembler = MassAssembler<TetTag, RealType, KeyType, cstone::GpuTag>;

template<typename RealType = float, typename KeyType = uint64_t>
using TetBCHandler = BoundaryConditionHandler<TetTag, RealType, KeyType, cstone::GpuTag>;

template<typename RealType = float, typename KeyType = uint64_t>
using TetCGSolver = ConjugateGradientSolver<RealType, KeyType, cstone::GpuTag>;

template<typename RealType = float, typename KeyType = uint64_t>
using TetBiCGSTABSolver = BiCGSTABSolver<RealType, KeyType, cstone::GpuTag>;

template<typename RealType = float, typename KeyType = uint64_t>
using TetGMRESSolver = GMRESSolver<RealType, KeyType, cstone::GpuTag>;

template<typename RealType = float, typename KeyType = uint64_t>
using TetSparseMatrix = SparseMatrix<KeyType, RealType, cstone::GpuTag>;

template<typename RealType = float, typename KeyType = uint64_t>
using TetUnstructuredDofHandler = UnstructuredDofHandler<TetTag, RealType, KeyType, cstone::GpuTag>;

} // namespace fem
} // namespace mars
