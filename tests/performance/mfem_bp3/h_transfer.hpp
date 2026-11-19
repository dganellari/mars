#pragma once
#include "mfem.hpp"
#include <memory>

namespace bp3_transfer
{
inline void add_assembled_level(mfem::ParFiniteElementSpaceHierarchy& hierarchy)
{
    auto& coarse = hierarchy.GetFinestFESpace();
    auto mesh = std::make_unique<mfem::ParMesh>(*coarse.GetParMesh());
    mesh->UniformRefinement();
    auto fine = std::make_unique<mfem::ParFiniteElementSpace>(mesh.get(), coarse.FEColl(),
                                                            coarse.GetVDim(), coarse.GetOrdering());
    // ANY_TYPE selects a host element loop in MFEM 4.8; ParCSR applies the same map on the device.
    mfem::OperatorHandle transfer(mfem::Operator::Hypre_ParCSR);
    fine->GetTrueTransferOperator(coarse, transfer);
    hierarchy.AddLevel(mesh.get(), fine.get(), transfer.Ptr(), true, true, true);
    mesh.release();
    fine.release();
    transfer.SetOperatorOwner(false);
}
} // namespace bp3_transfer
