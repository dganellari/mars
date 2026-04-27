#pragma once

// MARS Adaptive Mesh Refinement for the unstructured backend
//
// GPU-native AMR built on cornerstone-octree.
// Works with any element-based PDE solver (FEM, CVFEM, etc.)
//
// Components:
//   AmrOctree        — octree-based refinement level tracking + 1-irregular enforcement
//   HexRefiner       — GPU-native hex8 isotropic refinement (1 hex → 8 children)
//   ErrorIndicator   — gradient-based error indicator (one example; users can provide their own)
//   SolutionTransfer — interpolation of solution from old to new mesh
//   AmrManager       — orchestrates the full AMR cycle
//
// Usage:
//   AmrManager<HexTag, KeyType, RealType> amr(config);  // or TetTag
//   amr.initialize(meshFile, rank, numRanks);
//
//   while (amr.shouldContinue(stats)) {
//       // Solve on current mesh
//       auto& domain = amr.domain();
//       solve(domain, solution);
//
//       // Compute error indicator (user-defined)
//       auto d_error = computeMyError(domain, solution);
//
//       // Adapt mesh
//       DeviceVector<RealType> newSolution;
//       auto stats = amr.adaptMesh(d_error.data(), solution.data(), newSolution);
//       solution = std::move(newSolution);
//   }

#include "mars_amr_octree.hpp"
#include "mars_amr_hex_refine.hpp"
#include "mars_amr_tet_refine.hpp"
#include "mars_amr_error_indicator.hpp"
#include "mars_amr_solution_transfer.hpp"
#include "mars_amr_octree_native.hpp"
#include "mars_amr_manager.hpp"
