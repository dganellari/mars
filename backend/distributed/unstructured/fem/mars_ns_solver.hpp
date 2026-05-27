#pragma once
//
// Incompressible Navier-Stokes Chorin projection solver, header-only.
//
// Extracted from mars_amr_ns_projection.cu so multiple drivers (cavity,
// channel, TGV with periodic BCs) can share one solver. This is a verbatim
// lossless extraction: the cavity/channel driver behavior is preserved
// bit-for-bit. The Periodic BCKind, runNsStep split into sub-steps, and
// removeMean integration land in follow-up commits.
//
// Per timestep (current monolithic runNsStep):
//   1) PREDICTOR        (explicit advection + grad p, halo exchange)
//   2) IMPLICIT DIFFUSION (3 CG/Hypre solves for u,v,w)
//   3) PRESSURE POISSON   (K-path or DDT-path; pin removes constant mode)
//   4) CORRECTOR        (project velocity, update pressure)
//   5) DIVERGENCE DIAG  (per-step max |div u| for monitoring)
//
// BCKind: Cavity | Channel (the Periodic mode lands in a follow-up).
//
// Note on namespacing: kernels and helpers live at file scope (no enclosing
// namespace) to preserve the original symbol names. Drivers that include
// this header get the same global symbols they had pre-refactor.
//

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include "backend/distributed/unstructured/fem/mars_perf_counters.hpp"
#include "backend/distributed/unstructured/fem/mars_sparse_matrix.hpp"
#include "backend/distributed/unstructured/fem/mars_periodic_bc.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#ifdef MARS_ENABLE_HYPRE
#include "backend/distributed/unstructured/solvers/mars_hypre_pcg_solver.hpp"
#include "backend/distributed/unstructured/solvers/mars_hypre_gmres_solver.hpp"
#endif
#include "backend/distributed/unstructured/amr/mars_amr.hpp"

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>
#include <thrust/extrema.h>
#include <thrust/fill.h>
#include <thrust/system/cuda/execution_policy.h>

#include <memory>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

// Sync on lap so wall-clock reflects all-rank work, not rank-0's stragglers.
// Also live-prints each phase on rank 0 so long setup phases on large meshes
// (e.g. 15M-hex wing) show progress in real time instead of a silent block
// before the final summary.
struct PhaseTimer
{
    using clk = std::chrono::high_resolution_clock;
    clk::time_point t0;
    std::vector<std::pair<std::string, float>> phases;
    bool sync;
    int  liveRank = -1;       // rank that should print laps live (-1 = no live print)

    explicit PhaseTimer(bool sync_ = true) : sync(sync_) { reset(); }
    explicit PhaseTimer(bool sync_, int liveRank_) : sync(sync_), liveRank(liveRank_) { reset(); }

    void reset()
    {
        if (sync) { cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD); }
        t0 = clk::now();
    }

    void lap(const std::string& name)
    {
        if (sync) { cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD); }
        auto t1 = clk::now();
        float ms = std::chrono::duration<float, std::milli>(t1 - t0).count();
        phases.emplace_back(name, ms);
        t0 = t1;
        int worldRank = 0;
        if (liveRank >= 0)
        {
            MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
            if (worldRank == liveRank)
            {
                std::cout << "  [lap] " << std::left << std::setw(36) << name
                          << std::right << std::fixed << std::setprecision(2)
                          << std::setw(12) << ms << " ms\n" << std::flush;
            }
        }
    }

    void report(int rank, const std::string& header)
    {
        if (rank != 0) return;
        std::cout << "  [" << header << "] phase breakdown:\n";
        float total = 0;
        for (auto& [n, ms] : phases) total += ms;
        for (auto& [n, ms] : phases)
        {
            std::cout << "    " << std::left << std::setw(36) << n
                      << std::right << std::fixed << std::setprecision(2)
                      << std::setw(10) << ms << " ms ("
                      << std::setw(5) << std::setprecision(1)
                      << (100.0f * ms / std::max(total, 1e-3f)) << "%)\n";
        }
        std::cout << "    " << std::left << std::setw(36) << "TOTAL"
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(10) << total << " ms\n";
    }
};

// =============================================================================
// BC infrastructure: per-component velocity Dirichlet flags + per-node target
// values for the lid-driven cavity; pressure pin DOF id for the Poisson solve.
// =============================================================================

// Cavity BC: u=1 on top face (z=z_max), q=0 on the other 5 faces; v=w=0 on all
// 6 faces. Sets a per-node boundary mask AND per-node target values for each
// velocity component. The boundary mask is computed in DOF order (only owned
// DOFs need it; ghosts inherit through halo exchange of velocity).
//
// We compute per-NODE target values so the predictor/corrector can simply
// overwrite owned-boundary node slots; the implicit diffusion solve then
// inherits the same values via the lift-off RHS.
template<typename RealType>
__global__ void markCavityBCKernel(const RealType* nodeX,
                                   const RealType* nodeY,
                                   const RealType* nodeZ,
                                   const uint8_t* ownership,
                                   const int* nodeToDof,
                                   uint8_t* isBdryDof,        // size numOwnedDofs
                                   RealType* uTarget,         // size nodeCount
                                   RealType* vTarget,
                                   RealType* wTarget,
                                   size_t numNodes,
                                   RealType xmin, RealType xmax,
                                   RealType ymin, RealType ymax,
                                   RealType zmin, RealType zmax,
                                   RealType eps,
                                   RealType lidU)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;

    RealType x = nodeX[i], y = nodeY[i], z = nodeZ[i];
    bool onBdry = (fabs(x - xmin) < eps || fabs(x - xmax) < eps ||
                   fabs(y - ymin) < eps || fabs(y - ymax) < eps ||
                   fabs(z - zmin) < eps || fabs(z - zmax) < eps);
    bool onTop  = (fabs(z - zmax) < eps);
    // Pure-Dirichlet for velocity: only set targets on the actual boundary.
    // Interior nodes' targets are unused (overwritten by the solve).
    if (onBdry)
    {
        uTarget[i] = onTop ? lidU : RealType(0);
        vTarget[i] = RealType(0);
        wTarget[i] = RealType(0);
    }
    else
    {
        // Initialize interior to 0 so the predictor/corrector overwrite logic
        // doesn't read stale memory if the array happened to alias.
        uTarget[i] = RealType(0);
        vTarget[i] = RealType(0);
        wTarget[i] = RealType(0);
    }

    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    isBdryDof[dof] = onBdry ? 1 : 0;
}

// Channel flow BC marker.
//   Velocity:
//     Inflow face (x = xmin):  Dirichlet u = Uinf, v = w = 0
//     Outflow face (x = xmax): natural homogeneous-Neumann (do nothing — not
//                              marked in isBdryDof)
//     Tunnel walls (y=ymin/ymax, z=zmin/zmax): Dirichlet no-slip u = v = w = 0
//   Pressure: handled separately by markChannelPressureBCKernel (Dirichlet
//   p = 0 on the entire outflow face). The corner-pin used by the cavity
//   setup is skipped in channel mode — pinning a single corner fights the
//   natural inflow-to-outflow pressure gradient and a free null space lets
//   y-z symmetric modes grow.
template<typename RealType>
__global__ void markChannelBCKernel(const RealType* nodeX,
                                    const RealType* nodeY,
                                    const RealType* nodeZ,
                                    const uint8_t* ownership,
                                    const int* nodeToDof,
                                    uint8_t* isBdryDof,
                                    RealType* uTarget,
                                    RealType* vTarget,
                                    RealType* wTarget,
                                    size_t numNodes,
                                    RealType xmin, RealType xmax,
                                    RealType ymin, RealType ymax,
                                    RealType zmin, RealType zmax,
                                    RealType eps,
                                    RealType Uinf)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;

    RealType x = nodeX[i], y = nodeY[i], z = nodeZ[i];
    bool onInflow = (fabs(x - xmin) < eps);
    bool onWall   = (fabs(y - ymin) < eps || fabs(y - ymax) < eps ||
                     fabs(z - zmin) < eps || fabs(z - zmax) < eps);
    // Outflow (x=xmax) deliberately excluded from onBdry so the natural BC
    // (zero-grad / weak Neumann) applies through the matrix.
    bool onBdry = onInflow || onWall;

    if (onBdry)
    {
        uTarget[i] = onInflow ? Uinf : RealType(0);
        vTarget[i] = RealType(0);
        wTarget[i] = RealType(0);
    }
    else
    {
        uTarget[i] = RealType(0);
        vTarget[i] = RealType(0);
        wTarget[i] = RealType(0);
    }

    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    isBdryDof[dof] = onBdry ? 1 : 0;
}

// Channel pressure BC marker. Flags every owned node on the outflow face
// (x = xmax) as a pressure-Dirichlet DOF. The pressure-matrix rows for those
// DOFs become identity (row=0, diag=1) and the pressure RHS rows are zeroed
// (target value = 0) so phi[outflow] = 0 exactly. This replaces the
// single-point pin used by cavity and removes the null space without
// fighting the natural channel pressure gradient.
template<typename RealType>
__global__ void markChannelPressureBCKernel(const RealType* nodeX,
                                            const uint8_t* ownership,
                                            const int* nodeToDof,
                                            uint8_t* isPressureBdryDof,
                                            size_t numNodes,
                                            RealType xmax,
                                            RealType eps)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    bool onOutflow = fabs(nodeX[i] - xmax) < eps;
    isPressureBdryDof[dof] = onOutflow ? 1 : 0;
}

// Mark the single pressure-pin DOF. A pure-Dirichlet velocity problem leaves
// the pressure Poisson rank-deficient by one (constant mode). Pinning ONE node
// removes the singularity at minimal cost; we pick the lowest-rank node closest
// to (xmin, ymin, zmin) so the choice is deterministic across runs.
//
// Approach: each rank finds its own closest owned DOF to the corner, ties
// broken by lowest local DOF id. A single MPI_MINLOC then chooses the global
// winner (which rank holds the pin), and that rank flips the bit. Returns the
// winning rank in pinRankOut (only meaningful on rank 0).
template<typename RealType>
__global__ void findPressurePinCandidateKernel(const RealType* nodeX,
                                               const RealType* nodeY,
                                               const RealType* nodeZ,
                                               const uint8_t* ownership,
                                               const int* nodeToDof,
                                               RealType* d2PerDof,
                                               int* dofIdPerDof,
                                               size_t numNodes,
                                               int numOwnedDofs,
                                               RealType xPin, RealType yPin, RealType zPin)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;
    RealType dx = nodeX[i] - xPin;
    RealType dy = nodeY[i] - yPin;
    RealType dz = nodeZ[i] - zPin;
    d2PerDof[dof]    = dx * dx + dy * dy + dz * dz;
    dofIdPerDof[dof] = dof;
}

// Pressure pin: row=0, diag=1, RHS=0 on the chosen DOF (and only on the rank
// that owns it). All other rows have natural Neumann (no row touched).
template<typename RealType>
__global__ void enforcePinRowMatrixKernel(int pinDof,
                                          const int* rowPtr,
                                          const int* colInd,
                                          const int* diagPtr,
                                          RealType* values)
{
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t != 0) return;
    if (pinDof < 0) return;
    for (int j = rowPtr[pinDof]; j < rowPtr[pinDof + 1]; ++j)
        values[j] = RealType(0);
    int dp = diagPtr[pinDof];
    if (dp >= 0) values[dp] = RealType(1);
}

// Zero RHS at the pin (called every step before the pressure solve).
template<typename RealType>
__global__ void enforcePinRhsKernel(int pinDof, RealType* rhs)
{
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t != 0) return;
    if (pinDof >= 0) rhs[pinDof] = RealType(0);
}

// Symmetric pin: zero column pinDof in every owned row i != pinDof. The row
// kernel (enforcePinRowMatrixKernel) zeroes the pin's row + sets its diagonal
// to 1, but the original off-diagonal couplings A[i, pinDof] in OTHER rows
// remain, leaving the matrix non-symmetric. Hypre BoomerAMG's coarsening
// heuristic assumes near-symmetry and degrades on a one-sided pin (suspected
// contributor to the cube16 DDT+Hypre divergence). Pin value = 0, so no RHS
// lift is needed: we simply zero the column entries.
//
// One thread per owned row. Each row holds at most one column equal to pinDof
// (CSR rows have unique columns), so a linear scan + direct write (no atomic)
// is correct.
template<typename RealType>
__global__ void enforcePinColMatrixKernel(int pinDof,
                                          const int* rowPtr,
                                          const int* colInd,
                                          RealType* values,
                                          int numOwnedRows)
{
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= numOwnedRows) return;
    if (pinDof < 0) return;
    if (row == pinDof) return;  // row kernel already handled this row
    int rs = rowPtr[row];
    int re = rowPtr[row + 1];
    for (int j = rs; j < re; ++j)
    {
        if (colInd[j] == pinDof)
        {
            values[j] = RealType(0);
            return;
        }
    }
}

// Multi-DOF symmetric Dirichlet (channel outflow, wing outlet). For every
// owned non-Dirichlet row, zero every column whose dof is flagged in the
// boundary mask. Pin value = 0, no RHS lift.
template<typename RealType>
__global__ void enforceBcColMatrixKernel(const uint8_t* isBoundaryDof,
                                         const int* rowPtr,
                                         const int* colInd,
                                         RealType* values,
                                         int numOwnedRows)
{
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= numOwnedRows) return;
    if (isBoundaryDof[row]) return;  // already an identity row
    int rs = rowPtr[row];
    int re = rowPtr[row + 1];
    for (int j = rs; j < re; ++j)
    {
        int c = colInd[j];
        if (c >= 0 && c < numOwnedRows && isBoundaryDof[c])
            values[j] = RealType(0);
    }
}

// Zero RHS at every Dirichlet-pressure DOF (channel outflow face). Each
// thread handles one owned DOF.
template<typename RealType>
__global__ void enforcePressureBcRhsKernel(const uint8_t* isPressureBdryDof,
                                           RealType* rhs,
                                           int numDofs)
{
    int dof = blockIdx.x * blockDim.x + threadIdx.x;
    if (dof >= numDofs) return;
    if (isPressureBdryDof[dof]) rhs[dof] = RealType(0);
}

// Build a per-NODE pressure-Dirichlet mask from the existing per-DOF mask +
// optional single-pin DOF. Each owned node maps its dof through the existing
// boundary flag; ghost nodes get 0 (filled by halo exchange later).
template<typename RealType>
__global__ void buildPressureBdryNodeMaskKernel(const int* nodeToDof,
                                                const uint8_t* ownership,
                                                const uint8_t* isPressureBdryDof, // may be null
                                                int pinDof,                       // -1 if no pin
                                                uint8_t* isPressureBdryNode,
                                                size_t numNodes,
                                                int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    isPressureBdryNode[i] = 0;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;
    bool isPin  = (pinDof >= 0 && dof == pinDof);
    bool isMask = (isPressureBdryDof != nullptr && isPressureBdryDof[dof] != 0);
    isPressureBdryNode[i] = (isPin || isMask) ? uint8_t(1) : uint8_t(0);
}

// Symmetric column-zero for the assembled DDT CSR, using the per-NODE
// pressure-Dirichlet mask (halo-exchanged so ghosts see the owner's flag).
// For every owned row whose own node is NOT pressure-Dirichlet, scan its CSR
// neighbors; for any neighbor (local DOF c) whose corresponding LOCAL NODE is
// flagged Dirichlet, zero values[j]. This eliminates pin/mask couplings on
// EVERY rank (including ghost-copies of the pin on non-owning ranks), which
// is necessary for BoomerAMG to see a globally symmetric operator.
//
// Maps a column local DOF c back to a local node by walking the nodeToDof
// inverse via a precomputed dofToNode array (passed in as dofToNode[c]).
template<typename RealType>
__global__ void enforceBcColMatrixKernelByNode(const uint8_t* isPressureBdryNode,
                                               const int* nodeToDof,
                                               const int* dofToNode,
                                               int numTotalDofs,
                                               const int* rowPtr,
                                               const int* colInd,
                                               RealType* values,
                                               int numOwnedRows)
{
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= numOwnedRows) return;
    // Skip if THIS row is itself Dirichlet (already an identity row).
    int rowNode = dofToNode[row];
    if (rowNode >= 0 && isPressureBdryNode[rowNode]) return;
    int rs = rowPtr[row];
    int re = rowPtr[row + 1];
    for (int j = rs; j < re; ++j)
    {
        int c = colInd[j];
        if (c < 0 || c >= numTotalDofs) continue;
        int cnode = dofToNode[c];
        if (cnode >= 0 && isPressureBdryNode[cnode])
            values[j] = RealType(0);
    }
}

// Build the inverse map: dofToNode[d] = local node that has nodeToDof[node] = d.
// One thread per node. Multiple nodes can never share a DOF (in non-periodic
// modes), so this is a simple scatter.
__global__ void buildDofToNodeKernel(const int* nodeToDof,
                                     int* dofToNode,
                                     int numTotalDofs,
                                     size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    int dof = nodeToDof[i];
    if (dof >= 0 && dof < numTotalDofs) dofToNode[dof] = static_cast<int>(i);
}

// Velocity Dirichlet: zero the rows of boundary DOFs and put 1 on the diagonal.
// Applied ONCE to the (M/dt + nu K) matrix so all three components inherit the
// same row structure. Per-component target values are pushed through the RHS.
template<typename RealType>
__global__ void enforceBcMatrixKernel(const uint8_t* isBoundaryDof,
                                      const int* rowPtr,
                                      const int* colInd,
                                      const int* diagPtr,
                                      RealType* values,
                                      int numDofs)
{
    int dof = blockIdx.x * blockDim.x + threadIdx.x;
    if (dof >= numDofs || !isBoundaryDof[dof]) return;
    for (int j = rowPtr[dof]; j < rowPtr[dof + 1]; ++j)
        values[j] = RealType(0);
    int dp = diagPtr[dof];
    if (dp >= 0) values[dp] = RealType(1);
}

// Lift-off velocity Dirichlet on the RHS: rhs[dof] = qTarget[node(dof)] for
// owned boundary nodes. Combined with the matrix BC above, this implements
// q[boundary] = qTarget exactly.
template<typename RealType>
__global__ void enforceBcRhsFromTargetKernel(const uint8_t* isBoundaryDof,
                                             const RealType* qTarget,
                                             const int* nodeToDof,
                                             const uint8_t* ownership,
                                             RealType* rhs,
                                             size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    if (isBoundaryDof[dof]) rhs[dof] = qTarget[i];
}

// =============================================================================
// Lumped mass: per-node + reverse-halo pattern (same as B.2/B.3/B.4). Single
// element-volume scatter feeds both the diffusion matrix (added to its diag as
// M/dt) and the per-step predictor/corrector normalization. V_node matches the
// per-node scatter used by gradient/divergence so the projection step is
// algebraically consistent across operators.
// =============================================================================

template<typename KeyType, typename RealType>
__global__ void computeLumpedMassPerNodeKernel(const KeyType* c0, const KeyType* c1,
                                               const KeyType* c2, const KeyType* c3,
                                               const KeyType* c4, const KeyType* c5,
                                               const KeyType* c6, const KeyType* c7,
                                               const RealType* nodeX,
                                               const RealType* nodeY,
                                               const RealType* nodeZ,
                                               RealType* massNode,
                                               size_t startElem,
                                               size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;
    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};

    RealType x[8], y[8], z[8];
    for (int i = 0; i < 8; ++i)
    {
        x[i] = nodeX[n[i]];
        y[i] = nodeY[n[i]];
        z[i] = nodeZ[n[i]];
    }
    RealType xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0], zmin = z[0], zmax = z[0];
    for (int i = 1; i < 8; ++i)
    {
        if (x[i] < xmin) xmin = x[i]; if (x[i] > xmax) xmax = x[i];
        if (y[i] < ymin) ymin = y[i]; if (y[i] > ymax) ymax = y[i];
        if (z[i] < zmin) zmin = z[i]; if (z[i] > zmax) zmax = z[i];
    }
    RealType contrib = (xmax - xmin) * (ymax - ymin) * (zmax - zmin) * RealType(0.125);

    for (int i = 0; i < 8; ++i)
        atomicAdd(&massNode[n[i]], contrib);
}

// Owned-node mass slots -> per-DOF mass array (post reverse halo). atomicAdd
// because multiple nodes may map to the same DOF under periodic collapse
// (slave + master both index the master DOF); plain assignment would race.
template<typename RealType>
__global__ void gatherOwnedNodeMassToDofKernel(const RealType* massNode,
                                               const int* nodeToDof,
                                               const uint8_t* ownership,
                                               RealType* massDof,
                                               size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    atomicAdd(&massDof[dof], massNode[i]);
}

// Add M[i]/dt onto the diagonal of (nu K). Runs once at setup, before BC enforce.
template<typename RealType>
__global__ void addLumpedMassDiagonalKernel(const RealType* mass,
                                            const int* diagPtr,
                                            RealType invDt,
                                            RealType* values,
                                            int numOwnedDofs)
{
    int dof = blockIdx.x * blockDim.x + threadIdx.x;
    if (dof >= numOwnedDofs) return;
    int dp = diagPtr[dof];
    if (dp >= 0) values[dp] += mass[dof] * invDt;
}

// =============================================================================
// Advection scatter: per-node + reverse-halo. Computes net out-flux at each
// owned node for each velocity component q in {u, v, w}. Reused for all
// three predictor components by passing different q arrays.
//
// Three forms are available, selected at runtime via the `useSkewSymmetric` flag:
//
//   MARS_NS_ADV_UPWIND (default for cavity/channel; useSkewSymmetric=false):
//     q_face = q[iL] if mdot > 0 else q[iR]   (1st-order upwind)
//     flux   = mdot * q_face
//     Stable for any Re including high-Re, but injects/extracts KE proportional
//     to |mdot|. Fine when boundary dissipation (Dirichlet walls) absorbs it;
//     fatal for periodic NS because there is no wall to absorb the energy.
//
//   MARS_NS_ADV_SKEW (recommended for periodic; useSkewSymmetric=true):
//     Verstappen-Veldman discrete skew-symmetric form:
//       N_skew(u,q)_i = 1/2 (N_div + N_con)_i  =  N_div_i - 1/2 q_i (div u)_i
//     On an SCS-face flux scatter this collapses to the "3:1 stencil":
//       L gets  -0.25 * mdot * (3 q[L] + q[R])
//       R gets  +0.25 * mdot * (q[L] + 3 q[R])
//     Property: (N_skew q, q)_M = 0 exactly for ANY velocity field u (even
//     when div u != 0). Discrete kinetic energy is conserved by advection.
//     This is the canonical Nek5000 / Verstappen / Sandia-Nalu form.
//     Energy conservation is a SPATIAL property; works with any time
//     integrator (forward Euler still has O(dt) time-truncation drift).
//
// Sign convention matches the legacy upwind path: net flux out of L, into R.
// =============================================================================

template<typename KeyType, typename RealType>
__global__ void explicitAdvectionFluxScatterPerNodeKernel(const KeyType* c0, const KeyType* c1,
                                                          const KeyType* c2, const KeyType* c3,
                                                          const KeyType* c4, const KeyType* c5,
                                                          const KeyType* c6, const KeyType* c7,
                                                          const RealType* vx,
                                                          const RealType* vy,
                                                          const RealType* vz,
                                                          const RealType* q,
                                                          const RealType* areaVecX,
                                                          const RealType* areaVecY,
                                                          const RealType* areaVecZ,
                                                          RealType* dqdtNode,
                                                          size_t startElem,
                                                          size_t numLocal,
                                                          bool useSkewSymmetric)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;
    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};

    #pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        int nodeL = d_hexLRSCV[ip * 2];
        int nodeR = d_hexLRSCV[ip * 2 + 1];
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];

        RealType vfx = RealType(0.5) * (vx[iL] + vx[iR]);
        RealType vfy = RealType(0.5) * (vy[iL] + vy[iR]);
        RealType vfz = RealType(0.5) * (vz[iL] + vz[iR]);

        size_t off = e * 12 + ip;
        RealType mdot = vfx * areaVecX[off] + vfy * areaVecY[off] + vfz * areaVecZ[off];

        if (useSkewSymmetric)
        {
            // Verstappen-Veldman discrete skew-symmetric advection.
            //
            // Decomposition:
            //   N_skew(u,q) = 1/2 [ div(u q)  +  (u . grad) q ]
            //              = div(u q) - 1/2 q (div u)
            //
            // Per face f between (L, R), with mdot outflow-positive at L:
            //   divergence flux: L <- -mdot*0.5*(q_L + q_R) ; R <- +mdot*0.5*(q_L + q_R)
            //   per-node -1/2 q (div u): at L it adds -1/2 q_L * mdot
            //                            at R it adds +1/2 q_R * mdot
            // Sum per face:
            //   L:  -0.5 * mdot * (2 q_L + q_R)
            //   R:  +0.5 * mdot * (q_L + 2 q_R)
            //
            // Discrete property: sum_i q_i (dq/dt)_i = 0 for ANY u, even when
            // div u != 0. Discrete KE is conserved by advection. Forward-Euler
            // still has O(dt) time-truncation drift, but no exponential
            // injection. Spatial property -- works with any time integrator.
            RealType qL = q[iL];
            RealType qR = q[iR];
            atomicAdd(&dqdtNode[iL], -RealType(0.5) * mdot * (RealType(2) * qL + qR));
            atomicAdd(&dqdtNode[iR], +RealType(0.5) * mdot * (qL + RealType(2) * qR));
        }
        else
        {
            // 1st-order upwind (standard CVFEM, used by cavity/channel).
            RealType q_face = (mdot > RealType(0)) ? q[iL] : q[iR];
            RealType flux   = mdot * q_face;
            atomicAdd(&dqdtNode[iL], -flux);
            atomicAdd(&dqdtNode[iR], +flux);
        }
    }
}

// =============================================================================
// Gradient scatter (B.4 pattern): three components in a single pass per call.
// For each SCS face f with midpoint pressure p_f = 0.5*(p[iL]+p[iR]) and area
// vector A_f pointing L->R: +p_f*A_f to L, -p_f*A_f to R. Three separate
// buffers because CUDA has no native 3-vector atomicAdd.
// =============================================================================

template<typename KeyType, typename RealType>
__global__ void computeGradientPerNodeKernel(const KeyType* c0, const KeyType* c1,
                                             const KeyType* c2, const KeyType* c3,
                                             const KeyType* c4, const KeyType* c5,
                                             const KeyType* c6, const KeyType* c7,
                                             const RealType* p,
                                             const RealType* areaVecX,
                                             const RealType* areaVecY,
                                             const RealType* areaVecZ,
                                             RealType* gxAccNode,
                                             RealType* gyAccNode,
                                             RealType* gzAccNode,
                                             size_t startElem,
                                             size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;
    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};

    #pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType iL = n[d_hexLRSCV[ip * 2]];
        KeyType iR = n[d_hexLRSCV[ip * 2 + 1]];
        RealType pf = RealType(0.5) * (p[iL] + p[iR]);

        size_t off = e * 12 + ip;
        RealType cx = pf * areaVecX[off];
        RealType cy = pf * areaVecY[off];
        RealType cz = pf * areaVecZ[off];

        atomicAdd(&gxAccNode[iL], +cx);
        atomicAdd(&gyAccNode[iL], +cy);
        atomicAdd(&gzAccNode[iL], +cz);
        atomicAdd(&gxAccNode[iR], -cx);
        atomicAdd(&gyAccNode[iR], -cy);
        atomicAdd(&gzAccNode[iR], -cz);
    }
}

// Path A: literal transpose of the divergence scatter, so that the Chorin
// projection identity D u^{n+1} = D u** - (dt/rho) D D^T phi cancels exactly.
//
// The divergence kernel writes:
//   divAcc[L] += +0.5 * A_f . (u[L] + u[R])
//   divAcc[R] += -0.5 * A_f . (u[L] + u[R])
// Block-transposing the per-face 2x2 contribution gives:
//   gradAcc[L] += +0.5 * A_f * (p[L] - p[R])
//   gradAcc[R] += +0.5 * A_f * (p[L] - p[R])   (same sign as L, NOT opposite)
// That same-sign scatter of dp=0.5*(p[L]-p[R]) is what distinguishes this from
// computeGradientPerNodeKernel, which uses (p[L]+p[R]) with opposite signs and
// is the gradient one would get from independent integration by parts.
template<typename KeyType, typename RealType>
__global__ void applyDivTransposePerNodeKernel(const KeyType* c0, const KeyType* c1,
                                                const KeyType* c2, const KeyType* c3,
                                                const KeyType* c4, const KeyType* c5,
                                                const KeyType* c6, const KeyType* c7,
                                                const RealType* p,
                                                const RealType* areaVecX,
                                                const RealType* areaVecY,
                                                const RealType* areaVecZ,
                                                RealType* gxAccNode,
                                                RealType* gyAccNode,
                                                RealType* gzAccNode,
                                                size_t startElem,
                                                size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;
    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};

    #pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType iL = n[d_hexLRSCV[ip * 2]];
        KeyType iR = n[d_hexLRSCV[ip * 2 + 1]];
        RealType dp = RealType(0.5) * (p[iL] - p[iR]);

        size_t off = e * 12 + ip;
        RealType cx = dp * areaVecX[off];
        RealType cy = dp * areaVecY[off];
        RealType cz = dp * areaVecZ[off];

        atomicAdd(&gxAccNode[iL], cx);
        atomicAdd(&gyAccNode[iL], cy);
        atomicAdd(&gzAccNode[iL], cz);
        atomicAdd(&gxAccNode[iR], cx);
        atomicAdd(&gyAccNode[iR], cy);
        atomicAdd(&gzAccNode[iR], cz);
    }
}

// NOTE: this kernel is currently UNUSED. Per-element assembly of A = D M^-1 D^T
// is fundamentally insufficient because the per-node M^-1 step in
// applyDDTPerNode couples ANY pair of SCS faces (f, g) that share an
// intermediate node i -- regardless of which element they belong to. A
// per-element kernel only sees its OWN 12 faces, so the cross-element
// face-pair contributions to A[r, c] (where r is in element e, c in element
// e', and the pair shares an intermediate node) are silently dropped.
// Validated on 1 rank: assembled-CG converges to a different solution than
// matrix-free CG (which is bit-for-bit correct against the reference).
//
// Correct assembly would require a node-driven kernel that gathers the
// node-to-element-to-face adjacency, then enumerates (g, f) pairs per
// node across all touching elements. Left for later (see TGV_HANDOFF.md);
// production path is matrix-free CG.
//
// The pair-double-sum derivation is kept below for reference. With
// sigma_h(n) = +1 if n==L_h, -1 if n==R_h:
//   A[r,c] = 0.25 * sum_{face pair (g,f)} sigma_g(r) sigma_f(c) (A_g . A_f) *
//            sum_{i in {L_g,R_g} INTERSECT {L_f,R_f}} (1/V[i])
//
// Per-element implementation: enumerate the 8 nodes; for each node i with
// the 3 SCS faces incident on it, enumerate ordered pairs (g, f) of incident
// faces. The pair contributes to A[r, c] for the 2x2 grid (r in {L_g, R_g},
// c in {L_f, R_f}) scaled by 1/V[i]. ONLY captures within-element pairs --
// cross-element pairs sharing an intermediate i across element boundaries
// are NOT captured.
template<typename KeyType, typename RealType>
__global__ void assembleDDTPerElementKernel(const KeyType* c0, const KeyType* c1,
                                            const KeyType* c2, const KeyType* c3,
                                            const KeyType* c4, const KeyType* c5,
                                            const KeyType* c6, const KeyType* c7,
                                            const RealType* areaVecX,
                                            const RealType* areaVecY,
                                            const RealType* areaVecZ,
                                            const int* nodeToDof,
                                            const uint8_t* ownership,
                                            const RealType* lumpedMassNode, // per NODE (incl ghosts)
                                            const int* rowPtr,
                                            const int* colInd,
                                            int numOwnedDofs,
                                            RealType* values,
                                            size_t startElem,
                                            size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;
    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};

    // Per-node incidence: each hex corner touches 3 SCS faces (see d_hexLRSCV).
    // Indices into the 12-face local table. Listed in (face_idx, sign) pairs
    // where sign is sigma_face(node) = +1 if node is L, -1 if node is R.
    // Reading off d_hexLRSCV:
    //   0:(0,1), 1:(1,2), 2:(2,3), 3:(0,3),
    //   4:(4,5), 5:(5,6), 6:(6,7), 7:(4,7),
    //   8:(0,4), 9:(1,5), 10:(2,6), 11:(3,7)
    constexpr int nodeFaces[8][3] = { {0,3,8},  {0,1,9},  {1,2,10}, {2,3,11},
                                       {4,7,8},  {4,5,9},  {5,6,10}, {6,7,11} };

    // Cache the 12 area vectors locally so the inner double-loop reuses them.
    RealType Ax[12], Ay[12], Az[12];
    #pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        size_t off = e * 12 + ip;
        Ax[ip] = areaVecX[off];
        Ay[ip] = areaVecY[off];
        Az[ip] = areaVecZ[off];
    }

    // Iterate intermediate node i in {0..7}; the 1/V[i] factor multiplies the
    // pair contribution. V is read from per-NODE lumped mass (halo-correct).
    #pragma unroll
    for (int iLocal = 0; iLocal < 8; ++iLocal)
    {
        KeyType iGlob = n[iLocal];
        RealType vi = lumpedMassNode[iGlob];
        if (!(vi > RealType(0))) continue;
        RealType invVi = RealType(1) / vi;

        // For each ordered pair (g, f) of faces both incident to node i:
        //   - g contributes to rows r in {L_g, R_g}
        //   - f contributes to cols c in {L_f, R_f}
        // 3*3 = 9 face pairs per intermediate node, 8 nodes -> 72 pairs total
        // per element. Each pair writes a 2x2 block (4 atomic adds), giving
        // 288 atomic adds per element, vs the previous 96 (12 faces * 8).
        #pragma unroll
        for (int gi = 0; gi < 3; ++gi)
        {
            int gFace = nodeFaces[iLocal][gi];
            int gL = d_hexLRSCV[gFace * 2];
            int gR = d_hexLRSCV[gFace * 2 + 1];
            RealType agx = Ax[gFace], agy = Ay[gFace], agz = Az[gFace];

            #pragma unroll
            for (int fi = 0; fi < 3; ++fi)
            {
                int fFace = nodeFaces[iLocal][fi];
                int fL = d_hexLRSCV[fFace * 2];
                int fR = d_hexLRSCV[fFace * 2 + 1];
                RealType afx = Ax[fFace], afy = Ay[fFace], afz = Az[fFace];

                RealType dot   = agx*afx + agy*afy + agz*afz;
                RealType scale = RealType(0.25) * dot * invVi;

                // 2x2 block: rows = {L_g, R_g}, cols = {L_f, R_f}.
                // sigma_g(L_g) = +1, sigma_g(R_g) = -1.
                int   rIdx[2] = {gL, gR};
                int   cIdx[2] = {fL, fR};
                RealType sigR[2] = {RealType(+1), RealType(-1)};
                RealType sigC[2] = {RealType(+1), RealType(-1)};

                #pragma unroll
                for (int rr = 0; rr < 2; ++rr)
                {
                    KeyType nodeR = n[rIdx[rr]];
                    if (ownership[nodeR] != 1) continue;
                    int dofR = nodeToDof[nodeR];
                    if (dofR < 0 || dofR >= numOwnedDofs) continue;
                    int rs = rowPtr[dofR];
                    int re = rowPtr[dofR + 1];

                    #pragma unroll
                    for (int cc = 0; cc < 2; ++cc)
                    {
                        KeyType nodeC = n[cIdx[cc]];
                        int dofC = nodeToDof[nodeC];
                        if (dofC < 0) continue;
                        RealType v = sigR[rr] * sigC[cc] * scale;
                        fem::atomicAddSparseEntry(values, colInd, rs, re, dofC, v);
                    }
                }
            }
        }
    }
}

// Node-driven assembly of A = D M^-1 D^T into a CSR matrix sharing K's
// 27-NNZ full-element sparsity. Replicates matrix-free applyDDTPerNode
// EXACTLY, including cross-element face-pair contributions through the
// per-node M^-1 step.
//
// One CUDA thread per OWNED node i. For each i, we gather ALL SCS faces
// incident to i across every element that touches i (via the node->element
// CSR adjacency map d_nodeToElementOffsets/d_nodeToElementList). For each
// touching element e, we find which local corner of e is node i (search
// the 8 connectivity entries), then look up the 3 incident SCS faces of
// e at that corner and append each as a face record (other_node, area,
// sigma) to a per-thread list.
//
// Then we double-loop over the ordered (g, f) pairs of incident faces and
// contribute
//   delta_A[i, c] = 0.25 * sigma_g(i) * sigma_f(c) * (A_g . A_f) / V[i]
// for c in {i, other_f}. Row "other_g" is NOT written here -- it is owned
// by another node thread, which will fill that row from its own
// perspective. Each owned row is therefore written by exactly one thread.
//
// This is exactly the analytical expansion of the matrix-free chain
// D^T (scatter), M^-1 (per-node divide), D (scatter) applied to a unit
// vector e_c -- so the assembled matrix is bit-identical to the operator.
//
// Hex8 specific: max 8 elements per corner node * 3 SCS faces per element
// per corner = 24 faces in the worst case. We store them in a fixed-size
// register array (MAX_FACES = 32 for safety; nodes with more incidents --
// e.g., from non-conformal AMR -- would need a fallback, but plain hex
// meshes never exceed 24).
template<typename KeyType, typename RealType>
__global__ void assembleDDTPerNodeKernel(const KeyType* c0, const KeyType* c1,
                                         const KeyType* c2, const KeyType* c3,
                                         const KeyType* c4, const KeyType* c5,
                                         const KeyType* c6, const KeyType* c7,
                                         const RealType* areaVecX,
                                         const RealType* areaVecY,
                                         const RealType* areaVecZ,
                                         const KeyType* nodeToElemOffsets,
                                         const KeyType* nodeToElemList,
                                         const int* nodeToDof,
                                         const uint8_t* ownership,
                                         const RealType* lumpedMassNode,
                                         const int* rowPtr,
                                         const int* colInd,
                                         int numOwnedDofs,
                                         RealType* values,
                                         size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dofI = nodeToDof[i];
    if (dofI < 0 || dofI >= numOwnedDofs) return;

    RealType vi = lumpedMassNode[i];
    if (!(vi > RealType(0))) return;
    RealType invVi = RealType(1) / vi;

    // Per-node incidence: which 3 of the 12 local SCS faces touch each
    // local corner.
    constexpr int nodeFaces[8][3] = { {0,3,8},  {0,1,9},  {1,2,10}, {2,3,11},
                                       {4,7,8},  {4,5,9},  {5,6,10}, {6,7,11} };

    constexpr int MAX_FACES = 32;
    int      incCount = 0;
    KeyType  incOtherNode[MAX_FACES];
    int      incOtherDof[MAX_FACES];
    RealType incAx[MAX_FACES], incAy[MAX_FACES], incAz[MAX_FACES];
    int8_t   incSigma[MAX_FACES];

    KeyType eStart = nodeToElemOffsets[i];
    KeyType eEnd   = nodeToElemOffsets[i + 1];
    for (KeyType ep = eStart; ep < eEnd; ++ep)
    {
        KeyType e = nodeToElemList[ep];
        KeyType en[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};
        int iLocal = -1;
        #pragma unroll
        for (int k = 0; k < 8; ++k) if (en[k] == (KeyType)i) iLocal = k;
        if (iLocal < 0) continue;

        #pragma unroll
        for (int gi = 0; gi < 3; ++gi)
        {
            if (incCount >= MAX_FACES) break;
            int gFace = nodeFaces[iLocal][gi];
            int gL_local = d_hexLRSCV[gFace * 2];
            int gR_local = d_hexLRSCV[gFace * 2 + 1];
            int sigma = (iLocal == gL_local) ? +1 : -1;
            int otherLocal = (iLocal == gL_local) ? gR_local : gL_local;
            KeyType otherNode = en[otherLocal];
            size_t off = e * 12 + gFace;

            incOtherNode[incCount] = otherNode;
            incOtherDof[incCount]  = nodeToDof[otherNode];
            incAx[incCount]        = areaVecX[off];
            incAy[incCount]        = areaVecY[off];
            incAz[incCount]        = areaVecZ[off];
            incSigma[incCount]     = int8_t(sigma);
            incCount++;
        }
    }

    // Row r=i is owned by THIS thread; write that row directly.
    int rs_i = rowPtr[dofI];
    int re_i = rowPtr[dofI + 1];

    // Double-loop over face pairs (g, f) at intermediate node i.
    // The 2x2 contribution to A is at rows r in {i, other_g}, cols c in
    // {i, other_f} with signs sigma_g(r) * sigma_f(c).
    //   sigma_g(i)        = sg_i        sigma_g(other_g) = -sg_i
    //   sigma_f(i)        = sf_i        sigma_f(other_f) = -sf_i
    // We must write all four cells because no other node's thread will
    // ever visit this pair at intermediate i (other_g's thread iterates
    // pairs at intermediate other_g, not at i). Writes to "other" rows
    // gate on ownership[other_g] == 1.
    for (int g = 0; g < incCount; ++g)
    {
        RealType agx = incAx[g], agy = incAy[g], agz = incAz[g];
        int sg_i = incSigma[g];
        KeyType otherG_node = incOtherNode[g];
        // Look up row info for other_g if owned.
        int otherG_dof = -1;
        int rs_og = 0, re_og = 0;
        if (ownership[otherG_node] == 1)
        {
            otherG_dof = nodeToDof[otherG_node];
            if (otherG_dof >= 0 && otherG_dof < numOwnedDofs)
            {
                rs_og = rowPtr[otherG_dof];
                re_og = rowPtr[otherG_dof + 1];
            }
            else otherG_dof = -1;
        }

        for (int f = 0; f < incCount; ++f)
        {
            RealType dot   = agx * incAx[f] + agy * incAy[f] + agz * incAz[f];
            RealType scale = RealType(0.25) * dot * invVi;
            int sf_i = incSigma[f];
            int otherF_dof = incOtherDof[f];
            RealType base = RealType(sg_i) * RealType(sf_i) * scale;
            // (r=i,         c=i):         + base
            // (r=i,         c=other_f):   - base   (sigma_f(other_f)=-sf_i)
            // (r=other_g,   c=i):         - base   (sigma_g(other_g)=-sg_i)
            // (r=other_g,   c=other_f):   + base
            fem::atomicAddSparseEntry(values, colInd, rs_i, re_i, dofI,  base);
            if (otherF_dof >= 0)
                fem::atomicAddSparseEntry(values, colInd, rs_i, re_i, otherF_dof, -base);
            if (otherG_dof >= 0)
            {
                fem::atomicAddSparseEntry(values, colInd, rs_og, re_og, dofI, -base);
                if (otherF_dof >= 0)
                    fem::atomicAddSparseEntry(values, colInd, rs_og, re_og, otherF_dof, base);
            }
        }
    }
}

// Divide accumulator by V; write to per-node output (kept node-indexed for the
// predictor/corrector apply, which iterates over nodes).
template<typename RealType>
__global__ void normalizeGradientPerNodeKernel(const RealType* gxAccNode,
                                               const RealType* gyAccNode,
                                               const RealType* gzAccNode,
                                               const RealType* lumpedMass,
                                               const int* nodeToDof,
                                               const uint8_t* ownership,
                                               RealType* gx,
                                               RealType* gy,
                                               RealType* gz,
                                               size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    int dof = (ownership[i] == 1) ? nodeToDof[i] : -1;
    if (dof < 0) { gx[i] = gy[i] = gz[i] = RealType(0); return; }
    RealType invV = RealType(1) / lumpedMass[dof];
    gx[i] = gxAccNode[i] * invV;
    gy[i] = gyAccNode[i] * invV;
    gz[i] = gzAccNode[i] * invV;
}

// =============================================================================
// Divergence scatter (B.3 pattern): integrates SCS-face fluxes
// v_f . A_f into the L slot, -v_f . A_f into the R slot. Used in step 3 to
// build the pressure-Poisson RHS: b = (rho/dt) * div(u**) * V * sign? Below
// we keep the scatter in the un-normalized accumulator form because the RHS
// for the Poisson system is the integrated source f*V_i (lumped), which is
// exactly what the un-normalized accumulator already is when multiplied by
// the constant (rho/dt). Skipping the normalize-then-multiply-by-V cancels out.
// =============================================================================

template<typename KeyType, typename RealType>
__global__ void computeDivergencePerNodeKernel(const KeyType* c0, const KeyType* c1,
                                               const KeyType* c2, const KeyType* c3,
                                               const KeyType* c4, const KeyType* c5,
                                               const KeyType* c6, const KeyType* c7,
                                               const RealType* vx,
                                               const RealType* vy,
                                               const RealType* vz,
                                               const RealType* areaVecX,
                                               const RealType* areaVecY,
                                               const RealType* areaVecZ,
                                               RealType* divAccNode,
                                               size_t startElem,
                                               size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;
    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};

    #pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        int nodeL = d_hexLRSCV[ip * 2];
        int nodeR = d_hexLRSCV[ip * 2 + 1];
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];

        RealType vfx = RealType(0.5) * (vx[iL] + vx[iR]);
        RealType vfy = RealType(0.5) * (vy[iL] + vy[iR]);
        RealType vfz = RealType(0.5) * (vz[iL] + vz[iR]);

        size_t off = e * 12 + ip;
        RealType flow = vfx * areaVecX[off] + vfy * areaVecY[off] + vfz * areaVecZ[off];

        atomicAdd(&divAccNode[iL], +flow);
        atomicAdd(&divAccNode[iR], -flow);
    }
}

// =============================================================================
// Rhie-Chow corrected divergence kernel for periodic Q1-Q1 CVFEM. Same as
// computeDivergencePerNodeKernel but adds the Rhie-Chow pressure-velocity
// coupling correction at each SCS face:
//
//   mdot_RC = avg(u) . A - tau * (p_R - p_L) * |A|^2 / (A . dx)
//                        + tau * avg(grad_p) . A
//
// The middle term is the "compact" face-gradient stencil; the third is the
// "smooth" nodal-gradient average. Their difference IS the Rhie-Chow
// stabilization that suppresses pressure-velocity decoupling on co-located
// grids. The formula matches Nalu-Wind's MdotEdgeAlg and OpenFOAM's PISO
// pEqn.flux() construction.
//
// tau ~ dt/rho for explicit Chorin; for more accurate variants, use the
// inverse momentum-matrix diagonal 1/A_diag (the SIMPLE/PISO weight).
// =============================================================================

template<typename KeyType, typename RealType>
__global__ void computeDivergenceRhieChowKernel(const KeyType* c0, const KeyType* c1,
                                                const KeyType* c2, const KeyType* c3,
                                                const KeyType* c4, const KeyType* c5,
                                                const KeyType* c6, const KeyType* c7,
                                                const RealType* vx,
                                                const RealType* vy,
                                                const RealType* vz,
                                                const RealType* p,
                                                const RealType* gradPx,
                                                const RealType* gradPy,
                                                const RealType* gradPz,
                                                const RealType* nodeX,
                                                const RealType* nodeY,
                                                const RealType* nodeZ,
                                                const RealType* areaVecX,
                                                const RealType* areaVecY,
                                                const RealType* areaVecZ,
                                                RealType tauRC,
                                                RealType* divAccNode,
                                                size_t startElem,
                                                size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;
    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};

    #pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        int nodeL = d_hexLRSCV[ip * 2];
        int nodeR = d_hexLRSCV[ip * 2 + 1];
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];

        // Standard average face velocity.
        RealType vfx = RealType(0.5) * (vx[iL] + vx[iR]);
        RealType vfy = RealType(0.5) * (vy[iL] + vy[iR]);
        RealType vfz = RealType(0.5) * (vz[iL] + vz[iR]);

        size_t off = e * 12 + ip;
        RealType Ax = areaVecX[off];
        RealType Ay = areaVecY[off];
        RealType Az = areaVecZ[off];

        RealType flow = vfx * Ax + vfy * Ay + vfz * Az;

        // Node-pair geometry: dx = position(R) - position(L)
        RealType dx = nodeX[iR] - nodeX[iL];
        RealType dy = nodeY[iR] - nodeY[iL];
        RealType dz = nodeZ[iR] - nodeZ[iL];

        RealType asq  = Ax * Ax + Ay * Ay + Az * Az;
        RealType axdx = Ax * dx + Ay * dy + Az * dz;

        // CHORIN Rhie-Chow: compact-pressure-gradient term ONLY.
        //
        // In a Chorin projection, the predictor `u** = u^n + dt*adv - dt*grad(p^n)/rho`
        // has already subtracted the smooth nodal grad(p^n) from velocity, so
        // div(u**) already implicitly includes the smooth-gradient pressure
        // contribution. The Nalu/OpenFOAM RC formula `-compact + smooth` is
        // appropriate for SIMPLE/PISO loops where the predicted velocity
        // u* = HbyA does NOT yet contain grad p. Applying that form to a
        // post-predictor u** double-counts the smooth contribution and
        // produces catastrophic divergence (verified empirically).
        //
        // The Code_Saturne `arak`-only form and the OpenFOAM `phi = phiHbyA
        // - pEqn.flux()` pattern both use the compact term alone:
        //
        //   mdot_RC = avg(u**)*A - tau * (p_R - p_L) * |A|^2 / (A.dx)
        //
        // This subtracts only the compact-vs-smooth DISCREPANCY (which is
        // exactly the checkerboard-mode amplitude when p oscillates by sign
        // at adjacent nodes), leaving smooth pressure fields unaffected.
        if (fabs(axdx) > RealType(1e-30))
        {
            RealType compact = tauRC * (p[iR] - p[iL]) * asq / axdx;
            flow += -compact;
        }
        (void)gradPx; (void)gradPy; (void)gradPz;  // preserved for future variants

        atomicAdd(&divAccNode[iL], +flow);
        atomicAdd(&divAccNode[iR], -flow);
    }
}

// Normalized divergence (V_i^{-1} * sum), per node. Used for the div_max
// diagnostic only -- the Poisson RHS uses the un-normalized accumulator.
template<typename RealType>
__global__ void normalizeDivergencePerNodeKernel(const RealType* divAccNode,
                                                 const RealType* lumpedMass,
                                                 const int* nodeToDof,
                                                 const uint8_t* ownership,
                                                 RealType* divU,
                                                 size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) { divU[i] = RealType(0); return; }
    int dof = nodeToDof[i];
    if (dof < 0) { divU[i] = RealType(0); return; }
    divU[i] = divAccNode[i] / lumpedMass[dof];
}

// =============================================================================
// Predictor / corrector apply kernels.
// =============================================================================

// Predictor: q* = q^n + dt * dqdtNode / V - dt * (1/rho) * gradPnq on interior;
// on boundary, snap to qTarget (so the implicit RHS lifts it without drift).
// dqdtNode contains the per-node UN-normalized flux sum from the advection
// scatter; we divide by V here (no separate normalize step). gradPnq is the
// per-node gradient component (already V-normalized by the gradient kernel).
template<typename RealType>
__global__ void applyPredictorPerNodeKernel(RealType* qStar,
                                            const RealType* qN,
                                            const RealType* dqdtNode,
                                            const RealType* gradPnq,
                                            const RealType* lumpedMass,
                                            const RealType* /*lumpedMassNode*/,  // unused; kept for ABI
                                            const RealType* qTarget,
                                            const uint8_t* isBdryDof,
                                            const int* nodeToDof,
                                            const uint8_t* ownership,
                                            RealType dt,
                                            RealType invRho,
                                            size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;

    if (isBdryDof[dof])
    {
        qStar[i] = qTarget[i];
        return;
    }
    RealType V = lumpedMass[dof];
    qStar[i] = qN[i] + dt * dqdtNode[i] / V - dt * invRho * gradPnq[i];
}

// BDF2 + EXT2 predictor: 3-point backward time derivative + Adams-Bashforth-
// style extrapolation 2 N(u^n) - N(u^(n-1)) for advection.
//   u* = (4/3) u^n - (1/3) u^(n-1)
//        + (2 dt / 3) * [ (2 N_n - N_nm1)/V - invRho * grad p^n ]
// Implements the Karniadakis-Israeli-Orszag stencil also used by deal.II
// step-35, MFEM-Navier, and NekRS.
template<typename RealType>
__global__ void applyPredictorBdf2PerNodeKernel(RealType* qStar,
                                                const RealType* qN,
                                                const RealType* qNm1,
                                                const RealType* dqdtNode_n,
                                                const RealType* dqdtNode_nm1,
                                                const RealType* gradPnq,
                                                const RealType* lumpedMass,
                                                const RealType* qTarget,
                                                const uint8_t* isBdryDof,
                                                const int* nodeToDof,
                                                const uint8_t* ownership,
                                                RealType dt,
                                                RealType invRho,
                                                size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    if (isBdryDof[dof])
    {
        qStar[i] = qTarget[i];
        return;
    }
    RealType V       = lumpedMass[dof];
    RealType adv_ext = RealType(2) * dqdtNode_n[i] - dqdtNode_nm1[i];
    RealType coef    = RealType(2) * dt / RealType(3);
    qStar[i] = (RealType(4) / RealType(3)) * qN[i]
             - (RealType(1) / RealType(3)) * qNm1[i]
             + coef * (adv_ext / V - invRho * gradPnq[i]);
}

// Corrector: q^{n+1} = q** - (dt/rho) (grad phi)_q on interior; snap to
// qTarget on boundary (Dirichlet velocity unchanged by the corrector since
// grad(phi).n = 0 is enforced on impermeable walls -- here we just lock it).
template<typename RealType>
__global__ void applyCorrectorPerNodeKernel(RealType* q,
                                            const RealType* qStarStar,
                                            const RealType* gradPhiq,
                                            const RealType* qTarget,
                                            const uint8_t* isBdryDof,
                                            const int* nodeToDof,
                                            const uint8_t* ownership,
                                            RealType dt,
                                            RealType invRho,
                                            size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;

    if (isBdryDof[dof])
    {
        q[i] = qTarget[i];
        return;
    }
    q[i] = qStarStar[i] - dt * invRho * gradPhiq[i];
}

// Pressure update: p^{n+1} = p^n + phi on owned nodes (Chorin incremental).
// Ghosts are refreshed by an explicit halo exchange afterward.
template<typename RealType>
__global__ void updatePressureKernel(RealType* p,
                                     const RealType* phi,
                                     const int* nodeToDof,
                                     const uint8_t* ownership,
                                     size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    p[i] += phi[i];
}

// =============================================================================
// RHS builders for the implicit solves.
// =============================================================================

// Velocity implicit RHS: b[dof] = (M[dof]/dt) * qStar[node(dof)] on owned
// nodes. Same shape as mars_amr_advdiff's buildRhsFromTimeNKernel; called
// three times (one per component) with the appropriate qStar.
template<typename RealType>
__global__ void buildVelocityImplicitRhsKernel(const RealType* qStar,
                                               const int* nodeToDof,
                                               const uint8_t* ownership,
                                               const RealType* mass,
                                               RealType invDt,
                                               RealType* rhs,
                                               size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    rhs[dof] = mass[dof] * invDt * qStar[i];
}

// Pressure RHS from the un-normalized divergence accumulator: the per-node
// scatter sum IS the integrated flux (no extra /V).
//
// Sign convention: CvfemHexAssembler::assembleFull's K is the discrete operator
// for -Laplacian (see mars_amr_poisson / mars_cvfem_poisson: "Solve -Δu = f").
// The continuous projection equation is +Lap(phi) = (rho/dt) div(u**), so the
// discrete system is K phi = -(rho/dt) * div(u**)_lumped. divAccNode[i] is the
// lumped integral V_i * div(u**)_i, so rhs[dof] = -coef * divAccNode[i].
template<typename RealType>
__global__ void buildPressureRhsKernel(const RealType* divAccNode,
                                       const int* nodeToDof,
                                       const uint8_t* ownership,
                                       RealType coef,        // = rho / dt
                                       RealType* rhs,
                                       size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    rhs[dof] = -coef * divAccNode[i];
}

// DOF-indexed solver output -> per-node array. Reused for all velocity solves
// and the pressure solve.
template<typename RealType>
__global__ void scatterDofToNodeKernel(const RealType* sol,
                                       const int* nodeToDof,
                                       RealType* nodeOut,
                                       size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    int dof = nodeToDof[i];
    if (dof >= 0) nodeOut[i] = sol[dof];
}

// =============================================================================
// Solver kind + Stepper state. All persistent state lives here -- the per-step
// function reads / writes the device vectors in this struct and never allocates
// in the hot path (modulo per-call temporaries for the implicit solves).
// =============================================================================

enum class SolverKind { CG, Hypre };

// Krylov method hint for solveOneComponent. PCG is the default for velocity
// (Avel) and K-path pressure (Apre). GMRES is required for the DDT pressure
// operator (D M^-1 D^T), which is SPSD with a constant null space and is
// rejected by Hypre PCG (error 1). Only consulted when s.solverKind == Hypre;
// the in-house CG path ignores it.
//
// Env MARS_HYPRE_KRYLOV=pcg|gmres (DEBUG ONLY) overrides the caller's hint.
enum class KrylovHint { PCG, GMRES };

// Phase E: matrix-free pressure-Poisson choice.
//   K   -> existing CVFEM stiffness matrix K phi = -(rho/dt) D u**  (Galerkin
//          Laplacian; not discretely D D^T, so D u^{n+1} only decays slowly).
//   DDT -> matrix-free (D M^{-1} D^T) phi = (rho/dt) D u** ; gives exact
//          projection D u^{n+1} = 0 by construction because the corrector's
//          (1/rho) D^T phi cancels D u** algebraically through that operator.
// PressureSolveKind:
//   K   = (default, stable) solve K*phi = rhs with assembled stiffness from
//         assembleFull. Produces a stable, well-behaved cavity solve but
//         div_max decays slowly toward zero because K is the continuous
//         Laplacian discretization, not the discrete D D^T conjugate of
//         this code's divergence operator.
//   DDT = Matrix-free Galerkin Laplacian A = D M^{-1} D^T. The corrector's
//         grad operator uses M^{-1} D^T to make the projection algebraically
//         exact: div(u^{n+1}) = roundoff every step. ~3-4x slower per CG iter
//         vs K-path; iteration count similar on uniform meshes (~50-200 on
//         cube16). Production preconditioning (multigrid or Hypre on assembled
//         D D^T) is a separate work item.
enum class PressureSolveKind { K, DDT };

template<typename KeyType, typename RealType> struct NSStepper;

// Forward declaration so the setup-time MARS_DDT_PROBE_DIFF diagnostic can
// call applyDDTPerNode (defined later in this file).
template<typename KeyType, typename RealType>
void applyDDTPerNode(NSStepper<KeyType, RealType>& s,
                     cstone::DeviceVector<RealType>& phi,
                     cstone::DeviceVector<RealType>& outAcc,
                     cstone::DeviceVector<RealType>& gxAcc,
                     cstone::DeviceVector<RealType>& gyAcc,
                     cstone::DeviceVector<RealType>& gzAcc);

template<typename KeyType, typename RealType>
struct NSStepper
{
    using DomainT = ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>;

    DomainT& domain;
    SolverKind solverKind;
    int blockSize;
    int maxIter;
    RealType tolerance;
    int rank;
    int numRanks;
    PressureSolveKind pressureSolve = PressureSolveKind::K;
    // Path A (div^T gradient) was attempted but is unstable: the CVFEM stiffness
    // matrix K assembled by assembleFull discretizes the continuous Laplacian,
    // not D D^T from this scatter pattern. So phi solved against K is not the
    // phi that D^T would cancel divergence with -> projection diverges.
    // Default is the SCS gradient (legacy); the div^T path stays available
    // under --experimental-divT for future Path B (Rhie-Chow) work.
    bool useLegacyGradient = true;

    size_t nodeCount    = 0;
    size_t elementCount = 0;
    int numOwnedDofs    = 0;
    int numTotalDofs    = 0;
    int nnz             = 0;

    // DOF mapping + lumped mass + BC mask (per-OWNED-DOF for the mask).
    cstone::DeviceVector<int>     d_node_to_dof;
    cstone::DeviceVector<RealType> d_mass;       // per OWNED DOF
    cstone::DeviceVector<RealType> d_massNode;   // per NODE (owned + ghosts halo-exchanged); read by DDT assembler
    cstone::DeviceVector<uint8_t>  d_isBdryDof;
    // Pressure BC: in cavity mode, use a single corner pin (pressurePinDof).
    // In channel mode, use a Dirichlet mask (d_isPressureBdryDof) over the
    // entire outflow face. Exactly one mechanism is active per run.
    int pressurePinDof = -1;   // owned-DOF id on the owning rank; -1 elsewhere (or in channel mode)
    int pressurePinRank = 0;   // global rank that owns the pin
    cstone::DeviceVector<uint8_t> d_isPressureBdryDof;   // size numOwnedDofs; only used in channel mode
    // Per-NODE pressure-Dirichlet mask (size nodeCount, halo-exchanged so
    // ghost slots see the owner's flag). Includes the single corner pin
    // (cavity) AND every outflow-face DOF (channel/wing). Used by symmetric
    // column-zeroing on the assembled DDT CSR so non-owning ranks also zero
    // out their ghost-copy columns; without this, the cavity pin was only
    // enforced symmetrically on rank 0 and the channel mask only on each
    // owner, leaving the global matrix asymmetric and AMG misbehaving.
    cstone::DeviceVector<uint8_t> d_isPressureBdryNode;
    // Inverse of nodeToDof: dofToNode[dof] = local-node-index. Size numTotalDofs
    // (= nodeCount). Used by symmetric col-zero kernel to map a column local
    // DOF back to its local node, so the per-NODE mask can be consulted.
    cstone::DeviceVector<int> d_dofToNode;

    // Per-node velocity targets for the cavity BC. Boundary nodes hold the
    // prescribed value (u=1 on top, 0 elsewhere); interior nodes hold 0 (unused).
    cstone::DeviceVector<RealType> d_uTarget, d_vTarget, d_wTarget;

    // Geometry: SCS area vectors per (element, face). Same source used by the
    // implicit assembler, the advection scatter, the gradient, the divergence.
    cstone::DeviceVector<RealType> d_areaVec_x;
    cstone::DeviceVector<RealType> d_areaVec_y;
    cstone::DeviceVector<RealType> d_areaVec_z;

    // Two CSR matrices: velocity (M/dt + nu K) with Dirichlet velocity BC;
    // pressure (K) with pure Neumann + a single-DOF pin. Share rowPtr/colInd/
    // diagPtr (full 27-NNZ Laplacian sparsity required by assembleFull);
    // values are stored separately.
    cstone::DeviceVector<int> d_rowPtr;
    cstone::DeviceVector<int> d_colInd;
    cstone::DeviceVector<int> d_diagPtr;
    cstone::DeviceVector<RealType> d_valuesVel;
    cstone::DeviceVector<RealType> d_valuesPre;
    // Assembled D M^-1 D^T uses its OWN reduced (7-NNZ) sparsity that exactly
    // matches the SCS-face graph the assembler writes to. Sharing K's 27-NNZ
    // sparsity left ~20 explicit-zero entries per row that confused Hypre
    // BoomerAMG's coarsening (the stored zeros made the graph denser than the
    // operator). Separate CSR avoids that. Diagnosed via parallel research
    // agents during Phase E Stage 4 debug.
    cstone::DeviceVector<int> d_rowPtrDDT;
    cstone::DeviceVector<int> d_colIndDDT;
    cstone::DeviceVector<int> d_diagPtrDDT;
    int nnzDDT = 0;
    cstone::DeviceVector<RealType> d_valuesDDT;

    // Owned-row SparseMatrix wrappers consumed by CG/Hypre.
    using Matrix = SparseMatrix<int, RealType, cstone::GpuTag>;
    Matrix Avel;
    Matrix Apre;
    Matrix AddT;

    // Per-node solution fields. Sized nodeCount; ghost slots refreshed by
    // exchangeNodeHalo after each update.
    cstone::DeviceVector<RealType> d_u, d_v, d_w;            // velocity
    cstone::DeviceVector<RealType> d_uStar, d_vStar, d_wStar;        // after predictor
    cstone::DeviceVector<RealType> d_uStarStar, d_vStarStar, d_wStarStar;  // after diffusion
    cstone::DeviceVector<RealType> d_p;                       // pressure (n)
    cstone::DeviceVector<RealType> d_phi;                     // pressure correction
    // Stash of normalized div(u**) per node, captured during the pressure
    // solve, reused by the rotational pressure correction: p = p + phi - nu*div(u**).
    // Sized lazily in runPressureSolveStep.
    cstone::DeviceVector<RealType> d_divUStar;
    bool rotationalPressureCorrection = false;
    RealType nuCached = 0;  // remember nu so runCorrectorStep can apply -nu*div(u*)
    // Advection form: false = 1st-order upwind (default, stable at any Re; the
    // form used by cavity/channel). true = central skew-symmetric (discrete KE
    // conservation; needed for periodic NS where there is no wall dissipation).
    bool useSkewSymmetricAdvection = false;
    cstone::DeviceVector<RealType> d_gradPx, d_gradPy, d_gradPz;
    cstone::DeviceVector<RealType> d_gradPhix, d_gradPhiy, d_gradPhiz;

    // BDF2 / EXT2 (Karniadakis-Israeli-Orszag 1991) history.
    // ----------------------------------------------------------
    // Standard Chorin is first-order in time + forward-Euler explicit on
    // advection. For periodic NS this is unstable: every "fix" we tried
    // (Bochev-Dohrmann, Rhie-Chow, broadcasts, skew-symmetric stencil)
    // converges to the same ~step 100 blowup. The fundamental limit is
    // the time integrator. BDF2/EXT2 is what NekRS, MFEM-Navier, and
    // deal.II step-35 all use.
    //
    // Formula (per-component scalar, kinematic pressure):
    //   predictor (step n+1, n>=1):
    //     u* = (4/3) u^n  -  (1/3) u^(n-1)
    //          - (2 dt / 3) * [ N_ext + grad(p^n)/rho ]
    //     N_ext = 2 N(u^n) - N(u^(n-1))    // EXT2 extrapolation
    //   diffusion:  (3 M / (2 dt) + nu K) u** = (3 M / (2 dt)) u*
    //   pressure:   K phi = (3 rho / (2 dt)) div(u**)
    //   corrector:  u^(n+1) = u** - (2 dt / (3 rho)) grad(phi)
    //               p^(n+1) = p^n + phi    (incremental Chorin)
    //
    // Step 0 (bootstrap): no u^(n-1) yet -> fall back to BDF1 (the
    // existing Chorin formulas). At step 1, switch to BDF2 with the
    // u^0 stored as u^(n-1).
    bool useBdf2 = true;
    int  bdfStep = 0;     // 0 = use BDF1 (no history yet); >=1 = use BDF2
    // Velocity history (u at step n-1). Sized lazily by setupNSStepper.
    cstone::DeviceVector<RealType> d_u_nm1, d_v_nm1, d_w_nm1;
    // Advection accumulator history: N(u^n) and N(u^(n-1)) per node.
    // Promoted from step-local scratch to persistent so EXT2 can reuse.
    // Sized lazily by setupNSStepper.
    cstone::DeviceVector<RealType> d_advU_n,   d_advV_n,   d_advW_n;
    cstone::DeviceVector<RealType> d_advU_nm1, d_advV_nm1, d_advW_nm1;
    // Second velocity matrix for BDF2 (coefficient 3M/(2dt) + nu K).
    // Built once at setup alongside the BDF1 matrix; used from step 1 on.
    cstone::DeviceVector<RealType> d_valuesVel_bdf2;
    Matrix Avel_bdf2;

    // Per-step diagnostics
    int lastPressureIters = 0;
    int lastUIters = 0, lastVIters = 0, lastWIters = 0;
    RealType lastDivMax     = 0;  // |div(u^{n+1})| max -- post-corrector
    RealType lastDivRms     = 0;  // |div(u^{n+1})| RMS over interior owned DOFs; less ring-sensitive
    RealType lastDivMaxPre  = 0;  // |div(u**)|    max -- pre-corrector (= b magnitude / V scaled)
    RealType lastGradPRms   = 0;  // RMS of grad(p^n) magnitude  -- predictor input
    RealType lastGradPhiRms = 0;  // RMS of grad(phi)  magnitude -- corrector input

    // Hypre-specific persistent state (only filled when solverKind == Hypre).
    int64_t globalRowStart = 0;
    int64_t globalRowEnd = 0;
    int64_t numInteriorGlobal = 0;
#ifdef MARS_ENABLE_HYPRE
    thrust::device_vector<HYPRE_BigInt> d_localToGlobalDof;
#endif

    // Cached bbox for BC marking (constant on a fixed mesh).
    RealType xmin = 0, xmax = 1, ymin = 0, ymax = 1, zmin = 0, zmax = 1;
    RealType bboxEps = 0;
    RealType lidU = 1;

    // BC configuration. Cavity: Dirichlet u=lidU on top face, u=0 on others.
    // Channel: Dirichlet u=Uinf on x=xmin inflow, u=v=w=0 on y/z walls,
    // natural (Neumann) on x=xmax outflow.
    // Wing: per-side-set Dirichlet driven by Exodus side-sets. blade=no-slip
    //   (u=v=w=0), in=Uinf, out=Dirichlet p=0, other named side-sets get
    //   full Dirichlet u=Uinf (free-stream tunnel walls; true slip BCs are
    //   a follow-up). Side-set local-node lists must be populated before
    //   setupNSStepper(). NSStepper does not allocate or free them.
    // Periodic: triply periodic box. No Dirichlet on velocity. Pressure
    // null-space removed by global-mean subtraction every step. Requires
    // periodicMap to be set before setupNSStepper().
    enum class BCKind { Cavity, Channel, Wing, Periodic };
    BCKind bcKind = BCKind::Cavity;
    RealType Uinf = 1;   // channel/wing inflow speed; reuses lidU value via CLI

    // Wing-mode side-set local-node lists. Each vector holds owned-or-ghost
    // LOCAL node indices (size_t) that the driver mapped from global Exodus
    // node IDs via cstone's d_localToGlobalNodeMap_. Populated only when
    // bcKind == Wing. Names match the Exodus mesh side-set names.
    std::vector<int> wingBladeNodes;    // no-slip
    std::vector<int> wingInletNodes;    // u = Uinf
    std::vector<int> wingOutletNodes;   // Dirichlet p = 0
    std::vector<int> wingFarFieldNodes; // top/bottom/side/sym -- u = Uinf (full Dirichlet shortcut)

    // Optional. When bcKind == Periodic, must point to a built PeriodicMap.
    // Owned by the driver; NSStepper does not allocate or free it. The map
    // must be rebuilt by the driver after any AMR step (coords/node count
    // change) before the next runNsStep call.
    const PeriodicMap<KeyType, RealType>* periodicMap = nullptr;

    // Bochev-Dohrmann polynomial pressure projection stabilization.
    // Equal-order CG-FEM (Q1-Q1 here) does NOT satisfy the LBB / inf-sup
    // condition: the discrete pressure space admits spurious checkerboard
    // modes that the discrete gradient cannot detect. In domains with
    // Dirichlet pressure boundaries those modes are pinned by the BC; in
    // *periodic* domains they have no anchor and grow exponentially via the
    // grad(p) feedback in the predictor. Bochev & Dohrmann (2006, SISC) add
    // a consistent element-local stabilization term S = tau (q - Pi q, p - Pi p)
    // where Pi is the L2 projection onto element constants (P0 on Q1). For
    // Q1 hexes with lumped mass V_e/8 per corner, S^e_ij = tau (M^e_ij - V_e/64)
    // -- a small, symmetric element matrix added to the pressure system. This
    // controls the pressure null-space without artificially increasing the
    // velocity dissipation. tau = h^2 / (4 mu) is the canonical choice; we
    // expose it as a knob via stabPressureTau (set <= 0 to auto-pick).
    bool stabBochevDohrmann = false;
    RealType stabPressureTau = -1;
    // Rhie-Chow pressure-velocity coupling on the divergence operator.
    // The textbook CVFEM/co-located-FV stabilization (Nalu-Wind, OpenFOAM,
    // Fluent all use it). For Q1-Q1 on periodic boundaries this is the
    // correct fix for the inf-sup pathology; alternative FEM stabilizations
    // (Bochev-Dohrmann, PSPG, Codina) are FEM-canonical answers to the same
    // question but don't fit a face-flux discretization.
    bool useRhieChow      = false;
    RealType rhieChowTau  = -1;   // <=0 => auto-pick = dt/rho
};

// =============================================================================
// Helpers for the two-matrix Laplacian assembly. Both matrices want the same
// CVFEM Laplacian K (gamma=1), assembled by CvfemHexAssembler::assembleFull
// with zeroed advection inputs. We assemble K once into a temporary buffer,
// then copy it into d_valuesVel (and add M/dt + enforce velocity BC) and into
// d_valuesPre (and enforce single-pin BC) so the two matrices end up with
// different values but the same row/column structure.
// =============================================================================

template<typename KeyType, typename RealType>
void assembleLaplacian(NSStepper<KeyType, RealType>& s,
                       const cstone::DeviceVector<int>& /*unused*/,
                       cstone::DeviceVector<RealType>& d_valuesOut,
                       CvfemKernelVariant kernelVariant)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();

    using MatrixT = CSRMatrix<RealType>;
    MatrixT* d_matrix;
    cudaMalloc(&d_matrix, sizeof(MatrixT));
    MatrixT h_matrix{s.d_rowPtr.data(), s.d_colInd.data(), d_valuesOut.data(), s.d_diagPtr.data(),
                     s.numTotalDofs, s.nnz};
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixT), cudaMemcpyHostToDevice);

    // gamma=1 -> pure Laplacian; all advection inputs zero so assembler doesn't
    // accumulate an advective contribution. RHS unused.
    cstone::DeviceVector<RealType> d_gamma(s.nodeCount, RealType(1));
    cstone::DeviceVector<RealType> d_phi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_beta(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_gphi_x(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_gphi_y(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_gphi_z(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_mdot_zero(s.elementCount * 12, RealType(0));
    cstone::DeviceVector<RealType> d_rhs_unused(s.numTotalDofs, RealType(0));

    typename CvfemHexAssembler<KeyType, RealType>::Config config;
    config.blockSize = s.blockSize;
    config.variant   = kernelVariant;

    thrust::fill(thrust::device_pointer_cast(d_valuesOut.data()),
                 thrust::device_pointer_cast(d_valuesOut.data() + s.nnz),
                 RealType(0));

    CvfemHexAssembler<KeyType, RealType>::assembleFull(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), s.elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        d_gamma.data(), d_phi.data(), d_beta.data(),
        d_gphi_x.data(), d_gphi_y.data(), d_gphi_z.data(),
        d_mdot_zero.data(), s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
        s.d_node_to_dof.data(), d_nodeOwnership.data(), d_matrix, d_rhs_unused.data(), config);
    cudaDeviceSynchronize();
    cudaFree(d_matrix);
}

// Wrap a CSR (rowPtr/colInd/values) into a SparseMatrix and copy into A.
// numCols = numTotalDofs so ghost columns are addressable; numRows = owned only.
template<typename RealType>
void wrapIntoSparseMatrix(SparseMatrix<int, RealType, cstone::GpuTag>& A,
                          int numOwnedDofs, int numTotalDofs, int nnz,
                          const int* d_rowPtr, const int* d_colInd, const RealType* d_values)
{
    A.allocate(numOwnedDofs, numTotalDofs, nnz);
    cudaMemcpy(A.rowOffsetsPtr(), d_rowPtr, (numOwnedDofs + 1) * sizeof(int), cudaMemcpyDeviceToDevice);
    cudaMemcpy(A.colIndicesPtr(), d_colInd, nnz * sizeof(int),                cudaMemcpyDeviceToDevice);
    cudaMemcpy(A.valuesPtr(),     d_values, nnz * sizeof(RealType),           cudaMemcpyDeviceToDevice);
}

// =============================================================================
// Bochev-Dohrmann polynomial pressure projection stabilization for Q1 hexes.
// Adds, per element, the symmetric 8x8 matrix
//   S^e_ij = tau * V_e * (delta_ij/8 - 1/64)
// to the (already-assembled) pressure matrix. This is the lumped-mass form
// of the BD stabilizer with Pi(q) = element-mean. The matrix is symmetric,
// has row sum zero (consistent: kills constant-pressure mode), positive
// semidefinite, and adds element-local control of the discrete pressure
// null-space that equal-order CG-FEM cannot otherwise see.
//
// tau has units of [pressure] / [pressure-Poisson]. For (M/dt + nu*K) the
// natural choice is tau = h^2 / (4 * nu * rho); the agent docs use this
// scaling. We pass `tau` already-computed by the caller.
//
// The scatter pattern mirrors CvfemHexAssembler::assembleFull: each (i,j)
// pair of element corners writes to either the diagonal slot (if i and j
// reduce to the same DOF after periodic collapse) or finds the column
// position by binary search and atomicAdds there.
// =============================================================================

template<typename KeyType, typename RealType>
__global__ void assembleBDStabilizationKernel(const KeyType* c0, const KeyType* c1,
                                              const KeyType* c2, const KeyType* c3,
                                              const KeyType* c4, const KeyType* c5,
                                              const KeyType* c6, const KeyType* c7,
                                              const RealType* nodeX,
                                              const RealType* nodeY,
                                              const RealType* nodeZ,
                                              size_t numElements,
                                              const int* nodeToDof,
                                              const uint8_t* ownership,
                                              RealType tau,
                                              const int* rowPtr,
                                              const int* colInd,
                                              const int* diagPtr,
                                              RealType* values)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};

    // Element volume via box approximation: V_e ~ |x7-x0| * |y7-y0| * |z7-z0|.
    // For a uniform hex mesh this is exact; for distorted hexes it's the
    // determinant of the trilinear Jacobian at the element center, an O(1)
    // approximation that suffices for a stabilization parameter.
    RealType V_e;
    {
        RealType x0 = nodeX[n[0]], x6 = nodeX[n[6]];
        RealType y0 = nodeY[n[0]], y6 = nodeY[n[6]];
        RealType z0 = nodeZ[n[0]], z6 = nodeZ[n[6]];
        V_e = fabs((x6 - x0) * (y6 - y0) * (z6 - z0));
    }
    RealType diagCoef = tau * V_e * (RealType(1) / RealType(8) - RealType(1) / RealType(64));
    RealType offCoef  = tau * V_e * (-RealType(1) / RealType(64));

    int dofs[8];
    uint8_t own[8];
    #pragma unroll
    for (int i = 0; i < 8; ++i) { dofs[i] = nodeToDof[n[i]]; own[i] = ownership[n[i]]; }

    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        if (own[i] != 1) continue;
        int row_dof = dofs[i];
        if (row_dof < 0) continue;
        int row_start = rowPtr[row_dof];
        int row_end   = rowPtr[row_dof + 1];
        int diag_pos  = diagPtr[row_dof];

        #pragma unroll
        for (int j = 0; j < 8; ++j) {
            int col_dof = dofs[j];
            if (col_dof < 0) continue;
            RealType v = (i == j) ? diagCoef : offCoef;

            if (col_dof == row_dof) {
                atomicAdd(&values[diag_pos], v);
            } else {
                int lo = row_start, hi = row_end - 1;
                while (lo <= hi) {
                    int mid = (lo + hi) >> 1;
                    int col = colInd[mid];
                    if (col == col_dof) { atomicAdd(&values[mid], v); break; }
                    else if (col < col_dof) lo = mid + 1;
                    else                    hi = mid - 1;
                }
            }
        }
    }
}

// Helper: add BD stabilization to the pressure matrix s.d_valuesPre.
template<typename KeyType, typename RealType>
void addBochevDohrmannStab(NSStepper<KeyType, RealType>& s, RealType tau)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();

    size_t nElem = s.elementCount;
    int blockSize = s.blockSize;
    int grid = int((nElem + blockSize - 1) / blockSize);
    if (grid <= 0) return;

    assembleBDStabilizationKernel<KeyType, RealType><<<grid, blockSize>>>(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
        std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
        std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
        d_x.data(), d_y.data(), d_z.data(),
        nElem, s.d_node_to_dof.data(), d_nodeOwnership.data(),
        tau,
        s.d_rowPtr.data(), s.d_colInd.data(), s.d_diagPtr.data(),
        s.d_valuesPre.data());
    cudaDeviceSynchronize();
}

// =============================================================================
// Setup: builds DOF mapping, sparsity, two matrices, lumped mass, BC mask,
// pressure pin, area vectors. Runs ONCE before the time loop.
// =============================================================================

template<typename KeyType, typename RealType>
void setupNSStepper(NSStepper<KeyType, RealType>& s,
                    RealType nu,
                    RealType dt,
                    CvfemKernelVariant kernelVariant)
{
    // Live-print each setup lap on rank 0; the final report still summarizes.
    PhaseTimer pt(/*sync=*/true, /*liveRank=*/0);

    s.nodeCount    = s.domain.getNodeCount();
    s.elementCount = s.domain.getElementCount();

    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    s.domain.cacheNodeCoordinates();
    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();
    pt.lap("lazy domain prep");

    // DOF mapping: owned -> [0, numOwnedDofs); ghost -> [numOwnedDofs, nodeCount).
    s.d_node_to_dof.resize(s.nodeCount);
    s.numOwnedDofs = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(), s.d_node_to_dof.data(), s.nodeCount);
    s.numTotalDofs = static_cast<int>(s.nodeCount);

    // PERIODIC DOF COLLAPSE: slave nodes share their master's DOF index, and
    // we then COMPACT the DOF numbering so there are no orphan rows in the
    // sparsity. After this pass:
    //   - nodeToDof[slave] = nodeToDof[master]
    //   - all surviving DOF indices are dense in [0, nLive)
    //   - s.numOwnedDofs = nLive
    // The downstream sparsity builder / assembler / kernels read nodeToDof
    // unchanged, so they index through the compact DOF space transparently.
    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Periodic && s.periodicMap)
    {
        const int* d_partner = s.periodicMap->d_periodicPartner.data();
        int* d_n2d           = s.d_node_to_dof.data();
        size_t nNodes        = s.nodeCount;
        int numOwned         = s.numOwnedDofs;

        // Pass 1: redirect slave -> master DOF (using OLD master dof).
        cstone::DeviceVector<int> d_oldDof(s.d_node_to_dof);
        thrust::for_each(thrust::device,
                         thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(nNodes),
                         [d_partner, d_old = d_oldDof.data(), d_n2d]
                         __device__ (size_t i) {
                             int master = d_partner[i];
                             if (master >= 0)
                                 d_n2d[i] = d_old[master];
                         });
        cudaDeviceSynchronize();

        // Pass 2: mark live DOFs.
        cstone::DeviceVector<int> d_isLive(numOwned, 0);
        thrust::for_each(thrust::device,
                         thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(nNodes),
                         [d_n2d, d_live = d_isLive.data(), numOwned,
                          d_own = d_nodeOwnership.data()] __device__ (size_t i) {
                             if (d_own[i] != 1) return;
                             int dof = d_n2d[i];
                             if (dof >= 0 && dof < numOwned) d_live[dof] = 1;
                         });
        cudaDeviceSynchronize();

        // Pass 3: compact map: compact[old_dof] = prefix-sum of d_isLive.
        cstone::DeviceVector<int> d_compact(numOwned, -1);
        thrust::exclusive_scan(thrust::device, d_isLive.begin(), d_isLive.end(),
                               d_compact.begin());
        // For orphan DOFs we leave compact = scan value but they will be
        // unreachable from nodeToDof, so it doesn't matter. Mark them -1
        // explicitly for safety so any stale lookup is detected.
        thrust::transform(thrust::device,
                          d_isLive.begin(), d_isLive.end(),
                          d_compact.begin(), d_compact.begin(),
                          [] __device__ (int live, int idx) {
                              return live ? idx : -1;
                          });
        cudaDeviceSynchronize();

        int nLive = int(thrust::reduce(thrust::device,
                                        d_isLive.begin(), d_isLive.end(), 0));

        // Pass 4: rewrite nodeToDof through the compact map. Ghost nodes
        // (own != 1) keep their old DOF since they index into another rank's
        // (uncompacted) numbering -- but multi-rank periodic isn't supported
        // yet, so on >1 rank this would need cross-rank renumbering.
        thrust::for_each(thrust::device,
                         thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(nNodes),
                         [d_n2d, d_comp = d_compact.data(), numOwned,
                          d_own = d_nodeOwnership.data()] __device__ (size_t i) {
                             int dof = d_n2d[i];
                             if (dof >= 0 && dof < numOwned && d_own[i] == 1)
                                 d_n2d[i] = d_comp[dof];
                         });
        cudaDeviceSynchronize();

        s.numOwnedDofs = nLive;
        if (s.rank == 0)
        {
            int nOrphan = numOwned - nLive;
            std::cout << "  periodic DOF collapse: " << nOrphan << " slave DOFs collapsed -> "
                      << nLive << " active equations (was " << numOwned << ")\n";
        }
    }
    pt.lap("DOF mapping");

    // FULL 27-NNZ sparsity (same pattern used by both matrices). Building once
    // means the velocity and pressure matrices share row/col/diag pointers.
    s.d_rowPtr.resize(s.numTotalDofs + 1);
    s.d_diagPtr.resize(s.numTotalDofs);
    s.nnz = CvfemSparsityBuilder<KeyType>::buildFullSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), s.elementCount,
        s.d_node_to_dof.data(), s.numTotalDofs,
        s.d_rowPtr.data(), nullptr, nullptr, 0);
    s.d_colInd.resize(s.nnz);
    CvfemSparsityBuilder<KeyType>::buildFullSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), s.elementCount,
        s.d_node_to_dof.data(), s.numTotalDofs,
        s.d_rowPtr.data(), s.d_colInd.data(), s.d_diagPtr.data(), 0);
    s.d_valuesVel.resize(s.nnz);
    s.d_valuesPre.resize(s.nnz);
    pt.lap("sparsity build");

    // Area vectors (single source of truth for face geometry across operators).
    s.d_areaVec_x.resize(s.elementCount * 12);
    s.d_areaVec_y.resize(s.elementCount * 12);
    s.d_areaVec_z.resize(s.elementCount * 12);
    precomputeAreaVectorsGpu<KeyType, RealType>(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), s.elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data());
    pt.lap("area vectors");

    // Velocity matrix: nu * K assembled first (then we add M/dt and apply BC).
    assembleLaplacian<KeyType, RealType>(s, s.d_node_to_dof, s.d_valuesVel, kernelVariant);
    // The Laplacian above was assembled with gamma=1 -> raw K. Scale by nu for
    // the diffusion matrix.
    {
        const RealType nuScale = nu;
        RealType* valPtr = s.d_valuesVel.data();
        thrust::for_each(thrust::device,
                          thrust::counting_iterator<int>(0),
                          thrust::counting_iterator<int>(s.nnz),
                          [valPtr, nuScale] __device__ (int k) { valPtr[k] *= nuScale; });
        cudaDeviceSynchronize();
    }
    pt.lap("assembly nu*K (velocity)");

    // Pressure matrix: same K (no nu scale). Assemble into d_valuesPre.
    assembleLaplacian<KeyType, RealType>(s, s.d_node_to_dof, s.d_valuesPre, kernelVariant);
    pt.lap("assembly K (pressure)");

    // Lumped mass (per-node + reverse-halo). Identical to B.2/B.3/B.4 pattern.
    s.d_mass.resize(s.numOwnedDofs);
    thrust::fill(thrust::device_pointer_cast(s.d_mass.data()),
                 thrust::device_pointer_cast(s.d_mass.data() + s.numOwnedDofs),
                 RealType(0));
    // Persistent per-NODE lumped mass with halo-exchanged ghost values.
    // Needed by the assembled-DDT preconditioner builder so face entries
    // touching a ghost endpoint use the correct V_ghost (which lives on
    // another rank). Without this, every assembled row that has a ghost
    // neighbor is systematically wrong (entry magnitude ~halved), AMG sees
    // an unbalanced operator, and DDT+Hypre diverges on cube16 step 1.
    s.d_massNode.resize(s.nodeCount);
    thrust::fill(thrust::device_pointer_cast(s.d_massNode.data()),
                 thrust::device_pointer_cast(s.d_massNode.data() + s.nodeCount),
                 RealType(0));
    {
        size_t startElem = s.domain.startIndex();
        size_t numLocal  = s.domain.localElementCount();
        if (numLocal > 0)
        {
            int eBlocks = int((numLocal + s.blockSize - 1) / s.blockSize);
            computeLumpedMassPerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
                std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
                std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
                std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
                d_x.data(), d_y.data(), d_z.data(),
                s.d_massNode.data(), startElem, numLocal);
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(s.d_massNode);
        maybePeriodicSum<KeyType, RealType>(s, s.d_massNode);
        // Forward halo: ghost slots get owner-rank V values so DDT assembler
        // can read s.d_massNode[ghostNode] and get the correct V.
        s.domain.exchangeNodeHalo(s.d_massNode);
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        gatherOwnedNodeMassToDofKernel<RealType><<<nBlocks, s.blockSize>>>(
            s.d_massNode.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_mass.data(), s.nodeCount);
        cudaDeviceSynchronize();
    }
    pt.lap("lumped mass");

    // Assemble the matrix-free DDT operator into a CSR matrix with the
    // DEDICATED two-hop-through-shared-node sparsity. The operator
    // A = D M^-1 D^T couples pairs (r,c) of nodes that share ANY face-endpoint
    // i, INCLUDING the case where r and c live in different hex elements
    // that share only the corner i. buildFullSparsity only emits intra-element
    // pairs, so cross-element pairs that the assembler kernel writes get
    // silently dropped by atomicAddSparseEntry's binary search. The dedicated
    // buildDDTSparsity uses the same node->element CSR the assembler uses
    // and emits exactly the pairs the assembler writes; verified bit-exact
    // against applyDDTPerNode via MARS_DDT_PROBE_DIFF=1.
    const auto& d_n2eOff_sparsity  = s.domain.getNodeToElementOffsets();
    const auto& d_n2eList_sparsity = s.domain.getNodeToElementList();
    s.d_rowPtrDDT.resize(s.numTotalDofs + 1);
    s.d_diagPtrDDT.resize(s.numTotalDofs);
    s.nnzDDT = CvfemSparsityBuilder<KeyType>::buildDDTSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
        s.elementCount, s.nodeCount,
        d_n2eOff_sparsity.data(), d_n2eList_sparsity.data(),
        s.d_node_to_dof.data(), s.numTotalDofs,
        s.d_rowPtrDDT.data(), nullptr, nullptr, 0);
    s.d_colIndDDT.resize(s.nnzDDT);
    CvfemSparsityBuilder<KeyType>::buildDDTSparsity(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
        s.elementCount, s.nodeCount,
        d_n2eOff_sparsity.data(), d_n2eList_sparsity.data(),
        s.d_node_to_dof.data(), s.numTotalDofs,
        s.d_rowPtrDDT.data(), s.d_colIndDDT.data(), s.d_diagPtrDDT.data(), 0);
    s.d_valuesDDT.resize(s.nnzDDT);
    thrust::fill(thrust::device_pointer_cast(s.d_valuesDDT.data()),
                 thrust::device_pointer_cast(s.d_valuesDDT.data() + s.nnzDDT),
                 RealType(0));
    {
        // Node-driven assembly: one thread per OWNED node i, gather all
        // SCS faces incident to i via the node->element CSR (which spans
        // local + halo elements), enumerate all face pairs (g, f), and
        // contribute analytically. Captures cross-element face pairs that
        // a per-element kernel would miss because of the per-node M^-1
        // step in applyDDTPerNode (root cause documented above the kernel).
        const auto& d_n2eOff  = s.domain.getNodeToElementOffsets();
        const auto& d_n2eList = s.domain.getNodeToElementList();
        if (s.nodeCount > 0)
        {
            int nBlocks = int((s.nodeCount + s.blockSize - 1) / s.blockSize);
            assembleDDTPerNodeKernel<KeyType, RealType><<<nBlocks, s.blockSize>>>(
                std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
                std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
                std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
                std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                d_n2eOff.data(), d_n2eList.data(),
                s.d_node_to_dof.data(), d_nodeOwnership.data(),
                s.d_massNode.data(),
                s.d_rowPtrDDT.data(), s.d_colIndDDT.data(), s.numOwnedDofs,
                s.d_valuesDDT.data(), s.nodeCount);
            cudaDeviceSynchronize();
        }
        // (DDT-assembly NaN scan diagnostic removed; assembler verified clean.)
    }
    pt.lap("assembly D M^-1 D^T (node-driven, 27-NNZ)");

    // Optional diagonal shift on the DDT matrix: A := A + eps * I. Hypre PCG
    // strictly requires SPD. D M^-1 D^T is SPSD (constant mode in nullspace);
    // even with the pin BC, near-null modes survive in finite precision and
    // trigger PCG's pAp <= 0 check, returning error 1. A small diagonal shift
    // breaks the null mode (operator becomes A + eps*I with smallest eigenvalue
    // >= eps > 0) at the cost of perturbing the projection slightly. Default
    // 0 (no shift); set env MARS_DDT_DIAG_SHIFT=1e-10 to enable.
    {
        const char* shiftEnv = std::getenv("MARS_DDT_DIAG_SHIFT");
        RealType eps = (shiftEnv) ? RealType(std::atof(shiftEnv)) : RealType(0);
        if (eps > RealType(0))
        {
            int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
            const int* dpDdt = s.d_diagPtrDDT.data();
            RealType* vDdt   = s.d_valuesDDT.data();
            int n            = s.numOwnedDofs;
            thrust::for_each(thrust::device,
                thrust::counting_iterator<int>(0),
                thrust::counting_iterator<int>(n),
                [dpDdt, vDdt, eps] __device__ (int r) {
                    vDdt[dpDdt[r]] += eps;
                });
            cudaDeviceSynchronize();
            if (s.rank == 0)
                std::cout << "  [DDT] diagonal shift A := A + "
                          << std::scientific << eps << " * I (Hypre SPD regularization)\n"
                          << std::defaultfloat;
        }
    }

    // Add M/dt to the velocity matrix diagonal -> (M/dt + nu K).
    {
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        addLumpedMassDiagonalKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_mass.data(), s.d_diagPtr.data(), RealType(1) / dt,
            s.d_valuesVel.data(), s.numOwnedDofs);
        cudaDeviceSynchronize();
    }
    pt.lap("add M/dt (velocity)");

    // BDF2 velocity matrix: same K but with mass coefficient 3/(2 dt) instead
    // of 1/dt. Built as a copy of d_valuesVel BEFORE M/dt was added, then
    // M-diagonal increment of 3/(2 dt). Equivalently: start from current
    // d_valuesVel (which has M/dt added) and add an EXTRA 1/(2 dt) of mass
    // to upgrade the diagonal from M/dt to 3M/(2 dt).
    if (s.useBdf2)
    {
        s.d_valuesVel_bdf2.resize(s.d_valuesVel.size());
        thrust::copy(thrust::device,
                     s.d_valuesVel.begin(), s.d_valuesVel.end(),
                     s.d_valuesVel_bdf2.begin());
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        addLumpedMassDiagonalKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_mass.data(), s.d_diagPtr.data(), RealType(1) / (RealType(2) * dt),
            s.d_valuesVel_bdf2.data(), s.numOwnedDofs);
        cudaDeviceSynchronize();
        pt.lap("BDF2 velocity matrix (3M/(2dt) + nu K)");
    }

    // Global bounding box for BC marking. Reduce locally then MPI_Allreduce.
    auto xb = thrust::device_pointer_cast(d_x.data());
    auto yb = thrust::device_pointer_cast(d_y.data());
    auto zb = thrust::device_pointer_cast(d_z.data());
    RealType lxmin = thrust::reduce(thrust::device, xb, xb + s.nodeCount, RealType(1e30),  thrust::minimum<RealType>());
    RealType lxmax = thrust::reduce(thrust::device, xb, xb + s.nodeCount, RealType(-1e30), thrust::maximum<RealType>());
    RealType lymin = thrust::reduce(thrust::device, yb, yb + s.nodeCount, RealType(1e30),  thrust::minimum<RealType>());
    RealType lymax = thrust::reduce(thrust::device, yb, yb + s.nodeCount, RealType(-1e30), thrust::maximum<RealType>());
    RealType lzmin = thrust::reduce(thrust::device, zb, zb + s.nodeCount, RealType(1e30),  thrust::minimum<RealType>());
    RealType lzmax = thrust::reduce(thrust::device, zb, zb + s.nodeCount, RealType(-1e30), thrust::maximum<RealType>());

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&lxmin, &s.xmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lxmax, &s.xmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lymin, &s.ymin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lymax, &s.ymax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&lzmin, &s.zmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&lzmax, &s.zmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    s.bboxEps = RealType(1e-10) * std::max({s.xmax - s.xmin, s.ymax - s.ymin, s.zmax - s.zmin});
    pt.lap("global bbox");

    // BC mask + cavity target velocity.
    s.d_isBdryDof.resize(s.numOwnedDofs);
    thrust::fill(thrust::device_pointer_cast(s.d_isBdryDof.data()),
                 thrust::device_pointer_cast(s.d_isBdryDof.data() + s.numOwnedDofs),
                 uint8_t(0));
    s.d_uTarget.resize(s.nodeCount);
    s.d_vTarget.resize(s.nodeCount);
    s.d_wTarget.resize(s.nodeCount);
    {
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        using BCK = typename NSStepper<KeyType, RealType>::BCKind;
        if (s.bcKind == BCK::Cavity)
        {
            markCavityBCKernel<RealType><<<nBlocks, s.blockSize>>>(
                d_x.data(), d_y.data(), d_z.data(),
                d_nodeOwnership.data(), s.d_node_to_dof.data(),
                s.d_isBdryDof.data(),
                s.d_uTarget.data(), s.d_vTarget.data(), s.d_wTarget.data(),
                s.nodeCount, s.xmin, s.xmax, s.ymin, s.ymax, s.zmin, s.zmax,
                s.bboxEps, s.lidU);
        }
        else if (s.bcKind == BCK::Channel)
        {
            markChannelBCKernel<RealType><<<nBlocks, s.blockSize>>>(
                d_x.data(), d_y.data(), d_z.data(),
                d_nodeOwnership.data(), s.d_node_to_dof.data(),
                s.d_isBdryDof.data(),
                s.d_uTarget.data(), s.d_vTarget.data(), s.d_wTarget.data(),
                s.nodeCount, s.xmin, s.xmax, s.ymin, s.ymax, s.zmin, s.zmax,
                s.bboxEps, s.Uinf);
        }
        else if (s.bcKind == BCK::Wing)
        {
            // Wing BC: per-side-set Dirichlet driven by the local-node lists
            // populated by the driver. Build host arrays then copy to device:
            //   - isBdryDof[dof] = 1 for nodes on blade/inlet/farfield.
            //   - uTarget/vTarget/wTarget per side-set.
            // Pressure BC (outlet) is handled in the pressure-BC block below.
            std::vector<uint8_t> hostIsBdry(s.numOwnedDofs, 0);
            std::vector<RealType> hostUTgt(s.nodeCount, RealType(0));
            std::vector<RealType> hostVTgt(s.nodeCount, RealType(0));
            std::vector<RealType> hostWTgt(s.nodeCount, RealType(0));
            // Mirror node_to_dof + ownership to host so we can resolve DOF IDs.
            std::vector<int> hostNodeToDof(s.nodeCount, -1);
            std::vector<uint8_t> hostOwn(s.nodeCount, 0);
            thrust::copy(thrust::device_pointer_cast(s.d_node_to_dof.data()),
                         thrust::device_pointer_cast(s.d_node_to_dof.data() + s.nodeCount),
                         hostNodeToDof.begin());
            thrust::copy(thrust::device_pointer_cast(d_nodeOwnership.data()),
                         thrust::device_pointer_cast(d_nodeOwnership.data() + s.nodeCount),
                         hostOwn.begin());

            auto tag = [&] (const std::vector<int>& nodes, RealType u, RealType v, RealType w) {
                for (int li : nodes)
                {
                    if (li < 0 || (size_t)li >= s.nodeCount) continue;
                    // Target velocity is set on EVERY local copy (owned or ghost)
                    // so predictor/corrector see consistent values across rank
                    // boundaries.
                    hostUTgt[li] = u;
                    hostVTgt[li] = v;
                    hostWTgt[li] = w;
                    // Dirichlet mask is per OWNED DOF only.
                    if (hostOwn[li] != 1) continue;
                    int dof = hostNodeToDof[li];
                    if (dof < 0 || dof >= s.numOwnedDofs) continue;
                    hostIsBdry[dof] = 1;
                }
            };
            tag(s.wingBladeNodes,    RealType(0),     RealType(0), RealType(0));
            tag(s.wingInletNodes,    RealType(s.Uinf), RealType(0), RealType(0));
            tag(s.wingFarFieldNodes, RealType(s.Uinf), RealType(0), RealType(0));

            thrust::copy(hostIsBdry.begin(), hostIsBdry.end(),
                         thrust::device_pointer_cast(s.d_isBdryDof.data()));
            thrust::copy(hostUTgt.begin(), hostUTgt.end(),
                         thrust::device_pointer_cast(s.d_uTarget.data()));
            thrust::copy(hostVTgt.begin(), hostVTgt.end(),
                         thrust::device_pointer_cast(s.d_vTarget.data()));
            thrust::copy(hostWTgt.begin(), hostWTgt.end(),
                         thrust::device_pointer_cast(s.d_wTarget.data()));
            if (s.rank == 0)
            {
                std::cout << "  wing BC: " << s.wingBladeNodes.size()    << " blade nodes (no-slip), "
                          << s.wingInletNodes.size()    << " inlet nodes (u=" << s.Uinf << "), "
                          << s.wingFarFieldNodes.size() << " far-field nodes (u=" << s.Uinf << ")\n";
            }
        }
        // Periodic: no boundary DOFs (mask stays zero) and no per-node targets
        // (targets stay zero; predictor/corrector write the field on every node).
        cudaDeviceSynchronize();
    }
    // Sync ghost slots of the target arrays so per-node corrector/predictor
    // can read them on owned-boundary nodes even when those nodes are touched
    // by ghost elements on other ranks.
    s.domain.exchangeNodeHalo(s.d_uTarget);
    s.domain.exchangeNodeHalo(s.d_vTarget);
    s.domain.exchangeNodeHalo(s.d_wTarget);
    pt.lap("BC mark + target velocity");

    // Enforce velocity Dirichlet on the velocity matrix (row=0, diag=1).
    // Periodic: skip -- no DOFs are flagged so the kernel would be a no-op,
    // but we elide it explicitly to make the no-Dirichlet semantics obvious.
    if (s.bcKind != NSStepper<KeyType, RealType>::BCKind::Periodic)
    {
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        enforceBcMatrixKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_isBdryDof.data(), s.d_rowPtr.data(), s.d_colInd.data(),
            s.d_diagPtr.data(), s.d_valuesVel.data(), s.numOwnedDofs);
        if (s.useBdf2 && s.d_valuesVel_bdf2.size() > 0)
        {
            enforceBcMatrixKernel<RealType><<<dofBlocks, s.blockSize>>>(
                s.d_isBdryDof.data(), s.d_rowPtr.data(), s.d_colInd.data(),
                s.d_diagPtr.data(), s.d_valuesVel_bdf2.data(), s.numOwnedDofs);
        }
        cudaDeviceSynchronize();
    }
    pt.lap("velocity matrix BC");

    // Pressure null-space removal.
    //   Cavity: pin a single owned DOF closest to (xmin, ymin, zmin) and
    //   leave Neumann everywhere else. Works because pressure is free to be
    //   zero anywhere with pure-Dirichlet velocity BCs.
    //   Channel: pin the ENTIRE outflow face (x = xmax) to p = 0 by adding a
    //   per-owned-DOF mask, then turning matrix rows for those DOFs into
    //   identity. This matches the physically correct outflow condition for
    //   channel flow and avoids the single-corner-pin pathology that fought
    //   the natural inflow-to-outflow pressure gradient.
    using BCK = typename NSStepper<KeyType, RealType>::BCKind;
    if (s.bcKind == BCK::Periodic)
    {
        // Periodic: pressure has a 1D constant-mode null space that the solver
        // tolerates and the per-step removeMean call cleans up. No pin, no
        // outflow Dirichlet, no matrix-row enforcement on the pressure matrix.
        s.pressurePinDof  = -1;
        s.pressurePinRank = -1;
        if (s.rank == 0)
            std::cout << "  pressure BC: periodic (no pin; null-space removed each step via removeMean)\n";
    }
    else if (s.bcKind == BCK::Channel)
    {
        s.d_isPressureBdryDof.resize(s.numOwnedDofs);
        thrust::fill(thrust::device_pointer_cast(s.d_isPressureBdryDof.data()),
                     thrust::device_pointer_cast(s.d_isPressureBdryDof.data() + s.numOwnedDofs),
                     uint8_t(0));
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        markChannelPressureBCKernel<RealType><<<nBlocks, s.blockSize>>>(
            d_x.data(),
            d_nodeOwnership.data(), s.d_node_to_dof.data(),
            s.d_isPressureBdryDof.data(),
            s.nodeCount, s.xmax, s.bboxEps);
        cudaDeviceSynchronize();

        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        enforceBcMatrixKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_isPressureBdryDof.data(), s.d_rowPtr.data(), s.d_colInd.data(),
            s.d_diagPtr.data(), s.d_valuesPre.data(), s.numOwnedDofs);
        cudaDeviceSynchronize();

        s.pressurePinDof  = -1;     // disable corner-pin path
        s.pressurePinRank = -1;
        if (s.rank == 0)
            std::cout << "  pressure BC: Dirichlet p=0 on outflow face x=" << s.xmax << "\n";
    }
    else if (s.bcKind == BCK::Wing)
    {
        // Wing: pressure Dirichlet p=0 on the 'out' side-set (outlet face).
        s.d_isPressureBdryDof.resize(s.numOwnedDofs);
        thrust::fill(thrust::device_pointer_cast(s.d_isPressureBdryDof.data()),
                     thrust::device_pointer_cast(s.d_isPressureBdryDof.data() + s.numOwnedDofs),
                     uint8_t(0));
        // Host-side scatter from the outlet local-node list. Owned DOFs only.
        std::vector<uint8_t> hostMask(s.numOwnedDofs, 0);
        std::vector<int> hostNodeToDof(s.nodeCount, -1);
        std::vector<uint8_t> hostOwn(s.nodeCount, 0);
        thrust::copy(thrust::device_pointer_cast(s.d_node_to_dof.data()),
                     thrust::device_pointer_cast(s.d_node_to_dof.data() + s.nodeCount),
                     hostNodeToDof.begin());
        thrust::copy(thrust::device_pointer_cast(d_nodeOwnership.data()),
                     thrust::device_pointer_cast(d_nodeOwnership.data() + s.nodeCount),
                     hostOwn.begin());
        size_t ownedOutletCount = 0;
        for (int li : s.wingOutletNodes)
        {
            if (li < 0 || (size_t)li >= s.nodeCount) continue;
            if (hostOwn[li] != 1) continue;
            int dof = hostNodeToDof[li];
            if (dof < 0 || dof >= s.numOwnedDofs) continue;
            hostMask[dof] = 1;
            ++ownedOutletCount;
        }
        thrust::copy(hostMask.begin(), hostMask.end(),
                     thrust::device_pointer_cast(s.d_isPressureBdryDof.data()));

        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        enforceBcMatrixKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_isPressureBdryDof.data(), s.d_rowPtr.data(), s.d_colInd.data(),
            s.d_diagPtr.data(), s.d_valuesPre.data(), s.numOwnedDofs);
        cudaDeviceSynchronize();

        s.pressurePinDof  = -1;
        s.pressurePinRank = -1;
        if (s.rank == 0)
            std::cout << "  pressure BC: Dirichlet p=0 on outlet side-set ("
                      << ownedOutletCount << " owned DOFs on rank 0)\n";
    }
    else  // Cavity
    {
        cstone::DeviceVector<RealType> d_d2(s.numOwnedDofs, RealType(1e30));
        cstone::DeviceVector<int>      d_dofId(s.numOwnedDofs, -1);
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        findPressurePinCandidateKernel<RealType><<<nBlocks, s.blockSize>>>(
            d_x.data(), d_y.data(), d_z.data(),
            d_nodeOwnership.data(), s.d_node_to_dof.data(),
            d_d2.data(), d_dofId.data(),
            s.nodeCount, s.numOwnedDofs,
            s.xmin, s.ymin, s.zmin);
        cudaDeviceSynchronize();

        int localBestDof = -1;
        RealType localBestD2 = RealType(1e30);
        if (s.numOwnedDofs > 0)
        {
            auto dp = thrust::device_pointer_cast(d_d2.data());
            auto minIt = thrust::min_element(thrust::device, dp, dp + s.numOwnedDofs);
            size_t bestIdx = static_cast<size_t>(minIt - dp);
            thrust::copy(thrust::device_pointer_cast(d_d2.data() + bestIdx),
                         thrust::device_pointer_cast(d_d2.data() + bestIdx + 1),
                         &localBestD2);
            thrust::copy(thrust::device_pointer_cast(d_dofId.data() + bestIdx),
                         thrust::device_pointer_cast(d_dofId.data() + bestIdx + 1),
                         &localBestDof);
        }

        struct { double d; int r; } in{static_cast<double>(localBestD2), s.rank}, out{};
        MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        s.pressurePinRank = out.r;
        s.pressurePinDof  = (s.rank == out.r) ? localBestDof : -1;

        if (s.pressurePinDof >= 0)
        {
            enforcePinRowMatrixKernel<RealType><<<1, 1>>>(
                s.pressurePinDof, s.d_rowPtr.data(), s.d_colInd.data(),
                s.d_diagPtr.data(), s.d_valuesPre.data());
            cudaDeviceSynchronize();
        }
        if (s.rank == 0)
        {
            std::cout << "  pressure pin: rank=" << s.pressurePinRank
                      << ", dof=" << ((s.rank == s.pressurePinRank) ? s.pressurePinDof : -1)
                      << " (anchor near corner (" << s.xmin << "," << s.ymin << "," << s.zmin << "))\n";
        }
    }
    pt.lap("pressure BC");

    // Apply the same pressure-Dirichlet identity-row enforcement to the
    // assembled-DDT matrix so its rows for pinned outlet DOFs (and the cavity
    // single-pin) match the rule applyDDTPerNode uses at solve time. Without
    // this, AMG would see inconsistent rows and the preconditioner would
    // misbehave on boundary DOFs.
    // Build the dof->node inverse map + per-NODE pressure-Dirichlet mask
    // (halo-exchanged), used by the symmetric column-zeroing kernel below so
    // EVERY rank zeros its ghost-copy pin/Dirichlet columns. The earlier
    // approach only zeroed columns on the owning rank, leaving non-owning
    // ranks with stale couplings to the pin/mask DOFs -> globally asymmetric
    // matrix -> BoomerAMG produced wildly wrong phi (|gphi|=32 vs expected 5
    // on cube16-cavity at step 1).
    s.d_dofToNode.resize(s.numTotalDofs);
    thrust::fill(thrust::device_pointer_cast(s.d_dofToNode.data()),
                 thrust::device_pointer_cast(s.d_dofToNode.data() + s.numTotalDofs),
                 -1);
    {
        int nodeBlocks = int((s.nodeCount + s.blockSize - 1) / s.blockSize);
        buildDofToNodeKernel<<<nodeBlocks, s.blockSize>>>(
            s.d_node_to_dof.data(), s.d_dofToNode.data(),
            s.numTotalDofs, s.nodeCount);
        cudaDeviceSynchronize();
    }
    s.d_isPressureBdryNode.resize(s.nodeCount);
    thrust::fill(thrust::device_pointer_cast(s.d_isPressureBdryNode.data()),
                 thrust::device_pointer_cast(s.d_isPressureBdryNode.data() + s.nodeCount),
                 uint8_t(0));
    {
        int nodeBlocks = int((s.nodeCount + s.blockSize - 1) / s.blockSize);
        const uint8_t* maskPtr = (s.d_isPressureBdryDof.size() > 0)
                                 ? s.d_isPressureBdryDof.data() : nullptr;
        int pinDofLocal = s.pressurePinDof;  // -1 if this rank doesn't own pin
        buildPressureBdryNodeMaskKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            maskPtr, pinDofLocal,
            s.d_isPressureBdryNode.data(), s.nodeCount, s.numOwnedDofs);
        cudaDeviceSynchronize();
    }
    // Forward halo: ghost slots get owner-rank Dirichlet flag. cstone's
    // exchangeNodeHalo is locked to the domain's RealType (double), so we
    // exchange a double proxy then cast back. Proxy values are 0.0/1.0.
    {
        cstone::DeviceVector<RealType> d_maskProxy(s.nodeCount, RealType(0));
        thrust::transform(thrust::device,
                          thrust::device_pointer_cast(s.d_isPressureBdryNode.data()),
                          thrust::device_pointer_cast(s.d_isPressureBdryNode.data() + s.nodeCount),
                          thrust::device_pointer_cast(d_maskProxy.data()),
                          [] __device__ (uint8_t v) -> RealType { return v ? RealType(1) : RealType(0); });
        s.domain.exchangeNodeHalo(d_maskProxy);
        thrust::transform(thrust::device,
                          thrust::device_pointer_cast(d_maskProxy.data()),
                          thrust::device_pointer_cast(d_maskProxy.data() + s.nodeCount),
                          thrust::device_pointer_cast(s.d_isPressureBdryNode.data()),
                          [] __device__ (RealType v) -> uint8_t { return (v > RealType(0.5)) ? uint8_t(1) : uint8_t(0); });
        cudaDeviceSynchronize();
    }

    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Channel
        || s.bcKind == NSStepper<KeyType, RealType>::BCKind::Wing)
    {
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        // Row enforcement: zero each Dirichlet row, set diagonal to 1.
        enforceBcMatrixKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_isPressureBdryDof.data(), s.d_rowPtrDDT.data(), s.d_colIndDDT.data(),
            s.d_diagPtrDDT.data(), s.d_valuesDDT.data(), s.numOwnedDofs);
        // Symmetric column enforcement using per-NODE mask -- works on EVERY
        // rank (including ghost-copies of Dirichlet nodes).
        enforceBcColMatrixKernelByNode<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_isPressureBdryNode.data(),
            s.d_node_to_dof.data(), s.d_dofToNode.data(), s.numTotalDofs,
            s.d_rowPtrDDT.data(), s.d_colIndDDT.data(),
            s.d_valuesDDT.data(), s.numOwnedDofs);
        cudaDeviceSynchronize();
    }
    else if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Cavity
             && s.pressurePinDof >= 0)
    {
        // Row enforcement on pin row (one DOF, owning rank only).
        enforcePinRowMatrixKernel<RealType><<<1, 1>>>(
            s.pressurePinDof, s.d_rowPtrDDT.data(), s.d_colIndDDT.data(),
            s.d_diagPtrDDT.data(), s.d_valuesDDT.data());
    }
    // Re-enable cavity-mode column-zero on the DDT matrix. Hypre PCG requires
    // a symmetric matrix; without the col-clear, A is asymmetric (row 0 has
    // A[0,0]=1 with all other A[0,j]=0 from row-clear, but column 0 still has
    // A[r,0] != 0 from neighbour rows). PCG on an asymmetric matrix produces
    // garbage including NaN by step 2. The K-path matrix happens to work
    // without col-clear because it has different scaling; the DDT matrix
    // does not. With pin value b[pin]=0, symmetric treatment is just
    // row-clear + col-clear (no RHS modification needed).
    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Cavity)
    {
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        enforceBcColMatrixKernelByNode<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_isPressureBdryNode.data(),
            s.d_node_to_dof.data(), s.d_dofToNode.data(), s.numTotalDofs,
            s.d_rowPtrDDT.data(), s.d_colIndDDT.data(),
            s.d_valuesDDT.data(), s.numOwnedDofs);
        cudaDeviceSynchronize();
    }

    // Bochev-Dohrmann polynomial pressure projection stabilization. Added to
    // the assembled K-path pressure matrix (d_valuesPre, full 27-NNZ
    // sparsity). The DDT pressure solve also needs BD; we apply it
    // matrix-free in applyDDTPerNode() to avoid the sparsity mismatch (DDT
    // uses graph sparsity with only 7 NNZ/row, which cannot hold all 28
    // distinct corner-pair contributions per Q1 hex element).
    //
    // Q1-Q1 equal-order CG-FEM is LBB-unstable without BD on periodic
    // domains: the pressure checkerboard mode is unconstrained, grad(p)
    // explodes via predictor feedback, simulation blows up around step 30
    // regardless of dt or Re. BD adds the element-local symmetric PSD term
    //   S^e_ij = tau * V_e * (delta_ij/8 - 1/64)
    // which has row sum zero (kills constant pressure) and removes the
    // checkerboard null-space. Applied AFTER pin/outflow enforcement so
    // identity rows are preserved (BD is only enabled when bcKind == Periodic,
    // which has no pin and no outflow Dirichlet, so this ordering is moot in
    // practice but the assert is cheap).
    if (s.stabBochevDohrmann)
    {
        RealType tau = s.stabPressureTau;
        if (tau <= 0)
        {
            RealType bx = std::max(RealType(1e-30), s.xmax - s.xmin);
            RealType by = std::max(RealType(1e-30), s.ymax - s.ymin);
            RealType bz = std::max(RealType(1e-30), s.zmax - s.zmin);
            RealType h = std::cbrt((bx * by * bz) /
                                   RealType(std::max<size_t>(1, s.elementCount)));
            // Lumped-mass BD over-stabilizes by ~27x on the checkerboard mode
            // vs the consistent form (verified numerically: scratch/verify_bd.py).
            // Scale tau by 1/8 so we land between the Frobenius-match (~1/5)
            // and the eigenvalue-match (~1/3 to 1/27). Conservative but
            // physically sensible: kills the spurious mode strongly without
            // over-damping resolved pressure structures.
            tau = h * h / RealType(8);
        }
        addBochevDohrmannStab<KeyType, RealType>(s, tau);
        if (s.rank == 0)
            std::cout << "  pressure stabilization: Bochev-Dohrmann polynomial projection (K-path matrix + DDT matrix-free), tau=" << tau << "\n";
        pt.lap("BD pressure stabilization");
    }

    // Wrap into SparseMatrix for the solvers. Avel and Apre share K's 27-NNZ
    // layout; AddT uses its own reduced 7-NNZ layout (rowPtrDDT/colIndDDT).
    // Avel_bdf2 (when present) uses the SAME row layout as Avel but with the
    // 3M/(2 dt) mass coefficient instead of M/dt; used from step 1 of BDF2.
    wrapIntoSparseMatrix<RealType>(s.Avel, s.numOwnedDofs, s.numTotalDofs, s.nnz,
                                   s.d_rowPtr.data(), s.d_colInd.data(), s.d_valuesVel.data());
    if (s.useBdf2 && s.d_valuesVel_bdf2.size() > 0)
    {
        wrapIntoSparseMatrix<RealType>(s.Avel_bdf2, s.numOwnedDofs, s.numTotalDofs, s.nnz,
                                       s.d_rowPtr.data(), s.d_colInd.data(), s.d_valuesVel_bdf2.data());
    }
    wrapIntoSparseMatrix<RealType>(s.Apre, s.numOwnedDofs, s.numTotalDofs, s.nnz,
                                   s.d_rowPtr.data(), s.d_colInd.data(), s.d_valuesPre.data());
    wrapIntoSparseMatrix<RealType>(s.AddT, s.numOwnedDofs, s.numTotalDofs, s.nnzDDT,
                                   s.d_rowPtrDDT.data(), s.d_colIndDDT.data(), s.d_valuesDDT.data());
    pt.lap("SparseMatrix wrap (vel+pre+ddt)");

    // Optional diagnostic: env MARS_DDT_PROBE_DIFF=1 enables a one-shot
    // probe-and-compare check at setup time. For a few owned DOFs k, set
    // phi = e_k, apply BOTH the matrix-free applyDDTPerNode AND a manual
    // CSR SpMV of the assembled s.AddT, gather both back to per-DOF, and
    // print row-by-row diffs for the worst entries. If the assembled matrix
    // is bit-correct, all diffs print as ~1e-15; if it's wrong, we see
    // exactly which rows and how much. Used to debug the node-driven
    // assembler vs matrix-free identity.
    if (s.bcKind != NSStepper<KeyType, RealType>::BCKind::Periodic
        && std::getenv("MARS_DDT_PROBE_DIFF"))
    {
        // Choose probe DOFs scattered across the owned range. Skip dof 0
        // because in cavity mode it's the pressure pin and BC enforcement
        // (row-zero + col-zero) is INTENTIONALLY applied to the assembled
        // matrix but NOT to the matrix-free apply -- comparing them at the
        // pin would always show a false positive.
        std::vector<int> probeDofs;
        if (s.numOwnedDofs > 100) {
            probeDofs = { 50, 100, 500, 1500 };
        } else {
            for (int k = 1; k < std::min(s.numOwnedDofs, 5); ++k) probeDofs.push_back(k);
        }
        if (s.rank == 0)
            std::cout << "  [MARS_DDT_PROBE_DIFF] running " << probeDofs.size()
                      << " probe(s) of e_k through both apply-paths\n";

        cstone::DeviceVector<RealType> phiNode(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> yMfNode(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> xDof(s.numTotalDofs, RealType(0));
        cstone::DeviceVector<RealType> yAsmDof(s.numOwnedDofs, RealType(0));

        const int*     n2dPtr = s.d_node_to_dof.data();
        const uint8_t* ownPtr = d_nodeOwnership.data();
        size_t nNodes = s.nodeCount;
        size_t nTotal = s.numTotalDofs;
        size_t nOwned = s.numOwnedDofs;

        for (int targetDof : probeDofs)
        {
            if (targetDof >= s.numOwnedDofs) continue;
            // Reset
            cudaMemset(phiNode.data(), 0, nNodes * sizeof(RealType));
            cudaMemset(yMfNode.data(), 0, nNodes * sizeof(RealType));
            cudaMemset(xDof.data(),    0, nTotal * sizeof(RealType));
            cudaMemset(yAsmDof.data(), 0, nOwned * sizeof(RealType));

            // Set phiNode[node such that nodeToDof[node]==targetDof, owned] = 1.
            thrust::for_each(thrust::device,
                thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(nNodes),
                [n2dPtr, ownPtr, targetDof, p = phiNode.data()]
                __device__ (size_t k) {
                    if (ownPtr[k] == 1 && n2dPtr[k] == targetDof) p[k] = RealType(1);
                });
            // Set xDof[targetDof] = 1.
            thrust::for_each(thrust::device,
                thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(1),
                [targetDof, x = xDof.data()] __device__ (int) { x[targetDof] = RealType(1); });
            cudaDeviceSynchronize();

            // Forward halo of phiNode so D^T scatter sees correct ghost values.
            s.domain.exchangeNodeHalo(phiNode);

            // Matrix-free path
            applyDDTPerNode<KeyType, RealType>(s, phiNode, yMfNode, gx, gy, gz);

            // Manual CSR SpMV: yAsm[r] = sum_{k in rowPtr[r]..rowPtr[r+1]} A[r,k] * x[k]
            const int* rpPtr = s.d_rowPtrDDT.data();
            const int* ciPtr = s.d_colIndDDT.data();
            const RealType* vPtr = s.d_valuesDDT.data();
            const RealType* xPtr = xDof.data();
            RealType* yAsmPtr = yAsmDof.data();
            thrust::for_each(thrust::device,
                thrust::counting_iterator<int>(0), thrust::counting_iterator<int>((int)nOwned),
                [rpPtr, ciPtr, vPtr, xPtr, yAsmPtr] __device__ (int r) {
                    int rs = rpPtr[r], re = rpPtr[r + 1];
                    RealType acc = RealType(0);
                    for (int j = rs; j < re; ++j) acc += vPtr[j] * xPtr[ciPtr[j]];
                    yAsmPtr[r] = acc;
                });
            cudaDeviceSynchronize();

            // Pull both vectors back to host (small probe count, OK).
            std::vector<RealType> hMf(nNodes), hAsm(nOwned);
            std::vector<int> hN2D(nNodes);
            std::vector<uint8_t> hOwn(nNodes);
            thrust::copy(thrust::device_pointer_cast(yMfNode.data()),
                         thrust::device_pointer_cast(yMfNode.data() + nNodes), hMf.begin());
            thrust::copy(thrust::device_pointer_cast(yAsmDof.data()),
                         thrust::device_pointer_cast(yAsmDof.data() + nOwned), hAsm.begin());
            thrust::copy(thrust::device_pointer_cast(s.d_node_to_dof.data()),
                         thrust::device_pointer_cast(s.d_node_to_dof.data() + nNodes), hN2D.begin());
            thrust::copy(thrust::device_pointer_cast(d_nodeOwnership.data()),
                         thrust::device_pointer_cast(d_nodeOwnership.data() + nNodes), hOwn.begin());

            // Compare per OWNED dof. Print top-5 |diff|.
            std::vector<std::pair<int, RealType>> diffs;
            // gather yMf[node] -> per-dof
            std::vector<RealType> yMfDof(nOwned, RealType(0));
            for (size_t n = 0; n < nNodes; ++n) {
                if (hOwn[n] == 1 && hN2D[n] >= 0 && hN2D[n] < (int)nOwned)
                    yMfDof[hN2D[n]] = hMf[n];
            }
            for (int r = 0; r < (int)nOwned; ++r) {
                RealType d = hAsm[r] - yMfDof[r];
                if (std::abs(d) > 1e-12)
                    diffs.emplace_back(r, d);
            }
            std::sort(diffs.begin(), diffs.end(),
                      [](auto& a, auto& b){ return std::abs(a.second) > std::abs(b.second); });
            if (s.rank == 0) {
                RealType maxDiff = diffs.empty() ? RealType(0) : std::abs(diffs[0].second);
                std::cout << "    probe dof " << targetDof << " : "
                          << diffs.size() << " rows with |diff|>1e-12, maxDiff="
                          << std::scientific << std::setprecision(3) << maxDiff
                          << std::defaultfloat << "\n";
                for (size_t k = 0; k < std::min(diffs.size(), size_t(5)); ++k) {
                    int r = diffs[k].first;
                    std::cout << "      row " << r
                              << "  asm=" << std::scientific << std::setprecision(6) << hAsm[r]
                              << "  mf="  << yMfDof[r]
                              << "  diff="  << diffs[k].second
                              << std::defaultfloat << "\n";
                }
            }
        }
    }

#ifdef MARS_ENABLE_HYPRE
    if (s.solverKind == SolverKind::Hypre)
    {
        // Build contiguous global-DOF numbering (same recipe as
        // mars_amr_advdiff/mars_amr_pressure_poisson). Used by both matrix
        // solves; pre and vel share the same map since they share row layout.
        int64_t numOwnedLocal = static_cast<int64_t>(s.numOwnedDofs);
        int64_t exscan = 0;
        MPI_Exscan(&numOwnedLocal, &exscan, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
        if (s.rank == 0) exscan = 0;
        s.globalRowStart    = exscan;
        s.globalRowEnd      = s.globalRowStart + numOwnedLocal;
        MPI_Allreduce(&numOwnedLocal, &s.numInteriorGlobal, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);

        cstone::DeviceVector<RealType> d_nodeGlobalDof(s.nodeCount, RealType(-1));
        {
            const int* nodeToDofPtr = s.d_node_to_dof.data();
            const uint8_t* ownPtr   = d_nodeOwnership.data();
            RealType rs             = static_cast<RealType>(s.globalRowStart);
            thrust::for_each(thrust::device,
                              thrust::counting_iterator<size_t>(0),
                              thrust::counting_iterator<size_t>(s.nodeCount),
                              [nodeToDofPtr, ownPtr, rs,
                               out = d_nodeGlobalDof.data()] __device__(size_t i)
                              {
                                  if (ownPtr[i] == 1)
                                  {
                                      int dof = nodeToDofPtr[i];
                                      if (dof >= 0) out[i] = rs + static_cast<RealType>(dof);
                                  }
                              });
            s.domain.exchangeNodeHalo(d_nodeGlobalDof);
        }
        // Build the local-DOF -> global-DOF map directly on the device. The
        // per-node d_nodeGlobalDof has already been halo-exchanged above, so
        // ghost slots have the owner's global ID.
        s.d_localToGlobalDof.resize(s.numTotalDofs);
        thrust::fill(s.d_localToGlobalDof.begin(),
                     s.d_localToGlobalDof.end(),
                     HYPRE_BigInt(-1));
        {
            const int* nodeToDofPtr  = s.d_node_to_dof.data();
            const RealType* gPtr     = d_nodeGlobalDof.data();
            HYPRE_BigInt* outPtr     = thrust::raw_pointer_cast(s.d_localToGlobalDof.data());
            thrust::for_each(thrust::device,
                              thrust::counting_iterator<size_t>(0),
                              thrust::counting_iterator<size_t>(s.nodeCount),
                              [nodeToDofPtr, gPtr, outPtr] __device__(size_t i)
                              {
                                  int dof = nodeToDofPtr[i];
                                  if (dof >= 0)
                                  {
                                      RealType g = gPtr[i];
                                      if (g >= RealType(0))
                                          outPtr[dof] = static_cast<HYPRE_BigInt>(g);
                                  }
                              });
        }
        cudaDeviceSynchronize();
        pt.lap("Hypre global DOF map (on-device)");
    }
#endif

    // Allocate the per-step solution / staging buffers (all per-node, sized
    // nodeCount so ghost slots are addressable).
    auto zfill = [&](cstone::DeviceVector<RealType>& v) {
        v.resize(s.nodeCount);
        thrust::fill(thrust::device_pointer_cast(v.data()),
                     thrust::device_pointer_cast(v.data() + s.nodeCount),
                     RealType(0));
    };
    zfill(s.d_u); zfill(s.d_v); zfill(s.d_w);
    zfill(s.d_uStar); zfill(s.d_vStar); zfill(s.d_wStar);
    zfill(s.d_uStarStar); zfill(s.d_vStarStar); zfill(s.d_wStarStar);
    zfill(s.d_p); zfill(s.d_phi);
    zfill(s.d_gradPx); zfill(s.d_gradPy); zfill(s.d_gradPz);
    zfill(s.d_gradPhix); zfill(s.d_gradPhiy); zfill(s.d_gradPhiz);
    // BDF2/EXT2 history: velocity at n-1 and the persistent advection
    // accumulators for N(u^n) and N(u^(n-1)). bdfStep reset to 0; the very
    // first runNsStep call uses BDF1 (no history), step 1+ uses BDF2.
    if (s.useBdf2)
    {
        zfill(s.d_u_nm1); zfill(s.d_v_nm1); zfill(s.d_w_nm1);
        zfill(s.d_advU_n);   zfill(s.d_advV_n);   zfill(s.d_advW_n);
        zfill(s.d_advU_nm1); zfill(s.d_advV_nm1); zfill(s.d_advW_nm1);
        s.bdfStep = 0;
    }
    pt.lap("field allocation");

    pt.report(s.rank, "setup");
}

// =============================================================================
// Set the IC -- velocity = 0 inside, lid value on top face, pressure = 0.
// Boundary node slots already carry the right values via d_uTarget etc; we
// just initialize the actual velocity arrays from those targets so the very
// first predictor sees the BC correctly. Halo exchange syncs ghost slots.
// =============================================================================

template<typename RealType>
__global__ void initialConditionKernel(const uint8_t* isBdryDof,
                                       const int* nodeToDof,
                                       const uint8_t* ownership,
                                       const RealType* uTarget,
                                       const RealType* vTarget,
                                       const RealType* wTarget,
                                       RealType* u,
                                       RealType* v,
                                       RealType* w,
                                       size_t numNodes,
                                       // If nonzero, set interior u=interiorU
                                       // instead of u=0 (wing free-stream IC).
                                       RealType interiorU,
                                       RealType interiorV,
                                       RealType interiorW)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    // Default interior = (interiorU, interiorV, interiorW). For cavity that is
    // (0,0,0); for wing flow it is the free-stream (Uinf, 0, 0) so step 1's
    // predictor doesn't see a 1-cell-thick velocity discontinuity at the inlet
    // and far-field Dirichlet faces.
    RealType uu = interiorU, vv = interiorV, ww = interiorW;
    if (ownership[i] == 1)
    {
        int dof = nodeToDof[i];
        if (dof >= 0 && isBdryDof[dof])
        {
            uu = uTarget[i]; vv = vTarget[i]; ww = wTarget[i];
        }
    }
    u[i] = uu; v[i] = vv; w[i] = ww;
}

template<typename KeyType, typename RealType>
void applyInitialCondition(NSStepper<KeyType, RealType>& s)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    // Interior IC: cavity/channel default to (0,0,0). Wing starts at free-stream
    // (Uinf, 0, 0) so step 1's predictor doesn't see a velocity discontinuity
    // at the inlet/far-field boundaries.
    RealType iu = 0, iv = 0, iw = 0;
    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Wing)
    {
        iu = s.Uinf;
    }
    initialConditionKernel<RealType><<<nBlocks, s.blockSize>>>(
        s.d_isBdryDof.data(), s.d_node_to_dof.data(), d_nodeOwnership.data(),
        s.d_uTarget.data(), s.d_vTarget.data(), s.d_wTarget.data(),
        s.d_u.data(), s.d_v.data(), s.d_w.data(),
        s.nodeCount, iu, iv, iw);
    thrust::fill(thrust::device_pointer_cast(s.d_p.data()),
                 thrust::device_pointer_cast(s.d_p.data() + s.nodeCount),
                 RealType(0));
    cudaDeviceSynchronize();
    s.domain.exchangeNodeHalo(s.d_u);
    s.domain.exchangeNodeHalo(s.d_v);
    s.domain.exchangeNodeHalo(s.d_w);
    s.domain.exchangeNodeHalo(s.d_p);
}

// =============================================================================
// One projection-method step. d_u/v/w and d_p come in as the n-state and leave
// as the (n+1)-state. All halo exchanges between substeps are explicit.
// =============================================================================

// Helper: per-component velocity solve. b is built by the caller into b_rhs;
// the solve writes the per-DOF solution into x, which we then scatter back to
// the per-node array qOut. Returns iteration count.
template<typename KeyType, typename RealType>
int solveOneComponent(NSStepper<KeyType, RealType>& s,
                      cstone::DeviceVector<RealType>& b_rhs,
                      cstone::DeviceVector<RealType>& xVec,
                      cstone::DeviceVector<RealType>& qOut,
                      typename NSStepper<KeyType, RealType>::Matrix& A,
                      KrylovHint krylov = KrylovHint::PCG)
{
    bool converged = false;
    int iters = 0;

    if (s.solverKind == SolverKind::CG)
    {
        ConjugateGradientSolver<RealType, int, cstone::GpuTag> solver(s.maxIter, s.tolerance);
        solver.setVerbose(false);
        solver.setOwnedSize(s.numOwnedDofs);
        if (s.numRanks > 1)
        {
            const int* dofMapPtr = s.d_node_to_dof.data();
            auto& dom = s.domain;
            solver.setHaloExchangeCallback(
                [&dom, dofMapPtr](cstone::DeviceVector<RealType>& p)
                {
                    dom.exchangeNodeHalo(p, dofMapPtr);
                });
        }
        converged = solver.solve(A, b_rhs, xVec);
        iters = solver.getIterations();
    }
    else
    {
#ifdef MARS_ENABLE_HYPRE
        // Env MARS_HYPRE_PRECOND=jacobi switches BoomerAMG -> Jacobi for debug.
        // Env MARS_HYPRE_KRYLOV=pcg|gmres (DEBUG ONLY) overrides the caller's
        // hint. Production behavior: PCG default; DDT pressure call passes
        // GMRES explicitly.
        KrylovHint krylovEff = krylov;
        {
            const char* kEnv = std::getenv("MARS_HYPRE_KRYLOV");
            if (kEnv)
            {
                std::string v(kEnv);
                if      (v == "pcg")   krylovEff = KrylovHint::PCG;
                else if (v == "gmres") krylovEff = KrylovHint::GMRES;
            }
        }
        if (krylovEff == KrylovHint::PCG)
        {
            using HSol = mars::fem::HyprePCGSolver<RealType, int, cstone::GpuTag>;
            auto precond = HSol::BOOMERAMG;
            {
                const char* p = std::getenv("MARS_HYPRE_PRECOND");
                if (p && std::string(p) == "jacobi") precond = HSol::JACOBI;
            }
            HSol hypreSolver(MPI_COMM_WORLD, s.maxIter, s.tolerance, precond);
            hypreSolver.setVerbose(false);
            converged = hypreSolver.solve(
                A, b_rhs, xVec,
                static_cast<int>(s.globalRowStart), static_cast<int>(s.globalRowEnd),
                0, static_cast<int>(s.numInteriorGlobal),
                s.d_localToGlobalDof);
            iters = converged ? hypreSolver.getLastIterations() : -2;
            if (s.rank == 0 && hypreSolver.getLastFinalResidual() > s.tolerance * 10)
            {
                std::cout << "  [hypre-pcg] iters=" << hypreSolver.getLastIterations()
                          << " final_res=" << hypreSolver.getLastFinalResidual()
                          << " (tol=" << s.tolerance << ", looser than expected)\n";
            }
        }
        else  // KrylovHint::GMRES
        {
            using HSol = mars::fem::HypreGMRESSolver<RealType, int, cstone::GpuTag>;
            auto precond = HSol::BOOMERAMG;
            {
                const char* p = std::getenv("MARS_HYPRE_PRECOND");
                if (p && std::string(p) == "jacobi") precond = HSol::JACOBI;
            }
            HSol hypreSolver(MPI_COMM_WORLD, s.maxIter, s.tolerance, precond);
            hypreSolver.setVerbose(false);
            converged = hypreSolver.solve(
                A, b_rhs, xVec,
                static_cast<int>(s.globalRowStart), static_cast<int>(s.globalRowEnd),
                0, static_cast<int>(s.numInteriorGlobal),
                s.d_localToGlobalDof);
            iters = converged ? hypreSolver.getLastIterations() : -2;
            if (s.rank == 0 && hypreSolver.getLastFinalResidual() > s.tolerance * 10)
            {
                std::cout << "  [hypre-gmres] iters=" << hypreSolver.getLastIterations()
                          << " final_res=" << hypreSolver.getLastFinalResidual()
                          << " (tol=" << s.tolerance << ", looser than expected)\n";
            }
        }
#else
        if (s.rank == 0)
            std::cerr << "ERROR: --solver=hypre requested but MARS built without MARS_ENABLE_HYPRE\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    cudaDeviceSynchronize();

    // x -> qOut (per-node scatter).
    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    scatterDofToNodeKernel<RealType><<<nBlocks, s.blockSize>>>(
        xVec.data(), s.d_node_to_dof.data(), qOut.data(), s.nodeCount);
    cudaDeviceSynchronize();

    return converged ? iters : -2;
}

// =============================================================================
// Phase E: matrix-free pressure-Poisson SpMV  A phi = (D M^{-1} D^T) phi.
//
// D is the per-node divergence scatter, D^T its literal transpose. A is SPSD:
//   phi^T A phi = (D^T phi)^T M^{-1} (D^T phi) >= 0.
// The transpose scatter (0.5*(p[L]-p[R])*A_f same-sign to L and R) integrates
// to D^T p ~ -V grad(p), so M^{-1} D^T ~ -grad. The corrector
//   u^{n+1} = u** - (dt/rho) * (M^{-1} D^T phi)
// therefore adds + (dt/rho) grad(phi) to u**, and closing  D u^{n+1} = 0  gives
//   A phi = -(rho/dt) D u**   (RHS sign matches the K-path's buildPressureRhs).
// The K-path lacks this algebraic identity: K (Galerkin Laplacian from
// assembleFull) is NOT D M^{-1} D^T, so its phi only approximately cancels
// D u**. The DDT path cancels exactly, hence div_max ~ roundoff every step.
// =============================================================================

// out = A phi = (D M^{-1} D^T) phi. gxAcc/gyAcc/gzAcc are caller-allocated
// scratch reused across CG iterations (no allocation in the hot path).
// Mirrors mars_amr_ddt.cu's applyDDT (Stages 1+2 validated). Local addition:
// identity-row enforcement (out[i]=phi[i]) on Dirichlet pressure DOFs --
// either a single corner pin (cavity) or the d_isPressureBdryDof mask over
// the outflow face (channel). Combined with b[i]=0 at the same slots, phi is
// pinned to 0 there exactly.
template<typename KeyType, typename RealType>
void applyDDTPerNode(NSStepper<KeyType, RealType>& s,
                     cstone::DeviceVector<RealType>& phi,
                     cstone::DeviceVector<RealType>& outAcc,
                     cstone::DeviceVector<RealType>& gxAcc,
                     cstone::DeviceVector<RealType>& gyAcc,
                     cstone::DeviceVector<RealType>& gzAcc)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    auto c0 = std::get<0>(d_conn).data(); auto c1 = std::get<1>(d_conn).data();
    auto c2 = std::get<2>(d_conn).data(); auto c3 = std::get<3>(d_conn).data();
    auto c4 = std::get<4>(d_conn).data(); auto c5 = std::get<5>(d_conn).data();
    auto c6 = std::get<6>(d_conn).data(); auto c7 = std::get<7>(d_conn).data();
    const size_t startElem = s.domain.startIndex();
    const size_t numLocal  = s.domain.localElementCount();
    const int eBlocks    = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
    const int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;

    auto zeroVec = [&](cstone::DeviceVector<RealType>& v) {
        thrust::fill(thrust::device_pointer_cast(v.data()),
                     thrust::device_pointer_cast(v.data() + s.nodeCount), RealType(0));
    };

    // Step a: g = D^T phi (un-normalized per-node 3-vector accumulator).
    zeroVec(gxAcc); zeroVec(gyAcc); zeroVec(gzAcc);
    if (eBlocks > 0)
    {
        applyDivTransposePerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            c0, c1, c2, c3, c4, c5, c6, c7, phi.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            gxAcc.data(), gyAcc.data(), gzAcc.data(), startElem, numLocal);
        cudaDeviceSynchronize();
    }
    // Reverse-halo each component so the owner sees contributions from neighbor
    // ranks' local elements at corner-only-neighbor nodes.
    s.domain.reverseExchangeNodeHaloAdd(gxAcc);
    s.domain.reverseExchangeNodeHaloAdd(gyAcc);
    s.domain.reverseExchangeNodeHaloAdd(gzAcc);
    maybePeriodicSum<KeyType, RealType>(s, gxAcc);
    maybePeriodicSum<KeyType, RealType>(s, gyAcc);
    maybePeriodicSum<KeyType, RealType>(s, gzAcc);

    // Step b: g <- M^{-1} g (in-place; gather is same-index, safe).
    normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        gxAcc.data(), gyAcc.data(), gzAcc.data(), s.d_mass.data(),
        s.d_node_to_dof.data(), d_nodeOwnership.data(),
        gxAcc.data(), gyAcc.data(), gzAcc.data(), s.nodeCount);
    cudaDeviceSynchronize();
    // Forward halo: D-scatter in step c reads g[iL], g[iR] at ghost slots too.
    s.domain.exchangeNodeHalo(gxAcc);
    s.domain.exchangeNodeHalo(gyAcc);
    s.domain.exchangeNodeHalo(gzAcc);

    // Step c: out = D g (un-normalized per-node divergence accumulator).
    zeroVec(outAcc);
    if (eBlocks > 0)
    {
        computeDivergencePerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            c0, c1, c2, c3, c4, c5, c6, c7,
            gxAcc.data(), gyAcc.data(), gzAcc.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            outAcc.data(), startElem, numLocal);
        cudaDeviceSynchronize();
    }
    s.domain.reverseExchangeNodeHaloAdd(outAcc);
    maybePeriodicSum<KeyType, RealType>(s, outAcc);

    // NOTE: BD is intentionally NOT applied to the matrix-free DDT operator.
    // Adding S*phi to (D M^-1 D^T) phi breaks the algebraic identity
    // div(u^{n+1}) = 0 that DDT is built for: the corrector
    // u^{n+1} = u** - dt/rho * grad phi cancels exactly only when phi solves
    // (D M^-1 D^T) phi = source. BD on the DDT path produces a phi that does
    // NOT cancel divergence, and the residual is proportional to the BD
    // strength -- destroying the whole point of using DDT. For periodic
    // pressure runs use --pressure-solve=K, which is regularized by BD via
    // the assembled d_valuesPre matrix.
    if (false && s.stabBochevDohrmann && eBlocks > 0)
    {
        RealType tau = s.stabPressureTau;
        if (tau <= 0)
        {
            RealType bx = std::max(RealType(1e-30), s.xmax - s.xmin);
            RealType by = std::max(RealType(1e-30), s.ymax - s.ymin);
            RealType bz = std::max(RealType(1e-30), s.zmax - s.zmin);
            RealType h = std::cbrt((bx * by * bz) /
                                   RealType(std::max<size_t>(1, s.elementCount)));
            // Match the K-path auto-tau (lumped BD over-stabilizes by ~27x
            // on checkerboard vs consistent form, see scratch/verify_bd.py).
            tau = h * h / RealType(8);
        }
        const RealType* d_xNode = s.domain.getNodeX().data();
        const RealType* d_yNode = s.domain.getNodeY().data();
        const RealType* d_zNode = s.domain.getNodeZ().data();
        RealType* outPtr        = outAcc.data();
        const RealType* phiPtr  = phi.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(numLocal),
            [c0, c1, c2, c3, c4, c5, c6, c7, startElem,
             d_xNode, d_yNode, d_zNode,
             tau, outPtr, phiPtr,
             ownPtr = d_nodeOwnership.data(), n2d = s.d_node_to_dof.data()]
            __device__ (size_t k) {
                size_t e = startElem + k;
                KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};
                RealType x0 = d_xNode[n[0]], x6 = d_xNode[n[6]];
                RealType y0 = d_yNode[n[0]], y6 = d_yNode[n[6]];
                RealType z0 = d_zNode[n[0]], z6 = d_zNode[n[6]];
                RealType V_e = fabs((x6 - x0) * (y6 - y0) * (z6 - z0));
                // S^e_ij = tau V_e (delta_ij/8 - 1/64); algebra collapses
                // (S^e phi)_i to tau V_e (phi_i - meanPhi) / 8.
                RealType phiE[8];
                #pragma unroll
                for (int i = 0; i < 8; ++i) phiE[i] = phiPtr[n[i]];
                RealType meanPhi = RealType(0);
                #pragma unroll
                for (int i = 0; i < 8; ++i) meanPhi += phiE[i];
                meanPhi *= RealType(1)/RealType(8);
                // (S phi)_i = tau V_e (phi_i - meanPhi) / 8  -- equivalent
                // to diagC phi_i + offC sum_{j!=i} phi_j after algebra.
                #pragma unroll
                for (int i = 0; i < 8; ++i) {
                    if (ownPtr[n[i]] != 1) continue;
                    if (n2d[n[i]] < 0) continue;
                    RealType contrib = tau * V_e * (phiE[i] - meanPhi) / RealType(8);
                    atomicAdd(&outPtr[n[i]], contrib);
                }
            });
        cudaDeviceSynchronize();
        s.domain.reverseExchangeNodeHaloAdd(outAcc);
        maybePeriodicSum<KeyType, RealType>(s, outAcc);
    }

    // Identity-row enforcement on pinned pressure DOFs (cavity: single corner;
    // channel: outflow-face mask). out[i] = phi[i] makes the row act as identity.
    const int pinDof          = s.pressurePinDof;
    const uint8_t* maskPtr    = (s.d_isPressureBdryDof.size() > 0)
                                ? s.d_isPressureBdryDof.data() : nullptr;
    if (pinDof >= 0 || maskPtr != nullptr)
    {
        const int* dofPtr     = s.d_node_to_dof.data();
        const uint8_t* ownPtr = d_nodeOwnership.data();
        const RealType* phPtr = phi.data();
        RealType* outPtr      = outAcc.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [pinDof, maskPtr, dofPtr, ownPtr, phPtr, outPtr] __device__ (size_t i) {
                if (ownPtr[i] != 1) return;
                int dof = dofPtr[i];
                if (dof < 0) return;
                bool flagged = (pinDof >= 0 && dof == pinDof) ||
                               (maskPtr != nullptr && maskPtr[dof] != 0);
                if (flagged) outPtr[i] = phPtr[i];
            });
        cudaDeviceSynchronize();
    }
}

// Plain Euclidean dot over owned nodes (CG operates on un-normalized acc form).
template<typename RealType>
RealType ownedDot(const cstone::DeviceVector<RealType>& a,
                  const cstone::DeviceVector<RealType>& b,
                  const cstone::DeviceVector<int>& d_nodeToDof,
                  const uint8_t* d_ownership,
                  size_t numNodes)
{
    const RealType* aPtr = a.data();
    const RealType* bPtr = b.data();
    const int* dofPtr    = d_nodeToDof.data();
    RealType localSum = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(numNodes),
        [aPtr, bPtr, dofPtr, d_ownership] __device__ (size_t i) -> RealType {
            if (d_ownership[i] != 1 || dofPtr[i] < 0) return RealType(0);
            return aPtr[i] * bPtr[i];
        }, RealType(0), thrust::plus<RealType>());
    RealType globalSum = 0;
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&localSum, &globalSum, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    return globalSum;
}

// out = a + alpha * b on owned nodes only (ghosts untouched).
template<typename RealType>
__global__ void axpyOwnedKernel(RealType* out, const RealType* a, const RealType* b,
                                RealType alpha,
                                const int* nodeToDof, const uint8_t* ownership,
                                size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1 || nodeToDof[i] < 0) return;
    out[i] = a[i] + alpha * b[i];
}

// DDT pressure RHS: b = -coef * divAccNode  (coef = rho/dt).
// Sign comes from the projection algebra. D^T as written
// (applyDivTransposePerNodeKernel: same-sign 0.5*(p[L]-p[R])*A to both L,R)
// integrates to  D^T p ~ -V * grad(p)  so  M^{-1} D^T ~ -grad. The corrector
// u^{n+1} = u** - (dt/rho) * gradPhiq  with gradPhiq = M^{-1} D^T phi  then
// adds  +(dt/rho) grad(phi) to u**. Closing  D u^{n+1} = D u** + (dt/rho) A phi
// to zero requires  A phi = -(rho/dt) D u**, i.e. RHS has negative sign.
// Same magnitude/sign convention as the K-path's buildPressureRhsKernel.
template<typename RealType>
__global__ void buildPressureRhsDDTKernel(const RealType* divAccNode,
                                          const int* nodeToDof,
                                          const uint8_t* ownership,
                                          RealType coef, RealType* rhsNode,
                                          size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1 || nodeToDof[i] < 0) { rhsNode[i] = RealType(0); return; }
    // Confirmed by sign-discrimination test (mars_amr_ddt.cu --test=sign):
    // A models the standard -Laplacian (all three axes give (A*phi)/V_i ~ -2
    // for phi=coord^2). Closing the projection D u^{n+1}=0 with corrector
    // u = u** - (dt/rho) grad(phi) requires A phi = -(rho/dt) D u**.
    rhsNode[i] = -coef * divAccNode[i];
}

// Unpreconditioned matrix-free CG for (D M^{-1} D^T) phi = b. Returns iters
// at convergence, -2 on max-iter. Mirrors runCgTest from mars_amr_ddt.cu
// byte-for-byte (rho on <r,r>, convergence on sqrt(rho)=|r|, owned-only dots).
// Pin/mask handling: b[i]=0 at every Dirichlet-pressure DOF (cavity single
// pin or channel outflow mask), matching the SpMV's identity-row rule
// out[i]=phi[i] at those slots. A production version adds Jacobi/AMG.
template<typename KeyType, typename RealType>
int solvePressureDDT(NSStepper<KeyType, RealType>& s,
                     cstone::DeviceVector<RealType>& b_node,
                     cstone::DeviceVector<RealType>& phi_node)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;

    cstone::DeviceVector<RealType> r(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> p(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Ap(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));

    // Fresh start: phi = 0. With b[i]=0 at every pinned DOF and out[i]=phi[i]
    // in the SpMV, the pinned-row equation reads 0 = phi[i] => phi[i] stays 0.
    thrust::fill(thrust::device_pointer_cast(phi_node.data()),
                 thrust::device_pointer_cast(phi_node.data() + s.nodeCount), RealType(0));

    const int pinDof          = s.pressurePinDof;
    const uint8_t* maskPtr    = (s.d_isPressureBdryDof.size() > 0)
                                ? s.d_isPressureBdryDof.data() : nullptr;
    if (pinDof >= 0 || maskPtr != nullptr)
    {
        const int* dofPtr     = s.d_node_to_dof.data();
        const uint8_t* ownPtr = d_nodeOwnership.data();
        RealType* bPtr        = b_node.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [pinDof, maskPtr, dofPtr, ownPtr, bPtr] __device__ (size_t i) {
                if (ownPtr[i] != 1) return;
                int dof = dofPtr[i];
                if (dof < 0) return;
                bool flagged = (pinDof >= 0 && dof == pinDof) ||
                               (maskPtr != nullptr && maskPtr[dof] != 0);
                if (flagged) bPtr[i] = RealType(0);
            });
        cudaDeviceSynchronize();
    }

    // r = b - A*phi  (phi=0 so r = b on owned slots). p = r.
    thrust::copy(thrust::device_pointer_cast(b_node.data()),
                 thrust::device_pointer_cast(b_node.data() + s.nodeCount),
                 thrust::device_pointer_cast(r.data()));
    thrust::copy(thrust::device_pointer_cast(r.data()),
                 thrust::device_pointer_cast(r.data() + s.nodeCount),
                 thrust::device_pointer_cast(p.data()));

    RealType rho_old = ownedDot<RealType>(r, r, s.d_node_to_dof,
                                          d_nodeOwnership.data(), s.nodeCount);
    RealType r0_norm = std::sqrt(rho_old);
    if (r0_norm < std::numeric_limits<RealType>::min()) return 0;
    const RealType absTol = s.tolerance * r0_norm;

    int iters = -2;
    // Live progress every N CG iterations on rank 0 so large-mesh runs show
    // residual-decay rate instead of going silent for minutes. Each line
    // prints |r|/|r0| so the user can extrapolate iters-to-convergence.
    const int liveEvery = (s.maxIter >= 500) ? 100 : 25;
    for (int it = 0; it < s.maxIter; ++it)
    {
        // Ghost slots of p must be valid: A's D^T reads p[L],p[R] at every face.
        s.domain.exchangeNodeHalo(p);
        applyDDTPerNode<KeyType, RealType>(s, p, Ap, gx, gy, gz);

        RealType pAp = ownedDot<RealType>(p, Ap, s.d_node_to_dof,
                                          d_nodeOwnership.data(), s.nodeCount);
        if (pAp <= RealType(0)) { iters = -2; break; }
        RealType alpha = rho_old / pAp;

        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            phi_node.data(), phi_node.data(), p.data(), alpha,
            s.d_node_to_dof.data(), d_nodeOwnership.data(), s.nodeCount);
        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            r.data(), r.data(), Ap.data(), -alpha,
            s.d_node_to_dof.data(), d_nodeOwnership.data(), s.nodeCount);
        cudaDeviceSynchronize();

        RealType rho_new = ownedDot<RealType>(r, r, s.d_node_to_dof,
                                              d_nodeOwnership.data(), s.nodeCount);
        if (s.rank == 0 && (it == 0 || (it + 1) % liveEvery == 0))
        {
            std::cout << "    [cg-ddt] iter " << std::setw(6) << (it + 1)
                      << "  |r|/|r0| = " << std::scientific << std::setprecision(3)
                      << (std::sqrt(rho_new) / std::max(r0_norm, std::numeric_limits<RealType>::min()))
                      << std::defaultfloat << "\n" << std::flush;
        }
        if (std::sqrt(rho_new) < absTol) { iters = it + 1; break; }
        RealType beta = rho_new / rho_old;
        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            p.data(), r.data(), p.data(), beta,
            s.d_node_to_dof.data(), d_nodeOwnership.data(), s.nodeCount);
        cudaDeviceSynchronize();

        // (Inner-loop ortho-projection was tried here based on the
        // Nek5000/MFEM/PETSc pattern but caused immediate divergence in
        // this codebase. Most likely cause: the per-node CG vector `p`
        // includes ghost slots that are not yet halo-synced when the mean
        // is computed, biasing the projection. Disabled; outer-loop
        // broadcasts at all 4 sync points are the working path. See
        // TGV_HANDOFF.md for the full diagnosis and the next-session plan
        // to do this correctly (mass-weighted mean, owned-only, separate
        // halo sync per iter).)
        rho_old = rho_new;
    }

    s.domain.exchangeNodeHalo(phi_node);
    return iters;
}

// Negate three per-node vectors in place on owned DOFs only. Used when the
// gradient was produced by applyDivTransposePerNodeKernel (which integrates
// to -V*grad), so after normalize we have -grad; flipping the sign here makes
// the result a true +grad that predictor/corrector kernels (which subtract
// dt/rho * grad) can consume without their own sign awareness. Ghost slots
// are zeroed -- caller does a forward halo if needed.
template<typename RealType>
__global__ void negateThreeOwnedKernel(RealType* gx, RealType* gy, RealType* gz,
                                       const int* nodeToDof, const uint8_t* ownership,
                                       size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) { gx[i] = gy[i] = gz[i] = RealType(0); return; }
    int dof = nodeToDof[i];
    if (dof < 0) { gx[i] = gy[i] = gz[i] = RealType(0); return; }
    gx[i] = -gx[i];
    gy[i] = -gy[i];
    gz[i] = -gz[i];
}

// =============================================================================
// Periodic post-scatter helper. With the current DOF-collapse implementation
// (nodeToDof[slave] = nodeToDof[master]), slave and master nodes write into
// the same DOF in every per-DOF kernel. So no extra pair-summing is needed at
// the accumulator level -- the DOF mapping already takes care of it. This
// helper is a NO-OP and kept only so the call sites can stay in place across
// the four sub-steps.
// =============================================================================
template<typename KeyType, typename RealType>
inline void maybePeriodicSum(NSStepper<KeyType, RealType>& /*s*/,
                             cstone::DeviceVector<RealType>& /*d_field*/)
{
}

// =============================================================================
// Diagnostic helpers used by the sub-steps below. RMS of a 3-vector magnitude
// and max of |scalar| over owned interior DOFs (boundary skipped because lid
// Dirichlet picks up an artificial step). Defined above the step routines so
// they are visible at template instantiation time.
// =============================================================================
template<typename KeyType, typename RealType>
RealType rmsOwnedInterior3(NSStepper<KeyType, RealType>& s,
                            const cstone::DeviceVector<RealType>& d_gx,
                            const cstone::DeviceVector<RealType>& d_gy,
                            const cstone::DeviceVector<RealType>& d_gz)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const uint8_t* ownPtr = d_nodeOwnership.data();
    const int* dofPtr     = s.d_node_to_dof.data();
    const uint8_t* bnd    = s.d_isBdryDof.data();
    const uint8_t* pmask  = (s.d_isPressureBdryDof.size() > 0)
                            ? s.d_isPressureBdryDof.data() : nullptr;
    const RealType* gxP   = d_gx.data();
    const RealType* gyP   = d_gy.data();
    const RealType* gzP   = d_gz.data();
    RealType locSumSq = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount),
        [ownPtr, dofPtr, bnd, pmask, gxP, gyP, gzP] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1) return RealType(0);
            int dof = dofPtr[i];
            if (dof < 0 || bnd[dof])             return RealType(0);
            if (pmask != nullptr && pmask[dof])  return RealType(0);
            RealType gx = gxP[i], gy = gyP[i], gz = gzP[i];
            return gx*gx + gy*gy + gz*gz;
        }, RealType(0), thrust::plus<RealType>());
    long long locCnt = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount),
        [ownPtr, dofPtr, bnd, pmask] __device__ (size_t i) -> long long {
            if (ownPtr[i] != 1) return 0LL;
            int dof = dofPtr[i];
            if (dof < 0 || bnd[dof])             return 0LL;
            if (pmask != nullptr && pmask[dof])  return 0LL;
            return 1LL;
        }, 0LL, thrust::plus<long long>());
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    RealType gSum = 0;
    long long gCnt = 0;
    MPI_Allreduce(&locSumSq, &gSum, 1, mpiType,       MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locCnt,   &gCnt, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return (gCnt > 0) ? std::sqrt(gSum / RealType(gCnt)) : RealType(0);
}

// RMS of a scalar per-node field over interior owned DOFs (skip velocity-
// Dirichlet AND pressure-Dirichlet DOFs). Complements maxOwnedInteriorAbs:
// max is dominated by 1-cell ring next to Dirichlet pressure outflow; RMS
// averages it down so the bulk-interior divergence is visible.
template<typename KeyType, typename RealType>
RealType rmsOwnedInterior1(NSStepper<KeyType, RealType>& s,
                            const cstone::DeviceVector<RealType>& d_q)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const uint8_t* ownPtr = d_nodeOwnership.data();
    const int* dofPtr     = s.d_node_to_dof.data();
    const uint8_t* bnd    = s.d_isBdryDof.data();
    const uint8_t* pmask  = (s.d_isPressureBdryDof.size() > 0)
                            ? s.d_isPressureBdryDof.data() : nullptr;
    const RealType* qP    = d_q.data();
    RealType locSumSq = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount),
        [ownPtr, dofPtr, bnd, pmask, qP] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1) return RealType(0);
            int dof = dofPtr[i];
            if (dof < 0 || bnd[dof])             return RealType(0);
            if (pmask != nullptr && pmask[dof])  return RealType(0);
            RealType q = qP[i];
            return q * q;
        }, RealType(0), thrust::plus<RealType>());
    long long locCnt = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount),
        [ownPtr, dofPtr, bnd, pmask] __device__ (size_t i) -> long long {
            if (ownPtr[i] != 1) return 0LL;
            int dof = dofPtr[i];
            if (dof < 0 || bnd[dof])             return 0LL;
            if (pmask != nullptr && pmask[dof])  return 0LL;
            return 1LL;
        }, 0LL, thrust::plus<long long>());
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    RealType gSum = 0;
    long long gCnt = 0;
    MPI_Allreduce(&locSumSq, &gSum, 1, mpiType,       MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locCnt,   &gCnt, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return (gCnt > 0) ? std::sqrt(gSum / RealType(gCnt)) : RealType(0);
}

template<typename KeyType, typename RealType>
RealType maxOwnedInteriorAbs(NSStepper<KeyType, RealType>& s,
                             const cstone::DeviceVector<RealType>& d_q)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const uint8_t* ownPtr = d_nodeOwnership.data();
    const int* dofPtr     = s.d_node_to_dof.data();
    const uint8_t* bnd    = s.d_isBdryDof.data();
    const uint8_t* pmask  = (s.d_isPressureBdryDof.size() > 0)
                            ? s.d_isPressureBdryDof.data() : nullptr;
    const RealType* qP    = d_q.data();
    RealType locMax = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount),
        [ownPtr, dofPtr, bnd, pmask, qP] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1) return RealType(0);
            int dof = dofPtr[i];
            if (dof < 0 || bnd[dof])             return RealType(0);
            if (pmask != nullptr && pmask[dof])  return RealType(0);
            return fabs(qP[i]);
        }, RealType(0), thrust::maximum<RealType>());
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    RealType gMax = 0;
    MPI_Allreduce(&locMax, &gMax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    return gMax;
}

// =============================================================================
// runNsStep is composed of four sub-steps. They share NO state outside
// NSStepper; each can be invoked in isolation for unit testing (e.g. drive
// only runPressureSolveStep with a known velocity field to reproduce the
// DDT bug without setting up a full time loop).
//
// The sub-steps assume the previous one has left the relevant fields in a
// consistent (ghost-synced) state. runNsStep handles the chaining; isolated
// callers must do their own halo sync of inputs.
// =============================================================================

template<typename KeyType, typename RealType>
void runPredictorStep(NSStepper<KeyType, RealType>& s, RealType dt, RealType rho)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    auto c0 = std::get<0>(d_conn).data(); auto c1 = std::get<1>(d_conn).data();
    auto c2 = std::get<2>(d_conn).data(); auto c3 = std::get<3>(d_conn).data();
    auto c4 = std::get<4>(d_conn).data(); auto c5 = std::get<5>(d_conn).data();
    auto c6 = std::get<6>(d_conn).data(); auto c7 = std::get<7>(d_conn).data();
    const size_t startElem = s.domain.startIndex();
    const size_t numLocal  = s.domain.localElementCount();
    const int nodeBlocks   = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    const int eBlocks      = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
    const RealType invRho  = RealType(1) / rho;

    s.domain.exchangeNodeHalo(s.d_u);
    s.domain.exchangeNodeHalo(s.d_v);
    s.domain.exchangeNodeHalo(s.d_w);
    s.domain.exchangeNodeHalo(s.d_p);

    // Periodic: re-broadcast master->slave for u, v, w, p at the START of
    // each step (before the predictor's grad p and advection scatters read
    // them). Across steps, even after the post-corrector broadcast and the
    // halo exchanges, slave-node values can drift slightly from master from
    // floating-point reduction ordering and the per-node implicit-diffusion
    // solve. If grad(p) and the advection mdot at periodic faces read
    // mismatched slave/master values, the discrete operator is asymmetric
    // and pressure-velocity coupling drifts -- the classic "pressure
    // doubles every 2-3 steps while KE stays flat" signature. Found by
    // multi-agent DDT-path audit (round 7).
    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Periodic && s.periodicMap)
    {
        int blk = 256, grd = int((s.nodeCount + blk - 1) / blk);
        const int* d_partner = s.periodicMap->d_periodicPartner.data();
        size_t nN            = s.nodeCount;
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(d_partner, nN, s.d_u.data());
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(d_partner, nN, s.d_v.data());
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(d_partner, nN, s.d_w.data());
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(d_partner, nN, s.d_p.data());
        cudaDeviceSynchronize();
        s.domain.exchangeNodeHalo(s.d_u);
        s.domain.exchangeNodeHalo(s.d_v);
        s.domain.exchangeNodeHalo(s.d_w);
        s.domain.exchangeNodeHalo(s.d_p);
    }

    // grad p^n -> (d_gradPx, d_gradPy, d_gradPz). Default uses div^T (Path A);
    // --legacy-gradient falls back to the SCS-face gradient.
    {
        cstone::DeviceVector<RealType> d_gxAcc(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> d_gyAcc(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> d_gzAcc(s.nodeCount, RealType(0));
        if (eBlocks > 0)
        {
            if (s.useLegacyGradient)
            {
                computeGradientPerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    s.d_p.data(),
                    s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                    d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
                    startElem, numLocal);
            }
            else
            {
                applyDivTransposePerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    s.d_p.data(),
                    s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                    d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
                    startElem, numLocal);
            }
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(d_gxAcc);
        s.domain.reverseExchangeNodeHaloAdd(d_gyAcc);
        s.domain.reverseExchangeNodeHaloAdd(d_gzAcc);
        maybePeriodicSum<KeyType, RealType>(s, d_gxAcc);
        maybePeriodicSum<KeyType, RealType>(s, d_gyAcc);
        maybePeriodicSum<KeyType, RealType>(s, d_gzAcc);
        normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
            s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_gradPx.data(), s.d_gradPy.data(), s.d_gradPz.data(),
            s.nodeCount);
        cudaDeviceSynchronize();
        // applyDivTransposePerNodeKernel integrates to -V*grad (so M^{-1} D^T p
        // produces -grad p, not +grad p). The predictor kernel expects gradPnq
        // = +grad p (it computes u* = u^n - dt*(...+ grad(p)/rho - ...)). The
        // corrector applies the same flip at line 4135. Without this flip the
        // pressure gradient is added with the wrong sign in the predictor,
        // which drives the divergence to grow by ~10x per step regardless of
        // pressure stabilization or time integrator.
        if (!s.useLegacyGradient)
        {
            negateThreeOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                s.d_gradPx.data(), s.d_gradPy.data(), s.d_gradPz.data(),
                s.d_node_to_dof.data(), d_nodeOwnership.data(), s.nodeCount);
            cudaDeviceSynchronize();
        }
        s.lastGradPRms = rmsOwnedInterior3<KeyType, RealType>(
            s, s.d_gradPx, s.d_gradPy, s.d_gradPz);
    }

    // Per-component advection scatter + predictor apply. Each scatter writes
    // into a PERSISTENT advection-accumulator (s.d_advU_n etc.) so it can be
    // reused as N(u^(n-1)) at step n+1 via EXT2 extrapolation. Predictor
    // dispatches on s.bdfStep: 0 = BDF1 (Chorin), >=1 = BDF2/EXT2.
    auto runPredictor = [&] (cstone::DeviceVector<RealType>& qN,
                              cstone::DeviceVector<RealType>& qNm1,
                              cstone::DeviceVector<RealType>& qStar,
                              cstone::DeviceVector<RealType>& advN,
                              cstone::DeviceVector<RealType>& advNm1,
                              cstone::DeviceVector<RealType>& gradPnq,
                              cstone::DeviceVector<RealType>& qTarget)
    {
        // Compute current-step advection into advN.
        if (advN.size() != s.nodeCount) advN.resize(s.nodeCount);
        thrust::fill(thrust::device_pointer_cast(advN.data()),
                     thrust::device_pointer_cast(advN.data() + s.nodeCount),
                     RealType(0));
        if (eBlocks > 0)
        {
            explicitAdvectionFluxScatterPerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7,
                s.d_u.data(), s.d_v.data(), s.d_w.data(),
                qN.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                advN.data(), startElem, numLocal,
                s.useSkewSymmetricAdvection);
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(advN);
        maybePeriodicSum<KeyType, RealType>(s, advN);

        if (s.useBdf2 && s.bdfStep >= 1 && advNm1.size() == s.nodeCount)
        {
            applyPredictorBdf2PerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                qStar.data(), qN.data(), qNm1.data(),
                advN.data(), advNm1.data(),
                gradPnq.data(), s.d_mass.data(), qTarget.data(),
                s.d_isBdryDof.data(), s.d_node_to_dof.data(),
                d_nodeOwnership.data(),
                dt, invRho, s.nodeCount);
        }
        else
        {
            applyPredictorPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                qStar.data(), qN.data(), advN.data(),
                gradPnq.data(), s.d_mass.data(), s.d_massNode.data(), qTarget.data(),
                s.d_isBdryDof.data(), s.d_node_to_dof.data(),
                d_nodeOwnership.data(),
                dt, invRho, s.nodeCount);
        }
        cudaDeviceSynchronize();
    };

    runPredictor(s.d_u, s.d_u_nm1, s.d_uStar, s.d_advU_n, s.d_advU_nm1, s.d_gradPx, s.d_uTarget);
    runPredictor(s.d_v, s.d_v_nm1, s.d_vStar, s.d_advV_n, s.d_advV_nm1, s.d_gradPy, s.d_vTarget);
    runPredictor(s.d_w, s.d_w_nm1, s.d_wStar, s.d_advW_n, s.d_advW_nm1, s.d_gradPz, s.d_wTarget);

    // Sync ghosts of q* so the implicit RHS / CG warm-start read correct ghosts.
    s.domain.exchangeNodeHalo(s.d_uStar);
    s.domain.exchangeNodeHalo(s.d_vStar);
    s.domain.exchangeNodeHalo(s.d_wStar);

    // Periodic: broadcast master->slave so u* is identical on both sides of
    // the periodic boundary. Without this, the implicit diffusion RHS sees
    // a non-periodic u* and produces a non-periodic u**, propagating drift.
    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Periodic && s.periodicMap)
    {
        int blk = 256, grd = int((s.nodeCount + blk - 1) / blk);
        const int* dp = s.periodicMap->d_periodicPartner.data();
        size_t nN     = s.nodeCount;
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(dp, nN, s.d_uStar.data());
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(dp, nN, s.d_vStar.data());
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(dp, nN, s.d_wStar.data());
        cudaDeviceSynchronize();
        s.domain.exchangeNodeHalo(s.d_uStar);
        s.domain.exchangeNodeHalo(s.d_vStar);
        s.domain.exchangeNodeHalo(s.d_wStar);
    }
}

// -------------------------------------------------------------------------
// Step 2: IMPLICIT DIFFUSION. (M/dt + nu K) q** = (M/dt) q*  for each
// component. The matrix is FROZEN (assembled once at setup); per step we
// rebuild only the RHS. Velocity Dirichlet lift-off is applied to the RHS
// for each component using its own qTarget array.
// -------------------------------------------------------------------------
template<typename KeyType, typename RealType>
void runImplicitDiffusionStep(NSStepper<KeyType, RealType>& s, RealType dt)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const int nodeBlocks  = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    // BDF1 (step 0): mass coefficient = 1/dt, matrix = Avel.
    // BDF2 (step >= 1): mass coefficient = 3/(2*dt), matrix = Avel_bdf2.
    const bool bdf2Active = (s.useBdf2 && s.bdfStep >= 1 && s.d_valuesVel_bdf2.size() > 0);
    const RealType invDt  = bdf2Active ? (RealType(3) / (RealType(2) * dt))
                                       : (RealType(1) / dt);

    cstone::DeviceVector<RealType> b(s.numOwnedDofs);
    cstone::DeviceVector<RealType> xVec(s.numTotalDofs);

    auto runImplicit = [&] (cstone::DeviceVector<RealType>& qStar,
                             cstone::DeviceVector<RealType>& qStarStar,
                             cstone::DeviceVector<RealType>& qTarget) -> int
    {
        // RHS = (mass coef) * mass * qStar; coef = 1/dt for BDF1, 3/(2dt) for BDF2.
        thrust::fill(thrust::device_pointer_cast(b.data()),
                     thrust::device_pointer_cast(b.data() + s.numOwnedDofs),
                     RealType(0));
        buildVelocityImplicitRhsKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            qStar.data(), s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_mass.data(), invDt, b.data(), s.nodeCount);
        cudaDeviceSynchronize();

        enforceBcRhsFromTargetKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            s.d_isBdryDof.data(), qTarget.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            b.data(), s.nodeCount);
        cudaDeviceSynchronize();

        // Warm-start x with qStar scattered into DOF order.
        cudaMemset(xVec.data(), 0, s.numTotalDofs * sizeof(RealType));
        thrust::for_each(thrust::device,
                          thrust::counting_iterator<size_t>(0),
                          thrust::counting_iterator<size_t>(s.nodeCount),
                          [nodeToDof = s.d_node_to_dof.data(),
                           src       = qStar.data(),
                           xOut      = xVec.data()] __device__(size_t i) {
                              int dof = nodeToDof[i];
                              if (dof >= 0) xOut[dof] = src[i];
                          });
        cudaDeviceSynchronize();

        return solveOneComponent<KeyType, RealType>(
            s, b, xVec, qStarStar,
            bdf2Active ? s.Avel_bdf2 : s.Avel);
    };

    s.lastUIters = runImplicit(s.d_uStar, s.d_uStarStar, s.d_uTarget);
    s.lastVIters = runImplicit(s.d_vStar, s.d_vStarStar, s.d_vTarget);
    s.lastWIters = runImplicit(s.d_wStar, s.d_wStarStar, s.d_wTarget);

    // Sync ghosts of q** for the divergence scatter (face donors read q** at
    // ghost corners on rank boundaries).
    s.domain.exchangeNodeHalo(s.d_uStarStar);
    s.domain.exchangeNodeHalo(s.d_vStarStar);
    s.domain.exchangeNodeHalo(s.d_wStarStar);

    // Periodic: broadcast master->slave so u** is identical on both sides.
    // The divergence scatter reads u** at face corners; an asymmetric u**
    // produces an asymmetric div(u**) source, which biases the pressure
    // solve and feeds the runaway. See mars_periodic_bc.hpp + post-corrector
    // broadcast for the matching pattern at u^{n+1}.
    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Periodic && s.periodicMap)
    {
        int blk = 256, grd = int((s.nodeCount + blk - 1) / blk);
        const int* dp = s.periodicMap->d_periodicPartner.data();
        size_t nN     = s.nodeCount;
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(dp, nN, s.d_uStarStar.data());
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(dp, nN, s.d_vStarStar.data());
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(dp, nN, s.d_wStarStar.data());
        cudaDeviceSynchronize();
        s.domain.exchangeNodeHalo(s.d_uStarStar);
        s.domain.exchangeNodeHalo(s.d_vStarStar);
        s.domain.exchangeNodeHalo(s.d_wStarStar);
    }
}

// -------------------------------------------------------------------------
// Step 3: PRESSURE POISSON. K phi = (rho/dt) div(u**). The un-normalized
// divAccNode is exactly the integrated source f_i*V_i, so the lumped RHS is
// rhs[dof] = (rho/dt) * divAccNode[node(dof)] -- no separate
// normalize-then-multiply-by-V chain. Enforce the pin (single Dirichlet DOF
// at the chosen corner) on the RHS each step.
//
// Can be called in isolation with a pre-populated (s.d_uStarStar,
// s.d_vStarStar, s.d_wStarStar) and ghosts in sync. Useful for reproducing
// the DDT bug without the surrounding momentum solve.
// -------------------------------------------------------------------------
template<typename KeyType, typename RealType>
void runPressureSolveStep(NSStepper<KeyType, RealType>& s, RealType dt, RealType rho)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    auto c0 = std::get<0>(d_conn).data(); auto c1 = std::get<1>(d_conn).data();
    auto c2 = std::get<2>(d_conn).data(); auto c3 = std::get<3>(d_conn).data();
    auto c4 = std::get<4>(d_conn).data(); auto c5 = std::get<5>(d_conn).data();
    auto c6 = std::get<6>(d_conn).data(); auto c7 = std::get<7>(d_conn).data();
    const size_t startElem = s.domain.startIndex();
    const size_t numLocal  = s.domain.localElementCount();
    const int nodeBlocks   = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    const int eBlocks      = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
    const RealType invDt   = RealType(1) / dt;

    cstone::DeviceVector<RealType> b(s.numOwnedDofs);
    cstone::DeviceVector<RealType> xVec(s.numTotalDofs);

    cstone::DeviceVector<RealType> d_divAccNode(s.nodeCount, RealType(0));
    if (eBlocks > 0)
    {
        if (s.useRhieChow)
        {
            // Rhie-Chow corrected divergence: adds the pressure-velocity
            // coupling term that suppresses checkerboard on Q1-Q1 periodic
            // grids. This is the standard CVFEM/co-located-FV stabilization
            // used by Nalu-Wind and OpenFOAM. tau auto = dt/rho.
            RealType tauRC = s.rhieChowTau;
            if (tauRC <= 0) tauRC = dt / rho;
            const auto& d_x = s.domain.getNodeX();
            const auto& d_y = s.domain.getNodeY();
            const auto& d_z = s.domain.getNodeZ();
            computeDivergenceRhieChowKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7,
                s.d_uStarStar.data(), s.d_vStarStar.data(), s.d_wStarStar.data(),
                s.d_p.data(),
                s.d_gradPx.data(), s.d_gradPy.data(), s.d_gradPz.data(),
                d_x.data(), d_y.data(), d_z.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                tauRC,
                d_divAccNode.data(), startElem, numLocal);
        }
        else
        {
            computeDivergencePerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7,
                s.d_uStarStar.data(), s.d_vStarStar.data(), s.d_wStarStar.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                d_divAccNode.data(), startElem, numLocal);
        }
        cudaDeviceSynchronize();
    }
    s.domain.reverseExchangeNodeHaloAdd(d_divAccNode);
    maybePeriodicSum<KeyType, RealType>(s, d_divAccNode);

    // Diagnostic: |div(u**)| max BEFORE the corrector. This is the source term
    // the pressure solve must cancel. Compare with lastDivMax (post-corrector)
    // to see how much the projection actually reduces divergence.
    //
    // Also stash the normalized per-node divergence into s.d_divUStar so the
    // rotational pressure correction in runCorrectorStep can read it without
    // recomputing the scatter. When rotational correction is off the field is
    // simply unused.
    {
        if (s.d_divUStar.size() != s.nodeCount) s.d_divUStar.resize(s.nodeCount);
        thrust::fill(thrust::device_pointer_cast(s.d_divUStar.data()),
                     thrust::device_pointer_cast(s.d_divUStar.data() + s.nodeCount),
                     RealType(0));
        normalizeDivergencePerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_divAccNode.data(), s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_divUStar.data(), s.nodeCount);
        cudaDeviceSynchronize();
        s.lastDivMaxPre = maxOwnedInteriorAbs<KeyType, RealType>(s, s.d_divUStar);
    }

    if (s.pressureSolve == PressureSolveKind::K)
    {
        // K-path: Galerkin Laplacian (K is for -Δu=f, so RHS = -coef * divAcc).
        thrust::fill(thrust::device_pointer_cast(b.data()),
                     thrust::device_pointer_cast(b.data() + s.numOwnedDofs), RealType(0));
        RealType coef = rho * invDt;
        buildPressureRhsKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_divAccNode.data(), s.d_node_to_dof.data(),
            d_nodeOwnership.data(), coef, b.data(), s.nodeCount);
        cudaDeviceSynchronize();
        // Pressure BC RHS: cavity zeroes the single pin DOF; channel zeroes
        // every owned Dirichlet-pressure DOF on the outflow face.
        if (s.pressurePinDof >= 0)
        {
            enforcePinRhsKernel<RealType><<<1, 1>>>(s.pressurePinDof, b.data());
            cudaDeviceSynchronize();
        }
        if (s.d_isPressureBdryDof.size() > 0)
        {
            int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
            enforcePressureBcRhsKernel<RealType><<<dofBlocks, s.blockSize>>>(
                s.d_isPressureBdryDof.data(), b.data(), s.numOwnedDofs);
            cudaDeviceSynchronize();
        }
        // Periodic: project the RHS onto range(K) by subtracting its global
        // mean. The pure-Neumann Laplacian K has a constant-mode null space,
        // so K phi = b is solvable only when sum(b) = 0. Without this, CG hits
        // pAp <= 0 on the null-space direction and bails.
        if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Periodic && s.numOwnedDofs > 0)
        {
            auto bp = thrust::device_pointer_cast(b.data());
            RealType localSum = thrust::reduce(thrust::device, bp, bp + s.numOwnedDofs,
                                                RealType(0), thrust::plus<RealType>());
            RealType globalSum = 0;
            long long localN = s.numOwnedDofs, globalN = 0;
            MPI_Datatype mpiR = std::is_same<RealType, double>::value ? MPI_DOUBLE : MPI_FLOAT;
            MPI_Allreduce(&localSum, &globalSum, 1, mpiR,        MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&localN,   &globalN,   1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            if (globalN > 0)
            {
                RealType mean = globalSum / RealType(globalN);
                thrust::transform(thrust::device, bp, bp + s.numOwnedDofs, bp,
                                  [mean] __device__ (RealType v) { return v - mean; });
            }
        }
        // Fresh warm-start: phi is per-step correction, not cumulative.
        cudaMemset(xVec.data(), 0, s.numTotalDofs * sizeof(RealType));
        s.lastPressureIters = solveOneComponent<KeyType, RealType>(s, b, xVec, s.d_phi, s.Apre);
    }
    else  // DDT: (D M^{-1} D^T) phi = -(rho/dt) D u**
    {
        // DIAGNOSTIC: env MARS_DDT_USE_ASSEMBLED_CG=1 forces CG to use the
        // assembled s.AddT matrix. With the buildDDTSparsity dedicated
        // pattern (cross-element 2-hop pairs included) and the node-driven
        // assembleDDTPerNodeKernel, the assembled matrix equals the matrix-
        // free apply bit-for-bit -- so this diagnostic should give identical
        // numbers to matrix-free CG (and use the same iteration count).
        const char* useAsmEv = std::getenv("MARS_DDT_USE_ASSEMBLED_CG");
        bool useAssembledCG = (useAsmEv && std::string(useAsmEv) != "0")
                              && (s.solverKind == SolverKind::CG);
        if (s.solverKind == SolverKind::Hypre || useAssembledCG)
        {
            // DDT + (Hypre or assembled-CG): use the assembled DDT matrix
            // (s.AddT, built at setup). Build per-OWNED-DOF b exactly like the
            // K path does.
            thrust::fill(thrust::device_pointer_cast(b.data()),
                         thrust::device_pointer_cast(b.data() + s.numOwnedDofs),
                         RealType(0));
            RealType coef = rho * invDt;
            buildPressureRhsKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                d_divAccNode.data(), s.d_node_to_dof.data(),
                d_nodeOwnership.data(), coef, b.data(), s.nodeCount);
            cudaDeviceSynchronize();
            if (s.pressurePinDof >= 0)
            {
                enforcePinRhsKernel<RealType><<<1, 1>>>(s.pressurePinDof, b.data());
                cudaDeviceSynchronize();
            }
            if (s.d_isPressureBdryDof.size() > 0)
            {
                int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
                enforcePressureBcRhsKernel<RealType><<<dofBlocks, s.blockSize>>>(
                    s.d_isPressureBdryDof.data(), b.data(), s.numOwnedDofs);
                cudaDeviceSynchronize();
            }
            cudaMemset(xVec.data(), 0, s.numTotalDofs * sizeof(RealType));
            // DDT operator (D M^-1 D^T) is SPSD with constant null space; Hypre
            // PCG rejects it (error 1). GMRES tolerates it. Hint passes through
            // to solveOneComponent which selects HypreGMRESSolver wrapper.
            s.lastPressureIters = solveOneComponent<KeyType, RealType>(
                s, b, xVec, s.d_phi, s.AddT, KrylovHint::GMRES);
            // DIAGNOSTIC: sample WHOLE vectors via D2H, compute max-abs on
            // host. Lets us see whether the solution actually propagated.
            // Active only when env MARS_DDT_DIAG_AFTER_HYPRE is set.
            if (std::getenv("MARS_DDT_DIAG_AFTER_HYPRE") && s.rank == 0
                && s.numOwnedDofs > 0)
            {
                std::vector<RealType> h_b(s.numOwnedDofs), h_x(s.numOwnedDofs);
                std::vector<RealType> h_phi(s.nodeCount);
                cudaMemcpy(h_b.data(),   b.data(),       s.numOwnedDofs * sizeof(RealType), cudaMemcpyDeviceToHost);
                cudaMemcpy(h_x.data(),   xVec.data(),    s.numOwnedDofs * sizeof(RealType), cudaMemcpyDeviceToHost);
                cudaMemcpy(h_phi.data(), s.d_phi.data(), s.nodeCount   * sizeof(RealType), cudaMemcpyDeviceToHost);
                RealType bMax = 0, xMax = 0, phiMax = 0;
                int bMaxIdx = -1, xMaxIdx = -1;
                for (int i = 0; i < s.numOwnedDofs; ++i) {
                    if (std::abs(h_b[i]) > bMax) { bMax = std::abs(h_b[i]); bMaxIdx = i; }
                    if (std::abs(h_x[i]) > xMax) { xMax = std::abs(h_x[i]); xMaxIdx = i; }
                }
                for (size_t i = 0; i < s.nodeCount; ++i)
                    if (std::abs(h_phi[i]) > phiMax) phiMax = std::abs(h_phi[i]);
                std::cout << "  [DDT-diag] |b|_max=" << bMax << " (at dof " << bMaxIdx << ")"
                          << "  |xVec|_max=" << xMax << " (at dof " << xMaxIdx << ")"
                          << "  |d_phi|_max=" << phiMax
                          << "  iters=" << s.lastPressureIters << "\n";
            }
        }
        else
        {
            cstone::DeviceVector<RealType> d_bNode(s.nodeCount, RealType(0));
            RealType coef = rho * invDt;
            buildPressureRhsDDTKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                d_divAccNode.data(), s.d_node_to_dof.data(),
                d_nodeOwnership.data(), coef, d_bNode.data(), s.nodeCount);
            cudaDeviceSynchronize();
            // Periodic: D M^{-1} D^T is pure-Neumann with a constant null space.
            if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Periodic)
            {
                mars::fem::removeMean<RealType>(s.domain, d_bNode, MPI_COMM_WORLD);
            }
            s.lastPressureIters = solvePressureDDT<KeyType, RealType>(s, d_bNode, s.d_phi);
        }
    }

    // Periodic: pure-Neumann pressure has a constant-mode null space. Subtract
    // global mean so phi is uniquely defined (and the constant doesn't drift
    // step over step).
    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Periodic)
    {
        mars::fem::removeMean<RealType>(s.domain, s.d_phi, MPI_COMM_WORLD);
    }

    // Sync ghosts of phi -- needed by both the gradient (face donors) and the
    // pressure update on owned nodes that read ghost slots only for VTU.
    s.domain.exchangeNodeHalo(s.d_phi);

    // Periodic: broadcast phi master->slave so the corrector's grad(phi) reads
    // a periodic field. Without this, grad(phi) is asymmetric at the periodic
    // boundary, the corrector produces a non-divergence-free u^{n+1}, and the
    // next pressure source is correspondingly biased.
    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Periodic && s.periodicMap)
    {
        int blk = 256, grd = int((s.nodeCount + blk - 1) / blk);
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(
            s.periodicMap->d_periodicPartner.data(), s.nodeCount, s.d_phi.data());
        cudaDeviceSynchronize();
        s.domain.exchangeNodeHalo(s.d_phi);
    }
}

// -------------------------------------------------------------------------
// Step 4: CORRECTOR + pressure update + divergence diagnostic.
//   q^{n+1} = q** - (dt/rho) (grad phi)_q  on interior
//   q^{n+1} = qTarget                       on boundary (already correct)
//   p^{n+1} = p^n + phi
// grad(phi) computed once via the same B.4 kernel, applied per-component.
// Writes s.lastDivMax for monitoring.
// -------------------------------------------------------------------------
template<typename KeyType, typename RealType>
void runCorrectorStep(NSStepper<KeyType, RealType>& s, RealType dt, RealType rho)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    auto c0 = std::get<0>(d_conn).data(); auto c1 = std::get<1>(d_conn).data();
    auto c2 = std::get<2>(d_conn).data(); auto c3 = std::get<3>(d_conn).data();
    auto c4 = std::get<4>(d_conn).data(); auto c5 = std::get<5>(d_conn).data();
    auto c6 = std::get<6>(d_conn).data(); auto c7 = std::get<7>(d_conn).data();
    const size_t startElem = s.domain.startIndex();
    const size_t numLocal  = s.domain.localElementCount();
    const int nodeBlocks   = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    const int eBlocks      = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
    const RealType invRho  = RealType(1) / rho;
    // BDF2: corrector projection scaling matches the predictor's effective dt.
    // BDF1: u^{n+1} = u** - (dt/rho) grad(phi)
    // BDF2: u^{n+1} = u** - (2*dt/(3*rho)) grad(phi)
    const bool bdf2ActiveCorr = (s.useBdf2 && s.bdfStep >= 1 && s.d_valuesVel_bdf2.size() > 0);
    const RealType dtEff      = bdf2ActiveCorr ? (RealType(2) * dt / RealType(3)) : dt;

    {
        cstone::DeviceVector<RealType> d_gxAcc(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> d_gyAcc(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> d_gzAcc(s.nodeCount, RealType(0));
        // DDT mode MUST use the divT gradient (= literal D^T) to keep the
        // projection identity D u^{n+1} = 0 algebraically exact. K mode keeps
        // the legacy SCS gradient unless --experimental-divT overrides.
        const bool useDivT = (s.pressureSolve == PressureSolveKind::DDT) || !s.useLegacyGradient;
        if (eBlocks > 0)
        {
            if (!useDivT)
            {
                computeGradientPerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    s.d_phi.data(),
                    s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                    d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
                    startElem, numLocal);
            }
            else
            {
                applyDivTransposePerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    s.d_phi.data(),
                    s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                    d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
                    startElem, numLocal);
            }
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(d_gxAcc);
        s.domain.reverseExchangeNodeHaloAdd(d_gyAcc);
        s.domain.reverseExchangeNodeHaloAdd(d_gzAcc);
        maybePeriodicSum<KeyType, RealType>(s, d_gxAcc);
        maybePeriodicSum<KeyType, RealType>(s, d_gyAcc);
        maybePeriodicSum<KeyType, RealType>(s, d_gzAcc);
        normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
            s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_gradPhix.data(), s.d_gradPhiy.data(), s.d_gradPhiz.data(),
            s.nodeCount);
        cudaDeviceSynchronize();
        // applyDivTransposePerNodeKernel integrates to -V*grad (validator-
        // confirmed by mars_amr_ddt.cu --test=sign), so M^{-1} D^T phi = -grad.
        // The corrector kernel expects gradPhiq = +grad. Flip sign in place.
        // The SCS-gradient path produces +grad already, so no flip needed there.
        if (useDivT)
        {
            negateThreeOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                s.d_gradPhix.data(), s.d_gradPhiy.data(), s.d_gradPhiz.data(),
                s.d_node_to_dof.data(), d_nodeOwnership.data(), s.nodeCount);
            cudaDeviceSynchronize();
        }
        s.lastGradPhiRms = rmsOwnedInterior3<KeyType, RealType>(
            s, s.d_gradPhix, s.d_gradPhiy, s.d_gradPhiz);
    }

    auto runCorrector = [&] (cstone::DeviceVector<RealType>& qOut,
                              cstone::DeviceVector<RealType>& qStarStar,
                              cstone::DeviceVector<RealType>& gradPhiq,
                              cstone::DeviceVector<RealType>& qTarget)
    {
        applyCorrectorPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            qOut.data(), qStarStar.data(), gradPhiq.data(), qTarget.data(),
            s.d_isBdryDof.data(), s.d_node_to_dof.data(),
            d_nodeOwnership.data(), dtEff, invRho, s.nodeCount);
        cudaDeviceSynchronize();
    };

    runCorrector(s.d_u, s.d_uStarStar, s.d_gradPhix, s.d_uTarget);
    runCorrector(s.d_v, s.d_vStarStar, s.d_gradPhiy, s.d_vTarget);
    runCorrector(s.d_w, s.d_wStarStar, s.d_gradPhiz, s.d_wTarget);

    // Pressure update + ghost sync (next step's predictor needs ghost p^{n+1}).
    // Standard incremental Chorin: p^{n+1} = p^n + phi.
    updatePressureKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        s.d_p.data(), s.d_phi.data(),
        s.d_node_to_dof.data(), d_nodeOwnership.data(),
        s.nodeCount);
    cudaDeviceSynchronize();

    // Periodic: re-anchor pressure gauge after each cumulative incremental
    // update. The pressure equation has a constant-mode null space; phi has
    // its mean removed each step, but tiny floating-point residual means
    // accumulate through p^{n+1} = p^n + phi, biasing grad(p) in the next
    // predictor and producing an exponential pressure runaway (~1.5x per step
    // observed on cube16 TGV). Removing the mean of p here, once per step,
    // breaks the feedback.
    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Periodic)
    {
        mars::fem::removeMean<RealType>(s.domain, s.d_p, MPI_COMM_WORLD);
    }

    // Rotational pressure correction (Timmermans 1996, Guermond & Quartapelle 2000):
    //   p^{n+1} = p^n + phi - nu * div(u**)
    // The -nu*div(u**) term cures Chorin's O(dt^{1/2}) splitting error at boundaries
    // and periodic copies, raising pressure accuracy to O(dt^{3/2}). Strongly
    // recommended for periodic NS, where ordinary incremental Chorin is unstable
    // (no wall dissipation to absorb the splitting error).
    if (s.rotationalPressureCorrection && s.d_divUStar.size() == s.nodeCount)
    {
        const RealType nu = s.nuCached;
        RealType* pPtr            = s.d_p.data();
        const RealType* divPtr    = s.d_divUStar.data();
        const int* dofPtr2        = s.d_node_to_dof.data();
        const uint8_t* ownPtr2    = d_nodeOwnership.data();
        thrust::for_each(thrust::device,
                         thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(s.nodeCount),
                         [pPtr, divPtr, dofPtr2, ownPtr2, nu] __device__ (size_t i) {
                             if (ownPtr2[i] != 1) return;
                             int dof = dofPtr2[i];
                             if (dof < 0) return;
                             pPtr[i] -= nu * divPtr[i];
                         });
        cudaDeviceSynchronize();
    }

    s.domain.exchangeNodeHalo(s.d_u);
    s.domain.exchangeNodeHalo(s.d_v);
    s.domain.exchangeNodeHalo(s.d_w);
    s.domain.exchangeNodeHalo(s.d_p);

    // Periodic broadcast: per-NODE fields drift between slave and master
    // every step because the corrector and pressure update write each owned
    // node independently (DOF compaction makes the linear algebra unique but
    // not the per-node storage). Without this broadcast, slave-side grad(p)
    // in the next predictor reads a drifted pressure, generates an asymmetric
    // pressure-gradient force, and the system blows up around step 30.
    if (s.bcKind == NSStepper<KeyType, RealType>::BCKind::Periodic && s.periodicMap)
    {
        int blk = 256, grd = int((s.nodeCount + blk - 1) / blk);
        const int* d_partner = s.periodicMap->d_periodicPartner.data();
        size_t nN            = s.nodeCount;
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(d_partner, nN, s.d_u.data());
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(d_partner, nN, s.d_v.data());
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(d_partner, nN, s.d_w.data());
        mars::fem::periodicBroadcastKernel<<<grd, blk>>>(d_partner, nN, s.d_p.data());
        cudaDeviceSynchronize();
        s.domain.exchangeNodeHalo(s.d_u);
        s.domain.exchangeNodeHalo(s.d_v);
        s.domain.exchangeNodeHalo(s.d_w);
        s.domain.exchangeNodeHalo(s.d_p);
    }

    // -------------------------------------------------------------------------
    // div_max diagnostic over interior OWNED DOFs. Computes div(u^{n+1}); the
    // projection step's algebraic guarantee is div_max ~ roundoff every step.
    // -------------------------------------------------------------------------
    {
        cstone::DeviceVector<RealType> d_divAccPost(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> d_divNorm(s.nodeCount, RealType(0));
        if (eBlocks > 0)
        {
            computeDivergencePerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7,
                s.d_u.data(), s.d_v.data(), s.d_w.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                d_divAccPost.data(), startElem, numLocal);
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(d_divAccPost);
        maybePeriodicSum<KeyType, RealType>(s, d_divAccPost);
        normalizeDivergencePerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_divAccPost.data(), s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            d_divNorm.data(), s.nodeCount);
        cudaDeviceSynchronize();

        // Take max |div| over interior owned DOFs (skip boundary -- the boundary
        // divergence picks up the non-zero through-flux of u_lid and is not
        // expected to be zero pointwise). Also skip pressure-Dirichlet DOFs
        // (channel outflow mask): the matrix-free SpMV's identity-row at those
        // slots breaks the local D u = 0 identity by design, so leaving them
        // in the metric reports a meaningless cross-boundary flux that does
        // NOT indicate the interior projection is wrong.
        const uint8_t* ownPtr = d_nodeOwnership.data();
        const int* dofPtr     = s.d_node_to_dof.data();
        const uint8_t* bnd    = s.d_isBdryDof.data();
        const uint8_t* pmask  = (s.d_isPressureBdryDof.size() > 0)
                                ? s.d_isPressureBdryDof.data() : nullptr;
        const RealType* dPtr  = d_divNorm.data();
        RealType localMax = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [ownPtr, dofPtr, bnd, pmask, dPtr] __device__ (size_t i) -> RealType {
                if (ownPtr[i] != 1) return RealType(0);
                int dof = dofPtr[i];
                if (dof < 0)                       return RealType(0);
                if (bnd[dof])                      return RealType(0);
                if (pmask != nullptr && pmask[dof]) return RealType(0);
                return fabs(dPtr[i]);
            }, RealType(0), thrust::maximum<RealType>());

        auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
        RealType globalMax = 0;
        MPI_Allreduce(&localMax, &globalMax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
        s.lastDivMax = globalMax;
        // Bulk-interior RMS divergence: dominated by the body of the domain,
        // not the 1-cell ring next to Dirichlet pressure outflow that dominates
        // max. Channel-flow with a wide pressure mask shows RMS << max if the
        // projection is doing its job.
        s.lastDivRms = rmsOwnedInterior1<KeyType, RealType>(s, d_divNorm);
    }
}

// =============================================================================
// Top-level: chain the four sub-steps. nu is captured by the implicit matrix
// (assembled once at setup time), not used here directly -- but the parameter
// remains for API stability with the original signature.
// =============================================================================
// Max absolute value of an owned-portion device vector (skip NaN; report
// Inf-equivalent if found). Diagnostic only.
template<typename RealType>
inline RealType maxAbsOwned(const cstone::DeviceVector<RealType>& v, size_t n)
{
    if (n == 0) return RealType(0);
    auto bp = thrust::device_pointer_cast(v.data());
    auto mm = thrust::transform_reduce(
        thrust::device, bp, bp + n,
        [] __device__ (RealType x) -> RealType {
            return isnan(x) ? RealType(-1) : fabs(x);
        },
        RealType(0), thrust::maximum<RealType>());
    return mm;
}

// Kinetic energy over owned nodes weighted by lumped mass. Consistent
// with the projection step's L^2 inner product on the velocity space.
// Returns 0.5 * sum_(i in owned) m_i * (u_i^2 + v_i^2 + w_i^2).
template<typename KeyType, typename RealType>
inline RealType keOwned(NSStepper<KeyType, RealType>& s,
                        const cstone::DeviceVector<RealType>& du,
                        const cstone::DeviceVector<RealType>& dv,
                        const cstone::DeviceVector<RealType>& dw)
{
    const auto& d_own = s.domain.getNodeOwnershipMap();
    size_t n = s.nodeCount;
    RealType local = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(n),
        [own = d_own.data(), n2d = s.d_node_to_dof.data(),
         mass = s.d_mass.data(), nOwn = s.numOwnedDofs,
         u = du.data(), v = dv.data(), w = dw.data()] __device__ (size_t i) -> RealType {
            if (own[i] != 1) return RealType(0);
            int dof = n2d[i];
            if (dof < 0 || dof >= nOwn) return RealType(0);
            RealType uu = u[i], vv = v[i], ww = w[i];
            return RealType(0.5) * mass[dof] * (uu * uu + vv * vv + ww * ww);
        }, RealType(0), thrust::plus<RealType>());
    RealType global = 0;
    MPI_Datatype mpi_r = std::is_same<RealType, double>::value ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&local, &global, 1, mpi_r, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

// Sum of an owned-node device vector (e.g. divergence accumulator). Used to
// check the consistency condition int(div u) = 0 for periodic NS.
template<typename KeyType, typename RealType>
inline RealType sumOwned(NSStepper<KeyType, RealType>& s,
                         const cstone::DeviceVector<RealType>& v)
{
    const auto& d_own = s.domain.getNodeOwnershipMap();
    size_t n = s.nodeCount;
    RealType local = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(n),
        [own = d_own.data(), n2d = s.d_node_to_dof.data(),
         nOwn = s.numOwnedDofs, vp = v.data()] __device__ (size_t i) -> RealType {
            if (own[i] != 1) return RealType(0);
            int dof = n2d[i];
            if (dof < 0 || dof >= nOwn) return RealType(0);
            return vp[i];
        }, RealType(0), thrust::plus<RealType>());
    RealType global = 0;
    MPI_Datatype mpi_r = std::is_same<RealType, double>::value ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&local, &global, 1, mpi_r, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

// Compute div(u, v, w) lumped per node, return (max_abs, RMS) over owned
// interior nodes. Used to inspect divergence of intermediate velocity fields.
template<typename KeyType, typename RealType>
inline void divMaxAndRmsOwned(NSStepper<KeyType, RealType>& s,
                              const cstone::DeviceVector<RealType>& d_u,
                              const cstone::DeviceVector<RealType>& d_v,
                              const cstone::DeviceVector<RealType>& d_w,
                              RealType& outMax, RealType& outRms)
{
    const auto& d_own = s.domain.getNodeOwnershipMap();
    const auto& d_conn = s.domain.getElementToNodeConnectivity();
    auto c0 = std::get<0>(d_conn).data(); auto c1 = std::get<1>(d_conn).data();
    auto c2 = std::get<2>(d_conn).data(); auto c3 = std::get<3>(d_conn).data();
    auto c4 = std::get<4>(d_conn).data(); auto c5 = std::get<5>(d_conn).data();
    auto c6 = std::get<6>(d_conn).data(); auto c7 = std::get<7>(d_conn).data();
    size_t startElem = s.domain.startIndex();
    size_t numLocal  = s.domain.localElementCount();
    int eBlocks = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
    int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;

    cstone::DeviceVector<RealType> d_divAcc(s.nodeCount, RealType(0));
    if (eBlocks > 0) {
        computeDivergencePerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            c0, c1, c2, c3, c4, c5, c6, c7,
            d_u.data(), d_v.data(), d_w.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            d_divAcc.data(), startElem, numLocal);
        cudaDeviceSynchronize();
    }
    s.domain.reverseExchangeNodeHaloAdd(d_divAcc);
    cstone::DeviceVector<RealType> d_divNorm(s.nodeCount, RealType(0));
    normalizeDivergencePerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        d_divAcc.data(), s.d_mass.data(),
        s.d_node_to_dof.data(), d_own.data(),
        d_divNorm.data(), s.nodeCount);
    cudaDeviceSynchronize();
    outMax = maxOwnedInteriorAbs<KeyType, RealType>(s, d_divNorm);
    outRms = rmsOwnedInterior1<KeyType, RealType>(s, d_divNorm);
}

// =============================================================================
// Per-node vorticity magnitude |omega| = |curl(u)|. Three SCS-gradient passes
// (one per velocity component) give the 9-component velocity gradient tensor;
// the curl is then assembled per node. Uses the same scatter -> reverse-halo
// -> periodic-pair-sum -> normalize-by-mass pipeline as the predictor's
// grad(p), so it inherits the periodic and ghost handling automatically.
// d_omegaMag is resized to nodeCount and overwritten; owned-node values are
// the true vorticity magnitude; ghost slots are zeroed by the normalize step.
// =============================================================================
template<typename KeyType, typename RealType>
void computeVorticityMagnitudePerNode(NSStepper<KeyType, RealType>& s,
                                      cstone::DeviceVector<RealType>& d_u,
                                      cstone::DeviceVector<RealType>& d_v,
                                      cstone::DeviceVector<RealType>& d_w,
                                      cstone::DeviceVector<RealType>& d_omegaMag)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    auto c0 = std::get<0>(d_conn).data(); auto c1 = std::get<1>(d_conn).data();
    auto c2 = std::get<2>(d_conn).data(); auto c3 = std::get<3>(d_conn).data();
    auto c4 = std::get<4>(d_conn).data(); auto c5 = std::get<5>(d_conn).data();
    auto c6 = std::get<6>(d_conn).data(); auto c7 = std::get<7>(d_conn).data();
    const size_t startElem = s.domain.startIndex();
    const size_t numLocal  = s.domain.localElementCount();
    const int nodeBlocks   = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    const int eBlocks      = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;

    if (d_omegaMag.size() != s.nodeCount) d_omegaMag.resize(s.nodeCount);

    auto gradScalar = [&] (const RealType* qNode,
                           cstone::DeviceVector<RealType>& gx,
                           cstone::DeviceVector<RealType>& gy,
                           cstone::DeviceVector<RealType>& gz)
    {
        cstone::DeviceVector<RealType> gxAcc(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> gyAcc(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> gzAcc(s.nodeCount, RealType(0));
        if (eBlocks > 0)
        {
            computeGradientPerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7, qNode,
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                gxAcc.data(), gyAcc.data(), gzAcc.data(),
                startElem, numLocal);
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(gxAcc);
        s.domain.reverseExchangeNodeHaloAdd(gyAcc);
        s.domain.reverseExchangeNodeHaloAdd(gzAcc);
        maybePeriodicSum<KeyType, RealType>(s, gxAcc);
        maybePeriodicSum<KeyType, RealType>(s, gyAcc);
        maybePeriodicSum<KeyType, RealType>(s, gzAcc);
        normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            gxAcc.data(), gyAcc.data(), gzAcc.data(),
            s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            gx.data(), gy.data(), gz.data(),
            s.nodeCount);
        cudaDeviceSynchronize();
    };

    // Ghosts of u/v/w must be current before SCS scatters read q_f = 1/2(q[L]+q[R]).
    s.domain.exchangeNodeHalo(d_u);
    s.domain.exchangeNodeHalo(d_v);
    s.domain.exchangeNodeHalo(d_w);

    cstone::DeviceVector<RealType> dudx(s.nodeCount), dudy(s.nodeCount), dudz(s.nodeCount);
    cstone::DeviceVector<RealType> dvdx(s.nodeCount), dvdy(s.nodeCount), dvdz(s.nodeCount);
    cstone::DeviceVector<RealType> dwdx(s.nodeCount), dwdy(s.nodeCount), dwdz(s.nodeCount);
    gradScalar(d_u.data(), dudx, dudy, dudz);
    gradScalar(d_v.data(), dvdx, dvdy, dvdz);
    gradScalar(d_w.data(), dwdx, dwdy, dwdz);

    // omega = (dwdy - dvdz, dudz - dwdx, dvdx - dudy); |omega| per node.
    const RealType* aP = dudy.data(); const RealType* bP = dudz.data();
    const RealType* cP = dvdx.data(); const RealType* dP = dvdz.data();
    const RealType* eP = dwdx.data(); const RealType* fP = dwdy.data();
    RealType* outP = d_omegaMag.data();
    thrust::for_each(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount),
        [aP, bP, cP, dP, eP, fP, outP] __device__ (size_t i) {
            RealType wx = fP[i] - dP[i];
            RealType wy = bP[i] - eP[i];
            RealType wz = cP[i] - aP[i];
            outP[i] = sqrt(wx * wx + wy * wy + wz * wz);
        });
    cudaDeviceSynchronize();
}

// Set to >0 to dump per-phase diagnostics for the first N calls.
static int g_nsDebugStepsLeft = 0;

template<typename KeyType, typename RealType>
void runNsStep(NSStepper<KeyType, RealType>& s, RealType dt, RealType nu, RealType rho)
{
    s.nuCached = nu;
    bool dbg = (g_nsDebugStepsLeft > 0);

    // ENTRY: state of u^n, p^n before the step.
    if (dbg && s.rank == 0)
    {
        RealType KE_n      = keOwned<KeyType, RealType>(s, s.d_u, s.d_v, s.d_w);
        RealType div_n_max = 0, div_n_rms = 0;
        divMaxAndRmsOwned<KeyType, RealType>(s, s.d_u, s.d_v, s.d_w, div_n_max, div_n_rms);
        std::cout << "    [ns-dbg ENTRY] KE_n=" << std::scientific << std::setprecision(6) << KE_n
                  << " |u_n|max=" << maxAbsOwned(s.d_u, s.nodeCount)
                  << " |p_n|max=" << maxAbsOwned(s.d_p, s.nodeCount)
                  << " div(u_n)max=" << div_n_max
                  << " div(u_n)rms=" << div_n_rms << "\n";
    }

    runPredictorStep<KeyType, RealType>(s, dt, rho);
    if (dbg && s.rank == 0)
    {
        RealType KE_star      = keOwned<KeyType, RealType>(s, s.d_uStar, s.d_vStar, s.d_wStar);
        RealType div_star_max = 0, div_star_rms = 0;
        divMaxAndRmsOwned<KeyType, RealType>(s, s.d_uStar, s.d_vStar, s.d_wStar, div_star_max, div_star_rms);
        std::cout << "    [ns-dbg PRED ] KE*=" << KE_star
                  << " |u*|max="    << maxAbsOwned(s.d_uStar, s.nodeCount)
                  << " |gP|max="    << maxAbsOwned(s.d_gradPx, s.nodeCount)
                  << " div(u*)max=" << div_star_max
                  << " div(u*)rms=" << div_star_rms << "\n";
    }

    runImplicitDiffusionStep<KeyType, RealType>(s, dt);
    if (dbg && s.rank == 0)
    {
        RealType KE_ss      = keOwned<KeyType, RealType>(s, s.d_uStarStar, s.d_vStarStar, s.d_wStarStar);
        RealType div_ss_max = 0, div_ss_rms = 0;
        divMaxAndRmsOwned<KeyType, RealType>(s, s.d_uStarStar, s.d_vStarStar, s.d_wStarStar, div_ss_max, div_ss_rms);
        std::cout << "    [ns-dbg DIFF ] KE**=" << KE_ss
                  << " |u**|max="    << maxAbsOwned(s.d_uStarStar, s.nodeCount)
                  << " div(u**)max=" << div_ss_max
                  << " div(u**)rms=" << div_ss_rms
                  << " cg_uvw="      << s.lastUIters << "/" << s.lastVIters << "/" << s.lastWIters << "\n";
    }

    runPressureSolveStep<KeyType, RealType>(s, dt, rho);
    if (dbg && s.rank == 0)
    {
        std::cout << "    [ns-dbg POIS ] |phi|max=" << maxAbsOwned(s.d_phi, s.nodeCount)
                  << " cg_p=" << s.lastPressureIters << "\n";
    }

    // BDF2 velocity-history shuffle: snapshot u^n into u_{n-1} BEFORE the
    // corrector overwrites s.d_u with u^{n+1}. The advection-history copy
    // (advN -> advNm1) happens after the corrector since the advection slots
    // are only consumed by next step's predictor.
    if (s.useBdf2
        && s.d_u_nm1.size() == s.d_u.size()
        && s.d_v_nm1.size() == s.d_v.size()
        && s.d_w_nm1.size() == s.d_w.size())
    {
        thrust::copy(thrust::device, s.d_u.begin(), s.d_u.end(), s.d_u_nm1.begin());
        thrust::copy(thrust::device, s.d_v.begin(), s.d_v.end(), s.d_v_nm1.begin());
        thrust::copy(thrust::device, s.d_w.begin(), s.d_w.end(), s.d_w_nm1.begin());
    }

    runCorrectorStep<KeyType, RealType>(s, dt, rho);
    if (dbg && s.rank == 0)
    {
        RealType KE_np1 = keOwned<KeyType, RealType>(s, s.d_u, s.d_v, s.d_w);
        std::cout << "    [ns-dbg CORR ] KE_(n+1)=" << KE_np1
                  << " |u_(n+1)|max=" << maxAbsOwned(s.d_u, s.nodeCount)
                  << " |p_(n+1)|max=" << maxAbsOwned(s.d_p, s.nodeCount)
                  << " div_max=" << s.lastDivMax << "\n";
        std::cout << std::defaultfloat;
    }

    // BDF2 advection-history shuffle + bootstrap promotion.
    //   step 0  (BDF1): produces N(u^0) in advN. Save -> advNm1 for step 1.
    //   step >=1 (BDF2): produces N(u^n) in advN. Save -> advNm1 for step n+1.
    // After the shuffle advance bdfStep so the NEXT runNsStep call uses BDF2.
    if (s.useBdf2
        && s.d_advU_n.size() == s.nodeCount
        && s.d_advU_nm1.size() == s.nodeCount)
    {
        thrust::copy(thrust::device, s.d_advU_n.begin(), s.d_advU_n.end(), s.d_advU_nm1.begin());
        thrust::copy(thrust::device, s.d_advV_n.begin(), s.d_advV_n.end(), s.d_advV_nm1.begin());
        thrust::copy(thrust::device, s.d_advW_n.begin(), s.d_advW_n.end(), s.d_advW_nm1.begin());
        if (s.bdfStep < 1) s.bdfStep = 1;
    }

    if (g_nsDebugStepsLeft > 0) --g_nsDebugStepsLeft;
}

// =============================================================================
// L2 norm over owned interior+boundary DOFs (volume-weighted by lumped mass).
// Used for per-step |u|, |v|, |w|, |p| logging.
// =============================================================================

template<typename KeyType, typename RealType>
RealType computeWeightedL2Norm(NSStepper<KeyType, RealType>& s,
                               const cstone::DeviceVector<RealType>& d_q)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    cstone::DeviceVector<RealType> d_sq(s.numOwnedDofs, RealType(0));
    thrust::for_each(thrust::device,
                      thrust::counting_iterator<size_t>(0),
                      thrust::counting_iterator<size_t>(s.nodeCount),
                      [nodeToDof = s.d_node_to_dof.data(),
                       ownPtr    = d_nodeOwnership.data(),
                       mass      = s.d_mass.data(),
                       q         = d_q.data(),
                       out       = d_sq.data()] __device__(size_t i)
                      {
                          if (ownPtr[i] != 1) return;
                          int dof = nodeToDof[i];
                          if (dof < 0) return;
                          out[dof] = q[i] * q[i] * mass[dof];
                      });
    cudaDeviceSynchronize();
    auto sp = thrust::device_pointer_cast(d_sq.data());
    RealType localSum = thrust::reduce(thrust::device, sp, sp + s.numOwnedDofs, RealType(0));
    RealType globalSum = 0;
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&localSum, &globalSum, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    return std::sqrt(globalSum);
}
