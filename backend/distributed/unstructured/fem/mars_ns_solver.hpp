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
#include "backend/distributed/unstructured/fem/mars_element_traits.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_tet_area.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_tet_assembler.hpp"
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
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>
#include <thrust/extrema.h>
#include <thrust/fill.h>
#include <thrust/system/cuda/execution_policy.h>

#include <memory>
#include <array>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <string>
#include <cstdlib>
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
                                   RealType lidU,
                                   int numOwnedDofs)
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
    if (dof < 0 || dof >= numOwnedDofs) return;
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
                                    RealType Uinf,
                                    int numOwnedDofs)
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
    if (dof < 0 || dof >= numOwnedDofs) return;
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
                                            RealType eps,
                                            int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;
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

// Pin enforcement that PRESERVES the diagonal magnitude: zero the row's
// off-diagonals, keep A[pin,pin] at its assembled value (do NOT force it to
// 1). The DDT operator has tiny diagonals (~0.04-0.14); forcing the pin row
// to 1.0 makes it 10-25x its neighbors -- CG tolerates that spike, but
// BoomerAMG's strength-of-connection coarsening degrades and Hypre setup
// returns generic error 1. Keeping the native diagonal removes the spike so
// AMG sees a uniformly-scaled SPD row. Pin value = 0, so RHS needs no lift.
template<typename RealType>
__global__ void enforcePinRowKeepDiagKernel(int pinDof,
                                            const int* rowPtr,
                                            const int* colInd,
                                            const int* diagPtr,
                                            RealType* values)
{
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t != 0) return;
    if (pinDof < 0) return;
    int dp = diagPtr[pinDof];
    for (int j = rowPtr[pinDof]; j < rowPtr[pinDof + 1]; ++j)
        if (j != dp) values[j] = RealType(0);
    // Leave values[dp] as assembled. If the diagonal is somehow absent or
    // non-positive, fall back to 1 so the row stays well-posed.
    if (dp < 0) return;
    if (!(values[dp] > RealType(0))) values[dp] = RealType(1);
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

// Multi-DOF symmetric Dirichlet (channel outflow, pump outlet). For every
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

// Path-B helper: mark owned slave DOFs of cross-rank periodic pairs in a
// per-owned-DOF boolean mask. One thread per send-list entry. Writes are
// race-free because each slave node maps to a distinct owned DOF.
__global__ void markPeriodicXRSlaveDofKernel(const int* d_sendOwnedSlaveNodeIds,
                                             const int* d_nodeToDof,
                                             int numSendIds,
                                             int numOwnedDofs,
                                             uint8_t* d_mask)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numSendIds) return;
    int node = d_sendOwnedSlaveNodeIds[i];
    int dof  = d_nodeToDof[node];
    if (dof >= 0 && dof < numOwnedDofs) d_mask[dof] = 1;
}

// Path-B helper: zero the implicit-diffusion RHS at cross-rank slave DOFs so
// the Dirichlet-identity slave row solves trivially to x[slave] = 0. The
// post-solve crossRankPeriodicBroadcastDof then overwrites x[slave] with
// x[master_owned_D] via MPI.
template<typename RealType>
__global__ void zeroDofAtCrossRankSlavesKernel(const int* d_sendOwnedSlaveNodeIds,
                                               const int* d_nodeToDof,
                                               int numSendIds,
                                               int numOwnedDofs,
                                               RealType* b)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numSendIds) return;
    int node = d_sendOwnedSlaveNodeIds[i];
    int dof  = d_nodeToDof[node];
    if (dof >= 0 && dof < numOwnedDofs) b[dof] = RealType(0);
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
                                             size_t numNodes,
                                             int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;
    if (isBoundaryDof[dof]) rhs[dof] = qTarget[i];
}

// =============================================================================
// Lumped mass: per-node + reverse-halo pattern (same as B.2/B.3/B.4). Single
// element-volume scatter feeds both the diffusion matrix (added to its diag as
// M/dt) and the per-step predictor/corrector normalization. V_node matches the
// per-node scatter used by gradient/divergence so the projection step is
// algebraically consistent across operators.
// =============================================================================

// Element-tag dispatch for the SCS (sub-control-surface) left/right node indices.
// Hex uses d_hexLRSCV[24] (mars_cvfem_utils.hpp), tet uses d_tetLRSCV[12]
// (mars_element_traits.hpp). Both are __constant__ int arrays laid out [L,R,...].
template<typename ElementTag>
__device__ __forceinline__ void scsLR(int ip, int& L, int& R)
{
    if constexpr (std::is_same_v<ElementTag, HexTag>) { L = d_hexLRSCV[ip * 2]; R = d_hexLRSCV[ip * 2 + 1]; }
    else { L = d_tetLRSCV[ip * 2]; R = d_tetLRSCV[ip * 2 + 1]; }
}

// Host-side connectivity pointer gather. Returns NodesPerElem device pointers
// from the connectivity tuple so launch sites can dispatch on element type
// without spelling std::get<7> for tet (which has only 4 columns).
template<typename ElementTag, typename KeyType, typename ConnTuple>
inline std::array<const KeyType*, ElemTraits<ElementTag>::NodesPerElem> connPtrs(const ConnTuple& d_conn)
{
    std::array<const KeyType*, ElemTraits<ElementTag>::NodesPerElem> p{};
    if constexpr (std::is_same_v<ElementTag, HexTag>) {
        p = { std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
              std::get<4>(d_conn).data(), std::get<5>(d_conn).data(), std::get<6>(d_conn).data(), std::get<7>(d_conn).data() };
    } else {
        p = { std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(), std::get<3>(d_conn).data() };
    }
    return p;
}

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
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
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    RealType x[NPE], y[NPE], z[NPE];
    for (int i = 0; i < NPE; ++i)
    {
        x[i] = nodeX[n[i]];
        y[i] = nodeY[n[i]];
        z[i] = nodeZ[n[i]];
    }
    // Scatter to all corners (owned + ghost). The dispatcher passes only the
    // OWNED element range; ghost-slot contributions are then folded back to the
    // true owner by reverseExchangeNodeHaloAdd. (For cross-rank periodic pairs,
    // slave and master have DISTINCT SFC keys -- reverse-halo can't bridge them
    // automatically; a periodic-pair MPI exchange is required after this step.)
    if constexpr (std::is_same_v<ElementTag, HexTag>)
    {
        RealType xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0], zmin = z[0], zmax = z[0];
        for (int i = 1; i < 8; ++i)
        {
            if (x[i] < xmin) xmin = x[i]; if (x[i] > xmax) xmax = x[i];
            if (y[i] < ymin) ymin = y[i]; if (y[i] > ymax) ymax = y[i];
            if (z[i] < zmin) zmin = z[i]; if (z[i] > zmax) zmax = z[i];
        }
        RealType contrib = (xmax - xmin) * (ymax - ymin) * (zmax - zmin) * RealType(0.125);
        for (int i = 0; i < 8; ++i) atomicAdd(&massNode[n[i]], contrib);
    }
    else
    {
        // Tet: true cell volume V = |det[e1 e2 e3]| / 6, split equally to 4 nodes.
        RealType e1[3] = {x[1] - x[0], y[1] - y[0], z[1] - z[0]};
        RealType e2[3] = {x[2] - x[0], y[2] - y[0], z[2] - z[0]};
        RealType e3[3] = {x[3] - x[0], y[3] - y[0], z[3] - z[0]};
        RealType det = e1[0] * (e2[1] * e3[2] - e2[2] * e3[1])
                     - e1[1] * (e2[0] * e3[2] - e2[2] * e3[0])
                     + e1[2] * (e2[0] * e3[1] - e2[1] * e3[0]);
        RealType V = fabs(det) / RealType(6);
        RealType contrib = V / RealType(4);
        for (int i = 0; i < 4; ++i) atomicAdd(&massNode[n[i]], contrib);
    }
}

// Owned-node mass slots -> per-DOF mass array (post reverse halo). atomicAdd
// because multiple nodes may map to the same DOF under periodic collapse
// (slave + master both index the master DOF); plain assignment would race.
template<typename RealType>
__global__ void gatherOwnedNodeMassToDofKernel(const RealType* massNode,
                                               const int* nodeToDof,
                                               const uint8_t* ownership,
                                               RealType* massDof,
                                               size_t numNodes,
                                               int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;
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

// O(1)-per-row diagonal extraction using the precomputed d_diagPtr indirection.
// Cached at setup once (after BC enforcement + optional MARS_DDT_DIAG_SHIFT)
// and reused as the Jacobi preconditioner in solvePressureDDT every CG iter.
// Floor by 1e-30 so the 1/diag step in jacobiPrecondKernel cannot divide by
// zero on a (pathological) empty/zero-diag row.
template<typename RealType>
__global__ void extractDiagByDiagPtrKernel(const int* diagPtr,
                                           const RealType* values,
                                           RealType* diag,
                                           int numOwnedDofs)
{
    int dof = blockIdx.x * blockDim.x + threadIdx.x;
    if (dof >= numOwnedDofs) return;
    int dp = diagPtr[dof];
    RealType d = (dp >= 0) ? values[dp] : RealType(0);
    // Defensive: a non-positive diagonal would break Jacobi. Floor to keep CG
    // robust without changing the SpMV (we only divide r by this).
    if (!(d > RealType(1e-30))) d = RealType(1);
    diag[dof] = d;
}

// Matrix-free Jacobi diagonal for the tet D M^-1 D^T pressure operator.
// The assembled DDT matrix (and extractDiagByDiagPtrKernel above) is hex-only,
// so tet has no diagonal to precondition with -> CG stalls. This builds the
// diagonal directly from the same SCS area vectors and node lumped mass the
// matrix-free applyDDTPerNode uses.
//
// Derivation (verified against the 3 operator steps): set phi = e_i and keep
// each face's own self-term. For SCS face f=(L,R) with area A_f, applyDivT
// scatters dp=0.5*(p_L-p_R) * A_f to BOTH endpoints (same sign), normalize
// scales by 1/M, and computeDivergence forms 0.5*(g_L+g_R).A_f. The result:
// each incident face contributes  +0.25 * |A_f|^2 * (1/M_L + 1/M_R)  to the
// diagonal at i (sign^2 = 1, so strictly positive -> SPD-safe Jacobi). The
// per-face value is symmetric in L,R, so both endpoints receive the same
// contribution. Cross-face coupling (two distinct incident faces) is dropped;
// harmless for a preconditioner (only needs to be a good positive scaling).
//
// Mass MUST come from massNode (NODE-indexed, halo-filled) not the DOF-indexed
// d_mass: the neighbor endpoint can be a ghost with no DOF slot.
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
__global__ void computeTetDDTDiagonalKernel(const KeyType* c0, const KeyType* c1,
                                            const KeyType* c2, const KeyType* c3,
                                            const KeyType* c4, const KeyType* c5,
                                            const KeyType* c6, const KeyType* c7,
                                            const RealType* areaVecX,
                                            const RealType* areaVecY,
                                            const RealType* areaVecZ,
                                            const RealType* massNode,
                                            RealType* diagAccNode,
                                            size_t startElem,
                                            size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    constexpr int NSCS = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    #pragma unroll
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int nodeL, nodeR; scsLR<ElementTag>(ip, nodeL, nodeR);
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];

        size_t off = e * NSCS + ip;
        RealType ax = areaVecX[off], ay = areaVecY[off], az = areaVecZ[off];
        RealType A2 = ax*ax + ay*ay + az*az;

        RealType mL = massNode[iL], mR = massNode[iR];
        RealType invSum = RealType(0);
        if (mL > RealType(0)) invSum += RealType(1) / mL;
        if (mR > RealType(0)) invSum += RealType(1) / mR;
        RealType c = RealType(0.25) * A2 * invSum;

        atomicAdd(&diagAccNode[iL], c);
        atomicAdd(&diagAccNode[iR], c);
    }
}

// Fold owned-node diagonal accumulators into the per-DOF diagonal, flooring
// non-positive entries (mirrors extractDiagByDiagPtrKernel's robustness).
template<typename RealType>
__global__ void gatherTetDDTDiagToDofKernel(const RealType* diagAccNode,
                                            const int* nodeToDof,
                                            const uint8_t* ownership,
                                            RealType* diagDof,
                                            size_t numNodes,
                                            int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;
    atomicAdd(&diagDof[dof], diagAccNode[i]);
}

// Clipped Jacobi preconditioner: z[i] = r[i] / max(diag[dof(i)], epsClip)
// for owned nodes, 0 elsewhere. The clip floor is essential on wing-type
// meshes: high-aspect-ratio boundary-layer cells produce DDT diagonals
// down to ~4.6e-5 while interior diag ~12.3, a 5-orders-of-magnitude spread.
// Pure 1/diag would amplify residuals at the tiny-diag rows by ~250000x
// every iter, biasing the search direction and breaking pAp>0. The clip
// caps amplification at 1/epsClip while leaving the well-scaled bulk of
// the operator preconditioned correctly. epsClip=0 disables the clip
// (recovers the vanilla Jacobi).
template<typename RealType>
__global__ void jacobiPrecondNodeKernel(const RealType* r,
                                        const RealType* diag,
                                        const int* nodeToDof,
                                        const uint8_t* ownership,
                                        RealType epsClip,
                                        RealType* z,
                                        size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) { z[i] = RealType(0); return; }
    int dof = nodeToDof[i];
    if (dof < 0) { z[i] = RealType(0); return; }
    RealType d = diag[dof];
    if (d < epsClip) d = epsClip;
    z[i] = r[i] / d;
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
//   MARS_NS_ADV_BJ (2nd-order Barth-Jespersen limited upwind; tet-only):
//     Like upwind, but the upwind node value is linearly reconstructed to the
//     face with a limited node gradient:
//       up      = (mdot>0) ? L : R
//       q_face  = q[up] + phi[up] * grad_q[up] . (M - x_up)
//     M is the SCS-edge midpoint, so (M - x_up) = +-0.5*(xR - xL). phi in [0,1]
//     is the Barth-Jespersen limiter that keeps q_face monotone (no new
//     extrema vs the edge-neighbor min/max). phi=0 collapses to 1st-order
//     upwind; phi=1 is the unlimited 2nd-order reconstruction. This matches
//     the mesh developers' legacy advection. Gradient and phi are computed in
//     runPredictorStep and halo-exchanged before this kernel runs.
//
// advMode selects the form at runtime: 0=skew, 1=upwind, 2=Barth-Jespersen.
// Sign convention matches the legacy upwind path: net flux out of L, into R.
// =============================================================================

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
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
                                                          int advMode,
                                                          // Barth-Jespersen inputs (advMode==2 only). For the convected
                                                          // component q: limiter phi[node] in [0,1] and node gradient
                                                          // grad_q. nodeX/Y/Z give the edge midpoint for reconstruction.
                                                          // All nullptr / unread for skew and plain upwind.
                                                          const RealType* phiQ,
                                                          const RealType* gradQx,
                                                          const RealType* gradQy,
                                                          const RealType* gradQz,
                                                          const RealType* nodeX,
                                                          const RealType* nodeY,
                                                          const RealType* nodeZ)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    constexpr int NSCS = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    #pragma unroll
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int nodeL, nodeR; scsLR<ElementTag>(ip, nodeL, nodeR);
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];

        RealType vfx = RealType(0.5) * (vx[iL] + vx[iR]);
        RealType vfy = RealType(0.5) * (vy[iL] + vy[iR]);
        RealType vfz = RealType(0.5) * (vz[iL] + vz[iR]);

        size_t off = e * NSCS + ip;
        RealType mdot = vfx * areaVecX[off] + vfy * areaVecY[off] + vfz * areaVecZ[off];

        if (advMode == 0)
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
        else if (advMode == 1)
        {
            // 1st-order upwind (standard CVFEM, used by cavity/channel).
            RealType q_face = (mdot > RealType(0)) ? q[iL] : q[iR];
            RealType flux   = mdot * q_face;
            atomicAdd(&dqdtNode[iL], -flux);
            atomicAdd(&dqdtNode[iR], +flux);
        }
        else
        {
            // 2nd-order Barth-Jespersen limited upwind.
            // Reconstruct the upwind node value to the SCS-edge midpoint M and
            // limit the slope with the precomputed per-node phi. Midpoint
            // reconstruction (r = M - x_up = +-0.5*(xR-xL)) is a defensible
            // 2nd-order CVFEM target that needs only node coords; the full
            // dual-quad-centroid IP is a later refinement.
            KeyType up = (mdot > RealType(0)) ? iL : iR;
            RealType mx = RealType(0.5) * (nodeX[iL] + nodeX[iR]);
            RealType my = RealType(0.5) * (nodeY[iL] + nodeY[iR]);
            RealType mz = RealType(0.5) * (nodeZ[iL] + nodeZ[iR]);
            RealType rx = mx - nodeX[up];
            RealType ry = my - nodeY[up];
            RealType rz = mz - nodeZ[up];
            RealType recon = gradQx[up] * rx + gradQy[up] * ry + gradQz[up] * rz;
            RealType q_face = q[up] + phiQ[up] * recon;
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

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
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
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    constexpr int NSCS = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    #pragma unroll
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int nodeL, nodeR; scsLR<ElementTag>(ip, nodeL, nodeR);
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];
        RealType pf = RealType(0.5) * (p[iL] + p[iR]);

        size_t off = e * NSCS + ip;
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

// =============================================================================
// Barth-Jespersen limiter support kernels (used only when advScheme==
// BarthJespersen). Defined here, ABOVE runPredictorStep, so the predictor can
// launch them without a forward reference.
//
// CUDA has no native atomicMin/Max on double, so we do the standard atomicCAS
// loop on the 64-bit pattern. Monotone bit order of IEEE-754 doubles is NOT
// preserved across sign, so we compare the decoded double each iteration rather
// than comparing the raw bits. This is a limiter (not exactness-critical), but
// the CAS loop is exact anyway.
// =============================================================================
template<typename RealType>
__device__ __forceinline__ void atomicMinDouble(RealType* addr, RealType val)
{
    unsigned long long* a = reinterpret_cast<unsigned long long*>(addr);
    unsigned long long  old = *a, assumed;
    do {
        double cur = __longlong_as_double(static_cast<long long>(old));
        if (cur <= double(val)) break;               // already <= val, done
        assumed = old;
        old = atomicCAS(a, assumed,
                        static_cast<unsigned long long>(__double_as_longlong(double(val))));
    } while (assumed != old);
}

template<typename RealType>
__device__ __forceinline__ void atomicMaxDouble(RealType* addr, RealType val)
{
    unsigned long long* a = reinterpret_cast<unsigned long long*>(addr);
    unsigned long long  old = *a, assumed;
    do {
        double cur = __longlong_as_double(static_cast<long long>(old));
        if (cur >= double(val)) break;               // already >= val, done
        assumed = old;
        old = atomicCAS(a, assumed,
                        static_cast<unsigned long long>(__double_as_longlong(double(val))));
    } while (assumed != old);
}

// Seed qmin = qmax = q at every node before the neighbor min/max scatter, so a
// node with no scattered neighbor (shouldn't happen on a valid mesh, but is
// safe) still bounds itself.
template<typename RealType>
__global__ void bjSeedMinMaxKernel(const RealType* q, RealType* qmin, RealType* qmax, size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    qmin[i] = q[i];
    qmax[i] = q[i];
}

// Neighbor min/max over edge-connected nodes. For each SCS edge (iL,iR) push
// the neighbor's value into this node's running min/max. After this scatter
// qmin[i]/qmax[i] hold min/max of q over {i} U {edge neighbors of i} that are
// reachable from this rank's OWNED elements.
//
// Halo: run over OWNED elements only; an owned boundary node's off-rank
// neighbors are completed by reverseExchangeNodeHaloMin/Max in the caller (a
// MIN/MAX fold -- the SUM-based reverse fold would corrupt min/max), then a
// forward exchangeNodeHalo publishes the complete owner bounds to ghosts.
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
__global__ void bjNeighborMinMaxScatterKernel(const KeyType* c0, const KeyType* c1,
                                              const KeyType* c2, const KeyType* c3,
                                              const KeyType* c4, const KeyType* c5,
                                              const KeyType* c6, const KeyType* c7,
                                              const RealType* q,
                                              RealType* qmin, RealType* qmax,
                                              size_t startElem, size_t numElems)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numElems) return;
    size_t e = startElem + k;
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    constexpr int NSCS = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    #pragma unroll
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int nodeL, nodeR; scsLR<ElementTag>(ip, nodeL, nodeR);
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];
        RealType qL = q[iL];
        RealType qR = q[iR];
        atomicMinDouble<RealType>(&qmin[iL], qR);
        atomicMaxDouble<RealType>(&qmax[iL], qR);
        atomicMinDouble<RealType>(&qmin[iR], qL);
        atomicMaxDouble<RealType>(&qmax[iR], qL);
    }
}

// Barth-Jespersen limiter phi per node. For each face incident to a node, the
// reconstruction predicts a face delta Delta = grad_q[i] . (M - x_i). phi_face
// caps the slope so q[i]+phi*Delta does not exceed the neighbor bounds:
//   Delta > eps : phi = min(1, (qmax-q)/Delta)
//   Delta < -eps: phi = min(1, (qmin-q)/Delta)
//   else        : phi = 1
// phi[i] is the MIN of phi_face over all faces incident to i. Seed phi=1
// (bjSeedPhi) then atomicMin each face's phi into both endpoints.
template<typename RealType>
__global__ void bjSeedPhiKernel(RealType* phi, size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    phi[i] = RealType(1);
}

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
__global__ void bjLimiterScatterKernel(const KeyType* c0, const KeyType* c1,
                                       const KeyType* c2, const KeyType* c3,
                                       const KeyType* c4, const KeyType* c5,
                                       const KeyType* c6, const KeyType* c7,
                                       const RealType* q,
                                       const RealType* qmin, const RealType* qmax,
                                       const RealType* gradQx,
                                       const RealType* gradQy,
                                       const RealType* gradQz,
                                       const RealType* nodeX,
                                       const RealType* nodeY,
                                       const RealType* nodeZ,
                                       RealType* phi,
                                       size_t startElem, size_t numElems)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numElems) return;
    size_t e = startElem + k;
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    constexpr int NSCS = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    const RealType eps = RealType(1e-12);

    #pragma unroll
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int nodeL, nodeR; scsLR<ElementTag>(ip, nodeL, nodeR);
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];

        RealType mx = RealType(0.5) * (nodeX[iL] + nodeX[iR]);
        RealType my = RealType(0.5) * (nodeY[iL] + nodeY[iR]);
        RealType mz = RealType(0.5) * (nodeZ[iL] + nodeZ[iR]);

        // phi_face seen from each endpoint of this edge.
        #pragma unroll
        for (int side = 0; side < 2; ++side)
        {
            KeyType i = (side == 0) ? iL : iR;
            RealType rx = mx - nodeX[i];
            RealType ry = my - nodeY[i];
            RealType rz = mz - nodeZ[i];
            RealType delta = gradQx[i] * rx + gradQy[i] * ry + gradQz[i] * rz;
            RealType phi_face;
            if (delta > eps)
                phi_face = fmin(RealType(1), (qmax[i] - q[i]) / delta);
            else if (delta < -eps)
                phi_face = fmin(RealType(1), (qmin[i] - q[i]) / delta);
            else
                phi_face = RealType(1);
            // Clamp to [0,1]: (qmax-q) is >=0 and (qmin-q) is <=0 by
            // construction, so the ratio is >=0; the fmin caps at 1. Guard the
            // floor anyway against roundoff that could make qmax<q by an ulp.
            phi_face = fmax(RealType(0), phi_face);
            atomicMinDouble<RealType>(&phi[i], phi_face);
        }
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
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
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    constexpr int NSCS = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    #pragma unroll
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int nodeL, nodeR; scsLR<ElementTag>(ip, nodeL, nodeR);
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];
        RealType dp = RealType(0.5) * (p[iL] - p[iR]);

        size_t off = e * NSCS + ip;
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
                                               const RealType* lumpedMassNode,  // per-NODE, sized nodeCount (halo-synced)
                                               const int* /*nodeToDof*/,        // kept for ABI; unused
                                               const uint8_t* ownership,
                                               RealType* gx,
                                               RealType* gy,
                                               RealType* gz,
                                               size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) { gx[i] = gy[i] = gz[i] = RealType(0); return; }
    // Index per-NODE mass directly: works for cross-rank periodic slaves
    // whose nodeToDof was redirected to a ghost DOF (out of d_mass bounds).
    RealType m = lumpedMassNode[i];
    if (m == RealType(0)) { gx[i] = gy[i] = gz[i] = RealType(0); return; }
    RealType invV = RealType(1) / m;
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

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
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
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    constexpr int NSCS = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    #pragma unroll
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int nodeL, nodeR; scsLR<ElementTag>(ip, nodeL, nodeR);
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];

        RealType vfx = RealType(0.5) * (vx[iL] + vx[iR]);
        RealType vfy = RealType(0.5) * (vy[iL] + vy[iR]);
        RealType vfz = RealType(0.5) * (vz[iL] + vz[iR]);

        size_t off = e * NSCS + ip;
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

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
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
    constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;
    constexpr int NSCS = ElemTraits<ElementTag>::ScsPerElem;
    const KeyType* cc[8] = {c0, c1, c2, c3, c4, c5, c6, c7};
    KeyType n[NPE];
    for (int i = 0; i < NPE; ++i) n[i] = cc[i][e];

    #pragma unroll
    for (int ip = 0; ip < NSCS; ++ip)
    {
        int nodeL, nodeR; scsLR<ElementTag>(ip, nodeL, nodeR);
        KeyType iL = n[nodeL];
        KeyType iR = n[nodeR];

        // Standard average face velocity.
        RealType vfx = RealType(0.5) * (vx[iL] + vx[iR]);
        RealType vfy = RealType(0.5) * (vy[iL] + vy[iR]);
        RealType vfz = RealType(0.5) * (vz[iL] + vz[iR]);

        size_t off = e * NSCS + ip;
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
                                                 const RealType* lumpedMassNode,  // per-NODE, sized nodeCount
                                                 const int* /*nodeToDof*/,
                                                 const uint8_t* ownership,
                                                 RealType* divU,
                                                 size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) { divU[i] = RealType(0); return; }
    RealType m = lumpedMassNode[i];
    if (m == RealType(0)) { divU[i] = RealType(0); return; }
    divU[i] = divAccNode[i] / m;
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
                                            RealType bodyForce,
                                            size_t numNodes,
                                            int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;

    if (isBdryDof[dof])
    {
        qStar[i] = qTarget[i];
        return;
    }
    RealType V = lumpedMass[dof];
    // BDF1: u* = u^n + dt*(-adv/V - invRho*grad p + invRho*f) (interior only)
    qStar[i] = qN[i] + dt * dqdtNode[i] / V
                    - dt * invRho * gradPnq[i]
                    + dt * invRho * bodyForce;
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
                                                RealType bodyForce,
                                                size_t numNodes,
                                                int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;
    if (isBdryDof[dof])
    {
        qStar[i] = qTarget[i];
        return;
    }
    RealType V       = lumpedMass[dof];
    RealType adv_ext = RealType(2) * dqdtNode_n[i] - dqdtNode_nm1[i];
    RealType coef    = RealType(2) * dt / RealType(3);
    // BDF2 + EXT2 with body-force source: KIO stencil + coef*invRho*f.
    qStar[i] = (RealType(4) / RealType(3)) * qN[i]
             - (RealType(1) / RealType(3)) * qNm1[i]
             + coef * (adv_ext / V - invRho * gradPnq[i] + invRho * bodyForce);
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
                                            size_t numNodes,
                                            int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;

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
                                               size_t numNodes,
                                               int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    // dof can be >= numOwnedDofs when periodic DOF collapse redirected a
    // cross-rank slave onto the master's ghost DOF on this rank. The master
    // rank assembles the matching b[master_owned] from its periodic-image
    // halo elements, so skipping here is correct (no contribution lost).
    if (dof < 0 || dof >= numOwnedDofs) return;
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
                                       size_t numNodes,
                                       int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;
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

template<typename KeyType, typename RealType, typename ElementTag = HexTag> struct NSStepper;

// Forward declaration so the setup-time MARS_DDT_PROBE_DIFF diagnostic can
// call applyDDTPerNode (defined later in this file).
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void applyDDTPerNode(NSStepper<KeyType, RealType, ElementTag>& s,
                     cstone::DeviceVector<RealType>& phi,
                     cstone::DeviceVector<RealType>& outAcc,
                     cstone::DeviceVector<RealType>& gxAcc,
                     cstone::DeviceVector<RealType>& gyAcc,
                     cstone::DeviceVector<RealType>& gzAcc);

template<typename KeyType, typename RealType, typename ElementTag>
struct NSStepper
{
    using DomainT = ElementDomain<ElementTag, RealType, KeyType, cstone::GpuTag>;

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
    // Path-B mask for multi-rank periodic: marks owned slave DOFs whose master
    // lives on a remote rank. enforceBcMatrixKernel zeros those rows and sets
    // diag=1 (Dirichlet identity); the runImplicit RHS-zero step + post-solve
    // crossRankPeriodicBroadcastDof restore x[slave]:=x[master]. Empty on
    // single-rank or non-periodic runs.
    cstone::DeviceVector<uint8_t>  d_isPeriodicXRSlaveDof;
    cstone::DeviceVector<uint8_t>  d_isPeriodicXRSlaveNode;  // per-NODE Path-B mask, halo-exchanged so symmetric column-zero hits ghost copies too
    // Pressure BC: in cavity mode, use a single corner pin (pressurePinDof).
    // In channel mode, use a Dirichlet mask (d_isPressureBdryDof) over the
    // entire outflow face. Exactly one mechanism is active per run.
    int pressurePinDof = -1;   // owned-DOF id on the owning rank; -1 elsewhere (or in channel mode)
    int pressurePinRank = 0;   // global rank that owns the pin
    cstone::DeviceVector<uint8_t> d_isPressureBdryDof;   // size numOwnedDofs; only used in channel mode
    // Per-NODE pressure-Dirichlet mask (size nodeCount, halo-exchanged so
    // ghost slots see the owner's flag). Includes the single corner pin
    // (cavity) AND every outflow-face DOF (channel/pump). Used by symmetric
    // column-zeroing on the assembled DDT CSR so non-owning ranks also zero
    // out their ghost-copy columns; without this, the cavity pin was only
    // enforced symmetrically on rank 0 and the channel mask only on each
    // owner, leaving the global matrix asymmetric and AMG misbehaving.
    cstone::DeviceVector<uint8_t> d_isPressureBdryNode;
    // Per-NODE velocity-Dirichlet mask (size nodeCount, halo-exchanged so ghost
    // slots see the owner's flag). Mirrors d_isPressureBdryNode but driven by
    // the velocity BC predicate (d_isBdryDof). Used by symmetric column-zero
    // on d_valuesVel for Channel/Pump -- without it the velocity matrix is
    // row-cleared but not column-cleared, breaking SPD with many ~125k
    // Dirichlet DOFs + 5-OOM diag spread on anisotropic BL cells flips pAp
    // sign and breaks cstone CG (cg_iter_uvw=FAIL on any non-trivial RHS).
    cstone::DeviceVector<uint8_t> d_isBdryNode;
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
    // Cached diagonal of the BC-modified DDT operator (size numOwnedDofs),
    // used as the Jacobi preconditioner by solvePressureDDT. Extracted once
    // at setup after BC enforcement + optional shift.
    cstone::DeviceVector<RealType> d_diagDDT;
    // Floor for the clipped Jacobi preconditioner: z = r / max(diag, eps).
    // Default is 1e-3 * max(diag), computed at the setup-time cache build.
    // Env MARS_DDT_JACOBI_CLIP_FRAC overrides the 1e-3 fraction.
    RealType diagDDTEpsClip = RealType(0);

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
    // Advection form selector. The flux kernel branches on advScheme:
    //   Skew   = Verstappen central skew-symmetric (discrete KE conservation;
    //            needed for periodic NS where no wall dissipation exists).
    //   Upwind = 1st-order upwind (stable at any Re; cavity/channel default).
    //   BarthJespersen = 2nd-order upwind with a node gradient reconstruction
    //            limited by Barth-Jespersen, matching the mesh developers'
    //            legacy advection. Tet-only path built in this file.
    enum class AdvScheme { Skew, Upwind, BarthJespersen };
    AdvScheme advScheme = AdvScheme::Skew;
    cstone::DeviceVector<RealType> d_gradPx, d_gradPy, d_gradPz;
    cstone::DeviceVector<RealType> d_gradPhix, d_gradPhiy, d_gradPhiz;

    // Barth-Jespersen state (allocated lazily in runPredictorStep only when
    // advScheme==BarthJespersen; non-BJ runs never touch these). Per-component
    // node velocity gradients grad(u), grad(v), grad(w), each (1/V_i) sum_f
    // q_face*A_f, then halo-exchanged. Limiter phi and neighbor min/max are
    // recomputed per component into the scratch buffers below.
    cstone::DeviceVector<RealType> d_gradUx, d_gradUy, d_gradUz;
    cstone::DeviceVector<RealType> d_gradVx, d_gradVy, d_gradVz;
    cstone::DeviceVector<RealType> d_gradWx, d_gradWy, d_gradWz;
    // Per-component scratch reused across the 3 momentum components: neighbor
    // min/max over edge-connected nodes (seed = node's own value) and the
    // Barth-Jespersen limiter phi in [0,1]. Sized nodeCount, halo-exchanged.
    cstone::DeviceVector<RealType> d_bjQmin, d_bjQmax, d_bjPhi;

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
    RealType lastPressR0    = 0;  // initial absolute residual |r0| of the pressure CG this step
    RealType lastPressResid = 0;  // final relative residual |r|/|r0| at the pressure CG's exit
    int lastUIters = 0, lastVIters = 0, lastWIters = 0;
    RealType lastDivMax     = 0;  // |div(u^{n+1})| max -- post-corrector, PLAIN nodal operator
    RealType lastDivRms     = 0;  // |div(u^{n+1})| RMS over interior owned DOFs; less ring-sensitive
    RealType lastDivMaxPre  = 0;  // |div(u**)|    max -- pre-corrector (= b magnitude / V scaled)
    RealType lastDivRC      = 0;  // |div_RC(u^{n+1})| max -- the operator RC actually zeros (only set when useRhieChow)
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
    // Pump: per-side-set Dirichlet driven by Exodus side-sets. walls=no-slip
    //   (u=v=w=0), in=Uinf, out=Dirichlet p=0, other named side-sets get
    //   full Dirichlet u=Uinf (extra Dirichlet walls; true slip BCs are
    //   a follow-up). Side-set local-node lists must be populated before
    //   setupNSStepper(). NSStepper does not allocate or free them.
    // Periodic: triply periodic box. No Dirichlet on velocity. Pressure
    // null-space removed by global-mean subtraction every step. Requires
    // periodicMap to be set before setupNSStepper().
    enum class BCKind { Cavity, Channel, Pump, Periodic };
    BCKind bcKind = BCKind::Cavity;
    RealType Uinf = 1;   // channel/pump inflow speed; reuses lidU value via CLI
    // Inlet velocity direction (unit vector). Default +x for legacy channel/pump.
    // The pump driver sets this to the INWARD normal of the inlet side-set so
    // the prescribed Uinf is applied normal to the opening, not along a global
    // axis. Inlet velocity = Uinf * (inletDirX, inletDirY, inletDirZ).
    RealType inletDirX = 1, inletDirY = 0, inletDirZ = 0;
    // Mass-conserving outlet (pump). When outletU > 0, the Pump outlet nodes are
    // tagged as a velocity-Dirichlet OUTFLOW with velocity outletU along the
    // OUTWARD normal (outletDir), sized so the outlet flux removes the inlet
    // flux -> balances net through-flux so the projection is well-posed (a
    // velocity-inlet + pressure-only outlet leaves an unprojectable mode). The
    // outlet still gets p=0 for the pressure null-space. outletU<=0 keeps the
    // legacy natural-Neumann outflow.
    RealType outletU = -1;
    RealType outletDirX = 1, outletDirY = 0, outletDirZ = 0;
    bool outletDoNothing = true;  // do-nothing outlet (free velocity + p=0 face); false = mass-conserving velocity outlet

    // Constant streamwise (and orthogonal) momentum body force, added in the
    // predictor as +dt*f/rho (BDF1) and +(2dt/3)*f/rho (BDF2). Defaults to 0
    // to match Nalu-Wind / NekRS / MFEM-Navier convention (no force unless the
    // user opts in via CLI). Used for full-Dirichlet flows where a
    // uniform IC + full-Dirichlet BC has no transient otherwise; for cavity
    // the lid shears the flow and the body force stays 0. Activates per
    // component via --body-force-x|y|z=<val>.
    RealType bodyForceX = 0;
    RealType bodyForceY = 0;
    RealType bodyForceZ = 0;
    // IC velocity perturbation amplitude as a fraction of Uinf. Adds a small
    // deterministic per-node perturbation to (v, w) interior at t=0 so the
    // body force / advection has something to amplify. Without this on a
    // uniform-IC + full-Dirichlet domain, the system is a discrete
    // steady state and no flow develops -- exactly matches NekRS userdat2 /
    // Nalu-Wind ABLForcingMomentum boot pattern. Default 0 (opt-in via
    // --ic-perturb=<eps>). Standard NekRS value is ~1e-3.
    RealType icPerturbMag = 0;

    // Pump-mode side-set local-node lists. Each vector holds owned-or-ghost
    // LOCAL node indices (size_t) that the driver mapped from global Exodus
    // node IDs via cstone's d_localToGlobalNodeMap_. Populated only when
    // bcKind == Pump. Names match the Exodus mesh side-set names.
    std::vector<int> wallNodes;     // no-slip
    std::vector<int> inletNodes;    // u = Uinf
    std::vector<int> outletNodes;   // Dirichlet p = 0
    std::vector<int> extraNodes;    // top/bottom/side/sym -- u = Uinf (full Dirichlet shortcut)

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

    // Ownership map every kernel uses. Always the domain's cached map -- we no
    // longer demote periodic slaves (see setupNSStepper). Kept as a single
    // accessor so all kernels bind through one place.
    const cstone::DeviceVector<uint8_t>& ownershipMap() const {
        return domain.getNodeOwnershipMap();
    }
};

// =============================================================================
// Helpers for the two-matrix Laplacian assembly. Both matrices want the same
// CVFEM Laplacian K (gamma=1), assembled by CvfemHexAssembler::assembleFull
// with zeroed advection inputs. We assemble K once into a temporary buffer,
// then copy it into d_valuesVel (and add M/dt + enforce velocity BC) and into
// d_valuesPre (and enforce single-pin BC) so the two matrices end up with
// different values but the same row/column structure.
// =============================================================================

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void assembleLaplacian(NSStepper<KeyType, RealType, ElementTag>& s,
                       const cstone::DeviceVector<int>& /*unused*/,
                       cstone::DeviceVector<RealType>& d_valuesOut,
                       CvfemKernelVariant kernelVariant)
{
    const auto& d_nodeOwnership = s.ownershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();

    const uint8_t* effectiveOwnership = d_nodeOwnership.data();

    using MatrixT = CSRMatrix<RealType>;
    MatrixT* d_matrix;
    cudaMalloc(&d_matrix, sizeof(MatrixT));
    // numRows = numTotalDofs (matrix is sized for owned + ghost rows so ghost
    // column indices are addressable). numOwnedRows = numOwnedDofs tells the
    // kernel which rows are actually solved + halo-exchanged; scatters into
    // rows >= numOwnedRows would be silent leaks (e.g. periodic slave redirected
    // to a cross-rank master ghost).
    MatrixT h_matrix{s.d_rowPtr.data(), s.d_colInd.data(), d_valuesOut.data(), s.d_diagPtr.data(),
                     s.numTotalDofs, s.nnz, s.numOwnedDofs};
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixT), cudaMemcpyHostToDevice);

    // gamma=1 -> pure Laplacian; all advection inputs zero so assembler doesn't
    // accumulate an advective contribution. RHS unused.
    cstone::DeviceVector<RealType> d_gamma(s.nodeCount, RealType(1));
    cstone::DeviceVector<RealType> d_phi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_beta(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_gphi_x(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_gphi_y(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_gphi_z(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_mdot_zero(s.elementCount * ElemTraits<ElementTag>::ScsPerElem, RealType(0));
    cstone::DeviceVector<RealType> d_rhs_unused(s.numTotalDofs, RealType(0));

    thrust::fill(thrust::device_pointer_cast(d_valuesOut.data()),
                 thrust::device_pointer_cast(d_valuesOut.data() + s.nnz),
                 RealType(0));

    if constexpr (std::is_same_v<ElementTag, HexTag>)
    {
        typename CvfemHexAssembler<KeyType, RealType>::Config config;
        config.blockSize = s.blockSize;
        config.variant   = kernelVariant;

        CvfemHexAssembler<KeyType, RealType>::assembleFull(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
            std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), s.elementCount,
            d_x.data(), d_y.data(), d_z.data(),
            d_gamma.data(), d_phi.data(), d_beta.data(),
            d_gphi_x.data(), d_gphi_y.data(), d_gphi_z.data(),
            d_mdot_zero.data(), s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            s.d_node_to_dof.data(), effectiveOwnership, d_matrix, d_rhs_unused.data(), config);
    }
    else
    {
        typename CvfemTetAssembler<KeyType, RealType>::Config tetConfig;
        tetConfig.blockSize = s.blockSize;

        CvfemTetAssembler<KeyType, RealType>::assembleFull(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(), s.elementCount,
            d_x.data(), d_y.data(), d_z.data(),
            d_gamma.data(), d_phi.data(), d_beta.data(),
            d_gphi_x.data(), d_gphi_y.data(), d_gphi_z.data(),
            d_mdot_zero.data(),
            s.d_node_to_dof.data(), effectiveOwnership, d_matrix, d_rhs_unused.data(), tetConfig);
    }
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void addBochevDohrmannStab(NSStepper<KeyType, RealType, ElementTag>& s, RealType tau)
{
    const auto& d_nodeOwnership = s.ownershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();

    size_t nElem = s.elementCount;
    int blockSize = s.blockSize;
    int grid = int((nElem + blockSize - 1) / blockSize);
    if (grid <= 0) return;

    // BD stabilization is a Q1-hex-specific 8x8 scatter (std::get<7>, 8 corners).
    // Tet pressure stabilization is a separate follow-up; guard so the tet
    // instantiation never touches the 8-column connectivity.
    if constexpr (std::is_same_v<ElementTag, HexTag>)
    {
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
}

// =============================================================================
// Setup: builds DOF mapping, sparsity, two matrices, lumped mass, BC mask,
// pressure pin, area vectors. Runs ONCE before the time loop.
// =============================================================================

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void setupNSStepper(NSStepper<KeyType, RealType, ElementTag>& s,
                    RealType nu,
                    RealType dt,
                    CvfemKernelVariant kernelVariant)
{
    // Live-print each setup lap on rank 0; the final report still summarizes.
    PhaseTimer pt(/*sync=*/true, /*liveRank=*/0);

    s.nodeCount    = s.domain.getNodeCount();
    s.elementCount = s.domain.getElementCount();

    // Optional one-shot diagnostic: count halo elements per rank and how many
    // touch the periodic-min faces. Confirms cstone delivers periodic-image
    // halo elements (verified yes on cube16/4-rank: 1238 halo elements touch
    // min-faces). Set MARS_PERIODIC_HALO_DBG=1 to enable.
    // Hex-only: the lambda spells std::get<7> and a fixed n[8]; a runtime
    // `if` would still INSTANTIATE that on the 4-tuple tet domain. if constexpr
    // keeps it out of the tet compile entirely.
    if constexpr (std::is_same_v<ElementTag, HexTag>)
    if (std::getenv("MARS_PERIODIC_HALO_DBG") != nullptr)
    {
        size_t startE  = s.domain.startIndex();
        size_t endE    = s.domain.endIndex();
        size_t numOwn  = endE - startE;
        size_t numHalo = s.elementCount > numOwn ? s.elementCount - numOwn : 0;

        // Local bbox + Allreduce to get global box (s.xmin etc. not yet set).
        const auto& d_cn = s.domain.getElementToNodeConnectivity();
        const auto& d_x  = s.domain.getNodeX();
        const auto& d_y  = s.domain.getNodeY();
        const auto& d_z  = s.domain.getNodeZ();
        RealType lxmin = thrust::reduce(thrust::device, d_x.data(), d_x.data() + s.nodeCount, RealType(1e30),  thrust::minimum<RealType>());
        RealType lxmax = thrust::reduce(thrust::device, d_x.data(), d_x.data() + s.nodeCount, RealType(-1e30), thrust::maximum<RealType>());
        RealType lymin = thrust::reduce(thrust::device, d_y.data(), d_y.data() + s.nodeCount, RealType(1e30),  thrust::minimum<RealType>());
        RealType lymax = thrust::reduce(thrust::device, d_y.data(), d_y.data() + s.nodeCount, RealType(-1e30), thrust::maximum<RealType>());
        RealType lzmin = thrust::reduce(thrust::device, d_z.data(), d_z.data() + s.nodeCount, RealType(1e30),  thrust::minimum<RealType>());
        RealType lzmax = thrust::reduce(thrust::device, d_z.data(), d_z.data() + s.nodeCount, RealType(-1e30), thrust::maximum<RealType>());
        RealType xmin, xmax, ymin, ymax, zmin, zmax;
        MPI_Datatype mpiType = std::is_same<RealType, double>::value ? MPI_DOUBLE : MPI_FLOAT;
        MPI_Allreduce(&lxmin, &xmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&lxmax, &xmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&lymin, &ymin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&lymax, &ymax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&lzmin, &zmin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&lzmax, &zmax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
        RealType eps  = RealType(1e-3) * std::max({xmax - xmin, ymax - ymin, zmax - zmin});

        long long localOnMin = 0;
        if (numHalo > 0)
        {
            localOnMin = thrust::count_if(
                thrust::device,
                thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(s.elementCount),
                [startE, endE,
                 c0 = std::get<0>(d_cn).data(), c1 = std::get<1>(d_cn).data(),
                 c2 = std::get<2>(d_cn).data(), c3 = std::get<3>(d_cn).data(),
                 c4 = std::get<4>(d_cn).data(), c5 = std::get<5>(d_cn).data(),
                 c6 = std::get<6>(d_cn).data(), c7 = std::get<7>(d_cn).data(),
                 x = d_x.data(), y = d_y.data(), z = d_z.data(),
                 xmin, ymin, zmin, eps]
                __device__ (size_t e) -> bool {
                    if (e >= startE && e < endE) return false;
                    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e],
                                    c4[e], c5[e], c6[e], c7[e]};
                    for (int i = 0; i < 8; ++i)
                    {
                        if (fabs(double(x[n[i]] - xmin)) < double(eps) ||
                            fabs(double(y[n[i]] - ymin)) < double(eps) ||
                            fabs(double(z[n[i]] - zmin)) < double(eps))
                            return true;
                    }
                    return false;
                });
        }

        long long lOwn  = (long long)numOwn;
        long long lHalo = (long long)numHalo;
        long long gOwn = 0, gHalo = 0, gOnMin = 0;
        MPI_Reduce(&lOwn,       &gOwn,   1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&lHalo,      &gHalo,  1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&localOnMin, &gOnMin, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        if (s.rank == 0)
        {
            std::cout << "[DBG halo] sum owned=" << gOwn
                      << " sum halo=" << gHalo
                      << " sum halo-on-min-face=" << gOnMin
                      << "  (rank0 own=" << lOwn << " halo=" << lHalo
                      << " onMin=" << localOnMin << ")\n";
        }
    }

    // Ownership: normally the domain's cached map. For PERIODIC multi-rank we
    // Periodic slaves are NOT demoted. MARS's multi-rank matrix-assembly
    // invariant requires every owned DOF to receive all its element stiffness
    // locally (the assembler loops over cstone's halo-extended element set and
    // scatters owned rows only -- no cross-rank matrix coupling). With the
    // periodic cstone Box, the periodic-image element + master node arrive in
    // each rank's halo, so a periodic master's owned row is filled on its
    // owner rank exactly as for an interior node. Slaves keep ownership=1 and
    // the DOF collapse below identifies slave_dof with the (now locally
    // present, owned-or-halo) master_dof. Demoting slaves broke this invariant
    // (under-filled master rows -> velocity CG NaN), so we use the domain's
    // ownership map unchanged.
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

    // PERIODIC DOF COLLAPSE: a slave node shares its master's DOF index, then
    // we COMPACT the DOF numbering so there are no orphan rows. After this:
    //   - nodeToDof[slave] = nodeToDof[master]  (only for LOCAL-OWNED masters)
    //   - all surviving DOF indices are dense in [0, nLive)
    //   - s.numOwnedDofs = nLive
    //
    // MULTI-RANK: we ONLY collapse a slave onto its master when the master is
    // locally OWNED (own==1). Then the merged DOF is a real owned row on this
    // rank and the assembler fills it from the local halo-extended element set.
    //
    // When the master is a REMOTE node (present only as a periodic-halo ghost
    // here, own!=1), we do NOT collapse: the slave keeps its own owned DOF.
    // The slave and the remote master are the same periodic point but each
    // rank owns its own copy's row -- exactly how cstone treats any shared
    // node. cstone's periodic element halo fills each owned row locally, and
    // exchangeNodeHalo keeps the slave's value synced to the master's owner.
    // Collapsing across ranks would point the slave at a ghost DOF (out of the
    // owned range) and orphan its row -> velocity CG NaN (the bug we hit).
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic && s.periodicMap)
    {
        const int* d_partner = s.periodicMap->d_periodicPartner.data();
        const uint8_t* d_own = d_nodeOwnership.data();
        int* d_n2d           = s.d_node_to_dof.data();
        size_t nNodes        = s.nodeCount;
        int numOwned         = s.numOwnedDofs;

        // Pass 1: redirect slave -> master DOF ONLY for SAME-RANK pairs (both
        // the slave and the master are locally owned). Cross-rank pairs (only
        // one is owned, the other is a ghost via cstone's periodic-image halo)
        // are NOT collapsed. Doing so would break the matrix symmetry:
        //   - On rank A (slave owner, master is ghost): the slave row would
        //     have a column for the master-ghost, but rank D's master row
        //     would lose its column for the slave-ghost (because the redirect
        //     maps slave-ghost-DOF back onto master-owned-DOF, eliminating
        //     the off-diagonal slot).
        //
        // For cross-rank pairs, each owner keeps its own owned DOF and its
        // own row. The slave-master identity is enforced post-solve via
        // crossRankPeriodicBroadcastDof: x[slave_on_A] := x[master_owned_on_D].
        // The matrix stays SPD with paired (slave,master) and (master,slave)
        // off-diagonal entries on the two ranks.
        cstone::DeviceVector<int> d_oldDof(s.d_node_to_dof);
        thrust::for_each(thrust::device,
                         thrust::counting_iterator<size_t>(0),
                         thrust::counting_iterator<size_t>(nNodes),
                         [d_partner, d_old = d_oldDof.data(), d_n2d, d_own]
                         __device__ (size_t i) {
                             int master = d_partner[i];
                             // Same-rank pair only: BOTH the current node (slave)
                             // and the master must be locally owned.
                             if (master >= 0 && d_own[master] == 1 && d_own[i] == 1)
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

        // Pass 4: rewrite nodeToDof through the compact map.
        // - Owned nodes whose dof is in [0, numOwned): apply compact[].
        // - Owned slaves whose Pass-1 redirect put them at a ghost dof
        //   (>= numOwned, cross-rank periodic case): skip; the ghost dof is
        //   already in the post-compact ghost range [nLive, numTotalDofs)
        //   semantically because d_n2d ghost slots don't get renumbered and
        //   cstone halo exchange writes the ghost dof slot from the master
        //   owner's owned-dof slot (which IS in [0, nLive_remote) on the
        //   master rank). This works because we never reorder the ghost-dof
        //   slots themselves -- only the owned slots compact.
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
        {
            // Per-rank print: cross-rank slaves only show on the slave-owner rank,
            // not rank 0. Keep all ranks visible to verify Pass-1 redirect coverage.
            int nOrphan = numOwned - nLive;
            std::cout << "  [rank " << s.rank << "] periodic DOF collapse: "
                      << nOrphan << " slave DOFs collapsed -> "
                      << nLive << " active equations (was " << numOwned << ")\n";
        }
    }
    pt.lap("DOF mapping");

    // DIAGNOSTIC: owned-node partition health. Sum of per-rank owned DOFs across
    // ranks MUST equal the global unique-DOF count -- if it is LARGER, some nodes
    // are doubly-owned (counted on >1 rank), which double-counts them in the CG
    // inner products and corrupts the projection (multi-rank div != single-rank).
    {
        long long localOwned = s.numOwnedDofs;
        long long sumOwned = 0, maxOwned = 0;
        MPI_Allreduce(&localOwned, &sumOwned, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localOwned, &maxOwned, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        if (s.rank == 0)
            std::cout << "  [owned-node check] sum(numOwnedDofs over ranks)=" << sumOwned
                      << "  (this should EQUAL the global unique node/DOF count; "
                      << "larger => doubly-owned boundary nodes)\n";
    }

    // FULL 27-NNZ sparsity (same pattern used by both matrices). Building once
    // means the velocity and pressure matrices share row/col/diag pointers.
    s.d_rowPtr.resize(s.numTotalDofs + 1);
    s.d_diagPtr.resize(s.numTotalDofs);
    if constexpr (std::is_same_v<ElementTag, HexTag>)
    {
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
    }
    else
    {
        s.nnz = CvfemTetSparsityBuilder<KeyType>::buildFullSparsity(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(), s.elementCount,
            s.d_node_to_dof.data(), s.numTotalDofs,
            s.d_rowPtr.data(), nullptr, nullptr, 0);
        s.d_colInd.resize(s.nnz);
        CvfemTetSparsityBuilder<KeyType>::buildFullSparsity(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(), s.elementCount,
            s.d_node_to_dof.data(), s.numTotalDofs,
            s.d_rowPtr.data(), s.d_colInd.data(), s.d_diagPtr.data(), 0);
    }
    s.d_valuesVel.resize(s.nnz);
    s.d_valuesPre.resize(s.nnz);
    pt.lap("sparsity build");

    // Area vectors (single source of truth for face geometry across operators).
    s.d_areaVec_x.resize(s.elementCount * ElemTraits<ElementTag>::ScsPerElem);
    s.d_areaVec_y.resize(s.elementCount * ElemTraits<ElementTag>::ScsPerElem);
    s.d_areaVec_z.resize(s.elementCount * ElemTraits<ElementTag>::ScsPerElem);
    if constexpr (std::is_same_v<ElementTag, HexTag>)
    {
        precomputeAreaVectorsGpu<KeyType, RealType>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
            std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), s.elementCount,
            d_x.data(), d_y.data(), d_z.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data());
    }
    else
    {
        precomputeTetAreaVectorsGpu<KeyType, RealType>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(), s.elementCount,
            d_x.data(), d_y.data(), d_z.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data());
    }
    pt.lap("area vectors");

    // Velocity matrix: nu * K assembled first (then we add M/dt and apply BC).
    assembleLaplacian<KeyType, RealType, ElementTag>(s, s.d_node_to_dof, s.d_valuesVel, kernelVariant);
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

    // Cross-rank periodic seam: ship rank A's slave-row stiffness to the
    // master-owner rank and accumulate into the master row. cstone delivers
    // the periodic-image elements as halos, but the assembler's ownership
    // gate skips them on the master-owner rank (the slave-image node is a
    // ghost there). Without this MPI exchange the master row is missing
    // ~70% of its off-diagonal stiffness. One-shot at setup. No-op when
    // numRanks==1 or no cross-rank peers exist.
    // Pressure matrix: same K (no nu scale). Assemble into d_valuesPre.
    assembleLaplacian<KeyType, RealType, ElementTag>(s, s.d_node_to_dof, s.d_valuesPre, kernelVariant);
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
        // OWNED-only per-element scatter + reverseExchangeNodeHaloAdd: standard
        // pattern for cross-rank shared-FACE nodes. Periodic-pair cross-rank
        // contributions are NOT routed this way (slave/master have distinct
        // SFC keys, see domain.hpp comment around line 1607) -- a separate
        // periodic-pair MPI exchange runs below.
        size_t startElem = s.domain.startIndex();
        size_t numLocal  = s.domain.localElementCount();
        if (numLocal > 0)
        {
            int eBlocks = int((numLocal + s.blockSize - 1) / s.blockSize);
            auto cp = connPtrs<ElementTag, KeyType>(d_conn);
            if constexpr (std::is_same_v<ElementTag, HexTag>)
                computeLumpedMassPerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                    cp[0], cp[1], cp[2], cp[3], cp[4], cp[5], cp[6], cp[7],
                    d_x.data(), d_y.data(), d_z.data(),
                    s.d_massNode.data(), startElem, numLocal);
            else
                computeLumpedMassPerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                    cp[0], cp[1], cp[2], cp[3], nullptr, nullptr, nullptr, nullptr,
                    d_x.data(), d_y.data(), d_z.data(),
                    s.d_massNode.data(), startElem, numLocal);
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(s.d_massNode);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, s.d_massNode);
        // Forward halo: ghost slots get owner-rank V values so DDT assembler
        // can read s.d_massNode[ghostNode] and get the correct V.
        s.domain.exchangeNodeHalo(s.d_massNode);
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        gatherOwnedNodeMassToDofKernel<RealType><<<nBlocks, s.blockSize>>>(
            s.d_massNode.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_mass.data(), s.nodeCount, s.numOwnedDofs);
        cudaDeviceSynchronize();

        // DIAGNOSTIC: total control volume = sum of OWNED-node lumped mass over
        // ALL ranks. This MUST equal the single-rank total (the domain volume) --
        // it is rank-count-invariant. If multi-rank sum < single-rank, boundary
        // nodes are MISSING off-rank element volume => their 1/V_i is too large
        // => grad(p) and the corrector blow up there (seen as |gP|max ~ 1e4).
        {
            const uint8_t* ownPtr = d_nodeOwnership.data();
            const RealType* mP    = s.d_massNode.data();
            double localMass = thrust::transform_reduce(thrust::device,
                thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(s.nodeCount),
                [ownPtr, mP] __device__ (size_t i) -> double {
                    return (ownPtr[i] == 1) ? double(mP[i]) : 0.0;
                }, 0.0, thrust::plus<double>());
            double globalMass = 0;
            MPI_Allreduce(&localMass, &globalMass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (s.rank == 0)
                std::cout << "  [mass-sum check] sum(owned d_massNode over ranks)=" << std::scientific
                          << std::setprecision(8) << globalMass << std::defaultfloat
                          << "  (rank-count-INVARIANT = domain volume; if 4-rank < 1-rank, "
                          << "boundary nodes miss off-rank volume)\n";
        }
    }
    pt.lap("lumped mass");

    // Optional: print min/max/mean of d_mass + count of slots at "full" mass
    // (= V_e for the unit-mesh cube16 case) vs "half" mass. Distinguishes
    // periodic-pair sums "actually happening" from "missing".
    if (std::getenv("MARS_PERIODIC_HALO_DBG") != nullptr)
    {
        thrust::host_vector<RealType> h_mass(s.numOwnedDofs);
        thrust::copy(thrust::device_pointer_cast(s.d_mass.data()),
                     thrust::device_pointer_cast(s.d_mass.data() + s.numOwnedDofs),
                     h_mass.begin());
        double lmin = 1e30, lmax = -1e30, lsum = 0.0;
        for (size_t i = 0; i < h_mass.size(); ++i) {
            double v = double(h_mass[i]);
            if (v < lmin) lmin = v;
            if (v > lmax) lmax = v;
            lsum += v;
        }
        double gmin, gmax, gsum;
        long long lcount = (long long)h_mass.size(), gcount = 0;
        MPI_Reduce(&lmin, &gmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&lmax, &gmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&lsum, &gsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&lcount, &gcount, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        if (s.rank == 0) {
            std::cout << "[DBG mass] global owned-DOF count=" << gcount
                      << std::scientific << std::setprecision(6)
                      << " min=" << gmin
                      << " max=" << gmax
                      << " sum=" << gsum
                      << std::defaultfloat
                      << " (expected for cube16/unit-box: sum=1.0; full V_e ~ 2.44e-4; half V_e ~ 1.22e-4)\n";
        }
    }

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
    //
    // Hex-only: buildDDTSparsity / assembleDDTPerNodeKernel spell std::get<7>
    // and the fixed hex 8-corner / nodeFaces[8][3] incidence. The tet DDT
    // operator is a separate follow-up; if constexpr keeps the whole block
    // (sparsity, node-driven assembly, health check, diagonal shift) out of
    // the tet compile entirely.
    if constexpr (std::is_same_v<ElementTag, HexTag>)
    {
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
        // Post-assembly health check: every owned row must have a non-zero
        // diagonal entry (Hypre BoomerAMG strength-of-connection asserts on
        // |A[i,i]| > 0; missing or zero diag is a documented cause of
        // HYPRE_ERROR_GENERIC on setup). Run on rank 0 only since this is
        // a structural property -- all ranks should see the same answer for
        // their own owned rows. (Print only if a defect is found.)
        if (s.rank == 0)
        {
            std::vector<int> hRowPtr(s.numOwnedDofs + 1);
            cudaMemcpy(hRowPtr.data(), s.d_rowPtrDDT.data(),
                       (s.numOwnedDofs + 1) * sizeof(int), cudaMemcpyDeviceToHost);
            std::vector<int> hColInd(hRowPtr[s.numOwnedDofs]);
            std::vector<RealType> hVals(hRowPtr[s.numOwnedDofs]);
            cudaMemcpy(hColInd.data(), s.d_colIndDDT.data(),
                       hRowPtr[s.numOwnedDofs] * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(hVals.data(), s.d_valuesDDT.data(),
                       hRowPtr[s.numOwnedDofs] * sizeof(RealType), cudaMemcpyDeviceToHost);
            int emptyRows = 0, missingDiag = 0, zeroDiag = 0, negDiag = 0;
            int firstEmpty = -1, firstMissingDiag = -1;
            RealType minDiag = 1e30, maxDiag = -1e30;
            for (int r = 0; r < s.numOwnedDofs; ++r)
            {
                int rs = hRowPtr[r], re = hRowPtr[r + 1];
                if (re == rs) { emptyRows++; if (firstEmpty<0) firstEmpty=r; continue; }
                int diagIdx = -1;
                for (int k = rs; k < re; ++k) if (hColInd[k] == r) { diagIdx = k; break; }
                if (diagIdx < 0) {
                    missingDiag++;
                    if (firstMissingDiag<0) firstMissingDiag=r;
                    continue;
                }
                RealType d = hVals[diagIdx];
                if (d == RealType(0)) zeroDiag++;
                else if (d < RealType(0)) negDiag++;
                if (d < minDiag) minDiag = d;
                if (d > maxDiag) maxDiag = d;
            }
            std::cout << "  [DDT-health] emptyRows=" << emptyRows
                      << " missingDiag=" << missingDiag
                      << " zeroDiag=" << zeroDiag
                      << " negDiag=" << negDiag
                      << " diag range=[" << minDiag << ", " << maxDiag << "]\n";
            if (emptyRows > 0)
                std::cout << "    first empty row: " << firstEmpty << "\n";
            if (missingDiag > 0)
                std::cout << "    first missing-diag row: " << firstMissingDiag << "\n";
        }
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
    } // end if constexpr hex (DDT operator block)

    // Add M/dt to the velocity matrix diagonal -> (M/dt + nu K).
    {
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        addLumpedMassDiagonalKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_mass.data(), s.d_diagPtr.data(), RealType(1) / dt,
            s.d_valuesVel.data(), s.numOwnedDofs);
        cudaDeviceSynchronize();
    }
    pt.lap("add M/dt (velocity)");

    // DIAGNOSTIC: velocity diffusion matrix (nu K + M/dt) INTERIOR-row consistency
    // across ranks. K is a Laplacian -> its rows sum ~0, so each interior row of
    // (nu K + M/dt) sums to ~ mass[row]/dt -- a PHYSICAL, rank-count-INVARIANT
    // value. The matrix is assembled over OWNED+HALO elements; if a halo element
    // is double-counted (or an owned row's stencil differs across ranks), the
    // interior row-sum will DIFFER on 4 ranks vs 1 -> the diffusion solve injects
    // a bulk perturbation every step (suspected cause of div(u*) 17x on 4 ranks).
    // We report global max |rowsum - mass/dt| over interior owned rows: ~0 if
    // assembly is rank-consistent, large if not.
    {
        // Compare the velocity-matrix DIAGONAL (nu K_diag + mass/dt) across ranks
        // on a per-DOF basis. Uses d_diagPtr indirection -- the SAME proven-safe
        // indexing as addLumpedMassDiagonalKernel above (diagPtr[dof], guarded
        // dp>=0) -- so no rowPtr walk (that crashed: tet rowPtr layout is not a
        // plain owned-row CSR here). We report max and SUM of the K-diagonal
        // part (diag - mass/dt) over owned DOFs; both are rank-count-INVARIANT if
        // the owned+halo assembly is consistent. A rank-VARYING sum => an owned
        // row's stencil differs across ranks (suspected div(u*) 17x cause).
        const int* dptr      = s.d_diagPtr.data();
        const RealType* vv   = s.d_valuesVel.data();
        const RealType* mass = s.d_mass.data();
        RealType invdt = RealType(1) / dt;
        double localKsum = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(s.numOwnedDofs),
            [dptr, vv, mass, invdt] __device__ (int r) -> double {
                int dp = dptr[r];
                if (dp < 0) return 0.0;
                return double(vv[dp] - mass[r] * invdt);   // nu*K diagonal part
            }, 0.0, thrust::plus<double>());
        RealType localKmax = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(s.numOwnedDofs),
            [dptr, vv, mass, invdt] __device__ (int r) -> RealType {
                int dp = dptr[r];
                if (dp < 0) return RealType(0);
                return fabs(vv[dp] - mass[r] * invdt);
            }, RealType(0), thrust::maximum<RealType>());
        auto mt = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
        double gKsum = 0; RealType gKmax = 0;
        MPI_Allreduce(&localKsum, &gKsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&localKmax, &gKmax, 1, mt, MPI_MAX, MPI_COMM_WORLD);
        if (s.rank == 0)
            std::cout << "  [Avel-diag check] sum(nu*K_diag over owned)=" << std::scientific
                      << std::setprecision(8) << gKsum << "  max=" << gKmax << std::defaultfloat
                      << "  (BOTH rank-count-INVARIANT if owned+halo assembly is consistent; "
                      << "rank-varying => the bug)\n";
    }

    // Probe: dump master row stats on each rank that owns masters. Gated by
    // MARS_PERIODIC_AVEL_DBG=1. Distinguishes:
    //   (a) diag == expected single-rank value + sumAbsOff > 0 -> seam coupling
    //       made it into the master row. Bug is elsewhere.
    //   (b) diag << expected -> slave-side stiffness never reached master row.
    //       Either cstone didn't deliver the periodic-image element to this
    //       rank, OR the assembler on this rank skipped its scatter.
    if (std::getenv("MARS_PERIODIC_AVEL_DBG") != nullptr
        && s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
        && s.periodicMap != nullptr
        && !s.periodicMap->cross_.d_recvOwnedMasterIds_.empty())
    {
        const auto& xr = s.periodicMap->cross_;
        int firstMasterNode = -1;
        cudaMemcpy(&firstMasterNode, xr.d_recvOwnedMasterIds_.data(),
                   sizeof(int), cudaMemcpyDeviceToHost);
        if (firstMasterNode >= 0 && firstMasterNode < int(s.nodeCount))
        {
            int masterDof = -1;
            cudaMemcpy(&masterDof, s.d_node_to_dof.data() + firstMasterNode,
                       sizeof(int), cudaMemcpyDeviceToHost);
            if (masterDof >= 0 && masterDof < s.numOwnedDofs)
            {
                int rowOff[2] = {0, 0};
                cudaMemcpy(rowOff, s.d_rowPtr.data() + masterDof,
                           2 * sizeof(int), cudaMemcpyDeviceToHost);
                int rowLen = rowOff[1] - rowOff[0];
                std::vector<RealType> rowVals(rowLen);
                cudaMemcpy(rowVals.data(),
                           s.d_valuesVel.data() + rowOff[0],
                           rowLen * sizeof(RealType), cudaMemcpyDeviceToHost);
                int diagIdx = -1;
                cudaMemcpy(&diagIdx, s.d_diagPtr.data() + masterDof,
                           sizeof(int), cudaMemcpyDeviceToHost);
                double diag = 0.0, sumAbsOff = 0.0, maxAbsOff = 0.0;
                int nnzOff = 0;
                int diagLocal = diagIdx - rowOff[0];
                for (int j = 0; j < rowLen; ++j)
                {
                    double v = double(rowVals[j]);
                    if (j == diagLocal) { diag = v; }
                    else if (v != 0.0)
                    {
                        ++nnzOff;
                        sumAbsOff += std::fabs(v);
                        if (std::fabs(v) > maxAbsOff) maxAbsOff = std::fabs(v);
                    }
                }
                std::cerr << "[avel-probe rank " << s.rank
                          << "] masterNode=" << firstMasterNode
                          << " masterDof=" << masterDof
                          << " rowLen=" << rowLen
                          << " diag=" << diag
                          << " nnzOff=" << nnzOff
                          << " sumAbsOff=" << sumAbsOff
                          << " maxAbsOff=" << maxAbsOff
                          << " (expect interior diag ~ 2.45 nnzOff~17 sumAbs~0.04)"
                          << std::endl;
            }
        }
    }

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
        using BCK = typename NSStepper<KeyType, RealType, ElementTag>::BCKind;
        if (s.bcKind == BCK::Cavity)
        {
            markCavityBCKernel<RealType><<<nBlocks, s.blockSize>>>(
                d_x.data(), d_y.data(), d_z.data(),
                d_nodeOwnership.data(), s.d_node_to_dof.data(),
                s.d_isBdryDof.data(),
                s.d_uTarget.data(), s.d_vTarget.data(), s.d_wTarget.data(),
                s.nodeCount, s.xmin, s.xmax, s.ymin, s.ymax, s.zmin, s.zmax,
                s.bboxEps, s.lidU, s.numOwnedDofs);
        }
        else if (s.bcKind == BCK::Channel)
        {
            markChannelBCKernel<RealType><<<nBlocks, s.blockSize>>>(
                d_x.data(), d_y.data(), d_z.data(),
                d_nodeOwnership.data(), s.d_node_to_dof.data(),
                s.d_isBdryDof.data(),
                s.d_uTarget.data(), s.d_vTarget.data(), s.d_wTarget.data(),
                s.nodeCount, s.xmin, s.xmax, s.ymin, s.ymax, s.zmin, s.zmax,
                s.bboxEps, s.Uinf, s.numOwnedDofs);
        }
        else if (s.bcKind == BCK::Pump)
        {
            // BC: per-side-set Dirichlet driven by the local-node lists
            // populated by the driver. Build host arrays then copy to device:
            //   - isBdryDof[dof] = 1 for nodes on walls/inlet/extra.
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
            tag(s.wallNodes,    RealType(0),     RealType(0), RealType(0));
            tag(s.inletNodes,   RealType(s.Uinf * s.inletDirX),
                                RealType(s.Uinf * s.inletDirY),
                                RealType(s.Uinf * s.inletDirZ));
            tag(s.extraNodes,   RealType(s.Uinf), RealType(0), RealType(0));
            // Mass-conserving outlet: velocity-Dirichlet outflow along the
            // outward normal, sized to remove the inlet flux. Only when the
            // driver set outletU > 0; otherwise the outlet stays natural-Neumann
            // (legacy) and only gets the p=0 pressure Dirichlet below.
            if (s.outletU > RealType(0))
                tag(s.outletNodes, RealType(s.outletU * s.outletDirX),
                                       RealType(s.outletU * s.outletDirY),
                                       RealType(s.outletU * s.outletDirZ));

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
                std::cout << "  BC: " << s.wallNodes.size()  << " wall nodes (no-slip), "
                          << s.inletNodes.size() << " inlet nodes (u=" << s.Uinf << "), "
                          << s.extraNodes.size() << " extra nodes (u=" << s.Uinf << ")\n";
            }
        }
        // Periodic: no boundary DOFs (mask stays zero) and no per-node targets
        // (targets stay zero; predictor/corrector write the field on every node).
        cudaDeviceSynchronize();
    }

    // DIAGNOSTIC: global count of OWNED Dirichlet-velocity DOFs (d_isBdryDof==1)
    // and the owned-sum of |uTarget|^2. Both are rank-count-INVARIANT: a global
    // node is owned by exactly one rank, so summed over ranks the tagged-DOF
    // count and the target energy must MATCH the single-rank values. If 4-rank
    // differs from 1-rank, BC nodes are being dropped or double-counted across
    // ranks -- the direct cause of u* (and KE*) differing from step 1.
    {
        const uint8_t* bnd = s.d_isBdryDof.data();
        long long locBnd = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<int>(0),
            thrust::counting_iterator<int>(s.numOwnedDofs),
            [bnd] __device__ (int d) -> long long { return bnd[d] ? 1LL : 0LL; },
            0LL, thrust::plus<long long>());
        // owned-sum of |uTarget|^2 over Dirichlet nodes (positive => no cancel).
        const uint8_t* ownPtr = d_nodeOwnership.data();
        const int* dofPtr     = s.d_node_to_dof.data();
        const uint8_t* bnode  = (s.d_isBdryNode.size() == s.nodeCount) ? s.d_isBdryNode.data() : nullptr;
        const RealType* uT = s.d_uTarget.data();
        const RealType* vT = s.d_vTarget.data();
        const RealType* wT = s.d_wTarget.data();
        double locTgt = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [ownPtr, dofPtr, uT, vT, wT] __device__ (size_t i) -> double {
                if (ownPtr[i] != 1 || dofPtr[i] < 0) return 0.0;
                double u = uT[i], v = vT[i], w = wT[i];
                return u*u + v*v + w*w;
            }, 0.0, thrust::plus<double>());
        (void)bnode;
        long long gBnd = 0; double gTgt = 0;
        MPI_Allreduce(&locBnd, &gBnd, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locTgt, &gTgt, 1, MPI_DOUBLE,    MPI_SUM, MPI_COMM_WORLD);
        if (s.rank == 0)
            std::cout << "  [bc-count check] owned Dirichlet-vel DOFs (sum over ranks)=" << gBnd
                      << "  sum|uTarget|^2(owned)=" << std::scientific << std::setprecision(8) << gTgt
                      << std::defaultfloat << "  (BOTH rank-count-INVARIANT; differ 1 vs 4 => BC nodes dropped/doubled)\n";
    }

    // Sync ghost slots of the target arrays so per-node corrector/predictor
    // can read them on owned-boundary nodes even when those nodes are touched
    // by ghost elements on other ranks.
    s.domain.exchangeNodeHalo(s.d_uTarget);
    s.domain.exchangeNodeHalo(s.d_vTarget);
    s.domain.exchangeNodeHalo(s.d_wTarget);
    pt.lap("BC mark + target velocity");

    // d_dofToNode is needed for column-zero enforcement below. The pressure
    // DDT path builds it later (line ~3291) but we need it earlier for the
    // velocity periodic-XR column-zero step. Build once here; subsequent
    // re-builds in the DDT block are idempotent (same input, same output).
    if (s.d_dofToNode.size() != size_t(s.numTotalDofs))
    {
        s.d_dofToNode.resize(s.numTotalDofs);
        thrust::fill(thrust::device_pointer_cast(s.d_dofToNode.data()),
                     thrust::device_pointer_cast(s.d_dofToNode.data() + s.numTotalDofs),
                     -1);
        int nodeBlocks = int((s.nodeCount + s.blockSize - 1) / s.blockSize);
        buildDofToNodeKernel<<<nodeBlocks, s.blockSize>>>(
            s.d_node_to_dof.data(), s.d_dofToNode.data(),
            s.numTotalDofs, s.nodeCount);
        cudaDeviceSynchronize();
    }

    // Build the per-NODE velocity-Dirichlet mask + halo-exchange it so the
    // symmetric column-zero step on s.d_valuesVel below hits ghost copies of
    // Dirichlet nodes on non-owning ranks. Mirrors the DDT pattern around
    // line ~3429-3460 (d_isPressureBdryNode). Re-uses the generic
    // buildPressureBdryNodeMaskKernel: it takes any per-DOF mask + optional
    // pin and writes the per-NODE result.
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Channel
        || s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Pump)
    {
        s.d_isBdryNode.resize(s.nodeCount);
        thrust::fill(thrust::device_pointer_cast(s.d_isBdryNode.data()),
                     thrust::device_pointer_cast(s.d_isBdryNode.data() + s.nodeCount),
                     uint8_t(0));
        {
            int nodeBlocks = int((s.nodeCount + s.blockSize - 1) / s.blockSize);
            const uint8_t* maskPtr = (s.d_isBdryDof.size() > 0)
                                     ? s.d_isBdryDof.data() : nullptr;
            buildPressureBdryNodeMaskKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                s.d_node_to_dof.data(), d_nodeOwnership.data(),
                maskPtr, /*pinDof=*/-1,
                s.d_isBdryNode.data(), s.nodeCount, s.numOwnedDofs);
            cudaDeviceSynchronize();
        }
        // Forward halo via RealType proxy (cstone exchangeNodeHalo is locked
        // to the domain's RealType). Proxy values are 0.0/1.0; cast back via
        // >0.5 threshold. Same pattern as DDT pressure mask exchange.
        {
            cstone::DeviceVector<RealType> d_maskProxy(s.nodeCount, RealType(0));
            thrust::transform(thrust::device,
                              thrust::device_pointer_cast(s.d_isBdryNode.data()),
                              thrust::device_pointer_cast(s.d_isBdryNode.data() + s.nodeCount),
                              thrust::device_pointer_cast(d_maskProxy.data()),
                              [] __device__ (uint8_t v) -> RealType { return v ? RealType(1) : RealType(0); });
            s.domain.exchangeNodeHalo(d_maskProxy);
            thrust::transform(thrust::device,
                              thrust::device_pointer_cast(d_maskProxy.data()),
                              thrust::device_pointer_cast(d_maskProxy.data() + s.nodeCount),
                              thrust::device_pointer_cast(s.d_isBdryNode.data()),
                              [] __device__ (RealType v) -> uint8_t { return (v > RealType(0.5)) ? uint8_t(1) : uint8_t(0); });
            cudaDeviceSynchronize();
        }
        pt.lap("velocity per-node BC mask + halo");
    }

    // Fix Y: ship each cross-rank slave's full assembled velocity row to its
    // master-owner rank and accumulate into the master row. Runs BEFORE the
    // row-zero step below so the slave row still carries its original
    // off-diagonals. After this call, the master row on rank D holds the
    // fully merged equation; the slave row will become a Dirichlet identity
    // a few lines later (enforceBcMatrixKernel + d_isPeriodicXRSlaveDof).
    // Fix Y (cross-rank velocity row sum) is intentionally DISABLED.
    // The DOF-collapse Pass 1 above now redirects ALL slave nodeToDof entries
    // to their master's DOF, including cross-rank pairs where the master is a
    // ghost on this rank. After this, the slave's "row" doesn't exist as a
    // separate owned DOF -- the assembler scatters straight into the master's
    // row via the ghost DOF (cstone's periodic-image halo delivers the rank-A
    // slave-side elements onto rank D, where the master is owned). No row-sum
    // exchange needed; same flow as single-rank periodic.

    // Path-B mask: flag owned slave DOFs of cross-rank periodic pairs. The
    // sparsity builder didn't emit slave_row<->master_ghost_col edges on
    // slave-owner ranks (MARS_PERIODIC_XR_SYMCHECK proved A_local=0), so the
    // assembled slave row is structurally non-symmetric. Make it a Dirichlet
    // identity row instead; the post-solve crossRankPeriodicBroadcastDof
    // restores x[slave]:=x[master_owned_D] via MPI.
    s.d_isPeriodicXRSlaveDof.resize(s.numOwnedDofs);
    thrust::fill(thrust::device_pointer_cast(s.d_isPeriodicXRSlaveDof.data()),
                 thrust::device_pointer_cast(s.d_isPeriodicXRSlaveDof.data() + s.numOwnedDofs),
                 uint8_t(0));
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
        && s.periodicMap != nullptr
        && !s.periodicMap->cross_.d_sendOwnedSlaveIds_.empty())
    {
        int nSend = int(s.periodicMap->cross_.d_sendOwnedSlaveIds_.size());
        int blk = (nSend + s.blockSize - 1) / s.blockSize;
        markPeriodicXRSlaveDofKernel<<<blk, s.blockSize>>>(
            s.periodicMap->cross_.d_sendOwnedSlaveIds_.data(),
            s.d_node_to_dof.data(),
            nSend, s.numOwnedDofs,
            s.d_isPeriodicXRSlaveDof.data());
        cudaDeviceSynchronize();
    }

    // Per-NODE Path-B mask for symmetric column enforcement. The d_isPeriodicXRSlaveDof
    // mask above flags OWNED slave DOFs; this per-NODE version is halo-exchanged so
    // ghost-copies of slave nodes on neighboring ranks also carry the flag. That lets
    // enforceBcColMatrixKernelByNode zero column entries on EVERY rank that has an
    // edge to a slave column. Same pattern as d_isPressureBdryNode at line ~3302.
    s.d_isPeriodicXRSlaveNode.resize(s.nodeCount);
    thrust::fill(thrust::device_pointer_cast(s.d_isPeriodicXRSlaveNode.data()),
                 thrust::device_pointer_cast(s.d_isPeriodicXRSlaveNode.data() + s.nodeCount),
                 uint8_t(0));
    // Halo-exchange ALWAYS runs on every rank in Periodic mode, regardless of
    // whether THIS rank has local slaves: exchangeNodeHalo is an MPI collective.
    // Earlier version gated the entire block on !cross_.d_sendOwnedSlaveIds_.empty()
    // -- rank 0 (master-owner, no local slaves) skipped while ranks 1/2/3 entered
    // and called exchangeNodeHalo -> classic MPI deadlock. The scatter that
    // populates the mask is INSIDE its own !empty() guard so the H2D side of
    // the scatter doesn't index into a zero-sized array.
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
        && s.periodicMap != nullptr)
    {
        if (!s.periodicMap->cross_.d_sendOwnedSlaveIds_.empty())
        {
            int nSend = int(s.periodicMap->cross_.d_sendOwnedSlaveIds_.size());
            // Tiny scatter kernel inlined as a thrust::for_each: for each entry in
            // d_sendOwnedSlaveIds_, set d_isPeriodicXRSlaveNode[node] = 1.
            const int* d_sids = s.periodicMap->cross_.d_sendOwnedSlaveIds_.data();
            uint8_t* d_mask = s.d_isPeriodicXRSlaveNode.data();
            thrust::for_each(thrust::device,
                             thrust::counting_iterator<int>(0),
                             thrust::counting_iterator<int>(nSend),
                             [d_sids, d_mask] __device__ (int i) {
                                 d_mask[d_sids[i]] = 1;
                             });
            cudaDeviceSynchronize();
        }

        // Halo-exchange the mask via a RealType proxy (cstone exchangeNodeHalo is
        // locked to the domain's RealType). Proxy values are 0.0/1.0; threshold
        // 0.5 on the way back. Identical pattern to the pressure-DDT path. Runs
        // on ALL ranks (rank 0 will receive 1s in ghost slots from rank 1/2/3
        // slave-owners; rank 0's interior rows that reference those ghost columns
        // then get column-zeroed correctly).
        if (s.numRanks > 1)
        {
            cstone::DeviceVector<RealType> d_maskProxy(s.nodeCount, RealType(0));
            thrust::transform(thrust::device,
                              thrust::device_pointer_cast(s.d_isPeriodicXRSlaveNode.data()),
                              thrust::device_pointer_cast(s.d_isPeriodicXRSlaveNode.data() + s.nodeCount),
                              thrust::device_pointer_cast(d_maskProxy.data()),
                              [] __device__ (uint8_t v) -> RealType { return v ? RealType(1) : RealType(0); });
            s.domain.exchangeNodeHalo(d_maskProxy);
            thrust::transform(thrust::device,
                              thrust::device_pointer_cast(d_maskProxy.data()),
                              thrust::device_pointer_cast(d_maskProxy.data() + s.nodeCount),
                              thrust::device_pointer_cast(s.d_isPeriodicXRSlaveNode.data()),
                              [] __device__ (RealType v) -> uint8_t { return (v > RealType(0.5)) ? uint8_t(1) : uint8_t(0); });
            cudaDeviceSynchronize();
        }
    }

    // Probe 1: per-rank count of marked slave DOFs. Should be ~equal to
    // d_sendOwnedSlaveIds_.size() on slave-owner ranks, 0 on master-owner
    // (e.g. rank 0). Gated by MARS_PERIODIC_PATHB_DBG=1.
    if (std::getenv("MARS_PERIODIC_PATHB_DBG") != nullptr)
    {
        int sendListSize = s.periodicMap
            ? int(s.periodicMap->cross_.d_sendOwnedSlaveIds_.size()) : 0;
        long long maskOnes = thrust::count(
            thrust::device,
            thrust::device_pointer_cast(s.d_isPeriodicXRSlaveDof.data()),
            thrust::device_pointer_cast(s.d_isPeriodicXRSlaveDof.data() + s.numOwnedDofs),
            uint8_t(1));
        // Probe 1b: per-NODE mask count (post halo-exchange). Rank 0 has no
        // local slaves but should now carry GHOST flags from peer slave-owners,
        // so nodeMaskOnes > 0 on rank 0 too -- proves the halo exchange landed.
        long long nodeMaskOnes = thrust::count(
            thrust::device,
            thrust::device_pointer_cast(s.d_isPeriodicXRSlaveNode.data()),
            thrust::device_pointer_cast(s.d_isPeriodicXRSlaveNode.data() + s.nodeCount),
            uint8_t(1));
        std::cerr << "[pathb-dbg rank " << s.rank
                  << "] mask_ones=" << maskOnes
                  << " send_list_size=" << sendListSize
                  << " numOwnedDofs=" << s.numOwnedDofs
                  << " nodeMaskOnes=" << nodeMaskOnes
                  << " nodeCount=" << s.nodeCount << std::endl;
    }

    // Enforce velocity Dirichlet on the velocity matrix (row=0, diag=1).
    // - Cavity/channel/pump: zero rows + diag=1 for d_isBdryDof slots.
    // - Periodic: zero rows + diag=1 for d_isPeriodicXRSlaveDof slots
    //   (Path B: slave row becomes trivial identity, master row carries the
    //   physics; broadcast + post-solve patch handle the slave value).
    {
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        if (s.bcKind != NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic)
        {
            enforceBcMatrixKernel<RealType><<<dofBlocks, s.blockSize>>>(
                s.d_isBdryDof.data(), s.d_rowPtr.data(), s.d_colInd.data(),
                s.d_diagPtr.data(), s.d_valuesVel.data(), s.numOwnedDofs);
            if (s.useBdf2 && s.d_valuesVel_bdf2.size() > 0)
            {
                enforceBcMatrixKernel<RealType><<<dofBlocks, s.blockSize>>>(
                    s.d_isBdryDof.data(), s.d_rowPtr.data(), s.d_colInd.data(),
                    s.d_diagPtr.data(), s.d_valuesVel_bdf2.data(), s.numOwnedDofs);
            }
            // Symmetric column-zero for Pump/Channel. The cstone CG used for
            // the velocity solve is textbook symmetric PCG (alpha=rho/pAp);
            // the row-only clear above leaves the velocity matrix asymmetric.
            // With many Dirichlet DOFs (~125k + 5-OOM diag spread on
            // anisotropic boundary-layer cells) the asymmetric perturbation
            // flips pAp sign and breaks SPD: cg_iter_uvw=FAIL on any non-
            // trivial RHS, even at amplitude=1e-5 (validated this session).
            // Cavity is intentionally skipped: its Dirichlet set is sparse
            // enough that the existing asymmetric row-clear converges fine
            // (same rationale as DDT-cavity at lines ~3500-3508).
            //
            // NOTE: no Dirichlet "lift" yet. For zero qTarget (wall no-slip)
            // it isn't needed. For non-zero qTarget (inlet/extra u=Uinf)
            // this introduces an O(off-diag * Uinf) perturbation localized to
            // the 1-cell ring around those faces. Acceptable for the "BC
            // actually runs" milestone; promote to a proper lift kernel as
            // a follow-up after validating convergence.
            if ((s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Channel
                 || s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Pump)
                && s.d_isBdryNode.size() == static_cast<size_t>(s.nodeCount))
            {
                enforceBcColMatrixKernelByNode<RealType><<<dofBlocks, s.blockSize>>>(
                    s.d_isBdryNode.data(),
                    s.d_node_to_dof.data(), s.d_dofToNode.data(), s.numTotalDofs,
                    s.d_rowPtr.data(), s.d_colInd.data(),
                    s.d_valuesVel.data(), s.numOwnedDofs);
                if (s.useBdf2 && s.d_valuesVel_bdf2.size() > 0)
                {
                    enforceBcColMatrixKernelByNode<RealType><<<dofBlocks, s.blockSize>>>(
                        s.d_isBdryNode.data(),
                        s.d_node_to_dof.data(), s.d_dofToNode.data(), s.numTotalDofs,
                        s.d_rowPtr.data(), s.d_colInd.data(),
                        s.d_valuesVel_bdf2.data(), s.numOwnedDofs);
                }
                cudaDeviceSynchronize();
            }
        }
        else if (s.periodicMap != nullptr)
        {
            // Periodic, conditional-collapse design (pump-pattern): cross-rank
            // slaves keep their own owned DOF. Their rows are real physics
            // rows assembled in full from their owner's xmax-side elements.
            // DO NOT row-zero them -- that would identity-Dirichlet the
            // slave row to u[slave]=0 and kill the field.
            //
            // Slave-master identity is enforced at FIELD EXCHANGE time:
            // crossRankPeriodicBroadcastDof post-solve copies the master's
            // owned value into the slave's owned slot, and standard halo
            // exchange syncs the ghost copies.
        }
        cudaDeviceSynchronize();
    }
    pt.lap("velocity matrix BC");

    // Probe 2: confirm enforceBcMatrixKernel actually zeroed the slave row
    // and set diag=1 on Avel. Read back one cross-rank slave row's diagonal
    // value and the count of nonzero off-diagonal entries on that row.
    if (std::getenv("MARS_PERIODIC_PATHB_DBG") != nullptr
        && s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
        && s.periodicMap != nullptr
        && !s.periodicMap->cross_.d_sendOwnedSlaveIds_.empty())
    {
        const auto& xr = s.periodicMap->cross_;
        int firstSlaveNode = -1;
        {
            std::vector<int> tmp(1);
            cudaMemcpy(tmp.data(), xr.d_sendOwnedSlaveIds_.data(),
                       sizeof(int), cudaMemcpyDeviceToHost);
            firstSlaveNode = tmp[0];
        }
        int slaveDof = -1;
        if (firstSlaveNode >= 0 && firstSlaveNode < int(s.nodeCount))
        {
            std::vector<int> tmp(1);
            cudaMemcpy(tmp.data(),
                       s.d_node_to_dof.data() + firstSlaveNode,
                       sizeof(int), cudaMemcpyDeviceToHost);
            slaveDof = tmp[0];
        }
        if (slaveDof >= 0 && slaveDof < s.numOwnedDofs)
        {
            // Read this row of Avel from CSR.
            std::vector<int> rowOffsets(2);
            cudaMemcpy(rowOffsets.data(),
                       s.d_rowPtr.data() + slaveDof,
                       2 * sizeof(int), cudaMemcpyDeviceToHost);
            int rs = rowOffsets[0], re = rowOffsets[1];
            int nnzInRow = re - rs;
            std::vector<int>      rowCols(nnzInRow);
            std::vector<RealType> rowVals(nnzInRow);
            cudaMemcpy(rowCols.data(),
                       s.d_colInd.data() + rs,
                       nnzInRow * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(rowVals.data(),
                       s.d_valuesVel.data() + rs,
                       nnzInRow * sizeof(RealType), cudaMemcpyDeviceToHost);
            int nonzeroOff = 0;
            RealType diagVal = RealType(0), maxOff = RealType(0);
            for (int j = 0; j < nnzInRow; ++j)
            {
                if (rowCols[j] == slaveDof) diagVal = rowVals[j];
                else if (rowVals[j] != RealType(0))
                {
                    ++nonzeroOff;
                    if (std::abs(rowVals[j]) > std::abs(maxOff)) maxOff = rowVals[j];
                }
            }
            std::cerr << "[pathb-dbg rank " << s.rank
                      << "] firstSlaveDof=" << slaveDof
                      << " rowNnz=" << nnzInRow
                      << " diag=" << diagVal
                      << " nonzeroOff=" << nonzeroOff
                      << " maxOffVal=" << maxOff
                      << " (expect: diag=1, nonzeroOff=0 after Path-B mask)"
                      << std::endl;
        }
    }

    // Optional one-shot cross-rank periodic matrix-symmetry diagnostic.
    // For up to N cross-rank pairs, prints A[slave_row, master_ghost_col] on
    // the slave-owner rank vs A[master_row, slave_ghost_col] on the master
    // owner rank. Equal -> matrix is symmetric across the seam. Different
    // -> assembler is producing asymmetric cross-rank coupling. Gated by
    // MARS_PERIODIC_XR_SYMCHECK=1.
    if (std::getenv("MARS_PERIODIC_XR_SYMCHECK") != nullptr
        && s.numRanks > 1
        && s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
        && s.periodicMap != nullptr
        && !s.periodicMap->cross_.peers_.empty())
    {
        const auto& xr   = s.periodicMap->cross_;
        const int numPeers = int(xr.peers_.size());
        const int sendTotal = xr.sendOffsets_.empty() ? 0 : xr.sendOffsets_.back();
        const int recvTotal = xr.recvOffsets_.empty() ? 0 : xr.recvOffsets_.back();

        // Host copies of CSR + DOF map + partner.
        std::vector<int>      h_rowPtr(s.numTotalDofs + 1);
        std::vector<int>      h_colInd(s.nnz);
        std::vector<RealType> h_valsVel(s.nnz);
        std::vector<int>      h_n2d(s.nodeCount);
        std::vector<int>      h_partner(s.nodeCount);
        std::vector<int>      h_sendIds(sendTotal);
        std::vector<int>      h_recvIds(recvTotal);
        cudaMemcpy(h_rowPtr.data(), s.d_rowPtr.data(), (s.numTotalDofs+1)*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_colInd.data(), s.d_colInd.data(), s.nnz*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_valsVel.data(), s.d_valuesVel.data(), s.nnz*sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_n2d.data(), s.d_node_to_dof.data(), s.nodeCount*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_partner.data(), s.periodicMap->d_periodicPartner.data(), s.nodeCount*sizeof(int), cudaMemcpyDeviceToHost);
        if (sendTotal > 0)
            cudaMemcpy(h_sendIds.data(), xr.d_sendOwnedSlaveIds_.data(), sendTotal*sizeof(int), cudaMemcpyDeviceToHost);
        if (recvTotal > 0)
            cudaMemcpy(h_recvIds.data(), xr.d_recvOwnedMasterIds_.data(), recvTotal*sizeof(int), cudaMemcpyDeviceToHost);

        // For each peer bucket: we own slaves (sendIds), peer owns the master.
        // We exchange the off-diagonal entry value A[slave_row, master_ghost_col]
        // with the peer, who looks up A[master_row, slave_ghost_col] from its
        // own CSR. Use point-to-point pair-wise exchange via the xr.comm_.
        const int N_PROBE = 4;
        auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;

        for (int pi = 0; pi < numPeers; ++pi)
        {
            int peer = xr.peers_[pi];
            int sCnt = xr.sendOffsets_[pi+1] - xr.sendOffsets_[pi];
            int rCnt = xr.recvOffsets_[pi+1] - xr.recvOffsets_[pi];
            int nProbe = std::min(N_PROBE, std::max(sCnt, rCnt));

            // Slave-side: build our list of (our-A_value-at-slave-row-master-col).
            std::vector<RealType> sendVals(nProbe, RealType(0));
            for (int k = 0; k < std::min(nProbe, sCnt); ++k)
            {
                int slave_node = h_sendIds[xr.sendOffsets_[pi] + k];
                int master_ghost_node = h_partner[slave_node];
                int slave_dof = h_n2d[slave_node];
                int master_ghost_dof = (master_ghost_node >= 0)
                                       ? h_n2d[master_ghost_node] : -2;
                if (slave_dof < 0 || slave_dof >= s.numOwnedDofs || master_ghost_dof < 0) continue;
                int rs = h_rowPtr[slave_dof];
                int re = h_rowPtr[slave_dof + 1];
                for (int j = rs; j < re; ++j)
                {
                    if (h_colInd[j] == master_ghost_dof)
                    {
                        sendVals[k] = h_valsVel[j];
                        break;
                    }
                }
            }
            // Master-side: build our list of (our-A_value-at-master-row-slave-col)
            // where the slave is the one peer (sender) will be asking about.
            // We don't know which OUR-owned-master corresponds to which peer-slave
            // by k, so we send both halves and let each side print.
            std::vector<RealType> masterVals(nProbe, RealType(0));
            for (int k = 0; k < std::min(nProbe, rCnt); ++k)
            {
                int master_node = h_recvIds[xr.recvOffsets_[pi] + k];
                int master_dof = h_n2d[master_node];
                if (master_dof < 0 || master_dof >= s.numOwnedDofs) continue;
                // Walk the master row; pick the off-diagonal entry whose column is
                // a GHOST (>= numOwnedDofs) -- that's the slave ghost on this rank.
                int rs = h_rowPtr[master_dof];
                int re = h_rowPtr[master_dof + 1];
                for (int j = rs; j < re; ++j)
                {
                    int col = h_colInd[j];
                    if (col >= s.numOwnedDofs)
                    {
                        masterVals[k] = h_valsVel[j];
                        break;
                    }
                }
            }
            // Send our slave-side values; receive peer's master-side values for our slaves.
            std::vector<RealType> peerMasterVals(nProbe, RealType(0));
            MPI_Sendrecv(masterVals.data(), nProbe, mpiType, peer, 0xBEEF,
                         peerMasterVals.data(), nProbe, mpiType, peer, 0xBEEF,
                         xr.comm_, MPI_STATUS_IGNORE);
            // Print on the SLAVE-owner side (where sCnt > 0). Rank 0 typically
            // owns masters and has sCnt=0; the slave-owning ranks are the
            // meaningful side for this check.
            if (sCnt > 0)
            {
                std::cerr << "[xr-symcheck rank " << s.rank
                          << " (slave-owner) <-> peer " << peer
                          << " (master-owner)] sCnt=" << sCnt
                          << " rCnt=" << rCnt
                          << " nProbe=" << nProbe << std::endl;
                int nIter = std::min(nProbe, sCnt);
                for (int k = 0; k < nIter; ++k)
                {
                    int slave_node = h_sendIds[xr.sendOffsets_[pi] + k];
                    int master_ghost_node = h_partner[slave_node];
                    int slave_dof = h_n2d[slave_node];
                    int master_ghost_dof = (master_ghost_node >= 0)
                                           ? h_n2d[master_ghost_node] : -2;
                    std::cerr << "  k=" << k
                              << "  slave_node=" << slave_node
                              << " (dof=" << slave_dof << ")"
                              << "  master_ghost=" << master_ghost_node
                              << " (dof=" << master_ghost_dof << ")"
                              << "  A_local=" << sendVals[k]
                              << "  A_peer=" << peerMasterVals[k]
                              << "  diff=" << (sendVals[k] - peerMasterVals[k])
                              << std::endl;
                }
            }
        }
        MPI_Barrier(xr.comm_);
    }

    // Pressure null-space removal.
    //   Cavity: pin a single owned DOF closest to (xmin, ymin, zmin) and
    //   leave Neumann everywhere else. Works because pressure is free to be
    //   zero anywhere with pure-Dirichlet velocity BCs.
    //   Channel: pin the ENTIRE outflow face (x = xmax) to p = 0 by adding a
    //   per-owned-DOF mask, then turning matrix rows for those DOFs into
    //   identity. This matches the physically correct outflow condition for
    //   channel flow and avoids the single-corner-pin pathology that fought
    //   the natural inflow-to-outflow pressure gradient.
    using BCK = typename NSStepper<KeyType, RealType, ElementTag>::BCKind;
    if (s.bcKind == BCK::Periodic)
    {
        // Periodic: pressure is pure-Neumann with a 1D constant-mode null space.
        // Default: leave it Neumann and clean up with per-step removeMean.
        //
        // MARS_PERIODIC_PIN: pin ONE owned DOF (nearest the box corner) to break
        // the null space explicitly, exactly like the cavity. On >1 rank the
        // Neumann+removeMean route was observed to break down (the assembled
        // operator's numerical null vector is not exactly the constant across the
        // cross-rank seam, so the mean projection is inconsistent -> CG search
        // direction lands in the near-null space -> pAp<=0). A hard pin removes
        // the null direction outright and does not depend on seam completeness.
        // Same mechanism as the cavity pin below (findPressurePinCandidateKernel
        // + global MINLOC + enforcePinRowMatrixKernel).
        const char* pinEv = std::getenv("MARS_PERIODIC_PIN");
        bool periodicPin = (pinEv && std::string(pinEv) != "0");
        if (!periodicPin)
        {
            s.pressurePinDof  = -1;
            s.pressurePinRank = -1;
            if (s.rank == 0)
                std::cout << "  pressure BC: periodic (no pin; null-space removed each step via removeMean)\n";
        }
        else
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
                std::cout << "  pressure BC: periodic + single-DOF pin (MARS_PERIODIC_PIN): rank="
                          << s.pressurePinRank << ", dof="
                          << ((s.rank == s.pressurePinRank) ? s.pressurePinDof : -1)
                          << " (anchor near corner (" << s.xmin << "," << s.ymin << "," << s.zmin << "))\n";
        }
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
            s.nodeCount, s.xmax, s.bboxEps, s.numOwnedDofs);
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
    else if (s.bcKind == BCK::Pump)
    {
        // Pump: pressure Dirichlet p=0 on the 'out' side-set (outlet face).
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
        // When the outlet is a VELOCITY-Dirichlet outflow (mass-conserving pump,
        // outletU>0), prescribing p over the WHOLE outlet face too would
        // over-constrain it. The pressure then needs only ONE reference DOF, so
        // mask a single outlet DOF (the first owned one, on the lowest rank that
        // has one). When outletU<=0 (legacy pressure-outlet), mask the whole
        // outlet face as before (natural-Neumann velocity there).
        bool singlePin = (s.outletU > RealType(0));
        int  pinRankLocal = singlePin ? s.numRanks : -1;  // for the global argmin
        size_t ownedOutletCount = 0;
        int firstOutletDof = -1;
        for (int li : s.outletNodes)
        {
            if (li < 0 || (size_t)li >= s.nodeCount) continue;
            if (hostOwn[li] != 1) continue;
            int dof = hostNodeToDof[li];
            if (dof < 0 || dof >= s.numOwnedDofs) continue;
            if (firstOutletDof < 0) firstOutletDof = dof;
            if (!singlePin) { hostMask[dof] = 1; ++ownedOutletCount; }
        }
        if (singlePin)
        {
            // Pick the single global pin: lowest rank that owns an outlet DOF.
            int haveOutlet = (firstOutletDof >= 0) ? s.rank : s.numRanks;
            int pinRank = s.numRanks;
            MPI_Allreduce(&haveOutlet, &pinRank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            if (s.rank == pinRank && firstOutletDof >= 0)
            {
                hostMask[firstOutletDof] = 1;
                ownedOutletCount = 1;
            }
            (void)pinRankLocal;
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
            std::cout << "  pressure BC: " << (singlePin
                         ? "single p=0 pin (mass-conserving velocity outlet)"
                         : "Dirichlet p=0 on whole outlet side-set")
                      << " (" << ownedOutletCount << " masked DOFs on rank 0)\n";
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

    // These enforce Dirichlet BCs on the ASSEMBLED DDT matrix (s.d_*DDT). Those
    // buffers are only built inside the if constexpr(hex) DDT block above, so
    // nnzDDT==0 for tet / any run without an assembled DDT operator -- skip,
    // or the kernels dereference a null DDT rowPtr (illegal access). The K-path
    // velocity/pressure matrices get their own BC enforcement separately.
    if (s.nnzDDT > 0
        && (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Channel
            || s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Pump))
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
    else if (s.nnzDDT > 0
             && s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Cavity
             && s.pressurePinDof >= 0)
    {
        // Row enforcement on pin row, KEEPING the assembled diagonal magnitude
        // (not forcing it to 1). The DDT operator's diagonal range is ~0.04 to
        // 0.14; a pin row with diag=1 is 10-25x its neighbors, which CG
        // tolerates but BoomerAMG's coarsening does not -- Hypre setup returns
        // generic error 1 on that conditioning spike. Preserving the native
        // diagonal removes the spike so AMG sees a uniformly-scaled SPD row.
        enforcePinRowKeepDiagKernel<RealType><<<1, 1>>>(
            s.pressurePinDof, s.d_rowPtrDDT.data(), s.d_colIndDDT.data(),
            s.d_diagPtrDDT.data(), s.d_valuesDDT.data());
    }
    // Cavity DDT: column-clear stays disabled. A symmetric col-clear changed
    // neighboring rows' row-sum from ~0 to ~|A[r,pin]|, breaking AMG's
    // strength-of-connection metric (earlier attempt -> Hypre setup error 1).
    // K-path doesn't col-clear for cavity and works fine; matching that
    // asymmetric-pin treatment for DDT-cavity, now with a scale-matched pin
    // diagonal so AMG can coarsen the resulting SPD row uniformly. For
    // channel/pump the col-clear above is kept since those have an entire
    // Dirichlet face, not a single pin, and the row-sum drift is negligible
    // relative to the boundary mass.

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
        addBochevDohrmannStab<KeyType, RealType, ElementTag>(s, tau);
        if (s.rank == 0)
            std::cout << "  pressure stabilization: Bochev-Dohrmann polynomial projection (K-path matrix + DDT matrix-free), tau=" << tau << "\n";
        pt.lap("BD pressure stabilization");
    }

    // Cache the DDT diagonal (BC-modified, post-shift) for use as the Jacobi
    // preconditioner in solvePressureDDT. Single O(numOwnedDofs) kernel; the
    // diagonal stays valid for the whole run because BC enforcement and the
    // optional shift are setup-time-only modifications of s.d_valuesDDT.
    // Note: BD-periodic adds matrix-free terms inside applyDDTPerNode that are
    // NOT reflected in s.d_valuesDDT, so this cached diagonal will be slightly
    // off the actual operator diagonal on periodic runs. Acceptable: Jacobi
    // is a preconditioner, not the exact inverse; a small mismatch only slows
    // convergence, not correctness.
    if (s.numOwnedDofs > 0 && s.nnzDDT > 0)
    {
        s.d_diagDDT.resize(s.numOwnedDofs);
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        extractDiagByDiagPtrKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_diagPtrDDT.data(), s.d_valuesDDT.data(),
            s.d_diagDDT.data(), s.numOwnedDofs);
        cudaDeviceSynchronize();
        // Cavity-only override: the cavity pin is enforced in the matrix-free
        // applyDDTPerNode as an identity row (out[pin]=phi[pin], effective
        // diag=1.0), but the assembled s.d_valuesDDT keeps the pin's native
        // diagonal (~0.04-0.14) via enforcePinRowKeepDiagKernel. Without this
        // override the Jacobi step would compute z[pin] = r[pin]/0.04 = ~25x
        // amplification at the pin DOF every iter -- biasing the global
        // search direction and producing the |gphi|~5x divergence injection
        // observed earlier in this session before the loop-bound fix landed.
        // Channel/pump already get diag=1 from enforceBcMatrixKernel so they
        // need no override. Predicted Jacobi-PCG pump iter count: ~400-700
        // to 1e-8 (vs un-precond stall at 1e-3 after 20k).
        if (s.pressurePinDof >= 0
            && s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Cavity)
        {
            RealType one = RealType(1);
            cudaMemcpy(s.d_diagDDT.data() + s.pressurePinDof, &one,
                       sizeof(RealType), cudaMemcpyHostToDevice);
        }
    }
    else if (s.numOwnedDofs > 0 && s.pressureSolve == PressureSolveKind::DDT
             && std::is_same_v<ElementTag, TetTag>)
    {
        // Tet has no assembled DDT matrix (nnzDDT==0), so build the Jacobi
        // diagonal matrix-free from the SCS area vectors + node lumped mass,
        // matching the matrix-free applyDDTPerNode operator. Without this the
        // tet DDT CG runs unpreconditioned and stalls at the iter cap.
        //
        // HARD DEPENDENCY: this build reads s.d_areaVec_{x,y,z} (filled by
        // precomputeTetAreaVectorsGpu) and s.d_massNode (filled by the lumped-
        // mass step). Both run earlier in setupNSStepper, so they are populated
        // here. If a future refactor moves this block ahead of either, the
        // accumulator comes back all-zero -> floored to a constant 1.0 -> the
        // Jacobi preconditioner becomes a no-op (useJacobi stays true but CG
        // runs effectively un-preconditioned). Catch that at the source.
        if (s.d_areaVec_x.size() != s.elementCount * ElemTraits<ElementTag>::ScsPerElem
            || s.d_massNode.size() != s.nodeCount)
        {
            if (s.rank == 0)
                std::cout << "[cg-ddt] WARNING: tet DDT diagonal build ran before its "
                             "inputs were ready (areaVec/massNode unsized) -> Jacobi "
                             "preconditioner will be a useless constant. Move this block "
                             "AFTER the area-vector + lumped-mass steps.\n";
        }
        s.d_diagDDT.resize(s.numOwnedDofs);
        thrust::fill(thrust::device_pointer_cast(s.d_diagDDT.data()),
                     thrust::device_pointer_cast(s.d_diagDDT.data() + s.numOwnedDofs),
                     RealType(0));
        cstone::DeviceVector<RealType> d_diagAccNode(s.nodeCount, RealType(0));

        auto cp = connPtrs<ElementTag, KeyType>(d_conn);
        const KeyType* c0 = cp[0]; const KeyType* c1 = cp[1];
        const KeyType* c2 = cp[2]; const KeyType* c3 = cp[3];
        const size_t startElem = s.domain.startIndex();
        const size_t numLocal  = s.domain.localElementCount();
        const int eBlocks = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
        if (eBlocks > 0)
        {
            computeTetDDTDiagonalKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, nullptr, nullptr, nullptr, nullptr,
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                s.d_massNode.data(), d_diagAccNode.data(), startElem, numLocal);
            cudaDeviceSynchronize();
        }
        // Sum cross-rank / periodic incident-face contributions exactly as the
        // operator does for its accumulators.
        s.domain.reverseExchangeNodeHaloAdd(d_diagAccNode);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_diagAccNode);

        int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        gatherTetDDTDiagToDofKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_diagAccNode.data(), s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_diagDDT.data(), s.nodeCount, s.numOwnedDofs);
        cudaDeviceSynchronize();
        // Floor any zero/negative diagonal so the Jacobi 1/diag step is safe.
        thrust::transform(thrust::device,
            thrust::device_pointer_cast(s.d_diagDDT.data()),
            thrust::device_pointer_cast(s.d_diagDDT.data() + s.numOwnedDofs),
            thrust::device_pointer_cast(s.d_diagDDT.data()),
            [] __device__ (RealType d) { return (d > RealType(1e-30)) ? d : RealType(1); });
    }

    // Shared clipped-Jacobi epsilon floor (runs for whichever path filled
    // d_diagDDT -- hex assembled or tet matrix-free). Wing meshes with
    // anisotropic boundary-layer cells have a diag spread of 5+ orders of
    // magnitude; without a floor the Jacobi step amplifies tiny-diag rows by
    // ~1/diag, biasing the CG search direction. The floor caps amplification
    // at 1/eps. 1e-3 fraction is the PETSc/Hypre default; env override
    // MARS_DDT_JACOBI_CLIP_FRAC=F sets eps = F * max(diag).
    if (s.d_diagDDT.size() == static_cast<size_t>(s.numOwnedDofs) && s.numOwnedDofs > 0)
    {
        RealType localMax = thrust::reduce(
            thrust::device,
            thrust::device_pointer_cast(s.d_diagDDT.data()),
            thrust::device_pointer_cast(s.d_diagDDT.data() + s.numOwnedDofs),
            RealType(0), thrust::maximum<RealType>());
        RealType globalMax = localMax;
        if (s.numRanks > 1)
        {
            MPI_Datatype mpiR = std::is_same<RealType, double>::value
                                ? MPI_DOUBLE : MPI_FLOAT;
            MPI_Allreduce(&localMax, &globalMax, 1, mpiR, MPI_MAX, MPI_COMM_WORLD);
        }
        RealType frac = RealType(1e-3);
        const char* ev = std::getenv("MARS_DDT_JACOBI_CLIP_FRAC");
        if (ev) { double v = std::atof(ev); if (v > 0) frac = RealType(v); }
        s.diagDDTEpsClip = frac * globalMax;
        if (s.rank == 0)
        {
            auto oldFlags = std::cout.flags();
            std::cout << std::scientific << std::setprecision(4)
                      << "  [DDT-jacobi-clip] eps = " << s.diagDDTEpsClip
                      << " (= " << frac << " * max_diag " << globalMax
                      << "; floors tiny boundary-layer diagonals)\n";
            std::cout.flags(oldFlags);
        }
    }
    pt.lap("DDT diagonal cache (for Jacobi PCG)");

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
    // AddT only exists when the DDT operator was assembled (hex). For tet,
    // nnzDDT==0 and the d_*DDT buffers are null/empty -- wrapping them does a
    // cudaMemcpy on a null pointer (cudaErrorInvalidValue -> poisoned context).
    if (s.nnzDDT > 0)
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
    //
    // Hex-only: applyDDTPerNode is a no-op for tet and the assembled DDT CSR
    // (s.d_rowPtrDDT etc.) is only built in the hex branch above. if constexpr
    // keeps the whole probe out of the tet compile.
    if constexpr (std::is_same_v<ElementTag, HexTag>)
    if (s.bcKind != NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
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
            applyDDTPerNode<KeyType, RealType, ElementTag>(s, phiNode, yMfNode, gx, gy, gz);

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
                                       // instead of u=0 (pump free-stream IC).
                                       RealType interiorU,
                                       RealType interiorV,
                                       RealType interiorW)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    // Default interior = (interiorU, interiorV, interiorW). For cavity that is
    // (0,0,0); for pump flow it is the free-stream (Uinf, 0, 0) so step 1's
    // predictor doesn't see a 1-cell-thick velocity discontinuity at the inlet
    // and extra Dirichlet faces.
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

// Smooth low-mode perturbation on interior (v, w) so the wing system has
// a broken symmetry for body force / advection to amplify. Coordinate-based
// (sinusoidal) rather than per-DOF white noise: white-noise IC fails the
// implicit-diffusion CG on high-aspect-ratio boundary-layer cells because
// the noise wavelength is below the mesh's local resolution in the wall-
// normal direction, producing huge K*u_perturb and an ill-conditioned RHS.
// NekRS's userdat2 examples (eddy_uv, taylor_green) use exactly this form:
// few-mode trigonometric perturbation, eps*Uinf*sin(k*x)... amplitude.
// Owned nodes only; halo resyncs ghosts after.
template<typename RealType>
__global__ void perturbInteriorVWKernel(const uint8_t* isBdryDof,
                                        const int* nodeToDof,
                                        const uint8_t* ownership,
                                        const RealType* nodeX,
                                        const RealType* nodeY,
                                        const RealType* nodeZ,
                                        RealType* v,
                                        RealType* w,
                                        size_t numNodes,
                                        RealType amplitude,
                                        RealType kx, RealType ky, RealType kz,
                                        int numOwnedDofs)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0 || dof >= numOwnedDofs) return;
    if (isBdryDof[dof]) return;
    RealType x = nodeX[i], y = nodeY[i], z = nodeZ[i];
    // v: sin(kx*x) * cos(ky*y) * sin(kz*z) -- divergence-friendly, low-mode
    // w: cos(kx*x) * sin(ky*y) * cos(kz*z) -- distinct phase so v,w differ
    RealType vp = sin(kx * x) * cos(ky * y) * sin(kz * z);
    RealType wp = cos(kx * x) * sin(ky * y) * cos(kz * z);
    v[i] += amplitude * vp;
    w[i] += amplitude * wp;
}

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void applyInitialCondition(NSStepper<KeyType, RealType, ElementTag>& s)
{
    const auto& d_nodeOwnership = s.ownershipMap();
    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    // Interior IC: cavity/channel default to (0,0,0). Pump starts at free-stream
    // (Uinf, 0, 0) so step 1's predictor doesn't see a velocity discontinuity
    // at the inlet/extra Dirichlet boundaries.
    RealType iu = 0, iv = 0, iw = 0;
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Pump)
    {
        iu = s.Uinf;
    }
    initialConditionKernel<RealType><<<nBlocks, s.blockSize>>>(
        s.d_isBdryDof.data(), s.d_node_to_dof.data(), d_nodeOwnership.data(),
        s.d_uTarget.data(), s.d_vTarget.data(), s.d_wTarget.data(),
        s.d_u.data(), s.d_v.data(), s.d_w.data(),
        s.nodeCount, iu, iv, iw);
    // Optional interior perturbation on (v, w). Required for wing/tunnel runs
    // where the IC is uniform free-stream -- a discrete steady state of the
    // predictor unless something breaks the symmetry (NekRS userdat2 pattern).
    if (s.icPerturbMag > RealType(0))
    {
        RealType amp = s.icPerturbMag * s.Uinf;
        // Few-mode trigonometric perturbation (~3 waves across each bbox axis)
        // -- well-resolved on any mesh that resolves the geometry, regardless
        // of anisotropy. Two-pi/(L/3) = 6*pi/L wavenumber.
        const RealType twoPi = RealType(2) * RealType(3.14159265358979323846);
        RealType Lx = std::max(RealType(1e-30), s.xmax - s.xmin);
        RealType Ly = std::max(RealType(1e-30), s.ymax - s.ymin);
        RealType Lz = std::max(RealType(1e-30), s.zmax - s.zmin);
        RealType kx = twoPi * RealType(3) / Lx;
        RealType ky = twoPi * RealType(3) / Ly;
        RealType kz = twoPi * RealType(3) / Lz;
        const auto& d_x = s.domain.getNodeX();
        const auto& d_y = s.domain.getNodeY();
        const auto& d_z = s.domain.getNodeZ();
        perturbInteriorVWKernel<RealType><<<nBlocks, s.blockSize>>>(
            s.d_isBdryDof.data(), s.d_node_to_dof.data(), d_nodeOwnership.data(),
            d_x.data(), d_y.data(), d_z.data(),
            s.d_v.data(), s.d_w.data(), s.nodeCount, amp, kx, ky, kz,
            s.numOwnedDofs);
        if (s.rank == 0)
        {
            // Force scientific 4-sig-figs to override stream state leaked
            // from prior std::fixed/std::setprecision(2) used by pt.lap printer.
            auto oldFlags = std::cout.flags();
            std::cout << std::scientific << std::setprecision(4)
                      << "IC perturbation: (v, w) += eps*Uinf * smooth-mode-3, eps=" << s.icPerturbMag
                      << " Uinf=" << s.Uinf << " amplitude=" << amp
                      << " kx,ky,kz=(" << kx << "," << ky << "," << kz << ")"
                      << "  (NekRS userdat2 trig-mode pattern)\n";
            std::cout.flags(oldFlags);
        }
    }
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
int solveOneComponent(NSStepper<KeyType, RealType, ElementTag>& s,
                      cstone::DeviceVector<RealType>& b_rhs,
                      cstone::DeviceVector<RealType>& xVec,
                      cstone::DeviceVector<RealType>& qOut,
                      typename NSStepper<KeyType, RealType, ElementTag>::Matrix& A,
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
            // Periodic-aware halo exchange. The standard cstone exchange syncs
            // owners->ghosts for ordinary shared nodes. For multi-rank periodic
            // we additionally broadcast master->slave on the search vector p
            // every iteration, otherwise the periodic-image DOF columns that
            // the assembler put in the matrix reference stale values, A*p
            // produces garbage, and pAp<=0 -> CG breakdown -> NaN.
            const int* partnerPtr = (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
                                     && s.periodicMap)
                                    ? s.periodicMap->d_periodicPartner.data()
                                    : nullptr;
            size_t nNodes = s.nodeCount;
            int bs = s.blockSize;
            // Path B: pass the cross-rank-slave skip mask so the per-iter
            // broadcast does NOT overwrite p[slave_owned_dof]. Slave rows are
            // Dirichlet-identity (1 * x[slave] = 0), so we want p[slave] to
            // stay 0 throughout the iteration -- otherwise the broadcast
            // re-introduces the spurious coupling.
            const uint8_t* xrSkipPtr = (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
                                        && !s.d_isPeriodicXRSlaveDof.empty())
                                       ? s.d_isPeriodicXRSlaveDof.data() : nullptr;
            int nOwnedDofs = s.numOwnedDofs;
            // Probe 3: capture firstSlaveDof and a flag for the first 2 callback
            // invocations so we can verify p[firstSlaveDof] stays 0.
            const bool pathbDbg = (std::getenv("MARS_PERIODIC_PATHB_DBG") != nullptr);
            int firstSlaveDof = -1;
            if (pathbDbg && xrSkipPtr != nullptr)
            {
                // Find the first dof with mask==1 on host.
                std::vector<uint8_t> hMask(nOwnedDofs);
                cudaMemcpy(hMask.data(), xrSkipPtr,
                           nOwnedDofs * sizeof(uint8_t),
                           cudaMemcpyDeviceToHost);
                for (int i = 0; i < nOwnedDofs; ++i)
                    if (hMask[i]) { firstSlaveDof = i; break; }
            }
            std::shared_ptr<int> callCounter = std::make_shared<int>(0);
            int rankCapture = s.rank;
            solver.setHaloExchangeCallback(
                [&dom, dofMapPtr, partnerPtr, nNodes, bs, xrSkipPtr, nOwnedDofs,
                 pathbDbg, firstSlaveDof, callCounter, rankCapture]
                (cstone::DeviceVector<RealType>& p)
                {
                    // Probe 3a: read p[firstSlaveDof] BEFORE halo exchange/broadcast
                    RealType pBefore = RealType(0), pAfter = RealType(0);
                    if (pathbDbg && firstSlaveDof >= 0 && *callCounter < 3)
                    {
                        cudaMemcpy(&pBefore,
                                   p.data() + firstSlaveDof,
                                   sizeof(RealType), cudaMemcpyDeviceToHost);
                    }
                    dom.exchangeNodeHalo(p, dofMapPtr);
                    if (partnerPtr)
                    {
                        int grd = int((nNodes + bs - 1) / bs);
                        mars::fem::periodicBroadcastDofKernel<RealType>
                            <<<grd, bs>>>(partnerPtr, dofMapPtr, nNodes, p.data(),
                                          xrSkipPtr, nOwnedDofs);
                        cudaDeviceSynchronize();
                    }
                    // Probe 3b: read p[firstSlaveDof] AFTER halo+broadcast.
                    // If the skip-mask works, pAfter MUST equal pBefore (the
                    // broadcast skipped this slot). If pAfter != pBefore, the
                    // skip-mask is being ignored.
                    if (pathbDbg && firstSlaveDof >= 0 && *callCounter < 3)
                    {
                        cudaMemcpy(&pAfter,
                                   p.data() + firstSlaveDof,
                                   sizeof(RealType), cudaMemcpyDeviceToHost);
                        std::cerr << "[pathb-dbg rank " << rankCapture
                                  << " cgcall=" << *callCounter
                                  << "] firstSlaveDof=" << firstSlaveDof
                                  << " p_before=" << pBefore
                                  << " p_after=" << pAfter
                                  << " (expect equal)" << std::endl;
                    }
                    ++(*callCounter);
                });

            // Path B end-state: slave rows are now Dirichlet-identity rows
            // (A[slave,:]=0, A[slave,slave]=1, b[slave]=0). Master rows carry
            // the full physics via the master-ghost columns the assembler
            // already populated on the master-owner rank. The post-solve
            // crossRankPeriodicBroadcastDof on xVec restores x[slave]=x[master].
            // spmvPostCallback is not needed -- Ap[slave]=p[slave]=0, so the
            // slave's contribution to dot(p,Ap) is bit-exact zero.
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

    // Only scatter a CONVERGED solution into qOut. A failed Hypre solve
    // (error 1 on the near-singular DDT operator) can leave xVec holding
    // partial garbage; scattering it would poison qOut (= phi) and the
    // corrector would turn that into a NaN velocity field, cascading the
    // whole run. On failure we leave qOut untouched (caller warm-started it
    // to its previous value) and signal -2 so the step is recognizably bad.
    if (converged)
    {
        // Path B: before scattering xVec back to per-node qOut, fix up
        // cross-rank slave DOFs by an MPI broadcast from the master's owner
        // rank. Slave rows were Dirichlet-identity with b[slave]=0, so the
        // solve produced x[slave]=0; the broadcast overwrites that with the
        // master's converged value so the scatter writes the correct per-node
        // value. No-op on single-rank / non-periodic / empty cross-rank-pair.
        // Post-solve broadcast disabled: the new cross-rank DOF collapse points
        // slave nodeToDof at the master's ghost DOF, so scatterDofToNodeKernel
        // below writes x[master_ghost_dof] (synced by cstone halo from
        // master_owned on the master rank) into the slave's per-node slot
        // automatically. No MPI patch needed.
        if (false && s.numRanks > 1
            && s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
            && s.periodicMap != nullptr
            && !s.periodicMap->cross_.peers_.empty())
        {
            mars::fem::crossRankPeriodicBroadcastDof<KeyType, RealType>(
                *s.periodicMap, s.d_node_to_dof.data(), xVec);

            // MARS_PERIODIC_PATHB_CHECK=1: assert x[slave_owned_dof_A] equals
            // x[master_owned_dof_D] for up to N probed cross-rank pairs. Both
            // sides should see the same value after the broadcast.
            if (std::getenv("MARS_PERIODIC_PATHB_CHECK") != nullptr)
            {
                const auto& xr = s.periodicMap->cross_;
                const int sendTotal = xr.sendOffsets_.empty() ? 0 : xr.sendOffsets_.back();
                const int recvTotal = xr.recvOffsets_.empty() ? 0 : xr.recvOffsets_.back();
                const int N_PROBE = 4;

                std::vector<int> h_sendIds(sendTotal);
                std::vector<int> h_recvIds(recvTotal);
                std::vector<int> h_n2d(s.nodeCount);
                std::vector<RealType> h_xVec(s.numOwnedDofs);
                if (sendTotal > 0)
                    cudaMemcpy(h_sendIds.data(), xr.d_sendOwnedSlaveIds_.data(),
                               sendTotal*sizeof(int), cudaMemcpyDeviceToHost);
                if (recvTotal > 0)
                    cudaMemcpy(h_recvIds.data(), xr.d_recvOwnedMasterIds_.data(),
                               recvTotal*sizeof(int), cudaMemcpyDeviceToHost);
                cudaMemcpy(h_n2d.data(), s.d_node_to_dof.data(),
                           s.nodeCount*sizeof(int), cudaMemcpyDeviceToHost);
                cudaMemcpy(h_xVec.data(), xVec.data(),
                           s.numOwnedDofs*sizeof(RealType), cudaMemcpyDeviceToHost);

                auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
                for (size_t pi = 0; pi < xr.peers_.size(); ++pi)
                {
                    int peer = xr.peers_[pi];
                    int sCnt = xr.sendOffsets_[pi+1] - xr.sendOffsets_[pi];
                    int rCnt = xr.recvOffsets_[pi+1] - xr.recvOffsets_[pi];
                    int nProbe = std::min(N_PROBE, std::max(sCnt, rCnt));
                    std::vector<RealType> mineMaster(nProbe, RealType(0));
                    for (int k = 0; k < std::min(nProbe, rCnt); ++k)
                    {
                        int node = h_recvIds[xr.recvOffsets_[pi] + k];
                        int dof = h_n2d[node];
                        if (dof >= 0 && dof < s.numOwnedDofs) mineMaster[k] = h_xVec[dof];
                    }
                    std::vector<RealType> peerMaster(nProbe, RealType(0));
                    MPI_Sendrecv(mineMaster.data(), nProbe, mpiType, peer, 0xB0B0,
                                 peerMaster.data(), nProbe, mpiType, peer, 0xB0B0,
                                 xr.comm_, MPI_STATUS_IGNORE);
                    if (sCnt > 0)
                    {
                        int fails = 0;
                        for (int k = 0; k < std::min(nProbe, sCnt); ++k)
                        {
                            int snode = h_sendIds[xr.sendOffsets_[pi] + k];
                            int sdof  = h_n2d[snode];
                            RealType mine = (sdof >= 0 && sdof < s.numOwnedDofs)
                                            ? h_xVec[sdof] : RealType(0);
                            if (mine != peerMaster[k]) ++fails;
                        }
                        std::cerr << "[pathb-check rank " << s.rank
                                  << " <-> peer " << peer
                                  << "] " << (fails == 0 ? "PASS" : "FAIL")
                                  << " (" << fails << "/" << std::min(nProbe, sCnt)
                                  << " slave!=master)" << std::endl;
                    }
                }
                MPI_Barrier(xr.comm_);
            }
        }

        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        scatterDofToNodeKernel<RealType><<<nBlocks, s.blockSize>>>(
            xVec.data(), s.d_node_to_dof.data(), qOut.data(), s.nodeCount);
        cudaDeviceSynchronize();
    }

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
// Default arg for ElementTag is specified ONLY on the forward decl above
// (~line 1661); repeating it here is a C++ "redefinition of default
// argument" error.
template<typename KeyType, typename RealType, typename ElementTag>
void applyDDTPerNode(NSStepper<KeyType, RealType, ElementTag>& s,
                     cstone::DeviceVector<RealType>& phi,
                     cstone::DeviceVector<RealType>& outAcc,
                     cstone::DeviceVector<RealType>& gxAcc,
                     cstone::DeviceVector<RealType>& gyAcc,
                     cstone::DeviceVector<RealType>& gzAcc)
{
    // Matrix-free D M^-1 D^T operator. Element-generic: steps a/b/c go through
    // the templated applyDivTransposePerNodeKernel / computeDivergencePerNodeKernel
    // (which read c0..c3 only for tet), and the halo/periodic/normalize helpers
    // are element-agnostic. connPtrs<ElementTag> sets c4..c7 = nullptr for tet.
    // The only hex-specific code is the Bochev-Dohrmann lambda below (n[8],
    // n[6], c4..c7 deref), which is already dead (if(false)) AND now also
    // if constexpr(hex)-gated so it never compiles for tet.
    const auto& d_nodeOwnership = s.ownershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    // connPtrs gives NodesPerElem pointers; c4..c7 stay nullptr for tet (the
    // templated Phase-1 kernels read c0..c3 only when ElementTag==TetTag).
    auto cp = connPtrs<ElementTag, KeyType>(d_conn);
    const KeyType* c0 = cp[0]; const KeyType* c1 = cp[1];
    const KeyType* c2 = cp[2]; const KeyType* c3 = cp[3];
    const KeyType* c4 = nullptr; const KeyType* c5 = nullptr;
    const KeyType* c6 = nullptr; const KeyType* c7 = nullptr;
    if constexpr (std::is_same_v<ElementTag, HexTag>) { c4 = cp[4]; c5 = cp[5]; c6 = cp[6]; c7 = cp[7]; }
    // OWNED-only per-element scatter. Each rank owns its share of the elements
    // it sees; rank-boundary face contributions are picked up by the OWNER side
    // of the face and then routed to the ghost side via reverseExchangeNodeHaloAdd
    // immediately below. Looping owned+halo here would DOUBLE-COUNT shared faces
    // because the halo element on this rank is the same physical element the
    // owning rank already scatters, and reverseExchangeNodeHaloAdd then folds the
    // ghost-side contribution back. Restored to the proven-green pattern from
    // commit 5da5d6a (regression: div_max=1.65e-5 at cube16 cavity step 500).
    const size_t startElem = s.domain.startIndex();
    const size_t numLocal  = s.domain.localElementCount();
    const int eBlocks    = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
    const int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;

    auto zeroVec = [&](cstone::DeviceVector<RealType>& v) {
        thrust::fill(thrust::device_pointer_cast(v.data()),
                     thrust::device_pointer_cast(v.data() + s.nodeCount), RealType(0));
    };

    // Periodic prolongation P (input side): copy phi[master] -> phi[slave] for
    // every collapsed periodic pair BEFORE D^T reads phi. This is the MFEM/
    // deal.II P^T A P pattern: the operator is A_reduced = P^T A P, where P
    // spreads the single periodic DOF's value to both its node slots and P^T
    // (maybePeriodicSum on outAcc in step c) sums both slots' contributions back
    // onto the one DOF. Applying P here -- inside the matvec, on the operator
    // input -- makes the D^T read transpose-consistent with the D/P^T write, so
    // A is symmetric across the cross-rank seam. Doing this on the OPERATOR input
    // (not by mutating the CG iterate r/p between matvecs) is what keeps the
    // Krylov recurrence conjugate; the per-iteration broadcast it replaces is the
    // documented anti-pattern that stalled CG on >1 rank. No-op for Dirichlet/
    // pump (no periodic map). The slave's own incoming value is overwritten and
    // inert -- it never contributes to a dot product (ownedDot skips slaves).
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic && s.periodicMap)
    {
        const int* d_partner = s.periodicMap->d_periodicPartner.data();
        mars::fem::periodicBroadcastSameRankKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_partner, d_nodeOwnership.data(), s.nodeCount, phi.data());
        cudaDeviceSynchronize();
        if (s.numRanks > 1)
            mars::fem::crossRankPeriodicBroadcast<KeyType, RealType>(*s.periodicMap, phi);
        s.domain.exchangeNodeHalo(phi);
    }

    // Step a: g = D^T phi (un-normalized per-node 3-vector accumulator).
    zeroVec(gxAcc); zeroVec(gyAcc); zeroVec(gzAcc);
    if (eBlocks > 0)
    {
        applyDivTransposePerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
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
    maybePeriodicSum<KeyType, RealType, ElementTag>(s, gxAcc);
    maybePeriodicSum<KeyType, RealType, ElementTag>(s, gyAcc);
    maybePeriodicSum<KeyType, RealType, ElementTag>(s, gzAcc);

    // Step b: g <- M^{-1} g (in-place; gather is same-index, safe).
    normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        gxAcc.data(), gyAcc.data(), gzAcc.data(), s.d_massNode.data(),
        s.d_node_to_dof.data(), d_nodeOwnership.data(),
        gxAcc.data(), gyAcc.data(), gzAcc.data(), s.nodeCount);
    cudaDeviceSynchronize();
    // Forward halo: D-scatter in step c reads g[iL], g[iR] at ghost slots too.
    s.domain.exchangeNodeHalo(gxAcc);
    s.domain.exchangeNodeHalo(gyAcc);
    s.domain.exchangeNodeHalo(gzAcc);

    // Periodic g-consistency: maybePeriodicSum above merged the same-rank slave
    // accumulator onto the master and ZEROED the slave slot; normalize then left
    // g[slave]=0/mass=0 while g[master] carries M^-1 D^T phi. The forward halo
    // does NOT fix this (slave and master have distinct SFC keys). If step c
    // scatters with g[slave]=0 the D operator is ASYMMETRIC at the periodic
    // face and the projection identity div(u^{n+1})=0 never closes there.
    // Broadcast g master->slave (same-rank only) so both node slots of the one
    // collapsed DOF hold the same value -> symmetric D-scatter. No-op for
    // Dirichlet/pump (no periodic map) and for cross-rank slaves.
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic && s.periodicMap)
    {
        // Diagnostic toggles (MARS_DDT_SYMPROBE workflow): the master->slave
        // broadcast of g overwrites the slave slot, which is NOT self-adjoint;
        // inserting it unpaired between M^-1 and the D-scatter can break the
        // symmetry of A = D M^-1 D^T. These two env vars let us disable each
        // broadcast independently (without recompiling) to MEASURE which one
        // breaks adjointness via the symmetry probe. Default (unset) = both
        // broadcasts ON = current behavior. Read once, cache in a static bool.
        static const bool skipSrGbcast = (std::getenv("MARS_DDT_NO_SR_GBCAST") != nullptr);
        static const bool skipXrGbcast = (std::getenv("MARS_DDT_NO_XR_GBCAST") != nullptr);

        const int* d_partner = s.periodicMap->d_periodicPartner.data();
        if (!skipSrGbcast)
        {
            mars::fem::periodicBroadcastSameRankKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                d_partner, d_nodeOwnership.data(), s.nodeCount, gxAcc.data());
            mars::fem::periodicBroadcastSameRankKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                d_partner, d_nodeOwnership.data(), s.nodeCount, gyAcc.data());
            mars::fem::periodicBroadcastSameRankKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                d_partner, d_nodeOwnership.data(), s.nodeCount, gzAcc.data());
            cudaDeviceSynchronize();
        }
        // Cross-rank leg: refresh a cross-rank slave's OWNED g slot from the
        // master's owner rank. The forward halo above only touches ghosts, and
        // the same-rank kernel skipped this pair (master is a ghost here). The
        // subsequent D-scatter (step c) reads g[slave] and g[master] at the
        // seam; they must be equal or D is asymmetric there. No-op when peers_
        // empty (single-rank / no cross-rank pairs).
        if (s.numRanks > 1 && !skipXrGbcast)
        {
            mars::fem::crossRankPeriodicBroadcast<KeyType, RealType>(*s.periodicMap, gxAcc);
            mars::fem::crossRankPeriodicBroadcast<KeyType, RealType>(*s.periodicMap, gyAcc);
            mars::fem::crossRankPeriodicBroadcast<KeyType, RealType>(*s.periodicMap, gzAcc);
        }
    }

    // Step c: out = D g (un-normalized per-node divergence accumulator).
    zeroVec(outAcc);
    if (eBlocks > 0)
    {
        computeDivergencePerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
            c0, c1, c2, c3, c4, c5, c6, c7,
            gxAcc.data(), gyAcc.data(), gzAcc.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            outAcc.data(), startElem, numLocal);
        cudaDeviceSynchronize();
    }
    s.domain.reverseExchangeNodeHaloAdd(outAcc);
    // Diagnostic gate: MARS_DDT_NO_PERIODIC_SUM=1 disables the cross-rank
    // slave->master Ap sum on the SpMV result. Under the conditional-collapse
    // design, slave_on_A and master_on_D are separate owned DOFs representing
    // two distinct equations. The cross-rank slave->master sum was designed
    // for the OLD shared-DOF design and now merges two equations -> non-SPD.
    // Standard cstone reverseExchangeNodeHaloAdd already routes per-node
    // ghost accumulations to their owners by SFC key, so the cross-rank
    // physics SHOULD be captured by the reverseExchange alone -- as long as
    // the matrix-free operator emits the correct per-node face contributions
    // and the periodic-image element is iterated on both ranks.
    if (std::getenv("MARS_DDT_NO_PERIODIC_SUM") == nullptr)
    {
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, outAcc);
    }
    else if (std::getenv("MARS_DDT_PROBE") != nullptr && s.rank == 0)
    {
        // One-shot print so we see the gate is active.
        static bool printed = false;
        if (!printed) { std::cerr << "[ddt] periodic-sum DISABLED on Ap" << std::endl; printed = true; }
    }

    // NOTE: BD is intentionally NOT applied to the matrix-free DDT operator.
    // Adding S*phi to (D M^-1 D^T) phi breaks the algebraic identity
    // div(u^{n+1}) = 0 that DDT is built for: the corrector
    // u^{n+1} = u** - dt/rho * grad phi cancels exactly only when phi solves
    // (D M^-1 D^T) phi = source. BD on the DDT path produces a phi that does
    // NOT cancel divergence, and the residual is proportional to the BD
    // strength -- destroying the whole point of using DDT. For periodic
    // pressure runs use --pressure-solve=K, which is regularized by BD via
    // the assembled d_valuesPre matrix.
    // if constexpr(hex): the lambda reads n[8]={c0..c7} and n[6] (fixed hex
    // 8-corner stencil). Already dead via if(false); the if constexpr keeps the
    // c4..c7 / n[6] deref out of the tet compile entirely.
    if constexpr (std::is_same_v<ElementTag, HexTag>)
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
        // BD scatter: same gate as the main DDT sum above.
        if (std::getenv("MARS_DDT_NO_PERIODIC_SUM") == nullptr)
        {
            maybePeriodicSum<KeyType, RealType, ElementTag>(s, outAcc);
        }
    }

    // Identity-row enforcement on pinned pressure DOFs (cavity: single corner;
    // channel: outflow-face mask). out[i] = phi[i] makes the row act as identity.
    const int pinDof          = s.pressurePinDof;
    const uint8_t* maskPtr    = (s.d_isPressureBdryDof.size() > 0)
                                ? s.d_isPressureBdryDof.data() : nullptr;
    if (pinDof >= 0 || maskPtr != nullptr)
    {
        const int  numOwnedDofs = s.numOwnedDofs;
        const int* dofPtr       = s.d_node_to_dof.data();
        const uint8_t* ownPtr   = d_nodeOwnership.data();
        const RealType* phPtr   = phi.data();
        RealType* outPtr        = outAcc.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [pinDof, maskPtr, dofPtr, ownPtr, phPtr, outPtr, numOwnedDofs] __device__ (size_t i) {
                if (ownPtr[i] != 1) return;
                int dof = dofPtr[i];
                if (dof < 0 || dof >= numOwnedDofs) return;
                bool flagged = (pinDof >= 0 && dof == pinDof) ||
                               (maskPtr != nullptr && maskPtr[dof] != 0);
                if (flagged) outPtr[i] = phPtr[i];
            });
        cudaDeviceSynchronize();
    }
}

// Plain Euclidean dot over owned nodes (CG operates on un-normalized acc form).
// d_partner (optional): the periodic partner table. A periodic slave (partner>=0)
// ALIASES its master's DOF -- it is the SAME equation living on a second node.
// In the DDT path maybePeriodicSum merges the slave onto the master and zeroes/
// mirrors the slave slot, and bcastMasterToSlave keeps slave==master, so the pair
// is ONE DOF. Counting the slave would weight that DOF twice: by face multiplicity
// for a same-rank pair, and GLOBALLY (slave on this rank + master on the owner
// rank) for a cross-rank pair. The latter breaks SPD -> pAp<=0 -> CG breakdown.
// So we skip ALL slaves; the DOF is counted exactly once on the rank that owns its
// master (where ownership[master]==1). d_partner==nullptr (Dirichlet/pump, no
// collapse) reduces to the old owned-node dot exactly.
template<typename RealType>
RealType ownedDot(const cstone::DeviceVector<RealType>& a,
                  const cstone::DeviceVector<RealType>& b,
                  const cstone::DeviceVector<int>& d_nodeToDof,
                  const uint8_t* d_ownership,
                  size_t numNodes,
                  const int* d_partner = nullptr)
{
    const RealType* aPtr = a.data();
    const RealType* bPtr = b.data();
    const int* dofPtr    = d_nodeToDof.data();
    RealType localSum = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(numNodes),
        [aPtr, bPtr, dofPtr, d_ownership, d_partner] __device__ (size_t i) -> RealType {
            if (d_ownership[i] != 1 || dofPtr[i] < 0) return RealType(0);
            if (d_partner && d_partner[i] >= 0) return RealType(0);
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
int solvePressureDDT(NSStepper<KeyType, RealType, ElementTag>& s,
                     cstone::DeviceVector<RealType>& b_node,
                     cstone::DeviceVector<RealType>& phi_node)
{
    const auto& d_nodeOwnership = s.ownershipMap();
    const int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;

    cstone::DeviceVector<RealType> r(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> p(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Ap(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> z(s.nodeCount, RealType(0));   // M^-1 r
    cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));

    // Periodic DOF aliasing: on a periodic mesh setupNSStepper collapses each
    // same-rank slave node onto its master's DOF (node_to_dof[slave] ==
    // node_to_dof[master]) but slave and master remain TWO distinct NODES, both
    // owned. The CG below is NODE-indexed, so that one collapsed DOF occupies
    // two node slots. To make the CG consistent in DOF space we (i) count each
    // DOF once -- ownedDot skips same-rank slaves via partnerPtr -- and (ii)
    // keep the two slots identical -- bcastMasterToSlave copies master->slave
    // after every update of r/p/phi and inside the operator on g. partnerPtr is
    // null off the periodic path, so both are exact no-ops for Dirichlet/pump.
    const int* partnerPtr = (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
                             && s.periodicMap)
                            ? s.periodicMap->d_periodicPartner.data() : nullptr;
    auto bcastMasterToSlave = [&](cstone::DeviceVector<RealType>& v) {
        if (!partnerPtr) return;
        // (a) same-rank leg: collapsed slave node slot := its owned master's slot.
        mars::fem::periodicBroadcastSameRankKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            partnerPtr, d_nodeOwnership.data(), s.nodeCount, v.data());
        cudaDeviceSynchronize();
        // (b) cross-rank leg: a cross-rank slave is OWNED here but its master is
        // a ghost, so the same-rank kernel skipped it and the cstone halo (owners
        // ->ghosts) never touches an owned slot. Without this refresh its node
        // slot drifts from the master's value on the master-owner rank, and the
        // matrix-free D^T scatter reads asymmetric p across the periodic seam ->
        // seam divergence blows up. Mirrors maybePeriodicSum's cross-rank leg for
        // the SUM direction. No-op when peers_ empty (single-rank / no XR pairs).
        if (s.numRanks > 1)
        {
            mars::fem::crossRankPeriodicBroadcast<KeyType, RealType>(*s.periodicMap, v);
        }
    };

    // Jacobi preconditioning is ON by default. Two earlier session bugs masked
    // its correctness; both are now fixed:
    //   1. The per-element scatter loop-bound double-count (FIXED in this
    //      session): used to make the matrix-free operator disagree with the
    //      assembled s.d_valuesDDT at rank-boundary nodes, so the cached diag
    //      didn't match what CG saw.
    //   2. The cavity pin-row diag mismatch (FIXED in the setup-time cache
    //      build): matrix-free applyDDTPerNode enforces identity (diag=1) at
    //      the pin, but enforcePinRowKeepDiagKernel keeps the assembled value
    //      (~0.04) in s.d_valuesDDT. The cache build now overrides the pin
    //      slot to 1.0 so Jacobi z[pin] = r[pin]/1.0 matches the operator.
    // With both fixes in, Jacobi PCG should converge in ~30-60 iters on cube16
    // cavity (vs un-precond 104) and ~400-700 iters on the wing 15M-hex (vs
    // un-precond stalling at 1e-3 after 20k). Env MARS_DDT_NO_JACOBI=1 opts
    // OUT for diagnostic A/B against the un-precond baseline.
    bool useJacobi = (s.d_diagDDT.size() == static_cast<size_t>(s.numOwnedDofs)
                      && s.numOwnedDofs > 0);
    {
        const char* ev = std::getenv("MARS_DDT_NO_JACOBI");
        if (ev && std::string(ev) != "0") useJacobi = false;
    }
    if (s.rank == 0)
    {
        static bool printedOnce = false;
        if (!printedOnce)
        {
            std::cout << "[cg-ddt] preconditioner = "
                      << (useJacobi
                          ? (std::is_same_v<ElementTag, TetTag>
                             ? "Jacobi (diag from s.d_diagDDT; matrix-free area+mass)"
                             : "Jacobi (diag from s.d_diagDDT; extracted from assembled s.d_valuesDDT)")
                          : "NONE (un-precond)")
                      << "\n";
            if (useJacobi)
            {
                // One-shot diag range print so we see what we're dividing by.
                auto dp = thrust::device_pointer_cast(s.d_diagDDT.data());
                auto mm = thrust::minmax_element(thrust::device, dp,
                                                  dp + s.numOwnedDofs);
                RealType dmin = 0, dmax = 0;
                cudaMemcpy(&dmin, thrust::raw_pointer_cast(&*mm.first),  sizeof(RealType), cudaMemcpyDeviceToHost);
                cudaMemcpy(&dmax, thrust::raw_pointer_cast(&*mm.second), sizeof(RealType), cudaMemcpyDeviceToHost);
                std::cout << "[cg-ddt] cached diag range on rank 0: ["
                          << dmin << ", " << dmax << "]\n";
                // A constant diagonal carries no per-DOF scaling: Jacobi then
                // reduces to z = r / const, i.e. CG runs effectively
                // UN-preconditioned even though useJacobi==true. For tet this
                // means the matrix-free area+mass accumulator came back empty
                // (all zero -> floored to 1.0): check that d_areaVec and
                // d_massNode were populated before this build. Loud on purpose
                // so the "preconditioner present but useless" state is never
                // silent (it only shows up as a sluggish, high-iter solve).
                if (dmin == dmax)
                    std::cout << "[cg-ddt] WARNING: DDT Jacobi diagonal is CONSTANT ("
                              << dmin << ") -> no per-DOF preconditioning; CG will be slow. "
                              << (std::is_same_v<ElementTag, TetTag>
                                  ? "Tet matrix-free diag accumulator was empty/uniform."
                                  : "Assembled diag was empty/uniform.")
                              << "\n";
            }
            printedOnce = true;
        }
    }

    // ---- MARS_DDT_SYMPROBE: cross-rank operator-symmetry measurement ----
    // A = D M^-1 D^T must be symmetric for CG (pAp>0). On >1 rank with a
    // periodic seam the CG breaks (pAp<=0) even with Jacobi off, so A itself
    // is non-SPD across the seam. The suspect is the master->slave g-broadcast
    // (an overwrite is not self-adjoint) inserted between M^-1 and the
    // D-scatter inside applyDDTPerNode. This probe MEASURES adjointness:
    // for a cross-rank seam pair of DOFs (i,j) it computes A[i,j]=(A e_j)[i]
    // and A[j,i]=(A e_i)[j] and prints |A[i,j]-A[j,i]|. Symmetric => roundoff;
    // asymmetric => O(operator). Run it under the two NO_*_GBCAST toggles to
    // isolate which broadcast breaks symmetry. Runs ONCE (static guard), only
    // on the first solve, and is fully gated -- default behavior unchanged.
    if constexpr (std::is_same_v<ElementTag, HexTag>)
    {
        static bool symProbeDoneV2 = false;
        if (!symProbeDoneV2 && std::getenv("MARS_DDT_SYMPROBE") != nullptr
            && s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
            && s.periodicMap != nullptr)
        {
            symProbeDoneV2 = true;

            // Bilinear-form symmetry test (no pair/coupling guesswork): for the
            // self-adjoint operator A=D M^-1 D^T we must have <u,A v> == <v,A u>
            // for ANY u,v. Build two arbitrary OWNED vectors, apply A to each,
            // and compare the two cross inner products. The gap IS the asymmetry.
            // No RNG (unavailable on device-build); use deterministic index hashes.
            // Run under the NO_*_GBCAST toggles to see which broadcast breaks it.
            cstone::DeviceVector<RealType> uVec(s.nodeCount, RealType(0));
            cstone::DeviceVector<RealType> vVec(s.nodeCount, RealType(0));
            cstone::DeviceVector<RealType> AuVec(s.nodeCount, RealType(0));
            cstone::DeviceVector<RealType> AvVec(s.nodeCount, RealType(0));
            cstone::DeviceVector<RealType> g1(s.nodeCount, RealType(0));
            cstone::DeviceVector<RealType> g2(s.nodeCount, RealType(0));
            cstone::DeviceVector<RealType> g3(s.nodeCount, RealType(0));
            {
                RealType* uP = uVec.data();
                RealType* vP = vVec.data();
                const int*     dofP = s.d_node_to_dof.data();
                const uint8_t* ownP = d_nodeOwnership.data();
                const int*     parP = partnerPtr;
                // owned, non-slave nodes get a smooth deterministic value; the
                // exact values don't matter, only that u != v and both nonzero.
                thrust::for_each(thrust::device,
                    thrust::counting_iterator<size_t>(0),
                    thrust::counting_iterator<size_t>(s.nodeCount),
                    [uP, vP, dofP, ownP, parP] __device__ (size_t i) {
                        if (ownP[i] != 1 || dofP[i] < 0) return;
                        if (parP && parP[i] >= 0) return;
                        int d = dofP[i];
                        uP[i] = RealType(1) + RealType(d % 7);
                        vP[i] = RealType(1) + RealType((d * 3 + 1) % 11);
                    });
                cudaDeviceSynchronize();
            }
            // Av and Au through the SAME operator path the CG uses (incl. the
            // master->slave broadcasts + halo the toggles control).
            bcastMasterToSlave(vVec); s.domain.exchangeNodeHalo(vVec);
            applyDDTPerNode<KeyType, RealType, ElementTag>(s, vVec, AvVec, g1, g2, g3);
            bcastMasterToSlave(uVec); s.domain.exchangeNodeHalo(uVec);
            applyDDTPerNode<KeyType, RealType, ElementTag>(s, uVec, AuVec, g1, g2, g3);

            RealType uAv = ownedDot<RealType>(uVec, AvVec, s.d_node_to_dof,
                                              d_nodeOwnership.data(), s.nodeCount, partnerPtr);
            RealType vAu = ownedDot<RealType>(vVec, AuVec, s.d_node_to_dof,
                                              d_nodeOwnership.data(), s.nodeCount, partnerPtr);
            if (s.rank == 0)
            {
                RealType denom = std::max(std::abs(uAv), std::abs(vAu));
                RealType rel = (denom > RealType(0)) ? std::abs(uAv - vAu) / denom : RealType(0);
                std::cout << "  [SYMPROBE] <u,Av>=" << uAv << "  <v,Au>=" << vAu
                          << "  |diff|=" << std::abs(uAv - vAu)
                          << "  rel=" << rel
                          << "   (symmetric => rel~roundoff; asymmetric => rel~O(1))\n";
            }
        }
        // ---- legacy pair-probe (kept disabled; superseded by the bilinear test) ----
        static bool symProbeDone = false;
        if (false && !symProbeDone && std::getenv("MARS_DDT_SYMPROBE_PAIRS") != nullptr
            && s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
            && s.periodicMap != nullptr)
        {
            symProbeDone = true;

            // Step 1 (host, tiny): each rank picks up to 3 of its OWNED
            // cross-rank periodic masters i (those with a remote slave) and,
            // for each, an OWNED interior neighbor j that shares an element
            // with i (via the node->element CSR + element->node connectivity).
            // We gather (masterNode i, neighborNode j) pairs to host so the
            // collective apply below can build e_i / e_j on device. Only the
            // rank that owns the master proposes -- the others contribute 0 to
            // the MPI reduction, so no DOF is counted twice and none is missed.
            const auto& xr = s.periodicMap->cross_;
            std::vector<int> hMasters;
            const int wantPairs = 3;
            if (!xr.d_recvOwnedMasterIds_.empty())
            {
                int nMasters = std::min<int>(int(xr.d_recvOwnedMasterIds_.size()), wantPairs);
                hMasters.resize(nMasters);
                cudaMemcpy(hMasters.data(), xr.d_recvOwnedMasterIds_.data(),
                           nMasters * sizeof(int), cudaMemcpyDeviceToHost);
            }

            // Pull the small CSR + connectivity + dof/partner/ownership arrays
            // to host once to find neighbors. nodeCount-sized D2H, but the
            // probe is one-shot and diagnostic-only (no steady-state cost).
            const auto& d_n2eOff = s.domain.getNodeToElementOffsets();
            const auto& d_n2eLst = s.domain.getNodeToElementList();
            const auto& d_conn   = s.domain.getElementToNodeConnectivity();
            auto cp = connPtrs<ElementTag, KeyType>(d_conn);
            constexpr int NPE = ElemTraits<ElementTag>::NodesPerElem;

            std::vector<KeyType> hOff(s.nodeCount + 1);
            cudaMemcpy(hOff.data(), d_n2eOff.data(),
                       (s.nodeCount + 1) * sizeof(KeyType), cudaMemcpyDeviceToHost);
            std::vector<KeyType> hLst(hOff.empty() ? 0 : size_t(hOff.back()));
            if (!hLst.empty())
                cudaMemcpy(hLst.data(), d_n2eLst.data(),
                           hLst.size() * sizeof(KeyType), cudaMemcpyDeviceToHost);
            std::vector<int>     hN2D(s.nodeCount);
            std::vector<uint8_t> hOwn(s.nodeCount);
            std::vector<int>     hPartner(s.nodeCount, -1);
            cudaMemcpy(hN2D.data(), s.d_node_to_dof.data(),
                       s.nodeCount * sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(hOwn.data(), d_nodeOwnership.data(),
                       s.nodeCount * sizeof(uint8_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(hPartner.data(), s.periodicMap->d_periodicPartner.data(),
                       s.nodeCount * sizeof(int), cudaMemcpyDeviceToHost);
            // Element connectivity columns (NPE node-id arrays). Hex only here.
            std::vector<std::vector<KeyType>> hCol(NPE);
            for (int c = 0; c < NPE; ++c)
            {
                hCol[c].resize(s.elementCount);
                cudaMemcpy(hCol[c].data(), cp[c],
                           s.elementCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
            }

            // For each master i, scan its incident elements (node->element CSR)
            // and pick a neighbor node j off the SAME element via the
            // element->node connectivity. Sharing an element makes A[i,j] a real
            // nonzero entry of the 27-pt D M^-1 D^T stencil (A[i,j]!=0 iff i,j
            // share an element), so |A[i,j]-A[j,i]| actually measures asymmetry
            // instead of probing a structural zero. j is required owned with a
            // distinct valid owned DOF and partner<0 (a clean interior DOF, not a
            // periodic alias). Owned-on-this-rank => i,j are same-rank here; the
            // print labels that, and a future remotely-owned j would auto-label
            // cross-rank via the ownerOfJ reduction below.
            std::vector<std::pair<int,int>> pairsIJ;  // (masterNode i, neighbor j)
            for (int mi : hMasters)
            {
                if (mi < 0 || mi >= int(s.nodeCount)) continue;
                int dofI = hN2D[mi];
                if (dofI < 0 || dofI >= s.numOwnedDofs) continue;
                int found = -1;
                KeyType b = hOff[mi], e = hOff[mi + 1];
                for (KeyType t = b; t < e && found < 0; ++t)
                {
                    size_t elem = size_t(hLst[t]);
                    if (elem >= s.elementCount) continue;
                    for (int c = 0; c < NPE; ++c)
                    {
                        int nj = int(hCol[c][elem]);
                        if (nj < 0 || nj >= int(s.nodeCount)) continue;
                        if (nj == mi) continue;
                        if (hOwn[nj] != 1) continue;
                        if (hPartner[nj] >= 0) continue;            // skip slaves/aliases
                        int dofJ = hN2D[nj];
                        if (dofJ < 0 || dofJ >= s.numOwnedDofs) continue;
                        if (dofJ == dofI) continue;                 // distinct DOF
                        found = nj;
                        break;
                    }
                }
                if (found >= 0) pairsIJ.emplace_back(mi, found);
            }

            // Step 2: round-robin a GLOBAL list of pairs across ranks so the
            // collective applyDDTPerNode is called the same number of times on
            // every rank. Each rank announces how many pairs it has; we then
            // process up to wantPairs pairs total, one per (rank, slot).
            int myPairs = int(pairsIJ.size());
            std::vector<int> allCounts(s.numRanks, 0);
            MPI_Allgather(&myPairs, 1, MPI_INT, allCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

            // Device scratch for the unit vectors and the operator output. Sized
            // like the existing scaffold (phi over nodeCount, Ap over nodeCount).
            cstone::DeviceVector<RealType> phiE(s.nodeCount, RealType(0));
            cstone::DeviceVector<RealType> ApE(s.nodeCount, RealType(0));
            const int*     n2dPtr = s.d_node_to_dof.data();
            const uint8_t* ownPtr = d_nodeOwnership.data();

            // Set phiE = e_{dof} on the owner rank only (others leave it 0); the
            // halo exchange then fills any ghost copies of that node so the
            // collective D^T scatter sees a globally-consistent unit vector.
            auto setUnit = [&](int targetDof, bool ownerHere) {
                cudaMemset(phiE.data(), 0, s.nodeCount * sizeof(RealType));
                if (ownerHere)
                {
                    thrust::for_each(thrust::device,
                        thrust::counting_iterator<size_t>(0),
                        thrust::counting_iterator<size_t>(s.nodeCount),
                        [n2dPtr, ownPtr, targetDof, p = phiE.data()] __device__ (size_t k) {
                            if (ownPtr[k] == 1 && n2dPtr[k] == targetDof) p[k] = RealType(1);
                        });
                    cudaDeviceSynchronize();
                }
                s.domain.exchangeNodeHalo(phiE);
            };
            // Read (A e)[node with owned dof==targetDof] on the owner rank only;
            // returns 0 on non-owners so MPI_SUM picks exactly one contributor.
            auto readAt = [&](int targetDof, bool ownerHere) -> RealType {
                if (!ownerHere) return RealType(0);
                RealType acc = thrust::transform_reduce(thrust::device,
                    thrust::counting_iterator<size_t>(0),
                    thrust::counting_iterator<size_t>(s.nodeCount),
                    [n2dPtr, ownPtr, targetDof, a = ApE.data()] __device__ (size_t k) -> RealType {
                        return (ownPtr[k] == 1 && n2dPtr[k] == targetDof) ? a[k] : RealType(0);
                    }, RealType(0), thrust::plus<RealType>());
                return acc;
            };

            auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
            int printed = 0;
            for (int r = 0; r < s.numRanks && printed < wantPairs; ++r)
            {
                int slots = std::min(allCounts[r], wantPairs - printed);
                for (int sIdx = 0; sIdx < slots; ++sIdx, ++printed)
                {
                    bool ownerHere = (s.rank == r);
                    // Local DOF indices of this pair live ONLY on the owner rank
                    // r; other ranks pass -1 (their set/read self-select to 0).
                    int dofI = -1, dofJ = -1;
                    if (ownerHere)
                    {
                        dofI = hN2D[pairsIJ[sIdx].first];
                        dofJ = hN2D[pairsIJ[sIdx].second];
                    }

                    // A[j,i] = (A e_i)[j]: collective apply with phi=e_i.
                    setUnit(dofI, ownerHere);
                    applyDDTPerNode<KeyType, RealType, ElementTag>(s, phiE, ApE, gx, gy, gz);
                    RealType Aji_local = readAt(dofJ, ownerHere);

                    // A[i,j] = (A e_j)[i]: collective apply with phi=e_j.
                    setUnit(dofJ, ownerHere);
                    applyDDTPerNode<KeyType, RealType, ElementTag>(s, phiE, ApE, gx, gy, gz);
                    RealType Aij_local = readAt(dofI, ownerHere);

                    // A[i,i] = (A e_i)[i]: diagonal sanity. If this comes back 0
                    // the unit-vector/readAt machinery itself is broken (a node
                    // is always coupled to itself in the 27-pt DDT stencil), so a
                    // zero here, not just |Aij-Aji|=0, is the first thing to check:
                    // it means setUnit/exchange/applyDDTPerNode/readAt are wrong
                    // and the off-diagonal numbers are meaningless. A[i,i] must be
                    // > 0 (SPD diagonal). Same collective shape as the legs above.
                    setUnit(dofI, ownerHere);
                    applyDDTPerNode<KeyType, RealType, ElementTag>(s, phiE, ApE, gx, gy, gz);
                    RealType Aii_local = readAt(dofI, ownerHere);

                    // Each entry is owned by exactly one rank (owner of i for
                    // A[i,j]/A[i,i], owner of j for A[j,i]); everyone else packed
                    // 0, so SUM merges them with no double-count and no miss.
                    RealType Aji = 0, Aij = 0, Aii = 0;
                    MPI_Allreduce(&Aji_local, &Aji, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(&Aij_local, &Aij, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(&Aii_local, &Aii, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);

                    // Rank that owns j (the proposing rank r owns both i and j
                    // because the neighbor walk only accepts owned nodes; -1 from
                    // everyone else, MAX picks the real owner). Report same-rank
                    // vs cross-rank so a future selection that picks a remotely
                    // owned j is labelled correctly without further changes.
                    int ownerOfJ_local = ownerHere ? s.rank : -1;
                    int ownerOfJ = -1;
                    MPI_Allreduce(&ownerOfJ_local, &ownerOfJ, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                    if (s.rank == 0)
                    {
                        // dofI/dofJ are owner-local; report owner rank + locals
                        // for traceability. The signal is |Aij-Aji|; A[i,i]>0
                        // gates whether that signal is trustworthy at all.
                        int repI = -1, repJ = -1;
                        if (r == 0) { repI = dofI; repJ = dofJ; }
                        const char* loc = (ownerOfJ == r) ? "same-rank" : "cross-rank";
                        std::cout << "  [MARS_DDT_SYMPROBE] seam pair (ownerRank(i)=" << r
                                  << ", ownerRank(j)=" << ownerOfJ << ", " << loc
                                  << ", dofI=" << repI << ", dofJ=" << repJ << "): "
                                  << std::scientific << std::setprecision(6)
                                  << "A[i,i]=" << Aii
                                  << " A[i,j]=" << Aij << " A[j,i]=" << Aji
                                  << " |diff|=" << std::abs(Aij - Aji)
                                  << std::defaultfloat << "\n";
                    }
                }
            }
            if (printed == 0 && s.rank == 0)
                std::cout << "  [MARS_DDT_SYMPROBE] no cross-rank seam masters found "
                             "(single-rank or no cross-rank periodic pairs)\n";
        }
    }

    // Fresh start: phi = 0. With b[i]=0 at every pinned DOF and out[i]=phi[i]
    // in the SpMV, the pinned-row equation reads 0 = phi[i] => phi[i] stays 0.
    thrust::fill(thrust::device_pointer_cast(phi_node.data()),
                 thrust::device_pointer_cast(phi_node.data() + s.nodeCount), RealType(0));

    const int pinDof          = s.pressurePinDof;
    const uint8_t* maskPtr    = (s.d_isPressureBdryDof.size() > 0)
                                ? s.d_isPressureBdryDof.data() : nullptr;
    if (pinDof >= 0 || maskPtr != nullptr)
    {
        const int  numOwnedDofs = s.numOwnedDofs;
        const int* dofPtr       = s.d_node_to_dof.data();
        const uint8_t* ownPtr   = d_nodeOwnership.data();
        RealType* bPtr          = b_node.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [pinDof, maskPtr, dofPtr, ownPtr, bPtr, numOwnedDofs] __device__ (size_t i) {
                if (ownPtr[i] != 1) return;
                int dof = dofPtr[i];
                if (dof < 0 || dof >= numOwnedDofs) return;
                bool flagged = (pinDof >= 0 && dof == pinDof) ||
                               (maskPtr != nullptr && maskPtr[dof] != 0);
                if (flagged) bPtr[i] = RealType(0);
            });
        cudaDeviceSynchronize();
    }

    // r = b - A*phi  (phi=0 so r = b on owned slots). z = M^-1 r. p = z.
    thrust::copy(thrust::device_pointer_cast(b_node.data()),
                 thrust::device_pointer_cast(b_node.data() + s.nodeCount),
                 thrust::device_pointer_cast(r.data()));
    // NOTE: the periodic master<->slave identity is now applied INSIDE the
    // operator (P on phi at the top of applyDDTPerNode, P^T = maybePeriodicSum on
    // outAcc) -- the P^T A P pattern. The CG iterate r/p/phi is therefore NEVER
    // re-broadcast between matvecs (that was the anti-pattern that broke Krylov
    // conjugacy on >1 rank). The slave slot of r is left at b[slave]=0 by the RHS
    // assembly; it rides along inert because ownedDot skips slaves and the next
    // matvec's P re-establishes phi[slave]=phi[master].
    if (useJacobi)
    {
        jacobiPrecondNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            r.data(), s.d_diagDDT.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.diagDDTEpsClip, z.data(), s.nodeCount);
        cudaDeviceSynchronize();
    }
    else
    {
        thrust::copy(thrust::device_pointer_cast(r.data()),
                     thrust::device_pointer_cast(r.data() + s.nodeCount),
                     thrust::device_pointer_cast(z.data()));
    }
    thrust::copy(thrust::device_pointer_cast(z.data()),
                 thrust::device_pointer_cast(z.data() + s.nodeCount),
                 thrust::device_pointer_cast(p.data()));

    // rho = (r, z): preconditioned inner product drives alpha/beta.
    // Convergence check still uses the Euclidean residual |r| so the tolerance
    // semantics match the un-preconditioned baseline and the K-path numbers.
    RealType rho_old = ownedDot<RealType>(r, z, s.d_node_to_dof,
                                          d_nodeOwnership.data(), s.nodeCount, partnerPtr);
    RealType rr0     = ownedDot<RealType>(r, r, s.d_node_to_dof,
                                          d_nodeOwnership.data(), s.nodeCount, partnerPtr);
    RealType r0_norm = std::sqrt(rr0);
    // Diagnostic capture (every step, regardless of the gated per-iter print):
    // |r0| is the pressure RHS magnitude; lastPressResid starts at the not-yet-
    // iterated ratio 1 and is overwritten with the latest |r|/|r0| in the loop.
    s.lastPressR0    = r0_norm;
    s.lastPressResid = RealType(1);
    if (r0_norm < std::numeric_limits<RealType>::min()) { s.lastPressResid = RealType(0); return 0; }
    const RealType absTol = s.tolerance * r0_norm;

    int iters = -2;
    // Live progress every N CG iterations on rank 0 so large-mesh runs show
    // residual-decay rate instead of going silent for minutes. Each line
    // prints |r|/|r0| so the user can extrapolate iters-to-convergence.
    // liveEvery: print |r|/|r0| every N iters on rank 0. Default 100 for
    // maxIter>=500 (cube/wing scale), 25 for smaller runs. Env override
    // MARS_DDT_CG_PRINT_EVERY=N forces N (useful on wing to see whether the
    // solve is making progress without waiting for 100 SpMVs to elapse, which
    // at 4M owned DOFs on H100 can be tens of seconds per print interval).
    int liveEvery = (s.maxIter >= 500) ? 100 : 25;
    {
        const char* ev = std::getenv("MARS_DDT_CG_PRINT_EVERY");
        if (ev)
        {
            int v = std::atoi(ev);
            if (v > 0) liveEvery = v;
        }
    }
    for (int it = 0; it < s.maxIter; ++it)
    {
        // applyDDTPerNode applies P (periodic master->slave) and the ghost halo
        // to its input internally, so we pass p straight in -- no iterate mutation.
        applyDDTPerNode<KeyType, RealType, ElementTag>(s, p, Ap, gx, gy, gz);

        RealType pAp = ownedDot<RealType>(p, Ap, s.d_node_to_dof,
                                          d_nodeOwnership.data(), s.nodeCount, partnerPtr);
        // MARS_DDT_PAP_DEBUG: on the FIRST solve's FIRST iteration, dump the
        // actual scalars that decide pAp<=0. The PER-RANK local contribution
        // (before the Allreduce) is the key unknown: if one rank's local pAp is
        // a large negative or NaN/inf, the operator/vectors are inconsistent on
        // that rank's seam. Also dump |p|,|Ap|, min/max(Ap) over owned to see if
        // Ap carries a poison value. Gated + once -> zero cost in production.
        {
            static bool papDbgDone = false;
            if (!papDbgDone && std::getenv("MARS_DDT_PAP_DEBUG") != nullptr && it == 0)
            {
                papDbgDone = true;
                const RealType* pPtr  = p.data();
                const RealType* ApPtr = Ap.data();
                const int*      dofP  = s.d_node_to_dof.data();
                const uint8_t*  ownP  = d_nodeOwnership.data();
                const int*      parP  = partnerPtr;
                // local pAp (same gate as ownedDot, but keep the un-reduced value)
                RealType localPap = thrust::transform_reduce(thrust::device,
                    thrust::counting_iterator<size_t>(0),
                    thrust::counting_iterator<size_t>(s.nodeCount),
                    [pPtr, ApPtr, dofP, ownP, parP] __device__ (size_t i) -> RealType {
                        if (ownP[i] != 1 || dofP[i] < 0) return RealType(0);
                        if (parP && parP[i] >= 0) return RealType(0);
                        return pPtr[i] * ApPtr[i];
                    }, RealType(0), thrust::plus<RealType>());
                // local |p|^2 and |Ap|^2 over owned (same gate)
                RealType localPP = thrust::transform_reduce(thrust::device,
                    thrust::counting_iterator<size_t>(0),
                    thrust::counting_iterator<size_t>(s.nodeCount),
                    [pPtr, dofP, ownP, parP] __device__ (size_t i) -> RealType {
                        if (ownP[i] != 1 || dofP[i] < 0) return RealType(0);
                        if (parP && parP[i] >= 0) return RealType(0);
                        return pPtr[i] * pPtr[i];
                    }, RealType(0), thrust::plus<RealType>());
                RealType localApAp = thrust::transform_reduce(thrust::device,
                    thrust::counting_iterator<size_t>(0),
                    thrust::counting_iterator<size_t>(s.nodeCount),
                    [ApPtr, dofP, ownP, parP] __device__ (size_t i) -> RealType {
                        if (ownP[i] != 1 || dofP[i] < 0) return RealType(0);
                        if (parP && parP[i] >= 0) return RealType(0);
                        return ApPtr[i] * ApPtr[i];
                    }, RealType(0), thrust::plus<RealType>());
                cudaDeviceSynchronize();
                // gather per-rank local pAp to rank 0
                std::vector<double> allLocal(s.numRanks, 0.0);
                double myLocal = static_cast<double>(localPap);
                MPI_Gather(&myLocal, 1, MPI_DOUBLE, allLocal.data(), 1, MPI_DOUBLE,
                           0, MPI_COMM_WORLD);
                double gPP = 0, gApAp = 0, lPP = localPP, lApAp = localApAp;
                MPI_Allreduce(&lPP,   &gPP,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&lApAp, &gApAp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                if (s.rank == 0)
                {
                    std::cout << "  [pAp-dbg] pAp(global)=" << pAp
                              << "  rho_old(r,z)=" << rho_old
                              << "  rr0(r,r)=" << rr0
                              << "  |p|=" << std::sqrt(gPP)
                              << "  |Ap|=" << std::sqrt(gApAp) << "\n";
                    std::cout << "  [pAp-dbg] per-rank local pAp:";
                    for (int rr = 0; rr < s.numRanks; ++rr)
                        std::cout << " r" << rr << "=" << allLocal[rr];
                    std::cout << "\n";
                }
            }
        }
        if (pAp <= RealType(0)) { iters = -2; break; }
        RealType alpha = rho_old / pAp;

        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            phi_node.data(), phi_node.data(), p.data(), alpha,
            s.d_node_to_dof.data(), d_nodeOwnership.data(), s.nodeCount);
        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            r.data(), r.data(), Ap.data(), -alpha,
            s.d_node_to_dof.data(), d_nodeOwnership.data(), s.nodeCount);
        cudaDeviceSynchronize();
        // No iterate re-broadcast: r[slave] (and phi[slave]) ride along inert.
        // The slave is excluded from every ownedDot, and each matvec re-applies
        // P (phi[slave]:=phi[master]) inside applyDDTPerNode, so the stale slave
        // slots never enter a reduction or the operator. Re-broadcasting them
        // here is the documented anti-pattern (a projection between matvecs that
        // breaks CG conjugacy -> the >1-rank converge-then-diverge stall).

        // Convergence on the Euclidean residual; alpha/beta on the
        // preconditioned (r,z) inner product. Equivalent to un-precond CG when
        // z == r (the useJacobi=false fallback path).
        RealType rr_new  = ownedDot<RealType>(r, r, s.d_node_to_dof,
                                              d_nodeOwnership.data(), s.nodeCount, partnerPtr);
        // Keep the diagnostic resid current on every iter (the print below is
        // gated; this store is not, so lastPressResid is the true exit ratio).
        s.lastPressResid = std::sqrt(rr_new) / std::max(r0_norm, std::numeric_limits<RealType>::min());
        if (s.rank == 0 && (it == 0 || (it + 1) % liveEvery == 0))
        {
            std::cout << "    [cg-ddt] iter " << std::setw(6) << (it + 1)
                      << "  |r|/|r0| = " << std::scientific << std::setprecision(3)
                      << s.lastPressResid
                      << std::defaultfloat << "\n" << std::flush;
        }
        if (std::sqrt(rr_new) < absTol) { iters = it + 1; break; }
        // z = M^-1 r (or z := r if no preconditioner).
        if (useJacobi)
        {
            jacobiPrecondNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                r.data(), s.d_diagDDT.data(),
                s.d_node_to_dof.data(), d_nodeOwnership.data(),
                s.diagDDTEpsClip, z.data(), s.nodeCount);
            cudaDeviceSynchronize();
        }
        else
        {
            thrust::copy(thrust::device_pointer_cast(r.data()),
                         thrust::device_pointer_cast(r.data() + s.nodeCount),
                         thrust::device_pointer_cast(z.data()));
        }
        RealType rho_new = ownedDot<RealType>(r, z, s.d_node_to_dof,
                                              d_nodeOwnership.data(), s.nodeCount, partnerPtr);
        RealType beta    = rho_new / rho_old;
        // p = z + beta * p   (un-precond reduces to p = r + beta*p when z=r).
        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            p.data(), z.data(), p.data(), beta,
            s.d_node_to_dof.data(), d_nodeOwnership.data(), s.nodeCount);
        cudaDeviceSynchronize();
        // p is NOT re-broadcast: its slave slot rides along inert (excluded from
        // ownedDot; re-established by P inside the next matvec). Mutating the
        // search direction here breaks M-conjugacy -- the P^T A P structure keeps
        // the iterate consistent without touching it (MFEM/deal.II pattern).
        rho_old = rho_new;
    }

    // cg_p=-2 means the loop ran all maxIter without hitting absTol (iters stays
    // at its -2 init) -- a NON-convergence, NOT a pAp<=0 breakdown (that path
    // sets -2 too but we instrument pAp separately). Print the exit state so we
    // can tell "slow but decreasing" from "stalled" from "residual went NaN".
    if (iters == -2 && std::getenv("MARS_DDT_PAP_DEBUG") != nullptr && s.rank == 0)
    {
        std::cout << "  [cg-exit] ran maxIter=" << s.maxIter
                  << " without converging: final |r|/|r0|=" << s.lastPressResid
                  << "  absTol/|r0|=" << s.tolerance << "\n";
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
// Periodic post-scatter helper. Single-rank: gated intra-rank kernel is a
// no-op because the DOF collapse aliased nodeToDof[slave]=nodeToDof[master],
// so per-DOF accumulators already merged. Multi-rank periodic: the gated
// intra-rank kernel handles same-rank pairs (where the local-master-collapse
// hit); the cross-rank pair sum bridges pairs whose master is owned on a
// remote rank.
//
// Pre-condition: d_field is per-NODE (size = nodeCount). Caller has scattered
// owned-element contributions onto all 8 corners and (optionally) run
// reverseExchangeNodeHaloAdd to fold standard cross-rank ghost contributions
// back to owners.
// =============================================================================
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
inline void maybePeriodicSum(NSStepper<KeyType, RealType, ElementTag>& s,
                             cstone::DeviceVector<RealType>& d_field)
{
    if (s.bcKind != NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic) return;
    if (s.periodicMap == nullptr) return;

    const auto& d_own = s.domain.getNodeOwnershipMap();
    size_t nNodes = d_field.size();
    int blk = 256, grd = int((nNodes + blk - 1) / blk);

    // (a) Gated intra-rank pair sum. Same-rank pairs: master += slave (atomic),
    // slave=0. Cross-rank pairs (own[master]!=1): no-op.
    mars::fem::periodicPairSumKernel<RealType><<<grd, blk>>>(
        s.periodicMap->d_periodicPartner.data(),
        d_own.data(), nNodes, d_field.data());
    cudaDeviceSynchronize();

    // (b) Cross-rank pair sum. No-op when peers_ empty.
    if (s.numRanks > 1)
    {
        mars::fem::crossRankPeriodicPairSum<KeyType, RealType>(
            *s.periodicMap, d_field);
    }
}

// =============================================================================
// Diagnostic helpers used by the sub-steps below. RMS of a 3-vector magnitude
// and max of |scalar| over owned interior DOFs (boundary skipped because lid
// Dirichlet picks up an artificial step). Defined above the step routines so
// they are visible at template instantiation time.
// =============================================================================
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
RealType rmsOwnedInterior3(NSStepper<KeyType, RealType, ElementTag>& s,
                            const cstone::DeviceVector<RealType>& d_gx,
                            const cstone::DeviceVector<RealType>& d_gy,
                            const cstone::DeviceVector<RealType>& d_gz)
{
    const auto& d_nodeOwnership = s.ownershipMap();
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
RealType rmsOwnedInterior1(NSStepper<KeyType, RealType, ElementTag>& s,
                            const cstone::DeviceVector<RealType>& d_q)
{
    const auto& d_nodeOwnership = s.ownershipMap();
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

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
RealType maxOwnedInteriorAbs(NSStepper<KeyType, RealType, ElementTag>& s,
                             const cstone::DeviceVector<RealType>& d_q)
{
    const auto& d_nodeOwnership = s.ownershipMap();
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

// Set to >0 to dump per-phase diagnostics for the first N calls. Declared here
// (before the predictor/corrector) so all step functions can read it.
static int g_nsDebugStepsLeft = 0;

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void runPredictorStep(NSStepper<KeyType, RealType, ElementTag>& s, RealType dt, RealType rho)
{
    const auto& d_nodeOwnership = s.ownershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    // connPtrs gives NodesPerElem pointers; c4..c7 stay nullptr for tet (the
    // templated Phase-1 kernels read c0..c3 only when ElementTag==TetTag).
    auto cp = connPtrs<ElementTag, KeyType>(d_conn);
    const KeyType* c0 = cp[0]; const KeyType* c1 = cp[1];
    const KeyType* c2 = cp[2]; const KeyType* c3 = cp[3];
    const KeyType* c4 = nullptr; const KeyType* c5 = nullptr;
    const KeyType* c6 = nullptr; const KeyType* c7 = nullptr;
    if constexpr (std::is_same_v<ElementTag, HexTag>) { c4 = cp[4]; c5 = cp[5]; c6 = cp[6]; c7 = cp[7]; }
    // OWNED-only per-element scatter (proven-green pattern, commit 5da5d6a).
    // Looping owned+halo here would DOUBLE-COUNT shared rank-boundary faces
    // because the halo element on this rank is the same physical element the
    // owning rank already scatters; reverseExchangeNodeHaloAdd then folds the
    // ghost-side contribution back to owner, giving 2x at shared faces and 4x
    // at shared corners.
    const size_t startElem = s.domain.startIndex();
    const size_t numLocal  = s.domain.localElementCount();
    const int nodeBlocks   = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    const int eBlocks      = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
    const RealType invRho  = RealType(1) / rho;

    using AdvScheme = typename NSStepper<KeyType, RealType, ElementTag>::AdvScheme;
    const bool   bjMode  = (s.advScheme == AdvScheme::BarthJespersen);
    // Runtime mode passed to the flux kernel: 0=skew, 1=upwind, 2=Barth-Jespersen.
    const int    advMode = (s.advScheme == AdvScheme::Skew)   ? 0
                         : (s.advScheme == AdvScheme::Upwind) ? 1
                                                              : 2;

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
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic && s.periodicMap)
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
                computeGradientPerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    s.d_p.data(),
                    s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                    d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
                    startElem, numLocal);
            }
            else
            {
                applyDivTransposePerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
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
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_gxAcc);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_gyAcc);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_gzAcc);
        normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
            s.d_massNode.data(),
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
        s.lastGradPRms = rmsOwnedInterior3<KeyType, RealType, ElementTag>(
            s, s.d_gradPx, s.d_gradPy, s.d_gradPz);
    }

    // Node coords for the Barth-Jespersen midpoint reconstruction. Always valid
    // (read by the BJ flux branch and the limiter); cheap to fetch.
    const auto& d_bjNodeX = s.domain.getNodeX();
    const auto& d_bjNodeY = s.domain.getNodeY();
    const auto& d_bjNodeZ = s.domain.getNodeZ();

    // ---- Barth-Jespersen pre-pass: node velocity gradients ----
    // Only when advScheme==BarthJespersen; skew/upwind pay zero cost. Build the
    // per-component node gradient grad(q) = (1/V_i) sum_f q_face*A_f exactly like
    // the grad(p) block above (SCS-face scatter -> reverse-halo SUM -> periodic
    // sum -> normalize by lumped mass), then forward-exchange so ghosts carry the
    // owner gradient. The per-component min/max and limiter phi are computed
    // inside runPredictor (they depend on qN and this gradient).
    if (bjMode)
    {
        // Lazy-allocate the 9 gradient fields + 3 scratch buffers once.
        if (s.d_gradUx.size() != s.nodeCount)
        {
            s.d_gradUx.resize(s.nodeCount); s.d_gradUy.resize(s.nodeCount); s.d_gradUz.resize(s.nodeCount);
            s.d_gradVx.resize(s.nodeCount); s.d_gradVy.resize(s.nodeCount); s.d_gradVz.resize(s.nodeCount);
            s.d_gradWx.resize(s.nodeCount); s.d_gradWy.resize(s.nodeCount); s.d_gradWz.resize(s.nodeCount);
            s.d_bjQmin.resize(s.nodeCount); s.d_bjQmax.resize(s.nodeCount); s.d_bjPhi.resize(s.nodeCount);
        }

        // Green-Gauss node gradient grad(q) = (1/V_i) sum_f q_face*A_f, same
        // pattern as grad(p): SCS-face scatter -> reverse-halo SUM -> periodic
        // sum -> normalize by lumped mass -> forward-exchange so ghosts carry the
        // owner gradient. (An edge-stencil LSQ gradient was tried; it did NOT
        // change the multi-rank behaviour -- the blow-up was an explicit-BJ+EXT2
        // CFL limit, not a gradient-stencil issue -- so the cheaper Green-Gauss
        // form is kept.)
        auto computeNodeGradient = [&] (cstone::DeviceVector<RealType>& qField,
                                        cstone::DeviceVector<RealType>& gx,
                                        cstone::DeviceVector<RealType>& gy,
                                        cstone::DeviceVector<RealType>& gz)
        {
            cstone::DeviceVector<RealType> d_gxAcc(s.nodeCount, RealType(0));
            cstone::DeviceVector<RealType> d_gyAcc(s.nodeCount, RealType(0));
            cstone::DeviceVector<RealType> d_gzAcc(s.nodeCount, RealType(0));
            if (eBlocks > 0)
            {
                computeGradientPerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    qField.data(),
                    s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                    d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
                    startElem, numLocal);
                cudaDeviceSynchronize();
            }
            s.domain.reverseExchangeNodeHaloAdd(d_gxAcc);
            s.domain.reverseExchangeNodeHaloAdd(d_gyAcc);
            s.domain.reverseExchangeNodeHaloAdd(d_gzAcc);
            maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_gxAcc);
            maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_gyAcc);
            maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_gzAcc);
            normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
                s.d_massNode.data(),
                s.d_node_to_dof.data(), d_nodeOwnership.data(),
                gx.data(), gy.data(), gz.data(),
                s.nodeCount);
            cudaDeviceSynchronize();
            // Ghosts must carry the owner gradient: the BJ flux branch reads
            // grad_q[up] where up may be a ghost node on a rank boundary.
            s.domain.exchangeNodeHalo(gx);
            s.domain.exchangeNodeHalo(gy);
            s.domain.exchangeNodeHalo(gz);
        };

        computeNodeGradient(s.d_u, s.d_gradUx, s.d_gradUy, s.d_gradUz);
        computeNodeGradient(s.d_v, s.d_gradVx, s.d_gradVy, s.d_gradVz);
        computeNodeGradient(s.d_w, s.d_gradWx, s.d_gradWy, s.d_gradWz);
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
                              cstone::DeviceVector<RealType>& qTarget,
                              RealType bodyForce,
                              cstone::DeviceVector<RealType>& gradQx,
                              cstone::DeviceVector<RealType>& gradQy,
                              cstone::DeviceVector<RealType>& gradQz)
    {
        // Barth-Jespersen per-component limiter: neighbor min/max over edge
        // neighbors (owned scatter + reverse MIN/MAX fold for the cross-rank
        // ring), then phi. Must finish + be published to ghosts BEFORE the flux
        // kernel reads phi[up]/grad[up] (up may be a ghost on a rank boundary).
        if (bjMode)
        {
            bjSeedMinMaxKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                qN.data(), s.d_bjQmin.data(), s.d_bjQmax.data(), s.nodeCount);
            cudaDeviceSynchronize();
            // Scatter over OWNED elements only. An owned boundary node's full
            // edge-neighbor ring is completed across ranks by the reverse
            // min/max fold below -- exactly as the gradient/divergence use the
            // reverse SUM fold. (Looping the cstone element halo here would
            // STILL miss corner-only off-rank neighbors, so it is both
            // insufficient and a double-count risk; the reverse fold is the
            // correct mechanism.)
            if (eBlocks > 0)
            {
                bjNeighborMinMaxScatterKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    qN.data(), s.d_bjQmin.data(), s.d_bjQmax.data(),
                    startElem, numLocal);
                cudaDeviceSynchronize();
            }
            // Reduce each owner's bounds across all ranks holding the node as a
            // ghost (MIN for qmin, MAX for qmax) -- this is the cross-rank
            // completion the forward exchange cannot do -- THEN publish the
            // complete owner bounds to ghosts. Without this fold the limiter is
            // too permissive at rank boundaries and the reconstruction overshoots
            // (multi-rank blow-up).
            s.domain.reverseExchangeNodeHaloMin(s.d_bjQmin);
            s.domain.reverseExchangeNodeHaloMax(s.d_bjQmax);
            s.domain.exchangeNodeHalo(s.d_bjQmin);
            s.domain.exchangeNodeHalo(s.d_bjQmax);

            bjSeedPhiKernel<RealType><<<nodeBlocks, s.blockSize>>>(s.d_bjPhi.data(), s.nodeCount);
            cudaDeviceSynchronize();
            if (eBlocks > 0)
            {
                bjLimiterScatterKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    qN.data(), s.d_bjQmin.data(), s.d_bjQmax.data(),
                    gradQx.data(), gradQy.data(), gradQz.data(),
                    d_bjNodeX.data(), d_bjNodeY.data(), d_bjNodeZ.data(),
                    s.d_bjPhi.data(), startElem, numLocal);
                cudaDeviceSynchronize();
            }
            // phi is a per-node MIN over incident faces; an owned boundary node's
            // off-rank faces are folded in with reverse-MIN, then published.
            s.domain.reverseExchangeNodeHaloMin(s.d_bjPhi);
            s.domain.exchangeNodeHalo(s.d_bjPhi);
        }

        // Compute current-step advection into advN.
        if (advN.size() != s.nodeCount) advN.resize(s.nodeCount);
        thrust::fill(thrust::device_pointer_cast(advN.data()),
                     thrust::device_pointer_cast(advN.data() + s.nodeCount),
                     RealType(0));
        if (eBlocks > 0)
        {
            explicitAdvectionFluxScatterPerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7,
                s.d_u.data(), s.d_v.data(), s.d_w.data(),
                qN.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                advN.data(), startElem, numLocal,
                advMode,
                s.d_bjPhi.data(),
                gradQx.data(), gradQy.data(), gradQz.data(),
                d_bjNodeX.data(), d_bjNodeY.data(), d_bjNodeZ.data());
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(advN);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, advN);

        // DIAGNOSTIC: owned-sum of the (folded) advection term. This is a
        // rank-count-INVARIANT physical quantity if the advection scatter +
        // reverse-fold is cross-rank-consistent. If sum(advN[owned]) differs on
        // 1 vs 4 ranks at the same step, the advection -- the unmeasured
        // predictor input -- is where the cross-rank inconsistency (div(u*) 17x)
        // enters. Owned-only sum; compare 1-rank vs 4-rank at the same step.
        if (g_nsDebugStepsLeft > 0)
        {
            const uint8_t* ownPtr = s.ownershipMap().data();
            const int* dofPtr     = s.d_node_to_dof.data();
            const RealType* aP    = advN.data();
            double locSum = thrust::transform_reduce(thrust::device,
                thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(s.nodeCount),
                [ownPtr, dofPtr, aP] __device__ (size_t i) -> double {
                    if (ownPtr[i] != 1 || dofPtr[i] < 0) return 0.0;
                    return double(aP[i]);
                }, 0.0, thrust::plus<double>());
            double gSum = 0;
            MPI_Allreduce(&locSum, &gSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (s.rank == 0)
                std::cout << "    [advN-sum] sum(advN[owned])=" << std::scientific
                          << std::setprecision(8) << gSum << std::defaultfloat
                          << " (rank-INVARIANT if advection is cross-rank-consistent)\n";
        }

        if (s.useBdf2 && s.bdfStep >= 1 && advNm1.size() == s.nodeCount)
        {
            applyPredictorBdf2PerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                qStar.data(), qN.data(), qNm1.data(),
                advN.data(), advNm1.data(),
                gradPnq.data(), s.d_mass.data(), qTarget.data(),
                s.d_isBdryDof.data(), s.d_node_to_dof.data(),
                d_nodeOwnership.data(),
                dt, invRho, bodyForce, s.nodeCount, s.numOwnedDofs);
        }
        else
        {
            applyPredictorPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
                qStar.data(), qN.data(), advN.data(),
                gradPnq.data(), s.d_mass.data(), s.d_massNode.data(), qTarget.data(),
                s.d_isBdryDof.data(), s.d_node_to_dof.data(),
                d_nodeOwnership.data(),
                dt, invRho, bodyForce, s.nodeCount, s.numOwnedDofs);
        }
        cudaDeviceSynchronize();
    };

    runPredictor(s.d_u, s.d_u_nm1, s.d_uStar, s.d_advU_n, s.d_advU_nm1, s.d_gradPx, s.d_uTarget, s.bodyForceX,
                 s.d_gradUx, s.d_gradUy, s.d_gradUz);
    runPredictor(s.d_v, s.d_v_nm1, s.d_vStar, s.d_advV_n, s.d_advV_nm1, s.d_gradPy, s.d_vTarget, s.bodyForceY,
                 s.d_gradVx, s.d_gradVy, s.d_gradVz);
    runPredictor(s.d_w, s.d_w_nm1, s.d_wStar, s.d_advW_n, s.d_advW_nm1, s.d_gradPz, s.d_wTarget, s.bodyForceZ,
                 s.d_gradWx, s.d_gradWy, s.d_gradWz);

    // Sync ghosts of q* so the implicit RHS / CG warm-start read correct ghosts.
    s.domain.exchangeNodeHalo(s.d_uStar);
    s.domain.exchangeNodeHalo(s.d_vStar);
    s.domain.exchangeNodeHalo(s.d_wStar);

    // Periodic: broadcast master->slave so u* is identical on both sides of
    // the periodic boundary. Without this, the implicit diffusion RHS sees
    // a non-periodic u* and produces a non-periodic u**, propagating drift.
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic && s.periodicMap)
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void runImplicitDiffusionStep(NSStepper<KeyType, RealType, ElementTag>& s, RealType dt)
{
    const auto& d_nodeOwnership = s.ownershipMap();
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
            s.d_mass.data(), invDt, b.data(), s.nodeCount, s.numOwnedDofs);
        cudaDeviceSynchronize();

        enforceBcRhsFromTargetKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            s.d_isBdryDof.data(), qTarget.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            b.data(), s.nodeCount, s.numOwnedDofs);
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

        // Path B: cross-rank slave rows are now Dirichlet-identity. Zero
        // b[slave] AND x_warm[slave] so the initial residual r[slave] = 0
        // throughout CG (the Ap[slave] = 1*p[slave] term plus the
        // periodicBroadcastDofKernel skip-mask keep p[slave]=0 every iter,
        // so r[slave] stays 0). Without this, r[slave] = -nu*Kdiag*qStar
        // would freeze a non-zero floor that prevents ||r||/||b|| from ever
        // hitting tolerance. Post-solve crossRankPeriodicBroadcastDof in
        // solveOneComponent restores x[slave] := x[master_owned_D].
        if (s.numRanks > 1
            && s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
            && s.periodicMap != nullptr
            && !s.periodicMap->cross_.d_sendOwnedSlaveIds_.empty())
        {
            int nSend = int(s.periodicMap->cross_.d_sendOwnedSlaveIds_.size());
            int blk = (nSend + s.blockSize - 1) / s.blockSize;
            zeroDofAtCrossRankSlavesKernel<RealType><<<blk, s.blockSize>>>(
                s.periodicMap->cross_.d_sendOwnedSlaveIds_.data(),
                s.d_node_to_dof.data(),
                nSend, s.numOwnedDofs, b.data());
            // xVec is sized numTotalDofs but the slave dofs are in [0, numOwnedDofs)
            zeroDofAtCrossRankSlavesKernel<RealType><<<blk, s.blockSize>>>(
                s.periodicMap->cross_.d_sendOwnedSlaveIds_.data(),
                s.d_node_to_dof.data(),
                nSend, s.numOwnedDofs, xVec.data());
            cudaDeviceSynchronize();
        }

        return solveOneComponent<KeyType, RealType, ElementTag>(
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
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic && s.periodicMap)
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void runPressureSolveStep(NSStepper<KeyType, RealType, ElementTag>& s, RealType dt, RealType rho)
{
    const auto& d_nodeOwnership = s.ownershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    // connPtrs gives NodesPerElem pointers; c4..c7 stay nullptr for tet (the
    // templated Phase-1 kernels read c0..c3 only when ElementTag==TetTag).
    auto cp = connPtrs<ElementTag, KeyType>(d_conn);
    const KeyType* c0 = cp[0]; const KeyType* c1 = cp[1];
    const KeyType* c2 = cp[2]; const KeyType* c3 = cp[3];
    const KeyType* c4 = nullptr; const KeyType* c5 = nullptr;
    const KeyType* c6 = nullptr; const KeyType* c7 = nullptr;
    if constexpr (std::is_same_v<ElementTag, HexTag>) { c4 = cp[4]; c5 = cp[5]; c6 = cp[6]; c7 = cp[7]; }
    // OWNED-only per-element scatter (proven-green pattern, commit 5da5d6a).
    // Looping owned+halo here would DOUBLE-COUNT shared rank-boundary faces
    // because the halo element on this rank is the same physical element the
    // owning rank already scatters; reverseExchangeNodeHaloAdd then folds the
    // ghost-side contribution back to owner, giving 2x at shared faces and 4x
    // at shared corners.
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
        // Rhie-Chow corrected divergence: adds the compact pressure-velocity
        // coupling term that suppresses the checkerboard pressure mode on
        // Q1-Q1 co-located grids. Standard CVFEM/co-located-FV stabilization
        // (Nalu-Wind, OpenFOAM). Element-generic via scsLR<ElementTag> +
        // ElemTraits<ElementTag>::ScsPerElem. tau auto = dt/rho.
        bool useRC = s.useRhieChow;
        if (useRC)
        {
            RealType tauRC = s.rhieChowTau;
            if (tauRC <= 0) tauRC = dt / rho;
            const auto& d_x = s.domain.getNodeX();
            const auto& d_y = s.domain.getNodeY();
            const auto& d_z = s.domain.getNodeZ();
            computeDivergenceRhieChowKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
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
            computeDivergencePerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7,
                s.d_uStarStar.data(), s.d_vStarStar.data(), s.d_wStarStar.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                d_divAccNode.data(), startElem, numLocal);
        }
        cudaDeviceSynchronize();
    }
    s.domain.reverseExchangeNodeHaloAdd(d_divAccNode);
    maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_divAccNode);

    // Source-of-NaN trace: count non-finite divAccNode over ALL owned nodes and
    // over the owned-BOUNDARY subset separately. The interior diagnostic below
    // (maxOwnedInteriorAbs) excludes boundary DOFs, so a boundary NaN is
    // invisible to it -- this catches it. Gated behind env, free in production.
    if (std::getenv("MARS_DDT_RHS_TRACE"))
    {
        const uint8_t* ownPtr = s.domain.getNodeOwnershipMap().data();
        const int*     dofPtr = s.d_node_to_dof.data();
        const uint8_t* bndPtr = s.d_isBdryDof.data();
        const RealType* dP    = d_divAccNode.data();
        size_t nN = s.nodeCount;
        long long locBadAll = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(nN),
            [ownPtr, dofPtr, dP] __device__ (size_t i) -> long long {
                if (ownPtr[i] != 1 || dofPtr[i] < 0) return 0LL;
                return isfinite(static_cast<double>(dP[i])) ? 0LL : 1LL;
            }, 0LL, thrust::plus<long long>());
        long long locBadBnd = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(nN),
            [ownPtr, dofPtr, bndPtr, dP] __device__ (size_t i) -> long long {
                if (ownPtr[i] != 1) return 0LL;
                int dof = dofPtr[i];
                if (dof < 0 || !bndPtr[dof]) return 0LL;
                return isfinite(static_cast<double>(dP[i])) ? 0LL : 1LL;
            }, 0LL, thrust::plus<long long>());
        long long gBadAll = 0, gBadBnd = 0;
        MPI_Allreduce(&locBadAll, &gBadAll, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&locBadBnd, &gBadBnd, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (s.rank == 0)
            std::cout << "  [DDT-rhs] post-divexch divAccNode owned: nonfinite(all)="
                      << gBadAll << " nonfinite(boundary)=" << gBadBnd << "\n";
    }

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
            d_divAccNode.data(), s.d_massNode.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_divUStar.data(), s.nodeCount);
        cudaDeviceSynchronize();
        s.lastDivMaxPre = maxOwnedInteriorAbs<KeyType, RealType, ElementTag>(s, s.d_divUStar);
    }

    if (s.pressureSolve == PressureSolveKind::K)
    {
        // K-path: Galerkin Laplacian (K is for -Δu=f, so RHS = -coef * divAcc).
        thrust::fill(thrust::device_pointer_cast(b.data()),
                     thrust::device_pointer_cast(b.data() + s.numOwnedDofs), RealType(0));
        RealType coef = rho * invDt;
        buildPressureRhsKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_divAccNode.data(), s.d_node_to_dof.data(),
            d_nodeOwnership.data(), coef, b.data(), s.nodeCount, s.numOwnedDofs);
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
        if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic && s.numOwnedDofs > 0)
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
        s.lastPressureIters = solveOneComponent<KeyType, RealType, ElementTag>(s, b, xVec, s.d_phi, s.Apre);
    }
    else  // DDT: (D M^{-1} D^T) phi = -(rho/dt) D u**
    {
        // DDT pressure ALWAYS uses the matrix-free path (the else branch
        // below). The assembled DDT operator s.AddT holds POSITIVE off-
        // diagonals at boundary rows (D M^-1 D^T is not an M-matrix),
        // which BoomerAMG's classical strength-of-connection rejects --
        // every Hypre Setup attempt on s.AddT returned HYPRE_ERROR_GENERIC
        // (cg_iter_p=FAIL every step on cube16 cavity). The matrix-free
        // CG path with Jacobi preconditioning (see solvePressureDDT) is
        // the working production route. The Hypre wrappers still serve
        // velocity solves and the K-path pressure solve, both of which
        // are well-conditioned M-matrices where AMG works as designed.
        //
        // The diagnostic env MARS_DDT_USE_ASSEMBLED_CG=1 still runs the
        // assembled-CG path so we can verify the assembled matrix equals
        // the matrix-free apply bit-for-bit (useful when changing the
        // DDT assembler). It only fires when --solver=cg.
        const char* useAsmEv = std::getenv("MARS_DDT_USE_ASSEMBLED_CG");
        bool useAssembledCG = (useAsmEv && std::string(useAsmEv) != "0")
                              && (s.solverKind == SolverKind::CG);
        if (useAssembledCG)
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
                d_nodeOwnership.data(), coef, b.data(), s.nodeCount, s.numOwnedDofs);
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
            // RHS finiteness check. A non-finite entry here (e.g. a boundary
            // node whose divergence went bad upstream) would, if fed into the
            // SUM-Allreduce below, poison b on EVERY rank and make Hypre return
            // a generic error. The masked reduce localizes the source instead
            // of smearing NaN across the whole vector. Gated behind env so it
            // costs nothing in production.
            if (s.numOwnedDofs > 0 && std::getenv("MARS_DDT_RHS_TRACE"))
            {
                auto bp = thrust::device_pointer_cast(b.data());
                long long localBad = thrust::transform_reduce(thrust::device,
                    bp, bp + s.numOwnedDofs,
                    [] __device__ (RealType v) -> long long {
                        return isfinite(static_cast<double>(v)) ? 0LL : 1LL;
                    }, 0LL, thrust::plus<long long>());
                RealType localMax = thrust::transform_reduce(thrust::device,
                    bp, bp + s.numOwnedDofs,
                    [] __device__ (RealType v) -> RealType {
                        double d = static_cast<double>(v);
                        return isfinite(d) ? static_cast<RealType>(fabs(d)) : RealType(0);
                    }, RealType(0), thrust::maximum<RealType>());
                long long gBad = 0; RealType gMax = 0;
                MPI_Datatype mpiR = std::is_same<RealType, double>::value ? MPI_DOUBLE : MPI_FLOAT;
                MPI_Allreduce(&localBad, &gBad, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&localMax, &gMax, 1, mpiR, MPI_MAX, MPI_COMM_WORLD);
                if (s.rank == 0)
                    std::cout << "  [DDT-rhs] post-build owned-b: nonfinite=" << gBad
                              << " |b|max(finite)=" << gMax << "\n";
            }
            // Null-space projection: subtract mean(b)*1 so b lies in range(A).
            // Required ONLY for the pure-Neumann (periodic) operator, whose
            // constant null mode makes A phi = b solvable iff sum(b)=0. For
            // cavity (single Dirichlet pin) and channel/pump (Dirichlet outflow
            // face) the pin/face already removes the constant mode, so the K
            // path does NOT project there -- we match it exactly. Subtracting a
            // global mean from a pinned system is both unnecessary and, when an
            // upstream non-finite entry sneaks in, the mechanism that turned all
            // of b into NaN on every rank. The sum still uses a finiteness mask
            // as defense in depth.
            if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic
                && s.numOwnedDofs > 0)
            {
                auto bp = thrust::device_pointer_cast(b.data());
                RealType localSum = thrust::transform_reduce(thrust::device,
                    bp, bp + s.numOwnedDofs,
                    [] __device__ (RealType v) -> RealType {
                        double d = static_cast<double>(v);
                        return isfinite(d) ? static_cast<RealType>(d) : RealType(0);
                    }, RealType(0), thrust::plus<RealType>());
                RealType globalSum = 0;
                long long localN = s.numOwnedDofs, globalN = 0;
                MPI_Datatype mpiR = std::is_same<RealType, double>::value
                                    ? MPI_DOUBLE : MPI_FLOAT;
                MPI_Allreduce(&localSum, &globalSum, 1, mpiR, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&localN,   &globalN,   1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
                if (globalN > 0)
                {
                    RealType mean = globalSum / RealType(globalN);
                    thrust::transform(thrust::device, bp, bp + s.numOwnedDofs, bp,
                                      [mean] __device__ (RealType v) { return v - mean; });
                }
            }
            // Re-enforce the pin (RHS at pin must be 0).
            if (s.pressurePinDof >= 0)
            {
                enforcePinRhsKernel<RealType><<<1, 1>>>(s.pressurePinDof, b.data());
                cudaDeviceSynchronize();
            }
            cudaMemset(xVec.data(), 0, s.numTotalDofs * sizeof(RealType));
            // DDT operator (D M^-1 D^T) is SPSD with constant null space; Hypre
            // PCG rejects it (error 1). GMRES tolerates it. Hint passes through
            // to solveOneComponent which selects HypreGMRESSolver wrapper.
            s.lastPressureIters = solveOneComponent<KeyType, RealType, ElementTag>(
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
            // DOF-weighted mean (partner table) so collapsed periodic DOFs are
            // not over-counted -> RHS is exactly mean-zero in DOF space.
            if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic)
            {
                const int* partnerPtr = s.periodicMap
                                        ? s.periodicMap->d_periodicPartner.data() : nullptr;
                mars::fem::removeMean<RealType>(s.domain, d_bNode, MPI_COMM_WORLD, partnerPtr);
            }
            s.lastPressureIters = solvePressureDDT<KeyType, RealType, ElementTag>(s, d_bNode, s.d_phi);
        }
    }

    // Periodic: pure-Neumann pressure has a constant-mode null space. Subtract
    // the DOF-space mean so phi is uniquely defined (and the constant doesn't
    // drift step over step). Partner table => collapsed DOFs counted once.
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic)
    {
        const int* partnerPtr = s.periodicMap
                                ? s.periodicMap->d_periodicPartner.data() : nullptr;
        mars::fem::removeMean<RealType>(s.domain, s.d_phi, MPI_COMM_WORLD, partnerPtr);
    }

    // Sync ghosts of phi -- needed by both the gradient (face donors) and the
    // pressure update on owned nodes that read ghost slots only for VTU.
    s.domain.exchangeNodeHalo(s.d_phi);

    // Periodic: broadcast phi master->slave so the corrector's grad(phi) reads
    // a periodic field. Without this, grad(phi) is asymmetric at the periodic
    // boundary, the corrector produces a non-divergence-free u^{n+1}, and the
    // next pressure source is correspondingly biased.
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic && s.periodicMap)
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
// Divergence of the Rhie-Chow-CORRECTED face flux, not the plain nodal flux.
// This is the operator the pressure solve actually drives to zero when RC is on.
// The plain divMaxAndRmsOwned still "sees" the checkerboard mode that RC has
// decoupled from the flux, so it can stay O(1) even when RC works -- use THIS
// to tell whether RC is doing its job. tau matches the solve (s.rhieChowTau or
// dt/rho). Returns roundoff-level max when the projection is consistent.
template<typename KeyType, typename RealType, typename ElementTag>
inline void divMaxRhieChowOwned(NSStepper<KeyType, RealType, ElementTag>& s,
                                RealType dt, RealType rho,
                                RealType& outMax, RealType& outRms)
{
    const auto& d_own  = s.domain.getNodeOwnershipMap();
    const auto& d_conn = s.domain.getElementToNodeConnectivity();
    auto cp = connPtrs<ElementTag, KeyType>(d_conn);
    const KeyType* c0 = cp[0]; const KeyType* c1 = cp[1];
    const KeyType* c2 = cp[2]; const KeyType* c3 = cp[3];
    const KeyType* c4 = nullptr; const KeyType* c5 = nullptr;
    const KeyType* c6 = nullptr; const KeyType* c7 = nullptr;
    if constexpr (std::is_same_v<ElementTag, HexTag>) { c4 = cp[4]; c5 = cp[5]; c6 = cp[6]; c7 = cp[7]; }
    size_t startElem = s.domain.startIndex();
    size_t numLocal  = s.domain.localElementCount();
    int eBlocks = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
    int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;

    RealType tauRC = s.rhieChowTau;
    if (tauRC <= 0) tauRC = dt / rho;
    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();

    cstone::DeviceVector<RealType> d_divAcc(s.nodeCount, RealType(0));
    if (eBlocks > 0) {
        computeDivergenceRhieChowKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
            c0, c1, c2, c3, c4, c5, c6, c7,
            s.d_u.data(), s.d_v.data(), s.d_w.data(),
            s.d_p.data(),
            s.d_gradPx.data(), s.d_gradPy.data(), s.d_gradPz.data(),
            d_x.data(), d_y.data(), d_z.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            tauRC, d_divAcc.data(), startElem, numLocal);
        cudaDeviceSynchronize();
    }
    s.domain.reverseExchangeNodeHaloAdd(d_divAcc);
    cstone::DeviceVector<RealType> d_divNorm(s.nodeCount, RealType(0));
    normalizeDivergencePerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        d_divAcc.data(), s.d_massNode.data(),
        s.d_node_to_dof.data(), d_own.data(),
        d_divNorm.data(), s.nodeCount);
    cudaDeviceSynchronize();
    outMax = maxOwnedInteriorAbs<KeyType, RealType, ElementTag>(s, d_divNorm);
    outRms = rmsOwnedInterior1<KeyType, RealType, ElementTag>(s, d_divNorm);
}

//   p^{n+1} = p^n + phi
// grad(phi) computed once via the same B.4 kernel, applied per-component.
// Writes s.lastDivMax for monitoring.
// -------------------------------------------------------------------------
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void runCorrectorStep(NSStepper<KeyType, RealType, ElementTag>& s, RealType dt, RealType rho)
{
    const auto& d_nodeOwnership = s.ownershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    // connPtrs gives NodesPerElem pointers; c4..c7 stay nullptr for tet (the
    // templated Phase-1 kernels read c0..c3 only when ElementTag==TetTag).
    auto cp = connPtrs<ElementTag, KeyType>(d_conn);
    const KeyType* c0 = cp[0]; const KeyType* c1 = cp[1];
    const KeyType* c2 = cp[2]; const KeyType* c3 = cp[3];
    const KeyType* c4 = nullptr; const KeyType* c5 = nullptr;
    const KeyType* c6 = nullptr; const KeyType* c7 = nullptr;
    if constexpr (std::is_same_v<ElementTag, HexTag>) { c4 = cp[4]; c5 = cp[5]; c6 = cp[6]; c7 = cp[7]; }
    // OWNED-only per-element scatter (proven-green pattern, commit 5da5d6a).
    // Looping owned+halo here would DOUBLE-COUNT shared rank-boundary faces
    // because the halo element on this rank is the same physical element the
    // owning rank already scatters; reverseExchangeNodeHaloAdd then folds the
    // ghost-side contribution back to owner, giving 2x at shared faces and 4x
    // at shared corners.
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
                computeGradientPerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                    c0, c1, c2, c3, c4, c5, c6, c7,
                    s.d_phi.data(),
                    s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                    d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
                    startElem, numLocal);
            }
            else
            {
                applyDivTransposePerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
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
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_gxAcc);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_gyAcc);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_gzAcc);
        normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
            s.d_massNode.data(),
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
        s.lastGradPhiRms = rmsOwnedInterior3<KeyType, RealType, ElementTag>(
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
            d_nodeOwnership.data(), dtEff, invRho, s.nodeCount, s.numOwnedDofs);
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
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic)
    {
        // DOF-weighted (partner table): a same-rank periodic slave aliases its
        // master's DOF, so counting it would weight the re-anchor mean by the
        // face multiplicity and re-bias p every step -- the same gauge drift
        // this re-anchor is meant to kill. Counted once per DOF instead.
        const int* partnerPtr = s.periodicMap
                                ? s.periodicMap->d_periodicPartner.data() : nullptr;
        mars::fem::removeMean<RealType>(s.domain, s.d_p, MPI_COMM_WORLD, partnerPtr);
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
    if (s.bcKind == NSStepper<KeyType, RealType, ElementTag>::BCKind::Periodic && s.periodicMap)
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
            computeDivergencePerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7,
                s.d_u.data(), s.d_v.data(), s.d_w.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                d_divAccPost.data(), startElem, numLocal);
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(d_divAccPost);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, d_divAccPost);
        normalizeDivergencePerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_divAccPost.data(), s.d_massNode.data(),
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
        s.lastDivRms = rmsOwnedInterior1<KeyType, RealType, ElementTag>(s, d_divNorm);

        // DIAGNOSTIC: split div RMS into RANK-BOUNDARY-owned vs INTERIOR-owned
        // nodes. If the multi-rank divergence excess is boundary-localized, the
        // boundary RMS >> interior RMS => a cross-rank halo-completeness bug at
        // shared nodes. If both are inflated equally => a global issue. Builds a
        // per-node boundary mask from the node-halo topology's sendNodeIds_
        // (owned nodes shared with a peer). Runs only when g_nsDebugStepsLeft>0.
        if (s.numRanks > 1 && g_nsDebugStepsLeft > 0)
        {
            cstone::DeviceVector<uint8_t> d_isBndNode(s.nodeCount, uint8_t(0));
            const auto& topo = s.domain.getNodeHaloTopology();
            size_t nSend = topo.sendNodeIds_.size();
            if (nSend > 0)
            {
                const int* snd = topo.sendNodeIds_.data();
                uint8_t* bnd = d_isBndNode.data();
                size_t nn = s.nodeCount;
                thrust::for_each(thrust::device,
                    thrust::counting_iterator<size_t>(0),
                    thrust::counting_iterator<size_t>(nSend),
                    [snd, bnd, nn] __device__ (size_t k) {
                        int nid = snd[k];
                        if (nid >= 0 && size_t(nid) < nn) bnd[nid] = uint8_t(1);
                    });
                cudaDeviceSynchronize();
            }
            const uint8_t* ownPtr = s.ownershipMap().data();
            const int* dofPtr     = s.d_node_to_dof.data();
            const uint8_t* bndPtr = d_isBndNode.data();
            const RealType* dP    = d_divNorm.data();
            // (sumSq, count) for boundary and interior separately.
            auto accum = [&] (bool wantBoundary) {
                RealType locSq = thrust::transform_reduce(thrust::device,
                    thrust::counting_iterator<size_t>(0),
                    thrust::counting_iterator<size_t>(s.nodeCount),
                    [ownPtr, dofPtr, bndPtr, dP, wantBoundary] __device__ (size_t i) -> RealType {
                        if (ownPtr[i] != 1) return RealType(0);
                        if (dofPtr[i] < 0)  return RealType(0);
                        bool isB = (bndPtr[i] != 0);
                        if (isB != wantBoundary) return RealType(0);
                        RealType q = dP[i]; return q * q;
                    }, RealType(0), thrust::plus<RealType>());
                long long locCnt = thrust::transform_reduce(thrust::device,
                    thrust::counting_iterator<size_t>(0),
                    thrust::counting_iterator<size_t>(s.nodeCount),
                    [ownPtr, dofPtr, bndPtr, wantBoundary] __device__ (size_t i) -> long long {
                        if (ownPtr[i] != 1) return 0LL;
                        if (dofPtr[i] < 0)  return 0LL;
                        bool isB = (bndPtr[i] != 0);
                        return (isB == wantBoundary) ? 1LL : 0LL;
                    }, 0LL, thrust::plus<long long>());
                auto mt = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
                RealType gSq = 0; long long gCnt = 0;
                MPI_Allreduce(&locSq, &gSq, 1, mt, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&locCnt, &gCnt, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
                RealType rms = (gCnt > 0) ? std::sqrt(gSq / RealType(gCnt)) : RealType(0);
                return std::make_pair(rms, gCnt);
            };
            auto bRms = accum(true);
            auto iRms = accum(false);
            if (s.rank == 0)
                std::cout << "    [div-split] boundary RMS=" << std::scientific << std::setprecision(3)
                          << bRms.first << " (" << bRms.second << " nodes)  interior RMS="
                          << iRms.first << " (" << iRms.second << " nodes)  ratio="
                          << (iRms.first > 0 ? bRms.first / iRms.first : RealType(0))
                          << std::defaultfloat << "\n";
        }
    }
    // When RC is on, also report the divergence of the RC-corrected flux -- the
    // operator the pressure solve drives to zero. The plain lastDivMax above can
    // stay O(1) (it still resolves the checkerboard mode RC decoupled from the
    // flux), so lastDivRC is the honest "is RC working" signal.
    if (s.useRhieChow)
    {
        RealType rcMax = 0, rcRms = 0;
        divMaxRhieChowOwned<KeyType, RealType, ElementTag>(s, dt, rho, rcMax, rcRms);
        s.lastDivRC = rcMax;
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
inline RealType keOwned(NSStepper<KeyType, RealType, ElementTag>& s,
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
inline RealType sumOwned(NSStepper<KeyType, RealType, ElementTag>& s,
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
inline void divMaxAndRmsOwned(NSStepper<KeyType, RealType, ElementTag>& s,
                              const cstone::DeviceVector<RealType>& d_u,
                              const cstone::DeviceVector<RealType>& d_v,
                              const cstone::DeviceVector<RealType>& d_w,
                              RealType& outMax, RealType& outRms)
{
    const auto& d_own = s.domain.getNodeOwnershipMap();
    const auto& d_conn = s.domain.getElementToNodeConnectivity();
    // connPtrs gives NodesPerElem pointers; c4..c7 stay nullptr for tet (the
    // templated Phase-1 kernels read c0..c3 only when ElementTag==TetTag).
    auto cp = connPtrs<ElementTag, KeyType>(d_conn);
    const KeyType* c0 = cp[0]; const KeyType* c1 = cp[1];
    const KeyType* c2 = cp[2]; const KeyType* c3 = cp[3];
    const KeyType* c4 = nullptr; const KeyType* c5 = nullptr;
    const KeyType* c6 = nullptr; const KeyType* c7 = nullptr;
    if constexpr (std::is_same_v<ElementTag, HexTag>) { c4 = cp[4]; c5 = cp[5]; c6 = cp[6]; c7 = cp[7]; }
    // OWNED-only per-element scatter (proven-green pattern, commit 5da5d6a).
    // Looping owned+halo would double-count shared rank-boundary faces because
    // reverseExchangeNodeHaloAdd is called right after the kernel.
    size_t startElem = s.domain.startIndex();
    size_t numLocal  = s.domain.localElementCount();
    int eBlocks = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
    int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;

    cstone::DeviceVector<RealType> d_divAcc(s.nodeCount, RealType(0));
    if (eBlocks > 0) {
        computeDivergencePerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
            c0, c1, c2, c3, c4, c5, c6, c7,
            d_u.data(), d_v.data(), d_w.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            d_divAcc.data(), startElem, numLocal);
        cudaDeviceSynchronize();
    }
    s.domain.reverseExchangeNodeHaloAdd(d_divAcc);
    cstone::DeviceVector<RealType> d_divNorm(s.nodeCount, RealType(0));
    normalizeDivergencePerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        d_divAcc.data(), s.d_massNode.data(),
        s.d_node_to_dof.data(), d_own.data(),
        d_divNorm.data(), s.nodeCount);
    cudaDeviceSynchronize();
    outMax = maxOwnedInteriorAbs<KeyType, RealType, ElementTag>(s, d_divNorm);
    outRms = rmsOwnedInterior1<KeyType, RealType, ElementTag>(s, d_divNorm);
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
template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void computeVorticityMagnitudePerNode(NSStepper<KeyType, RealType, ElementTag>& s,
                                      cstone::DeviceVector<RealType>& d_u,
                                      cstone::DeviceVector<RealType>& d_v,
                                      cstone::DeviceVector<RealType>& d_w,
                                      cstone::DeviceVector<RealType>& d_omegaMag)
{
    const auto& d_nodeOwnership = s.ownershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    // connPtrs gives NodesPerElem pointers; c4..c7 stay nullptr for tet (the
    // templated Phase-1 kernels read c0..c3 only when ElementTag==TetTag).
    auto cp = connPtrs<ElementTag, KeyType>(d_conn);
    const KeyType* c0 = cp[0]; const KeyType* c1 = cp[1];
    const KeyType* c2 = cp[2]; const KeyType* c3 = cp[3];
    const KeyType* c4 = nullptr; const KeyType* c5 = nullptr;
    const KeyType* c6 = nullptr; const KeyType* c7 = nullptr;
    if constexpr (std::is_same_v<ElementTag, HexTag>) { c4 = cp[4]; c5 = cp[5]; c6 = cp[6]; c7 = cp[7]; }
    // OWNED-only per-element scatter (proven-green pattern, commit 5da5d6a).
    // Looping owned+halo would double-count shared rank-boundary faces because
    // reverseExchangeNodeHaloAdd is called right after the kernel.
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
            computeGradientPerNodeKernel<KeyType, RealType, ElementTag><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7, qNode,
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                gxAcc.data(), gyAcc.data(), gzAcc.data(),
                startElem, numLocal);
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(gxAcc);
        s.domain.reverseExchangeNodeHaloAdd(gyAcc);
        s.domain.reverseExchangeNodeHaloAdd(gzAcc);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, gxAcc);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, gyAcc);
        maybePeriodicSum<KeyType, RealType, ElementTag>(s, gzAcc);
        normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            gxAcc.data(), gyAcc.data(), gzAcc.data(),
            s.d_massNode.data(),
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

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
void runNsStep(NSStepper<KeyType, RealType, ElementTag>& s, RealType dt, RealType nu, RealType rho)
{
    s.nuCached = nu;
    bool dbg = (g_nsDebugStepsLeft > 0);

    // ENTRY: state of u^n, p^n before the step.
    // IMPORTANT: keOwned / divMaxAndRmsOwned / computeWeightedL2Norm all do an
    // MPI_Allreduce internally. They MUST be called on EVERY rank or the
    // collective deadlocks. Compute on all ranks; print only on rank 0.
    // (maxAbsOwned is a local thrust reduce -- no collective -- so it can stay
    // inside the rank-0 print.)
    if (dbg)
    {
        RealType KE_n      = keOwned<KeyType, RealType, ElementTag>(s, s.d_u, s.d_v, s.d_w);
        RealType div_n_max = 0, div_n_rms = 0;
        divMaxAndRmsOwned<KeyType, RealType, ElementTag>(s, s.d_u, s.d_v, s.d_w, div_n_max, div_n_rms);
        if (s.rank == 0)
            std::cout << "    [ns-dbg ENTRY] KE_n=" << std::scientific << std::setprecision(6) << KE_n
                      << " |u_n|max=" << maxAbsOwned(s.d_u, s.nodeCount)
                      << " |p_n|max=" << maxAbsOwned(s.d_p, s.nodeCount)
                      << " div(u_n)max=" << div_n_max
                      << " div(u_n)rms=" << div_n_rms << "\n";
    }

    runPredictorStep<KeyType, RealType, ElementTag>(s, dt, rho);
    if (dbg)
    {
        RealType KE_star      = keOwned<KeyType, RealType, ElementTag>(s, s.d_uStar, s.d_vStar, s.d_wStar);
        RealType div_star_max = 0, div_star_rms = 0;
        divMaxAndRmsOwned<KeyType, RealType, ElementTag>(s, s.d_uStar, s.d_vStar, s.d_wStar, div_star_max, div_star_rms);
        if (s.rank == 0)
            std::cout << "    [ns-dbg PRED ] KE*=" << KE_star
                      << " |u*|max="    << maxAbsOwned(s.d_uStar, s.nodeCount)
                      << " |gP|max="    << maxAbsOwned(s.d_gradPx, s.nodeCount)
                      << " div(u*)max=" << div_star_max
                      << " div(u*)rms=" << div_star_rms << "\n";
    }

    runImplicitDiffusionStep<KeyType, RealType, ElementTag>(s, dt);
    if (dbg)
    {
        RealType KE_ss      = keOwned<KeyType, RealType, ElementTag>(s, s.d_uStarStar, s.d_vStarStar, s.d_wStarStar);
        RealType div_ss_max = 0, div_ss_rms = 0;
        divMaxAndRmsOwned<KeyType, RealType, ElementTag>(s, s.d_uStarStar, s.d_vStarStar, s.d_wStarStar, div_ss_max, div_ss_rms);
        if (s.rank == 0)
            std::cout << "    [ns-dbg DIFF ] KE**=" << KE_ss
                      << " |u**|max="    << maxAbsOwned(s.d_uStarStar, s.nodeCount)
                      << " div(u**)max=" << div_ss_max
                      << " div(u**)rms=" << div_ss_rms
                      << " cg_uvw="      << s.lastUIters << "/" << s.lastVIters << "/" << s.lastWIters << "\n";
    }

    runPressureSolveStep<KeyType, RealType, ElementTag>(s, dt, rho);
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

    runCorrectorStep<KeyType, RealType, ElementTag>(s, dt, rho);
    if (dbg)
    {
        RealType KE_np1 = keOwned<KeyType, RealType, ElementTag>(s, s.d_u, s.d_v, s.d_w);   // collective: all ranks
        if (s.rank == 0)
        {
            std::cout << "    [ns-dbg CORR ] KE_(n+1)=" << KE_np1
                      << " |u_(n+1)|max=" << maxAbsOwned(s.d_u, s.nodeCount)
                      << " |p_(n+1)|max=" << maxAbsOwned(s.d_p, s.nodeCount)
                      << " div_max=" << s.lastDivMax << "\n";
            std::cout << std::defaultfloat;
        }
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

template<typename KeyType, typename RealType, typename ElementTag = HexTag>
RealType computeWeightedL2Norm(NSStepper<KeyType, RealType, ElementTag>& s,
                               const cstone::DeviceVector<RealType>& d_q)
{
    const auto& d_nodeOwnership = s.ownershipMap();
    cstone::DeviceVector<RealType> d_sq(s.numOwnedDofs, RealType(0));
    const int numOwnedDofs = s.numOwnedDofs;
    thrust::for_each(thrust::device,
                      thrust::counting_iterator<size_t>(0),
                      thrust::counting_iterator<size_t>(s.nodeCount),
                      [nodeToDof = s.d_node_to_dof.data(),
                       ownPtr    = d_nodeOwnership.data(),
                       mass      = s.d_mass.data(),
                       q         = d_q.data(),
                       out       = d_sq.data(),
                       numOwnedDofs] __device__(size_t i)
                      {
                          if (ownPtr[i] != 1) return;
                          int dof = nodeToDof[i];
                          if (dof < 0 || dof >= numOwnedDofs) return;
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
