// Incompressible Navier-Stokes driver, Chorin projection method (Phase C).
//
// Per step, given (u^n, v^n, w^n, p^n) on nodes:
//   1) PREDICTOR (explicit):
//        q* = q^n - dt * [ (u^n.grad) q^n + (1/rho) (grad p^n)_q ]
//      for each velocity component q in {u, v, w}. The advective derivative
//      reuses the per-node + reverse-halo upwind scatter from mars_amr_advdiff
//      (B.2); the pressure-gradient component comes from B.4's SCS gradient.
//   2) IMPLICIT DIFFUSION (per component):
//        (M/dt + nu K) q** = (M/dt) q*
//      Same matrix for all three components (constant in time), assembled ONCE.
//      Three CG/Hypre solves, one per component.
//   3) PRESSURE POISSON:
//        div u** = B.3_divergence(u**, v**, w**)
//        K phi   = (rho/dt) div u**
//      K is the SAME Laplacian (gamma=1), but with Neumann BC on all faces +
//      a single Dirichlet pin to fix the constant-mode null space. SEPARATE
//      matrix from the velocity (BC treatment differs).
//   4) CORRECTOR:
//        q^{n+1} = q** - (dt/rho) (grad phi)_q
//        p^{n+1} = p^n + phi      (incremental Chorin)
//
// Test cases (selected via --bc=):
//   cavity:  Lid-driven cavity, u=lidU on z=zmax, no-slip on the other 5
//            velocity faces. Pressure has Neumann everywhere with a single
//            corner DOF pinned to 0 to remove the constant-mode null space.
//   channel: Channel flow, Dirichlet u=Uinf on x=xmin inflow, no-slip on y/z
//            tunnel walls, natural (Neumann) velocity on x=xmax outflow.
//            Pressure has Dirichlet p=0 on the entire outflow face (this is
//            the physically correct condition for channel flow and removes
//            the null space without fighting the natural inflow-to-outflow
//            pressure gradient).
// Note that with the K-path pressure solve div_max decays slowly (not at
// roundoff) because K from assembleFull is not the discrete D D^T; see the
// PressureSolveKind comment for the experimental DDT alternative.
//
// Reuses the per-node + reverseExchangeNodeHaloAdd scatter pattern validated in
// B.2/B.3/B.4 across all five scatter kernels (mass, advection x3 components,
// gradient, divergence). All velocity/pressure halo exchanges happen between
// steps explicitly; ghost-side staleness is the most common projection-driver
// bug, so we exchange whenever a field changes.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include "backend/distributed/unstructured/fem/mars_perf_counters.hpp"
#include "backend/distributed/unstructured/fem/mars_sparse_matrix.hpp"
#include "backend/distributed/unstructured/solvers/mars_cg_solver.hpp"
#ifdef MARS_ENABLE_HYPRE
#include "backend/distributed/unstructured/solvers/mars_hypre_pcg_solver.hpp"
#endif
#include "backend/distributed/unstructured/amr/mars_amr.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"

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

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

// Sync on lap so wall-clock reflects all-rank work, not rank-0's stragglers.
struct PhaseTimer
{
    using clk = std::chrono::high_resolution_clock;
    clk::time_point t0;
    std::vector<std::pair<std::string, float>> phases;
    bool sync;

    explicit PhaseTimer(bool sync_ = true) : sync(sync_) { reset(); }

    void reset()
    {
        if (sync) { cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD); }
        t0 = clk::now();
    }

    void lap(const std::string& name)
    {
        if (sync) { cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD); }
        auto t1 = clk::now();
        phases.emplace_back(name, std::chrono::duration<float, std::milli>(t1 - t0).count());
        t0 = t1;
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

// Owned-node mass slots -> per-DOF mass array (post reverse halo).
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
    massDof[dof] = massNode[i];
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
// Advection scatter: per-node + reverse-halo upwind. Computes net out-flux at
// each owned node for each velocity component q in {u, v, w}. Reused for all
// three predictor components by passing different q arrays.
// (u.grad)q decomposes as div(u q) - q div(u); on a single time-step with the
// previous-step velocity, the upwind flux mdot * q_upwind is the conservative
// surrogate used by mars_amr_advdiff -- and div(u^n) is not exactly zero (the
// projection is one step behind), so we use the conservative flux form
// directly. Sign matches advdiff: -mdot * q_up on L, +mdot * q_up on R.
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
        RealType mdot = vfx * areaVecX[off] + vfy * areaVecY[off] + vfz * areaVecZ[off];

        RealType q_up = (mdot > RealType(0)) ? q[iL] : q[iR];
        RealType flux = mdot * q_up;

        atomicAdd(&dqdtNode[iL], -flux);
        atomicAdd(&dqdtNode[iR], +flux);
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
        // Boundary value enforced exactly -- the implicit RHS will lift it.
        qStar[i] = qTarget[i];
        return;
    }
    RealType V = lumpedMass[dof];
    // dt * (advection: divide flux sum by V) + dt * (-grad p^n / rho)
    qStar[i] = qN[i] + dt * dqdtNode[i] / V - dt * invRho * gradPnq[i];
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
//   DDT = (EXPERIMENTAL, currently DIVERGES) matrix-free CG on D M^{-1} D^T.
//         The algebra is sound on paper but the implementation cube16/Re=100
//         test diverges in ~25 steps. Diagnosis pending; needs isolation-test
//         of A phi against an analytical phi before re-enabling.
enum class PressureSolveKind { K, DDT };

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
    cstone::DeviceVector<RealType> d_mass;
    cstone::DeviceVector<uint8_t>  d_isBdryDof;
    // Pressure BC: in cavity mode, use a single corner pin (pressurePinDof).
    // In channel mode, use a Dirichlet mask (d_isPressureBdryDof) over the
    // entire outflow face. Exactly one mechanism is active per run.
    int pressurePinDof = -1;   // owned-DOF id on the owning rank; -1 elsewhere (or in channel mode)
    int pressurePinRank = 0;   // global rank that owns the pin
    cstone::DeviceVector<uint8_t> d_isPressureBdryDof;   // size numOwnedDofs; only used in channel mode

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
    // diagPtr (same Laplacian sparsity); values are stored separately.
    cstone::DeviceVector<int> d_rowPtr;
    cstone::DeviceVector<int> d_colInd;
    cstone::DeviceVector<int> d_diagPtr;
    cstone::DeviceVector<RealType> d_valuesVel;
    cstone::DeviceVector<RealType> d_valuesPre;

    // Owned-row SparseMatrix wrappers consumed by CG/Hypre.
    using Matrix = SparseMatrix<int, RealType, cstone::GpuTag>;
    Matrix Avel;
    Matrix Apre;

    // Per-node solution fields. Sized nodeCount; ghost slots refreshed by
    // exchangeNodeHalo after each update.
    cstone::DeviceVector<RealType> d_u, d_v, d_w;            // velocity
    cstone::DeviceVector<RealType> d_uStar, d_vStar, d_wStar;        // after predictor
    cstone::DeviceVector<RealType> d_uStarStar, d_vStarStar, d_wStarStar;  // after diffusion
    cstone::DeviceVector<RealType> d_p;                       // pressure (n)
    cstone::DeviceVector<RealType> d_phi;                     // pressure correction
    cstone::DeviceVector<RealType> d_gradPx, d_gradPy, d_gradPz;
    cstone::DeviceVector<RealType> d_gradPhix, d_gradPhiy, d_gradPhiz;

    // Per-step diagnostics
    int lastPressureIters = 0;
    int lastUIters = 0, lastVIters = 0, lastWIters = 0;
    RealType lastDivMax = 0;

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
    enum class BCKind { Cavity, Channel };
    BCKind bcKind = BCKind::Cavity;
    RealType Uinf = 1;   // channel inflow speed; reuses lidU value via CLI
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
// Setup: builds DOF mapping, sparsity, two matrices, lumped mass, BC mask,
// pressure pin, area vectors. Runs ONCE before the time loop.
// =============================================================================

template<typename KeyType, typename RealType>
void setupNSStepper(NSStepper<KeyType, RealType>& s,
                    RealType nu,
                    RealType dt,
                    CvfemKernelVariant kernelVariant)
{
    PhaseTimer pt;

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
    {
        cstone::DeviceVector<RealType> d_massNode(s.nodeCount, RealType(0));
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
                d_massNode.data(), startElem, numLocal);
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(d_massNode);
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        gatherOwnedNodeMassToDofKernel<RealType><<<nBlocks, s.blockSize>>>(
            d_massNode.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_mass.data(), s.nodeCount);
        cudaDeviceSynchronize();
    }
    pt.lap("lumped mass");

    // Add M/dt to the velocity matrix diagonal -> (M/dt + nu K).
    {
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        addLumpedMassDiagonalKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_mass.data(), s.d_diagPtr.data(), RealType(1) / dt,
            s.d_valuesVel.data(), s.numOwnedDofs);
        cudaDeviceSynchronize();
    }
    pt.lap("add M/dt (velocity)");

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
        else  // Channel
        {
            markChannelBCKernel<RealType><<<nBlocks, s.blockSize>>>(
                d_x.data(), d_y.data(), d_z.data(),
                d_nodeOwnership.data(), s.d_node_to_dof.data(),
                s.d_isBdryDof.data(),
                s.d_uTarget.data(), s.d_vTarget.data(), s.d_wTarget.data(),
                s.nodeCount, s.xmin, s.xmax, s.ymin, s.ymax, s.zmin, s.zmax,
                s.bboxEps, s.Uinf);
        }
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
    {
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        enforceBcMatrixKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_isBdryDof.data(), s.d_rowPtr.data(), s.d_colInd.data(),
            s.d_diagPtr.data(), s.d_valuesVel.data(), s.numOwnedDofs);
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
    if (s.bcKind == BCK::Channel)
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

    // Wrap both into SparseMatrix for the solvers.
    wrapIntoSparseMatrix<RealType>(s.Avel, s.numOwnedDofs, s.numTotalDofs, s.nnz,
                                   s.d_rowPtr.data(), s.d_colInd.data(), s.d_valuesVel.data());
    wrapIntoSparseMatrix<RealType>(s.Apre, s.numOwnedDofs, s.numTotalDofs, s.nnz,
                                   s.d_rowPtr.data(), s.d_colInd.data(), s.d_valuesPre.data());
    pt.lap("SparseMatrix wrap (vel+pre)");

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
                                       size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    // Default to 0; the halo exchange will refresh ghosts after this kernel.
    RealType uu = 0, vv = 0, ww = 0;
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
    initialConditionKernel<RealType><<<nBlocks, s.blockSize>>>(
        s.d_isBdryDof.data(), s.d_node_to_dof.data(), d_nodeOwnership.data(),
        s.d_uTarget.data(), s.d_vTarget.data(), s.d_wTarget.data(),
        s.d_u.data(), s.d_v.data(), s.d_w.data(),
        s.nodeCount);
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
                      typename NSStepper<KeyType, RealType>::Matrix& A)
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
        mars::fem::HyprePCGSolver<RealType, int, cstone::GpuTag>
            hypreSolver(MPI_COMM_WORLD, s.maxIter, s.tolerance,
                        mars::fem::HyprePCGSolver<RealType, int, cstone::GpuTag>::BOOMERAMG);
        hypreSolver.setVerbose(false);
        converged = hypreSolver.solve(
            A, b_rhs, xVec,
            static_cast<int>(s.globalRowStart), static_cast<int>(s.globalRowEnd),
            0, static_cast<int>(s.numInteriorGlobal),
            s.d_localToGlobalDof);
        // Hypre doesn't expose iteration count through the wrapper here -- use
        // -1 as a sentinel so the per-step log shows "hypre" instead of a count.
        iters = converged ? -1 : -2;
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

    // Step a: g = D^T phi, un-normalized per-node 3-vector.
    zeroVec(gxAcc); zeroVec(gyAcc); zeroVec(gzAcc);
    if (eBlocks > 0)
    {
        applyDivTransposePerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            c0, c1, c2, c3, c4, c5, c6, c7, phi.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            gxAcc.data(), gyAcc.data(), gzAcc.data(), startElem, numLocal);
        cudaDeviceSynchronize();
    }
    s.domain.reverseExchangeNodeHaloAdd(gxAcc);
    s.domain.reverseExchangeNodeHaloAdd(gyAcc);
    s.domain.reverseExchangeNodeHaloAdd(gzAcc);

    // Step b: g <- M^{-1} g  (in place; gather is same-index, safe).
    normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        gxAcc.data(), gyAcc.data(), gzAcc.data(), s.d_mass.data(),
        s.d_node_to_dof.data(), d_nodeOwnership.data(),
        gxAcc.data(), gyAcc.data(), gzAcc.data(), s.nodeCount);
    cudaDeviceSynchronize();
    s.domain.exchangeNodeHalo(gxAcc);
    s.domain.exchangeNodeHalo(gyAcc);
    s.domain.exchangeNodeHalo(gzAcc);

    // Step c: out = D g, un-normalized per-node divergence accumulator.
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

    // Pin row: out[pin] := phi[pin] kills the constant-mode null space.
    if (s.pressurePinDof >= 0)
    {
        const int pinDof      = s.pressurePinDof;
        const int* dofPtr     = s.d_node_to_dof.data();
        const uint8_t* ownPtr = d_nodeOwnership.data();
        const RealType* phPtr = phi.data();
        RealType* outPtr      = outAcc.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [pinDof, dofPtr, ownPtr, phPtr, outPtr] __device__ (size_t i) {
                if (ownPtr[i] == 1 && dofPtr[i] == pinDof) outPtr[i] = phPtr[i];
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
    rhsNode[i] = -coef * divAccNode[i];
}

// Unpreconditioned matrix-free CG for (D M^{-1} D^T) phi = b. Returns iters
// at convergence, -2 on max-iter. A production version adds Jacobi/AMG.
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

    // Fresh start: phi = 0. The pin row is phi[pin]=0, so the RHS slot also
    // gets zeroed (consistent with the SpMV's pin rule out[pin]=phi[pin]).
    thrust::fill(thrust::device_pointer_cast(phi_node.data()),
                 thrust::device_pointer_cast(phi_node.data() + s.nodeCount), RealType(0));
    if (s.pressurePinDof >= 0)
    {
        const int pinDof      = s.pressurePinDof;
        const int* dofPtr     = s.d_node_to_dof.data();
        const uint8_t* ownPtr = d_nodeOwnership.data();
        RealType* bPtr        = b_node.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [pinDof, dofPtr, ownPtr, bPtr] __device__ (size_t i) {
                if (ownPtr[i] == 1 && dofPtr[i] == pinDof) bPtr[i] = RealType(0);
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
    for (int it = 0; it < s.maxIter; ++it)
    {
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
        if (std::sqrt(rho_new) < absTol) { iters = it + 1; break; }
        RealType beta = rho_new / rho_old;
        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            p.data(), r.data(), p.data(), beta,
            s.d_node_to_dof.data(), d_nodeOwnership.data(), s.nodeCount);
        cudaDeviceSynchronize();
        rho_old = rho_new;
    }

    s.domain.exchangeNodeHalo(phi_node);
    return iters;
}

template<typename KeyType, typename RealType>
void runNsStep(NSStepper<KeyType, RealType>& s, RealType dt, RealType nu, RealType rho)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    auto c0 = std::get<0>(d_conn).data(); auto c1 = std::get<1>(d_conn).data();
    auto c2 = std::get<2>(d_conn).data(); auto c3 = std::get<3>(d_conn).data();
    auto c4 = std::get<4>(d_conn).data(); auto c5 = std::get<5>(d_conn).data();
    auto c6 = std::get<6>(d_conn).data(); auto c7 = std::get<7>(d_conn).data();

    const size_t startElem = s.domain.startIndex();
    const size_t numLocal  = s.domain.localElementCount();
    const int nodeBlocks  = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    const int eBlocks     = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
    const RealType invRho = RealType(1) / rho;
    const RealType invDt  = RealType(1) / dt;

    // -------------------------------------------------------------------------
    // Step 1: PREDICTOR. Compute (u^n . grad) q^n via the upwind flux scatter
    // for each component, then q* = q^n + dt * (advection / V) - dt * grad p^n / rho.
    //
    // Single grad(p^n) call (three components in one scatter) used for all
    // three predictors. Ghosts of u^n, v^n, w^n, p^n must be in sync because
    // face donors read them at ghost corners.
    // -------------------------------------------------------------------------
    s.domain.exchangeNodeHalo(s.d_u);
    s.domain.exchangeNodeHalo(s.d_v);
    s.domain.exchangeNodeHalo(s.d_w);
    s.domain.exchangeNodeHalo(s.d_p);

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
        normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
            s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_gradPx.data(), s.d_gradPy.data(), s.d_gradPz.data(),
            s.nodeCount);
        cudaDeviceSynchronize();
    }

    // Per-component advection scatter + predictor apply. Each scatter allocates
    // a per-node accumulator, runs the kernel, reverse-halos, then the
    // predictor reads (q^n, accumulator, gradPnq) to produce q*.
    auto runPredictor = [&] (cstone::DeviceVector<RealType>& qN,
                              cstone::DeviceVector<RealType>& qStar,
                              cstone::DeviceVector<RealType>& gradPnq,
                              cstone::DeviceVector<RealType>& qTarget)
    {
        cstone::DeviceVector<RealType> d_dqdtNode(s.nodeCount, RealType(0));
        if (eBlocks > 0)
        {
            explicitAdvectionFluxScatterPerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7,
                s.d_u.data(), s.d_v.data(), s.d_w.data(),
                qN.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                d_dqdtNode.data(), startElem, numLocal);
            cudaDeviceSynchronize();
        }
        s.domain.reverseExchangeNodeHaloAdd(d_dqdtNode);

        applyPredictorPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            qStar.data(), qN.data(), d_dqdtNode.data(),
            gradPnq.data(), s.d_mass.data(), qTarget.data(),
            s.d_isBdryDof.data(), s.d_node_to_dof.data(),
            d_nodeOwnership.data(),
            dt, invRho, s.nodeCount);
        cudaDeviceSynchronize();
    };

    runPredictor(s.d_u, s.d_uStar, s.d_gradPx, s.d_uTarget);
    runPredictor(s.d_v, s.d_vStar, s.d_gradPy, s.d_vTarget);
    runPredictor(s.d_w, s.d_wStar, s.d_gradPz, s.d_wTarget);

    // Sync ghosts of q* so the implicit RHS / CG warm-start read correct ghosts.
    s.domain.exchangeNodeHalo(s.d_uStar);
    s.domain.exchangeNodeHalo(s.d_vStar);
    s.domain.exchangeNodeHalo(s.d_wStar);

    // -------------------------------------------------------------------------
    // Step 2: IMPLICIT DIFFUSION. (M/dt + nu K) q** = (M/dt) q*  for each
    // component. The matrix is FROZEN (assembled once at setup); per step we
    // rebuild only the RHS. Velocity Dirichlet lift-off is applied to the RHS
    // for each component using its own qTarget array.
    // -------------------------------------------------------------------------
    cstone::DeviceVector<RealType> b(s.numOwnedDofs);
    cstone::DeviceVector<RealType> xVec(s.numTotalDofs);

    auto runImplicit = [&] (cstone::DeviceVector<RealType>& qStar,
                             cstone::DeviceVector<RealType>& qStarStar,
                             cstone::DeviceVector<RealType>& qTarget) -> int
    {
        // RHS = M/dt * qStar; then enforce Dirichlet (rhs[bdry] = qTarget[node]).
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

        return solveOneComponent<KeyType, RealType>(s, b, xVec, qStarStar, s.Avel);
    };

    s.lastUIters = runImplicit(s.d_uStar, s.d_uStarStar, s.d_uTarget);
    s.lastVIters = runImplicit(s.d_vStar, s.d_vStarStar, s.d_vTarget);
    s.lastWIters = runImplicit(s.d_wStar, s.d_wStarStar, s.d_wTarget);

    // Sync ghosts of q** for the divergence scatter (face donors read q** at
    // ghost corners on rank boundaries).
    s.domain.exchangeNodeHalo(s.d_uStarStar);
    s.domain.exchangeNodeHalo(s.d_vStarStar);
    s.domain.exchangeNodeHalo(s.d_wStarStar);

    // -------------------------------------------------------------------------
    // Step 3: PRESSURE POISSON. K phi = (rho/dt) div(u**). The
    // un-normalized divAccNode is exactly the integrated source f_i*V_i, so
    // the lumped RHS is rhs[dof] = (rho/dt) * divAccNode[node(dof)] -- no
    // separate normalize-then-multiply-by-V chain. Enforce the pin (single
    // Dirichlet DOF at the chosen corner) on the RHS each step.
    // -------------------------------------------------------------------------
    cstone::DeviceVector<RealType> d_divAccNode(s.nodeCount, RealType(0));
    if (eBlocks > 0)
    {
        computeDivergencePerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            c0, c1, c2, c3, c4, c5, c6, c7,
            s.d_uStarStar.data(), s.d_vStarStar.data(), s.d_wStarStar.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            d_divAccNode.data(), startElem, numLocal);
        cudaDeviceSynchronize();
    }
    s.domain.reverseExchangeNodeHaloAdd(d_divAccNode);

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
        // Fresh warm-start: phi is per-step correction, not cumulative.
        cudaMemset(xVec.data(), 0, s.numTotalDofs * sizeof(RealType));
        s.lastPressureIters = solveOneComponent<KeyType, RealType>(s, b, xVec, s.d_phi, s.Apre);
    }
    else  // DDT: matrix-free (D M^{-1} D^T) phi = -(rho/dt) D u**
    {
        cstone::DeviceVector<RealType> d_bNode(s.nodeCount, RealType(0));
        RealType coef = rho * invDt;
        buildPressureRhsDDTKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_divAccNode.data(), s.d_node_to_dof.data(),
            d_nodeOwnership.data(), coef, d_bNode.data(), s.nodeCount);
        cudaDeviceSynchronize();
        s.lastPressureIters = solvePressureDDT<KeyType, RealType>(s, d_bNode, s.d_phi);
    }

    // Sync ghosts of phi -- needed by both the gradient (face donors) and the
    // pressure update on owned nodes that read ghost slots only for VTU.
    s.domain.exchangeNodeHalo(s.d_phi);

    // -------------------------------------------------------------------------
    // Step 4: CORRECTOR + pressure update.
    //   q^{n+1} = q** - (dt/rho) (grad phi)_q  on interior
    //   q^{n+1} = qTarget                       on boundary (already correct)
    //   p^{n+1} = p^n + phi
    // grad(phi) computed once via the same B.4 kernel, applied per-component.
    // -------------------------------------------------------------------------
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
        normalizeGradientPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
            s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_gradPhix.data(), s.d_gradPhiy.data(), s.d_gradPhiz.data(),
            s.nodeCount);
        cudaDeviceSynchronize();
    }

    auto runCorrector = [&] (cstone::DeviceVector<RealType>& qOut,
                              cstone::DeviceVector<RealType>& qStarStar,
                              cstone::DeviceVector<RealType>& gradPhiq,
                              cstone::DeviceVector<RealType>& qTarget)
    {
        applyCorrectorPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            qOut.data(), qStarStar.data(), gradPhiq.data(), qTarget.data(),
            s.d_isBdryDof.data(), s.d_node_to_dof.data(),
            d_nodeOwnership.data(), dt, invRho, s.nodeCount);
        cudaDeviceSynchronize();
    };

    runCorrector(s.d_u, s.d_uStarStar, s.d_gradPhix, s.d_uTarget);
    runCorrector(s.d_v, s.d_vStarStar, s.d_gradPhiy, s.d_vTarget);
    runCorrector(s.d_w, s.d_wStarStar, s.d_gradPhiz, s.d_wTarget);

    // Pressure update + ghost sync (next step's predictor needs ghost p^{n+1}).
    updatePressureKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        s.d_p.data(), s.d_phi.data(),
        s.d_node_to_dof.data(), d_nodeOwnership.data(),
        s.nodeCount);
    cudaDeviceSynchronize();

    s.domain.exchangeNodeHalo(s.d_u);
    s.domain.exchangeNodeHalo(s.d_v);
    s.domain.exchangeNodeHalo(s.d_w);
    s.domain.exchangeNodeHalo(s.d_p);

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
        normalizeDivergencePerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_divAccPost.data(), s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            d_divNorm.data(), s.nodeCount);
        cudaDeviceSynchronize();

        // Take max |div| over interior owned DOFs (skip boundary -- the boundary
        // divergence picks up the non-zero through-flux of u_lid and is not
        // expected to be zero pointwise).
        const uint8_t* ownPtr = d_nodeOwnership.data();
        const int* dofPtr     = s.d_node_to_dof.data();
        const uint8_t* bnd    = s.d_isBdryDof.data();
        const RealType* dPtr  = d_divNorm.data();
        RealType localMax = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [ownPtr, dofPtr, bnd, dPtr] __device__ (size_t i) -> RealType {
                if (ownPtr[i] != 1) return RealType(0);
                int dof = dofPtr[i];
                if (dof < 0)            return RealType(0);
                if (bnd[dof])           return RealType(0);
                return fabs(dPtr[i]);
            }, RealType(0), thrust::maximum<RealType>());

        auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
        RealType globalMax = 0;
        MPI_Allreduce(&localMax, &globalMax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
        s.lastDivMax = globalMax;
    }
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

// =============================================================================
// main: parse args, init AMR/mesh, run setup, time loop, optionally write VTU.
// =============================================================================

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    std::string meshFile;
    CvfemKernelVariant kernelVariant = CvfemKernelVariant::Tensor;
    int blockSize    = 256;
    int bucketSize   = 64;
    int maxIter      = 1000;
    double tolerance = 1e-10;
    double rho       = 1.0;
    double nu        = 0.01;
    double dt        = 0.01;
    int numSteps     = 500;
    int vtuEvery     = 10;
    std::string vtuPrefix;
    SolverKind solverKind = SolverKind::CG;
    PressureSolveKind pressureSolve = PressureSolveKind::K;
    std::string bcKind = "cavity";
    double lidU = 1.0;
    bool useLegacyGradient = true;    // SCS gradient is the stable default; --experimental-divT swaps

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg.find("--mesh=") == 0)              meshFile = arg.substr(7);
        else if (arg.find("--kernel=") == 0)
        {
            std::string v = arg.substr(9);
            if      (v == "tensor")        kernelVariant = CvfemKernelVariant::Tensor;
            else if (v == "original")      kernelVariant = CvfemKernelVariant::Original;
            else if (v == "optimized")     kernelVariant = CvfemKernelVariant::Optimized;
            else if (v == "shmem")         kernelVariant = CvfemKernelVariant::Shmem;
            else if (v == "tensor_colored")kernelVariant = CvfemKernelVariant::TensorColored;
            else if (v == "tensor_aos")    kernelVariant = CvfemKernelVariant::TensorAoS;
            else if (v == "wmma_tensor")   kernelVariant = CvfemKernelVariant::WmmaTensor;
        }
        else if (arg.find("--block-size=") == 0)   blockSize  = std::stoi(arg.substr(13));
        else if (arg.find("--bucket-size=") == 0)  bucketSize = std::stoi(arg.substr(14));
        else if (arg.find("--max-iter=") == 0)     maxIter    = std::stoi(arg.substr(11));
        else if (arg.find("--tol=") == 0)          tolerance  = std::stod(arg.substr(6));
        else if (arg.find("--rho=") == 0)          rho        = std::stod(arg.substr(6));
        else if (arg.find("--nu=") == 0)           nu         = std::stod(arg.substr(5));
        else if (arg.find("--dt=") == 0)           dt         = std::stod(arg.substr(5));
        else if (arg.find("--num-steps=") == 0)    numSteps   = std::stoi(arg.substr(12));
        else if (arg.find("--vtu-every=") == 0)    vtuEvery   = std::stoi(arg.substr(12));
        else if (arg.find("--vtu-output=") == 0)   vtuPrefix  = arg.substr(13);
        else if (arg.find("--lid-u=") == 0)        lidU       = std::stod(arg.substr(8));
        else if (arg.find("--bc=") == 0)           bcKind     = arg.substr(5);
        else if (arg.find("--solver=") == 0)
        {
            std::string v = arg.substr(9);
            if      (v == "cg")    solverKind = SolverKind::CG;
            else if (v == "hypre") solverKind = SolverKind::Hypre;
            else
            {
                if (rank == 0) std::cerr << "Error: --solver must be cg or hypre, got '" << v << "'\n";
                MPI_Finalize();
                return 1;
            }
        }
        else if (arg == "--experimental-divT") useLegacyGradient = false;
        else if (arg.find("--pressure-solve=") == 0)
        {
            std::string v = arg.substr(17);
            if      (v == "K")   pressureSolve = PressureSolveKind::K;
            else if (v == "DDT") pressureSolve = PressureSolveKind::DDT;
            else
            {
                if (rank == 0) std::cerr << "Error: --pressure-solve must be K or DDT, got '" << v << "'\n";
                MPI_Finalize();
                return 1;
            }
        }
        else if (arg[0] != '-' && meshFile.empty()) meshFile = arg;
    }

    if (meshFile.empty())
    {
        if (rank == 0)
        {
            std::cout << "Usage: " << argv[0] << " --mesh=FILE [options]\n\n"
                      << "Chorin projection NS driver (Phase C).\n\n"
                      << "Options:\n"
                      << "  --mesh=FILE           Mesh file (.mesh or .exo) [REQUIRED]\n"
                      << "  --bc=cavity|channel   BC config (default: cavity)\n"
                      << "                          cavity:  lid u=lidU on top face, no-slip elsewhere\n"
                      << "                          channel: Dirichlet u=lidU inflow on x=xmin,\n"
                      << "                                   no-slip walls, natural outflow on x=xmax\n"
                      << "  --rho=VALUE           Density (default 1.0)\n"
                      << "  --nu=VALUE            Kinematic viscosity (default 0.01)\n"
                      << "  --dt=VALUE            Timestep (default 0.01)\n"
                      << "  --num-steps=N         Number of time steps (default 500)\n"
                      << "  --lid-u=VALUE         Lid velocity for cavity (default 1.0)\n"
                      << "  --kernel=VARIANT      tensor (default), optimized, shmem, wmma_tensor, ...\n"
                      << "  --solver=KIND         cg | hypre (default cg)\n"
                      << "  --pressure-solve=KIND K | DDT (default K, matrix-free Galerkin if DDT)\n"
                      << "  --vtu-output=PREFIX   Write per-step VTU/PVTU/PVD frames\n"
                      << "  --vtu-every=N         Write a frame every N steps (default 10)\n"
                      << "  --max-iter=N          Solver max iterations (default 1000)\n"
                      << "  --tol=VALUE           Solver tolerance (default 1e-10)\n"
                      << "  --bucket-size=N       Cornerstone bucket size (default 64)\n"
                      << "  --block-size=N        CUDA block size (default 256)\n";
        }
        MPI_Finalize();
        return 1;
    }

    if (bcKind != "cavity" && bcKind != "channel")
    {
        if (rank == 0)
            std::cerr << "Error: --bc must be 'cavity' or 'channel', got '" << bcKind << "'\n";
        MPI_Finalize();
        return 1;
    }
    if (bcKind == "channel" && pressureSolve == PressureSolveKind::DDT)
    {
        if (rank == 0)
            std::cerr << "Error: --pressure-solve=DDT is not supported with --bc=channel.\n"
                      << "       The DDT path's pin enforcement assumes a single corner-pin DOF;\n"
                      << "       channel mode uses a full outflow-face Dirichlet BC instead.\n"
                      << "       Use --pressure-solve=K (the default) with --bc=channel.\n";
        MPI_Finalize();
        return 1;
    }

    using KeyType  = uint64_t;
    using RealType = double;

    if (rank == 0)
    {
        const double Re = lidU * 1.0 / nu;  // L = 1 for the unit cube
        std::cout << "\n========================================\n";
        std::cout << "MARS NS Chorin projection driver (Phase C)\n";
        std::cout << "========================================\n";
        if (bcKind == "channel")
            std::cout << "BC:        channel flow, Uinf=" << lidU << " inflow on x=xmin, no-slip walls, natural outflow on x=xmax\n";
        else
            std::cout << "BC:        lid-driven cavity, u=" << lidU << " on top face\n";
        std::cout << "rho       = " << rho << "\n";
        std::cout << "nu        = " << nu  << "\n";
        std::cout << "dt        = " << dt  << "\n";
        std::cout << "Re (~L*U/nu, L=1) = " << Re << "\n";
        std::cout << "numSteps  = " << numSteps << " (T_final = " << numSteps * dt << ")\n";
        std::cout << "Mesh:    " << meshFile << "\n";
        std::cout << "Kernel:  " << CvfemHexAssembler<KeyType, RealType>::variantName(kernelVariant) << "\n";
        std::cout << "Solver:  " << (solverKind == SolverKind::CG ? "cg" : "hypre (PCG+BoomerAMG)") << "\n";
        std::cout << "Pressure solve: " << (pressureSolve == PressureSolveKind::K
                                            ? "K (CVFEM stiffness)"
                                            : "DDT (matrix-free D M^-1 D^T)") << "\n";
        std::cout << "MPI ranks: " << numRanks << "\n";
        std::cout << "========================================\n\n";
    }

    // AmrManager only for the mesh load + lazy halo/adj/coord build (maxLevels=0
    // means the mesh is frozen for this driver -- AMR-on-NS is out of scope for
    // Phase C).
    AmrManager<HexTag, KeyType, RealType>::Config amrConfig;
    amrConfig.maxLevels  = 0;
    amrConfig.blockSize  = blockSize;
    amrConfig.bucketSize = bucketSize;

    AmrManager<HexTag, KeyType, RealType> amr(amrConfig);
    amr.initialize(meshFile, rank, numRanks);
    auto initT = amr.initTimings();

    if (rank == 0)
    {
        std::cout << "Initial mesh:\n";
        std::cout << "  Elements: " << amr.domain().getElementCount() << "\n";
        std::cout << "  Nodes:    " << amr.domain().getNodeCount() << "\n";
        std::cout << "  Init breakdown (ms):\n";
        std::cout << "    domain (file read + sync):    " << std::fixed << initT.domainSyncTimeMs << "\n";
        std::cout << "    halo + node topology:         " << initT.haloTopoTimeMs   << "\n";
        std::cout << "    adjacency CSR (e2n + n2e):    " << initT.adjacencyTimeMs  << "\n";
        std::cout << "    coord cache (SFC decode):     " << initT.coordCacheTimeMs << "\n";
        std::cout << "    AMR octree state:             " << initT.octreeTimeMs     << "\n";
        std::cout << "    TOTAL:                        " << initT.totalMs          << "\n\n";
    }

    NSStepper<KeyType, RealType> s{amr.domain(), solverKind, blockSize, maxIter,
                                   RealType(tolerance), rank, numRanks};
    s.lidU = RealType(lidU);
    s.Uinf = RealType(lidU);   // reuse the --lid-u CLI flag as channel inflow speed
    s.bcKind = (bcKind == "channel") ? NSStepper<KeyType, RealType>::BCKind::Channel
                                      : NSStepper<KeyType, RealType>::BCKind::Cavity;
    s.useLegacyGradient = useLegacyGradient;
    s.pressureSolve = pressureSolve;
    setupNSStepper<KeyType, RealType>(s, RealType(nu), RealType(dt), kernelVariant);

    // IC: u = lid on top face, 0 elsewhere; p = 0.
    applyInitialCondition<KeyType, RealType>(s);

    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriterU;
    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriterP;
    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriterUmag;
    if (!vtuPrefix.empty())
    {
        // One writer per field: ParaView reads them as separate PVD timelines.
        vtuWriterU    = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix + "_u");
        vtuWriterP    = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix + "_p");
        vtuWriterUmag = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix + "_umag");
        if (rank == 0)
            std::cout << "VTU output enabled: " << vtuPrefix << "_{u,p,umag}_step*.pvtu\n\n";
    }

    auto writeVtuFrame = [&] (int step, double t)
    {
        if (!vtuWriterU) return;
        vtuWriterU->writeFrame(step, t, amr.domain(), s.d_u, "u");
        vtuWriterP->writeFrame(step, t, amr.domain(), s.d_p, "p");
        // |u| magnitude as a third field for visualization.
        cstone::DeviceVector<RealType> d_umag(s.nodeCount, RealType(0));
        const RealType* uPtr = s.d_u.data();
        const RealType* vPtr = s.d_v.data();
        const RealType* wPtr = s.d_w.data();
        RealType* mPtr       = d_umag.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [uPtr, vPtr, wPtr, mPtr] __device__ (size_t i) {
                RealType u = uPtr[i], v = vPtr[i], w = wPtr[i];
                mPtr[i] = sqrt(u * u + v * v + w * w);
            });
        cudaDeviceSynchronize();
        vtuWriterUmag->writeFrame(step, t, amr.domain(), d_umag, "umag");
    };

    // Step 0 report: norms of the IC (|u| is the lid contribution only).
    {
        RealType nu0 = computeWeightedL2Norm<KeyType, RealType>(s, s.d_u);
        RealType nv0 = computeWeightedL2Norm<KeyType, RealType>(s, s.d_v);
        RealType nw0 = computeWeightedL2Norm<KeyType, RealType>(s, s.d_w);
        RealType np0 = computeWeightedL2Norm<KeyType, RealType>(s, s.d_p);
        if (rank == 0)
        {
            std::cout << "Step " << std::setw(4) << 0 << ": "
                      << "t=" << std::fixed << std::setprecision(4) << 0.0
                      << " |u|=" << std::scientific << std::setprecision(3) << nu0
                      << " |v|=" << nv0
                      << " |w|=" << nw0
                      << " |p|=" << np0
                      << " div_max=0.000e+00 cg_iter_p=N/A\n"
                      << std::defaultfloat;
        }
        if (vtuWriterU) writeVtuFrame(0, 0.0);
    }

    auto totalStart = std::chrono::high_resolution_clock::now();

    for (int step = 1; step <= numSteps; ++step)
    {
        runNsStep<KeyType, RealType>(s, RealType(dt), RealType(nu), RealType(rho));

        double t = step * dt;

        RealType nuN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_u);
        RealType nvN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_v);
        RealType nwN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_w);
        RealType npN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_p);

        if (rank == 0)
        {
            std::cout << "Step " << std::setw(4) << step << ": "
                      << "t=" << std::fixed << std::setprecision(4) << t
                      << " |u|=" << std::scientific << std::setprecision(3) << nuN
                      << " |v|=" << nvN
                      << " |w|=" << nwN
                      << " |p|=" << npN
                      << " div_max=" << s.lastDivMax
                      << " cg_iter_p=";
            if (s.lastPressureIters == -1)       std::cout << "hypre";
            else if (s.lastPressureIters == -2)  std::cout << "FAIL";
            else                                  std::cout << s.lastPressureIters;
            // cg_iter_uvw: report the max across the three velocity solves so
            // a single number summarizes convergence (Hypre: -1; FAIL: -2).
            int vMax = std::max({s.lastUIters, s.lastVIters, s.lastWIters});
            std::cout << " cg_iter_uvw=";
            if (vMax == -1)      std::cout << "hypre";
            else if (vMax == -2) std::cout << "FAIL";
            else                  std::cout << vMax;
            std::cout << "\n" << std::defaultfloat;
        }

        if (vtuWriterU && (step % vtuEvery == 0 || step == numSteps))
        {
            writeVtuFrame(step, t);
            if (rank == 0)
                std::cout << "  VTU: wrote " << vtuPrefix << "_{u,p,umag}_step"
                          << std::setw(4) << std::setfill('0') << step << ".pvtu\n"
                          << std::setfill(' ');
        }
    }

    auto totalEnd = std::chrono::high_resolution_clock::now();
    float totalMs = std::chrono::duration<float, std::milli>(totalEnd - totalStart).count();

    if (rank == 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "NS projection run complete\n";
        std::cout << "  Time steps: " << numSteps << "\n";
        std::cout << "  Total wall: " << std::fixed << totalMs << " ms ("
                  << (totalMs / std::max(numSteps, 1)) << " ms/step)\n";
        std::cout << "========================================\n";
    }

    MPI_Finalize();
    return 0;
}
