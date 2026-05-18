// Advection-diffusion validator on a fixed cube mesh. Operator-split scheme:
//   1) explicit upwind advection: T* = T^n - dt * div(v T^n)_upwind
//   2) implicit backward-Euler diffusion: (M + dt*nu*K) T^{n+1} = M T*
// Matrix (M/dt + nu K) is assembled and BC-enforced ONCE in setup since the
// velocity doesn't enter it. Step 1 reduces to T*=T^n when omega=0 (recovers
// B.1 bit-for-bit); step 2 reduces to T^{n+1}=T* when nu=0 (pure-advection
// mass conservation test). Velocity is solid-body rotation about the z-axis:
//     v(x,y,z) = ( -Omega*(y - 0.5),  Omega*(x - 0.5),  0 ).
// Analytical reference (nu>0, optionally back-rotated by -Omega*t):
//     T_ana(x, t) = (sigma^2/(sigma^2+2 nu t))^(3/2)
//                   * exp(-|R(-Omega*t)(x-c)+c-x_c|^2 / (2 (sigma^2+2 nu t))).

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
#include <thrust/inner_product.h>
#include <thrust/transform.h>
#include <thrust/system/cuda/execution_policy.h>

#include <memory>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <limits>

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

// PhaseTimer: sync on lap so the wall-clock reflects all-rank work.
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
            std::cout << "    " << std::left << std::setw(28) << n
                      << std::right << std::fixed << std::setprecision(2)
                      << std::setw(10) << ms << " ms ("
                      << std::setw(5) << std::setprecision(1)
                      << (100.0f * ms / std::max(total, 1e-3f)) << "%)\n";
        }
        std::cout << "    " << std::left << std::setw(28) << "TOTAL"
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(10) << total << " ms\n";
    }
};

// Mark owned boundary DOFs (true on any of the 6 cube faces).
template<typename RealType>
__global__ void applyBCKernel(const RealType* nodeX,
                              const RealType* nodeY,
                              const RealType* nodeZ,
                              const uint8_t* ownership,
                              const int* nodeToDof,
                              uint8_t* isBoundaryDof,
                              size_t numNodes,
                              RealType xmin,
                              RealType xmax,
                              RealType ymin,
                              RealType ymax,
                              RealType zmin,
                              RealType zmax,
                              RealType eps)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;

    if (ownership[i] == 1)
    {
        int dof = nodeToDof[i];
        if (dof >= 0)
        {
            RealType x = nodeX[i], y = nodeY[i], z = nodeZ[i];
            bool onBdry = (fabs(x - xmin) < eps || fabs(x - xmax) < eps || fabs(y - ymin) < eps ||
                           fabs(y - ymax) < eps || fabs(z - zmin) < eps || fabs(z - zmax) < eps);
            isBoundaryDof[dof] = onBdry ? 1 : 0;
        }
    }
}

// Zero boundary rows, set diag=1. Applied once; matrix is frozen after that.
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

// Zero RHS on boundary DOFs (T=0 matches T_ana ~ 0 far from x_c).
template<typename RealType>
__global__ void enforceBcRhsKernel(const uint8_t* isBoundaryDof,
                                   RealType* rhs,
                                   int numDofs)
{
    int dof = blockIdx.x * blockDim.x + threadIdx.x;
    if (dof >= numDofs) return;
    if (isBoundaryDof[dof]) rhs[dof] = RealType(0);
}

// Legacy lumped-mass path: per-DOF scatter with ownership filter on targets.
// Under-counts at corner-only-neighbor partitions where the cstone halo misses
// element fans around a boundary node. Kept for --legacy-scatter comparisons.
template<typename KeyType, typename RealType>
__global__ void computeLumpedMassKernel(const KeyType* conn0,
                                        const KeyType* conn1,
                                        const KeyType* conn2,
                                        const KeyType* conn3,
                                        const KeyType* conn4,
                                        const KeyType* conn5,
                                        const KeyType* conn6,
                                        const KeyType* conn7,
                                        const RealType* nodeX,
                                        const RealType* nodeY,
                                        const RealType* nodeZ,
                                        const int* nodeToDof,
                                        const uint8_t* ownership,
                                        RealType* mass,
                                        size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType n[8] = {conn0[e], conn1[e], conn2[e], conn3[e],
                    conn4[e], conn5[e], conn6[e], conn7[e]};

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
    {
        if (ownership[n[i]] == 1)
        {
            int dof = nodeToDof[n[i]];
            if (dof >= 0) atomicAdd(&mass[dof], contrib);
        }
    }
}

// Per-node lumped-mass variant matching the per-node flux scatter: iterate
// LOCAL elements only, scatter V/8 to all 8 corner node slots with no
// ownership filter. Caller follows with reverseExchangeNodeHaloAdd + a gather
// to the per-DOF layout to close the corner-only-neighbor gap.
template<typename KeyType, typename RealType>
__global__ void computeLumpedMassPerNodeKernel(const KeyType* conn0,
                                               const KeyType* conn1,
                                               const KeyType* conn2,
                                               const KeyType* conn3,
                                               const KeyType* conn4,
                                               const KeyType* conn5,
                                               const KeyType* conn6,
                                               const KeyType* conn7,
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

    KeyType n[8] = {conn0[e], conn1[e], conn2[e], conn3[e],
                    conn4[e], conn5[e], conn6[e], conn7[e]};

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

// Add M[i]/dt onto the diagonal of (nu K). Run after assembleFull, before BC.
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

// Initial condition: 3D Gaussian sampled at every (owned + ghost) node.
template<typename RealType>
__global__ void gaussianInitialConditionKernel(const RealType* nodeX,
                                               const RealType* nodeY,
                                               const RealType* nodeZ,
                                               RealType* nodeT,
                                               size_t numNodes,
                                               RealType xc,
                                               RealType yc,
                                               RealType zc,
                                               RealType sigma)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    RealType dx = nodeX[i] - xc;
    RealType dy = nodeY[i] - yc;
    RealType dz = nodeZ[i] - zc;
    RealType r2 = dx * dx + dy * dy + dz * dz;
    nodeT[i] = exp(-r2 / (RealType(2) * sigma * sigma));
}

// b[dof] = M[dof]/dt * T[node(dof)]. Called with T=T* for the implicit RHS.
template<typename RealType>
__global__ void buildRhsFromTimeNKernel(const RealType* nodeT,
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
    rhs[dof] = mass[dof] * invDt * nodeT[i];
}

// Per-owned-DOF (T-T_ana)^2 * V contribution to ||err||_L2. theta=-omega*t
// back-rotates the sampler so a rigidly rotating Gaussian lines up.
template<typename RealType>
__global__ void l2ErrorKernel(const RealType* nodeX,
                              const RealType* nodeY,
                              const RealType* nodeZ,
                              const RealType* nodeT,
                              const int* nodeToDof,
                              const uint8_t* ownership,
                              const RealType* mass,
                              RealType* errSqPerDof,
                              size_t numNodes,
                              RealType xc,
                              RealType yc,
                              RealType zc,
                              RealType sigma,
                              RealType nu,
                              RealType t,
                              RealType theta,
                              RealType rotCx,
                              RealType rotCy)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;

    RealType c = cos(theta), s = sin(theta);
    RealType px = nodeX[i] - rotCx, py = nodeY[i] - rotCy;
    RealType xRot = rotCx + c * px - s * py;
    RealType yRot = rotCy + s * px + c * py;

    RealType s2 = sigma * sigma;
    RealType denom = s2 + RealType(2) * nu * t;
    RealType amp = pow(s2 / denom, RealType(1.5));
    RealType dx = xRot - xc;
    RealType dy = yRot - yc;
    RealType dz = nodeZ[i] - zc;
    RealType r2 = dx * dx + dy * dy + dz * dz;
    RealType Tana = amp * exp(-r2 / (RealType(2) * denom));

    RealType diff = nodeT[i] - Tana;
    errSqPerDof[dof] = diff * diff * mass[dof];
}

// Solid-body z-rotation around (cx,cy,*). Steady, evaluated once at setup.
template<typename RealType>
__global__ void setRotationVelocityKernel(const RealType* nodeX,
                                          const RealType* nodeY,
                                          RealType* vx,
                                          RealType* vy,
                                          RealType* vz,
                                          size_t numNodes,
                                          RealType omega,
                                          RealType cx,
                                          RealType cy)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    vx[i] = -omega * (nodeY[i] - cy);
    vy[i] =  omega * (nodeX[i] - cx);
    vz[i] = RealType(0);
}

// Per-node upwind flux scatter. The per-node accumulator (sized nodeCount, no
// ownership filter) over LOCAL elements only ensures each global face is
// scattered exactly once with equal-and-opposite contributions. The follow-up
// reverseExchangeNodeHaloAdd sums ghost-slot contributions back to owners,
// closing the corner-only-neighbor gap the legacy per-DOF + ownership-filter
// path silently drops.
template<typename KeyType, typename RealType>
__global__ void explicitAdvectionFluxScatterPerNodeKernel(const KeyType* conn0,
                                                          const KeyType* conn1,
                                                          const KeyType* conn2,
                                                          const KeyType* conn3,
                                                          const KeyType* conn4,
                                                          const KeyType* conn5,
                                                          const KeyType* conn6,
                                                          const KeyType* conn7,
                                                          const RealType* vx,
                                                          const RealType* vy,
                                                          const RealType* vz,
                                                          const RealType* T,
                                                          const RealType* areaVecX,
                                                          const RealType* areaVecY,
                                                          const RealType* areaVecZ,
                                                          RealType* dTdtNode,
                                                          size_t startElem,
                                                          size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;

    KeyType n[8] = {conn0[e], conn1[e], conn2[e], conn3[e],
                    conn4[e], conn5[e], conn6[e], conn7[e]};

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

        RealType T_up = (mdot > RealType(0)) ? T[iL] : T[iR];
        RealType flux = mdot * T_up;

        atomicAdd(&dTdtNode[iL], -flux);
        atomicAdd(&dTdtNode[iR], +flux);
    }
}

// Legacy per-DOF + ownership-filter scatter. Kept for --legacy-scatter.
template<typename KeyType, typename RealType>
__global__ void explicitAdvectionFluxScatterKernel(const KeyType* conn0,
                                                   const KeyType* conn1,
                                                   const KeyType* conn2,
                                                   const KeyType* conn3,
                                                   const KeyType* conn4,
                                                   const KeyType* conn5,
                                                   const KeyType* conn6,
                                                   const KeyType* conn7,
                                                   const RealType* vx,
                                                   const RealType* vy,
                                                   const RealType* vz,
                                                   const RealType* T,
                                                   const RealType* areaVecX,
                                                   const RealType* areaVecY,
                                                   const RealType* areaVecZ,
                                                   const int* nodeToDof,
                                                   const uint8_t* ownership,
                                                   RealType* dTdt,
                                                   size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType n[8] = {conn0[e], conn1[e], conn2[e], conn3[e],
                    conn4[e], conn5[e], conn6[e], conn7[e]};

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

        RealType T_up = (mdot > RealType(0)) ? T[iL] : T[iR];
        RealType flux = mdot * T_up;

        if (ownership[iL] == 1)
        {
            int dofL = nodeToDof[iL];
            if (dofL >= 0) atomicAdd(&dTdt[dofL], -flux);
        }
        if (ownership[iR] == 1)
        {
            int dofR = nodeToDof[iR];
            if (dofR >= 0) atomicAdd(&dTdt[dofR], +flux);
        }
    }
}

// T += dt * dTdtNode[i] / M[dof(i)] on owned nodes. Ghosts stale; caller
// re-exchanges before the implicit RHS build.
template<typename RealType>
__global__ void applyExplicitAdvectionPerNodeKernel(RealType* T,
                                                    const RealType* dTdtNode,
                                                    const RealType* lumpedMass,
                                                    const int* nodeToDof,
                                                    const uint8_t* ownership,
                                                    RealType dt,
                                                    size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    T[i] += dt * dTdtNode[i] / lumpedMass[dof];
}

// Legacy apply paired with the per-DOF dTdt accumulator.
template<typename RealType>
__global__ void applyExplicitAdvectionKernel(RealType* T,
                                             const RealType* dTdt,
                                             const RealType* lumpedMass,
                                             const int* nodeToDof,
                                             const uint8_t* ownership,
                                             RealType dt,
                                             size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    T[i] += dt * dTdt[dof] / lumpedMass[dof];
}

enum class SolverKind { CG, Hypre };

// Driver state persisted across steps. Built once; loop reuses everything.
template<typename KeyType, typename RealType>
struct AdvDiffStepper
{
    using DomainT = ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>;

    DomainT& domain;
    SolverKind solverKind;
    int blockSize;
    int maxIter;
    RealType tolerance;
    int rank;
    bool useLegacy = false;  // --legacy-scatter: per-DOF + ownership-filter scatter

    size_t nodeCount = 0;
    size_t elementCount = 0;
    int numOwnedDofs = 0;
    int numTotalDofs = 0;
    int nnz = 0;

    cstone::DeviceVector<int> d_node_to_dof;
    cstone::DeviceVector<int> d_rowPtr;
    cstone::DeviceVector<int> d_colInd;
    cstone::DeviceVector<int> d_diagPtr;
    cstone::DeviceVector<RealType> d_values;
    cstone::DeviceVector<RealType> d_mass;
    cstone::DeviceVector<uint8_t> d_isBdryDof;

    // Per-node steady velocity field, read by the per-step advection scatter.
    cstone::DeviceVector<RealType> d_vx;
    cstone::DeviceVector<RealType> d_vy;
    cstone::DeviceVector<RealType> d_vz;

    // Geometry, precomputed once; only the explicit scatter reads it.
    cstone::DeviceVector<RealType> d_areaVec_x;
    cstone::DeviceVector<RealType> d_areaVec_y;
    cstone::DeviceVector<RealType> d_areaVec_z;

    // Hypre-specific persistent state (only filled when solverKind == Hypre).
    int64_t globalRowStart = 0;
    int64_t globalRowEnd = 0;
    int64_t numInteriorGlobal = 0;
    std::vector<int64_t> localToGlobalDof;

    // Owned-row CSR wrapper consumed by the CG / Hypre solvers.
    using Matrix = SparseMatrix<int, RealType, cstone::GpuTag>;
    Matrix A;

    // Cached bbox for BC marking (constant on a fixed mesh).
    RealType xmin = 0, xmax = 1, ymin = 0, ymax = 1, zmin = 0, zmax = 1;
    RealType bboxEps = 0;
};

// One-time setup: implicit matrix (M/dt+nu K), lumped mass, BC mask, velocity,
// SCS area vectors, persistent CSR/DOF state. Frozen thereafter.
template<typename KeyType, typename RealType>
void setupAdvDiffStepper(AdvDiffStepper<KeyType, RealType>& s,
                         RealType nu,
                         RealType dt,
                         RealType omega,
                         RealType rotCx,
                         RealType rotCy,
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

    // Owned nodes -> [0, numOwnedDofs); ghosts -> [numOwnedDofs, nodeCount).
    s.d_node_to_dof.resize(s.nodeCount);
    s.numOwnedDofs = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(), s.d_node_to_dof.data(), s.nodeCount);
    s.numTotalDofs = static_cast<int>(s.nodeCount);
    pt.lap("DOF mapping");

    // Full 27-NNZ pattern, matching assembleFull's 8x8 element stencil.
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

    s.d_values.resize(s.nnz);
    thrust::fill(thrust::device_pointer_cast(s.d_values.data()),
                 thrust::device_pointer_cast(s.d_values.data() + s.nnz),
                 RealType(0));
    pt.lap("sparsity build");

    // Assembler config and CSR wrapper for the matrix.
    using MatrixT = CSRMatrix<RealType>;
    MatrixT* d_matrix;
    cudaMalloc(&d_matrix, sizeof(MatrixT));
    MatrixT h_matrix{s.d_rowPtr.data(), s.d_colInd.data(), s.d_values.data(), s.d_diagPtr.data(),
                     s.numTotalDofs, s.nnz};
    cudaMemcpy(d_matrix, &h_matrix, sizeof(MatrixT), cudaMemcpyHostToDevice);

    // gamma=nu, all advection inputs zero -> assembler produces pure nu*K.
    // RHS is built per step outside the assembler.
    cstone::DeviceVector<RealType> d_gamma(s.nodeCount, RealType(nu));
    cstone::DeviceVector<RealType> d_phi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_beta(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_grad_phi_x(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_grad_phi_y(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_grad_phi_z(s.nodeCount, RealType(0));

    s.d_areaVec_x.resize(s.elementCount * 12);
    s.d_areaVec_y.resize(s.elementCount * 12);
    s.d_areaVec_z.resize(s.elementCount * 12);
    precomputeAreaVectorsGpu<KeyType, RealType>(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), s.elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data());

    // Per-node rotation velocity (steady).
    s.d_vx.resize(s.nodeCount);
    s.d_vy.resize(s.nodeCount);
    s.d_vz.resize(s.nodeCount);
    {
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        setRotationVelocityKernel<RealType><<<nBlocks, s.blockSize>>>(
            d_x.data(), d_y.data(),
            s.d_vx.data(), s.d_vy.data(), s.d_vz.data(),
            s.nodeCount, omega, rotCx, rotCy);
        cudaDeviceSynchronize();
    }
    pt.lap("velocity + area vectors");

    cstone::DeviceVector<RealType> d_mdot_zero(s.elementCount * 12, RealType(0));
    cstone::DeviceVector<RealType> d_rhs_unused(s.numTotalDofs, RealType(0));

    typename CvfemHexAssembler<KeyType, RealType>::Config config;
    config.blockSize = s.blockSize;
    config.variant   = kernelVariant;

    CvfemHexAssembler<KeyType, RealType>::assembleFull(
        std::get<0>(d_conn).data(), std::get<1>(d_conn).data(), std::get<2>(d_conn).data(),
        std::get<3>(d_conn).data(), std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
        std::get<6>(d_conn).data(), std::get<7>(d_conn).data(), s.elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        d_gamma.data(), d_phi.data(), d_beta.data(),
        d_grad_phi_x.data(), d_grad_phi_y.data(), d_grad_phi_z.data(),
        d_mdot_zero.data(), s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
        s.d_node_to_dof.data(), d_nodeOwnership.data(), d_matrix, d_rhs_unused.data(), config);
    cudaDeviceSynchronize();
    pt.lap("assembly (nu*K)");

    // Lumped mass: default uses per-node + reverse-halo so M matches the per-
    // node flux scatter at corner-only-neighbor partitions; --legacy-scatter
    // keeps the per-DOF + ownership-filter path for side-by-side comparison.
    s.d_mass.resize(s.numOwnedDofs);
    thrust::fill(thrust::device_pointer_cast(s.d_mass.data()),
                 thrust::device_pointer_cast(s.d_mass.data() + s.numOwnedDofs),
                 RealType(0));
    if (s.useLegacy)
    {
        int eBlocks = (s.elementCount + s.blockSize - 1) / s.blockSize;
        computeLumpedMassKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            d_x.data(), d_y.data(), d_z.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_mass.data(), s.elementCount);
        cudaDeviceSynchronize();
    }
    else
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

    // (M/dt + nu K) = (nu K already in values) + (M/dt on diag).
    {
        int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
        addLumpedMassDiagonalKernel<RealType><<<dofBlocks, s.blockSize>>>(
            s.d_mass.data(), s.d_diagPtr.data(), RealType(1) / dt,
            s.d_values.data(), s.numOwnedDofs);
        cudaDeviceSynchronize();
    }
    pt.lap("add M/dt to diagonal");

    // Global bounding box (used for BC marking only).
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

    s.d_isBdryDof.resize(s.numOwnedDofs);
    thrust::fill(thrust::device_pointer_cast(s.d_isBdryDof.data()),
                 thrust::device_pointer_cast(s.d_isBdryDof.data() + s.numOwnedDofs),
                 uint8_t(0));
    int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    applyBCKernel<<<nodeBlocks, s.blockSize>>>(
        d_x.data(), d_y.data(), d_z.data(), d_nodeOwnership.data(),
        s.d_node_to_dof.data(), s.d_isBdryDof.data(), s.nodeCount,
        s.xmin, s.xmax, s.ymin, s.ymax, s.zmin, s.zmax, s.bboxEps);

    int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
    enforceBcMatrixKernel<RealType><<<dofBlocks, s.blockSize>>>(
        s.d_isBdryDof.data(), s.d_rowPtr.data(), s.d_colInd.data(), s.d_diagPtr.data(),
        s.d_values.data(), s.numOwnedDofs);
    cudaDeviceSynchronize();
    pt.lap("BC enforcement");

    // CSR wrapper for the solver. Frozen after this point.
    s.A.allocate(s.numOwnedDofs, s.numTotalDofs, s.nnz);
    cudaMemcpy(s.A.rowOffsetsPtr(), s.d_rowPtr.data(), (s.numOwnedDofs + 1) * sizeof(int),
               cudaMemcpyDeviceToDevice);
    cudaMemcpy(s.A.colIndicesPtr(), s.d_colInd.data(), s.nnz * sizeof(int),
               cudaMemcpyDeviceToDevice);
    cudaMemcpy(s.A.valuesPtr(),     s.d_values.data(), s.nnz * sizeof(RealType),
               cudaMemcpyDeviceToDevice);
    cudaFree(d_matrix);
    pt.lap("SparseMatrix wrap");

#ifdef MARS_ENABLE_HYPRE
    if (s.solverKind == SolverKind::Hypre)
    {
        // Build contiguous global owned-DOF numbering: rank r owns
        // [globalRowStart, globalRowEnd). Ghost local columns get filled by
        // a halo exchange of the per-node global DOF index.
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

        std::vector<RealType> h_nodeGlobalDof(s.nodeCount);
        thrust::copy(thrust::device_pointer_cast(d_nodeGlobalDof.data()),
                     thrust::device_pointer_cast(d_nodeGlobalDof.data() + s.nodeCount),
                     h_nodeGlobalDof.begin());
        std::vector<int> h_nodeToDof(s.nodeCount);
        thrust::copy(thrust::device_pointer_cast(s.d_node_to_dof.data()),
                     thrust::device_pointer_cast(s.d_node_to_dof.data() + s.nodeCount),
                     h_nodeToDof.begin());

        s.localToGlobalDof.assign(s.numTotalDofs, int64_t(-1));
        for (size_t i = 0; i < s.nodeCount; ++i)
        {
            int dof = h_nodeToDof[i];
            if (dof >= 0) s.localToGlobalDof[dof] = static_cast<int64_t>(h_nodeGlobalDof[i]);
        }
        pt.lap("Hypre global DOF map");
    }
#endif

    pt.report(s.rank, "setup");
}

// One operator-split step. d_nodeT comes in as T^n, leaves as T^{n+1}.
template<typename KeyType, typename RealType>
void runAdvDiffTimeStep(AdvDiffStepper<KeyType, RealType>& s,
                        RealType dt,
                        cstone::DeviceVector<RealType>& d_nodeT,
                        int& outConverged)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();

    // Step 1: explicit upwind advection (T^n -> T*). Face donors read T on
    // ghosts, so we sync ghosts on entry; the apply leaves ghosts stale and we
    // re-exchange before building the implicit RHS.
    s.domain.exchangeNodeHalo(d_nodeT);

    if (s.useLegacy)
    {
        cstone::DeviceVector<RealType> d_dTdt(s.numOwnedDofs, RealType(0));
        int eBlocks = (s.elementCount + s.blockSize - 1) / s.blockSize;
        explicitAdvectionFluxScatterKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            s.d_vx.data(), s.d_vy.data(), s.d_vz.data(),
            d_nodeT.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            d_dTdt.data(), s.elementCount);
        cudaDeviceSynchronize();

        int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        applyExplicitAdvectionKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_nodeT.data(), d_dTdt.data(), s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            dt, s.nodeCount);
        cudaDeviceSynchronize();
    }
    else
    {
        // Per-node accumulator + reverse halo closes the corner-only-neighbor
        // gap that the legacy per-DOF + ownership-filter scatter silently drops.
        cstone::DeviceVector<RealType> d_dTdtNode(s.nodeCount, RealType(0));

        size_t startElem = s.domain.startIndex();
        size_t numLocal  = s.domain.localElementCount();
        if (numLocal > 0)
        {
            int eBlocks = int((numLocal + s.blockSize - 1) / s.blockSize);
            explicitAdvectionFluxScatterPerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
                std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
                std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
                std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
                s.d_vx.data(), s.d_vy.data(), s.d_vz.data(),
                d_nodeT.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                d_dTdtNode.data(), startElem, numLocal);
            cudaDeviceSynchronize();
        }

        s.domain.reverseExchangeNodeHaloAdd(d_dTdtNode);

        int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        applyExplicitAdvectionPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            d_nodeT.data(), d_dTdtNode.data(), s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            dt, s.nodeCount);
        cudaDeviceSynchronize();
    }

    // Sync ghosts of T*: implicit RHS / CG warm-start need them.
    s.domain.exchangeNodeHalo(d_nodeT);

    // -------- Step 2: implicit backward-Euler diffusion (T* -> T^{n+1}) -----
    cstone::DeviceVector<RealType> b(s.numOwnedDofs, RealType(0));
    cstone::DeviceVector<RealType> x(s.numTotalDofs);

    int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    buildRhsFromTimeNKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        d_nodeT.data(), s.d_node_to_dof.data(), d_nodeOwnership.data(),
        s.d_mass.data(), RealType(1) / dt, b.data(), s.nodeCount);

    int dofBlocks = (s.numOwnedDofs + s.blockSize - 1) / s.blockSize;
    enforceBcRhsKernel<RealType><<<dofBlocks, s.blockSize>>>(
        s.d_isBdryDof.data(), b.data(), s.numOwnedDofs);
    cudaDeviceSynchronize();

    cudaMemset(x.data(), 0, s.numTotalDofs * sizeof(RealType));

    // Warm-start CG with T* scattered into DOF order.
    thrust::for_each(thrust::device,
                      thrust::counting_iterator<size_t>(0),
                      thrust::counting_iterator<size_t>(s.nodeCount),
                      [nodeToDof = s.d_node_to_dof.data(),
                       Tn        = d_nodeT.data(),
                       xOut      = x.data()] __device__(size_t i) {
                          int dof = nodeToDof[i];
                          if (dof >= 0) xOut[dof] = Tn[i];
                      });

    int nRanks = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

    bool converged = false;
    if (s.solverKind == SolverKind::CG)
    {
        ConjugateGradientSolver<RealType, int, cstone::GpuTag> solver(s.maxIter, s.tolerance);
        solver.setVerbose(false);
        solver.setOwnedSize(s.numOwnedDofs);
        if (nRanks > 1)
        {
            const int* dofMapPtr = s.d_node_to_dof.data();
            auto& dom = s.domain;
            solver.setHaloExchangeCallback(
                [&dom, dofMapPtr](cstone::DeviceVector<RealType>& p)
                {
                    dom.exchangeNodeHalo(p, dofMapPtr);
                });
        }
        converged = solver.solve(s.A, b, x);
    }
    else
    {
#ifdef MARS_ENABLE_HYPRE
        using IdxType = int64_t;
        mars::fem::HyprePCGSolver<RealType, int, cstone::GpuTag>
            hypreSolver(MPI_COMM_WORLD, s.maxIter, s.tolerance,
                        mars::fem::HyprePCGSolver<RealType, int, cstone::GpuTag>::BOOMERAMG);
        hypreSolver.setVerbose(false);
        converged = hypreSolver.template solve<IdxType>(
            s.A, b, x,
            static_cast<int>(s.globalRowStart), static_cast<int>(s.globalRowEnd),
            0, static_cast<int>(s.numInteriorGlobal),
            s.localToGlobalDof);
#else
        if (s.rank == 0)
            std::cerr << "ERROR: --solver=hypre requested but MARS was built without MARS_ENABLE_HYPRE\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    cudaDeviceSynchronize();

    // x (DOF-indexed, owned rows + ghost cols) -> d_nodeT (node-indexed).
    thrust::for_each(thrust::device,
                      thrust::counting_iterator<size_t>(0),
                      thrust::counting_iterator<size_t>(s.nodeCount),
                      [nodeToDof = s.d_node_to_dof.data(),
                       sol       = x.data(),
                       nodeSol   = d_nodeT.data()] __device__(size_t i)
                      {
                          int dof = nodeToDof[i];
                          if (dof >= 0) nodeSol[i] = sol[dof];
                      });

    s.domain.exchangeNodeHalo(d_nodeT);

    outConverged = converged ? 1 : 0;
}

// Global L2 error vs the (back-rotated) analytical Gaussian.
template<typename KeyType, typename RealType>
RealType computeL2ErrorVsGaussian(AdvDiffStepper<KeyType, RealType>& s,
                                  const cstone::DeviceVector<RealType>& d_nodeT,
                                  RealType xc,
                                  RealType yc,
                                  RealType zc,
                                  RealType sigma,
                                  RealType nu,
                                  RealType t,
                                  RealType theta,
                                  RealType rotCx,
                                  RealType rotCy)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();

    cstone::DeviceVector<RealType> d_errSq(s.numOwnedDofs, RealType(0));

    int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    l2ErrorKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        d_x.data(), d_y.data(), d_z.data(),
        d_nodeT.data(), s.d_node_to_dof.data(), d_nodeOwnership.data(),
        s.d_mass.data(), d_errSq.data(), s.nodeCount,
        xc, yc, zc, sigma, nu, t, theta, rotCx, rotCy);
    cudaDeviceSynchronize();

    auto ep = thrust::device_pointer_cast(d_errSq.data());
    RealType localSum = thrust::reduce(thrust::device, ep, ep + s.numOwnedDofs, RealType(0));
    RealType globalSum = 0;
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&localSum, &globalSum, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    return std::sqrt(globalSum);
}

// Global integrated mass (sum_i T_i V_i over owned DOFs). Conserved by the
// explicit advection step; decays only via the Dirichlet boundary leak.
template<typename KeyType, typename RealType>
RealType computeIntegratedMass(AdvDiffStepper<KeyType, RealType>& s,
                               const cstone::DeviceVector<RealType>& d_nodeT)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    cstone::DeviceVector<RealType> d_mv(s.numOwnedDofs, RealType(0));
    thrust::for_each(thrust::device,
                      thrust::counting_iterator<size_t>(0),
                      thrust::counting_iterator<size_t>(s.nodeCount),
                      [nodeToDof = s.d_node_to_dof.data(),
                       ownPtr    = d_nodeOwnership.data(),
                       mass      = s.d_mass.data(),
                       Tn        = d_nodeT.data(),
                       out       = d_mv.data()] __device__(size_t i)
                      {
                          if (ownPtr[i] != 1) return;
                          int dof = nodeToDof[i];
                          if (dof < 0) return;
                          out[dof] = Tn[i] * mass[dof];
                      });
    cudaDeviceSynchronize();

    auto sp = thrust::device_pointer_cast(d_mv.data());
    RealType localSum = thrust::reduce(thrust::device, sp, sp + s.numOwnedDofs, RealType(0));
    RealType globalSum = 0;
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&localSum, &globalSum, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    return globalSum;
}

// Global ||T||_2 over owned DOFs (volume-weighted by lumped mass).
template<typename KeyType, typename RealType>
RealType computeWeightedL2Norm(AdvDiffStepper<KeyType, RealType>& s,
                               const cstone::DeviceVector<RealType>& d_nodeT)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    cstone::DeviceVector<RealType> d_sq(s.numOwnedDofs, RealType(0));
    thrust::for_each(thrust::device,
                      thrust::counting_iterator<size_t>(0),
                      thrust::counting_iterator<size_t>(s.nodeCount),
                      [nodeToDof = s.d_node_to_dof.data(),
                       ownPtr    = d_nodeOwnership.data(),
                       mass      = s.d_mass.data(),
                       Tn        = d_nodeT.data(),
                       out       = d_sq.data()] __device__(size_t i)
                      {
                          if (ownPtr[i] != 1) return;
                          int dof = nodeToDof[i];
                          if (dof < 0) return;
                          out[dof] = Tn[i] * Tn[i] * mass[dof];
                      });
    cudaDeviceSynchronize();

    auto sp = thrust::device_pointer_cast(d_sq.data());
    RealType localSum = thrust::reduce(thrust::device, sp, sp + s.numOwnedDofs, RealType(0));
    RealType globalSum = 0;
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&localSum, &globalSum, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    return std::sqrt(globalSum);
}

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
    double nu        = 0.001;
    double dt        = 0.01;
    int numSteps     = 50;
    double sigma     = 0.05;
    double xc        = 0.5;
    double yc        = 0.5;
    double zc        = 0.5;
    double omega     = 0.0;  // solid-body rotation rate around the z-axis
    int vtuEvery     = 10;
    std::string vtuPrefix;
    SolverKind solverKind = SolverKind::CG;
    bool useLegacy = false;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg.find("--mesh=") == 0)
            meshFile = arg.substr(7);
        else if (arg.find("--kernel=") == 0)
        {
            std::string v = arg.substr(9);
            if (v == "tensor") kernelVariant = CvfemKernelVariant::Tensor;
            else if (v == "original") kernelVariant = CvfemKernelVariant::Original;
            else if (v == "optimized") kernelVariant = CvfemKernelVariant::Optimized;
            else if (v == "shmem") kernelVariant = CvfemKernelVariant::Shmem;
            else if (v == "tensor_colored") kernelVariant = CvfemKernelVariant::TensorColored;
            else if (v == "tensor_aos") kernelVariant = CvfemKernelVariant::TensorAoS;
            else if (v == "wmma_tensor") kernelVariant = CvfemKernelVariant::WmmaTensor;
        }
        else if (arg.find("--block-size=") == 0)
            blockSize = std::stoi(arg.substr(13));
        else if (arg.find("--bucket-size=") == 0)
            bucketSize = std::stoi(arg.substr(14));
        else if (arg.find("--max-iter=") == 0)
            maxIter = std::stoi(arg.substr(11));
        else if (arg.find("--tol=") == 0)
            tolerance = std::stod(arg.substr(6));
        else if (arg.find("--nu=") == 0)
            nu = std::stod(arg.substr(5));
        else if (arg.find("--dt=") == 0)
            dt = std::stod(arg.substr(5));
        else if (arg.find("--num-steps=") == 0)
            numSteps = std::stoi(arg.substr(12));
        else if (arg.find("--sigma=") == 0)
            sigma = std::stod(arg.substr(8));
        else if (arg.find("--center=") == 0)
        {
            std::string v = arg.substr(9);
            size_t c1 = v.find(','), c2 = (c1 == std::string::npos ? std::string::npos : v.find(',', c1 + 1));
            if (c1 == std::string::npos || c2 == std::string::npos)
            {
                if (rank == 0)
                    std::cerr << "Error: --center=x,y,z requires 3 comma-separated values\n";
                MPI_Finalize();
                return 1;
            }
            xc = std::stod(v.substr(0, c1));
            yc = std::stod(v.substr(c1 + 1, c2 - c1 - 1));
            zc = std::stod(v.substr(c2 + 1));
        }
        else if (arg.find("--omega=") == 0)
            omega = std::stod(arg.substr(8));
        else if (arg.find("--vtu-every=") == 0)
            vtuEvery = std::stoi(arg.substr(12));
        else if (arg.find("--vtu-output=") == 0)
            vtuPrefix = arg.substr(13);
        else if (arg.find("--solver=") == 0)
        {
            std::string v = arg.substr(9);
            if (v == "cg") solverKind = SolverKind::CG;
            else if (v == "hypre") solverKind = SolverKind::Hypre;
            else
            {
                if (rank == 0) std::cerr << "Error: --solver must be cg or hypre, got '" << v << "'\n";
                MPI_Finalize();
                return 1;
            }
        }
        else if (arg == "--legacy-scatter")
            useLegacy = true;
        else if (arg[0] != '-' && meshFile.empty())
            meshFile = arg;
    }

    if (meshFile.empty())
    {
        if (rank == 0)
        {
            std::cout << "Usage: " << argv[0] << " [options]\n\n";
            std::cout << "Options:\n";
            std::cout << "  --mesh=FILE          Mesh file (.mesh or .exo) [REQUIRED]\n";
            std::cout << "  --kernel=VARIANT     tensor (default), original, optimized, shmem, wmma_tensor\n";
            std::cout << "  --nu=VALUE           Diffusion coefficient (default: 0.001)\n";
            std::cout << "  --dt=VALUE           Timestep (default: 0.01)\n";
            std::cout << "  --num-steps=N        Number of time steps (default: 50)\n";
            std::cout << "  --sigma=VALUE        Gaussian width sigma (default: 0.05)\n";
            std::cout << "  --center=X,Y,Z       Gaussian center (default: 0.5,0.5,0.5)\n";
            std::cout << "  --omega=VALUE        Solid-body z-axis rotation rate (default: 0 = pure diffusion)\n";
            std::cout << "  --vtu-output=PREFIX  Write VTU/PVTU/PVD frames (default: off)\n";
            std::cout << "  --vtu-every=N        Write a frame every N steps (default: 10)\n";
            std::cout << "  --max-iter=N         Solver max iterations (default: 1000)\n";
            std::cout << "  --tol=VALUE          Solver tolerance (default: 1e-10)\n";
            std::cout << "  --bucket-size=N      Cornerstone bucket size (default: 64)\n";
            std::cout << "  --block-size=N       CUDA block size (default: 256)\n";
            std::cout << "  --solver=KIND        cg | hypre (default: cg; hypre requires MARS_ENABLE_HYPRE=ON)\n";
            std::cout << "  --legacy-scatter     Use the per-DOF + ownership-filter flux/mass scatter (bug demo)\n";
        }
        MPI_Finalize();
        return 1;
    }

    using KeyType  = uint64_t;
    using RealType = double;

    // Rotation axis is fixed to the cube centerline (0.5, 0.5). The Gaussian
    // center is independent: pick something off-axis (--center=0.7,0.5,0.5)
    // to actually see orbital motion.
    const double rotCx = 0.5;
    const double rotCy = 0.5;

    if (rank == 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "MARS Advection-Diffusion Validator (CVFEM, operator-split: explicit upwind + implicit BE)\n";
        std::cout << "========================================\n";
        std::cout << "Problem: dT/dt + v.grad(T) = nu * Laplacian(T), T = 0 on cube faces\n";
        std::cout << "         v = (-omega*(y-0.5), omega*(x-0.5), 0)\n";
        std::cout << "IC:      Gaussian, sigma=" << sigma
                  << ", center=(" << xc << "," << yc << "," << zc << ")\n";
        std::cout << "nu       = " << nu << "\n";
        std::cout << "omega    = " << omega
                  << (omega == 0.0 ? "  (pure diffusion)" : "  (rotation+diffusion)") << "\n";
        std::cout << "dt       = " << dt << "\n";
        std::cout << "numSteps = " << numSteps << " (T_final = " << numSteps * dt << ")\n";
        std::cout << "Mesh:    " << meshFile << "\n";
        std::cout << "Kernel:  " << CvfemHexAssembler<KeyType, RealType>::variantName(kernelVariant) << "\n";
        std::cout << "Solver:  " << (solverKind == SolverKind::CG ? "cg" : "hypre (PCG+BoomerAMG)") << "\n";
        std::cout << "MPI ranks: " << numRanks << "\n";
        std::cout << "Scatter: " << (useLegacy ? "LEGACY (per-DOF + ownership filter)"
                                                : "per-node + reverseExchangeNodeHaloAdd") << "\n";
        std::cout << "========================================\n\n";
    }

    // Reuse AmrManager only for the mesh load + lazy halo+adjacency+coord build.
    // maxLevels=0 means we never call adaptMesh; the mesh is frozen.
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

    AdvDiffStepper<KeyType, RealType> s{amr.domain(), solverKind, blockSize, maxIter,
                                        RealType(tolerance), rank, useLegacy};
    setupAdvDiffStepper<KeyType, RealType>(s, RealType(nu), RealType(dt),
                                           RealType(omega), RealType(rotCx), RealType(rotCy),
                                           kernelVariant);

    // Initial condition (per node, including ghosts so VTU output is clean
    // even before the first halo exchange).
    cstone::DeviceVector<RealType> d_nodeT(s.nodeCount, RealType(0));
    {
        const auto& d_x = amr.domain().getNodeX();
        const auto& d_y = amr.domain().getNodeY();
        const auto& d_z = amr.domain().getNodeZ();
        int nodeBlocks = (s.nodeCount + blockSize - 1) / blockSize;
        gaussianInitialConditionKernel<RealType><<<nodeBlocks, blockSize>>>(
            d_x.data(), d_y.data(), d_z.data(),
            d_nodeT.data(), s.nodeCount,
            RealType(xc), RealType(yc), RealType(zc), RealType(sigma));
        cudaDeviceSynchronize();
    }

    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriter;
    if (!vtuPrefix.empty())
    {
        vtuWriter = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix);
        if (rank == 0)
            std::cout << "VTU output enabled: " << vtuPrefix << "_step*.pvtu, "
                      << vtuPrefix << ".pvd\n\n";
    }

    // Step 0: report IC norm and error (should be ~0 since T_ana(t=0) = T_IC).
    {
        RealType tnorm = computeWeightedL2Norm<KeyType, RealType>(s, d_nodeT);
        RealType tmass = computeIntegratedMass<KeyType, RealType>(s, d_nodeT);
        RealType terr  = computeL2ErrorVsGaussian<KeyType, RealType>(
            s, d_nodeT, RealType(xc), RealType(yc), RealType(zc),
            RealType(sigma), RealType(nu), RealType(0),
            RealType(0), RealType(rotCx), RealType(rotCy));
        if (rank == 0)
        {
            std::cout << "Step " << std::setw(4) << 0 << ":  t=" << std::fixed << std::setprecision(4) << 0.0
                      << "  |T|=" << std::scientific << std::setprecision(3) << tnorm
                      << "  intT=" << tmass
                      << "  |err|=" << terr
                      << "  solver=initial\n"
                      << std::defaultfloat;
        }
        if (vtuWriter) vtuWriter->writeFrame(0, 0.0, amr.domain(), d_nodeT, "T");
    }

    auto totalStart = std::chrono::high_resolution_clock::now();

    for (int step = 1; step <= numSteps; ++step)
    {
        int converged = 0;
        runAdvDiffTimeStep<KeyType, RealType>(s, RealType(dt), d_nodeT, converged);

        double t = step * dt;
        // Back-rotate by -omega*t so a rigid rotation lands the reference on
        // top of the numerical solution. theta=0 for omega=0 (pure diffusion).
        double theta = -omega * t;

        RealType tnorm = computeWeightedL2Norm<KeyType, RealType>(s, d_nodeT);
        RealType tmass = computeIntegratedMass<KeyType, RealType>(s, d_nodeT);
        RealType terr  = computeL2ErrorVsGaussian<KeyType, RealType>(
            s, d_nodeT, RealType(xc), RealType(yc), RealType(zc),
            RealType(sigma), RealType(nu), RealType(t),
            RealType(theta), RealType(rotCx), RealType(rotCy));

        if (rank == 0)
        {
            std::cout << "Step " << std::setw(4) << step << ":  t=" << std::fixed << std::setprecision(4) << t
                      << "  |T|=" << std::scientific << std::setprecision(3) << tnorm
                      << "  intT=" << tmass
                      << "  |err|=" << terr
                      << "  solver=" << (converged ? "ok" : "FAIL")
                      << "\n" << std::defaultfloat;
        }

        if (vtuWriter && (step % vtuEvery == 0 || step == numSteps))
        {
            vtuWriter->writeFrame(step, t, amr.domain(), d_nodeT, "T");
            if (rank == 0)
                std::cout << "  VTU: wrote " << vtuPrefix << "_step"
                          << std::setw(4) << std::setfill('0') << step << ".pvtu\n"
                          << std::setfill(' ');
        }
    }

    auto totalEnd = std::chrono::high_resolution_clock::now();
    float totalMs = std::chrono::duration<float, std::milli>(totalEnd - totalStart).count();

    if (rank == 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "Advection-diffusion run complete\n";
        std::cout << "  Time steps: " << numSteps << "\n";
        std::cout << "  Total wall: " << std::fixed << totalMs << " ms ("
                  << (totalMs / std::max(numSteps, 1)) << " ms/step)\n";
        std::cout << "========================================\n";
    }

    MPI_Finalize();
    return 0;
}
