// Discrete divergence operator validator on a fixed cube mesh.
//
// Phase B.3 of the NS-projection plan. The pressure Poisson equation in the
// projection method is
//     div(grad p) = (1/dt) * div(u*),
// so we need a div(u) operator that is *consistent with the same CVFEM
// control-volume framework and SCS face layout* as the implicit Laplacian.
// Otherwise the projection cannot drive the discrete divergence to zero.
//
// Discrete divergence at owned node i with control volume V_i = lumpedMass[i]:
//     (div u)_i = (1/V_i) * sum over SCS faces f incident to i:  v_face_f . A_f
// where the sum runs over the same 12 SCS faces per hex used in the advection
// scatter, v_face = 0.5 (v[L] + v[R]) on a face connecting corners (L, R), and
// the area vector A_f follows precomputeAreaVectorsGpu's convention (pointing
// from L to R). Each face contributes +v.A to L's net outflow and -v.A to R's
// (this is the *opposite sign* of the advection scatter, which records flux
// *leaving* L as -flux on L).
//
// Two analytical test fields:
//   --field=rotation : v = (-omega*(y-cy), omega*(x-cx), 0)  ->  div = 0 exactly
//   --field=sinx     : v = (sin(pi*x), 0, 0)                 ->  div = pi*cos(pi*x)
//
// The rotation case verifies normalization and sign cancellation; the sinx
// case verifies sign convention and O(h) convergence away from boundaries.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/amr/mars_amr.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>
#include <thrust/system/cuda/execution_policy.h>

#include <memory>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

// Same pattern as mars_amr_advdiff: sync on lap so the wall-clock includes
// every rank's work, not just rank-0's stragglers.
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

enum class TestField { Rotation, Sinx };

// Per-element axis-aligned-box volume scattered to its 8 owned corner DOFs
// (each gets V/8). Same V_i = lumpedMass[i] used to normalize the divergence
// so the validator's V_i matches whatever projection solver will eventually use.
// Legacy variant (use with --use-legacy): iterates ALL elements with ownership
// filter, suffers the same cstone-halo-gap under-count as the legacy flux
// scatter at corner-only-neighbor partitions.
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

// Per-NODE mass variant that mirrors computeDivergencePerNodeKernel: iterates
// LOCAL elements only and scatters V/8 to all 8 corner node slots (no
// ownership filter). The caller follows with reverseExchangeNodeHaloAdd to
// sum ghost-slot contributions back to each node's owner. This keeps the
// lumped mass consistent with the per-node flux scatter so that single-rank
// and multi-rank produce the same V_node for every owned node.
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
    {
        atomicAdd(&massNode[n[i]], contrib);
    }
}

// Gather owned-node mass values into a per-DOF array (size = numOwnedDofs).
// After computeLumpedMassPerNodeKernel + reverseExchangeNodeHaloAdd, massNode
// holds the full lumped mass at every OWNED node's local slot. This kernel
// transposes that into the per-DOF layout the divergence-normalize kernel
// expects, keeping the public mass-array shape unchanged downstream.
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

// Per-node test velocity. Steady, evaluated once at setup.
template<typename RealType>
__global__ void setTestVelocityKernel(const RealType* nodeX,
                                      const RealType* nodeY,
                                      const RealType* /*nodeZ*/,
                                      RealType* vx,
                                      RealType* vy,
                                      RealType* vz,
                                      size_t numNodes,
                                      int fieldKind,    // 0 = rotation, 1 = sinx
                                      RealType omega,
                                      RealType cx,
                                      RealType cy)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;

    if (fieldKind == 0)
    {
        vx[i] = -omega * (nodeY[i] - cy);
        vy[i] =  omega * (nodeX[i] - cx);
        vz[i] = RealType(0);
    }
    else
    {
        constexpr RealType PI = RealType(3.14159265358979323846);
        vx[i] = sin(PI * nodeX[i]);
        vy[i] = RealType(0);
        vz[i] = RealType(0);
    }
}

// Divergence scatter. For each SCS face f of element e separating corners
// (iL, iR) given by d_hexLRSCV with area vector A pointing from L to R:
//     v_face = 0.5 (v[iL] + v[iR])
//     flow   = v_face . A_f   (signed contribution to L's net outflow)
//     divAcc[L] += +flow      (mass leaves L)
//     divAcc[R] += -flow      (mass enters R)
//
// This is the *opposite sign* of explicitAdvectionFluxScatterKernel, which
// records flow leaving L as -flux on L (because that kernel scatters dT/dt for
// a transported scalar, not a divergence). Magnitude per face is identical.
//
// Per-node accumulator variant: writes into a nodeCount-sized buffer indexed by
// LOCAL NODE ID (not DOF), with NO ownership filter. Iterates ONLY local
// elements [startElem, startElem+numLocal) -- not halos -- so each physical
// face is scattered exactly once globally (by the rank that owns its element).
// After scatter, ghost-slot contributions on each rank are summed back to the
// owner via reverseExchangeNodeHaloAdd so that nodes touched by elements on
// neighbor ranks get the missing remote contributions.
template<typename KeyType, typename RealType>
__global__ void computeDivergencePerNodeKernel(const KeyType* conn0,
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
        RealType flow = vfx * areaVecX[off] + vfy * areaVecY[off] + vfz * areaVecZ[off];

        // No ownership filter: scatter into both endpoints' node slots. The
        // reverseExchangeNodeHaloAdd later moves ghost-slot contributions back
        // to the rightful owners.
        atomicAdd(&divAccNode[iL], +flow);
        atomicAdd(&divAccNode[iR], -flow);
    }
}

// Legacy per-rank kernel (kept for reference and side-by-side comparison via
// --use-legacy). For an SCS face whose L is owned by rank A and R by rank B,
// this assumes both ranks visit the same physical face. When the cstone halo
// does NOT include a corner-only-neighbor element owned by the other rank,
// some contributions are missed and the global face sum is non-zero. See the
// header comment on computeDivergencePerNodeKernel for the conservative fix.
template<typename KeyType, typename RealType>
__global__ void computeDivergenceKernel(const KeyType* conn0,
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
                                        const RealType* areaVecX,
                                        const RealType* areaVecY,
                                        const RealType* areaVecZ,
                                        const int* nodeToDof,
                                        const uint8_t* ownership,
                                        RealType* divAcc,
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
        RealType flow = vfx * areaVecX[off] + vfy * areaVecY[off] + vfz * areaVecZ[off];

        if (ownership[iL] == 1)
        {
            int dofL = nodeToDof[iL];
            if (dofL >= 0) atomicAdd(&divAcc[dofL], +flow);
        }
        if (ownership[iR] == 1)
        {
            int dofR = nodeToDof[iR];
            if (dofR >= 0) atomicAdd(&divAcc[dofR], -flow);
        }
    }
}

// Divide accumulated net outflow by the lumped CV volume to get (div u)_i.
// Scatters the per-DOF result back to a per-node array so the VTU writer and
// the analytical-error kernel can see it node-indexed (and halo exchange can
// fill ghosts if needed for visualization).
template<typename RealType>
__global__ void applyDivergenceNormalizeKernel(const RealType* divAcc,
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
    divU[i] = divAcc[dof] / lumpedMass[dof];
}

// Per-node-indexed normalize. divAccNode is sized nodeCount and indexed by
// local node id (matches computeDivergencePerNodeKernel + the reverse halo).
// The lumped-mass array is still per-owned-DOF, so we look up via nodeToDof.
template<typename RealType>
__global__ void applyDivergenceNormalizePerNodeKernel(const RealType* divAccNode,
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

// ============================================================================
// Diagnostic kernels (gated by --diagnose). Help pin down cross-rank scatter
// inconsistencies: which #1-#4 cause from the task description is at fault.
// ============================================================================

// Diagnostic 1: per-owned-node SCS-face-touch count. Iterates ALL elements the
// rank sees (local + halo) and matches the LEGACY scatter's predication: only
// count faces where the rank-owns the endpoint (matches the legacy kernel).
// For interior cube16 nodes a single-rank run reports 12 touches per owned
// node (4 incident edges x 3 SCS faces per edge / 2 sharing); boundary nodes
// fewer. On 4 ranks the SAME global node should still report 12 (after
// summation across ranks since each rank counts only what IT sees). If on the
// owner-rank the count drops below the single-rank reference for the same
// physical position, then the owner's halo is missing element fans (cause #1).
template<typename KeyType>
__global__ void countFacesTouchedKernel(const KeyType* conn0,
                                        const KeyType* conn1,
                                        const KeyType* conn2,
                                        const KeyType* conn3,
                                        const KeyType* conn4,
                                        const KeyType* conn5,
                                        const KeyType* conn6,
                                        const KeyType* conn7,
                                        const uint8_t* ownership,
                                        int* faceTouchCountNode,
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

        if (ownership[iL] == 1) atomicAdd(&faceTouchCountNode[iL], 1);
        if (ownership[iR] == 1) atomicAdd(&faceTouchCountNode[iR], 1);
    }
}

// Diagnostic 2: per-NODE net flow (no ownership filter on scatter). Iterates
// LOCAL elements only [startElem, startElem+numLocal). Used to compare the
// owner-side partial sum (this rank's local-element contributions only)
// against the reverse-halo-summed value. The delta is exactly what the
// reverse halo recovers from neighbor ranks' local elements that touch the
// node here but aren't in this rank's halo.
template<typename KeyType, typename RealType>
__global__ void perNodeNetFlowKernel(const KeyType* conn0,
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
        RealType flow = vfx * areaVecX[off] + vfy * areaVecY[off] + vfz * areaVecZ[off];

        atomicAdd(&divAccNode[iL], +flow);
        atomicAdd(&divAccNode[iR], -flow);
    }
}

// Diagnostic 3: dump per-face contributions for a specific target local node
// from LOCAL elements only. Each rank dumps what its local elements scatter to
// the node. For the bug case, the per-rank dumps summed across ranks at the
// same physical node should give a zero-divergence signature for an interior
// node; deviations from zero localize the cause (#2 area-vector mismatch,
// #3 connectivity-order issue, #4 velocity-ghost staleness).
template<typename KeyType, typename RealType>
__global__ void dumpFaceContribsForNodeKernel(const KeyType* conn0,
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
                                              const RealType* areaVecX,
                                              const RealType* areaVecY,
                                              const RealType* areaVecZ,
                                              const uint8_t* ownership,
                                              KeyType targetLocalNode,
                                              int*       outCount,
                                              int*       outElem,
                                              int*       outIp,
                                              int*       outRoleIsL,   // 1 if target == iL, 0 if iR
                                              int*       outOwnL,
                                              int*       outOwnR,
                                              RealType*  outFlow,
                                              RealType*  outAx,
                                              RealType*  outAy,
                                              RealType*  outAz,
                                              size_t startElem,
                                              size_t numLocal,
                                              int maxOut)
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

        bool hit = (iL == targetLocalNode) || (iR == targetLocalNode);
        if (!hit) continue;

        RealType vfx = RealType(0.5) * (vx[iL] + vx[iR]);
        RealType vfy = RealType(0.5) * (vy[iL] + vy[iR]);
        RealType vfz = RealType(0.5) * (vz[iL] + vz[iR]);

        size_t off = e * 12 + ip;
        RealType flow = vfx * areaVecX[off] + vfy * areaVecY[off] + vfz * areaVecZ[off];

        int slot = atomicAdd(outCount, 1);
        if (slot < maxOut)
        {
            outElem[slot]    = int(e);
            outIp[slot]      = ip;
            outRoleIsL[slot] = (iL == targetLocalNode) ? 1 : 0;
            outOwnL[slot]    = int(ownership[iL]);
            outOwnR[slot]    = int(ownership[iR]);
            outFlow[slot]    = flow;
            outAx[slot]      = areaVecX[off];
            outAy[slot]      = areaVecY[off];
            outAz[slot]      = areaVecZ[off];
        }
    }
}

// Search every local node (owned or ghost) for the one whose physical position
// is closest to (tx, ty, tz). Returns (local_node_id, distance). Used to pick
// the same physical node on every rank for cross-rank dumps. Distance==1e30
// signals "no local node found".
template<typename RealType>
struct ProbeLocalNodeAt
{
    RealType targetX, targetY, targetZ;
    const RealType* nodeX;
    const RealType* nodeY;
    const RealType* nodeZ;

    __device__ thrust::pair<RealType, size_t> operator()(size_t i) const
    {
        RealType dx = nodeX[i] - targetX;
        RealType dy = nodeY[i] - targetY;
        RealType dz = nodeZ[i] - targetZ;
        return {dx*dx + dy*dy + dz*dz, i};
    }
};

template<typename RealType>
struct PairMinByFirst
{
    __device__ thrust::pair<RealType, size_t> operator()(
        const thrust::pair<RealType, size_t>& a,
        const thrust::pair<RealType, size_t>& b) const
    {
        return (a.first < b.first) ? a : b;
    }
};

// Per-owned-DOF (divU_num - divU_ana)^2 * V_node for the L2 error reduction.
// fieldKind selects the analytical reference.
//
// Also writes a per-DOF boundary mask. The CVFEM control-volume divergence
// at a boundary node intentionally measures the net outflow through the CV,
// which includes the boundary face flux. For a field like rigid-body rotation
// where v.n != 0 on the cube faces, that boundary CV term is non-zero even
// though div(v) = 0 analytically -- not an operator bug, just the definition
// of the discrete CV divergence at a domain boundary. So we separate boundary
// and interior contributions in the L2 norms downstream.
//
// Hardcoded unit cube [0,1]^3 with eps=1e-3; matches the cube16 driver use.
template<typename RealType>
__global__ void computeDivErrorKernel(const RealType* nodeX,
                                      const RealType* nodeY,
                                      const RealType* nodeZ,
                                      const RealType* divU,
                                      const int* nodeToDof,
                                      const uint8_t* ownership,
                                      const RealType* mass,
                                      RealType* errSqIntPerDof,
                                      RealType* anaSqIntPerDof,
                                      RealType* numSqIntPerDof,
                                      RealType* numSqBndPerDof,
                                      uint8_t* boundaryMaskPerDof,
                                      size_t numNodes,
                                      int fieldKind)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;

    RealType ana = RealType(0);
    if (fieldKind == 1)
    {
        constexpr RealType PI = RealType(3.14159265358979323846);
        ana = PI * cos(PI * nodeX[i]);
    }

    RealType num  = divU[i];
    RealType diff = num - ana;
    RealType m    = mass[dof];

    constexpr RealType eps = RealType(1e-3);
    RealType x = nodeX[i], y = nodeY[i], z = nodeZ[i];
    bool isBnd = (fabs(x) < eps) || (fabs(x - RealType(1)) < eps) ||
                 (fabs(y) < eps) || (fabs(y - RealType(1)) < eps) ||
                 (fabs(z) < eps) || (fabs(z - RealType(1)) < eps);

    boundaryMaskPerDof[dof] = isBnd ? uint8_t(1) : uint8_t(0);

    if (isBnd)
    {
        numSqBndPerDof[dof] = num * num * m;
    }
    else
    {
        errSqIntPerDof[dof] = diff * diff * m;
        anaSqIntPerDof[dof] = ana  * ana  * m;
        numSqIntPerDof[dof] = num  * num  * m;
    }
}

// Search owned interior nodes for the one closest to a target point; return
// the local divU value plus the L1 distance to the target. Used for the
// spot-check report ("at x=0.5 div: ~0"). Run as a thrust transform_reduce
// across owned DOFs and pick the rank-local minimum, then MPI_Allreduce with
// MPI_MINLOC to grab the global winner's value.
template<typename RealType>
struct ProbeAtX
{
    RealType targetX;
    RealType targetY;
    RealType targetZ;
    const RealType* nodeX;
    const RealType* nodeY;
    const RealType* nodeZ;
    const RealType* divU;
    const int* nodeToDof;
    const uint8_t* ownership;

    __device__ thrust::pair<RealType, RealType> operator()(size_t i) const
    {
        if (ownership[i] != 1) return {RealType(1e30), RealType(0)};
        int dof = nodeToDof[i];
        if (dof < 0) return {RealType(1e30), RealType(0)};
        RealType dx = nodeX[i] - targetX;
        RealType dy = nodeY[i] - targetY;
        RealType dz = nodeZ[i] - targetZ;
        RealType d2 = dx * dx + dy * dy + dz * dz;
        return {d2, divU[i]};
    }
};

template<typename RealType>
struct MinByFirst
{
    __device__ thrust::pair<RealType, RealType> operator()(
        const thrust::pair<RealType, RealType>& a,
        const thrust::pair<RealType, RealType>& b) const
    {
        return (a.first < b.first) ? a : b;
    }
};

// Returns the divU at the owned node closest to (tx, ty, tz). All ranks
// participate; the winning rank's value is broadcast via MPI_MINLOC.
template<typename RealType>
RealType probeDivAt(const cstone::DeviceVector<RealType>& d_x,
                    const cstone::DeviceVector<RealType>& d_y,
                    const cstone::DeviceVector<RealType>& d_z,
                    const cstone::DeviceVector<RealType>& d_divU,
                    const cstone::DeviceVector<int>& d_nodeToDof,
                    const cstone::DeviceVector<uint8_t>& d_ownership,
                    size_t numNodes,
                    RealType tx, RealType ty, RealType tz)
{
    ProbeAtX<RealType> op{tx, ty, tz,
                          d_x.data(), d_y.data(), d_z.data(),
                          d_divU.data(),
                          d_nodeToDof.data(), d_ownership.data()};

    auto init = thrust::make_pair(RealType(1e30), RealType(0));
    auto best = thrust::transform_reduce(
        thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(numNodes),
        op, init, MinByFirst<RealType>());

    // MPI_MINLOC: pick the (distance, rank) pair with smallest distance, then
    // bcast the corresponding divU value from that rank. Use double for the
    // wire type so MPI_DOUBLE_INT works regardless of RealType.
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    struct { double dist; int rank; } in, out;
    in.dist  = static_cast<double>(best.first);
    in.rank  = rank;
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    double divVal = static_cast<double>(best.second);
    MPI_Bcast(&divVal, 1, MPI_DOUBLE, out.rank, MPI_COMM_WORLD);
    return static_cast<RealType>(divVal);
}

// Driver state. Built once; the divergence is then computed once, no time loop.
template<typename KeyType, typename RealType>
struct DivStepper
{
    using DomainT = ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>;

    DomainT& domain;
    int blockSize;
    int rank;
    bool useLegacy = false;     // --use-legacy: per-DOF+ownership-filter scatter

    size_t nodeCount = 0;
    size_t elementCount = 0;
    int numOwnedDofs = 0;
    int numTotalDofs = 0;

    cstone::DeviceVector<int> d_node_to_dof;
    cstone::DeviceVector<RealType> d_mass;

    cstone::DeviceVector<RealType> d_vx;
    cstone::DeviceVector<RealType> d_vy;
    cstone::DeviceVector<RealType> d_vz;

    cstone::DeviceVector<RealType> d_areaVec_x;
    cstone::DeviceVector<RealType> d_areaVec_y;
    cstone::DeviceVector<RealType> d_areaVec_z;
};

template<typename KeyType, typename RealType>
void setupDivStepper(DivStepper<KeyType, RealType>& s,
                     TestField field,
                     RealType omega,
                     RealType rotCx,
                     RealType rotCy)
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

    s.d_node_to_dof.resize(s.nodeCount);
    s.numOwnedDofs = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(), s.d_node_to_dof.data(), s.nodeCount);
    s.numTotalDofs = static_cast<int>(s.nodeCount);
    pt.lap("DOF mapping");

    // Area vectors A_f for every (element, SCS face). Same data the assembler
    // and the advection scatter consume -- single source of truth for face geometry.
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

    s.d_vx.resize(s.nodeCount);
    s.d_vy.resize(s.nodeCount);
    s.d_vz.resize(s.nodeCount);
    {
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        setTestVelocityKernel<RealType><<<nBlocks, s.blockSize>>>(
            d_x.data(), d_y.data(), d_z.data(),
            s.d_vx.data(), s.d_vy.data(), s.d_vz.data(),
            s.nodeCount, (field == TestField::Rotation ? 0 : 1),
            omega, rotCx, rotCy);
        cudaDeviceSynchronize();
    }
    pt.lap("test velocity");

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
        // Build mass via the same per-node + reverse-halo pattern so that
        // every owned node's V matches the single-rank reference exactly,
        // and so div = sum_flow / V uses a V that includes contributions
        // from neighbor ranks' local elements at corner-only-neighbor nodes.
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

    pt.report(s.rank, "setup");
}

// Compute (div u)_i at every owned node. Writes node-indexed result into d_divU
// (zero on ghosts; caller can halo-exchange if needed for VTU).
//
// Default path: per-NODE accumulator + reverseExchangeNodeHaloAdd. Conservative
// regardless of cstone halo coverage. Each element scatters +/-flow into both
// L and R local-node slots (no ownership filter). After the scatter, the
// reverse halo SUMS ghost-slot contributions on every rank back into the
// owner's slot, so the owner sees contributions from ALL elements touching the
// node globally -- including corner-only-neighbor elements that aren't in this
// rank's cstone halo.
//
// --use-legacy: keeps the old per-DOF + ownership-filter scatter for side-by-
// side comparison. Demonstrates the original 9e-2 noise floor on 4 ranks.
template<typename KeyType, typename RealType>
void computeDivergence(DivStepper<KeyType, RealType>& s,
                       cstone::DeviceVector<RealType>& d_divU)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();

    // Velocity ghosts must be in sync before scatter: face donors that involve
    // ghost corners read v on those ghosts.
    s.domain.exchangeNodeHalo(s.d_vx);
    s.domain.exchangeNodeHalo(s.d_vy);
    s.domain.exchangeNodeHalo(s.d_vz);

    d_divU.resize(s.nodeCount);

    if (s.useLegacy)
    {
        cstone::DeviceVector<RealType> d_divAcc(s.numOwnedDofs, RealType(0));

        int eBlocks = (s.elementCount + s.blockSize - 1) / s.blockSize;
        computeDivergenceKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            s.d_vx.data(), s.d_vy.data(), s.d_vz.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            d_divAcc.data(), s.elementCount);
        cudaDeviceSynchronize();

        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        applyDivergenceNormalizeKernel<RealType><<<nBlocks, s.blockSize>>>(
            d_divAcc.data(), s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            d_divU.data(), s.nodeCount);
        cudaDeviceSynchronize();
    }
    else
    {
        // Per-node accumulator: sized nodeCount, indexed by local node id. The
        // scatter writes to BOTH owned and ghost slots (no ownership filter)
        // but iterates LOCAL elements only so each global face is scattered
        // exactly once globally (by the rank owning the element).
        cstone::DeviceVector<RealType> d_divAccNode(s.nodeCount, RealType(0));

        size_t startElem = s.domain.startIndex();
        size_t numLocal  = s.domain.localElementCount();

        int eBlocks = numLocal > 0 ? int((numLocal + s.blockSize - 1) / s.blockSize) : 0;
        if (eBlocks > 0)
        {
            computeDivergencePerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
                std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
                std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
                std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
                s.d_vx.data(), s.d_vy.data(), s.d_vz.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                d_divAccNode.data(), startElem, numLocal);
            cudaDeviceSynchronize();
        }

        // Sum ghost-slot contributions back into the owning rank's slot. This
        // is the closure step that makes the scatter globally conservative
        // even when cstone halos miss corner-only-neighbor elements.
        s.domain.reverseExchangeNodeHaloAdd(d_divAccNode);

        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        applyDivergenceNormalizePerNodeKernel<RealType><<<nBlocks, s.blockSize>>>(
            d_divAccNode.data(), s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            d_divU.data(), s.nodeCount);
        cudaDeviceSynchronize();
    }

    // Make ghost slots carry the owner-rank values so VTU output looks
    // continuous across rank boundaries.
    s.domain.exchangeNodeHalo(d_divU);
}

// ============================================================================
// Diagnostic driver: face-touch counts, per-rank net-flow comparison, and
// per-face contribution dump for one owned target node.
// ============================================================================
template<typename KeyType, typename RealType>
void runDivergenceDiagnostics(DivStepper<KeyType, RealType>& s,
                              RealType probeX, RealType probeY, RealType probeZ)
{
    int rank = s.rank;
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_x             = s.domain.getNodeX();
    const auto& d_y             = s.domain.getNodeY();
    const auto& d_z             = s.domain.getNodeZ();

    // Velocity ghosts in sync (mirrors what computeDivergence does first).
    s.domain.exchangeNodeHalo(s.d_vx);
    s.domain.exchangeNodeHalo(s.d_vy);
    s.domain.exchangeNodeHalo(s.d_vz);

    if (rank == 0)
    {
        std::cout << "----------------------------------------\n";
        std::cout << "[diagnose] cross-rank scatter audit\n";
        std::cout << "----------------------------------------\n";
    }

    // ----- Diagnostic 1: face-touch count per owned node ---------------------
    // Expected: every owned interior node sees 12 SCS-face touches. Boundary
    // owned nodes see fewer (6 on a face, 3 on an edge, 1 at a corner) because
    // fewer hex elements abut them. If a SINGLE-rank reference run reports 12
    // for some owned node and the multi-rank run reports <12 for the SAME
    // global node, halo coverage (cause #1) is at fault.
    cstone::DeviceVector<int> d_faceTouchNode(s.nodeCount, 0);
    {
        int eBlocks = (s.elementCount + s.blockSize - 1) / s.blockSize;
        countFacesTouchedKernel<KeyType><<<eBlocks, s.blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            d_nodeOwnership.data(),
            d_faceTouchNode.data(), s.elementCount);
        cudaDeviceSynchronize();
    }

    // Host-side histogram + min/max over owned nodes ONLY (touch=0 on ghosts).
    std::vector<int>     h_touch(s.nodeCount);
    std::vector<uint8_t> h_own(s.nodeCount);
    cudaMemcpy(h_touch.data(), d_faceTouchNode.data(), s.nodeCount * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_own.data(),   d_nodeOwnership.data(), s.nodeCount * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    int locMin = std::numeric_limits<int>::max();
    int locMax = 0;
    long long locSum = 0;
    int locOwned = 0;
    std::vector<int> hist(32, 0);
    for (size_t i = 0; i < s.nodeCount; ++i)
    {
        if (h_own[i] != 1) continue;
        int t = h_touch[i];
        if (t < locMin) locMin = t;
        if (t > locMax) locMax = t;
        locSum += t;
        ++locOwned;
        if (t < 32) hist[t]++;
    }
    if (locOwned == 0) locMin = 0;

    int gMin = 0, gMax = 0, gOwned = 0;
    long long gSum = 0;
    MPI_Allreduce(&locMin,   &gMin,   1, MPI_INT,       MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&locMax,   &gMax,   1, MPI_INT,       MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locSum,   &gSum,   1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locOwned, &gOwned, 1, MPI_INT,       MPI_SUM, MPI_COMM_WORLD);

    std::vector<int> gHist(32, 0);
    MPI_Allreduce(hist.data(), gHist.data(), 32, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "[D1] Face-touch count per owned node (globally aggregated):\n";
        std::cout << "      total owned DOFs:  " << gOwned  << "\n";
        std::cout << "      min:               " << gMin    << "\n";
        std::cout << "      max:               " << gMax    << "\n";
        std::cout << "      mean (sum/owned):  " << (gOwned ? double(gSum) / double(gOwned) : 0.0) << "\n";
        std::cout << "      histogram (touches : count_of_owned_nodes):\n";
        for (int t = 0; t < 32; ++t)
        {
            if (gHist[t] == 0) continue;
            std::cout << "        " << std::setw(3) << t << " : " << gHist[t] << "\n";
        }
        std::cout << "      Reference (single rank): every interior node = 12,\n";
        std::cout << "      face nodes = 6, edge nodes = 3, corner nodes = 1. If multi-rank\n";
        std::cout << "      histogram has counts at values BELOW the single-rank counts for\n";
        std::cout << "      the same node positions, halo coverage (cause #1) is the bug.\n\n";
    }

    // ----- Diagnostic 2: per-rank per-node net flow vs. owner-only-sum -------
    // Run the per-node scatter (without reverse halo) on LOCAL elements only,
    // then take owner-side slot values. Compare to: same scatter + reverse
    // halo (the fix). Reductions of |delta| over owned DOFs quantify how much
    // mass the reverse-halo recovers from neighbor ranks' local elements.
    size_t startElem = s.domain.startIndex();
    size_t numLocal  = s.domain.localElementCount();

    cstone::DeviceVector<RealType> d_divAccNode_noClose(s.nodeCount, RealType(0));
    if (numLocal > 0)
    {
        int eBlocks = int((numLocal + s.blockSize - 1) / s.blockSize);
        perNodeNetFlowKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            s.d_vx.data(), s.d_vy.data(), s.d_vz.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            d_divAccNode_noClose.data(), startElem, numLocal);
        cudaDeviceSynchronize();
    }

    cstone::DeviceVector<RealType> d_divAccNode_closed(s.nodeCount, RealType(0));
    if (numLocal > 0)
    {
        int eBlocks = int((numLocal + s.blockSize - 1) / s.blockSize);
        perNodeNetFlowKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            s.d_vx.data(), s.d_vy.data(), s.d_vz.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            d_divAccNode_closed.data(), startElem, numLocal);
        cudaDeviceSynchronize();
    }
    s.domain.reverseExchangeNodeHaloAdd(d_divAccNode_closed);

    // Sum |delta| = |closed - noClose| over owned slots. Non-zero -> the
    // reverse halo delivered missing contributions (i.e. cause #1 was real).
    std::vector<RealType> h_no(s.nodeCount), h_cl(s.nodeCount);
    cudaMemcpy(h_no.data(), d_divAccNode_noClose.data(), s.nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_cl.data(), d_divAccNode_closed.data(),  s.nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);

    RealType locDeltaL2 = 0, locDeltaMax = 0;
    int      locDeltaPos = 0;
    for (size_t i = 0; i < s.nodeCount; ++i)
    {
        if (h_own[i] != 1) continue;
        RealType d = std::fabs(h_cl[i] - h_no[i]);
        locDeltaL2  += d * d;
        if (d > locDeltaMax) locDeltaMax = d;
        if (d > RealType(1e-14)) ++locDeltaPos;
    }
    RealType gDeltaL2 = 0, gDeltaMax = 0;
    int      gDeltaPos = 0;
    auto mpiTypeD = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&locDeltaL2,  &gDeltaL2,  1, mpiTypeD, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locDeltaMax, &gDeltaMax, 1, mpiTypeD, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locDeltaPos, &gDeltaPos, 1, MPI_INT,  MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "[D2] |closed - noClose| reduction over OWNED nodes:\n";
        std::cout << "      L2 of delta:     " << std::scientific << std::setprecision(3) << std::sqrt(gDeltaL2) << "\n";
        std::cout << "      max of delta:    " << gDeltaMax << "\n";
        std::cout << "      owned w/ |d|>1e-14: " << gDeltaPos << " of " << gOwned << "\n";
        std::cout << "      delta = (reverse-halo contribution recovered per owned node).\n";
        std::cout << "      If max(delta) is O(1) and ~boundary-touching count of owned slots\n";
        std::cout << "      is non-zero, cause #1 (halo coverage) is the bug.\n";
        std::cout << "      If max(delta) is exactly zero everywhere, the bug is bit-level\n";
        std::cout << "      (cause #2/3/4) and per-face dump (D3) will localize it.\n\n";
        std::cout << std::defaultfloat;
    }

    // ----- Diagnostic 3: per-face dump for one target owned node -------------
    // Pick the owned node closest to (probeX, probeY, probeZ) on each rank.
    // Use MINLOC across ranks; the winning rank dumps its face contributions.
    // Same probe coords passed on 1-rank and 4-rank lets us compare contributions.
    struct Probe { double dist; int rank; };
    ProbeLocalNodeAt<RealType> op{probeX, probeY, probeZ,
                                  d_x.data(), d_y.data(), d_z.data()};
    auto init = thrust::make_pair(RealType(1e30), size_t(0));
    auto best = thrust::transform_reduce(
        thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount),
        op, init, PairMinByFirst<RealType>());

    // Check if best is owned; if not, set dist=1e30 so another rank wins.
    size_t bestLocal = best.second;
    double bestDist  = double(best.first);
    if (bestLocal >= s.nodeCount || h_own[bestLocal] != 1)
        bestDist = 1e30;

    Probe in{bestDist, rank}, out{};
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

    // Every rank dumps its face contributions to the SAME global node (whether
    // it's owned here or not). The winning rank publishes the target coords;
    // others find a local node at the same (xb, yb, zb). This is the per-face
    // audit -- if rank A and B both see the same physical face, their flow
    // values should sum (across L and R contributions) to zero.
    double targetXYZ[3];
    if (rank == out.rank)
    {
        std::vector<RealType> h_x(s.nodeCount), h_y(s.nodeCount), h_z(s.nodeCount);
        cudaMemcpy(h_x.data(), d_x.data(), s.nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_y.data(), d_y.data(), s.nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_z.data(), d_z.data(), s.nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
        targetXYZ[0] = double(h_x[bestLocal]);
        targetXYZ[1] = double(h_y[bestLocal]);
        targetXYZ[2] = double(h_z[bestLocal]);
    }
    MPI_Bcast(targetXYZ, 3, MPI_DOUBLE, out.rank, MPI_COMM_WORLD);

    // Find the local node closest to the broadcasted physical position on every rank.
    ProbeLocalNodeAt<RealType> op2{RealType(targetXYZ[0]), RealType(targetXYZ[1]), RealType(targetXYZ[2]),
                                   d_x.data(), d_y.data(), d_z.data()};
    auto best2 = thrust::transform_reduce(
        thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount),
        op2, init, PairMinByFirst<RealType>());
    KeyType targetLocal = KeyType(best2.second);
    RealType d2 = best2.first;

    constexpr int MAX_DUMP = 64;
    cstone::DeviceVector<int>      d_outCount(1, 0);
    cstone::DeviceVector<int>      d_outElem(MAX_DUMP, 0);
    cstone::DeviceVector<int>      d_outIp(MAX_DUMP, 0);
    cstone::DeviceVector<int>      d_outRoleIsL(MAX_DUMP, 0);
    cstone::DeviceVector<int>      d_outOwnL(MAX_DUMP, 0);
    cstone::DeviceVector<int>      d_outOwnR(MAX_DUMP, 0);
    cstone::DeviceVector<RealType> d_outFlow(MAX_DUMP, RealType(0));
    cstone::DeviceVector<RealType> d_outAx(MAX_DUMP, RealType(0));
    cstone::DeviceVector<RealType> d_outAy(MAX_DUMP, RealType(0));
    cstone::DeviceVector<RealType> d_outAz(MAX_DUMP, RealType(0));

    // Only dump if this rank actually has a node within ~1e-6 of the target
    // physical position (rules out far-away mismatches when the partition
    // is far from the target and no local node lies near).
    bool dumpThisRank = (double(d2) < 1e-6 * 1e-6);

    if (dumpThisRank && numLocal > 0)
    {
        int eBlocks = int((numLocal + s.blockSize - 1) / s.blockSize);
        dumpFaceContribsForNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            std::get<0>(d_conn).data(), std::get<1>(d_conn).data(),
            std::get<2>(d_conn).data(), std::get<3>(d_conn).data(),
            std::get<4>(d_conn).data(), std::get<5>(d_conn).data(),
            std::get<6>(d_conn).data(), std::get<7>(d_conn).data(),
            s.d_vx.data(), s.d_vy.data(), s.d_vz.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            d_nodeOwnership.data(),
            targetLocal,
            d_outCount.data(),
            d_outElem.data(), d_outIp.data(), d_outRoleIsL.data(),
            d_outOwnL.data(), d_outOwnR.data(),
            d_outFlow.data(),
            d_outAx.data(), d_outAy.data(), d_outAz.data(),
            startElem, numLocal, MAX_DUMP);
        cudaDeviceSynchronize();
    }

    int h_count = 0;
    if (dumpThisRank) cudaMemcpy(&h_count, d_outCount.data(), sizeof(int), cudaMemcpyDeviceToHost);
    int h_clamped = std::min(h_count, MAX_DUMP);

    std::vector<int>      hE(h_clamped), hIp(h_clamped), hRoleL(h_clamped), hOL(h_clamped), hOR(h_clamped);
    std::vector<RealType> hFlow(h_clamped), hAx(h_clamped), hAy(h_clamped), hAz(h_clamped);
    if (h_clamped > 0)
    {
        cudaMemcpy(hE.data(),     d_outElem.data(),    h_clamped * sizeof(int),      cudaMemcpyDeviceToHost);
        cudaMemcpy(hIp.data(),    d_outIp.data(),      h_clamped * sizeof(int),      cudaMemcpyDeviceToHost);
        cudaMemcpy(hRoleL.data(), d_outRoleIsL.data(), h_clamped * sizeof(int),      cudaMemcpyDeviceToHost);
        cudaMemcpy(hOL.data(),    d_outOwnL.data(),    h_clamped * sizeof(int),      cudaMemcpyDeviceToHost);
        cudaMemcpy(hOR.data(),    d_outOwnR.data(),    h_clamped * sizeof(int),      cudaMemcpyDeviceToHost);
        cudaMemcpy(hFlow.data(),  d_outFlow.data(),    h_clamped * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hAx.data(),    d_outAx.data(),      h_clamped * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hAy.data(),    d_outAy.data(),      h_clamped * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hAz.data(),    d_outAz.data(),      h_clamped * sizeof(RealType), cudaMemcpyDeviceToHost);
    }

    if (rank == 0)
    {
        std::cout << "[D3] Per-face contributions to the owned node nearest ("
                  << std::fixed << std::setprecision(4)
                  << targetXYZ[0] << ", " << targetXYZ[1] << ", " << targetXYZ[2] << ")\n";
        std::cout << "      Each rank prints its locally-seen faces touching this node.\n";
        std::cout << "      Per-face row format:\n";
        std::cout << "        rank  isOwner  e_idx  ip  role(L/R)  ownL  ownR  signed_flow  |A|\n";
        std::cout << "      For correctness the GLOBAL sum of signed_flow (sign by role:\n";
        std::cout << "      +flow if role=L, -flow if role=R) at this node must equal the\n";
        std::cout << "      surface flux through the node's CV. For an interior node and\n";
        std::cout << "      divergence-free v that is exactly zero.\n";
        std::cout << std::defaultfloat;
    }

    // Serial print per rank (rank 0 first, then 1, 2, ...).
    for (int r = 0; r < s.domain.numRanks(); ++r)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (r != rank) continue;
        if (!dumpThisRank)
        {
            std::cout << "      [rank " << r << "] no local node within 1e-6 of probe; skipped\n" << std::flush;
            continue;
        }
        bool isOwner = (rank == out.rank);
        std::cout << "      [rank " << r << "] localNode=" << targetLocal
                  << " isOwner=" << (isOwner ? 1 : 0)
                  << " dump_count=" << h_count
                  << (h_count > MAX_DUMP ? " (CLAMPED)" : "") << "\n";

        double localSignedSum = 0.0;
        for (int i = 0; i < h_clamped; ++i)
        {
            double s_flow = (hRoleL[i] == 1) ? hFlow[i] : -hFlow[i];
            localSignedSum += s_flow;
            double aMag = std::sqrt(double(hAx[i])*hAx[i] + double(hAy[i])*hAy[i] + double(hAz[i])*hAz[i]);
            std::cout << "        " << std::setw(3) << r << "  "
                      << std::setw(7) << (isOwner ? 1 : 0) << "  "
                      << std::setw(6) << hE[i] << "  "
                      << std::setw(2) << hIp[i] << "  "
                      << std::setw(7) << (hRoleL[i] ? "L" : "R") << "  "
                      << std::setw(4) << hOL[i] << "  "
                      << std::setw(4) << hOR[i] << "  "
                      << std::scientific << std::setprecision(6)
                      << std::setw(15) << s_flow << "  "
                      << std::setw(13) << aMag << "\n"
                      << std::defaultfloat;
        }
        std::cout << "      [rank " << r << "] local signed-flow sum: "
                  << std::scientific << std::setprecision(6) << localSignedSum
                  << std::defaultfloat << "\n" << std::flush;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) std::cout << "----------------------------------------\n\n";
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
    int blockSize    = 256;
    int bucketSize   = 64;
    TestField field  = TestField::Rotation;
    double omega     = 2.0;
    std::string vtuPrefix;
    bool diagnose    = false;      // --diagnose: run cross-rank scatter audit
    bool useLegacy   = false;      // --use-legacy: per-DOF+ownership-filter scatter
    double probeX    = 0.5;
    double probeY    = 0.5;
    double probeZ    = 0.5;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg.find("--mesh=") == 0)
            meshFile = arg.substr(7);
        else if (arg.find("--block-size=") == 0)
            blockSize = std::stoi(arg.substr(13));
        else if (arg.find("--bucket-size=") == 0)
            bucketSize = std::stoi(arg.substr(14));
        else if (arg.find("--field=") == 0)
        {
            std::string v = arg.substr(8);
            if (v == "rotation")    field = TestField::Rotation;
            else if (v == "sinx")   field = TestField::Sinx;
            else
            {
                if (rank == 0)
                    std::cerr << "Error: --field must be rotation or sinx, got '" << v << "'\n";
                MPI_Finalize();
                return 1;
            }
        }
        else if (arg.find("--omega=") == 0)
            omega = std::stod(arg.substr(8));
        else if (arg.find("--vtu-output=") == 0)
            vtuPrefix = arg.substr(13);
        else if (arg == "--diagnose")
            diagnose = true;
        else if (arg == "--use-legacy")
            useLegacy = true;
        else if (arg.find("--probe-x=") == 0)
            probeX = std::stod(arg.substr(10));
        else if (arg.find("--probe-y=") == 0)
            probeY = std::stod(arg.substr(10));
        else if (arg.find("--probe-z=") == 0)
            probeZ = std::stod(arg.substr(10));
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
            std::cout << "  --field=KIND         rotation (default) | sinx\n";
            std::cout << "  --omega=VALUE        Rotation rate for --field=rotation (default: 2.0)\n";
            std::cout << "  --vtu-output=PREFIX  Write a single-frame VTU/PVTU/PVD of divU (default: off)\n";
            std::cout << "  --bucket-size=N      Cornerstone bucket size (default: 64)\n";
            std::cout << "  --block-size=N       CUDA block size (default: 256)\n";
            std::cout << "  --diagnose           Run cross-rank scatter audit (face-touch count + per-face dump)\n";
            std::cout << "  --use-legacy         Use the original per-DOF + ownership-filter scatter (bug demo)\n";
            std::cout << "  --probe-x/y/z=VAL    Physical coordinate of the node to dump in --diagnose mode\n";
            std::cout << "                       (default 0.5,0.5,0.5 — the cube center; pick a node near a\n";
            std::cout << "                        rank-partition corner for cause-#1 evidence)\n";
        }
        MPI_Finalize();
        return 1;
    }

    using KeyType  = uint64_t;
    using RealType = double;

    const RealType rotCx = 0.5;
    const RealType rotCy = 0.5;

    if (rank == 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "MARS Divergence Operator Validator (CVFEM SCS, Phase B.3)\n";
        std::cout << "========================================\n";
        if (field == TestField::Rotation)
        {
            std::cout << "Field:    rotation, v = (-omega*(y-0.5), omega*(x-0.5), 0)\n";
            std::cout << "omega:    " << omega << "\n";
            std::cout << "Expected: div(v) = 0 exactly\n";
        }
        else
        {
            std::cout << "Field:    sinx, v = (sin(pi*x), 0, 0)\n";
            std::cout << "Expected: div(v) = pi*cos(pi*x)\n";
        }
        std::cout << "Mesh:     " << meshFile << "\n";
        std::cout << "MPI ranks: " << numRanks << "\n";
        std::cout << "========================================\n\n";
    }

    // maxLevels=0 so we never refine; AmrManager is here purely for the
    // mesh-load + lazy halo/adjacency/coord build path.
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

    DivStepper<KeyType, RealType> s{amr.domain(), blockSize, rank, useLegacy};
    setupDivStepper<KeyType, RealType>(s, field, RealType(omega), rotCx, rotCy);

    if (rank == 0)
    {
        std::cout << "Scatter mode: " << (useLegacy ? "LEGACY (per-DOF + ownership filter)"
                                                    : "PER-NODE + reverse halo (conservative)") << "\n\n";
    }

    if (diagnose)
    {
        runDivergenceDiagnostics<KeyType, RealType>(s, RealType(probeX), RealType(probeY), RealType(probeZ));
    }

    cstone::DeviceVector<RealType> d_divU;
    auto t0 = std::chrono::high_resolution_clock::now();
    computeDivergence<KeyType, RealType>(s, d_divU);
    cudaDeviceSynchronize();
    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = std::chrono::high_resolution_clock::now();
    float divMs = std::chrono::duration<float, std::milli>(t1 - t0).count();

    // -------- Reductions: |div u|_L2 split into interior vs boundary, plus
    // analytical L2 and L2 error on interior owned DOFs only. The CVFEM CV
    // divergence at boundary nodes captures the boundary face flux, so e.g.
    // for rigid rotation the boundary |div|_L2 is intentionally non-zero
    // (through-flow on the cube faces) and only the interior norm is a fair
    // operator test against div=0. --------------------------------
    cstone::DeviceVector<RealType> d_errSqInt(s.numOwnedDofs, RealType(0));
    cstone::DeviceVector<RealType> d_anaSqInt(s.numOwnedDofs, RealType(0));
    cstone::DeviceVector<RealType> d_numSqInt(s.numOwnedDofs, RealType(0));
    cstone::DeviceVector<RealType> d_numSqBnd(s.numOwnedDofs, RealType(0));
    cstone::DeviceVector<uint8_t>  d_bndMask (s.numOwnedDofs, uint8_t(0));
    {
        const auto& d_nodeOwn = amr.domain().getNodeOwnershipMap();
        const auto& d_x = amr.domain().getNodeX();
        const auto& d_y = amr.domain().getNodeY();
        const auto& d_z = amr.domain().getNodeZ();
        int nBlocks = (s.nodeCount + blockSize - 1) / blockSize;
        computeDivErrorKernel<RealType><<<nBlocks, blockSize>>>(
            d_x.data(), d_y.data(), d_z.data(),
            d_divU.data(), s.d_node_to_dof.data(), d_nodeOwn.data(),
            s.d_mass.data(),
            d_errSqInt.data(), d_anaSqInt.data(), d_numSqInt.data(),
            d_numSqBnd.data(), d_bndMask.data(),
            s.nodeCount, (field == TestField::Rotation ? 0 : 1));
        cudaDeviceSynchronize();
    }

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    auto sumLocal = [&](const cstone::DeviceVector<RealType>& v)
    {
        auto p = thrust::device_pointer_cast(v.data());
        return thrust::reduce(thrust::device, p, p + s.numOwnedDofs, RealType(0));
    };
    RealType locErrInt = sumLocal(d_errSqInt);
    RealType locAnaInt = sumLocal(d_anaSqInt);
    RealType locNumInt = sumLocal(d_numSqInt);
    RealType locNumBnd = sumLocal(d_numSqBnd);
    RealType gErrInt = 0, gAnaInt = 0, gNumInt = 0, gNumBnd = 0;
    MPI_Allreduce(&locErrInt, &gErrInt, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locAnaInt, &gAnaInt, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locNumInt, &gNumInt, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locNumBnd, &gNumBnd, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    RealType gNumAll = gNumInt + gNumBnd;

    // Per-owned-DOF |divU| max, split into interior and boundary.
    RealType locMaxInt = 0, locMaxBnd = 0;
    {
        const auto& d_nodeOwn = amr.domain().getNodeOwnershipMap();
        const uint8_t* ownPtr = d_nodeOwn.data();
        const int* dofPtr     = s.d_node_to_dof.data();
        const RealType* dPtr  = d_divU.data();
        const uint8_t* bndPtr = d_bndMask.data();
        locMaxInt = thrust::transform_reduce(
            thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [ownPtr, dofPtr, dPtr, bndPtr] __device__(size_t i) -> RealType {
                if (ownPtr[i] != 1) return RealType(0);
                int dof = dofPtr[i];
                if (dof < 0) return RealType(0);
                if (bndPtr[dof] != 0) return RealType(0);
                return fabs(dPtr[i]);
            },
            RealType(0), thrust::maximum<RealType>());
        locMaxBnd = thrust::transform_reduce(
            thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [ownPtr, dofPtr, dPtr, bndPtr] __device__(size_t i) -> RealType {
                if (ownPtr[i] != 1) return RealType(0);
                int dof = dofPtr[i];
                if (dof < 0) return RealType(0);
                if (bndPtr[dof] == 0) return RealType(0);
                return fabs(dPtr[i]);
            },
            RealType(0), thrust::maximum<RealType>());
    }
    RealType gMaxInt = 0, gMaxBnd = 0;
    MPI_Allreduce(&locMaxInt, &gMaxInt, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locMaxBnd, &gMaxBnd, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "Divergence computed in " << std::fixed << std::setprecision(2)
                  << divMs << " ms\n\n";
        std::cout << "========================================\n";
        if (field == TestField::Rotation)
            std::cout << "Test: rotation, omega=" << omega << "\n";
        else
            std::cout << "Test: sinx\n";
        std::cout << std::scientific << std::setprecision(3);
        std::cout << "  Numerical  |div u|_L2 (all owned):     " << std::sqrt(gNumAll)
                  << "   (includes boundary CV through-flow)\n";
        std::cout << "  Numerical  |div u|_L2 (interior only): " << std::sqrt(gNumInt)
                  << "   (true bulk divergence)\n";
        std::cout << "  Analytical |div u|_L2:                 " << std::sqrt(gAnaInt) << "\n";
        std::cout << "  L2 error (interior only):              " << std::sqrt(gErrInt) << "\n";
        std::cout << "  Max |div u| interior:                  " << gMaxInt << "\n";
        std::cout << "  Max |div u| boundary:                  " << gMaxBnd << "\n";
        std::cout << std::defaultfloat;
    }

    if (field == TestField::Sinx)
    {
        const auto& d_x = amr.domain().getNodeX();
        const auto& d_y = amr.domain().getNodeY();
        const auto& d_z = amr.domain().getNodeZ();
        const auto& d_nodeOwn = amr.domain().getNodeOwnershipMap();
        RealType d05 = probeDivAt<RealType>(d_x, d_y, d_z, d_divU,
                                            s.d_node_to_dof, d_nodeOwn, s.nodeCount,
                                            RealType(0.5), RealType(0.5), RealType(0.5));
        RealType d025 = probeDivAt<RealType>(d_x, d_y, d_z, d_divU,
                                             s.d_node_to_dof, d_nodeOwn, s.nodeCount,
                                             RealType(0.25), RealType(0.5), RealType(0.5));
        RealType d075 = probeDivAt<RealType>(d_x, d_y, d_z, d_divU,
                                             s.d_node_to_dof, d_nodeOwn, s.nodeCount,
                                             RealType(0.75), RealType(0.5), RealType(0.5));
        if (rank == 0)
        {
            constexpr double PI = 3.14159265358979323846;
            std::cout << std::fixed << std::setprecision(4);
            std::cout << "  At x=0.5  div: " << std::setw(8) << d05
                      << "   (expected " << std::setw(8) << (PI * std::cos(PI * 0.5)) << ")\n";
            std::cout << "  At x=0.25 div: " << std::setw(8) << d025
                      << "   (expected " << std::setw(8) << (PI * std::cos(PI * 0.25)) << ")\n";
            std::cout << "  At x=0.75 div: " << std::setw(8) << d075
                      << "   (expected " << std::setw(8) << (PI * std::cos(PI * 0.75)) << ")\n";
            std::cout << std::defaultfloat;
        }
    }
    if (rank == 0) std::cout << "========================================\n";

    if (!vtuPrefix.empty())
    {
        fem::VTUParallelWriter<KeyType, RealType> vtuWriter(vtuPrefix);
        vtuWriter.writeFrame(0, 0.0, amr.domain(), d_divU, "divU");
        if (rank == 0)
        {
            std::cout << "\nVTU written: " << vtuPrefix << "_step0000.pvtu, "
                      << vtuPrefix << ".pvd\n";
        }
    }

    MPI_Finalize();
    return 0;
}
