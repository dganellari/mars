// Galerkin Laplacian A = D M^{-1} D^T standalone validator (Phase E, Stages 1+2).
//
// The NS projection step solves a pressure-Poisson system with an operator
// formed algebraically as A = D M^{-1} D^T, where
//   D    : the per-node + reverse-halo discrete CVFEM divergence (Phase B.3,
//          same kernel as mars_amr_div.cu's computeDivergencePerNodeKernel).
//   M    : lumped mass (diagonal, per owned DOF) -- same per-node + reverse-halo
//          construction as mars_amr_div.cu, so V at every owned node matches
//          single-rank exactly.
//   D^T  : the *algebraic* transpose of D. Per face f connecting L,R with area
//          vector A_f pointing L -> R:
//            D^T phi accumulates  +0.5 * A_f * (p[L] - p[R])  into BOTH endpoints
//            with the SAME sign (NOT the opposite-sign scatter of the CVFEM
//            SCS gradient in mars_amr_grad.cu's computeGradientPerNodeKernel).
//            This is what applyDivTransposePerNodeKernel in mars_amr_ns_projection.cu
//            already implements; we copy that scatter pattern here byte for byte.
//
// A is symmetric (M diagonal => M^{-1} symmetric; sandwich preserves symmetry)
// and positive semi-definite. Its null space is the constant pressure mode
// (D^T 1 = 0 face-by-face since p[L]=p[R]=1 makes dp=0). A*phi at every owned
// node is left in *un-normalized* (V-scaled) form -- the same shape used by the
// NS driver's pressure-Poisson SpMV, so b vs Ax magnitudes match downstream.
//
// Stage 1 validation: NO wiring into the NS driver. Standalone tests:
//   --test=linear     A*phi for phi=a*x+b*y+c*z; expected interior |Aphi|=roundoff
//   --test=sign       A*phi for phi=x^2; the SIGNED interior value (A*phi)/V_i
//                     equals +2 or -2 -- isolates the global sign of A which
//                     PSD/CG tests cannot see (they accept both A and -A).
//   --test=symmetry   u^T A v == v^T A u to roundoff for random u,v
//   --test=nullspace  A*1 = 0 to roundoff (constant in null space)
//   --test=psd        phi_mean-removed random phi -> phi^T A phi > 0 strictly
//   --test=cg         un-preconditioned matrix-free CG converges to a sin*sin*sin
//                     reference field (range of A by construction)
//
// Stage 2 addition: Jacobi-preconditioned CG (PCG) on the same RHS.
//   --test=cg-jacobi  matrix-free PCG with M_jac = diag(A); same convergence
//                     check as Test 5a, expected to drop iteration count
//                     significantly (~10-20 vs 44 un-preconditioned on cube16/4-rank).
//
// diag(A) is computed via brute-force probe-vector method: numOwnedDofs separate
// A-applies on unit vectors. This is O(N^2) and only OK for cube16-sized cases
// (~5000 DOFs -> seconds). For production use this must be replaced by either
// graph-coloured probes (~27 applies on a uniform hex stencil) or a dedicated
// analytic per-face self-contribution kernel.
//
// All tests use the same per-node + reverse-halo conservation closure as
// mars_amr_div.cu and mars_amr_grad.cu so they're correct on multi-rank.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/amr/mars_amr.hpp"

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/inner_product.h>
#include <thrust/system/cuda/execution_policy.h>

#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>
#include <string>
#include <algorithm>

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

enum class TestKind { Linear, Sign, Symmetry, NullSpace, Psd, Cg, CgJacobi, All };

// Optional coordinate stretching applied to ALL nodes (owned + ghosts) after
// cstone has set up the domain. Used to test Jacobi PCG on non-uniform meshes
// where diag(A) actually varies. Stretch is applied to the x coordinate only,
// taking the original x in [xmin, xmax] and mapping it onto the same range
// with a chosen non-linear transform. Y and Z are left unchanged.
enum class StretchKind { None, SinX, Exp, Tanh };

// Per-node coordinate transform. Applied in-place on the cstone-cached device
// arrays before area vectors / lumped mass are computed, so everything
// downstream (D, D^T, M, diag(A)) automatically sees the stretched cells.
//
//   SinX:  x' = x + 0.1*(xmax-xmin)*sin(2*pi*(x-xmin)/(xmax-xmin))
//          ~3x compression at the antinodes, ~3x expansion at the nodes.
//   Exp:   x' = xmin + (xmax-xmin) * (exp(alpha*t) - 1) / (exp(alpha) - 1)
//          where t = (x-xmin)/(xmax-xmin) and alpha = 3.
//          Boundary-layer style: ~20x cell-size variation from xmin to xmax.
//   Tanh:  x' = xmin + 0.5*(xmax-xmin) * (1 + tanh(beta*(2*t-1))/tanh(beta))
//          where t = (x-xmin)/(xmax-xmin) and beta = 2.
//          Cluster near both ends (two boundary layers).
template<typename RealType>
__global__ void applyCoordStretchKernel(RealType* nodeX,
                                        size_t numNodes,
                                        int kind,
                                        RealType xmin,
                                        RealType xmax)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    RealType x = nodeX[i];
    RealType range = xmax - xmin;
    if (range <= RealType(0)) return;
    RealType t = (x - xmin) / range;
    RealType xNew = x;
    if (kind == 1) // SinX
    {
        const RealType amp = RealType(0.1);
        xNew = x + amp * range * sin(RealType(2 * M_PI) * t);
    }
    else if (kind == 2) // Exp
    {
        const RealType alpha = RealType(3);
        RealType eAlpha = exp(alpha);
        xNew = xmin + range * (exp(alpha * t) - RealType(1)) / (eAlpha - RealType(1));
    }
    else if (kind == 3) // Tanh
    {
        const RealType beta = RealType(2);
        RealType tB = tanh(beta);
        xNew = xmin + RealType(0.5) * range * (RealType(1) + tanh(beta * (RealType(2) * t - RealType(1))) / tB);
    }
    nodeX[i] = xNew;
}

// Sync on lap so the wall-clock includes every rank's work, not just rank-0's
// stragglers. Mirrors mars_amr_div.cu / mars_amr_grad.cu.
struct PhaseTimer
{
    using clk = std::chrono::high_resolution_clock;
    clk::time_point t0;
    std::vector<std::pair<std::string, float>> phases;
    explicit PhaseTimer() { reset(); }
    void reset() { cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD); t0 = clk::now(); }
    void lap(const std::string& name)
    {
        cudaDeviceSynchronize(); MPI_Barrier(MPI_COMM_WORLD);
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
            std::cout << "    " << std::left << std::setw(28) << n
                      << std::right << std::fixed << std::setprecision(2)
                      << std::setw(10) << ms << " ms ("
                      << std::setw(5) << std::setprecision(1)
                      << (100.0f * ms / std::max(total, 1e-3f)) << "%)\n";
        std::cout << "    " << std::left << std::setw(28) << "TOTAL"
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(10) << total << " ms\n";
    }
};

// =============================================================================
// Lumped mass: per-node + reverse-halo closure, gathered into per-owned-DOF
// array. Identical pattern to mars_amr_div.cu's mass build.
// =============================================================================

template<typename KeyType, typename RealType>
__global__ void computeLumpedMassPerNodeKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    RealType* massNode, size_t startElem, size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;

    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};
    RealType x[8], y[8], z[8];
    for (int i = 0; i < 8; ++i) { x[i] = nodeX[n[i]]; y[i] = nodeY[n[i]]; z[i] = nodeZ[n[i]]; }

    RealType xmin = x[0], xmax = x[0], ymin = y[0], ymax = y[0], zmin = z[0], zmax = z[0];
    for (int i = 1; i < 8; ++i) {
        xmin = fmin(xmin, x[i]); xmax = fmax(xmax, x[i]);
        ymin = fmin(ymin, y[i]); ymax = fmax(ymax, y[i]);
        zmin = fmin(zmin, z[i]); zmax = fmax(zmax, z[i]);
    }
    RealType contrib = (xmax - xmin) * (ymax - ymin) * (zmax - zmin) * RealType(0.125);

    for (int i = 0; i < 8; ++i) atomicAdd(&massNode[n[i]], contrib);
}

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

// =============================================================================
// D^T scatter (the *algebraic* transpose of D). For each SCS face of element e
// between corners (iL, iR) with area vector A_f pointing L -> R:
//   dp = 0.5 * (phi[iL] - phi[iR])
//   gAcc[iL] += dp * A_f
//   gAcc[iR] += dp * A_f      (SAME sign as L, NOT opposite)
// This is byte-for-byte applyDivTransposePerNodeKernel from
// mars_amr_ns_projection.cu. The block-transpose of D's per-face 2x2 acts on
// the *difference* of endpoints; the SCS gradient acts on the *sum* -- they
// are different operators. Use this one (it's the transpose of div).
// =============================================================================

template<typename KeyType, typename RealType>
__global__ void applyDivTransposePerNodeKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
    const RealType* phi,
    const RealType* areaVecX, const RealType* areaVecY, const RealType* areaVecZ,
    RealType* gxAccNode, RealType* gyAccNode, RealType* gzAccNode,
    size_t startElem, size_t numLocal)
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
        RealType dp = RealType(0.5) * (phi[iL] - phi[iR]);

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

// Divide accumulator by V on owned nodes; ghosts zeroed (will be filled by a
// subsequent forward halo exchange). In-place safe (same-index gather).
template<typename RealType>
__global__ void normalizeByMassPerNodeKernel(const RealType* gxAccNode,
                                             const RealType* gyAccNode,
                                             const RealType* gzAccNode,
                                             const RealType* lumpedMass,
                                             const int* nodeToDof,
                                             const uint8_t* ownership,
                                             RealType* gx, RealType* gy, RealType* gz,
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
// D scatter: integrated SCS-face flux. Same kernel as mars_amr_div.cu's
// computeDivergencePerNodeKernel. The result divAccNode is un-normalized
// (V_i-scaled) at every owned node after reverseExchangeNodeHaloAdd. We keep
// it un-normalized: A = D M^{-1} D^T leaves a V_i factor on every row, which
// is the shape the NS driver consumes (b = -coef * divAccNode matches it).
// =============================================================================

template<typename KeyType, typename RealType>
__global__ void computeDivergencePerNodeKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
    const RealType* vx, const RealType* vy, const RealType* vz,
    const RealType* areaVecX, const RealType* areaVecY, const RealType* areaVecZ,
    RealType* divAccNode, size_t startElem, size_t numLocal)
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
// Driver state. Built once at setup; the apply uses these device buffers.
// =============================================================================

template<typename KeyType, typename RealType>
struct DDTStepper
{
    using DomainT = ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>;
    DomainT& domain;
    int blockSize;
    int rank;
    StretchKind stretch = StretchKind::None;

    size_t nodeCount    = 0;
    size_t elementCount = 0;
    int numOwnedDofs    = 0;

    cstone::DeviceVector<int>      d_node_to_dof;
    cstone::DeviceVector<RealType> d_mass;
    cstone::DeviceVector<RealType> d_areaVec_x;
    cstone::DeviceVector<RealType> d_areaVec_y;
    cstone::DeviceVector<RealType> d_areaVec_z;
};

template<typename KeyType, typename RealType>
void setupStepper(DDTStepper<KeyType, RealType>& s)
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

    // Optional coordinate stretching, applied to ALL nodes (owned + ghosts)
    // before any face geometry is computed. cstone's partition stays as-is;
    // we just relabel cells with non-uniform sizes -- enough to demonstrate
    // that diag(A) varies and Jacobi PCG actually helps.
    if (s.stretch != StretchKind::None)
    {
        const auto& box = s.domain.getBoundingBox();
        RealType xmin = static_cast<RealType>(box.xmin());
        RealType xmax = static_cast<RealType>(box.xmax());
        int kindInt = static_cast<int>(s.stretch);
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        // d_x is a const ref to a DeviceVector whose buffer is mutable.
        RealType* nodeXPtr = const_cast<RealType*>(d_x.data());
        applyCoordStretchKernel<RealType><<<nBlocks, s.blockSize>>>(
            nodeXPtr, s.nodeCount, kindInt, xmin, xmax);
        cudaDeviceSynchronize();
    }
    pt.lap("lazy domain prep");

    s.d_node_to_dof.resize(s.nodeCount);
    s.numOwnedDofs = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(),
                                                 s.d_node_to_dof.data(), s.nodeCount);
    pt.lap("DOF mapping");

    auto c0 = std::get<0>(d_conn).data(); auto c1 = std::get<1>(d_conn).data();
    auto c2 = std::get<2>(d_conn).data(); auto c3 = std::get<3>(d_conn).data();
    auto c4 = std::get<4>(d_conn).data(); auto c5 = std::get<5>(d_conn).data();
    auto c6 = std::get<6>(d_conn).data(); auto c7 = std::get<7>(d_conn).data();

    // Same area-vector source as the implicit assembler, divergence, gradient
    // -- single source of truth for face geometry.
    s.d_areaVec_x.resize(s.elementCount * 12);
    s.d_areaVec_y.resize(s.elementCount * 12);
    s.d_areaVec_z.resize(s.elementCount * 12);
    precomputeAreaVectorsGpu<KeyType, RealType>(
        c0, c1, c2, c3, c4, c5, c6, c7, s.elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data());
    pt.lap("area vectors");

    // Per-node + reverse-halo lumped mass (matches mars_amr_div.cu).
    s.d_mass.resize(s.numOwnedDofs);
    thrust::fill(thrust::device_pointer_cast(s.d_mass.data()),
                 thrust::device_pointer_cast(s.d_mass.data() + s.numOwnedDofs), RealType(0));
    {
        cstone::DeviceVector<RealType> d_massNode(s.nodeCount, RealType(0));
        size_t startElem = s.domain.startIndex();
        size_t numLocal  = s.domain.localElementCount();
        if (numLocal > 0)
        {
            int eBlocks = int((numLocal + s.blockSize - 1) / s.blockSize);
            computeLumpedMassPerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7,
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

// =============================================================================
// A = D M^{-1} D^T applied to phi. phi must have valid ghost values at entry
// (caller exchanges). Output outAcc is un-normalized (V_i-scaled), per-node,
// reverse-halo-closed. gxAcc/gyAcc/gzAcc are scratch buffers (size nodeCount).
// =============================================================================
template<typename KeyType, typename RealType>
void applyDDT(DDTStepper<KeyType, RealType>& s,
              const cstone::DeviceVector<RealType>& phi,
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

    // Step b: g <- M^{-1} g (in-place; gather is same-index, safe).
    normalizeByMassPerNodeKernel<RealType><<<nodeBlocks, s.blockSize>>>(
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
}

// =============================================================================
// Helpers: owned-DOF dot product over per-node arrays, owned-DOF axpy.
// =============================================================================
template<typename RealType>
RealType ownedDot(const cstone::DeviceVector<RealType>& a,
                  const cstone::DeviceVector<RealType>& b,
                  const int* d_nodeToDof,
                  const uint8_t* d_ownership,
                  size_t numNodes)
{
    const RealType* aPtr = a.data();
    const RealType* bPtr = b.data();
    RealType localSum = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(numNodes),
        [aPtr, bPtr, d_nodeToDof, d_ownership] __device__ (size_t i) -> RealType {
            if (d_ownership[i] != 1 || d_nodeToDof[i] < 0) return RealType(0);
            return aPtr[i] * bPtr[i];
        }, RealType(0), thrust::plus<RealType>());
    RealType globalSum = 0;
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&localSum, &globalSum, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    return globalSum;
}

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

// Per-node sum over owned slots only (used for the "remove mean" projector in
// Test 4 and to verify Test 1 / Test 3).
template<typename RealType>
RealType ownedSum(const cstone::DeviceVector<RealType>& a,
                  const int* d_nodeToDof,
                  const uint8_t* d_ownership,
                  size_t numNodes)
{
    const RealType* aPtr = a.data();
    RealType localSum = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(numNodes),
        [aPtr, d_nodeToDof, d_ownership] __device__ (size_t i) -> RealType {
            if (d_ownership[i] != 1 || d_nodeToDof[i] < 0) return RealType(0);
            return aPtr[i];
        }, RealType(0), thrust::plus<RealType>());
    RealType globalSum = 0;
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&localSum, &globalSum, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    return globalSum;
}

long long ownedCount(const int* d_nodeToDof, const uint8_t* d_ownership, size_t numNodes)
{
    long long local = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(numNodes),
        [d_nodeToDof, d_ownership] __device__ (size_t i) -> long long {
            if (d_ownership[i] != 1 || d_nodeToDof[i] < 0) return 0;
            return 1;
        }, 0LL, thrust::plus<long long>());
    long long global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

// =============================================================================
// Boundary mask: unit cube [0,1]^3 with eps=1e-3. Same convention as mars_amr_div.
// =============================================================================
template<typename RealType>
__global__ void buildBoundaryMaskKernel(const RealType* nodeX,
                                        const RealType* nodeY,
                                        const RealType* nodeZ,
                                        uint8_t* boundaryNode,
                                        size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    constexpr RealType eps = RealType(1e-3);
    RealType x = nodeX[i], y = nodeY[i], z = nodeZ[i];
    bool isBnd = (fabs(x) < eps) || (fabs(x - RealType(1)) < eps) ||
                 (fabs(y) < eps) || (fabs(y - RealType(1)) < eps) ||
                 (fabs(z) < eps) || (fabs(z - RealType(1)) < eps);
    boundaryNode[i] = isBnd ? uint8_t(1) : uint8_t(0);
}

// =============================================================================
// Test field generators.
// =============================================================================

template<typename RealType>
__global__ void setLinearFieldKernel(const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
                                     RealType* phi, size_t numNodes,
                                     RealType a, RealType b, RealType c)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    phi[i] = a * nodeX[i] + b * nodeY[i] + c * nodeZ[i];
}

template<typename RealType>
__global__ void setConstantFieldKernel(RealType* phi, size_t numNodes, RealType val)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    phi[i] = val;
}

// Random-but-reproducible per-node field: uses a hash of (seed, rank, nodeId)
// so each rank produces independent values on its OWN nodes, then a halo
// exchange propagates owner values to ghost slots so A reads consistent phi.
template<typename RealType>
__global__ void setRandomFieldKernel(RealType* phi, size_t numNodes,
                                     const int* nodeToDof, const uint8_t* ownership,
                                     unsigned seed, int rank)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1 || nodeToDof[i] < 0) { phi[i] = RealType(0); return; }
    int dof = nodeToDof[i];
    // Cheap deterministic mix: 32-bit splittable hash on (seed, rank, dof).
    unsigned h = seed;
    h ^= (unsigned)(rank * 0x9E3779B1u);
    h ^= (unsigned)(dof  * 0x85EBCA77u);
    h = (h ^ (h >> 16)) * 0x7FEB352Du;
    h = (h ^ (h >> 15)) * 0x846CA68Bu;
    h =  h ^ (h >> 16);
    // Map to (-1, 1).
    RealType u = RealType(h) * RealType(1.0 / 4294967295.0);
    phi[i] = RealType(2) * u - RealType(1);
}

template<typename RealType>
__global__ void setSinSinSinKernel(const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
                                   RealType* phi, size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    constexpr RealType PI = RealType(3.14159265358979323846);
    phi[i] = sin(PI * nodeX[i]) * sin(PI * nodeY[i]) * sin(PI * nodeZ[i]);
}

// Quadratic phi = coord^2 on a chosen axis. Continuum Laplacian is +2
// everywhere, so the V-scaled discrete operator should produce |A*phi|/V_i ~
// +/-2 on every interior node. Used by runSignTest to fix the global sign of
// A; running it on x, y, z separately also catches per-axis sign anisotropy.
template<typename RealType>
__global__ void setQuadraticAxisFieldKernel(const RealType* nodeAxis, RealType* phi, size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    phi[i] = nodeAxis[i] * nodeAxis[i];
}

// Subtract a scalar from owned slots only.
template<typename RealType>
__global__ void subtractScalarOwnedKernel(RealType* phi, RealType c,
                                          const int* nodeToDof, const uint8_t* ownership,
                                          size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1 || nodeToDof[i] < 0) return;
    phi[i] -= c;
}

// =============================================================================
// Test 1: linear field. A*phi should be roundoff on every INTERIOR owned DOF.
// Reports max |Aphi| interior, max |Aphi| boundary, |Aphi|_L2 interior.
// =============================================================================

template<typename KeyType, typename RealType>
void runLinearTest(DDTStepper<KeyType, RealType>& s, RealType a, RealType b, RealType c)
{
    int rank = s.rank;
    if (rank == 0)
    {
        std::cout << "\n----- Test 1: linear field A*phi -----\n";
        std::cout << "  phi = " << a << "*x + " << b << "*y + " << c << "*z\n";
        std::cout << "  Expected: |A*phi| interior ~ roundoff; D^T phi -> (a,b,c)*V; D constant -> 0 interior.\n";
    }

    cstone::DeviceVector<RealType> phi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Aphi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));

    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();
    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    setLinearFieldKernel<RealType><<<nBlocks, s.blockSize>>>(
        d_x.data(), d_y.data(), d_z.data(), phi.data(), s.nodeCount, a, b, c);
    cudaDeviceSynchronize();
    // phi is a closed-form linear function of coords, so ghost values from the
    // kernel are already correct -- but exchange anyway for paranoia and to
    // match the apply-phase contract (phi must be ghost-consistent on entry).
    s.domain.exchangeNodeHalo(phi);

    applyDDT(s, phi, Aphi, gx, gy, gz);

    // Build per-node boundary mask once.
    cstone::DeviceVector<uint8_t> d_bnd(s.nodeCount, uint8_t(0));
    buildBoundaryMaskKernel<RealType><<<nBlocks, s.blockSize>>>(
        d_x.data(), d_y.data(), d_z.data(), d_bnd.data(), s.nodeCount);
    cudaDeviceSynchronize();

    const auto& d_own = s.domain.getNodeOwnershipMap();
    const uint8_t* ownPtr = d_own.data();
    const uint8_t* bndPtr = d_bnd.data();
    const int*     dofPtr = s.d_node_to_dof.data();
    const RealType* Aptr  = Aphi.data();

    // L2 (sum of squares) and max over interior/boundary owned slots.
    auto cBegin = thrust::counting_iterator<size_t>(0);
    auto cEnd   = thrust::counting_iterator<size_t>(s.nodeCount);

    RealType locL2sqInt = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, bndPtr, Aptr] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1 || dofPtr[i] < 0 || bndPtr[i] != 0) return RealType(0);
            return Aptr[i] * Aptr[i];
        }, RealType(0), thrust::plus<RealType>());
    RealType locMaxInt = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, bndPtr, Aptr] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1 || dofPtr[i] < 0 || bndPtr[i] != 0) return RealType(0);
            return fabs(Aptr[i]);
        }, RealType(0), thrust::maximum<RealType>());
    RealType locMaxBnd = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, bndPtr, Aptr] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1 || dofPtr[i] < 0 || bndPtr[i] == 0) return RealType(0);
            return fabs(Aptr[i]);
        }, RealType(0), thrust::maximum<RealType>());

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    RealType gL2sqInt = 0, gMaxInt = 0, gMaxBnd = 0;
    MPI_Allreduce(&locL2sqInt, &gL2sqInt, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locMaxInt,  &gMaxInt,  1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locMaxBnd,  &gMaxBnd,  1, mpiType, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << std::scientific << std::setprecision(3)
                  << "  |A*phi|_L2 interior:  " << std::sqrt(gL2sqInt) << "\n"
                  << "  max |A*phi| interior: " << gMaxInt << "\n"
                  << "  max |A*phi| boundary: " << gMaxBnd << "  (expected non-zero: boundary CV captures through-flow)\n";
        bool pass = (double(gMaxInt) < 1e-10);
        std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL") << " (interior max < 1e-10 in fp64)\n"
                  << std::defaultfloat;
    }
}

// =============================================================================
// Sign-discrimination test: phi = coord^2 has continuum Laplacian +2 everywhere.
// The V-scaled discrete operator (D M^{-1} D^T) phi should therefore produce
// either +2*V_i  (A models  +Laplacian, "positive Poisson")
// or     -2*V_i  (A models  -Laplacian, "standard SPD Poisson, like K-path")
// on every INTERIOR owned node.
//
// Why we need this: the PSD test only sees phi^T A phi > 0, which is true for
// both A and -A applied with a flipped CG operator -- it does not pin the sign
// of A relative to the corrector's grad(phi). coord^2 has a known SIGNED
// Laplacian so it isolates exactly that. We run on x, y, z separately so any
// per-axis sign anisotropy in A's D^T or D scatter would surface.
// =============================================================================

// Per-axis result block: signed extrema, mean, RMS deviation from +/-2.
template<typename RealType>
struct SignAxisResult
{
    long long cnt = 0;
    RealType minVal = 0;
    RealType maxVal = 0;
    RealType meanVal = 0;
    RealType rmsPos  = 0; // deviation from +2
    RealType rmsNeg  = 0; // deviation from -2
};

// Reduce one axis's (A*phi)/V over interior owned nodes against +/-2 targets.
template<typename KeyType, typename RealType>
SignAxisResult<RealType> runSignAxis(DDTStepper<KeyType, RealType>& s,
                                     const RealType* d_axis,
                                     const uint8_t* bndPtr)
{
    cstone::DeviceVector<RealType> phi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Aphi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));

    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    setQuadraticAxisFieldKernel<RealType><<<nBlocks, s.blockSize>>>(
        d_axis, phi.data(), s.nodeCount);
    cudaDeviceSynchronize();
    // Closed-form on coords, but honour the apply-phase contract.
    s.domain.exchangeNodeHalo(phi);

    applyDDT(s, phi, Aphi, gx, gy, gz);

    const auto& d_own = s.domain.getNodeOwnershipMap();
    const uint8_t* ownPtr  = d_own.data();
    const int*     dofPtr  = s.d_node_to_dof.data();
    const RealType* AphiPtr = Aphi.data();
    const RealType* massPtr = s.d_mass.data();

    auto cBegin = thrust::counting_iterator<size_t>(0);
    auto cEnd   = thrust::counting_iterator<size_t>(s.nodeCount);

    RealType locSum = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, bndPtr, AphiPtr, massPtr] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1 || dofPtr[i] < 0 || bndPtr[i] != 0) return RealType(0);
            return AphiPtr[i] / massPtr[dofPtr[i]];
        }, RealType(0), thrust::plus<RealType>());
    long long locCnt = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, bndPtr] __device__ (size_t i) -> long long {
            if (ownPtr[i] != 1 || dofPtr[i] < 0 || bndPtr[i] != 0) return 0LL;
            return 1LL;
        }, 0LL, thrust::plus<long long>());
    RealType locMin = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, bndPtr, AphiPtr, massPtr] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1 || dofPtr[i] < 0 || bndPtr[i] != 0) return RealType( 1e300);
            return AphiPtr[i] / massPtr[dofPtr[i]];
        }, RealType(1e300), thrust::minimum<RealType>());
    RealType locMax = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, bndPtr, AphiPtr, massPtr] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1 || dofPtr[i] < 0 || bndPtr[i] != 0) return RealType(-1e300);
            return AphiPtr[i] / massPtr[dofPtr[i]];
        }, RealType(-1e300), thrust::maximum<RealType>());
    RealType locL2sqPos = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, bndPtr, AphiPtr, massPtr] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1 || dofPtr[i] < 0 || bndPtr[i] != 0) return RealType(0);
            RealType v = AphiPtr[i] / massPtr[dofPtr[i]] - RealType(2);
            return v * v;
        }, RealType(0), thrust::plus<RealType>());
    RealType locL2sqNeg = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, bndPtr, AphiPtr, massPtr] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1 || dofPtr[i] < 0 || bndPtr[i] != 0) return RealType(0);
            RealType v = AphiPtr[i] / massPtr[dofPtr[i]] + RealType(2);
            return v * v;
        }, RealType(0), thrust::plus<RealType>());

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    SignAxisResult<RealType> r;
    RealType gSum = 0, gL2sqPos = 0, gL2sqNeg = 0;
    MPI_Allreduce(&locSum,     &gSum,     1, mpiType,        MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locMin,     &r.minVal, 1, mpiType,        MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&locMax,     &r.maxVal, 1, mpiType,        MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locCnt,     &r.cnt,    1, MPI_LONG_LONG,  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locL2sqPos, &gL2sqPos, 1, mpiType,        MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locL2sqNeg, &gL2sqNeg, 1, mpiType,        MPI_SUM, MPI_COMM_WORLD);
    r.meanVal = (r.cnt > 0) ? (gSum / RealType(r.cnt)) : RealType(0);
    r.rmsPos  = (r.cnt > 0) ? std::sqrt(gL2sqPos / RealType(r.cnt)) : RealType(0);
    r.rmsNeg  = (r.cnt > 0) ? std::sqrt(gL2sqNeg / RealType(r.cnt)) : RealType(0);
    return r;
}

template<typename KeyType, typename RealType>
void runSignTest(DDTStepper<KeyType, RealType>& s)
{
    int rank = s.rank;
    if (rank == 0)
    {
        std::cout << "\n----- Test 1b: sign-discrimination phi = x^2, y^2, z^2 -----\n";
        std::cout << "  Continuum Laplacian = +2 everywhere; expected (A*phi)/V_i = +/- 2 interior.\n";
    }

    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();
    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;

    cstone::DeviceVector<uint8_t> d_bnd(s.nodeCount, uint8_t(0));
    buildBoundaryMaskKernel<RealType><<<nBlocks, s.blockSize>>>(
        d_x.data(), d_y.data(), d_z.data(), d_bnd.data(), s.nodeCount);
    cudaDeviceSynchronize();

    auto rx = runSignAxis<KeyType, RealType>(s, d_x.data(), d_bnd.data());
    auto ry = runSignAxis<KeyType, RealType>(s, d_y.data(), d_bnd.data());
    auto rz = runSignAxis<KeyType, RealType>(s, d_z.data(), d_bnd.data());

    auto axisSign = [] (const SignAxisResult<RealType>& r) -> int {
        return (r.rmsPos < r.rmsNeg) ? +1 : -1;
    };
    auto axisLabel = [] (int s) { return s > 0 ? "+ (matches +2)" : "- (matches -2)"; };

    if (rank == 0)
    {
        auto report = [&] (const char* name, const SignAxisResult<RealType>& r) {
            std::cout << "  axis " << name << ":\n"
                      << std::scientific << std::setprecision(6)
                      << "    interior DOFs           " << r.cnt    << "\n"
                      << "    (A*phi/V) min interior  " << r.minVal << "\n"
                      << "    (A*phi/V) max interior  " << r.maxVal << "\n"
                      << "    (A*phi/V) mean interior " << r.meanVal << "\n"
                      << "    RMS deviation from +2   " << r.rmsPos  << "\n"
                      << "    RMS deviation from -2   " << r.rmsNeg  << "\n"
                      << "    inferred sign           " << axisLabel(axisSign(r)) << "\n"
                      << std::defaultfloat;
        };
        report("x", rx);
        report("y", ry);
        report("z", rz);

        int sx = axisSign(rx), sy = axisSign(ry), sz = axisSign(rz);
        bool consistent = (sx == sy) && (sy == sz);
        std::cout << "  per-axis signs (x,y,z): "
                  << (sx > 0 ? "+ " : "- ")
                  << (sy > 0 ? "+ " : "- ")
                  << (sz > 0 ? "+"  : "-") << "\n";
        std::cout << "  STATUS: " << (consistent ? "CONSISTENT" : "ANISOTROPIC -- BUG") << "\n";
        if (consistent)
        {
            const char* sign = (sx > 0) ? "+ (A models +Laplacian)"
                                        : "- (A models -Laplacian, K-path convention)";
            std::cout << "  INFERRED GLOBAL SIGN OF A: " << sign << "\n";
            if (sx > 0)
            {
                std::cout << "  IMPLIED NS RHS SIGN: b = +(rho/dt) * divAccNode\n"
                          << "                       (corrector: u = u** - (dt/rho) * grad(phi))\n";
            }
            else
            {
                std::cout << "  IMPLIED NS RHS SIGN: b = -(rho/dt) * divAccNode\n"
                          << "                       (corrector: u = u** - (dt/rho) * grad(phi))\n";
            }
        }
    }
}

// =============================================================================
// Test 2: symmetry. Random u, v -> u^T (A v) == v^T (A u) to roundoff.
// =============================================================================

template<typename KeyType, typename RealType>
void runSymmetryTest(DDTStepper<KeyType, RealType>& s, unsigned seed)
{
    int rank = s.rank;
    if (rank == 0)
    {
        std::cout << "\n----- Test 2: symmetry u^T A v == v^T A u -----\n";
        std::cout << "  seed=" << seed << "\n";
    }

    cstone::DeviceVector<RealType> u(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> v(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Av(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Au(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));

    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    const auto& d_own = s.domain.getNodeOwnershipMap();
    setRandomFieldKernel<RealType><<<nBlocks, s.blockSize>>>(
        u.data(), s.nodeCount, s.d_node_to_dof.data(), d_own.data(), seed,      rank);
    setRandomFieldKernel<RealType><<<nBlocks, s.blockSize>>>(
        v.data(), s.nodeCount, s.d_node_to_dof.data(), d_own.data(), seed + 1u, rank);
    cudaDeviceSynchronize();
    // Make ghost slots see owner-rank values -- A needs ghost phi to scatter
    // face contributions correctly.
    s.domain.exchangeNodeHalo(u);
    s.domain.exchangeNodeHalo(v);

    applyDDT(s, v, Av, gx, gy, gz);
    applyDDT(s, u, Au, gx, gy, gz);

    RealType s1 = ownedDot(u, Av, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType s2 = ownedDot(v, Au, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);

    if (rank == 0)
    {
        RealType diff = std::fabs(s1 - s2);
        RealType denom = std::max(std::fabs(s1), std::fabs(s2));
        RealType rel = (denom > RealType(0)) ? diff / denom : diff;
        std::cout << std::scientific << std::setprecision(6)
                  << "  s1 = u^T A v = " << s1 << "\n"
                  << "  s2 = v^T A u = " << s2 << "\n"
                  << "  |s1 - s2|     = " << diff << "\n"
                  << "  rel diff      = " << rel  << "\n";
        bool pass = (double(rel) < 1e-12);
        std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL") << " (rel < 1e-12 in fp64)\n"
                  << std::defaultfloat;
    }
}

// =============================================================================
// Test 3: null space. A*1 = 0 to roundoff.
// =============================================================================

template<typename KeyType, typename RealType>
void runNullSpaceTest(DDTStepper<KeyType, RealType>& s)
{
    int rank = s.rank;
    if (rank == 0)
    {
        std::cout << "\n----- Test 3: null space A*1 = 0 -----\n";
    }

    cstone::DeviceVector<RealType> phi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Aphi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));

    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    setConstantFieldKernel<RealType><<<nBlocks, s.blockSize>>>(phi.data(), s.nodeCount, RealType(1));
    cudaDeviceSynchronize();
    // phi=1 everywhere already so ghost halo exchange is a no-op, but keep
    // the contract consistent.
    s.domain.exchangeNodeHalo(phi);

    applyDDT(s, phi, Aphi, gx, gy, gz);

    const auto& d_own = s.domain.getNodeOwnershipMap();
    const uint8_t* ownPtr = d_own.data();
    const int*     dofPtr = s.d_node_to_dof.data();
    const RealType* Aptr  = Aphi.data();
    auto cBegin = thrust::counting_iterator<size_t>(0);
    auto cEnd   = thrust::counting_iterator<size_t>(s.nodeCount);
    RealType locMax = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, Aptr] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1 || dofPtr[i] < 0) return RealType(0);
            return fabs(Aptr[i]);
        }, RealType(0), thrust::maximum<RealType>());
    RealType locL2sq = thrust::transform_reduce(thrust::device, cBegin, cEnd,
        [ownPtr, dofPtr, Aptr] __device__ (size_t i) -> RealType {
            if (ownPtr[i] != 1 || dofPtr[i] < 0) return RealType(0);
            return Aptr[i] * Aptr[i];
        }, RealType(0), thrust::plus<RealType>());

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    RealType gMax = 0, gL2sq = 0;
    MPI_Allreduce(&locMax,  &gMax,  1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locL2sq, &gL2sq, 1, mpiType, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << std::scientific << std::setprecision(3)
                  << "  max |A*1| (all owned): " << gMax << "\n"
                  << "  |A*1|_L2 (all owned):  " << std::sqrt(gL2sq) << "\n";
        bool pass = (double(gMax) < 1e-12);
        std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL") << " (max < 1e-12 in fp64)\n"
                  << std::defaultfloat;
    }
}

// =============================================================================
// Test 4: positive semi-definiteness on orthogonal-to-constant subspace.
// =============================================================================

template<typename KeyType, typename RealType>
void runPsdTest(DDTStepper<KeyType, RealType>& s, unsigned seed)
{
    int rank = s.rank;
    if (rank == 0)
    {
        std::cout << "\n----- Test 4: phi^T A phi > 0 on null-space-orthogonal phi -----\n";
        std::cout << "  seed=" << seed << "\n";
    }

    cstone::DeviceVector<RealType> phi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Aphi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));

    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    const auto& d_own = s.domain.getNodeOwnershipMap();
    setRandomFieldKernel<RealType><<<nBlocks, s.blockSize>>>(
        phi.data(), s.nodeCount, s.d_node_to_dof.data(), d_own.data(), seed, rank);
    cudaDeviceSynchronize();

    // Project out the constant: phi <- phi - mean(phi).
    RealType sumPhi = ownedSum(phi, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    long long N = ownedCount(s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType meanPhi = (N > 0) ? sumPhi / RealType(N) : RealType(0);
    subtractScalarOwnedKernel<RealType><<<nBlocks, s.blockSize>>>(
        phi.data(), meanPhi, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    cudaDeviceSynchronize();
    s.domain.exchangeNodeHalo(phi);

    applyDDT(s, phi, Aphi, gx, gy, gz);
    RealType q = ownedDot(phi, Aphi, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);

    if (rank == 0)
    {
        std::cout << std::scientific << std::setprecision(6)
                  << "  phi^T A phi = " << q << "\n"
                  << "  mean(phi)   = " << meanPhi << "  (projected to ~0)\n";
        bool pass = (double(q) > 0.0);
        std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL") << " (strictly > 0)\n"
                  << std::defaultfloat;
    }
}

// =============================================================================
// diag(A) by brute-force probe (Stage 2). For each owned DOF j we apply A to
// e_j and read (A e_j)[j]. Every applyDDT is a global collective, so all ranks
// run in lockstep up to max(numOwnedDofs); ranks that ran out of local DOFs
// apply a zero probe. O(N^2) total -- OK for cube16, replace with coloured
// probes (~27 applies on uniform-hex) or an analytic kernel before scaling.
// =============================================================================

template<typename RealType>
__global__ void setProbeUnitKernel(RealType* phi, size_t numNodes,
                                   const int* nodeToDof, const uint8_t* ownership,
                                   int targetDof)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    RealType v = RealType(0);
    if (ownership[i] == 1 && nodeToDof[i] == targetDof) v = RealType(1);
    phi[i] = v;
}

// Pull (A e_j)[j] = Aphi[i] where node i is the owner of DOF j on this rank.
// Returns 0 on ranks that don't own j; caller does a global reduce.
template<typename RealType>
__global__ void extractDiagEntryKernel(const RealType* Aphi, size_t numNodes,
                                       const int* nodeToDof, const uint8_t* ownership,
                                       int targetDof, RealType* outScalar)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] == 1 && nodeToDof[i] == targetDof) *outScalar = Aphi[i];
}

template<typename KeyType, typename RealType>
void computeDiagAByProbe(DDTStepper<KeyType, RealType>& s,
                         cstone::DeviceVector<RealType>& d_diagA)
{
    int rank = s.rank;
    if (rank == 0)
        std::cout << "  Setup: computing diag(A) via probe (this is brute-force "
                     "O(N^2) -- OK for cube16, replace before scaling)\n";

    int locOwn = s.numOwnedDofs;
    int maxOwn = 0;
    MPI_Allreduce(&locOwn, &maxOwn, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    d_diagA.resize(s.numOwnedDofs);
    thrust::fill(thrust::device_pointer_cast(d_diagA.data()),
                 thrust::device_pointer_cast(d_diagA.data() + s.numOwnedDofs),
                 RealType(0));

    cstone::DeviceVector<RealType> phi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Aphi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> d_scalar(1, RealType(0));

    const auto& d_own = s.domain.getNodeOwnershipMap();
    const uint8_t* ownPtr = d_own.data();
    const int*     dofPtr = s.d_node_to_dof.data();
    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;

    // Buffer of per-DOF diagonals back on host -- we write into it as we go and
    // copy back to d_diagA at the end.
    std::vector<RealType> h_diag(s.numOwnedDofs, RealType(0));

    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;

    auto tStart = std::chrono::high_resolution_clock::now();
    int reportEvery = std::max(maxOwn / 20, 1);

    // Each rank picks targetDof from [0, numOwnedDofs) on its local DOF
    // numbering; ranks past their local count send a -1 (no-op probe).
    for (int probe = 0; probe < maxOwn; ++probe)
    {
        int targetDof = (probe < s.numOwnedDofs) ? probe : -1;

        // Set phi to e_{targetDof} on this rank (or 0 if we ran out of DOFs).
        setProbeUnitKernel<RealType><<<nBlocks, s.blockSize>>>(
            phi.data(), s.nodeCount, dofPtr, ownPtr, targetDof);
        cudaDeviceSynchronize();
        s.domain.exchangeNodeHalo(phi);
        applyDDT(s, phi, Aphi, gx, gy, gz);

        if (targetDof >= 0)
        {
            thrust::fill(thrust::device_pointer_cast(d_scalar.data()),
                         thrust::device_pointer_cast(d_scalar.data() + 1), RealType(0));
            extractDiagEntryKernel<RealType><<<nBlocks, s.blockSize>>>(
                Aphi.data(), s.nodeCount, dofPtr, ownPtr, targetDof, d_scalar.data());
            cudaDeviceSynchronize();
            RealType v = 0;
            cudaMemcpy(&v, d_scalar.data(), sizeof(RealType), cudaMemcpyDeviceToHost);
            h_diag[targetDof] = v;
        }

        if (rank == 0 && ((probe + 1) % reportEvery == 0 || probe + 1 == maxOwn))
        {
            auto tNow = std::chrono::high_resolution_clock::now();
            float ms = std::chrono::duration<float, std::milli>(tNow - tStart).count();
            std::cout << "    probe " << std::setw(6) << probe + 1 << " / " << maxOwn
                      << "  elapsed " << std::fixed << std::setprecision(0) << ms
                      << " ms\n" << std::defaultfloat;
        }
    }

    cudaMemcpy(d_diagA.data(), h_diag.data(),
               sizeof(RealType) * s.numOwnedDofs, cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();

    // Sanity: count non-positive diagonals (would flag an operator bug).
    long long locNeg = 0;
    RealType locMin = std::numeric_limits<RealType>::infinity(), locMax = 0;
    for (int j = 0; j < s.numOwnedDofs; ++j)
    {
        if (h_diag[j] <= RealType(0)) ++locNeg;
        if (h_diag[j] < locMin) locMin = h_diag[j];
        if (h_diag[j] > locMax) locMax = h_diag[j];
    }
    long long gNeg = 0;
    RealType gMin = 0, gMax = 0;
    MPI_Allreduce(&locNeg, &gNeg, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&locMin, &gMin, 1, mpiType, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&locMax, &gMax, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0)
    {
        auto tEnd = std::chrono::high_resolution_clock::now();
        float ms = std::chrono::duration<float, std::milli>(tEnd - tStart).count();
        std::cout << std::scientific << std::setprecision(3)
                  << "  diag(A) min = " << gMin << "  max = " << gMax
                  << "  ratio = " << (gMin > 0 ? gMax / gMin : RealType(0))
                  << "  non-positive entries (global): " << gNeg << "\n"
                  << std::fixed << std::setprecision(0)
                  << "  diag(A) build wall time: " << ms << " ms\n"
                  << std::defaultfloat;
    }
}

// Elementwise z = r / diag(A) on owned DOFs; ghost / non-owned set to 0.
template<typename RealType>
__global__ void jacobiPrecondKernel(const RealType* r, const RealType* diagA,
                                    const int* nodeToDof, const uint8_t* ownership,
                                    RealType* z, size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1 || nodeToDof[i] < 0) { z[i] = RealType(0); return; }
    int dof = nodeToDof[i];
    RealType d = diagA[dof];
    // Diag must be positive for SPSD A; guard against the pathological 0 case.
    z[i] = (d > RealType(0)) ? r[i] / d : r[i];
}

// =============================================================================
// Test 5: un-preconditioned matrix-free CG solving A*phi = b_ref where
// b_ref = A * phi_ref and phi_ref = sin(pi*x)*sin(pi*y)*sin(pi*z).
// phi_ref vanishes on the cube boundary and integrates to zero on the cube, so
// it is null-space-orthogonal and b_ref is in the range of A by construction
// -- no pin needed.
// =============================================================================

template<typename KeyType, typename RealType>
void runCgTest(DDTStepper<KeyType, RealType>& s, int maxIter, RealType tol)
{
    int rank = s.rank;
    if (rank == 0)
    {
        std::cout << "\n----- Test 5a: matrix-free CG (un-preconditioned), A phi = b_ref -----\n";
        std::cout << "  phi_ref = sin(pi*x)*sin(pi*y)*sin(pi*z); b_ref = A*phi_ref\n";
        std::cout << "  max iters=" << maxIter << ", tol=" << tol << "\n";
    }

    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();
    const auto& d_own = s.domain.getNodeOwnershipMap();

    cstone::DeviceVector<RealType> phiRef(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> b(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> phi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> r(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> p(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Ap(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));

    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    setSinSinSinKernel<RealType><<<nBlocks, s.blockSize>>>(
        d_x.data(), d_y.data(), d_z.data(), phiRef.data(), s.nodeCount);
    cudaDeviceSynchronize();
    s.domain.exchangeNodeHalo(phiRef);

    // b = A * phi_ref. Un-normalized accumulator form -- the same shape the
    // CG operates on throughout (Ap, r, p, phi are all per-node accumulators).
    applyDDT(s, phiRef, b, gx, gy, gz);

    // CG. phi = 0 initial guess so r0 = b - A*0 = b.
    thrust::fill(thrust::device_pointer_cast(phi.data()),
                 thrust::device_pointer_cast(phi.data() + s.nodeCount), RealType(0));
    thrust::copy(thrust::device_pointer_cast(b.data()),
                 thrust::device_pointer_cast(b.data() + s.nodeCount),
                 thrust::device_pointer_cast(r.data()));
    thrust::copy(thrust::device_pointer_cast(r.data()),
                 thrust::device_pointer_cast(r.data() + s.nodeCount),
                 thrust::device_pointer_cast(p.data()));

    RealType rho_old = ownedDot(r, r, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType r0_norm = std::sqrt(rho_old);
    if (r0_norm < std::numeric_limits<RealType>::min())
    {
        if (rank == 0) std::cout << "  trivial: |b| = 0; nothing to solve\n";
        return;
    }
    const RealType absTol = tol * r0_norm;

    if (rank == 0)
    {
        std::cout << std::scientific << std::setprecision(3)
                  << "  |r0| = |b| = " << r0_norm << "\n"
                  << "  absTol = " << absTol << "\n"
                  << std::defaultfloat;
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    int iters = -2;
    int nodeBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    for (int it = 0; it < maxIter; ++it)
    {
        // Ghost slots of p must be valid: A's D^T reads p[L], p[R] at every face.
        s.domain.exchangeNodeHalo(p);
        applyDDT(s, p, Ap, gx, gy, gz);

        RealType pAp = ownedDot(p, Ap, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        if (pAp <= RealType(0))
        {
            if (rank == 0)
                std::cout << "  WARNING: pAp = " << pAp << " <= 0 at iter " << it
                          << " -- A appears non-PSD on residual subspace (BUG)\n";
            iters = -2;
            break;
        }
        RealType alpha = rho_old / pAp;

        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            phi.data(), phi.data(), p.data(), alpha,
            s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            r.data(), r.data(), Ap.data(), -alpha,
            s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        cudaDeviceSynchronize();

        RealType rho_new = ownedDot(r, r, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        if (std::sqrt(rho_new) < absTol) { iters = it + 1; break; }
        RealType beta = rho_new / rho_old;
        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            p.data(), r.data(), p.data(), beta,
            s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        cudaDeviceSynchronize();
        rho_old = rho_new;

        if (rank == 0 && (it < 5 || (it % 50) == 0))
        {
            std::cout << "    iter " << std::setw(5) << it + 1
                      << "  |r| = " << std::scientific << std::setprecision(3)
                      << std::sqrt(rho_new) << std::defaultfloat << "\n";
        }
    }

    cudaDeviceSynchronize();
    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = std::chrono::high_resolution_clock::now();
    float ms = std::chrono::duration<float, std::milli>(t1 - t0).count();

    if (iters < 0) iters = maxIter;
    s.domain.exchangeNodeHalo(phi);

    // CG solves A phi = b up to a constant (A has a 1D null space). Compare
    // (phi - phi_ref) AFTER projecting out the constant from both (the constant
    // offset between them is undetermined by the equation).
    cstone::DeviceVector<RealType> diff(s.nodeCount, RealType(0));
    {
        const RealType* phPtr  = phi.data();
        const RealType* refPtr = phiRef.data();
        RealType* dPtr         = diff.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [phPtr, refPtr, dPtr] __device__ (size_t i) {
                dPtr[i] = phPtr[i] - refPtr[i];
            });
        cudaDeviceSynchronize();
    }
    RealType sumDiff = ownedSum(diff, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    long long N = ownedCount(s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType meanDiff = (N > 0) ? sumDiff / RealType(N) : RealType(0);
    subtractScalarOwnedKernel<RealType><<<nBlocks, s.blockSize>>>(
        diff.data(), meanDiff, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    cudaDeviceSynchronize();

    // |phi - phi_ref - mean|_L2 over owned, both raw and relative to |phi_ref|_L2.
    RealType errL2sq  = ownedDot(diff, diff, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType refL2sq  = ownedDot(phiRef, phiRef, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType locMaxErr = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount),
        [d_own_ptr = d_own.data(), dofPtr = s.d_node_to_dof.data(), dPtr = diff.data()] __device__ (size_t i) -> RealType {
            if (d_own_ptr[i] != 1 || dofPtr[i] < 0) return RealType(0);
            return fabs(dPtr[i]);
        }, RealType(0), thrust::maximum<RealType>());
    RealType gMaxErr = 0;
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&locMaxErr, &gMaxErr, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0)
    {
        RealType errL2 = std::sqrt(errL2sq);
        RealType refL2 = std::sqrt(refL2sq);
        RealType rel   = (refL2 > RealType(0)) ? errL2 / refL2 : errL2;
        std::cout << std::fixed << std::setprecision(2) << "  CG time: " << ms << " ms\n"
                  << "  iters:   " << iters
                  << (iters >= maxIter ? "  (MAX HIT)" : "") << "\n";
        std::cout << std::scientific << std::setprecision(3)
                  << "  |phi - phi_ref - mean|_L2:  " << errL2 << "\n"
                  << "  |phi_ref|_L2:               " << refL2 << "\n"
                  << "  relative L2 error:          " << rel   << "\n"
                  << "  max |phi - phi_ref - mean|: " << gMaxErr << "\n";
        // Stage-1 expectation: relative error in the O(h^2) range. For cube16
        // (h ~ 1/16) on the unit cube that's a few percent; CG should converge
        // well below 1% in a few hundred iterations as long as the operator is
        // correct. Use 5e-2 as a loose Stage-1 gate -- accuracy is dominated by
        // the discretization, not the solver.
        bool pass = (iters > 0 && iters < maxIter && double(rel) < 5e-2);
        std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL") << " (converged + rel < 5e-2)\n"
                  << std::defaultfloat;
    }
}

// =============================================================================
// Test 5b: Jacobi-PCG. Same RHS as Test 5a (b = A * sin*sin*sin) so iteration
// counts are directly comparable. Expected to drop from ~44 to ~10-20 on cube16.
// =============================================================================
template<typename KeyType, typename RealType>
void runCgJacobiTest(DDTStepper<KeyType, RealType>& s, int maxIter, RealType tol)
{
    int rank = s.rank;
    if (rank == 0)
        std::cout << "\n----- Test 5b: matrix-free PCG (Jacobi), A phi = b_ref -----\n"
                  << "  phi_ref = sin(pi*x)*sin(pi*y)*sin(pi*z); b_ref = A*phi_ref\n"
                  << "  M_jac = diag(A); max iters=" << maxIter << ", tol=" << tol << "\n";

    const auto& d_x = s.domain.getNodeX();
    const auto& d_y = s.domain.getNodeY();
    const auto& d_z = s.domain.getNodeZ();
    const auto& d_own = s.domain.getNodeOwnershipMap();

    cstone::DeviceVector<RealType> diagA;
    computeDiagAByProbe(s, diagA);

    cstone::DeviceVector<RealType> phiRef(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> b(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> phi(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> r(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> z(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> p(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> Ap(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gx(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gy(s.nodeCount, RealType(0));
    cstone::DeviceVector<RealType> gz(s.nodeCount, RealType(0));

    int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
    int nodeBlocks = nBlocks;

    setSinSinSinKernel<RealType><<<nBlocks, s.blockSize>>>(
        d_x.data(), d_y.data(), d_z.data(), phiRef.data(), s.nodeCount);
    cudaDeviceSynchronize();
    s.domain.exchangeNodeHalo(phiRef);
    applyDDT(s, phiRef, b, gx, gy, gz);

    // phi=0 initial guess => r0 = b; z0 = M_jac^{-1} r0; p0 = z0.
    thrust::fill(thrust::device_pointer_cast(phi.data()),
                 thrust::device_pointer_cast(phi.data() + s.nodeCount), RealType(0));
    thrust::copy(thrust::device_pointer_cast(b.data()),
                 thrust::device_pointer_cast(b.data() + s.nodeCount),
                 thrust::device_pointer_cast(r.data()));
    jacobiPrecondKernel<RealType><<<nodeBlocks, s.blockSize>>>(
        r.data(), diagA.data(), s.d_node_to_dof.data(), d_own.data(),
        z.data(), s.nodeCount);
    cudaDeviceSynchronize();
    thrust::copy(thrust::device_pointer_cast(z.data()),
                 thrust::device_pointer_cast(z.data() + s.nodeCount),
                 thrust::device_pointer_cast(p.data()));

    // Convergence check stays on the unpreconditioned |r| -- directly comparable to Test 5a.
    RealType rho_old = ownedDot(r, z, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType r0sq    = ownedDot(r, r, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType r0_norm = std::sqrt(r0sq);
    if (r0_norm < std::numeric_limits<RealType>::min())
    {
        if (rank == 0) std::cout << "  trivial: |b| = 0; nothing to solve\n";
        return;
    }
    const RealType absTol = tol * r0_norm;

    if (rank == 0)
        std::cout << std::scientific << std::setprecision(3)
                  << "  |r0| = |b| = " << r0_norm << "  <r0,z0> = " << rho_old
                  << "  absTol = " << absTol << "\n" << std::defaultfloat;

    auto t0 = std::chrono::high_resolution_clock::now();

    int iters = -2;
    for (int it = 0; it < maxIter; ++it)
    {
        s.domain.exchangeNodeHalo(p);
        applyDDT(s, p, Ap, gx, gy, gz);

        RealType pAp = ownedDot(p, Ap, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        if (pAp <= RealType(0))
        {
            if (rank == 0)
                std::cout << "  WARNING: pAp = " << pAp << " <= 0 at iter " << it
                          << " -- A appears non-PSD on residual subspace (BUG)\n";
            iters = -2;
            break;
        }
        RealType alpha = rho_old / pAp;

        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            phi.data(), phi.data(), p.data(), alpha,
            s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            r.data(), r.data(), Ap.data(), -alpha,
            s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        cudaDeviceSynchronize();

        RealType rNormSq = ownedDot(r, r, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        if (std::sqrt(rNormSq) < absTol) { iters = it + 1; break; }

        jacobiPrecondKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            r.data(), diagA.data(), s.d_node_to_dof.data(), d_own.data(),
            z.data(), s.nodeCount);
        cudaDeviceSynchronize();

        RealType rho_new = ownedDot(r, z, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        RealType beta    = rho_new / rho_old;
        axpyOwnedKernel<RealType><<<nodeBlocks, s.blockSize>>>(
            p.data(), z.data(), p.data(), beta,
            s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
        cudaDeviceSynchronize();
        rho_old = rho_new;

        if (rank == 0 && (it < 5 || (it % 10) == 0))
        {
            std::cout << "    iter " << std::setw(5) << it + 1
                      << "  |r| = " << std::scientific << std::setprecision(3)
                      << std::sqrt(rNormSq) << std::defaultfloat << "\n";
        }
    }

    cudaDeviceSynchronize();
    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = std::chrono::high_resolution_clock::now();
    float ms = std::chrono::duration<float, std::milli>(t1 - t0).count();

    if (iters < 0) iters = maxIter;
    s.domain.exchangeNodeHalo(phi);

    // Compare (phi - phi_ref) AFTER projecting out the constant (1D null space).
    cstone::DeviceVector<RealType> diff(s.nodeCount, RealType(0));
    {
        const RealType* phPtr  = phi.data();
        const RealType* refPtr = phiRef.data();
        RealType* dPtr         = diff.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [phPtr, refPtr, dPtr] __device__ (size_t i) {
                dPtr[i] = phPtr[i] - refPtr[i];
            });
        cudaDeviceSynchronize();
    }
    RealType sumDiff = ownedSum(diff, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    long long N = ownedCount(s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType meanDiff = (N > 0) ? sumDiff / RealType(N) : RealType(0);
    subtractScalarOwnedKernel<RealType><<<nBlocks, s.blockSize>>>(
        diff.data(), meanDiff, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    cudaDeviceSynchronize();

    RealType errL2sq  = ownedDot(diff, diff, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType refL2sq  = ownedDot(phiRef, phiRef, s.d_node_to_dof.data(), d_own.data(), s.nodeCount);
    RealType locMaxErr = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0),
        thrust::counting_iterator<size_t>(s.nodeCount),
        [d_own_ptr = d_own.data(), dofPtr = s.d_node_to_dof.data(), dPtr = diff.data()] __device__ (size_t i) -> RealType {
            if (d_own_ptr[i] != 1 || dofPtr[i] < 0) return RealType(0);
            return fabs(dPtr[i]);
        }, RealType(0), thrust::maximum<RealType>());
    RealType gMaxErr = 0;
    auto mpiType = std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&locMaxErr, &gMaxErr, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0)
    {
        RealType errL2 = std::sqrt(errL2sq);
        RealType refL2 = std::sqrt(refL2sq);
        RealType rel   = (refL2 > RealType(0)) ? errL2 / refL2 : errL2;
        std::cout << std::fixed << std::setprecision(2) << "  PCG solve time: " << ms << " ms\n"
                  << "  iters:   " << iters
                  << (iters >= maxIter ? "  (MAX HIT)" : "") << "\n";
        std::cout << std::scientific << std::setprecision(3)
                  << "  |phi - phi_ref - mean|_L2:  " << errL2 << "\n"
                  << "  |phi_ref|_L2:               " << refL2 << "\n"
                  << "  relative L2 error:          " << rel   << "\n"
                  << "  max |phi - phi_ref - mean|: " << gMaxErr << "\n";
        bool pass = (iters > 0 && iters < maxIter && double(rel) < 5e-2);
        std::cout << "  STATUS: " << (pass ? "PASS" : "FAIL") << " (converged + rel < 5e-2)\n"
                  << std::defaultfloat;
    }
}

// =============================================================================
// CLI parsing + dispatch.
// =============================================================================
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    std::string meshFile;
    int blockSize = 256, bucketSize = 64;
    TestKind test = TestKind::All;
    StretchKind stretch = StretchKind::None;
    double aCoef = 1.0, bCoef = 2.0, cCoef = 3.0;
    unsigned seed = 42;
    int    maxIter = 1000;
    double tol     = 1e-10;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if      (arg.find("--mesh=") == 0)        meshFile = arg.substr(7);
        else if (arg.find("--block-size=") == 0)  blockSize = std::stoi(arg.substr(13));
        else if (arg.find("--bucket-size=") == 0) bucketSize = std::stoi(arg.substr(14));
        else if (arg.find("--test=") == 0)
        {
            std::string v = arg.substr(7);
            if      (v == "linear")    test = TestKind::Linear;
            else if (v == "sign")      test = TestKind::Sign;
            else if (v == "symmetry")  test = TestKind::Symmetry;
            else if (v == "nullspace") test = TestKind::NullSpace;
            else if (v == "psd")       test = TestKind::Psd;
            else if (v == "cg")        test = TestKind::Cg;
            else if (v == "cg-jacobi") test = TestKind::CgJacobi;
            else if (v == "all")       test = TestKind::All;
            else
            {
                if (rank == 0)
                    std::cerr << "Error: --test must be one of "
                                 "linear|sign|symmetry|nullspace|psd|cg|cg-jacobi|all, got '"
                              << v << "'\n";
                MPI_Finalize();
                return 1;
            }
        }
        else if (arg.find("--a=") == 0)        aCoef = std::stod(arg.substr(4));
        else if (arg.find("--b=") == 0)        bCoef = std::stod(arg.substr(4));
        else if (arg.find("--c=") == 0)        cCoef = std::stod(arg.substr(4));
        else if (arg.find("--seed=") == 0)     seed = (unsigned)std::stoul(arg.substr(7));
        else if (arg.find("--max-iter=") == 0) maxIter = std::stoi(arg.substr(11));
        else if (arg.find("--tol=") == 0)      tol = std::stod(arg.substr(6));
        else if (arg.find("--stretch=") == 0)
        {
            std::string v = arg.substr(10);
            if      (v == "none") stretch = StretchKind::None;
            else if (v == "sinx") stretch = StretchKind::SinX;
            else if (v == "exp")  stretch = StretchKind::Exp;
            else if (v == "tanh") stretch = StretchKind::Tanh;
            else
            {
                if (rank == 0)
                    std::cerr << "Error: --stretch must be one of none|sinx|exp|tanh, got '" << v << "'\n";
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
            std::cout << "Usage: " << argv[0] << " --mesh=FILE [options]\n"
                      << "  --test=linear|sign|symmetry|nullspace|psd|cg|cg-jacobi|all (default: all)\n"
                      << "  --stretch=none|sinx|exp|tanh   coordinate stretch on x (default: none)\n"
                      << "                                 sinx: +/-10% sinusoidal pinch (~3x cell ratio)\n"
                      << "                                 exp:  alpha=3 exponential (~20x cell ratio)\n"
                      << "                                 tanh: beta=2 double-side clustering\n"
                      << "  --a=VAL --b=VAL --c=VAL   linear-field coefficients (defaults 1, 2, 3)\n"
                      << "  --seed=N                  RNG seed for symmetry/psd (default 42)\n"
                      << "  --max-iter=N --tol=VAL    CG / PCG controls (defaults 1000, 1e-10)\n"
                      << "  --block-size=N --bucket-size=N\n";
        }
        MPI_Finalize();
        return 1;
    }

    using KeyType  = uint64_t;
    using RealType = double;

    if (rank == 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "MARS A = D M^{-1} D^T Validator (Phase E, Stages 1+2)\n";
        std::cout << "========================================\n";
        std::cout << "Mesh:     " << meshFile << "\n";
        std::cout << "MPI ranks: " << numRanks << "\n";
        std::cout << "Tests:    ";
        switch (test)
        {
            case TestKind::Linear:    std::cout << "linear\n"; break;
            case TestKind::Sign:      std::cout << "sign\n"; break;
            case TestKind::Symmetry:  std::cout << "symmetry\n"; break;
            case TestKind::NullSpace: std::cout << "nullspace\n"; break;
            case TestKind::Psd:       std::cout << "psd\n"; break;
            case TestKind::Cg:        std::cout << "cg\n"; break;
            case TestKind::CgJacobi:  std::cout << "cg-jacobi\n"; break;
            case TestKind::All:       std::cout << "all\n"; break;
        }
        std::cout << "Stretch:  ";
        switch (stretch)
        {
            case StretchKind::None: std::cout << "none (uniform mesh)\n"; break;
            case StretchKind::SinX: std::cout << "sinx (sinusoidal pinch in x, ~3x cell-size variation)\n"; break;
            case StretchKind::Exp:  std::cout << "exp  (alpha=3 exponential in x, ~20x cell-size variation)\n"; break;
            case StretchKind::Tanh: std::cout << "tanh (beta=2 two-sided cluster in x)\n"; break;
        }
        std::cout << "========================================\n\n";
    }

    AmrManager<HexTag, KeyType, RealType>::Config amrConfig;
    amrConfig.maxLevels  = 0;
    amrConfig.blockSize  = blockSize;
    amrConfig.bucketSize = bucketSize;
    AmrManager<HexTag, KeyType, RealType> amr(amrConfig);
    amr.initialize(meshFile, rank, numRanks);
    auto initT = amr.initTimings();

    if (rank == 0)
    {
        std::cout << "Initial mesh: " << amr.domain().getElementCount() << " elements, "
                  << amr.domain().getNodeCount() << " nodes\n";
        std::cout << "  Init (ms): domain " << std::fixed << initT.domainSyncTimeMs
                  << ", halo+topo " << initT.haloTopoTimeMs
                  << ", adj " << initT.adjacencyTimeMs
                  << ", coords " << initT.coordCacheTimeMs
                  << ", octree " << initT.octreeTimeMs
                  << ", total " << initT.totalMs << "\n";
    }

    DDTStepper<KeyType, RealType> s{amr.domain(), blockSize, rank, stretch};
    setupStepper<KeyType, RealType>(s);

    if (rank == 0)
        std::cout << "Owned DOFs (rank 0): " << s.numOwnedDofs << "\n";

    if (test == TestKind::Linear    || test == TestKind::All)
        runLinearTest<KeyType, RealType>(s, RealType(aCoef), RealType(bCoef), RealType(cCoef));
    if (test == TestKind::Sign      || test == TestKind::All)
        runSignTest<KeyType, RealType>(s);
    if (test == TestKind::Symmetry  || test == TestKind::All)
        runSymmetryTest<KeyType, RealType>(s, seed);
    if (test == TestKind::NullSpace || test == TestKind::All)
        runNullSpaceTest<KeyType, RealType>(s);
    if (test == TestKind::Psd       || test == TestKind::All)
        runPsdTest<KeyType, RealType>(s, seed);
    if (test == TestKind::Cg        || test == TestKind::All)
        runCgTest<KeyType, RealType>(s, maxIter, RealType(tol));
    if (test == TestKind::CgJacobi  || test == TestKind::All)
        runCgJacobiTest<KeyType, RealType>(s, maxIter, RealType(tol));

    if (rank == 0) std::cout << "\n========================================\n";

    MPI_Finalize();
    return 0;
}
