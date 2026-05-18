// Discrete gradient operator validator on a fixed cube mesh.
//
// Phase B.4 of the NS-projection plan. The pressure-correction step
//   u^{n+1} = u* - dt * grad p
// needs a grad(p) operator consistent with the same CVFEM CV / SCS-face layout
// as the implicit Laplacian and the discrete divergence (B.3); otherwise the
// projection can't make div(u^{n+1}) = 0 discretely.
//
// At owned node i with V_i = lumpedMass[i]:
//   (grad p)_i = (1/V_i) * sum SCS faces f at i:  p_face * A_f
// p_face = 0.5*(p[iL]+p[iR]); A_f points L -> R; +p_face*A_f to L, -p_face*A_f
// to R. Dual of divergence: div sums v.A (scalar), grad sums p*A (vector).
// Same scatter structure, three components instead of one.
//
// Tests (see --field):
//   linear   : p = a*x+b*y+c*z              -> grad p = (a,b,c) exact
//   sin      : p = sin(pi*x)*sin(pi*y)*sin(pi*z)
//   constant : p = 1                         -> grad p = 0 exact
//
// Per-node accumulator + reverseExchangeNodeHaloAdd (once per component) is
// the conservative closure: gradient at every owned node picks up contributions
// from neighbor ranks' local elements that touch the node even when those
// elements are outside this rank's cstone halo.

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

// Sync on lap so the wall-clock includes every rank's work.
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

enum class TestField { Linear, Sin, Constant };

// Per-element axis-aligned-box volume / 8 (one slice per corner).
template<typename RealType>
__device__ inline RealType elementVolumeOverEight(const RealType* x, const RealType* y, const RealType* z)
{
    RealType xmin=x[0],xmax=x[0],ymin=y[0],ymax=y[0],zmin=z[0],zmax=z[0];
    for (int i = 1; i < 8; ++i) {
        xmin = fmin(xmin, x[i]); xmax = fmax(xmax, x[i]);
        ymin = fmin(ymin, y[i]); ymax = fmax(ymax, y[i]);
        zmin = fmin(zmin, z[i]); zmax = fmax(zmax, z[i]);
    }
    return (xmax - xmin) * (ymax - ymin) * (zmax - zmin) * RealType(0.125);
}

// Legacy lumped-mass scatter (--use-legacy): ALL elements + ownership filter.
// Suffers the cstone-halo-gap under-count at corner-only-neighbor partitions.
template<typename KeyType, typename RealType>
__global__ void computeLumpedMassKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    const int* nodeToDof, const uint8_t* ownership,
    RealType* mass, size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    KeyType n[8] = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};
    RealType x[8], y[8], z[8];
    for (int i = 0; i < 8; ++i) { x[i] = nodeX[n[i]]; y[i] = nodeY[n[i]]; z[i] = nodeZ[n[i]]; }
    RealType contrib = elementVolumeOverEight<RealType>(x, y, z);
    for (int i = 0; i < 8; ++i)
    {
        if (ownership[n[i]] != 1) continue;
        int dof = nodeToDof[n[i]];
        if (dof >= 0) atomicAdd(&mass[dof], contrib);
    }
}

// Per-NODE mass: LOCAL elements only, V/8 to each corner slot. Caller follows
// with reverseExchangeNodeHaloAdd. Keeps V consistent with the gradient scatter.
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
    RealType contrib = elementVolumeOverEight<RealType>(x, y, z);
    for (int i = 0; i < 8; ++i) atomicAdd(&massNode[n[i]], contrib);
}

// Owned-node mass -> per-DOF (the shape the normalize kernel expects).
template<typename RealType>
__global__ void gatherOwnedNodeMassToDofKernel(
    const RealType* massNode, const int* nodeToDof, const uint8_t* ownership,
    RealType* massDof, size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;
    massDof[dof] = massNode[i];
}

// Pressure set once at setup on every local node slot (owned + ghost) so the
// scatter can read p[iL], p[iR] for any face after a halo exchange.
template<typename RealType>
__global__ void setTestPressureKernel(
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    RealType* p, size_t numNodes, int fieldKind,
    RealType a, RealType b, RealType c)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    RealType x = nodeX[i], y = nodeY[i], z = nodeZ[i];
    if (fieldKind == 0)       p[i] = a * x + b * y + c * z;
    else if (fieldKind == 1)
    {
        constexpr RealType PI = RealType(3.14159265358979323846);
        p[i] = sin(PI * x) * sin(PI * y) * sin(PI * z);
    }
    else                      p[i] = RealType(1);
}

// Gradient scatter, per-NODE accumulator. For each SCS face f of element e
// between corners (iL, iR) with area vector A_f pointing L -> R:
//     p_face = 0.5*(p[iL]+p[iR]);  gradAcc[L] += +p_face*A_f;  gradAcc[R] -= p_face*A_f
// Three separate buffers because CUDA has no native 3-vector atomicAdd. LOCAL
// elements only -- each global face scatters once (from the rank owning its
// element); ghost-slot contributions are reclaimed by reverseExchangeNodeHaloAdd.
template<typename KeyType, typename RealType>
__global__ void computeGradientPerNodeKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
    const RealType* p,
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

// Legacy per-rank kernel (--use-legacy): misses contributions at corner-only-
// neighbor nodes when the cstone halo doesn't include the other rank's element
// fan (see mars_amr_div.cu for the long form of the bug).
template<typename KeyType, typename RealType>
__global__ void computeGradientKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
    const RealType* p,
    const RealType* areaVecX, const RealType* areaVecY, const RealType* areaVecZ,
    const int* nodeToDof, const uint8_t* ownership,
    RealType* gxAcc, RealType* gyAcc, RealType* gzAcc, size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;
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

        if (ownership[iL] == 1)
        {
            int dofL = nodeToDof[iL];
            if (dofL >= 0)
            {
                atomicAdd(&gxAcc[dofL], +cx);
                atomicAdd(&gyAcc[dofL], +cy);
                atomicAdd(&gzAcc[dofL], +cz);
            }
        }
        if (ownership[iR] == 1)
        {
            int dofR = nodeToDof[iR];
            if (dofR >= 0)
            {
                atomicAdd(&gxAcc[dofR], -cx);
                atomicAdd(&gyAcc[dofR], -cy);
                atomicAdd(&gzAcc[dofR], -cz);
            }
        }
    }
}

// Divide accumulator by V. Legacy variant uses per-DOF accumulator.
template<typename RealType>
__global__ void applyGradientNormalizeKernel(
    const RealType* gxAcc, const RealType* gyAcc, const RealType* gzAcc,
    const RealType* lumpedMass, const int* nodeToDof, const uint8_t* ownership,
    RealType* gx, RealType* gy, RealType* gz, size_t numNodes)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    int dof = (ownership[i] == 1) ? nodeToDof[i] : -1;
    if (dof < 0) { gx[i] = gy[i] = gz[i] = RealType(0); return; }
    RealType invV = RealType(1) / lumpedMass[dof];
    gx[i] = gxAcc[dof] * invV;
    gy[i] = gyAcc[dof] * invV;
    gz[i] = gzAcc[dof] * invV;
}

template<typename RealType>
__global__ void applyGradientNormalizePerNodeKernel(
    const RealType* gxAccNode, const RealType* gyAccNode, const RealType* gzAccNode,
    const RealType* lumpedMass, const int* nodeToDof, const uint8_t* ownership,
    RealType* gx, RealType* gy, RealType* gz, size_t numNodes)
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

// Per-owned-DOF |g - g_ana|^2 * V and |g_ana|^2 * V for L2 reductions. Split
// interior vs boundary because the CV gradient at a boundary node picks up
// the boundary face flux of p, biasing the value for non-zero-boundary fields.
// Hardcoded unit cube [0,1]^3 with eps=1e-3.
template<typename RealType>
__global__ void computeGradErrorKernel(
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    const RealType* gx, const RealType* gy, const RealType* gz,
    const int* nodeToDof, const uint8_t* ownership, const RealType* mass,
    RealType* errSqInt, RealType* anaSqInt, RealType* numSqInt, RealType* numSqBnd,
    uint8_t* boundaryMask, size_t numNodes, int fieldKind,
    RealType a, RealType b, RealType c)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    if (ownership[i] != 1) return;
    int dof = nodeToDof[i];
    if (dof < 0) return;

    RealType x = nodeX[i], y = nodeY[i], z = nodeZ[i];
    RealType ax = RealType(0), ay = RealType(0), az = RealType(0);
    if (fieldKind == 0) { ax = a; ay = b; az = c; }
    else if (fieldKind == 1)
    {
        constexpr RealType PI = RealType(3.14159265358979323846);
        RealType sx = sin(PI * x), cx = cos(PI * x);
        RealType sy = sin(PI * y), cy = cos(PI * y);
        RealType sz = sin(PI * z), cz = cos(PI * z);
        ax = PI * cx * sy * sz;
        ay = PI * sx * cy * sz;
        az = PI * sx * sy * cz;
    }

    RealType nx = gx[i], ny = gy[i], nz = gz[i];
    RealType dx = nx - ax, dy = ny - ay, dz = nz - az;
    RealType errSq = dx * dx + dy * dy + dz * dz;
    RealType anaSq = ax * ax + ay * ay + az * az;
    RealType numSq = nx * nx + ny * ny + nz * nz;
    RealType m = mass[dof];

    constexpr RealType eps = RealType(1e-3);
    bool isBnd = (fabs(x) < eps) || (fabs(x - RealType(1)) < eps) ||
                 (fabs(y) < eps) || (fabs(y - RealType(1)) < eps) ||
                 (fabs(z) < eps) || (fabs(z - RealType(1)) < eps);
    boundaryMask[dof] = isBnd ? uint8_t(1) : uint8_t(0);

    if (isBnd) { numSqBnd[dof] = numSq * m; }
    else
    {
        errSqInt[dof] = errSq * m;
        anaSqInt[dof] = anaSq * m;
        numSqInt[dof] = numSq * m;
    }
}

// Owned-node MINLOC probe -- closest node to (tx,ty,tz), returns its (gx,gy,gz).
template<typename RealType>
struct ProbeGradAt
{
    RealType targetX, targetY, targetZ;
    const RealType* nodeX; const RealType* nodeY; const RealType* nodeZ;
    const RealType* gx;    const RealType* gy;    const RealType* gz;
    const int* nodeToDof;  const uint8_t* ownership;

    struct Result { RealType d2; RealType gx; RealType gy; RealType gz; };

    __device__ Result operator()(size_t i) const
    {
        if (ownership[i] != 1) return {RealType(1e30), RealType(0), RealType(0), RealType(0)};
        int dof = nodeToDof[i];
        if (dof < 0)             return {RealType(1e30), RealType(0), RealType(0), RealType(0)};
        RealType dx = nodeX[i] - targetX;
        RealType dy = nodeY[i] - targetY;
        RealType dz = nodeZ[i] - targetZ;
        return {dx * dx + dy * dy + dz * dz, gx[i], gy[i], gz[i]};
    }
};

template<typename RealType>
struct MinByDist
{
    using R = typename ProbeGradAt<RealType>::Result;
    __device__ R operator()(const R& a, const R& b) const { return (a.d2 < b.d2) ? a : b; }
};

template<typename RealType>
void probeGradAt(const cstone::DeviceVector<RealType>& d_x,
                 const cstone::DeviceVector<RealType>& d_y,
                 const cstone::DeviceVector<RealType>& d_z,
                 const cstone::DeviceVector<RealType>& d_gx,
                 const cstone::DeviceVector<RealType>& d_gy,
                 const cstone::DeviceVector<RealType>& d_gz,
                 const cstone::DeviceVector<int>& d_nodeToDof,
                 const cstone::DeviceVector<uint8_t>& d_ownership,
                 size_t numNodes, RealType tx, RealType ty, RealType tz,
                 RealType& outGx, RealType& outGy, RealType& outGz)
{
    using R = typename ProbeGradAt<RealType>::Result;
    ProbeGradAt<RealType> op{tx, ty, tz,
                             d_x.data(), d_y.data(), d_z.data(),
                             d_gx.data(), d_gy.data(), d_gz.data(),
                             d_nodeToDof.data(), d_ownership.data()};
    R init{RealType(1e30), RealType(0), RealType(0), RealType(0)};
    R best = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(numNodes),
        op, init, MinByDist<RealType>());
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    struct { double dist; int rank; } in{double(best.d2), rank}, out{};
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    double v[3] = {double(best.gx), double(best.gy), double(best.gz)};
    MPI_Bcast(v, 3, MPI_DOUBLE, out.rank, MPI_COMM_WORLD);
    outGx = RealType(v[0]); outGy = RealType(v[1]); outGz = RealType(v[2]);
}

template<typename KeyType, typename RealType>
struct GradStepper
{
    using DomainT = ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>;
    DomainT& domain;
    int blockSize;
    int rank;
    bool useLegacy = false;

    size_t nodeCount = 0;
    size_t elementCount = 0;
    int numOwnedDofs = 0;

    cstone::DeviceVector<int> d_node_to_dof;
    cstone::DeviceVector<RealType> d_mass;
    cstone::DeviceVector<RealType> d_p;
    cstone::DeviceVector<RealType> d_areaVec_x;
    cstone::DeviceVector<RealType> d_areaVec_y;
    cstone::DeviceVector<RealType> d_areaVec_z;
};

template<typename KeyType, typename RealType>
void setupGradStepper(GradStepper<KeyType, RealType>& s, TestField field,
                      RealType a, RealType b, RealType c)
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
    s.numOwnedDofs = buildDofMappingGpu<KeyType>(d_nodeOwnership.data(),
                                                 s.d_node_to_dof.data(), s.nodeCount);
    pt.lap("DOF mapping");

    // Connectivity arrays used everywhere downstream.
    auto c0 = std::get<0>(d_conn).data(); auto c1 = std::get<1>(d_conn).data();
    auto c2 = std::get<2>(d_conn).data(); auto c3 = std::get<3>(d_conn).data();
    auto c4 = std::get<4>(d_conn).data(); auto c5 = std::get<5>(d_conn).data();
    auto c6 = std::get<6>(d_conn).data(); auto c7 = std::get<7>(d_conn).data();

    // Same area-vector source as the implicit assembler, advection scatter, divergence.
    s.d_areaVec_x.resize(s.elementCount * 12);
    s.d_areaVec_y.resize(s.elementCount * 12);
    s.d_areaVec_z.resize(s.elementCount * 12);
    precomputeAreaVectorsGpu<KeyType, RealType>(
        c0, c1, c2, c3, c4, c5, c6, c7, s.elementCount,
        d_x.data(), d_y.data(), d_z.data(),
        s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data());
    pt.lap("area vectors");

    s.d_p.resize(s.nodeCount);
    {
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        int fk = (field == TestField::Linear) ? 0 : (field == TestField::Sin) ? 1 : 2;
        setTestPressureKernel<RealType><<<nBlocks, s.blockSize>>>(
            d_x.data(), d_y.data(), d_z.data(),
            s.d_p.data(), s.nodeCount, fk, a, b, c);
        cudaDeviceSynchronize();
    }
    pt.lap("test pressure");

    s.d_mass.resize(s.numOwnedDofs);
    thrust::fill(thrust::device_pointer_cast(s.d_mass.data()),
                 thrust::device_pointer_cast(s.d_mass.data() + s.numOwnedDofs),
                 RealType(0));
    if (s.useLegacy)
    {
        int eBlocks = (s.elementCount + s.blockSize - 1) / s.blockSize;
        computeLumpedMassKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            c0, c1, c2, c3, c4, c5, c6, c7,
            d_x.data(), d_y.data(), d_z.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            s.d_mass.data(), s.elementCount);
        cudaDeviceSynchronize();
    }
    else
    {
        // Same per-node + reverse-halo closure as the gradient scatter --
        // matches V to the operator so single-rank == multi-rank per owned node.
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

// Default path: per-NODE accumulator + reverseExchangeNodeHaloAdd, called
// once per component (the method is templated on a single per-node array).
// --use-legacy: per-DOF + ownership-filter for diagnostic comparison.
template<typename KeyType, typename RealType>
void computeGradient(GradStepper<KeyType, RealType>& s,
                     cstone::DeviceVector<RealType>& d_gx,
                     cstone::DeviceVector<RealType>& d_gy,
                     cstone::DeviceVector<RealType>& d_gz)
{
    const auto& d_nodeOwnership = s.domain.getNodeOwnershipMap();
    const auto& d_conn          = s.domain.getElementToNodeConnectivity();

    // Pressure ghosts must be in sync -- face donors touching ghost corners read p there.
    s.domain.exchangeNodeHalo(s.d_p);

    d_gx.resize(s.nodeCount);
    d_gy.resize(s.nodeCount);
    d_gz.resize(s.nodeCount);

    auto c0 = std::get<0>(d_conn).data(); auto c1 = std::get<1>(d_conn).data();
    auto c2 = std::get<2>(d_conn).data(); auto c3 = std::get<3>(d_conn).data();
    auto c4 = std::get<4>(d_conn).data(); auto c5 = std::get<5>(d_conn).data();
    auto c6 = std::get<6>(d_conn).data(); auto c7 = std::get<7>(d_conn).data();

    if (s.useLegacy)
    {
        cstone::DeviceVector<RealType> d_gxAcc(s.numOwnedDofs, RealType(0));
        cstone::DeviceVector<RealType> d_gyAcc(s.numOwnedDofs, RealType(0));
        cstone::DeviceVector<RealType> d_gzAcc(s.numOwnedDofs, RealType(0));

        int eBlocks = (s.elementCount + s.blockSize - 1) / s.blockSize;
        computeGradientKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
            c0, c1, c2, c3, c4, c5, c6, c7,
            s.d_p.data(),
            s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(), s.elementCount);
        cudaDeviceSynchronize();

        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        applyGradientNormalizeKernel<RealType><<<nBlocks, s.blockSize>>>(
            d_gxAcc.data(), d_gyAcc.data(), d_gzAcc.data(),
            s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            d_gx.data(), d_gy.data(), d_gz.data(), s.nodeCount);
        cudaDeviceSynchronize();
    }
    else
    {
        cstone::DeviceVector<RealType> d_gxAccNode(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> d_gyAccNode(s.nodeCount, RealType(0));
        cstone::DeviceVector<RealType> d_gzAccNode(s.nodeCount, RealType(0));
        size_t startElem = s.domain.startIndex();
        size_t numLocal  = s.domain.localElementCount();
        if (numLocal > 0)
        {
            int eBlocks = int((numLocal + s.blockSize - 1) / s.blockSize);
            computeGradientPerNodeKernel<KeyType, RealType><<<eBlocks, s.blockSize>>>(
                c0, c1, c2, c3, c4, c5, c6, c7,
                s.d_p.data(),
                s.d_areaVec_x.data(), s.d_areaVec_y.data(), s.d_areaVec_z.data(),
                d_gxAccNode.data(), d_gyAccNode.data(), d_gzAccNode.data(),
                startElem, numLocal);
            cudaDeviceSynchronize();
        }
        // One reverse-halo per component -- the method templates on a single
        // per-node array. Each call sums ghost-slot contributions back to owner.
        s.domain.reverseExchangeNodeHaloAdd(d_gxAccNode);
        s.domain.reverseExchangeNodeHaloAdd(d_gyAccNode);
        s.domain.reverseExchangeNodeHaloAdd(d_gzAccNode);
        int nBlocks = (s.nodeCount + s.blockSize - 1) / s.blockSize;
        applyGradientNormalizePerNodeKernel<RealType><<<nBlocks, s.blockSize>>>(
            d_gxAccNode.data(), d_gyAccNode.data(), d_gzAccNode.data(),
            s.d_mass.data(),
            s.d_node_to_dof.data(), d_nodeOwnership.data(),
            d_gx.data(), d_gy.data(), d_gz.data(), s.nodeCount);
        cudaDeviceSynchronize();
    }

    // Fill ghosts with owner values so VTU output is continuous across ranks.
    s.domain.exchangeNodeHalo(d_gx);
    s.domain.exchangeNodeHalo(d_gy);
    s.domain.exchangeNodeHalo(d_gz);
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
    int blockSize = 256, bucketSize = 64;
    TestField field = TestField::Linear;
    double aCoef = 1.0, bCoef = 2.0, cCoef = 3.0;
    std::string vtuPrefix;
    bool useLegacy = false;
    double probeX = 0.5, probeY = 0.5, probeZ = 0.5;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg.find("--mesh=") == 0)              meshFile = arg.substr(7);
        else if (arg.find("--block-size=") == 0)   blockSize = std::stoi(arg.substr(13));
        else if (arg.find("--bucket-size=") == 0)  bucketSize = std::stoi(arg.substr(14));
        else if (arg.find("--field=") == 0)
        {
            std::string v = arg.substr(8);
            if      (v == "linear")   field = TestField::Linear;
            else if (v == "sin")      field = TestField::Sin;
            else if (v == "constant") field = TestField::Constant;
            else
            {
                if (rank == 0)
                    std::cerr << "Error: --field must be linear|sin|constant, got '" << v << "'\n";
                MPI_Finalize();
                return 1;
            }
        }
        else if (arg.find("--a=") == 0)            aCoef = std::stod(arg.substr(4));
        else if (arg.find("--b=") == 0)            bCoef = std::stod(arg.substr(4));
        else if (arg.find("--c=") == 0)            cCoef = std::stod(arg.substr(4));
        else if (arg.find("--vtu-output=") == 0)   vtuPrefix = arg.substr(13);
        else if (arg == "--use-legacy")            useLegacy = true;
        else if (arg.find("--probe-x=") == 0)      probeX = std::stod(arg.substr(10));
        else if (arg.find("--probe-y=") == 0)      probeY = std::stod(arg.substr(10));
        else if (arg.find("--probe-z=") == 0)      probeZ = std::stod(arg.substr(10));
        else if (arg[0] != '-' && meshFile.empty()) meshFile = arg;
    }

    if (meshFile.empty())
    {
        if (rank == 0)
        {
            std::cout << "Usage: " << argv[0] << " --mesh=FILE [options]\n"
                      << "  --field=linear|sin|constant  --a/b/c=VALUE (linear coefs)\n"
                      << "  --probe-x/y/z=VAL            --vtu-output=PREFIX\n"
                      << "  --bucket-size=N --block-size=N --use-legacy\n";
        }
        MPI_Finalize();
        return 1;
    }

    using KeyType  = uint64_t;
    using RealType = double;

    if (rank == 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "MARS Gradient Operator Validator (CVFEM SCS, Phase B.4)\n";
        std::cout << "========================================\n";
        if (field == TestField::Linear)
            std::cout << "Field: linear, p = " << aCoef << "*x + " << bCoef << "*y + " << cCoef
                      << "*z  =>  grad p = (" << aCoef << ", " << bCoef << ", " << cCoef << ")\n";
        else if (field == TestField::Sin)
            std::cout << "Field: sin, p = sin(pi*x)*sin(pi*y)*sin(pi*z)\n";
        else
            std::cout << "Field: constant, p = 1  =>  grad p = (0, 0, 0)\n";
        std::cout << "Mesh: " << meshFile << ", " << numRanks << " MPI ranks\n";
        std::cout << "========================================\n\n";
    }

    // maxLevels=0 -- AmrManager only used for mesh load + lazy halo + coords.
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
                  << ", total " << initT.totalMs << "\n\n";
    }

    GradStepper<KeyType, RealType> s{amr.domain(), blockSize, rank, useLegacy};
    setupGradStepper<KeyType, RealType>(s, field, RealType(aCoef), RealType(bCoef), RealType(cCoef));

    if (rank == 0)
    {
        std::cout << "Scatter mode: " << (useLegacy ? "LEGACY (per-DOF + ownership filter)"
                                                    : "PER-NODE + reverse halo (conservative)") << "\n\n";
    }

    cstone::DeviceVector<RealType> d_gx, d_gy, d_gz;
    auto t0 = std::chrono::high_resolution_clock::now();
    computeGradient<KeyType, RealType>(s, d_gx, d_gy, d_gz);
    cudaDeviceSynchronize();
    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = std::chrono::high_resolution_clock::now();
    float gradMs = std::chrono::duration<float, std::milli>(t1 - t0).count();

    int fieldKind = (field == TestField::Linear) ? 0 : (field == TestField::Sin) ? 1 : 2;

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
        computeGradErrorKernel<RealType><<<nBlocks, blockSize>>>(
            d_x.data(), d_y.data(), d_z.data(),
            d_gx.data(), d_gy.data(), d_gz.data(),
            s.d_node_to_dof.data(), d_nodeOwn.data(), s.d_mass.data(),
            d_errSqInt.data(), d_anaSqInt.data(), d_numSqInt.data(),
            d_numSqBnd.data(), d_bndMask.data(),
            s.nodeCount, fieldKind,
            RealType(aCoef), RealType(bCoef), RealType(cCoef));
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

    // Max |grad| (interior + boundary) and max |grad - ana| (interior).
    // Three single-value max-reductions; analytical recomputed inline.
    RealType locMaxInt = 0, locMaxBnd = 0, locMaxErrInt = 0;
    {
        const auto& d_nodeOwn = amr.domain().getNodeOwnershipMap();
        const auto& d_x = amr.domain().getNodeX();
        const auto& d_y = amr.domain().getNodeY();
        const auto& d_z = amr.domain().getNodeZ();
        const uint8_t* ownPtr = d_nodeOwn.data();
        const int* dofPtr     = s.d_node_to_dof.data();
        const RealType* gxPtr = d_gx.data(), *gyPtr = d_gy.data(), *gzPtr = d_gz.data();
        const RealType* xPtr  = d_x.data(),  *yPtr  = d_y.data(),  *zPtr  = d_z.data();
        const uint8_t* bndPtr = d_bndMask.data();
        const int fk = fieldKind;
        const RealType ka = RealType(aCoef), kb = RealType(bCoef), kc = RealType(cCoef);
        auto countBegin = thrust::counting_iterator<size_t>(0);
        auto countEnd   = thrust::counting_iterator<size_t>(s.nodeCount);

        locMaxInt = thrust::transform_reduce(thrust::device, countBegin, countEnd,
            [ownPtr, dofPtr, gxPtr, gyPtr, gzPtr, bndPtr] __device__ (size_t i) -> RealType {
                if (ownPtr[i] != 1) return RealType(0);
                int dof = dofPtr[i];
                if (dof < 0 || bndPtr[dof] != 0) return RealType(0);
                RealType gx = gxPtr[i], gy = gyPtr[i], gz = gzPtr[i];
                return sqrt(gx * gx + gy * gy + gz * gz);
            }, RealType(0), thrust::maximum<RealType>());
        locMaxBnd = thrust::transform_reduce(thrust::device, countBegin, countEnd,
            [ownPtr, dofPtr, gxPtr, gyPtr, gzPtr, bndPtr] __device__ (size_t i) -> RealType {
                if (ownPtr[i] != 1) return RealType(0);
                int dof = dofPtr[i];
                if (dof < 0 || bndPtr[dof] == 0) return RealType(0);
                RealType gx = gxPtr[i], gy = gyPtr[i], gz = gzPtr[i];
                return sqrt(gx * gx + gy * gy + gz * gz);
            }, RealType(0), thrust::maximum<RealType>());
        locMaxErrInt = thrust::transform_reduce(thrust::device, countBegin, countEnd,
            [ownPtr, dofPtr, gxPtr, gyPtr, gzPtr, xPtr, yPtr, zPtr, bndPtr, fk, ka, kb, kc] __device__ (size_t i) -> RealType {
                if (ownPtr[i] != 1) return RealType(0);
                int dof = dofPtr[i];
                if (dof < 0 || bndPtr[dof] != 0) return RealType(0);
                RealType ax = RealType(0), ay = RealType(0), az = RealType(0);
                if (fk == 0) { ax = ka; ay = kb; az = kc; }
                else if (fk == 1)
                {
                    constexpr RealType PI = RealType(3.14159265358979323846);
                    RealType sx = sin(PI * xPtr[i]), cx = cos(PI * xPtr[i]);
                    RealType sy = sin(PI * yPtr[i]), cy = cos(PI * yPtr[i]);
                    RealType sz = sin(PI * zPtr[i]), cz = cos(PI * zPtr[i]);
                    ax = PI * cx * sy * sz;
                    ay = PI * sx * cy * sz;
                    az = PI * sx * sy * cz;
                }
                RealType dx = gxPtr[i] - ax, dy = gyPtr[i] - ay, dz = gzPtr[i] - az;
                return sqrt(dx * dx + dy * dy + dz * dz);
            }, RealType(0), thrust::maximum<RealType>());
    }
    RealType gMaxInt = 0, gMaxBnd = 0, gMaxErrInt = 0;
    MPI_Allreduce(&locMaxInt,    &gMaxInt,    1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locMaxBnd,    &gMaxBnd,    1, mpiType, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&locMaxErrInt, &gMaxErrInt, 1, mpiType, MPI_MAX, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "Gradient computed in " << std::fixed << std::setprecision(2)
                  << gradMs << " ms\n\n";
        std::cout << "========================================\n";
        if (field == TestField::Linear)
            std::cout << "Test: linear, a=" << aCoef << ", b=" << bCoef << ", c=" << cCoef << "\n";
        else
            std::cout << "Test: " << (field == TestField::Sin ? "sin" : "constant") << "\n";
        std::cout << std::scientific << std::setprecision(3)
                  << "  Numerical  |grad p|_L2 (all owned):     " << std::sqrt(gNumAll)
                  << "   (includes boundary CV through-flux)\n"
                  << "  Numerical  |grad p|_L2 (interior only): " << std::sqrt(gNumInt) << "\n"
                  << "  Analytical |grad p|_L2 (interior only): " << std::sqrt(gAnaInt) << "\n"
                  << "  L2 error (interior only):               " << std::sqrt(gErrInt) << "\n"
                  << "  Max |grad p - analytical| interior:     " << gMaxErrInt << "\n"
                  << "  Max |grad p| interior:                  " << gMaxInt << "\n"
                  << "  Max |grad p| boundary:                  " << gMaxBnd << "\n"
                  << std::defaultfloat;
    }

    {
        // Probe: linear -> (a,b,c) to roundoff; sin -> spot accuracy; constant -> ~0.
        const auto& d_x = amr.domain().getNodeX();
        const auto& d_y = amr.domain().getNodeY();
        const auto& d_z = amr.domain().getNodeZ();
        const auto& d_nodeOwn = amr.domain().getNodeOwnershipMap();
        RealType pgx = 0, pgy = 0, pgz = 0;
        probeGradAt<RealType>(d_x, d_y, d_z, d_gx, d_gy, d_gz,
                              s.d_node_to_dof, d_nodeOwn, s.nodeCount,
                              RealType(probeX), RealType(probeY), RealType(probeZ),
                              pgx, pgy, pgz);
        if (rank == 0)
            std::cout << std::fixed << std::setprecision(4)
                      << "  Probe at (" << probeX << ", " << probeY << ", " << probeZ
                      << "): grad = (" << pgx << ", " << pgy << ", " << pgz << ")\n"
                      << std::defaultfloat;
    }
    if (rank == 0) std::cout << "========================================\n";

    if (!vtuPrefix.empty())
    {
        // VTU writer takes one scalar -- emit |grad p| (magnitude) for visualization.
        cstone::DeviceVector<RealType> d_gmag(s.nodeCount, RealType(0));
        const RealType* gxPtr = d_gx.data();
        const RealType* gyPtr = d_gy.data();
        const RealType* gzPtr = d_gz.data();
        RealType* mPtr        = d_gmag.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [gxPtr, gyPtr, gzPtr, mPtr] __device__(size_t i) {
                RealType gx = gxPtr[i], gy = gyPtr[i], gz = gzPtr[i];
                mPtr[i] = sqrt(gx * gx + gy * gy + gz * gz);
            });
        cudaDeviceSynchronize();
        fem::VTUParallelWriter<KeyType, RealType> vtuWriter(vtuPrefix);
        vtuWriter.writeFrame(0, 0.0, amr.domain(), d_gmag, "gradMag");
        if (rank == 0)
            std::cout << "\nVTU written: " << vtuPrefix << "_step0000.pvtu, "
                      << vtuPrefix << ".pvd\n";
    }

    MPI_Finalize();
    return 0;
}
