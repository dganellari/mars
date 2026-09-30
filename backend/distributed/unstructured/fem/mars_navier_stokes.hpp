#pragma once

// Incompressible Navier-Stokes on hex meshes: CVFEM, one GPU per MPI rank. The
// solver behind the Poiseuille, Taylor-Green and lid-driven cavity examples.
//
// Velocity and pressure live at the nodes (equal order, collocated). Each hex
// contributes 12 sub-control faces f (area vector A_f from node L to node R)
// and 8 sub-control volumes (lumped mass M). The flow through a face is a
// stored flux F_f, and continuity is the balance of these fluxes:
//   D_F   face flux divergence   (D_F F)_L += F_f, (D_F F)_R -= F_f, plus inlets and outlets
//   G     nodal gradient         G p = -M^-1 D^T p, (D^T p)_L = (D^T p)_R += (p_L - p_R) / 2 A_f
//   K     CVFEM Laplacian        -(K p)_L = sum_f (grad p . A)_f, compact face gradients
//   N     skew-symmetric advection by the face fluxes F
//   Q     zeroes the velocity at nodes where it is prescribed
//
// On equal-order nodes the plain average A . (u_L + u_R) / 2 leaves the
// pressure checkerboard invisible. The stabilized (Rhie-Chow) face flux
//   F = A . (u_L + u_R) / 2 - h [ (grad p . A)_f - A . (G p_L + G p_R) / 2 ],  h = dtEff / rho,
// adds the difference between the compact and the averaged nodal pressure
// gradient: O(h^2) for smooth p, large for a checkerboard.
//
// One time step, BDF2 with extrapolated advection (BDF1 on the first step),
// dtEff = 2 dt / 3 (dt with BDF1):
//   predictor   u*  = BDF history + dtEff M^-1 (N_ext(F^n) - rho^-1 M Q G p^n)
//   viscous     (M / dtEff + nu K) u** = M u* / dtEff,    u** = prescribed where fixed
//   projection  K phi = -(rho / dtEff) (D_F F** + openings),   F** the stabilized flux of u**, p^n
//   corrector   F = F** - h (grad phi . A),  u = u** - h Q G phi,  p += phi
// K is the operator of the flux correction, so D_F F = 0 to the solver
// tolerance on any number of ranks. The nodal velocity follows the flux up to
// the stabilization.
//
// Parallel layout. Each rank assembles its own elements into its node slots:
// owned nodes, ghost copies of other ranks' nodes and, on a periodic mesh,
// several slots for one periodic point. The DofSpace maps slots to the unknowns
// (P) and back (P^T). Time step kernels scatter over the local elements and
// restrict; the constant matrices are assembled per rank and reduced by Hypre
// to P^T A_local P, the way MFEM assembles. Both systems are solved with PCG
// and BoomerAMG, whose iteration count stays flat as the mesh and the number
// of GPUs grow.

#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_hex_kernel.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_dof_space.hpp"
#include "backend/distributed/unstructured/solvers/mars_hypre_amg_pcg_solver.hpp"

#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/tuple.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mpi.h>
#include <type_traits>
#include <vector>

namespace mars
{
namespace fem
{

// ---------------------------------------------------------------------------
// Problem description
// ---------------------------------------------------------------------------

// Conditions at one node, from its position. The solver copies the value at
// the DOF slot to every copy of the node, so a periodic image takes the
// conditions of its master.
template<typename RealType>
struct NodeCondition
{
    bool velocityFixed   = false; // velocity prescribed
    RealType velocity[3] = {0, 0, 0};
    bool pressureFixed   = false; // p = 0
};

// No boundary conditions (a fully periodic box).
template<typename RealType>
struct FreeNodes
{
    __device__ NodeCondition<RealType> operator()(RealType, RealType, RealType) const { return {}; }
};

// An inlet or an outlet: the element faces on the plane x[axis] = position.
// The elements lie on the side +1 (x[axis] > position) or -1; the outward
// normal points the other way. Fluid crosses these faces; everywhere else the
// boundary is closed (walls) or periodic.
template<typename RealType>
struct Opening
{
    int axis;
    RealType position;
    int side;
};

template<typename RealType>
struct Box
{
    RealType lo[3];
    RealType hi[3];
    RealType extent() const { return std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]}); }
};

// Global bounding box of the node coordinates.
template<typename Domain>
auto boundingBox(const Domain& domain)
{
    using RealType           = std::decay_t<decltype(*domain.getNodeX().data())>;
    const RealType* coord[3] = {domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data()};
    const size_t n           = domain.getNodeCount();
    const RealType inf       = std::numeric_limits<RealType>::infinity();
    Box<RealType> box;
    for (int d = 0; d < 3; ++d)
    {
        RealType lo = thrust::reduce(thrust::device, coord[d], coord[d] + n, inf, thrust::minimum<RealType>());
        RealType hi = thrust::reduce(thrust::device, coord[d], coord[d] + n, -inf, thrust::maximum<RealType>());
        MPI_Allreduce(&lo, &box.lo[d], 1, mpiDatatype<RealType>(), MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&hi, &box.hi[d], 1, mpiDatatype<RealType>(), MPI_MAX, MPI_COMM_WORLD);
    }
    return box;
}

// Elements [first, first + count); corner c of element e is node[c][e].
template<typename KeyType>
struct HexElements
{
    const KeyType* node[8];
    size_t first;
    size_t count;
};

// Up to three pointers, one per velocity component.
template<typename T>
struct Components
{
    T c[3];
};


// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

// det J of the trilinear map of a hex at reference point xi in [-1, 1]^3.
// Corner order: (-,-,-) (+,-,-) (+,+,-) (-,+,-) (-,-,+) (+,-,+) (+,+,+) (-,+,+).
template<typename RealType>
__host__ __device__ inline RealType hexJacobianDeterminant(const RealType x[8][3], const RealType xi[3])
{
    const int sign[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                            {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
    RealType J[3][3] = {};
    for (int n = 0; n < 8; ++n)
    {
        RealType f[3] = {1 + sign[n][0] * xi[0], 1 + sign[n][1] * xi[1], 1 + sign[n][2] * xi[2]};
        RealType dN[3] = {sign[n][0] * f[1] * f[2] / 8, sign[n][1] * f[0] * f[2] / 8, sign[n][2] * f[0] * f[1] / 8};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
                J[a][b] += x[n][a] * dN[b];
    }
    return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
           J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
}

// Volume of the sub-control volume of corner c: the image of the reference
// octant between the corner and the element centre. det J is quadratic in each
// reference coordinate, so the 2x2x2 Gauss rule on the octant is exact.
template<typename RealType>
__host__ __device__ inline RealType hexSubVolume(const RealType x[8][3], int c)
{
    const int sign[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                            {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
    const RealType g[2]  = {RealType(0.5) - RealType(0.5) / sqrt(RealType(3)),
                            RealType(0.5) + RealType(0.5) / sqrt(RealType(3))};
    RealType volume      = 0;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
            {
                RealType xi[3] = {sign[c][0] * g[i], sign[c][1] * g[j], sign[c][2] * g[k]};
                volume += hexJacobianDeterminant(x, xi);
            }
    return volume / 8; // octant volume 1, eight points
}

template<typename KeyType, typename RealType>
__global__ void nsMassKernel(HexElements<KeyType> hex, const RealType* x, const RealType* y, const RealType* z,
                             RealType* mass)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    RealType corner[8][3];
    for (int c = 0; c < 8; ++c)
    {
        KeyType n    = hex.node[c][e];
        corner[c][0] = x[n];
        corner[c][1] = y[n];
        corner[c][2] = z[n];
    }
    for (int c = 0; c < 8; ++c)
        atomicAdd(&mass[hex.node[c][e]], hexSubVolume(corner, c));
}

// Outward area of the faces on an opening, per node: each face node gets its
// quarter of the face (the sub-quad of node, edge midpoints and face centre).
// A face counts once, from the element on the opening's side.
template<typename KeyType, typename RealType>
__global__ void nsOpeningAreaKernel(HexElements<KeyType> hex, const RealType* x, const RealType* y, const RealType* z,
                                    Opening<RealType> opening, RealType tolerance, RealType* area)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e                 = hex.first + k;
    const RealType* coord[3] = {x, y, z};
    const int a = opening.axis, b = (a + 1) % 3, c = (a + 2) % 3;
    KeyType face[4];
    int onPlane = 0;
    for (int j = 0; j < 8; ++j)
    {
        KeyType n  = hex.node[j][e];
        RealType d = coord[a][n] - opening.position;
        if (fabs(d) <= tolerance)
        {
            if (onPlane < 4) face[onPlane] = n;
            ++onPlane;
        }
        else if (d * opening.side <= 0) return;
    }
    if (onPlane != 4) return;
    RealType cb = 0, cc = 0;
    for (int j = 0; j < 4; ++j)
    {
        cb += RealType(0.25) * coord[b][face[j]];
        cc += RealType(0.25) * coord[c][face[j]];
    }
    // The element's corner order does not give the face cycle; sort by angle.
    RealType angle[4];
    for (int j = 0; j < 4; ++j)
        angle[j] = atan2(coord[c][face[j]] - cc, coord[b][face[j]] - cb);
    for (int j = 1; j < 4; ++j)
        for (int i = j; i > 0 && angle[i] < angle[i - 1]; --i)
        {
            RealType t   = angle[i];
            angle[i]     = angle[i - 1];
            angle[i - 1] = t;
            KeyType f    = face[i];
            face[i]      = face[i - 1];
            face[i - 1]  = f;
        }
    for (int j = 0; j < 4; ++j)
    {
        KeyType n = face[j], next = face[(j + 1) % 4], prev = face[(j + 3) % 4];
        RealType nb    = RealType(0.5) * (coord[b][prev] - coord[b][next]);
        RealType nc    = RealType(0.5) * (coord[c][prev] - coord[c][next]);
        RealType share = RealType(0.5) * fabs((cb - coord[b][n]) * nc - (cc - coord[c][n]) * nb);
        atomicAdd(&area[n], -opening.side * share);
    }
}

// Outward face areas of an opening per node (component opening.axis), summed
// over the ranks. For measurements; the solver also folds periodic copies.
template<typename KeyType, typename RealType, typename Domain>
void openingFaceAreas(const Domain& domain, const Opening<RealType>& opening, RealType tolerance,
                      cstone::DeviceVector<RealType>& area)
{
    const size_t n = domain.getNodeCount();
    area.resize(n);
    cudaMemset(area.data(), 0, n * sizeof(RealType));
    const auto& conn = domain.getElementToNodeConnectivity();
    HexElements<KeyType> hex{{std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                              std::get<3>(conn).data(), std::get<4>(conn).data(), std::get<5>(conn).data(),
                              std::get<6>(conn).data(), std::get<7>(conn).data()},
                             domain.startIndex(),
                             domain.localElementCount()};
    if (hex.count > 0)
        nsOpeningAreaKernel<KeyType, RealType><<<int((hex.count + 255) / 256), 256>>>(
            hex, domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(), opening, tolerance,
            area.data());
    cudaCheckError();
    domain.reverseExchangeNodeHaloAdd(area);
    domain.exchangeNodeHalo(area);
}

// ---------------------------------------------------------------------------
// Setup kernels
// ---------------------------------------------------------------------------

template<typename RealType, class Rule>
__global__ void nsBoundaryKernel(size_t n, const RealType* x, const RealType* y, const RealType* z, Rule rule,
                                 RealType* flags, Components<RealType*> target)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    NodeCondition<RealType> node = rule(x[i], y[i], z[i]);
    // Both flags packed in one number, so the DofSpace can copy them to every copy of the node.
    flags[i] = RealType((node.velocityFixed ? 1 : 0) + (node.pressureFixed ? 2 : 0));
    for (int d = 0; d < 3; ++d)
        target.c[d][i] = node.velocityFixed ? node.velocity[d] : RealType(0);
}

template<typename RealType>
__global__ void nsUnpackFlagsKernel(size_t n, const RealType* flags, uint8_t* velocityFixed, uint8_t* pressureFixed)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int f            = int(flags[i] + RealType(0.5));
    velocityFixed[i] = f & 1;
    pressureFixed[i] = (f >> 1) & 1;
}

// Nodes on neither of the two planes z = lo, hi.
template<typename RealType>
struct NsOffPlanes
{
    RealType lo, hi, eps;
    __device__ bool operator()(RealType z) const { return fabs(z - lo) > eps && fabs(z - hi) > eps; }
};

template<typename RealType>
__global__ void nsGlobalIdKernel(size_t n, const uint8_t* isDof, const int* dofIndex, long long dofStart, RealType* id)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) id[i] = isDof[i] ? RealType(dofStart + dofIndex[i]) : RealType(-1);
}

// Rows of the prolongation P: slot localStart + i takes the value of DOF gid[i], per component.
template<typename RealType>
__global__ void nsProlongationTripletsKernel(size_t n, int comps, long long localStart, const RealType* gid,
                                             HYPRE_BigInt* rows, HYPRE_BigInt* cols, RealType* values)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    for (int d = 0; d < comps; ++d)
    {
        size_t k  = i * comps + d;
        rows[k]   = comps * (localStart + HYPRE_BigInt(i)) + d;
        cols[k]   = comps * HYPRE_BigInt(gid[i]) + d;
        values[k] = RealType(1);
    }
}

// With no p = 0 anywhere the pressure is only defined up to a constant: fix it
// at the DOF with global id 0.
template<typename RealType>
__global__ void nsPinFirstDofKernel(size_t n, const uint8_t* isDof, const RealType* gid, const uint8_t* pressureFixed,
                                    RealType* flag)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) flag[i] = pressureFixed[i] || (isDof[i] && gid[i] == RealType(0)) ? RealType(1) : RealType(0);
}

template<typename RealType>
__global__ void nsFlagFromRealKernel(size_t n, const RealType* flag, uint8_t* out)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) out[i] = flag[i] > RealType(0.5) ? 1 : 0;
}

template<typename KeyType, typename RealType>
__device__ inline void nsCorners(const HexElements<KeyType>& hex, size_t e, const RealType* x, const RealType* y,
                                 const RealType* z, KeyType n[8], double corner[8][3])
{
    for (int c = 0; c < 8; ++c)
    {
        n[c]         = hex.node[c][e];
        corner[c][0] = x[n[c]];
        corner[c][1] = y[n[c]];
        corner[c][2] = z[n[c]];
    }
}

// (grad phi . A) at the integration point of face ip is sum_n w[n] phi_n, with
// the shape derivatives of the library's CVFEM kernels.
template<typename RealType>
__device__ inline void nsFaceWeights(int ip, const double corner[8][3], const Components<const RealType*>& area,
                                     size_t f, double w[8])
{
    double dndx[8][3];
    computeShapeDerivatives(ip, corner, dndx);
    for (int m = 0; m < 8; ++m)
        w[m] = dndx[m][0] * area.c[0][f] + dndx[m][1] * area.c[1][f] + dndx[m][2] * area.c[2][f];
}

// The element matrix of scale * K as 64 triplets over local slots. Entries in
// a removed row or column are 0: the Galerkin product keeps them empty and the
// removed DOFs get identity rows. liftComps > 0 moves the removed columns of a
// free row to the right-hand side, -scale K_rc target_c, which keeps K symmetric.
template<typename KeyType, typename RealType>
__global__ void nsLaplacianTripletsKernel(HexElements<KeyType> hex, const RealType* x, const RealType* y,
                                          const RealType* z, Components<const RealType*> area, RealType scale,
                                          long long localStart, const uint8_t* removed,
                                          Components<const RealType*> target, Components<RealType*> lift,
                                          int liftComps, HYPRE_BigInt* rows, HYPRE_BigInt* cols, RealType* values)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    KeyType n[8];
    double corner[8][3];
    nsCorners(hex, hex.first + k, x, y, z, n, corner);
    double Ke[8][8] = {};
    for (int ip = 0; ip < 12; ++ip)
    {
        double w[8];
        nsFaceWeights(ip, corner, area, k * 12 + ip, w);
        int L = hexLRSCV[2 * ip], R = hexLRSCV[2 * ip + 1];
        for (int m = 0; m < 8; ++m)
        {
            Ke[L][m] -= w[m];
            Ke[R][m] += w[m];
        }
    }
    size_t out = k * 64;
    for (int a = 0; a < 8; ++a)
        for (int b = 0; b < 8; ++b, ++out)
        {
            KeyType r = n[a], c = n[b];
            bool drop = removed && (removed[r] || removed[c]);
            rows[out]   = localStart + HYPRE_BigInt(r);
            cols[out]   = localStart + HYPRE_BigInt(c);
            values[out] = drop ? RealType(0) : RealType(scale * Ke[a][b]);
            if (liftComps > 0 && !removed[r] && removed[c])
                for (int d = 0; d < liftComps; ++d)
                    atomicAdd(&lift.c[d][r], RealType(-scale * Ke[a][b]) * target.c[d][c]);
        }
}

// The diagonal c M over the local slots; zero in removed rows.
template<typename RealType>
__global__ void nsMassTripletsKernel(size_t n, long long localStart, const RealType* massLocal, const uint8_t* removed,
                                     RealType c, HYPRE_BigInt* rows, HYPRE_BigInt* cols, RealType* values)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    rows[i]   = localStart + HYPRE_BigInt(i);
    cols[i]   = rows[i];
    values[i] = removed[i] ? RealType(0) : c * massLocal[i];
}

template<typename RealType>
__global__ void nsTestVectorKernel(size_t n, const uint8_t* isDof, const uint8_t* pressureFixed, const RealType* gid,
                                   RealType* x)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    long long g = (long long)gid[i];
    x[i] = isDof[i] && !pressureFixed[i] ? RealType(1) + RealType((g * 7919) % 97) / RealType(97) : RealType(0);
}

// Slot values into DOF order.
template<typename T>
__global__ void nsToDofKernel(size_t n, const uint8_t* isDof, const int* dofIndex, const T* slot, T* dof)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n && isDof[i]) dof[dofIndex[i]] = slot[i];
}

// ---------------------------------------------------------------------------
// Time step kernels. Element loops run over the local elements, with the face
// data of element k at k * 12 + face; node loops over the DOF slots (isDof).
// ---------------------------------------------------------------------------

// The stabilized face flux F = A . (u_L + u_R) / 2 - h [ (grad p . A) + A . (g_L + g_R) / 2 ]
// with g = M^-1 D^T p = -G p the nodal gradient term, see the top of the file. The
// velocity and g have comps components; the compact gradient is always 3D (in
// planar flow p does not depend on z, so it adds nothing across z).
template<typename KeyType, typename RealType>
__global__ void nsFaceFluxKernel(HexElements<KeyType> hex, const RealType* x, const RealType* y, const RealType* z,
                                 Components<const RealType*> area, int comps, Components<const RealType*> u,
                                 const RealType* p, Components<const RealType*> g, RealType h, RealType* flux)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    KeyType n[8];
    double corner[8][3];
    nsCorners(hex, hex.first + k, x, y, z, n, corner);
    for (int ip = 0; ip < 12; ++ip)
    {
        size_t f = k * 12 + ip;
        double w[8];
        nsFaceWeights(ip, corner, area, f, w);
        KeyType L = n[hexLRSCV[2 * ip]], R = n[hexLRSCV[2 * ip + 1]];
        double average = 0, nodal = 0, compact = 0;
        for (int d = 0; d < comps; ++d)
        {
            average += 0.5 * (u.c[d][L] + u.c[d][R]) * area.c[d][f];
            nodal += 0.5 * (g.c[d][L] + g.c[d][R]) * area.c[d][f];
        }
        for (int m = 0; m < 8; ++m)
            compact += w[m] * p[n[m]];
        flux[f] = RealType(average - h * (compact + nodal));
    }
}

// flux_f += c (grad phi . A)_f
template<typename KeyType, typename RealType>
__global__ void nsFaceGradientKernel(HexElements<KeyType> hex, const RealType* x, const RealType* y,
                                     const RealType* z, Components<const RealType*> area, const RealType* phi,
                                     RealType c, RealType* flux)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    KeyType n[8];
    double corner[8][3];
    nsCorners(hex, hex.first + k, x, y, z, n, corner);
    for (int ip = 0; ip < 12; ++ip)
    {
        size_t f = k * 12 + ip;
        double w[8], gradient = 0;
        nsFaceWeights(ip, corner, area, f, w);
        for (int m = 0; m < 8; ++m)
            gradient += w[m] * phi[n[m]];
        flux[f] += RealType(c * gradient);
    }
}

// D_F F: the flux through a face leaves node L and enters node R.
template<typename KeyType, typename RealType>
__global__ void nsFaceDivergenceKernel(HexElements<KeyType> hex, const RealType* flux, RealType* div)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
#pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType L = hex.node[hexLRSCV[2 * ip]][e], R = hex.node[hexLRSCV[2 * ip + 1]][e];
        atomicAdd(&div[L], flux[k * 12 + ip]);
        atomicAdd(&div[R], -flux[k * 12 + ip]);
    }
}


// Skew-symmetric advection by the face fluxes: node L receives -F q_R / 2 and
// node R receives +F q_L / 2. Its kinetic energy production, sum_i q_i (N q)_i,
// cancels face by face.
template<typename KeyType, typename RealType>
__global__ void nsAdvectionKernel(HexElements<KeyType> hex, const RealType* flux, int comps,
                                  Components<const RealType*> u, Components<RealType*> adv)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
#pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType L = hex.node[hexLRSCV[2 * ip]][e], R = hex.node[hexLRSCV[2 * ip + 1]][e];
        RealType F = flux[k * 12 + ip];
        for (int d = 0; d < comps; ++d)
        {
            atomicAdd(&adv.c[d][L], -RealType(0.5) * F * u.c[d][R]);
            atomicAdd(&adv.c[d][R], RealType(0.5) * F * u.c[d][L]);
        }
    }
}

// D^T p, the transpose of the scatter above: (p_L - p_R) / 2 A_f to both ends of the face.
template<typename KeyType, typename RealType>
__global__ void nsDivergenceTransposeKernel(HexElements<KeyType> hex, Components<const RealType*> area, int comps,
                                            const RealType* p, Components<RealType*> g)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n[8];
    for (int c = 0; c < 8; ++c)
        n[c] = hex.node[c][e];
#pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType L = n[hexLRSCV[2 * ip]], R = n[hexLRSCV[2 * ip + 1]];
        size_t f    = k * 12 + ip;
        RealType dp = RealType(0.5) * (p[L] - p[R]);
        for (int d = 0; d < comps; ++d)
        {
            atomicAdd(&g.c[d][L], dp * area.c[d][f]);
            atomicAdd(&g.c[d][R], dp * area.c[d][f]);
        }
    }
}

// Flux out through the opening faces of node i: the prescribed velocity where
// it is fixed (inflow), the computed one elsewhere (outflow).
template<typename RealType>
__device__ inline RealType nsOpeningFlux(size_t i, int comps, Components<const RealType*> area, const uint8_t* fixed,
                                         Components<const RealType*> target, Components<const RealType*> u)
{
    RealType flux = 0;
    for (int d = 0; d < comps; ++d)
        flux += area.c[d][i] * (fixed[i] ? target.c[d][i] : u.c[d][i]);
    return flux;
}

template<typename RealType>
__global__ void nsAddOpeningFluxKernel(size_t n, const uint8_t* isDof, int comps, Components<const RealType*> area,
                                       const uint8_t* fixed, Components<const RealType*> target,
                                       Components<const RealType*> u, RealType* div)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n && isDof[i]) div[i] += nsOpeningFlux(i, comps, area, fixed, target, u);
}

// The skew form of an opening face with flux m is -m q / 2 at its node.
template<typename RealType>
__global__ void nsOpeningAdvectionKernel(size_t n, const uint8_t* isDof, int comps, Components<const RealType*> area,
                                         const uint8_t* fixed, Components<const RealType*> target,
                                         Components<const RealType*> u, Components<RealType*> adv)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n || !isDof[i]) return;
    RealType m = nsOpeningFlux(i, comps, area, fixed, target, u);
    for (int d = 0; d < comps; ++d)
        adv.c[d][i] -= RealType(0.5) * u.c[d][i] * m;
}

// g = Q M^-1 g at the DOF slots, 0 elsewhere. withQ == false keeps fixed nodes (gradients for output).
template<typename RealType>
__global__ void nsInverseMassKernel(size_t n, const uint8_t* isDof, const uint8_t* fixed, int withQ,
                                    const RealType* massDof, int comps, Components<RealType*> g)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    bool keep  = isDof[i] && !(withQ && fixed[i]);
    RealType s = keep ? RealType(1) / massDof[i] : RealType(0);
    for (int d = 0; d < comps; ++d)
        g.c[d][i] = keep ? g.c[d][i] * s : RealType(0);
}

// u* from u^n (and u^{n-1} for BDF2), the advection a = N(u^n) (extrapolated to
// 2 a - a_prev for BDF2), and g = Q M^-1 D^T p^n ~ -grad p^n.
template<typename RealType>
__global__ void nsPredictorKernel(size_t n, const uint8_t* isDof, const uint8_t* fixed, int bdf2, RealType dt,
                                  RealType invRho, const RealType* massDof, int comps,
                                  Components<const RealType*> target, Components<const RealType*> u,
                                  Components<const RealType*> um1, Components<const RealType*> adv,
                                  Components<const RealType*> advm1, Components<const RealType*> g,
                                  Components<RealType*> out)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n || !isDof[i]) return;
    RealType V = massDof[i];
    for (int d = 0; d < comps; ++d)
    {
        if (fixed[i]) out.c[d][i] = target.c[d][i];
        else if (bdf2)
        {
            RealType c = RealType(2) * dt / RealType(3);
            out.c[d][i] = RealType(4) / RealType(3) * u.c[d][i] - RealType(1) / RealType(3) * um1.c[d][i] +
                          c * ((RealType(2) * adv.c[d][i] - advm1.c[d][i]) / V + invRho * g.c[d][i]);
        }
        else out.c[d][i] = u.c[d][i] + dt * adv.c[d][i] / V + dt * invRho * g.c[d][i];
    }
}

// Right-hand side M q* / dtEff + lift and initial guess q* of one velocity
// component, per DOF. A fixed DOF has an identity row: its right-hand side is
// the prescribed value.
template<typename RealType>
__global__ void nsViscousRhsKernel(size_t n, const uint8_t* isDof, const int* dofIndex, const uint8_t* fixed,
                                   const RealType* massDof, RealType invDt, const RealType* qStar,
                                   const RealType* lift, const RealType* target, RealType* rhs, RealType* x)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n || !isDof[i]) return;
    int dof  = dofIndex[i];
    rhs[dof] = fixed[i] ? target[i] : massDof[i] * invDt * qStar[i] + lift[i];
    x[dof]   = qStar[i];
}

// K phi = -(rho / dtEff) (D F** + openings) per DOF; a p = 0 row is an identity row with value 0.
template<typename RealType>
__global__ void nsPressureRhsKernel(size_t n, const uint8_t* isDof, const int* dofIndex, const uint8_t* pressureFixed,
                                    RealType coef, const RealType* div, RealType* rhs)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n && isDof[i]) rhs[dofIndex[i]] = pressureFixed[i] ? RealType(0) : -coef * div[i];
}

template<typename RealType>
__global__ void nsFromDofKernel(size_t n, const uint8_t* isDof, const int* dofIndex, const RealType* x, RealType* q)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n && isDof[i]) q[i] = x[dofIndex[i]];
}

// u = u** + s g where free, the prescribed velocity where fixed, and p += phi.
template<typename RealType>
__global__ void nsCorrectorKernel(size_t n, const uint8_t* isDof, const uint8_t* fixed, RealType s, int comps,
                                  Components<const RealType*> target, Components<const RealType*> uss,
                                  Components<const RealType*> g, const RealType* phi, Components<RealType*> u,
                                  RealType* p)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n || !isDof[i]) return;
    for (int d = 0; d < comps; ++d)
        u.c[d][i] = fixed[i] ? target.c[d][i] : uss.c[d][i] + s * g.c[d][i];
    p[i] += phi[i];
}

template<typename RealType>
__global__ void nsApplyTargetKernel(size_t n, const uint8_t* fixed, int comps, Components<const RealType*> target,
                                    Components<RealType*> u)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n || !fixed[i]) return;
    for (int d = 0; d < comps; ++d)
        u.c[d][i] = target.c[d][i];
}

// ---------------------------------------------------------------------------
// Functors for thrust. They live at namespace scope: nvcc rejects private
// class members as template arguments of the kernels thrust instantiates.
// ---------------------------------------------------------------------------

template<typename RealType>
struct NsScaleBy
{
    RealType s;
    __device__ RealType operator()(RealType a) const { return a * s; }
};

struct NsToInt
{
    __device__ int operator()(uint8_t b) const { return b; }
};

template<typename RealType>
struct NsNegative
{
    __device__ bool operator()(RealType g) const { return g < 0; }
};

template<typename RealType>
struct NsFlaggedDof
{
    const uint8_t* isDof;
    const uint8_t* flag;
    __device__ bool operator()(size_t i) const { return isDof[i] && flag[i]; }
};

template<typename RealType>
struct NsGidOf
{
    const RealType* gid;
    __device__ HYPRE_BigInt operator()(size_t i) const { return HYPRE_BigInt(gid[i]); }
};

template<typename RealType>
struct NsMassNormTerm
{
    const uint8_t* isDof;
    const RealType* mass;
    const RealType* q;
    __device__ double operator()(size_t i) const { return isDof[i] ? double(mass[i]) * q[i] * q[i] : 0.0; }
};

template<typename RealType>
struct NsKineticTerm
{
    const uint8_t* isDof;
    const RealType* mass;
    int comps;
    Components<const RealType*> u;
    __device__ double operator()(size_t i) const
    {
        if (!isDof[i]) return 0.0;
        double s = 0;
        for (int d = 0; d < comps; ++d)
            s += double(u.c[d][i]) * u.c[d][i];
        return 0.5 * mass[i] * s;
    }
};

// |D_F F + openings| / M where the projection enforces continuity: DOFs without p = 0.
template<typename RealType>
struct NsContinuityTerm
{
    const uint8_t* isDof;
    const uint8_t* pressureFixed;
    const RealType* mass;
    const RealType* div;
    __device__ double operator()(size_t i) const
    {
        return isDof[i] && !pressureFixed[i] ? fabs(double(div[i]) / mass[i]) : 0.0;
    }
};

// (|assembled - matrix-free|, |matrix-free|) per DOF without p = 0.
template<typename RealType>
struct NsOperatorDifference
{
    const uint8_t* isDof;
    const uint8_t* pressureFixed;
    const int* dofIndex;
    const RealType* matFree;
    const RealType* assembled;
    __device__ thrust::tuple<double, double> operator()(size_t i) const
    {
        if (!isDof[i] || pressureFixed[i]) return thrust::make_tuple(0.0, 0.0);
        double ref = matFree[i];
        return thrust::make_tuple(fabs(assembled[dofIndex[i]] - ref), fabs(ref));
    }
};

struct NsMaxPair
{
    __device__ thrust::tuple<double, double> operator()(const thrust::tuple<double, double>& a,
                                                        const thrust::tuple<double, double>& b) const
    {
        return thrust::make_tuple(fmax(thrust::get<0>(a), thrust::get<0>(b)),
                                  fmax(thrust::get<1>(a), thrust::get<1>(b)));
    }
};

// Sums of the projection check, see NavierStokes::projectionReport.
struct NsProjectionSums
{
    double sum[7]{};
    double max[2]{};
    __host__ __device__ NsProjectionSums operator+(const NsProjectionSums& o) const
    {
        NsProjectionSums r;
        for (int k = 0; k < 7; ++k)
            r.sum[k] = sum[k] + o.sum[k];
        for (int k = 0; k < 2; ++k)
            r.max[k] = max[k] > o.max[k] ? max[k] : o.max[k];
        return r;
    }
};

template<typename RealType>
struct NsProjectionTerm
{
    const uint8_t* isDof;
    const uint8_t* pressureFixed;
    const uint8_t* fixed;
    const RealType* mass;
    const RealType* before; // D_F F** + openings
    const RealType* after;  // D_F F + openings
    const RealType* action; // K phi
    int comps;
    Components<const RealType*> area;
    Components<const RealType*> target;
    Components<const RealType*> u;
    RealType h; // dtEff / rho

    __device__ NsProjectionSums operator()(size_t i) const
    {
        NsProjectionSums s;
        if (!isDof[i]) return s;
        double V = mass[i], a = after[i], b = before[i], ap = action[i];
        if (!(V > 0) || !isfinite(V) || !isfinite(a) || !isfinite(b) || !isfinite(ap))
        {
            s.max[1] = 1;
            return s;
        }
        if (!pressureFixed[i])
        {
            double error = a - b - double(h) * ap;
            s.sum[0]     = b * b / V;
            s.sum[1]     = a * a / V;
            s.sum[2]     = error * error / V;
            s.sum[3]     = V;
            s.max[0]     = fabs(a) / V;
        }
        double flux = 0;
        for (int d = 0; d < comps; ++d)
            flux += double(area.c[d][i]) * (fixed[i] ? target.c[d][i] : u.c[d][i]);
        (fixed[i] ? s.sum[4] : s.sum[5]) = flux;
        s.sum[6] = a;
        return s;
    }
};

template<typename KeyType, typename RealType>
class NavierStokes
{
    static_assert(std::is_same_v<RealType, HYPRE_Complex>, "Hypre is built for another floating point type");

public:
    using Domain = ElementDomain<HexTag, RealType, KeyType, cstone::execution::Gpu>;
    using Vector = cstone::DeviceVector<RealType>;
    using Space  = DofSpace<KeyType, RealType, Domain>;
    using Map    = PeriodicMap<KeyType, RealType>;

    struct Params
    {
        RealType rho       = 1;
        RealType nu        = RealType(0.01);
        RealType dt        = RealType(0.01);
        bool planar        = false; // flow in x-y on one layer of elements: w = 0, two unknowns per node
        bool bdf2          = true;
        RealType tolerance = RealType(1e-10); // relative PCG tolerance of both solves
        int maxIter        = 1000;
        int blockSize      = 256;
    };

    // Milliseconds per stage, summed over the steps, and the summed AMG-PCG iterations.
    struct Timing
    {
        double predictor = 0, viscous = 0, pressure = 0, corrector = 0;
        double hypre = 0; // inside the Hypre solves, part of viscous and pressure
        long velocityIterations = 0, pressureIterations = 0;
    };

    // The last projection, checked. continuity = D_F F + openings per unit volume.
    struct ProjectionReport
    {
        double identity;        // |D_F F - D_F F** - h K phi| / |D_F F**|, at the solver tolerance
        double identityRms;     // the same, per unit volume
        double continuityRms;   // RMS of continuity over the rows the projection enforces
        double continuityMax;   // its maximum
        double inflow;          // flux through openings with prescribed velocity
        double outflow;         // flux through openings with computed velocity
        double balance;         // |inflow + outflow| / |inflow|
        double balanceIdentity; // total continuity against inflow + outflow
        bool finite;
    };

    // State per node slot. w stays 0 in planar mode.
    Vector u, v, w, p;

    // rule(x, y, z) gives the NodeCondition of every node; openings are the
    // inlets and outlets; periodic pairs the periodic nodes (nullptr: none).
    template<class Rule>
    NavierStokes(Domain& domain, const Params& params, Rule rule, const std::vector<Opening<RealType>>& openings,
                 const Map* periodic = nullptr)
        : domain_(domain)
        , prm_(params)
        , n_(domain.getNodeCount())
        , rank_(domain.rank())
        , comps_(params.planar ? 2 : 3)
        , space_(domain, periodic, params.blockSize)
    {
        domain_.cacheNodeCoordinates();
        const auto& conn = domain_.getElementToNodeConnectivity();
        hex_ = {{std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                 std::get<3>(conn).data(), std::get<4>(conn).data(), std::get<5>(conn).data(),
                 std::get<6>(conn).data(), std::get<7>(conn).data()},
                domain_.startIndex(),
                domain_.localElementCount()};
        box_ = boundingBox(domain_);

        for (Vector* f :
             {&u, &v, &w, &p, &phi_, &div_, &divStar_, &massLocal_, &massDof_, &gid_, &work_[0], &work_[1], &work_[2]})
            allocate(*f);
        for (int d = 0; d < 3; ++d)
            for (Vector* f : {&um1_[d], &star_[d], &sstar_[d], &adv_[d], &advm1_[d], &g_[d], &target_[d],
                              &openArea_[d], &lift_[d]})
                allocate(*f);
        flux_.resize(12 * hex_.count);

        if (prm_.planar) checkPlanarMesh();
        buildGeometry();
        applyConditions(rule);
        numberDofs();
        buildOpenings(openings);
        buildViscousSolvers();
        buildPressureSolver();
    }

    // After setting u, v, w, p: prescribed velocities, periodic copies and ghosts,
    // and the face fluxes of u; the next step is BDF1.
    void start()
    {
        nsApplyTargetKernel<RealType><<<grid(), bs()>>>(n_, fixed_.data(), comps_, cview(target_), view(vel()));
        cudaCheckError();
        for (Vector* q : velocity())
            space_.prolong(*q);
        space_.prolong(p);
        if (hex_.count > 0)
            nsFaceFluxKernel<KeyType, RealType><<<elemGrid(), bs()>>>(hex_, x(), y(), z(), cview(area_), comps_,
                                                                      cview(vel()), p.data(), cview(g_), RealType(0),
                                                                      flux_.data());
        cudaCheckError();
        steps_ = 0;
    }

    // One time step. Returns false (on every rank) if a linear solve did not converge.
    bool step()
    {
        const bool bdf2 = prm_.bdf2 && steps_ > 0;
        StageClock clock;
        predict(bdf2);
        timing_.predictor += clock.lap();
        if (!solveViscous(bdf2)) return false;
        timing_.viscous += clock.lap();
        if (!project(bdf2)) return false;
        timing_.pressure += clock.lap();
        correct(bdf2);
        timing_.corrector += clock.lap();
        lastBdf2_ = bdf2;
        ++steps_;
        return true;
    }

    // ||q||_M = sqrt(sum_i M_i q_i^2) over the DOFs.
    double norm(const Vector& q) const
    {
        return std::sqrt(dofSum(NsMassNormTerm<RealType>{space_.isDof(), massDof_.data(), q.data()}));
    }

    double kineticEnergy() const
    {
        return dofSum(NsKineticTerm<RealType>{space_.isDof(), massDof_.data(), comps_, cview(vel())});
    }

    // max |D_F F + openings| / M over the rows the projection enforces.
    double maxContinuity()
    {
        fluxDivergence(div_, vel());
        double local = thrust::transform_reduce(
            thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n_),
            NsContinuityTerm<RealType>{space_.isDof(), pressureFixed_.data(), massDof_.data(), div_.data()}, 0.0,
            thrust::maximum<double>());
        double global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return global;
    }

    // Checks the last step: the corrected fluxes must satisfy
    // D_F F = D_F F** + h K phi (h = dtEff / rho) with the matrix the pressure solve inverted.
    ProjectionReport projectionReport()
    {
        Vector& after  = work_[1];
        Vector& action = work_[2];
        fluxDivergence(after, vel());
        applyLaplacian(phi_, action);
        const RealType dtEff = lastBdf2_ ? RealType(2) * prm_.dt / RealType(3) : prm_.dt;
        NsProjectionTerm<RealType> term{space_.isDof(),   pressureFixed_.data(), fixed_.data(),
                                        massDof_.data(),  divStar_.data(),       after.data(),
                                        action.data(),    comps_,                cview(openArea_),
                                        cview(target_),   cview(vel()),          dtEff / prm_.rho};
        NsProjectionSums local = thrust::transform_reduce(
            thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n_), term,
            NsProjectionSums{}, thrust::plus<NsProjectionSums>());
        NsProjectionSums g;
        MPI_Allreduce(local.sum, g.sum, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(local.max, g.max, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        ProjectionReport r;
        double scale      = std::max(g.sum[0], 1e-60);
        double volume     = std::max(g.sum[3], 1e-60);
        r.identity        = std::sqrt(g.sum[2] / scale);
        r.identityRms     = std::sqrt(g.sum[2] / volume);
        r.continuityRms   = std::sqrt(g.sum[1] / volume);
        r.continuityMax   = g.max[0];
        r.inflow          = g.sum[4];
        r.outflow         = g.sum[5];
        r.balance         = std::abs(g.sum[4] + g.sum[5]) / std::max(std::abs(g.sum[4]), 1e-30);
        r.balanceIdentity = std::abs(g.sum[6] - g.sum[4] - g.sum[5]) /
                            std::max(std::abs(g.sum[4]) + std::abs(g.sum[5]), 1e-30);
        r.finite          = g.max[1] == 0 && std::isfinite(r.identity) && std::isfinite(r.continuityRms);
        return r;
    }

    // grad q per node slot (for output, e.g. vorticity). q must be current on every slot.
    void gradient(const Vector& q, Vector& gx, Vector& gy, Vector& gz)
    {
        gradientOf(q, false);
        Vector* out[3] = {&gx, &gy, &gz};
        for (int d = 0; d < 3; ++d)
        {
            allocate(*out[d]);
            if (d < comps_)
                thrust::transform(thrust::device, g_[d].data(), g_[d].data() + n_, out[d]->data(),
                                  NsScaleBy<RealType>{RealType(-1)});
            space_.prolong(*out[d]);
        }
    }

    const Vector& mass() const { return massDof_; } // lumped mass of each DOF, on every copy
    const Space& space() const { return space_; }
    const Box<RealType>& box() const { return box_; }
    const Timing& timing() const { return timing_; }
    void resetTiming() { timing_ = {}; }

    // Stage times of the slowest rank over the steps since the last resetTiming (collective).
    void printTiming(long steps) const
    {
        steps         = std::max(steps, 1L);
        const auto& t = timing_;
        double stages = t.predictor + t.viscous + t.pressure + t.corrector;
        double local[6] = {t.predictor, t.viscous, t.pressure, t.corrector, t.hypre, stages - t.hypre}, slowest[6];
        MPI_Allreduce(local, slowest, 6, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (rank_ != 0) return;
        int ranks = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &ranks);
        double total = slowest[0] + slowest[1] + slowest[2] + slowest[3];
        long long n  = globalDofs();
        std::cout << std::fixed << std::setprecision(3) << "[timing] ranks=" << ranks << " nodes=" << n
                  << " nodes/rank=" << n / ranks << " steps=" << steps << " ms/step: total=" << total / steps
                  << " predictor=" << slowest[0] / steps << " viscous=" << slowest[1] / steps
                  << " pressure=" << slowest[2] / steps << " corrector=" << slowest[3] / steps
                  << " | pressure_it/step=" << double(timing_.pressureIterations) / steps
                  << " ms/pressure_it="
                  << (timing_.pressureIterations > 0 ? slowest[2] / timing_.pressureIterations : 0.0)
                  << " velocity_it/step=" << double(timing_.velocityIterations) / steps
                  << " | hypre=" << slowest[4] / steps << " mars=" << slowest[5] / steps << "\n"
                  << std::defaultfloat;
    }
    int velocityIterations(int component) const { return velocityIters_[component]; }
    int pressureIterations() const { return pressureIters_; }
    long long globalDofs() const { return space_.numDofs(); }
    int components() const { return comps_; }

private:
    struct StageClock
    {
        double t;
        StageClock()
        {
            cudaDeviceSynchronize();
            t = MPI_Wtime();
        }
        double lap()
        {
            cudaDeviceSynchronize();
            double now = MPI_Wtime(), ms = 1e3 * (now - t);
            t          = now;
            return ms;
        }
    };

    std::array<Vector*, 3> velocity() { return {&u, &v, &w}; }
    Components<RealType*> vel() { return {{u.data(), v.data(), w.data()}}; }
    Components<const RealType*> vel() const { return {{u.data(), v.data(), w.data()}}; }
    static Components<RealType*> view(Components<RealType*> c) { return c; }
    static Components<RealType*> view(Vector (&f)[3]) { return {{f[0].data(), f[1].data(), f[2].data()}}; }
    static Components<const RealType*> cview(const Vector (&f)[3])
    {
        return {{f[0].data(), f[1].data(), f[2].data()}};
    }
    static Components<const RealType*> cview(Components<RealType*> c) { return {{c.c[0], c.c[1], c.c[2]}}; }
    static Components<const RealType*> cview(Components<const RealType*> c) { return c; }
    const RealType* x() const { return domain_.getNodeX().data(); }
    const RealType* y() const { return domain_.getNodeY().data(); }
    const RealType* z() const { return domain_.getNodeZ().data(); }

    // 1. u* from the history, the advection by F^n and the old pressure.
    void predict(bool bdf2)
    {
        scatterGradient(p);
        for (int d = 0; d < comps_; ++d)
            zero(adv_[d]);
        if (hex_.count > 0)
            nsAdvectionKernel<KeyType, RealType>
                <<<elemGrid(), bs()>>>(hex_, flux_.data(), comps_, cview(vel()), view(adv_));
        cudaCheckError();
        Vector* scattered[6];
        for (int d = 0; d < comps_; ++d)
        {
            scattered[d]          = &g_[d];
            scattered[comps_ + d] = &adv_[d];
        }
        space_.restrict(scattered, 2 * comps_);
        inverseMass(true);
        nsOpeningAdvectionKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), comps_, cview(openArea_),
                                                             fixed_.data(), cview(target_), cview(vel()),
                                                             view(adv_));
        nsPredictorKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), fixed_.data(), int(bdf2), prm_.dt,
                                                      RealType(1) / prm_.rho, massDof_.data(), comps_,
                                                      cview(target_), cview(vel()), cview(um1_), cview(adv_),
                                                      cview(advm1_), cview(g_), view(star_));
        cudaCheckError();
        for (int d = 0; d < comps_; ++d)
            adv_[d].swap(advm1_[d]);
    }

    // 2. (M / dtEff + nu K) u** = M u* / dtEff + lift per component, warm-started from u*.
    bool solveViscous(bool bdf2)
    {
        const RealType invDt      = bdf2 ? RealType(3) / (RealType(2) * prm_.dt) : RealType(1) / prm_.dt;
        HypreAmgPcgSolver& solver = bdf2 ? *viscousBdf2_ : *viscousBdf1_;
        for (int d = 0; d < comps_; ++d)
        {
            nsViscousRhsKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), dofIndex_.data(), fixed_.data(),
                                                           massDof_.data(), invDt, star_[d].data(), lift_[d].data(),
                                                           target_[d].data(), rhs_.data(), x_.data());
            cudaCheckError();
            StageClock clock;
            velocityIters_[d] = solver.solve(rhs_.data(), x_.data(), true);
            timing_.hypre += clock.lap();
            if (velocityIters_[d] < 0) return failed("velocity", solver);
            timing_.velocityIterations += velocityIters_[d];
            nsFromDofKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), dofIndex_.data(), x_.data(),
                                                        sstar_[d].data());
            cudaCheckError();
        }
        space_.prolong(fieldsOf(sstar_).data(), comps_);
        return true;
    }

    // 3. The stabilized fluxes F** of u** and p^n, then K phi = -(rho / dtEff) (D_F F** + openings).
    bool project(bool bdf2)
    {
        const RealType invDt = bdf2 ? RealType(3) / (RealType(2) * prm_.dt) : RealType(1) / prm_.dt;
        gradientOf(p, false);
        space_.prolong(fieldsOf(g_).data(), comps_);
        if (hex_.count > 0)
            nsFaceFluxKernel<KeyType, RealType><<<elemGrid(), bs()>>>(hex_, x(), y(), z(), cview(area_), comps_,
                                                                      cview(sstar_), p.data(), cview(g_),
                                                                      RealType(1) / (prm_.rho * invDt),
                                                                      flux_.data());
        cudaCheckError();
        fluxDivergence(divStar_, view(sstar_));
        nsPressureRhsKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), dofIndex_.data(), pressureFixed_.data(),
                                                        prm_.rho * invDt, divStar_.data(), rhs_.data());
        cudaCheckError();
        StageClock clock;
        pressureIters_ = pressure_->solve(rhs_.data(), x_.data());
        timing_.hypre += clock.lap();
        if (pressureIters_ < 0) return failed("pressure", *pressure_);
        timing_.pressureIterations += pressureIters_;
        nsFromDofKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), dofIndex_.data(), x_.data(), phi_.data());
        cudaCheckError();
        space_.prolong(phi_);
        return true;
    }

    // 4. F = F** - h (grad phi . A), u = u** - h Q G phi, p += phi. u^n becomes the BDF2 history.
    void correct(bool bdf2)
    {
        const RealType h = (bdf2 ? RealType(2) * prm_.dt / RealType(3) : prm_.dt) / prm_.rho;
        if (hex_.count > 0)
            nsFaceGradientKernel<KeyType, RealType>
                <<<elemGrid(), bs()>>>(hex_, x(), y(), z(), cview(area_), phi_.data(), -h, flux_.data());
        cudaCheckError();
        for (int d = 0; d < comps_; ++d)
            velocity()[d]->swap(um1_[d]);
        gradientOf(phi_, true);
        nsCorrectorKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), fixed_.data(), h, comps_, cview(target_),
                                                      cview(sstar_), cview(g_), phi_.data(), vel(), p.data());
        cudaCheckError();
        Vector* fields[4] = {velocity()[0], velocity()[1], velocity()[2], &p};
        fields[comps_]    = &p;
        space_.prolong(fields, comps_ + 1);
    }

    bool failed(const char* system, const HypreAmgPcgSolver& solver) const
    {
        if (rank_ == 0)
            std::cerr << "NavierStokes: the " << system << " solve did not converge: relative residual "
                      << std::scientific << solver.lastRelativeResidual() << std::defaultfloat << " after "
                      << solver.lastIterations() << " iterations\n";
        return false;
    }

    // g = M^-1 D^T q (= -G q, ~ -grad q) at the DOF slots, with Q if withQ. q must be current on every slot.
    void gradientOf(const Vector& q, bool withQ)
    {
        scatterGradient(q);
        space_.restrict(fieldsOf(g_).data(), comps_);
        inverseMass(withQ);
    }

    // g = D^T q over this rank's elements, before restrict.
    void scatterGradient(const Vector& q)
    {
        for (int d = 0; d < comps_; ++d)
            zero(g_[d]);
        if (hex_.count > 0)
            nsDivergenceTransposeKernel<KeyType, RealType>
                <<<elemGrid(), bs()>>>(hex_, cview(area_), comps_, q.data(), view(g_));
        cudaCheckError();
    }

    void inverseMass(bool withQ)
    {
        nsInverseMassKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), fixed_.data(), int(withQ),
                                                        massDof_.data(), comps_, view(g_));
        cudaCheckError();
    }

    static std::array<Vector*, 3> fieldsOf(Vector (&c)[3]) { return {&c[0], &c[1], &c[2]}; }

    // out = D_F F + openings at the DOF slots; the outlets carry the nodal velocity.
    void fluxDivergence(Vector& out, Components<RealType*> velocityField)
    {
        zero(out);
        if (hex_.count > 0)
            nsFaceDivergenceKernel<KeyType, RealType><<<elemGrid(), bs()>>>(hex_, flux_.data(), out.data());
        cudaCheckError();
        space_.restrict(out);
        nsAddOpeningFluxKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), comps_, cview(openArea_),
                                                           fixed_.data(), cview(target_), cview(velocityField),
                                                           out.data());
        cudaCheckError();
    }

    // out = K q at the DOF slots, matrix-free: -D_F of the compact face gradients.
    // q must be current on every slot.
    void applyLaplacian(const Vector& q, Vector& out)
    {
        Vector& faces = faceWork_;
        faces.resize(12 * hex_.count);
        zero(faces);
        zero(out);
        if (hex_.count > 0)
        {
            nsFaceGradientKernel<KeyType, RealType>
                <<<elemGrid(), bs()>>>(hex_, x(), y(), z(), cview(area_), q.data(), RealType(-1), faces.data());
            nsFaceDivergenceKernel<KeyType, RealType><<<elemGrid(), bs()>>>(hex_, faces.data(), out.data());
        }
        cudaCheckError();
        space_.restrict(out);
    }

    void checkPlanarMesh()
    {
        // One layer of elements between two z planes, so w = 0 and nothing depends on z.
        const RealType eps = RealType(1e-8) * box_.extent();
        long long bad =
            thrust::count_if(thrust::device, z(), z() + n_, NsOffPlanes<RealType>{box_.lo[2], box_.hi[2], eps});
        MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (bad > 0 || !(box_.hi[2] > box_.lo[2]))
        {
            if (rank_ == 0)
                std::cerr << "NavierStokes: planar flow needs one layer of elements between two z planes\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Sub-control face area vectors of the local elements, and the lumped mass:
    // local sub-volumes per slot (for assembly) and the DOF mass on every copy.
    void buildGeometry()
    {
        for (int d = 0; d < 3; ++d)
            area_[d].resize(12 * hex_.count);
        if (hex_.count > 0)
        {
            const KeyType* c[8];
            for (int i = 0; i < 8; ++i)
                c[i] = hex_.node[i] + hex_.first;
            precomputeAreaVectorsGpu<KeyType, RealType>(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], hex_.count,
                                                        x(), y(), z(), area_[0].data(), area_[1].data(),
                                                        area_[2].data());
            cudaCheckError();
            nsMassKernel<KeyType, RealType><<<elemGrid(), bs()>>>(hex_, x(), y(), z(), massLocal_.data());
        }
        cudaCheckError();
        copy(massLocal_, massDof_);
        space_.restrict(massDof_);
        space_.prolong(massDof_);
    }

    template<class Rule>
    void applyConditions(Rule rule)
    {
        Vector& flags = work_[0];
        nsBoundaryKernel<RealType, Rule>
            <<<grid(), bs()>>>(n_, x(), y(), z(), rule, flags.data(), view(target_));
        cudaCheckError();
        // The DOF slot decides; every copy of the node takes its conditions.
        space_.prolong(flags);
        for (int d = 0; d < 3; ++d)
            space_.prolong(target_[d]);
        fixed_.resize(n_);
        pressureFixed_.resize(n_);
        nsUnpackFlagsKernel<RealType><<<grid(), bs()>>>(n_, flags.data(), fixed_.data(), pressureFixed_.data());
        cudaCheckError();
    }

    // DOF slots get consecutive local indices; rank r owns the global ids
    // [dofStart_, dofStart_ + numOwned_). Every slot learns the global id of its
    // DOF. With no p = 0 anywhere, the DOF with id 0 gets p = 0.
    void numberDofs()
    {
        dofIndex_.resize(n_);
        thrust::transform(thrust::device, space_.isDof(), space_.isDof() + n_, dofIndex_.data(), NsToInt{});
        thrust::exclusive_scan(thrust::device, dofIndex_.data(), dofIndex_.data() + n_, dofIndex_.data());
        numOwned_ = int(thrust::count(thrust::device, space_.isDof(), space_.isDof() + n_, uint8_t(1)));

        long long owned = numOwned_, slots = (long long)n_;
        MPI_Exscan(&owned, &dofStart_, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Exscan(&slots, &localStart_, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (rank_ == 0) dofStart_ = localStart_ = 0;

        nsGlobalIdKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), dofIndex_.data(), dofStart_, gid_.data());
        cudaCheckError();
        space_.prolong(gid_);
        long long missing = thrust::count_if(thrust::device, gid_.data(), gid_.data() + n_, NsNegative<RealType>{});
        MPI_Allreduce(MPI_IN_PLACE, &missing, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (missing > 0)
        {
            if (rank_ == 0) std::cerr << "NavierStokes: " << missing << " node copies have no DOF\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        rhs_.resize(numOwned_);
        x_.resize(numOwned_);

        long long fixedPressure = thrust::count_if(thrust::device, thrust::counting_iterator<size_t>(0),
                                                   thrust::counting_iterator<size_t>(n_),
                                                   NsFlaggedDof<RealType>{space_.isDof(), pressureFixed_.data()});
        MPI_Allreduce(MPI_IN_PLACE, &fixedPressure, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (fixedPressure == 0)
        {
            Vector& flag = work_[0];
            nsPinFirstDofKernel<RealType>
                <<<grid(), bs()>>>(n_, space_.isDof(), gid_.data(), pressureFixed_.data(), flag.data());
            cudaCheckError();
            space_.prolong(flag);
            nsFlagFromRealKernel<RealType><<<grid(), bs()>>>(n_, flag.data(), pressureFixed_.data());
            cudaCheckError();
        }
    }

    void buildOpenings(const std::vector<Opening<RealType>>& openings)
    {
        const RealType tolerance = RealType(1e-6) * box_.extent();
        for (const auto& opening : openings)
            if (hex_.count > 0)
                nsOpeningAreaKernel<KeyType, RealType><<<elemGrid(), bs()>>>(hex_, x(), y(), z(), opening, tolerance,
                                                                             openArea_[opening.axis].data());
        cudaCheckError();
        for (int d = 0; d < 3; ++d)
        {
            space_.restrict(openArea_[d]);
            space_.prolong(openArea_[d]);
        }
    }

    // P: local slot i (all ranks' slots numbered consecutively) takes the value of
    // its DOF. comps > 1: one copy per velocity component.
    HypreMatrix prolongation(int comps)
    {
        const size_t count = n_ * comps;
        thrust::device_vector<HYPRE_BigInt> rows(count), cols(count);
        thrust::device_vector<RealType> values(count);
        nsProlongationTripletsKernel<RealType><<<grid(), bs()>>>(
            n_, comps, localStart_, gid_.data(), thrust::raw_pointer_cast(rows.data()),
            thrust::raw_pointer_cast(cols.data()), thrust::raw_pointer_cast(values.data()));
        cudaCheckError();
        return hypreAssemble(MPI_COMM_WORLD, comps * localStart_, comps * (localStart_ + (long long)n_),
                             comps * dofStart_, comps * (dofStart_ + numOwned_),
                             HypreCoo{thrust::raw_pointer_cast(rows.data()), thrust::raw_pointer_cast(cols.data()),
                                      thrust::raw_pointer_cast(values.data()), HYPRE_Int(count)});
    }

    // Global ids of the DOFs that satisfy flag.
    thrust::device_vector<HYPRE_BigInt> flaggedDofs(const uint8_t* flag)
    {
        thrust::device_vector<HYPRE_BigInt> ids(numOwned_);
        auto gids =
            thrust::make_transform_iterator(thrust::counting_iterator<size_t>(0), NsGidOf<RealType>{gid_.data()});
        size_t count = thrust::copy_if(thrust::device, gids, gids + n_, thrust::counting_iterator<size_t>(0),
                                       ids.begin(), NsFlaggedDof<RealType>{space_.isDof(), flag}) -
                       ids.begin();
        ids.resize(count);
        return ids;
    }

    // Local triplets of scale * K (element matrices, then the diagonal slots
    // c M if massCoefficient is given) with the removed slots dropped, reduced to
    // P^T A P plus identity rows at the removed DOFs.
    HypreMatrix assembleReduced(const HypreMatrix& P, RealType scale, const uint8_t* removed,
                                const RealType* massCoefficient, bool withLift)
    {
        const size_t countK = 64 * hex_.count, countM = massCoefficient ? n_ : 0, count = countK + countM;
        thrust::device_vector<HYPRE_BigInt> rows(std::max<size_t>(count, 1)), cols(std::max<size_t>(count, 1));
        thrust::device_vector<RealType> values(std::max<size_t>(count, 1));
        HYPRE_BigInt* r = thrust::raw_pointer_cast(rows.data());
        HYPRE_BigInt* c = thrust::raw_pointer_cast(cols.data());
        RealType* val   = thrust::raw_pointer_cast(values.data());
        if (hex_.count > 0)
            nsLaplacianTripletsKernel<KeyType, RealType><<<elemGrid(), bs()>>>(
                hex_, x(), y(), z(), cview(area_), scale, localStart_, removed, cview(target_), view(lift_),
                withLift ? comps_ : 0, r, c, val);
        if (massCoefficient)
            nsMassTripletsKernel<RealType><<<grid(), bs()>>>(n_, localStart_, massLocal_.data(), removed,
                                                             *massCoefficient, r + countK, c + countK, val + countK);
        cudaCheckError();
        HypreMatrix local = hypreAssemble(MPI_COMM_WORLD, localStart_, localStart_ + (long long)n_, localStart_,
                                          localStart_ + (long long)n_, HypreCoo{r, c, val, HYPRE_Int(count)});
        auto ids          = flaggedDofs(removed);
        return hypreAddIdentityRows(MPI_COMM_WORLD, hypreGalerkin(P, local), thrust::raw_pointer_cast(ids.data()),
                                    HYPRE_Int(ids.size()));
    }

    // a1 = M / dt + nu K and a2 = 3 M / (2 dt) + nu K, K the CVFEM Laplacian,
    // with the fixed velocity DOFs removed. The lift of the removed
    // columns is gathered while assembling a1 (it is the same for a2).
    void buildViscousSolvers()
    {
        HypreMatrix P      = prolongation(1);
        const RealType c1  = RealType(1) / prm_.dt;
        const RealType c2  = RealType(3) / (RealType(2) * prm_.dt);
        viscousBdf1_       = std::make_unique<HypreAmgPcgSolver>("MARS_VAMG");
        viscousBdf1_->setup(MPI_COMM_WORLD, assembleReduced(P, prm_.nu, fixed_.data(), &c1, true), prm_.tolerance,
                            prm_.maxIter);
        for (int d = 0; d < comps_; ++d)
            space_.restrict(lift_[d]);
        if (!prm_.bdf2) return;
        viscousBdf2_ = std::make_unique<HypreAmgPcgSolver>("MARS_VAMG");
        viscousBdf2_->setup(MPI_COMM_WORLD, assembleReduced(P, prm_.nu, fixed_.data(), &c2, false),
                            prm_.tolerance, prm_.maxIter);
    }

    // K over the DOFs with identity rows where p = 0; the operator of the flux correction.
    void buildPressureSolver()
    {
        HypreMatrix P = prolongation(1);
        pressure_     = std::make_unique<HypreAmgPcgSolver>("MARS_PAMG");
        pressure_->setup(MPI_COMM_WORLD, assembleReduced(P, RealType(1), pressureFixed_.data(), nullptr, false),
                         prm_.tolerance, prm_.maxIter);
        checkPressureOperator();
    }

    // K x from Hypre against the matrix-free K x of the time step, for a varied
    // x that is zero where p = 0. Stops if they differ.
    void checkPressureOperator()
    {
        Vector& xs      = work_[1];
        Vector& matFree = work_[2];
        nsTestVectorKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), pressureFixed_.data(), gid_.data(),
                                                       xs.data());
        cudaCheckError();
        space_.prolong(xs);
        applyLaplacian(xs, matFree);

        Vector xDof(numOwned_), assembled(numOwned_);
        nsToDofKernel<RealType><<<grid(), bs()>>>(n_, space_.isDof(), dofIndex_.data(), xs.data(), xDof.data());
        cudaCheckError();
        pressure_->apply(xDof.data(), assembled.data());
        NsOperatorDifference<RealType> diff{space_.isDof(), pressureFixed_.data(), dofIndex_.data(), matFree.data(),
                                            assembled.data()};
        thrust::tuple<double, double> local = thrust::transform_reduce(
            thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n_), diff,
            thrust::make_tuple(0.0, 0.0), NsMaxPair{});
        double loc[2] = {thrust::get<0>(local), thrust::get<1>(local)}, glob[2];
        MPI_Allreduce(loc, glob, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double relative = glob[1] > 0 ? glob[0] / glob[1] : glob[0];
        if (rank_ == 0)
            std::cout << "Pressure operator: assembled vs matrix-free, max |difference| / max |Kx| = "
                      << std::scientific << relative << std::defaultfloat << "\n";
        if (!(relative <= 1e-10))
        {
            if (rank_ == 0) std::cerr << "NavierStokes: the assembled pressure operator is wrong\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    template<class Term>
    double dofSum(Term term) const
    {
        double local = thrust::transform_reduce(thrust::device, thrust::counting_iterator<size_t>(0),
                                                thrust::counting_iterator<size_t>(n_), term, 0.0,
                                                thrust::plus<double>());
        double global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return global;
    }

    void allocate(Vector& a)
    {
        a.resize(n_);
        zero(a);
    }
    void zero(Vector& a) { cudaMemsetAsync(a.data(), 0, a.size() * sizeof(RealType)); }
    void copy(const Vector& from, Vector& to)
    {
        cudaMemcpyAsync(to.data(), from.data(), n_ * sizeof(RealType), cudaMemcpyDeviceToDevice);
    }
    // At least one block: a rank may own no elements, and an empty grid is a launch error.
    int bs() const { return prm_.blockSize; }
    int grid() const { return std::max(1, int((n_ + bs() - 1) / bs())); }
    int elemGrid() const { return std::max(1, int((hex_.count + bs() - 1) / bs())); }

    Domain& domain_;
    Params prm_;
    size_t n_;
    int rank_;
    int comps_;
    Space space_;
    HexElements<KeyType> hex_{};
    Box<RealType> box_{};

    int steps_            = 0;
    bool lastBdf2_        = false;
    int velocityIters_[3] = {0, 0, 0};
    int pressureIters_    = 0;
    Timing timing_;

    // Conditions per node slot
    cstone::DeviceVector<uint8_t> fixed_;         // velocity prescribed
    cstone::DeviceVector<uint8_t> pressureFixed_; // p = 0
    Vector target_[3];                            // prescribed velocity
    Vector openArea_[3];                          // outward area vector of the opening faces

    // DOFs
    int numOwned_         = 0;
    long long dofStart_   = 0;
    long long localStart_ = 0;
    cstone::DeviceVector<int> dofIndex_; // local DOF index of each DOF slot
    Vector gid_;                         // global DOF id of each slot's DOF

    // Geometry and solvers
    Vector area_[3]; // sub-control face area vectors, 12 per local element
    Vector flux_;    // stabilized face fluxes F, 12 per local element
    Vector massLocal_, massDof_;
    Vector lift_[3];
    std::unique_ptr<HypreAmgPcgSolver> viscousBdf1_, viscousBdf2_, pressure_;

    // Work vectors
    Vector um1_[3];   // u^{n-1}
    Vector star_[3];  // u*
    Vector sstar_[3]; // u**
    Vector adv_[3];   // N(u^n)
    Vector advm1_[3]; // N(u^{n-1})
    Vector g_[3];     // M^-1 D^T q
    Vector div_, divStar_, phi_;
    Vector work_[3], faceWork_;
    Vector rhs_, x_; // DOF order
};

} // namespace fem
} // namespace mars
