#pragma once

// Incompressible Navier-Stokes in a planar channel: the solver of the
// Poiseuille tutorial (docs/poiseuille_tutorial.md). CVFEM on Q1 hexes, one GPU
// per MPI rank.
//
// The mesh is one layer of axis-aligned hexes between two z planes, with
//   x = xmin          inlet    velocity (U, 0) prescribed
//   y = ymin, ymax    walls    no slip, velocity 0 (a wall wins where it meets the inlet)
//   x = xmax          outlet   p = 0, velocity free
// The flow is planar: w = 0 and nothing depends on z, so u and v are the unknowns.
//
// Operators, built from the sub-control faces f of the elements (area vector
// A_f from node L to node R) and the lumped mass M:
//   D     divergence    (D u)_L += A_f . (u_L + u_R) / 2 and (D u)_R -= the same,
//                       plus the flux through the inlet and outlet faces
//   D^T   its transpose; M^-1 D^T p approximates -grad p
//   Q     zeroes the velocity components the boundary conditions fix
//   K     CVFEM viscous stiffness, assembled
//   N     skew-symmetric advection
//
// One time step, BDF2 with extrapolated advection (BDF1 on the first step),
// dtEff = 2 dt / 3 (dt with BDF1):
//   predictor   u*  = BDF history + dtEff M^-1 (N_ext + Q D^T p / rho)
//   viscous     (M / dtEff + nu K) u** = M u* / dtEff,    u** = prescribed where fixed
//   projection  A phi = -(rho / dtEff) D u**,              A = D Q M^-1 D^T
//   corrector   u = u** + (dtEff / rho) Q M^-1 D^T phi,    p += phi
// The corrector applies the operator that A inverts, so D u = 0 holds to the
// solver tolerance on any number of ranks.
//
// Both matrices are constant in time. They are assembled once and solved with
// Hypre PCG + BoomerAMG, whose iteration count stays flat when the mesh and
// the number of GPUs grow.
//
// Parallel layout: cornerstone gives each rank its local elements and a halo
// of neighbour elements, and each node one owner rank. Element loops scatter
// into node slots; reverseExchangeNodeHaloAdd adds the ghost slots into their
// owners, and exchangeNodeHalo copies owner values back to the ghosts.

#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include "backend/distributed/unstructured/solvers/mars_hypre_amg_pcg_solver.hpp"

#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/reduce.h>
#include <thrust/remove.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/tuple.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <mpi.h>
#include <type_traits>

namespace mars
{
namespace fem
{

// Elements [first, first + count); corner c of element e is node node[c][e].
template<typename KeyType>
struct HexElementRange
{
    const KeyType* node[8];
    size_t first;
    size_t count;
};

template<typename RealType>
struct ChannelBox
{
    RealType lo[3];
    RealType hi[3];
};

// True if the corners are those of an axis-aligned box whose sides are all longer than tolerance.
template<typename RealType>
__host__ __device__ inline bool isAxisAlignedHex(const RealType corner[8][3], RealType tolerance)
{
    RealType lo[3], hi[3];
    for (int d = 0; d < 3; ++d) lo[d] = hi[d] = corner[0][d];
    for (int c = 1; c < 8; ++c)
        for (int d = 0; d < 3; ++d)
        {
            lo[d] = corner[c][d] < lo[d] ? corner[c][d] : lo[d];
            hi[d] = corner[c][d] > hi[d] ? corner[c][d] : hi[d];
        }
    for (int d = 0; d < 3; ++d)
        if (!(hi[d] - lo[d] > tolerance)) return false;
    unsigned seen = 0;
    for (int c = 0; c < 8; ++c)
    {
        int bits = 0;
        for (int d = 0; d < 3; ++d)
        {
            if (corner[c][d] - lo[d] <= tolerance) continue;
            if (!(hi[d] - corner[c][d] <= tolerance)) return false;
            bits |= 1 << d;
        }
        seen |= 1u << bits;
    }
    return seen == 255;
}

// ---------------------------------------------------------------------------
// Setup kernels
// ---------------------------------------------------------------------------

template<typename KeyType, typename RealType>
__global__ void channelMeshCheckKernel(HexElementRange<KeyType> hex, const RealType* x, const RealType* y,
                                       const RealType* z, RealType tolerance, int* bad)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    RealType corner[8][3];
    for (int c = 0; c < 8; ++c)
    {
        KeyType n  = hex.node[c][e];
        corner[c][0] = x[n];
        corner[c][1] = y[n];
        corner[c][2] = z[n];
    }
    if (!isAxisAlignedHex(corner, tolerance)) atomicAdd(bad, 1);
}

// Boundary conditions of every node slot, from the coordinates, so ghosts agree with their owners.
template<typename RealType>
__global__ void channelMarkBoundaryKernel(size_t n, const RealType* x, const RealType* y, ChannelBox<RealType> box,
                                          RealType eps, RealType inflow, uint8_t* fixed, uint8_t* pressureFixed,
                                          RealType* uTarget)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    bool inlet  = fabs(x[i] - box.lo[0]) < eps;
    bool outlet = fabs(x[i] - box.hi[0]) < eps;
    bool wall   = fabs(y[i] - box.lo[1]) < eps || fabs(y[i] - box.hi[1]) < eps;
    fixed[i]    = inlet || wall;
    uTarget[i]  = inlet && !wall ? inflow : RealType(0);
    // At an inlet-wall corner every edge neighbour has a fixed velocity, so the
    // pressure there acts on nothing: a zero column of A. Fixing it removes that null space.
    pressureFixed[i] = outlet || (inlet && wall);
}

// Each corner takes 1/8 of the element volume, exact for axis-aligned hexes.
template<typename KeyType, typename RealType>
__global__ void channelMassKernel(HexElementRange<KeyType> hex, const RealType* x, const RealType* y,
                                  const RealType* z, RealType* mass)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e     = hex.first + k;
    KeyType n0   = hex.node[0][e];
    RealType lo[3] = {x[n0], y[n0], z[n0]};
    RealType hi[3] = {lo[0], lo[1], lo[2]};
    for (int c = 1; c < 8; ++c)
    {
        KeyType n           = hex.node[c][e];
        const RealType p[3] = {x[n], y[n], z[n]};
        for (int d = 0; d < 3; ++d)
        {
            lo[d] = fmin(lo[d], p[d]);
            hi[d] = fmax(hi[d], p[d]);
        }
    }
    RealType share = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]) * RealType(0.125);
    for (int c = 0; c < 8; ++c)
        atomicAdd(&mass[hex.node[c][e]], share);
}

// Per-node share of the element faces on the plane x = x0, times scale. A face
// counts once, from the owned element on the given side of the plane. A node's
// share is its quarter of the face: the sub-quad of node, edge midpoints and
// face centre, the boundary face of its control volume.
template<typename KeyType, typename RealType>
__global__ void channelPlaneFaceAreaKernel(HexElementRange<KeyType> hex, const RealType* x, const RealType* y,
                                           const RealType* z, double x0, double tolerance, int side, double scale,
                                           RealType* area)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType face[4];
    int onPlane = 0;
    for (int c = 0; c < 8; ++c)
    {
        KeyType n = hex.node[c][e];
        double d  = double(x[n]) - x0;
        if (fabs(d) <= tolerance)
        {
            if (onPlane < 4) face[onPlane] = n;
            ++onPlane;
        }
        else if (d * side <= 0) return;
    }
    if (onPlane != 4) return;
    double cy = 0, cz = 0;
    for (int c = 0; c < 4; ++c)
    {
        cy += 0.25 * y[face[c]];
        cz += 0.25 * z[face[c]];
    }
    // The element's corner order does not give the face cycle; sort by angle.
    double angle[4];
    for (int c = 0; c < 4; ++c)
        angle[c] = atan2(z[face[c]] - cz, y[face[c]] - cy);
    for (int c = 1; c < 4; ++c)
        for (int j = c; j > 0 && angle[j] < angle[j - 1]; --j)
        {
            double a     = angle[j];
            angle[j]     = angle[j - 1];
            angle[j - 1] = a;
            KeyType f    = face[j];
            face[j]      = face[j - 1];
            face[j - 1]  = f;
        }
    for (int c = 0; c < 4; ++c)
    {
        KeyType p = face[c], next = face[(c + 1) % 4], prev = face[(c + 3) % 4];
        double ny = 0.5 * (y[prev] - y[next]), nz = 0.5 * (z[prev] - z[next]);
        double share = 0.5 * fabs((cy - y[p]) * nz - (cz - z[p]) * ny);
        atomicAdd(&area[p], RealType(scale * share));
    }
}

template<typename IndexType>
__global__ void channelDofToNodeKernel(size_t n, const IndexType* nodeToDof, IndexType* dofToNode)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) dofToNode[nodeToDof[i]] = IndexType(i);
}

// Global DOF number of each owned node, as a real number so the node halo can carry it.
template<typename RealType>
__global__ void channelGlobalIdKernel(size_t n, const uint8_t* owned, const int* nodeToDof, long long rowStart,
                                      RealType* id)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) id[i] = owned[i] == 1 ? RealType(rowStart + nodeToDof[i]) : RealType(-1);
}

template<typename RealType>
__global__ void channelToGlobalIdKernel(size_t n, const RealType* id, const int* nodeToDof, HYPRE_BigInt* nodeGid,
                                        HYPRE_BigInt* dofGid)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    nodeGid[i]           = HYPRE_BigInt(id[i]);
    dofGid[nodeToDof[i]] = HYPRE_BigInt(id[i]);
}

template<typename RealType>
__global__ void channelAddMassKernel(int rows, const int* diagPtr, const int* dofToNode, const RealType* mass,
                                     RealType c, RealType* a)
{
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row < rows) a[diagPtr[row]] += mass[dofToNode[row]] * c;
}

// u = U except where the velocity is prescribed.
template<typename RealType>
__global__ void channelInitialConditionKernel(size_t n, const uint8_t* fixed, const RealType* uTarget,
                                              RealType inflow, RealType* u)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n) u[i] = fixed[i] ? uTarget[i] : inflow;
}

// A varied pressure, zero on the pressure-Dirichlet nodes, to test the assembled operator.
template<typename RealType>
__global__ void channelTestVectorKernel(size_t n, const uint8_t* owned, const uint8_t* pressureFixed,
                                        const HYPRE_BigInt* gid, RealType* x)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    x[i] = owned[i] == 1 && !pressureFixed[i] ? RealType(1) + RealType((gid[i] * 7919) % 97) / RealType(97)
                                              : RealType(0);
}

template<typename RealType>
__global__ void channelToDofKernel(size_t n, const uint8_t* owned, const int* nodeToDof, const RealType* q,
                                   RealType* x)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n && owned[i] == 1) x[nodeToDof[i]] = q[i];
}

// Removes the fixed velocity DOFs from the two viscous matrices a1 and a2 (they
// differ on the diagonal only). A fixed row becomes an identity row. In a free
// row the fixed columns move to the right-hand side as the lift -a_rc U_c, which
// keeps the matrices symmetric; v is fixed to 0 and needs no lift.
template<typename RealType>
__global__ void channelEliminateFixedKernel(int rows, const int* rowPtr, const int* colInd, const int* diagPtr,
                                            const int* dofToNode, const uint8_t* fixed, const RealType* uTarget,
                                            RealType* a1, RealType* a2, RealType* liftU)
{
    int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= rows) return;
    if (fixed[dofToNode[row]])
    {
        for (int k = rowPtr[row]; k < rowPtr[row + 1]; ++k)
            a1[k] = a2[k] = RealType(0);
        a1[diagPtr[row]] = a2[diagPtr[row]] = RealType(1);
        liftU[row]                          = RealType(0);
        return;
    }
    RealType lift = 0;
    for (int k = rowPtr[row]; k < rowPtr[row + 1]; ++k)
    {
        int node = dofToNode[colInd[k]];
        if (!fixed[node]) continue;
        lift -= a1[k] * uTarget[node];
        a1[k] = a2[k] = RealType(0);
    }
    liftU[row] = lift;
}

// COO triplets of D and of D S, S = Q M^-1, for the owned rows. The loop covers
// local and halo elements: SFC ownership gives the owner of a node every element
// around it, so its row is complete without communication. Pressure-Dirichlet
// rows are left out (row -1), which also drops their columns from A = (D S) D^T.
// Velocity column 2 g + d is component d at the node with global DOF g.
template<typename KeyType, typename RealType>
__global__ void channelProjectionTripletsKernel(HexElementRange<KeyType> hex, const RealType* ax, const RealType* ay,
                                                const HYPRE_BigInt* gid, const RealType* mass, const uint8_t* fixed,
                                                const uint8_t* pressureFixed, const uint8_t* owned,
                                                HYPRE_BigInt* rows, HYPRE_BigInt* cols, RealType* valD,
                                                RealType* valDS)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n[8];
    for (int c = 0; c < 8; ++c) n[c] = hex.node[c][e];
    size_t out = k * 12 * 8; // per face: 2 rows x 2 columns x 2 components
    for (int ip = 0; ip < 12; ++ip)
    {
        const KeyType ends[2] = {n[hexLRSCV[2 * ip]], n[hexLRSCV[2 * ip + 1]]};
        size_t f              = e * 12 + ip;
        const RealType area[2] = {ax[f], ay[f]};
        for (int r = 0; r < 2; ++r)
        {
            KeyType row        = ends[r];
            HYPRE_BigInt rowId = owned[row] == 1 && !pressureFixed[row] ? gid[row] : HYPRE_BigInt(-1);
            RealType half      = r == 0 ? RealType(0.5) : RealType(-0.5);
            for (int c = 0; c < 2; ++c)
            {
                KeyType col    = ends[c];
                RealType scale = mass[col] == RealType(0) || fixed[col] ? RealType(0) : RealType(1) / mass[col];
                for (int d = 0; d < 2; ++d, ++out)
                {
                    rows[out]  = rowId;
                    cols[out]  = 2 * gid[col] + d;
                    valD[out]  = half * area[d];
                    valDS[out] = half * area[d] * scale;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Time step kernels. Element loops run over the local elements; the per-face
// area vectors are stored for all elements, at e * 12 + face.
// ---------------------------------------------------------------------------

// D u over the sub-control faces: the flux leaves node L and enters node R.
template<typename KeyType, typename RealType>
__global__ void channelDivergenceKernel(HexElementRange<KeyType> hex, const RealType* ax, const RealType* ay,
                                        const RealType* u, const RealType* v, RealType* div)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n[8];
    for (int c = 0; c < 8; ++c) n[c] = hex.node[c][e];
#pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType L = n[hexLRSCV[2 * ip]], R = n[hexLRSCV[2 * ip + 1]];
        size_t f      = e * 12 + ip;
        RealType flow = RealType(0.5) * (u[L] + u[R]) * ax[f] + RealType(0.5) * (v[L] + v[R]) * ay[f];
        atomicAdd(&div[L], flow);
        atomicAdd(&div[R], -flow);
    }
}

// D^T p, the transpose of the scatter above: (p_L - p_R) / 2 A_f to both ends of the face.
template<typename KeyType, typename RealType>
__global__ void channelDivergenceTransposeKernel(HexElementRange<KeyType> hex, const RealType* ax,
                                                 const RealType* ay, const RealType* p, RealType* gx, RealType* gy)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n[8];
    for (int c = 0; c < 8; ++c) n[c] = hex.node[c][e];
#pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType L = n[hexLRSCV[2 * ip]], R = n[hexLRSCV[2 * ip + 1]];
        size_t f    = e * 12 + ip;
        RealType dp = RealType(0.5) * (p[L] - p[R]);
        atomicAdd(&gx[L], dp * ax[f]);
        atomicAdd(&gy[L], dp * ay[f]);
        atomicAdd(&gx[R], dp * ax[f]);
        atomicAdd(&gy[R], dp * ay[f]);
    }
}

// Skew-symmetric advection of u and v by the face flux m = A_f . (u_L + u_R) / 2:
// node L receives -m q_R / 2 and node R receives +m q_L / 2. Its kinetic energy
// production, sum_i q_i (N q)_i, cancels face by face.
template<typename KeyType, typename RealType>
__global__ void channelAdvectionKernel(HexElementRange<KeyType> hex, const RealType* ax, const RealType* ay,
                                       const RealType* u, const RealType* v, RealType* au, RealType* av)
{
    size_t k = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (k >= hex.count) return;
    size_t e = hex.first + k;
    KeyType n[8];
    for (int c = 0; c < 8; ++c) n[c] = hex.node[c][e];
#pragma unroll
    for (int ip = 0; ip < 12; ++ip)
    {
        KeyType L = n[hexLRSCV[2 * ip]], R = n[hexLRSCV[2 * ip + 1]];
        size_t f   = e * 12 + ip;
        RealType m = RealType(0.5) * (u[L] + u[R]) * ax[f] + RealType(0.5) * (v[L] + v[R]) * ay[f];
        atomicAdd(&au[L], -RealType(0.5) * m * u[R]);
        atomicAdd(&au[R], RealType(0.5) * m * u[L]);
        atomicAdd(&av[L], -RealType(0.5) * m * v[R]);
        atomicAdd(&av[R], RealType(0.5) * m * v[L]);
    }
}

// Flux out through the inlet and outlet faces of node i: the inlet carries the
// prescribed velocity, the outlet the computed one. Inlet areas are negative
// (outward normal -x).
template<typename RealType>
__device__ inline RealType channelOpeningFlux(size_t i, const RealType* inletArea, const RealType* outletArea,
                                              const RealType* uTarget, const RealType* u)
{
    return inletArea[i] * uTarget[i] + outletArea[i] * u[i];
}

template<typename RealType>
__global__ void channelAddOpeningFluxKernel(size_t n, const uint8_t* owned, const RealType* inletArea,
                                            const RealType* outletArea, const RealType* uTarget, const RealType* u,
                                            RealType* div)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n && owned[i] == 1) div[i] += channelOpeningFlux(i, inletArea, outletArea, uTarget, u);
}

// The skew form of an opening face with flux m is -m q / 2 at its node.
template<typename RealType>
__global__ void channelOpeningAdvectionKernel(size_t n, const uint8_t* owned, const RealType* inletArea,
                                              const RealType* outletArea, const RealType* uTarget, const RealType* u,
                                              const RealType* v, RealType* au, RealType* av)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n || owned[i] != 1) return;
    RealType m = channelOpeningFlux(i, inletArea, outletArea, uTarget, u);
    au[i] -= RealType(0.5) * u[i] * m;
    av[i] -= RealType(0.5) * v[i] * m;
}

// g = Q M^-1 g at the owned nodes, 0 at the ghosts.
template<typename RealType>
__global__ void channelInverseMassKernel(size_t n, const uint8_t* owned, const uint8_t* fixed, const RealType* mass,
                                         RealType* gx, RealType* gy)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    bool keep = owned[i] == 1 && !fixed[i] && mass[i] != RealType(0);
    RealType s = keep ? RealType(1) / mass[i] : RealType(0);
    gx[i]      = keep ? gx[i] * s : RealType(0);
    gy[i]      = keep ? gy[i] * s : RealType(0);
}

// u* from u^n (and u^{n-1} for BDF2), the advection a = N(u^n) (extrapolated to
// 2 a - a_prev for BDF2), and g = Q M^-1 D^T p^n ~ -grad p^n.
template<typename RealType>
__global__ void channelPredictorKernel(size_t n, const uint8_t* owned, const uint8_t* fixed, int bdf2, RealType dt,
                                       RealType invRho, const RealType* mass, const RealType* uTarget,
                                       const RealType* u, const RealType* v, const RealType* um1, const RealType* vm1,
                                       const RealType* au, const RealType* av, const RealType* aum1,
                                       const RealType* avm1, const RealType* gx, const RealType* gy, RealType* us,
                                       RealType* vs)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n || owned[i] != 1) return;
    if (fixed[i])
    {
        us[i] = uTarget[i];
        vs[i] = RealType(0);
        return;
    }
    RealType V = mass[i];
    if (bdf2)
    {
        RealType c = RealType(2) * dt / RealType(3);
        us[i] = RealType(4) / RealType(3) * u[i] - RealType(1) / RealType(3) * um1[i] +
                c * ((RealType(2) * au[i] - aum1[i]) / V + invRho * gx[i]);
        vs[i] = RealType(4) / RealType(3) * v[i] - RealType(1) / RealType(3) * vm1[i] +
                c * ((RealType(2) * av[i] - avm1[i]) / V + invRho * gy[i]);
    }
    else
    {
        us[i] = u[i] + dt * au[i] / V + dt * invRho * gx[i];
        vs[i] = v[i] + dt * av[i] / V + dt * invRho * gy[i];
    }
}

// Right-hand side M q* / dtEff + lift and initial guess q* of one velocity
// component, per owned DOF. A fixed DOF has an identity row, so its right-hand
// side is the prescribed value. target == nullptr means 0 and lift == nullptr none.
template<typename RealType>
__global__ void channelViscousRhsKernel(size_t n, const uint8_t* owned, const int* nodeToDof, const uint8_t* fixed,
                                        const RealType* mass, RealType invDt, const RealType* qStar,
                                        const RealType* lift, const RealType* target, RealType* rhs, RealType* x)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n || owned[i] != 1) return;
    int dof = nodeToDof[i];
    if (fixed[i]) rhs[dof] = target ? target[i] : RealType(0);
    else rhs[dof] = mass[i] * invDt * qStar[i] + (lift ? lift[dof] : RealType(0));
    x[dof] = qStar[i];
}

// A phi = -(rho / dtEff) D u** per owned DOF; a pressure-Dirichlet row is an identity row with value 0.
template<typename RealType>
__global__ void channelPressureRhsKernel(size_t n, const uint8_t* owned, const int* nodeToDof,
                                         const uint8_t* pressureFixed, RealType coef, const RealType* div,
                                         RealType* rhs)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n && owned[i] == 1) rhs[nodeToDof[i]] = pressureFixed[i] ? RealType(0) : -coef * div[i];
}

template<typename RealType>
__global__ void channelFromDofKernel(size_t n, const uint8_t* owned, const int* nodeToDof, const RealType* x,
                                     RealType* q)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < n && owned[i] == 1) q[i] = x[nodeToDof[i]];
}

// u = u** + s g where free, the prescribed velocity where fixed, and p += phi.
template<typename RealType>
__global__ void channelCorrectorKernel(size_t n, const uint8_t* owned, const uint8_t* fixed, RealType s,
                                       const RealType* uTarget, const RealType* uss, const RealType* vss,
                                       const RealType* gx, const RealType* gy, const RealType* phi, RealType* u,
                                       RealType* v, RealType* p)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n || owned[i] != 1) return;
    u[i] = fixed[i] ? uTarget[i] : uss[i] + s * gx[i];
    v[i] = fixed[i] ? RealType(0) : vss[i] + s * gy[i];
    p[i] += phi[i];
}

// Functors for thrust. They live at namespace scope: nvcc rejects private
// class members as template arguments of the kernels thrust instantiates.

template<typename RealType>
struct ChannelOffPlanes
{
    RealType lo, hi, eps;
    __device__ bool operator()(RealType z) const { return fabs(z - lo) > eps && fabs(z - hi) > eps; }
};

template<typename RealType>
struct ChannelScaleBy
{
    RealType s;
    __device__ RealType operator()(RealType a) const { return a * s; }
};

struct ChannelNotOwnedRow
{
    template<class Tuple>
    __device__ bool operator()(const Tuple& t) const
    {
        return thrust::get<0>(t) < 0;
    }
};

struct ChannelNonZeroValue
{
    template<class Tuple>
    __device__ bool operator()(const Tuple& t) const
    {
        return thrust::get<2>(t) != 0;
    }
};

struct ChannelPressureFixedDof
{
    const int* dofToNode;
    const uint8_t* pressureFixed;
    __device__ bool operator()(int dof) const { return pressureFixed[dofToNode[dof]] != 0; }
};

// (|assembled - matrix-free|, |matrix-free|) per owned row that is not pressure-Dirichlet.
template<typename RealType>
struct ChannelOperatorDifference
{
    const uint8_t* owned;
    const uint8_t* pressureFixed;
    const int* nodeToDof;
    const RealType* matFree;
    const RealType* assembled;
    __device__ thrust::tuple<double, double> operator()(size_t i) const
    {
        if (owned[i] != 1 || pressureFixed[i]) return thrust::make_tuple(0.0, 0.0);
        double ref = matFree[i];
        return thrust::make_tuple(fabs(assembled[nodeToDof[i]] - ref), fabs(ref));
    }
};

struct ChannelMaxPair
{
    __device__ thrust::tuple<double, double> operator()(const thrust::tuple<double, double>& a,
                                                        const thrust::tuple<double, double>& b) const
    {
        return thrust::make_tuple(fmax(thrust::get<0>(a), thrust::get<0>(b)),
                                  fmax(thrust::get<1>(a), thrust::get<1>(b)));
    }
};

// ---------------------------------------------------------------------------
// Reductions over owned node slots
// ---------------------------------------------------------------------------

template<typename RealType>
struct ChannelMassNormTerm
{
    const uint8_t* owned;
    const RealType* mass;
    const RealType* q;
    __device__ double operator()(size_t i) const { return owned[i] == 1 ? double(q[i] * q[i] * mass[i]) : 0.0; }
};

// |D u| / M where the projection enforces continuity: owned, not pressure-Dirichlet.
template<typename RealType>
struct ChannelContinuityTerm
{
    const uint8_t* owned;
    const uint8_t* pressureFixed;
    const RealType* mass;
    const RealType* div;
    __device__ double operator()(size_t i) const
    {
        return owned[i] == 1 && !pressureFixed[i] ? fabs(double(div[i]) / mass[i]) : 0.0;
    }
};

// Sums of the projection check, see ChannelFlow::projectionReport.
struct ChannelProjectionSums
{
    double sum[7]{};
    double max[3]{};
    __host__ __device__ ChannelProjectionSums operator+(const ChannelProjectionSums& o) const
    {
        ChannelProjectionSums r;
        for (int k = 0; k < 7; ++k) r.sum[k] = sum[k] + o.sum[k];
        for (int k = 0; k < 3; ++k) r.max[k] = max[k] > o.max[k] ? max[k] : o.max[k];
        return r;
    }
};

template<typename RealType>
struct ChannelProjectionTerm
{
    const uint8_t* owned;
    const uint8_t* pressureFixed;
    const RealType* mass;
    const RealType* before; // D u** with openings
    const RealType* after;  // D u with openings
    const RealType* action; // A phi
    const RealType* inletArea;
    const RealType* outletArea;
    const RealType* uTarget;
    const RealType* u;
    RealType h; // dtEff / rho

    __device__ ChannelProjectionSums operator()(size_t i) const
    {
        ChannelProjectionSums s;
        if (owned[i] != 1) return s;
        double V = mass[i], a = after[i], b = before[i], ap = action[i];
        if (!(V > 0) || !isfinite(V) || !isfinite(a) || !isfinite(b) || !isfinite(ap) || !isfinite(double(u[i])))
        {
            s.max[2] = 1;
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
        else if (outletArea[i] == RealType(0)) s.max[1] = fabs(a) / V; // inlet-wall corner
        s.sum[4] = double(inletArea[i]) * uTarget[i];
        s.sum[5] = double(outletArea[i]) * u[i];
        s.sum[6] = a;
        return s;
    }
};

// Face areas of the plane x = x0 per node, see channelPlaneFaceAreaKernel.
template<typename KeyType, typename RealType, typename Domain>
void planeFaceAreas(const Domain& domain, double x0, double tolerance, int side, double scale,
                    cstone::DeviceVector<RealType>& area)
{
    const size_t n = domain.getNodeCount();
    area.resize(n);
    cudaMemset(area.data(), 0, n * sizeof(RealType));
    const auto& conn = domain.getElementToNodeConnectivity();
    HexElementRange<KeyType> local{{std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                                    std::get<3>(conn).data(), std::get<4>(conn).data(), std::get<5>(conn).data(),
                                    std::get<6>(conn).data(), std::get<7>(conn).data()},
                                   domain.startIndex(),
                                   domain.localElementCount()};
    if (local.count > 0)
        channelPlaneFaceAreaKernel<KeyType, RealType><<<int((local.count + 255) / 256), 256>>>(
            local, domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(), x0, tolerance, side,
            scale, area.data());
    cudaCheckError();
    domain.reverseExchangeNodeHaloAdd(area);
    domain.exchangeNodeHalo(area);
}

template<typename KeyType, typename RealType>
class ChannelFlow
{
    static_assert(std::is_same_v<RealType, HYPRE_Complex>, "Hypre is built for another floating point type");

public:
    using Domain = ElementDomain<HexTag, RealType, KeyType, cstone::execution::Gpu>;
    using Vector = cstone::DeviceVector<RealType>;

    struct Params
    {
        RealType rho       = 1;
        RealType nu        = RealType(0.01);
        RealType dt        = RealType(0.01);
        RealType inflow    = 1; // inlet velocity U
        bool bdf2          = true;
        RealType tolerance = RealType(1e-10); // relative PCG tolerance of both solves
        int maxIter        = 1000;
        int blockSize      = 256;
    };

    // Milliseconds per stage, summed over the steps, and the summed AMG-PCG iterations.
    struct Timing
    {
        double predictor = 0, viscous = 0, pressure = 0, corrector = 0;
        long velocityIterations = 0, pressureIterations = 0;
    };

    // The last projection, checked. continuity = D u + openings per unit volume,
    // over the rows the projection enforces.
    struct ProjectionReport
    {
        double identity;        // |D u - D u** - h A phi| / |D u**|, should be at the solver tolerance
        double identityRms;     // the same, per unit volume
        double continuityRms;   // RMS of the continuity residual
        double continuityMax;   // its maximum
        double cornerMax;       // continuity at the inlet-wall corners
        double inflow, outflow; // flux through inlet and outlet
        double balance;         // |inflow + outflow| / |inflow|
        double balanceIdentity; // total continuity against inflow + outflow
        bool ok;
    };

    // State per node slot (owned and ghost). w = 0.
    Vector u, v, p;

    ChannelFlow(Domain& domain, const Params& params)
        : domain_(domain)
        , prm_(params)
        , n_(domain.getNodeCount())
        , rank_(domain.rank())
    {
        domain_.cacheNodeCoordinates();
        const auto& conn = domain.getElementToNodeConnectivity();
        const KeyType* nodes[8] = {std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                                   std::get<3>(conn).data(), std::get<4>(conn).data(), std::get<5>(conn).data(),
                                   std::get<6>(conn).data(), std::get<7>(conn).data()};
        for (int c = 0; c < 8; ++c)
        {
            local_.node[c] = nodes[c];
            all_.node[c]   = nodes[c];
        }
        local_.first = domain.startIndex();
        local_.count = domain.localElementCount();
        all_.first   = 0;
        all_.count   = domain.getElementCount();

        for (Vector* f : {&u, &v, &p, &um1_, &vm1_, &us_, &vs_, &uss_, &vss_, &au_, &av_, &aum1_, &avm1_, &gx_,
                          &gy_, &div_, &phi_, &work_[0], &work_[1], &work_[2], &mass_, &uTarget_})
        {
            f->resize(n_);
            zero(*f);
        }
        checkMeshAndMarkBoundaries();
        numberDofs();
        buildGeometry();
        buildViscousSolvers();
        buildPressureSolver();
        setInitialCondition();
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

    // ||u||_M of the streamwise component: sqrt(sum_i M_i u_i^2).
    double streamwiseNorm() const
    {
        return std::sqrt(ownedSum(ChannelMassNormTerm<RealType>{owned(), mass_.data(), u.data()}));
    }

    // max |D u + openings| / M over the rows the projection enforces.
    double maxContinuity()
    {
        divergenceOf(u, v, div_, true);
        double local = thrust::transform_reduce(
            thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n_),
            ChannelContinuityTerm<RealType>{owned(), pressureFixed_.data(), mass_.data(), div_.data()}, 0.0,
            thrust::maximum<double>());
        double global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return global;
    }

    // Checks the last step: the corrected velocity must satisfy
    // D u = D u** + h A phi (h = dtEff / rho) with A the assembled operator the
    // pressure solve inverted, and the flux must balance.
    ProjectionReport projectionReport()
    {
        Vector& before = work_[0];
        Vector& after  = work_[1];
        Vector& action = work_[2];
        divergenceOf(uss_, vss_, before, true);
        divergenceOf(u, v, after, true);
        applyProjection(phi_, action);
        const RealType dtEff = lastBdf2_ ? RealType(2) * prm_.dt / RealType(3) : prm_.dt;
        ChannelProjectionTerm<RealType> term{owned(),           pressureFixed_.data(), mass_.data(),
                                             before.data(),     after.data(),          action.data(),
                                             inletArea_.data(), outletArea_.data(),    uTarget_.data(),
                                             u.data(),          dtEff / prm_.rho};
        ChannelProjectionSums local = thrust::transform_reduce(
            thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n_), term,
            ChannelProjectionSums{}, thrust::plus<ChannelProjectionSums>());
        ChannelProjectionSums g;
        MPI_Allreduce(local.sum, g.sum, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(local.max, g.max, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        ProjectionReport r;
        double scale      = std::max(g.sum[0], 1e-60);
        double volume     = std::max(g.sum[3], 1e-60);
        r.identity        = std::sqrt(g.sum[2] / scale);
        r.identityRms     = std::sqrt(g.sum[2] / volume);
        r.continuityRms   = std::sqrt(g.sum[1] / volume);
        r.continuityMax   = g.max[0];
        r.cornerMax       = g.max[1];
        r.inflow          = g.sum[4];
        r.outflow         = g.sum[5];
        r.balance         = std::abs(g.sum[4] + g.sum[5]) / std::max(std::abs(g.sum[4]), 1e-30);
        r.balanceIdentity = std::abs(g.sum[6] - g.sum[4] - g.sum[5]) /
                            std::max(std::abs(g.sum[4]) + std::abs(g.sum[5]), 1e-30);
        // When both divergences approach zero the relative identity loses its meaning; then the
        // absolute one must be at roundoff of the velocity gradient scale U / H.
        double identityFloor = 1e-10 * double(prm_.inflow) / double(box_.hi[1] - box_.lo[1]);
        r.ok = g.max[2] == 0 && std::isfinite(r.identity) && std::isfinite(r.continuityRms) &&
               std::isfinite(r.balanceIdentity) && !(r.identity > 1e-7 && r.identityRms > identityFloor) &&
               r.balanceIdentity <= 1e-10 && r.cornerMax <= 1e-8;
        return r;
    }

    const Vector& mass() const { return mass_; }
    const ChannelBox<RealType>& box() const { return box_; }
    const Timing& timing() const { return timing_; }
    void resetTiming() { timing_ = {}; }
    int velocityIterations(int component) const { return velocityIters_[component]; }
    int pressureIterations() const { return pressureIters_; }
    long long globalDofs() const { return numGlobal_; }

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

    // 1. u* from the history, the advection and the old pressure.
    void predict(bool bdf2)
    {
        gradientOf(p);
        zero(au_);
        zero(av_);
        if (local_.count > 0)
            channelAdvectionKernel<KeyType, RealType><<<elemGrid(), bs()>>>(local_, areaX_.data(), areaY_.data(),
                                                                            u.data(), v.data(), au_.data(), av_.data());
        cudaCheckError();
        domain_.reverseExchangeNodeHaloAdd(au_);
        domain_.reverseExchangeNodeHaloAdd(av_);
        channelOpeningAdvectionKernel<RealType><<<grid(), bs()>>>(n_, owned(), inletArea_.data(), outletArea_.data(),
                                                                  uTarget_.data(), u.data(), v.data(), au_.data(),
                                                                  av_.data());
        channelPredictorKernel<RealType><<<grid(), bs()>>>(
            n_, owned(), fixed_.data(), int(bdf2), prm_.dt, RealType(1) / prm_.rho, mass_.data(), uTarget_.data(),
            u.data(), v.data(), um1_.data(), vm1_.data(), au_.data(), av_.data(), aum1_.data(), avm1_.data(),
            gx_.data(), gy_.data(), us_.data(), vs_.data());
        cudaCheckError();
        au_.swap(aum1_);
        av_.swap(avm1_);
    }

    // 2. (M / dtEff + nu K) u** = M u* / dtEff + lift for u and v, warm-started from u*.
    bool solveViscous(bool bdf2)
    {
        const RealType invDt      = bdf2 ? RealType(3) / (RealType(2) * prm_.dt) : RealType(1) / prm_.dt;
        HypreAmgPcgSolver& solver = bdf2 ? *viscousBdf2_ : *viscousBdf1_;
        Vector* star[2]           = {&us_, &vs_};
        Vector* result[2]         = {&uss_, &vss_};
        const RealType* lift[2]   = {liftU_.data(), nullptr};
        const RealType* target[2] = {uTarget_.data(), nullptr};
        for (int c = 0; c < 2; ++c)
        {
            channelViscousRhsKernel<RealType><<<grid(), bs()>>>(n_, owned(), nodeToDof_.data(), fixed_.data(),
                                                                mass_.data(), invDt, star[c]->data(), lift[c],
                                                                target[c], rhs_.data(), x_.data());
            cudaCheckError();
            velocityIters_[c] = solver.solve(rhs_.data(), x_.data(), true);
            if (velocityIters_[c] < 0) return false;
            timing_.velocityIterations += velocityIters_[c];
            channelFromDofKernel<RealType><<<grid(), bs()>>>(n_, owned(), nodeToDof_.data(), x_.data(),
                                                             result[c]->data());
            cudaCheckError();
            domain_.exchangeNodeHalo(*result[c]);
        }
        return true;
    }

    // 3. A phi = -(rho / dtEff) (D u** + openings).
    bool project(bool bdf2)
    {
        const RealType invDt = bdf2 ? RealType(3) / (RealType(2) * prm_.dt) : RealType(1) / prm_.dt;
        divergenceOf(uss_, vss_, div_, true);
        channelPressureRhsKernel<RealType><<<grid(), bs()>>>(n_, owned(), nodeToDof_.data(), pressureFixed_.data(),
                                                             prm_.rho * invDt, div_.data(), rhs_.data());
        cudaCheckError();
        pressureIters_ = pressure_->solve(rhs_.data(), x_.data());
        if (pressureIters_ < 0) return false;
        timing_.pressureIterations += pressureIters_;
        channelFromDofKernel<RealType><<<grid(), bs()>>>(n_, owned(), nodeToDof_.data(), x_.data(), phi_.data());
        cudaCheckError();
        domain_.exchangeNodeHalo(phi_);
        return true;
    }

    // 4. u = u** + (dtEff / rho) Q M^-1 D^T phi, p += phi. u^n becomes the BDF2 history.
    void correct(bool bdf2)
    {
        const RealType dtEff = bdf2 ? RealType(2) * prm_.dt / RealType(3) : prm_.dt;
        u.swap(um1_);
        v.swap(vm1_);
        gradientOf(phi_);
        channelCorrectorKernel<RealType><<<grid(), bs()>>>(n_, owned(), fixed_.data(), dtEff / prm_.rho,
                                                           uTarget_.data(), uss_.data(), vss_.data(), gx_.data(),
                                                           gy_.data(), phi_.data(), u.data(), v.data(), p.data());
        cudaCheckError();
        domain_.exchangeNodeHalo(u);
        domain_.exchangeNodeHalo(v);
        domain_.exchangeNodeHalo(p);
    }

    // gx, gy = Q M^-1 D^T q (~ -grad q) at the owned nodes. q must be current on the ghosts.
    void gradientOf(const Vector& q)
    {
        zero(gx_);
        zero(gy_);
        if (local_.count > 0)
            channelDivergenceTransposeKernel<KeyType, RealType><<<elemGrid(), bs()>>>(
                local_, areaX_.data(), areaY_.data(), q.data(), gx_.data(), gy_.data());
        cudaCheckError();
        domain_.reverseExchangeNodeHaloAdd(gx_);
        domain_.reverseExchangeNodeHaloAdd(gy_);
        channelInverseMassKernel<RealType><<<grid(), bs()>>>(n_, owned(), fixed_.data(), mass_.data(), gx_.data(),
                                                             gy_.data());
        cudaCheckError();
    }

    // out = D (a, b) at the owned nodes, plus the opening fluxes if requested.
    void divergenceOf(const Vector& a, const Vector& b, Vector& out, bool openings)
    {
        zero(out);
        if (local_.count > 0)
            channelDivergenceKernel<KeyType, RealType><<<elemGrid(), bs()>>>(local_, areaX_.data(), areaY_.data(),
                                                                             a.data(), b.data(), out.data());
        cudaCheckError();
        domain_.reverseExchangeNodeHaloAdd(out);
        if (openings)
            channelAddOpeningFluxKernel<RealType><<<grid(), bs()>>>(n_, owned(), inletArea_.data(),
                                                                    outletArea_.data(), uTarget_.data(), a.data(),
                                                                    out.data());
        cudaCheckError();
    }

    // out = D Q M^-1 D^T x at the owned nodes, matrix-free, with the kernels of the time step.
    void applyProjection(const Vector& x, Vector& out)
    {
        gradientOf(x);
        domain_.exchangeNodeHalo(gx_);
        domain_.exchangeNodeHalo(gy_);
        divergenceOf(gx_, gy_, out, false);
    }

    void checkMeshAndMarkBoundaries()
    {
        const RealType* coord[3] = {domain_.getNodeX().data(), domain_.getNodeY().data(), domain_.getNodeZ().data()};
        const RealType inf       = std::numeric_limits<RealType>::infinity();
        for (int d = 0; d < 3; ++d)
        {
            RealType lo = thrust::reduce(thrust::device, coord[d], coord[d] + n_, inf, thrust::minimum<RealType>());
            RealType hi = thrust::reduce(thrust::device, coord[d], coord[d] + n_, -inf, thrust::maximum<RealType>());
            MPI_Allreduce(&lo, &box_.lo[d], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&hi, &box_.hi[d], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }
        // 1e-5 of the height: the first interior node of the graded test mesh is 1e-4 off the wall.
        const RealType eps = RealType(1e-5) * std::max(RealType(1), box_.hi[1] - box_.lo[1]);

        // One layer in z, every cell an axis-aligned box.
        RealType zlo = box_.lo[2], zhi = box_.hi[2];
        long long bad =
            thrust::count_if(thrust::device, coord[2], coord[2] + n_, ChannelOffPlanes<RealType>{zlo, zhi, eps});
        thrust::device_vector<int> badCells(1, 0);
        RealType extent = std::max({box_.hi[0] - box_.lo[0], box_.hi[1] - box_.lo[1], box_.hi[2] - box_.lo[2]});
        if (local_.count > 0)
            channelMeshCheckKernel<KeyType, RealType><<<elemGrid(), bs()>>>(
                local_, coord[0], coord[1], coord[2], RealType(1e-10) * extent,
                thrust::raw_pointer_cast(badCells.data()));
        cudaCheckError();
        bad += int(badCells[0]);
        long long badGlobal = 0;
        MPI_Allreduce(&bad, &badGlobal, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (badGlobal > 0 || !(zhi > zlo))
        {
            if (rank_ == 0)
                std::cerr << "ChannelFlow: the mesh must be one layer of axis-aligned hexes between two z planes\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        fixed_.resize(n_);
        pressureFixed_.resize(n_);
        channelMarkBoundaryKernel<RealType><<<grid(), bs()>>>(n_, coord[0], coord[1], box_, eps, prm_.inflow,
                                                              fixed_.data(), pressureFixed_.data(), uTarget_.data());
        cudaCheckError();
    }

    // Owned nodes are DOFs [0, numOwned_), ghosts follow. Hypre numbers the
    // owned DOFs of rank r contiguously from rowStart_.
    void numberDofs()
    {
        nodeToDof_.resize(n_);
        dofToNode_.resize(n_);
        numOwned_ = buildDofMappingGpu<KeyType>(owned(), nodeToDof_.data(), n_);
        channelDofToNodeKernel<int><<<grid(), bs()>>>(n_, nodeToDof_.data(), dofToNode_.data());
        cudaCheckError();

        long long ownedCount = numOwned_;
        MPI_Exscan(&ownedCount, &rowStart_, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (rank_ == 0) rowStart_ = 0;
        MPI_Allreduce(&ownedCount, &numGlobal_, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        Vector& id = work_[0];
        channelGlobalIdKernel<RealType><<<grid(), bs()>>>(n_, owned(), nodeToDof_.data(), rowStart_, id.data());
        cudaCheckError();
        domain_.exchangeNodeHalo(id);
        nodeGid_.resize(n_);
        dofGid_.resize(n_);
        channelToGlobalIdKernel<RealType><<<grid(), bs()>>>(n_, id.data(), nodeToDof_.data(),
                                                            thrust::raw_pointer_cast(nodeGid_.data()),
                                                            thrust::raw_pointer_cast(dofGid_.data()));
        cudaCheckError();
        long long missing = thrust::count(thrust::device, nodeGid_.begin(), nodeGid_.end(), HYPRE_BigInt(-1));
        MPI_Allreduce(MPI_IN_PLACE, &missing, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (missing > 0)
        {
            if (rank_ == 0) std::cerr << "ChannelFlow: " << missing << " ghost nodes have no owner\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        rhs_.resize(numOwned_);
        x_.resize(numOwned_);
    }

    // Sub-control face area vectors of all elements, the lumped mass, and the
    // outward face areas of the inlet (negative) and the outlet.
    void buildGeometry()
    {
        const size_t m = all_.count;
        areaX_.resize(12 * m);
        areaY_.resize(12 * m);
        areaZ_.resize(12 * m);
        const RealType* x = domain_.getNodeX().data();
        const RealType* y = domain_.getNodeY().data();
        const RealType* z = domain_.getNodeZ().data();
        if (m > 0)
            precomputeAreaVectorsGpu<KeyType, RealType>(all_.node[0], all_.node[1], all_.node[2], all_.node[3],
                                                        all_.node[4], all_.node[5], all_.node[6], all_.node[7], m, x,
                                                        y, z, areaX_.data(), areaY_.data(), areaZ_.data());
        cudaCheckError();
        if (local_.count > 0)
            channelMassKernel<KeyType, RealType><<<elemGrid(), bs()>>>(local_, x, y, z, mass_.data());
        cudaCheckError();
        domain_.reverseExchangeNodeHaloAdd(mass_);
        domain_.exchangeNodeHalo(mass_);

        const double eps = 1e-5 * std::max(1.0, double(box_.hi[1] - box_.lo[1]));
        planeFaceAreas<KeyType>(domain_, box_.lo[0], eps, +1, -1.0, inletArea_);
        planeFaceAreas<KeyType>(domain_, box_.hi[0], eps, -1, +1.0, outletArea_);
    }

    // CVFEM stiffness K (unit diffusivity, no advection) on the full 27-point
    // sparsity, then a1 = M / dt + nu K and a2 = 3 M / (2 dt) + nu K with the
    // fixed velocity DOFs removed. Each gets its own BoomerAMG hierarchy.
    void buildViscousSolvers()
    {
        const auto& n = all_.node;
        rowPtr_.resize(n_ + 1);
        diagPtr_.resize(n_);
        int nnz = CvfemSparsityBuilder<KeyType>::buildFullSparsity(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7],
                                                                     all_.count, nodeToDof_.data(), int(n_),
                                                                     rowPtr_.data(), nullptr, nullptr, 0);
        colInd_.resize(nnz);
        CvfemSparsityBuilder<KeyType>::buildFullSparsity(n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], all_.count,
                                                         nodeToDof_.data(), int(n_), rowPtr_.data(), colInd_.data(),
                                                         diagPtr_.data(), 0);
        Vector a1(nnz, RealType(0)), a2(nnz);
        {
            // gamma = 1 and zero advection turn the CVFEM assembler into the Laplacian.
            Vector ones(n_, RealType(1)), zeros(n_, RealType(0)), zeroFlux(12 * all_.count, RealType(0)),
                unusedRhs(n_, RealType(0));
            CSRMatrix<RealType> host{rowPtr_.data(), colInd_.data(), a1.data(), diagPtr_.data(), int(n_), nnz,
                                     numOwned_};
            thrust::device_vector<CSRMatrix<RealType>> matrix(1, host);
            typename CvfemHexAssembler<KeyType, RealType>::Config config;
            config.blockSize = bs();
            config.variant   = CvfemKernelVariant::Tensor;
            CvfemHexAssembler<KeyType, RealType>::assembleFull(
                n[0], n[1], n[2], n[3], n[4], n[5], n[6], n[7], all_.count, domain_.getNodeX().data(),
                domain_.getNodeY().data(), domain_.getNodeZ().data(), ones.data(), zeros.data(), zeros.data(),
                zeros.data(), zeros.data(), zeros.data(), zeroFlux.data(), areaX_.data(), areaY_.data(),
                areaZ_.data(), nodeToDof_.data(), owned(), thrust::raw_pointer_cast(matrix.data()), unusedRhs.data(),
                config);
            cudaCheckError();
            cudaDeviceSynchronize();
        }
        thrust::transform(thrust::device, a1.data(), a1.data() + nnz, a1.data(), ChannelScaleBy<RealType>{prm_.nu});
        const int rowGrid = std::max(1, (numOwned_ + bs() - 1) / bs());
        channelAddMassKernel<RealType><<<rowGrid, bs()>>>(numOwned_, diagPtr_.data(), dofToNode_.data(),
                                                          mass_.data(), RealType(1) / prm_.dt, a1.data());
        cudaMemcpy(a2.data(), a1.data(), nnz * sizeof(RealType), cudaMemcpyDeviceToDevice);
        channelAddMassKernel<RealType><<<rowGrid, bs()>>>(numOwned_, diagPtr_.data(), dofToNode_.data(),
                                                          mass_.data(), RealType(1) / (RealType(2) * prm_.dt),
                                                          a2.data());
        liftU_.resize(numOwned_);
        channelEliminateFixedKernel<RealType><<<rowGrid, bs()>>>(numOwned_, rowPtr_.data(), colInd_.data(),
                                                                 diagPtr_.data(), dofToNode_.data(), fixed_.data(),
                                                                 uTarget_.data(), a1.data(), a2.data(),
                                                                 liftU_.data());
        cudaCheckError();
        cudaDeviceSynchronize();

        auto build = [&](const Vector& values) {
            auto solver = std::make_unique<HypreAmgPcgSolver>("MARS_VAMG");
            solver->setupCsr(MPI_COMM_WORLD, HYPRE_BigInt(rowStart_), HYPRE_BigInt(rowStart_ + numOwned_),
                             rowPtr_.data(), colInd_.data(), values.data(), thrust::raw_pointer_cast(dofGid_.data()),
                             n_, prm_.tolerance, prm_.maxIter);
            return solver;
        };
        viscousBdf1_ = build(a1);
        if (prm_.bdf2) viscousBdf2_ = build(a2);
    }

    // Assembles A = (D S) D^T with Hypre and checks it against the matrix-free
    // operator the time step applies.
    void buildPressureSolver()
    {
        const size_t count = all_.count * 12 * 8;
        if (count > size_t(std::numeric_limits<HYPRE_Int>::max()))
        {
            std::cerr << "ChannelFlow: too many projection triplets on rank " << rank_ << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        thrust::device_vector<HYPRE_BigInt> rows(count), cols(count);
        thrust::device_vector<RealType> valD(count), valDS(count);
        if (all_.count > 0)
            channelProjectionTripletsKernel<KeyType, RealType><<<int((all_.count + bs() - 1) / bs()), bs()>>>(
                all_, areaX_.data(), areaY_.data(), thrust::raw_pointer_cast(nodeGid_.data()), mass_.data(),
                fixed_.data(), pressureFixed_.data(), owned(), thrust::raw_pointer_cast(rows.data()),
                thrust::raw_pointer_cast(cols.data()), thrust::raw_pointer_cast(valD.data()),
                thrust::raw_pointer_cast(valDS.data()));
        cudaCheckError();

        // Drop the rows other ranks own, then the zero entries of D S.
        auto all         = thrust::make_zip_iterator(thrust::make_tuple(rows.begin(), cols.begin(), valD.begin(),
                                                                         valDS.begin()));
        const size_t nnzD = thrust::remove_if(thrust::device, all, all + count, ChannelNotOwnedRow{}) - all;
        thrust::device_vector<HYPRE_BigInt> rowsS(nnzD), colsS(nnzD);
        thrust::device_vector<RealType> valS(nnzD);
        auto in           = thrust::make_zip_iterator(thrust::make_tuple(rows.begin(), cols.begin(), valDS.begin()));
        auto out          = thrust::make_zip_iterator(thrust::make_tuple(rowsS.begin(), colsS.begin(), valS.begin()));
        const size_t nnzS = thrust::copy_if(thrust::device, in, in + nnzD, out, ChannelNonZeroValue{}) - out;

        thrust::device_vector<HYPRE_BigInt> pinned(numOwned_);
        const size_t numPinned =
            thrust::copy_if(thrust::device, dofGid_.begin(), dofGid_.begin() + numOwned_,
                            thrust::counting_iterator<int>(0), pinned.begin(),
                            ChannelPressureFixedDof{dofToNode_.data(), pressureFixed_.data()}) -
            pinned.begin();

        using Coo = HypreAmgPcgSolver::Coo;
        pressure_ = std::make_unique<HypreAmgPcgSolver>("MARS_PAMG");
        pressure_->setupProjection(
            MPI_COMM_WORLD, HYPRE_BigInt(rowStart_), HYPRE_BigInt(rowStart_ + numOwned_), 2,
            Coo{thrust::raw_pointer_cast(rows.data()), thrust::raw_pointer_cast(cols.data()),
                thrust::raw_pointer_cast(valD.data()), HYPRE_Int(nnzD)},
            Coo{thrust::raw_pointer_cast(rowsS.data()), thrust::raw_pointer_cast(colsS.data()),
                thrust::raw_pointer_cast(valS.data()), HYPRE_Int(nnzS)},
            thrust::raw_pointer_cast(pinned.data()), HYPRE_Int(numPinned), prm_.tolerance, prm_.maxIter);
        checkPressureOperator();
    }

    // A x from Hypre against the matrix-free D Q M^-1 D^T x, for a varied x
    // that is zero on the pressure-Dirichlet nodes. Stops if they differ.
    void checkPressureOperator()
    {
        Vector& x       = work_[0];
        Vector& matFree = work_[1];
        channelTestVectorKernel<RealType><<<grid(), bs()>>>(n_, owned(), pressureFixed_.data(),
                                                            thrust::raw_pointer_cast(nodeGid_.data()), x.data());
        cudaCheckError();
        domain_.exchangeNodeHalo(x);
        applyProjection(x, matFree);

        thrust::device_vector<RealType> xDof(numOwned_), assembled(numOwned_);
        channelToDofKernel<RealType><<<grid(), bs()>>>(n_, owned(), nodeToDof_.data(), x.data(),
                                                       thrust::raw_pointer_cast(xDof.data()));
        cudaCheckError();
        pressure_->apply(thrust::raw_pointer_cast(xDof.data()), thrust::raw_pointer_cast(assembled.data()));
        ChannelOperatorDifference<RealType> diff{owned(), pressureFixed_.data(), nodeToDof_.data(), matFree.data(),
                                                 thrust::raw_pointer_cast(assembled.data())};
        thrust::tuple<double, double> local = thrust::transform_reduce(
            thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n_), diff,
            thrust::make_tuple(0.0, 0.0), ChannelMaxPair{});
        double loc[2] = {thrust::get<0>(local), thrust::get<1>(local)}, glob[2];
        MPI_Allreduce(loc, glob, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double relative = glob[1] > 0 ? glob[0] / glob[1] : glob[0];
        if (rank_ == 0)
            std::cout << "Pressure operator: assembled vs matrix-free, max |difference| / max |Ax| = "
                      << std::scientific << relative << std::defaultfloat << "\n";
        if (!(relative <= 1e-10))
        {
            if (rank_ == 0) std::cerr << "ChannelFlow: the assembled pressure operator is wrong\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Uniform flow u = U everywhere except on the walls, p = 0.
    void setInitialCondition()
    {
        channelInitialConditionKernel<RealType><<<grid(), bs()>>>(n_, fixed_.data(), uTarget_.data(), prm_.inflow,
                                                                  u.data());
        cudaCheckError();
        zero(v);
        zero(p);
        steps_ = 0;
    }

    template<class Term>
    double ownedSum(Term term) const
    {
        double local = thrust::transform_reduce(thrust::device, thrust::counting_iterator<size_t>(0),
                                                thrust::counting_iterator<size_t>(n_), term, 0.0,
                                                thrust::plus<double>());
        double global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return global;
    }

    const uint8_t* owned() const { return domain_.getNodeOwnershipMap().data(); }
    void zero(Vector& a) { cudaMemsetAsync(a.data(), 0, a.size() * sizeof(RealType)); }
    // At least one block: a rank may own no elements, and an empty grid is a launch error.
    int bs() const { return prm_.blockSize; }
    int grid() const { return std::max(1, int((n_ + bs() - 1) / bs())); }
    int elemGrid() const { return std::max(1, int((local_.count + bs() - 1) / bs())); }

    Domain& domain_;
    Params prm_;
    size_t n_;
    int rank_;
    HexElementRange<KeyType> local_{}; // this rank's elements
    HexElementRange<KeyType> all_{};   // local and halo elements
    ChannelBox<RealType> box_{};

    int steps_     = 0;
    bool lastBdf2_ = false;
    int velocityIters_[2] = {0, 0};
    int pressureIters_    = 0;
    Timing timing_;

    // Boundary conditions per node slot
    cstone::DeviceVector<uint8_t> fixed_;         // velocity prescribed
    cstone::DeviceVector<uint8_t> pressureFixed_; // p = 0
    Vector uTarget_;                              // prescribed u; v is 0
    Vector inletArea_, outletArea_;               // outward x-area of the opening faces per node

    // DOFs
    int numOwned_        = 0;
    long long rowStart_  = 0;
    long long numGlobal_ = 0;
    cstone::DeviceVector<int> nodeToDof_, dofToNode_;
    thrust::device_vector<HYPRE_BigInt> nodeGid_, dofGid_; // global DOF per node slot and per local DOF

    // Geometry and matrices
    Vector areaX_, areaY_, areaZ_; // sub-control face area vectors, 12 per element
    Vector mass_;
    cstone::DeviceVector<int> rowPtr_, colInd_, diagPtr_;
    Vector liftU_;
    std::unique_ptr<HypreAmgPcgSolver> viscousBdf1_, viscousBdf2_, pressure_;

    // Work vectors
    Vector um1_, vm1_;   // u^{n-1}
    Vector us_, vs_;     // u*
    Vector uss_, vss_;   // u**
    Vector au_, av_;     // N(u^n)
    Vector aum1_, avm1_; // N(u^{n-1})
    Vector gx_, gy_;     // Q M^-1 D^T q
    Vector div_, phi_;
    Vector work_[3];
    Vector rhs_, x_; // owned DOFs
};

} // namespace fem
} // namespace mars
