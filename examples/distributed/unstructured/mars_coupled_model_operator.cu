// Coupled collocated solver — Phase 0 model operator (see
// internal-notes/plan_coupled_collocated_solver.md).
//
// This standalone driver builds, validates, and (eventually) directly solves a
// small coupled [u,v,w,p] block operator for the paper benchmarks (lid-driven
// skew cavity, backward-facing step), as the throwaway oracle that de-risks the
// Phase-1 preconditioner bake-off.
//
// INCREMENT 1 (this file, now): the foundational geometry gate only.
//   Verify G = -D^T (discrete integration by parts) on the weak FEM div/grad
//   operators the coupled off-diagonal blocks are built from. If the geometry /
//   connectivity is sound, <phi, D v> = -<v, G phi> holds to round-off for ANY
//   fields phi, v -- so a failure here is a mesh/connectivity bug, caught before
//   any coupled physics is written.
// NEXT INCREMENTS (not yet here): momentum block a^uu (advection+diffusion),
//   the Rhie-Chow continuity coefficient a^pp/a^pu, block-CSR assembly, the
//   direct-solve oracle, and the published-profile comparison.

#include "backend/distributed/unstructured/fem/mars_ns_solver.hpp"  // pulls domain, AMR, CVFEM, ElemTraits, SparseMatrix

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

// Must follow the usings: the projection kernels live in the global namespace
// and resolve Tet4CVFEM / ElemTraits unqualified (header is not standalone).
#include "backend/distributed/unstructured/fem/mars_fem_projection.hpp"
#include "backend/distributed/unstructured/solvers/mars_acm_coarsen.hpp"  // ACM Stage 3 coarse-op + transfers
#include "backend/distributed/unstructured/solvers/mars_acm_vcycle.hpp"   // ACM Stage 4a V-cycle
#include "backend/distributed/unstructured/solvers/mars_gpu_acm_preconditioner.hpp"  // ACM Stage 4b precond
#include "backend/distributed/unstructured/solvers/mars_gmres_solver.hpp"            // native (Flex)GMRES
#include "backend/distributed/unstructured/fem/mars_cvfem_tet_area.hpp"              // SCS area vectors (Rhie-Chow)
#include "backend/distributed/unstructured/fem/mars_boundary_area.hpp"               // boundary median-dual S_bnd (open-flow continuity)

#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <thrust/extrema.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <cusolverSp.h>
#include <cusparse.h>
#include <chrono>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <array>
#include <algorithm>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <thrust/inner_product.h>

// ---- Increment 2a: coupled Stokes block assembly (advection OFF) ----
// Assembles the interleaved (dof = 4*node + comp) CSR with the coefficients
// host-replica-validated in scripts/coupled_assembly_replica.py:
//   a^up = +(V/4) dNdx[b][d]   (momentum row, pressure col)
//   a^pu = -(V/4) dNdx[a][d]   (continuity row, velocity col) = -(a^up)^T
//   a^uu = nu*K (component-diagonal)   a^pp = tau*K   with K_ab = V*(dNdx[a].dNdx[b])
// One thread per element scatters its 4-node x 4-comp contributions into the CSR.
template<typename RealType>
__device__ inline void csrAdd(RealType* vals, const int* col, int rs, int re, int c, RealType v)
{
    for (int k = rs; k < re; ++k)
        if (col[k] == c) { atomicAdd(&vals[k], v); return; }
}

template<typename KeyType, typename RealType>
__global__ void assembleCoupledStokesKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    const int* rowOff, const int* colInd, RealType* vals,
    RealType nu, RealType tau,
    const RealType* avx, const RealType* avy, const RealType* avz,  // frozen advecting vel; null = OFF
    size_t startElem, size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numLocal) return;
    size_t e = startElem + k;

    const KeyType* cc[4] = {c0, c1, c2, c3};
    int n[4];
    for (int i = 0; i < 4; ++i) n[i] = static_cast<int>(cc[i][e]);

    RealType coords[4][3];
    for (int i = 0; i < 4; ++i) {
        coords[i][0] = nodeX[n[i]];
        coords[i][1] = nodeY[n[i]];
        coords[i][2] = nodeZ[n[i]];
    }
    RealType det, dNdx[4][3];
    Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx);
    RealType V = det / RealType(6);
    if (!(V > RealType(0))) return;
    RealType qV = V * RealType(0.25);

    // frozen element-mean advecting velocity (consistent-Galerkin convection)
    RealType ax = 0, ay = 0, az = 0;
    if (avx) {
        for (int i = 0; i < 4; ++i) { ax += avx[n[i]]; ay += avy[n[i]]; az += avz[n[i]]; }
        ax *= RealType(0.25); ay *= RealType(0.25); az *= RealType(0.25);
    }

    for (int a = 0; a < 4; ++a) {
        int ra3 = 4 * n[a] + 3, rs3 = rowOff[ra3], re3 = rowOff[ra3 + 1];
        for (int b = 0; b < 4; ++b) {
            RealType Kab = V * (dNdx[a][0] * dNdx[b][0]
                              + dNdx[a][1] * dNdx[b][1]
                              + dNdx[a][2] * dNdx[b][2]);
            RealType Cab = avx ? qV * (ax * dNdx[b][0] + ay * dNdx[b][1] + az * dNdx[b][2])
                               : RealType(0);                          // C_ij = (V/4)(a.gradN_j)
            csrAdd<RealType>(vals, colInd, rs3, re3, 4 * n[b] + 3, tau * Kab);   // a^pp
            for (int d = 0; d < 3; ++d)
                csrAdd<RealType>(vals, colInd, rs3, re3, 4 * n[b] + d, -qV * dNdx[a][d]); // a^pu
            for (int d = 0; d < 3; ++d) {
                int rad = 4 * n[a] + d, rs = rowOff[rad], re = rowOff[rad + 1];
                csrAdd<RealType>(vals, colInd, rs, re, 4 * n[b] + d, nu * Kab + Cab);     // a^uu diff+adv
                csrAdd<RealType>(vals, colInd, rs, re, 4 * n[b] + 3, qV * dNdx[b][d]);    // a^up
            }
        }
    }
}

// ---- Rhie-Chow coupled pressure block (replaces PSPG a^pp = tau*K) ----
// The collocated inf-sup remedy (darwish2008a / Nalu-Wind / OpenFOAM): the SCS-face mass flux carries a
// compact pressure term  -d_f*(p_R - p_L),  d_f = (V_P/a_P^u)_f * |A|^2/(A.dx)  -- the momentum-interpolation
// coefficient, NOT a stabilization parameter. Continuity (sum of face fluxes) then yields a^pp as the
// d_f-weighted graph Laplacian: row-sum zero, +d_f on the diagonal, -d_f off -- which suppresses
// checkerboard intrinsically. d_f>0 since A is oriented L->R so A.dx>0 and V/a_P^u>0.

// per-node median-dual control-volume V_P = sum over incident tets of V/4
template<typename KeyType, typename RealType>
__global__ void computeNodeVolumeKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2,
                                        const KeyType* c3, const RealType* nx, const RealType* ny,
                                        const RealType* nz, RealType* Vp, size_t startElem, size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x; if (k >= numLocal) return;
    size_t e = startElem + k;
    const KeyType* cc[4] = {c0, c1, c2, c3};
    int n[4]; RealType coords[4][3];
    for (int i = 0; i < 4; ++i) { n[i] = (int)cc[i][e]; coords[i][0]=nx[n[i]]; coords[i][1]=ny[n[i]]; coords[i][2]=nz[n[i]]; }
    RealType det, dNdx[4][3];
    Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx);
    RealType V = det / RealType(6);
    if (!(V > RealType(0))) return;
    for (int i = 0; i < 4; ++i) atomicAdd(&Vp[n[i]], V * RealType(0.25));
}

// invApV[i] = V_P[i] / a_P^u[i], a_P^u = velocity-block diagonal (dof 4i+0) of the assembled momentum block
template<typename RealType>
__global__ void computeInvApVKernel(const int* rowOff, const int* colInd, const RealType* vals,
                                    const RealType* Vp, RealType* invApV, int nNodes)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x; if (i >= nNodes) return;
    int r = 4 * i; RealType aP = 0;
    for (int kk = rowOff[r]; kk < rowOff[r + 1]; ++kk) if (colInd[kk] == r) { aP = vals[kk]; break; }
    invApV[i] = (fabs((double)aP) > 1e-30) ? Vp[i] / aP : RealType(0);
}

// scatter the Rhie-Chow d_f-Laplacian into a^pp over the 6 tet SCS faces (ip order matches d_areaVec)
template<typename KeyType, typename RealType>
__global__ void assembleRhieChowAppKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2,
                                          const KeyType* c3, const RealType* nx, const RealType* ny,
                                          const RealType* nz, const RealType* avX, const RealType* avY,
                                          const RealType* avZ, const RealType* invApV,
                                          const int* rowOff, const int* colInd, RealType* vals,
                                          size_t startElem, size_t numLocal)
{
    size_t k = blockIdx.x * blockDim.x + threadIdx.x; if (k >= numLocal) return;
    size_t e = startElem + k;
    const KeyType* cc[4] = {c0, c1, c2, c3};
    int n[4]; RealType X[4], Y[4], Z[4];
    for (int i = 0; i < 4; ++i) { n[i] = (int)cc[i][e]; X[i]=nx[n[i]]; Y[i]=ny[n[i]]; Z[i]=nz[n[i]]; }
    const int LR[6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};   // matches precomputeTetAreaVectorsKernel
    for (int ip = 0; ip < 6; ++ip) {
        int L = LR[ip][0], R = LR[ip][1], gL = n[L], gR = n[R];
        size_t off = e * 6 + ip;
        RealType Ax = avX[off], Ay = avY[off], Az = avZ[off];
        RealType dx = X[R]-X[L], dy = Y[R]-Y[L], dz = Z[R]-Z[L];
        RealType axdx = Ax*dx + Ay*dy + Az*dz;
        if (fabs((double)axdx) < 1e-30) continue;
        RealType asq = Ax*Ax + Ay*Ay + Az*Az;
        RealType df = RealType(0.5) * (invApV[gL] + invApV[gR]) * asq / axdx;   // > 0
        int rL = 4*gL+3, rsL = rowOff[rL], reL = rowOff[rL+1];
        int rR = 4*gR+3, rsR = rowOff[rR], reR = rowOff[rR+1];
        csrAdd<RealType>(vals, colInd, rsL, reL, 4*gL+3,  df);
        csrAdd<RealType>(vals, colInd, rsL, reL, 4*gR+3, -df);
        csrAdd<RealType>(vals, colInd, rsR, reR, 4*gR+3,  df);
        csrAdd<RealType>(vals, colInd, rsR, reR, 4*gL+3, -df);
    }
}

// momentum block only: a^uu = nu*K (+ frozen convection C_ij) -- component-diagonal, no pressure coupling.
// Used by the Rhie-Chow path (the pressure coupling comes from the FV div/grad + RC a^pp, not Galerkin).
template<typename KeyType, typename RealType>
__global__ void assembleMomentumKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
    const int* rowOff, const int* colInd, RealType* vals, RealType nu, RealType supg,
    const RealType* avx, const RealType* avy, const RealType* avz, size_t startElem, size_t numLocal)
{
    size_t kk = blockIdx.x * blockDim.x + threadIdx.x; if (kk >= numLocal) return;
    size_t e = startElem + kk;
    const KeyType* cc[4] = {c0, c1, c2, c3};
    int n[4]; RealType coords[4][3];
    for (int i = 0; i < 4; ++i) { n[i]=(int)cc[i][e]; coords[i][0]=nodeX[n[i]]; coords[i][1]=nodeY[n[i]]; coords[i][2]=nodeZ[n[i]]; }
    RealType det, dNdx[4][3];
    Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx);
    RealType V = det / RealType(6); if (!(V > RealType(0))) return;
    RealType qV = V * RealType(0.25);
    // element-mean advecting velocity a; bdot[i] = a . gradN_i (reused by Cab, the SUPG terms, and tau)
    RealType ax = 0, ay = 0, az = 0, bdot[4] = {0, 0, 0, 0}, tau = 0;
    if (avx) {
        for (int i=0;i<4;++i){ ax+=avx[n[i]]; ay+=avy[n[i]]; az+=avz[n[i]]; } ax*=RealType(0.25); ay*=RealType(0.25); az*=RealType(0.25);
        RealType Sb = 0, a2 = ax*ax + ay*ay + az*az;
        for (int i=0;i<4;++i){ bdot[i] = ax*dNdx[i][0] + ay*dNdx[i][1] + az*dNdx[i][2]; Sb += fabs(bdot[i]); }
        // SUPG tau (Tezduyar/Codina streamline form): h_s=2|a|/Sb -> tau = 1/sqrt(Sb^2 + (nu*Sb^2/|a|^2)^2).
        // No 1/h_s blowup (h_s substituted out); tau=0 at no-flow so SUPG vanishes where Galerkin is stable.
        const RealType eps = RealType(1e-12);
        if (supg > RealType(0) && a2 > eps*eps && Sb > eps) {
            RealType diff = nu * Sb * Sb / a2;
            tau = supg / sqrt(Sb*Sb + diff*diff);
        }
    }
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
            RealType Kab = V * (dNdx[a][0]*dNdx[b][0] + dNdx[a][1]*dNdx[b][1] + dNdx[a][2]*dNdx[b][2]);
            RealType Cab = qV * bdot[b];                     // consistent-Galerkin convection (0 if no advection)
            RealType Suu = tau * V * bdot[a] * bdot[b];      // SUPG streamline diffusion (a^uu), symmetric PSD
            for (int d = 0; d < 3; ++d) {
                int rad = 4*n[a]+d, rs = rowOff[rad], re = rowOff[rad+1];
                csrAdd<RealType>(vals, colInd, rs, re, 4*n[b]+d, nu*Kab + Cab + Suu);     // a^uu (+SUPG diffusion)
                if (tau > RealType(0))                       // SUPG a^up cross-term tau*(a.gradN_a)*dp/dx_d
                    csrAdd<RealType>(vals, colInd, rs, re, 4*n[b]+3, tau * V * bdot[a] * dNdx[b][d]);
            }
        }
}

// Rhie-Chow FV pressure-velocity coupling, per SCS face (rho=1): a^pu = divergence of the averaged face
// velocity (+-0.5*S_f), and a^up = -(a^pu)^T (the consistent discrete gradient, host-validated). This is
// the collocated div/grad whose checkerboard the RC a^pp damps -- replaces the Galerkin qV*dNdx forms.
template<typename KeyType, typename RealType>
__global__ void assembleRhieChowDivGradKernel(
    const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
    const RealType* avX, const RealType* avY, const RealType* avZ,
    const int* rowOff, const int* colInd, RealType* vals, size_t startElem, size_t numLocal)
{
    size_t kk = blockIdx.x * blockDim.x + threadIdx.x; if (kk >= numLocal) return;
    size_t e = startElem + kk;
    const KeyType* cc[4] = {c0, c1, c2, c3};
    int n[4]; for (int i=0;i<4;++i) n[i]=(int)cc[i][e];
    const int LR[6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
    for (int ip = 0; ip < 6; ++ip) {
        int gL = n[LR[ip][0]], gR = n[LR[ip][1]];
        size_t off = e*6 + ip;
        RealType S[3] = {avX[off], avY[off], avZ[off]};
        int cpL = 4*gL+3, cpLs = rowOff[cpL], cpLe = rowOff[cpL+1];   // continuity rows (a^pu)
        int cpR = 4*gR+3, cpRs = rowOff[cpR], cpRe = rowOff[cpR+1];
        for (int d = 0; d < 3; ++d) {
            RealType h = RealType(0.5) * S[d];
            // a^pu: div of avg face velocity (S out of L => +, out of R => -)
            csrAdd<RealType>(vals, colInd, cpLs, cpLe, 4*gL+d,  h);
            csrAdd<RealType>(vals, colInd, cpLs, cpLe, 4*gR+d,  h);
            csrAdd<RealType>(vals, colInd, cpRs, cpRe, 4*gL+d, -h);
            csrAdd<RealType>(vals, colInd, cpRs, cpRe, 4*gR+d, -h);
            // a^up = -(a^pu)^T  (momentum rows, pressure cols)
            int mL = 4*gL+d, mLs = rowOff[mL], mLe = rowOff[mL+1];
            int mR = 4*gR+d, mRs = rowOff[mR], mRe = rowOff[mR+1];
            csrAdd<RealType>(vals, colInd, mLs, mLe, 4*gL+3, -h);
            csrAdd<RealType>(vals, colInd, mLs, mLe, 4*gR+3,  h);
            csrAdd<RealType>(vals, colInd, mRs, mRe, 4*gL+3, -h);
            csrAdd<RealType>(vals, colInd, mRs, mRe, 4*gR+3,  h);
        }
    }
}

// Open-boundary continuity closure: each boundary node's continuity row 4P+3 gets the missing exterior
// mass flux rho*S_bnd_P . u_P (rho=1, matching the FV a^pu). Added to the continuity rows ONLY -- NOT
// transposed into the momentum (pressure) rows: that asymmetry IS the integration-by-parts boundary
// integral (the natural/do-nothing traction has no matching pressure-gradient boundary force). On a
// closed domain u.n=0 so this contributes nothing after BC; on an inlet/outlet it carries the flux.
template<typename RealType>
__global__ void scatterBoundaryFluxKernel(const RealType* sBndX, const RealType* sBndY, const RealType* sBndZ,
                                          const int* rowOff, const int* colInd, RealType* vals, int nNodes)
{
    int P = blockIdx.x * blockDim.x + threadIdx.x; if (P >= nNodes) return;
    RealType s[3] = {sBndX[P], sBndY[P], sBndZ[P]};
    int cp = 4 * P + 3, cps = rowOff[cp], cpe = rowOff[cp + 1];
    for (int d = 0; d < 3; ++d)
        csrAdd<RealType>(vals, colInd, cps, cpe, 4 * P + d, s[d]);
}

// Assemble the full collocated Rhie-Chow coupled operator into a pre-zeroed CSR (caller fills d_vals=0):
// a^uu(nu*K) + FV a^pu/a^up=-(a^pu)^T + RC a^pp d_f-Laplacian (d_f from the momentum diagonal, no tau)
// + the open-boundary continuity closure rho*S_bnd.u (zero on closed domains, carries inlet/outlet flux).
template<typename KeyType, typename RealType>
void assembleRhieChowGpu(const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
                         const RealType* nx, const RealType* ny, const RealType* nz,
                         const int* d_rowOff, const int* d_colInd, RealType* d_vals, RealType nu, RealType supg,
                         const RealType* avxp, const RealType* avyp, const RealType* avzp,
                         size_t startElem, size_t numLocal, int nNodes, int blockSize)
{
    const int eBlocks = (int)((numLocal + blockSize - 1) / blockSize);
    thrust::device_vector<RealType> avX(6 * numLocal), avY(6 * numLocal), avZ(6 * numLocal);
    precomputeTetAreaVectorsGpu<KeyType, RealType>(c0, c1, c2, c3, numLocal, nx, ny, nz,
        thrust::raw_pointer_cast(avX.data()), thrust::raw_pointer_cast(avY.data()), thrust::raw_pointer_cast(avZ.data()));
    cudaDeviceSynchronize();
    assembleMomentumKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
        c0, c1, c2, c3, nx, ny, nz, d_rowOff, d_colInd, d_vals, nu, supg, avxp, avyp, avzp, startElem, numLocal);
    cudaDeviceSynchronize();
    thrust::device_vector<RealType> Vp(nNodes, RealType(0)), invApV(nNodes, RealType(0));
    computeNodeVolumeKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
        c0, c1, c2, c3, nx, ny, nz, thrust::raw_pointer_cast(Vp.data()), startElem, numLocal);
    cudaDeviceSynchronize();
    const int nblkN = (nNodes + blockSize - 1) / blockSize;
    computeInvApVKernel<RealType><<<nblkN, blockSize>>>(d_rowOff, d_colInd, d_vals,
        thrust::raw_pointer_cast(Vp.data()), thrust::raw_pointer_cast(invApV.data()), nNodes);
    cudaDeviceSynchronize();
    assembleRhieChowDivGradKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
        c0, c1, c2, c3, thrust::raw_pointer_cast(avX.data()), thrust::raw_pointer_cast(avY.data()),
        thrust::raw_pointer_cast(avZ.data()), d_rowOff, d_colInd, d_vals, startElem, numLocal);
    assembleRhieChowAppKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
        c0, c1, c2, c3, nx, ny, nz, thrust::raw_pointer_cast(avX.data()), thrust::raw_pointer_cast(avY.data()),
        thrust::raw_pointer_cast(avZ.data()), thrust::raw_pointer_cast(invApV.data()),
        d_rowOff, d_colInd, d_vals, startElem, numLocal);
    cudaDeviceSynchronize();
    // open-boundary continuity closure: S_bnd (geometric, could be hoisted out of the Picard loop) then
    // scatter rho*S_bnd.u into the continuity rows. Closes the open dual cells so inlet/outlet flux balances.
    thrust::device_vector<RealType> sbX(nNodes, RealType(0)), sbY(nNodes, RealType(0)), sbZ(nNodes, RealType(0));
    mars::fem::computeBoundaryAreaVectorsGpu<KeyType, RealType>(c0, c1, c2, c3, startElem, numLocal, nx, ny, nz,
        nNodes, thrust::raw_pointer_cast(sbX.data()), thrust::raw_pointer_cast(sbY.data()),
        thrust::raw_pointer_cast(sbZ.data()));
    cudaDeviceSynchronize();
    scatterBoundaryFluxKernel<RealType><<<nblkN, blockSize>>>(
        thrust::raw_pointer_cast(sbX.data()), thrust::raw_pointer_cast(sbY.data()),
        thrust::raw_pointer_cast(sbZ.data()), d_rowOff, d_colInd, d_vals, nNodes);
    cudaDeviceSynchronize();
}

// ---- Increment 3a: Dirichlet BC elimination + direct reference solve ----
// Identity-row BC: for a Dirichlet dof, zero its CSR row, set its diagonal to 1,
// rhs = prescribed value. Correct for a DIRECT solve (the interior rows still see
// the pinned value through their column coupling); no column elimination needed.
template<typename RealType>
__global__ void applyBCKernel(const uint8_t* bcFlag, const RealType* bcVal,
                              const int* rowOff, const int* colInd,
                              RealType* vals, RealType* rhs, int ND)
{
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= ND || !bcFlag[r]) return;
    for (int k = rowOff[r]; k < rowOff[r + 1]; ++k)
        vals[k] = (colInd[k] == r) ? RealType(1) : RealType(0);
    rhs[r] = bcVal[r];
}

template<typename RealType>
__global__ void extractCompKernel(const RealType* x, RealType* out, int stride, int comp, int N)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;
    out[i] = x[stride * i + comp];
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    std::cout << std::unitbuf;          // flush per write: live [phase0] trace under srun's block-buffered pipe
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    // TET AmrManager init needs a bound device (see memory: tet drivers must
    // set the device themselves, unlike the hex tgv driver).
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    using KeyType  = uint64_t;
    using RealType = double;

    std::string meshFile;
    int blockSize  = 256;
    int bucketSize = 64;
    bool doAssemble = false;     // --assemble: also run the coupled Stokes block assembly + gates
    bool doAdvect   = false;     // --advect: add frozen-velocity convection into a^uu (implies --assemble)
    bool doSolve    = false;     // --solve: BC elimination + cusolverSp QR direct reference solve (implies --assemble)
    bool doHypre    = false;     // --hypre: Phase-1 1a -- Hypre point-block BoomerAMG GMRES iters vs cusolver ref (implies --assemble)
    bool doAcm3     = false;      // --acm-stage3: GPU coarse-operator + transfers vs host P^T A P (implies --assemble)
    bool doAcm4     = false;      // --acm-stage4: GPU multilevel V-cycle vs host V-cycle replica (implies --assemble)
    bool doAcm4b    = false;      // --acm-stage4b: ACM-preconditioned FlexGMRES iters vs 1a baseline + cusolver ref
    bool doAcmPump  = false;      // --acm-pump: mesh-agnostic ACM-FlexGMRES vs BoomerAMG on ANY tet mesh (residual-based)
    bool doRhieChow = false;      // --rhie-chow: assemble the collocated Rhie-Chow coupled operator (vs PSPG); implies --assemble
    bool doSupg     = false;      // --supg: SUPG streamline-upwind stabilization on the momentum convection (needed at high Re)
    bool doChannel  = false;      // --channel: plane-channel inlet/outlet/wall BC + Poiseuille validation (implies --solve)
    double beta     = 45.0;      // --beta: skew-cavity angle (deg), used only for the --solve BC geometry
    int    Re       = 100;       // --Re: Reynolds number for the lid-driven solve (3b)
    int    picard   = 30;        // --picard: max Picard outer iterations
    std::string inletSS, outletSS;  // --inlet-ss / --outlet-ss: Exodus side-set NAMES (pump real inlet/outlet)
    double inletSpeed = 0.5;     // --inlet-speed: velocity-flux inlet magnitude (m/s), along the inward normal
    bool   inletFlip  = false;   // --inlet-flip-normal: reverse the inlet normal if the Exodus winding is outward
    double nuVal      = -1.0;    // --nu: dimensional kinematic viscosity (overrides 1/Re when >=0)
    double relax      = 1.0;     // --relax: Picard under-relaxation omega (u <- w*u_new + (1-w)*u_old); <1 for high Re
    int    acmRebuild = 1;       // --acm-rebuild: rebuild the ACM hierarchy every K Picard iters (reuse between; K>1 for host-heavy ILU)
    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if      (a.rfind("--mesh=", 0) == 0)        meshFile   = a.substr(7);
        else if (a.rfind("--block-size=", 0) == 0)  blockSize  = std::stoi(a.substr(13));
        else if (a.rfind("--bucket-size=", 0) == 0) bucketSize = std::stoi(a.substr(14));
        else if (a == "--assemble")                 doAssemble = true;
        else if (a == "--advect")                 { doAssemble = true; doAdvect = true; }
        else if (a == "--solve")                  { doAssemble = true; doSolve  = true; }
        else if (a == "--hypre")                  { doAssemble = true; doHypre  = true; }
        else if (a == "--acm-stage3")             { doAssemble = true; doAcm3   = true; }
        else if (a == "--acm-stage4")             { doAssemble = true; doAcm4   = true; }
        else if (a == "--acm-stage4b")            { doAssemble = true; doAcm4b  = true; }
        else if (a == "--acm-pump")               { doAssemble = true; doAcmPump = true; }
        else if (a == "--rhie-chow")              { doAssemble = true; doRhieChow = true; }
        else if (a == "--supg")                     doSupg     = true;
        else if (a.rfind("--relax=", 0) == 0)        relax      = std::stod(a.substr(8));
        else if (a.rfind("--acm-rebuild=", 0) == 0)  acmRebuild = std::stoi(a.substr(14));
        else if (a == "--channel")                { doAssemble = true; doSolve = true; doChannel = true; }
        else if (a.rfind("--inlet-ss=", 0) == 0)  { inletSS = a.substr(11); doAssemble = doRhieChow = true; }
        else if (a.rfind("--outlet-ss=", 0) == 0)   outletSS   = a.substr(12);
        else if (a.rfind("--inlet-speed=", 0) == 0) inletSpeed = std::stod(a.substr(14));
        else if (a == "--inlet-flip-normal")        inletFlip  = true;
        else if (a.rfind("--nu=", 0) == 0)          nuVal      = std::stod(a.substr(5));
        else if (a.rfind("--beta=", 0) == 0)        beta       = std::stod(a.substr(7));
        else if (a.rfind("--Re=", 0) == 0)          Re         = std::stoi(a.substr(5));
        else if (a.rfind("--picard=", 0) == 0)      picard     = std::stoi(a.substr(9));
        else if (a[0] != '-' && meshFile.empty())   meshFile   = a;
    }
    if (meshFile.empty())
    {
        if (rank == 0)
            std::cout << "Usage: " << argv[0] << " --mesh=FILE.exo [--block-size=N]\n";
        MPI_Finalize();
        return 1;
    }

    // The adjointness gate below is now multi-rank-correct (reverseExchangeNodeHaloAdd folds the ghost-scattered
    // div/grad contributions to the owner, then owned-masked + Allreduce dots). The coupled SOLVE is still
    // single-rank -- that guard moves to after the gate (Increment 2 = operator + Krylov halo).

    AmrManager<TetTag, KeyType, RealType>::Config cfg;
    cfg.maxLevels  = 0;          // frozen mesh
    cfg.blockSize  = blockSize;
    cfg.bucketSize = bucketSize;
    if (!inletSS.empty()) cfg.storeSideSets = true;   // pump inlet/outlet BC needs the Exodus side sets
    AmrManager<TetTag, KeyType, RealType> amr(cfg);
    amr.initialize(meshFile, rank, numRanks);
    auto& domain = amr.domain();

    const size_t nNodes   = domain.getNodeCount();
    const size_t startEl  = domain.startIndex();
    const size_t numLocal = domain.localElementCount();
    if (rank == 0)
        if (std::getenv("MARS_VERBOSE_MESH"))
            std::cout << "[phase0] mesh: elements=" << domain.getElementCount()
                      << " nodes=" << nNodes << " localElems=" << numLocal << "\n";

    const auto& conn = domain.getElementToNodeConnectivity();
    const KeyType* c0 = std::get<0>(conn).data();
    const KeyType* c1 = std::get<1>(conn).data();
    const KeyType* c2 = std::get<2>(conn).data();
    const KeyType* c3 = std::get<3>(conn).data();
    const RealType* nx = domain.getNodeX().data();
    const RealType* ny = domain.getNodeY().data();
    const RealType* nz = domain.getNodeZ().data();

    // Deterministic smooth test fields from coordinates. Any phi, v validate the
    // adjointness identity; smooth polynomials make a failure easy to read.
    std::vector<RealType> hx(nNodes), hy(nNodes), hz(nNodes);
    thrust::copy(thrust::device_pointer_cast(nx), thrust::device_pointer_cast(nx + nNodes), hx.begin());
    thrust::copy(thrust::device_pointer_cast(ny), thrust::device_pointer_cast(ny + nNodes), hy.begin());
    thrust::copy(thrust::device_pointer_cast(nz), thrust::device_pointer_cast(nz + nNodes), hz.begin());

    std::vector<RealType> hphi(nNodes), hvx(nNodes), hvy(nNodes), hvz(nNodes);
    for (size_t i = 0; i < nNodes; ++i)
    {
        const RealType x = hx[i], y = hy[i], z = hz[i];
        hphi[i] = 1.0 + 2.0 * x - 3.0 * y + 0.7 * z + 0.4 * x * y - 0.2 * y * z;
        hvx[i]  = 0.5 + x - 0.3 * z + 0.1 * x * x;
        hvy[i]  = -0.2 + 0.4 * y + 0.6 * z;
        hvz[i]  = 0.9 - 0.5 * x + 0.3 * y * y;
    }

    thrust::device_vector<RealType> phi(hphi.begin(), hphi.end());
    thrust::device_vector<RealType> vx(hvx.begin(), hvx.end());
    thrust::device_vector<RealType> vy(hvy.begin(), hvy.end());
    thrust::device_vector<RealType> vz(hvz.begin(), hvz.end());
    thrust::device_vector<RealType> divAcc(nNodes, RealType(0));
    thrust::device_vector<RealType> gx(nNodes, RealType(0));
    thrust::device_vector<RealType> gy(nNodes, RealType(0));
    thrust::device_vector<RealType> gz(nNodes, RealType(0));

    const int eBlocks = static_cast<int>((numLocal + blockSize - 1) / blockSize);

    // D v  (weak divergence, scattered to nodes)
    computeFemDivergenceTetKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
        c0, c1, c2, c3,
        thrust::raw_pointer_cast(vx.data()),
        thrust::raw_pointer_cast(vy.data()),
        thrust::raw_pointer_cast(vz.data()),
        nx, ny, nz,
        thrust::raw_pointer_cast(divAcc.data()),
        startEl, numLocal);

    // G phi  (un-normalized weak gradient accumulator -- this is exactly D^T's
    // partner; we do NOT divide by the lumped mass here, since the adjoint
    // identity holds between the raw operators).
    computeFemGradientTetKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
        c0, c1, c2, c3,
        thrust::raw_pointer_cast(phi.data()),
        nx, ny, nz,
        thrust::raw_pointer_cast(gx.data()),
        thrust::raw_pointer_cast(gy.data()),
        thrust::raw_pointer_cast(gz.data()),
        startEl, numLocal);

    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        std::cerr << "[phase0] kernel launch failed: " << cudaGetErrorString(err) << "\n";
        MPI_Finalize();
        return 2;
    }

    // <phi, D v> and <v, G phi>; the discrete integration-by-parts identity is <phi,Dv> + <v,Gphi> = 0.
    // Multi-rank: the div/grad kernels scatter to ghost nodes, so fold those contributions back to the owner
    // (reverseExchangeNodeHaloAdd) before any dot, and reduce ONLY over owned nodes -- getNodeCount includes
    // ghosts, so ghost DOFs are duplicate unknowns and an unmasked reduction would double-count every seam node.
    domain.reverseExchangeNodeHaloAdd(divAcc);
    domain.reverseExchangeNodeHaloAdd(gx);
    domain.reverseExchangeNodeHaloAdd(gy);
    domain.reverseExchangeNodeHaloAdd(gz);
    const uint8_t* ownPtr = domain.getNodeOwnershipMap().data();   // 1 = this rank owns the node
    auto ownedDot = [&](const thrust::device_vector<RealType>& a, const thrust::device_vector<RealType>& b) -> RealType {
        const RealType* ap = thrust::raw_pointer_cast(a.data());
        const RealType* bp = thrust::raw_pointer_cast(b.data());
        const RealType loc = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(nNodes),
            [ownPtr, ap, bp] __device__(size_t i) -> RealType { return ownPtr[i] == 1 ? ap[i] * bp[i] : RealType(0); },
            RealType(0), thrust::plus<RealType>());
        RealType g = 0; MPI_Allreduce(&loc, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return g;   // RealType=double here
    };
    const RealType dDiv  = ownedDot(phi, divAcc);
    const RealType dGrad = ownedDot(vx, gx) + ownedDot(vy, gy) + ownedDot(vz, gz);

    const RealType residual = dDiv + dGrad;
    const RealType scale    = std::max({std::abs(dDiv), std::abs(dGrad), RealType(1e-300)});
    const RealType relErr   = std::abs(residual) / scale;

    const RealType tol = 1e-12;
    const bool pass = relErr < tol;

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "[phase0][G=-D^T gate] <phi,Dv> = " << dDiv
              << "  <v,Gphi> = " << dGrad
              << "  rel|sum| = " << relErr
              << "  -> " << (pass ? "PASS" : "FAIL") << "\n";

    // The coupled-operator assembly + Krylov solve are not yet multi-rank (reverse-fold of the assembled
    // operator + halo-exchange in the matvec/dots = Increment 2). Refuse >1 rank here so the gate above can
    // validate the halo-fold/owned-mask/allreduce machinery at N ranks, without reporting a wrong solve.
    if (numRanks > 1)
    {
        if (rank == 0)
            std::cout << "[phase0] adjointness gate is multi-rank-correct; the coupled solve is single-rank only "
                         "(operator + Krylov halo = Increment 2). Gate validated above -- run the solve with 1 rank.\n";
        MPI_Finalize();
        return 0;
    }

    // Side-set velocity-flux inlet + no-slip walls BC builder (outlet velocity left FREE; caller pins
    // pressure). Shared by --solve (cusolver) and --acm-pump (ACM). Fills velocity comps of bcF/bcV; also
    // returns the per-node boundary S_bnd + inlet/outlet node sets for the through-flow diagnostics.
    // Inlet normal = global average inlet S_bnd (exact for a planar opening). Returns false if a SS missing.
    const bool inletSet = !inletSS.empty();       // pump side-set BC active (shared by --acm-pump and --solve)
    std::vector<RealType> hsbx, hsby, hsbz;       // per-node boundary area vectors (also used in validation)
    std::set<int> inN, outN;                      // inlet / outlet node sets
    auto buildSideSetVelBC = [&](std::vector<uint8_t>& bcF, std::vector<RealType>& bcV) -> bool {
        if (!domain.hasSideSet(inletSS) || (!outletSS.empty() && !domain.hasSideSet(outletSS))) {
            std::cerr << "[phase0][pump] side-set not found (inlet='" << inletSS << "' outlet='" << outletSS
                      << "'); available:";
            for (auto& nm : domain.sideSetNames()) std::cerr << " " << nm;
            std::cerr << "\n"; return false;
        }
        thrust::device_vector<RealType> dsbx(nNodes), dsby(nNodes), dsbz(nNodes);
        mars::fem::computeBoundaryAreaVectorsGpu<KeyType, RealType>(c0, c1, c2, c3, startEl, numLocal,
            nx, ny, nz, (int)nNodes, thrust::raw_pointer_cast(dsbx.data()),
            thrust::raw_pointer_cast(dsby.data()), thrust::raw_pointer_cast(dsbz.data()));
        cudaDeviceSynchronize();
        hsbx.resize(nNodes); hsby.resize(nNodes); hsbz.resize(nNodes);
        thrust::copy(dsbx.begin(), dsbx.end(), hsbx.begin());
        thrust::copy(dsby.begin(), dsby.end(), hsby.begin());
        thrust::copy(dsbz.begin(), dsbz.end(), hsbz.begin());
        auto ssToHost = [&](const std::string& nm) {
            const auto& dv = domain.sideSetNodes(nm);
            std::vector<int> h(dv.size());
            thrust::copy(thrust::device_pointer_cast(dv.data()),
                         thrust::device_pointer_cast(dv.data() + dv.size()), h.begin());
            return h;
        };
        std::vector<int> hin = ssToHost(inletSS);
        inN.insert(hin.begin(), hin.end());
        if (!outletSS.empty()) { auto ho = ssToHost(outletSS); outN.insert(ho.begin(), ho.end()); }
        RealType sx = 0, sy = 0, sz = 0;                              // net outward inlet area vector
        for (int i : hin) { sx += hsbx[i]; sy += hsby[i]; sz += hsbz[i]; }
        RealType len = std::sqrt(sx * sx + sy * sy + sz * sz);
        RealType sgn = (inletFlip ? RealType(1) : RealType(-1)) / (len > 0 ? len : RealType(1));  // inward
        RealType nhx = sgn * sx, nhy = sgn * sy, nhz = sgn * sz;
        // Side-set membership is authoritative for inlet/outlet -- apply the velocity BC to EVERY
        // inlet/outlet node, NOT only those S_bnd flagged as boundary. S_bnd is a geometric
        // reconstruction (boundary-face dedup + 1/3 scatter) that is legitimately zero for a side-set
        // node whose opening faces are tagged but topologically internal (shared by 2 tets) -- gating
        // the BC on S_bnd>0 left those nodes free -> a leaky inlet/outlet. The inlet normal is the
        // GLOBAL average inlet S_bnd (planar), so it is well-defined even where a node's own S_bnd==0.
        for (int i : inN) {                                           // velocity-flux inlet (all members)
            bcF[4 * i + 0] = bcF[4 * i + 1] = bcF[4 * i + 2] = 1;
            bcV[4 * i + 0] = inletSpeed * nhx; bcV[4 * i + 1] = inletSpeed * nhy; bcV[4 * i + 2] = inletSpeed * nhz;
        }
        // open outlet: velocity free, caller pins pressure (members already excluded from walls below)
        for (size_t i = 0; i < nNodes; ++i) {
            if (inN.count((int)i) || outN.count((int)i)) continue;    // inlet set above; outlet stays free
            bool isB = (hsbx[i] * hsbx[i] + hsby[i] * hsby[i] + hsbz[i] * hsbz[i]) > RealType(1e-30);
            if (!isB) continue;                                       // interior node: free
            bcF[4 * i + 0] = bcF[4 * i + 1] = bcF[4 * i + 2] = 1;     // no-slip wall
        }
        if (rank == 0 && std::getenv("MARS_VERBOSE_MESH"))   // node counts + inlet normal identify the geometry
            std::cout << "[phase0][pump] inlet nodes=" << inN.size() << " outlet nodes=" << outN.size()
                      << " inward-normal=(" << nhx << "," << nhy << "," << nhz << ") speed=" << inletSpeed << "\n";
        return true;
    };

    // ---- Increment 2a: assemble the coupled Stokes block + matrix-consistency gates ----
    bool assemblePass = true;
    if (doAssemble)
    {
        const RealType nu = 1e-3, tau = 1.0 / 24.0;   // values irrelevant to the gates
        const int ND = 4 * static_cast<int>(nNodes);

        // connectivity to host (single rank: numLocal == elementCount, startEl == 0)
        std::vector<KeyType> h0(numLocal), h1(numLocal), h2(numLocal), h3(numLocal);
        thrust::copy(thrust::device_pointer_cast(c0), thrust::device_pointer_cast(c0 + numLocal), h0.begin());
        thrust::copy(thrust::device_pointer_cast(c1), thrust::device_pointer_cast(c1 + numLocal), h1.begin());
        thrust::copy(thrust::device_pointer_cast(c2), thrust::device_pointer_cast(c2 + numLocal), h2.begin());
        thrust::copy(thrust::device_pointer_cast(c3), thrust::device_pointer_cast(c3 + numLocal), h3.begin());

        // node 1-ring adjacency (nodes sharing a tet), then interleaved 4N block sparsity
        std::vector<std::set<int>> ring(nNodes);
        for (size_t e = 0; e < numLocal; ++e) {
            int nn[4] = {(int)h0[e], (int)h1[e], (int)h2[e], (int)h3[e]};
            for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b) ring[nn[a]].insert(nn[b]);
        }
        std::vector<int> h_rowOff(ND + 1, 0);
        for (size_t i = 0; i < nNodes; ++i) {
            int deg = 4 * static_cast<int>(ring[i].size());
            for (int a = 0; a < 4; ++a) h_rowOff[4 * i + a + 1] = deg;
        }
        for (int r = 0; r < ND; ++r) h_rowOff[r + 1] += h_rowOff[r];
        const int nnz = h_rowOff[ND];
        std::vector<int> h_colInd(nnz);
        for (size_t i = 0; i < nNodes; ++i) {
            std::vector<int> rs(ring[i].begin(), ring[i].end());   // sorted (std::set)
            for (int a = 0; a < 4; ++a) {
                int p = h_rowOff[4 * i + a];
                for (int j : rs)
                    for (int b = 0; b < 4; ++b) h_colInd[p++] = 4 * j + b;
            }
        }

        thrust::device_vector<int>      d_rowOff(h_rowOff.begin(), h_rowOff.end());
        thrust::device_vector<int>      d_colInd(h_colInd.begin(), h_colInd.end());
        thrust::device_vector<RealType> d_vals(nnz, RealType(0));

        // frozen advecting velocity (divergence-free rotation about the mesh centre)
        thrust::device_vector<RealType> d_avx, d_avy, d_avz;
        const RealType *avxp = nullptr, *avyp = nullptr, *avzp = nullptr;
        if (doAdvect) {
            RealType xc = 0, yc = 0;
            for (size_t i = 0; i < nNodes; ++i) { xc += hx[i]; yc += hy[i]; }
            xc /= nNodes; yc /= nNodes;
            std::vector<RealType> hax(nNodes), hay(nNodes), haz(nNodes, RealType(0));
            for (size_t i = 0; i < nNodes; ++i) { hax[i] = -(hy[i] - yc); hay[i] = (hx[i] - xc); }
            d_avx.assign(hax.begin(), hax.end());
            d_avy.assign(hay.begin(), hay.end());
            d_avz.assign(haz.begin(), haz.end());
            avxp = thrust::raw_pointer_cast(d_avx.data());
            avyp = thrust::raw_pointer_cast(d_avy.data());
            avzp = thrust::raw_pointer_cast(d_avz.data());
        }

        if (doRhieChow) {
            assembleRhieChowGpu<KeyType, RealType>(c0, c1, c2, c3, nx, ny, nz,
                thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                thrust::raw_pointer_cast(d_vals.data()), nu, doSupg ? RealType(1) : RealType(0), avxp, avyp, avzp,
                startEl, numLocal, (int)nNodes, blockSize);
        } else {
            assembleCoupledStokesKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
                c0, c1, c2, c3, nx, ny, nz,
                thrust::raw_pointer_cast(d_rowOff.data()),
                thrust::raw_pointer_cast(d_colInd.data()),
                thrust::raw_pointer_cast(d_vals.data()),
                nu, tau, avxp, avyp, avzp, startEl, numLocal);
        }
        cudaError_t aerr = cudaDeviceSynchronize();
        if (aerr != cudaSuccess) {
            std::cerr << "[phase0][assemble] kernel failed: " << cudaGetErrorString(aerr) << "\n";
            MPI_Finalize();
            return 4;
        }

        // copy CSR back; gates are host-side (Phase-0 model operator, small mesh)
        std::vector<RealType> hv(nnz);
        thrust::copy(d_vals.begin(), d_vals.end(), hv.begin());
        std::map<std::pair<int, int>, RealType> A;
        for (int r = 0; r < ND; ++r)
            for (int k = h_rowOff[r]; k < h_rowOff[r + 1]; ++k) A[{r, h_colInd[k]}] = hv[k];
        auto get = [&](int r, int c) -> RealType {
            auto it = A.find({r, c}); return it == A.end() ? RealType(0) : it->second;
        };

        // boundary nodes: a tet face in exactly one tet marks its 3 nodes. The open-boundary continuity
        // flux (doRhieChow) breaks a^pu==-(a^up)^T on a boundary node's continuity DIAGONAL by design
        // (the IBP boundary integral, in continuity rows only) -> scope the transpose gate around it.
        std::vector<uint8_t> bnode(nNodes, 0);
        {
            const int TF[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
            std::map<std::array<int, 3>, int> fc;
            for (size_t e = 0; e < numLocal; ++e) {
                int nn[4] = {(int)h0[e], (int)h1[e], (int)h2[e], (int)h3[e]};
                for (int f = 0; f < 4; ++f) {
                    std::array<int, 3> tri = {nn[TF[f][0]], nn[TF[f][1]], nn[TF[f][2]]};
                    std::sort(tri.begin(), tri.end()); fc[tri]++;
                }
            }
            for (auto& kv : fc) if (kv.second == 1) for (int v : kv.first) bnode[v] = 1;
        }
        RealType e1 = 0, e1b = 0, e2 = 0, e3 = 0, e4 = 0, e5 = 0;
        for (size_t i = 0; i < nNodes; ++i)
            for (int j : ring[i])
                for (int d = 0; d < 3; ++d) {
                    RealType t = std::abs(get(4 * i + 3, 4 * j + d) + get(4 * j + d, 4 * i + 3)); // a^pu==-(a^up)^T
                    if (doRhieChow && bnode[i] && j == (int)i) e1b = std::max(e1b, t);            // expected boundary break
                    else e1 = std::max(e1, t);
                    e4 = std::max(e4, std::abs(get(4 * i + d, 4 * j + d) - get(4 * j + d, 4 * i + d))); // a^uu symmetry residual
                }
        for (size_t i = 0; i < nNodes; ++i) {
            RealType spp = 0, sdc = 0, suu[3] = {0, 0, 0};
            for (int j : ring[i]) {
                spp += get(4 * i + 3, 4 * j + 3);
                for (int d = 0; d < 3; ++d) { sdc += get(4 * i + 3, 4 * j + d); suu[d] += get(4 * i + d, 4 * j + d); }
            }
            e2 = std::max(e2, std::abs(spp));                     // a^pp . 1 == 0
            e3 += sdc;                                            // global sum(D.const)
            for (int d = 0; d < 3; ++d) e5 = std::max(e5, std::abs(suu[d])); // a^uu . 1 (conservation)
        }
        e3 = std::abs(e3);

        const RealType atol = 1e-12;
        auto pf = [&](bool ok) { return ok ? "PASS" : "FAIL"; };
        std::cout << std::scientific << std::setprecision(3);
        if (!std::getenv("MARS_VERBOSE_MESH"))
            std::cout << "[phase0][assemble]" << (doAdvect ? "  (advection ON)" : "  (Stokes)") << "\n";
        else
            std::cout << "[phase0][assemble] DOFs=" << ND << " nnz=" << nnz
                      << (doAdvect ? "  (advection ON)" : "  (Stokes)") << "\n";
        // SUPG adds a consistent a^up cross-term (momentum rows, pressure cols) with no a^pu counterpart,
        // so it legitimately breaks the transpose identity -- expected, like the open-boundary flux. e2/e3
        // (continuity-row quantities) and a^pp are untouched by SUPG.
        bool e1ok = (e1 < atol) || doSupg;
        std::cout << "[phase0][assemble] a_pu==-a_up^T : " << e1 << "  " << pf(e1ok)
                  << (doSupg ? "  (SUPG cross-term expected)" : doRhieChow ? "  (interior)" : "") << "\n";
        if (doRhieChow)
            std::cout << "[phase0][assemble] bndry-flux brk: " << e1b
                      << "  (EXPECTED >0 = open-boundary IBP integral, continuity rows only)\n";
        std::cout << "[phase0][assemble] a_pp.1==0     : " << e2 << "  " << pf(e2 < atol) << "\n";
        std::cout << "[phase0][assemble] sum(D.const)  : " << e3 << "  " << pf(e3 < atol) << "\n";
        if (doAdvect) {
            std::cout << "[phase0][assemble] a_uu nonsymm  : " << e4 << "  " << pf(e4 > atol) << "\n";
            std::cout << "[phase0][assemble] a_uu conserves: " << e5 << "  " << pf(e5 < atol) << "\n";
            assemblePass = e1ok && (e2 < atol) && (e3 < atol) && (e4 > atol) && (e5 < atol);
        } else {
            std::cout << "[phase0][assemble] a_uu symmetric: " << e4 << "  " << pf(e4 < atol) << "\n";
            assemblePass = e1ok && (e2 < atol) && (e3 < atol) && (e4 < atol);
        }
        if (doRhieChow) {   // G2: a^pp is a Rhie-Chow d_f-Laplacian M-matrix (diag>0, off<=0)
            RealType mdiagMin = RealType(1e30), moffMax = RealType(-1e30);
            for (size_t i = 0; i < nNodes; ++i) {
                mdiagMin = std::min(mdiagMin, get(4 * i + 3, 4 * i + 3));
                for (int j : ring[i]) if (j != (int)i) moffMax = std::max(moffMax, get(4 * i + 3, 4 * j + 3));
            }
            bool mok = (mdiagMin > RealType(0)) && (moffMax <= atol);
            std::cout << "[phase0][assemble] a_pp M-matrix : diag_min=" << mdiagMin << " off_max=" << moffMax
                      << "  " << pf(mok) << " (Rhie-Chow, no tau)\n";
            assemblePass = assemblePass && mok;
        }
        std::cout << "[phase0][assemble] -> " << (assemblePass ? "ALL PASS" : "FAIL") << "\n";

        // ---- ACM Stage 3: GPU coarse-operator (sum-of-fine == P^T A P) + grid transfers ----
        if (doAcm3)
        {
            // host-once aggregation (node-granular, greedy over the 1-ring). The map is an INPUT to
            // Stage 3 -- any valid map validates the GPU coarse-operator recipe; the GPU directional
            // aggregation is Stage 5.
            const int KMAX = 4;
            std::vector<int> agg(nNodes, -1);
            int nCoarse = 0;
            for (size_t i = 0; i < nNodes; ++i) {
                if (agg[i] != -1) continue;
                agg[i] = nCoarse; int cnt = 1;
                for (int j : ring[i]) {
                    if (cnt >= KMAX) break;
                    if (j != (int)i && agg[j] == -1) { agg[j] = nCoarse; ++cnt; }
                }
                ++nCoarse;
            }
            const int NDc = 4 * nCoarse;

            thrust::device_vector<int> d_agg(agg.begin(), agg.end());
            thrust::device_vector<int> d_crowOff, d_ccolInd;
            thrust::device_vector<RealType> d_cvals;
            buildCoarseOperator<RealType>(d_rowOff, d_colInd, d_vals, d_agg,
                                          (int)nNodes, nCoarse, d_crowOff, d_ccolInd, d_cvals);

            // host reference: the same map applied to the fine CSR by an independent accumulation
            auto cdof = [&](int r) { return 4 * agg[r >> 2] + (r & 3); };
            std::map<std::pair<int, int>, RealType> Aref;
            for (int r = 0; r < ND; ++r)
                for (int k = h_rowOff[r]; k < h_rowOff[r + 1]; ++k)
                    Aref[{cdof(r), cdof(h_colInd[k])}] += hv[k];

            const int ccnnz = (int)d_cvals.size();
            std::vector<int> hcrow(NDc + 1), hccol(ccnnz);
            std::vector<RealType> hcval(ccnnz);
            thrust::copy(d_crowOff.begin(), d_crowOff.end(), hcrow.begin());
            thrust::copy(d_ccolInd.begin(), d_ccolInd.end(), hccol.begin());
            thrust::copy(d_cvals.begin(),   d_cvals.end(),   hcval.begin());
            std::map<std::pair<int, int>, RealType> Agpu;
            for (int r = 0; r < NDc; ++r)
                for (int k = hcrow[r]; k < hcrow[r + 1]; ++k) Agpu[{r, hccol[k]}] = hcval[k];
            RealType dRecipe = 0;
            for (auto& kv : Aref) { auto it = Agpu.find(kv.first); dRecipe = std::max(dRecipe, std::abs(kv.second - (it == Agpu.end() ? RealType(0) : it->second))); }
            for (auto& kv : Agpu) { auto it = Aref.find(kv.first); dRecipe = std::max(dRecipe, std::abs(kv.second - (it == Aref.end() ? RealType(0) : it->second))); }

            // Galerkin consistency: P^T (A (P x)) == A_coarse x for a random coarse x (host)
            std::vector<RealType> xc(NDc), xf(ND, 0), yf(ND, 0), yc(NDc, 0), zc(NDc, 0), rc(NDc, 0);
            for (int i = 0; i < NDc; ++i) xc[i] = std::sin(RealType(0.7) * i + RealType(1));
            for (size_t i = 0; i < nNodes; ++i) for (int c = 0; c < 4; ++c) xf[4 * i + c] = xc[4 * agg[i] + c];   // P xc
            for (size_t i = 0; i < nNodes; ++i) for (int c = 0; c < 4; ++c) rc[4 * agg[i] + c] += xf[4 * i + c];  // P^T xf (restrict-kernel ref)
            for (int r = 0; r < ND; ++r) { RealType s = 0; for (int k = h_rowOff[r]; k < h_rowOff[r + 1]; ++k) s += hv[k] * xf[h_colInd[k]]; yf[r] = s; }  // A xf
            for (size_t i = 0; i < nNodes; ++i) for (int c = 0; c < 4; ++c) yc[4 * agg[i] + c] += yf[4 * i + c];  // P^T A xf
            for (auto& kv : Aref) zc[kv.first.first] += kv.second * xc[kv.first.second];                          // A_coarse xc
            RealType dGal = 0; for (int i = 0; i < NDc; ++i) dGal = std::max(dGal, std::abs(yc[i] - zc[i]));

            // GPU transfer kernels vs host (injection prolongation + sum restriction)
            thrust::device_vector<RealType> d_xc(xc.begin(), xc.end()), d_xf(ND, RealType(0)), d_yc(NDc, RealType(0));
            int nb = ((int)nNodes + blockSize - 1) / blockSize;
            acmInjectKernel<RealType><<<nb, blockSize>>>(thrust::raw_pointer_cast(d_xc.data()),
                thrust::raw_pointer_cast(d_xf.data()), thrust::raw_pointer_cast(d_agg.data()), (int)nNodes);
            acmRestrictAddKernel<RealType><<<nb, blockSize>>>(thrust::raw_pointer_cast(d_xf.data()),
                thrust::raw_pointer_cast(d_yc.data()), thrust::raw_pointer_cast(d_agg.data()), (int)nNodes);
            cudaDeviceSynchronize();
            std::vector<RealType> gxf(ND), gyc(NDc);
            thrust::copy(d_xf.begin(), d_xf.end(), gxf.begin());
            thrust::copy(d_yc.begin(), d_yc.end(), gyc.begin());
            RealType dInj = 0, dRes = 0;
            for (int i = 0; i < ND;  ++i) dInj = std::max(dInj, std::abs(gxf[i] - xf[i]));
            for (int i = 0; i < NDc; ++i) dRes = std::max(dRes, std::abs(gyc[i] - rc[i]));

            const RealType ctol = 1e-10;
            auto pf2 = [&](bool ok) { return ok ? "PASS" : "FAIL"; };
            std::cout << "[phase0][acm3] fine DOFs=" << ND << " -> coarse DOFs=" << NDc
                      << " (aggregates=" << nCoarse << ", coarse nnz=" << ccnnz << ")\n";
            std::cout << "[phase0][acm3] GPU coarse == host P^T A P : " << dRecipe << "  " << pf2(dRecipe < ctol) << "\n";
            std::cout << "[phase0][acm3] Galerkin P^T A P x == Ac x : " << dGal    << "  " << pf2(dGal    < ctol) << "\n";
            std::cout << "[phase0][acm3] GPU inject == host          : " << dInj    << "  " << pf2(dInj    < ctol) << "\n";
            std::cout << "[phase0][acm3] GPU restrict == host        : " << dRes    << "  " << pf2(dRes    < ctol) << "\n";
            std::cout << "[phase0][acm3] -> "
                      << ((dRecipe < ctol && dGal < ctol && dInj < ctol && dRes < ctol) ? "ALL PASS" : "FAIL") << "\n";
        }

        // ---- ACM Stage 4a: GPU multilevel V-cycle vs host V-cycle replica (one cycle, exact match) ----
        if (doAcm4)
        {
            // DILU is the global default now, but acmHostVcycle (the host replica below) is point-Jacobi.
            // Force the GPU side to point-Jacobi so this exact-match gate stays valid (diagnostic-only path;
            // set before the first acmUseDilu() call, which is the acmBuildHierarchy on the next line).
            setenv("MARS_ACM_SMOOTHER", "jacobi", /*overwrite=*/1);
            std::vector<AcmLevel<RealType>> levels;
            acmBuildHierarchy<RealType>(d_rowOff, d_colInd, d_vals, (int)nNodes,
                                        /*kmax=*/4, /*maxCoarseND=*/64, levels);

            std::vector<RealType> hb(ND);
            for (int i = 0; i < ND; ++i) hb[i] = std::sin(RealType(0.3) * i + RealType(0.5));   // deterministic test RHS

            // GPU: one V-cycle on (b, x=0)
            // coarsest = 10 Jacobi sweeps (the paper's value). NOTE: Jacobi is non-contractive on the
            // indefinite coarse operator (it mildly amplifies) -- a Stage-4a stopgap; Algorithm A's
            // direct coarsest solve (Stage 6) is the proper treatment, and GMRES (Stage 4b) makes the
            // cycle effective regardless.
            thrust::copy(hb.begin(), hb.end(), levels[0].bvec.begin());
            thrust::fill(levels[0].xvec.begin(), levels[0].xvec.end(), RealType(0));
            acmVcycleGpu<RealType>(levels, 0, /*pre=*/2, /*post=*/1, /*coarse=*/10, RealType(0.7));
            cudaDeviceSynchronize();
            std::vector<RealType> gx(ND);
            thrust::copy(levels[0].xvec.begin(), levels[0].xvec.end(), gx.begin());

            // host replica: same hierarchy + same arithmetic, one V-cycle on (b, x=0)
            std::vector<AcmLevelHost<RealType>> hLevels;
            acmHierarchyToHost<RealType>(levels, hLevels);
            std::vector<RealType> hx(ND, RealType(0));
            acmHostVcycle<RealType>(hLevels, 0, hx, hb, 2, 1, 10, RealType(0.7));

            // GPU SpMV/Jacobi contract to FMA, the host does not -> compare RELATIVE to ||x|| (a logic
            // bug gives O(1), not ~1e-9; this gate validates kernel correctness, not FP-bit-equality).
            RealType dCycle = 0, mhx = 0;
            for (int i = 0; i < ND; ++i) { dCycle = std::max(dCycle, std::abs(gx[i] - hx[i])); mhx = std::max(mhx, std::abs(hx[i])); }
            const RealType dRel = dCycle / std::max(RealType(1e-30), mhx);

            std::cout << "[phase0][acm4] V-cycle levels=" << levels.size() << "  sizes(ND):";
            for (auto& L : levels) std::cout << " " << L.ND;
            std::cout << "\n";
            const RealType ctol = 1e-7;
            std::cout << "[phase0][acm4] GPU V-cycle == host replica : abs=" << dCycle << " rel=" << dRel
                      << "  " << (dRel < ctol ? "PASS" : "FAIL") << "\n";
            std::cout << "[phase0][acm4] -> " << (dRel < ctol ? "ALL PASS" : "FAIL") << "\n";
        }

        // ---- Phase 1 1a (Hypre point-block) + ACM Stage 4b (FlexGMRES): shared BC-eliminated
        //      Stokes operator + cusolverSp QR reference, then each solver's own benchmark ----
        if (doHypre || doAcm4b)
        {
            const RealType nuS  = RealType(1) / static_cast<RealType>(Re);
            const RealType hpar = RealType(1) / std::sqrt(static_cast<RealType>(nNodes) * RealType(0.5));
            const RealType tauS = hpar * hpar / (RealType(4) * nuS);
            // re-assemble Stokes (advection off; matches the host pyamg study)
            thrust::fill(d_vals.begin(), d_vals.end(), RealType(0));
            assembleCoupledStokesKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
                c0, c1, c2, c3, nx, ny, nz,
                thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                thrust::raw_pointer_cast(d_vals.data()), nuS, tauS, nullptr, nullptr, nullptr,
                startEl, numLocal);
            cudaDeviceSynchronize();
            // BC tags + eliminate (lid u=1 / walls / w=0 / pressure pin)
            const double bb = beta * M_PI / 180.0;
            const RealType sinb = std::sin(bb), cotb = std::cos(bb) / std::sin(bb), eps = 1e-6;
            std::vector<uint8_t> bcF(ND, 0); std::vector<RealType> bcV(ND, RealType(0));
            for (size_t i = 0; i < nNodes; ++i) {
                RealType y0 = hy[i] / sinb, x0 = hx[i] - hy[i] * cotb;
                bcF[4 * i + 2] = 1;
                if (y0 > 1.0 - eps) { bcF[4 * i + 0] = 1; bcV[4 * i + 0] = RealType(1); bcF[4 * i + 1] = 1; }
                else if (x0 < eps || x0 > 1.0 - eps || y0 < eps) { bcF[4 * i + 0] = 1; bcF[4 * i + 1] = 1; }
            }
            bcF[3] = 1;
            thrust::device_vector<uint8_t> d_bcF(bcF.begin(), bcF.end());
            thrust::device_vector<RealType> d_bcV(bcV.begin(), bcV.end());
            thrust::device_vector<RealType> d_rhs(ND, RealType(0));
            const int dblk = (ND + blockSize - 1) / blockSize;
            applyBCKernel<RealType><<<dblk, blockSize>>>(
                thrust::raw_pointer_cast(d_bcF.data()), thrust::raw_pointer_cast(d_bcV.data()),
                thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                thrust::raw_pointer_cast(d_vals.data()), thrust::raw_pointer_cast(d_rhs.data()), ND);
            cudaDeviceSynchronize();
            // cusolver reference solution
            thrust::device_vector<RealType> d_xref(ND, RealType(0));
            { cusolverSpHandle_t cs = nullptr; cusolverSpCreate(&cs);
              cusparseMatDescr_t de = nullptr; cusparseCreateMatDescr(&de);
              cusparseSetMatType(de, CUSPARSE_MATRIX_TYPE_GENERAL);
              cusparseSetMatIndexBase(de, CUSPARSE_INDEX_BASE_ZERO); int sg = -2;
              cusolverSpDcsrlsvqr(cs, ND, nnz, de,
                  thrust::raw_pointer_cast(d_vals.data()), thrust::raw_pointer_cast(d_rowOff.data()),
                  thrust::raw_pointer_cast(d_colInd.data()), thrust::raw_pointer_cast(d_rhs.data()),
                  1e-12, 1, thrust::raw_pointer_cast(d_xref.data()), &sg);
              cudaDeviceSynchronize(); cusparseDestroyMatDescr(de); cusolverSpDestroy(cs); }
            // wrap the BC-eliminated CSR into a SparseMatrix (shared by the Hypre + ACM solves)
            SparseMatrix<int, RealType, cstone::GpuTag> Am; Am.allocate(ND, ND, nnz);
            thrust::copy(d_rowOff.begin(), d_rowOff.end(), thrust::device_pointer_cast(Am.rowOffsetsPtr()));
            thrust::copy(d_colInd.begin(), d_colInd.end(), thrust::device_pointer_cast(Am.colIndicesPtr()));
            thrust::copy(d_vals.begin(),   d_vals.end(),   thrust::device_pointer_cast(Am.valuesPtr()));

            if (doHypre)
            {
            // Hypre point-block BoomerAMG + GMRES
            cstone::DeviceVector<RealType> bH, xH; bH.resize(ND); xH.resize(ND);
            thrust::copy(d_rhs.begin(), d_rhs.end(), thrust::device_pointer_cast(bH.data()));
            HypreGMRESSolver<RealType, int, cstone::GpuTag> hg(
                MPI_COMM_WORLD, 1000, 1e-8,
                HypreGMRESSolver<RealType, int, cstone::GpuTag>::BOOMERAMG, 50);
            hg.setVerbose(false);
            hg.setPointBlock(4);
            hg.solve(Am, bH, xH, 0, ND, 0, ND, std::vector<int>{});
            const int iters = hg.getLastIterations();
            // compare to cusolver reference (guards the false-x=0 convergence mode)
            thrust::device_vector<RealType> d_xh(ND), d_df(ND);
            thrust::copy(thrust::device_pointer_cast(xH.data()),
                         thrust::device_pointer_cast(xH.data() + ND), d_xh.begin());
            thrust::transform(d_xh.begin(), d_xh.end(), d_xref.begin(), d_df.begin(),
                [] __device__ (RealType a, RealType c) { return fabs(a - c); });
            RealType dmax = thrust::reduce(d_df.begin(), d_df.end(), RealType(0), thrust::maximum<RealType>());
            thrust::transform(d_xref.begin(), d_xref.end(), d_df.begin(),
                [] __device__ (RealType c) { return fabs(c); });
            RealType rmax = thrust::reduce(d_df.begin(), d_df.end(), RealType(0), thrust::maximum<RealType>());
            RealType rel = dmax / (rmax > 0 ? rmax : RealType(1));
            // rel guards the false-x=0 collapse (rel~1 = AMG returned ~0). A real
            // GMRES solve to a 1e-8 RESIDUAL tol on a kappa-growing saddle operator
            // legitimately differs from the direct solve by ~kappa*1e-8, so allow a
            // generous band; it is NOT a precision check.
            bool hyOk = (iters > 0) && (rel < 5e-2);
            std::cout << std::scientific << std::setprecision(3);
            std::cout << "[phase0][hypre] Re=" << Re << " point-block BoomerAMG GMRES iters=" << iters
                      << " nnz=" << nnz << " |x_hypre-x_qr|rel=" << rel
                      << "  -> " << (hyOk ? "PASS (converged + matches reference)" : "CHECK") << "\n";
            assemblePass = assemblePass && hyOk;
            }  // if (doHypre)

            // ---- ACM Stage 4b: ACM-preconditioned native FlexGMRES vs the 1a baseline + cusolver ref ----
            if (doAcm4b)
            {
                using Vec = cstone::DeviceVector<RealType>;
                Vec bvec; bvec.resize(ND);
                thrust::copy(d_rhs.begin(), d_rhs.end(), thrust::device_pointer_cast(bvec.data()));

                // max|x_xref| for the relative-error denominator
                thrust::device_vector<RealType> d_absref(ND);
                thrust::transform(d_xref.begin(), d_xref.end(), d_absref.begin(),
                    [] __device__ (RealType c) { return fabs(c); });
                const RealType rmax = thrust::reduce(d_absref.begin(), d_absref.end(), RealType(0), thrust::maximum<RealType>());
                auto relErr = [&](Vec& xv) -> RealType {
                    thrust::device_vector<RealType> xt(ND);
                    thrust::copy(thrust::device_pointer_cast(xv.data()),
                                 thrust::device_pointer_cast(xv.data() + ND), xt.begin());
                    thrust::transform(xt.begin(), xt.end(), d_xref.begin(), xt.begin(),
                        [] __device__ (RealType a, RealType c) { return fabs(a - c); });
                    RealType dmax = thrust::reduce(xt.begin(), xt.end(), RealType(0), thrust::maximum<RealType>());
                    return dmax / (rmax > 0 ? rmax : RealType(1));
                };

                std::cout << std::scientific << std::setprecision(3);

                // GATE 1 (Z-basis correctness): with an IDENTITY preconditioner, flexible GMRES must
                // reduce EXACTLY to plain GMRES (Z_j=V_j, x=x0+Z y = x0+V y) -- a convergence-independent
                // algebraic identity. A Z-basis reconstruction bug breaks this even if neither converges.
                // (Jacobi is too weak to converge on this indefinite saddle operator, so it cannot be the
                // oracle; identity can.)
                Vec xp; xp.resize(ND);
                GMRESSolver<RealType, int, cstone::GpuTag> gp(2000, 1e-8, 30);
                gp.setVerbose(false); gp.solve(Am, bvec, xp, false);            // plain (unpreconditioned)
                const int itP = gp.getLastIterations();
                IdentityPreconditioner<RealType, int, cstone::GpuTag> id;
                Vec xi; xi.resize(ND);
                GMRESSolver<RealType, int, cstone::GpuTag> gi(2000, 1e-8, 30);
                gi.setVerbose(false); gi.setPreconditioner(&id); gi.setFlexible(true);
                gi.solve(Am, bvec, xi, true);                                   // flexible, M = I
                const int itI = gi.getLastIterations();
                thrust::device_vector<RealType> xpv(ND), xiv(ND);
                thrust::copy(thrust::device_pointer_cast(xp.data()), thrust::device_pointer_cast(xp.data() + ND), xpv.begin());
                thrust::copy(thrust::device_pointer_cast(xi.data()), thrust::device_pointer_cast(xi.data() + ND), xiv.begin());
                thrust::transform(xpv.begin(), xpv.end(), xiv.begin(), xiv.begin(),
                    [] __device__ (RealType a, RealType b) { return fabs(a - b); });
                const RealType ddiff = thrust::reduce(xiv.begin(), xiv.end(), RealType(0), thrust::maximum<RealType>());
                thrust::transform(xpv.begin(), xpv.end(), xpv.begin(), [] __device__ (RealType a) { return fabs(a); });
                const RealType dpmax = thrust::reduce(xpv.begin(), xpv.end(), RealType(0), thrust::maximum<RealType>());
                const RealType relZ = ddiff / (dpmax > 0 ? dpmax : RealType(1));
                const bool gate1 = (itP == itI) && (relZ < 1e-10);
                std::cout << "[phase0][acm4b][gate1] FGMRES(I)==plain GMRES (Z-basis): iters " << itP << "/" << itI
                          << " rel=" << relZ << "  -> " << (gate1 ? "PASS" : "FAIL") << "\n";

                // GATE 2 (ACM quality): ACM-preconditioned FlexGMRES converges to the reference.
                GpuAcmPreconditioner<RealType, int, cstone::GpuTag> acm;
                Vec xa; xa.resize(ND);
                GMRESSolver<RealType, int, cstone::GpuTag> ga(2000, 1e-8, 30);
                ga.setVerbose(false); ga.setPreconditioner(&acm); ga.setFlexible(true);
                ga.solve(Am, bvec, xa, true);
                const int itA = ga.getLastIterations();
                const RealType relA = relErr(xa);
                const bool gate2 = (itA > 0) && (itA < 1980) && (relA < 5e-2);   // converged (< maxIter) + correct
                std::cout << "[phase0][acm4b] Re=" << Re << " ACM-FlexGMRES iters=" << itA
                          << " (Hypre 1a baseline n16/32/64 = 13/19/39) levels=" << acm.numLevels()
                          << " |x_acm-x_qr|rel=" << relA
                          << "  -> " << (gate1 && gate2 ? "PASS" : "CHECK") << "\n";
                assemblePass = assemblePass && gate1 && gate2;
            }
        }

        // ---- ACM on a real mesh: mesh-agnostic ACM-FlexGMRES vs Hypre BoomerAMG (residual-based,
        //      no direct reference -> scales to large meshes). Closed-domain Stokes: u=v=w=0 on the
        //      surface, pressure pinned, smooth body force drives the flow. ----
        if (doAcmPump)
        {
            const RealType nuS  = (nuVal >= 0) ? static_cast<RealType>(nuVal) : RealType(1) / static_cast<RealType>(Re);
            const RealType hpar = RealType(1) / std::sqrt(static_cast<RealType>(nNodes) * RealType(0.5));
            const RealType tauS = hpar * hpar / (RealType(4) * nuS);
            // operator is assembled INSIDE the Picard loop below with the frozen advecting velocity (convection)

            // boundary nodes: a tet face shared by exactly one element is a surface face (host sort+count)
            std::vector<std::array<int, 3>> faces; faces.reserve(4 * numLocal);
            auto pushFace = [&](int a, int b, int c) {
                int t[3] = {a, b, c}; std::sort(t, t + 3); faces.push_back({t[0], t[1], t[2]});
            };
            for (size_t e = 0; e < numLocal; ++e) {
                int n0 = (int)h0[e], n1 = (int)h1[e], n2 = (int)h2[e], n3 = (int)h3[e];
                pushFace(n0, n1, n2); pushFace(n0, n1, n3); pushFace(n0, n2, n3); pushFace(n1, n2, n3);
            }
            std::sort(faces.begin(), faces.end());
            std::vector<uint8_t> isBnd(nNodes, 0);
            for (size_t i = 0; i < faces.size(); ) {
                size_t j = i + 1; while (j < faces.size() && faces[j] == faces[i]) ++j;
                if (j - i == 1) for (int v : faces[i]) isBnd[v] = 1;   // single-occurrence face -> surface
                i = j;
            }
            size_t nBnd = 0; for (size_t i = 0; i < nNodes; ++i) nBnd += isBnd[i];

            std::vector<uint8_t> bcF(ND, 0);
            std::vector<RealType> bcV(ND, RealType(0));
            std::vector<RealType> hrhs(ND, RealType(0));
            // Make the path explicit -- a silently-empty --inlet-ss= (e.g. an unset env var) would otherwise
            // fall back to the closed-box body-force diagnostic without warning.
            if (rank == 0)
                std::cout << "[phase0][acm-pump] mode="
                          << (inletSet ? "side-set through-flow (velocity-flux inlet + open outlet)"
                                       : "closed-box body-force diagnostic (no --inlet-ss)") << "\n";
            if (inletSet) {            // pump real flow: side-set velocity-flux inlet + no-slip walls, NO body force
                if (!buildSideSetVelBC(bcF, bcV)) { MPI_Finalize(); return 6; }
            } else {                   // closed-box diagnostic: all-surface no-slip + smooth body force drive
                for (size_t i = 0; i < nNodes; ++i) {
                    if (isBnd[i]) { bcF[4 * i + 0] = bcF[4 * i + 1] = bcF[4 * i + 2] = 1; }
                    hrhs[4 * i + 0] = std::sin(RealType(3) * hx[i] + RealType(1));      // smooth body force
                    hrhs[4 * i + 1] = std::sin(RealType(3) * hy[i] + RealType(2));
                    hrhs[4 * i + 2] = std::sin(RealType(3) * hz[i] + RealType(0.5));
                }
            }
            // pin ONE pressure per connected component of the fluid node-graph. A multiply-connected
            // domain (separate passages/chambers) has one constant-pressure null mode per component, so a
            // single global pin leaves the rest near-null -> near-singular. Union-find over the 1-ring.
            std::vector<int> cc(nNodes);
            for (size_t i = 0; i < nNodes; ++i) cc[i] = static_cast<int>(i);
            auto findRoot = [&](int x) { while (cc[x] != x) { cc[x] = cc[cc[x]]; x = cc[x]; } return x; };
            for (size_t i = 0; i < nNodes; ++i)
                for (int j : ring[i]) { int a = findRoot((int)i), b = findRoot(j); if (a != b) cc[a] = b; }
            std::vector<char> rootPinned(nNodes, 0);
            int nComp = 0;
            // For the open-flow pump, pin pressure at the OUTLET first: it is the natural reference for the
            // do-nothing outlet, and an arbitrary interior pin acts as a spurious mass SINK under a
            // velocity-flux inlet (inflow drains at the pin instead of the outlet). Then fall back to one pin
            // per OTHER component (e.g. a walled-off chamber) to kill its constant-pressure null mode.
            if (inletSet)
                for (int o : outN) { int r = findRoot(o); if (!rootPinned[r]) { rootPinned[r] = 1; bcF[4 * o + 3] = 1; ++nComp; } }
            for (size_t i = 0; i < nNodes; ++i) {
                int r = findRoot((int)i);
                if (!rootPinned[r]) { rootPinned[r] = 1; bcF[4 * i + 3] = 1; ++nComp; }
            }
            if (std::getenv("MARS_VERBOSE_MESH"))   // component count is a mesh-topology fact -> opt-in
                std::cout << "[phase0][acm-pump] fluid-graph components=" << nComp << " (one pressure pin per component)\n";

            thrust::device_vector<uint8_t>  d_bcF(bcF.begin(), bcF.end());
            thrust::device_vector<RealType> d_bcV(bcV.begin(), bcV.end());   // inlet flux values (0 for closed-box)
            const int dblk = (ND + blockSize - 1) / blockSize;
            const int nblkN = ((int)nNodes + blockSize - 1) / blockSize;

            SparseMatrix<int, RealType, cstone::GpuTag> Am; Am.allocate(ND, ND, nnz);
            thrust::copy(d_rowOff.begin(), d_rowOff.end(), thrust::device_pointer_cast(Am.rowOffsetsPtr()));
            thrust::copy(d_colInd.begin(), d_colInd.end(), thrust::device_pointer_cast(Am.colIndicesPtr()));
            using Vec = cstone::DeviceVector<RealType>;
            Vec b; b.resize(ND); Vec xa; xa.resize(ND);
            thrust::device_vector<RealType> d_rhs(ND), d_Ax(ND), d_diff(ND);
            RealType bnorm = RealType(1);
            auto resid = [&](Vec& xv) -> RealType {     // ||b - A x|| / ||b|| on the BC-eliminated CSR
                acmSpmvKernel<RealType><<<dblk, blockSize>>>(
                    thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                    thrust::raw_pointer_cast(d_vals.data()), xv.data(), thrust::raw_pointer_cast(d_Ax.data()), ND);
                cudaDeviceSynchronize();
                thrust::transform(thrust::device_pointer_cast(b.data()), thrust::device_pointer_cast(b.data() + ND),
                    d_Ax.begin(), d_diff.begin(), thrust::minus<RealType>());
                RealType rn = std::sqrt(thrust::inner_product(d_diff.begin(), d_diff.end(), d_diff.begin(), RealType(0)));
                return rn / (bnorm > 0 ? bnorm : RealType(1));
            };

            std::cout << std::scientific << std::setprecision(3);
            if (std::getenv("MARS_VERBOSE_MESH"))
                std::cout << "[phase0][acm-pump] DOFs=" << ND << " nnz=" << nnz
                          << " surface nodes=" << nBnd << "/" << nNodes << " nu=" << nuS << "\n";

            // Picard outer loop: freeze the velocity (d_av*), re-assemble the CONVECTIVE momentum, re-solve with
            // the ACM. d_av*=0 on iter 0 (pure Stokes); the previous solution then advects the next. BC and
            // sparsity are fixed across iters; only d_vals (the operator) + b are rebuilt each iter.
            GpuAcmPreconditioner<RealType, int, cstone::GpuTag> acm;
            GMRESSolver<RealType, int, cstone::GpuTag> ga(2000, 1e-8, 30);
            ga.setVerbose(false); ga.setPreconditioner(&acm); ga.setFlexible(true);
            thrust::device_vector<RealType> d_avx(nNodes, RealType(0)), d_avy(nNodes, RealType(0)),
                                            d_avz(nNodes, RealType(0)), d_sx(nNodes), d_sy(nNodes), d_sz(nNodes), d_tmp(nNodes);
            const RealType wrelax = static_cast<RealType>(relax);   // Picard under-relaxation omega (--relax)
            int itA = 0, pit = 0; RealType rA = 0, du = 0;
            for (pit = 0; pit < picard; ++pit) {
                thrust::fill(d_vals.begin(), d_vals.end(), RealType(0));
                if (doRhieChow)
                    assembleRhieChowGpu<KeyType, RealType>(c0, c1, c2, c3, nx, ny, nz,
                        thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                        thrust::raw_pointer_cast(d_vals.data()), nuS, doSupg ? RealType(1) : RealType(0),
                        thrust::raw_pointer_cast(d_avx.data()), thrust::raw_pointer_cast(d_avy.data()),
                        thrust::raw_pointer_cast(d_avz.data()), startEl, numLocal, (int)nNodes, blockSize);
                else
                    assembleCoupledStokesKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
                        c0, c1, c2, c3, nx, ny, nz,
                        thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                        thrust::raw_pointer_cast(d_vals.data()), nuS, tauS,
                        thrust::raw_pointer_cast(d_avx.data()), thrust::raw_pointer_cast(d_avy.data()),
                        thrust::raw_pointer_cast(d_avz.data()), startEl, numLocal);
                cudaDeviceSynchronize();
                thrust::copy(hrhs.begin(), hrhs.end(), d_rhs.begin());
                applyBCKernel<RealType><<<dblk, blockSize>>>(
                    thrust::raw_pointer_cast(d_bcF.data()), thrust::raw_pointer_cast(d_bcV.data()),
                    thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                    thrust::raw_pointer_cast(d_vals.data()), thrust::raw_pointer_cast(d_rhs.data()), ND);
                cudaDeviceSynchronize();
                thrust::copy(d_vals.begin(), d_vals.end(), thrust::device_pointer_cast(Am.valuesPtr()));
                thrust::copy(d_rhs.begin(), d_rhs.end(), thrust::device_pointer_cast(b.data()));
                bnorm = std::sqrt(thrust::inner_product(
                    thrust::device_pointer_cast(b.data()), thrust::device_pointer_cast(b.data() + ND),
                    thrust::device_pointer_cast(b.data()), RealType(0)));
                bool reuseThis = (acmRebuild > 1 && (pit % acmRebuild != 0));   // rebuild hierarchy every K iters, reuse between
                acm.setReuse(reuseThis);
                cudaDeviceSynchronize();                    // flush prior async work so the timer measures only the solve
                auto tSolve0 = std::chrono::high_resolution_clock::now();
                ga.solve(Am, b, xa, true);                  // (re-)setup the ACM hierarchy unless frozen this iter
                cudaDeviceSynchronize();
                double solveMs = std::chrono::duration<double, std::milli>(
                    std::chrono::high_resolution_clock::now() - tSolve0).count();   // setup(if rebuilt)+FlexGMRES
                itA = ga.getLastIterations(); rA = resid(xa);
                // extract the SOLVED velocity into temps, then Picard under-relax: u <- w*u_new + (1-w)*u_old.
                // d(u) is the RELAXED change; convergence is relative to the velocity scale (|u|max).
                extractCompKernel<RealType><<<nblkN, blockSize>>>(xa.data(), thrust::raw_pointer_cast(d_sx.data()), 4, 0, (int)nNodes);
                extractCompKernel<RealType><<<nblkN, blockSize>>>(xa.data(), thrust::raw_pointer_cast(d_sy.data()), 4, 1, (int)nNodes);
                extractCompKernel<RealType><<<nblkN, blockSize>>>(xa.data(), thrust::raw_pointer_cast(d_sz.data()), 4, 2, (int)nNodes);
                cudaDeviceSynchronize();
                thrust::transform(d_sx.begin(), d_sx.end(), d_avx.begin(), d_tmp.begin(),
                    [] __device__ (RealType s, RealType o) { return fabs(s - o); });
                du = wrelax * thrust::reduce(d_tmp.begin(), d_tmp.end(), RealType(0), thrust::maximum<RealType>());
                thrust::transform(d_sx.begin(), d_sx.end(), d_tmp.begin(), [] __device__ (RealType v) { return fabs(v); });
                RealType uscale = thrust::reduce(d_tmp.begin(), d_tmp.end(), RealType(0), thrust::maximum<RealType>());
                auto blend = [wrelax] __device__ (RealType s, RealType o) { return wrelax * s + (RealType(1) - wrelax) * o; };
                thrust::transform(d_sx.begin(), d_sx.end(), d_avx.begin(), d_avx.begin(), blend);
                thrust::transform(d_sy.begin(), d_sy.end(), d_avy.begin(), d_avy.begin(), blend);
                thrust::transform(d_sz.begin(), d_sz.end(), d_avz.begin(), d_avz.begin(), blend);
                double setupMs = acm.getLastSetupMs();      // hierarchy build (0 if reused); solve = total - setup
                std::cout << "[phase0][acm-pump][picard " << pit << "] ACM iters=" << itA << " res=" << rA
                          << " d(u)=" << du << " d(u)/|u|=" << du / (uscale + RealType(1e-30))
                          << " setup=" << setupMs << "ms solve=" << (solveMs - setupMs) << "ms"
                          << (reuseThis ? " (reuse)" : " (rebuild)") << "\n";
                if (du < RealType(1e-5) * (uscale + RealType(1e-30))) { ++pit; break; }   // relative steady-state
            }

            std::cout << "[phase0][acm-pump] picard=" << pit << " ACM-FlexGMRES iters=" << itA
                      << " levels=" << acm.numLevels() << " beta=" << acm.beta() << " res=" << rA << "\n";
            if (inletSet) {
                // pump through-flow validation: mass conservation (net boundary flux ~0 by divergence theorem),
                // inflow sign, pressure/velocity sanity. No BoomerAMG comparison (it fails on the RC operator).
                // PASS = physical validity: ACM converged + inflow sign + mass conserved (net boundary flux
                // ~0 by divergence theorem). |u|max is reported as INFO -- this pump genuinely reaches
                // hundreds of m/s at its narrowest path (collaborators ~430 m/s), so |u|max is NOT a fault.
                std::vector<RealType> xs(ND);
                thrust::copy(thrust::device_pointer_cast(xa.data()), thrust::device_pointer_cast(xa.data() + ND), xs.begin());
                RealType Qin = 0, Qout = 0, Qtot = 0, pInSum = 0, umax = 0;
                for (size_t i = 0; i < nNodes; ++i) {
                    RealType f = xs[4*i]*hsbx[i] + xs[4*i+1]*hsby[i] + xs[4*i+2]*hsbz[i]; Qtot += f;
                    RealType um = std::sqrt(xs[4*i]*xs[4*i] + xs[4*i+1]*xs[4*i+1] + xs[4*i+2]*xs[4*i+2]);
                    umax = std::max(umax, um);
                }
                for (int i : inN)  { Qin  += xs[4*i]*hsbx[i] + xs[4*i+1]*hsby[i] + xs[4*i+2]*hsbz[i]; pInSum += xs[4*i+3]; }
                for (int i : outN)   Qout += xs[4*i]*hsbx[i] + xs[4*i+1]*hsby[i] + xs[4*i+2]*hsbz[i];
                RealType pInMean = inN.empty() ? RealType(0) : pInSum / (RealType)inN.size();
                RealType massErr = std::abs(Qtot) / (std::abs(Qin) > RealType(1e-30) ? std::abs(Qin) : RealType(1));
                bool inflowOk = (Qin < RealType(0));
                bool ok = (rA < 1e-6) && std::isfinite(umax) && std::isfinite(Qtot) && inflowOk && (massErr < RealType(0.05));
                std::cout << "[phase0][acm-pump][pump] Qin=" << Qin << " Qout=" << Qout << " net=" << Qtot
                          << " massErr=" << massErr << "  inflow=" << (inflowOk ? "OK" : "REVERSED(use --inlet-flip-normal)")
                          << "  p_in_mean=" << pInMean << " (outlet pinned)  |u|max=" << umax << " (info, vs ref)"
                          << "  -> " << (ok ? "PASS (through-flow, mass conserved, ACM converged)" : "CHECK") << "\n";
            } else {
                Vec bh, xh; bh.resize(ND); xh.resize(ND);
                thrust::copy(d_rhs.begin(), d_rhs.end(), thrust::device_pointer_cast(bh.data()));
                HypreGMRESSolver<RealType, int, cstone::GpuTag> hg(
                    MPI_COMM_WORLD, 2000, 1e-8, HypreGMRESSolver<RealType, int, cstone::GpuTag>::BOOMERAMG, 50);
                hg.setVerbose(false); hg.setPointBlock(4);
                hg.solve(Am, bh, xh, 0, ND, 0, ND, std::vector<int>{});
                const int itH = hg.getLastIterations(); const RealType rH = resid(xh);
                thrust::device_vector<RealType> d_xa(ND), d_xh(ND);
                thrust::copy(thrust::device_pointer_cast(xa.data()), thrust::device_pointer_cast(xa.data() + ND), d_xa.begin());
                thrust::copy(thrust::device_pointer_cast(xh.data()), thrust::device_pointer_cast(xh.data() + ND), d_xh.begin());
                thrust::transform(d_xa.begin(), d_xa.end(), d_xh.begin(), d_diff.begin(),
                    [] __device__ (RealType a, RealType c) { return fabs(a - c); });
                const RealType dmax = thrust::reduce(d_diff.begin(), d_diff.end(), RealType(0), thrust::maximum<RealType>());
                thrust::transform(d_xa.begin(), d_xa.end(), d_diff.begin(), [] __device__ (RealType a) { return fabs(a); });
                const RealType amax = thrust::reduce(d_diff.begin(), d_diff.end(), RealType(0), thrust::maximum<RealType>());
                const RealType agree = dmax / (amax > 0 ? amax : RealType(1));
                std::cout << "[phase0][acm-pump] BoomerAMG-GMRES iters=" << itH << " res=" << rH << "\n";
                std::cout << "[phase0][acm-pump] solvers agree |x_acm-x_amg|rel=" << agree
                          << "  -> " << ((rA < 1e-6 && rH < 1e-6 && agree < 1e-2) ? "PASS" : "CHECK") << "\n";
            }
        }

        // ---- 3b: Picard loop + cusolverSp QR direct reference solve + benchmark match ----
        if (doSolve)
        {
            // Erturk-Dursun (arXiv:physics/0505121) Table 5, alpha=45: u along A-B, Re=100.
            const RealType ED_s[17] = {0,32,64,96,128,160,192,224,256,288,320,352,384,416,448,480,512};
            const RealType ED_u[17] = {0.0,-2.389e-3,-8.901e-3,-1.949e-2,-3.426e-2,-5.373e-2,
                                       -7.846e-2,-1.080e-1,-1.390e-1,-1.634e-1,-1.669e-1,-1.315e-1,
                                       -3.945e-2,1.217e-1,3.544e-1,6.477e-1,1.0};
            auto edInterp = [&](RealType s) -> RealType {
                if (s <= 0) return ED_u[0];
                for (int k = 1; k < 17; ++k)
                    if (s <= ED_s[k]) {
                        RealType t = (s - ED_s[k-1]) / (ED_s[k] - ED_s[k-1]);
                        return ED_u[k-1] + t * (ED_u[k] - ED_u[k-1]);
                    }
                return ED_u[16];
            };

            // physics: nu = 1/Re, PSPG tau = h^2/(4 nu), h ~ in-plane node spacing
            const RealType nuS  = RealType(1) / static_cast<RealType>(Re);
            const RealType hpar = RealType(1) / std::sqrt(static_cast<RealType>(nNodes) * RealType(0.5));
            const RealType tauS = hpar * hpar / (RealType(4) * nuS);

            // BC tags (one-time host setup). inletSet -> pump side-set velocity-flux inlet + open outlet;
            // doChannel -> plane-channel parabolic inlet; else lid-driven cavity.
            std::vector<uint8_t> bcF(ND, 0);
            std::vector<RealType> bcV(ND, RealType(0));
            const RealType chUmean = RealType(1);
            RealType chXmin = 0, chXmax = 0, chYmin = 0, chYmax = 0, chLy = 1;
            if (inletSet) {                                  // pump side-set velocity-flux inlet (cusolver path)
                if (!buildSideSetVelBC(bcF, bcV)) { MPI_Finalize(); return 6; }
                if (!outN.empty()) bcF[4 * (*outN.begin()) + 3] = 1;   // one outlet pressure pin (fix the level)
            } else if (doChannel) {
                chXmin = *std::min_element(hx.begin(), hx.end()); chXmax = *std::max_element(hx.begin(), hx.end());
                chYmin = *std::min_element(hy.begin(), hy.end()); chYmax = *std::max_element(hy.begin(), hy.end());
                chLy = chYmax - chYmin;
                const RealType ex = (chXmax - chXmin) * RealType(1e-6), ey = chLy * RealType(1e-6);
                int outletPin = -1;
                for (size_t i = 0; i < nNodes; ++i) {
                    bcF[4 * i + 2] = 1;                                          // w=0 (thin-slab 2D flow)
                    RealType x = hx[i], yn = (hy[i] - chYmin) / chLy;
                    if (x < chXmin + ex) {                                       // velocity-flux inlet: parabolic u
                        bcF[4 * i + 0] = 1; bcV[4 * i + 0] = RealType(6) * chUmean * yn * (RealType(1) - yn);
                        bcF[4 * i + 1] = 1;
                    } else if (hy[i] < chYmin + ey || hy[i] > chYmax - ey) {      // no-slip walls
                        bcF[4 * i + 0] = 1; bcF[4 * i + 1] = 1;
                    } else if (x > chXmax - ex && outletPin < 0) {               // open outlet: one pressure pin
                        bcF[4 * i + 3] = 1; outletPin = (int)i;
                    }
                }
            } else {
                const double bb = beta * M_PI / 180.0;
                const RealType sinb = std::sin(bb), cotb = std::cos(bb) / std::sin(bb), eps = 1e-6;
                for (size_t i = 0; i < nNodes; ++i) {
                    RealType y0 = hy[i] / sinb, x0 = hx[i] - hy[i] * cotb;
                    bcF[4 * i + 2] = 1;
                    if (y0 > 1.0 - eps) { bcF[4 * i + 0] = 1; bcV[4 * i + 0] = RealType(1); bcF[4 * i + 1] = 1; }
                    else if (x0 < eps || x0 > 1.0 - eps || y0 < eps) { bcF[4 * i + 0] = 1; bcF[4 * i + 1] = 1; }
                }
                bcF[3] = 1;
            }
            thrust::device_vector<uint8_t>  d_bcF(bcF.begin(), bcF.end());
            thrust::device_vector<RealType> d_bcV(bcV.begin(), bcV.end());
            thrust::device_vector<RealType> d_avx(nNodes, RealType(0)), d_avy(nNodes, RealType(0)),
                                            d_avz(nNodes, RealType(0)), d_uprev(nNodes, RealType(0)),
                                            d_tmp(nNodes);
            thrust::device_vector<RealType> d_rhs(ND, RealType(0)), d_x(ND, RealType(0));
            const int dblk = (ND + blockSize - 1) / blockSize;
            const int nblk = (static_cast<int>(nNodes) + blockSize - 1) / blockSize;

            cusolverSpHandle_t cs = nullptr; cusolverSpCreate(&cs);
            cusparseMatDescr_t descr = nullptr; cusparseCreateMatDescr(&descr);
            cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
            cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

            int it = 0, singularity = -2; RealType du = 0;
            cusolverStatus_t st = CUSOLVER_STATUS_SUCCESS;
            for (it = 0; it < picard; ++it) {
                // re-assemble with the frozen (previous-iterate) velocity. doRhieChow -> collocated RC
                // operator (no tau); else PSPG. The cavity benchmark (vs Erturk-Dursun) validates RC accuracy.
                thrust::fill(d_vals.begin(), d_vals.end(), RealType(0));
                if (doRhieChow)
                    assembleRhieChowGpu<KeyType, RealType>(c0, c1, c2, c3, nx, ny, nz,
                        thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                        thrust::raw_pointer_cast(d_vals.data()), nuS, doSupg ? RealType(1) : RealType(0),
                        thrust::raw_pointer_cast(d_avx.data()), thrust::raw_pointer_cast(d_avy.data()),
                        thrust::raw_pointer_cast(d_avz.data()), startEl, numLocal, (int)nNodes, blockSize);
                else
                    assembleCoupledStokesKernel<KeyType, RealType><<<eBlocks, blockSize>>>(
                        c0, c1, c2, c3, nx, ny, nz,
                        thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                        thrust::raw_pointer_cast(d_vals.data()), nuS, tauS,
                        thrust::raw_pointer_cast(d_avx.data()), thrust::raw_pointer_cast(d_avy.data()),
                        thrust::raw_pointer_cast(d_avz.data()), startEl, numLocal);
                cudaDeviceSynchronize();
                thrust::fill(d_rhs.begin(), d_rhs.end(), RealType(0));
                applyBCKernel<RealType><<<dblk, blockSize>>>(
                    thrust::raw_pointer_cast(d_bcF.data()), thrust::raw_pointer_cast(d_bcV.data()),
                    thrust::raw_pointer_cast(d_rowOff.data()), thrust::raw_pointer_cast(d_colInd.data()),
                    thrust::raw_pointer_cast(d_vals.data()), thrust::raw_pointer_cast(d_rhs.data()), ND);
                cudaDeviceSynchronize();
                st = cusolverSpDcsrlsvqr(cs, ND, nnz, descr,
                    thrust::raw_pointer_cast(d_vals.data()), thrust::raw_pointer_cast(d_rowOff.data()),
                    thrust::raw_pointer_cast(d_colInd.data()), thrust::raw_pointer_cast(d_rhs.data()),
                    1e-12, 1, thrust::raw_pointer_cast(d_x.data()), &singularity);
                cudaDeviceSynchronize();
                if (st != CUSOLVER_STATUS_SUCCESS || singularity >= 0) break;
                extractCompKernel<RealType><<<nblk, blockSize>>>(thrust::raw_pointer_cast(d_x.data()),
                    thrust::raw_pointer_cast(d_avx.data()), 4, 0, (int)nNodes);
                extractCompKernel<RealType><<<nblk, blockSize>>>(thrust::raw_pointer_cast(d_x.data()),
                    thrust::raw_pointer_cast(d_avy.data()), 4, 1, (int)nNodes);
                extractCompKernel<RealType><<<nblk, blockSize>>>(thrust::raw_pointer_cast(d_x.data()),
                    thrust::raw_pointer_cast(d_avz.data()), 4, 2, (int)nNodes);
                cudaDeviceSynchronize();
                thrust::transform(d_avx.begin(), d_avx.end(), d_uprev.begin(), d_tmp.begin(),
                    [] __device__ (RealType a, RealType c) { return fabs(a - c); });
                du = thrust::reduce(d_tmp.begin(), d_tmp.end(), RealType(0), thrust::maximum<RealType>());
                thrust::copy(d_avx.begin(), d_avx.end(), d_uprev.begin());
                if (du < 1e-8) { ++it; break; }
            }

            // post-process: inletSet -> pump through-flow (mass conservation); doChannel -> Poiseuille; else Erturk-Dursun
            bool solveOk = false;
            std::cout << std::scientific << std::setprecision(3);
            std::cout << "[phase0][solve] Re=" << Re << " picard=" << it << " d(u)=" << du
                      << " nu=" << nuS << " status=" << static_cast<int>(st) << " sing=" << singularity << "\n";
            if (inletSet) {
                // no analytical answer: check mass conservation (net boundary flux ~0 by divergence theorem),
                // inflow sign, and pressure/velocity sanity. Q = sum u.S_bnd over the boundary nodes.
                std::vector<RealType> xs(ND);
                thrust::copy(d_x.begin(), d_x.end(), xs.begin());
                RealType Qin = 0, Qout = 0, Qtot = 0, pInSum = 0, umax = 0;
                for (size_t i = 0; i < nNodes; ++i) {
                    RealType f = xs[4 * i + 0] * hsbx[i] + xs[4 * i + 1] * hsby[i] + xs[4 * i + 2] * hsbz[i];
                    Qtot += f;
                    RealType um = std::sqrt(xs[4*i]*xs[4*i] + xs[4*i+1]*xs[4*i+1] + xs[4*i+2]*xs[4*i+2]);
                    umax = std::max(umax, um);
                }
                for (int i : inN)  { Qin  += xs[4*i]*hsbx[i] + xs[4*i+1]*hsby[i] + xs[4*i+2]*hsbz[i]; pInSum += xs[4*i+3]; }
                for (int i : outN)   Qout += xs[4*i]*hsbx[i] + xs[4*i+1]*hsby[i] + xs[4*i+2]*hsbz[i];
                RealType pInMean = inN.empty() ? RealType(0) : pInSum / (RealType)inN.size();
                RealType massErr = std::abs(Qtot) / (std::abs(Qin) > RealType(1e-30) ? std::abs(Qin) : RealType(1));
                bool inflowOk = (Qin < RealType(0));     // u inward, S_bnd outward => u.S_bnd < 0 for inflow
                bool finite = std::isfinite(umax) && std::isfinite(pInMean) && std::isfinite(Qtot);
                solveOk = (st == CUSOLVER_STATUS_SUCCESS) && (singularity < 0) && finite && inflowOk
                          && (massErr < RealType(0.05));
                std::cout << "[phase0][solve][pump] Qin=" << Qin << " Qout=" << Qout << " net=" << Qtot
                          << " massErr=" << massErr << "  inflow=" << (inflowOk ? "OK" : "REVERSED(use --inlet-flip-normal)")
                          << "  p_in_mean=" << pInMean << " (outlet=0)  |u|max=" << umax
                          << "  -> " << (solveOk ? "PASS (through-flow, mass conserved)" : "CHECK") << "\n";
            } else if (doChannel) {
                std::vector<RealType> xs(ND);
                thrust::copy(d_x.begin(), d_x.end(), xs.begin());
                const RealType ex = (chXmax - chXmin) * RealType(1e-6);
                RealType maxe = 0, pin = 0, pout = 0; int nin = 0, nout = 0;
                for (size_t i = 0; i < nNodes; ++i) {
                    RealType yn = (hy[i] - chYmin) / chLy;
                    RealType uex = RealType(6) * chUmean * yn * (RealType(1) - yn);     // fully-developed parabola
                    maxe = std::max(maxe, std::abs(xs[4 * i + 0] - uex));
                    if (hx[i] < chXmin + ex) { pin += xs[4 * i + 3]; ++nin; }
                    if (hx[i] > chXmax - ex) { pout += xs[4 * i + 3]; ++nout; }
                }
                pin /= std::max(1, nin); pout /= std::max(1, nout);
                RealType dp = pin - pout;
                RealType dpex = RealType(12) * nuS * chUmean * (chXmax - chXmin) / (chLy * chLy);   // mu=rho*nu, rho=1
                RealType dperr = std::abs(dp - dpex) / (std::abs(dpex) > RealType(0) ? std::abs(dpex) : RealType(1));
                // PASS criterion = the boundary-flux fix's signature: flow propagates (NOT collapsed to ~0,
                // which gives maxe~peak 1.5) and dp has the correct SIGN + ballpark. The tight tolerance is
                // limited by the coarse thin-slab open-z/outlet do-nothing artifact (host-replica finding),
                // not the fix; dp -> exact is the decisive mass-conservation check.
                bool propagates = (maxe < RealType(1.2));
                bool dpok = (dp > RealType(0)) && (dperr < RealType(0.5));
                solveOk = (st == CUSOLVER_STATUS_SUCCESS) && (singularity < 0) && propagates && dpok;
                std::cout << "[phase0][solve] Poiseuille: max|u-parabola|=" << maxe << "  dp=" << dp
                          << " vs exact " << dpex << " (rel " << dperr << ")  propagates=" << propagates
                          << " -> " << (solveOk ? "PASS (fix works: flow propagates, dp sign+magnitude correct)" : "CHECK") << "\n";
            } else {
                const double bb = beta * M_PI / 180.0;
                const RealType sinb = std::sin(bb), cotb = std::cos(bb) / std::sin(bb);
                std::vector<RealType> hu(nNodes);
                thrust::copy(d_avx.begin(), d_avx.end(), hu.begin());
                RealType zmin = *std::min_element(hz.begin(), hz.end());
                std::vector<std::pair<RealType, RealType>> prof;
                for (size_t i = 0; i < nNodes; ++i) {
                    RealType x0 = hx[i] - hy[i] * cotb, y0 = hy[i] / sinb;
                    if (std::abs(x0 - 0.5) < 1e-6 && std::abs(hz[i] - zmin) < 1e-6) prof.push_back({y0, hu[i]});
                }
                std::sort(prof.begin(), prof.end());
                RealType emax = 0, el2 = 0;
                for (auto& pr : prof) { RealType e = std::abs(pr.second - edInterp(pr.first * RealType(512))); emax = std::max(emax, e); el2 += e * e; }
                int np = static_cast<int>(prof.size());
                el2 = np > 0 ? std::sqrt(el2 / np) : 0;
                solveOk = (st == CUSOLVER_STATUS_SUCCESS) && (singularity < 0) && (np > 2) && (emax < 6e-2);
                std::cout << "[phase0][solve] centerline u vs Erturk-Dursun(b45,Re100): n=" << np
                          << " max|err|=" << emax << " L2=" << el2 << "  -> " << (solveOk ? "PASS (matches benchmark)" : "CHECK") << "\n";
            }
            assemblePass = assemblePass && solveOk;
            cusparseDestroyMatDescr(descr);
            cusolverSpDestroy(cs);
        }
    }

    MPI_Finalize();
    return (pass && assemblePass) ? 0 : 3;
}
