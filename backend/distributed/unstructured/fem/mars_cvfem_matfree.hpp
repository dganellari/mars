#pragma once

#include "mars_cvfem_hex_kernel.hpp"

namespace mars {
namespace fem {

// Per-element 8x8 CVFEM diffusion stiffness (the SCS-flux "lhs" the tensor
// assembler builds), shared by the recompute-apply and the store-operator build
// so all paths form the SAME operator. Row-major lhs[i*8+j]. Diffusion only.
template<typename KeyType, typename RealType>
__device__ inline void cvfem_hex_diffusion_lhs(
    const RealType coords[8][3], const RealType gamma[8],
    const RealType* __restrict__ avx, const RealType* __restrict__ avy, const RealType* __restrict__ avz,
    size_t e, RealType lhs[64])
{
    #pragma unroll
    for (int i = 0; i < 64; ++i) lhs[i] = 0.0;
    #pragma unroll
    for (int ip = 0; ip < 12; ++ip) {
        int nodeL = hexLRSCV[ip*2], nodeR = hexLRSCV[ip*2+1];
        RealType areaVec[3] = {avx[e*12+ip], avy[e*12+ip], avz[e*12+ip]};
        RealType gamma_ip = 0.0;
        #pragma unroll
        for (int n = 0; n < 8; ++n) gamma_ip += hexShapeFcn[ip][n] * gamma[n];
        RealType dndx[8][3];
        computeShapeDerivatives(ip, coords, dndx);
        #pragma unroll
        for (int n = 0; n < 8; ++n) {
            RealType diff = -(gamma_ip * (dndx[n][0]*areaVec[0] + dndx[n][1]*areaVec[1] + dndx[n][2]*areaVec[2]));
            lhs[nodeL*8+n] += diff;
            lhs[nodeR*8+n] -= diff;
        }
    }
}

// (A) RECOMPUTE matrix-free: recompute lhs from geometry every matvec, apply,
// scatter. Lightest memory (just area vectors), heaviest compute. d_Au zeroed.
template<typename KeyType, typename RealType, int BlockSize = 256>
__global__ void __launch_bounds__(BlockSize)
cvfem_hex_matfree_apply_diffusion(
    const KeyType* __restrict__ d_c0,const KeyType* __restrict__ d_c1,const KeyType* __restrict__ d_c2,const KeyType* __restrict__ d_c3,
    const KeyType* __restrict__ d_c4,const KeyType* __restrict__ d_c5,const KeyType* __restrict__ d_c6,const KeyType* __restrict__ d_c7,
    size_t numElements,
    const RealType* __restrict__ d_x,const RealType* __restrict__ d_y,const RealType* __restrict__ d_z,
    const RealType* __restrict__ d_gamma,
    const RealType* __restrict__ d_avx,const RealType* __restrict__ d_avy,const RealType* __restrict__ d_avz,
    const int* __restrict__ d_node_to_dof,const uint8_t* __restrict__ d_ownership,int numOwnedRows,
    const RealType* __restrict__ d_u,RealType* __restrict__ d_Au)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    KeyType nd[8] = {d_c0[e],d_c1[e],d_c2[e],d_c3[e],d_c4[e],d_c5[e],d_c6[e],d_c7[e]};
    RealType coords[8][3], gamma[8], uloc[8]; int dofs[8]; uint8_t own[8];
    #pragma unroll
    for (int n = 0; n < 8; ++n) {
        KeyType k = nd[n];
        coords[n][0]=d_x[k]; coords[n][1]=d_y[k]; coords[n][2]=d_z[k]; gamma[n]=d_gamma[k];
        int dof=d_node_to_dof[k]; dofs[n]=dof; own[n]=d_ownership[k]; uloc[n]=(dof>=0)?d_u[dof]:RealType(0);
    }
    RealType lhs[64];
    cvfem_hex_diffusion_lhs<KeyType,RealType>(coords,gamma,d_avx,d_avy,d_avz,e,lhs);
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        if (own[i]==0 || dofs[i]<0 || dofs[i]>=numOwnedRows) continue;
        RealType yi = 0.0;
        #pragma unroll
        for (int j = 0; j < 8; ++j) yi += lhs[i*8+j]*uloc[j];
        atomicAdd(&d_Au[dofs[i]], yi);
    }
}

// (B) Nalu-style STORE-operator build: compute each element's 8x8 lhs ONCE into
// d_Ke[numElements*64] (row-major). Pay the shape-function cost once, not per matvec.
template<typename KeyType, typename RealType, int BlockSize = 256>
__global__ void __launch_bounds__(BlockSize)
cvfem_hex_matfree_build_lhs(
    const KeyType* __restrict__ d_c0,const KeyType* __restrict__ d_c1,const KeyType* __restrict__ d_c2,const KeyType* __restrict__ d_c3,
    const KeyType* __restrict__ d_c4,const KeyType* __restrict__ d_c5,const KeyType* __restrict__ d_c6,const KeyType* __restrict__ d_c7,
    size_t numElements,
    const RealType* __restrict__ d_x,const RealType* __restrict__ d_y,const RealType* __restrict__ d_z,
    const RealType* __restrict__ d_gamma,
    const RealType* __restrict__ d_avx,const RealType* __restrict__ d_avy,const RealType* __restrict__ d_avz,
    RealType* __restrict__ d_Ke)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    KeyType nd[8] = {d_c0[e],d_c1[e],d_c2[e],d_c3[e],d_c4[e],d_c5[e],d_c6[e],d_c7[e]};
    RealType coords[8][3], gamma[8];
    #pragma unroll
    for (int n = 0; n < 8; ++n) {
        KeyType k = nd[n];
        coords[n][0]=d_x[k]; coords[n][1]=d_y[k]; coords[n][2]=d_z[k]; gamma[n]=d_gamma[k];
    }
    RealType lhs[64];
    cvfem_hex_diffusion_lhs<KeyType,RealType>(coords,gamma,d_avx,d_avy,d_avz,e,lhs);
    #pragma unroll
    for (int i = 0; i < 64; ++i) d_Ke[e*64+i] = lhs[i];
}

// (B) Nalu-style STORE-operator apply (element-by-element): y_e = Ke_e * u_e,
// scatter. Reads the stored operator (regular/coalesced) -- no recompute, no
// global matrix. d_Au zeroed before launch.
template<typename KeyType, typename RealType, int BlockSize = 256>
__global__ void __launch_bounds__(BlockSize)
cvfem_hex_ebe_apply(
    const KeyType* __restrict__ d_c0,const KeyType* __restrict__ d_c1,const KeyType* __restrict__ d_c2,const KeyType* __restrict__ d_c3,
    const KeyType* __restrict__ d_c4,const KeyType* __restrict__ d_c5,const KeyType* __restrict__ d_c6,const KeyType* __restrict__ d_c7,
    size_t numElements, const RealType* __restrict__ d_Ke,
    const int* __restrict__ d_node_to_dof,const uint8_t* __restrict__ d_ownership,int numOwnedRows,
    const RealType* __restrict__ d_u,RealType* __restrict__ d_Au)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    KeyType nd[8] = {d_c0[e],d_c1[e],d_c2[e],d_c3[e],d_c4[e],d_c5[e],d_c6[e],d_c7[e]};
    int dofs[8]; uint8_t own[8]; RealType uloc[8];
    #pragma unroll
    for (int n = 0; n < 8; ++n) {
        int dof=d_node_to_dof[nd[n]]; dofs[n]=dof; own[n]=d_ownership[nd[n]]; uloc[n]=(dof>=0)?d_u[dof]:RealType(0);
    }
    const RealType* lhs = d_Ke + e*64;
    #pragma unroll
    for (int i = 0; i < 8; ++i) {
        if (own[i]==0 || dofs[i]<0 || dofs[i]>=numOwnedRows) continue;
        RealType yi = 0.0;
        #pragma unroll
        for (int j = 0; j < 8; ++j) yi += lhs[i*8+j]*uloc[j];
        atomicAdd(&d_Au[dofs[i]], yi);
    }
}

} // namespace fem
} // namespace mars
