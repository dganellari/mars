#pragma once

#include "mars_cvfem_kernel.hpp"     // Tet4CVFEM reference data
#include "mars_cvfem_hex_kernel.hpp" // mars::fem::CSRMatrix lives here
#include <cuda_runtime.h>

namespace mars {
namespace fem {

// Full CVFEM assembly kernel for linear tetrahedra: scatters the full 4x4
// local element matrix into the global CSR (no diagonal lumping). Sparsity
// pattern must contain all 16 (row, col) pairs per element.
//
// Diffusion here is classical linear-tet stiffness K_ij = volume * gamma_elem
// * (dNdx_i . dNdx_j) with constant dNdx and gamma at the centroid. That is
// NOT the same operator as mars_cvfem_tet_kernel_graph.hpp, which uses the
// per-SCS edge-area flux form (-gamma_ip * dNdx . areaVec). Only the SCS-edge
// advection (upwind + deferred correction) matches across the two kernels.
// mars_cvfem_tet_kernel_full_perip.hpp is the math twin of THIS kernel; the
// two are bit-identical at RealType=double.
template<typename KeyType, typename RealType>
__global__ void cvfem_tet_assembly_kernel_full(
    const KeyType* __restrict__ d_conn0,
    const KeyType* __restrict__ d_conn1,
    const KeyType* __restrict__ d_conn2,
    const KeyType* __restrict__ d_conn3,
    size_t numElements,
    const RealType* __restrict__ d_x,
    const RealType* __restrict__ d_y,
    const RealType* __restrict__ d_z,
    const RealType* __restrict__ d_gamma,
    const RealType* __restrict__ d_phi,
    const RealType* __restrict__ d_beta,
    const RealType* __restrict__ d_grad_phi_x,
    const RealType* __restrict__ d_grad_phi_y,
    const RealType* __restrict__ d_grad_phi_z,
    const RealType* __restrict__ d_mdot,           // [numElements * 6]
    const int* __restrict__ d_node_to_dof,
    const uint8_t* __restrict__ d_ownership,
    CSRMatrix<RealType>* matrix,
    RealType* __restrict__ d_rhs)
{
    int elemIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (elemIdx >= numElements) return;

    constexpr int nnodes = 4;
    constexpr int nscs   = 6;
    constexpr int ndim   = 3;

    KeyType nodes[nnodes];
    nodes[0] = d_conn0[elemIdx];
    nodes[1] = d_conn1[elemIdx];
    nodes[2] = d_conn2[elemIdx];
    nodes[3] = d_conn3[elemIdx];

    RealType coords[nnodes][ndim];
    RealType phi[nnodes], gamma[nnodes], beta[nnodes];
    RealType grad_phi[nnodes][ndim];
    for (int n = 0; n < nnodes; ++n) {
        coords[n][0] = d_x[nodes[n]];
        coords[n][1] = d_y[nodes[n]];
        coords[n][2] = d_z[nodes[n]];
        phi[n]      = d_phi[nodes[n]];
        gamma[n]    = d_gamma[nodes[n]];
        beta[n]     = d_beta[nodes[n]];
        grad_phi[n][0] = d_grad_phi_x[nodes[n]];
        grad_phi[n][1] = d_grad_phi_y[nodes[n]];
        grad_phi[n][2] = d_grad_phi_z[nodes[n]];
    }

    RealType det;
    RealType dNdx[nnodes][ndim];
    Tet4CVFEM::jacobian_and_dNdx<RealType>(coords, det, dNdx);

    RealType lhs[nnodes * nnodes] = {RealType(0)};
    RealType rhs[nnodes]          = {RealType(0)};

    // Linear tet stiffness from classical FEM:
    //   K_ij = volume * gamma_avg * (dNdx_i . dNdx_j)
    // dNdx is constant per element (linear basis), so the integral over the
    // element is just dNdx_i . dNdx_j * vol. Gamma is interpolated at the
    // element centroid (xi=eta=zeta=1/4 -> all N = 1/4) which is the standard
    // 1-point Gauss rule for linear tets.
    RealType volume = fabs(det) / RealType(6);
    RealType gamma_elem = RealType(0.25) * (gamma[0] + gamma[1] + gamma[2] + gamma[3]);

    for (int i = 0; i < nnodes; ++i) {
        for (int j = 0; j < nnodes; ++j) {
            RealType k_ij = RealType(0);
            for (int d = 0; d < ndim; ++d) k_ij += dNdx[i][d] * dNdx[j][d];
            lhs[i * nnodes + j] += volume * gamma_elem * k_ij;
        }
    }

    // Residual form: rhs -= K * phi  (assembleFull pattern, matches hex kernel).
    for (int i = 0; i < nnodes; ++i) {
        RealType r = RealType(0);
        for (int j = 0; j < nnodes; ++j) r += lhs[i * nnodes + j] * phi[j];
        rhs[i] -= r;
    }

    // Advection (upwind on edges, optional). Tet has 6 SCS edges; mdot is
    // per-(element,SCS). Lumped-style: each edge contributes a 2x2 block.
    for (int ip = 0; ip < nscs; ++ip) {
        RealType mdot = d_mdot[elemIdx * nscs + ip];
        if (mdot == RealType(0)) continue;

        int nodeL, nodeR;
        Tet4CVFEM::get_scs_nodes(ip, nodeL, nodeR);

        double xi, eta, zeta;
        Tet4CVFEM::scs_coords(ip, xi, eta, zeta);
        double N[nnodes];
        Tet4CVFEM::shape_fcn(xi, eta, zeta, N);

        RealType coords_ip[ndim] = {RealType(0), RealType(0), RealType(0)};
        for (int n = 0; n < nnodes; ++n)
            for (int d = 0; d < ndim; ++d) coords_ip[d] += N[n] * coords[n][d];

        RealType phi_upwind, beta_upwind, dcorr = RealType(0);
        int upw = (mdot > RealType(0)) ? nodeL : nodeR;
        phi_upwind  = phi[upw];
        beta_upwind = beta[upw];
        for (int d = 0; d < ndim; ++d) {
            RealType dx = coords_ip[d] - coords[upw][d];
            dcorr += grad_phi[upw][d] * dx;
        }
        dcorr *= beta_upwind;
        RealType adv_flux = mdot * (phi_upwind + dcorr);

        rhs[nodeL] -= adv_flux;
        rhs[nodeR] += adv_flux;

        RealType lhsfac_L = RealType(0.5) * (mdot + fabs(mdot));
        RealType lhsfac_R = RealType(0.5) * (mdot - fabs(mdot));
        lhs[nodeL * nnodes + nodeL] += lhsfac_L;
        lhs[nodeR * nnodes + nodeL] -= lhsfac_L;
        lhs[nodeL * nnodes + nodeR] += lhsfac_R;
        lhs[nodeR * nnodes + nodeR] -= lhsfac_R;
    }

    // Global scatter: full 4x4, no lumping. Sparsity must include all pairs.
    for (int i = 0; i < nnodes; ++i) {
        KeyType row_node = nodes[i];
        int row_dof      = d_node_to_dof[row_node];
        uint8_t ownership = d_ownership[row_node];
        if (ownership == 0 || row_dof < 0) continue;

        atomicAdd(&d_rhs[row_dof], rhs[i]);

        for (int j = 0; j < nnodes; ++j) {
            KeyType col_node = nodes[j];
            int col_dof      = d_node_to_dof[col_node];
            if (col_dof < 0) continue;
            matrix->addValue(row_dof, col_dof, lhs[i * nnodes + j]);
        }
    }
}

} // namespace fem
} // namespace mars
