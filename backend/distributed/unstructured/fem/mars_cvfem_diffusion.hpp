#pragma once

// The diffusion operator -∇·(γ∇u) by the control-volume finite element method (CVFEM) on hexahedra, the
// discretization of the Navier-Stokes solver, assembled into a SparseMatrix of the owned rows: the same
// call shape as TetStiffnessAssembler::assemble for the Galerkin method. Used by tutorial step 4.
#include "../domain.hpp"
#include "mars_cvfem_assembler.hpp"
#include "mars_cvfem_utils.hpp"
#include "mars_sparse_matrix.hpp"
#include "mars_sparsity_builder.hpp"
#include <thrust/execution_policy.h>
#include <thrust/fill.h>

namespace mars
{
namespace fem
{

// A: one row per owned DOF (d_nodeToDof from buildDofMapping: owned DOFs first), columns over all local
// DOFs. Row i is the diffusive flux out of node i's control volume, so A is positive definite once
// Dirichlet rows are set. Every element the rank holds contributes, which completes the owned rows.
template<typename RealType, typename KeyType>
void assembleCvfemDiffusion(const ElementDomain<HexTag, RealType, KeyType, cstone::execution::Gpu>& domain,
                            SparseMatrix<int, RealType, cstone::execution::Gpu>& A, const int* d_nodeToDof,
                            int numOwned, RealType gamma = 1)
{
    using Vector         = cstone::DeviceVector<RealType>;
    const auto& conn     = domain.getElementToNodeConnectivity();
    const size_t nodes   = domain.getNodeCount();
    const size_t elems   = domain.getElementCount();
    const int numDofs    = int(nodes);
    const KeyType* c[8]  = {std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(),
                            std::get<3>(conn).data(), std::get<4>(conn).data(), std::get<5>(conn).data(),
                            std::get<6>(conn).data(), std::get<7>(conn).data()};

    // the 27-point sparsity of every local DOF (first pass counts, second fills); owned rows come first
    cstone::DeviceVector<int> rowPtr(numDofs + 1), diagPtr(numDofs);
    const int nnzAll = CvfemSparsityBuilder<KeyType>::buildFullSparsity(
        c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], elems, d_nodeToDof, numDofs, rowPtr.data(), nullptr, nullptr, 0);
    cstone::DeviceVector<int> colInd(nnzAll);
    CvfemSparsityBuilder<KeyType>::buildFullSparsity(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], elems,
                                                      d_nodeToDof, numDofs, rowPtr.data(), colInd.data(),
                                                      diagPtr.data(), 0);

    // the assembly kernels read the matrix through a descriptor in device memory
    Vector values(nnzAll, RealType(0)), rhs(numDofs, RealType(0));
    CSRMatrix<RealType> h_desc{rowPtr.data(), colInd.data(), values.data(), diagPtr.data(), numDofs, nnzAll};
    h_desc.numOwnedRows = numOwned;
    cstone::DeviceVector<char> d_desc(sizeof(h_desc));
    cudaCheckCall(cudaMemcpy(d_desc.data(), &h_desc, sizeof(h_desc), cudaMemcpyHostToDevice));

    // the general advection-diffusion assembly with diffusion only: γ given, no advection (β = 0, no mass
    // flux), and u = 0 for the explicit terms, whose right-hand side is then zero and not used here
    Vector g(nodes, gamma), zeroNode(nodes, RealType(0)), zeroFace(elems * 12, RealType(0));
    Vector ax(elems * 12), ay(elems * 12), az(elems * 12);
    const RealType *x = domain.getNodeX().data(), *y = domain.getNodeY().data(), *z = domain.getNodeZ().data();
    precomputeAreaVectors<KeyType, RealType>(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], elems, x, y, z,
                                                ax.data(), ay.data(), az.data());
    typename CvfemHexAssembler<KeyType, RealType>::Config config;
    config.variant = CvfemKernelVariant::Tensor;
    CvfemHexAssembler<KeyType, RealType>::assembleFull(
        c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], elems, x, y, z, g.data(), zeroNode.data(), zeroNode.data(),
        zeroNode.data(), zeroNode.data(), zeroNode.data(), zeroFace.data(), ax.data(), ay.data(), az.data(),
        d_nodeToDof, domain.getNodeOwnershipMap().data(), reinterpret_cast<CSRMatrix<RealType>*>(d_desc.data()),
        rhs.data(), config);
    cudaCheckCall(cudaDeviceSynchronize());

    int nnz = 0; // of the owned rows
    cudaCheckCall(cudaMemcpy(&nnz, rowPtr.data() + numOwned, sizeof(int), cudaMemcpyDeviceToHost));
    A.allocate(numOwned, numDofs, nnz);
    cudaCheckCall(cudaMemcpy(A.rowOffsetsPtr(), rowPtr.data(), (numOwned + 1) * sizeof(int), cudaMemcpyDeviceToDevice));
    cudaCheckCall(cudaMemcpy(A.colIndicesPtr(), colInd.data(), size_t(nnz) * sizeof(int), cudaMemcpyDeviceToDevice));
    cudaCheckCall(cudaMemcpy(A.valuesPtr(), values.data(), size_t(nnz) * sizeof(RealType), cudaMemcpyDeviceToDevice));
}

// The control volume of every local node. Node i's part of a hexahedron is the octant of the reference
// cube at corner i, mapped to the element; on a distorted hexahedron the eight parts differ. |det J| of
// the trilinear map is at most quadratic in each reference coordinate, so 2x2x2 Gauss points inside the
// octant integrate it exactly. Complete for the owned nodes, whose elements the rank holds.
template<typename KeyType, typename RealType>
__global__ void cvfemControlVolumeKernel(const KeyType* c0, const KeyType* c1, const KeyType* c2, const KeyType* c3,
                                         const KeyType* c4, const KeyType* c5, const KeyType* c6, const KeyType* c7,
                                         const RealType* x, const RealType* y, const RealType* z, RealType* vol,
                                         size_t numElements)
{
    const size_t e = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (e >= numElements) return;
    const KeyType n[8]     = {c0[e], c1[e], c2[e], c3[e], c4[e], c5[e], c6[e], c7[e]};
    const RealType s[8][3] = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                              {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
    RealType p[8][3];
    for (int k = 0; k < 8; ++k)
    {
        p[k][0] = x[n[k]];
        p[k][1] = y[n[k]];
        p[k][2] = z[n[k]];
    }
    // Gauss points of [0, 1]; the octant of corner i is [0, s_i] in each reference direction
    const RealType g[2] = {RealType(0.5 - 0.5 / 1.7320508075688772), RealType(0.5 + 0.5 / 1.7320508075688772)};
    for (int i = 0; i < 8; ++i)
    {
        RealType v = 0;
        for (int q = 0; q < 8; ++q)
        {
            const RealType a = s[i][0] * g[q & 1], b = s[i][1] * g[(q >> 1) & 1], c = s[i][2] * g[q >> 2];
            RealType J[3][3] = {};
            for (int k = 0; k < 8; ++k)
            {
                const RealType dN[3] = {s[k][0] * (1 + s[k][1] * b) * (1 + s[k][2] * c) / 8,
                                        s[k][1] * (1 + s[k][0] * a) * (1 + s[k][2] * c) / 8,
                                        s[k][2] * (1 + s[k][0] * a) * (1 + s[k][1] * b) / 8};
                for (int r = 0; r < 3; ++r)
                    for (int t = 0; t < 3; ++t)
                        J[r][t] += p[k][r] * dN[t];
            }
            v += fabs(J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                      J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                      J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]));
        }
        // the octant has reference volume 1, so each of the 8 points weighs 1/8
        atomicAdd(&vol[n[i]], v / 8);
    }
}

template<typename RealType, typename KeyType>
void cvfemControlVolumes(const ElementDomain<HexTag, RealType, KeyType, cstone::execution::Gpu>& domain,
                         cstone::DeviceVector<RealType>& vol)
{
    const auto& conn   = domain.getElementToNodeConnectivity();
    const size_t elems = domain.getElementCount();
    vol.resize(domain.getNodeCount());
    thrust::fill(thrust::device, vol.data(), vol.data() + vol.size(), RealType(0));
    if (elems == 0) return;
    cvfemControlVolumeKernel<KeyType, RealType><<<unsigned((elems + 255) / 256), 256>>>(
        std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(), std::get<3>(conn).data(),
        std::get<4>(conn).data(), std::get<5>(conn).data(), std::get<6>(conn).data(), std::get<7>(conn).data(),
        domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(), vol.data(), elems);
    cudaCheckCall(cudaDeviceSynchronize());
}

} // namespace fem
} // namespace mars
