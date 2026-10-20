#pragma once

// Finite-element assembly with your own element kernel, on tetrahedra and hexahedra. You write a functor that
// computes one element's matrix and vector from its corner coordinates; assembleElements runs it on the GPU
// for every element the rank holds and adds the results into the owned rows of A and b, ready for the
// Krylov solvers. ElementValues gives the functor what it needs at each quadrature point, as deal.II's
// FEValues does: basis values, basis gradients in physical coordinates, the point, and weight times |det J|.
// Tutorial step 5 writes a Poisson operator with it.
#include "../domain.hpp"
#include "mars_reference_element.hpp"
#include "mars_sparse_matrix.hpp"
#include "mars_sparsity_builder.hpp"
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <type_traits>
#include <utility>

namespace mars
{
namespace fem
{

template<typename ElementTag, typename RealType>
struct ElementValues
{
    using Reference                    = ReferenceElement<ElementTag, RealType>;
    static constexpr int numNodes      = Reference::numNodes;
    static constexpr int numQuadPoints = Reference::numQuadPoints;

    RealType phi[numNodes];     // basis function values
    RealType grad[numNodes][3]; // basis function gradients in physical coordinates
    RealType x[3];              // the quadrature point in physical coordinates
    RealType JxW;               // quadrature weight times |det J|

    // the values at quadrature point q of the element with corner coordinates X
    __host__ __device__ void reinit(const RealType (&X)[numNodes][3], int q)
    {
        const auto qp = Reference::getQuadraturePoint(q);
        RealType cx[numNodes], cy[numNodes], cz[numNodes];
        for (int i = 0; i < numNodes; ++i)
        {
            cx[i] = X[i][0];
            cy[i] = X[i][1];
            cz[i] = X[i][2];
        }
        RealType J[3][3], Jinv[3][3];
        if constexpr (std::is_same_v<ElementTag, HexTag>)
            Reference::computeJacobian(cx, cy, cz, J, qp.xi, qp.eta, qp.zeta);
        else
            Reference::computeJacobian(cx, cy, cz, J);
        Reference::computeJacobianInverse(J, Jinv);
        JxW = qp.weight * fabs(Reference::computeJacobianDeterminant(J));

        x[0] = x[1] = x[2] = 0;
        for (int i = 0; i < numNodes; ++i)
        {
            phi[i] = Reference::evaluateBasis(i, qp.xi, qp.eta, qp.zeta);
            RealType g[3];
            if constexpr (std::is_same_v<ElementTag, HexTag>)
                Reference::evaluateGradient(i, g, qp.xi, qp.eta, qp.zeta);
            else
                Reference::evaluateGradient(i, g);
            Reference::transformGradient(Jinv, g, grad[i]);
            for (int k = 0; k < 3; ++k)
                x[k] += phi[i] * X[i][k];
        }
    }
};

template<typename KeyType, int N>
struct CornerColumns
{
    const KeyType* c[N];
};

// One thread per element: gather the corners, call the element kernel, add the owned rows.
template<int N, typename RealType, typename KeyType, typename ElementKernel>
__global__ void assembleElementsKernel(CornerColumns<KeyType, N> conn, size_t numElements, const RealType* x,
                                       const RealType* y, const RealType* z, const int* nodeToDof, int numOwned,
                                       const int* rowPtr, const int* colInd, RealType* values, RealType* rhs,
                                       ElementKernel kernel)
{
    const size_t e = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (e >= numElements) return;

    KeyType n[N];
    RealType X[N][3];
    for (int i = 0; i < N; ++i)
    {
        n[i]    = conn.c[i][e];
        X[i][0] = x[n[i]];
        X[i][1] = y[n[i]];
        X[i][2] = z[n[i]];
    }
    RealType K[N][N] = {}, F[N] = {};
    kernel(X, K, F);

    for (int i = 0; i < N; ++i)
    {
        const int row = nodeToDof[n[i]];
        if (row < 0 || row >= numOwned) continue; // a ghost row: its owner assembles it
        atomicAdd(&rhs[row], F[i]);
        for (int j = 0; j < N; ++j)
        {
            const int col = nodeToDof[n[j]];
            int lo = rowPtr[row], hi = rowPtr[row + 1];
            while (lo < hi) // the columns of a row are sorted
            {
                const int mid = (lo + hi) / 2;
                if (colInd[mid] < col) lo = mid + 1;
                else hi = mid;
            }
            if (lo < rowPtr[row + 1] && colInd[lo] == col) atomicAdd(&values[lo], K[i][j]);
        }
    }
}

// A: one row per owned DOF (d_nodeToDof from buildDofMapping: owned DOFs first), columns over all local
// DOFs; b: one entry per owned DOF. kernel(X, K, F) gets the corner coordinates X[N][3] and adds the element
// matrix into K[N][N] and the element vector into F[N] (both start at zero). Every element the rank holds
// contributes, which completes the owned rows.
template<typename ElementTag, typename RealType, typename KeyType, typename ElementKernel>
void assembleElements(const ElementDomain<ElementTag, RealType, KeyType, cstone::execution::Gpu>& domain,
                      SparseMatrix<int, RealType, cstone::execution::Gpu>& A, cstone::DeviceVector<RealType>& b,
                      const int* d_nodeToDof, int numOwned, ElementKernel kernel)
{
    constexpr int N    = ElementTag::NodesPerElement;
    const auto& conn   = domain.getElementToNodeConnectivity();
    const size_t elems = domain.getElementCount();
    const int numDofs  = int(domain.getNodeCount());
    CornerColumns<KeyType, N> c;
    [&]<std::size_t... I>(std::index_sequence<I...>) { ((c.c[I] = std::get<I>(conn).data()), ...); }
    (std::make_index_sequence<N>{});

    // the sparsity of every local DOF, owned rows first: each element couples all its corners. The first
    // pass counts the columns, the second fills them.
    cstone::DeviceVector<int> rowPtr(numDofs + 1), colInd;
    auto build = [&](int* cols) {
        if constexpr (std::is_same_v<ElementTag, HexTag>)
            return CvfemSparsityBuilder<KeyType>::buildFullSparsity(c.c[0], c.c[1], c.c[2], c.c[3], c.c[4], c.c[5],
                                                                     c.c[6], c.c[7], elems, d_nodeToDof, numDofs,
                                                                     rowPtr.data(), cols, nullptr, 0);
        else
        {
            static_assert(std::is_same_v<ElementTag, TetTag>, "assembleElements supports tetrahedra and hexahedra");
            return CvfemTetSparsityBuilder<KeyType>::buildFullSparsity(c.c[0], c.c[1], c.c[2], c.c[3], elems,
                                                                        d_nodeToDof, numDofs, rowPtr.data(), cols,
                                                                        nullptr, 0);
        }
    };
    colInd.resize(build(nullptr));
    build(colInd.data());

    int nnz = 0; // of the owned rows
    cudaCheckCall(cudaMemcpy(&nnz, rowPtr.data() + numOwned, sizeof(int), cudaMemcpyDeviceToHost));
    A.allocate(numOwned, numDofs, nnz);
    cudaCheckCall(cudaMemcpy(A.rowOffsetsPtr(), rowPtr.data(), (numOwned + 1) * sizeof(int), cudaMemcpyDeviceToDevice));
    cudaCheckCall(cudaMemcpy(A.colIndicesPtr(), colInd.data(), size_t(nnz) * sizeof(int), cudaMemcpyDeviceToDevice));
    A.zero();
    b.resize(numOwned);
    thrust::fill(thrust::device, b.data(), b.data() + numOwned, RealType(0));

    if (elems == 0) return;
    assembleElementsKernel<N, RealType, KeyType><<<unsigned((elems + 255) / 256), 256>>>(
        c, elems, domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(), d_nodeToDof,
        numOwned, A.rowOffsetsPtr(), A.colIndicesPtr(), A.valuesPtr(), b.data(), kernel);
    cudaCheckCall(cudaGetLastError());
    cudaCheckCall(cudaDeviceSynchronize());
}

} // namespace fem
} // namespace mars
