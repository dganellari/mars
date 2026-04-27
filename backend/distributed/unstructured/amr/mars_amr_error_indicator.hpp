#pragma once

#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/extrema.h>
#include <thrust/transform.h>
#include <cmath>
#include <mpi.h>

#include "cstone/cuda/cuda_utils.hpp"

namespace mars
{
namespace amr
{

// Hex8: gradient-jump error proxy via least-squares fit over 8 nodes
template<typename KeyType, typename RealType>
__global__ void computeGradientErrorHexKernel(const KeyType* conn0,
                                              const KeyType* conn1,
                                              const KeyType* conn2,
                                              const KeyType* conn3,
                                              const KeyType* conn4,
                                              const KeyType* conn5,
                                              const KeyType* conn6,
                                              const KeyType* conn7,
                                              const RealType* solution,
                                              const RealType* nodeX,
                                              const RealType* nodeY,
                                              const RealType* nodeZ,
                                              const int* nodeToDof,
                                              RealType* errorPerElement,
                                              size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType nodes[8] = {conn0[e], conn1[e], conn2[e], conn3[e], conn4[e], conn5[e], conn6[e], conn7[e]};

    RealType x[8], y[8], z[8], u[8];
    for (int i = 0; i < 8; ++i)
    {
        KeyType n = nodes[i];
        x[i]      = nodeX[n];
        y[i]      = nodeY[n];
        z[i]      = nodeZ[n];
        int dof   = nodeToDof[n];
        u[i]      = (dof >= 0) ? solution[dof] : RealType(0);
    }

    RealType cx = 0, cy = 0, cz = 0, uAvg = 0;
    for (int i = 0; i < 8; ++i)
    {
        cx += x[i]; cy += y[i]; cz += z[i]; uAvg += u[i];
    }
    cx *= RealType(0.125); cy *= RealType(0.125); cz *= RealType(0.125); uAvg *= RealType(0.125);

    // Least-squares gradient: A * grad = b, where A = sum(dx_i dx_i^T), b = sum(dx_i * du_i)
    RealType A00 = 0, A01 = 0, A02 = 0, A11 = 0, A12 = 0, A22 = 0;
    RealType b0 = 0, b1 = 0, b2 = 0;
    for (int i = 0; i < 8; ++i)
    {
        RealType dx = x[i] - cx, dy = y[i] - cy, dz = z[i] - cz, du = u[i] - uAvg;
        A00 += dx * dx; A01 += dx * dy; A02 += dx * dz;
        A11 += dy * dy; A12 += dy * dz; A22 += dz * dz;
        b0 += dx * du; b1 += dy * du; b2 += dz * du;
    }

    // Cramer's rule for the 3x3 normal equations
    RealType det = A00 * (A11 * A22 - A12 * A12) - A01 * (A01 * A22 - A12 * A02) + A02 * (A01 * A12 - A11 * A02);

    // Element-size threshold: scale by element extent so threshold tracks h^4
    // (det of A scales like h^4 for 8 corner-relative offsets in 3D)
    RealType hx_est = x[1] - x[0], hy_est = y[3] - y[0], hz_est = z[4] - z[0];
    RealType hScale = fabs(hx_est * hy_est * hz_est);
    RealType detTol = hScale * hScale * RealType(1e-12);
    if (detTol < RealType(1e-30)) detTol = RealType(1e-30);

    RealType gradX = 0, gradY = 0, gradZ = 0;
    if (fabs(det) > detTol)
    {
        RealType invDet = RealType(1) / det;
        gradX = ((A11 * A22 - A12 * A12) * b0 + (A02 * A12 - A01 * A22) * b1 + (A01 * A12 - A02 * A11) * b2) * invDet;
        gradY = ((A02 * A12 - A01 * A22) * b0 + (A00 * A22 - A02 * A02) * b1 + (A01 * A02 - A00 * A12) * b2) * invDet;
        gradZ = ((A01 * A12 - A02 * A11) * b0 + (A01 * A02 - A00 * A12) * b1 + (A00 * A11 - A01 * A01) * b2) * invDet;
    }

    RealType gradMag2 = gradX * gradX + gradY * gradY + gradZ * gradZ;
    // Guard: if Cramer produced inf/nan, skip this element (mark as zero error)
    if (!isfinite(gradMag2))
    {
        errorPerElement[e] = RealType(0);
        return;
    }
    RealType gradMag = sqrt(gradMag2);

    RealType h = cbrt(hScale);
    if (h < RealType(1e-30)) h = RealType(1e-15);

    RealType err = h * h * gradMag;
    errorPerElement[e] = isfinite(err) ? err : RealType(0);
}

// Tet4: exact gradient via Jacobian inverse (gradient is constant in linear tets)
template<typename KeyType, typename RealType>
__global__ void computeGradientErrorTetKernel(const KeyType* conn0,
                                              const KeyType* conn1,
                                              const KeyType* conn2,
                                              const KeyType* conn3,
                                              const RealType* solution,
                                              const RealType* nodeX,
                                              const RealType* nodeY,
                                              const RealType* nodeZ,
                                              const int* nodeToDof,
                                              RealType* errorPerElement,
                                              size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    KeyType n[4] = {conn0[e], conn1[e], conn2[e], conn3[e]};

    RealType x[4], y[4], z[4], u[4];
    for (int i = 0; i < 4; ++i)
    {
        x[i] = nodeX[n[i]];
        y[i] = nodeY[n[i]];
        z[i] = nodeZ[n[i]];
        int dof = nodeToDof[n[i]];
        u[i] = (dof >= 0) ? solution[dof] : RealType(0);
    }

    // Jacobian columns: edge vectors from node 0
    RealType dx1 = x[1] - x[0], dy1 = y[1] - y[0], dz1 = z[1] - z[0];
    RealType dx2 = x[2] - x[0], dy2 = y[2] - y[0], dz2 = z[2] - z[0];
    RealType dx3 = x[3] - x[0], dy3 = y[3] - y[0], dz3 = z[3] - z[0];

    // det(J) = 6 * volume
    RealType det = dx1 * (dy2 * dz3 - dy3 * dz2)
                 - dy1 * (dx2 * dz3 - dx3 * dz2)
                 + dz1 * (dx2 * dy3 - dx3 * dy2);

    RealType gradX = 0, gradY = 0, gradZ = 0;
    if (fabs(det) > RealType(1e-30))
    {
        RealType invDet = RealType(1) / det;
        RealType du1 = u[1] - u[0], du2 = u[2] - u[0], du3 = u[3] - u[0];

        // grad(u) = J^{-T} * [du1, du2, du3]
        gradX = ((dy2 * dz3 - dy3 * dz2) * du1 + (dy3 * dz1 - dy1 * dz3) * du2 + (dy1 * dz2 - dy2 * dz1) * du3) * invDet;
        gradY = ((dx3 * dz2 - dx2 * dz3) * du1 + (dx1 * dz3 - dx3 * dz1) * du2 + (dx2 * dz1 - dx1 * dz2) * du3) * invDet;
        gradZ = ((dx2 * dy3 - dx3 * dy2) * du1 + (dx3 * dy1 - dx1 * dy3) * du2 + (dx1 * dy2 - dx2 * dy1) * du3) * invDet;
    }

    // h = (6V)^(1/3), V = det/6
    RealType h = cbrt(fabs(det) / RealType(6));
    if (h < RealType(1e-30)) h = RealType(1e-15);

    RealType gradMag   = sqrt(gradX * gradX + gradY * gradY + gradZ * gradZ);
    errorPerElement[e] = h * h * gradMag;
}

template<typename RealType>
__global__ void markElementsKernel(const RealType* errorPerElement,
                                   RealType refineThreshold,
                                   RealType coarsenThreshold,
                                   uint8_t* markFlags, // 1 = refine, 2 = coarsen, 0 = keep
                                   size_t numElements)
{
    size_t e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= numElements) return;

    RealType err = errorPerElement[e];
    if (err > refineThreshold)
        markFlags[e] = 1;
    else if (err < coarsenThreshold)
        markFlags[e] = 2;
    else
        markFlags[e] = 0;
}

// Hex8 error indicator
template<typename KeyType, typename RealType>
class HexErrorIndicator
{
public:
    static cstone::DeviceVector<RealType> computeError(const KeyType* conn0, const KeyType* conn1,
                                                       const KeyType* conn2, const KeyType* conn3,
                                                       const KeyType* conn4, const KeyType* conn5,
                                                       const KeyType* conn6, const KeyType* conn7,
                                                       const RealType* solution,
                                                       const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
                                                       const int* nodeToDof,
                                                       size_t numElements, int blockSize = 256)
    {
        cstone::DeviceVector<RealType> d_error(numElements);
        int numBlocks = (numElements + blockSize - 1) / blockSize;
        computeGradientErrorHexKernel<KeyType, RealType><<<numBlocks, blockSize>>>(
            conn0, conn1, conn2, conn3, conn4, conn5, conn6, conn7,
            solution, nodeX, nodeY, nodeZ, nodeToDof, d_error.data(), numElements);
        cudaDeviceSynchronize();
        return d_error;
    }
};

// Tet4 error indicator
template<typename KeyType, typename RealType>
class TetErrorIndicator
{
public:
    static cstone::DeviceVector<RealType> computeError(const KeyType* conn0, const KeyType* conn1,
                                                       const KeyType* conn2, const KeyType* conn3,
                                                       const RealType* solution,
                                                       const RealType* nodeX, const RealType* nodeY, const RealType* nodeZ,
                                                       const int* nodeToDof,
                                                       size_t numElements, int blockSize = 256)
    {
        cstone::DeviceVector<RealType> d_error(numElements);
        int numBlocks = (numElements + blockSize - 1) / blockSize;
        computeGradientErrorTetKernel<KeyType, RealType><<<numBlocks, blockSize>>>(
            conn0, conn1, conn2, conn3,
            solution, nodeX, nodeY, nodeZ, nodeToDof, d_error.data(), numElements);
        cudaDeviceSynchronize();
        return d_error;
    }
};

// Shared utilities (element-type independent)
template<typename KeyType, typename RealType>
class ErrorIndicator
{
public:
    struct Config
    {
        // Doerfler marking: refine if error > refineFraction * maxError
        RealType refineFraction  = 0.3;
        RealType coarsenFraction = 0.03;
        int maxLevel = 5;
        int blockSize = 256;
    };

    static cstone::DeviceVector<uint8_t> markFromError(const cstone::DeviceVector<RealType>& d_error,
                                                       const Config& config = Config{})
    {
        size_t numElements = d_error.size();
        auto err_b         = thrust::device_pointer_cast(d_error.data());
        auto err_e         = thrust::device_pointer_cast(d_error.data() + numElements);
        auto maxIt         = thrust::max_element(thrust::device, err_b, err_e);
        RealType maxError  = 0;
        cudaMemcpy(&maxError, thrust::raw_pointer_cast(&*maxIt), sizeof(RealType), cudaMemcpyDeviceToHost);

        RealType refineThreshold  = config.refineFraction * maxError;
        RealType coarsenThreshold = config.coarsenFraction * maxError;

        cstone::DeviceVector<uint8_t> d_marks(numElements);
        int numBlocks = (numElements + config.blockSize - 1) / config.blockSize;
        markElementsKernel<<<numBlocks, config.blockSize>>>(
            d_error.data(), refineThreshold, coarsenThreshold, d_marks.data(), numElements);
        cudaDeviceSynchronize();
        return d_marks;
    }

    static RealType globalErrorNorm(cstone::DeviceVector<RealType>& d_error, MPI_Comm comm = MPI_COMM_WORLD)
    {
        auto eb = thrust::device_pointer_cast(d_error.data());
        auto ee = thrust::device_pointer_cast(d_error.data() + d_error.size());
        RealType localSum =
            thrust::transform_reduce(thrust::device, eb, ee,
                                     [] __device__(RealType e) -> RealType { return e * e; },
                                     RealType(0), thrust::plus<RealType>());

        RealType globalSum;
        MPI_Allreduce(&localSum, &globalSum, 1,
                       std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT,
                       MPI_SUM, comm);

        return std::sqrt(globalSum);
    }

    // Variant that sums only over the OWNED element range [startIdx, endIdx)
    // to avoid double-counting halo elements across ranks. Required when
    // d_error includes halos (e.g. computed via per-element kernel over the
    // full local element array).
    static RealType globalErrorNormOwned(cstone::DeviceVector<RealType>& d_error,
                                          size_t startIdx, size_t endIdx,
                                          MPI_Comm comm = MPI_COMM_WORLD)
    {
        auto eb = thrust::device_pointer_cast(d_error.data() + startIdx);
        auto ee = thrust::device_pointer_cast(d_error.data() + endIdx);
        RealType localSum =
            thrust::transform_reduce(thrust::device, eb, ee,
                                     [] __device__(RealType e) -> RealType { return e * e; },
                                     RealType(0), thrust::plus<RealType>());

        RealType globalSum;
        MPI_Allreduce(&localSum, &globalSum, 1,
                       std::is_same_v<RealType, double> ? MPI_DOUBLE : MPI_FLOAT,
                       MPI_SUM, comm);

        return std::sqrt(globalSum);
    }
};

} // namespace amr
} // namespace mars
