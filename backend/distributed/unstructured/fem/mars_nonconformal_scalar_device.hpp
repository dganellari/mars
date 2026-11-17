#pragma once

#include "../utils/mars_gpu_runtime.hpp"
#include "mars_nonconformal_scalar.hpp"
#include "mars_nonconformal_overlap.hpp"
#include <climits>

namespace mars::fem
{
namespace nonconformal_detail
{
template<class T>
__global__ void scalarStencilsKernel(int count,
                                     const InterfacePoint<T>* points,
                                     const InterfaceTetFace<T>* faces,
                                     int faceCount,
                                     const int* offsets,
                                     const int* candidates,
                                     int candidateCount,
                                     InterfacePairingOptions<T> options,
                                     InterfaceStencil<T>* output)
{
    auto index = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= static_cast<size_t>(count)) return;
    int i     = static_cast<int>(index);
    output[i] = {};
    if (offsets[i] < 0 || offsets[i + 1] < offsets[i] || offsets[i + 1] > candidateCount) return;
    output[i] = nonconformalScalarStencil(points[i], faces, faceCount, candidates ? candidates + offsets[i] : nullptr,
                                          offsets[i + 1] - offsets[i], options);
}
} // namespace nonconformal_detail

// Arrays and scratch belong to the caller and remain on the device. A successful
// launch does not imply successful pairing: inspect/reduce the per-point statuses.
template<class T>
cudaError_t buildNonconformalScalarStencils(int count,
                                            const InterfacePoint<T>* points,
                                            const InterfaceTetFace<T>* faces,
                                            int faceCount,
                                            const int* offsets,
                                            const int* candidates,
                                            int candidateCount,
                                            InterfacePairingOptions<T> options,
                                            InterfaceStencil<T>* output,
                                            cudaStream_t stream = 0)
{
    if (count < 0 || faceCount < 0 || candidateCount < 0) return cudaErrorInvalidValue;
    if (count == 0) return cudaSuccess;
    if (!points || !faces || !offsets || !output || (candidateCount && !candidates)) return cudaErrorInvalidValue;
    constexpr int block = 128;
    nonconformal_detail::scalarStencilsKernel<<<(count - 1) / block + 1, block, 0, stream>>>(
        count, points, faces, faceCount, offsets, candidates, candidateCount, options, output);
    return cudaGetLastError();
}
namespace nonconformal_detail
{
template<class T>
__global__ void conservativePairsKernel(int count,
                                        const InterfaceFacePair* pairs,
                                        const InterfaceTetFace<T>* faces,
                                        int faceCount,
                                        InterfacePairingOptions<T> options,
                                        InterfaceStencil<T>* out,
                                        int* counts,
                                        InterfaceStatus* status)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= size_t(count)) return;
    counts[i] = 0;
    status[i] = InterfaceStatus::invalid_input;
    auto p    = pairs[i];
    if (p.first < 0 || p.first >= faceCount || p.second < 0 || p.second >= faceCount) return;
    status[i] =
        nonconformalConservativePair(faces[p.first], faces[p.second], options, out + i * nonconformal_pair_capacity,
                                     nonconformal_pair_capacity, counts[i]);
}
} // namespace nonconformal_detail

// Output has count*72 slots. Only counts[i] entries of each pair are initialized;
// reduce statuses collectively before compacting or assembling any of them.
template<class T>
cudaError_t buildNonconformalConservativePairs(int count,
                                               const InterfaceFacePair* pairs,
                                               const InterfaceTetFace<T>* faces,
                                               int faceCount,
                                               InterfacePairingOptions<T> options,
                                               InterfaceStencil<T>* out,
                                               int* counts,
                                               InterfaceStatus* status,
                                               cudaStream_t stream = 0)
{
    if (count < 0 || count > INT_MAX / nonconformal_pair_capacity || faceCount < 0) return cudaErrorInvalidValue;
    if (count == 0) return cudaSuccess;
    if (!pairs || !faces || !out || !counts || !status) return cudaErrorInvalidValue;
    constexpr int block = 64;
    nonconformal_detail::conservativePairsKernel<<<(count - 1) / block + 1, block, 0, stream>>>(
        count, pairs, faces, faceCount, options, out, counts, status);
    return cudaGetLastError();
}
} // namespace mars::fem
