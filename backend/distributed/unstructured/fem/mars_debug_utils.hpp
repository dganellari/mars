#pragma once

#include <cstone/cuda/cuda_utils.hpp>
#include <iostream>

namespace mars
{
namespace fem
{
namespace debug
{

// Kernel to count NaNs and Infs in a device vector
template<typename RealType>
__global__ void checkVectorKernel(const RealType* data, size_t size, int* nanCount, int* infCount)
{
    int nans = 0;
    int infs = 0;
    for (size_t i = threadIdx.x + blockIdx.x * blockDim.x; i < size; i += blockDim.x * gridDim.x)
    {
        if (isnan(data[i]))
        {
            atomicAdd(&nans, 1);
        }
        if (isinf(data[i]))
        {
            atomicAdd(&infs, 1);
        }
    }
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        atomicAdd(nanCount, nans);
        atomicAdd(infCount, infs);
    }
}

// Host function to check a device vector for NaNs/Infs
template<typename VectorType>
void checkVector(const std::string& name, const VectorType& vec)
{
    if (vec.size() == 0) return;

    cstone::DeviceVector<int> counts(2);
    thrust::fill(thrust::device_pointer_cast(counts.data()), thrust::device_pointer_cast(counts.data() + 2), 0);

    int* d_counts = thrust::raw_pointer_cast(counts.data());
    checkVectorKernel<<<256, 256>>>(thrust::raw_pointer_cast(vec.data()), vec.size(), d_counts, d_counts + 1);
    cudaDeviceSynchronize();

    std::vector<int> h_counts(2);
    thrust::copy(thrust::device_pointer_cast(counts.data()), thrust::device_pointer_cast(counts.data() + 2), h_counts.begin());

    if (h_counts[0] > 0 || h_counts[1] > 0)
    {
        std::cout << "!!! WARNING: Vector '" << name << "' contains "
                  << h_counts[0] << " NaNs and "
                  << h_counts[1] << " Infs. !!!" << std::endl;
    }
}

// Host function to check a sparse matrix for NaNs/Infs
template<typename MatrixType>
void checkMatrix(const std::string& name, const MatrixType& A)
{
    checkVector(name + " values", A.values());
}

} // namespace debug
} // namespace fem
} // namespace mars
