#pragma once

// y = A x through cuSPARSE / hipSPARSE for the Krylov solvers, set up once per solve: the matrix and vector
// descriptors, the work buffer and, on AMD, rocSPARSE's CSR analysis (kept in the matrix descriptor) are
// made once instead of once per product, which on MI300A cost more than the products themselves. Each
// product only points the vector descriptors at its x and y; all x must have one size, all y another.
#include "../utils/mars_gpu_sparse_blas.hpp"
#include "cstone/cuda/device_vector.h"
#include <thrust/memory.h>
#include <type_traits>

namespace mars
{
namespace fem
{

template<typename RealType, typename Matrix, typename Vector>
class SpmvPlan
{
public:
    SpmvPlan(cusparseHandle_t handle, const Matrix& A, const Vector& x, Vector& y)
        : handle_(handle)
    {
        cusparseCreateCsr(&mat_, A.numRows(), A.numCols(), A.nnz(), (void*)A.rowOffsetsPtr(), (void*)A.colIndicesPtr(),
                          (void*)A.valuesPtr(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
                          valueType);
        cusparseCreateDnVec(&x_, x.size(), (void*)thrust::raw_pointer_cast(x.data()), valueType);
        cusparseCreateDnVec(&y_, y.size(), (void*)thrust::raw_pointer_cast(y.data()), valueType);
        size_t bytes = 0;
        cusparseSpMV_bufferSize(handle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &one_, mat_, x_, &zero_, y_, valueType,
                                CUSPARSE_SPMV_ALG_DEFAULT, &bytes);
        buffer_.resize(bytes);
    }
    SpmvPlan(const SpmvPlan&)            = delete;
    SpmvPlan& operator=(const SpmvPlan&) = delete;
    ~SpmvPlan()
    {
        cusparseDestroySpMat(mat_);
        cusparseDestroyDnVec(x_);
        cusparseDestroyDnVec(y_);
    }

    void apply(const Vector& x, Vector& y)
    {
        cusparseDnVecSetValues(x_, (void*)thrust::raw_pointer_cast(x.data()));
        cusparseDnVecSetValues(y_, (void*)thrust::raw_pointer_cast(y.data()));
        cusparseSpMV(handle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &one_, mat_, x_, &zero_, y_, valueType,
                     CUSPARSE_SPMV_ALG_DEFAULT, buffer_.data());
    }

private:
    static constexpr auto valueType = std::is_same_v<RealType, double> ? CUDA_R_64F : CUDA_R_32F;
    const RealType one_ = 1, zero_ = 0;
    cusparseHandle_t handle_;
    cusparseSpMatDescr_t mat_;
    cusparseDnVecDescr_t x_, y_;
    cstone::DeviceVector<char> buffer_;
};

} // namespace fem
} // namespace mars
