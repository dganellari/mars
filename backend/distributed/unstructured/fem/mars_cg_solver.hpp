#pragma once

#include "mars_sparse_matrix.hpp"
#include <cusparse.h>
#include <cublas_v2.h>
#include <iostream>
#include <cmath>

namespace mars
{
namespace fem
{

// Simple Conjugate Gradient solver for GPU
template<typename RealType, typename IndexType, typename AcceleratorTag>
class ConjugateGradientSolver
{
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    ConjugateGradientSolver(int maxIter = 1000, RealType tolerance = 1e-10)
        : maxIter_(maxIter)
        , tolerance_(tolerance)
        , verbose_(true)
    {
        // Initialize cuBLAS and cuSPARSE
        cublasCreate(&cublasHandle_);
        cusparseCreate(&cusparseHandle_);
    }

    ~ConjugateGradientSolver()
    {
        cublasDestroy(cublasHandle_);
        cusparseDestroy(cusparseHandle_);
    }

    // Solve Ax = b for x (initial guess in x)
    bool solve(const Matrix& A, const Vector& b, Vector& x)
    {
        size_t n = A.numRows();

        // Allocate work vectors
        Vector r(n);  // Residual
        Vector p(n);  // Search direction
        Vector Ap(n); // A * p

        // Initialize x to zero if not provided
        if (x.size() != n)
        {
            x.resize(n);
            thrust::fill(thrust::device_pointer_cast(x.data()), thrust::device_pointer_cast(x.data() + x.size()),
                         RealType(0));
        }

        // r = b - Ax
        spmv(A, x, r);       // r = Ax
        axpy(-1.0, r, b, r); // r = b - r

        // p = r
        thrust::copy(thrust::device_pointer_cast(r.data()), thrust::device_pointer_cast(r.data() + r.size()),
                     thrust::device_pointer_cast(p.data()));

        // rho = r^T r
        RealType rho  = dot(r, r);
        RealType rho0 = rho;

        if (verbose_) { std::cout << "CG iteration 0: residual = " << std::sqrt(rho) << std::endl; }

        // CG iterations
        for (int iter = 0; iter < maxIter_; ++iter)
        {
            // Ap = A * p
            spmv(A, p, Ap);

            // alpha = rho / (p^T Ap)
            RealType pAp   = dot(p, Ap);
            RealType alpha = rho / pAp;

            // x = x + alpha * p
            axpy(alpha, p, x, x);

            // r = r - alpha * Ap
            axpy(-alpha, Ap, r, r);

            // rho_new = r^T r
            RealType rho_new = dot(r, r);

            RealType residual = std::sqrt(rho_new / rho0);

            if (verbose_ && (iter % 10 == 0))
            {
                std::cout << "CG iteration " << iter + 1 << ": residual = " << residual << std::endl;
            }

            // Check convergence
            if (residual < tolerance_)
            {
                if (verbose_)
                {
                    std::cout << "CG converged in " << iter + 1 << " iterations, residual = " << residual << std::endl;
                }
                return true;
            }

            // beta = rho_new / rho
            RealType beta = rho_new / rho;

            // p = r + beta * p
            axpby(1.0, r, beta, p, p);

            rho = rho_new;
        }

        if (verbose_) { std::cout << "CG did not converge in " << maxIter_ << " iterations" << std::endl; }
        return false;
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }
    void setMaxIterations(int maxIter) { maxIter_ = maxIter; }
    void setTolerance(RealType tol) { tolerance_ = tol; }

private:
    // Sparse matrix-vector product: y = A * x
    void spmv(const Matrix& A, const Vector& x, Vector& y)
    {
        const RealType alpha = 1.0;
        const RealType beta  = 0.0;

        // Create matrix descriptor
        cusparseSpMatDescr_t matA;
        cusparseDnVecDescr_t vecX, vecY;

        if constexpr (std::is_same_v<RealType, double>)
        {
            cusparseCreateCsr(&matA, A.numRows(), A.numCols(), A.nnz(), (void*)A.rowOffsetsPtr(),
                              (void*)A.colIndicesPtr(), (void*)A.valuesPtr(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                              CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);

            cusparseCreateDnVec(&vecX, x.size(), (void*)thrust::raw_pointer_cast(x.data()), CUDA_R_64F);
            cusparseCreateDnVec(&vecY, y.size(), (void*)thrust::raw_pointer_cast(y.data()), CUDA_R_64F);

            size_t bufferSize = 0;
            cusparseSpMV_bufferSize(cusparseHandle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY,
                                    CUDA_R_64F, CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);

            void* buffer = nullptr;
            cudaMalloc(&buffer, bufferSize);

            cusparseSpMV(cusparseHandle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
                         CUSPARSE_SPMV_ALG_DEFAULT, buffer);

            cudaFree(buffer);
        }
        else
        {
            // Float version
            cusparseCreateCsr(&matA, A.numRows(), A.numCols(), A.nnz(), (void*)A.rowOffsetsPtr(),
                              (void*)A.colIndicesPtr(), (void*)A.valuesPtr(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                              CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);

            cusparseCreateDnVec(&vecX, x.size(), (void*)thrust::raw_pointer_cast(x.data()), CUDA_R_32F);
            cusparseCreateDnVec(&vecY, y.size(), (void*)thrust::raw_pointer_cast(y.data()), CUDA_R_32F);

            size_t bufferSize = 0;
            cusparseSpMV_bufferSize(cusparseHandle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY,
                                    CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);

            void* buffer = nullptr;
            cudaMalloc(&buffer, bufferSize);

            cusparseSpMV(cusparseHandle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY, CUDA_R_32F,
                         CUSPARSE_SPMV_ALG_DEFAULT, buffer);

            cudaFree(buffer);
        }

        cusparseDestroySpMat(matA);
        cusparseDestroyDnVec(vecX);
        cusparseDestroyDnVec(vecY);
    }

    // Dot product: result = x^T y
    RealType dot(const Vector& x, const Vector& y)
    {
        RealType result;
        if constexpr (std::is_same_v<RealType, double>)
        {
            cublasDdot(cublasHandle_, x.size(), thrust::raw_pointer_cast(x.data()), 1,
                       thrust::raw_pointer_cast(y.data()), 1, &result);
        }
        else
        {
            cublasSdot(cublasHandle_, x.size(), thrust::raw_pointer_cast(x.data()), 1,
                       thrust::raw_pointer_cast(y.data()), 1, &result);
        }
        return result;
    }

    // y = alpha * x + y
    void axpy(RealType alpha, const Vector& x, const Vector& y, Vector& result)
    {
        thrust::copy(thrust::device_pointer_cast(y.data()), thrust::device_pointer_cast(y.data() + y.size()),
                     thrust::device_pointer_cast(result.data()));
        if constexpr (std::is_same_v<RealType, double>)
        {
            cublasDaxpy(cublasHandle_, x.size(), &alpha, thrust::raw_pointer_cast(x.data()), 1,
                        thrust::raw_pointer_cast(result.data()), 1);
        }
        else
        {
            cublasSaxpy(cublasHandle_, x.size(), &alpha, thrust::raw_pointer_cast(x.data()), 1,
                        thrust::raw_pointer_cast(result.data()), 1);
        }
    }

    // result = alpha * x + beta * y
    void axpby(RealType alpha, const Vector& x, RealType beta, const Vector& y, Vector& result)
    {
        // result = alpha * x
        if constexpr (std::is_same_v<RealType, double>)
        {
            cublasDcopy(cublasHandle_, x.size(), thrust::raw_pointer_cast(x.data()), 1,
                        thrust::raw_pointer_cast(result.data()), 1);
            cublasDscal(cublasHandle_, x.size(), &alpha, thrust::raw_pointer_cast(result.data()), 1);
        }
        else
        {
            cublasScopy(cublasHandle_, x.size(), thrust::raw_pointer_cast(x.data()), 1,
                        thrust::raw_pointer_cast(result.data()), 1);
            cublasSscal(cublasHandle_, x.size(), &alpha, thrust::raw_pointer_cast(result.data()), 1);
        }

        // result += beta * y
        if constexpr (std::is_same_v<RealType, double>)
        {
            cublasDaxpy(cublasHandle_, y.size(), &beta, thrust::raw_pointer_cast(y.data()), 1,
                        thrust::raw_pointer_cast(result.data()), 1);
        }
        else
        {
            cublasSaxpy(cublasHandle_, y.size(), &beta, thrust::raw_pointer_cast(y.data()), 1,
                        thrust::raw_pointer_cast(result.data()), 1);
        }
    }

    int maxIter_;
    RealType tolerance_;
    bool verbose_;

    cublasHandle_t cublasHandle_;
    cusparseHandle_t cusparseHandle_;
};

} // namespace fem
} // namespace mars
