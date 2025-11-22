#pragma once

#include "../fem/mars_sparse_matrix.hpp"
#include <cusparse.h>
#include <cublas_v2.h>
#include <iostream>
#include <cmath>

namespace mars
{
namespace fem
{

// BiCGSTAB solver for GPU - works for non-symmetric systems
template<typename RealType, typename IndexType, typename AcceleratorTag>
class BiCGSTABSolver
{
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    BiCGSTABSolver(int maxIter = 1000, RealType tolerance = 1e-6)
        : maxIter_(maxIter)
        , tolerance_(tolerance)
        , verbose_(true)
    {
        // Initialize cuBLAS and cuSPARSE
        cublasCreate(&cublasHandle_);
        cusparseCreate(&cusparseHandle_);
    }

    ~BiCGSTABSolver()
    {
        cublasDestroy(cublasHandle_);
        cusparseDestroy(cusparseHandle_);
    }

    // Solve Ax = b for x
    bool solve(const Matrix& A, const Vector& b, Vector& x)
    {
        size_t n = A.numRows();

        // Allocate work vectors
        Vector r(n);      // Residual
        Vector r0(n);     // Initial residual (fixed)
        Vector p(n);      // Search direction
        Vector v(n);      // A*p
        Vector s(n);      // Intermediate residual
        Vector t(n);      // A*s

        // Initialize x to zero
        if (x.size() != n)
        {
            x.resize(n);
        }
        thrust::fill(thrust::device_pointer_cast(x.data()), 
                    thrust::device_pointer_cast(x.data() + x.size()),
                    RealType(0));

        // r = b - Ax (with x=0, r = b)
        thrust::copy(thrust::device_pointer_cast(b.data()), 
                    thrust::device_pointer_cast(b.data() + b.size()),
                    thrust::device_pointer_cast(r.data()));

        // r0 = r (arbitrary choice, could be random)
        thrust::copy(thrust::device_pointer_cast(r.data()), 
                    thrust::device_pointer_cast(r.data() + r.size()),
                    thrust::device_pointer_cast(r0.data()));

        // p = r
        thrust::copy(thrust::device_pointer_cast(r.data()), 
                    thrust::device_pointer_cast(r.data() + r.size()),
                    thrust::device_pointer_cast(p.data()));

        RealType rho = dot(r0, r);
        RealType rho_prev = rho;
        
        RealType b_norm = std::sqrt(dot(b, b));
        RealType r_norm = std::sqrt(dot(r, r));
        
        if (b_norm < 1e-14) {
            if (verbose_) {
                std::cout << "BiCGSTAB: RHS is zero\n";
            }
            return true;
        }

        if (verbose_) { 
            std::cout << "BiCGSTAB iteration 0: residual = " << r_norm / b_norm << std::endl; 
        }

        // BiCGSTAB iterations
        bool first = true;
        for (int iter = 0; iter < maxIter_; ++iter)
        {
            // v = A * p
            spmv(A, p, v);

            // alpha = rho / (r0^T * v)
            RealType r0v = dot(r0, v);
            
            if (std::abs(r0v) < 1e-30) {
                if (verbose_) {
                    std::cout << "BiCGSTAB: r0^T v = " << r0v << " too small, breakdown\n";
                }
                return false;
            }
            
            RealType alpha = rho / r0v;

            // s = r - alpha * v
            thrust::copy(thrust::device_pointer_cast(r.data()), 
                        thrust::device_pointer_cast(r.data() + r.size()),
                        thrust::device_pointer_cast(s.data()));
            axpy(-alpha, v, s, s);

            // Check if we can stop (s is small enough)
            RealType s_norm = std::sqrt(dot(s, s));
            if (s_norm / b_norm < tolerance_) {
                // x = x + alpha * p
                axpy(alpha, p, x, x);
                
                if (verbose_) {
                    std::cout << "BiCGSTAB converged in " << iter + 1 
                             << " iterations (early), residual = " << s_norm / b_norm << std::endl;
                }
                return true;
            }

            // t = A * s
            spmv(A, s, t);

            // omega = (t^T * s) / (t^T * t)
            RealType ts = dot(t, s);
            RealType tt = dot(t, t);
            
            if (std::abs(tt) < 1e-30) {
                if (verbose_) {
                    std::cout << "BiCGSTAB: t^T t = " << tt << " too small, breakdown\n";
                }
                return false;
            }
            
            RealType omega = ts / tt;

            // x = x + alpha * p + omega * s
            axpy(alpha, p, x, x);
            axpy(omega, s, x, x);

            // r = s - omega * t
            thrust::copy(thrust::device_pointer_cast(s.data()), 
                        thrust::device_pointer_cast(s.data() + s.size()),
                        thrust::device_pointer_cast(r.data()));
            axpy(-omega, t, r, r);

            // Check convergence
            r_norm = std::sqrt(dot(r, r));
            RealType residual = r_norm / b_norm;

            if (verbose_ && (iter % 10 == 0))
            {
                std::cout << "BiCGSTAB iteration " << iter + 1 << ": residual = " << residual << std::endl;
            }

            if (residual < tolerance_)
            {
                if (verbose_)
                {
                    std::cout << "BiCGSTAB converged in " << iter + 1 
                             << " iterations, residual = " << residual << std::endl;
                }
                return true;
            }

            // rho_new = r0^T * r
            RealType rho_new = dot(r0, r);
            
            if (std::abs(rho_new) < 1e-30) {
                if (verbose_) {
                    std::cout << "BiCGSTAB: rho = " << rho_new << " too small, breakdown\n";
                }
                return false;
            }

            if (!first) {
                // beta = (rho_new / rho) * (alpha / omega)
                RealType beta = (rho_new / rho) * (alpha / omega);

                // p = r + beta * (p - omega * v)
                Vector p_temp(n);
                thrust::copy(thrust::device_pointer_cast(p.data()), 
                            thrust::device_pointer_cast(p.data() + p.size()),
                            thrust::device_pointer_cast(p_temp.data()));
                axpy(-omega, v, p_temp, p_temp);  // p_temp = p - omega*v
                
                scale(beta, p_temp);               // p_temp = beta * p_temp
                axpy(1.0, p_temp, r, p);           // p = r + p_temp
            }

            first = false;
            rho = rho_new;
        }

        if (verbose_) { 
            std::cout << "BiCGSTAB did not converge in " << maxIter_ << " iterations" << std::endl; 
        }
        return false;
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }
    void setMaxIterations(int maxIter) { maxIter_ = maxIter; }
    void setTolerance(RealType tol) { tolerance_ = tol; }

private:
    int maxIter_;
    RealType tolerance_;
    bool verbose_;
    cublasHandle_t cublasHandle_;
    cusparseHandle_t cusparseHandle_;

    // Sparse matrix-vector product: y = A * x
    void spmv(const Matrix& A, const Vector& x, Vector& y)
    {
        const RealType alpha = 1.0;
        const RealType beta  = 0.0;

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

    // Dot product
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

    // result = y + alpha * x (overwrites result)
    void axpy(RealType alpha, const Vector& x, const Vector& y, Vector& result)
    {
        if (result.data() != y.data()) {
            thrust::copy(thrust::device_pointer_cast(y.data()), 
                        thrust::device_pointer_cast(y.data() + y.size()),
                        thrust::device_pointer_cast(result.data()));
        }
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

    // x = alpha * x
    void scale(RealType alpha, Vector& x)
    {
        if constexpr (std::is_same_v<RealType, double>)
        {
            cublasDscal(cublasHandle_, x.size(), &alpha, thrust::raw_pointer_cast(x.data()), 1);
        }
        else
        {
            cublasSscal(cublasHandle_, x.size(), &alpha, thrust::raw_pointer_cast(x.data()), 1);
        }
    }
};

} // namespace fem
} // namespace mars