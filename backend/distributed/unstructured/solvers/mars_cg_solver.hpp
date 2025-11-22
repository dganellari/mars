#pragma once

#include "../fem/mars_sparse_matrix.hpp"
#include <cusparse.h>
#include <cublas_v2.h>
#include <iostream>
#include <cmath>
#include <functional>

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
        , haloExchangeCallback_(nullptr)
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

    // Solve Ax = b for x with Jacobi preconditioning
    // Handles rectangular matrices (rows=owned DOFs, cols=owned+ghost DOFs)
    bool solve(const Matrix& A, const Vector& b, Vector& x)
    {
        size_t m = A.numRows();  // Number of owned DOFs (rows)
        size_t n = A.numCols();  // Number of local DOFs including ghosts (columns)

        // Extract diagonal for preconditioner (only for owned DOFs)
        Vector diag = A.getDiagonal();
        
        // Check for zero or very small diagonal entries
        thrust::host_vector<RealType> h_diag(diag.size());
        thrust::copy(thrust::device_pointer_cast(diag.data()),
                    thrust::device_pointer_cast(diag.data() + diag.size()),
                    h_diag.begin());
        
        for (size_t i = 0; i < h_diag.size(); ++i) {
            if (std::abs(h_diag[i]) < 1e-14) {
                h_diag[i] = 1.0;  // Avoid division by zero
            }
        }
        
        thrust::copy(h_diag.begin(), h_diag.end(),
                    thrust::device_pointer_cast(diag.data()));
        
        // Allocate work vectors
        // r, z, Ap are size m (owned DOFs only)
        // p, x are size n (owned + ghost DOFs for matrix compatibility)
        Vector r(m);   // Residual
        Vector z(m);   // Preconditioned residual
        Vector p(n);   // Search direction (includes ghost DOF slots)
        Vector Ap(m);  // A * p (result is owned DOFs only)

        // Initialize x to zero if not provided
        if (x.size() != n)
        {
            x.resize(n);
        }
        thrust::fill(thrust::device_pointer_cast(x.data()), 
                    thrust::device_pointer_cast(x.data() + x.size()),
                    RealType(0));

        // r = b - Ax  (with x = 0, this is just r = b)
        thrust::copy(thrust::device_pointer_cast(b.data()), 
                    thrust::device_pointer_cast(b.data() + b.size()),
                    thrust::device_pointer_cast(r.data()));

        // z = M^{-1} r (Jacobi: z_i = r_i / A_ii)
        jacobiPrecondition(diag, r, z);

        // p = z (copy to owned portion of p, ghost portion remains zero initially)
        thrust::copy(thrust::device_pointer_cast(z.data()), 
                    thrust::device_pointer_cast(z.data() + z.size()),
                    thrust::device_pointer_cast(p.data()));
        // Ghost portion of p will be filled by halo exchange before first SpMV

        // rho = r^T z
        RealType rho  = dot(r, z);
        RealType rho0 = rho;
        
        // Also compute norm of b for debugging
        RealType b_norm = std::sqrt(dot(b, b));
        
        // Handle case where b is nearly zero
        if (b_norm < 1e-14) {
            if (verbose_) {
                std::cout << "CG: RHS vector is nearly zero (||b|| = " << b_norm << "), solution is zero." << std::endl;
            }
            thrust::fill(thrust::device_pointer_cast(x.data()), 
                        thrust::device_pointer_cast(x.data() + x.size()),
                        RealType(0));
            return true;
        }
        
        if (rho0 < 1e-20) {
            rho0 = b_norm * b_norm;
        }

        // Compute norms for MFEM-style output
        RealType r_norm_0 = std::sqrt(dot(r, r));
        RealType z_norm_0 = std::sqrt(dot(z, z));

        if (verbose_) { 
            std::cout << "Iteration : 0,   ||r||_B = " << z_norm_0 << ",   ||r||_2 = " << r_norm_0 << std::endl; 
        }

        // PCG iterations
        for (int iter = 0; iter < maxIter_; ++iter)
        {
            // Call halo exchange before SpMV if callback is set
            if (haloExchangeCallback_) {
                haloExchangeCallback_(p);
            }
            
            // Ap = A * p
            spmv(A, p, Ap);

            // alpha = rho / (p^T Ap)
            RealType pAp   = dot(p, Ap);
            
            if (std::abs(pAp) < 1e-30) {
                if (verbose_) {
                    std::cout << "CG: p^T Ap = " << pAp << " is too small, stopping." << std::endl;
                }
                return false;
            }
            
            RealType alpha = rho / pAp;

            // x = x + alpha * p
            axpy(alpha, p, x, x);

            // r = r - alpha * Ap
            axpy(-alpha, Ap, r, r);

            // z = M^{-1} r
            jacobiPrecondition(diag, r, z);

            // rho_new = r^T z
            RealType rho_new = dot(r, z);

            // Compute norms for MFEM-style output
            RealType r_norm = std::sqrt(dot(r, r));
            RealType z_norm = std::sqrt(dot(z, z));

            if (verbose_)
            {
                std::cout << "Iteration : " << iter + 1 << ",   ||r||_B = " << z_norm << ",   ||r||_2 = " << r_norm << std::endl;
            }

            // Check convergence
            if (r_norm / b_norm < tolerance_)
            {
                if (verbose_)
                {
                    std::cout << "CG converged in " << iter + 1 << " iterations" << std::endl;
                    RealType avg_reduction = std::pow(r_norm / r_norm_0, 1.0 / (iter + 1));
                    std::cout << "Average reduction factor = " << avg_reduction << std::endl;
                }
                return true;
            }

            // beta = rho_new / rho
            RealType beta = rho_new / rho;

            // p = z + beta * p (only update owned portion, ghost updated by halo exchange)
            // p[0:m] = z + beta * p[0:m]
            axpbyPartial(1.0, z, beta, p, m);

            rho = rho_new;
        }

        if (verbose_) { std::cout << "CG did not converge in " << maxIter_ << " iterations" << std::endl; }
        return false;
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }
    void setMaxIterations(int maxIter) { maxIter_ = maxIter; }
    void setTolerance(RealType tol) { tolerance_ = tol; }
    
    // Set callback for halo exchange before SpMV (for distributed solves)
    void setHaloExchangeCallback(std::function<void(Vector&)> callback) {
        haloExchangeCallback_ = callback;
    }

    // Sparse matrix-vector product: y = A * x (public for debugging)
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

private:

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
    
    // Partial axpby: result[0:n] = alpha * x + beta * result[0:n]
    // Used when x is size n but result is larger (has ghost DOFs)
    void axpbyPartial(RealType alpha, const Vector& x, RealType beta, Vector& result, size_t n)
    {
        // Scale owned portion: result[0:n] *= beta
        if constexpr (std::is_same_v<RealType, double>)
        {
            cublasDscal(cublasHandle_, n, &beta, thrust::raw_pointer_cast(result.data()), 1);
            cublasDaxpy(cublasHandle_, n, &alpha, thrust::raw_pointer_cast(x.data()), 1,
                       thrust::raw_pointer_cast(result.data()), 1);
        }
        else
        {
            cublasSscal(cublasHandle_, n, &beta, thrust::raw_pointer_cast(result.data()), 1);
            cublasSaxpy(cublasHandle_, n, &alpha, thrust::raw_pointer_cast(x.data()), 1,
                       thrust::raw_pointer_cast(result.data()), 1);
        }
    }
    
    // Jacobi preconditioner: z = D^{-1} r
    void jacobiPrecondition(const Vector& diag, const Vector& r, Vector& z)
    {
        thrust::transform(thrust::device_pointer_cast(r.data()),
                         thrust::device_pointer_cast(r.data() + r.size()),
                         thrust::device_pointer_cast(diag.data()),
                         thrust::device_pointer_cast(z.data()),
                         thrust::divides<RealType>());
    }

    int maxIter_;
    RealType tolerance_;
    bool verbose_;
    
    std::function<void(Vector&)> haloExchangeCallback_;

    cublasHandle_t cublasHandle_;
    cusparseHandle_t cusparseHandle_;
};

} // namespace fem
} // namespace mars