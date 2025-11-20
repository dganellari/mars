#pragma once

#include "mars_sparse_matrix.hpp"
#include <cusparse.h>
#include <cublas_v2.h>
#include <iostream>
#include <cmath>
#include <vector>

namespace mars
{
namespace fem
{

// GMRES solver for GPU - more robust for ill-conditioned systems
template<typename RealType, typename IndexType, typename AcceleratorTag>
class GMRESSolver
{
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    GMRESSolver(int maxIter = 1000, RealType tolerance = 1e-6, int restart = 30)
        : maxIter_(maxIter)
        , tolerance_(tolerance)
        , restart_(restart)
        , verbose_(true)
    {
        cublasCreate(&cublasHandle_);
        cusparseCreate(&cusparseHandle_);
    }

    ~GMRESSolver()
    {
        cublasDestroy(cublasHandle_);
        cusparseDestroy(cusparseHandle_);
    }

    bool solve(const Matrix& A, const Vector& b, Vector& x, bool usePreconditioner = false)
    {
        size_t n = A.numRows();
        int m = restart_;

        // Initialize x to zero
        if (x.size() != n) {
            x.resize(n);
        }
        thrust::fill(thrust::device_pointer_cast(x.data()), 
                    thrust::device_pointer_cast(x.data() + x.size()),
                    RealType(0));

        // Setup preconditioner (Jacobi: M = diag(A))
        Vector M_inv;
        if (usePreconditioner) {
            M_inv = A.getDiagonal();
            // M_inv = 1 / M_inv, with safety for zeros
            thrust::transform(thrust::device_pointer_cast(M_inv.data()),
                            thrust::device_pointer_cast(M_inv.data() + M_inv.size()),
                            thrust::device_pointer_cast(M_inv.data()),
                            [] __device__ (RealType d) { return (std::abs(d) > 1e-14) ? RealType(1) / d : RealType(1); });
        }

        RealType b_norm = std::sqrt(dot(b, b));
        if (b_norm < 1e-14) {
            if (verbose_) {
                std::cout << "GMRES: RHS is zero\n";
            }
            return true;
        }

        // Allocate Krylov subspace basis vectors
        std::vector<Vector> V(m + 1);
        for (int i = 0; i <= m; ++i) {
            V[i].resize(n);
        }

        // Hessenberg matrix (host)
        std::vector<std::vector<RealType>> H(m + 1, std::vector<RealType>(m, 0.0));
        std::vector<RealType> s(m + 1);
        std::vector<RealType> cs(m);
        std::vector<RealType> sn(m);

        Vector r(n);
        Vector w(n);

        int totalIter = 0;

        // GMRES outer iterations
        for (int outer = 0; outer < maxIter_ / m; ++outer)
        {
            // r = b - A*x
            spmv(A, x, r);
            // r = b - r
            thrust::transform(thrust::device_pointer_cast(b.data()),
                             thrust::device_pointer_cast(b.data() + b.size()),
                             thrust::device_pointer_cast(r.data()),
                             thrust::device_pointer_cast(r.data()),
                             thrust::minus<RealType>());

            // Apply preconditioner: r = M^{-1} * r
            if (usePreconditioner) {
                thrust::transform(thrust::device_pointer_cast(r.data()),
                                thrust::device_pointer_cast(r.data() + r.size()),
                                thrust::device_pointer_cast(M_inv.data()),
                                thrust::device_pointer_cast(r.data()),
                                thrust::multiplies<RealType>());
            }

            RealType beta = std::sqrt(dot(r, r));
            
            if (verbose_ && totalIter % 10 == 0) {
                std::cout << "GMRES iteration " << totalIter << ": residual = " << beta / b_norm << std::endl;
            }

            if (beta / b_norm < tolerance_) {
                if (verbose_) {
                    std::cout << "GMRES converged in " << totalIter 
                             << " iterations, residual = " << beta / b_norm << std::endl;
                }
                return true;
            }

            // V[0] = r / beta
            thrust::copy(thrust::device_pointer_cast(r.data()),
                        thrust::device_pointer_cast(r.data() + r.size()),
                        thrust::device_pointer_cast(V[0].data()));
            scale(1.0 / beta, V[0]);

            // Initialize Givens rotation
            s[0] = beta;
            for (int i = 1; i <= m; ++i) {
                s[i] = 0.0;
            }

            // Arnoldi iteration
            int j;
            for (j = 0; j < m; ++j, ++totalIter)
            {
                // w = A * V[j]
                spmv(A, V[j], w);

                // Apply preconditioner: w = M^{-1} * w
                if (usePreconditioner) {
                    thrust::transform(thrust::device_pointer_cast(w.data()),
                                    thrust::device_pointer_cast(w.data() + w.size()),
                                    thrust::device_pointer_cast(M_inv.data()),
                                    thrust::device_pointer_cast(w.data()),
                                    thrust::multiplies<RealType>());
                }

                // Modified Gram-Schmidt orthogonalization
                for (int i = 0; i <= j; ++i) {
                    H[i][j] = dot(w, V[i]);
                    axpy(-H[i][j], V[i], w, w);  // w = w - H[i][j] * V[i]
                }

                H[j + 1][j] = std::sqrt(dot(w, w));

                if (std::abs(H[j + 1][j]) < 1e-14) {
                    // Lucky breakdown
                    m = j + 1;
                    break;
                }

                // V[j+1] = w / H[j+1][j]
                thrust::copy(thrust::device_pointer_cast(w.data()),
                            thrust::device_pointer_cast(w.data() + w.size()),
                            thrust::device_pointer_cast(V[j + 1].data()));
                scale(1.0 / H[j + 1][j], V[j + 1]);

                // Apply previous Givens rotations to H[:, j]
                for (int i = 0; i < j; ++i) {
                    RealType temp = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                    H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                    H[i][j] = temp;
                }

                // Compute new Givens rotation
                RealType h_norm = std::sqrt(H[j][j] * H[j][j] + H[j + 1][j] * H[j + 1][j]);
                cs[j] = H[j][j] / h_norm;
                sn[j] = H[j + 1][j] / h_norm;

                // Apply new Givens rotation
                H[j][j] = cs[j] * H[j][j] + sn[j] * H[j + 1][j];
                H[j + 1][j] = 0.0;

                // Update residual
                s[j + 1] = -sn[j] * s[j];
                s[j] = cs[j] * s[j];

                RealType residual = std::abs(s[j + 1]) / b_norm;

                if (verbose_ && totalIter % 10 == 0) {
                    std::cout << "GMRES iteration " << totalIter << ": residual = " << residual << std::endl;
                }

                if (residual < tolerance_) {
                    m = j + 1;
                    break;
                }
            }

            // Solve upper triangular system H * y = s
            std::vector<RealType> y(m);
            for (int i = m - 1; i >= 0; --i) {
                y[i] = s[i];
                for (int k = i + 1; k < m; ++k) {
                    y[i] -= H[i][k] * y[k];
                }
                y[i] /= H[i][i];
            }

            // Update solution: x = x + V * y
            for (int i = 0; i < m; ++i) {
                axpy(y[i], V[i], x, x);
            }

            if (std::abs(s[m]) / b_norm < tolerance_) {
                if (verbose_) {
                    std::cout << "GMRES converged in " << totalIter 
                             << " iterations, residual = " << std::abs(s[m]) / b_norm << std::endl;
                }
                return true;
            }
        }

        if (verbose_) {
            std::cout << "GMRES did not converge in " << maxIter_ << " iterations\n";
        }
        return false;
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }
    void setMaxIterations(int maxIter) { maxIter_ = maxIter; }
    void setTolerance(RealType tol) { tolerance_ = tol; }
    void setRestart(int restart) { restart_ = restart; }

private:
    int maxIter_;
    RealType tolerance_;
    int restart_;
    bool verbose_;
    cublasHandle_t cublasHandle_;
    cusparseHandle_t cusparseHandle_;

    void spmv(const Matrix& A, const Vector& x, Vector& y)
    {
        const RealType alpha = 1.0;
        const RealType beta  = 0.0;

        cusparseSpMatDescr_t matA;
        cusparseDnVecDescr_t vecX, vecY;

        if constexpr (std::is_same_v<RealType, double>) {
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
        } else {
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

    RealType dot(const Vector& x, const Vector& y)
    {
        RealType result;
        if constexpr (std::is_same_v<RealType, double>) {
            cublasDdot(cublasHandle_, x.size(), thrust::raw_pointer_cast(x.data()), 1,
                       thrust::raw_pointer_cast(y.data()), 1, &result);
        } else {
            cublasSdot(cublasHandle_, x.size(), thrust::raw_pointer_cast(x.data()), 1,
                       thrust::raw_pointer_cast(y.data()), 1, &result);
        }
        return result;
    }

    void axpy(RealType alpha, const Vector& x, const Vector& y, Vector& result)
    {
        if (result.data() != y.data()) {
            thrust::copy(thrust::device_pointer_cast(y.data()), 
                        thrust::device_pointer_cast(y.data() + y.size()),
                        thrust::device_pointer_cast(result.data()));
        }
        if constexpr (std::is_same_v<RealType, double>) {
            cublasDaxpy(cublasHandle_, x.size(), &alpha, thrust::raw_pointer_cast(x.data()), 1,
                        thrust::raw_pointer_cast(result.data()), 1);
        } else {
            cublasSaxpy(cublasHandle_, x.size(), &alpha, thrust::raw_pointer_cast(x.data()), 1,
                        thrust::raw_pointer_cast(result.data()), 1);
        }
    }

    void scale(RealType alpha, Vector& x)
    {
        if constexpr (std::is_same_v<RealType, double>) {
            cublasDscal(cublasHandle_, x.size(), &alpha, thrust::raw_pointer_cast(x.data()), 1);
        } else {
            cublasSscal(cublasHandle_, x.size(), &alpha, thrust::raw_pointer_cast(x.data()), 1);
        }
    }
};

} // namespace fem
} // namespace mars
