#pragma once

#include "../fem/mars_sparse_matrix.hpp"
#include "mars_cg_solver_with_preconditioner.hpp"   // Preconditioner base (pluggable + flexible GMRES)
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
    using Precond = Preconditioner<RealType, IndexType, AcceleratorTag>;

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
        if (spmvBuffer_) cudaFree(spmvBuffer_);
    }

    bool solve(const Matrix& A, const Vector& b, Vector& x, bool usePreconditioner = false)
    {
        // flexible GMRES path (Z-basis) for a pluggable / non-stationary preconditioner; the
        // stationary path below is left untouched so existing callers are unaffected.
        if (flexible_ && precond_) return solveFlexible(A, b, x);

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
            lastIters_ = 0;
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
                lastIters_ = totalIter;
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
                lastIters_ = totalIter;
                return true;
            }
        }

        lastIters_ = totalIter;
        if (verbose_) {
            std::cout << "GMRES did not converge in " << maxIter_ << " iterations\n";
        }
        return false;
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }
    void setMaxIterations(int maxIter) { maxIter_ = maxIter; }
    void setTolerance(RealType tol) { tolerance_ = tol; }
    void setRestart(int restart) { restart_ = restart; }
    void setPreconditioner(Precond* M) { precond_ = M; }   // opt-in pluggable preconditioner
    void setFlexible(bool f) { flexible_ = f; }            // FGMRES (Z-basis) -- needs a preconditioner set
    int getLastIterations() const { return lastIters_; }

private:
    int maxIter_;
    RealType tolerance_;
    int restart_;
    bool verbose_;
    cublasHandle_t cublasHandle_;
    cusparseHandle_t cusparseHandle_;
    void* spmvBuffer_ = nullptr;        // persistent cusparse SpMV workspace, grown on demand -> no per-iter cudaMalloc/Free
    size_t spmvBufferBytes_ = 0;
    Precond* precond_ = nullptr;   // not owned; nullptr -> stationary path (inline Jacobi / none)
    bool flexible_ = false;
    int lastIters_ = 0;

    // Flexible GMRES: store the preconditioned basis Z_j = M^-1 V_j, build the Krylov space on
    // A*Z_j, and reconstruct x = x0 + Z*y. This captures a varying/non-stationary preconditioner
    // (the ACM V-cycle) per-column, where stationary GMRES (fixed-M assumption) stagnates.
    bool solveFlexible(const Matrix& A, const Vector& b, Vector& x)
    {
        size_t n = A.numRows();
        int m = restart_;
        if (x.size() != n) x.resize(n);
        thrust::fill(thrust::device_pointer_cast(x.data()),
                    thrust::device_pointer_cast(x.data() + x.size()), RealType(0));

        precond_->setup(A);   // build / refresh M once

        RealType b_norm = std::sqrt(dot(b, b));
        if (b_norm < 1e-14) { lastIters_ = 0; if (verbose_) std::cout << "FGMRES: RHS is zero\n"; return true; }

        // V = orthonormal Krylov basis stored CONTIGUOUSLY (n x (m+1), column-major) so the Arnoldi
        // orthogonalization is batched cublas gemv (CGS2) -- O(1) syncs/iter, not the j+1 host-synced
        // dots of MGS (nsys: ~57% of the solve). Z (preconditioned basis) stays per-column (only the
        // O(m) reconstruction axpy touches it, not the inner loop).
        Vector Vmat; Vmat.resize(n * static_cast<size_t>(m + 1));
        std::vector<Vector> Z(m);
        for (int i = 0; i < m; ++i)  Z[i].resize(n);
        Vector hco(m + 1), gco(m + 1), vj(n);                  // CGS2 coeffs (device) + current-column scratch for apply()
        std::vector<RealType> hh(m + 1);                       // host copy of one column's coefficients
        auto col = [&](int i) { return thrust::raw_pointer_cast(Vmat.data()) + static_cast<size_t>(i) * n; };
        std::vector<std::vector<RealType>> H(m + 1, std::vector<RealType>(m, 0.0));
        std::vector<RealType> s(m + 1), cs(m), sn(m);
        Vector r(n), w(n);
        int totalIter = 0;

        for (int outer = 0; outer < maxIter_ / m; ++outer)
        {
            spmv(A, x, r);                                          // r = A x
            thrust::transform(thrust::device_pointer_cast(b.data()),
                             thrust::device_pointer_cast(b.data() + b.size()),
                             thrust::device_pointer_cast(r.data()),
                             thrust::device_pointer_cast(r.data()),
                             thrust::minus<RealType>());            // r = b - A x  (NOT preconditioned: flexible convention)
            RealType beta = std::sqrt(dot(r, r));
            if (beta / b_norm < tolerance_) {
                lastIters_ = totalIter;
                if (verbose_) std::cout << "FGMRES converged in " << totalIter << " iters, res=" << beta / b_norm << "\n";
                return true;
            }
            cudaMemcpyAsync(col(0), thrust::raw_pointer_cast(r.data()),
                            n * sizeof(RealType), cudaMemcpyDeviceToDevice, 0);     // V[:,0] = r
            colScale(static_cast<int>(n), RealType(1) / beta, col(0));              // V[:,0] /= beta
            s[0] = beta; for (int i = 1; i <= m; ++i) s[i] = 0.0;

            int j;
            for (j = 0; j < m; ++j, ++totalIter)
            {
                cudaMemcpyAsync(thrust::raw_pointer_cast(vj.data()), col(j),
                                n * sizeof(RealType), cudaMemcpyDeviceToDevice, 0); // col(j) -> Vector for apply()
                precond_->apply(vj, Z[j]);                         // Z[j] = M^-1 V[j]
                spmv(A, Z[j], w);                                  // w = A Z[j]
                // batched CGS2 vs V[0..j]: 4 gemv (async) + 1 D2H of the coefficients, not j+1 synced dots
                const int k = j + 1;
                RealType* Vp = col(0);
                RealType* wp = thrust::raw_pointer_cast(w.data());
                RealType* hp = thrust::raw_pointer_cast(hco.data());
                RealType* gp = thrust::raw_pointer_cast(gco.data());
                gemvT  (static_cast<int>(n), k, Vp, wp, hp);       // h  = V^T w
                gemvNsub(static_cast<int>(n), k, Vp, hp, wp);      // w -= V h
                gemvT  (static_cast<int>(n), k, Vp, wp, gp);       // g  = V^T w   (reorthogonalize for stability)
                gemvNsub(static_cast<int>(n), k, Vp, gp, wp);      // w -= V g
                axpyDev(k, RealType(1), gp, hp);                   // h += g  (total projection coefficients)
                cudaMemcpy(hh.data(), hp, k * sizeof(RealType), cudaMemcpyDeviceToHost);   // ONE D2H per iter
                for (int i = 0; i < k; ++i) H[i][j] = hh[i];
                H[j + 1][j] = std::sqrt(dot(w, w));
                if (std::abs(H[j + 1][j]) < 1e-14) { m = j + 1; break; }
                cudaMemcpyAsync(col(j + 1), wp, n * sizeof(RealType), cudaMemcpyDeviceToDevice, 0);  // V[:,j+1] = w
                colScale(static_cast<int>(n), RealType(1) / H[j + 1][j], col(j + 1));                // V[:,j+1] /= H[j+1][j]
                for (int i = 0; i < j; ++i) {
                    RealType t = cs[i] * H[i][j] + sn[i] * H[i + 1][j];
                    H[i + 1][j] = -sn[i] * H[i][j] + cs[i] * H[i + 1][j];
                    H[i][j] = t;
                }
                RealType hn = std::sqrt(H[j][j] * H[j][j] + H[j + 1][j] * H[j + 1][j]);
                cs[j] = H[j][j] / hn; sn[j] = H[j + 1][j] / hn;
                H[j][j] = cs[j] * H[j][j] + sn[j] * H[j + 1][j]; H[j + 1][j] = 0.0;
                s[j + 1] = -sn[j] * s[j]; s[j] = cs[j] * s[j];
                if (std::abs(s[j + 1]) / b_norm < tolerance_) { m = j + 1; break; }
            }

            std::vector<RealType> y(m);
            for (int i = m - 1; i >= 0; --i) { y[i] = s[i]; for (int k = i + 1; k < m; ++k) y[i] -= H[i][k] * y[k]; y[i] /= H[i][i]; }
            for (int i = 0; i < m; ++i) axpy(y[i], Z[i], x, x);    // x += Z*y  (flexible: accumulate the PRECONDITIONED basis)

            if (std::abs(s[m]) / b_norm < tolerance_) {
                lastIters_ = totalIter;
                if (verbose_) std::cout << "FGMRES converged in " << totalIter << " iters\n";
                return true;
            }
            m = restart_;   // restore restart depth for the next cycle (a break may have shrunk it)
        }
        lastIters_ = totalIter;
        if (verbose_) std::cout << "FGMRES did not converge in " << maxIter_ << " iters\n";
        return false;
    }

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
            if (bufferSize > spmvBufferBytes_) {   // grow once (matrix size is constant -> alloc on the first SpMV, then reused)
                if (spmvBuffer_) cudaFree(spmvBuffer_);
                cudaMalloc(&spmvBuffer_, bufferSize); spmvBufferBytes_ = bufferSize;
            }
            cusparseSpMV(cusparseHandle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY, CUDA_R_64F,
                         CUSPARSE_SPMV_ALG_DEFAULT, spmvBuffer_);
        } else {
            cusparseCreateCsr(&matA, A.numRows(), A.numCols(), A.nnz(), (void*)A.rowOffsetsPtr(),
                              (void*)A.colIndicesPtr(), (void*)A.valuesPtr(), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                              CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
            cusparseCreateDnVec(&vecX, x.size(), (void*)thrust::raw_pointer_cast(x.data()), CUDA_R_32F);
            cusparseCreateDnVec(&vecY, y.size(), (void*)thrust::raw_pointer_cast(y.data()), CUDA_R_32F);

            size_t bufferSize = 0;
            cusparseSpMV_bufferSize(cusparseHandle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY,
                                    CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);
            if (bufferSize > spmvBufferBytes_) {
                if (spmvBuffer_) cudaFree(spmvBuffer_);
                cudaMalloc(&spmvBuffer_, bufferSize); spmvBufferBytes_ = bufferSize;
            }
            cusparseSpMV(cusparseHandle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, vecX, &beta, vecY, CUDA_R_32F,
                         CUSPARSE_SPMV_ALG_DEFAULT, spmvBuffer_);
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

    // batched classical Gram-Schmidt (CGS2) primitives over the contiguous Krylov basis V (n x k, col-major).
    // One gemv replaces k host-synced dots, so the FGMRES Arnoldi is O(1) syncs/iter, not O(k) (nsys: the MGS
    // host-synced dots were ~57% of the solve). gemv output is on-device -> the gemv calls themselves DON'T sync
    // (host-pointer alpha/beta only); only the single D2H of the coefficient column + the norm dot sync.
    // ORDERING: all cublas here runs on the NULL (default) stream, so the per-iter gemv -> coefficient-D2H ->
    // norm-dot chain is correctly serialized. A future cublasSetStream would break that ordering -> add syncs.
    void gemvT(int n, int k, const RealType* V, const RealType* w, RealType* h)   // h = V^T w
    {
        const RealType one = 1, zero = 0;
        if constexpr (std::is_same_v<RealType, double>) cublasDgemv(cublasHandle_, CUBLAS_OP_T, n, k, &one, V, n, w, 1, &zero, h, 1);
        else                                            cublasSgemv(cublasHandle_, CUBLAS_OP_T, n, k, &one, V, n, w, 1, &zero, h, 1);
    }
    void gemvNsub(int n, int k, const RealType* V, const RealType* h, RealType* w) // w -= V h
    {
        const RealType neg = -1, one = 1;
        if constexpr (std::is_same_v<RealType, double>) cublasDgemv(cublasHandle_, CUBLAS_OP_N, n, k, &neg, V, n, h, 1, &one, w, 1);
        else                                            cublasSgemv(cublasHandle_, CUBLAS_OP_N, n, k, &neg, V, n, h, 1, &one, w, 1);
    }
    void axpyDev(int k, RealType a, const RealType* x, RealType* y)                // y += a*x (device, length k)
    {
        if constexpr (std::is_same_v<RealType, double>) cublasDaxpy(cublasHandle_, k, &a, x, 1, y, 1);
        else                                            cublasSaxpy(cublasHandle_, k, &a, x, 1, y, 1);
    }
    void colScale(int n, RealType a, RealType* x)                                 // x *= a (one basis column)
    {
        if constexpr (std::is_same_v<RealType, double>) cublasDscal(cublasHandle_, n, &a, x, 1);
        else                                            cublasSscal(cublasHandle_, n, &a, x, 1);
    }
};

} // namespace fem
} // namespace mars