#pragma once

#include "../fem/mars_sparse_matrix.hpp"
#include "../utils/mars_gpu_sparse_blas.hpp"
#include "mars_spmv_plan.hpp"
#include <iostream>
#include <cmath>
#include <functional>
#include <memory>

namespace mars
{
namespace fem
{

// Base preconditioner interface
template<typename RealType, typename IndexType, typename AcceleratorTag>
class Preconditioner {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    virtual ~Preconditioner() {}
    virtual void setup(const Matrix& A) = 0;
    virtual void apply(const Vector& r, Vector& z) = 0;
};

// Jacobi preconditioner (existing implementation)
template<typename RealType, typename IndexType, typename AcceleratorTag>
class JacobiPreconditioner : public Preconditioner<RealType, IndexType, AcceleratorTag> {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    void setup(const Matrix& A) override {
        diag_ = A.getDiagonal();

        // Check for zero or very small diagonal entries
        thrust::host_vector<RealType> h_diag(diag_.size());
        thrust::copy(thrust::device_pointer_cast(diag_.data()),
                    thrust::device_pointer_cast(diag_.data() + diag_.size()),
                    h_diag.begin());

        for (size_t i = 0; i < h_diag.size(); ++i) {
            if (std::abs(h_diag[i]) < 1e-14) {
                h_diag[i] = 1.0;  // Avoid division by zero
            }
        }

        thrust::copy(h_diag.begin(), h_diag.end(),
                    thrust::device_pointer_cast(diag_.data()));
    }

    void apply(const Vector& r, Vector& z) override {
        // z = r ./ diag
        thrust::transform(thrust::device_pointer_cast(r.data()),
                         thrust::device_pointer_cast(r.data() + r.size()),
                         thrust::device_pointer_cast(diag_.data()),
                         thrust::device_pointer_cast(z.data()),
                         thrust::divides<RealType>());
    }

private:
    Vector diag_;
};

// Identity preconditioner (z = r). Makes flexible GMRES reduce EXACTLY to plain GMRES, so it serves
// as a convergence-independent unit test of the FGMRES Z-basis reconstruction.
template<typename RealType, typename IndexType, typename AcceleratorTag>
class IdentityPreconditioner : public Preconditioner<RealType, IndexType, AcceleratorTag> {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    void setup(const Matrix&) override {}
    void apply(const Vector& r, Vector& z) override {
        if (z.size() != r.size()) z.resize(r.size());
        thrust::copy(thrust::device_pointer_cast(r.data()),
                     thrust::device_pointer_cast(r.data() + r.size()),
                     thrust::device_pointer_cast(z.data()));
    }
};

// Conjugate Gradient solver with pluggable preconditioners
template<typename RealType, typename IndexType, typename AcceleratorTag>
class PreconditionedConjugateGradientSolver
{
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;
    using Precond = Preconditioner<RealType, IndexType, AcceleratorTag>;

    PreconditionedConjugateGradientSolver(int maxIter = 1000, RealType tolerance = 1e-10)
        : maxIter_(maxIter)
        , tolerance_(tolerance)
        , verbose_(true)
        , haloExchangeCallback_(nullptr)
        , preconditioner_(std::make_shared<JacobiPreconditioner<RealType, IndexType, AcceleratorTag>>())
    {
        // Initialize cuBLAS and cuSPARSE
        cublasCreate(&cublasHandle_);
        cusparseCreate(&cusparseHandle_);
    }

    ~PreconditionedConjugateGradientSolver()
    {
        cublasDestroy(cublasHandle_);
        cusparseDestroy(cusparseHandle_);
    }

    // Set custom preconditioner
    void setPreconditioner(std::shared_ptr<Precond> precond) {
        preconditioner_ = precond;
    }

    // Solve Ax = b for x
    // Handles rectangular matrices (rows=owned DOFs, cols=owned+ghost DOFs)
    bool solve(const Matrix& A, const Vector& b, Vector& x)
    {
        size_t m = A.numRows();  // Number of owned DOFs (rows)
        size_t n = A.numCols();  // Number of local DOFs including ghosts (columns)

        // Setup preconditioner
        preconditioner_->setup(A);

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

        // z = M^{-1} r
        preconditioner_->apply(r, z);

        // p = z (copy to owned portion of p, ghost portion remains zero initially)
        thrust::copy(thrust::device_pointer_cast(z.data()),
                    thrust::device_pointer_cast(z.data() + z.size()),
                    thrust::device_pointer_cast(p.data()));
        // Ghost portion of p will be filled by halo exchange before first SpMV

        // rho = r^T z
        RealType rho  = dot(r, z);
        RealType rho0 = rho;

        // Also compute norm of b for debugging
        RealType norm_b = std::sqrt(dot(b, b));

        if (verbose_) {
            std::cout << "CG: Initial residual norm: " << std::sqrt(rho) << std::endl;
        }

        SpmvPlan<RealType, Matrix, Vector> spmvPlan(cusparseHandle_, A, p, Ap);
        int iter = 0;
        for (iter = 0; iter < maxIter_; ++iter)
        {
            // Perform halo exchange on p if needed
            if (haloExchangeCallback_)
            {
                haloExchangeCallback_(p);
            }

            // Ap = A * p
            spmvPlan.apply(p, Ap);

            // alpha = rho / (p^T * Ap)
            RealType pAp = dot(p, Ap);  // This handles ghost DOFs correctly
            RealType alpha = rho / pAp;

            // x = x + alpha * p
            axpy(alpha, p, x);

            // r = r - alpha * Ap
            axpy(-alpha, Ap, r);

            // z = M^{-1} r
            preconditioner_->apply(r, z);

            // rho_new = r^T z
            RealType rho_new = dot(r, z);

            // Check convergence
            RealType residual_norm = std::sqrt(rho_new);
            if (verbose_ && (iter % 10 == 0 || iter > maxIter_ - 10)) {
                std::cout << "CG iteration " << iter << ": residual = " << residual_norm
                         << " (relative: " << residual_norm / norm_b << ")" << std::endl;
            }

            if (residual_norm < tolerance_ * norm_b)
            {
                if (verbose_) {
                    std::cout << "CG converged in " << iter + 1 << " iterations" << std::endl;
                }
                return true;
            }

            // beta = rho_new / rho
            RealType beta = rho_new / rho;

            // p = z + beta * p
            // First update owned portion: p_owned = z + beta * p_owned
            thrust::transform(thrust::device_pointer_cast(z.data()),
                             thrust::device_pointer_cast(z.data() + m),
                             thrust::device_pointer_cast(p.data()),
                             thrust::device_pointer_cast(p.data()),
                             [beta] __host__ __device__ (RealType z_val, RealType p_val) {
                                 return z_val + beta * p_val;
                             });
            // Ghost portion of p will be updated in next halo exchange

            rho = rho_new;
        }

        if (verbose_) {
            std::cout << "CG did not converge in " << maxIter_ << " iterations" << std::endl;
        }
        return false;
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }
    void setHaloExchangeCallback(std::function<void(Vector&)> callback) {
        haloExchangeCallback_ = callback;
    }

private:

    // Helper functions using cuBLAS
    RealType dot(const Vector& a, const Vector& b) {
        RealType result;
        if constexpr (std::is_same_v<RealType, float>) {
            cublasSdot(cublasHandle_, a.size(), a.data(), 1, b.data(), 1, &result);
        } else {
            cublasDdot(cublasHandle_, a.size(), a.data(), 1, b.data(), 1, &result);
        }
        return result;
    }

    void axpy(RealType alpha, const Vector& x, Vector& y) {
        if constexpr (std::is_same_v<RealType, float>) {
            cublasSaxpy(cublasHandle_, x.size(), &alpha, x.data(), 1, y.data(), 1);
        } else {
            cublasDaxpy(cublasHandle_, x.size(), &alpha, x.data(), 1, y.data(), 1);
        }
    }

    int maxIter_;
    RealType tolerance_;
    bool verbose_;
    std::function<void(Vector&)> haloExchangeCallback_;
    std::shared_ptr<Precond> preconditioner_;

    cublasHandle_t cublasHandle_;
    cusparseHandle_t cusparseHandle_;
};

} // namespace fem
} // namespace mars