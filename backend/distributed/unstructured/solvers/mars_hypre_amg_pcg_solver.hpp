#pragma once

// Hypre PCG preconditioned by one BoomerAMG V-cycle, for SPD matrices that stay constant in
// time: the matrix and the AMG hierarchy are built once, and each solve only runs PCG.
//
// The matrix is built on the GPU from a few pieces: COO assembly of rows this rank owns,
// products and transposes of distributed matrices, and identity rows. The Galerkin product
// P^T A P turns a matrix over each rank's local node slots into the matrix over the DOFs.

#include "mars_hypre_pcg_solver.hpp"
#include <_hypre_parcsr_mv.h>
#include <thrust/device_vector.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>

namespace mars
{
namespace fem
{

inline bool hypreEnvFlag(const char* name, bool fallback)
{
    const char* option = std::getenv(name);
    return option ? std::string(option) == "1" : fallback;
}

// Hypre's own GPU kernels by default:
// - SpMV: the vendor path returned wrong products intermittently with Hypre 2.33 on the SIMPLE duct.
// - SpGEMM: cuSPARSE's default SpGEMM stopped the setup with "insufficient resources" on 16 and
//   64 GPUs (2M nodes per GPU); 4 GPUs with the same load passed.
// - GPU-aware MPI: without it Hypre copies every halo exchange through host memory, on every
//   matvec of every AMG level. MARS already passes device buffers to MPI in its own halos.
// MARS_HYPRE_SPMV_VENDOR=1, MARS_HYPRE_SPGEMM_VENDOR=1 and MARS_HYPRE_GPU_AWARE=0 switch back.
inline void hypreInitializeOnce()
{
    static HypreInitGuard guard;
    static bool configured = false;
    if (configured) return;
    configured = true;
    const bool spmvVendor   = hypreEnvFlag("MARS_HYPRE_SPMV_VENDOR", false);
    const bool spgemmVendor = hypreEnvFlag("MARS_HYPRE_SPGEMM_VENDOR", false);
    const bool gpuAwareMpi  = hypreEnvFlag("MARS_HYPRE_GPU_AWARE", true);
    const int buildGpuAware = hypre_GetGpuAwareMPI();
    HYPRE_SetSpMVUseVendor(spmvVendor ? 1 : 0);
    HYPRE_SetSpGemmUseVendor(spgemmVendor ? 1 : 0);
    HYPRE_SetGpuAwareMPI(gpuAwareMpi ? 1 : 0);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        std::printf("Hypre: SpMV %s, SpGEMM %s, GPU-aware MPI %s (build default %s)\n",
                    spmvVendor ? "vendor" : "hypre", spgemmVendor ? "vendor" : "hypre", gpuAwareMpi ? "on" : "off",
                    buildGpuAware ? "on" : "off");
}

// Collective: stops every rank if one rank fails.
inline void hypreCheck(MPI_Comm comm, bool localOk, const char* message)
{
    int bad = localOk ? 0 : 1, anyBad = 0;
    MPI_Allreduce(&bad, &anyBad, 1, MPI_INT, MPI_MAX, comm);
    if (!anyBad) return;
    int rank = 0;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) std::fprintf(stderr, "ERROR: Hypre: %s\n", message);
    MPI_Abort(comm, 1);
}

// Owns a distributed Hypre matrix.
class HypreMatrix
{
public:
    HypreMatrix() = default;
    explicit HypreMatrix(hypre_ParCSRMatrix* A)
        : A_(A)
    {
    }
    HypreMatrix(HypreMatrix&& other) noexcept
        : A_(std::exchange(other.A_, nullptr))
    {
    }
    HypreMatrix& operator=(HypreMatrix&& other) noexcept
    {
        std::swap(A_, other.A_);
        return *this;
    }
    HypreMatrix(const HypreMatrix&)            = delete;
    HypreMatrix& operator=(const HypreMatrix&) = delete;
    ~HypreMatrix()
    {
        if (A_) hypre_ParCSRMatrixDestroy(A_);
    }

    hypre_ParCSRMatrix* get() const { return A_; }
    hypre_ParCSRMatrix* release() { return std::exchange(A_, nullptr); }
    HYPRE_BigInt firstRow() const { return hypre_ParCSRMatrixFirstRowIndex(A_); }
    HYPRE_BigInt endRow() const { return hypre_ParCSRMatrixLastRowIndex(A_) + 1; }

private:
    hypre_ParCSRMatrix* A_ = nullptr;
};

struct HypreCoo
{
    const HYPRE_BigInt* rows;
    const HYPRE_BigInt* cols;
    const HYPRE_Complex* values;
    HYPRE_Int size;
};

// This rank's rows are [rowStart, rowEnd) and its share of the columns [colStart, colEnd).
// Every triplet must lie in an owned row: Hypre 2.33 crashed sending device triplets to
// other ranks. Duplicate entries are summed. All arrays are on the device.
inline HypreMatrix hypreAssemble(MPI_Comm comm, HYPRE_BigInt rowStart, HYPRE_BigInt rowEnd, HYPRE_BigInt colStart,
                                 HYPRE_BigInt colEnd, HypreCoo coo)
{
    hypreInitializeOnce();
    HYPRE_IJMatrix ij = nullptr;
    HYPRE_IJMatrixCreate(comm, rowStart, rowEnd - 1, colStart, colEnd - 1, &ij);
    HYPRE_IJMatrixSetObjectType(ij, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize_v2(ij, HYPRE_MEMORY_DEVICE);
    // ncols == nullptr: one entry per triplet.
    if (coo.size > 0) HYPRE_IJMatrixAddToValues(ij, coo.size, nullptr, coo.rows, coo.cols, coo.values);
    HYPRE_IJMatrixAssemble(ij);
    hypre_ParCSRMatrix* A = nullptr;
    HYPRE_IJMatrixGetObject(ij, reinterpret_cast<void**>(&A));
    // The IJ object owns A; keep a copy that outlives it.
    HypreMatrix copy(hypre_ParCSRMatrixClone(A, 1));
    HYPRE_IJMatrixDestroy(ij);
    hypreCheck(comm, HYPRE_GetError() == 0 && copy.get() != nullptr, "COO assembly failed");
    return copy;
}

inline HypreMatrix hypreTranspose(const HypreMatrix& A)
{
    hypre_ParCSRMatrix* T = nullptr;
    hypre_ParCSRMatrixTranspose(A.get(), &T, 1);
    return HypreMatrix(T);
}

inline HypreMatrix hypreMultiply(const HypreMatrix& A, const HypreMatrix& B)
{
    return HypreMatrix(hypre_ParCSRMatMat(A.get(), B.get()));
}

// P^T A P
inline HypreMatrix hypreGalerkin(const HypreMatrix& P, const HypreMatrix& A)
{
    return hypreMultiply(hypreTranspose(P), hypreMultiply(A, P));
}

// A + I on the given rows, global ids owned by this rank. The rows must be empty in A.
inline HypreMatrix hypreAddIdentityRows(MPI_Comm comm, HypreMatrix A, const HYPRE_BigInt* rows, HYPRE_Int count)
{
    long long total = count;
    MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_LONG_LONG, MPI_SUM, comm);
    if (total == 0) return A;
    thrust::device_vector<HYPRE_Complex> ones(count, HYPRE_Complex(1));
    HypreMatrix I = hypreAssemble(comm, A.firstRow(), A.endRow(), A.firstRow(), A.endRow(),
                                  HypreCoo{rows, rows, thrust::raw_pointer_cast(ones.data()), count});
    hypre_ParCSRMatrix* sum = nullptr;
    hypre_ParCSRMatrixAdd(1.0, A.get(), 1.0, I.get(), &sum);
    hypreCheck(comm, sum != nullptr, "adding identity rows failed");
    return HypreMatrix(sum);
}

class HypreAmgPcgSolver
{
public:
    // envPrefix selects the environment overrides of the AMG settings, e.g. MARS_PAMG_STRONG.
    explicit HypreAmgPcgSolver(std::string envPrefix)
        : envPrefix_(std::move(envPrefix))
    {
    }
    HypreAmgPcgSolver(const HypreAmgPcgSolver&)            = delete;
    HypreAmgPcgSolver& operator=(const HypreAmgPcgSolver&) = delete;
    ~HypreAmgPcgSolver() { destroy(); }

    // Takes the matrix and builds the AMG hierarchy. This rank's rows of A are its unknowns.
    void setup(MPI_Comm comm, HypreMatrix A, double tolerance, int maxIter)
    {
        hypreInitializeOnce();
        destroy();
        comm_            = comm;
        partitioning_[0] = A.firstRow();
        partitioning_[1] = A.endRow();
        rows_            = HYPRE_Int(partitioning_[1] - partitioning_[0]);
        A_               = std::move(A);
        trace("diagonal first");
#if defined(HYPRE_USING_GPU)
        // Relaxation reads the diagonal as the first entry of each row.
        hypre_CSRMatrixMoveDiagFirstDevice(hypre_ParCSRMatrixDiag(A_.get()));
#endif
        HYPRE_BigInt globalRows = hypre_ParCSRMatrixGlobalNumRows(A_.get());
        HYPRE_ParVectorCreate(comm_, globalRows, partitioning_, &f_);
        HYPRE_ParVectorCreate(comm_, globalRows, partitioning_, &u_);
        HYPRE_ParVectorInitialize(f_);
        HYPRE_ParVectorInitialize(u_);

        HYPRE_BoomerAMGCreate(&amg_);
        HYPRE_BoomerAMGSetCoarsenType(amg_, envInt("COARSEN", 8)); // PMIS
        HYPRE_BoomerAMGSetInterpType(amg_, envInt("INTERP", 6));   // extended+i
        HYPRE_BoomerAMGSetPMaxElmts(amg_, envInt("PMAX", 4));
        HYPRE_BoomerAMGSetAggNumLevels(amg_, envInt("AGG", 0));
        HYPRE_BoomerAMGSetStrongThreshold(amg_, envDouble("STRONG", 0.25));
        HYPRE_BoomerAMGSetRelaxType(amg_, 18); // l1-Jacobi: symmetric and parallel, so the V-cycle suits CG
        HYPRE_BoomerAMGSetRelaxOrder(amg_, 0);
        HYPRE_BoomerAMGSetNumSweeps(amg_, 1);
        // Direct elimination would break on a singular coarse operator; smooth there instead.
        HYPRE_BoomerAMGSetCycleRelaxType(amg_, 18, 3);
        HYPRE_BoomerAMGSetCycleNumSweeps(amg_, envInt("COARSE_SWEEPS", 4), 3);
        HYPRE_BoomerAMGSetKeepTranspose(amg_, 1);
        HYPRE_BoomerAMGSetMaxLevels(amg_, 25);
        HYPRE_BoomerAMGSetTol(amg_, 0.0);
        HYPRE_BoomerAMGSetMaxIter(amg_, 1);
        HYPRE_BoomerAMGSetPrintLevel(amg_, envInt("PRINT", 0));

        HYPRE_ParCSRPCGCreate(comm_, &pcg_);
        HYPRE_PCGSetTol(pcg_, tolerance);
        HYPRE_PCGSetAbsoluteTol(pcg_, 0.0);
        HYPRE_PCGSetTwoNorm(pcg_, 1); // ||r|| <= tol ||b||
        HYPRE_PCGSetMaxIter(pcg_, maxIter);
        HYPRE_PCGSetLogging(pcg_, 1);
        HYPRE_PCGSetPrintLevel(pcg_, 0);
        HYPRE_PCGSetPrecond(pcg_, reinterpret_cast<HYPRE_PtrToSolverFcn>(HYPRE_BoomerAMGSolve),
                            reinterpret_cast<HYPRE_PtrToSolverFcn>(HYPRE_BoomerAMGSetup), amg_);
        trace("PCG + BoomerAMG setup");
        HYPRE_Int error = HYPRE_ParCSRPCGSetup(pcg_, matrix(), f_, u_);
        hypreCheck(comm_, error == 0 && HYPRE_GetError() == 0, "PCG/BoomerAMG setup failed");
        trace("setup done");
    }

    // f and u hold this rank's unknowns on the device. u is the initial guess when
    // useInitialGuess is set, and is overwritten. Returns the PCG iterations, or -2 if the
    // relative tolerance ||b - A u|| <= tol ||b|| was not reached.
    int solve(const HYPRE_Complex* f, HYPRE_Complex* u, bool useInitialGuess = false)
    {
        copyIn(f, f_);
        HYPRE_Real rhsNorm2 = 0;
        HYPRE_ParVectorInnerProd(f_, f_, &rhsNorm2);
        if (rhsNorm2 == 0)
        {
            // Hypre returns x = 0 for a zero right-hand side without flagging convergence.
            HYPRE_ParVectorSetConstantValues(u_, 0.0);
            copyOut(u_, u);
            lastRelativeResidual_ = 0;
            lastIterations_       = 0;
            return 0;
        }
        if (useInitialGuess) copyIn(u, u_);
        else HYPRE_ParVectorSetConstantValues(u_, 0.0);
        trace("solve");
        HYPRE_Int error = HYPRE_ParCSRPCGSolve(pcg_, matrix(), f_, u_);
        trace("solve done");
        HYPRE_Int iterations = 0, converged = 0;
        HYPRE_PCGGetNumIterations(pcg_, &iterations);
        HYPRE_PCGGetConverged(pcg_, &converged);
        HYPRE_PCGGetFinalRelativeResidualNorm(pcg_, &lastRelativeResidual_);
        lastIterations_ = int(iterations);
        hypreCheck(comm_, (error & ~HYPRE_ERROR_CONV) == 0, "PCG solve failed");
        HYPRE_ClearAllErrors();
        copyOut(u_, u);
        return converged ? int(iterations) : -2;
    }

    // y = A x on this rank's unknowns.
    void apply(const HYPRE_Complex* x, HYPRE_Complex* y)
    {
        copyIn(x, f_);
        HYPRE_ParCSRMatrixMatvec(1.0, matrix(), f_, 0.0, u_);
        copyOut(u_, y);
    }

    double lastRelativeResidual() const { return lastRelativeResidual_; }
    int lastIterations() const { return lastIterations_; }

private:
    HYPRE_ParCSRMatrix matrix() const { return reinterpret_cast<HYPRE_ParCSRMatrix>(A_.get()); }

    void copyIn(const HYPRE_Complex* src, HYPRE_ParVector dst) const
    {
        auto* data = hypre_VectorData(hypre_ParVectorLocalVector(reinterpret_cast<hypre_ParVector*>(dst)));
        if (rows_ > 0) cudaMemcpy(data, src, rows_ * sizeof(HYPRE_Complex), cudaMemcpyDeviceToDevice);
    }

    void copyOut(HYPRE_ParVector src, HYPRE_Complex* dst) const
    {
        auto* data = hypre_VectorData(hypre_ParVectorLocalVector(reinterpret_cast<hypre_ParVector*>(src)));
        if (rows_ > 0) cudaMemcpy(dst, data, rows_ * sizeof(HYPRE_Complex), cudaMemcpyDeviceToDevice);
    }

    // MARS_AMG_TRACE=1: each setup and solve phase per rank, flushed, to locate a crash or hang.
    void trace(const char* phase) const
    {
        static const bool enabled = std::getenv("MARS_AMG_TRACE") != nullptr;
        if (!enabled) return;
        int rank = 0;
        MPI_Comm_rank(comm_, &rank);
        std::fprintf(stderr, "[amg-trace] %s rank %d: %s\n", envPrefix_.c_str(), rank, phase);
        std::fflush(stderr);
    }

    int envInt(const char* name, int fallback) const
    {
        const char* v = std::getenv((envPrefix_ + "_" + name).c_str());
        return v && *v ? std::atoi(v) : fallback;
    }

    double envDouble(const char* name, double fallback) const
    {
        const char* v = std::getenv((envPrefix_ + "_" + name).c_str());
        return v && *v ? std::atof(v) : fallback;
    }

    void destroy()
    {
        if (pcg_) HYPRE_ParCSRPCGDestroy(pcg_);
        if (amg_) HYPRE_BoomerAMGDestroy(amg_);
        if (f_) HYPRE_ParVectorDestroy(f_);
        if (u_) HYPRE_ParVectorDestroy(u_);
        pcg_ = amg_ = nullptr;
        f_ = u_ = nullptr;
        A_      = HypreMatrix();
    }

    std::string envPrefix_;
    MPI_Comm comm_                = MPI_COMM_WORLD;
    HYPRE_Int rows_               = 0;
    HYPRE_BigInt partitioning_[2] = {0, 0};
    HypreMatrix A_;
    HYPRE_ParVector f_           = nullptr;
    HYPRE_ParVector u_           = nullptr;
    HYPRE_Solver amg_            = nullptr;
    HYPRE_Solver pcg_            = nullptr;
    double lastRelativeResidual_ = 0;
    int lastIterations_          = 0;
};

} // namespace fem
} // namespace mars
