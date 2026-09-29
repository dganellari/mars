#pragma once

// Hypre PCG preconditioned by one BoomerAMG V-cycle, for SPD matrices that stay constant in
// time: the matrix and the AMG hierarchy are built once, and each solve only runs PCG.
//
// Two ways to build the matrix, both on the GPU:
//  - setupProjection: the exact projection operator A = (D S) D^T with S diagonal. The caller
//    gives D and D S as global COO triplets; rows may belong to other ranks and duplicates are
//    summed, so each rank adds the faces of its own elements. Pinned rows become identity rows.
//  - setupCsr: owned rows of a local CSR matrix whose columns are local DOFs (owned and ghost).

#include "mars_hypre_pcg_solver.hpp"
#include <_hypre_parcsr_mv.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/fill.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>

namespace mars
{
namespace fem
{

template<typename ValueType>
__global__ void zeroCsrRowsKernel(const HYPRE_BigInt* pinnedRows, HYPRE_Int numPinned, HYPRE_BigInt rowStart,
                                  const HYPRE_Int* diagI, ValueType* diagData, const HYPRE_Int* offdI,
                                  ValueType* offdData)
{
    HYPRE_Int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= numPinned) return;
    HYPRE_Int row = HYPRE_Int(pinnedRows[k] - rowStart);
    for (HYPRE_Int j = diagI[row]; j < diagI[row + 1]; ++j) diagData[j] = 0;
    if (offdI)
    {
        for (HYPRE_Int j = offdI[row]; j < offdI[row + 1]; ++j) offdData[j] = 0;
    }
}

// Global triplets of a local CSR block. A column without a global DOF gets row -1.
template<typename IndexType>
__global__ void localCsrToGlobalCooKernel(const IndexType* rowPtr, const IndexType* colInd, IndexType numRows,
                                          const HYPRE_BigInt* localToGlobal, size_t numLocalCols,
                                          HYPRE_BigInt rowStart, HYPRE_BigInt* rows, HYPRE_BigInt* cols)
{
    IndexType i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numRows) return;
    for (IndexType j = rowPtr[i]; j < rowPtr[i + 1]; ++j)
    {
        IndexType c    = colInd[j];
        HYPRE_BigInt g = (c >= 0 && size_t(c) < numLocalCols) ? localToGlobal[c] : HYPRE_BigInt(-1);
        rows[j]        = g < 0 ? HYPRE_BigInt(-1) : rowStart + i;
        cols[j]        = g;
    }
}

class HypreAmgPcgSolver
{
public:
    struct Coo
    {
        const HYPRE_BigInt* rows;
        const HYPRE_BigInt* cols;
        const HYPRE_Complex* values;
        HYPRE_Int size;
    };

    // envPrefix selects the environment overrides of the AMG settings, e.g. MARS_PAMG_STRONG.
    explicit HypreAmgPcgSolver(std::string envPrefix)
        : envPrefix_(std::move(envPrefix))
    {
    }
    HypreAmgPcgSolver(const HypreAmgPcgSolver&)            = delete;
    HypreAmgPcgSolver& operator=(const HypreAmgPcgSolver&) = delete;
    ~HypreAmgPcgSolver() { destroy(); }

    // Owned rows are [rowStart, rowEnd); the columns of D are [components * rowStart,
    // components * rowEnd). All arrays are on the device.
    void setupProjection(MPI_Comm comm, HYPRE_BigInt rowStart, HYPRE_BigInt rowEnd, int components,
                         Coo divergence, Coo scaledDivergence, const HYPRE_BigInt* pinnedRows,
                         HYPRE_Int numPinned, double tolerance, int maxIter)
    {
        begin(comm, rowStart, rowEnd);
        HYPRE_BigInt colStart = components * rowStart, colEnd = components * rowEnd;
        trace("assemble D");
        HYPRE_IJMatrix ijD = assemble(rowStart, rowEnd, colStart, colEnd, divergence);
        trace("assemble D S");
        HYPRE_IJMatrix ijDS = assemble(rowStart, rowEnd, colStart, colEnd, scaledDivergence);
        hypre_ParCSRMatrix *D = nullptr, *DS = nullptr, *DT = nullptr;
        HYPRE_IJMatrixGetObject(ijD, reinterpret_cast<void**>(&D));
        HYPRE_IJMatrixGetObject(ijDS, reinterpret_cast<void**>(&DS));
        trace("transpose D");
        hypre_ParCSRMatrixTranspose(D, &DT, 1);
        trace("product (D S) D^T");
        hypre_ParCSRMatrix* A0 = hypre_ParCSRMatMat(DS, DT);
        HYPRE_IJMatrixDestroy(ijD);
        HYPRE_IJMatrixDestroy(ijDS);
        hypre_ParCSRMatrixDestroy(DT);
        check(A0 != nullptr, "sparse product D S D^T failed");

        // A pinned row may still hold entries: zero it before adding the identity row.
        if (numPinned > 0)
        {
            hypre_CSRMatrix* diag = hypre_ParCSRMatrixDiag(A0);
            hypre_CSRMatrix* offd = hypre_ParCSRMatrixOffd(A0);
            bool haveOffd         = hypre_CSRMatrixNumCols(offd) > 0 && hypre_CSRMatrixI(offd);
            zeroCsrRowsKernel<HYPRE_Complex><<<(numPinned + 255) / 256, 256>>>(
                pinnedRows, numPinned, rowStart, hypre_CSRMatrixI(diag), hypre_CSRMatrixData(diag),
                haveOffd ? hypre_CSRMatrixI(offd) : nullptr, haveOffd ? hypre_CSRMatrixData(offd) : nullptr);
        }
        // Collective: every rank checks, also ranks without pinned rows.
        check(cudaGetLastError() == cudaSuccess, "pinned row elimination failed");
        thrust::device_vector<HYPRE_Complex> ones(numPinned, HYPRE_Complex(1));
        trace("pinned identity rows");
        HYPRE_IJMatrix ijI =
            assemble(rowStart, rowEnd, rowStart, rowEnd,
                     Coo{pinnedRows, pinnedRows, thrust::raw_pointer_cast(ones.data()), numPinned});
        hypre_ParCSRMatrix* I = nullptr;
        HYPRE_IJMatrixGetObject(ijI, reinterpret_cast<void**>(&I));
        trace("add identity rows");
        hypre_ParCSRMatrixAdd(1.0, A0, 1.0, I, &A_);
        HYPRE_IJMatrixDestroy(ijI);
        hypre_ParCSRMatrixDestroy(A0);
        check(A_ != nullptr, "adding the pinned identity rows failed");
        finish(tolerance, maxIter);
    }

    // Owned rows [rowStart, rowEnd) of a local CSR matrix; column c is the global DOF
    // localToGlobal[c]. Every entry must have a global column, so A is the caller's matrix.
    template<typename IndexType>
    void setupCsr(MPI_Comm comm, HYPRE_BigInt rowStart, HYPRE_BigInt rowEnd, const IndexType* rowPtr,
                  const IndexType* colInd, const HYPRE_Complex* values, const HYPRE_BigInt* localToGlobal,
                  size_t numLocalCols, double tolerance, int maxIter)
    {
        begin(comm, rowStart, rowEnd);
        trace("assemble CSR");
        IndexType nnz = 0;
        if (rows_ > 0) cudaMemcpy(&nnz, rowPtr + rows_, sizeof(IndexType), cudaMemcpyDeviceToHost);
        thrust::device_vector<HYPRE_BigInt> rows(nnz), cols(nnz);
        if (rows_ > 0)
        {
            localCsrToGlobalCooKernel<IndexType><<<(rows_ + 255) / 256, 256>>>(
                rowPtr, colInd, IndexType(rows_), localToGlobal, numLocalCols, rowStart,
                thrust::raw_pointer_cast(rows.data()), thrust::raw_pointer_cast(cols.data()));
        }
        check(cudaDeviceSynchronize() == cudaSuccess, "CSR conversion failed");
        check(thrust::count(rows.begin(), rows.end(), HYPRE_BigInt(-1)) == 0,
              "a matrix column has no global DOF");
        // The IJ matrix owns A_ and is kept until destroy().
        ij_ = assemble(rowStart, rowEnd, rowStart, rowEnd,
                       Coo{thrust::raw_pointer_cast(rows.data()), thrust::raw_pointer_cast(cols.data()), values,
                           HYPRE_Int(nnz)});
        HYPRE_IJMatrixGetObject(ij_, reinterpret_cast<void**>(&A_));
        check(A_ != nullptr, "matrix assembly failed");
        finish(tolerance, maxIter);
    }

    // f and u hold the owned DOF values on the device. u is the initial guess when
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
            return 0;
        }
        if (useInitialGuess) copyIn(u, u_);
        else HYPRE_ParVectorSetConstantValues(u_, 0.0);
        trace("solve");
        HYPRE_Int error = HYPRE_ParCSRPCGSolve(pcg_, reinterpret_cast<HYPRE_ParCSRMatrix>(A_), f_, u_);
        trace("solve done");
        HYPRE_Int iterations = 0, converged = 0;
        HYPRE_PCGGetNumIterations(pcg_, &iterations);
        HYPRE_PCGGetConverged(pcg_, &converged);
        HYPRE_PCGGetFinalRelativeResidualNorm(pcg_, &lastRelativeResidual_);
        check((error & ~HYPRE_ERROR_CONV) == 0, "PCG solve failed");
        HYPRE_ClearAllErrors();
        copyOut(u_, u);
        return converged ? int(iterations) : -2;
    }

    // y = A x on the owned DOFs.
    void apply(const HYPRE_Complex* x, HYPRE_Complex* y)
    {
        trace("matvec");
        copyIn(x, f_);
        HYPRE_ParCSRMatrixMatvec(1.0, reinterpret_cast<HYPRE_ParCSRMatrix>(A_), f_, 0.0, u_);
        copyOut(u_, y);
        trace("matvec done");
    }

    double lastRelativeResidual() const { return lastRelativeResidual_; }

private:
    void begin(MPI_Comm comm, HYPRE_BigInt rowStart, HYPRE_BigInt rowEnd)
    {
        static HypreInitGuard hypreInit;
        (void)hypreInit;
        destroy();
        comm_ = comm;
        rows_ = HYPRE_Int(rowEnd - rowStart);
        selectSpMV();
        partitioning_[0] = rowStart;
        partitioning_[1] = rowEnd;
    }

    void finish(double tolerance, int maxIter)
    {
        trace("diagonal first");
#if defined(HYPRE_USING_GPU)
        // Relaxation reads the diagonal as the first entry of each row.
        hypre_CSRMatrixMoveDiagFirstDevice(hypre_ParCSRMatrixDiag(A_));
#endif
        HYPRE_BigInt globalRows = hypre_ParCSRMatrixGlobalNumRows(A_);
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
        HYPRE_PCGSetTwoNorm(pcg_, 1); // ||r|| <= tol ||b||, as in the MARS CG
        HYPRE_PCGSetMaxIter(pcg_, maxIter);
        HYPRE_PCGSetLogging(pcg_, 1);
        HYPRE_PCGSetPrintLevel(pcg_, 0);
        HYPRE_PCGSetPrecond(pcg_, reinterpret_cast<HYPRE_PtrToSolverFcn>(HYPRE_BoomerAMGSolve),
                            reinterpret_cast<HYPRE_PtrToSolverFcn>(HYPRE_BoomerAMGSetup), amg_);
        trace("PCG + BoomerAMG setup");
        HYPRE_Int error = HYPRE_ParCSRPCGSetup(pcg_, reinterpret_cast<HYPRE_ParCSRMatrix>(A_), f_, u_);
        check(error == 0 && HYPRE_GetError() == 0, "PCG/BoomerAMG setup failed");
        trace("setup done");
    }

    HYPRE_IJMatrix assemble(HYPRE_BigInt rowStart, HYPRE_BigInt rowEnd, HYPRE_BigInt colStart,
                            HYPRE_BigInt colEnd, Coo coo)
    {
        HYPRE_IJMatrix ij = nullptr;
        HYPRE_IJMatrixCreate(comm_, rowStart, rowEnd - 1, colStart, colEnd - 1, &ij);
        HYPRE_IJMatrixSetObjectType(ij, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize_v2(ij, HYPRE_MEMORY_DEVICE);
        // ncols == nullptr: one entry per triplet. Off-rank rows are sent to their owners.
        if (coo.size > 0)
        {
            HYPRE_IJMatrixAddToValues(ij, coo.size, nullptr, coo.rows, coo.cols, coo.values);
        }
        HYPRE_IJMatrixAssemble(ij);
        check(HYPRE_GetError() == 0, "IJ assembly failed");
        return ij;
    }

    // The vendor SpMV path returned wrong products intermittently with Hypre 2.33 on the
    // SIMPLE duct; Hypre's own GPU SpMV did not. MARS_HYPRE_SPMV_VENDOR=1 opts back in.
    static void selectSpMV()
    {
        const char* option = std::getenv("MARS_HYPRE_SPMV_VENDOR");
        HYPRE_SetSpMVUseVendor(option && std::string(option) == "1" ? 1 : 0);
    }

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

    // MARS_PAMG_TRACE=1: each setup phase per rank, flushed, to locate a crash.
    void trace(const char* phase) const
    {
        static const bool enabled = std::getenv("MARS_PAMG_TRACE") != nullptr;
        if (!enabled) return;
        int rank = 0;
        MPI_Comm_rank(comm_, &rank);
        std::fprintf(stderr, "[amg-trace] %s rank %d: %s\n", envPrefix_.c_str(), rank, phase);
        std::fflush(stderr);
    }

    void check(bool localOk, const char* message) const
    {
        int bad = localOk ? 0 : 1, anyBad = 0;
        MPI_Allreduce(&bad, &anyBad, 1, MPI_INT, MPI_MAX, comm_);
        if (anyBad)
        {
            int rank = 0;
            MPI_Comm_rank(comm_, &rank);
            if (rank == 0) std::cerr << "ERROR: " << envPrefix_ << " solver: " << message << '\n';
            MPI_Abort(comm_, 1);
        }
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
        if (ij_) HYPRE_IJMatrixDestroy(ij_);
        else if (A_) hypre_ParCSRMatrixDestroy(A_);
        pcg_ = amg_ = nullptr;
        f_ = u_ = nullptr;
        A_       = nullptr;
        ij_      = nullptr;
    }

    std::string envPrefix_;
    MPI_Comm comm_                = MPI_COMM_WORLD;
    HYPRE_Int rows_               = 0;
    HYPRE_BigInt partitioning_[2] = {0, 0};
    HYPRE_IJMatrix ij_            = nullptr;
    hypre_ParCSRMatrix* A_        = nullptr;
    HYPRE_ParVector f_            = nullptr;
    HYPRE_ParVector u_            = nullptr;
    HYPRE_Solver amg_             = nullptr;
    HYPRE_Solver pcg_             = nullptr;
    double lastRelativeResidual_  = 0;
};

} // namespace fem
} // namespace mars
