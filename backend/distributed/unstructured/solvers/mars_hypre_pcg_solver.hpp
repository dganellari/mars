#pragma once

#include "../fem/mars_sparse_matrix.hpp"
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_utilities.h>
#include <mpi.h>
#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/scan.h>
#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <thrust/logical.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/functional.h>
#include <cmath>
#include <iostream>
#include <vector>

namespace mars {
namespace fem {

// One-time Hypre init/finalize for GPU execution.
// Constructed lazily on first solve() so MPI_Init is guaranteed to have run.
struct HypreInitGuard {
    HypreInitGuard() {
#if defined(HYPRE_RELEASE_NUMBER) && HYPRE_RELEASE_NUMBER >= 22000
        HYPRE_Initialize();
#else
        HYPRE_Init();
#endif
        HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE);
        HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE);
    }
    ~HypreInitGuard() {
#if defined(HYPRE_RELEASE_NUMBER) && HYPRE_RELEASE_NUMBER >= 22000
        HYPRE_Finalize();
#else
        HYPRE_Finalize();
#endif
    }
};

// Pass 1: per-row count valid entries + record whether the diagonal is present.
// One thread per owned row. Walks that row's CSR slice, bounds-checks each
// local column, looks up the global column, and counts valid entries.
template<typename IndexType, typename HBigInt>
__global__ void countValidPerRowKernel(
    const IndexType* d_rowOffsets,
    const IndexType* d_colIndicesLocal,
    const HBigInt*   d_localToGlobalDof,
    size_t           numLocalCols,
    HBigInt          ilower,
    IndexType        m,
    int*             d_perRowCount,
    int*             d_hasDiagonal)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= m) return;

    IndexType rowStart = d_rowOffsets[i];
    IndexType rowEnd   = d_rowOffsets[i + 1];
    HBigInt globalRow  = ilower + static_cast<HBigInt>(i);

    int valid = 0;
    int hasDiag = 0;
    for (IndexType j = rowStart; j < rowEnd; ++j) {
        IndexType localCol = d_colIndicesLocal[j];
        // signed-safe bound check
        if (localCol < 0) continue;
        size_t lc = static_cast<size_t>(localCol);
        if (lc >= numLocalCols) continue;
        HBigInt globalCol = d_localToGlobalDof[lc];
        if (globalCol < 0) continue;
        ++valid;
        if (globalCol == globalRow) hasDiag = 1;
    }
    d_perRowCount[i] = valid;
    d_hasDiagonal[i] = hasDiag;
}

// Pass 2: compact filtered (globalCol, value) entries into per-row contiguous slots.
// One thread per owned row. Uses the exclusive-scanned d_outOffsets as the
// write base, so the output is dense (CSR-shaped, but with width = perRowCount[i]).
template<typename IndexType, typename RealType, typename HBigInt, typename HReal>
__global__ void compactGlobalCsrKernel(
    const IndexType* d_rowOffsets,
    const IndexType* d_colIndicesLocal,
    const RealType*  d_valuesLocal,
    const HBigInt*   d_localToGlobalDof,
    size_t           numLocalCols,
    const int*       d_outOffsets,
    IndexType        m,
    HBigInt*         d_colsGlobalCompact,
    HReal*           d_valsCompact)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= m) return;

    IndexType rowStart = d_rowOffsets[i];
    IndexType rowEnd   = d_rowOffsets[i + 1];
    int writeBase      = d_outOffsets[i];
    int written        = 0;

    for (IndexType j = rowStart; j < rowEnd; ++j) {
        IndexType localCol = d_colIndicesLocal[j];
        if (localCol < 0) continue;
        size_t lc = static_cast<size_t>(localCol);
        if (lc >= numLocalCols) continue;
        HBigInt globalCol = d_localToGlobalDof[lc];
        if (globalCol < 0) continue;
        d_colsGlobalCompact[writeBase + written] = globalCol;
        d_valsCompact      [writeBase + written] = static_cast<HReal>(d_valuesLocal[j]);
        ++written;
    }
}

// Fill global row-index array [ilower, ilower+m) on device.
template<typename HBigInt>
__global__ void fillGlobalRowIndicesKernel(HBigInt* d_rows, HBigInt ilower, int m)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) d_rows[i] = ilower + static_cast<HBigInt>(i);
}

// NaN/Inf scan reduction predicate.
template<typename T>
struct IsNonFinite {
    __host__ __device__ bool operator()(T v) const {
#if defined(__CUDA_ARCH__)
        return !isfinite(static_cast<double>(v));
#else
        return !std::isfinite(static_cast<double>(v));
#endif
    }
};

// GPU-resident Hypre PCG + BoomerAMG solver.
// All matrix/RHS/solution data stays on device end-to-end.
template<typename RealType, typename IndexType, typename AcceleratorTag>
class HyprePCGSolver {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    enum PrecondType { BOOMERAMG, JACOBI };

    HyprePCGSolver(MPI_Comm comm = MPI_COMM_WORLD, int maxIter = 1000, RealType tolerance = 1e-6, PrecondType precondType = BOOMERAMG)
        : comm_(comm), maxIter_(maxIter), tolerance_(tolerance),
          solver_(nullptr), precond_(nullptr), A_hypre_(nullptr),
          b_hypre_(nullptr), x_hypre_(nullptr), verbose_(true), precondType_(precondType) {
    }

    ~HyprePCGSolver() {
        destroy();
    }

    bool solve(const Matrix& A, const Vector& b, Vector& x,
               IndexType globalDofStart = 0, IndexType globalDofEnd = 0) {
        return solve<IndexType>(A, b, x, globalDofStart, globalDofEnd,
                               globalDofStart, globalDofEnd, {});
    }

    // GPU-direct overload: caller provides a device-resident map of local DOF
    // index -> global DOF index. No H2D copy at all. Used by the NS/Poisson/
    // advection drivers where the map is built by a device kernel and never
    // touches the host.
    bool solve(const Matrix& A, const Vector& b, Vector& x,
               IndexType globalDofStart, IndexType globalDofEnd,
               IndexType globalColStart, IndexType globalColEnd,
               const thrust::device_vector<HYPRE_BigInt>& d_localToGlobalDof) {
        static HypreInitGuard g_hypreInit;
        (void)g_hypreInit;

        int rank;
        MPI_Comm_rank(comm_, &rank);
        if (verbose_) {
            std::cout << "Rank " << rank << ": Entering Hypre solve (device-map) with globalDofRange ["
                      << globalDofStart << ", " << globalDofEnd << ")" << std::endl;
        }

        destroy();
        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());
        if (globalDofEnd == 0) globalDofEnd = globalDofStart + m;
        globalDofStart_ = globalDofStart;
        globalDofEnd_   = globalDofEnd;

        // Take a device-side copy of the caller's map so the wrapper has stable
        // storage across the whole solve sequence (matrix setup, then solve).
        d_localToGlobalDof_ = d_localToGlobalDof;

        return solveImpl(A, b, x, globalColStart, globalColEnd);
    }

    // Legacy host-vector overload (compat shim). Uploads to device once and
    // delegates. New code should call the device-vector overload above.
    // globalDofStart/End: this rank's owned global row range [start, end)
    // globalColStart/End: total global interior column range
    // localToGlobalDof:   host vector mapping local col index -> global DOF
    template<typename KeyType>
    bool solve(const Matrix& A, const Vector& b, Vector& x,
               IndexType globalDofStart, IndexType globalDofEnd,
               IndexType globalColStart, IndexType globalColEnd,
               const std::vector<KeyType>& localToGlobalDof) {
        // Initialize Hypre exactly once per process, lazily, after MPI is up.
        static HypreInitGuard g_hypreInit;
        (void)g_hypreInit;

        int rank;
        MPI_Comm_rank(comm_, &rank);
        if (verbose_) {
            std::cout << "Rank " << rank << ": Entering Hypre solve with globalDofRange ["
                      << globalDofStart << ", " << globalDofEnd << ")" << std::endl;
        }

        destroy();

        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());

        if (globalDofEnd == 0) {
            globalDofEnd = globalDofStart + m;
        }
        globalDofStart_ = globalDofStart;
        globalDofEnd_   = globalDofEnd;

        // Upload local->global DOF map to device once.
        // For the no-mapping fallback case (identity = ilower + localCol), build
        // the identity map on the fly so the kernel always sees a valid array.
        d_localToGlobalDof_.clear();
        d_localToGlobalDof_.shrink_to_fit();
        if (!localToGlobalDof.empty()) {
            std::vector<HYPRE_BigInt> hLocalToGlobal(localToGlobalDof.size());
            for (size_t i = 0; i < localToGlobalDof.size(); ++i) {
                hLocalToGlobal[i] = static_cast<HYPRE_BigInt>(localToGlobalDof[i]);
            }
            d_localToGlobalDof_ = hLocalToGlobal;  // H2D
        } else {
            // Identity: localCol -> ilower + localCol over the owned range.
            std::vector<HYPRE_BigInt> identity(m);
            for (HYPRE_Int i = 0; i < m; ++i) {
                identity[i] = static_cast<HYPRE_BigInt>(globalDofStart_) + i;
            }
            d_localToGlobalDof_ = identity;
        }

        return solveImpl(A, b, x, globalColStart, globalColEnd);
    }

    void setVerbose(bool verbose) { verbose_ = verbose; }

    // The helpers below are conceptually private — they take internal Hypre
    // state and aren't meant to be called from outside — but nvcc rejects
    // extended __device__ lambdas inside private member functions, so we
    // leave them under the public access label. No external caller invokes
    // them.

    // Shared by both overloads. Assumes d_localToGlobalDof_, globalDofStart_/End_
    // are populated, destroy() has been called, and globalDofEnd defaulting has
    // been done.
    bool solveImpl(const Matrix& A, const Vector& b, Vector& x,
                   IndexType globalColStart, IndexType globalColEnd) {
        int rank;
        MPI_Comm_rank(comm_, &rank);
        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());

        setupHypreMatrix(A, globalColStart, globalColEnd);

        // RHS + initial guess: pure device path.
        HYPRE_BigInt ilower = static_cast<HYPRE_BigInt>(globalDofStart_);
        HYPRE_BigInt iupper = static_cast<HYPRE_BigInt>(globalDofEnd_ - 1);

        HYPRE_IJVectorCreate(comm_, ilower, iupper, &b_hypre_);
        HYPRE_IJVectorCreate(comm_, ilower, iupper, &x_hypre_);
        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(b_hypre_);
        HYPRE_IJVectorInitialize(x_hypre_);

        // Device-side global row indices [ilower, ilower+m).
        thrust::device_vector<HYPRE_BigInt> d_rowGlobal(m);
        {
            const int bs = 256;
            const int gs = (m + bs - 1) / bs;
            fillGlobalRowIndicesKernel<HYPRE_BigInt><<<gs, bs>>>(
                thrust::raw_pointer_cast(d_rowGlobal.data()), ilower, m);
        }

        // RHS NaN/Inf summary + min/max/sum, all on device.
        validateVector(b.data(), m, rank, "RHS");

        // Hypre RHS values must be HYPRE_Real; cast on device if RealType != HYPRE_Real.
        const HYPRE_Real* d_b_hypre = nullptr;
        thrust::device_vector<HYPRE_Real> d_b_cast;
        if constexpr (std::is_same_v<RealType, HYPRE_Real>) {
            d_b_hypre = reinterpret_cast<const HYPRE_Real*>(b.data());
        } else {
            d_b_cast.resize(m);
            thrust::copy(thrust::device_pointer_cast(b.data()),
                         thrust::device_pointer_cast(b.data() + m),
                         d_b_cast.begin());
            d_b_hypre = thrust::raw_pointer_cast(d_b_cast.data());
        }

        HYPRE_IJVectorSetValues(b_hypre_, m,
                                thrust::raw_pointer_cast(d_rowGlobal.data()),
                                d_b_hypre);

        // Initial guess.
        const HYPRE_Real* d_x_hypre = nullptr;
        thrust::device_vector<HYPRE_Real> d_x_cast;
        if constexpr (std::is_same_v<RealType, HYPRE_Real>) {
            d_x_hypre = reinterpret_cast<const HYPRE_Real*>(x.data());
        } else {
            d_x_cast.resize(m);
            thrust::copy(thrust::device_pointer_cast(x.data()),
                         thrust::device_pointer_cast(x.data() + m),
                         d_x_cast.begin());
            d_x_hypre = thrust::raw_pointer_cast(d_x_cast.data());
        }
        HYPRE_IJVectorSetValues(x_hypre_, m,
                                thrust::raw_pointer_cast(d_rowGlobal.data()),
                                d_x_hypre);

        HYPRE_IJVectorAssemble(b_hypre_);
        HYPRE_IJVectorAssemble(x_hypre_);

        if (verbose_) std::cout << "Rank " << rank << ": Vectors assembled, getting ParVector objects..." << std::endl;
        HYPRE_IJVectorGetObject(b_hypre_, (void**)&par_b_);
        HYPRE_IJVectorGetObject(x_hypre_, (void**)&par_x_);

        if (!par_b_ || !par_x_) {
            std::cerr << "Rank " << rank << ": Failed to get Hypre ParVector objects" << std::endl;
            return false;
        }

        if (verbose_ && rank == 0) std::cout << "Creating preconditioner..." << std::endl;

        if (precondType_ == BOOMERAMG) {
            if (verbose_ && rank == 0) std::cout << "Using BoomerAMG preconditioner (GPU)" << std::endl;
            HYPRE_BoomerAMGCreate(&precond_);
            HYPRE_BoomerAMGSetPrintLevel(precond_, 0);
        } else if (precondType_ == JACOBI) {
            if (verbose_ && rank == 0) std::cout << "Using Jacobi preconditioner" << std::endl;
            precond_ = (HYPRE_Solver) parcsr_A_;
        }

        if (verbose_ && rank == 0) std::cout << "Creating PCG solver..." << std::endl;
        HYPRE_ParCSRPCGCreate(comm_, &solver_);
        if (!solver_) {
            std::cerr << "Failed to create Hypre PCG solver" << std::endl;
            return false;
        }
        HYPRE_PCGSetMaxIter(solver_, maxIter_);
        HYPRE_PCGSetTol(solver_, tolerance_);
        HYPRE_PCGSetPrintLevel(solver_, verbose_ ? 2 : 0);

        if (precondType_ == BOOMERAMG && precond_) {
            if (verbose_ && rank == 0) std::cout << "Setting BoomerAMG preconditioner..." << std::endl;
            HYPRE_PCGSetPrecond(solver_,
                               (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSolve,
                               (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSetup,
                               precond_);
        } else if (precondType_ == JACOBI) {
            if (verbose_ && rank == 0) std::cout << "Setting Jacobi (diagonal) preconditioner..." << std::endl;
            HYPRE_PCGSetPrecond(solver_,
                               (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_ParCSRDiagScale,
                               (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_ParCSRDiagScaleSetup,
                               (HYPRE_Solver) parcsr_A_);
        }

        if (verbose_ && rank == 0) std::cout << "Setting up PCG solver..." << std::endl;
        MPI_Barrier(comm_);
        HYPRE_ParCSRPCGSetup(solver_, parcsr_A_, par_b_, par_x_);
        if (verbose_ && rank == 0) std::cout << "PCG setup complete, starting solve..." << std::endl;
        MPI_Barrier(comm_);
        HYPRE_ParCSRPCGSolve(solver_, parcsr_A_, par_b_, par_x_);
        if (verbose_ && rank == 0) std::cout << "PCG solve complete." << std::endl;

        // Read solution directly into device storage.
        if constexpr (std::is_same_v<RealType, HYPRE_Real>) {
            HYPRE_IJVectorGetValues(x_hypre_, m,
                                    thrust::raw_pointer_cast(d_rowGlobal.data()),
                                    reinterpret_cast<HYPRE_Real*>(x.data()));
        } else {
            d_x_cast.resize(m);
            HYPRE_IJVectorGetValues(x_hypre_, m,
                                    thrust::raw_pointer_cast(d_rowGlobal.data()),
                                    thrust::raw_pointer_cast(d_x_cast.data()));
            thrust::copy(d_x_cast.begin(), d_x_cast.end(),
                         thrust::device_pointer_cast(x.data()));
        }

        int    num_iterations = 0;
        double final_res_norm = 0.0;
        HYPRE_PCGGetNumIterations(solver_, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver_, &final_res_norm);

        if (verbose_) {
            std::cout << "Hypre PCG converged in " << num_iterations
                      << " iterations, final residual: " << final_res_norm << std::endl;
        }

        return (final_res_norm < tolerance_);
    }

    // Build A_hypre_ entirely on the device:
    //   1) count valid entries per row + diagonal presence
    //   2) exclusive-scan to per-row write offsets (total = filtered nnz)
    //   3) compact (globalCol, value) into dense per-row slots
    //   4) one HYPRE_IJMatrixSetValues call with all device pointers
    void setupHypreMatrix(const Matrix& A,
                          IndexType globalColStart, IndexType globalColEnd) {
        int rank;
        MPI_Comm_rank(comm_, &rank);
        HYPRE_Int m       = static_cast<HYPRE_Int>(A.numRows());
        HYPRE_BigInt ilower = static_cast<HYPRE_BigInt>(globalDofStart_);
        HYPRE_BigInt iupper = static_cast<HYPRE_BigInt>(globalDofEnd_ - 1);

        // Column partitioning must match row partitioning so x is distributed like b.
        // Hypre handles off-processor (ghost) columns in ParCSR automatically.
        (void)globalColStart;
        (void)globalColEnd;

        HYPRE_IJMatrixCreate(comm_, ilower, iupper, ilower, iupper, &A_hypre_);
        HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(A_hypre_);

        // Quick NaN/Inf scan on raw values (single allreduce-style reduction).
        bool hasNaN = thrust::any_of(thrust::device_pointer_cast(A.valuesPtr()),
                                     thrust::device_pointer_cast(A.valuesPtr() + A.nnz()),
                                     IsNonFinite<RealType>());
        if (hasNaN) {
            std::cerr << "Rank " << rank << ": ERROR - matrix contains NaN/Inf values!" << std::endl;
            return;
        }

        const size_t numLocalCols = d_localToGlobalDof_.size();
        const HYPRE_Int nnzLocal  = static_cast<HYPRE_Int>(A.nnz());

        // Pass 1: per-row valid count + diagonal flag.
        thrust::device_vector<int> d_perRowCount(m, 0);
        thrust::device_vector<int> d_hasDiagonal(m, 0);
        {
            const int bs = 256;
            const int gs = (m + bs - 1) / bs;
            countValidPerRowKernel<IndexType, HYPRE_BigInt><<<gs, bs>>>(
                A.rowOffsetsPtr(),
                A.colIndicesPtr(),
                thrust::raw_pointer_cast(d_localToGlobalDof_.data()),
                numLocalCols,
                ilower,
                m,
                thrust::raw_pointer_cast(d_perRowCount.data()),
                thrust::raw_pointer_cast(d_hasDiagonal.data()));
        }

        // Exclusive scan -> per-row write offsets; total filtered nnz = last + count[m-1].
        thrust::device_vector<int> d_outOffsets(m + 1, 0);
        thrust::exclusive_scan(d_perRowCount.begin(), d_perRowCount.end(),
                               d_outOffsets.begin());
        // total = scan[m-1] + count[m-1]
        int totalFiltered = 0;
        int lastOffset = 0, lastCount = 0;
        if (m > 0) {
            cudaMemcpy(&lastOffset,
                       thrust::raw_pointer_cast(d_outOffsets.data()) + (m - 1),
                       sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(&lastCount,
                       thrust::raw_pointer_cast(d_perRowCount.data()) + (m - 1),
                       sizeof(int), cudaMemcpyDeviceToHost);
        }
        totalFiltered = lastOffset + lastCount;

        // Pass 2: compact.
        thrust::device_vector<HYPRE_BigInt> d_colsGlobalCompact(totalFiltered);
        thrust::device_vector<HYPRE_Real>   d_valsCompact(totalFiltered);
        {
            const int bs = 256;
            const int gs = (m + bs - 1) / bs;
            compactGlobalCsrKernel<IndexType, RealType, HYPRE_BigInt, HYPRE_Real><<<gs, bs>>>(
                A.rowOffsetsPtr(),
                A.colIndicesPtr(),
                A.valuesPtr(),
                thrust::raw_pointer_cast(d_localToGlobalDof_.data()),
                numLocalCols,
                thrust::raw_pointer_cast(d_outOffsets.data()),
                m,
                thrust::raw_pointer_cast(d_colsGlobalCompact.data()),
                thrust::raw_pointer_cast(d_valsCompact.data()));
        }

        // Device-side row indices for the SetValues call.
        thrust::device_vector<HYPRE_BigInt> d_rows(m);
        {
            const int bs = 256;
            const int gs = (m + bs - 1) / bs;
            fillGlobalRowIndicesKernel<HYPRE_BigInt><<<gs, bs>>>(
                thrust::raw_pointer_cast(d_rows.data()), ilower, m);
        }

        // Validation summaries (small reductions, single int copy each).
        int emptyRows = static_cast<int>(
            thrust::count(d_perRowCount.begin(), d_perRowCount.end(), 0));
        int rowsMissingDiagonal = thrust::transform_reduce(
            thrust::make_counting_iterator(0),
            thrust::make_counting_iterator(static_cast<int>(m)),
            [d_perRowCount_ptr = thrust::raw_pointer_cast(d_perRowCount.data()),
             d_hasDiagonal_ptr = thrust::raw_pointer_cast(d_hasDiagonal.data())]
            __device__ (int i) -> int {
                return (d_perRowCount_ptr[i] > 0 && d_hasDiagonal_ptr[i] == 0) ? 1 : 0;
            },
            0,
            thrust::plus<int>());

        // Global col range (cheap reduction; only on non-zero rows would be exact,
        // but min/max over a possibly-empty compact array still gives a useful summary).
        HYPRE_BigInt minGlobalCol = 0, maxGlobalCol = -1;
        if (totalFiltered > 0) {
            auto mm = thrust::minmax_element(d_colsGlobalCompact.begin(),
                                             d_colsGlobalCompact.end());
            minGlobalCol = *mm.first;
            maxGlobalCol = *mm.second;
        }

        if (verbose_ || emptyRows > 0 || rowsMissingDiagonal > 0) {
            std::cout << "Rank " << rank << ": Matrix filtering: " << nnzLocal
                      << " original entries -> " << totalFiltered << " valid entries ("
                      << (nnzLocal - totalFiltered) << " filtered)" << std::endl;
            std::cout << "Rank " << rank << ": Global column range: [" << minGlobalCol
                      << ", " << maxGlobalCol << "], global row range: [" << ilower
                      << ", " << iupper << "]" << std::endl;
        }
        if (emptyRows > 0) {
            std::cout << "Rank " << rank << ": WARNING - " << emptyRows
                      << " rows became empty after filtering!" << std::endl;
        }
        if (rowsMissingDiagonal > 0) {
            std::cout << "Rank " << rank << ": WARNING - " << rowsMissingDiagonal
                      << " rows missing diagonal after filtering!" << std::endl;
        }

        // One device-pointer SetValues call: HYPRE_MEMORY_DEVICE has been set globally,
        // so Hypre reads ncols/rows/cols/values directly from device memory.
        HYPRE_IJMatrixSetValues(A_hypre_, m,
                                thrust::raw_pointer_cast(d_perRowCount.data()),
                                thrust::raw_pointer_cast(d_rows.data()),
                                thrust::raw_pointer_cast(d_colsGlobalCompact.data()),
                                thrust::raw_pointer_cast(d_valsCompact.data()));

        HYPRE_IJMatrixAssemble(A_hypre_);
        if (verbose_) std::cout << "Rank " << rank << ": Matrix assembled, getting ParCSR object..." << std::endl;
        HYPRE_IJMatrixGetObject(A_hypre_, (void**)&parcsr_A_);

        if (!parcsr_A_) {
            std::cerr << "Rank " << rank << ": Failed to get Hypre ParCSR matrix object" << std::endl;
            return;
        }
    }

    // Single thrust reduction: NaN/Inf check + min/max/sum, printed once.
    void validateVector(const RealType* d_ptr, HYPRE_Int n, int rank, const char* label) {
        if (n <= 0) return;
        bool hasNaN = thrust::any_of(thrust::device_pointer_cast(d_ptr),
                                     thrust::device_pointer_cast(d_ptr + n),
                                     IsNonFinite<RealType>());
        if (verbose_ || hasNaN) {
            // For min/max/sum we restrict to finite values by masking. Cheap pass.
            double sumVal = thrust::transform_reduce(
                thrust::device_pointer_cast(d_ptr),
                thrust::device_pointer_cast(d_ptr + n),
                [] __device__ (RealType v) -> double {
                    return isfinite(static_cast<double>(v)) ? static_cast<double>(v) : 0.0;
                },
                0.0,
                thrust::plus<double>());
            double minVal = thrust::transform_reduce(
                thrust::device_pointer_cast(d_ptr),
                thrust::device_pointer_cast(d_ptr + n),
                [] __device__ (RealType v) -> double {
                    return isfinite(static_cast<double>(v)) ? static_cast<double>(v) : 1e100;
                },
                1e100,
                thrust::minimum<double>());
            double maxVal = thrust::transform_reduce(
                thrust::device_pointer_cast(d_ptr),
                thrust::device_pointer_cast(d_ptr + n),
                [] __device__ (RealType v) -> double {
                    return isfinite(static_cast<double>(v)) ? static_cast<double>(v) : -1e100;
                },
                -1e100,
                thrust::maximum<double>());
            std::cout << "Rank " << rank << ": Before Hypre - " << label << " sum=" << sumVal
                      << ", range=[" << minVal << ", " << maxVal << "], NaN/Inf="
                      << (hasNaN ? "YES" : "no") << " /" << n << std::endl;
            if (hasNaN) {
                std::cout << "Rank " << rank << " ERROR: " << label << " contains NaN/Inf!" << std::endl;
            }
        }
    }

    void destroy() {
        if (solver_) {
            HYPRE_ParCSRPCGDestroy(solver_);
            solver_ = nullptr;
        }
        if (precond_ && precondType_ == BOOMERAMG) {
            HYPRE_BoomerAMGDestroy(precond_);
            precond_ = nullptr;
        }
        if (A_hypre_) {
            HYPRE_IJMatrixDestroy(A_hypre_);
            A_hypre_ = nullptr;
        }
        if (b_hypre_) {
            HYPRE_IJVectorDestroy(b_hypre_);
            b_hypre_ = nullptr;
        }
        if (x_hypre_) {
            HYPRE_IJVectorDestroy(x_hypre_);
            x_hypre_ = nullptr;
        }
    }

private:
    MPI_Comm comm_;
    int maxIter_;
    RealType tolerance_;
    bool verbose_;
    PrecondType precondType_;

    HYPRE_Solver solver_;
    HYPRE_Solver precond_;

    HYPRE_IJMatrix A_hypre_;
    HYPRE_ParCSRMatrix parcsr_A_;

    HYPRE_IJVector b_hypre_;
    HYPRE_IJVector x_hypre_;
    HYPRE_ParVector par_b_;
    HYPRE_ParVector par_x_;

    IndexType globalDofStart_;
    IndexType globalDofEnd_;

    // Device-resident local->global DOF map, uploaded once per solve.
    thrust::device_vector<HYPRE_BigInt> d_localToGlobalDof_;
};

} // namespace fem
} // namespace mars
