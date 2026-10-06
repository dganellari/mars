#pragma once

// Hypre PCG preconditioned by one BoomerAMG V-cycle, for SPD matrices that stay constant in
// time: the matrix and the AMG hierarchy are built once, and each solve only runs PCG.
//
// hypreFromEntries builds the distributed matrix on the GPU from (local row, global column,
// value) entries; the Navier-Stokes solver forms them with DofSpace::restrictMatrix.

#include "mars_hypre_pcg_solver.hpp"
#include <_hypre_parcsr_mv.h>
#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/partition.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/transform.h>
#include <thrust/unique.h>
#include <cstdio>
#include <cstdlib>
#include <limits>
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

inline size_t hypreEnvMiB(const char* name, size_t fallback)
{
    const char* option = std::getenv(name);
    return (option ? size_t(std::strtoull(option, nullptr, 10)) : fallback) << 20;
}

// Hypre's own GPU kernels by default:
// - SpMV: the vendor path returned wrong products intermittently with Hypre 2.33 on the SIMPLE duct.
// - SpGEMM: cuSPARSE's default SpGEMM stopped the setup with "insufficient resources" on 16 and
//   64 GPUs (2M nodes per GPU); 4 GPUs with the same load passed.
// - GPU-aware MPI: without it Hypre copies every halo exchange through host memory, on every
//   matvec of every AMG level. MARS already passes device buffers to MPI in its own halos.
// MARS_HYPRE_SPMV_VENDOR=1, MARS_HYPRE_SPGEMM_VENDOR=1 and MARS_HYPRE_GPU_AWARE=0 switch back.
// MARS_HYPRE_POOL_MAX_MIB and MARS_HYPRE_POOL_CACHE_MIB set the device pool limits below.
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

    // Hypre's CUB device pool (when Hypre is built with it) by default rounds every allocation
    // up to the next power of 8 and keeps every freed block, so setup memory never goes back to
    // the GPU; TGV with 8M nodes per GPU ran out of memory on 8 GPUs. Pool only blocks up to
    // 128 MiB, rounded to a power of 2, and cache at most 2 GiB; larger arrays are allocated
    // exactly and freed at once. This must come before Hypre's first device allocation.
    const size_t poolMax   = hypreEnvMiB("MARS_HYPRE_POOL_MAX_MIB", 128);
    const size_t poolCache = hypreEnvMiB("MARS_HYPRE_POOL_CACHE_MIB", 2048);
    int maxBin             = 0;
    while ((size_t(2) << maxBin) <= poolMax)
        ++maxBin;
    HYPRE_SetGPUMemoryPoolSize(2, 1, maxBin, poolCache);
    char memory[96];
#if defined(HYPRE_USING_UMPIRE_DEVICE)
    std::snprintf(memory, sizeof(memory), "Umpire pool");
#elif defined(HYPRE_USING_CUDA) && defined(HYPRE_USING_DEVICE_POOL)
    std::snprintf(memory, sizeof(memory), "CUB pool, blocks up to %zu MiB, at most %zu MiB cached",
                  (size_t(1) << maxBin) >> 20, poolCache >> 20);
#elif defined(HYPRE_USING_CUDA) && defined(HYPRE_USING_DEVICE_MALLOC_ASYNC)
    std::snprintf(memory, sizeof(memory), "cudaMallocAsync");
#else
    std::snprintf(memory, sizeof(memory), "no pool");
#endif

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
        std::printf("Hypre: SpMV %s, SpGEMM %s, GPU-aware MPI %s (build default %s), device memory: %s\n",
                    spmvVendor ? "vendor" : "hypre", spgemmVendor ? "vendor" : "hypre", gpuAwareMpi ? "on" : "off",
                    buildGpuAware ? "on" : "off", memory);
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

struct HypreRowColumnKey
{
    __host__ __device__ unsigned long long operator()(int row, unsigned long long column) const
    {
        return ((unsigned long long)row << 32) | column;
    }
};

struct HypreKeyRow
{
    __host__ __device__ int operator()(unsigned long long key) const { return int(key >> 32); }
};

// Compressed column ids [lo, hi) are this rank's own columns.
struct HypreKeyInBlock
{
    unsigned long long lo, hi;
    __host__ __device__ bool operator()(const thrust::tuple<unsigned long long, HYPRE_Complex>& entry) const
    {
        unsigned long long c = thrust::get<0>(entry) & 0xffffffffull;
        return c >= lo && c < hi;
    }
};

struct HypreDiagColumn
{
    const long long* columns;
    long long first;
    __host__ __device__ HYPRE_Int operator()(unsigned long long key) const
    {
        return HYPRE_Int(columns[key & 0xffffffffull] - first);
    }
};

struct HypreOffdColumn
{
    unsigned long long lo, width;
    __host__ __device__ HYPRE_Int operator()(unsigned long long key) const
    {
        unsigned long long c = key & 0xffffffffull;
        return HYPRE_Int(c < lo ? c : c - width);
    }
};

// The matrix over this rank's rows [firstRow, firstRow + rows) from its entries (row = local
// row, col = global column, duplicates summed), in the layout Hypre keeps: this rank's own
// columns (diag, local ids) and the others (offd, ids into a sorted column map). Two radix
// sorts of 64-bit keys do the work: the columns compressed to this rank's sorted column set,
// then (row << 32 | column). The input vectors are consumed.
inline HypreMatrix hypreFromEntries(MPI_Comm comm, long long firstRow, HYPRE_Int rows, long long globalRows,
                                    thrust::device_vector<int>& row, thrust::device_vector<long long>& col,
                                    thrust::device_vector<HYPRE_Complex>& val)
{
    hypreInitializeOnce();
    // Rows are numbered in 64 bits here; a Hypre built with 32-bit global indices would wrap
    // them silently above 2^31 - 1.
    hypreCheck(comm, globalRows <= (long long)std::numeric_limits<HYPRE_BigInt>::max(),
               "more global rows than this Hypre build can index; build Hypre with --enable-mixedint");
    const size_t n = row.size();
    thrust::device_vector<long long> columns(col);
    thrust::sort(thrust::device, columns.begin(), columns.end());
    columns.resize(thrust::unique(thrust::device, columns.begin(), columns.end()) - columns.begin());

    thrust::device_vector<unsigned long long> keys(n);
    thrust::lower_bound(thrust::device, columns.begin(), columns.end(), col.begin(), col.end(), keys.begin());
    thrust::transform(thrust::device, row.begin(), row.end(), keys.begin(), keys.begin(), HypreRowColumnKey{});
    thrust::device_vector<int>().swap(row);
    thrust::device_vector<long long>().swap(col);
    thrust::sort_by_key(thrust::device, keys.begin(), keys.end(), val.begin());
    thrust::device_vector<unsigned long long> entryKey(n);
    thrust::device_vector<HYPRE_Complex> entryVal(n);
    const size_t nnz = thrust::reduce_by_key(thrust::device, keys.begin(), keys.end(), val.begin(), entryKey.begin(),
                                             entryVal.begin())
                           .first -
                       entryKey.begin();
    thrust::device_vector<unsigned long long>().swap(keys);
    thrust::device_vector<HYPRE_Complex>().swap(val);

    // Own columns first, rows kept in order within each part.
    const unsigned long long lo =
        thrust::lower_bound(thrust::device, columns.begin(), columns.end(), (long long)firstRow) - columns.begin();
    const unsigned long long hi =
        thrust::lower_bound(thrust::device, columns.begin(), columns.end(), (long long)(firstRow + rows)) -
        columns.begin();
    auto entries = thrust::make_zip_iterator(thrust::make_tuple(entryKey.begin(), entryVal.begin()));
    const size_t nnzDiag =
        thrust::stable_partition(thrust::device, entries, entries + nnz, HypreKeyInBlock{lo, hi}) - entries;
    const size_t nnzOffd = nnz - nnzDiag;
    const HYPRE_Int colsOffd = HYPRE_Int(columns.size() - (hi - lo));

    HYPRE_BigInt starts[2] = {HYPRE_BigInt(firstRow), HYPRE_BigInt(firstRow + rows)};
    hypre_ParCSRMatrix* A  = hypre_ParCSRMatrixCreate(comm, HYPRE_BigInt(globalRows), HYPRE_BigInt(globalRows), starts,
                                                      starts, colsOffd, HYPRE_Int(nnzDiag), HYPRE_Int(nnzOffd));
    hypre_ParCSRMatrixInitialize_v2(A, HYPRE_MEMORY_DEVICE);
    hypre_CSRMatrix* diag = hypre_ParCSRMatrixDiag(A);
    hypre_CSRMatrix* offd = hypre_ParCSRMatrixOffd(A);
    auto rowOf            = thrust::make_transform_iterator(entryKey.begin(), HypreKeyRow{});
    auto rowIds           = thrust::counting_iterator<int>(0);
    thrust::lower_bound(thrust::device, rowOf, rowOf + nnzDiag, rowIds, rowIds + rows + 1,
                        thrust::device_pointer_cast(hypre_CSRMatrixI(diag)));
    thrust::lower_bound(thrust::device, rowOf + nnzDiag, rowOf + nnz, rowIds, rowIds + rows + 1,
                        thrust::device_pointer_cast(hypre_CSRMatrixI(offd)));
    const long long* columnIds = thrust::raw_pointer_cast(columns.data());
    if (nnzDiag > 0)
    {
        thrust::transform(thrust::device, entryKey.begin(), entryKey.begin() + nnzDiag,
                          thrust::device_pointer_cast(hypre_CSRMatrixJ(diag)), HypreDiagColumn{columnIds, firstRow});
        thrust::copy(thrust::device, entryVal.begin(), entryVal.begin() + nnzDiag,
                     thrust::device_pointer_cast(hypre_CSRMatrixData(diag)));
    }
    if (nnzOffd > 0)
    {
        thrust::transform(thrust::device, entryKey.begin() + nnzDiag, entryKey.begin() + nnz,
                          thrust::device_pointer_cast(hypre_CSRMatrixJ(offd)), HypreOffdColumn{lo, hi - lo});
        thrust::copy(thrust::device, entryVal.begin() + nnzDiag, entryVal.begin() + nnz,
                     thrust::device_pointer_cast(hypre_CSRMatrixData(offd)));
    }
    if (colsOffd > 0)
    {
        thrust::device_vector<HYPRE_BigInt> colMap(colsOffd);
        auto next = thrust::copy(thrust::device, columns.begin(), columns.begin() + lo, colMap.begin());
        thrust::copy(thrust::device, columns.begin() + hi, columns.end(), next);
        cudaMemcpy(hypre_ParCSRMatrixColMapOffd(A), thrust::raw_pointer_cast(colMap.data()),
                   colsOffd * sizeof(HYPRE_BigInt), cudaMemcpyDeviceToHost);
    }
    hypre_ParCSRMatrixSetNumNonzeros(A);
    hypreCheck(comm, HYPRE_GetError() == 0, "matrix setup from entries failed");
    return HypreMatrix(A);
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
        HYPRE_BoomerAMGSetAggInterpType(amg_, envInt("AGG_INTERP", 6)); // 5-8 have GPU kernels; 6: 2-stage ext+i
        HYPRE_BoomerAMGSetStrongThreshold(amg_, envDouble("STRONG", 0.25));
        HYPRE_BoomerAMGSetRelaxType(amg_, 18); // l1-Jacobi: symmetric and parallel, so the V-cycle suits CG
        HYPRE_BoomerAMGSetRelaxOrder(amg_, 0);
        HYPRE_BoomerAMGSetNumSweeps(amg_, 1);
        // Direct elimination would break on a singular coarse operator; smooth there instead.
        HYPRE_BoomerAMGSetCycleRelaxType(amg_, 18, 3);
        HYPRE_BoomerAMGSetCycleNumSweeps(amg_, envInt("COARSE_SWEEPS", 4), 3);
        HYPRE_BoomerAMGSetKeepTranspose(amg_, 1);
        HYPRE_BoomerAMGSetMaxLevels(amg_, envInt("MAX_LEVELS", 25)); // 1: l1-Jacobi PCG, no hierarchy
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
