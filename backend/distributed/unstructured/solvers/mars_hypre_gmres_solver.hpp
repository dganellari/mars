#pragma once

// All shared utilities (HypreInitGuard, countValidPerRowKernel,
// compactGlobalCsrKernel, fillGlobalRowIndicesKernel, IsNonFinite, the
// HYPRE_/thrust/mpi headers, the namespace declarations) come from the PCG
// header. Including it here means we don't redefine those symbols.
#include "mars_hypre_pcg_solver.hpp"
#if !defined(HYPRE_RELEASE_NUMBER) || HYPRE_RELEASE_NUMBER < 30000
// Older releases export ParVectorAxpy but declare it only in this header.
#include <_hypre_parcsr_mv.h>
#endif
#include "mars_solver_profile.hpp"
#include <limits>
#include <stdexcept>
#include <thrust/gather.h>

namespace mars {
namespace fem {

// Keep the source slot as well as the global column: later solves only gather values.
template<typename IndexType>
__global__ void pack_hypre_graph_kernel(const IndexType* offsets, const IndexType* columns,
    const HYPRE_BigInt* map, size_t map_size, const int* packed_offsets, int rows,
    HYPRE_BigInt* packed_columns, IndexType* source_slots)
{
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= rows) return;
    int slot = packed_offsets[row];
    for (IndexType j = offsets[row]; j < offsets[row + 1]; ++j) {
        const IndexType col = columns[j];
        if (col < 0 || static_cast<size_t>(col) >= map_size || map[col] < 0) continue;
        packed_columns[slot] = map[col];
        source_slots[slot++] = j;
    }
}

// GPU-resident Hypre GMRES + BoomerAMG solver.
// API-identical drop-in for HyprePCGSolver. Uses Hypre's restarted GMRES
// (HYPRE_ParCSRGMRES) instead of PCG, which removes the strict SPD
// requirement: works for SPSD pressure-Poisson with constant null space
// (D M^-1 D^T) where PCG returns error 1.
//
// Memory cost: GMRES(k) restart length stores ~k+5 Krylov vectors of size
// numOwnedDofs. Default k=30 -> 35 vectors of 8 bytes per rank for FP64.
// At wing scale (15M DOFs, 4 ranks): 35*3.75M*8 = ~1 GB per rank. Fine.
template<typename RealType, typename IndexType, typename AcceleratorTag>
class HypreGMRESSolver {
public:
    using Matrix = SparseMatrix<IndexType, RealType, AcceleratorTag>;
    using Vector = typename mars::VectorSelector<RealType, AcceleratorTag>::type;

    enum PrecondType { BOOMERAMG, JACOBI };

    void set_profile(SolverProfile* profile) { profile_ = profile; }

    // Use ||b-Ax||/||b|| instead of unit-dependent x/b magnitude heuristics.
    // A zero RHS uses an absolute residual. The caller must fix any nullspace.
    void enable_true_residual_check() {
        true_residual_check_ = true;
        mixed_residual_check_ = false;
    }
    // Accept the caller's ||b-Ax|| <= atol + rtol*||b|| without changing
    // the tighter Krylov target. Hypre can stop before reaching that target.
    void enable_true_residual_check(double absolute_tolerance, double relative_tolerance) {
        if (!std::isfinite(absolute_tolerance) || !std::isfinite(relative_tolerance)
            || absolute_tolerance < 0 || relative_tolerance < 0
            || (absolute_tolerance == 0 && relative_tolerance == 0))
            throw std::runtime_error("invalid true residual tolerances");
        true_residual_check_ = true;
        mixed_residual_check_ = true;
        residual_absolute_tolerance_ = absolute_tolerance;
        residual_relative_tolerance_ = relative_tolerance;
    }

    // The caller must invalidate before changing matrix/map contents in place.
    // The outlet owner bounds reuse to a single frozen physical step.
    void enable_reuse(bool amg_cycle = false) {
        destroy();
        fixed_graph_updates_ = false;
        reuse_enabled_ = true;
        amg_cycle_ = amg_cycle;
    }
    // Device-map overload and nonempty owned partitions only. Invalidate before
    // changing CSR/map contents in place. Hypre may still rebuild its inner GPU CSR.
    // All numeric values and the AMG hierarchy are refreshed on every solve.
    void enable_fixed_graph_updates() {
        destroy();
        reuse_enabled_ = false;
        amg_cycle_ = false;
        fixed_graph_updates_ = true;
    }
    void invalidate_setup() { destroy(); }
    int get_setup_count() const { return setup_count_; }
    int get_graph_build_count() const { return graph_build_count_; }
    int get_numeric_update_count() const { return numeric_update_count_; }

    struct SolveTiming {
        double prepare_seconds = 0, packing_seconds = 0, setup_seconds = 0;
        double solve_seconds = 0, finish_seconds = 0;
    };
    // Local wall API times without synchronization; packing is part of prepare.
    void enable_timing(bool enabled = true) { timing_enabled_ = enabled; }
    const SolveTiming& get_last_timing() const { return last_timing_; }

    HypreGMRESSolver(MPI_Comm comm = MPI_COMM_WORLD, int maxIter = 1000, RealType tolerance = 1e-6,
                     PrecondType precondType = BOOMERAMG, int kDim = 30)
        : comm_(comm), maxIter_(maxIter), tolerance_(tolerance),
          solver_(nullptr), precond_(nullptr), A_hypre_(nullptr),
          b_hypre_(nullptr), x_hypre_(nullptr), verbose_(true), precondType_(precondType),
          kDim_(kDim) {
        // Plain GMRES vs FlexGMRES. A 1-V-cycle BoomerAMG preconditioner IS a
        // fixed linear operator, so plain GMRES (HYPRE_GMRESSetPrecond) is valid
        // AND is the most battle-tested GPU path. FlexGMRES (varying precond) was
        // tried as the default but its precond-apply appeared to NO-OP on our GPU
        // Hypre build (AMG returned x=0 with a denormal ~1e-315 residual), while
        // plain GMRES+AMG is the standard. Default OFF (plain GMRES). FlexGMRES
        // only matters if the precond varies (>1 AMG sweep w/ inner tol); enable
        // with MARS_HYPRE_FLEXGMRES=1.
        useFlexGmres_ = false;
        const char* fev = std::getenv("MARS_HYPRE_FLEXGMRES");
        if (fev) useFlexGmres_ = (std::string(fev) != "0");
    }

    ~HypreGMRESSolver() {
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
        last_timing_ = {};
        prepare_start_ = wall_stamp();
        if (profile_) profile_start_ = profile_->stamp();
        static HypreInitGuard g_hypreInit;
        (void)g_hypreInit;

        int rank;
        MPI_Comm_rank(comm_, &rank);
        if (verbose_) {
            std::cout << "Rank " << rank << ": Entering Hypre solve (device-map) with globalDofRange ["
                      << globalDofStart << ", " << globalDofEnd << ")" << std::endl;
        }

        HYPRE_Int m = static_cast<HYPRE_Int>(A.numRows());
        if (globalDofEnd == 0) globalDofEnd = globalDofStart + m;
        if (reuse_enabled_ || fixed_graph_updates_) {
            const bool same = prepared_ && matrix_ == &A
                && (fixed_graph_updates_ || values_ == A.valuesPtr())
                && row_offsets_ == A.rowOffsetsPtr() && columns_ == A.colIndicesPtr()
                && rows_ == A.numRows() && cols_ == A.numCols() && nnz_ == A.nnz()
                && map_ == thrust::raw_pointer_cast(d_localToGlobalDof.data())
                && map_size_ == d_localToGlobalDof.size()
                && globalDofStart_ == globalDofStart && globalDofEnd_ == globalDofEnd
                && global_col_start_ == globalColStart && global_col_end_ == globalColEnd;
            int rebuild = same ? 0 : 1, any_rebuild = 0;
            MPI_Allreduce(&rebuild, &any_rebuild, 1, MPI_INT, MPI_MAX, comm_);
            if (!any_rebuild) {
                if (fixed_graph_updates_) {
                    update_matrix_values(A);
                    if (!update_vectors(b, x)) return false;
                    if (profile_) profile_start_ = profile_->lap(SolverProfile::Prepare, profile_start_);
                    finish_prepare();
                    const double start = wall_stamp();
                    const HYPRE_Int error = useFlexGmres_
                        ? HYPRE_ParCSRFlexGMRESSetup(solver_, parcsr_A_, par_b_, par_x_)
                        : HYPRE_ParCSRGMRESSetup(solver_, parcsr_A_, par_b_, par_x_);
                    if (timing_enabled_) last_timing_.setup_seconds = wall_stamp() - start;
                    if (profile_) profile_start_ = profile_->lap(SolverProfile::Setup, profile_start_);
                    require_reuse(error == 0 && HYPRE_GetError() == 0, "GMRES/AMG refresh failed");
                    ++setup_count_;
                    return solve_vectors(b, x);
                }
                if (!update_vectors(b, x)) return false;
                if (!vectors_changed_) {
                    if (profile_) profile_->lap(SolverProfile::Prepare, profile_start_);
                    finish_prepare();
                    return solve_vectors(b, x);
                }
            }
        }
        destroy();
        globalDofStart_ = globalDofStart;
        globalDofEnd_   = globalDofEnd;
        matrix_ = &A; values_ = A.valuesPtr(); row_offsets_ = A.rowOffsetsPtr();
        columns_ = A.colIndicesPtr(); rows_ = A.numRows(); cols_ = A.numCols(); nnz_ = A.nnz();
        map_ = thrust::raw_pointer_cast(d_localToGlobalDof.data());
        map_size_ = d_localToGlobalDof.size();
        global_col_start_ = globalColStart; global_col_end_ = globalColEnd;

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
        last_timing_ = {};
        prepare_start_ = wall_stamp();
        if (fixed_graph_updates_)
            require_reuse(false, "fixed graph updates require the device-map overload");
        if (profile_) profile_start_ = profile_->stamp();
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

    // No-op preconditioner-setup shim. When AMG was already set up on a SEPARATE
    // matrix K (precondMatrix_ path), we register this as the precond "setup" so
    // HYPRE_GMRESSetup does NOT rebuild the AMG hierarchy on the solved matrix A.
    // The real solve still calls HYPRE_BoomerAMGSolve with the K-built hierarchy.
    static HYPRE_Int noopPrecondSetup(HYPRE_Solver, HYPRE_Matrix,
                                      HYPRE_Vector, HYPRE_Vector) { return 0; }

    // Set a SEPARATE matrix K to build the BoomerAMG preconditioner on, instead
    // of the solved matrix A. K must share A's partition (same row/col global
    // range and same local->global DOF map). Pass nullptr to revert to the
    // classic single-matrix path. See precondMatrix_ for the rationale.
    void setPrecondMatrix(const Matrix* K) {
        if (precondMatrix_ != K) invalidate_setup();
        precondMatrix_ = K;
    }

    // Env-override helpers for the BoomerAMG knobs. Return the default when the
    // var is unset or unparseable, so a typo never silently disables AMG tuning.
    static int getEnvInt(const char* name, int defVal) {
        const char* v = std::getenv(name);
        if (!v || !*v) return defVal;
        char* end = nullptr;
        long parsed = std::strtol(v, &end, 10);
        return (end && end != v) ? static_cast<int>(parsed) : defVal;
    }
    static double getEnvDouble(const char* name, double defVal) {
        const char* v = std::getenv(name);
        if (!v || !*v) return defVal;
        char* end = nullptr;
        double parsed = std::strtod(v, &end);
        return (end && end != v) ? parsed : defVal;
    }

    // Iteration count and final relative residual from the most recent solve.
    // Hypre fills these via HYPRE_GMRESGetNumIterations / GetFinalRelativeResidualNorm
    // at end of solve. Driver code reads these for per-step diagnostics so a
    // silently-non-converged AMG run (returns success but residual stagnates
    // above tol) is visible in the step log instead of buried in Hypre's own
    // print output.
    int    getLastIterations()    const { return lastNumIters_; }

    // Opt-in systems / point-block AMG: with N>1 and node-major interleaved DOFs
    // (dof = N*node + comp), BoomerAMG auto-generates dof_func[i] = i % N, so a
    // single SetNumFunctions(N) call enables point-block coarsening. Default N=1
    // (scalar) -> no behavior change for existing callers.
    void   setPointBlock(int n) {
        const int next = n > 1 ? n : 1;
        if (pointBlock_ != next) invalidate_setup();
        pointBlock_ = next;
    }
    // Keep Hypre's default unless a caller supplies a coarse-level relaxation.
    void setAMGCoarseRelaxType(int type) {
        if (coarseRelaxType_ != type) invalidate_setup();
        coarseRelaxType_ = type;
    }
    double getLastFinalResidual() const { return lastFinalRes_; }
    double getLastSolutionMax()   const { return lastSolutionMax_; }
    bool   lastReturnedNullSolution() const { return nullSolutionReturned_; }

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

        if (reuse_enabled_ || fixed_graph_updates_ || true_residual_check_) {
            const bool bad = A.nnz() > 0 && thrust::any_of(thrust::device_pointer_cast(A.valuesPtr()),
                thrust::device_pointer_cast(A.valuesPtr() + A.nnz()), IsNonFinite<RealType>());
            require_reuse(!bad && globalDofEnd_ - globalDofStart_ == m
                          && (!fixed_graph_updates_ || (m > 0 && !precondMatrix_)),
                          "nonfinite matrix, inconsistent partition, or empty partition/separate K in fixed graph mode");
            HYPRE_ClearAllErrors();
        }
        setupHypreMatrix(A, globalColStart, globalColEnd);
        if (reuse_enabled_ || fixed_graph_updates_ || true_residual_check_)
            require_reuse(parcsr_A_ != nullptr && HYPRE_GetError() == 0, "matrix preparation failed");

        // RHS + initial guess: pure device path.
        HYPRE_BigInt ilower = static_cast<HYPRE_BigInt>(globalDofStart_);
        HYPRE_BigInt iupper = static_cast<HYPRE_BigInt>(globalDofEnd_ - 1);

        HYPRE_IJVectorCreate(comm_, ilower, iupper, &b_hypre_);
        HYPRE_IJVectorCreate(comm_, ilower, iupper, &x_hypre_);
        HYPRE_IJVectorSetObjectType(b_hypre_, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(x_hypre_, HYPRE_PARCSR);
        if (reuse_enabled_ || fixed_graph_updates_ || true_residual_check_)
            require_reuse(b_hypre_ && x_hypre_ && HYPRE_GetError() == 0, "vector creation failed");

        // Device-side global row indices [ilower, ilower+m).
        auto& d_rowGlobal = d_row_global_;
        d_rowGlobal.resize(m);
        if (m > 0) {
            const int bs = 256;
            const int gs = (m + bs - 1) / bs;
            fillGlobalRowIndicesKernel<HYPRE_BigInt><<<gs, bs>>>(
                thrust::raw_pointer_cast(d_rowGlobal.data()), ilower, m);
        }

        if (!update_vectors(b, x)) return false;

        if (profile_) profile_start_ = profile_->lap(SolverProfile::Prepare, profile_start_);
        finish_prepare();
        const double setup_start = wall_stamp();

        if (verbose_ && rank == 0) std::cout << "Creating preconditioner..." << std::endl;

        if (precondType_ == BOOMERAMG) {
            if (verbose_ && rank == 0) std::cout << "Using BoomerAMG preconditioner (GPU)" << std::endl;
            HYPRE_BoomerAMGCreate(&precond_);
            if (reuse_enabled_ || fixed_graph_updates_ || true_residual_check_)
                require_reuse(precond_ && HYPRE_GetError() == 0, "AMG creation failed");
            // Env var MARS_HYPRE_VERBOSE=1 turns on per-iter Hypre prints so a
            // stalled AMG solve shows its residual history without rebuilding.
            const char* ev = std::getenv("MARS_HYPRE_VERBOSE");
            int hyprePrintLevel = (ev && std::string(ev) != "0") ? 3 : 0;
            HYPRE_BoomerAMGSetPrintLevel(precond_, hyprePrintLevel);
            // GPU-required + tuned-for-3D-pressure-Poisson parameters. The
            // assembled DDT operator (D M^-1 D^T + tau*PSPG, single-pin anchored)
            // is a near-singular SPD Laplacian-like matrix. Prior defaults
            // (RelaxType=18 l1-Jacobi, StrongThr=0.5, AggNumLevels=1) are WEAK
            // for this operator: l1-Jacobi smoothing barely damps the
            // low-frequency error, the 0.5 threshold drops too many couplings in
            // 3D so coarsening loses the operator's connectivity, and aggressive
            // coarsening on an irregular tet stencil produces a coarse operator
            // whose near-constant null mode dominates the V-cycle output. Under
            // FlexGMRES that makes the FIRST preconditioned direction collapse
            // onto the (near-)null space: the relative residual ||r||/||b|| drops
            // below tol in ~2 iters while the actual solution x stays ~0 (the
            // false x=0 "convergence"). Hybrid symmetric Gauss-Seidel (RelaxType
            // 6/8) is the standard strong smoother for pressure Poisson and
            // damps the low modes the V-cycle must remove; StrongThr=0.25 is the
            // 3D-Poisson value; AggNumLevels=0 keeps the full operator
            // connectivity so the coarse grid still represents the constant mode
            // correctly. All env-overridable so they can be swept without a
            // rebuild.
            //   CoarsenType=8  = PMIS (mandatory on GPU; cheaper than CLJP)
            //   InterpType=6   = ext+i (long-range, works with PMIS)
            //   RelaxType=18   = l1-Jacobi. MANDATORY on GPU: the GS-family
            //                   smoothers (hybrid GS=3/6, l1-hybrid-SSOR=8)
            //                   become non-symmetric / NO-OP on GPU (the cuda
            //                   Hypre build does not implement them as a real
            //                   sweep), so AMG returns x=0 with a garbage ~1e-315
            //                   residual on the first cycle. The working PCG
            //                   wrapper uses 18 for exactly this reason (not just
            //                   PCG symmetry). Earlier 8 here was the cause of the
            //                   pump's AMG x=0. Override via MARS_AMG_RELAX.
            //   RelaxOrder=0   = lexicographic (mandatory on GPU)
            //   KeepTranspose=1= avoid SpMTV on GPU
            //   StrongThr=0.25 = 3D-Poisson value (0.5 is for anisotropic)
            //   PMaxElmts=4    = truncate P (GPU memory + perf)
            //   AggNumLevels=0 = no aggressive coarsening (keeps tet connectivity)
            int   amgCoarsen   = getEnvInt   ("MARS_AMG_COARSEN",   8);
            int   amgInterp    = getEnvInt   ("MARS_AMG_INTERP",    6);
            int   amgRelax     = getEnvInt   ("MARS_AMG_RELAX",    18);
            int   amgRelaxOrder= getEnvInt   ("MARS_AMG_RELAXORDER",0);
            double amgStrong   = getEnvDouble("MARS_AMG_STRONG",    0.25);
            int   amgPMax      = getEnvInt   ("MARS_AMG_PMAX",      4);
            int   amgAgg       = getEnvInt   ("MARS_AMG_AGG",       0);
            int   amgSweeps    = getEnvInt   ("MARS_AMG_SWEEPS",    2);  // l1-Jacobi(18) is weaker than SSOR; 2 sweeps restore smoothing
            if (verbose_ && rank == 0) {
                std::cout << "  [BoomerAMG] coarsen=" << amgCoarsen
                          << " interp=" << amgInterp << " relax=" << amgRelax
                          << " strong=" << amgStrong << " pmax=" << amgPMax
                          << " agg=" << amgAgg << " sweeps=" << amgSweeps
                          << std::endl;
            }
            HYPRE_BoomerAMGSetCoarsenType(precond_, amgCoarsen);
            HYPRE_BoomerAMGSetInterpType(precond_, amgInterp);
            HYPRE_BoomerAMGSetRelaxType(precond_, amgRelax);
            if (coarseRelaxType_ >= 0)
                HYPRE_BoomerAMGSetCycleRelaxType(precond_, coarseRelaxType_, 3);
            HYPRE_BoomerAMGSetRelaxOrder(precond_, amgRelaxOrder);
            HYPRE_BoomerAMGSetKeepTranspose(precond_, 1);
            HYPRE_BoomerAMGSetStrongThreshold(precond_, amgStrong);
            HYPRE_BoomerAMGSetPMaxElmts(precond_, amgPMax);
            HYPRE_BoomerAMGSetAggNumLevels(precond_, amgAgg);
            HYPRE_BoomerAMGSetNumSweeps(precond_, amgSweeps);
            // Systems / point-block AMG: coarsen the N-DOF nodal block together
            // (Hypre auto-builds dof_func = i % N for node-major interleaved DOFs).
            if (pointBlock_ > 1) {
                HYPRE_BoomerAMGSetNumFunctions(precond_, pointBlock_);
                if (verbose_ && rank == 0)
                    std::cout << "  [BoomerAMG] point-block NumFunctions=" << pointBlock_ << std::endl;
            }
            HYPRE_BoomerAMGSetMaxLevels(precond_, 25);
            HYPRE_BoomerAMGSetMinCoarseSize(precond_, 32);  // avoid too-small coarse grid on small problems
            HYPRE_BoomerAMGSetMaxCoarseSize(precond_, 128); // upper bound on coarsest direct solve
            HYPRE_BoomerAMGSetTol(precond_, 0.0);           // GMRES controls outer tol
            HYPRE_BoomerAMGSetMaxIter(precond_, 1);         // 1 V-cycle per GMRES iter
            // Note: an earlier attempt called HYPRE_BoomerAMGSetInterpVectors
            // with the constant-of-ones to declare the near-null mode of the
            // pressure-Poisson. That call is REJECTED with HYPRE_ERROR_GENERIC
            // by hypre/parcsr_ls/par_amg_setup.c:506 unless paired with
            // SetNodal(1) + SetInterpVecVariant(2) + SetInterpVecQMax(...) +
            // SetSmoothInterpVectors(0) (cf. MFEM hypre.cpp:5152-5160).
            // Since the assembled DDT operator also fails BoomerAMG for a
            // separate reason -- its boundary rows hold POSITIVE off-diagonals
            // (non-M-matrix), which BoomerAMG's classical strength-of-
            // connection drops -- we instead route DDT pressure through the
            // matrix-free CG path (with Jacobi preconditioning) entirely, and
            // keep this wrapper for the K-path / velocity solves where the
            // matrix IS an M-matrix and AMG works as designed.
            // MARS_HYPRE_VERBOSE=1 also enables BoomerAMG setup-phase prints
            // (level structure, complexity, row sums per level).
            {
                const char* ev = std::getenv("MARS_HYPRE_VERBOSE");
                if (ev && std::string(ev) != "0") {
                    HYPRE_BoomerAMGSetPrintLevel(precond_, 3);
                }
            }
            // Schur-complement preconditioning: build the V-cycle on the separate
            // K (AMG-friendly Galerkin stiffness) instead of the solved A (the
            // AMG-hostile Gram operator A_fem). We run BoomerAMGSetup HERE on
            // parcsr_K_, then register AMG with a no-op setup shim below so
            // HYPRE_GMRESSetup does NOT rebuild the hierarchy on parcsr_A_. K
            // shares A_fem's partition exactly, so the same d_localToGlobalDof_
            // and global range apply.
            if (precondMatrix_ != nullptr) {
                buildParCsr(*precondMatrix_, globalColStart, globalColEnd,
                            K_hypre_, parcsr_K_);
                if (!parcsr_K_) {
                    std::cerr << "[HypreGMRES] precond matrix K build failed; "
                                 "falling back to AMG-on-A\n";
                } else {
                    if (verbose_ && rank == 0)
                        std::cout << "  [BoomerAMG] building hierarchy on separate "
                                     "K preconditioner matrix\n";
                    HYPRE_Int kSetupErr = HYPRE_BoomerAMGSetup(
                        precond_, parcsr_K_, par_b_, par_x_);
                    if (kSetupErr != 0 && rank == 0) {
                        std::cerr << "[HypreGMRES] BoomerAMGSetup on K returned "
                                  << kSetupErr << "\n";
                        HYPRE_ClearAllErrors();
                    }
                }
            }
        } else if (precondType_ == JACOBI) {
            if (verbose_ && rank == 0) std::cout << "Using Jacobi preconditioner" << std::endl;
            precond_ = (HYPRE_Solver) parcsr_A_;
        }

        if (amg_cycle_) {
            require_reuse(precondType_ == BOOMERAMG && precond_ && !precondMatrix_,
                          "direct cycle requires BoomerAMG on Apre");
            const HYPRE_Int error = HYPRE_BoomerAMGSetup(precond_, parcsr_A_, par_b_, par_x_);
            if (timing_enabled_) last_timing_.setup_seconds = wall_stamp() - setup_start;
            if (profile_) profile_->lap(SolverProfile::Setup, profile_start_);
            require_reuse(error == 0 && HYPRE_GetError() == 0, "AMG setup failed");
            ++setup_count_;
            prepared_ = true;
            return solve_vectors(b, x);
        }

        const char* krylovName = useFlexGmres_ ? "FlexGMRES" : "GMRES";
        if (verbose_ && rank == 0) std::cout << "Creating " << krylovName << " solver..." << std::endl;
        if (useFlexGmres_) HYPRE_ParCSRFlexGMRESCreate(comm_, &solver_);
        else               HYPRE_ParCSRGMRESCreate(comm_, &solver_);
        if (reuse_enabled_ || fixed_graph_updates_ || true_residual_check_)
            require_reuse(solver_ && HYPRE_GetError() == 0, "GMRES creation failed");
        if (!solver_) {
            std::cerr << "Failed to create Hypre " << krylovName << " solver" << std::endl;
            return false;
        }
        // These handles have different layouts; every API must match the solver.
        (useFlexGmres_ ? HYPRE_FlexGMRESSetMaxIter : HYPRE_GMRESSetMaxIter)(solver_, maxIter_);
        (useFlexGmres_ ? HYPRE_FlexGMRESSetTol : HYPRE_GMRESSetTol)(solver_, tolerance_);
        (useFlexGmres_ ? HYPRE_FlexGMRESSetKDim : HYPRE_GMRESSetKDim)(solver_, kDim_);
        // Floor on iterations. For the near-singular DDT pressure operator (one
        // constant null mode, single pin), BoomerAMG's first V-cycle can map the
        // initial residual almost entirely into the null space, dropping the
        // RELATIVE residual ||r||/||b|| below tol in ~2 iters while the actual
        // solution x is still ~0 (the false x=0 "convergence"). A small minimum
        // iteration count forces GMRES to keep building the Krylov space past
        // that first deceptive cycle so a real x emerges. Env-overridable.
        {
            int minIt = getEnvInt("MARS_HYPRE_MINITER", 3);
            if (minIt > 0)
                (useFlexGmres_ ? HYPRE_FlexGMRESSetMinIter : HYPRE_GMRESSetMinIter)(solver_, minIt);
        }
        // Optional absolute stopping floor: Hypre uses max(atol, rtol*||b||).
        // Default 0 keeps the relative target; the wrapper checks acceptance below.
        {
            double absTol = getEnvDouble("MARS_HYPRE_ABSTOL", 0.0);
            if (absTol > 0.0)
                (useFlexGmres_ ? HYPRE_FlexGMRESSetAbsoluteTol : HYPRE_GMRESSetAbsoluteTol)(solver_, absTol);
        }
        // print level: 0 silent, 2 per-iter residuals. Env MARS_HYPRE_VERBOSE=1.
        {
            const char* ev = std::getenv("MARS_HYPRE_VERBOSE");
            int gmresPrint = (verbose_ || (ev && std::string(ev) != "0")) ? 2 : 0;
            (useFlexGmres_ ? HYPRE_FlexGMRESSetPrintLevel : HYPRE_GMRESSetPrintLevel)(solver_, gmresPrint);
        }

        if (precondType_ == BOOMERAMG && precond_) {
            if (verbose_ && rank == 0) std::cout << "Setting BoomerAMG preconditioner..." << std::endl;
            // When AMG was pre-built on the separate K (parcsr_K_), register a
            // no-op setup so GMRESSetup does not rebuild the hierarchy on A.
            const bool kPrecondReady = (precondMatrix_ != nullptr && parcsr_K_ != nullptr);
            auto amgSetupFn = kPrecondReady
                ? (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))noopPrecondSetup
                : (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSetup;
            if (useFlexGmres_)
                HYPRE_FlexGMRESSetPrecond(solver_,
                                   (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSolve,
                                   amgSetupFn,
                                   precond_);
            else
                HYPRE_GMRESSetPrecond(solver_,
                                   (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSolve,
                                   amgSetupFn,
                                   precond_);
        } else if (precondType_ == JACOBI) {
            if (verbose_ && rank == 0) std::cout << "Setting Jacobi (diagonal) preconditioner..." << std::endl;
            if (useFlexGmres_)
                HYPRE_FlexGMRESSetPrecond(solver_,
                                   (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_ParCSRDiagScale,
                                   (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_ParCSRDiagScaleSetup,
                                   (HYPRE_Solver) parcsr_A_);
            else
                HYPRE_GMRESSetPrecond(solver_,
                                   (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_ParCSRDiagScale,
                                   (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_ParCSRDiagScaleSetup,
                                   (HYPRE_Solver) parcsr_A_);
        }

        if (verbose_ && rank == 0) std::cout << "Setting up " << krylovName << " solver..." << std::endl;
        if (!fixed_graph_updates_) MPI_Barrier(comm_);
        HYPRE_Int setup_err = useFlexGmres_
            ? HYPRE_ParCSRFlexGMRESSetup(solver_, parcsr_A_, par_b_, par_x_)
            : HYPRE_ParCSRGMRESSetup(solver_, parcsr_A_, par_b_, par_x_);
        if (timing_enabled_) last_timing_.setup_seconds = wall_stamp() - setup_start;
        if (profile_) profile_start_ = profile_->lap(SolverProfile::Setup, profile_start_);
        if (reuse_enabled_ || fixed_graph_updates_ || true_residual_check_)
            require_reuse(setup_err == 0 && HYPRE_GetError() == 0, "GMRES/AMG setup failed");
        if (setup_err != 0 && rank == 0) {
            std::cerr << "[HypreGMRES] Setup returned error " << setup_err
                      << " (HYPRE_GetError=" << HYPRE_GetError() << ")\n";
            char errbuf[256];
            HYPRE_DescribeError(HYPRE_GetError(), errbuf);
            std::cerr << "[HypreGMRES] " << errbuf << "\n";
            HYPRE_ClearAllErrors();
        }
        if (verbose_ && rank == 0) std::cout << "GMRES setup complete, starting solve..." << std::endl;

        if (setup_err == 0) ++setup_count_;
        prepared_ = (reuse_enabled_ || fixed_graph_updates_) && setup_err == 0;
        return solve_vectors(b, x);
    }

    // Reinitialize the same IJ vectors; their partition and storage stay fixed.
    bool update_vectors(const Vector& b, Vector& x) {
        int rank = 0;
        MPI_Comm_rank(comm_, &rank);
        const HYPRE_Int m = globalDofEnd_ - globalDofStart_;
        auto& d_rowGlobal = d_row_global_;
        const auto previous_b = par_b_, previous_x = par_x_;
        if (reuse_enabled_) {
            require_reuse(b.size() >= size_t(m) && x.size() >= size_t(m), "undersized RHS or solution");
            const bool bad_b = m > 0 && thrust::any_of(thrust::device_pointer_cast(b.data()),
                thrust::device_pointer_cast(b.data() + m), IsNonFinite<RealType>());
            const bool bad_x = m > 0 && thrust::any_of(thrust::device_pointer_cast(x.data()),
                thrust::device_pointer_cast(x.data() + m), IsNonFinite<RealType>());
            require_reuse(!bad_b && !bad_x, "nonfinite RHS or initial guess");
            HYPRE_ClearAllErrors();
        } else if (fixed_graph_updates_) {
            const bool sized = b.size() >= size_t(m) && x.size() >= size_t(m);
            const bool bad_b = sized && m > 0 && thrust::any_of(thrust::device_pointer_cast(b.data()),
                thrust::device_pointer_cast(b.data() + m), IsNonFinite<RealType>());
            const bool bad_x = sized && m > 0 && thrust::any_of(thrust::device_pointer_cast(x.data()),
                thrust::device_pointer_cast(x.data() + m), IsNonFinite<RealType>());
            require_reuse(sized && !bad_b && !bad_x, "undersized or nonfinite RHS/initial guess");
            HYPRE_ClearAllErrors();
        }
        HYPRE_IJVectorInitialize(b_hypre_);
        HYPRE_IJVectorInitialize(x_hypre_);
        // RHS NaN/Inf summary + min/max/sum, all on device.
        if (!fixed_graph_updates_ || verbose_) validateVector(b.data(), m, rank, "RHS");

        // Hypre RHS values must be HYPRE_Real; cast on device if RealType != HYPRE_Real.
        const HYPRE_Real* d_b_hypre = nullptr;
        auto& d_b_cast = d_b_cast_;
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
        auto& d_x_cast = d_x_cast_;
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

        if (reuse_enabled_) {
            require_reuse(par_b_ && par_x_ && HYPRE_GetError() == 0, "vector update failed");
            int changed = prepared_ && (previous_b != par_b_ || previous_x != par_x_);
            int any_changed = 0;
            MPI_Allreduce(&changed, &any_changed, 1, MPI_INT, MPI_MAX, comm_);
            vectors_changed_ = any_changed != 0;
        } else if (fixed_graph_updates_ || true_residual_check_) {
            require_reuse(par_b_ && par_x_ && HYPRE_GetError() == 0, "vector update failed");
        }

        if (!par_b_ || !par_x_) {
            std::cerr << "Rank " << rank << ": Failed to get Hypre ParVector objects" << std::endl;
            return false;
        }

        return true;
    }

    bool solve_vectors(const Vector& b, Vector& x) {
        int rank = 0;
        MPI_Comm_rank(comm_, &rank);
        const HYPRE_Int m = globalDofEnd_ - globalDofStart_;
        auto& d_rowGlobal = d_row_global_;
        auto& d_x_cast = d_x_cast_;
        if (profile_) profile_start_ = profile_->stamp();
        const double solve_start = wall_stamp();
        if (!fixed_graph_updates_) MPI_Barrier(comm_);
        HYPRE_Int solve_err = amg_cycle_
            ? HYPRE_BoomerAMGSolve(precond_, parcsr_A_, par_b_, par_x_)
            : useFlexGmres_
            ? HYPRE_ParCSRFlexGMRESSolve(solver_, parcsr_A_, par_b_, par_x_)
            : HYPRE_ParCSRGMRESSolve(solver_, parcsr_A_, par_b_, par_x_);
        if (timing_enabled_) last_timing_.solve_seconds = wall_stamp() - solve_start;
        const double finish_start = wall_stamp();
        if (profile_) profile_start_ = profile_->lap(SolverProfile::Solve, profile_start_);
        if (reuse_enabled_ || fixed_graph_updates_ || true_residual_check_) {
            require_reuse((amg_cycle_ ? solve_err : (solve_err & ~HYPRE_ERROR_CONV)) == 0,
                          "Hypre apply failed");
            HYPRE_ClearAllErrors();
        }
        if (solve_err != 0 && rank == 0) {
            std::cerr << "[HypreGMRES] Solve returned error " << solve_err << '\n';
            char errbuf[256];
            HYPRE_DescribeError(solve_err, errbuf);
            std::cerr << "[HypreGMRES] " << errbuf << "\n";
            HYPRE_ClearAllErrors();
        }
        if (verbose_ && rank == 0) std::cout << "GMRES solve complete." << std::endl;

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

        if (reuse_enabled_ || fixed_graph_updates_ || true_residual_check_) {
            const bool bad = m > 0 && thrust::any_of(thrust::device_pointer_cast(x.data()),
                thrust::device_pointer_cast(x.data() + m), IsNonFinite<RealType>());
            require_reuse(!bad && HYPRE_GetError() == 0 && cudaGetLastError() == cudaSuccess,
                          "nonfinite result, extraction error, or CUDA failure");
        }
        if (amg_cycle_) {
            const double local[2] = {
                m > 0 ? thrust::transform_reduce(thrust::device_pointer_cast(x.data()),
                    thrust::device_pointer_cast(x.data() + m),
                    [] __device__(RealType v) -> double { return fabs(double(v)); },
                    0.0, thrust::maximum<double>()) : 0.0,
                m > 0 ? thrust::transform_reduce(thrust::device_pointer_cast(b.data()),
                    thrust::device_pointer_cast(b.data() + m),
                    [] __device__(RealType v) -> double { return fabs(double(v)); },
                    0.0, thrust::maximum<double>()) : 0.0};
            double global[2] = {};
            MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MAX, comm_);
            lastNumIters_ = 1;
            lastFinalRes_ = std::numeric_limits<double>::quiet_NaN();
            lastSolutionMax_ = global[0];
            nullSolutionReturned_ = global[1] > 0 && global[0] == 0;
            if (profile_) profile_->lap(SolverProfile::Finish, profile_start_);
            if (timing_enabled_) last_timing_.finish_seconds = wall_stamp() - finish_start;
            // A cycle is a preconditioner action, not a converged linear solve.
            return !nullSolutionReturned_;
        }

        int    num_iterations = 0;
        double final_res_norm = 0.0;
        (useFlexGmres_ ? HYPRE_FlexGMRESGetNumIterations : HYPRE_GMRESGetNumIterations)(solver_, &num_iterations);
        (useFlexGmres_ ? HYPRE_FlexGMRESGetFinalRelativeResidualNorm : HYPRE_GMRESGetFinalRelativeResidualNorm)(solver_, &final_res_norm);
        const double reported_res_norm = final_res_norm;
        // Hypre's zero-residual early return can leave the previous solve's norm.
        if (true_residual_check_ || (fixed_graph_updates_ && num_iterations == 0))
            final_res_norm = true_relative_residual();
        lastNumIters_ = num_iterations;
        lastFinalRes_ = final_res_norm;

        // Keep legacy magnitude guards for projection callers. Anchored SIMPLE
        // systems use the explicit residual: pressure and mass RHS have different
        // units, and changing the matrix scale changes x/b without harming a solve.
        {
            double localXmax = (m > 0)
                ? thrust::transform_reduce(
                      thrust::device_pointer_cast(x.data()),
                      thrust::device_pointer_cast(x.data() + m),
                      [] __device__ (RealType v) -> double { return fabs((double)v); },
                      0.0, thrust::maximum<double>())
                : 0.0;
            double localBmax = (m > 0)
                ? thrust::transform_reduce(
                      thrust::device_pointer_cast(b.data()),
                      thrust::device_pointer_cast(b.data() + m),
                      [] __device__ (RealType v) -> double { return fabs((double)v); },
                      0.0, thrust::maximum<double>())
                : 0.0;
            double gXmax = 0.0, gBmax = 0.0;
            if (fixed_graph_updates_) {
                const double local[2] = {localXmax, localBmax};
                double global[2] = {};
                MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MAX, comm_);
                gXmax = global[0]; gBmax = global[1];
            } else {
                MPI_Allreduce(&localXmax, &gXmax, 1, MPI_DOUBLE, MPI_MAX, comm_);
                MPI_Allreduce(&localBmax, &gBmax, 1, MPI_DOUBLE, MPI_MAX, comm_);
            }
            lastSolutionMax_ = gXmax;
            double nullRatio = getEnvDouble("MARS_HYPRE_NULLX_RATIO", 1e-12);
            nullSolutionReturned_ = !true_residual_check_ && (gBmax > 0.0 && gXmax < nullRatio * gBmax);
            if (nullSolutionReturned_ && rank == 0) {
                std::cout << "[HypreGMRES] WARNING: near-zero solution (iters="
                          << num_iterations << ", rel_res=" << final_res_norm
                          << ") but |x|inf=" << gXmax << " << |b|inf=" << gBmax
                          << " -- reporting NOT converged; cause is undetermined.\n";
            }
            double maxXRatio = getEnvDouble("MARS_HYPRE_MAXX_RATIO", 1e6);
            if (!true_residual_check_ && gBmax > 0.0 && gXmax > maxXRatio * gBmax) {
                nullSolutionReturned_ = true;
                if (rank == 0)
                    std::cout << "[HypreGMRES] WARNING: |x|inf=" << gXmax
                              << " >> |b|inf=" << gBmax << " (ratio>" << maxXRatio
                              << ") -- under-resolved/null-contaminated solution; "
                              << "reporting NOT converged.\n";
            }
        }

        const double residual_limit = residual_absolute_tolerance_
                                    + residual_relative_tolerance_ * last_rhs_norm_;
        const bool residual_ok = mixed_residual_check_
            ? std::isfinite(residual_limit) && last_absolute_residual_ <= residual_limit
            : final_res_norm < tolerance_;
        const bool converged = residual_ok && !nullSolutionReturned_;
        if (true_residual_check_ && !converged) {
            // Read the solver's work vector without replacing our explicit check.
            // After an early stop it need not contain the final b-Ax.
            HYPRE_Real work_residual2 = std::numeric_limits<HYPRE_Real>::quiet_NaN();
            if (num_iterations > 0) {
                HYPRE_ParVector work_residual = nullptr;
                (useFlexGmres_ ? HYPRE_ParCSRFlexGMRESGetResidual : HYPRE_ParCSRGMRESGetResidual)(solver_, &work_residual);
                require_reuse(work_residual && HYPRE_GetError() == 0, "Krylov residual lookup failed");
                HYPRE_ParVectorInnerProd(work_residual, work_residual, &work_residual2);
                require_reuse(HYPRE_GetError() == 0, "Krylov residual norm failed");
            }
            if (rank == 0) {
                std::cerr << "[HypreGMRES] rejected: backend=" << (useFlexGmres_ ? "FlexGMRES" : "GMRES")
                          << " iterations=" << num_iterations
                          << '/' << maxIter_ << " restart=" << kDim_
                          << " reported_relative=" << reported_res_norm
                          << " true_relative_or_absolute=" << final_res_norm
                          << " tolerance=" << tolerance_ << " solve_error=" << solve_err
                          << " krylov_work_norm=" << std::sqrt(work_residual2)
                          << " rhs_norm=" << last_rhs_norm_;
                if (mixed_residual_check_)
                    std::cerr << " absolute_residual=" << last_absolute_residual_
                              << " acceptance_limit=" << residual_limit;
                std::cerr << '\n';
            }
        }
        if (verbose_) {
            std::cout << "Hypre GMRES " << (converged ? "converged" : "did not converge")
                      << " in " << num_iterations
                      << " iterations, final residual: " << final_res_norm << std::endl;
        }

        // Legacy callers also retain their magnitude safeguards.
        if (profile_) profile_->lap(SolverProfile::Finish, profile_start_);
        if (timing_enabled_) last_timing_.finish_seconds = wall_stamp() - finish_start;
        return converged;
    }

    double true_relative_residual() {
        if (!r_hypre_) {
            HYPRE_IJVectorCreate(comm_, static_cast<HYPRE_BigInt>(globalDofStart_),
                                 static_cast<HYPRE_BigInt>(globalDofEnd_ - 1), &r_hypre_);
            HYPRE_IJVectorSetObjectType(r_hypre_, HYPRE_PARCSR);
            HYPRE_IJVectorInitialize(r_hypre_);
            HYPRE_IJVectorAssemble(r_hypre_);
            HYPRE_IJVectorGetObject(r_hypre_, reinterpret_cast<void**>(&par_r_));
            require_reuse(par_r_ && HYPRE_GetError() == 0, "residual vector creation failed");
        }
        // Keep b-Ax on Hypre's compute stream without a runtime memcpy before SpMV.
        HYPRE_Int error = HYPRE_ParCSRMatrixMatvec(-1.0, parcsr_A_, par_x_, 0.0, par_r_);
        error |= HYPRE_ParVectorAxpy(1.0, par_b_, par_r_);
        HYPRE_Real residual2 = 0, rhs2 = 0;
        error |= HYPRE_ParVectorInnerProd(par_r_, par_r_, &residual2);
        error |= HYPRE_ParVectorInnerProd(par_b_, par_b_, &rhs2);
        require_reuse(error == 0 && HYPRE_GetError() == 0 && std::isfinite(residual2) && std::isfinite(rhs2)
                      && residual2 >= 0 && rhs2 >= 0,
                      "true residual evaluation failed");
        last_absolute_residual_ = std::sqrt(residual2);
        last_rhs_norm_ = std::sqrt(rhs2);
        return std::sqrt(rhs2 > 0 ? residual2 / rhs2 : residual2);
    }

    // Build A_hypre_ entirely on the device:
    //   1) count valid entries per row + diagonal presence
    //   2) exclusive-scan to per-row write offsets (total = filtered nnz)
    //   3) compact (globalCol, value) into dense per-row slots
    //   4) one HYPRE_IJMatrixSetValues call with all device pointers
    void setupHypreMatrix(const Matrix& A,
                          IndexType globalColStart, IndexType globalColEnd) {
        if (fixed_graph_updates_) {
            build_fixed_graph(A);
            update_matrix_values(A);
            return;
        }
        // Default target: the solved operator's handles.
        buildParCsr(A, globalColStart, globalColEnd, A_hypre_, parcsr_A_);
    }

    void build_fixed_graph(const Matrix& A) {
        const double start = wall_stamp();
        const int m = static_cast<int>(A.numRows());
        const HYPRE_BigInt lower = globalDofStart_, upper = globalDofEnd_ - 1;
        d_packed_counts_.resize(m);
        d_graph_diagonal_.resize(m);
        d_packed_offsets_.resize(m);
        d_row_global_.resize(m);
        if (m > 0) {
            countValidPerRowKernel<IndexType, HYPRE_BigInt><<<(m + 255) / 256, 256>>>(
                A.rowOffsetsPtr(), A.colIndicesPtr(),
                thrust::raw_pointer_cast(d_localToGlobalDof_.data()), d_localToGlobalDof_.size(),
                lower, m, thrust::raw_pointer_cast(d_packed_counts_.data()),
                thrust::raw_pointer_cast(d_graph_diagonal_.data()));
            fillGlobalRowIndicesKernel<HYPRE_BigInt><<<(m + 255) / 256, 256>>>(
                thrust::raw_pointer_cast(d_row_global_.data()), lower, m);
        }
        thrust::exclusive_scan(d_packed_counts_.begin(), d_packed_counts_.end(), d_packed_offsets_.begin());
        require_reuse(thrust::count(d_packed_counts_.begin(), d_packed_counts_.end(), 0) == 0
                      && thrust::count(d_graph_diagonal_.begin(), d_graph_diagonal_.end(), 0) == 0,
                      "fixed graph contains an empty row or a missing diagonal");
        // Hypre's insertion API needs a host count for allocating the packed arrays.
        const int count = thrust::reduce(d_packed_counts_.begin(), d_packed_counts_.end(), 0);
        d_packed_columns_.resize(count);
        d_source_slots_.resize(count);
        d_packed_values_.resize(count);
        if (m > 0) {
            pack_hypre_graph_kernel<IndexType><<<(m + 255) / 256, 256>>>(
                A.rowOffsetsPtr(), A.colIndicesPtr(),
                thrust::raw_pointer_cast(d_localToGlobalDof_.data()), d_localToGlobalDof_.size(),
                thrust::raw_pointer_cast(d_packed_offsets_.data()), m,
                thrust::raw_pointer_cast(d_packed_columns_.data()),
                thrust::raw_pointer_cast(d_source_slots_.data()));
        }
        require_reuse(cudaGetLastError() == cudaSuccess, "graph packing CUDA launch failed");
        HYPRE_IJMatrixCreate(comm_, lower, upper, lower, upper, &A_hypre_);
        HYPRE_IJMatrixSetObjectType(A_hypre_, HYPRE_PARCSR);
        require_reuse(A_hypre_ && HYPRE_GetError() == 0, "matrix creation failed");
        ++graph_build_count_;
        if (timing_enabled_) last_timing_.packing_seconds += wall_stamp() - start;
    }

    void update_matrix_values(const Matrix& A) {
        const double start = wall_stamp();
        const bool bad = A.nnz() > 0 && thrust::any_of(thrust::device_pointer_cast(A.valuesPtr()),
            thrust::device_pointer_cast(A.valuesPtr() + A.nnz()), IsNonFinite<RealType>());
        HYPRE_ClearAllErrors();
        thrust::gather(d_source_slots_.begin(), d_source_slots_.end(),
                       thrust::device_pointer_cast(A.valuesPtr()), d_packed_values_.begin());
        require_reuse(!bad && cudaGetLastError() == cudaSuccess, "nonfinite matrix or numeric packing failure");
        if (timing_enabled_) last_timing_.packing_seconds += wall_stamp() - start;
        const auto previous = parcsr_A_;
        HYPRE_IJMatrixInitialize(A_hypre_);
        // Write every graph entry, including entries that changed to zero.
        HYPRE_IJMatrixSetValues(A_hypre_, static_cast<HYPRE_Int>(A.numRows()),
            thrust::raw_pointer_cast(d_packed_counts_.data()),
            thrust::raw_pointer_cast(d_row_global_.data()),
            thrust::raw_pointer_cast(d_packed_columns_.data()),
            thrust::raw_pointer_cast(d_packed_values_.data()));
        HYPRE_IJMatrixAssemble(A_hypre_);
        HYPRE_IJMatrixGetObject(A_hypre_, reinterpret_cast<void**>(&parcsr_A_));
        require_reuse(parcsr_A_ && (!previous || previous == parcsr_A_) && HYPRE_GetError() == 0,
                      "matrix refresh failed or changed the ParCSR object");
        ++numeric_update_count_;
    }

    // Build a ParCSR from a device CSR into the GIVEN handles, using the same
    // partition (globalDofStart_/End_, d_localToGlobalDof_) as the solved matrix.
    // Used both for the solved operator A and, on the K-preconditioner path, for
    // the separate precond matrix K (which shares A_fem's partition exactly).
    void buildParCsr(const Matrix& A,
                     IndexType globalColStart, IndexType globalColEnd,
                     HYPRE_IJMatrix& ij_out, HYPRE_ParCSRMatrix& parcsr_out) {
        const double packing_start = wall_stamp();
        int rank;
        MPI_Comm_rank(comm_, &rank);
        HYPRE_Int m       = static_cast<HYPRE_Int>(A.numRows());
        HYPRE_BigInt ilower = static_cast<HYPRE_BigInt>(globalDofStart_);
        HYPRE_BigInt iupper = static_cast<HYPRE_BigInt>(globalDofEnd_ - 1);

        // Column partitioning must match row partitioning so x is distributed like b.
        // Hypre handles off-processor (ghost) columns in ParCSR automatically.
        (void)globalColStart;
        (void)globalColEnd;

        HYPRE_IJMatrixCreate(comm_, ilower, iupper, ilower, iupper, &ij_out);
        if (reuse_enabled_ || fixed_graph_updates_ || true_residual_check_)
            require_reuse(ij_out && HYPRE_GetError() == 0, "matrix creation failed");
        HYPRE_IJMatrixSetObjectType(ij_out, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(ij_out);

        // Quick NaN/Inf scan on raw values (single allreduce-style reduction).
        bool hasNaN = A.nnz() > 0 && thrust::any_of(thrust::device_pointer_cast(A.valuesPtr()),
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
        if (m > 0) {
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
        if (m > 0) {
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
        if (m > 0) {
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
        if (timing_enabled_) last_timing_.packing_seconds += wall_stamp() - packing_start;
        HYPRE_IJMatrixSetValues(ij_out, m,
                                thrust::raw_pointer_cast(d_perRowCount.data()),
                                thrust::raw_pointer_cast(d_rows.data()),
                                thrust::raw_pointer_cast(d_colsGlobalCompact.data()),
                                thrust::raw_pointer_cast(d_valsCompact.data()));

        HYPRE_IJMatrixAssemble(ij_out);
        if (verbose_) std::cout << "Rank " << rank << ": Matrix assembled, getting ParCSR object..." << std::endl;
        HYPRE_IJMatrixGetObject(ij_out, (void**)&parcsr_out);

        if (!parcsr_out) {
            std::cerr << "Rank " << rank << ": Failed to get Hypre ParCSR matrix object" << std::endl;
            return;
        }
        ++graph_build_count_;
        ++numeric_update_count_;
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
        prepared_ = false;
        if (solver_) {
            if (useFlexGmres_) HYPRE_ParCSRFlexGMRESDestroy(solver_);
            else               HYPRE_ParCSRGMRESDestroy(solver_);
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
        if (K_hypre_) {
            HYPRE_IJMatrixDestroy(K_hypre_);
            K_hypre_ = nullptr;
            parcsr_K_ = nullptr;
        }
        if (b_hypre_) {
            HYPRE_IJVectorDestroy(b_hypre_);
            b_hypre_ = nullptr;
        }
        if (x_hypre_) {
            HYPRE_IJVectorDestroy(x_hypre_);
            x_hypre_ = nullptr;
        }
        if (r_hypre_) {
            HYPRE_IJVectorDestroy(r_hypre_);
            r_hypre_ = nullptr;
        }
        parcsr_A_ = nullptr;
        par_b_ = nullptr;
        par_x_ = nullptr;
        par_r_ = nullptr;
    }

    void require_reuse(bool local_ok, const char* message) const {
        int bad = local_ok ? 0 : 1, any_bad = 0;
        MPI_Allreduce(&bad, &any_bad, 1, MPI_INT, MPI_MAX, comm_);
        if (any_bad) {
            int rank = 0;
            MPI_Comm_rank(comm_, &rank);
            if (rank == 0) std::cerr << "ERROR: prepared Hypre: " << message << '\n';
            MPI_Abort(comm_, 1);
            std::abort();
        }
    }

private:
    double wall_stamp() const { return timing_enabled_ ? MPI_Wtime() : 0; }
    void finish_prepare() {
        if (timing_enabled_) last_timing_.prepare_seconds = wall_stamp() - prepare_start_;
    }
    bool fixed_graph_updates_ = false, timing_enabled_ = false, true_residual_check_ = false;
    bool mixed_residual_check_ = false;
    double residual_absolute_tolerance_ = 0, residual_relative_tolerance_ = 0;
    int graph_build_count_ = 0, numeric_update_count_ = 0;
    SolveTiming last_timing_;
    double prepare_start_ = 0;
    thrust::device_vector<int> d_packed_counts_, d_packed_offsets_, d_graph_diagonal_;
    thrust::device_vector<IndexType> d_source_slots_;
    thrust::device_vector<HYPRE_BigInt> d_packed_columns_;
    thrust::device_vector<HYPRE_Real> d_packed_values_;
    bool reuse_enabled_ = false, amg_cycle_ = false, prepared_ = false, vectors_changed_ = false;
    int setup_count_ = 0;
    const Matrix* matrix_ = nullptr;
    const void *values_ = nullptr, *row_offsets_ = nullptr, *columns_ = nullptr, *map_ = nullptr;
    size_t map_size_ = 0;
    IndexType rows_ = 0, cols_ = 0, nnz_ = 0, global_col_start_ = 0, global_col_end_ = 0;
    thrust::device_vector<HYPRE_BigInt> d_row_global_;
    thrust::device_vector<HYPRE_Real> d_b_cast_, d_x_cast_;
    SolverProfile* profile_ = nullptr;
    double profile_start_ = 0;
    MPI_Comm comm_;
    int maxIter_;
    RealType tolerance_;
    bool verbose_;
    PrecondType precondType_;

    // Per-solve diagnostics (filled by solveImpl after Hypre returns).
    int    lastNumIters_ = 0;
    double lastFinalRes_ = 0.0;
    double last_absolute_residual_ = 0, last_rhs_norm_ = 0;
    double lastSolutionMax_ = 0.0;        // ||x||inf of the returned solution
    bool   nullSolutionReturned_ = false; // converged-but-x~0 (null-mode) flag

    HYPRE_Solver solver_;
    HYPRE_Solver precond_;

    HYPRE_IJMatrix A_hypre_;
    HYPRE_ParCSRMatrix parcsr_A_ = nullptr;

    // Optional SEPARATE preconditioner matrix K (the AMG-friendly Galerkin
    // stiffness). When set, BoomerAMG is built on parcsr_K_ while GMRES still
    // solves parcsr_A_ (the Gram operator A_fem). This is the Schur-complement
    // preconditioning fix: A_fem projects exactly but is AMG-hostile, K is
    // spectrally equivalent and AMG coarsens it well. Both share A_fem's
    // partition. nullptr -> classic single-matrix path (AMG on A itself).
    const Matrix*      precondMatrix_ = nullptr;
    HYPRE_IJMatrix     K_hypre_  = nullptr;
    HYPRE_ParCSRMatrix parcsr_K_ = nullptr;

    HYPRE_IJVector b_hypre_;
    HYPRE_IJVector x_hypre_;
    HYPRE_IJVector r_hypre_ = nullptr;
    HYPRE_ParVector par_b_ = nullptr;
    HYPRE_ParVector par_x_ = nullptr;
    HYPRE_ParVector par_r_ = nullptr;

    IndexType globalDofStart_;
    IndexType globalDofEnd_;

    // Device-resident map, copied when the setup is rebuilt.
    thrust::device_vector<HYPRE_BigInt> d_localToGlobalDof_;

    // GMRES(k) restart length.
    int kDim_;
    bool useFlexGmres_ = false;  // FlexGMRES (varying precond) vs plain GMRES
    int  pointBlock_   = 1;      // >1: systems/point-block BoomerAMG via SetNumFunctions
    int  coarseRelaxType_ = -1;
};

} // namespace fem
} // namespace mars
