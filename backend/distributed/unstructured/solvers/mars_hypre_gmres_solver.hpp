#pragma once

// All shared utilities (HypreInitGuard, countValidPerRowKernel,
// compactGlobalCsrKernel, fillGlobalRowIndicesKernel, IsNonFinite, the
// HYPRE_/thrust/mpi headers, the namespace declarations) come from the PCG
// header. Including it here means we don't redefine those symbols.
#include "mars_hypre_pcg_solver.hpp"

namespace mars {
namespace fem {


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
            HYPRE_BoomerAMGSetRelaxOrder(precond_, amgRelaxOrder);
            HYPRE_BoomerAMGSetKeepTranspose(precond_, 1);
            HYPRE_BoomerAMGSetStrongThreshold(precond_, amgStrong);
            HYPRE_BoomerAMGSetPMaxElmts(precond_, amgPMax);
            HYPRE_BoomerAMGSetAggNumLevels(precond_, amgAgg);
            HYPRE_BoomerAMGSetNumSweeps(precond_, amgSweeps);
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
        } else if (precondType_ == JACOBI) {
            if (verbose_ && rank == 0) std::cout << "Using Jacobi preconditioner" << std::endl;
            precond_ = (HYPRE_Solver) parcsr_A_;
        }

        const char* krylovName = useFlexGmres_ ? "FlexGMRES" : "GMRES";
        if (verbose_ && rank == 0) std::cout << "Creating " << krylovName << " solver..." << std::endl;
        if (useFlexGmres_) HYPRE_ParCSRFlexGMRESCreate(comm_, &solver_);
        else               HYPRE_ParCSRGMRESCreate(comm_, &solver_);
        if (!solver_) {
            std::cerr << "Failed to create Hypre " << krylovName << " solver" << std::endl;
            return false;
        }
        // FlexGMRES and GMRES share the HYPRE_GMRES* setter/getter symbols
        // (FlexGMRES is a GMRES variant in Hypre); only Create/Setup/Solve/
        // Destroy and SetPrecond differ. So the Set* calls below are common.
        HYPRE_GMRESSetMaxIter(solver_, maxIter_);
        HYPRE_GMRESSetTol(solver_, tolerance_);
        HYPRE_GMRESSetKDim(solver_, kDim_);  // restart length
        // Floor on iterations. For the near-singular DDT pressure operator (one
        // constant null mode, single pin), BoomerAMG's first V-cycle can map the
        // initial residual almost entirely into the null space, dropping the
        // RELATIVE residual ||r||/||b|| below tol in ~2 iters while the actual
        // solution x is still ~0 (the false x=0 "convergence"). A small minimum
        // iteration count forces GMRES to keep building the Krylov space past
        // that first deceptive cycle so a real x emerges. Env-overridable.
        {
            int minIt = getEnvInt("MARS_HYPRE_MINITER", 3);
            if (minIt > 0) HYPRE_GMRESSetMinIter(solver_, minIt);
        }
        // Absolute residual floor. A pure relative test can be satisfied by a
        // tiny ||b|| whose mass sits in the null mode; pairing it with an
        // absolute tol means convergence requires the TRUE residual to be small,
        // not just relatively small versus a near-null RHS. Default 0 (disabled)
        // keeps legacy behavior; set MARS_HYPRE_ABSTOL>0 to engage.
        {
            double absTol = getEnvDouble("MARS_HYPRE_ABSTOL", 0.0);
            if (absTol > 0.0) HYPRE_GMRESSetAbsoluteTol(solver_, absTol);
        }
        // print level: 0 silent, 2 per-iter residuals. Env MARS_HYPRE_VERBOSE=1.
        {
            const char* ev = std::getenv("MARS_HYPRE_VERBOSE");
            int gmresPrint = (verbose_ || (ev && std::string(ev) != "0")) ? 2 : 0;
            HYPRE_GMRESSetPrintLevel(solver_, gmresPrint);
        }

        if (precondType_ == BOOMERAMG && precond_) {
            if (verbose_ && rank == 0) std::cout << "Setting BoomerAMG preconditioner..." << std::endl;
            if (useFlexGmres_)
                HYPRE_FlexGMRESSetPrecond(solver_,
                                   (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSolve,
                                   (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSetup,
                                   precond_);
            else
                HYPRE_GMRESSetPrecond(solver_,
                                   (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSolve,
                                   (HYPRE_Int (*)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector))HYPRE_BoomerAMGSetup,
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
        MPI_Barrier(comm_);
        HYPRE_Int setup_err = useFlexGmres_
            ? HYPRE_ParCSRFlexGMRESSetup(solver_, parcsr_A_, par_b_, par_x_)
            : HYPRE_ParCSRGMRESSetup(solver_, parcsr_A_, par_b_, par_x_);
        if (setup_err != 0 && rank == 0) {
            std::cerr << "[HypreGMRES] Setup returned error " << setup_err
                      << " (HYPRE_GetError=" << HYPRE_GetError() << ")\n";
            char errbuf[256];
            HYPRE_DescribeError(HYPRE_GetError(), errbuf);
            std::cerr << "[HypreGMRES] " << errbuf << "\n";
            HYPRE_ClearAllErrors();
        }
        if (verbose_ && rank == 0) std::cout << "GMRES setup complete, starting solve..." << std::endl;

        MPI_Barrier(comm_);
        HYPRE_Int solve_err = useFlexGmres_
            ? HYPRE_ParCSRFlexGMRESSolve(solver_, parcsr_A_, par_b_, par_x_)
            : HYPRE_ParCSRGMRESSolve(solver_, parcsr_A_, par_b_, par_x_);
        if (solve_err != 0 && rank == 0) {
            std::cerr << "[HypreGMRES] Solve returned error " << solve_err
                      << " (HYPRE_GetError=" << HYPRE_GetError() << ")\n";
            char errbuf[256];
            HYPRE_DescribeError(HYPRE_GetError(), errbuf);
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

        int    num_iterations = 0;
        double final_res_norm = 0.0;
        HYPRE_GMRESGetNumIterations(solver_, &num_iterations);
        HYPRE_GMRESGetFinalRelativeResidualNorm(solver_, &final_res_norm);
        lastNumIters_ = num_iterations;
        lastFinalRes_ = final_res_norm;

        // False-convergence guard against the near-singular DDT operator. When
        // BoomerAMG collapses the first cycle onto the constant null mode,
        // FlexGMRES can report converged=true with a tiny relative residual while
        // the returned x is essentially zero. Detect that here: compute ||x||inf
        // and ||b||inf (global), and if Hypre claims convergence but x is ~0 while
        // b is not, the "solution" is the null vector -- report NON-converged so
        // the caller does not scatter a zero phi as if it were a real projection.
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
            MPI_Allreduce(&localXmax, &gXmax, 1, MPI_DOUBLE, MPI_MAX, comm_);
            MPI_Allreduce(&localBmax, &gBmax, 1, MPI_DOUBLE, MPI_MAX, comm_);
            lastSolutionMax_ = gXmax;
            // x is "null" if it is many orders below b. The default 1e-12 only
            // fires on an essentially-exact x=0 (the observed failure mode), so
            // it never rejects a physically meaningful phi. Raise the ratio via
            // MARS_HYPRE_NULLX_RATIO if a borderline near-null x appears.
            double nullRatio = getEnvDouble("MARS_HYPRE_NULLX_RATIO", 1e-12);
            nullSolutionReturned_ = (gBmax > 0.0 && gXmax < nullRatio * gBmax);
            if (nullSolutionReturned_ && rank == 0) {
                std::cout << "[HypreGMRES] WARNING: solver reported converged (iters="
                          << num_iterations << ", rel_res=" << final_res_norm
                          << ") but |x|inf=" << gXmax << " << |b|inf=" << gBmax
                          << " -- this is the null-mode false convergence; "
                          << "reporting NOT converged.\n";
            }
            // UPPER-BOUND reject: a |x| that is HUGE relative to |b| is an
            // under-resolved near-null-space solution (the observed 1e7..1e9 phi
            // from a stalled solve on a poorly-conditioned operator). Its small
            // residual can slip a loose acceptance gate, but its gradient blows up
            // the corrector. Reject it. Default ceiling 1e6 (an O(1)-conditioned
            // pressure solve has |phi| ~ |b|/diag, not 1e6x|b|); MARS_HYPRE_MAXX_RATIO.
            double maxXRatio = getEnvDouble("MARS_HYPRE_MAXX_RATIO", 1e6);
            if (gBmax > 0.0 && gXmax > maxXRatio * gBmax) {
                nullSolutionReturned_ = true;
                if (rank == 0)
                    std::cout << "[HypreGMRES] WARNING: |x|inf=" << gXmax
                              << " >> |b|inf=" << gBmax << " (ratio>" << maxXRatio
                              << ") -- under-resolved/null-contaminated solution; "
                              << "reporting NOT converged.\n";
            }
        }

        if (verbose_) {
            std::cout << "Hypre GMRES converged in " << num_iterations
                      << " iterations, final residual: " << final_res_norm << std::endl;
        }

        // Real convergence requires BOTH the residual test AND a non-null x.
        return (final_res_norm < tolerance_) && !nullSolutionReturned_;
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

    // Per-solve diagnostics (filled by solveImpl after Hypre returns).
    int    lastNumIters_ = 0;
    double lastFinalRes_ = 0.0;
    double lastSolutionMax_ = 0.0;        // ||x||inf of the returned solution
    bool   nullSolutionReturned_ = false; // converged-but-x~0 (null-mode) flag

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

    // GMRES(k) restart length.
    int kDim_;
    bool useFlexGmres_ = false;  // FlexGMRES (varying precond) vs plain GMRES
};

} // namespace fem
} // namespace mars
