// BP3 with MFEM: the baseline for mars_marsir_bp --spec bp3 --mg --amg. Built from MFEM's
// examples/ex26p.cpp (order-refined GeometricMultigrid, partial assembly, Chebyshev
// smoothing, BoomerAMG on p = 1) and set to the same discrete problem as the MARS run:
//
// - the unit cube as ne^3 trilinear hexahedra, vertices moved by x += a sin(pi x) sin(pi y)
//   sin(pi z) (1,1,1) (a vertex-only mesh, so the elements stay trilinear);
// - H1 GLL nodal basis of degree p, q = p + 2 Gauss points per axis on every level;
// - b = A u with u = sin(pi x) sin(pi y) sin(pi z) at the nodes, zero on the boundary, so
//   the solve must return u; x starts at 0.
//
// Same space, quadrature and geometry, so A and b equal the MARS ones up to rounding.
// The preconditioners differ in details that are MFEM's own: Chebyshev on [0.3, 1.2]
// lmax of D^-1 A (MARS [1.1/15, 1.1] lmax), 10 power iterations (MARS 30), hypre's own
// coarsest level. The V-cycle has the same shape: degrees p -> 3 -> 1, order-2 Chebyshev
// before and after (2k+1 = 5 fine applies), BoomerAMG with --amg-cycles V-cycles (PMIS,
// ext+i, Pmax 4, theta 0.25, l1-Jacobi: MFEM's GPU defaults, also the MARS settings).
// CG stops on sqrt(<r,z>/<r0,z0>) like FCG, so compare the printed true residual
// ||B - A X|| / ||B|| (each free DoF once, the same quantity MARS prints).
//
// --amg-settings mars (default): Hypre and BoomerAMG as MARS sets them, so both codes run the
// same coarse solve: Hypre's own SpMV and SpGEMM kernels (MFEM leaves cuSPARSE's SpMV on CUDA),
// GPU-aware MPI (this Hypre build defaults to host copies), the device pool of
// fem::hypreInitializeOnce, and on the coarsest level 4 l1-Jacobi sweeps instead of a direct
// solve (cellwise::P1Amg). --amg-settings mfem: MFEM's and Hypre's defaults.
//
// The serial mesh is (ne / 2^r)^3, split into a Cartesian rank grid, then refined r times
// in parallel (each rank holds the serial mesh, so it stays at most 32^3).
//
// --solves k runs the solve k times from X = 0 and prints each time; the results reported are
// those of the last solve (MFEM allocates device memory on first use, so the first solve also
// pays for that). --profile brackets the last solve with cudaProfilerStart / Stop (nsys
// --capture-range=cudaProfilerApi).
//
// Run: mfem_bp3 --p 3 --ne 64 [--deform 0.05] [--tol 1e-10] [--maxit 2000] [--pc mg|jacobi]
//               [--amg-cycles 3] [--amg-settings mars|mfem] [--reps 20] [--solves 1] [--profile]
//               [--device cuda]
//               [--no-gpu-aware-mpi]

#include "mfem.hpp"
#ifdef MFEM_USE_CUDA
#include <cuda_profiler_api.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <vector>

using namespace mfem;

namespace {

double g_deform = 0.05;

void deform(const Vector& x, Vector& y)
{
    const double bump = g_deform * std::sin(M_PI * x(0)) * std::sin(M_PI * x(1)) * std::sin(M_PI * x(2));
    y.SetSize(3);
    for (int d = 0; d < 3; ++d) y(d) = x(d) + bump;
}

double exact(const Vector& x) { return std::sin(M_PI * x(0)) * std::sin(M_PI * x(1)) * std::sin(M_PI * x(2)); }

// Gauss-Legendre with p + 2 points per axis: exact to order 2 (p + 2) - 1.
const IntegrationRule& rule(int p) { return IntRules.Get(Geometry::CUBE, 2 * (p + 2) - 1); }

// ex26p's DiffusionMultigrid with the BP3 rule on every level, FULL assembly on the
// coarse level (on the device) and BoomerAMG alone as the coarse solve.
class Bp3Multigrid : public GeometricMultigrid
{
public:
    Bp3Multigrid(ParFiniteElementSpaceHierarchy& fespaces, Array<int>& ess_bdr, int cheb_order, int amg_cycles,
                 bool mars_amg)
        : GeometricMultigrid(fespaces, ess_bdr)
    {
        ParFiniteElementSpace& coarse = fespaces.GetFESpaceAtLevel(0);
        form(coarse, AssemblyLevel::FULL);
        // FULL assembly returns a new matrix owned by the handle (ex26p's HypreParMatrix&
        // overload would leave a dangling reference); take it over.
        OperatorPtr matrix(Operator::Hypre_ParCSR);
        bfs[0]->FormSystemMatrix(*essentialTrueDofs[0], matrix);
        matrix.SetOperatorOwner(false);
        coarse_matrix_ = matrix.As<HypreParMatrix>();
        amg_ = new HypreBoomerAMG(*coarse_matrix_);
        amg_->SetPrintLevel(0);
        amg_->SetMaxIter(amg_cycles);
        if (mars_amg) {
            const HYPRE_Solver s = *amg_;
            HYPRE_BoomerAMGSetRelaxOrder(s, 0);
            HYPRE_BoomerAMGSetMaxCoarseSize(s, 9);
            HYPRE_BoomerAMGSetCycleRelaxType(s, 18, 3);
            HYPRE_BoomerAMGSetCycleNumSweeps(s, 4, 3);
        }
        // BoomerAMG sets up on its first Mult; do it here so the setup time holds it.
        Vector zero(coarse.GetTrueVSize()), out(coarse.GetTrueVSize());
        zero.UseDevice(true);
        out.UseDevice(true);
        zero = 0.0;
        amg_->Setup(zero, out);
        AddLevel(coarse_matrix_, amg_, true, true);

        for (int level = 1; level < fespaces.GetNumLevels(); ++level) {
            ParFiniteElementSpace& fes = fespaces.GetFESpaceAtLevel(level);
            const Array<int>& ess = *essentialTrueDofs[level];
            form(fes, AssemblyLevel::PARTIAL);
            OperatorPtr op;
            op.SetType(Operator::ANY_TYPE);
            bfs.Last()->FormSystemMatrix(ess, op);
            op.SetOperatorOwner(false);
            Vector diag(fes.GetTrueVSize());
            bfs.Last()->AssembleDiagonal(diag);
            Solver* smoother = new OperatorChebyshevSmoother(*op, diag, ess, cheb_order, fes.GetParMesh()->GetComm());
            AddLevel(op.Ptr(), smoother, true, true);
        }
    }

private:
    void form(ParFiniteElementSpace& fes, AssemblyLevel level)
    {
        auto* f = new ParBilinearForm(&fes);
        f->SetAssemblyLevel(level);
        if (level == AssemblyLevel::FULL) f->EnableSparseMatrixSorting(Device::IsEnabled());
        auto* integ = new DiffusionIntegrator(one_);
        integ->SetIntRule(&rule(fes.GetMaxElementOrder()));
        f->AddDomainIntegrator(integ);
        f->Assemble();
        bfs.Append(f);
    }

    ConstantCoefficient one_{1.0};
    HypreParMatrix* coarse_matrix_ = nullptr;   // owned by the Multigrid (AddLevel)
    HypreBoomerAMG* amg_ = nullptr;             // owned by the Multigrid (AddLevel)
};

// sqrt(<r,z>/<r0,z0>) after each iteration (CG passes <r,z> itself), as MARS prints it.
class History : public IterativeSolverMonitor
{
public:
    void MonitorResidual(int it, real_t norm, const Vector&, bool final) override
    {
        if (final) return;
        if (it == 0) {
            nom0_ = norm;
            rel.clear();
        } else {
            rel.push_back(std::sqrt(norm / nom0_));
        }
    }
    std::vector<double> rel;

private:
    double nom0_ = 1.0;
};

double all_max(double v)
{
    MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return v;
}

}  // namespace

int main(int argc, char* argv[])
{
    Mpi::Init(argc, argv);
    const int rank = Mpi::WorldRank(), ranks = Mpi::WorldSize();
    const bool root = rank == 0;
    Hypre::Init();

    int p = 3, ne = 32, maxit = 2000, cheb_order = 2, amg_cycles = 3, reps = 20, solves = 1;
    bool profile = false;
    double tol = 1e-10;
    const char* pc = "mg";
    const char* amg_settings = "mars";
    const char* device_config = "cuda";
    bool gpu_aware_mpi = true;
    OptionsParser args(argc, argv);
    args.AddOption(&p, "-p", "--p", "Degree of the finest level: 1, 3 or 7.");
    args.AddOption(&ne, "-ne", "--ne", "Elements per axis.");
    args.AddOption(&g_deform, "-a", "--deform", "Amplitude of the vertex deformation.");
    args.AddOption(&tol, "-t", "--tol", "CG relative tolerance on sqrt(<r,z>).");
    args.AddOption(&maxit, "-it", "--maxit", "CG iteration limit.");
    args.AddOption(&pc, "-pc", "--pc", "Preconditioner: mg (V-cycle) or jacobi.");
    args.AddOption(&cheb_order, "-cheb", "--cheb-order", "Chebyshev order on each multigrid level.");
    args.AddOption(&amg_cycles, "-amg", "--amg-cycles", "BoomerAMG V-cycles per coarse solve.");
    args.AddOption(&amg_settings, "-amgs", "--amg-settings", "Hypre and BoomerAMG settings: mars or mfem.");
    args.AddOption(&reps, "-reps", "--reps", "Operator applications to time.");
    args.AddOption(&solves, "-s", "--solves", "Solves to run; the last is reported.");
    args.AddOption(&profile, "-prof", "--profile", "-no-prof", "--no-profile",
                   "cudaProfilerStart / Stop around the last solve.");
    args.AddOption(&device_config, "-d", "--device", "MFEM device configuration.");
    args.AddOption(&gpu_aware_mpi, "-gam", "--gpu-aware-mpi", "-no-gam", "--no-gpu-aware-mpi",
                   "Pass device buffers to MPI.");
    args.Parse();
    const bool mg = !std::strcmp(pc, "mg");
    const bool mars_amg = !std::strcmp(amg_settings, "mars");
    if (!args.Good() || (p != 1 && p != 3 && p != 7) || (!mg && std::strcmp(pc, "jacobi")) || solves < 1 ||
        (!mars_amg && std::strcmp(amg_settings, "mfem"))) {
        if (root) args.PrintUsage(std::cout);
        return 1;
    }

    if (mars_amg) {
        // Before Hypre's first device allocation (the pool size).
        HYPRE_SetSpMVUseVendor(0);
        HYPRE_SetSpGemmUseVendor(0);
        HYPRE_SetGpuAwareMPI(1);
        int max_bin = 0;
        while ((size_t(2) << max_bin) <= (size_t(128) << 20)) ++max_bin;
        HYPRE_SetGPUMemoryPoolSize(2, 1, max_bin, size_t(2048) << 20);
    }

    // One GPU per rank, by the rank's index on its node.
    int local = 0;
    {
        MPI_Comm node;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node);
        MPI_Comm_rank(node, &local);
        MPI_Comm_free(&node);
    }
    int gpus = 1;
#ifdef MFEM_USE_CUDA
    gpus = std::max(1, CuGetDeviceCount());
#endif
    Device device(device_config, local % gpus);
    if (gpu_aware_mpi) Device::SetGPUAwareMPI(true);
    if (root) device.Print();

    // Serial mesh of at most 32^3 elements, r parallel refinements to ne^3.
    int r = 0;
    while ((ne >> r) > 32 && (ne >> r) % 2 == 0) ++r;
    const int ne0 = ne >> r;
    int dims[3] = {0, 0, 0};
    MPI_Dims_create(ranks, 3, dims);
    if (ne0 < dims[0] || ne0 < dims[1] || ne0 < dims[2]) {
        if (root) printf("serial mesh %d^3 has fewer elements per axis than ranks (%d x %d x %d)\n", ne0, dims[0],
                         dims[1], dims[2]);
        return 1;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    tic_toc.Clear();
    tic_toc.Start();
    std::unique_ptr<ParMesh> pmesh;
    {
        Mesh serial = Mesh::MakeCartesian3D(ne0, ne0, ne0, Element::HEXAHEDRON, 1.0, 1.0, 1.0, false);
        // Partition by element centers before the deformation, so each rank owns a box.
        std::unique_ptr<int[]> part(serial.CartesianPartitioning(dims));
        pmesh.reset(new ParMesh(MPI_COMM_WORLD, serial, part.get()));
    }
    for (int l = 0; l < r; ++l) pmesh->UniformRefinement();
    pmesh->Transform(deform);

    // Degrees 1 [-> 3] [-> 7]; the hierarchy owns the mesh and the spaces.
    std::vector<std::unique_ptr<FiniteElementCollection>> collections;
    collections.emplace_back(new H1_FECollection(1, 3));
    auto* coarse_fes = new ParFiniteElementSpace(pmesh.get(), collections.back().get());
    ParFiniteElementSpaceHierarchy fespaces(pmesh.release(), coarse_fes, true, true);
    for (int deg : {3, 7}) {
        if (deg > p) break;
        collections.emplace_back(new H1_FECollection(deg, 3));
        fespaces.AddOrderRefinedLevel(collections.back().get());
    }
    ParFiniteElementSpace& fes = fespaces.GetFinestFESpace();
    const HYPRE_BigInt unique = fes.GlobalTrueVSize();

    Array<int> ess_bdr(fes.GetParMesh()->bdr_attributes.Max());
    ess_bdr = 1;
    Array<int> ess;
    fes.GetEssentialTrueDofs(ess_bdr, ess);

    // The finest operator with the boundary rows eliminated, and the preconditioner.
    std::unique_ptr<Bp3Multigrid> multigrid;
    std::unique_ptr<ParBilinearForm> fine_form;
    std::unique_ptr<OperatorJacobiSmoother> jacobi;
    ConstantCoefficient one(1.0);
    OperatorPtr fine_op;
    Operator* A = nullptr;
    Solver* M = nullptr;
    if (mg) {
        multigrid.reset(new Bp3Multigrid(fespaces, ess_bdr, cheb_order, amg_cycles, mars_amg));
        multigrid->SetCycleType(Multigrid::CycleType::VCYCLE, 1, 1);
        A = multigrid->GetOperatorAtFinestLevel();
        M = multigrid.get();
    } else {
        fine_form.reset(new ParBilinearForm(&fes));
        fine_form->SetAssemblyLevel(AssemblyLevel::PARTIAL);
        auto* integ = new DiffusionIntegrator(one);
        integ->SetIntRule(&rule(p));
        fine_form->AddDomainIntegrator(integ);
        fine_form->Assemble();
        fine_op.SetType(Operator::ANY_TYPE);
        fine_form->FormSystemMatrix(ess, fine_op);
        A = fine_op.Ptr();
        Vector diag(fes.GetTrueVSize());
        fine_form->AssembleDiagonal(diag);
        jacobi.reset(new OperatorJacobiSmoother(diag, ess));
        M = jacobi.get();
    }
    MFEM_DEVICE_SYNC;
    MPI_Barrier(MPI_COMM_WORLD);
    tic_toc.Stop();
    const double setup_ms = all_max(tic_toc.RealTime() * 1e3);

    // u at the nodes (the trilinear image of the GLL points), exactly zero on the boundary.
    ParGridFunction ugf(&fes);
    FunctionCoefficient u_coef(exact);
    ugf.ProjectCoefficient(u_coef);
    const int n = fes.GetTrueVSize();
    Vector U(n), B(n), X(n), Y(n);
    U.UseDevice(true);
    B.UseDevice(true);
    X.UseDevice(true);
    Y.UseDevice(true);
    ugf.GetTrueDofs(U);
    U.SetSubVector(ess, 0.0);
    A->Mult(U, B);
    X = 0.0;

    CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(tol);
    cg.SetAbsTol(0.0);
    cg.SetMaxIter(maxit);
    cg.SetPrintLevel(IterativeSolver::PrintLevel().None());
    cg.SetOperator(*A);
    cg.SetPreconditioner(*M);
    History history;
    cg.SetMonitor(history);
    std::vector<double> solve_times;
    for (int k = 0; k < solves; ++k) {
        const bool last = k + 1 == solves;
        X = 0.0;
        MFEM_DEVICE_SYNC;
        MPI_Barrier(MPI_COMM_WORLD);
#ifdef MFEM_USE_CUDA
        if (last && profile) cudaProfilerStart();
#endif
        tic_toc.Clear();
        tic_toc.Start();
        cg.Mult(B, X);
        MFEM_DEVICE_SYNC;
        tic_toc.Stop();
#ifdef MFEM_USE_CUDA
        if (last && profile) cudaProfilerStop();
#endif
        solve_times.push_back(all_max(tic_toc.RealTime() * 1e3));
    }
    const double solve_ms = solve_times.back();
    const int iterations = cg.GetNumIterations();

    // True residual over the true DoFs (boundary rows: B = X = 0 there).
    A->Mult(X, Y);
    Y -= B;
    Y.Neg();
    const double true_rel = std::sqrt(InnerProduct(MPI_COMM_WORLD, Y, Y) / InnerProduct(MPI_COMM_WORLD, B, B));
    Y = X;
    Y -= U;
    const double err = GlobalLpNorm(infinity(), Y.Normlinf(), MPI_COMM_WORLD);
    const double umax = GlobalLpNorm(infinity(), U.Normlinf(), MPI_COMM_WORLD);

    // The assembled operator: E-vector gather, element kernel, scatter, the shared-DoF exchange.
    A->Mult(U, Y);
    MFEM_DEVICE_SYNC;
    MPI_Barrier(MPI_COMM_WORLD);
    tic_toc.Clear();
    tic_toc.Start();
    for (int i = 0; i < reps; ++i) A->Mult(U, Y);
    MFEM_DEVICE_SYNC;
    tic_toc.Stop();
    const double op_ms = all_max(tic_toc.RealTime() * 1e3 / reps);

    const bool pass = cg.GetConverged() && err <= 1e-8 * umax;
    if (root) {
        printf("MFEM %s BP3, p=%d, q=%d, %d^3 elements (%lld), %lld unique DoFs, deform %.3f\n", MFEM_VERSION_STRING,
               p, p + 2, ne, (long long)ne * ne * ne, (long long)unique, g_deform);
        printf("  %d ranks (%d x %d x %d), serial mesh %d^3 refined %d times, GPU-aware MPI %s\n", ranks, dims[0],
               dims[1], dims[2], ne0, r, gpu_aware_mpi ? "on" : "off");
        if (mg)
            printf("  multigrid p=%d%s%s, Chebyshev order %d V(1,1), coarse BoomerAMG x%d; setup %.1f ms\n", p,
                   p == 7 ? " -> 3" : "", p > 1 ? " -> 1" : "", cheb_order, amg_cycles, setup_ms);
        else
            printf("  Jacobi; setup %.1f ms\n", setup_ms);
        printf("  Hypre: %s\n", mars_amg ? "MARS's settings (SpMV and SpGEMM hypre, GPU-aware MPI, pool 128 MiB "
                                            "blocks / 2 GiB; coarsest level 4 l1-Jacobi sweeps)"
                                          : "MFEM's and Hypre's defaults");
        printf("  CG + %s: %d iterations, sqrt(<r,z>/<r0,z0>) = %.3e, true residual %.3e\n", mg ? "V-cycle" : "Jacobi",
               iterations, cg.GetFinalRelNorm(), true_rel);
        for (int k : {1, 5, 10, 20, 50, 100})
            if (k <= (int)history.rel.size()) printf("  it %3d  %.6e\n", k, history.rel[k - 1]);
        printf("  max |x - u| / max |u| = %.3e  %s\n", err / umax, pass ? "PASS" : "FAIL");
        printf("  operator     %8.3f ms  %6.2f GDoF/s unique  (%s)\n", op_ms, unique / (op_ms * 1e6),
               mg && p == 1 ? "assembled matrix, HypreParMatrix" : "partial assembly, the assembled operator");
        if (solves > 1) {
            printf("  solves");
            for (double t : solve_times) printf(" %.1f", t);
            printf(" ms (the last is reported)\n");
        }
        printf("  CG iteration %8.3f ms  (solve %.1f ms, + setup %.0f ms)\n", iterations ? solve_ms / iterations : 0.0,
               solve_ms, setup_ms);
    }
    return pass ? 0 : 1;
}
