// Taylor-Green vortex in a triply periodic box: the MARS periodic
// incompressible Navier-Stokes example (tutorial: docs/periodic_tgv_tutorial.md).
//
// main() follows the steps of the tutorial:
//   1. mesh and domain   a hex mesh of the box, distributed over the ranks by
//                        cornerstone; the box is periodic, so every rank also
//                        receives the elements across the opposite faces
//   2. periodic DOFs     each max-face node belongs to its min-face master;
//                        a periodic point is one unknown for every field
//   3. solver            the same NavierStokes solver as the Poiseuille and
//                        cavity examples; its matrices are assembled on the
//                        periodic unknowns and solved with BoomerAMG
//   4. time loop         BDF2 incremental projection; div u = 0 every step
//   5. output            kinetic energy against the viscous decay, VTU frames
//
// Initial condition, k = 2 pi / L on the box [lo, hi]^3 with L = hi - lo:
//   u =  V0 sin(kx) cos(ky) cos(kz)
//   v = -V0 cos(kx) sin(ky) cos(kz)
//   w =  0
//   p =  rho V0^2 / 16 (cos 2kx + cos 2ky)(cos 2kz + 2)
//
// Example (unit box, low Reynolds number, same result on any rank count):
//   mpirun -np 4 ./mars_tgv --mesh=cube16 --box-lo=0 --box-hi=1 --nu=0.05 --dt=1e-4 --num-steps=300

#include "backend/distributed/unstructured/fem/mars_navier_stokes.hpp"
#include "backend/distributed/unstructured/amr/mars_amr.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

using KeyType  = uint64_t;
using RealType = double;
using Solver   = mars::fem::NavierStokes<KeyType, RealType>;
using Amr      = mars::amr::AmrManager<mars::HexTag, KeyType, RealType>;
using Domain   = Solver::Domain;

// Nine device pointers, one per velocity gradient component.
template<typename T>
struct Components9
{
    const T* c[9];
};

// =============================================================================
// Command line
// =============================================================================

struct Options
{
    std::string mesh;
    std::string vtuPrefix;
    RealType V0          = 1;
    RealType rho         = 1;
    RealType nu          = RealType(1) / RealType(1600);
    RealType dt          = RealType(1e-3);
    RealType boxLo       = 0;
    RealType boxHi       = RealType(2 * M_PI);
    RealType tolerance   = RealType(1e-10);
    int numSteps         = 1000;
    int reportEvery      = 10;
    int vtuEvery         = 20;
    int maxIter          = 1000;
};

void printUsage()
{
    std::cout << "Usage: mars_tgv --mesh=DIR [options]\n"
                 "  --mesh=DIR            hex mesh of the periodic box (required)\n"
                 "  --box-lo=X --box-hi=X periodic box [lo, hi]^3 (default [0, 2 pi])\n"
                 "  --V0=X --rho=X        velocity amplitude, density (default 1, 1)\n"
                 "  --nu=X | --Re=X       kinematic viscosity, or nu = 1/Re (default 1/1600)\n"
                 "  --dt=X --num-steps=N  time step and step count (default 1e-3, 1000)\n"
                 "  --tol=X --max-iter=N  AMG-PCG tolerance and iteration cap (default 1e-10, 1000)\n"
                 "  --report-every=N      energy report interval (default 10)\n"
                 "  --vtu-output=PREFIX   write PVTU frames; --vtu-every=N (default 20)\n";
}

// Returns false when the program should stop (help, or a bad option).
bool parseOptions(int argc, char** argv, int rank, Options& o, int& exitCode)
{
    exitCode = 0;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        auto take = [&arg](const char* key, auto& out) {
            const std::string prefix = std::string("--") + key + "=";
            if (arg.rfind(prefix, 0) != 0) return false;
            std::istringstream(arg.substr(prefix.size())) >> out;
            return true;
        };
        RealType Re = 0;
        if (arg == "--help" || arg == "-h")
        {
            if (rank == 0) printUsage();
            return false;
        }
        if (take("Re", Re))
        {
            o.nu = RealType(1) / Re;
            continue;
        }
        bool known = take("mesh", o.mesh) || take("vtu-output", o.vtuPrefix) || take("V0", o.V0) ||
                     take("rho", o.rho) || take("nu", o.nu) || take("dt", o.dt) || take("box-lo", o.boxLo) ||
                     take("box-hi", o.boxHi) || take("tol", o.tolerance) || take("num-steps", o.numSteps) ||
                     take("report-every", o.reportEvery) || take("vtu-every", o.vtuEvery) ||
                     take("max-iter", o.maxIter);
        if (!known)
        {
            if (rank == 0)
            {
                std::cerr << "Unknown option: " << arg << "\n";
                printUsage();
            }
            exitCode = 1;
            return false;
        }
    }
    if (o.mesh.empty())
    {
        if (rank == 0) std::cerr << "Error: --mesh=DIR is required\n";
        exitCode = 1;
        return false;
    }
    return true;
}

// =============================================================================
// Initial condition
// =============================================================================

template<typename T>
__global__ void tgvInitialConditionKernel(const T* x, const T* y, const T* z, size_t n, T lo, T k, T V0, T rho,
                                          T* u, T* v, T* w, T* p)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    // One wavelength per box length, so opposite faces get identical values.
    T X = k * (x[i] - lo), Y = k * (y[i] - lo), Z = k * (z[i] - lo);
    u[i] = V0 * sin(X) * cos(Y) * cos(Z);
    v[i] = -V0 * cos(X) * sin(Y) * cos(Z);
    w[i] = T(0);
    p[i] = rho * V0 * V0 / T(16) * (cos(T(2) * X) + cos(T(2) * Y)) * (cos(T(2) * Z) + T(2));
}

void setInitialCondition(Solver& solver, const Domain& domain, const Options& o)
{
    const size_t n = domain.getNodeCount();
    const RealType k = RealType(2 * M_PI) / (o.boxHi - o.boxLo);
    tgvInitialConditionKernel<RealType><<<int((n + 255) / 256), 256>>>(
        domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(), n, o.boxLo, k, o.V0, o.rho,
        solver.u.data(), solver.v.data(), solver.w.data(), solver.p.data());
    cudaCheckError();
    solver.start();
}

// =============================================================================
// Validation output (kept apart from the solver steps in main)
//
// At low Reynolds number the vortex decays like the Stokes solution: the
// initial field is an eigenfunction of the Laplacian with eigenvalue -3 k^2,
// so KE(t) ~ KE(0) exp(-6 nu k^2 t). KE / KE_Stokes stays close to 1 on a
// resolved mesh (the nonlinear transfer is weak there); at Re = 1600 it does
// not, and the ratio then only shows the transition. div is max |D_F F| / M,
// the continuity of the face fluxes: roundoff level on every rank count.
// =============================================================================

struct EnergyReport
{
    RealType ke0      = 0;
    RealType decay    = 0;   // 6 nu k^2
    RealType maxDiv   = 0;
    int rank          = 0;

    EnergyReport(Solver& solver, const Options& o, int rank_)
        : rank(rank_)
    {
        const RealType k = RealType(2 * M_PI) / (o.boxHi - o.boxLo);
        ke0              = solver.kineticEnergy();
        decay            = 6 * o.nu * k * k;
    }

    void print(Solver& solver, int step, RealType t)
    {
        RealType ke  = solver.kineticEnergy();
        RealType div = solver.maxContinuity();
        maxDiv       = std::max(maxDiv, div);
        if (rank != 0) return;
        std::cout << "Step " << std::setw(6) << step << "  t=" << std::fixed << std::setprecision(5) << t
                  << std::scientific << std::setprecision(10) << "  KE=" << ke
                  << "  KE/KE_Stokes=" << std::fixed << std::setprecision(8) << ke / (ke0 * std::exp(-decay * t))
                  << std::scientific << std::setprecision(3) << "  div=" << div << "  amg(u,v,w,p)="
                  << solver.velocityIterations(0) << "/" << solver.velocityIterations(1) << "/"
                  << solver.velocityIterations(2) << "/" << solver.pressureIterations() << "\n"
                  << std::defaultfloat;
    }

    void summary(Solver& solver, RealType t)
    {
        RealType ke = solver.kineticEnergy();
        if (rank != 0) return;
        std::cout << std::scientific << std::setprecision(10) << "TGV final: t=" << t << "  KE=" << ke
                  << "  KE_Stokes=" << ke0 * std::exp(-decay * t) << std::fixed << std::setprecision(8)
                  << "  KE/KE_Stokes=" << ke / (ke0 * std::exp(-decay * t)) << std::scientific << std::setprecision(3)
                  << "  max_div=" << maxDiv << "\n"
                  << std::defaultfloat;
    }
};

// =============================================================================
// VTU output: velocity, pressure, vorticity magnitude
// =============================================================================

// |curl u| from the velocity gradients: g[3 c + d] = d u_c / d x_d.
template<typename T>
__global__ void vorticityKernel(size_t n, Components9<T> g, T* out)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    T ox = g.c[7][i] - g.c[5][i]; // dw/dy - dv/dz
    T oy = g.c[2][i] - g.c[6][i]; // du/dz - dw/dx
    T oz = g.c[3][i] - g.c[1][i]; // dv/dx - du/dy
    out[i] = sqrt(ox * ox + oy * oy + oz * oz);
}

struct FrameWriter
{
    std::unique_ptr<mars::fem::VTUParallelWriter<KeyType, RealType>> writer;
    cstone::DeviceVector<RealType> grad[9], omega;

    explicit FrameWriter(const std::string& prefix)
    {
        if (!prefix.empty()) writer = std::make_unique<mars::fem::VTUParallelWriter<KeyType, RealType>>(prefix);
    }

    void write(Solver& solver, const Domain& domain, int step, RealType t)
    {
        if (!writer) return;
        const cstone::DeviceVector<RealType>* velocity[3] = {&solver.u, &solver.v, &solver.w};
        for (int c = 0; c < 3; ++c)
            solver.gradient(*velocity[c], grad[3 * c], grad[3 * c + 1], grad[3 * c + 2]);
        const size_t n = domain.getNodeCount();
        omega.resize(n);
        Components9<RealType> g;
        for (int k = 0; k < 9; ++k)
            g.c[k] = grad[k].data();
        vorticityKernel<RealType><<<int((n + 255) / 256), 256>>>(n, g, omega.data());
        cudaCheckError();
        using FD = typename mars::fem::VTUParallelWriter<KeyType, RealType>::FieldDesc;
        std::vector<FD> fields = {{"u", FD::Kind::PointScalar, &solver.u, nullptr, nullptr},
                                  {"v", FD::Kind::PointScalar, &solver.v, nullptr, nullptr},
                                  {"w", FD::Kind::PointScalar, &solver.w, nullptr, nullptr},
                                  {"p", FD::Kind::PointScalar, &solver.p, nullptr, nullptr},
                                  {"omega", FD::Kind::PointScalar, &omega, nullptr, nullptr},
                                  {"velocity", FD::Kind::PointVector3, &solver.u, &solver.v, &solver.w}};
        writer->writeMultiFieldFrame(step, t, domain, fields);
    }
};

// Marks the nodes on the max faces. The solver's DofSpace finds each one's master by key,
// on whatever rank owns it, so no cross-rank pair table is needed.
void pairPeriodicNodes(const Domain& domain, mars::fem::PeriodicMap<KeyType, RealType>& map, const Options& o,
                      RealType faceEps)
{
    mars::fem::buildPeriodicMap<KeyType, RealType>(domain, map, o.boxLo, o.boxHi, o.boxLo, o.boxHi, o.boxLo, o.boxHi,
                                                   faceEps);
}

// =============================================================================
// main: the tutorial path
// =============================================================================

// Steps 1 to 5. The domain, the periodic map and the solver hold MPI
// resources, so they must be destroyed before MPI_Finalize: they live here.
int runTgv(const Options& opt, int rank, int numRanks)
{
    // 1. Mesh and domain. periodicAxesMask = 7 makes the cornerstone box
    //    periodic in x, y and z: node coordinates stay real, and each rank's
    //    halo also holds the elements on the other side of every periodic face.
    //    The AMR manager only loads the mesh here (maxLevels = 0): the solver
    //    does not constrain the hanging nodes that refinement would leave.
    Amr::Config amrConfig;
    amrConfig.maxLevels = 0;
    Amr amr(amrConfig);
    amr.initialize(opt.mesh, rank, numRanks, /*periodicAxesMask=*/7, opt.boxLo, opt.boxHi);
    Domain& domain = amr.domain();

    // 2. Periodic DOFs. Every node on a max face (x, y or z = hi) is marked; its
    //    master on the min faces has the same key with those axes moved to lo,
    //    on whatever rank owns it. The solver keeps one unknown per periodic point.
    mars::fem::PeriodicMap<KeyType, RealType> periodicMap;
    pairPeriodicNodes(domain, periodicMap, opt, RealType(1e-6) * (opt.boxHi - opt.boxLo));

    // 3. Solver. A periodic box has no boundary conditions and no openings; the
    //    periodic map makes every operator act on one unknown per periodic point.
    Solver::Params params;
    params.nu        = opt.nu;
    params.rho       = opt.rho;
    params.dt        = opt.dt;
    params.maxIter   = opt.maxIter;
    params.tolerance = opt.tolerance;
    Solver solver(domain, params, mars::fem::FreeNodes<RealType>{}, std::vector<mars::fem::Opening<RealType>>{},
                  &periodicMap);

    if (rank == 0)
        std::cout << "TGV: ranks=" << numRanks << "  periodic DOFs=" << solver.globalDofs()
                  << "  (the same on every rank count)\n"
                  << "     nu=" << opt.nu << " dt=" << opt.dt << " steps=" << opt.numSteps << "\n";

    setInitialCondition(solver, domain, opt);

    // 4. Time loop.
    EnergyReport report(solver, opt, rank);
    FrameWriter frames(opt.vtuPrefix);
    report.print(solver, 0, 0);
    frames.write(solver, domain, 0, 0);

    auto wallStart = std::chrono::steady_clock::now();
    for (int step = 1; step <= opt.numSteps; ++step)
    {
        if (!solver.step())
        {
            if (rank == 0) std::cerr << "Step " << step << ": a linear solve did not converge, stopping\n";
            return 1;
        }
        if (step == 1 && opt.numSteps > 1) solver.resetTiming();

        // 5. Output.
        RealType t = step * opt.dt;
        if (step % opt.reportEvery == 0 || step == opt.numSteps) report.print(solver, step, t);
        if (step % opt.vtuEvery == 0 || step == opt.numSteps) frames.write(solver, domain, step, t);
    }
    double wallMs =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - wallStart).count();

    report.summary(solver, opt.numSteps * opt.dt);
    solver.printTiming(opt.numSteps - 1); // step 1 (first-use allocations, BDF1) is left out
    if (rank == 0)
        std::cout << "Wall time " << std::fixed << std::setprecision(1) << wallMs << " ms, "
                  << wallMs / std::max(opt.numSteps, 1) << " ms/step\n";
    return 0;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    mars::abortAllRanksOnUncaughtException();
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    Options opt;
    int exitCode = 0;
    if (parseOptions(argc, argv, rank, opt, exitCode)) exitCode = runTgv(opt, rank, numRanks);
    MPI_Finalize();
    return exitCode;
}
