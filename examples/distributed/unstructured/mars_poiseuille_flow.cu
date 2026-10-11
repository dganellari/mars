// Poiseuille flow in a plane channel: the MARS incompressible Navier-Stokes
// example (tutorial: docs/poiseuille_tutorial.md).
//
// Uniform inflow U at x = xmin, no-slip walls at y = ymin and y = ymax, and a
// pressure outlet at x = xmax. Downstream of the entrance the flow develops into
// the parabola
//   u(y) = 1.5 U (1 - eta^2),   eta = (2 y - ymin - ymax) / H,
// held by the pressure gradient -dp/dx = 12 rho nu U / H^2.
//
// main() follows the steps of the tutorial:
//   1. mesh and domain   read a mesh, or generate the channel on every rank;
//                        cornerstone distributes the elements over the ranks
//   2. solver            the boundary conditions as a function of the node
//                        position, and the inlet and outlet; the solver assembles
//                        the viscous and pressure matrices and builds BoomerAMG
//   3. time loop         BDF2 projection steps, D u = 0 after every step
//   4. result            timing, and the profile against the parabola
//                        (mars_poiseuille_validation.hpp)
//
// Examples:
//   mpirun -np 4 ./mars_poiseuille_flow --mesh=poiseuille_hex_14k_elem.e --num-steps=1500
//   mpirun -np 64 ./mars_poiseuille_flow --cells=8192,2048 --num-steps=20      (scaling)

#include "mars_poiseuille_validation.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"

#include <thrust/fill.h>

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <string>
#include <vector>

using poiseuille::Domain;
using poiseuille::KeyType;
using poiseuille::RealType;
using poiseuille::Solver;

// =============================================================================
// Command line
// =============================================================================

struct Options
{
    std::string mesh;
    std::array<size_t, 2> cells{0, 0}; // generated channel [0,10] x [0,1] x [0,0.06], one cell in z
    double yGrading = 0;               // clusters the generated y spacing at the walls, in [0, 1)
    double inflow   = 1;               // inlet velocity U
    Solver::Params params;
    int numSteps    = 1000;
    int reportEvery = 1;
    std::string vtuPrefix;
    int vtuEvery   = 50;
    int bucketSize = 64;
    poiseuille::Options validation;
};

void printUsage()
{
    std::cout << "Usage: mars_poiseuille_flow --mesh=FILE | --cells=NX,NY [options]\n"
                 "  --mesh=FILE           one layer of axis-aligned hexes (.e / .exo / binary mesh directory)\n"
                 "  --cells=NX,NY         generate the channel [0,10] x [0,1] x [0,0.06] instead\n"
                 "  --y-grading=S         cluster the generated y spacing at the walls, S in [0,1) (default 0)\n"
                 "  --uinf=X --rho=X --nu=X  inflow velocity, density, viscosity (default 1, 1, 0.01)\n"
                 "  --dt=X --num-steps=N  time step and step count (default 0.01, 1000)\n"
                 "  --bdf1                first-order time stepping (default BDF2)\n"
                 "  --tol=X --max-iter=N  relative AMG-PCG tolerance and iteration cap (default 1e-10, 1000)\n"
                 "  --report-every=N      progress line interval (default 1)\n"
                 "  --vtu-output=PREFIX   write PVTU frames; --vtu-every=N (default 50)\n"
                 "  --block-size=N --bucket-size=N  CUDA block size, cornerstone bucket size (default 256, 64)\n";
    poiseuille::printOptions();
}

// Returns false when the program should stop (help, or a bad option).
bool parseOptions(int argc, char** argv, int rank, Options& o, int& exitCode)
{
    exitCode = 0;
    bool bad = false;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        auto take = [&arg, &bad](const char* key, auto& out) {
            const std::string prefix = std::string("--") + key + "=";
            if (arg.rfind(prefix, 0) != 0) return false;
            std::istringstream in(arg.substr(prefix.size()));
            bad = bad || !(in >> out) || !in.eof();
            return true;
        };
        if (arg == "--help" || arg == "-h")
        {
            if (rank == 0) printUsage();
            return false;
        }
        std::string cells;
        if (take("cells", cells))
        {
            std::replace(cells.begin(), cells.end(), ',', ' ');
            std::istringstream in(cells);
            bad = bad || !(in >> o.cells[0] >> o.cells[1]);
            continue;
        }
        if (arg == "--bdf1")
        {
            o.params.bdf2 = false;
            continue;
        }
        Solver::Params& p = o.params;
        bool known = take("mesh", o.mesh) || take("y-grading", o.yGrading) || take("uinf", o.inflow) ||
                     take("rho", p.rho) || take("nu", p.nu) || take("dt", p.dt) || take("num-steps", o.numSteps) ||
                     take("tol", p.tolerance) || take("max-iter", p.maxIter) || take("report-every", o.reportEvery) ||
                     take("vtu-output", o.vtuPrefix) || take("vtu-every", o.vtuEvery) ||
                     take("block-size", p.blockSize) || take("bucket-size", o.bucketSize) ||
                     poiseuille::parseOption(arg, o.validation, bad);
        if (!known)
        {
            if (rank == 0) std::cerr << "Unknown option: " << arg << "\n";
            bad = true;
        }
    }
    const Solver::Params& p = o.params;
    bool generated          = o.cells[0] > 0 && o.cells[1] > 0;
    bad = bad || o.mesh.empty() == !generated || !(o.yGrading >= 0 && o.yGrading < 1) || !(o.inflow > 0) ||
          !(p.rho > 0) || !(p.nu > 0) || !(p.dt > 0) || !(p.tolerance > 0) || p.maxIter <= 0 || o.numSteps <= 0 ||
          o.reportEvery <= 0 || o.vtuEvery <= 0;
    if (bad)
    {
        if (rank == 0) printUsage();
        exitCode = 1;
        return false;
    }
    return true;
}

// =============================================================================
// 1. Mesh and domain
// =============================================================================

// storeOriginalCoords keeps the exact node coordinates; otherwise they are
// decoded from the space-filling-curve keys.
std::unique_ptr<Domain> makeDomain(const Options& o, int rank, int numRanks)
{
    if (!o.mesh.empty()) return std::make_unique<Domain>(o.mesh, rank, numRanks, true, o.bucketSize);

    // Every rank generates its own brick of the channel: no mesh file, so no
    // rank reads or broadcasts the whole mesh.
    [[maybe_unused]] auto [nodes, elements, x, y, z, conn] = mars::generateBoxElementPartition<RealType, KeyType>(
        {o.cells[0], o.cells[1], 1}, {0.0, 0.0, 0.0}, {10.0, 1.0, 0.06}, o.yGrading, rank, numRanks);
    Domain::HostCoordsTuple coords{std::move(x), std::move(y), std::move(z)};
    Domain::HostConnectivityTuple connectivity{std::move(conn[0]), std::move(conn[1]), std::move(conn[2]),
                                               std::move(conn[3]), std::move(conn[4]), std::move(conn[5]),
                                               std::move(conn[6]), std::move(conn[7])};
    return std::make_unique<Domain>(coords, connectivity, rank, numRanks, o.bucketSize, true);
}

// =============================================================================
// Output
// =============================================================================

// PREFIX.pvd with the point fields u, v and p (scripts/render_poiseuille.py).
struct FrameWriter
{
    std::unique_ptr<mars::fem::VTUParallelWriter<KeyType, RealType>> writer;

    explicit FrameWriter(const std::string& prefix)
    {
        if (!prefix.empty()) writer = std::make_unique<mars::fem::VTUParallelWriter<KeyType, RealType>>(prefix);
    }
    void write(Solver& s, const Domain& domain, int step, double t)
    {
        if (!writer) return;
        using FD = mars::fem::VTUParallelWriter<KeyType, RealType>::FieldDesc;
        std::vector<FD> fields{{"u", FD::Kind::PointScalar, &s.u, nullptr, nullptr},
                               {"v", FD::Kind::PointScalar, &s.v, nullptr, nullptr},
                               {"p", FD::Kind::PointScalar, &s.p, nullptr, nullptr}};
        writer->writeMultiFieldFrame(step, t, domain, fields);
    }
};

// =============================================================================
// main: the tutorial path
// =============================================================================

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
    if (!parseOptions(argc, argv, rank, opt, exitCode))
    {
        MPI_Finalize();
        return exitCode;
    }

    // 1. Mesh and domain.
    auto domain = makeDomain(opt, rank, numRanks);
    {
        // 2. Boundary conditions and solver. The conditions are a function of the
        //    node position, evaluated on the GPU for every node: walls on the y
        //    faces (they win at the corners), inflow U at x = xmin, p = 0 at x = xmax.
        //    The z faces need nothing: planar flow keeps w = 0 there. Fluid enters
        //    and leaves through the two openings.
        //    Hypre objects live inside this scope: they must be destroyed before MPI_Finalize.
        const auto box     = mars::fem::boundingBox(*domain);
        const RealType eps = RealType(1e-5) * std::max(RealType(1), box.hi[1] - box.lo[1]);
        const RealType U   = opt.inflow;
        auto channel = [box, eps, U] __device__(RealType x, RealType y, RealType) -> mars::fem::NodeCondition<RealType> {
            mars::fem::NodeCondition<RealType> node;
            bool wall          = fabs(y - box.lo[1]) < eps || fabs(y - box.hi[1]) < eps;
            bool inlet         = fabs(x - box.lo[0]) < eps;
            node.velocityFixed = wall || inlet;
            node.velocity[0]   = inlet && !wall ? U : RealType(0);
            node.pressureFixed = fabs(x - box.hi[0]) < eps;
            return node;
        };
        std::vector<mars::fem::Opening<RealType>> openings{{0, box.lo[0], +1}, {0, box.hi[0], -1}};
        opt.params.planar = true;
        Solver solver(*domain, opt.params, channel, openings);

        // Start from uniform flow; start() puts the walls and the inlet at their values.
        thrust::fill(thrust::device, solver.u.data(), solver.u.data() + solver.u.size(), U);
        solver.start();
        if (rank == 0)
            std::cout << "Poiseuille channel: " << solver.globalDofs() << " nodes on " << numRanks
                      << " ranks, Re = U H / nu = " << U * (box.hi[1] - box.lo[1]) / opt.params.nu << "\n";

        // Diagnostics for large runs: the exchange pattern alone, before any time step.
        if (const char* reps = std::getenv("MARS_EXCHANGE_BENCH")) solver.benchmarkExchanges(std::atoi(reps));

        FrameWriter frames(opt.vtuPrefix);
        poiseuille::Monitor monitor(opt.validation, solver, opt.numSteps, U);
        frames.write(solver, *domain, 0, 0.0);
        monitor.afterStep(solver, *domain, 0, 0.0);

        // 3. Time loop.
        for (int step = 1; step <= opt.numSteps; ++step)
        {
            if (!solver.step())
            {
                if (rank == 0) std::cerr << "Step " << step << ": a linear solve did not converge\n";
                exitCode = 1;
                break;
            }
            if (step == 1 && opt.numSteps > 1) solver.resetTiming();
            const double t = step * opt.params.dt;
            if (step % opt.reportEvery == 0 || step == opt.numSteps)
            {
                double norm       = solver.norm(solver.u);
                double continuity = solver.maxContinuity();
                if (rank == 0)
                    std::cout << "Step " << std::setw(6) << step << "  t=" << std::fixed << std::setprecision(4) << t
                              << std::scientific << std::setprecision(10) << "  |u|_M=" << norm
                              << std::setprecision(3) << "  continuity=" << continuity
                              << "  amg(u,v,p)=" << solver.velocityIterations(0) << "/"
                              << solver.velocityIterations(1) << "/" << solver.pressureIterations() << "\n"
                              << std::defaultfloat;
            }
            if (step % opt.vtuEvery == 0 || step == opt.numSteps) frames.write(solver, *domain, step, t);
            monitor.afterStep(solver, *domain, step, t);
        }

        // 4. Result.
        if (exitCode == 0)
        {
            solver.printTiming(opt.numSteps - 1); // step 1 (first-use allocations, BDF1) is left out
            exitCode = monitor.finish(solver, *domain, opt.vtuPrefix, opt.params.rho, opt.params.nu);
        }
    }
    domain.reset();
    MPI_Finalize();
    return exitCode;
}
