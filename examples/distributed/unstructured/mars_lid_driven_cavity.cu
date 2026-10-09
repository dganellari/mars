// Lid-driven cavity: fluid in a closed box, set in motion by its top wall
// sliding with velocity U in x. The classic closed-flow test of an
// incompressible solver, and the smallest MARS Navier-Stokes example.
//
// main() follows three steps:
//   1. mesh and domain   read a hex mesh of the box, or generate a cube on every
//                        rank; cornerstone distributes the elements
//   2. solver            no-slip walls everywhere, the lid moving; no inlet or
//                        outlet, so the pressure is only defined up to a
//                        constant, which the solver removes
//   3. time loop         kinetic energy and continuity of the flow spinning up
//
// It runs on the same NavierStokes solver as the Poiseuille and Taylor-Green
// examples (fem/mars_navier_stokes.hpp).
//
// Example (Re = U L / nu = 100 in the unit cube):
//   mpirun -np 4 ./mars_lid_driven_cavity --cells=32 --nu=0.01 --dt=0.005 --num-steps=400

#include "backend/distributed/unstructured/fem/mars_navier_stokes.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"

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
using Domain   = Solver::Domain;

struct Options
{
    std::string mesh;
    size_t cells    = 0; // generated unit cube with cells^3 hexes
    RealType lidU   = 1;
    Solver::Params params;
    int numSteps    = 100;
    int reportEvery = 10;
    std::string vtuPrefix;
    int vtuEvery   = 50;
    int bucketSize = 64;
};

void printUsage()
{
    std::cout << "Usage: mars_lid_driven_cavity --mesh=FILE | --cells=N [options]\n"
                 "  --mesh=FILE           hex mesh of a box\n"
                 "  --cells=N             generate the unit cube with N^3 hexes instead\n"
                 "  --lid-u=X             lid velocity in x (default 1)\n"
                 "  --rho=X --nu=X        density, kinematic viscosity (default 1, 0.01)\n"
                 "  --dt=X --num-steps=N  time step and step count (default 0.01, 100)\n"
                 "  --tol=X --max-iter=N  relative AMG-PCG tolerance and iteration cap (default 1e-10, 1000)\n"
                 "  --report-every=N      progress line interval (default 10)\n"
                 "  --vtu-output=PREFIX   write PVTU frames; --vtu-every=N (default 50)\n";
}

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
        Solver::Params& p = o.params;
        bool known = take("mesh", o.mesh) || take("cells", o.cells) || take("lid-u", o.lidU) ||
                     take("rho", p.rho) || take("nu", p.nu) || take("dt", p.dt) || take("num-steps", o.numSteps) ||
                     take("tol", p.tolerance) || take("max-iter", p.maxIter) || take("report-every", o.reportEvery) ||
                     take("vtu-output", o.vtuPrefix) || take("vtu-every", o.vtuEvery) ||
                     take("bucket-size", o.bucketSize);
        if (!known)
        {
            if (rank == 0) std::cerr << "Unknown option: " << arg << "\n";
            bad = true;
        }
    }
    const Solver::Params& p = o.params;
    bad = bad || o.mesh.empty() == (o.cells == 0) || !(p.rho > 0) || !(p.nu > 0) || !(p.dt > 0) ||
          !(p.tolerance > 0) || p.maxIter <= 0 || o.numSteps <= 0 || o.reportEvery <= 0 || o.vtuEvery <= 0;
    if (bad)
    {
        if (rank == 0) printUsage();
        exitCode = 1;
        return false;
    }
    return true;
}

// storeOriginalCoords keeps the exact node coordinates; otherwise they are
// decoded from the space-filling-curve keys.
std::unique_ptr<Domain> makeDomain(const Options& o, int rank, int numRanks)
{
    if (!o.mesh.empty()) return std::make_unique<Domain>(o.mesh, rank, numRanks, true, o.bucketSize);
    [[maybe_unused]] auto [nodes, elements, x, y, z, conn] = mars::generateBoxElementPartition<RealType, KeyType>(
        {o.cells, o.cells, o.cells}, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 0.0, rank, numRanks);
    Domain::HostCoordsTuple coords{std::move(x), std::move(y), std::move(z)};
    Domain::HostConnectivityTuple connectivity{std::move(conn[0]), std::move(conn[1]), std::move(conn[2]),
                                               std::move(conn[3]), std::move(conn[4]), std::move(conn[5]),
                                               std::move(conn[6]), std::move(conn[7])};
    return std::make_unique<Domain>(coords, connectivity, rank, numRanks, o.bucketSize, true);
}

// Steps 1 to 3. The domain and the solver hold MPI and Hypre resources, so they
// must be destroyed before MPI_Finalize: they live here.
int runCavity(Options opt, int rank, int numRanks)
{
    // 1. Mesh and domain.
    auto domain = makeDomain(opt, rank, numRanks);

    // 2. Boundary conditions and solver. Every boundary node is a wall; the lid
    //    at z = zmax slides in x. The side walls win on the lid edges.
    const auto box     = mars::fem::boundingBox(*domain);
    const RealType eps = RealType(1e-8) * box.extent();
    const RealType U   = opt.lidU;
    auto cavity = [box, eps, U] __device__(RealType x, RealType y, RealType z) -> mars::fem::NodeCondition<RealType> {
        mars::fem::NodeCondition<RealType> node;
        bool side = fabs(x - box.lo[0]) < eps || fabs(x - box.hi[0]) < eps || fabs(y - box.lo[1]) < eps ||
                    fabs(y - box.hi[1]) < eps;
        bool bottom        = fabs(z - box.lo[2]) < eps;
        bool lid           = fabs(z - box.hi[2]) < eps && !side;
        node.velocityFixed = side || bottom || lid;
        node.velocity[0]   = lid ? U : RealType(0);
        return node;
    };
    Solver solver(*domain, opt.params, cavity, {});
    solver.start(); // fluid at rest, lid moving
    if (rank == 0)
        std::cout << "Lid-driven cavity: " << solver.globalDofs() << " nodes on " << numRanks
                  << " ranks, Re = U L / nu = " << U * (box.hi[0] - box.lo[0]) / opt.params.nu << "\n";

    std::unique_ptr<mars::fem::VTUParallelWriter<KeyType, RealType>> frames;
    if (!opt.vtuPrefix.empty())
        frames = std::make_unique<mars::fem::VTUParallelWriter<KeyType, RealType>>(opt.vtuPrefix);
    auto write = [&](int step, double t) {
        if (!frames) return;
        using FD = mars::fem::VTUParallelWriter<KeyType, RealType>::FieldDesc;
        std::vector<FD> fields{{"velocity", FD::Kind::PointVector3, &solver.u, &solver.v, &solver.w},
                               {"p", FD::Kind::PointScalar, &solver.p, nullptr, nullptr}};
        frames->writeMultiFieldFrame(step, t, *domain, fields);
    };
    write(0, 0.0);

    // 3. Time loop.
    for (int step = 1; step <= opt.numSteps; ++step)
    {
        if (!solver.step())
        {
            if (rank == 0) std::cerr << "Step " << step << ": a linear solve did not converge\n";
            return 1;
        }
        if (step == 1 && opt.numSteps > 1) solver.resetTiming();
        const double t = step * opt.params.dt;
        if (step % opt.reportEvery == 0 || step == opt.numSteps)
        {
            double ke         = solver.kineticEnergy();
            double continuity = solver.maxContinuity();
            if (rank == 0)
                std::cout << "Step " << std::setw(6) << step << "  t=" << std::fixed << std::setprecision(4) << t
                          << std::scientific << std::setprecision(10) << "  KE=" << ke << std::setprecision(3)
                          << "  continuity=" << continuity << "  amg(u,v,w,p)=" << solver.velocityIterations(0)
                          << "/" << solver.velocityIterations(1) << "/" << solver.velocityIterations(2) << "/"
                          << solver.pressureIterations() << "\n"
                          << std::defaultfloat;
        }
        if (step % opt.vtuEvery == 0 || step == opt.numSteps) write(step, t);
    }

    solver.printTiming(opt.numSteps - 1); // step 1 (first-use allocations, BDF1) is left out
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
    if (parseOptions(argc, argv, rank, opt, exitCode)) exitCode = runCavity(opt, rank, numRanks);
    MPI_Finalize();
    return exitCode;
}
