// Poiseuille flow validation driver.
//
// Channel flow: uniform inflow u=Uinf on x=xmin, no-slip on the lateral walls,
// free-velocity outflow on x=xmax with pressure pinned to 0 there (the FLUYA
// reference config from report.html). The long channel lets the flow develop so
// the downstream cross-section relaxes to the analytic parabolic profile:
//   u(y) = U_max * (1 - (2y/H)^2),  U_max = 1.5 * Uinf
// (Hagen-Poiseuille, plane channel).
//
// Uses the CHANNEL FORK of the NS solver (mars_ns_channel_solver.hpp), whose
// one delta is the balanced opening-flux source ported from the pump fork: the
// interior SCS divergence never sees the boundary opening faces, so without the
// source the prescribed inlet is invisible to the pressure solve and the
// channel stays dead (verified: |u| frozen at the inlet-sliver norm, div_max
// pinned, no downstream propagation over t=7).

#include "backend/distributed/unstructured/fem/mars_ns_channel_solver.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"
#include "backend/distributed/unstructured/utils/mars_generate_cube.hpp"

#include <cmath>
#include <vector>
#include <algorithm>
#include <charconv>
#include <fstream>

#include <thrust/copy.h>
#include <thrust/pair.h>

// Per-node share of the element faces on the plane x = plane_x, times `scale`.
// A face counts only from the owned element whose other four nodes lie on
// `side` of the plane, so a face shared by two elements (or two ranks) counts
// once. The reverse halo completes the owned nodes; the forward halo copies
// them to the ghosts. A node's share is the y-z area of its sub-quad (node,
// edge midpoints, face centre), the CVFEM boundary sub-face; for the
// rectangular faces of this mesh it is a quarter of the face.
template<class KeyType, class RealType, class Domain>
void plane_face_areas(const Domain& domain, double plane_x, double tolerance, int side,
                      double scale, cstone::DeviceVector<RealType>& area)
{
    const size_t n = domain.getNodeCount();
    area.resize(n);
    thrust::fill(thrust::device_pointer_cast(area.data()),
                 thrust::device_pointer_cast(area.data() + n), RealType(0));
    const auto cp = connPtrs<HexTag, KeyType>(domain.getElementToNodeConnectivity());
    const size_t first = domain.startIndex();
    const size_t count = domain.localElementCount();
    thrust::for_each(thrust::device,
        thrust::counting_iterator<size_t>(first), thrust::counting_iterator<size_t>(first + count),
        [c0=cp[0], c1=cp[1], c2=cp[2], c3=cp[3], c4=cp[4], c5=cp[5], c6=cp[6], c7=cp[7],
         x=domain.getNodeX().data(), y=domain.getNodeY().data(), z=domain.getNodeZ().data(),
         a=area.data(), plane_x, tolerance, side, scale] __device__ (size_t e) {
            const KeyType nodes[8] = {c0[e],c1[e],c2[e],c3[e],c4[e],c5[e],c6[e],c7[e]};
            KeyType face[4];
            int on_plane = 0;
            for (int j = 0; j < 8; ++j)
            {
                const double d = double(x[nodes[j]]) - plane_x;
                if (fabs(d) <= tolerance) { if (on_plane < 4) face[on_plane] = nodes[j]; ++on_plane; }
                else if (d * side <= 0) return;
            }
            if (on_plane != 4) return;
            double cy = 0, cz = 0;
            for (int k = 0; k < 4; ++k) { cy += 0.25 * y[face[k]]; cz += 0.25 * z[face[k]]; }
            // Local node order does not give the face cycle; sort by angle.
            double angle[4];
            for (int k = 0; k < 4; ++k) angle[k] = atan2(z[face[k]] - cz, y[face[k]] - cy);
            for (int k = 1; k < 4; ++k)
                for (int j = k; j > 0 && angle[j] < angle[j - 1]; --j)
                {
                    double t = angle[j]; angle[j] = angle[j - 1]; angle[j - 1] = t;
                    KeyType f = face[j]; face[j] = face[j - 1]; face[j - 1] = f;
                }
            for (int k = 0; k < 4; ++k)
            {
                const KeyType p = face[k], next = face[(k + 1) % 4], prev = face[(k + 3) % 4];
                const double ny = 0.5 * (y[prev] - y[next]), nz = 0.5 * (z[prev] - z[next]);
                const double share = 0.5 * fabs((cy - y[p]) * nz - (cz - z[p]) * ny);
                atomicAdd(&a[p], RealType(scale * share));
            }
        });
    if (cudaDeviceSynchronize() != cudaSuccess) MPI_Abort(MPI_COMM_WORLD, 1);
    domain.reverseExchangeNodeHaloAdd(area);
    domain.exchangeNodeHalo(area);
}

// x of the node plane nearest `target`, the same on every rank: global minimum
// of (|x - target|, x), so a tie resolves to the smaller x.
template<class RealType>
double nearest_node_plane(const cstone::DeviceVector<RealType>& d_x, size_t n, double target)
{
    using Candidate = thrust::pair<double, double>;
    const double inf = std::numeric_limits<double>::infinity();
    const Candidate local = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
        [x=d_x.data(), target] __device__ (size_t i) -> Candidate {
            return Candidate(fabs(double(x[i]) - target), double(x[i]));
        }, Candidate(inf, inf),
        [] __device__ (const Candidate& a, const Candidate& b) -> Candidate { return b < a ? b : a; });
    double distance = inf, plane = inf;
    MPI_Allreduce(&local.first, &distance, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    const double candidate = local.first == distance ? local.second : inf;
    MPI_Allreduce(&candidate, &plane, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return plane;
}

template<class Function>
struct OwnedTerm
{
    const uint8_t* ownership;
    Function f;
    __device__ double operator()(size_t i) const { return ownership[i] == 1 ? f(i) : 0.0; }
};

// Global sum over owned nodes of f(i).
template<class Function>
double owned_sum(const uint8_t* ownership, size_t n, Function f)
{
    const double local = thrust::transform_reduce(thrust::device,
        thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
        OwnedTerm<Function>{ownership, f}, 0.0, thrust::plus<double>());
    double global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

// Analytic plane-Poiseuille profile, normalized cross-section coordinate
// eta in [-1, 1] (eta=0 center, |eta|=1 walls). U_max = 1.5 * U_avg.
struct PoiseuilleProfile
{
    double uMax;
    double wallLo;   // cross-section coordinate at one wall
    double wallHi;   // cross-section coordinate at the other wall

    __host__ __device__ double eta(double s) const
    {
        double mid  = 0.5 * (wallLo + wallHi);
        double half = 0.5 * (wallHi - wallLo);
        return half > 0 ? (s - mid) / half : 0.0;
    }
    __host__ __device__ double analytic(double s) const
    {
        double e = eta(s);
        return uMax * (1.0 - e * e);
    }
};

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    using KeyType  = uint64_t;
    using RealType = double;

    std::string meshFile;
    // --cells=NX,NY generates the public channel's box, one cell through z, on
    // every rank without a mesh file (scaling runs).
    std::array<size_t, 2> cells{0, 0};
    double yGrading = 0.0;
    CvfemKernelVariant kernelVariant = CvfemKernelVariant::Tensor;
    int    blockSize  = 256;
    int    bucketSize = 64;
    int    maxIter    = 1000;
    double tolerance  = 1e-10;
    double rho        = 1.0;
    double nu         = 0.01;
    double dt         = 0.01;
    int    numSteps   = 1000;
    int    vtuEvery   = 50;
    std::string vtuPrefix;
    std::string comparisonPrefix;
    std::string comparisonEveryText;
    int    comparisonEvery = 50;
    bool   comparisonRequested = false;
    bool   comparisonEveryRequested = false;
    SolverKind        solverKind    = SolverKind::CG;
    // The planar pressure action includes the velocity constraint Q:
    // D Q M^-1 D^T, with symmetric pressure elimination.
    PressureSolveKind pressureSolve = PressureSolveKind::DDT;
    double Uinf       = 1.0;

    // Drive mode. Default is INLET -- the FLUYA reference config (report.html):
    // velocity-Dirichlet inlet u=Uinf, static-pressure outlet (p=0, velocity
    // FREE), no-slip walls, symmetry on the thin z-faces. The free pressure
    // outlet lets mass balance through the pressure field and the parabola
    // develops downstream. (--drive=bodyforce is the textbook periodic-style
    // alternative but on this free-in/out mesh the projection cancels the mean
    // flow -> collapses; kept for study, needs periodic-x to work.)
    std::string driveMode = "inlet";       // inlet | bodyforce
    double bodyForceX = -1.0;              // <0 = auto (target U_max = Uinf)

    // Outlet cross-section probe: nodes within xTol of profileX. Default -1
    // means "auto" = pick a plane 90% down the channel from xmin to xmax.
    double profileX    = -1.0;
    double profileXTol = -1.0;   // auto = one element width
    // Cross-section axis the parabola varies over: "y" or "z" (the wall-normal
    // direction). The streamwise component is always u (x). Default y.
    std::string crossAxis = "y";
    // Seed interior u=Uinf at the IC so the long channel starts full (default
    // ON). Without it the rest-IC inlet front takes many flow-throughs to fill
    // an 11-unit channel and the profile never develops in a short run.
    bool seedInterior = true;
    // Balanced opening-flux source (the pump-fork fix; channel-fork solver).
    // Default ON in inlet mode -- without it the inlet is invisible to the
    // pressure solve. --no-opening-flux-source gives the A/B baseline.
    bool openingFluxSource = true;
    bool planar_projection = false;
    bool pressureAmg = false;
    bool velocityAmg = false;
    bool forceBdf1 = false;   // --bdf1: 1st-order time stepping (disable BDF2)
    // Regression-check mode: exit 1 unless the final profile RMS is below
    // rmsTol (the FLUYA reference tolerance) AND every interior-flux ratio is
    // within fluxTol of 1. Drives the ctest entry.
    bool   checkMode = false;
    double rmsTol    = 6e-3;
    double fluxTol   = 0.10;
    double steady_tol = 1e-6;
    double continuity_tol = 1e-6;

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if      (arg.find("--mesh=") == 0)         meshFile   = arg.substr(7);
        else if (arg.find("--cells=") == 0)
        {
            const std::string v = arg.substr(8);
            const size_t comma = v.find(',');
            if (comma != std::string::npos)
                cells = {size_t(std::stoull(v.substr(0, comma))), size_t(std::stoull(v.substr(comma + 1)))};
        }
        else if (arg.find("--y-grading=") == 0)    yGrading   = std::stod(arg.substr(12));
        else if (arg.find("--nu=") == 0)           nu         = std::stod(arg.substr(5));
        else if (arg.find("--dt=") == 0)           dt         = std::stod(arg.substr(5));
        else if (arg.find("--rho=") == 0)          rho        = std::stod(arg.substr(6));
        else if (arg.find("--num-steps=") == 0)    numSteps   = std::stoi(arg.substr(12));
        else if (arg.find("--vtu-every=") == 0)    vtuEvery   = std::stoi(arg.substr(12));
        else if (arg.find("--vtu-output=") == 0)   vtuPrefix  = arg.substr(13);
        else if (arg.find("--comparison-output=") == 0)
        {
            comparisonRequested = true;
            comparisonPrefix = arg.substr(20);
        }
        else if (arg.find("--comparison-every=") == 0)
        {
            comparisonEveryRequested = true;
            comparisonEveryText = arg.substr(19);
        }
        else if (arg.find("--lid-u=") == 0)        Uinf       = std::stod(arg.substr(8));
        else if (arg.find("--uinf=") == 0)         Uinf       = std::stod(arg.substr(7));
        else if (arg.find("--profile-x=") == 0)    profileX   = std::stod(arg.substr(12));
        else if (arg.find("--profile-xtol=") == 0) profileXTol= std::stod(arg.substr(15));
        else if (arg.find("--cross-axis=") == 0)   crossAxis  = arg.substr(13);
        else if (arg.find("--drive=") == 0)        driveMode  = arg.substr(8);
        else if (arg.find("--body-force-x=") == 0) bodyForceX = std::stod(arg.substr(15));
        else if (arg == "--no-seed-interior")      seedInterior = false;
        else if (arg == "--no-opening-flux-source") openingFluxSource = false;
        else if (arg == "--planar-ddt")            planar_projection = true;
        else if (arg == "--pressure-amg")          pressureAmg = true;
        else if (arg == "--velocity-amg")          velocityAmg = true;
        else if (arg == "--bdf1")                  forceBdf1 = true;
        else if (arg == "--check")                 checkMode = true;
        else if (arg.find("--rms-tol=") == 0)      rmsTol  = std::stod(arg.substr(10));
        else if (arg.find("--flux-tol=") == 0)     fluxTol = std::stod(arg.substr(11));
        else if (arg.find("--steady-tol=") == 0)   steady_tol = std::stod(arg.substr(13));
        else if (arg.find("--continuity-tol=") == 0) continuity_tol = std::stod(arg.substr(17));
        else if (arg.find("--block-size=") == 0)   blockSize  = std::stoi(arg.substr(13));
        else if (arg.find("--bucket-size=") == 0)  bucketSize = std::stoi(arg.substr(14));
        else if (arg.find("--max-iter=") == 0)     maxIter    = std::stoi(arg.substr(11));
        else if (arg.find("--tol=") == 0)          tolerance  = std::stod(arg.substr(6));
        else if (arg == "--kernel=wmma_tensor")    kernelVariant = CvfemKernelVariant::WmmaTensor;
        else if (arg.find("--solver=") == 0)
        {
            std::string v = arg.substr(9);
            if      (v == "cg")    solverKind = SolverKind::CG;
            else if (v == "hypre") solverKind = SolverKind::Hypre;
        }
        else if (arg.find("--pressure-solve=") == 0)
        {
            std::string v = arg.substr(17);
            if      (v == "K")   pressureSolve = PressureSolveKind::K;
            else if (v == "DDT") pressureSolve = PressureSolveKind::DDT;
        }
        else if (arg[0] != '-' && meshFile.empty()) meshFile = arg;
    }

    if (!(std::isfinite(dt) && dt > 0 && std::isfinite(nu) && nu > 0
          && std::isfinite(rho) && rho > 0 && std::isfinite(Uinf) && Uinf > 0
          && std::isfinite(tolerance) && tolerance > 0 && maxIter > 0 && numSteps > 0
          && std::isfinite(rmsTol) && rmsTol > 0 && std::isfinite(fluxTol) && fluxTol > 0
          && std::isfinite(steady_tol) && steady_tol > 0
          && std::isfinite(continuity_tol) && continuity_tol > 0))
    {
        if (rank == 0) std::cerr << "ERROR: require positive finite physical parameters and tolerances\n";
        MPI_Finalize();
        return 1;
    }

    int comparisonError = 0;
    if (comparisonRequested && comparisonPrefix.empty()) comparisonError |= 1;
    if (comparisonEveryRequested)
    {
        int parsed = 0;
        auto result = std::from_chars(comparisonEveryText.data(),
                                      comparisonEveryText.data() + comparisonEveryText.size(), parsed);
        if (result.ec != std::errc{} ||
            result.ptr != comparisonEveryText.data() + comparisonEveryText.size() || parsed <= 0)
            comparisonError |= 2;
        else
            comparisonEvery = parsed;
    }
    if (comparisonEveryRequested && !comparisonRequested) comparisonError |= 4;
    int comparisonErrorGlobal = 0;
    MPI_Allreduce(&comparisonError, &comparisonErrorGlobal, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD);
    if (comparisonErrorGlobal)
    {
        if (rank == 0)
        {
            std::cerr << "Invalid comparison export options:";
            if (comparisonErrorGlobal & 1) std::cerr << " --comparison-output needs a nonempty prefix;";
            if (comparisonErrorGlobal & 2) std::cerr << " --comparison-every needs a positive integer;";
            if (comparisonErrorGlobal & 4) std::cerr << " --comparison-every requires --comparison-output;";
            std::cerr << "\n";
        }
        MPI_Finalize();
        return 1;
    }

    const bool generated = cells[0] > 0 && cells[1] > 0;
    if (meshFile.empty() == !generated || !(yGrading >= 0 && yGrading < 1))
    {
        if (rank == 0)
            std::cout << "Usage: " << argv[0] << " --mesh=FILE [options]\n\n"
                      << "Poiseuille channel flow validation (parabolic profile).\n\n"
                      << "  --mesh=FILE         Mesh file (.exo / .mesh); or\n"
                      << "  --cells=NX,NY       Generate the channel [0,10]x[0,1]x[0,0.06], one cell in z\n"
                      << "  --y-grading=S       Cluster generated y spacing at the walls, S in [0,1) (default 0)\n"
                      << "  --drive=bodyforce|inlet  Flow driver (default inlet)\n"
                      << "                      bodyforce: constant streamwise force G, no-slip y-walls,\n"
                      << "                                 U_max=G*H^2/(8 rho nu). inlet: prescribed velocity\n"
                      << "                                 inlet+free-velocity pressure outlet\n"
                      << "  --body-force-x=G    Body force value (default: auto = G to hit U_max=Uinf)\n"
                      << "  --uinf=VALUE        Target U_max (bodyforce) / inflow speed (inlet); default 1.0\n"
                      << "  --nu=VALUE          Kinematic viscosity (default 0.01)\n"
                      << "  --dt=VALUE          Timestep (default 0.01)\n"
                      << "  --rho=VALUE         Density (default 1.0)\n"
                      << "  --num-steps=N       Time steps (default 1000)\n"
                      << "  --cross-axis=y|z    Wall-normal axis the parabola varies over (default y)\n"
                      << "  --no-seed-interior  Start from rest (default seeds interior with mean flow)\n"
                      << "  --no-opening-flux-source  Disable the balanced opening-flux source (A/B baseline;\n"
                      << "                      without it the inlet is invisible to the pressure solve)\n"
                      << "  --check             Regression mode: exit 1 unless RMS < rms-tol and flux\n"
                      << "                      ratios within flux-tol of 1; planar mode also requires continuity and steadiness\n"
                      << "  --rms-tol=VAL       RMS pass threshold (default 6e-3, the reference tol)\n"
                      << "  --flux-tol=VAL      Relative flux-ratio tolerance (default 0.10)\n"
                      << "  --profile-x=X       Outlet probe plane x (default: 90% down the channel)\n"
                      << "  --profile-xtol=TOL  Half-width of the probe plane (default: one element)\n"
                      << "  --solver=cg|hypre   Linear solver (default cg)\n"
                      << "  --planar-ddt        Consistent xy projection for a one-layer rectilinear channel (inlet + CG/DDT)\n"
                      << "  --pressure-amg      Assemble the DDT pressure operator and solve it with Hypre PCG + BoomerAMG\n"
                      << "  --velocity-amg      Solve the implicit velocity systems with Hypre PCG + BoomerAMG\n"
                      << "  --steady-tol=TOL    Final 20-step velocity change / Uinf (default 1e-6)\n"
                      << "  --continuity-tol=TOL  Final planar continuity RMS*H/U and relative boundary balance (default 1e-6)\n"
                      << "  --pressure-solve=K|DDT  Pressure operator (default DDT; K FAILs on this channel)\n"
                      << "  --vtu-output=PREFIX Write VTU/PVTU/PVD frames\n"
                      << "  --vtu-every=N       Frame every N steps (default 50)\n"
                      << "  --comparison-output=PREFIX  Export public Poiseuille u,v,w,p frames\n"
                      << "  --comparison-every=N       Comparison frame every N steps (default 50)\n"
                      << "  --tol=VALUE         Solver tolerance (default 1e-10)\n"
                      << "  --max-iter=N        Solver max iters (default 1000)\n";
        MPI_Finalize();
        return 1;
    }

    if (rank == 0)
    {
        const double Re = Uinf * 1.0 / nu;
        std::cout << "\n========================================\n";
        std::cout << "MARS Poiseuille channel-flow validation\n";
        std::cout << "========================================\n";
        std::cout << "Drive:     " << driveMode
                  << (driveMode == "bodyforce"
                          ? " (constant streamwise G, no-slip y-walls)"
                          : " (prescribed inlet u=Uinf, mass-conserving outlet)") << "\n";
        std::cout << "rho       = " << rho << "\n";
        std::cout << "nu        = " << nu  << "\n";
        std::cout << "dt        = " << dt  << "\n";
        std::cout << "Re (~L*U/nu, L=1) = " << Re << "\n";
        std::cout << "numSteps  = " << numSteps << " (T_final = " << numSteps * dt << ")\n";
        if (generated)
            std::cout << "Mesh:    generated " << cells[0] << " x " << cells[1] << " x 1 channel, y grading "
                      << yGrading << "\n";
        else
            std::cout << "Mesh:    " << meshFile << "\n";
        std::cout << "Cross-axis (wall-normal): " << crossAxis << "\n";
        std::cout << "Pressure solve: " << (pressureSolve == PressureSolveKind::K ? "K" : "DDT") << "\n";
        std::cout << "MPI ranks: " << numRanks << "\n";
        std::cout << "========================================\n\n";
    }

    AmrManager<HexTag, KeyType, RealType>::Config amrConfig;
    amrConfig.maxLevels  = 0;
    amrConfig.blockSize  = blockSize;
    amrConfig.bucketSize = bucketSize;

    AmrManager<HexTag, KeyType, RealType> amr(amrConfig);
    if (generated)
    {
        using Domain = AmrManager<HexTag, KeyType, RealType>::Domain;
        [[maybe_unused]] auto [nodes, elements, x, y, z, conn] = generateBoxElementPartition<RealType, KeyType>(
            {cells[0], cells[1], 1}, {0.0, 0.0, 0.0}, {10.0, 1.0, 0.06}, yGrading, rank, numRanks);
        typename Domain::HostCoordsTuple coords{std::move(x), std::move(y), std::move(z)};
        typename Domain::HostConnectivityTuple connectivity{
            std::move(conn[0]), std::move(conn[1]), std::move(conn[2]), std::move(conn[3]),
            std::move(conn[4]), std::move(conn[5]), std::move(conn[6]), std::move(conn[7])};
        amr.initialize(std::make_unique<Domain>(coords, connectivity, rank, numRanks, bucketSize, true),
                       rank, numRanks);
    }
    else
        amr.initialize(meshFile, rank, numRanks);

    if (rank == 0)
    {
        std::cout << "Mesh: " << amr.domain().getElementCount() << " elements, "
                  << amr.domain().getNodeCount() << " nodes\n\n";
    }

    NSStepper<KeyType, RealType> s{amr.domain(), solverKind, blockSize, maxIter,
                                   RealType(tolerance), rank, numRanks};
    s.lidU   = RealType(Uinf);
    s.Uinf   = RealType(Uinf);
    s.useLegacyGradient = true;
    s.pressureSolve     = pressureSolve;
    if (forceBdf1) s.useBdf2 = false;   // 1st-order time stepping

    const bool useBodyForce = (driveMode == "bodyforce");
    double targetUMax = Uinf;   // inlet: 1.5*Uinf below; bodyforce: from G

    // Node coordinates + global bbox, needed BEFORE setup: inlet mode feeds the
    // BCs through the solver's Pump node-list path so every derived BC object
    // (halo-synced per-node mask, matrix row/col elimination, pressure pin) is
    // built inside setupNSStepper from the REAL masks. The first multirank
    // attempt rebuilt d_isBdryDof after setup: the derived state kept the
    // original thin-z all-Dirichlet marking, owner and ghost ranks disagreed,
    // and the run blew up at step 1 (single rank has no ghosts -> worked).
    amr.domain().getLocalToGlobalSfcMap();
    const size_t nNodes = amr.domain().getNodeCount();
    std::vector<RealType> hx(nNodes), hy(nNodes), hz(nNodes);
    {
        const auto& d_x = amr.domain().getNodeX();
        const auto& d_y = amr.domain().getNodeY();
        const auto& d_z = amr.domain().getNodeZ();
        cudaMemcpy(hx.data(), d_x.data(), nNodes * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hy.data(), d_y.data(), nNodes * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hz.data(), d_z.data(), nNodes * sizeof(RealType), cudaMemcpyDeviceToHost);
    }
    RealType xMinL = *std::min_element(hx.begin(), hx.end());
    RealType xMaxL = *std::max_element(hx.begin(), hx.end());
    RealType yMinL = *std::min_element(hy.begin(), hy.end());
    RealType yMaxL = *std::max_element(hy.begin(), hy.end());
    RealType zMinL = *std::min_element(hz.begin(), hz.end());
    RealType zMaxL = *std::max_element(hz.begin(), hz.end());
    RealType xMin = xMinL, xMax = xMaxL, yMin = yMinL, yMax = yMaxL;
    RealType zMin = zMinL, zMax = zMaxL;
    MPI_Allreduce(&xMinL, &xMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&xMaxL, &xMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&yMinL, &yMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&yMaxL, &yMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&zMinL, &zMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&zMaxL, &zMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    // 1e-5, NOT 1e-4: the first interior node sits 1e-4 off the wall; eps=1e-4
    // float-rounded the TOP first interior node into the wall set.
    RealType eps = 1e-5 * std::max(RealType(1), yMax - yMin);
    double H = static_cast<double>(yMax - yMin);   // channel height
    if (planar_projection)
    {
        int invalid = useBodyForce || crossAxis != "y" || !openingFluxSource
                      || solverKind != SolverKind::CG || pressureSolve != PressureSolveKind::DDT;
        invalid |= !(xMax > xMin && yMax > yMin && zMax > zMin);
        for (size_t i = 0; i < nNodes; ++i)
            invalid |= std::abs(hz[i] - zMin) > eps && std::abs(hz[i] - zMax) > eps;
        const auto& d_x = amr.domain().getNodeX();
        const auto& d_y = amr.domain().getNodeY();
        const auto& d_z = amr.domain().getNodeZ();
        const auto cp = connPtrs<HexTag, KeyType>(amr.domain().getElementToNodeConnectivity());
        const size_t first = amr.domain().startIndex();
        const size_t count = amr.domain().localElementCount();
        const RealType cell_tolerance = RealType(1e-10) * std::max({xMax-xMin, yMax-yMin, zMax-zMin});
        int invalid_cells = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<size_t>(first), thrust::counting_iterator<size_t>(first + count),
            [c0=cp[0], c1=cp[1], c2=cp[2], c3=cp[3], c4=cp[4], c5=cp[5], c6=cp[6], c7=cp[7],
             x=d_x.data(), y=d_y.data(), z=d_z.data(), cell_tolerance] __device__ (size_t e) -> int {
                const KeyType n[8] = {c0[e],c1[e],c2[e],c3[e],c4[e],c5[e],c6[e],c7[e]};
                RealType coords[8][3];
                for (int j=0;j<8;++j) {coords[j][0]=x[n[j]];coords[j][1]=y[n[j]];coords[j][2]=z[n[j]];}
                return channel_rectilinear_hex(coords, cell_tolerance) ? 0 : 1;
            }, 0, thrust::maximum<int>());
        invalid |= invalid_cells;
        int any_invalid = 0;
        MPI_Allreduce(&invalid, &any_invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (any_invalid)
        {
            if (rank == 0) std::cerr << "--planar-ddt requires inlet, y walls, CG/DDT, openings and rectilinear cells in one z layer\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
        s.planar_projection = true;
        s.useLegacyGradient = false;
    }
    s.pressure_amg = pressureAmg;
    s.velocity_amg = velocityAmg;

    if (useBodyForce)
    {
        // Study mode (single-rank only): shared Channel marking at setup, then
        // a post-setup wall-only rebuild below -- known NOT multirank-safe.
        s.bcKind = NSStepper<KeyType, RealType>::BCKind::Channel;
    }
    else
    {
        // INLET mode (FLUYA reference): walls=no-slip, inlet u=Uinf along +x,
        // outlet velocity stays natural-Neumann (outletU<0 default) and the
        // whole outlet plane gets the p=0 pressure Dirichlet inside setup --
        // exactly the reference outlet. z faces + interior free. Corner policy:
        // inlet wins (corner nodes go in inletNodes only).
        s.bcKind = NSStepper<KeyType, RealType>::BCKind::Pump;
        s.inletDirX = 1; s.inletDirY = 0; s.inletDirZ = 0;
        for (size_t i = 0; i < nNodes; ++i)
        {
            bool onInflow  = std::abs(hx[i] - xMin) < eps;
            bool onYWall   = std::abs(hy[i] - yMin) < eps || std::abs(hy[i] - yMax) < eps;
            bool onOutflow = std::abs(hx[i] - xMax) < eps;
            if (planar_projection && onYWall) s.wallNodes.push_back(int(i));
            else if (onInflow) s.inletNodes.push_back(int(i));
            else if (onYWall) s.wallNodes.push_back(int(i));
            if (planar_projection && onInflow && onYWall)
                s.pressure_null_nodes.push_back(int(i));
            if (onOutflow)    s.outletNodes.push_back(int(i));
        }
        targetUMax = 1.5 * Uinf;
        if (rank == 0)
            std::cout << "Drive: INLET (FLUYA-ref via pump BC path): walls=no-slip, "
                      << "inlet u=" << Uinf << ", outlet natural + p=0"
                      << (planar_projection ? ", planar symmetry w=0; walls win at corners\n" : ", z free\n");
    }

    setupNSStepper<KeyType, RealType>(s, RealType(nu), RealType(dt), kernelVariant);

    if (useBodyForce)
    {
        const auto& d_ownership = s.ownershipMap();
        size_t n = s.nodeCount;
        std::vector<int>     h_n2d(n);
        std::vector<uint8_t> h_own(n);
        cudaMemcpy(h_n2d.data(), s.d_node_to_dof.data(), n * sizeof(int),     cudaMemcpyDeviceToHost);
        cudaMemcpy(h_own.data(), d_ownership.data(),     n * sizeof(uint8_t), cudaMemcpyDeviceToHost);

        std::vector<uint8_t>  hostIsBdry(s.numOwnedDofs, 0);
        std::vector<RealType> hostU(n, 0), hostV(n, 0), hostW(n, 0);
        long nWall = 0;
        for (size_t i = 0; i < n; ++i)
        {
            bool onYWall = std::abs(hy[i] - yMin) < eps || std::abs(hy[i] - yMax) < eps;
            if (!onYWall) continue;                 // x-planes + z-faces + interior free
            hostU[i] = 0; hostV[i] = 0; hostW[i] = 0;
            if (h_own[i] != 1) continue;
            int dof = h_n2d[i];
            if (dof < 0 || dof >= s.numOwnedDofs) continue;
            hostIsBdry[dof] = 1;
            ++nWall;
        }
        double G = (bodyForceX >= 0) ? bodyForceX : 8.0 * rho * nu * targetUMax / (H * H);
        s.bodyForceX = RealType(G);
        targetUMax = G * H * H / (8.0 * rho * nu);
        long nWallG = nWall;
        MPI_Allreduce(&nWall, &nWallG, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (rank == 0)
            std::cout << "Drive: BODY FORCE G=" << G << " (x), H=" << H
                      << " -> analytic U_max=" << targetUMax << "; y-walls=" << nWallG << "\n";
        thrust::copy(hostIsBdry.begin(), hostIsBdry.end(),
                     thrust::device_pointer_cast(s.d_isBdryDof.data()));
        thrust::copy(hostU.begin(), hostU.end(), thrust::device_pointer_cast(s.d_uTarget.data()));
        thrust::copy(hostV.begin(), hostV.end(), thrust::device_pointer_cast(s.d_vTarget.data()));
        thrust::copy(hostW.begin(), hostW.end(), thrust::device_pointer_cast(s.d_wTarget.data()));
        cudaDeviceSynchronize();
    }
    else if (openingFluxSource)
    {
        // Balanced opening-flux source: per-node OUTWARD x-areas of the inlet
        // and outlet faces; the outlet is rescaled per step (oScale) so the net
        // source is machine-zero. The openings are x planes, so y/z stay zero.
        const size_t n = s.nodeCount;
        plane_face_areas<KeyType>(amr.domain(), double(xMin), double(eps), +1, -1.0, s.d_openInAreaX);
        plane_face_areas<KeyType>(amr.domain(), double(xMax), double(eps), -1, +1.0, s.d_openOutAreaX);
        for (auto* d : {&s.d_openInAreaY, &s.d_openInAreaZ, &s.d_openOutAreaY, &s.d_openOutAreaZ})
        {
            d->resize(n);
            thrust::fill(thrust::device_pointer_cast(d->data()),
                         thrust::device_pointer_cast(d->data() + n), RealType(0));
        }
        const uint8_t* own = s.ownershipMap().data();
        const double Ain  = owned_sum(own, n, [a=s.d_openInAreaX.data()] __device__ (size_t i) -> double { return -double(a[i]); });
        const double Aout = owned_sum(own, n, [a=s.d_openOutAreaX.data()] __device__ (size_t i) -> double { return double(a[i]); });

        RealType srcOutletU = (Aout > 0) ? RealType(Uinf * Ain / Aout) : RealType(Uinf);
        s.openInletVel[0]  = RealType(Uinf);
        s.openOutletVel[0] = srcOutletU;
        s.useOpeningFluxSource = true;
        if (rank == 0 && planar_projection)
            std::cout << "Opening flux: nodal inlet targets and solved outlet velocity; no outlet rescaling"
                      << std::setprecision(17) << " Ain=" << Ain << " Aout=" << Aout
                      << std::setprecision(6) << '\n';
        else if (rank == 0)
            std::cout << "Opening-flux source: Ain=" << Ain << " Aout=" << Aout
                      << " Qin=" << Uinf * Ain << " outletU=" << srcOutletU
                      << " (net zeroed per step via oScale)\n";
    }

    applyInitialCondition<KeyType, RealType>(s);

    // Seed interior streamwise velocity (default ON) so the channel starts with
    // flow and relaxes to the parabola quickly instead of growing from rest.
    // Body-force mode seeds the MEAN velocity (2/3 U_max); inlet mode seeds Uinf.
    // With the BC rebuilt above, interior nodes are genuinely non-boundary, so
    // this takes effect.
    // Pump-path IC already seeds interior u=Uinf in inlet mode; the explicit
    // seed is only for the body-force study path.
    RealType uSeed = RealType(2.0 / 3.0 * targetUMax);
    if (seedInterior && useBodyForce)
    {
        const auto& d_ownership = s.ownershipMap();
        size_t n = s.nodeCount;
        std::vector<uint8_t> h_bdry(s.numOwnedDofs);
        std::vector<int>     h_n2d(n);
        std::vector<uint8_t> h_own(n);
        std::vector<RealType> hu(n);
        cudaMemcpy(h_bdry.data(), s.d_isBdryDof.data(),  s.numOwnedDofs * sizeof(uint8_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_n2d.data(),  s.d_node_to_dof.data(), n * sizeof(int),     cudaMemcpyDeviceToHost);
        cudaMemcpy(h_own.data(),  d_ownership.data(),     n * sizeof(uint8_t), cudaMemcpyDeviceToHost);
        cudaMemcpy(hu.data(),     s.d_u.data(),           n * sizeof(RealType), cudaMemcpyDeviceToHost);
        for (size_t i = 0; i < n; ++i)
        {
            if (h_own[i] != 1) continue;
            int dof = h_n2d[i];
            if (dof < 0 || dof >= s.numOwnedDofs) continue;
            if (h_bdry[dof]) continue;        // leave wall/inlet Dirichlet values
            hu[i] = uSeed;
        }
        cudaMemcpy(s.d_u.data(), hu.data(), n * sizeof(RealType), cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();
        s.domain.exchangeNodeHalo(s.d_u);
        if (rank == 0)
            std::cout << "Interior IC: seeded u=" << uSeed << " (channel starts with flow)\n";
    }

    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriterU;
    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriterP;
    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriterUmag;
    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> comparisonWriter;
    if (!vtuPrefix.empty())
    {
        vtuWriterU    = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix + "_u");
        vtuWriterP    = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix + "_p");
        vtuWriterUmag = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix + "_umag");
        if (rank == 0)
            std::cout << "VTU output enabled: " << vtuPrefix << "_{u,p,umag}_step*.pvtu\n\n";
    }
    if (comparisonRequested)
    {
        comparisonWriter = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(
            comparisonPrefix, true);
        if (rank == 0)
            std::cout << "Comparison output enabled: " << comparisonPrefix
                      << "_step*.pvtu (u,v,w,p; every " << comparisonEvery << " steps)\n";
    }

    auto writeVtuFrame = [&] (int step, double t)
    {
        if (!vtuWriterU) return;
        vtuWriterU->writeFrame(step, t, amr.domain(), s.d_u, "u");
        vtuWriterP->writeFrame(step, t, amr.domain(), s.d_p, "p");
        cstone::DeviceVector<RealType> d_umag(s.nodeCount, RealType(0));
        const RealType* uPtr = s.d_u.data();
        const RealType* vPtr = s.d_v.data();
        const RealType* wPtr = s.d_w.data();
        RealType* mPtr       = d_umag.data();
        thrust::for_each(thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(s.nodeCount),
            [uPtr, vPtr, wPtr, mPtr] __device__ (size_t i) {
                RealType u = uPtr[i], v = vPtr[i], w = wPtr[i];
                mPtr[i] = sqrt(u * u + v * v + w * w);
            });
        cudaDeviceSynchronize();
        vtuWriterUmag->writeFrame(step, t, amr.domain(), d_umag, "umag");
    };

    auto writeComparisonFrame = [&] (int step, double t)
    {
        if (!comparisonWriter) return;
        using FD = typename fem::VTUParallelWriter<KeyType, RealType>::FieldDesc;
        std::vector<FD> fields{
            {"u", FD::Kind::PointScalar, &s.d_u, nullptr, nullptr},
            {"v", FD::Kind::PointScalar, &s.d_v, nullptr, nullptr},
            {"w", FD::Kind::PointScalar, &s.d_w, nullptr, nullptr},
            {"p", FD::Kind::PointScalar, &s.d_p, nullptr, nullptr}
        };
        comparisonWriter->writeMultiFieldFrame(step, t, amr.domain(), fields);
    };

    if (vtuWriterU) writeVtuFrame(0, 0.0);
    if (comparisonWriter) writeComparisonFrame(0, 0.0);

    // Step-0 norm: confirms the seeded IC actually populated d_u before any step.
    {
        RealType nu0 = computeWeightedL2Norm<KeyType, RealType>(s, s.d_u);
        if (rank == 0)
            std::cout << "Step    0: t=0.0000 |u|=" << std::scientific
                      << std::setprecision(3) << nu0 << " (IC)\n" << std::defaultfloat;
    }

    // A small timestep can hide drift in a one-step change. Compare a fixed window.
    const int steady_window = std::min(20, numSteps - 1);
    std::array<cstone::DeviceVector<RealType>, 3> d_steady_start;
    if (checkMode && planar_projection)
        for (auto& component : d_steady_start) component.resize(s.nodeCount);
    for (int step = 1; step <= numSteps; ++step)
    {
        if (checkMode && planar_projection && step == numSteps - steady_window + 1)
        {
            thrust::copy(thrust::device, s.d_u.begin(), s.d_u.end(), d_steady_start[0].begin());
            thrust::copy(thrust::device, s.d_v.begin(), s.d_v.end(), d_steady_start[1].begin());
            thrust::copy(thrust::device, s.d_w.begin(), s.d_w.end(), d_steady_start[2].begin());
        }
        s.check_projection = checkMode && step == numSteps;
        runNsStep<KeyType, RealType>(s, RealType(dt), RealType(nu), RealType(rho));
        // Step 1 carries first-use allocations and the BDF1 start; time steps 2..N.
        if (step == 1 && numSteps > 1) s.step_timing = {};
        double t = step * dt;

        RealType nuN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_u);
        if (rank == 0)
        {
            std::cout << "Step " << std::setw(4) << step
                      << ": t=" << std::fixed << std::setprecision(4) << t
                      << " |u|=" << std::scientific << std::setprecision(3) << nuN
                      << " div_max=" << s.lastDivMax
                      << " cg_iter_p=";
            if (s.lastPressureIters == -1)      std::cout << "hypre";
            else if (s.lastPressureIters == -2) std::cout << "FAIL";
            else                                std::cout << s.lastPressureIters;
            std::cout << "\n" << std::defaultfloat;
        }

        if (vtuWriterU && (step % vtuEvery == 0 || step == numSteps))
            writeVtuFrame(step, t);
        if (comparisonWriter && (step % comparisonEvery == 0 || step == numSteps))
            writeComparisonFrame(step, t);
    }

    {
        // Stage times are the slowest rank's; iteration counts are global.
        const auto& t = s.step_timing;
        const double steps = numSteps > 1 ? numSteps - 1 : 1;
        double local[5] = {t.total, t.predictor, t.diffusion, t.pressure, t.corrector}, slowest[5];
        MPI_Allreduce(local, slowest, 5, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        long owned = s.numOwnedDofs, nodes = 0;
        MPI_Allreduce(&owned, &nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (rank == 0)
            std::cout << "[timing] ranks=" << numRanks << " nodes=" << nodes << " nodes/rank=" << nodes / numRanks
                      << " steps=" << long(steps) << std::fixed << std::setprecision(3)
                      << " ms/step: total=" << slowest[0] / steps << " predictor=" << slowest[1] / steps
                      << " diffusion=" << slowest[2] / steps << " pressure=" << slowest[3] / steps
                      << " corrector=" << slowest[4] / steps
                      << " | pressure_it/step=" << t.pressure_iterations / steps
                      << " ms/pressure_it=" << (t.pressure_iterations > 0 ? slowest[3] / t.pressure_iterations : 0.0)
                      << " velocity_it/step=" << t.velocity_iterations / steps << '\n' << std::defaultfloat;
    }

    double steady_change = std::numeric_limits<double>::infinity();
    if (checkMode && planar_projection && steady_window > 0)
    {
        const uint8_t* own = s.ownershipMap().data();
        const double difference = owned_sum(own, s.nodeCount,
            [u=s.d_u.data(), v=s.d_v.data(), w=s.d_w.data(), mass=s.d_massNode.data(),
             u0=d_steady_start[0].data(), v0=d_steady_start[1].data(),
             w0=d_steady_start[2].data()] __device__ (size_t i) -> double {
                const double du=u[i]-u0[i], dv=v[i]-v0[i], dw=w[i]-w0[i];
                return mass[i]*(du*du+dv*dv+dw*dw);
            });
        const double volume = owned_sum(own, s.nodeCount,
            [mass=s.d_massNode.data()] __device__ (size_t i) -> double { return mass[i]; });
        if (!(volume > 0) || cudaGetLastError() != cudaSuccess) MPI_Abort(MPI_COMM_WORLD, 1);
        steady_change = std::sqrt(difference / volume) / Uinf;
        if (rank == 0)
            std::cout << "[channel-steady] window=" << steady_window << " velocity_change/U="
                      << std::scientific << std::setprecision(8) << steady_change << '\n';
    }

    double finalRms     = -1.0;      // raw RMS vs analytic (reference tol 6e-3)
    double fluxRatio[3] = {0, 0, 0}; // Q(25/50/75%) / Q(inlet)

    // -------------------------------------------------------------------------
    // Outlet-plane validation against the analytic parabolic profile. Sums run
    // on the device over owned nodes and are then reduced over ranks, so every
    // node counts once on any rank count.
    // -------------------------------------------------------------------------
    {
        const size_t n = s.nodeCount;
        const uint8_t* own = s.ownershipMap().data();
        const auto& d_x = amr.domain().getNodeX();
        const RealType* node_x = d_x.data();
        const RealType* node_s = (crossAxis == "z") ? amr.domain().getNodeZ().data()
                                                    : amr.domain().getNodeY().data();
        const RealType* vel_u = s.d_u.data();
        const RealType* pres  = s.d_p.data();
        const double snapTol = 1e-9 * std::max(1.0, double(xMax - xMin));

        double xProbe = (profileX >= 0) ? profileX : (xMin + 0.9 * (xMax - xMin));
        double xTol   = (profileXTol > 0) ? profileXTol : 0.02 * (xMax - xMin);

        // Wall extents of the probe slab.
        const double inf = std::numeric_limits<double>::infinity();
        const double sLoL = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
            [node_x, node_s, xProbe, xTol, inf] __device__ (size_t i) -> double {
                return fabs(double(node_x[i]) - xProbe) < xTol ? double(node_s[i]) : inf;
            }, inf, thrust::minimum<double>());
        const double sHiL = thrust::transform_reduce(thrust::device,
            thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
            [node_x, node_s, xProbe, xTol, inf] __device__ (size_t i) -> double {
                return fabs(double(node_x[i]) - xProbe) < xTol ? double(node_s[i]) : -inf;
            }, -inf, thrust::maximum<double>());
        double sLo = sLoL, sHi = sHiL;
        MPI_Allreduce(&sLoL, &sLo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&sHiL, &sHi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        PoiseuilleProfile prof{targetUMax, sLo, sHi};

        const double sumSq = owned_sum(own, n,
            [node_x, node_s, vel_u, xProbe, xTol, prof] __device__ (size_t i) -> double {
                if (!(fabs(double(node_x[i]) - xProbe) < xTol)) return 0.0;
                const double err = double(vel_u[i]) - prof.analytic(double(node_s[i]));
                return err * err;
            });
        const long cnt = long(owned_sum(own, n, [node_x, xProbe, xTol] __device__ (size_t i) -> double {
            return fabs(double(node_x[i]) - xProbe) < xTol ? 1.0 : 0.0;
        }));

        double rms = cnt > 0 ? std::sqrt(sumSq / cnt) : -1.0;
        finalRms = rms;

        // Profile CSV for plotting (scripts/plot_poiseuille_profile.py):
        // cross-coordinate, solved u, analytic u at ONE fixed x (the node plane
        // nearest xProbe) -- the report's figure samples a single station, so
        // the CSV must too (the RMS above uses a slab, which is fine for the
        // norm but would smear the dot plot). Each rank sends the owned nodes
        // of that plane to rank 0.
        const double xSnap = nearest_node_plane(d_x, n, xProbe);
        cstone::DeviceVector<int> d_ids(n);
        const int localCount = int(thrust::copy_if(thrust::device,
            thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(int(n)), d_ids.data(),
            [own, node_x, xSnap, snapTol] __device__ (int i) -> bool {
                return own[i] == 1 && fabs(double(node_x[i]) - xSnap) < snapTol;
            }) - d_ids.data());
        cstone::DeviceVector<double> d_plane(2 * size_t(localCount));
        thrust::for_each(thrust::device,
            thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(localCount),
            [ids=d_ids.data(), node_s, vel_u, out=d_plane.data()] __device__ (int k) {
                out[2 * k]     = node_s[ids[k]];
                out[2 * k + 1] = vel_u[ids[k]];
            });
        std::vector<double> h_plane(2 * size_t(localCount));
        cudaMemcpy(h_plane.data(), d_plane.data(), h_plane.size() * sizeof(double), cudaMemcpyDeviceToHost);
        const int sendCount = int(h_plane.size());
        std::vector<int> counts(numRanks, 0), displs(numRanks, 0);
        MPI_Gather(&sendCount, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
        for (int r = 1; r < numRanks; ++r) displs[r] = displs[r - 1] + counts[r - 1];
        std::vector<double> allPlane(rank == 0 ? size_t(displs.back() + counts.back()) : 0);
        MPI_Gatherv(h_plane.data(), sendCount, MPI_DOUBLE, allPlane.data(), counts.data(),
                    displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            std::vector<std::pair<double,double>> pts;
            for (size_t k = 0; k + 1 < allPlane.size(); k += 2)
                pts.emplace_back(allPlane[k], allPlane[k + 1]);
            std::sort(pts.begin(), pts.end());
            std::string csvName = (vtuPrefix.empty() ? std::string("poiseuille") : vtuPrefix)
                                  + "_profile.csv";
            std::ofstream csv(csvName);
            csv << "y,u_solved,u_analytic\n";
            for (auto& p : pts)
                csv << p.first << "," << p.second << "," << prof.analytic(p.first) << "\n";
            std::cout << "Profile CSV written: " << csvName
                      << "  (single plane x=" << xSnap << ", " << pts.size() << " nodes)\n";

            // G from the solved VELOCITY: quadratic fit over the core
            // (u > 0.5 U_max -- avoids the near-wall overshoot layer), then
            // G = -mu * u'' = -2 a mu. Report this independently of solved p.
            double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, b0 = 0, b1 = 0, b2 = 0;
            for (auto& pr : pts)
            {
                if (pr.second <= 0.5 * targetUMax) continue;
                double yv = pr.first, uv = pr.second;
                double y2 = yv * yv;
                s0 += 1;  s1 += yv;  s2 += y2;  s3 += y2 * yv;  s4 += y2 * y2;
                b0 += uv; b1 += yv * uv; b2 += y2 * uv;
            }
            if (s0 >= 5)
            {
                // Solve the 3x3 normal equations [s4 s3 s2; s3 s2 s1; s2 s1 s0]
                // (a b c)^T = (b2 b1 b0)^T by Cramer's rule.
                auto det3 = [](double a11, double a12, double a13,
                               double a21, double a22, double a23,
                               double a31, double a32, double a33) {
                    return a11 * (a22 * a33 - a23 * a32)
                         - a12 * (a21 * a33 - a23 * a31)
                         + a13 * (a21 * a32 - a22 * a31);
                };
                double D  = det3(s4, s3, s2,  s3, s2, s1,  s2, s1, s0);
                double Da = det3(b2, s3, s2,  b1, s2, s1,  b0, s1, s0);
                double aFit = (D != 0) ? Da / D : 0;
                double mu = rho * nu;
                double Gfit = -2.0 * aFit * mu;
                double Gexact = 8.0 * mu * targetUMax;   // h=1-normalized below
                double Hloc = sHi - sLo;
                if (Hloc > 0) Gexact = 8.0 * mu * targetUMax / (Hloc * Hloc);
                std::cout << "G from velocity (core parabola fit, " << long(s0) << " nodes): "
                          << std::scientific << std::setprecision(4) << Gfit
                          << "   exact " << Gexact << "   ("
                          << std::fixed << std::setprecision(1)
                          << (Gexact != 0 ? 100.0 * Gfit / Gexact : 0.0) << "% of exact)\n";
            }
        }

        if (rank == 0)
        {
            std::cout << "\n========================================\n";
            std::cout << "Poiseuille profile validation\n";
            std::cout << "  probe plane x = " << std::fixed << std::setprecision(4) << xProbe
                      << " (+/- " << xTol << "),  channel x in [" << xMin << ", " << xMax << "]\n";
            std::cout << "  wall extents on plane: [" << sLo << ", " << sHi << "]  (axis " << crossAxis << ")\n";
            std::cout << "  U_max analytic = " << prof.uMax
                      << (useBodyForce ? "  (= G*H^2/(8 rho nu))" : "  (= 1.5*Uinf)") << "\n";
            std::cout << "  probe nodes    = " << cnt << "\n";
            std::cout << "  RMS error      = " << std::scientific << std::setprecision(6) << rms
                      << "  (normalized: " << (rms / prof.uMax) << ")\n";
            std::cout << "========================================\n";
        }

        // Interior-flux probe (the honest through-flow diagnostic from the pump
        // fork): Q(x*) at 25/50/75% of the channel from the SOLVED interior
        // velocity, vs Qin at the inlet plane. Cannot be faked by BC values --
        // a dead channel shows ratio ~0 here no matter what the BCs claim.
        {
            // Flux through the node plane nearest xTarget, u weighted by the
            // face areas of that plane.
            cstone::DeviceVector<RealType> d_area;
            auto planeFlux = [&](double xTarget) -> double
            {
                const double xPlane = nearest_node_plane(d_x, n, xTarget);
                const int side = xPlane < double(xMax) - snapTol ? +1 : -1;
                plane_face_areas<KeyType>(amr.domain(), xPlane, snapTol, side, 1.0, d_area);
                return owned_sum(own, n, [a=d_area.data(), vel_u] __device__ (size_t i) -> double {
                    return double(vel_u[i]) * double(a[i]);
                });
            };

            double qIn = planeFlux(double(xMin));
            const double fracs[3] = {0.25, 0.50, 0.75};
            double qAt[3];
            for (int k = 0; k < 3; ++k)
            {
                qAt[k] = planeFlux(double(xMin) + fracs[k] * double(xMax - xMin));
                fluxRatio[k] = (std::abs(qIn) > 0) ? qAt[k] / qIn : 0.0;
            }
            if (rank == 0)
            {
                std::cout << "Interior-flux probe (solved u, face-lumped):\n";
                std::cout << "  Q(inlet) = " << std::scientific << std::setprecision(4) << qIn << "\n";
                for (int k = 0; k < 3; ++k)
                    std::cout << "  Q(" << std::fixed << std::setprecision(0) << fracs[k] * 100
                              << "%) = " << std::scientific << std::setprecision(4) << qAt[k]
                              << "  ratio=" << std::fixed << std::setprecision(3)
                              << fluxRatio[k] << "\n";
                std::cout << "========================================\n";
            }
        }

        // Wikipedia plane-Poiseuille closure check: the exact solution is
        // u(y) = G/(2 mu) y(h-y) with G = -dp/dx = 8 mu U_max/h^2. Measure the
        // SOLVED pressure gradient between two developed stations (60%, 90%)
        // and compare -- validates the pressure field, not just the velocity.
        //
        // Use the moving core for this diagnostic. An empty sample is invalid,
        // not a zero pressure gradient.
        {
            auto planeMeanP = [&](double xTarget) -> double
            {
                const double xPlane = nearest_node_plane(d_x, n, xTarget);
                const double uCore = 0.5 * targetUMax;
                const double sum = owned_sum(own, n,
                    [node_x, vel_u, pres, xPlane, snapTol, uCore] __device__ (size_t i) -> double {
                        const bool core = fabs(double(node_x[i]) - xPlane) < snapTol && double(vel_u[i]) > uCore;
                        return core ? double(pres[i]) : 0.0;
                    });
                const double count = owned_sum(own, n,
                    [node_x, vel_u, xPlane, snapTol, uCore] __device__ (size_t i) -> double {
                        const bool core = fabs(double(node_x[i]) - xPlane) < snapTol && double(vel_u[i]) > uCore;
                        return core ? 1.0 : 0.0;
                    });
                return count > 0 ? sum / count : std::numeric_limits<double>::quiet_NaN();
            };
            double xA = double(xMin) + 0.60 * double(xMax - xMin);
            double xB = double(xMin) + 0.90 * double(xMax - xMin);
            double pA = planeMeanP(xA);
            double pB = planeMeanP(xB);
            double H      = sHi - sLo;
            double Gmeas  = (pA - pB) / (xB - xA);
            double Gexact = (H > 0) ? 8.0 * rho * nu * targetUMax / (H * H) : 0.0;
            if (rank == 0)
                std::cout << "Pressure-gradient from solved p (moving-core station means):\n"
                          << "  -dp/dx = " << std::scientific << std::setprecision(4) << Gmeas
                          << "   physical G = " << Gexact << "\n"
                          << "========================================\n";
        }
    }

    // Regression verdict. finalRms and fluxRatio are identical on all ranks
    // (both come from Allreduced sums), so every rank computes the same code.
    int exitCode = 0;
    if (checkMode)
    {
        bool rmsOk  = (finalRms >= 0.0) && (finalRms < rmsTol);
        bool fluxOk = true;
        for (int k = 0; k < 3; ++k)
            fluxOk = fluxOk && (std::abs(fluxRatio[k] - 1.0) <= fluxTol);
        const bool steady_ok = !planar_projection ||
            (std::isfinite(steady_change) && steady_change <= steady_tol);
        const double scaled_continuity = s.channel_continuity_rms * H / Uinf;
        const bool continuity_ok = !planar_projection || channel_state_converged(
            scaled_continuity, s.channel_boundary_balance, 0.0, continuity_tol, steady_tol);
        bool pass = rmsOk && fluxOk && (!planar_projection || channel_state_converged(
            scaled_continuity, s.channel_boundary_balance, steady_change, continuity_tol, steady_tol));
        if (rank == 0 && planar_projection)
            std::cout << "[channel-validation] continuity*H/U=" << std::scientific << scaled_continuity
                      << " boundary_balance=" << s.channel_boundary_balance
                      << " tolerance=" << continuity_tol
                      << " steady_change/U=" << steady_change << " tolerance=" << steady_tol << '\n';
        if (rank == 0)
            std::cout << "VALIDATION " << (pass ? "PASS" : "FAIL")
                      << ": RMS=" << std::scientific << std::setprecision(3) << finalRms
                      << (rmsOk ? " < " : " >= ") << rmsTol
                      << ", flux ratios " << std::fixed << std::setprecision(3)
                      << fluxRatio[0] << "/" << fluxRatio[1] << "/" << fluxRatio[2]
                      << (fluxOk ? " within " : " OUTSIDE ") << "1 +/- " << fluxTol
                      << ", steady=" << (steady_ok ? "PASS" : "FAIL")
                      << ", continuity=" << (continuity_ok ? "PASS" : "FAIL") << "\n";
        exitCode = pass ? 0 : 1;
    }

#ifdef MARS_ENABLE_HYPRE
    // Hypre objects must be destroyed while MPI is still initialized.
    s.pressureAmgSolver.reset();
    s.velocityAmgSolver.reset();
    s.velocityAmgSolverBdf2.reset();
#endif
    MPI_Finalize();
    return exitCode;
}
