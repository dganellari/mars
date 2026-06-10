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

#include <cmath>
#include <vector>
#include <algorithm>

// Voronoi-lumped per-node areas of a y-z cross-section plane of this quasi-2D
// mesh (one cell thick in z, two z node-planes): per node, (Voronoi y-interval)
// x dz/2. Sums to the exact face area H*dz over a full plane. ys holds the y
// coordinate of each plane node (parallel to the caller's node-id list).
// Single-rank assumption: the local node set must cover the full plane.
static std::vector<double> planeVoronoiAreas(const std::vector<double>& ys, double dz)
{
    std::vector<double> yu(ys);
    std::sort(yu.begin(), yu.end());
    yu.erase(std::unique(yu.begin(), yu.end(),
                         [](double a, double b) { return std::abs(a - b) < 1e-12; }),
             yu.end());
    std::vector<double> areas(ys.size(), 0.0);
    if (yu.size() < 2) return areas;
    for (size_t k = 0; k < ys.size(); ++k)
    {
        size_t j = std::lower_bound(yu.begin(), yu.end(), ys[k] - 1e-12) - yu.begin();
        if (j >= yu.size()) j = yu.size() - 1;
        double lo = (j == 0) ? yu.front() : 0.5 * (yu[j - 1] + yu[j]);
        double hi = (j + 1 >= yu.size()) ? yu.back() : 0.5 * (yu[j] + yu[j + 1]);
        areas[k] = (hi - lo) * dz * 0.5;
    }
    return areas;
}

// Analytic plane-Poiseuille profile, normalized cross-section coordinate
// eta in [-1, 1] (eta=0 center, |eta|=1 walls). U_max = 1.5 * U_avg.
struct PoiseuilleProfile
{
    double uMax;
    double wallLo;   // cross-section coordinate at one wall
    double wallHi;   // cross-section coordinate at the other wall

    double eta(double s) const
    {
        double mid  = 0.5 * (wallLo + wallHi);
        double half = 0.5 * (wallHi - wallLo);
        return half > 0 ? (s - mid) / half : 0.0;
    }
    double analytic(double s) const
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
    SolverKind        solverKind    = SolverKind::CG;
    // DDT (matrix-free D M^-1 D^T + Jacobi-PCG) is the validated channel
    // pressure path; the K-path on this channel hits pAp<=0 (cg_iter_p=FAIL).
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

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if      (arg.find("--mesh=") == 0)         meshFile   = arg.substr(7);
        else if (arg.find("--nu=") == 0)           nu         = std::stod(arg.substr(5));
        else if (arg.find("--dt=") == 0)           dt         = std::stod(arg.substr(5));
        else if (arg.find("--rho=") == 0)          rho        = std::stod(arg.substr(6));
        else if (arg.find("--num-steps=") == 0)    numSteps   = std::stoi(arg.substr(12));
        else if (arg.find("--vtu-every=") == 0)    vtuEvery   = std::stoi(arg.substr(12));
        else if (arg.find("--vtu-output=") == 0)   vtuPrefix  = arg.substr(13);
        else if (arg.find("--lid-u=") == 0)        Uinf       = std::stod(arg.substr(8));
        else if (arg.find("--uinf=") == 0)         Uinf       = std::stod(arg.substr(7));
        else if (arg.find("--profile-x=") == 0)    profileX   = std::stod(arg.substr(12));
        else if (arg.find("--profile-xtol=") == 0) profileXTol= std::stod(arg.substr(15));
        else if (arg.find("--cross-axis=") == 0)   crossAxis  = arg.substr(13);
        else if (arg.find("--drive=") == 0)        driveMode  = arg.substr(8);
        else if (arg.find("--body-force-x=") == 0) bodyForceX = std::stod(arg.substr(15));
        else if (arg == "--no-seed-interior")      seedInterior = false;
        else if (arg == "--no-opening-flux-source") openingFluxSource = false;
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

    if (meshFile.empty())
    {
        if (rank == 0)
            std::cout << "Usage: " << argv[0] << " --mesh=FILE [options]\n\n"
                      << "Poiseuille channel flow validation (parabolic profile).\n\n"
                      << "  --mesh=FILE         Mesh file (.exo / .mesh) [REQUIRED]\n"
                      << "  --drive=bodyforce|inlet  Flow driver (default bodyforce; textbook plane Poiseuille)\n"
                      << "                      bodyforce: constant streamwise force G, no-slip y-walls,\n"
                      << "                                 U_max=G*H^2/(8 rho nu). inlet: prescribed velocity\n"
                      << "                                 inlet+mass-conserving outlet (study; does not sustain)\n"
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
                      << "  --profile-x=X       Outlet probe plane x (default: 90% down the channel)\n"
                      << "  --profile-xtol=TOL  Half-width of the probe plane (default: one element)\n"
                      << "  --solver=cg|hypre   Linear solver (default cg)\n"
                      << "  --pressure-solve=K|DDT  Pressure operator (default DDT; K FAILs on this channel)\n"
                      << "  --vtu-output=PREFIX Write VTU/PVTU/PVD frames\n"
                      << "  --vtu-every=N       Frame every N steps (default 50)\n"
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
    s.bcKind = NSStepper<KeyType, RealType>::BCKind::Channel;
    s.useLegacyGradient = true;
    s.pressureSolve     = pressureSolve;

    setupNSStepper<KeyType, RealType>(s, RealType(nu), RealType(dt), kernelVariant);

    // Body-force drive (textbook plane Poiseuille). The channel height H is the
    // y-extent; centerline U_max = G*H^2/(8*rho*nu), so to hit a target U_max
    // (default Uinf) set G = 8*rho*nu*U_max/H^2. Computed after we know the
    // y-bbox below; stash the target here.
    const bool useBodyForce = (driveMode == "bodyforce");
    double targetUMax = Uinf;     // body-force mode aims for this centerline speed
    //
    // The shared Channel marker treats ALL six bbox faces' tangential ones as
    // walls -- including z=zmin/zmax. This mesh is ONE element thick in z (two
    // node-planes, bbox z extent ~0.066), so EVERY node lies on a z-face and
    // gets pinned no-slip -> the whole velocity field is Dirichlet -> flow is
    // frozen and div(u)=0 every step (confirmed: bc-count=all 30000 DOFs).
    //
    // Plane Poiseuille is inherently 2D: z is the extrusion (slip/symmetry)
    // direction, NOT a wall. So we rebuild the velocity BC here, marking:
    //   - inflow  x=xmin : u=(Uinf,0,0) Dirichlet
    //   - outflow x=xmax : u=(u_out,0,0) Dirichlet, mass-conserving (see below)
    //   - y-walls y=ymin,ymax : no-slip (0,0,0)
    // and leaving z-faces + interior FREE. The w-component is left unconstrained
    // (free-slip z); for a 1-cell-thick mesh w stays ~0 by symmetry anyway.
    //
    // MASS-CONSERVING OUTLET (the load-bearing fix): a velocity-inlet with a
    // natural (Neumann) outlet leaves the net discrete boundary flux nonzero --
    // the inlet injects Uinf*A_in that the divergence operator (interior SCS
    // faces only, no boundary-face term) cannot route to any sink, so the
    // pressure RHS b = -(rho/dt) div(u**) has a component OUTSIDE range(A).
    // CG then stalls (residual floors at the out-of-range size) and div_max
    // scales ~1/dt (the fixed geometric source amplified by rho/dt). Prescribing
    // u_out = Uinf * (N_in/N_out) on the outflow makes the net boundary flux
    // zero by construction, so the RHS is consistent and the projection
    // converges. Keep the p=0 pressure Dirichlet on x=xmax (removes the
    // constant mode). NOTE: a flat u_out over-constrains the exit, so the
    // parabola must develop UPSTREAM -- the profile probe stays at 90% down
    // the channel (its default), never at the outlet plane.
    {
        const auto& d_x = amr.domain().getNodeX();
        const auto& d_y = amr.domain().getNodeY();
        const auto& d_z = amr.domain().getNodeZ();
        const auto& d_ownership = s.ownershipMap();
        size_t n = s.nodeCount;

        std::vector<RealType> hx(n), hy(n), hz(n);
        std::vector<int>      h_n2d(n);
        std::vector<uint8_t>  h_own(n);
        cudaMemcpy(hx.data(),    d_x.data(),              n * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hy.data(),    d_y.data(),              n * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hz.data(),    d_z.data(),              n * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_n2d.data(), s.d_node_to_dof.data(), n * sizeof(int),      cudaMemcpyDeviceToHost);
        cudaMemcpy(h_own.data(), d_ownership.data(),     n * sizeof(uint8_t),  cudaMemcpyDeviceToHost);

        // Global bbox for face detection (single-rank: local == global).
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
        RealType eps = 1e-4 * std::max(RealType(1), yMax - yMin);
        double H = static_cast<double>(yMax - yMin);   // channel height

        std::vector<uint8_t>  hostIsBdry(s.numOwnedDofs, 0);
        std::vector<RealType> hostU(n, 0), hostV(n, 0), hostW(n, 0);

        if (useBodyForce)
        {
            // BODY-FORCE mode: ONLY the y-walls are no-slip Dirichlet. The
            // x-planes (inflow/outflow) are left FREE -- the flow is driven by
            // the streamwise body force G, not by a prescribed inlet, so there
            // is no net-boundary-flux trap. The outflow p=0 pressure Dirichlet
            // (set by setup) anchors the constant pressure mode.
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
            // G = 8 rho nu U_max / H^2 to hit the target centerline speed.
            double G = (bodyForceX >= 0)
                ? bodyForceX
                : 8.0 * rho * nu * targetUMax / (H * H);
            s.bodyForceX = RealType(G);
            // Record the achieved analytic U_max for the validation print.
            targetUMax = G * H * H / (8.0 * rho * nu);

            long nWallG = nWall;
            MPI_Allreduce(&nWall, &nWallG, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
            if (rank == 0)
                std::cout << "Drive: BODY FORCE G=" << G << " (x), H=" << H
                          << " -> analytic U_max=G*H^2/(8 rho nu)=" << targetUMax << "\n"
                          << "Quasi-2D channel BC: y-walls=" << nWallG
                          << " nodes (no-slip), x-planes + z-faces FREE\n";
        }
        else
        {
            // INLET mode -- matches the FLUYA reference (report.html input.yml)
            // that produces the validated parabola: inlet=velocity-Dirichlet
            // u=(Uinf,0,0), outlet=static-pressure p=0 with velocity FREE,
            // wall=no-slip, front_back(z)=symmetry/free. The OUTLET VELOCITY IS
            // NOT PINNED -- it develops freely, and the outflow p=0 Dirichlet
            // (applied by the shared channel pressure marker, setup) both anchors
            // the pressure null space and lets mass balance through the pressure
            // field. This is the reference config; the parabola develops upstream
            // of the free outlet, so the profile probe stays at 90% down.
            long nInflow = 0, nWall = 0;
            for (size_t i = 0; i < n; ++i)
            {
                bool onInflow = std::abs(hx[i] - xMin) < eps;
                bool onYWall  = std::abs(hy[i] - yMin) < eps || std::abs(hy[i] - yMax) < eps;
                // Outlet (x=xMax) is deliberately NOT tagged -> free velocity.
                if (!onInflow && !onYWall) continue;     // outlet + z-faces + interior free
                RealType uT = onInflow ? RealType(Uinf) : RealType(0);  // inlet wins corners
                hostU[i] = uT; hostV[i] = 0; hostW[i] = 0;
                if (h_own[i] != 1) continue;
                int dof = h_n2d[i];
                if (dof < 0 || dof >= s.numOwnedDofs) continue;
                hostIsBdry[dof] = 1;
                if (onInflow) ++nInflow; else ++nWall;
            }
            targetUMax = 1.5 * Uinf;   // developed parabola: U_max = 1.5 * mean = 1.5 Uinf

            long nInflowG = nInflow, nWallG = nWall;
            MPI_Allreduce(&nInflow, &nInflowG, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&nWall,   &nWallG,   1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
            if (rank == 0)
                std::cout << "Drive: INLET (FLUYA-ref) inflow=" << nInflowG << " (u=" << Uinf << "), "
                          << "outlet=pressure p=0 (velocity FREE), "
                          << "y-walls=" << nWallG << " (no-slip), z=symmetry/free\n";

            // Balanced opening-flux source wiring: per-node OUTWARD area
            // vectors for the two opening planes (inlet outward = -x, outlet
            // outward = +x), Voronoi-lumped from the plane node coordinates.
            // The prescribed outlet speed is sized for the area ratio; the
            // solver's per-step oScale then makes the net source exactly zero.
            if (openingFluxSource)
            {
                std::vector<size_t> inIds, outIds;
                std::vector<double> inYs, outYs;
                for (size_t i = 0; i < n; ++i)
                {
                    if (std::abs(hx[i] - xMin) < eps)
                    { inIds.push_back(i); inYs.push_back(double(hy[i])); }
                    else if (std::abs(hx[i] - xMax) < eps)
                    { outIds.push_back(i); outYs.push_back(double(hy[i])); }
                }
                double dz = double(zMax - zMin);
                auto inAreas  = planeVoronoiAreas(inYs, dz);
                auto outAreas = planeVoronoiAreas(outYs, dz);

                std::vector<RealType> aInX(n, 0), aZero(n, 0), aOutX(n, 0);
                double Ain = 0, Aout = 0;
                for (size_t k = 0; k < inIds.size(); ++k)
                { aInX[inIds[k]] = RealType(-inAreas[k]); Ain += inAreas[k]; }
                for (size_t k = 0; k < outIds.size(); ++k)
                { aOutX[outIds[k]] = RealType(+outAreas[k]); Aout += outAreas[k]; }

                auto upload = [&](cstone::DeviceVector<RealType>& d,
                                  const std::vector<RealType>& h)
                {
                    d.resize(n);
                    thrust::copy(h.begin(), h.end(), thrust::device_pointer_cast(d.data()));
                };
                upload(s.d_openInAreaX, aInX);   upload(s.d_openInAreaY, aZero);
                upload(s.d_openInAreaZ, aZero);
                upload(s.d_openOutAreaX, aOutX); upload(s.d_openOutAreaY, aZero);
                upload(s.d_openOutAreaZ, aZero);

                RealType outletU = (Aout > 0) ? RealType(Uinf * Ain / Aout) : RealType(Uinf);
                s.openInletVel[0]  = RealType(Uinf);
                s.openOutletVel[0] = outletU;
                s.useOpeningFluxSource = true;

                if (rank == 0)
                {
                    if (numRanks > 1)
                        std::cout << "WARNING: opening-flux Voronoi lumping assumes the full plane "
                                     "is rank-local; multi-rank areas may be wrong at partition cuts\n";
                    std::cout << "Opening-flux source: Ain=" << Ain << " Aout=" << Aout
                              << " Qin=" << Uinf * Ain << " outletU=" << outletU
                              << " (net zeroed per step via oScale)\n";
                }
            }
        }

        thrust::copy(hostIsBdry.begin(), hostIsBdry.end(),
                     thrust::device_pointer_cast(s.d_isBdryDof.data()));
        thrust::copy(hostU.begin(), hostU.end(), thrust::device_pointer_cast(s.d_uTarget.data()));
        thrust::copy(hostV.begin(), hostV.end(), thrust::device_pointer_cast(s.d_vTarget.data()));
        thrust::copy(hostW.begin(), hostW.end(), thrust::device_pointer_cast(s.d_wTarget.data()));
        cudaDeviceSynchronize();
    }

    applyInitialCondition<KeyType, RealType>(s);

    // Seed interior streamwise velocity (default ON) so the channel starts with
    // flow and relaxes to the parabola quickly instead of growing from rest.
    // Body-force mode seeds the MEAN velocity (2/3 U_max); inlet mode seeds Uinf.
    // With the BC rebuilt above, interior nodes are genuinely non-boundary, so
    // this takes effect.
    RealType uSeed = useBodyForce ? RealType(2.0 / 3.0 * targetUMax) : RealType(Uinf);
    if (seedInterior)
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
    if (!vtuPrefix.empty())
    {
        vtuWriterU    = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix + "_u");
        vtuWriterP    = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix + "_p");
        vtuWriterUmag = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix + "_umag");
        if (rank == 0)
            std::cout << "VTU output enabled: " << vtuPrefix << "_{u,p,umag}_step*.pvtu\n\n";
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

    if (vtuWriterU) writeVtuFrame(0, 0.0);

    // Step-0 norm: confirms the seeded IC actually populated d_u before any step.
    {
        RealType nu0 = computeWeightedL2Norm<KeyType, RealType>(s, s.d_u);
        if (rank == 0)
            std::cout << "Step    0: t=0.0000 |u|=" << std::scientific
                      << std::setprecision(3) << nu0 << " (IC)\n" << std::defaultfloat;
    }

    for (int step = 1; step <= numSteps; ++step)
    {
        runNsStep<KeyType, RealType>(s, RealType(dt), RealType(nu), RealType(rho));
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
    }

    // -------------------------------------------------------------------------
    // Outlet-plane validation against the analytic parabolic profile.
    // Pull u and node coords to host, gather the probe-plane nodes, fit the
    // wall extents from the data, and report RMS error vs the analytic parabola.
    // -------------------------------------------------------------------------
    {
        const auto& d_x = amr.domain().getNodeX();
        const auto& d_y = amr.domain().getNodeY();
        const auto& d_z = amr.domain().getNodeZ();
        size_t n = s.nodeCount;

        std::vector<RealType> h_u(n), h_x(n), h_y(n), h_z(n);
        cudaMemcpy(h_u.data(), s.d_u.data(), n * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_x.data(), d_x.data(),   n * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_y.data(), d_y.data(),   n * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_z.data(), d_z.data(),   n * sizeof(RealType), cudaMemcpyDeviceToHost);

        // Local x-extent -> reduce to global for the auto probe plane.
        RealType xMinL = *std::min_element(h_x.begin(), h_x.end());
        RealType xMaxL = *std::max_element(h_x.begin(), h_x.end());
        RealType xMin = xMinL, xMax = xMaxL;
        MPI_Allreduce(&xMinL, &xMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&xMaxL, &xMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        double xProbe = (profileX >= 0) ? profileX : (xMin + 0.9 * (xMax - xMin));
        double xTol   = (profileXTol > 0) ? profileXTol : 0.02 * (xMax - xMin);

        const std::vector<RealType>& cross = (crossAxis == "z") ? h_z : h_y;

        std::vector<std::pair<double,double>> plane;   // (cross-coord, u)
        for (size_t i = 0; i < n; ++i)
            if (std::abs(h_x[i] - xProbe) < xTol)
                plane.emplace_back(cross[i], h_u[i]);

        // Wall extents from the gathered plane (reduce min/max across ranks).
        double sLoL =  1e300, sHiL = -1e300;
        for (auto& p : plane) { sLoL = std::min(sLoL, p.first); sHiL = std::max(sHiL, p.first); }
        double sLo = sLoL, sHi = sHiL;
        MPI_Allreduce(&sLoL, &sLo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&sHiL, &sHi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        PoiseuilleProfile prof{targetUMax, sLo, sHi};

        // Local sum of squared error + node count, reduce to global RMS.
        double sumSqL = 0.0;
        long   cntL   = 0;
        for (auto& p : plane)
        {
            double err = p.second - prof.analytic(p.first);
            sumSqL += err * err;
            ++cntL;
        }
        double sumSq = 0.0; long cnt = 0;
        MPI_Allreduce(&sumSqL, &sumSq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&cntL,   &cnt,   1, MPI_LONG,   MPI_SUM, MPI_COMM_WORLD);

        double rms = cnt > 0 ? std::sqrt(sumSq / cnt) : -1.0;

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
            RealType zMinL2 = *std::min_element(h_z.begin(), h_z.end());
            RealType zMaxL2 = *std::max_element(h_z.begin(), h_z.end());
            RealType zMin2 = zMinL2, zMax2 = zMaxL2;
            MPI_Allreduce(&zMinL2, &zMin2, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(&zMaxL2, &zMax2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            double dz = double(zMax2 - zMin2);

            // Flux through the mesh x-plane nearest to xTarget: snap to the
            // closest actual node-plane, then Voronoi-lump that plane's nodes.
            auto planeFlux = [&](double xTarget) -> double
            {
                double xSnap = 1e300;
                for (size_t i = 0; i < n; ++i)
                    if (std::abs(h_x[i] - xTarget) < std::abs(xSnap - xTarget)) xSnap = h_x[i];
                double snapTol = 1e-9 * std::max(1.0, double(xMax - xMin));
                std::vector<double> ys, us;
                for (size_t i = 0; i < n; ++i)
                    if (std::abs(h_x[i] - xSnap) < snapTol)
                    { ys.push_back(double(h_y[i])); us.push_back(double(h_u[i])); }
                auto areas = planeVoronoiAreas(ys, dz);
                double q = 0;
                for (size_t k = 0; k < us.size(); ++k) q += us[k] * areas[k];
                return q;
            };

            double qIn = planeFlux(double(xMin));
            if (rank == 0)
            {
                std::cout << "Interior-flux probe (solved u, Voronoi-lumped):\n";
                std::cout << "  Q(inlet) = " << std::scientific << std::setprecision(4) << qIn << "\n";
                for (double f : {0.25, 0.50, 0.75})
                {
                    double q = planeFlux(double(xMin) + f * double(xMax - xMin));
                    std::cout << "  Q(" << std::fixed << std::setprecision(0) << f * 100
                              << "%) = " << std::scientific << std::setprecision(4) << q
                              << "  ratio=" << std::fixed << std::setprecision(3)
                              << (std::abs(qIn) > 0 ? q / qIn : 0.0) << "\n";
                }
                std::cout << "========================================\n";
            }
        }
    }

    MPI_Finalize();
    return 0;
}
