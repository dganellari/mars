// Incompressible Navier-Stokes driver, Chorin projection method (Phase C).
//
// Per step, given (u^n, v^n, w^n, p^n) on nodes:
//   1) PREDICTOR (explicit):
//        q* = q^n - dt * [ (u^n.grad) q^n + (1/rho) (grad p^n)_q ]
//      for each velocity component q in {u, v, w}. The advective derivative
//      reuses the per-node + reverse-halo upwind scatter from mars_amr_advdiff
//      (B.2); the pressure-gradient component comes from B.4's SCS gradient.
//   2) IMPLICIT DIFFUSION (per component):
//        (M/dt + nu K) q** = (M/dt) q*
//      Same matrix for all three components (constant in time), assembled ONCE.
//      Three CG/Hypre solves, one per component.
//   3) PRESSURE POISSON:
//        div u** = B.3_divergence(u**, v**, w**)
//        K phi   = (rho/dt) div u**
//      K is the SAME Laplacian (gamma=1), but with Neumann BC on all faces +
//      a single Dirichlet pin to fix the constant-mode null space. SEPARATE
//      matrix from the velocity (BC treatment differs).
//   4) CORRECTOR:
//        q^{n+1} = q** - (dt/rho) (grad phi)_q
//        p^{n+1} = p^n + phi      (incremental Chorin)
//
// Test cases (selected via --bc=):
//   cavity:  Lid-driven cavity, u=lidU on z=zmax, no-slip on the other 5
//            velocity faces. Pressure has Neumann everywhere with a single
//            corner DOF pinned to 0 to remove the constant-mode null space.
//   channel: Channel flow, Dirichlet u=Uinf on x=xmin inflow, no-slip on y/z
//            tunnel walls, natural (Neumann) velocity on x=xmax outflow.
//            Pressure has Dirichlet p=0 on the entire outflow face (this is
//            the physically correct condition for channel flow and removes
//            the null space without fighting the natural inflow-to-outflow
//            pressure gradient).
// Note that with the K-path pressure solve div_max decays slowly (not at
// roundoff) because K from assembleFull is not the discrete D D^T; see the
// PressureSolveKind comment for the experimental DDT alternative.
//
// Reuses the per-node + reverseExchangeNodeHaloAdd scatter pattern validated in
// B.2/B.3/B.4 across all five scatter kernels (mass, advection x3 components,
// gradient, divergence). All velocity/pressure halo exchanges happen between
// steps explicitly; ghost-side staleness is the most common projection-driver
// bug, so we exchange whenever a field changes.


// Solver kernels, NSStepper struct, setup/runNsStep all live in the
// shared header so multiple drivers can reuse them.
#include "backend/distributed/unstructured/fem/mars_ns_solver.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"


// =============================================================================
// main: parse args, init AMR/mesh, run setup, time loop, optionally write VTU.
// =============================================================================

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    std::string meshFile;
    CvfemKernelVariant kernelVariant = CvfemKernelVariant::Tensor;
    int blockSize    = 256;
    int bucketSize   = 64;
    int maxIter      = 1000;
    double tolerance = 1e-10;
    double rho       = 1.0;
    double nu        = 0.01;
    double dt        = 0.01;
    int numSteps     = 500;
    int vtuEvery     = 10;
    std::string vtuPrefix;
    SolverKind solverKind = SolverKind::CG;
    PressureSolveKind pressureSolve = PressureSolveKind::K;
    std::string bcKind = "cavity";
    double lidU = 1.0;
    bool useLegacyGradient = true;    // SCS gradient is the stable default; --experimental-divT swaps

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg.find("--mesh=") == 0)              meshFile = arg.substr(7);
        else if (arg.find("--kernel=") == 0)
        {
            std::string v = arg.substr(9);
            if      (v == "tensor")        kernelVariant = CvfemKernelVariant::Tensor;
            else if (v == "original")      kernelVariant = CvfemKernelVariant::Original;
            else if (v == "optimized")     kernelVariant = CvfemKernelVariant::Optimized;
            else if (v == "shmem")         kernelVariant = CvfemKernelVariant::Shmem;
            else if (v == "tensor_colored")kernelVariant = CvfemKernelVariant::TensorColored;
            else if (v == "tensor_aos")    kernelVariant = CvfemKernelVariant::TensorAoS;
            else if (v == "wmma_tensor")   kernelVariant = CvfemKernelVariant::WmmaTensor;
        }
        else if (arg.find("--block-size=") == 0)   blockSize  = std::stoi(arg.substr(13));
        else if (arg.find("--bucket-size=") == 0)  bucketSize = std::stoi(arg.substr(14));
        else if (arg.find("--max-iter=") == 0)     maxIter    = std::stoi(arg.substr(11));
        else if (arg.find("--tol=") == 0)          tolerance  = std::stod(arg.substr(6));
        else if (arg.find("--rho=") == 0)          rho        = std::stod(arg.substr(6));
        else if (arg.find("--nu=") == 0)           nu         = std::stod(arg.substr(5));
        else if (arg.find("--dt=") == 0)           dt         = std::stod(arg.substr(5));
        else if (arg.find("--num-steps=") == 0)    numSteps   = std::stoi(arg.substr(12));
        else if (arg.find("--vtu-every=") == 0)    vtuEvery   = std::stoi(arg.substr(12));
        else if (arg.find("--vtu-output=") == 0)   vtuPrefix  = arg.substr(13);
        else if (arg.find("--lid-u=") == 0)        lidU       = std::stod(arg.substr(8));
        else if (arg.find("--bc=") == 0)           bcKind     = arg.substr(5);
        else if (arg.find("--solver=") == 0)
        {
            std::string v = arg.substr(9);
            if      (v == "cg")    solverKind = SolverKind::CG;
            else if (v == "hypre") solverKind = SolverKind::Hypre;
            else
            {
                if (rank == 0) std::cerr << "Error: --solver must be cg or hypre, got '" << v << "'\n";
                MPI_Finalize();
                return 1;
            }
        }
        else if (arg == "--experimental-divT") useLegacyGradient = false;
        else if (arg.find("--pressure-solve=") == 0)
        {
            std::string v = arg.substr(17);
            if      (v == "K")   pressureSolve = PressureSolveKind::K;
            else if (v == "DDT") pressureSolve = PressureSolveKind::DDT;
            else
            {
                if (rank == 0) std::cerr << "Error: --pressure-solve must be K or DDT, got '" << v << "'\n";
                MPI_Finalize();
                return 1;
            }
        }
        else if (arg[0] != '-' && meshFile.empty()) meshFile = arg;
    }

    if (meshFile.empty())
    {
        if (rank == 0)
        {
            std::cout << "Usage: " << argv[0] << " --mesh=FILE [options]\n\n"
                      << "Chorin projection NS driver (Phase C).\n\n"
                      << "Options:\n"
                      << "  --mesh=FILE           Mesh file (.mesh or .exo) [REQUIRED]\n"
                      << "  --bc=cavity|channel   BC config (default: cavity)\n"
                      << "                          cavity:  lid u=lidU on top face, no-slip elsewhere\n"
                      << "                          channel: Dirichlet u=lidU inflow on x=xmin,\n"
                      << "                                   no-slip walls, natural outflow on x=xmax\n"
                      << "  --rho=VALUE           Density (default 1.0)\n"
                      << "  --nu=VALUE            Kinematic viscosity (default 0.01)\n"
                      << "  --dt=VALUE            Timestep (default 0.01)\n"
                      << "  --num-steps=N         Number of time steps (default 500)\n"
                      << "  --lid-u=VALUE         Lid velocity for cavity (default 1.0)\n"
                      << "  --kernel=VARIANT      tensor (default), optimized, shmem, wmma_tensor, ...\n"
                      << "  --solver=KIND         cg | hypre (default cg)\n"
                      << "  --pressure-solve=KIND K | DDT (default K, matrix-free Galerkin if DDT)\n"
                      << "  --vtu-output=PREFIX   Write per-step VTU/PVTU/PVD frames\n"
                      << "  --vtu-every=N         Write a frame every N steps (default 10)\n"
                      << "  --max-iter=N          Solver max iterations (default 1000)\n"
                      << "  --tol=VALUE           Solver tolerance (default 1e-10)\n"
                      << "  --bucket-size=N       Cornerstone bucket size (default 64)\n"
                      << "  --block-size=N        CUDA block size (default 256)\n";
        }
        MPI_Finalize();
        return 1;
    }

    if (bcKind != "cavity" && bcKind != "channel")
    {
        if (rank == 0)
            std::cerr << "Error: --bc must be 'cavity' or 'channel', got '" << bcKind << "'\n";
        MPI_Finalize();
        return 1;
    }
    using KeyType  = uint64_t;
    using RealType = double;

    if (rank == 0)
    {
        const double Re = lidU * 1.0 / nu;  // L = 1 for the unit cube
        std::cout << "\n========================================\n";
        std::cout << "MARS NS Chorin projection driver (Phase C)\n";
        std::cout << "========================================\n";
        if (bcKind == "channel")
            std::cout << "BC:        channel flow, Uinf=" << lidU << " inflow on x=xmin, no-slip walls, natural outflow on x=xmax\n";
        else
            std::cout << "BC:        lid-driven cavity, u=" << lidU << " on top face\n";
        std::cout << "rho       = " << rho << "\n";
        std::cout << "nu        = " << nu  << "\n";
        std::cout << "dt        = " << dt  << "\n";
        std::cout << "Re (~L*U/nu, L=1) = " << Re << "\n";
        std::cout << "numSteps  = " << numSteps << " (T_final = " << numSteps * dt << ")\n";
        std::cout << "Mesh:    " << meshFile << "\n";
        std::cout << "Kernel:  " << CvfemHexAssembler<KeyType, RealType>::variantName(kernelVariant) << "\n";
        std::cout << "Solver:  " << (solverKind == SolverKind::CG ? "cg" : "hypre (PCG+BoomerAMG)") << "\n";
        std::cout << "Pressure solve: " << (pressureSolve == PressureSolveKind::K
                                            ? "K (CVFEM stiffness)"
                                            : "DDT (matrix-free D M^-1 D^T)") << "\n";
        std::cout << "MPI ranks: " << numRanks << "\n";
        std::cout << "========================================\n\n";
    }

    // AmrManager only for the mesh load + lazy halo/adj/coord build (maxLevels=0
    // means the mesh is frozen for this driver -- AMR-on-NS is out of scope for
    // Phase C).
    AmrManager<HexTag, KeyType, RealType>::Config amrConfig;
    amrConfig.maxLevels  = 0;
    amrConfig.blockSize  = blockSize;
    amrConfig.bucketSize = bucketSize;

    AmrManager<HexTag, KeyType, RealType> amr(amrConfig);
    amr.initialize(meshFile, rank, numRanks);
    auto initT = amr.initTimings();

    if (rank == 0)
    {
        std::cout << "Initial mesh:\n";
        std::cout << "  Elements: " << amr.domain().getElementCount() << "\n";
        std::cout << "  Nodes:    " << amr.domain().getNodeCount() << "\n";
        std::cout << "  Init breakdown (ms):\n";
        std::cout << "    domain (file read + sync):    " << std::fixed << initT.domainSyncTimeMs << "\n";
        std::cout << "    halo + node topology:         " << initT.haloTopoTimeMs   << "\n";
        std::cout << "    adjacency CSR (e2n + n2e):    " << initT.adjacencyTimeMs  << "\n";
        std::cout << "    coord cache (SFC decode):     " << initT.coordCacheTimeMs << "\n";
        std::cout << "    AMR octree state:             " << initT.octreeTimeMs     << "\n";
        std::cout << "    TOTAL:                        " << initT.totalMs          << "\n\n";
    }

    NSStepper<KeyType, RealType> s{amr.domain(), solverKind, blockSize, maxIter,
                                   RealType(tolerance), rank, numRanks};
    s.lidU = RealType(lidU);
    s.Uinf = RealType(lidU);   // reuse the --lid-u CLI flag as channel inflow speed
    s.bcKind = (bcKind == "channel") ? NSStepper<KeyType, RealType>::BCKind::Channel
                                      : NSStepper<KeyType, RealType>::BCKind::Cavity;
    s.useLegacyGradient = useLegacyGradient;
    s.pressureSolve = pressureSolve;
    setupNSStepper<KeyType, RealType>(s, RealType(nu), RealType(dt), kernelVariant);

    // IC: u = lid on top face, 0 elsewhere; p = 0.
    applyInitialCondition<KeyType, RealType>(s);

    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriterU;
    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriterP;
    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vtuWriterUmag;
    if (!vtuPrefix.empty())
    {
        // One writer per field: ParaView reads them as separate PVD timelines.
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
        // |u| magnitude as a third field for visualization.
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

    // Step 0 report: norms of the IC (|u| is the lid contribution only).
    {
        RealType nu0 = computeWeightedL2Norm<KeyType, RealType>(s, s.d_u);
        RealType nv0 = computeWeightedL2Norm<KeyType, RealType>(s, s.d_v);
        RealType nw0 = computeWeightedL2Norm<KeyType, RealType>(s, s.d_w);
        RealType np0 = computeWeightedL2Norm<KeyType, RealType>(s, s.d_p);
        if (rank == 0)
        {
            std::cout << "Step " << std::setw(4) << 0 << ": "
                      << "t=" << std::fixed << std::setprecision(4) << 0.0
                      << " |u|=" << std::scientific << std::setprecision(3) << nu0
                      << " |v|=" << nv0
                      << " |w|=" << nw0
                      << " |p|=" << np0
                      << " div_max=0.000e+00 cg_iter_p=N/A\n"
                      << std::defaultfloat;
        }
        if (vtuWriterU) writeVtuFrame(0, 0.0);
    }

    auto totalStart = std::chrono::high_resolution_clock::now();

    for (int step = 1; step <= numSteps; ++step)
    {
        runNsStep<KeyType, RealType>(s, RealType(dt), RealType(nu), RealType(rho));

        double t = step * dt;

        RealType nuN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_u);
        RealType nvN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_v);
        RealType nwN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_w);
        RealType npN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_p);

        if (rank == 0)
        {
            std::cout << "Step " << std::setw(4) << step << ": "
                      << "t=" << std::fixed << std::setprecision(4) << t
                      << " |u|=" << std::scientific << std::setprecision(3) << nuN
                      << " |v|=" << nvN
                      << " |w|=" << nwN
                      << " |p|=" << npN
                      << " div_max=" << s.lastDivMax
                      << " div_rms=" << s.lastDivRms
                      << " div_pre=" << s.lastDivMaxPre
                      << " |gp|=" << s.lastGradPRms
                      << " |gphi|=" << s.lastGradPhiRms
                      << " cg_iter_p=";
            if (s.lastPressureIters == -1)       std::cout << "hypre";
            else if (s.lastPressureIters == -2)  std::cout << "FAIL";
            else                                  std::cout << s.lastPressureIters;
            // cg_iter_uvw: report the max across the three velocity solves so
            // a single number summarizes convergence (Hypre: -1; FAIL: -2).
            int vMax = std::max({s.lastUIters, s.lastVIters, s.lastWIters});
            std::cout << " cg_iter_uvw=";
            if (vMax == -1)      std::cout << "hypre";
            else if (vMax == -2) std::cout << "FAIL";
            else                  std::cout << vMax;
            std::cout << "\n" << std::defaultfloat;
        }

        if (vtuWriterU && (step % vtuEvery == 0 || step == numSteps))
        {
            writeVtuFrame(step, t);
            if (rank == 0)
                std::cout << "  VTU: wrote " << vtuPrefix << "_{u,p,umag}_step"
                          << std::setw(4) << std::setfill('0') << step << ".pvtu\n"
                          << std::setfill(' ');
        }
    }

    auto totalEnd = std::chrono::high_resolution_clock::now();
    float totalMs = std::chrono::duration<float, std::milli>(totalEnd - totalStart).count();

    if (rank == 0)
    {
        std::cout << "\n========================================\n";
        std::cout << "NS projection run complete\n";
        std::cout << "  Time steps: " << numSteps << "\n";
        std::cout << "  Total wall: " << std::fixed << totalMs << " ms ("
                  << (totalMs / std::max(numSteps, 1)) << " ms/step)\n";
        std::cout << "========================================\n";
    }

    MPI_Finalize();
    return 0;
}
