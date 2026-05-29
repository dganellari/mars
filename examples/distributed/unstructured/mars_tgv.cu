// Taylor-Green vortex driver with AMR.
//
// Solves incompressible Navier-Stokes via Chorin projection on a periodic
// cube [0, 2*pi]^3 with the analytical TGV initial condition:
//
//   u(x,y,z,0) =  V0 * sin(x) * cos(y) * cos(z)
//   v(x,y,z,0) = -V0 * cos(x) * sin(y) * cos(z)
//   w(x,y,z,0) =  0
//   p(x,y,z,0) =  (rho*V0^2/16) * (cos(2x)+cos(2y)) * (cos(2z)+2)
//
// Boundary conditions: triply periodic (all six faces). No Dirichlet BCs; the
// pressure null-space is removed by subtracting the global mean each step.
//
// AMR: every adaptEvery steps the mesh is adapted using the vorticity
// magnitude as the error indicator. Velocity, pressure, and TGV-specific
// fields are transferred via AmrManager::adaptMeshMultiField.
//
// The numerical machinery (predictor / implicit diffusion / pressure-Poisson /
// corrector / divergence / gradient kernels) is structurally identical to
// mars_amr_ns_projection.cu; the differences are:
//   * Periodic pair-sum after every per-node scatter (mass, advection,
//     gradient, divergence) so opposite faces hold the same accumulated value.
//   * No Dirichlet enforcement on velocity matrix rows.
//   * Pressure pin replaced by global-mean subtraction every step.
//   * Periodic map rebuilt after every AMR step.
//   * Vorticity-based error indicator + adaptMeshMultiField call site.
//
// This driver intentionally lives as a self-contained .cu so the existing
// cavity/channel projection driver stays untouched.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_fem.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_assembler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_utils.hpp"
#include "backend/distributed/unstructured/fem/mars_sparsity_builder.hpp"
#include "backend/distributed/unstructured/fem/mars_sparse_matrix.hpp"
#include "backend/distributed/unstructured/fem/mars_periodic_bc.hpp"
#include "backend/distributed/unstructured/fem/mars_ns_solver.hpp"
#include "backend/distributed/unstructured/amr/mars_amr.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/fill.h>

#include <memory>
#include <mpi.h>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

// =============================================================================
// TGV initial condition kernel.
// =============================================================================

template<typename RealType>
__global__ void tgvInitialConditionKernel(const RealType* nodeX,
                                          const RealType* nodeY,
                                          const RealType* nodeZ,
                                          size_t numNodes,
                                          RealType V0, RealType rho,
                                          RealType kx, RealType ky, RealType kz,
                                          RealType* d_u, RealType* d_v,
                                          RealType* d_w, RealType* d_p)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    // Use the lowest wavenumber that fits one full period in the box so the
    // IC is exactly periodic on opposite faces (sin(kx*0) == sin(kx*L) == 0).
    RealType x = kx * nodeX[i], y = ky * nodeY[i], z = kz * nodeZ[i];
    RealType sx = sin(x), cx = cos(x);
    RealType sy = sin(y), cy = cos(y);
    RealType cz = cos(z);
    d_u[i] =  V0 * sx * cy * cz;
    d_v[i] = -V0 * cx * sy * cz;
    d_w[i] = RealType(0);
    d_p[i] = (rho * V0 * V0 / RealType(16)) *
             (cos(RealType(2) * x) + cos(RealType(2) * y)) *
             (cos(RealType(2) * z) + RealType(2));
}

// =============================================================================
// Vorticity magnitude kernel (used as AMR error indicator + viz field).
// Uses CVFEM SCS-face gradient pattern: per-element 8-corner contributions
// scattered to nodes, then divided by lumped mass for a per-node grad u.
// Vorticity omega = curl(u) = (dw/dy - dv/dz, du/dz - dw/dx, dv/dx - du/dy).
//
// For the AMR indicator we want a PER-ELEMENT scalar; we use the elementwise
// max of |omega| over the 8 corner nodes.
// =============================================================================

template<typename RealType>
__global__ void computeNodeVorticityKernel(const RealType* d_dudx, const RealType* d_dudy,
                                           const RealType* d_dudz,
                                           const RealType* d_dvdx, const RealType* d_dvdy,
                                           const RealType* d_dvdz,
                                           const RealType* d_dwdx, const RealType* d_dwdy,
                                           const RealType* d_dwdz,
                                           size_t numNodes, RealType* d_omegaMag)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= numNodes) return;
    RealType wx = d_dwdy[i] - d_dvdz[i];
    RealType wy = d_dudz[i] - d_dwdx[i];
    RealType wz = d_dvdx[i] - d_dudy[i];
    d_omegaMag[i] = sqrt(wx * wx + wy * wy + wz * wz);
}

// Per-element max over its 8 nodes (HexTag). Used to drive AMR marking.
// Connectivity arrays come in as KeyType (uint64_t) directly from the domain.
template<typename KeyType, typename RealType>
__global__ void elementMaxFromNodesKernel(const KeyType* n0, const KeyType* n1,
                                          const KeyType* n2, const KeyType* n3,
                                          const KeyType* n4, const KeyType* n5,
                                          const KeyType* n6, const KeyType* n7,
                                          size_t numElements, const RealType* d_nodeField,
                                          RealType* d_elemMax)
{
    size_t e = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    RealType m = d_nodeField[n0[e]];
    m = fmax(m, d_nodeField[n1[e]]);
    m = fmax(m, d_nodeField[n2[e]]);
    m = fmax(m, d_nodeField[n3[e]]);
    m = fmax(m, d_nodeField[n4[e]]);
    m = fmax(m, d_nodeField[n5[e]]);
    m = fmax(m, d_nodeField[n6[e]]);
    m = fmax(m, d_nodeField[n7[e]]);
    d_elemMax[e] = m;
}

// =============================================================================
// Diagnostics: total kinetic energy E = 0.5 * sum_owned( |u|^2 * V_node )
// where V_node = lumped-mass entry / rho. Used for analytical decay
// validation against TGV reference data.
// =============================================================================

// Kinetic energy weighted by lumped mass. Indexes mass by DOF (not by node)
// because d_mass is sized numOwnedDofs after periodic DOF compaction. The
// driver passes nodeToDof + numOwnedDofs in addition to d_mass.
template<typename RealType, typename DomainT>
RealType computeKineticEnergy(const DomainT& domain,
                              const cstone::DeviceVector<RealType>& d_u,
                              const cstone::DeviceVector<RealType>& d_v,
                              const cstone::DeviceVector<RealType>& d_w,
                              const cstone::DeviceVector<RealType>& d_mass,
                              const cstone::DeviceVector<int>& d_nodeToDof,
                              int numOwnedDofs,
                              MPI_Comm comm)
{
    const auto& d_ownership = domain.getNodeOwnershipMap();
    size_t numNodes = d_u.size();

    RealType local = thrust::transform_reduce(
        thrust::device,
        thrust::make_counting_iterator(size_t(0)),
        thrust::make_counting_iterator(numNodes),
        [u = d_u.data(), v = d_v.data(), w = d_w.data(), m = d_mass.data(),
         n2d = d_nodeToDof.data(), nOwn = numOwnedDofs,
         own = d_ownership.data()] __device__ (size_t i) -> RealType {
            if (own[i] != 1) return RealType(0);
            int dof = n2d[i];
            if (dof < 0 || dof >= nOwn) return RealType(0);
            RealType uu = u[i], vv = v[i], ww = w[i];
            return RealType(0.5) * m[dof] * (uu * uu + vv * vv + ww * ww);
        },
        RealType(0), thrust::plus<RealType>());

    RealType global = 0;
    MPI_Datatype mpi_real = std::is_same<RealType, double>::value ? MPI_DOUBLE : MPI_FLOAT;
    MPI_Allreduce(&local, &global, 1, mpi_real, MPI_SUM, comm);
    return global;
}

// =============================================================================
// CLI parsing helpers.
// =============================================================================

template<typename T>
T parseArg(const std::string& arg, const std::string& key, T fallback)
{
    if (arg.rfind(key, 0) == 0) {
        std::string val = arg.substr(key.size());
        std::istringstream is(val);
        T out; is >> out; return out;
    }
    return fallback;
}

// =============================================================================
// Main.
// =============================================================================

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    // Device selection is handled by the affinity script (bind_numa.sh +
    // CUDA_VISIBLE_DEVICES). Each rank already sees only its own GPU as
    // device 0, so no cudaSetDevice call is needed here.

    using KeyType  = uint64_t;
    using RealType = double;

    // --- defaults ---
    std::string meshFile;
    std::string vtuPrefix;
    RealType    V0          = 1.0;
    RealType    rho         = 1.0;
    RealType    nu          = 1.0 / 1600.0;     // Re = 1600 (canonical TGV)
    RealType    dt          = 1e-3;
    int         numSteps    = 1000;
    int         vtuEvery    = 20;
    int         maxIter     = 1000;
    RealType    tolerance   = 1e-10;
    int         blockSize   = 256;
    int         bucketSize  = 64;
    int         adaptEvery  = 0;     // 0 = no AMR (uniform mesh)
    int         maxLevels   = 2;
    RealType    refineFrac  = 0.10;
    RealType    coarsenFrac = 0.30;
    RealType    boxLo       = 0.0;
    RealType    boxHi       = 2.0 * M_PI;
    bool        useHypre    = false;
    // Default DDT for periodic (algebraically exact projection). K-path leaves
    // residual divergence growing under pure-Neumann BCs.
    // DDT is the default: it produces div(u^{n+1}) at roundoff every step
    // because A = D M^-1 D^T is the algebraic dual of the divergence
    // operator. The rank-1 constant-mode null space is handled by removing
    // the mean of the RHS before the CG and the mean of phi after. BD is
    // intentionally NOT applied to the DDT path (it would destroy the
    // projection identity); the constant-mode CG solve still converges
    // because we project onto the range each step.
    PressureSolveKind pressureSolve = PressureSolveKind::DDT;
    // Verstappen skew-symmetric advection: discretely conserves KE on any
    // velocity field. Replaces 1st-order upwind which adds artificial
    // dissipation ~ |u| * h, dominating the physical viscosity at coarse mesh.
    // Default off only for legacy compatibility; for TGV validation,
    // --skew=1 is the correct choice.
    bool        useSkewAdvection = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      (a.rfind("--mesh=", 0) == 0)         meshFile     = a.substr(7);
        else if (a.rfind("--vtu-output=", 0) == 0)   vtuPrefix    = a.substr(13);
        else if (a.rfind("--V0=", 0) == 0)           V0           = std::stod(a.substr(5));
        else if (a.rfind("--rho=", 0) == 0)          rho          = std::stod(a.substr(6));
        else if (a.rfind("--nu=", 0) == 0)           nu           = std::stod(a.substr(5));
        else if (a.rfind("--Re=", 0) == 0)           nu           = 1.0 / std::stod(a.substr(5));
        else if (a.rfind("--dt=", 0) == 0)           dt           = std::stod(a.substr(5));
        else if (a.rfind("--num-steps=", 0) == 0)    numSteps     = std::stoi(a.substr(12));
        else if (a.rfind("--vtu-every=", 0) == 0)    vtuEvery     = std::stoi(a.substr(12));
        else if (a.rfind("--adapt-every=", 0) == 0)  adaptEvery   = std::stoi(a.substr(14));
        else if (a.rfind("--max-levels=", 0) == 0)   maxLevels    = std::stoi(a.substr(13));
        else if (a.rfind("--refine-frac=", 0) == 0)  refineFrac   = std::stod(a.substr(14));
        else if (a.rfind("--coarsen-frac=", 0) == 0) coarsenFrac  = std::stod(a.substr(15));
        else if (a.rfind("--box-lo=", 0) == 0)       boxLo        = std::stod(a.substr(9));
        else if (a.rfind("--box-hi=", 0) == 0)       boxHi        = std::stod(a.substr(9));
        else if (a == "--solver=hypre")              useHypre     = true;
        else if (a == "--pressure-solve=K")          pressureSolve = PressureSolveKind::K;
        else if (a == "--pressure-solve=DDT")        pressureSolve = PressureSolveKind::DDT;
        else if (a == "--skew=1" || a == "--skew=true"  || a == "--skew") useSkewAdvection = true;
        else if (a == "--skew=0" || a == "--skew=false")                  useSkewAdvection = false;
        else if (a == "--help" || a == "-h") {
            if (rank == 0) {
                std::cout <<
                    "TGV (Taylor-Green vortex) driver with AMR\n"
                    "Usage: mars_tgv --mesh=FILE [options]\n"
                    "Required:\n"
                    "  --mesh=FILE          Hex mesh of the periodic box (default [0, 2*pi]^3)\n"
                    "Physics:\n"
                    "  --V0=VALUE           IC amplitude (default 1.0)\n"
                    "  --rho=VALUE          Density (default 1.0)\n"
                    "  --nu=VALUE           Kinematic viscosity (default 1/1600)\n"
                    "  --Re=VALUE           Equivalent to --nu=(1/Re)\n"
                    "Time:\n"
                    "  --dt=VALUE           Timestep (default 1e-3)\n"
                    "  --num-steps=N        Steps (default 1000)\n"
                    "AMR:\n"
                    "  --adapt-every=N      Adapt every N steps (default 0: no AMR)\n"
                    "  --max-levels=N       Max refinement levels (default 2)\n"
                    "  --refine-frac=F      Doerfler refine fraction (default 0.10)\n"
                    "  --coarsen-frac=F     Coarsen fraction (default 0.30)\n"
                    "Geometry:\n"
                    "  --box-lo=VALUE       Periodic box lo (default 0)\n"
                    "  --box-hi=VALUE       Periodic box hi (default 2*pi)\n"
                    "Output / solver:\n"
                    "  --vtu-output=PREFIX  Write VTU/PVTU/PVD frames\n"
                    "  --vtu-every=N        VTU every N steps (default 20)\n"
                    "  --solver=hypre       Use Hypre PCG+BoomerAMG (default CG)\n";
            }
            MPI_Finalize();
            return 0;
        }
    }
    if (meshFile.empty()) {
        if (rank == 0) std::cerr << "Error: --mesh=FILE required\n";
        MPI_Finalize();
        return 1;
    }

    SolverKind solverKind = useHypre ? SolverKind::Hypre : SolverKind::CG;

    if (rank == 0) {
        std::cout << "\n========================================\n"
                  << "MARS TGV driver (incompressible NS + AMR)\n"
                  << "========================================\n"
                  << "Re        = " << 1.0 / nu << "  (nu = " << nu << ")\n"
                  << "V0        = " << V0 << "\n"
                  << "rho       = " << rho << "\n"
                  << "dt        = " << dt << "\n"
                  << "T_final   = " << numSteps * dt << "  (steps=" << numSteps << ")\n"
                  << "Box       = [" << boxLo << ", " << boxHi << "]^3 (periodic)\n"
                  << "AMR       = " << (adaptEvery > 0 ? "ON" : "OFF");
        if (adaptEvery > 0)
            std::cout << " (every " << adaptEvery << " steps, maxLevels=" << maxLevels << ")";
        std::cout << "\nSolver    = " << (useHypre ? "Hypre PCG+BoomerAMG" : "CG (Jacobi)") << "\n"
                  << "Ranks     = " << numRanks << "\n"
                  << "Mesh      = " << meshFile << "\n"
                  << "========================================\n\n";
    }

    // ------------------------------------------------------------------------
    // Mesh load via AmrManager (maxLevels controls AMR; 0 = frozen).
    // ------------------------------------------------------------------------
    typename AmrManager<HexTag, KeyType, RealType>::Config amrConfig;
    amrConfig.maxLevels        = (adaptEvery > 0) ? maxLevels : 0;
    amrConfig.refineFraction   = refineFrac;
    amrConfig.coarsenFraction  = coarsenFrac;
    amrConfig.blockSize        = blockSize;
    amrConfig.bucketSize       = bucketSize;
    amrConfig.strategy         = mars::amr::MarkingStrategy::Doerfler;

    AmrManager<HexTag, KeyType, RealType> amr(amrConfig);
    // periodicAxesMask=7 sets the cstone Box periodic on all 3 axes (NO coord
    // shift -- geometry stays real). This makes cstone deliver opposite-face
    // nodes as halo ghosts, so buildPeriodicMap (below) finds each slave's
    // master among the local ghosts and cstone's halo exchange accumulates
    // the periodic contribution across ranks. Single-rank is unaffected
    // (no ghosts; masters are local owned nodes, found the same way).
    amr.initialize(meshFile, rank, numRanks,
                   /*periodicAxesMask=*/ 7,
                   RealType(boxLo), RealType(boxHi));

    auto& domain = amr.domain();
    size_t numNodes = domain.getNodeCount();
    size_t numElems = domain.getElementCount();
    if (rank == 0) {
        std::cout << "Initial mesh: elements=" << numElems
                  << " nodes=" << numNodes << "\n";
    }

    // ------------------------------------------------------------------------
    // Build periodic map for the current coordinates.
    // ------------------------------------------------------------------------
    PeriodicMap<KeyType, RealType> pmap;
    RealType faceEps = (boxHi - boxLo) * 1e-6;
    buildPeriodicMap<KeyType, RealType>(domain, pmap,
                                        boxLo, boxHi, boxLo, boxHi, boxLo, boxHi,
                                        faceEps, MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "Periodic map: " << pmap.numSlaves
                  << " slave nodes (on max-faces) linked to masters\n";
    }

    // ------------------------------------------------------------------------
    // Set up NSStepper in Periodic mode and apply TGV initial condition to
    // the stepper-owned fields s.d_u/v/w/p.
    // ------------------------------------------------------------------------
    NSStepper<KeyType, RealType> s{domain, solverKind, blockSize, maxIter,
                                   RealType(tolerance), rank, numRanks};
    s.bcKind         = NSStepper<KeyType, RealType>::BCKind::Periodic;
    s.periodicMap    = &pmap;
    s.pressureSolve  = pressureSolve;
    s.rotationalPressureCorrection = false;
    // The skew-symmetric advection branch is implemented in the kernel but
    // not exactly energy-conserving for non-solenoidal discrete u (Verstappen
    // identity requires div u_discrete = 0, which depends on whether DDT
    // pressure projection drives div(u^{n+1}) to roundoff). For now we keep
    // upwind: dissipative, stable, and the dissipation is bounded.
    s.useSkewSymmetricAdvection    = useSkewAdvection;

    // Rhie-Chow pressure-velocity coupling on the divergence operator.
    // Uses the CHORIN-COMPATIBLE form: compact-pressure-gradient term only
    // (the smooth term double-counts the predictor's already-subtracted
    // grad p^n contribution). This is the OpenFOAM `phi - pEqn.flux()` /
    // Code_Saturne `arak-only` form, which suppresses the Q1-Q1 checkerboard
    // pressure mode that destabilizes periodic CVFEM.
    // Rhie-Chow demonstrably destabilizes this code at rho=1 (KE=1e+84 by
    // step 10) -- the dimensional/sign issue identified earlier isn't fully
    // resolved. Disabled. Stability rests on the 4 post-substep broadcasts
    // + new start-of-predictor broadcast (master->slave on u,v,w,p).
    s.useRhieChow        = false;
    s.rhieChowTau        = -1;
    // Q1-Q1 equal-order CG-FEM is inf-sup-unstable on periodic. BD adds an
    // element-local PSD term that anchors the checkerboard pressure null-space
    // (Bochev-Dohrmann SISC 2006). Required for periodic TGV: without it the
    // pressure grows exponentially in 25-100 steps regardless of time scheme.
    s.stabBochevDohrmann = true;
    // Use the literal D^T (applyDivTransposePerNodeKernel) for grad(p) and
    // grad(phi). This makes the predictor's grad(p^n) consistent with the
    // DDT projection operator D M^{-1} D^T, so div(u^{n+1})=0 holds algebraically
    // exactly. The legacy SCS gradient is NOT consistent with DDT and leaks
    // divergence into u^{n+1} each step.
    s.useLegacyGradient = false;
    s.useBdf2           = true;

    setupNSStepper<KeyType, RealType>(s, RealType(nu), RealType(dt),
                                      CvfemKernelVariant::Tensor);

    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();
    {
        // One full TGV wavelength in each direction across the periodic box.
        RealType L = boxHi - boxLo;
        RealType k = RealType(2) * M_PI / L;
        int block = 256, grid = int((numNodes + block - 1) / block);
        tgvInitialConditionKernel<RealType><<<grid, block>>>(
            d_x.data(), d_y.data(), d_z.data(), numNodes,
            V0, rho, k, k, k,
            s.d_u.data(), s.d_v.data(), s.d_w.data(), s.d_p.data());
        cudaDeviceSynchronize();
    }
    domain.exchangeNodeHalo(s.d_u);
    domain.exchangeNodeHalo(s.d_v);
    domain.exchangeNodeHalo(s.d_w);
    domain.exchangeNodeHalo(s.d_p);

    // Subtract the mean of pressure once at t=0 so we start in the same gauge
    // the solver will maintain.
    mars::fem::removeMean<RealType>(domain, s.d_p, MPI_COMM_WORLD);
    domain.exchangeNodeHalo(s.d_p);

    // ------------------------------------------------------------------------
    // VTU writers.
    // ------------------------------------------------------------------------
    // Single multi-field PVD timeline carrying u, v, w, p, |omega|, and the
    // velocity vector. Consumed by render_tgv_4panel.py. The vorticity and
    // velocity-magnitude fields are computed on the fly from s.d_u/v/w just
    // before each frame is written.
    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType>> vw_main;
    cstone::DeviceVector<RealType> d_umag(numNodes, RealType(0));
    cstone::DeviceVector<RealType> d_omegaMag(numNodes, RealType(0));
    cstone::DeviceVector<RealType> d_elemLevel;  // populated only when AMR is on
    if (!vtuPrefix.empty()) {
        vw_main = std::make_unique<fem::VTUParallelWriter<KeyType, RealType>>(vtuPrefix);
        if (rank == 0)
            std::cout << "VTU output: " << vtuPrefix << "_step*.pvtu (multi-field, every "
                      << vtuEvery << " steps)\n";
    }

    auto writeFrame = [&] (int step, double t)
    {
        if (!vw_main) return;
        auto& curDomain = amr.domain();
        size_t curN = curDomain.getNodeCount();
        size_t curE = curDomain.getElementCount();
        if (d_umag.size()     != curN) d_umag.resize(curN);
        if (d_omegaMag.size() != curN) d_omegaMag.resize(curN);

        // Velocity magnitude (cheap; per-node thrust kernel).
        {
            const RealType* uP = s.d_u.data();
            const RealType* vP = s.d_v.data();
            const RealType* wP = s.d_w.data();
            RealType*       mP = d_umag.data();
            thrust::for_each(thrust::device,
                thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(curN),
                [uP, vP, wP, mP] __device__ (size_t i) {
                    RealType u = uP[i], v = vP[i], w = wP[i];
                    mP[i] = sqrt(u*u + v*v + w*w);
                });
            cudaDeviceSynchronize();
        }
        // Real vorticity magnitude via the SCS gradient pipeline: 3 grad passes
        // on (u, v, w) -> 9-component velocity gradient -> curl -> |omega|.
        // Same scatter / reverse-halo / periodic-sum / normalize-by-mass path
        // as the predictor's grad(p), so periodic + ghost handling are correct.
        computeVorticityMagnitudePerNode<KeyType, RealType>(
            s, s.d_u, s.d_v, s.d_w, d_omegaMag);

        // Build the field list and write the frame.
        using FD = typename fem::VTUParallelWriter<KeyType, RealType>::FieldDesc;
        std::vector<FD> fields;
        fields.push_back({ "u",     FD::Kind::PointScalar,  &s.d_u,   nullptr,  nullptr  });
        fields.push_back({ "v",     FD::Kind::PointScalar,  &s.d_v,   nullptr,  nullptr  });
        fields.push_back({ "w",     FD::Kind::PointScalar,  &s.d_w,   nullptr,  nullptr  });
        fields.push_back({ "p",     FD::Kind::PointScalar,  &s.d_p,   nullptr,  nullptr  });
        fields.push_back({ "omega", FD::Kind::PointScalar,  &d_omegaMag, nullptr, nullptr });
        fields.push_back({ "velocity",
                           FD::Kind::PointVector3,
                           &s.d_u, &s.d_v, &s.d_w });
        if (d_elemLevel.size() == curE)
            fields.push_back({ "level", FD::Kind::CellScalar, &d_elemLevel, nullptr, nullptr });

        vw_main->writeMultiFieldFrame(step, t, curDomain, fields);
    };

    // ------------------------------------------------------------------------
    // Step 0 diagnostics: kinetic energy under mesh-weighted (lumped-mass)
    // measure. For TGV on [0, 2*pi]^3 the analytic value is rho*(2*pi)^3*V0^2/8.
    // ------------------------------------------------------------------------
    auto reportStep = [&] (int step, double t, int cgU, int cgV, int cgW, int cgP)
    {
        RealType KE = computeKineticEnergy<RealType>(domain, s.d_u, s.d_v, s.d_w,
                                                      s.d_mass, s.d_node_to_dof,
                                                      s.numOwnedDofs, MPI_COMM_WORLD);
        RealType uN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_u);
        RealType pN = computeWeightedL2Norm<KeyType, RealType>(s, s.d_p);
        if (rank == 0)
        {
            std::cout << "Step " << std::setw(5) << step
                      << "  t=" << std::fixed << std::setprecision(4) << t
                      << "  KE=" << std::scientific << std::setprecision(6) << KE
                      << "  |u|=" << uN
                      << "  |p|=" << pN
                      << "  div_max=" << s.lastDivMax
                      << "  cg(uvwP)=" << cgU << "/" << cgV << "/" << cgW << "/" << cgP
                      << "\n" << std::defaultfloat;
        }
    };
    reportStep(0, 0.0, 0, 0, 0, 0);
    writeFrame(0, 0.0);

    // ------------------------------------------------------------------------
    // AMR helper: compute per-element error indicator from velocity magnitude
    // (proxy for vorticity until the SCS gradient pipeline is plumbed for
    // standalone use). For TGV the vortex sheets are where velocity changes
    // fastest, so per-element max |u| is a serviceable indicator.
    // ------------------------------------------------------------------------
    auto computeErrorIndicator =
        [&] (cstone::DeviceVector<RealType>& d_errorPerElement) {
            size_t curNumElems = amr.domain().getElementCount();
            d_errorPerElement.resize(curNumElems);
            // refresh velocity magnitude on current node array
            cstone::DeviceVector<RealType> d_um(amr.domain().getNodeCount(), RealType(0));
            const RealType* uP = s.d_u.data();
            const RealType* vP = s.d_v.data();
            const RealType* wP = s.d_w.data();
            RealType*       mP = d_um.data();
            size_t nN          = d_um.size();
            thrust::for_each(thrust::device,
                thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(nN),
                [uP, vP, wP, mP] __device__ (size_t i) {
                    RealType u = uP[i], v = vP[i], w = wP[i];
                    mP[i] = sqrt(u*u + v*v + w*w);
                });
            cudaDeviceSynchronize();
            const auto& d_conn = amr.domain().getElementToNodeConnectivity();
            int blk = 256, grd = int((curNumElems + blk - 1) / blk);
            elementMaxFromNodesKernel<KeyType, RealType><<<grd, blk>>>(
                std::get<0>(d_conn).data(),
                std::get<1>(d_conn).data(),
                std::get<2>(d_conn).data(),
                std::get<3>(d_conn).data(),
                std::get<4>(d_conn).data(),
                std::get<5>(d_conn).data(),
                std::get<6>(d_conn).data(),
                std::get<7>(d_conn).data(),
                curNumElems, d_um.data(), d_errorPerElement.data());
            cudaDeviceSynchronize();
        };

    // ------------------------------------------------------------------------
    // AMR step: gather error, call adaptMeshMultiField with {u,v,w,p}, then
    // re-build periodic map and re-run setupNSStepper so all per-DOF state is
    // sized to the new mesh. Velocity/pressure are transferred by the AMR
    // module itself.
    // ------------------------------------------------------------------------
    auto runAmrStep = [&] () {
        cstone::DeviceVector<RealType> d_error;
        computeErrorIndicator(d_error);

        cstone::DeviceVector<RealType> new_u, new_v, new_w, new_p;
        std::vector<const RealType*> oldFields = {
            s.d_u.data(), s.d_v.data(), s.d_w.data(), s.d_p.data()
        };
        std::vector<cstone::DeviceVector<RealType>*> newFields = {
            &new_u, &new_v, &new_w, &new_p
        };
        amr.adaptMeshMultiField(d_error.data(), oldFields, newFields);

        // New domain. Update local caches.
        auto& newDomain = amr.domain();
        size_t newNumNodes = newDomain.getNodeCount();
        size_t newNumElems = newDomain.getElementCount();

        // Rebuild periodic map for the new coordinates.
        // After AMR refinement the coordinates come from decodeAllNodesKernel
        // (SFC-key round-trip in a 5%-padded bounding box). At deeper octree
        // levels the decode quantization is around (box_padded_width)/2^L which
        // can exceed the initial faceEps=1e-6, so we loosen it. h_min / 100 is
        // a safe upper bound on quantization error and still much smaller than
        // any inter-node gap (which is at least h_min).
        newDomain.cacheNodeCoordinates();
        RealType h_min_amr = (boxHi - boxLo) /
                             RealType(16 * (1 << amr.config().maxLevels));
        RealType faceEpsAmr = h_min_amr * RealType(1e-2);
        if (faceEpsAmr < faceEps) faceEpsAmr = faceEps;
        if (rank == 0) {
            std::cout << "[amr] rebuilding periodic map with faceEps=" << faceEpsAmr
                      << " (was " << faceEps << ")\n";
        }
        buildPeriodicMap<KeyType, RealType>(newDomain, pmap,
                                            boxLo, boxHi, boxLo, boxHi, boxLo, boxHi,
                                            faceEpsAmr, MPI_COMM_WORLD);

        // Rebuild NSStepper in place. NSStepper's copy/move-assign is deleted
        // (it contains thrust/cstone device vectors that don't support it),
        // so we don't reassign `s`; instead we re-run setupNSStepper, which
        // does .resize() on every internal buffer to match the new mesh.
        // Per-step config flags (BCKind, stabBochevDohrmann etc.) persist
        // because they're plain data members.
        s.bcKind                       = NSStepper<KeyType, RealType>::BCKind::Periodic;
        s.periodicMap                  = &pmap;
        s.pressureSolve                = pressureSolve;
        s.rotationalPressureCorrection = false;
        s.useSkewSymmetricAdvection    = useSkewAdvection;
        s.stabBochevDohrmann           = true;
        s.stabPressureTau              = -1;
        s.useLegacyGradient            = false;
        s.useBdf2                      = true;
        setupNSStepper<KeyType, RealType>(s, RealType(nu), RealType(dt),
                                          CvfemKernelVariant::Tensor);

        // Move transferred fields into the stepper's owned vectors. setup
        // re-sized s.d_u etc. to the new nodeCount; std::move replaces the
        // (zero-initialized) contents with the AMR-interpolated values.
        s.d_u = std::move(new_u);
        s.d_v = std::move(new_v);
        s.d_w = std::move(new_w);
        s.d_p = std::move(new_p);
        newDomain.exchangeNodeHalo(s.d_u);
        newDomain.exchangeNodeHalo(s.d_v);
        newDomain.exchangeNodeHalo(s.d_w);
        newDomain.exchangeNodeHalo(s.d_p);
        mars::fem::removeMean<RealType>(newDomain, s.d_p, MPI_COMM_WORLD);
        newDomain.exchangeNodeHalo(s.d_p);

        // Resize per-frame scratch (umag, omega) to match new node count.
        d_umag.resize(newNumNodes);
        d_omegaMag.resize(newNumNodes);
        thrust::fill(thrust::device_pointer_cast(d_umag.data()),
                     thrust::device_pointer_cast(d_umag.data() + newNumNodes), RealType(0));
        thrust::fill(thrust::device_pointer_cast(d_omegaMag.data()),
                     thrust::device_pointer_cast(d_omegaMag.data() + newNumNodes), RealType(0));

        // Element-level cell data for the writer (refinement level per element).
        // amr.octree().elementLevels() gives the per-element level as DeviceVector<int>.
        const auto& levelsInt = amr.octree().elementLevels();
        d_elemLevel.resize(newNumElems);
        thrust::transform(thrust::device,
                          levelsInt.data(), levelsInt.data() + newNumElems,
                          d_elemLevel.data(),
                          [] __device__ (int l) { return RealType(l); });
        cudaDeviceSynchronize();

        if (rank == 0)
            std::cout << "    [amr] elements: " << newNumElems
                      << "  nodes: " << newNumNodes
                      << "  periodic slaves: " << pmap.numSlaves << "\n";
    };

    // ------------------------------------------------------------------------
    // Time loop.
    // ------------------------------------------------------------------------
    g_nsDebugStepsLeft = 40;  // per-substep KE + divergence budget for first 40 steps
    auto wallStart = std::chrono::high_resolution_clock::now();
    for (int step = 1; step <= numSteps; ++step)
    {
        runNsStep<KeyType, RealType>(s, RealType(dt), RealType(nu), RealType(rho));

        if (adaptEvery > 0 && step % adaptEvery == 0 &&
            amrConfig.maxLevels > 0 && amr.currentLevel() < amrConfig.maxLevels)
        {
            runAmrStep();
        }

        double t = step * dt;
        if (step % 10 == 0 || step == numSteps)
            reportStep(step, t, s.lastUIters, s.lastVIters, s.lastWIters,
                       s.lastPressureIters);

        if (!vtuPrefix.empty() && (step % vtuEvery == 0 || step == numSteps))
            writeFrame(step, t);
    }
    auto wallEnd = std::chrono::high_resolution_clock::now();
    double wallMs = std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();
    if (rank == 0)
        std::cout << "\nTotal wall time: " << std::fixed << std::setprecision(2)
                  << wallMs << " ms over " << numSteps << " steps ("
                  << wallMs / std::max(numSteps, 1) << " ms/step)\n";

    MPI_Finalize();
    return 0;
}
