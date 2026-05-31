// Pump fluid-only driver: incompressible Navier-Stokes (Chorin projection) on a
// tetrahedral pump mesh, for comparison against a legacy non-MARS solver.
//
// Internal flow, named Exodus side-sets (pass the names at runtime):
//   --inlet-ss=NAME   : Dirichlet u = inletU * inward face normal
//   --outlet-ss=NAME  : Dirichlet p = 0 (outflow)
//   every other side-set : no-slip wall (u = 0)
//
// This reuses the stepper's per-side-set BC machinery (Dirichlet node
// lists + outlet pressure mask) with pump semantics: walls go into the
// no-slip wall bucket and there is NO extra bucket (an internal flow
// has no free-stream boundary). The element type is TetTag, so this is also
// the first driver to instantiate NSStepper<KeyType,RealType,TetTag>.
//
// Pressure: matrix-free DDT (D M^-1 D^T) + Jacobi-preconditioned CG. On tets
// the K-path gradient and the SCS divergence are inconsistent (near-orthogonal),
// so DDT is the consistent projection that actually zeros the divergence.
// Advection: skew-symmetric (default), 1st-order upwind, or 2nd-order
// Barth-Jespersen limited (--bj, matches the legacy reference scheme).
//
// Non-dimensional by default (rho=1, pick nu for the target Re); pass --rho /
// --nu for a dimensional (e.g. water) setup. All BCs and properties are CLI-
// driven so a case can be matched to the legacy reference without recompiling.

#include "backend/distributed/unstructured/fem/mars_ns_solver.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"
#include "backend/distributed/unstructured/utils/mars_read_exodus_mesh.hpp"

#include <unordered_map>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace mars;
using namespace mars::fem;
using namespace mars::amr;

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    // GPU selection is the launch environment's job (bind_numa.sh /
    // CUDA_VISIBLE_DEVICES), never the code's. The driver does not call
    // cudaSetDevice.

    using KeyType  = uint64_t;
    using RealType = double;

    // --- defaults ---
    std::string meshFile;
    std::string vtuPrefix;
    std::string bcMode   = "pump";    // "pump" (Exodus side-sets) or "cavity" (tet NS validation)
    std::string outletMode = "do-nothing";  // "do-nothing" (free vel + p=0 face) or "mass-conserving" (vel outflow + pin)
    // Advection scheme: "skew" (KE-conserving, default), "upwind" (1st-order),
    // or "barth-jespersen" (2nd-order limited, matches mesh developers' legacy).
    std::string advScheme  = "skew";
    bool        useRhieChow = false;   // compact RC is geometrically unsafe on tets (blows up at every tau); --rhie-chow to force on
    RealType    rhieTau    = -1;       // RC strength; <=0 -> auto dt/rho. --rhie-tau= to sweep
    std::string inletSS;   // inlet side-set name (required for --bc=pump; pass --inlet-ss=)
    std::string outletSS;  // outlet side-set name (required for --bc=pump; pass --outlet-ss=)
    RealType    inletU   = 0.5;        // inlet speed (m/s, or non-dim)
    RealType    rho      = 1.0;
    RealType    nu       = 1.0e-2;      // pick for target Re; non-dim default
    double      reqRe    = -1;          // if >0, nu is set from Re after L known
    RealType    dt       = 1.0e-3;
    int         numSteps = 200;
    int         vtuEvery = 20;
    int         maxIter  = 2000;
    RealType    tolerance = 1e-8;
    int         blockSize  = 256;
    int         bucketSize = 64;
    RealType    icPerturb  = 0.0;

    for (int i = 1; i < argc; ++i)
    {
        std::string a = argv[i];
        if      (a.rfind("--mesh=", 0) == 0)         meshFile  = a.substr(7);
        else if (a.rfind("--vtu-output=", 0) == 0)   vtuPrefix = a.substr(13);
        else if (a == "--bc=cavity")                 bcMode    = "cavity";  // tet NS validation: lid-driven cavity, no side-sets
        else if (a == "--outlet=do-nothing")         outletMode = "do-nothing";    // free velocity + p=0 face (default)
        else if (a == "--outlet=mass-conserving")    outletMode = "mass-conserving"; // velocity-Dirichlet outflow + single pin
        else if (a == "--upwind")                    advScheme  = "upwind";  // 1st-order upwind
        else if (a == "--skew")                      advScheme  = "skew";    // skew-symmetric (default)
        else if (a == "--bj")                        advScheme  = "barth-jespersen"; // 2nd-order limited
        else if (a.rfind("--advection=", 0) == 0)    advScheme  = a.substr(12); // skew|upwind|barth-jespersen
        else if (a == "--no-rhie")                   useRhieChow = false; // plain Galerkin divergence (checkerboard-prone)
        else if (a == "--rhie-chow")                 useRhieChow = true;
        else if (a.rfind("--rhie-tau=", 0) == 0)     rhieTau   = std::stod(a.substr(11));
        else if (a.rfind("--inlet-ss=", 0) == 0)     inletSS   = a.substr(11);
        else if (a.rfind("--outlet-ss=", 0) == 0)    outletSS  = a.substr(12);
        else if (a.rfind("--inlet-velocity=", 0) == 0) inletU  = std::stod(a.substr(17));
        else if (a.rfind("--rho=", 0) == 0)          rho       = std::stod(a.substr(6));
        else if (a.rfind("--nu=", 0) == 0)         { nu = std::stod(a.substr(5)); reqRe = -1; }
        else if (a.rfind("--Re=", 0) == 0)           reqRe = std::stod(a.substr(5)); // nu set after L is known
        else if (a.rfind("--dt=", 0) == 0)           dt        = std::stod(a.substr(5));
        else if (a.rfind("--num-steps=", 0) == 0)    numSteps  = std::stoi(a.substr(12));
        else if (a.rfind("--vtu-every=", 0) == 0)    vtuEvery  = std::stoi(a.substr(12));
        else if (a.rfind("--max-iter=", 0) == 0)     maxIter   = std::stoi(a.substr(11));
        else if (a.rfind("--tol=", 0) == 0)          tolerance = std::stod(a.substr(6));
        else if (a.rfind("--ic-perturb=", 0) == 0)   icPerturb = std::stod(a.substr(13));
        else if (a == "--help" || a == "-h")
        {
            if (rank == 0)
            {
                std::cout <<
                    "Pump fluid-only NS driver (tet mesh, Exodus side-sets)\n"
                    "Usage: mars_pump --mesh=FILE.exo [options]\n"
                    "  --inlet-ss=NAME      inlet side-set name (required for pump BC)\n"
                    "  --outlet-ss=NAME     outlet side-set name (required for pump BC)\n"
                    "  --inlet-velocity=V   inlet speed along x (default 0.5)\n"
                    "  --advection=NAME     skew (default) | upwind | barth-jespersen (--bj)\n"
                    "  --rho=V --nu=V       fluid properties (default rho=1, nu=1e-2)\n"
                    "  --Re=V               set nu = inletU / Re (L=1 non-dim)\n"
                    "  --dt=V --num-steps=N time stepping (default 1e-3, 200)\n"
                    "  --vtu-output=PREFIX --vtu-every=N   VTU/PVTU output\n"
                    "  --ic-perturb=F       interior IC perturbation (break symmetry)\n";
            }
            MPI_Finalize();
            return 0;
        }
    }
    if (meshFile.empty())
    {
        if (rank == 0) std::cerr << "Error: --mesh=FILE.exo required\n";
        MPI_Finalize();
        return 1;
    }
    // Pump mode needs the inlet/outlet side-set names (no baked-in defaults).
    if (bcMode == "pump" && (inletSS.empty() || outletSS.empty()))
    {
        if (rank == 0)
            std::cerr << "Error: pump BC needs --inlet-ss=NAME and --outlet-ss=NAME "
                         "(or use --bc=cavity)\n";
        MPI_Finalize();
        return 1;
    }

    if (rank == 0)
    {
        std::cout << "\n========================================\n"
                  << "MARS tet NS (Chorin projection, DDT pressure)\n"
                  << "  mode      = " << bcMode << "\n"
                  << "========================================\n"
                  << "Mesh        = " << meshFile << "\n";
        if (bcMode == "cavity")
            std::cout << "BC          = lid-driven cavity (lid u = " << inletU << ")\n";
        else
            std::cout << "Inlet  SS   = " << inletSS  << "  (u = " << inletU << ")\n"
                      << "Outlet SS   = " << outletSS << "  (p = 0)\n"
                      << "Walls       = all other side-sets (no-slip)\n";
        std::cout << "rho         = " << rho << "   nu = " << nu << "\n"
                  << "dt          = " << dt << "   steps = " << numSteps << "\n"
                  << "advection   = " << (advScheme == "skew" ? "skew-symmetric"
                                          : advScheme == "upwind" ? "1st-order upwind"
                                          : "2nd-order Barth-Jespersen") << "\n"
                  << "Rhie-Chow   = " << (useRhieChow ? "ON" : "OFF")
                  << (useRhieChow && rhieTau > 0 ? "  (tau=" : "  (tau=auto dt/rho")
                  << (useRhieChow && rhieTau > 0 ? std::to_string(rhieTau) + ")" : ")") << "\n"
                  << "MPI ranks   = " << numRanks << "\n"
                  << "========================================\n\n";
    }

    // Mesh load via AmrManager (maxLevels=0 -> frozen mesh; no AMR this driver).
    AmrManager<TetTag, KeyType, RealType>::Config amrConfig;
    amrConfig.maxLevels  = 0;
    amrConfig.blockSize  = blockSize;
    amrConfig.bucketSize = bucketSize;

    AmrManager<TetTag, KeyType, RealType> amr(amrConfig);
    amr.initialize(meshFile, rank, numRanks);

    if (rank == 0)
    {
        std::cout << "Initial mesh: elements=" << amr.domain().getElementCount()
                  << " nodes=" << amr.domain().getNodeCount() << "\n";
    }

    // Geometry length scale (bbox diagonal). If the user requested a Reynolds
    // number, set nu = U * L / Re now that L is known (the pump mesh is in
    // physical units, so a fixed nu gives an arbitrary Re -- e.g. nu=1e-2 on an
    // 8cm geometry is Re~4, near-Stokes. Re is the physical control knob).
    const auto& box = amr.domain().getBoundingBox();
    double Lx = double(box.xmax() - box.xmin());
    double Ly = double(box.ymax() - box.ymin());
    double Lz = double(box.zmax() - box.zmin());
    double Lscale = std::sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
    if (reqRe > 0 && Lscale > 0)
    {
        nu = RealType(double(inletU) * Lscale / reqRe);
        if (rank == 0)
            std::cout << "  Re=" << reqRe << " requested -> nu = U*L/Re = "
                      << nu << "  (U=" << inletU << ", L=" << Lscale << ")\n";
    }

    NSStepper<KeyType, RealType, TetTag> s{amr.domain(), SolverKind::CG, blockSize,
                                           maxIter, RealType(tolerance), rank, numRanks};
    s.Uinf          = RealType(inletU);
    s.icPerturbMag  = RealType(icPerturb);
    // DDT (matrix-free D M^-1 D^T): the corrector applies the literal D^T, so the
    // projection cancels divergence algebraically -> div(u^{n+1}) at roundoff. The
    // K-path solves phi against the Galerkin stiffness but corrects with the SCS
    // gradient; on unstructured tets those operators are nearly orthogonal, so the
    // K-path projection FAILS (div_max blows up). DDT forces useDivT=true.
    s.pressureSolve = PressureSolveKind::DDT;
    s.useLegacyGradient = false;
    // Through-flow has no wall dissipation to absorb advective noise (unlike the
    // closed cavity), so 1st-order upwind + BDF2 blows up after a few flow-
    // throughs. Skew-symmetric advection discretely conserves KE (no exponential
    // energy injection) and is the right pairing with BDF2 (the NekRS
    // combination). Default for the pump. --bj selects 2nd-order Barth-Jespersen
    // limited advection (apples-to-apples with the mesh developers' legacy code).
    using PumpAdv = NSStepper<KeyType, RealType, TetTag>::AdvScheme;
    s.advScheme = (advScheme == "upwind")          ? PumpAdv::Upwind
                : (advScheme == "barth-jespersen") ? PumpAdv::BarthJespersen
                                                   : PumpAdv::Skew;
    // Compact Rhie-Chow on the divergence operator: couples odd/even pressure
    // nodes so the projection can see (and kill) the checkerboard mode that
    // leaves div*L/U stuck. tau auto = dt/rho. --no-rhie for an A/B comparison.
    s.useRhieChow = useRhieChow;
    s.rhieChowTau = rhieTau;   // <=0 -> kernel falls back to dt/rho
    // Cavity = lid-driven cavity (geometric BC, no side-sets) -- a controlled
    // tet-NS validation case with a known flow pattern. Pump = per-side-set
    // BC machinery (Dirichlet velocity lists + outlet pressure mask).
    bool cavityMode = (bcMode == "cavity");
    s.bcKind = cavityMode
        ? NSStepper<KeyType, RealType, TetTag>::BCKind::Cavity
        : NSStepper<KeyType, RealType, TetTag>::BCKind::Pump;
    s.lidU = RealType(inletU);   // cavity lid speed reuses the inlet-velocity flag
    // Outlet treatment: do-nothing (free velocity + p=0 face, default, stable) or
    // mass-conserving (velocity-Dirichlet outflow + single pin). The latter
    // over-constrains velocity and was seen to blow up after several flow-throughs.
    s.outletDoNothing = (outletMode != "mass-conserving");

    // -------- Resolve side-sets to local node lists (pump semantics) --------
    if (!cavityMode)
    {
        ExodusSideSets ss = readExodusSideSetsTet4(meshFile, rank);

        const auto& d_l2g = amr.domain().getLocalToGlobalNodeMap();
        std::vector<KeyType> hostL2G(d_l2g.size());
        thrust::copy(thrust::device_pointer_cast(d_l2g.data()),
                     thrust::device_pointer_cast(d_l2g.data() + d_l2g.size()),
                     hostL2G.begin());
        std::unordered_map<uint64_t, int> g2l;
        g2l.reserve(hostL2G.size() * 2);
        for (size_t li = 0; li < hostL2G.size(); ++li)
            g2l.emplace(static_cast<uint64_t>(hostL2G[li]), static_cast<int>(li));

        auto resolve = [&] (const std::vector<uint64_t>& globalIds) {
            std::vector<int> local;
            local.reserve(globalIds.size());
            for (uint64_t g : globalIds)
            {
                auto it = g2l.find(g);
                if (it != g2l.end()) local.push_back(it->second);
            }
            return local;
        };

        // Inlet + outlet by name; every other side-set is a no-slip wall.
        std::vector<int> wallNodes;
        for (const auto& [name, gids] : ss.nodesByName)
        {
            if (name == inletSS || name == outletSS) continue;
            auto r = resolve(gids);
            wallNodes.insert(wallNodes.end(), r.begin(), r.end());
        }
        std::sort(wallNodes.begin(), wallNodes.end());
        wallNodes.erase(std::unique(wallNodes.begin(), wallNodes.end()), wallNodes.end());

        // On a name miss, dump the mesh's side-set names with explicit [..]
        // delimiters and lengths so a mismatch (trailing space, case, hidden
        // char) is visible without guessing. Prints the names already in the
        // mesh -- the user passed these on their own command line.
        auto dumpAvailable = [&] () {
            if (rank != 0) return;
            std::cerr << "  side-sets present in mesh (name shown in [brackets], with length):\n";
            for (const auto& [name, gids] : ss.nodesByName)
                std::cerr << "    [" << name << "]  len=" << name.size()
                          << "  nodes=" << gids.size() << "\n";
        };

        std::vector<int> inletNodes, outletNodes;
        {
            auto it = ss.nodesByName.find(inletSS);
            if (it != ss.nodesByName.end()) inletNodes = resolve(it->second);
            else if (rank == 0)
            {
                std::cerr << "WARNING: inlet side-set [" << inletSS << "] (len="
                          << inletSS.size() << ") not found in mesh\n";
                dumpAvailable();
            }
            it = ss.nodesByName.find(outletSS);
            if (it != ss.nodesByName.end()) outletNodes = resolve(it->second);
            else if (rank == 0)
            {
                std::cerr << "WARNING: outlet side-set [" << outletSS << "] (len="
                          << outletSS.size() << ") not found in mesh\n";
                dumpAvailable();
            }
        }

        // Inlet/outlet win over walls where a node is shared on an edge.
        std::vector<int> special = inletNodes;
        special.insert(special.end(), outletNodes.begin(), outletNodes.end());
        std::sort(special.begin(), special.end());
        special.erase(std::unique(special.begin(), special.end()), special.end());
        std::vector<int> wallFiltered;
        wallFiltered.reserve(wallNodes.size());
        std::set_difference(wallNodes.begin(), wallNodes.end(),
                            special.begin(), special.end(),
                            std::back_inserter(wallFiltered));

        // Map onto the BC buckets: walls=no-slip, in=inlet, out=p0,
        // extra empty (internal flow has none).
        s.wallNodes    = std::move(wallFiltered);
        s.inletNodes   = std::move(inletNodes);
        s.outletNodes  = std::move(outletNodes);
        s.extraNodes.clear();

        if (rank == 0)
        {
            std::cout << "  pump side-sets resolved on rank 0:\n"
                      << "    walls (no-slip): " << s.wallNodes.size()   << " local nodes\n"
                      << "    inlet (u set):   " << s.inletNodes.size()  << " local nodes\n"
                      << "    outlet (p=0):    " << s.outletNodes.size() << " local nodes\n";
        }

        // FOUND-vs-RESOLVED diagnostic: the Exodus side-set has G global node
        // IDs (rank 0 reads all); each rank resolves the ones it owns/ghosts via
        // g2l. If sum-of-resolved across ranks << G, the global->local map is
        // dropping nodes (likely an Exodus-original vs cstone-SFC node-ID
        // mismatch) -- the real reason the inlet/outlet are only a few nodes.
        {
            auto found = [&](const std::string& nm) -> long long {
                auto it = ss.nodesByName.find(nm);
                return it == ss.nodesByName.end() ? 0 : (long long)it->second.size();
            };
            long long gInletFound  = found(inletSS);   // same on every rank (rank0 read)
            long long gOutletFound = found(outletSS);
            long long li = (long long)s.inletNodes.size();
            long long lo = (long long)s.outletNodes.size();
            long long gInletRes = 0, gOutletRes = 0;
            MPI_Allreduce(&li, &gInletRes,  1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&lo, &gOutletRes, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            if (rank == 0)
                std::cout << "  [ss-resolve] inlet:  found " << gInletFound
                          << " global nodes, resolved " << gInletRes << " (summed over ranks)\n"
                          << "  [ss-resolve] outlet: found " << gOutletFound
                          << " global nodes, resolved " << gOutletRes << "\n"
                          << "  (resolved << found => global->local node-ID mismatch is dropping BC nodes)\n";
        }

        // -------- Inlet/outlet face geometry --------
        // Compute the global area-weighted OUTWARD normal vector (sum of each
        // triangle's 0.5*(p1-p0)x(p2-p0); Exodus winding is outward) and its
        // magnitude (= opening area) for a side-set. Each rank sums the faces it
        // can resolve to local coords; the MPI sum is partition-independent.
        amr.domain().cacheNodeCoordinates();
        const auto& d_x = amr.domain().getNodeX();
        const auto& d_y = amr.domain().getNodeY();
        const auto& d_z = amr.domain().getNodeZ();
        std::vector<RealType> hx(d_x.size()), hy(d_y.size()), hz(d_z.size());
        thrust::copy(thrust::device_pointer_cast(d_x.data()),
                     thrust::device_pointer_cast(d_x.data() + d_x.size()), hx.begin());
        thrust::copy(thrust::device_pointer_cast(d_y.data()),
                     thrust::device_pointer_cast(d_y.data() + d_y.size()), hy.begin());
        thrust::copy(thrust::device_pointer_cast(d_z.data()),
                     thrust::device_pointer_cast(d_z.data() + d_z.size()), hz.begin());

        auto faceAreaVec = [&](const std::string& nm, double out[3]) -> double {
            double nx = 0, ny = 0, nz = 0;
            auto tit = ss.trianglesByName.find(nm);
            if (tit != ss.trianglesByName.end())
                for (size_t f = 0; f + 2 < tit->second.size(); f += 3)
                {
                    auto a = g2l.find(tit->second[f]);
                    auto b = g2l.find(tit->second[f + 1]);
                    auto c = g2l.find(tit->second[f + 2]);
                    if (a == g2l.end() || b == g2l.end() || c == g2l.end()) continue;
                    int ia = a->second, ib = b->second, ic = c->second;
                    double e1x = hx[ib]-hx[ia], e1y = hy[ib]-hy[ia], e1z = hz[ib]-hz[ia];
                    double e2x = hx[ic]-hx[ia], e2y = hy[ic]-hy[ia], e2z = hz[ic]-hz[ia];
                    nx += 0.5*(e1y*e2z - e1z*e2y);
                    ny += 0.5*(e1z*e2x - e1x*e2z);
                    nz += 0.5*(e1x*e2y - e1y*e2x);
                }
            double l[3] = {nx, ny, nz};
            MPI_Allreduce(l, out, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            return std::sqrt(out[0]*out[0] + out[1]*out[1] + out[2]*out[2]);
        };

        double aIn[3], aOut[3];
        double areaIn  = faceAreaVec(inletSS,  aIn);
        double areaOut = faceAreaVec(outletSS, aOut);

        // Inlet velocity along the INWARD normal (-outward), magnitude Uinf.
        if (areaIn > 1e-30)
        {
            s.inletDirX = RealType(-aIn[0] / areaIn);
            s.inletDirY = RealType(-aIn[1] / areaIn);
            s.inletDirZ = RealType(-aIn[2] / areaIn);
        }
        // Mass-conserving outlet: a velocity-inlet + pressure-only outlet leaves
        // an unprojectable net flux (the inlet injects U*A_in with no balanced
        // sink). Make the outlet a Dirichlet OUTFLOW that removes exactly the
        // inlet volume flux: U_out = (U_in * A_in) / A_out, along the OUTWARD
        // normal. Then div(u**) is in range of the projection and flow can
        // establish a clean inlet->outlet path. The outlet keeps p=0 for the
        // pressure null-space (set in the stepper).
        if (areaOut > 1e-30)
        {
            double Uout = (double(inletU) * areaIn) / areaOut;
            s.outletU    = RealType(Uout);
            s.outletDirX = RealType(aOut[0] / areaOut);   // outward
            s.outletDirY = RealType(aOut[1] / areaOut);
            s.outletDirZ = RealType(aOut[2] / areaOut);
        }
        if (rank == 0)
            std::cout << "    inlet:  area=" << areaIn  << "  inflow dir=("
                      << s.inletDirX << "," << s.inletDirY << "," << s.inletDirZ << ")  Q=" << (double(inletU)*areaIn) << "\n"
                      << "    outlet: area=" << areaOut << "  U_out=" << s.outletU
                      << "  outflow dir=(" << s.outletDirX << "," << s.outletDirY << "," << s.outletDirZ << ")\n";

        // Fail fast on an ill-posed setup: if the inlet or outlet side-set
        // matched 0 nodes on EVERY rank, the velocity inflow / pressure
        // Dirichlet is missing, the pressure system is singular, and setup
        // would later crash deep in a GPU kernel with a cryptic error
        // (out-of-bounds -> poisoned context -> "invalid device ordinal").
        // Catch it here with a clear message and the names that were passed.
        {
            long long li = (long long)s.inletNodes.size();
            long long lo = (long long)s.outletNodes.size();
            long long gi = 0, go = 0;
            MPI_Allreduce(&li, &gi, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&lo, &go, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            if (gi == 0 || go == 0)
            {
                if (rank == 0)
                {
                    std::cerr << "\nERROR: ill-posed pump BCs -- aborting before GPU setup.\n";
                    if (gi == 0)
                        std::cerr << "  inlet side-set '" << inletSS
                                  << "' matched 0 nodes on ANY rank (need --inlet-ss=<exact Exodus name>).\n";
                    if (go == 0)
                        std::cerr << "  outlet side-set '" << outletSS
                                  << "' matched 0 nodes on ANY rank (need --outlet-ss=<exact Exodus name>).\n";
                    std::cerr << "  Check the names against `ncdump -v ss_names <mesh>`.\n" << std::endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Finalize();
                return 1;
            }
        }
    }

    setupNSStepper<KeyType, RealType, TetTag>(s, RealType(nu), RealType(dt),
                                              CvfemKernelVariant::Tensor);

    // -------- VTU output --------
    std::unique_ptr<fem::VTUParallelWriter<KeyType, RealType, TetTag>> vw;
    if (!vtuPrefix.empty())
        vw = std::make_unique<fem::VTUParallelWriter<KeyType, RealType, TetTag>>(vtuPrefix);

    auto writeFrame = [&] (int step, double t) {
        if (!vw) return;
        auto& dom = amr.domain();
        using FD = typename fem::VTUParallelWriter<KeyType, RealType, TetTag>::FieldDesc;
        std::vector<FD> fields;
        fields.push_back({ "u", FD::Kind::PointScalar, &s.d_u, nullptr, nullptr });
        fields.push_back({ "v", FD::Kind::PointScalar, &s.d_v, nullptr, nullptr });
        fields.push_back({ "w", FD::Kind::PointScalar, &s.d_w, nullptr, nullptr });
        fields.push_back({ "p", FD::Kind::PointScalar, &s.d_p, nullptr, nullptr });
        fields.push_back({ "velocity", FD::Kind::PointVector3, &s.d_u, &s.d_v, &s.d_w });
        vw->writeMultiFieldFrame(step, t, dom, fields);
    };
    writeFrame(0, 0.0);

    // Total fluid volume (sum of owned lumped masses) and the geometry length
    // scale (bbox diagonal). The pump mesh is in physical units and may be
    // sub-mm, so the raw integral norms |u|=sqrt(integral u^2 dV) and
    // div=flux/V are tiny/huge and unreadable. Report SCALE-INDEPENDENT
    // diagnostics: RMS velocity = sqrt(integral u^2 dV / V_total) ~ O(U), and a
    // non-dimensional divergence div_max * L / U.
    double Vlocal = thrust::reduce(
        thrust::device,
        thrust::device_pointer_cast(s.d_mass.data()),
        thrust::device_pointer_cast(s.d_mass.data() + s.numOwnedDofs),
        double(0), thrust::plus<double>());
    double Vtotal = 0;
    MPI_Allreduce(&Vlocal, &Vtotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // Lscale already computed from the domain bbox before setup.
    if (rank == 0)
        std::cout << "  geometry: V_total=" << Vtotal << "  L(bbox diag)=" << Lscale
                  << "  Re=U*L/nu=" << (inletU * Lscale / nu) << "\n";

    // -------- Time loop --------
    // Flow-through time = L/U: the time for fluid to traverse the geometry once.
    // Developed flow needs SEVERAL flow-throughs; report progress in those units
    // plus the steady-state residual d(u_rms) so convergence is visible.
    double tFlow = (inletU > 0 && Lscale > 0) ? (Lscale / double(inletU)) : 1.0;
    double prevURms = 0.0;
    auto wallStart = std::chrono::high_resolution_clock::now();
    for (int step = 1; step <= numSteps; ++step)
    {
        runNsStep<KeyType, RealType, TetTag>(s, RealType(dt), RealType(nu), RealType(rho));

        if (step % 10 == 0 || step == numSteps)
        {
            RealType uN = computeWeightedL2Norm<KeyType, RealType, TetTag>(s, s.d_u);
            RealType pN = computeWeightedL2Norm<KeyType, RealType, TetTag>(s, s.d_p);
            // Scale-independent: RMS velocity ~ O(U), non-dim divergence.
            double uRms = (Vtotal > 0) ? double(uN) / std::sqrt(Vtotal) : 0.0;
            double uFrac = (inletU > 0) ? uRms / double(inletU) : uRms;       // u_rms as a fraction of inlet U
            double dURms = uRms - prevURms;                                   // steady-state residual
            prevURms = uRms;
            // Peak INTERIOR velocity (excludes the pinned inlet/outlet/wall DOFs):
            // does the inlet jet propagate into the domain? u_max/U ~ O(1) near
            // the jet means flow IS entering even if the volume-average u_rms is
            // tiny (small port into a large chamber -> bulk is nearly stagnant).
            RealType umx = maxOwnedInteriorAbs<KeyType, RealType, TetTag>(s, s.d_u);
            RealType vmx = maxOwnedInteriorAbs<KeyType, RealType, TetTag>(s, s.d_v);
            RealType wmx = maxOwnedInteriorAbs<KeyType, RealType, TetTag>(s, s.d_w);
            double uMax = std::sqrt(double(umx)*double(umx) + double(vmx)*double(vmx) + double(wmx)*double(wmx));
            double divND = (inletU > 0 && Lscale > 0)
                           ? double(s.lastDivMax) * Lscale / inletU : double(s.lastDivMax);
            // RC-flux divergence (the operator RC actually zeros); only meaningful
            // with --rhie-chow. Plain div*L/U stays high by construction.
            double divRCnd = (inletU > 0 && Lscale > 0)
                           ? double(s.lastDivRC) * Lscale / inletU : double(s.lastDivRC);
            double flowThroughs = (step * dt) / tFlow;
            if (rank == 0)
            {
                std::cout << "Step " << std::setw(5) << step
                          << "  ft=" << std::fixed << std::setprecision(2) << flowThroughs  // flow-throughs elapsed
                          << "  u_rms=" << std::scientific << std::setprecision(3) << uRms
                          << "  u_max=" << std::fixed << std::setprecision(3) << uMax        // peak interior speed
                          << "  uMax/U=" << std::setprecision(2) << (inletU > 0 ? uMax/double(inletU) : uMax)
                          << "  d(u_rms)=" << std::scientific << std::setprecision(2) << dURms
                          << "  div*L/U=" << std::fixed << std::setprecision(2) << divND;
                if (useRhieChow)
                    std::cout << "  divRC*L/U=" << std::fixed << std::setprecision(2) << divRCnd;
                std::cout << "  cg_p=" << s.lastPressureIters
                          << "  cg_uvw=" << s.lastUIters
                          << "\n" << std::defaultfloat;
            }
        }
        if (!vtuPrefix.empty() && (step % vtuEvery == 0 || step == numSteps))
            writeFrame(step, step * dt);
    }
    auto wallEnd = std::chrono::high_resolution_clock::now();
    double wallMs = std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();
    if (rank == 0)
        std::cout << "\nPump run complete: " << numSteps << " steps, "
                  << std::fixed << std::setprecision(1) << wallMs << " ms ("
                  << wallMs / std::max(numSteps, 1) << " ms/step)\n";

    MPI_Finalize();
    return 0;
}
