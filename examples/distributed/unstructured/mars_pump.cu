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

// Pump uses a FORKED NS solver pinned to the verified-good pump state (commit
// 8bb0402, 2026-06-02: 1-rank==4-rank bit-identical, stable jet to ft=26). The
// shared mars_ns_solver.hpp was heavily modified by the TGV multi-rank periodic
// work (55 commits, +1587 lines) which regressed the pump's DDT pressure solve
// (cg_p diverges). Keeping a separate header decouples them; reconcile to one
// solver later. TGV keeps using mars_ns_solver.hpp unchanged.
#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"
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
#include <cstdlib>

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
    bool        useBdf2     = true;    // --bdf1 forces BDF1/Chorin (1st-order time, more stable for explicit advection)
    bool        useRhieChow = false;   // compact RC is geometrically unsafe on tets (blows up at every tau); --rhie-chow to force on
    RealType    rhieTau    = -1;       // RC strength; <=0 -> auto dt/rho. --rhie-tau= to sweep
    bool        useVMSStab = false;    // Nalu-Wind VMS pressure stab (NOT Rhie-Chow); --vms-stab. EXPERIMENTAL, validate.
    std::string inletSS;   // inlet side-set name (required for --bc=pump; pass --inlet-ss=)
    std::string outletSS;  // outlet side-set name (required for --bc=pump; pass --outlet-ss=)
    // Inlet velocity follows each node's OWN local surface normal (area-weighted
    // from the faces touching it) instead of one global averaged direction. On a
    // curved inlet the per-node form keeps the prescribed velocity normal to the
    // surface, matching the reference solver. --no-inlet-pernode-normal falls
    // back to the old single global normal.
    bool inletPernodeNormal = true;
    // Flip the inlet normal sign. The Exodus face winding determines whether the
    // computed area-normal is inward or outward; the default negates to inward.
    // If the vectors come out OUTWARD on your mesh, pass --inlet-flip-normal.
    bool inletFlipNormal = false;
    // Pump default: start from REST and let the inlet drive the flow. The legacy
    // uniform (Uinf,0,0) IC seeds spurious +x flow in both tanks; --pump-uniform-ic
    // restores it for comparison.
    bool pumpUniformIC = false;
    RealType    inletU   = 0.5;        // inlet speed (m/s)
    RealType    rho      = 1000.0;     // water (kg/m^3)
    RealType    nu       = 1.0e-6;     // water kinematic viscosity (m^2/s) = mu/rho
    double      reqRe    = -1;          // legacy: if >0, override nu from a bbox-based Re
    RealType    dt       = 1.0e-3;     // initial/fixed dt; capped by --cfl if set
    double      cflMax   = -1;         // >0 enables adaptive dt: cap advective CFL (uMax*dt/dx) at this
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
        else if (a == "--bdf1")                      useBdf2    = false; // diagnostic: BDF1/Chorin instead of BDF2/EXT2
        else if (a.rfind("--advection=", 0) == 0)    advScheme  = a.substr(12); // skew|upwind|barth-jespersen
        else if (a == "--no-rhie")                   useRhieChow = false; // plain Galerkin divergence (checkerboard-prone)
        else if (a == "--rhie-chow")                 useRhieChow = true;
        else if (a.rfind("--rhie-tau=", 0) == 0)     rhieTau   = std::stod(a.substr(11));
        else if (a == "--vms-stab")                  useVMSStab = true;  // Nalu VMS pressure stab (tet-only, experimental)
        else if (a.rfind("--inlet-ss=", 0) == 0)     inletSS   = a.substr(11);
        else if (a.rfind("--outlet-ss=", 0) == 0)    outletSS  = a.substr(12);
        else if (a == "--inlet-pernode-normal")      inletPernodeNormal = true;
        else if (a == "--no-inlet-pernode-normal")   inletPernodeNormal = false; // global averaged normal (legacy)
        else if (a == "--inlet-flip-normal")         inletFlipNormal = true;     // flip if vectors come out outward
        else if (a == "--pump-uniform-ic")           pumpUniformIC = true;       // legacy free-stream IC (default: start from rest)
        else if (a.rfind("--inlet-velocity=", 0) == 0) inletU  = std::stod(a.substr(17));
        else if (a.rfind("--rho=", 0) == 0)          rho       = std::stod(a.substr(6));
        else if (a.rfind("--nu=", 0) == 0)         { nu = std::stod(a.substr(5)); reqRe = -1; }
        else if (a.rfind("--Re=", 0) == 0)           reqRe = std::stod(a.substr(5)); // nu set after L is known
        else if (a.rfind("--dt=", 0) == 0)           dt        = std::stod(a.substr(5));
        else if (a.rfind("--cfl=", 0) == 0)          cflMax    = std::stod(a.substr(6)); // adaptive dt: cap advective CFL
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
                    "  --inlet-flip-normal  flip inlet normal sign (use if vectors come out OUTWARD)\n"
                    "  --outlet-ss=NAME     outlet side-set name (required for pump BC)\n"
                    "  --inlet-velocity=V   inlet speed along x (default 0.5)\n"
                    "  --no-inlet-pernode-normal  use one global inlet normal (default: per-node)\n"
                    "  --advection=NAME     skew (default) | upwind | barth-jespersen (--bj)\n"
                    "  --rho=V --nu=V       physical fluid properties (default water: rho=1000, nu=1e-6)\n"
                    "  --Re=V               LEGACY: override nu = inletU*L_bbox/Re. L_bbox is the\n"
                    "                       whole-geometry diagonal, NOT the passage scale, so this Re\n"
                    "                       does NOT match a physically-defined Re. Prefer --nu/--rho.\n"
                    "  --dt=V --num-steps=N time stepping (default 1e-3, 200)\n"
                    "  --cfl=C              adaptive dt: cap advective CFL at C (~0.5 for BJ+BDF2)\n"
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
                  << "VMS-stab    = " << (useVMSStab ? "ON (Nalu, tet-only, EXPERIMENTAL)" : "OFF") << "\n"
                  << "MPI ranks   = " << numRanks << "\n"
                  << "========================================\n\n";
    }

    // Mesh load via AmrManager (maxLevels=0 -> frozen mesh; no AMR this driver).
    AmrManager<TetTag, KeyType, RealType>::Config amrConfig;
    amrConfig.maxLevels  = 0;
    amrConfig.blockSize  = blockSize;
    amrConfig.bucketSize = bucketSize;
    // Keep the Exodus side sets on-device so the pump BC node lists come from the
    // domain (defaults false -> sideSetNodes() would be empty without this).
    amrConfig.storeSideSets = true;

    AmrManager<TetTag, KeyType, RealType> amr(amrConfig);
    amr.initialize(meshFile, rank, numRanks);

    if (rank == 0)
    {
        std::cout << "Initial mesh: elements=" << amr.domain().getElementCount()
                  << " nodes=" << amr.domain().getNodeCount() << "\n";
    }

    // Geometry length scale (bbox diagonal), used for diagnostics (tFlow, div*L/U).
    const auto& box = amr.domain().getBoundingBox();
    double Lx = double(box.xmax() - box.xmin());
    double Ly = double(box.ymax() - box.ymin());
    double Lz = double(box.zmax() - box.zmin());
    double Lscale = std::sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
    // LEGACY --Re override: nu is normally the physical value (--nu, default water
    // 1e-6). --Re forces nu = U*L_bbox/Re, but L_bbox is the whole-geometry diagonal,
    // not the passage scale -- this Re is NOT a physical Reynolds number and makes the
    // fluid far too viscous (Re=100 on a real pump -> ~50x too viscous -> velocity is
    // capped far below the true flow). Kept only for old non-dimensional runs.
    if (reqRe > 0 && Lscale > 0)
    {
        nu = RealType(double(inletU) * Lscale / reqRe);
        if (rank == 0)
            std::cout << "  WARNING --Re=" << reqRe << " (legacy) overrides physical nu -> nu = U*L_bbox/Re = "
                      << nu << " (U=" << inletU << ", L_bbox=" << Lscale << "). This is NOT a physical Re; "
                      << "use --nu for a real fluid.\n";
    }

    NSStepper<KeyType, RealType, TetTag> s{amr.domain(), SolverKind::CG, blockSize,
                                           maxIter, RealType(tolerance), rank, numRanks};
    s.Uinf          = RealType(inletU);
    s.pumpZeroIC    = !pumpUniformIC;   // default: start from rest
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
    s.useBdf2 = useBdf2;
    s.useRhieChow = useRhieChow;
    s.useVMSStab  = useVMSStab;   // Nalu VMS pressure stabilization (tet-only)
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
        // The face-normal/area path below still needs the reader's per-triangle
        // coords (the domain side-set API stores node lists only, no faces yet),
        // so keep the reader read alive for that path and the global-found
        // diagnostics. The node LISTS (wall/inlet/outlet buckets) now come from
        // the domain's on-device side sets instead of an in-driver resolve.
        ExodusSideSets ss = readExodusSideSetsTet4(meshFile, rank);

        // Copy a domain side-set node list (device) to a host std::vector<int>.
        auto toHost = [](const cstone::DeviceVector<int>& d) {
            std::vector<int> h(d.size());
            if (!h.empty())
                thrust::copy(thrust::device_pointer_cast(d.data()),
                             thrust::device_pointer_cast(d.data() + d.size()), h.begin());
            return h;
        };

        // Inlet + outlet by name; every other side-set is a no-slip wall. Node
        // lists are pre-resolved on-device by the domain (opt-in storeSideSets).
        std::vector<int> wallNodes;
        for (const auto& name : amr.domain().sideSetNames())
        {
            if (name == inletSS || name == outletSS) continue;
            auto r = toHost(amr.domain().sideSetNodes(name));
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
            for (const auto& name : amr.domain().sideSetNames())
                std::cerr << "    [" << name << "]  len=" << name.size()
                          << "  nodes=" << amr.domain().sideSetNodeKeys(name).size() << "\n";
        };

        std::vector<int> inletNodes, outletNodes;
        {
            if (amr.domain().hasSideSet(inletSS)) inletNodes = toHost(amr.domain().sideSetNodes(inletSS));
            else if (rank == 0)
            {
                std::cerr << "WARNING: inlet side-set [" << inletSS << "] (len="
                          << inletSS.size() << ") not found in mesh\n";
                dumpAvailable();
            }
            if (amr.domain().hasSideSet(outletSS)) outletNodes = toHost(amr.domain().sideSetNodes(outletSS));
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
        // IDs (rank 0 reads all); each rank resolves the ones it has locally via
        // the coord->SFC-key lookup. Summed-resolved across ranks should be >= G
        // (shared nodes resolve on >1 rank). If it is << G, nodes are being
        // dropped (SFC-key mismatch) -- watch this after the index-space fix.
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

        // [diag] Per-rank inlet ownership audit (MARS_PUMP_BC_DEBUG=1). Confirms
        // the multi-rank BC failure mode: each inlet node resolves as a GHOST on
        // the ranks that hold it and is OWNED by none -> tag()'s ownership gate
        // marks zero inlet Dirichlet DOFs. For each rank print:
        //   coords  = inlet coords the reader handed this rank (global, all ranks)
        //   resolved= inlet nodes this rank found in its SFC map (owned or ghost)
        //   ownedOf = of those, how many this rank actually OWNS (hostOwn==1)
        // If sum(ownedOf) over ranks < #global-inlet-nodes, some inlet nodes are
        // owned by no rank (the bug). Healthy: sum(ownedOf) == #inlet-nodes.
        if (std::getenv("MARS_PUMP_BC_DEBUG"))
        {
            const auto& d_own = amr.domain().getNodeOwnershipMap();
            std::vector<uint8_t> hostOwn(amr.domain().getNodeCount(), 0);
            thrust::copy(thrust::device_pointer_cast(d_own.data()),
                         thrust::device_pointer_cast(d_own.data() + d_own.size()),
                         hostOwn.begin());
            long long coords   = (long long)ss.nodeCoordsByName[inletSS].size();
            long long resolved = (long long)s.inletNodes.size();
            long long ownedOf  = 0;
            for (int li : s.inletNodes)
                if (li >= 0 && (size_t)li < hostOwn.size() && hostOwn[li] == 1) ++ownedOf;
            // Serialize the per-rank lines so they do not interleave.
            for (int r = 0; r < numRanks; ++r)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                if (r == rank)
                    std::cerr << "[bc-debug] rank " << rank << " inlet: coords=" << coords
                              << " resolved=" << resolved << " ownedOf=" << ownedOf << "\n";
            }
            MPI_Barrier(MPI_COMM_WORLD);
            long long gOwnedOf = 0;
            MPI_Allreduce(&ownedOf, &gOwnedOf, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            if (rank == 0)
                std::cerr << "[bc-debug] inlet owned-of-resolved summed over ranks = "
                          << gOwnedOf << " (must equal the #global inlet nodes; 0 => globally unowned)\n";
        }

        // -------- Inlet/outlet face geometry --------
        // Compute the global area-weighted OUTWARD normal vector (sum of each
        // triangle's 0.5*(p1-p0)x(p2-p0); Exodus winding is outward) and its
        // magnitude (= opening area) for a side-set. Each rank sums the faces it
        // OWNS (first node owned here); the MPI sum is partition-independent.
        // Coords come from the reader's per-triangle coords, so no node-coord
        // cache / host copy is needed here.

        // Inlet/outlet face normal, summed over faces and MPI-reduced. Compute
        // each triangle's area-normal directly from the reader's per-triangle
        // coords (triangleCoordsByName, 3 (x,y,z) per face). To avoid counting a
        // face once per rank under the Allreduce(SUM), only sum a face if its
        // first node is LOCAL to and OWNED by this rank (resolve coord->SFC local,
        // check ownership) -- so each face is counted by exactly one rank.
        std::vector<uint8_t> hostOwnFA(amr.domain().getNodeCount(), 0);
        {
            const auto& d_own = amr.domain().getNodeOwnershipMap();
            thrust::copy(thrust::device_pointer_cast(d_own.data()),
                         thrust::device_pointer_cast(d_own.data() + d_own.size()),
                         hostOwnFA.begin());
        }
        auto faceAreaVec = [&](const std::string& nm, double out[3]) -> double {
            double nx = 0, ny = 0, nz = 0;
            auto tit = ss.triangleCoordsByName.find(nm);
            if (tit != ss.triangleCoordsByName.end())
            {
                // Resolve all triangle nodes to SFC local ids in ONE batch via
                // the shared domain resolver (returns -1 for nodes not local to
                // this rank, preserving alignment). Count each face only on the
                // rank that OWNS its first node, so the MPI sum is once-per-face.
                std::vector<int> triLocal =
                    amr.domain().resolveSideSetNodesToLocalKeepMisses(tit->second);
                for (size_t f = 0; f + 2 < tit->second.size(); f += 3)
                {
                    int la = triLocal[f];
                    if (la < 0 || hostOwnFA[la] != 1) continue;   // not owned here
                    const auto& A = tit->second[f];
                    const auto& B = tit->second[f + 1];
                    const auto& C = tit->second[f + 2];
                    double e1x = B[0]-A[0], e1y = B[1]-A[1], e1z = B[2]-A[2];
                    double e2x = C[0]-A[0], e2y = C[1]-A[1], e2z = C[2]-A[2];
                    nx += 0.5*(e1y*e2z - e1z*e2y);
                    ny += 0.5*(e1z*e2x - e1x*e2z);
                    nz += 0.5*(e1x*e2y - e1y*e2x);
                }
            }
            double l[3] = {nx, ny, nz};
            MPI_Allreduce(l, out, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            return std::sqrt(out[0]*out[0] + out[1]*out[1] + out[2]*out[2]);
        };

        double aIn[3], aOut[3];
        double areaIn  = faceAreaVec(inletSS,  aIn);
        double areaOut = faceAreaVec(outletSS, aOut);

        // Inlet velocity along the INWARD normal (-outward), magnitude Uinf.
        // sgn flips the whole convention if the mesh's face winding makes the
        // default come out outward (--inlet-flip-normal).
        const double sgn = inletFlipNormal ? +1.0 : -1.0;
        if (areaIn > 1e-30)
        {
            s.inletDirX = RealType(sgn * aIn[0] / areaIn);
            s.inletDirY = RealType(sgn * aIn[1] / areaIn);
            s.inletDirZ = RealType(sgn * aIn[2] / areaIn);
        }

        // -------- Per-node inlet inward normals --------
        // The single global inletDir makes every inlet vector parallel, which on
        // a curved inlet looks non-normal. Instead give each inlet node its OWN
        // normal: the area-weighted sum of the inlet faces touching it. The
        // accumulation is keyed by runtime local node id (same id space as
        // s.inletNodes, which sideSetNodes() returns) and folded across ranks
        // with the standard node-halo primitives, so a node whose incident faces
        // span ranks still gets the full sum before normalizing.
        if (inletPernodeNormal && areaIn > 1e-30)
        {
            const size_t nNodes = amr.domain().getNodeCount();
            std::vector<RealType> h_nx(nNodes, RealType(0)),
                                  h_ny(nNodes, RealType(0)),
                                  h_nz(nNodes, RealType(0));
            auto tit = ss.triangleCoordsByName.find(inletSS);
            if (tit != ss.triangleCoordsByName.end())
            {
                // Resolve every triangle vertex to its local id (KeepMisses keeps
                // 1:1 alignment with the coord array). Scatter each face exactly
                // ONCE globally -- on the rank that OWNS its first node -- the
                // same once-per-face gate faceAreaVec uses. reverseExchangeNodeHaloAdd
                // then routes any contribution that landed on a ghost endpoint to
                // that node's true owner; the forward exchange publishes the
                // complete owner sum back to all ghosts. This mirrors the lumped-
                // mass scatter pattern (owned-only scatter + reverse-add + forward).
                std::vector<int> triLocal =
                    amr.domain().resolveSideSetNodesToLocalKeepMisses(tit->second);
                for (size_t f = 0; f + 2 < tit->second.size(); f += 3)
                {
                    int la = triLocal[f];
                    if (la < 0 || hostOwnFA[la] != 1) continue;  // count once: owner of first node
                    const auto& A = tit->second[f];
                    const auto& B = tit->second[f + 1];
                    const auto& C = tit->second[f + 2];
                    double e1x = B[0]-A[0], e1y = B[1]-A[1], e1z = B[2]-A[2];
                    double e2x = C[0]-A[0], e2y = C[1]-A[1], e2z = C[2]-A[2];
                    double ax = 0.5*(e1y*e2z - e1z*e2y);
                    double ay = 0.5*(e1z*e2x - e1x*e2z);
                    double az = 0.5*(e1x*e2y - e1y*e2x);
                    for (int j = 0; j < 3; ++j)
                    {
                        int li = triLocal[f + j];
                        if (li < 0 || (size_t)li >= nNodes) continue;
                        h_nx[li] += RealType(ax);
                        h_ny[li] += RealType(ay);
                        h_nz[li] += RealType(az);
                    }
                }
            }

            // Fold cross-rank partials onto owners, then publish complete sums
            // back to ghosts (mirrors every other nodal-accumulate path here).
            cstone::DeviceVector<RealType> d_nx(nNodes), d_ny(nNodes), d_nz(nNodes);
            cudaMemcpy(d_nx.data(), h_nx.data(), nNodes*sizeof(RealType), cudaMemcpyHostToDevice);
            cudaMemcpy(d_ny.data(), h_ny.data(), nNodes*sizeof(RealType), cudaMemcpyHostToDevice);
            cudaMemcpy(d_nz.data(), h_nz.data(), nNodes*sizeof(RealType), cudaMemcpyHostToDevice);
            amr.domain().reverseExchangeNodeHaloAdd(d_nx);
            amr.domain().reverseExchangeNodeHaloAdd(d_ny);
            amr.domain().reverseExchangeNodeHaloAdd(d_nz);
            amr.domain().exchangeNodeHalo(d_nx);
            amr.domain().exchangeNodeHalo(d_ny);
            amr.domain().exchangeNodeHalo(d_nz);
            cudaMemcpy(h_nx.data(), d_nx.data(), nNodes*sizeof(RealType), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_ny.data(), d_ny.data(), nNodes*sizeof(RealType), cudaMemcpyDeviceToHost);
            cudaMemcpy(h_nz.data(), d_nz.data(), nNodes*sizeof(RealType), cudaMemcpyDeviceToHost);

            // Per inlet node: normalize then flip to INWARD (negate), matching the
            // -aIn/areaIn convention of the global normal. Degenerate sums fall
            // back to the global inlet direction so no node is left unset.
            s.inletDirXPerNode.resize(s.inletNodes.size());
            s.inletDirYPerNode.resize(s.inletNodes.size());
            s.inletDirZPerNode.resize(s.inletNodes.size());
            for (size_t k = 0; k < s.inletNodes.size(); ++k)
            {
                int li = s.inletNodes[k];
                double nx = 0, ny = 0, nz = 0;
                if (li >= 0 && (size_t)li < nNodes)
                {
                    nx = h_nx[li]; ny = h_ny[li]; nz = h_nz[li];
                }
                double len = std::sqrt(nx*nx + ny*ny + nz*nz);
                if (len > 1e-30)
                {
                    s.inletDirXPerNode[k] = RealType(sgn * nx / len);
                    s.inletDirYPerNode[k] = RealType(sgn * ny / len);
                    s.inletDirZPerNode[k] = RealType(sgn * nz / len);
                }
                else
                {
                    s.inletDirXPerNode[k] = s.inletDirX;
                    s.inletDirYPerNode[k] = s.inletDirY;
                    s.inletDirZPerNode[k] = s.inletDirZ;
                }
            }
            if (rank == 0)
                std::cout << "    inlet:  per-node inward normals ON ("
                          << s.inletNodes.size() << " local nodes; --no-inlet-pernode-normal to disable)\n";
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

        // MARS_BC_OWN_DEBUG=1: per-rank confirm of the resolve-vs-own relation
        // for the inlet/outlet side-sets. For each side-set, on THIS rank:
        //   coords   = # side-set node coords the reader gave us (per-rank
        //              complete -> same on every rank)
        //   resolved = # of those coords whose SFC key is in this rank's local
        //              map (held owned-or-ghost here)
        //   ownedRes = # of resolved nodes whose corrected ownership flag == 1
        //              (this rank is the authoritative owner)
        // Global SUM(ownedRes) MUST equal the # of distinct side-set nodes (a
        // global node is owned by exactly one rank). If it is 0 while resolved
        // sum > 0, the owner of every inlet node was demoted to ghost (the
        // Step-2 over-yield orphan) -- every holder sees it only as a ghost.
        if (std::getenv("MARS_BC_OWN_DEBUG"))
        {
            auto ownDbg = [&](const std::string& nm, const std::vector<int>& resolved) {
                auto cit = ss.nodeCoordsByName.find(nm);
                long long lc = (cit == ss.nodeCoordsByName.end()) ? 0
                                                                  : (long long)cit->second.size();
                long long lr = (long long)resolved.size();
                long long lo = 0;
                for (int li : resolved)
                    if (li >= 0 && (size_t)li < hostOwnFA.size() && hostOwnFA[li] == 1) ++lo;
                long long gc = 0, gr = 0, go = 0;
                MPI_Allreduce(&lc, &gc, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
                MPI_Allreduce(&lr, &gr, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(&lo, &go, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
                std::cout << "[bc-own-debug r" << rank << "] " << nm
                          << ": coords=" << lc << " resolved=" << lr
                          << " ownedRes=" << lo << std::endl;
                std::cout.flush();
                MPI_Barrier(MPI_COMM_WORLD);
                if (rank == 0)
                    std::cout << "[bc-own-debug GLOBAL] " << nm
                              << ": found=" << gc << " resolved(sum)=" << gr
                              << " ownedRes(sum)=" << go
                              << "  (ownedRes(sum) must == found; 0 => owner demoted to ghost)"
                              << std::endl;
            };
            ownDbg(inletSS,  s.inletNodes);
            ownDbg(outletSS, s.outletNodes);
        }

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

    // Per-node BC tag so the side sets survive into the VTU (Exodus side sets
    // are not a VTK concept; we carry membership as a point scalar instead).
    // 0=interior, 1=inlet, 2=outlet. In ParaView, Threshold on bc_tag to extract
    // exactly the inlet/outlet nodes -- the by-value equivalent of Extract Block.
    cstone::DeviceVector<RealType> d_bcTag;
    if (vw) {
        std::vector<RealType> h_tag(s.d_u.size(), RealType(0));
        for (int li : s.inletNodes)  if (li >= 0 && size_t(li) < h_tag.size()) h_tag[li] = RealType(1);
        for (int li : s.outletNodes) if (li >= 0 && size_t(li) < h_tag.size()) h_tag[li] = RealType(2);
        d_bcTag.resize(h_tag.size());
        cudaMemcpy(d_bcTag.data(), h_tag.data(), h_tag.size() * sizeof(RealType),
                   cudaMemcpyHostToDevice);
    }

    // side_set_id carries EVERY named Exodus side set (1,2,3,... in sideSetNames
    // order), not just inlet/outlet, so ParaView can Threshold on side_set_id==N
    // to extract any set -- the by-value equivalent of Exodus Extract Block.
    // Sourced from the domain's stored sets (commit 9a83b0d).
    cstone::DeviceVector<RealType> d_ssId;
    std::vector<std::string> ssIdNames;  // index i -> side_set_id (i+1); printed once below
    if (vw) {
        std::vector<RealType> h_id(s.d_u.size(), RealType(0));
        const auto names = amr.domain().sideSetNames();
        for (size_t k = 0; k < names.size(); ++k) {
            const auto& d_loc = amr.domain().sideSetNodes(names[k]);
            std::vector<int> h_loc(d_loc.size());
            if (!h_loc.empty())
                thrust::copy(thrust::device_pointer_cast(d_loc.data()),
                             thrust::device_pointer_cast(d_loc.data() + d_loc.size()),
                             h_loc.begin());
            for (int li : h_loc) if (li >= 0 && size_t(li) < h_id.size()) h_id[li] = RealType(k + 1);
            ssIdNames.push_back(names[k]);
        }
        d_ssId.resize(h_id.size());
        cudaMemcpy(d_ssId.data(), h_id.data(), h_id.size() * sizeof(RealType),
                   cudaMemcpyHostToDevice);
        if (rank == 0 && !ssIdNames.empty()) {
            std::cout << "  side_set_id legend (Threshold on side_set_id==N):\n";
            for (size_t k = 0; k < ssIdNames.size(); ++k)
                std::cout << "    " << (k + 1) << " = [" << ssIdNames[k] << "]\n";
        }
    }

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
        fields.push_back({ "bc_tag", FD::Kind::PointScalar, &d_bcTag, nullptr, nullptr });
        fields.push_back({ "side_set_id", FD::Kind::PointScalar, &d_ssId, nullptr, nullptr });
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

    // Adaptive-dt (--cfl) setup. Explicit BJ+EXT2 advection has a real CFL limit:
    // as the jet develops, uMax*dt/dx crosses the EXT2 stability bound and the
    // run diverges. Cap the advective CFL by shrinking dt. dx is the mean cell
    // size estimated from the bbox diagonal and the GLOBAL element count:
    // dx ~ L_bbox / nElem^(1/3). Uses Lscale (always > 0) not V_total (which is
    // tiny at physical cm-scale and underflows), so dx is robust at any scale.
    unsigned long long nElemLocal = amr.domain().getElementCount();
    unsigned long long nElemGlobal = 0;
    MPI_Allreduce(&nElemLocal, &nElemGlobal, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    double dxMean = (nElemGlobal > 0 && Lscale > 0)
                    ? double(Lscale) / std::cbrt(double(nElemGlobal)) : double(Lscale);
    double dtInit = dt;                         // never grow above the user's dt
    double simTime = 0.0;                        // accumulated time (dt may vary)
    bool   adaptDt = (cflMax > 0 && dxMean > 0);
    if (rank == 0 && adaptDt)
        std::cout << "  adaptive dt: CFL<=" << cflMax << "  dx_mean=" << std::scientific
                  << std::setprecision(3) << dxMean << "  dt0=" << dt
                  << "  cap=" << dtInit << std::defaultfloat << "\n";

    // Enable the per-step solver diagnostics for the first steps -- lets us
    // compare |gP|max (pressure-gradient magnitude) and div between 1 and N
    // ranks. Enabled for ALL rank counts so single vs multi is comparable.
    g_nsDebugStepsLeft = 40;

    auto wallStart = std::chrono::high_resolution_clock::now();
    for (int step = 1; step <= numSteps; ++step)
    {
        // Adapt dt to the current jet strength before stepping. uMax over the
        // whole field; cap dt at cflMax*dx/uMax, never above dtInit. BDF2 uses
        // CONSTANT-dt coefficients ((4/3)u^n-(1/3)u^{n-1}); a dt that changes
        // fast breaks that assumption and injects instability -- so limit BOTH
        // grow and SHRINK to <=5%/step. A genuine CFL drop is reached over a
        // few steps instead of one big jump, keeping BDF2 well-posed.
        if (adaptDt)
        {
            RealType ux = maxOwnedInteriorAbs<KeyType, RealType, TetTag>(s, s.d_u);
            RealType uy = maxOwnedInteriorAbs<KeyType, RealType, TetTag>(s, s.d_v);
            RealType uz = maxOwnedInteriorAbs<KeyType, RealType, TetTag>(s, s.d_w);
            double um = std::sqrt(double(ux)*double(ux) + double(uy)*double(uy) + double(uz)*double(uz));
            double dtCfl = (um > 0) ? cflMax * dxMean / um : dtInit;
            double dtNew = std::min(dtCfl, dtInit);
            dtNew = std::min(dtNew, 1.05 * dt);     // grow <=5%/step (BDF2)
            dtNew = std::max(dtNew, 0.95 * dt);     // shrink <=5%/step (BDF2)
            dtNew = std::max(dtNew, 1e-6 * dtInit); // floor: never collapse dt to ~0
            dt = dtNew;
        }
        runNsStep<KeyType, RealType, TetTag>(s, RealType(dt), RealType(nu), RealType(rho));
        simTime += dt;

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
            double flowThroughs = simTime / tFlow;     // accumulated time (dt may vary)
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
                if (adaptDt)
                    std::cout << "  dt=" << std::scientific << std::setprecision(2) << dt;
                std::cout << "  cg_p=" << s.lastPressureIters
                          << "  pres_r0=" << std::scientific << std::setprecision(3) << s.lastPressR0
                          << "  pres_res=" << std::scientific << std::setprecision(3) << s.lastPressResid
                          << "  cg_uvw=" << s.lastUIters
                          << "\n" << std::defaultfloat;
            }
        }
        if (!vtuPrefix.empty() && (step % vtuEvery == 0 || step == numSteps))
            writeFrame(step, simTime);
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
