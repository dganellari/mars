#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"

#ifndef MARS_ENABLE_HYPRE
#error "Face flux checks require a GPU Hypre build of MARS"
#endif

// The face-flux projection (useFaceFlux) on public Tet4 and Hex8 boxes, inflow at x0 (Uinf along +x),
// a do-nothing pressure outlet at x1, no-slip walls; Hypre, divT gradient. The first pressure solve
// runs the operator self-check, which covers the outlet rows (the corrector's own gradient on the
// outlets' boundary faces) and aborts on a mismatch.
//   {tet,hex}_{ddt,k}_{scs,source}
//            assembled DDT or K; the openings' flux through the SCS term (useOpeningFluxSource =
//            false: where the outlet flux used to enter twice on DDT, and the inlet's not at all on K)
//            or the prescribed source (true, as ldp_fsi runs). Every step: the pressure solve
//            succeeds; each control volume balances (lastFaceFluxDivMax, the cached balance); and,
//            independently of it, the inflow and the outflow measured on the face flux itself (the
//            net interior flux out of the inlet plane's control volumes, and into the outlet
//            plane's) and on the solved velocity at the eligible outlets (u . A_b) all equal Uinf
//            times the inlet's area;
//   failed   a deliberately failed pressure solve returns -2, and the convergence check rejects it;
//   reject   K with the legacy gradient, CG, MARS_HYPRE_USE_DDT, and no outlet without a pin:
//            std::invalid_argument on every rank, nothing changed.

using Key  = uint64_t;
using Real = double;

void require(bool pass, const char* message)
{
    int bad = pass ? 0 : 1, global = 0;
    MPI_Allreduce(&bad, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global) throw std::runtime_error(message);
}
template<class T>
std::vector<T> host(const cstone::DeviceVector<T>& d)
{
    std::vector<T> h(d.size());
    if (!h.empty()) cudaCheckCall(cudaMemcpy(h.data(), d.data(), h.size() * sizeof(T), cudaMemcpyDeviceToHost));
    return h;
}
template<class T>
void upload(const std::vector<T>& h, cstone::DeviceVector<T>& d)
{
    d.resize(h.size());
    if (!h.empty()) cudaCheckCall(cudaMemcpy(d.data(), h.data(), h.size() * sizeof(T), cudaMemcpyHostToDevice));
}

// The box [x0, x0 + 1] x [y0, y0 + ly] x [z0, z0 + lz] of cells^3 cubes: five tets each (the
// tests/outlet_pressure pattern) or one hex with warped interior nodes; elements dealt round-robin.
struct Box
{
    static constexpr double x0 = 0.123456789, y0 = .234567891, z0 = .345678912, ly = .7, lz = .6;
    int cells;
    explicit Box(int cells_) : cells(cells_) {}
    std::array<double, 3> point(int x, int y, int z, bool warp) const
    {
        const bool interior = warp && x > 0 && x < cells && y > 0 && y < cells && z > 0 && z < cells;
        const double h      = 1. / cells;
        const double wx = interior ? .12 * h * std::sin(2.1 * x + 1.3 * y + .7 * z) : 0,
                     wy = interior ? .12 * h * std::cos(1.7 * x - .9 * y + 1.1 * z) : 0,
                     wz = interior ? .12 * h * std::sin(.8 * x + 1.9 * y - 1.4 * z) : 0;
        return {x0 + x * h + wx, y0 + ly * (y * h + wy), z0 + lz * (z * h + wz)};
    }
    Key index(int x, int y, int z) const { return Key(x + (cells + 1) * (y + (cells + 1) * z)); }
};
double determinant(std::array<double, 3> a, std::array<double, 3> b, std::array<double, 3> c, std::array<double, 3> d)
{
    for (int k = 0; k < 3; ++k)
    {
        b[k] -= a[k];
        c[k] -= a[k];
        d[k] -= a[k];
    }
    return b[0] * (c[1] * d[2] - c[2] * d[1]) - b[1] * (c[0] * d[2] - c[2] * d[0]) + b[2] * (c[0] * d[1] - c[1] * d[0]);
}
template<class Tag>
std::unique_ptr<mars::ElementDomain<Tag, Real, Key, cstone::execution::Gpu>> make_domain(const Box& box, int rank, int ranks)
{
    using Domain   = mars::ElementDomain<Tag, Real, Key, cstone::execution::Gpu>;
    constexpr bool hex = std::is_same_v<Tag, mars::HexTag>;
    typename Domain::HostCoordsTuple xyz;
    typename Domain::HostConnectivityTuple conn;
    const int n = box.cells;
    std::vector<std::array<double, 3>> points;
    for (int z = 0; z <= n; ++z)
        for (int y = 0; y <= n; ++y)
            for (int x = 0; x <= n; ++x)
            {
                points.push_back(box.point(x, y, z, hex));
                std::get<0>(xyz).push_back(points.back()[0]);
                std::get<1>(xyz).push_back(points.back()[1]);
                std::get<2>(xyz).push_back(points.back()[2]);
            }
    int element = 0;
    for (int z = 0; z < n; ++z)
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x)
            {
                if constexpr (hex)
                {
                    const Key c[8] = {box.index(x, y, z),         box.index(x + 1, y, z),
                                      box.index(x + 1, y + 1, z), box.index(x, y + 1, z),
                                      box.index(x, y, z + 1),     box.index(x + 1, y, z + 1),
                                      box.index(x + 1, y + 1, z + 1), box.index(x, y + 1, z + 1)};
                    if (element++ % ranks != rank) continue;
                    std::get<0>(conn).push_back(c[0]);
                    std::get<1>(conn).push_back(c[1]);
                    std::get<2>(conn).push_back(c[2]);
                    std::get<3>(conn).push_back(c[3]);
                    std::get<4>(conn).push_back(c[4]);
                    std::get<5>(conn).push_back(c[5]);
                    std::get<6>(conn).push_back(c[6]);
                    std::get<7>(conn).push_back(c[7]);
                }
                else
                {
                    const int pattern[5][4] = {{1, 2, 4, 7}, {0, 1, 2, 4}, {3, 1, 2, 7}, {5, 1, 4, 7}, {6, 2, 4, 7}};
                    for (auto& cell : pattern)
                    {
                        std::array<Key, 4> nodes;
                        for (int a = 0; a < 4; ++a)
                        {
                            const int c = cell[a] ^ ((x + y + z) & 1);
                            nodes[a]    = box.index(x + (c & 1), y + ((c >> 1) & 1), z + ((c >> 2) & 1));
                        }
                        if (determinant(points[nodes[0]], points[nodes[1]], points[nodes[2]], points[nodes[3]]) < 0)
                            std::swap(nodes[0], nodes[1]);
                        if (element++ % ranks != rank) continue;
                        std::get<0>(conn).push_back(nodes[0]);
                        std::get<1>(conn).push_back(nodes[1]);
                        std::get<2>(conn).push_back(nodes[2]);
                        std::get<3>(conn).push_back(nodes[3]);
                    }
                }
            }
    return std::make_unique<Domain>(xyz, conn, rank, ranks, 8, false, 1);
}

// Boundaries, and the openings' outward area per node (owner-complete): the rectangular dual of the
// face grid, a quarter cell per adjacent cell, so each opening's areas sum to its area.
template<class Stepper, class Domain>
void openings(Stepper& s, const Domain& domain, const Box& box, bool outlet_area)
{
    // The domain rebuilds node coordinates from the SFC keys (quantized to ~1e-6): planes within 1e-5.
    const auto x = host(domain.getNodeX());
    const auto y = host(domain.getNodeY());
    const auto z = host(domain.getNodeZ());
    const double hy = Box::ly / box.cells, hz = Box::lz / box.cells, x1 = Box::x0 + 1;
    std::vector<double> ax(x.size(), 0), zero(x.size(), 0), bx(x.size(), 0);
    for (size_t i = 0; i < x.size(); ++i)
    {
        const int j = int(std::lround((y[i] - Box::y0) / hy)), k = int(std::lround((z[i] - Box::z0) / hz));
        const double w = (j == 0 || j == box.cells ? .5 : 1.) * (k == 0 || k == box.cells ? .5 : 1.) * hy * hz;
        if (std::abs(x[i] - Box::x0) < 1e-5)
        {
            s.inletNodes.push_back(int(i));
            ax[i] = -w;
        }
        else if (std::abs(x[i] - x1) < 1e-5)
        {
            s.outletNodes.push_back(int(i));
            bx[i] = outlet_area ? w : 0;
        }
        else if (std::abs(y[i] - Box::y0) < 1e-5 || std::abs(y[i] - (Box::y0 + Box::ly)) < 1e-5 ||
                 std::abs(z[i] - Box::z0) < 1e-5 || std::abs(z[i] - (Box::z0 + Box::lz)) < 1e-5)
            s.wallNodes.push_back(int(i));
    }
    upload(ax, s.d_inletAreaVecX);
    upload(zero, s.d_inletAreaVecY);
    upload(zero, s.d_inletAreaVecZ);
    upload(bx, s.d_outletAreaVecX);
    upload(zero, s.d_outletAreaVecY);
    upload(zero, s.d_outletAreaVecZ);
}

// The last pressure solve succeeded: solveOneComponent returns -2 on failure, else its iterations,
// below the cap. (The cap alone accepts -2.)
template<class Stepper>
bool pressure_solve_succeeded(const Stepper& s)
{
    return s.lastPressureIters >= 0 && s.lastPressureIters < s.maxIter;
}

template<class Tag>
struct Case
{
    using Domain  = mars::ElementDomain<Tag, Real, Key, cstone::execution::Gpu>;
    using Stepper = NSStepper<Key, Real, Tag>;
    Box box{6};
    std::unique_ptr<Domain> domain;
    std::unique_ptr<Stepper> s;
    Case(int rank, int ranks, PressureSolveKind solve, bool source, bool outlet_area = true)
    {
        domain               = make_domain<Tag>(box, rank, ranks);
        s                    = std::make_unique<Stepper>(*domain, SolverKind::Hypre, 256, 2000, 1e-10, rank, ranks);
        s->Uinf              = 1;
        s->uinfBase          = 1;
        s->pumpZeroIC        = true;
        s->bcKind            = Stepper::BCKind::Pump;
        s->pressureSolve     = solve;
        s->useLegacyGradient = false;
        s->outletDoNothing   = true;
        s->useDirichletLift  = true;
        s->useOpeningFluxSource = source;
        s->useFaceFlux       = true;
        openings(*s, *domain, box, outlet_area);
        setupNSStepper(*s, .01, .001, CvfemKernelVariant::Tensor);
    }
    // The openings' fluxes over all ranks: the inlet (u . A_in at the owned inlet nodes) and the
    // eligible outlets (u . A_b).
    std::array<double, 2> opening_fluxes() const
    {
        const auto own = host(s->ownershipMap());
        const auto out = host(s->d_faceFluxOutlet);
        const auto u = host(s->d_u);
        const auto v = host(s->d_v);
        const auto w = host(s->d_w);
        const auto ix = host(s->d_inletAreaVecX);
        const auto ox = host(s->d_outletAreaVecX);
        std::array<double, 2> q{0, 0};
        for (size_t i = 0; i < u.size(); ++i)
        {
            if (own[i] != 1) continue;
            q[0] += u[i] * ix[i];
            if (i < out.size() && out[i]) q[1] += u[i] * ox[i];
        }
        MPI_Allreduce(MPI_IN_PLACE, q.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return q;
    }
    // On the face flux itself, independently of the openings' terms and the cached balance: the net
    // interior flux out of the inlet plane's control volumes (their shared faces cancel), and into the
    // outlet plane's. faceFluxDivergence is owner-complete.
    std::array<double, 2> plane_fluxes() const
    {
        cstone::DeviceVector<double> d;
        faceFluxDivergence<Key, Real, Tag>(*s, d);
        const auto div = host(d);
        const auto own = host(s->ownershipMap());
        const auto x   = host(domain->getNodeX());
        std::array<double, 2> q{0, 0};
        for (size_t i = 0; i < div.size(); ++i)
        {
            if (own[i] != 1) continue;
            if (std::abs(x[i] - Box::x0) < 1e-5) q[0] += div[i];
            if (std::abs(x[i] - (Box::x0 + 1)) < 1e-5) q[1] -= div[i];
        }
        MPI_Allreduce(MPI_IN_PLACE, q.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return q;
    }
};

template<class Tag>
void check_conservation(int rank, int ranks, PressureSolveKind solve, bool source)
{
    Case<Tag> c(rank, ranks, solve, source);
    const double inflow = 1. * Box::ly * Box::lz;  // Uinf times the inlet's area
    for (int step = 0; step < 4; ++step)
    {
        runNsStep(*c.s, .001, .01, 1.);
        require(pressure_solve_succeeded(*c.s), "The pressure solve did not converge");
        const auto q     = c.opening_fluxes();
        const auto plane = c.plane_fluxes();
        if (rank == 0)
            std::cout << std::scientific << "step " << step << ": max |CV balance| " << c.s->lastFaceFluxDivMax
                      << ", sum " << c.s->lastFaceFluxDivSum << "; inflow " << -q[0] << " (prescribed), "
                      << plane[0] << " (face flux out of the inlet plane); outflow " << q[1] << " (u . A_b), "
                      << plane[1] << " (face flux into the outlet plane); pressure iterations "
                      << c.s->lastPressureIters << std::defaultfloat << "\n";
        // Every control volume, openings included, balances to the solve's tolerance.
        require(c.s->lastFaceFluxDivMax <= 1e-7 * inflow, "A control volume does not balance");
        // The inlet carries Uinf times its area in, and the outlets carry it out.
        require(std::abs(q[0] + inflow) <= 1e-12 * inflow, "The inlet flux is not Uinf times its area");
        require(std::abs(q[0] + q[1]) <= 1e-7 * inflow, "The openings' fluxes do not sum to zero");
        // Independently of the openings' terms: the face flux carries the inflow out of the inlet plane
        // and into the outlet plane.
        require(std::abs(plane[0] - inflow) <= 1e-7 * inflow, "The face flux does not carry the inflow off the inlet");
        require(std::abs(plane[1] - inflow) <= 1e-7 * inflow, "The face flux does not carry the inflow to the outlet");
    }
}

void rejected_without_mutation(NSStepper<Key, Real, mars::TetTag>& s)
{
    const auto p = host(s.d_p);
    const auto u = host(s.d_u);
    const auto v = host(s.d_v);
    const auto w = host(s.d_w);
    const int history = s.bdfStep;
    bool rejected     = false;
    try
    {
        runNsStep(s, .001, .01, 1.);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "An unsupported face-flux configuration was accepted on a rank");
    require(p == host(s.d_p) && u == host(s.d_u) && v == host(s.d_v) && w == host(s.d_w) && s.bdfStep == history,
            "A rejected step changed the state");
}
// One GMRES iteration and an acceptance residual nothing reaches: the pressure solve fails.
void check_failed_solve(int rank, int ranks)
{
    Case<mars::TetTag> c(rank, ranks, PressureSolveKind::DDT, true);
    c.s->maxIter   = 1;
    c.s->tolerance = 1e-300;
    setenv("MARS_HYPRE_ACCEPT_RES", "1e-300", 1);
    runNsStep(*c.s, .001, .01, 1.);
    unsetenv("MARS_HYPRE_ACCEPT_RES");
    if (rank == 0) std::cout << "failed: pressure solve returned " << c.s->lastPressureIters << "\n";
    require(c.s->lastPressureIters == -2, "The forced pressure failure did not return -2");
    require(!pressure_solve_succeeded(*c.s), "The convergence check accepted a failed pressure solve");
}

void check_rejections(int rank, int ranks)
{
    const bool inject = rank == ranks - 1;  // on one rank only: the refusal must be collective
    {
        Case<mars::TetTag> c(rank, ranks, PressureSolveKind::DDT, true);
        auto& s = *c.s;
        if (inject)
        {
            s.pressureSolve     = PressureSolveKind::K;
            s.useLegacyGradient = true;
        }
        rejected_without_mutation(s);
        s.pressureSolve     = PressureSolveKind::DDT;
        s.useLegacyGradient = false;
        if (inject) s.solverKind = SolverKind::CG;
        rejected_without_mutation(s);
        s.solverKind = SolverKind::Hypre;
        if (inject) setenv("MARS_HYPRE_USE_DDT", "1", 1);
        rejected_without_mutation(s);
        unsetenv("MARS_HYPRE_USE_DDT");
        // The restored configuration steps.
        runNsStep(s, .001, .01, 1.);
        require(pressure_solve_succeeded(s), "The pressure solve did not converge");
    }
    // No eligible outlet (zero outlet areas) and no pin: no pressure gauge.
    Case<mars::TetTag> c(rank, ranks, PressureSolveKind::DDT, true, false);
    rejected_without_mutation(*c.s);
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank, ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    try
    {
        MPI_Comm local;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &local);
        int local_rank, devices;
        MPI_Comm_rank(local, &local_rank);
        MPI_Comm_free(&local);
        cudaCheckCall(cudaGetDeviceCount(&devices));
        require(devices > 0 && argc == 2, "Expected a GPU and a check mode");
        cudaCheckCall(cudaSetDevice(local_rank % devices));
        mars::requestFullNodeHaloExchange();
        unsetenv("MARS_HYPRE_USE_DDT");
        const std::string mode(argv[1]);
        const auto DDT = PressureSolveKind::DDT, K = PressureSolveKind::K;
        if (mode == "tet_ddt_scs")
            check_conservation<mars::TetTag>(rank, ranks, DDT, false);
        else if (mode == "tet_ddt_source")
            check_conservation<mars::TetTag>(rank, ranks, DDT, true);
        else if (mode == "tet_k_scs")
            check_conservation<mars::TetTag>(rank, ranks, K, false);
        else if (mode == "tet_k_source")
            check_conservation<mars::TetTag>(rank, ranks, K, true);
        else if (mode == "hex_ddt_scs")
            check_conservation<mars::HexTag>(rank, ranks, DDT, false);
        else if (mode == "hex_ddt_source")
            check_conservation<mars::HexTag>(rank, ranks, DDT, true);
        else if (mode == "hex_k_scs")
            check_conservation<mars::HexTag>(rank, ranks, K, false);
        else if (mode == "hex_k_source")
            check_conservation<mars::HexTag>(rank, ranks, K, true);
        else if (mode == "failed")
            check_failed_solve(rank, ranks);
        else if (mode == "reject")
            check_rejections(rank, ranks);
        else
            require(false, "Unknown mode");
        cudaCheckCall(cudaDeviceSynchronize());
        if (rank == 0) std::cout << "PASS: face flux " << argv[1] << " ranks=" << ranks << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
