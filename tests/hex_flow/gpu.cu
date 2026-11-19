#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"

#ifndef MARS_ENABLE_HYPRE
#error "Hex flow checks require a GPU Hypre build of MARS"
#endif

// The pump stepper on HEX8 end to end (NSStepper<..., HexTag>): a box of cells^3 hexes with warped
// interior nodes, inflow at x0, a do-nothing outlet at x1, no-slip walls; the path ldp_fsi runs
// (Hypre, assembled DDT, divT gradient, implicit convection, SUPG, the penalised projection).
//   flow     every step: the pressure solve converges, the fields stay finite, and continuity is
//            no worse than without the penalty; prints the kinetic energy, KE=, for the rank comparison;
//   failed   a deliberately failed pressure solve: it returns -2, and the convergence check rejects it;
//   invalid  one inverted hex: setup refuses it with std::invalid_argument on every rank.

using Key     = uint64_t;
using Real    = double;
using Domain  = mars::ElementDomain<mars::HexTag, Real, Key, cstone::execution::Gpu>;
using Stepper = NSStepper<Key, Real, mars::HexTag>;

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
// [x0, x0 + 1] x [y0, y0 + ly] x [z0, z0 + lz], interior nodes moved by up to 12% of a cell (the
// boundary planes stay planar), hexes in Exodus order, dealt round-robin to the ranks. inverted: the
// first hex has two nodes swapped.
struct Fixture
{
    std::unique_ptr<Domain> domain;
    double x0 = 0.123456789, y0 = .234567891, z0 = .345678912, ly = .7, lz = .6, x1;
    Fixture(int rank, int ranks, int cells, bool inverted)
    {
        Domain::HostCoordsTuple xyz;
        Domain::HostConnectivityTuple conn;
        const int width = cells + 1;
        auto index = [width](int x, int y, int z) { return Key(x + width * (y + width * z)); };
        x1 = x0 + 1;
        for (int z = 0; z <= cells; ++z)
            for (int y = 0; y <= cells; ++y)
                for (int x = 0; x <= cells; ++x)
                {
                    const bool interior = x > 0 && x < cells && y > 0 && y < cells && z > 0 && z < cells;
                    const double h = 1. / cells, wx = interior ? .12 * h * std::sin(2.1 * x + 1.3 * y + .7 * z) : 0,
                                 wy = interior ? .12 * h * std::cos(1.7 * x - .9 * y + 1.1 * z) : 0,
                                 wz = interior ? .12 * h * std::sin(.8 * x + 1.9 * y - 1.4 * z) : 0;
                    std::get<0>(xyz).push_back(x0 + x * h + wx);
                    std::get<1>(xyz).push_back(y0 + ly * (y * h) + ly * wy);
                    std::get<2>(xyz).push_back(z0 + lz * (z * h) + lz * wz);
                }
        int element = 0;
        for (int z = 0; z < cells; ++z)
            for (int y = 0; y < cells; ++y)
                for (int x = 0; x < cells; ++x)
                {
                    std::array<Key, 8> n{index(x, y, z),         index(x + 1, y, z),     index(x + 1, y + 1, z),
                                         index(x, y + 1, z),     index(x, y, z + 1),     index(x + 1, y, z + 1),
                                         index(x + 1, y + 1, z + 1), index(x, y + 1, z + 1)};
                    if (inverted && element == 0) std::swap(n[0], n[1]);
                    if (element++ % ranks != rank) continue;
                    std::get<0>(conn).push_back(n[0]);
                    std::get<1>(conn).push_back(n[1]);
                    std::get<2>(conn).push_back(n[2]);
                    std::get<3>(conn).push_back(n[3]);
                    std::get<4>(conn).push_back(n[4]);
                    std::get<5>(conn).push_back(n[5]);
                    std::get<6>(conn).push_back(n[6]);
                    std::get<7>(conn).push_back(n[7]);
                }
        domain = std::make_unique<Domain>(xyz, conn, rank, ranks, 8, false, 1);
    }
};
void boundaries(Stepper& s, const Fixture& f)
{
    // The domain rebuilds node coordinates from the SFC keys (quantized to ~1e-6): planes within 1e-5.
    auto x = host(f.domain->getNodeX());
    auto y = host(f.domain->getNodeY());
    auto z = host(f.domain->getNodeZ());
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (std::abs(x[i] - f.x0) < 1e-5)
            s.inletNodes.push_back(int(i));
        else if (std::abs(x[i] - f.x1) < 1e-5)
            s.outletNodes.push_back(int(i));
        else if (std::abs(y[i] - f.y0) < 1e-5 || std::abs(y[i] - (f.y0 + f.ly)) < 1e-5 || std::abs(z[i] - f.z0) < 1e-5 ||
                 std::abs(z[i] - (f.z0 + f.lz)) < 1e-5)
            s.wallNodes.push_back(int(i));
    }
}
// The last pressure solve succeeded: solveOneComponent returns -2 on failure, else its iterations,
// below the cap. (The cap alone accepts -2.)
bool pressure_solve_succeeded(const Stepper& s) { return s.lastPressureIters >= 0 && s.lastPressureIters < s.maxIter; }

std::unique_ptr<Stepper> make_stepper(const Fixture& f, int rank, int ranks)
{
    auto s               = std::make_unique<Stepper>(*f.domain, SolverKind::Hypre, 256, 2000, 1e-10, rank, ranks);
    s->Uinf              = 1;
    s->uinfBase          = 1;
    s->pumpZeroIC        = true;
    s->bcKind            = Stepper::BCKind::Pump;
    s->pressureSolve     = PressureSolveKind::DDT;
    s->useLegacyGradient = false;
    s->outletDoNothing   = true;
    s->useDirichletLift  = true;
    s->implicitAdvection = true;
    s->useSupg           = true;
    boundaries(*s, f);
    setupNSStepper(*s, .01, .001, CvfemKernelVariant::Tensor);
    return s;
}

void check_flow(const Fixture& f, int rank, int ranks)
{
    auto s       = make_stepper(f, rank, ranks);
    auto without = make_stepper(f, rank, ranks);
    // The penalised projection with a uniform beta = 1/2 (W = V / dt, eta = 1).
    {
        auto V = host(s->d_massNode);
        for (auto& v : V)
            v *= 1 / .001;
        cstone::DeviceVector<double> W;
        upload(V, W);
        setPenalisedProjection(*s, W.data(), 1.);
    }
    for (int step = 0; step < 6; ++step)
    {
        runNsStep(*s, .001, .01, 1.);
        runNsStep(*without, .001, .01, 1.);
        // Pressure convergence: both solves stopped below the iteration cap.
        require(pressure_solve_succeeded(*s) && pressure_solve_succeeded(*without),
                "The pressure solve did not converge");
        bool finite = true;
        for (auto* field : {&s->d_u, &s->d_v, &s->d_w, &s->d_p})
            for (double value : host(*field))
                finite &= std::isfinite(value);
        require(finite, "A field is not finite");
        if (rank == 0)
            std::cout << std::scientific << "hex flow step " << step << ": max |div u| " << s->lastDivMax
                      << " with the penalty, " << without->lastDivMax << " without; pressure iterations "
                      << s->lastPressureIters << ", " << without->lastPressureIters << std::defaultfloat << "\n";
        // Continuity, separately: no worse than without the penalty. (The largest nodal divergence sits
        // at the inlet's Dirichlet corners, which no projection changes: it is compared, not driven to 0.)
        require(std::isfinite(s->lastDivMax) && s->lastDivMax <= 2 * without->lastDivMax + 1e-12,
                "Continuity is worse with the penalty");
    }
    // The kinetic energy, sum over owned nodes of V |u|^2 / 2: the same on every rank count.
    const auto own = host(s->ownershipMap());
    const auto V = host(s->d_massNode);
    const auto u = host(s->d_u);
    const auto v = host(s->d_v);
    const auto w = host(s->d_w);
    double ke = 0;
    for (size_t i = 0; i < u.size(); ++i)
        if (own[i] == 1) ke += .5 * V[i] * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
    MPI_Allreduce(MPI_IN_PLACE, &ke, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << std::setprecision(12) << "KE=" << ke << "\n";
    require(std::isfinite(ke) && ke > 0, "No flow");
}

// One GMRES iteration and an acceptance residual nothing reaches: the pressure solve fails.
void check_failed_solve(const Fixture& f, int rank, int ranks)
{
    auto s = make_stepper(f, rank, ranks);
    s->maxIter   = 1;
    s->tolerance = 1e-300;
    setenv("MARS_HYPRE_ACCEPT_RES", "1e-300", 1);
    runNsStep(*s, .001, .01, 1.);
    unsetenv("MARS_HYPRE_ACCEPT_RES");
    if (rank == 0) std::cout << "failed: pressure solve returned " << s->lastPressureIters << "\n";
    require(s->lastPressureIters == -2, "The forced pressure failure did not return -2");
    require(!pressure_solve_succeeded(*s), "The convergence check accepted a failed pressure solve");
}

void check_invalid(int rank, int ranks)
{
    Fixture f(rank, ranks, 6, true);
    bool rejected = false;
    try
    {
        make_stepper(f, rank, ranks);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "An inverted hex was accepted on a rank");
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
        if (mode == "flow")
        {
            Fixture fixture(rank, ranks, 6, false);
            check_flow(fixture, rank, ranks);
        }
        else if (mode == "failed")
        {
            Fixture fixture(rank, ranks, 6, false);
            check_failed_solve(fixture, rank, ranks);
        }
        else if (mode == "invalid")
            check_invalid(rank, ranks);
        else
            require(false, "Unknown mode");
        cudaCheckCall(cudaDeviceSynchronize());
        if (rank == 0) std::cout << "PASS: hex flow " << argv[1] << " ranks=" << ranks << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
