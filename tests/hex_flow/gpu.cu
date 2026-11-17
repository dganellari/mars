#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"

#ifndef MARS_ENABLE_HYPRE
#error "Hex flow checks require a GPU Hypre build of MARS"
#endif

// The pump stepper on HEX8 end to end (NSStepper<..., HexTag>): a box of cells^3 hexes with warped
// interior nodes, inflow at x0, a do-nothing outlet at x1, no-slip walls; the path ldp_fsi runs
// (Hypre, assembled DDT, divT gradient, implicit convection, SUPG, the penalised projection).
//   flow     every step: the pressure solve converges, the fields stay finite, and the projection
//            removes most of the divergence; prints the kinetic energy, KE=, for the rank comparison;
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
    auto x = host(f.domain->getNodeX()), y = host(f.domain->getNodeY()), z = host(f.domain->getNodeZ());
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (std::abs(x[i] - f.x0) < 1e-9)
            s.inletNodes.push_back(int(i));
        else if (std::abs(x[i] - f.x1) < 1e-9)
            s.outletNodes.push_back(int(i));
        else if (std::abs(y[i] - f.y0) < 1e-9 || std::abs(y[i] - (f.y0 + f.ly)) < 1e-9 || std::abs(z[i] - f.z0) < 1e-9 ||
                 std::abs(z[i] - (f.z0 + f.lz)) < 1e-9)
            s.wallNodes.push_back(int(i));
    }
}
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
    auto s = make_stepper(f, rank, ranks);
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
        require(s->lastPressureIters < s->maxIter, "The pressure solve did not converge");
        bool finite = true;
        for (auto* field : {&s->d_u, &s->d_v, &s->d_w, &s->d_p})
            for (double value : host(*field))
                finite &= std::isfinite(value);
        require(finite, "A field is not finite");
        if (step > 0)
            require(s->lastDivMax <= .1 * s->lastDivMaxPre,
                    "The projection did not remove most of the divergence");
    }
    // The kinetic energy, sum over owned nodes of V |u|^2 / 2: the same on every rank count.
    const auto own = host(s->ownershipMap()), V = host(s->d_massNode);
    const auto u = host(s->d_u), v = host(s->d_v), w = host(s->d_w);
    double ke = 0;
    for (size_t i = 0; i < u.size(); ++i)
        if (own[i] == 1) ke += .5 * V[i] * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
    MPI_Allreduce(MPI_IN_PLACE, &ke, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << std::setprecision(12) << "hex flow: |div u| " << s->lastDivMax << " (before the projection "
                  << s->lastDivMaxPre << "), pressure iterations " << s->lastPressureIters << "\nKE=" << ke << "\n";
    require(std::isfinite(ke) && ke > 0, "No flow");
}

void check_invalid(int rank, int ranks)
{
    Fixture f(rank, ranks, 4, true);
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
            Fixture fixture(rank, ranks, 4, false);
            check_flow(fixture, rank, ranks);
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
