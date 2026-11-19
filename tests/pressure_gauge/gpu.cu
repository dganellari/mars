#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"

#ifndef MARS_ENABLE_HYPRE
#error "Pressure gauge checks require a GPU Hypre build of MARS"
#endif

// The pressure increment's mean (NSStepper::pressureIncrementMean) on the public Tet4 box with
// inflow, a prescribed outlet pressure T = 10 (outletPressureLift: Dirichlet rows) and the
// opening-flux source on: the configuration whose phi mean the pump removed (af88fd3d), so the stored
// pressure drifted off its prescribed value. Hypre, assembled DDT, 5 steps.
//   auto      the default: the outlet keeps p = T (to round-off), the pressure being gauged;
//   always    the old behaviour: the outlet's pressure moves off T (the shift the mean removal makes);
//   velocity  Auto and Always give the same velocity (grad phi does not see the shift), to the
//             pressure solve's tolerance; their pressures differ by a constant.

using Key     = uint64_t;
using Real    = double;
using Domain  = mars::ElementDomain<mars::TetTag, Real, Key, cstone::execution::Gpu>;
using Stepper = NSStepper<Key, Real, mars::TetTag>;
using Point   = std::array<double, 3>;

void require(bool pass, const char* message)
{
    int bad = pass ? 0 : 1, global = 0;
    MPI_Allreduce(&bad, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global) throw std::runtime_error(message);
}
template<class T>
std::vector<T> host(const T* d, size_t n)
{
    std::vector<T> h(n);
    if (n) cudaCheckCall(cudaMemcpy(h.data(), d, n * sizeof(T), cudaMemcpyDeviceToHost));
    return h;
}
template<class T>
std::vector<T> host(const cstone::DeviceVector<T>& d)
{
    return host(d.data(), d.size());
}
template<class T>
void upload(const std::vector<T>& h, cstone::DeviceVector<T>& d)
{
    d.resize(h.size());
    if (!h.empty()) cudaCheckCall(cudaMemcpy(d.data(), h.data(), h.size() * sizeof(T), cudaMemcpyHostToDevice));
}
double determinant(Point a, Point b, Point c, Point d)
{
    for (int k = 0; k < 3; ++k)
    {
        b[k] -= a[k];
        c[k] -= a[k];
        d[k] -= a[k];
    }
    return b[0] * (c[1] * d[2] - c[2] * d[1]) - b[1] * (c[0] * d[2] - c[2] * d[0]) + b[2] * (c[0] * d[1] - c[1] * d[0]);
}
// The box of tests/outlet_pressure: cells^3 cubes of five tets, x in [x0, x0 + 1] (planes of nodes
// 1 / cells apart), elements dealt round-robin to the ranks.
struct Fixture
{
    std::unique_ptr<Domain> domain;
    int cells;
    double x0 = 0.123456789, x1 = 1.253456789;
    Fixture(int rank, int ranks, int cells_ = 4) : cells(cells_)
    {
        Domain::HostCoordsTuple xyz;
        Domain::HostConnectivityTuple conn;
        std::vector<Point> points;
        const int width = cells + 1;
        auto index = [width](int x, int y, int z) { return Key(x + width * (y + width * z)); };
        for (int z = 0; z <= cells; ++z)
            for (int y = 0; y <= cells; ++y)
                for (int x = 0; x <= cells; ++x)
                {
                    Point p{x0 + double(x) / cells, .234567891 + .257 * y * (4. / cells),
                            .345678912 + .243 * z * (4. / cells)};
                    points.push_back(p);
                    std::get<0>(xyz).push_back(p[0]);
                    std::get<1>(xyz).push_back(p[1]);
                    std::get<2>(xyz).push_back(p[2]);
                }
        x1                      = points[cells][0];
        const int pattern[5][4] = {{1, 2, 4, 7}, {0, 1, 2, 4}, {3, 1, 2, 7}, {5, 1, 4, 7}, {6, 2, 4, 7}};
        int element             = 0;
        for (int z = 0; z < cells; ++z)
            for (int y = 0; y < cells; ++y)
                for (int x = 0; x < cells; ++x)
                    for (auto& cell : pattern)
                    {
                        std::array<Key, 4> nodes;
                        for (int a = 0; a < 4; ++a)
                        {
                            int c    = cell[a] ^ ((x + y + z) & 1);
                            nodes[a] = index(x + (c & 1), y + ((c >> 1) & 1), z + ((c >> 2) & 1));
                        }
                        double det =
                            determinant(points[nodes[0]], points[nodes[1]], points[nodes[2]], points[nodes[3]]);
                        if (det < 0) std::swap(nodes[0], nodes[1]);
                        if (element++ % ranks != rank) continue;
                        std::get<0>(conn).push_back(nodes[0]);
                        std::get<1>(conn).push_back(nodes[1]);
                        std::get<2>(conn).push_back(nodes[2]);
                        std::get<3>(conn).push_back(nodes[3]);
                    }
        domain = std::make_unique<Domain>(xyz, conn, rank, ranks, 8, false, 1);
    }
    // The plane of nodes of x (0 at the inlet, cells at the outlet).
    int plane(double x) const { return int(std::lround((x - x0) * cells)); }
};
void boundaries(Stepper& s, const Fixture& f)
{
    auto x = host(f.domain->getNodeX());
    auto y = host(f.domain->getNodeY());
    auto z = host(f.domain->getNodeZ());
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (std::abs(x[i] - f.x0) < 1e-5)
            s.inletNodes.push_back(int(i));
        else if (std::abs(x[i] - f.x1) < 1e-5)
            s.outletNodes.push_back(int(i));
        else if (std::abs(y[i] - .234567891) < 1e-5 || std::abs(y[i] - (.234567891 + 4 * .257)) < 1e-5 ||
                 std::abs(z[i] - .345678912) < 1e-5 || std::abs(z[i] - (.345678912 + 4 * .243)) < 1e-5)
            s.wallNodes.push_back(int(i));
    }
}

constexpr double TARGET = 10;

std::unique_ptr<Stepper> make_stepper(const Fixture& f, int rank, int ranks, Stepper::PressureIncrementMean mode)
{
    auto s                   = std::make_unique<Stepper>(*f.domain, SolverKind::Hypre, 256, 2000, 1e-10, rank, ranks);
    s->Uinf                  = 1;
    s->uinfBase              = 1;
    s->pumpZeroIC            = true;
    s->bcKind                = Stepper::BCKind::Pump;
    s->pressureSolve         = PressureSolveKind::DDT;
    s->useLegacyGradient     = false;
    s->outletDoNothing       = true;
    s->useDirichletLift      = true;
    s->useOpeningFluxSource  = true;
    s->pressureIncrementMean = mode;
    boundaries(*s, f);
    setupNSStepper(*s, .1, .001, CvfemKernelVariant::Tensor);
    // The outlet's prescribed pressure on its Dirichlet rows.
    s->outletPressureLift = true;
    const auto mask       = host(s->d_isPressureBdryDof);
    std::vector<double> targets(s->numOwnedDofs, NAN);
    for (size_t i = 0; i < targets.size(); ++i)
        if (mask[i]) targets[i] = TARGET;
    upload(targets, s->d_pPhiTargetDof);
    return s;
}

// The owned outlet nodes' pressures (the Dirichlet rows): min and max over all ranks.
std::array<double, 2> outlet_pressure(const Stepper& s)
{
    const auto own  = host(s.ownershipMap());
    const auto n2d  = host(s.d_node_to_dof);
    const auto mask = host(s.d_isPressureBdryDof);
    const auto p    = host(s.d_p);
    double lo = 1e300, hi = -1e300;
    for (size_t i = 0; i < p.size(); ++i)
    {
        const int d = n2d[i];
        if (own[i] != 1 || d < 0 || d >= s.numOwnedDofs || !mask[d]) continue;
        lo = std::min(lo, p[i]);
        hi = std::max(hi, p[i]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &lo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &hi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return {lo, hi};
}

void steps(Stepper& s)
{
    for (int k = 0; k < 5; ++k)
    {
        runNsStep(s, .001, .1, 1.);
        require(s.lastPressureIters >= 0 && s.lastPressureIters < s.maxIter, "The pressure solve did not converge");
    }
}

void check_auto(const Fixture& f, int rank, int ranks)
{
    auto s = make_stepper(f, rank, ranks, Stepper::PressureIncrementMean::Auto);
    require(pressureHasGauge(*s), "The prescribed outlet pressure is not seen as a gauge");
    steps(*s);
    const auto p = outlet_pressure(*s);
    if (rank == 0) std::cout << std::scientific << "auto: outlet p in [" << p[0] << ", " << p[1] << "], T = " << TARGET << "\n";
    require(std::abs(p[0] - TARGET) < 1e-9 * TARGET && std::abs(p[1] - TARGET) < 1e-9 * TARGET,
            "Auto moved the outlet's pressure off its prescribed value");
}

void check_always(const Fixture& f, int rank, int ranks)
{
    auto s = make_stepper(f, rank, ranks, Stepper::PressureIncrementMean::Always);
    steps(*s);
    const auto p = outlet_pressure(*s);
    if (rank == 0) std::cout << std::scientific << "always: outlet p in [" << p[0] << ", " << p[1] << "], T = " << TARGET << "\n";
    require(std::abs(p[1] - TARGET) > 1e-2 * TARGET, "Always did not shift the outlet's pressure (the old behaviour)");
}

void check_velocity(const Fixture& f, int rank, int ranks)
{
    auto a = make_stepper(f, rank, ranks, Stepper::PressureIncrementMean::Auto);
    auto b = make_stepper(f, rank, ranks, Stepper::PressureIncrementMean::Always);
    steps(*a);
    steps(*b);
    double err = 0, scale = 0, lo = 1e300, hi = -1e300;
    for (auto field : {&Stepper::d_u, &Stepper::d_v, &Stepper::d_w})
    {
        const auto x = host((*a).*field);
        const auto y = host((*b).*field);
        for (size_t i = 0; i < x.size(); ++i)
        {
            err   = std::max(err, std::abs(x[i] - y[i]));
            scale = std::max(scale, std::abs(y[i]));
        }
    }
    const auto own = host(a->ownershipMap());
    const auto pa  = host(a->d_p);
    const auto pb  = host(b->d_p);
    for (size_t i = 0; i < pa.size(); ++i)
    {
        if (own[i] != 1) continue;
        lo = std::min(lo, pa[i] - pb[i]);
        hi = std::max(hi, pa[i] - pb[i]);
    }
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &scale, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &lo, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &hi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << std::scientific << "velocity: max |u_auto - u_always| / max |u| " << err / scale
                  << "; p_auto - p_always in [" << lo << ", " << hi << "]\n";
    require(scale > 0 && err <= 1e-8 * scale, "Auto and Always give different velocities");
    require(hi - lo <= 1e-6 * TARGET, "Auto and Always differ by more than a constant in pressure");
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
        Fixture fixture(rank, ranks, 6);
        if (mode == "auto")
            check_auto(fixture, rank, ranks);
        else if (mode == "always")
            check_always(fixture, rank, ranks);
        else if (mode == "velocity")
            check_velocity(fixture, rank, ranks);
        else
            require(false, "Unknown mode");
        cudaCheckCall(cudaDeviceSynchronize());
        if (rank == 0) std::cout << "PASS: pressure gauge " << argv[1] << " ranks=" << ranks << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
