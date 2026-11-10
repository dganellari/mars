#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"

#ifndef MARS_ENABLE_HYPRE
#error "Fluid contract pressure checks require a GPU Hypre build of MARS"
#endif

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
struct Fixture
{
    std::unique_ptr<Domain> domain;
    double x0 = 0.123456789, x1 = 1.253456789;
    Fixture(int rank, int ranks, int cells = 4)
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
};
void boundaries(Stepper& s, const Fixture& f)
{
    auto x = host(f.domain->getNodeX()), y = host(f.domain->getNodeY()), z = host(f.domain->getNodeZ());
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

void rejected_without_mutation(Stepper& s)
{
    const auto p = host(s.d_p), u = host(s.d_u), v = host(s.d_v), w = host(s.d_w);
    const int history = s.bdfStep;
    const auto nu     = s.nuCached;
    bool rejected     = false;
    try
    {
        runNsStep(s, .001, .37, 1.);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "Unsupported pressure configuration was accepted on a rank");
    require(p == host(s.d_p) && u == host(s.d_u) && v == host(s.d_v) && w == host(s.d_w) && s.bdfStep == history &&
                s.nuCached == nu,
            "Rejected step mutated state");
}
void check_rejections(Stepper& s, int rank, int ranks)
{
    s.outletPressureLift = true;
    upload(std::vector<double>(s.numOwnedDofs, 0), s.d_pPhiTargetDof);
    const bool inject = rank == ranks - 1;
    if (inject) s.d_pPhiTargetDof.resize(s.numOwnedDofs + 1);
    rejected_without_mutation(s);
    upload(std::vector<double>(s.numOwnedDofs, 0), s.d_pPhiTargetDof);
    auto mask = host(s.d_isPressureBdryDof);
    int owner = ranks;
    for (auto value : mask)
        if (value) owner = rank;
    MPI_Allreduce(MPI_IN_PLACE, &owner, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    require(owner < ranks, "Fixture has no pressure target");
    auto target = host(s.d_pPhiTargetDof);
    if (rank == owner)
        for (size_t i = 0; i < mask.size(); ++i)
            if (mask[i])
            {
                target[i] = NAN;
                break;
            }
    upload(target, s.d_pPhiTargetDof);
    rejected_without_mutation(s);
    upload(std::vector<double>(s.numOwnedDofs, 0), s.d_pPhiTargetDof);
    if (inject) s.solverKind = SolverKind::CG;
    rejected_without_mutation(s);
    s.solverKind = SolverKind::Hypre;
    if (inject) s.pressureSolve = PressureSolveKind::K;
    rejected_without_mutation(s);
    s.pressureSolve = PressureSolveKind::DDT;
    if (inject) s.nCorrectors = 2;
    rejected_without_mutation(s);
    s.nCorrectors = 1;
    for (bool* option : {&s.rotationalPressureCorrection, &s.useFemProjection, &s.fluxNeumann})
    {
        if (inject) *option = true;
        rejected_without_mutation(s);
        *option = false;
    }
    if (inject) s.outletDoNothing = false;
    rejected_without_mutation(s);
    s.outletDoNothing = true;
    if (inject) s.outletBeta = 0;
    rejected_without_mutation(s);
    s.outletBeta = -1;
    if (inject) s.pumpDp = 1;
    rejected_without_mutation(s);
    s.pumpDp = 0;
    if (inject) setenv("MARS_HYPRE_USE_DDT", "1", 1);
    rejected_without_mutation(s);
    unsetenv("MARS_HYPRE_USE_DDT");
}
void check_pressure(Stepper& s)
{
    require(!s.outletPressureLift, "Prescribed pressure must default off");
    runNsStep(s, .001, .1, 1.);
    for (auto* field : {&s.d_u, &s.d_v, &s.d_w, &s.d_p})
    {
        bool zero = true;
        for (double value : host(*field))
            zero &= value == 0;
        require(zero, "Default quiescent state changed");
    }
    s.outletPressureLift = true;
    // The override must not force a nonsymmetric row system through PCG.
    setenv("MARS_HYPRE_KRYLOV", "pcg", 1);
    for (double target : {2., -1., 0.})
    {
        auto mask = host(s.d_isPressureBdryDof);
        std::vector<double> targets(s.numOwnedDofs, NAN);
        for (size_t i = 0; i < targets.size(); ++i)
            if (mask[i]) targets[i] = target;
        upload(targets, s.d_pPhiTargetDof);
        runNsStep(s, .001, .1, 1.);
        bool correct = true;
        for (double value : host(s.d_p))
            correct &= std::isfinite(value) && std::abs(value - target) < 1e-7;
        require(correct, "Target did not propagate into the interior or return to zero");
    }
    unsetenv("MARS_HYPRE_KRYLOV");
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
        require(devices > 0 && argc == 2, "Expected GPU and check mode");
        cudaCheckCall(cudaSetDevice(local_rank % devices));
        mars::requestFullNodeHaloExchange();
        unsetenv("MARS_HYPRE_USE_DDT");
        const std::string mode(argv[1]);
        require(mode == "pressure" || mode == "reject" || mode == "refined", "Unknown mode");
        Fixture fixture(rank, ranks, mode == "refined" ? 8 : 4);
        Stepper s{*fixture.domain, SolverKind::Hypre, 256, 2000, 1e-10, rank, ranks};
        s.Uinf              = 0;
        s.pumpZeroIC        = true;
        s.bcKind            = Stepper::BCKind::Pump;
        s.pressureSolve     = PressureSolveKind::DDT;
        s.useLegacyGradient = false;
        s.outletDoNothing   = true;
        s.useDirichletLift  = true;
        boundaries(s, fixture);
        setupNSStepper(s, .1, .001, CvfemKernelVariant::Tensor);
        if (mode == "reject")
            check_rejections(s, rank, ranks);
        else
            check_pressure(s);
        cudaCheckCall(cudaDeviceSynchronize());
        if (rank == 0) std::cout << "PASS: outlet pressure " << argv[1] << " ranks=" << ranks << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
