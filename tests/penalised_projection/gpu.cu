#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"

#ifndef MARS_ENABLE_HYPRE
#error "Penalised projection checks require a GPU Hypre build of MARS"
#endif

// The penalised projection (setPenalisedProjection) on public Tet4 boxes, through the pump stepper:
//   baseline  W = 0 leaves the setup operator and the steps unchanged;
//   varying   piecewise W in slabs: every interior plane carries the flux of the layered operator;
//   update    W = w V gives a uniform beta = c / (c + w / eta): K becomes beta K0, after W, eta, dt change;
//   bdf2      the same at BDF2 startup: dtEff = dt on the first step, 2 dt / 3 after it;
//   flow      with inflow, the pressure solve converges, and continuity holds as well as without the
//             penalty (two separate checks);
//   reject    invalid W, eta, dt or configuration: std::invalid_argument on every rank, nothing changed.

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
// The pump path ldp_fsi uses: Hypre's assembled operator, DDT, divT gradient, do-nothing outlet.
std::unique_ptr<Stepper> make_stepper(const Fixture& f, int rank, int ranks, double inflow, bool bdf2)
{
    auto s               = std::make_unique<Stepper>(*f.domain, SolverKind::Hypre, 256, 2000, 1e-10, rank, ranks);
    s->Uinf              = inflow;
    s->uinfBase          = inflow;
    s->pumpZeroIC        = true;
    s->bcKind            = Stepper::BCKind::Pump;
    s->pressureSolve     = PressureSolveKind::DDT;
    s->useLegacyGradient = false;
    s->outletDoNothing   = true;
    s->useDirichletLift  = true;
    s->useBdf2           = bdf2;
    boundaries(*s, f);
    setupNSStepper(*s, .1, .001, CvfemKernelVariant::Tensor);
    return s;
}
std::vector<double> matrix(const Stepper& s) { return host(s.Apre.valuesPtr(), size_t(s.nnz)); }
// The last pressure solve succeeded: solveOneComponent returns -2 on failure, else its iterations,
// below the cap. (The cap alone accepts -2.)
bool pressure_solve_succeeded(const Stepper& s) { return s.lastPressureIters >= 0 && s.lastPressureIters < s.maxIter; }

// max |A - b K0| over the rows the projection assembles (owned, not boundary or pin rows), relative
// to max |K0|.
double scaled_operator_error(const Stepper& s, const std::vector<double>& K0, double b)
{
    const auto A = matrix(s);
    const auto row = host(s.d_rowPtr);
    const auto bdry = host(s.d_isPressureBdryDof);
    double err = 0, scale = 0;
    for (int r = 0; r < s.numOwnedDofs; ++r)
    {
        if (bdry[r] || r == s.pressurePinDof) continue;
        for (int k = row[r]; k < row[r + 1]; ++k)
        {
            err   = std::max(err, std::abs(A[k] - b * K0[k]));
            scale = std::max(scale, std::abs(K0[k]));
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &scale, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return err / scale;
}
// W = w V (the nodal volume): beta = c / (c + w invEta) on every node.
cstone::DeviceVector<double> proportional_W(const Stepper& s, double w)
{
    auto V = host(s.d_massNode);
    for (auto& v : V)
        v *= w;
    cstone::DeviceVector<double> d;
    upload(V, d);
    return d;
}
double uniform_beta(double dtEff, double w, double invEta) { return (1 / dtEff) / (1 / dtEff + w * invEta); }

void check_baseline(const Fixture& f, int rank, int ranks)
{
    // The operator: W = 0 refreshes it to the setup values.
    {
        auto s        = make_stepper(f, rank, ranks, 0, false);
        const auto K0 = matrix(*s);
        cstone::DeviceVector<double> W;
        upload(std::vector<double>(s->nodeCount, 0), W);
        setPenalisedProjection(*s, W.data(), 1e6);
        runNsStep(*s, .001, .1, 1.);
        require(scaled_operator_error(*s, K0, 1) < 1e-14, "W = 0 changed the pressure operator");
    }
    // The steps: with inflow, W = 0 gives the fields of a stepper without the penalty.
    auto a = make_stepper(f, rank, ranks, 1, true);
    auto b = make_stepper(f, rank, ranks, 1, true);
    cstone::DeviceVector<double> W;
    upload(std::vector<double>(a->nodeCount, 0), W);
    setPenalisedProjection(*a, W.data(), 1e6);
    for (int step = 0; step < 3; ++step)
    {
        runNsStep(*a, .001, .1, 1.);
        runNsStep(*b, .001, .1, 1.);
        require(pressure_solve_succeeded(*a) && pressure_solve_succeeded(*b), "The pressure solve did not converge");
    }
    double err = 0, scale = 0;
    for (auto field : {&Stepper::d_u, &Stepper::d_v, &Stepper::d_w, &Stepper::d_p})
    {
        const auto x = host((*a).*field);
        const auto y = host((*b).*field);
        for (size_t i = 0; i < x.size(); ++i)
        {
            err   = std::max(err, std::abs(x[i] - y[i]));
            scale = std::max(scale, std::abs(y[i]));
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &scale, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << std::scientific << "baseline: max |fields(W = 0) - fields(no penalty)| / max |fields| " << err / scale
                  << std::defaultfloat << "\n";
    // The operators agree to round-off (checked above); the fields then agree to the pressure solve's
    // tolerance (1e-10, relative), not to round-off: 6e-11 measured on GH200.
    require(scale > 0 && err <= 1e-8 * scale, "W = 0 changed the steps");
}

// Layered W, piecewise constant per plane of nodes, with jumps. For phi = phi(x) with slope g_L in
// layer L (between planes L and L + 1), the rows of an interior plane p sum to
//   sum_(i in p) (K phi)_i = (g_(p-1) G_(p-1) - g_p G_p) / h,   G_L = sum_(tets of L) gamma_e V_e,
// gamma_e the mean of the tet's nodal beta (the Galerkin operator's element coefficient). MARS's
// CVFEM K satisfies the Galerkin identity to ~1e-5 on this mesh even with a uniform beta (8.9e-6 on
// the 4^3 box), so the tolerance is 1e-4. That it still tells weightings apart is checked: the same
// identity with gamma_e the largest nodal beta must miss by more than 1e-2.
void check_varying(const Fixture& f, int rank, int ranks)
{
    auto s           = make_stepper(f, rank, ranks, 0, false);
    const int cells  = f.cells;
    const auto x = host(f.domain->getNodeX());
    const auto y = host(f.domain->getNodeY());
    const auto z = host(f.domain->getNodeZ());
    const double eta = 1e-2, h = 1. / cells;
    std::vector<double> w_plane(cells + 1);
    for (int p = 0; p <= cells; ++p)
        w_plane[p] = (p % 2 ? 40. : 3.) * (1 + p);
    {
        auto V = host(s->d_massNode);
        std::vector<double> W(s->nodeCount);
        for (size_t i = 0; i < W.size(); ++i)
            W[i] = w_plane[f.plane(x[i])] * V[i];
        cstone::DeviceVector<double> dW;
        upload(W, dW);
        setPenalisedProjection(*s, dW.data(), 1 / eta);
        runNsStep(*s, .001, .1, 1.);
    }
    const auto beta = host(s->d_penaltyBeta);
    // G_L over this rank's owned elements (each element counted once).
    std::vector<double> G(cells, 0), Gmax(cells, 0);
    {
        const auto& conn = f.domain->getElementToNodeConnectivity();
        std::array<std::vector<Key>, 4> c{host(std::get<0>(conn)), host(std::get<1>(conn)), host(std::get<2>(conn)),
                                          host(std::get<3>(conn))};
        const size_t start = f.domain->startIndex(), end = start + f.domain->localElementCount();
        for (size_t e = start; e < end; ++e)
        {
            Point v[4];
            double gamma = 0, gammaMax = 0;
            int lowest   = cells;
            for (int a = 0; a < 4; ++a)
            {
                const auto n = c[a][e];
                v[a]         = Point{x[n], y[n], z[n]};
                gamma += beta[n] / 4;
                gammaMax = std::max(gammaMax, beta[n]);
                lowest = std::min(lowest, f.plane(x[n]));
            }
            const double volume = std::abs(determinant(v[0], v[1], v[2], v[3])) / 6;
            G[lowest] += gamma * volume;
            Gmax[lowest] += gammaMax * volume;
        }
        MPI_Allreduce(MPI_IN_PLACE, G.data(), cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, Gmax.data(), cells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    // phi(x) with a slope per layer; the sums of the owned rows of each plane.
    std::vector<double> slope(cells), phi(cells + 1, 0);
    for (int L = 0; L < cells; ++L)
    {
        slope[L]   = 1 + 0.37 * L * L;
        phi[L + 1] = phi[L] + slope[L] * h;
    }
    const auto A = matrix(*s);
    const auto row = host(s->d_rowPtr);
    const auto col = host(s->d_colInd);
    const auto dof = host(s->d_node_to_dof);
    std::vector<int> dof_plane(s->numTotalDofs, -1);
    for (size_t i = 0; i < dof.size(); ++i)
        if (dof[i] >= 0) dof_plane[dof[i]] = f.plane(x[i]);
    std::vector<double> sum(cells + 1, 0);
    for (size_t i = 0; i < dof.size(); ++i)
    {
        const int r = dof[i];
        if (r < 0 || r >= s->numOwnedDofs) continue;
        double value = 0;
        for (int k = row[r]; k < row[r + 1]; ++k)
            value += A[k] * phi[dof_plane[col[k]]];
        sum[f.plane(x[i])] += value;
    }
    MPI_Allreduce(MPI_IN_PLACE, sum.data(), cells + 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // K's sign convention, from its interior diagonal (the identity holds for +stiffness).
    double diagonal = 0;
    {
        const auto diag = host(s->d_diagPtr);
        const auto bdry = host(s->d_isPressureBdryDof);
        for (int r = 0; r < s->numOwnedDofs; ++r)
            if (!bdry[r] && r != s->pressurePinDof && diag[r] >= 0) diagonal += A[diag[r]];
        MPI_Allreduce(MPI_IN_PLACE, &diagonal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    const double sign = diagonal > 0 ? 1. : -1.;
    double err = 0, errMax = 0, scale = 0;
    for (int p = 1; p < cells; ++p)
    {
        const double expected = sign * (slope[p - 1] * G[p - 1] - slope[p] * G[p]) / h;
        const double wrong    = sign * (slope[p - 1] * Gmax[p - 1] - slope[p] * Gmax[p]) / h;
        err                   = std::max(err, std::abs(sum[p] - expected));
        errMax                = std::max(errMax, std::abs(sum[p] - wrong));
        scale                 = std::max(scale, std::abs(expected));
    }
    if (rank == 0)
        std::cout << std::scientific << "varying: plane balance error / largest balance " << err / scale
                  << " (with the largest nodal beta per tet instead of the mean: " << errMax / scale << ")"
                  << std::defaultfloat << "\n";
    require(scale > 0 && err <= 1e-4 * scale, "Layered penalty: a plane's flux differs from the layered operator's");
    require(errMax > 1e-2 * scale, "The plane balance does not tell the element coefficient apart");
}

void check_update(const Fixture& f, int rank, int ranks)
{
    auto s        = make_stepper(f, rank, ranks, 0, false);
    const auto K0 = matrix(*s);
    const struct
    {
        double w, eta, dt;
    } cases[] = {{1., 1e-2, .001}, {5., 1e-2, .001}, {5., 1e-3, .001}, {5., 1e-3, .002}, {0., 1e-3, .002}};
    for (const auto& c : cases)
    {
        auto W = proportional_W(*s, c.w);
        setPenalisedProjection(*s, W.data(), 1 / c.eta);
        runNsStep(*s, c.dt, .1, 1.);
        require(scaled_operator_error(*s, K0, uniform_beta(c.dt, c.w, 1 / c.eta)) < 1e-12,
                "The operator did not follow a change of W, eta or dt");
    }
    // A change of dt alone (the same W and eta) rebuilds it too.
    runNsStep(*s, .004, .1, 1.);
    require(scaled_operator_error(*s, K0, uniform_beta(.004, 0., 1e3)) < 1e-12, "dt alone did not rebuild");
}

void check_bdf2(const Fixture& f, int rank, int ranks)
{
    auto s        = make_stepper(f, rank, ranks, 0, true);
    const auto K0 = matrix(*s);
    const double w = 3, invEta = 1e2, dt = .001;
    auto W = proportional_W(*s, w);
    setPenalisedProjection(*s, W.data(), invEta);
    runNsStep(*s, dt, .1, 1.);
    require(scaled_operator_error(*s, K0, uniform_beta(dt, w, invEta)) < 1e-12, "BDF2's first step is not at dt");
    require(s->bdfStep >= 1, "BDF2 did not start");
    for (int step = 0; step < 2; ++step)
    {
        runNsStep(*s, dt, .1, 1.);
        require(scaled_operator_error(*s, K0, uniform_beta(2 * dt / 3, w, invEta)) < 1e-12,
                "BDF2's later steps are not at 2 dt / 3");
    }
}

void check_flow(const Fixture& f, int rank, int ranks)
{
    auto with = make_stepper(f, rank, ranks, 1, true);
    auto without = make_stepper(f, rank, ranks, 1, true);
    auto W = proportional_W(*with, 1);
    setPenalisedProjection(*with, W.data(), 1e3);  // beta = 0.5 at dt = 1e-3
    for (int step = 0; step < 4; ++step)
    {
        runNsStep(*with, .001, .1, 1.);
        runNsStep(*without, .001, .1, 1.);
        // Pressure convergence: the solve stopped below the iteration cap.
        require(pressure_solve_succeeded(*with) && pressure_solve_succeeded(*without),
                "The pressure solve did not converge");
    }
    // Continuity: the post-corrector divergence is no worse than without the penalty.
    if (rank == 0)
        std::cout << "flow: |div u| max " << with->lastDivMax << " with the penalty, " << without->lastDivMax
                  << " without; pressure iterations " << with->lastPressureIters << ", " << without->lastPressureIters
                  << "\n";
    require(std::isfinite(with->lastDivMax) && with->lastDivMax <= 2 * without->lastDivMax + 1e-12,
            "Continuity is worse with the penalty");
}

void rejected_without_mutation(Stepper& s, double dt)
{
    const auto p = host(s.d_p);
    const auto u = host(s.d_u);
    const auto v = host(s.d_v);
    const auto w = host(s.d_w);
    const auto A = matrix(s);
    const int history = s.bdfStep;
    const auto nu     = s.nuCached;
    bool rejected     = false;
    try
    {
        runNsStep(s, dt, .37, 1.);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "An invalid penalised projection was accepted on a rank");
    require(p == host(s.d_p) && u == host(s.d_u) && v == host(s.d_v) && w == host(s.d_w) && A == matrix(s) &&
                s.bdfStep == history && s.nuCached == nu,
            "A rejected step changed the state");
}
void setter_rejected(Stepper& s, const double* W, double invEta)
{
    const auto W0     = host(s.d_penaltyW);
    const auto inv0   = s.penaltyInvEta;
    const bool on0    = s.penalisedProjection;
    bool rejected     = false;
    try
    {
        setPenalisedProjection(s, W, invEta);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "setPenalisedProjection accepted invalid input on a rank");
    require(W0 == host(s.d_penaltyW) && inv0 == s.penaltyInvEta && on0 == s.penalisedProjection,
            "A rejected setPenalisedProjection changed the state");
}
void check_rejections(const Fixture& f, int rank, int ranks)
{
    auto s            = make_stepper(f, rank, ranks, 0, false);
    const bool inject = rank == ranks - 1;  // invalid on one rank only: the refusal must be collective
    auto valid        = proportional_W(*s, 1);
    // The setter: W NaN or negative, invEta NaN, negative or infinite.
    for (double bad : {double(NAN), -0.5})
    {
        auto W = host(valid);
        if (inject && !W.empty()) W[W.size() / 2] = bad;
        cstone::DeviceVector<double> dW;
        upload(W, dW);
        setter_rejected(*s, dW.data(), 1.);
    }
    for (double invEta : {double(NAN), -1., double(INFINITY)})
        setter_rejected(*s, valid.data(), inject ? invEta : 1.);
    // The step: dt zero, negative or NaN; W changed behind the setter's back; a wrong size.
    setPenalisedProjection(*s, valid.data(), 1.);
    for (double dt : {0., -.001, double(NAN)})
        rejected_without_mutation(*s, dt);
    {
        auto W = host(s->d_penaltyW);
        if (inject && !W.empty()) W[0] = -1;
        upload(W, s->d_penaltyW);
        rejected_without_mutation(*s, .001);
        setPenalisedProjection(*s, valid.data(), 1.);
    }
    if (inject) s->d_penaltyW.resize(s->nodeCount + 1);
    rejected_without_mutation(*s, .001);
    setPenalisedProjection(*s, valid.data(), 1.);
    // Unsupported configurations.
    if (inject) s->solverKind = SolverKind::CG;
    rejected_without_mutation(*s, .001);
    s->solverKind = SolverKind::Hypre;
    for (bool* option : {&s->useFemProjection, &s->useRhieChow, &s->useVMSStab, &s->useRcImplicit, &s->useRcOnly})
    {
        if (inject) *option = true;
        rejected_without_mutation(*s, .001);
        *option = false;
    }
    if (inject) s->outletBeta = 0;
    rejected_without_mutation(*s, .001);
    s->outletBeta = -1;
    if (inject) setenv("MARS_HYPRE_USE_DDT", "1", 1);
    rejected_without_mutation(*s, .001);
    unsetenv("MARS_HYPRE_USE_DDT");
    // And back to a valid state, which steps.
    runNsStep(*s, .001, .1, 1.);
}

// A deliberately failed pressure solve, the velocity's solves succeeding (so the pressure's RHS, the
// inflow's divergence, is not 0): a step of 1e-8 (the velocity matrix is then its mass: its solves
// converge in a few iterations), at most 5 iterations, Jacobi in place of BoomerAMG for the pressure
// and an acceptance residual nothing reaches. The solve returns -2, and pressure_solve_succeeded, the
// flow and baseline checks' test, rejects it.
void check_failed_solve(const Fixture& f, int rank, int ranks)
{
    auto s = make_stepper(f, rank, ranks, 1, false);
    auto W = proportional_W(*s, 1);
    setPenalisedProjection(*s, W.data(), 1e3);
    s->maxIter = 5;
    setenv("MARS_HYPRE_PRECOND", "jacobi", 1);
    setenv("MARS_HYPRE_ACCEPT_RES", "1e-300", 1);
    runNsStep(*s, 1e-8, .1, 1.);
    unsetenv("MARS_HYPRE_PRECOND");
    unsetenv("MARS_HYPRE_ACCEPT_RES");
    if (rank == 0)
        std::cout << "failed: velocity iterations " << s->lastUIters << " " << s->lastVIters << " " << s->lastWIters << "\n";
    require(s->lastUIters >= 0 && s->lastVIters >= 0 && s->lastWIters >= 0,
            "The forced failure's velocity solves did not succeed");
    if (rank == 0) std::cout << "failed: pressure solve returned " << s->lastPressureIters << "\n";
    require(s->lastPressureIters == -2, "The forced pressure failure did not return -2");
    require(!pressure_solve_succeeded(*s), "The convergence check accepted a failed pressure solve");
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
        // 6^3 cubes: on 4^3, BoomerAMG under GMRES fails the first solve with inflow (Hypre error 1,
        // with or without the penalty); 6^3 and 8^3 converge in about 12 iterations at 1, 2 and 4 ranks.
        Fixture fixture(rank, ranks, 6);
        if (mode == "baseline")
            check_baseline(fixture, rank, ranks);
        else if (mode == "varying")
            check_varying(fixture, rank, ranks);
        else if (mode == "update")
            check_update(fixture, rank, ranks);
        else if (mode == "bdf2")
            check_bdf2(fixture, rank, ranks);
        else if (mode == "flow")
            check_flow(fixture, rank, ranks);
        else if (mode == "failed")
            check_failed_solve(fixture, rank, ranks);
        else if (mode == "reject")
            check_rejections(fixture, rank, ranks);
        else
            require(false, "Unknown mode");
        cudaCheckCall(cudaDeviceSynchronize());
        if (rank == 0) std::cout << "PASS: penalised projection " << argv[1] << " ranks=" << ranks << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
