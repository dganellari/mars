#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"

#ifndef MARS_ENABLE_HYPRE
#error "Slip-wall checks require a GPU Hypre build of MARS"
#endif

// Free-slip walls (NSStepper::slipNodes) on the 2D Taylor-Green vortex, an exact Navier-Stokes
// solution in [0, pi]^2 with free-slip walls on its four sides:
//   u = sin x cos y e^(-2 nu t),  v = -cos x sin y e^(-2 nu t),  p = (cos 2x + cos 2y) / 4 e^(-4 nu t)
// (rho 1), in a slab two cells thick, free slip on its faces too (w = 0). Every face is a slip wall:
// no openings, the pressure pinned at one interior node (fluxNeumann). Hypre, assembled DDT, divT,
// implicit advection, the face-flux projection. (Not the nodal projection: in a closed box its
// pressure solve does not converge, GMRES stalling with Jacobi and BoomerAMG returning a
// null-space-contaminated solution; the interior pin evidently does not reach the operator it solves.)
//   {tet,hex}_ff         the face-flux projection and its convection (the stored flux initialised
//                        from the initial velocity, initFaceFluxFromVelocity). On 16^2 and 32^2
//                        cells, every step: the pressure
//                        solve succeeds, the fields are finite, and the normal component is exactly 0
//                        on every slip plane. At the end, the velocity's relative L2 error against
//                        the exact solution is small and falls with h at least as h^1.5; the kinetic
//                        energy follows e^(-4 nu t).
//   reject               slip nodes with BCKind::Channel: std::invalid_argument on every rank.

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
bool pressure_solve_succeeded(int iters, int cap) { return iters >= 0 && iters < cap; }

constexpr double NU = 0.1;

// The slab [0, pi] x [0, pi] x [0, 2 pi / n] of n x n x 2 cubes: one Hex8 each, or five Tet4
// (alternating, as tests/outlet_pressure); elements dealt round-robin to the ranks.
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
std::unique_ptr<mars::ElementDomain<Tag, Real, Key, cstone::execution::Gpu>> make_domain(int n, int rank, int ranks)
{
    using Domain   = mars::ElementDomain<Tag, Real, Key, cstone::execution::Gpu>;
    constexpr bool hex = std::is_same_v<Tag, mars::HexTag>;
    typename Domain::HostCoordsTuple xyz;
    typename Domain::HostConnectivityTuple conn;
    const int nz = 2;
    const double h = M_PI / n;
    auto index     = [n](int x, int y, int z) { return Key(x + (n + 1) * (y + (n + 1) * z)); };
    std::vector<std::array<double, 3>> points;
    for (int z = 0; z <= nz; ++z)
        for (int y = 0; y <= n; ++y)
            for (int x = 0; x <= n; ++x)
            {
                points.push_back({x * h, y * h, z * h});
                std::get<0>(xyz).push_back(x * h);
                std::get<1>(xyz).push_back(y * h);
                std::get<2>(xyz).push_back(z * h);
            }
    int element = 0;
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x)
            {
                if constexpr (hex)
                {
                    const Key c[8] = {index(x, y, z),         index(x + 1, y, z),         index(x + 1, y + 1, z),
                                      index(x, y + 1, z),     index(x, y, z + 1),         index(x + 1, y, z + 1),
                                      index(x + 1, y + 1, z + 1), index(x, y + 1, z + 1)};
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
                            nodes[a]    = index(x + (c & 1), y + ((c >> 1) & 1), z + ((c >> 2) & 1));
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

template<class Tag>
struct Slab
{
    using Domain  = mars::ElementDomain<Tag, Real, Key, cstone::execution::Gpu>;
    using Stepper = NSStepper<Key, Real, Tag>;
    int n;
    std::unique_ptr<Domain> domain;
    std::unique_ptr<Stepper> s;
    Slab(int n_, int rank, int ranks, bool faceFlux, typename Stepper::BCKind bc = Stepper::BCKind::Pump) : n(n_)
    {
        domain               = make_domain<Tag>(n, rank, ranks);
        s                    = std::make_unique<Stepper>(*domain, SolverKind::Hypre, 256, 2000, 1e-10, rank, ranks);
        s->Uinf              = 0;
        s->uinfBase          = 0;
        s->pumpZeroIC        = true;
        s->bcKind            = bc;
        s->pressureSolve     = PressureSolveKind::DDT;
        s->useLegacyGradient = false;
        s->outletDoNothing   = true;
        s->useDirichletLift  = true;
        s->implicitAdvection = true;
        s->fluxNeumann       = true;  // a closed box: one interior pressure pin
        s->useFaceFlux       = faceFlux;
        s->faceFluxConvection = faceFlux;
        // Slip planes, within 1e-5 (the domain's coordinates are SFC-quantized, ~1e-6).
        const auto x = host(domain->getNodeX()), y = host(domain->getNodeY()), z = host(domain->getNodeZ());
        const double zmax = 2 * M_PI / n;
        for (size_t i = 0; i < x.size(); ++i)
        {
            if (std::abs(x[i]) < 1e-5 || std::abs(x[i] - M_PI) < 1e-5) s->slipNodes[0].push_back(int(i));
            if (std::abs(y[i]) < 1e-5 || std::abs(y[i] - M_PI) < 1e-5) s->slipNodes[1].push_back(int(i));
            if (std::abs(z[i]) < 1e-5 || std::abs(z[i] - zmax) < 1e-5) s->slipNodes[2].push_back(int(i));
        }
        setupNSStepper(*s, NU, .002, CvfemKernelVariant::Tensor);
        // The exact field at t = 0 on every local node, and the face flux from it.
        std::vector<double> u(x.size()), v(x.size()), w(x.size(), 0.0), p(x.size());
        for (size_t i = 0; i < x.size(); ++i)
        {
            u[i] = std::sin(x[i]) * std::cos(y[i]);
            v[i] = -std::cos(x[i]) * std::sin(y[i]);
            p[i] = .25 * (std::cos(2 * x[i]) + std::cos(2 * y[i]));
        }
        upload(u, s->d_u);
        upload(v, s->d_v);
        upload(w, s->d_w);
        upload(p, s->d_p);
        initFaceFluxFromVelocity(*s);
    }
    // The largest |normal component| over the slip nodes (owned and ghost), over all ranks.
    double slip_violation() const
    {
        const cstone::DeviceVector<Real>* q[3] = {&s->d_u, &s->d_v, &s->d_w};
        double worst = 0;
        for (int a = 0; a < 3; ++a)
        {
            const auto m = host(s->d_isSlipNode[a]);
            const auto c = host(*q[a]);
            for (size_t i = 0; i < m.size(); ++i)
                if (m[i]) worst = std::max(worst, std::abs(c[i]));
        }
        MPI_Allreduce(MPI_IN_PLACE, &worst, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return worst;
    }
    // The relative L2 error of the velocity against the exact one at t, lumped over the owned nodes,
    // and the kinetic energy's ratio to the exact one.
    std::array<double, 2> errors(double t) const
    {
        const auto own = host(s->ownershipMap());
        const auto V = host(s->d_massNode);
        const auto x = host(domain->getNodeX()), y = host(domain->getNodeY());
        const auto u = host(s->d_u), v = host(s->d_v), w = host(s->d_w);
        const double decay = std::exp(-2 * NU * t);
        double sums[3] = {0, 0, 0};  // |u - u_ex|^2, |u_ex|^2, |u|^2 (each times V)
        for (size_t i = 0; i < u.size(); ++i)
        {
            if (own[i] != 1) continue;
            const double ue = std::sin(x[i]) * std::cos(y[i]) * decay, ve = -std::cos(x[i]) * std::sin(y[i]) * decay;
            sums[0] += V[i] * ((u[i] - ue) * (u[i] - ue) + (v[i] - ve) * (v[i] - ve) + w[i] * w[i]);
            sums[1] += V[i] * (ue * ue + ve * ve);
            sums[2] += V[i] * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
        }
        MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return {std::sqrt(sums[0] / sums[1]), sums[2] / sums[1]};
    }
};

template<class Tag>
double run(int n, int rank, int ranks, bool faceFlux)
{
    Slab<Tag> c(n, rank, ranks, faceFlux);
    const double dt = .002;
    const int steps = 25;
    for (int step = 0; step < steps; ++step)
    {
        runNsStep(*c.s, dt, NU, 1.);
        require(pressure_solve_succeeded(c.s->lastPressureIters, c.s->maxIter), "The pressure solve did not converge");
        require(c.slip_violation() == 0, "A slip node has a normal velocity");
        bool finite = true;
        for (auto* f : {&c.s->d_u, &c.s->d_v, &c.s->d_w, &c.s->d_p})
            for (double q : host(*f))
                finite &= std::isfinite(q);
        require(finite, "A field is not finite");
    }
    const auto e = c.errors(steps * dt);
    if (rank == 0)
        std::cout << std::scientific << "n " << n << ": relative L2 velocity error " << e[0] << ", KE / KE_exact " << e[1]
                  << std::defaultfloat << "\n";
    require(std::abs(e[1] - 1) < 2e-2, "The kinetic energy does not follow the exact decay");
    return e[0];
}

template<class Tag>
void check_convergence(int rank, int ranks, bool faceFlux)
{
    const double coarse = run<Tag>(16, rank, ranks, faceFlux);
    const double fine   = run<Tag>(32, rank, ranks, faceFlux);
    const double order  = std::log(coarse / fine) / std::log(2.);
    if (rank == 0) std::cout << "observed order " << order << "\n";
    require(coarse < 5e-2, "The 16^2 velocity error is too large");
    require(order > 1.5, "The velocity error does not fall with h");
}

void check_reject(int rank, int ranks)
{
    bool rejected = false;
    try
    {
        Slab<mars::HexTag> c(8, rank, ranks, false, NSStepper<Key, Real, mars::HexTag>::BCKind::Channel);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "Slip walls were accepted without the pump boundary conditions");
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
        if (mode == "tet_ff")
            check_convergence<mars::TetTag>(rank, ranks, true);
        else if (mode == "hex_ff")
            check_convergence<mars::HexTag>(rank, ranks, true);
        else if (mode == "reject")
            check_reject(rank, ranks);
        else
            require(false, "Unknown mode");
        cudaCheckCall(cudaDeviceSynchronize());
        if (rank == 0) std::cout << "PASS: slip walls " << argv[1] << " ranks=" << ranks << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
