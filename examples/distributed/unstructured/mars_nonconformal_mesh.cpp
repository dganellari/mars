#include "mars_nonconformal_scalar_runtime.hpp"
#include <bit>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "mars_nonconformal_mesh_input.hpp"

struct VolumeStencil
{
    const Vec3<double>* x;
    const int* cells;
    const double* k;
    InterfaceStencil<double>* points;
    int64_t first, last;
    MARS_FS_FN void operator()(int p) const
    {
        int e = p / 4, a = p % 4;
        auto& s = points[p];
        s       = {};
        s.row   = cells[p];
        if (s.row < first || s.row >= last) return;
        InterfaceTetFace<double> f{};
        f.opposite    = 3;
        f.diffusivity = k[e];
        for (int j = 0; j < 4; ++j)
        {
            f.dof[j] = cells[4 * e + j];
            f.x[j]   = x[f.dof[j]];
        }
        nonconformal_detail::Geometry<double> g;
        s.status = nonconformal_detail::geometry(f, g);
        if (s.status != InterfaceStatus::valid) return;
        double volume =
            std::abs(dot(sub(f.x[1], f.x[0]), nonconformal_detail::cross(sub(f.x[2], f.x[0]), sub(f.x[3], f.x[0])))) /
            6.;
        for (int j = 0; j < 8; ++j)
        {
            s.column[j] = f.dof[j % 4];
            s.value[j]  = j < 4 ? k[e] * volume * dot(g.gradient[a], g.gradient[j]) : 0.;
        }
    }
};
struct OwnedPoint
{
    int64_t first, last;
    MARS_FS_FN bool operator()(const InterfaceStencil<double>& s) const { return s.row >= first && s.row < last; }
};
static Graph
assemble_mesh(const Mesh& m, const std::vector<int64_t>& rows, const Buffer<InterfaceStencil<double>>& interface)
{
    Buffer<Vec3<double>> coords(m.x.begin(), m.x.end());
    Buffer<int> cells(m.cells.begin(), m.cells.end());
    Buffer<double> k(m.k.begin(), m.k.end());
    Buffer<InterfaceStencil<double>> all(m.cells.size()), points(m.cells.size() + interface.size());
    int64_t first = int64_t(m.x.size()) * rankId / rankCount, last = int64_t(m.x.size()) * (rankId + 1) / rankCount;
    each<device>(int(all.size()),
                 VolumeStencil{pointer(coords), pointer(cells), pointer(k), pointer(all), first, last});
    completed();
#if defined(__CUDACC__) || defined(__HIPCC__)
    auto end     = thrust::copy_if(all.begin(), all.end(), points.begin(), OwnedPoint{first, last});
    size_t count = end - points.begin();
    thrust::copy(interface.begin(), interface.end(), end);
#else
    auto end     = std::copy_if(all.begin(), all.end(), points.begin(), OwnedPoint{first, last});
    size_t count = end - points.begin();
    std::copy(interface.begin(), interface.end(), end);
#endif
    points.resize(count + interface.size());
    Graph base, out;
    base.rowIds = Buffer<int64_t>(rows.begin(), rows.end());
    base.offsets.resize(rows.size() + 1, 0);
#if defined(__CUDACC__) || defined(__HIPCC__)
    check(buildNonconformalScalarGraph(base.view(), int(points.size()), pointer(points), out),
          "volume/interface graph");
    Buffer<int> valid(rows.size(), 0);
    check(assembleNonconformalScalarRows(out.view(), int(points.size()), pointer(points), pointer(out.pointOffsets),
                                         pointer(out.pointIds), pointer(valid)) == cudaSuccess,
          "scalar assembly launch");
    completed();
    check(thrust::count(valid.begin(), valid.end(), 1) == int(rows.size()), "scalar rows accepted");
#else
    check(buildNonconformalScalarGraphHost(base.view(), int(points.size()), pointer(points), out),
          "volume/interface graph");
    check(assembleNonconformalScalarRowsHost(out.view(), int(points.size()), pointer(points), pointer(out.pointOffsets),
                                             pointer(out.pointIds)),
          "scalar rows accepted");
#endif
    return out;
}

static void solve_mesh(const Mesh& m, double tolerance, const std::string& output, bool analytic)
{
    auto rows = ownedRows(m);
    auto f    = extract(m);
    NonconformalInterfacePlan<double, device> plan(MPI_COMM_WORLD);
    check(plan.build(f.faces, f.ids, f.owners, rows, peers(), tolerance), "interface coverage and pairing");
    auto a                  = assemble_mesh(m, rows, plan.stencils());
    Buffer<double> original = a.values;
    constexpr int restart   = 100;
    ChannelOps ops(a, m, restart);
    if (rankId == 0) std::puts("Interface paired and scalar system assembled; solving diffusion.");
    auto result = outlet_fgmres(ops, restart, 4000, 1e-11, std::numeric_limits<double>::epsilon());
    completed();
    check(result.converged, "mesh scalar solve converged");
    // Refresh solution ghosts once before independent residuals and output.
    ops.apply(ops.solution(), ops.work());
    completed();
    auto values = host<double>(original), x = host<double>(ops.x);
    auto offsets = host<int>(a.offsets), columns = host<int>(ops.columns);
    auto global_columns = host<int64_t>(a.columns);
    std::map<int64_t, int> slot;
    for (size_t r = 0; r < rows.size(); ++r)
        slot[rows[r]] = int(r);
    for (size_t j = 0; j < columns.size(); ++j)
        slot[global_columns[j]] = columns[j];
    long double res2 = 0, b2 = 0, reactions[2]{}, exact2 = 0;
    double error = 0;
    for (size_t r = 0; r < rows.size(); ++r)
    {
        int64_t id      = rows[r];
        long double raw = 0, exact = 0;
        for (int j = offsets[r]; j < offsets[r + 1]; ++j)
        {
            raw += static_cast<long double>(values[j]) * x[columns[j]];
            if (analytic) exact += static_cast<long double>(values[j]) * (m.x[global_columns[j]].x / 2.);
        }
        int bc = m.boundary[id];
        if (bc)
        {
            reactions[bc - 1] += raw;
            raw = x[r] - double(bc == 2);
            b2 += double(bc == 2);
        }
        if (analytic)
        {
            if (bc) exact = m.x[id].x / 2. - double(bc == 2);
            exact2 += exact * exact;
            error = std::max(error, std::abs(x[r] - m.x[id].x / 2.));
        }
        res2 += raw * raw;
    }
    long double flux[2]{};
    for (const auto& s : host<InterfaceStencil<double>>(plan.stencils()))
        for (int j = 0; j < 8; ++j)
            flux[m.side[s.row]] += static_cast<long double>(s.value[j]) * x[slot.at(s.column[j])];
    long double local[]{res2, b2, reactions[0], reactions[1], flux[0], flux[1], exact2}, global[7]{};
    MPI_Allreduce(local, global, 7, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double max_error = 0;
    MPI_Allreduce(&error, &max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    double relative = std::sqrt(double(global[0] / global[1]));
    check(std::isfinite(relative) && relative < 2e-11, "independent mesh residual");
    double flux_scale = std::max({1., std::abs(double(global[4])), std::abs(double(global[5]))});
    check(std::abs(global[4] + global[5]) < 2e-11 * flux_scale, "shared interface flux balance");
    check(std::abs(global[2] + global[3]) < 2e-8 * flux_scale, "outer boundary balance");
    if (analytic)
    {
        check(std::sqrt(double(global[6] / global[1])) < 1e-12, "u=x/2 satisfies this discrete boundary problem");
        check(max_error < 2e-8, "linear channel solution");
        check(std::abs(global[2] + .5) < 2e-8 && std::abs(global[3] - .5) < 2e-8 && std::abs(global[4] + .5) < 2e-8 &&
                  std::abs(global[5] - .5) < 2e-8,
              "unit-area channel transmitted flux");
    }
    if (!output.empty())
    {
        std::ofstream out(output + ".rank" + std::to_string(rankId) + ".csv");
        out << "node_id,x,y,z,scalar\n" << std::setprecision(17);
        for (size_t r = 0; r < rows.size(); ++r)
        {
            auto id = rows[r];
            out << id + 1 << ',' << m.x[id].x << ',' << m.x[id].y << ',' << m.x[id].z << ',' << x[r] << '\n';
        }
        out.close();
        check(bool(out), "owned field output");
    }
    if (rankId == 0)
    {
        std::printf("PASS: stationary non-conformal mesh scalar solve ranks=%d iterations=%d independent_relative=%.3e "
                    "interface_balance=%.3e boundary_balance=%.3e\n",
                    rankCount, result.iterations, relative, double(std::abs(global[4] + global[5])),
                    double(std::abs(global[2] + global[3])));
        if (analytic) std::printf("PASS: u=x/2 max_error=%.3e; analytic boundary and interface flux\n", max_error);
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankId);
    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
    try
    {
        std::string input, output;
        double tolerance = 1e-12;
        bool analytic    = false;
        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];
            if (arg == "--mesh" && i + 1 < argc)
                input = argv[++i];
            else if (arg == "--output" && i + 1 < argc)
                output = argv[++i];
            else if (arg == "--distance-tolerance" && i + 1 < argc)
            {
                std::string value = argv[++i];
                size_t used       = 0;
                tolerance         = std::stod(value, &used);
                if (used != value.size()) throw std::runtime_error("invalid distance tolerance");
            }
            else if (arg == "--check-linear-channel")
                analytic = true;
            else
                throw std::runtime_error("usage: --mesh prepared.mnc [--output prefix] [--distance-tolerance value] "
                                         "[--check-linear-channel]");
        }
        if (input.empty() || !std::isfinite(tolerance) || tolerance <= 0)
            throw std::runtime_error("mesh and positive finite distance tolerance required");
        auto mesh = read_mesh(input);
        solve_mesh(mesh, tolerance, output, analytic);
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "rank %d mesh scalar example: %s\n", rankId, error.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
    return 0;
}
