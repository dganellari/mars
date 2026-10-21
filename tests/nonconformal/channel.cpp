#include "examples/distributed/unstructured/mars_nonconformal_scalar_runtime.hpp"
#include "channel_mesh.hpp"
struct VolumeRows
{
    const Vec3<double>* x;
    const int* cells;
    const double* k;
    const int* offsets;
    const int* incident;
    InterfaceCsrView<double> matrix;
    MARS_FS_FN void operator()(int r) const
    {
        int64_t row = matrix.rowIds[r];
        for (int entry = offsets[r]; entry < offsets[r + 1]; ++entry)
        {
            int e = incident[entry], a = -1;
            InterfaceTetFace<double> f{};
            f.diffusivity = k[e];
            f.opposite    = 3;
            for (int j = 0; j < 4; ++j)
            {
                int id   = cells[4 * e + j];
                f.x[j]   = x[id];
                f.dof[j] = id;
                if (id == row) a = j;
            }
            nonconformal_detail::Geometry<double> g;
            nonconformal_detail::geometry(f, g);
            double volume = std::abs(dot(sub(f.x[1], f.x[0]),
                                         nonconformal_detail::cross(sub(f.x[2], f.x[0]), sub(f.x[3], f.x[0])))) /
                            6;
            for (int j = 0; j < 4; ++j)
            {
                int slot =
                    nonconformal_detail::columnSlot(matrix.columns, matrix.offsets[r], matrix.offsets[r + 1], f.dof[j]);
                matrix.values[slot] += k[e] * volume * dot(g.gradient[a], g.gradient[j]);
            }
        }
    }
};
static Graph assemble(const Mesh& m, const std::vector<int64_t>& rows, const Buffer<InterfaceStencil<double>>& stencils)
{
    std::vector<int> offsets{0}, incOffsets{0}, incident;
    std::vector<int64_t> columns;
    for (auto row : rows)
    {
        std::set<int64_t> col{row};
        for (int e = 0; e < int(m.elements.size()); ++e)
        {
            bool has = false;
            for (int j = 0; j < 4; ++j)
                has = has || m.cells[4 * e + j] == row;
            if (has)
            {
                incident.push_back(e);
                for (int j = 0; j < 4; ++j)
                    col.insert(m.cells[4 * e + j]);
            }
        }
        columns.insert(columns.end(), col.begin(), col.end());
        offsets.push_back(int(columns.size()));
        incOffsets.push_back(int(incident.size()));
    }
    Graph base, out;
    base.rowIds  = Buffer<int64_t>(rows.begin(), rows.end());
    base.columns = Buffer<int64_t>(columns.begin(), columns.end());
    base.offsets = Buffer<int>(offsets.begin(), offsets.end());
    base.values.resize(columns.size(), 0.);
    Buffer<Vec3<double>> coords(m.x.begin(), m.x.end());
    Buffer<int> cells(m.cells.begin(), m.cells.end()), io(incOffsets.begin(), incOffsets.end()),
        ie(incident.begin(), incident.end());
    Buffer<double> k(m.k.begin(), m.k.end());
    each<device>(int(rows.size()),
                 VolumeRows{pointer(coords), pointer(cells), pointer(k), pointer(io), pointer(ie), base.view()});
    completed();
#if defined(__CUDACC__) || defined(__HIPCC__)
    check(buildNonconformalScalarGraph(base.view(), int(stencils.size()), pointer(stencils), out),
          "interface graph union");
    Buffer<int> row_valid(rows.size(), 0);
    check(assembleNonconformalScalarRows(out.view(), int(stencils.size()), pointer(stencils), pointer(out.pointOffsets),
                                         pointer(out.pointIds), pointer(row_valid)) == cudaSuccess,
          "interface row assembly launch");
    completed();
    // A successful launch does not mean every row accepted its stencil list.
    check(thrust::count(row_valid.begin(), row_valid.end(), 1) == int(rows.size()),
          "interface row assembly accepted on all ranks");
#else
    check(buildNonconformalScalarGraphHost(base.view(), int(stencils.size()), pointer(stencils), out),
          "interface graph union");
    check(assembleNonconformalScalarRowsHost(out.view(), int(stencils.size()), pointer(stencils),
                                             pointer(out.pointOffsets), pointer(out.pointIds)),
          "interface row assembly");
#endif
    completed();
    return out;
}
static void rejections()
{
    auto m    = mesh(1, 1.);
    auto rows = ownedRows(m);
    NonconformalInterfacePlan<double, device> plan(MPI_COMM_WORLD);
    auto trial = [&](const Mesh& test, bool expected, const char* name)
    {
        auto f = extract(test);
        check(plan.build(f.faces, f.ids, f.owners, rows, peers(), 1e-12) == expected, name);
        if (!expected) check(plan.stencils().empty(), "rejected plan clears old rows");
    };
    trial(m, true, "complete interface accepted");
    auto hole    = m;
    auto remove  = std::find_if(hole.tags.begin(), hole.tags.end(), [](auto t) { return t.side == 1; });
    auto removed = *remove;
    hole.tags.erase(remove);
    trial(hole, false, "missing coverage rejected");
    auto duplicate = [&](Mesh& a, InterfaceFaceTag source)
    {
        int e = int(a.elements.size());
        a.elements.push_back(e);
        a.owners.push_back(e % rankCount);
        a.k.push_back(a.k[source.element]);
        for (int j = 0; j < 4; ++j)
            a.cells.push_back(a.cells[4 * source.element + j]);
        source.element = e;
        a.tags.push_back(source);
    };
    auto repeated = m;
    duplicate(repeated, removed);
    trial(repeated, false, "duplicate coverage rejected");
    auto compensated = hole;
    auto other       = *std::find_if(hole.tags.begin(), hole.tags.end(), [](auto t) { return t.side == 1; });
    duplicate(compensated, other);
    trial(compensated, false, "equal-area hole and duplicate rejected");
    auto gap = m;
    for (size_t i = 0; i < gap.x.size(); ++i)
        if (gap.side[i]) gap.x[i].x += 1e-5;
    trial(gap, false, "separated interface rejected");
    trial(m, true, "valid rebuild after rejection");
    auto f = extract(m);
    check(!plan.build(f.faces, f.ids, f.owners, rows, peers(), rankId == 0 ? -1. : 1e-12),
          "rank-local bad tolerance rejected collectively");
    auto idle = m;
    std::fill(idle.owners.begin(), idle.owners.end(), 0);
    auto fi = extract(idle);
    if (rankId)
    {
        fi.faces.clear();
        fi.ids.clear();
        fi.owners.clear();
    }
    check(plan.build(fi.faces, fi.ids, fi.owners, rankId == 0 ? m.dofs : std::vector<int64_t>{}, peers(), 1e-12),
          "empty interface ranks participate");
    if (rankId == 0)
        check(!plan.stencils().empty(), "active interface has rows");
    else
        check(plan.stencils().empty(), "idle interface has no rows");
    if (rankId == 0) std::puts("PASS: search rejects holes, duplicate coverage, compensated holes and separated faces");
}
static void run(int n, double rightK)
{
    auto m    = mesh(n, rightK);
    auto rows = ownedRows(m);
    auto f    = extract(m);
    NonconformalInterfacePlan<double, device> plan(MPI_COMM_WORLD);
    check(plan.build(f.faces, f.ids, f.owners, rows, peers(), 1e-12), "automatic interface search and coverage");
    auto a = assemble(m, rows, plan.stencils());
    // Preserve unconstrained rows for the independent boundary-reaction balance.
    auto rawOffsets       = host<int>(a.offsets);
    auto rawColumns       = host<int64_t>(a.columns);
    auto rawValues        = host<double>(a.values);
    constexpr int restart = 60;
    ChannelOps ops(a, m, restart);
    auto result = outlet_fgmres(ops, restart, 2000, 1e-11, std::numeric_limits<double>::epsilon());
    completed();
    check(result.converged, "channel FGMRES converged");
    // Field downloads below are test diagnostics, never part of interface setup or the solver.
    auto localX = host<double>(ops.x);
    std::vector<double> partial(m.x.size(), 0.), solution(m.x.size());
    for (size_t i = 0; i < rows.size(); ++i)
        partial[rows[i]] = localX[i];
    MPI_Allreduce(partial.data(), solution.data(), int(solution.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double q = 1. / (1. + 1. / rightK), error = 0., reactionLeft = 0., reactionRight = 0., balance = 0.;
    long double residualSquared = 0., rhsSquared = 0.;
    for (size_t r = 0; r < rows.size(); ++r)
    {
        auto id  = rows[r];
        double x = m.x[id].x, exact = m.side[id] ? q + q * x / rightK : q * (x + 1.);
        error           = std::max(error, std::abs(solution[id] - exact));
        double residual = 0;
        for (int j = rawOffsets[r]; j < rawOffsets[r + 1]; ++j)
            residual += rawValues[j] * solution[rawColumns[j]];
        if (x == -1.) reactionLeft += residual;
        if (x == 1.) reactionRight += residual;
        balance += residual;
        long double independent = 0.;
        for (int j = rawOffsets[r]; j < rawOffsets[r + 1]; ++j)
            independent += static_cast<long double>(rawValues[j]) * solution[rawColumns[j]];
        if (std::abs(x) == 1.)
        {
            double rhs  = x > 0 ? 1. : 0.;
            independent = solution[id] - rhs;
            rhsSquared += rhs * rhs;
        }
        residualSquared += independent * independent;
    }
    double flux[2]{0., 0.};
    auto hs = host<InterfaceStencil<double>>(plan.stencils());
    for (const auto& s : hs)
    {
        double value = 0.;
        for (int j = 0; j < 8; ++j)
            value += s.value[j] * solution[s.column[j]];
        flux[m.side[s.row]] += value;
    }
    double local[]{error, reactionLeft, reactionRight, balance, flux[0], flux[1]}, global[6]{};
    MPI_Allreduce(local, global, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double maxError = 0.;
    MPI_Allreduce(&error, &maxError, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    long double sums[]{residualSquared, rhsSquared}, totals[2]{};
    MPI_Allreduce(sums, totals, 2, MPI_LONG_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double independentRelative = std::sqrt(totals[0] / totals[1]);
    check(std::isfinite(independentRelative) && independentRelative < 2e-11, "independent constrained residual");
    check(maxError < 2e-8, "piecewise affine channel solution");
    check(std::abs(global[1] + q) < 2e-8 && std::abs(global[2] - q) < 2e-8, "analytic outer boundary reactions");
    check(std::abs(global[3]) < 2e-8 && std::abs(global[4] + global[5]) < 2e-12, "global and interface flux balance");
    check(std::abs(global[4] + q) < 2e-8 && std::abs(global[5] - q) < 2e-8, "analytic transmitted interface flux");
    if (rankId == 0)
        std::printf("PASS ranks=%d n=%d k_right=%.0f: channel iterations=%d relative=%.3e max_error=%.3e "
                    "interface_balance=%.3e\n",
                    rankCount, n, rightK, result.iterations, result.relative_residual, maxError,
                    std::abs(global[4] + global[5]));
}
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankId);
    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
    // The launch binding script selects the visible GPU; do not override it here.
    {
        rejections();
        for (int n : {1, 2})
            for (double k : {1., 4.})
                run(n, k);
    }
    if (rankId == 0) std::puts("PASS: stationary non-conformal scalar channel checks");
    MPI_Finalize();
    return 0;
}
