#include "fixture.hpp"
#include "overlap.hpp"
#include "backend/distributed/unstructured/fem/mars_nonconformal_donor_exchange.hpp"
#include "backend/distributed/unstructured/fem/mars_nonconformal_scalar_graph.hpp"
#include <map>
#include <numeric>

static int rankId, rankCount;
static void check(bool ok, const char* message)
{
    if (!ok)
    {
        std::fprintf(stderr, "rank %d FAIL: %s\n", rankId, message);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
static void close(double a, double b, const char* message)
{
    check(std::isfinite(a) && std::isfinite(b) && std::abs(a - b) < 2e-12 * (1 + std::abs(b)), message);
}
static std::vector<Stencil> build(const std::vector<Point>& p,
                                  const std::vector<Face>& faces,
                                  const std::vector<int>& offsets,
                                  const std::vector<int>& candidates,
                                  [[maybe_unused]] const Face* deviceFaces = nullptr)
{
    std::vector<Stencil> result(p.size());
    InterfacePairingOptions<double> opts{1e-12};
#if defined(__CUDACC__) || defined(__HIPCC__)
    thrust::device_vector<Point> dp(p);
    thrust::device_vector<Face> df;
    if (!deviceFaces)
    {
        df          = faces;
        deviceFaces = thrust::raw_pointer_cast(df.data());
    }
    thrust::device_vector<int> d_off(offsets), d_ids(candidates);
    thrust::device_vector<Stencil> d_result(p.size());
    check(buildNonconformalScalarStencils(int(p.size()), thrust::raw_pointer_cast(dp.data()), deviceFaces,
                                          int(faces.size()), thrust::raw_pointer_cast(d_off.data()),
                                          thrust::raw_pointer_cast(d_ids.data()), int(candidates.size()), opts,
                                          thrust::raw_pointer_cast(d_result.data())) == cudaSuccess,
          "stencil launch");
    check(cudaDeviceSynchronize() == cudaSuccess, "stencil completion");
    thrust::copy(d_result.begin(), d_result.end(), result.begin());
#else
    check(buildNonconformalScalarStencilsHost(int(p.size()), p.data(), faces.data(), int(faces.size()), offsets.data(),
                                              candidates.data(), int(candidates.size()), opts, result.data()),
          "host stencils");
#endif
    for (const auto& s : result)
        check(s.status == InterfaceStatus::valid, "paired point");
    return result;
}
static void assemblyRejections()
{
    Stencil s;
    s.status    = InterfaceStatus::valid;
    s.row       = 7;
    int64_t row = 7, columns[8];
    double values[8];
    int offsets[]{0, 8}, po[]{0, 1}, pi[]{0};
    for (int j = 0; j < 8; ++j)
    {
        columns[j]  = 10 + j;
        s.column[j] = 10 + j;
        s.value[j]  = j - 3.5;
        values[j]   = 1.;
    }
    InterfaceCsrView<double> view{1, 8, &row, offsets, columns, values};
    check(assembleNonconformalScalarRowsHost(view, 1, &s, po, pi), "add to existing matrix");
    for (int j = 0; j < 8; ++j)
        close(values[j], 1. + s.value[j], "additive coefficient");
    check(assembleNonconformalScalarRowsHost(view, 1, &s, po, pi), "repeat assembly");
    for (int j = 0; j < 8; ++j)
        close(values[j], 1. + 2 * s.value[j], "repeat additive coefficient");
    auto bad   = s;
    bad.status = InterfaceStatus::no_donor;
    check(!assembleNonconformalScalarRowsHost(view, 1, &bad, po, pi), "unpaired stencil rejected");
    bad     = s;
    bad.row = 8;
    check(!assembleNonconformalScalarRowsHost(view, 1, &bad, po, pi), "foreign row rejected");
    bad          = s;
    bad.value[3] = std::numeric_limits<double>::quiet_NaN();
    check(!assembleNonconformalScalarRowsHost(view, 1, &bad, po, pi), "nonfinite coefficient rejected");
    pi[0] = -1;
    check(!assembleNonconformalScalarRowsHost(view, 1, &s, po, pi), "invalid point index rejected");
    pi[0]      = 0;
    columns[1] = columns[0];
    check(!assembleNonconformalScalarRowsHost(view, 1, &s, po, pi), "duplicate CSR column rejected");
    columns[1] = 11;
    offsets[1] = 9;
    check(!assembleNonconformalScalarRowsHost(view, 1, &s, po, pi), "invalid CSR span rejected");
    offsets[1] = 8;
    for (int j = 0; j < 8; ++j)
        close(values[j], 1. + 2 * s.value[j], "rejections preserve row");
}

static void graphChecks()
{
    std::vector<int64_t> rows{2, 7, 9}, columns{2, 99, 9};
    std::vector<int> offsets{0, 2, 2, 3};
    std::vector<double> values{3, 4, 5};
    std::vector<Stencil> stencils(3);
    for (int p = 0; p < 3; ++p)
    {
        stencils[p].status = InterfaceStatus::valid;
        stencils[p].row    = p == 1 ? 2 : 7;
        int64_t cols[]{2, 7, 9, 99, 2, 7, 9, p == 1 ? 99 : INT64_C(1000000000000)};
        for (int j = 0; j < 8; ++j)
        {
            stencils[p].column[j] = cols[j];
            stencils[p].value[j]  = j - 3.5;
        }
    }
    NonconformalScalarGraphHost<double> host;
#if defined(__CUDACC__) || defined(__HIPCC__)
    NonconformalScalarGraph<double> device;
    cudaStream_t stream;
    check(cudaStreamCreate(&stream) == cudaSuccess, "graph setup stream");
#endif
    auto exercise = [&](bool expected)
    {
        InterfaceCsrView<double> base{int(rows.size()), int(columns.size()), rows.data(),
                                      offsets.data(),   columns.data(),      values.data()};
        auto old = host;
        check(buildNonconformalScalarGraphHost(base, int(stencils.size()), stencils.data(), host) == expected,
              "graph validation result");
        if (!expected)
            check(host.rowIds == old.rowIds && host.columns == old.columns && host.offsets == old.offsets &&
                      host.pointOffsets == old.pointOffsets && host.pointIds == old.pointIds &&
                      host.values == old.values,
                  "rejected graph preserves previous output");
#if defined(__CUDACC__) || defined(__HIPCC__)
        thrust::device_vector<int64_t> dr(rows), dc(columns);
        thrust::device_vector<int> d_offsets(offsets);
        thrust::device_vector<double> dv(values);
        thrust::device_vector<Stencil> ds(stencils);
        InterfaceCsrView<double> db{int(rows.size()),
                                    int(columns.size()),
                                    thrust::raw_pointer_cast(dr.data()),
                                    thrust::raw_pointer_cast(d_offsets.data()),
                                    thrust::raw_pointer_cast(dc.data()),
                                    thrust::raw_pointer_cast(dv.data())};
        check(buildNonconformalScalarGraph(db, int(stencils.size()), thrust::raw_pointer_cast(ds.data()), device,
                                           stream) == expected,
              "device graph validation result");
        auto equal = [](const auto& d, const auto& h)
        {
            using T                    = typename std::decay_t<decltype(h)>::value_type;
            thrust::host_vector<T> got = d;
            return std::vector<T>(got.begin(), got.end()) == h;
        };
        check(equal(device.rowIds, host.rowIds) && equal(device.columns, host.columns) &&
                  equal(device.offsets, host.offsets) && equal(device.pointOffsets, host.pointOffsets) &&
                  equal(device.pointIds, host.pointIds) && equal(device.values, host.values),
              "device graph and preservation match host");
#endif
    };
    exercise(true);
    check(host.columns == std::vector<int64_t>({2, 7, 9, 99, 2, 7, 9, 99, INT64_C(1000000000000), 9}) &&
              host.offsets == std::vector<int>({0, 4, 9, 10}) && host.pointOffsets == std::vector<int>({0, 1, 3, 3}) &&
              host.pointIds == std::vector<int>({1, 0, 2}) &&
              host.values == std::vector<double>({3, 0, 0, 4, 0, 0, 0, 0, 0, 5}),
          "graph union, duplicate removal, int64 columns and base values");
    auto assembled = host;
    auto expected  = host.values;
    for (const auto& stencil : stencils)
    {
        auto row = std::find(host.rowIds.begin(), host.rowIds.end(), stencil.row) - host.rowIds.begin();
        for (int j = 0; j < 8; ++j)
            for (int k = host.offsets[row]; k < host.offsets[row + 1]; ++k)
                if (host.columns[k] == stencil.column[j]) expected[k] += stencil.value[j];
    }
    check(assembleNonconformalScalarRowsHost(assembled.view(), int(stencils.size()), stencils.data(),
                                             assembled.pointOffsets.data(), assembled.pointIds.data()),
          "assemble into generated graph");
    check(assembled.values == expected, "volume and interface entries retained together");
#if defined(__CUDACC__) || defined(__HIPCC__)
    auto assembledDevice = device;
    thrust::device_vector<Stencil> ds(stencils);
    thrust::device_vector<int> valid(rows.size());
    check(assembleNonconformalScalarRows(assembledDevice.view(), int(stencils.size()),
                                         thrust::raw_pointer_cast(ds.data()),
                                         thrust::raw_pointer_cast(assembledDevice.pointOffsets.data()),
                                         thrust::raw_pointer_cast(assembledDevice.pointIds.data()),
                                         thrust::raw_pointer_cast(valid.data()), stream) == cudaSuccess,
          "device assembly into generated graph");
    check(cudaStreamSynchronize(stream) == cudaSuccess, "generated graph assembly completion");
    thrust::host_vector<double> assembledValues = assembledDevice.values;
    thrust::host_vector<int> accepted           = valid;
    for (int v : accepted)
        check(v == 1, "generated graph row accepted");
    check(std::vector<double>(assembledValues.begin(), assembledValues.end()) == expected,
          "device volume and interface entries retained together");
#endif
    stencils[0].row = 8;
    exercise(false);
    stencils[0].row       = 7;
    stencils[0].column[0] = -1;
    exercise(false);
    stencils[0].column[0] = 2;
    stencils[0].status    = InterfaceStatus::no_donor;
    exercise(false);
    stencils[0].status   = InterfaceStatus::valid;
    stencils[0].value[0] = std::numeric_limits<double>::quiet_NaN();
    exercise(false);
    stencils[0].value[0] = -3.5;
    rows[1]              = 2;
    exercise(false);
    rows[1]    = 7;
    offsets[1] = 4;
    exercise(false);
    offsets[1] = 2;
    offsets[0] = 1;
    exercise(false);
    offsets[0] = 0;
    columns[1] = 2;
    exercise(false);
    columns[1] = 99;
    values[0]  = std::numeric_limits<double>::infinity();
    exercise(false);
    values[0] = 3;
    stencils.clear();
    exercise(true);
    check(host.columns == std::vector<int64_t>({2, 99, 7, 9}) && host.values == std::vector<double>({3, 4, 0, 5}),
          "empty interface preserves base and adds missing diagonal");
    rows.clear();
    columns.clear();
    values.clear();
    offsets = {0};
    exercise(true);
    check(host.offsets == std::vector<int>({0}) && host.pointOffsets == std::vector<int>({0}) && host.columns.empty() &&
              host.pointIds.empty(),
          "empty graph on idle rank");
    InterfaceCsrView<double> invalid{INT_MAX, 0, nullptr, nullptr, nullptr, nullptr};
    check(!buildNonconformalScalarGraphHost(invalid, 0, static_cast<const Stencil*>(nullptr), host),
          "oversized graph rejected before pointer access");
#if defined(__CUDACC__) || defined(__HIPCC__)
    check(!buildNonconformalScalarGraph(invalid, 0, static_cast<const Stencil*>(nullptr), device),
          "oversized device graph rejected before pointer access");
    check(cudaStreamDestroy(stream) == cudaSuccess, "graph stream destroy");
#endif
}

static void conservativeDistributed(const std::vector<Face>& faces, int active, int round)
{
    std::vector<InterfaceFacePair> pairs;
    if (!faces.empty())
    {
        int coarse =
            int(std::find_if(faces.begin(), faces.end(), [](const Face& f) { return f.id == 10; }) - faces.begin());
        for (int i = 0; i < int(faces.size()); ++i)
            if (faces[i].side == 1) pairs.push_back({coarse, i});
        std::sort(pairs.begin(), pairs.end(), [&](auto a, auto b) { return faces[a.second].id < faces[b.second].id; });
    }
    auto all = conservativeFixture(faces, pairs, 1e-12,
#if defined(__CUDACC__) || defined(__HIPCC__)
                                   true
#else
                                   false
#endif
    );
    std::vector<Stencil> stencils;
    for (auto s : all)
        if (s.row % active == rankId) stencils.push_back(s);
    std::vector<int64_t> rows;
    if (rankId < active)
        for (int64_t i = rankId; i < 11; i += active)
            rows.push_back(i);
    std::vector<int> offsets(rows.size() + 1, 0);
    InterfaceCsrView<double> base{int(rows.size()), 0, rows.data(), offsets.data(), nullptr, nullptr};
    NonconformalScalarGraphHost<double> graph;
    check(buildNonconformalScalarGraphHost(base, int(stencils.size()), stencils.data(), graph),
          "conservative host graph");
    auto host = InterfaceCsrView<double>{int(rows.size()),     int(graph.columns.size()), rows.data(),
                                         graph.offsets.data(), graph.columns.data(),      graph.values.data()};
    check(assembleNonconformalScalarRowsHost(host, int(stencils.size()), stencils.data(), graph.pointOffsets.data(),
                                             graph.pointIds.data()),
          "conservative host assembly");
#if defined(__CUDACC__) || defined(__HIPCC__)
    thrust::device_vector<int64_t> dr(rows);
    thrust::device_vector<int> off(offsets), valid(rows.size());
    thrust::device_vector<Stencil> ds(stencils);
    InterfaceCsrView<double> db{
        int(rows.size()), 0,      thrust::raw_pointer_cast(dr.data()), thrust::raw_pointer_cast(off.data()),
        nullptr,          nullptr};
    NonconformalScalarGraph<double> dg;
    check(buildNonconformalScalarGraph(db, int(stencils.size()), thrust::raw_pointer_cast(ds.data()), dg),
          "conservative device graph");
    check(assembleNonconformalScalarRows(dg.view(), int(stencils.size()), thrust::raw_pointer_cast(ds.data()),
                                         thrust::raw_pointer_cast(dg.pointOffsets.data()),
                                         thrust::raw_pointer_cast(dg.pointIds.data()),
                                         thrust::raw_pointer_cast(valid.data())) == cudaSuccess,
          "conservative device assembly");
    check(cudaDeviceSynchronize() == cudaSuccess, "conservative assembly completion");
    thrust::host_vector<int> ok = valid, ho = dg.offsets;
    thrust::host_vector<int64_t> hc = dg.columns;
    thrust::host_vector<double> hv  = dg.values;
    for (int v : ok)
        check(v == 1, "conservative valid row");
    check(std::vector<int>(ho.begin(), ho.end()) == graph.offsets &&
              std::vector<int64_t>(hc.begin(), hc.end()) == graph.columns,
          "conservative graph parity");
    check(hv.size() == graph.values.size(), "conservative coefficient count");
    for (size_t i = 0; i < hv.size(); ++i)
        close(hv[i], graph.values[i], "conservative coefficient parity");
    graph.values.assign(hv.begin(), hv.end());
#endif
    double column[11]{}, global[11]{}, field[11]{}, flux[2]{}, total[2]{};
    for (const auto& f : faces)
        for (int j = 0; j < 4; ++j)
            field[f.dof[j]] = affine(f.x[j]);
    for (size_t i = 0; i < rows.size(); ++i)
    {
        double null = 0;
        for (int j = graph.offsets[i]; j < graph.offsets[i + 1]; ++j)
        {
            auto col = graph.columns[j];
            column[col] += graph.values[j];
            null += graph.values[j];
            flux[rows[i] < 4 ? 0 : 1] += graph.values[j] * field[col];
        }
        close(null, 0, "conservative row constant nullspace");
    }
    MPI_Allreduce(column, global, 11, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (double v : global)
        close(v, 0, "conservative assembled basis balance");
    MPI_Allreduce(flux, total, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    close(total[0], -2 * (2. + round), "conservative affine coarse flux");
    close(total[1], 2 * (2. + round), "conservative affine fine flux");
    if (rankId == 0)
        std::printf("PASS ranks=%d round=%d: shared overlap, owned CSR, all-basis flux balance\n", rankCount,
                    round + 1);
}

static void run()
{
    graphChecks();
    assemblyRejections();
    int active         = std::min(rankCount, 3);
    auto expectedFaces = fixture();
    std::vector<int> peers;
    if (rankId < active)
        for (int i = 0; i < active; ++i)
            if (i != rankId) peers.push_back(i);
    std::vector<uint64_t> ids;
    std::vector<int> owners, order;
    if (rankId < active)
    {
        for (int i = 0; i < 5; ++i)
            order.push_back((i + rankId) % 5);
        if (rankId % 2) std::reverse(order.begin(), order.end());
        for (int i : order)
        {
            ids.push_back(expectedFaces[i].id);
            owners.push_back(i % active);
        }
    }
    NonconformalDonorExchange<double> exchange(MPI_COMM_WORLD);
    check(exchange.buildFromHost(ids, owners, peers), "build donor plan");
    for (int round = 0; round < 2; ++round)
    {
        for (auto& f : expectedFaces)
            f.diffusivity = 2.0 + round;
        std::vector<Face> faces(ids.size());
        for (size_t i = 0; i < ids.size(); ++i)
            if (owners[i] == rankId) faces[i] = expectedFaces[order[i]];
#if defined(__CUDACC__) || defined(__HIPCC__)
        thrust::device_vector<Face> d_faces(faces);
        check(exchange.forward(thrust::raw_pointer_cast(d_faces.data()), d_faces.size()), "device donor exchange");
        thrust::copy(d_faces.begin(), d_faces.end(), faces.begin());
#else
        check(exchange.forwardHost(faces), "host donor exchange");
#endif
        for (size_t i = 0; i < faces.size(); ++i)
        {
            const auto& expected = expectedFaces[order[i]];
            check(faces[i].id == expected.id && faces[i].side == expected.side, "face identity");
            close(faces[i].diffusivity, expected.diffusivity, "donor numeric refresh");
            for (int j = 0; j < 4; ++j)
            {
                check(faces[i].dof[j] == expected.dof[j], "donor global DOF");
                close(faces[i].x[j].x, expected.x[j].x, "donor x");
                close(faces[i].x[j].y, expected.x[j].y, "donor y");
                close(faces[i].x[j].z, expected.x[j].z, "donor z");
            }
        }
        conservativeDistributed(faces, active, round);
        std::vector<Point> p;
        std::vector<int> offsets{0}, candidates;
        if (rankId < active)
            for (auto point : points(expectedFaces))
            {
                if (expectedFaces[point.face].dof[point.rowCorner] % active != rankId) continue;
                int oldFace = point.face;
                point.face  = int(std::find(order.begin(), order.end(), oldFace) - order.begin());
                p.push_back(point);
                for (size_t j = 0; j < faces.size(); ++j)
                    if (faces[j].side != faces[point.face].side) candidates.push_back(int(j));
                offsets.push_back(int(candidates.size()));
            }
#if defined(__CUDACC__) || defined(__HIPCC__)
        auto stencils = build(p, faces, offsets, candidates, thrust::raw_pointer_cast(d_faces.data()));
#else
        auto stencils = build(p, faces, offsets, candidates);
#endif
        std::vector<Stencil> oracle(p.size());
        std::vector<Face> oracleFaces;
        for (int i : order)
            oracleFaces.push_back(expectedFaces[i]);
        check(buildNonconformalScalarStencilsHost(int(p.size()), p.data(), oracleFaces.data(), int(oracleFaces.size()),
                                                  offsets.data(), candidates.data(), int(candidates.size()),
                                                  InterfacePairingOptions<double>{1e-12}, oracle.data()),
              "single-rank oracle");
        std::vector<int64_t> rows, columns;
        std::vector<int> rowOffsets{0}, pointOffsets{0}, pointIds;
        std::vector<double> expected, values;
        if (rankId < active)
            for (int64_t row = rankId; row < 11; row += active)
            {
                rows.push_back(row);
                std::map<int64_t, double> entry{{row, 0.0}};
                for (size_t i = 0; i < oracle.size(); ++i)
                    if (oracle[i].row == row)
                    {
                        pointIds.push_back(int(i));
                        for (int j = 0; j < 8; ++j)
                            entry[oracle[i].column[j]] += oracle[i].value[j];
                    }
                for (const auto& e : entry)
                {
                    columns.push_back(e.first);
                    expected.push_back(e.second);
                }
                rowOffsets.push_back(int(columns.size()));
                pointOffsets.push_back(int(pointIds.size()));
            }
        std::vector<int> baseOffsets(rows.size() + 1);
        std::iota(baseOffsets.begin(), baseOffsets.end(), 0);
        std::vector<double> baseValues(rows.size(), 0.0);
        InterfaceCsrView<double> base{int(rows.size()),   int(rows.size()), rows.data(),
                                      baseOffsets.data(), rows.data(),      baseValues.data()};
        NonconformalScalarGraphHost<double> graph;
        check(buildNonconformalScalarGraphHost(base, int(stencils.size()), stencils.data(), graph), "host graph build");
        check(graph.rowIds == rows && graph.columns == columns && graph.offsets == rowOffsets &&
                  graph.pointOffsets == pointOffsets && graph.pointIds == pointIds,
              "graph matches independent map");
        values = graph.values;
        InterfaceCsrView<double> view{int(rows.size()),  int(columns.size()), rows.data(),
                                      rowOffsets.data(), columns.data(),      values.data()};
#if defined(__CUDACC__) || defined(__HIPCC__)
        thrust::device_vector<int64_t> d_rows(rows);
        thrust::device_vector<int> d_baseOffsets(baseOffsets), d_valid(rows.size());
        thrust::device_vector<double> d_baseValues(baseValues);
        thrust::device_vector<Stencil> d_stencils(stencils);
        InterfaceCsrView<double> d_base{int(rows.size()),
                                        int(rows.size()),
                                        thrust::raw_pointer_cast(d_rows.data()),
                                        thrust::raw_pointer_cast(d_baseOffsets.data()),
                                        thrust::raw_pointer_cast(d_rows.data()),
                                        thrust::raw_pointer_cast(d_baseValues.data())};
        NonconformalScalarGraph<double> d_graph;
        check(buildNonconformalScalarGraph(d_base, int(stencils.size()), thrust::raw_pointer_cast(d_stencils.data()),
                                           d_graph),
              "device graph build");
        auto device                               = d_graph.view();
        auto& d_po                                = d_graph.pointOffsets;
        auto& d_pi                                = d_graph.pointIds;
        auto& d_values                            = d_graph.values;
        thrust::host_vector<int64_t> graphColumns = d_graph.columns;
        thrust::host_vector<int> graphOffsets = d_graph.offsets, graphPO = d_po, graphPI = d_pi;
        check(std::vector<int64_t>(graphColumns.begin(), graphColumns.end()) == columns &&
                  std::vector<int>(graphOffsets.begin(), graphOffsets.end()) == rowOffsets &&
                  std::vector<int>(graphPO.begin(), graphPO.end()) == pointOffsets &&
                  std::vector<int>(graphPI.begin(), graphPI.end()) == pointIds,
              "device graph matches independent map");
        check(assembleNonconformalScalarRows(device, int(stencils.size()), thrust::raw_pointer_cast(d_stencils.data()),
                                             thrust::raw_pointer_cast(d_po.data()),
                                             thrust::raw_pointer_cast(d_pi.data()),
                                             thrust::raw_pointer_cast(d_valid.data())) == cudaSuccess,
              "assembly launch");
        check(cudaDeviceSynchronize() == cudaSuccess, "assembly completion");
        thrust::host_vector<int> valid = d_valid;
        for (int v : valid)
            check(v == 1, "device CSR row accepted");
        thrust::copy(d_values.begin(), d_values.end(), values.begin());
#else
        check(assembleNonconformalScalarRowsHost(view, int(stencils.size()), stencils.data(), pointOffsets.data(),
                                                 pointIds.data()),
              "host CSR assembly");
#endif
        for (size_t k = 0; k < values.size(); ++k)
            close(values[k], expected[k], "owned CSR coefficient");
        double field[11]{};
        for (const auto& f : expectedFaces)
            for (int j = 0; j < 4; ++j)
                field[f.dof[j]] = affine(f.x[j]);
        double flux[2]{};
        for (size_t row = 0; row < rows.size(); ++row)
        {
            double constant = 0, action = 0;
            for (int k = rowOffsets[row]; k < rowOffsets[row + 1]; ++k)
            {
                constant += values[k];
                action += values[k] * field[columns[k]];
            }
            close(constant, 0, "assembled constant nullspace");
            flux[rows[row] < 4 ? 0 : 1] += action;
        }
        double global[2];
        MPI_Allreduce(flux, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        close(global[0], -2 * (2.0 + round), "coarse-side affine flux");
        close(global[1], 2 * (2.0 + round), "fine-side affine flux");
        if (rankId == 0)
            std::printf("PASS ranks=%d round=%d: donor refresh, generated CSR, constant and affine flux\n", rankCount,
                        round + 1);

        if (!stencils.empty())
        {
            auto invalid              = stencils;
            invalid.front().column[0] = 9999;
            auto sentinel             = values;
            check(!assembleNonconformalScalarRowsHost(view, int(invalid.size()), invalid.data(), pointOffsets.data(),
                                                      pointIds.data()),
                  "missing remote column rejected");
            // Other rows may have been assembled; only the rejected row must remain untouched.
            auto r = std::find(rows.begin(), rows.end(), invalid.front().row) - rows.begin();
            for (int k = rowOffsets[r]; k < rowOffsets[r + 1]; ++k)
                close(values[k], sentinel[k], "rejected row unchanged");
#if defined(__CUDACC__) || defined(__HIPCC__)
            d_stencils    = invalid;
            d_values      = sentinel;
            device.values = thrust::raw_pointer_cast(d_values.data());
            check(assembleNonconformalScalarRows(
                      device, int(invalid.size()), thrust::raw_pointer_cast(d_stencils.data()),
                      thrust::raw_pointer_cast(d_po.data()), thrust::raw_pointer_cast(d_pi.data()),
                      thrust::raw_pointer_cast(d_valid.data())) == cudaSuccess,
                  "invalid graph assembly launch");
            check(cudaDeviceSynchronize() == cudaSuccess, "invalid graph completion");
            valid = d_valid;
            check(valid[r] == 0, "device missing column rejected");
            thrust::copy(d_values.begin(), d_values.end(), values.begin());
            for (int k = rowOffsets[r]; k < rowOffsets[r + 1]; ++k)
                close(values[k], sentinel[k], "device rejected row unchanged");
#endif
        }
    }
    std::vector<Face> wrongSize(ids.size());
    if (rankId == 0) wrongSize.push_back({});
    check(!exchange.forwardHost(wrongSize), "short/long field rejected collectively");
#if defined(__CUDACC__) || defined(__HIPCC__)
    thrust::device_vector<Face> wrongDevice(wrongSize);
    check(!exchange.forward(thrust::raw_pointer_cast(wrongDevice.data()), wrongDevice.size()),
          "device capacity rejected collectively");
#endif
    auto badOwners = owners;
    if (rankId == 0) badOwners[0] = rankCount;
    check(!exchange.buildFromHost(ids, badOwners, peers), "invalid owner rejected collectively");
    auto badIds = ids;
    if (rankId == 0) badIds[1] = badIds[0];
    check(!exchange.buildFromHost(badIds, owners, peers), "duplicate identity rejected collectively");
    if (active > 1)
    {
        auto asymmetric = peers;
        // Keep valid owner lists out of this check so the reciprocal-peer validation is exercised.
        if (rankId == 0) asymmetric.clear();
        check(!exchange.buildFromHost({}, {}, asymmetric), "one-sided peers rejected collectively");
        badIds = ids;
        if (rankId == 0)
            for (size_t i = 0; i < owners.size(); ++i)
                if (owners[i] != rankId)
                {
                    badIds[i] = 999;
                    break;
                }
        check(!exchange.buildFromHost(badIds, owners, peers), "unresolved donor rejected collectively");
    }
    check(exchange.buildFromHost({}, {}, {}), "empty plan rebuild");
    std::vector<Face> empty;
    check(exchange.forwardHost(empty), "empty exchange");
    if (rankId == 0) std::puts("PASS: stationary non-conformal distributed scalar checks");
}
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankId);
    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
    try
    {
        run();
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "rank %d: %s\n", rankId, error.what());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    MPI_Finalize();
}
