#include "hypre_host_shim.hpp"
// Inspect resource identity without exposing production test hooks.
#define private public
#include "host_gmres.hpp"
#undef private

using Solver = mars::fem::HypreGMRESSolver<double, int, mars::HostTestTag>;
using Matrix = Solver::Matrix;

void check(bool okay, const char* message) {
    if (!okay) throw std::runtime_error(message);
}

template<class Function> void expect_failure(Function function) {
    bool rejected = false;
    try { function(); }
    catch (const std::runtime_error&) { rejected = true; }
    check(rejected, "invalid input was accepted");
}

void make_graph(Matrix& matrix, int n) {
    matrix.column_count = n + 1;
    matrix.offsets = {0};
    matrix.columns.clear();
    for (int row = 0; row < n; ++row) {
        // Reversed local numbering tests the supplied map. Invalid entries must
        // stay excluded regardless of their large numeric values.
        matrix.columns.push_back(n - 1 - row);
        if (row > 0) matrix.columns.push_back(n - row);
        if (row + 1 < n) matrix.columns.push_back(n - 2 - row);
        matrix.columns.push_back(-1);
        matrix.columns.push_back(n);
        matrix.columns.push_back(n + 7);
        matrix.offsets.push_back(matrix.columns.size());
    }
    matrix.values.resize(matrix.columns.size());
}

void fill_system(Matrix& matrix, const std::vector<HYPRE_BigInt>& map, int epoch,
                 std::vector<double>& rhs, std::vector<double>& truth) {
    const int n = matrix.numRows();
    rhs.assign(n, 0);
    truth.resize(n);
    for (int i = 0; i < n; ++i) truth[i] = 2 + std::sin(0.11 * i + 0.3 * epoch);
    for (int row = 0; row < n; ++row) {
        for (int slot = matrix.offsets[row]; slot < matrix.offsets[row + 1]; ++slot) {
            const int local = matrix.columns[slot];
            const int col = local >= 0 && local < int(map.size()) ? map[local] : -1;
            double value = 1e12;
            if (col >= 0) {
                value = col == row ? 4.0 + 0.01 * row + 0.8 * epoch
                    : (epoch == 1 || (epoch == 3 && row % 3 == 0)) ? 0.0
                    : col < row ? -0.5 - 0.1 * epoch : -1.0;
                rhs[row] += value * truth[col];
            }
            matrix.values[slot] = value;
        }
    }
}

void verify(const Matrix& matrix, const std::vector<HYPRE_BigInt>& map,
            const std::vector<double>& b, const std::vector<double>& x,
            const std::vector<double>& truth) {
    double residual2 = 0, rhs2 = 0, error = 0;
    for (int row = 0; row < matrix.numRows(); ++row) {
        double residual = b[row];
        for (int slot = matrix.offsets[row]; slot < matrix.offsets[row + 1]; ++slot) {
            const int local = matrix.columns[slot];
            if (local >= 0 && local < int(map.size()) && map[local] >= 0)
                residual -= matrix.values[slot] * x[map[local]];
        }
        residual2 += residual * residual;
        rhs2 += b[row] * b[row];
        error = std::max(error, std::abs(x[row] - truth[row]));
    }
    check(std::sqrt(residual2 / rhs2) < 2e-10, "independent true residual failed");
    check(error < 2e-9, "manufactured solution error failed");
}

void scaled_systems() {
    Matrix matrix;
    make_graph(matrix, 32);
    std::vector<HYPRE_BigInt> map(33);
    for (int i=0;i<32;++i) map[i]=31-i;
    map.back()=-1;
    for (bool cached : {false,true}) {
        Solver solver(0,300,1e-10);
        solver.setVerbose(false);
        solver.setAMGCoarseRelaxType(18);
        solver.enable_true_residual_check();
        if (cached) solver.enable_fixed_graph_updates();
        for (double scale : {1e-9,1.0,1e15}) {
            std::vector<double> b,truth,x(32,0);
            fill_system(matrix,map,0,b,truth);
            for (double& value:matrix.values) value*=scale;
            for (double& value:b) value*=scale;
            check(solver.solve(matrix,b,x,0,32,0,32,map), "scaled system rejected");
            verify(matrix,map,b,x,truth);
            check(solver.getLastFinalResidual()<1e-10, "explicit residual not reported");
        }
        std::vector<double> b(32,0),truth,x(32,0);
        check(solver.solve(matrix,b,x,0,32,0,32,map), "explicit zero RHS rejected");
        check(solver.getLastFinalResidual()==0, "explicit zero residual was stale");
        fill_system(matrix,map,0,b,truth);
        Solver stalled(0,1,1e-14,Solver::JACOBI,1);
        stalled.setVerbose(false);
        stalled.enable_true_residual_check();
        if (cached) stalled.enable_fixed_graph_updates();
        check(!stalled.solve(matrix,b,x,0,32,0,32,map), "unconverged solve accepted");
        check(stalled.getLastFinalResidual()>1e-14, "failed solve hid its true residual");
    }
}

int main() {
    setenv("MARS_HYPRE_MINITER", "0", 1);
    unsetenv("MARS_HYPRE_FLEXGMRES");
    try {
        scaled_systems();
        Matrix matrix;
        make_graph(matrix, 160);
        std::vector<HYPRE_BigInt> map(161);
        for (int i = 0; i < 160; ++i) map[i] = 159 - i;
        map.back() = -1;
        std::vector<double> b, truth, x(160, 0), baseline_x(160, 0);
        Solver persistent(0, 300, 1e-10), baseline(0, 300, 1e-10);
        persistent.setVerbose(false);
        baseline.setVerbose(false);
        baseline.enable_timing();
        persistent.enable_fixed_graph_updates();
        persistent.enable_timing();
        void *ij = nullptr, *parcsr = nullptr, *solver = nullptr, *amg = nullptr;
        void *rhs = nullptr, *solution = nullptr;
        const void *packed = nullptr, *slots = nullptr;
        for (int epoch = 0; epoch < 8; ++epoch) {
            fill_system(matrix, map, epoch, b, truth);
            std::fill(x.begin(), x.end(), 0.1 * epoch);
            std::fill(baseline_x.begin(), baseline_x.end(), 0.1 * epoch);
            const int barriers = host_barriers;
            const int reductions = host_reductions;
            check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "persistent solve failed");
            check(host_barriers == barriers, "fixed graph path added a barrier");
            if (epoch > 0) check(host_reductions - reductions == 9, "steady update collectives changed");
            verify(matrix, map, b, x, truth);
            check(baseline.solve(matrix, b, baseline_x, 0, 160, 0, 160, map), "baseline solve failed");
            check(host_barriers == barriers + 2, "legacy barrier behavior changed");
            verify(matrix, map, b, baseline_x, truth);
            check(baseline.get_last_timing().packing_seconds > 0, "baseline packing time absent");
            for (int i = 0; i < 160; ++i)
                check(std::abs(x[i] - baseline_x[i]) < 2e-9, "fresh/persistent field mismatch");
            if (epoch == 0) {
                ij = persistent.A_hypre_; parcsr = persistent.parcsr_A_;
                solver = persistent.solver_; amg = persistent.precond_;
                rhs = persistent.b_hypre_; solution = persistent.x_hypre_;
                packed = persistent.d_packed_values_.data(); slots = persistent.d_source_slots_.data();
            }
            check(ij == persistent.A_hypre_ && parcsr == persistent.parcsr_A_
                && solver == persistent.solver_ && amg == persistent.precond_
                && rhs == persistent.b_hypre_ && solution == persistent.x_hypre_
                && packed == persistent.d_packed_values_.data() && slots == persistent.d_source_slots_.data(),
                "fixed graph resources changed identity");
            check(persistent.get_graph_build_count() == 1, "graph unexpectedly rebuilt");
            check(persistent.get_numeric_update_count() == epoch + 1, "numeric refresh skipped");
            check(persistent.get_setup_count() == epoch + 1, "numeric AMG setup skipped");
            const auto& timing = persistent.get_last_timing();
            check(timing.prepare_seconds >= timing.packing_seconds && timing.packing_seconds > 0
                && timing.setup_seconds > 0 && timing.solve_seconds > 0 && timing.finish_seconds > 0,
                "timing scopes invalid");
        }
        check(baseline.get_graph_build_count() == 8, "baseline did not rebuild each time");

        // Value-buffer relocation alone does not invalidate a fixed graph.
        auto relocated_values = matrix.values;
        matrix.values.swap(relocated_values);
        fill_system(matrix, map, 8, b, truth);
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "relocated values failed");
        verify(matrix, map, b, x, truth);
        check(persistent.get_graph_build_count() == 1, "value relocation rebuilt graph");

        // Moving graph storage is detected even when the dimensions stay fixed.
        auto relocated_columns = matrix.columns;
        matrix.columns.swap(relocated_columns);
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "graph relocation failed");
        check(persistent.get_graph_build_count() == 2, "graph relocation was not detected");
        verify(matrix, map, b, x, truth);

        // In-place map edits require explicit invalidation under the API contract.
        std::swap(map[0], map[1]);
        persistent.invalidate_setup();
        fill_system(matrix, map, 9, b, truth);
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "invalidated map failed");
        verify(matrix, map, b, x, truth);
        check(persistent.get_graph_build_count() == 3, "explicit invalidation did not rebuild");
        persistent.setPointBlock(2);
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "configuration invalidation failed");
        check(persistent.get_graph_build_count() == 4, "configuration did not rebuild");
        verify(matrix, map, b, x, truth);

        // Frozen reuse keeps its old contract and setup count.
        persistent.setPointBlock(1);
        persistent.enable_reuse();
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "frozen setup failed");
        const int frozen_setups = persistent.get_setup_count();
        for (double& value : b) value *= 2;
        for (double& value : truth) value *= 2;
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "frozen solve failed");
        verify(matrix, map, b, x, truth);
        check(persistent.get_setup_count() == frozen_setups, "frozen mode rebuilt setup");
        persistent.enable_fixed_graph_updates();
        fill_system(matrix, map, 10, b, truth);
        check(persistent.solve(matrix, b, x, 0, 160, 0, 160, map), "mode transition failed");
        verify(matrix, map, b, x, truth);

        // New partitions and maps rebuild the graph and vector ranges together.
        const int graphs_before_resize = persistent.get_graph_build_count();
        make_graph(matrix, 192);
        map.resize(193);
        for (int i = 0; i < 192; ++i) map[i] = 191 - i;
        map.back() = -1;
        fill_system(matrix, map, 11, b, truth);
        x.assign(192, 0);
        check(persistent.solve(matrix, b, x, 0, 192, 0, 192, map), "partition resize failed");
        verify(matrix, map, b, x, truth);
        check(persistent.get_graph_build_count() == graphs_before_resize + 1,
              "partition resize did not rebuild");

        const double saved = matrix.values[0];
        matrix.values[0] = std::numeric_limits<double>::quiet_NaN();
        expect_failure([&] { persistent.solve(matrix, b, x, 0, 192, 0, 192, map); });
        matrix.values[0] = saved;
        const double saved_b = b[0];
        b[0] = std::numeric_limits<double>::infinity();
        expect_failure([&] { persistent.solve(matrix, b, x, 0, 192, 0, 192, map); });
        b[0] = saved_b;
        auto short_rhs = b;
        short_rhs.resize(2);
        expect_failure([&] { persistent.solve(matrix, short_rhs, x, 0, 192, 0, 192, map); });
        expect_failure([&] { persistent.solve<HYPRE_BigInt>(matrix, b, x, 0, 192, 0, 192, map); });
        persistent.setPrecondMatrix(&matrix);
        expect_failure([&] { persistent.solve(matrix, b, x, 0, 192, 0, 192, map); });
        persistent.setPrecondMatrix(nullptr);
        check(persistent.solve(matrix, b, x, 0, 192, 0, 192, map), "post-rejection rebuild failed");
        verify(matrix, map, b, x, truth);

        // The same fixed-graph path also supports the diagonal preconditioner.
        Solver jacobi(0, 300, 1e-10, Solver::JACOBI);
        jacobi.setVerbose(false);
        jacobi.enable_fixed_graph_updates();
        for (int epoch = 0; epoch < 3; ++epoch) {
            fill_system(matrix, map, epoch, b, truth);
            std::fill(x.begin(), x.end(), 0);
            check(jacobi.solve(matrix, b, x, 0, 192, 0, 192, map), "Jacobi refresh failed");
            verify(matrix, map, b, x, truth);
        }
        check(jacobi.get_setup_count() == 3 && jacobi.get_graph_build_count() == 1,
              "Jacobi lifecycle failed");
        b.assign(192, 0);
        truth.assign(192, 0);
        x.assign(192, 1);
        std::vector<double> fresh_zero_x(192, 1);
        Solver fresh_jacobi(0, 300, 1e-10, Solver::JACOBI);
        fresh_jacobi.setVerbose(false);
        const bool refreshed_zero = jacobi.solve(matrix, b, x, 0, 192, 0, 192, map);
        const bool fresh_zero = fresh_jacobi.solve(matrix, b, fresh_zero_x, 0, 192, 0, 192, map);
        check(refreshed_zero == fresh_zero, "zero RHS warm-guess acceptance changed");
        for (int i = 0; i < 192; ++i)
            check(std::abs(x[i] - fresh_zero_x[i]) < 2e-9, "zero RHS warm-guess parity failed");
        x.assign(192, 0);
        check(jacobi.solve(matrix, b, x, 0, 192, 0, 192, map), "zero RHS refresh failed");
        for (double value : x) check(std::abs(value) < 1e-9, "zero RHS retained stale solution");
        check(jacobi.getLastIterations() == 0 && jacobi.getLastFinalResidual() == 0,
              "zero-iteration result retained previous residual");
        const auto residual_scratch = jacobi.r_hypre_;
        fill_system(matrix, map, 1, b, truth);
        x = truth;
        check(jacobi.solve(matrix, b, x, 0, 192, 0, 192, map), "exact nonzero initial guess failed");
        verify(matrix, map, b, x, truth);
        check(jacobi.getLastIterations() == 0 && jacobi.r_hypre_ == residual_scratch,
              "zero-iteration residual scratch was not reused");
        Matrix empty;
        make_graph(empty, 0);
        const std::vector<HYPRE_BigInt> empty_map;
        std::vector<double> empty_rhs, empty_x;
        expect_failure([&] { jacobi.solve(empty, empty_rhs, empty_x, 0, 0, 0, 0, empty_map); });
        std::cout << "PASS: scale-independent acceptance and rejected stalled solves, refreshed values/zeros, true residuals, fresh-solve parity, stable resources, "
                     "relocation, resize, invalidation, frozen/Jacobi modes, rejected bad inputs, "
                     "timing and collective counts\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
