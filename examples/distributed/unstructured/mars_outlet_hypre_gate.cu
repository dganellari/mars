// GPT/Codex, 2026-09-12: actual Hypre reuse on an analytic distributed matrix.
#include "backend/distributed/unstructured/fem/mars_ns_pump_solver.hpp"
#include <numeric>

namespace {
using Solver = mars::fem::HypreGMRESSolver<double, int, cstone::GpuTag>;
using Matrix = Solver::Matrix;
using Vector = Solver::Vector;

void require(bool ok, const char* message)
{
    int bad = ok ? 0 : 1, global = 0;
    MPI_Allreduce(&bad, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (global) {
        std::cerr << "FAIL: " << message << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        std::abort();
    }
}

template<class VectorType, class T>
void upload(VectorType& destination, const std::vector<T>& source)
{
    destination.resize(source.size());
    if (!source.empty())
        require(cudaMemcpy(thrust::raw_pointer_cast(destination.data()), source.data(),
            source.size()*sizeof(T), cudaMemcpyHostToDevice) == cudaSuccess, "fixture upload");
}

std::vector<double> download(const Vector& source)
{
    std::vector<double> result(source.size());
    require(cudaMemcpy(result.data(), source.data(), result.size()*sizeof(double),
                       cudaMemcpyDeviceToHost) == cudaSuccess, "fixture download");
    return result;
}

void gate(int rank, int ranks)
{
    constexpr int global_rows = 512;
    const int n = global_rows/ranks, start = rank*n, end = start+n;
    Matrix matrix;
    matrix.allocate(n, global_rows, 3*n);
    std::vector<int> rows(n+1), columns(3*n);
    std::vector<double> values(3*n), rhs(n), exact(n);
    std::vector<HYPRE_BigInt> mapping(global_rows);
    std::iota(mapping.begin(), mapping.end(), HYPRE_BigInt(0));
    for (int i = 0; i < n; ++i) {
        rows[i] = 3*i;
        columns[3*i] = start+i;
        columns[3*i+1] = (start+i+global_rows-1)%global_rows;
        columns[3*i+2] = (start+i+1)%global_rows;
        values[3*i] = 4;
        values[3*i+1] = -1;
        values[3*i+2] = -.5;
    }
    rows[n] = 3*n;
    upload(matrix.rowOffsets(), rows);
    upload(matrix.colIndices(), columns);
    upload(matrix.values(), values);
    thrust::device_vector<HYPRE_BigInt> device_map;
    upload(device_map, mapping);
    Vector b, fresh_x, cached_x, cycle_x;
    auto fill_rhs = [&](int sample) {
        auto field = [sample](int i) {
            return sample == 1 ? -2*(1 + .2*std::sin(.1*i)) : 1 + .2*std::sin(.1*i + sample);
        };
        for (int i = 0; i < n; ++i) {
            exact[i] = field(start+i);
            rhs[i] = values[3*i]*exact[i] - field(columns[3*i+1]) - .5*field(columns[3*i+2]);
        }
        upload(b, rhs);
        const std::vector<double> zero(n, 0);
        upload(fresh_x, zero); upload(cached_x, zero); upload(cycle_x, zero);
    };
    Solver cached(MPI_COMM_WORLD, 500, 1e-11), cycle(MPI_COMM_WORLD, 500, 1e-11);
    cached.setVerbose(false); cycle.setVerbose(false);
    cached.enable_reuse(); cycle.enable_reuse(true);
    auto solve = [&](Solver& solver, Vector& x) {
        return solver.solve(matrix, b, x, start, end, 0, global_rows, device_map);
    };
    std::vector<double> previous_action;
    for (int sample = 0; sample < 3; ++sample) {
        if (sample == 2) {
            // In-place matrix changes require collective explicit invalidation.
            for (int i = 0; i < n; ++i) values[3*i] = 5;
            upload(matrix.values(), values);
            cached.invalidate_setup(); cycle.invalidate_setup();
        }
        fill_rhs(sample);
        Solver fresh(MPI_COMM_WORLD, 500, 1e-11), fresh_cycle(MPI_COMM_WORLD, 500, 1e-11);
        fresh.setVerbose(false); fresh_cycle.setVerbose(false);
        fresh_cycle.enable_reuse(true);
        require(solve(fresh, fresh_x), "fresh GMRES convergence");
        require(solve(cached, cached_x), "reused GMRES convergence");
        const auto reference = download(fresh_x), reused = download(cached_x);
        double error = 0;
        for (int i = 0; i < n; ++i)
            error = std::max(error, std::max(std::abs(reused[i]-reference[i]), std::abs(reused[i]-exact[i])));
        require(error < 1e-8, "reused solve agrees with fresh solve and analytic solution");
        require(solve(cycle, cycle_x), "cached AMG application");
        const auto first = download(cycle_x);
        if (sample == 1) {
            error = 0;
            for (int i = 0; i < n; ++i) error = std::max(error, std::abs(first[i]+2*previous_action[i]));
            require(error < 1e-10, "cached AMG uses the new scaled RHS");
        }
        previous_action = first;
        const std::vector<double> zero(n, 0);
        upload(cycle_x, zero);
        require(solve(cycle, cycle_x), "repeated AMG application");
        const auto second = download(cycle_x);
        upload(fresh_x, zero);
        require(solve(fresh_cycle, fresh_x), "fresh AMG application");
        const auto fresh_action = download(fresh_x);
        error = 0;
        for (int i = 0; i < n; ++i)
            error = std::max(error, std::abs(first[i]-second[i]));
        require(error < 1e-10, "cached AMG action is repeatable from a zero initial guess");
        // Fresh coarsening need not choose the same hierarchy. Check its action independently.
        double local_error = 0, local_rhs = 0;
        for (int i = 0; i < n; ++i) {
            local_error += (fresh_action[i]-exact[i])*(fresh_action[i]-exact[i]);
            local_rhs += exact[i]*exact[i];
        }
        double sums[2] = {local_error, local_rhs}, global_sums[2] = {};
        MPI_Allreduce(sums, global_sums, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        require(std::isfinite(global_sums[0]) && global_sums[0] < global_sums[1],
                "fresh AMG action improves this diagonally dominant fixture");
        require(cached.get_setup_count() == (sample == 2 ? 2 : 1)
             && cycle.get_setup_count() == (sample == 2 ? 2 : 1), "setup reuse and invalidation counts");
    }
    // A changed device map allocation triggers the collective metadata check.
    auto alternate_map = device_map;
    fill_rhs(3);
    const auto& selected_map = rank == 0 ? alternate_map : device_map;
    require(cached.solve(matrix, b, cached_x, start, end, 0, global_rows, selected_map), "map-storage rebuild");
    require(cached.get_setup_count() == 3, "map-storage change invalidates setup");
    require(cudaDeviceSynchronize() == cudaSuccess, "asynchronous CUDA completion");
    if (rank == 0) std::cout << "PASS: prepared Hypre and AMG-cycle gate ranks=" << ranks << '\n';
}
} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, ranks = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ranks);
    require(ranks == 1 || ranks == 2 || ranks == 4, "use 1/2/4 ranks");
    // Match the existing outlet gate: the launcher controls GPU visibility.
    require(cudaFree(nullptr) == cudaSuccess, "initialize visible GPU");
    gate(rank, ranks);
    MPI_Finalize();
}
