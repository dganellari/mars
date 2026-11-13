#pragma once

// Diagnostics for the Galerkin driver; every dot product on unknowns counts owned DOFs once.
#include "backend/distributed/unstructured/amr/mars_amr_tet_audit_math.hpp"

namespace amr_audit
{
using Vector = cstone::DeviceVector<double>;

template<class F>
double sum(size_t n, F f)
{
    double value = mars::transformReduce(thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n), f,
                                         0.0, thrust::plus<double>());
    cudaCheckCall(cudaDeviceSynchronize());
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return value;
}

inline void require(bool valid, const char* check, int level)
{
    int bad = !valid;
    MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (!bad) return;
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) std::fprintf(stderr, "[amr-algebra-audit] FAIL level=%d check=%s\n", level, check);
    MPI_Abort(MPI_COMM_WORLD, 78);
    std::abort();
}

__host__ __device__ inline double random_value(uint64_t key, uint64_t salt)
{
    return double(mars::amr::audit_mix(key ^ salt) >> 11) * 0x1.0p-52 - 1.0;
}

template<int N, class Space>
void check_prolongation(Space& space,
                        Corners<KeyType, N> conn,
                        size_t first,
                        size_t own,
                        size_t n,
                        const KeyType* keys,
                        const double* x,
                        const double* y,
                        const double* z,
                        double cell_size,
                        int rank,
                        int level)
{
    const auto* synced = space.isSynced();
    const auto* dof    = space.isDof();
    double missing =
        sum(own * N, [=] __host__ __device__(size_t i) { return synced[conn.c[i % N][first + i / N]] ? 0.0 : 1.0; });
    require(missing == 0, "owned-corner coverage", level);
    Vector px(n), y_slot(n), restricted(n);
    double* px_data = px.data();
    double* y_data  = y_slot.data();
    thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                     [=] __host__ __device__(size_t i)
                     {
                         px_data[i] = dof[i] ? random_value(keys[i], 17) : 0;
                         y_data[i]  = synced[i] ? random_value(keys[i], uint64_t(rank + 1) * 91) : 0;
                     });
    space.prolong(px);
    if (n) cudaCheckCall(cudaMemcpy(restricted.data(), y_data, n * sizeof(double), cudaMemcpyDeviceToDevice));
    space.restrict(restricted);
    const double* pt_y = restricted.data();
    double lhs         = sum(n, [=] __host__ __device__(size_t i) { return synced[i] ? px_data[i] * y_data[i] : 0; });
    double rhs = sum(n, [=] __host__ __device__(size_t i) { return dof[i] ? random_value(keys[i], 17) * pt_y[i] : 0; });
    double scale = sum(n, [=] __host__ __device__(size_t i) { return synced[i] ? fabs(px_data[i] * y_data[i]) : 0; });
    double error = std::abs(lhs - rhs) / std::max(scale, 1e-300);
    if (rank == 0)
        std::printf("[amr-algebra-audit] level=%d adjoint_relative=%.17g scale=%.17g\n", level, error, scale);
    require(std::isfinite(error) && error <= 1e-10, "P/P-transpose adjoint", level);

    for (int basis = 0; basis < 4; ++basis)
    {
        const double* coordinate = basis == 1 ? x : basis == 2 ? y : z;
        thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                          px.data(),
                          [=] __host__ __device__(size_t i) { return dof[i] ? (basis ? coordinate[i] : 1.0) : 0; });
        space.prolong(px);
        double maximum = mars::transformReduce(
            thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
            [=] __host__ __device__(size_t i)
            {
                if (!synced[i]) return 0.0;
                return isfinite(px_data[i]) ? fabs(px_data[i] - (basis ? coordinate[i] : 1.0)) : INFINITY;
            },
            0.0, thrust::maximum<double>());
        MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        const double tolerance = basis ? 8 * cell_size : 1e-12;
        if (rank == 0)
            std::printf("[amr-algebra-audit] level=%d basis=%d max=%.17g tolerance=%.17g\n", level, basis, maximum,
                        tolerance);
        require(std::isfinite(maximum) && maximum <= tolerance, "affine basis prolongation", level);
    }
}

template<class Space>
void apply_slots(Space& space,
                 size_t n,
                 const int* index,
                 const double* removed,
                 const int* row_ptr,
                 const int* columns,
                 const double* values,
                 const double* packed,
                 double* result)
{
    const auto* dof = space.isDof();
    Vector input(n), output(n);
    thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                      input.data(),
                      [=] __host__ __device__(size_t s) { return dof[s] && !removed[s] ? packed[index[s]] : 0; });
    space.prolong(input);
    const double* u = input.data();
    thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                      output.data(),
                      [=] __host__ __device__(size_t s)
                      {
                          double a = 0;
                          for (int k = row_ptr[s]; k < row_ptr[s + 1]; ++k)
                              a += values[k] * u[columns[k]];
                          return a;
                      });
    space.restrict(output);
    const double* a = output.data();
    thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                     [=] __host__ __device__(size_t s)
                     {
                         if (dof[s]) result[index[s]] = removed[s] ? packed[index[s]] : a[s];
                     });
}

inline double difference(size_t n, const double* a, const double* b)
{
    double diff  = sum(n,
                       [=] __host__ __device__(size_t i)
                       {
                          double d = a[i] - b[i];
                          return d * d;
                      });
    double scale = sum(n, [=] __host__ __device__(size_t i) { return a[i] * a[i] + b[i] * b[i]; });
    return std::sqrt(diff / std::max(scale, 1e-300));
}

template<class Space>
void check_operator(Space& space,
                    HypreAmgPcgSolver& solver,
                    size_t n,
                    size_t owned,
                    const KeyType* keys,
                    const int* index,
                    const double* removed,
                    const int* row_ptr,
                    const int* columns,
                    const double* values,
                    const double* rhs,
                    const double* x,
                    const double* y,
                    const double* z,
                    Problem p,
                    int rank,
                    int level)
{
    const auto* dof = space.isDof();
    Vector first(owned), first_action(owned), second(owned), action(owned), reference(owned);
    for (int trial = 0; trial < 3; ++trial)
    {
        double* input  = trial == 0 ? first.data() : second.data();
        double* actual = trial == 0 ? first_action.data() : action.data();
        thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                         [=] __host__ __device__(size_t s)
                         {
                             if (dof[s])
                                 input[index[s]] = removed[s]   ? 0
                                                   : trial == 2 ? p.u(x[s], y[s], z[s])
                                                                : random_value(keys[s], uint64_t(47 + trial * 29));
                         });
        apply_slots(space, n, index, removed, row_ptr, columns, values, input, reference.data());
        cudaCheckCall(cudaDeviceSynchronize());
        solver.apply(input, actual);
        cudaCheckCall(cudaDeviceSynchronize());
        require(HYPRE_GetError() == 0, "Hypre matrix action", level);
        double error = difference(owned, actual, reference.data());
        if (rank == 0)
            std::printf("[amr-algebra-audit] level=%d trial=%d operator_relative=%.17g\n", level, trial, error);
        require(std::isfinite(error) && error <= 1e-10, "assembled P-transpose A P action", level);
        if (trial == 0)
        {
            double energy = sum(owned, [=] __host__ __device__(size_t i) { return input[i] * actual[i]; });
            require(std::isfinite(energy) && energy > 0, "positive energy", level);
        }
        if (trial == 1)
        {
            const double* u  = first.data();
            const double* au = first_action.data();
            double left      = sum(owned, [=] __host__ __device__(size_t i) { return u[i] * actual[i]; });
            double right     = sum(owned, [=] __host__ __device__(size_t i) { return input[i] * au[i]; });
            double scale     = sum(owned, [=] __host__ __device__(size_t i)
                                   { return fabs(u[i] * actual[i]) + fabs(input[i] * au[i]); });
            double symmetry  = std::abs(left - right) / std::max(scale, 1e-300);
            if (rank == 0) std::printf("[amr-algebra-audit] level=%d symmetry_relative=%.17g\n", level, symmetry);
            require(std::isfinite(symmetry) && symmetry <= 1e-10, "matrix symmetry", level);
        }
        if (trial == 2)
        {
            double assembled = difference(owned, actual, rhs);
            double composed  = difference(owned, reference.data(), rhs);
            if (rank == 0)
                std::printf("[amr-algebra-audit] level=%d exact_residual_assembled=%.17g "
                            "exact_residual_composed=%.17g\n",
                            level, assembled, composed);
            require(std::isfinite(assembled) && std::isfinite(composed), "finite exact-field residual", level);
        }
    }
}

template<class Space>
void check_solution(Space& space,
                    HypreAmgPcgSolver& solver,
                    size_t n,
                    size_t owned,
                    const int* index,
                    const double* removed,
                    const int* row_ptr,
                    const int* columns,
                    const double* values,
                    const double* rhs,
                    const double* solution,
                    double tolerance,
                    int rank,
                    int level)
{
    Vector actual(owned), composed(owned);
    apply_slots(space, n, index, removed, row_ptr, columns, values, solution, composed.data());
    cudaCheckCall(cudaDeviceSynchronize());
    solver.apply(solution, actual.data());
    cudaCheckCall(cudaDeviceSynchronize());
    require(HYPRE_GetError() == 0, "Hypre residual action", level);
    double norm2  = sum(owned, [=] __host__ __device__(size_t i) { return rhs[i] * rhs[i]; });
    auto residual = [&](const double* product)
    {
        return std::sqrt(sum(owned,
                             [=] __host__ __device__(size_t i)
                             {
                                 double r = rhs[i] - product[i];
                                 return r * r;
                             }) /
                         (norm2 > 0 ? norm2 : 1));
    };
    double a = residual(actual.data()), b = residual(composed.data());
    if (rank == 0)
        std::printf("[amr-algebra-audit] level=%d true_residual_assembled=%.17g "
                    "true_residual_composed=%.17g\n",
                    level, a, b);
    require(std::isfinite(a) && std::isfinite(b) && a <= 10 * tolerance && b <= 10 * tolerance,
            "true solution residual", level);
}
} // namespace amr_audit
