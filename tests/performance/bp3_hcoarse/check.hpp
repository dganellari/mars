#pragma once

#include "backend/distributed/unstructured/solvers/mars_cellwise_h_coarse.hpp"
#include <thrust/inner_product.h>

namespace mars::cellwise::h_coarse {

inline void report_memory(MPI_Comm comm)
{
    size_t free_bytes = 0, total_bytes = 0;
    MARS_CELLWISE_CK(cudaMemGetInfo(&free_bytes, &total_bytes));
    unsigned long long used = total_bytes - free_bytes, free = free_bytes, max_used = 0, min_free = 0;
    MARS_CELLWISE_MPI(MPI_Allreduce(&used, &max_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, comm));
    MARS_CELLWISE_MPI(MPI_Allreduce(&free, &min_free, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, comm));
    int rank = 0;
    MARS_CELLWISE_MPI(MPI_Comm_rank(comm, &rank));
    if (rank == 0) std::printf("BP3 h memory after setup: max used %llu bytes, min free %llu bytes; not peak\n", max_used, min_free);
}

// Coordinate-based values also check the topology's physical boundary independently.
__global__ void check_values_kernel(double* v, cw::Block box, unsigned seed, int kind)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= box.elements() * 8) return;
    const long long e = t / 8;
    const int x = box.ox + int(e / ((long long)box.ny * box.nz));
    const int y = box.oy + int(e / box.nz % box.ny), z = box.oz + int(e % box.nz);
    const int a = x + int(t % 8) / 4, b = y + int(t % 4) / 2, c = z + int(t % 2);
    const bool boundary = a == 0 || b == 0 || c == 0 || a == box.NX || b == box.NY || c == box.NZ;
    const unsigned long long key = kind == 1 ? ((long long)a * (box.NY + 1) + b) * (box.NZ + 1) + c
        : (((long long)x * box.NY + y) * box.NZ + z) * 8 + t % 8;
    unsigned h = unsigned(key) ^ seed;
    h ^= h >> 16; h *= 0x7feb352dU; h ^= h >> 15; h *= 0x846ca68bU; h ^= h >> 16;
    v[t] = kind == 3 ? 1.0 : boundary ? 0.0 : kind == 2 ? 1.0 : double(h) / 4294967296.0 - 0.5;
}

__global__ void check_divide_kernel(double* v, const double* count, long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) v[t] /= count[t];
}

struct Check {
    MPI_Comm comm;
    cudaStream_t stream;
    double worst = 0;

    double dot(const thrust::device_vector<double>& a, const thrust::device_vector<double>& b)
    {
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
        const double local = thrust::inner_product(thrust::device, a.begin(), a.end(), b.begin(), 0.0);
        double global = 0;
        MARS_CELLWISE_MPI(MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, comm));
        mars::fem::hypreCheck(comm, std::isfinite(global), "nonfinite h validation dot");
        return global;
    }

    void require(double error, const char* message)
    {
        mars::fem::hypreCheck(comm, std::isfinite(error) && error < 5e-11, message);
        worst = std::max(worst, error);
    }

    void values(thrust::device_vector<double>& a, const cw::Block& box, unsigned seed, int kind)
    {
        check_values_kernel<<<blocks(a.size()), cw::kThreads, 0, stream>>>(raw(a), box, seed, kind);
        MARS_CELLWISE_CK(cudaGetLastError());
    }

    void same(const thrust::device_vector<double>& a, const thrust::device_vector<double>& b, const char* message)
    {
        thrust::device_vector<double> difference(a.size());
        MARS_CELLWISE_CK(cudaDeviceSynchronize());
        cw::pmg_detail::axpby_kernel<<<blocks(a.size()), cw::kThreads, 0, stream>>>(
            raw(difference), raw(a), raw(b), -1.0, a.size());
        require(std::sqrt(dot(difference, difference) / std::max(1e-300, dot(a, a))), message);
    }

    void continuous(Level& level, const cw::Block& box, const thrust::device_vector<double>& x)
    {
        thrust::device_vector<double> ones(x.size()), count(x.size()), sum(x.size());
        MARS_CELLWISE_CK(cudaDeviceSynchronize());
        values(ones, box, 0, 3);
        cw::dss<1>(raw(ones), raw(count), *level.tables, level.halo, stream);
        cw::dss<1>(raw(x), raw(sum), *level.tables, level.halo, stream);
        check_divide_kernel<<<blocks(x.size()), cw::kThreads, 0, stream>>>(raw(sum), raw(count), x.size());
        same(x, sum, "h correction has unequal shared copies");
    }
};

inline void validate(Hierarchy& hierarchy, const cw::Decomposition& dec, cudaStream_t stream)
{
    Check test{dec.comm, stream};
    auto box = dec.local;
    for (std::size_t k = 0; k < hierarchy.size(); ++k) {
        Level& level = hierarchy.level(k);
        thrust::device_vector<double> mask(level.x.size()), expected(level.x.size());
        MARS_CELLWISE_CK(cudaDeviceSynchronize());
        test.values(mask, box, 0, 3);
        cw::zero_dirichlet<1>(raw(mask), *level.tables, stream);
        test.values(expected, box, 0, 2);
        test.same(expected, mask, "h topology boundary differs from global coordinates");
        if (k + 1 < hierarchy.size()) {
            Level& child = hierarchy.level(k + 1);
            thrust::device_vector<double> c(child.x.size()), c2(child.x.size()), r(child.x.size()),
                rap(child.x.size()), ac(child.x.size()), ac2(child.x.size()),
                f(level.x.size()), pc(level.x.size()), ap(level.x.size());
            MARS_CELLWISE_CK(cudaDeviceSynchronize());
            test.values(c, box.coarse(), 17, 1);
            test.values(c2, box.coarse(), 199, 1);
            test.values(f, box, 71, 0);
            prolong_kernel<false><<<blocks(pc.size()), cw::kThreads, 0, stream>>>(raw(c), raw(pc), level.grid);
            restrict_kernel<<<blocks(r.size()), cw::kThreads, 0, stream>>>(raw(f), raw(r), level.grid);
            test.require(std::abs(test.dot(f, pc) - test.dot(r, c)) /
                std::sqrt(test.dot(f, f) * test.dot(pc, pc)), "h transfers are not adjoint");
            test.continuous(level, box, pc);
            level.apply(raw(pc), raw(ap), stream);
            restrict_kernel<<<blocks(rap.size()), cw::kThreads, 0, stream>>>(raw(ap), raw(rap), level.grid);
            child.apply(raw(c), raw(ac), stream);
            child.apply(raw(c2), raw(ac2), stream);
            test.same(rap, ac, "h cached matrix differs from transpose(P) A P");
            test.require(std::abs(test.dot(c2, ac) - test.dot(c, ac2)) /
                std::sqrt(test.dot(c2, c2) * test.dot(ac, ac)), "h Galerkin matrix is not symmetric");
            mars::fem::hypreCheck(dec.comm, test.dot(c, ac) > 0, "h Galerkin energy is not positive");
        }
        box = box.coarse();
    }
    Level& fine = hierarchy.level(0);
    thrust::device_vector<double> a(fine.x.size()), b(a.size()), x(a.size()), y(a.size()), again(a.size()),
                                  sum(a.size()), combined(a.size());
    MARS_CELLWISE_CK(cudaDeviceSynchronize());
    test.values(a, dec.local, 81, 0);
    test.values(b, dec.local, 123, 0);
    hierarchy.solve(raw(a), raw(x));
    hierarchy.solve(raw(b), raw(y));
    hierarchy.solve(raw(a), raw(again));
    test.same(x, again, "h cycle retained an earlier correction");
    cw::zero_dirichlet<1>(raw(again), *fine.tables, stream);
    test.same(x, again, "h correction violates the Dirichlet boundary");
    const double ax = test.dot(a, x), by = test.dot(b, y);
    mars::fem::hypreCheck(dec.comm, ax > 0 && by > 0, "h cycle has nonpositive curvature");
    test.require(std::abs(test.dot(a, y) - test.dot(b, x)) / std::sqrt(ax * by), "h cycle is not symmetric");
    test.continuous(fine, dec.local, x);
    cw::pmg_detail::axpby_kernel<<<blocks(a.size()), cw::kThreads, 0, stream>>>(raw(sum), raw(a), raw(b), 2.0, a.size());
    hierarchy.solve(raw(sum), raw(combined));
    cw::pmg_detail::axpby_kernel<<<blocks(a.size()), cw::kThreads, 0, stream>>>(raw(sum), raw(x), raw(y), 2.0, a.size());
    test.same(sum, combined, "h cycle is not linear");
    MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
    if (dec.rank == 0) std::printf("PASS h validation: levels=%zu, adjoint/Galerkin/symmetry/linearity/continuity/repeat worst %.6e\n",
                                  hierarchy.size() - 1, test.worst);
}

} // namespace mars::cellwise::h_coarse
