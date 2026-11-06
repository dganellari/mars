#pragma once

#include "mars_cellwise_h_transfer.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_amg.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_pmg.hpp"
#include <thrust/extrema.h>
#include <memory>

namespace mars::cellwise::h_coarse {

namespace cw = mars::cellwise;
using Apply = cw::P1Amg::Apply1;

template<class T> T* raw(thrust::device_vector<T>& v) { return thrust::raw_pointer_cast(v.data()); }
template<class T> const T* raw(const thrust::device_vector<T>& v) { return thrust::raw_pointer_cast(v.data()); }
inline unsigned blocks(long long n) { return unsigned((n + cw::kThreads - 1) / cw::kThreads); }

template<bool Add>
__global__ void prolong_kernel(const double* coarse, double* fine, Grid grid)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < grid.elements() * 8) {
        const double v = prolong_value(coarse, grid, t);
        fine[t] = Add ? fine[t] + v : v;
    }
}

__global__ void restrict_kernel(const double* fine, double* coarse, Grid grid)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < grid.coarse().elements() * 8) coarse[t] = restrict_value(fine, grid, t);
}

// Store columns in the same layout used by P1Amg's element assembly.
__global__ void matrix_apply_kernel(const double* matrix, const double* u, double* y, long long elements)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= elements * 8) return;
    double value = 0.0;
    for (int j = 0; j < 8; ++j) value += matrix[(long long)j * elements * 8 + t] * u[t / 8 * 8 + j];
    y[t] = value;
}

__global__ void bound_column_kernel(const double* column, double* diagonal, double* row_abs, long long n, int j)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) {
        row_abs[t] += fabs(column[t]);
        if (t % 8 == j) diagonal[t] = column[t];
    }
}

__global__ void bound_ratio_kernel(const double* row_abs, double* ratio, const double* diagonal,
                                    cw::TopologyView topology)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= topology.local * 8) return;
    const int corner = cw::amg_detail::corner_of(int(t % 8));
    const bool boundary = topology.vert_dirichlet[topology.vert_of[t / 8 * 8 + corner]];
    ratio[t] = boundary ? 0.0 : diagonal[t] > 0.0 && isfinite(diagonal[t]) && isfinite(row_abs[t])
        ? row_abs[t] / diagonal[t] : INFINITY;
}

__global__ void topology_keys_kernel(const unsigned long long* gid, long long count, int ne,
                                      unsigned long long* const* keys, int* const* ids)
{
    const long long e = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= count) return;
    const long long g = (long long)gid[e], side = ne + 1;
    for (int corner = 0; corner < 8; ++corner) {
        int x, y, z;
        cw::corner_bits(corner, x, y, z);
        const long long v = ((g / ((long long)ne * ne) + x) * side + g / ne % ne + y) * side + g % ne + z;
        keys[corner][e] = (unsigned long long)v;
        ids[corner][e] = int(v);
    }
}

inline cw::UnstructuredTopology topology(const cw::Decomposition& dec,
                                         thrust::device_vector<unsigned long long>& gid,
                                         thrust::device_vector<int>& owner, cudaStream_t stream)
{
    cw::block_elements(dec, gid, owner);
    thrust::device_vector<unsigned long long> keys[8];
    thrust::device_vector<int> ids[8];
    std::vector<unsigned long long*> key_ptrs(8);
    std::vector<int*> id_ptrs(8);
    for (int c = 0; c < 8; ++c) {
        keys[c].resize(gid.size()); ids[c].resize(gid.size());
        key_ptrs[c] = raw(keys[c]); id_ptrs[c] = raw(ids[c]);
    }
    thrust::device_vector<unsigned long long*> dk(key_ptrs);
    thrust::device_vector<int*> di(id_ptrs);
    MARS_CELLWISE_CK(cudaDeviceSynchronize()); // Setup buffers were initialized on the default stream.
    topology_keys_kernel<<<blocks(gid.size()), cw::kThreads, 0, stream>>>(raw(gid), gid.size(), dec.local.NX, raw(dk), raw(di));
    MARS_CELLWISE_CK(cudaGetLastError());
    MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
    std::vector<const unsigned long long*> ck(key_ptrs.begin(), key_ptrs.end());
    std::vector<const int*> ci(id_ptrs.begin(), id_ptrs.end());
    thrust::device_vector<const unsigned long long*> dck(ck);
    thrust::device_vector<const int*> dci(ci);
    return cw::build_topology(raw(dck), raw(dci), gid.size(), dec.local.elements());
}

struct Level {
    Grid grid;
    std::unique_ptr<cw::UnstructuredTopology> owned_topology;
    std::unique_ptr<cw::UnstructuredHalo> owned_halo;
    const cw::UnstructuredTopology* tables;
    cw::UnstructuredHalo* halo;
    Apply apply;
    thrust::device_vector<double> matrix, b, x, residual, y, direction, diagonal;
    double bound = 0.0;

    Level(Grid grid, const cw::UnstructuredTopology* tables, cw::UnstructuredHalo* halo, Apply apply)
        : grid(grid), tables(tables), halo(halo), apply(std::move(apply)), b(grid.elements() * 8),
          x(b.size()), residual(b.size()), y(b.size()), direction(b.size()), diagonal(b.size())
    { MARS_CELLWISE_CK(cudaDeviceSynchronize()); }

    void prepare_smoother(MPI_Comm comm, cudaStream_t stream)
    {
        const long long n = grid.elements() * 8;
        MARS_CELLWISE_CK(cudaMemsetAsync(raw(residual), 0, n * sizeof(double), stream));
        for (int j = 0; j < 8; ++j) {
            cw::amg_detail::unit_kernel<<<blocks(n), cw::kThreads, 0, stream>>>(raw(x), grid.elements(), j);
            apply(raw(x), raw(y), stream);
            bound_column_kernel<<<blocks(n), cw::kThreads, 0, stream>>>(raw(y), raw(b), raw(residual), n, j);
        }
        MARS_CELLWISE_CK(cudaGetLastError());
        cw::dss<1>(raw(b), raw(diagonal), *tables, halo, stream);
        cw::dss<1>(raw(residual), raw(y), *tables, halo, stream);
        bound_ratio_kernel<<<blocks(n), cw::kThreads, 0, stream>>>(raw(y), raw(residual), raw(diagonal), cw::view(*tables));
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
        const double local = *thrust::max_element(thrust::device, residual.begin(), residual.end());
        MARS_CELLWISE_MPI(MPI_Allreduce(&local, &bound, 1, MPI_DOUBLE, MPI_MAX, comm));
        mars::fem::hypreCheck(comm, std::isfinite(bound) && bound > 0.0, "invalid geometric smoother bound");
    }
};

class Hierarchy {
public:
    Hierarchy(const cw::Decomposition& dec, const cw::UnstructuredTopology& fine_topology,
              cw::UnstructuredHalo* fine_halo, const Apply& fine_apply, int depth, int cycles,
              const cw::AmgSettings& settings, cudaStream_t stream)
        : comm_(dec.comm), stream_(stream)
    {
        const auto& box = dec.local;
        mars::fem::hypreCheck(comm_, depth >= 1 && depth <= 3, "h depth must be between one and three");
        const int factor = 1 << depth;
        mars::fem::hypreCheck(comm_, box.nx % factor == 0 && box.ny % factor == 0 &&
                             box.nz % factor == 0 && box.nx / factor >= 2 && box.ny / factor >= 2 && box.nz / factor >= 2,
                             "h-coarsening needs aligned blocks with at least two cells per axis on the last level");
        levels_.emplace_back(new Level({box.nx, box.ny, box.nz}, &fine_topology, fine_halo, fine_apply));
        for (int k = 1; k <= depth; ++k) {
            Level& parent = *levels_.back();
            parent.prepare_smoother(comm_, stream_);
            const int ne = box.NX >> k;
            cw::Decomposition coarse_dec(comm_, ne, ne, ne);
            thrust::device_vector<unsigned long long> gid;
            thrust::device_vector<int> owner;
            auto tables = std::make_unique<cw::UnstructuredTopology>(topology(coarse_dec, gid, owner, stream_));
            std::unique_ptr<cw::UnstructuredHalo> halo;
            if (dec.size > 1) halo.reset(new cw::UnstructuredHalo(*tables, raw(gid), raw(owner), comm_, 1));
            auto child = std::make_unique<Level>(parent.grid.coarse(), tables.get(), halo.get(), Apply{});
            child->owned_topology = std::move(tables);
            child->owned_halo = std::move(halo);
            const long long coarse_n = child->grid.elements() * 8;
            child->matrix.resize(coarse_n * 8);
            MARS_CELLWISE_CK(cudaDeviceSynchronize());
            for (int j = 0; j < 8; ++j) {
                cw::amg_detail::unit_kernel<<<blocks(coarse_n), cw::kThreads, 0, stream_>>>(raw(child->x), child->grid.elements(), j);
                prolong_kernel<false><<<blocks(parent.x.size()), cw::kThreads, 0, stream_>>>(raw(child->x), raw(parent.x), parent.grid);
                parent.apply(raw(parent.x), raw(parent.y), stream_);
                restrict_kernel<<<blocks(coarse_n), cw::kThreads, 0, stream_>>>(raw(parent.y), raw(child->matrix) + j * coarse_n, parent.grid);
            }
            MARS_CELLWISE_CK(cudaGetLastError());
            Level* ptr = child.get();
            child->apply = [ptr](const double* in, double* out, cudaStream_t s) {
                matrix_apply_kernel<<<blocks(ptr->grid.elements() * 8), cw::kThreads, 0, s>>>(raw(ptr->matrix), in, out, ptr->grid.elements());
                MARS_CELLWISE_CK(cudaGetLastError());
            };
            levels_.push_back(std::move(child));
        }
        Level& last = *levels_.back();
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
        bottom_.reset(new cw::P1Amg(*last.tables, last.apply, last.halo, comm_, stream_, cycles, settings));
        const long long side = (box.NX >> depth) - 1;
        mars::fem::hypreCheck(comm_, bottom_->unknowns() == side * side * side, "unexpected geometric coarse row count");
        if (dec.rank == 0) {
            std::printf("BP3 geometric levels: %d, global cells per axis", depth);
            for (int k = 0; k <= depth; ++k) std::printf(" %d", box.NX >> k);
            std::printf("; bottom rows %lld; smoother bounds", bottom_->unknowns());
            for (int k = 0; k < depth; ++k) std::printf(" %.6g", levels_[k]->bound);
            std::printf("\n");
        }
    }

    void solve(const double* rhs, double* correction) { cycle(0, rhs, correction); }
    std::size_t size() const { return levels_.size(); }
    Level& level(std::size_t k) { return *levels_.at(k); }

private:
    template<bool First>
    void step(Level& level, const double* rhs, const double* y, double* x, double c1, double c2)
    {
        cw::run_unstructured_pass(level.halo, level.grid.elements(), rhs, y, stream_,
            [&](const cw::ElementList& set, const double* ghost, int) {
                static const int grid = cw::resident_grid(cw::unstructured_chebyshev_kernel<1, First>);
                cw::unstructured_chebyshev_kernel<1, First><<<grid, cw::kThreads, 0, stream_>>>(
                    rhs, y, ghost, raw(level.diagonal), raw(level.direction), x, c1, c2, cw::view(*level.tables), set);
                MARS_CELLWISE_CK(cudaGetLastError());
                return grid;
            });
    }

    void smooth(Level& level, const double* rhs, double* x, bool zero)
    {
        cw::chebyshev_steps(level.bound, 2, [&](int i, double c1, double c2) {
            const double* y = nullptr;
            if (i || !zero) { level.apply(x, raw(level.y), stream_); y = raw(level.y); }
            if (i == 0) step<true>(level, rhs, y, x, c1, c2);
            else step<false>(level, rhs, y, x, c1, c2);
        });
    }

    void cycle(std::size_t index, const double* rhs, double* x)
    {
        if (index + 1 == levels_.size()) { bottom_->solve(rhs, x); return; }
        Level& level = *levels_[index];
        Level& child = *levels_[index + 1];
        smooth(level, rhs, x, true);
        level.apply(x, raw(level.y), stream_);
        cw::pmg_detail::axpby_kernel<<<blocks(level.y.size()), cw::kThreads, 0, stream_>>>(
            raw(level.residual), rhs, raw(level.y), -1.0, level.y.size());
        cw::zero_dirichlet<1>(raw(level.residual), *level.tables, stream_);
        restrict_kernel<<<blocks(child.b.size()), cw::kThreads, 0, stream_>>>(raw(level.residual), raw(child.b), level.grid);
        MARS_CELLWISE_CK(cudaGetLastError());
        cycle(index + 1, raw(child.b), raw(child.x));
        prolong_kernel<true><<<blocks(level.y.size()), cw::kThreads, 0, stream_>>>(raw(child.x), x, level.grid);
        MARS_CELLWISE_CK(cudaGetLastError());
        smooth(level, rhs, x, false);
    }

    MPI_Comm comm_;
    cudaStream_t stream_;
    std::vector<std::unique_ptr<Level>> levels_;
    std::unique_ptr<cw::P1Amg> bottom_; // Destroy before the topology and halo it borrows.
};

} // namespace mars::cellwise::h_coarse
