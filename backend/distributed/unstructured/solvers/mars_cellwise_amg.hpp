#pragma once
// The p = 1 coarse solve of the cell-wise p-multigrid (mars_cellwise_pmg.hpp): BoomerAMG
// V-cycles on the assembled p = 1 operator. Needs a Hypre build.
//
// Each vertex belongs to the rank that holds its counting copy (the canonically first).
// That rank numbers its vertices densely after the lower ranks' (one MPI_Exscan) and
// writes number + 1 into its counting copies, so one DSS gives every copy on every rank
// its vertex's global number; Dirichlet vertices are not numbered and read -1. The same
// with the owner rank. The matrix comes from the 8 x 8 element matrices (the p = 1
// operator on unit vectors): entries of rows owned elsewhere go to their owners over the
// NBX sparseExchange, then hypreFromEntries builds the ParCSR matrix on the GPU. In a
// V-cycle the right-hand side is the DSS of the raw residual read at the counting copies,
// and the solution returns to every copy through the counting copies and a second DSS.

// domain.hpp first: the Hypre headers use its VectorSelector without including it.
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/mars_sfc_ownership.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_operator.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"
#include "backend/distributed/unstructured/solvers/mars_hypre_amg_pcg_solver.hpp"

#include <cuda_runtime.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/tuple.h>

#include <functional>
#include <span>
#include <vector>

namespace mars {
namespace cellwise {

namespace amg_detail {

using topo_detail::blocks;

// The corner of the p = 1 local node l (node (a, b, c) = (a * 2 + b) * 2 + c).
__device__ inline int corner_of(int l) { return corner_from_bits(l >> 2, (l >> 1) & 1, l & 1); }

// 1 for the copies of the vertices this rank numbers: counting copies, not Dirichlet.
__global__ void owned_kernel(TopologyView T, int* __restrict__ flag)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= T.local * 8) return;
    const long long e = t / 8;
    const int c = corner_of((int)(t % 8));
    flag[t] = ((T.counting[e] >> (kVertexBit + c)) & 1u) && !T.vert_dirichlet[T.vert_of[e * 8 + c]];
}

// v[copy[k]] = base + k + 1 (v is zero elsewhere).
__global__ void number_kernel(const int* __restrict__ copy, long long rows, long long base, double* __restrict__ v)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k < rows) v[copy[k]] = (double)(base + k + 1);
}

// v[copy[k]] = rank + 1 (v is zero elsewhere).
__global__ void owner_kernel(const int* __restrict__ copy, long long rows, int rank, double* __restrict__ v)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k < rows) v[copy[k]] = (double)(rank + 1);
}

__global__ void unit_kernel(double* __restrict__ v, long long E, int j)   // 1 at local node j of every element
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < E * 8) v[t] = t % 8 == j ? 1.0 : 0.0;
}

// One item per (element, row node i, column node j): entry A_e[i][j] at (number of i,
// number of j). dest: -1 for a row of this rank, the owner rank otherwise, -2 to drop
// (a Dirichlet row or column).
__global__ void entries_kernel(const double* __restrict__ number, const double* __restrict__ owner,
                               const double* __restrict__ columns, long long E, int rank, int* __restrict__ dest,
                               long long* __restrict__ row, long long* __restrict__ col, double* __restrict__ val)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * 64) return;
    const long long e = t / 64;
    const int i = (int)(t % 64) / 8, j = (int)(t % 8);
    const long long r = (long long)number[e * 8 + i] - 1, c = (long long)number[e * 8 + j] - 1;
    const int o = (int)owner[e * 8 + i] - 1;
    dest[t] = r < 0 || c < 0 ? -2 : (o == rank ? -1 : o);
    row[t] = r;
    col[t] = c;
    val[t] = columns[(long long)j * E * 8 + e * 8 + i];
}

// Functors and kernels rather than device lambdas: nvcc allows none inside a constructor.
struct NonZero {
    __device__ bool operator()(int f) const { return f != 0; }
};

// The remote entries as (row, column, value) doubles.
__global__ void pack_remote_kernel(const long long* __restrict__ row, const long long* __restrict__ col,
                                   const double* __restrict__ val, long long n, double* __restrict__ out)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    out[3 * k] = (double)row[k];
    out[3 * k + 1] = (double)col[k];
    out[3 * k + 2] = val[k];
}

// This rank's rows: its own entries, then the received ones; rows relative to first.
__global__ void local_entries_kernel(const long long* __restrict__ row, const long long* __restrict__ col,
                                     const double* __restrict__ val, long long mine, const double* __restrict__ got,
                                     long long n, long long first, int* __restrict__ orow,
                                     long long* __restrict__ ocol, double* __restrict__ oval)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    if (k < mine) {
        orow[k] = (int)(row[k] - first);
        ocol[k] = col[k];
        oval[k] = val[k];
    } else {
        const long long m = k - mine;
        orow[k] = (int)((long long)got[3 * m] - first);
        ocol[k] = (long long)got[3 * m + 1];
        oval[k] = got[3 * m + 2];
    }
}

__global__ void gather_kernel(const double* __restrict__ y, const int* __restrict__ copy, long long rows,
                              double* __restrict__ f)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k < rows) f[k] = y[copy[k]];
}

__global__ void scatter_kernel(const double* __restrict__ u, const int* __restrict__ copy, long long rows,
                               double* __restrict__ v)
{
    const long long k = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (k < rows) v[copy[k]] = u[k];
}

}  // namespace amg_detail

// The BoomerAMG options that set how deep and how wide the hierarchy is; the defaults are
// those of HypreAmgPcgSolver. All run on the GPU in Hypre 2.33: PMIS coarsening, and
// aggressive levels with 2-stage extended+e interpolation (Hypre's device default).
struct AmgSettings {
    double strong_threshold = 0.25;
    int aggressive_levels = 0;
    int max_coarse_size = 9;   // coarsening stops below this many rows
    bool print = false;        // Hypre's setup statistics (levels, complexities)
    // Groups of this many consecutive ranks hand their rows to the group's first rank,
    // which alone runs BoomerAMG: fewer, larger ranks per level (1: every rank solves).
    int ranks_per_solve = 1;
};

class P1Amg {
public:
    // y = A u at p = 1 on the local elements, for the element matrices.
    using Apply1 = std::function<void(const double*, double*, cudaStream_t)>;

    // G1, M1: the p = 1 CVFEM metric and matrices; halo1: the p = 1 exchange (null on one
    // rank); cycles: V-cycles per coarse solve.
    P1Amg(const UnstructuredTopology& T, const double* G1, const CvfemMatrices<1>& M1, UnstructuredHalo* halo1,
          MPI_Comm comm, cudaStream_t stream, int cycles = 1)
        : P1Amg(
              T,
              [G1, M1, E = T.local](const double* u, double* y, cudaStream_t s) { cvfem_apply<1>(u, y, G1, E, M1, s); },
              halo1, comm, stream, cycles)
    {
    }

    // Any p = 1 operator (e.g. the BP3 Galerkin form), symmetric for BoomerAMG's V-cycle
    // to be a symmetric preconditioner.
    P1Amg(const UnstructuredTopology& T, const Apply1& apply1, UnstructuredHalo* halo1, MPI_Comm comm,
          cudaStream_t stream, int cycles = 1, const AmgSettings& settings = {})
        : T_(T), halo_(halo1), stream_(stream), y_(T.local * 8), v_(T.local * 8)
    {
        using namespace amg_detail;
        int rank = 0;
        MARS_CELLWISE_MPI(MPI_Comm_rank(comm, &rank));
        const TopologyView tv = view(T);
        const long long E = T.local, n = E * 8;

        // The vertices this rank numbers, and their numbers.
        thrust::device_vector<int> flag(n, 0);
        if (n > 0) owned_kernel<<<blocks(n), kThreads, 0, stream_>>>(tv, raw(flag));
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
        copy_.resize(n);
        const auto ids = thrust::counting_iterator<int>(0);
        rows_ = thrust::copy_if(thrust::device, ids, ids + n, flag.begin(), copy_.begin(), NonZero{}) - copy_.begin();
        copy_.resize(rows_);
        // Rows are numbered group by group, and inside a group rank by rank, so each
        // leader's rows are contiguous and in the order of the gather.
        grouped_ = settings.ranks_per_solve > 1;
        const int leader_rank = rank / settings.ranks_per_solve * settings.ranks_per_solve;
        MARS_CELLWISE_MPI(MPI_Comm_split(comm, rank / settings.ranks_per_solve, rank, &group_));
        int g = 0, gsize = 1;
        MARS_CELLWISE_MPI(MPI_Comm_rank(group_, &g));
        MARS_CELLWISE_MPI(MPI_Comm_size(group_, &gsize));
        leader_ = g == 0;
        MARS_CELLWISE_MPI(MPI_Comm_split(comm, leader_ ? 0 : MPI_UNDEFINED, rank, &solvers_));
        long long rows = rows_, in_group = 0, group_rows = 0, group_first = 0, total = 0;
        MARS_CELLWISE_MPI(MPI_Exscan(&rows, &in_group, 1, MPI_LONG_LONG, MPI_SUM, group_));
        if (g == 0) in_group = 0;
        MARS_CELLWISE_MPI(MPI_Allreduce(&rows, &group_rows, 1, MPI_LONG_LONG, MPI_SUM, group_));
        if (leader_) {
            int s = 0;
            MARS_CELLWISE_MPI(MPI_Comm_rank(solvers_, &s));
            MARS_CELLWISE_MPI(MPI_Exscan(&group_rows, &group_first, 1, MPI_LONG_LONG, MPI_SUM, solvers_));
            if (s == 0) group_first = 0;
        }
        MARS_CELLWISE_MPI(MPI_Bcast(&group_first, 1, MPI_LONG_LONG, 0, group_));
        MARS_CELLWISE_MPI(MPI_Allreduce(&rows, &total, 1, MPI_LONG_LONG, MPI_SUM, comm));
        first_ = group_first + in_group;
        solve_first_ = group_first;
        solve_rows_ = leader_ ? group_rows : 0;
        // The leader's counts and offsets of the gather and scatter.
        counts_.assign(leader_ ? gsize : 0, 0);
        displs_.assign(leader_ ? gsize : 0, 0);
        const int rows32 = (int)rows_;
        MARS_CELLWISE_MPI(MPI_Gather(&rows32, 1, MPI_INT, counts_.data(), 1, MPI_INT, 0, group_));
        for (int i = 1; i < (int)counts_.size(); ++i) displs_[i] = displs_[i - 1] + counts_[i - 1];
        if (grouped_) own_.resize(rows_);

        // Every copy's vertex number and owner, through the counting copies and a DSS.
        thrust::device_vector<double> number(n), owner(n);
        MARS_CELLWISE_CK(cudaMemsetAsync(raw(y_), 0, n * sizeof(double), stream_));
        if (rows_ > 0) number_kernel<<<blocks(rows_), kThreads, 0, stream_>>>(raw(copy_), rows_, first_, raw(y_));
        dss<1>(raw(y_), raw(number), T_, halo_, stream_);
        MARS_CELLWISE_CK(cudaMemsetAsync(raw(y_), 0, n * sizeof(double), stream_));
        if (rows_ > 0) owner_kernel<<<blocks(rows_), kThreads, 0, stream_>>>(raw(copy_), rows_, leader_rank, raw(y_));
        dss<1>(raw(y_), raw(owner), T_, halo_, stream_);

        // The element matrices, one column per unit vector.
        thrust::device_vector<double> columns(8 * n);
        for (int j = 0; j < 8; ++j) {
            if (n > 0) unit_kernel<<<blocks(n), kThreads, 0, stream_>>>(raw(v_), E, j);
            apply1(raw(v_), raw(columns) + (long long)j * n, stream_);
        }

        // Entries: this rank's rows directly, the others' sent to their owners (the leaders).
        thrust::device_vector<int> dest(E * 64);
        thrust::device_vector<long long> row(E * 64), col(E * 64);
        thrust::device_vector<double> val(E * 64);
        if (E > 0)
            entries_kernel<<<blocks(E * 64), kThreads, 0, stream_>>>(raw(number), raw(owner), raw(columns), E, rank,
                                                                     raw(dest), raw(row), raw(col), raw(val));
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));
        thrust::sort_by_key(thrust::device, dest.begin(), dest.end(),
                            thrust::make_zip_iterator(thrust::make_tuple(row.begin(), col.begin(), val.begin())));
        const long long dropped = thrust::count(thrust::device, dest.begin(), dest.end(), -2);
        const long long mine = thrust::count(thrust::device, dest.begin(), dest.end(), -1);
        const long long remote = E * 64 - dropped - mine;
        // Remote entries as (row, column, value) doubles, grouped by owner (dest >= 0 is sorted last).
        thrust::device_vector<double> send(3 * remote);
        const long long r0 = dropped + mine;
        if (remote > 0)
            pack_remote_kernel<<<blocks(remote), kThreads>>>(raw(row) + r0, raw(col) + r0, raw(val) + r0, remote,
                                                             raw(send));
        MARS_CELLWISE_CK(cudaGetLastError());
        thrust::device_vector<int> peers(remote), counts(remote);
        const auto pe = thrust::reduce_by_key(thrust::device, dest.begin() + r0, dest.end(),
                                              thrust::constant_iterator<int>(3), peers.begin(), counts.begin());
        const long long np = pe.first - peers.begin();
        std::vector<int> h_peers(np), h_counts(np);
        thrust::copy(peers.begin(), peers.begin() + np, h_peers.begin());
        thrust::copy(counts.begin(), counts.begin() + np, h_counts.begin());
        MARS_CELLWISE_CK(cudaDeviceSynchronize());
        // Only rows on rank boundaries travel, so the received entries land on the host.
        std::vector<double> h_recv;
        constexpr int kTagEntries = 0x3350;   // clear of the halo (0x3300) and the domain tables (0x3340)
        sparseExchange<double>(std::span<const int>(h_peers), std::span<const int>(h_counts), raw(send), h_recv,
                               kTagEntries, comm);
        const thrust::device_vector<double> recv(h_recv.begin(), h_recv.end());
        const long long got = (long long)recv.size() / 3;

        thrust::device_vector<int> lrow(mine + got);
        thrust::device_vector<long long> lcol(mine + got);
        thrust::device_vector<HYPRE_Complex> lval(mine + got);
        if (mine + got > 0)
            local_entries_kernel<<<blocks(mine + got), kThreads>>>(raw(row) + dropped, raw(col) + dropped,
                                                                   raw(val) + dropped, mine, raw(recv), mine + got,
                                                                   solve_first_, raw(lrow), raw(lcol), raw(lval));
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaDeviceSynchronize());
        rows_global_ = total;
        if (!leader_) return;
        fem::hypreCheck(solvers_, solve_rows_ <= 0x7fffffffLL, "cell-wise p = 1 BoomerAMG: too many rows on one rank");
        A_ = fem::hypreFromEntries(solvers_, solve_first_, (HYPRE_Int)solve_rows_, total, lrow, lcol, lval);
#if defined(HYPRE_USING_GPU)
        hypre_CSRMatrixMoveDiagFirstDevice(hypre_ParCSRMatrixDiag(A_.get()));   // relaxation reads it first
#endif
        HYPRE_BigInt part[2] = {(HYPRE_BigInt)solve_first_, (HYPRE_BigInt)(solve_first_ + solve_rows_)};
        HYPRE_ParVectorCreate(solvers_, (HYPRE_BigInt)total, part, &f_);
        HYPRE_ParVectorCreate(solvers_, (HYPRE_BigInt)total, part, &u_);
        HYPRE_ParVectorInitialize(f_);
        HYPRE_ParVectorInitialize(u_);
        // The settings of HypreAmgPcgSolver (PMIS, extended+i, l1-Jacobi), as `cycles` V-cycles.
        HYPRE_BoomerAMGCreate(&amg_);
        HYPRE_BoomerAMGSetCoarsenType(amg_, 8);
        HYPRE_BoomerAMGSetInterpType(amg_, 6);
        HYPRE_BoomerAMGSetPMaxElmts(amg_, 4);
        HYPRE_BoomerAMGSetStrongThreshold(amg_, settings.strong_threshold);
        if (settings.aggressive_levels > 0) {
            HYPRE_BoomerAMGSetAggNumLevels(amg_, settings.aggressive_levels);
            HYPRE_BoomerAMGSetAggInterpType(amg_, 7);
        }
        HYPRE_BoomerAMGSetMaxCoarseSize(amg_, settings.max_coarse_size);
        HYPRE_BoomerAMGSetRelaxType(amg_, 18);
        HYPRE_BoomerAMGSetRelaxOrder(amg_, 0);
        HYPRE_BoomerAMGSetNumSweeps(amg_, 1);
        HYPRE_BoomerAMGSetCycleRelaxType(amg_, 18, 3);   // smooth on the coarsest level, no direct solve
        HYPRE_BoomerAMGSetCycleNumSweeps(amg_, 4, 3);
        HYPRE_BoomerAMGSetMaxLevels(amg_, 25);
        HYPRE_BoomerAMGSetTol(amg_, 0.0);
        HYPRE_BoomerAMGSetMaxIter(amg_, cycles);
        HYPRE_BoomerAMGSetPrintLevel(amg_, settings.print ? 1 : 0);
        const HYPRE_Int error = HYPRE_BoomerAMGSetup(amg_, matrix(), f_, u_);
        fem::hypreCheck(solvers_, error == 0 && HYPRE_GetError() == 0, "cell-wise p = 1 BoomerAMG setup failed");
    }
    ~P1Amg()
    {
        if (amg_) HYPRE_BoomerAMGDestroy(amg_);
        if (f_) HYPRE_ParVectorDestroy(f_);
        if (u_) HYPRE_ParVectorDestroy(u_);
        int finalized = 0;   // a caller may destroy this after MPI_Finalize
        MPI_Finalized(&finalized);
        if (finalized) return;
        if (group_ != MPI_COMM_NULL) MPI_Comm_free(&group_);
        if (solvers_ != MPI_COMM_NULL) MPI_Comm_free(&solvers_);
    }
    P1Amg(const P1Amg&) = delete;
    P1Amg& operator=(const P1Amg&) = delete;

    long long unknowns() const { return rows_global_; }

    // x1 = the coarse correction of the raw p = 1 residual b1, on every copy.
    void solve(const double* b1, double* x1) const
    {
        using namespace amg_detail;
        dss<1>(b1, raw(y_), T_, halo_, stream_);   // the assembled residual on every copy
        // In groups the rows go through own_ to the leader and back; alone, straight into Hypre's vectors.
        double* f = grouped_ ? raw(own_) : data(f_);
        if (rows_ > 0) gather_kernel<<<blocks(rows_), kThreads, 0, stream_>>>(raw(y_), raw(copy_), rows_, f);
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream_));   // MPI and Hypre (on its own stream) read it next
        if (grouped_)
            MARS_CELLWISE_MPI(MPI_Gatherv(f, (int)rows_, MPI_DOUBLE, leader_ ? data(f_) : nullptr, counts_.data(),
                                          displs_.data(), MPI_DOUBLE, 0, group_));
        if (leader_) {
            HYPRE_ParVectorSetConstantValues(u_, 0.0);
            HYPRE_BoomerAMGSolve(amg_, matrix(), f_, u_);
            MARS_CELLWISE_CK(cudaDeviceSynchronize());
        }
        double* u = grouped_ ? raw(own_) : data(u_);
        if (grouped_)
            MARS_CELLWISE_MPI(MPI_Scatterv(leader_ ? data(u_) : nullptr, counts_.data(), displs_.data(), MPI_DOUBLE, u,
                                           (int)rows_, MPI_DOUBLE, 0, group_));
        MARS_CELLWISE_CK(cudaMemsetAsync(raw(v_), 0, v_.size() * sizeof(double), stream_));
        if (rows_ > 0) scatter_kernel<<<blocks(rows_), kThreads, 0, stream_>>>(u, raw(copy_), rows_, raw(v_));
        MARS_CELLWISE_CK(cudaGetLastError());
        dss<1>(raw(v_), x1, T_, halo_, stream_);
    }

private:
    template <typename T>
    static T* raw(thrust::device_vector<T>& v) { return thrust::raw_pointer_cast(v.data()); }
    template <typename T>
    static const T* raw(const thrust::device_vector<T>& v) { return thrust::raw_pointer_cast(v.data()); }
    static double* data(HYPRE_ParVector v)
    {
        return hypre_VectorData(hypre_ParVectorLocalVector(reinterpret_cast<hypre_ParVector*>(v)));
    }
    HYPRE_ParCSRMatrix matrix() const { return reinterpret_cast<HYPRE_ParCSRMatrix>(A_.get()); }

    const UnstructuredTopology& T_;
    UnstructuredHalo* halo_;
    cudaStream_t stream_;
    mutable thrust::device_vector<double> y_, v_;   // scratch of the p = 1 size
    thrust::device_vector<int> copy_;                // the counting copy of every numbered vertex
    long long rows_ = 0, first_ = 0, rows_global_ = 0;
    // Groups: this rank's group and the leaders' communicator (MPI_COMM_NULL off the leaders),
    // the leader's rows, and its gather counts and offsets per group rank.
    bool grouped_ = false, leader_ = true;
    MPI_Comm group_ = MPI_COMM_NULL, solvers_ = MPI_COMM_NULL;
    long long solve_first_ = 0, solve_rows_ = 0;
    std::vector<int> counts_, displs_;
    mutable thrust::device_vector<double> own_;   // this rank's rows on their way to and from the leader
    fem::HypreMatrix A_;
    HYPRE_ParVector f_ = nullptr, u_ = nullptr;
    HYPRE_Solver amg_ = nullptr;
};

}  // namespace cellwise
}  // namespace mars
