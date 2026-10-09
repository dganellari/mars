#pragma once
// Layout of element-local (cell-wise) vectors on one GPU and across GPUs.
//
// One structured block of hexahedra at p = 7, split over a grid of MPI ranks. Every
// element stores its own copy of each node it touches. Sharing, copy counts and the
// Dirichlet boundary follow the GLOBAL position of an element, so a node on a rank
// boundary is shared exactly as it would be on one GPU.
//
// The DSS gathers every copy of a node in the order of the dimensionally split
// cascade (x pairs, then y pairs of those, then z). Copies held by other ranks come
// from a ghost exchange: each rank sends the raw shared layer facing each of its up
// to 26 neighbours (a face plane, an edge line or a corner node per element) and the
// gather reads them in place of the missing neighbour elements. The result is
// bit-identical to the one-GPU DSS of the same global vector (checked in
// marsir-mlir/test/cellwise_distributed_ref.py). Message sizes follow from the block
// layout on both sides, so sends and receives match by construction.

#include <cuda_runtime.h>
#include <mpi.h>

#include <cstdio>
#include <cstdlib>
#include <vector>

namespace mars {
namespace cellwise {

constexpr int kN = 8, kNN = 64, kN3 = 512;   // nodes per edge, face, element (p = 7)
constexpr int kThreads = 256;                // threads per block of every kernel here

#define MARS_CELLWISE_CK(call)                                                   \
    do {                                                                         \
        cudaError_t e_ = (call);                                                 \
        if (e_ != cudaSuccess) {                                                 \
            fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e_),  \
                    __FILE__, __LINE__);                                         \
            std::abort();                                                        \
        }                                                                        \
    } while (0)

#define MARS_CELLWISE_MPI(call)                                                  \
    do {                                                                         \
        int e_ = (call);                                                         \
        if (e_ != MPI_SUCCESS) {                                                 \
            char s_[MPI_MAX_ERROR_STRING];                                       \
            int l_ = 0;                                                          \
            MPI_Error_string(e_, s_, &l_);                                       \
            fprintf(stderr, "MPI error %s at %s:%d\n", s_, __FILE__, __LINE__);  \
            MPI_Abort(MPI_COMM_WORLD, 1);                                        \
        }                                                                        \
    } while (0)

// A rank's part of a structured block: nx * ny * nz local elements at global element
// offset (ox, oy, oz) of an NX * NY * NZ block. Local element (ex, ey, ez) is
// ((ex * ny) + ey) * nz + ez; its node (a, b, c) is (a * 8 + b) * 8 + c with axis a
// along x.
struct Block {
    int nx = 0, ny = 0, nz = 0;
    int ox = 0, oy = 0, oz = 0;
    int NX = 0, NY = 0, NZ = 0;
    Block() = default;
    __host__ __device__ Block(int x, int y, int z) : nx(x), ny(y), nz(z), NX(x), NY(y), NZ(z) {}
    __host__ __device__ Block(int x, int y, int z, int px, int py, int pz, int gx, int gy, int gz)
        : nx(x), ny(y), nz(z), ox(px), oy(py), oz(pz), NX(gx), NY(gy), NZ(gz)
    {
    }
    __host__ __device__ long long elements() const { return (long long)nx * ny * nz; }
    __host__ __device__ long long values() const { return elements() * kN3; }
    // The same rank's part one multigrid level down: every dimension halved.
    __host__ __device__ Block coarse() const
    {
        return Block(nx / 2, ny / 2, nz / 2, ox / 2, oy / 2, oz / 2, NX / 2, NY / 2, NZ / 2);
    }
    __host__ __device__ bool halvable() const
    {
        return nx % 2 == 0 && ny % 2 == 0 && nz % 2 == 0;
    }
};

struct Node {
    long long e, t;   // local element, and value index e * kN3 + local node index
    int a, b, c;      // local node coordinates
    int ex, ey, ez;   // local element coordinates
};

// Along one axis, node coordinate a of the element at global index g among n: -1
// when the node is shared with element g - 1, +1 with g + 1, 0 when only this
// element has it.
__host__ __device__ inline int share_side(int a, int g, int n)
{
    if (a == 0) return g > 0 ? -1 : 0;
    if (a == kN - 1) return g < n - 1 ? 1 : 0;
    return 0;
}

__device__ inline bool on_outer_boundary(const Node& nd, const Block& blk)
{
    auto outer = [](int a, int g, int n) { return (a == 0 && g == 0) || (a == kN - 1 && g == n - 1); };
    return outer(nd.a, nd.ex + blk.ox, blk.NX) || outer(nd.b, nd.ey + blk.oy, blk.NY) ||
           outer(nd.c, nd.ez + blk.oz, blk.NZ);
}

// 1 / (number of copies): a power of two, so scaling by it is exact.
__device__ inline double weight(const Node& nd, const Block& blk)
{
    double w = 1.0;
    if (share_side(nd.a, nd.ex + blk.ox, blk.NX)) w *= 0.5;
    if (share_side(nd.b, nd.ey + blk.oy, blk.NY)) w *= 0.5;
    if (share_side(nd.c, nd.ez + blk.oz, blk.NZ)) w *= 0.5;
    return w;
}

// The local elements a kernel pass visits: all of them, the box [lo, hi) whose copies
// are all on this rank, or the shell around that box, which needs ghost copies. The
// shell is six slabs (x below and above the box, then y, then z, each inside the
// previous ones), so every shell element is visited once.
struct ElementSet {
    enum Kind { kAll, kBox, kShell };
    Kind kind = kAll;
    int lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};

    __host__ __device__ void slab(int s, const Block& b, int (&from)[3], int (&to)[3]) const
    {
        const int n[3] = {b.nx, b.ny, b.nz};
        const int axis = s / 2;
        for (int k = 0; k < 3; ++k) {
            if (k < axis) { from[k] = lo[k]; to[k] = hi[k]; }
            else if (k > axis) { from[k] = 0; to[k] = n[k]; }
            else if (s % 2 == 0) { from[k] = 0; to[k] = lo[k]; }
            else { from[k] = hi[k]; to[k] = n[k]; }
        }
    }
    __host__ __device__ static long long volume(const int (&from)[3], const int (&to)[3])
    {
        long long v = 1;
        for (int k = 0; k < 3; ++k) v *= to[k] > from[k] ? to[k] - from[k] : 0;
        return v;
    }
    __host__ __device__ long long count(const Block& b) const
    {
        if (kind == kAll) return b.elements();
        if (kind == kBox) return volume(lo, hi);
        long long c = 0;
        for (int s = 0; s < 6; ++s) {
            int from[3], to[3];
            slab(s, b, from, to);
            c += volume(from, to);
        }
        return c;
    }
    __host__ __device__ void element(long long k, const Block& b, int& ex, int& ey, int& ez) const
    {
        int from[3] = {0, 0, 0}, to[3] = {b.nx, b.ny, b.nz};
        if (kind == kBox) {
            for (int i = 0; i < 3; ++i) { from[i] = lo[i]; to[i] = hi[i]; }
        } else if (kind == kShell) {
            for (int s = 0; s < 6; ++s) {
                slab(s, b, from, to);
                const long long v = volume(from, to);
                if (k < v) break;
                k -= v;
            }
        }
        const int dy = to[1] - from[1], dz = to[2] - from[2];
        ez = from[2] + (int)(k % dz);
        ey = from[1] + (int)((k / dz) % dy);
        ex = from[0] + (int)(k / ((long long)dz * dy));
    }
};

// Calls f(node) for each node of the set's elements that this thread handles. Each
// thread block owns a contiguous range of the set, so neighbouring elements are read
// by the same block close in time. A box (all elements, or the interior) is walked
// by stepping the coordinates, without a division per element; the shell (SHELL)
// looks each element up in its slabs.
template <bool SHELL = false, typename F>
__device__ inline void for_each_node(const Block& blk, const ElementSet& set, F&& f)
{
    const long long count = set.count(blk);
    const long long chunk = (count + gridDim.x - 1) / gridDim.x;
    const long long begin = (long long)blockIdx.x * chunk;
    const long long end = begin + chunk < count ? begin + chunk : count;
    if (begin >= end) return;
    Node nd;
    const bool all = set.kind == ElementSet::kAll;
    const int lo[3] = {all ? 0 : set.lo[0], all ? 0 : set.lo[1], all ? 0 : set.lo[2]};
    const int hi[3] = {all ? blk.nx : set.hi[0], all ? blk.ny : set.hi[1], all ? blk.nz : set.hi[2]};
    if (!SHELL) set.element(begin, blk, nd.ex, nd.ey, nd.ez);
    for (long long k = begin; k < end; ++k) {
        if (SHELL) set.element(k, blk, nd.ex, nd.ey, nd.ez);
        nd.e = ((long long)nd.ex * blk.ny + nd.ey) * blk.nz + nd.ez;
        for (int l = threadIdx.x; l < kN3; l += kThreads) {
            nd.t = nd.e * kN3 + l;
            nd.a = l / kNN;
            nd.b = (l / kN) % kN;
            nd.c = l % kN;
            f(nd);
        }
        if (!SHELL && ++nd.ez == hi[2]) {
            nd.ez = lo[2];
            if (++nd.ey == hi[1]) {
                nd.ey = lo[1];
                ++nd.ex;
            }
        }
    }
}

template <typename F>
__device__ inline void for_each_node(const Block& blk, F&& f)
{
    for_each_node<false>(blk, ElementSet{}, f);
}

// Neighbour directions d = (dx, dy, dz) in {-1, 0, 1}^3 as slots 0..26, slot 13 being
// the rank itself; the opposite direction is slot 26 - s.
__host__ __device__ inline int direction_slot(int dx, int dy, int dz)
{
    return (dx + 1) * 9 + (dy + 1) * 3 + dz + 1;
}

// Where a rank finds the copies of its shared nodes held by other ranks: one region
// per neighbour direction, holding that neighbour's shared layer. A region spans the
// axes where the direction is 0: elements in (x, y, z) order, times the 8 nodes per
// such axis in (a, b, c) order. offset[slot] < 0 when there is no neighbour.
struct GhostView {
    const double* data = nullptr;
    long long offset[27] = {};
};

// Copy (a, b, c) of the element at local coordinates (x, y, z); a coordinate one past
// either end of the local block names the neighbouring rank's element, read from the
// ghost regions. value(i) gives local value i.
template <typename Value>
__device__ inline double fetch_with_ghosts(const Value& value, const GhostView& g, const Block& blk,
                                           int x, int y, int z, int a, int b, int c)
{
    const int dx = x < 0 ? -1 : (x >= blk.nx ? 1 : 0);
    const int dy = y < 0 ? -1 : (y >= blk.ny ? 1 : 0);
    const int dz = z < 0 ? -1 : (z >= blk.nz ? 1 : 0);
    if (dx == 0 && dy == 0 && dz == 0)
        return value((((long long)x * blk.ny + y) * blk.nz + z) * kN3 + (a * kN + b) * kN + c);
    long long elem = 0, node = 0, nodes = 1;
    if (dx == 0) { elem = elem * blk.nx + x; node = node * kN + a; nodes *= kN; }
    if (dy == 0) { elem = elem * blk.ny + y; node = node * kN + b; nodes *= kN; }
    if (dz == 0) { elem = elem * blk.nz + z; node = node * kN + c; nodes *= kN; }
    return g.data[g.offset[direction_slot(dx, dy, dz)] + elem * nodes + node];
}

// Sum of every copy of the node, in the cascade's order: x pairs first, then y pairs
// of those sums, then z pairs. On a shared axis the pair is (low element, high
// element); the low one holds the node on its last plane, the high one on plane 0.
// fetch(x, y, z, a, b, c) returns copy (a, b, c) of local element (x, y, z).
template <typename Fetch>
__device__ inline double gather_sum_at(const Fetch& fetch, const Node& nd, const Block& blk)
{
    const int sx = share_side(nd.a, nd.ex + blk.ox, blk.NX);
    const int sy = share_side(nd.b, nd.ey + blk.oy, blk.NY);
    const int sz = share_side(nd.c, nd.ez + blk.oz, blk.NZ);
    const int lx = nd.ex - (sx < 0), ly = nd.ey - (sy < 0), lz = nd.ez - (sz < 0);
    const int a0 = sx ? kN - 1 : nd.a, b0 = sy ? kN - 1 : nd.b, c0 = sz ? kN - 1 : nd.c;
    double zsum = 0.0;
#pragma unroll
    for (int k = 0; k < 2; ++k) {
        if (k == 1 && !sz) break;
        double ysum = 0.0;
#pragma unroll
        for (int j = 0; j < 2; ++j) {
            if (j == 1 && !sy) break;
            const int bb = j ? 0 : b0, cc = k ? 0 : c0;
            double xsum = fetch(lx, ly + j, lz + k, a0, bb, cc);
            if (sx) xsum += fetch(lx + 1, ly + j, lz + k, 0, bb, cc);
            ysum = j == 0 ? xsum : ysum + xsum;
        }
        zsum = k == 0 ? ysum : zsum + ysum;
    }
    return zsum;
}

// The gather of value(i) over this rank's copies only (SHELL = false: every copy must
// be local), or over local and ghost copies (SHELL = true).
template <bool SHELL, typename Value>
__device__ inline double gather_value(const Value& value, const GhostView& g, const Node& nd,
                                      const Block& blk)
{
    if constexpr (SHELL) {
        return gather_sum_at(
            [&](int x, int y, int z, int a, int b, int c) {
                return fetch_with_ghosts(value, g, blk, x, y, z, a, b, c);
            },
            nd, blk);
    } else {
        return gather_sum_at(
            [&](int x, int y, int z, int a, int b, int c) {
                return value((((long long)x * blk.ny + y) * blk.nz + z) * kN3 + (a * kN + b) * kN + c);
            },
            nd, blk);
    }
}

// The ranks as a PX * PY * PZ grid, rank (rx * PY + ry) * PZ + rz, each owning an
// equal sub-block of the global NX * NY * NZ block.
struct Decomposition {
    MPI_Comm comm = MPI_COMM_SELF;
    int rank = 0, size = 1;
    int P[3] = {1, 1, 1}, coord[3] = {0, 0, 0};
    int neighbour[27];   // rank in each direction, -1 for none (and for slot 13)
    Block local;

    Decomposition(MPI_Comm c, int NX, int NY, int NZ) : comm(c)
    {
        MARS_CELLWISE_MPI(MPI_Comm_rank(comm, &rank));
        MARS_CELLWISE_MPI(MPI_Comm_size(comm, &size));
        int dims[3] = {0, 0, 0};
        MARS_CELLWISE_MPI(MPI_Dims_create(size, 3, dims));
        const int N[3] = {NX, NY, NZ};
        for (int k = 0; k < 3; ++k) {
            P[k] = dims[k];
            if (N[k] % P[k] != 0) {
                if (rank == 0)
                    fprintf(stderr, "decomposition: %d elements along axis %d do not split over %d ranks\n",
                            N[k], k, P[k]);
                MPI_Abort(comm, 1);
            }
        }
        coord[2] = rank % P[2];
        coord[1] = (rank / P[2]) % P[1];
        coord[0] = rank / (P[1] * P[2]);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    const int q[3] = {coord[0] + dx, coord[1] + dy, coord[2] + dz};
                    const bool inside = q[0] >= 0 && q[0] < P[0] && q[1] >= 0 && q[1] < P[1] &&
                                        q[2] >= 0 && q[2] < P[2];
                    const bool self = dx == 0 && dy == 0 && dz == 0;
                    neighbour[direction_slot(dx, dy, dz)] =
                        inside && !self ? (q[0] * P[1] + q[1]) * P[2] + q[2] : -1;
                }
        const int n[3] = {NX / P[0], NY / P[1], NZ / P[2]};
        local = Block(n[0], n[1], n[2], coord[0] * n[0], coord[1] * n[1], coord[2] * n[2], NX, NY, NZ);
    }
};

// Region sizes and offsets of the 27 direction slots for one block layout: shared by
// the pack kernel (send side) and the GhostView (receive side).
struct RegionTable {
    long long offset[27];
    long long count[27];
};

// Packs value = b - q (b alone when q is null) of the shared layer facing each
// neighbour: one thread per packed value, all regions in one launch.
__global__ void __launch_bounds__(kThreads)
pack_kernel(const double* __restrict__ b, const double* __restrict__ q, double* __restrict__ out,
            RegionTable table, long long total, Block blk)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= total) return;
    int s = 0;
    while (t >= table.offset[s] + table.count[s]) ++s;
    const int d[3] = {s / 9 - 1, (s / 3) % 3 - 1, s % 3 - 1};
    const int n[3] = {blk.nx, blk.ny, blk.nz};
    long long i = t - table.offset[s];
    int nodes = 1;
    for (int k = 0; k < 3; ++k)
        if (d[k] == 0) nodes *= kN;
    long long elem = i / nodes, node = i % nodes;
    int e[3], l[3];
    for (int k = 2; k >= 0; --k) {   // decode the last axis first
        if (d[k] == 0) {
            e[k] = (int)(elem % n[k]);
            elem /= n[k];
            l[k] = (int)(node % kN);
            node /= kN;
        } else {
            e[k] = d[k] > 0 ? n[k] - 1 : 0;
            l[k] = d[k] > 0 ? kN - 1 : 0;
        }
    }
    const long long v = (((long long)e[0] * blk.ny + e[1]) * blk.nz + e[2]) * kN3 + (l[0] * kN + l[1]) * kN + l[2];
    out[t] = q ? b[v] - q[v] : b[v];
}

// The ghost exchange for one block layout (one multigrid level). start() packs and
// posts the messages; the caller runs the elements that need no ghosts; finish()
// waits. All on device buffers: GPU-aware MPI (on Cray MPICH,
// MPICH_GPU_SUPPORT_ENABLED=1).
class Halo {
public:
    Halo(const Decomposition& dec, const Block& blk) : comm_(dec.comm), blk_(blk)
    {
        const int n[3] = {blk.nx, blk.ny, blk.nz};
        long long total = 0;
        for (int s = 0; s < 27; ++s) {
            neighbour_[s] = dec.neighbour[s];
            table_.offset[s] = total;
            table_.count[s] = 0;
            if (neighbour_[s] < 0) continue;
            const int d[3] = {s / 9 - 1, (s / 3) % 3 - 1, s % 3 - 1};
            long long c = 1;
            for (int k = 0; k < 3; ++k)
                if (d[k] == 0) c *= (long long)n[k] * kN;
            table_.count[s] = c;
            total += c;
        }
        total_ = total;
        for (int s = 0; s < 27; ++s) view_.offset[s] = neighbour_[s] < 0 ? -1 : table_.offset[s];
        // A neighbour below along an axis means the elements at the low end of that
        // axis need ghosts; likewise above.
        for (int k = 0; k < 3; ++k) {
            const int below = neighbour_[direction_slot(k == 0 ? -1 : 0, k == 1 ? -1 : 0, k == 2 ? -1 : 0)];
            const int above = neighbour_[direction_slot(k == 0 ? 1 : 0, k == 1 ? 1 : 0, k == 2 ? 1 : 0)];
            interior_.lo[k] = below >= 0 ? 1 : 0;
            interior_.hi[k] = n[k] - (above >= 0 ? 1 : 0);
            if (interior_.lo[k] > n[k]) interior_.lo[k] = n[k];
            if (interior_.hi[k] < interior_.lo[k]) interior_.hi[k] = interior_.lo[k];
        }
        interior_.kind = ElementSet::kBox;
        shell_ = interior_;
        shell_.kind = ElementSet::kShell;
        requests_.reserve(2 * 26);
        if (total_ > 0) {
            MARS_CELLWISE_CK(cudaMalloc(&send_, total_ * sizeof(double)));
            MARS_CELLWISE_CK(cudaMalloc(&recv_, total_ * sizeof(double)));
        }
        view_.data = recv_;
    }
    ~Halo()
    {
        cudaFree(send_);
        cudaFree(recv_);
    }
    Halo(const Halo&) = delete;
    Halo& operator=(const Halo&) = delete;

    bool active() const { return total_ > 0; }
    const GhostView& view() const { return view_; }
    const ElementSet& interior() const { return interior_; }
    const ElementSet& shell() const { return shell_; }

    // Region s goes to the neighbour in direction s, which files it as its region
    // for the opposite direction: the tag carries the sender's slot.
    void start(const double* b, const double* q, cudaStream_t stream)
    {
        pack_kernel<<<(unsigned)((total_ + kThreads - 1) / kThreads), kThreads, 0, stream>>>(
            b, q, send_, table_, total_, blk_);
        MARS_CELLWISE_CK(cudaGetLastError());
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));   // MPI reads send_ next
        requests_.clear();
        for (int s = 0; s < 27; ++s) {
            if (neighbour_[s] < 0) continue;
            requests_.emplace_back();
            MARS_CELLWISE_MPI(MPI_Irecv(recv_ + table_.offset[s], (int)table_.count[s], MPI_DOUBLE,
                                        neighbour_[s], kTagBase + (26 - s), comm_, &requests_.back()));
        }
        for (int s = 0; s < 27; ++s) {
            if (neighbour_[s] < 0) continue;
            requests_.emplace_back();
            MARS_CELLWISE_MPI(MPI_Isend(send_ + table_.offset[s], (int)table_.count[s], MPI_DOUBLE,
                                        neighbour_[s], kTagBase + s, comm_, &requests_.back()));
        }
    }
    void finish()
    {
        MARS_CELLWISE_MPI(MPI_Waitall((int)requests_.size(), requests_.data(), MPI_STATUSES_IGNORE));
    }

private:
    static constexpr int kTagBase = 0x3100;   // clear of the other MARS halo tags
    MPI_Comm comm_;
    Block blk_;
    int neighbour_[27];
    RegionTable table_;
    long long total_ = 0;
    double *send_ = nullptr, *recv_ = nullptr;
    GhostView view_;
    ElementSet interior_, shell_;
    std::vector<MPI_Request> requests_;
};

// Runs a gather pass over the local block. On one rank (or without neighbours) that
// is one launch over all elements. Otherwise: start the exchange of value = b - q,
// launch the elements whose copies are all local while the messages travel, wait,
// then launch the shell. launch(set, shell, first) returns its number of thread
// blocks; `first` is how many blocks earlier launches used (for partial sums). The
// result is the total number of blocks.
template <typename Launch>
int run_gather_pass(Halo* halo, const double* b, const double* q, cudaStream_t stream,
                    Launch&& launch)
{
    if (!halo || !halo->active()) return launch(ElementSet{}, false, 0);
    halo->start(b, q, stream);
    const int first = launch(halo->interior(), false, 0);
    halo->finish();
    return first + launch(halo->shell(), true, first);
}

}  // namespace cellwise
}  // namespace mars
