#pragma once
#include "backend/distributed/unstructured/fem/mars_nonconformal_mesh_faces.hpp"
#include "backend/distributed/unstructured/fem/mars_nonconformal_search.hpp"
#include "backend/distributed/unstructured/fem/mars_outlet_fgmres.hpp"
#include <array>
#include <cstdio>
#include <numeric>
#include <set>
#if defined(__CUDACC__) || defined(__HIPCC__)
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#endif
using namespace mars::fem;
using namespace mars::fem::nonconformal_search_detail;
#if defined(__CUDACC__) || defined(__HIPCC__)
constexpr bool device = true;
using Graph           = NonconformalScalarGraph<double>;
#else
constexpr bool device = false;
using Graph           = NonconformalScalarGraphHost<double>;
#endif
template<class T>
using Buffer = Array<T, device>;
static int rankId, rankCount;
static void check(bool ok, const char* message)
{
    int local = ok, global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (!global)
    {
        if (!ok) std::fprintf(stderr, "rank %d FAIL: %s\n", rankId, message);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
template<class T>
static std::vector<T> host(const Buffer<T>& a)
{
    std::vector<T> b(a.size());
#if defined(__CUDACC__) || defined(__HIPCC__)
    thrust::copy(a.begin(), a.end(), b.begin());
#else
    std::copy(a.begin(), a.end(), b.begin());
#endif
    return b;
}
static void completed()
{
#if defined(__CUDACC__) || defined(__HIPCC__)
    bool ok = cudaGetLastError() == cudaSuccess;
    ok      = (cudaDeviceSynchronize() == cudaSuccess) && ok;
    check(ok, "device completion");
#endif
}
struct Mesh
{
    std::vector<Vec3<double>> x;
    std::vector<int> cells, owners;
    std::vector<int64_t> dofs;
    std::vector<uint64_t> elements;
    std::vector<double> k;
    std::vector<InterfaceFaceTag> tags, outer;
    std::vector<int> side, boundary;
    int globalNodes  = 0;
    bool partitioned = false;
    std::vector<int> candidatePeers;
};
static int owner(int64_t id, int n)
{
    for (int r = 0; r < rankCount; ++r)
        if (id < int64_t(n) * (r + 1) / rankCount) return r;
    return -1;
}
static std::vector<int> peers()
{
    std::vector<int> p;
    for (int r = 0; r < rankCount; ++r)
        if (r != rankId) p.push_back(r);
    return p;
}
static int global_node_count(const Mesh& m) { return m.globalNodes ? m.globalNodes : int(m.x.size()); }
static std::vector<int> peers(const Mesh& m) { return m.partitioned ? m.candidatePeers : peers(); }
static std::vector<int64_t> ownedRows(const Mesh& m)
{
    std::vector<int64_t> rows;
    for (int i = int(int64_t(global_node_count(m)) * rankId / rankCount);
         i < int(int64_t(global_node_count(m)) * (rankId + 1) / rankCount); ++i)
        rows.push_back(i);
    return rows;
}
struct Extracted
{
    Buffer<InterfaceTetFace<double>> faces;
    std::vector<uint64_t> ids;
    std::vector<int> owners;
};
static Extracted extract(const Mesh& m)
{
    std::vector<InterfaceFaceTag> tags, outer;
    for (auto tag : m.tags)
    {
        bool need = m.owners[tag.element] == rankId;
        for (int j = 0; j < 4; ++j)
            if (j != tag.opposite)
                need = need || owner(m.dofs[m.cells[4 * tag.element + j]], global_node_count(m)) == rankId;
        if (need) tags.push_back(tag);
    }
    Buffer<Vec3<double>> x(m.x.begin(), m.x.end());
    Buffer<int> cells(m.cells.begin(), m.cells.end()), owners(m.owners.begin(), m.owners.end());
    Buffer<int64_t> dofs(m.dofs.begin(), m.dofs.end());
    Buffer<uint64_t> elements(m.elements.begin(), m.elements.end());
    Buffer<double> k(m.k.begin(), m.k.end());
    Buffer<InterfaceFaceTag> dtags(tags.begin(), tags.end());
    Extracted out;
    out.faces.resize(tags.size());
    Buffer<uint64_t> ids(tags.size());
    Buffer<int> own(tags.size()), valid(tags.size());
    InterfaceTetMeshView<double> view{int(m.x.size()), int(m.elements.size()), pointer(x),      pointer(cells),
                                      pointer(dofs),   pointer(elements),      pointer(owners), pointer(k)};
    each<device>(int(tags.size()), ExtractInterfaceFace<double>{view, pointer(dtags), pointer(out.faces), pointer(ids),
                                                                pointer(own), pointer(valid)});
    completed();
    auto hv = host<int>(valid);
    check(std::find(hv.begin(), hv.end(), 0) == hv.end(), "mesh side-set extraction");
    // GhostRegistry currently consumes host identities; face geometry stays on device.
    out.ids    = host<uint64_t>(ids);
    out.owners = host<int>(own);
    return out;
}
struct Dirichlet
{
    InterfaceCsrView<double> a;
    const Vec3<double>* x;
    const int* boundaryCodes;
    double* b;
    double* inverseDiagonal;
    MARS_FS_FN void operator()(int r) const
    {
        auto row      = a.rowIds[r];
        bool boundary = boundaryCodes ? boundaryCodes[row] != 0 : std::abs(x[row].x) == 1.;
        b[r]          = boundaryCodes ? double(boundaryCodes[row] == 2) : double(boundary && x[row].x > 0);
        for (int j = a.offsets[r]; j < a.offsets[r + 1]; ++j)
        {
            if (boundary) a.values[j] = a.columns[j] == row ? 1. : 0.;
            if (a.columns[j] == row) inverseDiagonal[r] = 1. / a.values[j];
        }
    }
};
struct VectorOp
{
    double* y;
    const double* x;
    double a;
    int mode;
    MARS_FS_FN void operator()(int i) const
    {
        if (mode == 0)
            y[i] = 0;
        else if (mode == 1)
            y[i] = x[i];
        else if (mode == 2)
            y[i] *= a;
        else if (mode == 3)
            y[i] += a * x[i];
        else
            y[i] = a;
    }
};
struct Multiply
{
    double* y;
    const double* x;
    const double* diagonal;
    MARS_FS_FN void operator()(int i) const { y[i] = x[i] * diagonal[i]; }
};
struct InvalidDiagonal
{
    MARS_FS_FN bool operator()(double value) const { return !nonconformal_detail::finite(value) || value == 0; }
};
struct Dot
{
    const double* x;
    const double* y;
    MARS_FS_FN double operator()(int i) const { return x[i] * y[i]; }
};
struct Matvec
{
    const int* offsets;
    const int* columns;
    const double* values;
    const double* x;
    double* y;
    MARS_FS_FN void operator()(int r) const
    {
        double s = 0;
        for (int j = offsets[r]; j < offsets[r + 1]; ++j)
            s += values[j] * x[columns[j]];
        y[r] = s;
    }
};
// Example operations share the FGMRES controller and donor exchange.
// Only owned rows enter global reductions.
class ChannelOps
{
public:
    ChannelOps(Graph& a, const Mesh& m, int restart, int components = 1, bool imposeBoundary = true)
        : halo(MPI_COMM_WORLD)
        , n(int(a.rowIds.size()))
    {
        check((components == 1 || components == 3) && (!imposeBoundary || components == 1), "linear component layout");
        auto rows = host<int64_t>(a.rowIds), cols = host<int64_t>(a.columns);
        std::vector<uint64_t> ids(rows.begin(), rows.end());
        std::set<int64_t> ghosts;
        for (auto c : cols)
            if (owner(c / components, global_node_count(m)) != rankId) ghosts.insert(c);
        ids.insert(ids.end(), ghosts.begin(), ghosts.end());
        width = int(ids.size());
        std::map<int64_t, int> slots;
        std::vector<int> owners;
        for (int i = 0; i < width; ++i)
        {
            slots[ids[i]] = i;
            owners.push_back(owner(ids[i] / components, global_node_count(m)));
        }
        std::vector<int> local;
        for (auto c : cols)
            local.push_back(slots.at(c));
        columns = Buffer<int>(local.begin(), local.end());
        check(halo.buildFromHost(ids, owners, peers(m)), "solution halo setup");
        offsets = &a.offsets;
        values  = &a.values;
        x.resize(width);
        b.resize(width);
        r.resize(width);
        w.resize(width);
        diagonal.resize(n);
        basisValues.resize(size_t(restart + 1) * width);
        directions.resize(size_t(restart) * width);
        if (!imposeBoundary) return;
        Buffer<Vec3<double>> coords(m.x.begin(), m.x.end());
        Buffer<int> codes(m.boundary.begin(), m.boundary.end());
        each<device>(n, Dirichlet{a.view(), pointer(coords), codes.empty() ? nullptr : pointer(codes), pointer(b),
                                  pointer(diagonal)});
        completed();
#if defined(__CUDACC__) || defined(__HIPCC__)
        check(thrust::count_if(diagonal.begin(), diagonal.end(), InvalidDiagonal{}) == 0, "finite Jacobi diagonal");
#else
        check(std::none_of(diagonal.begin(), diagonal.end(), InvalidDiagonal{}), "finite Jacobi diagonal");
#endif
    }
    double* solution() { return pointer(x); }
    double* residual() { return pointer(r); }
    double* work() { return pointer(w); }
    double* rhs() { return pointer(b); }
    double* basis(int i) { return pointer(basisValues) + size_t(i) * width; }
    double* direction(int i) { return pointer(directions) + size_t(i) * width; }
    void zero(double* y) { each<device>(width, VectorOp{y, nullptr, 0, 0}); }
    void copy(const double* a, double* y) { each<device>(n, VectorOp{y, a, 0, 1}); }
    void scale(double a, double* y) { each<device>(n, VectorOp{y, nullptr, a, 2}); }
    void axpy(double a, const double* v, double* y) { each<device>(n, VectorOp{y, v, a, 3}); }
    double localDot(const double* a, const double* b)
    {
#if defined(__CUDACC__) || defined(__HIPCC__)
        return thrust::transform_reduce(thrust::device, thrust::counting_iterator<int>(0),
                                        thrust::counting_iterator<int>(n), Dot{a, b}, 0., thrust::plus<double>());
#else
        double s = 0;
        for (int i = 0; i < n; ++i)
            s += a[i] * b[i];
        return s;
#endif
    }
    double norm(const double* a)
    {
        double local = localDot(a, a), global = 0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return std::sqrt(global);
    }
    bool precondition(const double* a, double* z)
    {
        each<device>(n, Multiply{z, a, pointer(diagonal)});
        return true;
    }
    void apply(double* a, double* y)
    {
#if defined(__CUDACC__) || defined(__HIPCC__)
        check(halo.forward(a, width), "solution halo refresh");
#else
        // The host registry exposes a vector API; this is the host oracle only.
        std::vector<double> field(a, a + width);
        check(halo.forwardHost(field), "solution halo refresh");
        std::copy(field.begin() + n, field.end(), a + n);
#endif
        each<device>(n, Matvec{pointer(*offsets), pointer(columns), pointer(*values), a, y});
    }
    void orthogonalize(double* y, int count, double* h)
    {
        std::fill(h, h + count, 0.);
        std::vector<double> local(count), global(count);
        for (int pass = 0; pass < 2; ++pass)
        {
            for (int j = 0; j < count; ++j)
                local[j] = localDot(y, basis(j));
            MPI_Allreduce(local.data(), global.data(), count, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            for (int j = 0; j < count; ++j)
            {
                h[j] += global[j];
                axpy(-global[j], basis(j), y);
            }
        }
    }
    void report(int, double) {}
    NonconformalDonorExchange<double, double> halo;
    int n, width;
    Buffer<double> x, b, r, w, diagonal, basisValues, directions;
    Buffer<int> columns;
    Buffer<int>* offsets;
    Buffer<double>* values;
};
