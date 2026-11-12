// Poisson on an adapted non-conforming mesh: Galerkin Q1 on hexahedra, P1 on tetrahedra.
//
// Every level refines the elements near a point and keeps the others, so nodes hang in the middle of the edges
// (and faces) of the coarser elements next to the refined ones. HangingNodes finds them and DofSpace gives each
// the mean of the nodes it hangs on, which keeps the field continuous: the unknowns are the regular nodes, the
// matrix is P^T A P (DofSpace::restrictMatrix) with the Dirichlet values eliminated symmetrically, and PCG with
// BoomerAMG solves it.
//
//   --problem=linear  u = 1 + x + 2y + 3z. The discrete space holds u, so the solution equals it at every node,
//                     hanging ones included, up to where the SFC keys put the nodes (1e-6 here): a direct check of
//                     the constraints, which would be off by 1e-2 and more.
//   --problem=peak    u = exp(-alpha |x - c|^2) with the refinement around c: the L2 error at each level.
//
// The mesh must be the unit cube (its faces are the Dirichlet boundary). Example, three refinements on 4 ranks:
//   srun -n 4 ./mars_amr_galerkin --mesh=cube16.mesh --element=hex --problem=linear --levels=3

#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/amr/mars_amr.hpp"
#include "backend/distributed/unstructured/fem/mars_dof_space.hpp"
#include "backend/distributed/unstructured/fem/mars_element_assembly.hpp"
#include "backend/distributed/unstructured/solvers/mars_hypre_amg_pcg_solver.hpp"
#include "backend/distributed/unstructured/utils/mars_gpu_reduce.hpp"

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/functional.h>
#include <thrust/find.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>

using namespace mars;
using namespace mars::fem;
using KeyType  = uint64_t;
using RealType = double;

struct Problem
{
    bool linear = true;
    double alpha = 100, c[3] = {0.3, 0.4, 0.45};

    __host__ __device__ double r2(double x, double y, double z) const
    {
        return (x - c[0]) * (x - c[0]) + (y - c[1]) * (y - c[1]) + (z - c[2]) * (z - c[2]);
    }
    __host__ __device__ double u(double x, double y, double z) const
    {
        return linear ? 1 + x + 2 * y + 3 * z : exp(-alpha * r2(x, y, z));
    }
    // f = -laplace u
    __host__ __device__ double f(double x, double y, double z) const
    {
        if (linear) return 0;
        const double s = r2(x, y, z);
        return (6 * alpha - 4 * alpha * alpha * s) * exp(-alpha * s);
    }
};

struct Options
{
    std::string mesh, element = "hex", marking = "doerfler";
    Problem problem;
    int levels      = 3;
    double radius   = 0.25; // refine the elements whose centroid lies this close to c
    double tol      = 1e-10;
    int maxIter     = 500;
    int bucketSize  = 64;
    bool measure_adapt = false;
    bool diagnose_keys = false;
    bool mesh_only = false;
    bool audit_mesh = false;
    bool audit_algebra = false;
    bool diagnose   = false; // per-rank errors by kind of slot, before and after the solve
};

template<typename KeyT, int N>
struct Corners
{
    const KeyT* c[N];
};

#include "mars_amr_galerkin_audit.hpp"

// Every corner pair of the own elements as a (row << 32 | col) key, then the diagonal of every slot: an owned node
// whose elements are all halos here still needs its row.
template<int N>
__global__ void patternKernel(Corners<KeyType, N> conn, size_t first, size_t own, size_t numSlots,
                              unsigned long long* keys)
{
    const size_t t = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (t < own)
    {
        const size_t e = first + t;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                keys[(t * N + i) * N + j] =
                    (unsigned long long)conn.c[i][e] << 32 | (unsigned long long)conn.c[j][e];
    }
    else if (t < own + numSlots)
    {
        const unsigned long long s = t - own;
        keys[own * N * N + s]      = s << 32 | s;
    }
}

struct PatternCol
{
    __host__ __device__ int operator()(unsigned long long k) const { return int(k & 0xffffffffULL); }
};
struct PatternRow
{
    __host__ __device__ unsigned long long operator()(long long r) const { return (unsigned long long)r << 32; }
};

template<int N>
__device__ inline void cornersOf(const Corners<KeyType, N>& conn, size_t e, const double* x, const double* y,
                                 const double* z, int (&n)[N], double (&X)[N][3])
{
    for (int i = 0; i < N; ++i)
    {
        n[i]    = int(conn.c[i][e]);
        X[i][0] = x[n[i]];
        X[i][1] = y[n[i]];
        X[i][2] = z[n[i]];
    }
}

// Element stiffness and load of the own elements, added into the slot CSR and the slot load vector
template<typename ElementTag, int N>
__global__ void assembleKernel(Corners<KeyType, N> conn, size_t first, size_t own, const double* x, const double* y,
                               const double* z, Problem p, const int* rowPtr, const int* cols, double* values,
                               double* load)
{
    const size_t t = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (t >= own) return;
    int n[N];
    double X[N][3];
    cornersOf(conn, first + t, x, y, z, n, X);
    double K[N][N] = {}, F[N] = {};
    ElementValues<ElementTag, double> v;
    for (int q = 0; q < v.numQuadPoints; ++q)
    {
        v.reinit(X, q);
        const double fq = p.f(v.x[0], v.x[1], v.x[2]);
        for (int i = 0; i < N; ++i)
        {
            F[i] += fq * v.phi[i] * v.JxW;
            for (int j = 0; j < N; ++j)
                K[i][j] += (v.grad[i][0] * v.grad[j][0] + v.grad[i][1] * v.grad[j][1] + v.grad[i][2] * v.grad[j][2]) *
                           v.JxW;
        }
    }
    for (int i = 0; i < N; ++i)
    {
        atomicAdd(&load[n[i]], F[i]);
        for (int j = 0; j < N; ++j)
        {
            int lo = rowPtr[n[i]], hi = rowPtr[n[i] + 1];
            while (lo < hi) // sorted columns
            {
                const int mid = (lo + hi) / 2;
                if (cols[mid] < n[j]) lo = mid + 1;
                else hi = mid;
            }
            atomicAdd(&values[lo], K[i][j]);
        }
    }
}

// (u_h - u)^2 integrated over each own element
template<typename ElementTag, int N>
__global__ void l2ErrorKernel(Corners<KeyType, N> conn, size_t first, size_t own, const double* x, const double* y,
                              const double* z, Problem p, const double* uh, double* err)
{
    const size_t t = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (t >= own) return;
    int n[N];
    double X[N][3];
    cornersOf(conn, first + t, x, y, z, n, X);
    ElementValues<ElementTag, double> v;
    double sum = 0;
    for (int q = 0; q < v.numQuadPoints; ++q)
    {
        v.reinit(X, q);
        double value = 0;
        for (int i = 0; i < N; ++i)
            value += uh[n[i]] * v.phi[i];
        const double d = value - p.u(v.x[0], v.x[1], v.x[2]);
        sum += d * d * v.JxW;
    }
    err[t] = sum;
}

// 1 for elements whose centroid lies within radius of c, else 0: the marking of the next level
template<int N>
__global__ void nearCentreKernel(Corners<KeyType, N> conn, size_t numElements, const double* x, const double* y,
                                 const double* z, Problem p, double radius, double* indicator)
{
    const size_t e = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (e >= numElements) return;
    double m[3] = {0, 0, 0};
    for (int i = 0; i < N; ++i)
    {
        const KeyType s = conn.c[i][e];
        m[0] += x[s] / N;
        m[1] += y[s] / N;
        m[2] += z[s] / N;
    }
    indicator[e] = p.r2(m[0], m[1], m[2]) < radius * radius ? 1.0 : 0.0;
}

struct OnBox
{
    const double *x, *y, *z;
    const uint8_t* isDof;
    double lo[3], hi[3], eps;
    __host__ __device__ double operator()(size_t s) const
    {
        if (!isDof[s]) return 0;
        const double v[3] = {x[s], y[s], z[s]};
        bool on           = false;
        for (int d = 0; d < 3; ++d)
            on |= fabs(v[d] - lo[d]) < eps || fabs(v[d] - hi[d]) < eps;
        return on ? 1.0 : 0.0;
    }
};

struct NodalError
{
    const double *x, *y, *z, *uh;
    const uint8_t* synced;
    Problem p;
    __host__ __device__ double operator()(size_t s) const
    {
        return synced[s] ? fabs(uh[s] - p.u(x[s], y[s], z[s])) : 0.0;
    }
};

// --diagnose: |field - u| at the slots of one kind (0 DOF, 1 copy, 2 hanging) that prolong keeps current, else 0
struct SlotError
{
    const double *x, *y, *z, *field;
    const uint8_t *isDof, *synced, *hanging;
    Problem p;
    int kind;
    __host__ __device__ double operator()(size_t s) const
    {
        const int k = hanging[s] ? 2 : isDof[s] ? 0 : 1;
        return synced[s] && k == kind ? fabs(field[s] - p.u(x[s], y[s], z[s])) : 0.0;
    }
};
struct SlotAbove
{
    SlotError e;
    double tol;
    __host__ __device__ bool operator()(size_t s) const { return e(s) > tol; }
};
struct SlotEquals
{
    SlotError e;
    double value;
    __host__ __device__ bool operator()(size_t s) const { return e(s) == value; }
};

// Per rank and kind of slot: the largest error, how many slots exceed 1e-4, and where the worst one is
void diagnose(const char* stage, int level, int rank, size_t n, SlotError e)
{
    const char* names[3] = {"DOF", "copy", "hanging"};
    for (int kind = 0; kind < 3; ++kind)
    {
        e.kind            = kind;
        auto first        = thrust::counting_iterator<size_t>(0);
        const double max  = n ? mars::transformReduce(first, first + n, e, 0.0, thrust::maximum<double>()) : 0.0;
        const long above  = n ? long(mars::countIf(first, first + n, SlotAbove{e, 1e-4})) : 0;
        double where[3]   = {0, 0, 0};
        if (above > 0)
        {
            const size_t s = *thrust::find_if(thrust::device, first, first + n, SlotEquals{e, max});
            cudaMemcpy(&where[0], e.x + s, sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(&where[1], e.y + s, sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(&where[2], e.z + s, sizeof(double), cudaMemcpyDeviceToHost);
        }
        std::printf("[diag] %s level=%d rank=%d %-7s max=%.3e above1e-4=%ld worst=(%.6f, %.6f, %.6f)\n", stage, level,
                    rank, names[kind], max, above, where[0], where[1], where[2]);
    }
    std::fflush(stdout);
}

// Opt-in check of the node identities before the solve changes any fields.
__global__ void diagnoseKeysKernel(const KeyType* keys, size_t n, amr::KeyedPosition<KeyType, RealType> position,
                                  const uint8_t* owned, int rank, int level, int stage, unsigned* bad)
{
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i >= n) return;
    RealType x, y, z;
    position(keys[i], x, y, z);
    const KeyType encoded = position.encode(x, y, z);
    if (encoded != keys[i] && atomicAdd(bad, 1u) < 8)
        printf("[amr-key-invalid] rank=%d level=%d stage=%d slot=%llu owned=%u key=%llu encoded=%llu\n",
               rank, level, stage, (unsigned long long)i, unsigned(owned[i]),
               (unsigned long long)keys[i], (unsigned long long)encoded);
}

struct PhaseClock
{
    std::vector<std::pair<std::string, double>> phases;
    std::chrono::steady_clock::time_point t0 = now();
    static std::chrono::steady_clock::time_point now()
    {
        cudaDeviceSynchronize();
        MPI_Barrier(MPI_COMM_WORLD);
        return std::chrono::steady_clock::now();
    }
    void lap(const char* name)
    {
        auto t1 = now();
        phases.emplace_back(name, std::chrono::duration<double, std::milli>(t1 - t0).count());
        t0 = t1;
    }
};

template<typename ElementTag>
int run(const Options& opt, int rank, int numRanks)
{
    constexpr int N   = ElementTag::NodesPerElement;
    using Manager     = amr::AmrManager<ElementTag, KeyType, RealType>;
    using Domain      = typename Manager::Domain;
    using Space       = DofSpace<KeyType, RealType, Domain>;
    using Vector      = cstone::DeviceVector<double>;
    const Problem& p  = opt.problem;

    typename Manager::Config config;
    config.strategy       = opt.marking == "octree" ? amr::MarkingStrategy::OctreeNative :
                                                    amr::MarkingStrategy::Doerfler;
    config.refineFraction = 0.5; // the indicator is 0 or 1
    config.maxLevels      = 64;  // the level gate is not used; the loop below stops after opt.levels
    config.bucketSize     = opt.bucketSize;
    Manager amr(config);
    amr.initialize(opt.mesh, rank, numRanks);

    bool exact = true;
    Vector u;
    cstone::DeviceVector<unsigned> key_failures(opt.diagnose_keys ? 1 : 0);
    for (int level = 0; level <= opt.levels; ++level)
    {
        PhaseClock clock;
        Domain& domain = amr.domain();
        domain.cacheNodeCoordinates();
        const size_t n = domain.getNodeCount(), numElements = domain.getElementCount();
        const size_t first = domain.startIndex(), own = domain.localElementCount();
        const double *x = domain.getNodeX().data(), *y = domain.getNodeY().data(), *z = domain.getNodeZ().data();
        const auto& conn = domain.getElementToNodeConnectivity();
        Corners<KeyType, N> c;
        [&]<std::size_t... I>(std::index_sequence<I...>) { ((c.c[I] = std::get<I>(conn).data()), ...); }
        (std::make_index_sequence<N>{});
        clock.lap("domain");
        auto check_keys = [&](int stage) {
            if (!opt.diagnose_keys) return;
            cudaCheckCall(cudaMemsetAsync(key_failures.data(), 0, sizeof(unsigned)));
            if (n)
                diagnoseKeysKernel<<<unsigned((n + 255) / 256), 256>>>(
                    domain.getLocalToGlobalSfcMap().data(), n, {domain.getBoundingBox()},
                    domain.getNodeOwnershipMap().data(), rank, level, stage, key_failures.data());
            cudaCheckError();
            unsigned bad = 0;
            cudaCheckCall(cudaMemcpy(&bad, key_failures.data(), sizeof(bad), cudaMemcpyDeviceToHost));
            MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            if (rank == 0) std::printf("[amr-key-check] level=%d stage=%d invalid=%u\n", level, stage, bad);
            if (bad) MPI_Abort(MPI_COMM_WORLD, 79);
        };
        check_keys(0); // rebuilt domain

        auto adapt = [&](bool transfer) {
            Vector indicator(numElements > 0 ? numElements : 1, 0.0);
            if (numElements > 0)
                nearCentreKernel<N><<<unsigned((numElements + 255) / 256), 256>>>(
                    c, numElements, x, y, z, p, opt.radius, indicator.data());
            cudaCheckError();
            Vector next;
            if (opt.measure_adapt)
            {
                checkGpuErrors(cudaDeviceSynchronize());
                MPI_Barrier(MPI_COMM_WORLD);
            }
            const auto adapt_start = std::chrono::steady_clock::now();
            auto stats = transfer ? amr.adaptMesh(indicator.data(), u.data(), next) :
                                    amr.adaptMeshMultiField(indicator.data(), {}, {});
            if (rank == 0)
                std::printf("    adapt: %zu elements refined, %.2f ms\n", size_t(stats.elementsRefined), stats.totalTimeMs);
            if (transfer) u = std::move(next);
            if (opt.measure_adapt)
            {
                checkGpuErrors(cudaDeviceSynchronize());
                double elapsed = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - adapt_start).count();
                MPI_Allreduce(MPI_IN_PLACE, &elapsed, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                if (rank == 0) std::printf("[amr-adapt-wall] level=%d ms=%.9g\n", level, elapsed);
            }
        };
        if (opt.mesh_only)
        {
            uint64_t elements = own;
            MPI_Allreduce(MPI_IN_PLACE, &elements, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
            if (rank == 0) std::printf("[amr-mesh-only] level=%d elements=%llu\n", level, (unsigned long long)elements);
            if (level < opt.levels) adapt(false);
            continue;
        }

        HangingNodes<KeyType, RealType> hanging(domain);
        clock.lap("hanging nodes");
        Space space(domain, nullptr, &hanging);
        const uint8_t* isDof = space.isDof();
        clock.lap("DOF space");
        check_keys(1);
        if (opt.audit_algebra)
        {
            const auto& box = domain.getBoundingBox();
            const auto bits = box.getBoxDimBits(cstone::maxTreeLevel<KeyType>{});
            const double cell_size = std::max({box.lx() / double(uint64_t(1) << bits[0]),
                box.ly() / double(uint64_t(1) << bits[1]), box.lz() / double(uint64_t(1) << bits[2])});
            amr_audit::check_prolongation(space, c, first, own, n, domain.getLocalToGlobalSfcMap().data(),
                                          x, y, z, cell_size, rank, level);
        }


        thrust::device_vector<uint8_t> hangingSlot;
        SlotError slotError{x, y, z, nullptr, isDof, space.isSynced(), nullptr, p, 0};
        if (opt.diagnose)
        {
            // P of the exact values: the constraints, copies and extra values alone, no solve
            hangingSlot.assign(n, 0);
            auto hangingSlots = thrust::device_pointer_cast(hanging.slots());
            thrust::fill(thrust::device, thrust::make_permutation_iterator(hangingSlot.begin(), hangingSlots),
                         thrust::make_permutation_iterator(hangingSlot.begin(), hangingSlots + hanging.size()),
                         uint8_t(1));
            slotError.hanging = thrust::raw_pointer_cast(hangingSlot.data());
            Vector exactAtDofs(n);
            thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                              exactAtDofs.data(), [=] __host__ __device__(size_t s) {
                                  return isDof[s] ? p.u(x[s], y[s], z[s]) : 0.0;
                              });
            space.prolong(exactAtDofs);
            slotError.field = exactAtDofs.data();
            diagnose("prolong", level, rank, n, slotError);
            clock.lap("diagnose");
        }

        // DOF numbering: owned DOFs of rank r follow those of ranks < r
        thrust::device_vector<int> dofIndex(n);
        thrust::transform(thrust::device, isDof, isDof + n, dofIndex.begin(), thrust::identity<int>());
        thrust::exclusive_scan(thrust::device, dofIndex.begin(), dofIndex.end(), dofIndex.begin(), 0);
        const long long numOwned = mars::count(isDof, isDof + n, uint8_t(1));
        long long dofStart       = 0;
        MPI_Exscan(&numOwned, &dofStart, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (rank == 0) dofStart = 0;
        const int* idx = thrust::raw_pointer_cast(dofIndex.data());

        // The global id of every slot's DOF, sent to the copies as a field
        Vector gid(n);
        thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                          gid.data(), [=] __host__ __device__(size_t s) {
                              return isDof[s] ? double(dofStart + idx[s]) : -1.0;
                          });
        space.prolong(gid);
        thrust::device_vector<long long> slotGid(n);
        thrust::device_vector<int> dofRow(n);
        const double* g = gid.data();
        thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                          slotGid.begin(), [=] __host__ __device__(size_t s) { return llround(g[s]); });
        thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                          dofRow.begin(), [=] __host__ __device__(size_t s) { return isDof[s] ? idx[s] : -1; });

        // Slot CSR over the own elements
        const size_t numKeys = own * N * N + n;
        thrust::device_vector<unsigned long long> keys(numKeys);
        patternKernel<N><<<unsigned((own + n + 255) / 256), 256>>>(c, first, own, n,
                                                                    thrust::raw_pointer_cast(keys.data()));
        cudaCheckError();
        thrust::sort(thrust::device, keys.begin(), keys.end());
        const size_t nnz = thrust::unique(thrust::device, keys.begin(), keys.end()) - keys.begin();
        thrust::device_vector<int> rowPtr(n + 1), cols(nnz);
        thrust::transform(thrust::device, keys.begin(), keys.begin() + nnz, cols.begin(), PatternCol{});
        auto rowKeys = thrust::make_transform_iterator(thrust::counting_iterator<long long>(0), PatternRow{});
        thrust::lower_bound(thrust::device, keys.begin(), keys.begin() + nnz, rowKeys, rowKeys + n + 1,
                            rowPtr.begin());
        thrust::device_vector<double> values(nnz, 0.0);
        Vector load(n, 0.0);
        if (own > 0)
            assembleKernel<ElementTag, N><<<unsigned((own + 255) / 256), 256>>>(
                c, first, own, x, y, z, p, thrust::raw_pointer_cast(rowPtr.data()),
                thrust::raw_pointer_cast(cols.data()), thrust::raw_pointer_cast(values.data()), load.data());
        cudaCheckError();
        clock.lap("assembly");

        // Dirichlet DOFs on the faces of the unit cube, and u_D: the boundary values there, zero elsewhere
        double lo[3], hi[3];
        const double* xyz[3] = {x, y, z};
        for (int d = 0; d < 3; ++d)
        {
            auto b = thrust::device_pointer_cast(xyz[d]);
            lo[d]  = mars::reduce(b, b + n, 1e30, thrust::minimum<double>());
            hi[d]  = mars::reduce(b, b + n, -1e30, thrust::maximum<double>());
        }
        MPI_Allreduce(MPI_IN_PLACE, lo, 3, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, hi, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        OnBox onBox{x, y, z, isDof, {lo[0], lo[1], lo[2]}, {hi[0], hi[1], hi[2]}, 1e-8};
        Vector removed(n), uD(n);
        thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                          removed.data(), onBox);
        const double* rm = removed.data();
        thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                          uD.data(), [=] __host__ __device__(size_t s) {
                              return rm[s] != 0 ? p.u(x[s], y[s], z[s]) : 0.0;
                          });
        space.prolong(uD);

        // b = P^T (F - A P u_D) on the free DOFs, 0 on the Dirichlet ones
        Vector residual(n);
        {
            const int* rp = thrust::raw_pointer_cast(rowPtr.data());
            const int* cc = thrust::raw_pointer_cast(cols.data());
            const double* vv = thrust::raw_pointer_cast(values.data());
            const double *f = load.data(), *ud = uD.data();
            thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0),
                              thrust::counting_iterator<size_t>(n), residual.data(),
                              [=] __host__ __device__(size_t s) {
                                  double r = f[s];
                                  for (int j = rp[s]; j < rp[s + 1]; ++j)
                                      r -= vv[j] * ud[cc[j]];
                                  return r;
                              });
        }
        space.restrict(residual);
        thrust::device_vector<double> rhs(numOwned), solution(numOwned, 0.0);
        {
            const double* r = residual.data();
            double* b       = thrust::raw_pointer_cast(rhs.data());
            thrust::for_each(thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(n),
                             [=] __host__ __device__(size_t s) {
                                 if (isDof[s]) b[idx[s]] = rm[s] != 0 ? 0.0 : r[s];
                             });
        }

        thrust::device_vector<int> rows;
        thrust::device_vector<long long> colGids;
        thrust::device_vector<double> entries;
        space.restrictMatrix(thrust::raw_pointer_cast(rowPtr.data()), thrust::raw_pointer_cast(cols.data()),
                             thrust::raw_pointer_cast(values.data()), thrust::raw_pointer_cast(slotGid.data()),
                             thrust::raw_pointer_cast(dofRow.data()), rows, colGids, entries, removed.data());
        HypreMatrix A = hypreFromEntries(MPI_COMM_WORLD, dofStart, HYPRE_Int(numOwned), space.numDofs(), rows,
                                         colGids, entries);
        clock.lap("P^T A P");
        check_keys(2);
        HypreAmgPcgSolver solver("MARS_AMR_AMG");
        solver.setup(MPI_COMM_WORLD, std::move(A), opt.tol, opt.maxIter);
        clock.lap("AMG setup");
        if (opt.audit_algebra)
            amr_audit::check_operator(space, solver, n, numOwned, domain.getLocalToGlobalSfcMap().data(), idx,
                removed.data(), thrust::raw_pointer_cast(rowPtr.data()), thrust::raw_pointer_cast(cols.data()),
                thrust::raw_pointer_cast(values.data()), thrust::raw_pointer_cast(rhs.data()), x, y, z, p, rank, level);

        const int iterations = solver.solve(thrust::raw_pointer_cast(rhs.data()),
                                            thrust::raw_pointer_cast(solution.data()));
        clock.lap("PCG");
        if (opt.audit_algebra)
            amr_audit::check_solution(space, solver, n, numOwned, idx, removed.data(),
                thrust::raw_pointer_cast(rowPtr.data()), thrust::raw_pointer_cast(cols.data()),
                thrust::raw_pointer_cast(values.data()), thrust::raw_pointer_cast(rhs.data()),
                thrust::raw_pointer_cast(solution.data()), opt.tol, rank, level);


        // u = u_D + solution on the DOFs, then on every slot
        u.resize(n);
        {
            const double* sol = thrust::raw_pointer_cast(solution.data());
            const double* ud  = uD.data();
            thrust::transform(thrust::device, thrust::counting_iterator<size_t>(0),
                              thrust::counting_iterator<size_t>(n), u.data(), [=] __host__ __device__(size_t s) {
                                  return isDof[s] ? ud[s] + sol[idx[s]] : 0.0;
                              });
        }
        space.prolong(u);
        if (opt.diagnose)
        {
            slotError.field = u.data();
            diagnose("solution", level, rank, n, slotError);
        }

        double maxErr = mars::transformReduce(thrust::counting_iterator<size_t>(0),
                                              thrust::counting_iterator<size_t>(n),
                                              NodalError{x, y, z, u.data(), space.isSynced(), p}, 0.0,
                                              thrust::maximum<double>());
        Vector elementErr(own > 0 ? own : 1, 0.0);
        if (own > 0)
            l2ErrorKernel<ElementTag, N><<<unsigned((own + 255) / 256), 256>>>(c, first, own, x, y, z, p, u.data(),
                                                                              elementErr.data());
        cudaCheckError();
        double l2 = own > 0 ? mars::reduce(thrust::device_pointer_cast(elementErr.data()),
                                           thrust::device_pointer_cast(elementErr.data()) + own, 0.0,
                                           thrust::plus<double>())
                            : 0.0;
        long long elements = (long long)own;
        MPI_Allreduce(MPI_IN_PLACE, &maxErr, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &l2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &elements, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        l2 = std::sqrt(l2);
        clock.lap("errors");
        check_keys(3);
        if (iterations < 0 || !std::isfinite(l2) || !std::isfinite(maxErr)) exact = false;
        // The linear field is in the discrete space up to where the nodes sit: positions come from SFC keys, about
        // 2^-21 of the box apart, so a hanging node is up to a key cell off the mean of the nodes it hangs on (a few
        // 1e-6 here). A missing or wrong constraint is off by about h |grad u|, 1e-2 and more.
        if (p.linear && maxErr > 1e-4) exact = false;

        if (rank == 0)
        {
            std::printf("[amr-galerkin] level=%d elements=%lld dofs=%lld hanging=%lld iterations=%d relres=%.3e "
                        "l2=%.6e max=%.6e\n",
                        level, elements, space.numDofs(), hanging.globalCount(), iterations,
                        solver.lastRelativeResidual(), l2, maxErr);
            for (const auto& [name, ms] : clock.phases)
                std::printf("    %-16s %10.2f ms\n", name.c_str(), ms);
            std::fflush(stdout);
        }

        if (level == opt.levels) break;
        adapt(true);
    }
    if (rank == 0) std::printf("[amr-galerkin] %s\n", exact ? "PASS" : "FAIL");
    return exact ? 0 : 1;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) cudaSetDevice(rank % deviceCount);

    char host[MPI_MAX_PROCESSOR_NAME], pci[32];
    int host_length = 0, device = 0;
    MPI_Get_processor_name(host, &host_length);
    checkGpuErrors(cudaGetDevice(&device));
#if defined(__HIPCC__)
    checkGpuErrors(hipDeviceGetPCIBusId(pci, sizeof(pci), device));
#else
    checkGpuErrors(cudaDeviceGetPCIBusId(pci, sizeof(pci), device));
#endif
    std::printf("[amr-device] rank=%d host=%s pci=%s\n", rank, host, pci);
    std::fflush(stdout);

    Options opt;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        auto value          = [&](const char* key) { return a.substr(std::string(key).size()); };
        if (a.rfind("--mesh=", 0) == 0) opt.mesh = value("--mesh=");
        else if (a.rfind("--element=", 0) == 0) opt.element = value("--element=");
        else if (a == "--problem=linear") opt.problem.linear = true;
        else if (a == "--problem=peak") opt.problem.linear = false;
        else if (a.rfind("--levels=", 0) == 0) opt.levels = std::stoi(value("--levels="));
        else if (a.rfind("--radius=", 0) == 0) opt.radius = std::stod(value("--radius="));
        else if (a.rfind("--alpha=", 0) == 0) opt.problem.alpha = std::stod(value("--alpha="));
        else if (a.rfind("--tol=", 0) == 0) opt.tol = std::stod(value("--tol="));
        else if (a.rfind("--max-iter=", 0) == 0) opt.maxIter = std::stoi(value("--max-iter="));
        else if (a.rfind("--bucket-size=", 0) == 0) opt.bucketSize = std::stoi(value("--bucket-size="));
        else if (a.rfind("--marking=", 0) == 0) opt.marking = value("--marking=");
        else if (a == "--measure-adapt") opt.measure_adapt = true;
        else if (a == "--diagnose-keys") opt.diagnose_keys = true;
        else if (a == "--diagnose") opt.diagnose = true;
        else if (a == "--mesh-only") opt.mesh_only = true;
        else if (a == "--audit-mesh") opt.audit_mesh = true;
        else if (a == "--audit-algebra") opt.audit_algebra = true;
        else
        {
            if (rank == 0) std::fprintf(stderr, "unknown option %s\n", a.c_str());
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
    }
    if ((opt.marking != "doerfler" && opt.marking != "octree") || opt.mesh.empty() || (opt.element != "hex" && opt.element != "tet") ||
        (opt.audit_mesh && opt.element != "tet") || (opt.mesh_only && opt.audit_algebra))
    {
        if (rank == 0)
            std::fprintf(stderr, "usage: %s --mesh=DIR --element=hex|tet [--problem=linear|peak] [--levels=N] "
                                 "[--radius=R] [--alpha=A] [--tol=T] [--max-iter=N] [--bucket-size=N] [--diagnose] [--mesh-only] [--audit-mesh (tet)] [--audit-algebra] [--marking=doerfler|octree] [--measure-adapt] [--diagnose-keys]\n",
                         argv[0]);
        MPI_Finalize();
        return 2;
    }
    if (opt.audit_mesh) setenv("MARS_AMR_AUDIT_TET", "1", 1);
    if (rank == 0)
        std::printf("[amr-galerkin] %s, %s, %d levels, %d ranks, refinement within %g of (%g, %g, %g)\n",
                    opt.element.c_str(), opt.problem.linear ? "u = 1 + x + 2y + 3z" : "peak", opt.levels, numRanks,
                    opt.radius, opt.problem.c[0], opt.problem.c[1], opt.problem.c[2]);

    int status = opt.element == "hex" ? run<HexTag>(opt, rank, numRanks) : run<TetTag>(opt, rank, numRanks);
    MPI_Finalize();
    return status;
}
