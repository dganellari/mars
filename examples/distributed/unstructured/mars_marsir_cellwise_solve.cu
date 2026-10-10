// Element-local (cell-wise) solve of the high-order CVFEM Laplacian on one or more
// GPUs: the mixed form of the plan from internal-notes/FlexibleCG.pdf. The operator is the
// kernel MARSIR generates from high-level MLIR (a PTX file, loaded at run time);
// the DSS and the Krylov solver are hand-written
// (backend/distributed/unstructured/solvers/mars_cellwise_krylov.hpp). No assembled
// vector, index map, atomic or coloring appears anywhere.
//
// Mesh: the unit cube as ne^3 trilinear hexahedra, optionally deformed by
// x += a sin(pi x) sin(pi y) sin(pi z) (1,1,1), which keeps the boundary in place and
// gives every element its own metric. The manufactured solution
// u = sin(pi x) sin(pi y) sin(pi z) at the nodes gives the RHS b = A u, so the
// discrete solution is u itself and the error measures the solver alone.
//
// All per-element work runs on the GPU: corners, metric, diagonal, u, b. On several
// MPI ranks each rank owns a sub-block of the ne^3 block (one GPU per rank, chosen by
// the rank on its node); the run needs GPU-aware MPI (on Cray MPICH,
// MPICH_GPU_SUPPORT_ENABLED=1). Any rank count gives the same iterations as one rank,
// up to the rounding of the global sums.
//
// Run:  mars_marsir_cellwise_solve --ptx <hl_full_p7_sm90.ptx> [--ne 8] [--deform 0]
//                                  [--tol 1e-10] [--maxit 500] [--reps 10]
//                                  [--mg] [--pre 0] [--post 3]
// --unstructured runs the same cube through the unstructured tables
// (solvers/mars_cellwise_topology.hpp) instead of the structured block, on any number of
// ranks (each with the ghost elements around its sub-block), and --rotate gives every
// element a random proper rotation of its local frame (Jacobi only). The discrete
// problem does not change, so the iterations must match the structured run.
// --mg preconditions with the geometric multigrid V-cycle instead of DSS + Jacobi
// (ne must be a power of two), with --pre / --post Chebyshev steps around the coarse
// correction; its history must match marsir-mlir/test/cellwise_multigrid_ref.py
// with the same options.
// With the same --ne and --deform, the residual history must match
// marsir-mlir/test/cellwise_krylov_ref.py to about 4 digits (BiCGStab amplifies the
// rounding of the different summation orders).

#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"
#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_operator.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_multigrid.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"

#include <cuda_runtime.h>
#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

using namespace mars;
using mars::cellwise::kN;
using mars::cellwise::kN3;
using mars::cellwise::kNN;

namespace {

constexpr int kP = 7;
constexpr double kPi = 3.14159265358979323846;
constexpr int kGElem = 3 * kP * 3 * kNN;   // metric doubles per element, component-major

__constant__ int c_rot[24][9];   // proper rotations of the cube; 0 is the identity

__device__ unsigned long long mix(unsigned long long z)
{
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

// Rotation of the local frame of global element g (lattice order of the whole block):
// local corner c sits at the block corner whose signs are R s_c. By global index, so
// every element has the same frame on any number of ranks.
__device__ int rotation_of(long long g, bool rotate) { return rotate ? (int)(mix(g) % 24) : 0; }

__device__ void corner_offset(int rot, int c, int (&b)[3])
{
    int cb[3];
    cellwise::corner_bits(c, cb[0], cb[1], cb[2]);
    for (int i = 0; i < 3; ++i) {
        int s = 0;
        for (int j = 0; j < 3; ++j) s += c_rot[rot][i * 3 + j] * (2 * cb[j] - 1);
        b[i] = s > 0;
    }
}

// Vertex (i, j, k) of the global vertex lattice of block `blk` on the unit cube,
// deformed in place. Coarser levels see the same deformed points.
__device__ void vertex(int i, int j, int k, const cellwise::Block& blk, double deform, double x[3])
{
    x[0] = (double)i / blk.NX;
    x[1] = (double)j / blk.NY;
    x[2] = (double)k / blk.NZ;
    const double bump = deform * sin(kPi * x[0]) * sin(kPi * x[1]) * sin(kPi * x[2]);
    x[0] += bump;
    x[1] += bump;
    x[2] += bump;
}

// The 8 corners of local element e in its own corner order (c_hexCornerRef), with the
// element's frame rotated when `rotate`.
__device__ void element_corners(long long e, const cellwise::Block& blk, double deform, bool rotate,
                                double corners[8][3])
{
    const int ez = (int)(e % blk.nz) + blk.oz;
    const int ey = (int)((e / blk.nz) % blk.ny) + blk.oy;
    const int ex = (int)(e / ((long long)blk.nz * blk.ny)) + blk.ox;
    const int rot = rotation_of(((long long)ex * blk.NY + ey) * blk.NZ + ez, rotate);
    for (int c = 0; c < 8; ++c) {
        int b[3];
        corner_offset(rot, c, b);
        vertex(ex + b[0], ey + b[1], ez + b[2], blk, deform, corners[c]);
    }
}

// Corner keys and local ids of the unstructured tables for the block's elements
// elem[0 .. E) (global lattice order): global lattice vertex ids.
__global__ void mesh_kernel(cellwise::Block blk, bool rotate, const unsigned long long* elem, long long E,
                            unsigned long long* const* key, int* const* lid)
{
    const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= E) return;
    const long long g = (long long)elem[i];
    const int ez = (int)(g % blk.NZ), ey = (int)((g / blk.NZ) % blk.NY);
    const int ex = (int)(g / ((long long)blk.NZ * blk.NY));
    const int rot = rotation_of(g, rotate);
    for (int c = 0; c < 8; ++c) {
        int b[3];
        corner_offset(rot, c, b);
        const long long v = (((long long)(ex + b[0]) * (blk.NY + 1)) + ey + b[1]) * (blk.NZ + 1) + ez + b[2];
        key[c][i] = (unsigned long long)v;
        lid[c][i] = (int)v;
    }
}

// One thread per metric point (e, dir, l, s, r): MARS's device metric, written in
// the component-major layout the MARSIR kernel reads, [dir][face][comp][s][r].
__global__ void metric_kernel(double* __restrict__ G, cellwise::Block blk, double deform, bool rotate)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long points = 3LL * kP * kNN;
    if (t >= blk.elements() * points) return;
    const long long e = t / points;
    int q = (int)(t % points);
    const int r = q % kN;
    q /= kN;
    const int s = q % kN;
    q /= kN;
    const int l = q % kP;
    const int dir = q / kP;
    double corners[8][3];
    element_corners(e, blk, deform, rotate, corners);
    double g[3];
    fem::ho_cvfem_metric_point<kP>(corners, dir, l, s, r, g);
    double* out = G + e * kGElem + (dir * kP + l) * 3 * kNN + s * kN + r;
    out[0] = g[0];
    out[kNN] = g[1];
    out[2 * kNN] = g[2];
}

// One thread per node: the diagonal of the element operator in closed form. For a
// unit vector at node (A, B, C) in direction-d coordinates, only face A-1 (plus) and
// face A (minus) see it, and each face's W (x) W smoothing returns
//   W[B][B] W[C][C] g2 Dtil[l][A]
//   + W[B][B] Btil[l][A] sum_r W[C][r] g0[B][r] D[r][C]
//   + W[C][C] Btil[l][A] sum_s W[B][s] g1[s][C] D[s][B]
// at the node itself. Checked against the assembled element matrix (2.9e-16).
__global__ void diagonal_kernel(const double* __restrict__ G, double* __restrict__ diag,
                                long long E)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * kN3) return;
    const long long e = t / kN3;
    const int node = (int)(t % kN3);
    const int a = node / kNN, b = (node / kN) % kN, c = node % kN;
    const double* Ge = G + e * kGElem;
    double sum = 0.0;
    for (int dir = 0; dir < 3; ++dir) {
        const int A = dir == 0 ? a : dir == 1 ? b : c;
        const int B = dir == 0 ? b : a;
        const int C = dir == 2 ? b : c;
        for (int side = 0; side < 2; ++side) {
            const int l = side == 0 ? A - 1 : A;
            if (l < 0 || l >= kP) continue;
            const double* g0 = Ge + (dir * kP + l) * 3 * kNN;
            const double* g1 = g0 + kNN;
            const double* g2 = g0 + 2 * kNN;
            const double bt = fem::c_Btil[l * kN + A];
            double term = fem::c_W[B * kN + B] * fem::c_W[C * kN + C] * g2[B * kN + C] *
                          fem::c_Dtil[l * kN + A];
            double t0 = 0.0, t1 = 0.0;
            for (int q = 0; q < kN; ++q) {
                t0 += fem::c_W[C * kN + q] * g0[B * kN + q] * fem::c_D[q * kN + C];
                t1 += fem::c_W[B * kN + q] * g1[q * kN + C] * fem::c_D[q * kN + B];
            }
            term += fem::c_W[B * kN + B] * bt * t0 + fem::c_W[C * kN + C] * bt * t1;
            sum += side == 0 ? term : -term;
        }
    }
    diag[t] = sum;
}

// u = sin(pi x) sin(pi y) sin(pi z) at every element-local node, through the
// element's trilinear map.
__global__ void exact_kernel(double* __restrict__ u, cellwise::Block blk, double deform, bool rotate)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= blk.values()) return;
    const long long e = t / kN3;
    const int node = (int)(t % kN3);
    const double rf[3] = {fem::c_zeta[node / kNN], fem::c_zeta[(node / kN) % kN],
                          fem::c_zeta[node % kN]};
    double corners[8][3];
    element_corners(e, blk, deform, rotate, corners);
    double x[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 8; ++k) {
        double w = 0.125;
        for (int d = 0; d < 3; ++d) w *= 1.0 + fem::c_hexCornerRef[k][d] * rf[d];
        for (int d = 0; d < 3; ++d) x[d] += w * corners[k][d];
    }
    u[t] = sin(kPi * x[0]) * sin(kPi * x[1]) * sin(kPi * x[2]);
}

__global__ void difference_kernel(double* __restrict__ out, const double* __restrict__ a,
                                  const double* __restrict__ b, long long n)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) out[t] = a[t] - b[t];
}

unsigned blocks_for(long long n, int threads) { return (unsigned)((n + threads - 1) / threads); }

double* device_array(long long n)
{
    double* p = nullptr;
    MARS_CELLWISE_CK(cudaMalloc(&p, n * sizeof(double)));
    return p;
}

double* device_copy(const std::vector<double>& h)
{
    double* p = device_array((long long)h.size());
    MARS_CELLWISE_CK(cudaMemcpy(p, h.data(), h.size() * sizeof(double), cudaMemcpyHostToDevice));
    return p;
}

}  // namespace

int main(int argc, char** argv)
{
    MARS_CELLWISE_MPI(MPI_Init(&argc, &argv));
    int rank = 0;
    MARS_CELLWISE_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    // One GPU per rank, by the rank's index on its node, before any CUDA object exists
    // (the reductions and the PTX module bind to the current device when built).
    {
        MPI_Comm node;
        MARS_CELLWISE_MPI(MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank,
                                              MPI_INFO_NULL, &node));
        int local = 0, devices = 0;
        MARS_CELLWISE_MPI(MPI_Comm_rank(node, &local));
        MARS_CELLWISE_MPI(MPI_Comm_free(&node));
        MARS_CELLWISE_CK(cudaGetDeviceCount(&devices));
        MARS_CELLWISE_CK(cudaSetDevice(local % devices));
        MARS_CELLWISE_CK(cudaFree(nullptr));
    }

    std::string ptx_path;
    int ne = 8, max_iterations = 500, reps = 10;
    double deform = 0.0, tol = 1e-10;
    bool use_mg = false, unstructured = false, rotate = false;
    int pre = 0, post = 3;
    const char* error = nullptr;
    for (int i = 1; i < argc && !error; i += 2) {
        if (!strcmp(argv[i], "--mg")) { use_mg = true; --i; continue; }
        if (!strcmp(argv[i], "--unstructured")) { unstructured = true; --i; continue; }
        if (!strcmp(argv[i], "--rotate")) { rotate = unstructured = true; --i; continue; }
        if (i + 1 >= argc) { error = "an option needs a value"; break; }
        if (!strcmp(argv[i], "--ptx")) ptx_path = argv[i + 1];
        else if (!strcmp(argv[i], "--ne")) ne = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--deform")) deform = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "--tol")) tol = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "--maxit")) max_iterations = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--reps")) reps = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--pre")) pre = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--post")) post = atoi(argv[i + 1]);
        else error = "unknown option";
    }
    if (!error && ptx_path.empty()) error = "--ptx is required";
    if (!error && use_mg && (ne & (ne - 1))) error = "--mg needs ne to be a power of two";
    if (!error && unstructured && use_mg) error = "--unstructured runs with Jacobi only for now";
    if (error) {
        if (rank == 0)
            fprintf(stderr, "%s\nusage: %s --ptx <hl_full_p7_sm90.ptx> [--ne N] [--deform a] [--tol t] "
                            "[--maxit k] [--reps r] [--mg] [--pre k] [--post k] [--unstructured] "
                            "[--rotate]\n", error, argv[0]);
        MPI_Finalize();
        return 1;
    }

    // The 1D basis matrices (at most 64 numbers each) are built once on the host,
    // as for every MARS high-order kernel; all per-element work below is on the GPU.
    const fem::HoCvfemOperators ops = fem::buildHoCvfemOperators(kP);
    MARS_CELLWISE_CK(fem::ho_cvfem_upload_operators(kP, ops.Btil.data(), ops.Dtil.data(),
                                                    ops.D.data(), ops.W.data(), ops.xi.data(),
                                                    ops.zeta.data()));
    // The MARSIR kernel takes Btil/Dtil with the face dimension padded to 8 rows.
    std::vector<double> btil8(kNN, 0.0), dtil8(kNN, 0.0);
    std::copy(ops.Btil.begin(), ops.Btil.end(), btil8.begin());
    std::copy(ops.Dtil.begin(), ops.Dtil.end(), dtil8.begin());
    double* d_btil = device_copy(btil8);
    double* d_dtil = device_copy(dtil8);
    double* d_w = device_copy(ops.W);
    double* d_d = device_copy(ops.D);
    // Pageable uploads may still be in flight when the calls return, and the solver's
    // stream below does not wait for the default stream.
    MARS_CELLWISE_CK(cudaDeviceSynchronize());

    const cellwise::Decomposition dec(MPI_COMM_WORLD, ne, ne, ne);
    const cellwise::Block blk = dec.local;
    // The ghost exchange and the coarse gather pass device pointers to MPI. MPICH (Cray
    // MPICH on Alps) only accepts them with MPICH_GPU_SUPPORT_ENABLED=1; without it the
    // first exchange dies with an unreadable error, so stop here instead.
    if (dec.size > 1) {
        char version[MPI_MAX_LIBRARY_VERSION_STRING];
        int len = 0;
        MARS_CELLWISE_MPI(MPI_Get_library_version(version, &len));
        const char* gpu = getenv("MPICH_GPU_SUPPORT_ENABLED");
        if (strstr(version, "MPICH") && !(gpu && !strcmp(gpu, "1"))) {
            if (rank == 0)
                fprintf(stderr, "MPICH_GPU_SUPPORT_ENABLED=1 is required on several ranks (GPU-aware MPI)\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    // The 24 proper rotations of the cube (signed permutation matrices, determinant +1),
    // identity first.
    {
        std::vector<int> rot;
        const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
        for (const auto& pm : perms)
            for (int sg = 0; sg < 8; ++sg) {
                int M[9] = {0};
                for (int i = 0; i < 3; ++i) M[i * 3 + pm[i]] = (sg >> i & 1) ? -1 : 1;
                const int det = M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6]) +
                                M[2] * (M[3] * M[7] - M[4] * M[6]);
                if (det == 1) rot.insert(rot.end(), M, M + 9);
            }
        MARS_CELLWISE_CK(cudaMemcpyToSymbol(c_rot, rot.data(), rot.size() * sizeof(int)));
        MARS_CELLWISE_CK(cudaDeviceSynchronize());
    }
    const long long E = blk.elements(), n = blk.values();
    const long long unique = (long long)(ne * kP + 1) * (ne * kP + 1) * (ne * kP + 1);
    constexpr int threads = 256;
    cudaStream_t stream;
    MARS_CELLWISE_CK(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));

    double* d_G = device_array(E * kGElem);
    double* d_diag = device_array(n);
    double* d_uex = device_array(n);
    double* d_b = device_array(n);
    double* d_x = device_array(n);
    std::unique_ptr<cellwise::Halo> halo(dec.size > 1 && !unstructured ? new cellwise::Halo(dec, blk) : nullptr);
    std::unique_ptr<cellwise::UnstructuredHalo> uhalo;

    metric_kernel<<<blocks_for(E * 3LL * kP * kNN, threads), threads, 0, stream>>>(d_G, blk, deform, rotate);
    diagonal_kernel<<<blocks_for(n, threads), threads, 0, stream>>>(d_G, d_b, E);   // d_b: scratch until b = A u
    exact_kernel<<<blocks_for(n, threads), threads, 0, stream>>>(d_uex, blk, deform, rotate);
    MARS_CELLWISE_CK(cudaGetLastError());
    // The unstructured tables of the same cube: this rank's elements, then the ghost
    // elements around them; corner keys = global lattice vertex ids.
    cellwise::UnstructuredTopology topo;
    const cellwise::Block flat(E, 1, 1);   // unstructured iteration space: E elements in a row
    if (unstructured) {
        thrust::device_vector<unsigned long long> elem;
        thrust::device_vector<int> owner;
        cellwise::block_elements(dec, elem, owner);
        const long long ET = (long long)elem.size();
        thrust::device_vector<unsigned long long> key[8];
        thrust::device_vector<int> lid[8];
        std::vector<unsigned long long*> kp(8);
        std::vector<int*> lp(8);
        for (int c = 0; c < 8; ++c) {
            key[c].resize(ET);
            lid[c].resize(ET);
            kp[c] = thrust::raw_pointer_cast(key[c].data());
            lp[c] = thrust::raw_pointer_cast(lid[c].data());
        }
        thrust::device_vector<unsigned long long*> d_kp(kp.begin(), kp.end());
        thrust::device_vector<int*> d_lp(lp.begin(), lp.end());
        mesh_kernel<<<blocks_for(ET, threads), threads, 0, stream>>>(blk, rotate, thrust::raw_pointer_cast(elem.data()),
                                                                     ET, thrust::raw_pointer_cast(d_kp.data()),
                                                                     thrust::raw_pointer_cast(d_lp.data()));
        MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
        const std::vector<const unsigned long long*> ckp(kp.begin(), kp.end());
        const std::vector<const int*> clp(lp.begin(), lp.end());
        thrust::device_vector<const unsigned long long*> d_ckp(ckp.begin(), ckp.end());
        thrust::device_vector<const int*> d_clp(clp.begin(), clp.end());
        topo = cellwise::build_topology(thrust::raw_pointer_cast(d_ckp.data()), thrust::raw_pointer_cast(d_clp.data()),
                                        ET, E);
        if (dec.size > 1)
            uhalo.reset(new cellwise::UnstructuredHalo(topo, thrust::raw_pointer_cast(elem.data()),
                                                       thrust::raw_pointer_cast(owner.data()), MPI_COMM_WORLD));
        cellwise::dss(d_b, d_diag, topo, uhalo.get(), stream);   // assembled diagonal on every copy
    } else {
        cellwise::dss(d_b, d_diag, blk, halo.get(), stream);
    }
    const cellwise::TopologyView tview = cellwise::view(topo);
    const cellwise::CountingWeight counting{thrust::raw_pointer_cast(topo.counting.data())};

    const marsir::LaplacianPtx op(ptx_path, kP);
    auto apply = [&](const double* u, double* y) {
        op.apply(u, d_btil, d_dtil, d_w, d_d, d_G, y, E, stream);
    };
    apply(d_uex, d_b);                   // b = A u, unassembled

    // Multigrid levels: the metric and the unassembled diagonal of any block (a rank's
    // part or, on rank 0, a whole coarse block), made on the GPU.
    auto apply_on = [&](const double* u, double* y, const double* G, long long elements,
                        cudaStream_t s) { op.apply(u, d_btil, d_dtil, d_w, d_d, G, y, elements, s); };
    auto make_level = [&](const cellwise::Block& b) {
        double* g = device_array(b.elements() * kGElem);
        double* raw = device_array(b.values());
        metric_kernel<<<blocks_for(b.elements() * 3LL * kP * kNN, threads), threads, 0, stream>>>(g, b, deform, false);
        diagonal_kernel<<<blocks_for(b.values(), threads), threads, 0, stream>>>(g, raw, b.elements());
        MARS_CELLWISE_CK(cudaGetLastError());
        return std::pair<double*, double*>(g, raw);
    };
    using Multigrid = cellwise::Multigrid<decltype(apply_on), decltype(make_level)>;
    std::unique_ptr<Multigrid> mg;
    float setup_ms = 0.0f;
    cudaEvent_t t0, t1;
    MARS_CELLWISE_CK(cudaEventCreate(&t0));
    MARS_CELLWISE_CK(cudaEventCreate(&t1));
    if (use_mg) {
        MARS_CELLWISE_CK(cudaEventRecord(t0, stream));
        mg.reset(new Multigrid(apply_on, make_level, {blk, d_G, d_diag}, dec, ops.zeta.data(), pre, post,
                               stream));
        MARS_CELLWISE_CK(cudaEventRecord(t1, stream));
        MARS_CELLWISE_CK(cudaEventSynchronize(t1));
        MARS_CELLWISE_CK(cudaEventElapsedTime(&setup_ms, t0, t1));
    }
    MARS_CELLWISE_CK(cudaStreamSynchronize(stream));
    MARS_CELLWISE_MPI(MPI_Barrier(MPI_COMM_WORLD));

    MARS_CELLWISE_CK(cudaEventRecord(t0, stream));
    const cellwise::SolveResult res =
        unstructured ? cellwise::bicgstab(apply, cellwise::UnstructuredJacobi{tview, d_diag, uhalo.get()}, counting,
                                          flat, d_b, d_x, MPI_COMM_WORLD, tol, max_iterations, stream)
        : use_mg     ? cellwise::bicgstab(apply, *mg, blk, d_b, d_x, MPI_COMM_WORLD, tol, max_iterations, stream)
                     : cellwise::bicgstab(apply, cellwise::JacobiPreconditioner{d_diag, blk, halo.get()}, blk,
                                          d_b, d_x, MPI_COMM_WORLD, tol, max_iterations, stream);
    MARS_CELLWISE_CK(cudaEventRecord(t1, stream));
    MARS_CELLWISE_CK(cudaEventSynchronize(t1));
    float solve_ms = 0.0f;
    MARS_CELLWISE_CK(cudaEventElapsedTime(&solve_ms, t0, t1));

    // Relative error against the manufactured solution, in the weighted (= assembled) norm.
    double rel_err = 0.0;
    {
        cellwise::Reduction red(MPI_COMM_WORLD);
        difference_kernel<<<blocks_for(n, threads), threads, 0, stream>>>(d_b, d_x, d_uex, n);
        if (unstructured) cellwise::wdot(d_b, d_b, flat, counting, red, cellwise::kRR, stream);
        else cellwise::wdot(d_b, d_b, blk, red, cellwise::kRR, stream);
        const double err = cellwise::read_norm(red, stream);
        if (unstructured) cellwise::wdot(d_uex, d_uex, flat, counting, red, cellwise::kRR, stream);
        else cellwise::wdot(d_uex, d_uex, blk, red, cellwise::kRR, stream);
        rel_err = err / cellwise::read_norm(red, stream);
    }

    if (rank == 0) {
        printf("cell-wise BiCGStab, p=%d, %d^3 elements (%lld), %lld unique DoFs, deform %.3f, %s\n",
               kP, ne, (long long)ne * ne * ne, unique, deform,
               unstructured ? (rotate ? "unstructured tables, rotated frames, DSS + Jacobi" : "unstructured tables, DSS + Jacobi")
               : use_mg     ? "multigrid V-cycle" : "DSS + Jacobi");
        printf("  %d rank(s) as %d x %d x %d, %d x %d x %d elements each\n", dec.size, dec.P[0], dec.P[1],
               dec.P[2], blk.nx, blk.ny, blk.nz);
        if (use_mg) {
            printf("  multigrid (%d,%d) setup %.1f ms; lambda_max(P_J A) per level:", pre, post, setup_ms);
            for (double l : mg->lambda_max()) printf(" %.4f", l);
            printf("\n");
        }
        printf("  iterations %d, ||P r0||_w = %.3e, ||P r||_w = %.3e\n", res.iterations, res.residual0,
               res.residual);
        for (int k : {1, 5, 10, 20, 30, 40})
            if (k <= (int)res.history.size())
                printf("  it %3d  ||P r||_w = %.3e\n", k, res.history[k - 1]);
        printf("  ||x - u||_w / ||u||_w = %.3e  %s\n", rel_err, rel_err < 1e-8 ? "PASS" : "FAIL");
    }

    // Timing of the building blocks and of a whole iteration (rank 0's view; the
    // exchanges and global sums keep the ranks in step).
    auto time_ms = [&](auto&& f) {
        f();
        MARS_CELLWISE_CK(cudaEventRecord(t0, stream));
        for (int r = 0; r < reps; ++r) f();
        MARS_CELLWISE_CK(cudaEventRecord(t1, stream));
        MARS_CELLWISE_CK(cudaEventSynchronize(t1));
        float ms = 0.0f;
        MARS_CELLWISE_CK(cudaEventElapsedTime(&ms, t0, t1));
        return ms / reps;
    };
    // Bandwidth counts each element-local value read and written once (8 B each), plus
    // the diagonal read in the preconditioner: what a single pass must move at least.
    const double vec_gb = n * 8.0 / 1e9;
    const float op_ms = time_ms([&] { apply(d_uex, d_b); });
    const float dss_ms = unstructured ? time_ms([&] { cellwise::dss(d_b, d_x, topo, uhalo.get(), stream); })
                                      : time_ms([&] { cellwise::dss(d_b, d_x, blk, halo.get(), stream); });
    const float pre_ms = unstructured
        ? time_ms([&] {
              cellwise::UnstructuredJacobi{tview, d_diag, uhalo.get()}.operator()<false, false>(d_b, d_x, nullptr,
                                                                                               nullptr, stream);
          })
        : time_ms([&] { cellwise::precondition(d_b, d_x, d_diag, blk, halo.get(), stream); });
    const float vc_ms = use_mg ? time_ms([&] { mg->vcycle(0, d_b, d_x, stream); }) : 0.0f;
    const double it_ms = res.iterations ? solve_ms / res.iterations : 0.0;
    if (rank == 0) {
        printf("  operator     %8.3f ms  %7.1f ns/elem  %6.2f GDoF/s (unique, per rank)\n", op_ms,
               op_ms * 1e6 / E, unique / dec.size / (op_ms * 1e-3) / 1e9);
        printf("  DSS gather   %8.3f ms  %5.2f TB/s%s\n", dss_ms, 2 * vec_gb / dss_ms,
               dec.size > 1 ? "  (with the ghost exchange)" : "");
        printf("  precond      %8.3f ms  %5.2f TB/s  (gather + Jacobi, one pass)\n", pre_ms, 3 * vec_gb / pre_ms);
        if (use_mg) printf("  V-cycle      %8.3f ms  (%.1f operator calls' worth)\n", vc_ms, vc_ms / op_ms);
        printf("  iteration    %8.3f ms  (2 operator, 2 precond + dots, 3 updates)\n", it_ms);
        printf("  solve        %8.3f ms\n", solve_ms);
    }

    mg.reset();
    halo.reset();
    uhalo.reset();
    for (double* p : {d_btil, d_dtil, d_w, d_d, d_G, d_diag, d_uex, d_b, d_x})
        cudaFree(p);
    MARS_CELLWISE_CK(cudaStreamDestroy(stream));
    MPI_Finalize();
    return rel_err < 1e-8 ? 0 : 1;
}
