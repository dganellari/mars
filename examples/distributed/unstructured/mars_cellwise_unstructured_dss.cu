// The unstructured element-local DSS (backend/distributed/unstructured/solvers/
// mars_cellwise_topology.hpp) on a cube of ne^3 hexahedra, plain or with every element
// in a random proper rotation of its local frame (so all 8 face orientations and
// reversed edges occur).
//
// One rank checks:
//   1. the DSS equals a scatter-add over global node ids (to rounding);
//   2. every copy of a node holds the same bits;
//   3. the rotated cube gives the same bits per global node as the plain one: the input
//      depends only on (element, global node), and the summation order only on the
//      elements' corner keys, so the local frames must not matter;
//   4. on the plain cube it agrees with the structured DSS of mars_cellwise_layout.hpp
//      (to rounding: the structured version sums in the cascade order).
// Then times the unstructured DSS against the structured one.
//
// Several ranks (one GPU each) split the rotated cube into sub-blocks; each rank holds
// its sub-block plus the ghost elements around it and exchanges the shared nodes. Every
// local copy must hold the same bits as a one-domain DSS of the whole cube, which each
// rank also computes. Then times the DSS with the exchange.
//
// Run: mars_cellwise_unstructured_dss [--ne 8] [--reps 20]
// (several ranks need GPU-aware MPI: on Cray MPICH, MPICH_GPU_SUPPORT_ENABLED=1)

#include "backend/distributed/unstructured/solvers/mars_cellwise_krylov.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_topology.hpp"

#include <cuda_runtime.h>
#include <mpi.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

using namespace mars;
using namespace mars::cellwise;

namespace {

__constant__ int c_rot[24][9];   // the 24 proper rotations of the cube, row-major

__device__ unsigned long long mix(unsigned long long z)
{
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

__device__ int rotation_of(long long e, bool rotate) { return rotate ? (int)(mix(e) % 24) : 0; }

// A point at local reference position x (node coordinates 0..p) of an element with
// rotation r sits at the block's reference position R x (about the element centre).
__device__ void rotate_node(int r, const int (&x)[3], int (&y)[3])
{
    for (int i = 0; i < 3; ++i) {
        int s = 0;
        for (int j = 0; j < 3; ++j) s += c_rot[r][i * 3 + j] * (2 * x[j] - kP);
        y[i] = (s + kP) / 2;
    }
}

// Element i of a mesh is the cube's element elem[i] (i itself when elem is null), in
// lattice order of the ne^3 cube.
__device__ long long cube_element(const unsigned long long* elem, long long i) { return elem ? (long long)elem[i] : i; }

// Corner keys and local ids: vertex ((ex,ey,ez) + bits) of the (ne+1)^3 lattice.
__global__ void mesh_kernel(int ne, bool rotate, const unsigned long long* elem, long long E,
                            unsigned long long* const* key, int* const* lid)
{
    const long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= E) return;
    const long long e = cube_element(elem, i);
    const int ez = (int)(e % ne), ey = (int)((e / ne) % ne), ex = (int)(e / ((long long)ne * ne));
    const int r = rotation_of(e, rotate);
    for (int c = 0; c < 8; ++c) {
        int b[3];
        corner_bits(c, b[0], b[1], b[2]);
        const int x[3] = {b[0] * kP, b[1] * kP, b[2] * kP};
        int y[3];
        rotate_node(r, x, y);
        const long long v = (((long long)(ex + y[0] / kP) * (ne + 1)) + ey + y[1] / kP) * (ne + 1) + ez + y[2] / kP;
        key[c][i] = (unsigned long long)v;
        lid[c][i] = (int)v;
    }
}

// Global node id of every copy of the first E elements, and the local index the same
// node has in the plain (unrotated) element.
__global__ void node_ids_kernel(int ne, bool rotate, const unsigned long long* elem, long long E, long long* gid,
                                int* plain_local)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * kN3) return;
    const long long e = cube_element(elem, t / kN3);
    const int l = (int)(t % kN3);
    const int ez = (int)(e % ne), ey = (int)((e / ne) % ne), ex = (int)(e / ((long long)ne * ne));
    const int x[3] = {l / kNN, (l / kN) % kN, l % kN};
    int y[3];
    rotate_node(rotation_of(e, rotate), x, y);
    const long long NG = (long long)ne * kP + 1;
    gid[t] = (((long long)ex * kP + y[0]) * NG + ey * kP + y[1]) * NG + ez * kP + y[2];
    plain_local[t] = node_at(y[0], y[1], y[2]);
}

// The input depends only on (cube element, global node).
__global__ void input_kernel(const unsigned long long* elem, const long long* gid, long long n, double* v)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    const unsigned long long e = (unsigned long long)cube_element(elem, t / kN3);
    const unsigned long long h = mix(e * 0x100000001b3ULL ^ (unsigned long long)gid[t]);
    v[t] = (double)(h >> 11) / 9007199254740992.0 - 0.5;
}

// A rank's local copy against the one-domain copy of the same cube element and node.
__global__ void ranks_vs_one_kernel(const double* out, const unsigned long long* elem, long long n,
                                    const double* one, int* bitdiff)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    if (__double_as_longlong(out[t]) != __double_as_longlong(one[(long long)elem[t / kN3] * kN3 + t % kN3]))
        atomicAdd(bitdiff, 1);
}

__global__ void scatter_add_kernel(const double* v, const long long* gid, long long n, double* acc)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) atomicAdd(&acc[gid[t]], v[t]);
}

// Racy on purpose: if every copy of a node holds the same bits, any winner is right.
__global__ void scatter_any_kernel(const double* v, const long long* gid, long long n, double* g)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) g[gid[t]] = v[t];
}

// |out - acc[gid]| / scale, and the number of copies that differ from g[gid] in any bit.
__global__ void compare_kernel(const double* out, const long long* gid, const double* acc, const double* g,
                               long long n, double* err, int* bitdiff)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    err[t] = fabs(out[t] - acc[gid[t]]);
    if (out[t] != g[gid[t]]) atomicAdd(bitdiff, 1);
}

// Rotated vs plain: same element, same global node, the copies must match bit for bit.
__global__ void rotated_vs_plain_kernel(const double* rot, const double* plain, const int* plain_local,
                                        long long n, int* bitdiff)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    if (rot[t] != plain[(t / kN3) * kN3 + plain_local[t]]) atomicAdd(bitdiff, 1);
}

__global__ void abs_diff_kernel(const double* a, const double* b, long long n, double* d)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < n) d[t] = fabs(a[t] - b[t]);
}

unsigned grid_of(long long n) { return (unsigned)((n + kThreads - 1) / kThreads); }

template <typename T>
T* raw(thrust::device_vector<T>& v) { return thrust::raw_pointer_cast(v.data()); }

struct Mesh {
    thrust::device_vector<unsigned long long> key[8];
    thrust::device_vector<int> lid[8];
    thrust::device_vector<long long> gid;
    thrust::device_vector<int> plain_local;
    UnstructuredTopology topo;
};

// Tables over E elements (cube elements elem[0 .. E), or the whole cube when elem is
// null), the first `local` of them this rank's; node ids for the local ones.
void build(Mesh& m, int ne, bool rotate, const unsigned long long* elem, long long E, long long local)
{
    std::vector<unsigned long long*> kp(8);
    std::vector<int*> lp(8);
    for (int c = 0; c < 8; ++c) {
        m.key[c].resize(E);
        m.lid[c].resize(E);
        kp[c] = raw(m.key[c]);
        lp[c] = raw(m.lid[c]);
    }
    thrust::device_vector<unsigned long long*> d_kp(kp.begin(), kp.end());
    thrust::device_vector<int*> d_lp(lp.begin(), lp.end());
    mesh_kernel<<<grid_of(E), kThreads>>>(ne, rotate, elem, E, raw(d_kp), raw(d_lp));
    m.gid.resize(local * kN3);
    m.plain_local.resize(local * kN3);
    node_ids_kernel<<<grid_of(local * kN3), kThreads>>>(ne, rotate, elem, local, raw(m.gid), raw(m.plain_local));
    MARS_CELLWISE_CK(cudaGetLastError());
    std::vector<const unsigned long long*> ckp(kp.begin(), kp.end());
    std::vector<const int*> clp(lp.begin(), lp.end());
    thrust::device_vector<const unsigned long long*> d_ckp(ckp.begin(), ckp.end());
    thrust::device_vector<const int*> d_clp(clp.begin(), clp.end());
    m.topo = build_topology(raw(d_ckp), raw(d_clp), E, local);
}

template <typename F>
float time_ms(int reps, F&& f)
{
    cudaEvent_t t0, t1;
    MARS_CELLWISE_CK(cudaEventCreate(&t0));
    MARS_CELLWISE_CK(cudaEventCreate(&t1));
    f();
    MARS_CELLWISE_CK(cudaEventRecord(t0));
    for (int r = 0; r < reps; ++r) f();
    MARS_CELLWISE_CK(cudaEventRecord(t1));
    MARS_CELLWISE_CK(cudaEventSynchronize(t1));
    float ms = 0.0f;
    MARS_CELLWISE_CK(cudaEventElapsedTime(&ms, t0, t1));
    MARS_CELLWISE_CK(cudaEventDestroy(t0));
    MARS_CELLWISE_CK(cudaEventDestroy(t1));
    return ms / reps;
}

bool one_rank(int ne, int reps)
{
    const long long E = (long long)ne * ne * ne, n = E * kN3, NG = (long long)ne * kP + 1;
    Mesh plain, rotated;
    build(plain, ne, false, nullptr, E, E);
    build(rotated, ne, true, nullptr, E, E);
    const UnstructuredTopology& tr = rotated.topo;
    printf("unstructured cell-wise DSS, %d^3 elements, %lld unique nodes\n", ne, NG * NG * NG);
    printf("  tables: %lld face pairs, %lld single faces, %lld edges, %lld vertices\n", tr.face_pairs,
           tr.single_faces, tr.edges(), tr.vertices());

    bool ok = true;
    thrust::device_vector<double> in(n), out_plain(n), out_rot(n), acc(NG * NG * NG), g(NG * NG * NG), err(n);
    thrust::device_vector<int> bitdiff(1);
    for (int pass = 0; pass < 2; ++pass) {
        Mesh& m = pass ? rotated : plain;
        thrust::device_vector<double>& out = pass ? out_rot : out_plain;
        input_kernel<<<grid_of(n), kThreads>>>(nullptr, raw(m.gid), n, raw(in));
        dss(raw(in), raw(out), m.topo);
        thrust::fill(acc.begin(), acc.end(), 0.0);
        scatter_add_kernel<<<grid_of(n), kThreads>>>(raw(in), raw(m.gid), n, raw(acc));
        scatter_any_kernel<<<grid_of(n), kThreads>>>(raw(out), raw(m.gid), n, raw(g));
        bitdiff[0] = 0;
        compare_kernel<<<grid_of(n), kThreads>>>(raw(out), raw(m.gid), raw(acc), raw(g), n, raw(err), raw(bitdiff));
        MARS_CELLWISE_CK(cudaDeviceSynchronize());
        const double e = *thrust::max_element(err.begin(), err.end());
        const int bd = bitdiff[0];
        printf("  %s: vs scatter-add %.1e, copies differing from each other %d\n", pass ? "rotated" : "plain  ",
               e, bd);
        ok &= e < 1e-13 && bd == 0;
    }
    bitdiff[0] = 0;
    rotated_vs_plain_kernel<<<grid_of(n), kThreads>>>(raw(out_rot), raw(out_plain), raw(rotated.plain_local), n,
                                                      raw(bitdiff));
    MARS_CELLWISE_CK(cudaDeviceSynchronize());
    printf("  rotated vs plain, per global node: %d copies differ in any bit\n", (int)bitdiff[0]);
    ok &= bitdiff[0] == 0;

    // The structured DSS on the plain cube (same element order and frames).
    thrust::device_vector<double> out_struct(n), diff(n);
    input_kernel<<<grid_of(n), kThreads>>>(nullptr, raw(plain.gid), n, raw(in));
    cellwise::dss(raw(in), raw(out_struct), Block(ne, ne, ne), nullptr);
    abs_diff_kernel<<<grid_of(n), kThreads>>>(raw(out_struct), raw(out_plain), n, raw(diff));
    MARS_CELLWISE_CK(cudaDeviceSynchronize());
    const double sd = *thrust::max_element(diff.begin(), diff.end());
    printf("  plain vs structured DSS: max |difference| %.1e (summation orders differ)\n", sd);
    ok &= sd < 1e-13;

    const double unique = (double)NG * NG * NG, gb = n * 16.0 / 1e9;   // one read + one write per value
    const float ms_struct = time_ms(reps, [&] { cellwise::dss(raw(in), raw(out_struct), Block(ne, ne, ne), nullptr); });
    const float ms_plain = time_ms(reps, [&] { dss(raw(in), raw(out_plain), plain.topo); });
    const float ms_rot = time_ms(reps, [&] { dss(raw(in), raw(out_rot), rotated.topo); });
    auto row = [&](const char* name, float ms) {
        printf("  %-28s %8.3f ms  %6.2f GDoF/s  %5.2f TB/s\n", name, ms, unique / ms / 1e6, gb / ms);
    };
    row("structured gather", ms_struct);
    row("unstructured, plain", ms_plain);
    row("unstructured, rotated", ms_rot);
    return ok;
}

bool several_ranks(int ne, int reps, MPI_Comm comm)
{
    const Decomposition dec(comm, ne, ne, ne);
    thrust::device_vector<unsigned long long> elem;
    thrust::device_vector<int> owner;
    const long long L = block_elements(dec, elem, owner), E = (long long)elem.size(), n = L * kN3;
    Mesh m;
    build(m, ne, true, raw(elem), E, L);
    UnstructuredHalo halo(m.topo, raw(elem), raw(owner), comm);
    thrust::device_vector<double> in(n), out(n);
    input_kernel<<<grid_of(n), kThreads>>>(raw(elem), raw(m.gid), n, raw(in));
    dss(raw(in), raw(out), m.topo, &halo);

    // The same cube as one domain on this rank.
    const long long EC = (long long)ne * ne * ne;
    Mesh one;
    build(one, ne, true, nullptr, EC, EC);
    thrust::device_vector<double> in_one(EC * kN3), out_one(EC * kN3);
    input_kernel<<<grid_of(EC * kN3), kThreads>>>(nullptr, raw(one.gid), EC * kN3, raw(in_one));
    dss(raw(in_one), raw(out_one), one.topo);
    thrust::device_vector<int> bitdiff(1, 0);
    ranks_vs_one_kernel<<<grid_of(n), kThreads>>>(raw(out), raw(elem), n, raw(out_one), raw(bitdiff));
    MARS_CELLWISE_CK(cudaDeviceSynchronize());
    int bd = bitdiff[0];
    MARS_CELLWISE_MPI(MPI_Allreduce(MPI_IN_PLACE, &bd, 1, MPI_INT, MPI_SUM, comm));

    MARS_CELLWISE_MPI(MPI_Barrier(comm));
    const float ms = time_ms(reps, [&] { dss(raw(in), raw(out), m.topo, &halo); });
    if (dec.rank == 0) {
        printf("unstructured cell-wise DSS on %d ranks (%d x %d x %d), %d^3 elements, rotated frames\n", dec.size,
               dec.P[0], dec.P[1], dec.P[2], ne);
        printf("  rank 0: %lld local + %lld ghost elements (%lld read ghosts), %lld values sent to %d ranks\n", L,
               E - L, halo.boundary().count, halo.values_sent(), halo.peers());
        printf("  several ranks vs one domain: %d copies differ in any bit\n", bd);
        printf("  DSS with the exchange  %8.3f ms (rank 0)\n", ms);
    }
    return bd == 0;
}

}  // namespace

int main(int argc, char** argv)
{
    MARS_CELLWISE_MPI(MPI_Init(&argc, &argv));
    int rank = 0, size = 1;
    MARS_CELLWISE_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
    MARS_CELLWISE_MPI(MPI_Comm_size(MPI_COMM_WORLD, &size));
    // One GPU per rank, by the rank's index on its node.
    {
        MPI_Comm node;
        MARS_CELLWISE_MPI(MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node));
        int local = 0, devices = 0;
        MARS_CELLWISE_MPI(MPI_Comm_rank(node, &local));
        MARS_CELLWISE_MPI(MPI_Comm_free(&node));
        MARS_CELLWISE_CK(cudaGetDeviceCount(&devices));
        MARS_CELLWISE_CK(cudaSetDevice(local % devices));
    }
    int ne = 8, reps = 20;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "--ne")) ne = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--reps")) reps = atoi(argv[i + 1]);
        else { if (rank == 0) fprintf(stderr, "unknown option %s\n", argv[i]); MPI_Finalize(); return 1; }
    }
    if (size > 1) {   // device pointers go to MPI: MPICH needs GPU support switched on
        char version[MPI_MAX_LIBRARY_VERSION_STRING];
        int len = 0;
        MARS_CELLWISE_MPI(MPI_Get_library_version(version, &len));
        const char* gpu = getenv("MPICH_GPU_SUPPORT_ENABLED");
        if (strstr(version, "MPICH") && !(gpu && !strcmp(gpu, "1"))) {
            if (rank == 0) fprintf(stderr, "MPICH_GPU_SUPPORT_ENABLED=1 is required on several ranks (GPU-aware MPI)\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    // The 24 proper rotations: signed permutation matrices of determinant +1.
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
    // Identity first, so rotation 0 means "plain".
    MARS_CELLWISE_CK(cudaMemcpyToSymbol(c_rot, rot.data(), rot.size() * sizeof(int)));
    MARS_CELLWISE_CK(cudaDeviceSynchronize());

    const bool ok = size == 1 ? one_rank(ne, reps) : several_ranks(ne, reps, MPI_COMM_WORLD);
    if (rank == 0) printf("UNSTRUCTURED DSS: %s\n", ok ? "PASS" : "FAIL");
    MPI_Finalize();
    return ok ? 0 : 1;
}
