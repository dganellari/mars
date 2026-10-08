// Element-local (cell-wise) solve of the high-order CVFEM Laplacian on one GPU: the
// mixed form of the plan from internal-notes/FlexibleCG.pdf. The operator is the
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
// All per-element work runs on the GPU: corners, metric, diagonal, u, b.
//
// Run:  mars_marsir_cellwise_solve --ptx <hl_full_p7_sm90.ptx> [--ne 8] [--deform 0]
//                                  [--tol 1e-10] [--maxit 500] [--reps 10] [--mg]
// --mg preconditions with the geometric multigrid V-cycle instead of DSS + Jacobi
// (ne must be a power of two); its history must match
// marsir-mlir/test/cellwise_multigrid_ref.py.
// With the same --ne and --deform, the residual history must match
// marsir-mlir/test/cellwise_krylov_ref.py to about 4 digits (BiCGStab amplifies the
// rounding of the different summation orders).

#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"
#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_operator.hpp"
#include "backend/distributed/unstructured/solvers/mars_cellwise_multigrid.hpp"

#include <cuda_runtime.h>

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

// Vertex (i, j, k) of the global vertex lattice, deformed in place.
__device__ void vertex(int i, int j, int k, int ne, double deform, double x[3])
{
    const double h = 1.0 / ne;
    x[0] = i * h;
    x[1] = j * h;
    x[2] = k * h;
    const double bump = deform * sin(kPi * x[0]) * sin(kPi * x[1]) * sin(kPi * x[2]);
    x[0] += bump;
    x[1] += bump;
    x[2] += bump;
}

// The 8 corners of element e in the order of c_hexCornerRef (sign -1 = low side).
__device__ void element_corners(long long e, int ne, double deform, double corners[8][3])
{
    const int ez = (int)(e % ne);
    const int ey = (int)((e / ne) % ne);
    const int ex = (int)(e / ((long long)ne * ne));
    for (int c = 0; c < 8; ++c)
        vertex(ex + (fem::c_hexCornerRef[c][0] > 0), ey + (fem::c_hexCornerRef[c][1] > 0),
               ez + (fem::c_hexCornerRef[c][2] > 0), ne, deform, corners[c]);
}

// One thread per metric point (e, dir, l, s, r): MARS's device metric, written in
// the component-major layout the MARSIR kernel reads, [dir][face][comp][s][r].
__global__ void metric_kernel(double* __restrict__ G, long long E, int ne, double deform)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    const long long points = 3LL * kP * kNN;
    if (t >= E * points) return;
    const long long e = t / points;
    int q = (int)(t % points);
    const int r = q % kN;
    q /= kN;
    const int s = q % kN;
    q /= kN;
    const int l = q % kP;
    const int dir = q / kP;
    double corners[8][3];
    element_corners(e, ne, deform, corners);
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
__global__ void exact_kernel(double* __restrict__ u, long long E, int ne, double deform)
{
    const long long t = (long long)blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= E * kN3) return;
    const long long e = t / kN3;
    const int node = (int)(t % kN3);
    const double rf[3] = {fem::c_zeta[node / kNN], fem::c_zeta[(node / kN) % kN],
                          fem::c_zeta[node % kN]};
    double corners[8][3];
    element_corners(e, ne, deform, corners);
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
    std::string ptx_path;
    int ne = 8, max_iterations = 500, reps = 10;
    double deform = 0.0, tol = 1e-10;
    bool use_mg = false;
    for (int i = 1; i < argc; i += 2) {
        if (!strcmp(argv[i], "--mg")) { use_mg = true; --i; continue; }
        if (i + 1 >= argc) { fprintf(stderr, "option %s needs a value\n", argv[i]); return 1; }
        if (!strcmp(argv[i], "--ptx")) ptx_path = argv[i + 1];
        else if (!strcmp(argv[i], "--ne")) ne = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--deform")) deform = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "--tol")) tol = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "--maxit")) max_iterations = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "--reps")) reps = atoi(argv[i + 1]);
        else { fprintf(stderr, "unknown option %s\n", argv[i]); return 1; }
    }
    if (ptx_path.empty()) {
        fprintf(stderr, "usage: %s --ptx <hl_full_p7_sm90.ptx> [--ne N] [--deform a] "
                        "[--tol t] [--maxit k] [--reps r] [--mg]\n", argv[0]);
        return 1;
    }
    if (use_mg && (ne & (ne - 1))) {
        fprintf(stderr, "--mg needs ne to be a power of two\n");
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

    const cellwise::Block blk{ne, ne, ne};
    const long long E = blk.elements(), n = blk.values();
    const long long unique = (long long)(ne * kP + 1) * (ne * kP + 1) * (ne * kP + 1);
    constexpr int threads = 256;

    double* d_G = device_array(E * kGElem);
    double* d_diag = device_array(n);
    double* d_uex = device_array(n);
    double* d_b = device_array(n);
    double* d_x = device_array(n);

    metric_kernel<<<blocks_for(E * 3LL * kP * kNN, threads), threads>>>(d_G, E, ne, deform);
    diagonal_kernel<<<blocks_for(n, threads), threads>>>(d_G, d_b, E);   // d_b: scratch until b = A u
    exact_kernel<<<blocks_for(n, threads), threads>>>(d_uex, E, ne, deform);
    MARS_CELLWISE_CK(cudaGetLastError());
    cellwise::dss(d_b, d_diag, blk);     // assembled diagonal on every copy

    const marsir::LaplacianPtx op(ptx_path, kP);
    auto apply = [&](const double* u, double* y) {
        op.apply(u, d_btil, d_dtil, d_w, d_d, d_G, y, E);
    };
    apply(d_uex, d_b);                   // b = A u, unassembled

    // Multigrid levels: the same mesh with every block dimension halved, down to one
    // element; each level gets its own metric and assembled diagonal, on the GPU.
    auto apply_on = [&](const double* u, double* y, const double* G, long long elements) {
        op.apply(u, d_btil, d_dtil, d_w, d_d, G, y, elements);
    };
    std::vector<double*> level_arrays;
    std::unique_ptr<cellwise::Multigrid<decltype(apply_on)>> mg;
    float setup_ms = 0.0f;
    cudaEvent_t t0, t1;
    MARS_CELLWISE_CK(cudaEventCreate(&t0));
    MARS_CELLWISE_CK(cudaEventCreate(&t1));
    if (use_mg) {
        MARS_CELLWISE_CK(cudaEventRecord(t0));
        std::vector<cellwise::MultigridLevel> levels{{blk, d_G, d_diag}};
        for (int m = ne / 2; m >= 1; m /= 2) {
            const cellwise::Block bl{m, m, m};
            const long long em = bl.elements(), nm = bl.values();
            double* g = device_array(em * kGElem);
            double* dg = device_array(nm);
            metric_kernel<<<blocks_for(em * 3LL * kP * kNN, threads), threads>>>(g, em, m, deform);
            diagonal_kernel<<<blocks_for(nm, threads), threads>>>(g, d_x, em);   // d_x: scratch
            cellwise::dss(d_x, dg, bl);
            levels.push_back({bl, g, dg});
            level_arrays.push_back(g);
            level_arrays.push_back(dg);
        }
        mg = std::make_unique<cellwise::Multigrid<decltype(apply_on)>>(apply_on, levels,
                                                                      ops.zeta.data());
        MARS_CELLWISE_CK(cudaEventRecord(t1));
        MARS_CELLWISE_CK(cudaEventSynchronize(t1));
        MARS_CELLWISE_CK(cudaEventElapsedTime(&setup_ms, t0, t1));
    }
    MARS_CELLWISE_CK(cudaDeviceSynchronize());

    MARS_CELLWISE_CK(cudaEventRecord(t0));
    const cellwise::SolveResult res =
        use_mg ? cellwise::bicgstab(apply, *mg, blk, d_b, d_x, tol, max_iterations)
               : cellwise::bicgstab(apply, cellwise::JacobiPreconditioner{d_diag, blk}, blk, d_b,
                                    d_x, tol, max_iterations);
    MARS_CELLWISE_CK(cudaEventRecord(t1));
    MARS_CELLWISE_CK(cudaEventSynchronize(t1));
    float solve_ms = 0.0f;
    MARS_CELLWISE_CK(cudaEventElapsedTime(&solve_ms, t0, t1));

    // Relative error against the manufactured solution, in the weighted (= assembled) norm.
    double rel_err = 0.0;
    {
        cellwise::Reduction red;
        difference_kernel<<<blocks_for(n, threads), threads>>>(d_b, d_x, d_uex, n);
        cellwise::wdot(d_b, d_b, blk, red, cellwise::kRR);
        const double err = cellwise::read_norm(red, 0);
        cellwise::wdot(d_uex, d_uex, blk, red, cellwise::kRR);
        rel_err = err / cellwise::read_norm(red, 0);
    }

    printf("cell-wise BiCGStab, p=%d, %d^3 elements (%lld), %lld unique DoFs, deform %.3f, %s\n",
           kP, ne, E, unique, deform, use_mg ? "multigrid V-cycle" : "DSS + Jacobi");
    if (use_mg) {
        printf("  multigrid setup %.1f ms; lambda_max(P_J A) per level:", setup_ms);
        for (int l = 0; (ne >> l) > 1; ++l) printf(" %.4f", mg->lambda_max(l));
        printf("\n");
    }
    printf("  iterations %d, ||P r0||_w = %.3e, ||P r||_w = %.3e\n", res.iterations,
           res.residual0, res.residual);
    for (int k : {1, 5, 10, 20, 30, 40})
        if (k <= (int)res.history.size())
            printf("  it %3d  ||P r||_w = %.3e\n", k, res.history[k - 1]);
    printf("  ||x - u||_w / ||u||_w = %.3e  %s\n", rel_err, rel_err < 1e-8 ? "PASS" : "FAIL");

    // Timing of the three building blocks and of a whole iteration.
    auto time_ms = [&](auto&& f) {
        f();
        MARS_CELLWISE_CK(cudaEventRecord(t0));
        for (int r = 0; r < reps; ++r) f();
        MARS_CELLWISE_CK(cudaEventRecord(t1));
        MARS_CELLWISE_CK(cudaEventSynchronize(t1));
        float ms = 0.0f;
        MARS_CELLWISE_CK(cudaEventElapsedTime(&ms, t0, t1));
        return ms / reps;
    };
    // Bandwidth counts each element-local value read and written once (8 B each), plus
    // the diagonal read in the preconditioner: what a single pass must move at least.
    const double vec_gb = n * 8.0 / 1e9;
    const float op_ms = time_ms([&] { apply(d_uex, d_b); });
    const float dss_ms = time_ms([&] { cellwise::dss(d_b, d_x, blk); });
    const float cascade_ms = time_ms([&] { cellwise::dss_cascade(d_x, blk); });
    const float pre_ms = time_ms([&] { cellwise::precondition(d_b, d_x, d_diag, blk); });
    const double it_ms = res.iterations ? solve_ms / res.iterations : 0.0;
    printf("  operator     %8.3f ms  %7.1f ns/elem  %6.2f GDoF/s (unique)\n", op_ms,
           op_ms * 1e6 / E, unique / (op_ms * 1e-3) / 1e9);
    printf("  DSS gather   %8.3f ms  %6.2f GDoF/s (unique)  %5.2f TB/s\n", dss_ms,
           unique / (dss_ms * 1e-3) / 1e9, 2 * vec_gb / dss_ms);
    printf("  DSS cascade  %8.3f ms  %6.2f GDoF/s (unique)  (in place, 3 passes)\n", cascade_ms,
           unique / (cascade_ms * 1e-3) / 1e9);
    printf("  precond      %8.3f ms  %5.2f TB/s  (gather + Jacobi, one pass)\n", pre_ms,
           3 * vec_gb / pre_ms);
    if (use_mg) {
        const float vc_ms = time_ms([&] { mg->vcycle(0, d_b, d_x, 0); });
        printf("  V-cycle      %8.3f ms  (%.1f operator calls' worth)\n", vc_ms, vc_ms / op_ms);
    }
    printf("  iteration    %8.3f ms  (2 operator, 2 precond + dots, 3 updates)\n", it_ms);

    mg.reset();
    for (double* p : level_arrays) cudaFree(p);
    for (double* p : {d_btil, d_dtil, d_w, d_d, d_G, d_diag, d_uex, d_b, d_x})
        cudaFree(p);
    return rel_err < 1e-8 ? 0 : 1;
}
