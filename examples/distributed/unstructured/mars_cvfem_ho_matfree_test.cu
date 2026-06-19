// GPU correctness + perf gate for the optimized HO-CVFEM matrix-free diffusion
// apply (mars_cvfem_ho_matfree.hpp). Single rank, in-memory structured cube --
// reuses the host patch-test setup (HODofHandler + the same elemDof/coord
// convention) so any divergence localizes to the GPU kernel, not the DOF map.
//
// Three gates per order p in {1,2,4}:
//   (A) ELEMENT bit-exactness: GPU PerPoint metric kernel must reproduce the host
//       computeElementMetric layout, and the GPU apply on one element must match
//       host applyHoCvfemElement (general cross-term path) to <1e-12 relative
//       (only FP reduction-order differs).
//   (B) PATCH A*1 == 0 everywhere (constant nullspace), assembled over the cube.
//   (C) PATCH A*linear == 0 at INTERIOR DOFs (consistency).
// Then a timing loop reports MDOF/s and an estimated effective GB/s per p.
//
// d_y MUST be zeroed before every apply (scatter is additive); operators are
// uploaded once per p via ho_cvfem_upload_operators.

#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"

#include <cuda_runtime.h>
#include <array>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace mars::fem;

// One row of the matrix-free order sweep (throughput + memory per order).
struct SweepRow { int p; long dofs; double ms, mdofs, gbs, mfBpd, asmBpd; };
static std::vector<SweepRow> g_sweep;

// variadic so template-argument commas (e.g. <double, P>) inside the call don't
// get split into separate macro arguments by the preprocessor.
#define CK(...) do { cudaError_t e_=(__VA_ARGS__); if(e_!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return false; } } while(0)

// Build the structured-cube mesh (E^3 hexes) + HODofHandler, identical to the
// host patch test. Returns corners (per element, 8x xyz) for the metric kernel.
struct CubeMesh {
    HODofHandler dh;
    std::vector<std::array<int,3>> ijk;
    long nDof;
    size_t nEl;
    int n, N3;
    std::vector<double> h_corners;   // [nEl*8*3]
};

static CubeMesh buildCube(const HoCvfemOperators& op, int P, int E)
{
    CubeMesh m;
    const int n = P + 1, N3 = n * n * n;
    auto cg = [&](int x, int y, int z) { return (x*(E+1)+y)*(E+1)+z; };
    std::vector<std::array<int,8>> ec;
    for (int ex=0; ex<E; ++ex) for (int ey=0; ey<E; ++ey) for (int ez=0; ez<E; ++ez) {
        ec.push_back({cg(ex,ey,ez),cg(ex+1,ey,ez),cg(ex+1,ey+1,ez),cg(ex,ey+1,ez),
                      cg(ex,ey,ez+1),cg(ex+1,ey,ez+1),cg(ex+1,ey+1,ez+1),cg(ex,ey+1,ez+1)});
        m.ijk.push_back({ex,ey,ez});
    }
    m.dh.build(ec, long(E+1)*(E+1)*(E+1), P);
    m.nDof = m.dh.numDof; m.nEl = ec.size(); m.n = n; m.N3 = N3;

    // Physical corner coords (unit cube, h = 1/E). Hex corner order matches the
    // host hexCornerRef / handler convention.
    const double h = 1.0 / E;
    static const int corncoord[8][3] =
        {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    m.h_corners.resize(m.nEl * 24);
    for (size_t e = 0; e < m.nEl; ++e)
        for (int c = 0; c < 8; ++c) {
            m.h_corners[e*24 + c*3 + 0] = (m.ijk[e][0] + corncoord[c][0]) * h;
            m.h_corners[e*24 + c*3 + 1] = (m.ijk[e][1] + corncoord[c][1]) * h;
            m.h_corners[e*24 + c*3 + 2] = (m.ijk[e][2] + corncoord[c][2]) * h;
        }
    return m;
}

// One full GPU assembled apply: zero d_y, launch. Templated on P so the launcher
// picks the per-P block/elems defaults.
template<int P>
static bool gpuApply(const double* d_u, double* d_y, const int* d_elemDof,
                     const double* d_G, size_t nEl, long nDof)
{
    CK(cudaMemset(d_y, 0, sizeof(double) * nDof));
    CK(ho_cvfem_apply_launch<double, P>(d_u, d_y, d_elemDof, d_G, nEl));
    CK(cudaDeviceSynchronize());
    return true;
}

template<int P>
static bool runOrder(int E)
{
    auto op = buildHoCvfemOperators(P);
    const int n = P + 1, N3 = n * n * n;
    CubeMesh m = buildCube(op, P, E);
    const size_t nEl = m.nEl; const long nDof = m.nDof;

    // Upload reference operators (constant memory) and xi/zeta for the metric.
    CK(ho_cvfem_upload_operators(P, op.Btil.data(), op.Dtil.data(),
                                 op.D.data(), op.W.data(), op.xi.data(), op.zeta.data()));

    // Device buffers.
    int*    d_elemDof; double* d_corners; double* d_G;
    double* d_u; double* d_y;
    const size_t gLen = nEl * (size_t)(3 * P * n * n) * 3;
    CK(cudaMalloc(&d_elemDof, sizeof(int)    * nEl * N3));
    CK(cudaMalloc(&d_corners, sizeof(double) * nEl * 24));
    CK(cudaMalloc(&d_G,       sizeof(double) * gLen));
    CK(cudaMalloc(&d_u,       sizeof(double) * nDof));
    CK(cudaMalloc(&d_y,       sizeof(double) * nDof));
    CK(cudaMemcpy(d_elemDof, m.dh.elemDof.data(), sizeof(int) * nEl * N3, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_corners, m.h_corners.data(),  sizeof(double) * nEl * 24, cudaMemcpyHostToDevice));

    // Build the PerPoint metric on device.
    CK(ho_cvfem_metric_perpoint_launch<double, P>(d_corners, d_G, nEl));
    CK(cudaDeviceSynchronize());

    // ---- Gate (A.1): GPU metric == host computeElementMetric, element 0. ----
    double corners0[8][3];
    for (int c=0;c<8;++c) for (int d=0;d<3;++d) corners0[c][d] = m.h_corners[c*3+d];
    auto Ghost0 = computeElementMetric(op, corners0);
    std::vector<double> Ggpu0(Ghost0.size()*3);
    CK(cudaMemcpy(Ggpu0.data(), d_G, sizeof(double)*Ghost0.size()*3, cudaMemcpyDeviceToHost));
    double metricErr = 0;
    for (size_t i=0;i<Ghost0.size();++i) for (int c=0;c<3;++c)
        metricErr = std::max(metricErr, std::abs(Ggpu0[i*3+c] - Ghost0[i][c]));

    // ---- Gate (A.2): single-element apply bit-exactness vs host general path. ----
    // Random u over the cube; compare element 0's contribution. We run the host
    // applyHoCvfemElement on element 0 with its own gathered u and host metric,
    // and compare to the GPU per-element output. To isolate one element we apply
    // the GPU kernel with a 1-element elemDof identity map.
    std::mt19937 rng(12345); std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<double> uGlob(nDof);
    for (long i=0;i<nDof;++i) uGlob[i] = dist(rng);

    // host element-0 reference
    std::vector<double> ul(N3), yl_host;
    for (int l=0;l<N3;++l) ul[l] = uGlob[m.dh.elemDof[0*N3 + l]];
    applyHoCvfemElement(op, Ghost0, ul, yl_host);

    // GPU single element: identity elemDof [0..N3), metric = element 0's G.
    int* d_edof1; double* d_u1; double* d_y1;
    CK(cudaMalloc(&d_edof1, sizeof(int)*N3));
    CK(cudaMalloc(&d_u1, sizeof(double)*N3));
    CK(cudaMalloc(&d_y1, sizeof(double)*N3));
    std::vector<int> eid(N3); for (int l=0;l<N3;++l) eid[l]=l;
    CK(cudaMemcpy(d_edof1, eid.data(), sizeof(int)*N3, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_u1, ul.data(), sizeof(double)*N3, cudaMemcpyHostToDevice));
    CK(cudaMemset(d_y1, 0, sizeof(double)*N3));
    CK(ho_cvfem_apply_launch<double, P>(d_u1, d_y1, d_edof1, d_G, 1));
    CK(cudaDeviceSynchronize());
    std::vector<double> yl_gpu(N3);
    CK(cudaMemcpy(yl_gpu.data(), d_y1, sizeof(double)*N3, cudaMemcpyDeviceToHost));
    double ynorm=0, elemErr=0;
    for (int l=0;l<N3;++l) { ynorm = std::max(ynorm, std::abs(yl_host[l]));
        elemErr = std::max(elemErr, std::abs(yl_gpu[l]-yl_host[l])); }
    double elemRel = elemErr / (ynorm > 0 ? ynorm : 1.0);

    // ---- Gate (B): A*1 == 0 everywhere. ----
    std::vector<double> ones(nDof, 1.0);
    CK(cudaMemcpy(d_u, ones.data(), sizeof(double)*nDof, cudaMemcpyHostToDevice));
    if (!gpuApply<P>(d_u, d_y, d_elemDof, d_G, nEl, nDof)) return false;
    std::vector<double> y1(nDof);
    CK(cudaMemcpy(y1.data(), d_y, sizeof(double)*nDof, cudaMemcpyDeviceToHost));
    double nullMax=0; for (long i=0;i<nDof;++i) nullMax = std::max(nullMax, std::abs(y1[i]));

    // ---- Gate (C): A*linear == 0 at interior DOFs. ----
    // Build dofX + boundary flag exactly as the host patch test.
    std::vector<double> dofX(nDof,0.0); std::vector<uint8_t> bdry(nDof,0);
    const double h = 1.0/E;
    auto coord = [&](int e_ijk, int loc){ return (e_ijk + 0.5*(op.zeta[loc]+1.0))*h; };
    for (size_t e=0;e<nEl;++e)
        for (int i=0;i<n;++i) for (int j=0;j<n;++j) for (int k=0;k<n;++k) {
            int dof = m.dh.elemDof[e*N3 + i*n*n+j*n+k];
            double px=coord(m.ijk[e][0],i), py=coord(m.ijk[e][1],j), pz=coord(m.ijk[e][2],k);
            dofX[dof]=px;
            bdry[dof]=(px<1e-12||px>1-1e-12||py<1e-12||py>1-1e-12||pz<1e-12||pz>1-1e-12)?1:0;
        }
    CK(cudaMemcpy(d_u, dofX.data(), sizeof(double)*nDof, cudaMemcpyHostToDevice));
    if (!gpuApply<P>(d_u, d_y, d_elemDof, d_G, nEl, nDof)) return false;
    std::vector<double> yx(nDof);
    CK(cudaMemcpy(yx.data(), d_y, sizeof(double)*nDof, cudaMemcpyDeviceToHost));
    double interMax=0; for (long i=0;i<nDof;++i) if(!bdry[i]) interMax=std::max(interMax,std::abs(yx[i]));

    // ---- Timing loop (random u). MDOF/s + effective GB/s estimate. ----
    CK(cudaMemcpy(d_u, uGlob.data(), sizeof(double)*nDof, cudaMemcpyHostToDevice));
    const int warm=5, iters=100;
    for (int it=0; it<warm; ++it) {
        CK(cudaMemset(d_y, 0, sizeof(double)*nDof));
        CK(ho_cvfem_apply_launch<double, P>(d_u, d_y, d_elemDof, d_G, nEl));
    }
    CK(cudaDeviceSynchronize());
    cudaEvent_t t0,t1; cudaEventCreate(&t0); cudaEventCreate(&t1);
    cudaEventRecord(t0);
    for (int it=0; it<iters; ++it) {
        CK(cudaMemset(d_y, 0, sizeof(double)*nDof));
        CK(ho_cvfem_apply_launch<double, P>(d_u, d_y, d_elemDof, d_G, nEl));
    }
    cudaEventRecord(t1); CK(cudaEventSynchronize(t1));
    float ms=0; cudaEventElapsedTime(&ms, t0, t1);
    double sPerApply = (ms/1e3)/iters;
    double mdofs = (double)nDof / sPerApply / 1e6;
    // Traffic estimate (DRAM lower bound; EXCLUDES smem round-trips, which the
    // ncu prior says dominate). gather u (N3 doubles/elem), scatter y (N3/elem),
    // metric G read (3*p*n*n*3 doubles/elem PerPoint), and the elemDof read
    // (N3 ints = N3*4 bytes/elem -- real DRAM, ~25% of bytes at p=1, omitting it
    // understates the floor).
    double bytesPerElem = (double)(2*N3 + 3*P*n*n*3) * 8.0 + (double)N3 * 4.0;
    double gbs = bytesPerElem * nEl / sPerApply / 1e9;
    cudaEventDestroy(t0); cudaEventDestroy(t1);

    bool okA = (metricErr < 1e-12) && (elemRel < 1e-12);
    bool okB = (nullMax < 1e-9);
    bool okC = (interMax < 1e-9);
    bool ok = okA && okB && okC;
    printf("p=%d E=%d nDof=%ld nEl=%zu | metricErr=%.2e elemRel=%.2e | A*1=%.2e A*lin(int)=%.2e | %s\n",
           P, E, nDof, nEl, metricErr, elemRel, nullMax, interMax, ok?"PASS":"FAIL");
    printf("    perf: %.3f ms/apply | %.1f MDOF/s | %.1f GB/s (useful traffic est)\n",
           sPerApply*1e3, mdofs, gbs);
    // record the sweep row: matrix-free bytes/DOF (what the matvec reads) vs the
    // assembled CSR it avoids ((2p+1)^3 nnz/row x 12B = value+colInd).
    double mfBpd  = bytesPerElem * (double)nEl / (double)nDof;
    double asmBpd = std::pow(2.0*P+1.0, 3.0) * 12.0;
    g_sweep.push_back({P, nDof, sPerApply*1e3, mdofs, gbs, mfBpd, asmBpd});

    cudaFree(d_elemDof); cudaFree(d_corners); cudaFree(d_G);
    cudaFree(d_u); cudaFree(d_y);
    cudaFree(d_edof1); cudaFree(d_u1); cudaFree(d_y1);
    return ok;
}

// ---- Gate (D): NON-UNIFORM (per-element sheared) metric + apply. ----
// The uniform cube cannot catch a per-element metric base-offset error (every
// element shares the same G) nor distinguish the cross terms g0/g1 from zero.
// Here each element gets a distinct straight-sided shear so detJ and g0/g1 are
// nonzero and element-distinct. We compare the GPU PerPoint metric AND the GPU
// assembled apply element-by-element against the host computeElementMetric /
// applyHoCvfemElement. This is what actually validates the PerPoint metric port
// the header advertises for curved/sheared hexes, and it would fail loudly on
// the factor-of-3 per-element stride bug.
template<int P>
static bool runShearGate(int E)
{
    auto op = buildHoCvfemOperators(P);
    const int n = P + 1, N3 = n * n * n;
    CubeMesh m = buildCube(op, P, E);
    const size_t nEl = m.nEl; const long nDof = m.nDof;

    // Per-element shear: x += a*z, y += b*z with element-distinct (a,b). Keeps
    // hexes straight-sided but makes J (and thus G) non-diagonal and per-element.
    for (size_t e = 0; e < nEl; ++e) {
        double a = 0.13 + 0.01 * (double)(e % 7);
        double b = 0.07 + 0.01 * (double)(e % 5);
        for (int c = 0; c < 8; ++c) {
            double z = m.h_corners[e*24 + c*3 + 2];
            m.h_corners[e*24 + c*3 + 0] += a * z;
            m.h_corners[e*24 + c*3 + 1] += b * z;
        }
    }

    CK(ho_cvfem_upload_operators(P, op.Btil.data(), op.Dtil.data(),
                                 op.D.data(), op.W.data(), op.xi.data(), op.zeta.data()));

    int* d_elemDof; double* d_corners; double* d_G; double* d_u; double* d_y;
    const size_t gLen = nEl * (size_t)(3 * P * n * n) * 3;
    CK(cudaMalloc(&d_elemDof, sizeof(int)    * nEl * N3));
    CK(cudaMalloc(&d_corners, sizeof(double) * nEl * 24));
    CK(cudaMalloc(&d_G,       sizeof(double) * gLen));
    CK(cudaMalloc(&d_u,       sizeof(double) * nDof));
    CK(cudaMalloc(&d_y,       sizeof(double) * nDof));
    CK(cudaMemcpy(d_elemDof, m.dh.elemDof.data(), sizeof(int) * nEl * N3, cudaMemcpyHostToDevice));
    CK(cudaMemcpy(d_corners, m.h_corners.data(),  sizeof(double) * nEl * 24, cudaMemcpyHostToDevice));
    CK(ho_cvfem_metric_perpoint_launch<double, P>(d_corners, d_G, nEl));
    CK(cudaDeviceSynchronize());

    // Metric check on element 0 AND the last element (per-element stride probe).
    const size_t perElem = (size_t)(3 * P * n * n);
    std::vector<double> Gall(gLen);
    CK(cudaMemcpy(Gall.data(), d_G, sizeof(double)*gLen, cudaMemcpyDeviceToHost));
    double metricErr = 0;
    for (size_t e : {(size_t)0, nEl - 1}) {
        double cor[8][3];
        for (int c=0;c<8;++c) for (int d=0;d<3;++d) cor[c][d] = m.h_corners[e*24 + c*3 + d];
        auto Gh = computeElementMetric(op, cor);
        for (size_t i=0;i<Gh.size();++i) for (int c=0;c<3;++c)
            metricErr = std::max(metricErr, std::abs(Gall[(e*perElem + i)*3 + c] - Gh[i][c]));
    }

    // Assembled apply vs host loop over applyHoCvfemElement (each element's own G).
    std::mt19937 rng(987); std::uniform_real_distribution<double> dist(-1.0,1.0);
    std::vector<double> uGlob(nDof); for (long i=0;i<nDof;++i) uGlob[i]=dist(rng);
    std::vector<double> yHost(nDof, 0.0);
    for (size_t e=0;e<nEl;++e) {
        double cor[8][3];
        for (int c=0;c<8;++c) for (int d=0;d<3;++d) cor[c][d] = m.h_corners[e*24 + c*3 + d];
        auto Gh = computeElementMetric(op, cor);
        std::vector<double> ul(N3), yl;
        for (int l=0;l<N3;++l) { int dof=m.dh.elemDof[e*N3+l]; ul[l]= dof>=0 ? uGlob[dof] : 0.0; }
        applyHoCvfemElement(op, Gh, ul, yl);
        for (int l=0;l<N3;++l) { int dof=m.dh.elemDof[e*N3+l]; if(dof>=0) yHost[dof]+=yl[l]; }
    }
    CK(cudaMemcpy(d_u, uGlob.data(), sizeof(double)*nDof, cudaMemcpyHostToDevice));
    CK(cudaMemset(d_y, 0, sizeof(double)*nDof));
    CK(ho_cvfem_apply_launch<double, P>(d_u, d_y, d_elemDof, d_G, nEl));
    CK(cudaDeviceSynchronize());
    std::vector<double> yGpu(nDof);
    CK(cudaMemcpy(yGpu.data(), d_y, sizeof(double)*nDof, cudaMemcpyDeviceToHost));
    double ynorm=0, applyErr=0;
    for (long i=0;i<nDof;++i){ ynorm=std::max(ynorm,std::abs(yHost[i]));
        applyErr=std::max(applyErr,std::abs(yGpu[i]-yHost[i])); }
    double applyRel = applyErr / (ynorm>0?ynorm:1.0);

    bool ok = (metricErr < 1e-12) && (applyRel < 1e-11);
    printf("p=%d E=%d (SHEAR) nEl=%zu | metricErr=%.2e applyRel=%.2e | %s\n",
           P, E, nEl, metricErr, applyRel, ok?"PASS":"FAIL");

    cudaFree(d_elemDof); cudaFree(d_corners); cudaFree(d_G); cudaFree(d_u); cudaFree(d_y);
    return ok;
}

int main()
{
    int dev=0; cudaGetDeviceCount(&dev); if (dev>0) cudaSetDevice(0);
    bool ok = true;
    // Larger meshes so the timing loop saturates the GPU (132 SMs on H100); the
    // tiny-mesh numbers are launch/latency-bound and not representative. The
    // correctness gates A/B/C are still cheap at these sizes (element-0 + max
    // reductions, no per-element host loop).
    // Full order sweep p=1..7. E sized per p for GPU saturation + memory fit
    // (DOFs=(E*p+1)^3, metric<150MB each). Structured cube == block-mesh geometry.
    ok &= runOrder<1>(48);
    ok &= runOrder<2>(40);
    ok &= runOrder<3>(28);
    ok &= runOrder<4>(24);
    ok &= runOrder<5>(20);
    ok &= runOrder<6>(18);
    ok &= runOrder<7>(16);
    // Non-uniform metric gates (catch per-element stride + cross-term port).
    ok &= runShearGate<1>(4);
    ok &= runShearGate<2>(3);
    ok &= runShearGate<4>(2);

    printf("\n=== MATRIX-FREE ORDER SWEEP (throughput + memory; structured cube = block-mesh geometry) ===\n");
    printf("   p |    DOFs    | ms/apply | MDOF/s |  GB/s | matfree B/DOF | assembled B/DOF | asm/mf\n");
    for (auto& r : g_sweep)
        printf("  %2d | %9ld | %8.3f | %6.0f | %5.0f | %13.1f | %15.1f | %6.1fx\n",
               r.p, r.dofs, r.ms, r.mdofs, r.gbs, r.mfBpd, r.asmBpd, r.asmBpd / r.mfBpd);

    printf("\nHO-CVFEM matrix-free GPU gate: %s\n", ok?"PASS":"FAIL");
    return ok?0:1;
}
