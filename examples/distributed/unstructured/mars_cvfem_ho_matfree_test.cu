// GPU correctness + perf gate for the optimized HO-CVFEM matrix-free diffusion
// apply (mars_cvfem_ho_matfree.hpp). Single rank, in-memory structured cube --
// same elemDof/coord convention as the host patch test, so any divergence
// localizes to the GPU kernel and not the DOF map. The numbering itself is built
// on the DEVICE (buildGpu) and stays there -- the apply reads elemDof straight
// out of HoOwnershipDeviceData, with no host build and no H2D.
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
//
// --dof-self-check (or MARS_HO_DOF_SELF_CHECK=1) adds a numbering gate in front:
// host HODofHandler::build() vs the single-rank GPU buildGpu(), compared on the
// permutation-invariant quantities. Off by default; the gates below are unchanged.
// --dof-self-check=E (or MARS_HO_DOF_SELF_CHECK_E=E) sets the cube size. E=4 is
// the cheap correctness cube; a large E is what shows the setup-time win, since
// the host numbering goes through std::map and the GPU path through a radix
// dedup. Both build times are printed per p.

#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler_gpu.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_matfree.hpp"

#include <cuda_runtime.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
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
    // The numbering is built on the device and STAYS there: own.elemDof is what
    // the apply reads, so there is no H2D upload of it anywhere below.
    HoOwnershipDeviceData own;
    std::vector<std::array<int,3>> ijk;
    long nDof;
    size_t nEl;
    int n, N3;
    std::vector<double> h_corners;   // [nEl*8*3]
};

static void makeCubeCorners(int E, std::vector<std::array<int,8>>& ec,
                            std::vector<std::array<int,3>>& ijk)
{
    auto cg = [&](int x, int y, int z) { return (x*(E+1)+y)*(E+1)+z; };
    for (int ex=0; ex<E; ++ex) for (int ey=0; ey<E; ++ey) for (int ez=0; ez<E; ++ez) {
        ec.push_back({cg(ex,ey,ez),cg(ex+1,ey,ez),cg(ex+1,ey+1,ez),cg(ex,ey+1,ez),
                      cg(ex,ey,ez+1),cg(ex+1,ey,ez+1),cg(ex+1,ey+1,ez+1),cg(ex,ey+1,ez+1)});
        ijk.push_back({ex,ey,ez});
    }
}

static CubeMesh buildCube(const HoCvfemOperators& op, int P, int E)
{
    CubeMesh m;
    const int n = P + 1, N3 = n * n * n;
    std::vector<std::array<int,8>> ec;
    makeCubeCorners(E, ec, m.ijk);
    buildGpu(m.dh, ec, long(E+1)*(E+1)*(E+1), P, &m.own);
    m.nDof = m.dh.numDof; m.nEl = ec.size(); m.n = n; m.N3 = N3;
    // TEST ONLY: the gates below score the device apply against a host reference,
    // which indexes elemDof on the host. Passing keepOwn deliberately skips that
    // download, so take it once here. Nothing in the apply path needs it.
    m.dh.elemDof.resize((size_t)m.nEl * N3);
    thrust::copy(m.own.elemDof.begin(), m.own.elemDof.end(), m.dh.elemDof.begin());

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
    double* d_corners; double* d_G;
    const int* d_elemDof = thrust::raw_pointer_cast(m.own.elemDof.data());
    double* d_u; double* d_y;
    const size_t gLen = nEl * (size_t)(3 * P * n * n) * 3;
    CK(cudaMalloc(&d_corners, sizeof(double) * nEl * 24));
    CK(cudaMalloc(&d_G,       sizeof(double) * gLen));
    CK(cudaMalloc(&d_u,       sizeof(double) * nDof));
    CK(cudaMalloc(&d_y,       sizeof(double) * nDof));
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

    cudaFree(d_corners); cudaFree(d_G);
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

    double* d_corners; double* d_G; double* d_u; double* d_y;
    const int* d_elemDof = thrust::raw_pointer_cast(m.own.elemDof.data());
    const size_t gLen = nEl * (size_t)(3 * P * n * n) * 3;
    CK(cudaMalloc(&d_corners, sizeof(double) * nEl * 24));
    CK(cudaMalloc(&d_G,       sizeof(double) * gLen));
    CK(cudaMalloc(&d_u,       sizeof(double) * nDof));
    CK(cudaMalloc(&d_y,       sizeof(double) * nDof));
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

    cudaFree(d_corners); cudaFree(d_G); cudaFree(d_u); cudaFree(d_y);
    return ok;
}

// ---- Opt-in gate: single-rank host build() vs the GPU buildGpu(). ----
// --dof-self-check (or MARS_HO_DOF_SELF_CHECK=1). Additive: the A/B/C/shear gates
// below run exactly as before with or without it.
//
// The GPU numbers edges/faces in sorted-key order, the host in std::map insertion
// order, so the two DOF numberings are a PERMUTATION of each other by construction.
// Only permutation-invariant properties can be compared:
//   (1) numDof / nEdge / nFace,
//   (2) the multiset of canonical DofKeys,
//   (3) the elemDof identification classes: the (element,local-node) slots that
//       share a DOF must be the same set on both sides. This is the strongest of
//       the three -- it proves ONE bijection relates the two maps, i.e. the two
//       spaces have identical continuity, not just the same keys somewhere.
// (1)+(2) are what mars_ho_dist_apply_test --self-check compares for the
// distributed path; (3) is available here because there is a single numbering.
static bool checkDofNumbering(int P, int E)
{
    std::vector<std::array<int,8>> ec;
    std::vector<std::array<int,3>> ijk;
    makeCubeCorners(E, ec, ijk);
    const long nCornerNodes = (long)(E+1)*(E+1)*(E+1);
    const long nElem = (long)ec.size();
    const int  n = P + 1, N3 = n*n*n;

    // The degenerate single-rank configuration buildGpu feeds the shared numbering
    // core, spelled out here so the host oracle sees the identical inputs.
    std::vector<long> cornerGid(nCornerNodes);
    for (long i = 0; i < nCornerNodes; ++i) cornerGid[i] = i;   // global id == local id
    const std::vector<int>     cornerOwner(nCornerNodes, 0);
    const std::vector<uint8_t> sharedCorner(nCornerNodes, 0);
    const std::vector<int>     elemOwner(nElem, 0);

    HODofHandler dofPlain, dofHost, dofGpu;
    dofPlain.build(ec, nCornerNodes, P);
    using Clk = std::chrono::steady_clock;
    // buildDistributed() calls build() and then tags the DofKeys build() alone does
    // not produce. Comparing dofPlain to dofHost verifies that rather than assuming
    // it, so the oracle really is build().
    const auto t0 = Clk::now();
    dofHost.buildDistributed(ec, nCornerNodes, P, cornerGid, cornerOwner, elemOwner, 0, sharedCorner);
    const auto t1 = Clk::now();
    buildGpu(dofGpu, ec, nCornerNodes, P);
    cudaDeviceSynchronize();
    const auto t2 = Clk::now();
    const double hostS = std::chrono::duration<double>(t1 - t0).count();
    const double gpuS  = std::chrono::duration<double>(t2 - t1).count();

    const bool oracleOk = (dofPlain.numDof == dofHost.numDof) && (dofPlain.elemDof == dofHost.elemDof);
    const bool countOk  = (dofHost.numDof == dofGpu.numDof) && (dofHost.nEdge == dofGpu.nEdge) &&
                          (dofHost.nFace == dofGpu.nFace);

    // (2) DofKey multiset.
    long keyBad = -1;
    if (dofHost.dofKey.size() == dofGpu.dofKey.size()) {
        auto packed = [](const HODofHandler::DofKey& k) {
            return std::array<long,6>{ (long)k.kind, k.g0, k.g1, k.g2, k.g3, (long)k.pos }; };
        std::vector<std::array<long,6>> hk(dofHost.dofKey.size()), gk(dofGpu.dofKey.size());
        for (size_t i = 0; i < hk.size(); ++i) hk[i] = packed(dofHost.dofKey[i]);
        for (size_t i = 0; i < gk.size(); ++i) gk[i] = packed(dofGpu.dofKey[i]);
        std::sort(hk.begin(), hk.end());
        std::sort(gk.begin(), gk.end());
        keyBad = 0;
        for (size_t i = 0; i < hk.size(); ++i) if (hk[i] != gk[i]) ++keyBad;
    }

    // (3) elemDof identification classes: host dof <-> GPU dof must be a bijection.
    long permBad = 0, unmapped = 0;
    bool identity = true;
    if (countOk && dofHost.elemDof.size() == dofGpu.elemDof.size()) {
        std::vector<int> h2g(dofHost.numDof, -1), g2h(dofGpu.numDof, -1);
        for (size_t s = 0; s < dofHost.elemDof.size(); ++s) {
            const int hd = dofHost.elemDof[s], gd = dofGpu.elemDof[s];
            if (hd != gd) identity = false;
            if (h2g[hd] < 0) h2g[hd] = gd; else if (h2g[hd] != gd) ++permBad;
            if (g2h[gd] < 0) g2h[gd] = hd; else if (g2h[gd] != hd) ++permBad;
        }
        for (long d = 0; d < dofHost.numDof; ++d) if (h2g[d] < 0) ++unmapped;
    } else {
        permBad = -1;
    }

    // With no shared corners nothing may be flagged shared or on a boundary, on
    // either side -- a direct check that the degenerate inputs landed as intended.
    long hShared = 0, gShared = 0, hBnd = 0, gBnd = 0;
    for (auto v : dofHost.dofShared)   hShared += v;
    for (auto v : dofGpu.dofShared)    gShared += v;
    for (auto v : dofHost.dofBoundary) hBnd += v;
    for (auto v : dofGpu.dofBoundary)  gBnd += v;
    const bool flagOk = (hShared == 0 && gShared == 0 && hBnd == 0 && gBnd == 0);

    // Single rank -> every DOF is owned by rank 0 on both sides.
    bool ownOk = ((long)dofHost.dofOwner.size() == dofHost.numDof) &&
                 ((long)dofGpu.dofOwner.size() == dofGpu.numDof);
    for (int o : dofHost.dofOwner) if (o != 0) ownOk = false;
    for (int o : dofGpu.dofOwner)  if (o != 0) ownOk = false;

    const bool ok = oracleOk && countOk && keyBad == 0 && permBad == 0 &&
                    unmapped == 0 && flagOk && ownOk;
    printf("[dof-self-check] p=%d E=%d nEl=%ld N3=%d | numDof h=%ld g=%ld nEdge h=%ld g=%ld nFace h=%ld g=%ld"
           " | build host %.3fs gpu %.3fs (%.2fx)\n",
           P, E, nElem, N3, dofHost.numDof, dofGpu.numDof, dofHost.nEdge, dofGpu.nEdge, dofHost.nFace, dofGpu.nFace,
           hostS, gpuS, gpuS > 0.0 ? hostS / gpuS : 0.0);
    printf("                 oracle=%s keyMismatch=%ld permMismatch=%ld unmapped=%ld flags=%s owner=%s perm=%s | %s\n",
           oracleOk ? "ok" : "BAD", keyBad, permBad, unmapped,
           flagOk ? "ok" : "BAD", ownOk ? "ok" : "BAD",
           identity ? "identity" : "permuted", ok ? "PASS" : "FAIL");
    return ok;
}

int main(int argc, char** argv)
{
    int dev=0; cudaGetDeviceCount(&dev); if (dev>0) cudaSetDevice(0);
    bool ok = true;

    // Opt-in only: proves host build() and the single-rank GPU buildGpu() number the
    // same space before the operator gates run on the host numbering.
    bool dofSelfCheck = std::getenv("MARS_HO_DOF_SELF_CHECK") != nullptr;
    int  selfCheckE = 4;   // cells per side; 4 is the cheap correctness cube
    if (const char* ev = std::getenv("MARS_HO_DOF_SELF_CHECK_E")) {
        const int v = std::atoi(ev);
        if (v > 0) { selfCheckE = v; dofSelfCheck = true; }
    }
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--dof-self-check") == 0) {
            dofSelfCheck = true;
        } else if (std::strncmp(argv[i], "--dof-self-check=", 17) == 0) {
            dofSelfCheck = true;
            const int v = std::atoi(argv[i] + 17);
            if (v > 0) selfCheckE = v;
        }
    }
    if (dofSelfCheck) {
        printf("=== HO DOF numbering self-check: host build() vs GPU buildGpu() (E=%d) ===\n",
               selfCheckE);
        for (int p = 1; p <= 7; ++p) ok &= checkDofNumbering(p, selfCheckE);
        printf("\n");
    }
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
