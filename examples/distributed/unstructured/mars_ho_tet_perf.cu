// Single-GPU performance sweep for the high-order TET operators -- the tet counterpart of the
// hex HO CVFEM talk numbers (throughput per GPU vs order, bytes per DOF).
//   1. tet CVFEM smem sum-fac kernel (block/element)  -- parity-gated vs host, then timed
//   2. tet CVFEM dense correctness kernel             -- baseline
//   3. tet Galerkin stored-G scalar kernel            -- the stored-metric contrast (PA-style)
// CG semantics on a sorted Kuhn mesh (merged collapsed DOF, gather/scatter + atomics).
//
//   srun -N1 -n1 ./examples/distributed/unstructured/mars_ho_tet_perf --ncells=12 --p=4 --reps=50
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_tet.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_laplacian_tet.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler_tet.hpp"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <cmath>
using namespace mars::fem;

#define CK(x) do{ cudaError_t e_=(x); if(e_!=cudaSuccess){ \
    std::printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); std::exit(1);} }while(0)

template <typename T> static T* up(const std::vector<T>& h){ T* d; CK(cudaMalloc(&d,h.size()*sizeof(T)));
    CK(cudaMemcpy(d,h.data(),h.size()*sizeof(T),cudaMemcpyHostToDevice)); return d; }

// GEMM pipeline for the universal-matrix path: gather element blocks, V = K_stack * U (cuBLAS,
// column-major), then per-element coefficient combine + scatter.
template <typename RealType>
__global__ void tet_gather_kernel(const RealType* __restrict__ u, const int* __restrict__ dof,
                                  RealType* __restrict__ U, size_t total) {
    const size_t i = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (i < total) U[i] = u[dof[i]];
}
template <typename RealType, int NN>
__global__ void tet_combine_scatter(const RealType* __restrict__ V, const int* __restrict__ dof,
                                    const RealType* __restrict__ geom, RealType* __restrict__ y, size_t nE) {
    const size_t w = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (w >= nE*(size_t)NN) return;
    const size_t e = w / NN; const int i = (int)(w % NN);
    const RealType* MMt = geom + e*10; const RealType detB = geom[e*10+9];
    const RealType* Ve = V + e*(size_t)6*NN;
    const RealType acc = detB*(MMt[0]*Ve[0*NN+i] + MMt[1]*Ve[1*NN+i] + MMt[2]*Ve[2*NN+i]
                             + MMt[4]*Ve[3*NN+i] + MMt[5]*Ve[4*NN+i] + MMt[8]*Ve[5*NN+i]);
    atomicAdd(&y[dof[w]], acc);
}

template <int P>
static void run_order(int ncells, int reps) {
    constexpr int M = (P + 3) / 2;
    auto o = buildHoCvfemTetOps<double>(P);
    const int n = o.n, NN = o.NN;

    // ---- parity gate on a small mesh: smem kernel vs host sum-fac assembled ----
    {
        std::vector<std::array<double,3>> coords; std::vector<int> ec;
        buildKuhnTetMesh(2, coords, ec);
        const size_t nE = ec.size()/4;
        std::vector<double> Zd(o.Z.begin(), o.Z.end());
        HoCvfemTetDofHandler dh; dh.build(ec, coords, Zd);
        const long nDof = dh.numDof;
        std::vector<double> geom(nE*10);
        for (size_t e=0;e<nE;++e){ double c4[4][3];
            for (int c=0;c<4;++c) for (int d=0;d<3;++d) c4[c][d]=coords[dh.sortedCorners[e*4+c]][d];
            double MMt[9], detB; computeTetGeom<double>(c4, MMt, detB);
            for (int i=0;i<9;++i) geom[e*10+i]=MMt[i]; geom[e*10+9]=detB; }
        std::vector<double> u(nDof), yref(nDof,0.0), ue(NN), ye(NN);
        for (long i=0;i<nDof;++i) u[i]=std::sin(0.9*i+0.4);
        for (size_t e=0;e<nE;++e){
            for (int i=0;i<NN;++i) ue[i]=u[dh.elemDof[e*NN+i]];
            applyHoCvfemTetElement<double>(o, &geom[e*10], geom[e*10+9], ue.data(), ye.data());
            for (int i=0;i<NN;++i) yref[dh.elemDof[e*NN+i]] += ye[i]; }
        std::vector<double> hxs(o.xs.begin(),o.xs.end()), hyq(o.yq.begin(),o.yq.end()), hwq(o.wqv.begin(),o.wqv.end()),
                            hLs(o.Ls.begin(),o.Ls.end()), hDs(o.Ds.begin(),o.Ds.end()),
                            hLq(o.Lq.begin(),o.Lq.end()), hDq(o.Dq.begin(),o.Dq.end());
        double *d_xs=up(hxs),*d_yq=up(hyq),*d_wq=up(hwq),*d_Ls=up(hLs),*d_Ds=up(hDs),*d_Lq=up(hLq),*d_Dq=up(hDq),
               *d_geom=up(geom),*d_u=up(u),*d_y; int *d_ed=up(dh.elemDof);
        CK(cudaMalloc(&d_y,nDof*sizeof(double))); CK(cudaMemset(d_y,0,nDof*sizeof(double)));
        ho_cvfem_tet_apply_cg_smem<double,P,M><<<(unsigned)nE,256>>>(d_u,d_y,d_ed,d_xs,d_yq,d_wq,d_Ls,d_Ds,d_Lq,d_Dq,d_geom,nE);
        CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
        std::vector<double> yg(nDof); CK(cudaMemcpy(yg.data(),d_y,nDof*sizeof(double),cudaMemcpyDeviceToHost));
        double err=0; for (long i=0;i<nDof;++i) err=std::max(err,std::fabs(yg[i]-yref[i]));
        std::printf(" p=%d  smem-kernel parity vs host: %.2e %s\n", P, err, err<1e-10?"PASS":"FAIL");
        // universal-matrix (affine fast path) parity
        auto K6 = buildCvfemTetRefMatrices<double>(o);
        double* d_K = up(K6);
        CK(cudaMemset(d_y,0,nDof*sizeof(double)));
        ho_cvfem_tet_apply_cg_univ<double,P><<<(unsigned)nE,128>>>(d_u,d_y,d_ed,d_K,d_geom,nE);
        CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
        CK(cudaMemcpy(yg.data(),d_y,nDof*sizeof(double),cudaMemcpyDeviceToHost));
        double erru=0; for (long i=0;i<nDof;++i) erru=std::max(erru,std::fabs(yg[i]-yref[i]));
        std::printf(" p=%d  univ-kernel parity vs host: %.2e %s\n", P, erru, erru<1e-10?"PASS":"FAIL");
        cudaFree(d_K);
        cudaFree(d_xs);cudaFree(d_yq);cudaFree(d_wq);cudaFree(d_Ls);cudaFree(d_Ds);cudaFree(d_Lq);cudaFree(d_Dq);
        cudaFree(d_geom);cudaFree(d_u);cudaFree(d_y);cudaFree(d_ed);
        if (err >= 1e-10 || erru >= 1e-10) return;
    }

    // ---- timing mesh ----
    std::vector<std::array<double,3>> coords; std::vector<int> ec;
    buildKuhnTetMesh(ncells, coords, ec);
    const size_t nE = ec.size()/4;
    std::vector<double> Zd(o.Z.begin(), o.Z.end());
    HoCvfemTetDofHandler dh; dh.build(ec, coords, Zd);
    const long nDof = dh.numDof;
    std::vector<double> geom(nE*10);
    for (size_t e=0;e<nE;++e){ double c4[4][3];
        for (int c=0;c<4;++c) for (int d=0;d<3;++d) c4[c][d]=coords[dh.sortedCorners[e*4+c]][d];
        double MMt[9], detB; computeTetGeom<double>(c4, MMt, detB);
        for (int i=0;i<9;++i) geom[e*10+i]=MMt[i]; geom[e*10+9]=detB; }

    std::vector<double> hxs(o.xs.begin(),o.xs.end()), hyq(o.yq.begin(),o.yq.end()), hwq(o.wqv.begin(),o.wqv.end()),
                        hLs(o.Ls.begin(),o.Ls.end()), hDs(o.Ds.begin(),o.Ds.end()),
                        hLq(o.Lq.begin(),o.Lq.end()), hDq(o.Dq.begin(),o.Dq.end());
    double *d_xs=up(hxs),*d_yq=up(hyq),*d_wq=up(hwq),*d_Ls=up(hLs),*d_Ds=up(hDs),*d_Lq=up(hLq),*d_Dq=up(hDq),
           *d_geom=up(geom);
    int *d_ed=up(dh.elemDof);
    std::vector<double> u(nDof); for (long i=0;i<nDof;++i) u[i]=std::sin(0.9*i+0.4);
    double *d_u=up(u), *d_y; CK(cudaMalloc(&d_y,nDof*sizeof(double)));

    auto timeit=[&](auto&& launch, int rep)->double{
        launch(); CK(cudaDeviceSynchronize());                       // warmup
        cudaEvent_t t0,t1; cudaEventCreate(&t0); cudaEventCreate(&t1);
        cudaEventRecord(t0);
        for (int r=0;r<rep;++r) launch();
        cudaEventRecord(t1); CK(cudaEventSynchronize(t1));
        float ms; cudaEventElapsedTime(&ms,t0,t1);
        cudaEventDestroy(t0); cudaEventDestroy(t1);
        return ms/rep*1e-3;    // seconds per apply
    };

    // 1. CVFEM smem sum-fac
    const double t_smem = timeit([&]{
        cudaMemsetAsync(d_y,0,nDof*sizeof(double));
        ho_cvfem_tet_apply_cg_smem<double,P,M><<<(unsigned)nE,256>>>(d_u,d_y,d_ed,d_xs,d_yq,d_wq,d_Ls,d_Ds,d_Lq,d_Dq,d_geom,nE);
    }, reps);
    // 2. CVFEM universal-matrix via cuBLAS: V = K_stack U (column-major), combine + scatter.
    constexpr int NNC = (P+1)*(P+1)*(P+1);
    auto K6 = buildCvfemTetRefMatrices<double>(o);
    std::vector<double> Kst((size_t)6*NN*NN);
    for (int c=0;c<6;++c) for (int i=0;i<NN;++i) for (int j=0;j<NN;++j)
        Kst[(size_t)j*(6*NN) + (size_t)c*NN + i] = K6[((size_t)c*NN+i)*NN + j];
    double* d_Kst = up(Kst);
    double *d_U, *d_V;
    const size_t tot = nE*(size_t)NN;
    CK(cudaMalloc(&d_U, tot*sizeof(double)));
    CK(cudaMalloc(&d_V, tot*6*sizeof(double)));
    cublasHandle_t cbh; cublasCreate(&cbh);
    const double onec = 1.0, zeroc = 0.0;
    auto gemmApply=[&]{
        cudaMemsetAsync(d_y,0,nDof*sizeof(double));
        tet_gather_kernel<double><<<(unsigned)((tot+255)/256),256>>>(d_u,d_ed,d_U,tot);
        cublasDgemm(cbh, CUBLAS_OP_N, CUBLAS_OP_N, 6*NN, (int)nE, NN,
                    &onec, d_Kst, 6*NN, d_U, NN, &zeroc, d_V, 6*NN);
        tet_combine_scatter<double,NNC><<<(unsigned)((tot+255)/256),256>>>(d_V,d_ed,d_geom,d_y,nE);
    };
    // parity vs the (host-gated) smem kernel on this mesh
    {
        std::vector<double> ys(nDof), yg2(nDof);
        CK(cudaMemset(d_y,0,nDof*sizeof(double)));
        ho_cvfem_tet_apply_cg_smem<double,P,M><<<(unsigned)nE,256>>>(d_u,d_y,d_ed,d_xs,d_yq,d_wq,d_Ls,d_Ds,d_Lq,d_Dq,d_geom,nE);
        CK(cudaDeviceSynchronize());
        CK(cudaMemcpy(ys.data(),d_y,nDof*sizeof(double),cudaMemcpyDeviceToHost));
        gemmApply(); CK(cudaDeviceSynchronize());
        CK(cudaMemcpy(yg2.data(),d_y,nDof*sizeof(double),cudaMemcpyDeviceToHost));
        double e2=0; for (long i=0;i<nDof;++i) e2=std::max(e2,std::fabs(ys[i]-yg2[i]));
        std::printf(" p=%d  gemm parity vs smem: %.2e %s\n", P, e2, e2<1e-10?"PASS":"FAIL");
    }
    const double t_gemm = timeit(gemmApply, reps);
    // 3. CVFEM dense (baseline correctness kernel)
    const double t_dense = timeit([&]{
        cudaMemsetAsync(d_y,0,nDof*sizeof(double));
        ho_cvfem_tet_apply_cg<double,P,M><<<(unsigned)((nE+63)/64),64>>>(d_u,d_y,d_ed,d_xs,d_yq,d_wq,d_Ls,d_Ds,d_Lq,d_Dq,d_geom,nE);
    }, reps>10?10:reps);

    const double bytesPerDofCvfem = (nE*10.0*8 + nE*(double)NN*4 + 2.0*8*nDof) / nDof;
    std::printf(" p=%d | elems=%7zu | nDof=%9ld | gemm: %8.3f ms %7.2f GDOF/s | smem: %8.3f ms %5.2f GDOF/s | dense: %8.3f ms | %5.1f B/DOF\n",
                P, nE, nDof, t_gemm*1e3, nDof/t_gemm*1e-9, t_smem*1e3, nDof/t_smem*1e-9, t_dense*1e3, bytesPerDofCvfem);
    cublasDestroy(cbh);
    cudaFree(d_Kst); cudaFree(d_U); cudaFree(d_V);

    cudaFree(d_xs);cudaFree(d_yq);cudaFree(d_wq);cudaFree(d_Ls);cudaFree(d_Ds);cudaFree(d_Lq);cudaFree(d_Dq);
    cudaFree(d_geom);cudaFree(d_ed);cudaFree(d_u);cudaFree(d_y);

    // 3. Galerkin stored-G scalar kernel (modal blocks, PA-style stored metric contrast)
    {
        auto og = buildHoTetOps<double>(P);
        const int W = og.W, W3 = W*W*W, nq = og.nq, nMod = og.nModes;
        const size_t nEl = nE;
        double c4[4][3] = {{0.1,0.2,0.3},{2.1,0.3,0.4},{0.6,1.7,0.1},{0.4,0.5,1.5}};
        auto Gc1 = computeHoTetMetric<double>(c4, og);
        std::vector<double> Gall(nEl*(size_t)nq*6);
        for (size_t e=0;e<nEl;++e) std::memcpy(&Gall[e*(size_t)nq*6], Gc1.data(), Gc1.size()*sizeof(double));
        std::vector<double> ub(nEl*(size_t)W3, 0.0);
        for (size_t e=0;e<nEl;++e) for (int m=0;m<nMod;++m)
            ub[e*W3 + ((size_t)og.modes[m][0]*W+og.modes[m][1])*W+og.modes[m][2]] = std::sin(0.9*m+0.4);
        double *g_u=up(ub), *g_G=up(Gall), *g_A=up(og.A), *g_Ad=up(og.Ad), *g_B=up(og.B), *g_Bd=up(og.Bd),
               *g_C=up(og.C), *g_Cd=up(og.Cd), *g_y;
        CK(cudaMalloc(&g_y, ub.size()*sizeof(double)));
        const double t_gal = timeit([&]{
            ho_laplacian_tet_apply_scalar<double,P><<<(unsigned)((nEl+127)/128),128>>>(g_u,g_y,g_A,g_Ad,g_B,g_Bd,g_C,g_Cd,g_G,nEl);
        }, reps);
        const double dofGal = nEl*(double)nMod;
        const double bpdGal = (nEl*(double)nq*6*8 + 2.0*8*nEl*W3) / dofGal;
        std::printf("      |                                 | Galerkin stored-G: %8.3f ms  %7.2f GDOF/s | %5.1f B/DOF\n",
                    t_gal*1e3, dofGal/t_gal*1e-9, bpdGal);
        cudaFree(g_u);cudaFree(g_G);cudaFree(g_A);cudaFree(g_Ad);cudaFree(g_B);cudaFree(g_Bd);
        cudaFree(g_C);cudaFree(g_Cd);cudaFree(g_y);
    }
}

int main(int argc, char** argv) {
    int ncells=12, p=0, reps=50;
    for (int i=1;i<argc;++i){
        if (!std::strncmp(argv[i],"--ncells=",9)) ncells=std::atoi(argv[i]+9);
        if (!std::strncmp(argv[i],"--p=",4))      p=std::atoi(argv[i]+4);
        if (!std::strncmp(argv[i],"--reps=",7))   reps=std::atoi(argv[i]+7);
    }
    CK(cudaSetDevice(0));
    std::printf("=== HO tet perf sweep (ncells=%d, reps=%d) ===\n", ncells, reps);
    if (p==0 || p==2) run_order<2>(ncells,reps);
    if (p==0 || p==3) run_order<3>(ncells,reps);
    if (p==0 || p==4) run_order<4>(ncells,reps);
    if (p==0 || p==5) run_order<5>(ncells,reps);
    return 0;
}
