// Stage-0 GPU gate for the high-order tet Laplacian (collapsed sum-factorization).
// Replicates one tet across numElements, applies the device kernel, and checks it against the
// host reference applyHoLaplacianTetElement (the Mac-validated oracle) plus A.1=0. Tet twin of
// mars_ho_laplacian_test.cu. Single rank, no MPI -- pure element-local kernel parity.
#include "backend/distributed/unstructured/fem/mars_ho_laplacian_tet.hpp"
#include <cuda_runtime.h>
#include <cstdio>
#include <vector>
#include <cmath>
using namespace mars::fem;

#define CK(x) do{ cudaError_t e_=(x); if(e_!=cudaSuccess){ \
    std::printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); return false; } }while(0)

template <typename T> static T* up(const std::vector<T>& h){ T* d; cudaMalloc(&d,h.size()*sizeof(T));
    cudaMemcpy(d,h.data(),h.size()*sizeof(T),cudaMemcpyHostToDevice); return d; }

template <int P>
static bool run_order(int numElements) {
    auto o = buildHoTetOps<double>(P);
    const int W=o.W, nq=o.nq, W3=W*W*W, M=o.nModes;
    double corners[4][3]={{0.1,0.2,0.3},{2.1,0.3,0.4},{0.6,1.7,0.1},{0.4,0.5,1.5}};
    auto Gc = computeHoTetMetric<double>(corners, o);

    // per-element random modal field + host-reference output (the oracle)
    std::vector<double> u((size_t)numElements*W3, 0.0), yhost((size_t)numElements*W3, 0.0);
    for (int e=0;e<numElements;++e) {
        for (int m=0;m<M;++m){ auto&mm=o.modes[m]; u[(size_t)e*W3 + ((size_t)mm[0]*W+mm[1])*W+mm[2]] = std::sin(0.9*m+0.4+0.013*e); }
        applyHoLaplacianTetElement<double>(o, Gc.data(), &u[(size_t)e*W3], &yhost[(size_t)e*W3]);
    }
    std::vector<double> Gall((size_t)numElements*nq*6);
    for (int e=0;e<numElements;++e) for (size_t i=0;i<(size_t)nq*6;++i) Gall[(size_t)e*nq*6+i]=Gc[i];

    double *d_u=up(u), *d_A=up(o.A), *d_Ad=up(o.Ad), *d_B=up(o.B), *d_Bd=up(o.Bd),
           *d_C=up(o.C), *d_Cd=up(o.Cd), *d_G=up(Gall), *d_y=nullptr;
    CK(cudaMalloc(&d_y, u.size()*sizeof(double)));

    const int bs=128, nb=(numElements+bs-1)/bs;
    ho_laplacian_tet_apply_scalar<double,P><<<nb,bs>>>(d_u,d_y,d_A,d_Ad,d_B,d_Bd,d_C,d_Cd,d_G,numElements);
    CK(cudaGetLastError()); CK(cudaDeviceSynchronize());

    std::vector<double> ygpu(u.size()); CK(cudaMemcpy(ygpu.data(),d_y,u.size()*sizeof(double),cudaMemcpyDeviceToHost));
    double err=0; for (int e=0;e<numElements;++e) for (int m=0;m<M;++m){ size_t x=(size_t)e*W3+((size_t)o.modes[m][0]*W+o.modes[m][1])*W+o.modes[m][2]; err=std::max(err,std::fabs(ygpu[x]-yhost[x])); }

    // A.1 = 0 : element 0, constant mode
    std::vector<double> one((size_t)numElements*W3,0.0); for(int e=0;e<numElements;++e) one[(size_t)e*W3]=1.0;
    double* d_one=up(one); ho_laplacian_tet_apply_scalar<double,P><<<nb,bs>>>(d_one,d_y,d_A,d_Ad,d_B,d_Bd,d_C,d_Cd,d_G,numElements);
    CK(cudaDeviceSynchronize()); CK(cudaMemcpy(ygpu.data(),d_y,u.size()*sizeof(double),cudaMemcpyDeviceToHost));
    double a1=0; for(int m=0;m<M;++m) a1=std::max(a1,std::fabs(ygpu[((size_t)o.modes[m][0]*W+o.modes[m][1])*W+o.modes[m][2]]));

    cudaFree(d_u);cudaFree(d_y);cudaFree(d_A);cudaFree(d_Ad);cudaFree(d_B);cudaFree(d_Bd);cudaFree(d_C);cudaFree(d_Cd);cudaFree(d_G);cudaFree(d_one);
    bool ok = (err<1e-10) && (a1<1e-9);
    std::printf(" p=%d | modes=%3d | elems=%d | max|gpu-host|=%.2e | A.1=%.2e | %s\n", P, M, numElements, err, a1, ok?"PASS":"FAIL");
    return ok;
}

int main() {
    std::printf("=== HO tet Laplacian Stage-0 GPU parity (collapsed sum-fac) ===\n");
    bool ok = true;
    ok &= run_order<2>(2000);
    ok &= run_order<3>(2000);
    ok &= run_order<4>(1000);
    ok &= run_order<5>(400);
    std::printf(ok ? "\n==== ALL TET HO GPU PARITY CHECKS PASS ====\n" : "\n==== FAILURE ====\n");
    return ok ? 0 : 1;
}
