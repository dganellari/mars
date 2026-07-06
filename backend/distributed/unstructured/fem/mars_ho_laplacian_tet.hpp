#ifndef MARS_HO_LAPLACIAN_TET_HPP
#define MARS_HO_LAPLACIAN_TET_HPP

// High-order matrix-free Laplacian on tetrahedra via collapsed-coordinate (Duffy) sum-factorization.
// Tet twin of mars_ho_laplacian.hpp (hex). Element operator y = B^T G_c B u, PKD modal basis + the
// collapse metric. This is the exactly-consistent, symmetric Galerkin operator validated bit-exact
// against tet_ho_ref rung 9 -- NOT the CVFEM flux balance (which is only approximately consistent on
// the rational collapse metric, see TET_SUMFAC_HANDOFF.md).
//
// Modal DOF: (p,q,r) with p+q+r<=P, laid out u[(p*W+q)*W+r], W=P+1. The (W^3 - nModes) unused slots
// are a small waste kept for index simplicity; a compact packing is a later optimization.

#include "mars_ho_tet_basis.hpp"
#include <vector>
#include <array>
#include <cmath>

namespace mars {
namespace fem {

// Per-quad-point symmetric metric G_c (6 comps {G00,G01,G02,G11,G12,G22}, already * quad weight W_q),
// from the affine element corners and the collapse Jacobian J_c. Layout: Gc[point*6 + c].
template <typename RealType>
std::vector<RealType> computeHoTetMetric(const RealType corners[4][3], const HoTetOps<RealType>& o) {
    double B[9];
    for (int d = 0; d < 3; ++d) {
        B[d*3+0] = corners[1][d] - corners[0][d];
        B[d*3+1] = corners[2][d] - corners[0][d];
        B[d*3+2] = corners[3][d] - corners[0][d];
    }
    const double det = B[0]*(B[4]*B[8]-B[5]*B[7]) - B[1]*(B[3]*B[8]-B[5]*B[6]) + B[2]*(B[3]*B[7]-B[4]*B[6]);
    const double detB = std::fabs(det), inv = 1.0 / det;
    double M[9];   // Minv = B^-1 via cofactors
    M[0]=(B[4]*B[8]-B[5]*B[7])*inv; M[1]=(B[2]*B[7]-B[1]*B[8])*inv; M[2]=(B[1]*B[5]-B[2]*B[4])*inv;
    M[3]=(B[5]*B[6]-B[3]*B[8])*inv; M[4]=(B[0]*B[8]-B[2]*B[6])*inv; M[5]=(B[2]*B[3]-B[0]*B[5])*inv;
    M[6]=(B[3]*B[7]-B[4]*B[6])*inv; M[7]=(B[1]*B[6]-B[0]*B[7])*inv; M[8]=(B[0]*B[4]-B[1]*B[3])*inv;
    double MMt[9];
    for (int a=0;a<3;++a) for (int b=0;b<3;++b) { double v=0; for(int c=0;c<3;++c) v+=M[a*3+c]*M[b*3+c]; MMt[a*3+b]=v; }

    const int n = o.n;
    std::vector<RealType> Gc((size_t)o.nq * 6);
    for (int i=0;i<n;++i) for (int j=0;j<n;++j) for (int k=0;k<n;++k) {
        const double a=0.5*(1+o.ga.x[i]), b=0.5*(1+o.gb.x[j]), c=0.5*(1+o.gc.x[k]);
        const double Wq=(0.5*o.ga.w[i])*(0.25*o.gb.w[j])*(0.125*o.gc.w[k])*detB;
        const double r=a*(1-b)*(1-c), s=b*(1-c);
        const double D1=1.0-s-c, D2=1.0-c;
        const double Jc[9]={ 1.0/D1, r/(D1*D1), r/(D1*D1),  0.0, 1.0/D2, s/(D2*D2),  0.0,0.0,1.0 };
        double JM[9]; for(int p=0;p<3;++p)for(int t=0;t<3;++t){ double v=0; for(int u=0;u<3;++u) v+=Jc[p*3+u]*MMt[u*3+t]; JM[p*3+t]=v; }
        double G[9];  for(int p=0;p<3;++p)for(int t=0;t<3;++t){ double v=0; for(int u=0;u<3;++u) v+=JM[p*3+u]*Jc[t*3+u]; G[p*3+t]=Wq*v; }
        const size_t x=((size_t)i*n+j)*n+k;
        Gc[x*6+0]=G[0]; Gc[x*6+1]=G[1]; Gc[x*6+2]=G[2]; Gc[x*6+3]=G[4]; Gc[x*6+4]=G[5]; Gc[x*6+5]=G[8];
    }
    return Gc;
}

// y = B^T G_c B u. u,y sized W^3 (ragged-indexed). Host reference; the CUDA kernel mirrors this.
template <typename RealType>
void applyHoLaplacianTetElement(const HoTetOps<RealType>& o, const RealType* Gc,
                                const RealType* u, RealType* y) {
    const int P=o.P, W=o.W, n=o.n, nq=o.nq;
    const RealType *A=o.A.data(),*Ad=o.Ad.data(),*B=o.B.data(),*Bd=o.Bd.data(),*C=o.C.data(),*Cd=o.Cd.data();
    std::vector<std::array<RealType,3>> grad(nq,{0,0,0}), flux(nq,{0,0,0});

    auto sweep=[&](const RealType* Af,const RealType* Bf,const RealType* Cf,int comp){
        std::vector<RealType> f1((size_t)W*W*n,0), f2((size_t)W*n*n,0);
        for(int p=0;p<=P;++p)for(int q=0;q<=P-p;++q)for(int k=0;k<n;++k){ RealType v=0;
            for(int r=0;r<=P-p-q;++r) v+=u[((size_t)p*W+q)*W+r]*Cf[(((size_t)p*W+q)*W+r)*n+k]; f1[((size_t)p*W+q)*n+k]=v; }
        for(int p=0;p<=P;++p)for(int j=0;j<n;++j)for(int k=0;k<n;++k){ RealType v=0;
            for(int q=0;q<=P-p;++q) v+=Bf[((size_t)p*W+q)*n+j]*f1[((size_t)p*W+q)*n+k]; f2[((size_t)p*n+j)*n+k]=v; }
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){ RealType v=0;
            for(int p=0;p<=P;++p) v+=Af[(size_t)p*n+i]*f2[((size_t)p*n+j)*n+k]; grad[((size_t)i*n+j)*n+k][comp]=v; } };
    sweep(Ad,B,C,0); sweep(A,Bd,C,1); sweep(A,B,Cd,2);

    for(int x=0;x<nq;++x){ const RealType* g=grad[x].data(); const RealType* G=&Gc[(size_t)x*6];
        flux[x][0]=G[0]*g[0]+G[1]*g[1]+G[2]*g[2];
        flux[x][1]=G[1]*g[0]+G[3]*g[1]+G[4]*g[2];
        flux[x][2]=G[2]*g[0]+G[4]*g[1]+G[5]*g[2]; }

    for(size_t I=0;I<(size_t)W*W*W;++I) y[I]=0;
    auto tsweep=[&](const RealType* Af,const RealType* Bf,const RealType* Cf,int comp){
        std::vector<RealType> h1((size_t)W*n*n,0), h2((size_t)W*W*n,0);
        for(int p=0;p<=P;++p)for(int j=0;j<n;++j)for(int k=0;k<n;++k){ RealType v=0;
            for(int i=0;i<n;++i) v+=Af[(size_t)p*n+i]*flux[((size_t)i*n+j)*n+k][comp]; h1[((size_t)p*n+j)*n+k]=v; }
        for(int p=0;p<=P;++p)for(int q=0;q<=P-p;++q)for(int k=0;k<n;++k){ RealType v=0;
            for(int j=0;j<n;++j) v+=Bf[((size_t)p*W+q)*n+j]*h1[((size_t)p*n+j)*n+k]; h2[((size_t)p*W+q)*n+k]=v; }
        for(int p=0;p<=P;++p)for(int q=0;q<=P-p;++q)for(int r=0;r<=P-p-q;++r){ RealType v=0;
            for(int k=0;k<n;++k) v+=Cf[(((size_t)p*W+q)*W+r)*n+k]*h2[((size_t)p*W+q)*n+k]; y[((size_t)p*W+q)*W+r]+=v; } };
    tsweep(Ad,B,C,0); tsweep(A,Bd,C,1); tsweep(A,B,Cd,2);
}

// Nodal element apply for C0 assembly: y_nodal = V^-T ( K_modal ( V^-1 u_nodal ) ).
// u,y sized Np (= nModes); the ragged W^3 modal layout is internal.
template <typename RealType>
void applyHoLaplacianTetElementNodal(const HoTetOps<RealType>& o, const HoTetNodal<RealType>& nd,
                                     const RealType* Gc, const RealType* u, RealType* y) {
    const int Np = nd.Np, W = o.W;
    std::vector<RealType> um((size_t)W*W*W, 0), ym((size_t)W*W*W, 0), uh(Np), yh(Np);
    for (int m = 0; m < Np; ++m) { RealType v = 0;
        for (int nn = 0; nn < Np; ++nn) v += nd.Vinv[(size_t)m*Np+nn] * u[nn];
        uh[m] = v; }
    for (int m = 0; m < Np; ++m) { auto& e = o.modes[m]; um[((size_t)e[0]*W+e[1])*W+e[2]] = uh[m]; }
    applyHoLaplacianTetElement<RealType>(o, Gc, um.data(), ym.data());
    for (int m = 0; m < Np; ++m) { auto& e = o.modes[m]; yh[m] = ym[((size_t)e[0]*W+e[1])*W+e[2]]; }
    for (int nn = 0; nn < Np; ++nn) { RealType v = 0;
        for (int m = 0; m < Np; ++m) v += nd.Vinv[(size_t)m*Np+nn] * yh[m];   // Vinv^T
        y[nn] = v; }
}

// ---------------------------------------------------------------------------------------------------
// Device kernel: one thread per element, scalar reference (tet twin of hex ho_laplacian_apply_scalar).
// A direct transcription of applyHoLaplacianTetElement above -- correct by construction, to be
// GPU-validated on Alps. Local temporaries sized at compile time from P; heavy but for correctness,
// not speed (the shared-memory / multi-element optimized kernel is a follow-on).
#if defined(USE_CUDA)

template <typename RealType, int P>
__device__ inline void ho_tet_fwd_sweep(const RealType* __restrict__ u,
        const RealType* __restrict__ Af, const RealType* __restrict__ Bf, const RealType* __restrict__ Cf,
        int comp, RealType grad[][3]) {
    constexpr int n=P+1, W=P+1, NF1=W*W*n, NF2=W*n*n;
    RealType f1[NF1], f2[NF2];
    for (int p=0;p<=P;++p) for (int q=0;q<=P-p;++q) for (int k=0;k<n;++k) { RealType v=0;
        for (int r=0;r<=P-p-q;++r) v += u[((size_t)p*W+q)*W+r]*Cf[(((size_t)p*W+q)*W+r)*n+k]; f1[((size_t)p*W+q)*n+k]=v; }
    for (int p=0;p<=P;++p) for (int j=0;j<n;++j) for (int k=0;k<n;++k) { RealType v=0;
        for (int q=0;q<=P-p;++q) v += Bf[((size_t)p*W+q)*n+j]*f1[((size_t)p*W+q)*n+k]; f2[((size_t)p*n+j)*n+k]=v; }
    for (int i=0;i<n;++i) for (int j=0;j<n;++j) for (int k=0;k<n;++k) { RealType v=0;
        for (int p=0;p<=P;++p) v += Af[(size_t)p*n+i]*f2[((size_t)p*n+j)*n+k]; grad[((size_t)i*n+j)*n+k][comp]=v; }
}

template <typename RealType, int P>
__device__ inline void ho_tet_bwd_sweep(const RealType flux[][3],
        const RealType* __restrict__ Af, const RealType* __restrict__ Bf, const RealType* __restrict__ Cf,
        int comp, RealType* __restrict__ y) {
    constexpr int n=P+1, W=P+1, NH1=W*n*n, NH2=W*W*n;
    RealType h1[NH1], h2[NH2];
    for (int p=0;p<=P;++p) for (int j=0;j<n;++j) for (int k=0;k<n;++k) { RealType v=0;
        for (int i=0;i<n;++i) v += Af[(size_t)p*n+i]*flux[((size_t)i*n+j)*n+k][comp]; h1[((size_t)p*n+j)*n+k]=v; }
    for (int p=0;p<=P;++p) for (int q=0;q<=P-p;++q) for (int k=0;k<n;++k) { RealType v=0;
        for (int j=0;j<n;++j) v += Bf[((size_t)p*W+q)*n+j]*h1[((size_t)p*n+j)*n+k]; h2[((size_t)p*W+q)*n+k]=v; }
    for (int p=0;p<=P;++p) for (int q=0;q<=P-p;++q) for (int r=0;r<=P-p-q;++r) { RealType v=0;
        for (int k=0;k<n;++k) v += Cf[(((size_t)p*W+q)*W+r)*n+k]*h2[((size_t)p*W+q)*n+k]; y[((size_t)p*W+q)*W+r]+=v; }
}

// d_u,d_y: [numElements * W^3]; factor tables as built by buildHoTetOps; d_G: [numElements * nq * 6].
template <typename RealType, int P>
__global__ void ho_laplacian_tet_apply_scalar(
        const RealType* __restrict__ d_u, RealType* __restrict__ d_y,
        const RealType* __restrict__ d_A, const RealType* __restrict__ d_Ad,
        const RealType* __restrict__ d_B, const RealType* __restrict__ d_Bd,
        const RealType* __restrict__ d_C, const RealType* __restrict__ d_Cd,
        const RealType* __restrict__ d_G, size_t numElements) {
    constexpr int n=P+1, W=P+1, nq=n*n*n, W3=W*W*W;
    const size_t e = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    const RealType* u = d_u + e*W3;
    RealType* y = d_y + e*W3;
    const RealType* G = d_G + e*(size_t)nq*6;

    RealType grad[nq][3], flux[nq][3];
    ho_tet_fwd_sweep<RealType,P>(u, d_Ad, d_B,  d_C,  0, grad);
    ho_tet_fwd_sweep<RealType,P>(u, d_A,  d_Bd, d_C,  1, grad);
    ho_tet_fwd_sweep<RealType,P>(u, d_A,  d_B,  d_Cd, 2, grad);
    for (int x=0;x<nq;++x) { const RealType* g=grad[x]; const RealType* Gp=&G[(size_t)x*6];
        flux[x][0]=Gp[0]*g[0]+Gp[1]*g[1]+Gp[2]*g[2];
        flux[x][1]=Gp[1]*g[0]+Gp[3]*g[1]+Gp[4]*g[2];
        flux[x][2]=Gp[2]*g[0]+Gp[4]*g[1]+Gp[5]*g[2]; }
    for (int I=0;I<W3;++I) y[I]=0;
    ho_tet_bwd_sweep<RealType,P>(flux, d_Ad, d_B,  d_C,  0, y);
    ho_tet_bwd_sweep<RealType,P>(flux, d_A,  d_Bd, d_C,  1, y);
    ho_tet_bwd_sweep<RealType,P>(flux, d_A,  d_B,  d_Cd, 2, y);
}

// CG-style nodal apply over a C0 mesh: gather nodal DOF via d_elemDof, transform to modal (Vinv),
// sum-fac apply, transform back (Vinv^T), atomic scatter. NMODES = (P+1)(P+2)(P+3)/6.
// d_modes: [NMODES*3] the (p,q,r) list in HoTetOps order.
template <typename RealType, int P>
__global__ void ho_laplacian_tet_apply_cg(
        const RealType* __restrict__ d_u_dof, RealType* __restrict__ d_y_dof,
        const int* __restrict__ d_elemDof,
        const RealType* __restrict__ d_Vinv, const int* __restrict__ d_modes,
        const RealType* __restrict__ d_A, const RealType* __restrict__ d_Ad,
        const RealType* __restrict__ d_B, const RealType* __restrict__ d_Bd,
        const RealType* __restrict__ d_C, const RealType* __restrict__ d_Cd,
        const RealType* __restrict__ d_G, size_t numElements) {
    constexpr int n = P+1, W = P+1, nq = n*n*n, W3 = W*W*W;
    constexpr int Np = (P+1)*(P+2)*(P+3)/6;
    const size_t e = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (e >= numElements) return;
    const int* dof = d_elemDof + e*Np;
    const RealType* G = d_G + e*(size_t)nq*6;

    RealType ue[Np], uh[Np], um[W3], ym[W3];
    for (int i = 0; i < Np; ++i) ue[i] = d_u_dof[dof[i]];
    for (int m = 0; m < Np; ++m) { RealType v = 0;
        for (int i = 0; i < Np; ++i) v += d_Vinv[(size_t)m*Np+i]*ue[i]; uh[m] = v; }
    for (int I = 0; I < W3; ++I) um[I] = 0;
    for (int m = 0; m < Np; ++m)
        um[((size_t)d_modes[m*3]*W + d_modes[m*3+1])*W + d_modes[m*3+2]] = uh[m];

    RealType grad[nq][3], flux[nq][3];
    ho_tet_fwd_sweep<RealType,P>(um, d_Ad, d_B,  d_C,  0, grad);
    ho_tet_fwd_sweep<RealType,P>(um, d_A,  d_Bd, d_C,  1, grad);
    ho_tet_fwd_sweep<RealType,P>(um, d_A,  d_B,  d_Cd, 2, grad);
    for (int x = 0; x < nq; ++x) { const RealType* g = grad[x]; const RealType* Gp = &G[(size_t)x*6];
        flux[x][0] = Gp[0]*g[0]+Gp[1]*g[1]+Gp[2]*g[2];
        flux[x][1] = Gp[1]*g[0]+Gp[3]*g[1]+Gp[4]*g[2];
        flux[x][2] = Gp[2]*g[0]+Gp[4]*g[1]+Gp[5]*g[2]; }
    for (int I = 0; I < W3; ++I) ym[I] = 0;
    ho_tet_bwd_sweep<RealType,P>(flux, d_Ad, d_B,  d_C,  0, ym);
    ho_tet_bwd_sweep<RealType,P>(flux, d_A,  d_Bd, d_C,  1, ym);
    ho_tet_bwd_sweep<RealType,P>(flux, d_A,  d_B,  d_Cd, 2, ym);

    for (int m = 0; m < Np; ++m)
        uh[m] = ym[((size_t)d_modes[m*3]*W + d_modes[m*3+1])*W + d_modes[m*3+2]];
    for (int i = 0; i < Np; ++i) { RealType v = 0;
        for (int m = 0; m < Np; ++m) v += d_Vinv[(size_t)m*Np+i]*uh[m];   // Vinv^T
        atomicAdd(&d_y_dof[dof[i]], v); }
}

#endif  // USE_CUDA

}  // namespace fem
}  // namespace mars
#endif
