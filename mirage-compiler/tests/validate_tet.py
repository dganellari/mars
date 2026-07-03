#!/usr/bin/env python3
"""Local proof of the Mirage tet collapsed sum-factorization apply.

Emits mars::mirage::tet_laplacian_apply_host (the generated B^T G_c B sweep) and
compares it against the DENSE per-mode PKD oracle on the same collapsed tet, for
p=2..5. Reproduces rung 9 (see mirage-compiler/TET_SUMFAC_HANDOFF.md):
  generated == dense (~1e-14),  A.1=0 exactly,  symmetric (~1e-14).

Requires the reference headers (jacobi/gauss_jacobi/basis) from the sibling
tet_ho_ref reference; point MIRAGE_TET_ORACLE at that dir (default: this session's
scratchpad copy). Run from mirage-compiler/:  python3 tests/validate_tet.py
"""

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GEN = os.path.join(ROOT, "generated")
ORACLE = os.environ.get(
    "MIRAGE_TET_ORACLE",
    "/private/tmp/claude-501/-Users-gandanie-scratch-santis-mars/"
    "8ed75b86-c40b-4adf-9eeb-9a07dee15bc2/scratchpad/tet_ho_ref")

sys.path.insert(0, ROOT)
from mirage.backends import tet_galerkin          # noqa: E402

# Driver: build the PKD tables + collapse metric + a DENSE oracle (all from the
# reference library), then diff the GENERATED apply against the dense operator.
DRIVER = r"""
#include "gauss_jacobi.hpp"
#include "jacobi.hpp"
#include "basis.hpp"                 // mat_inverse
#include "tet_laplacian_apply_host.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
using namespace tetref;

static double det3(const double* M){ return M[0]*(M[4]*M[8]-M[5]*M[7])-M[1]*(M[3]*M[8]-M[5]*M[6])+M[2]*(M[3]*M[7]-M[4]*M[6]); }

int main() {
    const double V[4][3] = {{0.1,0.2,0.3},{2.1,0.3,0.4},{0.6,1.7,0.1},{0.4,0.5,1.5}};
    double Bm[9]; for(int d=0;d<3;++d){ Bm[d*3+0]=V[1][d]-V[0][d]; Bm[d*3+1]=V[2][d]-V[0][d]; Bm[d*3+2]=V[3][d]-V[0][d]; }
    const double detB = std::fabs(det3(Bm));
    std::vector<double> Minv, Bv(Bm,Bm+9); mat_inverse(Bv,3,Minv);
    double MMt[9]; for(int a=0;a<3;++a) for(int b=0;b<3;++b){ double s=0; for(int c=0;c<3;++c) s+=Minv[a*3+c]*Minv[b*3+c]; MMt[a*3+b]=s; }

    std::printf(" p | modes | max|gen - dense|    |  A.1=0   | symmetry\n");
    std::printf("---+-------+---------------------+----------+---------\n");
    int fails = 0;
    for (int D = 2; D <= 5; ++D) {
        const int n = D+1, W_ = D+1, NQ = n*n*n;
        Quad ga=gauss_jacobi(n,0.0,0.0), gb=gauss_jacobi(n,1.0,0.0), gc=gauss_jacobi(n,2.0,0.0);
        std::vector<std::array<int,3>> modes;
        for(int p=0;p<=D;++p)for(int q=0;q<=D-p;++q)for(int r=0;r<=D-p-q;++r) modes.push_back({p,q,r});
        const int Mn = (int)modes.size();

        std::vector<double> A((size_t)W_*n), Ad((size_t)W_*n), B((size_t)W_*W_*n), Bd((size_t)W_*W_*n),
                            C((size_t)W_*W_*W_*n), Cd((size_t)W_*W_*W_*n);
        for(int p=0;p<=D;++p) for(int i=0;i<n;++i){ A[(size_t)p*n+i]=jacobi(p,0.0,0.0,ga.x[i]); Ad[(size_t)p*n+i]=2.0*djacobi(p,0.0,0.0,ga.x[i]); }
        for(int p=0;p<=D;++p) for(int q=0;q<=D-p;++q) for(int j=0;j<n;++j){ double x=gb.x[j],w=0.5*(1-x);
            double J=jacobi(q,2.0*p+1,0.0,x), dJ=djacobi(q,2.0*p+1,0.0,x);
            B[((size_t)p*W_+q)*n+j]=std::pow(w,p)*J;
            Bd[((size_t)p*W_+q)*n+j]=((p>0)?-(double)p*std::pow(w,p-1)*J:0.0)+std::pow(w,p)*2.0*dJ; }
        for(int p=0;p<=D;++p) for(int q=0;q<=D-p;++q) for(int r=0;r<=D-p-q;++r) for(int k=0;k<n;++k){ int e=p+q; double x=gc.x[k],w=0.5*(1-x);
            double J=jacobi(r,2.0*e+2,0.0,x), dJ=djacobi(r,2.0*e+2,0.0,x);
            C[(((size_t)p*W_+q)*W_+r)*n+k]=std::pow(w,e)*J;
            Cd[(((size_t)p*W_+q)*W_+r)*n+k]=((e>0)?-(double)e*std::pow(w,e-1)*J:0.0)+std::pow(w,e)*2.0*dJ; }

        std::vector<std::array<double,9>> Gc(NQ);
        for(int i=0;i<n;++i) for(int j=0;j<n;++j) for(int k=0;k<n;++k){
            const double a=0.5*(1+ga.x[i]), b=0.5*(1+gb.x[j]), c=0.5*(1+gc.x[k]);
            const double wa=0.5*ga.w[i], wb=0.25*gb.w[j], wc=0.125*gc.w[k];
            const double Wq=wa*wb*wc*detB;
            const double r=a*(1-b)*(1-c), s=b*(1-c);
            const double D1=1.0-s-c, D2=1.0-c;
            const double Jc[9]={ 1.0/D1, r/(D1*D1), r/(D1*D1),   0.0, 1.0/D2, s/(D2*D2),   0.0, 0.0, 1.0 };
            double JM[9]; for(int p=0;p<3;++p)for(int t=0;t<3;++t){ double v=0; for(int u=0;u<3;++u) v+=Jc[p*3+u]*MMt[u*3+t]; JM[p*3+t]=v; }
            double G[9];  for(int p=0;p<3;++p)for(int t=0;t<3;++t){ double v=0; for(int u=0;u<3;++u) v+=JM[p*3+u]*Jc[t*3+u]; G[p*3+t]=Wq*v; }
            for(int m=0;m<9;++m) Gc[((size_t)i*n+j)*n+k][m]=G[m];
        }

        // DENSE per-mode oracle (grad_x psi_m) -- the ground truth.
        std::vector<double> GX((size_t)Mn*NQ*3, 0.0);
        for(int m=0;m<Mn;++m){ int p=modes[m][0],q=modes[m][1],r=modes[m][2];
            for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k){ size_t x=((size_t)i*n+j)*n+k;
                double gc0=Ad[(size_t)p*n+i]*B[((size_t)p*W_+q)*n+j]*C[(((size_t)p*W_+q)*W_+r)*n+k];
                double gc1=A[(size_t)p*n+i]*Bd[((size_t)p*W_+q)*n+j]*C[(((size_t)p*W_+q)*W_+r)*n+k];
                double gc2=A[(size_t)p*n+i]*B[((size_t)p*W_+q)*n+j]*Cd[(((size_t)p*W_+q)*W_+r)*n+k];
                const double a=0.5*(1+ga.x[i]), b=0.5*(1+gb.x[j]), c=0.5*(1+gc.x[k]);
                const double rr=a*(1-b)*(1-c), ss=b*(1-c); const double D1=1.0-ss-c,D2=1.0-c;
                const double Jc[9]={1.0/D1,rr/(D1*D1),rr/(D1*D1), 0.0,1.0/D2,ss/(D2*D2), 0.0,0.0,1.0};
                double gcv[3]={gc0,gc1,gc2}, gref[3], gx[3];
                for(int d=0;d<3;++d){ double v=0; for(int u=0;u<3;++u) v+=Jc[u*3+d]*gcv[u]; gref[d]=v; }
                for(int d=0;d<3;++d){ double v=0; for(int u=0;u<3;++u) v+=Minv[u*3+d]*gref[u]; gx[d]=v; }
                GX[((size_t)m*NQ+x)*3+0]=gx[0]; GX[((size_t)m*NQ+x)*3+1]=gx[1]; GX[((size_t)m*NQ+x)*3+2]=gx[2];
            } }
        std::vector<double> Wq(NQ);
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)for(int k=0;k<n;++k) Wq[((size_t)i*n+j)*n+k]=(0.5*ga.w[i])*(0.25*gb.w[j])*(0.125*gc.w[k])*detB;
        auto Aden=[&](const std::vector<double>& u, std::vector<double>& y){
            std::vector<std::array<double,3>> gxu(NQ,{0,0,0});
            for(int m=0;m<Mn;++m){ int pp=modes[m][0],qq=modes[m][1],rr2=modes[m][2]; double um=u[((size_t)pp*W_+qq)*W_+rr2];
                for(int x=0;x<NQ;++x) for(int d=0;d<3;++d) gxu[x][d]+=um*GX[((size_t)m*NQ+x)*3+d]; }
            y.assign((size_t)W_*W_*W_,0.0);
            for(int m=0;m<Mn;++m){ int pp=modes[m][0],qq=modes[m][1],rr2=modes[m][2]; double s=0;
                for(int x=0;x<NQ;++x){ double d=0; for(int c=0;c<3;++c) d+=GX[((size_t)m*NQ+x)*3+c]*gxu[x][c]; s+=Wq[x]*d; }
                y[((size_t)pp*W_+qq)*W_+rr2]=s; }
        };

        // GENERATED apply from Mirage.
        auto Agen=[&](const std::vector<double>& u, std::vector<double>& y){
            mars::mirage::tet_laplacian_apply_host(D, A, Ad, B, Bd, C, Cd, Gc, u, y);
        };

        std::vector<double> u((size_t)W_*W_*W_,0.0), y1((size_t)W_*W_*W_), y2((size_t)W_*W_*W_);
        for(int m=0;m<Mn;++m){ auto& e=modes[m]; u[((size_t)e[0]*W_+e[1])*W_+e[2]]=std::sin(0.9*m+0.4); }
        Agen(u,y1); Aden(u,y2);
        double err=0; for(int m=0;m<Mn;++m){ auto& e=modes[m]; size_t x=((size_t)e[0]*W_+e[1])*W_+e[2]; err=std::max(err,std::fabs(y1[x]-y2[x])); }

        std::vector<double> one((size_t)W_*W_*W_,0.0), yo((size_t)W_*W_*W_); one[0]=1.0; Agen(one,yo);
        double a1=0; for(int m=0;m<Mn;++m){ auto& e=modes[m]; a1=std::max(a1,std::fabs(yo[((size_t)e[0]*W_+e[1])*W_+e[2]])); }

        std::vector<double> v((size_t)W_*W_*W_,0.0), Av((size_t)W_*W_*W_), Au((size_t)W_*W_*W_);
        for(int m=0;m<Mn;++m){ auto& e=modes[m]; v[((size_t)e[0]*W_+e[1])*W_+e[2]]=std::cos(0.5*m+1.1); }
        Agen(v,Av); Agen(u,Au); double uAv=0,vAu=0;
        for(int m=0;m<Mn;++m){ auto& e=modes[m]; size_t x=((size_t)e[0]*W_+e[1])*W_+e[2]; uAv+=u[x]*Av[x]; vAu+=v[x]*Au[x]; }
        double sym=std::fabs(uAv-vAu);
        bool ok = (err<1e-11)&&(a1<1e-12)&&(sym<1e-11);
        std::printf(" %d | %5d |      %.2e      | %.1e  | %.1e  %s\n", D, Mn, err, a1, sym, ok?"OK":"FAIL");
        if(!ok) ++fails;
    }
    std::printf(fails? "\n%d FAILED\n":"\nALL PASS (generated tet apply == dense oracle)\n", fails);
    return fails?1:0;
}
"""


def main():
    if not os.path.isdir(ORACLE):
        print(f"ERROR: tet reference not found at {ORACLE}\n"
              f"Set MIRAGE_TET_ORACLE to the tet_ho_ref dir.", file=sys.stderr)
        sys.exit(2)
    os.makedirs(GEN, exist_ok=True)
    with open(os.path.join(GEN, "tet_laplacian_apply_host.hpp"), "w") as f:
        f.write(tet_galerkin.emit("tet_laplacian"))

    drv = os.path.join(GEN, "validate_tet_driver.cpp")
    with open(drv, "w") as f:
        f.write(DRIVER)
    exe = os.path.join(GEN, "validate_tet_driver")
    cc = os.environ.get("CXX", "c++")
    cmd = [cc, "-std=c++17", "-O2", f"-I{ORACLE}", f"-I{GEN}", drv, "-o", exe]
    print("compile:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("run:", exe)
    sys.exit(subprocess.run([exe]).returncode)


if __name__ == "__main__":
    main()
