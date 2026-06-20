// Stage B validation (host-only): the general per-element metric + non-orthogonal
// cross terms must keep the operator a CONSISTENT Laplacian on a DISTORTED hex
// mesh (sheared + stretched => non-orthogonal, variable Jacobian). Patch test:
//   A * 1          = 0   everywhere
//   A * (x_phys)   = 0   at interior DOFs   (linear physical field annihilated)
// The affine cube never exercised cross terms or a varying J; this does.

#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp"

#include <array>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace mars::fem;

// smooth distortion of [0,1]^3 -> non-orthogonal (shear) + variable-J (quadratic)
static void distort(double X, double Y, double Z, double& x, double& y, double& z) {
    x = X + 0.15*Y + 0.05*X*(1-X);
    y = Y + 0.10*Z + 0.05*Y*(1-Y);
    z = Z + 0.05*Z*(1-Z);
}

static void trilinearPos(const double C[8][3], const double r[3], double out[3]) {
    const auto& S = hexCornerRef();
    out[0]=out[1]=out[2]=0;
    for (int c=0;c<8;++c) {
        double N = 0.125*(1+S[c][0]*r[0])*(1+S[c][1]*r[1])*(1+S[c][2]*r[2]);
        for (int a=0;a<3;++a) out[a]+=N*C[c][a];
    }
}

static bool patch(int P, int E)
{
    auto op = buildHoCvfemOperators(P);
    const int n=P+1, N3=n*n*n;
    const double h=1.0/E;
    auto cg = [&](int x,int y,int z){ return (x*(E+1)+y)*(E+1)+z; };

    std::vector<std::array<int,8>> ec; std::vector<std::array<int,3>> ijk;
    for (int ex=0;ex<E;++ex) for (int ey=0;ey<E;++ey) for (int ez=0;ez<E;++ez){
        ec.push_back({cg(ex,ey,ez),cg(ex+1,ey,ez),cg(ex+1,ey+1,ez),cg(ex,ey+1,ez),
                      cg(ex,ey,ez+1),cg(ex+1,ey,ez+1),cg(ex+1,ey+1,ez+1),cg(ex,ey+1,ez+1)});
        ijk.push_back({ex,ey,ez});
    }
    HODofHandler dh; dh.build(ec, long(E+1)*(E+1)*(E+1), P);
    long nDof = dh.numDof; size_t nEl = ec.size();

    // base corner offsets (unit ref-cube order) in element-local units
    const auto& S = hexCornerRef();

    std::vector<double> dofX(nDof,0.0); std::vector<uint8_t> bdry(nDof,0);
    std::vector<std::vector<std::array<double,3>>> Gelem(nEl);

    for (size_t e=0;e<nEl;++e) {
        // base + distorted corner coords of this element
        double baseC[8][3], distC[8][3];
        for (int c=0;c<8;++c) {
            double X = (ijk[e][0] + (S[c][0]+1)/2.0)*h;
            double Y = (ijk[e][1] + (S[c][1]+1)/2.0)*h;
            double Z = (ijk[e][2] + (S[c][2]+1)/2.0)*h;
            baseC[c][0]=X; baseC[c][1]=Y; baseC[c][2]=Z;
            distort(X,Y,Z, distC[c][0],distC[c][1],distC[c][2]);
        }
        Gelem[e] = computeElementMetric(op, distC);
        for (int i=0;i<n;++i) for (int j=0;j<n;++j) for (int k=0;k<n;++k) {
            int dof = dh.elemDof[e*N3 + i*n*n+j*n+k];
            double r[3]={op.zeta[i],op.zeta[j],op.zeta[k]}, bp[3], dp[3];
            trilinearPos(baseC, r, bp); trilinearPos(distC, r, dp);
            dofX[dof] = dp[0];
            bdry[dof] = (bp[0]<1e-12||bp[0]>1-1e-12||bp[1]<1e-12||bp[1]>1-1e-12||bp[2]<1e-12||bp[2]>1-1e-12)?1:0;
        }
    }

    auto apply = [&](const std::vector<double>& u){
        std::vector<double> y(nDof,0.0), ul(N3), yl;
        for (size_t e=0;e<nEl;++e){
            for (int l=0;l<N3;++l) ul[l]=u[dh.elemDof[e*N3+l]];
            applyHoCvfemElement(op, Gelem[e], ul, yl);
            for (int l=0;l<N3;++l) y[dh.elemDof[e*N3+l]] += yl[l];
        }
        return y;
    };

    std::vector<double> ones(nDof,1.0);
    auto y1=apply(ones); double nullMax=0; for(double v:y1) nullMax=std::max(nullMax,std::abs(v));
    auto yx=apply(dofX); double interMax=0; for(long p=0;p<nDof;++p) if(!bdry[p]) interMax=std::max(interMax,std::abs(yx[p]));

    bool ok=(nullMax<1e-9)&&(interMax<1e-9);
    printf("p=%d E=%d (distorted) | nDof=%ld | A*1=%.2e | A*linear interior=%.2e | %s\n",
           P,E,nDof,nullMax,interMax,ok?"PASS":"FAIL");
    return ok;
}

int main()
{
    bool ok=true;
    ok &= patch(1,4);
    ok &= patch(2,3);
    ok &= patch(4,2);
    printf("Stage B (HO-CVFEM general metric + cross terms, distorted mesh): %s\n", ok?"PASS":"FAIL");
    return ok?0:1;
}
