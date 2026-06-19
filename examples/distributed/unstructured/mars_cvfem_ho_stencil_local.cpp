// Local p=1 stencil gate: dump the high-order CVFEM (Knaus Alg 2) element
// stiffness K_e (8x8) on a unit cube and compare it to the analytic trilinear-hex
// (Q1) Galerkin stiffness. Purpose: understand the measured 64% structural
// difference vs MARS's SCS-flux CvfemHexAssembler -- is OUR operator the textbook
// stencil (=> the other assembler is the outlier), or a distinct CVFEM stencil?
// Host-only; same MARS hex corner ordering as the kernels.

#include "backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp"

#include <array>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace mars::fem;

// analytic Q1 Galerkin stiffness on the unit cube [0,1]^3, MARS corner ordering
// (corner c at reference signs hexCornerRef[c]). K_ij = int_[0,1]^3 grad Ni . grad Nj.
static std::array<double,64> q1Galerkin()
{
    const auto& S = hexCornerRef();
    // corner physical coords on [0,1]^3
    double X[8][3];
    for (int c=0;c<8;++c) for (int d=0;d<3;++d) X[c][d] = (S[c][d]+1)*0.5;
    // trilinear shape grad at physical point p (unit cube): Ni = prod (p_d or 1-p_d)
    auto gradN = [&](int c, const double p[3], double g[3]) {
        double f[3], df[3];
        for (int d=0;d<3;++d) { if (X[c][d] > 0.5) { f[d]=p[d];     df[d]=1.0; }
                                else                { f[d]=1.0-p[d]; df[d]=-1.0; } }
        g[0]=df[0]*f[1]*f[2]; g[1]=f[0]*df[1]*f[2]; g[2]=f[0]*f[1]*df[2];
    };
    // 2-pt Gauss on [0,1] (exact for the quadratic integrand)
    const double gp[2]={0.5-0.5/std::sqrt(3.0), 0.5+0.5/std::sqrt(3.0)}, gw=0.5;
    std::array<double,64> K{}; K.fill(0.0);
    for (int a=0;a<2;++a) for (int b=0;b<2;++b) for (int cc=0;cc<2;++cc) {
        double p[3]={gp[a],gp[b],gp[cc]}, w=gw*gw*gw;
        double g[8][3]; for (int c=0;c<8;++c) gradN(c,p,g[c]);
        for (int i=0;i<8;++i) for (int j=0;j<8;++j)
            K[i*8+j] += w*(g[i][0]*g[j][0]+g[i][1]*g[j][1]+g[i][2]*g[j][2]);
    }
    return K;
}

int main()
{
    auto op = buildHoCvfemOperators(1);
    const auto& S = hexCornerRef();

    // unit-cube corners, MARS order
    double corners[8][3];
    for (int c=0;c<8;++c) for (int d=0;d<3;++d) corners[c][d] = (S[c][d]+1)*0.5;
    auto G = computeElementMetric(op, corners);

    // OUR K_e: apply to the 8 unit vectors. Local node l = i*4+j*2+k; corner c maps
    // to l via the same (i,j,k)=(S+1)/2 rule, so build K in CORNER order to align
    // with the Q1 matrix.
    int cornerToLocal[8];
    for (int c=0;c<8;++c) {
        int i=(S[c][0]+1)/2, j=(S[c][1]+1)/2, k=(S[c][2]+1)/2;
        cornerToLocal[c] = i*4 + j*2 + k;
    }
    std::array<double,64> Kours{}; Kours.fill(0.0);
    for (int cj=0;cj<8;++cj) {
        std::vector<double> e(8,0.0), y; e[cornerToLocal[cj]] = 1.0;
        applyHoCvfemElement(op, G, e, y);
        for (int ci=0;ci<8;++ci) Kours[ci*8+cj] = y[cornerToLocal[ci]];
    }

    auto Kq1 = q1Galerkin();

    // row sums (nullspace), symmetry, and the comparison
    auto rowsum = [&](const std::array<double,64>& K, int i){ double s=0; for(int j=0;j<8;++j) s+=K[i*8+j]; return s; };
    double ourNull=0, q1Null=0, ourSym=0;
    for (int i=0;i<8;++i){ ourNull=std::max(ourNull,std::abs(rowsum(Kours,i))); q1Null=std::max(q1Null,std::abs(rowsum(Kq1,i))); }
    for (int i=0;i<8;++i) for (int j=0;j<8;++j) ourSym=std::max(ourSym,std::abs(Kours[i*8+j]-Kours[j*8+i]));

    auto dump = [&](const char* name, const std::array<double,64>& K){
        printf("\n%s (8x8):\n", name);
        for (int i=0;i<8;++i){ printf("  "); for (int j=0;j<8;++j) printf("%8.4f ", K[i*8+j]); printf("\n"); }
    };
    dump("OUR Knaus CVFEM K_e", Kours);
    dump("analytic Q1 Galerkin K_e", Kq1);

    // best constant scale aligning ours to Q1, then residual
    double num=0, den=0;
    for (int e=0;e<64;++e){ num+=Kours[e]*Kq1[e]; den+=Kours[e]*Kours[e]; }
    double scale = den>0 ? num/den : 0.0;
    double maxAbsQ1=0, maxDiff=0;
    for (int e=0;e<64;++e){ maxAbsQ1=std::max(maxAbsQ1,std::abs(Kq1[e])); maxDiff=std::max(maxDiff,std::abs(scale*Kours[e]-Kq1[e])); }

    printf("\n-- summary --\n");
    printf("  OUR  : rowsum(null)=%.2e  sym=%.2e\n", ourNull, ourSym);
    printf("  Q1   : rowsum(null)=%.2e\n", q1Null);
    printf("  best scale (ours->Q1) = %.6f\n", scale);
    printf("  max|scale*ours - Q1| / max|Q1| = %.3e  [%s]\n",
           maxDiff/(maxAbsQ1>0?maxAbsQ1:1.0),
           maxDiff/(maxAbsQ1>0?maxAbsQ1:1.0) < 1e-10 ? "OUR == Q1 (up to scale)" : "OUR != Q1 (distinct stencil)");
    return 0;
}
