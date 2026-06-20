// Paper Table 1 replica (host): manufactured-solution convergence of the HO-CVFEM
// diffusion operator. Solve  -div(grad u) = f  on the unit cube with
//   u_exact = sin(pi x) sin(pi y) sin(pi z)   (=> homogeneous Dirichlet, u=0 on dOmega)
//   f = 3 pi^2 u_exact
// using the consistent CVFEM stiffness (Knaus Alg 2) and consistent control-volume
// load b = M f (Alg 1). Small dense solve (the operator is non-symmetric at p>=2,
// so a direct LU, not CG). Report L2 error + convergence rate vs DOFs per order.
// The paper reports ~p+1 order; this checks our operator reproduces it.

#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_mass.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace mars::fem;
static const double PI = std::acos(-1.0);

// dense LU solve with partial pivoting (Ax=b, A row-major n x n; A,b clobbered)
static void solveDense(std::vector<double>& A, std::vector<double>& b, int n, std::vector<double>& x)
{
    for (int k = 0; k < n; ++k) {
        int mx = k; for (int i = k+1; i < n; ++i) if (std::abs(A[i*n+k]) > std::abs(A[mx*n+k])) mx = i;
        if (mx != k) { for (int j = k; j < n; ++j) std::swap(A[k*n+j], A[mx*n+j]); std::swap(b[k], b[mx]); }
        double akk = A[k*n+k];
        for (int i = k+1; i < n; ++i) {
            double fac = A[i*n+k] / akk;
            if (fac == 0.0) continue;
            for (int j = k; j < n; ++j) A[i*n+j] -= fac * A[k*n+j];
            b[i] -= fac * b[k];
        }
    }
    x.assign(n, 0.0);
    for (int i = n-1; i >= 0; --i) { double s = b[i]; for (int j = i+1; j < n; ++j) s -= A[i*n+j]*x[j]; x[i] = s / A[i*n+i]; }
}

// true global L2 error ||u_h - u_exact||_L2 / ||u_exact||_L2 via Gauss quadrature
// of the FE interpolant error (the paper's metric, gives the optimal p+1 rate).
template<int P>
static double globalL2(const HoCvfemOperators& op, const HODofHandler& dh,
                       const std::vector<std::array<int,3>>& ijk, int E,
                       const std::vector<double>& uh)
{
    const int n=P+1, N3=n*n*n, nq=P+2;
    const double h=1.0/E, detJ=std::pow(h*0.5,3.0);
    std::vector<double> gx, gw; gaussLegendre(nq, gx, gw);
    std::vector<std::vector<double>> L(nq);   // L[q][i] = l_i(gx[q]) on GLL nodes
    for (int q=0;q<nq;++q){ std::vector<double> val,der; lagrangeValDeriv(op.zeta, gx[q], val, der); L[q]=val; }
    double err2=0, nrm2=0;
    for (size_t e=0;e<ijk.size();++e) {
        std::vector<double> ul(N3);
        for (int l=0;l<N3;++l) ul[l]=uh[dh.elemDof[e*N3+l]];
        for (int qi=0;qi<nq;++qi) for (int qj=0;qj<nq;++qj) for (int qk=0;qk<nq;++qk) {
            double w=gw[qi]*gw[qj]*gw[qk], uhq=0;
            for (int i=0;i<n;++i) for (int j=0;j<n;++j) for (int k=0;k<n;++k)
                uhq += L[qi][i]*L[qj][j]*L[qk][k]*ul[(i*n+j)*n+k];
            double px=(ijk[e][0]+0.5*(gx[qi]+1))*h, py=(ijk[e][1]+0.5*(gx[qj]+1))*h, pz=(ijk[e][2]+0.5*(gx[qk]+1))*h;
            double ue=std::sin(PI*px)*std::sin(PI*py)*std::sin(PI*pz), d=uhq-ue;
            err2 += d*d*w*detJ; nrm2 += ue*ue*w*detJ;
        }
    }
    return std::sqrt(err2/(nrm2>0?nrm2:1.0));
}

template<int P>
static double runConv(int E, long& outDof)
{
    auto op = buildHoCvfemOperators(P);
    const int n = P+1, N3 = n*n*n;
    const double h = 1.0 / E;
    auto cg = [&](int x,int y,int z){ return (x*(E+1)+y)*(E+1)+z; };
    std::vector<std::array<int,8>> ec; std::vector<std::array<int,3>> ijk;
    for (int ex=0;ex<E;++ex) for (int ey=0;ey<E;++ey) for (int ez=0;ez<E;++ez){
        ec.push_back({cg(ex,ey,ez),cg(ex+1,ey,ez),cg(ex+1,ey+1,ez),cg(ex,ey+1,ez),
                      cg(ex,ey,ez+1),cg(ex+1,ey,ez+1),cg(ex+1,ey+1,ez+1),cg(ex,ey+1,ez+1)});
        ijk.push_back({ex,ey,ez});
    }
    HODofHandler dh; dh.build(ec, long(E+1)*(E+1)*(E+1), P);
    const long nDof = dh.numDof; const size_t nEl = ec.size();

    // DOF coords (GLL positions), exact solution, source, boundary flag
    std::vector<double> uex(nDof,0.0), f(nDof,0.0); std::vector<uint8_t> bdry(nDof,0);
    auto coord = [&](int e_ijk, int loc){ return (e_ijk + 0.5*(op.zeta[loc]+1.0))*h; };
    for (size_t e=0;e<nEl;++e)
        for (int i=0;i<n;++i) for (int j=0;j<n;++j) for (int k=0;k<n;++k) {
            long dof = dh.elemDof[e*N3 + i*n*n+j*n+k];
            double px=coord(ijk[e][0],i), py=coord(ijk[e][1],j), pz=coord(ijk[e][2],k);
            double u = std::sin(PI*px)*std::sin(PI*py)*std::sin(PI*pz);
            uex[dof]=u; f[dof]=3.0*PI*PI*u;
            if (px<1e-12||px>1-1e-12||py<1e-12||py>1-1e-12||pz<1e-12||pz>1-1e-12) bdry[dof]=1;
        }

    // assemble dense stiffness A, consistent load b = M f, lumped volume V = M 1
    std::vector<double> A((size_t)nDof*nDof, 0.0), b(nDof,0.0), V(nDof,0.0);
    const double detJ = std::pow(h*0.5, 3.0);
    const int (&S)[8][3] = hexCornerRef();
    for (size_t e=0;e<nEl;++e) {
        double corners[8][3];
        for (int c=0;c<8;++c) {
            corners[c][0]=(ijk[e][0]+(S[c][0]+1)*0.5)*h;
            corners[c][1]=(ijk[e][1]+(S[c][1]+1)*0.5)*h;
            corners[c][2]=(ijk[e][2]+(S[c][2]+1)*0.5)*h;
        }
        auto G = computeElementMetric(op, corners);
        // K_e columns -> dense A
        std::vector<double> el(N3,0.0), y;
        for (int J=0;J<N3;++J) {
            el.assign(N3,0.0); el[J]=1.0;
            applyHoCvfemElement(op, G, el, y);
            long gJ = dh.elemDof[e*N3+J];
            for (int I=0;I<N3;++I) A[(size_t)dh.elemDof[e*N3+I]*nDof + gJ] += y[I];
        }
        // load + lumped volume
        std::vector<double> fe(N3), ones(N3,1.0), ym;
        for (int l=0;l<N3;++l) fe[l]=f[dh.elemDof[e*N3+l]];
        applyMassCube(op, detJ, fe, ym);   for (int l=0;l<N3;++l) b[dh.elemDof[e*N3+l]] += ym[l];
        applyMassCube(op, detJ, ones, ym); for (int l=0;l<N3;++l) V[dh.elemDof[e*N3+l]] += ym[l];
    }

    // homogeneous Dirichlet: boundary rows -> identity, rhs -> u_exact (=0)
    for (long d=0; d<nDof; ++d) if (bdry[d]) {
        for (long j=0;j<nDof;++j) A[(size_t)d*nDof+j]=0.0;
        A[(size_t)d*nDof+d]=1.0; b[d]=uex[d];
    }

    std::vector<double> u;
    solveDense(A, b, (int)nDof, u);

    outDof = nDof;
    return globalL2<P>(op, dh, ijk, E, u);   // true integrated L2 (paper metric)
}

template<int P>
static void sweep(const std::vector<int>& Es)
{
    printf("\n p=%d :  #DOFs   |   L2 error   | rate\n", P);
    double prevErr=0; long prevDof=0;
    for (size_t i=0;i<Es.size();++i) {
        long dof; double err = runConv<P>(Es[i], dof);
        double rate = (i==0) ? 0.0 :
            std::log(prevErr/err) / std::log(std::cbrt((double)dof/(double)prevDof));
        if (i==0) printf("       %8ld | %.4e |   -\n", dof, err);
        else      printf("       %8ld | %.4e | %.2f\n", dof, err, rate);
        prevErr=err; prevDof=dof;
    }
}

int main()
{
    printf("HO-CVFEM manufactured-solution convergence (paper Table 1 replica)\n");
    printf("u = sin(pi x)sin(pi y)sin(pi z), homogeneous Dirichlet; expect ~p+1 rate\n");
    sweep<1>({4,6,8,12});
    sweep<2>({2,3,4,6});
    sweep<3>({2,3,4});
    sweep<4>({1,2,3});
    return 0;
}
