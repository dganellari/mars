// Stage C consistency patch test (host-only): assemble the high-order CVFEM
// diffusion operator on a structured cube via the topological HO DOF map, and
// verify it is a CONSISTENT Laplacian:
//   A * 1            = 0   (everywhere)
//   A * (linear x)   = 0   at INTERIOR DOFs
// This is the decisive check -- it holds for a consistent operator whether or
// not it is symmetric, so it tells us if the p>=2 asymmetry is legitimate CVFEM
// (Petrov-Galerkin) or an Alg-2 bug. coeff=1 per element (consistency is
// scaling-invariant). Uses the same elemDof convention as the apply
// (l = i*n*n + j*n + k).

#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_apply.hpp"

#include <array>
#include <cstdio>
#include <cmath>
#include <vector>

using namespace mars::fem;

static bool patch(int P, int E)
{
    auto op = buildHoCvfemOperators(P);
    const int n = P + 1, N3 = n * n * n;

    // structured cube corner connectivity (hex order matching the handler)
    auto cg = [&](int x, int y, int z) { return (x*(E+1)+y)*(E+1)+z; };
    std::vector<std::array<int,8>> ec; std::vector<std::array<int,3>> ijk;
    for (int ex=0; ex<E; ++ex) for (int ey=0; ey<E; ++ey) for (int ez=0; ez<E; ++ez) {
        ec.push_back({cg(ex,ey,ez),cg(ex+1,ey,ez),cg(ex+1,ey+1,ez),cg(ex,ey+1,ez),
                      cg(ex,ey,ez+1),cg(ex+1,ey,ez+1),cg(ex+1,ey+1,ez+1),cg(ex,ey+1,ez+1)});
        ijk.push_back({ex,ey,ez});
    }
    HODofHandler dh; dh.build(ec, long(E+1)*(E+1)*(E+1), P);
    long nDof = dh.numDof; size_t nEl = ec.size();
    const double h = 1.0 / E;

    // physical x-coord + boundary flag per DOF (GLL node positions)
    std::vector<double> dofX(nDof, 0.0); std::vector<uint8_t> bdry(nDof, 0);
    auto coord = [&](int e_ijk, int loc){ return (e_ijk + 0.5*(op.zeta[loc]+1.0))*h; };
    for (size_t e=0; e<nEl; ++e)
        for (int i=0;i<n;++i) for (int j=0;j<n;++j) for (int k=0;k<n;++k) {
            int dof = dh.elemDof[e*N3 + i*n*n+j*n+k];
            double px=coord(ijk[e][0],i), py=coord(ijk[e][1],j), pz=coord(ijk[e][2],k);
            dofX[dof]=px;
            bdry[dof] = (px<1e-12||px>1-1e-12||py<1e-12||py>1-1e-12||pz<1e-12||pz>1-1e-12)?1:0;
        }

    auto apply = [&](const std::vector<double>& u){
        std::vector<double> y(nDof, 0.0), ul(N3), yl;
        for (size_t e=0;e<nEl;++e){
            for (int l=0;l<N3;++l) ul[l]=u[dh.elemDof[e*N3+l]];
            applyHoCvfemElementCube(op, 1.0, ul, yl);
            for (int l=0;l<N3;++l) y[dh.elemDof[e*N3+l]] += yl[l];
        }
        return y;
    };

    std::vector<double> ones(nDof,1.0);
    auto y1 = apply(ones);
    double nullMax=0; for (double v:y1) nullMax=std::max(nullMax,std::abs(v));

    auto yx = apply(dofX);
    double interMax=0; for (long p=0;p<nDof;++p) if(!bdry[p]) interMax=std::max(interMax,std::abs(yx[p]));

    bool ok = (nullMax<1e-9) && (interMax<1e-9);
    printf("p=%d E=%d | nDof=%ld | A*1=%.2e | A*linear interior=%.2e | %s\n",
           P, E, nDof, nullMax, interMax, ok?"PASS":"FAIL");
    return ok;
}

int main()
{
    bool ok = true;
    ok &= patch(1, 4);
    ok &= patch(2, 3);
    ok &= patch(4, 2);
    printf("Stage C consistency (HO-CVFEM patch test): %s\n", ok?"PASS":"FAIL");
    return ok?0:1;
}
