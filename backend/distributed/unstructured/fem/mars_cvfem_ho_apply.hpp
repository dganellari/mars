#pragma once

// Stage C (host reference): element diffusion-stiffness apply y = K_e u for a
// tensor-product hex, following Knaus Alg 2. This is the correctness ground
// truth the GPU kernel will be ported against. Diffusion only; here specialized
// to an affine (orthogonal) cube element so the metric is a constant diagonal
// (no cross terms) -- enough to validate the sum-factorization STRUCTURE
// (constant-nullspace, symmetry on the orthogonal cube, consistency). The exact
// metric/area scaling is folded into `coeff`; structural correctness is
// scaling-invariant, and absolute scaling is pinned later vs --kernel=tensor.
//
// Alg 2 per direction `dir` (normal index), at each subcontrol-surface l:
//   deriv = Dtil . u        (normal derivative GLL->SCS face)
//   flux  = coeff * deriv    (diffusion flux; orthogonal cube -> no cross terms)
//   integrate flux over the face: W in both tangential directions
//   distribute (Deltatil): node l -= , node l+1 +=

#include "mars_cvfem_ho_basis.hpp"
#include <vector>
#include <array>

namespace mars {
namespace fem {

// Reference-cube corner signs, matching generate_hex_cube hex order
// (0:(-1,-1,-1) 1:(1,-1,-1) 2:(1,1,-1) 3:(-1,1,-1) 4:(-1,-1,1) 5:(1,-1,1) 6:(1,1,1) 7:(-1,1,1)).
inline const int (&hexCornerRef())[8][3] {
    static const int s[8][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
    return s;
}

// Trilinear Jacobian J[a][b] = d x_a / d r_b at reference point r, from 8 corners.
inline void trilinearJacobian(const double corners[8][3], const double r[3], double J[3][3])
{
    const auto& S = hexCornerRef();
    for (int a=0;a<3;++a) for (int b=0;b<3;++b) J[a][b]=0.0;
    for (int c=0;c<8;++c) {
        double sh[3];   // d N_c / d r_b
        sh[0]=0.125*S[c][0]*(1+S[c][1]*r[1])*(1+S[c][2]*r[2]);
        sh[1]=0.125*(1+S[c][0]*r[0])*S[c][1]*(1+S[c][2]*r[2]);
        sh[2]=0.125*(1+S[c][0]*r[0])*(1+S[c][1]*r[1])*S[c][2];
        for (int a=0;a<3;++a) for (int b=0;b<3;++b) J[a][b]+=sh[b]*corners[c][a];
    }
}

// Per-element metric for the general apply: at each (dir, SCS face l, tangential
// s,r) point [normal = Gauss xi_l, tangential = GLL zeta_s,zeta_r], store the
// 3-vector G = detJ (JJ^T)^{-1} e_dir, decomposed into (normal, tang1, tang2)
// components used by Alg 2. Layout: Gout[((dir*p + l)*n + s)*n + r][0..2] with
// [2]=normal, [1]=tang1 (s-axis), [0]=tang2 (r-axis).
inline std::vector<std::array<double,3>>
computeElementMetric(const HoCvfemOperators& op, const double corners[8][3])
{
    const int p=op.p, n=p+1;
    // tangential axes for each direction (matches the idx() mapping in the apply)
    const int t1axis[3]={1,0,0}, t2axis[3]={2,2,1};
    std::vector<std::array<double,3>> G((size_t)3*p*n*n);
    for (int dir=0; dir<3; ++dir)
        for (int l=0;l<p;++l) for (int s=0;s<n;++s) for (int r=0;r<n;++r) {
            double rf[3];
            rf[dir]        = op.xi[l];
            rf[t1axis[dir]]= op.zeta[s];
            rf[t2axis[dir]]= op.zeta[r];
            double J[3][3]; trilinearJacobian(corners, rf, J);
            double det = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
                       - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
                       + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
            // J^{-1} (cofactor / det)
            double Ji[3][3];
            Ji[0][0]=(J[1][1]*J[2][2]-J[1][2]*J[2][1])/det; Ji[0][1]=(J[0][2]*J[2][1]-J[0][1]*J[2][2])/det; Ji[0][2]=(J[0][1]*J[1][2]-J[0][2]*J[1][1])/det;
            Ji[1][0]=(J[1][2]*J[2][0]-J[1][0]*J[2][2])/det; Ji[1][1]=(J[0][0]*J[2][2]-J[0][2]*J[2][0])/det; Ji[1][2]=(J[0][2]*J[1][0]-J[0][0]*J[1][2])/det;
            Ji[2][0]=(J[1][0]*J[2][1]-J[1][1]*J[2][0])/det; Ji[2][1]=(J[0][1]*J[2][0]-J[0][0]*J[2][1])/det; Ji[2][2]=(J[0][0]*J[1][1]-J[0][1]*J[1][0])/det;
            // isoparametric Laplacian metric = detJ * J^{-1} J^{-T};
            // Gvec = metric * e_dir = detJ * (J^{-1}J^{-T})[:,dir]
            double gvec[3];
            for (int a=0;a<3;++a) { double v=0; for(int k=0;k<3;++k) v+=Ji[a][k]*Ji[dir][k]; gvec[a]=det*v; }
            auto& g = G[((size_t)(dir*p + l)*n + s)*n + r];
            g[2]=gvec[dir]; g[1]=gvec[t1axis[dir]]; g[0]=gvec[t2axis[dir]];
        }
    return G;
}

// General element diffusion apply (full Alg 2 with non-orthogonal cross terms),
// using the per-element metric G from computeElementMetric. Works for arbitrary
// straight-sided hexes. The affine cube is the special case G[0]=G[1]=0.
inline void applyHoCvfemElement(const HoCvfemOperators& op,
                                const std::vector<std::array<double,3>>& G,
                                const std::vector<double>& u, std::vector<double>& y)
{
    const int p=op.p, n=p+1, n3=n*n*n;
    y.assign(n3, 0.0);
    auto idx = [&](int dir,int nrm,int t1,int t2)->int {
        if (dir==0) return nrm*n*n + t1*n + t2;
        if (dir==1) return t1*n*n + nrm*n + t2;
        return t1*n*n + t2*n + nrm;
    };
    for (int dir=0; dir<3; ++dir)
        for (int l=0;l<p;++l) {
            std::vector<double> interp(n*n,0.0), deriv(n*n,0.0), flux(n*n,0.0), tmp(n*n,0.0), intf(n*n,0.0);
            for (int s=0;s<n;++s) for (int r=0;r<n;++r) {
                double bi=0, di=0;
                for (int q=0;q<n;++q) { bi+=op.Btil[l*n+q]*u[idx(dir,q,s,r)]; di+=op.Dtil[l*n+q]*u[idx(dir,q,s,r)]; }
                interp[s*n+r]=bi; deriv[s*n+r]=di;
            }
            for (int s=0;s<n;++s) for (int r=0;r<n;++r) {
                const auto& g = G[((size_t)(dir*p + l)*n + s)*n + r];
                double dt2=0, dt1=0;                       // tangential derivatives of interp
                for (int q=0;q<n;++q) dt2 += op.D[r*n+q]*interp[s*n+q];   // along t2 (r-axis)
                for (int q=0;q<n;++q) dt1 += op.D[s*n+q]*interp[q*n+r];   // along t1 (s-axis)
                flux[s*n+r] = g[2]*deriv[s*n+r] + g[0]*dt2 + g[1]*dt1;
            }
            for (int s=0;s<n;++s) for (int r=0;r<n;++r) { double v=0; for(int q=0;q<n;++q) v+=op.W[r*n+q]*flux[s*n+q]; tmp[s*n+r]=v; }
            for (int s=0;s<n;++s) for (int r=0;r<n;++r) { double v=0; for(int q=0;q<n;++q) v+=op.W[s*n+q]*tmp[q*n+r]; intf[s*n+r]=v; }
            for (int s=0;s<n;++s) for (int r=0;r<n;++r) {
                y[idx(dir,l,  s,r)] -= intf[s*n+r];
                y[idx(dir,l+1,s,r)] += intf[s*n+r];
            }
        }
}

inline void applyHoCvfemElementCube(const HoCvfemOperators& op, double coeff,
                                    const std::vector<double>& u, std::vector<double>& y)
{
    const int p = op.p, n = p + 1, n3 = n * n * n;
    y.assign(n3, 0.0);

    // (a,b,c) flat index with `nrm` in the normal slot for `dir`, (t1,t2) tangential.
    auto idx = [&](int dir, int nrm, int t1, int t2) -> int {
        if (dir == 0) return nrm * n * n + t1 * n + t2;
        if (dir == 1) return t1 * n * n + nrm * n + t2;
        return t1 * n * n + t2 * n + nrm;
    };

    for (int dir = 0; dir < 3; ++dir)
        for (int l = 0; l < p; ++l)               // subcontrol-surface faces
        {
            std::vector<double> flux(n * n, 0.0), w1(n * n, 0.0), intf(n * n, 0.0);
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {
                    double deriv = 0;
                    for (int q = 0; q < n; ++q) deriv += op.Dtil[l * n + q] * u[idx(dir, q, s, r)];
                    flux[s * n + r] = coeff * deriv;
                }
            // tangential integration (W in each tangential direction)
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {
                    double v = 0; for (int q = 0; q < n; ++q) v += op.W[r * n + q] * flux[s * n + q];
                    w1[s * n + r] = v;
                }
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {
                    double v = 0; for (int q = 0; q < n; ++q) v += op.W[s * n + q] * w1[q * n + r];
                    intf[s * n + r] = v;
                }
            // distribute to the two nodes bounding the face (Deltatil -1/+1)
            for (int s = 0; s < n; ++s)
                for (int r = 0; r < n; ++r) {
                    y[idx(dir, l,     s, r)] -= intf[s * n + r];
                    y[idx(dir, l + 1, s, r)] += intf[s * n + r];
                }
        }
}

} // namespace fem
} // namespace mars
