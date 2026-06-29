// Stage 2: validate the ORIENTATION-AWARE HO DOF handler on a mesh where the two
// hexes sharing a face genuinely disagree on the face parameterization.
//
// We take the structured ExExE cube, then relabel EACH element's 8 local corners
// by a RANDOM cube rotation (one of the 24 proper rotations of the cube). The
// geometry stays a valid axis-aligned cube -- only WHICH physical corner is local
// 0..7 changes -- so neighbors now present each shared face in different dihedral
// orientations. Connectivity is built from the relabeled corners, and the metric
// is computed from the actual edge vectors (general Jacobian, not the
// assume-c0c1-is-x shortcut), so a wrong corner labeling cannot be hidden by a
// shortcut metric.
//
// With the orientation-BLIND face pos this FAILS the continuity gate at p>=3
// (face-interior DOFs are scattered to mismatched slots). With the canonical
// face pos (hex_face_canonical_pos) it PASSES. p=2 has a single interior face DOF
// (no permutation freedom) and is a passing sanity baseline either way.
//
// Gates (same as mars_ho_dof_test.cu):
//   1. numDof == (E*p+1)^3        (entity enumeration completeness, no double-count)
//   2. nullspace  A*1 = 0
//   3. symmetry   u^T A v = v^T A u
//   4. continuity A*(physical x) = 0 at interior DOFs   <-- the orientation gate

#include "backend/distributed/unstructured/fem/mars_spectral_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_laplacian.hpp"

#include <cuda_runtime.h>
#include <array>
#include <cstdio>
#include <cmath>
#include <random>
#include <vector>

using namespace mars::fem;

#define CK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); return false; } } while(0)

// Unit-cube positions of the 8 hex corners, in the handler's corner convention.
static const int CPX[8] = {0,1,1,0,0,1,1,0};
static const int CPY[8] = {0,0,1,1,0,0,1,1};
static const int CPZ[8] = {0,0,0,0,1,1,1,1};

// The 24 proper cube rotations as corner permutations: rot[c] = the OLD local
// corner that becomes the NEW local slot c. Built from signed-permutation 3x3
// matrices with det = +1 (reflections excluded -> geometry stays a real cube,
// not a mirror image). Each rotation maps a corner's centered position to another
// corner's position; we look up which old corner lands on each new slot.
static std::vector<std::array<int,8>> buildCubeRotations()
{
    std::vector<std::array<int,8>> rots;
    static const int PERM[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    for (int sx = -1; sx <= 1; sx += 2)
    for (int sy = -1; sy <= 1; sy += 2)
    for (int sz = -1; sz <= 1; sz += 2)
    for (int pm = 0; pm < 6; ++pm)
    {
        const int sgn[3] = {sx, sy, sz};
        // R: column j is sgn[j] * e_{PERM[pm][j]}.
        int R[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        for (int j = 0; j < 3; ++j) R[PERM[pm][j]][j] = sgn[j];
        const int det = R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1])
                      - R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0])
                      + R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]);
        if (det < 0) continue;                    // proper rotations only

        std::array<int,8> p8{};
        for (int c = 0; c < 8; ++c)
        {
            const int u[3] = {2*CPX[c]-1, 2*CPY[c]-1, 2*CPZ[c]-1};   // centered {-1,+1}
            const int rx = R[0][0]*u[0]+R[0][1]*u[1]+R[0][2]*u[2];
            const int ry = R[1][0]*u[0]+R[1][1]*u[1]+R[1][2]*u[2];
            const int rz = R[2][0]*u[0]+R[2][1]*u[1]+R[2][2]*u[2];
            const int bx = (rx > 0), by = (ry > 0), bz = (rz > 0);
            int oc = -1;
            for (int q = 0; q < 8; ++q)
                if (CPX[q]==bx && CPY[q]==by && CPZ[q]==bz) { oc = q; break; }
            p8[c] = oc;
        }
        rots.push_back(p8);
    }
    return rots;
}

template<int P>
bool runStage2(int E)
{
    constexpr int n  = P + 1;
    constexpr int N3 = n * n * n;
    using Real = double;
    const Real h = 1.0 / E;

    auto cgid = [&](int cx, int cy, int cz) { return (cx * (E + 1) + cy) * (E + 1) + cz; };

    // Structured connectivity, then a random cube rotation applied per element.
    auto rots = buildCubeRotations();
    std::mt19937 rng(2026);
    std::uniform_int_distribution<int> pickRot(0, (int)rots.size() - 1);

    std::vector<std::array<int,8>> elemCorners;
    elemCorners.reserve(size_t(E) * E * E);
    for (int ex = 0; ex < E; ++ex)
    for (int ey = 0; ey < E; ++ey)
    for (int ez = 0; ez < E; ++ez)
    {
        std::array<int,8> base = {
            cgid(ex,ey,ez),     cgid(ex+1,ey,ez),   cgid(ex+1,ey+1,ez),   cgid(ex,ey+1,ez),
            cgid(ex,ey,ez+1),   cgid(ex+1,ey,ez+1), cgid(ex+1,ey+1,ez+1), cgid(ex,ey+1,ez+1) };
        const auto& r = rots[pickRot(rng)];
        std::array<int,8> relabeled{};
        for (int c = 0; c < 8; ++c) relabeled[c] = base[r[c]];
        elemCorners.push_back(relabeled);
    }
    long nCorner = long(E + 1) * (E + 1) * (E + 1);

    HODofHandler dh;
    dh.build(elemCorners, nCorner, P);
    long nDof = dh.numDof;
    long expect = long(E * P + 1) * (E * P + 1) * (E * P + 1);

    std::vector<double> x, w, D;
    gllBasis(P, x, w, D);
    std::vector<Real> Dr(D.begin(), D.end());

    // Recover a corner's grid position from its global id (inverse of cgid).
    auto idToXYZ = [&](int gid, Real& X, Real& Y, Real& Z) {
        int cz = gid % (E + 1);
        int cy = (gid / (E + 1)) % (E + 1);
        int cx = gid / ((E + 1) * (E + 1));
        X = cx * h; Y = cy * h; Z = cz * h;
    };

    size_t nEl = elemCorners.size();
    std::vector<Real>    G(nEl * N3 * 6, 0.0);
    std::vector<Real>    dofX(nDof, 0.0);
    std::vector<uint8_t> isBdry(nDof, 0);

    // Per-element general metric from the 3 edge vectors of the relabeled cube:
    //   e1 = X[c1]-X[c0], e2 = X[c3]-X[c0], e3 = X[c4]-X[c0]   (ref axes span [-1,1])
    //   J = [e1 e2 e3]/2 ;  G = w * (J^-1 J^-T) * |det J|   per GLL node.
    // Under a cube relabeling J is an orthogonal scaled (signed-permutation)
    // matrix, so G stays the correct diagonal -- but we compute it generally so a
    // wrong Jacobian (or a wrong corner labeling) cannot pass.
    for (size_t e = 0; e < nEl; ++e)
    {
        const auto& c = elemCorners[e];
        Real Xc[8][3];
        for (int cc = 0; cc < 8; ++cc) idToXYZ(c[cc], Xc[cc][0], Xc[cc][1], Xc[cc][2]);

        Real J[3][3];
        for (int d = 0; d < 3; ++d) {
            J[d][0] = 0.5 * (Xc[1][d] - Xc[0][d]);   // e1/2
            J[d][1] = 0.5 * (Xc[3][d] - Xc[0][d]);   // e2/2
            J[d][2] = 0.5 * (Xc[4][d] - Xc[0][d]);   // e3/2
        }
        Real det = J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])
                 - J[0][1]*(J[1][0]*J[2][2]-J[1][2]*J[2][0])
                 + J[0][2]*(J[1][0]*J[2][1]-J[1][1]*J[2][0]);
        Real adet = std::abs(det);
        Real Ji[3][3];
        Ji[0][0]=(J[1][1]*J[2][2]-J[1][2]*J[2][1])/det;
        Ji[0][1]=(J[0][2]*J[2][1]-J[0][1]*J[2][2])/det;
        Ji[0][2]=(J[0][1]*J[1][2]-J[0][2]*J[1][1])/det;
        Ji[1][0]=(J[1][2]*J[2][0]-J[1][0]*J[2][2])/det;
        Ji[1][1]=(J[0][0]*J[2][2]-J[0][2]*J[2][0])/det;
        Ji[1][2]=(J[0][2]*J[1][0]-J[0][0]*J[1][2])/det;
        Ji[2][0]=(J[1][0]*J[2][1]-J[1][1]*J[2][0])/det;
        Ji[2][1]=(J[0][1]*J[2][0]-J[0][0]*J[2][1])/det;
        Ji[2][2]=(J[0][0]*J[1][1]-J[0][1]*J[1][0])/det;
        Real M[3][3];
        for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            Real s = 0; for (int kk = 0; kk < 3; ++kk) s += Ji[a][kk]*Ji[b][kk];
            M[a][b] = s;                              // J^-1 J^-T
        }

        for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
        for (int k = 0; k < n; ++k) {
            int pidx = i*n*n + j*n + k;
            Real wp = w[i]*w[j]*w[k]*adet;
            Real* g = &G[(e*N3 + pidx) * 6];
            g[0]=wp*M[0][0]; g[1]=wp*M[0][1]; g[2]=wp*M[0][2];
            g[3]=wp*M[1][1]; g[4]=wp*M[1][2]; g[5]=wp*M[2][2];

            // Physical coord via trilinear map of the relabeled corners.
            Real a0=0.5*(1-x[i]), a1=0.5*(1+x[i]);
            Real b0=0.5*(1-x[j]), b1=0.5*(1+x[j]);
            Real c0=0.5*(1-x[k]), c1=0.5*(1+x[k]);
            Real px=0, py=0, pz=0;
            for (int cc = 0; cc < 8; ++cc) {
                Real sh = (CPX[cc]?a1:a0)*(CPY[cc]?b1:b0)*(CPZ[cc]?c1:c0);
                px += sh*Xc[cc][0]; py += sh*Xc[cc][1]; pz += sh*Xc[cc][2];
            }
            int dof = dh.elemDof[e*N3 + pidx];
            dofX[dof] = px;
            bool b = (px<1e-9||px>1-1e-9||py<1e-9||py>1-1e-9||pz<1e-9||pz>1-1e-9);
            isBdry[dof] = b ? 1 : 0;
        }
    }

    int *dElemDof; Real *dD, *dG, *dU, *dY;
    CK(cudaMalloc(&dElemDof, dh.elemDof.size()*sizeof(int)));
    CK(cudaMalloc(&dD, Dr.size()*sizeof(Real)));
    CK(cudaMalloc(&dG, G.size()*sizeof(Real)));
    CK(cudaMalloc(&dU, nDof*sizeof(Real)));
    CK(cudaMalloc(&dY, nDof*sizeof(Real)));
    CK(cudaMemcpy(dElemDof, dh.elemDof.data(), dh.elemDof.size()*sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dD, Dr.data(), Dr.size()*sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dG, G.data(), G.size()*sizeof(Real), cudaMemcpyHostToDevice));

    int block = 128, grid = int((nEl + block - 1)/block);
    auto apply = [&](const std::vector<Real>& u, std::vector<Real>& y) -> bool {
        CK(cudaMemcpy(dU, u.data(), nDof*sizeof(Real), cudaMemcpyHostToDevice));
        CK(cudaMemset(dY, 0, nDof*sizeof(Real)));
        ho_laplacian_apply_cg<Real, P><<<grid, block>>>(dU, dY, dElemDof, dD, dG, nEl);
        CK(cudaDeviceSynchronize());
        y.resize(nDof);
        CK(cudaMemcpy(y.data(), dY, nDof*sizeof(Real), cudaMemcpyDeviceToHost));
        return true;
    };

    std::vector<Real> ones(nDof, 1.0), y;
    if (!apply(ones, y)) return false;
    Real nullMax = 0; for (Real v : y) nullMax = std::max(nullMax, std::abs(v));

    std::mt19937 r2(7); std::uniform_real_distribution<Real> uni(-1,1);
    std::vector<Real> uu(nDof), vv(nDof), Au, Av;
    for (auto& z : uu) z = uni(r2);
    for (auto& z : vv) z = uni(r2);
    if (!apply(uu, Au)) return false;
    if (!apply(vv, Av)) return false;
    Real uAv=0, vAu=0; for (long p=0;p<nDof;++p){ uAv+=uu[p]*Av[p]; vAu+=vv[p]*Au[p]; }
    Real symErr = std::abs(uAv-vAu)/(std::abs(uAv)+1.0);

    std::vector<Real> Ax;
    if (!apply(dofX, Ax)) return false;
    Real interMax = 0;
    for (long p=0;p<nDof;++p) if (!isBdry[p]) interMax = std::max(interMax, std::abs(Ax[p]));

    bool ok = (nDof==expect) && (nullMax<1e-9) && (symErr<1e-10) && (interMax<1e-9);
    printf("P=%d E=%d | nDof=%ld (expect %ld) | nEdge=%ld nFace=%ld | null=%.2e | sym=%.2e | continuity=%.2e | %s\n",
           P, E, nDof, expect, dh.nEdge, dh.nFace, nullMax, symErr, interMax, ok?"PASS":"FAIL");

    cudaFree(dElemDof); cudaFree(dD); cudaFree(dG); cudaFree(dU); cudaFree(dY);
    return ok;
}

int main()
{
    bool ok = true;
    // p>=3 cases carry interior-face permutation freedom: these catch a wrong
    // orientation table. p=2 has a single interior face DOF -> passing baseline.
    ok &= runStage2<2>(4);
    ok &= runStage2<4>(3);
    ok &= runStage2<7>(2);
    printf("Stage 2 (orientation-aware DOF handler, relabeled-cube validation): %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
