#ifndef MARS_HO_DOF_HANDLER_TET_HPP
#define MARS_HO_DOF_HANDLER_TET_HPP

// C0 nodal DOF numbering for high-order tetrahedra. Tet counterpart of HODofHandler (hex).
// Each nodal DOF is keyed by its barycentric signature against the GLOBAL vertex ids:
// sorted [(gid, weight)] pairs with zero weights dropped. Shared vertices/edges/faces then match
// automatically in any element orientation -- no edge-direction or face-rotation tables needed.
// Host build; the resulting elemDof feeds the device gather/scatter.

#include "mars_ho_tet_basis.hpp"
#include <algorithm>
#include <map>
#include <vector>
#include <array>
#include <utility>
#include <cstdint>

namespace mars {
namespace fem {

struct HoTetDofHandler {
    int P = 0, Np = 0;
    long numDof = 0;
    std::vector<int> elemDof;                       // [nElem * Np] local node -> global DOF
    std::vector<std::array<double,3>> dofPos;       // physical position per DOF (for gates/BCs)

    // elemCorners: [nElem*4] global vertex ids; coords: [nVerts][3].
    template <typename RealType>
    void build(const std::vector<int>& elemCorners, const std::vector<std::array<double,3>>& coords,
               const HoTetNodal<RealType>& nd) {
        P = nd.P; Np = nd.Np;
        const size_t nElem = elemCorners.size() / 4;
        elemDof.assign(nElem * Np, -1);
        dofPos.clear();
        std::map<std::vector<std::pair<int,int>>, int> dofOf;
        for (size_t e = 0; e < nElem; ++e) {
            const int g[4] = { elemCorners[e*4+0], elemCorners[e*4+1], elemCorners[e*4+2], elemCorners[e*4+3] };
            for (int nn = 0; nn < Np; ++nn) {
                const int i = nd.bary[nn][0], j = nd.bary[nn][1], k = nd.bary[nn][2];
                const int w[4] = { P - i - j - k, i, j, k };   // weights vs corners 0..3
                std::vector<std::pair<int,int>> key;
                for (int c = 0; c < 4; ++c) if (w[c] > 0) key.push_back({g[c], w[c]});
                std::sort(key.begin(), key.end());
                auto it = dofOf.find(key);
                int id;
                if (it == dofOf.end()) {
                    id = (int)numDof++;
                    dofOf.emplace(key, id);
                    std::array<double,3> x = {0,0,0};
                    for (int c = 0; c < 4; ++c) for (int d = 0; d < 3; ++d)
                        x[d] += (w[c] / (double)P) * coords[g[c]][d];
                    dofPos.push_back(x);
                } else id = it->second;
                elemDof[e*Np + nn] = id;
            }
        }
    }
};

// C0 DOF numbering for the CVFEM collapsed GLL tensor grid. Nodes are keyed by their barycentric
// weights against the element's SORTED global vertex ids, quantized to 2^32 (GLL positions are
// irrational; both sides of a shared face compute identical doubles, so quantization is safe).
// Coincident collapsed nodes (the degenerate apex/edge images) get identical keys and are MERGED
// into one DOF. Requires ascending-by-gid corner ordering per element (sortedCorners output) --
// that ordering makes every shared face induce the same canonical 2D Duffy grid from both sides.
struct HoCvfemTetDofHandler {
    int P = 0, n = 0, NN = 0;
    long numDof = 0;
    std::vector<int> elemDof;                       // [nElem * NN] tensor node -> global DOF
    std::vector<int> sortedCorners;                 // [nElem * 4] ascending-gid corner order
    std::vector<std::array<double,3>> dofPos;       // physical position per DOF

    void build(const std::vector<int>& elemCorners, const std::vector<std::array<double,3>>& coords,
               const std::vector<double>& Z /* n GLL nodes on [0,1] */) {
        n = (int)Z.size(); P = n - 1; NN = n*n*n;
        const size_t nElem = elemCorners.size() / 4;
        elemDof.assign(nElem * NN, -1);
        sortedCorners.assign(nElem * 4, -1);
        dofPos.clear(); numDof = 0;
        std::map<std::vector<std::pair<int,long long>>, int> dofOf;
        const double Q = 4294967296.0;   // 2^32 quantization
        for (size_t e = 0; e < nElem; ++e) {
            int g[4] = { elemCorners[e*4+0], elemCorners[e*4+1], elemCorners[e*4+2], elemCorners[e*4+3] };
            std::sort(g, g + 4);
            for (int c = 0; c < 4; ++c) sortedCorners[e*4+c] = g[c];
            for (int ia = 0; ia < n; ++ia) for (int ib = 0; ib < n; ++ib) for (int ic = 0; ic < n; ++ic) {
                const double a = Z[ia], b = Z[ib], cc = Z[ic];
                const double r = a*(1-b)*(1-cc), s = b*(1-cc), t = cc;
                const double w[4] = { 1.0-r-s-t, r, s, t };
                std::vector<std::pair<int,long long>> key;
                for (int c = 0; c < 4; ++c) if (w[c] > 1e-9) key.push_back({g[c], (long long)std::llround(w[c]*Q)});
                std::sort(key.begin(), key.end());
                auto it = dofOf.find(key);
                int id;
                if (it == dofOf.end()) {
                    id = (int)numDof++;
                    dofOf.emplace(key, id);
                    std::array<double,3> x = {0,0,0};
                    for (int c = 0; c < 4; ++c) for (int d = 0; d < 3; ++d) x[d] += w[c]*coords[g[c]][d];
                    dofPos.push_back(x);
                } else id = it->second;
                elemDof[e*NN + ((size_t)ia*n+ib)*n+ic] = id;
            }
        }
    }
};

// Kuhn triangulation of an ncells^3 cube grid: 6 tets per cell (all sharing the cell diagonal),
// face-to-face conforming across cells. Test-mesh utility for the tet HO gates.
inline void buildKuhnTetMesh(int ncells, std::vector<std::array<double,3>>& coords,
                             std::vector<int>& elemCorners) {
    const int N = ncells + 1;
    coords.clear(); elemCorners.clear();
    for (int z = 0; z < N; ++z) for (int y = 0; y < N; ++y) for (int x = 0; x < N; ++x)
        coords.push_back({(double)x, (double)y, (double)z});
    auto gid = [&](int x, int y, int z) { return (z*N + y)*N + x; };
    const int perms[6][3] = {{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
    for (int z = 0; z < ncells; ++z) for (int y = 0; y < ncells; ++y) for (int x = 0; x < ncells; ++x)
        for (int p = 0; p < 6; ++p) {
            int c[3] = {x, y, z};
            elemCorners.push_back(gid(c[0], c[1], c[2]));
            for (int step = 0; step < 3; ++step) {
                c[perms[p][step]] += 1;
                elemCorners.push_back(gid(c[0], c[1], c[2]));
            }
        }
}

}  // namespace fem
}  // namespace mars
#endif
