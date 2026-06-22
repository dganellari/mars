#pragma once

// Topological high-order DOF numbering for tensor-product hex spectral elements.
// Builds the element-local -> global DOF map (elemDof[e*N3 + l]) that the CG
// matrix-free apply (ho_laplacian_apply_cg) gathers/scatters through. Continuity
// is established by keying shared entities (edges, faces) on the GLOBAL CORNER
// IDS they connect -- never on high-order node coordinates (which would hit the
// SFC quantization wall at p=7). Corner continuity comes from the mesh's
// existing corner-node numbering.
//
// DOF layout (global):
//   [0, nCorner)                              corner DOFs (= corner node ids)
//   [nCorner, +nEdge*(p-1))                   edge-interior DOFs
//   [+, +nFace*(p-1)^2)                       face-interior DOFs
//   [+, +nElem*(p-1)^3)                       element-interior DOFs
//
// Orientation:
//   - edges: (p-1) interior DOFs ordered from the lower-global-id endpoint to
//     the higher. An element traversing the edge the other way reverses. (done)
//   - faces: (p-1)^2 interior DOFs are placed in a canonical 2D frame derived
//     from the 4 global corner ids (hex_face_canonical_pos), so the two hexes
//     sharing a face agree on each interior node's slot under any of the 8
//     dihedral relabelings -- continuity holds for general unstructured meshes.
//
// Hex corner convention (matches generate_hex_cube.py / CVFEM connectivity):
//   0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0) 4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)

#include "mars_hex_face_orientations.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <vector>

namespace mars {
namespace fem {

// Local corner index from the (i==p, j==p, k==p) bits.
inline int hexCornerIndex(int bi, int bj, int bk)
{
    static const int tbl[2][2][2] = {
        // [bk][bj][bi]
        { {0, 1}, {3, 2} },   // bk=0
        { {4, 5}, {7, 6} },   // bk=1
    };
    return tbl[bk][bj][bi];
}

class HODofHandler
{
public:
    int  P = 0, n = 0, N3 = 0;
    long nCorner = 0, nEdge = 0, nFace = 0, nElem = 0;
    long numDof  = 0;
    std::vector<int> elemDof;   // [nElem * N3] -> global DOF id

    // elemCorners[e] = 8 global corner-node ids in hex order.
    void build(const std::vector<std::array<int, 8>>& elemCorners, long numCornerNodes, int order)
    {
        P  = order;
        n  = P + 1;
        N3 = n * n * n;
        nElem   = (long)elemCorners.size();
        nCorner = numCornerNodes;
        const int pm1 = P - 1;                 // interior count per edge dim

        // Hex topology in the corner convention above.
        static const int EDGES[12][2] = {
            {0,1},{1,2},{2,3},{3,0}, {4,5},{5,6},{6,7},{7,4}, {0,4},{1,5},{2,6},{3,7}
        };
        static const int FACES[6][4] = {
            {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, {3,2,6,7}, {1,2,6,5}, {0,3,7,4}
        };

        // --- Pass 1: enumerate unique edges & faces by sorted global-corner key ---
        std::map<std::array<int,2>, int> edgeId;
        std::map<std::array<int,4>, int> faceId;
        for (long e = 0; e < nElem; ++e)
        {
            const auto& c = elemCorners[e];
            for (int ed = 0; ed < 12; ++ed)
            {
                std::array<int,2> k{ c[EDGES[ed][0]], c[EDGES[ed][1]] };
                if (k[0] > k[1]) std::swap(k[0], k[1]);
                edgeId.emplace(k, (int)edgeId.size());
            }
            for (int f = 0; f < 6; ++f)
            {
                std::array<int,4> k{ c[FACES[f][0]], c[FACES[f][1]], c[FACES[f][2]], c[FACES[f][3]] };
                std::sort(k.begin(), k.end());
                faceId.emplace(k, (int)faceId.size());
            }
        }
        nEdge = (long)edgeId.size();
        nFace = (long)faceId.size();

        const long edgeBase = nCorner;
        const long faceBase = edgeBase + nEdge * pm1;
        const long intrBase = faceBase + nFace * (long)pm1 * pm1;
        numDof = intrBase + nElem * (long)pm1 * pm1 * pm1;

        // --- Pass 2: fill elemDof ---
        elemDof.assign(nElem * N3, -1);
        for (long e = 0; e < nElem; ++e)
        {
            const auto& c = elemCorners[e];
            for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
            {
                int l = i * n * n + j * n + k;
                int onI = (i == 0 || i == P), onJ = (j == 0 || j == P), onK = (k == 0 || k == P);
                int cnt = onI + onJ + onK;
                long dof = -1;

                if (cnt == 3)
                {
                    dof = c[hexCornerIndex(i == P, j == P, k == P)];
                }
                else if (cnt == 2)
                {
                    // edge: one index varies (the non-on one)
                    int vary = onI ? (onJ ? 2 : 1) : 0;        // index that varies
                    int t    = (vary == 0) ? i : (vary == 1) ? j : k;   // 1..p-1
                    // endpoints: vary=0 -> corner, vary=p -> corner
                    int bi = (vary == 0) ? 0 : (i == P), bj = (vary == 1) ? 0 : (j == P), bk = (vary == 2) ? 0 : (k == P);
                    int cA = c[hexCornerIndex(bi, bj, bk)];
                    bi = (vary == 0) ? 1 : (i == P); bj = (vary == 1) ? 1 : (j == P); bk = (vary == 2) ? 1 : (k == P);
                    int cB = c[hexCornerIndex(bi, bj, bk)];
                    std::array<int,2> key{ cA, cB };
                    int rev = 0;
                    if (key[0] > key[1]) { std::swap(key[0], key[1]); rev = 1; }
                    int eid = edgeId.at(key);
                    int pos = rev ? (pm1 - t) : (t - 1);       // canonical low->high
                    dof = edgeBase + (long)eid * pm1 + pos;
                }
                else if (cnt == 1)
                {
                    // face: one index on a boundary, two vary
                    int on = onI ? 0 : onJ ? 1 : 2;
                    int v1, v2;                                 // the two varying values (1..p-1)
                    if (on == 0)      { v1 = j; v2 = k; }
                    else if (on == 1) { v1 = i; v2 = k; }
                    else              { v1 = i; v2 = j; }
                    const int* fc = FACES[on == 0 ? (i == P ? 4 : 5) : on == 1 ? (j == P ? 3 : 2) : (k == P ? 1 : 0)];
                    std::array<int,4> key{ c[fc[0]], c[fc[1]], c[fc[2]], c[fc[3]] };
                    std::sort(key.begin(), key.end());
                    int fid = faceId.at(key);
                    // Canonical orientation: both neighbors of a shared face elect
                    // the same origin/s/t frame from the 4 global ids, so this pos
                    // agrees across the dihedral relabeling (continuity at p>=3).
                    int pos = hex_face_canonical_pos(P, v1, v2, c[fc[0]], c[fc[1]], c[fc[2]], c[fc[3]]);
                    dof = faceBase + (long)fid * pm1 * pm1 + pos;
                }
                else
                {
                    // interior: element-local block
                    int pos = ((i - 1) * pm1 + (j - 1)) * pm1 + (k - 1);
                    dof = intrBase + e * (long)pm1 * pm1 * pm1 + pos;
                }
                elemDof[e * N3 + l] = (int)dof;
            }
        }
    }

    // Canonical cross-rank identity of a high-order DOF (for the multi-rank halo).
    // kind: 0 corner, 1 edge, 2 face, 3 interior. g* = sorted global defining-corner
    // ids (unused = -1); pos = index within the edge/face/interior block. Two ranks
    // that share a DOF compute the SAME key, so the halo matches by it.
    struct DofKey { int kind; long g0, g1, g2, g3; int pos; };

    std::vector<int>    dofOwner;   // [numDof] owning rank (filled by buildDistributed)
    std::vector<DofKey> dofKey;     // [numDof] canonical cross-rank key
    std::vector<uint8_t> dofShared; // [numDof] 1 if edge/face DOF on a rank boundary
                                    // (provisional owner = myRank; resolve via
                                    //  resolveHoDofOwnership before building the halo)
    std::vector<uint8_t> dofBoundary; // [numDof] 1 if ANY shared DOF (corner/edge/face on
                                      // a rank boundary). Superset of dofShared (adds shared
                                      // corners). Only these can be exchanged -> the halo
                                      // keys/maps over dofBoundary, NOT all numDof (else a
                                      // numDof-sized std::map OOMs the host at scale).

    // Multi-rank numbering. Calls build() for the LOCAL dense elemDof, then tags each
    // DOF with owner + canonical key. Ownership:
    //   corner   -> cornerOwner[lc]  (consistent cstone P1 node ownership; the working
    //               p=1 case -- the owner always contains the node)
    //   interior -> elemOwner[e]     (element-local, never shared)
    //   edge/face-> PROVISIONAL myRank. The lowest-global-corner heuristic is WRONG
    //               here: that corner's owner may have the corner but not the edge/face,
    //               orphaning the DOF. resolveHoDofOwnership() (a peer key exchange)
    //               sets the real owner = min rank among the ranks that actually hold
    //               it. dofShared[d]=1 flags edge/face DOF whose defining corners are
    //               all P1-shared -> candidates for that resolution.
    //   elemCorners[e] : 8 LOCAL corner indices (hex order)  -- for numbering
    //   cornerGid[lc]  : LOCAL corner index -> GLOBAL id      -- for identity + key
    //   cornerOwner[lc]: LOCAL corner index -> owning rank    -- P1 node ownership
    //   elemOwner[e]   : owning rank of element e
    //   sharedCorner[lc]: 1 if LOCAL corner lc is on a rank boundary (P1 send U recv)
    void buildDistributed(const std::vector<std::array<int,8>>& elemCorners,
                          long numCornerNodes, int order,
                          const std::vector<long>&    cornerGid,
                          const std::vector<int>&     cornerOwner,
                          const std::vector<int>&     elemOwner,
                          int                          myRank,
                          const std::vector<uint8_t>&  sharedCorner)
    {
        build(elemCorners, numCornerNodes, order);
        dofOwner.assign(numDof, -1);
        dofKey.assign(numDof, DofKey{-1,-1,-1,-1,-1,-1});
        dofShared.assign(numDof, 0);
        dofBoundary.assign(numDof, 0);

        static const int FACES[6][4] = {
            {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, {3,2,6,7}, {1,2,6,5}, {0,3,7,4} };
        const int pm1 = P - 1;

        for (long e = 0; e < nElem; ++e)
        {
            const auto& c = elemCorners[e];
            for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
            {
                int l   = i*n*n + j*n + k;
                int dof = elemDof[e*N3 + l];
                int onI = (i==0||i==P), onJ=(j==0||j==P), onK=(k==0||k==P);
                int cnt = onI + onJ + onK;
                DofKey key{-1,-1,-1,-1,-1,-1}; int owner = -1; int shared = 0; int boundary = 0;

                if (cnt == 3) {                       // corner -- consistent P1 ownership
                    int lc = c[hexCornerIndex(i==P, j==P, k==P)];
                    owner    = cornerOwner[lc];
                    boundary = sharedCorner[lc];      // shared corner -> in the halo (but owner stays P1)
                    key      = DofKey{0, cornerGid[lc], -1,-1,-1, 0};
                } else if (cnt == 2) {                // edge -- provisional, resolve later
                    int vary = onI ? (onJ?2:1) : 0;
                    int t    = (vary==0)?i:(vary==1)?j:k;
                    int bi=(vary==0)?0:(i==P), bj=(vary==1)?0:(j==P), bk=(vary==2)?0:(k==P);
                    int cA = c[hexCornerIndex(bi,bj,bk)];
                    bi=(vary==0)?1:(i==P); bj=(vary==1)?1:(j==P); bk=(vary==2)?1:(k==P);
                    int cB = c[hexCornerIndex(bi,bj,bk)];
                    long gA = cornerGid[cA], gB = cornerGid[cB];
                    long klo = (gA<=gB)?gA:gB, khi = (gA<=gB)?gB:gA;
                    int  pos = (gA<=gB) ? (t-1) : (pm1-t);   // canonical low->high
                    owner  = myRank;
                    shared = (sharedCorner[cA] && sharedCorner[cB]) ? 1 : 0;
                    key    = DofKey{1, klo, khi, -1, -1, pos};
                } else if (cnt == 1) {                // face -- provisional, resolve later
                    int on = onI?0:onJ?1:2;
                    int v1,v2; if(on==0){v1=j;v2=k;} else if(on==1){v1=i;v2=k;} else {v1=i;v2=j;}
                    const int* fc = FACES[on==0?(i==P?4:5):on==1?(j==P?3:2):(k==P?1:0)];
                    int  lcc[4] = { c[fc[0]], c[fc[1]], c[fc[2]], c[fc[3]] };
                    long s[4]   = { cornerGid[lcc[0]], cornerGid[lcc[1]], cornerGid[lcc[2]], cornerGid[lcc[3]] };
                    std::sort(s, s+4);
                    int  pos = hex_face_canonical_pos(P, v1, v2, lcc[0], lcc[1], lcc[2], lcc[3]);
                    owner  = myRank;
                    shared = (sharedCorner[lcc[0]] && sharedCorner[lcc[1]] &&
                              sharedCorner[lcc[2]] && sharedCorner[lcc[3]]) ? 1 : 0;
                    key    = DofKey{2, s[0], s[1], s[2], s[3], pos};
                } else {                              // interior (element-local, no halo)
                    int pos = ((i-1)*pm1 + (j-1))*pm1 + (k-1);
                    owner = elemOwner[e];
                    key   = DofKey{3, (long)e, -1,-1,-1, pos};
                }
                dofOwner[dof]    = owner;
                dofKey[dof]      = key;
                dofShared[dof]   = (uint8_t)shared;            // edge/face -> ownership resolution
                dofBoundary[dof] = (uint8_t)(boundary | shared); // any shared DOF -> halo keys/maps
            }
        }
    }
};

} // namespace fem
} // namespace mars
