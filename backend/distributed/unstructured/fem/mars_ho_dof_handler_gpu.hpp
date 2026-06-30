#pragma once

// GPU-native build of the SAME high-order distributed DOF numbering as
// HODofHandler::buildDistributed (mars_ho_dof_handler.hpp). Flag-gated, optional
// path; the host buildDistributed stays the default and is left untouched.
//
// WHY: at ~10M elements/rank the host numbering spends minutes in std::map edge/
// face dedup over ~180M (12 edge + 6 face) entries/rank. The dedup is a sort +
// unique + binary search -- a textbook thrust pipeline -- so we move it (and the
// per-DOF tagging) onto the GPU.
//
// WHAT WE REPRODUCE EXACTLY (must match the host convention bit-for-bit):
//   - the four DOF kinds (corner cnt==3, edge cnt==2, face cnt==1, interior cnt==0)
//     and the (i,j,k)->l index l = i*n*n + j*n + k,
//   - the hex EDGES[12] / FACES[6] tables and hexCornerIndex(),
//   - the edge pos = canonical low->high (t-1 or pm1-t by global-id order),
//   - the face pos via hex_face_canonical_pos() (orientation-invariant),
//   - the interior pos = ((i-1)*pm1 + (j-1))*pm1 + (k-1),
//   - the canonical DofKey {kind, sorted global corner ids, pos},
//   - dofOwner/dofShared/dofBoundary rules from buildDistributed.
//
// WHAT MAY LEGITIMATELY DIFFER FROM HOST: the LOCAL edge/face ids (host assigns
// them in std::map insertion order; we assign them in sorted-key order). That
// permutes the local DOF ids inside elemDof/dofOwner/dofKey/... but the cross-rank
// identity is the DofKey, never the local id. So the correct equivalence is: same
// numDof/nEdge/nFace, and same MULTISET of DofKeys. The host driver self-check
// (--self-check) verifies exactly this, plus the end-to-end A.1 apply gate.
//
// MEMORY (see writeup): the transient sort buffers dominate setup HBM. Per owned
// element we touch 18 entity keys (12 edge + 6 face) * 16B for the emit arrays,
// plus the per-local-dof arrays (~8 longs + ints per DOF). These compete with the
// operator metric d_G for HBM, so this path has a LOWER per-GPU DOF ceiling than
// host numbering. All transients are freed before returning so the apply's d_G
// allocation sees them gone.
//
// Guarded by USE_CUDA; only the .cu driver includes this. Pure-host translation
// units never see the device code.

#include "mars_ho_dof_handler.hpp"
#include "mars_hex_face_orientations.hpp"

#ifdef USE_CUDA

#include <cstdint>
#include <vector>
#include <utility>

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/binary_search.h>
#include <thrust/copy.h>
#include <thrust/transform.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/for_each.h>

namespace mars {
namespace fem {

// Packed sorted-corner key for an edge (2 local corners) or face (4 local corners).
// Local corner indices fit in 32 bits (node count per rank << 2^31). We pack into
// two 64-bit words so a single thrust::sort over the pair sorts entities by their
// full sorted-corner tuple -- identical entities become adjacent for unique().
//   edge: hi = (c0<<32)|c1 (c0<=c1), lo = 0
//   face: hi = (c0<<32)|c1, lo = (c2<<32)|c3 (c0<=c1<=c2<=c3)
struct HoEntityKey {
    uint64_t hi;
    uint64_t lo;
    __host__ __device__ bool operator<(const HoEntityKey& o) const {
        return hi < o.hi || (hi == o.hi && lo < o.lo);
    }
    __host__ __device__ bool operator==(const HoEntityKey& o) const {
        return hi == o.hi && lo == o.lo;
    }
};

namespace ho_gpu_detail {

// Device-side hex topology, mirroring the host EDGES/FACES tables exactly.
// Returned from accessor functions (not __constant__ globals) so the header is
// ODR-safe if ever included in more than one translation unit.
__host__ __device__ inline int edgeCorner(int ed, int end) {
    const int E[12][2] = {
        {0,1},{1,2},{2,3},{3,0}, {4,5},{5,6},{6,7},{7,4}, {0,4},{1,5},{2,6},{3,7} };
    return E[ed][end];
}
__host__ __device__ inline int faceCorner(int f, int slot) {
    const int F[6][4] = {
        {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, {3,2,6,7}, {1,2,6,5}, {0,3,7,4} };
    return F[f][slot];
}

__host__ __device__ inline int hexCornerIndexDev(int bi, int bj, int bk) {
    // Same table as host hexCornerIndex (tbl[bk][bj][bi]).
    const int tbl[2][2][2] = { { {0,1},{3,2} }, { {4,5},{7,6} } };
    return tbl[bk][bj][bi];
}

// Device copy of hex_face_canonical_pos. The included host version is plain host
// inline; the numbering kernel needs a __device__ form. Logic is line-for-line
// identical so positions match the host path exactly.
__host__ __device__ inline int hexFaceCanonicalPosDev(int p, int v1, int v2,
                                                      long g0, long g1, long g2, long g3) {
    const int pm1 = p - 1;
    const long g[4] = {g0, g1, g2, g3};
    int o = 0;
    for (int a = 1; a < 4; ++a) if (g[a] < g[o]) o = a;
    const int nprev = (o + 3) & 3;
    const int nnext = (o + 1) & 3;
    int sN, tN;
    if (g[nnext] < g[nprev]) { sN = nnext; tN = nprev; }
    else                     { sN = nprev; tN = nnext; }
    const int bitX[4] = {0, 1, 1, 0};
    const int bitY[4] = {0, 0, 1, 1};
    const int obx = bitX[o];
    const int oby = bitY[o];
    const int r1 = (obx == 0) ? v1 : (p - v1);
    const int r2 = (oby == 0) ? v2 : (p - v2);
    const int s = ((bitX[sN] - obx) != 0) ? r1 : r2;
    const int t = ((bitX[tN] - obx) != 0) ? r1 : r2;
    return (s - 1) * pm1 + (t - 1);
}

// Edge key as a single uint64 (2 local corners, sorted) -> thrust::sort takes the
// RADIX path (not the custom-struct merge_sort that crashed at ~116M keys).
__host__ __device__ inline uint64_t edgeKeyU64(int a, int b) {
    int c0 = a, c1 = b;
    if (c0 > c1) { int t = c0; c0 = c1; c1 = t; }
    return ((uint64_t)(uint32_t)c0 << 32) | (uint32_t)c1;
}

__host__ __device__ inline HoEntityKey faceKey(int a, int b, int c, int d) {
    int s[4] = {a, b, c, d};
    // insertion sort of 4 -- tiny and branch-cheap on device
    for (int i = 1; i < 4; ++i) {
        int v = s[i], j = i - 1;
        while (j >= 0 && s[j] > v) { s[j+1] = s[j]; --j; }
        s[j+1] = v;
    }
    HoEntityKey k;
    k.hi = ((uint64_t)(uint32_t)s[0] << 32) | (uint32_t)s[1];
    k.lo = ((uint64_t)(uint32_t)s[2] << 32) | (uint32_t)s[3];
    return k;
}

} // namespace ho_gpu_detail

// Per-DOF device columns kept resident for the ALL-GPU resolveHoDofOwnershipGpu (the
// resolve's exact inputs/output). Normally buildDistributedGpu frees these at scope exit;
// when the caller passes a non-null HoOwnershipDeviceData* it MOVES them out instead so the
// device-resident resolve can run with no host round-trip. d_dofOwner is mutated in place
// by the resolve; copy it back to dof.dofOwner afterwards.
struct HoOwnershipDeviceData {
    long numDof = 0;
    thrust::device_vector<long> dofKind, dofG0, dofG1, dofG2, dofG3;
    thrust::device_vector<int>  dofPos, dofShared, dofOwner;
    // dofBoundary kept too so the device halo build (HoHalo::buildDevice) keys only the
    // boundary DOF straight from device -- same columns, no re-upload. uint8 lane.
    thrust::device_vector<uint8_t> dofBoundary;
    // elemDof [nElem*N3] kept DEVICE-resident so the apply consumes it directly -- no
    // host download here + no H2D re-upload there. The apply needs a d_elemDof of exactly
    // this size regardless, so keeping it is not extra peak HBM. The caller MOVES it out
    // before freeing the rest of this struct (the key columns are huge; elemDof is not).
    thrust::device_vector<int>  elemDof;
};

// GPU build of the distributed HO numbering. Same inputs as buildDistributed;
// writes the SAME host members (elemDof, dofOwner, dofKey, dofShared, dofBoundary,
// numDof, nEdge, nFace, ...) by pulling the device results back, so the downstream
// host resolveHoDofOwnership + HoHalo::build work unchanged.
//
// Kept as a free function taking the handler by reference (not a member) to avoid
// touching the host class -- per the constraint "add alongside, don't change".
//
// If keepOwn != null, the per-DOF device columns the resolve needs are MOVED into it
// (kept resident) for the all-GPU resolveHoDofOwnershipGpu; otherwise they free here.
//
// CORE: all numbering inputs are ALREADY device-resident here (the SoA int corner
// buffer d_corners_ptr[c*nElem+e], and per-corner gid/owner/shared/elemOwner pointers).
// The two public wrappers below feed it either from host std::vectors (host-input path,
// for fallback/self-check) or straight from the domain's device columns (device-input
// path, zero H2D). Keeping ONE numbering body guarantees both paths number bit-for-bit
// identically -> the A.1 oracle is shared.
inline void buildDistributedGpuCore(HODofHandler&  dof,
                                    const int*     d_corners_ptr,   // SoA [c*nElem+e], LOCAL ids
                                    long           nElem,
                                    long           numCornerNodes,
                                    int            order,
                                    const long*    d_gid_ptr,       // per LOCAL corner id
                                    const int*     d_cowner_ptr,
                                    const uint8_t* d_scorner_ptr,
                                    const int*     d_eowner_ptr,    // per OWNED element
                                    int            myRank,
                                    HoOwnershipDeviceData* keepOwn)
{
    using namespace ho_gpu_detail;

    const int P   = order;
    const int n   = P + 1;
    const int N3  = n * n * n;
    const int pm1 = P - 1;
    const long nCorner = numCornerNodes;

    dof.P = P; dof.n = n; dof.N3 = N3;
    dof.nElem = nElem; dof.nCorner = nCorner;

    // STAGES 1+2 use RADIX-friendly uint64 sorts, NOT a custom 16-byte struct sort.
    // The original thrust::sort over HoEntityKey took cub's MERGE_SORT path and hit
    // cudaErrorIllegalAddress at ~116M keys; uint64 keys take the robust RADIX path.
    // Edge = 2 corners -> one uint64. Face = 4 corners (>64 bits) -> two uint64
    // columns sorted by a STABLE 2-pass radix (by lo, then by hi). unique()/lower_bound
    // are compactions / per-thread searches, not parallel sorts, so they keep the
    // struct form for faces and are unaffected.
    long nEdge = 0, nFace = 0;
    thrust::device_vector<uint64_t>    d_uEdge;   // unique sorted edge keys (uint64)
    thrust::device_vector<HoEntityKey> d_uFace;   // unique sorted face keys (struct, for lower_bound)

    // ---- STAGE 1: edges -> uint64 radix sort + unique ----
    {
        thrust::device_vector<uint64_t> d_edgeAll((size_t)nElem * 12);
        uint64_t* p = thrust::raw_pointer_cast(d_edgeAll.data());
        thrust::for_each(
            thrust::counting_iterator<long>(0),
            thrust::counting_iterator<long>(nElem),
            [d_corners_ptr, p, nElem] __device__ (long e) {
                int c[8];
                #pragma unroll
                for (int x = 0; x < 8; ++x) c[x] = d_corners_ptr[(long)x * nElem + e];
                #pragma unroll
                for (int ed = 0; ed < 12; ++ed)
                    p[e * 12 + ed] = edgeKeyU64(c[edgeCorner(ed,0)], c[edgeCorner(ed,1)]);
            });
        thrust::sort(d_edgeAll.begin(), d_edgeAll.end());   // RADIX (uint64)
        d_uEdge.resize(d_edgeAll.size());
        auto uend = thrust::unique_copy(d_edgeAll.begin(), d_edgeAll.end(), d_uEdge.begin());
        nEdge = (long)(uend - d_uEdge.begin());
        d_uEdge.resize(nEdge);
    } // d_edgeAll freed here

    // ---- STAGE 2: faces -> two uint64 columns, stable 2-pass radix, zip-unique ----
    {
        const size_t Nf = (size_t)nElem * 6;
        thrust::device_vector<uint64_t> d_hi(Nf), d_lo(Nf);
        uint64_t* phi = thrust::raw_pointer_cast(d_hi.data());
        uint64_t* plo = thrust::raw_pointer_cast(d_lo.data());
        thrust::for_each(
            thrust::counting_iterator<long>(0),
            thrust::counting_iterator<long>(nElem),
            [d_corners_ptr, phi, plo, nElem] __device__ (long e) {
                int c[8];
                #pragma unroll
                for (int x = 0; x < 8; ++x) c[x] = d_corners_ptr[(long)x * nElem + e];
                #pragma unroll
                for (int f = 0; f < 6; ++f) {
                    HoEntityKey k = faceKey(c[faceCorner(f,0)], c[faceCorner(f,1)],
                                            c[faceCorner(f,2)], c[faceCorner(f,3)]);
                    phi[e * 6 + f] = k.hi;
                    plo[e * 6 + f] = k.lo;
                }
            });
        // sort by (hi,lo): radix on lo, then STABLE radix on hi (equal-hi keep lo order)
        thrust::stable_sort_by_key(d_lo.begin(), d_lo.end(), d_hi.begin());
        thrust::stable_sort_by_key(d_hi.begin(), d_hi.end(), d_lo.begin());
        // dedup the sorted (hi,lo) pairs via a zip iterator (tuple==, a compaction not a sort)
        auto zbeg = thrust::make_zip_iterator(thrust::make_tuple(d_hi.begin(), d_lo.begin()));
        auto zend = thrust::make_zip_iterator(thrust::make_tuple(d_hi.end(),   d_lo.end()));
        auto znew = thrust::unique(zbeg, zend);
        nFace = (long)(znew - zbeg);
        // pack the unique columns into HoEntityKey for the stage-3 per-thread lower_bound
        d_uFace.resize(nFace);
        HoEntityKey*    puf  = thrust::raw_pointer_cast(d_uFace.data());
        const uint64_t* cphi = thrust::raw_pointer_cast(d_hi.data());
        const uint64_t* cplo = thrust::raw_pointer_cast(d_lo.data());
        thrust::for_each(
            thrust::counting_iterator<long>(0),
            thrust::counting_iterator<long>(nFace),
            [puf, cphi, cplo] __device__ (long i) { puf[i].hi = cphi[i]; puf[i].lo = cplo[i]; });
    } // d_hi/d_lo freed here

    dof.nEdge = nEdge;
    dof.nFace = nFace;

    const long edgeBase = nCorner;
    const long faceBase = edgeBase + nEdge * pm1;
    const long intrBase = faceBase + nFace * (long)pm1 * pm1;
    const long numDof   = intrBase + nElem * (long)pm1 * pm1 * pm1;
    dof.numDof = numDof;

    const uint64_t*    d_uEdge_ptr = thrust::raw_pointer_cast(d_uEdge.data());
    const HoEntityKey* d_uFace_ptr = thrust::raw_pointer_cast(d_uFace.data());

    // ================================================================
    // STAGE 3: numbering kernel.
    // For each (e, i, j, k) -> l, classify the node and compute its local DOF id
    // (corner = corner local id; edge/face id via lower_bound over the unique key
    // arrays; interior = element block). Write elemDof and SCATTER the per-DOF
    // owner/key/shared into the dof-indexed arrays. boundary is an OR across all
    // elements touching the DOF (atomicOr on an int lane).
    // ================================================================
    thrust::device_vector<int> d_elemDof((size_t)nElem * N3, -1);
    int* d_elemDof_ptr = thrust::raw_pointer_cast(d_elemDof.data());

    // per-DOF outputs (device), pulled back to host members at the end.
    thrust::device_vector<int>  d_dofOwner(numDof, -1);
    thrust::device_vector<long> d_dofKind(numDof, -1);
    thrust::device_vector<long> d_dofG0(numDof, -1);
    thrust::device_vector<long> d_dofG1(numDof, -1);
    thrust::device_vector<long> d_dofG2(numDof, -1);
    thrust::device_vector<long> d_dofG3(numDof, -1);
    thrust::device_vector<int>  d_dofPos(numDof, -1);
    thrust::device_vector<int>  d_dofShared(numDof, 0);   // int lane for atomicOr
    thrust::device_vector<int>  d_dofBoundary(numDof, 0); // int lane for atomicOr

    int*  p_dofOwner    = thrust::raw_pointer_cast(d_dofOwner.data());
    long* p_dofKind     = thrust::raw_pointer_cast(d_dofKind.data());
    long* p_dofG0       = thrust::raw_pointer_cast(d_dofG0.data());
    long* p_dofG1       = thrust::raw_pointer_cast(d_dofG1.data());
    long* p_dofG2       = thrust::raw_pointer_cast(d_dofG2.data());
    long* p_dofG3       = thrust::raw_pointer_cast(d_dofG3.data());
    int*  p_dofPos      = thrust::raw_pointer_cast(d_dofPos.data());
    int*  p_dofShared   = thrust::raw_pointer_cast(d_dofShared.data());
    int*  p_dofBoundary = thrust::raw_pointer_cast(d_dofBoundary.data());

    thrust::for_each(
        thrust::counting_iterator<long>(0),
        thrust::counting_iterator<long>(nElem),
        [=] __device__ (long e) {
            int c[8];
            #pragma unroll
            for (int x = 0; x < 8; ++x) c[x] = d_corners_ptr[(long)x * nElem + e];

            for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
            {
                int l   = i*n*n + j*n + k;
                int onI = (i==0||i==P), onJ=(j==0||j==P), onK=(k==0||k==P);
                int cnt = onI + onJ + onK;

                long dofId = -1;
                int  kind = -1, owner = -1, shared = 0, boundary = 0, pos = 0;
                long gg0 = -1, gg1 = -1, gg2 = -1, gg3 = -1;

                if (cnt == 3) {                       // corner -- P1 ownership
                    int lc = c[hexCornerIndexDev(i==P, j==P, k==P)];
                    dofId    = lc;                    // corner DOF id = corner local id
                    kind     = 0;
                    owner    = d_cowner_ptr[lc];
                    boundary = d_scorner_ptr[lc];
                    gg0      = d_gid_ptr[lc];
                    pos      = 0;
                } else if (cnt == 2) {                // edge -- provisional, resolve later
                    int vary = onI ? (onJ?2:1) : 0;
                    int t    = (vary==0)?i:(vary==1)?j:k;
                    int bi=(vary==0)?0:(i==P), bj=(vary==1)?0:(j==P), bk=(vary==2)?0:(k==P);
                    int cA = c[hexCornerIndexDev(bi,bj,bk)];
                    bi=(vary==0)?1:(i==P); bj=(vary==1)?1:(j==P); bk=(vary==2)?1:(k==P);
                    int cB = c[hexCornerIndexDev(bi,bj,bk)];
                    // local edge id by binary search over unique edge keys (keyed on
                    // sorted LOCAL corners, just like the host std::map<array<int,2>>)
                    uint64_t ek = edgeKeyU64(cA, cB);
                    long eid = thrust::lower_bound(thrust::seq, d_uEdge_ptr, d_uEdge_ptr + nEdge, ek)
                               - d_uEdge_ptr;
                    long gA = d_gid_ptr[cA], gB = d_gid_ptr[cB];
                    long klo = (gA<=gB)?gA:gB, khi = (gA<=gB)?gB:gA;
                    int  posCanon = (gA<=gB) ? (t-1) : (pm1-t);   // canonical low->high
                    dofId    = edgeBase + eid * pm1 + posCanon;
                    kind     = 1;
                    owner    = myRank;
                    shared   = (d_scorner_ptr[cA] && d_scorner_ptr[cB]) ? 1 : 0;
                    gg0 = klo; gg1 = khi; pos = posCanon;
                } else if (cnt == 1) {                // face -- provisional, resolve later
                    int on = onI?0:onJ?1:2;
                    int v1,v2; if(on==0){v1=j;v2=k;} else if(on==1){v1=i;v2=k;} else {v1=i;v2=j;}
                    int faceSel = on==0?(i==P?4:5):on==1?(j==P?3:2):(k==P?1:0);
                    int lcc0 = c[faceCorner(faceSel,0)], lcc1 = c[faceCorner(faceSel,1)];
                    int lcc2 = c[faceCorner(faceSel,2)], lcc3 = c[faceCorner(faceSel,3)];
                    // local face id by binary search over unique face keys (keyed on
                    // sorted LOCAL corners, just like the host std::map<array<int,4>>)
                    HoEntityKey fk = faceKey(lcc0, lcc1, lcc2, lcc3);
                    long fid = thrust::lower_bound(thrust::seq, d_uFace_ptr, d_uFace_ptr + nFace, fk)
                               - d_uFace_ptr;
                    long s0 = d_gid_ptr[lcc0], s1 = d_gid_ptr[lcc1];
                    long s2 = d_gid_ptr[lcc2], s3 = d_gid_ptr[lcc3];
                    // sort the 4 GLOBAL ids for the key (matches host std::sort(s,s+4))
                    long sg[4] = {s0,s1,s2,s3};
                    for (int a=1;a<4;++a){ long v=sg[a]; int b=a-1; while(b>=0&&sg[b]>v){sg[b+1]=sg[b];--b;} sg[b+1]=v; }
                    // canonical pos from LOCAL corner indices -- EXACTLY as host L256
                    // (host passes lcc[] = local ids, not global). Local ids are
                    // SFC-key-monotonic so the argmin/ordering matches the global
                    // frame; passing locals reproduces the host bit-for-bit and does
                    // not rely on that monotonicity holding.
                    int posCanon = hexFaceCanonicalPosDev(P, v1, v2, lcc0, lcc1, lcc2, lcc3);
                    dofId   = faceBase + fid * (long)pm1 * pm1 + posCanon;
                    kind    = 2;
                    owner   = myRank;
                    shared  = (d_scorner_ptr[lcc0] && d_scorner_ptr[lcc1] &&
                               d_scorner_ptr[lcc2] && d_scorner_ptr[lcc3]) ? 1 : 0;
                    gg0=sg[0]; gg1=sg[1]; gg2=sg[2]; gg3=sg[3]; pos = posCanon;
                } else {                              // interior -- element-local
                    int posCanon = ((i-1)*pm1 + (j-1))*pm1 + (k-1);
                    dofId   = intrBase + e * (long)pm1 * pm1 * pm1 + posCanon;
                    kind    = 3;
                    owner   = d_eowner_ptr[e];
                    gg0 = e; pos = posCanon;
                }

                d_elemDof_ptr[e * N3 + l] = (int)dofId;

                // Scatter the per-DOF identity. Many elements write the same dof;
                // owner/kind/key/pos are identical from every writer (functions of
                // global ids only), so a plain store is race-safe (same value). The
                // share/boundary flags are an OR across writers -> atomicOr.
                p_dofOwner[dofId] = owner;
                p_dofKind[dofId]  = kind;
                p_dofG0[dofId]    = gg0;
                p_dofG1[dofId]    = gg1;
                p_dofG2[dofId]    = gg2;
                p_dofG3[dofId]    = gg3;
                p_dofPos[dofId]   = pos;
                // dofBoundary = boundary | shared (host L269): corners carry the share
                // flag in 'boundary'; edge/face carry it in 'shared'.
                if (shared)   atomicOr(&p_dofShared[dofId], 1);
                if (boundary || shared) atomicOr(&p_dofBoundary[dofId], 1);
            }
        });

    // ---- pull device results into the host members ----
    // elemDof: download to the host member ONLY when no device handoff is requested
    // (host fallback / self-check upload it via H2D). With keepOwn the apply reads the
    // device copy directly -> skip the download (and the later H2D) entirely.
    if (!keepOwn) {
        dof.elemDof.resize((size_t)nElem * N3);
        thrust::copy(d_elemDof.begin(), d_elemDof.end(), dof.elemDof.begin());
    }

    dof.dofOwner.resize(numDof);
    thrust::copy(d_dofOwner.begin(), d_dofOwner.end(), dof.dofOwner.begin());

    // dofKey/dofShared/dofBoundary host members feed ONLY the host resolve + host halo build.
    // The device path (keepOwn) runs the resolve + buildDevice off the device columns moved
    // into keepOwn below, so skip this per-DOF D2H + host rebuild entirely -- it was dead work
    // (8 per-DOF copies + a numDof host loop) on the full-GPU path.
    if (!keepOwn) {
        std::vector<long> h_kind(numDof), h_g0(numDof), h_g1(numDof), h_g2(numDof), h_g3(numDof);
        std::vector<int>  h_pos(numDof), h_shared(numDof), h_boundary(numDof);
        thrust::copy(d_dofKind.begin(),     d_dofKind.end(),     h_kind.begin());
        thrust::copy(d_dofG0.begin(),       d_dofG0.end(),       h_g0.begin());
        thrust::copy(d_dofG1.begin(),       d_dofG1.end(),       h_g1.begin());
        thrust::copy(d_dofG2.begin(),       d_dofG2.end(),       h_g2.begin());
        thrust::copy(d_dofG3.begin(),       d_dofG3.end(),       h_g3.begin());
        thrust::copy(d_dofPos.begin(),      d_dofPos.end(),      h_pos.begin());
        thrust::copy(d_dofShared.begin(),   d_dofShared.end(),   h_shared.begin());
        thrust::copy(d_dofBoundary.begin(), d_dofBoundary.end(), h_boundary.begin());

        dof.dofKey.resize(numDof);
        dof.dofShared.resize(numDof);
        dof.dofBoundary.resize(numDof);
        for (long d = 0; d < numDof; ++d) {
            dof.dofKey[d] = HODofHandler::DofKey{ (int)h_kind[d], h_g0[d], h_g1[d],
                                                  h_g2[d], h_g3[d], h_pos[d] };
            dof.dofShared[d]   = (uint8_t)(h_shared[d]   ? 1 : 0);
            dof.dofBoundary[d] = (uint8_t)(h_boundary[d] ? 1 : 0);
        }
    }

    // Hand the resolve's device columns to the caller (move, no copy) so the all-GPU
    // resolve runs without re-uploading them. They stay alive in keepOwn; everything else
    // below frees as usual. d_dofShared/d_dofOwner are int lanes -- the resolve's inputs.
    if (keepOwn) {
        keepOwn->numDof    = numDof;
        keepOwn->elemDof   = std::move(d_elemDof);   // device-resident for the apply (no H2D)
        keepOwn->dofKind   = std::move(d_dofKind);
        keepOwn->dofG0     = std::move(d_dofG0);
        keepOwn->dofG1     = std::move(d_dofG1);
        keepOwn->dofG2     = std::move(d_dofG2);
        keepOwn->dofG3     = std::move(d_dofG3);
        keepOwn->dofPos    = std::move(d_dofPos);
        keepOwn->dofShared = std::move(d_dofShared);
        keepOwn->dofOwner  = std::move(d_dofOwner);
        // d_dofBoundary is an int atomicOr lane; the halo build wants a uint8 mask -> shrink
        // it here (one transform, ~numDof bytes) so HoHalo::buildDevice keys boundary DOF
        // from device with no host round-trip.
        keepOwn->dofBoundary.resize(numDof);
        thrust::transform(d_dofBoundary.begin(), d_dofBoundary.end(), keepOwn->dofBoundary.begin(),
                          [] __device__ (int v) { return (uint8_t)(v ? 1 : 0); });
    }
    // All device transients (d_uEdge/d_uFace/d_elemDof/d_dof*) destruct here as the
    // thrust::device_vectors leave scope -> HBM freed before the caller's apply
    // allocates d_G. (Those moved into keepOwn outlive this scope by design.)
}

// HOST-INPUT wrapper: same signature/behaviour as before. Uploads the host corner/gid/
// owner/shared/elemOwner vectors, packs the SoA int corner buffer, then runs the core.
// Kept for MARS_HO_HOST + the self-check and any host-vector caller. UNCHANGED numbering.
inline void buildDistributedGpu(HODofHandler&                         dof,
                                const std::vector<std::array<int,8>>& elemCorners,
                                long numCornerNodes, int order,
                                const std::vector<long>&     cornerGid,
                                const std::vector<int>&      cornerOwner,
                                const std::vector<int>&      elemOwner,
                                int                          myRank,
                                const std::vector<uint8_t>&  sharedCorner,
                                HoOwnershipDeviceData*       keepOwn = nullptr)
{
    const long nElem = (long)elemCorners.size();
    // SoA layout d_corners[c*nElem + e] keeps coalesced per-corner reads (same as before).
    std::vector<int> h_corners((size_t)nElem * 8);
    for (long e = 0; e < nElem; ++e)
        for (int c = 0; c < 8; ++c)
            h_corners[(size_t)c * nElem + e] = elemCorners[e][c];
    thrust::device_vector<int>     d_corners(h_corners.begin(), h_corners.end());
    thrust::device_vector<long>    d_cornerGid(cornerGid.begin(), cornerGid.end());
    thrust::device_vector<int>     d_cornerOwner(cornerOwner.begin(), cornerOwner.end());
    thrust::device_vector<uint8_t> d_sharedCorner(sharedCorner.begin(), sharedCorner.end());
    thrust::device_vector<int>     d_elemOwner(elemOwner.begin(), elemOwner.end());
    buildDistributedGpuCore(dof, thrust::raw_pointer_cast(d_corners.data()), nElem,
                            numCornerNodes, order,
                            thrust::raw_pointer_cast(d_cornerGid.data()),
                            thrust::raw_pointer_cast(d_cornerOwner.data()),
                            thrust::raw_pointer_cast(d_sharedCorner.data()),
                            thrust::raw_pointer_cast(d_elemOwner.data()),
                            myRank, keepOwn);
}

// Pack the SoA int corner buffer [c*nElem+e] from the domain's 8 KeyType connectivity
// columns (full element range) over the OWNED slice [startE, startE+nElem). LOCAL node
// ids fit in 32 bits, so the KeyType->int narrowing is exact. This replaces the host
// h_corners pack + its H2D on the device-input path.
template<typename KeyType>
__global__ void ho_pack_corner_soa_kernel(
    const KeyType* __restrict__ c0, const KeyType* __restrict__ c1,
    const KeyType* __restrict__ c2, const KeyType* __restrict__ c3,
    const KeyType* __restrict__ c4, const KeyType* __restrict__ c5,
    const KeyType* __restrict__ c6, const KeyType* __restrict__ c7,
    int* __restrict__ d_corners, long startElem, long nElem)
{
    const long e = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nElem) return;
    const long ge = startElem + e;
    d_corners[0L * nElem + e] = (int)c0[ge];
    d_corners[1L * nElem + e] = (int)c1[ge];
    d_corners[2L * nElem + e] = (int)c2[ge];
    d_corners[3L * nElem + e] = (int)c3[ge];
    d_corners[4L * nElem + e] = (int)c4[ge];
    d_corners[5L * nElem + e] = (int)c5[ge];
    d_corners[6L * nElem + e] = (int)c6[ge];
    d_corners[7L * nElem + e] = (int)c7[ge];
}

// DEVICE-INPUT wrapper: every numbering input is already device-resident (the domain's
// connectivity columns + the per-corner gid/owner/shared columns built by device kernels
// in extractDistDof). Packs the SoA int corner buffer ON DEVICE from the connectivity,
// then runs the SAME core -> identical numbering, zero H2D. conn columns are over the
// FULL element range; the owned slice starts at startElem.
template<typename KeyType>
inline void buildDistributedGpuDevice(HODofHandler&  dof,
                                      const KeyType* d_conn[8],
                                      long           startElem,
                                      long           nElem,
                                      long           numCornerNodes,
                                      int            order,
                                      const long*    d_cornerGid,    // per LOCAL corner id
                                      const int*     d_cornerOwner,
                                      const uint8_t* d_sharedCorner,
                                      const int*     d_elemOwner,    // per OWNED element
                                      int            myRank,
                                      HoOwnershipDeviceData* keepOwn = nullptr)
{
    thrust::device_vector<int> d_corners((size_t)nElem * 8);
    const int block = 256;
    const long grid = (nElem + block - 1) / block;
    if (grid > 0)
        ho_pack_corner_soa_kernel<KeyType><<<(unsigned)grid, block>>>(
            d_conn[0], d_conn[1], d_conn[2], d_conn[3],
            d_conn[4], d_conn[5], d_conn[6], d_conn[7],
            thrust::raw_pointer_cast(d_corners.data()), startElem, nElem);
    buildDistributedGpuCore(dof, thrust::raw_pointer_cast(d_corners.data()), nElem,
                            numCornerNodes, order, d_cornerGid, d_cornerOwner,
                            d_sharedCorner, d_elemOwner, myRank, keepOwn);
}

} // namespace fem
} // namespace mars

#endif // USE_CUDA
