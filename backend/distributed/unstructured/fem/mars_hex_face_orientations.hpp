#pragma once

// Canonical face-interior DOF position for tensor-product hex spectral elements.
//
// Two hexes that share a quad face usually parameterize that face in different
// orientations: one of the 8 dihedral symmetries of a square (4 rotations x 2
// reflections). The (p-1)^2 face-INTERIOR DOFs must land in the same global
// slot from both sides, using ONLY the four global corner-node ids of the face
// (no coordinates, no communication). The sorted-4-corner key already gives an
// orientation-invariant face IDENTITY; this function gives the orientation-
// invariant POSITION inside that face.
//
// Why a canonical frame is needed: an element stores a face as the cyclic corner
// order FACES[f] = (c0,c1,c2,c3), with local varying indices v1 along edge
// c0->c1 and v2 along edge c0->c3. The two sides traverse the same physical quad
// in different (possibly mirrored) senses, so their local (v1,v2) frames are
// rotated/reflected relative to each other. Indexing pos by the element's own
// (v1,v2) scatters the SAME physical interior node to DIFFERENT slots -> broken
// C0 continuity at p>=3 (at p=2 there is a single interior node, so no defect).
//
// The fix elects an absolute 2D frame from the global ids alone:
//   origin  = the face corner with the SMALLEST global id,
//   s-axis  = toward its edge-neighbor with the smaller global id,
//   t-axis  = toward its edge-neighbor with the larger global id.
// Because origin/s/t are functions of the four shared global ids only, both
// neighbors pick the SAME physical corner as origin and the SAME physical edges
// as s and t. Each element then maps its own (v1,v2) into (s,t): a reflection at
// the origin (r = p - v when the origin sits at the high end of that local axis)
// plus an axis selection (which of v1/v2 feeds s) absorbs all 8 dihedral cases.
// pos = (s-1)*(p-1) + (t-1) is then identical on both sides.
//
// The mapping is a bijection of the (p-1)^2 interior nodes onto [0,(p-1)^2) for
// every orientation, and yields one common pos per physical node across all 8
// orientations (verified exhaustively for p=2..8). Assumes a conforming mesh
// with globally-unique corner ids, so the argmin origin is tie-free.

namespace mars {
namespace fem {

// Local (x,y) bits of the 4 cyclic face slots in FACES[f] order:
//   slot0=(0,0) slot1=(1,0) slot2=(1,1) slot3=(0,1).
// v1 runs along the x-bit edge (slot0->slot1), v2 along the y-bit edge
// (slot0->slot3); this regular alignment holds for all 6 hex faces.
inline constexpr int kFaceSlotBitX[4] = {0, 1, 1, 0};
inline constexpr int kFaceSlotBitY[4] = {0, 0, 1, 1};

// Canonical interior position for a face node with local varying indices
// (v1,v2), each in [1,p-1], given the face's four global corner ids in cyclic
// FACES[f] order. Returns pos in [0,(p-1)^2). dof = faceBase + fid*(p-1)^2 + pos
// (fid still from the existing sorted-4-corner key, already orientation-invariant).
inline int hex_face_canonical_pos(int p, int v1, int v2,
                                  int g0, int g1, int g2, int g3)
{
    const int pm1 = p - 1;
    const int g[4] = {g0, g1, g2, g3};

    // origin = cyclic slot of the smallest global id (unique for a valid mesh).
    int o = 0;
    for (int a = 1; a < 4; ++a)
        if (g[a] < g[o]) o = a;

    // origin's two cyclic edge-neighbors; the (o+2) slot is the diagonal.
    const int nprev = (o + 3) & 3;
    const int nnext = (o + 1) & 3;

    // s-axis toward the smaller-id neighbor, t-axis toward the larger.
    int sN, tN;
    if (g[nnext] < g[nprev]) { sN = nnext; tN = nprev; }
    else                     { sN = nprev; tN = nnext; }

    // Node distance from the origin corner along the two local axes. When the
    // origin sits at the high end of a local axis (bit==1), that axis is
    // reversed, so the distance is p-v instead of v (this absorbs reflections).
    const int obx = kFaceSlotBitX[o];
    const int oby = kFaceSlotBitY[o];
    const int r1 = (obx == 0) ? v1 : (p - v1);
    const int r2 = (oby == 0) ? v2 : (p - v2);

    // Pick which local axis feeds s and which feeds t: the origin->sN edge moves
    // along v1 iff it changes the x-bit. The complementary edge feeds t. This
    // axis selection swaps under a mirrored cycle, absorbing the remaining
    // dihedral freedom.
    const int s = ((kFaceSlotBitX[sN] - obx) != 0) ? r1 : r2;
    const int t = ((kFaceSlotBitX[tN] - obx) != 0) ? r1 : r2;

    return (s - 1) * pm1 + (t - 1);
}

} // namespace fem
} // namespace mars
