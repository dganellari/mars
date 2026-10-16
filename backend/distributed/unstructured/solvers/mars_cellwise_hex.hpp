#pragma once
// Local entities of a p = 7 hexahedron in its own frame, shared by every cell-wise
// kernel and table builder. No CUDA or MPI dependency, so it is unit-tested on the host
// against marsir-mlir/test/cellwise_unstructured_ref.py.

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_CW_HD __host__ __device__
#else
#define MARS_CW_HD
#endif

namespace mars {
namespace cellwise {

constexpr int kN = 8, kNN = 64, kN3 = 512;   // nodes per edge, face, element (p = 7)
constexpr int kP = kN - 1;   // polynomial degree (7)

// ---- local entities of a hex in its own frame (VTK corner order) --------------------
// Corner c has signs (x, y, z) with corner 0 at (-,-,-), 1 (+,-,-), 2 (+,+,-), 3 (-,+,-),
// 4-7 the same at z = +. Face f: axis f / 2 at coordinate (f % 2) * p, face coordinates
// (i, j) over the other two axes in increasing order, corner slot s = i-bit + 2 j-bit.
// Edge k: axis k / 4 varies (position t from coordinate 0 to p), the other two axes fixed
// at ((k / 2) % 2) * p and (k % 2) * p, in increasing axis order.
// The node functions take the degree P as a template argument (p = 7 by default) for the
// p-multigrid levels; node (a, b, c) of degree P is (a * (P + 1) + b) * (P + 1) + c.

template <int P = kP>
MARS_CW_HD inline int node_at(int x, int y, int z) { return (x * (P + 1) + y) * (P + 1) + z; }

MARS_CW_HD inline int corner_from_bits(int bx, int by, int bz)
{
    return 4 * bz + (by ? (bx ? 2 : 3) : (bx ? 1 : 0));
}

MARS_CW_HD inline void corner_bits(int c, int& bx, int& by, int& bz)
{
    bz = c >> 2;
    const int q = c & 3;
    bx = (q == 1 || q == 2);
    by = (q >= 2);
}

template <int P = kP>
MARS_CW_HD inline int corner_node(int c)
{
    int bx, by, bz;
    corner_bits(c, bx, by, bz);
    return node_at<P>(bx * P, by * P, bz * P);
}

MARS_CW_HD inline void other_axes(int axis, int& o0, int& o1)
{
    o0 = axis == 0 ? 1 : 0;
    o1 = axis == 2 ? 1 : 2;
}

template <int P = kP>
MARS_CW_HD inline int face_node(int f, int i, int j)
{
    const int axis = f / 2;
    int x[3], o0, o1;
    other_axes(axis, o0, o1);
    x[axis] = (f % 2) * P;
    x[o0] = i;
    x[o1] = j;
    return node_at<P>(x[0], x[1], x[2]);
}

MARS_CW_HD inline int face_corner(int f, int s)
{
    const int axis = f / 2;
    int b[3], o0, o1;
    other_axes(axis, o0, o1);
    b[axis] = f % 2;
    b[o0] = s & 1;
    b[o1] = s >> 1;
    return corner_from_bits(b[0], b[1], b[2]);
}

template <int P = kP>
MARS_CW_HD inline int edge_node(int k, int t)
{
    const int axis = k / 4;
    int x[3], o0, o1;
    other_axes(axis, o0, o1);
    x[axis] = t;
    x[o0] = ((k / 2) % 2) * P;
    x[o1] = (k % 2) * P;
    return node_at<P>(x[0], x[1], x[2]);
}

MARS_CW_HD inline int edge_corner(int k, int end)
{
    const int axis = k / 4;
    int b[3], o0, o1;
    other_axes(axis, o0, o1);
    b[axis] = end;
    b[o0] = (k / 2) % 2;
    b[o1] = k % 2;
    return corner_from_bits(b[0], b[1], b[2]);
}

// The 4 edges of face f (as local edge indices) and its 4 corners.
MARS_CW_HD inline int face_edge(int f, int which)
{
    // Edges of face f run along its two in-face axes, at the face's coordinate on the
    // normal axis and at 0 or p on the other in-face axis.
    const int axis = f / 2, side = f % 2;
    int o0, o1;
    other_axes(axis, o0, o1);
    const int along = which < 2 ? o0 : o1;      // axis the edge varies along
    const int across = which < 2 ? o1 : o0;     // in-face axis it is fixed on
    int fixed[3] = {0, 0, 0};
    fixed[axis] = side;
    fixed[across] = which % 2;
    int a0, a1;
    other_axes(along, a0, a1);
    return along * 4 + fixed[a0] * 2 + fixed[a1];
}

// Canonical face frame from the global keys at the 4 corner slots: origin at the
// smallest key, first axis toward the smaller of its two neighbours (the convention of
// hexFaceCanonicalPosDev). Code = origin slot | (first axis is i) << 2.
MARS_CW_HD inline int face_frame_code(const unsigned long long (&key)[4])
{
    int o = 0;
    for (int s = 1; s < 4; ++s)
        if (key[s] < key[o]) o = s;
    const bool first_is_i = key[o ^ 1] < key[o ^ 2];
    return o | (first_is_i ? 4 : 0);
}

template <int P = kP>
MARS_CW_HD inline void face_to_canonical(int code, int i, int j, int& I, int& J)
{
    const int di = (code & 1) ? P - i : i, dj = (code & 2) ? P - j : j;
    if (code & 4) { I = di; J = dj; }
    else { I = dj; J = di; }
}

template <int P = kP>
MARS_CW_HD inline void canonical_to_face(int code, int I, int J, int& i, int& j)
{
    const int di = (code & 4) ? I : J, dj = (code & 4) ? J : I;
    i = (code & 1) ? P - di : di;
    j = (code & 2) ? P - dj : dj;
}

// Entity bits of an element's 32-bit masks (counting copies, Dirichlet entities): face f
// at bit f, edge k at kEdgeBit + k, corner c at kVertexBit + c, interior nodes at
// kInteriorBit. The gather kernel loads each entity's descriptor on the lane with the
// same number, so one ballot gives the Dirichlet mask.
constexpr int kEdgeBit = 8, kVertexBit = 20, kInteriorBit = 31;

// The local entity that holds node (a, b, c): kind 0 interior, 1 face, 2 edge, 3 corner,
// and its local number (f, k or c).
struct LocalEntity {
    int kind, index;
};

template <int P = kP>
MARS_CW_HD inline LocalEntity local_entity(int a, int b, int c)
{
    const int x[3] = {a, b, c};
    int on = 0;
    for (int k = 0; k < 3; ++k) on += (x[k] == 0 || x[k] == P);
    if (on == 0) return {0, 0};
    if (on == 3) return {3, corner_from_bits(x[0] == P, x[1] == P, x[2] == P)};
    int axis = 0, o0, o1;
    if (on == 1) {
        while (x[axis] != 0 && x[axis] != P) ++axis;
        return {1, axis * 2 + (x[axis] == P)};
    }
    while (x[axis] == 0 || x[axis] == P) ++axis;   // the edge runs along the free axis
    other_axes(axis, o0, o1);
    return {2, axis * 4 + (x[o0] == P) * 2 + (x[o1] == P)};
}

MARS_CW_HD inline int entity_bit(const LocalEntity& le)
{
    if (le.kind == 0) return kInteriorBit;
    return (le.kind == 1 ? 0 : (le.kind == 2 ? kEdgeBit : kVertexBit)) + le.index;
}

}  // namespace cellwise
}  // namespace mars
