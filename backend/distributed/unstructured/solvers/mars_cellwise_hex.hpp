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

MARS_CW_HD inline int node_at(int x, int y, int z) { return (x * kN + y) * kN + z; }

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

MARS_CW_HD inline int corner_node(int c)
{
    int bx, by, bz;
    corner_bits(c, bx, by, bz);
    return node_at(bx * kP, by * kP, bz * kP);
}

MARS_CW_HD inline void other_axes(int axis, int& o0, int& o1)
{
    o0 = axis == 0 ? 1 : 0;
    o1 = axis == 2 ? 1 : 2;
}

MARS_CW_HD inline int face_node(int f, int i, int j)
{
    const int axis = f / 2;
    int x[3], o0, o1;
    other_axes(axis, o0, o1);
    x[axis] = (f % 2) * kP;
    x[o0] = i;
    x[o1] = j;
    return node_at(x[0], x[1], x[2]);
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

MARS_CW_HD inline int edge_node(int k, int t)
{
    const int axis = k / 4;
    int x[3], o0, o1;
    other_axes(axis, o0, o1);
    x[axis] = t;
    x[o0] = ((k / 2) % 2) * kP;
    x[o1] = (k % 2) * kP;
    return node_at(x[0], x[1], x[2]);
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

MARS_CW_HD inline void face_to_canonical(int code, int i, int j, int& I, int& J)
{
    const int di = (code & 1) ? kP - i : i, dj = (code & 2) ? kP - j : j;
    if (code & 4) { I = di; J = dj; }
    else { I = dj; J = di; }
}

MARS_CW_HD inline void canonical_to_face(int code, int I, int J, int& i, int& j)
{
    const int di = (code & 4) ? I : J, dj = (code & 4) ? J : I;
    i = (code & 1) ? kP - di : di;
    j = (code & 2) ? kP - dj : dj;
}

}  // namespace cellwise
}  // namespace mars
