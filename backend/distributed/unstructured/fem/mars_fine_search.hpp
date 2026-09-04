#ifndef MARS_FINE_SEARCH_HPP
#define MARS_FINE_SEARCH_HPP

// Fine search: project a point onto a P1 triangle and return where it landed.
//
// The second half of the interface pairing. coarseSearch narrows a slave integration point down to
// a handful of candidate master faces by box overlap; this decides WHICH of them actually owns the
// point and WHERE on it, which is what the DG-IP flux needs to interpolate the opposing state.
//
// Closed form, no iteration. A P1 tet face is a planar triangle, so the projection is exact: drop
// the point onto the plane, solve the 2x2 for the barycentric coordinates, and if it lands outside
// the triangle clamp to the nearest edge or vertex. That last part matters -- on a curved or
// mismatched interface the exact projection often falls just outside every candidate, and without
// clamping such a point silently finds no donor.
//
// Device-first, per the MARS convention: the plain names are the device ones, `Host` marks the
// reference. The projection itself is a __host__ __device__ free function, so both paths run
// literally the same arithmetic and the gate compares them exactly.

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_FS_CUDA 1
#include <thrust/device_vector.h>
#define MARS_FS_FN __host__ __device__
#else
#define MARS_FS_FN
#endif

#include <cmath>
#include <cstdint>

namespace mars
{
namespace fem
{

template<typename T>
struct Vec3
{
    T x, y, z;
};

template<typename T>
MARS_FS_FN inline Vec3<T> sub(const Vec3<T>& a, const Vec3<T>& b)
{
    return Vec3<T>{a.x - b.x, a.y - b.y, a.z - b.z};
}
template<typename T>
MARS_FS_FN inline T dot(const Vec3<T>& a, const Vec3<T>& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Where a point projected onto a triangle: the barycentric weights of the landing point and how
// far the original point was from it. `inside` is false when the point had to be clamped to an
// edge or vertex -- the weights are still valid and still sum to 1, they just sit on the boundary.
template<typename T>
struct Projection
{
    T    w0, w1, w2;
    T    distance;
    bool inside;
};

// Closest point on triangle (A,B,C) to P, in barycentric form.
//
// The seven-region test from Ericson, Real-Time Collision Detection ch.5: check the three vertex
// Voronoi regions, then the three edge regions, and only then the interior. Written branch-by-
// branch rather than as a projection followed by a clamp, because clamping barycentric
// coordinates independently does NOT give the closest point once the triangle is obtuse.
template<typename T>
MARS_FS_FN Projection<T> projectPointToTriangle(const Vec3<T>& p, const Vec3<T>& a,
                                                const Vec3<T>& b, const Vec3<T>& c)
{
    const Vec3<T> ab = sub(b, a), ac = sub(c, a), ap = sub(p, a);
    const T d1 = dot(ab, ap), d2 = dot(ac, ap);

    Projection<T> r;
    r.inside = false;
    T u, v, w;

    if (d1 <= T(0) && d2 <= T(0)) { u = 1; v = 0; w = 0; }            // vertex A
    else
    {
        const Vec3<T> bp = sub(p, b);
        const T d3 = dot(ab, bp), d4 = dot(ac, bp);
        if (d3 >= T(0) && d4 <= d3) { u = 0; v = 1; w = 0; }          // vertex B
        else
        {
            const Vec3<T> cp = sub(p, c);
            const T d5 = dot(ab, cp), d6 = dot(ac, cp);
            if (d6 >= T(0) && d5 <= d6) { u = 0; v = 0; w = 1; }      // vertex C
            else
            {
                const T vc = d1 * d4 - d3 * d2;
                if (vc <= T(0) && d1 >= T(0) && d3 <= T(0))           // edge AB
                {
                    const T t = d1 / (d1 - d3);
                    u = 1 - t; v = t; w = 0;
                }
                else
                {
                    const T vb = d5 * d2 - d1 * d6;
                    if (vb <= T(0) && d2 >= T(0) && d6 <= T(0))       // edge AC
                    {
                        const T t = d2 / (d2 - d6);
                        u = 1 - t; v = 0; w = t;
                    }
                    else
                    {
                        const T va = d3 * d6 - d5 * d4;
                        if (va <= T(0) && (d4 - d3) >= T(0) && (d5 - d6) >= T(0))   // edge BC
                        {
                            const T t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
                            u = 0; v = 1 - t; w = t;
                        }
                        else                                           // interior
                        {
                            const T den = T(1) / (va + vb + vc);
                            v = vb * den;
                            w = vc * den;
                            u = T(1) - v - w;
                            r.inside = true;
                        }
                    }
                }
            }
        }
    }

    r.w0 = u; r.w1 = v; r.w2 = w;
    const Vec3<T> q{u * a.x + v * b.x + w * c.x, u * a.y + v * b.y + w * c.y,
                    u * a.z + v * b.z + w * c.z};
    const Vec3<T> d = sub(p, q);
    // std::sqrt is constexpr-callable from device code under --expt-relaxed-constexpr, which the
    // build sets; plain sqrt would resolve to the C library on the host and to the device
    // intrinsic in a kernel, which is the same thing here but less obvious.
    r.distance = std::sqrt(dot(d, d));
    return r;
}

// One resolved pairing: which candidate face won, and where on it the point sits.
template<typename T>
struct DonorHit
{
    uint64_t face;        // id of the winning candidate
    T        w0, w1, w2;  // barycentric weights on that face
    T        distance;
    bool     inside;
    bool     found;
};

// Pick the nearest candidate for one point. Candidates come from coarseSearch; `faceIds` and the
// three corner arrays are indexed together. A point is allowed to land outside every candidate --
// the nearest clamped projection still wins, which is what keeps a slightly mismatched interface
// from dropping integration points.
template<typename T>
MARS_FS_FN DonorHit<T> selectDonor(const Vec3<T>& p, const uint64_t* faceIds, const Vec3<T>* cornerA,
                                   const Vec3<T>* cornerB, const Vec3<T>* cornerC, int nCand)
{
    DonorHit<T> best;
    best.found = false;
    best.face  = 0;
    best.w0 = best.w1 = best.w2 = T(0);
    best.distance = T(0);
    best.inside   = false;
    for (int i = 0; i < nCand; ++i)
    {
        const Projection<T> pr = projectPointToTriangle(p, cornerA[i], cornerB[i], cornerC[i]);
        if (!best.found || pr.distance < best.distance)
        {
            best.found    = true;
            best.face     = faceIds[i];
            best.w0       = pr.w0;
            best.w1       = pr.w1;
            best.w2       = pr.w2;
            best.distance = pr.distance;
            best.inside   = pr.inside;
        }
    }
    return best;
}

}  // namespace fem
}  // namespace mars

#endif  // MARS_FINE_SEARCH_HPP
