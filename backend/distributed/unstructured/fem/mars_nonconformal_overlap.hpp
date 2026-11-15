#pragma once

#include "mars_nonconformal_scalar.hpp"

namespace mars::fem
{
// Each triangular face has six median-dual triangles. A pair of such triangles
// has one convex overlap, contributing one row on each side.
inline constexpr int nonconformal_pair_capacity = 72;
struct InterfaceFacePair
{
    int first, second;
};

// One physical overlap, oriented out of side zero. The two rows use the same
// quadrature point and opposite signs; do not recompute a reverse-side flux.
template<class T>
struct InterfacePatch
{
    int64_t row[2]{-1, -1}, column[8]{};
    T shape[8]{}, area{}, inverseHeight[2]{};
    Vec3<T> gradient[8]{}, normal{};
};

namespace nonconformal_overlap_detail
{
template<class T>
struct Point2
{
    T x, y;
};
template<class T>
MARS_FS_FN T cross2(Point2<T> a, Point2<T> b, Point2<T> c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}
template<class T>
MARS_FS_FN Point2<T> mean(Point2<T> a, Point2<T> b)
{
    return {(a.x + b.x) / T(2), (a.y + b.y) / T(2)};
}
template<class T>
MARS_FS_FN void patch(const Point2<T>* face, int row, int edge, Point2<T>* out)
{
    out[0] = face[row];
    out[1] = mean(face[row], face[(row + edge) % 3]);
    out[2] = {(face[0].x + face[1].x + face[2].x) / T(3), (face[0].y + face[1].y + face[2].y) / T(3)};
    if (cross2(out[0], out[1], out[2]) < T(0))
    {
        auto tmp = out[1];
        out[1]   = out[2];
        out[2]   = tmp;
    }
}
// Work in a face-local orthonormal plane, so areas are physical and large
// absolute mesh coordinates do not enter polygon moments.
template<class T>
MARS_FS_FN bool overlap(const Point2<T>* a, const Point2<T>* b, T& area, Point2<T>& center)
{
    Point2<T> polygon[8], next[8];
    int size = 3;
    for (int i = 0; i < 3; ++i)
        polygon[i] = a[i];
    for (int edge = 0; edge < 3 && size; ++edge)
    {
        int count = 0;
        auto prev = polygon[size - 1];
        T dp      = cross2(b[edge], b[(edge + 1) % 3], prev);
        for (int i = 0; i < size; ++i)
        {
            auto cur = polygon[i];
            T dc     = cross2(b[edge], b[(edge + 1) % 3], cur);
            if ((dp >= T(0)) != (dc >= T(0)))
            {
                T t = dp / (dp - dc);
                if (count == 8 || !nonconformal_detail::finite(t)) return false;
                next[count++] = {prev.x + t * (cur.x - prev.x), prev.y + t * (cur.y - prev.y)};
            }
            if (dc >= T(0))
            {
                if (count == 8) return false;
                next[count++] = cur;
            }
            prev = cur;
            dp   = dc;
        }
        size = count;
        for (int i = 0; i < size; ++i)
            polygon[i] = next[i];
    }
    area   = T(0);
    center = {T(0), T(0)};
    for (int i = 1; i + 1 < size; ++i)
    {
        T left   = (polygon[i].x - polygon[0].x) * (polygon[i + 1].y - polygon[0].y);
        T right  = (polygon[i].y - polygon[0].y) * (polygon[i + 1].x - polygon[0].x);
        T w      = (left - right) / T(2);
        T extent = T(0);
        for (int k = 0; k < 3; ++k)
            extent += std::abs(a[k].x - a[0].x) + std::abs(a[k].y - a[0].y) + std::abs(b[k].x - a[0].x) +
                      std::abs(b[k].y - a[0].y);
        T bound = T(32) * std::numeric_limits<T>::epsilon() * extent * extent;
        if (!nonconformal_detail::finite(w) || w < -bound) return false;
        // Clipping can leave collinear vertices a few ulps out of order.
        if (w < T(0)) w = T(0);
        T total = area + w;
        if (w > T(0))
        {
            Point2<T> c{(polygon[0].x + polygon[i].x + polygon[i + 1].x) / T(3),
                        (polygon[0].y + polygon[i].y + polygon[i + 1].y) / T(3)};
            center.x += (w / total) * (c.x - center.x);
            center.y += (w / total) * (c.y - center.y);
        }
        area = total;
    }
    return nonconformal_detail::finite(area) && nonconformal_detail::finite(center.x) &&
           nonconformal_detail::finite(center.y);
}
} // namespace nonconformal_overlap_detail

// Steady scalar diffusion: q.n = -k grad(u).n. Use the existing average flux
// and jump penalty, but integrate once on the intersection of the two median
// dual patches. Every emitted row has an exact opposite on the other side.
// Caller supplies a unique face pair in a consistent order on every rank and
// discards all entries unless status is valid. Nonoverlapping faces emit none.
template<class T>
MARS_FS_FN InterfaceStatus nonconformalConservativePair(const InterfaceTetFace<T>& a,
                                                        const InterfaceTetFace<T>& b,
                                                        InterfacePairingOptions<T> options,
                                                        InterfaceStencil<T>* out,
                                                        int capacity,
                                                        int& count,
                                                        InterfacePatch<T>* patches = nullptr)
{
    using namespace nonconformal_detail;
    using namespace nonconformal_overlap_detail;
    count = 0;
    if (!out || capacity < nonconformal_pair_capacity || a.id == b.id || a.side == b.side ||
        !finite(options.distanceTolerance) || options.distanceTolerance < T(0) || !finite(options.normalTolerance) ||
        options.normalTolerance < T(0) || options.normalTolerance >= T(1))
        return InterfaceStatus::invalid_input;
    Geometry<T> ga, gb;
    if (geometry(a, ga) != InterfaceStatus::valid || geometry(b, gb) != InterfaceStatus::valid)
        return InterfaceStatus::invalid_geometry;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (a.dof[i] == b.dof[j]) return InterfaceStatus::invalid_input;
    if (T(1) + dot(ga.normal, gb.normal) > options.normalTolerance) return InterfaceStatus::invalid_geometry;
    auto origin = a.x[ga.corner[0]];
    auto u      = sub(a.x[ga.corner[1]], origin);
    u           = scale(u, T(1) / std::sqrt(dot(u, u)));
    auto v      = cross(ga.normal, u);
    Point2<T> fa[3], fb[3];
    for (int i = 0; i < 3; ++i)
    {
        auto da = sub(a.x[ga.corner[i]], origin), db = sub(b.x[gb.corner[i]], origin);
        if (std::abs(dot(db, ga.normal)) > options.distanceTolerance) return InterfaceStatus::invalid_geometry;
        fa[i] = {dot(da, u), dot(da, v)};
        fb[i] = {dot(db, u), dot(db, v)};
    }
    InterfaceTetFace<T> faces[]{a, b};
    int candidate = 1, written = 0;
    for (int ra = 0; ra < 3; ++ra)
        for (int ea = 1; ea < 3; ++ea)
            for (int rb = 0; rb < 3; ++rb)
                for (int eb = 1; eb < 3; ++eb)
                {
                    Point2<T> pa[3], pb[3], center;
                    patch(fa, ra, ea, pa);
                    patch(fb, rb, eb, pb);
                    T area;
                    if (!overlap(pa, pb, area, center)) return InterfaceStatus::invalid_geometry;
                    if (area == T(0)) continue;
                    auto x   = add(origin, add(scale(u, center.x), scale(v, center.y)));
                    auto hit = projectPointToTriangle(x, a.x[ga.corner[0]], a.x[ga.corner[1]], a.x[ga.corner[2]]);
                    InterfacePoint<T> point{0, ga.corner[ra], {hit.w0, hit.w1, hit.w2}, area};
                    auto s = nonconformalScalarStencil(point, faces, 2, &candidate, 1, options);
                    if (s.status != InterfaceStatus::valid) return s.status;
                    if (patches)
                    {
                        auto& p            = patches[written / 2];
                        p                  = {};
                        p.row[0]           = s.row;
                        p.row[1]           = b.dof[gb.corner[rb]];
                        p.area             = area;
                        p.normal           = ga.normal;
                        p.inverseHeight[0] = ga.inverseHeight;
                        p.inverseHeight[1] = gb.inverseHeight;
                        for (int j = 0; j < 4; ++j)
                        {
                            p.column[j]       = a.dof[j];
                            p.column[4 + j]   = b.dof[j];
                            p.gradient[j]     = ga.gradient[j];
                            p.gradient[4 + j] = gb.gradient[j];
                        }
                        auto bhit = projectPointToTriangle(x, b.x[gb.corner[0]], b.x[gb.corner[1]], b.x[gb.corner[2]]);
                        T wa[]{hit.w0, hit.w1, hit.w2}, wb[]{bhit.w0, bhit.w1, bhit.w2};
                        T sa = 0, sb = 0;
                        for (int j = 0; j < 3; ++j)
                        {
                            if (wa[j] < -T(128) * std::numeric_limits<T>::epsilon() ||
                                wb[j] < -T(128) * std::numeric_limits<T>::epsilon())
                                return InterfaceStatus::invalid_geometry;
                            wa[j] = wa[j] < 0 ? T(0) : wa[j];
                            wb[j] = wb[j] < 0 ? T(0) : wb[j];
                            sa += wa[j];
                            sb += wb[j];
                        }
                        for (int j = 0; j < 3; ++j)
                        {
                            p.shape[ga.corner[j]]     = wa[j] / sa;
                            p.shape[4 + gb.corner[j]] = wb[j] / sb;
                        }
                    }
                    out[written++] = s;
                    s.row          = b.dof[gb.corner[rb]];
                    s.donor        = a.id;
                    for (int j = 0; j < 8; ++j)
                        s.value[j] = -s.value[j];
                    out[written++] = s;
                }
    count = written;
    return InterfaceStatus::valid;
}
} // namespace mars::fem
