#pragma once

#include "mars_fine_search.hpp"
#include <cstdint>
#include <limits>

namespace mars::fem
{

// IDs identify separate interface sides even when their coordinates coincide.
// Donor faces and their four DOFs must already be available locally or as ghosts.
template<class T>
struct InterfaceTetFace
{
    Vec3<T> x[4];
    int64_t dof[4];
    uint64_t id;
    int side;
    int opposite; // face consists of the other three corners
    T diffusivity;
};

template<class T>
struct InterfacePoint
{
    int face;
    int rowCorner;
    T bary[3]; // ascending local corner order, excluding opposite
    T area;    // physical quadrature area, including the CVFEM test weight
};

enum class InterfaceStatus : int
{
    valid,
    invalid_input,
    invalid_geometry,
    no_donor
};

template<class T>
struct InterfaceStencil
{
    InterfaceStatus status = InterfaceStatus::invalid_input;
    uint64_t donor         = 0;
    int64_t row            = -1;
    int64_t column[8]{};
    T value[8]{};
};

template<class T>
struct InterfacePairingOptions
{
    T distanceTolerance; // physical length; caller must choose from mesh scale
    T normalTolerance = T(64) * std::numeric_limits<T>::epsilon(); // allowed 1 + n_current dot n_donor
};

namespace nonconformal_detail
{
template<class T>
MARS_FS_FN bool finite(T x)
{
    return x == x && x <= std::numeric_limits<T>::max() && x >= -std::numeric_limits<T>::max();
}
template<class T>
MARS_FS_FN Vec3<T> cross(Vec3<T> a, Vec3<T> b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template<class T>
MARS_FS_FN Vec3<T> scale(Vec3<T> a, T b)
{
    return {a.x * b, a.y * b, a.z * b};
}
template<class T>
MARS_FS_FN Vec3<T> add(Vec3<T> a, Vec3<T> b)
{
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

template<class T>
struct Geometry
{
    Vec3<T> gradient[4], normal;
    int corner[3];
    T inverseHeight;
};

template<class T>
MARS_FS_FN InterfaceStatus geometry(const InterfaceTetFace<T>& f, Geometry<T>& g)
{
    if (f.opposite < 0 || f.opposite >= 4 || !finite(f.diffusivity) || f.diffusivity <= T(0))
        return InterfaceStatus::invalid_input;
    T length2 = 0;
    for (int i = 0; i < 4; ++i)
    {
        if (f.dof[i] < 0 || !finite(f.x[i].x) || !finite(f.x[i].y) || !finite(f.x[i].z))
            return InterfaceStatus::invalid_input;
        for (int j = 0; j < i; ++j)
        {
            if (f.dof[i] == f.dof[j]) return InterfaceStatus::invalid_input;
            auto d = sub(f.x[i], f.x[j]);
            T l    = dot(d, d);
            if (l > length2) length2 = l;
        }
    }
    T length = std::sqrt(length2);
    if (!finite(length) || length <= T(0)) return InterfaceStatus::invalid_geometry;
    // Normalize first so the degeneracy criterion is independent of mesh units.
    auto a = scale(sub(f.x[1], f.x[0]), T(1) / length);
    auto b = scale(sub(f.x[2], f.x[0]), T(1) / length);
    auto c = scale(sub(f.x[3], f.x[0]), T(1) / length);
    T det  = dot(a, cross(b, c));
    if (!finite(det) || std::abs(det) <= T(64) * std::numeric_limits<T>::epsilon())
        return InterfaceStatus::invalid_geometry;
    T factor      = (T(1) / det) / length;
    g.gradient[1] = scale(cross(b, c), factor);
    g.gradient[2] = scale(cross(c, a), factor);
    g.gradient[3] = scale(cross(a, b), factor);
    g.gradient[0] = scale(add(add(g.gradient[1], g.gradient[2]), g.gradient[3]), T(-1));
    int n         = 0;
    for (int i = 0; i < 4; ++i)
    {
        if (!finite(g.gradient[i].x) || !finite(g.gradient[i].y) || !finite(g.gradient[i].z))
            return InterfaceStatus::invalid_geometry;
        if (i != f.opposite) g.corner[n++] = i;
    }
    g.inverseHeight = std::sqrt(dot(g.gradient[f.opposite], g.gradient[f.opposite]));
    if (!finite(g.inverseHeight) || g.inverseHeight <= T(0)) return InterfaceStatus::invalid_geometry;
    // N_opposite grows into the element; the outward face normal points the other way.
    g.normal = scale(g.gradient[f.opposite], -T(1) / g.inverseHeight);
    return InterfaceStatus::valid;
}
} // namespace nonconformal_detail

// One current-side CVFEM row, with all four columns from each tetrahedron.
// Independent passes reproduce affine fluxes but need not conserve arbitrary fields
// on nonmatching faces. Use nonconformalConservativePair for shared-overlap assembly.
// Candidates are face indices from coarse search, with remote donors imported beforehand.
template<class T>
MARS_FS_FN InterfaceStencil<T> nonconformalScalarStencil(const InterfacePoint<T>& point,
                                                         const InterfaceTetFace<T>* faces,
                                                         int faceCount,
                                                         const int* candidates,
                                                         int candidateCount,
                                                         InterfacePairingOptions<T> options)
{
    using namespace nonconformal_detail;
    InterfaceStencil<T> out;
    if (!faces || faceCount < 1 || point.face < 0 || point.face >= faceCount || candidateCount < 0 ||
        (candidateCount && !candidates) || !finite(point.area) || point.area <= T(0) ||
        !finite(options.distanceTolerance) || options.distanceTolerance < T(0) || !finite(options.normalTolerance) ||
        options.normalTolerance < T(0) || options.normalTolerance >= T(1))
        return out;
    const auto& a = faces[point.face];
    Geometry<T> ga;
    out.status = geometry(a, ga);
    if (out.status != InterfaceStatus::valid) return out;
    out.status = InterfaceStatus::invalid_input;
    if (point.rowCorner < 0 || point.rowCorner >= 4 || point.rowCorner == a.opposite) return out;
    T sum = 0, na[4]{};
    Vec3<T> x{};
    for (int j = 0; j < 3; ++j)
    {
        if (!finite(point.bary[j]) || point.bary[j] < T(0) || point.bary[j] > T(1)) return out;
        sum += point.bary[j];
        na[ga.corner[j]] = point.bary[j];
        x                = add(x, scale(a.x[ga.corner[j]], point.bary[j]));
    }
    if (std::abs(sum - T(1)) > T(64) * std::numeric_limits<T>::epsilon()) return out;
    int best       = -1;
    T bestDistance = std::numeric_limits<T>::max(), nb[4]{};
    Geometry<T> gb;
    for (int k = 0; k < candidateCount; ++k)
    {
        int index = candidates[k];
        if (index < 0 || index >= faceCount) return out;
        const auto& b = faces[index];
        if (b.side == a.side) continue;
        Geometry<T> g;
        auto status = geometry(b, g);
        if (status != InterfaceStatus::valid)
        {
            out.status = status;
            return out;
        }
        // Separate sides must keep distinct equation IDs, including coincident nodes.
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (a.dof[i] == b.dof[j]) return out;
        if (T(1) + dot(ga.normal, g.normal) > options.normalTolerance) continue;
        auto hit = projectPointToTriangle(x, b.x[g.corner[0]], b.x[g.corner[1]], b.x[g.corner[2]]);
        // Distance to the clamped triangle also rejects points outside its footprint.
        if (!finite(hit.distance) || hit.distance > options.distanceTolerance) continue;
        if (best < 0 || hit.distance < bestDistance || (hit.distance == bestDistance && b.id < faces[best].id))
        {
            best         = index;
            bestDistance = hit.distance;
            gb           = g;
            for (int j = 0; j < 4; ++j)
                nb[j] = 0;
            nb[g.corner[0]] = hit.w0;
            nb[g.corner[1]] = hit.w1;
            nb[g.corner[2]] = hit.w2;
        }
    }
    if (best < 0)
    {
        out.status = InterfaceStatus::no_donor;
        return out;
    }
    const auto& b = faces[best];
    T penalty     = T(0.5) * (a.diffusivity * ga.inverseHeight + b.diffusivity * gb.inverseHeight);
    for (int j = 0; j < 4; ++j)
    {
        out.column[j]     = a.dof[j];
        out.column[4 + j] = b.dof[j];
        out.value[j]      = point.area * (-T(0.5) * a.diffusivity * dot(ga.gradient[j], ga.normal) + penalty * na[j]);
        out.value[4 + j]  = point.area * (T(0.5) * b.diffusivity * dot(gb.gradient[j], gb.normal) - penalty * nb[j]);
        if (!finite(out.value[j]) || !finite(out.value[4 + j]))
        {
            out = InterfaceStencil<T>{};
            return out;
        }
    }
    out.status = InterfaceStatus::valid;
    out.row    = a.dof[point.rowCorner];
    out.donor  = b.id;
    return out;
}

// Host reference for a CSR candidate list. Failures remain per-point and must be
// reduced collectively before any solver consumes the resulting sparse rows.
template<class T>
bool buildNonconformalScalarStencilsHost(int count,
                                         const InterfacePoint<T>* points,
                                         const InterfaceTetFace<T>* faces,
                                         int faceCount,
                                         const int* offsets,
                                         const int* candidates,
                                         int candidateCount,
                                         InterfacePairingOptions<T> options,
                                         InterfaceStencil<T>* output)
{
    if (count < 0 || faceCount < 0 || candidateCount < 0) return false;
    if (count == 0) return true;
    if (!points || !faces || !offsets || !output || (candidateCount && !candidates)) return false;
    for (int i = 0; i < count; ++i)
    {
        output[i] = {};
        if (offsets[i] < 0 || offsets[i + 1] < offsets[i] || offsets[i + 1] > candidateCount) continue;
        output[i] =
            nonconformalScalarStencil(points[i], faces, faceCount, candidates ? candidates + offsets[i] : nullptr,
                                      offsets[i + 1] - offsets[i], options);
    }
    return true;
}

} // namespace mars::fem
