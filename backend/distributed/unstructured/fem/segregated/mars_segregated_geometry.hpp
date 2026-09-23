#pragma once
#include <cmath>
#include <limits>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_GEOMETRY_HD __host__ __device__
#else
#define MARS_GEOMETRY_HD
#endif

namespace mars::segregated {

struct TetBoundaryFace { int element, ordinal; };

template<class RealType> struct TetGeometry {
    RealType volume{}, gradient[12]{}, area[18]{};
};

MARS_GEOMETRY_HD inline int tet_edge_node(int sample, int end)
{
    constexpr int nodes[] = {0,1, 1,2, 0,2, 0,3, 1,3, 2,3};
    return nodes[2*sample+end];
}

MARS_GEOMETRY_HD inline int tet_face_node(int face, int node)
{
    constexpr int nodes[] = {0,1,3, 1,2,3, 0,3,2, 0,2,1};
    return nodes[3*face+node];
}

MARS_GEOMETRY_HD inline int tet_opposite_node(int face)
{
    constexpr int nodes[] = {2,0,1,3};
    return nodes[face];
}

template<class RealType>
MARS_GEOMETRY_HD inline bool geometry_finite(RealType value)
{
    constexpr RealType bound = std::numeric_limits<RealType>::max();
    return value >= -bound && value <= bound;
}

// Positive Tet4 orientation is required; silently flipping nodes would change sample identities.
template<class RealType>
MARS_GEOMETRY_HD inline bool tet_geometry(const RealType* coordinates, TetGeometry<RealType>& out)
{
    out = {};
    RealType a[3], b[3], c[3], cross[3][3];
    for (int i = 0; i < 12; ++i) if (!geometry_finite(coordinates[i])) return false;
    for (int j = 0; j < 3; ++j) {
        a[j] = coordinates[3+j]-coordinates[j];
        b[j] = coordinates[6+j]-coordinates[j];
        c[j] = coordinates[9+j]-coordinates[j];
    }
    for (int j = 0; j < 3; ++j) {
        const int k = (j+1)%3, l = (j+2)%3;
        cross[0][j] = b[k]*c[l]-b[l]*c[k];
        cross[1][j] = c[k]*a[l]-c[l]*a[k];
        cross[2][j] = a[k]*b[l]-a[l]*b[k];
    }
    const RealType determinant = a[0]*cross[0][0]+a[1]*cross[0][1]+a[2]*cross[0][2];
    if (!(determinant > 0) || !geometry_finite(determinant)) return false;
    out.volume = determinant/RealType(6);
    if (!(out.volume > 0)) return false;
    for (int j = 0; j < 3; ++j) {
        for (int n = 1; n < 4; ++n) out.gradient[3*n+j] = cross[n-1][j]/determinant;
        out.gradient[j] = -out.gradient[3+j]-out.gradient[6+j]-out.gradient[9+j];
    }
    // Affine image of the median-dual quadrilateral: A_LR = V/4 (grad N_R - grad N_L).
    for (int s = 0; s < 6; ++s)
        for (int j = 0; j < 3; ++j)
            out.area[3*s+j] = (out.volume/RealType(4))
                *(out.gradient[3*tet_edge_node(s,1)+j]-out.gradient[3*tet_edge_node(s,0)+j]);
    for (RealType value : out.gradient) if (!geometry_finite(value)) return false;
    for (RealType value : out.area) if (!geometry_finite(value)) return false;
    return true;
}

template<class RealType>
MARS_GEOMETRY_HD inline void tet_sample_shape(int sample, bool shifted, RealType* shape)
{
    const int left = tet_edge_node(sample,0), right = tet_edge_node(sample,1);
    for (int n = 1; n < 4; ++n)
        shape[n] = shifted ? ((n == left || n == right) ? RealType(0.5) : RealType(0))
                           : ((n == left || n == right) ? RealType(13)/36 : RealType(5)/36);
    shape[0] = RealType(1)-shape[1]-shape[2]-shape[3];
}

template<class RealType>
MARS_GEOMETRY_HD inline void tri_sample_shape(int sample, bool shifted, RealType* shape)
{
    for (int n = 1; n < 3; ++n)
        shape[n] = shifted ? (n == sample ? RealType(1) : RealType(0))
                           : (n == sample ? RealType(11)/18 : RealType(7)/36);
    shape[0] = RealType(1)-shape[1]-shape[2];
}

template<class RealType>
MARS_GEOMETRY_HD inline void tet_boundary_area(const TetGeometry<RealType>& geometry, int face, RealType* area)
{
    // Each of the three Tri3 samples receives one third of the outward face vector.
    for (int j = 0; j < 3; ++j)
        area[j] = -geometry.volume*geometry.gradient[3*tet_opposite_node(face)+j];
}

// Component-major gradients per node; numerator has units field * area.
template<int Components, class RealType>
MARS_GEOMETRY_HD inline void tet_gradient_numerator(const TetGeometry<RealType>& geometry,
    const RealType* field, bool shifted, bool incremental, RealType* numerator)
{
    for (int i = 0; i < 4*Components*3; ++i) numerator[i] = 0;
    for (int s = 0; s < 6; ++s) {
        RealType shape[4]; tet_sample_shape(s, shifted, shape);
        const int left = tet_edge_node(s,0), right = tet_edge_node(s,1);
        for (int c = 0; c < Components; ++c) {
            RealType value = 0;
            for (int n = 0; n < 4; ++n) value += shape[n]*field[n*Components+c];
            const RealType l = value-(incremental ? field[left*Components+c] : RealType(0));
            const RealType r = value-(incremental ? field[right*Components+c] : RealType(0));
            for (int j = 0; j < 3; ++j) {
                numerator[(left*Components+c)*3+j] += l*geometry.area[3*s+j];
                numerator[(right*Components+c)*3+j] -= r*geometry.area[3*s+j];
            }
        }
    }
}

template<int Components, class RealType>
MARS_GEOMETRY_HD inline void tri_gradient_numerator(const RealType* area,
    const RealType* field, bool shifted, bool incremental, RealType* numerator)
{
    for (int s = 0; s < 3; ++s) {
        RealType shape[3]; tri_sample_shape(s, shifted, shape);
        for (int c = 0; c < Components; ++c) {
            RealType value = 0;
            for (int n = 0; n < 3; ++n) value += shape[n]*field[n*Components+c];
            value -= incremental ? field[s*Components+c] : RealType(0);
            for (int j = 0; j < 3; ++j) numerator[(s*Components+c)*3+j] = value*area[j];
        }
    }
}

template<class RealType>
MARS_GEOMETRY_HD inline bool finish_gradient(RealType numerator, RealType volume, RealType previous,
    RealType relaxation, bool initialized, RealType& result)
{
    if (!(volume > 0) || !geometry_finite(volume) || !geometry_finite(numerator)
        || !(relaxation >= 0 && relaxation <= 1)) return false;
    const RealType fresh = numerator/volume;
    // First reconstruction and the unrelaxed pressure increment must not read stale history.
    result = (!initialized || relaxation == RealType(1)) ? fresh
        : (RealType(1)-relaxation)*previous+relaxation*fresh;
    return geometry_finite(result);
}
} // namespace mars::segregated
#undef MARS_GEOMETRY_HD
