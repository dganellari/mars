#pragma once
#include "mars_segregated_node.hpp"
#include "mars_segregated_geometry.hpp"
#include "mars_segregated_tet_interior.hpp"
#include "mars_segregated_boundary.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_ASSEMBLY_HD __host__ __device__
#else
#define MARS_ASSEMBLY_HD
#endif

namespace mars::segregated {

struct BoundaryAssemblyInput {
    int element, face, nodes[4];
    BoundaryInput values;
};

MARS_ASSEMBLY_HD inline void native_interior(TetInteriorInput& x, const TetGeometry<double>& g,
    const int* nodes, const double* velocity_gradient, const double* pressure_gradient,
    const double* influence)
{
    for (int s = 0; s < 6; ++s) {
        for (int end = 0; end < 2; ++end) x.edges[2*s+end] = tet_edge_node(s,end);
        tet_sample_shape(s,false,x.velocity_shape+4*s);
        tet_sample_shape(s,false,x.coordinate_shape+4*s);
        for (int j = 0; j < 12; ++j) x.shape_gradient[12*s+j] = g.gradient[j];
        for (int j = 0; j < 3; ++j) x.area[3*s+j] = g.area[3*s+j];
    }
    for (int n = 0; n < 4; ++n) {
        for (int j = 0; j < 9; ++j) x.velocity_gradient[9*n+j] = velocity_gradient[9*nodes[n]+j];
        for (int j = 0; j < 3; ++j) {
            x.pressure_gradient[3*n+j] = pressure_gradient[3*nodes[n]+j];
            if (x.stage == 0) x.influence_lhs[3*n+j] = x.influence_rhs[3*n+j] = influence[3*nodes[n]+j];
        }
    }
}

MARS_ASSEMBLY_HD inline bool native_boundary(BoundaryInput& x, const BoundaryAssemblyInput& input,
    const TetGeometry<double>& g, const int* element_nodes, const double* pressure_gradient,
    const double* influence)
{
    const int count = x.stage == 5 ? 3 : 4;
    for (int s = 0; s < 3; ++s) {
        tet_boundary_area(g,input.face,x.area+3*s);
        int nearest = -1;
        for (int f = 0; f < 3; ++f) if (x.face_nodes[f] == x.nearest[s]) nearest = f;
        if (nearest < 0) return false;
        tri_sample_shape(nearest,false,x.shape+3*s);
    }
    for (int n = 0; n < count; ++n) {
        int local = -1;
        for (int k = 0; k < 4; ++k) if (element_nodes[k] == input.nodes[n]) local = k;
        if (local < 0) return false;
        x.bc_multiplier[n] = 1;
        for (int f = 0; f < 3; ++f) if (x.face_nodes[f] == n) x.bc_multiplier[n] = 0;
        for (int j = 0; j < 3; ++j) {
            x.pressure_gradient[3*n+j] = pressure_gradient[3*input.nodes[n]+j];
            for (int s = 0; s < 3; ++s) x.gradient[12*s+3*n+j] = g.gradient[3*local+j];
        }
    }
    if (x.stage == 1)
        for (int f = 0; f < 3; ++f)
            for (int j = 0; j < 3; ++j)
                x.influence_lhs[3*f+j] = x.influence_rhs[3*f+j] = influence[3*input.nodes[x.face_nodes[f]]+j];
    return true;
}

// Sorted node adjacency; each entry holds a row-major Components x Components block.
template<int Components> struct BlockCsrView {
    int nodes;
    const int *offsets, *columns;
    double *values, *rhs;
};

template<int Components>
MARS_ASSEMBLY_HD inline int block_position(BlockCsrView<Components> a, int row, int column)
{
    if (row < 0 || column < 0 || row >= a.nodes || column >= a.nodes) return -1;
    int lo = a.offsets[row], hi = a.offsets[row+1];
    const int end = hi;
    while (lo < hi) {
        const int mid = lo+(hi-lo)/2;
        if (a.columns[mid] < column) lo = mid+1; else hi = mid;
    }
    return lo < end && a.columns[lo] == column ? lo : -1;
}

MARS_ASSEMBLY_HD inline void assembly_add(double* address, double value)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    atomicAdd(address,value);
#else
    *address += value;
#endif
}

template<int Components>
MARS_ASSEMBLY_HD inline bool scatter_block(BlockCsrView<Components> a, const int* nodes,
    int count, const double* lhs, const double* rhs)
{
    if (count < 1 || count > 4) return false;
    int positions[16];
    // Check every entry, including zeros: a missing pressure anchor must never disappear silently.
    for (int r = 0; r < count; ++r)
        for (int c = 0; c < count; ++c) {
            positions[count*r+c] = block_position(a,nodes[r],nodes[c]);
            if (positions[count*r+c] < 0) return false;
        }
    const int width = count*Components;
    for (int r = 0; r < count; ++r)
        for (int i = 0; i < Components; ++i) {
            assembly_add(a.rhs+Components*nodes[r]+i,rhs[Components*r+i]);
            for (int c = 0; c < count; ++c)
                for (int j = 0; j < Components; ++j)
                    assembly_add(a.values+Components*Components*positions[count*r+c]+Components*i+j,
                                 lhs[width*(Components*r+i)+Components*c+j]);
        }
    return true;
}

MARS_ASSEMBLY_HD inline bool finish_momentum_row(BlockCsrView<3> a, int row,
    double volume, double alpha, double boundary_factor, bool consistent,
    double* influence, double* influence_tilde)
{
    if (!(alpha > 0 && alpha <= 1) || !(boundary_factor > 0 && boundary_factor <= 1)) return false;
    const int diagonal = block_position(a,row,row);
    if (diagonal < 0) return false;
    relax_momentum_diagonal(a.values+9*diagonal,alpha);
    if (!momentum_influence(volume,a.values+9*a.offsets[row],a.offsets[row+1]-a.offsets[row],
                            diagonal-a.offsets[row],consistent,influence,influence_tilde)) return false;
    relax_boundary_rhs(a.rhs+3*row,boundary_factor);
    return true;
}
} // namespace mars::segregated
#undef MARS_ASSEMBLY_HD
