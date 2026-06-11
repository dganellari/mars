#pragma once

#include "backend/distributed/unstructured/domain.hpp"  // HexTag, TetTag
#include <cuda_runtime.h>
#include <type_traits>

namespace mars {
namespace fem {

// Compile-time element geometry constants, keyed on the element tag. The NS
// solver and assemblers are being templated on ElementTag; these traits
// replace the literal 8 / 12 / 3 that were scattered through the hex-only
// code. Anything that is a pure count lives here; anything that needs a
// per-element formula (volume, area vectors) stays in a tag-dispatched
// function, not a constant.
//
// facesPerCorner is 3 for BOTH hex and tet, but that is a coincidence of the
// node valence inside one element, NOT a shared table -- the (face, sign)
// CONTENTS differ. Do not treat the tables as interchangeable.
template<typename ElementTag>
struct ElemTraits;

template<>
struct ElemTraits<HexTag> {
    static constexpr int NodesPerElem  = 8;
    static constexpr int ScsPerElem    = 12;   // sub-control surfaces (one per edge)
    static constexpr int FacesPerCorner = 3;   // SCS faces incident to each corner
    // Worst-case SCS faces touching one node across all incident elements,
    // for the DDT per-node assembly register buffer. Hex: ~8 elems/corner * 3.
    static constexpr int MaxFacesPerNode = 32;
};

template<>
struct ElemTraits<TetTag> {
    static constexpr int NodesPerElem  = 4;
    static constexpr int ScsPerElem    = 6;
    static constexpr int FacesPerCorner = 3;   // each tet node is an endpoint of 3 of the 6 edges
    // Tet meshes have much higher vertex valence than hex (~20-24 elems/node
    // is common), so the hex bound of 32 silently truncates. Sized generously;
    // the DDT per-node kernel must guard against overflow regardless.
    static constexpr int MaxFacesPerNode = 128;
};

// Device-resident SCS left/right node-pair tables. d_hexLRSCV already lives in
// mars_cvfem_utils.hpp; d_tetLRSCV is added here so the templated scatter
// kernels can pick the right table by tag. Layout matches Tet4CVFEM's host
// lr_pairs: {0,1, 1,2, 0,2, 0,3, 1,3, 2,3}.
__device__ __constant__ int d_tetLRSCV[12] = {
    0, 1, 1, 2, 0, 2, 0, 3, 1, 3, 2, 3
};

} // namespace fem
} // namespace mars
