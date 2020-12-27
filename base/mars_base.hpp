#ifndef MARS_BASE_HPP
#define MARS_BASE_HPP

#include <cassert>
#include "mars_config.hpp"

/**
 * M.A.R.S Mesh Adaptive Refinement for Simplical meshes
 */
namespace mars {
    using Real = double;
    using Integer = long;
    static const long INVALID_INDEX = -1;

    enum DofOrient : int { xDir = 0, yDir= 1, zDir= 2 };

    enum DofLabel : int { lAll = 0, lVolume = 1, lCorner = 2, lFace = 3 };

    enum ElementType : int {
        // 1D
        Edge2 = 2,
        // 2D
        Quad4 = 4,
        Quad8 = 8,
        Quad9 = 9,
        // 3D
        Hex8 = 8,
        Hex20 = 20,
        Hex27 = 27,
        // possible 4D
        // NonSimplex4D,
        InvalidElem = -1
    };
}  // namespace mars

#endif  // MARS_BASE_HPP
