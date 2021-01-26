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

    enum DofOrient : int { xDir = 0, yDir = 1, zDir = 2 };

    enum DofLabel : int { lNone = -1, lAll = 0, lVolume = 1, lCorner = 2, lFace = 3 };

    enum SLabel : int { Diagonal = 0, Right = 1, Left = 2, Up = 3, Down = 4 };

    enum SSXLabel : int {
        PRight = 1,
        PLeft = 2,
        VXRight = 5,
        VXLeft = 6,
        VXUp = 7,
        VXDown = 8,
        /* VYUpRight = 12,
        VYUpLeft = 11,
        VYDownRight = 10,
        VYDownLeft = 9 */
    };

    enum SSYLabel : int {
        PUp = 3,
        PDown = 4,
        VYRight = 5,
        VYLeft = 6,
        VYUp = 7,
        VYDown = 8,
        /* VXUpRight = 12,
        VXUpLeft = 11,
        VXDownRight = 10,
        VXDownLeft = 9 */
    };

    enum SSXFLabel : int {
        VYUpRight = 3,
        VYUpLeft = 2,
        VYDownRight = 1,
        VYDownLeft = 0
    };

    enum SSYFLabel : int {
        VXUpRight = 3,
        VXUpLeft = 2,
        VXDownRight = 1,
        VXDownLeft = 0
    };

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
