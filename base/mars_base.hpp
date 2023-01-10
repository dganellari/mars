#ifndef MARS_BASE_HPP
#define MARS_BASE_HPP

#include <cassert>
#include "mars_config.hpp"

/**
 * M.A.R.S Mesh Adaptive Refinement for Simplical meshes
 */
namespace mars {
    using Real = double;
#ifdef _WIN32
    using Integer = int;
#else
    using Integer = long;
#endif
    using Unsigned = unsigned long;
    static constexpr long INVALID_INDEX = -1;

    enum DofOrient : int { xDir = 0, yDir = 1, zDir = 2 };

    // 15 = 1111 binary which is enough if we have only 4 other label representations.
    enum DofLabel : int { lNone = 0, lAll = 15, lVolume = 1, lCorner = 2, lFace = 4, lEdge = 8 };

    enum SLabel : int { Diagonal = 0, Right = 1, Left = 2, Up = 3, Down = 4 };

    enum SSOLabel : int {
        CornerRight = 1,
        CornerLeft = 2,
        VolumeUp = 3,
        VolumeDown = 4,
        FaceRight = 5,
        FaceLeft = 6,
        FaceUp = 7,
        FaceDown = 8,
        FaceUpRight = 12,
        FaceUpLeft = 11,
        FaceDownRight = 10,
        FaceDownLeft = 9
    };

    enum SSXLabel : int {
        VolumeXRight = 1,
        VolumeXLeft = 2,
        CornerXUp = 3,
        CornerXDown = 4,
        FaceXRight = 5,
        FaceXLeft = 6,
        FaceXUp = 7,
        FaceXDown = 8,
        FaceYUpRight = 12,
        FaceYUpLeft = 11,
        FaceYDownRight = 10,
        FaceYDownLeft = 9
    };

    enum SSYLabel : int {
        CornerYRight = 1,
        CornerYLeft = 2,
        VolumeYUp = 3,
        VolumeYDown = 4,
        FaceYRight = 5,
        FaceYLeft = 6,
        FaceYUp = 7,
        FaceYDown = 8,
        FaceXUpRight = 12,
        FaceXUpLeft = 11,
        FaceXDownRight = 10,
        FaceXDownLeft = 9
    };
    /*
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
        }; */

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
