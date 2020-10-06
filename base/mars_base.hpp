#ifndef MARS_BASE_HPP
#define MARS_BASE_HPP

#include "mars_config.hpp"
#include <cassert>

/**
 * M.A.R.S Mesh Adaptive Refinement for Simplical meshes
 */
namespace mars
{
using Real = double;
using Integer = long;
static const long INVALID_INDEX = -1;

enum ElementType : int
{
    //1D
    Edge2 = 2,
    //2D
    Quad4 = 4,
    Quad8 = 8,
    Quad9 = 9,
    //3D
    Hex8 = 8,
    Hex20 = 20,
    Hex27 = 27,
    //possible 4D
    //NonSimplex4D,
    InvalidElem = -1
};
} // namespace mars

#endif //MARS_BASE_HPP
