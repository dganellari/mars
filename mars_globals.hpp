#ifndef MARS_GLOBALS_HPP
#define MARS_GLOBALS_HPP

#include <vector>
#include "mars_config.hpp"
#include "mars_base.hpp"

#ifdef WITH_KOKKOS
#include <Kokkos_Core.hpp>
#define MARS_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION 
#else
#define MARS_INLINE_FUNCTION inline
#endif



namespace mars {
    constexpr int hex_n_sides = 6; // 6 faces in total for the hex27.
    constexpr int hex_n_nodes = 27; // 27 nodes for the hex27.
    constexpr int hex_side_n_nodes = 9; // 9 nodes per face for the hex27.


    //FIXME not its place here
    inline void add_side(std::vector<Integer>& side, const Integer a, const Integer b,
            const Integer index) {

        if (a != 1 && b != 1) { //add only nodes which are not mid faces or mid edges
            side.push_back(index);
        } else if (a == 1 && b == 1) { // then add only mid faces
            side.push_back(index);
        }
    }

    MARS_INLINE_FUNCTION
    Integer index(const Integer xDim, const Integer yDim, const Integer i,
            const Integer j, const Integer k) {
        //return k+ (2*zDim +1) * (j + i* (2*yDim + 1));
        return i + (2 * xDim + 1) * (j + k * (2 * yDim + 1));
    }

    //host and device function used for the serial version as well (without kokkos).
    template<typename T>
    MARS_INLINE_FUNCTION
    void build_hex27(T&& nodes, const Integer xDim,
            const Integer yDim, const int i, const int j, const int k) {

        nodes[0] = index(xDim, yDim, i, j, k);
        nodes[1] = index(xDim, yDim, i + 2, j, k);
        nodes[2] = index(xDim, yDim, i + 2, j + 2, k);
        nodes[3] = index(xDim, yDim, i, j + 2, k);
        nodes[4] = index(xDim, yDim, i, j, k + 2);
        nodes[5] = index(xDim, yDim, i + 2, j, k + 2);
        nodes[6] = index(xDim, yDim, i + 2, j + 2, k + 2);
        nodes[7] = index(xDim, yDim, i, j + 2, k + 2);
        nodes[8] = index(xDim, yDim, i + 1, j, k);
        nodes[9] = index(xDim, yDim, i + 2, j + 1, k);
        nodes[10] = index(xDim, yDim, i + 1, j + 2, k);
        nodes[11] = index(xDim, yDim, i, j + 1, k);
        nodes[12] = index(xDim, yDim, i, j, k + 1);
        nodes[13] = index(xDim, yDim, i + 2, j, k + 1);
        nodes[14] = index(xDim, yDim, i + 2, j + 2, k + 1);
        nodes[15] = index(xDim, yDim, i, j + 2, k + 1);
        nodes[16] = index(xDim, yDim, i + 1, j, k + 2);
        nodes[17] = index(xDim, yDim, i + 2, j + 1, k + 2);
        nodes[18] = index(xDim, yDim, i + 1, j + 2, k + 2);
        nodes[19] = index(xDim, yDim, i, j + 1, k + 2);
        nodes[20] = index(xDim, yDim, i + 1, j + 1, k);
        nodes[21] = index(xDim, yDim, i + 1, j, k + 1);
        nodes[22] = index(xDim, yDim, i + 2, j + 1, k + 1);
        nodes[23] = index(xDim, yDim, i + 1, j + 2, k + 1);
        nodes[24] = index(xDim, yDim, i, j + 1, k + 1);
        nodes[25] = index(xDim, yDim, i + 1, j + 1, k + 2);
        nodes[26] = index(xDim, yDim, i + 1, j + 1, k + 1);
    }

    //libmesh method to map the sides to nodes.
    const std::vector<std::vector<unsigned int>> hex_side_nodes{ { 0, 3, 2,
            1, 11, 10, 9, 8, 20 }, // Side 0
            { 0, 1, 5, 4, 8, 13, 16, 12, 21 }, // Side 1
            { 1, 2, 6, 5, 9, 14, 17, 13, 22 }, // Side 2
            { 2, 3, 7, 6, 10, 15, 18, 14, 23 }, // Side 3
            { 3, 0, 4, 7, 11, 12, 19, 15, 24 }, // Side 4
            { 4, 5, 6, 7, 16, 17, 18, 19, 25 }  // Side 5
    };

}



#endif //MARS_GLOBALS_HPP
