#ifndef MARS_MORTON_CODE_HPP
#define MARS_MORTON_CODE_HPP

#include "mars_base.hpp"
#include "mars_config.hpp"
#include "mars_globals.hpp"

namespace mars {

    MARS_INLINE_FUNCTION
    Unsigned interleave_2D(Unsigned x) {
        x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
        x = (x | (x << 8)) & 0x00FF00FF00FF00FF;
        x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0F;
        x = (x | (x << 2)) & 0x3333333333333333;
        x = (x | (x << 1)) & 0x5555555555555555;
        return x;
    }

    MARS_INLINE_FUNCTION
    Unsigned encode_morton_2D(Unsigned x, Unsigned y) { return (interleave_2D(y) << 1) | interleave_2D(x); }

    MARS_INLINE_FUNCTION
    Unsigned compact_2D(Unsigned x) {
        x &= 0x5555555555555555;
        x = (x ^ (x >> 1)) & 0x3333333333333333;  // x = --fe --dc --ba --98 --76 --54 --32 --10
        x = (x ^ (x >> 2)) & 0x0F0F0F0F0F0F0F0F;
        x = (x ^ (x >> 4)) & 0x00FF00FF00FF00FF;
        x = (x ^ (x >> 8)) & 0x0000FFFF0000FFFF;
        x = (x ^ (x >> 16)) & 0x00000000FFFFFFFF;
        return x;
    }

    MARS_INLINE_FUNCTION
    Unsigned decode_morton_2DX(Unsigned code) { return compact_2D(code >> 0); }

    MARS_INLINE_FUNCTION
    Unsigned decode_morton_2DY(Unsigned code) { return compact_2D(code >> 1); }

    // method to seperate bits from a given integer 3 positions apart
    MARS_INLINE_FUNCTION Unsigned interleave_3D(Unsigned a) {
        Unsigned x = a & 0x1fffff;
        x = (x | x << 32) & 0x1f00000000ffff;
        x = (x | x << 16) & 0x1f0000ff0000ff;
        x = (x | x << 8) & 0x100f00f00f00f00f;
        x = (x | x << 4) & 0x10c30c30c30c30c3;
        x = (x | x << 2) & 0x1249249249249249;
        return x;
    }

    MARS_INLINE_FUNCTION
    Unsigned encode_morton_3D(Unsigned x, Unsigned y, Unsigned z) {
        Unsigned sfc = 0;
        sfc |= interleave_3D(x) | interleave_3D(y) << 1 | interleave_3D(z) << 2;
        return sfc;
    }

    MARS_INLINE_FUNCTION Unsigned compact_3D(Unsigned n) {
        n &= 0x1249249249249249;
        n = (n ^ (n >> 2)) & 0x30c30c30c30c30c3;
        n = (n ^ (n >> 4)) & 0xf00f00f00f00f00f;
        n = (n ^ (n >> 8)) & 0x00ff0000ff0000ff;
        n = (n ^ (n >> 16)) & 0x00ff00000000ffff;
        n = (n ^ (n >> 32)) & 0x1fffff;
        return n;
    }

    MARS_INLINE_FUNCTION
    Unsigned decode_morton_3DX(Unsigned code) { return compact_3D(code >> 0); }

    MARS_INLINE_FUNCTION
    Unsigned decode_morton_3DY(Unsigned code) { return compact_3D(code >> 1); }

    MARS_INLINE_FUNCTION
    Unsigned decode_morton_3DZ(Unsigned code) { return compact_3D(code >> 2); }

    /******************** morton z curve from fgiesen.wordpress.com with slight
     * modifications*****************************/

    // "Insert" a 0 bit after each of the 16 low bits of x
    MARS_INLINE_FUNCTION
    unsigned int interleave_by1(unsigned int x) {
        x &= 0x0000ffff;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
        x = (x ^ (x << 8)) & 0x00ff00ff;  // x = ---- ---- fedc ba98 ---- ---- 7654 3210
        x = (x ^ (x << 4)) & 0x0f0f0f0f;  // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
        x = (x ^ (x << 2)) & 0x33333333;  // x = --fe --dc --ba --98 --76 --54 --32 --10
        x = (x ^ (x << 1)) & 0x55555555;  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
        return x;
    }

    // "Insert" two 0 bits after each of the 10 low bits of x
    // 10bit number and 2 interleaving bits (for 32 bit) x<=1024.
    MARS_INLINE_FUNCTION
    unsigned int interleave_by2(unsigned int x) {
        x &= 0x000003ff;                   // x = ---- ---- ---- ---- ---- --98 7654 3210
        x = (x ^ (x << 16)) & 0xff0000ff;  // x = ---- --98 ---- ---- ---- ---- 7654 3210
        x = (x ^ (x << 8)) & 0x0300f00f;   // x = ---- --98 ---- ---- 7654 ---- ---- 3210
        x = (x ^ (x << 4)) & 0x030c30c3;   // x = ---- --98 ---- 76-- --54 ---- 32-- --10
        x = (x ^ (x << 2)) & 0x09249249;   // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
        return x;
    }

    /* MARS_INLINE_FUNCTION
    unsigned int encode_morton_2D(unsigned int x, unsigned int y)
    {
        return (interleave_by1(y) << 1) | interleave_by1(x);
    } */

    /* MARS_INLINE_FUNCTION
    unsigned int encode_morton_3D(unsigned int x, unsigned int y, unsigned int z)
    {
        return (interleave_by2(z) << 2) | (interleave_by2(y) << 1) | interleave_by2(x);
    } */

    // Inverse of Part1By1 - "delete" all odd-indexed bits
    MARS_INLINE_FUNCTION
    unsigned int recall_by1(unsigned int x) {
        x &= 0x55555555;                  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
        x = (x ^ (x >> 1)) & 0x33333333;  // x = --fe --dc --ba --98 --76 --54 --32 --10
        x = (x ^ (x >> 2)) & 0x0f0f0f0f;  // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
        x = (x ^ (x >> 4)) & 0x00ff00ff;  // x = ---- ---- fedc ba98 ---- ---- 7654 3210
        x = (x ^ (x >> 8)) & 0x0000ffff;  // x = ---- ---- ---- ---- fedc ba98 7654 3210
        return x;
    }

    // Inverse of Part1By2 - "delete" all bits not at positions divisible by 3
    MARS_INLINE_FUNCTION
    unsigned int recall_by2(unsigned int x) {
        x &= 0x09249249;                   // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
        x = (x ^ (x >> 2)) & 0x030c30c3;   // x = ---- --98 ---- 76-- --54 ---- 32-- --10
        x = (x ^ (x >> 4)) & 0x0300f00f;   // x = ---- --98 ---- ---- 7654 ---- ---- 3210
        x = (x ^ (x >> 8)) & 0xff0000ff;   // x = ---- --98 ---- ---- ---- ---- 7654 3210
        x = (x ^ (x >> 16)) & 0x000003ff;  // x = ---- ---- ---- ---- ---- --98 7654 3210
        return x;
    }

    /* MARS_INLINE_FUNCTION
    unsigned int decode_morton_2DX(unsigned int code)
    {
        return recall_by1(code >> 0);
    }

    MARS_INLINE_FUNCTION
    unsigned int decode_morton_2DY(unsigned int code)
    {
        return recall_by1(code >> 1);
    }

    MARS_INLINE_FUNCTION
    unsigned int decode_morton_3DX(unsigned int code)
    {
        return recall_by2(code >> 0);
    }

    MARS_INLINE_FUNCTION
    unsigned int decode_morton_3DY(unsigned int code)
    {
        return recall_by2(code >> 1);
    }

    MARS_INLINE_FUNCTION
    unsigned int decode_morton_3DZ(unsigned int code)
    {
        return recall_by2(code >> 2);
    } */

}  // namespace mars
#endif
