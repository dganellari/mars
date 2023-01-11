#ifndef MARS_MORTON_CODE_HPP
#define MARS_MORTON_CODE_HPP

#include "mars_distributed_octant.hpp"
#include "mars_sfc_utils.hpp"

namespace mars {

    /* -------------------------------------64 bit ------------------------------------- */
    MARS_INLINE_FUNCTION
    Unsigned interleave_2D(Unsigned x) {
        x = (x | (x << 16u)) & 0x0000FFFF0000FFFFlu;
        x = (x | (x << 8u)) & 0x00FF00FF00FF00FFlu;
        x = (x | (x << 4u)) & 0x0F0F0F0F0F0F0F0Flu;
        x = (x | (x << 2u)) & 0x3333333333333333lu;
        x = (x | (x << 1u)) & 0x5555555555555555lu;
        return x;
    }

    MARS_INLINE_FUNCTION
    Unsigned compact_2D(Unsigned x) {
        x &= 0x5555555555555555lu;
        x = (x ^ (x >> 1u)) & 0x3333333333333333lu;  // x = --fe --dc --ba --98 --76 --54 --32 --10
        x = (x ^ (x >> 2u)) & 0x0F0F0F0F0F0F0F0Flu;
        x = (x ^ (x >> 4u)) & 0x00FF00FF00FF00FFlu;
        x = (x ^ (x >> 8u)) & 0x0000FFFF0000FFFFlu;
        x = (x ^ (x >> 16u)) & 0x00000000FFFFFFFFlu;
        return x;
    }

    MARS_INLINE_FUNCTION
    Unsigned decode_morton_2DX(Unsigned code) { return compact_2D(code >> 0); }

    MARS_INLINE_FUNCTION
    Unsigned decode_morton_2DY(Unsigned code) { return compact_2D(code >> 1); }

    // method to seperate bits from a given integer 3 positions apart
    MARS_INLINE_FUNCTION Unsigned interleave_3D(Unsigned a) {
        Unsigned x = a & 0x1fffffu;
        x = (x | x << 32u) & 0x1f00000000fffflu;
        x = (x | x << 16u) & 0x1f0000ff0000fflu;
        x = (x | x << 8u) & 0x100f00f00f00f00flu;
        x = (x | x << 4u) & 0x10c30c30c30c30c3lu;
        x = (x | x << 2u) & 0x1249249249249249lu;
        return x;
    }

    MARS_INLINE_FUNCTION Unsigned compact_3D(Unsigned n) {
        n &= 0x1249249249249249lu;
        n = (n ^ (n >> 2u)) & 0x30c30c30c30c30c3lu;
        n = (n ^ (n >> 4u)) & 0xf00f00f00f00f00flu;
        n = (n ^ (n >> 8u)) & 0x00ff0000ff0000fflu;
        n = (n ^ (n >> 16u)) & 0x00ff00000000fffflu;
        n = (n ^ (n >> 32u)) & 0x1ffffflu;
        return n;
    }

    MARS_INLINE_FUNCTION
    Unsigned decode_morton_3DX(Unsigned code) { return compact_3D(code >> 0); }

    MARS_INLINE_FUNCTION
    Unsigned decode_morton_3DY(Unsigned code) { return compact_3D(code >> 1); }

    MARS_INLINE_FUNCTION
    Unsigned decode_morton_3DZ(Unsigned code) { return compact_3D(code >> 2); }

    /* -------------------------------------32 bit ------------------------------------- */

    // "Insert" a 0 bit after each of the 16 low bits of x
    MARS_INLINE_FUNCTION
    unsigned interleave_2D(unsigned x) {
        x &= 0x0000ffffu;                   // x = ---- ---- ---- ---- fedc ba98 7654 3210
        x = (x ^ (x << 8u)) & 0x00ff00ffu;  // x = ---- ---- fedc ba98 ---- ---- 7654 3210
        x = (x ^ (x << 4u)) & 0x0f0f0f0fu;  // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
        x = (x ^ (x << 2u)) & 0x33333333u;  // x = --fe --dc --ba --98 --76 --54 --32 --10
        x = (x ^ (x << 1u)) & 0x55555555u;  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
        return x;
    }

    // "Insert" two 0 bits after each of the 10 low bits of x
    // 10bit number and 2 interleaving bits (for 32 bit) x<=1024.
    MARS_INLINE_FUNCTION
    unsigned interleave_3D(unsigned x) {
        x &= 0x000003ffu;                    // x = ---- ---- ---- ---- ---- --98 7654 3210
        x = (x ^ (x << 16u)) & 0xff0000ffu;  // x = ---- --98 ---- ---- ---- ---- 7654 3210
        x = (x ^ (x << 8u)) & 0x0300f00fu;   // x = ---- --98 ---- ---- 7654 ---- ---- 3210
        x = (x ^ (x << 4u)) & 0x030c30c3u;   // x = ---- --98 ---- 76-- --54 ---- 32-- --10
        x = (x ^ (x << 2u)) & 0x09249249u;   // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
        return x;
    }

    // Inverse of Part1By1 - "delete" all odd-indexed bits
    MARS_INLINE_FUNCTION
    unsigned compact_2D(unsigned x) {
        x &= 0x55555555u;                   // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
        x = (x ^ (x >> 1u)) & 0x33333333u;  // x = --fe --dc --ba --98 --76 --54 --32 --10
        x = (x ^ (x >> 2u)) & 0x0f0f0f0fu;  // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
        x = (x ^ (x >> 4u)) & 0x00ff00ffu;  // x = ---- ---- fedc ba98 ---- ---- 7654 3210
        x = (x ^ (x >> 8u)) & 0x0000ffffu;  // x = ---- ---- ---- ---- fedc ba98 7654 3210
        return x;
    }

    // Inverse of Part1By2 - "delete" all bits not at positions divisible by 3
    MARS_INLINE_FUNCTION
    unsigned compact_3D(unsigned x) {
        x &= 0x09249249u;                    // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
        x = (x ^ (x >> 2u)) & 0x030c30c3u;   // x = ---- --98 ---- 76-- --54 ---- 32-- --10
        x = (x ^ (x >> 4u)) & 0x0300f00fu;   // x = ---- --98 ---- ---- 7654 ---- ---- 3210
        x = (x ^ (x >> 8u)) & 0xff0000ffu;   // x = ---- --98 ---- ---- ---- ---- 7654 3210
        x = (x ^ (x >> 16u)) & 0x000003ffu;  // x = ---- ---- ---- ---- ---- --98 7654 3210
        return x;
    }

    MARS_INLINE_FUNCTION
    unsigned decode_morton_2DX(unsigned code) { return compact_2D(code >> 0); }

    MARS_INLINE_FUNCTION
    unsigned decode_morton_2DY(unsigned code) { return compact_2D(code >> 1); }

    MARS_INLINE_FUNCTION
    unsigned decode_morton_3DX(unsigned code) { return compact_3D(code >> 0); }

    MARS_INLINE_FUNCTION
    unsigned decode_morton_3DY(unsigned code) { return compact_3D(code >> 1); }

    MARS_INLINE_FUNCTION
    unsigned decode_morton_3DZ(unsigned code) { return compact_3D(code >> 2); }

    /* -------------------------3D Morton encoding/decoding in 32- and 64-bit using the magic number method---- */

    template <typename KeyType = Unsigned>
    MARS_INLINE_FUNCTION std::enable_if_t<std::is_unsigned_v<KeyType>, KeyType> encode_morton_2D(unsigned x,
                                                                                                 unsigned y) {
        return (interleave_2D(KeyType(y)) << 1) | interleave_2D(KeyType(x));
    }

    template <typename KeyType = Unsigned>
    MARS_INLINE_FUNCTION std::enable_if_t<std::is_unsigned_v<KeyType>, KeyType> encode_morton_3D(unsigned x,
                                                                                                 unsigned y,
                                                                                                 unsigned z) {
        KeyType sfc = 0;
        sfc |= interleave_3D(KeyType(x)) | interleave_3D(KeyType(y)) << 1 | interleave_3D(KeyType(z)) << 2;
        return sfc;
    }

    template <typename KeyType = Unsigned>
    MARS_INLINE_FUNCTION Octant decode_morton_2D(KeyType code) {
        return Octant(decode_morton_2DX(code), decode_morton_2DY(code));
    }

    template <typename KeyType = Unsigned>
    MARS_INLINE_FUNCTION Octant decode_morton_3D(KeyType code) {
        return Octant(decode_morton_3DX(code), decode_morton_3DY(code), decode_morton_3DZ(code));
    }

}  // namespace mars
#endif
