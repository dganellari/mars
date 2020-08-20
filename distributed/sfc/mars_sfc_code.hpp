#ifndef MARS_SFC_CODE_HPP
#define MARS_SFC_CODE_HPP

namespace mars
{

/* Magicbitcode morton z curve from fgiesen.wordpress.com with slight modifications*/

// "Insert" a 0 bit after each of the 16 low bits of x
MARS_INLINE_FUNCTION
unsigned int interleave_by1(unsigned int x)
{
    x &= 0x0000ffff;                 // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x << 8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x << 4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x << 2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x << 1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    return x;
}

// "Insert" two 0 bits after each of the 10 low bits of x
MARS_INLINE_FUNCTION
unsigned int interleave_by2(unsigned int x)
{
    x &= 0x000003ff;                  // x = ---- ---- ---- ---- ---- --98 7654 3210
    x = (x ^ (x << 16)) & 0xff0000ff; // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x ^ (x << 8)) & 0x0300f00f;  // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x ^ (x << 4)) & 0x030c30c3;  // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x ^ (x << 2)) & 0x09249249;  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    return x;
}

MARS_INLINE_FUNCTION
unsigned int encode_morton_2D(unsigned int x, unsigned int y)
{
    return (interleave_by1(y) << 1) | interleave_by1(x);
}

MARS_INLINE_FUNCTION
unsigned int encode_morton_3D(unsigned int x, unsigned int y, unsigned int z)
{
    return (interleave_by2(z) << 2) | (interleave_by2(y) << 1) | interleave_by2(x);
}

// Inverse of Part1By1 - "delete" all odd-indexed bits
MARS_INLINE_FUNCTION
unsigned int recall_by1(unsigned int x)
{
    x &= 0x55555555;                 // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x = (x ^ (x >> 1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >> 2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >> 4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >> 8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
    return x;
}

// Inverse of Part1By2 - "delete" all bits not at positions divisible by 3
MARS_INLINE_FUNCTION
unsigned int recall_by2(unsigned int x)
{
    x &= 0x09249249;                  // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    x = (x ^ (x >> 2)) & 0x030c30c3;  // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x ^ (x >> 4)) & 0x0300f00f;  // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x ^ (x >> 8)) & 0xff0000ff;  // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x ^ (x >> 16)) & 0x000003ff; // x = ---- ---- ---- ---- ---- --98 7654 3210
    return x;
}

MARS_INLINE_FUNCTION
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
}

} // namespace mars
#endif 