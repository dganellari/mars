/*! @file
 * @brief  3D Hilbert encoding/decoding in 32- and 64-bit
 *
 * @author Daniel Ganellari <ganellari.d@gmail.com>
 *
 * The 2D Implementation based the Hilbert's Curve chapter from Hackers delight book:
 * https://learning.oreilly.com/library/view/hackers-delight-second/9780133084993.
 *
 * The 3D Hilbert implementation is based on the following paper:
 *
 * Yohei Miki, Masayuki Umemura
 * GOTHIC: Gravitational oct-tree code accelerated by hierarchical time step controlling
 * https://doi.org/10.1016/j.newast.2016.10.007
 *
 * implemented on the sph-exa code base: https://github.com/unibas-dmi-hpc/SPH-EXA.
 */

#ifndef MARS_HILBERT_CODE_HPP
#define MARS_HILBERT_CODE_HPP

#include "mars_base.hpp"
#include "mars_config.hpp"
#include "mars_globals.hpp"

namespace mars {

template<class KeyType>
struct maxTreeLevel
{
};

template<>
struct maxTreeLevel<unsigned> : stl::integral_constant<unsigned, 10>
{
};

template<>
struct maxTreeLevel<unsigned long long> : stl::integral_constant<unsigned, 21>
{
};
template<>
struct maxTreeLevel<Unsigned> : stl::integral_constant<unsigned, 21>
{
};

//! @brief maximum integer coordinate
template<class KeyType>
struct maxCoord : stl::integral_constant<unsigned, (1u << maxTreeLevel<KeyType>{})>
{
};

template<typename KeyType = Unsigned>
MARS_INLINE_FUNCTION
void decode_hilbert_2D(KeyType s, Unsigned *xp, Unsigned *yp) noexcept {


   auto order = maxTreeLevel<KeyType>{};

   s = s | (0x55555555 << 2 * order); // Pad s on left with 01
   const Unsigned sr = (s >> 1) & 0x55555555;  // (no change) groups.
   Unsigned cs = ((s & 0x55555555) + sr) // Compute complement &
        ^ 0x55555555;           // swap info in two-bit
                                // groups.
   // Parallel prefix xor op to propagate both complement
   // and swap info together from left to right (there is
   // no step "cs ^= cs >> 1", so in effect it computes
   // two independent parallel prefix operations on two
   // interleaved sets of sixteen bits).

   cs = cs ^ (cs >> 2);
   cs = cs ^ (cs >> 4);
   cs = cs ^ (cs >> 8);
   cs = cs ^ (cs >> 16);
   const Unsigned swap = cs & 0x55555555;      // Separate the swap and
   const Unsigned comp = (cs >> 1) & 0x55555555;  // complement bits.

   Unsigned t = (s & swap) ^ comp;       // Calculate x and y in
   s = s ^ sr ^ t ^ (t << 1);   // the odd & even bit
                                // positions, resp.
   s = s & ((1 << 2 * order) - 1);    // Clear out any junk
                                // on the left (unpad).

   // Now "unshuffle" to separate the x and y bits.

   t = (s ^ (s >> 1)) & 0x22222222; s = s ^ t ^ (t << 1);
   t = (s ^ (s >> 2)) & 0x0C0C0C0C; s = s ^ t ^ (t << 2);
   t = (s ^ (s >> 4)) & 0x00F000F0; s = s ^ t ^ (t << 4);
   t = (s ^ (s >> 8)) & 0x0000FF00; s = s ^ t ^ (t << 8);

   *xp = s >> 16;              // Assign the two halves
   *yp = s & 0xFFFF;           // of t to x and y.
}

template<typename KeyType = Unsigned>
MARS_INLINE_FUNCTION
std::enable_if_t<std::is_unsigned_v<KeyType>, KeyType> encode_hilbert_2D(Unsigned x, Unsigned y) noexcept {

   int i, xi, yi;
   Unsigned temp;
   KeyType s = 0;

   for (i = maxTreeLevel<KeyType>{} - 1; i >= 0; i--) {
      xi = (x >> i) & 1;            // Get bit i of x.
      yi = (y >> i) & 1;           // Get bit i of y.

      if (yi == 0) {
         temp = x;                 // Swap x and y and,
         x = y^(-xi);              // if xi = 1,
         y = temp^(-xi);           // complement them.
      }
      s = 4*s + 2*xi + (xi^yi);    // Append two bits to s.
   }
   return s;
}

template<typename KeyType = Unsigned>
MARS_INLINE_FUNCTION
std::enable_if_t<std::is_unsigned_v<KeyType>, KeyType> hilbert_s_from_xy(Unsigned x, Unsigned y) {

   int i;
   Unsigned state, s, row;

   state = 0;                         // Initialize.
   s = 0;

   for (i = maxTreeLevel<KeyType>{} - 1; i >= 0; i--) {
      row = 4*state | 2*((x >> i) & 1) | (y >> i) & 1;
      s = (s << 2) | (0x361E9CB4 >> 2*row) & 3;
      state = (0x8FE65831 >> 2*row) & 3;
   }
   return s;
}

//---------------3D------------------------
//
template<typename KeyType = Unsigned>
MARS_INLINE_FUNCTION
Octant decode_hilbert_3D(KeyType key) noexcept
{
    Unsigned px = 0;
    Unsigned py = 0;
    Unsigned pz = 0;

    for (Unsigned level = 0; level < maxTreeLevel<KeyType>; ++level)
    {
        Unsigned octant   = (key >> (3 * level)) & 7u;
        const Unsigned xi = octant >> 2u;
        const Unsigned yi = (octant >> 1u) & 1u;
        const Unsigned zi = octant & 1u;

        if (yi ^ zi)
        {
            // cyclic rotation
            Unsigned pt = px;
            px          = pz;
            pz          = py;
            py          = pt;
        }
        else if ((!xi & !yi & !zi) || (xi & yi & zi))
        {
            // swap x and z
            Unsigned pt = px;
            px          = pz;
            pz          = pt;
        }

        // turn px, py and pz
        Unsigned mask = (1 << level) - 1;
        px ^= mask & (-(xi & (yi | zi)));
        py ^= mask & (-((xi & ((!yi) | (!zi))) | ((!xi) & yi & zi)));
        pz ^= mask & (-((xi & (!yi) & (!zi)) | (yi & zi)));

        // append 1 bit to the positions
        px |= (xi << level);
        py |= ((xi ^ yi) << level);
        pz |= ((yi ^ zi) << level);
    }

    return Octant(px, py, pz);
}


#if defined(__CUDACC__) || defined(__HIPCC__)
__device__ static unsigned mortonToHilbertDevice[8] = {0, 1, 3, 2, 7, 6, 4, 5};
#endif

/*! @brief compute the Hilbert key for a 3D point of integer coordinates
 *
 * @tparam     KeyType   32- or 64-bit unsigned integer
 * @param[in]  px,py,pz  input coordinates in [0:2^maxTreeLevel<KeyType>{}]
 * @return               the Hilbert key
 */
MARS_INLINE_FUNCTION
encode_hilbert_3D(unsigned px, unsigned py, unsigned pz) noexcept
{
    assert(px < (1u << maxTreeLevel<KeyType>{}));
    assert(py < (1u << maxTreeLevel<KeyType>{}));
    assert(pz < (1u << maxTreeLevel<KeyType>{}));

#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
    constexpr unsigned mortonToHilbert[8] = {0, 1, 3, 2, 7, 6, 4, 5};
#endif

    KeyType key = 0;

    for (int level = maxTreeLevel<KeyType>{} - 1; level >= 0; --level)
    {
        unsigned xi = (px >> level) & 1u;
        unsigned yi = (py >> level) & 1u;
        unsigned zi = (pz >> level) & 1u;

        // append 3 bits to the key
        unsigned octant = (xi << 2) | (yi << 1) | zi;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        key = (key << 3) + mortonToHilbertDevice[octant];
#else
        key = (key << 3) + mortonToHilbert[octant];
#endif

        // turn px, py and pz
        px ^= -(xi & ((!yi) | zi));
        py ^= -((xi & (yi | zi)) | (yi & (!zi)));
        pz ^= -((xi & (!yi) & (!zi)) | (yi & (!zi)));

        if (zi)
        {
            // cyclic rotation
            unsigned pt = px;
            px          = py;
            py          = pz;
            pz          = pt;
        }
        else if (!yi)
        {
            // swap x and z
            unsigned pt = px;
            px          = pz;
            pz          = pt;
        }
    }

    return key;
}

}  // namespace mars
#endif
