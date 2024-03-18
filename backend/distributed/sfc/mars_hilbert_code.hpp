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

#include "mars_distributed_octant.hpp"
#include "mars_sfc_utils.hpp"

namespace mars {

    //! @brief inverse function of iHilbert_2D 32 bit only up to oder 16 but works at constant time.
    MARS_INLINE_FUNCTION Octant decode_hilbert_2D(unsigned s) noexcept {
        unsigned order = maxTreeLevel<unsigned>{};

        s = s | (0x55555555 << 2 * order);  // Pad s on left with 01

        const unsigned sr = (s >> 1) & 0x55555555;           // (no change) groups.
        unsigned cs = ((s & 0x55555555) + sr) ^ 0x55555555;  // Compute complement & swap info in two-bit groups.
        // Parallel prefix xor op to propagate both complement
        // and swap info together from left to right (there is
        // no step "cs ^= cs >> 1", so in effect it computes
        // two independent parallel prefix operations on two
        // interleaved sets of sixteen bits).
        cs = cs ^ (cs >> 2);
        cs = cs ^ (cs >> 4);
        cs = cs ^ (cs >> 8);
        cs = cs ^ (cs >> 16);
        const unsigned swap = cs & 0x55555555;         // Separate the swap and
        const unsigned comp = (cs >> 1) & 0x55555555;  // complement bits.

        unsigned t = (s & swap) ^ comp;  // Calculate x and y in
        s = s ^ sr ^ t ^ (t << 1);       // the odd & even bit positions, resp.
        s = s & ((1 << 2 * order) - 1);  // Clear out any junk on the left (unpad).

        // Now "unshuffle" to separate the x and y bits.

        t = (s ^ (s >> 1)) & 0x22222222;
        s = s ^ t ^ (t << 1);
        t = (s ^ (s >> 2)) & 0x0C0C0C0C;
        s = s ^ t ^ (t << 2);
        t = (s ^ (s >> 4)) & 0x00F000F0;
        s = s ^ t ^ (t << 4);
        t = (s ^ (s >> 8)) & 0x0000FF00;
        s = s ^ t ^ (t << 8);

        unsigned xp = s >> 16;     // Assign the two halves
        unsigned yp = s & 0xFFFF;  // of t to x and y.

        return Octant(xp, yp);
    }

    // Lam and Shapiro inverse function of hilbert
    template <class KeyType>
    MARS_INLINE_FUNCTION Octant decode_hilbert_2D(KeyType key) noexcept {
        unsigned sa, sb;
        unsigned x = 0, y = 0, temp = 0;
        unsigned order = maxTreeLevel<KeyType>{};

        for (unsigned level = 0; level < 2 * order; level += 2) {
            // Get bit level+1 of key.
            sa = (key >> (level + 1)) & 1;
            // Get bit level of key.
            sb = (key >> level) & 1;
            if ((sa ^ sb) == 0) {
                // If sa,sb = 00 or 11,
                temp = x;
                // swap x and y,
                x = y ^ (-sa);
                // and if sa = 1,
                y = temp ^ (-sa);
                // complement them.
            }
            x = (x >> 1) | (sa << 31);         // Prepend sa to x and
            y = (y >> 1) | ((sa ^ sb) << 31);  // (sa ^ sb) to y.
        }
        unsigned px = x >> (32 - order);
        // Right-adjust x and y
        unsigned py = y >> (32 - order);
        // and return them to
        return Octant(px, py);
    }

    template <typename KeyType = Unsigned>
    MARS_INLINE_FUNCTION std::enable_if_t<std::is_unsigned_v<KeyType>, KeyType> encode_hilbert_2D(unsigned x,
                                                                                                  unsigned y) noexcept {
        assert(x < (1u << maxTreeLevel<KeyType>{}));
        assert(y < (1u << maxTreeLevel<KeyType>{}));

        unsigned xi, yi;
        Unsigned temp;
        KeyType s = 0;

        for (int i = maxTreeLevel<KeyType>{} - 1; i >= 0; i--) {
            xi = (x >> i) & 1;  // Get bit i of x.
            yi = (y >> i) & 1;  // Get bit i of y.

            if (yi == 0) {
                temp = x;          // Swap x and y and,
                x = y ^ (-xi);     // if xi = 1,
                y = temp ^ (-xi);  // complement them.
            }
            s = 4 * s + 2 * xi + (xi ^ yi);  // Append two bits to s.
        }
        return s;
    }

    template <typename KeyType = Unsigned>
    MARS_INLINE_FUNCTION std::enable_if_t<std::is_unsigned_v<KeyType>, KeyType> hilbert_s_from_xy(Unsigned x,
                                                                                                  Unsigned y) {
        int i;
        Unsigned state, s, row;

        state = 0;  // Initialize.
        s = 0;

        for (i = maxTreeLevel<KeyType>{} - 1; i >= 0; i--) {
            row = 4 * state | 2 * ((x >> i) & 1) | (y >> i) & 1;
            s = (s << 2) | (0x361E9CB4 >> 2 * row) & 3;
            state = (0x8FE65831 >> 2 * row) & 3;
        }
        return s;
    }

    //---------------3D------------------------
    //
    template <typename KeyType = Unsigned>
    MARS_INLINE_FUNCTION Octant decode_hilbert_3D(KeyType key) noexcept {
        unsigned px = 0;
        unsigned py = 0;
        unsigned pz = 0;

        for (unsigned level = 0; level < maxTreeLevel<KeyType>{}; ++level) {
            unsigned octant = (key >> (3 * level)) & 7u;
            const unsigned xi = octant >> 2u;
            const unsigned yi = (octant >> 1u) & 1u;
            const unsigned zi = octant & 1u;

            if (yi ^ zi) {
                // cyclic rotation
                unsigned pt = px;
                px = pz;
                pz = py;
                py = pt;
            } else if ((!xi & !yi & !zi) || (xi & yi & zi)) {
                // swap x and z
                unsigned pt = px;
                px = pz;
                pz = pt;
            }

            // turn px, py and pz
            unsigned mask = (1 << level) - 1;
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
    template <typename KeyType = Unsigned>
    MARS_INLINE_FUNCTION std::enable_if_t<std::is_unsigned_v<KeyType>, KeyType>
    encode_hilbert_3D(unsigned px, unsigned py, unsigned pz) noexcept {
        assert(px < (1u << maxTreeLevel<KeyType>{}));
        assert(py < (1u << maxTreeLevel<KeyType>{}));
        assert(pz < (1u << maxTreeLevel<KeyType>{}));

#if !defined(__CUDA_ARCH__) && !defined(__HIP_DEVICE_COMPILE__)
        constexpr unsigned mortonToHilbert[8] = {0, 1, 3, 2, 7, 6, 4, 5};
#endif

        KeyType key = 0;

        for (int level = maxTreeLevel<KeyType>{} - 1; level >= 0; --level) {
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

            if (zi) {
                // cyclic rotation
                unsigned pt = px;
                px = py;
                py = pz;
                pz = pt;
            } else if (!yi) {
                // swap x and z
                unsigned pt = px;
                px = pz;
                pz = pt;
            }
        }

        return key;
    }

}  // namespace mars
#endif
