/*
 * MIT License
 *
 * Copyright (c) 2021 CSCS, ETH Zurich
 *               2021 University of Basel
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef MARS_SFC_UTILS_HPP
#define MARS_SFC_UTILS_HPP


#include "mars_base.hpp"
#include "mars_config.hpp"
#include "mars_globals.hpp"

namespace mars {

template<class KeyType>
struct maxTreeLevel
{
};

template <>
struct maxTreeLevel<unsigned> : integral_constant<unsigned, 10> {};

template <>
struct maxTreeLevel<unsigned long> : integral_constant<unsigned, 21> {};

template <>
struct maxTreeLevel<unsigned long long> : integral_constant<unsigned, 21> {};

template<>
struct maxTreeLevel<int> : integral_constant<unsigned, 10> {};

//! @brief maximum integer coordinate
template <class KeyType>
struct maxCoord : integral_constant<unsigned, (1u << maxTreeLevel<KeyType>{})> {};

/*! @brief add (binary) zeros behind a prefix
 *
 * Allows comparisons, such as number of leading common bits (cpr)
 * of the prefix with SFC codes.
 *
 * @tparam KeyType  32- or 64-bit unsigned integer type
 * @param prefix    the bit pattern
 * @param length    number of bits in the prefix
 * @return          prefix padded out with zeros
 *
 * Examples:
 *  pad(0b011u,  3) -> 0b00011 << 27
 *  pad(0b011ul, 3) -> 0b0011ul << 60
 *
 *  i.e. @p length plus the number of zeros added adds up to 30 for 32-bit integers
 *  or 63 for 64-bit integers, because these are the numbers of usable bits in SFC codes.
 */
template<class KeyType>
constexpr KeyType pad(KeyType prefix, int length)
{
    return prefix << (3 * maxTreeLevel<KeyType>{} - length);
}

/*! @brief extract the n-th octal digit from an SFC key, starting from the most significant
 *
 * @tparam KeyType   32- or 64-bit unsigned integer type
 * @param code       Input SFC key code
 * @param position   Which digit place to extract. Return values will be meaningful for
 *                   @p position in [1:11] for 32-bit keys and in [1:22] for 64-bit keys and
 *                   will be zero otherwise, but a value of 0 for @p position can also be specified
 *                   to detect whether the 31st or 63rd bit for the last cornerstone is non-zero.
 *                   (The last cornerstone has a value of nodeRange<KeyType>(0) = 2^31 or 2^63)
 * @return           The value of the digit at place @p position
 *
 * The position argument correspondence to octal digit places has been chosen such that
 * octalDigit(code, pos) returns the octant at octree division level pos.
 */
template<class KeyType>
MARS_INLINE_FUNCTION constexpr unsigned octalDigit(KeyType code, unsigned position)
{
    return (code >> (3u * (maxTreeLevel<KeyType>{} - position))) & 7u;
}

/*! @brief compute the maximum range of an octree node at a given subdivision level
 *
 * @tparam KeyType    32- or 64-bit unsigned integer type
 * @param  treeLevel  octree subdivision level
 * @return            the range
 *
 * At treeLevel 0, the range is the entire 30 or 63 bits used in the SFC code.
 * After that, the range decreases by 3 bits for each level.
 *
 */
template<class KeyType>
MARS_INLINE_FUNCTION constexpr KeyType nodeRange(unsigned treeLevel)
{
    assert(treeLevel <= maxTreeLevel<KeyType>{});
    unsigned shifts = maxTreeLevel<KeyType>{} - treeLevel;

    return KeyType(1ul << (3u * shifts));
}


}  // namespace mars
#endif
