
#ifndef MARS_SFC_CODE_HPP
#define MARS_SFC_CODE_HPP

/*! @file
 * @brief  SFC encoding/decoding in 32- and 64-bit
 *
 * @author Daniel Ganellari
 *
 * Common interface to Morton and Hilbert keys based on strong C++ types
 * based on the https://github.com/unibas-dmi-hpc/SPH-EXA
 */

#include "mars_hilbert_code.hpp"
#include "mars_morton_code.hpp"

namespace mars {

    //! @brief Strong type for Morton keys
    template <class IntegerType>
    using MortonKey = StrongType<IntegerType, struct MortonKeyTag>;

    //! @brief Strong type for Hilbert keys
    template <class IntegerType>
    using HilbertKey = StrongType<IntegerType, struct HilbertKeyTag>;

    //! @brief Meta function to detect Morton key types
    template <class KeyType>
    struct IsMorton : std::bool_constant<std::is_same_v<KeyType, MortonKey<typename KeyType::ValueType>>> {};

    //! @brief Meta function to detect Hilbert key types
    template <class KeyType>
    struct IsHilbert : std::bool_constant<std::is_same_v<KeyType, HilbertKey<typename KeyType::ValueType>>> {};

    //! @brief Key encode overload for Morton keys
    template <class KeyType>
    MARS_INLINE_FUNCTION std::enable_if_t<IsMorton<KeyType>{}, KeyType> encode_sfc_2D(unsigned ix, unsigned iy) {
        return KeyType{encode_morton_2D<typename KeyType::ValueType>(ix, iy)};
    }

    //! @brief Key encode overload for Hilbert keys
    template <class KeyType>
    MARS_INLINE_FUNCTION std::enable_if_t<IsHilbert<KeyType>{}, KeyType> encode_sfc_2D(unsigned ix, unsigned iy) {
        return KeyType{encode_hilbert_2D<typename KeyType::ValueType>(ix, iy)};
    }

    //! @brief decode a Morton key
    template <class KeyType>
    MARS_INLINE_FUNCTION std::enable_if_t<IsMorton<KeyType>{}, Octant> decode_sfc_2D(KeyType key) {
        return decode_morton_2D(key);
    }

    //! @brief decode a Hilbert key
    template <class KeyType>
    MARS_INLINE_FUNCTION std::enable_if_t<IsHilbert<KeyType>{}, Octant> decode_sfc_2D(KeyType key) {
        return decode_hilbert_2D<typename KeyType::ValueType>(key);
    }

}  // namespace mars
#endif
