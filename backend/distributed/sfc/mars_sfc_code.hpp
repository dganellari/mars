
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


#include "mars_morton_code.hpp"
#include "mars_hilbert_code.hpp"

namespace mars {

//! @brief Strong type for Morton keys
template<class IntegerType>
using MortonKey = StrongType<IntegerType, struct MortonKeyTag>;

//! @brief Strong type for Hilbert keys
template<class IntegerType>
using HilbertKey = StrongType<IntegerType, struct HilbertKeyTag>;

//! @brief Meta function to detect Morton key types
template<class KeyType>
struct IsMorton : std::bool_constant<std::is_same_v<KeyType, MortonKey<typename KeyType::ValueType>>>
{
};

//! @brief Meta function to detect Hilbert key types
template<class KeyType>
struct IsHilbert : std::bool_constant<std::is_same_v<KeyType, HilbertKey<typename KeyType::ValueType>>>
{
};

//! @brief Key encode overload for Morton keys
template<class KeyType>
MARS_INLINE_FUNCTION std::enable_if_t<IsMorton<KeyType>{}, KeyType> iSfcKey(unsigned ix, unsigned iy, unsigned iz)
{
    return KeyType{iMorton<typename KeyType::ValueType>(ix, iy, iz)};
}

//! @brief Key encode overload for Hilbert keys
template<class KeyType>
MARS_INLINE_FUNCTION std::enable_if_t<IsHilbert<KeyType>{}, KeyType> iSfcKey(unsigned ix, unsigned iy, unsigned iz)
{
    return KeyType{iHilbert<typename KeyType::ValueType>(ix, iy, iz)};
}

//! @brief decode a Morton key
template<class KeyType>
MARS_INLINE_FUNCTION std::enable_if_t<IsMorton<KeyType>{}, util::tuple<unsigned, unsigned, unsigned>>
decodeSfc(KeyType key)
{
    return decodeMorton<typename KeyType::ValueType>(key);
}

//! @brief decode a Hilbert key
template<class KeyType>
MARS_INLINE_FUNCTION inline std::enable_if_t<IsHilbert<KeyType>{}, util::tuple<unsigned, unsigned, unsigned>>
decodeSfc(KeyType key)
{
    return decodeHilbert<typename KeyType::ValueType>(key);
}

}  // namespace mars
#endif
