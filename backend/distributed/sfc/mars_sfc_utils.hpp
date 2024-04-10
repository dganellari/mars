
#ifndef MARS_SFC_UTILS_HPP
#define MARS_SFC_UTILS_HPP

/*! @file
 * @brief  SFC encoding/decoding utilities and definitions
 *
 * @author Daniel Ganellari
 */

#include <_types/_uint64_t.h>
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
struct maxTreeLevel<uint64_t> : integral_constant<unsigned, 21> {};

template <>
struct maxTreeLevel<unsigned long> : integral_constant<unsigned, 21> {};

//! @brief maximum integer coordinate
template <class KeyType>
struct maxCoord : integral_constant<unsigned, (1u << maxTreeLevel<KeyType>{})> {};

}  // namespace mars
#endif
