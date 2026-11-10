#pragma once
#include "backend/distributed/unstructured/domain.hpp"
namespace mars
{
namespace fem
{
template<class Tag>
struct ElemTraits;
template<>
struct ElemTraits<TetTag>
{
    static constexpr int NodesPerElem = 4;
};
template<>
struct ElemTraits<HexTag>
{
    static constexpr int NodesPerElem = 8;
};
} // namespace fem
} // namespace mars
