#pragma once

#include <cmath>
#include <cstone/sfc/sfc.hpp>

namespace mars::amr
{
// The point an eighth of a cell inside the cell of cornerstone's encode (sfc3D) that holds a key: it encodes back to
// the key, and the edge midpoints and face centres of such points encode the same in any summation order. The
// decoded node coordinates (decodeSfcToPhysical) do not always encode back. Non-conforming refinement computes the
// new nodes from these points and snaps them the same way, in a box kept from the first mesh, so a node that a later
// level emits again (the coarse side of a hanging node) gets the key and position it had.
template<typename KeyType, typename RealType>
struct KeyedPosition
{
    cstone::Box<RealType> box;

    HOST_DEVICE_FUN static RealType axis(unsigned i, unsigned bits, RealType lo, RealType inverseLength)
    {
        const RealType m = RealType(1u << bits) * inverseLength; // the scale of sfc3D
        return (RealType(i) + ceil(lo * m) + RealType(0.125)) / m;
    }

    HOST_DEVICE_FUN void operator()(KeyType key, RealType& x, RealType& y, RealType& z) const
    {
        const auto bits   = box.getBoxDimBits(cstone::maxTreeLevel<KeyType>{});
        auto [ix, iy, iz] = cstone::decodeSfc(cstone::sfcKey(key), bits);
        x                 = axis(ix, bits[0], box.xmin(), box.ilx());
        y                 = axis(iy, bits[1], box.ymin(), box.ily());
        z                 = axis(iz, bits[2], box.zmin(), box.ilz());
    }

    HOST_DEVICE_FUN KeyType encode(RealType x, RealType y, RealType z) const
    {
        return cstone::sfc3D<cstone::SfcKind<KeyType>>(x, y, z, box);
    }
};

} // namespace mars::amr
