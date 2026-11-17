#ifndef MARS_DIST_OCTANT_UTILS_HPP
#define MARS_DIST_OCTANT_UTILS_HPP

#include "mars_sfc_code.hpp"
#include "mars_distributed_octant.hpp"

namespace mars {

    template <Integer Type, class KeyType = MortonKey<Unsigned>>
    MARS_INLINE_FUNCTION Octant get_octant_from_sfc(const typename KeyType::ValueType gl_index) {
        Octant ref_octant;

        switch (Type) {
            case ElementType::Quad4: {
                ref_octant = decode_sfc_2D<KeyType>(gl_index);
                break;
            }
            case ElementType::Hex8: {
                ref_octant = decode_sfc_3D<KeyType>(gl_index);
                break;
            }
        }
        return ref_octant;
    }

    template <Integer Type, class KeyType = MortonKey<Unsigned>>
    MARS_INLINE_FUNCTION typename KeyType::ValueType get_sfc_from_octant(const Octant &o) {
        typename KeyType::ValueType enc_oc = 0;

        switch (Type) {
            case ElementType::Quad4: {
                enc_oc = encode_sfc_2D<KeyType>(o.x, o.y);
                break;
            }
            case ElementType::Hex8: {
                enc_oc = encode_sfc_3D<KeyType>(o.x, o.y, o.z);
                break;
            }
        }

        return enc_oc;
    }

    //-1 for all boundary. 0 left, 1 right, 2 down, 3 up and 4 and 5 for z dim.
    template <Integer Type, class KeyType = MortonKey<Unsigned>>
    MARS_INLINE_FUNCTION bool is_boundary_sfc(const typename KeyType::ValueType sfc,
                                              const int xdim,
                                              const int ydim,
                                              const int zdim,
                                              const Integer Face = -1) {
        Octant o = get_octant_from_sfc<Type, KeyType>(sfc);
        return o.is_boundary<Type>(xdim, ydim, zdim, Face);
    }

    inline Integer find_map_side(const std::map<std::string, Integer> &side_map, const std::string side) {
        auto search = side_map.find(side);
        if (side_map.end() != search) {
            return search->second;
        }

        std::cerr
            << "An invalid side named: " << side
            << "has been provided for the boundary iterator. Valid definitions are: all, left, right, bottom, top (2D) "
               "back and front (3D)!"
            << std::endl;
        std::exit(1);
    }

    // SFC based face numbering needs the following if top bottom face nr are differnet for 2D and 3D.
    template <Integer Type>
    std::enable_if_t<Type == ElementType::Quad4, Integer> map_side_to_value(const std::string side) {
        std::map<std::string, Integer> side_map{{"all", -1}, {"left", 0}, {"right", 1}, {"bottom", 2}, {"top", 3}};
        return find_map_side(side_map, side);
    }

    template <Integer Type>
    std::enable_if_t<Type == ElementType::Hex8, Integer> map_side_to_value(const std::string side) {
        std::map<std::string, Integer> side_map{
            {"all", -1}, {"left", 0}, {"right", 1}, {"bottom", 2}, {"top", 3}, {"front", 4}, {"back", 5}};
        return find_map_side(side_map, side);
    }

    template <Integer Type, class KeyType = MortonKey<Unsigned>>
    MARS_INLINE_FUNCTION void get_vertex_coordinates_from_sfc(const typename KeyType::ValueType sfc,
                                                              double *point,
                                                              const Integer xDim,
                                                              const Integer yDim,
                                                              const Integer zDim) {
        Octant o = get_octant_from_sfc<Type, KeyType>(sfc);
        o.get_vertex_coordinates<Type>(point, xDim, yDim, zDim);
    }

}  // namespace mars

#endif  // MARS_DIST_OCTANT_HPP
