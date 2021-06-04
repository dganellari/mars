#ifndef MARS_DIST_OCTANT_HPP
#define MARS_DIST_OCTANT_HPP

#include <map>
#include <type_traits>
#include "mars_base.hpp"
#include "mars_globals.hpp"
#include "mars_sfc_code.hpp"

namespace mars {
    struct Octant {
        Integer x;
        Integer y;
        Integer z;

        bool valid;
        bool share_boundary;

        Octant() = default;

        MARS_INLINE_FUNCTION
        Octant(Integer i, Integer j) : x(i), y(j), z(0) {
            valid = true;
            share_boundary = false;
        }

        MARS_INLINE_FUNCTION
        Octant(Integer i, Integer j, Integer k) : x(i), y(j), z(k) {
            valid = true;
            share_boundary = false;
        }

        MARS_INLINE_FUNCTION
        bool shares_boundary() { return share_boundary; }

        MARS_INLINE_FUNCTION
        void set_shares_boundary() { share_boundary = true; }

        MARS_INLINE_FUNCTION
        bool is_valid() const { return valid; }

        MARS_INLINE_FUNCTION
        void set_invalid() { valid = false; }

        MARS_INLINE_FUNCTION
        bool is_equals(Octant o) const { return (x == o.x && y == o.y && z == o.z); }

        template <Integer Type>
        MARS_INLINE_FUNCTION bool is_boundary(const int xdim, const int ydim, const int zdim, const Integer Face = -1) {
            switch (Face) {
                case 0: {
                    return (x == 0);
                }
                case 1: {
                    return (x == xdim);
                }
                case 2: {
                    return (y == 0);
                }
                case 3: {
                    return (y == ydim);
                }
                case 4: {
                    return (z == 0);
                }
                case 5: {
                    return (z == zdim);
                }
                case -1: {
                    return boundary<Type>(xdim, ydim, zdim);
                }
                default: {
                    printf("Invalid face number template argument!\n");
                    return false;
                }
            }
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION bool boundary(const int xdim, const int ydim, const int zdim) {
            switch (Type) {
                case ElementType::Quad4: {
                    return (x == 0 || y == 0 || x == xdim || y == ydim);
                }
                case ElementType::Hex8: {
                    return (x == 0 || y == 0 || z == 0 || x == xdim || y == ydim || z == zdim);
                }
                default: {
                    printf("the element type is not valid\n");
                    return false;
                }
            }
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION Integer get_global_index(const int xDim, const int yDim) {
            switch (Type) {
                case ElementType::Quad4: {
                    return x + (xDim + 1) * y;
                }
                case ElementType::Hex8: {
                    return elem_index(x, y, z, xDim, yDim);
                }
                default: {
                    printf("The element type is not valid\n");
                    return -1;
                }
            }
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION void get_vertex_coordinates(double *point,
                                                         const Integer xDim,
                                                         const Integer yDim,
                                                         const Integer zDim) const {
            switch (Type) {
                case ElementType::Quad4: {
                    point[0] = static_cast<Real>(x) / static_cast<Real>(xDim);
                    point[1] = static_cast<Real>(y) / static_cast<Real>(yDim);
                    /* point[2] = 0.; */
                    break;
                }
                case ElementType::Hex8: {
                    point[0] = static_cast<Real>(x) / static_cast<Real>(xDim);
                    point[1] = static_cast<Real>(y) / static_cast<Real>(yDim);
                    point[2] = static_cast<Real>(z) / static_cast<Real>(zDim);
                    break;
                }
            }
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION void validate_dof_nbh(const Integer xDim,
                                                   const Integer yDim,
                                                   const Integer zDim,
                                                   const bool periodic) {
            switch (Type) {
                case ElementType::Quad4: {
                    if (x < 0 || y < 0 || x > xDim || y > yDim) {
                        set_invalid();
                    }
                    break;
                }
                case ElementType::Hex8: {
                    if (x < 0 || y < 0 || z < 0 || x > xDim || y > yDim || z > zDim) {
                        set_invalid();
                    }
                    break;
                }
                default: {
                    printf("Could not validate: The element type is not valid\n");
                    break;
                }
            }
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION void validate_nbh(const Integer xDim,
                                               const Integer yDim,
                                               const Integer zDim,
                                               const bool periodic) {
            switch (Type) {
                case ElementType::Quad4: {
                    if (x < 0 || y < 0 || x >= xDim || y >= yDim) {
                        // check only xdim and ydim and not for -1 because we do not want to process the same face twice
                        if (periodic) {
                            x = (xDim + x) % xDim;
                            y = (yDim + y) % yDim;
                        } else {
                            set_invalid();
                        }
                        set_shares_boundary();
                    }
                    break;
                }
                case ElementType::Hex8: {
                    if (x < 0 || y < 0 || z < 0 || x >= xDim || y >= yDim || z >= zDim) {
                        if (periodic) {
                            x = (xDim + x) % xDim;
                            y = (yDim + y) % yDim;
                            z = (zDim + z) % zDim;
                        } else {
                            set_invalid();
                        }
                        set_shares_boundary();
                    }
                    break;
                }
                default: {
                    printf("Could not validate: The element type is not valid\n");
                    break;
                }
            }
        }

        MARS_INLINE_FUNCTION Integer get_edge_direction(const int edge) const { return edge / 4; }

        MARS_INLINE_FUNCTION Octant sfc_edge_start(const int edge) const {
            assert(0 <= edge && edge < 12);
            auto direction = get_edge_direction(edge);
            auto e = edge % 4;
            Integer x_ = -1, y_ = -1, z_ = -1;

            switch (direction) {
                case 0:
                    x_ = x;
                    y_ = y + (e & 1);
                    z_ = z + e / 2;
                    break;
                case 1:
                    x_ = x + (e & 1);
                    y_ = y;
                    z_ = z + e / 2;
                    break;
                case 2:
                    x_ = x + (e & 1);
                    y_ = y + e / 2;
                    z_ = z;
                    break;
                default:
                    printf("The element type is not valid\n");
                    break;
            }
            return Octant(x_, y_, z_);
        }

        // Only a 3D functionality defined for hex8 elements only.
        template <Integer Type>
        MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Octant> edge_start(const int edge,
                                                                                            const Integer xDim,
                                                                                            const Integer yDim,
                                                                                            const Integer zDim,
                                                                                            const bool periodic) const {
            Octant nbh = sfc_edge_start(edge);
            nbh.validate_nbh<Type>(xDim, yDim, zDim, periodic);
            return nbh;
        }
        MARS_INLINE_FUNCTION Octant sfc_edge_nbh(const int edge) const {
            // adapted from the p4est corner neighbor for the mesh generation
            assert(0 <= edge && edge < 12);
            auto direction = get_edge_direction(edge);
            Integer x_ = -1, y_ = -1, z_ = -1;

            switch (direction) {
                case 0:
                    x_ = x;
                    y_ = y + (2 * (edge & 1) - 1);
                    z_ = z + ((edge & 2) - 1);
                    break;
                case 1:
                    x_ = x + (2 * (edge & 1) - 1);
                    y_ = y;
                    z_ = z + ((edge & 2) - 1);
                    break;
                case 2:
                    x_ = x + (2 * (edge & 1) - 1);
                    y_ = y + ((edge & 2) - 1);
                    z_ = z;
                    break;
                default:
                    printf("The element type is not valid\n");
                    break;
            }
            return Octant(x_, y_, z_);
        }

        // Only a 3D functionality defined for hex8 elements only.
        template <Integer Type>
        MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Octant> edge_nbh(const int edge,
                                                                                          const Integer xDim,
                                                                                          const Integer yDim,
                                                                                          const Integer zDim,
                                                                                          const bool periodic) const {
            Octant nbh = sfc_edge_nbh(edge);
            nbh.validate_nbh<Type>(xDim, yDim, zDim, periodic);
            return nbh;
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION Octant sfc_face_nbh(const int face) const {
            switch (Type) {
                case ElementType::Quad4: {
                    // adapted from the p4est corner neighbor for the mesh generation
                    const Integer x_ = x + ((face == 0) ? -1 : (face == 1) ? 1 : 0);
                    const Integer y_ = y + ((face == 2) ? -1 : (face == 3) ? 1 : 0);

                    return Octant(x_, y_);
                }
                case ElementType::Hex8: {
                    // adapted from the p4est corner neighbor for the mesh generation
                    const Integer x_ = x + ((face == 0) ? -1 : (face == 1) ? 1 : 0);
                    const Integer y_ = y + ((face == 2) ? -1 : (face == 3) ? 1 : 0);
                    const Integer z_ = z + ((face == 4) ? -1 : (face == 5) ? 1 : 0);

                    return Octant(x_, y_, z_);
                }
                default: {
                    printf("The element type is not valid\n");
                    return Octant();
                }
            }
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION Octant face_nbh(const int face,
                                             const Integer xDim,
                                             const Integer yDim,
                                             const Integer zDim,
                                             const bool periodic) const {
            Octant nbh = sfc_face_nbh<Type>(face);
            nbh.validate_nbh<Type>(xDim, yDim, zDim, periodic);
            return nbh;
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION Octant sfc_corner_nbh(const int corner) const {
            switch (Type) {
                case ElementType::Quad4: {
                    // adapted from the p4est corner neighbor for the mesh generation
                    const Integer x_ = x + 2 * (corner & 1) - 1;
                    const Integer y_ = y + (corner & 2) - 1;

                    return Octant(x_, y_);
                }
                case ElementType::Hex8: {
                    // adapted from the p4est corner neighbor for the mesh generation
                    const Integer x_ = x + 2 * (corner & 1) - 1;
                    const Integer y_ = y + (corner & 2) - 1;
                    const Integer z_ = z + (corner & 4) / 2 - 1;

                    return Octant(x_, y_, z_);
                }
                default: {
                    printf("The element type is not valid\n");
                    return Octant();
                }
            }
        }

        template <Integer Type>
        MARS_INLINE_FUNCTION Octant corner_nbh(const int corner,
                                               const Integer xDim,
                                               const Integer yDim,
                                               const Integer zDim,
                                               const bool periodic) const {
            Octant nbh = sfc_corner_nbh<Type>(corner);
            nbh.validate_nbh<Type>(xDim, yDim, zDim, periodic);
            return nbh;
        }

        struct Depth {
            Integer x;
            Integer y;
            Integer z;

            Depth(Integer xd, Integer yd, Integer zd) : x(xd), y(yd), z(zd) {}

            // This gives a depth one one ring neighbors since the largest you are substracting from the x is 1.
            template <Integer Type>
            static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, Depth> set_depth(
                const Integer depth) {
                Depth d(depth, depth, depth);
                return d;
            }

            // This is the 2D case. Meaning that the z coordinate is skipped.
            template <Integer Type>
            static MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Quad4, Depth> set_depth(
                const Integer depth) {
                Depth d(depth, depth, 1);
                return d;
            }
        };

        template <Integer Type>
        MARS_INLINE_FUNCTION void one_ring_corner_nbhs(Octant one_ring[Type],
                                                       const Integer xDim,
                                                       const Integer yDim,
                                                       const Integer zDim,
                                                       const bool periodic) const {
            const Integer ring_depth = 2;
            auto depth = Depth::set_depth<Type>(ring_depth);

            for (int k = 0; k < depth.z; k++) {
                for (int j = 0; j < depth.y; j++) {
                    for (int i = 0; i < depth.x; i++) {
                        Octant nbh(x - i, y - j, z - k);
                        nbh.validate_nbh<Type>(xDim, yDim, zDim, periodic);
                        one_ring[k * depth.y * depth.x + j * depth.x + i] = nbh;
                    }
                }
            }
        }

        template <Integer Type, typename F>
        MARS_INLINE_FUNCTION void one_ring_corner_nbhs(F f,
                                                       const Integer xDim,
                                                       const Integer yDim,
                                                       const Integer zDim,
                                                       const bool periodic) const {
            const Integer ring_depth = 2;
            auto depth = Depth::set_depth<Type>(ring_depth);

            for (int k = 0; k < depth.z; k++) {
                for (int j = 0; j < depth.y; j++) {
                    for (int i = 0; i < depth.x; i++) {
                        Octant nbh(x - i, y - j, z - k);
                        nbh.validate_nbh<Type>(xDim, yDim, zDim, periodic);
                        f(nbh);
                    }
                }
            }
        }

        template <Integer Type, typename F>
        MARS_INLINE_FUNCTION std::enable_if_t<Type == ElementType::Hex8, void> one_ring_edge_nbhs(
            F f,
            const Integer direction,
            const Integer xDim,
            const Integer yDim,
            const Integer zDim,
            const bool periodic) const {
            const Integer ring_depth = 2;
            auto depth = Depth::set_depth<Type>(ring_depth);

            for (int j = 0; j < depth.y; j++) {
                for (int i = 0; i < depth.x; i++) {
                    Octant nbh = get_edge_ring_nbh(i, j, direction);
                    nbh.validate_nbh<Type>(xDim, yDim, zDim, periodic);
                    f(nbh);
                }
            }
        }

        MARS_INLINE_FUNCTION Octant get_edge_ring_nbh(const int i, const int j, const Integer direction) const {
            Integer x_ = -1, y_ = -1, z_ = -1;

            switch (direction) {
                case 0:
                    x_ = x;
                    y_ = y - i;
                    z_ = z - j;
                    break;
                case 1:
                    x_ = x - i;
                    y_ = y;
                    z_ = z - j;
                    break;
                case 2:
                    x_ = x - i;
                    y_ = y - j;
                    z_ = z;
                    break;
                default:
                    printf("The element type is not valid\n");
                    break;
            }
            return Octant(x_, y_, z_);
        }
    };

    template <Integer Type>
    MARS_INLINE_FUNCTION Octant get_octant_from_sfc(const Integer gl_index) {
        Octant ref_octant;

        switch (Type) {
            case ElementType::Quad4: {
                const Integer i = decode_morton_2DX(gl_index);
                const Integer j = decode_morton_2DY(gl_index);

                ref_octant = Octant(i, j);
                break;
            }
            case ElementType::Hex8: {
                const Integer i = decode_morton_3DX(gl_index);
                const Integer j = decode_morton_3DY(gl_index);
                const Integer k = decode_morton_3DZ(gl_index);

                ref_octant = Octant(i, j, k);
                break;
            }
        }
        return ref_octant;
    }

    template <Integer Type>
    MARS_INLINE_FUNCTION Integer get_sfc_from_octant(const Octant &o) {
        Integer enc_oc = -1;

        switch (Type) {
            case ElementType::Quad4: {
                enc_oc = encode_morton_2D(o.x, o.y);
                break;
            }
            case ElementType::Hex8: {
                enc_oc = encode_morton_3D(o.x, o.y, o.z);
                break;
            }
        }

        return enc_oc;
    }

    //-1 for all boundary. 0 left, 1 right, 2 down, 3 up and 4 and 5 for z dim.
    template <Integer Type>
    MARS_INLINE_FUNCTION bool is_boundary_sfc(const Integer sfc,
                                              const int xdim,
                                              const int ydim,
                                              const int zdim,
                                              const Integer Face = -1) {
        Octant o = get_octant_from_sfc<Type>(sfc);
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
            {"all", -1}, {"left", 0}, {"right", 1}, {"front", 2}, {" back ", 3}, {" bottom ", 4}, {" top ", 5}};
        return find_map_side(side_map, side);
    }

    template <Integer Type>
    MARS_INLINE_FUNCTION void get_vertex_coordinates_from_sfc(const Integer sfc,
                                                              double *point,
                                                              const Integer xDim,
                                                              const Integer yDim,
                                                              const Integer zDim) {
        Octant o = get_octant_from_sfc<Type>(sfc);
        o.get_vertex_coordinates<Type>(point, xDim, yDim, zDim);
    }

}  // namespace mars

#endif  // MARS_DIST_OCTANT_HPP
