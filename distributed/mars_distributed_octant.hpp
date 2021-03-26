#ifndef MARS_DIST_OCTANT_HPP
#define MARS_DIST_OCTANT_HPP

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

        template <Integer Type, Integer Face = -1>
        MARS_INLINE_FUNCTION bool is_boundary(const int xdim, const int ydim, const int zdim) {
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

        MARS_INLINE_FUNCTION Octant sfc_edge_nbh(const int edge) const {
            // adapted from the p4est corner neighbor for the mesh generation
            assert(0 <= edge && edge < 12);
            auto direction = edge / 4;
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

        template <Integer Type, Integer ManifoldDim>
        MARS_INLINE_FUNCTION void one_ring_nbh(Octant one_ring[Type],
                                               const Integer xDim,
                                               const Integer yDim,
                                               const Integer zDim,
                                               const bool periodic) const {
            for (int i = 0; i < ManifoldDim; i++) {
                for (int j = 0; j < ManifoldDim; j++)
                /* for (int z = 0; z < ManifoldDim; z++) // 3D case*/
                {
                    Octant nbh(x - i, y - j);
                    nbh.validate_nbh<Type>(xDim, yDim, zDim, periodic);
                    one_ring[i * ManifoldDim + j] = nbh;
                }
            }
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
    template <Integer Type, Integer Face = -1>
    MARS_INLINE_FUNCTION bool is_boundary_sfc(const Integer sfc, const int xdim, const int ydim, const int zdim) {
        Octant o = get_octant_from_sfc<Type>(sfc);
        return o.is_boundary<Type, Face>(xdim, ydim, zdim);
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
