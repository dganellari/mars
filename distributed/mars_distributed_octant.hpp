#include "mars_base.hpp"
#include "mars_globals.hpp"
#include "mars_sfc_code.hpp"

namespace mars
{
struct Octant
{
    Integer x;
    Integer y;
    Integer z;

    bool valid;
    bool share_boundary;

    Octant() = default;

    MARS_INLINE_FUNCTION
    Octant(Integer i, Integer j) : x(i), y(j), z(0)
    {
        valid = true;
        share_boundary = false;
    }

    MARS_INLINE_FUNCTION
    Octant(Integer i, Integer j, Integer k) : x(i), y(j), z(k)
    {
        valid = true;
        share_boundary = false;
    }

    MARS_INLINE_FUNCTION
    bool shares_boundary()
    {
        return share_boundary;
    }

    MARS_INLINE_FUNCTION
    void set_shares_boundary()
    {
        share_boundary= true;
    }

    MARS_INLINE_FUNCTION
    bool is_valid() const
    {
        return valid;
    }

    MARS_INLINE_FUNCTION
    void set_invalid()
    {
        valid = false;
    }

    /* template <Integer Type>
    MARS_INLINE_FUNCTION bool shares_boundary_side(const int xDim, const int yDim, const int zDim)
    {
        switch (Type)
        {
        case ElementType::Quad4:
        {
            return (x == -1 || y == -1 || x == xDim || y == yDim);
        }
        case ElementType::Hex8:
        {
            return (x == -1 || y == -1 || z == -1 || x == xDim || y == yDim || z == zDim);
        }
        default:
        {
            printf("The element type is not valid\n");
            return false;
        }
        }
    }
 */
    template <Integer Type>
    MARS_INLINE_FUNCTION
        Integer
        get_global_index(const int xDim, const int yDim)
    {
        switch (Type)
        {
        case ElementType::Quad4:
        {
            return x + (xDim + 1) * y;
        }
        case ElementType::Hex8:
        {
            return elem_index(x, y, z, xDim, yDim);
        }
        default:
        {
            printf("The element type is not valid\n");
            return -1;
        }
        }
    }

    template <Integer Type>
    MARS_INLINE_FUNCTION void get_vertex_coordinates(double *point, const Integer xDim, const Integer yDim, const Integer zDim)
    {
        switch (Type)
        {
        case ElementType::Quad4:
        {
            point[0] = static_cast<Real>(x) / static_cast<Real>(xDim);
            point[1] = static_cast<Real>(y) / static_cast<Real>(yDim);
            point[2] = 0.;
            break;
        }
        case ElementType::Hex8:
        {
            point[0] = static_cast<Real>(x) / static_cast<Real>(xDim);
            point[1] = static_cast<Real>(y) / static_cast<Real>(yDim);
            point[2] = static_cast<Real>(z) / static_cast<Real>(zDim);
            break;
        }
        }
    }

    template <Integer Type>
    MARS_INLINE_FUNCTION Octant sfc_face_nbh(const int face)
    {
        switch (Type)
        {
        case ElementType::Quad4:
        {
            //adapted from the p4est corner neighbor for the mesh generation
            const Integer x_ = x + ((face == 0) ? -1 : (face == 1) ? 1 : 0);
            const Integer y_ = y + ((face == 2) ? -1 : (face == 3) ? 1 : 0);

            return Octant(x_, y_);
        }
        case ElementType::Hex8:
        {

            //adapted from the p4est corner neighbor for the mesh generation
            const Integer x_ = x + ((face == 0) ? -1 : (face == 1) ? 1 : 0);
            const Integer y_ = y + ((face == 2) ? -1 : (face == 3) ? 1 : 0);
            const Integer z_ = z + ((face == 4) ? -1 : (face == 5) ? 1 : 0);

            return Octant(x_, y_, z_);
        }
        default:
        {
            printf("The element type is not valid\n");
            return Octant();
        }
        }
    }

    template <Integer Type>
    MARS_INLINE_FUNCTION void validate_nbh(const Integer xDim, const Integer yDim, const Integer zDim, const bool periodic)
    {
        switch (Type)
        {
        case ElementType::Quad4:
        {
            if (x < 0 || y < 0 || x >= xDim || y >= yDim)
            {
                //check only xdim and ydim and not for -1 because we do not want to process the same face twice
                if (periodic)
                {
                    x = (xDim + x) % xDim;
                    y = (yDim + y) % yDim;
                }
                else
                {
                    set_invalid();
                }
                set_shares_boundary();
            }
            break;
        }
        case ElementType::Hex8:
        {

            if (x < 0 || y < 0 || z < 0 || x >= xDim || y >= yDim || z >= zDim)
            {
                if (periodic)
                {
                    x = (xDim + x) % xDim;
                    y = (yDim + y) % yDim;
                    z = (zDim + z) % zDim;
                }
                else
                {
                    set_invalid();
                }
                set_shares_boundary();
            }
            break;
        }
        default:
        {
            printf("Could not validate: The element type is not valid\n");
            break;
        }
        }
    }


    template <Integer Type>
    MARS_INLINE_FUNCTION Octant face_nbh(const int face, const Integer xDim, const Integer yDim, const Integer zDim, const bool periodic)
    {
        Octant nbh = sfc_face_nbh<Type>(face);
        nbh.validate_nbh<Type>(xDim, yDim, zDim, periodic);
        return nbh;
    }


    template <Integer Type>
    MARS_INLINE_FUNCTION Octant sfc_corner_nbh(const int corner)
    {
        switch (Type)
        {
        case ElementType::Quad4:
        {
            //adapted from the p4est corner neighbor for the mesh generation
            const Integer x_ = x + 2 * (corner & 1) - 1;
            const Integer y_ = y + (corner & 2) - 1;

            return Octant(x_, y_);
        }
        case ElementType::Hex8:
        {
            //adapted from the p4est corner neighbor for the mesh generation
            const Integer x_ = x + 2 * (corner & 1) - 1;
            const Integer y_ = y + (corner & 2) - 1;
            const Integer z_ = z + (corner & 4) / 2 - 1;

            return Octant(x_, y_, z_);
        }
        default:
        {
            printf("The element type is not valid\n");
            return Octant();
        }
        }
    }

    template <Integer Type>
    MARS_INLINE_FUNCTION Octant corner_nbh(const int corner, const Integer xDim, const Integer yDim, const Integer zDim, const bool periodic)
    {
        Octant nbh = sfc_corner_nbh<Type>(corner);
        nbh.validate_nbh<Type>(xDim, yDim, zDim, periodic);
        return nbh;
    }
};

template <Integer Type>
MARS_INLINE_FUNCTION Octant get_octant_from_sfc(const Integer gl_index)
{
    Octant ref_octant;

    switch (Type)
    {
    case ElementType::Quad4:
    {
        const Integer i = decode_morton_2DX(gl_index);
        const Integer j = decode_morton_2DY(gl_index);

        ref_octant = Octant(i, j);
        break;
    }
    case ElementType::Hex8:
    {
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
MARS_INLINE_FUNCTION Integer get_sfc_from_octant(const Octant &o)
{
    Integer enc_oc = -1;

    switch (Type)
    {
    case ElementType::Quad4:
    {
        enc_oc = encode_morton_2D(o.x, o.y);
        break;
    }
    case ElementType::Hex8:
    {
        enc_oc = encode_morton_3D(o.x, o.y, o.z);
        break;
    }
    }

    return enc_oc;
}

template <Integer Type>
MARS_INLINE_FUNCTION void get_vertex_coordinates_from_sfc(const Integer gl_index, double *point, const Integer xDim, const Integer yDim, const Integer zDim)
{
    switch (Type)
    {
    case ElementType::Quad4:
    {
        point[0] = static_cast<Real>(decode_morton_2DX(gl_index)) / static_cast<Real>(xDim);
        point[1] = static_cast<Real>(decode_morton_2DY(gl_index)) / static_cast<Real>(yDim);
        point[2] = 0.;
        break;
    }
    case ElementType::Hex8:
    {
        point[0] = static_cast<Real>(decode_morton_3DX(gl_index)) / static_cast<Real>(xDim);
        point[1] = static_cast<Real>(decode_morton_3DY(gl_index)) / static_cast<Real>(yDim);
        point[2] = static_cast<Real>(decode_morton_3DZ(gl_index)) / static_cast<Real>(zDim);
        break;
    }
    }
}

} // namespace mars
