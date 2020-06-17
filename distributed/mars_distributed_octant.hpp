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

    MARS_INLINE_FUNCTION
    Octant() = default;

    MARS_INLINE_FUNCTION
    Octant(Integer i, Integer j) : x(i), y(j), z(0)
    {
        valid = true;
    }

    MARS_INLINE_FUNCTION
    Octant(Integer i, Integer j, Integer k) : x(i), y(j), z(k)
    {
        valid = true;
    }

    MARS_INLINE_FUNCTION
    bool is_valid()
    {
        return valid;
    }

    MARS_INLINE_FUNCTION
    void set_invalid()
    {
        valid = false;
    }

    template <Integer Type>
    MARS_INLINE_FUNCTION
    bool shares_boundary_side(const int xDim, const int yDim, const int zDim)
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
        }
    }

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
MARS_INLINE_FUNCTION void get_vertex_coordinates_from_sfc(const Integer gl_index, double* point, const Integer xDim, const Integer yDim, const Integer zDim)
{
    switch (Type)
    {
    case ElementType::Quad4:
    {
        point[0] = static_cast<Real>(decode_morton_2DX(gl_index)) / static_cast<Real>(xDim);
        point[1] = static_cast<Real>(decode_morton_2DY(gl_index)) / static_cast<Real>(yDim);
        point[2]=0.;
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

template <Integer Type>
MARS_INLINE_FUNCTION Octant face_nbh(const Octant &ref_octant, const int face,
                                            const Integer xDim, const Integer yDim, const Integer zDim)
{

    switch (Type)
    {
    case ElementType::Quad4:
    {
        //adapted from the p4est corner neighbor for the mesh generation
        const Integer x = ref_octant.x + ((face == 0) ? -1 : (face == 1) ? 1 : 0);
        const Integer y = ref_octant.y + ((face == 2) ? -1 : (face == 3) ? 1 : 0);

        Octant o(x, y);
        if (o.x < 0 || o.y < 0 || o.x >= xDim || o.y >= yDim)
            o.set_invalid();

        return o;
    }
    case ElementType::Hex8:
    {

        //adapted from the p4est corner neighbor for the mesh generation
        const Integer x = ref_octant.x + ((face == 0) ? -1 : (face == 1) ? 1 : 0);
        const Integer y = ref_octant.y + ((face == 2) ? -1 : (face == 3) ? 1 : 0);
        const Integer z = ref_octant.z + ((face == 4) ? -1 : (face == 5) ? 1 : 0);

        Octant o(x, y, z);
        if (o.x < 0 || o.y < 0 || o.z < 0 || o.x >= xDim || o.y >= yDim || o.z >= zDim)
            o.set_invalid();
        return o;
    }
    }
}

template <Integer Type>
MARS_INLINE_FUNCTION Octant corner_nbh(const Octant &ref_octant, const int corner,
                                              const Integer xDim, const Integer yDim, const Integer zDim)
{

    switch (Type)
    {
    case ElementType::Quad4:
    {
        //adapted from the p4est corner neighbor for the mesh generation
        const Integer x = ref_octant.x + 2 * (corner & 1) - 1;
        const Integer y = ref_octant.y + (corner & 2) - 1;

        Octant o(x, y);
        if (o.x < 0 || o.y < 0 || o.x >= xDim || o.y >= yDim)
            o.set_invalid();

        return o;
    }
    case ElementType::Hex8:
    {
        //adapted from the p4est corner neighbor for the mesh generation
        const Integer x = ref_octant.x + 2 * (corner & 1) - 1;
        ;
        const Integer y = ref_octant.y + (corner & 2) - 1;
        const Integer z = ref_octant.z + (corner & 4) / 2 - 1;

        Octant o(x, y, z);
        if (o.x < 0 || o.y < 0 || o.z < 0 || o.x >= xDim || o.y >= yDim || o.z >= zDim)
            o.set_invalid();

        return o;
    }
    }
}

} // namespace mars
