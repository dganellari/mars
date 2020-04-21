#include "mars_base.hpp"
#include "mars_globals.hpp"
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
    Integer get_global_index(const int xDim, const int yDim)
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
};
} // namespace mars
