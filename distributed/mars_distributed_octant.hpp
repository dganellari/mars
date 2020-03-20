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
    Octant(Integer i, Integer j) : x(i), y(j)
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
};
} // namespace mars