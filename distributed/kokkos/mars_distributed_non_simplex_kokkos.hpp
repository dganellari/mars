#ifndef MARS_DIST_NON_SIMPLEX_KOKKOS_HPP
#define MARS_DIST_NON_SIMPLEX_KOKKOS_HPP

#include "mars_base.hpp"
#include "mars_imesh_kokkos.hpp"

#include <array>
#include <vector>

#include <cmath>
#include <cassert>
#include <initializer_list>
#include <algorithm>

#include "mars_sub_view.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

    template<Integer Type>
    class NonSimplex<Type, DistributedImplementation> : public ParallelIElem {
    public:

        static constexpr Integer ElemType = Type;

        /* SubView<Integer, Type> nodes; */
        //SubView<Integer, Type> user_data;

        Integer elem_id = INVALID_INDEX;
       // Integer parent_id = INVALID_INDEX;
        //Integer block = INVALID_INDEX;

        bool is_ghost;

        MARS_INLINE_FUNCTION NonSimplex() {}

        MARS_INLINE_FUNCTION NonSimplex(bool ghost) : is_ghost(ghost)
        {
        }

        MARS_INLINE_FUNCTION
        void set_elem_id(const Integer idx)
        {
            elem_id = idx;
        }

        MARS_INLINE_FUNCTION
        void set_ghost()
        {
            is_ghost = true;
        }


        /* MARS_INLINE_FUNCTION NonSimplex(SubView<Integer, Type> n) :
        	nodes(n)
        {
        } */

        MARS_INLINE_FUNCTION Integer n_nodes() const override
        {
            return Type;
        }

        /* MARS_INLINE_FUNCTION Integer node(const Integer idx) const override
        {
             assert(idx < Type);
             return nodes[idx];
        } */

        MARS_INLINE_FUNCTION Integer type() const override {
            return Type;
        }
    };

    template<Integer Type>
    class Side<Type, DistributedImplementation> : public NonSimplex<Type, DistributedImplementation>
    {
        Integer face;

        MARS_INLINE_FUNCTION
        Integer get_face()
        {
            return face;
        }

        MARS_INLINE_FUNCTION
        void set_face(const Integer f)
        {
            face =f;
        }
    };

    template <Integer Type>
    using FaceSide = Side<Type / 2, DistributedImplementation>;

    template <Integer Type>
    class Face
    {
    public:
        MARS_INLINE_FUNCTION
        Face(const int direction)
        {
            sides[0].set_face(2 * direction + 1);
            sides[1].set_face(2 * direction);
            set_direction(direction);
        }

        MARS_INLINE_FUNCTION
        Integer get_direction()
        {
            return direction;
        }

        MARS_INLINE_FUNCTION
        void set_direction(const Integer d)
        {
            direction = d;
        }

        MARS_INLINE_FUNCTION
        FaceSide<Type> * get_sides()
        {
            return sides;
        }

        MARS_INLINE_FUNCTION
        Integer get_second_side_face()
        {
            return sides[1].get_face();
        }

        MARS_INLINE_FUNCTION
        Integer get_first_side_face()
        {
            return sides[0].get_face();
        }

    private:

        Integer direction;
        FaceSide<Type> sides[2];
    };

    template<Integer Type>
    MARS_INLINE_FUNCTION constexpr static Integer n_nodes(const NonSimplex<Type, DistributedImplementation> &)
    {
        return Type;
    }

	using Quad4DElem = NonSimplex<ElementType::Quad4, DistributedImplementation>;
	using Hex8DElem  = NonSimplex<ElementType::Hex8, DistributedImplementation>;

	using Quad4DSide= Side<ElementType::Quad4, DistributedImplementation>;
	using Hex8DSide= Side<ElementType::Hex8, DistributedImplementation>;

}

#endif //MARS_NONSIMPLEX_HPP
