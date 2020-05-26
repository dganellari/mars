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
    class NonSimplex<Type, DistributedImplementation> final : public ParallelIElem {
    public:

        static constexpr Integer ElemType = Type;

    	SubView<Integer, Type> nodes;
        //SubView<Integer, Type> user_data;

        Integer id = INVALID_INDEX;
       // Integer parent_id = INVALID_INDEX;
        //Integer block = INVALID_INDEX;

        MARS_INLINE_FUNCTION NonSimplex() {}

        MARS_INLINE_FUNCTION NonSimplex(SubView<Integer, Type> n) :
        	nodes(n)
        {
    	}

        MARS_INLINE_FUNCTION Integer n_nodes() const override 
        { 
            return Type; 
        }

        MARS_INLINE_FUNCTION Integer node(const Integer idx) const override 
        {
             assert(idx < Type); 
             return nodes[idx]; 
        }

        MARS_INLINE_FUNCTION Integer type() const override {
            return Type;
        }        
    };

    
    template<Integer Type>
    MARS_INLINE_FUNCTION constexpr static Integer n_nodes(const NonSimplex<Type, DistributedImplementation> &)
    {
        return Type;
    }

	using Quad4DElem = NonSimplex<ElementType::Quad4, DistributedImplementation>;
	using Hex8DElem  = NonSimplex<ElementType::Hex8, DistributedImplementation>;
}

#endif //MARS_NONSIMPLEX_HPP
