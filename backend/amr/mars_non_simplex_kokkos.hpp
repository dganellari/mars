#ifndef MARS_NON_SIMPLEX_KOKKOS_HPP
#define MARS_NON_SIMPLEX_KOKKOS_HPP

#include "mars_base.hpp"
#include "mars_imesh_kokkos.hpp"

#include <array>
#include <vector>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <initializer_list>

#include "mars_sub_view.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

    template <Integer Type>
    class NonSimplex<Type, KokkosImplementation> final : public ParallelIElem {
    public:
        static constexpr Integer ElemType = Type;
        static constexpr Integer NNodes = Type;
        static constexpr Integer Dim = (Type == 4) ? 2 : 3;

        SubView<Integer, Type> nodes;
        // SubView<Integer,2> children; //TODO: templatize for the number of children based onthe select algorithm
        SubView<Integer, Type> side_tags;

        Integer id = INVALID_INDEX;
        Integer parent_id = INVALID_INDEX;
        Integer block = INVALID_INDEX;

        MARS_INLINE_FUNCTION NonSimplex() {}

        MARS_INLINE_FUNCTION NonSimplex(SubView<Integer, Type> n) : nodes(n) {}

        MARS_INLINE_FUNCTION Integer get_block() const  // override
        {
            return block;
        }

        MARS_INLINE_FUNCTION void set_block(const Integer block_id)  // override
        {
            block = block_id;
        }

        MARS_INLINE_FUNCTION Integer n_nodes() const override { return Type; }

        MARS_INLINE_FUNCTION Integer node(const Integer idx) const override {
            assert(idx < Type);
            return nodes[idx];
        }

        MARS_INLINE_FUNCTION Integer type() const override { return Type; }
    };

    template <Integer Type>
    MARS_INLINE_FUNCTION constexpr static Integer n_nodes(const NonSimplex<Type, KokkosImplementation> &) {
        return Type;
    }

    using Quad4PElem = mars::NonSimplex<ElementType::Quad4, KokkosImplementation>;
    using Hex8PElem = mars::NonSimplex<ElementType::Hex8, KokkosImplementation>;
}  // namespace mars

#endif  // MARS_NONSIMPLEX_HPP
