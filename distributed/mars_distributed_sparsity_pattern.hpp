#ifndef GENERATION_MARS_DISTRIBUTED_SP_HPP_
#define GENERATION_MARS_DISTRIBUTED_SP_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_finite_element.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    template <class Finite>
    class SparsityPattern {
    public:

        MARS_INLINE_FUNCTION
        SparsityPattern(Finite f) : finite(f) {}

        MARS_INLINE_FUNCTION
        const graph_type get_sparsity_pattern() const { return sparsity_pattern; }


    private:
        graph_type sparsity_pattern;
        Finite finite_item;
    };
}  // namespace mars

#endif
#endif

#endif
