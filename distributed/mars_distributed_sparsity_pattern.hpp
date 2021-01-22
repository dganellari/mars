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

    template <typename... Stencil>
    class SparsityPattern {
    public:
        MARS_INLINE_FUNCTION
        SparsityPattern(Stencil... f) : stencils(std::make_tuple(f...)) {}

        MARS_INLINE_FUNCTION
        const graph_type get_sparsity_pattern() const { return sparsity_pattern; }

        using stencil_tuple = std::tuple<Stencil...>;

        class CountUnique {
            CountUnique(ViewVectorType<Integer> c) : counter(c) {}

            MARS_INLINE_FUNCTION
            template <typename H>
            conditional_count(const Integer local_dof, H &handler, Integer &count) {
                if (local_dof > -1 && handler.get_label(local_dof) != DofLabel::lCorner &&
                    !handler.is_boundary(local_dof)) {
                    count++;
                }
            }

            MARS_INLINE_FUNCTION
            template <typename S>
            void count_unique(S st) {
                st.iterate(MARS_LAMBDA(const Integer stencil_index) {
                    Integer count = 0;
                    for (int i = 0; i < st.get_stencil().get_length(); i++) {
                        const Integer local_dof = st.get_value(stencil_index, i);
                        conditional_count(local_dof, st.get_dof_handler(), count);
                    }

                    for (int i = 0; i < st.get_face_length(); i++) {
                        // get the local dof of the i-th index within thelement
                        const Integer local_dof = st.get_face_value(stencil_index, i);
                        conditional_count(local_dof, st.get_dof_handler(), count);
                    }

                    for (int i = 0; i < st.get_corner_length(); i++) {
                        // get the local dof of the i-th index within thelement
                        const Integer local_dof = st.get_corner_value(stencil_index, i);
                        conditional_count(local_dof, st.get_dof_handler(), count);
                    }

                    const Integer dof = st.get_value(stencil_index, 0);
                    //TODO:check if the correct way to do it.
                    counter(st.get_dof_handler().sfc_to_global(dof)) = count;
                });
            }

            template <typename S>
            MARS_INLINE_FUNCTION void operator()(S stencil) const {
                count_unique(stencil);
            }

            ViewVectorType<Integer> counter;
        };

        void exclude_dof(const Integer label) {

        }

        template <Integer... dataidx>
        void generate_pattern() {
            //TODO: make a global owned dof size for all handlers.
            ViewVectorType<Integer> counter("counter", st.get_dof_handler().get_owned_dof_size());
            expand_tuple<CountUnique, stencil_tuple, dataidx...>(CountUnique(counter), stencils);
            scan(counter);
        }

    private:
        graph_type sparsity_pattern;
        stencil_tuple stencils;
    };
}  // namespace mars

#endif
#endif

#endif
