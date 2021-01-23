#ifndef GENERATION_MARS_DISTRIBUTED_SP_HPP_
#define GENERATION_MARS_DISTRIBUTED_SP_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include "mars_distributed_finite_element.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    /* template <class Finite>
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
 */
    template <Integer ExLabel = DofLabel::lNone, typename... Stencil>
    class SparsityPattern {
    public:
        MARS_INLINE_FUNCTION
        SparsityPattern(Stencil... f) : stencils(std::make_tuple(f...)) { generate_pattern(); }

        MARS_INLINE_FUNCTION
        const graph_type get_sparsity_pattern() const { return sparsity_pattern; }

        using stencil_tuple = std::tuple<Stencil...>;

        struct CountUnique {
            CountUnique(ViewVectorType<Integer> c) : counter(c) {}

            template <typename H>
            MARS_INLINE_FUNCTION void conditional_count(const Integer local_dof,const H &handler, Integer &count) const {
                if (local_dof > -1 && handler.get_label(local_dof) != ExLabel &&
                    !handler.template is_boundary<H::ElemType>(local_dof)) {
                    count++;
                }
            }

            template <typename S>
            MARS_INLINE_FUNCTION void count_dof(const S &st, const Integer stencil_index, Integer &count) const {
                for (int i = 0; i < st.get_length(); i++) {
                    const Integer local_dof = st.get_value(stencil_index, i);
                    conditional_count(local_dof, st.get_dof_handler(), count);
                }
            }

            template <typename DofHandler>
            MARS_INLINE_FUNCTION void count_dof(const StokesStencil<DOFHandler> &st,
                                                const Integer stencil_index,
                                                Integer &count) const {
                for (int i = 0; i < st.get_length(); i++) {
                    const Integer local_dof = st.get_value(stencil_index, i);
                    conditional_count(local_dof, st.get_dof_handler(), count);
                }

                for (int i = 0; i < st.get_face_length(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = st.get_face_value(stencil_index, i);
                    conditional_count(local_dof, st.get_dof_handler(), count);
                }
            }

            template <typename DofHandler>
            MARS_INLINE_FUNCTION void count_dof(const FullStokesStencil<DOFHandler> &st,
                                                const Integer stencil_index,
                                                Integer &count) const {
                for (int i = 0; i < st.get_length(); i++) {
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
            }

            template <typename S>
            MARS_INLINE_FUNCTION void count_unique(S st) const {
                st.iterate(MARS_LAMBDA(const Integer stencil_index) {
                    Integer count = 0;

                    count_dof(st, stencil_index, count);

                    const Integer dof = st.get_value(stencil_index, 0);
                    counter(st.get_dof_handler().local_to_global(dof)) = count;
                });
            }

            template <typename S>
            MARS_INLINE_FUNCTION void operator()(S &stencil, size_t I) const {
                count_unique<S>(stencil);
            }

            ViewVectorType<Integer> counter;
        };

        template <Integer... dataidx>
        void generate_pattern() {
            auto handler = std::get<0>(stencils).get_dof_handler();
            ViewVectorType<Integer> counter("counter", handler.get_global_dof_size());
            expand_tuple<CountUnique, stencil_tuple, dataidx...>(CountUnique(counter), stencils);
            /* scan(counter); */

            Kokkos::parallel_for(
                "test", handler.get_global_dof_size(), MARS_LAMBDA(const Integer i) {
                    printf("i: %li, count: %li \n", i, counter(i));
                });
        }

    private:
        graph_type sparsity_pattern;
        stencil_tuple stencils;
    };
}  // namespace mars

#endif
#endif

#endif
