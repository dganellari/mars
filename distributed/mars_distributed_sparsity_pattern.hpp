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
    template <Integer ExLabel = DofLabel::lNone, typename... ST>
    class SparsityPattern {
    public:
        MARS_INLINE_FUNCTION
        SparsityPattern(ST... f) : stencils(std::make_tuple(f...)) { generate_pattern(); }

        MARS_INLINE_FUNCTION
        const graph_type get_sparsity_pattern() const { return sparsity_pattern; }

        using stencil_tuple = std::tuple<ST...>;

        template <typename H>
        static MARS_INLINE_FUNCTION bool conditional(const Integer local_dof, const H &handler) {
            return (local_dof > -1 && handler.get_label(local_dof) != ExLabel &&
                    !handler.template is_boundary<H::ElemType>(local_dof));
        }

        struct CountUniqueDofs {
            CountUniqueDofs(ViewVectorType<Integer> c) : counter(c) {}

            template <typename S>
            MARS_INLINE_FUNCTION void count_dof(const S &st, const Integer stencil_index, Integer &count) const {
                for (int i = 0; i < st.get_length(); i++) {
                    const Integer local_dof = st.get_value(stencil_index, i);
                    if (conditional(local_dof, st.get_dof_handler())) {
                        count++;
                    }
                }

                for (int i = 0; i < st.get_face_length(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = st.get_face_value(stencil_index, i);
                    if (conditional(local_dof, st.get_dof_handler())) {
                        count++;
                    }
                }

                for (int i = 0; i < st.get_corner_length(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = st.get_corner_value(stencil_index, i);
                    if (conditional(local_dof, st.get_dof_handler())) {
                        count++;
                    }
                }
            }

            template <typename S>
            MARS_INLINE_FUNCTION void count_unique(S st) const {
                st.iterate(MARS_LAMBDA(const Integer stencil_index) {
                    Integer count = 0;

                    count_dof(st, stencil_index, count);

                    const Integer dof = st.get_value(stencil_index, 0);
                    const Integer total =
                        st.get_dof_handler().template is_boundary<S::DHandler::ElemType>(dof) ? 0 : count;
                    counter(st.get_dof_handler().local_to_global(dof)) = total;
                });
            }

            template <typename S>
            MARS_INLINE_FUNCTION void operator()(S &stencil, size_t I) const {
                count_unique<S>(stencil);
            }

            ViewVectorType<Integer> counter;
        };

        struct InsertSortedDofs {
            InsertSortedDofs(ViewVectorType<Integer> r, ViewVectorType<Integer> c) : row_ptr(r), col_idx(c) {}

            MARS_INLINE_FUNCTION
            void insert_sorted(const ViewVectorType<Integer> &col,
                               const Integer index,
                               const Integer value,
                               Integer &count) const {
                Integer i = 0;
                while (value > col(index + i) && i < count) {
                    i++;
                }

                for (int j = count; j > i; --j) {
                    col(index + j) = col(index + j - 1);
                }
                col(index + i) = value;
                ++count;
            }

            template <typename S>
            MARS_INLINE_FUNCTION void insert_dof(const S &st, const Integer stencil_index, Integer dof) const {
                Integer count = 0;

                Integer index = row_ptr(st.get_dof_handler().local_to_global(dof));

                for (int i = 0; i < st.get_length(); i++) {
                    const Integer local_dof = st.get_value(stencil_index, i);
                    if (conditional(local_dof, st.get_dof_handler())) {
                        const Integer global = st.get_dof_handler().local_to_global(local_dof);
                        insert_sorted(col_idx, index, global, count);
                    }
                }

                for (int i = 0; i < st.get_face_length(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = st.get_face_value(stencil_index, i);
                    if (conditional(local_dof, st.get_dof_handler())) {
                        const Integer global = st.get_dof_handler().local_to_global(local_dof);
                        insert_sorted(col_idx, index, global, count);
                    }
                }

                for (int i = 0; i < st.get_corner_length(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = st.get_corner_value(stencil_index, i);
                    if (conditional(local_dof, st.get_dof_handler())) {
                        const Integer global = st.get_dof_handler().local_to_global(local_dof);
                        insert_sorted(col_idx, index, global, count);
                    }
                }
            }

            template <typename S>
            MARS_INLINE_FUNCTION void insert_sorted(S st) const {
                st.iterate(MARS_LAMBDA(const Integer stencil_index) {
                    const Integer dof = st.get_value(stencil_index, 0);
                    if (!st.get_dof_handler().template is_boundary<S::DHandler::ElemType>(dof)) {
                        insert_dof(st, stencil_index, dof);
                    }
                });
            }

            template <typename S>
            MARS_INLINE_FUNCTION void operator()(S &stencil, size_t I) const {
                insert_sorted<S>(stencil);
            }

            ViewVectorType<Integer> row_ptr;
            ViewVectorType<Integer> col_idx;
        };

        template <Integer... dataidx>
        void generate_pattern() {
            auto handler = std::get<0>(stencils).get_dof_handler();
            auto global_size = handler.get_global_dof_size();

            ViewVectorType<Integer> counter("counter", global_size);
            expand_tuple<CountUniqueDofs, stencil_tuple, dataidx...>(CountUniqueDofs(counter), stencils);
            ViewVectorType<Integer> row_ptr("counter_scan", global_size + 1);
            incl_excl_scan(0, global_size, counter, row_ptr);

            auto ss = subview(row_ptr, global_size);
            auto h_ss = create_mirror_view(ss);
            deep_copy(h_ss, ss);

            ViewVectorType<Integer> col_idx("ColdIdx", h_ss());

            expand_tuple<InsertSortedDofs, stencil_tuple, dataidx...>(InsertSortedDofs(row_ptr, col_idx), stencils);

            /* segmented_sort(col_idx);
            sparsity_pattern = graph_type(row_ptr, col_idx); */

            Kokkos::parallel_for(
                "sp:", global_size, MARS_LAMBDA(const Integer i) {
                    for (int j = row_ptr(i); j < row_ptr(i + 1); j++) printf("spi: %li, global: %li\n", i, col_idx(j));
                });
            /*
                        Kokkos::parallel_for(
                            "test", global_size, MARS_LAMBDA(const Integer i) {
                                printf("i: %li, count: %li scan: %li\n", i, counter(i), row_ptr(i));
                            }); */
        }

    private:
        graph_type sparsity_pattern;
        stencil_tuple stencils;
    };
}  // namespace mars

#endif
#endif

#endif
