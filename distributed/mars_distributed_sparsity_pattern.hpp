#ifndef GENERATION_MARS_DISTRIBUTED_SP_HPP_
#define GENERATION_MARS_DISTRIBUTED_SP_HPP_

#ifdef WITH_MPI
#ifdef WITH_KOKKOS
#include <sstream>
#include "mars_distributed_finite_element.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    /* template <class Finite>
    class SparsityPattern {
    public:
        MARS_INLINE_FUNCTION
        SparsityPattern(Finite f) : finite(f) {}

        MARS_INLINE_FUNCTION
        const crs_graph get_sparsity_pattern() const { return sparsity_pattern; }

    private:
        crs_graph sparsity_pattern;
        Finite finite_item;
    };
 */
    template <typename V, typename SHandler, typename... ST>
    class SparsityPattern {
    public:
        using stencil_tuple = std::tuple<ST...>;

        using Scalar = V;
        using Ordinal = default_lno_t;
        using Offset = default_size_type;
        using Layout = default_layout;

        using device_type = typename Kokkos::Device<Kokkos::DefaultExecutionSpace, KokkosSpace>;

        using crs_matrix = typename KokkosSparse::CrsMatrix<Scalar, Ordinal, device_type, void, Offset>;

        using crs_graph = typename crs_matrix::staticcrsgraph_type;

        using crs_row = typename crs_graph::row_map_type::non_const_type;
        using crs_col = typename crs_graph::entries_type::non_const_type;
        using crs_value = typename crs_matrix::values_type::non_const_type;

        MARS_INLINE_FUNCTION
        SparsityPattern(SHandler h) : dof_handler(h) {}

        /* struct CompareLabel {
            MARS_INLINE_FUNCTION
            CompareLabel(ViewObject<bool> r, Integer l) : result(r), label(l) {}

            template <typename S>
            MARS_INLINE_FUNCTION bool operator()(S &stencil, size_t I) const {
                if (stencil.get_label() == label) result(0) = (result(0) || true);
            }

            ViewObject<bool> result;
            Integer label;
        };

        template <Integer... dataidx>
        static MARS_INLINE_FUNCTION bool expand_conditional(stencil_tuple sts, const Integer label) {
            ViewObject<bool> result("result");
            expand_tuple<CompareLabel, stencil_tuple, dataidx...>(CompareLabel(result, label), sts);
            return result(0);
        } */

        template <typename H>
        static MARS_INLINE_FUNCTION bool conditional(const Integer local_dof, const H &handler) {
            /* return (local_dof > -1 && expand_conditional(sts, handler.get_label(local_dof))); */
            return (local_dof > -1 && (SHandler::dofLabel & handler.get_label(local_dof)));
        }

        struct CountUniqueDofs {
            CountUniqueDofs(ViewVectorType<Integer> c, SHandler h) : counter(c), dhandler(h) {}

            template <typename S>
            MARS_INLINE_FUNCTION void count_dof(const S &st, const Integer stencil_index, Integer &count) const {
                for (int i = 0; i < st.get_length(); i++) {
                    const Integer local_dof = st.get_value(stencil_index, i);
                    if (conditional(local_dof, dhandler)) {
                        count++;
                    }
                }

                for (int i = 0; i < st.get_face_length(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = st.get_face_value(stencil_index, i);
                    if (conditional(local_dof, dhandler)) {
                        count++;
                    }
                }

                for (int i = 0; i < st.get_corner_length(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = st.get_corner_value(stencil_index, i);
                    if (conditional(local_dof, dhandler)) {
                        count++;
                    }
                }
            }

            // TODO:Specialize for stokes as in the next method below since normally the boundary is excluded.
            // In this case the total would be 0 for the boundary and in the specialization 1.
            /* template <typename S> */
            MARS_INLINE_FUNCTION Integer get_total(const Integer dof, Integer &count) const {
                return dhandler.template is_boundary<SHandler::ElemType>(dof) ? 1 : count;
                /* return handler.template is_boundary<SHandler::ElemType>(dof) ? 0 : count; */
            }

            /* template <>
            MARS_INLINE_FUNCTION Integer get_total<StokesStencil<DHandler>>(const StokesStencil<DHandler> &st) const {
                return dhandler.template is_boundary<SHandler::ElemType>(dof) ? 1 : count;
            } */

            template <typename S>
            MARS_INLINE_FUNCTION void count_unique(S st) const {
                st.iterate(MARS_LAMBDA(const Integer stencil_index) {
                    Integer count = 0;

                    count_dof(st, stencil_index, count);

                    const Integer dof = st.get_value(stencil_index, 0);
                    const Integer total = get_total(dof, count);
                    counter(dhandler.local_to_global(dof)) = total;
                });
            }

            template <typename S>
            MARS_INLINE_FUNCTION void operator()(S &stencil, size_t I) const {
                count_unique<S>(stencil);
            }

            ViewVectorType<Integer> counter;
            SHandler dhandler;
        };

        struct InsertSortedDofs {
            InsertSortedDofs(crs_row r, crs_col c, SHandler h) : row_ptr(r), col_idx(c), dhandler(h) {}

            MARS_INLINE_FUNCTION
            void insert_sorted(const crs_col &col, const Integer index, const Integer value, Integer &count) const {
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

                Integer index = row_ptr(dhandler.local_to_global(dof));

                for (int i = 0; i < st.get_length(); i++) {
                    const Integer local_dof = st.get_value(stencil_index, i);
                    if (conditional(local_dof, dhandler)) {
                        const Integer global = dhandler.local_to_global(local_dof);
                        insert_sorted(col_idx, index, global, count);
                    }
                }

                for (int i = 0; i < st.get_face_length(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = st.get_face_value(stencil_index, i);
                    if (conditional(local_dof, dhandler)) {
                        const Integer global = dhandler.local_to_global(local_dof);
                        insert_sorted(col_idx, index, global, count);
                    }
                }

                for (int i = 0; i < st.get_corner_length(); i++) {
                    // get the local dof of the i-th index within thelement
                    const Integer local_dof = st.get_corner_value(stencil_index, i);
                    if (conditional(local_dof, dhandler)) {
                        const Integer global = dhandler.local_to_global(local_dof);
                        insert_sorted(col_idx, index, global, count);
                    }
                }
            }

            /* TODO:specialize for stokes because normally the boundary would be excluded. */
            template <typename S>
            MARS_INLINE_FUNCTION void insert_sorted(S st) const {
                st.iterate(MARS_LAMBDA(const Integer stencil_index) {
                    const Integer dof = st.get_value(stencil_index, 0);
                    if (dhandler.template is_boundary<SHandler::ElemType>(dof)) {
                        Integer global = dhandler.local_to_global(dof);
                        Integer index = row_ptr(global);
                        col_idx(index) = global;
                    } else {
                        insert_dof(st, stencil_index, dof);
                    }
                });
            }

            template <typename S>
            MARS_INLINE_FUNCTION void operator()(S &stencil, size_t I) const {
                insert_sorted<S>(stencil);
            }

            crs_row row_ptr;
            crs_col col_idx;
            SHandler dhandler;
        };

        template <Integer... dataidx>
        void build_pattern(ST... s) {
            stencil_tuple stencils(std::make_tuple(s...));

            auto global_size = get_dof_handler().get_global_dof_size();
            printf("global_size %li\n", global_size);

            ViewVectorType<Integer> counter("counter", global_size);
            expand_tuple<CountUniqueDofs, stencil_tuple, dataidx...>(CountUniqueDofs(counter, get_dof_handler()),
                                                                     stencils);
            crs_row row_ptr("counter_scan", global_size + 1);
            incl_excl_scan(0, global_size, counter, row_ptr);

            auto ss = subview(row_ptr, global_size);
            auto h_ss = create_mirror_view(ss);
            deep_copy(h_ss, ss);

            crs_col col_idx("ColdIdx", h_ss());
            crs_value values("values", h_ss());

            expand_tuple<InsertSortedDofs, stencil_tuple, dataidx...>(
                InsertSortedDofs(row_ptr, col_idx, get_dof_handler()), stencils);

            /* TODO: */
            /* segmented_sort(col_idx);*/

            sparsity_pattern = crs_graph(col_idx, row_ptr);

            matrix = crs_matrix("crs_matrix", global_size, values, sparsity_pattern);

            /* Kokkos::parallel_for(
                "sp:", global_size, MARS_LAMBDA(const Integer i) {
                    for (int j = row_ptr(i); j < row_ptr(i + 1); j++) printf("spi: %li, global: %li\n", i, col_idx(j));
                });

            Kokkos::parallel_for(
                "test", global_size, MARS_LAMBDA(const Integer i) {
                    printf("i: %li, count: %li scan: %li\n", i, counter(i), row_ptr(i));
                }); */
        }

        MARS_INLINE_FUNCTION
        const crs_graph get_sparsity_pattern() const { return sparsity_pattern; }

        MARS_INLINE_FUNCTION
        const Integer get_col(const Integer index) const { return get_sparsity_pattern().entries(index); }

        MARS_INLINE_FUNCTION
        const crs_col get_col() const { return get_sparsity_pattern().entries; }

        MARS_INLINE_FUNCTION
        const Integer get_row_map(const Integer index) const { return get_sparsity_pattern().row_map(index); }

        MARS_INLINE_FUNCTION
        const crs_row get_row_map() { return get_sparsity_pattern().row_map; }

        MARS_INLINE_FUNCTION
        const Integer get_col_index_from_global(const Integer row, const Integer col) const {
            const Integer row_idx = sparsity_pattern.row_map(row);
            const Integer next_row_idx = sparsity_pattern.row_map(row + 1);

            const Integer i = binary_search(sparsity_pattern.entries, row_idx, next_row_idx, col);
            return i;
        }

        MARS_INLINE_FUNCTION
        const Integer get_col_index(const Integer local_row, const Integer local_col) const {
            auto row = get_dof_handler().local_to_global(local_row);
            auto col = get_dof_handler().local_to_global(local_col);

            return get_col_index_from_global(row, col);
        }
        MARS_INLINE_FUNCTION
        const void set_value(const Integer index, const V val) const { matrix.values(index) = val; }

        MARS_INLINE_FUNCTION
        const void set_value(const Integer row, const Integer col, const V val) const {
            const Integer index = get_col_index(row, col);
            if (index > -1) matrix.values(index) = val;
        }

        MARS_INLINE_FUNCTION
        const void set_value_from_global(const Integer row, const Integer col, const V val) const {
            const Integer index = get_col_index_from_global(row, col);
            if (index > -1) matrix.values(index) = val;
        }

        MARS_INLINE_FUNCTION
        const V get_value(const Integer index) const { return matrix.values(index); }

        MARS_INLINE_FUNCTION
        const V get_value(const Integer row, const Integer col) const { return matrix.values(get_col_index(row, col)); }

        MARS_INLINE_FUNCTION
        const V get_value_from_global(const Integer row, const Integer col) const {
            return matrix.values(get_col_index_from_global(row, col));
        }

        MARS_INLINE_FUNCTION
        const Integer get_num_rows() const { return matrix.numRows(); }

        /* ***************** print ******************************************** */

        bool write(const std::string &file_path) {
            auto proc = mars::rank(get_dof_handler().get_context());
            auto path = file_path + "_" + std::to_string(proc);

            std::cout << "Writing SparsityPattern to " << path << " file." << std::endl;

            std::ofstream os;
            os.open(path.c_str());
            if (!os.good()) {
                os.close();
                return false;
            }

            auto row_h = Kokkos::create_mirror_view(sparsity_pattern.row_map);
            auto col_h = Kokkos::create_mirror_view(sparsity_pattern.entries);
            auto val_h = Kokkos::create_mirror_view(matrix.values);
            Kokkos::deep_copy(row_h, sparsity_pattern.row_map);
            Kokkos::deep_copy(col_h, sparsity_pattern.entries);
            Kokkos::deep_copy(val_h, matrix.values);

            os << "Number of values: " << matrix.nnz() << " Number of Cols: " << matrix.numCols()
               << " Number of Rows: " << matrix.numRows() << "\n";

            for (int i = 0; i < matrix.nnz(); ++i) {
                os << matrix.values(i) << " ";
            }
            os << "\n";

            for (int i = 0; i < matrix.nnz(); ++i) {
                os << sparsity_pattern.entries(i) << " ";
            }
            os << "\n";

            for (int i = 0; i < matrix.numRows() + 1; ++i) {
                os << sparsity_pattern.row_map(i) << " ";
            }
            os << "\n";

            os.close();
            return true;
        }
        /* ********************************************************************************* */

        MARS_INLINE_FUNCTION
        const SHandler &get_dof_handler() const { return dof_handler; }

    private:
        crs_graph sparsity_pattern;
        crs_matrix matrix;

        SHandler dof_handler;
    };
}  // namespace mars

#endif
#endif

#endif
