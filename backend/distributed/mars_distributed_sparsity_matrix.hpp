#ifndef GENERATION_MARS_DISTRIBUTED_SM_HPP_
#define GENERATION_MARS_DISTRIBUTED_SM_HPP_

#ifdef MARS_ENABLE_MPI
#ifdef MARS_ENABLE_KOKKOS_KERNELS
#include "mars_distributed_sparsity_pattern.hpp"

namespace mars {

    class ISparsityMatrix {
    public:
        virtual MARS_INLINE_FUNCTION ~ISparsityMatrix() {}
    };
    // for default local/global ordinals, you can have for instantiation:
    // LO=default_lno_t
    // GO=default_size_type
    template <typename SP>
    class SparsityMatrix : public ISparsityMatrix {
    public:
        using V = typename SP::Scalar;
        using Ordinal = typename SP::Ordinal;
        using Offset = typename SP::Offset;
        using Layout = default_layout;

        using device_type = typename Kokkos::Device<Kokkos::DefaultExecutionSpace, KokkosSpace>;

        using crs_matrix = typename KokkosSparse::CrsMatrix<V, Ordinal, device_type, void, Offset>;

        using crs_graph = typename crs_matrix::staticcrsgraph_type;

        using crs_row = typename crs_graph::row_map_type::non_const_type;
        using crs_col = typename crs_graph::entries_type::non_const_type;
        using col_index_type = typename crs_graph::entries_type::value_type;
        using crs_value = typename crs_matrix::values_type::non_const_type;

        using DofHandler = typename SP::DofHandler;
        using SparsityPattern = SP;

        MARS_INLINE_FUNCTION
        SparsityMatrix(SparsityPattern sp) : sparsity_pattern(sp) {}

        MARS_INLINE_FUNCTION
        SparsityMatrix(SparsityPattern sp, crs_matrix mat) : sparsity_pattern(sp), matrix(mat) {}

        crs_matrix new_crs_matrix() const { return get_sparsity_pattern().new_crs_matrix(); }

        void build_crs_matrix() { matrix = new_crs_matrix(); }

        MARS_INLINE_FUNCTION
        void matrix_apply_constraints(const Integer row, const V value) const {
            get_sparsity_pattern().matrix_apply_constraints(row, matrix, value);
        }

        MARS_INLINE_FUNCTION
        void set_value(const Integer index, const V val) const { matrix.values(index) = val; }

        MARS_INLINE_FUNCTION
        void set_value(const Integer row, const Integer col, const V val) const {
            const Integer index = get_sparsity_pattern().get_col_index(row, col);
            assert(index > -1);
            if (index > -1) matrix.values(index) = val;
        }

        MARS_INLINE_FUNCTION
        void atomic_add_value(const Integer row, const Integer col, const V val) const {
            const Integer index = get_sparsity_pattern().get_col_index(row, col);
            assert(index > -1);
            if (index > -1) {
                Kokkos::atomic_fetch_add(&matrix.values(index), val);
            }
        }

        template <typename C>
        MARS_INLINE_FUNCTION void atomic_add_value(const Integer row,
                                                   const Integer col,
                                                   const V val,
                                                   C &crs_matrix) const {
            const Integer index = get_sparsity_pattern().get_col_index(row, col);
            assert(index > -1);
            if (index > -1) {
                Kokkos::atomic_fetch_add(&crs_matrix.values(index), val);
            }
        }

        MARS_INLINE_FUNCTION
        void set_value_from_global(const Integer row, const Integer col, const V val) const {
            const Integer index = get_sparsity_pattern().get_col_index_from_global(row, col);
            if (index > -1) matrix.values(index) = val;
        }

        MARS_INLINE_FUNCTION
        const V get_value(const Integer index) const { return matrix.values(index); }

        MARS_INLINE_FUNCTION
        const V get_value(const Integer row, const Integer col) const {
            return matrix.values(get_sparsity_pattern().get_col_index(row, col));
        }

        MARS_INLINE_FUNCTION
        const V get_value_from_global(const Integer row, const Integer col) const {
            return matrix.values(get_sparsity_pattern().get_col_index_from_global(row, col));
        }

        void print_sparsity_pattern() const {
            const Integer size = get_num_rows();
            auto sp = *this;
            Kokkos::parallel_for(
                "for", size, MARS_LAMBDA(const int row) {
                    const Integer start = sp.get_sparsity_pattern().get_row_map(row);
                    const Integer end = sp.get_sparsity_pattern().get_row_map(row + 1);

                    // print only if end - start > 0. Otherwise segfaults.
                    // The row index is not a global index of the current process!
                    for (int i = start; i < end; ++i) {
                        auto value = sp.get_value(i);
                        auto col = sp.get_sparsity_pattern().get_col(i);

                        const Integer local_dof = sp.get_dof_handler().get_owned_dof(row);
                        const Integer global_row = sp.get_dof_handler().local_to_global(local_dof);

                        auto base_col = sp.get_dof_handler().compute_base(col);
                        auto base_row = sp.get_dof_handler().compute_base(global_row);
                        printf("SP - Row_Dof: %li - %li, base_row:%li, col_dof: %li, base_col: %li, value: %lf\n",
                               row,
                               global_row,
                               base_row,
                               col,
                               base_col,
                               value);
                    }
                });
        }

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

            auto row_h = Kokkos::create_mirror_view(get_crs_graph().row_map);
            auto col_h = Kokkos::create_mirror_view(get_crs_graph().entries);
            auto val_h = Kokkos::create_mirror_view(matrix.values);
            Kokkos::deep_copy(row_h, get_crs_graph().row_map);
            Kokkos::deep_copy(col_h, get_crs_graph().entries);
            Kokkos::deep_copy(val_h, matrix.values);

            os << "Number of values: " << matrix.nnz() << " Number of Cols: " << matrix.numCols()
               << " Number of Rows: " << matrix.numRows() << "\n";

            for (int i = 0; i < matrix.nnz(); ++i) {
                os << val_h(i) << " ";
            }
            os << "\n";

            for (int i = 0; i < matrix.nnz(); ++i) {
                os << col_h(i) << " ";
            }
            os << "\n";

            for (int i = 0; i < matrix.numRows() + 1; ++i) {
                os << row_h(i) << " ";
            }
            os << "\n";

            os.close();
            return true;
        }

        MARS_INLINE_FUNCTION
        const Integer get_nnz() const { return matrix.nnz(); }

        MARS_INLINE_FUNCTION
        const Integer get_num_cols() const { return matrix.numCols(); }

        MARS_INLINE_FUNCTION
        const Integer get_num_rows() const { return matrix.numRows(); }

        MARS_INLINE_FUNCTION
        const crs_matrix get_crs_matrix() const { return matrix; }

        MARS_INLINE_FUNCTION
        const crs_graph get_crs_graph() const { return sparsity_pattern.get_sparsity_pattern(); }

        MARS_INLINE_FUNCTION
        const SparsityPattern get_sparsity_pattern() const { return sparsity_pattern; }

        MARS_INLINE_FUNCTION
        const DofHandler &get_dof_handler() const { return sparsity_pattern.get_dof_handler(); }

    private:
        SparsityPattern sparsity_pattern;
        crs_matrix matrix;
    };
}  // namespace mars

#endif
#endif

#endif
