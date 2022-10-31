#ifndef GENERATION_MARS_DISTRIBUTED_SP_HPP_
#define GENERATION_MARS_DISTRIBUTED_SP_HPP_

#ifdef MARS_ENABLE_MPI
#ifdef MARS_ENABLE_KOKKOS_KERNELS
#include <sstream>
#include "mars_distributed_finite_element.hpp"
#include "mars_distributed_stencil.hpp"

namespace mars {

    class ISparsityPattern {
    public:
        virtual MARS_INLINE_FUNCTION ~ISparsityPattern() {}
    };

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
    // for default local/global ordinals, you can have for instantiation:
    // LO=default_lno_t
    // GO=default_size_type
    template <typename V, typename LO, typename GO, typename SHandler, typename OS = GO>
    class SparsityPattern : public ISparsityPattern {
    public:
        using Scalar = V;
        using Ordinal = LO;
        using Offset = OS;
        using Layout = default_layout;

        using device_type = typename Kokkos::Device<Kokkos::DefaultExecutionSpace, KokkosSpace>;

        using crs_matrix = typename KokkosSparse::CrsMatrix<Scalar, Ordinal, device_type, void, Offset>;

        using crs_graph = typename crs_matrix::staticcrsgraph_type;

        using crs_row = typename crs_graph::row_map_type::non_const_type;
        using crs_col = typename crs_graph::entries_type::non_const_type;
        using col_index_type = typename crs_graph::entries_type::value_type;
        using crs_value = typename crs_matrix::values_type::non_const_type;

        using DofHandler = SHandler;

        static constexpr Integer Block = DofHandler::Block;

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

            template <typename H>
            MARS_INLINE_FUNCTION Integer get_total(const Integer dof, Integer &count, const H &handler) const {
                return handler.template is_boundary<SHandler::ElemType>(dof) ? 1 : count;
            }

            // TODO:Specialize for normal stencils as in the method below since normally the boundary is excluded.
            /* template <>
            MARS_INLINE_FUNCTION Integer
            get_total<DofHandler<SHandler::Mesh, SHandler::Degree>>(const StokesStencil<DHandler> &st) const {
                return dhandler.template is_boundary<SHandler::ElemType>(dof) ? 0 : count;
            } */

            template <typename S>
            MARS_INLINE_FUNCTION void count_unique(S st) const {
                auto handler = dhandler;
                auto cntr = counter;
                st.iterate(MARS_LAMBDA(const Integer stencil_index) {
                    Integer count = 0;

                    for (int i = 0; i < st.get_length(); i++) {
                        const Integer local_dof = st.get_value(stencil_index, i);
                        if (conditional(local_dof, handler)) {
                            count++;
                        }
                    }
                    const Integer dof = st.get_value(stencil_index, 0);
                    const Integer total = get_total(dof, count, handler);
                    cntr(handler.local_to_owned_index(dof)) = total;
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

            /* TODO:specialize for stokes because normally the boundary would be excluded. */
            template <typename S>
            MARS_INLINE_FUNCTION void insert_sorted(S st) const {
                auto handler = dhandler;
                auto rp = row_ptr;
                auto ci = col_idx;
                st.iterate(MARS_LAMBDA(const Integer stencil_index) {
                    const Integer dof = st.get_value(stencil_index, 0);
                    Integer index = rp(handler.local_to_owned_index(dof));

                    if (handler.template is_boundary<SHandler::ElemType>(dof)) {
                        Integer global = handler.local_to_global(dof);
                        ci(index) = global;
                    } else {
                        Integer count = 0;

                        for (int i = 0; i < st.get_length(); i++) {
                            const Integer local_dof = st.get_value(stencil_index, i);
                            if (conditional(local_dof, handler)) {
                                const Integer global = handler.local_to_global(local_dof);
                                insert_sorted(ci, index, global, count);
                            }
                        }
                    }
                });
            }

            template <typename S>
            MARS_INLINE_FUNCTION void operator()(S &stencil, size_t I) const {
                insert_sorted<S>(stencil);
                /* //To be used with kokkos radix sort. Currently slower than the insert_sorted routine.
                insert<S>(stencil); */
            }

            /* TODO:specialize for stokes because normally the boundary would be excluded. */
            template <typename S>
            MARS_INLINE_FUNCTION void insert(S st) const {
                auto handler = dhandler;
                auto rp = row_ptr;
                auto ci = col_idx;
                st.iterate(MARS_LAMBDA(const Integer stencil_index) {
                    const Integer dof = st.get_value(stencil_index, 0);
                    Integer index = rp(handler.local_to_owned_index(dof));

                    if (handler.template is_boundary<SHandler::ElemType>(dof)) {
                        Integer global = handler.local_to_global(dof);
                        ci(index) = global;
                    } else {
                        Integer count = 0;

                        for (int i = 0; i < st.get_length(); i++) {
                            const Integer local_dof = st.get_value(stencil_index, i);
                            if (conditional(local_dof, handler)) {
                                const Integer global = handler.local_to_global(local_dof);
                                ci(index + count) = global;
                                ++count;
                            }
                        }
                    }
                });
            }

            crs_row row_ptr;
            crs_col col_idx;
            SHandler dhandler;
        };

        template <typename... ST>
        void build_pattern(ST... s) {
            using stencil_tuple = std::tuple<ST...>;

            stencil_tuple stencils(std::make_tuple(s...));

            auto global_size = get_dof_handler().get_global_dof_size();
            auto owned_size = get_dof_handler().get_owned_dof_size();

            printf("Global Dof size: %li, Owned Dof Size: %li\n", global_size, owned_size);

            ViewVectorType<Integer> counter("counter", owned_size);
            expand_tuple<CountUniqueDofs, stencil_tuple>(CountUniqueDofs(counter, get_dof_handler()), stencils);
            crs_row row_ptr("counter_scan", owned_size + 1);
            incl_excl_scan(0, owned_size, counter, row_ptr);

            auto ss = subview(row_ptr, owned_size);
            auto h_ss = create_mirror_view(ss);
            deep_copy(h_ss, ss);

            crs_col col_idx("ColIdxST", h_ss());
            expand_tuple<InsertSortedDofs, stencil_tuple>(InsertSortedDofs(row_ptr, col_idx, get_dof_handler()),
                                                          stencils);

            /* Alternatively: use the following kokkos segmented radix sort.
            For that switching to the sort function instead of insert_sorted in the InsertSortedDofs struct is needed.
            KokkosKernels::Impl::sort_crs_graph<Kokkos::DefaultExecutionSpace, crs_row, crs_col>(row_ptr, col_idx); */

            sparsity_pattern = crs_graph(col_idx, row_ptr);
        }

        // Finite Element Sparsity pattern creation

        // Unique number of dofs on theelements that share a volume, corner, face or edge.
        template <Integer LL>
        constexpr Integer label_based_node_count() const {
            constexpr Integer base = 2 * (SHandler::Degree + 1) - 1;
            constexpr Integer height = power((SHandler::Degree + 1), SHandler::ManifoldDim - 1);
            switch (LL) {
                case DofLabel::lVolume: {
                    return SHandler::elem_dofs;
                }
                case DofLabel::lCorner: {
                    return power(base, SHandler::ManifoldDim);
                }
                case DofLabel::lFace: {
                    return base * height;
                }
                case DofLabel::lEdge: {
                    // height * length * width/2
                    return base * base * height;
                }
                default: {
                    printf("Invalid Label!\n");
                    return 0;
                }
            }
        }

        MARS_INLINE_FUNCTION
        void insert_sorted_unique_fe(const ViewMatrixType<Integer> &col,
                                     const Integer row,
                                     const Integer value,
                                     Integer &count) const {
            Integer i = 0;
            while (value > col(row, i) && i < count) {
                i++;
            }

            if (value == col(row, i)) return;

            for (int j = count; j > i; --j) {
                col(row, j) = col(row, j - 1);
            }
            col(row, i) = value;
            ++count;
        }

        MARS_INLINE_FUNCTION
        bool is_unique(const ViewMatrixType<Integer> &col,
                       const Integer row,
                       const Integer value,
                       const Integer count) const {
            auto index = binary_search(&col(row, 0), 0, count, value);
            return (index == INVALID_INDEX);
        }

        template <Integer L = SHandler::dofLabel, class FE>
        ViewMatrixType<Integer> generate_fe_node_to_node_matrix(const FE &fe,
                                                                const ViewVectorType<Integer> &counter,
                                                                ViewVectorType<Integer> &locally_owned_dofs) {
            auto handler = get_dof_handler();
            const auto block = handler.get_block();

            auto node_to_element = fe.template build_node_element_dof_map<L>(locally_owned_dofs);
            auto owned_size = node_to_element.extent(0);

            Integer node_max_size = block * label_based_node_count<L>();
            ViewMatrixType<Integer> ntn("count_node_view", owned_size, node_max_size);
            Kokkos::parallel_for(
                "Count_nodes", owned_size, MARS_LAMBDA(const Integer i) {
                    for (int j = 0; j < node_max_size; j++) {
                        ntn(i, j) = INVALID_INDEX;
                    }
                });

            /* compact_owned_dofs<L>(get_dof_handler(), locally_owned_dofs); */
            assert(owned_size == block * locally_owned_dofs.extent(0));

            auto el_max_size = fe.template label_based_element_count<L>();
            Kokkos::parallel_for(
                "Count_nodes", owned_size, MARS_LAMBDA(const Integer i) {
                    /* auto owned_dof = handler.get_owned_dof(i); */
                    /* auto label = handler.get_label(owned_dof); */

                    Integer count = 0;
                    for (int j = 0; j < el_max_size; j++) {
                        auto elem_index = node_to_element(i, j);

                        if (fe.is_valid(elem_index)) {
                            for (int k = 0; k < fe.get_fe_size(); k++) {
                                auto local_dof = fe.get_elem_local_dof(elem_index, k);
                                if (local_dof > INVALID_INDEX) {
                                    const Integer global = handler.local_to_global(local_dof);
                                    insert_sorted_unique_fe(ntn, i, global, count);
                                }
                            }
                        }
                    }

                    auto base_index = handler.compute_base(i);
                    auto local_owned_dof = locally_owned_dofs(base_index);
                    auto comp_index = handler.compute_component(i);

                    auto block_owned = handler.compute_block_index(local_owned_dof, comp_index);
                    auto local_owned_index = handler.local_to_owned_dof(block_owned);
                    counter(local_owned_index) = count;
                });

            /* print_node_to_element<L, FE>(fe, locally_owned_dofs, node_to_element); */
            return ntn;
        }

        template <Integer L = SHandler::dofLabel, class FE>
        void print_node_to_element(const FE &fe,
                                   const ViewVectorType<Integer> &locally_owned_dofs,
                                   const ViewMatrixType<Integer> &node_to_element) {
            auto owned_size = node_to_element.extent(0);
            auto handler = get_dof_handler();
            auto el_max_size = fe.template label_based_element_count<L>();
            Kokkos::parallel_for(
                "print_node_elem", owned_size, MARS_LAMBDA(const Integer i) {
                    auto owned_dof = locally_owned_dofs(i);
                    auto od = handler.get_octant_from_local(owned_dof);
                    auto label = handler.get_label(owned_dof);
                    auto gid = handler.local_to_global(owned_dof);
                    for (int j = 0; j < el_max_size; j++) {
                        auto elem_index = node_to_element(i, j);

                        if (fe.is_valid(elem_index)) {
                            auto elem_sfc = fe.get_elem_sfc(elem_index);
                            auto o = handler.get_mesh_manager().get_mesh()->octant_from_sfc(elem_sfc);
                            printf("Node: %li -  %li, Label: %li, [%li, %li, %li] - octant: [%li, %li, %li]\n",
                                   owned_dof,
                                   gid,
                                   label,
                                   od.x,
                                   od.y,
                                   od.z,
                                   o.x,
                                   o.y,
                                   o.z);
                        }
                    }
                });
        }

        struct GenColIdxFromNodeToNodeTuple {
            /* Possible optmization by putting ntn to the shared memory.
            Another way is to parallelize through col indices instead of owned dof indices,
            to end up with coalesced writes instead of the current coalesced reads. */
            template <typename M, typename S>
            void generate_col_idx_from_node_to_node(const M &ntn, const S &lod) const {
                auto handler = dhandler;
                auto r = row;
                auto c = col;
                /* auto owned_size = lod.extent(0); */
                auto owned_size = ntn.extent(0);
                assert(owned_size == handler.get_block() * lod.extent(0));
                Kokkos::parallel_for(
                    "generate columns from node to node connectivity", owned_size, MARS_LAMBDA(const Integer i) {
                        auto local_owned_dof = lod(handler.compute_base(i));
                        auto comp_index = handler.compute_component(i);
                        auto block_owned = handler.compute_block_index(local_owned_dof, comp_index);
                        auto owned_index = handler.local_to_owned_dof(block_owned);

                        auto index = r(owned_index);
                        auto count = r(owned_index + 1) - index;
                        for (int j = 0; j < count; ++j) {
                            c(index + j) = ntn(i, j);
                        }
                    });
            }

            template <typename M, typename S>
            MARS_INLINE_FUNCTION void operator()(const M &ntn, const S &lod) const {
                generate_col_idx_from_node_to_node(ntn, lod);
            }

            MARS_INLINE_FUNCTION
            GenColIdxFromNodeToNodeTuple(SHandler dh, crs_col c, crs_row r) : dhandler(dh), col(c), row(r) {}

            SHandler dhandler;
            crs_col col;
            crs_row row;
        };

        /* template <typename N>
        void generate_col_idx_from_node_to_node(crs_col col, const crs_row row, const N ntn) {
            auto owned_size = get_dof_handler().get_owned_dof_size();
            Kokkos::parallel_for(
                "generate columsn from node to node connectivity", owned_size, MARS_LAMBDA(const Integer i) {
                    auto count = row(i + 1) - row(i);
                    auto index = row(i);
                    for (int j = 0; j < count; ++j) {
                        col(index + j) = ntn(i, j);
                    }
                });
        } */

        /* Split and interated in a split mode is done for reducing the memory footprint.
        When done like this different size of matrix are needed for volume, corner face edge,
        reducing in this way the matrix padding. */
        template <class FE, Integer ET = SHandler::ElemType>
        std::enable_if_t<ET == ElementType::Quad4, void> build_all_node_to_node(
            const FE &fe,
            const ViewVectorType<Integer> &counter) {
            ViewVectorType<Integer> vlod, flod, clod;

            auto vntn = generate_fe_node_to_node_matrix<DofLabel::lVolume>(fe, counter, vlod);
            auto fntn = generate_fe_node_to_node_matrix<DofLabel::lFace>(fe, counter, flod);
            auto cntn = generate_fe_node_to_node_matrix<DofLabel::lCorner>(fe, counter, clod);

            auto ntn_tuple = std::make_tuple(vntn, fntn, cntn);
            auto lod_tuple = std::make_tuple(vlod, flod, clod);

            create_sparsity_pattern(ntn_tuple, lod_tuple, counter);
        }

        template <class FE, Integer ET = SHandler::ElemType>
        std::enable_if_t<ET == ElementType::Hex8, void> build_all_node_to_node(const FE &fe,
                                                                               const ViewVectorType<Integer> &counter) {
            ViewVectorType<Integer> vlod, flod, clod, elod;

            auto vntn = generate_fe_node_to_node_matrix<DofLabel::lVolume>(fe, counter, vlod);
            auto fntn = generate_fe_node_to_node_matrix<DofLabel::lFace>(fe, counter, flod);
            auto cntn = generate_fe_node_to_node_matrix<DofLabel::lCorner>(fe, counter, clod);
            auto entn = generate_fe_node_to_node_matrix<DofLabel::lEdge>(fe, counter, elod);

            auto ntn_tuple = std::make_tuple(vntn, fntn, cntn, entn);
            auto lod_tuple = std::make_tuple(vlod, flod, clod, elod);

            create_sparsity_pattern(ntn_tuple, lod_tuple, counter);
        }

        template <class M, class S>
        void create_sparsity_pattern(M &ntn_tuple, S &lod_tuple, const ViewVectorType<Integer> &counter) {
            auto handler = get_dof_handler();

            auto global_size = get_dof_handler().get_global_dof_size();
            auto owned_size = handler.get_owned_dof_size();
            assert(owned_size == counter.extent(0));

            printf("Global Dof size: %li, Owned Dof Size: %li\n", global_size, owned_size);

            crs_row row_ptr("counter_scan", owned_size + 1);
            incl_excl_scan(0, owned_size, counter, row_ptr);

            /* Kokkos::parallel_for(
                "Count_nodes", owned_size, MARS_LAMBDA(const Integer i) {
                    auto owned_dof = handler.get_owned_dof(i);
                    auto gid = handler.local_to_global(owned_dof);
                    printf("Node: %li  scan: %li : - ", gid, row_ptr(i));
                    printf("\n");
                }); */

            auto ss = subview(row_ptr, owned_size);
            auto h_ss = create_mirror_view(ss);
            deep_copy(h_ss, ss);

            crs_col col_idx("ColIdxFe", h_ss());

            expand_tuple<GenColIdxFromNodeToNodeTuple, M, S>(
                GenColIdxFromNodeToNodeTuple(get_dof_handler(), col_idx, row_ptr), ntn_tuple, lod_tuple);

            sparsity_pattern = crs_graph(col_idx, row_ptr);
        }

        crs_matrix new_crs_matrix() const {
            auto global_size = get_dof_handler().get_global_dof_size();
            crs_value values("values", get_nnz());
            return crs_matrix("crs_matrix", global_size, values, sparsity_pattern);
        }

        /* template <Integer... Label>
        void build_pattern(FEDofMap<SHandler, Label>... fe) { */

        template <class FE>
        void build_pattern(const FE &fe) {
            /* using fe_tuple = std::tuple<ST...>; */
            /* fe_tuple fes(std::make_tuple(fe...)); */

            auto counter_size = get_dof_handler().get_owned_dof_size();
            ViewVectorType<Integer> counter("count_node_view", counter_size);
            build_all_node_to_node(fe, counter);
        }

        MARS_INLINE_FUNCTION
        void matrix_apply_constraints(const Integer row, crs_matrix m, const Scalar value) const {
            const Integer diag_row = get_dof_handler().local_to_owned_index(row);
            const col_index_type diag_col = get_dof_handler().local_to_global(row);

            const Integer row_idx = sparsity_pattern.row_map(diag_row);
            const Integer next_row_idx = sparsity_pattern.row_map(diag_row + 1) - 1;

            const Integer col_index = binary_search(sparsity_pattern.entries.data(), row_idx, next_row_idx, diag_col);

            for (int i = row_idx; i <= next_row_idx; ++i) {
                if (i == col_index) {
                    m.values(i) = value;
                } else {
                    m.values(i) = 0;
                }
            }
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
        const Integer get_col_index_from_global(const Integer row, const col_index_type col) const {
            const Integer row_idx = sparsity_pattern.row_map(row);
            const Integer next_row_idx = sparsity_pattern.row_map(row + 1) - 1;

            const Integer i = binary_search(sparsity_pattern.entries.data(), row_idx, next_row_idx, col);
            /* if (i == -1) printf("indexrow: %li, %li, %li, %li\n", row, col, row_idx, next_row_idx); */
            return i;
        }

        MARS_INLINE_FUNCTION
        const Integer get_col_index(const Integer local_row, const Integer local_col) const {
            // translate the row to the owned index and the col to the global index first.
            auto row = get_dof_handler().local_to_owned_index(local_row);
            auto col = get_dof_handler().local_to_global(local_col);

            return get_col_index_from_global(row, col);
        }

        MARS_INLINE_FUNCTION
        const Integer get_nnz() const { return get_sparsity_pattern().entries.extent(0); }

        MARS_INLINE_FUNCTION
        const SHandler &get_dof_handler() const { return dof_handler; }

    private:
        crs_graph sparsity_pattern;
        SHandler dof_handler;
    };
}  // namespace mars

#endif
#endif

#endif
