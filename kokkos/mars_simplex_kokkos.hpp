#ifndef MARS_SIMPLEX_KOKKOS_HPP
#define MARS_SIMPLEX_KOKKOS_HPP

#include "mars_base.hpp"
#include "mars_imesh_kokkos.hpp"
// #include "mars_matrix.hpp"
// #include "mars_stream.hpp"

#include <array>
#include <iostream>
#include <ostream>
#include <vector>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <initializer_list>

#include "mars_static_math_kokkos.hpp"
#include "mars_sub_view.hpp"
#include "mars_utils_kokkos.hpp"

namespace mars {

    template <Integer Dim, Integer ManifoldDim>
    class Simplex<Dim, ManifoldDim, KokkosImplementation> final : public ParallelIElem {
    public:
        static constexpr Integer ManifoldDim_ = ManifoldDim;

        static constexpr Integer ElemType = ManifoldDim_ + 1;
        static constexpr Integer NNodes = ManifoldDim_ + 1;

        using Comb = Combinations<ManifoldDim + 1, 2, KokkosImplementation>;

        ViewMatrixTextureC<Integer, Comb::value, 2> combs;

        SubView<Integer, ManifoldDim + 1> nodes;
        SubView<Integer, 2> children;  // TODO: templatize for the number of children based onthe select algorithm
        SubView<Integer, ManifoldDim + 1> side_tags;

        Integer id = INVALID_INDEX;
        // It will a view of one element basically. Think of it as pointer to the parent view.
        // It will actually get a subview.
        ViewVectorType<Integer> parent;
        Integer block = INVALID_INDEX;

        MARS_INLINE_FUNCTION
        void set_view_parent(const ViewVectorType<Integer> p) { parent = p; }

        MARS_INLINE_FUNCTION
        void set_parent_id(const Integer pid) { parent(id) = pid; }

        MARS_INLINE_FUNCTION
        Integer get_parent_id() const { return parent(id); }

        MARS_INLINE_FUNCTION Simplex() {}

        MARS_INLINE_FUNCTION Simplex(SubView<Integer, ManifoldDim + 1> n) : nodes(n) {}

        MARS_INLINE_FUNCTION Simplex(SubView<Integer, ManifoldDim + 1> n,
                                     const ViewMatrixTextureC<Integer, Comb::value, 2> &cmbs)
            : nodes(n), combs(cmbs) {}

        MARS_INLINE_FUNCTION Simplex(SubView<Integer, ManifoldDim + 1> n, SubView<Integer, 2> ch)
            : nodes(n), children(ch) {}

        MARS_INLINE_FUNCTION Simplex(SubView<Integer, ManifoldDim + 1> n,
                                     SubView<Integer, 2> ch,
                                     SubView<Integer, ManifoldDim + 1> st)
            : nodes(n), children(ch), side_tags(st) {}

        /*inline void get_nodes(std::vector<Integer> &nodes_copy) const override
        {
            nodes_copy.resize(nodes.size());
            std::copy(std::begin(nodes), std::end(nodes), std::begin(nodes_copy));
        }*/

        MARS_INLINE_FUNCTION Integer get_block() const  // override
        {
            return block;
        }

        MARS_INLINE_FUNCTION void set_block(const Integer block_id)  // override
        {
            block = block_id;
        }

        MARS_INLINE_FUNCTION Integer n_nodes() const override { return ManifoldDim + 1; }

        MARS_INLINE_FUNCTION Integer node(const Integer idx) const override {
            assert(idx < ManifoldDim + 1);
            return nodes[idx];
        }

        MARS_INLINE_FUNCTION Integer type() const override { return ManifoldDim; }

        // MARS_INLINE_FUNCTION static std::vector<Vector<Real, Dim>> &ref() {
        //     static const Integer N = ManifoldDim + 1;

        //     static std::vector<Vector<Real, Dim>> ref_;

        //     if (ref_.empty()) {
        //         ref_.resize(N);

        //         ref_[0] = Vector<Real, Dim>().zero();

        //         for (Integer i = 0; i < ManifoldDim; ++i) {
        //             ref_[i + 1] = Vector<Real, Dim>().zero();
        //             ref_[i + 1](i) = 1.;
        //         }
        //     }

        //     return ref_;
        // }

        MARS_INLINE_FUNCTION
        void edge(const Integer &edge_num,
                  Integer &v1,
                  Integer &v2) const { /*Integer vs[2];
                                       Comb::instance().choose(edge_num, nodes, vs);*/
            assert(combs(edge_num, 0) < ManifoldDim + 1);
            assert(combs(edge_num, 0) >= 0);
            v1 = nodes[combs(edge_num, 0)];
            v2 = nodes[combs(edge_num, 1)];
        }

        /*  MARS_INLINE_FUNCTION
          void side(const Integer &side_num,
                    Simplex<Dim, ManifoldDim-1,KokkosImplementation> &side) const
          {
              Combinations<ManifoldDim+1, ManifoldDim, KokkosImplementation>::choose(side_num, nodes, side.nodes);
          }*/
    };

    template <Integer Dim, Integer ManifoldDim>
    MARS_INLINE_FUNCTION constexpr static Integer n_nodes(const Simplex<Dim, ManifoldDim, KokkosImplementation> &) {
        return ManifoldDim + 1;
    }

    template <Integer Dim, Integer ManifoldDim>
    MARS_INLINE_FUNCTION constexpr static Integer n_sides(const Simplex<Dim, ManifoldDim, KokkosImplementation> &) {
        return ManifoldDim + 1;
    }

    template <Integer Dim, Integer ManifoldDim>
    MARS_INLINE_FUNCTION constexpr static Integer n_dims(const Simplex<Dim, ManifoldDim, KokkosImplementation> &) {
        return Dim;
    }

    template <Integer Dim, Integer ManifoldDim>
    MARS_INLINE_FUNCTION constexpr static Integer n_edges(const Simplex<Dim, ManifoldDim, KokkosImplementation> &) {
        return Combinations<ManifoldDim + 1, 2, KokkosImplementation>::value;
    }

}  // namespace mars

#endif  // MARS_SIMPLEX_HPP
