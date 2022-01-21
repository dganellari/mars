#ifndef MARS_I_MESH_KOKKOS_HPP
#define MARS_I_MESH_KOKKOS_HPP

#include <algorithm>
#include <array>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>
#include "mars_point.hpp"
#include "mars_sub_view.hpp"

namespace mars {

    class ParallelIElem {
    public:
        virtual MARS_INLINE_FUNCTION ~ParallelIElem() {}
        // virtual void get_nodes(std::vector<Integer> &nodes) const = 0;
        virtual MARS_INLINE_FUNCTION Integer n_nodes() const = 0;
        virtual MARS_INLINE_FUNCTION Integer node(const Integer idx) const = 0;
        /*   virtual MARS_INLINE_FUNCTION Integer get_block() const = 0;
          virtual MARS_INLINE_FUNCTION void set_block(const Integer block_id) = 0; */

        // for the moment it just returns the simplex type (e.g. 2 for triangles)
        virtual MARS_INLINE_FUNCTION Integer type() const = 0;
    };

    /*MARS_INLINE_FUNCTION Integer n_nodes(const ParallelIElem &elem)
    {
            return elem.n_nodes();
    }*/

    class DimensionlessIMesh {
    public:
        virtual MARS_INLINE_FUNCTION ~DimensionlessIMesh() {}
    };

    template <Integer Dim_>
    class ParallelIMesh : public DimensionlessIMesh {
    public:
        static const Integer Dim = Dim_;
        using Elem = mars::ParallelIElem;

        virtual MARS_INLINE_FUNCTION ~ParallelIMesh() {}
        // virtual void points(const Integer id, std::vector<Point> &pts) const = 0;
        /*virtual Elem &elem(const Integer id) = 0;
        virtual const Elem &elem(const Integer id) const = 0;
*/
        /*virtual MARS_INLINE_FUNCTION Elem elem(const Integer id) = 0;
        virtual MARS_INLINE_FUNCTION const Elem elem(const Integer id) const = 0;*/
        virtual MARS_INLINE_FUNCTION bool is_active(const Integer id) const = 0;
        virtual MARS_INLINE_FUNCTION Integer n_nodes() const = 0;
        virtual MARS_INLINE_FUNCTION Integer n_elements() const = 0;
        // virtual MARS_INLINE_FUNCTION Integer n_active_elements() const = 0;
        // virtual Integer add_point(const Point<Real,Dim> &point) = 0;
        virtual MARS_INLINE_FUNCTION Point<Real, Dim> point(const Integer i) = 0;
        virtual MARS_INLINE_FUNCTION const Point<Real, Dim> point(const Integer i) const = 0;
        // virtual Integer add_elem(const IElem &elem) = 0;
        // virtual Integer add_elem(const std::vector<Integer> &nodes, const int row) = 0;
        // virtual const std::vector<Point> &points() const = 0;
        // virtual ViewMatrixType<Integer> get_view_elements() const =0;

        virtual MARS_INLINE_FUNCTION Integer type() const = 0;
        virtual void reserve(const std::size_t n_elements, const std::size_t n_points) = 0;
    };

}  // namespace mars

#endif  // MARS_I_MESH_HPP
