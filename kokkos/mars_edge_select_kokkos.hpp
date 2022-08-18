#ifndef MARS_EDGE_SELECT_KOKKOS_HPP
#define MARS_EDGE_SELECT_KOKKOS_HPP

#include "mars_base.hpp"
#ifdef WITH_KOKKOS
#include "mars_edge_kokkos.hpp"
#endif
#include "mars_fwd.hpp"
#include "mars_globals.hpp"

#include <ostream>

namespace mars {

    template <class Mesh>
    class EdgeSelect<Mesh, KokkosImplementation> {
    public:
        using Edge = typename Mesh::Edge;

        static const Integer Dim = Mesh::Dim;
        static const Integer ManifoldDim = Mesh::ManifoldDim;

        virtual MARS_INLINE_FUNCTION ~EdgeSelect() {}

        virtual MARS_INLINE_FUNCTION Integer select(const Mesh &mesh, const Integer element_id) const {
            // first edge selected
            return 0;
        }

        virtual MARS_INLINE_FUNCTION Integer select(const Mesh &mesh,
                                                    const Edge &neighbor_edge,
                                                    const Integer element_id) const {
            // first edge selected
            return 0;
        }

        // virtual MARS_INLINE_FUNCTION Integer select(const Mesh &mesh,
        //                                             const Integer element_id,
        //                                             const EdgeElementMap &edge_element_map) const {
        //     // first edge selected
        //     return 0;
        // }

        virtual MARS_INLINE_FUNCTION Integer stable_select(const Mesh &mesh, const Integer element_id) const {
            // first edge selected
            return 0;
        }

        virtual bool repair_element() { return true; }

        virtual bool can_refine(const Mesh &mesh, const Integer element_id) const { return true; }

        virtual void reorder_edge(const Mesh &mesh, const Integer element_id, Integer &v1, Integer &v2) const {}

        // virtual void edge_refined(const Mesh &mesh, const EdgeElementMap &eem, const Edge &edge) {
        //     //
        // }

        virtual void update(const Mesh &mesh) {}

        virtual void element_refined(const Mesh &mesh,
                                     const Integer element_id,
                                     const Edge &edge,
                                     const Integer local_midpoint_id) {
            // do smth
        }

        virtual MARS_INLINE_FUNCTION bool is_recursive() const { return false; }

        virtual void describe(std::ostream &os) const {}

        virtual std::string name() const = 0;
    };

}  // namespace mars

#endif  // MARS_EDGE_SELECT_HPP
