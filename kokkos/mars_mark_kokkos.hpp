#ifndef MARS_MARK_KOKKOS_HPP_
#define MARS_MARK_KOKKOS_HPP_

#include "mars_base.hpp"

#include "mars_mesh_kokkos.hpp"

namespace mars {

    template <class Mesh_>
    MARS_INLINE_FUNCTION static bool in_hypersphere(const Mesh_* mesh,
                                                    const Integer element_id,
                                                    const ViewVectorTypeC<Real, Mesh_::Dim> center,
                                                    const Real radius) {
        bool inside = false;
        bool outside = false;

        if (mesh->is_active(element_id))  // TODO: remove if by compacting on active elems.
        {
            auto e = mesh->elem(element_id);

            for (Integer i = 0; i < Mesh_::ManifoldDim + 1; ++i) {
                auto dir = mesh->point(e.nodes[i]) - center;
                auto d = dir.norm();

                if (d < radius) {
                    inside = true;
                } else if (d > radius) {
                    outside = true;
                } else if (abs(d) < 1e-16) {
                    inside = true;
                    outside = true;
                    break;
                }
            }
        }

        // return inside && outside;
        return inside;
    }

    template <class Mesh_>
    struct Hypersphere {
        static constexpr Integer Dim = Mesh_::Dim;
        static constexpr Integer ManifoldDim = Mesh_::ManifoldDim;

        Mesh_* mesh;
        ViewVectorTypeC<Real, Dim> center;
        Real radius;

        ViewVectorType<Integer> count;

        Hypersphere(Mesh_* ms, ViewVectorTypeC<Real, Dim> c, Real rd, ViewVectorType<Integer> pc)
            : mesh(ms), center(c), radius(rd), count(pc) {}

        Hypersphere(Mesh_* ms, ViewVectorTypeC<Real, Dim> c, Real rd) : mesh(ms), center(c), radius(rd) {}

        MARS_INLINE_FUNCTION
        void operator()(const int element_id) const {
            Integer hyper_count = 0;

            if (in_hypersphere(mesh, element_id, center, radius)) {
                ++hyper_count;
            }

            count(element_id + 1) = hyper_count;
            //+1 for leaving the first cell 0 and performing an inclusive scan on the rest
            // to have both exclusive and inclusive and the total on the last cell.
        }
    };

    template <class Mesh_>
    struct ScatterHypersphere : Hypersphere<Mesh_> {
        static constexpr Integer Dim = Mesh_::Dim;
        static constexpr Integer ManifoldDim = Mesh_::ManifoldDim;

        ViewVectorType<Integer> elements;
        ViewVectorType<Integer> scan;

        ScatterHypersphere(Mesh_* ms,
                           ViewVectorTypeC<Real, Dim> c,
                           Real rd,
                           ViewVectorType<Integer> sc,
                           ViewVectorType<Integer> elems)
            : Hypersphere<Mesh_>(ms, c, rd), scan(sc), elements(elems) {}

        MARS_INLINE_FUNCTION
        void operator()(const int element_id) const {
            if (in_hypersphere(this->mesh, element_id, this->center, this->radius)) {
                elements(scan(element_id)) = element_id;
            }
        }
    };

    template <class Mesh_>
    ViewVectorType<Integer> mark_hypersphere_for_refinement(Mesh_* mesh,
                                                            const ViewVectorTypeC<Real, Mesh_::Dim> center,
                                                            const Real radius,
                                                            const Integer nr_elements) {
        using namespace Kokkos;

        Timer timer1;

        ViewVectorType<Integer> count_ = ViewVectorType<Integer>("count_hypers", nr_elements + 1);

        parallel_for(nr_elements, Hypersphere<Mesh_>(mesh, center, radius, count_));

        double time1 = timer1.seconds();
        // std::cout << "Hypersphere Count took: " << time1 << " seconds." << std::endl;

        Timer timer2;

        inclusive_scan(1, nr_elements + 1, count_);

        double time2 = timer2.seconds();
        // std::cout << "Inclusive_scan took: " << time2 << " seconds." << std::endl;

        auto index_subview = subview(count_, nr_elements);
        auto h_ic = create_mirror_view(index_subview);

        // Deep copy device view to host view.
        deep_copy(h_ic, index_subview);

        ViewVectorType<Integer> elements = ViewVectorType<Integer>("hyper_elements", h_ic());

        Timer timer3;

        parallel_for(nr_elements, ScatterHypersphere<Mesh_>(mesh, center, radius, count_, elements));

        double time3 = timer3.seconds();
        // std::cout << "ScatterHypersphere took: " << time3 << " seconds." << std::endl;

        fence();

        return elements;
    }

    // mesh is a device pointer and can not be used here for getting the nr_elements.
    inline ViewVectorType<Integer> mark_active(const ViewVectorType<bool> active) {
        using namespace Kokkos;

        const Integer nr_elements = active.extent(0);
        ViewVectorType<Integer> scan = ViewVectorType<Integer>("count_hypers", nr_elements + 1);

        incl_excl_scan(0, nr_elements, active, scan);

        auto index_subview = subview(scan, nr_elements);
        auto h_ic = create_mirror_view(index_subview);

        // Deep copy device view to host view.
        deep_copy(h_ic, index_subview);

        Timer timer3;

        ViewVectorType<Integer> elements = ViewVectorType<Integer>("hyper_elements", h_ic());

        parallel_for(
            nr_elements, KOKKOS_LAMBDA(const Integer element_id) {
                if (active(element_id)) {
                    elements(scan(element_id)) = element_id;
                }
            });

        double time3 = timer3.seconds();

        return elements;
    }

    template <class Mesh_>
    bool compact(const Mesh_* mesh, ViewVectorType<Integer>& elements) {
        using namespace Kokkos;

        Timer timer3;

        const Integer nr_elements = elements.extent(0);
        // printf("scanI nr_elements: %li\n", nr_elements);

        ViewVectorType<bool> predicate("predicate", nr_elements);

        // build predicate
        parallel_for(
            nr_elements, KOKKOS_LAMBDA(const Integer i) {
                if (mesh->is_active(elements(i))) {
                    predicate(i) = true;
                }
            });

        ViewVectorType<Integer> scan = ViewVectorType<Integer>("count_hypers", nr_elements + 1);
        // scan the predicate for getting the addresses for compaction
        incl_excl_scan(0, nr_elements, predicate, scan);

        auto index_subview = subview(scan, nr_elements);
        auto h_ic = create_mirror_view(index_subview);

        // Deep copy device view to host view.
        deep_copy(h_ic, index_subview);
        // std::cout << "Scan for compact count (new nr_elements): " << h_ic(0)<< std::endl;

        ViewVectorType<Integer> el = ViewVectorType<Integer>("hyper_elements", h_ic());
        // compact using the addresses in scan
        parallel_for(
            nr_elements, KOKKOS_LAMBDA(const Integer i) {
                if (predicate(i)) {
                    el(scan(i)) = elements(i);
                }
            });

        elements = el;

        double time3 = timer3.seconds();
        // std::cout << "compact elements took: " << time3 << " seconds." << std::endl;

        return elements.extent(0) > 0;
    }

    template <typename T, Integer Dim>
    void fill_view(ViewVectorTypeC<T, Dim> point, const T value) {
        using namespace Kokkos;

        typename ViewVectorTypeC<T, Dim>::HostMirror h_pt = create_mirror_view(point);

        for (Integer i = 0; i < Dim; ++i) h_pt(i) = value;

        deep_copy(point, h_pt);
    }

}  // namespace mars
#endif /* GENERATION_MARS_UTILS_KOKKOS_HPP_ */
