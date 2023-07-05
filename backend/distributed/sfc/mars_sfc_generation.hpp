#ifndef MARS_SFC_GENERATION_HPP
#define MARS_SFC_GENERATION_HPP

#include "mars_sfc_code.hpp"

namespace mars {

    template <Integer Type>
    MARS_INLINE_FUNCTION Integer compute_all_range(const Integer x, const Integer y, const Integer z) {
        switch (Type) {
            case ElementType::Quad4: {
                return encode_morton_2D(x + 1, y + 1);
            }
            case ElementType::Hex8: {
                return encode_morton_3D(x + 1, y + 1, z + 1);
            }
            default: {
                return INVALID_INDEX;
            }
        }
    }

    template <Integer Type>
    class SFC {
    public:
        inline void compact_elements(const ViewVectorType<Integer> &sfc_to_local,
                                     const ViewVectorType<bool> all_elements) {
            using namespace Kokkos;

            exclusive_bool_scan(0, get_all_range(), sfc_to_local, all_elements);

            auto sfc_subview = subview(sfc_to_local, get_all_range() - 1);
            auto elm_subview = subview(all_elements, get_all_range() - 1);
            auto h_sfc = create_mirror_view(sfc_subview);
            auto h_elm = create_mirror_view(elm_subview);
            deep_copy(h_sfc, sfc_subview);
            deep_copy(h_elm, elm_subview);

            const Integer elem_size = h_sfc() + h_elm();
            reserve_elements(elem_size);

            // otherwise kokkos lambda will not work with CUDA
            ViewVectorType<Integer> el = elements_;
            parallel_for(
                get_all_range(), KOKKOS_LAMBDA(const Integer i) {
                    if (all_elements(i) == 1) {
                        Integer k = sfc_to_local(i);
                        el(k) = i;
                        // elements_(k) =i; It will not work with CUDA. this.elements_ is a host pointer.
                    }
                });
        }

        inline void compact_elements(const ViewVectorType<Integer> &sfc_to_local,
                                     const ViewVectorType<bool> all_elements,
                                     const ViewVectorType<Integer> all_labels,
                                     ViewVectorType<Integer> &elb) {
            using namespace Kokkos;
            elb = ViewVectorType<Integer>("compacted_elb", get_elem_size());
            // otherwise kokkos lambda will not work with CUDA
            parallel_for(
                get_all_range(), KOKKOS_LAMBDA(const Integer i) {
                    if (all_elements(i) == 1) {
                        Integer k = sfc_to_local(i);
                        elb(k) = all_labels(i);
                        // elements_(k) =i; It will not work with CUDA. this.elements_ is a host pointer.
                    }
                });
        }

        inline void compact_element_and_labels(const ViewVectorType<Integer> &sfc_to_local,
                                               const ViewVectorType<bool> all_elements,
                                               const ViewVectorType<Integer> all_labels) {
            compact_elements(sfc_to_local, all_elements);  // this should come first!
            /* compact_element_labels(all_elements, all_labels); */
            compact_elements(sfc_to_local, all_elements, all_labels, element_labels_);
        }

        void init_element_orientations(const Integer n_elements) {
            // reserve the view then init
            element_orientations_ = ViewVectorType<Integer>("Orientation", get_elem_size());

            // remove the shallow copy when switching to c++17.
            // Currently needed for the cuda code.
            auto orient = element_orientations_;
            Kokkos::parallel_for(
                "init_orientation", get_elem_size(), MARS_LAMBDA(const Integer i) { orient(i) = DofOrient::yDir; });
        }

        void reserve_element_labels(const Integer n_elements) {
            element_labels_ = ViewVectorType<Integer>("labels", n_elements);
        }

        void reserve_elements(const Integer n_elements) {
            elements_size_ = n_elements;
            elements_ = ViewVectorType<Integer>("morton_code", n_elements);
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_element_orientations() const  // override
        {
            return element_orientations_;
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_element_labels() const  // override
        {
            return element_labels_;
        }

        MARS_INLINE_FUNCTION
        const ViewVectorType<Integer> &get_view_elements() const  // override
        {
            return elements_;
        }

        MARS_INLINE_FUNCTION
        Integer get_sfc(const Integer i) const  // override
        {
            return elements_(i);
        }

        MARS_INLINE_FUNCTION
        Integer get_orientation(const Integer i) const  // override
        {
            return element_orientations_(i);
        }

        MARS_INLINE_FUNCTION
        void set_label(const Integer i, const Integer label) const  // override
        {
            element_labels_(i) = label;
        }

        MARS_INLINE_FUNCTION
        Integer get_label(const Integer i) const  // override
        {
            return element_labels_(i);
        }

        void set_elem_size(const Integer size) { elements_size_ = size; }

        MARS_INLINE_FUNCTION
        const Integer get_elem_size() const { return elements_size_; }

        MARS_INLINE_FUNCTION
        const Integer get_all_range() const { return all_range_; }

        SFC() = default;

        SFC(const Integer x, const Integer y, const Integer z) : xDim(x), yDim(y), zDim(z) {
            all_range_ = compute_all_range<Type>(xDim, yDim, zDim);
        }

        MARS_INLINE_FUNCTION
        const UnorderedMap<Integer, Integer> &get_sfc_to_local_map() const { return sfc_to_local_map_; }

        void generate_sfc_to_local_map() {
            sfc_to_local_map_ = UnorderedMap<Integer, Integer>(get_elem_size());
            // copies needed because of cuda lambda functions. Issue: Copied as *this.{} host pointer.
            auto element_view = elements_;
            auto sfc_map = sfc_to_local_map_;
            Kokkos::parallel_for(
                "generate_sfc_to_local_map", get_elem_size(), MARS_LAMBDA(const Integer i) {
                    sfc_map.insert(element_view(i), i);
                });
        }

        MARS_INLINE_FUNCTION
        Integer get_XDim() const { return xDim; }

        MARS_INLINE_FUNCTION
        Integer get_YDim() const { return yDim; }

        MARS_INLINE_FUNCTION
        Integer get_ZDim() const { return zDim; }

    private:
        ViewVectorType<Integer> elements_;
        ViewVectorType<Integer> element_labels_;
        ViewVectorType<Integer> element_orientations_;
        Integer elements_size_;

        UnorderedMap<Integer, Integer> sfc_to_local_map_;
        Integer all_range_;

        Integer xDim, yDim, zDim;
    };
}  // namespace mars
#endif
