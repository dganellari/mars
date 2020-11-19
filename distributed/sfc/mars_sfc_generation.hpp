#ifndef MARS_SFC_GENERATION_HPP
#define MARS_SFC_GENERATION_HPP

#include "mars_sfc_code.hpp"

namespace mars {
    template <Integer Type>
    class SFC {
    public:
        inline void compact_elements(const ViewVectorType<bool> all_elements) {
            using namespace Kokkos;

            exclusive_bool_scan(0, get_all_range(), get_view_sfc_to_local(), all_elements);

            auto sfc_subview = subview(get_view_sfc_to_local(), get_all_range() - 1);
            auto elm_subview = subview(all_elements, get_all_range() - 1);
            auto h_sfc = create_mirror_view(sfc_subview);
            auto h_elm = create_mirror_view(elm_subview);
            deep_copy(h_sfc, sfc_subview);
            deep_copy(h_elm, elm_subview);

            const Integer elem_size = h_sfc() + h_elm();
            reserve_elements(elem_size);

            // otherwise kokkos lambda will not work with CUDA
            ViewVectorType<Integer> el = elements_;
            ViewVectorType<Integer> sfc_to_local = sfc_to_local_;
            parallel_for(
                get_all_range(), KOKKOS_LAMBDA(const Integer i) {
                    if (all_elements(i) == 1) {
                        Integer k = sfc_to_local(i);
                        el(k) = i;
                        // elements_(k) =i; It will not work with CUDA. this.elements_ is a host pointer.
                    }
                });
        }

        inline void compact_element_labels(const ViewVectorType<bool> all_elements,
                                           const ViewVectorType<bool> all_labels) {
            using namespace Kokkos;

            reserve_element_labels(get_elem_size());

            // otherwise kokkos lambda will not work with CUDA
            ViewVectorType<Integer> elb = element_labels_;
            ViewVectorType<Integer> sfc_to_local = sfc_to_local_;
            parallel_for(
                get_all_range(), KOKKOS_LAMBDA(const Integer i) {
                    if (all_elements(i) == 1) {
                        Integer k = sfc_to_local(i);
                        elb(k) = all_labels(i);
                        // elements_(k) =i; It will not work with CUDA. this.elements_ is a host pointer.
                    }
                });
        }

        inline void compact_element_and_labels(const ViewVectorType<bool> all_elements,
                                               const ViewVectorType<bool> all_labels) {
            compact_elements(all_elements);  // this should come first!
            compact_element_labels(all_elements, all_labels);
        }

        struct GenerateSFC {
            ViewVectorType<bool> predicate;
            Integer xDim;
            Integer yDim;
            Integer zDim;

            GenerateSFC(ViewVectorType<bool> el, Integer xdm, Integer ydm) : predicate(el), xDim(xdm), yDim(ydm) {}
            GenerateSFC(ViewVectorType<bool> el, Integer xdm, Integer ydm, Integer zdm)
                : predicate(el), xDim(xdm), yDim(ydm), zDim(zdm) {}
            KOKKOS_INLINE_FUNCTION
            void operator()(int j, int i) const {
                // set to true only those elements from the vector that are generated.
                // in this way the array is already sorted and you just compact it using scan which is much faster in
                // parallel.
                assert(encode_morton_2D(i, j) < encode_morton_2D(xDim, yDim));
                if (encode_morton_2D(i, j) >= encode_morton_2D(xDim, yDim)) {
                    const Integer chunk_size = xDim * yDim;
                    printf("You have reached the mesh genration limit size. Can not generate mesh %li elements\n",
                           chunk_size);
                    exit(1);
                }

                predicate(encode_morton_2D(i, j)) = 1;
            }

            KOKKOS_INLINE_FUNCTION
            void operator()(int k, int j, int i) const {
                // set to true only those elements from the vector that are generated.
                // in this way the array is already sorted and you just compact it using scan which is much faster in
                // parallel.
                assert(encode_morton_3D(i, j, k) < encode_morton_3D(xDim, yDim, zDim));
                if (encode_morton_3D(i, j, k) >= encode_morton_3D(xDim, yDim, zDim)) {
                    const Integer chunk_size = xDim * yDim * zDim;
                    printf("You have reached the mesh genration limit size. Can not generate mesh %li elements\n",
                           chunk_size);
                    exit(1);
                }

                predicate(encode_morton_3D(i, j, k)) = 1;
            }
        };

        inline bool generate_sfc() {
            using namespace Kokkos;

            ViewVectorType<bool> all_elements("predicate", get_all_range());
            Integer n__anchor_nodes;

            switch (Type) {
                case ElementType::Quad4: {
                    assert(xDim != 0);
                    assert(yDim != 0);
                    assert(zDim == 0);

                    n__anchor_nodes = xDim * yDim;

                    parallel_for(MDRangePolicy<Rank<2>>({0, 0}, {yDim, xDim}), GenerateSFC(all_elements, xDim, yDim));
                }
                case ElementType::Hex8: {
                    assert(xDim != 0);
                    assert(yDim != 0);
                    assert(zDim != 0);

                    n__anchor_nodes = xDim * yDim * zDim;
                    parallel_for(MDRangePolicy<Rank<3>>({0, 0, 0}, {zDim, yDim, xDim}),
                                 GenerateSFC(all_elements, xDim, yDim, zDim));
                }
                default: {
                    return false;
                }
            }
            // compacting the 1 and 0 array and inserting the "true" index of the all elements
            // which is the correct morton code leaving the sfc elements array sorted.
            compact_elements(all_elements);
            assert(n__anchor_nodes == get_elem_size());

            return true;
        }

        bool generate_sfc_elements() {
            bool gen_sfc = generate_sfc();
            if (!gen_sfc) {
                std::cerr << "Not implemented for other dimensions yet" << std::endl;
            }
            return (gen_sfc);
        }

        void reserve_element_labels(const Integer n_elements) {
            element_labels_ = ViewVectorType<Integer>("labels", n_elements);
        }

        void reserve_elements(const Integer n_elements) {
            elements_size_ = n_elements;
            elements_ = ViewVectorType<Integer>("morton_code", n_elements);
        }

        const ViewVectorType<Integer> &get_view_element_labels() const  // override
        {
            return element_labels_;
        }

        const ViewVectorType<Integer> &get_view_elements() const  // override
        {
            return elements_;
        }

        void set_elem_size(const Integer size) { elements_size_ = size; }

        const Integer get_elem_size() const { return elements_size_; }

        const ViewVectorType<Integer> &get_view_sfc_to_local() const  // override
        {
            return sfc_to_local_;
        }

        Integer get_all_range() { return all_range_; }

        MARS_INLINE_FUNCTION
        SFC() = default;

        MARS_INLINE_FUNCTION
        SFC(const int x, const int y, const int z) : xDim(x), yDim(y), zDim(z) {
            switch (Type) {
                case ElementType::Quad4: {
                    all_range_ = encode_morton_2D(xDim + 1, yDim + 1);
                    break;
                }
                case ElementType::Hex8: {
                    all_range_ = encode_morton_3D(xDim + 1, yDim + 1, zDim + 1);
                    break;
                }
            }

            sfc_to_local_ = ViewVectorType<Integer>("sfc_to_local_", all_range_);
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
        Integer elements_size_;

        ViewVectorType<Integer> sfc_to_local_;
        Integer all_range_;

        Integer xDim, yDim, zDim;
    };

}  // namespace mars
#endif
