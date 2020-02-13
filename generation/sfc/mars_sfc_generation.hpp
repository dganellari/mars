#ifndef MARS_SFC_GENERATION_HPP
#define MARS_SFC_GENERATION_HPP

#include "mars_sfc_code.hpp"

namespace mars
{

class SFC
{
public:
    //add point functor
    struct GenerateSFC
    {

        ViewVectorType<bool> predicate;
        Integer xDim;
        Integer yDim;
        Integer zDim;

        GenerateSFC(ViewVectorType<bool> el, Integer xdm, Integer ydm) : predicate(el), xDim(xdm), yDim(ydm)
        {
        }
        GenerateSFC(ViewVectorType<bool> el, Integer xdm, Integer ydm, Integer zdm) : predicate(el), xDim(xdm), yDim(ydm), zDim(zdm)
        {
        }
        KOKKOS_INLINE_FUNCTION
        void operator()(int j, int i) const
        {
            // set to true only those elements from the vector that are generated.
            // in this way the array is already sorted and you just compact it using scan which is much faster in parallel.
            assert(encode_morton_2D(i, j) < encode_morton_2D(xDim, yDim));
            predicate(encode_morton_2D(i, j)) = 1;
        }

        /*       KOKKOS_INLINE_FUNCTION
        void operator()(int z, int y, int x) const
        {
        } */
    };

    inline bool generate_sfc(const int xDim, const int yDim, const int zDim,
                             Integer type)
    {
        using namespace Kokkos;

        switch (type)
        {
        case ElementType::Quad4:
        {
            assert(xDim != 0);
            assert(yDim != 0);
            assert(zDim == 0);

            const unsigned int n__anchor_nodes = xDim * yDim;
            reserve_elements(n__anchor_nodes);

            //calculate all range before compaction to avoid sorting.
            unsigned int allrange = encode_morton_2D(xDim, yDim); //TODO : check if enough. Test with xdim != ydim.
            ViewVectorType<bool> all_elements("predicate", allrange);
            ViewVectorType<unsigned int> scan_indices("scan_indices", allrange);

            parallel_for(
                MDRangePolicy<Rank<2>>({0, 0}, {yDim, xDim}),
                GenerateSFC(all_elements, xDim, yDim));

            //compacting the 1 and 0 array and inserting the "true" index of the all elements
            //which is the correct morton code leaving the sfc elements array sorted.
            compact_elements(scan_indices, all_elements, allrange);

            return true;
        }
        default:
        {
            return false;
        }
        }
    }

    inline void compact_elements(const ViewVectorType<unsigned int> scan_indices,
                                 const ViewVectorType<bool> all_elements, const unsigned int size)
    {
        using namespace Kokkos;

        Timer timer;

        exclusive_bool_scan(0, size, scan_indices, all_elements);

        //otherwise kokkos lambda will not work with CUDA
        ViewVectorType<unsigned int> tmp = elements_; 

        parallel_for(
            size, KOKKOS_LAMBDA(const unsigned int i) {
                if (all_elements(i) == 1)
                {
                    unsigned int k = scan_indices(i);
                    tmp(k) = i;
                    //elements_(k) =i; It will not work with CUDA. this.elements_ is a host pointer.
                }
            });
    }

    template <Integer Type>
    bool generate_sfc_elements(const Integer xDim, const Integer yDim, const Integer zDim)
    {
        bool gen_sfc = generate_sfc(xDim, yDim, zDim, Type);
        if (!gen_sfc)
        {
            std::cerr << "Not implemented for other dimensions yet" << std::endl;
        }
        return (gen_sfc);
    }

    void reserve_elements(const unsigned int n_elements)
    {
        elements_size_ = n_elements;
        elements_ = ViewVectorType<unsigned int>("morton_code", n_elements);
    }

    const ViewVectorType<unsigned int> &get_view_elements() const //override
    {
        return elements_;
    }

    unsigned int get_elem_size()
    {
        return elements_size_;
    }

private:
    ViewVectorType<unsigned int> elements_;
    unsigned int elements_size_;
};

} // namespace mars
#endif //MARS_MESH_HPP
