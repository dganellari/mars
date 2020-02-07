#ifndef MARS_SFC_GENERATION_HPP
#define MARS_SFC_GENERATION_HPP

#include "mars_mesh_kokkos.hpp"
#include "mars_sfc_code.hpp"

namespace mars
{

class SFC
{
    //add point functor
    struct GenerateSFC
    {

        ViewVectorType<bool> elems;
        Integer xDim;
        Integer yDim;
        Integer zDim;

        GenerateSFC(ViewVectorType<bool> el, Integer xdm, Integer ydm) : 
		elems(el), xDim(xdm), yDim(ydm)
        {
        }

        GenerateSFC(ViewVectorType<bool> el, Integer xdm, Integer ydm, Integer zdm) : 
		elems(el), xDim(xdm), yDim(ydm), zDim(zdm)
        {
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(int j, int i) const
        {
 			// set to true only those elements from the vector that are generated.
			// in this way the array is already sorted and you just compact it using scan which is much faster in parallel.
			elems(encode_morton_2D(i, j)) = 1;
        }

  /*       KOKKOS_INLINE_FUNCTION
        void operator()(int z, int y, int x) const
        {

            int index = (xDim + 1) * (yDim + 1) * z + (xDim + 1) * y + x;

            points(index, 0) = static_cast<Real>(x) / static_cast<Real>(xDim);
            points(index, 1) = static_cast<Real>(y) / static_cast<Real>(yDim);
            points(index, 2) = static_cast<Real>(z) / static_cast<Real>(zDim);
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

            const int n__anchor_nodes = xDim * yDim;
            reserve_elements(n__anchor_nodes);

            //calculate all range before compaction to avoid sorting.
            unsigned int allrange = encode_morton_2D(xDim, yDim); //TODO : check if enough. Test with xdim != ydim.
            ViewVectorType<bool> all_elements_ = ViewVectorType<bool>("allelems", allrange);

            parallel_for(
                MDRangePolicy<Rank<2>>({0, 0}, {yDim, xDim}),
                GenerateSFC(all_elements_, xDim, yDim));

			//compacting the 1 and 0 array and inserting the "true" index of the all elements 
			//which is the correct morton code leaving the sfc elements array sorted.
			compact_elements(all_elements_, allrange);

            return true;
        }
        default:
        {
            return false;
        }
        }
    }

    inline void compact_elements(ViewVectorType<bool> &all_elements_, const int size)
    {
		using namespace Kokkos;

		Timer timer;

		exclusive_scan(0, size, all_elements_);

		parallel_for(size, KOKKOS_LAMBDA(uint32_t i)
		{
			if(all_elements_(i) == 1)
			{
				Integer k = all_elements_(i);
				elements_(k) = i;
			}
		});
    }

    template <Integer Type>
    bool generate_sfc_elements(const Integer xDim, const Integer yDim, const Integer zDim)
    {

        assert(Dim <= 3);
        assert(ManifoldDim <= Dim);

        bool gen_sfc = generate_sfc(xDim, yDim, zDim, Type);

        if (!gen_sfc)
            std::cerr << "Not implemented for other dimensions yet" << std::endl;

        return (gen_sfc);
    }

    void reserve_elements(const std::size_t n_elements)
    {
        elements_size_ = n_elements;
        elements_ = ViewVectorType<Integer>("elems", n_elements);
    }

private:
    ViewVectorType<Integer> elements_;
    Integer elements_size_;
}

} // namespace mars
#endif //MARS_MESH_HPP
