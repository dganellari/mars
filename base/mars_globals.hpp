#ifndef MARS_GLOBALS_HPP
#define MARS_GLOBALS_HPP

#include <vector>
#include "mars_config.hpp"
#include "mars_base.hpp"
#include <numeric>

#ifdef WITH_KOKKOS
	#define MARS_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION 
	#define MARS_LAMBDA KOKKOS_LAMBDA
	#include <Kokkos_Core.hpp>
	#include <Kokkos_UnorderedMap.hpp>
	#include "mars_device_vector.hpp"
#else
	#define MARS_INLINE_FUNCTION inline
    #define MARS_LAMBDA [=]
#endif

namespace mars {

    constexpr int hex_n_sides = 6; // 6 faces in total for the hex27.
    constexpr int hex_n_nodes = 27; // 27 nodes for the hex27.
    constexpr int hex_side_n_nodes = 9; // 9 nodes per face for the hex27.


    //FIXME not its place here
    inline void add_side(std::vector<Integer>& side, const Integer a, const Integer b,
            const Integer index)
    {
        if (a != 1 && b != 1) { //add only nodes which are not mid faces or mid edges
            side.push_back(index);
        } else if (a == 1 && b == 1) { // then add only mid faces
            side.push_back(index);
        }
    }

    MARS_INLINE_FUNCTION
    Integer index(const Integer xDim, const Integer yDim, const Integer i,
            const Integer j, const Integer k)
    {
        //return k+ (2*zDim +1) * (j + i* (2*yDim + 1));
        return i + (2 * xDim + 1) * (j + k * (2 * yDim + 1));
    }

	MARS_INLINE_FUNCTION
	Integer elem_index(const Integer i, const Integer j, const Integer k,
		const Integer xDim, const Integer yDim)
	{
		return i + (xDim+1)*(j + k*(yDim+1));
	}

	//host and device function used for the serial version as well (without kokkos).
    template<typename T>
    MARS_INLINE_FUNCTION
    void build_hex27(T&& nodes, const Integer xDim,
            const Integer yDim, const int i, const int j, const int k)
    {

        nodes[0] = index(xDim, yDim, i, j, k);
        nodes[1] = index(xDim, yDim, i + 2, j, k);
        nodes[2] = index(xDim, yDim, i + 2, j + 2, k);
        nodes[3] = index(xDim, yDim, i, j + 2, k);
        nodes[4] = index(xDim, yDim, i, j, k + 2);
        nodes[5] = index(xDim, yDim, i + 2, j, k + 2);
        nodes[6] = index(xDim, yDim, i + 2, j + 2, k + 2);
        nodes[7] = index(xDim, yDim, i, j + 2, k + 2);
        nodes[8] = index(xDim, yDim, i + 1, j, k);
        nodes[9] = index(xDim, yDim, i + 2, j + 1, k);
        nodes[10] = index(xDim, yDim, i + 1, j + 2, k);
        nodes[11] = index(xDim, yDim, i, j + 1, k);
        nodes[12] = index(xDim, yDim, i, j, k + 1);
        nodes[13] = index(xDim, yDim, i + 2, j, k + 1);
        nodes[14] = index(xDim, yDim, i + 2, j + 2, k + 1);
        nodes[15] = index(xDim, yDim, i, j + 2, k + 1);
        nodes[16] = index(xDim, yDim, i + 1, j, k + 2);
        nodes[17] = index(xDim, yDim, i + 2, j + 1, k + 2);
        nodes[18] = index(xDim, yDim, i + 1, j + 2, k + 2);
        nodes[19] = index(xDim, yDim, i, j + 1, k + 2);
        nodes[20] = index(xDim, yDim, i + 1, j + 1, k);
        nodes[21] = index(xDim, yDim, i + 1, j, k + 1);
        nodes[22] = index(xDim, yDim, i + 2, j + 1, k + 1);
        nodes[23] = index(xDim, yDim, i + 1, j + 2, k + 1);
        nodes[24] = index(xDim, yDim, i, j + 1, k + 1);
        nodes[25] = index(xDim, yDim, i + 1, j + 1, k + 2);
        nodes[26] = index(xDim, yDim, i + 1, j + 1, k + 1);
    }

    //libmesh method to map the sides to nodes.
    const std::vector<std::vector<unsigned int>> hex_side_nodes{ { 0, 3, 2,
            1, 11, 10, 9, 8, 20 }, // Side 0
            { 0, 1, 5, 4, 8, 13, 16, 12, 21 }, // Side 1
            { 1, 2, 6, 5, 9, 14, 17, 13, 22 }, // Side 2
            { 2, 3, 7, 6, 10, 15, 18, 14, 23 }, // Side 3
            { 3, 0, 4, 7, 11, 12, 19, 15, 24 }, // Side 4
            { 4, 5, 6, 7, 16, 17, 18, 19, 25 }  // Side 5
    };


	MARS_INLINE_FUNCTION
	void swap(Integer* a, Integer* b)
	{
		int t = *a;
		*a = *b;
		*b = t;
	}

	template<typename T>
	MARS_INLINE_FUNCTION const T&
	min(const T& a, const T& b)
	{
		if (b < a)
		return b;

		return a;
	}

	template<typename T>
	MARS_INLINE_FUNCTION const T&
	abs(const T& a)
	{
		if (a<0)
		return -a;

		return a;
	}

    template <typename T>
    MARS_INLINE_FUNCTION constexpr Integer power(T base, T exp) noexcept
    {
        return (exp == 0 ? 1 : base * power(base, exp - 1));
    }

    template <typename T>
    MARS_INLINE_FUNCTION constexpr Integer power_of_2(T exp)
    {
        return (exp == 0 ? 1 : 2 * power_of_2(exp - 1));
    }

    // returns the prefix sum of C
    template <typename C>
    C make_scan_index(C const &c)
    {
        static_assert(
            std::is_integral<typename C::value_type>::value,
            "make_index only applies to integral types");

        C out(c.size() + 1);
        out[0] = 0;
        std::partial_sum(c.begin(), c.end(), out.begin() + 1);
        return out;
    }

    // returns the prefix sum of C into a mirror view
    template <typename C>
    void make_scan_index_mirror(const ViewVectorType<Integer>::HostMirror &out, C const &c)
    {
        static_assert(
            std::is_integral<typename C::value_type>::value,
            "make_index only applies to integral types");

        out(0) = 0;
        std::partial_sum(c.begin(), c.end(), out.data() + 1);
    }

#ifdef WITH_KOKKOS

    template <typename T, Integer N>
	MARS_INLINE_FUNCTION int find_pivot(TempArray<T, N> &in, int start, int end)
	{
		int pivot = in[end]; // pivot
		int i = (start - 1); // Index of smaller element

		for (int j = start; j <= end - 1; j++)
		{
			// If current element is smaller than the pivot
			if (in[j] < pivot)
			{
				i++; // increment index of smaller element
				swap(&in[i], &in[j]);
			}
		}
		swap(&in[i + 1], &in[end]);
		return (i + 1);
	}

	template <typename T, Integer N>
	MARS_INLINE_FUNCTION void quick_sort(TempArray<T, N> &in, int start, int end)
	{
		if (start < end)
		{
			int pivot = find_pivot(in, start, end);

			quick_sort(in, start, pivot - 1);
			quick_sort(in, pivot + 1, end);
		}
	}

	template <typename T>
	MARS_INLINE_FUNCTION void quick_sort(TempArray<T, 2> &in, const int start,
										const int end)
	{
		if (start < end)
		{
			if (in[end] < in[start])
				swap(&in[start], &in[end]);
		}
	}

#endif
}

#endif //MARS_GLOBALS_HPP
