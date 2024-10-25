#ifndef MARS_GLOBALS_HPP
#define MARS_GLOBALS_HPP

#include <numeric>
#include <vector>
#include "mars_base.hpp"
#include "mars_config.hpp"

#ifdef MARS_ENABLE_KOKKOS
#define MARS_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
#define MARS_LAMBDA KOKKOS_LAMBDA
#define MARS_CLASS_LAMBDA KOKKOS_CLASS_LAMBDA
#include <Kokkos_Core.hpp>
#include <Kokkos_UnorderedMap.hpp>
#include "mars_device_vector.hpp"
#else
#define MARS_INLINE_FUNCTION inline
#define MARS_LAMBDA [=]
#endif

namespace mars {

template<typename T, T v>
struct integral_constant
{
    static constexpr T value = v;
    typedef T value_type;
    typedef integral_constant<T, v> type;

    MARS_INLINE_FUNCTION
    constexpr operator value_type() const noexcept { return value; } // NOLINT
};

/*! @brief A template to create structs as a type-safe version to using declarations
 * based on the https://github.com/unibas-dmi-hpc/SPH-EXA
 *
 * Used in public API functions where a distinction between different
 * arguments of the same underlying type is desired. This provides a type-safe
 * version to using declarations. Instead of naming a type alias, the name
 * is used to define a struct that inherits from StrongType<T>, where T is
 * the underlying type.
 *
 * Due to the T() conversion and assignment from T,
 * an instance of StrongType<T> struct behaves essentially like an actual T, while construction
 * from T is disabled. This makes it impossible to pass a T as a function parameter
 * of type StrongType<T>.
 */
template<class T, class Phantom>
struct StrongType
{
    using ValueType [[maybe_unused]] = T;

    //! default ctor
    constexpr MARS_INLINE_FUNCTION StrongType()
        : value_{}
    {
    }
    //! construction from the underlying type T, implicit conversions disabled
    explicit constexpr MARS_INLINE_FUNCTION StrongType(T v)
        : value_(std::move(v))
    {
    }

    //! assignment from T
    constexpr MARS_INLINE_FUNCTION StrongType& operator=(T v)
    {
        value_ = std::move(v);
        return *this;
    }

    //! conversion to T
    constexpr MARS_INLINE_FUNCTION operator T() const { return value_; } // NOLINT

    //! access the underlying value
    constexpr MARS_INLINE_FUNCTION T value() const { return value_; }

private:
    T value_;
};

/*! @brief StrongType equality comparison
 *
 * Requires that both T and Phantom template parameters match.
 * For the case where a comparison between StrongTypes with matching T, but differing Phantom
 * parameters is desired, the underlying value attribute should be compared instead
 */
template<class T, class Phantom>
constexpr MARS_INLINE_FUNCTION bool operator==(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() == rhs.value();
}

//! @brief comparison function <
template<class T, class Phantom>
constexpr MARS_INLINE_FUNCTION bool operator<(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() < rhs.value();
}

//! @brief comparison function >
template<class T, class Phantom>
constexpr MARS_INLINE_FUNCTION bool operator>(const StrongType<T, Phantom>& lhs, const StrongType<T, Phantom>& rhs)
{
    return lhs.value() > rhs.value();
}

//! @brief addition
template<class T, class Phantom>
constexpr MARS_INLINE_FUNCTION StrongType<T, Phantom> operator+(const StrongType<T, Phantom>& lhs,
                                                           const StrongType<T, Phantom>& rhs)
{
    return StrongType<T, Phantom>(lhs.value() + rhs.value());
}

//! @brief subtraction
template<class T, class Phantom>
constexpr MARS_INLINE_FUNCTION StrongType<T, Phantom> operator-(const StrongType<T, Phantom>& lhs,
                                                           const StrongType<T, Phantom>& rhs)
{
    return StrongType<T, Phantom>(lhs.value() - rhs.value());
}

/* --------------------------------------------------------------------------------------------------------- */

    constexpr int hex_n_sides = 6;       // 6 faces in total for the hex27.
    constexpr int hex_n_nodes = 27;      // 27 nodes for the hex27.
    constexpr int hex_side_n_nodes = 9;  // 9 nodes per face for the hex27.

    // FIXME not its place here
    inline void add_side(std::vector<Integer>& side, const Integer a, const Integer b, const Integer index) {
        if (a != 1 && b != 1) {  // add only nodes which are not mid faces or mid edges
            side.push_back(index);
        } else if (a == 1 && b == 1) {  // then add only mid faces
            side.push_back(index);
        }
    }

    MARS_INLINE_FUNCTION
    Integer index(const Integer xDim, const Integer yDim, const Integer i, const Integer j, const Integer k) {
        // return k+ (2*zDim +1) * (j + i* (2*yDim + 1));
        return i + (2 * xDim + 1) * (j + k * (2 * yDim + 1));
    }

    MARS_INLINE_FUNCTION
    Integer elem_index(const Integer i, const Integer j, const Integer k, const Integer xDim, const Integer yDim) {
        return i + (xDim + 1) * (j + k * (yDim + 1));
    }

    // host and device function used for the serial version as well (without kokkos).
    template <typename T>
    MARS_INLINE_FUNCTION void build_hex27(T&& nodes,
                                          const Integer xDim,
                                          const Integer yDim,
                                          const int i,
                                          const int j,
                                          const int k) {
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

    // libmesh method to map the sides to nodes.
    const std::vector<std::vector<unsigned int>> hex_side_nodes{
        {0, 3, 2, 1, 11, 10, 9, 8, 20},    // Side 0
        {0, 1, 5, 4, 8, 13, 16, 12, 21},   // Side 1
        {1, 2, 6, 5, 9, 14, 17, 13, 22},   // Side 2
        {2, 3, 7, 6, 10, 15, 18, 14, 23},  // Side 3
        {3, 0, 4, 7, 11, 12, 19, 15, 24},  // Side 4
        {4, 5, 6, 7, 16, 17, 18, 19, 25}   // Side 5
    };

    MARS_INLINE_FUNCTION
    void swap(Integer* a, Integer* b) {
        int t = *a;
        *a = *b;
        *b = t;
    }

    //! @brief This does what you think it does
    template <class T>
    MARS_INLINE_FUNCTION constexpr const T& min(const T& a, const T& b) {
        if (b < a) return b;
        return a;
    }

    //! @brief This does what you think it does
    template <class T>
    MARS_INLINE_FUNCTION constexpr const T& max(const T& a, const T& b) {
        if (a < b) return b;
        return a;
    }

    //! @brief the std version is not constexpr, this here requires two's complement
    template <class T>
    MARS_INLINE_FUNCTION constexpr std::enable_if_t<std::is_signed_v<T>, T> abs(T a) {
#ifdef __CUDA_ARCH__
        if constexpr (std::is_same_v<T, int>) {
            return ::abs(a);
        } else {
            return ::labs(a);
        }
#else
        T mask = a >> (sizeof(T) * 8 - 1);
        return (a ^ mask) - mask;
#endif
    }


    /* template <typename T>
    MARS_INLINE_FUNCTION const T abs(const T& a) {
        if (a < 0) return -a;

        return a;
    } */

    template <typename T>
    MARS_INLINE_FUNCTION constexpr Integer power(T base, T exp) noexcept {
        return (exp == 0 ? 1 : base * power(base, exp - 1));
    }

    template <typename T>
    MARS_INLINE_FUNCTION constexpr Integer power_of_2(T exp) {
        return (exp == 0 ? 1 : 2 * power_of_2(exp - 1));
    }

    // returns the prefix sum of C
    template <typename C>
    C make_scan_index(C const& c) {
        static_assert(std::is_integral<typename C::value_type>::value, "make_index only applies to integral types");

        C out(c.size() + 1);
        out[0] = 0;
        std::partial_sum(c.begin(), c.end(), out.begin() + 1);
        return out;
    }

#ifdef MARS_ENABLE_KOKKOS

    template <typename T, Integer N>
    MARS_INLINE_FUNCTION int find_pivot(TempArray<T, N>& in, int start, int end) {
        int pivot = in[end];  // pivot
        int i = (start - 1);  // Index of smaller element

        for (int j = start; j <= end - 1; j++) {
            // If current element is smaller than the pivot
            if (in[j] < pivot) {
                i++;  // increment index of smaller element
                swap(&in[i], &in[j]);
            }
        }
        swap(&in[i + 1], &in[end]);
        return (i + 1);
    }

    template <typename T, Integer N>
    MARS_INLINE_FUNCTION void quick_sort(TempArray<T, N>& in, int start, int end) {
        if (start < end) {
            int pivot = find_pivot(in, start, end);

            quick_sort(in, start, pivot - 1);
            quick_sort(in, pivot + 1, end);
        }
    }

    template <typename T>
    MARS_INLINE_FUNCTION void quick_sort(TempArray<T, 2>& in, const int start, const int end) {
        if (start < end) {
            if (in[end] < in[start]) swap(&in[start], &in[end]);
        }
    }

#endif
}  // namespace mars

#endif  // MARS_GLOBALS_HPP
