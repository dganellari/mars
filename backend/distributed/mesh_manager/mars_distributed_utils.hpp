#ifndef MARS_DISTRIBUTED_UTILS_HPP
#define MARS_DISTRIBUTED_UTILS_HPP

#include <ostream>
#include <tuple>
#include <type_traits>
#include "mars_globals.hpp"

#ifdef MARS_ENABLE_KOKKOS_KERNELS
#include "Kokkos_ArithTraits.hpp"
#ifdef MARS_ENABLE_KOKKOS
#include "Kokkos_Atomic.hpp"
#include "Kokkos_Macros.hpp"
#include "mars_utils_kokkos.hpp"
#if KOKKOS_VERSION >= 30500
#include "impl/Kokkos_Atomic_Generic.hpp"
#endif
#endif

#ifdef MARS_ENABLE_MPI
#include <mpi.h>
#endif

namespace mars {

    inline void Abort() {
        int error_code = -1;
#ifdef MARS_ENABLE_MPI
        MPI_Abort(MPI_COMM_WORLD, error_code);
#else
        exit(error_code);
#endif
    }

    inline void Abort(const std::string& message) {
        std::printf("%s\n", message.c_str());
        Abort();
    }

    /******************** Tuple utils *****************************/

    /* getting the index of a type in a variadic template definition */

    template <class Tuple, class T, std::size_t Index = 0>
    struct TypeIdx;

    template <std::size_t Index, bool Valid>
    struct TypeIdxTest : public std::integral_constant<std::size_t, Index> {};

    template <std::size_t Index>
    struct TypeIdxTest<Index, false> {
        static_assert(Index < 0, "Type not found in the tuple!");
    };

    template <class Head, class T, std::size_t Index>
    struct TypeIdx<std::tuple<Head>, T, Index> : public TypeIdxTest<Index, std::is_same<Head, T>::value> {};

    template <class Head, class... Rest, class T, std::size_t Index>
    struct TypeIdx<std::tuple<Head, Rest...>, T, Index>
        : public std::conditional<std::is_same<Head, T>::value,
                                  std::integral_constant<std::size_t, Index>,
                                  TypeIdx<std::tuple<Rest...>, T, Index + 1>>::type {};

    /* ------------------------------------------------------------------------------------------------ */
    // getting the  nth type of a variadic.
    template <std::size_t I, typename... T>
    struct NthType;

    // recursive case
    template <std::size_t I, typename Head, typename... Tail>
    struct NthType<I, Head, Tail...> {
        typedef typename NthType<I - 1, Tail...>::type type;
    };

    template <class Head, class... Tail>
    struct NthType<0, Head, Tail...> {
        typedef Head type;
    };

    // getting the  nth type of a variadic.
    template <std::size_t I, Integer... T>
    struct NthValue;

    // recursive case
    template <std::size_t I, Integer Head, Integer... Tail>
    struct NthValue<I, Head, Tail...> {
        static constexpr Integer value = NthValue<I - 1, Tail...>::value;
    };

    template <Integer Head, Integer... Tail>
    struct NthValue<0, Head, Tail...> {
        static constexpr Integer value = Head;
    };

    /* ------------------------------------------------------------------------------------------- */
    // getting the type of a tuple element using a tuple instead. The same as std::tuple_element
    template <std::size_t I, class T>
    struct tuple_el;

    // recursive case
    template <std::size_t I, class Head, class... Tail>
    struct tuple_el<I, std::tuple<Head, Tail...>> : tuple_el<I - 1, std::tuple<Tail...>> {};

    template <class Head, class... Tail>
    struct tuple_el<0, std::tuple<Head, Tail...>> {
        typedef Head type;
    };

    template <std::size_t I = 0, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type print_tuple(std::tuple<Tp...>& t) {}

    template <std::size_t I = 0, typename... Tp>
        inline typename std::enable_if < I<sizeof...(Tp), void>::type print_tuple(std::tuple<Tp...>& t) {
        std::cout << std::get<I>(t) << std::endl;
        print_tuple<I + 1, Tp...>(t);
    }

    template <std::size_t I = 0, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type reserve_view_tuple(std::tuple<Tp...>& t,
                                                                                      const int size,
                                                                                      const std::string desc) {}

    template <std::size_t I = 0, typename... Tp>
        inline typename std::enable_if <
        I<sizeof...(Tp), void>::type reserve_view_tuple(std::tuple<Tp...>& t, const int size, const std::string desc) {
        std::get<I>(t) = typename std::tuple_element<I, typename std::decay<decltype(t)>::type>::type(
            desc + std::to_string(I), size);
        reserve_view_tuple<I + 1, Tp...>(t, size, desc);
    }

    /************* generic tuple expansion using a functor to apply a template function to each tuple ********/

    // backwards expansion of a tuple from N-0
    template <typename F, size_t Idx, typename... Vs>
    typename std::enable_if<Idx == 0, void>::type apply_each_element_impl(const F& fct, std::tuple<Vs...>& tuple) {
        fct(std::get<0>(tuple));
    }

    template <typename F, size_t Idx, typename... Vs>
    typename std::enable_if<Idx != 0, void>::type apply_each_element_impl(const F& fct, std::tuple<Vs...>& tuple) {
        fct(std::get<Idx>(tuple));
        apply_each_element_impl<F, Idx - 1, Vs...>(fct, tuple);
    }

    template <typename F, typename... Vs>
    void apply_each_element(const F& fct, std::tuple<Vs...>& tuple) {
        apply_each_element_impl<F, sizeof...(Vs) - 1, Vs...>(fct, tuple);
    }

    template <int I, class... Ts>
    auto get_nth_value(Ts&&... ts) -> decltype(std::get<I>(std::forward_as_tuple(ts...))) {
        return std::get<I>(std::forward_as_tuple(ts...));
    }

    template <std::size_t I = 0, std::size_t J, typename F>
    inline typename std::enable_if<I == J, void>::type for_each_tuple_elem(const F& f) {}

    template <std::size_t I = 0, std::size_t J, typename F>
        inline typename std::enable_if < I<J, void>::type for_each_tuple_elem(const F& f) {
        f(I);
        for_each_tuple_elem<I + 1, J, F>(f);
    }

    // forwards expansion of a tuple from 0-N
    template <typename F, std::size_t I = 0, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type apply_impl(const F& f, std::tuple<Tp...>& t) {}

    template <typename F, std::size_t I = 0, typename... Tp>
        inline typename std::enable_if < I<sizeof...(Tp), void>::type apply_impl(const F& f, std::tuple<Tp...>& t) {
        f(std::get<I>(t), I);
        apply_impl<F, I + 1, Tp...>(f, t);
    }

    /* forwards expansion of a tuple from 0-N */
    template <typename F, std::size_t I = 0, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type apply_impl(const F& f,
                                                                              std::tuple<Tp...>& t,
                                                                              std::tuple<Tp...>& v) {}

    template <typename F, std::size_t I = 0, typename... Tp>
        inline typename std::enable_if <
        I<sizeof...(Tp), void>::type apply_impl(const F& f, std::tuple<Tp...>& t, std::tuple<Tp...>& v) {
        f(std::get<I>(t), std::get<I>(v));
        apply_impl<F, I + 1, Tp...>(f, t, v);
    }

    /* forwards expansion of a tuple from 0-N */
    template <typename F, std::size_t I = 0, typename... Tp>
    inline typename std::enable_if<I == sizeof...(Tp), void>::type apply_impl(const F& f,
                                                                              std::tuple<ViewMatrixType<Tp>...>& t,
                                                                              std::tuple<ViewVectorType<Tp>...>& v) {}

    template <typename F, std::size_t I = 0, typename... Tp>
        inline typename std::enable_if < I<sizeof...(Tp), void>::type apply_impl(const F& f,
                                                                                 std::tuple<ViewMatrixType<Tp>...>& t,
                                                                                 std::tuple<ViewVectorType<Tp>...>& v) {
        f(std::get<I>(t), std::get<I>(v));
        apply_impl<F, I + 1, Tp...>(f, t, v);
    }

    template <typename F, Integer I = 0, Integer... Args, typename... Tp>
    typename std::enable_if<I == sizeof...(Args), void>::type for_each_arg(const F& f,
                                                                                                std::tuple<Tp...>& t) {}

    template <typename F, Integer I = 0, Integer... Args, typename... Tp>
        typename std::enable_if <
        I<sizeof...(Args), void>::type for_each_arg(const F& f, std::tuple<Tp...>& t) {
        constexpr Integer dataIndex = NthValue<I, Args...>::value;

        f(std::get<dataIndex>(t), dataIndex);
        for_each_arg<F, I + 1, Args...>(f, t);
    }

    template <typename F, Integer I = 0, Integer... Args, typename... Tp>
    typename std::enable_if<I == sizeof...(Args), void>::type for_each_arg(const F& f,
                                                                                                std::tuple<Tp...>& t,
                                                                                                std::tuple<Tp...>& v) {}

    template <typename F, Integer I = 0, Integer... Args, typename... Tp>
        typename std::enable_if <
        I<sizeof...(Args), void>::type for_each_arg(const F& f,
                                                                         std::tuple<Tp...>& t,
                                                                         std::tuple<Tp...>& v) {
        constexpr Integer dataIndex = NthValue<I, Args...>::value;

        f(std::get<dataIndex>(t), std::get<dataIndex>(v));
        for_each_arg<F, I + 1, Args...>(f, t, v);
    }

    template <typename F, Integer I = 0, Integer... Args, typename... Tp>
    typename std::enable_if<I == sizeof...(Args), void>::type
    for_each_arg(const F& f, std::tuple<ViewMatrixType<Tp>...>& t, std::tuple<ViewVectorType<Tp>...>& v) {}

    template <typename F, Integer I = 0, Integer... Args, typename... Tp>
        typename std::enable_if <
        I<sizeof...(Args), void>::type for_each_arg(const F& f,
                                                                         std::tuple<ViewMatrixType<Tp>...>& t,
                                                                         std::tuple<ViewVectorType<Tp>...>& v) {
        constexpr Integer dataIndex = NthValue<I, Args...>::value;

        f(std::get<dataIndex>(t), std::get<dataIndex>(v));
        for_each_arg<F, I + 1, Args...>(f, t, v);
    }

    template <typename F, typename T, Integer... dataidx>
    static void expand_tuple(const F& f, T& t) {
        if (sizeof...(dataidx) == 0) {
            apply_impl(f, t);
        } else {
            for_each_arg<F, 0, dataidx...>(f, t);
        }
    }

    template <typename F, typename T, Integer... dataidx>
    static void expand_tuple(const F& f, T& t, T& v) {
        if (sizeof...(dataidx) == 0) {
            apply_impl(f, t, v);
        } else {
            for_each_arg<F, 0, dataidx...>(f, t, v);
        }
    }

    template <typename F, typename M, typename T, Integer... dataidx>
    static void expand_tuple(const F& f, M& t, T& v) {
        if (sizeof...(dataidx) == 0) {
            apply_impl(f, t, v);
        } else {
            for_each_arg<F, 0, dataidx...>(f, t, v);
        }
    }

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

}  // namespace mars
#endif
#endif  // MARS_DISTRIBUTED_UTILS_HPP
