#include <tuple>
#include <type_traits>

namespace mars
{

/******************** Tuple utils *****************************/

/* getting the index of a type in a variadic template definition */

template <class Tuple, class T, std::size_t Index = 0>
struct TypeIdx;

template <std::size_t Index, bool Valid>
struct TypeIdxTest
    : public std::integral_constant<std::size_t, Index>
{
};

template <std::size_t Index>
struct TypeIdxTest<Index, false>
{
    static_assert(Index == -1, "Type not found in the tuple!");
};

template <class Head, class T, std::size_t Index>
struct TypeIdx<std::tuple<Head>, T, Index>
    : public TypeIdxTest<Index, std::is_same<Head, T>::value>
{
};

template <class Head, class... Rest, class T, std::size_t Index>
struct TypeIdx<std::tuple<Head, Rest...>, T, Index>
    : public std::conditional<std::is_same<Head, T>::value,
                              std::integral_constant<std::size_t, Index>,
                              TypeIdx<std::tuple<Rest...>, T, Index + 1>>::type
{
};

/* ------------------------------------------------------------------------------------------------ */
//getting the  nth type of a variadic.
template <std::size_t I, typename... T>
struct NthType;

// recursive case
template <std::size_t I, typename Head, typename... Tail>
struct NthType<I, Head, Tail...>
{
    typedef typename NthType<I-1, Tail...>::type type;
};

template <class Head, class... Tail>
struct NthType<0, Head, Tail...>
{
    typedef Head type;
};

/* ------------------------------------------------------------------------------------------- */
//getting the type of a tuple element using a tuple instead. The same as std::tuple_element
template <std::size_t I, class T>
struct tuple_el;

// recursive case
template <std::size_t I, class Head, class... Tail>
struct tuple_el<I, std::tuple<Head, Tail...>>
    : tuple_el<I - 1, std::tuple<Tail...>>
{
};

template <class Head, class... Tail>
struct tuple_el<0, std::tuple<Head, Tail...>>
{
    typedef Head type;
};

template <std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
print_tuple(std::tuple<Tp...> &t) {}

template <std::size_t I = 0, typename... Tp>
inline typename std::enable_if <
I<sizeof...(Tp), void>::type print_tuple(std::tuple<Tp...> &t)
{
    std::cout << std::get<I>(t) << std::endl;
    print_tuple<I + 1, Tp...>(t);
}

template <std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
reserve_view_tuple(std::tuple<Tp...> &t, const int size, const std::string desc) {}

template <std::size_t I = 0, typename... Tp>
 inline typename std::enable_if <
I<sizeof...(Tp), void>::type reserve_view_tuple(std::tuple<Tp...> &t, const int size, const std::string desc)
{
    std::get<I>(t) = typename std::tuple_element<I, typename std::decay<decltype(t)>::type>::type(desc + std::to_string(I), size);
    reserve_view_tuple<I + 1, Tp...>(t, size, desc);
}



/************* generic tuple expansion using a functor to apply a template function to each tuple ********/

//backwards expansion of a tuple from N-0
template <typename F, size_t Idx, typename... Vs>
typename std::enable_if<Idx == 0, void>::type
apply_each_element_impl(const F &fct, std::tuple<Vs...> &tuple)
{
    fct(std::get<0>(tuple));
}

template <typename F, size_t Idx, typename... Vs>
typename std::enable_if<Idx != 0, void>::type
apply_each_element_impl(const F &fct, std::tuple<Vs...> &tuple)
{
    fct(std::get<Idx>(tuple));
    apply_each_element_impl<F, Idx - 1, Vs...>(fct, tuple);
}

template <typename F, typename... Vs>
void apply_each_element(const F &fct, std::tuple<Vs...> &tuple)
{
    apply_each_element_impl<F, sizeof...(Vs) - 1, Vs...>(fct, tuple);
}

template <int I, class... Ts>
auto get_nth_value(Ts&&... ts) -> decltype(std::get<I>(std::forward_as_tuple(ts...)))
{
  return std::get<I>(std::forward_as_tuple(ts...));
}


template <std::size_t I=0, std::size_t J, typename F>
inline typename std::enable_if<I == J, void>::type
for_each_tuple_elem(const F &f) {}

template <std::size_t I=0, std::size_t J, typename F>
inline typename std::enable_if <
I<J, void>::type for_each_tuple_elem(const F &f)
{
    f(I);
    for_each_tuple_elem<I + 1, J, F>(f);
}


//forwards expansion of a tuple from 0-N
template <typename F, std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
apply_impl(const F &f, std::tuple<Tp...> &t) {}

template <typename F, std::size_t I = 0, typename... Tp>
inline typename std::enable_if <
I<sizeof...(Tp), void>::type apply_impl(const F &f, std::tuple<Tp...> &t)
{
    f(std::get<I>(t), I);
    apply_impl<F, I + 1, Tp...>(f, t);
}


/* forwards expansion of a tuple from 0-N */
template <typename F, std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
apply_impl(const F &f, std::tuple<Tp...> &t, std::tuple<Tp...> &v) {}

template <typename F, std::size_t I = 0, typename... Tp>
inline typename std::enable_if <
I<sizeof...(Tp), void>::type apply_impl(const F &f, std::tuple<Tp...> &t, std::tuple<Tp...> &v)
{
    f(std::get<I>(t), std::get<I>(v));
    apply_impl<F, I + 1, Tp...>(f, t, v);
}

/* *************************************************************************************** */

/* other utils */

struct resize_view_functor
{
    resize_view_functor(std::string d, size_t s) : _desc(d), _size(s)  {}
    template <typename ElementType>
    void operator()(ElementType &el, std::size_t I) const
    {
        el = ElementType(_desc + std::to_string(I), _size);
    }

    std::string _desc;
    size_t _size;
};

struct resize_functor
{
    resize_functor(size_t s) : _size(s) {}
    template <typename ElementType>
    void operator()(ElementType &el, std::size_t I) const
    {
        el = ElementType(_size);
    }
    size_t _size;
};

struct print_functor
{
    template <typename ElementType>
    void operator()(ElementType &el) const
    {
        std::cout << el << std::endl;
    }
};

//binary search on some view.
template <typename T>
MARS_INLINE_FUNCTION Integer find_owner_processor(const ViewVectorType<T> view,
                                                  const T enc_oc, const int offset, Integer guess)
{
    const int last_index = view.extent(0) / offset - 1;
    int first_proc = 0;
    int last_proc = last_index;

    //special implementation of the binary search considering found if an element is between current and next proc value.
    while (first_proc <= last_proc && first_proc != last_index)
    {
        T current = view(offset * guess);
        T next = view(offset * (guess + 1));

        if (enc_oc >= current && enc_oc < next)
        {
            return guess;
        }
        else if (enc_oc < current)
        {
            last_proc = guess - 1;
            guess = (first_proc + last_proc + 1) / 2;
        }
        else if (enc_oc >= next)
        {
            first_proc = guess + 1;
            guess = (first_proc + last_proc) / 2;
        }
    }

    return -1;
}

} //namespace mars
