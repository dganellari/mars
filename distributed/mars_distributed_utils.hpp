#include "tuple"

namespace mars
{

/* Tuple utils */

//getting the type of a tuple element. The same as std::tuple_element
template <std::size_t I, class T>
struct tuple_el;

// recursive case
template <std::size_t I, class Head, class... Tail>
struct tuple_el<I, std::tuple<Head, Tail...>>
    : tuple_el<I - 1, std::tuple<Tail...>>
{
};

// base case
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

/* generic tuple expansion using a functor to apply a template function to each tuple */

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


//forwards expansion of a tuple from 0-N
template <typename F, std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
apply_impl(const F &f, std::tuple<Tp...> &t) {}

template <typename F, std::size_t I = 0, typename... Tp>
inline typename std::enable_if <
I<sizeof...(Tp), void>::type apply_impl(const F &f, std::tuple<Tp...> &t)
{
    f(std::get<I>(t));
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
    apply_impl<F, I + 1, Tp...>(f, t);
}

/* other utils */

struct resize_view_functor
{
    resize_view_functor(std::string d, size_t s) : _desc(d), _size(s)  {}
    template <typename ElementType>
    void operator()(ElementType &el) const
    {
        el = ElementType(_desc, _size);
    }

    std::string _desc;
    size_t _size;
};

struct resize_functor
{
    resize_functor(size_t s) : _size(s) {}
    template <typename ElementType>
    void operator()(ElementType &el) const
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

} //namespace mars
