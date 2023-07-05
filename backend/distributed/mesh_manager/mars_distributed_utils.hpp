#ifndef MARS_DISTRIBUTED_UTILS_HPP
#define MARS_DISTRIBUTED_UTILS_HPP

#include <ostream>
#include <tuple>
#include <type_traits>
#include "mars_globals.hpp"

#ifdef MARS_ENABLE_KOKKOS

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

    /* ***************************general utils************************************************** */

#ifdef MARS_ENABLE_KOKKOS

    template <typename T, Integer Dim>
    void fill_view_vector(ViewVectorTextureC<T, Dim> v, const T value[Dim]) {
        using namespace Kokkos;

        typename ViewVectorTextureC<T, Dim>::HostMirror h_pt = create_mirror_view(v);

        for (Integer i = 0; i < Dim; ++i) h_pt(i) = value[i];

        deep_copy(v, h_pt);
    }

    template <typename T, Integer XDim, Integer YDim>
    void fill_view_matrix(ViewMatrixTextureC<T, XDim, YDim> m, const T value[XDim][YDim]) {
        using namespace Kokkos;

        typename ViewMatrixTextureC<T, XDim, YDim>::HostMirror h_pt = create_mirror_view(m);

        for (Integer i = 0; i < XDim; ++i) {
            for (Integer j = 0; j < YDim; ++j) {
                h_pt(i, j) = value[i][j];
            }
        }

        deep_copy(m, h_pt);
    }

    template <typename T>
    void print_view(const ViewVectorType<T>& view) {
        using namespace Kokkos;

        parallel_for(
            "print view", view.extent(0), KOKKOS_LAMBDA(const Integer i) {
                std::cout << "i: " << i << " value: " << view(i) << std::endl;
            });
    }

    struct resize_view_functor {
        resize_view_functor(std::string d, size_t s) : _desc(d), _size(s) {}
        template <typename ElementType>
        void operator()(ElementType& el, std::size_t I) const {
            el = ElementType(_desc + std::to_string(I), _size);
        }

        std::string _desc;
        size_t _size;
    };

    struct resize_functor {
        resize_functor(size_t s) : _size(s) {}
        template <typename ElementType>
        void operator()(ElementType& el, std::size_t I) const {
            el = ElementType(_size);
        }
        size_t _size;
    };

    struct print_functor {
        template <typename ElementType>
        void operator()(ElementType& el) const {
            std::cout << el << std::endl;
        }
    };

    /* special implementation of the binary search considering found
    if an element is between current and next proc value. */
    template <typename T>
    MARS_INLINE_FUNCTION Integer
    find_owner_processor(const ViewVectorType<T> view, const T enc_oc, const int offset, Integer guess) {
        const int last_index = view.extent(0) / offset - 1;
        int first_proc = 0;
        int last_proc = last_index;

        while (first_proc <= last_proc && first_proc != last_index) {
            T current = view(offset * guess);
            T next = view(offset * (guess + 1));

            if (enc_oc >= current && enc_oc < next) {
                return guess;
            } else if (enc_oc < current) {
                last_proc = guess - 1;
                guess = (first_proc + last_proc + 1) / 2;
            } else if (enc_oc >= next) {
                first_proc = guess + 1;
                guess = (first_proc + last_proc) / 2;
            }
        }

        return -1;
    }

    // standard binary search
    template <typename T>
    MARS_INLINE_FUNCTION Integer binary_search(const T* view, Integer f, Integer l, const T enc_oc) {
        while (f <= l) {
            Integer guess = (l + f) / 2;

            T current = *(view + guess);

            if (enc_oc == current) {
                return guess;
            } else if (enc_oc < current) {
                l = guess - 1;
            } else {
                f = guess + 1;
            }
        }

        return -1;
    }

//Kokkos way of doing Abs Max Atomics : https://github.com/kokkos/kokkos/pull/5816/files
//Trick the atomic add into doing max abs fetch operator
  template <class Scalar>
  struct AbsMaxHelper {
    Scalar value;

    KOKKOS_FUNCTION AbsMaxHelper& operator+=(AbsMaxHelper const& rhs) {
      Scalar lhs_abs_value = Kokkos::abs(value);
      Scalar rhs_abs_value = Kokkos::abs(rhs.value);
      value = lhs_abs_value > rhs_abs_value ? lhs_abs_value : rhs_abs_value;
      return *this;
    }

    KOKKOS_FUNCTION AbsMaxHelper operator+(AbsMaxHelper const& rhs) const {
      AbsMaxHelper ret = *this;
      ret += rhs;
      return ret;
    }
  };

    // max plus functor
    template <typename T>
    class MaxPlus {
    public:
        // Kokkos reduction functors need the value_type typedef.
        // This is the type of the result of the reduction.
        typedef T value_type;

        MaxPlus(const ViewVectorType<T> x) : x_(x) {}

        // This is helpful for determining the right index type,
        // especially if you expect to need a 64-bit index.
        typedef typename ViewVectorType<T>::size_type size_type;

        KOKKOS_INLINE_FUNCTION void operator()(const size_type i,
                                               value_type& update) const {  // max-plus semiring equivalent of "plus"
            if (update < x_(i)) {
                update = x_(i);
            }
        }

        // "Join" intermediate results from different threads.
        // This should normally implement the same reduction
        // operation as operator() above.
#if (KOKKOS_VERSION >= 40000)
        KOKKOS_INLINE_FUNCTION void join(
            value_type& dst,
            const value_type& src) const {  // max-plus semiring equivalent of "plus"
#else
        // Kokkos forces us to have the input values being declared volatile.
        KOKKOS_INLINE_FUNCTION void join(
            volatile value_type& dst,
            const volatile value_type& src) const {  // max-plus semiring equivalent of "plus"
#endif
            if (dst < src) {
                dst = src;
            }
        }

        // Tell each thread how to initialize its reduction result.
        KOKKOS_INLINE_FUNCTION void init(value_type& dst) const {  // The identity under max is -Inf.
            dst = Kokkos::reduction_identity<value_type>::max();
        }

    private:
        ViewVectorType<T> x_;
    };

    //**********************************dof handler utils******************************************

    // The staggered grid implementation forbids unifying these two methods.
    // They work for both general and staggered dof handlers (unified for handlers).
    // However not possible to unify for owned and local. Either one or the other.
    template <Integer Label, typename H>
    ViewVectorType<Integer> compact_owned_dofs(const H& dof_handler, ViewVectorType<Integer>& locally_owned_dofs) {
        using namespace Kokkos;

        const Integer local_size = dof_handler.template get_owned_dof_size<1>();

        ViewVectorType<bool> dof_predicate("label_dof_predicate", local_size);
        Kokkos::parallel_for(
            "separateowneddoflabels", local_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer local = dof_handler.template get_owned_dof<1>(i);
                if (dof_handler.template get_label<1>(local) & Label) {
                    dof_predicate(i) = 1;
                }
            });

        /* perform a scan on the dof predicate*/
        ViewVectorType<Integer> owned_dof_map("owned_dof_scan", local_size + 1);
        incl_excl_scan(0, local_size, dof_predicate, owned_dof_map);

        auto vol_subview = subview(owned_dof_map, local_size);
        auto h_vs = create_mirror_view(vol_subview);
        // Deep copy device view to host view.
        deep_copy(h_vs, vol_subview);

        locally_owned_dofs = ViewVectorType<Integer>("locally_owned_dofs", h_vs());

        parallel_for(
            local_size, KOKKOS_LAMBDA(const Integer i) {
                if (dof_predicate(i) == 1) {
                    Integer vindex = owned_dof_map(i);
                    const Integer local = dof_handler.template get_owned_dof<1>(i);
                    locally_owned_dofs(vindex) = local;
                }
            });

        return owned_dof_map;
    }

    template <Integer Label, typename H>
    ViewVectorType<Integer> compact_local_dofs(const H& dof_handler, ViewVectorType<Integer>& local_dofs) {
        using namespace Kokkos;

        const Integer local_size = dof_handler.get_base_dof_size();

        ViewVectorType<bool> dof_predicate("label_dof_predicate", local_size);
        Kokkos::parallel_for(
            "separatelocaldoflabels", local_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer local = dof_handler.get_local_dof(i);
                if (dof_handler.template get_label<1>(local) & Label) {
                    dof_predicate(i) = 1;
                }
            });

        /* perform a scan on the dof predicate*/
        ViewVectorType<Integer> local_dof_map("local_dof_scan", local_size + 1);
        incl_excl_scan(0, local_size, dof_predicate, local_dof_map);

        auto vol_subview = subview(local_dof_map, local_size);
        auto h_vs = create_mirror_view(vol_subview);
        // Deep copy device view to host view.
        deep_copy(h_vs, vol_subview);

        local_dofs = ViewVectorType<Integer>("local_dofs", h_vs());

        /* Compact the predicate into the volume and face dofs views */
        parallel_for(
            local_size, KOKKOS_LAMBDA(const Integer i) {
                if (dof_predicate(i) == 1) {
                    Integer vindex = local_dof_map(i);
                    const Integer local = dof_handler.get_local_dof(i);
                    local_dofs(vindex) = local;
                }
            });

        return local_dof_map;
    }

    template <typename DH, typename F>
    ViewVectorType<bool> build_sfc_to_local_predicate(DH dofhandler,
                                                      F f,
                                                      const Integer local_size,
                                                      const ViewVectorType<Integer> in) {
        ViewVectorType<bool> predicate("separated_predicate", local_size);

        Kokkos::parallel_for(
            "separatedoflabelss", local_size, KOKKOS_LAMBDA(const Integer i) {
                const Integer sfc = in(i);
                /* if (f(dofhandler.sfc_to_local(sfc))) { */
                if (f(sfc)) {
                    predicate(i) = 1;
                }
            });

        return predicate;
    }

    template <typename DH, typename F>
    inline ViewVectorType<bool> compact_sfc_to_local(DH dofhandler,
                                                     F f,
                                                     const ViewVectorType<Integer> in,
                                                     ViewVectorType<Integer>& out) {
        using namespace Kokkos;

        auto local_size = in.extent(0);
        auto in_predicate = build_sfc_to_local_predicate<DH, F>(dofhandler, f, local_size, in);

        /* perform a scan on the separated dof predicate*/
        auto bscan = ViewVectorType<Integer>("boudnary_separated_scan", local_size + 1);
        incl_excl_scan(0, local_size, in_predicate, bscan);

        auto vol_subview = subview(bscan, local_size);
        auto h_vs = create_mirror_view(vol_subview);
        // Deep copy device view to host view.
        deep_copy(h_vs, vol_subview);

        out = ViewVectorType<Integer>("local_separated_dofs", h_vs());
        /* ViewVectorType<Integer> bvds = out; */

        /* Compact the predicate into the separated and face dofs views */
        parallel_for(
            local_size, KOKKOS_LAMBDA(const Integer i) {
                if (in_predicate(i) == 1) {
                    Integer vindex = bscan(i);
                    out(vindex) = dofhandler.sfc_to_local(in(i));
                }
            });
        return in_predicate;
    }

    inline ViewVectorType<Integer> count_sfc_to_local(const ViewVectorType<Integer> in_scan,
                                                      const ViewVectorType<bool> in_predicate) {
        using namespace Kokkos;

        const Integer count_size = in_scan.extent(0) - 1;

        ViewVectorType<Integer> count("count_to_local", count_size);
        ViewVectorType<Integer> scan("scan_to_local", count_size + 1);
        /* Compact the predicate into the separated and face dofs views */
        parallel_for(
            in_predicate.extent(0), KOKKOS_LAMBDA(const Integer i) {
                if (in_predicate(i) == 1) {
                    Integer proc = find_owner_processor(in_scan, i, 1, 0);
                    atomic_increment(&count(proc));
                }
            });

        incl_excl_scan(0, count_size, count, scan);

        return scan;
    }
#endif
}  // namespace mars
#endif
#endif  // MARS_DISTRIBUTED_UTILS_HPP
