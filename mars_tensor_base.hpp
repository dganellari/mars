#ifndef MARS_TENSOR_BASE_HPP
#define MARS_TENSOR_BASE_HPP 
#include <array>
#include <initializer_list>
#include <cmath>
#include <iostream>
#include "mars_static_math.hpp"

namespace mars {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// TensorBase is used to define also static constexpr matrices and vectors, without compiler problems
    //// To have more informations regarding the eventual errors, visit:
    //// https://stackoverflow.com/questions/57131334/how-to-initialise-a-constexpr-matrix-of-type-matrixt-rows-cols
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template <typename, typename>
    class TensorBase;

    template <typename T, Integer ... Is>
    class TensorBase<T, std::index_sequence<Is...>>
    {
    protected:
        std::array<T, sizeof...(Is)> values{};

    public:
        constexpr TensorBase (getTypeSequence<T, Is> ... vals)
        : values{{vals...}}
        {}

        constexpr TensorBase (std::array<T, sizeof...(Is)> const & a)
        : values{a}
        {}

        constexpr TensorBase (std::array<T, sizeof...(Is)> && a)
        : values{std::move(a)}
        {}
 
        constexpr TensorBase () = default;

        ~TensorBase() = default;

        constexpr TensorBase (TensorBase const &) = default;
        constexpr TensorBase (TensorBase &&) = default;

        constexpr TensorBase & operator= (TensorBase const &) = default;
        constexpr TensorBase & operator= (TensorBase &&) = default;
    };

}


#endif
