#ifndef MARS_SubView_HPP
#define MARS_SubView_HPP

#include <array>
#include <cmath>
#include <initializer_list>
#include <iostream>

#include "mars_base.hpp"

#include "mars_utils_kokkos.hpp"

namespace mars {

    template <typename T, Integer Dim>
    class SubView {
    public:
        // Keep the allocation alive and make the handle safe to capture in a device kernel.
        ViewMatrixType<T> view;
        Integer index;
        bool valid;

        MARS_INLINE_FUNCTION SubView() : index(0), valid(true) {}

        MARS_INLINE_FUNCTION SubView(const ViewMatrixType<T> *v, Integer id) : view(*v), index(id), valid(true) {}

        MARS_INLINE_FUNCTION SubView(const ViewMatrixType<T> &v, Integer id) : view(v), index(id), valid(true) {}

        MARS_INLINE_FUNCTION void set_valid(bool v = true) { valid = v; }

        MARS_INLINE_FUNCTION bool is_valid() const { return valid; }

        MARS_INLINE_FUNCTION SubView operator+=(const SubView &right) {
            for (Integer i = 0; i < Dim; ++i) {
                view(index, i) += right.view(right.index, i);
            }

            return *this;
        }

        MARS_INLINE_FUNCTION SubView operator*=(const T &right) {
            for (Integer i = 0; i < Dim; ++i) {
                view(index, i) *= right;
            }

            return *this;
        }

        MARS_INLINE_FUNCTION SubView operator/=(const T &right) {
            for (Integer i = 0; i < Dim; ++i) {
                view(index, i) /= right;
            }

            return *this;
        }

        MARS_INLINE_FUNCTION T &operator()(const Integer i) {
            assert(i < Dim);
            return view(index, i);
        }

        MARS_INLINE_FUNCTION const T &operator()(const Integer i) const {
            assert(i < Dim);
            return view(index, i);
        }

        MARS_INLINE_FUNCTION T &operator[](const Integer i) {
            assert(i < Dim);
            return view(index, i);
        }

        MARS_INLINE_FUNCTION const T &operator[](const Integer i) const {
            assert(i < Dim);
            return view(index, i);
        }

        void describe(std::ostream &os) const {
            for (Integer i = 0; i < Dim; ++i) {
                os << view(index, i) << " ";
            }

            os << "\n";
        }

        friend std::ostream &operator<<(std::ostream &os, const SubView &v) {
            v.describe(os);
            return os;
        }

        MARS_INLINE_FUNCTION T squared_norm() const {
            T sqn = view(index, 0) * view(index, 0);
            for (Integer i = 1; i < Dim; ++i) {
                sqn += view(index, i) * view(index, i);
            }

            return sqn;
        }

        MARS_INLINE_FUNCTION T norm() const { return std::sqrt(squared_norm()); }

        MARS_INLINE_FUNCTION SubView &normalize() {
            T len = norm();

            for (Integer i = 0; i < Dim; ++i) {
                view(index, i) /= len;
            }

            return *this;
        }

        MARS_INLINE_FUNCTION SubView &zero() {
            for (Integer i = 0; i < Dim; ++i) {
                view(index, i) = 0;
            }

            return *this;
        }

        MARS_INLINE_FUNCTION SubView &set(const T value) {
            for (Integer i = 0; i < Dim; ++i) {
                view(index, i) = value;
            }

            return *this;
        }

        MARS_INLINE_FUNCTION Integer size() { return Dim; }
    };

    template <typename T, Integer Dim>
    MARS_INLINE_FUNCTION T dot(const SubView<T, Dim> &left, const SubView<T, Dim> &right) {
        T ret = 0.;
        for (Integer d = 0; d < Dim; ++d) {
            ret += left.view(left.index, d) * right.view(right.index, d);
        }

        return ret;
    }
}  // namespace mars

#endif  // MARS_SUBVIEW_HPP
