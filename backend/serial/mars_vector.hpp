#ifndef MARS_VECTOR_HPP
#define MARS_VECTOR_HPP

#include <array>
#include <cmath>
#include <initializer_list>
#include <iostream>

namespace mars {

    template <typename T, Integer Dim>
    class Vector {
    public:
        std::array<T, Dim> values;
        Vector() {}

        Vector(std::initializer_list<T> values) {
            assert(values.size() == Dim);
            std::copy(std::begin(values), std::begin(values) + Dim, std::begin(this->values));
        }

        friend Vector operator*(const T &factor, const Vector &v) {
            Vector ret;
            for (Integer i = 0; i < Dim; ++i) {
                ret(i) = factor * v(i);
            }

            return ret;
        }

        friend Vector operator/(const Vector &v, const T &factor) {
            Vector ret;
            for (Integer i = 0; i < Dim; ++i) {
                ret(i) = v(i) / factor;
            }

            return ret;
        }

        Vector operator-(const Vector &right) const {
            Vector ret;
            for (Integer i = 0; i < Dim; ++i) {
                ret(i) = (*this)(i)-right(i);
            }

            return ret;
        }

        bool operator<(const Vector &right) const {
            for (Integer i = 0; i < Dim; ++i) {
                if ((*this)(i) < right(i)) return true;
                if ((*this)(i) > right(i)) return false;
            }

            return false;
        }

        Vector operator+(const Vector &right) const {
            Vector ret;
            for (Integer i = 0; i < Dim; ++i) {
                ret(i) = (*this)(i) + right(i);
            }

            return ret;
        }

        Vector operator+=(const Vector &right) {
            for (Integer i = 0; i < Dim; ++i) {
                (*this)(i) += right(i);
            }

            return *this;
        }

        Vector operator*=(const T &right) {
            for (Integer i = 0; i < Dim; ++i) {
                (*this)(i) *= right;
            }

            return *this;
        }

        Vector operator/=(const T &right) {
            for (Integer i = 0; i < Dim; ++i) {
                (*this)(i) /= right;
            }

            return *this;
        }

        Vector operator*(const Vector &right) const {
            Vector ret;
            for (Integer i = 0; i < Dim; ++i) {
                ret(i) = (*this)(i)*right(i);
            }

            return ret;
        }

        inline T &operator()(const Integer i) {
            assert(i < Dim);
            return values[i];
        }

        inline const T &operator()(const Integer i) const {
            assert(i < Dim);
            return values[i];
        }

        inline T &operator[](const Integer i) {
            assert(i < Dim);
            return values[i];
        }

        inline const T &operator[](const Integer i) const {
            assert(i < Dim);
            return values[i];
        }

        void describe(std::ostream &os) const {
            for (Integer i = 0; i < Dim; ++i) {
                os << (*this)(i) << " ";
            }

            os << "\n";
        }

        friend std::ostream &operator<<(std::ostream &os, const Vector &v) {
            v.describe(os);
            return os;
        }

        inline T squared_norm() const {
            T sqn = (*this)(0) * (*this)(0);
            for (Integer i = 1; i < Dim; ++i) {
                sqn += (*this)(i) * (*this)(i);
            }

            return sqn;
        }

        inline T squared_distance(const Vector &right) const { return ((*this) - right).squared_norm(); }

        inline T norm() const { return std::sqrt(squared_norm()); }

        Vector &normalize() {
            T len = norm();

            for (Integer i = 0; i < Dim; ++i) {
                (*this)(i) /= len;
            }

            return *this;
        }

        inline Vector &zero() {
            std::fill(begin(values), end(values), 0.);
            return *this;
        }

        inline Vector &set(const Real value) {
            std::fill(begin(values), end(values), value);
            return *this;
        }
    };

    template <typename T, Integer Dim>
    inline T dot(const Vector<T, Dim> &left, const Vector<T, Dim> &right) {
        T ret = 0.;
        for (Integer d = 0; d < Dim; ++d) {
            ret += left(d) * right(d);
        }

        return ret;
    }
}  // namespace mars

#endif  // MARS_VECTOR_HPP
