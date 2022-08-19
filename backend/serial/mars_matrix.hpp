#ifndef MARS_MATRIX_HPP
#define MARS_MATRIX_HPP

#include "mars_base.hpp"
#include "mars_vector.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <iostream>

namespace mars {

    template <typename T, Integer Rows, Integer Cols>
    class Matrix {
    public:
        Matrix() {}

        Matrix(std::initializer_list<T> values) {
            assert(values.size() == Rows * Cols);
            std::copy(std::begin(values), std::begin(values) + (Rows * Cols), std::begin(this->values));
        }

        inline constexpr static Integer rows() { return Rows; }
        inline constexpr static Integer cols() { return Cols; }

        inline T &operator()(const Integer i, const Integer j) {
            assert(i < Rows);
            assert(j < Cols);
            return values[i * cols() + j];
        }

        inline const T &operator()(const Integer i, const Integer j) const {
            assert(i < Rows);
            assert(j < Cols);
            return values[i * cols() + j];
        }

        inline void col(const Integer c, const Vector<T, Rows> &v) {
            assert(c < Cols);

            for (Integer d = 0; d < Rows; ++d) {
                (*this)(d, c) = v(d);
            }
        }

        inline void zero() { std::fill(begin(values), end(values), 0.); }

        void describe(std::ostream &os) const {
            for (Integer i = 0; i < Rows; ++i) {
                for (Integer j = 0; j < Cols; ++j) {
                    os << (*this)(i, j) << " ";
                }
                os << "\n";
            }

            os << "\n";
        }

        friend std::ostream &operator<<(std::ostream &os, const Matrix &m) {
            m.describe(os);
            return os;
        }

        template <Integer OtherCols>
        inline Matrix<T, Rows, OtherCols> operator*(const Matrix<T, Cols, OtherCols> &other) const {
            Matrix<T, Rows, OtherCols> ret;
            ret.zero();

            for (Integer i = 0; i < Rows; ++i) {
                for (Integer j = 0; j < Cols; ++j) {
                    for (Integer k = 0; k < OtherCols; ++k) {
                        ret(i, k) += (*this)(i, j) * other(j, k);
                    }
                }
            }

            return ret;
        }

        std::array<T, Rows * Cols> values;
    };

    template <typename T, Integer N>
    inline void m_minor(const Integer cof_i,
                        const Integer cof_j,
                        const Matrix<T, N, N> &mat,
                        Matrix<T, N - 1, N - 1> &m) {
        Integer i_offset = 0;
        for (Integer i = 0; i < N; ++i) {
            if (i == cof_i) {
                i_offset = -1;
                continue;
            }

            Integer j_offset = 0;

            for (Integer j = 0; j < N; ++j) {
                if (j == cof_j) {
                    j_offset = -1;
                    continue;
                }

                m(i + i_offset, j + j_offset) = mat(i, j);
            }
        }
    }

    template <typename T, Integer N>
    inline T det_aux(const Matrix<T, N, N> &m) {
        static_assert(N < 7, "max size is 6");

        std::array<Integer, N> nnz;
        std::fill(std::begin(nnz), std::end(nnz), 0);

        for (Integer i = 0; i < N; ++i) {
            for (Integer j = 0; j < N; ++j) {
                nnz[i] += m(i, j) != 0.;
            }
        }

        Integer row = 0;
        for (Integer i = 0; i < N; ++i) {
            if (nnz[row] > nnz[i]) {
                row = i;
            }
        }

        assert(row < N);

        if (nnz[row] == 0) return 0.;

        Integer ret = 0.;
        Matrix<T, N - 1, N - 1> mij;
        for (Integer j = 0; j < N; ++j) {
            const Real coff = m(row, j);
            if (coff == 0.) continue;

            const Real sign = (j % 2) == 0 ? 1 : -1;
            m_minor(row, j, m, mij);
            ret += sign * coff * det(mij);
        }

        return ret;
    }

    template <typename T, Integer N>
    inline T det(const Matrix<T, N, N> &m) {
        return det_aux(m);
    }

    template <typename T, Integer Rows, Integer Cols>
    inline T det(const Matrix<T, Rows, Cols> &m) {
        if (Rows > Cols) {
            Matrix<T, Rows, Rows> m_square;
            for (Integer i = 0; i < Rows; ++i) {
                for (Integer j = 0; j < Cols; ++j) {
                    m_square(i, j) = m(i, j);
                }

                for (Integer j = Cols; j < Rows; ++j) {
                    m_square(i, j) = i == j;
                }
            }

            return det(m_square);
        } else {
            Matrix<T, Cols, Cols> m_square;

            for (Integer j = 0; j < Cols; ++j) {
                for (Integer i = 0; i < Rows; ++i) {
                    m_square(i, j) = m(i, j);
                }

                for (Integer i = Rows; i < Cols; ++i) {
                    m_square(i, j) = i == j;
                }
            }

            return det(m_square);
        }
    }

    template <typename T, Integer Rows, Integer Cols>
    Matrix<T, Cols, Rows> transpose(const Matrix<T, Rows, Cols> &mat) {
        Matrix<T, Cols, Rows> ret;
        for (Integer i = 0; i < Rows; ++i) {
            for (Integer j = 0; j < Cols; ++j) {
                ret(j, i) = mat(i, j);
            }
        }

        return ret;
    }

    template <typename T>
    inline T det(const Matrix<T, 1, 1> &m) {
        return m(0, 0);
    }

    template <typename T>
    inline T det(const Matrix<T, 2, 2> &m) {
        return m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1);
    }

    template <typename T>
    inline T det(const Matrix<T, 3, 3> &m) {
        return m(0, 0) * m(1, 1) * m(2, 2) + m(0, 1) * m(1, 2) * m(2, 0) + m(0, 2) * m(1, 0) * m(2, 1) -
               m(0, 0) * m(1, 2) * m(2, 1) - m(0, 1) * m(1, 0) * m(2, 2) - m(0, 2) * m(1, 1) * m(2, 0);
    }

    template <typename T>
    inline T det(const Matrix<T, 4, 4> &m) {
        const auto m00 = m(0, 0);
        const auto m01 = m(0, 1);
        const auto m02 = m(0, 2);
        const auto m03 = m(0, 3);

        const auto m10 = m(1, 0);
        const auto m11 = m(1, 1);
        const auto m12 = m(1, 2);
        const auto m13 = m(1, 3);

        const auto m20 = m(2, 0);
        const auto m21 = m(2, 1);
        const auto m22 = m(2, 2);
        const auto m23 = m(2, 3);

        const auto m30 = m(3, 0);
        const auto m31 = m(3, 1);
        const auto m32 = m(3, 2);
        const auto m33 = m(3, 3);

        return (m00 == 0. ? 0.
                          : (m00 * det(Matrix<T, 3, 3>({m11, m12, m13, m21, m22, m23, m31, m32, m33

                                   })))) -
               (m01 == 0. ? 0. : m01 * det(Matrix<T, 3, 3>({m10, m12, m13, m20, m22, m23, m30, m32, m33}))) +
               (m02 == 0. ? 0.
                          : (m02 * det(Matrix<T, 3, 3>({m10, m11, m13, m20, m21, m23, m30, m31, m33

                                   })))) -
               (m03 == 0. ? 0.
                          : (m03 * det(Matrix<T, 3, 3>({m10, m11, m12, m20, m21, m22, m30, m31, m32

                                   }))));
    }
}  // namespace mars

#endif  // MARS_MATRIX_HPP
