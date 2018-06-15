#ifndef MARS_SIMPLEX_HPP
#define MARS_SIMPLEX_HPP

#include <array>
#include <vector>
#include <ostream>
#include <iostream>

#include <cmath>
#include <cassert>
#include <initializer_list>
#include <algorithm>

/**
 * M.A.R.S Mesh Adaptive Refinement for Simplical meshes
 */
namespace mars {
    using Real    = double;
    using Integer = int;
    static const int INVALID_INDEX = -1;

    template<Integer N>
    class Factorial {
    public:
        static const Integer value = Factorial<N-1>::value * N;
    };

    template<>
    class Factorial<1> {
    public:
        static const Integer value = 1;
    };

    template<>
    class Factorial<0> {
    public:
        static const Integer value = 0;
    };
    
    template<typename T, Integer Dim>
    class Vector {
    public:
        std::array<T, Dim> values;
        Vector() {}

        Vector(std::initializer_list<T> values)
        {
            std::copy(std::begin(values), std::end(values), std::begin(this->values));
        }
        
        Vector operator-(const Vector &right) const
        {
            Vector ret;
            for(Integer i = 0; i < Dim; ++i) {
                ret(i) = (*this)(i) - right(i);
            }
            
            return ret;
        }
        
        Vector operator+(const Vector &right) const
        {
            Vector ret;
            for(Integer i = 0; i < Dim; ++i) {
                ret(i) = (*this)(i) + right(i);
            }
            
            return ret;
        }

        Vector operator*(const Vector &right) const
        {
            Vector ret;
            for(Integer i = 0; i < Dim; ++i) {
                ret(i) = (*this)(i) * right(i);
            }
            
            return ret;
        }
        
        inline T &operator()(const Integer i)
        {
            assert(i < Dim);
            return values[i];
        }
        
        inline const T &operator()(const Integer i) const
        {
            assert(i < Dim);
            return values[i];
        }
        
        inline T &operator[](const Integer i)
        {
            assert(i < Dim);
            return values[i];
        }
        
        inline const T &operator[](const Integer i) const
        {
            assert(i < Dim);
            return values[i];
        }


        void describe(std::ostream &os) const
        {
            for(Integer i = 0; i < Dim; ++i) {
                os << (*this)(i) << " ";
            }

            os << "\n";
        }

        friend std::ostream &operator<<(std::ostream &os, const Vector &v)
        {
            v.describe(os);
            return os;
        }

        inline T norm() const
        {
            T sqn = (*this)(0) * (*this)(0);
            for(Integer i = 1; i < Dim; ++i)
            {
                sqn += (*this)(i) * (*this)(i);
            }

            return std::sqrt(sqn);
        }

        Vector &normalize() {
            T len = norm();

            for(Integer i = 0; i < Dim; ++i) {
                (*this)(i) /= len;
            }

            return *this;
        }

        inline void zero()
        {
            std::fill(begin(values), end(values), 0.);
        }
    };
    
    template<typename T, Integer Rows, Integer Cols>
    class Matrix {
    public:
    	Matrix() {}

        Matrix(std::initializer_list<T> values)
        {
            std::copy(std::begin(values), std::end(values), std::begin(this->values));
        }
        
        inline constexpr static Integer rows() { return Rows; }
        inline constexpr static Integer cols() { return Cols; }
        
        inline T &operator()(const Integer i, const Integer j)
        {
            assert(i < Rows);
            assert(j < Cols);
            return values[i*cols() + j];
        }
        
        inline const T &operator()(const Integer i, const Integer j) const
        {
            assert(i < Rows);
            assert(j < Cols);
            return values[i*cols() + j];
        }
        
        inline void col(const Integer c, const Vector<T, Rows> &v)
        {
            assert(c < Cols);

            for(Integer d = 0; d < Rows; ++d) {
                (*this)(d, c) = v(d);
            }
        }
            
        inline void zero()
        {
            std::fill(begin(values), end(values), 0.);
        }


        void describe(std::ostream &os) const
        {
            for(Integer i = 0; i < Rows; ++i) {
                for(Integer j = 0; j < Cols; ++j) {
                    os << (*this)(i, j) << " ";
                }
                os << "\n";
            }

            os << "\n";
        }

        friend std::ostream &operator<<(std::ostream &os, const Matrix &m)
        {
            m.describe(os);
            return os;
        }

        std::array<T, Rows * Cols> values;
    };

    template<typename T, Integer Rows, Integer Cols>
    inline T det(const Matrix<T, Rows, Cols> &m)
    {
        if(Rows > Cols) {
            Matrix<T, Rows, Rows> m_square;
            for(Integer i = 0; i < Rows; ++i) {

                for(Integer j = 0; j < Cols; ++j) {
                    m_square(i, j) = m(i, j);
                }

                for(Integer j = Cols; j < Rows; ++j) {
                    m_square(i, j) = i == j;
                }
            }

            return det(m_square);
        } else {
            Matrix<T, Cols, Cols> m_square;

            for(Integer j = 0; j < Cols; ++j) {
                for(Integer i = 0; i < Rows; ++i) {
                    m_square(i, j) = m(i, j);
                }

                for(Integer i = Rows; i < Cols; ++i) {
                    m_square(i, j) = i == j;
                }
            }

            return det(m_square);
        }

        
    }
    
    template<typename T>
    inline T det(const Matrix<T, 1, 1> &m)
    {
        return m(0, 0);
    }
    
    template<typename T>
    inline T det(const Matrix<T, 2, 2> &m)
    {
        return m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1);
    }
    
    template<typename T>
    inline T det(const Matrix<T, 3, 3> &m)
    {
        return
        m(0,0) * m(1,1) * m(2,2)  +
        m(0,1) * m(1,2) * m(2,0)  +
        m(0,2) * m(1,0) * m(2,1)  -
        m(0,0) * m(1,2) * m(2,1)  -
        m(0,1) * m(1,0) * m(2,2)  -
        m(0,2) * m(1,1) * m(2,0);
    }
    
    template<typename T>
    inline T det(const Matrix<T, 4, 4> &m)
    {
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
        
        return (m00 == 0. ? 0. : (
        m00 * det(
            Matrix<T, 3, 3>({
                m11, m12, m13,
                m21, m22, m23,
                m31, m32, m33
            
        }))))
        -
        (m01 == 0. ? 0. :
         m01 * det(
            Matrix<T, 3, 3>({
                m10, m12, m13,
                m20, m22, m23,
                m30, m32, m33
        })))
        +
        (m02 == 0. ? 0. : (
         m02 * det(
            Matrix<T, 3, 3>({
                m10, m11, m13,
                m20, m21, m23,
                m30, m31, m33
            
        }))))
        -
        (m03 == 0. ? 0. : (
         m03 * det(
            Matrix<T, 3, 3>({
                m10, m11, m12,
                m20, m21, m22,
                m30, m31, m32
            
        }))));
    }
    
    template<Integer Dim, Integer ManifoldDim>
    class Simplex {};
    
    template<Integer Dim>
    using Node        = Simplex<Dim, 0>;
    
    template<Integer Dim>
    using Line        = Simplex<Dim, 1>;
    
    template<Integer Dim>
    using Triangle    = Simplex<Dim, 2>;
    
    template<Integer Dim>
    using Tetrahedron = Simplex<Dim, 3>;
    
    template<Integer Dim>
    using Pentatope   = Simplex<Dim, 4>;

    using Node1        = Node<1>;
    using Line1        = Line<1>;
    using Vector1r     = Vector<Real, 1>;

    using Node2        = Node<2>;
    using Line2        = Line<2>;
    using Triangle2    = Triangle<2>;
    using Vector2r     = Vector<Real, 2>;

    using Node3        = Node<3>;
    using Line3        = Line<3>;
    using Triangle3    = Triangle<3>;
    using Tetrahedron3 = Tetrahedron<3>;
    using Vector3r     = Vector<Real, 3>;
    
    using Node4 	   = Node<4>;
    using Line4        = Line<4>;
    using Triangle4    = Triangle<4>;
    using Tetrahedron4 = Tetrahedron<4>;
    using Pentatope4   = Pentatope<4>;
    using Vector4r     = Vector<Real, 4>;

    
    template<Integer Dim>
    class Simplex<Dim, 0> {
    public:
        Integer id;
    };
    
    template<Integer Dim>
    class Simplex<Dim, 1> {
    public:
        std::array<Integer, 2> nodes;
        Integer id;
    };
    
    template<Integer Dim>
    class Simplex<Dim, 2> {
    public:
        std::array<Integer, 3> sides;
        std::array<Integer, 3> nodes;
        Integer id;
    };
    
    template<Integer Dim>
    class Simplex<Dim, 3> {
    public:
        std::array<Integer, 4> sides;
        std::array<Integer, 4> nodes;
        Integer id;
    };
    
    template<Integer Dim>
    class Simplex<Dim, 4> {
    public:
        std::array<Integer, 5> sides;
        std::array<Integer, 5> nodes;
        Integer id;
    };
    
    template<Integer Dim, Integer ManifoldDim>
    inline constexpr static Integer n_nodes(const Simplex<Dim, ManifoldDim> &)
    {
        return ManifoldDim + 1;
    }
    
    template<Integer Dim, Integer ManifoldDim>
    inline constexpr static Integer n_sides(const Simplex<Dim, ManifoldDim> &)
    {
        return ManifoldDim + 1;
    }
    
    template<Integer Dim, Integer ManifoldDim>
    inline constexpr static Integer n_dims(const Simplex<Dim, ManifoldDim> &)
    {
        return Dim;
    }
    
    template<Integer Dim, Integer ManifoldDim>
    inline void jacobian(const Simplex<Dim, ManifoldDim>  &simplex,
                         const std::vector<Vector<Real, Dim>> &points,
                         Matrix<Real, Dim, ManifoldDim> &J)
    {
        static_assert(Dim >= ManifoldDim, "Dim must be greater or equal ManifoldDim");
        
        auto n = n_nodes(simplex);
        
        Vector<Real, Dim> v0 = points[simplex.nodes[0]];
        
        for(Integer i = 1; i < n; ++i) {
            const auto &vi = points[simplex.nodes[i]];
            J.col(i-1, vi - v0);
        }
    }

    template<Integer Dim, Integer ManifoldDim>
    bool check_and_fix_jac(Matrix<Real, Dim, ManifoldDim> &J)
    {
        Integer n_zero_rows = 0;
        for(Integer i = 0; i < Dim; ++i) {

            bool is_zero_row = true;
            for(Integer j = 0; j < ManifoldDim; ++j) {
                if(std::abs(J(i, j)) != 0.) {
                    is_zero_row = false;
                    break;
                }
            }

            if(is_zero_row) {
                ++n_zero_rows;
                
                if(i < ManifoldDim) {
                    J(i, i) = 1.;
                }
            }
        }

        return n_zero_rows <= (Dim - ManifoldDim);
    }
    
    template<Integer Dim, Integer ManifoldDim>
    inline Real volume(const Simplex<Dim, ManifoldDim>  &simplex,
                       const std::vector<Vector<Real, Dim>> &points)
    {
        static_assert(Dim >= ManifoldDim, "Dim must be greater or equal ManifoldDim");
        
        Matrix<Real, Dim, ManifoldDim> J;
        jacobian(simplex, points, J);

        //this hack does not work (must find submanifold plane and compute volume there)
        if(!check_and_fix_jac(J)) {
            return 0.;
        }

        auto ref_vol = (1./Factorial<ManifoldDim>::value);
        return ref_vol * det(J);
    }

    inline Vector4r normal(
    	const Vector4r &x0,
    	const Vector4r &x1,
    	const Vector4r &x2,
    	const Vector4r &x3) {

      	Vector4r ret;

      	//row 1
      	const Real m10 = x1(0) - x0(0);
      	const Real m11 = x1(1) - x0(1);
      	const Real m12 = x1(2) - x0(2);
      	const Real m13 = x1(3) - x0(3);

      	//row 2
      	const Real m20 = x2(0) - x0(0);
      	const Real m21 = x2(1) - x0(1);
      	const Real m22 = x2(2) - x0(2);
      	const Real m23 = x2(3) - x0(3);

      	//row 3
      	const Real m30 = x3(0) - x0(0);
      	const Real m31 = x3(1) - x0(1);
      	const Real m32 = x3(2) - x0(2);
      	const Real m33 = x3(3) - x0(3);

      	ret[0] = det(
            Matrix<Real, 3, 3>({
                m11, m12, m13,
                m21, m22, m23,
                m31, m32, m33
        }));
        
        ret[1] = -det(
            Matrix<Real, 3, 3>({
                m10, m12, m13,
                m20, m22, m23,
                m30, m32, m33
        }));
        
        ret[2] = det(
            Matrix<Real, 3, 3>({
                m10, m11, m13,
                m20, m21, m23,
                m30, m31, m33
        }));
        
        ret[3] = -det(
            Matrix<Real, 3, 3>({
                m10, m11, m12,
                m20, m21, m22,
                m30, m31, m32
        }));

        ret.normalize();

        return ret;
    }

    template<Integer Dim, Integer ManifoldDim>
    inline Vector<Real, Dim> normal(
        const Simplex<Dim, ManifoldDim>      &simplex,
        const std::vector<Vector<Real, Dim>> &points,
        const bool apply_normalization = true)
    {
        static_assert(Dim >= ManifoldDim, "simplex must be embedded in R^Dim");

        Matrix<Real, ManifoldDim, Dim> m;
        Matrix<Real, ManifoldDim, Dim-1> minor;
        m.zero();

        auto x0 = points[simplex.nodes[0]];
        for(Integer i = 0; i < ManifoldDim; ++i) {
            auto xi = points[simplex.nodes[i + 1]];

            for(Integer j = 0; j < Dim; ++j) {
                m(i, j) = xi(j) - x0(j);
            }
        }

        Vector<Real, Dim> ret;
        for(Integer d = 0; d < Dim; ++d) {
            for(Integer i = 0; i < ManifoldDim; ++i) {
                Integer k = 0;
                for(Integer j = 0; j < Dim; ++j) {
                    if(j == d) { continue; }
                    minor(i, k) = m(i, j);
                    k++;
                }
            }

            auto det_minor = det(minor);
            ret[d] = ((d & 1) == 0 ? 1. : -1.) * det_minor;
        }

        if(apply_normalization) {
            ret.normalize();
        }
        return ret;
    }

    inline Vector4r normal(
        const Vector4r &x0,
        const Vector4r &x1,
        const Vector4r &x2
    ) 
    {
        return Vector4r();
    }
}

#endif //MARS_SIMPLEX_HPP
