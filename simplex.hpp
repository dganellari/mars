#ifndef MARS_SIMPLEX_HPP
#define MARS_SIMPLEX_HPP

#include "base.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "static_math.hpp"

#include <array>
#include <vector>
#include <ostream>
#include <iostream>

#include <cmath>
#include <cassert>
#include <initializer_list>
#include <algorithm>

namespace mars {

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

        inline static std::vector<Vector<Real, Dim>> &ref()
        {
            static std::vector<Vector<Real, Dim>> ref_(1, Vector<Real, Dim>().zero());
            return ref_;
        }
    };
    
    template<Integer Dim>
    class Simplex<Dim, 1> {
    public:
        std::array<Integer, 2> nodes;
        Integer id;

        inline static std::vector<Vector<Real, Dim>> &ref()
        {
            static std::vector<Vector<Real, Dim>> ref_;
            if(ref_.empty()) {
                ref_.resize(2);
                ref_[0] = Vector<Real, Dim>().zero();
                ref_[1] = Vector<Real, Dim>().zero();
                ref_[1](0) = 1.;
            }
            return ref_;
        }
    };
    
    template<Integer Dim>
    class Simplex<Dim, 2> {
    public:
        std::array<Integer, 3> sides;
        std::array<Integer, 3> nodes;
        Integer id;


        inline static std::vector<Vector<Real, Dim>> &ref()
        {
            static std::vector<Vector<Real, Dim>> ref_;
            if(ref_.empty()) {
                ref_.resize(3);
                ref_[0] = Vector<Real, Dim>().zero();
                ref_[1] = Vector<Real, Dim>().zero();
                ref_[1](0) = 1.;

                ref_[2] = Vector<Real, Dim>().zero();
                ref_[2](1) = 1.;
            }

            return ref_;
        }
    };
    
    template<Integer Dim>
    class Simplex<Dim, 3> {
    public:
        std::array<Integer, 4> sides;
        std::array<Integer, 4> nodes;
        Integer id;


        inline static std::vector<Vector<Real, Dim>> &ref()
        {
            static std::vector<Vector<Real, Dim>> ref_;
            if(ref_.empty()) {
                ref_.resize(4);
                ref_[0] = Vector<Real, Dim>().zero();
               
                ref_[1] = Vector<Real, Dim>().zero();
                ref_[1](0) = 1.;

                ref_[2] = Vector<Real, Dim>().zero();
                ref_[2](1) = 1.;

                ref_[3] = Vector<Real, Dim>().zero();
                ref_[3](2) = 1.;
            }
            
            return ref_;
        }
    };
    
    template<Integer Dim>
    class Simplex<Dim, 4> {
    public:
        std::array<Integer, 5> sides;
        std::array<Integer, 5> nodes;
        Integer id;


        inline static std::vector<Vector<Real, Dim>> &ref()
        {
            static std::vector<Vector<Real, Dim>> ref_;
            if(ref_.empty()) {
                ref_.resize(5);
                ref_[0] = Vector<Real, Dim>().zero();
               
                ref_[1] = Vector<Real, Dim>().zero();
                ref_[1](0) = 1.;

                ref_[2] = Vector<Real, Dim>().zero();
                ref_[2](1) = 1.;

                ref_[3] = Vector<Real, Dim>().zero();
                ref_[3](2) = 1.;

                ref_[4] = Vector<Real, Dim>().zero();
                ref_[4](3) = 1.;
            }
            
            return ref_;
        }
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

    template<Integer ManifoldDim>
    inline Integer midpoint_index(
        const Integer i,
        const Integer j)
    {
        const auto ip1 = i + 1;
        const auto jp1 = j + 1;
        return ((ip1 - 1) * (ManifoldDim - (ip1/2.)) + jp1 + ManifoldDim) - 1;
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
        // if(!check_and_fix_jac(J)) {
        //     return 0.;
        // }

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


    template<Integer ManifoldDim>
    inline void red_refinement_interpolator(
        Matrix<Real, 
               ManifoldDim + 1 + Combinations<ManifoldDim + 1, 2>::value,
               ManifoldDim + 1> &interp)
    {
        interp.zero();

        for(Integer i = 0; i < ManifoldDim + 1; ++i) {
            interp(i, i) = 1.;
        }

        for(Integer i = 0; i < ManifoldDim + 1; ++i) {
            for(Integer j = i + 1; j < ManifoldDim + 1; ++j) {
                Integer offset = midpoint_index<ManifoldDim>(i, j); 
                interp(offset, i) = 0.5;
                interp(offset, j) = 0.5;
            }
        }
    }
}

#endif //MARS_SIMPLEX_HPP
