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
        Integer id = INVALID_INDEX;
        Integer parent_id = INVALID_INDEX;

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
        Integer id = INVALID_INDEX;
        Integer parent_id = INVALID_INDEX;

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
        // std::array<Integer, 3> sides;
        std::array<Integer, 3> nodes;
        Integer id = INVALID_INDEX;
        Integer parent_id = INVALID_INDEX;


        void side(
            const Integer &side_num,
            Simplex<Dim, 1> &side) const
        {
            side.nodes[0] = nodes[side_num];
            side.nodes[1] = nodes[side_num == 2? 0 : (side_num + 1)];
        }

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
        // std::array<Integer, 4> sides;
        std::array<Integer, 4> nodes;
        Integer id = INVALID_INDEX;
        Integer parent_id = INVALID_INDEX;

        void side(
            const Integer &side_num,
            Simplex<Dim, 2> &side) const
        {
            switch(side_num) {
                case 0:
                {   
                    side.nodes[0] = nodes[0];
                    side.nodes[1] = nodes[2];
                    side.nodes[2] = nodes[1];
                    break;
                }

                case 1:
                {
                    side.nodes[0] = nodes[0];
                    side.nodes[1] = nodes[3];
                    side.nodes[2] = nodes[2];
                    break;
                }

                case 2:
                {
                    side.nodes[0] = nodes[0];
                    side.nodes[1] = nodes[1];
                    side.nodes[2] = nodes[3];
                    break;
                }

                case 3:
                {
                    side.nodes[0] = nodes[1];
                    side.nodes[1] = nodes[2];
                    side.nodes[2] = nodes[3];
                    break;
                }

                default:
                {
                    assert(false);
                    break;
                }

            }
        }


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
        // std::array<Integer, 5> sides;
        std::array<Integer, 5> nodes;
        Integer id = INVALID_INDEX;
        Integer parent_id = INVALID_INDEX;


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

        void side(
            const Integer &side_num,
            Simplex<Dim, 3> &side) const
        {
            switch(side_num)
            {
                case 0:
                {
                    side.nodes[0] = nodes[0];
                    side.nodes[1] = nodes[1];
                    side.nodes[2] = nodes[2];
                    side.nodes[3] = nodes[3];
                    break;
                }

                case 1:
                {   
                    side.nodes[0] = nodes[0];
                    side.nodes[1] = nodes[1];
                    side.nodes[2] = nodes[3];
                    side.nodes[3] = nodes[4];
                    break;
                }

                case 2:
                {   
                    side.nodes[0] = nodes[0];
                    side.nodes[1] = nodes[3];
                    side.nodes[2] = nodes[2];
                    side.nodes[3] = nodes[4];
                    break;
                }

                case 3:
                {   
                    side.nodes[0] = nodes[0];
                    side.nodes[1] = nodes[1];
                    side.nodes[2] = nodes[4];
                    side.nodes[3] = nodes[2];
                    break;
                }

                case 4:
                {   
                    side.nodes[0] = nodes[4];
                    side.nodes[1] = nodes[1];
                    side.nodes[2] = nodes[3];
                    side.nodes[3] = nodes[2];
                    break;
                }

                default: {
                    assert(false);
                    break;
                }
            }
        }
    };

    template<Integer ManifoldDim>
    class NSubSimplices {
    public:
        static const Integer value = Power<2, ManifoldDim>::value;
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

    template<Integer Dim>
    inline void fixed_red_refinement(std::array<Simplex<Dim, 1>, 2> &sub_simplices)
    {
        sub_simplices[0].nodes = {0, 2};
        sub_simplices[1].nodes = {2, 1};
    }

    template<Integer Dim>
    inline void fixed_red_refinement(std::array<Simplex<Dim, 2>, 4> &sub_simplices)
    {
        sub_simplices[0].nodes = {0, 3, 4};
        sub_simplices[1].nodes = {1, 5, 3};
        sub_simplices[2].nodes = {2, 4, 5};
        sub_simplices[3].nodes = {3, 5, 4};
    }

    template<Integer Dim>
    inline void fixed_red_refinement(std::array<Simplex<Dim, 3>, 8> &sub_simplices)
    {
        //corner tets
        sub_simplices[0].nodes = {0, 4, 5, 6};
        sub_simplices[1].nodes = {4, 1, 7, 8};
        sub_simplices[2].nodes = {5, 7, 2, 9};
        sub_simplices[3].nodes = {6, 8, 9, 3};

        //octahedron tets
        sub_simplices[4].nodes = {4, 7, 5, 8};
        sub_simplices[5].nodes = {4, 8, 5, 6};
        sub_simplices[6].nodes = {6, 8, 5, 9};
        sub_simplices[7].nodes = {7, 6, 5, 9};
    }

    inline void fixed_red_refinement(std::array<Simplex<4, 4>, 16> &sub_simplices)
    {
        sub_simplices[0].nodes  = {0, 5, 6, 7, 8};
        sub_simplices[1].nodes  = {5, 1, 9, 10, 11};
        
        sub_simplices[2].nodes  = {6, 9, 2, 12, 13};
        sub_simplices[3].nodes  = {7, 10, 12, 3, 14};
        
        sub_simplices[4].nodes  = {8, 11, 13, 14, 4}; 
        //wrong
        // sub_simplices[5].nodes  = {5, 8, 9, 12, 14};
        sub_simplices[5].nodes  = {5, 8, 9, 14, 12};

        //wrong
        // sub_simplices[6].nodes  = {5, 7, 8, 12, 14};
        sub_simplices[6].nodes  = {5, 7, 8, 14, 12};
        sub_simplices[7].nodes  = {6, 8, 9, 12, 13};

        //wrong
        // sub_simplices[8].nodes  = {5, 6, 7, 8, 12};
        sub_simplices[8].nodes  = {5, 6, 7, 12, 8};
        sub_simplices[9].nodes  = {8, 9, 12, 14, 13};

        sub_simplices[10].nodes = {5, 8, 9, 11, 14};
        sub_simplices[11].nodes = {5, 7, 10, 12, 14};

        //wrong
        // sub_simplices[12].nodes = {5, 9, 10, 12, 14};
        sub_simplices[12].nodes = {5, 9, 10, 14, 12};
        //wrong
        // sub_simplices[13].nodes = {5, 6, 8, 9, 12};
        sub_simplices[13].nodes = {5, 6, 8, 12, 9};

        sub_simplices[14].nodes = {8, 9, 11, 13, 14};
        //wrong
        // sub_simplices[15].nodes = {5, 9, 10, 12, 14};
        sub_simplices[15].nodes = {5, 9, 10, 14, 12};
    }

    inline void non_cyclic_fixing(std::array<Simplex<4, 4>, 16> &sub_simplices)
    {
        sub_simplices[0].nodes  = {0, 5, 7, 8, 9};
        sub_simplices[1].nodes  = {5, 1, 9, 10, 11};
       
        sub_simplices[2].nodes  = {6, 9, 2, 12, 13};
        sub_simplices[3].nodes  = {7, 10, 12, 3, 14};
       
        sub_simplices[4].nodes  = {8, 11, 13, 14, 4}; 
        //wrong
        sub_simplices[5].nodes  = {7, 8, 10, 12, 14};

        sub_simplices[6].nodes  = {6, 8, 9, 12, 13};
        sub_simplices[7].nodes  = {6, 7, 8, 9, 12};

        //wrong
        sub_simplices[8].nodes  = {8, 9, 12, 13, 14};
        //wrong
        sub_simplices[9].nodes  = {5, 8, 9, 10, 11};

        //wrong
        sub_simplices[10].nodes = {5, 7, 8, 9, 10};
        sub_simplices[11].nodes = {8, 9, 10, 11, 14};

        sub_simplices[12].nodes = {5, 6, 7, 8, 9};
        sub_simplices[13].nodes = {8, 9, 11, 13, 14};

        //wrong
        sub_simplices[14].nodes = {8, 9, 10, 12, 14};
        sub_simplices[15].nodes = {7, 8, 9, 10, 12};
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


    template<Integer Dim, Integer ManifoldDim>
    inline Real unsigned_volume(const Simplex<Dim, ManifoldDim>  &simplex,
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
        
        if(Dim > ManifoldDim) {
          Matrix<Real, ManifoldDim, ManifoldDim> JtJ = transpose(J) * J;
          return ref_vol * std::sqrt(det(JtJ));
        } else {
            return ref_vol * std::abs(det(J));
        }
    }

    template<Integer Dim, Integer ManifoldDim>
    inline Vector<Real, Dim> barycenter(
        const Simplex<Dim, ManifoldDim>      &simplex,
        const std::vector<Vector<Real, Dim>> &points)
    {
        Vector<Real, Dim> b;
        b.zero();

        for(Integer i = 0; i < n_nodes(simplex); ++i) {
            b += points[simplex.nodes[i]];
        }

        b /= n_nodes(simplex);
        return b;
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
    using SimplexInterpolator = Matrix< Real,
                                        ManifoldDim + 1 + Combinations<ManifoldDim + 1, 2>::value,
                                        ManifoldDim + 1>; 


    template<Integer ManifoldDim>
    inline void red_refinement_interpolator(
        SimplexInterpolator<ManifoldDim> &interp)
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

    template<Integer Dim, Integer ManifoldDim>
    void refine_points(
        const std::vector<Vector<Real, Dim>> &parent_points,
        const SimplexInterpolator<ManifoldDim> &interp,
        std::vector<Vector<Real, Dim>> &children_points)
    {
        children_points.resize(interp.rows());
        for(Integer i = 0; i < interp.rows(); ++i) {
            children_points[i].zero();
            for(Integer j = 0; j < interp.cols(); ++j) {
                children_points[i] += interp(i, j) * parent_points[j];
            }
        }
    }

    template<Integer Dim, Integer ManifoldDim, Integer NSubs>
    void red_refinement(
        const Simplex<Dim, ManifoldDim> &parent,
        const std::vector<Vector<Real, Dim>> &parent_points,
        std::array<Simplex<Dim, ManifoldDim>, NSubs> &children, 
        std::vector<Vector<Real, Dim>> &children_points,
        SimplexInterpolator<ManifoldDim> &interp)
    {
        red_refinement_interpolator<ManifoldDim>(interp);
        refine_points<Dim, ManifoldDim>(parent_points, interp, children_points);
        fixed_red_refinement(children);

        for(auto &c : children) {
            c.parent_id = parent.id;
        }
    }
}

#endif //MARS_SIMPLEX_HPP
