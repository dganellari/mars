#ifndef MARS_JACOBIAN_HPP
#define MARS_JACOBIAN_HPP
#include "mars_simplex.hpp"
#include "mars_general_form.hpp"

namespace mars{




template<typename Elem>
class FiniteElem;

template<Integer Dim, Integer ManifoldDim>
class FiniteElem<Simplex<Dim, ManifoldDim>>
{
  public:


  FiniteElem(const Mesh<Dim,ManifoldDim>&mesh):
  already_set_(false),
  mesh_ptr_(std::make_shared<Mesh<Dim,ManifoldDim>>(mesh))
  {}
  
  constexpr void init_det(){detJ_=det(J_);}

  constexpr void init(const Integer id)
   {
   if(already_set_==false)
    {
      points_.resize(ManifoldDim+1);

      already_set_=true;
    }
   const auto& points=mesh_ptr_->points();
   elem_id_=id;
   elem_=mesh_ptr_->elem(id);
   n_sides_=n_sides(elem_);
   auto n = n_nodes(elem_);
   points_[0] = points[elem_.nodes[0]];
   for(Integer i = 1; i < n; ++i) 
      points_[i] = points[elem_.nodes[i]];
   jacobian(elem_,points,J_);
   volume_=unsigned_volume(elem_,points_);
   init_det();
   }

  constexpr void init_boundary(const Integer side_id)
  {
    // init must be called before init the boundary
    side_id_=side_id;
    elem_.side(side_id,side_);
    side_volume_=unsigned_volume(side_,points_);
  }

  constexpr auto & operator()()const {return J_;}
  constexpr auto & get_det()   const {return detJ_;}
  constexpr auto & elem_id()   const {return elem_id_;}
  constexpr auto & side_id()   const {return side_id_;}
  constexpr auto & n_side()   const {return n_sides_;}
  constexpr auto & side()   const {return side_;}
  constexpr auto & side_volume() const {return side_volume_;}
  constexpr auto & points() const {return points_;}
  constexpr auto & volume() const {return volume_;}

  private:  
  Integer n_sides_;
  Simplex<Dim, ManifoldDim-1> side_;
  Real side_volume_;
  Integer side_id_;
  bool already_set_;
  Integer elem_id_;
  Simplex<Dim, ManifoldDim> elem_;
  Real volume_;
  std::vector<Vector<Real,Dim>> points_;
  Matrix<Real, Dim, ManifoldDim> J_;
  Real detJ_;
  std::shared_ptr<Mesh<Dim,ManifoldDim>> mesh_ptr_;

};



}
#endif