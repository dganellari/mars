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
  static constexpr auto faces_combs=combinations_generate<ManifoldDim+1,ManifoldDim>();

  FiniteElem(const Mesh<Dim,ManifoldDim>&mesh):
  already_set_(false),
  mesh_ptr_(std::make_shared<Mesh<Dim,ManifoldDim>>(mesh))
  {}

  FiniteElem(const std::shared_ptr<Mesh<Dim,ManifoldDim>>&mesh_ptr):
  already_set_(false),
  mesh_ptr_(mesh_ptr)
  {}

  constexpr void init_det(){detJ_=det(J_);}
  constexpr void init_side_det(){det_side_J_=det(side_J_);}

  constexpr void init(const Integer id)
   {
   if(already_set_==false)
    {
      points_.resize(ManifoldDim+1);
      side_points_.resize(ManifoldDim);

      already_set_=true;
    }
   const auto& points=mesh_ptr_->points();
   elem_id_=id;
   elem_=mesh_ptr_->elem(id);
   side_tags_=elem_.side_tags;
   n_sides_=n_sides(elem_);
   auto n = n_nodes(elem_);
   points_[0] = points[elem_.nodes[0]];
   for(Integer i = 1; i < n; ++i) 
      points_[i] = points[elem_.nodes[i]];
   affine_transformation(elem_,points,J_,v0_);
   // jacobian(elem_,points,J_);   
   init_det();
   volume_=unsigned_volume(elem_,detJ_);
   }
  constexpr void init(const Integer id, const Integer level)
   {
    level_=level;
    init(id);
   }
  constexpr void init_boundary(const Integer side_id,const bool also_inside=false)
  {
    // init must be called before init the boundary
    side_id_=side_id;
    side_tag_=side_tags_[side_id];
    if(side_tag_!=INVALID_INDEX || also_inside)
    {
     const auto& points=mesh_ptr_->points();
     elem_.side_sorted(side_id,side_);
     // elem_.side(side_id,side_);
     // side_points_[0] = points_[faces_combs[side_id][0]];
     // for(Integer i = 1; i < ManifoldDim; ++i) 
     //    side_points_[i] = points_[faces_combs[side_id][i]];
    // std::cout<<"side nodes="<<std::endl;
    //  for(Integer i = 0; i < side_.nodes.size(); ++i) 
    //     std::cout<<side_.nodes[i]<<std::endl;

    // std::cout<<"points_"<<std::endl;
    //  for(Integer i = 0; i < points_.size(); ++i) 
    //     std::cout<<points_[i]<<std::endl;


     side_points_[0] = points[side_.nodes[0] ];
     for(Integer i = 1; i < ManifoldDim; ++i) 
        side_points_[i] = points[side_.nodes[i] ];

     // std::cout<<"side points="<<std::endl;
     // for(Integer i = 0; i < ManifoldDim; ++i) 
     //    std::cout<<side_points_[i]<<std::endl;

      // std::cout<<"side_id="<<side_id<<std::endl;
      // std::cout<<"side_points_=" <<std::endl;
      // for(Integer ii=0;ii<side_points_.size();ii++)
      //   std::cout<<side_points_[ii]<<std::endl;
      // jacobian(side_,points_,side_J_);
      affine_transformation(side_,points,side_J_,v0_);
      // std::cout<<"init boundary"<<std::endl;
      // for(Integer k = 0; k < side_.nodes.size(); ++k)
      // std::cout<<side_.nodes[k]<<" "<<std::endl;
      // std::cout<<std::endl;

      // std::cout<<"side_J_"<<side_J_<<std::endl;
      // std::cout<<"v0_"<<v0_<<std::endl;
      // for(Integer k = 0; k < points_.size(); ++k) 
      // std::cout<<"points_="<<points_[k]<<std::endl;
      // std::cout<<"side_J_="<<side_J_<<std::endl;
      // std::cout<<"v0="<<v0_<<std::endl;
      // std::cout<<"side_J_="<<std::endl;
      // std::cout<<side_J_ <<std::endl;

      init_side_det();
      // std::cout<<"det_side_J_="<<std::endl;
      // std::cout<<det_side_J_ <<std::endl;

      side_volume_=unsigned_volume(side_,det_side_J_);
      // std::cout<<"side_volume_="<<std::endl;
      // std::cout<<side_volume_ <<std::endl;

      // std::cout<<"detJ_="<<std::endl;
      // std::cout<<detJ_ <<std::endl;
      // std::cout<<"volume_="<<std::endl;
      // std::cout<<volume_ <<std::endl;
    }

  }
  constexpr inline void init_inv_jac() {J_inv_=inverse(J_);}

  constexpr auto & inv_jac()const {return J_inv_;}

  constexpr auto & operator()()const {return J_;}
  
  constexpr auto & jac()const {return J_;}

  constexpr auto & jac_side()const {return side_J_;}

  constexpr auto & v0()const {return v0_;}

  inline auto & mesh_ptr()const {return mesh_ptr_;}

  // constexpr auto Matrix<Real,Dim,1>& transform_point(const Matrix<Real,Dim,1>& Point)const 
  // {

  //   return v0_;
  // }

  constexpr void transform_point(Vector<Real,Dim>& vec_transformed_point, const Matrix<Real,Dim,1>& Point) const
  {
    for(Integer i=0; i<Dim;i++)
     {
       vec_transformed_point[i]=v0_[i];
       for(Integer j=0;j<J_.cols();j++)
        vec_transformed_point[i]+=J_(i,j)*Point[j];
     }  

    // return mat_transformed_point_;
  }

  constexpr void transform_point(Matrix<Real,1,Dim>& mat_transformed_point, const Matrix<Real,1,Dim>& Point) const
  {
    for(Integer i=0; i<Dim;i++)
     {
       mat_transformed_point(0,i)=v0_[i];
       for(Integer j=0;j<J_.cols();j++)
        mat_transformed_point(0,i)+=J_(i,j)*Point(0,j);
     }  
  }


  constexpr void transform_point(Matrix<Real,Dim,1>& mat_transformed_point, const Matrix<Real,1,Dim>& Point) const
  {
    // std::cout<<"v0="<<v0_<<std::endl;
    // std::cout<<"J_="<<J_<<std::endl;

    for(Integer i=0; i<Dim;i++)
     {
       mat_transformed_point[i]=v0_[i];
       for(Integer j=0;j<J_.cols();j++)
        mat_transformed_point[i]+=J_(i,j)*Point(0,j);
     }  
  }


  template<Integer DimAux>
  constexpr void transform_point(Matrix<Real,DimAux,Dim>& mat_transformed_point, const Matrix<Real,DimAux,Dim>& Point) const
  {
    
    // std::cout<<"transform_point" <<std::endl;

    // std::cout<<"J_" <<std::endl;
    // std::cout<<J_ <<std::endl;

    // std::cout<<"v0_" <<std::endl;
    // std::cout<<v0_ <<std::endl;

    // std::cout<<"Point" <<std::endl;
    // std::cout<<Point <<std::endl;

    for(Integer s=0;s<DimAux;s++)
    {
    for(Integer i=0; i<J_.rows();i++)
     {
       mat_transformed_point(s,i)=v0_[i];
       for(Integer j=0;j<Dim;j++)
        mat_transformed_point(s,i)+=J_(i,j)*Point(s,j);
     }
    }
  
  }

  template<Integer DimAux,Integer DimFrom>
  constexpr void side_transform_point(Matrix<Real,DimAux,Dim>& mat_transformed_point, const Matrix<Real,DimAux,DimFrom>& Point) const
  {

    // std::cout<<"side_J_"<<std::endl;
    // std::cout<<side_J_<<std::endl;
    // std::cout<<"Point"<<std::endl;
    // std::cout<<Point<<std::endl;
    for(Integer s=0;s<DimAux;s++)
    {
      // std::cout<<"s="<<s<<std::endl;
    for(Integer i=0; i<Dim;i++)
     {
       mat_transformed_point(s,i)=v0_[i];
       for(Integer j=0;j<side_J_.cols();j++)
        mat_transformed_point(s,i)+=side_J_(i,j)*Point(s,j);
     }
    }
  
  }
  
  template<Integer Rows>
  inline void points_to_reference_points(Matrix<Real,Rows,Dim>& ref,const Matrix<Real,Rows,Dim>& conf)const
  {
 

    std::cout<<"v0_"<<std::endl;
    std::cout<<v0_<<std::endl;
    std::cout<<"J_inv_"<<std::endl;
    std::cout<<J_inv_<<std::endl;
   for(Integer i=0;i<Rows;i++)
    {
      for(Integer k=0;k<Dim;k++)
      { 
        ref(i,k)=J_inv_(k,0)*(conf(i,0)-v0_[0]);
        for(Integer j=1;j<Dim;j++)
        {
                ref(i,k)+=J_inv_(k,j)*(conf(i,j)-v0_[j]);
        }

      }


    }
  }



  const auto& level()   const {return level_;}


  constexpr auto get_det()   const {return detJ_;}

  constexpr auto get_det_side()   const {return det_side_J_;}

  template<bool VolumeDet>
  constexpr std::enable_if_t<(VolumeDet==true),Real> 
  get_det()   const {return detJ_;}

  template<bool VolumeDet>
  constexpr std::enable_if_t<(VolumeDet==false),Real> 
  get_det()   const {return det_side_J_;}

  constexpr auto & elem_id()   const {return elem_id_;}
  constexpr auto & side_id()   const {return side_id_;}
  constexpr auto & side_tag()   const {return side_tag_;}
  constexpr auto & side_tags()   const {return side_tags_;}
  constexpr auto & n_side()   const {return n_sides_;}
  constexpr auto & side()   const {return side_;}
  constexpr auto & side_volume() const {return side_volume_;}
  constexpr auto & points() const {return points_;}
  constexpr auto & volume() const {return volume_;}

  constexpr bool is_on_boundary()   const 
  {
    for(std::size_t i=0;i<n_sides_;i++)
      if(side_tags_[i]!=INVALID_INDEX)
        return true;
    return false;
  }
  constexpr bool is_side_on_boundary()   const 
  {
      if(side_tag_==INVALID_INDEX)
        return false;
      else
        return true;
  }

  private:  
  bool already_set_;

  Integer elem_id_;
  Simplex<Dim, ManifoldDim> elem_;
  std::vector<Vector<Real,Dim>> points_;
  Matrix<Real, Dim, ManifoldDim> J_;
  Matrix<Real, Dim, ManifoldDim> J_inv_;
  Real detJ_;
  Real volume_;
  // Matrix<Real,Dim,1> mat_transformed_point_;
  // Vector<Real,Dim> vec_transformed_point_;

  
  Integer n_sides_;
  Integer side_id_;
  Integer side_tag_;
  std::array<Integer, ManifoldDim+1> side_tags_;
  Simplex<Dim, ManifoldDim-1> side_;
  std::vector<Vector<Real,Dim>> side_points_;
  Matrix<Real, Dim, ManifoldDim-1> side_J_;
  Real det_side_J_;
  Real side_volume_;
  Vector<Real, Dim> v0_;
  std::shared_ptr<Mesh<Dim,ManifoldDim>> mesh_ptr_;

  Integer level_;



};



}
#endif