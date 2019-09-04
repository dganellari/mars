




#ifndef MARS_SHAPE_FUNCTION_HPP
#define MARS_SHAPE_FUNCTION_HPP

#include "mars_simplex.hpp"
#include "mars_base_elementfunctionspace.hpp"
#include "mars_vector.hpp"
#include "mars_fqpexpressions.hpp"
#include "mars_quadrature_rules.hpp"
#include "mars_operators.hpp"
#include "mars_referencemap.hpp"
#include "mars_shape_function_coefficients.hpp"

namespace mars{





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////// class SignedNormal: for each face f of the element e                             ///////////////////
////////////// normal[e][f] returns the normal, outward[e][f] returns the sign of the normal    ///////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// returns the sign of val (1,0,-1)
// template <typename T> int sgn(T val) {
//     return (T(0) < val) - (val < T(0));
// }
// returns the order of the indeces of the sorted v
template <typename T>
T argsort(const T &v) {
  // initialize original index locations
  T idx;
  std::iota(idx.begin(), idx.end(), 0);
  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}
// returns v, resorted based on index 
template <typename T, typename S>
T sort_by_index(const T &v,const S& index) {
  T sortv;
  assert(index.size()==v.size() && "sort_by_index: v and index must have the same length, otherwise last elements in v are not initialized");
  for(Integer ii=0;ii<index.size();ii++)
      sortv[ii]=v[index[ii]];
  return sortv;
}
// returns the sub array of v of indices = index
template <typename S,std::size_t SDim, typename T,std::size_t TDim>
std::array<S,TDim> sub_array(const std::array<S,SDim> &v,const std::array<T,TDim>& index) {

  static_assert(TDim<=SDim,"in sub_array the length of the index vector must be smaller than the one of the vector");
  std::array<S,TDim> subvec;
  for(Integer ii=0;ii<TDim;ii++)
      subvec[ii]=v[index[ii]];

  return subvec;
}




template<typename Elem>
class SignedNormal;


template<Integer Dim, Integer ManifoldDim>
class SignedNormal<Simplex<Dim,ManifoldDim>>{
private:
      std::vector< Vector< Vector< Real, Dim> , ManifoldDim + 1 > >  normal_;
      std::vector< Vector<Real, ManifoldDim + 1 > > outward_;
public:

      std::vector< Vector< Vector< Real, Dim> , ManifoldDim + 1 > > operator () () const { return normal_; };
      // const Vector< Vector< Real, Dim> , ManifoldDim + 1 >& normal(const Integer elem_id) const {return normal_[elem_id];}; 
      std::vector< Vector<Real, ManifoldDim + 1 > > sign() const {return outward_;}; 
      const Vector<Real, ManifoldDim + 1 >& sign(const Integer elem_id) const {return outward_[elem_id];}; 
      
      template< typename MeshT>
      void init(const MeshT& mesh)
       {
      using Elem=typename MeshT::Elem;
      static_assert(Dim==ManifoldDim, "SignedNormal: computing internal normals of a mesh makes sense only if Dim==ManifoldDim");   
      Integer n_elements=mesh.n_elements();
      Integer facenodes_carray[ManifoldDim];
      std::array<Integer,ManifoldDim+1> elemnodes_local;
      std::array<Integer,ManifoldDim+1> elemnodes_global;
      std::array<Integer,ManifoldDim+1> elemnodes_local_sort;

      std::array< Integer, ManifoldDim> facenodes_local;
      std::array< Integer, ManifoldDim> facenodes_local_sort;
      std::array<Integer,ManifoldDim> facenodes_global;
      std::array<Integer,ManifoldDim> facenodes_global_sort;

      std::array<Vector<Real,Dim>, ManifoldDim+1> points;
      Vector<Vector<Real,Dim>, ManifoldDim> facepoints;
      Vector<Vector<Real,Dim>, ManifoldDim+1> elempoints_sort;

      Vector<Real,Dim> facepoint_mean;
      Vector<Real,Dim> elempoint_mean;

      // we initialize the simplex
      Simplex<Dim,ManifoldDim> simplex_elem;
      for(Integer nn=0;nn<ManifoldDim+1;nn++)
          simplex_elem.nodes[nn]=nn;

      normal_.resize(n_elements);
      outward_.resize(n_elements);


      // elemnodes_local is the array{ManifoldDim,ManifoldDim-1,...,1,0}
      for(Integer nn=0;nn<ManifoldDim+1 ;nn++)
          elemnodes_local[nn]=nn;
      std::reverse(std::begin(elemnodes_local), std::end(elemnodes_local));

      // loop on all the elements
      for(Integer ee=0;ee<mesh.n_elements();ee++)
        {
        Elem elem=mesh.elem(ee);
        elemnodes_global=elem.nodes;
        
        for(Integer mm=0;mm<points.size();mm++)
            points[mm]=mesh.point(elemnodes_global[mm]);
         
        // loop on all the faces of the simplex
        // for Simplex={0,1,2,3}, the order is: 0-1-2, 0-1-3, 0-2-3, 1-2-3
        for(Integer mm=0;mm<ManifoldDim+1;mm++)
           {
            // facenodes_local is a std::array containing the local nodes of the face mm
            Combinations<ManifoldDim + 1,ManifoldDim>::generate(mm,facenodes_carray);
            std::copy(std::begin(facenodes_carray), std::end(facenodes_carray), std::begin(facenodes_local));  
            // we reorder the local nodes based on the sorting order of the corresponding global nodes
            facenodes_global=sub_array(elemnodes_global,facenodes_local);
            facenodes_global_sort=argsort(facenodes_global);
            facenodes_local_sort=sort_by_index(facenodes_local,facenodes_global_sort);

            // in elemnodes_local_sort we have facenodes_local_sort and the remaining node, elemnodes_local[mm] 
            for(Integer nn=0;nn<ManifoldDim;nn++)
                elemnodes_local_sort[nn]=facenodes_local_sort[nn];

            elemnodes_local_sort[ManifoldDim]=elemnodes_local[mm];

            for(Integer nn=0;nn<ManifoldDim+1;nn++)
               elempoints_sort[nn]=points[elemnodes_local_sort[nn]];

            // we create the face midpoint, the element midpoint and their difference
            for(Integer nn=0;nn<ManifoldDim;nn++)
              facepoints[nn]=points[facenodes_local[nn]];

            facepoint_mean = facepoints.Tmean();
            elempoint_mean=elempoints_sort.Tmean();
            auto diff=facepoint_mean-elempoint_mean;

            auto n = normal(simplex_elem, elempoints_sort);
            outward_[ee][mm]=Sign(dot(diff,n)); 
            // save the normal and if its outward or inward 
            // if the face is on the boundary, then it must be outward  
            if(elem.side_tags[mm]>=0 && outward_[ee][mm]==-1)
            {
            // if the normal ins inward, force it to be outward
                normal_[ee][mm]=-n;
                outward_[ee][mm]=1;
            }
            else
            {
            normal_[ee][mm]=n;
            }
            }
            }
       };
            
      SignedNormal(){}

      template< typename MeshT>
      SignedNormal(const MeshT& mesh)
      {init(mesh);};

      template< typename MeshT>
      void print(const MeshT& mesh)
      {  
        using Elem=typename MeshT::Elem;
        Vector<Vector<Real,Dim>, ManifoldDim+1> elempoints_sort; 
        std::array<Integer,ManifoldDim+1> nodes;
        std::vector<Vector<Real,Dim>> points;
        Integer facenodes_carray[ManifoldDim];
        std::array< Integer, ManifoldDim> facenodes_local;

        for(Integer ee=0;ee<normal_.size();ee++)
        {
            std::cout<<std::endl<<"ELEMENT ID == "<<ee<<" WITH NODES"<<std::endl;
            Elem elem=mesh.elem(ee);
            nodes=elem.nodes;
            for(Integer mm=0;mm<ManifoldDim+1;mm++)
               std::cout<< nodes[mm] <<" ";
            std::cout<<std::endl;

            for(Integer mm=0;mm<ManifoldDim+1;mm++)
                {
                  Combinations<ManifoldDim + 1,ManifoldDim>::generate(mm,facenodes_carray);
                  std::copy(std::begin(facenodes_carray), std::end(facenodes_carray), std::begin(facenodes_local));    
                  std::cout<<"facenodes_local== "<<std::endl;
                  for(Integer ii=0;ii<ManifoldDim;ii++)
                      std::cout<<facenodes_local[ii]<<" ";           
                  std::cout<<"normal== ";
                  normal_[ee][mm].describe(std::cout);
                  std::cout<<"outward/inward(+1/-1)== "<<outward_[ee][mm]<<std::endl;
                }
           }
      };


};











































































template<typename Elem,typename BaseFunctionSpace,typename Operator, typename QuadratureRule>
class BaseShapeFunctionOperatorDependent;





template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace>
class BaseShapeFunctionOperatorDependent<Elem,BaseFunctionSpace,IdentityOperator,QuadratureRule>

{ 
  public:
  using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
  static constexpr Integer Dim=Elem::Dim;
  static constexpr Integer ManifoldDim=Elem::ManifoldDim;
  static constexpr Integer NComponents=BaseFunctionSpace::NComponents;
  static constexpr Integer ShapeFunctionDim1=FunctionSpace::ShapeFunctionDim1;
  static constexpr Integer ShapeFunctionDim2=FunctionSpace::ShapeFunctionDim2;
  static constexpr Integer NQPoints=QuadratureRule::NQPoints;
  static constexpr Integer Ntot=FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,BaseFunctionSpace>>::value;
  static constexpr Integer Ndofs=Ntot/NComponents;
  using single_type   = Matrix<Real, ShapeFunctionDim1, ShapeFunctionDim2>;
  using vector_single_type   = Vector<single_type,Ndofs>;
  using tot_type= Matrix<Real, ShapeFunctionDim1 * NComponents, ShapeFunctionDim2>;
  using type= FQPValues<tot_type,NQPoints,Ntot>;
  using Point = Vector<Real,Dim>;
  using QP = Matrix<Real,NQPoints,Dim>;
  using qp_points_type=typename QuadratureRule::qp_points_type;
  virtual ~BaseShapeFunctionOperatorDependent(){};
  
  virtual void
  value(const Point& point, vector_single_type& func_ )   =0;

  const FQPValues<single_type,NQPoints,Ndofs>& reference()const{return reference_func_values_;}
  const type& function()const{return func_values_;}


  template<Integer N=NComponents>
  typename std::enable_if< (1<N),const type& >::type function ()const {return func_values_;}

  template<Integer N=NComponents>
  typename std::enable_if< (1==N),const type& >::type function ()const {return component_func_values_;}

  // const type& function()const{return func_values_;}



  // compute the reference shape function for just one component
  void init_reference()
  {
   qp_points_=QuadratureRule::qp_points;
   for(Integer qp=0;qp<NQPoints;qp++)
    {
    qp_points_.get_row(qp,qp_point_);
    value(qp_point_,func_);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
        reference_func_values_[n_dof][qp]=func_[n_dof];     
     }
    }
  };



  // compute the actual shape function for just one component
  template<typename Mapping>
  void init_component(const Mapping& J)
  {
   map_.init(J);
   const auto& mapping=map_();
   for(Integer qp=0;qp<NQPoints;qp++)
    {
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
        {
          component_func_values_[n_dof][qp]= mapping * reference_func_values_[n_dof][qp];
        }
    }
  }

  template<typename Mapping>
  void init_component(const Mapping& J, const Vector<Real,Ndofs> &alpha)
  {
   map_.init(J);
   const auto& mapping=map_();
   for(Integer qp=0;qp<NQPoints;qp++)
    {
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
        {
          component_func_values_[n_dof][qp]=alpha[n_dof] * mapping * reference_func_values_[n_dof][qp];
        }
    }
  }

  template<typename Jacobian,Integer N=NComponents>
  typename std::enable_if< 1==N,void>::type
  init2(const Jacobian& J)
  {
  init_component(J);
  };

  template<typename Jacobian,Integer N=NComponents>
  typename std::enable_if< 1<N,void>::type
  init2(const Jacobian& J)
  {
  init_component(J);
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
          n_tot_=n_dof * NComponents +  n_comp ;
          n_=n_comp*ShapeFunctionDim1;
          for(Integer qp=0;qp<NQPoints;qp++)
          {             
            func_values_[n_tot_][qp].zero();
            assign(func_values_[n_tot_][qp],component_func_values_[n_dof][qp],n_,0);
          }
                 
      }
     }
  };


  template<typename Jacobian>
  void init(const Jacobian& J)
  {
  map_.init(J);
  const auto& mapping=map_();
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {

          n_tot_=n_dof * NComponents +  n_comp ;
          n_=n_comp*ShapeFunctionDim1;
          for(Integer qp=0;qp<NQPoints;qp++)
          {             
            func_values_[n_tot_][qp].zero();
            func_tmp_= mapping * reference_func_values_[n_dof][qp];
            assign(func_values_[n_tot_][qp],func_tmp_,n_,0);
          }
                 
      }
     }
  };

  template<typename Jacobian>
  void init(const Jacobian& J, const Vector<Real,Ndofs> &alpha)
  {
  map_.init(J);
  const auto& mapping=map_();
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {

          n_tot_=n_dof * NComponents +  n_comp ;
          n_=n_comp*ShapeFunctionDim1;
          for(Integer qp=0;qp<NQPoints;qp++)
          {             
            func_values_[n_tot_][qp].zero();
            func_tmp_=alpha[n_dof] * mapping * reference_func_values_[n_dof][qp];
            assign(func_values_[n_tot_][qp],func_tmp_,n_,0);
          }
                 
      }
     }
  };
  BaseShapeFunctionOperatorDependent()
  {};


  private: 
      vector_single_type func_;
      Point qp_point_;
      FQPValues<single_type,NQPoints,Ndofs> reference_func_values_;
      FQPValues<single_type,NQPoints,Ndofs> component_func_values_;
      type func_values_;
      single_type func_tmp_;
      QuadratureRule quadrature_;
      qp_points_type qp_points_;
      MapFromReference5<IdentityOperator,Elem,BaseFunctionSpace::FEFamily> map_;
      Integer n_tot_,n_;
};




template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace>
class BaseShapeFunctionOperatorDependent<Elem,BaseFunctionSpace,GradientOperator,QuadratureRule>
{ 
  public:
  using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
  static constexpr Integer Dim=Elem::Dim;
  static constexpr Integer ManifoldDim=Elem::ManifoldDim;
  static constexpr Integer NComponents=BaseFunctionSpace::NComponents;
  static constexpr Integer ShapeFunctionDim1=FunctionSpace::ShapeFunctionDim1;
  static constexpr Integer ShapeFunctionDim2=FunctionSpace::ShapeFunctionDim2;
  static constexpr Integer NQPoints=QuadratureRule::NQPoints;
  static constexpr Integer Ntot=FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,BaseFunctionSpace>>::value;
  static constexpr Integer Ndofs=Ntot/NComponents;
  using single_type   = Matrix<Real, ShapeFunctionDim1 , ShapeFunctionDim2 * Dim >;
  using vector_single_type   = Vector<single_type,Ndofs>;
  using tot_type= Matrix<Real, ShapeFunctionDim1 * NComponents, ShapeFunctionDim2 * Dim >;
  using type= FQPValues<tot_type,NQPoints,Ntot>;
  using Point = Vector<Real,Dim>;
  using QP = Matrix<Real,NQPoints,Dim>;
  using qp_points_type=typename QuadratureRule::qp_points_type;
  virtual ~BaseShapeFunctionOperatorDependent(){};
 
  virtual void value(const Point& point,    Vector<single_type,Ndofs>& func_grad_)=0;

 
  const FQPValues<single_type,NQPoints,Ndofs>& reference()const{return reference_grad_values_;}
  const type& function()const{return grad_values_;}


  void init_reference()
  {
   qp_points_=QuadratureRule::qp_points;
   for(Integer qp=0;qp<NQPoints;qp++)
    {
    qp_points_.get_row(qp,qp_point_);
    value(qp_point_,grad_);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      reference_grad_values_[n_dof][qp]=grad_[n_dof];           
     }
    }
  };



  template<typename Mapping>
  void init(const Mapping& J)
  {
   map_.init(J);
   const auto& mapping=map_();
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        n_tot_=n_dof * NComponents +  n_comp ;
        n_=n_comp * ShapeFunctionDim1;
        for(Integer qp=0;qp<NQPoints;qp++)
        {
          grad_values_[n_tot_][qp].zero();  
          grad_tmp_= contract(mapping, reference_grad_values_[n_dof][qp]);
          assign(grad_values_[n_tot_][qp],grad_tmp_,n_,0);    
        }
      }
     }
  };

  template<typename Mapping>
  void init(const Mapping& J,const Vector<Real,Ndofs> &alpha)
  {
   map_.init(J);
   const auto& mapping=map_();
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        n_tot_=n_dof * NComponents +  n_comp ;
        n_=n_comp * ShapeFunctionDim1;
        for(Integer qp=0;qp<NQPoints;qp++)
        {
          grad_values_[n_tot_][qp].zero();  
          grad_tmp_=alpha[n_dof]*contract(mapping, reference_grad_values_[n_dof][qp]);
          assign(grad_values_[n_tot_][qp],grad_tmp_,n_,0);    
        }
      }
     }
  };

  BaseShapeFunctionOperatorDependent(){};

  private: 
      Vector<single_type,Ndofs> grad_;
      Point qp_point_;
      FQPValues<single_type,NQPoints,Ndofs> reference_grad_values_;
      type grad_values_;
      single_type grad_tmp_;
      QuadratureRule quadrature_;
      qp_points_type qp_points_;
      Contraction contract;
      MapFromReference5<GradientOperator,Elem,BaseFunctionSpace::FEFamily> map_;
      Integer n_tot_,n_;
};

























template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace>
class BaseShapeFunctionOperatorDependent<Elem,BaseFunctionSpace,DivergenceOperator,QuadratureRule>
{ 
  public:
  using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
  static constexpr Integer Dim=Elem::Dim;
  static constexpr Integer ManifoldDim=Elem::ManifoldDim;
  static constexpr Integer NComponents=BaseFunctionSpace::NComponents;
  static constexpr Integer NQPoints=QuadratureRule::NQPoints;
  static constexpr Integer Ntot=FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,BaseFunctionSpace>>::value;
  static constexpr Integer Ndofs=Ntot/NComponents;
  using single_type   = Matrix<Real,1,1>;
  using vector_single_type   = Vector<single_type,Ndofs>;
  using tot_type= Matrix<Real, NComponents,1 >;
  using type= FQPValues<tot_type,NQPoints,Ntot>;
  using Point = Vector<Real,Dim>;
  using QP = Matrix<Real,NQPoints,Dim>;
  using qp_points_type=typename QuadratureRule::qp_points_type;
  virtual ~BaseShapeFunctionOperatorDependent(){};
 
  virtual void value(const Point& point,    Vector<single_type,Ndofs>& func_div_)=0;

 
  const FQPValues<single_type,NQPoints,Ndofs>& reference()const{return reference_div_values_;}
  const type& function()const {return div_values_;}


  void init_reference()
  {
   qp_points_=quadrature_.qp_points();
   for(Integer qp=0;qp<NQPoints;qp++)
    {
    qp_points_.get_row(qp,qp_point_);
    value(qp_point_,div_);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      reference_div_values_[n_dof][qp]=div_[n_dof];           
     }
    }
  };

  template<typename Jacobian>
  void init(const Jacobian& J)
  {
   map_.init(J);
   const auto& mapping=map_();
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     { 
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        n_tot_=n_dof * NComponents +  n_comp ;
        n_=n_comp;
        for(Integer qp=0;qp<NQPoints;qp++)
        {
          div_values_[n_tot_][qp].zero();  
          div_tmp_= contract(mapping, reference_div_values_[n_dof][qp]);           
          assign(div_values_[n_tot_][qp],div_tmp_,n_,0);   
        }
      }
     }
  };

  template<typename Jacobian>
  void init(const Jacobian& J,const Vector<Real,Ndofs> &alpha)
  {
   map_.init(J);
   const auto& mapping=map_();
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     { 
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        n_tot_=n_dof * NComponents +  n_comp ;
        n_=n_comp;
        for(Integer qp=0;qp<NQPoints;qp++)
        {
          div_values_[n_tot_][qp].zero();  
          div_tmp_=alpha[n_dof]*contract(mapping, reference_div_values_[n_dof][qp]);           
          assign(div_values_[n_tot_][qp],div_tmp_,n_,0);   
        }
      }
     }
  };

  BaseShapeFunctionOperatorDependent(){};

  private: 
      Vector<single_type,Ndofs> div_;
      Point qp_point_;
      FQPValues<single_type,NQPoints,Ndofs> reference_div_values_;
      type div_values_;
      single_type div_tmp_;
      QuadratureRule quadrature_;
      qp_points_type qp_points_;
      Contraction contract;
      MapFromReference5<DivergenceOperator,Elem,BaseFunctionSpace::FEFamily> map_;
      Integer n_tot_,n_;
};
















template<typename Elem,Integer FEFamily,Integer Order, typename ConstInput, typename ShapeFunctionCoefficient>
void shape_function_coefficients_init(const ConstInput& mesh_ptr,ShapeFunctionCoefficient& coeff);


  template<typename Elem,typename Operator,Integer FEFamily,Integer Order,typename Output,typename Point>
constexpr void value(const Point& point,Output& output);

  template<> 
    constexpr void value<Simplex<2,2>, IdentityOperator, LagrangeFE, 1>
    (const Vector<Real,2>& point, Vector<Matrix<Real, 1, 1>,3> & func)
    {
       Vector<Matrix<Real, 1, 1>,3> func2((1. - point[0] - point[1]), // 1 in (0,0)
                                           point[0],                  // 1 in (1,0)
                                           point[1]);                 // 1 in (0,1)
       func=func2; 
     }

  template<> 
     constexpr void value<Simplex<2,2>, GradientOperator, LagrangeFE, 1>
     (const Vector<Real,2>& point, Vector<Vector<Real, 2>,3> & func)
     {
      const auto& xi=point[0];
      const auto& eta=point[1];
      const Real zeta = 1. - xi - eta;
      Vector<Vector<Real, 2>,3> func2({-1,-1},
       {+1, 0},
       { 0,+1});     
      func=func2;
    }

    //shape_function_coefficients_init for: Simplex<2,2>, LagrangeFE, 1
    // do nothing: lagrange shape functions do not need any coefficient


  template<> 
    constexpr void value<Simplex<2,2>, IdentityOperator, LagrangeFE, 2>
    (const Vector<Real,2>& point, Vector<Matrix<Real, 1, 1>,6> & func)
    {
      const auto& xi=point[0];
      const auto& eta=point[1];
      const Real zeta = 1. - xi - eta;
        Vector<Matrix<Real, 1, 1>,6> func2(2.*zeta*(zeta-0.5), // 1 in (0,0)
                                           2.*xi*(xi-0.5),     // 1 in (1,0)
                                           2.*eta*(eta-0.5),   // 1 in (0,1)
                                           4.*zeta*xi,         // 1 in (0.5,0)
                                           4.*eta*zeta,        // 1 in (0,0.5)
                                           4.*xi*eta);         // 1 in (0.5,0.5)
        func=func2;
      }




  template<> 
      constexpr void value<Simplex<2,2>, GradientOperator, LagrangeFE, 2>
      (const Vector<Real,2>& point, Vector<Vector<Real, 2>,6> & func)
      {
        const auto& xi=point[0];
        const auto& eta=point[1];
        const Real zeta = 1. - xi - eta;
        const Real dxi_dxi    = 1.;
        const Real deta_deta  = 1.;
        const Real dzeta_dxi  = -1.;
        const Real dzeta_deta = -1.;
        Vector<Vector<Real, 2>,6> func2(
          {2.*zeta*dzeta_dxi  + 2*dzeta_dxi *(zeta-0.5), 2.*zeta*dzeta_deta + 2*dzeta_deta*(zeta-0.5)}, 
          {2.*xi*dxi_dxi  + 2.*dxi_dxi *(xi-0.5),0},
          {0,2.*eta*deta_deta + 2.*deta_deta*(eta-0.5)}, 
          {4.*zeta*dxi_dxi+4.*dzeta_dxi*xi,4.*dzeta_deta*xi},
          {4.*dxi_dxi*eta,4.*xi*deta_deta},
          {4.*eta*dzeta_dxi,4.*eta*dzeta_deta + 4.*deta_deta*zeta}); 

        func=func2;
      }





    //shape_function_coefficients_init for: Simplex<2,2>, LagrangeFE, 2
    // do nothing: lagrange shape functions do not need any coefficient

  template<> 
      constexpr void value<Simplex<2,2>, IdentityOperator, RaviartThomasFE, 0>
      (const Vector<Real,2>& point, Vector<Matrix<Real, 2, 1>,3> & func)
      {
       const auto& xi=point[0];
       const auto& eta=point[1];
       Vector<Matrix<Real, 2, 1>,3> func2{{xi,eta-1},{xi-1,eta},{xi,eta}};
       func=func2;
     }

  template<> 
     constexpr void value<Simplex<2,2>, DivergenceOperator, RaviartThomasFE, 0>
     (const Vector<Real,2>& point, Vector<Matrix<Real, 1, 1>,3> & func)
     {
       Vector<Matrix<Real, 1, 1>,3> func2{{2},{2},{2}};
       func=func2;
     } 


  template<>
    void shape_function_coefficients_init<Simplex<2,2>, RaviartThomasFE, 0>
    (const Vector<Real, 3 >& outward,Vector<Real, 3 >& coeff)
    {
     coeff[0]=outward[0];
     coeff[1]=outward[1];
     coeff[2]=outward[2];
    }

  template<> 
     constexpr void value<Simplex<2,2>, IdentityOperator, RaviartThomasFE, 1>
     (const Vector<Real,2>& point, Vector<Matrix<Real, 2, 1>,8> & func)
     {
       const auto& xi=point[0];
       const auto& eta=point[1];

       Vector<Matrix<Real, 2, 1>,8> func2{
        {(1. - xi - eta)*xi,(1. - xi - eta)*(eta-1)},   // 0 in (1,0), (0,1), non-zero normal on edge0
        {xi*xi,xi*(eta-1)},                             // 0 in (0,0), (0,1), non-zero normal on edge0
        {(1. - xi - eta)*(xi-1),(1. - xi - eta)*(eta)}, // 0 in (1,0), (0,1), non-zero normal on edge1
        {eta*(xi-1),eta*eta},                           // 0 in (0,0), (1,0), non-zero normal on edge1
        {xi*xi,xi*eta},                                 // 0 in (0,0), (0,1), non-zero normal on edge2
        {eta*xi,eta*eta},                               // 0 in (0,0), (1,0), non-zero normal on edge2
        {eta*xi,eta*(eta-1)},                           // normal 0 on all edges, element-dof
        {xi*(xi-1),xi*eta}                              // normal 0 on all edges, element-dof
      };
      func=func2;
    }


  template<> 
    constexpr void value<Simplex<2,2>, DivergenceOperator, RaviartThomasFE, 1>
    (const Vector<Real,2>& point, Vector<Matrix<Real, 1, 1>,8> & func)
    {
     const auto& xi=point[0];
     const auto& eta=point[1];
     Vector<Matrix<Real, 1, 1>,8> func2{
      {3*(1-xi-eta)},
      {3*xi},
      {3*(1-xi-eta)},
      {3*eta},
      {3*xi},
      {3*eta},
      {3*eta},
      {3*xi}
    };
    func=func2;
  } 

  template<>
    void shape_function_coefficients_init<Simplex<2,2>, RaviartThomasFE, 1>
    (const Vector<Real, 3 >& outward,Vector<Real, 8 >& coeff)
    {
     coeff[0]=outward[0];
     coeff[1]=outward[0];
     coeff[2]=outward[1];
     coeff[3]=outward[1];
     coeff[4]=outward[2];
     coeff[5]=outward[2];
     coeff[6]=outward[0];
     coeff[7]=outward[1];
    }

  template<typename Elem,typename Operator, Integer FEFamily,Integer Order,typename single_type,Integer Ndofs,
  Integer NQPoints,Integer Dim>
  constexpr const Vector<Vector<single_type,NQPoints>,Ndofs> reference_shape_function_init(const Matrix<Real,NQPoints,Dim>&qp_points)
  {
    Vector<Vector<single_type,NQPoints>,Ndofs> v;
    Vector<Real,Dim> qp_point;
    Vector<single_type,Ndofs> func;
    for(Integer qp=0;qp<NQPoints;qp++)
    {
     qp_point=qp_points.get_row(qp);
             // func=value<Elem,Operator,FEFamily,Order,single_type,Ndofs>(qp_point);
     value<Elem,Operator,FEFamily,Order>(qp_point,func);
     for(Integer n_dof = 0; n_dof < Ndofs; ++n_dof) {
      const_cast<single_type&>
      (static_cast<const std::array<single_type,NQPoints>& >
       ((static_cast<const std::array<Vector<single_type,NQPoints>,Ndofs>& >(v())[n_dof] )())[qp])=
      static_cast<const std::array<single_type,Ndofs>& >(func())[n_dof];
    }
  }
  return v;
 };


  template<typename T,Integer NQPoints,Integer Ndofs>
  constexpr FQPValues<T,NQPoints,Ndofs> 
  weighted_reference_shape_function_init(const FQPValues<T,NQPoints,Ndofs>&input,const Array<Real,NQPoints>& weights)
  {
    FQPValues<T,NQPoints,Ndofs> output;

    for(Integer qp=0;qp<NQPoints;qp++)
    {
     for(Integer n_dof = 0; n_dof < Ndofs; ++n_dof) {


      const_cast<T&>
      (static_cast<const Vector<T,NQPoints>& >
       ((static_cast<const FQPValues<T,NQPoints,Ndofs>& >(output)[n_dof] ))[qp])
      // output[n_dof][qp]
      =
      static_cast<const Array<Real,NQPoints>& >(weights)[qp]*
      const_cast<T&>
      (static_cast<const Vector<T,NQPoints>& >
       ((static_cast<const FQPValues<T,NQPoints,Ndofs>& >(input)[n_dof] ))[qp]);
    }
  }
  return output;
 };




template<typename FunctionSpace,typename Operator>
class SingleTypeShapeFunction;



template<typename FunctionSpace>
class SingleTypeShapeFunction<FunctionSpace,IdentityOperator>
{
public:
  static constexpr Integer NComponents=FunctionSpace::NComponents;
  static constexpr Integer ShapeFunctionDim1=FunctionSpace::ShapeFunctionDim1;
  static constexpr Integer ShapeFunctionDim2=FunctionSpace::ShapeFunctionDim2;  
  using SingleType=Matrix<Real, ShapeFunctionDim1, ShapeFunctionDim2>;
  using TotType= Matrix<Real, ShapeFunctionDim1 * NComponents, ShapeFunctionDim2>;
};
template<typename FunctionSpace>
class SingleTypeShapeFunction<FunctionSpace,GradientOperator>
{
public:
  static constexpr Integer NComponents=FunctionSpace::NComponents;
  static constexpr Integer ShapeFunctionDim1=FunctionSpace::ShapeFunctionDim1;
  static constexpr Integer ShapeFunctionDim2=FunctionSpace::ShapeFunctionDim2; 
  static constexpr Integer Dim=FunctionSpace::Elem::Dim;
  using SingleType= typename
  std::conditional_t<(1==ShapeFunctionDim1 && 1==ShapeFunctionDim2),
                      Vector<Real, Dim >,
                      Matrix<Real, ShapeFunctionDim1 , ShapeFunctionDim2 * Dim>
                  >;
  using TotType= Matrix<Real, ShapeFunctionDim1 * NComponents, ShapeFunctionDim2 * Dim >;
};


template<typename FunctionSpace>
class SingleTypeShapeFunction<FunctionSpace,DivergenceOperator>
{
public:
  static constexpr Integer NComponents=FunctionSpace::NComponents;
  static constexpr Integer ShapeFunctionDim1=1;
  static constexpr Integer ShapeFunctionDim2=1; 
  using SingleType=Matrix<Real,1,1>;
  using TotType= Matrix<Real, NComponents,1 >;   
};










template< typename Elem_,typename BaseFunctionSpace, typename Operator_, typename QuadratureRule>
class ShapeFunctionDependent 
{ 
  public:
  using Elem=Elem_;
  using Operator=Operator_;

  using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
  static constexpr Integer Dim=Elem::Dim;
  static constexpr Integer ManifoldDim=Elem::ManifoldDim;
  static constexpr Integer NComponents=BaseFunctionSpace::NComponents;
  static constexpr Integer NQPoints=QuadratureRule::NQPoints;
  static constexpr Integer Ntot=FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,BaseFunctionSpace>>::value;
  static constexpr Integer Ndofs=Ntot/NComponents;
  static constexpr Integer Order=BaseFunctionSpace::Order;
  static constexpr Integer FEFamily=BaseFunctionSpace::FEFamily;
  
  using SingleType   = typename SingleTypeShapeFunction<FunctionSpace,Operator>::SingleType;
  using VectorSingleType   = Vector<SingleType,Ndofs>;
  using tot_type= typename SingleTypeShapeFunction<FunctionSpace,Operator>::TotType;
  using type= FQPValues<tot_type,NQPoints,Ntot>;
  using Point = Vector<Real,Dim>;
  using QP = Matrix<Real,NQPoints,Dim>;
  using qp_points_type=typename QuadratureRule::qp_points_type;
  using Map=MapFromReference5<Operator,Elem,BaseFunctionSpace::FEFamily>;

  static constexpr Integer ShapeFunctionDim1=SingleTypeShapeFunction<FunctionSpace,Operator>::ShapeFunctionDim1;
  static constexpr Integer ShapeFunctionDim2=SingleTypeShapeFunction<FunctionSpace,Operator>::ShapeFunctionDim2;

  static constexpr FQPValues<SingleType,NQPoints,Ndofs>  
  reference_values{reference_shape_function_init<Elem,Operator,FEFamily,Order,SingleType,Ndofs>(QuadratureRule::qp_points)};
  static constexpr FQPValues<SingleType,NQPoints,Ndofs>
  weighted_reference_values{  weighted_reference_shape_function_init(reference_values,QuadratureRule::qp_sqrt_abs_weights)};

  constexpr const type& eval()const{return func_values_;}


   void init()
    {
    const auto& map=(*map_ptr);
    const auto& mapping=map();
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
       {
        for(Integer n_comp=0;n_comp<NComponents;n_comp++)
        {

            n_tot_=n_dof * NComponents +  n_comp ;
            n_=n_comp*ShapeFunctionDim1;
            for(Integer qp=0;qp<NQPoints;qp++)
            {             
                func_values_[n_tot_][qp].zero();
              func_tmp_=  mapping * weighted_reference_values[n_dof][qp];
              assign(func_values_[n_tot_][qp],func_tmp_,n_,0);
             }
                   
        }
       }
    }

   void init(const Vector<Real,Ndofs> &alpha)
    {
    const auto& map=(*map_ptr);
    const auto& mapping=map();
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
       {
        for(Integer n_comp=0;n_comp<NComponents;n_comp++)
        {

            n_tot_=n_dof * NComponents +  n_comp ;
            n_=n_comp*ShapeFunctionDim1;
            for(Integer qp=0;qp<NQPoints;qp++)
            {             
              func_values_[n_tot_][qp].zero();
              func_tmp_=alpha[n_dof] * mapping * weighted_reference_values[n_dof][qp];             
              assign(func_values_[n_tot_][qp],func_tmp_,n_,0);
         }
                   
        }
       }
    };

    constexpr void init_map(const Map& map){map_ptr=std::make_shared<Map>(map);}

    ShapeFunctionDependent(const Map& map):
    map_ptr(std::make_shared<Map>(map))
  {}

  ShapeFunctionDependent(){}

 const auto& map()const{return (*map_ptr);}

  private: 
      SingleType func_tmp_;
      VectorSingleType func_;
      Point qp_point_;
      FQPValues<SingleType,NQPoints,Ndofs> component_func_values_;
      type func_values_;
      std::shared_ptr<Map> map_ptr;
      Integer n_tot_;
      Integer n_;  
};

 






template<typename Elem,typename BaseFunctionSpace,typename Operator, typename QuadratureRule>
constexpr FQPValues<typename ShapeFunctionDependent<Elem,BaseFunctionSpace,Operator,QuadratureRule>::SingleType,
                    ShapeFunctionDependent<Elem,BaseFunctionSpace,Operator,QuadratureRule>::NQPoints,
                    ShapeFunctionDependent<Elem,BaseFunctionSpace,Operator,QuadratureRule>::Ndofs> 
ShapeFunctionDependent<Elem,BaseFunctionSpace,Operator,QuadratureRule>::reference_values;

template<typename Elem,typename BaseFunctionSpace,typename Operator, typename QuadratureRule>
constexpr FQPValues<typename ShapeFunctionDependent<Elem,BaseFunctionSpace,Operator,QuadratureRule>::SingleType,
                    ShapeFunctionDependent<Elem,BaseFunctionSpace,Operator,QuadratureRule>::NQPoints,
                    ShapeFunctionDependent<Elem,BaseFunctionSpace,Operator,QuadratureRule>::Ndofs> 
ShapeFunctionDependent<Elem,BaseFunctionSpace,Operator,QuadratureRule>::weighted_reference_values;





}



#endif






