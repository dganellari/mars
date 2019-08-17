




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
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
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
            outward_[ee][mm]=sgn(dot(diff,n)); 
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






























template<Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents>
class BaseShapeFunctionOperator2
{public:

  virtual ~BaseShapeFunctionOperator2(){};

  virtual const Matrix<Real,Ndofs,ShapeFunctionDim>
  phi_single_component(const IdentityOperator& o, const Vector<Real,Dim>& point )=0;

  virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
  phi_single_component(const GradientOperator&o, const Vector<Real,Dim>& point)=0;

  virtual const Vector<Real,Ndofs>
  phi_single_component(const DivergenceOperator&o, const Vector<Real,Dim>& qp_point )=0;



  template<Integer NQPoints>
  Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPoints>,Ndofs>
  operator()(const IdentityOperator&o,const Matrix<Real,NQPoints,Dim>& qp_points)
  {
   Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPoints>,Ndofs> f;
   //Vector<Real,ShapeFunctionDim> row_tmp;
   Vector<Real,Dim> qp_point;

   for(Integer qp=0;qp<NQPoints;qp++)
    {
    qp_points.get_row(qp,qp_point);
    auto phi_single=phi_single_component(o,qp_point);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer ii=0;ii<ShapeFunctionDim;ii++)
        f[n_dof][qp][ii]=phi_single(n_dof,ii);     
     }
    }
   std::cout<<"identity operator"<<std::endl;
   std::cout<<f<<std::endl;
   
   return f;
  };

  template<Integer NQPoints, typename Mapping>
  Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPoints>,NComponents>,Ndofs>
  operator()(const IdentityOperator&o,
       const Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPoints>,Ndofs>& reference_phi,
       const Mapping& mapping, 
       const Vector<Real,Ndofs> &alpha)
  {
   Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPoints>,NComponents>,Ndofs> result;
   Vector<Real,ShapeFunctionDim> row;
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        for(Integer qp=0;qp<NQPoints;qp++)
        { 
          result[n_dof][n_comp][qp].zero();
          row=alpha[n_dof] * mapping*reference_phi[n_dof][qp];
          result[n_dof][n_comp][qp].row(n_comp,row);
        }        
      }
     }
   return result;
  };

  template<Integer NQPoints>
  Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs>
  operator()(const GradientOperator& o, const Matrix<Real,NQPoints,Dim>& qp_points)
  {
   Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs> f;
   Vector<Real,Dim> qp_point;

   for(Integer qp=0;qp<NQPoints;qp++)
    {
    qp_points.get_row(qp,qp_point);
    auto grad_phi_single=phi_single_component(o,qp_point);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      auto grad_phi=grad_phi_single[n_dof];
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        f[n_dof][n_comp][qp].zero();
        for(Integer ii=0;ii<ShapeFunctionDim;ii++)
          {
            for(Integer jj=0;jj<Dim;jj++)
            {
              f[n_dof][n_comp][qp](n_comp*ShapeFunctionDim+ii,jj)=grad_phi(ii,jj);
            }
          }
      }
     }
    }

   return f;
  };


  template<Integer NQPoints, typename Mapping>
  Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs>
  operator()(const GradientOperator& o,
       const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs>& reference_grad_phi,
       const Mapping& mapping, 
       const Vector<Real,Ndofs> &alpha)
  {
   Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs> result;
   Vector<Real,Dim> row;
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        for(Integer qp=0;qp<NQPoints;qp++)
        {
          result[n_dof][n_comp][qp].zero();
          reference_grad_phi[n_dof][n_comp][qp].get_row(n_comp,row); 
          // std::cout<<" firs row ="<<std::endl;
          // reference_grad_phi[n_dof][n_comp][qp].describe(std::cout);     
          // row.describe(std::cout);     
          // std::cout<<" mapping ="<<std::endl;
          // mapping.describe(std::cout);    
          row=alpha[n_dof]*mapping*row;      
          result[n_dof][n_comp][qp].row(n_comp,row);
          // std::cout<<" result ="<<std::endl;
          // result[n_dof][n_comp][qp].describe(std::cout);        
        }
      }
     }
   return result;
  };

  template<Integer NQPoints>
  Vector<Vector< Real, NQPoints>,Ndofs>
  operator()(const DivergenceOperator&o,const Matrix<Real,NQPoints,Dim>& qp_points)
  {
        Vector<Vector< Real, NQPoints>,Ndofs> div;
        Vector<Real,Dim> qp_point;

        for(Integer qp=0;qp<NQPoints;qp++)
        {
          qp_points.get_row(qp,qp_point);
          auto div_phi_single=phi_single_component(o,qp_point);
              for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
                   div[n_dof][qp]=div_phi_single[n_dof];
        }
        return div;};

  template<Integer NQPoints,typename Mapping>
  Vector<Vector<Vector< Matrix<Real,NComponents,1>, NQPoints>,NComponents>,Ndofs>
  operator()(const DivergenceOperator&o,
       const Vector<Vector< Real, NQPoints>,Ndofs>& divphi_reference,
       const Mapping& mapping,
       const Vector<Real,Ndofs> &alpha)
  {
        Vector<Vector<Vector< Matrix<Real,NComponents,1>, NQPoints>,NComponents>,Ndofs> divphi;
        for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
         for(Integer n_comp=0;n_comp<NComponents;n_comp++)
          for(Integer qp=0;qp<NQPoints;qp++)
            { divphi[n_dof][n_comp][qp].zero();
              divphi[n_dof][n_comp][qp](n_comp,0)=alpha[n_dof] * mapping * divphi_reference[n_dof][qp];
            } 
    return divphi;};

};





template<typename Elem,typename BaseFunctionSpace>
class ShapeFunctionOperator2;



template<Integer NComponents_>
class ShapeFunctionOperator2<Simplex<2,2>, Lagrange1<NComponents_> > : 
public BaseShapeFunctionOperator2<3,2,2,1,NComponents_>
{
public:
  virtual const Matrix<Real,3,1>
  phi_single_component(const IdentityOperator& o,const Vector<Real,2>& point)
       {const auto& xi=point[0];
        const auto& eta=point[1];
        const Real zeta = 1. - xi - eta;
        const Matrix<Real,3,1> shapefunction{zeta,xi,eta};
        return shapefunction;};

  virtual const Vector< Matrix<Real,1,2>, 3>
  phi_single_component(const GradientOperator& o, const Vector<Real,2>& point)
  {
    Vector< Matrix<Real,1,2>, 3> grad{{-1,-1},{1,0},{0,1}};
    return grad;
  }

  virtual const Vector<Real,3>
  phi_single_component(const DivergenceOperator&o, const Vector<Real,2>& qp_point )
  {Vector<Real,3> div; return div;};


};


template<Integer NComponents_>
class ShapeFunctionOperator2<Simplex<2,2>, Lagrange2<NComponents_> > : 
public BaseShapeFunctionOperator2<6,2,2,1,NComponents_>
{
public:
      virtual const Matrix<Real,6,1>
      phi_single_component(const IdentityOperator& o, const Vector<Real,2>& point) 
      {
          const auto& xi=point[0];
          const auto& eta=point[1];
          const Real zeta = 1. - xi - eta;
          Matrix<Real,6,1> shape_function{2.*zeta*(zeta-0.5),
                                        2.*xi*(xi-0.5),
                                        2.*eta*(eta-0.5),
                                        4.*zeta*xi,
                                        4.*xi*eta, 
                                        4.*eta*zeta };
          return shape_function;
      };

  virtual const Vector< Matrix<Real,1,2>, 6>
  phi_single_component(const GradientOperator& o, const Vector<Real,2>& point)
  {
   const auto& xi=point[0];
   const auto& eta=point[1];
   const Real zeta = 1. - xi - eta;
   Vector< Matrix<Real,1,2>, 6> grad{
                                     { - 4 * (1 - xi - eta) + 1, - 4 * (1 - xi - eta) + 1},
                                     { 4 * xi - 1              , 0                       },
                                     { 0                       , 4 * eta -1              },
                                     { 4 - 8 * xi - 4 * eta    , - 4 * xi                },
                                     { 4 * eta                 , 4 * xi                  },
                                     { -4 * eta                , 4 - 4 * xi - 8 * eta    } };
  return grad;
  }


  virtual const Vector<Real,6>
  phi_single_component(const DivergenceOperator&o, const Vector<Real,2>& qp_point )
  {Vector<Real,6> div; return div;};



};







template<Integer NComponents_>
class ShapeFunctionOperator2<Simplex<2,2>, RT0<NComponents_> > : 
public BaseShapeFunctionOperator2<3,2,2,2,NComponents_>
{
public:
  virtual const  Matrix<Real,3,2>
  phi_single_component(const IdentityOperator&o, const Vector<Real,2>& point)
       {const auto& xi=point[0];
        const auto& eta=point[1];
        const Matrix<Real,3,2> shapefunction{xi,eta-1,xi-1,eta,xi,eta};
        return shapefunction;};

  virtual const Vector< Matrix<Real,2,2>, 3>
  phi_single_component(const GradientOperator&o, const Vector<Real,2>& point)
  {Vector< Matrix<Real,2,2>, 3> grad; return grad;};

  virtual const Vector<Real,3>
  phi_single_component(const DivergenceOperator& o, const Vector<Real,2>& qp_point )
  {
   Vector<Real,3> tryme;
     tryme[0]=2;
     tryme[1]=2;
     tryme[2]=2;    
   return tryme;
  };

};




































































template<typename QuadratureRule, Integer Ndofs, Integer Dim, Integer ManifoldDim, 
         Integer ShapeFunctionDim1=1, Integer ShapeFunctionDim2=1, Integer NComponents=1>
class BaseShapeFunctionOperator4
{ 
  public:
  static constexpr Integer NQPoints=QuadratureRule::NQPoints;
  static constexpr Integer Ntot = Ndofs * NComponents ;
  using FuncType   = Matrix<Real, ShapeFunctionDim1, ShapeFunctionDim2>;
  using TotFuncType= Matrix<Real, ShapeFunctionDim1 * NComponents, ShapeFunctionDim2>;
  using GradType   = Matrix<Real, ShapeFunctionDim1 , ShapeFunctionDim2 * Dim >;
  using TotGradType= Matrix<Real, ShapeFunctionDim1 * NComponents, ShapeFunctionDim2 * Dim >;
  using FunctionType= FQPValues<TotFuncType,NQPoints,Ntot>;
  using GradientType= FQPValues<TotGradType,NQPoints,Ntot>;
  using type= FunctionType;

  using DivType = Real;
  using TotDivType = Matrix<Real,NComponents,1>;
  using Point = Vector<Real,Dim>;
  using QP = Matrix<Real,NQPoints,Dim>;
  

  virtual ~BaseShapeFunctionOperator4(){};
  
  virtual void
  value(const IdentityOperator& o,  const Point& point,    Vector<FuncType,Ndofs>& func_ )   =0;

  virtual void
  value(const GradientOperator&o,   const Point& point,    Vector<GradType,Ndofs>& func_grad_)=0;

  virtual void
  value(const DivergenceOperator&o, const Point& qp_point, Vector<DivType,Ndofs>& func_div_)  =0;

  FQPValues<FuncType,NQPoints,Ndofs>& reference(const IdentityOperator&o){return reference_func_values_;}
  FQPValues<GradType,NQPoints,Ndofs>& reference(const GradientOperator&o){return reference_grad_values_;}
  FunctionType & function(const IdentityOperator&o){return func_values_;}
  GradientType& function(const GradientOperator&o){return grad_values_;}


  void operator()(const QP & qp_points,const IdentityOperator&o)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
    {
    qp_points.get_row(qp,qp_point_);
    value(o,qp_point_,func_);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
        reference_func_values_[n_dof][qp]=func_[n_dof];     
     }
    }
  };


  template<typename Mapping>
  void operator()(const IdentityOperator&o,
                  const Mapping& mapping, 
                  const Vector<Real,Ndofs> &alpha)
  {
  Integer n_tot,n1;
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {

          n_tot=n_dof * NComponents +  n_comp ;
          n1=n_comp*ShapeFunctionDim1;
          for(Integer qp=0;qp<NQPoints;qp++)
          {             
            func_values_[n_tot][qp].zero();
            func_tmp_=alpha[n_dof] * mapping * reference_func_values_[n_dof][qp];
            assign(func_values_[n_tot][qp],func_tmp_,n1,0);
          }
                 
      }
     }
  };

  void operator()(const QP & qp_points,const GradientOperator& o)
  {

   for(Integer qp=0;qp<NQPoints;qp++)
    {
    qp_points.get_row(qp,qp_point_);
    value(o,qp_point_,grad_);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      reference_grad_values_[n_dof][qp]=grad_[n_dof];           
     }
    }
  };

  template<typename Mapping>
  void operator()(const GradientOperator& o,
                  const Mapping& mapping, 
                  const Vector<Real,Ndofs> &alpha)
  {
    Integer n_tot,n1;
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        n_tot=n_dof * NComponents +  n_comp ;
        n1=n_comp * ShapeFunctionDim1;
        for(Integer qp=0;qp<NQPoints;qp++)
        {
          grad_values_[n_tot][qp].zero();  
          grad_tmp_=alpha[n_dof]*contract(mapping, reference_grad_values_[n_dof][qp]);
          assign(grad_values_[n_tot][qp],grad_tmp_,n1,0);    
        }
      }
     }
  };

  template<typename Operator1, typename Operator2, typename...Operators>
  void operator()(const QP & qp_points,const Operator1&o1, const Operator2& o2, const Operators&...operators)
  {
    (*this)(qp_points,o1);
    (*this)(qp_points,o2,operators...);
  };

  void operator()(const QP & qp_points)
  {
    (*this)(qp_points,IdentityOperator(),GradientOperator());
  };

  template<typename Operator1, typename Mapping1, typename...OperatorsAndMappings>
  void operator()(const Operator1&o1, const Mapping1& m1, const OperatorsAndMappings&...oandm)//,const Vector<Real,Ndofs> &alpha=1.0 )
  {
    (*this)(o1,m1);//,alpha);
    (*this)(oandm...);//,alpha);
  };

  BaseShapeFunctionOperator4(){};

  private: 
  Vector<FuncType,Ndofs> func_;
  Vector<GradType,Ndofs> grad_;
  Vector<DivType,Ndofs>  div_;
  Point qp_point_;

  FQPValues<FuncType,NQPoints,Ndofs> reference_func_values_;
  FQPValues<GradType,NQPoints,Ndofs> reference_grad_values_;

  FunctionType func_values_;
  GradientType grad_values_;
  FuncType func_tmp_;
  GradType grad_tmp_;
  Contraction contract;

};




template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace>
class ShapeFunctionOperator4;



template<typename QuadratureRule, Integer NComponents>
class ShapeFunctionOperator4<QuadratureRule,Simplex<2,2>, Lagrange1<NComponents> > : 
public BaseShapeFunctionOperator4<QuadratureRule,3,2,2,1,1,NComponents>
{
public:
  static constexpr Integer Ndofs=3;
  static constexpr Integer Dim=2;
  static constexpr Integer ManifoldDim=2;
  static constexpr Integer ShapeFunctionDim1=1;
  static constexpr Integer ShapeFunctionDim2=1;
  static constexpr Integer Ntot = Ndofs * NComponents ;
  static constexpr Integer NQPoints=QuadratureRule::NQPoints;
  using FuncType   = Matrix<Real, ShapeFunctionDim1, ShapeFunctionDim2>;
  using GradType   = Matrix<Real, ShapeFunctionDim1, ShapeFunctionDim2 * Dim>;
  using DivType = Real;
  using Point = Vector<Real,Dim>;
  using QP = Matrix<Real,NQPoints,Dim>;

  ShapeFunctionOperator4(){};
  ShapeFunctionOperator4(const QP & qp_points){(*this)(qp_points);};
 virtual void value(const IdentityOperator& o, const Point& point, Vector<FuncType,Ndofs>& func )
       {func[0](0,0)=point[0];
        func[1](0,0)=point[1];
        func[2](0,0)=1. - point[0] - point[1];};
  ;

  virtual void value(const GradientOperator&o, const Point& point, Vector<GradType,Ndofs>& func_grad)
   {
    func_grad[0](0,0)=-1;  func_grad[0](0,1)=-1; 
    func_grad[1](0,0)=+1;  func_grad[1](0,1)= 0; 
    func_grad[2](0,0)= 0;  func_grad[2](0,1)=+1; 
  } 

  virtual void value(const DivergenceOperator&o, const Point& qp_point, Vector<DivType,Ndofs>& func_div_)
  {
    std::cout<<" Lagrange1 divergence not implemented"<<std::endl;
  }


};



















template<typename Operator, typename QuadratureRule, Integer Ndofs, Integer Dim, Integer ManifoldDim, 
         Integer ShapeFunctionDim1=1, Integer ShapeFunctionDim2=1, Integer NComponents=1>
class BaseShapeFunctionOperator;





template<typename QuadratureRule, Integer Ndofs, Integer Dim, Integer ManifoldDim, 
         Integer ShapeFunctionDim1, Integer ShapeFunctionDim2, Integer NComponents>
class BaseShapeFunctionOperator<IdentityOperator, QuadratureRule, Ndofs, Dim, ManifoldDim, 
                                ShapeFunctionDim1, ShapeFunctionDim2, NComponents>

{ 
  public:
  static constexpr Integer NQPoints=QuadratureRule::NQPoints;
  static constexpr Integer Ntot = Ndofs * NComponents ;
  using FuncType   = Matrix<Real, ShapeFunctionDim1, ShapeFunctionDim2>;
  using TotFuncType= Matrix<Real, ShapeFunctionDim1 * NComponents, ShapeFunctionDim2>;
  using FunctionType= FQPValues<TotFuncType,NQPoints,Ntot>;
  using type= FunctionType;


  using Point = Vector<Real,Dim>;
  using QP = Matrix<Real,NQPoints,Dim>;
  

  virtual ~BaseShapeFunctionOperator(){};
  
  virtual void
  value(const IdentityOperator& o,  const Point& point,    Vector<FuncType,Ndofs>& func_ )   =0;

  FQPValues<FuncType,NQPoints,Ndofs>& reference(const IdentityOperator&o){return reference_func_values_;}
  FunctionType & function(const IdentityOperator&o){return func_values_;}


  void init(const QP & qp_points,const IdentityOperator&o)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
    {
    qp_points.get_row(qp,qp_point_);
    value(o,qp_point_,func_);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
        reference_func_values_[n_dof][qp]=func_[n_dof];     
     }
    }
  };


  template<typename Mapping>
  void operator()(const IdentityOperator&o,
                  const Mapping& mapping, 
                  const Vector<Real,Ndofs> &alpha)
  {
  Integer n_tot,n1;
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {

          n_tot=n_dof * NComponents +  n_comp ;
          n1=n_comp*ShapeFunctionDim1;
          for(Integer qp=0;qp<NQPoints;qp++)
          {             
            func_values_[n_tot][qp].zero();
            func_tmp_=alpha[n_dof] * mapping * reference_func_values_[n_dof][qp];
            assign(func_values_[n_tot][qp],func_tmp_,n1,0);
          }
                 
      }
     }
  };


  BaseShapeFunctionOperator(){};

  private: 
  Vector<FuncType,Ndofs> func_;
  Point qp_point_;

  FQPValues<FuncType,NQPoints,Ndofs> reference_func_values_;

  FunctionType func_values_;
  FuncType func_tmp_;
  Contraction contract;

};










template<typename QuadratureRule, Integer Ndofs, Integer Dim, Integer ManifoldDim, 
         Integer ShapeFunctionDim1, Integer ShapeFunctionDim2, Integer NComponents>
class BaseShapeFunctionOperator<GradientOperator, QuadratureRule, Ndofs, Dim, ManifoldDim, 
                                ShapeFunctionDim1, ShapeFunctionDim2, NComponents>
{ 
  public:
  static constexpr Integer NQPoints=QuadratureRule::NQPoints;
  static constexpr Integer Ntot = Ndofs * NComponents ;
  using GradType   = Matrix<Real, ShapeFunctionDim1 , ShapeFunctionDim2 * Dim >;
  using TotGradType= Matrix<Real, ShapeFunctionDim1 * NComponents, ShapeFunctionDim2 * Dim >;
  using GradientType= FQPValues<TotGradType,NQPoints,Ntot>;
  using type= GradientType;

  using Point = Vector<Real,Dim>;
  using QP = Matrix<Real,NQPoints,Dim>;
  

  virtual ~BaseShapeFunctionOperator(){};
 
  virtual void
  value(const GradientOperator&o,   const Point& point,    Vector<GradType,Ndofs>& func_grad_)=0;

 
  FQPValues<GradType,NQPoints,Ndofs>& reference(const GradientOperator&o){return reference_grad_values_;}
  GradientType& function(const GradientOperator&o){return grad_values_;}


  void init(const QP & qp_points,const GradientOperator& o)
  {

   for(Integer qp=0;qp<NQPoints;qp++)
    {
    qp_points.get_row(qp,qp_point_);
    value(o,qp_point_,grad_);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      reference_grad_values_[n_dof][qp]=grad_[n_dof];           
     }
    }
  };

  template<typename Mapping>
  void operator()(const GradientOperator& o,
                  const Mapping& mapping, 
                  const Vector<Real,Ndofs> &alpha)
  {
    Integer n_tot,n1;
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        n_tot=n_dof * NComponents +  n_comp ;
        n1=n_comp * ShapeFunctionDim1;
        for(Integer qp=0;qp<NQPoints;qp++)
        {
          grad_values_[n_tot][qp].zero();  
          grad_tmp_=alpha[n_dof]*contract(mapping, reference_grad_values_[n_dof][qp]);
          assign(grad_values_[n_tot][qp],grad_tmp_,n1,0);    
        }
      }
     }
  };

  BaseShapeFunctionOperator(){};

  private: 
  Vector<GradType,Ndofs> grad_;
  Point qp_point_;

  FQPValues<GradType,NQPoints,Ndofs> reference_grad_values_;

  GradientType grad_values_;
  GradType grad_tmp_;
  Contraction contract;

};






template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace>
class ShapeFunctionOperator;



template<typename QuadratureRule, Integer NComponents>
class ShapeFunctionOperator<QuadratureRule,Simplex<2,2>, Lagrange1<NComponents> > : 
public BaseShapeFunctionOperator<IdentityOperator, QuadratureRule,3,2,2,1,1,NComponents>,
public BaseShapeFunctionOperator<GradientOperator, QuadratureRule,3,2,2,1,1,NComponents>
{
public:

  static constexpr Integer Ndofs=3;
  static constexpr Integer Dim=2;
  static constexpr Integer ManifoldDim=2;
  static constexpr Integer ShapeFunctionDim1=1;
  static constexpr Integer ShapeFunctionDim2=1;
  static constexpr Integer Ntot = Ndofs * NComponents ;
  static constexpr Integer NQPoints=QuadratureRule::NQPoints;
  using Elem=Simplex<2,2>;
  using BaseFunctionSpace=Lagrange1<NComponents>;
  template<typename Operator>
  using Base=BaseShapeFunctionOperator<Operator,QuadratureRule,Ndofs,Dim,ManifoldDim,ShapeFunctionDim1,ShapeFunctionDim2,NComponents>;


  using Base<IdentityOperator>::reference;
  using Base<GradientOperator>::reference;
  using Base<IdentityOperator>::function;
  using Base<GradientOperator>::function;
  using Base<IdentityOperator>::operator();
  using Base<GradientOperator>::operator();




  using FuncType   = Matrix<Real, ShapeFunctionDim1, ShapeFunctionDim2>;
  using GradType   = Matrix<Real, ShapeFunctionDim1, ShapeFunctionDim2 * Dim>;
  using DivType = Real;
  using Point = Vector<Real,Dim>;
  using QP = Matrix<Real,NQPoints,Dim>;


 virtual void value(const IdentityOperator& o, const Point& point, Vector<FuncType,Ndofs>& func )
       {func[0](0,0)=point[0];
        func[1](0,0)=point[1];
        func[2](0,0)=1. - point[0] - point[1];};
  ;

  virtual void value(const GradientOperator&o, const Point& point, Vector<GradType,Ndofs>& func_grad)
   {
    func_grad[0](0,0)=-1;  func_grad[0](0,1)=-1; 
    func_grad[1](0,0)=+1;  func_grad[1](0,1)= 0; 
    func_grad[2](0,0)= 0;  func_grad[2](0,1)=+1; 
  } 

  template<typename Operator, typename...Operators>
  typename std::enable_if< 0==sizeof...(Operators),void >::type
  init(const QP & qp_points,const Operator&o, const Operators&...operators)
  {Base<Operator>::init(qp_points,o);}


  template<typename Operator, typename...Operators>
  typename std::enable_if< 0<sizeof...(Operators),void >::type
  init(const QP & qp_points,const Operator&o, const Operators&...operators)
  {Base<Operator>::init(qp_points,o);init(qp_points,operators...);}


  template<typename Jacobian, typename Operator, typename...Operators>
  typename std::enable_if< 0==sizeof...(Operators),void>::type 
  operator()(const Jacobian& J, 
             const Vector<Real,Ndofs> &alpha,
             const Operator& o
                  )
  {
  map.init(J,o);
  Base<Operator>::operator()(o,map(o),alpha);
  };

  template<typename Jacobian, typename Operator, typename...Operators>
  typename std::enable_if< 0<sizeof...(Operators),void>::type 
  operator()(const Jacobian& J, 
             const Vector<Real,Ndofs> &alpha,
             const Operator& o,
             const Operators&...operators)
  {
  map.init(J,o);
  Base<Operator>::operator()(o,map(o),alpha);
  operator()(J,alpha,operators...);
  };
 
  ShapeFunctionOperator(){};

  ShapeFunctionOperator(const Integer& add_const):
  range_({add_const,Ntot-1+add_const})
  {}

  template<typename...Operators>
  ShapeFunctionOperator(const QP & qp_points, const Integer& add_const=0,const Operators&...operators):
  range_({add_const,Ntot-1+add_const})
  {init(qp_points,operators...);};

  const std::array<Integer,2>& range(){return range_;};
  inline void range(const Integer& add_const){ range_[0]=add_const;range_[1]=add_const+Ntot-1;};
  inline const Integer& range1(){return range_[1];};

 private:
   MapFromReference4<Elem,BaseFunctionSpace> map;
   std::array<Integer,2> range_;
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













template<typename Elem,typename BaseFunctionSpace,typename Operator,typename QuadratureRule>
class ShapeFunctionOperatorDependent;



template<typename QuadratureRule, Integer NComponents>
class ShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange1<NComponents>,IdentityOperator,QuadratureRule > : 
public BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange1<NComponents>,IdentityOperator,QuadratureRule >  
{
public:
  using Base=BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange1<NComponents>,IdentityOperator,QuadratureRule >;
  using vector_single_type   = typename Base::vector_single_type;
  using Point = typename Base::Point;
  

  virtual void value(const Point& point, vector_single_type& func )
       {func[0](0,0)=point[0];
        func[1](0,0)=point[1];
        func[2](0,0)=1. - point[0] - point[1];};

  ShapeFunctionOperatorDependent()
  {}


};



template<typename QuadratureRule, Integer NComponents>
class ShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange2<NComponents>,IdentityOperator,QuadratureRule > : 
public BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange2<NComponents>,IdentityOperator,QuadratureRule >  
{
public:
  using Base=BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange2<NComponents>,IdentityOperator,QuadratureRule >;
  using vector_single_type   = typename Base::vector_single_type;
  using Point = typename Base::Point;
  

  virtual void value(const Point& point, vector_single_type& func )
       {
        const auto& xi=point[0];
        const auto& eta=point[1];
        const Real zeta = 1. - xi - eta;
        func[0](0,0)=2.*zeta*(zeta-0.5);
        func[1](0,0)=2.*xi*(xi-0.5);
        func[2](0,0)=2.*eta*(eta-0.5);
        func[3](0,0)=4.*zeta*xi;
        func[4](0,0)=4.*xi*eta;
        func[5](0,0)=4.*eta*zeta;};
  ShapeFunctionOperatorDependent()
  {}
};


template<typename QuadratureRule, Integer NComponents>
class ShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange3<NComponents>,IdentityOperator,QuadratureRule > : 
public BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange3<NComponents>,IdentityOperator,QuadratureRule >  
{
public:
  using Base=BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange3<NComponents>,IdentityOperator,QuadratureRule >;
  using vector_single_type   = typename Base::vector_single_type;
  using Point = typename Base::Point;
  

  virtual void value(const Point& point, vector_single_type& func )
       {};
  ShapeFunctionOperatorDependent()
  {}
};




template<typename QuadratureRule, Integer NComponents>
class ShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange1<NComponents>,GradientOperator,QuadratureRule > : 
public BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange1<NComponents>,GradientOperator,QuadratureRule >  
{
public:
  using Base=BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange1<NComponents>,GradientOperator,QuadratureRule >;
  using vector_single_type   = typename Base::vector_single_type;
  using Point = typename Base::Point;
  
  virtual void value(const Point& point, vector_single_type& func_grad)
   {
    func_grad[0](0,0)=-1;  func_grad[0](0,1)=-1; 
    func_grad[1](0,0)=+1;  func_grad[1](0,1)= 0; 
    func_grad[2](0,0)= 0;  func_grad[2](0,1)=+1; 
  } 
  ShapeFunctionOperatorDependent()
  {}
};


template<typename QuadratureRule, Integer NComponents>
class ShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange2<NComponents>,GradientOperator,QuadratureRule > : 
public BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange2<NComponents>,GradientOperator,QuadratureRule >  
{
public:
  using Base=BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange2<NComponents>,GradientOperator,QuadratureRule >;
  using vector_single_type   = typename Base::vector_single_type;
  using Point = typename Base::Point;
  
  virtual void value(const Point& point, vector_single_type& func_grad)
   {
   //TODO, FIXME
   } 
  ShapeFunctionOperatorDependent()
  {}
};

template<typename QuadratureRule, Integer NComponents>
class ShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange3<NComponents>,GradientOperator,QuadratureRule > : 
public BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange3<NComponents>,GradientOperator,QuadratureRule >  
{
public:
  using Base=BaseShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange3<NComponents>,GradientOperator,QuadratureRule >;
  using vector_single_type   = typename Base::vector_single_type;
  using Point = typename Base::Point;
  
  virtual void value(const Point& point, vector_single_type& func_grad)
   {
   //TODO, FIXME
   } 
  ShapeFunctionOperatorDependent()
  {}
};

template<typename QuadratureRule, Integer NComponents>
class ShapeFunctionOperatorDependent<Simplex<2,2>, RT0<NComponents>,IdentityOperator,QuadratureRule > : 
public BaseShapeFunctionOperatorDependent<Simplex<2,2>, RT0<NComponents>,IdentityOperator,QuadratureRule >  
{
public:
  using Base=BaseShapeFunctionOperatorDependent<Simplex<2,2>, RT0<NComponents>,IdentityOperator,QuadratureRule >  ;
  using vector_single_type   = typename Base::vector_single_type;
  using Point = typename Base::Point;
  
  virtual void value(const Point& point, vector_single_type& func)
   {
    const auto& xi=point[0];
    const auto& eta=point[1];
    func[0](0,0)=xi;   func[0](1,0)=eta-1;  
    func[1](0,0)=xi-1; func[1](1,0)=eta;
    func[2](0,0)=xi;   func[2](1,0)=eta;
  } 

  ShapeFunctionOperatorDependent()
  {}


};

template<typename QuadratureRule, Integer NComponents>
class ShapeFunctionOperatorDependent<Simplex<2,2>, RT0<NComponents>,DivergenceOperator,QuadratureRule > : 
public BaseShapeFunctionOperatorDependent<Simplex<2,2>, RT0<NComponents>,DivergenceOperator,QuadratureRule >  
{
public:
  using Base=BaseShapeFunctionOperatorDependent<Simplex<2,2>, RT0<NComponents>,DivergenceOperator,QuadratureRule>  ;
  using vector_single_type   = typename Base::vector_single_type;
  using Point = typename Base::Point;
  
  virtual void value(const Point& point, vector_single_type& func)
   {
    func[0](0,0)=2;  
    func[1](0,0)=2;
    func[2](0,0)=2;
  } 

  ShapeFunctionOperatorDependent()
  {}
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
              func_tmp_=  mapping * reference_values[n_dof][qp];
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
              func_tmp_=alpha[n_dof] * mapping * reference_values[n_dof][qp];             
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


















































































// template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace>
// class Shapeprova
// : public ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>, 
//   public FQPExpression<Shapeprova<QuadratureRule,  Elem, BaseFunctionSpace>,
//   typename ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::TotFuncType,
//   ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::NQPoints,
//   ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Ntot,
//   ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Dim>
// { 
// public:
//   using SFO=ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>;
//   using Map=MapFromReference4<Elem,BaseFunctionSpace>;
//   using S=typename ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>::TotFuncType;
//   static constexpr Integer NQPoints=ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::NQPoints;
//   static constexpr Integer Ntot=ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Ntot;
//   static constexpr Integer  Dim= ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Dim;
//   using Point = Vector<Real,Dim> ;

//   FQPValues<S,NQPoints,Ntot>& apply(const Matrix<Real,NQPoints,Dim> & qp_points)
//      {
//       value_=sfo_->function(Operator::id())();
//       return value_;
//      }

//   Shapeprova(const SFO& sfo, const Map& map):
//     sfo_(std::make_shared<SFO>(sfo)),
//     map_(std::make_shared<Map>(map))
//     {}; 

//    std::shared_ptr<SFO> shape()const{return sfo_;};
//    std::shared_ptr<Map> map()const{return map_;};
//    FQPValues<S,NQPoints,Ntot> value()const{return value_;};
// private:
//   std::shared_ptr<SFO> sfo_;
//   std::shared_ptr<Map> map_;
//   FQPValues<S,NQPoints,Ntot> value_;
// };





// template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace>
// class GradientExpression<QuadratureRule,Elem,BaseFunctionSpace>
// : public ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>, 
//   public FQPExpression<GradientExpression<QuadratureRule,  Elem, BaseFunctionSpace>,
//   typename ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::TotGradType,
//   ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::NQPoints,
//   ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Ntot,
//   ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Dim>
// { 
// public:
//   using SFO=ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>;
//   using Map=MapFromReference4<Elem,BaseFunctionSpace>;
//   using S=typename ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>::TotGradType;
//   static constexpr Integer NQPoints=ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::NQPoints;
//   static constexpr Integer Ntot=ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Ntot;
//   static constexpr Integer  Dim= ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Dim;
//   using Point = Vector<Real,Dim> ;

//   FQPValues<S,NQPoints,Ntot>& apply(const Matrix<Real,NQPoints,Dim> & qp_points)
//      {
//       value_= sp_.shape()->function(Operator::grad())();
//       return value_;
//      }

//   GradientExpression(const Shapeprova<QuadratureRule,Elem,BaseFunctionSpace>& sp):
//     sp_(sp)
//     {}; 
// private:
//   const Shapeprova<QuadratureRule,Elem,BaseFunctionSpace>& sp_;
//   FQPValues<S,NQPoints,Ntot> value_;
// };

// template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace>
// GradientExpression<QuadratureRule,Elem,BaseFunctionSpace> GradProva(const Shapeprova<QuadratureRule,Elem,BaseFunctionSpace>& sp)
// {return GradientExpression<QuadratureRule,Elem,BaseFunctionSpace>(sp);};

// class FQPExpression1: public FQPExpression<FQPExpression1,Matrix<Real,1,2>,6,2,2>
// {
//  public:   
//      static constexpr Integer NQPoints=6;
//      static constexpr Integer NComponents=2;
//      static constexpr Integer Dim=2;

//      using Point = Vector<Real,Dim> ;
//      using T=Matrix<Real,1,2>;
//      using S=Vector<T,NComponents>;

//      S& val(const Point& point)
//      {
//       vec_mat_[0](0,0)=point[0];vec_mat_[0](0,1)=point[0];
//       vec_mat_[1](0,0)=point[1];vec_mat_[1](0,1)=point[1];
//       return vec_mat_;
//      }

//      FQPValues<T,NQPoints,NComponents>& apply(const Matrix<Real,NQPoints,Dim> & qp_points)
//      {
//        for(Integer qp=0;qp<NQPoints;qp++)
//           { 
//             qp_points.get_row(qp,row_);
//             const auto& row_comp=val(row_);
//             for(Integer n_comp=0; n_comp<NComponents; n_comp++)
//                 {
//                   fqp_[n_comp][qp]=row_comp[n_comp];
//                 }
//           }

//       return fqp_;}

//   protected:
//     S vec_mat_;
//     FQPValues<T,NQPoints,NComponents> fqp_;
// };


// class QPExpression1: public QPExpression<QPExpression1,Matrix<Real,3,1>,6,2>
// {
//  public: 
//      static constexpr Integer Dim=2;
//      using Point = Vector<Real,Dim> ;
//       Matrix<Real,3,1>& apply(const Point& point)
//      {
//       mat_(0,0)=point[0];
//       mat_(1,0)=point[1];
//       mat_(2,0)=point[0]+point[1];
//       return mat_;
//      }
//   protected:
//     Matrix<Real,3,1> mat_;
// };


// class QPExpression2: public QPExpression<QPExpression2,Matrix<Real,1,3>,6,2>
// {
//  public: 
//      static constexpr Integer Dim=2;
//      using Point = Vector<Real,Dim> ;
//       Matrix<Real,1,3>& apply(const Point& point)
//      {
//       mat_(0,0)=point[0];
//       mat_(0,1)=point[1];
//       mat_(0,2)=point[0]+point[1];
//       return mat_;
//      }
//   protected:
//     Matrix<Real,1,3> mat_;
// };



// class QPExpression3: public QPExpression<QPExpression3,Matrix<Real,2,1>,6,2>
// {
//  public: 
//      static constexpr Integer Dim=2;
//      using T=Matrix<Real,2,1>;
//      using Point = Vector<Real,Dim> ;
//      T& apply(const Point& point)
//      {
//       mat_(0,0)=point[0];
//       mat_(1,0)=point[1];
//       return mat_;
//      }
//   protected:
//   T mat_;
// };



















template<typename PhiTrial,typename PhiTest,Integer NQPoints> 
Vector<Real, NQPoints> BilinearFormQPvalues(const Vector<PhiTrial,NQPoints>& phi_trial,
                                             const Vector<PhiTest,NQPoints>& phi_test,
                                             const Real& alpha=1.0 )
{
Vector<Real, NQPoints> result;
Contraction contr;
 for(Integer qp=0;qp<NQPoints;qp++)
 {
  result[qp]=alpha * contr(phi_trial[qp],phi_test[qp]);
 }
 return result;
}


template<typename Operation1,typename Operation2, typename MatTrial,typename PhiTrial,typename MatTest,typename PhiTest,Integer NQPoints> 
Vector<Real, NQPoints> BilinearFormQPvalues(const Vector<PhiTrial,NQPoints>& phi_trial,
                                             const Vector<PhiTest,NQPoints>& phi_test,
                                             const Vector<MatTrial,NQPoints>& A,
                                             const Real& alpha=1.0 )
{
Vector<Real, NQPoints> result;
Operation1 prod1;
Contraction contr;
 for(Integer qp=0;qp<NQPoints;qp++)
 {
  auto e1=prod1(A[qp],phi_trial[qp]);
  auto e3=contr(e1,phi_test[qp]);
  result[qp]=e3;
 }
 return result;
}

template<typename Operation1,typename Operation2, typename MatTrial,typename PhiTrial,typename MatTest,typename PhiTest,Integer NQPoints> 
Vector<Real, NQPoints> BilinearFormQPvalues(const Vector<PhiTrial,NQPoints>& phi_trial,
                                             const Vector<PhiTest,NQPoints>& phi_test,
                                             const Vector<MatTrial,NQPoints>& A, 
                                             const Vector<MatTest,NQPoints>& B,
                                             const Real& alpha=1.0 )
{
Vector<Real, NQPoints> result;
Operation1 prod1;
Operation2 prod2;
Contraction contr;
 for(Integer qp=0;qp<NQPoints;qp++)
 {
  auto e1=prod1(A[qp],phi_trial[qp]);
  auto e2=prod2(B[qp],phi_test[qp]);
  auto e3=contr(e1,e2);
  result[qp]=alpha * e3;
 }
 return result;
}

template< Integer NQPoints, typename Jacobian, typename Trial, typename Test, typename MapTrial, typename MapTest, 
          typename ReferenceTrial,typename ReferenceTest, typename AlphaTrial, typename AlphaTest,
          Integer Ndofs_trial, Integer Ndofs_test, Integer NComponents_trial, Integer NComponents_test,
          typename MatrixA, typename MatrixB, typename Matrix, typename QPpoints, typename Volume >
void localmatrix(const Jacobian& J, const Trial& trial, const Test& test, const MapTrial& map_trial,const MapTest& map_test,
                 const ReferenceTrial& reference_trial, const ReferenceTest& reference_test,
                 const AlphaTrial& alpha_trial, const AlphaTest& alpha_test,
                 const MatrixA& A,const MatrixB& B, const QPpoints& qp_weights,const Volume& volume, const Matrix& mat)
{
  //volume=unsigned_volume(simplex,points);
  const auto& mapping_trial=map_trial(J);
  const auto& mapping_test=map_test(J);
  const auto& element_trial=trial.phiN(reference_trial,mapping_trial,alpha_trial);
  const auto& element_test=test.phiN(reference_test,mapping_test,alpha_test);

  for(Integer n_dof_trial=0;n_dof_trial<Ndofs_trial;n_dof_trial++)
    for(Integer n_comp_trial=0;n_comp_trial<NComponents_trial;n_comp_trial++)
     for(Integer n_dof_test=0;n_dof_test<Ndofs_test;n_dof_test++)
      for(Integer n_comp_test=0;n_comp_test<NComponents_test;n_comp_test++)
        {const auto& i=n_comp_trial * NComponents_trial+ n_dof_trial;
         const auto& j=n_comp_test * NComponents_test + n_dof_test;
         mat(i,j)=1;//dot(vec,qp_weights)*volume;
        }
}













template<int N, int... Rest>
struct Array_impl {
    static constexpr auto& value = Array_impl<N - 1, N, Rest...>::value;
};

template<int... Rest>
struct Array_impl<0, Rest...> {
    static constexpr int value[] = { 0, Rest... };
};

template<int... Rest>
constexpr int Array_impl<0, Rest...>::value[];

template<int N>
struct Array {
    static_assert(N >= 0, "N must be at least 0");

    static constexpr auto& value = Array_impl<N>::value;

    Array() = delete;
    Array(const Array&) = delete;
    Array(Array&&) = delete;
};







template<Integer I,Integer J,Integer M, Integer N>
constexpr Integer MatrixIndexElem()
{ static_assert(I<M, " index I >= M in compile time matrix"); 
  static_assert(J<N, " index J >= N in compile time matrix"); 
  return (J + I * N);};

template<Integer S,Integer N,Integer Rows>
class IndexTo;

template<Integer N, Integer Rows>
class IndexTo<0,N,Rows>
{ public:
  static const Integer value=N/Rows;
};

template<Integer N,Integer Rows>
class IndexTo<1,N,Rows>
{ public:
  static const Integer value=Modulo<N,Rows>::value;
};




template<Integer S,typename OperatorType,typename ShapeType,typename QPType, Integer Rows,Integer M,Integer N=0>//Integer I, Integer J>//std::tuple_size<OperatorType>::value,Integer N=0>
struct ReferenceShapeFunctionType
{  


      using operatortype=typename std::tuple_element<N,OperatorType>::type;
      static constexpr Integer shape_index=IndexTo<S,N,Rows>::value;
      operatortype operation;
      using shapetype=typename std::tuple_element<shape_index,ShapeType>::type;
      shapetype shape;
      QPType qp_points;

      using singletype=decltype(shape(operation,qp_points));

       using rest = typename ReferenceShapeFunctionType<S,OperatorType,ShapeType,QPType,Rows,M,N+1>::type;
       using type = decltype(std::tuple_cat(std::declval< std::tuple<singletype> >(), std::declval< rest >() ) );

      template<Integer J>
      using tupletype=typename std::tuple_element<J,type>::type;

};

template<Integer S,typename OperatorType,typename ShapeType,typename QPType, Integer Rows,Integer M>
struct ReferenceShapeFunctionType<S,OperatorType,ShapeType,QPType,Rows,M,M-1>
{  

      static constexpr Integer N=M-1;
      using operatortype=typename std::tuple_element<N,OperatorType>::type;
      static constexpr Integer shape_index=IndexTo<S,N,Rows>::value;
      operatortype operation;
      using shapetype=typename std::tuple_element<shape_index,ShapeType>::type;
      shapetype shape;
      QPType qp_points;
      using singletype = decltype(shape(operation,qp_points));
      using type = typename std::tuple<singletype>;
      // using rest = typename ReferenceShapeFunctionType<OperatorType,ShapeType,QPType,M,N+1>::type;
      // using type = decltype(std::tuple_cat(std::declval< std::tuple<singletype> >(), std::declval< rest >() ) );
      template<Integer J>
      using tupletype=type;

};







template<Integer S,typename ReferenceShapeFunctionType,typename Operator,typename Tupleshapetype_trial,typename QPpoints, Integer M,Integer Rows,Integer N=0>//Integer I, Integer J>//std::tuple_size<OperatorType>::value,Integer N=0>
typename std::enable_if<N==(M-1), void >::type
ReferenceShapeFunction( ReferenceShapeFunctionType& tuple, const Operator& operators_trial,const Tupleshapetype_trial& tupleshape, const QPpoints& qp_points )
{  

constexpr Integer index=IndexTo<S,N,Rows>::value;
using operatortype=typename std::tuple_element<N,Operator>::type;
using shapefunctiontype=typename std::tuple_element<index,Tupleshapetype_trial>::type;

shapefunctiontype shape;
operatortype operation;
std::cout<<"ReferenceShapeFunction=="<<std::endl;
std::cout<<shape(operation,qp_points)<<std::endl;
std::get<N>(tuple)=shape(operation,qp_points);
};

template<Integer S,typename ReferenceShapeFunctionType,typename Operator,typename Tupleshapetype_trial,typename QPpoints,Integer M,Integer Rows,Integer N=0>
typename std::enable_if<N<(M-1), void >::type 
ReferenceShapeFunction( ReferenceShapeFunctionType& tuple, const Operator& operators_trial,const Tupleshapetype_trial& tupleshape, const QPpoints& qp_points )
{  
constexpr Integer index=IndexTo<S,N,Rows>::value;
using operatortype=typename std::tuple_element<N,Operator>::type;
using shapefunctiontype=typename std::tuple_element<index,Tupleshapetype_trial>::type;

shapefunctiontype shape;
operatortype operation;
std::get<N>(tuple)=shape(operation,qp_points);
std::cout<<"ReferenceShapeFunction=="<<std::endl;
std::cout<<shape(operation,qp_points)<<std::endl;
ReferenceShapeFunction<S,ReferenceShapeFunctionType,Operator,Tupleshapetype_trial,QPpoints,M,Rows,N+1>(tuple,operators_trial,tupleshape,qp_points);
};









template<typename TrialShapeFunction, typename TestShapeFunction,typename Mat,typename QPWeight, Integer Ndofs_trial, Integer NComponents_trial,Integer Ndofs_test, Integer NComponents_test>
void local_matrix(const TrialShapeFunction& trial,const TestShapeFunction& test,Mat& mat,const QPWeight& qp_weights,Real volume )
{
    for(Integer n_dof_trial=0;n_dof_trial<Ndofs_trial;n_dof_trial++)
      for(Integer n_comp_trial=0;n_comp_trial<NComponents_trial;n_comp_trial++)
       for(Integer n_dof_test=0;n_dof_test<Ndofs_test;n_dof_test++)
        for(Integer n_comp_test=0;n_comp_test<NComponents_test;n_comp_test++)
               {const auto& i=n_comp_trial * NComponents_trial+ n_dof_trial;
                const auto& j=n_comp_test * NComponents_test + n_dof_test;
                const auto& vec=BilinearFormQPvalues(trial[n_dof_trial][n_comp_trial],test[n_dof_test][n_comp_test]);
                std::cout<<"(n_dof_trial,n_comp_trial,n_dof_test,n_comp_test)=("<<n_dof_trial<<","<<n_comp_trial<<","<<n_dof_test<<","<<n_comp_test<<")"<<std::endl;
                vec.describe(std::cout);
                mat(i,j)=dot(vec,qp_weights)*volume;};

}





template<Integer I, Integer J,Integer Rows, Integer Cols>
class space_loop
{ public:

  template<typename TupleShapeTrial,typename TupleShapeTest,
             typename TupleReferenceTrial,typename TupleReferenceTest,
             typename TupleMapTrial, typename TupleMapTest,
             typename TupleOperatorTrial, typename TupleOperatorTest,
             typename Jacobian, typename QPpoints>
  static const void loop(const TupleShapeTrial& tuple_shape_trial,const TupleShapeTest& tuple_shape_test,
                  const TupleReferenceTrial& tuple_reference_trial,const TupleReferenceTest& tuple_reference_test,
                  const TupleMapTrial& tuple_map_trial, const TupleMapTest& tuple_map_test,
                  const TupleOperatorTrial& tuple_operator_trial, const TupleOperatorTest& tuple_operator_test,
                  const Jacobian& jacob, const QPpoints& qp_points)

  { 
    constexpr Integer N=MatrixIndexElem<I,J,Rows,Cols>();
    std::cout<<"-----------------------------(I,J)= "<<I<<", "<<J<<std::endl;
    std::cout<<"-----------------------------N= "<<N<<std::endl;
    

    const auto& operator_trial=std::get<N>(tuple_operator_trial);
    const auto& operator_test =std::get<N>(tuple_operator_test);


    auto shape_trial=std::get<J>(tuple_shape_trial);
    auto shape_test =std::get<I>(tuple_shape_test);

    const auto& reference_shape_trial=std::get<N>(tuple_reference_trial);//(operator_trial,qp_points);
    const auto& reference_shape_test =std::get<N>(tuple_reference_test);//(operator_test,qp_points);

    // auto trial=shape_trial(operator_trial,qp_points);
    // auto test=shape_test(operator_test,qp_points);
    // Vector<Real,2> punto{0.0,1.0};
    // auto ss=shape_trial.phi_single_component(operator_trial,punto);

    auto map_trial=std::get<J>(tuple_map_trial);
    auto map_test =std::get<I>(tuple_map_test);

    const auto& mapping_trial=map_trial(operator_trial,jacob);
    const auto& mapping_test =map_test (operator_test,jacob);

    // std::cout<<"shape_trial.phi_single_component= "<<std::endl;
    // ss.describe(std::cout);
    // std::cout<<"qp_points= "<<std::endl;
    // qp_points.describe(std::cout);
    // std::cout<<"mapping_test= "<<std::endl;
    // std::cout<<mapping_test<<std::endl;
    // // mapping_trial.describe(std::cout);
    // std::cout<<"jacob= "<<std::endl;
    // jacob.describe(std::cout);
    // std::cout<<"reference_shape_test= "<<std::endl;
    // reference_shape_test.describe(std::cout);
    auto element_trial=shape_trial(operator_trial,reference_shape_trial,mapping_trial);
    auto element_test =shape_test (operator_test, reference_shape_test, mapping_test);
    // local_matrix(element_trial,element_test,mat,qp_weights,volume );
  };

};









template<typename Elem, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
struct TupleShapeFunctionType
{
      using rest = typename TupleShapeFunctionType<Elem, BaseFunctionSpaces...>::type; 
      using ens = ShapeFunctionOperator2<Elem,BaseFunctionSpace>;
      using type = decltype(std::tuple_cat(std::declval< std::tuple<ens> >(), std::declval< rest >() ) );
      template<Integer N>
      using singletype=typename std::tuple_element<N,type>::type;
};



template<typename Elem,typename BaseFunctionSpace>
struct TupleShapeFunctionType<Elem,BaseFunctionSpace>
{    using type = typename std::tuple<ShapeFunctionOperator2<Elem,BaseFunctionSpace>>;
     template<Integer N>
     using singletype=type;

};


template<typename Elem, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
struct TupleMapType
{
      using rest = typename TupleMapType<Elem, BaseFunctionSpaces...>::type; 
      using ens = MapFromReference2<Elem,BaseFunctionSpace>;
      using type = decltype(std::tuple_cat(std::declval< std::tuple<ens> >(), std::declval< rest >() ) );
      template<Integer N>
      using singletype=typename std::tuple_element<N,type>::type;
};



template<typename Elem,typename BaseFunctionSpace>
struct TupleMapType<Elem,BaseFunctionSpace>
{    using type = typename std::tuple<MapFromReference2<Elem,BaseFunctionSpace>>;
     template<Integer N>
     using singletype=type;

};





template<typename Type,typename...Types>
struct VariadicTupleType
{
      using rest = typename VariadicTupleType<Types...>::type; 
      using type = decltype(std::tuple_cat(std::declval< std::tuple<Type> >(), std::declval< rest >() ) );
      template<Integer N>
      using singletype=typename std::tuple_element<N,type>::type;
};



template<typename Type>
struct VariadicTupleType<Type>
{
     using type = typename std::tuple<Type>;
      template<Integer N>
      using singletype=Type;

};







template<typename Type>
struct VariadicTupleTPrint
{
      template<Integer N>
      static void print()
     {
      using Access= typename Type::template singletype<N>;
      std::cout<<"access=="<<Access::FEFamily<<std::endl;
      };
     
};











template<Integer M, Integer N, Integer I=0, Integer J=0, 
         typename TupleShapeTrial,typename TupleShapeTest,
         typename TupleReferenceTrial,typename TupleReferenceTest,
         typename TupleMapTrial, typename TupleMapTest,
         typename TupleOperatorTrial, typename TupleOperatorTest,
         typename Jacobian, typename QPpoints>
typename std::enable_if< I<=(M-1) && J<(N-1), void>::type
MatrixLoop       (const TupleShapeTrial& tuple_shape_trial,const TupleShapeTest& tuple_shape_test,
                  const TupleReferenceTrial& tuple_reference_trial,const TupleReferenceTest& tuple_reference_test,
                  const TupleMapTrial& tuple_map_trial, const TupleMapTest& tuple_map_test,
                  const TupleOperatorTrial& tuple_operator_trial, const TupleOperatorTest& tuple_operator_test,
                  const Jacobian& jacob, const QPpoints& qp_points);

template<Integer M, Integer N, Integer I=0, Integer J=0,  
         typename TupleShapeTrial,typename TupleShapeTest,
         typename TupleReferenceTrial,typename TupleReferenceTest,
         typename TupleMapTrial, typename TupleMapTest,
         typename TupleOperatorTrial, typename TupleOperatorTest,
         typename Jacobian, typename QPpoints>
typename std::enable_if< I==(M-1) && J==(N-1), void>::type
MatrixLoop       (const TupleShapeTrial& tuple_shape_trial,const TupleShapeTest& tuple_shape_test,
                  const TupleReferenceTrial& tuple_reference_trial,const TupleReferenceTest& tuple_reference_test,
                  const TupleMapTrial& tuple_map_trial, const TupleMapTest& tuple_map_test,
                  const TupleOperatorTrial& tuple_operator_trial, const TupleOperatorTest& tuple_operator_test,
                  const Jacobian& jacob, const QPpoints& qp_points)
{std::cout<<"(M,N)=("<<M<<", "<<N<<")----(I,J)=("<<I<<", "<<J<<")"<<std::endl; 
 std::cout<<"Modulo="<<Modulo<MatrixIndexElem<I,J,M,N>(),N>::value <<std::endl;
 std::cout<<"Division="<<MatrixIndexElem<I,J,M,N>()/N <<std::endl;
 std::cout<<MatrixIndexElem<I,J,M,N>()<<std::endl;
 space_loop<I,J,M,N>::loop(tuple_shape_trial,tuple_shape_test,
                  tuple_reference_trial,tuple_reference_test,
                  tuple_map_trial,tuple_map_test,
                  tuple_operator_trial,tuple_operator_test,
                  jacob,qp_points);

 
};

template<Integer M, Integer N, Integer I=0, Integer J=0,
         typename TupleShapeTrial,typename TupleShapeTest,
         typename TupleReferenceTrial,typename TupleReferenceTest,
         typename TupleMapTrial, typename TupleMapTest,
         typename TupleOperatorTrial, typename TupleOperatorTest,
         typename Jacobian, typename QPpoints>
typename std::enable_if< I<(M-1) && J==(N-1), void>::type
MatrixLoop       (const TupleShapeTrial& tuple_shape_trial,const TupleShapeTest& tuple_shape_test,
                  const TupleReferenceTrial& tuple_reference_trial,const TupleReferenceTest& tuple_reference_test,
                  const TupleMapTrial& tuple_map_trial, const TupleMapTest& tuple_map_test,
                  const TupleOperatorTrial& tuple_operator_trial, const TupleOperatorTest& tuple_operator_test,
                  const Jacobian& jacob, const QPpoints& qp_points)
{std::cout<<"(M,N)=("<<M<<", "<<N<<")----(I,J)=("<<I<<", "<<J<<")"<<std::endl;
 std::cout<<"Modulo="<<Modulo<MatrixIndexElem<I,J,M,N>(),N>::value <<std::endl;
 std::cout<<"Division="<<MatrixIndexElem<I,J,M,N>()/N <<std::endl;
std::cout<<MatrixIndexElem<I,J,M,N>()<<std::endl;
 space_loop<I,J,M,N>::loop(tuple_shape_trial,tuple_shape_test,
                  tuple_reference_trial,tuple_reference_test,
                  tuple_map_trial,tuple_map_test,
                  tuple_operator_trial,tuple_operator_test,
                  jacob,qp_points);


 MatrixLoop<M,N,I+1,0>(tuple_shape_trial,tuple_shape_test,
                            tuple_reference_trial,tuple_reference_test,
                            tuple_map_trial,tuple_map_test,
                            tuple_operator_trial,tuple_operator_test,
                            jacob,qp_points);
};


template<Integer M, Integer N, Integer I, Integer J, 
         typename TupleShapeTrial,typename TupleShapeTest,
         typename TupleReferenceTrial,typename TupleReferenceTest,
         typename TupleMapTrial, typename TupleMapTest,
         typename TupleOperatorTrial, typename TupleOperatorTest,
         typename Jacobian, typename QPpoints>
typename std::enable_if< I<=(M-1) && J<(N-1), void>::type
MatrixLoop //(Operator operators_trial,TupleShapeFunction tupleshape, QPPoints qp_points, const Output& output, Parameters...parameters)//const T &tuple,A& array)
                 (const TupleShapeTrial& tuple_shape_trial,const TupleShapeTest& tuple_shape_test,
                  const TupleReferenceTrial& tuple_reference_trial,const TupleReferenceTest& tuple_reference_test,
                  const TupleMapTrial& tuple_map_trial, const TupleMapTest& tuple_map_test,
                  const TupleOperatorTrial& tuple_operator_trial, const TupleOperatorTest& tuple_operator_test,
                  const Jacobian& jacob, const QPpoints& qp_points)
{     std::cout<<"(M,N)=("<<M<<", "<<N<<")----(I,J)=("<<I<<", "<<J<<")"<<std::endl;
      std::cout<<"Modulo="<<Modulo<MatrixIndexElem<I,J,M,N>(),N>::value <<std::endl;
      std::cout<<"Division="<<MatrixIndexElem<I,J,M,N>()/N <<std::endl;
      std::cout<<MatrixIndexElem<I,J,M,N>()<<std::endl;


 space_loop<I,J,M,N>::loop(tuple_shape_trial,tuple_shape_test,
                  tuple_reference_trial,tuple_reference_test,
                  tuple_map_trial,tuple_map_test,
                  tuple_operator_trial,tuple_operator_test,
                  jacob,qp_points);


      MatrixLoop<M,N,I,J+1>(tuple_shape_trial,tuple_shape_test,
                            tuple_reference_trial,tuple_reference_test,
                            tuple_map_trial,tuple_map_test,
                            tuple_operator_trial,tuple_operator_test,
                            jacob,qp_points);
};



















class Mother{
public:
  virtual~Mother(){};
  virtual std::shared_ptr<Mother> get(const Integer i)=0;
  virtual void print()=0;
};

class Son1: public Mother
{
public:
  Son1(const Integer& i): value(i){}
  virtual void print(){std::cout<<value<<std::endl;}
  virtual std::shared_ptr<Mother> get(const Integer i){return std::make_shared<Son1>(*this);}
private:
  Integer  value;
};

class Son2: public Mother
{
 public:
  Son2():children() {};

  virtual void print(){for(Integer ii=0;ii< children.size();ii++)
                        children[ii]->print();}

  virtual std::shared_ptr<Mother> get(const Integer i){ return children[i];};

     void add(const std::shared_ptr<Mother>& kid)
    {
        children.push_back(kid);
    } 

private:
  std::vector<std::shared_ptr<Mother>> children;
};



template<typename...Args>
class Tuple
{
  public:

  template<Integer N>
  using type=typename std::tuple_element<N, std::tuple<Args...> >::type ;


  Tuple(const Args&...args): tuple_(std::make_tuple(args...)){}

  template<Integer N>
  const type<N>& 
  get() const {static_assert(N>=0,"Mars Tuple has only non negative components");
        static_assert(N<sizeof...(Args),"Mars Tuple exceeds maximum dimension");
        return std::get<N>(tuple_);}

  template<Integer N>
  type<N>& 
  get(){static_assert(N>=0,"Mars Tuple has only non negative components");
        static_assert(N<sizeof...(Args),"Mars Tuple exceeds maximum dimension");
        return std::get<N>(tuple_);}




template<Integer N,Integer...M>
class GetType
{ public:
  using type = typename type<N>:: template GetType<M...>::type;};


template<Integer N>
class GetType<N>
{ public: 
   using type =type<N>;};

  protected:
     std::tuple<Args...> tuple_;
};





template<typename...Args>
Tuple<Args...> TupleOf(const Args&...args){return Tuple<Args...>(args...);};
















template <Integer QPOrder, typename TrialFunction,typename TestFunction, typename OperatorTrial, typename OperatorTest>//,typename MeshT >
void BilinearFormIntegrator(const TrialFunction& trialspace,const TestFunction& testspace)
{









using Elem=typename TrialFunction::Elem;
constexpr Integer Dim=Elem::Dim;
constexpr Integer ManifoldDim=Elem::ManifoldDim;  
constexpr Integer Npoints=Elem::Npoints; 

constexpr Integer Nsubspaces_trial=TrialFunction::Nsubspaces;
constexpr Integer Nsubspaces_test =TestFunction::Nsubspaces;
constexpr Integer NtimesN=Nsubspaces_trial*Nsubspaces_test;

// quadrature rule
constexpr Integer NQPoints=GaussPoints<Elem,QPOrder>::NQPoints; 
GaussPoints<Elem,QPOrder> gauss;
const auto& qp_points=GaussPoints<Elem,QPOrder>::qp_points;
const auto& qp_weights=GaussPoints<Elem,QPOrder>::qp_weights;



 using operator_tuple_trial= VariadicTupleType<GradientOperator,IdentityOperator,IdentityOperator,IdentityOperator>;
 using operator_tuple_test = VariadicTupleType<GradientOperator,IdentityOperator,IdentityOperator,IdentityOperator>;

 using operator_tuple_type_trial=typename operator_tuple_trial::type;
 using operator_tuple_type_test =typename operator_tuple_test::type;
 
 operator_tuple_type_trial operators_trial;
 operator_tuple_type_test  operators_test;


  // using vartuple=VariadicTupleType<Lagrange1<2>,Lagrange2<1>,RT0<1>>;
  // std::cout<<"FEFamily="<<vartuple::singletype<0>::FEFamily<<std::endl;
  // std::cout<<"FEFamily="<<vartuple::singletype<1>::FEFamily<<std::endl;
  // std::cout<<"FEFamily="<<vartuple::singletype<2>::FEFamily<<std::endl;
// std::cout<<"FEFamily="<<VariadicTupleType<RT0<1>,Lagrange2<1>>::type3<1>::FEFamily<<std::endl;;
 // VariadicTupleTPrint<vartuple>::print<0>();
 using Tuplemaptype_trial=typename TupleMapType<Simplex<Dim,ManifoldDim>,Lagrange1<2>,Lagrange2<1>>::type;
 using Tuplemaptype_test =typename TupleMapType<Simplex<Dim,ManifoldDim>,Lagrange1<2>,Lagrange2<1>>::type;

 Tuplemaptype_trial tuplemap_trial;
 Tuplemaptype_test  tuplemap_test;

 auto map_trial0=std::get<0>(tuplemap_trial);
 auto map_trial1=std::get<1>(tuplemap_trial);
 // auto map_trial2=std::get<2>(tuplemap_trial);

 using Tupleshapetype_trial=typename TupleShapeFunctionType<Simplex<Dim,ManifoldDim>,Lagrange1<2>,Lagrange2<1>>::type;
 using Tupleshapetype_test =typename TupleShapeFunctionType<Simplex<Dim,ManifoldDim>,Lagrange1<2>,Lagrange2<1>>::type;

 Tupleshapetype_trial tupleshape_trial;
 Tupleshapetype_test  tupleshape_test;


 auto trial0=std::get<0>(tupleshape_trial);
 auto trial1=std::get<1>(tupleshape_trial);

GradientOperator grad;
IdentityOperator identity;
auto reference_trial=trial0(grad,qp_points);


Vector<Real, 3 > alpha_trial{1,1,1};



using QPType=typename GaussPoints<Elem,QPOrder>::qp_points_type;
using ReferenceShapeTrial = ReferenceShapeFunctionType<1,operator_tuple_type_trial,Tupleshapetype_trial,QPType,Nsubspaces_trial,NtimesN>;
using ReferenceShapeTest  = ReferenceShapeFunctionType<0,operator_tuple_type_test, Tupleshapetype_test, QPType,Nsubspaces_test, NtimesN>;

using reference_shape_type_trial=typename ReferenceShapeTrial::type;
using reference_shape_type_test =typename ReferenceShapeTest ::type;



// using reference_shape_trial0= typename std::tuple_element<0,reference_shape_type_trial>::type;// ReferenceShapeTrial ::template tupletype<0>;
// using reference_shape_trial1= typename ReferenceShapeTrial ::template tupletype<1>;
// using reference_shape_trial2= typename ReferenceShapeTrial ::template tupletype<2>;
// using reference_shape_trial3= typename ReferenceShapeTrial ::template tupletype<3>;


// reference_shape_trial0 proviamoci0=trial0(grad,qp_points);
// reference_shape_trial1 proviamoci1=trial1(identity,qp_points);
// reference_shape_trial2 proviamoci2=trial0(identity,qp_points);
// reference_shape_trial3 proviamoci3=trial1(identity,qp_points);

// std::cout<<"proviamoci0"<<std::endl;
// proviamoci0.describe(std::cout);
// std::cout<<"proviamoci1"<<std::endl;
// proviamoci1.describe(std::cout);
// std::cout<<"proviamoci2"<<std::endl;
// proviamoci2.describe(std::cout);
// std::cout<<"proviamoci3"<<std::endl;
// proviamoci3.describe(std::cout);


reference_shape_type_trial reference_shape_trial;
reference_shape_type_test  reference_shape_test;

ReferenceShapeFunction<1,reference_shape_type_trial,operator_tuple_type_trial,Tupleshapetype_trial,QPType, NtimesN, Nsubspaces_trial>(reference_shape_trial,operators_trial,tupleshape_trial,qp_points );
ReferenceShapeFunction<0,reference_shape_type_test, operator_tuple_type_test, Tupleshapetype_test, QPType, NtimesN, Nsubspaces_test> (reference_shape_test, operators_test, tupleshape_test, qp_points );

std::cout<<"ReferenceShapetrial0"<<std::endl;
std::get<0>(reference_shape_trial).describe(std::cout);
std::cout<<"ReferenceShapetrial1"<<std::endl;
std::get<1>(reference_shape_trial).describe(std::cout);
std::cout<<"ReferenceShapetrial2"<<std::endl;
std::get<2>(reference_shape_trial).describe(std::cout);
std::cout<<"ReferenceShapetrial3"<<std::endl;
std::get<3>(reference_shape_trial).describe(std::cout);

std::cout<<"ReferenceShapetest0"<<std::endl;
std::get<0>(reference_shape_trial).describe(std::cout);
std::cout<<"ReferenceShapetest1"<<std::endl;
std::get<1>(reference_shape_trial).describe(std::cout);
std::cout<<"ReferenceShapetest2"<<std::endl;
std::get<2>(reference_shape_trial).describe(std::cout);
std::cout<<"ReferenceShapetest3"<<std::endl;
std::get<3>(reference_shape_trial).describe(std::cout);


// Vector<Real, Ndofs_test > alpha_test=1;
// trial.phiN(grad,reference_trial,mapping_trial,alpha_trial);
// VariadicTupleType<Lagrange2<1>> tupletypes;

//  using tuplee= typename decltype(std::get<0>(tupletypes::type));
// std::cout<<"FEFamily="<<tupletypes::FEFamily<<std::endl;
// std::cout<<"FEFamily="<<std::get<1>(tupletypes)::FEFamily<<std::endl;
const auto& mesh=trialspace.mesh();
// using trial_function=Elem2FunctionSpace<TrialFunction>;
// using test_function=Elem2FunctionSpace<TestFunction>;

// 


// HO UN MAP_TRIAL, MAP_TEST PER OGNI COMPONENTE DELLA MATRICE. 
// MAP_TEST: LO SPAZIO HA LO SPAZIO COSTANTE LUNGO LE RIGHE
// MAP_TRIAL: LO SPAZIO HA LO SPAZIO COSTANTE LUNGO LE COLONNE
// GLI OPERATORI SONO DIVERSI PER TUTTE LE COMPONENTI



 // using TrialFunction= ElemFunctionSpace<Simplex<2,2>, Lagrange1<1>>;
 // using TestFunction= ElemFunctionSpace<Simplex<2,2>, Lagrange2<2>>;

// // trial function infos
// constexpr Integer NsubspacesTrial = TrialFunction::Nsubspaces;
// constexpr Integer NsubspacesTest = TestFunction::Nsubspaces;


// constexpr Integer FEFamily_trial=TrialFunction::FEFamily;
// constexpr Integer NComponents_trial=TrialFunction::NComponents;
// constexpr Integer ShapeFunctionDim_trial=ShapeFunctionOperator<Elem, trial_function>::ShapeFunctionDim;
// constexpr Integer Ndofs_trial=ShapeFunctionOperator<Elem, trial_function>::Ndofs;
// ShapeFunctionOperator<OperatorTrial, Elem, trial_function> trial;
// MapFromReference<OperatorTrial,Elem,FEFamily_trial,NComponents_trial> map_trial;

// // test function infos
// constexpr Integer FEFamily_test=TestFunction::FEFamily;
// constexpr Integer NComponents_test=TestFunction::NComponents;
// constexpr Integer ShapeFunctionDim_test=ShapeFunctionOperator<Elem, test_function>::ShapeFunctionDim;
// constexpr Integer Ndofs_test=ShapeFunctionOperator<Elem, test_function>::Ndofs;
// ShapeFunctionOperator<OperatorTest,Elem, test_function> test;
// MapFromReference<OperatorTest,Elem,FEFamily_test,NComponents_test> map_test;


auto sn=SignedNormal<Elem>(*mesh);
auto normal_sign= sn.sign(); 
// const auto& n_elements=mesh->n_elements();





// std::vector<Vector<Real,Dim>> points(Npoints);
// Matrix<Real, Dim, ManifoldDim> J;
// Real volume;
// Elem elem;


// for(Integer ii=0;ii<Npoints;ii++)
//    elem.nodes[ii]=ii;





// auto reference_trial=trial.phi(qp_points);
// auto reference_test=test.phi(qp_points);


// Vector<Matrix<Real, NComponents_trial,ShapeFunctionDim_trial*NComponents_trial>,NQPoints> A;
// Vector<Matrix<Real, 2,NComponents_test>,NQPoints> B;

// for(Integer qp=0;qp<NQPoints;qp++)
// {
//   A[qp]=2.0;
//   B[qp]=10.0;
// }
// Matrix<Real, Ndofs_test * NComponents_test, Ndofs_trial * NComponents_trial > mat;

//   auto elemnodes_global=mesh.elem(0).nodes;
//   for(Integer mm=0;mm<points.size();mm++)
//      points[mm]=mesh.point(elemnodes_global[mm]);

// std::cout<<"  points  "<<std::endl;
// points[0][0]=0;
// points[0][1]=0;
// points[1][0]=2;
// points[1][1]=3;
// points[2][0]=1;
// points[2][1]=4;

// for(Integer nn=0;nn<points.size();nn++)
// {
//   points[nn].describe(std::cout);
// }



//   jacobian(elem,points,J);



// std::cout<<"BEGIN MATRIX LOOP"<<std::endl;

// MatrixLoop<2,2>(tupleshape_trial,tupleshape_test,
//                 reference_shape_trial,reference_shape_test,
//                 tuplemap_trial,tuplemap_test,
//                 operators_trial,operators_test,
//                 J,qp_points);

// std::cout<<"END MATRIX LOOP"<<std::endl;



  // const auto& mapping_trial=map_trial(grad,J);
  // mapping_trial.describe(std::cout);
//   volume=unsigned_volume(elem,points);
//   alpha_test=normal_sign[0];
//   const auto& mapping_test=map_test(J);
//   const auto& element_test=test.phiN(reference_test,mapping_test,alpha_test);
  
//     for(Integer n_dof_trial=0;n_dof_trial<Ndofs_trial;n_dof_trial++)
//       for(Integer n_comp_trial=0;n_comp_trial<NComponents_trial;n_comp_trial++)
//        for(Integer n_dof_test=0;n_dof_test<Ndofs_test;n_dof_test++)
//         for(Integer n_comp_test=0;n_comp_test<NComponents_test;n_comp_test++)
//                {const auto& i=n_comp_trial * NComponents_trial+ n_dof_trial;
//                 const auto& j=n_comp_test * NComponents_test + n_dof_test;
//                 const auto& vec=BilinearFormQPvalues(element_test[n_dof_trial][n_comp_trial],element_test[n_dof_test][n_comp_test]);
//                 std::cout<<"(n_dof_trial,n_comp_trial,n_dof_test,n_comp_test)=("<<n_dof_trial<<","<<n_comp_trial<<","<<n_dof_test<<","<<n_comp_test<<")"<<std::endl;
//                 vec.describe(std::cout);
//                 mat(i,j)=dot(vec,qp_weights)*volume;};

// std::cout<<"  alpha_test  "<<std::endl;
// alpha_test.describe(std::cout);
// std::cout<<" RT mass matrix "<<std::endl;
// mat.describe(std::cout);
// std::cout<<" RT mass matrix "<<std::endl;
// element_test.describe(std::cout);
// std::cout<<"  volume "<<volume<<std::endl;
// std::cout<<"  qp_weights "<<volume<<std::endl;
// qp_weights.describe(std::cout);

//   for(Integer elem_iter=0 ; elem_iter < 1;elem_iter++)//n_elements ; elem_iter++)
//   {
//     // auto elemnodes_global=mesh->elem(elem_iter).nodes;
//     // for(Integer mm=0;mm<points.size();mm++)
//     //    points[mm]=mesh->point(elemnodes_global[mm]);
//     // jacobian(elem,points,J);
//     std::cout<<IndexTo<0,0,2>::value<<std::endl;
//     std::cout<<IndexTo<0,1,2>::value<<std::endl;
//     std::cout<<IndexTo<0,2,2>::value<<std::endl;
//     std::cout<<IndexTo<0,3,2>::value<<std::endl;

//     std::cout<<IndexTo<1,0,2>::value<<std::endl;
//     std::cout<<IndexTo<1,1,2>::value<<std::endl;
//     std::cout<<IndexTo<1,2,2>::value<<std::endl;
//     std::cout<<IndexTo<1,3,2>::value<<std::endl;
//     volume=unsigned_volume(elem,points);
//     const auto& mapping_trial0=map_trial0(grad,J);
//     const auto& mapping_trial1=map_trial1(identity,J);

// const auto reference_shape_trial0=std::get<0>(reference_shape_trial);

// Vector<Vector<Vector<double, 1>, 6>, 3> reference_shape_trial05=std::get<2>(reference_shape_trial);
// const auto& mapping_trial05=map_trial0(identity,J);

// const auto element_trial0=trial0(grad,reference_shape_trial0,mapping_trial0);
// const auto element_trial05=trial0(identity,reference_shape_trial05,mapping_trial05);

// const auto element_trial1=trial1(identity,std::get<1>(reference_shape_trial),mapping_trial1);
// // const auto element_trial2=trial2(identity,std::get<2>(reference_shape_trial),mapping_trial0);
// // const auto element_trial3=trial3(identity,std::get<3>(reference_shape_trial),mapping_trial0);
// auto naltro=trial0.phi_single_component(std::get<0>(operators_trial),points[0]);
// // std::cout<<"naltro--"<<nm<<std::endl;
// // std::cout<<naltro<<std::endl;
// // std::cout<<"mapping_trial05--"<<nm<<std::endl;
// // std::cout<<mapping_trial05<<std::endl;
// // std::cout<<"reference_shape_trial05--"<<nm<<std::endl;
// // std::cout<<reference_shape_trial05<<std::endl;
// // std::cout<<"jacobian="<<std::endl;
// // J.describe(std::cout);
// // std::cout<<"detJ="<<std::endl;
// // std::cout<<det(J)<<std::endl;
// // std::cout<<"element_trial0="<<std::endl;
// // element_trial0.describe(std::cout);
// // std::cout<<"element_trial1="<<std::endl;
// // element_trial1.describe(std::cout);

// //     const auto& mapping_trial=map_trial(J);
// //     const auto& mapping_test=map_test(J);
// //     const auto& element_trial=trial.phiN(reference_trial,mapping_trial,alpha_trial);
// //     const auto& element_test=test.phiN(reference_test,mapping_test,alpha_test);

// //     for(Integer n_dof_trial=0;n_dof_trial<Ndofs_trial;n_dof_trial++)
// //       for(Integer n_comp_trial=0;n_comp_trial<NComponents_trial;n_comp_trial++)
// //        for(Integer n_dof_test=0;n_dof_test<Ndofs_test;n_dof_test++)
// //         for(Integer n_comp_test=0;n_comp_test<NComponents_test;n_comp_test++)
// //           {const auto& i=n_comp_trial * NComponents_trial+ n_dof_trial;
// //            const auto& j=n_comp_test * NComponents_test + n_dof_test;
// //            // mat(i,j)=dot(vec,qp_weights)*volume;
// //           }

//   }

// ShapeFunctionOperator3<NQPoints, Elem, Lagrange1<1>> s;
// s(identity,qp_points);
// s(grad,qp_points);
// auto ref_identity=s.reference(identity);

// auto ref_grad=s.reference(grad);

// std::cout<<ref_identity()<<std::endl;
// // std::cout<<ref_grad()<<std::endl;

// ShapeFunctionOperator4<GaussPoints<Elem,QPOrder>, Elem, Lagrange1<4>> s4;
// s4(qp_points);
// // s4(qp_points,identity);
// // s4(qp_points,grad);
// // s4(qp_points,identity,identity,grad);
// std::cout<<s4.reference(identity)()<<std::endl;
// const auto& mapping_identity_lagrange1=map_trial0(identity,J);
// const auto& mapping_grad_lagrange1=map_trial0(grad,J);
// // s4(identity,mapping_identity_lagrange1);
// // s4(grad,mapping_grad_lagrange1);
// s4(identity,mapping_identity_lagrange1,grad,mapping_grad_lagrange1);
// std::cout<<s4.function(identity)()<<std::endl;
// std::cout<<s4.function(grad)()<<std::endl;


// Matrix<Real,2,2> mat2{0,1,2,3};


// s4(grad,qp_points);

// using TTT=Matrix<Real,2,1>;
// using SSS=Matrix<Real,2,2>;
// using Tvec=Vector<Real,3>;
// TTT ttt{1,1};
// TTT ttt2{2,2};
// SSS sss{3,3,3,3};
// Tvec tvec1=4;
// Tvec tvec2=5;
// Vector<TTT,2> initvec0=ttt;
// Vector<TTT,6> initvec=ttt;
// Vector<TTT,6> initvec2=ttt2;

// Vector<Tvec,2> init1=tvec1;
// Vector<Tvec,2> init2=tvec2;

// Vector<Vector<TTT,2>,3> initvecvec0=initvec0;
// Vector<Vector<TTT,6>,3> initvecvec=initvec;
// Vector<Vector<TTT,6>,3> initvecvec2=initvec2;

// Vector<Vector<Tvec,2>,3> init_vec=init1;
// Vector<Vector<Tvec,2>,3> init_vec2=init2;


// FQPValues<Tvec,2,3> fqpvalue1(init_vec);
// FQPValues<Tvec,2,3> fqpvalue2(init_vec2);


// FQPValues<TTT,2,3> fqpval0(initvecvec0);
// FQPValues<TTT,6,3> fqpval(initvecvec);
// FQPValues<TTT,6,3> fqpval2(initvecvec2);
// QPValues<SSS,6> qpval(sss);
// QPValues<SSS,2> qpval0(sss);

// std::cout<<"fqpval()"<<std::endl;
// std::cout<<fqpval()<<std::endl;
// std::cout<<"qpval()"<<std::endl;
// std::cout<<qpval()<<std::endl;

// auto result=qpval*fqpval+fqpval;
// std::cout<<"result()"<<std::endl;

// std::cout<<result()<<std::endl;

// // FQPValues<TTT,3,6> result2=qpval*fqpval+fqpval;
// // std::cout<<fqpval()<<std::endl;
// // std::cout<<qpval()<<std::endl;
// // std::cout<<result()<<std::endl;
// // std::cout<<result2()<<std::endl;




//     MapFromReference4<Elem,Lagrange1<1>> map3;
//     using QP=typename GaussPoints<Elem,QPOrder>::qp_points_type;
//     using QuadratureRule=GaussPoints<Elem,QPOrder>;
//     using BaseFunctionSpace=Lagrange1<1>;
//     using ShapeSpace=ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>;
//     ShapeSpace s;
//     map3.init(J,Operator::id());
//     map3.init(J,Operator::grad());
//     const auto& o=Operator::id();
//     const auto& o_grad=Operator::grad();
//     // const auto& mapping=map3(o);
//     s(qp_points);
//     s(o,map3(o));
//     s(o_grad,map3(o_grad));
  






//   FQPExpression1 fqpp1;
//   QPExpression1  qpp1;
//   QPExpression2  qpp2;
//   QPExpression3  qpp3;
//   Shapeprova<GaussPoints<Elem,QPOrder>, Elem, Lagrange1<1>> phi(s,map3);

//   std::cout<<"fqpp1"<<std::endl;
//   std::cout<<fqpp1.eval(qp_points)<<std::endl;
//   std::cout<<"qpp1"<<std::endl;
//   std::cout<<qpp1.eval(qp_points)<<std::endl;
//   auto fqpp2=qpp1*fqpp1;
//   std::cout<<"fqpp2"<<std::endl;
//   std::cout<<fqpp2.eval(qp_points)<<std::endl;
  
//   std::cout<<"qpp2"<<std::endl;
//   std::cout<<qpp2.eval(qp_points)<<std::endl;
//   auto minusqpp2=-qpp2/3.4;
//   auto plusqpp2=+qpp2;
//   auto sumqqp2 = -qpp2+2*qpp2*2;
//   auto minussumqqp2 = -qpp2-qpp2;
//   std::cout<<"minusqpp2"<<std::endl;
//   std::cout<<minusqpp2.eval(qp_points)<<std::endl;
//   std::cout<<"plusqpp2"<<std::endl;
//   std::cout<<plusqpp2.eval(qp_points)<<std::endl;
//   std::cout<<"sumqqp2"<<std::endl;
//   std::cout<<sumqqp2.eval(qp_points)<<std::endl;
//   std::cout<<"minussumqqp2"<<std::endl;
//   std::cout<<minussumqqp2.eval(qp_points)<<std::endl;



//   auto phitwo=phi+phi;
//   auto phitwominus=phi/4.0-phi/2;
//   auto phithree=phitwominus+phitwo;
//   std::cout<<"phi"<<std::endl;
//   std::cout<<phi.eval(qp_points)<<std::endl;

//   std::cout<<"phitwo"<<std::endl;
//   std::cout<<phitwo.eval(qp_points)<<std::endl;

//   std::cout<<"phitwominus"<<std::endl;
//   std::cout<<phitwominus.eval(qp_points)<<std::endl;

//   std::cout<<"phithree"<<std::endl;
//   std::cout<<phithree.eval(qp_points)<<std::endl;


//   auto fqpp3=-qpp2*qpp1*(-fqpp1);
//   std::cout<<"fqpp3"<<std::endl;
//   std::cout<<fqpp3.eval(qp_points)<<std::endl;
//   auto fqpp3equal=qpp2*qpp1*(fqpp1);
//   std::cout<<"fqpp3equal"<<std::endl;
//   std::cout<<fqpp3equal.eval(qp_points)<<std::endl;


//   std::cout<<"phi"<<std::endl;
//   std::cout<<phi.eval(qp_points)<<std::endl;

//   auto phi2= qpp1*phi;
//   std::cout<<"phi2"<<std::endl;
//   std::cout<<phi2.eval(qp_points)<<std::endl;
  
//   auto phi3=+2*GradProva(phi)*3;
//   std::cout<<"+phi3"<<std::endl;
//   std::cout<<phi3.eval(qp_points)<<std::endl;
//   auto phi3minus=-GradProva(phi);
//   std::cout<<"-Grad(phi)"<<std::endl;
//   std::cout<<phi3minus.eval(qp_points)<<std::endl;


// Matrix<Real,2,2> mat_der1{1,2,3,4};
// Matrix<Real,2,2> mat_der2{5,4,3,2};

// auto qpval3=-qpval*3.0+2.5*qpval;
// auto phi4=qpp3*GradProva(phi); 
// std::cout<<"qpp3*Grad(phi)"<<phi4.eval(qp_points)<<std::endl;

// qpval3+=qpval;
// qpval3/=0.5;
// qpval3-=2*qpval;
// std::cout<<"qpval"<<qpval<<std::endl;

// std::cout<<"qpval3"<<qpval3<<std::endl;
// auto fqpval3=fqpval*2-2*fqpval2;
// std::cout<<"fqpval3"<<fqpval3<<std::endl;

// auto fqpvalue3=fqpvalue1*2-fqpvalue2/4;
// std::cout<<"fqpvalue3"<<fqpvalue3<<std::endl;
//  fqpvalue3=-fqpvalue3;
//  std::cout<<"fqpvalue3"<<fqpvalue3<<std::endl;

//  fqpvalue3=+fqpvalue3;
//  std::cout<<"fqpvalue3"<<fqpvalue3<<std::endl;

//  fqpvalue3=2*fqpvalue3;
//  std::cout<<"fqpvalue3"<<fqpvalue3<<std::endl;

//  fqpvalue3=fqpvalue3*2;
//  std::cout<<"fqpvalue3"<<fqpvalue3<<std::endl;

//  fqpvalue3=fqpvalue3/2;
// std::cout<<"fqpvalue3"<<fqpvalue3<<std::endl;


//  auto phi3two= phi3.eval(qp_points)*phi3.eval(qp_points);
//  // auto integral=Integral<QuadratureRule>(phi3,phi3);
//  std::cout<<"phi3two"<<phi3two<<std::endl;
//  // std::cout<<"integral"<<integral.eval(qp_points)<<std::endl;


//  Matrix<Real,2,2> matrix0{0.5,0.5,0.5,0.5};
//  Matrix<Real,2,3> matrix1{1.,1.,1.,1.,1.,1.};
//  Matrix<Real,3,2> matrix2{2.,2.,2.,2.,2.,2.};
//  Matrix<Real,2,2> matrix3{3,3,3,3};
//  Matrix<Real,2,2> matrix4;

// Expression2Matrix<Real,2,2> exprmat0(matrix0);
// Expression2Matrix<Real,2,3> exprmat1(matrix1);
// Expression2Matrix<Real,3,2> exprmat2(matrix2);
// Expression2Matrix<Real,2,2> exprmat3(matrix3);
// auto exprmat=  -2*exprmat0*exprmat1*exprmat2 ;// * exprmat1 * exprmat0+(exprmat0-exprmat1)*exprmat0;
// std::cout<<exprmat0()<<std::endl;

// Evaluation<decltype(exprmat)> eval(exprmat);
// Matrix<Real,2,2> res =eval.apply();
// std::cout<<res<<std::endl;

// Expression2QPValues<SSS,2> qp1(qpval0);
// Expression2FQPValues<TTT,2,3> fqp1(fqpval0);

// auto exprqp=(-qp1+2*qp1)*qp1+2.*qp1/0.5;
// Evaluation<decltype(exprqp)> evalqp(exprqp);
// auto resqp =evalqp.apply();
// std::cout<<resqp<<std::endl;

// auto exprfqp=qp1*(-fqp1+2.*fqp1)+2.*fqp1/0.5;
// Evaluation<decltype(exprfqp)> evalfqp(exprfqp);
// auto resfqp =evalfqp.apply();
// std::cout<<resfqp<<std::endl;


// Expression2MatrixVar<Real,2,2> exprmatvar0;
// Expression2MatrixVar<Real,3,2> exprmatvar1;
// auto exprvar= - exprmatvar1 * exprmatvar0;
// auto exprvar2= -2*exprmatvar1*exprmatvar0/4+exprmatvar1*exprmatvar0*exprmatvar0;

// Evaluation<decltype(exprmatvar1)> evalvar0(exprmatvar1);
// Evaluation<decltype(exprmatvar0)> evalvar1(exprmatvar0);
// Evaluation<decltype(exprvar)> evalvar(exprvar);

// auto resvar0 =evalvar0.apply(qp_points);
// auto resvar1 =evalvar1.apply(qp_points);
// auto resvar =evalvar.apply(qp_points);
// std::cout<<resvar0<<std::endl;
// std::cout<<resvar1<<std::endl;
// std::cout<<resvar<<std::endl;

// auto lagr=ShapeFunctionOperator<GaussPoints<Elem,QPOrder>, Elem, Lagrange1<4>>();
// ShapeFunctionOperator<GaussPoints<Elem,QPOrder>, Elem, Lagrange1<4>> lagrangesf;
// lagrangesf.init(qp_points,Operator::id(),Operator::grad());
// lagrangesf(J,Vector<Real,3>(1),Operator::id(),Operator::grad());//,Vector<Real,3>(1.0));

// std::cout<<"lagrange_ identity"<<lagrangesf.function(identity)()<<std::endl;
// std::cout<<"lagrange_grad"<<lagrangesf.function(grad)()<<std::endl;


// OperatorType<Matrix<Real,2,2>>::type ess;
// std::cout<<"QuadratureOrder="<<QuadratureOrder<decltype(exprvar2)>::value<<std::endl;


// TupleCatType< std::tuple<int,float>, std::tuple<int>, std::tuple<char> > test{1,1.0f,2, 'a'};

// std::cout<<std::get<3>(test)<<std::endl;
// auto tuplecat=std::tuple_cat(std::tuple<int, char>(1,'a'),std::tuple<int, char,int>(1,'a',2),std::tuple<int>(1));
//     RemoveTupleDuplicates<std::tuple<int, char, char, int, char>> a(1,'a');
//     static_assert(std::is_same<RemoveTupleDuplicates<std::tuple<int, char>>, std::tuple<int, char>>::value, "ok");
//     static_assert(std::is_same<RemoveTupleDuplicates<std::tuple<int, char, int>>, std::tuple<int, char>>::value, "ok");
//     static_assert(std::is_same<RemoveTupleDuplicates<std::tuple<int, int, char>>, std::tuple<int, char>>::value, "ok");
//     static_assert(std::is_same<RemoveTupleDuplicates<std::tuple<int, char, char, int, char>>, std::tuple<int, char>>::value, "ok");

    
//    TupleRemoveType<int,std::tuple<char, float,int,int> >::type removed('a',2.3);
//    std::cout<<std::get<0>(removed)<<std::endl;
//    std::cout<<std::get<1>(removed)<<std::endl;

//    TupleRemovesingleType<int,std::tuple<int,char, int,float,int,int> > removed2('a',2,2.3,1,3);
//    std::cout<<std::tuple_size<decltype(removed2)>::value<<std::endl;
//    std::cout<<std::get<0>(removed2)<<std::endl;
//    std::cout<<std::get<1>(removed2)<<std::endl;
//    std::cout<<std::get<2>(removed2)<<std::endl;
//    std::cout<<std::get<3>(removed2)<<std::endl;
//    std::cout<<std::get<4>(removed2)<<std::endl;

// ShapeFunctionExpression<1,2,Lagrange1<1>> sfe1; 
// GradientShapeFunctionExpression<0,2,RT0<1>> sfe2;
// GradientShapeFunctionExpression<1,2,Lagrange1<1>> sfe3;
// ShapeFunctionExpression<0,2,RT0<2>> sfe4; 

// auto sfen=sfe1+sfe2+sfe3+sfe4;

// OperatorTupleType<decltype(sfen)>::type ec;




// constexpr Integer NComponents=4;
// ShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange1<NComponents>,IdentityOperator,QuadratureRule >  sfod;
// sfod.init_reference();
// sfod.init(J);
// sfod.function();
// std::cout<<"sfod.function()=="<<sfod.function()<<std::endl;

// ShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange1<NComponents>,GradientOperator,QuadratureRule >  sfod_grad;
// sfod_grad.init_reference();
// sfod_grad.init(J);
// sfod_grad.function();
// std::cout<<"sfod_grad.function()=="<<sfod_grad.function()<<std::endl;


// ShapeFunctionOperatorDependent<Simplex<2,2>, RT0<NComponents>,IdentityOperator,QuadratureRule>  sfod_rt;
// sfod_rt.init_reference();
// sfod_rt.init(J);
// sfod_rt.function();
// std::cout<<"sfod_rt.function()=="<<sfod_rt.function()<<std::endl;

// ShapeFunctionOperatorDependent<Simplex<2,2>, RT0<NComponents>,DivergenceOperator,QuadratureRule >  sfod_rt_div;
// sfod_rt_div.init_reference();
// sfod_rt_div.init(J);
// sfod_rt_div.function();
// std::cout<<"sfod_rt_div.function()=="<<sfod_rt_div.function()<<std::endl;
 

};



















}



#endif





