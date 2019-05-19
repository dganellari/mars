// Ho una classe ElemIntegral< Elem, QuadratureRule, Order> che ha un metodo integral che prende una funzione f,
// e calcola l'integrale di f sull'elemento (di riferimento?), quindi:
// ElemIntegral< Elem, QuadratureRule, Order>::integral(f)
// oppure f
// ElemIntegral< Elem, QuadratureRule, Order>::integral(f,actual_elem)
// per calcolare l'integrale sull'elemento attuale.

// Quindi, integral prende i pesi e li moltiplica per f valutata nei punti di quadratura.

// Devo dare f come input generico. Quindi devo poter creare f(p)=(A+1)(f_1(p))*f_2(p)*f_3(p).
// Oppure f=(A+1)(f_1) * f_2 * f_3



// class Function;

// class FunctionToIntegrate
// {

//   private:
//          Function f_;

//     FunctionToIntegrate operator * (const FunctionToIntegrate& rhs) const
//     {
//         FunctionToIntegrate p;
//         p.f_=this->f_ * rhs.f_;
//         p.a = this;
//         p.b = &rhs;
//         return p;
//     }

// }

// function(Point p)
// {return {this->f_(p)* rhs.f_(p);}


// template<typename T>
// {


  
// }

//   class gradient() prende una funzione di base generica  ene ritorna il gradiente

//   grad(phi)(p)
//   cosi posso fare:
//   gradient(shape_function)            
//   oppure
//   A= (I+div)(shape_function)
//   dove definisco I(shape_function)=shape_function e div(shape_function)=div_shape_function







#ifndef MARS_SHAPE_FUNCTION_HPP
#define MARS_SHAPE_FUNCTION_HPP

#include "mars_simplex.hpp"
#include "mars_base_elementfunctionspace.hpp"
#include "mars_vector.hpp"

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
      std::vector< Vector<Integer, ManifoldDim + 1 > > outward_;
public:

      std::vector< Vector< Vector< Real, Dim> , ManifoldDim + 1 > > operator () () const { return normal_; };
      std::vector< Vector<Integer, ManifoldDim + 1 > > sign() const {return outward_;}; 

      template< typename MeshT>
      SignedNormal(MeshT mesh)
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

            // save the normal and if its outward or inward 
            auto n = normal(simplex_elem, elempoints_sort);
            normal_[ee][mm]=n;
            outward_[ee][mm]=sgn(dot(diff,n));
              }
            }
       };

      template< typename MeshT>
      void print(MeshT mesh)
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





























template<typename TensorType>
class DotProduct;


template<>
class DotProduct<Real>
{

};

template<Integer Dim>
class DotProduct<Vector<Real,Dim>>
{

 public:
 static const Real compute(const Vector<Real,Dim>& v1, const Vector<Real,Dim>& v2, const double& alpha=1)
 {
  Real result=0.0;
  for(Integer nn=0;nn<Dim;nn++)
  {
   result+=v1[nn]*v2[nn];
  }
  return result*alpha;
 };
};

template<Integer Rows, Integer Cols>
class DotProduct<Matrix<Real,Rows,Cols>>
{

 public:


 static const Real compute(const Matrix<Real,Rows,Cols>& m1, const Matrix<Real,Rows,Cols>& m2, const Real& alpha=1)
 {
  Real result=0.0;
  for(Integer nn=0;nn<Rows;nn++)
    {
      for(Integer mm=0;mm<Cols;mm++)
      {
       result+=m1(nn,mm)*m2(nn,mm);
      }
    }
  return result*alpha;
 }

  static const Real compute(const Matrix<Real,Rows,Cols>& m1, const Matrix<Real,Rows,Cols>& m2, const Matrix<Real,Rows,Cols>& m3)
 {
    Real result=0.0;
    for(Integer nn=0;nn<Rows;nn++)
      {
        for(Integer mm=0;mm<Cols;mm++)
        {
         result+=m1(nn,mm)*m2(nn,mm)*m3(nn,mm);
        }
      } 
    return result;
   }
};




template<typename TensorType, Integer NComponents=1>
class QPVector;


template<Integer Rows>
class QPVector<Vector<Real,Rows>,1>
{
 public:
 // compute (U,V)
 static const Vector<Real,Rows> compute(const Vector<Real,Rows>& U, 
                                        const Vector<Real,Rows>& V)
 {
  Vector<Real,Rows> result;
  for(Integer mm=0;mm<Rows;mm++)
       result[mm]=U[mm]*V[mm];
  return result;
 }
};


template<Integer Rows, Integer Cols, Integer NComponents>
class QPVector<Matrix<Real,Rows,Cols>,NComponents>
{

 public:

 // compute (U,V)
 static const Vector<Real,Rows> compute(const Matrix<Real,Rows,Cols>& U, 
                                        const Matrix<Real,Rows,Cols>& V)
 {
  Vector<Real,Rows> result;

  for(Integer mm=0;mm<Rows;mm++)
    {
      result[mm]=0;
      for(Integer nn=0;nn<Cols;nn++)
      {
       result[mm]+=U(mm,nn)*V(mm,nn);
      }
    }
  return result;
 }


 // compute (U,V)
 static const Vector<Real,Rows> compute(const Real& alpha, 
                                        const Matrix<Real,Rows,Cols>& U,         
                                        const Matrix<Real,Rows,Cols>& V)
 {
  Vector<Real,Rows> result;

  for(Integer mm=0;mm<Rows;mm++)
    {
      result[mm]=0;
      for(Integer nn=0;nn<Cols;nn++)
      {
       result[mm]+=U(mm,nn)*V(mm,nn);
      }
      result[mm]*=alpha;
    }
  return result;
 }

 // compute (A U, V),A=const
 static const Vector<Real,Rows> compute(const Matrix<Real,Cols,Cols>& A, const Matrix<Real,Rows,Cols>& U, 
                                                                         const Matrix<Real,Rows,Cols>& V)
 {
  Vector<Real,Rows> result;
  Vector<Real,Cols> rowU;
  Vector<Real,Cols> rowV;
  for(Integer mm=0;mm<Rows;mm++)
    {
      result[mm]=0;
      U.get_row(mm,rowU);
      V.get_row(mm,rowV);
      rowU=A*rowU;
      for(Integer nn=0;nn<Cols;nn++)
      {
       result[mm]+=rowU(nn)*rowV(nn);
      }
    }
  return result;
 }

 // compute (A U, V),A=const
 static const Vector<Real,Rows> compute(                                 const Matrix<Real,Rows,Cols>& U, 
                                        const Matrix<Real,Cols,Cols>& B, const Matrix<Real,Rows,Cols>& V)
 {
  Vector<Real,Rows> result;
  Vector<Real,Cols> rowU;
  Vector<Real,Cols> rowV;
  for(Integer mm=0;mm<Rows;mm++)
    {
      result[mm]=0;
      U.get_row(mm,rowU);
      V.get_row(mm,rowV);
      rowV=B*rowV;
      for(Integer nn=0;nn<Cols;nn++)
      {
       result[mm]+=rowU(nn)*rowV(nn);
      }
    }
  return result;
 }

 // compute (A U,B V),A=const, B =const
 static const Vector<Real,Rows> compute(const Matrix<Real,Cols,Cols>& A, const Matrix<Real,Rows,Cols>& U, 
                                        const Matrix<Real,Cols,Cols>& B, const Matrix<Real,Rows,Cols>& V)
 {
  Vector<Real,Rows> result;
  Vector<Real,Cols> rowU;
  Vector<Real,Cols> rowV;
  for(Integer mm=0;mm<Rows;mm++)
    {
      result[mm]=0;
      U.get_row(mm,rowU);
      V.get_row(mm,rowV);
      rowU=A*rowU;
      rowV=B*rowV;
      for(Integer nn=0;nn<Cols;nn++)
      {
       result[mm]+=rowU(nn)*rowV(nn);
      }
    }
  return result;
 }

 // compute (A U,B V),A=const, B =const

 static const Matrix< Vector<Real,Rows>, NComponents,NComponents> compute(
                                        const Vector< Matrix<Real,NComponents,NComponents>,Rows>& A,
                                        const Vector< Matrix<Real,NComponents,NComponents>,Rows>& B)
 {
  Matrix< Vector<Real,Rows>, NComponents,NComponents> mat;
  for(Integer qp=0;qp<Rows;qp++)
      {
       auto mat_qp=transpose(A[qp])*B[qp];
       for(Integer ii=0; ii< NComponents; ii++)
           for(Integer jj=0; jj< NComponents; jj++)
              mat(ii,jj)[qp]=mat_qp(ii,jj); 
      }
  return mat;
 }

};


template<Integer Rows>
using QPVecV=QPVector<Vector<Real,Rows>,1>;

template<Integer Rows,Integer Cols,Integer NComponents=1>
using QPVecM=QPVector<Matrix<Real,Rows,Cols>,NComponents>;

















template<typename Elem,Integer Order>
class GaussPoints;

template< Integer Dim >
class GaussPoints< Simplex<Dim,2> , 1>
{
private:
Matrix<Real,1,2> qp_points_;
Vector<Real,1> weights_;

public:
  static constexpr Integer NQPpoints=1;
  const Matrix<Real,1,2> qp_points()const {return qp_points_;};
  const Vector<Real,1> weights()const {return weights_;};

  GaussPoints<Simplex<Dim,2>,1>():
  qp_points_({0.33333333333333, 0.33333333333333}),
  weights_({1})
   {}; 

};



template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,2>
{
private:
  Matrix<Real,3,2>  qp_points_;
  Vector<Real,3> weights_;

public:
  static constexpr Integer NQPpoints=3;
  const Matrix<Real,3,2> qp_points()const {return qp_points_;};
  const Vector<Real,3> weights()const {return weights_;};
  GaussPoints<Simplex<Dim,2>,2>():
  qp_points_({0.16666666666667, 0.16666666666667,
              0.16666666666667, 0.66666666666667,
              0.16666666666667, 0.16666666666667}),
  weights_({0.33333333333333,0.33333333333333,0.33333333333333})
  {}
};




template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,3>
{
private:
  Matrix<Real,4,2>  qp_points_;
  Vector<Real,4> weights_;

public:
  static constexpr Integer NQPpoints=4;
  const Matrix<Real,4,2> qp_points()const {return qp_points_;};
  const Vector<Real,4> weights()const {return weights_;};
  GaussPoints<Simplex<Dim,2>,3>():
  qp_points_({0.33333333333333, 0.33333333333333,
              0.20000000000000, 0.20000000000000,
              0.20000000000000, 0.60000000000000,
              0.60000000000000, 0.20000000000000}),
  weights_({ -0.56250000000000,0.52083333333333,0.52083333333333,0.52083333333333})
  {}
};

template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,4>
{
private:
  Matrix<Real,6,2>  qp_points_;
  Vector<Real,6> weights_;

public:
  static constexpr Integer NQPpoints=6;
  const Matrix<Real,6,2> qp_points()const {return qp_points_;};
  const Vector<Real,6> weights()const {return weights_;};
  GaussPoints<Simplex<Dim,2>,4>():
  qp_points_({0.44594849091597, 0.44594849091597,
              0.44594849091597, 0.10810301816807,
              0.10810301816807, 0.44594849091597, 
              0.09157621350977, 0.09157621350977, 
              0.09157621350977, 0.81684757298046,
              0.81684757298046, 0.09157621350977} ),
  weights_({0.22338158967801, 
            0.22338158967801,
            0.22338158967801,
            0.10995174365532,
            0.10995174365532, 
            0.10995174365532})
  {}
};


template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,5>
{
private:
  Matrix<Real,7,2>  qp_points_;
  Vector<Real,7> weights_;

public:
  static constexpr Integer NQPpoints=7;
  const Matrix<Real,7,2> qp_points()const {return qp_points_;};
  const Vector<Real,7> weights()const {return weights_;};
  GaussPoints<Simplex<Dim,2>,5>():
  qp_points_(   {0.33333333333333, 0.33333333333333, 
                 0.47014206410511, 0.47014206410511,
                 0.47014206410511, 0.05971587178977,
                 0.05971587178977, 0.47014206410511, 
                 0.10128650732346, 0.10128650732346, 
                 0.10128650732346, 0.79742698535309, 
                 0.79742698535309, 0.10128650732346}  ),
  weights_({0.22500000000000, 
            0.13239415278851, 
            0.13239415278851, 
            0.13239415278851, 
            0.12593918054483, 
            0.12593918054483, 
            0.12593918054483 })
  {}
};





template<Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents=1>
class DivPhi;

template<Integer Dim, Integer ManifoldDim,Integer NComponents>
class DivPhi<Dim,ManifoldDim,1,NComponents>
{
public:
using type = Matrix<Real,1,1>;
  // this class is for not vector shape function. 
  // this call is ambiguous for 1D elements. avoid them.
};


template<Integer Dim,Integer NComponents>
class DivPhi<Dim,Dim,Dim,NComponents>
{
  public:
  using type = Matrix<Real,NComponents,1>;
  static_assert(Dim>1,"Vector function spaces need to be at least in 2D (use lagrange elements for 1D, they are equivalent)");
  // this class is for not vector shape function. 
};

template<Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents>
using DivPhiType=typename DivPhi<Dim, ManifoldDim, ShapeFunctionDim, NComponents>::type;










template<Integer Ndofs_, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents=1>
class BaseShapeFunction
{public:
  static constexpr Integer Ndofs=Ndofs_;
  virtual ~BaseShapeFunction(){};

  virtual const Matrix<Real,Ndofs,ShapeFunctionDim>
  phi_single_component(const Vector<Real,Dim>& point )=0;

  virtual const Matrix<Real,Ndofs,ShapeFunctionDim>
  phi_single_component(const Vector<Real,Dim>& point,const std::vector<Vector<Real, Dim>>& elem_points)=0;


  virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
  grad_phi_single_component(const Vector<Real,Dim>& point)=0;


  virtual const Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs>
  div_phi(const Vector<Real,Dim>& qp_point )=0;

  template<Integer NQPpoints>
  const Vector<Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>, NQPpoints>,Ndofs>
  div_phiN(const Matrix<Real,NQPpoints,Dim>& qp_points)
  {
        Vector<Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>, NQPpoints>,Ndofs> div;
        Vector<Real,Dim> qp_point;
        for(Integer qp=0;qp<NQPpoints;qp++)
        {
          qp_points.get_row(qp,qp_point);
          auto ecco=div_phi(qp_point);
              for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
                  div[n_dof][qp]=ecco[n_dof];
        }

        return div;};






  template<Integer NQPpoints>
  const Vector<Matrix<Real,NQPpoints,ShapeFunctionDim>,Ndofs>
  reference_shape_function(const Matrix<Real,NQPpoints,Dim>& qp_points)
  {
        Vector<Matrix<Real,NQPpoints,ShapeFunctionDim>,Ndofs> shapefunction;
        Vector<Real,Dim> point;
        for(Integer nn=0;nn<Ndofs;nn++)
          {
            for(Integer qp=0;qp<NQPpoints;qp++)
            {
              qp_points.get_row(qp,point);
              auto& sf=this->phi_single_component(point);
              for(Integer dim=0;dim<ShapeFunctionDim;dim++)
              shapefunction[nn](qp,dim)=sf(nn,dim);  
            }            
          }

        return shapefunction;};

  template<Integer NQPpoints>
  const Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPpoints>,NComponents>,Ndofs>
  phi(const Matrix<Real,NQPpoints,Dim>& qp_points)
  {
   Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPpoints>,NComponents>,Ndofs> f;
   Vector<Real,ShapeFunctionDim> row_tmp;
   Vector<Real,Dim> qp_point;

   for(Integer qp=0;qp<NQPpoints;qp++)
    {
    qp_points.get_row(qp,qp_point);
    auto phi_single=phi_single_component(qp_point);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      phi_single.get_row(n_dof,row_tmp);
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        f[n_dof][n_comp][qp].zero();
        f[n_dof][n_comp][qp].row(n_comp,row_tmp);
      }
     }
    }

   return f;
  };


  template<Integer NQPpoints>
  const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPpoints>,NComponents>,Ndofs>
  grad_phi(const Matrix<Real,NQPpoints,Dim>& qp_points)
  {
   Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPpoints>,NComponents>,Ndofs> f;
   Vector<Real,Dim> qp_point;

   for(Integer qp=0;qp<NQPpoints;qp++)
    {
    qp_points.get_row(qp,qp_point);
    auto grad_phi_single=grad_phi_single_component(qp_point);
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


  template<Integer NQPpoints, typename Mapping>
  const Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPpoints>,NComponents>,Ndofs>
  phiN(const Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPpoints>,NComponents>,Ndofs>& reference_phi,
       const Mapping& mapping, 
       const Vector<Real,Ndofs> &alpha=1.0)
  {
   Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPpoints>,NComponents>,Ndofs> result;
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        for(Integer qp=0;qp<NQPpoints;qp++)
        {
          result[n_dof][n_comp][qp]= alpha[n_dof] * mapping * reference_phi[n_dof][n_comp] [qp];
        }
        
      }
     }
    

   return result;
  };


  template<Integer NQPpoints, typename Mapping>
  const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPpoints>,NComponents>,Ndofs>
  grad_phiN(const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPpoints>,NComponents>,Ndofs>& reference_grad_phi,
       const Mapping& mapping, 
       const Vector<Real,Ndofs> &alpha=1.0)
  {
   Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPpoints>,NComponents>,Ndofs> result;
   Vector<Real,Dim> row;
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        for(Integer qp=0;qp<NQPpoints;qp++)
        {
          reference_grad_phi[n_dof][n_comp][qp].get_row(n_comp,row); 
          std::cout<<" firs row ="<<std::endl;
          reference_grad_phi[n_dof][n_comp][qp].describe(std::cout);     
          row.describe(std::cout);     
          std::cout<<" mapping ="<<std::endl;
          mapping.describe(std::cout);    
          row=alpha[n_dof]*mapping*row;      
          result[n_dof][n_comp][qp].row(n_comp,row);
          std::cout<<" result ="<<std::endl;
          result[n_dof][n_comp][qp].describe(std::cout);        
        }
      }
     }
   return result;
  };



template<Integer NQPpoints, typename Mapping>
const Vector< Vector< Vector< Matrix<Real,NComponents,ShapeFunctionDim>, NQPpoints>,NComponents>, Ndofs>
shape_function(const Vector<Matrix<Real,NQPpoints,ShapeFunctionDim>,Ndofs> &reference_shape_function,
          const Mapping &mapping,
          const Real& alpha = 1)
{
 Vector< Vector< Vector< Matrix<Real,NComponents,ShapeFunctionDim>, NQPpoints>,NComponents>, Ndofs> res;
 Vector<Real, ShapeFunctionDim> row;
 Matrix<Real, NComponents,ShapeFunctionDim> shape;

 for(Integer n_dof=0;n_dof<Ndofs; n_dof++)
  {
    auto& shape_n_dof = reference_shape_function[n_dof];
    for(Integer n_comp=0;n_comp<NComponents; n_comp++)
    {

      for(Integer qp=0;qp<NQPpoints; qp++)
           {
            shape.zero();
            shape_n_dof.get_row(qp,row);
            shape.row(n_comp,row);
            //shape.describe(std::cout);
            auto new_shape=mapping*shape*alpha;
            //std::cout<<"(n_dof, n_comp, qp)=( "<<n_dof<<", "<<n_comp<<", "<<qp<<" )"<<std::endl;
            //new_shape.describe(std::cout);
            res[n_dof][n_comp][qp]= new_shape;
          }
    }
  }
  return res;
};



};



























template<typename Elem,typename BaseFunctionSpace>
class ShapeFunction;








template<Integer NComponents_>
class ShapeFunction<Simplex<2,2>, Lagrange1<NComponents_> > : 
public BaseShapeFunction<3,2,2,1,NComponents_>
{
public:
  static constexpr Integer Dim=2;
  static constexpr Integer ManifoldDim=2;
  static constexpr Integer ShapeFunctionDim=1;
  static constexpr Integer NComponents=NComponents_;
  static constexpr Integer Ndofs=3;

  virtual const Matrix<Real,3,1>
  phi_single_component(const Vector<Real,2>& point)
       {const auto& xi=point[0];
        const auto& eta=point[1];
        const Real zeta = 1. - xi - eta;
        const Matrix<Real,3,1> shapefunction{zeta,xi,eta};
        return shapefunction;};


  virtual const Matrix<Real,3,1>
  phi_single_component(const Vector<Real,2>& point,const std::vector<Vector<Real, 2>>& elem_points)
  {
    Matrix<Real,3,1> mat;
    return mat;
    // do stuff
  };



  virtual const Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs>
  div_phi(const Vector<Real,Dim>& qp_point )
  {
    static_assert(NComponents==ManifoldDim, "divergence for non vector shape functions requires: NComponents==ManifoldDim");
    static_assert(NComponents==Dim, "divergence for non vector shape functions requires: NComponents==Dim");
  Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs> tryme=1; 

    tryme[0](0,0)=1.0;
    tryme[1](0,0)=1.0;
    tryme[2](0,0)=-2.0;

  return tryme;};



  virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
  grad_phi_single_component(const Vector<Real,Dim>& point)
  {
    Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs> grad{{-1,-1},{1,0},{0,1}};
    return grad;
  }

};


template<Integer NComponents_>
class ShapeFunction<Simplex<2,2>, Lagrange2<NComponents_> > : 
public BaseShapeFunction<6,2,2,1,NComponents_>
{
public:
  static constexpr Integer Ndofs=6;
  static constexpr Integer Dim=2;
  static constexpr Integer ManifoldDim=2;
  static constexpr Integer ShapeFunctionDim=1;
  static constexpr Integer NComponents=NComponents_;
  

      virtual const Matrix<Real,6,1>
      phi_single_component(const Vector<Real,2>& point) 
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

  virtual const Matrix<Real,6,1>
  phi_single_component(const Vector<Real,2>& point,const std::vector<Vector<Real, 2>>& elem_points)
  {
     Matrix<Real,6,1> shape_function;
     return shape_function;
    // do stuff
  };

  virtual const Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs>
  div_phi(const Vector<Real,Dim>& qp_point )
  {
    Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs> tryme; 
    tryme[0](0,0)=(4.0+4.0);
    tryme[1](0,0)=(4.0);
    tryme[2](0,0)=(4.0);
    tryme[3](0,0)=(-8.0);
    tryme[4](0,0)=(4.0);
    tryme[5](0,0)=(-8.0);
    return tryme;
 };



  virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
  grad_phi_single_component(const Vector<Real,Dim>& point)
  {
   const auto& xi=point[0];
   const auto& eta=point[1];
   const Real zeta = 1. - xi - eta;
   Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs> grad{
                                     { - 4 * (1 - xi - eta) + 1, - 4 * (1 - xi - eta) + 1},
                                     { 4 * xi - 1              , 0                       },
                                     { 0                       , 4 * eta -1              },
                                     { 4 - 8 * xi - 4 * eta    , - 4 * xi                },
                                     { 4 * eta                 , 4 * xi                  },
                                     { -4 * eta                , 4 - 4 * xi - 8 * eta    } };
  return grad;
  }

};





template<Integer NComponents_>
class ShapeFunction<Simplex<2,2>, RT0<NComponents_> > : 
public BaseShapeFunction<3,2,2,2,NComponents_>
{
public:
  static constexpr Integer Ndofs=3;
  static constexpr Integer Dim=2;
  static constexpr Integer ManifoldDim=2;
  static constexpr Integer ShapeFunctionDim=2;
  static constexpr Integer NComponents=NComponents_;
  

  virtual const  Matrix<Real,3,2>
  phi_single_component(const Vector<Real,2>& point)
       {const auto& xi=point[0];
        const auto& eta=point[1];
        const Matrix<Real,3,2> shapefunction{xi,eta-1,xi-1,eta,xi,eta};
        return shapefunction;};


  virtual const Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs>
  div_phi(const Vector<Real,Dim>& qp_point )
  {
   Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs> tryme;
   for(Integer n_comp=0;n_comp<NComponents;n_comp++)
   {
     tryme[0](n_comp,0)=2;
     tryme[1](n_comp,0)=2;
     tryme[0](n_comp,0)=2;    
   }
   return tryme;
  };



  virtual const Matrix<Real,3,2>
  phi_single_component(const Vector<Real,2>& point,const std::vector<Vector<Real, 2>>& elem_points)
  {


    const Matrix<Real,3,2> shapefunction;
    return shapefunction;
    // do stuff
  };

 virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
 grad_phi_single_component(const Vector<Real,2>& point)
       {
        //static_assert(Dim==ShapeFunctionDim, "Grad(RT): shape function dim and space dim must be the same")
        Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs> grad{ {1,0, 0,1},{1,0, 0,1},{1,0, 0,1}} ;
        return grad;
        assert("gradient is not defined for RT elements");};



};



























template<typename Elem,Integer FEFamily, Integer Order, Integer Continuity, Integer NComponents=1>
class MapFromReference;

template<Integer Dim, Integer ManifoldDim,Integer Order,Integer Continuity,Integer NComponents>
class MapFromReference<Simplex<Dim,ManifoldDim>,LagrangeFE,Order,Continuity,NComponents>
{
 private:
 public:

        Matrix<Real, Dim, ManifoldDim>  map_grad
        (const Simplex<Dim,ManifoldDim>& simplex,
         const std::vector<Vector<Real, Dim>> &points)
         {
          static_assert(Dim==ManifoldDim,"Dim==ManifoldDim for inverting jacobian");
          Matrix<Real, Dim, ManifoldDim> J;
          jacobian(simplex,points,J); 
          auto inv=  inverse(J);
          return inv;
         }

        Real  map
        (const Simplex<Dim,ManifoldDim>& simplex,
         const std::vector<Vector<Real, Dim>> &points)
         {
          return 1.0;
         }


        template<Integer Ndofs, Integer ShapeFunctionDim>
        Matrix<Real,Ndofs,ShapeFunctionDim> operator() 
                        (const Matrix<Real,Ndofs,ShapeFunctionDim> reference_shape_function,
                         const std::vector<Vector<Real, Dim>> &points,
                         const Integer& sign_shape_function=1)
                   {
                    return reference_shape_function;
                   }

};


template<Integer Dim, Integer ManifoldDim,Integer Order,Integer Continuity,Integer NComponents>
class MapFromReference<Simplex<Dim,ManifoldDim>,RaviartThomasFE,Order,Continuity,NComponents>
{
 private:
 public:
        
        Matrix<Real, Dim, ManifoldDim>  map(const Simplex<Dim,ManifoldDim>& simplex,
                                            const std::vector<Vector<Real, Dim>> &points)
         {
          Matrix<Real, Dim, ManifoldDim> mapping;
          jacobian(simplex,points,mapping);
          mapping/=det(mapping);
          return mapping;
         }

        Real  map_div(const Simplex<Dim,ManifoldDim>& simplex,
                      const std::vector<Vector<Real, Dim>> &points)
         {
          Real mapping_div;
          Matrix<Real, Dim, ManifoldDim> J;
          jacobian(simplex,points,J);
          mapping_div=1.0/det(J);
          return mapping_div;
         }

        inline const Matrix<Real, Dim, ManifoldDim> map(const Matrix<Real, Dim, ManifoldDim>& J)
         {Matrix<Real, Dim, ManifoldDim> mapping=J;
          mapping/=det(mapping);
          return mapping;}

        inline const Real  map_div (const Matrix<Real, Dim, ManifoldDim>& J){return 1.0/det(J);}


        template<Integer Ndofs, Integer ShapeFunctionDim>
        Matrix<Real,Ndofs,ShapeFunctionDim> operator() 
                        (const Matrix<Real,Ndofs,ShapeFunctionDim> reference_shape_function,
                         const std::vector<Vector<Real, Dim>> &points,
                         const Integer& sign_shape_function=1)
                   {
                    static_assert(ShapeFunctionDim==ManifoldDim, "MapFromReference:  ShapeFunctionDim must be equal to ManifoldDim");

                    Simplex<Dim,ManifoldDim> simplex;
                    Matrix<Real,Ndofs,ShapeFunctionDim> result;
                    // // WHY DO WE NEED TO DO THIS EVERYTIME?
                    for(Integer ii=0;ii<ManifoldDim+1;ii++)
                        simplex.nodes[ii]=ii;
                    auto mapping=map(simplex,points);

                    // jacobian(simplex,points,mapping);
                    // auto detmapping=det(mapping);
                    // std::cout<<"J=="<<std::endl;
                    // mapping.describe(std::cout);

                    // mapping/=(detmapping*sign_shape_function);
                    // std::cout<<"J/detJ=="<<std::endl;
                    // mapping.describe(std::cout);


                    Vector<Real,ShapeFunctionDim> row;
                    for(Integer ii=0;ii<Ndofs;ii++)
                    {
                      reference_shape_function.get_row(ii,row);
                      result.row(ii,mapping*row);
                    }
                    // std::cout<<"reference_shape_function=="<<std::endl;
                    // reference_shape_function.describe(std::cout);
                    // std::cout<<"result=="<<std::endl;
                    // result.describe(std::cout);


                    return result;

//                      //return ;
                   }

};



template<typename T, Integer Rows1,Integer Cols1,Integer Rows2,Integer Cols2>
Matrix<Real, Rows1*Rows2, Cols1*Cols2 > tensorproduct(const Matrix<T, Rows1, Cols1>& mat1,const Matrix<T, Rows2, Cols2>& mat2 )
{
 Matrix<Real, Rows1*Rows2, Cols1*Cols2 > mat;

  for(Integer i1=0;i1<Rows1;i1++)
  {
    for(Integer j1=0;j1<Cols1;j1++)
    {
      for(Integer i2=0;i2<Rows2;i2++)
      {
       for(Integer j2=0;j2<Cols2;j2++)
        {
          mat(i1*Rows2+i2 ,j1*Cols2+j2)= DotProduct<T>::compute(mat1(i1,j1),mat2(i2,j2));
        }
      }
    }
  }
  return mat;
}













template <Integer QPOrder, typename TrialFunction,typename TestFunction, typename MeshT >
void MassIntegrator(const MeshT& mesh)
{
using trial_function=Elem2FunctionSpace<TrialFunction>;
using test_function=Elem2FunctionSpace<TestFunction>;

using Elem=typename TrialFunction::Elem;
constexpr Integer Dim=Elem::Dim;
constexpr Integer ManifoldDim=Elem::ManifoldDim;  

// trial function infos
constexpr Integer FEFamily_trial=TrialFunction::FEFamily;
constexpr Integer Order_trial=TrialFunction::Order;
constexpr Integer Continuity_trial=TrialFunction::Continuity;
constexpr Integer NComponents_trial=TrialFunction::NComponents;
constexpr Integer ShapeFunctionDim_trial=ShapeFunction<Elem, trial_function>::ShapeFunctionDim;
constexpr Integer Ndofs_trial=ShapeFunction<Elem, trial_function>::Ndofs;
ShapeFunction<Elem, trial_function> trial;

// test function infos
constexpr Integer FEFamily_test=TestFunction::FEFamily;
constexpr Integer Order_test=TestFunction::Order;
constexpr Integer Continuity_test=TestFunction::Continuity;
constexpr Integer NComponents_test=TestFunction::NComponents;
constexpr Integer ShapeFunctionDim_test=ShapeFunction<Elem, test_function>::ShapeFunctionDim;
constexpr Integer Ndofs_test=ShapeFunction<Elem, test_function>::Ndofs;
ShapeFunction<Elem, test_function> test;


MapFromReference<Simplex<Dim,ManifoldDim>,FEFamily_trial,Order_trial,Continuity_trial,NComponents_trial> map_trial;
MapFromReference<Simplex<Dim,ManifoldDim>,FEFamily_test,Order_test,Continuity_test,NComponents_test> map_test;


// quadrature rule
constexpr Integer NQPpoints=GaussPoints<Elem,QPOrder>::NQPpoints; 
GaussPoints<Elem,QPOrder> gauss;
auto qp_points=gauss.qp_points();
//Matrix<Real,Ndofs*NComponents,Ndofs*NComponents> mass;
std::vector<Vector<Real,Dim>> points(ManifoldDim+1);
Simplex<Dim,ManifoldDim> simplex;
for(Integer ii=0;ii<ManifoldDim+1;ii++)
   simplex.nodes[ii]=ii;



//Vector<Matrix<Real,NQPpoints,ShapeFunctionDim_trial>,Ndofs_trial> 


auto elemnodes_global=mesh.elem(0).nodes;
  for(Integer mm=0;mm<points.size();mm++)
     points[mm]=mesh.point(elemnodes_global[mm]);

auto mapping_trial=map_trial.map(simplex,points);
//auto mapping_grad_trial=map_trial.map_grad(simplex,points);

// std::cout<<"mapping_grad_trial=="<<std::endl;
// mapping_grad_trial.describe(std::cout);
//auto mapping_grad_trial=map_trial.map_grad(simplex,points);
 Matrix<Real,3,3> matt=1.0;//{1.0,1.0,1.0,0.0,4.0,1.0,0.0,0.0,5.0};
Vector<Real,3> vec=2.0;
 std::cout<<"MATT=="<<std::endl;
 vec.describe(std::cout);

 matt.describe(std::cout);
// Matrix<Real,3,3>  matt2;
// matt2=inverse(matt);
// std::cout<<"matt2*=4=="<<std::endl;
// matt2*=4;
// matt2.describe(std::cout);
// std::cout<<"matt2/=4=="<<std::endl;
// matt2/=4;
// matt2.describe(std::cout);

// std::cout<<"matt2+=matt=="<<std::endl;
// matt2+=matt;
// matt2.describe(std::cout);

// std::cout<<"matt2=(matt2/4)=="<<std::endl;
// matt2=(matt2/4);
// matt2.describe(std::cout);

// std::cout<<"matt2=(matt2*4)=="<<std::endl;
// matt2=(matt2*4);
// matt2.describe(std::cout);


auto shape_trial=trial.reference_shape_function(qp_points);
Vector<Real,Dim > row;
auto reference_trial=trial.phi(qp_points);
auto grad_reference_trial=trial.grad_phi(qp_points);
auto element_trial=trial.phiN(reference_trial,mapping_trial);
auto div_phi_trial=trial.div_phiN(qp_points);
std::cout<<"div_phi_trial.describe=="<<std::endl;

div_phi_trial.describe(std::cout);
// auto element_grad_trial=trial.grad_phiN(grad_reference_trial,mapping_grad_trial);

// std::cout<<"shape_trial.describe=="<<std::endl;
// shape_trial.describe(std::cout);
// std::cout<<"phi_trial.describe=="<<std::endl;
// reference_trial.describe(std::cout);
// std::cout<<"element_trial.describe=="<<std::endl;
// element_trial.describe(std::cout);
std::cout<<"grad_reference_trial.describe=="<<std::endl;
grad_reference_trial.describe(std::cout);
//Vector<Matrix<Real,NQPpoints,ShapeFunctionDim_test>,Ndofs_test>
auto shape_test=test.reference_shape_function(qp_points);

// vector for each dof of vector for each qp point of the shape_function (with all components) in qp
Vector< Vector< Vector< Matrix<Real,NComponents_trial,ShapeFunctionDim_trial>, NQPpoints>,NComponents_trial>, Ndofs_trial> sf_trial;
//Vector<Matrix<Real,NQPpoints,ShapeFunctionDim_test>,Ndofs_test*NComponents_test> sf_test;





   for(Integer nn_trial=0;nn_trial<Ndofs_trial;nn_trial++)
    for(Integer cc_trial=0;cc_trial<NComponents_trial;cc_trial++)
       {}//sf_trial[nn_trial][cc_trial]



for(Integer ee=0;ee<1;ee++)//mesh.n_elements();ee++)
{
  Elem elem=mesh.elem(ee);
  auto elemnodes_global=elem.nodes;
  for(Integer mm=0;mm<points.size();mm++)
    {
     points[mm]=mesh.point(elemnodes_global[mm]);
     std::cout<<"points="<<std::endl;
     //points[mm].describe(std::cout);
    }

  auto mapping_trial=map_trial.map(simplex,points);
  //std::cout<<"mapping_trial="<<std::endl;
 // mapping_trial.describe(std::cout);

  
  auto elemenet_shape_trial=trial.shape_function(shape_trial,mapping_trial);
  for(Integer ii=0;ii<Ndofs_test;ii++)
      {
       auto phi_test=map_test(shape_test[ii],points);
         for(Integer jj=0;jj<Ndofs_trial;jj++)
            {
              auto phi_trial=map_trial(shape_trial[jj],points);
              //auto& qpvec=QPVecM<NQPpoints,ShapeFunctionDim_trial>::compute( phi_test[ii], phi_trial[jj]);
            }
      }
}


// std::cout<<" shape_trial"<<std::endl;
// for(Integer ii=0;ii<shape_trial.size();ii++)
//      shape_trial[ii].describe(std::cout);
// std::cout<<" shape_test"<<std::endl;
// for(Integer ii=0;ii<shape_test.size();ii++)
//      shape_test[ii].describe(std::cout);

std::cout<<"----------------"<<std::endl;


//constexpr Integer NComponents=ShapeFunction<Elem, trialfun>::NComponents;
// Matrix<Real,NComponents,NComponents> matrix_a{1.0};
// Matrix<Real,NComponents,NComponents> matrix_b{1.0};


//Matrix<Real,NComponents_trial,NComponents_trial> matrix_a{1.0};
//Matrix<Real,NComponents_trial,NComponents_trial> matrix_a{1.0,0.0,0.0,1.0};
//Matrix<Real,NComponents_trial,NComponents_trial> matrix_b{1.0};
//Matrix<Real,NComponents_trial,NComponents_trial> matrix_b{1.0,0.0,0.0,1.0};

Matrix<Real,Ndofs_test,Ndofs_trial> mass;

std::cout<<"Ndofs_test,   Ndofs_trial======="<<Ndofs_test<<", "<<Ndofs_trial<<std::endl;

// Vector<decltype(matrix_a),NQPpoints> mata(matrix_a);
// Vector<decltype(matrix_b),NQPpoints> matb(matrix_b);
//std::cout<<" mata"<<std::endl;

// for(Integer qp=0;qp<NQPpoints;qp++)
//     mata[qp].describe(std::cout);

// Matrix< Vector<Real,NQPpoints>, NComponents_trial,NComponents_trial> qpmat=QPVecM<NQPpoints,ShapeFunctionDim_trial,NComponents_trial>::compute(mata,mata);
Matrix<Vector<Real,NQPpoints> ,Ndofs_test,Ndofs_trial> mat_mass;


// vector, NQPpoints long, whose components are the matrix A evaluated in different qp points
//Vector< Matrix<Real,NComponents,NComponents> , NQPpoints> matA;
 for(Integer nn_test=0;nn_test<Ndofs_test;nn_test++)
  for(Integer cc_test=0;cc_test<NComponents_test;cc_test++)
   for(Integer nn_trial=0;nn_trial<Ndofs_trial;nn_trial++)
    for(Integer cc_trial=0;cc_trial<NComponents_trial;cc_trial++)
        {
        }

for(Integer nn=0;nn<Ndofs_test;nn++)
   for(Integer mm=0;mm<Ndofs_trial;mm++)
   {

    // Quello che voglio fare alla fine e' una combinazione di matrici. Per H1 ad esempio:
    // (A gradu, gradv) + (B u,v) = (weights,  A GradUGradV+BUV)
    // (A divS,B divT)
    // dove BUV e' il vettore 
    // un elemento e' sempre il prodotto scalare tra wieghts e un'altro vettore lungo NQPpoints
    // quindi se ho il prodotto tra due Matrix(n_qppoints,Dim) mat, devo restituire
    // Vector(n_qppoints ) vec, dove vec(qp)=sum( mat1(qp,ii)*mat2(qp,ii))

    
    auto& qpvec=QPVecM<NQPpoints,ShapeFunctionDim_trial>::compute( shape_test[nn], shape_trial[mm]);
    std::cout<<" (nn,mm)==("<<nn<<", "<<mm<<")"<<std::endl;
    mat_mass(nn,mm)=QPVecV<NQPpoints>::compute(qpvec,gauss.weights());
    //mass(nn,mm)=DotProduct<Vector<Real,NQPpoints>>::compute(qpvec,gauss.weights());
   }
// std::cout<<"mass"<<std::endl;
// mass.describe(std::cout);
// std::cout<<"mat_mass"<<std::endl;
// mat_mass.describe(std::cout);
// std::cout<<"qpmat"<<std::endl;
// qpmat.describe(std::cout);
// std::cout<<"tensor product"<<std::endl;
// auto mattensor=tensorproduct(mat_mass,qpmat);
// mattensor.describe(std::cout);

// auto sn=SignedNormal<ManifoldDim>(mesh);

// sn.print();




    // Triangle2 tri2;
    // tri2.nodes[0] = 0;
    // tri2.nodes[1] = 1;
    // tri2.nodes[2] = 2;


    // Tetrahedron4 tet; 
    // tet.nodes[0] = 0;
    // tet.nodes[1] = 1;
    // tet.nodes[2] = 2;
    // tet.nodes[3] = 3;

    // std::vector<Vector2r> points2(
    // {
      
      
    //    { 1., 0. },
    //    { 0., 0. },
    //    { 0., 1. }
    // });

    // std::vector<Vector4r> points4(
    // {
    //   { 2., 0., 0., 0. },
    //   { 0., 2., 0., 0. },
    //   { 0., 0., 2., 0. },
    //   { 0., 0., 0., 2. },
    //   { 0., 0., 0., 0. }
      
      
      
      
    // });
    // auto n2 = normal(tri2, points2);
    // auto n4 = normal(tet, points4);

    // n2.describe(std::cout);
    // n4.describe(std::cout);
    // std::cout << n4 << std::endl;

};


// template< Integer Dim>
// shape_function<Simplex<Dim,2>, RaviartThomasFE, 0 >(const std::array<double,2>& point, std::array<double,3>)
// {
// const auto& xi=point[0];
// const auto& eta=point[1];
// const double zeta = 1. - xi - eta;

// shape_function[0]={xi,eta-1};
// shape_function[0]={xi-1,eta};
// shape_function[0]={xi,eta};
// }





}



#endif