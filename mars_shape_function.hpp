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


class IdentityOperator;
class DivergenceOperator;
class GradientOperator;

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
      std::vector< Vector<Real, ManifoldDim + 1 > > sign() const {return outward_;}; 

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


















class Multiply
{
 public:
 template<typename T,Integer Rows1,Integer Cols1,Integer Rows2, Integer Cols2>
 Matrix<T, Rows1,Cols2> operator()(const Matrix<T, Rows1,Cols1> &A,const Matrix<T, Rows2,Cols2> &B)
 {
       static_assert(Cols1==Rows2,"Multiply: matrix dimension must agree, i.e. Cols1==Rows2 ");
       return A*B;
 }

 template<typename T,Integer Rows,Integer Cols>
 Matrix<T, Rows,Cols> operator()(const Matrix<T, Rows,Cols> &A,const T &B)
 {
       return A*B;
 }

 template<typename T,Integer Rows,Integer Cols>
 Matrix<T, Rows,Cols> operator()(const T &A,const Matrix<T, Rows,Cols> &B)
 {
       return A*B;
 }

 template<typename T,Integer Rows,Integer Cols>
 Vector<T,Rows> operator()(const Matrix<T, Rows,Cols> &A,const Vector<T,Cols> &B)
 {
       return A*B;
 };

 template<typename T,Integer Rows,Integer Cols>
 Vector<T,Rows> operator()(const Matrix<T, Rows,Cols> &A,const Matrix<T,Cols,1> &B)
 {
       Vector<T,Cols> vec;
       B.get_col(0,vec);
       return A*vec;
 };

};

class Contraction
{
public:

 template<typename T,Integer Rows,Integer Cols>
 T operator()(const Matrix<T, Rows,Cols> &A,const Matrix<T, Rows,Cols> &B)
 {
       T result=0;
       for(Integer i = 0; i < Rows; ++i) 
         for(Integer j = 0; j < Cols; ++j) 
             result += A(i, j) * B(i,j);
       return result;
 }

 template<typename T>
 T operator()(const Matrix<T, 1,1> &A,const T &B)
 {
       return A(0,0)*B;
 }

 template<typename T>
 T operator()(const T &A,const Matrix<T, 1,1> &B)
 {
       return A*B(0,0);
 }

 template<typename T,Integer Dim>
 T operator()(const Vector<T, Dim> &A,const Vector<T, Dim> &B)
 {
       T result=0;
       for(Integer i = 0; i < Dim; ++i) 
             result += A[i] * B[i];
       return result;
 }

 template<typename T>
 T operator()(const T &A,const T &B)
 {
       return A*B;
 }

 template<typename T,Integer Rows>
 T operator()(const Matrix<T, Rows,1> &A,const Vector<T, Rows> &B)
 {
       T result=0;
       for(Integer i = 0; i < Rows; ++i) 
             result += A(i, 0) * B[i];
       return result;
 }

 template<typename T,Integer Cols>
 T operator()(const Matrix<T, 1,Cols> &A,const Vector<T, Cols> &B)
 {
       T result=0;
       for(Integer i = 0; i < Cols; ++i) 
             result += A(0,i) * B[i];
       return result;
 }

 template<typename T,Integer Rows>
 T operator()(const Vector<T, Rows> &B,const Matrix<T, Rows,1> &A)
 {
       T result=0;
       for(Integer i = 0; i < Rows; ++i) 
             result += A(i, 0) * B[i];
       return result;
 }

 template<typename T,Integer Cols>
 T operator()(const Vector<T, Cols> &B,const Matrix<T, 1,Cols> &A)
 {
       T result=0;
       for(Integer i = 0; i < Cols; ++i) 
             result += A(0,i) * B[i];
       return result;
 }

};





template<typename Operation>
class Product{
 public:
 template<typename Input1,typename Input2, typename Output>
 Output operator()(Input1 i1,Input2 i2){return Operation(i1,i2);};
};


// template<>
// class GeneralProduct<Contraction>
// {
// public:
//  template<typename Input1,typename Input2, typename Output>
//  Output operator()(Input1 i1,Input2 i2)
//  {return Contraction(i1,i2);};

// };





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
Vector<Real,1> qp_weights_;

public:
  static constexpr Integer NQPpoints=1;
  const Matrix<Real,1,2> qp_points()const {return qp_points_;};
  const Vector<Real,1> qp_weights()const {return qp_weights_;};

  GaussPoints<Simplex<Dim,2>,1>():
  qp_points_({0.33333333333333, 0.33333333333333}),
  qp_weights_({1})
   {}; 

};



template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,2>
{
private:
  Matrix<Real,3,2>  qp_points_;
  Vector<Real,3> qp_weights_;

public:
  static constexpr Integer NQPpoints=3;
  const Matrix<Real,3,2> qp_points()const {return qp_points_;};
  const Vector<Real,3> qp_weights()const {return qp_weights_;};
  GaussPoints<Simplex<Dim,2>,2>():
  qp_points_({0.16666666666667, 0.16666666666667,
              0.16666666666667, 0.66666666666667,
              0.16666666666667, 0.16666666666667}),
  qp_weights_({0.33333333333333,0.33333333333333,0.33333333333333})
  {}
};




template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,3>
{
private:
  Matrix<Real,4,2>  qp_points_;
  Vector<Real,4> qp_weights_;

public:
  static constexpr Integer NQPpoints=4;
  const Matrix<Real,4,2> qp_points()const {return qp_points_;};
  const Vector<Real,4> qp_weights()const {return qp_weights_;};
  GaussPoints<Simplex<Dim,2>,3>():
  qp_points_({0.33333333333333, 0.33333333333333,
              0.20000000000000, 0.20000000000000,
              0.20000000000000, 0.60000000000000,
              0.60000000000000, 0.20000000000000}),
  qp_weights_({ -0.56250000000000,0.52083333333333,0.52083333333333,0.52083333333333})
  {}
};

template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,4>
{
private:
  Matrix<Real,6,2>  qp_points_;
  Vector<Real,6> qp_weights_;

public:
  static constexpr Integer NQPpoints=6;
  const Matrix<Real,6,2> qp_points()const {return qp_points_;};
  const Vector<Real,6> qp_weights()const {return qp_weights_;};
  GaussPoints<Simplex<Dim,2>,4>():
  qp_points_({0.44594849091597, 0.44594849091597,
              0.44594849091597, 0.10810301816807,
              0.10810301816807, 0.44594849091597, 
              0.09157621350977, 0.09157621350977, 
              0.09157621350977, 0.81684757298046,
              0.81684757298046, 0.09157621350977} ),
  qp_weights_({0.22338158967801, 
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
  Vector<Real,7> qp_weights_;

public:
  static constexpr Integer NQPpoints=7;
  const Matrix<Real,7,2> qp_points()const {return qp_points_;};
  const Vector<Real,7> qp_weights()const {return qp_weights_;};
  GaussPoints<Simplex<Dim,2>,5>():
  qp_points_(   {0.33333333333333, 0.33333333333333, 
                 0.47014206410511, 0.47014206410511,
                 0.47014206410511, 0.05971587178977,
                 0.05971587178977, 0.47014206410511, 
                 0.10128650732346, 0.10128650732346, 
                 0.10128650732346, 0.79742698535309, 
                 0.79742698535309, 0.10128650732346}  ),
  qp_weights_({0.22500000000000, 
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
























template<Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents=1>
class DivPhiOperator;

template<Integer Dim, Integer ManifoldDim,Integer NComponents>
class DivPhiOperator<Dim,ManifoldDim,1,NComponents>
{
public:
using type = Real;//Matrix<Real,1,1>;
  // this class is for not vector shape function. 
  // this call is ambiguous for 1D elements. avoid them.
};


template<Integer Dim,Integer NComponents>
class DivPhiOperator<Dim,Dim,Dim,NComponents>
{
  public:
  using type = Matrix<Real,1,1>;
  static_assert(Dim>1,"Vector function spaces need to be at least in 2D (use lagrange elements for 1D, they are equivalent)");
  // this class is for not vector shape function. 
};

template<Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents>
using DivPhiOperatorType=typename DivPhiOperator<Dim, ManifoldDim, ShapeFunctionDim, NComponents>::type;




template<typename Operator, Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents=1>
class BaseShapeFunctionOperator;

template<Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents>
class BaseShapeFunctionOperator<IdentityOperator,Ndofs,Dim,ManifoldDim,ShapeFunctionDim,NComponents>
{public:

  virtual ~BaseShapeFunctionOperator(){};

  virtual const Matrix<Real,Ndofs,ShapeFunctionDim>
  phi_single_component(const Vector<Real,Dim>& point )=0;

  template<Integer NQPpoints>
  const Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPpoints>,Ndofs>
  phi(const Matrix<Real,NQPpoints,Dim>& qp_points)
  {
   Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPpoints>,Ndofs> f;
   //Vector<Real,ShapeFunctionDim> row_tmp;
   Vector<Real,Dim> qp_point;

   for(Integer qp=0;qp<NQPpoints;qp++)
    {
    qp_points.get_row(qp,qp_point);
    auto phi_single=phi_single_component(qp_point);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer ii=0;ii<ShapeFunctionDim;ii++)
        f[n_dof][qp][ii]=phi_single(n_dof,ii);     
     }
    }
   return f;
  };

  template<Integer NQPpoints, typename Mapping>
  const Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPpoints>,NComponents>,Ndofs>
  phiN(const Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPpoints>,Ndofs>& reference_phi,
       const Mapping& mapping, 
       const Vector<Real,Ndofs> &alpha=1.0)
  {
   Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPpoints>,NComponents>,Ndofs> result;
   Vector<Real,ShapeFunctionDim> row;
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        for(Integer qp=0;qp<NQPpoints;qp++)
        { 
          result[n_dof][n_comp][qp].zero();
          row=alpha[n_dof] * mapping*reference_phi[n_dof][qp];
          result[n_dof][n_comp][qp].row(n_comp,row);
        }        
      }
     }
   return result;
  };

};


template<Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents>
class BaseShapeFunctionOperator<DivergenceOperator,Ndofs,Dim,ManifoldDim,ShapeFunctionDim,NComponents>
{public:

  virtual ~BaseShapeFunctionOperator(){};

  virtual const Vector<Real,Ndofs>
  phi_single_component(const Vector<Real,Dim>& qp_point )=0;


  template<Integer NQPpoints>
  const Vector<Vector< Real, NQPpoints>,Ndofs>
  phi(const Matrix<Real,NQPpoints,Dim>& qp_points)
  {
        Vector<Vector< Real, NQPpoints>,Ndofs> div;
        Vector<Real,Dim> qp_point;

        for(Integer qp=0;qp<NQPpoints;qp++)
        {
          qp_points.get_row(qp,qp_point);
          auto div_phi_single=phi_single_component(qp_point);
              for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
                   div[n_dof][qp]=div_phi_single[n_dof];
        }
        return div;};

  template<Integer NQPpoints,typename Mapping>
  const Vector<Vector<Vector< Matrix<Real,NComponents,1>, NQPpoints>,NComponents>,Ndofs>
  phiN(const Vector<Vector< Real, NQPpoints>,Ndofs>& divphi_reference,
       const Mapping& mapping,
       const Vector<Real,Ndofs> &alpha=1.0)
  {
        Vector<Vector<Vector< Matrix<Real,NComponents,1>, NQPpoints>,NComponents>,Ndofs> divphi;
        for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
         for(Integer n_comp=0;n_comp<NComponents;n_comp++)
          for(Integer qp=0;qp<NQPpoints;qp++)
            { divphi[n_dof][n_comp][qp].zero();
              divphi[n_dof][n_comp][qp](n_comp,0)=alpha[n_dof] * mapping * divphi_reference[n_dof][qp];
            } 
    return divphi;};


};







template<Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents>
class BaseShapeFunctionOperator<GradientOperator,Ndofs,Dim,ManifoldDim,ShapeFunctionDim,NComponents>
{public:

  virtual ~BaseShapeFunctionOperator(){};

 virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
 phi_single_component(const Vector<Real,Dim>& point)=0;

  template<Integer NQPpoints>
  const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPpoints>,NComponents>,Ndofs>
  phi(const Matrix<Real,NQPpoints,Dim>& qp_points)
  {
   Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPpoints>,NComponents>,Ndofs> f;
   Vector<Real,Dim> qp_point;

   for(Integer qp=0;qp<NQPpoints;qp++)
    {
    qp_points.get_row(qp,qp_point);
    auto grad_phi_single=phi_single_component(qp_point);
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
  const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPpoints>,NComponents>,Ndofs>
  phiN(const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPpoints>,NComponents>,Ndofs>& reference_grad_phi,
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


};




template<typename Operator, typename Elem,typename BaseFunctionSpace>
class ShapeFunctionOperator;



template<Integer NComponents_>
class ShapeFunctionOperator<IdentityOperator,Simplex<2,2>, Lagrange1<NComponents_> > : 
public BaseShapeFunctionOperator<IdentityOperator,3,2,2,1,NComponents_>
{
public:
  virtual const Matrix<Real,3,1>
  phi_single_component(const Vector<Real,2>& point)
       {const auto& xi=point[0];
        const auto& eta=point[1];
        const Real zeta = 1. - xi - eta;
        const Matrix<Real,3,1> shapefunction{zeta,xi,eta};
        return shapefunction;};
};


// template<Integer NComponents>
// class ShapeFunctionOperator<DivergenceOperator,Simplex<2,2>, Lagrange1<NComponents> > : 
// public BaseShapeFunctionOperator<DivergenceOperator,3,2,2,1,NComponents>
// {
// public:
//   virtual const Vector<DivPhiOperator<2,2,1,NComponents>,3>
//   phi_single_component(const Vector<Real,2>& qp_point )
//   {
//   Vector<DivPhiType<2,2,1,NComponents>,3> tryme; 

//     tryme[0](0,0)=1.0;
//     tryme[1](0,0)=1.0;
//     tryme[2](0,0)=-2.0;

//   return tryme;};

// };

template<Integer NComponents>
class ShapeFunctionOperator<GradientOperator,Simplex<2,2>, Lagrange1<NComponents> > : 
public BaseShapeFunctionOperator<GradientOperator,3,2,2,1,NComponents>
{
public:

  virtual const Vector< Matrix<Real,1,2>, 3>
  phi_single_component(const Vector<Real,2>& point)
  {
    Vector< Matrix<Real,1,2>, 3> grad{{-1,-1},{1,0},{0,1}};
    return grad;
  }

};












template<Integer NComponents_>
class ShapeFunctionOperator<IdentityOperator, Simplex<2,2>, Lagrange2<NComponents_> > : 
public BaseShapeFunctionOperator<IdentityOperator,6,2,2,1,NComponents_>
{
public:
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
};

// template<Integer NComponents>
// class ShapeFunctionOperator<DivergenceOperator, Simplex<2,2>, Lagrange2<NComponents> > : 
// public BaseShapeFunctionOperator<DivergenceOperator,6,2,2,1,NComponents>
// {
// public:
//   virtual const Vector<DivPhiType<2,2,1,NComponents>,6>
//   phi_single_component(const Vector<Real,2>& qp_point )
//   {
//     //////////////// E' SBAGLIATO, CI VUOLE DIPENDENZA DALLE VARIABILI ////////////////
//     Vector<DivPhiType<2,2,1,NComponents>,6> tryme; 
//     tryme[0](0,0)=(4.0+4.0);
//     tryme[1](0,0)=(4.0);
//     tryme[2](0,0)=(4.0);
//     tryme[3](0,0)=(-8.0);
//     tryme[4](0,0)=(4.0);
//     tryme[5](0,0)=(-8.0);
//     return tryme;
//  };
// };


template<Integer NComponents>
class ShapeFunctionOperator<GradientOperator, Simplex<2,2>, Lagrange2<NComponents> > : 
public BaseShapeFunctionOperator<GradientOperator,6,2,2,1,NComponents>
{
public:

  virtual const Vector< Matrix<Real,1,2>, 6>
  phi_single_component(const Vector<Real,2>& point)
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

};





template<Integer NComponents_>
class ShapeFunctionOperator<IdentityOperator, Simplex<2,2>, RT0<NComponents_> > : 
public BaseShapeFunctionOperator<IdentityOperator,3,2,2,2,NComponents_>
{
public:
  virtual const  Matrix<Real,3,2>
  phi_single_component(const Vector<Real,2>& point)
       {const auto& xi=point[0];
        const auto& eta=point[1];
        const Matrix<Real,3,2> shapefunction{xi,eta-1,xi-1,eta,xi,eta};
        return shapefunction;};
};


template<Integer NComponents>
class ShapeFunctionOperator<DivergenceOperator, Simplex<2,2>, RT0<NComponents> > : 
public BaseShapeFunctionOperator<DivergenceOperator,3,2,2,2,NComponents>
{
public:

  virtual const Vector<Real,3>
  phi_single_component(const Vector<Real,2>& qp_point )
  {
   Vector<Real,3> tryme;
     tryme[0]=2;
     tryme[1]=2;
     tryme[2]=2;    
   return tryme;
  };
};

//  virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
//  grad_phi_single_component(const Vector<Real,2>& point)
//        {
//         //static_assert(Dim==ShapeFunctionDim, "Grad(RT): shape function dim and space dim must be the same")
//         Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs> grad{ {1,0, 0,1},{1,0, 0,1},{1,0, 0,1}} ;
//         return grad;
//         assert("gradient is not defined for RT elements");};
// };




















































template<Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents=1>
class BaseShapeFunction
{public:

  virtual ~BaseShapeFunction(){};

  virtual const Matrix<Real,Ndofs,ShapeFunctionDim>
  phi_single_component(const Vector<Real,Dim>& point )=0;


  virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
  grad_phi_single_component(const Vector<Real,Dim>& point)=0;


  virtual const Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs>
  div_phi_single_component(const Vector<Real,Dim>& qp_point )=0;

  template<Integer NQPpoints>
  const Vector<Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>, NQPpoints>,Ndofs>
  div_phi(const Matrix<Real,NQPpoints,Dim>& qp_points)
  {
        Vector<Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>, NQPpoints>,Ndofs> div;
        Vector<Real,Dim> qp_point;
        std::cout<<"div_phi Ndofs=="<<Ndofs<<std::endl;
        std::cout<<"div_phi NQPpoints=="<<NQPpoints<<std::endl;

        for(Integer qp=0;qp<NQPpoints;qp++)
        {
          qp_points.get_row(qp,qp_point);
          auto div_phi_single=div_phi_single_component(qp_point);
              for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
                div[n_dof][qp]=div_phi_single[n_dof];
        }

        return div;};

  template<Integer NQPpoints,typename Mapping>
  const Vector<Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>, NQPpoints>,Ndofs>
  div_phiN(const Vector<Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>, NQPpoints>,Ndofs>& divphi_reference,
           const Mapping& mapping,
           const Vector<Real,Ndofs> &alpha=1.0)
  {
        Vector<Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>, NQPpoints>,Ndofs> divphi;
        for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
         for(Integer qp=0;qp<NQPpoints;qp++)
            {
              divphi[n_dof][qp]=alpha[n_dof] * mapping * divphi_reference[n_dof][qp];
            }


    return divphi;};



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

  template<Integer NQPpoints>
  const Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPpoints>,Ndofs>
  phi(const Matrix<Real,NQPpoints,Dim>& qp_points)
  {
   Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPpoints>,Ndofs> f;
   //Vector<Real,ShapeFunctionDim> row_tmp;
   Vector<Real,Dim> qp_point;

   for(Integer qp=0;qp<NQPpoints;qp++)
    {
    qp_points.get_row(qp,qp_point);
    auto phi_single=phi_single_component(qp_point);
    for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      //phi_single.get_row(n_dof,row_tmp);
      // for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      // {
        //f[n_dof][n_comp][qp].zero();
      for(Integer ii=0;ii<ShapeFunctionDim;ii++)
        f[n_dof][qp][ii]=phi_single(n_dof,ii);
      
     }
    }

   return f;
  };

  template<Integer NQPpoints, typename Mapping>
  const Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPpoints>,NComponents>,Ndofs>
  phiN(const Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPpoints>,Ndofs>& reference_phi,
       const Mapping& mapping, 
       const Vector<Real,Ndofs> &alpha=1.0)
  {
   Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPpoints>,NComponents>,Ndofs> result;
   Vector<Real,ShapeFunctionDim> row;
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
     {
      for(Integer n_comp=0;n_comp<NComponents;n_comp++)
      {
        for(Integer qp=0;qp<NQPpoints;qp++)
        { 
          result[n_dof][n_comp][qp].zero();
          row=alpha[n_dof] * mapping*reference_phi[n_dof][qp];
          result[n_dof][n_comp][qp].row(n_comp,row);
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


};



























template<typename Elem,typename BaseFunctionSpace>
class ShapeFunction;








template<Integer NComponents_>
class ShapeFunction<Simplex<2,2>, Lagrange1<NComponents_> > : public BaseShapeFunction<3,2,2,1,NComponents_>
{
public:
  static constexpr Integer Ndofs=3;
  static constexpr Integer Dim=2;
  static constexpr Integer ManifoldDim=2;
  static constexpr Integer ShapeFunctionDim=1;
  static constexpr Integer NComponents=NComponents_;


  virtual const Matrix<Real,3,1>
  phi_single_component(const Vector<Real,2>& point)
       {const auto& xi=point[0];
        const auto& eta=point[1];
        const Real zeta = 1. - xi - eta;
        const Matrix<Real,3,1> shapefunction{zeta,xi,eta};
        return shapefunction;};



  virtual const Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs>
  div_phi_single_component(const Vector<Real,Dim>& qp_point )
  {
    assert(NComponents==ManifoldDim && "divergence for non vector shape functions requires: NComponents==ManifoldDim");
    assert(NComponents==Dim && "divergence for non vector shape functions requires: NComponents==Dim");
  Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs> tryme; 

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


  virtual const Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs>
  div_phi_single_component(const Vector<Real,Dim>& qp_point )
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
  div_phi_single_component(const Vector<Real,Dim>& qp_point )
  {
   Vector<DivPhiType<Dim,ManifoldDim,ShapeFunctionDim,NComponents>,Ndofs> tryme;
   for(Integer n_comp=0;n_comp<NComponents;n_comp++)
   {
     tryme[0](n_comp,0)=2;
     tryme[1](n_comp,0)=2;
     tryme[2](n_comp,0)=2;    
   }
   return tryme;
  };


 virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
 grad_phi_single_component(const Vector<Real,2>& point)
       {
        //static_assert(Dim==ShapeFunctionDim, "Grad(RT): shape function dim and space dim must be the same")
        Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs> grad{ {1,0, 0,1},{1,0, 0,1},{1,0, 0,1}} ;
        return grad;
        assert("gradient is not defined for RT elements");};
};

























template<typename Operator,typename Elem,Integer FEFamily, Integer NComponents=1>
class MapFromReference;



template<Integer Dim, Integer ManifoldDim,Integer NComponents>
class MapFromReference<IdentityOperator, Simplex<Dim,ManifoldDim>,LagrangeFE,NComponents>
{
 private:
 public:
        inline const Real  operator()(const Matrix<Real, Dim, ManifoldDim>& J){return 1.0;}
};

template<Integer Dim, Integer ManifoldDim,Integer NComponents>
class MapFromReference<GradientOperator, Simplex<Dim,ManifoldDim>,LagrangeFE,NComponents>
{
 private:
 public:
        inline const Matrix<Real, Dim, ManifoldDim>  operator() (const Matrix<Real, Dim, ManifoldDim>& J)
         {static_assert(Dim==ManifoldDim,"Dim==ManifoldDim for inverting jacobian");
          auto inv=  inverse(J);
          return inv;}
};



template<Integer Dim, Integer ManifoldDim,Integer NComponents>
class MapFromReference<IdentityOperator, Simplex<Dim,ManifoldDim>,RaviartThomasFE,NComponents>
{
 private:
 public:
        inline const Matrix<Real, Dim, ManifoldDim> operator()(const Matrix<Real, Dim, ManifoldDim>& J)
         {Matrix<Real, Dim, ManifoldDim> mapping=J;
          mapping/=det(mapping);
          return mapping;}
};

template<Integer Dim, Integer ManifoldDim,Integer NComponents>
class MapFromReference<DivergenceOperator, Simplex<Dim,ManifoldDim>,RaviartThomasFE,NComponents>
{
 private:
 public:
       inline const Real  operator() (const Matrix<Real, Dim, ManifoldDim>& J){return 1.0/det(J);}
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





template<Integer NQPpoints, typename Operation1,typename Operation2, typename MatTrial,typename PhiTrial,typename MatTest,typename PhiTest> 
Vector<Real, NQPpoints> BilinearFormQPvalues(const MatTrial& A, const PhiTrial& phi_trial,const MatTest& B, const PhiTest& phi_test)
{
 Vector<Real, NQPpoints> result;
 // auto left=A*phi_trial;
 // auto right=B*phi_test;

Operation1 prod1;
Operation2 prod2;
Contraction contr;
 for(Integer qp=0;qp<NQPpoints;qp++)
 {
  // std::cout<<"A"<<std::endl;
  // A.describe(std::cout);
  // std::cout<<"B"<<std::endl;
  // B.describe(std::cout);
  // std::cout<<"phi_trial"<<std::endl;
  // phi_trial.describe(std::cout);
  // std::cout<<"phi_test"<<std::endl;
  // phi_test.describe(std::cout);
  // std::cout<<"phi_test[qp]"<<std::endl;
  // std::cout<<phi_test[qp]<<std::endl;
  auto e1=prod1(A,phi_trial[qp]);
  auto e2=prod2(B,phi_test[qp]);
  // std::cout<<"left"<<std::endl;
  // e1.describe(std::cout);
  // std::cout<<"right"<<std::endl;
  //std::cout<<e2<<std::endl;
  // e2.describe(std::cout);
  auto e3=contr(e1,e2);
  // std::cout<<"result"<<std::endl;
  //std::cout<<e3<<std::endl;
  result[qp]=e3;
 }

 return result;

}


template <Integer QPOrder, typename TrialFunction,typename TestFunction, typename OperatorTrial, typename OperatorTest,typename MeshT >
void BilinearFormIntegrator(const MeshT& mesh)
{

using trial_function=Elem2FunctionSpace<TrialFunction>;
using test_function=Elem2FunctionSpace<TestFunction>;

using Elem=typename TrialFunction::Elem;
constexpr Integer Dim=Elem::Dim;
constexpr Integer ManifoldDim=Elem::ManifoldDim;  

// trial function infos
constexpr Integer FEFamily_trial=TrialFunction::FEFamily;
constexpr Integer NComponents_trial=TrialFunction::NComponents;
constexpr Integer ShapeFunctionDim_trial=ShapeFunction<Elem, trial_function>::ShapeFunctionDim;
constexpr Integer Ndofs_trial=ShapeFunction<Elem, trial_function>::Ndofs;
ShapeFunctionOperator<OperatorTrial, Elem, trial_function> trial;
MapFromReference<OperatorTrial,Simplex<Dim,ManifoldDim>,FEFamily_trial,NComponents_trial> map_trial;

// test function infos
constexpr Integer FEFamily_test=TestFunction::FEFamily;
constexpr Integer NComponents_test=TestFunction::NComponents;
constexpr Integer ShapeFunctionDim_test=ShapeFunction<Elem, test_function>::ShapeFunctionDim;
constexpr Integer Ndofs_test=ShapeFunction<Elem, test_function>::Ndofs;
ShapeFunctionOperator<OperatorTest,Elem, test_function> test;
MapFromReference<OperatorTest,Simplex<Dim,ManifoldDim>,FEFamily_test,NComponents_test> map_test;


auto sn=SignedNormal<Simplex<Dim,ManifoldDim>>(mesh);
auto normal_sign= sn.sign(); 
const auto& n_elements=mesh.n_elements();




// quadrature rule
constexpr Integer NQPpoints=GaussPoints<Elem,QPOrder>::NQPpoints; 
GaussPoints<Elem,QPOrder> gauss;
auto qp_points=gauss.qp_points();
auto qp_weights=gauss.qp_weights();

std::vector<Vector<Real,Dim>> points(ManifoldDim+1);
Matrix<Real, Dim, ManifoldDim> J;
Real volume;
Simplex<Dim,ManifoldDim> simplex;


for(Integer ii=0;ii<ManifoldDim+1;ii++)
   simplex.nodes[ii]=ii;


Vector<Real, Ndofs_trial > alpha_trial=1;
Vector<Real, Ndofs_test > alpha_test=1;


auto reference_trial=trial.phi(qp_points);
auto reference_test=test.phi(qp_points);


Matrix<Real, NComponents_trial,ShapeFunctionDim_trial*NComponents_trial> A=2.0;
Matrix<Real, 2,NComponents_test> B=10.0;

Matrix<Real, Ndofs_test * NComponents_test, Ndofs_trial * NComponents_trial > mat;


for(Integer elem_iter=0 ; elem_iter < n_elements ; elem_iter++)
{
  auto elemnodes_global=mesh.elem(elem_iter).nodes;
  for(Integer mm=0;mm<points.size();mm++)
     points[mm]=mesh.point(elemnodes_global[mm]);
  jacobian(simplex,points,J);
  volume=unsigned_volume(simplex,points);
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
         const auto& vec=BilinearFormQPvalues<NQPpoints,Multiply,Multiply>(A,element_trial[n_dof_trial][n_comp_trial],B,element_test[n_dof_test][n_comp_test]);
         mat(i,j)=dot(vec,qp_weights)*volume;
        }

}

};


}



#endif