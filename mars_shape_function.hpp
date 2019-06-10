




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

 template<typename T,Integer Rows,Integer Cols>
 Matrix<T,1,Rows> operator()(const Matrix<T, Rows,Cols> &A,const Matrix<T, 1,Cols> &B)
 {
       Matrix<T,1,Rows> result=0;
       for(Integer i = 0; i < Rows; ++i) 
         for(Integer j = 0; j < Cols; ++j) 
             result(0,i) += A(i, j) * B(0,j);
       return result;
 }
 template<typename T,Integer Rows,Integer Cols>
 Matrix<T,Rows,1> operator()(const Matrix<T, Rows,Cols> &A,const Matrix<T, Cols,1> &B)
 {
       Matrix<T,Rows,1> result=0;
       for(Integer i = 0; i < Rows; ++i) 
         for(Integer j = 0; j < Cols; ++j) 
             result(i,0) += A(i, j) * B(j,0);
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

 Real operator()(const Real &A,const Real &B)
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





template<typename T,Integer NQPoints,Integer NComponents>
class FQPValues;


template<typename T,Integer NQPoints>
class QPValues
{
public:


  using type=  Vector<T,NQPoints>;
  using FunctionType=  Vector<T,NQPoints>;
  using GradientType=  Vector<T,NQPoints>;
  using DivergenceType=  Vector<T,NQPoints>;
  using CurlType=  Vector<T,NQPoints>;
  QPValues(): qpvalues_(){};
  QPValues(const type& v): qpvalues_(v){};
  inline type operator()()const {return qpvalues_;};

  inline T &operator[](const Integer i)
      {
          assert(i < NQPoints);
          return qpvalues_[i];
      }
  inline const T &operator[](const Integer i)const
      {
          assert(i < NQPoints);
          return qpvalues_[i];
      }



private:
  type qpvalues_;

};




template<Integer Rows,Integer Cols,Integer NQPoints>
class QPValues<Matrix<Real,Rows,Cols>,NQPoints>
{
public:


  using T= Matrix<Real,Rows,Cols>;
  using type=  Vector<T,NQPoints>;
  using FunctionType=  Vector<T,NQPoints>;
  using GradientType=  Vector<T,NQPoints>;
  using DivergenceType=  Vector<T,NQPoints>;
  using CurlType=  Vector<T,NQPoints>;
  QPValues(): qpvalues_(){};
  QPValues(const type& v): qpvalues_(v){};
  inline type operator()()const {return qpvalues_;};

  inline T &operator[](const Integer i)
      {
          assert(i < NQPoints);
          return qpvalues_[i];
      }
  inline const T &operator[](const Integer i)const
      {
          assert(i < NQPoints);
          return qpvalues_[i];
      }

      template<Integer OtherDim, Integer NComponents>
      FQPValues<Matrix<Real,Rows,OtherDim>,NQPoints,NComponents> operator*( 
                                 const FQPValues<Matrix<Real,Cols,OtherDim>,NQPoints,NComponents>& fqpvalue)
      { 
        FQPValues<Matrix<Real,Rows,OtherDim>,NQPoints,NComponents> result;
        for(Integer i = 0; i < NComponents; ++i) 
        {
          for(Integer j = 0; j < NQPoints; ++j)
            {
            result[i][j]= qpvalues_[i]*fqpvalue[i][j];
            }          
        }       
        return result;
      };

      // tensor contraction tensor: Real= contract(M<Rows,Cols>, M<Rows,Cols>),contract(Vec<Rows>, Vec<Rows>),contract(Real,Real)
      template<Integer OtherCols>
      QPValues<Matrix<Real,Rows,OtherCols>,NQPoints> operator*(const QPValues<Matrix<Real,Cols,OtherCols>,NQPoints>& other)
      { 
        QPValues<Matrix<Real,Rows,OtherCols>,NQPoints> result;
        Contraction contract;
          for(Integer j = 0; j < NQPoints; ++j)
            {
            result[j]= contract(qpvalues_[j],other[j]);
            }               
        return result;
      };


      void describe(std::ostream &os) const
      {
          for(Integer i = 0; i < NQPoints; ++i) {
              os << (*this)[i] << " ";
          }

          os << "\n";
      }

      friend std::ostream &operator<<(std::ostream &os, const QPValues &v)
      {
          v.describe(os);
          return os;
      }

private:
  type qpvalues_;

};





// template<typename T,Integer NQPoints,Integer NComponents>
// class FQPValuesBase
// {
// public:

//       using type= Vector<Vector<T,NQPoints>,NComponents>;

//       FQPValuesBase(): fqpvalues_(){};
//       FQPValuesBase(const type& v): fqpvalues_(v){};
//       inline type operator()(){return fqpvalues_;};
//       inline Vector<T,NQPoints> operator()(const Integer& i)const{return fqpvalues_[i];};
//       inline void operator()(const Integer& i,const Vector<T,NQPoints>& u){fqpvalues_[i]=u;};
      
//       inline Vector<T,NQPoints> &operator[](const Integer i)
//       {
//           assert(i < NComponents);
//           return fqpvalues_[i];
//       }
//       inline const Vector<T,NQPoints> &operator[](const Integer i)const
//       {
//           assert(i < NComponents);
//           return fqpvalues_[i];
//       }
// protected:
//       type fqpvalues_;
//       Vector<T,NQPoints> tmp1_;

// };






// FQPVALUES REAL
template<Integer NQPoints,Integer NComponents>
class FQPValues<Real,NQPoints,NComponents> //: public FQPValuesBase<Real,NQPoints,NComponents>
{
public:
  using T=Real;
  using type= Vector<Vector<T,NQPoints>,NComponents>;
      FQPValues(): fqpvalues_(){};
      FQPValues(const type& v): fqpvalues_(v){};

      // FQPValues(): FQPValuesBase<T,NQPoints,NComponents>(){};
      // FQPValues(const type& v): FQPValuesBase<T,NQPoints,NComponents>(v){};
      inline type operator()()const{return fqpvalues_;};
      inline Vector<T,NQPoints> operator()(const Integer& i)const{return fqpvalues_[i];};
      inline void operator()(const Integer& i,const Vector<T,NQPoints>& u){fqpvalues_[i]=u;};   

      inline Vector<T,NQPoints> &operator[](const Integer i)
      {
          assert(i < NComponents);
          return fqpvalues_[i];
      }
      inline const Vector<T,NQPoints> &operator[](const Integer i)const
      {
          assert(i < NComponents);
          return fqpvalues_[i];
      }

      // equal
      inline FQPValues& operator = (const FQPValues &u)
      {            
        (*this).fqpvalues_=u();
        return *this;
      } 

      // unary add
      FQPValues operator+()
      {
        FQPValues result;
        result.fqpvalues_=+(*this).fqpvalues_;
        
        return result;
      };
      // binary add
      FQPValues operator+(const FQPValues &fqpvalue) const
      {
        FQPValues result;
        result.fqpvalues_=(*this).fqpvalues_+fqpvalue.fqpvalues_;
        
        return result;
      };

      // unary minus
      FQPValues operator-()
      {
        FQPValues result;
        result.fqpvalues_=-(*this).fqpvalues_;
        
        return result;
      };
      // binary minus
      FQPValues operator-(const FQPValues &fqpvalue) const
      {
        FQPValues result;
        result.fqpvalues_=(*this).fqpvalues_-fqpvalue.fqpvalues_;
        
        return result;
      };

      // right scalar multiply
      FQPValues operator*(const Real &value) const
      { 
        FQPValues result;
        result.fqpvalues_=(*this).fqpvalues_*value;
        return result;
      };
      // left scalar multiply
      friend FQPValues operator*(const Real &value, const FQPValues& F) 
      { 
        FQPValues result;
        result=value*F;
        return result;
      };

      // tensor contraction tensor: Real= contract(M<Rows,Cols>, M<Rows,Cols>),contract(Vec<Rows>, Vec<Rows>),contract(Real,Real)
      template<typename T>
      friend FQPValues operator*(const QPValues<T,NQPoints> &qpvalue,const FQPValues<T,NQPoints,NComponents>& fqpvalue)
      { 
        FQPValues result;
        Contraction contract;
        const auto& tmp=qpvalue();
        for(Integer i = 0; i < NComponents; ++i) 
        {
          for(Integer j = 0; j < NQPoints; ++j)
            {
            result[i][j]= contract(tmp[i],fqpvalue[i][j]);
            }          
        }       
        return result;
      };

      void describe(std::ostream &os) const
      {
        for(Integer n_comp=0;n_comp<NComponents;n_comp++)
          for(Integer qp = 0; qp < NQPoints; ++qp) {
              os << (*this)[n_comp][qp] << " ";
          }

          os << "\n";
      }

      friend std::ostream &operator<<(std::ostream &os, const FQPValues &v)
      {
          v.describe(os);
          return os;
      }

protected:
      type fqpvalues_;
      Vector<T,NQPoints> tmp1_;
};

// FQPVALUES VECTOR
template<Integer Rows,Integer NQPoints,Integer NComponents>
class FQPValues<Vector<Real,Rows>,NQPoints,NComponents> //: public FQPValuesBase<Vector<Real,Rows>,NQPoints,NComponents>
{
public:
  using T=Vector<Real,Rows>;
  using type= Vector<Vector<T,NQPoints>,NComponents>;

       FQPValues(): fqpvalues_(){};
      FQPValues(const type& v): fqpvalues_(v){};
      
      // FQPValues(): FQPValuesBase<T,NQPoints,NComponents>(){};
      // FQPValues(const type& v): FQPValuesBase<T,NQPoints,NComponents>(v){};
      inline type operator()()const{return fqpvalues_;};
      inline Vector<T,NQPoints> operator()(const Integer& i)const{return fqpvalues_[i];};
      inline void operator()(const Integer& i,const Vector<T,NQPoints>& u){fqpvalues_[i]=u;};     
      inline Vector<T,NQPoints> &operator[](const Integer i)
      {
          assert(i < NComponents);
          return fqpvalues_[i];
      }
      inline const Vector<T,NQPoints> &operator[](const Integer i)const
      {
          assert(i < NComponents);
          return fqpvalues_[i];
      }

            // equal
      inline FQPValues& operator = (const FQPValues &u)
      {            
        (*this).fqpvalues_=u();
        return *this;
      } 

      // unary add
      FQPValues operator+()
      {
        FQPValues result;
        result.fqpvalues_=+(*this).fqpvalues_;
        
        return result;
      };
      // binary add
      FQPValues operator+(const FQPValues &fqpvalue) const
      {
        FQPValues result;
        result.fqpvalues_=(*this).fqpvalues_+fqpvalue.fqpvalues_;
        
        return result;
      };

      // unary minus
      FQPValues operator-()
      {
        FQPValues result;
        result.fqpvalues_=-(*this).fqpvalues_;
        
        return result;
      };
      // binary minus
      FQPValues operator-(const FQPValues &fqpvalue) const
      {
        FQPValues result;
        result.fqpvalues_=(*this).fqpvalues_-fqpvalue.fqpvalues_;
        
        return result;
      };

      // right scalar multiply
      FQPValues operator*(const Real &value) const
      { 
        FQPValues result;
        result.fqpvalues_=(*this).fqpvalues_*value;
        return result;
      };
      // left scalar multiply
      friend FQPValues operator*(const Real &value, const FQPValues& F) 
      { 
        FQPValues result;
        result=value*F;
        return result;
      };

      // matrix times matrix multiply: M<Rows,Cols>= M<Rows,Cols2> * M<Rows2,Cols>
      template<Integer Rows2,Integer Cols2>
      friend FQPValues operator*(const QPValues<Matrix<Real,Rows2,Cols2>,NQPoints> &qpvalue,const FQPValues<T,NQPoints,NComponents>& fqpvalue)
      { 
        assert(Cols2==Rows);
        FQPValues result;
        const auto& tmp=qpvalue();
        for(Integer i = 0; i < NComponents; ++i) 
        {
          for(Integer j = 0; j < NQPoints; ++j)
            {
            result[i][j]= tmp[i]*fqpvalue[i][j];
            }          
        }       
        return result;
      };

      void describe(std::ostream &os) const
      {
        for(Integer n_comp=0;n_comp<NComponents;n_comp++)
          for(Integer qp = 0; qp < NQPoints; ++qp) {
              os << (*this)[n_comp][qp] << " ";
          }

          os << "\n";
      }

      friend std::ostream &operator<<(std::ostream &os, const FQPValues &v)
      {
          v.describe(os);
          return os;
      }

protected:
      type fqpvalues_;
      Vector<T,NQPoints> tmp1_;
};

// FQPVALUES MATRIX
template<Integer Rows,Integer Cols,Integer NQPoints,Integer NComponents>
class FQPValues<Matrix<Real,Rows,Cols>,NQPoints,NComponents> //: public FQPValuesBase<Matrix<Real,Rows,Cols>,NQPoints,NComponents>
{

public:
  using T=Matrix<Real,Rows,Cols>;
  using type= Vector<Vector<T,NQPoints>,NComponents>;

       FQPValues(): fqpvalues_(){};
      FQPValues(const type& v): fqpvalues_(v){};
      
      // FQPValues(): FQPValuesBase<T,NQPoints,NComponents>(){};
      // FQPValues(const type& v): FQPValuesBase<T,NQPoints,NComponents>(v){};
      inline type operator()()const{return fqpvalues_;};
      inline Vector<T,NQPoints> operator()(const Integer& i)const{return fqpvalues_[i];};
      inline void operator()(const Integer& i,const Vector<T,NQPoints>& u){fqpvalues_[i]=u;};     
       inline Vector<T,NQPoints> &operator[](const Integer i)
      {
          assert(i < NComponents);
          return fqpvalues_[i];
      }
      inline const Vector<T,NQPoints> &operator[](const Integer i)const
      {
          assert(i < NComponents);
          return fqpvalues_[i];
      }

           // equal
      inline FQPValues& operator = (const FQPValues &u)
      {            
        (*this).fqpvalues_=u();
        return *this;
      } 

      // unary add
      FQPValues operator+()
      {
        FQPValues result;
        result.fqpvalues_=+(*this).fqpvalues_;
        
        return result;
      };
      // binary add
      FQPValues operator+(const FQPValues &fqpvalue) const
      {
        FQPValues result;
        result.fqpvalues_=(*this).fqpvalues_+fqpvalue.fqpvalues_;
        
        return result;
      };

      // unary minus
      FQPValues operator-()
      {
        FQPValues result;
        result.fqpvalues_=-(*this).fqpvalues_;
        
        return result;
      };
      // binary minus
      FQPValues operator-(const FQPValues &fqpvalue) const
      {
        FQPValues result;
        result.fqpvalues_=(*this).fqpvalues_-fqpvalue.fqpvalues_;
        
        return result;
      };

      // right scalar multiply
      FQPValues operator*(const Real &value) const
      { 
        FQPValues result;
        result.fqpvalues_=(*this).fqpvalues_*value;
        return result;
      };
      // left scalar multiply
      friend FQPValues operator*(const Real &value, const FQPValues& F) 
      { 
        FQPValues result;
        result=F*value;
        return result;
      };

      // matrix times matrix multiply: M<Rows,Cols>= M<Rows,Cols2> * M<Rows2,Cols>
      friend FQPValues operator*(const QPValues<Real,NQPoints> &qpvalue, const FQPValues<Matrix<Real,Rows,Cols>,NQPoints,NComponents>& fqpvalue)
      { 
        Contraction contract;
        FQPValues result;
        const auto& tmp=qpvalue();
        for(Integer i = 0; i < NComponents; ++i) 
        {
          for(Integer j = 0; j < NQPoints; ++j)
            {
            result[i][j]= tmp[i]*fqpvalue[i][j];
            }          
        }       
        return result;
      };

      void describe(std::ostream &os) const
      {
        for(Integer n_comp=0;n_comp<NComponents;n_comp++)
          for(Integer qp = 0; qp < NQPoints; ++qp) {
              os << (*this)[n_comp][qp] << " ";
          }

          os << "\n";
      };

      friend std::ostream &operator<<(std::ostream &os, const FQPValues &v)
      {
          v.describe(os);
          return os;
      };


protected:
      type fqpvalues_;
      Vector<T,NQPoints> tmp1_;

};












template<typename Operation>
class Product{
 public:
 template<typename Input1,typename Input2, typename Output>
 Output operator()(Input1 i1,Input2 i2){return Operation(i1,i2);};
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















template<typename Elem_,Integer Order_>
class BaseQuadrature
{
 using Elem=Elem_;
 static constexpr Integer Order=Order;
};

template<typename Elem,Integer Order>
class GaussPoints;

template< Integer Dim >
class GaussPoints< Simplex<Dim,2> , 1>:BaseQuadrature<Simplex<Dim,2>,1>
{
private:
 Matrix<Real,1,2> qp_points_;
 Vector<Real,1> qp_weights_;
public:
  using qp_points_type=Matrix<Real,1,2>;
  static constexpr Integer NQPoints=1;
  const Matrix<Real,1,2> qp_points()const {return qp_points_;};
  const Vector<Real,1> qp_weights()const {return qp_weights_;};

  GaussPoints<Simplex<Dim,2>,1>():
  qp_points_({0.33333333333333, 0.33333333333333}),
  qp_weights_({1})
   {}; 

};



template<Integer Dim>
class GaussPoints<Simplex<Dim,2>,2>:BaseQuadrature<Simplex<Dim,2>,2>
{
private:
  Matrix<Real,3,2>  qp_points_;
  Vector<Real,3> qp_weights_;
public:
  using qp_points_type=Matrix<Real,3,2>;
  static constexpr Integer NQPoints=3;
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
class GaussPoints<Simplex<Dim,2>,3>:BaseQuadrature<Simplex<Dim,2>,3>
{
private:
  Matrix<Real,4,2>  qp_points_;
  Vector<Real,4> qp_weights_;
public:
  using qp_points_type=Matrix<Real,4,2>;
  static constexpr Integer NQPoints=4;
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
class GaussPoints<Simplex<Dim,2>,4>:BaseQuadrature<Simplex<Dim,2>,4>
{
private:
  Matrix<Real,6,2>  qp_points_;
  Vector<Real,6> qp_weights_;
public:
  using qp_points_type=Matrix<Real,6,2>;
  static constexpr Integer NQPoints=6;
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
class GaussPoints<Simplex<Dim,2>,5>:BaseQuadrature<Simplex<Dim,2>,5>
{
private:
  Matrix<Real,7,2>  qp_points_;
  Vector<Real,7> qp_weights_;
public:
  using qp_points_type=decltype(qp_points_);
  static constexpr Integer NQPoints=7;
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












// Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents=1>
// template< typename Elem,typename BaseFunctionSpace>
// class BaseShapeFunction2;//: public ElemFunctionSpace<Elem,BaseFunctionSpace>;

// template<typename Elem, typename BaseFunctionSpace>
// class BaseShapeFunction2;

// template<Integer Dim, Integer ManifoldDim, typename BaseFunctionSpace>
// class BaseShapeFunction2<Simplex<Dim,ManifoldDim>,BaseFunctionSpace>:
// public ElementFunctionSpace<Simplex<Dim,ManifoldDim>,
//                             BaseFunctionSpace::FEFamily,BaseFunctionSpace::Order,BaseFunctionSpace::Continuity,BaseFunctionSpace::NComponents>
// {
// };


// template<Integer Dim, Integer ManifoldDim, Integer Continuity,Integer NComponents>
// class BaseShapeFunction2<Simplex<Dim,ManifoldDim>>:
// public ElementFunctionSpace<Simplex<Dim,ManifoldDim>,LagrangeFE,1,Continuity,NComponents>
// {
// };











// template<typename Operator, Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents=1>
// class BaseShapeFunctionOperator;








// template<Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents>
// class BaseShapeFunctionOperator<IdentityOperator,Ndofs,Dim,ManifoldDim,ShapeFunctionDim,NComponents>
// {public:

//   virtual ~BaseShapeFunctionOperator(){};

//   virtual const Matrix<Real,Ndofs,ShapeFunctionDim>
//   phi_single_component(const Vector<Real,Dim>& point )=0;

//   template<Integer NQPoints>
//   Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPoints>,Ndofs>
//   phi(const Matrix<Real,NQPoints,Dim>& qp_points)
//   {
//    Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPoints>,Ndofs> f;
//    //Vector<Real,ShapeFunctionDim> row_tmp;
//    Vector<Real,Dim> qp_point;

//    for(Integer qp=0;qp<NQPoints;qp++)
//     {
//     qp_points.get_row(qp,qp_point);
//     auto phi_single=phi_single_component(qp_point);
//     for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//      {
//       for(Integer ii=0;ii<ShapeFunctionDim;ii++)
//         f[n_dof][qp][ii]=phi_single(n_dof,ii);     
//      }
//     }
//    return f;
//   };

//   template<Integer NQPoints, typename Mapping>
//   Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPoints>,NComponents>,Ndofs>
//   phiN(const Vector<Vector<Vector<Real,ShapeFunctionDim>,NQPoints>,Ndofs>& reference_phi,
//        const Mapping& mapping, 
//        const Vector<Real,Ndofs> &alpha=1.0)
//   {
//    Vector<Vector<Vector<Matrix<Real,NComponents,ShapeFunctionDim>,NQPoints>,NComponents>,Ndofs> result;
//    Vector<Real,ShapeFunctionDim> row;
//   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//      {
//       for(Integer n_comp=0;n_comp<NComponents;n_comp++)
//       {
//         for(Integer qp=0;qp<NQPoints;qp++)
//         { 
//           result[n_dof][n_comp][qp].zero();
//           row=alpha[n_dof] * mapping*reference_phi[n_dof][qp];
//           result[n_dof][n_comp][qp].row(n_comp,row);
//         }        
//       }
//      }
//    return result;
//   };

// };


// template<Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents>
// class BaseShapeFunctionOperator<DivergenceOperator,Ndofs,Dim,ManifoldDim,ShapeFunctionDim,NComponents>
// {public:

//   virtual ~BaseShapeFunctionOperator(){};

//   virtual const Vector<Real,Ndofs>
//   phi_single_component(const Vector<Real,Dim>& qp_point )=0;


//   template<Integer NQPoints>
//   const Vector<Vector< Real, NQPoints>,Ndofs>
//   phi(const Matrix<Real,NQPoints,Dim>& qp_points)
//   {
//         Vector<Vector< Real, NQPoints>,Ndofs> div;
//         Vector<Real,Dim> qp_point;

//         for(Integer qp=0;qp<NQPoints;qp++)
//         {
//           qp_points.get_row(qp,qp_point);
//           auto div_phi_single=phi_single_component(qp_point);
//               for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//                    div[n_dof][qp]=div_phi_single[n_dof];
//         }
//         return div;};

//   template<Integer NQPoints,typename Mapping>
//   const Vector<Vector<Vector< Matrix<Real,NComponents,1>, NQPoints>,NComponents>,Ndofs>
//   phiN(const Vector<Vector< Real, NQPoints>,Ndofs>& divphi_reference,
//        const Mapping& mapping,
//        const Vector<Real,Ndofs> &alpha=1.0)
//   {
//         Vector<Vector<Vector< Matrix<Real,NComponents,1>, NQPoints>,NComponents>,Ndofs> divphi;
//         for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//          for(Integer n_comp=0;n_comp<NComponents;n_comp++)
//           for(Integer qp=0;qp<NQPoints;qp++)
//             { divphi[n_dof][n_comp][qp].zero();
//               divphi[n_dof][n_comp][qp](n_comp,0)=alpha[n_dof] * mapping * divphi_reference[n_dof][qp];
//             } 
//     return divphi;};


// };







// template<Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents>
// class BaseShapeFunctionOperator<GradientOperator,Ndofs,Dim,ManifoldDim,ShapeFunctionDim,NComponents>
// {public:

//   virtual ~BaseShapeFunctionOperator(){};

//  virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
//  phi_single_component(const Vector<Real,Dim>& point)=0;

//   template<Integer NQPoints>
//   const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs>
//   phi(const Matrix<Real,NQPoints,Dim>& qp_points)
//   {
//    Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs> f;
//    Vector<Real,Dim> qp_point;

//    for(Integer qp=0;qp<NQPoints;qp++)
//     {
//     qp_points.get_row(qp,qp_point);
//     auto grad_phi_single=phi_single_component(qp_point);
//     for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//      {
//       auto grad_phi=grad_phi_single[n_dof];
//       for(Integer n_comp=0;n_comp<NComponents;n_comp++)
//       {
//         f[n_dof][n_comp][qp].zero();
//         for(Integer ii=0;ii<ShapeFunctionDim;ii++)
//           {
//             for(Integer jj=0;jj<Dim;jj++)
//             {
//               f[n_dof][n_comp][qp](n_comp*ShapeFunctionDim+ii,jj)=grad_phi(ii,jj);
//             }
//           }
//       }
//      }
//     }

//    return f;
//   };

//   template<Integer NQPoints, typename Mapping>
//   const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs>
//   phiN(const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs>& reference_grad_phi,
//        const Mapping& mapping, 
//        const Vector<Real,Ndofs> &alpha=1.0)
//   {
//    Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs> result;
//    Vector<Real,Dim> row;
//     for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//      {
//       for(Integer n_comp=0;n_comp<NComponents;n_comp++)
//       {
//         for(Integer qp=0;qp<NQPoints;qp++)
//         {
//           result[n_dof][n_comp][qp].zero();
//           reference_grad_phi[n_dof][n_comp][qp].get_row(n_comp,row); 
//           // std::cout<<" firs row ="<<std::endl;
//           // reference_grad_phi[n_dof][n_comp][qp].describe(std::cout);     
//           // row.describe(std::cout);     
//           // std::cout<<" mapping ="<<std::endl;
//           // mapping.describe(std::cout);    
//           row=alpha[n_dof]*mapping*row;      
//           result[n_dof][n_comp][qp].row(n_comp,row);
//           // std::cout<<" result ="<<std::endl;
//           // result[n_dof][n_comp][qp].describe(std::cout);        
//         }
//       }
//      }
//    return result;
//   };


// };




// template<typename Operator, typename Elem,typename BaseFunctionSpace>
// class ShapeFunctionOperator;



// template<Integer NComponents_>
// class ShapeFunctionOperator<IdentityOperator,Simplex<2,2>, Lagrange1<NComponents_> > : 
// public BaseShapeFunctionOperator<IdentityOperator,3,2,2,1,NComponents_>
// {
// public:
//   virtual const Matrix<Real,3,1>
//   phi_single_component(const Vector<Real,2>& point)
//        {const auto& xi=point[0];
//         const auto& eta=point[1];
//         const Real zeta = 1. - xi - eta;
//         const Matrix<Real,3,1> shapefunction{zeta,xi,eta};
//         return shapefunction;};
// };


// // template<Integer NComponents>
// // class ShapeFunctionOperator<DivergenceOperator,Simplex<2,2>, Lagrange1<NComponents> > : 
// // public BaseShapeFunctionOperator<DivergenceOperator,3,2,2,1,NComponents>
// // {
// // public:
// //   virtual const Vector<DivPhiOperator<2,2,1,NComponents>,3>
// //   phi_single_component(const Vector<Real,2>& qp_point )
// //   {
// //   Vector<DivPhiType<2,2,1,NComponents>,3> tryme; 

// //     tryme[0](0,0)=1.0;
// //     tryme[1](0,0)=1.0;
// //     tryme[2](0,0)=-2.0;

// //   return tryme;};

// // };

// template<Integer NComponents>
// class ShapeFunctionOperator<GradientOperator,Simplex<2,2>, Lagrange1<NComponents> > : 
// public BaseShapeFunctionOperator<GradientOperator,3,2,2,1,NComponents>
// {
// public:

//   virtual const Vector< Matrix<Real,1,2>, 3>
//   phi_single_component(const Vector<Real,2>& point)
//   {
//     Vector< Matrix<Real,1,2>, 3> grad{{-1,-1},{1,0},{0,1}};
//     return grad;
//   }

// };












// template<Integer NComponents_>
// class ShapeFunctionOperator<IdentityOperator, Simplex<2,2>, Lagrange2<NComponents_> > : 
// public BaseShapeFunctionOperator<IdentityOperator,6,2,2,1,NComponents_>
// {
// public:
//       virtual const Matrix<Real,6,1>
//       phi_single_component(const Vector<Real,2>& point) 
//       {
//           const auto& xi=point[0];
//           const auto& eta=point[1];
//           const Real zeta = 1. - xi - eta;
//           Matrix<Real,6,1> shape_function{2.*zeta*(zeta-0.5),
//                                         2.*xi*(xi-0.5),
//                                         2.*eta*(eta-0.5),
//                                         4.*zeta*xi,
//                                         4.*xi*eta, 
//                                         4.*eta*zeta };
//           return shape_function;
//       };
// };

// // template<Integer NComponents>
// // class ShapeFunctionOperator<DivergenceOperator, Simplex<2,2>, Lagrange2<NComponents> > : 
// // public BaseShapeFunctionOperator<DivergenceOperator,6,2,2,1,NComponents>
// // {
// // public:
// //   virtual const Vector<DivPhiType<2,2,1,NComponents>,6>
// //   phi_single_component(const Vector<Real,2>& qp_point )
// //   {
// //     //////////////// E' SBAGLIATO, CI VUOLE DIPENDENZA DALLE VARIABILI ////////////////
// //     Vector<DivPhiType<2,2,1,NComponents>,6> tryme; 
// //     tryme[0](0,0)=(4.0+4.0);
// //     tryme[1](0,0)=(4.0);
// //     tryme[2](0,0)=(4.0);
// //     tryme[3](0,0)=(-8.0);
// //     tryme[4](0,0)=(4.0);
// //     tryme[5](0,0)=(-8.0);
// //     return tryme;
// //  };
// // };


// template<Integer NComponents>
// class ShapeFunctionOperator<GradientOperator, Simplex<2,2>, Lagrange2<NComponents> > : 
// public BaseShapeFunctionOperator<GradientOperator,6,2,2,1,NComponents>
// {
// public:

//   virtual const Vector< Matrix<Real,1,2>, 6>
//   phi_single_component(const Vector<Real,2>& point)
//   {
//    const auto& xi=point[0];
//    const auto& eta=point[1];
//    const Real zeta = 1. - xi - eta;
//    Vector< Matrix<Real,1,2>, 6> grad{
//                                      { - 4 * (1 - xi - eta) + 1, - 4 * (1 - xi - eta) + 1},
//                                      { 4 * xi - 1              , 0                       },
//                                      { 0                       , 4 * eta -1              },
//                                      { 4 - 8 * xi - 4 * eta    , - 4 * xi                },
//                                      { 4 * eta                 , 4 * xi                  },
//                                      { -4 * eta                , 4 - 4 * xi - 8 * eta    } };
//   return grad;
//   }

// };





// template<Integer NComponents_>
// class ShapeFunctionOperator<IdentityOperator, Simplex<2,2>, RT0<NComponents_> > : 
// public BaseShapeFunctionOperator<IdentityOperator,3,2,2,2,NComponents_>
// {
// public:
//   virtual const  Matrix<Real,3,2>
//   phi_single_component(const Vector<Real,2>& point)
//        {const auto& xi=point[0];
//         const auto& eta=point[1];
//         const Matrix<Real,3,2> shapefunction{xi,eta-1,xi-1,eta,xi,eta};
//         return shapefunction;};
// };


// template<Integer NComponents>
// class ShapeFunctionOperator<DivergenceOperator, Simplex<2,2>, RT0<NComponents> > : 
// public BaseShapeFunctionOperator<DivergenceOperator,3,2,2,2,NComponents>
// {
// public:

//   virtual const Vector<Real,3>
//   phi_single_component(const Vector<Real,2>& qp_point )
//   {
//    Vector<Real,3> tryme;
//      tryme[0]=2;
//      tryme[1]=2;
//      tryme[2]=2;    
//    return tryme;
//   };
// };

//  virtual const Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs>
//  grad_phi_single_component(const Vector<Real,2>& point)
//        {
//         //static_assert(Dim==ShapeFunctionDim, "Grad(RT): shape function dim and space dim must be the same")
//         Vector< Matrix<Real,ShapeFunctionDim,Dim>, Ndofs> grad{ {1,0, 0,1},{1,0, 0,1},{1,0, 0,1}} ;
//         return grad;
//         assert("gradient is not defined for RT elements");};
// };


















































































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
       const Vector<Real,Ndofs> &alpha=1.0)
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
       const Vector<Real,Ndofs> &alpha=1.0)
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
       const Vector<Real,Ndofs> &alpha=1.0)
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





















template<typename Elem,typename BaseFunctionSpace>
class MapFromReference2;



template<Integer Dim, Integer ManifoldDim, Integer Order, Integer Continuity, Integer NComponents>
class MapFromReference2<Simplex<Dim,ManifoldDim>,BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents>>
{
 private:
 public:
 using Jacobian=Matrix<Real, Dim, ManifoldDim>;


        inline const Real  operator()(const IdentityOperator& o,const Jacobian& J){return 1.0;}

        inline const Jacobian  operator() (const GradientOperator& o, const Jacobian& J)
         {static_assert(Dim==ManifoldDim,"Dim==ManifoldDim for inverting jacobian");
          auto inv=  inverse(J);
          return inv;}

        inline const Real  operator()(const IdentityOperator& o,const Jacobian& J, Real& map){map=1.0;}

        inline const Jacobian  operator() (const GradientOperator& o, const Jacobian& J, Jacobian& map)
         {static_assert(Dim==ManifoldDim,"Dim==ManifoldDim for inverting jacobian");
          map=  inverse(J);}
};


template<Integer Dim, Integer ManifoldDim, Integer Order, Integer Continuity, Integer NComponents>
class MapFromReference2<Simplex<Dim,ManifoldDim>,BaseFunctionSpace<RaviartThomasFE,Order,Continuity,NComponents>>
{
 private:
 public:
 using Jacobian=Matrix<Real, Dim, ManifoldDim>;
 
        inline const Jacobian operator()(const IdentityOperator&o, const Matrix<Real, Dim, ManifoldDim>& J)
         {Jacobian mapping=J;
          mapping/=det(mapping);
          return mapping;}

         
         inline const Real operator()(const DivergenceOperator&o, const Matrix<Real, Dim, ManifoldDim>& J)
         {return 1.0/det(J);}


        
        inline const Jacobian operator()(const IdentityOperator&o, const Jacobian& J, Jacobian& map)
         {map=J;
          map/=det(map);}

         
         inline const Real operator()(const DivergenceOperator&o, const Jacobian& J,Real& map)
         {map= 1.0/det(J);}

};

















template<typename Elem,typename BaseFunctionSpace>
class MapFromReference3;



template<Integer Dim, Integer ManifoldDim, Integer Order, Integer Continuity, Integer NComponents>
class MapFromReference3<Simplex<Dim,ManifoldDim>,BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents>>
{
 public:
 using Jacobian=Matrix<Real, Dim, ManifoldDim>;


        inline void  init(const IdentityOperator& o,const Jacobian& J){id_=1.0;}

        inline void  init(const GradientOperator& o, const Jacobian& J){grad_=  inverse(J);}

        inline const Real&  operator()(const IdentityOperator& o)const{return id_;}

        inline const Jacobian&  operator() (const GradientOperator& o)const{return grad_;}

 private:
  Real id_;
  Jacobian grad_; 
};


template<Integer Dim, Integer ManifoldDim, Integer Order, Integer Continuity, Integer NComponents>
class MapFromReference3<Simplex<Dim,ManifoldDim>,BaseFunctionSpace<RaviartThomasFE,Order,Continuity,NComponents>>
{

 public:
 using Jacobian=Matrix<Real, Dim, ManifoldDim>;
     
         inline void init(const IdentityOperator&o, const Jacobian& J) {id_=J;id_/=det(id_);}
        
         inline void init(const DivergenceOperator&o, const Jacobian& J){div_= 1.0/det(J);}

         inline const Jacobian&  operator()(const IdentityOperator& o)const {return id_;}

         inline const Real&  operator() (const DivergenceOperator& o)const{return div_;}


 private:
 Jacobian id_;
 Real div_;
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
                  const Vector<Real,Ndofs> &alpha=1.0)
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
                  const Vector<Real,Ndofs> &alpha=1.0)
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































// Base Function: f(x) = x
template<typename T,typename...Parameters>
class ExpressionNew
{
    protected:
    using type=typename T::type;
    using FunctionType=typename T::FunctionType;
    using GradientType=typename T::GradientType;

    struct Implementation
    {
    virtual ~Implementation() {};

    virtual type evaluate(const Parameters &... params)
    {type filler; 
    return filler;};

    virtual GradientType evaluate(const GradientOperator&o, const Parameters &... params)
    {GradientType filler; 
    return filler;};
    };
    
    protected:
    std::shared_ptr<Implementation> self_;

    public:
    ExpressionNew()
    :   self_(std::make_shared<Implementation>())
    {}

    ExpressionNew(const std::shared_ptr<Implementation>& self)
    :   self_(self)
    {}

    type operator () (const Parameters &... params) const 
    { std::cout<<"ExpressionNew standard"<<std::endl;
      return self_->evaluate(params...); }
    GradientType operator () (const GradientOperator&o,const Parameters &... params) const 
    { return self_->evaluate(o,params...); }
    
};


// Unary Function: u(-f(x))
template<typename T,typename...Parameters>
class ExpressionUnaryMinusNew : public ExpressionNew<T,Parameters...>
{
public:

    using type=typename ExpressionNew<T,Parameters...>::type;
    protected:
    struct Implementation : ExpressionNew<T,Parameters...>::Implementation
    {
        ExpressionNew<T,Parameters...> f;

        Implementation(const ExpressionNew<T,Parameters...>& f1):   
        f(f1)
        {};

        virtual type evaluate(const Parameters&... params) override
        { 

          auto provo=f(params...);
          std::cout<<"provo"<<std::endl;
          std::cout<<provo()<<std::endl;
          return -f(params...); 
        }
    };

    public:
    ExpressionUnaryMinusNew(const ExpressionNew<T,Parameters...>& f):   
    ExpressionNew<T,Parameters...>(std::make_shared<Implementation>(f))
    {}
};
template<typename T,typename...Parameters>
inline ExpressionUnaryMinusNew<T,Parameters...> operator - 
(const ExpressionNew<T,Parameters...>& f) { return ExpressionUnaryMinusNew<T,Parameters...>(f); }




template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace,typename QP>
class ExpressionShapeNew: public ExpressionNew<ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>,QP>
{ 
public:
  using T=ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>;
  using ExprT=ExpressionNew<T,QP>;
  using Map=MapFromReference3<Elem,BaseFunctionSpace>;
  using FunctionType=typename T::FunctionType;
  using GradientType=typename T::GradientType;
  using type=FunctionType;
     protected:
    struct Implementation: ExprT::Implementation
    {
        std::shared_ptr<T> t_ptr;
        std::shared_ptr<Map> map_ptr;

        Implementation(const T& t,const Map& map):   
        t_ptr(std::make_shared<T>(t)),
        map_ptr(std::make_shared<Map>(map))
        {};

        virtual FunctionType 
        evaluate(const QP& qp_points) override
        { 
          return (*t_ptr).function(Operator::id())(); 
        };

        virtual GradientType 
        evaluate(const GradientOperator&o, const QP& qp_points) override
        { 
          return (*t_ptr).function(o)(); 
        };


    protected:
        ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace> shape;
    };  
public: 



    ExpressionShapeNew():
    ExprT
    (std::make_shared<Implementation>(ExprT()))    
    {};
    ExpressionShapeNew(const T& t, const Map& map):
    ExprT
    (std::make_shared<Implementation>(t,map))
    {};

};


template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace,typename QP>
class ExpressionShapeGradientNew: public ExpressionNew<ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>,QP>
{ 
public:
  using T=ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>;
  using ExprT=ExpressionNew<T,QP>;
  using Map=MapFromReference3<Elem,BaseFunctionSpace>;
  using FunctionType=typename T::FunctionType;
  using GradientType=typename T::GradientType;
  using type=GradientType;
     protected:
    struct Implementation: ExprT::Implementation
    {
        std::shared_ptr<T> t_ptr;
        std::shared_ptr<Map> map_ptr;

        Implementation(const T& t,const Map& map):   
        t_ptr(std::make_shared<T>(t)),
        map_ptr(std::make_shared<Map>(map))
        {};

        // virtual FunctionType 
        // evaluate(const QP& qp_points) override
        // { 
        //   return (*t_ptr).function(Operator::id())(); 
        // };

        virtual GradientType 
        evaluate(const GradientOperator&o, const QP& qp_points) override
        { 
          return (*t_ptr).function(o)(); 
        };


    protected:
        ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace> shape;
    };  
public: 



    ExpressionShapeGradientNew():
    ExprT
    (std::make_shared<Implementation>(ExprT()))    
    {};
    ExpressionShapeGradientNew(const T& t, const Map& map):
    ExprT
    (std::make_shared<Implementation>(t,map))
    {};

};






// Grad Function: grad(f(x))
template<typename T,typename...Parameters>
class GradExpressionNew : public ExpressionNew<T,Parameters...>
{
public:
    using type=typename T::GradientType;

    protected:
    struct Implementation : ExpressionNew<T,Parameters...>::Implementation
    {
        // using FunctionType=typename T::FunctionType;
        // using GradientType=typename T::GradientType;
        ExpressionNew<T,Parameters...> f;


        Implementation(const ExpressionNew<T,Parameters...>& f1):   
        f(f1)
        {};

        virtual type evaluate(const GradientOperator& o, const Parameters&... params) override
        { 
          return f(o,params...); 
        }
    };

 
    public:
    GradExpressionNew(const ExpressionNew<T,Parameters...>& f)
    :   ExpressionNew<T,Parameters...>(std::make_shared<Implementation>(f))
    {}
   type operator () (const Parameters &... params) const 
    { 
      std::cout<<"son qua"<<std::endl;
      auto ee=(*this).self_->evaluate(Operator::grad(),params...);
      std::cout<<ee()<<std::endl;
      return (*this).self_->evaluate(Operator::grad(),params...); }

};


template<typename T,typename...Parameters>
inline GradExpressionNew<T,Parameters...> GradNew 
(const ExpressionNew<T,Parameters...>& f) { return GradExpressionNew<T,Parameters...>(f); }
















 template<class Derived>
    class OperatorExpression {
    public:


        Derived& derived() { return static_cast<Derived &>(*this); }

        // virtual ~Expression() { }

        ///@return itself as its Derived type
        inline constexpr const Derived &derived() const { 
            return static_cast<const Derived &>(*this); 
        }

        ///@return itself as its Derived type
        inline operator Derived &() {
            return derived();
        }

        ///@return itself as its Derived type
        inline constexpr  operator const Derived &() const {
            return derived();
        }

        ///@return a string with the name of the class
        virtual std::string getClass() const {
            return "Expression";
        }
    };





template<typename...Inputs>
class Eval;

template<typename Left, typename Right>
class OperatorMultiply; //: OperatorExpression<OperatorMultiply<Left,Right>>



template<typename Left, typename Right>
class MultiplyProva;

class FQPprova1;

class QPprova1;

class QPprova2;

template<typename Derived>
class ExpressionProva;

template<typename Derived, typename S,Integer NQPoints, Integer Dim>
class QPprova;

template<typename Derived, typename S,Integer NQPoints, Integer NComponents, Integer Dim>
class FQPprova;

template<typename Operator>
class OperatorType;



template<typename Derived>
class ExpressionProva
{public:
 Derived& derived() { return static_cast<Derived &>(*this); }};



template<typename Derived, typename S,Integer NQPoints, Integer NComponents, Integer Dim>
class FQPprova:
public ExpressionProva<FQPprova<Derived, S, NQPoints, NComponents,Dim>>
// public OperatorExpression<FQPprova<Derived, S, NQPoints, NComponents, Dim>>
{
     public:
     using Point = Vector<Real,Dim> ;
     using type= FQPValues < S, NQPoints,NComponents>;
      Derived& derived() { return static_cast<Derived &>(*this); }
      const Derived& derived()const { return static_cast<const Derived &>(*this); }

      // Vector<S,NComponents> & apply(const Point& point);

      FQPValues<S,NQPoints,NComponents>& apply(const Matrix<Real,NQPoints,Dim> & qp_points);
      // inline FQPValues<S,NQPoints,NComponents> eval(const Matrix<Real,NQPoints,Dim> & qp_points)
      // {
      //  for(Integer qp=0;qp<NQPoints;qp++)
      //     { 
      //       qp_points.get_row(qp,row_);
      //       const auto& row_comp=derived().apply(row_);
      //       for(Integer n_comp=0; n_comp<NComponents; n_comp++)
      //           {
      //             value_[n_comp][qp]=row_comp[n_comp];
      //           }
      //     }
      //     return value_;
      // };

      inline FQPValues<S,NQPoints,NComponents> eval(const Matrix<Real,NQPoints,Dim> & qp_points)
      {
       // for(Integer qp=0;qp<NQPoints;qp++)
       //    { 
       //      qp_points.get_row(qp,row_);
       //      const auto& row_comp=derived().prova(row_);
       //      for(Integer n_comp=0; n_comp<NComponents; n_comp++)
       //          {
       //            value_[n_comp][qp]=row_comp[n_comp];
       //          }
       //    }
          return derived().apply(qp_points);
      };


    protected:
      FQPValues<S,NQPoints,NComponents> value_; 
      Point row_; 
};










template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace,typename QP>
class Shapeprova
: public ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>, 
  public FQPprova<Shapeprova<QuadratureRule,  Elem, BaseFunctionSpace,QP>,
  typename ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::TotFuncType,
  ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::NQPoints,
  ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Ntot,
  ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Dim>
{ 
public:
  using SFO=ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>;
  using Map=MapFromReference3<Elem,BaseFunctionSpace>;
  using S=typename ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>::TotFuncType;
  static constexpr Integer NQPoints=ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::NQPoints;
  static constexpr Integer Ntot=ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Ntot;
  static constexpr Integer  Dim= ShapeFunctionOperator4<QuadratureRule,  Elem, BaseFunctionSpace>::Dim;
  using Point = Vector<Real,Dim> ;

  FQPValues<S,NQPoints,Ntot>& apply(const Matrix<Real,NQPoints,Dim> & qp_points)
     {
      value_=sfo_->function(Operator::id())();
      return value_;
      // vec_mat_[0](0,0)=point[0];vec_mat_[0](0,1)=point[0];
      // // vec_mat_[0](1,0)=point[0];vec_mat_[0](1,1)=point[0];

      // vec_mat_[1](0,0)=point[1];vec_mat_[1](0,1)=point[1];
      // // vec_mat_[1](1,0)=point[1];vec_mat_[1](1,1)=point[1];
      // return vec_mat_;
     }




  Shapeprova(const SFO& sfo, const Map& map):
    sfo_(std::make_shared<SFO>(sfo)),
    map_(std::make_shared<Map>(map))
    {};

//      protected:
//     struct Implementation: ExprT::Implementation
//     {
//         using FunctionType=typename T::FunctionType;
//         using GradientType=typename T::GradientType;
//         std::shared_ptr<T> t_ptr;
//         std::shared_ptr<Map> map_ptr;

//         Implementation(const T& t,const Map& map):   
//         t_ptr(std::make_shared<T>(t)),
//         map_ptr(std::make_shared<Map>(map))
//         {};

//         virtual FunctionType 
//         evaluate(const IdentityOperator& o,const QP& qp_points) override
//         { 
//           // const auto& mapping=(*map_ptr)(o);
//           return (*t_ptr).function(o)(); 
//         };
//         virtual GradientType 
//         evaluate(const GradientOperator& o,const QP& qp_points) override
//         { 
//           // const auto& mapping=(*map_ptr)(o);
//           return (*t_ptr).function(o)(); 
//         };

//     protected:
//         ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace> shape;
//     };  
private:
  std::shared_ptr<SFO> sfo_;
  std::shared_ptr<Map> map_;
  FQPValues<S,NQPoints,Ntot> value_;
};








class FQPprova1: public FQPprova<FQPprova1,Matrix<Real,1,2>,6,2,2>,
                 public ExpressionProva<FQPprova1>
{
 public: 
     
     static constexpr Integer NQPoints=6;
     static constexpr Integer NComponents=2;
     static constexpr Integer Dim=2;

     using Point = Vector<Real,Dim> ;
     using T=Matrix<Real,1,2>;
     using S=Vector<T,NComponents>;

     // S& apply(const Point& point)
     // {
     //  vec_mat_[0](0,0)=point[0];vec_mat_[0](0,1)=point[0];
     //  vec_mat_[1](0,0)=point[1];vec_mat_[1](0,1)=point[1];
     //  return vec_mat_;
     // }

     S& val(const Point& point)
     {
      vec_mat_[0](0,0)=point[0];vec_mat_[0](0,1)=point[0];
      vec_mat_[1](0,0)=point[1];vec_mat_[1](0,1)=point[1];
      return vec_mat_;
     }

     FQPValues<T,NQPoints,NComponents>& apply(const Matrix<Real,NQPoints,Dim> & qp_points)
     {
       for(Integer qp=0;qp<NQPoints;qp++)
          { 
            qp_points.get_row(qp,row_);
            const auto& row_comp=val(row_);
            for(Integer n_comp=0; n_comp<NComponents; n_comp++)
                {
                  fqp_[n_comp][qp]=row_comp[n_comp];
                }
          }

      return fqp_;}

  protected:
    S vec_mat_;
    FQPValues<T,NQPoints,NComponents> fqp_;
};


template<typename Derived, typename S,Integer NQPoints, Integer Dim>
class QPprova: 
public ExpressionProva<QPprova<Derived, S, NQPoints, Dim>>
{
     public:
     using Point = Vector<Real,Dim> ;
     using type= QPValues < S, NQPoints>;
      Derived& derived() { return static_cast<Derived &>(*this); }
      const Derived& derived()const { return static_cast<const Derived &>(*this); }
      S& apply(const Point& point);

      inline QPValues<S,NQPoints> eval(const Matrix<Real,NQPoints,Dim> & qp_points)
      {
       for(Integer qp=0;qp<NQPoints;qp++)
          {
            qp_points.get_row(qp,row_);
            value_[qp]=derived().apply(row_);
          }

        return value_;
      };

    protected:
      QPValues<S,NQPoints> value_; 
      Point row_; 
};

class QPprova1: public QPprova<QPprova1,Matrix<Real,3,1>,6,2>,
                public ExpressionProva<QPprova1>
{
 public: 
     static constexpr Integer Dim=2;
     using Point = Vector<Real,Dim> ;
      Matrix<Real,3,1>& apply(const Point& point)
     {
      mat_(0,0)=point[0];
      mat_(1,0)=point[1];
      mat_(2,0)=point[0]+point[1];
      return mat_;
     }
  protected:
    Matrix<Real,3,1> mat_;
};


class QPprova2: public QPprova<QPprova2,Matrix<Real,1,3>,6,2>,
                public ExpressionProva<QPprova2>
{
 public: 
     static constexpr Integer Dim=2;
     using Point = Vector<Real,Dim> ;
      Matrix<Real,1,3>& apply(const Point& point)
     {
      mat_(0,0)=point[0];
      mat_(0,1)=point[1];
      mat_(0,2)=point[0]+point[1];
      return mat_;
     }
  protected:
    Matrix<Real,1,3> mat_;
};








template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight>
class OperatorType< MultiplyProva< Matrix<T, RowsLeft, Common>, Matrix<T, Common, ColsRight> > >
{ public:
  using type=Matrix<T, RowsLeft, ColsRight>;
};

template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight, Integer NQPoints, Integer NComponents>
class OperatorType< MultiplyProva< QPValues < Matrix<T, RowsLeft, Common>, NQPoints> ,
                                   FQPValues< Matrix<T, Common, ColsRight>, NQPoints, NComponents> > >
{ public:
   using type=FQPValues< Matrix<T, RowsLeft, ColsRight>, NQPoints, NComponents>;
};

template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight, Integer NQPoints>
class OperatorType< MultiplyProva< QPValues< Matrix<T, RowsLeft, Common>, NQPoints> ,
                                   QPValues< Matrix<T, Common, ColsRight>, NQPoints> > >
{ public:
   using type=QPValues< Matrix<T, RowsLeft, ColsRight>, NQPoints>;
};



template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
class OperatorType< MultiplyProva< QPprova < DerivedLeft, TLeft, NQPoints,Dim> ,
                                   FQPprova< DerivedRight,TRight, NQPoints, NComponents,Dim> > >
{ public:


   using type= FQPValues<typename OperatorType<MultiplyProva<TLeft,TRight>>::type , NQPoints, NComponents>;
};

template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
class OperatorType< MultiplyProva< QPprova< DerivedLeft, TLeft,  NQPoints,Dim> ,
                                   QPprova< DerivedRight,TRight, NQPoints,Dim> > >
{ public:


   using type= QPValues<typename OperatorType<MultiplyProva<TLeft,TRight>>::type , NQPoints>;
};
// template<typename Derived1, typename Derived2,
//         typename T,Integer Rows1,Integer Cols1, Integer Rows2,Integer Cols2,
//         Integer NQPoints, Integer NComponents, Integer Dim>
// class MultiplyProva<QPprova<Derived1,Matrix<T,Rows1,Cols1>,NQPoints,Dim>,
//                     FQPprova<Derived2,Matrix<T,Rows2,Cols2>,NQPoints,NComponents,Dim>>
                    // : 
                    // public FQPprova<Derived2,Matrix<Real,Rows1,Cols2>,NQPoints,NComponents,Dim>

                    // : 
                    // public FQPprova<MultiplyProva<QPprova<Derived1,Matrix<T,Rows1,Cols1>,NQPoints,Dim>,
                    //             FQPprova<Derived2,Matrix<T,Rows2,Cols2>,NQPoints,NComponents,Dim>>,
                    //             Matrix<Real,Rows1,Cols2>,NQPoints,NComponents,Dim>
// template<typename Left, typename Right>
// class MultiplyProva: public ExpressionProva<MultiplyProva<Left,Right>>

// {
//   public:
//     MultiplyProva(const Left& left, const Right&right): 
//     left_(left),
//     right_(right)
//     {};

//     template<typename QP>
//     typename OperatorType<MultiplyProva<Left,Right>>::type eval(const QP & qp_points)
//     {return left_.eval(qp_points)*right_.eval(qp_points);};

//   private:
//   Left left_;
//   Right right_;
// };


// template<typename Left, typename Right>
// MultiplyProva<Left,Right>
// operator*(const Left& left, const Right&right)// (ExpressionProva<Left> const& u, ExpressionProva<Right> const& v) 
// {
// return MultiplyProva<Left,Right>(left,right);//static_cast<const Left>(u), static_cast<const Right>(v));
// }


template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
class MultiplyProva< QPprova <DerivedLeft, TLeft, NQPoints,Dim>, 
                     QPprova <DerivedRight,TRight,NQPoints, Dim> >

: public QPprova<MultiplyProva<QPprova <DerivedLeft, TLeft, NQPoints,Dim>,
                               QPprova <DerivedRight,TRight,NQPoints,Dim> >,
                typename OperatorType<MultiplyProva<typename TLeft::type,typename TRight::type>>::type,
                NQPoints, Dim >
{
  public:
    using Left= QPprova<DerivedLeft, TLeft, NQPoints,Dim>;
    using Right=QPprova<DerivedRight,TRight,NQPoints,Dim>;
    MultiplyProva(const Left& left, const Right&right): 
    left_(left.derived()),
    right_(right.derived())
    {};

    template<typename QP>
    typename OperatorType< MultiplyProva<typename Left::type, typename Right::type>>::type  eval(const QP & qp_points)
    {return left_.eval(qp_points)*right_.eval(qp_points);};

  private:
  DerivedLeft left_;
  DerivedRight right_;
};


template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
class MultiplyProva< QPprova<DerivedLeft, TLeft, NQPoints,Dim>, 
                     QPprova<DerivedRight,TRight,NQPoints,Dim> >
operator*(const QPprova<DerivedLeft, TLeft, NQPoints,Dim>& left, 
          const QPprova<DerivedRight,TRight,NQPoints, Dim>&right)
{
 return MultiplyProva<QPprova<DerivedLeft, TLeft, NQPoints,Dim>, 
                      QPprova<DerivedRight,TRight,NQPoints,Dim> >(left,right);//static_cast<const Left>(u), static_cast<const Right>(v));
}



template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
class MultiplyProva< QPprova <DerivedLeft, TLeft, NQPoints,Dim>, 
                     FQPprova<DerivedRight,TRight,NQPoints,NComponents, Dim> >

: public FQPprova<MultiplyProva<QPprova <DerivedLeft, TLeft, NQPoints,Dim>,
                              FQPprova<DerivedRight,TRight,NQPoints,NComponents, Dim> >,
                typename OperatorType<MultiplyProva<typename TLeft::type,typename TRight::type>>::type,
                NQPoints, NComponents, Dim >
{
  public:
    using Left=QPprova <DerivedLeft, TLeft, NQPoints,Dim>;
    using Right=FQPprova<DerivedRight,TRight,NQPoints,NComponents, Dim>;
    MultiplyProva(const Left& left, const Right&right): 
    left_(left.derived()),
    right_(right.derived())
    {};

    template<typename QP>
    typename OperatorType< MultiplyProva<typename Left::type, typename Right::type>>::type  eval(const QP & qp_points)
    {return left_.eval(qp_points)*right_.eval(qp_points);};

  private:
  DerivedLeft left_;
  DerivedRight right_;
};

template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
class MultiplyProva< QPprova <DerivedLeft, TLeft, NQPoints,Dim>, 
                     FQPprova<DerivedRight,TRight,NQPoints,NComponents, Dim> >
operator*(const QPprova <DerivedLeft, TLeft, NQPoints,Dim>& left, 
          const FQPprova<DerivedRight,TRight,NQPoints,NComponents, Dim>&right)
{
return MultiplyProva<QPprova <DerivedLeft, TLeft, NQPoints,Dim>, FQPprova<DerivedRight,TRight,NQPoints,NComponents, Dim> >(left,right);//static_cast<const Left>(u), static_cast<const Right>(v));
}

// template<typename Derived1, typename Derived2,
//         typename T,Integer Rows1,Integer Cols1, Integer Rows2,Integer Cols2,
//         Integer NQPoints, Integer NComponents, Integer Dim>
// inline MultiplyProva<QPprova<Derived1,Matrix<T,Rows1,Cols1>,NQPoints,Dim>,
//                      FQPprova<Derived2,Matrix<T,Rows2,Cols2>,NQPoints,NComponents,Dim> > operator *
// (const QPprova<Derived1,Matrix<T,Rows1,Cols1>,NQPoints,Dim>& left, 
//  const FQPprova<Derived2,Matrix<T,Rows2,Cols2>,NQPoints,NComponents,Dim> & right) 
// { return MultiplyProva<QPprova<Derived1,Matrix<T,Rows1,Cols1>,NQPoints,Dim>,
//                        FQPprova<Derived2,Matrix<T,Rows2,Cols2>,NQPoints,NComponents,Dim> >
//                        (left,right);};























template<typename S,Integer NQPoints, Integer NComponents, Integer Dim>
class FQPExpression: public OperatorExpression<FQPExpression<S, NQPoints, NComponents, Dim>>
{
     public:
     using Point = Vector<Real,Dim> ;

      virtual ~FQPExpression(){};
      FQPExpression(){};
      virtual Vector<S,NComponents> & apply(const Point& point)=0;

      inline FQPValues<S,NQPoints,NComponents>& eval(const Matrix<Real,NQPoints,Dim> & qp_points)
      {
       for(Integer qp=0;qp<NQPoints;qp++)
          { 
            qp_points.get_row(qp,row_);
            const auto& row_comp=apply(row_);
            for(Integer n_comp=0; n_comp<NComponents; n_comp++)
                {
                  value_[n_comp][qp]=row_comp[n_comp];
                }
          }
          return value_;
      };

    protected:
      FQPValues<S,NQPoints,NComponents> value_; 
      Point row_; 
};

template<typename S,Integer NQPoints, Integer Dim>
class QPExpression: public OperatorExpression<QPExpression<S, NQPoints, Dim>>
{
     public:
     using Point = Vector<Real,Dim> ;

      virtual ~QPExpression(){};
      QPExpression(){};
      virtual S& apply(const Point& point)=0;

      inline QPValues<S,NQPoints>& eval(const Matrix<Real,NQPoints,Dim> & qp_points)
      {
       for(Integer qp=0;qp<NQPoints;qp++)
          {
            qp_points.get_row(qp,row_);
            value_[qp]=apply(row_);
          }

        return value_;
      };

    protected:
      QPValues<S,NQPoints> value_; 
      Point row_; 
};


template<typename T,Integer Rows1,Integer Cols1, Integer Rows2,Integer Cols2,Integer NQPoints, Integer NComponents, Integer Dim>
class OperatorMultiply<QPExpression<Matrix<T,Rows1,Cols1>,NQPoints,Dim>,
                       FQPExpression<Matrix<T,Rows2,Cols2>,NQPoints,NComponents,Dim>>: 
                       public FQPExpression<Matrix<T,Rows1,Cols2>,NQPoints,NComponents,Dim>
{
  public:
    using Result=Matrix<T,Rows1,Cols2>;
    template<typename Left, typename Right>
    OperatorMultiply(const Left& left, const Right&right): 
    left_(std::make_shared<Left>(left)),
    right_(std::make_shared<Right>(right))
    {};
 
    virtual Vector<Matrix<T,Rows1,Cols2>,NComponents> & apply(const Vector<T,Dim>& point)
    { 
      return vec_mat_; };

    FQPValues<Result,NQPoints,NComponents> eval(const Matrix<Real,NQPoints,Dim> & qp_points)
    {

      auto e1=left_->eval(qp_points);
      std::cout<<"ci arrivo left"<<e1<<std::endl;
      auto e2=right_->eval(qp_points);
      std::cout<<"ci arrivo right"<<e2<<std::endl;

     return left_->eval(qp_points)*right_->eval(qp_points);
    };

  private:
  std::shared_ptr<QPExpression<Matrix<T,Rows1,Cols1>,NQPoints,Dim>> left_;
  std::shared_ptr<FQPExpression<Matrix<T,Rows2,Cols2>,NQPoints,NComponents,Dim>> right_;
  Vector<Matrix<T,Rows1,Cols2>,NComponents> vec_mat_;
};



class QP1: public QPExpression<Matrix<Real,3,1>,6,2>
{
 public: 
     static constexpr Integer Dim=2;
     using Point = Vector<Real,Dim> ;
     virtual Matrix<Real,3,1>& apply(const Point& point) override
     {
      mat_(0,0)=point[0];
      mat_(1,0)=point[1];
      mat_(2,0)=point[0]+point[1];
      return mat_;
     }
  protected:
    Matrix<Real,3,1> mat_;
};



class QP2: public QPExpression<Matrix<Real,3,3>,6,2>
{
 public: 
     static constexpr Integer Dim=2;
     using Point = Vector<Real,Dim> ;
     virtual Matrix<Real,3,3>& apply(const Point& point) override
     {
      mat_(0,0)=point[0];mat_(0,1)=point[0];                  mat_(0,2)=point[0];
      mat_(1,0)=point[1];mat_(1,1)=point[1];                  mat_(1,2)=point[0];
      mat_(2,0)=point[0]+point[1];mat_(2,1)=point[1]+point[0];mat_(2,2)=point[0];
      return mat_;
     }
  protected:
    Matrix<Real,3,3> mat_;
};


class FQP1: public FQPExpression<Matrix<Real,1,2>,6,2,2>
{
 public: 
     static constexpr Integer Dim=2;
     using Point = Vector<Real,Dim> ;
     using S=Vector<Matrix<Real,1,2>,2>;
     virtual S& apply(const Point& point) override
     {
      vec_mat_[0](0,0)=point[0];vec_mat_[0](0,1)=point[0];
      // vec_mat_[0](1,0)=point[0];vec_mat_[0](1,1)=point[0];

      vec_mat_[1](0,0)=point[1];vec_mat_[1](0,1)=point[1];
      // vec_mat_[1](1,0)=point[1];vec_mat_[1](1,1)=point[1];
      return vec_mat_;
     }
  protected:
    S vec_mat_;
};
























// Base Function: f(x) = x
template<typename T,typename...Parameters>
class ExpressionT
{
    protected:
    using FunctionType=typename T::FunctionType;
    using GradientType=typename T::GradientType;
    struct Implementation
    {
    virtual ~Implementation() {};

    virtual T evaluate(const Parameters &... params)
    {T filler; 
    return filler;};

    virtual FunctionType evaluate(const IdentityOperator& o,const Parameters &... params)
    {FunctionType filler; 
    return filler;};
    virtual GradientType evaluate(const GradientOperator& o,const Parameters &... params)
    {GradientType filler; 
     return filler;};
    // virtual typename T::DivergenceType evaluate(const DivergenceOperator& o,const Parameters &... params){};
    };

    public:
    ExpressionT()
    :   self_(std::make_shared<Implementation>())
    {}

    ExpressionT(const std::shared_ptr<Implementation>& self)
    :   self_(self)
    {}

    T operator () (const Parameters &... params) const 
    { return self_->evaluate(params...); }

    FunctionType operator () (const IdentityOperator& o, const Parameters &... params) const 
    { return self_->evaluate(o,params...); }

    GradientType operator () (const GradientOperator& o, const Parameters &... params) const 
    { return self_->evaluate(o,params...); }
    // typename T::GradientType operator () (const GradientOperator& o, const Parameters &... params) const { return self_->evaluate(o,params...); }
    // typename T::DivergenceType operator () (const DivergenceOperator& o, const Parameters &... params) const { return self_->evaluate(o,params...); }

    private:
    std::shared_ptr<Implementation> self_;
};

// Unary Function: u(-f(x))
template<typename T,typename...Parameters>
class ExpressionUnaryMinusT : public ExpressionT<T,Parameters...>
{
public:


    protected:
    struct Implementation : ExpressionT<T,Parameters...>::Implementation
    {
        using FunctionType=typename T::FunctionType;
        ExpressionT<T,Parameters...> f;


        Implementation(const ExpressionT<T,Parameters...>& f1):   
        f(f1)
        {};

        virtual FunctionType evaluate(const IdentityOperator& o,const Parameters&... params) override
        { 
          return -f(o,params...); 
        }

    };

    public:
    ExpressionUnaryMinusT(const ExpressionT<T,Parameters...>& f)
    :   ExpressionT<T,Parameters...>(std::make_shared<Implementation>(f))
    {}
};
template<typename T,typename...Parameters>
inline ExpressionUnaryMinusT<T,Parameters...> operator - 
(const ExpressionT<T,Parameters...>& f) { return ExpressionUnaryMinusT<T,Parameters...>(f); }






// Grad Function: grad(f(x))
template<typename T,typename...Parameters>
class GradExpressionT : public ExpressionT<T,Parameters...>
{
public:
    using FunctionType=typename T::GradientType;

    protected:
    struct Implementation : ExpressionT<T,Parameters...>::Implementation
    {
        // using FunctionType=typename T::FunctionType;
        // using GradientType=typename T::GradientType;
        ExpressionT<T,Parameters...> f;


        Implementation(const ExpressionT<T,Parameters...>& f1):   
        f(f1)
        {};

        virtual FunctionType evaluate(const GradientOperator& o, const Parameters&... params) override
        { const auto& boh=f(o,params...);
          std::cout<<"grad f(o,params...)"<<boh()<<std::endl;
          return f(o,params...); 
        }
    };

    public:
    GradExpressionT(const ExpressionT<T,Parameters...>& f)
    :   ExpressionT<T,Parameters...>(std::make_shared<Implementation>(f))
    {}
};


template<typename T,typename...Parameters>
inline GradExpressionT<T,Parameters...> Grad 
(const ExpressionT<T,Parameters...>& f) { return GradExpressionT<T,Parameters...>(f); }




template<typename S,Integer NQPoints,typename QP>
class ExpressionMatrixFunction: public ExpressionT<QPValues<S,NQPoints>,QP>
{ 
  using T=QPValues<S,NQPoints>;
  using FunctionType=typename T::type;
  using ExprT=ExpressionT<T,QP>;
  
     protected:
    struct Implementation: ExpressionT<T,QP>::Implementation
    {
        Implementation()   
        {};


        template<typename Point>
        S& value(const Point& point)
        {
          s_tmp_(0,0)=point[0]; s_tmp_(0,1)=point[1];
          s_tmp_(1,0)=point[1]; s_tmp_(1,1)=point[0]+point[1];
          return s_tmp_; 
        }

        virtual T 
        evaluate(const QP& qp_points) override
        { 

          for(Integer qp=0;qp<NQPoints;qp++)
            {
              qp_points.get_row(qp,row);
              t_tmp_[qp]=value(row);

            }
          return t_tmp_;
        };
    protected:
      S s_tmp_;
      T t_tmp_;
      Vector<Real,QP::Cols> row;
    };  
public: 

    ExpressionMatrixFunction():
    ExprT(std::make_shared<Implementation>())    
    {};



};


template<Integer NQPoints,typename QP>
class ExpressionRealFunction: public ExpressionT<QPValues<Real,NQPoints>,QP>
{ 
  using S=Real;
  using T=QPValues<S,NQPoints>;
  using FunctionType=typename T::type;
  using ExprT=ExpressionT<T,QP>;
  
     protected:
    struct Implementation: ExpressionT<T,QP>::Implementation
    {
        Implementation()   
        {};


        template<typename Point>
        S& value(const Point& point)
        {
          s_tmp_=3.1415*point[0]+1.0*point[1];
          return s_tmp_; 
        }

        virtual T 
        evaluate(const QP& qp_points) override
        { 

          for(Integer qp=0;qp<NQPoints;qp++)
            {
              qp_points.get_row(qp,row);
              t_tmp_[qp]=value(row);

            }
          return t_tmp_;
        };
    protected:
      S s_tmp_;
      T t_tmp_;
      Vector<Real,QP::Cols> row;
    };  
public: 

    ExpressionRealFunction():
    ExprT(std::make_shared<Implementation>())    
    {};



};



// Binary Function: u(f(x) * g(x))
template<typename T,typename...Parameters>
class ExpressionBinaryMultiplyT :  public ExpressionT<T,Parameters...>
{
    protected:
    struct Implementation : ExpressionT<T,Parameters...>::Implementation
    {
        ExpressionT<T,Parameters...> f;
        ExpressionT<T,Parameters...> g;
        Implementation(const ExpressionT<T,Parameters...>& f1,const ExpressionT<T,Parameters...>& g1)
        :   f(f1), g(g1)
        {};

        virtual T evaluate(const Parameters &...params) override { return f(params...) * g(params...); }
    };

    public:
    ExpressionBinaryMultiplyT<T,Parameters...>(const ExpressionT<T,Parameters...>& f,const ExpressionT<T,Parameters...>& g)
    :   ExpressionT<T,Parameters...>(std::make_shared<Implementation>(f, g))
    {}
};









template<typename QuadratureRule, typename Elem,typename BaseFunctionSpace,typename QP>
class ExpressionShape: public ExpressionT<ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>,QP>
{ 
  using T=ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>;
  using ExprT=ExpressionT<T,QP>;
  using Map=MapFromReference3<Elem,BaseFunctionSpace>;
     protected:
    struct Implementation: ExprT::Implementation
    {
        using FunctionType=typename T::FunctionType;
        using GradientType=typename T::GradientType;
        std::shared_ptr<T> t_ptr;
        std::shared_ptr<Map> map_ptr;

        Implementation(const T& t,const Map& map):   
        t_ptr(std::make_shared<T>(t)),
        map_ptr(std::make_shared<Map>(map))
        {};

        virtual FunctionType 
        evaluate(const IdentityOperator& o,const QP& qp_points) override
        { 
          // const auto& mapping=(*map_ptr)(o);
          return (*t_ptr).function(o)(); 
        };
        virtual GradientType 
        evaluate(const GradientOperator& o,const QP& qp_points) override
        { 
          // const auto& mapping=(*map_ptr)(o);
          return (*t_ptr).function(o)(); 
        };

    protected:
        ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace> shape;
    };  
public: 

    ExpressionShape():
    ExprT
    (std::make_shared<Implementation>(ExprT()))    
    {};
    ExpressionShape(const T& t, const Map& map):
    ExprT
    (std::make_shared<Implementation>(t,map))
    {};

};





template<typename QuadratureRule, typename Elem, typename BaseFunctionSpace,typename QP>
class ExpressionBinaryMultiplyT<ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>,QP>: 
public ExpressionT<ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>,QP>
{

  using S=ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>;
  using ExprS=ExpressionT<S,QP>;
  using FunctionType=typename S::FunctionType;

    protected:
    template<typename S,Integer NQPoints>
    struct Implementation: ExprS::Implementation
    {   
        // using ExprT=ExpressionT<T,QP>;
        using ExprT=ExpressionT<QPValues<S,NQPoints>,QP>;
        ExprT expr_t;
        ExprS expr_s;
        Implementation(const ExprT& t,const ExprS& s)
        :   expr_t(t),expr_s(s)
        {};

        virtual FunctionType evaluate(const IdentityOperator& o,const QP &qp_points) override 
        {return expr_t(qp_points)*expr_s(o,qp_points); }
    };

    

    public:
    template<typename S,Integer NQPoints>
    ExpressionBinaryMultiplyT(const ExpressionT<QPValues<S,NQPoints>,QP>& t,const ExprS& s): 
    ExprS(std::make_shared<Implementation<S,NQPoints>>(t,s)) {};
};



template<typename QuadratureRule, typename Elem, typename BaseFunctionSpace,typename QP, typename S, Integer NQPoints>
inline ExpressionBinaryMultiplyT<ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>,QP> operator *
(const ExpressionT<QPValues<S,NQPoints>,QP>& t, 
 const ExpressionT<ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>,QP> & s) 
{ return ExpressionBinaryMultiplyT<ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>,QP>(t,s);};



// template<Integer NQPoints, Integer Ndofs, Integer Dim, Integer ManifoldDim, Integer ShapeFunctionDim,Integer NComponents>
// class BaseShapeFunctionOperator3
// { 
//   public:
//   /// maybe use Matrix<Real,ShapefunctionDim1,ShapeFunctionDim2>
//   using FuncType=Vector<Real,ShapeFunctionDim>;
//   using GradType=Matrix<Real,ShapeFunctionDim*NComponents,Dim>;
//   using DivType=Real;
//   using Point=Vector<Real,Dim>;
//   using QP=Matrix<Real,NQPoints,Dim>;
//   static constexpr Integer Ntot=Ndofs*NComponents;




//   public:



//   virtual ~BaseShapeFunctionOperator3(){};
  
//   virtual void
//   value(const IdentityOperator& o,  const Point& point,    Vector<FuncType,Ndofs>& func_ )   =0;

//   virtual void
//   value(const GradientOperator&o,   const Point& point,    Vector<GradType,Ndofs>& func_grad_)=0;

//   virtual void
//   value(const DivergenceOperator&o, const Point& qp_point, Vector<DivType,Ndofs>& func_div_)  =0;



//   FQPValues<FuncType,NQPoints,Ndofs>& reference(const IdentityOperator&o){return reference_func_values_;}
//   FQPValues<GradType,NQPoints,Ntot>& reference(const GradientOperator&o){return reference_grad_values_;}

//   FQPValues<FuncType,NQPoints,Ntot>& operator()(const IdentityOperator&o){return func_values_;}

//   void operator()(const IdentityOperator&o,const QP & qp_points)//,FQPValues<FuncType,NQPoints,Ndofs>& reference_func_values_)
//   {
//    for(Integer qp=0;qp<NQPoints;qp++)
//     {
//     qp_points.get_row(qp,qp_point_);
//     value(o,qp_point_,func_);
//     for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//      {
//       for(Integer ii=0;ii<ShapeFunctionDim;ii++)
//         reference_func_values_[n_dof][qp][ii]=func_[n_dof][ii];     
//      }
//     }
//   };



//   template<Integer NQPoints, typename Mapping>
//   void operator()(FQPValues<FuncType,NQPoints,Ntot>& func_values_,
//                   const FQPValues<FuncType,NQPoints,Ndofs>& reference_func_values_,
//                   const IdentityOperator&o,
//                   const Mapping& mapping, 
//                   const Vector<Real,Ndofs> &alpha=1.0)
//   {
//   Integer n_tot;
//   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//      {
//       for(Integer n_comp=0;n_comp<NComponents;n_comp++)
//       {
//         for(Integer qp=0;qp<NQPoints;qp++)
//         { 
//           n_tot=n_dof * NComponents + n_comp;
//           func_values_[n_tot][qp].zero();
//           func_tmp_=alpha[n_dof] * mapping*reference_func_values_[n_dof][qp];
//           func_values_[n_tot][qp].row(n_comp,func_tmp_);
//         }        
//       }
//      }
//   };





//   template<Integer NQPoints>
//   void operator()(const GradientOperator& o, const Matrix<Real,NQPoints,Dim>& qp_points)
//   {
//    Vector<Real,Dim> qp_point;

//    for(Integer qp=0;qp<NQPoints;qp++)
//     {
//     qp_points.get_row(qp,qp_point_);
//     value(o,qp_point_,grad_);
//     for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
//      {
//       reference_grad_values_[n_dof][qp]=grad_[n_dof];           
//      }
//     }
//   };


  // template<Integer NQPoints, typename Mapping>
  // Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs>
  // operator()(const GradientOperator& o,
  //      const Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs>& reference_grad_phi,
  //      const Mapping& mapping, 
  //      const Vector<Real,Ndofs> &alpha=1.0)
  // {
  //  Vector<Vector<Vector<Matrix<Real,ShapeFunctionDim*NComponents,Dim>,NQPoints>,NComponents>,Ndofs> result;
  //  Vector<Real,Dim> row;
  //   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
  //    {
  //     for(Integer n_comp=0;n_comp<NComponents;n_comp++)
  //     {
  //       for(Integer qp=0;qp<NQPoints;qp++)
  //       {
  //         result[n_dof][n_comp][qp].zero();
  //         reference_grad_phi[n_dof][n_comp][qp].get_row(n_comp,row); 
  //         // std::cout<<" firs row ="<<std::endl;
  //         // reference_grad_phi[n_dof][n_comp][qp].describe(std::cout);     
  //         // row.describe(std::cout);     
  //         // std::cout<<" mapping ="<<std::endl;
  //         // mapping.describe(std::cout);    
  //         row=alpha[n_dof]*mapping*row;      
  //         result[n_dof][n_comp][qp].row(n_comp,row);
  //         // std::cout<<" result ="<<std::endl;
  //         // result[n_dof][n_comp][qp].describe(std::cout);        
  //       }
  //     }
  //    }
  //  return result;
  // };

  // template<Integer NQPoints>
  // Vector<Vector< Real, NQPoints>,Ndofs>
  // operator()(const DivergenceOperator&o,const Matrix<Real,NQPoints,Dim>& qp_points)
  // {
  //       Vector<Vector< Real, NQPoints>,Ndofs> div;
  //       Vector<Real,Dim> qp_point;

  //       for(Integer qp=0;qp<NQPoints;qp++)
  //       {
  //         qp_points.get_row(qp,qp_point);
  //         auto div_phi_single=phi_single_component(o,qp_point);
  //             for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
  //                  div[n_dof][qp]=div_phi_single[n_dof];
  //       }
  //       return div;};

  // template<Integer NQPoints,typename Mapping>
  // Vector<Vector<Vector< Matrix<Real,NComponents,1>, NQPoints>,NComponents>,Ndofs>
  // operator()(const DivergenceOperator&o,
  //      const Vector<Vector< Real, NQPoints>,Ndofs>& divphi_reference,
  //      const Mapping& mapping,
  //      const Vector<Real,Ndofs> &alpha=1.0)
  // {
  //       Vector<Vector<Vector< Matrix<Real,NComponents,1>, NQPoints>,NComponents>,Ndofs> divphi;
  //       for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
  //        for(Integer n_comp=0;n_comp<NComponents;n_comp++)
  //         for(Integer qp=0;qp<NQPoints;qp++)
  //           { divphi[n_dof][n_comp][qp].zero();
  //             divphi[n_dof][n_comp][qp](n_comp,0)=alpha[n_dof] * mapping * divphi_reference[n_dof][qp];
  //           } 
  //   return divphi;};


//   private: 
//   Vector<FuncType,Ndofs> func_;
//   Vector<GradType,Ndofs> grad_;
//   Vector<DivType,Ndofs>  div_;
//   Point qp_point_;

//   FQPValues<FuncType,NQPoints,Ndofs> reference_func_values_;
//   FQPValues<GradType,NQPoints,Ndofs> reference_grad_values_;

//   FQPValues<FuncType,NQPoints,Ntot> func_values_;
//   FQPValues<GradType,NQPoints,Ntot> grad_values_;
//   FuncType func_tmp_;
// };











// template<Integer NQPoints, typename Elem,typename BaseFunctionSpace>
// class ShapeFunctionOperator3;



// template<Integer NQPoints, Integer NComponents>
// class ShapeFunctionOperator3<NQPoints,Simplex<2,2>, Lagrange1<NComponents> > : 
// public BaseShapeFunctionOperator3<NQPoints,3,2,2,1,NComponents>
// {
// public:
//   static constexpr Integer Ndofs=3;
//   static constexpr Integer Dim=2;
//   static constexpr Integer ManifoldDim=2;
//   static constexpr Integer ShapeFunctionDim=1;
//   using FuncType=Vector<Real,ShapeFunctionDim>;
//   using GradType=Matrix<Real,ShapeFunctionDim,Dim>;
//   using DivType=Real;
//   using Point=Vector<Real,Dim>;
//   using QP=Matrix<Real,NQPoints,Dim>;



//  virtual void value(const IdentityOperator& o, const Point& point, Vector<FuncType,Ndofs>& func )
//        {func[0][0]=point[0];
//         func[1][0]=point[1];
//         func[2][0]=1. - point[0] - point[1];};
//   ;

//   virtual void value(const GradientOperator&o, const Point& point, Vector<GradType,Ndofs>& func_grad)
//    {
//     func_grad[0](0,0)=-1;  func_grad[0](0,1)=-1; 
//     func_grad[1](0,0)=+1;  func_grad[1](0,1)= 0; 
//     func_grad[2](0,0)= 0;  func_grad[2](0,1)=+1; 
//   } 

//   virtual void value(const DivergenceOperator&o, const Point& qp_point, Vector<DivType,Ndofs>& func_div_)
//   {
//     std::cout<<" Lagrange1 divergence not implemented"<<std::endl;
//   }


// };


































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
         //const auto& vec=BilinearFormQPvalues<NQPoints,Multiply,Multiply>(A,element_trial[n_dof_trial][n_comp_trial],B,element_test[n_dof_test][n_comp_test]);
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
  virtual std::shared_ptr<Mother> operator[](const Integer i)=0;
  virtual void print()=0;
};

class Son1: public Mother
{
public:
   Son1(const Integer& i): value(i){};
  virtual void print(){std::cout<<value<<std::endl;}
  virtual std::shared_ptr<Mother> operator[](const Integer i)
  {assert(i==0 && "kids have only 1 component");
    return NULL;};
private:
  Integer  value;
};

class Son2: public Mother
{
 public:
  Son2():children() {};

  virtual void print(){for(Integer ii=0;ii< children.size();ii++)
                        children[ii]->print();}

  virtual std::shared_ptr<Mother> operator[](const Integer i){return children[i];};

     void add(const std::shared_ptr<Mother>& kid)
    {
        children.push_back(kid);
    } 

private:
  std::vector <std::shared_ptr< Mother > > children;
};





template <class Self>
struct Base
{
 Self& self() { return static_cast<Self&>(*this); }

  int algorithm(int x)
  {
    self().prologue();
    if (x > 42)
      x = self().downsize(x);
    x = self().crunch(x);
    self().epilogue();
    return x;
  }

  void prologue()
  {}

  int downsize(int x)
  { return x % 42; }

  int crunch(int x)
  { return -x; }

  void epilogue()
  {}


  Real var(const Real& p){var_=p; return var_;};
protected:
  Real var_;


};

struct Derived : Base<Derived>
{
  int algorithm(int x)
  {
 std::cout << "We're almost done!\n";

    return x;
  }

  int downsize(int x)
  {
    while (x > 42) x /= 2;
    return x;
  }

  void epilogue()
  { std::cout << "We're done!\n"; }

  Real var(const Real& p){var_=p; return var_+1;};
};





template <Integer QPOrder, typename TrialFunction,typename TestFunction, typename OperatorTrial, typename OperatorTest>//,typename MeshT >
void BilinearFormIntegrator(const TrialFunction& trialspace,const TestFunction& testspace)
{



std::shared_ptr<Son1> son1=std::make_shared<Son1>(1);
std::shared_ptr<Son1> son2=std::make_shared<Son1>(2);
std::shared_ptr<Son1> son3=std::make_shared<Son1>(3);

std::shared_ptr<Son1> daughter1=std::make_shared<Son1>(4);
std::shared_ptr<Son1> daughter2=std::make_shared<Son1>(5);


std::shared_ptr<Son2> mom=std::make_shared<Son2>();
std::shared_ptr<Son2> uncle1=std::make_shared<Son2>();
std::shared_ptr<Son2> uncle2=std::make_shared<Son2>();
std::shared_ptr<Son2> grandma=std::make_shared<Son2>();

son1->print();
son2->print();
son3->print();
std::cout<<"mom"<<std::endl;
mom->add(son1);
mom->add(son2);
std::cout<<"uncle1"<<std::endl;

uncle1->add(daughter1);
uncle1->add(daughter2);
std::cout<<"grandma"<<std::endl;

grandma->add(mom);
grandma->add(uncle1);

auto mom1=grandma->operator[](0);
auto kid1=mom1->operator[](0);
kid1->print();
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
const auto& qp_points=gauss.qp_points();
const auto& qp_weights=gauss.qp_weights();


 using multiplication_tuple_trial= VariadicTupleType<Contraction,Multiply,Multiply,Contraction>;
 using multiplication_tuple_test= VariadicTupleType<Contraction,Multiply,Multiply,Contraction>;


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


Vector<Real, 3 > alpha_trial=1;



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
const auto& n_elements=mesh->n_elements();





std::vector<Vector<Real,Dim>> points(Npoints);
Matrix<Real, Dim, ManifoldDim> J;
Real volume;
Elem elem;


for(Integer ii=0;ii<Npoints;ii++)
   elem.nodes[ii]=ii;





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
points[0][0]=0;
points[0][1]=0;
points[1][0]=2;
points[1][1]=3;
points[2][0]=1;
points[2][1]=4;

for(Integer nn=0;nn<points.size();nn++)
{
  points[nn].describe(std::cout);
}



  jacobian(elem,points,J);



std::cout<<"BEGIN MATRIX LOOP"<<std::endl;

MatrixLoop<2,2>(tupleshape_trial,tupleshape_test,
                reference_shape_trial,reference_shape_test,
                tuplemap_trial,tuplemap_test,
                operators_trial,operators_test,
                J,qp_points);

std::cout<<"END MATRIX LOOP"<<std::endl;



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

  for(Integer elem_iter=0 ; elem_iter < 1;elem_iter++)//n_elements ; elem_iter++)
  {
    // auto elemnodes_global=mesh->elem(elem_iter).nodes;
    // for(Integer mm=0;mm<points.size();mm++)
    //    points[mm]=mesh->point(elemnodes_global[mm]);
    // jacobian(elem,points,J);
    std::cout<<IndexTo<0,0,2>::value<<std::endl;
    std::cout<<IndexTo<0,1,2>::value<<std::endl;
    std::cout<<IndexTo<0,2,2>::value<<std::endl;
    std::cout<<IndexTo<0,3,2>::value<<std::endl;

    std::cout<<IndexTo<1,0,2>::value<<std::endl;
    std::cout<<IndexTo<1,1,2>::value<<std::endl;
    std::cout<<IndexTo<1,2,2>::value<<std::endl;
    std::cout<<IndexTo<1,3,2>::value<<std::endl;
    volume=unsigned_volume(elem,points);
    const auto& mapping_trial0=map_trial0(grad,J);
    const auto& mapping_trial1=map_trial1(identity,J);

const auto reference_shape_trial0=std::get<0>(reference_shape_trial);

Vector<Vector<Vector<double, 1>, 6>, 3> reference_shape_trial05=std::get<2>(reference_shape_trial);
const auto& mapping_trial05=map_trial0(identity,J);

const auto element_trial0=trial0(grad,reference_shape_trial0,mapping_trial0);
const auto element_trial05=trial0(identity,reference_shape_trial05,mapping_trial05);

const auto element_trial1=trial1(identity,std::get<1>(reference_shape_trial),mapping_trial1);
// const auto element_trial2=trial2(identity,std::get<2>(reference_shape_trial),mapping_trial0);
// const auto element_trial3=trial3(identity,std::get<3>(reference_shape_trial),mapping_trial0);
auto naltro=trial0.phi_single_component(std::get<0>(operators_trial),points[0]);
// std::cout<<"naltro--"<<nm<<std::endl;
// std::cout<<naltro<<std::endl;
// std::cout<<"mapping_trial05--"<<nm<<std::endl;
// std::cout<<mapping_trial05<<std::endl;
// std::cout<<"reference_shape_trial05--"<<nm<<std::endl;
// std::cout<<reference_shape_trial05<<std::endl;
// std::cout<<"jacobian="<<std::endl;
// J.describe(std::cout);
// std::cout<<"detJ="<<std::endl;
// std::cout<<det(J)<<std::endl;
// std::cout<<"element_trial0="<<std::endl;
// element_trial0.describe(std::cout);
// std::cout<<"element_trial1="<<std::endl;
// element_trial1.describe(std::cout);

//     const auto& mapping_trial=map_trial(J);
//     const auto& mapping_test=map_test(J);
//     const auto& element_trial=trial.phiN(reference_trial,mapping_trial,alpha_trial);
//     const auto& element_test=test.phiN(reference_test,mapping_test,alpha_test);

//     for(Integer n_dof_trial=0;n_dof_trial<Ndofs_trial;n_dof_trial++)
//       for(Integer n_comp_trial=0;n_comp_trial<NComponents_trial;n_comp_trial++)
//        for(Integer n_dof_test=0;n_dof_test<Ndofs_test;n_dof_test++)
//         for(Integer n_comp_test=0;n_comp_test<NComponents_test;n_comp_test++)
//           {const auto& i=n_comp_trial * NComponents_trial+ n_dof_trial;
//            const auto& j=n_comp_test * NComponents_test + n_dof_test;
//            // mat(i,j)=dot(vec,qp_weights)*volume;
//           }

  }

// ShapeFunctionOperator3<NQPoints, Elem, Lagrange1<1>> s;
// s(identity,qp_points);
// s(grad,qp_points);
// auto ref_identity=s.reference(identity);

// auto ref_grad=s.reference(grad);

// std::cout<<ref_identity()<<std::endl;
// std::cout<<ref_grad()<<std::endl;

ShapeFunctionOperator4<GaussPoints<Elem,QPOrder>, Elem, Lagrange1<4>> s4;
s4(qp_points);
// s4(qp_points,identity);
// s4(qp_points,grad);
// s4(qp_points,identity,identity,grad);
std::cout<<s4.reference(identity)()<<std::endl;
const auto& mapping_identity_lagrange1=map_trial0(identity,J);
const auto& mapping_grad_lagrange1=map_trial0(grad,J);
// s4(identity,mapping_identity_lagrange1);
// s4(grad,mapping_grad_lagrange1);
s4(identity,mapping_identity_lagrange1,grad,mapping_grad_lagrange1);
std::cout<<s4.function(identity)()<<std::endl;
std::cout<<s4.function(grad)()<<std::endl;


Matrix<Real,5,5> mat1=0;
Matrix<Real,2,2> mat2{0,1,2,3};

assign(mat1,mat2,2,3);
std::cout<<mat1<<std::endl;
// s4(grad,qp_points);

using TTT=Matrix<Real,2,1>;
using SSS=Matrix<Real,2,2>;
TTT ttt=1;
SSS sss=3;
Vector<TTT,6> initvec=ttt;
Vector<Vector<TTT,6>,3> initvecvec=initvec;
FQPValues<TTT,6,3> fqpval(initvecvec);
QPValues<SSS,6> qpval(sss);

std::cout<<fqpval()<<std::endl;
std::cout<<qpval()<<std::endl;

auto result=qpval*fqpval+fqpval;
std::cout<<result()<<std::endl;

// FQPValues<TTT,3,6> result2=qpval*fqpval+fqpval;
// std::cout<<fqpval()<<std::endl;
// std::cout<<qpval()<<std::endl;
// std::cout<<result()<<std::endl;
// std::cout<<result2()<<std::endl;


grandma->print();


    MapFromReference3<Elem,Lagrange1<1>> map3;
    using QP=typename GaussPoints<Elem,QPOrder>::qp_points_type;
    using QuadratureRule=GaussPoints<Elem,QPOrder>;
    using BaseFunctionSpace=Lagrange1<1>;
    using ShapeSpace=ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>;
    ShapeSpace s;
    map3.init(Operator::id(),J);
    map3.init(Operator::grad(),J);
    const auto& o=Operator::id();
    const auto& o_grad=Operator::grad();
    // const auto& mapping=map3(o);
    s(qp_points);
    s(o,map3(o));
    s(o_grad,map3(o_grad));
    // ExpressionShape<GaussPoints<Elem,QPOrder>, Elem, BaseFunctionSpace,QP> a(s,map3);
    // ExpressionMatrixFunction<Matrix<Real,2,2>,NQPoints,QP> z0;

    // ExpressionRealFunction<NQPoints,QP> g0;

    // auto grada=-Grad(a);
    // auto grada0=grada(Operator::grad(),qp_points);
    // FQPValues<Matrix<Real, 2,1>, 6, 6> pr1;
    // QPValues<Real, 6> pr2;
    
    // QPValues<Matrix<Real, 4,3>, 6> pr3;
    // FQPValues<Matrix<Real, 3,2>, 6, 6> pr4;
    // auto fff=pr2*pr1;
    // auto es=pr3*pr4;
    // std::cout<<"pr3()"<<std::endl;
    // std::cout<<pr3()[0]<<std::endl;
    // std::cout<<"pr4()"<<std::endl;
    // std::cout<<pr4()[0][0]<<std::endl;
    // std::cout<<"es"<<std::endl;
    // std::cout<<es()<<std::endl;
    // Matrix<Real, 3,1> m3=3;
    // Matrix<Real, 1,2> m4=4;
    // auto m5=m3*m4;
    // std::cout<<"m5="<<m5<<std::endl;
    // auto m=g0*(g0*a);//ExpressionBinaryMultiplyT<ShapeFunctionOperator4< QuadratureRule,  Elem, BaseFunctionSpace>,QP>(z0,a);
    // auto r=z0*a;
    // auto r0=r(Operator::id(),qp_points);
    // auto m0=m(Operator::id(),qp_points);

    // auto g1=g0(qp_points);

    // auto a1=a(Operator::id(),qp_points);
    // auto ga=g1*a1;
    // auto ga2=g0(qp_points)*a(Operator::id(),qp_points);
    // auto r0=r(Operator::id(),qp_points);
    // auto m0=m(Operator::id(),qp_points);
    // auto a0=a(Operator::id(),qp_points);
    // auto b=-a;
    // auto b0=b(Operator::id(),qp_points);
    // auto e=z0(qp_points);

    // //Operator::identity;
    //  std::cout<<"so acca"<<std::endl;
    //  std::cout<<a0()<<std::endl;
    //  std::cout<<" z0=========="<<std::endl;
    //  std::cout<<e()<<std::endl;
     // std::cout<<" r begin =========="<<std::endl;
     // std::cout<<r0()<<std::endl;
     // std::cout<<" m begin =========="<<std::endl;
     // std::cout<<m0()<<std::endl;
     // std::cout<<" grad(a) begin =========="<<std::endl;
     // std::cout<<grada0()<<std::endl;


     // ExpressionShapeNew<GaussPoints<Elem,QPOrder>, Elem, Lagrange1<2>,QP> a_new(s,map3);
     // ExpressionShapeGradientNew<GaussPoints<Elem,QPOrder>, Elem, Lagrange1<2>,QP> a_gradnew(s,map3);
     // // auto b_new=-a_new;
     // // auto b_new0=b_new(qp_points);
     // // auto grada_new=GradNew(a_new);
     // // auto grada0_new=grada_new(qp_points);
     // auto gradb_new=-GradNew(a_new);
     // auto gradb0_new=gradb_new(qp_points);
     // auto gradc_new=a_gradnew;
     // auto gradc0_new=gradc_new(qp_points);
     // std::cout<<" grada0_new begin =========="<<std::endl;
     // std::cout<<grada0_new()<<std::endl;
     //  std::cout<<" gradb0_new begin =========="<<std::endl;
     // std::cout<<gradb0_new()<<std::endl; 
     //  std::cout<<" gradc_new begin =========="<<std::endl;
     // std::cout<<gradc0_new()<<std::endl;    
      // std::cout<<" b_new0() begin =========="<<std::endl;
     // std::cout<<b_new0()<<std::endl;    



using QPexpr= QPExpression<Matrix<Real,3,1>,6,2>;
using FQPexpr= FQPExpression<Matrix<Real,1,2>,6,2,2>;
using QPexpr2= QPExpression<Matrix<Real,3,3>,6,2>;
using FQPexpr2= FQPExpression<Matrix<Real,3,2>,6,2,2>;

QP1 qp1;
// qp1(qp_points);
// auto qp1val=qp1();
// std::cout<<"qp1val"<<qp1val()<<std::endl;

QP2 qp2;
FQP1 fqp1;


auto prodqp=OperatorMultiply<QPexpr,FQPexpr>(qp1,fqp1);
auto prodqp2= OperatorMultiply<QPexpr2,FQPexpr2>(qp2,prodqp);
auto culo=prodqp.eval(qp_points);
auto qulo= qp2.eval(qp_points)* prodqp.eval(qp_points);
std::cout<<"culo"<<culo()<<std::endl;
auto culo2=prodqp2.eval(qp_points);
std::cout<<"qulo2"<<qulo()<<std::endl;





  FQPprova1 fqpp1;
  QPprova1  qpp1;
  QPprova2  qpp2;
  Shapeprova<GaussPoints<Elem,QPOrder>, Elem, Lagrange1<1>,QP> sfprova(s,map3);

  std::cout<<"fqpp1"<<std::endl;
  std::cout<<fqpp1.eval(qp_points)<<std::endl;
  std::cout<<"qpp1"<<std::endl;
  std::cout<<qpp1.eval(qp_points)<<std::endl;
  auto fqpp2=qpp1*fqpp1;
  std::cout<<"fqpp2"<<std::endl;
  std::cout<<fqpp2.eval(qp_points)<<std::endl;
  
    std::cout<<"qpp2"<<std::endl;
  std::cout<<qpp2.eval(qp_points)<<std::endl;

 Derived d;
  std::cout << d.algorithm(420)<<std::endl;
    std::cout << d.var(420)<<std::endl;
  
  auto fqpp3=qpp2*qpp1*fqpp1;
  std::cout<<"fqpp3"<<std::endl;
  std::cout<<fqpp3.eval(qp_points)<<std::endl;

  std::cout<<"sfprova"<<std::endl;
  std::cout<<sfprova.eval(qp_points)<<std::endl;

  auto sfprova2= qpp1*sfprova;
  std::cout<<"sfprova2"<<std::endl;
  std::cout<<sfprova2.eval(qp_points)<<std::endl;

};



















}



#endif






