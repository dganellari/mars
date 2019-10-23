#ifndef MARS_OPERATOR_APPLY_HPP
#define MARS_OPERATOR_APPLY_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"

namespace mars {

template<typename T>
class OperatorApply;

template<typename Left,typename Right>
class Addition;

template<typename Left,typename Right>
class Subtraction;

template<typename T,Integer Rows,Integer Cols>
class OperatorApply<Addition< Matrix<T,Rows,Cols>,Matrix<T,Rows,Cols> > >
{       
 public: 
 using Left=Matrix<T,Rows,Cols>;
 using Right=Matrix<T,Rows,Cols>;
 // using Output=OperatorType<Addition<Left,Right> >;
 template<typename Output>
 inline static void apply(      Output& A,
                          const Left& B,
                          const Right& C)
  {

    for(Integer ii=0;ii<Left::Rows;ii++)
    for(Integer jj=0;jj<Left::Cols;jj++)
           A(ii,jj)=B(ii,jj)+C(ii,jj);
        // AdditionMatrixAndTransposedAux<Left,Right>::apply(A,B,C);
  };
};

// template<typename T,Integer Rows,Integer Cols>
// class OperatorApply<Subtraction< Matrix<T,Rows,Cols>,Matrix<T,Rows,Cols> > >
// {       
//  public: 
//  using Left=Matrix<T,Rows,Cols>;
//  using Right=Matrix<T,Rows,Cols>;
//  template<typename Output>
//  inline static void apply(      Output& A,
//                           const Left& B,
//                           const Right& C)
//   {

//     for(Integer ii=0;ii<Left::Rows;ii++)
//     for(Integer jj=0;jj<Left::Cols;jj++)
//            A(ii,jj)=B(ii,jj)-C(ii,jj);
//   };
// };




template<template<class>class Unary, typename U,Integer NQPoints>
class OperatorApply<Unary<QPValues<U,NQPoints>>>
{       
 public:               
 template<typename V>
 inline static void apply(      QPValues<V,NQPoints>& A,
                          const QPValues<U,NQPoints>& B)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
        {
          Unary<U>::apply(A[qp],B[qp]);
        }
  };
};

template<template<class>class Unary, typename U,Integer NQPoints,Integer Ndofs>
class OperatorApply<Unary<FQPValues<U,NQPoints,Ndofs>>>
{       
 public:               
 template<typename V>
 inline static void apply(      FQPValues<V,NQPoints,Ndofs>& A,
                          const FQPValues<U,NQPoints,Ndofs>& B)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
    {
       Unary<U>::apply(A(n_dof,qp),B(n_dof,qp));
    }
  };

};










template<template<class,class>class Binary, typename U,typename V, Integer NQPoints>
class OperatorApply<Binary<QPValues<U,NQPoints>,QPValues<V,NQPoints>>> 
{       
 public:               
 inline static void apply(      QPValues<OperatorType<Binary<U,V>>,NQPoints>& A,
                          const QPValues<U,NQPoints>& B,
                          const QPValues<V,NQPoints>& C)
  {
   for(Integer qp=0;qp<NQPoints;qp++)
    {
      Binary<U,V>::apply(A[qp],B[qp],C[qp]);
    }    
    
  };

};

template<template<class,class>class Binary, typename U,typename V, Integer NQPoints,Integer Ndofs>
class OperatorApply<Binary<FQPValues<U,NQPoints,Ndofs>,FQPValues<V,NQPoints,Ndofs>>> 
{       
 public:               
 inline static void apply(      FQPValues<OperatorType<Binary<U,V>>,NQPoints,Ndofs>& A,
                          const FQPValues<U,NQPoints,Ndofs>& B,
                          const FQPValues<V,NQPoints,Ndofs>& C)
  {
  for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
   for(Integer qp=0;qp<NQPoints;qp++)
    {
       Binary<U,V>::apply(A(n_dof,qp),B(n_dof,qp),C(n_dof,qp));
    }
  };

};





template<template<class,class>class Binary, typename S, typename T,Integer NQPoints,Integer Ndofs>
class OperatorApply<Binary<QPValues<S,NQPoints>,FQPValues<T,NQPoints,Ndofs>>>
                      
{    
 public:            
 inline static void apply(      FQPValues<OperatorType<Binary<S,T>>,NQPoints,Ndofs>& A,
                          const QPValues<S,NQPoints>& B,
                          const FQPValues<T,NQPoints,Ndofs>& C)
  {
   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
    for(Integer qp=0;qp<NQPoints;qp++)
       {
        Binary<S,T>::apply(A(n_dof,qp),B[qp],C(n_dof,qp));
       }
  };
};

template<template<class,class>class Binary,typename S, typename T,Integer NQPoints,Integer Ndofs>
class OperatorApply<Binary<FQPValues<S,NQPoints,Ndofs>,QPValues<T,NQPoints>>> 
{    
 public:          
 inline static void apply(      FQPValues<OperatorType<Binary<S,T>>,NQPoints,Ndofs>& A,
                          const FQPValues<S,NQPoints,Ndofs>& B,
                          const QPValues<T,NQPoints>& C)
  {
   for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
    for(Integer qp=0;qp<NQPoints;qp++)
       {
        Binary<S,T>::apply(A(n_dof,qp),B(n_dof,qp),C[qp]);
       }
  };
};

// template<>
// class Assignment<Real> 
// {
//  public:       
//  inline static constexpr void apply(Real& A,const Real& B)
//   {
//            A=B;
//   };
// };



// template<typename T,Integer Rows,Integer Cols>
// class Assignment<Matrix<T,Rows,Cols>> 
// {
//  public:       
 
//  inline static constexpr void apply(Matrix<T,Rows,Cols>& A,const Matrix<T,Rows,Cols>& B)
//   {
//    for(Integer ii=0;ii<Rows;ii++)
//     for(Integer jj=0;jj<Cols;jj++)
//       Assignment<T>::apply(A(ii,jj),B(ii,jj));
//   };

//  inline static constexpr void apply(Transposed<Matrix<T,Rows,Cols>>& A,const Matrix<T,Rows,Cols>& B)
//   {
//   // update pointer to the non-transposed matrix of the transposed matrix A with B itself
//     // std::cout<<"Assignment<Transposed<Matrix<T,Rows,Cols>>> "<<B<<std::endl;

//        A(B);
//     // std::cout<<A<<std::endl;
//   };



// };


// template<typename T,Integer Rows,Integer Cols>
// class Assignment<Transposed<Matrix<T,Cols,Rows>>> 
// {
//  public:       
 
//   inline static constexpr void apply(Matrix<T,Rows,Cols>& A,const Transposed<Matrix<T,Cols,Rows>>& B)
//   {
//    for(Integer ii=0;ii<Rows;ii++)
//     for(Integer jj=0;jj<Cols;jj++)
//       Assignment<T>::apply(A(ii,jj),B(jj,ii));
//   };

//  inline static constexpr void apply(Transposed<Matrix<T,Cols,Rows>>& A,const Transposed<Matrix<T,Cols,Rows>>& B)
//   {
//   // update pointer to the non-transposed matrix of the transposed matrix A 
//   // with the not transposed matrix of B
//    A(B());
//    };


// };


// template<typename S,Integer NQPoints>
// class Assignment<QPValues<S,NQPoints>> 
// {
//  public:       
//  template<typename T>
//  inline static constexpr void apply(QPValues<T,NQPoints>& A,const QPValues<S,NQPoints>& B)
//   {
//     OperatorApply<Assignment<QPValues<T,NQPoints>>>::apply(A,B);

//     // T ok1(1);
//     // S ok2(2);
//    // for(Integer ii=0;ii<NQPoints;ii++)
//    //    {

//    //      Assign(A[ii],B[ii]);

//         // Assignment<S>::apply(A[ii],B[ii]);

//         // Assignment<T>::apply(A[ii],B[ii]);
//         std::cout<<"Assignment<QPValues<T,NQPoints>> "<<A<<" "<<B<<std::endl;
//       // }
//   };







// };

// template<typename T,Integer NQPoints,Integer Ndofs>
// class Assignment<FQPValues<T,NQPoints,Ndofs>> 
// {
//  public:       
//  template<typename S>
//  inline static constexpr void apply(FQPValues<S,NQPoints,Ndofs>& A,const FQPValues<T,NQPoints,Ndofs>& B)
//   {
//     OperatorApply<Assignment<FQPValues<T,NQPoints,Ndofs>>>::apply(A,B);
//    // for(Integer ii=0;ii<Ndofs;ii++)
//    //  for(Integer jj=0;jj<NQPoints;jj++)
//    //    Assign(A[ii][jj],B[ii][jj]);
//   };
// };


}
#endif