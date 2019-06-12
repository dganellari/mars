#ifndef MARS_OPERATORS_HPP
#define MARS_OPERATORS_HPP
#include "mars_base.hpp"

namespace mars {

///////////////////////////////////////////////////////////////
/////////////////////    MAIN OPERATORS    ////////////////////
///////////////////////////////////////////////////////////////
template<typename T>
class UnaryPlus;

template<typename T>
class UnaryMinus;

template<typename Left, typename Right>
class Multiply;

template<typename Left, typename Right>
class Divide;

template<typename Left, typename Right>
class Subtraction;

template<typename Left, typename Right>
class Addition;

template<typename Derived, typename S,Integer NQPoints, Integer Dim>
class QPprova;

template<typename Derived, typename S,Integer NQPoints, Integer NComponents, Integer Dim>
class FQPprova;

template<typename T,Integer NQPoints,Integer NComponents>
class FQPValues;

template<typename T,Integer NQPoints>
class QPValues;




template<typename Operator>
class OperatorType;


// type(T) = type(+T)
template<typename T>
class OperatorType< UnaryPlus< T > >
{ public:
  using type=T;
};


// type(T) = type(-T)
template<typename T>
class OperatorType< UnaryMinus< T > >
{ public:
  using type=T;
};


// type(A)=type(B) = type(A+B)
template<typename Left, typename Right>
class OperatorType< Addition< Left,Right > >
{ public:
  static_assert(std::is_same<Left,Right>::value, " In Addition, Left and Right types must be equal");
  using type=Left;
};


// type(A)=type(B) = type(A-B)
template<typename Left, typename Right>
class OperatorType< Subtraction< Left,Right > >
{ public:
  static_assert(std::is_same<Left,Right>::value, " In Subtraction, Left and Right types must be equal");
  using type=Left;
};

// type(T) = type(T/REAL)
template<typename T>
class OperatorType< Divide<T, Real > >
{ public:
  using type=T;
};

// type(VECTOR) = type(VECTOR * REAL)
template<typename T,Integer Dim>
class OperatorType< Multiply< Vector<T, Dim>, Real > >
{ public:
  using type=Vector<T, Dim>;
};


// type(MATRIX) = type(MATRIX * REAL)
template<typename T, Integer RowsLeft,Integer ColsLeft>
class OperatorType< Multiply< Matrix<T, RowsLeft, ColsLeft>, Real > >
{ public:
  using type=Matrix<T, RowsLeft, ColsLeft>;
};


// type(VECTOR) = type(MATRIX * VECTOR)
template<typename T, Integer RowsLeft,Integer Common>
class OperatorType< Multiply< Matrix<T, RowsLeft, Common>, Vector<T, Common> > >
{ public:
  using type=Vector<T, RowsLeft>;
};

// type(MATRIX) = type(MATRIX * MATRIX)

template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight>
class OperatorType< Multiply< Matrix<T, RowsLeft, Common>, Matrix<T, Common, ColsRight> > >
{ public:
  using type=Matrix<T, RowsLeft, ColsRight>;
};

template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight, Integer NQPoints, Integer NComponents>
class OperatorType< Multiply< QPValues < Matrix<T, RowsLeft, Common>, NQPoints> ,
                                   FQPValues< Matrix<T, Common, ColsRight>, NQPoints, NComponents> > >
{ public:
   using type=FQPValues< Matrix<T, RowsLeft, ColsRight>, NQPoints, NComponents>;
};

template<typename T, Integer RowsLeft,Integer Common,Integer ColsRight, Integer NQPoints>
class OperatorType< Multiply< QPValues< Matrix<T, RowsLeft, Common>, NQPoints> ,
                                   QPValues< Matrix<T, Common, ColsRight>, NQPoints> > >
{ public:
   using type=QPValues< Matrix<T, RowsLeft, ColsRight>, NQPoints>;
};


template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer NComponents,Integer Dim>
class OperatorType< Multiply< QPprova < DerivedLeft, TLeft, NQPoints,Dim> ,
                                   FQPprova< DerivedRight,TRight, NQPoints, NComponents,Dim> > >
{ public:


   using type= FQPValues<typename OperatorType<Multiply<TLeft,TRight>>::type , NQPoints, NComponents>;
};

template< typename DerivedLeft,typename DerivedRight, typename  TLeft,typename TRight, Integer NQPoints,Integer Dim>
class OperatorType< Multiply< QPprova< DerivedLeft, TLeft,  NQPoints,Dim> ,
                                   QPprova< DerivedRight,TRight, NQPoints,Dim> > >
{ public:


   using type= QPValues<typename OperatorType<Multiply<TLeft,TRight>>::type , NQPoints>;
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
}

#endif //MARS_OPERATORS_HPP
