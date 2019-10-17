#ifndef MARS_OPERATOR_TYPE_HPP
#define MARS_OPERATOR_TYPE_HPP
#include "mars_base.hpp"

namespace mars {



class IdentityOperator;

class GradientOperator;

class DivergenceOperator;

class CurlOperator;

template<typename...Ts>
class CompositeOperator;

template<typename ConstType,typename...Inputs>
class ConstantTensor;

template<typename FullSpace,Integer N,typename Operator_,typename FuncType>
class Function;


template<typename MixedSpace, Integer N, typename OperatorType>
class Test;

template<typename MixedSpace, Integer N, typename OperatorType>
class Trial;


template<typename Left,typename Right,Integer QR>
class L2DotProductIntegral;

template<typename Derived>
class Expression;

template<typename T>
class Transposed;

template<typename...Parameters>
class UnaryPlus;

template<typename...Parameters>
class UnaryMinus;

template<typename...Parameters>
class Multiplication;

template<typename...Parameters>
class Contraction2;

template<typename...Parameters>
class Division;

template<typename...Parameters>
class Subtraction;

template<typename...Parameters>
class Addition;



template<typename T,Integer NQPoints>
class QPValues;

template<typename T,Integer NQPoints,Integer Ndofs>
class FQPValues;


template<typename...Operators>
class OperatorTypeHelper;


template<typename...Inputs>
using OperatorType=typename OperatorTypeHelper<Inputs...>::type;




///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////// UnaryMinus for Matrix, QPValues, FQPValues and their transposes                             /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, Integer Rows, Integer Cols,typename...Ts>
class OperatorTypeHelper<UnaryMinus< Matrix<T,Cols,Rows> >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols,typename...Ts>
class OperatorTypeHelper<UnaryMinus< Transposed< Matrix<T,Cols,Rows> > >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols,Integer NQPoints, typename...Ts>
class OperatorTypeHelper<UnaryMinus< QPValues<Transposed< Matrix<T,Cols,Rows>>,NQPoints>  >, Ts...>
{ public:
  using type=QPValues< Matrix<T,Rows,Cols>,NQPoints>;
};

template<typename T, Integer Rows, Integer Cols,Integer NQPoints,Integer Ndofs, typename...Ts>
class OperatorTypeHelper<UnaryMinus< FQPValues<Transposed< Matrix<T,Cols,Rows>>,NQPoints,Ndofs> >, Ts...>
{ public:
  using type=FQPValues< Matrix<T,Rows,Cols>,NQPoints,Ndofs >;
};




///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////// Addition for Matrix, QPValues, FQPValues and their transposes combinations                  /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                             A + B                                             /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Addition<  Matrix<T,Rows,Cols>, Matrix<T,Rows,Cols> >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Addition< QPValues< Matrix<T,Rows,Cols>,NQPoints>, QPValues<Matrix<T,Rows,Cols>,NQPoints> >, Ts...>
{ public:
  using type=QPValues<Matrix<T,Rows,Cols>,NQPoints>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper< Addition< FQPValues< Matrix<T,Rows,Cols>,NQPoints,Ndofs>, FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs> >, Ts...>
{ public:
  using type=FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs>;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                             A + Transpose(B)                                  /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Addition<  Matrix<T,Rows,Cols>, Transposed< Matrix<T,Cols,Rows> > >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Addition< QPValues< Matrix<T,Rows,Cols>,NQPoints>, QPValues<Transposed< Matrix<T,Cols,Rows> >,NQPoints> >, Ts...>
{ public:
  using type=QPValues<Matrix<T,Rows,Cols>,NQPoints>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper< Addition< FQPValues< Matrix<T,Rows,Cols>,NQPoints,Ndofs>, FQPValues<Transposed< Matrix<T,Cols,Rows> >,NQPoints,Ndofs> >, Ts...>
{ public:
  using type=FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs>;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                       Transpose(A) + B                                        /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Addition< Transposed< Matrix<T,Cols,Rows> >, Matrix<T,Rows,Cols> >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Addition< QPValues<Transposed<Matrix<T,Cols,Rows>>,NQPoints>, QPValues<Matrix<T,Rows,Cols>,NQPoints> >, Ts...>
{ public:
  using type=QPValues<Matrix<T,Rows,Cols>,NQPoints>;
};

template<typename T, Integer Rows, Integer Cols,Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper< Addition< FQPValues<Transposed<Matrix<T,Cols,Rows>>,NQPoints,Ndofs> , FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs> >, Ts...>
{ public:
  using type=FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs>;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                Transpose(A) + Transpose(B)                                    /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Addition< Transposed< Matrix<T,Cols,Rows> >, Transposed< Matrix<T,Cols,Rows> > >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Addition< QPValues<Transposed< Matrix<T,Cols,Rows> >,NQPoints>, QPValues< Transposed< Matrix<T,Cols,Rows> >,NQPoints> >, Ts...>
{ public:
  using type=QPValues<Matrix<T,Rows,Cols>,NQPoints>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper< Addition< FQPValues<Transposed< Matrix<T,Cols,Rows> >,NQPoints,Ndofs> , FQPValues<Transposed< Matrix<T,Cols,Rows> >,NQPoints,Ndofs> >, Ts...>
{ public:
  using type=FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs>;
};











///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////// Subtraction for Matrix, QPValues, FQPValues and their transposes combinations               /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                             A - B                                             /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Subtraction<  Matrix<T,Rows,Cols>, Matrix<T,Rows,Cols> >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Subtraction< QPValues< Matrix<T,Rows,Cols>,NQPoints>, QPValues<Matrix<T,Rows,Cols>,NQPoints> >, Ts...>
{ public:
  using type=QPValues<Matrix<T,Rows,Cols>,NQPoints>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper< Subtraction< FQPValues< Matrix<T,Rows,Cols>,NQPoints,Ndofs>, FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs> >, Ts...>
{ public:
  using type=FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs>;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                             A - Transpose(B)                                  /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Subtraction<  Matrix<T,Rows,Cols>, Transposed< Matrix<T,Cols,Rows> > >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Subtraction< QPValues< Matrix<T,Rows,Cols>,NQPoints>, QPValues<Transposed< Matrix<T,Cols,Rows> >,NQPoints> >, Ts...>
{ public:
  using type=QPValues<Matrix<T,Rows,Cols>,NQPoints>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper< Subtraction< FQPValues< Matrix<T,Rows,Cols>,NQPoints,Ndofs>, FQPValues<Transposed< Matrix<T,Cols,Rows> >,NQPoints,Ndofs> >, Ts...>
{ public:
  using type=FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs>;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                       Transpose(A) - B                                        /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Subtraction< Transposed< Matrix<T,Cols,Rows> >, Matrix<T,Rows,Cols> >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Subtraction< QPValues<Transposed<Matrix<T,Cols,Rows>>,NQPoints>, QPValues<Matrix<T,Rows,Cols>,NQPoints> >, Ts...>
{ public:
  using type=QPValues<Matrix<T,Rows,Cols>,NQPoints>;
};

template<typename T, Integer Rows, Integer Cols,Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper< Subtraction< FQPValues<Transposed<Matrix<T,Cols,Rows>>,NQPoints,Ndofs> , FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs> >, Ts...>
{ public:
  using type=FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs>;
};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                Transpose(A) - Transpose(B)                                    /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Subtraction< Transposed< Matrix<T,Cols,Rows> >, Transposed< Matrix<T,Cols,Rows> > >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Subtraction< QPValues<Transposed< Matrix<T,Cols,Rows> >,NQPoints>, QPValues< Transposed< Matrix<T,Cols,Rows> >,NQPoints> >, Ts...>
{ public:
  using type=QPValues<Matrix<T,Rows,Cols>,NQPoints>;
};

template<typename T, Integer Rows, Integer Cols, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper< Subtraction< FQPValues<Transposed< Matrix<T,Cols,Rows> >,NQPoints,Ndofs> , FQPValues<Transposed< Matrix<T,Cols,Rows> >,NQPoints,Ndofs> >, Ts...>
{ public:
  using type=FQPValues<Matrix<T,Rows,Cols>,NQPoints,Ndofs>;
};







///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////// Multiplication between Matrices ( (M1 X N1) x (M2 x N2)) or a matrix and a scalar ( 1 x 1)  /////////
///////// and also their transposes combinations                                                      /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T,Integer Rows1,Integer Cols1,Integer Rows2, Integer Cols2, typename...Ts>
class OperatorTypeHelper<Multiplication<Matrix<T,Rows1,Cols1>, 
                                        Matrix<T,Rows2,Cols2> >, Ts...>
{
  public:
  using type=
    // if (1 x 1) x (1 x 1)
    typename
    std::conditional< (Rows1==1 && Cols1==1 && Rows2==1 && Cols2==1),
                       Matrix<T,1,1>,
                       // if (1 x 1) x (M x N)
                       typename
                       std::conditional< (Rows1==1 && Cols1==1),
                                          Matrix<T,Rows2,Cols2>,
                                          // if (M x N) x (1 x 1)
                                          typename
                                          std::conditional< (Rows2==1 && Cols2==1),
                                                             Matrix<T,Rows1,Cols1>,
                                                             Matrix<T,Rows1,Cols2>
                                                          >::type
                                        >::type
                    >::type;
};

template<typename T,Integer Rows1,Integer Cols1,Integer Rows2, Integer Cols2, typename...Ts>
class OperatorTypeHelper<Multiplication<Matrix<T,Rows1,Cols1>, 
                         Transposed<Matrix<T,Cols2,Rows2>> >, Ts...>
{
  public:
  using type=typename OperatorTypeHelper<Multiplication<Matrix<T,Rows1,Cols1>,
                                                        Matrix<T,Rows2,Cols2>>,Ts...>::type;
};



template<typename T,Integer Rows1,Integer Cols1,Integer Rows2, Integer Cols2, typename...Ts>
class OperatorTypeHelper<Multiplication<Transposed<Matrix<T,Cols1,Rows1>>, 
                                        Matrix<T,Rows2,Cols2> >, Ts...>
{
  public:
  using type=typename OperatorTypeHelper<Multiplication<Matrix<T,Rows1,Cols1>,
                                                        Matrix<T,Rows2,Cols2>>,Ts...>::type;
};

template<typename T,Integer Rows1,Integer Cols1,Integer Rows2, Integer Cols2, typename...Ts>
class OperatorTypeHelper<Multiplication<Transposed<Matrix<T,Cols1,Rows1>>, 
                                        Transposed<Matrix<T,Cols2,Rows2>> >, Ts...>
{
  public:
  using type=typename OperatorTypeHelper<Multiplication<Matrix<T,Rows1,Cols1>,
                                                        Matrix<T,Rows2,Cols2>>,Ts...>::type;
};









///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////// Division between a matrix(M X N)and a scalar ( 1 x 1)                                       /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<Integer Rows,Integer Cols, typename...Ts>
class OperatorTypeHelper<Division<Matrix<Real,Rows,Cols>, 
                                  Matrix<Real,1,1> >, Ts...>
{
  public:
  using type=Matrix<Real,Rows,Cols>;
};

template<Integer Rows,Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper<Division<QPValues< Matrix<Real,Rows,Cols>, NQPoints>, 
                                  QPValues< Matrix<Real,1,1>, NQPoints> >, Ts...>
{
  public:
  using type=QPValues<Matrix<Real,Rows,Cols>,NQPoints>;
};

template<Integer Rows,Integer Cols, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper<Division<FQPValues< Matrix<Real,Rows,Cols>, NQPoints, Ndofs>, 
                                  FQPValues< Matrix<Real,1,1>, NQPoints, Ndofs> >, Ts...>
{
  public:
  using type=FQPValues<Matrix<Real,Rows,Cols>,NQPoints,Ndofs>;
};









template<Integer Rows,Integer Cols, typename...Ts>
class OperatorTypeHelper<Division<Transposed<Matrix<Real,Cols,Rows>>, 
                                  Matrix<Real,1,1> >, Ts...>
{
  public:
  using type=Matrix<Real,Rows,Cols>;
};


template<Integer Rows,Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper<Division<QPValues< Transposed<Matrix<Real,Cols,Rows>>, NQPoints>, 
                                  QPValues< Matrix<Real,1,1>, NQPoints> >, Ts...>
{
  public:
  using type=QPValues<Matrix<Real,Rows,Cols>,NQPoints>;
};

template<Integer Rows,Integer Cols, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper<Division<FQPValues< Transposed<Matrix<Real,Cols,Rows>>, NQPoints, Ndofs>, 
                                  FQPValues< Matrix<Real,1,1>, NQPoints, Ndofs> >, Ts...>
{
  public:
  using type=FQPValues<Matrix<Real,Rows,Cols>,NQPoints,Ndofs>;
};








template<Integer Rows,Integer Cols, typename...Ts>
class OperatorTypeHelper<Division<Matrix<Real,Rows,Cols>, 
                                  Transposed<Matrix<Real,1,1>> >, Ts...>
{
  public:
  using type=Matrix<Real,Rows,Cols>;
};

template<Integer Rows,Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper<Division<QPValues< Matrix<Real,Rows,Cols>, NQPoints>, 
                                  QPValues< Transposed<Matrix<Real,1,1>>, NQPoints> >, Ts...>
{
  public:
  using type=QPValues<Matrix<Real,Rows,Cols>,NQPoints>;
};

template<Integer Rows,Integer Cols, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper<Division<FQPValues< Matrix<Real,Rows,Cols>, NQPoints, Ndofs>, 
                                  FQPValues< Transposed<Matrix<Real,1,1>>, NQPoints, Ndofs> >, Ts...>
{
  public:
  using type=FQPValues<Matrix<Real,Rows,Cols>,NQPoints,Ndofs>;
};













template<Integer Rows,Integer Cols, typename...Ts>
class OperatorTypeHelper<Division<Transposed<Matrix<Real,Cols,Rows>>, 
                                  Transposed<Matrix<Real,1,1>> >, Ts...>
{
  public:
  using type=Matrix<Real,Rows,Cols>;
};

template<Integer Rows,Integer Cols, Integer NQPoints, typename...Ts>
class OperatorTypeHelper<Division<QPValues< Transposed<Matrix<Real,Cols,Rows>>, NQPoints>, 
                                  QPValues< Transposed<Matrix<Real,1,1>>, NQPoints> >, Ts...>
{
  public:
  using type=QPValues<Matrix<Real,Rows,Cols>,NQPoints>;
};

template<Integer Rows,Integer Cols, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper<Division<FQPValues< Transposed<Matrix<Real,Cols,Rows>>, NQPoints, Ndofs>, 
                                  FQPValues< Transposed<Matrix<Real,1,1>>, NQPoints, Ndofs> >, Ts...>
{
  public:
  using type=FQPValues<Matrix<Real,Rows,Cols>,NQPoints,Ndofs>;
};



















template<typename T, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< UnaryPlus< QPValues<T,NQPoints>>, Ts...>
{
  public:
  using type=QPValues<typename OperatorTypeHelper<UnaryPlus<T,Ts...>>::type ,NQPoints>;
  // using type=QPValues<typename OperatorTypeHelper<T,Ts...>::type ,NQPoints>;

};


template<typename T, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper< UnaryPlus< FQPValues<T,NQPoints,Ndofs>>, Ts...>
{
  public:
  using type=FQPValues<typename OperatorTypeHelper<UnaryPlus<T,Ts...>>::type ,NQPoints,Ndofs>;
  // using type=FQPValues<typename OperatorTypeHelper<T,Ts...>::type ,NQPoints,Ndofs>;

};


template<typename T, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< UnaryMinus< QPValues<T,NQPoints>>, Ts...>
{
  public:
  using type=QPValues<typename OperatorTypeHelper<UnaryMinus<T,Ts...>>::type ,NQPoints>;
};

template<typename T, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper< UnaryMinus< FQPValues<T,NQPoints,Ndofs>>, Ts...>
{
  public:
  using type=FQPValues<typename OperatorTypeHelper<UnaryMinus<T,Ts...>>::type ,NQPoints,Ndofs>;
};











// type(T) = type(+T)
template<typename T, typename...Ts>
class OperatorTypeHelper<T,Ts...>
{ public:
  using type=T;
};




template<typename T,Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper<Multiplication<Matrix<T,Rows,Cols>, 
                                        Vector<T,Cols> >, Ts...>
{
  public:
  using type=Vector<T,Rows>;
};




template<typename T, typename S, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Multiplication< QPValues<T,NQPoints>, 
                                          QPValues<S,NQPoints> >, Ts...>
{
  public:
  using type=QPValues<typename OperatorTypeHelper<Multiplication<T,S>,Ts...>::type ,NQPoints>;
};

template<typename T, typename S, Integer NQPoints, typename...Ts>
class OperatorTypeHelper< Division< QPValues<T,NQPoints>, 
                                    QPValues<S,NQPoints> >, Ts...>
{
  public:
  using type=QPValues<typename OperatorTypeHelper<Division<T,S>,Ts...>::type ,NQPoints>;
};



template<typename T, typename S, Integer NQPoints,Integer NComponents, typename...Ts>
class OperatorTypeHelper<Multiplication< QPValues<S,NQPoints>, 
                                         FQPValues<T,NQPoints,NComponents> >, Ts...>
{
  public:
  using type=FQPValues<typename OperatorTypeHelper<Multiplication<S,T>,Ts...>::type ,NQPoints,NComponents>;
};

template<typename T, typename S, Integer NQPoints,Integer NComponents, typename...Ts>
class OperatorTypeHelper<Division< QPValues<S,NQPoints>, 
                                   FQPValues<T,NQPoints,NComponents> >, Ts...>
{
  public:
  using type=FQPValues<typename OperatorTypeHelper<Division<S,T>,Ts...>::type ,NQPoints,NComponents>;
};

template<typename T, typename S, Integer NQPoints,Integer NComponents, typename...Ts>
class OperatorTypeHelper<Multiplication< FQPValues<S,NQPoints,NComponents>,
                                         QPValues<T,NQPoints> >, Ts...>
{
  public:
  using type=FQPValues<typename OperatorTypeHelper<Multiplication<S,T>,Ts...>::type ,NQPoints,NComponents>;
};

template<typename T, typename S, Integer NQPoints,Integer NComponents, typename...Ts>
class OperatorTypeHelper<Division< FQPValues<S,NQPoints,NComponents>,
                                   QPValues<T,NQPoints> >, Ts...>
{
  public:
  using type=FQPValues<typename OperatorTypeHelper<Division<S,T>,Ts...>::type ,NQPoints,NComponents>;
};


template<typename T, Integer NQPoints, typename...Ts>
class OperatorTypeHelper<Transposed<QPValues<T,NQPoints>>,Ts...>
{ public:
  using type=QPValues< typename OperatorTypeHelper<Transposed<T>,Ts...>::type ,NQPoints>;
};

template<typename T, Integer NQPoints, Integer Ndofs, typename...Ts>
class OperatorTypeHelper<Transposed<FQPValues<T,NQPoints,Ndofs>>,Ts...>
{ public:
  using type=FQPValues< typename OperatorTypeHelper<Transposed<T>,Ts...>::type ,NQPoints,Ndofs>;
};




template<typename T, typename...Ts>
class OperatorTypeHelper<Transposed<Transposed<T>>,Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};





// type(T) = type(+T)
template<typename T, typename...Ts>
class OperatorTypeHelper<UnaryPlus<Expression<T> >, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};





// type(T) = type(-T)
template<typename T, typename...Ts>
class OperatorTypeHelper<UnaryMinus<Expression<T>>, Ts...>
{ public:
  using type=typename OperatorTypeHelper<T,Ts...>::type;
};


template<typename T, typename...Ts>//template<class > class Unary, 
class OperatorTypeHelper<UnaryMinus<Expression<UnaryPlus<Expression<T>>>>,Ts...>
{ public:
  using type=typename OperatorTypeHelper<UnaryMinus<Expression<T>>,Ts...>::type;
};


// type(A)=type(B) = type(A+B)
template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Addition< Expression<Left>, Expression<Right > >, Ts...>
{ public:
  using type1=typename OperatorTypeHelper<Left,Ts...>::type;
  using type2=typename OperatorTypeHelper<Right,Ts...>::type;
  // static_assert(IsAddable<typename OperatorTypeHelper<Left>::type,
  //                         typename OperatorTypeHelper<Right>::type
  //                        >::value, " In Addition, Left and Right types must be equal");
  using type=typename OperatorTypeHelper<Addition<type1,type2>,Ts...>::type;
};


// type(A)=type(B) = type(A-B)
template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Subtraction< Expression<Left>, Expression<Right > >, Ts...>
{ public:
  // static_assert(IsSame<typename OperatorTypeHelper<Left>::type,
  //                      typename OperatorTypeHelper<Right>::type
  //                     >::value, " In Subtraction, Left and Right types must be equal");
  using type1=typename OperatorTypeHelper<Left,Ts...>::type;
  using type2=typename OperatorTypeHelper<Right,Ts...>::type;
  using type=typename OperatorTypeHelper<Subtraction<Left,Right>,Ts...>::type;
};

template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Multiplication< Expression<Left>, Expression<Right > >, Ts...>
{ public:
  // static_assert(IsSame<typename OperatorTypeHelper<Left>::type,
  //                      typename OperatorTypeHelper<Right>::type
  //                     >::value, " In Subtraction, Left and Right types must be equal");

  using LeftType=typename OperatorTypeHelper<Left,Ts...>::type;
  using RightType=typename OperatorTypeHelper<Right,Ts...>::type;
  using type=typename OperatorTypeHelper<Multiplication<LeftType,RightType>,Ts...>::type;
};


template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Division< Expression<Left>, Expression<Right > >, Ts...>
{ public:
  using LeftType=typename OperatorTypeHelper<Left,Ts...>::type;
  using RightType=typename OperatorTypeHelper<Right,Ts...>::type;
  using type=typename OperatorTypeHelper<Division<LeftType,RightType>,Ts...>::type;
};





///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////// UnaryPlus for Matrix, QPValues, FQPValues and their transposes                              /////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, Integer Rows, Integer Cols,typename...Ts>
class OperatorTypeHelper<UnaryPlus< Matrix<T,Cols,Rows> >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols,typename...Ts>
class OperatorTypeHelper<UnaryPlus< Transposed< Matrix<T,Cols,Rows> > >, Ts...>
{ public:
  using type=Matrix<T,Rows,Cols>;
};

template<typename T, Integer Rows, Integer Cols,Integer NQPoints, typename...Ts>
class OperatorTypeHelper<UnaryPlus< Transposed< QPValues< Matrix<T,Cols,Rows>,NQPoints> > >, Ts...>
{ public:
  using type=QPValues< Matrix<T,Rows,Cols>,NQPoints>;
};

template<typename T, Integer Rows, Integer Cols,Integer NQPoints,Integer Ndofs, typename...Ts>
class OperatorTypeHelper<UnaryPlus< Transposed< FQPValues< Matrix<T,Cols,Rows>,NQPoints,Ndofs> > >, Ts...>
{ public:
  using type=FQPValues< Matrix<T,Rows,Cols>,NQPoints,Ndofs >;
};






































// type(A)=type(B) = type(A+B)
template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Addition< Left, Right >, Ts...>
{ public:
  // we add UnaryPlus, so that we can deal with transposed or not equally
  // using type1=typename OperatorTypeHelper<UnaryPlus<Left>,Ts...>::type;
  // using type2=typename OperatorTypeHelper<UnaryPlus<Right>,Ts...>::type;
  // using type=typename OperatorTypeHelper<Left,Ts...>::type;
 using subtypeleft=OperatorType<Left,Ts...>;
 using subtyperight=OperatorType<Right,Ts...>;
 using type=OperatorType<Addition<subtypeleft,subtyperight>,Ts...>;

  // using type=typename OperatorTypeHelper<UnaryPlus<Left>,Ts...>::type;

};

// type(A)=type(B) = type(A+B)
template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Subtraction< Left, Right >, Ts...>
{ public:
  // we add UnaryPlus, so that we can deal with transposed or not equally
  // using type1=typename OperatorTypeHelper<UnaryPlus<Left>,Ts...>::type;
  // using type2=typename OperatorTypeHelper<UnaryPlus<Right>,Ts...>::type;
  // using type=typename OperatorTypeHelper<Left,Ts...>::type;
 using subtypeleft=OperatorType<Left,Ts...>;
 using subtyperight=OperatorType<Right,Ts...>;
 using type=OperatorType<Subtraction<subtypeleft,subtyperight>,Ts...>;


  // using type=typename OperatorTypeHelper<UnaryPlus<Left>,Ts...>::type;
};













// template<typename T, typename...Ts>
// class OperatorTypeHelper<UnaryPlus<T>, Ts...>
// { public:
//   using type=typename OperatorTypeHelper<T,Ts...>::type;
// };


// type(T) = type(-T)
// template<typename T, typename...Ts>
// class OperatorTypeHelper<UnaryMinus<T>, Ts...>
// { public:
//   using type=typename OperatorTypeHelper<T,Ts...>::type;
// };


















// template<typename Left, typename Right, typename...Ts>
// class OperatorTypeHelper< Addition<Left,Right>, Ts...>
// { public:
//   using type=typename OperatorTypeHelper<Left,Ts...>::type;
// };


// // type(A)=type(B) = type(A-B)
// template<typename Left, typename Right, typename...Ts>
// class OperatorTypeHelper< Subtraction< Left, Right >, Ts...>
// { public:
//   using type=typename OperatorTypeHelper<Left,Ts...>::type;
// };

template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Multiplication< Left, Right >, Ts...>
{ public:
  using LeftType=typename OperatorTypeHelper<Left,Ts...>::type;
  using RightType=typename OperatorTypeHelper<Right,Ts...>::type;
  using type=typename OperatorTypeHelper<Multiplication<LeftType,RightType>,Ts...>::type;
};


template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Division< Left, Right >, Ts...>
{ public:
  using LeftType=typename OperatorTypeHelper<Left,Ts...>::type;
  using RightType=typename OperatorTypeHelper<Right,Ts...>::type;
  using type=typename OperatorTypeHelper<Division<LeftType,RightType>,Ts...>::type;
};










template<typename ConstType, typename...Inputs, typename...Ts>
class OperatorTypeHelper<Transposed<Expression<ConstantTensor<ConstType,Inputs...>>>,Ts...>
{ public:
  using type_tmp=typename OperatorTypeHelper<ConstantTensor<ConstType,Inputs...>,Ts...>::type;
  using type=typename OperatorTypeHelper<Transposed<type_tmp>,Ts...>::type;
};

template<typename FuncType,Integer N,typename OperatorType,typename FullSpace, typename...Ts>
class OperatorTypeHelper<Transposed<Expression<Function<FullSpace,N,OperatorType,FuncType>>>,Ts...>
{ public:
  using type_tmp=typename OperatorTypeHelper<Function<FullSpace,N,OperatorType,FuncType>,Ts...>::type;
  using type=typename OperatorTypeHelper<Transposed<type_tmp>,Ts...>::type;
};

template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename OperatorType,typename QRule>
class OperatorTypeHelper<Transposed<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,QRule >
{ public:
 static_assert((IsSame<TestOrTrial<MixedSpace,N,OperatorType>,Test<MixedSpace,N,OperatorType>>::value ||
                IsSame<TestOrTrial<MixedSpace,N,OperatorType>,Trial<MixedSpace,N,OperatorType>>::value )
               && "In OperatorTypeHelper<Transposed<TestOrTrial>>,TestOrTrial=Test or Trial ");  
  using tmptype= TestOrTrial<MixedSpace,N,OperatorType>;
  using FunctionSpace=typename tmptype::FunctionSpace;
  using FunctionSpaces=GetType<typename tmptype::UniqueElementFunctionSpacesTupleType,tmptype::value>;
  using Elem=typename FunctionSpaces::Elem;
  using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
  using Operator=typename tmptype::Operator; 
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using type=typename OperatorTypeHelper<Transposed<typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QRule>::type>,QRule>::type;
};




template<template<class...> class Nnary, typename...Tss, typename...Ts>
class OperatorTypeHelper<Transposed<Expression<Nnary<Tss...>>>,Ts...>
{ public:
  using type_tmp=typename OperatorTypeHelper<Nnary<Tss...>,Ts...>::type;
  using type=typename OperatorTypeHelper<Transposed<type_tmp>,Ts...>::type;
};

template<typename T, typename...Ts>
class OperatorTypeHelper<Transposed<Expression<T>>,Ts...>
{ public:
  using type_tmp=typename OperatorTypeHelper<Expression<T>,Ts...>::type;
  using type=typename OperatorTypeHelper<Transposed<type_tmp>,Ts...>::type;
};


// template<template<class > class Unary, typename T, typename...Ts>
// class OperatorTypeHelper<Unary<Transposed<T>>,Ts...>
// { public:
//   using type=typename OperatorTypeHelper<Unary<T>,Ts...>::type;
// };

template<typename T, typename...Ts>//template<class > class Unary, 
class OperatorTypeHelper<UnaryPlus<Expression<Transposed<Expression<T>>>>,Ts...>
{ public:
  using type_tmp=typename OperatorTypeHelper<UnaryPlus<Expression<T>>,Ts...>::type;
  using type=typename OperatorTypeHelper<UnaryPlus<Transposed<type_tmp>>,Ts...>::type;
};
template<typename T, typename...Ts>//template<class > class Unary, 
class OperatorTypeHelper<UnaryMinus<Expression<Transposed<Expression<T>>>>,Ts...>
{ public:
  using type_tmp=typename OperatorTypeHelper<UnaryMinus<Expression<T>>,Ts...>::type;
  using type=typename OperatorTypeHelper<UnaryPlus<Transposed<type_tmp>>,Ts...>::type;
};

template<typename Left,typename Right, typename...Ts>
class OperatorTypeHelper<Transposed<Expression<Addition<Expression<Left>,Expression<Right>>>>,Ts...>
{ public:
  using type_tmp=typename OperatorTypeHelper<Addition<Expression<Left>,Expression<Right>>,Ts...>::type;
  using type=typename OperatorTypeHelper<Transposed<type_tmp>,Ts...>::type;
};
template<typename Left,typename Right, typename...Ts>
class OperatorTypeHelper<Transposed<Expression<Subtraction<Expression<Left>,Expression<Right>>>>,Ts...>
{ public:
  using type_tmp=typename OperatorTypeHelper<Subtraction<Expression<Left>,Expression<Right>>,Ts...>::type;
  using type=typename OperatorTypeHelper<Transposed<type_tmp>,Ts...>::type;
};
template<typename Left,typename Right, typename...Ts>
class OperatorTypeHelper<Transposed<Expression<Multiplication<Expression<Left>,Expression<Right>>>>,Ts...>
{ public:
  using type_tmp=typename OperatorTypeHelper<Multiplication<Expression<Left>,Expression<Right>>,Ts...>::type;
  using type=typename OperatorTypeHelper<Transposed<type_tmp>,Ts...>::type;
};
template<typename Left,typename Right, typename...Ts>
class OperatorTypeHelper<Transposed<Expression<Division<Expression<Left>,Expression<Right>>>>,Ts...>
{ public:
  using type_tmp=typename OperatorTypeHelper<Division<Expression<Left>,Expression<Right>>,Ts...>::type;
  using type=typename OperatorTypeHelper<Transposed<type_tmp>,Ts...>::type;
};




// template<template<class,class> class Binary, typename Left,typename Right, typename...Ts>
// class OperatorTypeHelper< Transposed<Expression<Binary<Expression<Left>,Expression<Right> > > >, Ts...>
// { public:
//   using type=typename OperatorTypeHelper< Binary<Expression<Transposed<Expression<Left>>>,Expression<Transposed<Expression<Right>>>>,Ts...>::type;
// };

// template<typename Left,typename Right, typename...Ts>
// class OperatorTypeHelper< Transposed<Expression<Multiplication<Expression<Left>,Expression<Right> > > >, Ts...>
// { public:
//   using type=typename OperatorTypeHelper< Multiplication<Expression<Transposed<Expression<Left>>>,Expression<Transposed<Expression<Right>>>>,Ts...>::type;
// };


// template<typename T, typename...Ts>
// class OperatorTypeHelper<Transposed<T>,Ts...>
// { public:
//   using type=typename OperatorTypeHelper<T,Ts...>::type;
// };

// template<typename ConstType, typename...Inputs, typename...Ts>
// class OperatorTypeHelper<Transposed<Expression<ConstantTensor<ConstType,Inputs...>>>,Ts...>
// { public:
//   using type_tmp=typename OperatorTypeHelper<ConstantTensor<ConstType,Inputs...>,Ts...>::type;
//   using type=typename OperatorTypeHelper<Transposed<type_tmp>,Ts...>::type;
// };

// template<typename FuncType,Integer N,typename OperatorType,typename FullSpace, typename...Ts>
// class OperatorTypeHelper<Transposed<Expression<Function<FullSpace,N,OperatorType,FuncType>>>,Ts...>
// { public:
//   using type_tmp=typename OperatorTypeHelper<Function<FullSpace,N,OperatorType,FuncType>,Ts...>::type;
//   using type=typename OperatorTypeHelper<Transposed<type_tmp>,Ts...>::type;
// };


// template<typename T, typename...Ts>
// class OperatorTypeHelper<Transposed<Expression<T>>,Ts...>
// { public:
//   using type=Transposed<typename OperatorTypeHelper<T,Ts...>::type>;
// };








template<typename ConstType,typename...Inputs, typename Elem,typename BaseFunctionSpace, typename QuadratureRule_>
class OperatorTypeHelper< CompositeOperator<Expression< ConstantTensor<ConstType,Inputs...> > >,
                          Elem,
                          BaseFunctionSpace,
                          QuadratureRule_>
{ public:

  // type is simply T, but we must create QPValues<T,NQPoints>, so we use qpvalues<NQPoints>
  using type= typename ConstantTensor<ConstType,Inputs...>::template qptype<QuadratureRule_::NQPoints>;
};

template<typename FullSpace,Integer N,typename Operator_,typename FuncType, typename Elem,typename BaseFunctionSpace, typename QuadratureRule_>
class OperatorTypeHelper< CompositeOperator<Expression< Function<FullSpace,N,Operator_,FuncType> > >,
                          Elem,
                          BaseFunctionSpace,
                          QuadratureRule_>
{ public:

  // TODO FIXME
  // using type= typename ConstantTensor<ConstType,Inputs...>::type;
};

template<typename T, typename Elem,typename BaseFunctionSpace, typename QuadratureRule_>
class OperatorTypeHelper< CompositeOperator<Expression< T > >,
                          Elem,
                          BaseFunctionSpace,
                          QuadratureRule_>
{ public:

  // TODO FIXME
  using type= typename ShapeFunction<Elem,BaseFunctionSpace,T,QuadratureRule_>::type;
};



 template<typename T, typename Elem,typename BaseFunctionSpace, typename QuadratureRule_>
class OperatorTypeHelper< CompositeOperator<Expression< UnaryPlus< Expression< T > > > >,
                          Elem,
                          BaseFunctionSpace,
                          QuadratureRule_>
{ public:

  using type= typename OperatorTypeHelper<CompositeOperator<Expression<T>>,
                                              Elem,
                                              BaseFunctionSpace,
                                              QuadratureRule_
                                              >::type;
};

 template<typename T, typename Elem,typename BaseFunctionSpace, typename QuadratureRule_>
class OperatorTypeHelper< CompositeOperator<Expression< UnaryMinus< Expression< T > > > >,
                          Elem,
                          BaseFunctionSpace,
                          QuadratureRule_>
{ public:

  using type= typename OperatorTypeHelper<CompositeOperator<Expression<T>>,
                                              Elem,
                                              BaseFunctionSpace,
                                              QuadratureRule_
                                              >::type;
};




 template<typename Left, typename Right, typename Elem,typename BaseFunctionSpace, typename QuadratureRule_>
class OperatorTypeHelper< CompositeOperator<Expression< Division< Expression< Left >,Expression< Right > > > >,
                          Elem,
                          BaseFunctionSpace,
                          QuadratureRule_>
{ public:
  using typeLeft= typename OperatorTypeHelper<CompositeOperator<Expression<Left>>,
                                              Elem,
                                              BaseFunctionSpace,
                                              QuadratureRule_
                                              >::type;
  using typeRight= typename OperatorTypeHelper<CompositeOperator<Expression<Right>>,
                                              Elem,
                                              BaseFunctionSpace,
                                              QuadratureRule_
                                              >::type;
  using type=typename OperatorTypeHelper< Division<Expression<typeLeft>,Expression<typeRight>>
                                        >::type;
};

 template<typename Left, typename Right, typename Elem,typename BaseFunctionSpace, typename QuadratureRule_>
class OperatorTypeHelper< CompositeOperator<Expression< Multiplication< Expression< Left >,Expression< Right > > > >,
                          Elem,
                          BaseFunctionSpace,
                          QuadratureRule_>
{ public:
  using typeLeft= typename OperatorTypeHelper<CompositeOperator<Expression<Left>>,
                                              Elem,
                                              BaseFunctionSpace,
                                              QuadratureRule_
                                              >::type;
  using typeRight= typename OperatorTypeHelper<CompositeOperator<Expression<Right>>,
                                              Elem,
                                              BaseFunctionSpace,
                                              QuadratureRule_
                                              >::type;
  using type=typename OperatorTypeHelper< Multiplication<Expression<typeLeft>,Expression<typeRight>>
                                        >::type;
};


 template<typename Left, typename Right, typename Elem,typename BaseFunctionSpace, typename QuadratureRule_>
class OperatorTypeHelper< CompositeOperator<Expression< Addition< Expression< Left >,Expression< Right > > > >,
                          Elem,
                          BaseFunctionSpace,
                          QuadratureRule_>
{ public:
  using typeLeft= typename OperatorTypeHelper<CompositeOperator<Expression<Left>>,
                                              Elem,
                                              BaseFunctionSpace,
                                              QuadratureRule_
                                              >::type;
  using typeRight= typename OperatorTypeHelper<CompositeOperator<Expression<Right>>,
                                              Elem,
                                              BaseFunctionSpace,
                                              QuadratureRule_
                                              >::type;
  using type=typename OperatorTypeHelper< Addition<Expression<typeLeft>,Expression<typeRight>>
                                        >::type;
};

 template<typename Left, typename Right, typename Elem,typename BaseFunctionSpace, typename QuadratureRule_>
class OperatorTypeHelper< CompositeOperator<Expression< Subtraction< Expression< Left >,Expression< Right > > > >,
                          Elem,
                          BaseFunctionSpace,
                          QuadratureRule_>
{ public:
  using typeLeft= typename OperatorTypeHelper<CompositeOperator<Expression<Left>>,
                                              Elem,
                                              BaseFunctionSpace,
                                              QuadratureRule_
                                              >::type;
  using typeRight= typename OperatorTypeHelper<CompositeOperator<Expression<Right>>,
                                              Elem,
                                              BaseFunctionSpace,
                                              QuadratureRule_
                                              >::type;
  using type=typename OperatorTypeHelper< Subtraction<Expression<typeLeft>,Expression<typeRight>>
                                        >::type;
};




template<typename...Ts>
class OperatorTypeHelper< Contraction2< Expression<Real>, Expression<Real > >, Ts...>
{ public:
  using type=Real;
};

template<typename T, Integer Dim, typename...Ts>
class OperatorTypeHelper< Contraction2< Expression<Vector<T,Dim>>, Expression<Vector<T,Dim>> >, Ts...>
{ public:
  using type=T;
};

template<typename T, Integer Rows, Integer Cols, typename...Ts>
class OperatorTypeHelper< Contraction2< Expression<Matrix<T,Rows,Cols>>, Expression<Matrix<T,Rows,Cols>> >, Ts...>
{ public:
  using type=T;
};



template<typename Left, typename Right, typename...Ts>
class OperatorTypeHelper< Contraction2< Expression<Left>, Expression<Right > >, Ts...>
{ public:
  using typeleft=typename OperatorTypeHelper<Left,Ts...>::type;
  using typeright=typename OperatorTypeHelper<Right,Ts...>::type;
  using type=typename OperatorTypeHelper<Contraction2<typeleft,typeright>,Ts...>::type;
};






























template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename OperatorType,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,OperatorType>,QRule >
{ public:
 static_assert((IsSame<TestOrTrial<MixedSpace,N,OperatorType>,Test<MixedSpace,N,OperatorType>>::value ||
                IsSame<TestOrTrial<MixedSpace,N,OperatorType>,Trial<MixedSpace,N,OperatorType>>::value )
               && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  using tmptype= TestOrTrial<MixedSpace,N,OperatorType>;
  using FunctionSpace=typename tmptype::FunctionSpace;
  using FunctionSpaces=GetType<typename tmptype::UniqueElementFunctionSpacesTupleType,tmptype::value>;
  using Elem=typename FunctionSpaces::Elem;
  using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
  using Operator=typename tmptype::Operator; 
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using type=typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QRule>::type;
};


template<typename FuncType,Integer N,typename OperatorType,typename FullSpace, typename QRule>
class OperatorTypeHelper<Function<FullSpace,N,OperatorType,FuncType>,QRule >
{ public:
  using tmptype=Function<FullSpace,N,OperatorType,FuncType>;
  using FunctionSpace=typename tmptype::FunctionSpace;
  using FunctionSpaces=GetType<typename tmptype::UniqueElementFunctionSpacesTupleType,tmptype::value>;
  using Elem=typename FunctionSpaces::Elem;
  using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
  using Operator=typename tmptype::Operator; 
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using type=typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QRule>::qpvalues_type;
};

template<typename ConstType,typename...Inputs, typename QRule>
class OperatorTypeHelper<ConstantTensor<ConstType,Inputs...>,QRule >
{ public:
  using Constant=ConstantTensor<ConstType,Inputs...>;
  using tmptype=typename Constant::type;
  using type=QPValues<tmptype,QRule::NQPoints>;
};







template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename T,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< T >>>,QRule >
{ public:
  using Operator=CompositeOperator<Expression<UnaryPlus<Expression<T>>>>;
  using tmptype= TestOrTrial<MixedSpace,N,Operator>;
  static_assert((IsSame<tmptype,Test<MixedSpace,N,Operator>>::value ||
                 IsSame<tmptype,Trial<MixedSpace,N,Operator>>::value )
                && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using type=typename OperatorTypeHelper<TestOrTrial<MixedSpace,N,T>,QRule>::type;
};



template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename FuncType,Integer M,typename OperatorType,typename FullSpace,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< Function<FullSpace,M,OperatorType,FuncType> >>>,QRule >
{ public:
  using tmptype=Function<FullSpace,M,OperatorType,FuncType>;
  using FunctionSpace=typename tmptype::FunctionSpace;
  using FunctionSpaces=GetType<typename tmptype::UniqueElementFunctionSpacesTupleType,tmptype::value>;
  using Elem=typename FunctionSpaces::Elem;
  using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
  using Operator=typename tmptype::Operator; 
  using type=typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QRule>::qpvalues_type;
};

template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename FuncType,Integer M,typename OperatorType,typename FullSpace,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Transposed<Expression< Function<FullSpace,M,OperatorType,FuncType> > >>>>,QRule >
{ public:
  using tmptype=Function<FullSpace,M,OperatorType,FuncType>;
  using FunctionSpace=typename tmptype::FunctionSpace;
  using FunctionSpaces=GetType<typename tmptype::UniqueElementFunctionSpacesTupleType,tmptype::value>;
  using Elem=typename FunctionSpaces::Elem;
  using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
  using Operator=typename tmptype::Operator; 
  using typetmp=typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QRule>::qpvalues_type;
  using type=typename OperatorTypeHelper<Transposed<typetmp>,QRule>::type;
};


template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename ConstType,typename...Inputs,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< ConstantTensor<ConstType,Inputs...> >>>,QRule >
{ public:
  using type=typename OperatorTypeHelper<ConstantTensor<ConstType,Inputs...>,QRule>::type;
};

template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename ConstType,typename...Inputs,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Transposed<Expression< ConstantTensor<ConstType,Inputs...> >>>>>,QRule >
{ public:
  using typetmp=typename OperatorTypeHelper<ConstantTensor<ConstType,Inputs...>,QRule>::type;
  using type=typename OperatorTypeHelper<Transposed<typetmp>,QRule>::type;

};



template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename T,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< Transposed<Expression<T>> >>>,QRule >
{ public:
  using Operator=CompositeOperator<Expression<Transposed<Expression<T>>>>;
  using tmptype= TestOrTrial<MixedSpace,N,Operator>;
 static_assert((IsSame<tmptype,Test<MixedSpace,N,Operator>>::value ||
                IsSame<tmptype,Trial<MixedSpace,N,Operator>>::value )
               && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using type=typename OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<T>>>>::type;
};


template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename T,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< UnaryPlus<Expression<T>> >>>,QRule >
{ public:
  using Operator=CompositeOperator<Expression<UnaryPlus<Expression<T>>>>;
  using tmptype= TestOrTrial<MixedSpace,N,Operator>;
 static_assert((IsSame<tmptype,Test<MixedSpace,N,Operator>>::value ||
                IsSame<tmptype,Trial<MixedSpace,N,Operator>>::value )
               && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using type=typename OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<T>>>>::type;
};

template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename T,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< UnaryMinus<Expression<T>> >>>,QRule >
{ public:
  using Operator=CompositeOperator<Expression<UnaryMinus<Expression<T>>>>;
  using tmptype= TestOrTrial<MixedSpace,N,Operator>;
 static_assert((IsSame<tmptype,Test<MixedSpace,N,Operator>>::value ||
                IsSame<tmptype,Trial<MixedSpace,N,Operator>>::value )
               && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using type=typename OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<T>>>>::type;
};


template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename Left, typename Right,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< Addition<Expression<Left>,Expression<Right>> >>>,QRule >
{ public:
  using Operator=CompositeOperator<Expression<Addition<Expression<Left>,Expression<Right>>>>;
  using tmptype= TestOrTrial<MixedSpace,N,Operator>;
 static_assert((IsSame<tmptype,Test<MixedSpace,N,Operator>>::value ||
                IsSame<tmptype,Trial<MixedSpace,N,Operator>>::value )
               && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using TestOrTrialLeft=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Left>>>;
  using TestOrTrialRight=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Right>>>;
  using type=typename OperatorTypeHelper<Addition<Expression<TestOrTrialLeft>,Expression<TestOrTrialRight>>,QRule>::type;//typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QRule>::type;
};

template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename Left, typename Right,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< Subtraction<Expression<Left>,Expression<Right>> >>>,QRule >
{ public:
  using Operator=CompositeOperator<Expression<Subtraction<Expression<Left>,Expression<Right>>>>;
  using tmptype= TestOrTrial<MixedSpace,N,Operator>;
 static_assert((IsSame<tmptype,Test<MixedSpace,N,Operator>>::value ||
                IsSame<tmptype,Trial<MixedSpace,N,Operator>>::value )
               && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using TestOrTrialLeft=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Left>>>;
  using TestOrTrialRight=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Right>>>;
  using type=typename OperatorTypeHelper<Subtraction<Expression<TestOrTrialLeft>,Expression<TestOrTrialRight>>,QRule>::type;//typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QRule>::type;
};

template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename Left, typename Right,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< Multiplication<Expression<Left>,Expression<Right>> >>>,QRule >
{ public:
  using Operator=CompositeOperator<Expression<Multiplication<Expression<Left>,Expression<Right>>>>;
  using tmptype= TestOrTrial<MixedSpace,N,Operator>;
 static_assert((IsSame<tmptype,Test<MixedSpace,N,Operator>>::value ||
                IsSame<tmptype,Trial<MixedSpace,N,Operator>>::value )
               && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using TestOrTrialLeft=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Left>>>;
  using TestOrTrialRight=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Right>>>;
  using type=typename OperatorTypeHelper<Multiplication<Expression<TestOrTrialLeft>,Expression<TestOrTrialRight>>,QRule>::type;//typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QRule>::type;
};


template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename Left, typename Right,typename QRule>
class OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< Division<Expression<Left>,Expression<Right>> >>>,QRule >
{ public:
  using Operator=CompositeOperator<Expression<Division<Expression<Left>,Expression<Right>>>>;
  using tmptype= TestOrTrial<MixedSpace,N,Operator>;
 static_assert((IsSame<tmptype,Test<MixedSpace,N,Operator>>::value ||
                IsSame<tmptype,Trial<MixedSpace,N,Operator>>::value )
               && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using TestOrTrialLeft=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Left>>>;
  using TestOrTrialRight=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Right>>>;
  using type=typename OperatorTypeHelper<Division<Expression<TestOrTrialLeft>,Expression<TestOrTrialRight>>,QRule>::type;//typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QRule>::type;
};



template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, 
         Integer N, typename T,typename QRule>
class OperatorTypeHelper<Transposed<Expression<TestOrTrial<MixedSpace,N,CompositeOperator<Expression< T >>>>>,QRule >
{ public:
  using Operator=CompositeOperator<Expression<UnaryPlus<Expression<T>>>>;
  using tmptype= TestOrTrial<MixedSpace,N,Operator>;
  static_assert((IsSame<tmptype,Test<MixedSpace,N,Operator>>::value ||
                 IsSame<tmptype,Trial<MixedSpace,N,Operator>>::value )
                && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  // TODO FIX ME SUBTYPE AND NOT TYPE CHECK (now in fwpvalues<T,M,N>, tot_type=T)
  using typetmp=typename OperatorTypeHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<T>>>,QRule>::type;
  using type=typename OperatorTypeHelper<Transposed<typetmp>,QRule>::type;
};


template<typename Left,typename Right,Integer QR>
class OperatorTypeHelper<L2DotProductIntegral<Left,Right,QR>,Number<1>>
{ public:

   using L2=L2DotProductIntegral<Left,Right,QR>;
   using FunctionSpace=typename L2::FunctionSpace;
   using TestTrialNumbers=typename L2::TestTrialNumbers;
   static constexpr Integer TestNumber=GetType<TestTrialNumbers>::value;
   static constexpr Integer TestN=FunctionSpace::Nelem_dofs_array[TestNumber];
   using type=Vector<Real,TestN >;
 };

template<typename Left,typename Right,Integer QR>
class OperatorTypeHelper<L2DotProductIntegral<Left,Right,QR>,Number<2>>
{ public:

   using L2=L2DotProductIntegral<Left,Right,QR>;
   using FunctionSpace=typename L2::FunctionSpace;
   using TestTrialNumbers=typename L2::TestTrialNumbers;
   static constexpr Integer TestNumber=GetType<TestTrialNumbers,0>::value;
   static constexpr Integer TrialNumber=GetType<TestTrialNumbers,1>::value;
   static constexpr Integer TestN=FunctionSpace::Nelem_dofs_array[TestNumber];
   static constexpr Integer TrialN=FunctionSpace::Nelem_dofs_array[TrialNumber];
   using type=Matrix<Real,TestN,TrialN >;//typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QRule>::type;
};







































template<typename...Ts>
class FormOfCompositeOperatorAuxType;


template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Operator>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,Operator>,Expr>
{
public:
  using T=TestOrTrial<MixedSpace,N,Operator>;
  using FunctionSpace=typename T::FunctionSpace;
  using FromElementFunctionSpacesToFirstSpaceTupleType=typename FunctionSpace::FromElementFunctionSpacesToFirstSpaceTupleType;
  static constexpr Integer FirstSpace=GetType<FromElementFunctionSpacesToFirstSpaceTupleType,T::value>::value;
  using type=Test<MixedSpace,FirstSpace,Operator>;
};

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Operator>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Operator>
{
public:
  using T=TestOrTrial<MixedSpace,N,Operator>;
  using FunctionSpace=typename T::FunctionSpace;
  using FromElementFunctionSpacesToFirstSpaceTupleType=typename FunctionSpace::FromElementFunctionSpacesToFirstSpaceTupleType;
  static constexpr Integer FirstSpace=GetType<FromElementFunctionSpacesToFirstSpaceTupleType,T::value>::value;
  using type=Test<MixedSpace,FirstSpace,Operator>;
};

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename ConstType,typename...Inputs>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,ConstantTensor<ConstType,Inputs...>>
{
public:
  using type=ConstantTensor<ConstType,Inputs...>;
};

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename FullSpace, Integer M,typename Operator,typename FuncType>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Function<FullSpace,M,Operator,FuncType>>
{
public:
  using type=Function<FullSpace,M,Operator,FuncType>;
};

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename T>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,UnaryPlus<Expression<T>>>
{
public:
  using TT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,T>::type;
  using type=UnaryPlus<Expression<TT>>;
};

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename T>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,UnaryMinus<Expression<T>>>
{
public:
  using TT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,T>::type;
  using type=UnaryMinus<Expression<TT>>;
};

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Addition<Expression<Left>,Expression<Right>>>
{
public:
  using LeftT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Left>::type;
  using RightT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Right>::type;
  using type=Addition<Expression<LeftT>,Expression<RightT>>;
};

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Subtraction<Expression<Left>,Expression<Right>>>
{
public:
  using LeftT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Left>::type;
  using RightT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Right>::type;
  using type=Subtraction<Expression<LeftT>,Expression<RightT>>;
};

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Multiplication<Expression<Left>,Expression<Right>>>
{
public:
  using LeftT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Left>::type;
  using RightT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Right>::type;
  using type=Multiplication<Expression<LeftT>,Expression<RightT>>;
};

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Division<Expression<Left>,Expression<Right>>>
{
public:
  using LeftT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Left>::type;
  using RightT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Right>::type;
  using type=Division<Expression<LeftT>,Expression<RightT>>;
};


template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename T>
class FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Transposed<Expression<T>>>
{
public:
  using TT=typename FormOfCompositeOperatorAuxType<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>,T>::type;
  using type=Transposed<Expression<TT>>;

};

template<typename...Ts>
class FormOfCompositeOperatorType;

template<template<class,Integer,class > class TestOrTrial_, typename MixedSpace,Integer N,typename Expr>
class FormOfCompositeOperatorType<TestOrTrial_< MixedSpace,N,CompositeOperator<Expression<Expr>> >>
{
public:
  using type=typename FormOfCompositeOperatorAuxType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Expr>>>,Expr>::type;
};









// template<typename FullSpace, typename Eval,Integer M>
// class TupleOfCompositeTensorAuxAuxAux;


// template<typename FullSpace, typename Expr,typename QuadratureRule,Integer M>
// class TupleOfCompositeTensorAuxAuxAux<FullSpace,Evaluation<Expression<Expr>,QuadratureRule>,M>
// {
//   public:

// // Test<MixedSpace,M,CompositeOperator<Expression<Expr>>>
//   // using type=std::tuple<typename OperatorType<Expr,QuadratureRule>::type>;
// // 
//   using type=std::tuple<typename OperatorType<Test<FullSpace,M,CompositeOperator<Expression<Expr>>>,QuadratureRule>::type>;
// };

// template<typename FullSpace, Integer M>
// class TupleOfCompositeTensorAuxAuxAux<FullSpace,std::tuple<>,M>
// {
//   public:
//   using type=std::tuple<>;
// };


// template<typename FullSpace, typename Tuple,typename BuildTuple,Integer M,Integer Nmax,Integer N>
// class TupleOfCompositeTensorAuxAux;

// template<typename FullSpace, typename BuildTuple,Integer M,Integer Nmax,Integer N>
// class TupleOfCompositeTensorAuxAux<FullSpace,std::tuple<>,BuildTuple,M,Nmax,N>
// {public:
//  using type=BuildTuple;};


// template<typename FullSpace, typename Tuple,typename BuildTuple,Integer M,Integer Nmax>
// class TupleOfCompositeTensorAuxAux<FullSpace,Tuple,BuildTuple,M,Nmax,Nmax>
// {
//  public:
//   using TupleNth=GetType<Tuple,Nmax> ;
//   using single_type= typename TupleOfCompositeTensorAuxAuxAux<FullSpace,TupleNth,M>::type ;
//   using ChangeType=TupleCatType<GetType<BuildTuple,M>,single_type>;
//   using type=TupleChangeType<M, ChangeType, BuildTuple> ;

// };


// template<typename FullSpace, typename Tuple,typename BuildTuple,Integer M,Integer Nmax,Integer N>
// class TupleOfCompositeTensorAuxAux
// {
//   public:
//   using TupleNth=GetType<Tuple,N> ;
//   using single_type= typename TupleOfCompositeTensorAuxAuxAux<FullSpace,TupleNth,M>::type ;
//   using ChangeType=TupleCatType<GetType<BuildTuple,M>,single_type>;
//   using BuildTupleNew=TupleChangeType<M, ChangeType, BuildTuple> ;
//   using type=typename TupleOfCompositeTensorAuxAux<FullSpace,Tuple,BuildTupleNew,M,Nmax,N+1>::type;

// };







// template<typename FullSpace,typename Tuple,typename BuildTuple,Integer Nmax,Integer N>
// class TupleOfCompositeTensorAux;

// template<typename FullSpace,typename Tuple,typename BuildTuple,Integer Nmax>
// class TupleOfCompositeTensorAux<FullSpace,Tuple,BuildTuple,Nmax,Nmax>
// {
//  public:
//   using TupleNth=GetType<Tuple,Nmax> ;
//   static constexpr Integer Nmax_aux=TupleTypeSize<TupleNth>::value-1;
//   using type=typename TupleOfCompositeTensorAuxAux<FullSpace,TupleNth,BuildTuple,Nmax,Nmax,0>::type;//
// };


// template<typename FullSpace, typename Tuple,typename BuildTuple,Integer Nmax,Integer N>
// class TupleOfCompositeTensorAux
// {
//   public:
//   using TupleNth=GetType<Tuple,N> ;
//   static constexpr Integer Nmax_aux=TupleTypeSize<TupleNth>::value-1;
//   using NewBuildTuple=typename TupleOfCompositeTensorAuxAux<FullSpace,TupleNth,BuildTuple,N,Nmax,0>::type;
//   using type=typename TupleOfCompositeTensorAux<FullSpace,Tuple,NewBuildTuple,Nmax,N+1>::type;
// };

// template<typename FullSpace, typename Tuple>
// class TupleOfCompositeTensorHelper
// {
//  public:
//   static constexpr Integer Nmax=TupleTypeSize<Tuple>::value;
//   using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
//   using type=typename TupleOfCompositeTensorAux<FullSpace,Tuple,emptytuple,Nmax,0>::type;
// };

// template<typename FullSpace,typename Tuple>
// using TupleOfCompositeTensor=typename TupleOfCompositeTensorHelper<FullSpace,Tuple>::type;










template<typename ...Ts>
class TupleOfCombinationFunctionsAuxAux;

template<typename QuadratureRule, typename Tuple,typename TupleTensor,typename Other>
class TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Other>
{ 
 public:
 using type=Tuple ;
 using type_tensor=TupleTensor;
};



template<typename QuadratureRule, typename Tuple,typename TupleTensor,
         template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N, typename Expr>
class TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Expr>>>>
{
  public:
  using Test=Test<MixedSpace,N,CompositeOperator<Expression<Expr>>>;
  using CompositeType=typename FormOfCompositeOperatorType<Test>::type;
  using EvaluationCompositeType=Evaluation<Expression<CompositeType>,QuadratureRule>;
  using TupleNth=GetType<Tuple,Test::value>;

  static constexpr Integer IsNotAlreadyThere=IsNegative(TypeToTupleElementPosition<EvaluationCompositeType,TupleNth>::value);

  using ChangeType=typename std::conditional<IsNotAlreadyThere,TupleCatType<TupleNth,std::tuple<EvaluationCompositeType>>,TupleNth>::type;

  using type= TupleChangeType<Test::value, ChangeType, Tuple> ;

  using TensorType=OperatorType<Test,QuadratureRule>;
  using TupleTensorNth=GetType<TupleTensor,Test::value>;
  using ChangeTensorType=typename std::conditional<IsNotAlreadyThere,TupleCatType<TupleTensorNth,std::tuple<TensorType>>,TupleTensorNth>::type;
  using type_tensor=TupleChangeType<Test::value, ChangeTensorType, TupleTensor> ;

};





template<typename QuadratureRule, typename Tuple,typename TupleTensor,typename T>
class TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,UnaryPlus<Expression<T>>>
{
  public:
  using type=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,T>::type;
  using type_tensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,T>::type_tensor;
};

template<typename QuadratureRule, typename Tuple,typename TupleTensor,typename T>
class TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,UnaryMinus<Expression<T>>>
{
  public:
  using type=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,T>::type;
  using type_tensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,T>::type_tensor;
};

template<typename QuadratureRule, typename Tuple,typename TupleTensor,typename Left,typename Right>
class TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Addition<Expression<Left>,Expression<Right>>>
{
  public:
  using TupleNew=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Left>::type;
  using TupleNewTensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Left>::type_tensor;
  using type=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,TupleNew,TupleNewTensor,Right>::type;
  using type_tensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,TupleNew,TupleNewTensor,Right>::type_tensor;
};

template<typename QuadratureRule, typename Tuple,typename TupleTensor,typename Left,typename Right>
class TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Subtraction<Expression<Left>,Expression<Right>>>
{
  public:
  using TupleNew=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Left>::type;
  using TupleNewTensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Left>::type_tensor;
  using type=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,TupleNew,TupleNewTensor,Right>::type;
  using type_tensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,TupleNew,TupleNewTensor,Right>::type_tensor;

};

template<typename QuadratureRule, typename Tuple,typename TupleTensor,typename Left,typename Right>
class TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Multiplication<Expression<Left>,Expression<Right>>>
{
  public:
  using TupleNew=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Left>::type;
  using TupleNewTensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Left>::type_tensor;
  using type=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,TupleNew,TupleNewTensor,Right>::type;
  using type_tensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,TupleNew,TupleNewTensor,Right>::type_tensor;

};


template<typename QuadratureRule, typename Tuple,typename TupleTensor,typename Left,typename Right>
class TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Division<Expression<Left>,Expression<Right>>>
{
  public:
  using TupleNew=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Left>::type;
  using TupleNewTensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Left>::type_tensor;
  using type=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,TupleNew,TupleNewTensor,Right>::type;
  using type_tensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,TupleNew,TupleNewTensor,Right>::type_tensor;

};

template<typename QuadratureRule, typename Tuple,typename TupleTensor,typename T>
class TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,Transposed<Expression<T>>>
{
  public:
  using type=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,T>::type;
  using type_tensor=typename TupleOfCombinationFunctionsAuxAux<QuadratureRule,Tuple,TupleTensor,T>::type_tensor;
};


template<typename ...Ts>
class TupleOfCombinationFunctionsAux;




template<typename Tuple, typename TupleTensor,typename Left,typename Right,Integer QR>
class TupleOfCombinationFunctionsAux<Tuple,TupleTensor,L2DotProductIntegral<Left,Right,QR>>
{
  public:
  using L2=L2DotProductIntegral<Left,Right,QR>;
  using QRule=typename L2::QRule;
  using TupleNew=typename TupleOfCombinationFunctionsAuxAux<QRule,Tuple,TupleTensor,Left>::type;
  using TupleNewTensor=typename TupleOfCombinationFunctionsAuxAux<QRule,Tuple,TupleTensor,Left>::type_tensor;
  using type=typename TupleOfCombinationFunctionsAuxAux<QRule,TupleNew,TupleNewTensor,Right>::type;
  using type_tensor=typename TupleOfCombinationFunctionsAuxAux<QRule,TupleNew,TupleNewTensor,Right>::type_tensor;
};


template<typename Tuple, typename TupleTensor,typename Left,typename Right>
class TupleOfCombinationFunctionsAux<Tuple,TupleTensor,Addition<Expression<Left>,Expression<Right>>>
{
  public:
  using TupleNew=typename TupleOfCombinationFunctionsAux<Tuple,TupleTensor,Left>::type;
  using TupleNewTensor=typename TupleOfCombinationFunctionsAux<Tuple,TupleTensor,Left>::type_tensor;
  using type=typename TupleOfCombinationFunctionsAux<TupleNew,TupleNewTensor,Right>::type;
  using type_tensor=typename TupleOfCombinationFunctionsAux<TupleNew,TupleNewTensor,Right>::type_tensor;
};




template<typename ...Ts>
class TupleOfCombinationFunctionsForm;


template<typename Tuple,typename TupleTensor,typename Form>
class TupleOfCombinationFunctionsForm<Tuple,TupleTensor,Form>
{
  public:
  using type=typename TupleOfCombinationFunctionsAux<Tuple,TupleTensor,Form>::type;
  using type_tensor=typename TupleOfCombinationFunctionsAux<Tuple,TupleTensor,Form>::type_tensor;
};

template<typename Tuple,typename TupleTensor,typename Form,typename...Forms>
class TupleOfCombinationFunctionsForm<Tuple,TupleTensor,Form,Forms...>
{
  public:
  using TupleNew=typename TupleOfCombinationFunctionsAux<Tuple,TupleTensor,Form>::type;
  using TupleTensorNew=typename TupleOfCombinationFunctionsAux<Tuple,TupleTensor,Form>::type_tensor;
  using type=typename TupleOfCombinationFunctionsForm<TupleNew,TupleTensorNew,Forms...>::type;
  using type_tensor=typename TupleOfCombinationFunctionsForm<TupleNew,TupleTensorNew,Forms...>::type_tensor;
};



template<Integer Nmax,typename Form,typename...Forms>
class TupleOfCombinationFunctions
{
 public:
  using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
  using type=typename TupleOfCombinationFunctionsForm<emptytuple,emptytuple,Form,Forms...>::type;
  using type_tensor=typename TupleOfCombinationFunctionsForm<emptytuple,emptytuple,Form,Forms...>::type_tensor;
};


}



#endif