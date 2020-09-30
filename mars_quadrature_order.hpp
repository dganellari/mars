#ifndef MARS_QUADRATURE_ORDER_HPP
#define MARS_QUADRATURE_ORDER_HPP

// #include "mars_constant.hpp"
#include "mars_base.hpp"
#include "mars_elem_to_sub_elem.hpp"
#include "mars_operators.hpp"

namespace mars{


template<typename...T>
class QuadratureOrder;


template<typename Elem, Integer QuadratureRuleType,Integer ActualOrder>
class CheckMaxQuadratureOrder
{
public:
  static constexpr Integer value=Min(MaxOrder<Elem,QuadratureRuleType>::value,ActualOrder);
};

template<Integer N,typename OperatorType,typename FullSpace,typename FuncType,typename...Args>
class QuadratureOrder<Function<FullSpace,N,OperatorType,FuncType>,Args...>
{ public:
  using type=Function<FullSpace,N,OperatorType,FuncType>;
  using Operator=typename type::Operator;
  using UniqueElementFunctionSpacesTupleType=typename type::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,type::value>;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;

};

template<typename ConstType,typename...Inputs,typename...Args>
class QuadratureOrder<ConstantTensor<ConstType,Inputs...>,Args...>
{ public:
  static constexpr Integer value=0;
};

template<typename MixedSpace,Integer N_, typename OperatorKind,typename...Args>
class QuadratureOrder<Test<MixedSpace,N_,OperatorKind>,Args...>
{ public:
  static constexpr Integer N=Test<MixedSpace,N_,OperatorKind>::value;
  using Test_=Test<MixedSpace,N,OperatorKind>;
  using UniqueElementFunctionSpacesTupleType=typename Test_::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename Test_::Operator;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;
};

template<typename MixedSpace,Integer N_, typename OperatorKind,typename...Args>
class QuadratureOrder<Trial<MixedSpace,N_,OperatorKind>,Args...>
{ public:
  static constexpr Integer N=Trial<MixedSpace,N_,OperatorKind>::value;
  using Trial_=Trial<MixedSpace,N,OperatorKind>;
  using UniqueElementFunctionSpacesTupleType=typename Trial_::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename Trial_::Operator;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;
};

template<template<class>class Unary,typename T, typename...Args>
class QuadratureOrder< Unary< Expression<T> >,Args... >
{ public:
  static constexpr Integer value=QuadratureOrder<T,Args...>::value;
};

// order(T) = order(+T)
template<typename T, typename...Args>
class QuadratureOrder< UnaryPlus< Expression<T> >,Args... >
{ public:
  static constexpr Integer value=QuadratureOrder<T>::value;
};

// order(T) = order(-T)
template<typename T, typename...Args>
class QuadratureOrder< UnaryMinus< Expression<T> >,Args... >
{ public:
  static constexpr Integer value=QuadratureOrder<T,Args...>::value;
};

// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right, typename...Args>
class QuadratureOrder< Addition< Expression<Left>,Expression<Right> >,Args... >
{ public:
  static constexpr Integer value=Max(QuadratureOrder<Left,Args...>::value,
                                     QuadratureOrder<Right,Args...>::value);
};

// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right, typename...Args>
class QuadratureOrder< Subtraction< Expression<Left>,Expression<Right> >,Args... >
{ public:
  static constexpr Integer value=Max(QuadratureOrder<Left,Args...>::value,
                                     QuadratureOrder<Right,Args...>::value);
};

// order(T) = order(T/REAL)
template<typename Left, typename Right, typename...Args>
class QuadratureOrder< Division<Expression<Left>, Expression<Right> >,Args... >
{ public:
  static constexpr Integer value=QuadratureOrder<Left,Args...>::value-QuadratureOrder<Right,Args...>::value;
};

// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right, typename...Args>
class QuadratureOrder< Multiplication< Expression<Left>, Expression<Right> >,Args... >
{ public:

  static constexpr Integer value=QuadratureOrder<Left,Args...>::value + QuadratureOrder<Right,Args...>::value;
};

template<typename Left, typename Right, typename...Args>
class QuadratureOrder< InnerProduct< Expression<Left>, Expression<Right> >,Args... >
{ public:

  static constexpr Integer value=QuadratureOrder<Left,Args...>::value + 
                                 QuadratureOrder<Right,Args...>::value;
};



template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<IdentityOperator,BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order;
};

template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<TraceOperator, BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order;
};

template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<TraceGradientOperator, BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order-1;
};

template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<GradientOperator, BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order-1;
};


template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<IdentityOperator, BaseFunctionSpace<RaviartThomasFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order+1;
};

template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<TraceOperator, BaseFunctionSpace<RaviartThomasFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order;
};


template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<DivergenceOperator, BaseFunctionSpace<RaviartThomasFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order;
};


// template<typename T,Integer FEFamily, Integer Order,Integer Continuity, Integer NComponents>
// class QuadratureOrder<Operator,BaseFunctionSpace<FEFamily,Order,Continuity,NComponents> >
// { public:
//   static constexpr Integer value=QuadratureOrder<IdentityOperator,BaseFunctionSpace<FEFamily,Order,Continuity,NComponents>>::value;
// };


template<typename Operator,typename Elem, Integer FEFamily, Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<Operator,ElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=QuadratureOrder<IdentityOperator,BaseFunctionSpace<FEFamily,Order,Continuity,NComponents>>::value;
};





template<template<class>class Unary, typename MixedSpace,Integer N_, typename Expr>
class QuadratureOrder<Trial<MixedSpace,N_,CompositeOperator<Expression<Unary<Expression<Expr>>>>> >
{ public:
  static constexpr Integer value=
  QuadratureOrder<Trial<MixedSpace,N_,CompositeOperator<Expression<Expr>>>>::value;
};


template<template<class,class>class Binary, typename MixedSpace,Integer N_, typename Left,typename Right>
class QuadratureOrder<Trial<MixedSpace,N_,CompositeOperator<Expression<Binary<Expression<Left>,Expression<Right>>>>> >
{ public:

  static constexpr Integer value=
  QuadratureOrder<Binary<Expression<Trial<MixedSpace,N_,CompositeOperator<Expression<Left>>>>,
                         Expression<Trial<MixedSpace,N_,CompositeOperator<Expression<Right>>>>>>::value;
};



template<template<class>class Unary, typename MixedSpace,Integer N_, typename Expr>
class QuadratureOrder<Test<MixedSpace,N_,CompositeOperator<Expression<Unary<Expression<Expr>>>>> >
{ public:
  static constexpr Integer value=
  QuadratureOrder<Test<MixedSpace,N_,CompositeOperator<Expression<Expr>>>>::value;
};


template<template<class,class>class Binary, typename MixedSpace,Integer N_, typename Left,typename Right>
class QuadratureOrder<Test<MixedSpace,N_,CompositeOperator<Expression<Binary<Expression<Left>,Expression<Right>>>>> >
{ public:

  static constexpr Integer value=
  QuadratureOrder<Binary<Expression<Test<MixedSpace,N_,CompositeOperator<Expression<Left>>>>,
                         Expression<Test<MixedSpace,N_,CompositeOperator<Expression<Right>>>>>>::value;
};



template<typename MixedSpace,Integer N_, typename Expr>
class QuadratureOrder<Trial<MixedSpace,N_,CompositeOperator<Expression<Expr>>> >
{ public:
  static constexpr Integer N=Trial<MixedSpace,N_,Expression<Expr>>::value;
  using Trial_=Trial<MixedSpace,N,Expression<Expr>>;
  using UniqueElementFunctionSpacesTupleType=typename Trial_::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename Trial_::Operator;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Expr,BaseFunctionSpace>::value;
};


template<typename MixedSpace,Integer N_, typename Expr>
class QuadratureOrder<Test<MixedSpace,N_,CompositeOperator<Expression<Expr>>> >
{ public:
  static constexpr Integer N=Test<MixedSpace,N_,Expression<Expr>>::value;
  using Test_=Test<MixedSpace,N,Expression<Expr>>;
  using UniqueElementFunctionSpacesTupleType=typename Test_::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename Test_::Operator;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Expr,BaseFunctionSpace>::value;
};



// /////////////////////////////////////////////////////////////////////////////
// //// VolumeOrSurfaceElem return the volume/surface Element if true/false ////
// /////////////////////////////////////////////////////////////////////////////

// template<typename Elem, bool=true>
// class VolumeOrSurfaceElem;

// template<template<Integer Dim,Integer ManifoldDim>class Elem_,Integer Dim, Integer ManifoldDim>
// class VolumeOrSurfaceElem<Elem_<Dim,ManifoldDim>,true>
// {
//  public:
//   using type=Elem_<Dim,ManifoldDim>; 
// };

// template<template<Integer Dim,Integer ManifoldDim>class Elem_,Integer Dim, Integer ManifoldDim>
// class VolumeOrSurfaceElem<Elem_<Dim,ManifoldDim>,false>
// {
//  public:
//   using type=Elem_<Dim,ManifoldDim-1>; 
// };



}


#endif