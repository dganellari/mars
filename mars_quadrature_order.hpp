#ifndef MARS_QUADRATURE_ORDER_HPP
#define MARS_QUADRATURE_ORDER_HPP

// #include "mars_constant.hpp"
#include "mars_base.hpp"


namespace mars{

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

template<typename T>
class Transposed;


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
  using Test=Test<MixedSpace,N,OperatorKind>;
  using UniqueElementFunctionSpacesTupleType=typename Test::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename Test::Operator;
  // using BaseFunctionSpace=GetType<BaseFunctionSpaceAndElement,1>;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;
};

template<typename MixedSpace,Integer N_, typename OperatorKind,typename...Args>
class QuadratureOrder<Trial<MixedSpace,N_,OperatorKind>,Args...>
{ public:
  static constexpr Integer N=Trial<MixedSpace,N_,OperatorKind>::value;
  using Trial=Trial<MixedSpace,N,OperatorKind>;
  using UniqueElementFunctionSpacesTupleType=typename Trial::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename Trial::Operator;
  // using BaseFunctionSpace=GetType<BaseFunctionSpaceAndElement,1>;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;
};


// order(T) = order(+T)
template<typename T>
class QuadratureOrder< UnaryPlus< Expression<T> > >
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
class QuadratureOrder< Contraction2< Expression<Left>, Expression<Right> >,Args... >
{ public:

  static constexpr Integer value=QuadratureOrder<Left,Args...>::value + 
                                 QuadratureOrder<Right,Args...>::value;
};


template<typename T, typename...Args>
class QuadratureOrder< Transposed< Expression<T> >,Args... >
{ public:
  static constexpr Integer value=QuadratureOrder<T,Args...>::value;
};


template<Integer Order,Integer Continuity, Integer NComponents>
class QuadratureOrder<IdentityOperator,BaseFunctionSpace<LagrangeFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order;
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
class QuadratureOrder<DivergenceOperator, BaseFunctionSpace<RaviartThomasFE,Order,Continuity,NComponents> >
{ public:
  static constexpr Integer value=Order;
};


template<typename MixedSpace,Integer N_, typename Expr>
class QuadratureOrder<Trial<MixedSpace,N_,CompositeOperator<Expression<Expr>>> >
{ public:
  static constexpr Integer N=Trial<MixedSpace,N_,Expression<Expr>>::value;
  using Trial=Trial<MixedSpace,N,Expression<Expr>>;
  using UniqueElementFunctionSpacesTupleType=typename Trial::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename Trial::Operator;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Expr,BaseFunctionSpace>::value;
};

template<typename MixedSpace,Integer N_, typename Expr>
class QuadratureOrder<Test<MixedSpace,N_,CompositeOperator<Expression<Expr>>> >
{ public:
  static constexpr Integer N=Test<MixedSpace,N_,Expression<Expr>>::value;
  using Test=Test<MixedSpace,N,Expression<Expr>>;
  using UniqueElementFunctionSpacesTupleType=typename Test::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename Test::Operator;
  using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Expr,BaseFunctionSpace>::value;
};


// template<template<class,Integer,class> class TestOrTrial_,typename MixedSpace,Integer N_, typename Expr>
// class QuadratureOrder<TestOrTrial_<MixedSpace,N_,CompositeOperator<Expression<Expr>>> >
// { public:
//   static constexpr Integer N=TestOrTrial_<MixedSpace,N,Expression<Expr>>::value;
//   using TestOrTrial=TestOrTrial_<MixedSpace,N,Expression<Expr>>;
//   using UniqueElementFunctionSpacesTupleType=typename TestOrTrial::UniqueElementFunctionSpacesTupleType;
//   using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
//   using Operator=typename TestOrTrial::Operator;
//   using BaseFunctionSpace=Elem2FunctionSpace<BaseFunctionSpaceAndElement>;
//   static constexpr Integer value=QuadratureOrder<Expr,BaseFunctionSpace>::value;
// };

// template<template<class,Integer,class> class TestOrTrial_,typename MixedSpace,Integer N_, 
//          typename FullSpace,Integer M,typename Operator_,typename FuncType>
// class QuadratureOrder<TestOrTrial_<MixedSpace,N_,CompositeOperator<Expression<Function<FullSpace,M,Operator_,FuncType>>>> >
// { public:
//   static constexpr Integer value=QuadratureOrder<Function<FullSpace,M,Operator_,FuncType>>::value;
// };





}


#endif