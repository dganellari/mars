#ifndef MARS_EVALUATION_TEST_OR_TRIAL_HPP
#define MARS_EVALUATION_TEST_OR_TRIAL_HPP
#include "mars_base.hpp"
#include "mars_test_function.hpp"
#include "mars_trial_function.hpp"


namespace mars {




template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, Integer N, typename Operator_,typename...OtherTemplateArguments>
class Evaluation<Expression<TestOrTrial<MixedSpace,N,Operator_>>,OtherTemplateArguments...>
{
 public:
 using type= TestOrTrial<MixedSpace,N,Operator_>;

 static_assert((IsSame<TestOrTrial<MixedSpace,N,Operator_>,Test<MixedSpace,N,Operator_>>::value ||
                IsSame<TestOrTrial<MixedSpace,N,Operator_>,Trial<MixedSpace,N,Operator_>>::value )
               && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");
 using FunctionSpace=typename type::FunctionSpace;
 using FunctionSpaces=GetType<typename type::UniqueElementFunctionSpacesTupleType,type::value>;
 using Elem=typename FunctionSpaces::Elem;
 using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
 using Operator=typename type::Operator;
 using value_type=OperatorType<type,OtherTemplateArguments...>;

 Evaluation(){};

 Evaluation(const type& expr):
 eval_(expr)
 {};
 
 //  template<typename...Forms, typename FiniteElem, typename...Inputs>
 // constexpr void apply(value_type& value,const FiniteElem& J, const ShapeFunctions2<Forms...>& shape_functions)
 template<typename...Args, typename FiniteElem,typename...Inputs>
 constexpr void apply(value_type& value,const FiniteElem& J, const std::tuple<Args...>& shape_functions,const Inputs&...inputs)

 {
  using TupleOfTupleShapeFunction=std::tuple<Args...>;
  using tuple_type=GetType<TupleOfTupleShapeFunction,type::value>;
  const auto& tuple=tuple_get<type::value>(shape_functions);
    // const auto& tuple=tuple_get<type::value>(shape_functions());
  constexpr Integer M=TypeToTupleElementPosition<ShapeFunction<Elem,BaseFunctionSpace,Operator,OtherTemplateArguments...>,tuple_type>::value;
    // std::cout<<"befor Evaluation<Expression<TestOrTrial<MixedSpace,N,Operator_> "<<value<<std::endl;

 // Number<M> e5e(5);

 //  tuple_type ok(1);
 //  decltype(tuple) ee(3);
  Assign(value,tuple_get<M>(tuple).eval());




  // Assignment<value_type>::apply(value,tuple_get<M>(tuple).eval());
  std::cout<<"after Evaluation<Expression<TestOrTrial<MixedSpace,N,Operator_> "<<value<<std::endl;

 }

private:
 
 type eval_;
};









template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, Integer N, typename Expr,typename...OtherTemplateArguments>
class Evaluation<Expression<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>>,OtherTemplateArguments...>
{
 public:
 using Operator=CompositeOperator<Expression<Expr>>;
 using type= TestOrTrial<MixedSpace,N,Operator>;

 static_assert((IsSame<type,Test<MixedSpace,N,Operator>>::value ||
                IsSame<type,Trial<MixedSpace,N,Operator>>::value )
               && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,Composite<Operator>>>>,TestOrTrial=Test or Trial ");
 using FunctionSpace=typename type::FunctionSpace;
 using FunctionSpaces=GetType<typename type::UniqueElementFunctionSpacesTupleType,type::value>;
 using Elem=typename FunctionSpaces::Elem;
 using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
 // using value_type=OperatorType<type,OtherTemplateArguments...>;

 Evaluation(){};

 Evaluation(const type& expr):
 eval_(expr)
 {};
 
 //  template<typename ValueType,typename...Forms, typename FiniteElem, typename...Inputs>
 // constexpr void apply(ValueType& value,const FiniteElem& J, const ShapeFunctions2<Forms...>& shape_functions)
  template<typename ValueType, typename FiniteElem, typename...Args1,typename...Args2,typename...Args3>
 constexpr void apply(ValueType& value,const FiniteElem& J, const std::tuple<Args1...>& tuple_shape_functions,const std::tuple<Args2...>&tuple_tensor,const std::tuple<Args3...>&tuple_composite_shapes)

 {
 using single_type=Evaluation<Expression<typename FormOfCompositeOperatorType<type>::type >,OtherTemplateArguments...>;

 using TupleOfTupleCompositeShapeFunctionEval=std::tuple<Args3...>;
 using tuple_type=GetType<TupleOfTupleCompositeShapeFunctionEval,type::value>;
 // using single_type=Evaluation<Expression<typename FormOfCompositeOperatorType<type>::type >,OtherTemplateArguments...>;
 // const auto& tuple=tuple_get<type::value>(tuple_tensor);
 constexpr Integer M=TypeToTupleElementPosition<single_type,tuple_type>::value;

 auto& single_tensor=tuple_get<type::value,M >(tuple_tensor);
 Assign(value,single_tensor);


 // Assignment<ValueType>::apply(value,single_tensor);
   // std::cout<<"single_tensor= "<<single_tensor<<std::endl;

  // single_type ok1(1);
  // TupleOfTupleCompositeShapeFunctionEval ok2(3);
  std::cout<<"Evaluation<Expression<TestOrTrial<MixedSpace,N,Operator_> "<<value<<std::endl;


  }

private:
 
 type eval_;
};


}
#endif