#ifndef MARS_EVALUATION_UNARY_HPP
#define MARS_EVALUATION_UNARY_HPP
#include "mars_base.hpp"
#include "mars_operators.hpp"


namespace mars {

template<typename...Parameters>
class Evaluation;

template<template<class>class Unary,typename Derived,typename...OtherTemplateArguments>
class Evaluation<Expression<Unary<Expression<Derived> > >,OtherTemplateArguments...>
{
 public:
 using type=Unary<Expression<Derived>>;
 using Eval=Evaluation<Expression<Derived>,OtherTemplateArguments...>;
 using subtype= OperatorType<Derived,OtherTemplateArguments...>;
 using outputsubtype= OperatorType<type,OtherTemplateArguments...>;

 Evaluation(){};
 
 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_(Eval(expr_()))
 {};
 
 template<typename...Inputs>
 void apply(outputsubtype& output,const Inputs&...inputs)
 {

  std::cout<<"pre eval Evaluation<Expression<UnaryOperator "<<std::endl;
  eval_.apply(derived_value_,inputs...);
  OperatorApply<Unary<subtype>>::apply(output,derived_value_);
  std::cout<<"after Unary Evaluation<Expression<UnaryOperator "<<std::endl;
  std::cout<<"derived_value_="<<derived_value_<<std::endl;
  std::cout<<"output="<<output<<std::endl;
 }
       auto  expression()     {return expr_;}
 const auto& expression()const{return expr_;}

private:

 type expr_;
 Eval eval_;
 subtype derived_value_;
};




}
#endif

