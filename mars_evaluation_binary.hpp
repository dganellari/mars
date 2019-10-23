#ifndef MARS_EVALUATION_BINARY_HPP
#define MARS_EVALUATION_BINARY_HPP
#include "mars_base.hpp"
#include "mars_operators.hpp"


namespace mars {

template<typename...Parameters>
class Evaluation;


template<template<class,class>class Binary, typename DerivedRight,typename DerivedLeft,typename...OtherTemplateArguments>
class Evaluation<Expression<Binary< Expression<DerivedLeft>  ,  
                                      Expression<DerivedRight>>>,
                 OtherTemplateArguments...>
{
 public:
 using type=Binary<Expression<DerivedLeft>,Expression<DerivedRight>>;
 using EvalLeft=Evaluation<Expression<DerivedLeft>,OtherTemplateArguments...>;
 using EvalRight=Evaluation<Expression<DerivedRight>,OtherTemplateArguments...>;
 using subtypeleft=OperatorType<DerivedLeft,OtherTemplateArguments...>;
 using subtyperight=OperatorType<DerivedRight,OtherTemplateArguments...>;

 Evaluation(){};
 

 Evaluation(const Expression<type>& expr,OtherTemplateArguments&...args):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left(),args...)),
 eval_right_(EvalRight(expr_.right(),args...))
 {};

 Evaluation(const Expression<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(EvalRight(expr_.right()))
 {};

 template<typename subtype,typename...OtherTemplateArguments2,typename...Inputs>
 void apply(subtype& output,const Inputs&...inputs)
 {
  
  eval_left_.apply(left_value_,inputs...);
  eval_right_.apply(right_value_,inputs...);
  OperatorApply<Binary<subtypeleft,subtyperight>>::apply(output,left_value_,right_value_);

  std::cout<<"left_value= "<<left_value_<<std::endl;
  std::cout<<"right_value_= "<<right_value_<<std::endl;
  std::cout<<"Evaluation<Expression<Binary "<<output<<std::endl;

 }


       auto  expression()     {return expr_;}
 const auto& expression()const{return expr_;}

private:

 type expr_;
 subtypeleft left_value_;
 subtyperight right_value_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
};


}
#endif

