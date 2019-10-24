#ifndef MARS_EVALUATION_CONSTANT_TENSOR_HPP
#define MARS_EVALUATION_CONSTANT_TENSOR_HPP
#include "mars_base.hpp"
#include "mars_constant.hpp"


namespace mars {

template<typename...Parameters>
class Evaluation;

template<typename ConstType,typename...Inputs,typename QuadratureRule>
class Evaluation<Expression<ConstantTensor<ConstType,Inputs...>>,QuadratureRule>
{
 public:
 using type=ConstantTensor<ConstType,Inputs...>;
 using Input=typename type::type;
 using value_type=OperatorType<type,QuadratureRule>;


 Evaluation(){};

 Evaluation(const type& expr):
 eval_(expr)
 { };



 
 template<typename...Args, typename Jacobian,typename...OtherInputs>
 constexpr void apply(value_type& value,const Jacobian& J, const std::tuple<Args...>& tuple_of_tuple,const OtherInputs&...inputs)
 {
  // TODO here we copy the static static_value  into value, but it is useless. better to specialize evaluation
  // Assign(value,eval_.template eval<QuadratureRule::NQpoints>());
  // value_type oo(4,4,4,44,4,4);
  // decltype(eval_.template qp_eval<QuadratureRule::NQPoints>()) rrr(4,4,4,44,4,4);
  Assignment<value_type>::apply(value,eval_.template qp_eval<QuadratureRule::NQPoints>());
  std::cout<<"Evaluation<Expression<ConstantTensor "<<value<<std::endl;



 }


private:
 type eval_;
};



}
#endif

