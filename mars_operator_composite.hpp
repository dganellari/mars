#ifndef MARS_OPERATOR_COMPOSITE_HPP
#define MARS_OPERATOR_COMPOSITE_HPP
#include "mars_base.hpp"
#include "mars_matrix.hpp"
#include "mars_test_function.hpp"
#include "mars_trial_function.hpp"

namespace mars {


template<typename Expr>
class CompositeOperator<Expression<Expr>>
{
public:
  CompositeOperator(const Expr& expr):expr_(expr){}


template<typename MixedSpace,Integer N> 
constexpr Test<MixedSpace,N,CompositeOperator> 
operator()(const Test<MixedSpace,N,IdentityOperator>& t)
{return Test<MixedSpace,N,CompositeOperator> (t.spaces_ptr(),expr_);}

template<typename MixedSpace,Integer N> 
constexpr Trial<MixedSpace,N,CompositeOperator> 
operator()(const Trial<MixedSpace,N,IdentityOperator>& t)
{return Trial<MixedSpace,N,CompositeOperator> (t.spaces_ptr(),expr_);}



constexpr auto& composite_operator()const{return expr_;}
private:
  Expr expr_;
};


template<typename Expr>
constexpr auto NewOperator(const Expr& expr){return CompositeOperator<Expression<Expr>>(expr);}


}
#endif

