#ifndef MARS_OPERATORS_HPP
#define MARS_OPERATORS_HPP
#include "mars_base.hpp"
#include "mars_tuple_utilities.hpp"
#include "mars_quadrature_rules.hpp"
#include "mars_expression.hpp"
#include "mars_operator_type.hpp"
#include "mars_operator_tuple_type.hpp"
#include "mars_quadrature_order.hpp"
#include "mars_matrix_static_assert.hpp"

#include "mars_operator_addition.hpp"
#include "mars_operator_apply.hpp"
#include "mars_operator_assignment.hpp"
#include "mars_operator_inner_product.hpp"
#include "mars_operator_division.hpp"
#include "mars_operator_multiplication.hpp"
#include "mars_operator_subtraction.hpp"
#include "mars_operator_trace.hpp"
#include "mars_operator_transposed.hpp"
#include "mars_operator_unary_minus.hpp"
#include "mars_operator_unary_plus.hpp"


namespace mars {








template<typename...Ts>
class CompositeOperator;

template<typename ConstType,typename...Inputs>
class ConstantTensor;

template<typename FullSpace,Integer N,typename Operator_,typename FuncType>
class Function;


////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////---------- OPERATOR IDs ---------- ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class IdentityOperator: public Expression<IdentityOperator>  
{public: IdentityOperator(){};}; 
class DivergenceOperator: public Expression<DivergenceOperator>
{public: DivergenceOperator(){};};

class GradientOperator : public Expression<GradientOperator> 
{public: GradientOperator(){}};
class CurlOperator: public Expression<CurlOperator>  
{public: CurlOperator(){}};
template<Integer N_>
class ComponentOperator: public Expression<ComponentOperator<N_>>  
 {public: static constexpr Integer N=N_;
  ComponentOperator(){}};







template<typename Left,typename Right,Integer QR>
class L2DotProductIntegral;



























// template<typename Left,typename Right>
// class InnerProduct;


// template< typename DerivedLeft,typename DerivedRight>
// class InnerProduct< Expression <DerivedLeft>, Expression <DerivedRight> > 
// : public Expression< InnerProduct< Expression <DerivedLeft>, Expression <DerivedRight> > >
// {
//   public:
//     using Left=DerivedLeft;
//     using Right=DerivedRight;

//     InnerProduct(const Expression<DerivedLeft>& left, const Expression<DerivedRight>&right)
//     : 
//     left_(left.derived()),
//     right_(right.derived())
//     {};
//     const constexpr DerivedLeft& left()const{return left_;};
//     const constexpr DerivedRight& right()const{return right_;};
//   private:
//   DerivedLeft left_;
//   DerivedRight right_;
// };

// template< typename DerivedLeft,typename DerivedRight>
// class InnerProduct< Expression <DerivedLeft>, Expression <DerivedRight> >
// Inner(const Expression<DerivedLeft>&left, const Expression<DerivedRight>&right)
// {return InnerProduct< Expression <DerivedLeft>, Expression <DerivedRight> >(left,right);}








template<typename T, typename...Ts>
class MultipleAdditionHelper;

template<typename T, typename...Ts>
using MultipleAddition=typename MultipleAdditionHelper<T,Ts...>::type;

template<typename T>
class MultipleAdditionHelper<T>
{
  public:
    using type=T;
};


template<typename T, typename...Ts>
class MultipleAdditionHelper
{
  public:
    using type=Addition<Expression<T>,Expression<typename MultipleAdditionHelper<Ts...>::type>>;
};




// template<typename MeshT, typename Left,typename Right,Integer QR, typename Form,typename...OtherTemplateArguments>
// constexpr auto Eval(const L2DotProductIntegral<MeshT,Left,Right,QR>& t,const OtherTemplateArguments&...ts)
// {return Evaluation< Expression<remove_all_t<decltype(t)>>,OtherTemplateArguments...>(t,ts...);}
template<typename Left,typename Right,Integer QR, typename Form,typename...OtherTemplateArguments>
constexpr auto Eval(const L2DotProductIntegral<Left,Right,QR>& t,const OtherTemplateArguments&...ts)
{return Evaluation< Expression<remove_all_t<decltype(t)>>,OtherTemplateArguments...>(t,ts...);}


template<typename...OtherTemplateArguments, typename T>
constexpr auto Eval(const T& t){return Evaluation< Expression<remove_all_t<decltype(t)>>,OtherTemplateArguments...>(t);}

// first parameter is an expression, while the others input are utils 
template<typename...OtherTemplateArguments,typename T,typename ...Ts>
constexpr auto Eval(const T& t,const Ts&...ts){return Evaluation< Expression<remove_all_t<decltype(t)>>,
                                                                              remove_all_t<decltype(ts)>...,OtherTemplateArguments... >(t,ts...);}

template<typename...OtherTemplateArguments,typename T,typename ...Ts>
constexpr auto Eval(const T& t, Ts&...ts){return Evaluation< Expression<remove_all_t<decltype(t)>>,
                                                                         remove_all_t<decltype(ts)>...,OtherTemplateArguments... >(t,ts...);}


template<typename...Args>
class GeneralForm;

template<typename Form,typename...Forms>
class ShapeFunctions2;

template<typename Elem>
class Jacobian;

template<typename...Ts>
class TupleOfTestTrialPairsNumbers;

template<typename TupleOfPairsNumbers, typename Form>
class TupleOfL2Products;

template<typename ...Ts>
class TupleOfEvaluationExpressionOfTypesHelper;



template<typename...Ts>
class EvaluationOfL2Inners;

template<typename Form,typename GeneralForm_, typename...GeneralForms>
class Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions2<GeneralForm_,GeneralForms...>>
{
 public:
 using type= GeneralForm<Form>;
 using FunctionSpace= typename type::FunctionSpace;
 using TupleFunctionSpace=typename type::TupleFunctionSpace;
 using ShapeFunctions=ShapeFunctions2<GeneralForm_,GeneralForms...>;
 using Shapes=typename ShapeFunctions::TupleOfTupleShapeFunction;
 using TupleOfPairsNumbers=typename GeneralForm<Form>::TupleOfPairsNumbers;
 using L2Products=typename TupleOfL2Products< TupleOfPairsNumbers, Form >::type;
 using EvaluationOfL2Inners=EvaluationOfL2Inners<Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions>>;
 const Form& operator()()const{return general_form_();};

 Evaluation(const GeneralForm<Form>& general_form,ShapeFunctions& shapesform)
 :
 general_form_(general_form)
 ,
 shapesform_(shapesform)
 ,
 eval_inners_(general_form_,shapesform)
 {}
 
 template<typename Elem>
 constexpr void apply(const Jacobian<Elem>& J)
 {
  int token_mat=1;
  eval_inners_.apply(token_mat,J);
 }

 private:
 const GeneralForm<Form>& general_form_;
 ShapeFunctions& shapesform_;
 EvaluationOfL2Inners eval_inners_;
};


template<typename Form,typename FullSpace, typename GeneralForm_, typename...GeneralForms>
constexpr auto Eval(const GeneralForm<Form>& form, ShapeFunctions2<GeneralForm_,GeneralForms...>& shapes)
{return Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions2<GeneralForm_,GeneralForms...> >(form,shapes);}











































































































template<typename...Ts>
class CompositeOperator;

template<typename MixedSpace, Integer N, typename OperatorType>
class Test;

template<typename MixedSpace, Integer N, typename OperatorType>
class Trial;

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
