#ifndef MARS_EVALUATION_GENERAL_FORM_HPP
#define MARS_EVALUATION_GENERAL_FORM_HPP
#include "mars_base.hpp"
#include "mars_general_form.hpp"
#include "mars_operators.hpp"


namespace mars {


template<typename...Ts>
class EvaluationOfL2Inners;

template<typename Form,typename GeneralForm_, typename...GeneralForms>
class Evaluation<Expression<GeneralForm<Form>>,ShapeFunctionsCollection<GeneralForm_,GeneralForms...>>
{
 public:
 using type= GeneralForm<Form>;
 using FunctionSpace= typename type::FunctionSpace;
 using TupleFunctionSpace=typename type::TupleFunctionSpace;
 using ShapeFunctions=ShapeFunctionsCollection<GeneralForm_,GeneralForms...>;
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
constexpr auto Eval(const GeneralForm<Form>& form, ShapeFunctionsCollection<GeneralForm_,GeneralForms...>& shapes)
{return Evaluation<Expression<GeneralForm<Form>>,ShapeFunctionsCollection<GeneralForm_,GeneralForms...> >(form,shapes);}




}
#endif

