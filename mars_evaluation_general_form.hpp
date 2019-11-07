#ifndef MARS_EVALUATION_GENERAL_FORM_HPP
#define MARS_EVALUATION_GENERAL_FORM_HPP
#include "mars_base.hpp"
#include "mars_general_form.hpp"
#include "mars_operators.hpp"
#include "mars_evaluation_general_form_utils.hpp"

namespace mars {


// template<typename...Ts>
// class EvaluationOfL2Inners;

// template<typename Form,typename GeneralForm_, typename...GeneralForms>
// class Evaluation<Expression<GeneralForm<Form>>,ShapeFunctionsCollection<GeneralForm_,GeneralForms...>>
// {
//  public:
//  using type= GeneralForm<Form>;
//  using FunctionSpace= typename type::FunctionSpace;
//  using TupleFunctionSpace=typename type::TupleFunctionSpace;
//  using ShapeFunctions=ShapeFunctionsCollection<GeneralForm_,GeneralForms...>;
//  // using Shapes=typename ShapeFunctions::TupleOfTupleShapeFunctionVolumetric;
//  using TupleOfPairsNumbers=typename GeneralForm<Form>::TupleOfPairsNumbers;
//  using L2Products=typename TupleOfL2Products2<0, TupleOfPairsNumbers, Form >::type;
//  using EvaluationOfL2Inners=EvaluationOfL2Inners<Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions>>;
//  const Form& operator()()const{return general_form_();};

//  Evaluation(const GeneralForm<Form>& general_form,ShapeFunctions& shapesform)
//  :
//  general_form_(general_form)
//  ,
//  shapesform_(shapesform)
//  ,
//  eval_inners_(general_form_,shapesform)
//  {}
 
//  template<typename Elem>
//  constexpr void apply(const FiniteElem<Elem>& J)
//  {
//   int token_mat=1;
//   eval_inners_.apply(token_mat,J);
//  }

//  private:
//  const GeneralForm<Form>& general_form_;
//  ShapeFunctions& shapesform_;
//  EvaluationOfL2Inners eval_inners_;
// };


template<Integer H,typename...Ts>
class EvaluationOfL2Inners;


template<typename ShapeFunctions_>
class Evaluation<Expression<EmptyExpression>,ShapeFunctions_>
{
public:
 template<typename...Ts>
 Evaluation(){}
 template<typename Elem>
 constexpr void apply(const FiniteElem<Elem>& J)
 {}
};

template<typename Form,typename GeneralForm_, typename...GeneralForms>
class Evaluation<Expression<GeneralForm<Form>>,ShapeFunctionsCollection<GeneralForm_,GeneralForms...>>
{
 public:
 using type= GeneralForm<Form>;
 using FunctionSpace= typename type::FunctionSpace;
 using TupleFunctionSpace=typename type::TupleFunctionSpace;
 using ShapeFunctions=ShapeFunctionsCollection<GeneralForm_,GeneralForms...>;
 // using Shapes=typename ShapeFunctions::TupleOfTupleShapeFunctionVolumetric;
 using TupleOfPairsNumbers=typename GeneralForm<Form>::TupleOfPairsNumbers;
 template<Integer H>
 using L2Products=typename TupleOfL2Products2<H, TupleOfPairsNumbers, Form >::type;
 // using L2ProductsSurface=typename TupleOfL2Products2<1, TupleOfPairsNumbers, Form >::type;
 using EvaluationOfL2InnersVolume=EvaluationOfL2Inners<0,Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions>>;
 using EvaluationOfL2InnersSurface=EvaluationOfL2Inners<1,Evaluation<Expression<GeneralForm<Form>>,ShapeFunctions>>;

 const Form& operator()()const{return general_form_();};

 Evaluation(const GeneralForm<Form>& general_form,ShapeFunctions& shapesform)
 :
 general_form_(general_form)
 ,
 shapesform_(shapesform)
 ,
 eval_inners_volume_(general_form_,shapesform)
 ,
 eval_inners_surface_(general_form_,shapesform)
 {}
 
 template<typename Output,typename Elem>
 constexpr void apply(Output& token_mat,const FiniteElem<Elem>& J)
 {
  eval_inners_volume_.apply(token_mat,J,general_form_.spaces_ptr()->dofmap2());
 }

 template<typename Output, typename Elem>
 constexpr void apply_boundary(Output& token_mat,const FiniteElem<Elem>& FE )
 {
  eval_inners_surface_.apply_boundary(token_mat,FE,general_form_.spaces_ptr()->dofmap2());
 }
 private:
 const GeneralForm<Form>& general_form_;
 ShapeFunctions& shapesform_;
 EvaluationOfL2InnersVolume eval_inners_volume_;
 EvaluationOfL2InnersSurface eval_inners_surface_;
};



template<typename Form,typename FullSpace, typename GeneralForm_, typename...GeneralForms>
constexpr auto Eval(const GeneralForm<Form>& form, ShapeFunctionsCollection<GeneralForm_,GeneralForms...>& shapes)
{return Evaluation<Expression<GeneralForm<Form>>,ShapeFunctionsCollection<GeneralForm_,GeneralForms...> >(form,shapes);}




}
#endif

