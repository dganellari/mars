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
 using DofsDM= typename FunctionSpace::DofsDM;
 using ElemDofMap= typename DofsDM::ElemDofMap;
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

 // Evaluation(const Form& form,ShapeFunctions& shapesform)
 // :
 // general_form_(form)
 // ,
 // shapesform_(shapesform)
 // ,
 // eval_inners_volume_(general_form_,shapesform)
 // ,
 // eval_inners_surface_(general_form_,shapesform)
 // {}


 template<typename Output,typename Elem>
 constexpr void apply(Output& token_mat, FiniteElem<Elem>& FE)
 {
    // // std::cout<<" eval general form volume begin"<<std::endl;
  // eval_inners_volume_.apply(token_mat,FE,general_form_.spaces_ptr()->dofmap2());
eval_inners_volume_.apply(token_mat,FE,general_form_.spaces_ptr()->dofsdofmap());
    
// auto dofsdofmap=general_form_.spaces_ptr()->dofsdofmap();
// // DofsDM::DofMap feef(5);
// // DofsDM fe(5);
// // typename DofsDM::ElemDofMap fefe(5,4,5,6,7);
// // Form eeee(5,6);
// // decltype(general_form_.spaces_ptr()) ddeede(5);
// // decltype(dofsdofmap) fee(5); 
// // FunctionSpace mm(5);
// // for(Integer i=0;i<dofsdofmap.space_dofs().size();i++)
// // {
// //     // std::cout<<"===== space="<<i<<std::endl;
// //     for(Integer j=0;j<dofsdofmap.space_dofs_size(i);j++)
// //     {
// //       // for(Integer k=0;k<sd3[i]->operator[](j).size();k++)
// //          // std::cout<<dofsdofmap.space_dofs_get(i,j)<<std::endl;
// //     }

// // }

//     // std::cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||-----dofmap of space="<<0<<std::endl;

// for(Integer i=0;i<dofsdofmap.template dofmap_size<0>();i++)
// {
//          // std::cout<<dofsdofmap.n_dofs()<<std::endl;
//          // std::cout<<dofsdofmap.template dofmap_get<0>(i);
//          // std::cout<<std::endl;
// }
//     // std::cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||-----dofmap of space="<<1<<std::endl;

// for(Integer i=0;i<dofsdofmap.template dofmap_size<1>();i++)
// {
//          // std::cout<<dofsdofmap. template dofmap_get<1>(i)<<" ";
//          // std::cout<<std::endl;
// }

//     // std::cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||-----dofmap of space="<<2<<std::endl;

// for(Integer i=0;i<dofsdofmap.template dofmap_size<2>();i++)
// {
//          // std::cout<<dofsdofmap.template dofmap_get<2>(i)<<" ";
//          // std::cout<<std::endl;
// }
//     // std::cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||-----dofmap of space="<<3<<std::endl;

// for(Integer i=0;i<dofsdofmap.template dofmap_size<3>();i++)
// {
//          // std::cout<<dofsdofmap.template dofmap_get<3>(i)<<" ";
//          // std::cout<<std::endl;
// }


// // std::cout<<"-----dofmap of space="<<dofsdofmap.n_dofs()<<std::endl;
// // std::cout<<"-----dofmap of space="<<dofsdofmap.cumulative_n_dofs()<<std::endl;

  // // std::cout<<" eval general form volume end"<<std::endl;
 }

 template<typename Output, typename Elem>
 constexpr void apply_boundary(Output& token_mat, FiniteElem<Elem>& FE )
 {
  // eval_inners_surface_.apply_boundary(token_mat,FE,general_form_.spaces_ptr()->dofmap2());
  eval_inners_surface_.apply_boundary(token_mat,FE,general_form_.spaces_ptr()->dofsdofmap());
 }
 private:
 const GeneralForm<Form>& general_form_;
 ShapeFunctions& shapesform_;
 EvaluationOfL2InnersVolume eval_inners_volume_;
 EvaluationOfL2InnersSurface eval_inners_surface_;
};



// template<typename Form,typename FullSpace, typename GeneralForm_, typename...GeneralForms>
// constexpr auto Eval(const GeneralForm<Form>& form, ShapeFunctionsCollection<GeneralForm_,GeneralForms...>& shapes)
// {return Evaluation<Expression<GeneralForm<Form>>,ShapeFunctionsCollection<GeneralForm_,GeneralForms...> >(form,shapes);}



template<typename Form, typename GeneralForm_, typename...GeneralForms>
constexpr auto Eval(const GeneralForm<Form>& form, ShapeFunctionsCollection<GeneralForm_,GeneralForms...>& shapes)
{

    return Evaluation<Expression<GeneralForm<Form>>,ShapeFunctionsCollection<GeneralForm_,GeneralForms...> >(form,shapes);}
}


#endif

