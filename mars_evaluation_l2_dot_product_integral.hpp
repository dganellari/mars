
#ifndef MARS_EVALUATION_L2_DOT_PRODUCT_INTEGRAL_HPP
#define MARS_EVALUATION_L2_DOT_PRODUCT_INTEGRAL_HPP
#include "mars_l2_dot_product_integral.hpp"
#include "mars_evaluation_general_form_utils.hpp"
namespace mars {






template<typename Left_,typename Right_,bool VolumeIntegral, Integer QR, typename...Forms>
class Evaluation<Expression<L2DotProductIntegral<Left_,Right_,VolumeIntegral,QR>>, ShapeFunctionsCollection<Forms...>>
{
 public:
 using Left=Left_;
 using Right=Right_;
 using type= L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
 using QRule=typename type ::QRule;
 using TestTrialSpaces=typename type::TestTrialSpaces;
 using subtype= OperatorType<type,GetType<typename type::form>>;
 static constexpr bool PositiveWeights= IsPositive(QRule::qp_weights);
 
 Evaluation(const type& expr, ShapeFunctionsCollection<Forms...>& shape_functions):
 expr_(expr.left(),expr.right())
 ,
 shape_functions_(shape_functions)
 ,
 local_tensor_(expr)
 ,
 label_(expr.label())
 {};
 

 template<bool VolumeIntegralAux,typename Elem, typename...DofMaps>
 std::enable_if_t<(VolumeIntegralAux==true),void>
  apply_aux(subtype& mat, const FiniteElem<Elem>& J, const DofMaps&...dofmaps)
 {

  // changed todo fixme
   // local_tensor_.apply(mat,J,shape_functions_());
  std::cout<<"pre Evaluation<Expression<L2DotProductIntegral local tensor="<<std::endl;
  

  local_tensor_.apply(mat,J,shape_functions_, dofmaps...);//(),shape_functions_.composite_tensor(),shape_functions_.composite_shapes());
  std::cout<<"after Evaluation<Expression<L2DotProductIntegral local tensor="<<std::endl;
  std::cout<<mat<<std::endl;
 }

 template<bool VolumeIntegralAux,typename Elem, typename...DofMaps>
 std::enable_if_t<(VolumeIntegralAux==false),void>
  apply_aux(subtype& mat, const FiniteElem<Elem>& J, const DofMaps&...dofmaps)
 {

  // changed todo fixme
   // local_tensor_.apply(mat,J,shape_functions_());
  std::cout<<"L2DotProductIntegral LABEL="<<label_<<std::endl;
  std::cout<<"side="<<J.side_id()<<std::endl;
  std::cout<<"tag="<<J.side_tag()<<std::endl;
  
  if(J.side_tag()==label_)
   { 
    std::cout<<"VALUTO IL BOUNDARY"<<std::endl;
    local_tensor_.apply(mat,J,shape_functions_, dofmaps...);//(),shape_functions_.composite_tensor(),shape_functions_.composite_shapes());
    }
  std::cout<<"after BOUNDAR Evaluation<Expression<L2DotProductIntegral local tensor="<<std::endl;
  std::cout<<mat<<std::endl;
 }


 template<typename...Inputs>
 void apply(subtype& mat, const Inputs&...inputs)
 {apply_aux<VolumeIntegral>(mat,inputs...);}

 auto operator()(){return local_tensor_;}

       auto  expression(){return expr_;}
 const auto& expression()const{return expr_;}



private:
 type expr_;
 ShapeFunctionsCollection<Forms...>& shape_functions_;
 LocalTensor<PositiveWeights,TestTrialSpaces,type,GetType<typename type::form>> local_tensor_;
 Integer label_;
};


template<typename Left,typename Right,bool VolumeIntegral,Integer QR, typename Form,typename...OtherTemplateArguments>
constexpr auto Eval(const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& t,const OtherTemplateArguments&...ts)
{return Evaluation< Expression<remove_all_t<decltype(t)>>,OtherTemplateArguments...>(t,ts...);}



}
#endif
