
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
 expr_(expr.left(),expr.right(),expr.label())
 ,
 shape_functions_(shape_functions)
 ,
 local_tensor_(expr)
 ,
 label_(expr.label())
 {};
 

 template<bool VolumeIntegralAux,typename Elem, typename...DofMaps>
 std::enable_if_t<(VolumeIntegralAux==true),void>
  apply_aux(subtype& mat, FiniteElem<Elem>& J, const DofMaps&...dofmaps)
 {

  // changed todo fixme
   // local_tensor_.apply(mat,J,shape_functions_());
  // // std::cout<<"pre Evaluation<Expression<L2DotProductIntegral local tensor="<<std::endl;
  
    std::cout<<"Evaluation<Expression<L2DotProductIntegral--------------VOLUME"<<std::endl;

  local_tensor_.apply(mat,J,shape_functions_, dofmaps...);//(),shape_functions_.composite_tensor(),shape_functions_.composite_shapes());
  // // std::cout<<"after Evaluation<Expression<L2DotProductIntegral local tensor="<<std::endl;
  // // std::cout<<mat<<std::endl;
 }

 template<bool VolumeIntegralAux,typename Elem, typename...DofMaps>
 std::enable_if_t<(VolumeIntegralAux==false),void>
  apply_aux(subtype& mat, FiniteElem<Elem>& J, const DofMaps&...dofmaps)
 {

  // changed todo fixme
   // local_tensor_.apply(mat,J,shape_functions_());
  std::cout<<"L2DotProductIntegral LABEL="<<label_<<std::endl;
  std::cout<<"side="<<J.side_id()<<std::endl;
  std::cout<<"tag="<<J.side_tag()<<std::endl;
  std::cout<<"label_="<<label_<<std::endl;
  std::cout<<mat<<std::endl;
  std::cout<<"not evaluated yet"<<mat<<std::endl;
  // // 

  if(J.side_tag()==label_)
   { 
    std::cout<<"Evaluation<Expression<L2DotProductIntegral-----------------VALUTO IL BOUNDARY "<<label_<<std::endl;
    local_tensor_.apply(mat,J,shape_functions_, dofmaps...);//(),shape_functions_.composite_tensor(),shape_functions_.composite_shapes());
    std::cout<<"after BOUNDAR Evaluation<Expression<L2DotProductIntegral local tensor="<<std::endl;
    std::cout<<mat<<std::endl;
    }
    else
    {

      for(Integer i=0;i<mat.rows();i++)
        for(Integer j=0;j<mat.cols();j++)
          mat(i,j)=0.0;
    }

 }


 template<typename...Inputs>
 void apply(subtype& mat, Inputs&...inputs)
 {apply_aux<VolumeIntegral>(mat,inputs...);
 // // std::cout<<"after apply Evaluation<Expression<L2DotProductIntegral="<<std::endl;
 }

 auto operator()(){return local_tensor_;}

       auto  expression(){return expr_;}
 const auto& expression()const{return expr_;}



private:
 type expr_;
 ShapeFunctionsCollection<Forms...>& shape_functions_;
 // LocalTensor<PositiveWeights,TestTrialSpaces,type,GetType<typename type::form>> local_tensor_;

 LocalTensor<true,TestTrialSpaces,type,GetType<typename type::form>> local_tensor_;
 Integer label_;
};


template<typename Left,typename Right,bool VolumeIntegral,Integer QR, typename Form,typename...OtherTemplateArguments>
constexpr auto Eval(const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& t,const OtherTemplateArguments&...ts)
{  
  // // std::cout<<" <<<<<<<<<eval l2 dot product>>>>>>>>>>>> "<<std::endl;
  return Evaluation< Expression<remove_all_t<decltype(t)>>,OtherTemplateArguments...>(t,ts...);}



}
#endif
