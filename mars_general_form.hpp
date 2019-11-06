
#ifndef MARS_GENERAL_FORM_HPP
#define MARS_GENERAL_FORM_HPP
#include "mars_base.hpp"
#include "mars_l2_dot_product_integral.hpp"

namespace mars {


template<typename...Args>
class GeneralForm;

template<typename Form_>
class GeneralForm<Form_>
{

  public:
  using Form=Form_;
  template<typename T,Integer N>
  class KindType;

  template<typename T>
  class KindType<T,0>
  {
  public:
    using type=typename T::TupleFunctionSpace;
  };


 template<typename T>
  class KindType<T,1>
  {
  public:
    using type=typename T::UniqueElementFunctionSpacesTupleType;
  };


 template<typename T>
  class KindType<T,2>
  {
  public:
    using type=typename T::form;
  };

  template<typename FormTmp, Integer Kind>
  class ClassHelper;

  template<typename Left,typename Right, bool VolumeIntegral, Integer QR, Integer Kind>
  class ClassHelper<L2DotProductIntegral<Left,Right,VolumeIntegral,QR>,Kind  >
  {
  public:
   using type=RemoveTupleDuplicates< typename KindType<L2DotProductIntegral<Left,Right,VolumeIntegral,QR>,Kind>::type>;
  }; 

  template<typename Left1,typename Right1, typename Left2,typename Right2, bool VolumeIntegral1,bool VolumeIntegral2,Integer QR1, Integer QR2, Integer Kind>
  class ClassHelper<Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,
                             Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>> >,Kind  >
  {
  public:
   using type=RemoveTupleDuplicates< TupleCatType< typename KindType<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>,Kind>::type,
                                                   typename KindType<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>,Kind>::type >>;
  }; 


  template<typename Left1,typename Right1, typename Left2,typename Right2, bool VolumeIntegral1,bool VolumeIntegral2,Integer QR1, Integer QR2, Integer Kind>
  class ClassHelper<Subtraction<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,
                                Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>> >,Kind  >
  {
  public:
   using type=RemoveTupleDuplicates< TupleCatType< typename KindType<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>,Kind>::type,
                                                   typename KindType<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>,Kind>::type >>;
  }; 
 
  template<typename Left1,typename Right1, typename Right, bool VolumeIntegral1,Integer QR1, Integer Kind>
  class ClassHelper<Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,Expression<Right> >,Kind  >
  {
  public:
   using type=RemoveTupleDuplicates< TupleCatType< typename KindType<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>,Kind>::type,
                                                   typename ClassHelper<Right,Kind>::type >>;
  };  

 
  template<typename Left1,typename Right1, typename Right, bool VolumeIntegral1,Integer QR1, Integer Kind>
  class ClassHelper<Subtraction<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,Expression<Right> >,Kind  >
  {
  public:
   using type=RemoveTupleDuplicates< TupleCatType< typename KindType<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>,Kind>::type,
                                                   typename ClassHelper<Right,Kind>::type >>;
  };  

  
  template<typename Left,typename Left1,typename Right1, bool VolumeIntegral1,Integer QR1, Integer Kind>
  class ClassHelper<Addition<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>> >,Kind  >
  {
  public:
   using type=RemoveTupleDuplicates< TupleCatType< typename ClassHelper<Left,Kind>::type,
                                                   typename KindType<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>,Kind>::type  >>;
  }; 
 

  template<typename Left,typename Left1,typename Right1, bool VolumeIntegral1,Integer QR1, Integer Kind>
  class ClassHelper<Subtraction<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>> >,Kind  >
  {
  public:
   using type=RemoveTupleDuplicates< TupleCatType< typename ClassHelper<Left,Kind>::type,
                                                   typename KindType<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>,Kind>::type  >>;
  } ; 

  template<typename Left,typename Right, Integer Kind>
  class ClassHelper<Addition<Expression<Left>,Expression<Right> >,Kind >
  {
  public:
   using type=RemoveTupleDuplicates< TupleCatType< typename ClassHelper<Left,Kind>::type,
                                                   typename ClassHelper<Right,Kind>::type  >>;
  };   

  template<typename Left,typename Right, Integer Kind>
  class ClassHelper<Subtraction<Expression<Left>,Expression<Right> >,Kind  >
  {
  public:
   using type=RemoveTupleDuplicates< TupleCatType< typename ClassHelper<Left,Kind>::type,
                                                   typename ClassHelper<Right,Kind>::type  >>;
  };   


  template<Integer Kind>
  using type=typename ClassHelper<Form,Kind>::type;

  using TupleFunctionSpace=typename ClassHelper<Form,0>::type;      

  using FunctionSpace=GetType<TupleFunctionSpace,0>;  

  using UniqueElementFunctionSpacesTupleType=typename ClassHelper<Form,1>::type;      

  using form=typename ClassHelper<Form,2>::type;      

  using TupleOfPairsNumbers=BubbleSortTupleOfPairsNumbers<typename TupleOfTestTrialPairsNumbers<Form>::type>;

    GeneralForm(const Form& form)
    : 
    form_(form),
    spaces_ptr_(find_spaces_ptr<FunctionSpace>(form))
    {};

    // GeneralForm(const Form& form,const FullSpace& space)
    // : 
    // form_(form),
    // spaces_ptr_(std::make_shared<FullSpace>(space))
    // {};

    const Form& operator()()const{return form_;};
          Form& operator()()     {return form_;};

    // const std::shared_ptr<FullSpace>& spaces_ptr()const{return spaces_ptr_;};
    //       std::shared_ptr<FullSpace>& spaces_ptr()     {return spaces_ptr_;};

    auto spaces_ptr()     {return spaces_ptr_;}
    auto spaces_ptr()const{return spaces_ptr_;}
 
  private:
  Form form_;
  std::shared_ptr<FunctionSpace> spaces_ptr_;
};
   


template<typename Form>
constexpr auto general_form(const Form& form){return GeneralForm<Form>(form);}






}
#endif