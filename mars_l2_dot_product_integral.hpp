#ifndef MARS_L2_DOT_PRODUCT_INTEGRAL_HPP
#define MARS_L2_DOT_PRODUCT_INTEGRAL_HPP
#include "mars_base.hpp"
#include "mars_l2_dot_product_integral_utils.hpp"
#include "mars_quadrature_order.hpp"
namespace mars {



template<typename Left_,typename Right_,bool VolumeIntegral=true,Integer QR=GaussianQuadrature>
class L2DotProductIntegral: 
public Expression<L2DotProductIntegral<Left_,Right_,VolumeIntegral,QR>>
{  
   public:
    using Left=Left_;
    using Right=Right_;
    using type=InnerProduct<Expression <Left>, Expression <Right> > ;
    using TestOrTrialLeft= IsTestOrTrial<Left>;
    using TestOrTrialRight= IsTestOrTrial<Right>;
    using OperatorLeft=typename TestOrTrialLeft::Operator;
    using OperatorRight=typename TestOrTrialRight::Operator;
    static constexpr Integer leftN=TestOrTrialLeft::number;
    static constexpr Integer rightN=TestOrTrialRight::number;
    
   static_assert(IsVolumeOrSurfaceIntegral<VolumeIntegral,type>::value && "In Volume integrals, no trace operator can occur");
    

    static constexpr Integer TestOrTrialLeftValue =GetType<typename TestOrTrialLeft::type,0>::value;
    static constexpr Integer TestOrTrialRightValue =GetType<typename TestOrTrialRight::type,0>::value;

    using Elem=Choose<typename TestOrTrialLeft::Elem,typename TestOrTrialRight::Elem>;
    using QuadratureElem=typename VolumeOrSurfaceElem<Elem,VolumeIntegral>::type;
    
    static constexpr Integer Order=CheckMaxQuadratureOrder<QuadratureElem,QR,QuadratureOrder<type>::value+1>::value; 
    using QRule=typename QuadratureRule<QR>:: template rule<QuadratureElem,Order>;

    using TestOrTrialLeftType=GetType<typename TestOrTrialLeft::type>;
    using TestOrTrialRightType=GetType<typename TestOrTrialRight::type>;

    using form= std::tuple<typename TypeOfForm<TestOrTrialLeftType,TestOrTrialRightType>::type >;
    using TestTrialNumbers=typename FormTestTrialNumbers<GetType<form,0>::value,TestOrTrialLeftValue,TestOrTrialRightValue,leftN,rightN>::type;

    using UniqueElementFunctionSpacesTupleType=GetType<RemoveTupleDuplicates< TupleCatType< typename TestOrTrialLeft::UniqueElementFunctionSpacesTupleType,
                                                                                            typename TestOrTrialRight::UniqueElementFunctionSpacesTupleType  >>,0>;               
    using TupleFunctionSpace=RemoveTupleDuplicates< TupleCatType< typename IsTestOrTrial<Left>::TupleFunctionSpace,
                                                                  typename IsTestOrTrial<Right>::TupleFunctionSpace  >>;               
    
    using FunctionSpace=GetType<TupleFunctionSpace,0>;

    using TupleOfSpaces=typename FunctionSpace::TupleOfSpaces;


    template<typename TupleOfSpacesAux,typename TestTrialNumbersAux,Integer N>
    class ClassAux;

    template<typename TupleOfSpacesAux,typename TestTrialNumbersAux>
    class ClassAux<TupleOfSpacesAux,TestTrialNumbersAux,0>
    {
    public:
      using type=std::tuple<Real>;
    };

    template<typename TupleOfSpacesAux,typename TestTrialNumbersAux>
    class ClassAux<TupleOfSpacesAux,TestTrialNumbersAux,1>
    {
    public:
      using type_test= typename std::conditional< (-1==GetType<TestTrialNumbersAux,0>::value),
                                          std::tuple<>,
                                          GetType<TupleOfSpacesAux,GetType<TestTrialNumbersAux,0>::value>>::type;
      using type=std::tuple<type_test>;
    };

    template<typename TupleOfSpacesAux,typename TestTrialNumbersAux>
    class ClassAux<TupleOfSpacesAux,TestTrialNumbersAux,2>
    {
    public:
      using type_test= typename std::conditional< (-1==GetType<TestTrialNumbersAux,0>::value),
                                          std::tuple<>,
                                          GetType<TupleOfSpacesAux,GetType<TestTrialNumbersAux,0>::value>>::type;
      using type_trial=typename std::conditional< (-1==GetType<TestTrialNumbersAux,1>::value),
                                          std::tuple<>,
                                          GetType<TupleOfSpacesAux,GetType<TestTrialNumbersAux,1>::value>>::type;
      using type=std::tuple<type_test,type_trial>;
    };

   using TestTrialSpaces=typename ClassAux<TupleOfSpaces,TestTrialNumbers,GetType<form,0>::value>::type;
   
    const Left&  left() const{return left_;};
    const Right& right()const{return right_;};

    L2DotProductIntegral(const Expression<Left>& left,const Expression<Right>& right,const Integer label=-666)
    :
    left_(left.derived()),
    right_(right.derived()),
    product_(Inner(left,right)),
    label_(label),
    spaces_ptr_(find_spaces_ptr<FunctionSpace>(left_*right_))
    {}
     

    L2DotProductIntegral(const Expression<L2DotProductIntegral<Left,Right,VolumeIntegral,QR>>& l2prod)
    :
    left_(l2prod.derived().left()),
    right_(l2prod.derived().right()),
    product_(Inner(left_,right_)),
    label_(l2prod.label()),
    spaces_ptr_(l2prod.spaces_ptr())//find_spaces_ptr<FunctionSpace>(left_*right_))
    {}

     auto spaces_ptr()     {return spaces_ptr_;}
     auto spaces_ptr()const{return spaces_ptr_;}


     auto label()const{return label_;}
 

     auto operator()(){return L2DotProductIntegral<Left,Right,VolumeIntegral,QR>(left_,right_);}
     const auto operator()()const{return L2DotProductIntegral<Left,Right,VolumeIntegral,QR>(left_,right_);}
  private:
    Left left_;
    Right right_;
    type product_;
    Integer label_;
    std::shared_ptr<FunctionSpace> spaces_ptr_;
};



template<bool VolumeIntegral,Integer QR,typename Left,typename Right>
auto
L2InnerProduct(const Expression<Left>& left,const Expression<Right>& right)
{return L2DotProductIntegral<Left,Right,VolumeIntegral,QR>(left,right);}

template<bool VolumeIntegral,Integer QR,typename Left,typename Right>
auto
L2InnerProduct(const Expression<Left>& left,const Expression<Right>& right, const Integer label)
{return L2DotProductIntegral<Left,Right,VolumeIntegral,QR>(left,right,label);}









template<typename Left,typename Right,bool VolumeIntegral=true,Integer QR=GaussianQuadrature>
auto
L2Inner(const Expression<Left>& left,const Expression<Right>& right)
{return L2DotProductIntegral<Left,Right,VolumeIntegral,QR>(left,right);}

template<typename Left,typename Right,Integer QR=GaussianQuadrature>
auto
surface_integral(const Integer label,const Expression<Left>& left,const Expression<Right>& right)
{return L2DotProductIntegral<Left,Right,false,QR>(left,right,label);}

// template<typename Left,typename Right,Integer QR=GaussianQuadrature>
// auto
// surface_integral(const Expression<Left>& left,const Expression<Right>& right)
// { assert(1&&"Add a label for the surface integral");
//   return L2DotProductIntegral<Left,Right,false,QR>(left,right,-1);}


// template<typename Left,typename Right>
// auto
// L2Inner(const Expression<Left>& left,const Expression<Right>& right, const Integer label)
// {return L2DotProductIntegral<Left,Right>(left,right,label);}


template<typename Left,typename Right,bool VolumeIntegral=true,Integer QR=GaussianQuadrature>
constexpr auto
L2Inner(const Expression<UnaryMinus<Expression<Left>>>& left,const Expression<UnaryMinus<Expression<Right>>>& right)
{return L2DotProductIntegral<Left,Right,VolumeIntegral,QR>(left.derived().derived(),right.derived().derived());}

template<typename Left,typename Right,bool VolumeIntegral=true,Integer QR=GaussianQuadrature>
constexpr auto
L2Inner(const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Left>>> >>& left,const Expression<Right>& right)
{return L2DotProductIntegral<Left,Right,VolumeIntegral,QR>(left.derived().derived().derived(),right);}


template<typename Left,typename Right,bool VolumeIntegral=true,Integer QR=GaussianQuadrature>
constexpr auto
L2Inner(const Expression<Left>& left,const  Expression<UnaryMinus<Expression<UnaryMinus<Expression<Right>>> >>& right)
{return L2DotProductIntegral<Left,Right,VolumeIntegral,QR>(left,right.derived().derived().derived());}

template<typename Left,typename Right,bool VolumeIntegral=true,Integer QR=GaussianQuadrature>
constexpr auto
L2Inner(const Expression<UnaryMinus<Expression<Left>>>& left,const  Expression<UnaryMinus<Expression<UnaryMinus<Expression<Right>>> >>& right)
{return L2DotProductIntegral<UnaryMinus<Expression<Left>>,Right,VolumeIntegral,QR>(left.derived(),right.derived().derived().derived());}

template<typename Left,typename Right,bool VolumeIntegral=true,Integer QR=GaussianQuadrature>
constexpr auto
L2Inner(const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Left>>> >>& left,const Expression<UnaryMinus<Expression<Right>>>& right)
{return L2DotProductIntegral<Left,UnaryMinus<Expression<Right>>,VolumeIntegral,QR>(left.derived().derived().derived(),right.derived());}

template<typename Left,typename Right,bool VolumeIntegral=true,Integer QR=GaussianQuadrature>
constexpr auto
L2Inner(const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Left>>>>>& left,const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Right>>> >>& right)
{return L2DotProductIntegral<Left,Right,VolumeIntegral,QR>(left.derived().derived().derived(),right.derived().derived().derived());}




template<typename Left2,typename Right2, typename Left>
constexpr auto
L2Inner(const Expression<Left>& left,const Expression<Addition<Expression<Left2>,Expression<Right2>>>& right)
{return 
  Addition<
  Expression<decltype(L2Inner(left.derived(),right.derived().left()))>,
  Expression<decltype(L2Inner(left.derived(),right.derived().right()))>>
              (L2Inner(left.derived(),right.derived().left()),
               L2Inner(left.derived(),right.derived().right()) );}

template<typename Left2,typename Right2, typename Left>
constexpr auto
L2Inner(const Expression<Left>& left,const Expression<Subtraction<Expression<Left2>,Expression<Right2>>>& right)
{return Addition<
  Expression<decltype(L2Inner(left.derived(),right.derived().left()))>,
  Expression<decltype(L2Inner(-left.derived(),right.derived().right()))>>
              (L2Inner(left.derived(),right.derived().left()),
               L2Inner(-left.derived(),right.derived().right()) );}







template<typename Left1,typename Right1, typename Right>
constexpr auto
L2Inner(const Expression<Addition<Expression<Left1>,Expression<Right1>>>& left,const Expression<Right>& right)
{return Addition<
  Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
  Expression<decltype(L2Inner(left.derived().right(),right.derived()))>>
  (L2Inner(left.derived().left(),right.derived()),
   L2Inner(left.derived().right(),right.derived()) );}

template<typename Left1,typename Right1, typename Right>
constexpr auto
L2Inner(const Expression<Subtraction<Expression<Left1>,Expression<Right1>>>& left,const Expression<Right>& right)
{return Addition<
  Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
  Expression<decltype(L2Inner(-left.derived().right(),right.derived()))>>
  (L2Inner(left.derived().left(),right.derived()),
   L2Inner(-left.derived().right(),right.derived()) );}


template<typename Left1,typename Right1,typename Left2, typename Right2>
constexpr auto
L2Inner(const Expression<Addition<Expression<Left1>,Expression<Right1>>>& left,
        const Expression<Addition<Expression<Left2>,Expression<Right2>>>& right)
{return 
Addition<Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
         Expression<decltype(L2Inner(left.derived().right(),right.derived()))>
         >
  (
    L2Inner(left.derived().left(),right.derived()),
    L2Inner(left.derived().right(),right.derived())                
  )
  ;
}


template<typename Left1,typename Right1,typename Left2, typename Right2>
constexpr auto
L2Inner(const Expression<Addition<Expression<Left1>,Expression<Right1>>>& left,
        const Expression<Subtraction<Expression<Left2>,Expression<Right2>>>& right)
{return 
Addition<Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
         Expression<decltype(L2Inner(left.derived().right(),right.derived()))>
         >
  (
    L2Inner(left.derived().left(),right.derived()),
    L2Inner(left.derived().right(),right.derived())                
  )
  ;
}

template<typename Left1,typename Right1,typename Left2, typename Right2>
constexpr auto
L2Inner(const Expression<Subtraction<Expression<Left1>,Expression<Right1>>>& left,
        const Expression<Subtraction<Expression<Left2>,Expression<Right2>>>& right)
{return 
Addition<Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
         Expression<decltype(L2Inner(-left.derived().right(),right.derived()))>
         >
  (
    L2Inner(left.derived().left(),right.derived()),
    L2Inner(-left.derived().right(),right.derived())                
  )
  ;
}


template<typename Left1,typename Right1,typename Left2, typename Right2>
constexpr auto
L2Inner(const Expression<Subtraction<Expression<Left1>,Expression<Right1>>>& left,
        const Expression<Addition<Expression<Left2>,Expression<Right2>>>& right)
{return 
Addition<Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
         Expression<decltype(L2Inner(-left.derived().right(),right.derived()))>
         >
  (
    L2Inner(left.derived().left(),right.derived()),
    L2Inner(-left.derived().right(),right.derived())                
  )
  ;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////  + L2DotProductIntegral
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< typename Left,typename Right,bool VolumeIntegral,Integer QR>
class L2DotProductIntegral<Left,Right,VolumeIntegral,QR >
operator+(const Expression<L2DotProductIntegral<Left,Right,VolumeIntegral,QR>>&l2prod)
{return L2InnerProduct<VolumeIntegral,QR>(l2prod.derived().left(),l2prod.derived().right(),l2prod.derived().label());}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////  - L2DotProductIntegral
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< typename Left,typename Right,bool VolumeIntegral,Integer QR>
class L2DotProductIntegral<UnaryMinus<Expression<Left>>,Right,VolumeIntegral,QR >
operator-(const Expression<L2DotProductIntegral<Left,Right,VolumeIntegral,QR>>&l2prod)
{return L2InnerProduct<VolumeIntegral,QR>(-l2prod.derived().left(),l2prod.derived().right());}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////  Exr - L2DotProductIntegral
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< typename Left1,typename Left2,typename Right2, bool VolumeIntegral2, Integer QR2>
class Addition< Expression <Left1>, 
                Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,VolumeIntegral2,QR2>> >
operator-(const Expression<Left1>&left, 
          const Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>>&right)
{return left+L2InnerProduct<VolumeIntegral2,QR2>(-right.derived().left(),right.derived().right(),
                                                right.derived().label());}




template< typename Left1,typename Left2,typename Right2, bool VolumeIntegral2, Integer QR2>
class Addition< Expression <Left1>, 
                Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>> >
operator-(const Expression<Left1>&left, 
          const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,VolumeIntegral2,QR2>>&right)
{return left+L2InnerProduct<VolumeIntegral2,QR2>(right.derived().left().derived(),right.derived().right(),
                                                 right.derived().label());}



template< typename Left1,typename Left2,typename Right2, bool VolumeIntegral2,Integer QR2>
class Addition< Expression<Left1>, 
                Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>> >
operator-(const Expression<Left1>&left, 
          const Expression<L2DotProductIntegral<Left2,UnaryMinus<Expression<Right2>>,VolumeIntegral2,QR2>>&right)
{return left+L2InnerProduct<VolumeIntegral2,QR2>(right.derived().left(),right.derived().right().derived(),
                                                 right.derived().label());}




template< typename Left1,typename Left2,typename Right2, bool VolumeIntegral2,Integer QR2>
class Addition< Expression<Left1>, 
                Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,VolumeIntegral2,QR2>> >
operator-(const Expression<Left1>&left, 
          const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,UnaryMinus<Expression<Right2>>,VolumeIntegral2,QR2>>&right)
{return left+L2InnerProduct<VolumeIntegral2,QR2>(right.derived().left(),right.derived().right().derived(),
                                                 right.derived().label());}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////  L2DotProductIntegral - Exr
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// template< typename Left1,typename Right1, Integer QR, typename Right2>
// class Addition< Expression<L2DotProductIntegral<Left1,Right1,QR>>, Expression<UnaryMinus<Expression<Right2>>> >
// operator-(const Expression<L2DotProductIntegral<Left1,Right1,QR>>&left, 
//           const Expression<Right2>&right)
// {return Addition<Expression<L2DotProductIntegral<UnaryMinus<Expression<Left1>>,Right1,QR>>,
//                  Expression <UnaryMinus<Expression<Right2>>> >
//                  (left,-right);}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //////  L2DotProductIntegral - L2DotProductIntegral
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename Left1,typename Right1, Integer QR1, typename Left2,typename Right2, Integer QR2>
// class Addition< Expression<L2DotProductIntegral<Left1,Right1,QR1>>, 
//                 Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>> >
// operator-(const Expression<L2DotProductIntegral<Left1,Right1,QR1>>&left, 
//           const Expression<L2DotProductIntegral<Left2,Right2,QR2>>&right)
// {return left+L2Inner(-right.derived().left(),right.derived().right());}
//   // return Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//   //                Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>> >
//   //                (left,L2Inner(-right.derived().left(),right.derived().right()));}





// template< typename Left1,typename Right1, Integer QR1, typename Left2,typename Right2, Integer QR2>
// class Addition< Expression<L2DotProductIntegral<Left1,Right1,QR1>>, 
//                 Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// operator-(const Expression<L2DotProductIntegral<Left1,Right1,QR1>>&left, 
//           const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>>&right)
// {
// return left+L2Inner(right.derived().left().derived(),right.derived().right());
//   // return Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//   //                Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
//   //                (left,L2Inner(right.derived().left().derived(),right.derived().right()));
//                }

// template< typename Left1,typename Left2,typename Right2, Integer QR2>
// class Addition< Expression<Left1>, 
//                 Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// operator-(const Expression<Left1>&left, 
//           const Expression<L2DotProductIntegral<Left2,UnaryMinus<Expression<Right2>>,QR2>>&right)
// {return left+L2Inner(right.derived().left(),right.derived().right().derived());
//   // Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//   //                Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
//   //                (left,L2Inner(right.derived().left().derived(),right.derived().right()));
//                }


// template< typename Left1,typename Left2,typename Right2, Integer QR2>
// class Addition< Expression<Left1>, 
//                 Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>> >
// operator-(const Expression<Left1>&left, 
//           const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,UnaryMinus<Expression<Right2>>,QR2>>&right)
// {return left+L2Inner(right.derived().left(),right.derived().right().derived());
//   // Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//   //                Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
//   //                (left,L2Inner(right.derived().left().derived(),right.derived().right()));
//                }






//   Addition<Expression<Left1>,
//                  Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>> >
//                  (left,L2Inner(right.derived().left().derived(),right.derived().right().derived()));}



// template< typename Left1,typename Right1, Integer QR1, typename Left2,typename Right2, Integer QR2>
// class Addition< Expression<L2DotProductIntegral<Left1,Right1,QR1>>, 
//                 Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// operator-(const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left1>>,Right1,QR1>>&left, 
//           const Expression<L2DotProductIntegral<Left2,Right2,QR2>>&right)
// {return Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//                  Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
//                  (left,L2Inner(right.derived().left(),right.derived().right().derived()));}




// // template<typename Left1,typename Right1, typename Right>
// // Addition<Expression<L2DotProductIntegral<MeshT,Left1,Right>>,
// //           Expression<L2DotProductIntegral<MeshT,Right1,Right>>>
// // L2Inner(const Expression<Addition<Expression<Left1>,Expression<Right1>>>& left,const Expression<Right>& right, const Integer& label)
// // {return Addition<
// //   Expression<L2DotProductIntegral<MeshT,Left1,Right>>,
// //   Expression<L2DotProductIntegral<MeshT,Right1,Right>>>(L2Inner(left.derived().left(),right,label),
// //                                                         L2Inner(left.derived().right(),right,label) );}

template<Integer H,typename T>
class ExtractFormType;

template< typename Left, typename Right,bool VolumeIntegral,Integer QR>
class OperatorTypeHelper<ExtractFormType<-1, L2DotProductIntegral<Left,Right,VolumeIntegral,QR> >  >
{ public:
  using type=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
};

template< typename Left, typename Right, Integer QR>
class OperatorTypeHelper<ExtractFormType<0, L2DotProductIntegral<Left,Right,true,QR> >  >
{ public:
  using type=L2DotProductIntegral<Left,Right,true,QR>;
};

template< typename Left, typename Right, Integer QR>
class OperatorTypeHelper<ExtractFormType<0, L2DotProductIntegral<Left,Right,false,QR> >  >
{ public:
  using type=EmptyExpression;
};

template< typename Left, typename Right, Integer QR>
class OperatorTypeHelper<ExtractFormType<1, L2DotProductIntegral<Left,Right,false,QR> >  >
{ public:
  using type=L2DotProductIntegral<Left,Right,false,QR>;
};

template< typename Left, typename Right, Integer QR>
class OperatorTypeHelper<ExtractFormType<1, L2DotProductIntegral<Left,Right,true,QR> >  >
{ public:
  using type=EmptyExpression;
};


template<Integer H>
class OperatorTypeHelper<ExtractFormType<H, Addition<Expression<EmptyExpression>,Expression<EmptyExpression> > >  >
{ public:
  using type=EmptyExpression;
};

template<Integer H,typename Left>
class OperatorTypeHelper<ExtractFormType<H, Addition<Expression<Left>,Expression<EmptyExpression> > >  >
{ public:
  using type=Left;
};

template<Integer H,typename Right>
class OperatorTypeHelper<ExtractFormType<H, Addition<Expression<EmptyExpression>,Expression<Right> > >  >
{ public:
  using type=Right;
};


template<Integer H, typename Left, typename Right>
class OperatorTypeHelper<ExtractFormType<H, Addition<Expression<Left>,Expression<Right> > >  >
{ public:
  using LeftT=OperatorType<ExtractFormType<H,Left>>;
  using RightT=OperatorType<ExtractFormType<H,Right>>;
  using type= std::conditional_t<(IsSame<LeftT,EmptyExpression>::value && IsSame<RightT,EmptyExpression>::value),
                                 EmptyExpression,
                                 std::conditional_t<IsSame<LeftT,EmptyExpression>::value,
                                                    RightT,
                                                    std::conditional_t<IsSame<RightT,EmptyExpression>::value,
                                                                       LeftT,
                                                                       Addition<Expression<LeftT>,Expression<RightT>>
                                                                      >
                                                   >
                                >;
};



template<Integer H, typename Left,typename Right,bool VolumeIntegral,Integer QR>
constexpr std::enable_if_t<H==-1,L2DotProductIntegral<Left,Right,VolumeIntegral,QR>>
ExtractForm(const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod)
{return l2prod;}

template<Integer H, typename Left,typename Right,Integer QR>
constexpr std::enable_if_t<H==0,L2DotProductIntegral<Left,Right,true,QR>>
ExtractForm(const L2DotProductIntegral<Left,Right,true,QR>& l2prod)
{return l2prod;}

template<Integer H, typename Left,typename Right,Integer QR>
constexpr std::enable_if_t<H==0,EmptyExpression>
ExtractForm(const L2DotProductIntegral<Left,Right,false,QR>& l2prod)
{return EmptyExpression();}


template<Integer H, typename Left,typename Right,Integer QR>
constexpr std::enable_if_t<H==1,EmptyExpression>
ExtractForm(const L2DotProductIntegral<Left,Right,true,QR>& l2prod)
{return EmptyExpression();;}

template<Integer H, typename Left,typename Right,Integer QR>
constexpr std::enable_if_t<H==1,L2DotProductIntegral<Left,Right,false,QR>>
ExtractForm(const L2DotProductIntegral<Left,Right,false,QR>& l2prod)
{return l2prod;}

template<Integer H>
constexpr auto
ExtractForm(const EmptyExpression& t)
{return t;}



template<Integer H,typename T>
constexpr auto
ExtractForm(const UnaryPlus<Expression<T>>& t)
{return ExtractForm<H>(t);}

template<Integer H,typename T>
constexpr auto
ExtractForm(const UnaryMinus<Expression<T>>& t)
{return ExtractForm<H>(t);}


template<Integer H>
constexpr auto
ExtractForm(const Addition<Expression<EmptyExpression>,Expression<EmptyExpression>>& sum)
{return EmptyExpression();}

template<Integer H,typename Left>
constexpr OperatorType<ExtractFormType<H,Addition<Expression<Left>,Expression<EmptyExpression>>>>
ExtractForm(const Addition<Expression<Left>,Expression<EmptyExpression>>& sum)
{return (ExtractForm<H>(sum.left()));}

template<Integer H,typename Right>
constexpr OperatorType<ExtractFormType<H,Addition<Expression<EmptyExpression>,Expression<Right>>>>
ExtractForm(const Addition<Expression<EmptyExpression>,Expression<Right>>& sum)
{return (ExtractForm<H>(sum.right()));}

template<Integer H,typename Left,typename Right>
constexpr OperatorType<ExtractFormType<H,Addition<Expression<Left>,Expression<Right>>>>
ExtractForm(const Addition<Expression<Left>,Expression<Right>>& sum)
{return ExtractForm<H>(ExtractForm<H>(sum.left())+ExtractForm<H>(sum.right()));}



}
#endif
