#ifndef MARS_L2_DOT_PRODUCT_INTEGRAL_HPP
#define MARS_L2_DOT_PRODUCT_INTEGRAL_HPP
#include "mars_base.hpp"
#include "mars_l2_dot_product_integral_utils.hpp"

namespace mars {



template<typename Left_,typename Right_,Integer QR=GaussianQuadrature>
class L2DotProductIntegral: 
public Expression<L2DotProductIntegral<Left_,Right_,QR>>
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
   
   static_assert(IsVolumeOrSurfaceIntegral<type>::volume && "In Volume integrals, no trace operator can occur");
    

    static constexpr Integer TestOrTrialLeftValue =GetType<typename TestOrTrialLeft::type,0>::value;
    static constexpr Integer TestOrTrialRightValue =GetType<typename TestOrTrialRight::type,0>::value;

    using Elem=Choose<typename TestOrTrialLeft::Elem,typename TestOrTrialRight::Elem>;
    static constexpr Integer Order=CheckMaxQuadratureOrder<Elem,QR,QuadratureOrder<type>::value+1>::value; 
    using QRule=typename QuadratureRule<QR>:: template rule<Elem,Order>;


    using form= std::tuple<typename TypeOfForm<GetType<typename IsTestOrTrial<Left>::type,0>,
                                               GetType<typename IsTestOrTrial<Right>::type,0>
                            >::type >;
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
    label_(label)
    {}
     

    L2DotProductIntegral(const Expression<L2DotProductIntegral<Left,Right,QR>>& l2prod,const Integer label=-666)
    :
    left_(l2prod.derived().left()),
    right_(l2prod.derived().right()),
    product_(Inner(left_,right_)),
    label_(label)
    {}


     auto operator()(){return L2DotProductIntegral<Left,Right,QR>(left_,right_);}
     const auto operator()()const{return L2DotProductIntegral<Left,Right,QR>(left_,right_);}
  private:
    Left left_;
    Right right_;
    type product_;
    Integer label_;
};






template<typename Left,typename Right>
auto
L2Inner(const Expression<Left>& left,const Expression<Right>& right)
{return L2DotProductIntegral<Left,Right>(left,right);}

// template<typename Left,typename Right>
// auto
// L2Inner(const Expression<Left>& left,const Expression<Right>& right, const Integer label)
// {return L2DotProductIntegral<Left,Right>(left,right,label);}


template<typename Left,typename Right>
constexpr auto
L2Inner(const Expression<UnaryMinus<Expression<Left>>>& left,const Expression<UnaryMinus<Expression<Right>>>& right)
{return L2DotProductIntegral<Left,Right>(left.derived().derived(),right.derived().derived());}

template<typename Left,typename Right>
constexpr auto
L2Inner(const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Left>>> >>& left,const Expression<Right>& right)
{return L2DotProductIntegral<Left,Right>(left.derived().derived().derived(),right);}


template<typename Left,typename Right>
constexpr auto
L2Inner(const Expression<Left>& left,const  Expression<UnaryMinus<Expression<UnaryMinus<Expression<Right>>> >>& right)
{return L2DotProductIntegral<Left,Right>(left,right.derived().derived().derived());}

template<typename Left,typename Right>
constexpr auto
L2Inner(const Expression<UnaryMinus<Expression<Left>>>& left,const  Expression<UnaryMinus<Expression<UnaryMinus<Expression<Right>>> >>& right)
{return L2DotProductIntegral<UnaryMinus<Expression<Left>>,Right>(left.derived(),right.derived().derived().derived());}

template<typename Left,typename Right>
constexpr auto
L2Inner(const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Left>>> >>& left,const Expression<UnaryMinus<Expression<Right>>>& right)
{return L2DotProductIntegral<Left,UnaryMinus<Expression<Right>>>(left.derived().derived().derived(),right.derived());}

template<typename Left,typename Right>
constexpr auto
L2Inner(const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Left>>>>>& left,const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Right>>> >>& right)
{return L2DotProductIntegral<Left,Right>(left.derived().derived().derived(),right.derived().derived().derived());}




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
template< typename Left,typename Right, Integer QR>
class L2DotProductIntegral<Left,Right,QR >
operator+(const Expression<L2DotProductIntegral<Left,Right,QR>>&l2prod)
{return L2Inner(l2prod.derived().left(),l2prod.derived().right());}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////  - L2DotProductIntegral
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< typename Left,typename Right, Integer QR>
class L2DotProductIntegral<UnaryMinus<Expression<Left>>,Right,QR >
operator-(const Expression<L2DotProductIntegral<Left,Right,QR>>&l2prod)
{return L2Inner(-l2prod.derived().left(),l2prod.derived().right());}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////  Exr - L2DotProductIntegral
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< typename Left1,typename Left2,typename Right2, Integer QR>
class Addition< Expression <Left1>, 
                Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR>> >
operator-(const Expression<Left1>&left, 
          const Expression<L2DotProductIntegral<Left2,Right2,QR>>&right)
{return left+L2Inner(-right.derived().left(),right.derived().right());}




template< typename Left1,typename Left2,typename Right2, Integer QR>
class Addition< Expression <Left1>, 
                Expression<L2DotProductIntegral<Left2,Right2,QR>> >
operator-(const Expression<Left1>&left, 
          const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR>>&right)
{return left+L2Inner(right.derived().left().derived(),right.derived().right());}



template< typename Left1,typename Left2,typename Right2, Integer QR2>
class Addition< Expression<Left1>, 
                Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
operator-(const Expression<Left1>&left, 
          const Expression<L2DotProductIntegral<Left2,UnaryMinus<Expression<Right2>>,QR2>>&right)
{return left+L2Inner(right.derived().left(),right.derived().right().derived());}




template< typename Left1,typename Left2,typename Right2, Integer QR2>
class Addition< Expression<Left1>, 
                Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>> >
operator-(const Expression<Left1>&left, 
          const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,UnaryMinus<Expression<Right2>>,QR2>>&right)
{return left+L2Inner(right.derived().left(),right.derived().right().derived());}


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






}
#endif
