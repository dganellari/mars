#ifndef MARS_L2_DOT_PRODUCT_INTEGRAL_UTILS_HPP
#define MARS_L2_DOT_PRODUCT_INTEGRAL_UTILS_HPP
#include "mars_base.hpp"
#include "mars_constant.hpp"
#include "mars_function.hpp"
#include "mars_operators.hpp"
#include "mars_trial_function.hpp"
#include "mars_test_function.hpp"

namespace mars {



template<typename T>
class IsVolumeOrSurfaceIntegralHelper
{
public:
  static constexpr bool surface=true;
  static constexpr bool volume=true;
};

template<template<class,Integer,class>class TestOrTrial, typename MixedSpace,Integer N>
class IsVolumeOrSurfaceIntegralHelper<TestOrTrial<MixedSpace,N,TraceOperator>>{
public:
  static constexpr bool surface=true;
  static constexpr bool volume=false;
};

// template<template<class,Integer,class>class TestOrTrial, typename MixedSpace,Integer N, typename OperatorType>
// class IsVolumeOrSurfaceIntegralHelper<TestOrTrial<MixedSpace,N,OperatorType>>{
// public:
//   static constexpr bool surface=true;
//   static constexpr bool volume=true;
// };

template<template<class,Integer,class>class TestOrTrial, typename MixedSpace,Integer N, typename OperatorType>
class IsVolumeOrSurfaceIntegralHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<OperatorType>>>>
{
public:
  using T=TestOrTrial<MixedSpace,N,OperatorType>;
  static constexpr bool surface=IsVolumeOrSurfaceIntegralHelper<T>::surface;
  static constexpr bool volume=IsVolumeOrSurfaceIntegralHelper<T>::volume;
};


template<template<class,Integer,class>class TestOrTrial,template<class>class Unary, typename MixedSpace,Integer N, typename OperatorType>
class IsVolumeOrSurfaceIntegralHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Unary<Expression<OperatorType>>>>> >{
public:
  using T=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<OperatorType>>>;
  static constexpr bool surface=IsVolumeOrSurfaceIntegralHelper<T>::surface;
  static constexpr bool volume=IsVolumeOrSurfaceIntegralHelper<T>::volume;
};

template<template<class,Integer,class>class TestOrTrial,template<class,class>class Binary, typename MixedSpace,Integer N, typename Left,typename Right>
class IsVolumeOrSurfaceIntegralHelper<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Binary<Expression<Left>,Expression<Right>>>>> >{
public:
  using LeftT=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Left>>>;
  using RightT=TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Right>>>;

  static constexpr bool surface=IsVolumeOrSurfaceIntegralHelper<LeftT>::surface*IsVolumeOrSurfaceIntegralHelper<RightT>::surface;
  static constexpr bool volume=IsVolumeOrSurfaceIntegralHelper<LeftT>::volume*IsVolumeOrSurfaceIntegralHelper<RightT>::volume;
};

template<template<class>class Unary, typename Type>
class IsVolumeOrSurfaceIntegralHelper< Unary<Expression<Type>> >
{
public:
  static constexpr bool surface=IsVolumeOrSurfaceIntegralHelper<Type>::surface;
  static constexpr bool volume=IsVolumeOrSurfaceIntegralHelper<Type>::volume;
};


template<template<class,class>class Binary, typename Left,typename Right>
class IsVolumeOrSurfaceIntegralHelper< Binary<Expression<Left>,Expression<Right>> >
{
public:
  static constexpr bool surface=IsVolumeOrSurfaceIntegralHelper<Left>::surface*IsVolumeOrSurfaceIntegralHelper<Right>::surface;
  static constexpr bool volume=IsVolumeOrSurfaceIntegralHelper<Left>::volume*IsVolumeOrSurfaceIntegralHelper<Right>::volume;
};


template<typename Left,typename Right>
class IsVolumeOrSurfaceIntegralHelper<InnerProduct<Expression <Left>, Expression <Right> >  >
{
public:
  static constexpr bool surface=IsVolumeOrSurfaceIntegralHelper<Left>::surface*IsVolumeOrSurfaceIntegralHelper<Right>::surface;
  static constexpr bool volume=IsVolumeOrSurfaceIntegralHelper<Left>::volume*IsVolumeOrSurfaceIntegralHelper<Right>::volume;
};


template<bool VolumeIntegral, typename T>
class IsVolumeOrSurfaceIntegralAux;

template<typename T>
class IsVolumeOrSurfaceIntegralAux<true,T>
{
public:
  static constexpr bool value=IsVolumeOrSurfaceIntegralHelper<T>::volume;
};

template<typename T>
class IsVolumeOrSurfaceIntegralAux<false,T>
{
public:
  static constexpr bool value=IsVolumeOrSurfaceIntegralHelper<T>::surface;
};

template<bool VolumeIntegral, typename T>
class IsVolumeOrSurfaceIntegral
{
  public:
   static constexpr bool value=IsVolumeOrSurfaceIntegralAux<VolumeIntegral,T>::value;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// We associate Number<2>,Number<1>,Number<0> respectively to a X=Trial/Test/OtherFunction
///// We define IsTestOrTrial<X> which is used for each term Left or Right in L2Inner<Left,Right>
///// The possible results are:
///// 1)IsTestOrTrial<X> == std::tuple<Number<0>> (if only OtherFunction is there, it is not touched)
///// 2)IsTestOrTrial<X> == std::tuple<Number<1>> (if OtherFunction is there, Number<0> is removed and only Number<1> survives)
///// 3)IsTestOrTrial<X> == std::tuple<Number<2>> (if OtherFunction is there, Number<0> is removed and only Number<2> survives)
///// 4)Error with static assert (if more than one Number<1>/Number<2> is present or if both Number<1> and Number<2> are present)
/////   This assert means that we can have only: a function or exactly one trial or exactly one test
///// Then is check the form of L2inner, bu using TypeOfForm<Number<M>,Number<N>>:
///// 0) 0 form if:
/////  Left==Right==std::tuple<Number<0>>
///// 1) 1 form if: 
/////              Left==std::tuple<Number<0>> and  Right==std::tuple<Number<1>> (function,test)
/////              Left==std::tuple<Number<1>> and  Right==std::tuple<Number<0>> (test,function)
///// 2) 2 form if Left==Right==Number<0>
/////              Left==std::tuple<Number<1>> and  Right==std::tuple<Number<2>> (test,trial)
/////              Left==std::tuple<Number<2>> and  Right==std::tuple<Number<1>> (trial,test)
///// 3) everything else: error
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
class IsTestOrTrial;

// template<typename T>
// class IsTestOrTrial{
// public:
//   // using Elem=EmptyClass;
//   // using Operator=std::tuple<>;
//   // using TupleFunctionSpace=std::tuple<>;
//   // using UniqueElementFunctionSpacesTupleType=std::tuple<>;
//   using type=std::tuple<Number<-1>>;
//   // static constexpr Integer value=-1;
//   // static constexpr Integer number=-1;
// };


template<typename ConstType,typename...Inputs>
class IsTestOrTrial< ConstantTensor<ConstType,Inputs...>>{
public:
  using Elem=EmptyClass;
  using Operator=std::tuple<>;
  using TupleFunctionSpace=std::tuple<>;
  using UniqueElementFunctionSpacesTupleType=std::tuple<>;
  using type=std::tuple<Number<0>>;
  static constexpr Integer value=-1;
  static constexpr Integer number=-1;
};
template<typename ConstType,typename...Inputs>
class IsTestOrTrial<const ConstantTensor<ConstType,Inputs...>>
: public IsTestOrTrial< ConstantTensor<ConstType,Inputs...>>
{};


template<typename FuncType,typename FullSpace, Integer N, typename Operator_>
class IsTestOrTrial<Function<FullSpace,N,Operator_,FuncType>>{
public:
  using Elem=typename FullSpace::Elem;
  using Operator=std::tuple<Operator_>;
  using TupleFunctionSpace=std::tuple<>;
  using UniqueElementFunctionSpacesTupleType=std::tuple<>;
  using type=std::tuple<Number<0>>;
  static constexpr Integer value=-1;
  static constexpr Integer number=-1;
};

template<typename MixedSpace,Integer N, typename OperatorType>
class IsTestOrTrial<Test<MixedSpace,N,OperatorType>>{
public:
  using Elem=typename MixedSpace::Elem;
  using Operator=std::tuple<OperatorType>;
  using TupleFunctionSpace=std::tuple<MixedSpace>;
  using UniqueElementFunctionSpacesTupleType=std::tuple<typename MixedSpace::UniqueElementFunctionSpacesTupleType>;
  using type=std::tuple<Number<1>>;
  static constexpr Integer value=Test<MixedSpace,N,OperatorType>::value;
  static constexpr Integer number=Test<MixedSpace,N,OperatorType>::number;
};

template<typename MixedSpace,Integer N, typename OperatorType>
class IsTestOrTrial<Trial<MixedSpace,N,OperatorType>>{
public:
  using Elem=typename MixedSpace::Elem;
  using Operator=std::tuple<OperatorType>;
  using TupleFunctionSpace=std::tuple<MixedSpace>;
  using UniqueElementFunctionSpacesTupleType=std::tuple<typename MixedSpace::UniqueElementFunctionSpacesTupleType>;
  using type=std::tuple<Number<2>>;
  static constexpr Integer value=Trial<MixedSpace,N,OperatorType>::value;
  static constexpr Integer number=Trial<MixedSpace,N,OperatorType>::number;
};

template <typename T,typename ... Types>
class TupleTypeSize;




template<template<class>class Unary, typename Type>
class IsTestOrTrial< Unary<Expression<Type>> >
{
public:
  using Elem=typename IsTestOrTrial<Type>::Elem;

  using Operator=TupleCatType<typename IsTestOrTrial<Type>::Operator>;
  using TupleFunctionSpace=typename IsTestOrTrial<Type>::TupleFunctionSpace;
  using UniqueElementFunctionSpacesTupleType=typename IsTestOrTrial<Type>::UniqueElementFunctionSpacesTupleType;
  // WE REMOVE THE ZEROS, BUT WE SHOULD NOT REMOVE ALL ZEROS
  // if it has only one element, do not remove duplicates
  using type=TupleRemoveNumber0<typename IsTestOrTrial<Type>::type>;
  // using type=TupleRemoveType<Number<0>,typename IsTestOrTrial<Type>::type>;
  static_assert(TupleTypeSize<type>::value<2," In Unary<Type>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
  static constexpr Integer number= Heaviside(IsTestOrTrial<Type>::number);
};


template<template<class,class>class Binary, typename Left, typename Right>
class IsTestOrTrial<Binary<Expression<Left>,Expression<Right> > >
{
public:

  using Elem=Choose<typename IsTestOrTrial<Left>::Elem,typename IsTestOrTrial<Right>::Elem>;
  using Operator=TupleCatType<typename IsTestOrTrial<Left>::Operator,
                              typename IsTestOrTrial<Right>::Operator >;
  using TupleFunctionSpace=TupleCatType<typename IsTestOrTrial<Left>::TupleFunctionSpace,
                                        typename IsTestOrTrial<Right>::TupleFunctionSpace >;

  using UniqueElementFunctionSpacesTupleType=TupleCatType<typename IsTestOrTrial<Left>::UniqueElementFunctionSpacesTupleType,
                                                          typename IsTestOrTrial<Right>::UniqueElementFunctionSpacesTupleType >;
  using tmp_type=TupleRemoveType<Number<0>,TupleCatType<typename IsTestOrTrial<Left>::type,typename IsTestOrTrial<Right>::type >> ;
  using type=typename std::conditional<IsSame<tmp_type,std::tuple<>>::value, 
                                       std::tuple<Number<0>>,
                                       tmp_type>::type;   
  static constexpr Integer number= Heaviside(IsTestOrTrial<Left>::number)+Heaviside(IsTestOrTrial<Right>::number);
  static_assert(TupleTypeSize<type>::value<2," In Binary<Left,Right>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
};








template<typename Left,typename Right>
class TypeOfForm;


template<Integer M,Integer N>
class TypeOfForm<Number<M>,Number<N>>
{
public:
  using type=void;
  static_assert(0*Number<N>::value,"L2inner: the form is neither a 0-form(function,function), 1-form(function/test,test/function) or 2-form (test/trial,trial/test), where the function is neither a test nor a trial");
};



template<>
class TypeOfForm<Number<0>,Number<0>>
{
  public:
    using type=Number<0>; 
};

template<>
class TypeOfForm<Number<0>,Number<1>>
{
  public:
    using type=Number<1>; 
};

template<>
class TypeOfForm<Number<1>,Number<0>>
{
  public:
    using type=Number<1>; 
};

template<>
class TypeOfForm<Number<1>,Number<2>>
{
  public:
    using type=Number<2>; 
};


template<>
class TypeOfForm<Number<2>,Number<1>>
{
  public:
    using type=Number<2>; 
};






// template<Integer N,typename...Args >
// auto MakeTest(const FunctionSpace<Args...>& W)
// {return Test<FunctionSpace<Args...>,N>(W);}


// template<Integer N,typename...Args >
// auto MakeTrial(const FunctionSpace<Args...>& W)
// {return Trial<FunctionSpace<Args...>,N>(W);}









// template<Integer N,typename...Args >
// auto MakeTest(const MixedSpace<Args...>& W)
// {return Test<MixedSpace<Args...>,N>(W);}


// template<Integer N,typename...Args >
// auto MakeTrial(const MixedSpace<Args...>& W)
// {return Trial<MixedSpace<Args...>,N>(W);}


// template<Integer N,typename...Args >
// auto MakeTest(const FullSpace<Args...>& W)
// {return Test<FullSpace<Args...>,N>(W);}


// template<Integer N,typename...Args >
// auto MakeTrial(const FullSpace<Args...>& W)
// {return Trial<FullSpace<Args...>,N>(W);}



// template<typename MixedSpace,Integer N>
// Trial<MixedSpace,N,GradientOperator> 
// Grad(const Trial<MixedSpace,N,IdentityOperator>& t)
// {return Trial<MixedSpace,N,GradientOperator> (t.spaces_ptr());}

// template<typename MixedSpace,Integer N>
// Test<MixedSpace,N,GradientOperator> 
// Grad(const Test<MixedSpace,N,IdentityOperator>& t)
// {return Test<MixedSpace,N,GradientOperator> (t.spaces_ptr());}

// template<typename MixedSpace,Integer N>
// Trial<MixedSpace,N,DivergenceOperator> 
// Div(const Trial<MixedSpace,N,IdentityOperator>& t)
// {return Trial<MixedSpace,N,DivergenceOperator> (t.spaces_ptr());}

// template<typename MixedSpace,Integer N>
// Test<MixedSpace,N,DivergenceOperator> 
// Div(const Test<MixedSpace,N,IdentityOperator>& t)
// {return Test<MixedSpace,N,DivergenceOperator> (t.spaces_ptr());}


// template<typename MixedSpace,Integer N>
// Trial<MixedSpace,N,CurlOperator> 
// Curl(const Trial<MixedSpace,N,IdentityOperator>& t)
// {return Trial<MixedSpace,N,CurlOperator> (t.spaces_ptr());}

// template<typename MixedSpace,Integer N>
// Test<MixedSpace,N,CurlOperator> 
// Curl(const Test<MixedSpace,N,IdentityOperator>& t)
// {return Test<MixedSpace,N,CurlOperator> (t.spaces_ptr());}



// template<typename MixedSpace,Integer N>
// Trial<MixedSpace,N,TraceOperator> 
// Trace(const Trial<MixedSpace,N,IdentityOperator>& t)
// {return Trial<MixedSpace,N,TraceOperator> (t.spaces_ptr());}

// template<typename MixedSpace,Integer N>
// Test<MixedSpace,N,TraceOperator> 
// Trace(const Test<MixedSpace,N,IdentityOperator>& t)
// {return Test<MixedSpace,N,TraceOperator> (t.spaces_ptr());}






template<Integer...N>
class FormTestTrialNumbers;

template<Integer LeftSpaceNumber,Integer RightSpaceNumber>
class FormTestTrialNumbers<2,2,1,LeftSpaceNumber,RightSpaceNumber>
{
public:
  using type=std::tuple<Number<RightSpaceNumber>,Number<LeftSpaceNumber>>;
};

template<Integer LeftSpaceNumber,Integer RightSpaceNumber>
class FormTestTrialNumbers<2,1,2,LeftSpaceNumber,RightSpaceNumber>
{
public:
  using type=std::tuple<Number<LeftSpaceNumber>,Number<RightSpaceNumber>>;
};

template<Integer LeftSpaceNumber,Integer RightSpaceNumber>
class FormTestTrialNumbers<1,0,1,LeftSpaceNumber,RightSpaceNumber>
{
public:
  using type=std::tuple<Number<RightSpaceNumber>>;
};

template<Integer LeftSpaceNumber,Integer RightSpaceNumber>
class FormTestTrialNumbers<1,1,0,LeftSpaceNumber,RightSpaceNumber>
{
public:
  using type=std::tuple<Number<LeftSpaceNumber>>;
};



















template<Integer H,typename T1>
class TupleOfTestTrialPairsNumbersAuxAuxAux;

template<Integer H,typename Left,typename Right,Integer QR, bool VolumeIntegral1>
class TupleOfTestTrialPairsNumbersAuxAuxAux<H,L2DotProductIntegral<Left,Right,VolumeIntegral1,QR>>
{
public:
  using type=EmptyExpression;
};


template<typename Left,typename Right,Integer QR, bool VolumeIntegral1>
class TupleOfTestTrialPairsNumbersAuxAuxAux<-1,L2DotProductIntegral<Left,Right,VolumeIntegral1,QR>>
{
public:
  using type=L2DotProductIntegral<Left,Right,VolumeIntegral1,QR>;
};

template<typename Left,typename Right,Integer QR>
class TupleOfTestTrialPairsNumbersAuxAuxAux<0,L2DotProductIntegral<Left,Right,true,QR>>
{
public:
  using type=L2DotProductIntegral<Left,Right,true,QR>;
};

template<typename Left,typename Right,Integer QR>
class TupleOfTestTrialPairsNumbersAuxAuxAux<1,L2DotProductIntegral<Left,Right,false,QR>>
{
public:
  using type=L2DotProductIntegral<Left,Right,false,QR>;
};




template<Integer H,typename T1,typename T2>
class TupleOfTestTrialPairsNumbersAuxAux2;


template<Integer H,typename T>
class TupleOfTestTrialPairsNumbersAuxAux2<H, T, Expression<EmptyExpression> >
{
 public:
  using type=EmptyExpression;
};



template<Integer H,typename T, typename Left,typename Right,Integer QR, bool VolumeIntegral1>
class TupleOfTestTrialPairsNumbersAuxAux2<H,T, Expression<L2DotProductIntegral<Left,Right,VolumeIntegral1,QR>> >
{
 public:
  using L2=L2DotProductIntegral<Left,Right,VolumeIntegral1,QR>;
  using tmptype=typename TupleOfTestTrialPairsNumbersAuxAuxAux<H,L2>::type;
  using type=typename std::conditional<IsSame<T,typename L2::TestTrialNumbers>::value, tmptype, EmptyExpression>::type;
};


template<Integer H,typename T1, typename T2>
class TupleOfTestTrialPairsNumbersAuxAux2<H, T1, Expression<Addition<Expression<T2>,Expression<EmptyExpression> > >>
{
 public:
  using type=typename TupleOfTestTrialPairsNumbersAuxAux2<H, T1,Expression<T2>>::type;
};


template<Integer H,typename T1, typename T2>
class TupleOfTestTrialPairsNumbersAuxAux2<H, T1, Expression<Addition<Expression<EmptyExpression>, Expression<T2> > >>
{
 public:
  using type=typename TupleOfTestTrialPairsNumbersAuxAux2<H, T1,Expression<T2>>::type;
};

template<Integer H,typename T>
class TupleOfTestTrialPairsNumbersAuxAux2<H, T, Expression<Addition<Expression<EmptyExpression>, Expression<EmptyExpression> > >>
{
 public:
  using type=EmptyExpression;
};




template<Integer H,typename T, typename Left,typename Right>
class TupleOfTestTrialPairsNumbersAuxAux2<H,T,Expression<Addition<Expression<Left>,
                                                               Expression<Right> >>>
{
 public:
  using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux2<H,T,Expression<Left>>::type>,
                      Expression<typename TupleOfTestTrialPairsNumbersAuxAux2<H,T,Expression<Right>>::type>>;
};



template<Integer H,typename...Ts>
class TupleOfTestTrialPairsNumbersAux2;

template<Integer H,typename T,typename Left,typename Right,Integer QR,bool VolumeIntegral>
class TupleOfTestTrialPairsNumbersAux2<H,T,L2DotProductIntegral<Left,Right,VolumeIntegral,QR> >
{
 public:
  using S=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux2<H,T,Expression<S>>::type;

};

template<Integer H,typename T, typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1,
                     typename Left2,typename Right2,Integer QR2,bool VolumeIntegral2>
class TupleOfTestTrialPairsNumbersAux2<H,T,Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,
                                                   Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>>>>
{
 public:
  using Left=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;
  using Right=L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux2<H,T,
                          Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux2<H,T,Left>::type>,
                                              Expression<typename TupleOfTestTrialPairsNumbersAux2<H,T,Right>::type>>>
                            >::type;
};

template<Integer H,typename T,typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1,typename Right>
class TupleOfTestTrialPairsNumbersAux2<H,T, Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,Expression<Right > > >
{
 public:
  using Left=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux2<H,T,
                          Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux2<H,T,Left>::type>,
                                                Expression<typename TupleOfTestTrialPairsNumbersAux2<H,T,Right>::type>>>
                             >::type;

};

template<Integer H, typename T, typename Left, typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1>
class TupleOfTestTrialPairsNumbersAux2<H, T, Addition<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1> > > >
{
 public:
  using Right=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux2<H,T,
                             Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux2<H,T,Left>::type>,
                                                 Expression<typename TupleOfTestTrialPairsNumbersAux2<H,T,Right>::type>>>
                             >::type;
};

template<Integer H, typename T, typename Left,typename Right>
class TupleOfTestTrialPairsNumbersAux2<H, T, Addition<Expression<Left>,Expression<Right > > >
{
public:

  using type=typename TupleOfTestTrialPairsNumbersAuxAux2<H,T,
                            Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux2<H, T,Left>::type>,
                                                Expression<typename TupleOfTestTrialPairsNumbersAux2<H, T,Right>::type>>>
                             >::type;
};





























template<typename T1,typename T2>
class TupleOfTestTrialPairsNumbersAuxAux;

template<typename T>
class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<EmptyExpression> >
{
 public:
  using type=EmptyExpression;
};


template<typename T, typename Left,typename Right,Integer QR, bool VolumeIntegral1>
class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<L2DotProductIntegral<Left,Right,VolumeIntegral1,QR>> >
{
 public:
  using L2=L2DotProductIntegral<Left,Right,VolumeIntegral1,QR>;
  using type=typename std::conditional<IsSame<T,typename L2::TestTrialNumbers>::value, L2, EmptyExpression>::type;
};




template<typename T1, typename T2>
class TupleOfTestTrialPairsNumbersAuxAux<T1, Expression<Addition<Expression<T2>,Expression<EmptyExpression> > >>
{
 public:
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T1,Expression<T2>>::type;
};


template<typename T1, typename T2>
class TupleOfTestTrialPairsNumbersAuxAux<T1, Expression<Addition<Expression<EmptyExpression>, Expression<T2> > >>
{
 public:
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T1,Expression<T2>>::type;
};

template<typename T>
class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<Addition<Expression<EmptyExpression>, Expression<EmptyExpression> > >>
{
 public:
  using type=EmptyExpression;
};



// template<typename T, typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1,
//                      typename Left2,typename Right2,Integer QR2,bool VolumeIntegral2>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,
//                                                                Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>>>>>
// {
//  public:
//   using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>>::type>,
//                        Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>>>::type>>;
// };


// template<typename T, typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1,
//                      typename Right2>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,
//                                                                Expression<Right2>>>>
// {
//  public:
//   using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>>::type>,
//                        Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Right2>>::type>>;
// };
// template<typename T, typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,
//                                                                Expression<EmptyExpression>>>>
// {
//  public:
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>>::type;
                       
// };



// template<typename T, typename Left1,
//                      typename Left2,typename Right2,Integer QR2,bool VolumeIntegral2>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<Left1>,
//                                                                Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>>>>>
// {
//  public:
//   using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Left1>>::type>,
//                        Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>>>::type>>;
// };

// template<typename T, typename Left2,typename Right2,Integer QR2,bool VolumeIntegral2>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<EmptyExpression>,
//                                                                Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>>>>>
// {
//  public:
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>>>::type;
// };


template<typename T, typename Left,typename Right>
class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<Left>,
                                                               Expression<Right> >>>
{
 public:
  using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Left>>::type>,
                      Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Right>>::type>>;
};



template<typename...Ts>
class TupleOfTestTrialPairsNumbersAux;

template<typename T,typename Left,typename Right,Integer QR,bool VolumeIntegral>
class TupleOfTestTrialPairsNumbersAux<T,L2DotProductIntegral<Left,Right,VolumeIntegral,QR> >
{
 public:
  using S=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<S>>::type;

};

template<typename T, typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1,
                     typename Left2,typename Right2,Integer QR2,bool VolumeIntegral2>
class TupleOfTestTrialPairsNumbersAux<T,Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,
                        Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>>>>
{
 public:
  using Left=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;
  using Right=L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
                          Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
                                                Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
                            >::type;
};

template<typename T,typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1,typename Right>
class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,Expression<Right > > >
{
 public:
  using Left=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
                          Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
                                                Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
                             >::type;

};

template<typename T, typename Left, typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1>
class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1> > > >
{
 public:
  using Right=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
                             Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
                                                 Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
                             >::type;
};

template<typename T, typename Left,typename Right>
class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<Left>,Expression<Right > > >
{
public:

  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
                            Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
                                                Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
                             >::type;
};


/////////////////////////////////////////////////////////////////////////////////////////////////
///// It build the tuple of the pairs test/trial used in the form.
///// Example: X=[P1,P0], u,v in P1, r,s in P0
///// L= int u * v - int u* s -int r * v
///// TupleOfTestTrialPairsNumbers= {(0,0),(0,1),(1,0)}
///// where (1,1) is missing because there is no int r * s
/////////////////////////////////////////////////////////////////////////////////////////////////
template<typename...Ts>
class TupleOfTestTrialPairsNumbers;

template<typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1>
class TupleOfTestTrialPairsNumbers<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>
{
 public:
  using T=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;

  using type=std::tuple<typename T::TestTrialNumbers>;

};

template<typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1,
         typename Left2,typename Right2,Integer QR2,bool VolumeIntegral2>
class TupleOfTestTrialPairsNumbers<Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,
                        Expression<L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>>>>
{
 public:
  using Left=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;
  using Right=L2DotProductIntegral<Left2,Right2,VolumeIntegral2,QR2>;

  using type=RemoveTupleDuplicates<std::tuple<typename Left::TestTrialNumbers, typename Right::TestTrialNumbers>>;

};

template<typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1,typename Right>
class TupleOfTestTrialPairsNumbers<Addition<Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>>,Expression<Right > > >
{
 public:
  using Left=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;
  using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type >>;
};

template<typename Left,typename Left1,typename Right1,Integer QR1,bool VolumeIntegral1>
class TupleOfTestTrialPairsNumbers<Addition<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1> > > >
{
 public:
  using Right=L2DotProductIntegral<Left1,Right1,VolumeIntegral1,QR1>;
  using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type >>;
};

template<typename Left,typename Right>
class TupleOfTestTrialPairsNumbers<Addition<Expression<Left>,Expression<Right > > >
{
public:
  using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type>>;
};



template<typename TupleOfPairsNumbers, typename Form,Integer Nmax,Integer N>
class TupleOfL2ProductsHelper;

template<typename TupleOfPairsNumbers, typename Form,Integer Nmax>
class TupleOfL2ProductsHelper<TupleOfPairsNumbers,Form,Nmax,Nmax>
{
 public:
 using type=std::tuple<typename TupleOfTestTrialPairsNumbersAux<GetType<TupleOfPairsNumbers,Nmax>,Form>::type>;
};



template<typename TupleOfPairsNumbers, typename Form,Integer Nmax,Integer N>
class TupleOfL2ProductsHelper
{
 public:
 using type=TupleCatType<std::tuple<typename TupleOfTestTrialPairsNumbersAux<GetType<TupleOfPairsNumbers,N>,Form>::type> , 
                         typename TupleOfL2ProductsHelper<TupleOfPairsNumbers,Form,Nmax,N+1>::type >;
};

template<typename TupleOfPairsNumbers, typename Form>
class TupleOfL2Products
{
public:
using type=typename TupleOfL2ProductsHelper<TupleOfPairsNumbers,Form,TupleTypeSize<TupleOfPairsNumbers>::value-1,0>::type;


};


















template<Integer H, typename TupleOfPairsNumbers, typename Form,Integer Nmax,Integer N>
class TupleOfL2ProductsHelper2;

template<Integer H, typename TupleOfPairsNumbers, typename Form,Integer Nmax>
class TupleOfL2ProductsHelper2<H,TupleOfPairsNumbers,Form,Nmax,Nmax>
{
 public:
 using type=std::tuple<typename TupleOfTestTrialPairsNumbersAux2<H,GetType<TupleOfPairsNumbers,Nmax>,Form>::type>;
};



template<Integer H, typename TupleOfPairsNumbers, typename Form,Integer Nmax,Integer N>
class TupleOfL2ProductsHelper2
{
 public:
 using type=TupleCatType<std::tuple<typename TupleOfTestTrialPairsNumbersAux2<H,GetType<TupleOfPairsNumbers,N>,Form>::type> , 
                         typename TupleOfL2ProductsHelper2<H,TupleOfPairsNumbers,Form,Nmax,N+1>::type >;
};

template<Integer H, typename TupleOfPairsNumbers, typename Form>
class TupleOfL2Products2
{
public:
using type=typename TupleOfL2ProductsHelper2<H,TupleOfPairsNumbers,Form,TupleTypeSize<TupleOfPairsNumbers>::value-1,0>::type;


};




template<typename FunctionSpace, typename T>
std::shared_ptr<FunctionSpace> find_spaces_ptr(const T& t)
{
  std::cout<<"is null"<<std::endl;
 return std::shared_ptr<FunctionSpace>(nullptr);
}


template<typename FunctionSpace, template<class,Integer,class>class TestOrTrial_,typename MixedSpace,Integer N,typename OperatorType>
std::shared_ptr<FunctionSpace> find_spaces_ptr(const TestOrTrial_<MixedSpace,N,OperatorType>& t)
{
  std::cout<<"is not null"<<std::endl;
 return t.spaces_ptr();
}


template<typename FunctionSpace, typename Left, typename Right, bool VolumeIntegral, Integer QR>
std::shared_ptr<FunctionSpace> find_spaces_ptr(const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& t)
{
 return t.spaces_ptr();
}


template<typename FunctionSpace, template<class>class Unary, typename T>
std::shared_ptr<FunctionSpace> find_spaces_ptr(const Unary<Expression<T>>& t)
{
 return find_spaces_ptr<FunctionSpace>(t.derived());
}


// std::enable_if_t<>
// find_spaces_ptr_aux(const )
// {

// }
template<typename FunctionSpace, template<class,class>class Binary, typename Left,typename Right>
std::shared_ptr<FunctionSpace> find_spaces_ptr(const Binary<Expression<Left>,Expression<Right>>& t)
{
 auto p1=find_spaces_ptr<FunctionSpace>(t.left());
 auto p2=find_spaces_ptr<FunctionSpace>(t.right());
 return (p1==nullptr && p2==nullptr) ? nullptr : ((p1 != nullptr) ? p1 : p2);
 // if(tmp==NULL)
 //  return find_spaces_ptr(t.right());
 // else
 //  return tmp;
}


}
#endif