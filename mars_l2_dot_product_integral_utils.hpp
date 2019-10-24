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
class IsVolumeOrSurfaceIntegral
{
public:
  static constexpr bool surface=true;
  static constexpr bool volume=true;
};

template<template<class,Integer,class>class TestOrTrial, typename MixedSpace,Integer N>
class IsVolumeOrSurfaceIntegral<TestOrTrial<MixedSpace,N,TraceOperator>>{
public:
  static constexpr bool surface=true;
  static constexpr bool volume=false;
};

template<template<class,Integer,class>class TestOrTrial, typename MixedSpace,Integer N, typename OperatorType>
class IsVolumeOrSurfaceIntegral<TestOrTrial<MixedSpace,N,OperatorType>>{
public:
  static constexpr bool surface=false;
  static constexpr bool volume=true;
};

template<template<class>class Unary, typename Type>
class IsVolumeOrSurfaceIntegral< Unary<Expression<Type>> >
{
public:
  static constexpr bool surface=IsVolumeOrSurfaceIntegral<Type>::surface;
  static constexpr bool volume=IsVolumeOrSurfaceIntegral<Type>::volume;
};


template<template<class,class>class Binary, typename Left,typename Right>
class IsVolumeOrSurfaceIntegral< Binary<Expression<Left>,Expression<Right>> >
{
public:
  static constexpr bool surface=IsVolumeOrSurfaceIntegral<Left>::surface*IsVolumeOrSurfaceIntegral<Right>::surface;
  static constexpr bool volume=IsVolumeOrSurfaceIntegral<Left>::volume*IsVolumeOrSurfaceIntegral<Right>::volume;
};


template<typename Left,typename Right>
class IsVolumeOrSurfaceIntegral<InnerProduct<Expression <Left>, Expression <Right> >  >
{
public:
  static constexpr bool surface=IsVolumeOrSurfaceIntegral<Left>::surface*IsVolumeOrSurfaceIntegral<Right>::surface;
  static constexpr bool volume=IsVolumeOrSurfaceIntegral<Left>::volume*IsVolumeOrSurfaceIntegral<Right>::volume;
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























template<typename T1,typename T2>
class TupleOfTestTrialPairsNumbersAuxAux;

template<typename T>
class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<std::tuple<>> >
{
 public:
  using type=std::tuple<>;
};



// template<typename T, typename MeshT, typename Left,typename Right,Integer QR>
// class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<L2DotProductIntegral<MeshT,Left,Right,QR>> >
// {
//  public:
//   using L2=L2DotProductIntegral<MeshT,Left,Right,QR>;
//   using type=typename std::conditional<IsSame<T,typename L2::TestTrialNumbers>::value, L2, std::tuple<>>::type;
// };
template<typename T, typename Left,typename Right,Integer QR>
class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<L2DotProductIntegral<Left,Right,QR>> >
{
 public:
  using L2=L2DotProductIntegral<Left,Right,QR>;
  using type=typename std::conditional<IsSame<T,typename L2::TestTrialNumbers>::value, L2, std::tuple<>>::type;
};




template<typename T1, typename T2>
class TupleOfTestTrialPairsNumbersAuxAux<T1, Expression<Addition<Expression<T2>,Expression<std::tuple<>> > >>
{
 public:
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T1,Expression<T2>>::type;
};


template<typename T1, typename T2>
class TupleOfTestTrialPairsNumbersAuxAux<T1, Expression<Addition<Expression<std::tuple<>>, Expression<T2> > >>
{
 public:
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T1,Expression<T2>>::type;
};

template<typename T>
class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<Addition<Expression<std::tuple<>>, Expression<std::tuple<>> > >>
{
 public:
  using type=std::tuple<>;
};


// template<typename T, typename MeshT1, typename Left1,typename Right1,Integer QR1,
//                      typename MeshT2, typename Left2,typename Right2,Integer QR2>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<MeshT1,Left1,Right1,QR1>>,
//                                                                  Expression<L2DotProductIntegral<MeshT2,Left2,Right2,QR2>>>>>
// {
//  public:
//   using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<MeshT1,Left1,Right1,QR1>>>::type>,
//                        Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<MeshT2,Left2,Right2,QR2>>>::type>>;
// };
template<typename T, typename Left1,typename Right1,Integer QR1,
                     typename Left2,typename Right2,Integer QR2>
class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
                                                               Expression<L2DotProductIntegral<Left2,Right2,QR2>>>>>
{
 public:
  using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left1,Right1,QR1>>>::type>,
                       Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left2,Right2,QR2>>>::type>>;
};


template<typename T, typename Left1,typename Right1,Integer QR1,
                     typename Right2>
class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
                                                               Expression<Right2>>>>
{
 public:
  using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left1,Right1,QR1>>>::type>,
                       Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Right2>>::type>>;
};
template<typename T, typename Left1,typename Right1,Integer QR1>
class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
                                                               Expression<std::tuple<>>>>>
{
 public:
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left1,Right1,QR1>>>::type;
                       
};



template<typename T, typename Left1,
                     typename Left2,typename Right2,Integer QR2>
class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<Left1>,
                                                               Expression<L2DotProductIntegral<Left2,Right2,QR2>>>>>
{
 public:
  using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Left1>>::type>,
                       Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left2,Right2,QR2>>>::type>>;
};

template<typename T, typename Left2,typename Right2,Integer QR2>
class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<std::tuple<>>,
                                                               Expression<L2DotProductIntegral<Left2,Right2,QR2>>>>>
{
 public:
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left2,Right2,QR2>>>::type;
};


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





// template<typename T,typename MeshT, typename Left,typename Right,Integer QR>
// class TupleOfTestTrialPairsNumbersAux<T,L2DotProductIntegral<MeshT,Left,Right,QR> >
// {
//  public:
//   using S=L2DotProductIntegral<MeshT,Left,Right,QR>;
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<S>>::type;

// };
template<typename T,typename Left,typename Right,Integer QR>
class TupleOfTestTrialPairsNumbersAux<T,L2DotProductIntegral<Left,Right,QR> >
{
 public:
  using S=L2DotProductIntegral<Left,Right,QR>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<S>>::type;

};

// template<typename T, typename MeshT1, typename Left1,typename Right1,Integer QR1,
//                      typename MeshT2, typename Left2,typename Right2,Integer QR2>
// class TupleOfTestTrialPairsNumbersAux<T,Addition<Expression<L2DotProductIntegral<MeshT1,Left1,Right1,QR1>>,
//                         Expression<L2DotProductIntegral<MeshT2,Left2,Right2,QR2>>>>
// {
//  public:
//   using Left=L2DotProductIntegral<MeshT1,Left1,Right1,QR1>;
//   using Right=L2DotProductIntegral<MeshT2,Left2,Right2,QR2>;
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
//                           Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
//                                                 Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
//                             >::type;
// };
template<typename T, typename Left1,typename Right1,Integer QR1,
                     typename Left2,typename Right2,Integer QR2>
class TupleOfTestTrialPairsNumbersAux<T,Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
                        Expression<L2DotProductIntegral<Left2,Right2,QR2>>>>
{
 public:
  using Left=L2DotProductIntegral<Left1,Right1,QR1>;
  using Right=L2DotProductIntegral<Left2,Right2,QR2>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
                          Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
                                                Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
                            >::type;
};


// template<typename T,typename MeshT, typename Left1,typename Right1,Integer QR,typename Right>
// class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR>>,Expression<Right > > >
// {
//  public:
//   using Left=L2DotProductIntegral<MeshT,Left1,Right1,QR>;
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
//                           Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
//                                                 Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
//                              >::type;

// };

template<typename T,typename Left1,typename Right1,Integer QR,typename Right>
class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<L2DotProductIntegral<Left1,Right1,QR>>,Expression<Right > > >
{
 public:
  using Left=L2DotProductIntegral<Left1,Right1,QR>;
  using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
                          Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
                                                Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
                             >::type;

};

// template<typename T, typename Left,typename MeshT, typename Left1,typename Right1,Integer QR>
// class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<Left>,Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR> > > >
// {
//  public:
//   using Right=L2DotProductIntegral<MeshT,Left1,Right1,QR>;
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
//                              Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
//                                                    Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
//                              >::type;
// };
template<typename T, typename Left, typename Left1,typename Right1,Integer QR>
class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1,QR> > > >
{
 public:
  using Right=L2DotProductIntegral<Left1,Right1,QR>;
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



template<typename...Ts>
class TupleOfTestTrialPairsNumbers;

// template<typename MeshT1, typename Left1,typename Right1,Integer QR1>
// class TupleOfTestTrialPairsNumbers<L2DotProductIntegral<MeshT1,Left1,Right1,QR1>>
// {
//  public:
//   using T=L2DotProductIntegral<MeshT1,Left1,Right1,QR1>;

//   using type=std::tuple<typename T::TestTrialNumbers>;

// };

template<typename Left1,typename Right1,Integer QR1>
class TupleOfTestTrialPairsNumbers<L2DotProductIntegral<Left1,Right1,QR1>>
{
 public:
  using T=L2DotProductIntegral<Left1,Right1,QR1>;

  using type=std::tuple<typename T::TestTrialNumbers>;

};

// template<typename MeshT1, typename Left1,typename Right1,Integer QR1,
//          typename MeshT2, typename Left2,typename Right2,Integer QR2>
// class TupleOfTestTrialPairsNumbers<Addition<Expression<L2DotProductIntegral<MeshT1,Left1,Right1,QR1>>,
//                         Expression<L2DotProductIntegral<MeshT2,Left2,Right2,QR2>>>>
// {
//  public:
//   using Left=L2DotProductIntegral<MeshT1,Left1,Right1,QR1>;
//   using Right=L2DotProductIntegral<MeshT2,Left2,Right2,QR2>;

//   using type=RemoveTupleDuplicates<std::tuple<typename Left::TestTrialNumbers, typename Right::TestTrialNumbers>>;

// };

template<typename Left1,typename Right1,Integer QR1,
         typename Left2,typename Right2,Integer QR2>
class TupleOfTestTrialPairsNumbers<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
                        Expression<L2DotProductIntegral<Left2,Right2,QR2>>>>
{
 public:
  using Left=L2DotProductIntegral<Left1,Right1,QR1>;
  using Right=L2DotProductIntegral<Left2,Right2,QR2>;

  using type=RemoveTupleDuplicates<std::tuple<typename Left::TestTrialNumbers, typename Right::TestTrialNumbers>>;

};

// template<typename MeshT, typename Left1,typename Right1,Integer QR,typename Right>
// class TupleOfTestTrialPairsNumbers<Addition<Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR>>,Expression<Right > > >
// {
//  public:
//   using Left=L2DotProductIntegral<MeshT,Left1,Right1,QR>;
//   using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type >>;
// };

template<typename Left1,typename Right1,Integer QR,typename Right>
class TupleOfTestTrialPairsNumbers<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR>>,Expression<Right > > >
{
 public:
  using Left=L2DotProductIntegral<Left1,Right1,QR>;
  using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type >>;
};

// template<typename Left,typename MeshT, typename Left1,typename Right1,Integer QR>
// class TupleOfTestTrialPairsNumbers<Addition<Expression<Left>,Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR> > > >
// {
//  public:
//   using Right=L2DotProductIntegral<MeshT,Left1,Right1,QR>;
//   using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type >>;
// };

template<typename Left,typename Left1,typename Right1,Integer QR>
class TupleOfTestTrialPairsNumbers<Addition<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1,QR> > > >
{
 public:
  using Right=L2DotProductIntegral<Left1,Right1,QR>;
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
using type=  typename TupleOfL2ProductsHelper<TupleOfPairsNumbers,Form, TupleTypeSize<TupleOfPairsNumbers>::value-1,0>::type;


};





}
#endif