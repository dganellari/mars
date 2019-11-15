#ifndef MARS_OPERATOR_TUPLE_TYPE_HPP
#define MARS_OPERATOR_TUPLE_TYPE_HPP

#include "mars_tuple_utilities.hpp"
#include "mars_base.hpp"
#include "mars_operators.hpp"
#include "mars_operator_differential.hpp"
namespace mars{


template<typename...Ts>
class CompositeOperator;

template<typename ConstType,typename...Inputs>
class ConstantTensor;

template<typename FullSpace,Integer N,typename Operator_,typename FuncType>
class Function;

template<typename MixedSpace, Integer N, typename OperatorType>
class Test;

template<typename MixedSpace, Integer N, typename OperatorType>
class Trial;

class TraceOperator;

// H=-1: both volume and surface integrals
// H=0: only volume integral
// H=1: only surface integral
template<typename T,Integer H=-1>
class OperatorAndQuadratureTupleType;


template<template<class>class Unary, typename T,Integer H>
    class OperatorAndQuadratureTupleType< Unary< Expression<T> >, H >
    { public:
      using type=typename OperatorAndQuadratureTupleType<T,H>::type;
      using composite_type=typename OperatorAndQuadratureTupleType<T,H>::composite_type;
  };


template<template<class,class>class Binary,typename Left, typename Right,Integer H>
  class OperatorAndQuadratureTupleType< Binary< Expression<Left>,Expression<Right> >, H >
  { public:
      using LeftT=typename OperatorAndQuadratureTupleType<Left,H>::type;
      using RightT=typename OperatorAndQuadratureTupleType<Right,H>::type;
      using type=RemoveTupleOfTupleDuplicates< LeftT, RightT >;
      using CompositeLeftT=typename OperatorAndQuadratureTupleType<Left,H>::composite_type;
      using CompositeRightT=typename OperatorAndQuadratureTupleType<Right,H>::composite_type;
      using composite_type=RemoveTupleOfTupleDuplicates<  CompositeLeftT, CompositeRightT >;
  };

template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename OperatorType,Integer H>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,OperatorType>,H >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,OperatorType>;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,OperatorType>>::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,OperatorType>>::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");

      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using single_type=std::tuple<std::tuple< OperatorType,std::tuple<> >>;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using type=TupleChangeType<TestOrTrial::value,single_type,emptytuple>;
      using composite_type=emptytuple;
  };


template<typename FuncType, typename MixedSpace, Integer N,typename OperatorType,Integer H>
  class OperatorAndQuadratureTupleType<Function<MixedSpace,N,OperatorType,FuncType>,H >
  { public:
      using Func=Function<MixedSpace,N,OperatorType,FuncType>;
      static constexpr Integer Nmax= Func::Nmax;
      using single_type=std::tuple<std::tuple< OperatorType,std::tuple<> >>;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using type=TupleChangeType<Func::value,single_type,emptytuple>;
      using composite_type=emptytuple;
  };

template<typename ConstType, typename...Inputs,Integer H>
  class OperatorAndQuadratureTupleType<ConstantTensor<ConstType,Inputs...>,H >
  { public:
      using type=std::tuple<std::tuple<>>;
      using composite_type=std::tuple<std::tuple<>>;
  };


template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< TraceOperator > > >,0 >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< TraceOperator > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator< Expression<TraceOperator > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator< Expression<TraceOperator > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using type=TupleOfType<Nmax,std::tuple<> > ;
      using composite_type=type;

  };


template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< TraceOperator > > >,1 >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< TraceOperator > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator< Expression<TraceOperator > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator< Expression<TraceOperator > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using single_type=std::tuple<std::tuple< TraceOperator,std::tuple<> >>;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using type=TupleChangeType<TestOrTrial::value,single_type,emptytuple>;
      using single_composite_type=std::tuple<std::tuple< CompositeOperator<Expression< TraceOperator > >,std::tuple<> >>;
      using composite_type=TupleChangeType<TestOrTrial::value,single_composite_type,emptytuple>;

  };

template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename Expr,Integer H>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< Expr > > >,H >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< Expr > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator< Expression<Expr > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator< Expression<Expr > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using single_type=std::tuple<std::tuple< Expr,std::tuple<> >>;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using type=TupleChangeType<TestOrTrial::value,single_type,emptytuple>;
      using single_composite_type=std::tuple<std::tuple< CompositeOperator<Expression< Expr > >,std::tuple<> >>;
      using composite_type=TupleChangeType<TestOrTrial::value,single_composite_type,emptytuple>;

  };


template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename ConstType,typename...Inputs,Integer H>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< ConstantTensor<ConstType,Inputs...> > > >,H >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< ConstantTensor<ConstType,Inputs...> > > >;
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using type=std::tuple<std::tuple<>>;
      using composite_type=emptytuple;

  };



template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,
         typename FullSpace, Integer M,typename Operator_,typename FuncType,Integer H>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< Function<FullSpace,M,Operator_,FuncType> > > >,H >
  { public:

      using Func=Function<FullSpace,M,Operator_,FuncType>;
      static constexpr Integer Nmax= Func::Nmax;
      using single_type=std::tuple<std::tuple< Operator_,std::tuple<> >>;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using type=TupleChangeType<Func::value,single_type,emptytuple>;
      using composite_type=emptytuple;
  };


template<template<class> class UnaryOperator,template<class,Integer,class > class TestOrTrial_,
         typename MixedSpace, Integer N,typename Type,Integer H>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<UnaryOperator< Expression<Type>> > > >, H >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< UnaryOperator< Expression<Type>> > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator< Expression<UnaryOperator< Expression<Type>> > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator< Expression<UnaryOperator< Expression<Type>> > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");

      using TestOrTrialType=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Type>>>;
      using type=typename OperatorAndQuadratureTupleType<TestOrTrialType,H>::type;
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using single_composite_type=std::tuple<std::tuple< CompositeOperator<Expression< UnaryOperator< Expression<Type>> > >,std::tuple<> >>;
      using composite_type=TupleChangeType<TestOrTrial::value,single_composite_type,emptytuple>;

  };

template<template<class,class>class Binary, template<class,Integer,class > class TestOrTrial_,
         typename MixedSpace, Integer N,typename Left,typename Right,Integer H>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Binary< Expression<Left>,Expression<Right> > > > >, H >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< Binary<Expression<Left>,Expression<Right>> > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator< Expression<Binary< Expression<Left>,Expression<Right> > > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator< Expression<Binary< Expression<Left>,Expression<Right> > > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");

      using TestOrTrialLeft=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Left>>>;
      using TestOrTrialRight=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Right>>>;
      using LeftT=typename OperatorAndQuadratureTupleType<TestOrTrialLeft,H>::type;
      using RightT=typename OperatorAndQuadratureTupleType<TestOrTrialRight,H>::type;
      using type=RemoveTupleOfTupleDuplicates< LeftT, RightT >;
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using single_composite_type=std::tuple<std::tuple< CompositeOperator<Expression< Binary< Expression<Left>,Expression<Right> > > >,std::tuple<> >>;
      using composite_type=TupleChangeType<TestOrTrial::value,single_composite_type,emptytuple>;

  };


template<typename Left,typename Right, bool VolumeIntegral,Integer QR, Integer H>
class OperatorAndQuadratureTupleType<L2DotProductIntegral<Left,Right,VolumeIntegral,QR>,H>
{ 
public:
  using L2prod=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
  using QRule=typename L2prod::QRule;
  using Operatortuple=typename OperatorAndQuadratureTupleType<typename L2prod::type,H>::type;
  using type=TupleOfTupleChangeType<1,QRule, Operatortuple>;
  using CompositeOperatortuple=typename OperatorAndQuadratureTupleType<typename L2prod::type,H>::composite_type;
  using composite_type=TupleOfTupleChangeType<1,QRule, CompositeOperatortuple>;
} 
; 


template<typename Left,typename Right,Integer QR>
class OperatorAndQuadratureTupleType<L2DotProductIntegral<Left,Right,false,QR>,0>
{ 
public:
  using L2prod=L2DotProductIntegral<Left,Right,false,QR>;
  using UniqueElementFunctionSpacesTupleType=typename L2prod::UniqueElementFunctionSpacesTupleType;
  static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
  using type=TupleOfType<Nmax,std::tuple<> > ;
  using composite_type=type;  
} 
; 





template<typename Left,typename Right,Integer QR>
class OperatorAndQuadratureTupleType<L2DotProductIntegral<Left,Right,true,QR>,1>
{ 
public:
  using L2prod=L2DotProductIntegral<Left,Right,true,QR>;
  using UniqueElementFunctionSpacesTupleType=typename L2prod::UniqueElementFunctionSpacesTupleType;
  static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
  using type=TupleOfType<Nmax,std::tuple<> > ;
  using composite_type=type;  
} 
; 


// template<typename T>
//   class OperatorAndQuadratureTupleType< Transposed< Expression<Transposed<T> >>>
//   { public:

//       using type=typename OperatorAndQuadratureTupleType<T>::type;
//       using composite_type=typename OperatorAndQuadratureTupleType<T>::composite_type;
//   };

template<typename Tuple,typename DirichletBc,typename...DirichletBCs>
class DirichletBCMapsHelper;

template<typename Tuple,typename DirichletBC>
class DirichletBCMapsHelper<Tuple,DirichletBC>
{
public:
 static constexpr Integer map_value=DirichletBC::map_value;
 static constexpr Integer value=DirichletBC::value;
 using FunctionSpace=typename DirichletBC::FunctionSpace;
 using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
 using Elem=typename FunctionSpace::Elem;
 using TraceElem=TraceOf<Elem>; 
 static constexpr Integer FEFamily=GetType<UniqueElementFunctionSpacesTupleType,value>::FEFamily;
 using tmp_type=MapFromReference<TraceOperator,TraceElem, FEFamily>;
 using type=TupleChangeType<map_value,tmp_type,Tuple>;
};

//ShapeFunction<Elem,BaseFunctionSpace<-10,1, 1, 1>,TraceOperator,GaussPoints<TraceElem, 1> > >

template<typename Tuple,typename DirichletBC,typename...DirichletBCs>
class DirichletBCMapsHelper
{
public:
 static constexpr Integer map_value=DirichletBC::map_value;
 static constexpr Integer value=DirichletBC::value;
 using FunctionSpace=typename DirichletBC::FunctionSpace;
 using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
 using Elem=typename FunctionSpace::Elem;
 using TraceElem=TraceOf<Elem>; 
 static constexpr Integer FEFamily=GetType<UniqueElementFunctionSpacesTupleType,value>::FEFamily;
 using tmp_type=MapFromReference<TraceOperator,TraceElem,FEFamily>;
 using tmp_tuple=TupleChangeType<map_value,tmp_type,Tuple>;
 using type=typename DirichletBCMapsHelper<tmp_tuple,DirichletBCs...>::type;
};

template<typename DirichletBC, typename...DirichletBCs>
class DirichletBCMapsAux
{
public:
 using FunctionSpace=typename DirichletBC::FunctionSpace;
 using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
 using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;
 static constexpr Integer Nmax=MaxNumberInTuple<SpacesToUniqueFEFamilies>::value+1;
 using emptytuple=TupleOfTypeTCreate<std::tuple<>,Nmax>;
 using type=typename DirichletBCMapsHelper<emptytuple,DirichletBC,DirichletBCs...>::type;
};




template<typename...DirichletBCs>
using DirichletBCMapsCollection=typename DirichletBCMapsAux<DirichletBCs...>::type;









}
#endif