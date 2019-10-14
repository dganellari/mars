#ifndef MARS_OPERATOR_TUPLE_TYPE_HPP
#define MARS_OPERATOR_TUPLE_TYPE_HPP

#include "mars_tuple_utilities.hpp"
#include "mars_base.hpp"


namespace mars{

class IdentityOperator;

class GradientOperator;

class DivergenceOperator;

class CurlOperator;

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


template<typename Left,typename Right,Integer QR>
class L2DotProductIntegral;


template<typename...Parameters>
class UnaryPlus;

template<typename...Parameters>
class UnaryMinus;

template<typename...Parameters>
class Multiplication;

template<typename...Parameters>
class Contraction2;

template<typename...Parameters>
class Division;

template<typename...Parameters>
class Subtraction;

template<typename...Parameters>
class Addition;

template<typename T>
class Transposed;


template<typename...T>
    class OperatorAndQuadratureTupleType;


// order(T) = order(+T)
template<typename T>
    class OperatorAndQuadratureTupleType< UnaryPlus< Expression<T> > >
    { public:
      using type=typename OperatorAndQuadratureTupleType<T>::type;
      using composite_type=typename OperatorAndQuadratureTupleType<T>::composite_type;
  };


// order(T) = order(-T)
template<typename T>
  class OperatorAndQuadratureTupleType< UnaryMinus< Expression<T> > >
  { public:
      using type=typename OperatorAndQuadratureTupleType<T>::type;
      using composite_type=typename OperatorAndQuadratureTupleType<T>::composite_type;
  };



// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
  class OperatorAndQuadratureTupleType< Addition< Expression<Left>,Expression<Right> > >
  { public:
      using LeftT=typename OperatorAndQuadratureTupleType<Left>::type;
      using RightT=typename OperatorAndQuadratureTupleType<Right>::type;
      using type=RemoveTupleOfTupleDuplicates< LeftT, RightT >;
      using CompositeLeftT=typename OperatorAndQuadratureTupleType<Left>::composite_type;
      using CompositeRightT=typename OperatorAndQuadratureTupleType<Right>::composite_type;
      using composite_type=RemoveTupleOfTupleDuplicates<  CompositeLeftT, CompositeRightT >;
  };


// order(A+B) =max(order(A),order(B))
template<typename Left, typename Right>
  class OperatorAndQuadratureTupleType< Subtraction< Expression<Left>,Expression<Right> > >
  { public:
      using LeftT=typename OperatorAndQuadratureTupleType<Left>::type;
      using RightT=typename OperatorAndQuadratureTupleType<Right>::type;
      using type=RemoveTupleOfTupleDuplicates<  LeftT, RightT >;
      using CompositeLeftT=typename OperatorAndQuadratureTupleType<Left>::composite_type;
      using CompositeRightT=typename OperatorAndQuadratureTupleType<Right>::composite_type;
      using composite_type=RemoveTupleOfTupleDuplicates<  CompositeLeftT, CompositeRightT >;
  };



template<typename Left, typename Right>
  class OperatorAndQuadratureTupleType< Multiplication< Expression<Left>, Expression<Right> > >
  { public:
      using LeftT=typename OperatorAndQuadratureTupleType<Left>::type;
      using RightT=typename OperatorAndQuadratureTupleType<Right>::type;
      using type=RemoveTupleOfTupleDuplicates<  LeftT, RightT >;
      using CompositeLeftT=typename OperatorAndQuadratureTupleType<Left>::composite_type;
      using CompositeRightT=typename OperatorAndQuadratureTupleType<Right>::composite_type;
      using composite_type=RemoveTupleOfTupleDuplicates<  CompositeLeftT, CompositeRightT >;
  };

template<typename Left, typename Right>
  class OperatorAndQuadratureTupleType< Division< Expression<Left>, Expression<Right> > >
  { public:
      using LeftT=typename OperatorAndQuadratureTupleType<Left>::type;
      using RightT=typename OperatorAndQuadratureTupleType<Right>::type;
      using type=RemoveTupleOfTupleDuplicates<  LeftT, RightT >;
      using CompositeLeftT=typename OperatorAndQuadratureTupleType<Left>::composite_type;
      using CompositeRightT=typename OperatorAndQuadratureTupleType<Right>::composite_type;
      using composite_type=RemoveTupleOfTupleDuplicates<  CompositeLeftT, CompositeRightT >;
  };


// order(Left*Right) = order(Left) * order(Right)
template<typename Left, typename Right>
  class OperatorAndQuadratureTupleType< Contraction2< Expression<Left>, Expression<Right> > >
  { public:
      using LeftT=typename OperatorAndQuadratureTupleType<Left>::type;
      using RightT=typename OperatorAndQuadratureTupleType<Right>::type;
      using type=RemoveTupleOfTupleDuplicates<  LeftT, RightT >;
      using CompositeLeftT=typename OperatorAndQuadratureTupleType<Left>::composite_type;
      using CompositeRightT=typename OperatorAndQuadratureTupleType<Right>::composite_type;
      using composite_type=RemoveTupleOfTupleDuplicates<  CompositeLeftT, CompositeRightT >;
  };

template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename OperatorType>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,OperatorType> >
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



template<typename FuncType, typename MixedSpace, Integer N,typename OperatorType>
  class OperatorAndQuadratureTupleType<Function<MixedSpace,N,OperatorType,FuncType> >
  { public:
      using Func=Function<MixedSpace,N,OperatorType,FuncType>;
      static constexpr Integer Nmax= Func::Nmax;
      using single_type=std::tuple<std::tuple< OperatorType,std::tuple<> >>;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using type=TupleChangeType<Func::value,single_type,emptytuple>;
      using composite_type=emptytuple;
  };


template<typename ConstType, typename...Inputs>
  class OperatorAndQuadratureTupleType<ConstantTensor<ConstType,Inputs...> >
  { public:
      // using type=TupleOfType<3,std::tuple<> > ;
      // using composite_type=TupleOfType<3,std::tuple<> > ;
      using type=std::tuple<std::tuple<>>;
      using composite_type=std::tuple<std::tuple<>>;
  };


template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename Expr>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< Expr > > > >
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


template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename ConstType,typename...Inputs>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< ConstantTensor<ConstType,Inputs...> > > > >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< ConstantTensor<ConstType,Inputs...> > > >;
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using type=std::tuple<std::tuple<>>;
      using composite_type=emptytuple;

  };



template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,
         typename FullSpace, Integer M,typename Operator_,typename FuncType>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< Function<FullSpace,M,Operator_,FuncType> > > > >
  { public:

      using Func=Function<FullSpace,M,Operator_,FuncType>;
      static constexpr Integer Nmax= Func::Nmax;
      using single_type=std::tuple<std::tuple< Operator_,std::tuple<> >>;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using type=TupleChangeType<Func::value,single_type,emptytuple>;
      using composite_type=emptytuple;
  };

template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename Type>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<UnaryPlus< Expression<Type>> > > > >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< UnaryPlus< Expression<Type>> > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator< Expression<UnaryPlus< Expression<Type>> > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator< Expression<UnaryPlus< Expression<Type>> > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");

      using TestOrTrialType=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Type>>>;
      using type=typename OperatorAndQuadratureTupleType<TestOrTrialType>::type;
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using single_composite_type=std::tuple<std::tuple< CompositeOperator<Expression< UnaryPlus< Expression<Type>> > >,std::tuple<> >>;
      using composite_type=TupleChangeType<TestOrTrial::value,single_composite_type,emptytuple>;

  };

template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename Type>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<UnaryMinus< Expression<Type>> > > > >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< UnaryMinus< Expression<Type>> > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator< Expression<UnaryMinus< Expression<Type>> > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator< Expression<UnaryMinus< Expression<Type>> > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");

      using TestOrTrialType=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Type>>>;
      using type=typename OperatorAndQuadratureTupleType<TestOrTrialType>::type;
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using single_composite_type=std::tuple<std::tuple< CompositeOperator<Expression< UnaryMinus< Expression<Type>> > >,std::tuple<> >>;
      using composite_type=TupleChangeType<TestOrTrial::value,single_composite_type,emptytuple>;
  };

template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename Left,typename Right>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Addition< Expression<Left>,Expression<Right> > > > > >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< Addition<Expression<Left>,Expression<Right>> > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator< Expression<Addition< Expression<Left>,Expression<Right> > > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator< Expression<Addition< Expression<Left>,Expression<Right> > > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");

      using TestOrTrialLeft=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Left>>>;
      using TestOrTrialRight=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Right>>>;
      using LeftT=typename OperatorAndQuadratureTupleType<TestOrTrialLeft>::type;
      using RightT=typename OperatorAndQuadratureTupleType<TestOrTrialRight>::type;
      using type=RemoveTupleOfTupleDuplicates< LeftT, RightT >;
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using single_composite_type=std::tuple<std::tuple< CompositeOperator<Expression< Addition< Expression<Left>,Expression<Right> > > >,std::tuple<> >>;
      using composite_type=TupleChangeType<TestOrTrial::value,single_composite_type,emptytuple>;

  };

template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename Left,typename Right>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Subtraction< Expression<Left>,Expression<Right> > > > > >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< Subtraction<Expression<Left>,Expression<Right>> > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator< Expression<Subtraction< Expression<Left>,Expression<Right> > > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator< Expression<Subtraction< Expression<Left>,Expression<Right> > > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");

      using TestOrTrialLeft=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Left>>>;
      using TestOrTrialRight=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Right>>>;
      using LeftT=typename OperatorAndQuadratureTupleType<TestOrTrialLeft>::type;
      using RightT=typename OperatorAndQuadratureTupleType<TestOrTrialRight>::type;
      using type=RemoveTupleOfTupleDuplicates< LeftT, RightT >;
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using single_composite_type=std::tuple<std::tuple< CompositeOperator<Expression< Subtraction< Expression<Left>,Expression<Right> > > >,std::tuple<> >>;
      using composite_type=TupleChangeType<TestOrTrial::value,single_composite_type,emptytuple>;

  };

template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename Left,typename Right>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Multiplication< Expression<Left>,Expression<Right> > > > > >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< Multiplication<Expression<Left>,Expression<Right>> > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator<Expression<Multiplication< Expression<Left>,Expression<Right> > > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator<Expression<Multiplication< Expression<Left>,Expression<Right> > > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");

      using TestOrTrialLeft=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Left>>>;
      using TestOrTrialRight=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Right>>>;
      using LeftT=typename OperatorAndQuadratureTupleType<TestOrTrialLeft>::type;
      using RightT=typename OperatorAndQuadratureTupleType<TestOrTrialRight>::type;
      using type=RemoveTupleOfTupleDuplicates< LeftT, RightT >;
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using single_composite_type=std::tuple<std::tuple< CompositeOperator<Expression< Multiplication< Expression<Left>,Expression<Right> > > >,std::tuple<> >>;
      using composite_type=TupleChangeType<TestOrTrial::value,single_composite_type,emptytuple>;

  };


template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename Left,typename Right>
  class OperatorAndQuadratureTupleType<TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Division< Expression<Left>,Expression<Right> > > > > >
  { public:
      using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression< Division<Expression<Left>,Expression<Right>> > > >;
      static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,CompositeOperator<Expression<Division< Expression<Left>,Expression<Right> > > > > >::value ||
        IsSame<TestOrTrial,Trial<MixedSpace,N,CompositeOperator<Expression<Division< Expression<Left>,Expression<Right> > > > > >::value )
      && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");

      using TestOrTrialLeft=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Left>>>;
      using TestOrTrialRight=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Right>>>;
      using LeftT=typename OperatorAndQuadratureTupleType<TestOrTrialLeft>::type;
      using RightT=typename OperatorAndQuadratureTupleType<TestOrTrialRight>::type;
      using type=RemoveTupleOfTupleDuplicates< LeftT, RightT >;
      static constexpr Integer Nmax= TestOrTrial::Nmax;
      using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
      using single_composite_type=std::tuple<std::tuple< CompositeOperator<Expression< Division< Expression<Left>,Expression<Right> > > >,std::tuple<> >>;
      using composite_type=TupleChangeType<TestOrTrial::value,single_composite_type,emptytuple>;
  };




template<typename Left,typename Right,Integer QR>
class OperatorAndQuadratureTupleType<L2DotProductIntegral<Left,Right,QR>>
{ 
public:
  using L2prod=L2DotProductIntegral<Left,Right,QR>;
  using QRule=typename L2prod::QRule;
  using Operatortuple=typename OperatorAndQuadratureTupleType<typename L2prod::type>::type;
  using type=TupleOfTupleChangeType<1,QRule, Operatortuple>;
  using CompositeOperatortuple=typename OperatorAndQuadratureTupleType<typename L2prod::type>::composite_type;
  using composite_type=TupleOfTupleChangeType<1,QRule, CompositeOperatortuple>;
} 
; 



// order(T) = order(+T)
template<typename T>
    class OperatorAndQuadratureTupleType< Transposed< Expression<T> > >
    { public:
      using type=typename OperatorAndQuadratureTupleType<T>::type;
      using composite_type=typename OperatorAndQuadratureTupleType<T>::composite_type;
  };






}
#endif