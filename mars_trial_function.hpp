#ifndef MARS_TRIAL_FUNCTION_HPP
#define MARS_TRIAL_FUNCTION_HPP
#include "mars_base.hpp"
#include "mars_operators.hpp"

namespace mars {
template<typename MixedSpace, Integer N, typename OperatorType>
class Trial;

template<typename MixedSpace,Integer N, typename OperatorType=IdentityOperator>
class Trial: public Expression<Trial<MixedSpace,N,OperatorType>>
{ public:
  using FunctionSpace=MixedSpace;
  using Elem=typename FunctionSpace::Elem;
  static constexpr Integer Dim=Elem::Dim;
  static constexpr Integer ManifoldDim=Elem::ManifoldDim;
  using MeshT=Mesh<Dim,ManifoldDim>;
  using UniqueElementFunctionSpacesTupleType=typename MixedSpace::UniqueElementFunctionSpacesTupleType;
  static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
  static constexpr Integer number=N;
  static constexpr Integer value=GetType<typename MixedSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
  using Operator=OperatorType;

  Trial(const std::shared_ptr<MixedSpace>& W_ptr):
  spaces_ptr_(W_ptr)
  {}

  Trial(const MixedSpace& W):
  spaces_ptr_(std::make_shared<MixedSpace>(W))
  {}

  inline auto spaces_ptr()const {return spaces_ptr_;};

  private:
  std::shared_ptr<MixedSpace> spaces_ptr_;
};

template<typename MixedSpace, Integer N, typename Expr>
class Trial<MixedSpace,N,CompositeOperator<Expression<Expr>>>: 
public Expression<Trial< MixedSpace,N,CompositeOperator<Expression<Expr>> > >
{
public:
  using FunctionSpace=MixedSpace;
  using Elem=typename FunctionSpace::Elem;
  static constexpr Integer Dim=Elem::Dim;
  static constexpr Integer ManifoldDim=Elem::ManifoldDim;
  using MeshT=Mesh<Dim,ManifoldDim>;
  using UniqueElementFunctionSpacesTupleType=typename MixedSpace::UniqueElementFunctionSpacesTupleType;
  static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
  static constexpr Integer number=N;
  static constexpr Integer value=GetType<typename MixedSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
  using Operator=CompositeOperator<Expression<Expr>>;

  Trial(const std::shared_ptr<MixedSpace>& W_ptr, const Operator& op):
  spaces_ptr_(W_ptr),
  op_(op)
  {}

  Trial(const MixedSpace& W, const Operator& op):
  spaces_ptr_(std::make_shared<MixedSpace>(W)),
  op_(op)
  {}

  inline auto spaces_ptr()const {return spaces_ptr_;};
  inline constexpr auto composite_operator()const {return op_;};

  private:
  std::shared_ptr<MixedSpace> spaces_ptr_;
  Operator op_;
};





template<Integer N,typename...Args >
auto MakeTrial(const FunctionSpace<Args...>& W)
{return Trial<FunctionSpace<Args...>,N>(W);}

template<Integer N,typename...Args >
auto MakeTrial(const std::shared_ptr<FunctionSpace<Args...>>& W_ptr)
{return Trial<FunctionSpace<Args...>,N>(W_ptr);}


template<Integer N,typename...Args >
auto MakeTrial(const MixedSpace<Args...>& W)
{return Trial<MixedSpace<Args...>,N>(W);}

template<Integer N,typename...Args >
auto MakeTrial(const std::shared_ptr<MixedSpace<Args...>>& W_ptr)
{return Trial<MixedSpace<Args...>,N>(W_ptr);}


template<Integer N,typename...Args >
auto MakeTrial(const FullSpace<Args...>& W)
{return Trial<FullSpace<Args...>,N>(W);}

template<Integer N,typename...Args >
auto MakeTrial(const std::shared_ptr<FullSpace<Args...>>& W_ptr)
{return Trial<FullSpace<Args...>,N>(W_ptr);}




template<typename MixedSpace,Integer N>
Trial<MixedSpace,N,GradientOperator> 
Grad(const Trial<MixedSpace,N,IdentityOperator>& t)
{return Trial<MixedSpace,N,GradientOperator> (t.spaces_ptr());}

template<typename MixedSpace,Integer N>
Trial<MixedSpace,N,DivergenceOperator> 
Div(const Trial<MixedSpace,N,IdentityOperator>& t)
{return Trial<MixedSpace,N,DivergenceOperator> (t.spaces_ptr());}

template<typename MixedSpace,Integer N>
Trial<MixedSpace,N,CurlOperator> 
Curl(const Trial<MixedSpace,N,IdentityOperator>& t)
{return Trial<MixedSpace,N,CurlOperator> (t.spaces_ptr());}

template<typename MixedSpace,Integer N>
Trial<MixedSpace,N,TraceOperator> 
Trace(const Trial<MixedSpace,N,IdentityOperator>& t)
{return Trial<MixedSpace,N,TraceOperator> (t.spaces_ptr());}

template<typename MixedSpace,Integer N>
Trial<MixedSpace,N,TraceGradientOperator> 
TraceGrad(const Trial<MixedSpace,N,IdentityOperator>& t)
{return Trial<MixedSpace,N,TraceGradientOperator> (t.spaces_ptr());}





} 
#endif