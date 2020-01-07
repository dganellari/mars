#ifndef MARS_TEST_FUNCTION_HPP
#define MARS_TEST_FUNCTION_HPP
#include "mars_base.hpp"
#include "mars_operators.hpp"

namespace mars {


template<typename MixedSpace, Integer N, typename OperatorType>
class Test;


template<typename MeshT, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
class FunctionSpace; 

template<typename...Args>
class MixedSpace; 

template<typename...Args>
class FullSpace;


template<typename MixedSpace, Integer N, typename OperatorType=IdentityOperator>
class Test: public Expression<Test<MixedSpace,N,OperatorType>>
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
  using Operator=OperatorType;

  Test(const std::shared_ptr<MixedSpace>& W_ptr):
  spaces_ptr_(W_ptr)
  {}

  Test(const MixedSpace& W):
  spaces_ptr_(std::make_shared<MixedSpace>(W))
  {}

  inline auto spaces_ptr()const {return spaces_ptr_;};

  private:
  std::shared_ptr<MixedSpace> spaces_ptr_;
};


template<typename MixedSpace, Integer N, typename Expr>
class Test<MixedSpace,N,CompositeOperator<Expression<Expr>>>: 
public Expression<Test< MixedSpace,N,CompositeOperator<Expression<Expr>> > >
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

  Test(const std::shared_ptr<MixedSpace>& W_ptr, const Operator& op):
  spaces_ptr_(W_ptr),
  op_(op)
  {}

  Test(const MixedSpace& W, const Operator& op):
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
auto MakeTest(const FunctionSpace<Args...>& W)
{return Test<FunctionSpace<Args...>,N>(W);}

template<Integer N,typename...Args >
auto MakeTest(const std::shared_ptr<FunctionSpace<Args...>>& W_ptr)
{return Test<FunctionSpace<Args...>,N>(W_ptr);}




template<Integer N,typename...Args >
auto MakeTest(const MixedSpace<Args...>& W)
{return Test<MixedSpace<Args...>,N>(W);}

template<Integer N,typename...Args >
auto MakeTest(const std::shared_ptr<MixedSpace<Args...>>& W_ptr)
{return Test<MixedSpace<Args...>,N>(W_ptr);}




template<Integer N,typename...Args >
auto MakeTest(const FullSpace<Args...>& W)
{return Test<FullSpace<Args...>,N>(W);}

template<Integer N,typename...Args >
auto MakeTest(const std::shared_ptr<FullSpace<Args...>>& W_ptr)
{return Test<FullSpace<Args...>,N>(W_ptr);}




template<typename MixedSpace,Integer N>
Test<MixedSpace,N,GradientOperator> 
Grad(const Test<MixedSpace,N,IdentityOperator>& t)
{return Test<MixedSpace,N,GradientOperator> (t.spaces_ptr());}

template<typename MixedSpace,Integer N>
Test<MixedSpace,N,DivergenceOperator> 
Div(const Test<MixedSpace,N,IdentityOperator>& t)
{return Test<MixedSpace,N,DivergenceOperator> (t.spaces_ptr());}

template<typename MixedSpace,Integer N>
Test<MixedSpace,N,CurlOperator> 
Curl(const Test<MixedSpace,N,IdentityOperator>& t)
{return Test<MixedSpace,N,CurlOperator> (t.spaces_ptr());}

template<typename MixedSpace,Integer N>
Test<MixedSpace,N,TraceOperator> 
Trace(const Test<MixedSpace,N,IdentityOperator>& t)
{return Test<MixedSpace,N,TraceOperator> (t.spaces_ptr());}




}
#endif