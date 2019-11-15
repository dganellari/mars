#ifndef MARS_OPERATOR_DIFFERENTIAL_HPP
#define MARS_OPERATOR_DIFFERENTIAL_HPP
#include "mars_base.hpp"

namespace mars {


class IdentityOperator: public Expression<IdentityOperator>  
{public: IdentityOperator(){};}; 
class DivergenceOperator: public Expression<DivergenceOperator>
{public: DivergenceOperator(){};};

class GradientOperator : public Expression<GradientOperator> 
{public: GradientOperator(){}};

class CurlOperator: public Expression<CurlOperator>  
{public: CurlOperator(){}};

class TraceOperator: public Expression<TraceOperator>  
{public: TraceOperator(){}};

template<Integer N_>
class ComponentOperator: public Expression<ComponentOperator<N_>>  
 {public: static constexpr Integer N=N_;
  ComponentOperator(){}};

 template<typename...Ts>
 class TraceElemHelper;

 template<template<Integer,Integer> class Elem_,Integer Dim,Integer ManifoldDim>
 class TraceElemHelper<Elem_<Dim,ManifoldDim>>
 {
  public:
    using type=Elem_<Dim,ManifoldDim-1>;
 };

 template<class Elem>
 using TraceOf=typename TraceElemHelper<Elem>::type;


}
#endif