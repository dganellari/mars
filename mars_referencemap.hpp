
#ifndef MARS_REFERENCE_MAP_HPP
#define MARS_REFERENCE_MAP_HPP

#include "mars_simplex.hpp"
#include "mars_base_elementfunctionspace.hpp"
#include "mars_vector.hpp"
#include "mars_matrix.hpp"
#include "mars_operators.hpp"
#include "mars_jacobian.hpp" 
namespace mars{



template<typename Operator,typename Elem,Integer FEFamily> 
class MapFromReference5;



template<Integer Dim, Integer ManifoldDim> 
class MapFromReference5<IdentityOperator, Simplex<Dim,ManifoldDim>,LagrangeFE> 
{
 public:
 using Operator=IdentityOperator;
 using Elem=Simplex<Dim,ManifoldDim>;
 // using type=Real;
        inline constexpr void  init(const Jacobian<Simplex<Dim,ManifoldDim>>& J){id_=1.0;}
        inline constexpr const auto&  operator()()const{return id_;}
 private:
  Real id_;
};






template<Integer Dim, Integer ManifoldDim> 
class MapFromReference5<GradientOperator, Simplex<Dim,ManifoldDim>,LagrangeFE> 
{
 public:
 using Operator=GradientOperator;
 using Elem=Simplex<Dim,ManifoldDim>;
 // using type=Matrix<Real, Dim, ManifoldDim>;
        inline constexpr void  init(const Jacobian<Simplex<Dim,ManifoldDim>>& J){grad_=  inverse(J());}
        inline constexpr const auto&  operator() ()const{return grad_;}
 private:
  Matrix<Real, Dim, ManifoldDim> grad_; 
};



template<Integer Dim, Integer ManifoldDim> 
class MapFromReference5<IdentityOperator,Simplex<Dim,ManifoldDim>,RaviartThomasFE> 
{
 public:
 using Operator=IdentityOperator;
 using Elem=Simplex<Dim,ManifoldDim>;
 // using type=Matrix<Real, Dim, ManifoldDim>;
         inline constexpr void init(const Jacobian<Simplex<Dim,ManifoldDim>> &J){id_=J();id_/=J.get_det();}
         inline constexpr const auto&  operator()()const {return id_;}
 private:
 Matrix<Real, Dim, ManifoldDim> id_;
};

template<Integer Dim, Integer ManifoldDim>
class MapFromReference5<DivergenceOperator,Simplex<Dim,ManifoldDim>,RaviartThomasFE> 
{
 public:
 using Operator=DivergenceOperator;
 using Elem=Simplex<Dim,ManifoldDim>;
 // using type=Real;
         inline constexpr void init(const Jacobian<Simplex<Dim,ManifoldDim>> &J){div_= J.get_det();}
         inline constexpr const auto&  operator()()const{return div_;}
 private:
 Real div_;
};






template<template<class>class Unary,typename T, typename Elem_, Integer FEFamily>
class MapFromReference5<Unary<Expression<T>>,Elem_,FEFamily> 
{
 public:
 using Operator=T;
 using Elem=Elem_;
 using Map=MapFromReference5<T,Elem_,FEFamily>;
         inline constexpr void init(const Jacobian<Elem> &J){map_.init(J);}
         inline constexpr const auto&  operator()()const{return map_;}
 private:
 Map map_;
};

// template<typename T, typename Elem_, Integer FEFamily>
// class MapFromReference5<TraceOperator<Expression<T>>,Elem_,FEFamily> 
// {
//  public:
//  using Operator=T;
//  using Elem=Elem_;
//  using Map=MapFromReference5<T,Elem_,FEFamily>;
//          inline constexpr void init(const Jacobian<Elem> &J){map_.init(J);}
//          inline constexpr const auto&  operator()()const{return map_;}
//  private:
//  Map map_;
// };




















template<typename GeneralForm,typename...GeneralForms>
class ReferenceMaps2
{
public:
 
  using FunctionSpace=typename GeneralForm::FunctionSpace;
  using Form=MultipleAddition<typename GeneralForm::Form,typename GeneralForms::Form...>;
  using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
  using TupleOperatorsAndQuadrature= typename OperatorAndQuadratureTupleType<Form>::type;
  using TupleCompositeOperatorsAndQuadrature= typename OperatorAndQuadratureTupleType<Form>::composite_type;
  using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadrature>;
  using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;

  using Map=MapOperatorTupleOfTuple<TupleOfTupleNoQuadrature,UniqueElementFunctionSpacesTupleType>;
  using UniqueMapping=UniqueMap<SpacesToUniqueFEFamilies,Map> ;


  template<Integer Nmax,Integer N,typename Elem, typename Map>
  typename std::enable_if_t< (N>Nmax),void> 
  unique_mapping_init_aux_aux(Map& map, const Jacobian<Elem> &J)
  {}


  template<Integer Nmax,Integer N,typename Elem, typename Map>
  typename std::enable_if_t< (N<=Nmax),void> 
  unique_mapping_init_aux_aux(Map& map, const Jacobian<Elem> &J)

  {
    unique_mapping_init_aux_aux<Nmax,N+1>(map,J);
    std::get<N>(map).init(J);
  }

  template<Integer Nmax,Integer N,typename Elem>
  typename std::enable_if_t< (N>Nmax),void> 
  unique_mapping_init_aux(const Jacobian<Elem> &J)
  {}
 
  template<Integer Nmax,Integer N,typename Elem>
  typename std::enable_if_t< (N<=Nmax),void> 
  unique_mapping_init_aux(const Jacobian<Elem> &J)
  {
    unique_mapping_init_aux_aux<TupleTypeSize<GetType<UniqueMapping,N>>::value-1,0>(std::get<N>(tuple_maps_),J);
    unique_mapping_init_aux<Nmax,N+1>(J);
  }


  template<typename Elem>
  constexpr void init(const Jacobian<Elem> &J)
  {
   unique_mapping_init_aux<TupleTypeSize<UniqueMapping>::value-1,0>(J);
  }

   const auto& operator()()const{return tuple_maps_;}
         auto& operator()()     {return tuple_maps_;}
         
   template<Integer...Ns>
   const auto& get()const{return tuple_get<Ns...>(tuple_maps_)();}

   template<Integer...Ns>
         auto& get()     {return tuple_get<Ns...>(tuple_maps_)();}
   
private:
  UniqueMapping tuple_maps_;
};


template<typename ConstFormReference,typename...ConstFormReferences>
constexpr auto reference_maps2(const ConstFormReference& form,const ConstFormReferences&...forms)
{
 return ReferenceMaps2<ConstFormReference,ConstFormReferences...>();

}


 // using Form=typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type;
 // typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type,
 //                       typename std::remove_const<typename std::remove_reference<ConstFormReferences>::type>::type...>();  


}
#endif
