
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
class MapFromReference;



template<Integer Dim, Integer ManifoldDim> 
class MapFromReference<IdentityOperator, Simplex<Dim,ManifoldDim>,LagrangeFE> 
{
 public:
 using Operator=IdentityOperator;
 using Elem=Simplex<Dim,ManifoldDim>;
 // using type=Real;
        // if we have functions defined on manifolddim-d, d>0
        template<Integer ManifoldDim2>
        inline constexpr void  init(const FiniteElem<Simplex<Dim,ManifoldDim2>>& FE){id_=1.0;}
        inline constexpr const auto&  operator()()const{return id_;}
 private:
  Real id_;
};

template<Integer Dim, Integer ManifoldDim> 
class MapFromReference<TraceOperator, Simplex<Dim,ManifoldDim>,LagrangeFE> 
{
 public:
 using Operator=TraceOperator;
 using Elem=Simplex<Dim,ManifoldDim+1>;
 // todo fixme. check if it is correct to put +1
 // indeed the element is Simplex<Dim,ManifoldDim+1> and we pass FiniteElem<Simplex<Dim,ManifoldDim+1>>

        inline constexpr void  init(const FiniteElem<Simplex<Dim,ManifoldDim+1>>& FE){id_=1.0;}
        inline constexpr const auto&  operator()()const{return id_;}
 private:
  Real id_;
};




template<Integer Dim, Integer ManifoldDim> 
class MapFromReference<GradientOperator, Simplex<Dim,ManifoldDim>,LagrangeFE> 
{
 public:
 using Operator=GradientOperator;
 using Elem=Simplex<Dim,ManifoldDim>;
 // using type=Matrix<Real, Dim, ManifoldDim>;
        inline constexpr void  init(const FiniteElem<Simplex<Dim,ManifoldDim>>& FE){grad_=  inverse(FE());}
        inline constexpr const auto&  operator() ()const{return grad_;}
 private:
  Matrix<Real, Dim, ManifoldDim> grad_; 
};



template<Integer Dim, Integer ManifoldDim> 
class MapFromReference<IdentityOperator,Simplex<Dim,ManifoldDim>,RaviartThomasFE> 
{
 public:
 using Operator=IdentityOperator;
 using Elem=Simplex<Dim,ManifoldDim>;
 // using type=Matrix<Real, Dim, ManifoldDim>;
         inline constexpr void init(const FiniteElem<Simplex<Dim,ManifoldDim>> &FE){id_=FE();id_/=FE.get_det();}
         inline constexpr const auto&  operator()()const {return id_;}
 private:
 Matrix<Real, Dim, ManifoldDim> id_;
};

template<Integer Dim, Integer ManifoldDim>
class MapFromReference<DivergenceOperator,Simplex<Dim,ManifoldDim>,RaviartThomasFE> 
{
 public:
 using Operator=DivergenceOperator;
 using Elem=Simplex<Dim,ManifoldDim>;
         inline constexpr void init(const FiniteElem<Simplex<Dim,ManifoldDim>> &FE){div_= FE.get_det();}
         inline constexpr const auto&  operator()()const{return div_;}
 private:
 Real div_;
};



/////////////////////////////////////////////////////////////////////////////////////
///// For RaviartThomasFE of every order, we have that:                     /////////
///// Trace(RT_k)= RT_k*n|_{boundary}= (RT0 * n) * P_k |_{boundary}         /////////
/////            = 1/|boundary| * P_k |_{boundary}                          /////////
///// where:                                                                /////////
///// p_k is the polynomial of order k on the boundary and                  /////////
//// |boundary|is the measure of the boundary                               /////////
///// Since LagrangeFE have a map=1, here the map is given by 1/|boundary|  /////////
///// We assume however that the normal is outward                          /////////
/////////////////////////////////////////////////////////////////////////////////////

template<Integer Dim, Integer ManifoldDim>
class MapFromReference<TraceOperator,Simplex<Dim,ManifoldDim>,RaviartThomasFE> 
{
 public:
 using Operator=TraceOperator;
 using Elem=Simplex<Dim,ManifoldDim+1>;

 // todo fixme. check if it is correct to put +1
 // indeed the element is Simplex<Dim,ManifoldDim+1> and we pass FiniteElem<Simplex<Dim,ManifoldDim+1>>
         inline constexpr void init(const FiniteElem<Simplex<Dim,ManifoldDim+1>> &FE)
         {flux_inv_= 1.0/((ManifoldDim+1)*FE.volume());}
         inline constexpr const auto&  operator()()const{return flux_inv_;}
 private:
 Real flux_inv_;
};



template<template<class>class Unary,typename T, typename Elem_, Integer FEFamily>
class MapFromReference<Unary<Expression<T>>,Elem_,FEFamily> 
{
 public:
 using Operator=T;
 using Elem=Elem_;
 using Map=MapFromReference<T,Elem_,FEFamily>;
         inline constexpr void init(const FiniteElem<Elem> &FE){map_.init(FE);}
         inline constexpr const auto&  operator()()const{return map_;}
 private:
 Map map_;
};

// template<typename T, typename Elem_, Integer FEFamily>
// class MapFromReference<TraceOperator<Expression<T>>,Elem_,FEFamily> 
// {
//  public:
//  using Operator=T;
//  using Elem=Elem_;
//  using Map=MapFromReference<T,Elem_,FEFamily>;
//          inline constexpr void init(const FiniteElem<Elem> &J){map_.init(J);}
//          inline constexpr const auto&  operator()()const{return map_;}
//  private:
//  Map map_;
// };




















template<typename GeneralForm,typename...GeneralForms>
class MapFromReferenceCollection
{
public:
 
  using FunctionSpace=typename GeneralForm::FunctionSpace;
  using Form=MultipleAddition<typename GeneralForm::Form,typename GeneralForms::Form...>;
  using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
  using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;

  using TupleOperatorsAndQuadratureVolumetric= typename OperatorAndQuadratureTupleType<Form,0>::type;
  using TupleCompositeOperatorsAndQuadratureVolumetric= typename OperatorAndQuadratureTupleType<Form,0>::composite_type;

  using TupleOperatorsAndQuadratureSurface= typename OperatorAndQuadratureTupleType<Form,1>::type;
  using TupleCompositeOperatorsAndQuadratureSurface= typename OperatorAndQuadratureTupleType<Form,1>::composite_type;

  // H=0, we consider only volume shape functions

  using TupleOfTupleNoQuadratureVolumetric=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadratureVolumetric>;
  using TupleOfTupleNoQuadratureSurface=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadratureSurface>;

  using MapVolumetric=MapOperatorTupleOfTuple<TupleOfTupleNoQuadratureVolumetric,UniqueElementFunctionSpacesTupleType>;
  using MapSurface=MapOperatorTupleOfTuple<TupleOfTupleNoQuadratureSurface,UniqueElementFunctionSpacesTupleType>;

  using UniqueMappingVolumetric=UniqueMap<SpacesToUniqueFEFamilies,MapVolumetric> ;
  using UniqueMappingSurface=UniqueMap<SpacesToUniqueFEFamilies,MapSurface> ;


  // template<Integer Nmax,Integer N,typename Elem, typename Map>
  // typename std::enable_if_t< (IsSame<Map,std::tuple<>>::value),void> 
  // unique_mapping_init_aux_aux(Map& map, const FiniteElem<Elem> &FE)
  // {}

  template<Integer Nmax,Integer N,typename Elem, typename Map>
  typename std::enable_if_t< ( (N>Nmax) || (IsSame<Map,std::tuple<>>::value)),void> 
  unique_mapping_init_aux_aux(Map& map, const FiniteElem<Elem> &FE)
  {}


  template<Integer Nmax,Integer N,typename Elem, typename Map>
  typename std::enable_if_t< ((N<=Nmax) && (IsDifferent<Map,std::tuple<>>::value)),void> 
  unique_mapping_init_aux_aux(Map& map, const FiniteElem<Elem> &FE)

  {
    unique_mapping_init_aux_aux<Nmax,N+1>(map,FE);
    std::get<N>(map).init(FE);
  }

  template<typename KindOfMapping, Integer Nmax,Integer N,typename Tuple,typename Elem>
  typename std::enable_if_t< (N>Nmax),void> 
  unique_mapping_init_aux(Tuple&tuple, const FiniteElem<Elem> &FE)
  {}
 
  template<typename KindOfMapping, Integer Nmax,Integer N,typename Tuple,typename Elem>
  typename std::enable_if_t< (N<=Nmax),void> 
  unique_mapping_init_aux(Tuple&tuple, const FiniteElem<Elem> &FE)
  {
    unique_mapping_init_aux_aux<TupleTypeSize<GetType<KindOfMapping,N>>::value-1,0>(std::get<N>(tuple),FE);
    unique_mapping_init_aux<KindOfMapping,Nmax,N+1>(tuple,FE);
  }


  template<typename Elem>
  constexpr void init(const FiniteElem<Elem> &FE)
  {
   unique_mapping_init_aux<UniqueMappingVolumetric,TupleTypeSize<UniqueMappingVolumetric>::value-1,0>(tuple_maps_volumetric_,FE);
  }

  template<typename Elem>
  constexpr void init_boundary(const FiniteElem<Elem> &FE)
  {
   unique_mapping_init_aux<UniqueMappingSurface,TupleTypeSize<UniqueMappingSurface>::value-1,0>(tuple_maps_surface_,FE);
  }

   // const auto& operator()()const{return tuple_maps_volumetric_;}
   //       auto& operator()()     {return tuple_maps_volumetric_;}

   const auto& volumetric_map()const{return tuple_maps_volumetric_;}
         auto& volumetric_map()     {return tuple_maps_volumetric_;}

   const auto& surface_map()const{return tuple_maps_surface_;}
         auto& surface_map()     {return tuple_maps_surface_;}

      
   template<Integer...Ns>
   const auto& get()const{return tuple_get<Ns...>(tuple_maps_volumetric_)();}

   template<Integer...Ns>
         auto& get()     {return tuple_get<Ns...>(tuple_maps_volumetric_)();}
   
private:
  UniqueMappingVolumetric tuple_maps_volumetric_;
  UniqueMappingSurface tuple_maps_surface_;
};


template<typename ConstFormReference,typename...ConstFormReferences>
constexpr auto reference_maps2(const ConstFormReference& form,const ConstFormReferences&...forms)
{
 return MapFromReferenceCollection<ConstFormReference,ConstFormReferences...>();

}


 // using Form=typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type;
 // typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type,
 //                       typename std::remove_const<typename std::remove_reference<ConstFormReferences>::type>::type...>();  


}
#endif
