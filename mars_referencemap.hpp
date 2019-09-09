
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
 // using Jacobian=Matrix<Real, Dim, ManifoldDim>;
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
        // inline void  init(const Jacobian& J){grad_=  inverse(J);}
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
 // using Jacobian=Matrix<Real, Dim, ManifoldDim>;
     
         // inline void init(const Jacobian& J) {id_=J;id_/=det(id_);}
         inline constexpr void init(const Jacobian<Simplex<Dim,ManifoldDim>> &J){id_=J();id_/=J.det();}

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
 // using Jacobian=Matrix<Real, Dim, ManifoldDim>;   
         inline constexpr void init(const Jacobian<Simplex<Dim,ManifoldDim>> &J){div_= J.det();}
         // inline void init(const Jacobian& J){div_= 1.0/det(J);}
         inline constexpr const auto&  operator()()const{return div_;}
 private:
 Real div_;
};


  // template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim, typename...Maps>
  // typename std::enable_if_t<(N>Nmax),void> 
  // unique_mapping_init_aux_aux(std::tuple<Maps...>& tuple, const Matrix<Real,Dim,ManifoldDim> &J)
  // {}

  // template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim, typename...Maps>
  // typename std::enable_if_t< (N<=Nmax),void> 
  // unique_mapping_init_aux_aux(std::tuple<Maps...>& tuple, const Matrix<Real,Dim,ManifoldDim> &J)
  // {
  //   unique_mapping_init_aux_aux<Nmax,N+1>(tuple,J);
  //   std::get<N>(tuple).init(J);
  // }

  // template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim,typename...TuplesOfMaps>
  // typename std::enable_if_t< (N>Nmax),void> 
  // unique_mapping_init_aux(std::tuple<TuplesOfMaps...>& tuple, const Matrix<Real,Dim,ManifoldDim> &J)
  // {}

  // template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim,typename...TuplesOfMaps>
  // typename std::enable_if_t<N<=Nmax,void> 
  // unique_mapping_init_aux(std::tuple<TuplesOfMaps...>& tuple, const Matrix<Real,Dim,ManifoldDim> &J)
  // {
  //   constexpr Integer Nmax_aux=TupleTypeSize<GetType<std::tuple<TuplesOfMaps...>,N>>::value-1;
  //   unique_mapping_init_aux_aux<Nmax_aux,0>(std::get<N>(tuple),J);
  //   unique_mapping_init_aux<Nmax,N+1>(tuple,J);
  // }



  // template<Integer Dim, Integer ManifoldDim, typename...TuplesOfMaps>
  // constexpr void init(std::tuple<TuplesOfMaps...>& tuple, const Matrix<Real,Dim,ManifoldDim> &J)
  // {
  //  constexpr Integer Nmax=TupleTypeSize<std::tuple<TuplesOfMaps...>>::value-1;
  //  unique_mapping_init_aux<Nmax,0>(tuple,J);
  // }

// template<typename ConstFormReference>
// class ReferenceMaps
// {
// public:
//   using Form=typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type;
//   using UniqueElementFunctionSpacesTupleType=typename Form::UniqueElementFunctionSpacesTupleType;
//   using TupleOperatorsAndQuadrature= typename OperatorTupleType<Form>::type;
//   using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadrature>;
//   using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies<UniqueElementFunctionSpacesTupleType>;
//   using Map=MapOperatorTupleOfTuple<TupleOfTupleNoQuadrature,UniqueElementFunctionSpacesTupleType>;
//   using UniqueMapping=UniqueMap<SpacesToUniqueFEFamilies,Map> ;



//   template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim, typename Map>
//   typename std::enable_if_t<(N>Nmax),void> 
//   unique_mapping_init_aux_aux(Map& map, const Matrix<Real,Dim,ManifoldDim> &J)
//   {}

//   template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim, typename Map>
//   typename std::enable_if_t< (N<=Nmax),void> 
//   unique_mapping_init_aux_aux(Map& map, const Matrix<Real,Dim,ManifoldDim> &J)
//   {
//     unique_mapping_init_aux_aux<Nmax,N+1>(map,J);
//     std::get<N>(map).init(J);
//   }

//   template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim>
//   typename std::enable_if_t< (N>Nmax),void> 
//   unique_mapping_init_aux(const Matrix<Real,Dim,ManifoldDim> &J)
//   {}

//   template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim>
//   typename std::enable_if_t<(N<=Nmax),void> 
//   unique_mapping_init_aux(const Matrix<Real,Dim,ManifoldDim> &J)
//   {
//     unique_mapping_init_aux_aux<TupleTypeSize<GetType<UniqueMapping,N>>::value-1,0>(std::get<N>(tuple_maps_),J);
//     unique_mapping_init_aux<Nmax,N+1>(J);
//   }



//   template<Integer Dim, Integer ManifoldDim>
//   constexpr void init(const Matrix<Real,Dim,ManifoldDim> &J)
//   {
//    unique_mapping_init_aux<TupleTypeSize<UniqueMapping>::value-1,0>(J);
//   }

//    const auto& operator()()const{return tuple_maps_;}
//          auto& operator()()     {return tuple_maps_;}
//    template<Integer...Ns>
//    const auto& get()const{return tuple_get<Ns...>(tuple_maps_)();}

//    template<Integer...Ns>
//          auto& get()     {return tuple_get<Ns...>(tuple_maps_)();}

// private:
//   UniqueMapping tuple_maps_;
// };


// template<typename ConstFormReference>
// constexpr auto reference_maps(const ConstFormReference& form)
// {using Form=typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type;
//  return ReferenceMaps<Form>();  }




























template<typename ConstGeneralFormReference>
class ReferenceMaps2
{
public:
  using GeneralForm=typename std::remove_const<typename std::remove_reference<ConstGeneralFormReference>::type>::type;
  using Form=typename GeneralForm::Form;
  using FunctionSpace=typename GeneralForm::FunctionSpace;
  using UniqueElementFunctionSpacesTupleType=typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
  using TupleOperatorsAndQuadrature= typename OperatorTupleType<Form>::type;
  using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadrature>;
  using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies<UniqueElementFunctionSpacesTupleType>;
  using Map=MapOperatorTupleOfTuple<TupleOfTupleNoQuadrature,UniqueElementFunctionSpacesTupleType>;
  using UniqueMapping=UniqueMap<SpacesToUniqueFEFamilies,Map> ;



  // template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim, typename Map>
  // typename std::enable_if_t<(N>Nmax),void> 
  // unique_mapping_init_aux_aux(Map& map, const Matrix<Real,Dim,ManifoldDim> &J)
  template<Integer Nmax,Integer N,typename Elem, typename Map>
  typename std::enable_if_t< (N>Nmax),void> 
  unique_mapping_init_aux_aux(Map& map, const Jacobian<Elem> &J)
  {}

  // template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim, typename Map>
  // typename std::enable_if_t< (N<=Nmax),void> 
  // unique_mapping_init_aux_aux(Map& map, const Matrix<Real,Dim,ManifoldDim> &J)
  template<Integer Nmax,Integer N,typename Elem, typename Map>
  typename std::enable_if_t< (N<=Nmax),void> 
  unique_mapping_init_aux_aux(Map& map, const Jacobian<Elem> &J)

  {
    unique_mapping_init_aux_aux<Nmax,N+1>(map,J);
    std::get<N>(map).init(J);
  }

  // template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim>
  // typename std::enable_if_t< (N>Nmax),void> 
  // unique_mapping_init_aux(const Matrix<Real,Dim,ManifoldDim> &J)
  template<Integer Nmax,Integer N,typename Elem>
  typename std::enable_if_t< (N>Nmax),void> 
  unique_mapping_init_aux(const Jacobian<Elem> &J)
  {}

  // template<Integer Nmax,Integer N,Integer Dim, Integer ManifoldDim>
  // typename std::enable_if_t<(N<=Nmax),void> 
  template<Integer Nmax,Integer N,typename Elem>
  typename std::enable_if_t< (N<=Nmax),void> 
  unique_mapping_init_aux(const Jacobian<Elem> &J)
  // unique_mapping_init_aux(const Matrix<Real,Dim,ManifoldDim> &J)
  {
    unique_mapping_init_aux_aux<TupleTypeSize<GetType<UniqueMapping,N>>::value-1,0>(std::get<N>(tuple_maps_),J);
    unique_mapping_init_aux<Nmax,N+1>(J);
  }



  // template<Integer Dim, Integer ManifoldDim>
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


template<typename ConstFormReference>
constexpr auto reference_maps2(const ConstFormReference& form)
{using Form=typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type;
 return ReferenceMaps2<Form>();  }


}
#endif
