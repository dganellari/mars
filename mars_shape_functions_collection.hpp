#ifndef MARS_SHAPE_FUNCTION_COLLECTION_HPP
#define MARS_SHAPE_FUNCTION_COLLECTION_HPP
#include "mars_general_form.hpp"


namespace mars {





template<typename GeneralForm_, typename...GeneralForms_>
class ShapeFunctionsCollection
{
 public:
  using GeneralForm=GeneralForm_;
  using Elem=typename GeneralForm::FunctionSpace::Elem;
  using Form=MultipleAddition<typename GeneralForm_::Form,typename GeneralForms_::Form...> ;  
  using UniqueElementFunctionSpacesTupleType=typename GeneralForm::UniqueElementFunctionSpacesTupleType;
  using TupleOperatorsAndQuadrature= typename OperatorAndQuadratureTupleType<Form>::type;
  using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadrature>;
  using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;

  using Map=MapOperatorTupleOfTuple<TupleOfTupleNoQuadrature,UniqueElementFunctionSpacesTupleType>;
  using UniqueMapping=UniqueMap<SpacesToUniqueFEFamilies,Map> ;
  using TupleOfTupleShapeFunction=TupleOfTupleShapeFunctionType2<UniqueElementFunctionSpacesTupleType,TupleOperatorsAndQuadrature>;
  
  using TupleCompositeOperatorsAndQuadrature= typename OperatorAndQuadratureTupleType<Form>::composite_type;
  // using TupleOfTupleShapeFunctionCombinations=TupleOfTupleShapeFunctionType3<UniqueElementFunctionSpacesTupleType,
  //                                                                            TupleCompositeOperatorsAndQuadrature>;

  using MapTupleNumbers=MapTupleInit<TupleOfTupleShapeFunction, SpacesToUniqueFEFamilies, UniqueMapping>;
  
  using TupleOfTupleCompositeShapeFunction=typename TupleOfCombinationFunctions<GeneralForm::FunctionSpace::Nuniquesubspaces,Form>::type;
  using TupleOfTupleCompositeShapeFunctionTensor=typename TupleOfCombinationFunctions<GeneralForm::FunctionSpace::Nuniquesubspaces,Form>::type_tensor;

   template<Integer M,Integer Nmax_aux,Integer N,typename Tuple,typename Map>
   constexpr typename std::enable_if< (N>Nmax_aux) ,void>::type 
   init_map_aux_aux(Tuple& t,const Map& maps){}

   template<Integer M,Integer Nmax_aux,Integer N,typename Tuple,typename Map>
   constexpr typename std::enable_if< (N<=Nmax_aux) ,void>::type 
   init_map_aux_aux(Tuple& t,const Map& maps) 
   {
    auto& t_nth=std::get<N>(t); 
    auto& map_nth=std::get<GetType<MapTupleNumbers,M,N>::value>(maps); 
    t_nth.init_map(map_nth);
    init_map_aux_aux<M,Nmax_aux,N+1>(t,maps);
    }



   template<Integer Nmax_aux,Integer N,typename Map>
   constexpr typename std::enable_if< (N>Nmax_aux) ,void>::type 
   init_map_aux(const Map& maps){}

   template<Integer Nmax_aux,Integer N,typename Map>
   constexpr typename std::enable_if< (N<=Nmax_aux) ,void>::type 
   init_map_aux(const Map& maps) 
   {
    auto& t_nth=std::get<N>(tuple); 
    auto& map_nth=std::get<GetType<SpacesToUniqueFEFamilies,N>::value>(maps); 
    init_map_aux_aux<N,TupleTypeSize<decltype(t_nth)>::value-1,0>(t_nth,map_nth);
    init_map_aux<Nmax_aux,N+1>(maps);
    }


   template<typename...Args>
   constexpr void init_map(const MapFromReferenceCollection<Args...>& maps) 
   {init_map_aux<TupleTypeSize<MapTupleNumbers>::value-1,0>(maps());}



  template<Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
  constexpr typename std::enable_if_t< FEFamily==LagrangeFE,void> 
  shape_function_init_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficient<Args...> &coefficients)
  {
    std::cout<<M<<" "<<FEFamily<<" "<<Order<<std::endl;
    shape.init();
  }


  template<Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
  constexpr typename std::enable_if_t<FEFamily==RaviartThomasFE,void> 
  shape_function_init_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficient<Args...> &coefficients)
  {
    shape.init(tuple_get<M>(coefficients()));
  }




  template<Integer M, Integer Nmax,Integer N,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t<(N>Nmax),void> 
  shape_function_init_aux_aux(Tuple& tuple, const ShapeFunctionCoefficient<Args...> &coefficients)
  {}

  template<Integer M,Integer Nmax,Integer N,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t< (N<=Nmax),void> 
  shape_function_init_aux_aux(Tuple& tuple, const ShapeFunctionCoefficient<Args...> &coefficients)
  {

    using Shape=GetType<TupleOfTupleShapeFunction,M,N>;
    using Elem=typename Shape::Elem;
    constexpr Integer FEFamily=Shape::FEFamily;
    constexpr Integer Order=Shape::Order;
    shape_function_init_aux_aux_aux<M,Elem,FEFamily,Order>(tuple_get<N>(tuple),coefficients);
    shape_function_init_aux_aux<M,Nmax,N+1>(tuple,coefficients);
  }

  template<Integer Nmax,Integer N,typename...Args>
  constexpr typename std::enable_if_t< (N>Nmax),void> 
  shape_function_init_aux(const ShapeFunctionCoefficient<Args...> &coefficients)
  {}

  template<Integer Nmax,Integer N,typename...Args>
  constexpr typename std::enable_if_t<N<=Nmax,void> 
  shape_function_init_aux(const ShapeFunctionCoefficient<Args...> &coefficients)
  {
    constexpr Integer Nmax_aux=TupleTypeSize<GetType<TupleOfTupleShapeFunction,N>>::value-1;
    shape_function_init_aux_aux<N,Nmax_aux,0>(tuple_get<N>(tuple),coefficients);
    shape_function_init_aux<Nmax,N+1>(coefficients);
  }



  //template<typename Coefficients>
  constexpr void init_shape_functions()//const Coefficients& shape_coefficients)
  {
   constexpr Integer Nmax=TupleTypeSize<TupleOfTupleShapeFunction>::value-1;
   shape_function_init_aux<Nmax,0>(coeffs_);
   }





















  template<Integer Nmax,Integer N,typename...Args1,typename...Args2>
  constexpr typename std::enable_if_t<(N>Nmax),void> 
  init_composite_shape_functions_aux_aux(std::tuple<Args1...>& tuple_tensor,const Jacobian<Elem> &J, std::tuple<Args2...>& tuple_composite)
  {
    std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;

  }

  template<Integer Nmax,Integer N,typename...Args1>
  constexpr void init_composite_shape_functions_aux_aux(std::tuple<>& tuple_tensor,const Jacobian<Elem> &J, std::tuple<Args1...>& tuple_composite)
  {
    std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;
  }

  template<Integer Nmax,Integer N,typename...Args1>
  constexpr void init_composite_shape_functions_aux_aux(std::tuple<Args1...>& tuple_tensor,const Jacobian<Elem> &J, std::tuple<>& tuple_composite)
  {
    std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;

  }

  template<Integer Nmax,Integer N>
  constexpr void init_composite_shape_functions_aux_aux(std::tuple<>& tuple_tensor,const Jacobian<Elem> &J, std::tuple<>& tuple_composite)
  {
    std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;

  }

  template<Integer Nmax,Integer N,typename...Args1,typename...Args2>
  constexpr typename std::enable_if_t< (N<=Nmax),void> 
  init_composite_shape_functions_aux_aux(std::tuple<Args1...>& tuple_tensor,const Jacobian<Elem> &J, std::tuple<Args2...>& tuple_composite)
  {
    
    auto& composite=tuple_get<N>(tuple_composite);
    auto& tensor=tuple_get<N>(tuple_tensor);
    // decltype(composite) a1(1);
    // decltype(tensor) a2(1);
    // decltype(tuple_tensor)d344(5);
        // TupleOfTupleCompositeShapeFunction ee4e(6);

    // TupleOfTupleCompositeShapeFunctionTensor eee(6);
    // Number<N> ee3(5);
    // Number<Nmax> ee(5);
    std::cout<<"aux aux Nmax,N="<<Nmax<<", "<<N<<std::endl;

    // todo fixme scommentami
    composite.apply(tensor,J,tuple);
    init_composite_shape_functions_aux_aux<Nmax,N+1>(tuple_tensor,J,tuple_composite);
  }


  template<Integer Nmax,Integer N>
  constexpr typename std::enable_if_t< (N>Nmax),void> 
  init_composite_shape_functions_aux(const Jacobian<Elem> &J)
  {}

  template<Integer Nmax,Integer N>
  constexpr typename std::enable_if_t<N<=Nmax,void> 
  init_composite_shape_functions_aux(const Jacobian<Elem> &J)
  {
    auto& tensor=tuple_get<N>(tuple_tensors_);
    const auto& composite=tuple_get<N>(tuple_composite_);
    constexpr Integer Nmax_aux=TupleTypeSize<decltype(composite)>::value-1;
    constexpr Integer Nmax_au2=TupleTypeSize<decltype(tensor)>::value-1;

    // Number<Nmax_aux> rr(6);

    // Number<Nmax_aux> e(6);
    // Number<Nmax_au2> r5e(6);
    // decltype(tensor) rf(56);
    // decltype(composite) r4f(56);
    // std::cout<<"aux Nmax,N="<<Nmax<<", "<<N<<std::endl;
    // std::cout<<"Nmax_aux="<<TupleTypeSize<GetType<TupleOfTupleCompositeShapeFunction,N>>::value<<" "<<Nmax_aux<<std::endl;
    // std::cout<<"Nmax_aux="<<TupleTypeSize<decltype(tensor)>::value-1<<" "<<Nmax_aux<<std::endl;
    // std::cout<<"Nmax_aux="<<TupleTypeSize<decltype(composite)>::value-1<<" "<<Nmax_aux<<std::endl;

    init_composite_shape_functions_aux_aux<Nmax_aux,0>(tuple_get<N>(tuple_tensors_),J,tuple_get<N>(tuple_composite_));
    init_composite_shape_functions_aux<Nmax,N+1>(J);
  }

  //template<typename Coefficients>
  constexpr void init_composite_shape_functions(const Jacobian<Elem> &J)
  {
    std::cout<<"init_composite_shape_functions"<<std::endl;
   constexpr Integer Nmax=TupleTypeSize<TupleOfTupleCompositeShapeFunction>::value-1;
   init_composite_shape_functions_aux<Nmax,0>(J);
   }



  constexpr void init(const Jacobian<Elem>&J)
  {
    // TODO CHECK: SINCE WE HAVE ALREADY HAVE REFERENCES TO MAPS AND COEFFS, do we have to init_map?
    // INIT MAP: takes the corresponding map and compute the shape function in the actual element
    std::cout<<"init maps"<<std::endl;
    init_map(maps_);
    // init_shape_functions: takes also the coefficients which multiply the actual shape functions
    std::cout<<"init_shape_functions"<<std::endl;
    init_shape_functions();
    std::cout<<"init_composite_shape_functions"<<std::endl;
    init_composite_shape_functions(J);
   }




   constexpr       auto& operator()()      {return tuple;}
   constexpr const auto& operator()()const {return tuple;}


   constexpr       auto & composite_shapes()      {return tuple_composite_;}
   constexpr const auto & composite_shapes()const {return tuple_composite_;}

   constexpr       auto & composite_tensor()      {return tuple_tensors_;}
   constexpr const auto & composite_tensor()const {return tuple_tensors_;}

   template<Integer...Ns>
   constexpr const auto& get()const{return tuple_get<Ns...>(tuple);}

   template<Integer...Ns>
         auto& get()     {return tuple_get<Ns...>(tuple);}

   template<Integer...Ns>
   constexpr const auto& value()const{return tuple_get<Ns...>(tuple).eval();}

   template<Integer...Ns>
         auto& value()     {return tuple_get<Ns...>(tuple).eval();}


 ShapeFunctionsCollection(ShapeFunctionCoefficient<GeneralForm_,GeneralForms_...>&coeffs,
                 MapFromReferenceCollection<GeneralForm_,GeneralForms_...>& maps,
                 const GeneralForm_& form,const GeneralForms_&...forms):
 coeffs_(coeffs),
 maps_(maps)
 ,
 tuple_composite_(build_tuple_of_combination_functions<GeneralForm::FunctionSpace::Nuniquesubspaces>(form(),forms()...))
 { }


private:
   ShapeFunctionCoefficient<GeneralForm_,GeneralForms_...> & coeffs_;
   MapFromReferenceCollection<GeneralForm_,GeneralForms_...> & maps_;
   TupleOfTupleShapeFunction tuple;
   TupleOfTupleCompositeShapeFunction tuple_composite_;
   TupleOfTupleCompositeShapeFunctionTensor tuple_tensors_;
};


template<typename Form,typename...Forms>
constexpr auto shape_functions(ShapeFunctionCoefficient<Form,Forms...>&coeffs, MapFromReferenceCollection<Form,Forms...>&maps,const Form& form,const Forms&...forms)
{
  //using Form=typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type;
 return ShapeFunctionsCollection<Form,Forms...>(coeffs,maps,form,forms...);  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////// For explanation on how and why this works, check:
//////// https://stackoverflow.com/questions/1005476/how-to-detect-whether-there-is-a-specific-member-variable-in-class
//////// We check wheter T has a variable named is_test_or_trial
//////// Such a variable is used only by Test or Trial and respectively set to 0 or 1
//////// If in a linear or bilinear form, other quantities which are neither tests nor trials appear, then we can say so
//////// because they do not have such is_test_or_trial
//////// Therefore we can control how many is_test_or_trial appear and their respective values
//////// In a 0-form, no is_test_or_trial appear
//////// In a 1-form, is_test_or_trial=0 (test)
//////// In a 2-form, is_test_or_trial=0 and 1 (test and trial)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// template <typename T, typename = Integer>
// struct ITestOrTrial : std::false_type { };

// template <typename T>
// struct ITestOrTrial <T,  decltype((void) T::is_test_or_trial, static_cast<decltype(T::is_test_or_trial)>(T::is_test_or_trial) )> : std::true_type { };



}
#endif