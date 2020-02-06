#ifndef MARS_SHAPE_FUNCTION_COLLECTION_HPP
#define MARS_SHAPE_FUNCTION_COLLECTION_HPP
#include "mars_general_form.hpp"
#include "mars_dirichlet_bc.hpp"


namespace mars {


template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Operator,typename Expr>
constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,Operator>& t, 
                                              const Expr& expr)
{
    using T=TestOrTrial<MixedSpace,N,Operator>;
    using FunctionSpace=typename T::FunctionSpace;
    using FromElementFunctionSpacesToFirstSpaceTupleType=typename FunctionSpace::FromElementFunctionSpacesToFirstSpaceTupleType;
    constexpr Integer FirstSpace=GetType<FromElementFunctionSpacesToFirstSpaceTupleType,T::value>::value;
    return TestOrTrial<MixedSpace,FirstSpace,Operator>(t.spaces_ptr());
}

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Operator>
constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
                                              const Operator& op)
{
    using T=TestOrTrial<MixedSpace,N,Operator>;
    using FunctionSpace=typename T::FunctionSpace;
    using FromElementFunctionSpacesToFirstSpaceTupleType=typename FunctionSpace::FromElementFunctionSpacesToFirstSpaceTupleType;
    constexpr Integer FirstSpace=GetType<FromElementFunctionSpacesToFirstSpaceTupleType,T::value>::value;
    return TestOrTrial<MixedSpace,FirstSpace,Operator>(t.spaces_ptr());
}

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Operator>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Transposed<Expression<Operator>>& expr)
// {
//   auto e=form_of_composite_operator_aux(t,expr.derived());
//   // decltype(expr.derived()) eee(6);
//   return Transpose(e);

//   // using T=TestOrTrial<MixedSpace,N,Operator>;
//   // using FunctionSpace=typename T::FunctionSpace;
//   // using FromElementFunctionSpacesToFirstSpaceTupleType=typename FunctionSpace::FromElementFunctionSpacesToFirstSpaceTupleType;
//   // constexpr Integer FirstSpace=GetType<FromElementFunctionSpacesToFirstSpaceTupleType,T::value>::value;  
//   // return Transpose(TestOrTrial<MixedSpace,FirstSpace,Operator>(t.spaces_ptr()));
// }

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename ConstType,typename...Inputs>
constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
                                              const ConstantTensor<ConstType,Inputs...>& constant)
{return constant;}

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename ConstType,typename...Inputs>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Transposed<Expression<ConstantTensor<ConstType,Inputs...>>>& transposed_constant)
// {
//   return transposed_constant;
// }

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename ConstType,typename...Inputs>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const TraceOperator<Expression<ConstantTensor<ConstType,Inputs...>>>& trace_constant)
// {
//   return trace_constant;
// }


template<template<class,Integer,class > class TestOrTrial, template<class>class Unary,
typename MixedSpace,Integer N,typename Expr,typename ConstType,typename...Inputs>
constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
                                              const Unary<Expression<ConstantTensor<ConstType,Inputs...>>>& unary_operator_applied_to_constant)
{
    return unary_operator_applied_to_constant;
}

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename FullSpace, Integer M,typename Operator,typename FuncType>
constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
                                              const Function<FullSpace,M,Operator,FuncType>& func)
{return func;}

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename FullSpace, Integer M,typename Operator,typename FuncType>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Transposed<Expression<Function<FullSpace,M,Operator,FuncType>>>& transposed_func)
// {return transposed_func;}


template<template<class,Integer,class > class TestOrTrial, template<class>class Unary,
typename MixedSpace,Integer N,typename Expr,typename FullSpace, Integer M,typename Operator,typename FuncType>
constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
                                              const Unary<Expression<Function<FullSpace,M,Operator,FuncType>>>& unary_operator_applied_to_func)
{return unary_operator_applied_to_func;}












template<template<class,Integer,class > class TestOrTrial, template<class>class Unary,
typename MixedSpace,Integer N,typename Expr,typename T>
constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
                                              const Unary<Expression<T>>& expr)
{
    auto e=form_of_composite_operator_aux(t,expr.derived());
    // decltype(expr.derived()) eee(6);
    return Unary<Expression<decltype(e)>>(e);
}

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename T>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const UnaryPlus<Expression<T>>& expr)
// {
//   auto e=form_of_composite_operator_aux(t,expr.derived());
//   // decltype(expr.derived()) eee(6);
//   return +e;
// }

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename T>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const UnaryMinus<Expression<T>>& expr)
// {
//   auto e=form_of_composite_operator_aux(t,expr.derived());
//   return -e;
// }



template<template<class,Integer,class > class TestOrTrial, template<class,class>class Binary,
typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
                                              const Binary<Expression<Left>,Expression<Right>>& expr)
{
    auto left=form_of_composite_operator_aux(t,expr.left());
    auto right=form_of_composite_operator_aux(t,expr.right());
    
    return Binary<Expression<decltype(left)>,Expression<decltype(right)>>(left,right);
}


// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Addition<Expression<Left>,Expression<Right>>& expr)
// {
//   auto left=form_of_composite_operator_aux(t,expr.left());
//   auto right=form_of_composite_operator_aux(t,expr.right());

//   return left+right;
// }

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Subtraction<Expression<Left>,Expression<Right>>& expr)
// {
//   auto left=form_of_composite_operator_aux(t,expr.left());
//   auto right=form_of_composite_operator_aux(t,expr.right());

//   return left-right;
// }

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Multiplication<Expression<Left>,Expression<Right>>& expr)
// {
//   auto left=form_of_composite_operator_aux(t,expr.left());
//   auto right=form_of_composite_operator_aux(t,expr.right());

//   return left*right;
// }

// template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr,typename Left,typename Right>
// constexpr auto form_of_composite_operator_aux(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t, 
//                                               const Division<Expression<Left>,Expression<Right>>& expr)
// {
//   auto left=form_of_composite_operator_aux(t,expr.left());
//   auto right=form_of_composite_operator_aux(t,expr.right());

//   return left/right;
// }


template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Expr>
constexpr auto form_of_composite_operator(const TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>& t)
{
    return form_of_composite_operator_aux(t,t.composite_operator().composite_operator());
}

template<template<class,Integer,class > class TestOrTrial, typename MixedSpace,Integer N,typename Operator>
constexpr auto form_of_composite_operator(const TestOrTrial<MixedSpace,N,Operator>& t)
{
    return t;
}





















template<typename QuadratureRule, typename Tuple,typename Other>
constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, const Other& other)
{
 return tuple;
}

// template<typename QuadratureRule, typename Tuple,typename ConstType,typename...Inputs>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, const ConstantTensor<ConstType,Inputs...>& other)
// {
//  return tuple;
// }


// template<typename QuadratureRule, typename Tuple,typename FullSpace, Integer M,typename Operator,typename FuncType>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, const Function<FullSpace,M,Operator,FuncType>& other)
// {
//  return tuple;
// }

// template<typename QuadratureRule, typename Tuple,typename MixedSpace, Integer N, typename Operator>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, const Trial<MixedSpace,N,Operator>& other)
// {
//  return tuple;
// }
// template<typename QuadratureRule, typename Tuple,typename MixedSpace, Integer N, typename Operator>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, const Test<MixedSpace,N,Operator>& other)
// {
//  return tuple;
// }

template<typename QuadratureRule, typename Tuple,
         template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N, typename Expr>
constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple, 
                           const TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Expr>>>& testortrial)
{
 // check if already exists test or trial of the same input with quadrature ruel
  using TestOrTrial=TestOrTrial_<MixedSpace,N,CompositeOperator<Expression<Expr>>>;
  auto composite_testortrial=
       form_of_composite_operator(Test<MixedSpace,N,CompositeOperator<Expression<Expr>>>(testortrial.spaces_ptr(),testortrial.composite_operator()));
  auto eval_composite=Evaluation<Expression<decltype(composite_testortrial)>,QuadratureRule>(composite_testortrial);
  auto tuple_nth=tuple_get<TestOrTrial::value>(tuple);
  
  // decltype(composite_testortrial) eee(5);
  return tuple_change_element<TestOrTrial::value>(tuple,tuple_cat_unique(tuple_nth,decltype(eval_composite)(eval_composite)));
}









template<typename QuadratureRule, template<class>class Unary, typename Tuple,typename T>
constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
                           const Unary<Expression<T>>& unary)
{
  return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,unary.derived());
}

// template<typename QuadratureRule, typename Tuple,typename T>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
//                            const UnaryPlus<Expression<T>>& plusexpr)
// {
//   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,plusexpr.derived());
// }

// template<typename QuadratureRule, typename Tuple,typename T>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
//                            const UnaryMinus<Expression<T>>& minusexpr)
// {
//   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,minusexpr.derived());
// }


template<typename QuadratureRule, template<class,class>class Binary, typename Tuple,typename Left,typename Right>
constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
                           const Binary<Expression<Left>,Expression<Right>>& binary)
{
  auto tuple_new=build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,binary.left());
  return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple_new,binary.right());
}


// template<typename QuadratureRule, typename Tuple,typename Left,typename Right>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
//                            const Addition<Expression<Left>,Expression<Right>>& addition)
// {
//   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,addition.left());
//   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple_new,addition.right());
// }

// template<typename QuadratureRule, typename Tuple,typename Left,typename Right>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
//                            const Subtraction<Expression<Left>,Expression<Right>>& subtraction)
// {
//   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,subtraction.left());
//   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple_new,subtraction.right());
// }

// template<typename QuadratureRule, typename Tuple,typename Left,typename Right>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
//                            const Multiplication<Expression<Left>,Expression<Right>>& multiplication)
// {
//   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,multiplication.left());
//   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple_new,multiplication.right());
// }

// template<typename QuadratureRule, typename Tuple,typename Left,typename Right>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
//                            const Division<Expression<Left>,Expression<Right>>& division)
// {
//   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,division.left());
//   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple_new,division.right());
// }


// template<typename QuadratureRule, typename Tuple,typename T>
// constexpr auto build_tuple_of_combination_functions_aux_aux(const Tuple& tuple,
//                            const Transposed<Expression<T>>& transposed_expr)
// {
//   return build_tuple_of_combination_functions_aux_aux<QuadratureRule>(tuple,transposed_expr.derived());
// }


template<typename Tuple, typename Left,typename Right, Integer QR>
constexpr auto build_tuple_of_combination_functions_surface_aux(const Tuple& tuple,
                                                               const L2DotProductIntegral<Left,Right,true,QR>& l2prod)
{return tuple;}

template<typename Tuple, typename Left,typename Right, Integer QR>
constexpr auto build_tuple_of_combination_functions_surface_aux(const Tuple& tuple,
                                                        const L2DotProductIntegral<Left,Right,false,QR>& l2prod)
{
  using L2=L2DotProductIntegral<Left,Right,false,QR>;
  using QRule=typename L2::QRule;
  auto tuple_new=build_tuple_of_combination_functions_aux_aux<QRule>(tuple,l2prod.left());
  return build_tuple_of_combination_functions_aux_aux<QRule>(tuple_new,l2prod.right());
}


template<typename Tuple, typename Left,typename Right, Integer QR>
constexpr auto build_tuple_of_combination_functions_volume_aux(const Tuple& tuple,
                                                        const L2DotProductIntegral<Left,Right,false,QR>& l2prod)
{return tuple;}


template<typename Tuple, typename Left,typename Right, Integer QR>
constexpr auto build_tuple_of_combination_functions_volume_aux(const Tuple& tuple,
                                                        const L2DotProductIntegral<Left,Right,true,QR>& l2prod)
{
  using L2=L2DotProductIntegral<Left,Right,true,QR>;
  using QRule=typename L2::QRule;
  auto tuple_new=build_tuple_of_combination_functions_aux_aux<QRule>(tuple,l2prod.left());
  return build_tuple_of_combination_functions_aux_aux<QRule>(tuple_new,l2prod.right());
}







// template<typename Tuple, typename Left,typename Right,bool VolumeIntegral, Integer QR>
// constexpr auto build_tuple_of_combination_functions_aux(const Tuple& tuple,
//                                                         const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod)
// {
//   return build_tuple_of_combination_functions_volume_aux(tuple,l2prod);
// }



// template<typename Tuple, typename Left,typename Right,bool VolumeIntegral, Integer QR>
// constexpr auto build_tuple_of_combination_functions_aux(const Tuple& tuple,
//                                                         const L2DotProductIntegral<Left,Right,VolumeIntegral,QR>& l2prod)
// {
//   using L2=L2DotProductIntegral<Left,Right,VolumeIntegral,QR>;
//   using QRule=typename L2::QRule;
//   auto tuple_new=build_tuple_of_combination_functions_aux_aux<QRule>(tuple,l2prod.left());
//   return build_tuple_of_combination_functions_aux_aux<QRule>(tuple_new,l2prod.right());
// }


// template<typename Tuple, typename Left,typename Right>
// constexpr auto build_tuple_of_combination_functions_aux(const Tuple& tuple,
//                           const Addition<Expression<Left>,Expression<Right>>& addition)
// {
//   auto tuple_new=build_tuple_of_combination_functions_aux(tuple,addition.left());
//   return build_tuple_of_combination_functions_aux(tuple_new,addition.right());
// }


// template<template<class,class>class AdditionOrSubtraction,typename Tuple, typename Left,typename Right>
// constexpr auto build_tuple_of_combination_functions_aux(const Tuple& tuple,
//                           const AdditionOrSubtraction<Expression<Left>,Expression<Right>>& addition)
// {
//   auto tuple_new=build_tuple_of_combination_functions_aux(tuple,addition.left());
//   return build_tuple_of_combination_functions_aux(tuple_new,addition.right());
// }




// template<typename Tuple,typename Form>
// constexpr auto build_tuple_of_combination_functions_form(const Tuple& tuple,const Form& form)
// {
//   return build_tuple_of_combination_functions_aux(tuple,form);
// }

// template<typename Tuple,typename Form,typename...Forms>
// constexpr auto build_tuple_of_combination_functions_form(const Tuple& tuple,const Form& form, const Forms&...forms)
// {
//   auto tuple_new= build_tuple_of_combination_functions_aux(tuple,form);
//   return build_tuple_of_combination_functions_form(tuple_new,forms...);
// }


// template<Integer Nmax,typename Form,typename...Forms>
// constexpr auto build_tuple_of_combination_functions(const Form& form, const Forms&...forms)
// {
//   using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
//   return build_tuple_of_combination_functions_form(emptytuple(),form,forms...);
// }















template<typename Tuple, typename Left,typename Right>
constexpr auto build_tuple_of_combination_functions_surface_aux(const Tuple& tuple,
                          const Addition<Expression<Left>,Expression<Right>>& addition)
{
  auto tuple_new=build_tuple_of_combination_functions_surface_aux(tuple,addition.left());
  return build_tuple_of_combination_functions_surface_aux(tuple_new,addition.right());
}


template<typename Tuple,typename Form>
constexpr auto build_tuple_of_combination_functions_surface_form(const Tuple& tuple,const Form& form)
{
  return build_tuple_of_combination_functions_surface_aux(tuple,form);
}

template<typename Tuple,typename Form,typename...Forms>
constexpr auto build_tuple_of_combination_functions_surface_form(const Tuple& tuple,const Form& form, const Forms&...forms)
{
  auto tuple_new= build_tuple_of_combination_functions_surface_aux(tuple,form);
  return build_tuple_of_combination_functions_surface_form(tuple_new,forms...);
}

template<Integer Nmax,typename Form,typename...Forms>
constexpr auto build_tuple_of_combination_functions_surface(const Form& form, const Forms&...forms)
{
  using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
  return build_tuple_of_combination_functions_surface_form(emptytuple(),form,forms...);
}







template<typename Tuple, typename Left,typename Right>
constexpr auto build_tuple_of_combination_functions_volume_aux(const Tuple& tuple,
                          const Addition<Expression<Left>,Expression<Right>>& addition)
{
  auto tuple_new=build_tuple_of_combination_functions_volume_aux(tuple,addition.left());
  return build_tuple_of_combination_functions_volume_aux(tuple_new,addition.right());
}


template<typename Tuple,typename Form>
constexpr auto build_tuple_of_combination_functions_volume_form(const Tuple& tuple,const Form& form)
{
  return build_tuple_of_combination_functions_volume_aux(tuple,form);
}

template<typename Tuple,typename Form,typename...Forms>
constexpr auto build_tuple_of_combination_functions_volume_form(const Tuple& tuple,const Form& form, const Forms&...forms)
{
  auto tuple_new= build_tuple_of_combination_functions_volume_aux(tuple,form);
  return build_tuple_of_combination_functions_volume_form(tuple_new,forms...);
}

template<Integer Nmax,typename Form,typename...Forms>
constexpr auto build_tuple_of_combination_functions_volume(const Form& form, const Forms&...forms)
{
  using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
  return build_tuple_of_combination_functions_volume_form(emptytuple(),form,forms...);
}





template<typename...Ts>
class ShapeFunctionsCollection;



template<typename Form_, typename...Forms_>
class ShapeFunctionsCollection<GeneralForm<Form_>,GeneralForm<Forms_>...>
{
 public:
  using GenericForm=GeneralForm<Form_>;
  using Elem=typename GenericForm::FunctionSpace::Elem;
  using Form=MultipleAddition<typename GeneralForm<Form_>::Form,typename GeneralForm<Forms_>::Form...> ;  
  using VolumetricForm=OperatorType<ExtractFormType<0,Form>>;
  using SurfaceForm=OperatorType<ExtractFormType<1,Form>>;
  
  using UniqueElementFunctionSpacesTupleType=typename GenericForm::UniqueElementFunctionSpacesTupleType;
  // using TupleOperatorsAndQuadrature= typename OperatorAndQuadratureTupleType<Form>::type;
  // using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadrature>;
  // H=0, we consider only volume shape functions
  using TupleOperatorsAndQuadrature= typename OperatorAndQuadratureTupleType<Form>::type;
  using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadrature>;

  using TupleOperatorsAndQuadratureVolumetric= typename OperatorAndQuadratureTupleType<Form,0>::type;

  using TupleOperatorsAndQuadratureSurface= typename OperatorAndQuadratureTupleType<Form,1>::type;



  using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;

  using Map=MapOperatorTupleOfTuple<TupleOfTupleNoQuadrature,UniqueElementFunctionSpacesTupleType>;
  using UniqueMapping=UniqueMap<SpacesToUniqueFEFamilies,Map> ;
  using TupleOfTupleShapeFunctionVolumetric=TupleOfTupleShapeFunctionType2<UniqueElementFunctionSpacesTupleType,TupleOperatorsAndQuadratureVolumetric>;
  using TupleOfTupleShapeFunctionSurface=TupleOfTupleShapeFunctionType2<UniqueElementFunctionSpacesTupleType,TupleOperatorsAndQuadratureSurface>;
  
  // using TupleCompositeOperatorsAndQuadrature= typename OperatorAndQuadratureTupleType<Form>::composite_type;
  // H=0, we consider only volume shape functions
  using TupleCompositeOperatorsAndQuadratureVolumetric= typename OperatorAndQuadratureTupleType<Form,0>::composite_type;
  using MapTupleNumbersVolumetric=MapTupleInit<TupleOfTupleShapeFunctionVolumetric, SpacesToUniqueFEFamilies, UniqueMapping>;
  using TupleOfTupleCompositeShapeFunctionVolumetric=typename TupleOfCombinationFunctions<GenericForm::FunctionSpace::Nuniquesubspaces,VolumetricForm>::type;
  using TupleOfTupleCompositeShapeFunctionTensorVolumetric=typename TupleOfCombinationFunctions<GenericForm::FunctionSpace::Nuniquesubspaces,VolumetricForm>::type_tensor;

  using TupleCompositeOperatorsAndQuadratureSurface= typename OperatorAndQuadratureTupleType<Form,1>::composite_type;
  using MapTupleNumbersSurface=MapTupleInit<TupleOfTupleShapeFunctionSurface, SpacesToUniqueFEFamilies, UniqueMapping>;
  using TupleOfTupleCompositeShapeFunctionSurface=typename TupleOfCombinationFunctions<GenericForm::FunctionSpace::Nuniquesubspaces,SurfaceForm>::type;
  using TupleOfTupleCompositeShapeFunctionTensorSurface=typename TupleOfCombinationFunctions<GenericForm::FunctionSpace::Nuniquesubspaces,SurfaceForm>::type_tensor;
 
  


  using CoefficientsCollection=ShapeFunctionCoefficientsCollection<GeneralForm<Form_>,GeneralForm<Forms_>...>;
  using MapCollection=MapFromReferenceCollection<GeneralForm<Form_>,GeneralForm<Forms_>...>;




   template<Integer M,Integer Nmax_aux,Integer N,typename Tuple,typename Map>
   constexpr typename std::enable_if< (N>Nmax_aux || IsSame<Tuple,std::tuple<>>::value) ,void>::type 
   init_map_aux_aux(Tuple& t,const Map& maps){}



   template<Integer M,Integer Nmax_aux,Integer N,typename Tuple,typename Map>
   constexpr typename std::enable_if< ((N<=Nmax_aux) && IsDifferent<Tuple,std::tuple<>>::value) ,void>::type 
   init_map_aux_aux(Tuple& t,const Map& maps) 
   {
    auto& t_nth=std::get<N>(t); 
    using T=remove_all_t<decltype(t_nth )>;
    using MMM=MapFromReference<typename T::Operator,typename T::QuadratureElem,T::FEFamily>;
    auto& map_nth=std::get<TypeToTupleElementPosition<MMM,Map>::value>(maps); 
    // auto& map_nth=std::get<GetType<MapTupleNumbersVolumetric,M,N>::value>(maps); 
    t_nth.init_map(map_nth);
    init_map_aux_aux<M,Nmax_aux,N+1>(t,maps);
    }



   template<Integer Nmax_aux,Integer N,typename Tuple,typename Map>
   constexpr typename std::enable_if< (N>Nmax_aux) ,void>::type 
   init_map_aux(Tuple& tuple,const Map& maps){}

   template<Integer Nmax_aux,Integer N,typename Tuple,typename Map>
   constexpr typename std::enable_if< (N<=Nmax_aux) ,void>::type 
   init_map_aux(Tuple& tuple, const Map& maps) 
   {
    auto& t_nth=std::get<N>(tuple); 
    auto& map_nth=std::get<GetType<SpacesToUniqueFEFamilies,N>::value>(maps); 
    init_map_aux_aux<N,TupleTypeSize<decltype(t_nth)>::value-1,0>(t_nth,map_nth);
    init_map_aux<Nmax_aux,N+1>(tuple,maps);
    }


   template<typename MapTupleNumbers,typename...Args>
   constexpr void init_map_volume(const MapFromReferenceCollection<Args...>& maps) 
   {init_map_aux<TupleTypeSize<MapTupleNumbers>::value-1,0>(tuple_volume_,maps.volumetric_map());}

   template<typename MapTupleNumbers,typename...Args>
   constexpr void init_map_surface(const MapFromReferenceCollection<Args...>& maps) 
   {init_map_aux<TupleTypeSize<MapTupleNumbers>::value-1,0>(tuple_surface_,maps.surface_map());}






  template<Integer N,Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
  constexpr typename std::enable_if_t< IsSame<Shape,std::tuple<>>::value,void> 
  shape_function_init_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficientsCollection<Args...> &coefficients,FiniteElem<Elem>&FE)
  {}


  template<Integer N,Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
  constexpr typename std::enable_if_t< FEFamily==LagrangeFE,void> 
  shape_function_init_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficientsCollection<Args...> &coefficients,FiniteElem<Elem>&FE)
  {
    // std::cout<<" M = "<< M<<" FEFamily = "<<FEFamily<<" Order = "<<Order<<std::endl;
    tuple_get<N>(shape).init(FE);
    // std::cout<<" map = "<< tuple_get<N>(shape).map()()<<std::endl;
    // std::cout<<" reference_values = "<< tuple_get<N>(shape).reference_values<<std::endl;
    // std::cout<<tuple_get<N>(shape).eval()<<std::endl;
  }


  template<Integer N,Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
  constexpr typename std::enable_if_t<FEFamily==RaviartThomasFE,void> 
  shape_function_init_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficientsCollection<Args...> &coefficients,FiniteElem<Elem>&FE)
  {
    // std::cout<<" M = "<< M<<" FEFamily = "<<FEFamily<<" Order = "<<Order<<std::endl;
    tuple_get<N>(shape).init(tuple_get<M>(coefficients()),FE);
    // std::cout<<" map = "<< tuple_get<N>(shape).map()()<<std::endl;
    // std::cout<<" reference_values = "<< tuple_get<N>(shape).reference_values<<std::endl;
    // std::cout<<tuple_get<N>(shape).eval()<<std::endl;
    // std::cout<<"coefficients="<<std::endl;
    // std::cout<<tuple_get<M>(coefficients())<<std::endl;
  }



  template<typename Shape,Integer N,Integer M,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t< IsSame<Shape,std::tuple<>>::value,void> 
  shape_function_init_aux_aux_helper(Tuple& tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients,FiniteElem<Elem>&FE)
  {}


  template<typename Shape,Integer N,Integer M,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t< IsDifferent<Shape,std::tuple<>>::value,void> 
  shape_function_init_aux_aux_helper(Tuple& tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients,FiniteElem<Elem>&FE)
  {
    using Elem=typename Shape::Elem;
    constexpr Integer FEFamily=Shape::FEFamily;
    constexpr Integer Order=Shape::Order;
    shape_function_init_aux_aux_aux<N,M,Elem,FEFamily,Order>(tuple,coefficients,FE);
  }



  template<Integer M, Integer Nmax,Integer N,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t<(N>Nmax),void> 
  shape_function_init_aux_aux(Tuple& tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients,FiniteElem<Elem>&FE)
  {}

  template<Integer M,Integer Nmax,Integer N,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t< (N<=Nmax),void> 
  shape_function_init_aux_aux(Tuple& tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients,FiniteElem<Elem>&FE)
  {

    using Shape=GetType<TupleOfTupleShapeFunctionVolumetric,M,N>;
    // using Elem=typename Shape::Elem;
    // constexpr Integer FEFamily=Shape::FEFamily;
    // constexpr Integer Order=Shape::Order;
    // shape_function_init_aux_aux_aux<N,M,Elem,FEFamily,Order>(tuple,coefficients);
    shape_function_init_aux_aux_helper<Shape,N,M>(tuple,coefficients,FE);
    shape_function_init_aux_aux<M,Nmax,N+1>(tuple,coefficients,FE);
  }

  template<typename TupleOfTupleShapeFunction,Integer Nmax,Integer N,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t< (N>Nmax),void> 
  shape_function_init_aux(Tuple&tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients,FiniteElem<Elem>&FE)
  {}

  template<typename TupleOfTupleShapeFunction,Integer Nmax,Integer N,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t<N<=Nmax,void> 
  shape_function_init_aux(Tuple&tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients,FiniteElem<Elem>&FE)
  {
    constexpr Integer Nmax_aux=TupleTypeSize<GetType<TupleOfTupleShapeFunctionVolumetric,N>>::value-1;
    shape_function_init_aux_aux<N,Nmax_aux,0>(tuple_get<N>(tuple),coefficients,FE);
    shape_function_init_aux<TupleOfTupleShapeFunction,Nmax,N+1>(tuple,coefficients,FE);
  }



  //template<typename Coefficients>
  constexpr void init_shape_functions_volume(FiniteElem<Elem>&FE)//const Coefficients& shape_coefficients)
  {
   constexpr Integer Nmax=TupleTypeSize<TupleOfTupleShapeFunctionVolumetric>::value-1;
   shape_function_init_aux<TupleOfTupleShapeFunctionVolumetric,Nmax,0>(tuple_volume_,coeffs_,FE);
   }













  template<Integer N,Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
  constexpr typename std::enable_if_t< IsSame<Shape,std::tuple<>>::value,void> 
  shape_function_init_surface_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficientsCollection<Args...> &coefficients)
  {}


  template<Integer N,Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
  constexpr typename std::enable_if_t< FEFamily==LagrangeFE,void> 
  shape_function_init_surface_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficientsCollection<Args...> &coefficients, FiniteElem<Elem> &FE)
  {
    // std::cout<<"M,N= "<<M<<" "<<N<<std::endl;
    // from reference to actual element: use rectangular jacobian!!!
    tuple_get<N>(shape).init(FE.side_id(),FE);
    // std::cout<<tuple_get<N>(shape).eval()<<std::endl;
  }


  template<Integer N,Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
  constexpr typename std::enable_if_t<FEFamily==RaviartThomasFE,void> 
  shape_function_init_surface_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficientsCollection<Args...> &coefficients, FiniteElem<Elem> &FE)
  {
    // std::cout<<" shape_function_init_surface_aux_aux_aux "<<std::endl;
    // std::cout<<" "<<M<<" "<<N<<std::endl;
    std::cout<<tuple_get<M>(coefficients())<<std::endl;
    std::cout<<FE.side_id()<<std::endl;
    tuple_get<N>(shape).init(tuple_get<M>(coefficients()), FE.side_id(),FE);
  }



  template<Integer M, Integer Nmax,Integer N,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t<(N>Nmax || IsSame<Tuple,std::tuple<>>::value),void> 
  shape_function_init_surface_aux_aux(Tuple& tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients, FiniteElem<Elem> &FE)
  {}

  template<Integer M,Integer Nmax,Integer N,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t< (N<=Nmax && IsDifferent<Tuple,std::tuple<>>::value),void> 
  shape_function_init_surface_aux_aux(Tuple& tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients, FiniteElem<Elem> &FE)
  {

    using Shape=GetType<TupleOfTupleShapeFunctionSurface,M,N>;
    using Elem=typename Shape::Elem;
    constexpr Integer FEFamily=Shape::FEFamily;
    constexpr Integer Order=Shape::Order;
    // std::cout<<" shape_function_init_surface_aux_aux "<<N<<std::endl;
    shape_function_init_surface_aux_aux_aux<N,M,Elem,FEFamily,Order>(tuple,coefficients,FE);
    shape_function_init_surface_aux_aux<M,Nmax,N+1>(tuple,coefficients,FE);
  }


  template<typename TupleOfTupleShapeFunction,Integer Nmax,Integer N,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t< (N>Nmax),void> 
  shape_function_init_surface_aux(Tuple&tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients, FiniteElem<Elem> &FE)
  {}

  template<typename TupleOfTupleShapeFunction,Integer Nmax,Integer N,typename Tuple,typename...Args>
  constexpr typename std::enable_if_t<N<=Nmax,void> 
  shape_function_init_surface_aux(Tuple&tuple, const ShapeFunctionCoefficientsCollection<Args...> &coefficients, FiniteElem<Elem> &FE)
  {
    constexpr Integer Nmax_aux=TupleTypeSize<GetType<TupleOfTupleShapeFunctionSurface,N>>::value-1;
    // std::cout<<" shape_function_init_surface_aux "<<N<<std::endl;
    shape_function_init_surface_aux_aux<N,Nmax_aux,0>(tuple_get<N>(tuple),coefficients,FE);
    shape_function_init_surface_aux<TupleOfTupleShapeFunction,Nmax,N+1>(tuple,coefficients,FE);
  }



  //template<typename Coefficients>
  constexpr void init_shape_functions_surface(FiniteElem<Elem> &FE)//const Coefficients& shape_coefficients)
  {
   constexpr Integer Nmax=TupleTypeSize<TupleOfTupleShapeFunctionSurface>::value-1;
   shape_function_init_surface_aux<TupleOfTupleShapeFunctionSurface,Nmax,0>(tuple_surface_,coeffs_,FE);
   }

















  template<Integer Nmax,Integer N,typename Tuple,typename...Args1,typename...Args2>
  constexpr typename std::enable_if_t<(N>Nmax),void> 
  init_composite_shape_functions_aux_aux(Tuple& tuple, std::tuple<Args1...>& tuple_tensor, FiniteElem<Elem> &J, std::tuple<Args2...>& tuple_composite)
  {
    // std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;
  }

  template<Integer Nmax,Integer N,typename Tuple,typename...Args1>
  constexpr void init_composite_shape_functions_aux_aux(Tuple& tuple, std::tuple<>& tuple_tensor, FiniteElem<Elem> &J, std::tuple<Args1...>& tuple_composite)
  {
    // std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;
  }

  template<Integer Nmax,Integer N,typename Tuple,typename...Args1>
  constexpr void init_composite_shape_functions_aux_aux(Tuple& tuple, std::tuple<Args1...>& tuple_tensor, FiniteElem<Elem> &J, std::tuple<>& tuple_composite)
  {
    // std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;

  }

  template<Integer Nmax,Integer N,typename Tuple>
  constexpr void init_composite_shape_functions_aux_aux(Tuple& tuple, std::tuple<>& tuple_tensor, FiniteElem<Elem> &J, std::tuple<>& tuple_composite)
  {
    // std::cout<<"init_composite_shape_functions_aux_aux Nmax,N="<<Nmax<<", "<<N<<std::endl;
  }

  template<Integer Nmax,Integer N,typename Tuple,typename...Args1,typename...Args2>
  constexpr typename std::enable_if_t< (N<=Nmax),void> 
  init_composite_shape_functions_aux_aux(Tuple& tuple, std::tuple<Args1...>& tuple_tensor, FiniteElem<Elem> &J, std::tuple<Args2...>& tuple_composite)
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
    // std::cout<<"init_composite_shape_functions_aux_aux Nmax,N="<<Nmax<<", "<<N<<std::endl;

    // todo fixme scommentami
    composite.apply(tensor,J,tuple);
    init_composite_shape_functions_aux_aux<Nmax,N+1>(tuple,tuple_tensor,J,tuple_composite);
  }


  template<Integer Nmax,Integer N,typename Tuple, typename TupleTensor,typename TupleComposite>
  constexpr typename std::enable_if_t< (N>Nmax),void> 
  init_composite_shape_functions_aux(Tuple&tuple, TupleTensor&tuple_tensor,TupleComposite& tuple_composite, FiniteElem<Elem> &J)
  {}

  template<Integer Nmax,Integer N,typename Tuple, typename TupleTensor,typename TupleComposite>
  constexpr typename std::enable_if_t<N<=Nmax,void> 
  init_composite_shape_functions_aux(Tuple&tuple,TupleTensor&tuple_tensor,TupleComposite& tuple_composite, FiniteElem<Elem> &J)
  {
    auto& tensor=tuple_get<N>(tuple_tensor);
    const auto& composite=tuple_get<N>(tuple_composite);
    constexpr Integer Nmax_aux=TupleTypeSize<decltype(composite)>::value-1;
    constexpr Integer Nmax_au2=TupleTypeSize<decltype(tensor)>::value-1;
    init_composite_shape_functions_aux_aux<Nmax_aux,0>(tuple,tuple_get<N>(tuple_tensor),J,tuple_get<N>(tuple_composite));
    init_composite_shape_functions_aux<Nmax,N+1>(tuple,tuple_tensor,tuple_composite,J);
  }

  constexpr void init_composite_shape_functions_volume( FiniteElem<Elem> &J)
  {
    // std::cout<<"init_composite_shape_functions"<<std::endl;
   constexpr Integer Nmax=TupleTypeSize<TupleOfTupleCompositeShapeFunctionVolumetric>::value-1;
   init_composite_shape_functions_aux<Nmax,0>(tuple_volume_,tuple_tensors_volume_,tuple_composite_volume_,J);
   }


  constexpr void init_composite_shape_functions_surface( FiniteElem<Elem> &J)
  {
    // std::cout<<"init_composite_shape_functions"<<std::endl;
   constexpr Integer Nmax=TupleTypeSize<TupleOfTupleCompositeShapeFunctionSurface>::value-1;
   init_composite_shape_functions_aux<Nmax,0>(tuple_surface_,tuple_tensors_surface_,tuple_composite_surface_,J);
   }





  constexpr void init(FiniteElem<Elem>&FE)
  {
    // TODO CHECK: SINCE WE HAVE ALREADY HAVE REFERENCES TO MAPS AND COEFFS, do we have to init_map?
    // INIT MAP: takes the corresponding map and compute the shape function in the actual element
    // std::cout<<"init maps"<<std::endl;
    init_map_volume<MapTupleNumbersVolumetric>(maps_);
    // init_shape_functions: takes also the coefficients which multiply the actual shape functions
    // std::cout<<"init_shape_functions"<<std::endl;
    init_shape_functions_volume(FE);
    // std::cout<<"init_composite_shape_functions"<<std::endl;
    init_composite_shape_functions_volume(FE);
   }

  constexpr void init_boundary(FiniteElem<Elem>&FE)
  {
    // TODO CHECK: SINCE WE HAVE ALREADY HAVE REFERENCES TO MAPS AND COEFFS, do we have to init_map?
    // INIT MAP: takes the corresponding map and compute the shape function in the actual element
    // std::cout<<"init_boundary init maps"<<std::endl;
    init_map_surface<MapTupleNumbersSurface>(maps_);
    // init_shape_functions: takes also the coefficients which multiply the actual shape functions
    // std::cout<<"init_boundary init_shape_functions"<<std::endl;
    init_shape_functions_surface(FE);
    // std::cout<<"init_boundary init_composite_shape_functions"<<std::endl;
    init_composite_shape_functions_surface(FE);
   }


    
   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==true,TupleOfTupleShapeFunctionVolumetric &> 
   tuple(){return tuple_volume_;}

   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==false,TupleOfTupleShapeFunctionSurface &> 
   tuple(){return tuple_surface_;}

   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==true,const TupleOfTupleShapeFunctionVolumetric &> 
   tuple()const{return tuple_volume_;}

   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==false,const TupleOfTupleShapeFunctionSurface &> 
   tuple()const{return tuple_surface_;}




   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==true,TupleOfTupleCompositeShapeFunctionVolumetric &> 
   composite_shapes(){return tuple_composite_volume_;}

   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==false,TupleOfTupleCompositeShapeFunctionSurface &> 
   composite_shapes(){return tuple_composite_surface_;}

   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==true,const TupleOfTupleCompositeShapeFunctionVolumetric &> 
   composite_shapes()const{return tuple_composite_volume_;}

   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==false,const TupleOfTupleCompositeShapeFunctionSurface &> 
   composite_shapes()const{return tuple_composite_surface_;}



   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==true,TupleOfTupleCompositeShapeFunctionTensorVolumetric &> 
   composite_tensor(){return tuple_tensors_volume_;}

   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==false,TupleOfTupleCompositeShapeFunctionTensorSurface &> 
   composite_tensor(){return tuple_tensors_surface_;}

   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==true,const TupleOfTupleCompositeShapeFunctionTensorVolumetric &> 
   composite_tensor()const{return tuple_tensors_volume_;}

   template<bool VolumeIntegral> 
   constexpr std::enable_if_t<VolumeIntegral==false,const TupleOfTupleCompositeShapeFunctionTensorSurface &> 
   composite_tensor()const{return tuple_tensors_surface_;}
   // constexpr       auto& volumetric_tuple()      {return tuple_volume_;}
   // constexpr const auto& volumetric_tuple()const {return tuple_volume_;}

   // constexpr       auto& surface_tuple()      {return tuple_surface_;}
   // constexpr const auto& surface_tuple()const {return tuple_surface_;}


   // constexpr       auto& operator()()      {return tuple;}
   // constexpr const auto& operator()()const {return tuple;}


   // constexpr       auto & composite_shapes_volume()      {return tuple_composite_volume_;}
   // constexpr const auto & composite_shapes_volume()const {return tuple_composite_volume_;}

   // constexpr       auto & composite_shapes_surface()      {return tuple_composite_surface_;}
   // constexpr const auto & composite_shapes_surface()const {return tuple_composite_surface_;}

   // constexpr       auto & composite_tensor_volume()      {return tuple_tensors_volume_;}
   // constexpr const auto & composite_tensor_volume()const {return tuple_tensors_volume_;}


   template<Integer...Ns>
   constexpr const auto& get_volume()const{return tuple_get<Ns...>(tuple_volume_);}

   template<Integer...Ns>
         auto& get_volume()     {return tuple_get<Ns...>(tuple_volume_);}

   template<Integer...Ns>
   constexpr const auto& value_volume()const{return tuple_get<Ns...>(tuple_volume_).eval();}

   template<Integer...Ns>
         auto& value_volume()     {return tuple_get<Ns...>(tuple_volume_).eval();}


   template<Integer...Ns>
   constexpr const auto& get_surface()const{return tuple_get<Ns...>(tuple_surface_);}

   template<Integer...Ns>
         auto& get_surface()     {return tuple_get<Ns...>(tuple_surface_);}

   template<Integer...Ns>
   constexpr const auto& value_surface()const{return tuple_get<Ns...>(tuple_surface_).eval();}

   template<Integer...Ns>
         auto& value_surface()     {return tuple_get<Ns...>(tuple_surface_).eval();}



 ShapeFunctionsCollection(ShapeFunctionCoefficientsCollection<GeneralForm<Form_>,GeneralForm<Forms_>...>&coeffs,
                 MapFromReferenceCollection<GeneralForm<Form_>,GeneralForm<Forms_>...>& maps,
                 const GeneralForm<Form_>& form,const GeneralForm<Forms_>&...forms)
 :
 coeffs_(coeffs),
 maps_(maps)
 ,
 tuple_composite_volume_(build_tuple_of_combination_functions_volume<GenericForm::FunctionSpace::Nuniquesubspaces>(form(),forms()...)),
 tuple_composite_surface_(build_tuple_of_combination_functions_surface<GenericForm::FunctionSpace::Nuniquesubspaces>(form(),forms()...))
 { }


private:
   CoefficientsCollection & coeffs_;
   MapCollection & maps_;
   TupleOfTupleShapeFunctionVolumetric tuple_volume_;
   TupleOfTupleCompositeShapeFunctionVolumetric tuple_composite_volume_;
   TupleOfTupleCompositeShapeFunctionTensorVolumetric tuple_tensors_volume_;

   TupleOfTupleShapeFunctionSurface tuple_surface_;
   TupleOfTupleCompositeShapeFunctionSurface tuple_composite_surface_;
   TupleOfTupleCompositeShapeFunctionTensorSurface tuple_tensors_surface_;
};


template<typename...Forms_>
constexpr auto shape_functions(ShapeFunctionCoefficientsCollection<GeneralForm<Forms_>...>&coeffs, 
                               MapFromReferenceCollection<GeneralForm<Forms_>...>&maps,
                               const GeneralForm<Forms_>&...forms)
{
  //using Form=typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type;
 return ShapeFunctionsCollection<GeneralForm<Forms_>...>(coeffs,maps,forms...);  }
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




template<typename...BCs>
class ShapeFunctionsCollection//<BCs...>
{
 public:

};

}
#endif