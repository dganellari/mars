/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                             ////////
////// We define the FunctionSpace class:                                                                                ////////
////// 1) It takes a mesh and 1 or more FunctionSpaces (Lagrange1<2>, RT0<1>...)                                         ////////                                        
////// 2) It builds the dofmap: a vector (long n_elements), whose component is the array of all the dofs of the element  ////////
////// 2) dofmap(space_id,elem_id) returns the dofs of element elem_id corresponding to the space space_id               ////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MARS_FunctionSpace_HPP
#define MARS_FunctionSpace_HPP

#include "mars_base.hpp"
#include "mars_elementfunctionspace.hpp"
#include "mars_functionspace_dofmap.hpp"
#include "mars_shape_function.hpp"
#include "mars_tuple_utilities.hpp"
#include "mars_operators.hpp"
#include "mars_quadrature_rules.hpp"
#include "mars_vector.hpp"
#include "mars_number.hpp"
namespace mars{



template<typename MeshT, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
class FunctionSpace 
{
public:

      using Elem= typename MeshT::Elem;

      static constexpr Integer Nsubspaces=1+sizeof...(BaseFunctionSpaces);

      static constexpr Integer Nelem_dofs=DofsPerElemNums<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;

      using DofMapType=std::vector<std::array<Integer, Nelem_dofs>>;
      using OffSetType=std::array<std::vector<Integer>, Nsubspaces>;
      using SpacesDofsArrayType=std::array<std::vector<std::vector<Integer>>, Nsubspaces>;
      using SpacesInfosArrayType=std::array<std::array<Integer,4>,Nsubspaces>;
      using ElemsTupleType=std::tuple<Elem>;
      using ElementFunctionSpacesTupleType=std::tuple<std::tuple<Elem,BaseFunctionSpace>,std::tuple<Elem,BaseFunctionSpaces>...>;
      using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates<ElementFunctionSpacesTupleType>;
      using SpacesToUniqueFEFamily=SpacesToUniqueFEFamilies<UniqueElementFunctionSpacesTupleType>;
      using FromElementFunctionSpacesToUniqueNumbersTupleType=TupleAllToUniqueMap<ElementFunctionSpacesTupleType,UniqueElementFunctionSpacesTupleType>;
      inline Integer n_subspaces()const{return Nsubspaces;};

      inline const Integer& components (const Integer& space_id)const{return space_infos_[space_id][3];};

      inline const Integer& n_elem_dofs()const{return Nelem_dofs;};

      inline Integer n_elem_dofs(const Integer& space_id)const{
                                  const auto& os=offset_[space_id];
                                  const auto size=os[os.size()-1]-os[0];
                                  return size;}

      inline Integer n_elem_dofs(const Integer& space_id,const Integer& component_id)const{
                                  const auto& size=n_elem_dofs(space_id);
                                  return (size/space_infos_[space_id][3]);}


      inline const Integer& n_dofs()const{return n_dofs_;};

      inline const Integer n_dofs(const Integer& space_id,const Integer& component_id)const
                                 {return space_dofs_[space_id][component_id].size(); };

      inline const DofMapType& dofmap()const{return dofmap_;};

      inline void  dofmap(const DofMapType& dm)const{dm=dofmap_;};


      inline const std::array<Integer, Nelem_dofs>& dofmap(const Integer& elem_id)const
                         {return dofmap_[elem_id];};

      inline void  dofmap(const Integer& elem_id, const std::array<Integer, Nelem_dofs> & elem_dm)const
                         {elem_dm=dofmap_[elem_id];};


      inline std::vector<Integer> 
                   dofmap(const Integer& space_id,const Integer& elem_id)const{
                        const auto& os=offset_[space_id];
                        const auto& size=n_elem_dofs(space_id);
                        std::vector<Integer> output(size);
                        for(Integer nn=0;nn<size;nn++)
                             output[nn]=dofmap_[elem_id][nn+os[0]];
                        return output;};


      inline std::vector<Integer> 
                   dofmap(const Integer& space_id,const Integer& component_id,const Integer& elem_id)const{
                        const auto& comp=components(space_id);
                        const auto& os=offset_[space_id];
                        const auto& size= n_elem_dofs(space_id);
                        std::vector<Integer> output(size/comp);
                        space_infos_[space_id][3];
                        Integer mm=0;
                        for(Integer nn=component_id;nn<size;nn=nn+comp)
                             {output[mm]=dofmap_[elem_id][nn+os[0]];mm++;}
                        return output;};


      inline void dofmap(const Integer& space_id,const Integer& component_id,const Integer& elem_id, const std::vector<Integer>& elem_space_dm)const{
                        const auto& comp=components(space_id);
                        const auto& os=offset_[space_id];
                        const auto& size= n_elem_dofs(space_id);
                        std::vector<Integer> output(size/comp);
                        space_infos_[space_id][3];
                        Integer mm=0;
                        elem_space_dm.resize(size/comp);
                        for(Integer nn=component_id;nn<size;nn=nn+comp)
                             {elem_space_dm[mm]=dofmap_[elem_id][nn+os[0]];mm++;}
                        };

      inline const std::array<std::vector<Integer>, Nsubspaces>& offset() const {return offset_;};

      inline void  offset(const std::array<std::vector<Integer>, Nsubspaces> &os)const {os=offset_;};


      inline const std::vector<Integer>& offset(const Integer& space_id)const{return offset_[space_id];};

      inline void offset(Integer space_id, const std::vector<Integer>& space_os)const {space_os=offset_[space_id];};


      inline const std::vector<Integer>& space_dofs(const Integer& space_id,const Integer& component_id) const
                                         {return space_dofs_[space_id][component_id];};


      inline void space_dofs(const Integer& space_id, const Integer& component_id,std::vector<Integer>& spacedofs)const
                            {spacedofs.resize(n_dofs(space_id,component_id));
                             spacedofs=space_dofs_[space_id][component_id];};

      inline const SpacesInfosArrayType& space_info()const{return space_infos_;};

      inline std::shared_ptr< MeshT > mesh()const {return mesh_;};
       
      FunctionSpace(const MeshT& mesh,const Integer dof_count_start=0):
      mesh_(std::make_shared< MeshT >(mesh))
      {
      function_space_info<0,Nsubspaces,BaseFunctionSpace,BaseFunctionSpaces...>(space_infos_);
      dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap_,offset_,n_dofs_,space_infos_,space_dofs_);     
      };

private:
   
      std::shared_ptr< MeshT > mesh_;
      Integer n_dofs_;
      DofMapType dofmap_;
      OffSetType offset_;
      SpacesDofsArrayType space_dofs_;
      SpacesInfosArrayType space_infos_;
      ElementFunctionSpacesTupleType shape_functions_;

};





template<typename MeshT, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
class TraceSpace 
{

};





template<typename...Args>
class MixedSpace 
{
public:
  using DofMapType=std::tuple<typename Args::DofMapType...>;
  using ElementFunctionSpacesTupleType=TupleCatType<typename Args::ElementFunctionSpacesTupleType...>;
  using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates<TupleCatType<typename Args::UniqueElementFunctionSpacesTupleType...>>;
  using SpacesToUniqueFEFamily=SpacesToUniqueFEFamilies<UniqueElementFunctionSpacesTupleType>;
  using ElemsTupleType=RemoveTupleDuplicates<TupleCatType<typename Args::ElemsTupleType...>>;
  using FromElementFunctionSpacesToUniqueNumbersTupleType=TupleAllToUniqueMap<ElementFunctionSpacesTupleType,UniqueElementFunctionSpacesTupleType>;
  static constexpr Integer Nelem_dofs=Sum(Args::Nelem_dofs...);
  inline const Integer& n_dofs()const{return n_dofs_;}; 

  inline const DofMapType& dofmap()const{return dofmap_;};
  
  template<Integer...Ns>
  inline constexpr const auto& dofmap()const{return tuple_get<Ns...>(dofmap_);};


template<typename OtherArg,typename...OtherArgs>
  class tot_subspaces
  { public: static constexpr Integer value= OtherArg::Nsubspaces+tot_subspaces<OtherArgs...>::value;};

template<typename OtherArg>
  class tot_subspaces<OtherArg>
  { public:  static constexpr Integer value= OtherArg::Nsubspaces;};


template<Integer N>
  typename std::enable_if< 0==N, Integer >::type
  tot_n_dofs()
  { const auto& tmp=std::get<N>(spaces_);
    return tmp->n_dofs();}


template<Integer N>
    typename std::enable_if< 0<N, Integer >::type
    tot_n_dofs()
    { 
      static_assert(N>0, " tuple cannot have negative components");
      const auto& tmp=std::get<N>(spaces_);
      return tmp->n_dofs()+tot_n_dofs<N-1>();}



template<typename OtherArg,typename...OtherArgs>
      struct tuple_type
      {
        using rest = typename tuple_type<OtherArgs...>::type; 
        using tuple_ens=std::tuple<OtherArg>;
        using type = decltype( std::tuple_cat( std::declval< tuple_ens >(), std::declval< rest >() ) );
      };


template<typename OtherArg>
      struct tuple_type<OtherArg>
      {
       using type = typename std::tuple<OtherArg>;
     };


template<Integer N,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< 0==sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
     tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     { std::cout<<"N=="<<N<<". "<<tot_n_dofs<N-1>()<<std::endl;
     return std::tuple<OtherArg>(add_costant(otherarg,tot_n_dofs<N-1>()));}


template<Integer N,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< 0<N && 0<sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
     tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     {std::cout<<"N=="<<N<<". "<<tot_n_dofs<N-1>()<<std::endl;
     return std::tuple_cat(std::tuple<OtherArg>( add_costant(otherarg,tot_n_dofs<N-1>()) ),
       tuple_make<N+1,OtherArgs...>(otherargs...));}

template<Integer N,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< 0==N && 0<sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
     tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     {std::cout<<"N=="<<N<<". "<<0<<std::endl;
     return std::tuple_cat(std::tuple<OtherArg>( otherarg ),
       tuple_make<N+1,OtherArgs...>(otherargs...));}

     MixedSpace(const Args&...args):
     spaces_(std::make_tuple(std::make_shared<Args>(args)...)),
     n_dofs_(tot_n_dofs<sizeof...(Args)-1>()),
     dofmap_(tuple_make<0,typename Args::DofMapType...>(args.dofmap()...))
     {}

     static constexpr Integer Nsubspaces=tot_subspaces<Args...>::value;
   private:
    std::tuple<std::shared_ptr<Args>...> spaces_;
    Integer n_dofs_;
    DofMapType dofmap_;
  };

template<typename...Args>
MixedSpace<Args...> MixedFunctionSpace(const Args&...args){return MixedSpace<Args...>(args...);};


template<typename MixedSpace,Integer N, typename OperatorType=IdentityOperator>
class Trial: public Expression2<Trial<MixedSpace,N,OperatorType>>
{ public:
  using FunctionSpace=MixedSpace;
  using UniqueElementFunctionSpacesTupleType=typename MixedSpace::UniqueElementFunctionSpacesTupleType;
  static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
  static constexpr Integer value=N;
  using Operator=OperatorType;
};

template<typename MixedSpace, Integer N, typename OperatorType=IdentityOperator>
class Test: public Expression2<Test<MixedSpace,N,OperatorType>>
{
public:
  using FunctionSpace=MixedSpace;
  using UniqueElementFunctionSpacesTupleType=typename MixedSpace::UniqueElementFunctionSpacesTupleType;
  static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
  static constexpr Integer value=N;
  using Operator=OperatorType;
};



template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, Integer N, typename OperatorType>
class Evaluation<Expression2<TestOrTrial<MixedSpace,N,OperatorType>>>
{
 public:
 using type= TestOrTrial<MixedSpace,N,OperatorType>;

 static_assert((IsSame<TestOrTrial<MixedSpace,N,OperatorType>,Test<MixedSpace,N,OperatorType>>::value ||
                IsSame<TestOrTrial<MixedSpace,N,OperatorType>,Trial<MixedSpace,N,OperatorType>>::value )
               && "In Evaluation<Expression2<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");
 using FunctionSpace=typename type::FunctionSpace;
 using FunctionSpaces=GetType<typename type::UniqueElementFunctionSpacesTupleType,type::value>;
 using Elem=GetType<FunctionSpaces,0>;
 using BaseFunctionSpace=GetType<FunctionSpaces,1>;
 using Operator=typename type::Operator;
 template<typename QRule>
 using value_type=typename ShapeFunctionDependent<Elem,BaseFunctionSpace,Operator,QRule>::type;

 template<typename ...Ts> 
 Evaluation(const type& expr):
 eval_(expr)
 {};
 
 template<typename QRule, typename ...Ts>
 constexpr void apply(const std::tuple<Ts...>& tuple1,value_type<QRule>& value)
 {
  using tuple_type=GetType<typename std::tuple<Ts...>,type::value>;
  const auto& tuple=std::get<type::value>(tuple1);
  constexpr Integer M=TypeToTupleElementPosition<ShapeFunctionDependent<Elem,BaseFunctionSpace,Operator,QRule>,tuple_type>::value;
  value=std::get<M>(tuple);
  // return value_;
 }
private:

 type eval_;
};












template<typename Derived>
class Evaluation<Expression2<UnaryPlus2< Expression2<Derived> > > >
{
 public:
 using type=UnaryPlus2<Expression2<Derived>>;
 using Eval=Evaluation<Expression2<Derived>>;

 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_(Eval(expr_()))
 {};
 
 template<Integer Rows,Integer Cols>
 void apply(const Matrix<Real,Rows,Cols>& mat)
 {
  eval_.apply(mat);
 }
private:

 type expr_;
 Eval eval_;
};


template<typename Derived>
class Evaluation<Expression2<UnaryMinus2< Expression2<Derived> > > >
{
 public:
 using type=UnaryMinus2<Expression2<Derived>>;
 using Eval=Evaluation<Expression2<Derived>>;

 Evaluation(){};
 
 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_(Eval(expr_()))
 {};
 
 template<Integer Rows,Integer Cols>
 void apply(const Matrix<Real,Rows,Cols>& mat)
 {
  eval_.apply(mat);
 }
private:

 type expr_;
 Eval eval_;
};




template<typename DerivedRight,typename DerivedLeft>
class Evaluation<Expression2<Addition2< Expression2<DerivedLeft>  ,  
                                        Expression2<DerivedRight>
                                        >>>
{
 public:
 using type=Addition2<Expression2<DerivedLeft>,Expression2<DerivedRight>>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>>;
 using EvalRight=Evaluation<Expression2<DerivedRight>>;

 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(EvalRight(expr_.right()))
{};
 
 template<Integer Rows,Integer Cols>
 void apply(const Matrix<Real,Rows,Cols>& mat)
 {
  eval_left_.apply(mat);
  eval_right_.apply(mat); 
 }
private:

 type expr_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
};

template<typename DerivedRight,typename DerivedLeft>
class Evaluation<Expression2<Subtraction2< Expression2<DerivedLeft>  ,  
                                        Expression2<DerivedRight>
                                        >>>
{
 public:
 using type=Subtraction2<Expression2<DerivedLeft>,Expression2<DerivedRight>>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>>;
 using EvalRight=Evaluation<Expression2<DerivedRight>>;

 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(EvalRight(expr_.right()))
{};
 
 template<Integer Rows,Integer Cols>
 void apply(const Matrix<Real,Rows,Cols>& mat)
 {
  eval_left_.apply(mat);
  eval_right_.apply(mat); 
 }
private:

 type expr_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
};


template<typename DerivedRight>
class Evaluation<Expression2<Multiply2< Real, Expression2<DerivedRight> >>>                                     
{
 public:
 using type=Multiply2<Real,Expression2<DerivedRight>>;
 using EvalRight=Evaluation<Expression2<DerivedRight>>;
 
 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(expr_.left()),
 eval_right_(EvalRight(expr_.right()))
{};
 
 template<Integer Rows,Integer Cols>
 void apply(const Matrix<Real,Rows,Cols>& mat)
 {
  // eval_left_.apply(mat);
  // eval_right_.apply(mat); 
 }
private:

 type expr_;
 Real eval_left_;
 EvalRight eval_right_;
};




template< typename DerivedLeft>
class Evaluation<Expression2<Multiply2< Expression2<DerivedLeft>,
                                        Real
                                        >>>                                       
{
 public:
 using type=Multiply2<Expression2<DerivedLeft>,Real>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>>;
 
 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(expr_.right())
{};
 
 template<Integer Rows,Integer Cols>
 void apply(const Matrix<Real,Rows,Cols>& mat)
 {
  // eval_left_.apply(mat);
  // eval_right_.apply(mat); 
 }
private:

 type expr_;
 EvalLeft eval_left_;
 Real eval_right_;
};

template< typename DerivedLeft>
class Evaluation<Expression2<Divide2< Expression2<DerivedLeft>,
                                        Real
                                        >>>                                       
{
 public:
 using type=Divide2<Expression2<DerivedLeft>,Real>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>>;
 
 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left())),
 eval_right_(expr_.right())
{};
 
 template<Integer Rows,Integer Cols>
 void apply(const Matrix<Real,Rows,Cols>& mat)
 {
  // eval_left_.apply(mat);
  // eval_right_.apply(mat); 
 }
private:

 type expr_;
 EvalLeft eval_left_;
 Real eval_right_;
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// We associate Number<2>,Number<1>,Number<0> respectively to a X=Trial/Test/OtherFunction
///// We define IsTestOrTrial<X> which is used for each term Left or Right in L2Inner<Left,Right>
///// The possible results are:
///// 1)IsTestOrTrial<X> == std::tuple<Number<0>> (if only OtherFunction is there, it is not touched)
///// 2)IsTestOrTrial<X> == std::tuple<Number<1>> (if OtherFunction is there, Number<0> is removed and only Number<1> survives)
///// 3)IsTestOrTrial<X> == std::tuple<Number<2>> (if OtherFunction is there, Number<0> is removed and only Number<2> survives)
///// 4)Error with static assert (if more than one Number<1>/Number<2> is present or if both Number<1> and Number<2> are present)
/////   This assert means that we can have only: a function or exactly one trial or exactly one test
///// Then is check the form of L2inner, bu using TypeOfForm<Number<M>,Number<N>>:
///// 0) 0 form if:
/////  Left==Right==std::tuple<Number<0>>
///// 1) 1 form if: 
/////              Left==std::tuple<Number<0>> and  Right==std::tuple<Number<1>> (function,test)
/////              Left==std::tuple<Number<1>> and  Right==std::tuple<Number<0>> (test,function)
///// 2) 2 form if Left==Right==Number<0>
/////              Left==std::tuple<Number<1>> and  Right==std::tuple<Number<2>> (test,trial)
/////              Left==std::tuple<Number<2>> and  Right==std::tuple<Number<1>> (trial,test)
///// 3) everything else: error
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
class IsTestOrTrial{
public:
  using TupleFunctionSpace=std::tuple<>;
  using UniqueElementFunctionSpacesTupleType=std::tuple<>;
  using type=std::tuple<Number<0>>;
};

template<typename MixedSpace,Integer N, typename OperatorType>
class IsTestOrTrial<Test<MixedSpace,N,OperatorType>>{
public:
  using TupleFunctionSpace=std::tuple<MixedSpace>;
  using UniqueElementFunctionSpacesTupleType=std::tuple<typename MixedSpace::UniqueElementFunctionSpacesTupleType>;
  using type=std::tuple<Number<1>>;
};

template<typename MixedSpace,Integer N, typename OperatorType>
class IsTestOrTrial<Trial<MixedSpace,N,OperatorType>>{
public:
  using TupleFunctionSpace=std::tuple<MixedSpace>;
  using UniqueElementFunctionSpacesTupleType=std::tuple<typename MixedSpace::UniqueElementFunctionSpacesTupleType>;
  using type=std::tuple<Number<2>>;
};


template <typename T,typename ... Types>
class TupleTypeSize;


template<typename Left, typename Right>
class IsTestOrTrial< Multiply2<Expression2<Left>,Expression2<Right> > >
{
public:
  using TupleFunctionSpace=TupleCatType<typename IsTestOrTrial<Left>::TupleFunctionSpace,
                                        typename IsTestOrTrial<Right>::TupleFunctionSpace >;

  using UniqueElementFunctionSpacesTupleType=TupleCatType<typename IsTestOrTrial<Left>::UniqueElementFunctionSpacesTupleType,
                                                          typename IsTestOrTrial<Right>::UniqueElementFunctionSpacesTupleType >;
  using type=TupleRemoveType<Number<0>,TupleCatType<typename IsTestOrTrial<Left>::type,typename IsTestOrTrial<Right>::type >> ;
  static_assert(TupleTypeSize<type>::value<2," In Multiply<Left,Right>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
};

template<typename Left>
class IsTestOrTrial< Multiply2<Expression2<Left>,Real > >
{
public:
  using TupleFunctionSpace=typename IsTestOrTrial<Left>::TupleFunctionSpace;
  using UniqueElementFunctionSpacesTupleType=typename IsTestOrTrial<Left>::UniqueElementFunctionSpacesTupleType;
  using type=typename IsTestOrTrial<Left>::type;
  static_assert(TupleTypeSize<type>::value<2," In Multiply<Left,Real>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
};

template<typename Right>
class IsTestOrTrial< Multiply2<Real,Expression2<Right> > >
{
public:
  using TupleFunctionSpace=typename IsTestOrTrial<Right>::TupleFunctionSpace;
  using UniqueElementFunctionSpacesTupleType=typename IsTestOrTrial<Right>::UniqueElementFunctionSpacesTupleType;
  using type=typename IsTestOrTrial<Right>::type;
  static_assert(TupleTypeSize<type>::value<2," In Multiply<Real,Right>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");

};


template<typename Left>
class IsTestOrTrial< Divide2<Expression2<Left>,Real > >
{
public:
  using TupleFunctionSpace=typename IsTestOrTrial<Left>::TupleFunctionSpace;
  using UniqueElementFunctionSpacesTupleType=typename IsTestOrTrial<Left>::UniqueElementFunctionSpacesTupleType;
  using type=typename IsTestOrTrial<Left>::type;
  static_assert(TupleTypeSize<type>::value<2," In Divide<Left,Real>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
};


template<typename Left, typename Right>
class IsTestOrTrial< Addition2<Expression2<Left>,Expression2<Right> > >
{
public:
  using TupleFunctionSpace=TupleCatType<typename IsTestOrTrial<Left>::TupleFunctionSpace,
                                   typename IsTestOrTrial<Right>::TupleFunctionSpace >;

  using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates<TupleCatType<typename IsTestOrTrial<Left>::UniqueElementFunctionSpacesTupleType,
                                                                                typename IsTestOrTrial<Right>::UniqueElementFunctionSpacesTupleType >>;
  using type=RemoveTupleDuplicates<TupleRemoveType<Number<0>,TupleCatType<typename IsTestOrTrial<Left>::type,typename IsTestOrTrial<Right>::type >>> ;
  static_assert(TupleTypeSize<type>::value<2," In Addition<Left,Right>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
};

template<typename Left, typename Right>
class IsTestOrTrial< Subtraction2<Expression2<Left>,Expression2<Right> > >
{
public:
  using TupleFunctionSpace=TupleCatType<typename IsTestOrTrial<Left>::TupleFunctionSpace,
                                   typename IsTestOrTrial<Right>::TupleFunctionSpace >;

  using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates<TupleCatType<typename IsTestOrTrial<Left>::UniqueElementFunctionSpacesTupleType,
                                                                                typename IsTestOrTrial<Right>::UniqueElementFunctionSpacesTupleType >>;
  using type=RemoveTupleDuplicates<TupleRemoveType<Number<0>,TupleCatType<typename IsTestOrTrial<Left>::type,typename IsTestOrTrial<Right>::type > >>;
  static_assert(TupleTypeSize<type>::value<2," In Subtraction<Left,Right>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");

};



template<typename Type>
class IsTestOrTrial< UnaryPlus2<Expression2<Type>> >
{
public:
  using TupleFunctionSpace=typename IsTestOrTrial<Type>::TupleFunctionSpace;
  using UniqueElementFunctionSpacesTupleType=typename IsTestOrTrial<Type>::UniqueElementFunctionSpacesTupleType;
  using type=TupleRemoveType<Number<0>,typename IsTestOrTrial<Type>::type>;
  static_assert(TupleTypeSize<type>::value<2," In UnaryPlus<Type>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");

};

template<typename Type>
class IsTestOrTrial< UnaryMinus2<Expression2<Type>> >
{
public:
  using TupleFunctionSpace=typename IsTestOrTrial<Type>::TupleFunctionSpace;
  using UniqueElementFunctionSpacesTupleType=typename IsTestOrTrial<Type>::UniqueElementFunctionSpacesTupleType;
  using type=TupleRemoveType<Number<0>,typename IsTestOrTrial<Type>::type>;
  static_assert(TupleTypeSize<type>::value<2," In UnaryMinus<Type>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");

};

template<typename Left,typename Right>
class TypeOfForm;


template<Integer M,Integer N>
class TypeOfForm<Number<M>,Number<N>>
{
public:
  using type=void;
  static_assert(0*Number<N>::value,"L2inner: the form is neither a 0-form(function,function), 1-form(function/test,test/function) or 2-form (test/trial,trial/test), where the function is neither a test nor a trial");
};


template<>
class TypeOfForm<Number<0>,Number<0>>
{
  public:
    using type=Number<0>; 
};

template<>
class TypeOfForm<Number<0>,Number<1>>
{
  public:
    using type=Number<1>; 
};

template<>
class TypeOfForm<Number<1>,Number<0>>
{
  public:
    using type=Number<1>; 
};

template<>
class TypeOfForm<Number<1>,Number<2>>
{
  public:
    using type=Number<2>; 
};


template<>
class TypeOfForm<Number<2>,Number<1>>
{
  public:
    using type=Number<2>; 
};


// template<typename MixedSpace, Integer N,typename OperatorType>
// class OperatorTupleType<Test<MixedSpace,N,OperatorType> >
// { public:
//   using Test=Test<MixedSpace,N,OperatorType>;
//   static constexpr Integer Nmax= Test::Nmax;
//   using single_type=std::tuple<std::tuple< OperatorType,std::tuple<> >>;
//   using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
//   using type=TupleChangeType<N,single_type,emptytuple>;
// };

// template<typename MixedSpace, Integer N,typename OperatorType>
// class OperatorTupleType<Trial<MixedSpace,N,OperatorType> >
// { public:
//   using Trial=Trial<MixedSpace,N,OperatorType>;
//   static constexpr Integer Nmax= Trial::Nmax;
//   using single_type=std::tuple<std::tuple< OperatorType,std::tuple<> >>;
//   using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
//   using type=TupleChangeType<N,single_type,emptytuple>;
// };


template<template<class,Integer,class > class TestOrTrial_,typename MixedSpace, Integer N,typename OperatorType>
class OperatorTupleType<TestOrTrial_<MixedSpace,N,OperatorType> >
{ public:
  using TestOrTrial=TestOrTrial_<MixedSpace,N,OperatorType>;
  static_assert((IsSame<TestOrTrial,Test<MixedSpace,N,OperatorType>>::value ||
                IsSame<TestOrTrial,Trial<MixedSpace,N,OperatorType>>::value )
               && "In Evaluation<Expression2<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");

  static constexpr Integer Nmax= TestOrTrial::Nmax;
  using single_type=std::tuple<std::tuple< OperatorType,std::tuple<> >>;
  using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
  using type=TupleChangeType<N,single_type,emptytuple>;
};




template<Integer N,typename...Args >
Test<        MixedSpace<Args...>,//MixedSpace<Args...>::UniqueElementFunctionSpacesTupleType,
             GetType<typename MixedSpace<Args...>::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value>
MakeTest(const MixedSpace<Args...>& W)
{return Test<MixedSpace<Args...>,//typename MixedSpace<Args...>::UniqueElementFunctionSpacesTupleType
        GetType<typename MixedSpace<Args...>::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value>();}


template<Integer N,typename...Args >
Trial<        MixedSpace<Args...>,//typename MixedSpace<Args...>::UniqueElementFunctionSpacesTupleType,
              GetType<typename MixedSpace<Args...>::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value> 
MakeTrial(const MixedSpace<Args...>& W)
{return Trial<MixedSpace<Args...>,//typename MixedSpace<Args...>::UniqueElementFunctionSpacesTupleType
        GetType<typename MixedSpace<Args...>::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value>();}



template<template<class,Integer,class>class TestOrTrial, typename MixedSpace,Integer N>
TestOrTrial<MixedSpace,N,GradientOperator> 
Grad(const TestOrTrial<MixedSpace,N,IdentityOperator>& W)
{
  static_assert((IsSame<TestOrTrial<MixedSpace,N,IdentityOperator>,Test<MixedSpace,N,IdentityOperator>>::value ||
                 IsSame<TestOrTrial<MixedSpace,N,IdentityOperator>,Trial<MixedSpace,N,IdentityOperator>>::value )
                 && "In Evaluation<Expression2<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");  
  return TestOrTrial<MixedSpace,N,GradientOperator> ();}


template<typename MixedSpace,Integer N>
Trial<MixedSpace,N,DivergenceOperator> 
Div(const Trial<MixedSpace,N,IdentityOperator>& W)
{return Trial<MixedSpace,N,DivergenceOperator> ();}

template<typename MixedSpace,Integer N>
Test<MixedSpace,N,DivergenceOperator> 
Div(const Test<MixedSpace,N,IdentityOperator>& W)
{return Test<MixedSpace,N,DivergenceOperator> ();}


template<typename MixedSpace,Integer N>
Trial<MixedSpace,N,CurlOperator> 
Curl(const Trial<MixedSpace,N,IdentityOperator>& W)
{return Trial<MixedSpace,N,CurlOperator> ();}

template<typename MixedSpace,Integer N>
Test<MixedSpace,N,CurlOperator> 
Curl(const Test<MixedSpace,N,IdentityOperator>& W)
{return Test<MixedSpace,N,CurlOperator> ();}



template<typename MixedSpace,Integer N, typename OperatorKind>
class QuadratureOrder<Test<MixedSpace,N,OperatorKind> >
{ public:
  using Test=Test<MixedSpace,N,OperatorKind>;
  using UniqueElementFunctionSpacesTupleType=typename Test::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename Test::Operator;
  using BaseFunctionSpace=GetType<BaseFunctionSpaceAndElement,1>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;
};

template<typename MixedSpace,Integer N, typename OperatorKind>
class QuadratureOrder<Trial<MixedSpace,N,OperatorKind> >
{ public:
  using Trial=Trial<MixedSpace,N,OperatorKind>;
  using UniqueElementFunctionSpacesTupleType=typename Trial::UniqueElementFunctionSpacesTupleType;
  using BaseFunctionSpaceAndElement=GetType<UniqueElementFunctionSpacesTupleType,N>;
  using Operator=typename Trial::Operator;
  using BaseFunctionSpace=GetType<BaseFunctionSpaceAndElement,1>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;
};





template<typename MeshT, typename Left,typename Right,Integer QR=GaussianQuadrature>
class L2DotProductIntegral: public Expression2<L2DotProductIntegral<MeshT,Left,Right,QR>>
{  
   public:
    using type=Contraction2<Expression2 <Left>, Expression2 <Right> > ;
    using Elem=typename MeshT::Elem;
    static constexpr Integer Order=CheckMaxQuadratureOrder<Elem,QR,QuadratureOrder<type>::value+1>::value; 
    using QRule=typename QuadratureRule<QR>:: template rule<Elem,Order>;
    using form= std::tuple<typename TypeOfForm<GetType<typename IsTestOrTrial<Left>::type,0>,
                                    GetType<typename IsTestOrTrial<Right>::type,0>>::type >;
    using UniqueElementFunctionSpacesTupleType=GetType<RemoveTupleDuplicates< TupleCatType< typename IsTestOrTrial<Left>::UniqueElementFunctionSpacesTupleType,
                                                                                            typename IsTestOrTrial<Right>::UniqueElementFunctionSpacesTupleType  >>,0>;               
    using TupleFunctionSpace=RemoveTupleDuplicates< TupleCatType< typename IsTestOrTrial<Left>::TupleFunctionSpace,
                                                             typename IsTestOrTrial<Right>::TupleFunctionSpace  >>;               
    using FunctionSpace=GetType<TupleFunctionSpace,0>;
    L2DotProductIntegral(const MeshT& mesh,const Expression2<Left>& left,const Expression2<Right>& right,const Integer label=-666):
    mesh_(mesh),
    left_(left),
    right_(right),
    product_(Inner(left,right)),
    label_(label)
    {}
     
     const Left&  left() const{return left_;};
     const Right& right()const{return right_;};
  private:
    MeshT mesh_;
    Left left_;
    Right right_;
    type product_;
    Integer label_;
};



template<typename MeshT, typename Left,typename Right,Integer QR>
class OperatorTupleType<L2DotProductIntegral<MeshT,Left,Right,QR>>
{ 
public:
  using L2prod=L2DotProductIntegral<MeshT,Left,Right,QR>;
  using QRule=typename L2prod::QRule;
  using Operatortuple=typename OperatorTupleType<typename L2prod::type>::type;
  using type=TupleOfTupleChangeType<1,QRule, Operatortuple>;;
} 
;                               
 




template<typename MeshT, typename Left,typename Right>
L2DotProductIntegral<MeshT,Left,Right>
L2Inner(const MeshT& mesh,const Expression2<Left>& left,const Expression2<Right>& right)
{return L2DotProductIntegral<MeshT,Left,Right>(mesh,left,right);}




template<typename MeshT, typename Left,typename Right>
L2DotProductIntegral<MeshT,Left,Right>
L2Inner(const MeshT& mesh,const Expression2<Left>& left,const Expression2<Right>& right, const Integer label)
{return L2DotProductIntegral<MeshT,Left,Right>(mesh,left,right,label);}











template<typename DerivedLeft, typename MeshT, typename Left,typename Right>
class Addition2< Expression2<DerivedLeft>  ,  Expression2<L2DotProductIntegral<MeshT,Left,Right>>>
: public Expression2< Addition2< Expression2<DerivedLeft>  ,  Expression2<L2DotProductIntegral<MeshT,Left,Right>>  > > 
{

  public:
  using TupleFunctionSpace=RemoveTupleDuplicates< TupleCatType< typename L2DotProductIntegral<MeshT,Left,Right>::TupleFunctionSpace,
                                                           typename DerivedLeft::TupleFunctionSpace  >>;               

  using FunctionSpace=GetType<TupleFunctionSpace,0>;

  using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::UniqueElementFunctionSpacesTupleType,
                                                                                 typename DerivedLeft::UniqueElementFunctionSpacesTupleType > >;
  using form =RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::form,
                                                  typename DerivedLeft::form  > >;
    Addition2(const Expression2<DerivedLeft>& left, const Expression2<L2DotProductIntegral<MeshT,Left,Right>>&right)
    : 
    left_(left.derived()),
    right_(right.derived())
    {};
    const DerivedLeft & left()const{return left_;};
    const L2DotProductIntegral<MeshT,Left,Right>& right()const{return right_;};
  private:
  DerivedLeft left_;
  L2DotProductIntegral<MeshT,Left,Right> right_;
};


template<typename DerivedLeft, typename MeshT, typename Left,typename Right>
class Subtraction2< Expression2<DerivedLeft>  ,  Expression2<L2DotProductIntegral<MeshT,Left,Right>>>
: public Expression2< Subtraction2< Expression2<DerivedLeft>  ,  Expression2<L2DotProductIntegral<MeshT,Left,Right>>  > > 
{

  public:
  using TupleFunctionSpace=RemoveTupleDuplicates< TupleCatType< typename L2DotProductIntegral<MeshT,Left,Right>::TupleFunctionSpace,
                                                           typename DerivedLeft::TupleFunctionSpace  >>;               

  using FunctionSpace=GetType<TupleFunctionSpace,0>;

  using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::UniqueElementFunctionSpacesTupleType,
                                                                                 typename DerivedLeft::UniqueElementFunctionSpacesTupleType > >;
  using form =RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::form,
                                                  typename DerivedLeft::form  > >;
    Subtraction2(const Expression2<DerivedLeft>& left, const Expression2<L2DotProductIntegral<MeshT,Left,Right>>&right)
    : 
    left_(left.derived()),
    right_(right.derived())
    {};
    const DerivedLeft & left()const{return left_;};
    const L2DotProductIntegral<MeshT,Left,Right>& right()const{return right_;};
  private:
  DerivedLeft left_;
  L2DotProductIntegral<MeshT,Left,Right> right_;
};





template<typename Form>
class ShapeFunctions;

template<typename MeshT, typename Left,typename Right,Integer QR, typename Form>
class Evaluation<Expression2<L2DotProductIntegral<MeshT,Left,Right,QR>>, ShapeFunctions<Form>>
{
 public:
 using type= L2DotProductIntegral<MeshT,Left,Right,QR>;
 using EvalLeft=Evaluation<Expression2<Left>>;
 using EvalRight=Evaluation<Expression2<Right>>;
 Evaluation(){};
 
 Evaluation(const type& expr, ShapeFunctions<Form>& shape_functions):
 eval_(expr),
 eval_left_(EvalLeft(eval_.left())),
 eval_right_(EvalRight(eval_.right())),
 shape_functions_(shape_functions)
 {};
 
 template<typename QRule,typename ...Ts>
 void apply_aux(const std::tuple<Ts...>& tuple)
 {
  eval_left_.template apply<QRule>(tuple);
  eval_right_.template apply<QRule>(tuple);
 }

 template<typename ...Ts>
 void apply(const std::tuple<Ts...>& tuple)
 {
  apply_aux<typename type::QRule>(tuple);
 }

 template<Integer Rows,Integer Cols>
 void apply(const Matrix<Real,Rows,Cols>& mat)
 {
  // eval_left_.apply(mat);
  // eval_right_.apply(mat); 
  // apply_aux<typename type::QRule>(tuple);
 }
 
private:
 type eval_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
 ShapeFunctions<Form>& shape_functions_;
};


template<typename MeshT, typename Left,typename Right,typename DerivedLeft, typename Form>
class Evaluation<Expression2<Addition2< Expression2<DerivedLeft>  ,  
                                        Expression2<L2DotProductIntegral<MeshT,Left,Right>>>>,
                 ShapeFunctions<Form>>
{
 public:
 using type=Addition2<Expression2<DerivedLeft>,Expression2<L2DotProductIntegral<MeshT,Left,Right>>>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>,ShapeFunctions<Form>>;
 using EvalRight=Evaluation<Expression2<L2DotProductIntegral<MeshT,Left,Right>>,ShapeFunctions<Form>>;
 Evaluation(){};
 

 Evaluation(const Expression2<type>& expr, ShapeFunctions<Form>& shape_functions):
 expr_(expr.derived()),
 eval_left_(EvalLeft(expr_.left(),shape_functions)),
 eval_right_(EvalRight(expr_.right(),shape_functions)),
 shape_functions_(shape_functions)
{};
 
 template<Integer Rows,Integer Cols>
 void apply(const Matrix<Real,Rows,Cols>& mat)
 {
  eval_left_.apply(mat);
  eval_right_.apply(mat); 
 }
private:

 type expr_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
 ShapeFunctions<Form>& shape_functions_;
};

template<typename MeshT, typename Left,typename Right,typename DerivedLeft, typename Form>
class Evaluation<Expression2<Subtraction2< Expression2<DerivedLeft>  ,  
                                        Expression2<L2DotProductIntegral<MeshT,Left,Right>>>>,
                 ShapeFunctions<Form>>{
 public:
 using type=Subtraction2<Expression2<DerivedLeft>,Expression2<L2DotProductIntegral<MeshT,Left,Right>>>;
 using EvalLeft=Evaluation<Expression2<DerivedLeft>,ShapeFunctions<Form>>;
 using EvalRight=Evaluation<Expression2<L2DotProductIntegral<MeshT,Left,Right>>,ShapeFunctions<Form>>;
 Evaluation(){};
 

 Evaluation(const type& expr, ShapeFunctions<Form>& shape_functions):
 expr_(expr),
 eval_left_(EvalLeft(expr_.left(),shape_functions)),
 eval_right_(EvalRight(expr_.right(),shape_functions)),
 shape_functions_(shape_functions)
 {};
 
 template<Integer Rows,Integer Cols>
 void apply(const Matrix<Real,Rows,Cols>& mat)
 {
  eval_left_.apply(mat);
  eval_right_.apply(mat); 
 }
private:

 type expr_;
 EvalLeft eval_left_;
 EvalRight eval_right_;
 ShapeFunctions<Form>& shape_functions_;
};
























template<typename Form>
class ShapeFunctions
{
 public:
  using UniqueElementFunctionSpacesTupleType=typename Form::UniqueElementFunctionSpacesTupleType;
  using TupleOperatorsAndQuadrature= typename OperatorTupleType<Form>::type;
  using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadrature>;
  using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies<UniqueElementFunctionSpacesTupleType>;
  using Map=MapOperatorTupleOfTuple<TupleOfTupleNoQuadrature,UniqueElementFunctionSpacesTupleType>;
  using UniqueMapping=UniqueMap<SpacesToUniqueFEFamilies,Map> ;
  using TupleOfTupleShapeFunction=TupleOfTupleShapeFunctionType<UniqueElementFunctionSpacesTupleType,TupleOperatorsAndQuadrature>;
  using MapTupleNumbers=MapTupleInit<TupleOfTupleShapeFunction, SpacesToUniqueFEFamilies, UniqueMapping>;
 
  
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


   template<typename ConstFormReference>
   constexpr void init_map(const ReferenceMaps<ConstFormReference>& maps) 
   {init_map_aux<TupleTypeSize<MapTupleNumbers>::value-1,0>(maps());}


//////////////////////////////////////////////////



  template<Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
  constexpr typename std::enable_if_t< FEFamily==LagrangeFE,void> 
  shape_function_init_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficient<Args...> &coefficients)
  {
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



  template<typename Coefficients>
  constexpr void init(const Coefficients& shape_coefficients)
  {
   constexpr Integer Nmax=TupleTypeSize<TupleOfTupleShapeFunction>::value-1;
   shape_function_init_aux<Nmax,0>(shape_coefficients);
   }


  template<typename ConstFormReference, typename Coefficients>
  constexpr void init(const ReferenceMaps<ConstFormReference>& maps,const Coefficients& shape_coefficients)
  {
    init_map(maps);
    init(shape_coefficients);
   }




   constexpr       auto& operator()()      {return tuple;}
   constexpr const auto& operator()()const {return tuple;}


   template<Integer...Ns>
   constexpr const auto& get()const{return tuple_get<Ns...>(tuple);}

   template<Integer...Ns>
         auto& get()     {return tuple_get<Ns...>(tuple);}

   template<Integer...Ns>
   constexpr const auto& value()const{return tuple_get<Ns...>(tuple).eval();}

   template<Integer...Ns>
         auto& value()     {return tuple_get<Ns...>(tuple).eval();}




private:
 TupleOfTupleShapeFunction tuple;
};


template<typename ConstFormReference>
constexpr auto shape_functions(const ConstFormReference& form)
{using Form=typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type;
 return ShapeFunctions<Form>();  }


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
// struct IsTestOrTrial : std::false_type { };

// template <typename T>
// struct IsTestOrTrial <T,  decltype((void) T::is_test_or_trial, static_cast<decltype(T::is_test_or_trial)>(T::is_test_or_trial) )> : std::true_type { };




}


#endif