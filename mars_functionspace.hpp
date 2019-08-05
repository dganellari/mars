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
  using ElemsTupleType=RemoveTupleDuplicates<TupleCatType<typename Args::ElemsTupleType...>>;
  using FromElementFunctionSpacesToUniqueNumbersTupleType=TupleAllToUniqueMap<ElementFunctionSpacesTupleType,UniqueElementFunctionSpacesTupleType>;
  inline const Integer& n_dofs()const{return n_dofs_;}; 

  inline const DofMapType& dofmap()const{return dofmap_;};


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


template<typename BaseFunctionSpacesType,Integer N, typename OperatorType=IdentityOperator>
class Trial: public Expression2<Trial<BaseFunctionSpacesType,N,OperatorType>>
{ public:
  static constexpr Integer Nmax=TupleTypeSize<BaseFunctionSpacesType>::value;
  static constexpr Integer value=N;
  using BaseFunctionSpaces=BaseFunctionSpacesType;
  using Operator=OperatorType;
};

template<typename BaseFunctionSpacesType, Integer N, typename OperatorType=IdentityOperator>
class Test: public Expression2<Test<BaseFunctionSpacesType,N,OperatorType>>
{
public:
  static constexpr Integer Nmax=TupleTypeSize<BaseFunctionSpacesType>::value;
  static constexpr Integer value=N;
  using BaseFunctionSpaces=BaseFunctionSpacesType;
  using Operator=OperatorType;
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
  using UniqueBaseFunctionSpaces=std::tuple<>;
  using type=std::tuple<Number<0>>;
};

template<typename BaseFunctionSpacesType,Integer N, typename OperatorType>
class IsTestOrTrial<Test<BaseFunctionSpacesType,N,OperatorType>>{
public:
  using UniqueBaseFunctionSpaces=std::tuple<BaseFunctionSpacesType>;
  using type=std::tuple<Number<1>>;
};

template<typename BaseFunctionSpacesType,Integer N, typename OperatorType>
class IsTestOrTrial<Trial<BaseFunctionSpacesType,N,OperatorType>>{
public:
  using UniqueBaseFunctionSpaces=std::tuple<BaseFunctionSpacesType>;
  using type=std::tuple<Number<2>>;
};


template <typename T,typename ... Types>
class TupleTypeSize;


template<typename Left, typename Right>
class IsTestOrTrial< Multiply2<Expression2<Left>,Expression2<Right> > >
{
public:
  using UniqueBaseFunctionSpaces=TupleCatType<typename IsTestOrTrial<Left>::UniqueBaseFunctionSpaces,typename IsTestOrTrial<Right>::UniqueBaseFunctionSpaces >;
  using type=TupleRemoveType<Number<0>,TupleCatType<typename IsTestOrTrial<Left>::type,typename IsTestOrTrial<Right>::type >> ;
  static_assert(TupleTypeSize<type>::value<2," In Multiply<Left,Right>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
};

template<typename Left>
class IsTestOrTrial< Multiply2<Expression2<Left>,Real > >
{
public:
  using UniqueBaseFunctionSpaces=typename IsTestOrTrial<Left>::UniqueBaseFunctionSpaces;
  using type=typename IsTestOrTrial<Left>::type;
  static_assert(TupleTypeSize<type>::value<2," In Multiply<Left,Real>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
};

template<typename Right>
class IsTestOrTrial< Multiply2<Real,Expression2<Right> > >
{
public:
  using UniqueBaseFunctionSpaces=typename IsTestOrTrial<Right>::UniqueBaseFunctionSpaces;
  using type=typename IsTestOrTrial<Right>::type;
  static_assert(TupleTypeSize<type>::value<2," In Multiply<Real,Right>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");

};


template<typename Left>
class IsTestOrTrial< Divide2<Expression2<Left>,Real > >
{
public:
  using UniqueBaseFunctionSpaces=typename IsTestOrTrial<Left>::UniqueBaseFunctionSpaces;
  using type=typename IsTestOrTrial<Left>::type;
  static_assert(TupleTypeSize<type>::value<2," In Divide<Left,Real>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
};


template<typename Left, typename Right>
class IsTestOrTrial< Addition2<Expression2<Left>,Expression2<Right> > >
{
public:
  using UniqueBaseFunctionSpaces=RemoveTupleDuplicates<TupleCatType<typename IsTestOrTrial<Left>::UniqueBaseFunctionSpaces,typename IsTestOrTrial<Right>::UniqueBaseFunctionSpaces >>;
  using type=RemoveTupleDuplicates<TupleRemoveType<Number<0>,TupleCatType<typename IsTestOrTrial<Left>::type,typename IsTestOrTrial<Right>::type >>> ;
  static_assert(TupleTypeSize<type>::value<2," In Addition<Left,Right>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
};

template<typename Left, typename Right>
class IsTestOrTrial< Subtraction2<Expression2<Left>,Expression2<Right> > >
{
public:
  using UniqueBaseFunctionSpaces=RemoveTupleDuplicates<TupleCatType<typename IsTestOrTrial<Left>::UniqueBaseFunctionSpaces,typename IsTestOrTrial<Right>::UniqueBaseFunctionSpaces >>;
  using type=RemoveTupleDuplicates<TupleRemoveType<Number<0>,TupleCatType<typename IsTestOrTrial<Left>::type,typename IsTestOrTrial<Right>::type > >>;
  static_assert(TupleTypeSize<type>::value<2," In Subtraction<Left,Right>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");

};



template<typename Type>
class IsTestOrTrial< UnaryPlus2<Expression2<Type>> >
{
public:
  using UniqueBaseFunctionSpaces=typename IsTestOrTrial<Type>::UniqueBaseFunctionSpaces;
  using type=TupleRemoveType<Number<0>,typename IsTestOrTrial<Type>::type>;
  static_assert(TupleTypeSize<type>::value<2," In UnaryPlus<Type>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");

};

template<typename Type>
class IsTestOrTrial< UnaryMinus2<Expression2<Type>> >
{
public:
  using UniqueBaseFunctionSpaces=typename IsTestOrTrial<Type>::UniqueBaseFunctionSpaces;
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



template<typename BaseFunctionSpacesType, Integer N,typename OperatorType>
class OperatorTupleType<Test<BaseFunctionSpacesType,N,OperatorType> >
{ public:
  using Test=Test<BaseFunctionSpacesType,N,OperatorType>;
  static constexpr Integer Nmax= Test::Nmax;
  using single_type=std::tuple<std::tuple< OperatorType,std::tuple<> >>;
  using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
  using type=TupleChangeType<N,single_type,emptytuple>;
};

template<typename BaseFunctionSpacesType, Integer N,typename OperatorType>
class OperatorTupleType<Trial<BaseFunctionSpacesType,N,OperatorType> >
{ public:
  using Trial=Trial<BaseFunctionSpacesType,N,OperatorType>;
  static constexpr Integer Nmax= Trial::Nmax;
  using single_type=std::tuple<std::tuple< OperatorType,std::tuple<> >>;
  using emptytuple=TupleOfType<Nmax,std::tuple<> > ;
  using type=TupleChangeType<N,single_type,emptytuple>;
};






template<Integer N,typename...Args >
Test<typename MixedSpace<Args...>::UniqueElementFunctionSpacesTupleType,
             GetType<N,typename MixedSpace<Args...>::FromElementFunctionSpacesToUniqueNumbersTupleType>::value>
MakeTest(const MixedSpace<Args...>& W)
{return Test<typename MixedSpace<Args...>::UniqueElementFunctionSpacesTupleType,
        GetType<N,typename MixedSpace<Args...>::FromElementFunctionSpacesToUniqueNumbersTupleType>::value>();}


template<Integer N,typename...Args >
Trial<          typename MixedSpace<Args...>::UniqueElementFunctionSpacesTupleType,
              GetType<N,typename MixedSpace<Args...>::FromElementFunctionSpacesToUniqueNumbersTupleType>::value> 
MakeTrial(const MixedSpace<Args...>& W)
{return Trial<typename MixedSpace<Args...>::UniqueElementFunctionSpacesTupleType,
        GetType<N,typename MixedSpace<Args...>::FromElementFunctionSpacesToUniqueNumbersTupleType>::value>();}


template<typename BaseFunctionSpacesType,Integer N>
Trial<BaseFunctionSpacesType,N,GradientOperator> 
Grad(const Trial<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Trial<BaseFunctionSpacesType,N,GradientOperator> ();}

template<typename BaseFunctionSpacesType,Integer N>
Test<BaseFunctionSpacesType,N,GradientOperator> 
Grad(const Test<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Test<BaseFunctionSpacesType,N,GradientOperator> ();}


template<typename BaseFunctionSpacesType,Integer N>
Trial<BaseFunctionSpacesType,N,DivergenceOperator> 
Div(const Trial<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Trial<BaseFunctionSpacesType,N,DivergenceOperator> ();}

template<typename BaseFunctionSpacesType,Integer N>
Test<BaseFunctionSpacesType,N,DivergenceOperator> 
Div(const Test<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Test<BaseFunctionSpacesType,N,DivergenceOperator> ();}


template<typename BaseFunctionSpacesType,Integer N>
Trial<BaseFunctionSpacesType,N,CurlOperator> 
Curl(const Trial<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Trial<BaseFunctionSpacesType,N,CurlOperator> ();}

template<typename BaseFunctionSpacesType,Integer N>
Test<BaseFunctionSpacesType,N,CurlOperator> 
Curl(const Test<BaseFunctionSpacesType,N,IdentityOperator>& W)
{return Test<BaseFunctionSpacesType,N,CurlOperator> ();}



template<typename BaseFunctionSpacesType,Integer N, typename OperatorKind>
class QuadratureOrder<Test<BaseFunctionSpacesType,N,OperatorKind> >
{ public:

  using BaseFunctionSpaceAndElement=GetType<N,BaseFunctionSpacesType>;
  using Operator=typename Test<BaseFunctionSpacesType,N,OperatorKind>::Operator;
  using BaseFunctionSpace=GetType<1,BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;
};

template<typename BaseFunctionSpacesType,Integer N, typename OperatorKind>
class QuadratureOrder<Trial<BaseFunctionSpacesType,N,OperatorKind> >
{ public:
  using BaseFunctionSpaceAndElement=GetType<N,BaseFunctionSpacesType>;
  using Operator=typename Trial<BaseFunctionSpacesType,N,OperatorKind>::Operator;
  using BaseFunctionSpace=GetType<1,BaseFunctionSpaceAndElement>;
  static constexpr Integer value=QuadratureOrder<Operator,BaseFunctionSpace>::value;
};





template<typename MeshT, typename Left,typename Right,Integer QR=GaussianQuadrature>
class L2DotProductIntegral: public Expression2<L2DotProductIntegral<MeshT,Left,Right>>
{  
   public:
    using type=Contraction2<Expression2 <Left>, Expression2 <Right> > ;
    using Elem=typename MeshT::Elem;
    static constexpr Integer Order=CheckMaxQuadratureOrder<Elem,QR,QuadratureOrder<type>::value+1>::value; 
    using QRule=typename QuadratureRule<QR>:: template rule<Elem,Order>;
    using form= std::tuple<typename TypeOfForm<GetType<0,typename IsTestOrTrial<Left>::type>,
                                    GetType<0,typename IsTestOrTrial<Right>::type>>::type >;
    using UniqueBaseFunctionSpaces=GetType<0,RemoveTupleDuplicates< TupleCatType< typename IsTestOrTrial<Left>::UniqueBaseFunctionSpaces,typename IsTestOrTrial<Right>::UniqueBaseFunctionSpaces >>>;               
    L2DotProductIntegral(const MeshT& mesh,const Expression2<Left>& left,const Expression2<Right>& right ):
    mesh_(mesh),
    left_(left),
    right_(right),
    product_(Inner(left,right))
    {}

  private:
    MeshT mesh_;
    Left left_;
    Right right_;
    type product_;
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















template<typename DerivedLeft, typename MeshT, typename Left,typename Right>
class Addition2< Expression2<DerivedLeft>  ,  Expression2<L2DotProductIntegral<MeshT,Left,Right>>>
: public Expression2< Addition2< Expression2<DerivedLeft>  ,  Expression2<L2DotProductIntegral<MeshT,Left,Right>>  > > 
{

  public:
  using UniqueBaseFunctionSpaces=RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::UniqueBaseFunctionSpaces,
                                                  typename DerivedLeft::UniqueBaseFunctionSpaces > >;
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
  using UniqueBaseFunctionSpaces=RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::UniqueBaseFunctionSpaces,
                                                  typename DerivedLeft::UniqueBaseFunctionSpaces > >;
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


// template<typename DerivedLeft, typename MeshT, typename Left,typename Right>
// class Multiply2< Expression2<DerivedLeft>  ,  Expression2<L2DotProductIntegral<MeshT,Left,Right>>>
// : public Expression2< Multiply2< Expression2<DerivedLeft>  ,  Expression2<L2DotProductIntegral<MeshT,Left,Right>>  > > 
// {

//   public:
//   using L2prod=L2DotProductIntegral<MeshT,Left,Right>;
//   static_assert(0*GetType<0,typename L2prod::form>::value, "You cannot multiply a L2 inner product: you can only add or subtract two or more inner products");
//   using UniqueBaseFunctionSpaces=RemoveTupleDuplicates< TupleCatType<typename L2prod::UniqueBaseFunctionSpaces,
//                                                   typename DerivedLeft::UniqueBaseFunctionSpaces > >;
//   using form =RemoveTupleDuplicates< TupleCatType<typename L2prod::form,
//                                                   typename DerivedLeft::form  > >;
//     Multiply2(const Expression2<DerivedLeft>& left, const Expression2<L2prod>&right)
//     : 
//     left_(left.derived()),
//     right_(right.derived())
//     {};
//     const DerivedLeft & left()const{return left_;};
//     const L2prod& right()const{return right_;};
//   private:
//   DerivedLeft left_;
//   L2prod right_;
// };

// template<typename DerivedLeft, typename MeshT, typename Left,typename Right>
// class Multiply2<  Expression2<L2DotProductIntegral<MeshT,Left,Right>>, Expression2<DerivedLeft>>
// : public Expression2< Multiply2<  Expression2<L2DotProductIntegral<MeshT,Left,Right>> ,Expression2<DerivedLeft> > > 
// {

//   public:
//   static_assert(0*L2DotProductIntegral<MeshT,Left,Right>::form::value, "You cannot multiply a L2 inner product: you can only add or subtract two or more inner products");
//   using UniqueBaseFunctionSpaces=RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::UniqueBaseFunctionSpaces,
//                                                   typename DerivedLeft::UniqueBaseFunctionSpaces > >;
//   using form =RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::form,
//                                                   typename DerivedLeft::form  > >;
//     Multiply2(const Expression2<L2DotProductIntegral<MeshT,Left,Right>>& left, const Expression2<DerivedLeft>&right)
//     : 
//     left_(left.derived()),
//     right_(right.derived())
//     {};
//     const L2DotProductIntegral<MeshT,Left,Right>& left()const{return left_;};
//     const DerivedLeft & right()const{return right_;};
//   private:
//   L2DotProductIntegral<MeshT,Left,Right> left_;
//   DerivedLeft  right_;
// };


template<typename ConstFormReference>
void Assembly(const ConstFormReference& form)
{
  using Form=typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type;
  using UniqueBaseFunctionSpaces=typename Form::UniqueBaseFunctionSpaces;
  using TupleOperatorsAndQuadrature= typename OperatorTupleType<Form>::type;
  using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadrature>;
  using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies<UniqueBaseFunctionSpaces>;
  using Map=MapOperatorTupleOfTuple<TupleOfTupleNoQuadrature,UniqueBaseFunctionSpaces>;
  using UniqueMapping=UniqueMap<SpacesToUniqueFEFamilies,Map> ;
  using TupleOfTupleShapeFunction=TupleOfTupleShapeFunctionType<UniqueBaseFunctionSpaces,TupleOperatorsAndQuadrature>;
  using MapTupleNumbers=MapTupleInit<TupleOfTupleShapeFunction, SpacesToUniqueFEFamilies, UniqueMapping>;
 
  TupleOfTupleShapeFunction stuple;
  UniqueMapping mapping;
  init_map_aux< SpacesToUniqueFEFamilies,MapTupleNumbers,TupleTypeSize<MapTupleNumbers>::value-1,0>(stuple,mapping);
}


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