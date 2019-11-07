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
#include "mars_evaluation.hpp"
#include "mars_vector.hpp"
#include "mars_number.hpp"



namespace mars{




template<typename MeshT, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
class FunctionSpace 
{
public:

      using Elem= typename MeshT::Elem;
      using TupleOfSpaces=std::tuple<ElemFunctionSpace<Elem,BaseFunctionSpace>,ElemFunctionSpace<Elem,BaseFunctionSpaces>... >;
      static constexpr Integer Nsubspaces=1+sizeof...(BaseFunctionSpaces);

      static constexpr Integer Nelem_dofs=DofsPerElemNums<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;
      static constexpr Array<Integer,Nsubspaces> Nelem_dofs_array=
      concat(Array<Integer,1>{DofsPerElemNums<Elem,BaseFunctionSpace>::value},Array<Integer,1>{DofsPerElemNums<Elem,BaseFunctionSpaces>::value}...);
      using DofMapType=std::vector<std::array<Integer, Nelem_dofs>>;
      using DofMapType2=std::tuple<std::vector<Array<Integer, DofsPerElemNums<Elem,BaseFunctionSpace>::value>>,
                                   std::vector<Array<Integer, DofsPerElemNums<Elem,BaseFunctionSpaces>::value>>
                                  ...>;
      static constexpr auto faces_dofs_array=std::tuple_cat(std::make_tuple(trace_dofs<ElemFunctionSpace<Elem,BaseFunctionSpace>>()),
                                                             std::make_tuple(trace_dofs<ElemFunctionSpace<Elem,BaseFunctionSpaces>>())...);
      static constexpr Array<std::size_t,Nsubspaces> Nfaces_dofs_array{trace_dofs<ElemFunctionSpace<Elem,BaseFunctionSpace>>()[0].size(),trace_dofs<ElemFunctionSpace<Elem,BaseFunctionSpaces>>()[0].size()...};


      using OffSetType=std::array<std::vector<Integer>, Nsubspaces>;
      using SpacesDofsArrayType=std::array<std::vector<std::vector<Integer>>, Nsubspaces>;
      using SpacesInfosArrayType=std::array<std::array<Integer,4>,Nsubspaces>;
      using ElemsTupleType=std::tuple<Elem>;
      // using ElementFunctionSpacesTupleType=std::tuple<std::tuple<Elem,BaseFunctionSpace>,std::tuple<Elem,BaseFunctionSpaces>...>;
      // using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates<ElementFunctionSpacesTupleType>;
      // using SpacesToUniqueFEFamily=SpacesToUniqueFEFamilies<UniqueElementFunctionSpacesTupleType>;
      // using FromElementFunctionSpacesToUniqueNumbersTupleType=TupleAllToUniqueMap<ElementFunctionSpacesTupleType,UniqueElementFunctionSpacesTupleType>;

      using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates<TupleOfSpaces>;
      using SpacesToUniqueFEFamily=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;
      using FromElementFunctionSpacesToUniqueNumbersTupleType=
            TupleAllToUniqueMap<TupleOfSpaces,UniqueElementFunctionSpacesTupleType>;
      using FromElementFunctionSpacesToFirstSpaceTupleType=
            FromElementFunctionSpacesToFirstSpaceTupleType<FromElementFunctionSpacesToUniqueNumbersTupleType>;

      static constexpr Integer Nuniquesubspaces=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;


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

      inline const DofMapType2& dofmap2()const{return dofmap2_;};

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

      inline auto mesh_ptr()const {return mesh_ptr_;};
       
      FunctionSpace(const MeshT& mesh,const Integer dof_count_start=0):
      mesh_ptr_(std::make_shared< MeshT >(mesh))
      {
      function_space_info<0,Nsubspaces,BaseFunctionSpace,BaseFunctionSpaces...>(space_infos_);
      // dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap_,offset_,n_dofs_,space_infos_,space_dofs_);     
      dofmap_fespace2<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap2_,offset_,n_dofs_,space_infos_,space_dofs_);     
      };

private:
   
      std::shared_ptr< MeshT > mesh_ptr_;
      Integer n_dofs_;
      DofMapType dofmap_;
      DofMapType2 dofmap2_;
      OffSetType offset_;
      SpacesDofsArrayType space_dofs_;
      SpacesInfosArrayType space_infos_;
      // ElementFunctionSpacesTupleType shape_functions_;

};


template<typename MeshT, typename BaseFunctionSpace>
using AuxFunctionSpace=FunctionSpace<MeshT,BaseFunctionSpace>;




template<typename...Args>
class MixedSpace; 

template<typename Arg,typename...Args>
class MixedSpace<Arg,Args...>
{
public:
  using Elem=typename Arg::Elem;
  using MeshT=Mesh<Elem::Dim,Elem::ManifoldDim>;
  using TupleOfSpaces=TupleCatType<typename Arg::TupleOfSpaces,typename Args::TupleOfSpaces...>;
  using DofMapType=std::tuple<typename Arg::DofMapType,typename Args::DofMapType...>;
  using DofMapType2=TupleCatType<typename Arg::DofMapType2,typename Args::DofMapType2...>;
  // using ElementFunctionSpacesTupleType=TupleCatType<typename Arg::ElementFunctionSpacesTupleType,typename Args::ElementFunctionSpacesTupleType...>;
  // using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates<TupleCatType<typename Arg::UniqueElementFunctionSpacesTupleType,typename Args::UniqueElementFunctionSpacesTupleType...>>;
  // using SpacesToUniqueFEFamily=SpacesToUniqueFEFamilies<UniqueElementFunctionSpacesTupleType>;
  // using ElemsTupleType=RemoveTupleDuplicates<TupleCatType<typename Arg::ElemsTupleType,typename Args::ElemsTupleType...>>;
  // using FromElementFunctionSpacesToUniqueNumbersTupleType=TupleAllToUniqueMap<ElementFunctionSpacesTupleType,UniqueElementFunctionSpacesTupleType>;



  using UniqueElementFunctionSpacesTupleType=
  RemoveTupleDuplicates<TupleCatType<typename Arg::UniqueElementFunctionSpacesTupleType,typename Args::UniqueElementFunctionSpacesTupleType...>>;
  using SpacesToUniqueFEFamily=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;
  using FromElementFunctionSpacesToUniqueNumbersTupleType=
        TupleAllToUniqueMap<TupleOfSpaces,UniqueElementFunctionSpacesTupleType>;
  using FromElementFunctionSpacesToFirstSpaceTupleType=
        FromElementFunctionSpacesToFirstSpaceTupleType<FromElementFunctionSpacesToUniqueNumbersTupleType>;

  static constexpr Integer Nuniquesubspaces=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
  // using UniqueElementFunctionSpacesTupleType2=
  // RemoveTupleDuplicates<TupleCatType<typename Arg::UniqueElementFunctionSpacesTupleType2,typename Args::UniqueElementFunctionSpacesTupleType2...>>;
  // using SpacesToUniqueFEFamily2=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType2>;
  // using FromElementFunctionSpacesToUniqueNumbersTupleType2=
  // TupleAllToUniqueMap<TupleOfSpaces,UniqueElementFunctionSpacesTupleType2>;





 //  template<Integer N, typename T,typename...Ts>
 //  class FunctionSpacesCount;

 //  template<Integer N, typename...Ts>
 //  class FunctionSpacesCount<N,std::tuple<FunctionSpace<Ts...>>>
 //  {
 //   public: 
 //   using T=FunctionSpace<Ts...>;
 //   using type=TupleOfTypeTCreate<Number<N>,T::Nsubspaces>;     
 //   static constexpr Integer value= N+1;
 //   };

 //  template<Integer N, typename...Ts>
 //  class FunctionSpacesCount<N,std::tuple<MixedSpace<Ts...>>>
 //  {
 //   public: 
 //   using type=typename FunctionSpacesCount<N,typename MixedSpace<Ts...>::Spaces>::type;
 //   static constexpr Integer value= FunctionSpacesCount<N,typename MixedSpace<Ts...>::Spaces>::value;
 //  };

 //  template<Integer N, typename ...Ts1, typename...Ts>
 //  class FunctionSpacesCount<N,std::tuple<FunctionSpace<Ts1...>,Ts...>>
 //  {
 //   public: 
 //   using T=FunctionSpace<Ts1...>;
 //   using singletype=TupleOfTypeTCreate<Number<N>,T::Nsubspaces>;
 //   using rest = typename FunctionSpacesCount<N+1,std::tuple<Ts...>>::type;
 //   using type=TupleCatType<singletype,rest>; 
 //   static constexpr Integer value=FunctionSpacesCount<N+1,std::tuple<Ts...>>::value;
 // };

 //  template<Integer N, typename ...Ts1, typename...Ts>
 //  class FunctionSpacesCount<N,std::tuple<MixedSpace<Ts1...>,Ts...>>
 //  {
 //   public: 
 //   using T=MixedSpace<Ts1...>;
 //   using singletype=typename FunctionSpacesCount<N,typename T::Spaces>::type;
 //   static constexpr Integer value_tmp= FunctionSpacesCount<N,typename T::Spaces>::value;
 //   using rest = typename FunctionSpacesCount<value_tmp,std::tuple<Ts...>>::type;
 //   static constexpr Integer value= FunctionSpacesCount<value_tmp,std::tuple<Ts...>>::value;
 //   using type=TupleCatType<singletype,rest>;
 //  };

  // using FromSpacesToFunctionSpaces=typename FunctionSpacesCount<0, Spaces>::type;

  inline const Integer& n_dofs()const{return n_dofs_;}; 

  inline const DofMapType& dofmap()const{return dofmap_;};
  
  template<Integer...Ns>
  inline constexpr const auto& dofmap()const{return tuple_get<Ns...>(dofmap_);};

      inline const DofMapType2& dofmap2()const{return dofmap2_;};

      
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
     typename std::enable_if< (N==0&&0==sizeof...(OtherArgs)),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
     tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     { 
     return std::tuple<OtherArg>(add_costant(otherarg,0));}

template<Integer N,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< (0<N&&0==sizeof...(OtherArgs)),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
     tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     { 
     return std::tuple<OtherArg>(add_costant(otherarg,tot_n_dofs<N-1>()));}


template<Integer N,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< 0<N && 0<sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
     tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     {
     return std::tuple_cat(std::tuple<OtherArg>( add_costant(otherarg,tot_n_dofs<N-1>()) ),
       tuple_make<N+1,OtherArgs...>(otherargs...));}

template<Integer N,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< 0==N && 0<sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
     tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     {
     return std::tuple_cat(std::tuple<OtherArg>( otherarg ),
       tuple_make<N+1,OtherArgs...>(otherargs...));}

     MixedSpace(const Arg& arg,const Args&...args):
     mesh_ptr_(arg.mesh_ptr()),
     spaces_(std::make_tuple(std::make_shared<Arg>(arg),std::make_shared<Args>(args)...)),
     n_dofs_(tot_n_dofs<1+sizeof...(Args)-1>()),
     dofmap_(tuple_make<0,typename Arg::DofMapType,typename Args::DofMapType...>(arg.dofmap(),args.dofmap()...))
     ,
     dofmap2_(std::tuple_cat(arg.dofmap2(),args.dofmap2()...))//tuple_make<0,typename Arg::DofMapType2,typename Args::DofMapType2...>(arg.dofmap2(),args.dofmap2()...))
     {}


     inline auto mesh_ptr()const {return mesh_ptr_;};

     constexpr const auto& spaces(){return spaces_;}
     static constexpr Integer Nsubspaces=tot_subspaces<Arg,Args...>::value;
     static constexpr Integer Nelem_dofs=Sum(Arg::Nelem_dofs, Args::Nelem_dofs...);
     static constexpr Array<Integer,Nsubspaces> Nelem_dofs_array=concat(Arg::Nelem_dofs_array,Args::Nelem_dofs_array...);

     static constexpr auto faces_dofs_array=std::tuple_cat(Arg ::faces_dofs_array,
                                                            Args::faces_dofs_array...);
     static constexpr auto Nfaces_dofs_array=concat(Arg::Nfaces_dofs_array,Args::Nfaces_dofs_array...);
     
   private:
    std::shared_ptr< MeshT > mesh_ptr_;
    std::tuple<std::shared_ptr<Arg>,std::shared_ptr<Args>...> spaces_;
    Integer n_dofs_;
    DofMapType dofmap_;
    DofMapType2 dofmap2_;
};

template<typename...Args>
MixedSpace<Args...> MixedFunctionSpace(const Args&...args){return MixedSpace<Args...>(args...);};













template<typename...Args>
class AuxMixedSpace; 

template<typename MeshT,typename Arg,typename...Args>
class AuxMixedSpace<AuxFunctionSpace<MeshT,Arg>,AuxFunctionSpace<MeshT,Args>...>
{
public:
  using Elem=typename AuxFunctionSpace<MeshT,Arg>::Elem;
  using TupleOfSpaces=TupleCatType<typename AuxFunctionSpace<MeshT,Arg>::TupleOfSpaces,
  typename AuxFunctionSpace<MeshT,Args>::TupleOfSpaces...>;
  using DofMapType=std::tuple<typename AuxFunctionSpace<MeshT,Arg>::DofMapType,
  typename AuxFunctionSpace<MeshT,Args>::DofMapType...>;

  using DofMapType2=TupleCatType<typename AuxFunctionSpace<MeshT,Arg>:: DofMapType2,
                                 typename AuxFunctionSpace<MeshT,Args>::DofMapType2...>;


  using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates
  <TupleCatType<typename AuxFunctionSpace<MeshT,Arg>::UniqueElementFunctionSpacesTupleType,
  typename AuxFunctionSpace<MeshT,Args>::UniqueElementFunctionSpacesTupleType...>>;
  using SpacesToUniqueFEFamily=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;
  using FromElementFunctionSpacesToUniqueNumbersTupleType=
        TupleAllToUniqueMap<TupleOfSpaces,UniqueElementFunctionSpacesTupleType>;
  using FromElementFunctionSpacesToFirstSpaceTupleType=
        FromElementFunctionSpacesToFirstSpaceTupleType<FromElementFunctionSpacesToUniqueNumbersTupleType>;

  static constexpr Integer Nuniquesubspaces=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;

  AuxMixedSpace(const AuxFunctionSpace<MeshT,Arg>& arg,const AuxFunctionSpace<MeshT,Args>&...args):
  mesh_ptr_(arg.mesh_ptr()),
  spaces_(std::make_tuple(std::make_shared<AuxFunctionSpace<MeshT,Arg>>(arg),
   std::make_shared<AuxFunctionSpace<MeshT,Args>>(args)...)),
  dofmap_(arg.dofmap(),args.dofmap()...)
  ,
  dofmap2_(arg.dofmap2(),args.dofmap2()...)
  {}

  inline auto mesh_ptr()const {return mesh_ptr_;};

  inline const DofMapType& dofmap()const{return dofmap_;};

  inline const DofMapType2& dofmap2()const{return dofmap2_;};
  
     template<Integer...Ns>
  inline constexpr const auto& dofmap()const{return tuple_get<Ns...>(dofmap_);};
  constexpr const auto& spaces(){return spaces_;}

  static constexpr Integer Nsubspaces=TupleTypeSize<TupleOfSpaces>::value;
  static constexpr Array<Integer,Nsubspaces> Nelem_dofs_array=
  concat(AuxFunctionSpace<MeshT,Arg>::Nelem_dofs_array,
    AuxFunctionSpace<MeshT,Args>::Nelem_dofs_array...);

  static constexpr auto faces_dofs_array=std::tuple_cat(AuxFunctionSpace<MeshT,Arg>::faces_dofs_array,
                                                        AuxFunctionSpace<MeshT,Args>::faces_dofs_array...);

  static constexpr auto Nfaces_dofs_array=concat(AuxFunctionSpace<MeshT,Arg>::Nfaces_dofs_array,
                                                 AuxFunctionSpace<MeshT,Args>::Nfaces_dofs_array...);
     

private:
  std::shared_ptr< MeshT > mesh_ptr_;
  std::tuple<std::shared_ptr<AuxFunctionSpace<MeshT,Arg>>,
  std::shared_ptr<AuxFunctionSpace<MeshT,Args>>...> spaces_;
  DofMapType dofmap_;
  DofMapType2 dofmap2_;
};

template<typename MeshT,typename...Args>
auto AuxFunctionSpacesBuild(const AuxFunctionSpace<MeshT,Args>&...args){return AuxMixedSpace<AuxFunctionSpace<MeshT,Args>...>(args...);};





template<typename ...Spaces>
class FullSpace;

template<typename...Args>
class FullSpace<MixedSpace<Args...>>
{
public:
  using FunctionSpace=MixedSpace<Args...>;
  using Elem=typename FunctionSpace::Elem;
  using MeshT=Mesh<Elem::Dim,Elem::ManifoldDim>;
  using TupleOfSpaces=typename FunctionSpace::TupleOfSpaces;
  using DofMapType=std::tuple<typename FunctionSpace::DofMapType>;

  using DofMapType2=std::tuple<typename FunctionSpace::DofMapType2>;

  using UniqueElementFunctionSpacesTupleType= typename FunctionSpace::UniqueElementFunctionSpacesTupleType;
  using SpacesToUniqueFEFamily=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;
  using FromElementFunctionSpacesToUniqueNumbersTupleType=
        TupleAllToUniqueMap<TupleOfSpaces,UniqueElementFunctionSpacesTupleType>;
  using FromElementFunctionSpacesToFirstSpaceTupleType=
        FromElementFunctionSpacesToFirstSpaceTupleType<FromElementFunctionSpacesToUniqueNumbersTupleType>;

  static constexpr Integer Nuniquesubspaces=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;

  static constexpr Integer TrialSpaceSize=FunctionSpace::Nsubspaces;
  static constexpr Integer Nsubspaces=FunctionSpace::Nsubspaces;
  static constexpr Array<Integer,Nsubspaces> Nelem_dofs_array=FunctionSpace::Nelem_dofs_array;
  static constexpr auto faces_dofs_array=std::tuple_cat(Args::faces_dofs_array...);
  static constexpr auto Nfaces_dofs_array=concat(Args::Nfaces_dofs_array...);


  FullSpace(const FunctionSpace& W):
  spaces_ptr_(std::make_shared<FunctionSpace>(W))
  {}

  inline auto spaces_ptr(){return spaces_ptr_;}

  inline auto& dofmap2()const{return spaces_ptr_->dofmap2();};
private:
std::shared_ptr<FunctionSpace> spaces_ptr_;
};


template<typename...Args,typename...AuxArgs>
class FullSpace<MixedSpace<Args...>,AuxMixedSpace<AuxArgs...>>
{
public:
  using FunctionSpace=MixedSpace<Args...>;
  using AuxFunctionSpace=AuxMixedSpace<AuxArgs...>;
  using Elem=typename FunctionSpace::Elem;
  using MeshT=Mesh<Elem::Dim,Elem::ManifoldDim>;
  using TupleOfSpaces=TupleCatType<typename FunctionSpace::TupleOfSpaces,
                                   typename AuxFunctionSpace::TupleOfSpaces>;
  using DofMapType=std::tuple<typename FunctionSpace::DofMapType,
                              typename AuxFunctionSpace::DofMapType>;

  using DofMapType2=std::tuple<typename FunctionSpace::DofMapType2,
                              typename AuxFunctionSpace::DofMapType2>;


  using UniqueElementFunctionSpacesTupleType= RemoveTupleDuplicates<TupleCatType<
                                     typename FunctionSpace::UniqueElementFunctionSpacesTupleType,
                                     typename AuxFunctionSpace::UniqueElementFunctionSpacesTupleType>>;
  using SpacesToUniqueFEFamily=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;
  using FromElementFunctionSpacesToUniqueNumbersTupleType=
        TupleAllToUniqueMap<TupleOfSpaces,UniqueElementFunctionSpacesTupleType>;
  using FromElementFunctionSpacesToFirstSpaceTupleType=
        FromElementFunctionSpacesToFirstSpaceTupleType<FromElementFunctionSpacesToUniqueNumbersTupleType>;


  static constexpr Integer Nuniquesubspaces=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;

  static constexpr Integer TrialSpaceSize=FunctionSpace::Nsubspaces;
  static constexpr Integer Nsubspaces=FunctionSpace::Nsubspaces + AuxFunctionSpace::Nsubspaces;
  static constexpr Array<Integer,Nsubspaces> Nelem_dofs_array=concat(FunctionSpace::Nelem_dofs_array,
                                                                     AuxFunctionSpace::Nelem_dofs_array);

  static constexpr auto faces_dofs_array=std::tuple_cat(MixedSpace<Args...>::faces_dofs_array,
                                                        AuxMixedSpace<AuxArgs...>::faces_dofs_array);
  static constexpr auto Nfaces_dofs_array=concat(MixedSpace<Args...>::Nfaces_dofs_array,
                                                 AuxMixedSpace<AuxArgs...>::Nfaces_dofs_array);

  FullSpace(const FunctionSpace& W1, const AuxFunctionSpace& W2):
  spaces_ptr_(std::make_shared<FunctionSpace>(W1)),
  aux_spaces_ptr_(std::make_shared<AuxFunctionSpace>(W2))
  {}
  inline       auto  spaces_ptr()     {return spaces_ptr_;}
  inline const auto& spaces_ptr()const{return spaces_ptr_;}
  inline       auto  aux_spaces_ptr()     {return aux_spaces_ptr_;}
  inline const auto& aux_spaces_ptr()const{return aux_spaces_ptr_;}
  
  inline auto& dofmap2()const{return spaces_ptr_->dofmap2();};
private:
std::shared_ptr<FunctionSpace> spaces_ptr_;
std::shared_ptr<AuxFunctionSpace> aux_spaces_ptr_;
};

template<typename...Args>
auto FullSpaceBuild(const MixedSpace<Args...>& W)
    {return FullSpace<MixedSpace<Args...>>(W);}

template<typename...Args,typename...AuxArgs>
auto FullSpaceBuild(const MixedSpace<Args...>& W1,const AuxMixedSpace<AuxArgs...>& W2)
    {return FullSpace<MixedSpace<Args...>,AuxMixedSpace<AuxArgs...>>(W1,W2);}









// template<typename MixedSpace, Integer N, typename OperatorType>
// class Test;

// template<typename MixedSpace, Integer N, typename OperatorType>
// class Trial;

// template<typename MixedSpace,Integer N, typename OperatorType=IdentityOperator>
// class Trial: public Expression<Trial<MixedSpace,N,OperatorType>>
// { public:
//   using FunctionSpace=MixedSpace;
//   using Elem=typename FunctionSpace::Elem;
//   static constexpr Integer Dim=Elem::Dim;
//   static constexpr Integer ManifoldDim=Elem::ManifoldDim;
//   using MeshT=Mesh<Dim,ManifoldDim>;
//   using UniqueElementFunctionSpacesTupleType=typename MixedSpace::UniqueElementFunctionSpacesTupleType;
//   static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
//   static constexpr Integer number=N;
//   static constexpr Integer value=GetType<typename MixedSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
//   using Operator=OperatorType;

//   Trial(const std::shared_ptr<MixedSpace>& W_ptr):
//   spaces_ptr_(W_ptr)
//   {}

//   Trial(const MixedSpace& W):
//   spaces_ptr_(std::make_shared<MixedSpace>(W))
//   {}

//   inline auto spaces_ptr()const {return spaces_ptr_;};

//   private:
//   std::shared_ptr<MixedSpace> spaces_ptr_;
// };

// template<typename MixedSpace, Integer N, typename OperatorType=IdentityOperator>
// class Test: public Expression<Test<MixedSpace,N,OperatorType>>
// {
// public:
//   using FunctionSpace=MixedSpace;
//   using Elem=typename FunctionSpace::Elem;
//   static constexpr Integer Dim=Elem::Dim;
//   static constexpr Integer ManifoldDim=Elem::ManifoldDim;
//   using MeshT=Mesh<Dim,ManifoldDim>;
//   using UniqueElementFunctionSpacesTupleType=typename MixedSpace::UniqueElementFunctionSpacesTupleType;
//   static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
//   static constexpr Integer number=N;//GetType<typename MixedSpace::FromSpacesToFunctionSpaces,N>::value;
//   static constexpr Integer value=GetType<typename MixedSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
//   using Operator=OperatorType;

//   Test(const std::shared_ptr<MixedSpace>& W_ptr):
//   spaces_ptr_(W_ptr)
//   {}

//   Test(const MixedSpace& W):
//   spaces_ptr_(std::make_shared<MixedSpace>(W))
//   {}

//   inline auto spaces_ptr()const {return spaces_ptr_;};

//   private:
//   std::shared_ptr<MixedSpace> spaces_ptr_;
// };


// template<typename MixedSpace, Integer N, typename Expr>
// class Trial<MixedSpace,N,CompositeOperator<Expression<Expr>>>: 
// public Expression<Trial< MixedSpace,N,CompositeOperator<Expression<Expr>> > >
// {
// public:
//   using FunctionSpace=MixedSpace;
//   using Elem=typename FunctionSpace::Elem;
//   static constexpr Integer Dim=Elem::Dim;
//   static constexpr Integer ManifoldDim=Elem::ManifoldDim;
//   using MeshT=Mesh<Dim,ManifoldDim>;
//   using UniqueElementFunctionSpacesTupleType=typename MixedSpace::UniqueElementFunctionSpacesTupleType;
//   static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
//   static constexpr Integer number=N;
//   static constexpr Integer value=GetType<typename MixedSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
//   using Operator=CompositeOperator<Expression<Expr>>;

//   Trial(const std::shared_ptr<MixedSpace>& W_ptr, const Operator& op):
//   spaces_ptr_(W_ptr),
//   op_(op)
//   {}

//   Trial(const MixedSpace& W, const Operator& op):
//   spaces_ptr_(std::make_shared<MixedSpace>(W)),
//   op_(op)
//   {}

//   inline auto spaces_ptr()const {return spaces_ptr_;};
//   inline constexpr auto composite_operator()const {return op_;};

//   private:
//   std::shared_ptr<MixedSpace> spaces_ptr_;
//   Operator op_;
// };



// template<typename MixedSpace, Integer N, typename Expr>
// class Test<MixedSpace,N,CompositeOperator<Expression<Expr>>>: 
// public Expression<Test< MixedSpace,N,CompositeOperator<Expression<Expr>> > >
// {
// public:
//   using FunctionSpace=MixedSpace;
//   using Elem=typename FunctionSpace::Elem;
//   static constexpr Integer Dim=Elem::Dim;
//   static constexpr Integer ManifoldDim=Elem::ManifoldDim;
//   using MeshT=Mesh<Dim,ManifoldDim>;
//   using UniqueElementFunctionSpacesTupleType=typename MixedSpace::UniqueElementFunctionSpacesTupleType;
//   static constexpr Integer Nmax=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
//   static constexpr Integer number=N;
//   static constexpr Integer value=GetType<typename MixedSpace::FromElementFunctionSpacesToUniqueNumbersTupleType,N>::value;
//   using Operator=CompositeOperator<Expression<Expr>>;

//   Test(const std::shared_ptr<MixedSpace>& W_ptr, const Operator& op):
//   spaces_ptr_(W_ptr),
//   op_(op)
//   {}

//   Test(const MixedSpace& W, const Operator& op):
//   spaces_ptr_(std::make_shared<MixedSpace>(W)),
//   op_(op)
//   {}

//   inline auto spaces_ptr()const {return spaces_ptr_;};
//   inline constexpr auto composite_operator()const {return op_;};
  

//   private:
//   std::shared_ptr<MixedSpace> spaces_ptr_;
//   Operator op_;
// };






































// template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, Integer N, typename Operator_,typename...OtherTemplateArguments>
// class Evaluation<Expression<TestOrTrial<MixedSpace,N,Operator_>>,OtherTemplateArguments...>
// {
//  public:
//  using type= TestOrTrial<MixedSpace,N,Operator_>;

//  static_assert((IsSame<TestOrTrial<MixedSpace,N,Operator_>,Test<MixedSpace,N,Operator_>>::value ||
//                 IsSame<TestOrTrial<MixedSpace,N,Operator_>,Trial<MixedSpace,N,Operator_>>::value )
//                && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,OperatorType>>>,TestOrTrial=Test or Trial ");
//  using FunctionSpace=typename type::FunctionSpace;
//  using FunctionSpaces=GetType<typename type::UniqueElementFunctionSpacesTupleType,type::value>;
//  using Elem=typename FunctionSpaces::Elem;
//  using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
//  using Operator=typename type::Operator;
//  using value_type=OperatorType<type,OtherTemplateArguments...>;

//  Evaluation(){};

//  Evaluation(const type& expr):
//  eval_(expr)
//  {};
 
//  //  template<typename...Forms, typename FiniteElem, typename...Inputs>
//  // constexpr void apply(value_type& value,const FiniteElem& J, const ShapeFunctionsCollection<Forms...>& shape_functions)
//  template<typename...Args, typename FiniteElem,typename...Inputs>
//  constexpr void apply(value_type& value,const FiniteElem& J, const std::tuple<Args...>& shape_functions,const Inputs&...inputs)

//  {
//   using TupleOfTupleShapeFunction=std::tuple<Args...>;
//   using tuple_type=GetType<TupleOfTupleShapeFunction,type::value>;
//   const auto& tuple=tuple_get<type::value>(shape_functions);
//     // const auto& tuple=tuple_get<type::value>(shape_functions());
//   constexpr Integer M=TypeToTupleElementPosition<ShapeFunction<Elem,BaseFunctionSpace,Operator,OtherTemplateArguments...>,tuple_type>::value;
//     // std::cout<<"befor Evaluation<Expression<TestOrTrial<MixedSpace,N,Operator_> "<<value<<std::endl;

//  // Number<M> e5e(5);

//  //  tuple_type ok(1);
//  //  decltype(tuple) ee(3);
//   Assign(value,tuple_get<M>(tuple).eval());




//   // Assignment<value_type>::apply(value,tuple_get<M>(tuple).eval());
//   std::cout<<"after Evaluation<Expression<TestOrTrial<MixedSpace,N,Operator_> "<<value<<std::endl;

//  }

// private:
 
//  type eval_;
// };









// template<template<class,Integer,class > class TestOrTrial,  typename MixedSpace, Integer N, typename Expr,typename...OtherTemplateArguments>
// class Evaluation<Expression<TestOrTrial<MixedSpace,N,CompositeOperator<Expression<Expr>>>>,OtherTemplateArguments...>
// {
//  public:
//  using Operator=CompositeOperator<Expression<Expr>>;
//  using type= TestOrTrial<MixedSpace,N,Operator>;

//  static_assert((IsSame<type,Test<MixedSpace,N,Operator>>::value ||
//                 IsSame<type,Trial<MixedSpace,N,Operator>>::value )
//                && "In Evaluation<Expression<TestOrTrial<MixedSpace,N,Composite<Operator>>>>,TestOrTrial=Test or Trial ");
//  using FunctionSpace=typename type::FunctionSpace;
//  using FunctionSpaces=GetType<typename type::UniqueElementFunctionSpacesTupleType,type::value>;
//  using Elem=typename FunctionSpaces::Elem;
//  using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
//  // using value_type=OperatorType<type,OtherTemplateArguments...>;

//  Evaluation(){};

//  Evaluation(const type& expr):
//  eval_(expr)
//  {};
 
//  //  template<typename ValueType,typename...Forms, typename FiniteElem, typename...Inputs>
//  // constexpr void apply(ValueType& value,const FiniteElem& J, const ShapeFunctionsCollection<Forms...>& shape_functions)
//   template<typename ValueType, typename FiniteElem, typename...Args1,typename...Args2,typename...Args3>
//  constexpr void apply(ValueType& value,const FiniteElem& J, const std::tuple<Args1...>& tuple_shape_functions,const std::tuple<Args2...>&tuple_tensor,const std::tuple<Args3...>&tuple_composite_shapes)

//  {
//  using single_type=Evaluation<Expression<typename FormOfCompositeOperatorType<type>::type >,OtherTemplateArguments...>;

//  using TupleOfTupleCompositeShapeFunctionEval=std::tuple<Args3...>;
//  using tuple_type=GetType<TupleOfTupleCompositeShapeFunctionEval,type::value>;
//  // using single_type=Evaluation<Expression<typename FormOfCompositeOperatorType<type>::type >,OtherTemplateArguments...>;
//  // const auto& tuple=tuple_get<type::value>(tuple_tensor);
//  constexpr Integer M=TypeToTupleElementPosition<single_type,tuple_type>::value;

//  auto& single_tensor=tuple_get<type::value,M >(tuple_tensor);
//  Assign(value,single_tensor);


//  // Assignment<ValueType>::apply(value,single_tensor);
//    // std::cout<<"single_tensor= "<<single_tensor<<std::endl;

//   // single_type ok1(1);
//   // TupleOfTupleCompositeShapeFunctionEval ok2(3);
//   std::cout<<"Evaluation<Expression<TestOrTrial<MixedSpace,N,Operator_> "<<value<<std::endl;


//   }

// private:
 
//  type eval_;
// };












 



























// template<typename FuncType,Integer N,typename Operator_,typename FullSpace,typename...OtherTemplateArguments>
// class Evaluation<Expression<Function<FullSpace,N,Operator_,FuncType>>,OtherTemplateArguments...>
// {
//  public:
//  using type=Function<FullSpace,N,Operator_,FuncType>;
//  using FunctionSpace=typename type::FunctionSpace;
//  using FunctionSpaces=GetType<typename type::UniqueElementFunctionSpacesTupleType,type::value>;
//  using Elem=typename FunctionSpaces::Elem;
//  using BaseFunctionSpace=Elem2FunctionSpace<FunctionSpaces>;
//  using Operator=typename type::Operator;
//  using value_type=OperatorType<type,OtherTemplateArguments...>;
//  Evaluation(){};

//  // template<typename ...Ts> 
//  Evaluation(const type& expr):
//  eval_(expr)
//  {};

//  template<typename Shapes>
//  void compute(value_type& value, type&eval,const FiniteElem<Elem>&J, const Shapes& shapes )
//  {
//   using type1=typename Shapes::type; //Vector<Vector<T,NQpoints>,Ndofs>
//   using type2=typename Shapes::subtype; //Vector<T,NQpoints>
//   using type3=typename Shapes::subtype::subtype; // T 
  
//    // std::cout<<"pre eval.local_dofs_update(J)"<<std::endl;
//    eval.local_dofs_update(J);
//    // std::cout<<"post eval.local_dofs_update(J)"<<std::endl;
//    const auto& local_dofs=eval.local_dofs();
//     // loop on qp points
//     for(Integer ii=0;ii< type2::Dim;ii++)
//     {
//       for(Integer mm=0;mm<type3::Rows;mm++)
//         for(Integer nn=0;nn<type3::Cols;nn++)
//            value[ii](mm,nn)=local_dofs[0]*shapes[0][ii](mm,nn);
//     // loop on dofs
//      for(Integer jj=1;jj< type1::Dim;jj++)
//       // loop on components of subtype (is a matrix)
//       for(Integer mm=0;mm<type3::Rows;mm++)
//         for(Integer nn=0;nn<type3::Cols;nn++)
//            value[ii](mm,nn)+=local_dofs[jj]*shapes[jj][ii](mm,nn);

//     }


//  }
 
//  // template<typename...Forms, typename FiniteElem,typename...Inputs>
//  // constexpr void apply(value_type& value,const FiniteElem& J, const typename ShapeFunctionsCollection<Forms...>::TupleOfTupleShapeFunction& shape_functions,const Inputs&...inputs)
//  template<typename...Args, typename FiniteElem,typename...Inputs>
//  constexpr void apply(value_type& value,const FiniteElem& J, const std::tuple<Args...>& shape_functions,const Inputs&...inputs)

//  {
//   using TupleOfTupleShapeFunction=std::tuple<Args...>;
//   using tuple_type=GetType<TupleOfTupleShapeFunction,type::value>;
//   // const auto& tuple=tuple_get<type::value>(shape_functions());
//   const auto& tuple=tuple_get<type::value>(shape_functions);

//   // using tuple_type=GetType<std::tuple<Args...>,type::value>;
//   // const auto& tuple=tuple_get<type::value>(tuple_of_tuple);
//   constexpr Integer M=TypeToTupleElementPosition<ShapeFunction<Elem,BaseFunctionSpace,Operator,OtherTemplateArguments...>,tuple_type>::value;
//   // FIXME


//   const auto& shapes=tuple_get<M>(tuple).eval();
//   compute(value,eval_,J,shapes);

//   // decltype(tuple) oo(5);
//   // Number<M> owo(5);
//   // ShapeFunction<Elem,BaseFunctionSpace,Operator,OtherTemplateArguments...> rra(6);
//   // tuple_type ee(5);
//   // TupleOfTupleShapeFunction r5(6);
//   // std::tuple<Args...> kj(6);
//   std::cout<<"Evaluation<Expression<Function "<<value<<std::endl;


//  }


// private:
//  type eval_;
// };












// template<typename ConstType,typename...Inputs,typename QuadratureRule>
// class Evaluation<Expression<ConstantTensor<ConstType,Inputs...>>,QuadratureRule>
// {
//  public:
//  using type=ConstantTensor<ConstType,Inputs...>;
//  using Input=typename type::type;
//  using value_type=OperatorType<type,QuadratureRule>;


//  Evaluation(){};

//  Evaluation(const type& expr):
//  eval_(expr)
//  { };



 
//  template<typename...Args, typename FiniteElem,typename...OtherInputs>
//  constexpr void apply(value_type& value,const FiniteElem& J, const std::tuple<Args...>& tuple_of_tuple,const OtherInputs&...inputs)
//  {
//   // TODO here we copy the static static_value  into value, but it is useless. better to specialize evaluation
//   // Assign(value,eval_.template eval<QuadratureRule::NQpoints>());
//   // value_type oo(4,4,4,44,4,4);
//   // decltype(eval_.template qp_eval<QuadratureRule::NQPoints>()) rrr(4,4,4,44,4,4);
//   Assignment<value_type>::apply(value,eval_.template qp_eval<QuadratureRule::NQPoints>());
//   std::cout<<"Evaluation<Expression<ConstantTensor "<<value<<std::endl;



//  }


// private:
//  type eval_;
// };







// template<typename T>
// class IsVolumeOrSurfaceIntegral
// {
// public:
//   static constexpr bool surface=true;
//   static constexpr bool volume=true;
// };

// template<template<class,Integer,class>class TestOrTrial, typename MixedSpace,Integer N>
// class IsVolumeOrSurfaceIntegral<TestOrTrial<MixedSpace,N,TraceOperator>>{
// public:
//   static constexpr bool surface=true;
//   static constexpr bool volume=false;
// };

// template<template<class,Integer,class>class TestOrTrial, typename MixedSpace,Integer N, typename OperatorType>
// class IsVolumeOrSurfaceIntegral<TestOrTrial<MixedSpace,N,OperatorType>>{
// public:
//   static constexpr bool surface=false;
//   static constexpr bool volume=true;
// };

// template<template<class>class Unary, typename Type>
// class IsVolumeOrSurfaceIntegral< Unary<Expression<Type>> >
// {
// public:
//   static constexpr bool surface=IsVolumeOrSurfaceIntegral<Type>::surface;
//   static constexpr bool volume=IsVolumeOrSurfaceIntegral<Type>::volume;
// };


// template<template<class,class>class Binary, typename Left,typename Right>
// class IsVolumeOrSurfaceIntegral< Binary<Expression<Left>,Expression<Right>> >
// {
// public:
//   static constexpr bool surface=IsVolumeOrSurfaceIntegral<Left>::surface*IsVolumeOrSurfaceIntegral<Right>::surface;
//   static constexpr bool volume=IsVolumeOrSurfaceIntegral<Left>::volume*IsVolumeOrSurfaceIntegral<Right>::volume;
// };


// template<typename Left,typename Right>
// class IsVolumeOrSurfaceIntegral<InnerProduct<Expression <Left>, Expression <Right> >  >
// {
// public:
//   static constexpr bool surface=IsVolumeOrSurfaceIntegral<Left>::surface*IsVolumeOrSurfaceIntegral<Right>::surface;
//   static constexpr bool volume=IsVolumeOrSurfaceIntegral<Left>::volume*IsVolumeOrSurfaceIntegral<Right>::volume;
// };

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

// template<typename T>
// class IsTestOrTrial;

// // template<typename T>
// // class IsTestOrTrial{
// // public:
// //   // using Elem=EmptyClass;
// //   // using Operator=std::tuple<>;
// //   // using TupleFunctionSpace=std::tuple<>;
// //   // using UniqueElementFunctionSpacesTupleType=std::tuple<>;
// //   using type=std::tuple<Number<-1>>;
// //   // static constexpr Integer value=-1;
// //   // static constexpr Integer number=-1;
// // };


// template<typename ConstType,typename...Inputs>
// class IsTestOrTrial< ConstantTensor<ConstType,Inputs...>>{
// public:
//   using Elem=EmptyClass;
//   using Operator=std::tuple<>;
//   using TupleFunctionSpace=std::tuple<>;
//   using UniqueElementFunctionSpacesTupleType=std::tuple<>;
//   using type=std::tuple<Number<0>>;
//   static constexpr Integer value=-1;
//   static constexpr Integer number=-1;
// };
// template<typename ConstType,typename...Inputs>
// class IsTestOrTrial<const ConstantTensor<ConstType,Inputs...>>
// : public IsTestOrTrial< ConstantTensor<ConstType,Inputs...>>
// {};


// template<typename FuncType,typename FullSpace, Integer N, typename Operator_>
// class IsTestOrTrial<Function<FullSpace,N,Operator_,FuncType>>{
// public:
//   using Elem=typename FullSpace::Elem;
//   using Operator=std::tuple<Operator_>;
//   using TupleFunctionSpace=std::tuple<>;
//   using UniqueElementFunctionSpacesTupleType=std::tuple<>;
//   using type=std::tuple<Number<0>>;
//   static constexpr Integer value=-1;
//   static constexpr Integer number=-1;
// };

// template<typename MixedSpace,Integer N, typename OperatorType>
// class IsTestOrTrial<Test<MixedSpace,N,OperatorType>>{
// public:
//   using Elem=typename MixedSpace::Elem;
//   using Operator=std::tuple<OperatorType>;
//   using TupleFunctionSpace=std::tuple<MixedSpace>;
//   using UniqueElementFunctionSpacesTupleType=std::tuple<typename MixedSpace::UniqueElementFunctionSpacesTupleType>;
//   using type=std::tuple<Number<1>>;
//   static constexpr Integer value=Test<MixedSpace,N,OperatorType>::value;
//   static constexpr Integer number=Test<MixedSpace,N,OperatorType>::number;
// };

// template<typename MixedSpace,Integer N, typename OperatorType>
// class IsTestOrTrial<Trial<MixedSpace,N,OperatorType>>{
// public:
//   using Elem=typename MixedSpace::Elem;
//   using Operator=std::tuple<OperatorType>;
//   using TupleFunctionSpace=std::tuple<MixedSpace>;
//   using UniqueElementFunctionSpacesTupleType=std::tuple<typename MixedSpace::UniqueElementFunctionSpacesTupleType>;
//   using type=std::tuple<Number<2>>;
//   static constexpr Integer value=Trial<MixedSpace,N,OperatorType>::value;
//   static constexpr Integer number=Trial<MixedSpace,N,OperatorType>::number;
// };

// template <typename T,typename ... Types>
// class TupleTypeSize;




// template<template<class>class Unary, typename Type>
// class IsTestOrTrial< Unary<Expression<Type>> >
// {
// public:
//   using Elem=typename IsTestOrTrial<Type>::Elem;

//   using Operator=TupleCatType<typename IsTestOrTrial<Type>::Operator>;
//   using TupleFunctionSpace=typename IsTestOrTrial<Type>::TupleFunctionSpace;
//   using UniqueElementFunctionSpacesTupleType=typename IsTestOrTrial<Type>::UniqueElementFunctionSpacesTupleType;
//   // WE REMOVE THE ZEROS, BUT WE SHOULD NOT REMOVE ALL ZEROS
//   // if it has only one element, do not remove duplicates
//   using type=TupleRemoveNumber0<typename IsTestOrTrial<Type>::type>;
//   // using type=TupleRemoveType<Number<0>,typename IsTestOrTrial<Type>::type>;
//   static_assert(TupleTypeSize<type>::value<2," In Unary<Type>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
//   static constexpr Integer number= Heaviside(IsTestOrTrial<Type>::number);
// };


// template<template<class,class>class Binary, typename Left, typename Right>
// class IsTestOrTrial<Binary<Expression<Left>,Expression<Right> > >
// {
// public:

//   using Elem=Choose<typename IsTestOrTrial<Left>::Elem,typename IsTestOrTrial<Right>::Elem>;
//   using Operator=TupleCatType<typename IsTestOrTrial<Left>::Operator,
//                               typename IsTestOrTrial<Right>::Operator >;
//   using TupleFunctionSpace=TupleCatType<typename IsTestOrTrial<Left>::TupleFunctionSpace,
//                                         typename IsTestOrTrial<Right>::TupleFunctionSpace >;

//   using UniqueElementFunctionSpacesTupleType=TupleCatType<typename IsTestOrTrial<Left>::UniqueElementFunctionSpacesTupleType,
//                                                           typename IsTestOrTrial<Right>::UniqueElementFunctionSpacesTupleType >;
//   using tmp_type=TupleRemoveType<Number<0>,TupleCatType<typename IsTestOrTrial<Left>::type,typename IsTestOrTrial<Right>::type >> ;
//   using type=typename std::conditional<IsSame<tmp_type,std::tuple<>>::value, 
//                                        std::tuple<Number<0>>,
//                                        tmp_type>::type;   
//   static constexpr Integer number= Heaviside(IsTestOrTrial<Left>::number)+Heaviside(IsTestOrTrial<Right>::number);
//   static_assert(TupleTypeSize<type>::value<2," In Binary<Left,Right>, cannot have more than one test/more than one trial/one ore more trials and one or more tests");
// };








// template<typename Left,typename Right>
// class TypeOfForm;


// template<Integer M,Integer N>
// class TypeOfForm<Number<M>,Number<N>>
// {
// public:
//   using type=void;
//   static_assert(0*Number<N>::value,"L2inner: the form is neither a 0-form(function,function), 1-form(function/test,test/function) or 2-form (test/trial,trial/test), where the function is neither a test nor a trial");
// };



// template<>
// class TypeOfForm<Number<0>,Number<0>>
// {
//   public:
//     using type=Number<0>; 
// };

// template<>
// class TypeOfForm<Number<0>,Number<1>>
// {
//   public:
//     using type=Number<1>; 
// };

// template<>
// class TypeOfForm<Number<1>,Number<0>>
// {
//   public:
//     using type=Number<1>; 
// };

// template<>
// class TypeOfForm<Number<1>,Number<2>>
// {
//   public:
//     using type=Number<2>; 
// };


// template<>
// class TypeOfForm<Number<2>,Number<1>>
// {
//   public:
//     using type=Number<2>; 
// };






// template<Integer N,typename...Args >
// auto MakeTest(const FunctionSpace<Args...>& W)
// {return Test<FunctionSpace<Args...>,N>(W);}


// template<Integer N,typename...Args >
// auto MakeTrial(const FunctionSpace<Args...>& W)
// {return Trial<FunctionSpace<Args...>,N>(W);}









// template<Integer N,typename...Args >
// auto MakeTest(const MixedSpace<Args...>& W)
// {return Test<MixedSpace<Args...>,N>(W);}


// template<Integer N,typename...Args >
// auto MakeTrial(const MixedSpace<Args...>& W)
// {return Trial<MixedSpace<Args...>,N>(W);}


// template<Integer N,typename...Args >
// auto MakeTest(const FullSpace<Args...>& W)
// {return Test<FullSpace<Args...>,N>(W);}


// template<Integer N,typename...Args >
// auto MakeTrial(const FullSpace<Args...>& W)
// {return Trial<FullSpace<Args...>,N>(W);}



// template<typename MixedSpace,Integer N>
// Trial<MixedSpace,N,GradientOperator> 
// Grad(const Trial<MixedSpace,N,IdentityOperator>& t)
// {return Trial<MixedSpace,N,GradientOperator> (t.spaces_ptr());}

// template<typename MixedSpace,Integer N>
// Test<MixedSpace,N,GradientOperator> 
// Grad(const Test<MixedSpace,N,IdentityOperator>& t)
// {return Test<MixedSpace,N,GradientOperator> (t.spaces_ptr());}

// template<typename MixedSpace,Integer N>
// Trial<MixedSpace,N,DivergenceOperator> 
// Div(const Trial<MixedSpace,N,IdentityOperator>& t)
// {return Trial<MixedSpace,N,DivergenceOperator> (t.spaces_ptr());}

// template<typename MixedSpace,Integer N>
// Test<MixedSpace,N,DivergenceOperator> 
// Div(const Test<MixedSpace,N,IdentityOperator>& t)
// {return Test<MixedSpace,N,DivergenceOperator> (t.spaces_ptr());}


// template<typename MixedSpace,Integer N>
// Trial<MixedSpace,N,CurlOperator> 
// Curl(const Trial<MixedSpace,N,IdentityOperator>& t)
// {return Trial<MixedSpace,N,CurlOperator> (t.spaces_ptr());}

// template<typename MixedSpace,Integer N>
// Test<MixedSpace,N,CurlOperator> 
// Curl(const Test<MixedSpace,N,IdentityOperator>& t)
// {return Test<MixedSpace,N,CurlOperator> (t.spaces_ptr());}



// template<typename MixedSpace,Integer N>
// Trial<MixedSpace,N,TraceOperator> 
// Trace(const Trial<MixedSpace,N,IdentityOperator>& t)
// {return Trial<MixedSpace,N,TraceOperator> (t.spaces_ptr());}

// template<typename MixedSpace,Integer N>
// Test<MixedSpace,N,TraceOperator> 
// Trace(const Test<MixedSpace,N,IdentityOperator>& t)
// {return Test<MixedSpace,N,TraceOperator> (t.spaces_ptr());}






// template<Integer...N>
// class FormTestTrialNumbers;

// template<Integer LeftSpaceNumber,Integer RightSpaceNumber>
// class FormTestTrialNumbers<2,2,1,LeftSpaceNumber,RightSpaceNumber>
// {
// public:
//   using type=std::tuple<Number<RightSpaceNumber>,Number<LeftSpaceNumber>>;
// };

// template<Integer LeftSpaceNumber,Integer RightSpaceNumber>
// class FormTestTrialNumbers<2,1,2,LeftSpaceNumber,RightSpaceNumber>
// {
// public:
//   using type=std::tuple<Number<LeftSpaceNumber>,Number<RightSpaceNumber>>;
// };

// template<Integer LeftSpaceNumber,Integer RightSpaceNumber>
// class FormTestTrialNumbers<1,0,1,LeftSpaceNumber,RightSpaceNumber>
// {
// public:
//   using type=std::tuple<Number<RightSpaceNumber>>;
// };

// template<Integer LeftSpaceNumber,Integer RightSpaceNumber>
// class FormTestTrialNumbers<1,1,0,LeftSpaceNumber,RightSpaceNumber>
// {
// public:
//   using type=std::tuple<Number<LeftSpaceNumber>>;
// };























// template<typename T1,typename T2>
// class TupleOfTestTrialPairsNumbersAuxAux;

// template<typename T>
// class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<std::tuple<>> >
// {
//  public:
//   using type=std::tuple<>;
// };



// // template<typename T, typename MeshT, typename Left,typename Right,Integer QR>
// // class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<L2DotProductIntegral<MeshT,Left,Right,QR>> >
// // {
// //  public:
// //   using L2=L2DotProductIntegral<MeshT,Left,Right,QR>;
// //   using type=typename std::conditional<IsSame<T,typename L2::TestTrialNumbers>::value, L2, std::tuple<>>::type;
// // };
// template<typename T, typename Left,typename Right,Integer QR>
// class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<L2DotProductIntegral<Left,Right,QR>> >
// {
//  public:
//   using L2=L2DotProductIntegral<Left,Right,QR>;
//   using type=typename std::conditional<IsSame<T,typename L2::TestTrialNumbers>::value, L2, std::tuple<>>::type;
// };




// template<typename T1, typename T2>
// class TupleOfTestTrialPairsNumbersAuxAux<T1, Expression<Addition<Expression<T2>,Expression<std::tuple<>> > >>
// {
//  public:
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T1,Expression<T2>>::type;
// };


// template<typename T1, typename T2>
// class TupleOfTestTrialPairsNumbersAuxAux<T1, Expression<Addition<Expression<std::tuple<>>, Expression<T2> > >>
// {
//  public:
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T1,Expression<T2>>::type;
// };

// template<typename T>
// class TupleOfTestTrialPairsNumbersAuxAux<T, Expression<Addition<Expression<std::tuple<>>, Expression<std::tuple<>> > >>
// {
//  public:
//   using type=std::tuple<>;
// };


// // template<typename T, typename MeshT1, typename Left1,typename Right1,Integer QR1,
// //                      typename MeshT2, typename Left2,typename Right2,Integer QR2>
// // class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<MeshT1,Left1,Right1,QR1>>,
// //                                                                  Expression<L2DotProductIntegral<MeshT2,Left2,Right2,QR2>>>>>
// // {
// //  public:
// //   using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<MeshT1,Left1,Right1,QR1>>>::type>,
// //                        Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<MeshT2,Left2,Right2,QR2>>>::type>>;
// // };
// template<typename T, typename Left1,typename Right1,Integer QR1,
//                      typename Left2,typename Right2,Integer QR2>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//                                                                Expression<L2DotProductIntegral<Left2,Right2,QR2>>>>>
// {
//  public:
//   using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left1,Right1,QR1>>>::type>,
//                        Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left2,Right2,QR2>>>::type>>;
// };


// template<typename T, typename Left1,typename Right1,Integer QR1,
//                      typename Right2>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//                                                                Expression<Right2>>>>
// {
//  public:
//   using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left1,Right1,QR1>>>::type>,
//                        Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Right2>>::type>>;
// };
// template<typename T, typename Left1,typename Right1,Integer QR1>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//                                                                Expression<std::tuple<>>>>>
// {
//  public:
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left1,Right1,QR1>>>::type;
                       
// };



// template<typename T, typename Left1,
//                      typename Left2,typename Right2,Integer QR2>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<Left1>,
//                                                                Expression<L2DotProductIntegral<Left2,Right2,QR2>>>>>
// {
//  public:
//   using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Left1>>::type>,
//                        Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left2,Right2,QR2>>>::type>>;
// };

// template<typename T, typename Left2,typename Right2,Integer QR2>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<std::tuple<>>,
//                                                                Expression<L2DotProductIntegral<Left2,Right2,QR2>>>>>
// {
//  public:
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<L2DotProductIntegral<Left2,Right2,QR2>>>::type;
// };


// template<typename T, typename Left,typename Right>
// class TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Addition<Expression<Left>,
//                                                                Expression<Right> >>>
// {
//  public:
//   using type=Addition<Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Left>>::type>,
//                       Expression<typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<Right>>::type>>;
// };



// template<typename...Ts>
// class TupleOfTestTrialPairsNumbersAux;

// // template<typename T,typename MeshT, typename Left,typename Right,Integer QR>
// // class TupleOfTestTrialPairsNumbersAux<T,L2DotProductIntegral<MeshT,Left,Right,QR> >
// // {
// //  public:
// //   using S=L2DotProductIntegral<MeshT,Left,Right,QR>;
// //   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<S>>::type;

// // };
// template<typename T,typename Left,typename Right,Integer QR>
// class TupleOfTestTrialPairsNumbersAux<T,L2DotProductIntegral<Left,Right,QR> >
// {
//  public:
//   using S=L2DotProductIntegral<Left,Right,QR>;
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,Expression<S>>::type;

// };

// // template<typename T, typename MeshT1, typename Left1,typename Right1,Integer QR1,
// //                      typename MeshT2, typename Left2,typename Right2,Integer QR2>
// // class TupleOfTestTrialPairsNumbersAux<T,Addition<Expression<L2DotProductIntegral<MeshT1,Left1,Right1,QR1>>,
// //                         Expression<L2DotProductIntegral<MeshT2,Left2,Right2,QR2>>>>
// // {
// //  public:
// //   using Left=L2DotProductIntegral<MeshT1,Left1,Right1,QR1>;
// //   using Right=L2DotProductIntegral<MeshT2,Left2,Right2,QR2>;
// //   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
// //                           Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
// //                                                 Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
// //                             >::type;
// // };
// template<typename T, typename Left1,typename Right1,Integer QR1,
//                      typename Left2,typename Right2,Integer QR2>
// class TupleOfTestTrialPairsNumbersAux<T,Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//                         Expression<L2DotProductIntegral<Left2,Right2,QR2>>>>
// {
//  public:
//   using Left=L2DotProductIntegral<Left1,Right1,QR1>;
//   using Right=L2DotProductIntegral<Left2,Right2,QR2>;
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
//                           Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
//                                                 Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
//                             >::type;
// };


// // template<typename T,typename MeshT, typename Left1,typename Right1,Integer QR,typename Right>
// // class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR>>,Expression<Right > > >
// // {
// //  public:
// //   using Left=L2DotProductIntegral<MeshT,Left1,Right1,QR>;
// //   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
// //                           Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
// //                                                 Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
// //                              >::type;

// // };

// template<typename T,typename Left1,typename Right1,Integer QR,typename Right>
// class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<L2DotProductIntegral<Left1,Right1,QR>>,Expression<Right > > >
// {
//  public:
//   using Left=L2DotProductIntegral<Left1,Right1,QR>;
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
//                           Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
//                                                 Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
//                              >::type;

// };

// // template<typename T, typename Left,typename MeshT, typename Left1,typename Right1,Integer QR>
// // class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<Left>,Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR> > > >
// // {
// //  public:
// //   using Right=L2DotProductIntegral<MeshT,Left1,Right1,QR>;
// //   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
// //                              Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
// //                                                    Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
// //                              >::type;
// // };
// template<typename T, typename Left, typename Left1,typename Right1,Integer QR>
// class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1,QR> > > >
// {
//  public:
//   using Right=L2DotProductIntegral<Left1,Right1,QR>;
//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
//                              Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
//                                                  Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
//                              >::type;
// };

// template<typename T, typename Left,typename Right>
// class TupleOfTestTrialPairsNumbersAux<T, Addition<Expression<Left>,Expression<Right > > >
// {
// public:

//   using type=typename TupleOfTestTrialPairsNumbersAuxAux<T,
//                             Expression<Addition<Expression<typename TupleOfTestTrialPairsNumbersAux<T,Left>::type>,
//                                                 Expression<typename TupleOfTestTrialPairsNumbersAux<T,Right>::type>>>
//                              >::type;
// };



// template<typename...Ts>
// class TupleOfTestTrialPairsNumbers;

// // template<typename MeshT1, typename Left1,typename Right1,Integer QR1>
// // class TupleOfTestTrialPairsNumbers<L2DotProductIntegral<MeshT1,Left1,Right1,QR1>>
// // {
// //  public:
// //   using T=L2DotProductIntegral<MeshT1,Left1,Right1,QR1>;

// //   using type=std::tuple<typename T::TestTrialNumbers>;

// // };

// template<typename Left1,typename Right1,Integer QR1>
// class TupleOfTestTrialPairsNumbers<L2DotProductIntegral<Left1,Right1,QR1>>
// {
//  public:
//   using T=L2DotProductIntegral<Left1,Right1,QR1>;

//   using type=std::tuple<typename T::TestTrialNumbers>;

// };

// // template<typename MeshT1, typename Left1,typename Right1,Integer QR1,
// //          typename MeshT2, typename Left2,typename Right2,Integer QR2>
// // class TupleOfTestTrialPairsNumbers<Addition<Expression<L2DotProductIntegral<MeshT1,Left1,Right1,QR1>>,
// //                         Expression<L2DotProductIntegral<MeshT2,Left2,Right2,QR2>>>>
// // {
// //  public:
// //   using Left=L2DotProductIntegral<MeshT1,Left1,Right1,QR1>;
// //   using Right=L2DotProductIntegral<MeshT2,Left2,Right2,QR2>;

// //   using type=RemoveTupleDuplicates<std::tuple<typename Left::TestTrialNumbers, typename Right::TestTrialNumbers>>;

// // };

// template<typename Left1,typename Right1,Integer QR1,
//          typename Left2,typename Right2,Integer QR2>
// class TupleOfTestTrialPairsNumbers<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
//                         Expression<L2DotProductIntegral<Left2,Right2,QR2>>>>
// {
//  public:
//   using Left=L2DotProductIntegral<Left1,Right1,QR1>;
//   using Right=L2DotProductIntegral<Left2,Right2,QR2>;

//   using type=RemoveTupleDuplicates<std::tuple<typename Left::TestTrialNumbers, typename Right::TestTrialNumbers>>;

// };

// // template<typename MeshT, typename Left1,typename Right1,Integer QR,typename Right>
// // class TupleOfTestTrialPairsNumbers<Addition<Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR>>,Expression<Right > > >
// // {
// //  public:
// //   using Left=L2DotProductIntegral<MeshT,Left1,Right1,QR>;
// //   using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type >>;
// // };

// template<typename Left1,typename Right1,Integer QR,typename Right>
// class TupleOfTestTrialPairsNumbers<Addition<Expression<L2DotProductIntegral<Left1,Right1,QR>>,Expression<Right > > >
// {
//  public:
//   using Left=L2DotProductIntegral<Left1,Right1,QR>;
//   using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type >>;
// };

// // template<typename Left,typename MeshT, typename Left1,typename Right1,Integer QR>
// // class TupleOfTestTrialPairsNumbers<Addition<Expression<Left>,Expression<L2DotProductIntegral<MeshT,Left1,Right1,QR> > > >
// // {
// //  public:
// //   using Right=L2DotProductIntegral<MeshT,Left1,Right1,QR>;
// //   using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type >>;
// // };

// template<typename Left,typename Left1,typename Right1,Integer QR>
// class TupleOfTestTrialPairsNumbers<Addition<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1,QR> > > >
// {
//  public:
//   using Right=L2DotProductIntegral<Left1,Right1,QR>;
//   using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type >>;
// };

// template<typename Left,typename Right>
// class TupleOfTestTrialPairsNumbers<Addition<Expression<Left>,Expression<Right > > >
// {
// public:
//   using type=RemoveTupleDuplicates<TupleCatType<typename TupleOfTestTrialPairsNumbers<Left>::type,typename TupleOfTestTrialPairsNumbers<Right>::type>>;
// };




// template<typename TupleOfPairsNumbers, typename Form,Integer Nmax,Integer N>
// class TupleOfL2ProductsHelper;

// template<typename TupleOfPairsNumbers, typename Form,Integer Nmax>
// class TupleOfL2ProductsHelper<TupleOfPairsNumbers,Form,Nmax,Nmax>
// {
//  public:
//  using type=std::tuple<typename TupleOfTestTrialPairsNumbersAux<GetType<TupleOfPairsNumbers,Nmax>,Form>::type>;
// };



// template<typename TupleOfPairsNumbers, typename Form,Integer Nmax,Integer N>
// class TupleOfL2ProductsHelper
// {
//  public:
//  using type=TupleCatType<std::tuple<typename TupleOfTestTrialPairsNumbersAux<GetType<TupleOfPairsNumbers,N>,Form>::type> , 
//                          typename TupleOfL2ProductsHelper<TupleOfPairsNumbers,Form,Nmax,N+1>::type >;
// };

// template<typename TupleOfPairsNumbers, typename Form>
// class TupleOfL2Products
// {
// public:
// using type=  typename TupleOfL2ProductsHelper<TupleOfPairsNumbers,Form, TupleTypeSize<TupleOfPairsNumbers>::value-1,0>::type;


// };




// template<typename Left_,typename Right_,Integer QR=GaussianQuadrature>
// class L2DotProductIntegral: 
// public Expression<L2DotProductIntegral<Left_,Right_,QR>>
// {  
//    public:
//     using Left=Left_;
//     using Right=Right_;
//     using type=InnerProduct<Expression <Left>, Expression <Right> > ;
//     using TestOrTrialLeft= IsTestOrTrial<Left>;
//     using TestOrTrialRight= IsTestOrTrial<Right>;
//     using OperatorLeft=typename TestOrTrialLeft::Operator;
//     using OperatorRight=typename TestOrTrialRight::Operator;
//     static constexpr Integer leftN=TestOrTrialLeft::number;
//     static constexpr Integer rightN=TestOrTrialRight::number;
   
//    static_assert(IsVolumeOrSurfaceIntegral<type>::volume && "In Volume integrals, no trace operator can occur");
    

//     static constexpr Integer TestOrTrialLeftValue =GetType<typename TestOrTrialLeft::type,0>::value;
//     static constexpr Integer TestOrTrialRightValue =GetType<typename TestOrTrialRight::type,0>::value;

//     using Elem=Choose<typename TestOrTrialLeft::Elem,typename TestOrTrialRight::Elem>;
//     static constexpr Integer Order=CheckMaxQuadratureOrder<Elem,QR,QuadratureOrder<type>::value+1>::value; 
//     using QRule=typename QuadratureRule<QR>:: template rule<Elem,Order>;


//     using form= std::tuple<typename TypeOfForm<GetType<typename IsTestOrTrial<Left>::type,0>,
//                                                GetType<typename IsTestOrTrial<Right>::type,0>
//                             >::type >;
//     using TestTrialNumbers=typename FormTestTrialNumbers<GetType<form,0>::value,TestOrTrialLeftValue,TestOrTrialRightValue,leftN,rightN>::type;

//     using UniqueElementFunctionSpacesTupleType=GetType<RemoveTupleDuplicates< TupleCatType< typename TestOrTrialLeft::UniqueElementFunctionSpacesTupleType,
//                                                                                             typename TestOrTrialRight::UniqueElementFunctionSpacesTupleType  >>,0>;               
//     using TupleFunctionSpace=RemoveTupleDuplicates< TupleCatType< typename IsTestOrTrial<Left>::TupleFunctionSpace,
//                                                              typename IsTestOrTrial<Right>::TupleFunctionSpace  >>;               
    
//     using FunctionSpace=GetType<TupleFunctionSpace,0>;

//     using TupleOfSpaces=typename FunctionSpace::TupleOfSpaces;


//     template<typename TupleOfSpacesAux,typename TestTrialNumbersAux,Integer N>
//     class ClassAux;

//     template<typename TupleOfSpacesAux,typename TestTrialNumbersAux>
//     class ClassAux<TupleOfSpacesAux,TestTrialNumbersAux,1>
//     {
//     public:
//       using type_test= typename std::conditional< (-1==GetType<TestTrialNumbersAux,0>::value),
//                                           std::tuple<>,
//                                           GetType<TupleOfSpacesAux,GetType<TestTrialNumbersAux,0>::value>>::type;
//       using type=std::tuple<type_test>;
//     };

//     template<typename TupleOfSpacesAux,typename TestTrialNumbersAux>
//     class ClassAux<TupleOfSpacesAux,TestTrialNumbersAux,2>
//     {
//     public:
//       using type_test= typename std::conditional< (-1==GetType<TestTrialNumbersAux,0>::value),
//                                           std::tuple<>,
//                                           GetType<TupleOfSpacesAux,GetType<TestTrialNumbersAux,0>::value>>::type;
//       using type_trial=typename std::conditional< (-1==GetType<TestTrialNumbersAux,1>::value),
//                                           std::tuple<>,
//                                           GetType<TupleOfSpacesAux,GetType<TestTrialNumbersAux,1>::value>>::type;
//       using type=std::tuple<type_test,type_trial>;
//     };

//    using TestTrialSpaces=typename ClassAux<TupleOfSpaces,TestTrialNumbers,GetType<form,0>::value>::type;
   
//     const Left&  left() const{return left_;};
//     const Right& right()const{return right_;};

//     L2DotProductIntegral(const Expression<Left>& left,const Expression<Right>& right,const Integer label=-666)
//     :
//     left_(left.derived()),
//     right_(right.derived()),
//     product_(Inner(left,right)),
//     label_(label)
//     {}
     

//     L2DotProductIntegral(const Expression<L2DotProductIntegral<Left,Right,QR>>& l2prod,const Integer label=-666)
//     :
//     left_(l2prod.derived().left()),
//     right_(l2prod.derived().right()),
//     product_(Inner(left_,right_)),
//     label_(label)
//     {}


//      auto operator()(){return L2DotProductIntegral<Left,Right,QR>(left_,right_);}
//      const auto operator()()const{return L2DotProductIntegral<Left,Right,QR>(left_,right_);}
//   private:
//     Left left_;
//     Right right_;
//     type product_;
//     Integer label_;
// };













//////////////////////////////////////////////////////////////////////////////////////////////////////
////// This class is used to add the corresponding FunctionSpace of a Function (aka Type_), 
////// which is neither Test nor Trial, but which is projected onto a finite element space,
////// to the tuple containing all the necessary FunctionSpaces of the form
////// N_ is the counter of Function already added to Tuple_ (we start from 0)
////// Given Type_, for sure a new function has been found, so N=N_+1
////// Then we search if the corresponding FunctionSpace is already present in Tuple_
////// (Example Tuple_=tuple<Lagrange1,Lagrange2>, Type_=RT0)
////// If not, we add it to the tuple and we say that the Function position is N
////// If it is, Tuple is unchanged and the Function position is actually position
//////////////////////////////////////////////////////////////////////////////////////////////////////





























































































































                               






// template<typename Left,typename Right>
// auto
// L2Inner(const Expression<Left>& left,const Expression<Right>& right)
// {return L2DotProductIntegral<Left,Right>(left,right);}

// // template<typename Left,typename Right>
// // auto
// // L2Inner(const Expression<Left>& left,const Expression<Right>& right, const Integer label)
// // {return L2DotProductIntegral<Left,Right>(left,right,label);}


// template<typename Left,typename Right>
// constexpr auto
// L2Inner(const Expression<UnaryMinus<Expression<Left>>>& left,const Expression<UnaryMinus<Expression<Right>>>& right)
// {return L2DotProductIntegral<Left,Right>(left.derived().derived(),right.derived().derived());}

// template<typename Left,typename Right>
// constexpr auto
// L2Inner(const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Left>>> >>& left,const Expression<Right>& right)
// {return L2DotProductIntegral<Left,Right>(left.derived().derived().derived(),right);}


// template<typename Left,typename Right>
// constexpr auto
// L2Inner(const Expression<Left>& left,const  Expression<UnaryMinus<Expression<UnaryMinus<Expression<Right>>> >>& right)
// {return L2DotProductIntegral<Left,Right>(left,right.derived().derived().derived());}

// template<typename Left,typename Right>
// constexpr auto
// L2Inner(const Expression<UnaryMinus<Expression<Left>>>& left,const  Expression<UnaryMinus<Expression<UnaryMinus<Expression<Right>>> >>& right)
// {return L2DotProductIntegral<UnaryMinus<Expression<Left>>,Right>(left.derived(),right.derived().derived().derived());}

// template<typename Left,typename Right>
// constexpr auto
// L2Inner(const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Left>>> >>& left,const Expression<UnaryMinus<Expression<Right>>>& right)
// {return L2DotProductIntegral<Left,UnaryMinus<Expression<Right>>>(left.derived().derived().derived(),right.derived());}

// template<typename Left,typename Right>
// constexpr auto
// L2Inner(const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Left>>>>>& left,const Expression<UnaryMinus<Expression<UnaryMinus<Expression<Right>>> >>& right)
// {return L2DotProductIntegral<Left,Right>(left.derived().derived().derived(),right.derived().derived().derived());}




// template<typename Left2,typename Right2, typename Left>
// constexpr auto
// L2Inner(const Expression<Left>& left,const Expression<Addition<Expression<Left2>,Expression<Right2>>>& right)
// {return 
//   Addition<
//   Expression<decltype(L2Inner(left.derived(),right.derived().left()))>,
//   Expression<decltype(L2Inner(left.derived(),right.derived().right()))>>
//               (L2Inner(left.derived(),right.derived().left()),
//                L2Inner(left.derived(),right.derived().right()) );}

// template<typename Left2,typename Right2, typename Left>
// constexpr auto
// L2Inner(const Expression<Left>& left,const Expression<Subtraction<Expression<Left2>,Expression<Right2>>>& right)
// {return Addition<
//   Expression<decltype(L2Inner(left.derived(),right.derived().left()))>,
//   Expression<decltype(L2Inner(-left.derived(),right.derived().right()))>>
//               (L2Inner(left.derived(),right.derived().left()),
//                L2Inner(-left.derived(),right.derived().right()) );}







// template<typename Left1,typename Right1, typename Right>
// constexpr auto
// L2Inner(const Expression<Addition<Expression<Left1>,Expression<Right1>>>& left,const Expression<Right>& right)
// {return Addition<
//   Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
//   Expression<decltype(L2Inner(left.derived().right(),right.derived()))>>
//   (L2Inner(left.derived().left(),right.derived()),
//    L2Inner(left.derived().right(),right.derived()) );}

// template<typename Left1,typename Right1, typename Right>
// constexpr auto
// L2Inner(const Expression<Subtraction<Expression<Left1>,Expression<Right1>>>& left,const Expression<Right>& right)
// {return Addition<
//   Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
//   Expression<decltype(L2Inner(-left.derived().right(),right.derived()))>>
//   (L2Inner(left.derived().left(),right.derived()),
//    L2Inner(-left.derived().right(),right.derived()) );}


// template<typename Left1,typename Right1,typename Left2, typename Right2>
// constexpr auto
// L2Inner(const Expression<Addition<Expression<Left1>,Expression<Right1>>>& left,
//         const Expression<Addition<Expression<Left2>,Expression<Right2>>>& right)
// {return 
// Addition<Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
//          Expression<decltype(L2Inner(left.derived().right(),right.derived()))>
//          >
//   (
//     L2Inner(left.derived().left(),right.derived()),
//     L2Inner(left.derived().right(),right.derived())                
//   )
//   ;
// }


// template<typename Left1,typename Right1,typename Left2, typename Right2>
// constexpr auto
// L2Inner(const Expression<Addition<Expression<Left1>,Expression<Right1>>>& left,
//         const Expression<Subtraction<Expression<Left2>,Expression<Right2>>>& right)
// {return 
// Addition<Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
//          Expression<decltype(L2Inner(left.derived().right(),right.derived()))>
//          >
//   (
//     L2Inner(left.derived().left(),right.derived()),
//     L2Inner(left.derived().right(),right.derived())                
//   )
//   ;
// }

// template<typename Left1,typename Right1,typename Left2, typename Right2>
// constexpr auto
// L2Inner(const Expression<Subtraction<Expression<Left1>,Expression<Right1>>>& left,
//         const Expression<Subtraction<Expression<Left2>,Expression<Right2>>>& right)
// {return 
// Addition<Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
//          Expression<decltype(L2Inner(-left.derived().right(),right.derived()))>
//          >
//   (
//     L2Inner(left.derived().left(),right.derived()),
//     L2Inner(-left.derived().right(),right.derived())                
//   )
//   ;
// }


// template<typename Left1,typename Right1,typename Left2, typename Right2>
// constexpr auto
// L2Inner(const Expression<Subtraction<Expression<Left1>,Expression<Right1>>>& left,
//         const Expression<Addition<Expression<Left2>,Expression<Right2>>>& right)
// {return 
// Addition<Expression<decltype(L2Inner(left.derived().left(),right.derived()))>,
//          Expression<decltype(L2Inner(-left.derived().right(),right.derived()))>
//          >
//   (
//     L2Inner(left.derived().left(),right.derived()),
//     L2Inner(-left.derived().right(),right.derived())                
//   )
//   ;
// }

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //////  + L2DotProductIntegral
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename Left,typename Right, Integer QR>
// class L2DotProductIntegral<Left,Right,QR >
// operator+(const Expression<L2DotProductIntegral<Left,Right,QR>>&l2prod)
// {return L2Inner(l2prod.derived().left(),l2prod.derived().right());}



// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //////  - L2DotProductIntegral
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename Left,typename Right, Integer QR>
// class L2DotProductIntegral<UnaryMinus<Expression<Left>>,Right,QR >
// operator-(const Expression<L2DotProductIntegral<Left,Right,QR>>&l2prod)
// {return L2Inner(-l2prod.derived().left(),l2prod.derived().right());}



// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //////  Exr - L2DotProductIntegral
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// template< typename Left1,typename Left2,typename Right2, Integer QR>
// class Addition< Expression <Left1>, 
//                 Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR>> >
// operator-(const Expression<Left1>&left, 
//           const Expression<L2DotProductIntegral<Left2,Right2,QR>>&right)
// {return left+L2Inner(-right.derived().left(),right.derived().right());}




// template< typename Left1,typename Left2,typename Right2, Integer QR>
// class Addition< Expression <Left1>, 
//                 Expression<L2DotProductIntegral<Left2,Right2,QR>> >
// operator-(const Expression<Left1>&left, 
//           const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR>>&right)
// {return left+L2Inner(right.derived().left().derived(),right.derived().right());}



// template< typename Left1,typename Left2,typename Right2, Integer QR2>
// class Addition< Expression<Left1>, 
//                 Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// operator-(const Expression<Left1>&left, 
//           const Expression<L2DotProductIntegral<Left2,UnaryMinus<Expression<Right2>>,QR2>>&right)
// {return left+L2Inner(right.derived().left(),right.derived().right().derived());}




// template< typename Left1,typename Left2,typename Right2, Integer QR2>
// class Addition< Expression<Left1>, 
//                 Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>> >
// operator-(const Expression<Left1>&left, 
//           const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,UnaryMinus<Expression<Right2>>,QR2>>&right)
// {return left+L2Inner(right.derived().left(),right.derived().right().derived());}


// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //////  L2DotProductIntegral - Exr
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// // template< typename Left1,typename Right1, Integer QR, typename Right2>
// // class Addition< Expression<L2DotProductIntegral<Left1,Right1,QR>>, Expression<UnaryMinus<Expression<Right2>>> >
// // operator-(const Expression<L2DotProductIntegral<Left1,Right1,QR>>&left, 
// //           const Expression<Right2>&right)
// // {return Addition<Expression<L2DotProductIntegral<UnaryMinus<Expression<Left1>>,Right1,QR>>,
// //                  Expression <UnaryMinus<Expression<Right2>>> >
// //                  (left,-right);}

// // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// // //////  L2DotProductIntegral - L2DotProductIntegral
// // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// // template< typename Left1,typename Right1, Integer QR1, typename Left2,typename Right2, Integer QR2>
// // class Addition< Expression<L2DotProductIntegral<Left1,Right1,QR1>>, 
// //                 Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>> >
// // operator-(const Expression<L2DotProductIntegral<Left1,Right1,QR1>>&left, 
// //           const Expression<L2DotProductIntegral<Left2,Right2,QR2>>&right)
// // {return left+L2Inner(-right.derived().left(),right.derived().right());}
// //   // return Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
// //   //                Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>> >
// //   //                (left,L2Inner(-right.derived().left(),right.derived().right()));}





// // template< typename Left1,typename Right1, Integer QR1, typename Left2,typename Right2, Integer QR2>
// // class Addition< Expression<L2DotProductIntegral<Left1,Right1,QR1>>, 
// //                 Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// // operator-(const Expression<L2DotProductIntegral<Left1,Right1,QR1>>&left, 
// //           const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>>&right)
// // {
// // return left+L2Inner(right.derived().left().derived(),right.derived().right());
// //   // return Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
// //   //                Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// //   //                (left,L2Inner(right.derived().left().derived(),right.derived().right()));
// //                }

// // template< typename Left1,typename Left2,typename Right2, Integer QR2>
// // class Addition< Expression<Left1>, 
// //                 Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// // operator-(const Expression<Left1>&left, 
// //           const Expression<L2DotProductIntegral<Left2,UnaryMinus<Expression<Right2>>,QR2>>&right)
// // {return left+L2Inner(right.derived().left(),right.derived().right().derived());
// //   // Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
// //   //                Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// //   //                (left,L2Inner(right.derived().left().derived(),right.derived().right()));
// //                }


// // template< typename Left1,typename Left2,typename Right2, Integer QR2>
// // class Addition< Expression<Left1>, 
// //                 Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>> >
// // operator-(const Expression<Left1>&left, 
// //           const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,UnaryMinus<Expression<Right2>>,QR2>>&right)
// // {return left+L2Inner(right.derived().left(),right.derived().right().derived());
// //   // Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
// //   //                Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// //   //                (left,L2Inner(right.derived().left().derived(),right.derived().right()));
// //                }






// //   Addition<Expression<Left1>,
// //                  Expression<L2DotProductIntegral<UnaryMinus<Expression<Left2>>,Right2,QR2>> >
// //                  (left,L2Inner(right.derived().left().derived(),right.derived().right().derived()));}



// // template< typename Left1,typename Right1, Integer QR1, typename Left2,typename Right2, Integer QR2>
// // class Addition< Expression<L2DotProductIntegral<Left1,Right1,QR1>>, 
// //                 Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// // operator-(const Expression<L2DotProductIntegral<UnaryMinus<Expression<Left1>>,Right1,QR1>>&left, 
// //           const Expression<L2DotProductIntegral<Left2,Right2,QR2>>&right)
// // {return Addition<Expression<L2DotProductIntegral<Left1,Right1,QR1>>,
// //                  Expression<L2DotProductIntegral<Left2,Right2,QR2>> >
// //                  (left,L2Inner(right.derived().left(),right.derived().right().derived()));}




// // // template<typename Left1,typename Right1, typename Right>
// // // Addition<Expression<L2DotProductIntegral<MeshT,Left1,Right>>,
// // //           Expression<L2DotProductIntegral<MeshT,Right1,Right>>>
// // // L2Inner(const Expression<Addition<Expression<Left1>,Expression<Right1>>>& left,const Expression<Right>& right, const Integer& label)
// // // {return Addition<
// // //   Expression<L2DotProductIntegral<MeshT,Left1,Right>>,
// // //   Expression<L2DotProductIntegral<MeshT,Right1,Right>>>(L2Inner(left.derived().left(),right,label),
// // //                                                         L2Inner(left.derived().right(),right,label) );}












/////////////// GUARDA QUI TODO FIX LA HO COMMENTATA MA BOH
  // template<typename Left,typename Left1,typename Right1, Integer Kind>
  // class ClassHelper<Subtraction<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1>> >,Kind  >
  // {
  // public:
  //  using type=RemoveTupleDuplicates< TupleCatType< typename ClassHelper<Left,Kind>::type,
  //                                                  typename KindType<L2DotProductIntegral<Left1,Right1>,Kind>::type  >>;
  // } ; 

























// template<typename DerivedLeft, typename MeshT, typename Left,typename Right>
// class Addition< Expression<DerivedLeft>  ,  Expression<L2DotProductIntegral<MeshT,Left,Right>>>
// : public Expression< Addition< Expression<DerivedLeft>  ,  Expression<L2DotProductIntegral<MeshT,Left,Right>>  > > 
// {

//   public:
//   using TupleFunctionSpace=RemoveTupleDuplicates< TupleCatType< typename L2DotProductIntegral<MeshT,Left,Right>::TupleFunctionSpace,
//                                                            typename DerivedLeft::TupleFunctionSpace  >>;               

//   using FunctionSpace=GetType<TupleFunctionSpace,0>;

//   using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::UniqueElementFunctionSpacesTupleType,
//                                                                                  typename DerivedLeft::UniqueElementFunctionSpacesTupleType > >;
//   using form =RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::form,
//                                                   typename DerivedLeft::form  > >;
//     Addition(const Expression<DerivedLeft>& left, const Expression<L2DotProductIntegral<MeshT,Left,Right>>&right)
//     : 
//     left_(left.derived()),
//     right_(right.derived())
//     {};
//     const DerivedLeft & left()const{return left_;};
//     const L2DotProductIntegral<MeshT,Left,Right>& right()const{return right_;};
//   private:
//   DerivedLeft left_;
//   L2DotProductIntegral<MeshT,Left,Right> right_;
// };


// template<typename DerivedLeft, typename MeshT, typename Left,typename Right>
// class Subtraction< Expression<DerivedLeft>  ,  Expression<L2DotProductIntegral<MeshT,Left,Right>>>
// : public Expression< Subtraction< Expression<DerivedLeft>  ,  Expression<L2DotProductIntegral<MeshT,Left,Right>>  > > 
// {

//   public:
//   using TupleFunctionSpace=RemoveTupleDuplicates< TupleCatType< typename L2DotProductIntegral<MeshT,Left,Right>::TupleFunctionSpace,
//                                                            typename DerivedLeft::TupleFunctionSpace  >>;               

//   using FunctionSpace=GetType<TupleFunctionSpace,0>;

//   using UniqueElementFunctionSpacesTupleType=RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::UniqueElementFunctionSpacesTupleType,
//                                                                                  typename DerivedLeft::UniqueElementFunctionSpacesTupleType > >;
//   using form =RemoveTupleDuplicates< TupleCatType<typename L2DotProductIntegral<MeshT,Left,Right>::form,
//                                                   typename DerivedLeft::form  > >;
//     Subtraction(const Expression<DerivedLeft>& left, const Expression<L2DotProductIntegral<MeshT,Left,Right>>&right)
//     : 
//     left_(left.derived()),
//     right_(right.derived())
//     {};
//     const DerivedLeft & left()const{return left_;};
//     const L2DotProductIntegral<MeshT,Left,Right>& right()const{return right_;};
//   private:
//   DerivedLeft left_;
//   L2DotProductIntegral<MeshT,Left,Right> right_;
// };



 


// template<typename Form_>
// class GeneralForm<Form_>
// {

//   public:
//   using Form=Form_;
//   template<typename T,Integer N>
//   class KindType;

//   template<typename T>
//   class KindType<T,0>
//   {
//   public:
//     using type=typename T::TupleFunctionSpace;
//   };


//  template<typename T>
//   class KindType<T,1>
//   {
//   public:
//     using type=typename T::UniqueElementFunctionSpacesTupleType;
//   };


//  template<typename T>
//   class KindType<T,2>
//   {
//   public:
//     using type=typename T::form;
//   };

//   template<typename FormTmp, Integer Kind>
//   class ClassHelper;

//   template<typename Left,typename Right, Integer Kind>
//   class ClassHelper<L2DotProductIntegral<Left,Right>,Kind  >
//   {
//   public:
//    using type=RemoveTupleDuplicates< typename KindType<L2DotProductIntegral<Left,Right>,Kind>::type>;
//   }; 

//   template<typename Left1,typename Right1, typename Left2,typename Right2, Integer Kind>
//   class ClassHelper<Addition<Expression<L2DotProductIntegral<Left1,Right1>>,
//                              Expression<L2DotProductIntegral<Left2,Right2>> >,Kind  >
//   {
//   public:
//    using type=RemoveTupleDuplicates< TupleCatType< typename KindType<L2DotProductIntegral<Left1,Right1>,Kind>::type,
//                                                    typename KindType<L2DotProductIntegral<Left2,Right2>,Kind>::type >>;
//   }; 


//   template<typename Left1,typename Right1, typename Left2,typename Right2, Integer Kind>
//   class ClassHelper<Subtraction<Expression<L2DotProductIntegral<Left1,Right1>>,
//                                 Expression<L2DotProductIntegral<Left2,Right2>> >,Kind  >
//   {
//   public:
//    using type=RemoveTupleDuplicates< TupleCatType< typename KindType<L2DotProductIntegral<Left1,Right1>,Kind>::type,
//                                                    typename KindType<L2DotProductIntegral<Left2,Right2>,Kind>::type >>;
//   }; 
 
//   template<typename Left1,typename Right1, typename Right, Integer Kind>
//   class ClassHelper<Addition<Expression<L2DotProductIntegral<Left1,Right1>>,Expression<Right> >,Kind  >
//   {
//   public:
//    using type=RemoveTupleDuplicates< TupleCatType< typename KindType<L2DotProductIntegral<Left1,Right1>,Kind>::type,
//                                                    typename ClassHelper<Right,Kind>::type >>;
//   };  

 
//   template<typename Left1,typename Right1, typename Right, Integer Kind>
//   class ClassHelper<Subtraction<Expression<L2DotProductIntegral<Left1,Right1>>,Expression<Right> >,Kind  >
//   {
//   public:
//    using type=RemoveTupleDuplicates< TupleCatType< typename KindType<L2DotProductIntegral<Left1,Right1>,Kind>::type,
//                                                    typename ClassHelper<Right,Kind>::type >>;
//   };  

  
//   template<typename Left,typename Left1,typename Right1, Integer Kind>
//   class ClassHelper<Addition<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1>> >,Kind  >
//   {
//   public:
//    using type=RemoveTupleDuplicates< TupleCatType< typename ClassHelper<Left,Kind>::type,
//                                                    typename KindType<L2DotProductIntegral<Left1,Right1>,Kind>::type  >>;
//   }; 
 

//   template<typename Left,typename Left1,typename Right1, Integer Kind>
//   class ClassHelper<Subtraction<Expression<Left>,Expression<L2DotProductIntegral<Left1,Right1>> >,Kind  >
//   {
//   public:
//    using type=RemoveTupleDuplicates< TupleCatType< typename ClassHelper<Left,Kind>::type,
//                                                    typename KindType<L2DotProductIntegral<Left1,Right1>,Kind>::type  >>;
//   } ; 

//   template<typename Left,typename Right, Integer Kind>
//   class ClassHelper<Addition<Expression<Left>,Expression<Right> >,Kind >
//   {
//   public:
//    using type=RemoveTupleDuplicates< TupleCatType< typename ClassHelper<Left,Kind>::type,
//                                                    typename ClassHelper<Right,Kind>::type  >>;
//   };   

//   template<typename Left,typename Right, Integer Kind>
//   class ClassHelper<Subtraction<Expression<Left>,Expression<Right> >,Kind  >
//   {
//   public:
//    using type=RemoveTupleDuplicates< TupleCatType< typename ClassHelper<Left,Kind>::type,
//                                                    typename ClassHelper<Right,Kind>::type  >>;
//   };   


//   template<Integer Kind>
//   using type=typename ClassHelper<Form,Kind>::type;

//   using TupleFunctionSpace=typename ClassHelper<Form,0>::type;      

//   using FunctionSpace=GetType<TupleFunctionSpace,0>;  

//   using UniqueElementFunctionSpacesTupleType=typename ClassHelper<Form,1>::type;      

//   using form=typename ClassHelper<Form,2>::type;      

//   using TupleOfPairsNumbers=BubbleSortTupleOfPairsNumbers<typename TupleOfTestTrialPairsNumbers<Form>::type>;

//     GeneralForm(const Form& form)
//     : 
//     form_(form)
//     {};

//     // GeneralForm(const Form& form,const FullSpace& space)
//     // : 
//     // form_(form),
//     // spaces_ptr_(std::make_shared<FullSpace>(space))
//     // {};

//     const Form& operator()()const{return form_;};
//           Form& operator()()     {return form_;};

//     // const std::shared_ptr<FullSpace>& spaces_ptr()const{return spaces_ptr_;};
//     //       std::shared_ptr<FullSpace>& spaces_ptr()     {return spaces_ptr_;};

//   private:
//   Form form_;
//   // std::shared_ptr<FullSpace> spaces_ptr_;
// };
   


// template<typename Form>
// constexpr auto general_form(const Form& form){return GeneralForm<Form>(form);}









































// template<typename Left_,typename Right_,Integer QR, typename...Forms>
// class Evaluation<Expression<L2DotProductIntegral<Left_,Right_,QR>>, ShapeFunctionsCollection<Forms...>>
// {
//  public:
//  using Left=Left_;
//  using Right=Right_;
//  using type= L2DotProductIntegral<Left,Right,QR>;
//  using QRule=typename type ::QRule;
//  using TestTrialSpaces=typename type::TestTrialSpaces;
//  using subtype= OperatorType<type,GetType<typename type::form>>;
//  static constexpr bool PositiveWeights= IsPositive(QRule::qp_weights);
 
//  Evaluation(const type& expr, ShapeFunctionsCollection<Forms...>& shape_functions):
//  expr_(expr.left(),expr.right())
//  ,
//  shape_functions_(shape_functions)
//  ,
//  local_tensor_(expr)
//  {};
 
//  template<typename Elem>
//  void apply_aux(subtype& mat, const FiniteElem<Elem>& J)
//  {

//   // changed todo fixme
//    // local_tensor_.apply(mat,J,shape_functions_());
//   std::cout<<"pre Evaluation<Expression<L2DotProductIntegral local tensor="<<std::endl;
  

//   local_tensor_.apply(mat,J,shape_functions_);//(),shape_functions_.composite_tensor(),shape_functions_.composite_shapes());
//   std::cout<<"after Evaluation<Expression<L2DotProductIntegral local tensor="<<std::endl;
//  }


//  template<typename...Inputs>
//  void apply(subtype& mat, const Inputs&...inputs)
//  {apply_aux(mat,inputs...);}

//  auto operator()(){return local_tensor_;}

//        auto  expression(){return expr_;}
//  const auto& expression()const{return expr_;}



// private:
//  type expr_;
//  ShapeFunctionsCollection<Forms...>& shape_functions_;
// };

























































// template<typename GeneralForm_, typename...GeneralForms_>
// class ShapeFunctionsCollection
// {
//  public:
//   using GeneralForm=GeneralForm_;
//   using Elem=typename GeneralForm::FunctionSpace::Elem;
//   using Form=MultipleAddition<typename GeneralForm_::Form,typename GeneralForms_::Form...> ;  
//   using UniqueElementFunctionSpacesTupleType=typename GeneralForm::UniqueElementFunctionSpacesTupleType;
//   using TupleOperatorsAndQuadrature= typename OperatorAndQuadratureTupleType<Form>::type;
//   using TupleOfTupleNoQuadrature=TupleOfTupleRemoveQuadrature<TupleOperatorsAndQuadrature>;
//   using SpacesToUniqueFEFamilies=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;

//   using Map=MapOperatorTupleOfTuple<TupleOfTupleNoQuadrature,UniqueElementFunctionSpacesTupleType>;
//   using UniqueMapping=UniqueMap<SpacesToUniqueFEFamilies,Map> ;
//   using TupleOfTupleShapeFunction=TupleOfTupleShapeFunctionType2<UniqueElementFunctionSpacesTupleType,TupleOperatorsAndQuadrature>;
  
//   using TupleCompositeOperatorsAndQuadrature= typename OperatorAndQuadratureTupleType<Form>::composite_type;
//   // using TupleOfTupleShapeFunctionCombinations=TupleOfTupleShapeFunctionType3<UniqueElementFunctionSpacesTupleType,
//   //                                                                            TupleCompositeOperatorsAndQuadrature>;

//   using MapTupleNumbers=MapTupleInit<TupleOfTupleShapeFunction, SpacesToUniqueFEFamilies, UniqueMapping>;
  
//   using TupleOfTupleCompositeShapeFunction=typename TupleOfCombinationFunctions<GeneralForm::FunctionSpace::Nuniquesubspaces,Form>::type;
//   using TupleOfTupleCompositeShapeFunctionTensor=typename TupleOfCombinationFunctions<GeneralForm::FunctionSpace::Nuniquesubspaces,Form>::type_tensor;

//    template<Integer M,Integer Nmax_aux,Integer N,typename Tuple,typename Map>
//    constexpr typename std::enable_if< (N>Nmax_aux) ,void>::type 
//    init_map_aux_aux(Tuple& t,const Map& maps){}

//    template<Integer M,Integer Nmax_aux,Integer N,typename Tuple,typename Map>
//    constexpr typename std::enable_if< (N<=Nmax_aux) ,void>::type 
//    init_map_aux_aux(Tuple& t,const Map& maps) 
//    {
//     auto& t_nth=std::get<N>(t); 
//     auto& map_nth=std::get<GetType<MapTupleNumbers,M,N>::value>(maps); 
//     t_nth.init_map(map_nth);
//     init_map_aux_aux<M,Nmax_aux,N+1>(t,maps);
//     }



//    template<Integer Nmax_aux,Integer N,typename Map>
//    constexpr typename std::enable_if< (N>Nmax_aux) ,void>::type 
//    init_map_aux(const Map& maps){}

//    template<Integer Nmax_aux,Integer N,typename Map>
//    constexpr typename std::enable_if< (N<=Nmax_aux) ,void>::type 
//    init_map_aux(const Map& maps) 
//    {
//     auto& t_nth=std::get<N>(tuple); 
//     auto& map_nth=std::get<GetType<SpacesToUniqueFEFamilies,N>::value>(maps); 
//     init_map_aux_aux<N,TupleTypeSize<decltype(t_nth)>::value-1,0>(t_nth,map_nth);
//     init_map_aux<Nmax_aux,N+1>(maps);
//     }


//    template<typename...Args>
//    constexpr void init_map(const ReferenceMaps2<Args...>& maps) 
//    {init_map_aux<TupleTypeSize<MapTupleNumbers>::value-1,0>(maps());}



//   template<Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
//   constexpr typename std::enable_if_t< FEFamily==LagrangeFE,void> 
//   shape_function_init_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficient<Args...> &coefficients)
//   {
//     std::cout<<M<<" "<<FEFamily<<" "<<Order<<std::endl;
//     shape.init();
//   }


//   template<Integer M, typename Elem, Integer FEFamily,Integer Order,typename Shape,typename...Args>
//   constexpr typename std::enable_if_t<FEFamily==RaviartThomasFE,void> 
//   shape_function_init_aux_aux_aux(Shape& shape, const ShapeFunctionCoefficient<Args...> &coefficients)
//   {
//     shape.init(tuple_get<M>(coefficients()));
//   }




//   template<Integer M, Integer Nmax,Integer N,typename Tuple,typename...Args>
//   constexpr typename std::enable_if_t<(N>Nmax),void> 
//   shape_function_init_aux_aux(Tuple& tuple, const ShapeFunctionCoefficient<Args...> &coefficients)
//   {}

//   template<Integer M,Integer Nmax,Integer N,typename Tuple,typename...Args>
//   constexpr typename std::enable_if_t< (N<=Nmax),void> 
//   shape_function_init_aux_aux(Tuple& tuple, const ShapeFunctionCoefficient<Args...> &coefficients)
//   {

//     using Shape=GetType<TupleOfTupleShapeFunction,M,N>;
//     using Elem=typename Shape::Elem;
//     constexpr Integer FEFamily=Shape::FEFamily;
//     constexpr Integer Order=Shape::Order;
//     shape_function_init_aux_aux_aux<M,Elem,FEFamily,Order>(tuple_get<N>(tuple),coefficients);
//     shape_function_init_aux_aux<M,Nmax,N+1>(tuple,coefficients);
//   }

//   template<Integer Nmax,Integer N,typename...Args>
//   constexpr typename std::enable_if_t< (N>Nmax),void> 
//   shape_function_init_aux(const ShapeFunctionCoefficient<Args...> &coefficients)
//   {}

//   template<Integer Nmax,Integer N,typename...Args>
//   constexpr typename std::enable_if_t<N<=Nmax,void> 
//   shape_function_init_aux(const ShapeFunctionCoefficient<Args...> &coefficients)
//   {
//     constexpr Integer Nmax_aux=TupleTypeSize<GetType<TupleOfTupleShapeFunction,N>>::value-1;
//     shape_function_init_aux_aux<N,Nmax_aux,0>(tuple_get<N>(tuple),coefficients);
//     shape_function_init_aux<Nmax,N+1>(coefficients);
//   }



//   //template<typename Coefficients>
//   constexpr void init_shape_functions()//const Coefficients& shape_coefficients)
//   {
//    constexpr Integer Nmax=TupleTypeSize<TupleOfTupleShapeFunction>::value-1;
//    shape_function_init_aux<Nmax,0>(coeffs_);
//    }





















//   template<Integer Nmax,Integer N,typename...Args1,typename...Args2>
//   constexpr typename std::enable_if_t<(N>Nmax),void> 
//   init_composite_shape_functions_aux_aux(std::tuple<Args1...>& tuple_tensor,const FiniteElem<Elem> &J, std::tuple<Args2...>& tuple_composite)
//   {
//     std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;

//   }

//   template<Integer Nmax,Integer N,typename...Args1>
//   constexpr void init_composite_shape_functions_aux_aux(std::tuple<>& tuple_tensor,const FiniteElem<Elem> &J, std::tuple<Args1...>& tuple_composite)
//   {
//     std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;
//   }

//   template<Integer Nmax,Integer N,typename...Args1>
//   constexpr void init_composite_shape_functions_aux_aux(std::tuple<Args1...>& tuple_tensor,const FiniteElem<Elem> &J, std::tuple<>& tuple_composite)
//   {
//     std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;

//   }

//   template<Integer Nmax,Integer N>
//   constexpr void init_composite_shape_functions_aux_aux(std::tuple<>& tuple_tensor,const FiniteElem<Elem> &J, std::tuple<>& tuple_composite)
//   {
//     std::cout<<"aux aux void Nmax,N="<<Nmax<<", "<<N<<std::endl;

//   }

//   template<Integer Nmax,Integer N,typename...Args1,typename...Args2>
//   constexpr typename std::enable_if_t< (N<=Nmax),void> 
//   init_composite_shape_functions_aux_aux(std::tuple<Args1...>& tuple_tensor,const FiniteElem<Elem> &J, std::tuple<Args2...>& tuple_composite)
//   {
    
//     auto& composite=tuple_get<N>(tuple_composite);
//     auto& tensor=tuple_get<N>(tuple_tensor);
//     // decltype(composite) a1(1);
//     // decltype(tensor) a2(1);
//     // decltype(tuple_tensor)d344(5);
//         // TupleOfTupleCompositeShapeFunction ee4e(6);

//     // TupleOfTupleCompositeShapeFunctionTensor eee(6);
//     // Number<N> ee3(5);
//     // Number<Nmax> ee(5);
//     std::cout<<"aux aux Nmax,N="<<Nmax<<", "<<N<<std::endl;

//     // todo fixme scommentami
//     composite.apply(tensor,J,tuple);
//     init_composite_shape_functions_aux_aux<Nmax,N+1>(tuple_tensor,J,tuple_composite);
//   }


//   template<Integer Nmax,Integer N>
//   constexpr typename std::enable_if_t< (N>Nmax),void> 
//   init_composite_shape_functions_aux(const FiniteElem<Elem> &J)
//   {}

//   template<Integer Nmax,Integer N>
//   constexpr typename std::enable_if_t<N<=Nmax,void> 
//   init_composite_shape_functions_aux(const FiniteElem<Elem> &J)
//   {
//     auto& tensor=tuple_get<N>(tuple_tensors_);
//     const auto& composite=tuple_get<N>(tuple_composite_);
//     constexpr Integer Nmax_aux=TupleTypeSize<decltype(composite)>::value-1;
//     constexpr Integer Nmax_au2=TupleTypeSize<decltype(tensor)>::value-1;

//     // Number<Nmax_aux> rr(6);

//     // Number<Nmax_aux> e(6);
//     // Number<Nmax_au2> r5e(6);
//     // decltype(tensor) rf(56);
//     // decltype(composite) r4f(56);
//     // std::cout<<"aux Nmax,N="<<Nmax<<", "<<N<<std::endl;
//     // std::cout<<"Nmax_aux="<<TupleTypeSize<GetType<TupleOfTupleCompositeShapeFunction,N>>::value<<" "<<Nmax_aux<<std::endl;
//     // std::cout<<"Nmax_aux="<<TupleTypeSize<decltype(tensor)>::value-1<<" "<<Nmax_aux<<std::endl;
//     // std::cout<<"Nmax_aux="<<TupleTypeSize<decltype(composite)>::value-1<<" "<<Nmax_aux<<std::endl;

//     init_composite_shape_functions_aux_aux<Nmax_aux,0>(tuple_get<N>(tuple_tensors_),J,tuple_get<N>(tuple_composite_));
//     init_composite_shape_functions_aux<Nmax,N+1>(J);
//   }

//   //template<typename Coefficients>
//   constexpr void init_composite_shape_functions(const FiniteElem<Elem> &J)
//   {
//     std::cout<<"init_composite_shape_functions"<<std::endl;
//    constexpr Integer Nmax=TupleTypeSize<TupleOfTupleCompositeShapeFunction>::value-1;
//    init_composite_shape_functions_aux<Nmax,0>(J);
//    }



//   constexpr void init(const FiniteElem<Elem>&J)
//   {
//     // TODO CHECK: SINCE WE HAVE ALREADY HAVE REFERENCES TO MAPS AND COEFFS, do we have to init_map?
//     // INIT MAP: takes the corresponding map and compute the shape function in the actual element
//     std::cout<<"init maps"<<std::endl;
//     init_map(maps_);
//     // init_shape_functions: takes also the coefficients which multiply the actual shape functions
//     std::cout<<"init_shape_functions"<<std::endl;
//     init_shape_functions();
//     std::cout<<"init_composite_shape_functions"<<std::endl;
//     init_composite_shape_functions(J);
//    }




//    constexpr       auto& operator()()      {return tuple;}
//    constexpr const auto& operator()()const {return tuple;}


//    constexpr       auto & composite_shapes()      {return tuple_composite_;}
//    constexpr const auto & composite_shapes()const {return tuple_composite_;}

//    constexpr       auto & composite_tensor()      {return tuple_tensors_;}
//    constexpr const auto & composite_tensor()const {return tuple_tensors_;}

//    template<Integer...Ns>
//    constexpr const auto& get()const{return tuple_get<Ns...>(tuple);}

//    template<Integer...Ns>
//          auto& get()     {return tuple_get<Ns...>(tuple);}

//    template<Integer...Ns>
//    constexpr const auto& value()const{return tuple_get<Ns...>(tuple).eval();}

//    template<Integer...Ns>
//          auto& value()     {return tuple_get<Ns...>(tuple).eval();}


//  ShapeFunctionsCollection(ShapeFunctionCoefficient<GeneralForm_,GeneralForms_...>&coeffs,
//                  ReferenceMaps2<GeneralForm_,GeneralForms_...>& maps,
//                  const GeneralForm_& form,const GeneralForms_&...forms):
//  coeffs_(coeffs),
//  maps_(maps)
//  ,
//  tuple_composite_(build_tuple_of_combination_functions<GeneralForm::FunctionSpace::Nuniquesubspaces>(form(),forms()...))
//  { }


// private:
//    ShapeFunctionCoefficient<GeneralForm_,GeneralForms_...> & coeffs_;
//    ReferenceMaps2<GeneralForm_,GeneralForms_...> & maps_;
//    TupleOfTupleShapeFunction tuple;
//    TupleOfTupleCompositeShapeFunction tuple_composite_;
//    TupleOfTupleCompositeShapeFunctionTensor tuple_tensors_;
// };


// template<typename Form,typename...Forms>
// constexpr auto shape_functions2(ShapeFunctionCoefficient<Form,Forms...>&coeffs, ReferenceMaps2<Form,Forms...>&maps,const Form& form,const Forms&...forms)
// {
//   //using Form=typename std::remove_const<typename std::remove_reference<ConstFormReference>::type>::type;
//  return ShapeFunctionsCollection<Form,Forms...>(coeffs,maps,form,forms...);  }
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //////// For explanation on how and why this works, check:
// //////// https://stackoverflow.com/questions/1005476/how-to-detect-whether-there-is-a-specific-member-variable-in-class
// //////// We check wheter T has a variable named is_test_or_trial
// //////// Such a variable is used only by Test or Trial and respectively set to 0 or 1
// //////// If in a linear or bilinear form, other quantities which are neither tests nor trials appear, then we can say so
// //////// because they do not have such is_test_or_trial
// //////// Therefore we can control how many is_test_or_trial appear and their respective values
// //////// In a 0-form, no is_test_or_trial appear
// //////// In a 1-form, is_test_or_trial=0 (test)
// //////// In a 2-form, is_test_or_trial=0 and 1 (test and trial)
// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// // template <typename T, typename = Integer>
// // struct ITestOrTrial : std::false_type { };

// // template <typename T>
// // struct ITestOrTrial <T,  decltype((void) T::is_test_or_trial, static_cast<decltype(T::is_test_or_trial)>(T::is_test_or_trial) )> : std::true_type { };




}


#endif