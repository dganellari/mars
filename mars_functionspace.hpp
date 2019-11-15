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
      static constexpr auto faces_dofs_array=std::tuple_cat(std::make_tuple(TraceDofs<ElemFunctionSpace<Elem,BaseFunctionSpace>>::dofs()),
                                                             std::make_tuple(TraceDofs<ElemFunctionSpace<Elem,BaseFunctionSpaces>>::dofs())...);
      static constexpr Array<std::size_t,Nsubspaces> Nfaces_dofs_array{TraceDofs<ElemFunctionSpace<Elem,BaseFunctionSpace>>::dofs()[0].size(),TraceDofs<ElemFunctionSpace<Elem,BaseFunctionSpaces>>::dofs()[0].size()...};


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

      // inline const DofMapType& dofmap()const{return dofmap_;};

      inline const DofMapType2& dofmap2()const{return dofmap2_;};

      // inline void  dofmap(const DofMapType& dm)const{dm=dofmap_;};


      // inline const std::array<Integer, Nelem_dofs>& dofmap(const Integer& elem_id)const
      //                    {return dofmap_[elem_id];};

      // inline void  dofmap(const Integer& elem_id, const std::array<Integer, Nelem_dofs> & elem_dm)const
      //                    {elem_dm=dofmap_[elem_id];};


      // inline std::vector<Integer> 
      //              dofmap(const Integer& space_id,const Integer& elem_id)const{
      //                   const auto& os=offset_[space_id];
      //                   const auto& size=n_elem_dofs(space_id);
      //                   std::vector<Integer> output(size);
      //                   for(Integer nn=0;nn<size;nn++)
      //                        output[nn]=dofmap_[elem_id][nn+os[0]];
      //                   return output;};


      // inline std::vector<Integer> 
      //              dofmap(const Integer& space_id,const Integer& component_id,const Integer& elem_id)const{
      //                   const auto& comp=components(space_id);
      //                   const auto& os=offset_[space_id];
      //                   const auto& size= n_elem_dofs(space_id);
      //                   std::vector<Integer> output(size/comp);
      //                   space_infos_[space_id][3];
      //                   Integer mm=0;
      //                   for(Integer nn=component_id;nn<size;nn=nn+comp)
      //                        {output[mm]=dofmap_[elem_id][nn+os[0]];mm++;}
      //                   return output;};


      // inline void dofmap(const Integer& space_id,const Integer& component_id,const Integer& elem_id, const std::vector<Integer>& elem_space_dm)const{
      //                   const auto& comp=components(space_id);
      //                   const auto& os=offset_[space_id];
      //                   const auto& size= n_elem_dofs(space_id);
      //                   std::vector<Integer> output(size/comp);
      //                   space_infos_[space_id][3];
      //                   Integer mm=0;
      //                   elem_space_dm.resize(size/comp);
      //                   for(Integer nn=component_id;nn<size;nn=nn+comp)
      //                        {elem_space_dm[mm]=dofmap_[elem_id][nn+os[0]];mm++;}
      //                   };

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

  // inline const DofMapType& dofmap()const{return dofmap_;};
  
  // template<Integer...Ns>
  // inline constexpr const auto& dofmap()const{return tuple_get<Ns...>(dofmap_);};

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
      struct tuple_type2
      {
        // using rest = typename tuple_type<OtherArgs...>::type; 
        // using tuple_ens=std::tuple<OtherArg>;
        using type = TupleCatType<OtherArg,OtherArgs...>;
      };


template<typename OtherArg>
      struct tuple_type2<OtherArg>
      {
       using type = OtherArg;
     };

// template<typename OtherArg,typename...OtherArgs>
//       struct tuple_type
//       {
//         using rest = typename tuple_type<OtherArgs...>::type; 
//         using tuple_ens=std::tuple<OtherArg>;
//         using type = decltype( std::tuple_cat( std::declval< tuple_ens >(), std::declval< rest >() ) );
//       };


// template<typename OtherArg>
//       struct tuple_type<OtherArg>
//       {
//        using type = typename std::tuple<OtherArg>;
//      };

// template<Integer N,typename OtherArg, typename...OtherArgs>
//      typename std::enable_if< (N==0&&0==sizeof...(OtherArgs)),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
//      cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
//      { 
//      return std::tuple<OtherArg>(add_costant(otherarg,0));}

// template<Integer N,typename OtherArg, typename...OtherArgs>
//      typename std::enable_if< (0<N&&0==sizeof...(OtherArgs)),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
//      cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
//      { 
//      return std::tuple<OtherArg>(add_costant(otherarg,tot_n_dofs<N-1>()));}


// template<Integer N,typename OtherArg, typename...OtherArgs>
//      typename std::enable_if< 0<N && 0<sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
//      cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
//      {
//      return std::tuple_cat(std::tuple<OtherArg>( add_costant(otherarg,tot_n_dofs<N-1>()) ),
//        cumulative_tuple_make<N+1,OtherArgs...>(otherargs...));}

// template<Integer N,typename OtherArg, typename...OtherArgs>
//      typename std::enable_if< 0==N && 0<sizeof...(OtherArgs),typename tuple_type<OtherArg,OtherArgs...>::type >::type  
//      cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
//      {
//      return std::tuple_cat(std::tuple<OtherArg>( otherarg ),
//        cumulative_tuple_make<N+1,OtherArgs...>(otherargs...));}


template<Integer N,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< (N==0&&0==sizeof...(OtherArgs)),typename tuple_type2<OtherArg,OtherArgs...>::type >::type  
     cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     { 
     return otherarg;//add_costant(otherarg,0);
   }

template<Integer N,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< (0<N&&0==sizeof...(OtherArgs)),typename tuple_type2<OtherArg,OtherArgs...>::type >::type  
     cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     { 
      std::cout<<"addin="<<tot_n_dofs<N-1>()<<std::endl;
     return // otherarg;//
         add_costant(otherarg,tot_n_dofs<N-1>());
     }

template<Integer N,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< 0<N && 0<sizeof...(OtherArgs),typename tuple_type2<OtherArg,OtherArgs...>::type >::type  
     cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     {
      std::cout<<"addin="<<tot_n_dofs<N-1>()<<std::endl;
     return std::tuple_cat(add_costant(otherarg,tot_n_dofs<N-1>()) ,
       cumulative_tuple_make<N+1,OtherArgs...>(otherargs...));
   }

template<Integer N,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< 0==N && 0<sizeof...(OtherArgs),typename tuple_type2<OtherArg,OtherArgs...>::type >::type  
     cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     {
     return std::tuple_cat( otherarg ,
       cumulative_tuple_make<N+1,OtherArgs...>(otherargs...));
   }





     MixedSpace(const Arg& arg,const Args&...args):
     mesh_ptr_(arg.mesh_ptr()),
     spaces_(std::make_tuple(std::make_shared<Arg>(arg),std::make_shared<Args>(args)...)),
     n_dofs_(tot_n_dofs<1+sizeof...(Args)-1>()),
    // dofmap_(cumulative_tuple_make<0,typename Arg::DofMapType,typename Args::DofMapType...>(arg.dofmap(),args.dofmap()...))
 

      // dofmap2_(cumulative_tuple_make<0,Arg,Args...>(0,arg,args...))

     dofmap2_(cumulative_tuple_make<0,typename Arg::DofMapType2,typename Args::DofMapType2...>(arg.dofmap2(),args.dofmap2()...))
     // ,
     // dofmap2_(std::tuple_cat(arg.dofmap2(),args.dofmap2()...))//cumulative_tuple_make<0,typename Arg::DofMapType2,typename Args::DofMapType2...>(arg.dofmap2(),args.dofmap2()...))
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
  // dofmap_(arg.dofmap(),args.dofmap()...)
  // ,
  dofmap2_(arg.dofmap2(),args.dofmap2()...)
  {}

  inline auto mesh_ptr()const {return mesh_ptr_;};

  // inline const DofMapType& dofmap()const{return dofmap_;};

  inline const DofMapType2& dofmap2()const{return dofmap2_;};
  
  //    template<Integer...Ns>
  // inline constexpr const auto& dofmap()const{return tuple_get<Ns...>(dofmap_);};
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
  
  inline auto mesh_ptr(){return spaces_ptr_->mesh_ptr();}

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
  
  inline auto mesh_ptr(){return spaces_ptr_->mesh_ptr();}

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



}


#endif