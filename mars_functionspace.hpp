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




template<typename...Ts>
class SingleSpaceDofsDofMap;

template<typename Elem,typename...BaseFunctionSpaces>
class SingleSpaceDofsDofMap<Elem,BaseFunctionSpaces...>
{
public:
 static constexpr Integer Nsubspaces=sizeof...(BaseFunctionSpaces);
 static constexpr auto NsubspacesArray=Array<Integer,1>{Nsubspaces};
 using ElemDofMap=std::tuple<Array<Integer, DofsPerElemNums<Elem,BaseFunctionSpaces>::value>...>;
 using DofMap=std::tuple<std::shared_ptr<std::vector<Array<Integer, DofsPerElemNums<Elem,BaseFunctionSpaces>::value>>>...>;
 using SpacesDofsArray=Array<std::shared_ptr<std::vector<Integer>>, Nsubspaces>;
 static constexpr auto NLocalDofsArray=Array<Integer,Nsubspaces>{DofsPerElemNums<Elem,BaseFunctionSpaces>::value...};
 static constexpr auto NLocalDofs=Sum(DofsPerElemNums<Elem,BaseFunctionSpaces>::value...);
 inline auto& space_dofs() const {return space_dofs_;};
 inline auto& space_dofs()       {return space_dofs_;};
 inline auto& dofmap() const {return dofmap_;};
 inline auto& dofmap()       {return dofmap_;};
 inline auto& n_dofs()       {return n_dofs_;}
 inline auto& n_dofs() const {return n_dofs_;}

 inline auto& level_n_dofs_array()       {return level_n_dofs_array_;}
 inline auto& level_n_dofs_array() const {return level_n_dofs_array_;}

 template<Integer N=0>
 inline std::enable_if_t<(N>=Nsubspaces),void>
 init_dofmap_aux(const Integer n_elements){};


 template<Integer N=0>
 inline std::enable_if_t<(N<Nsubspaces),void>
 init_dofmap_aux(const Integer n_elements)
 {
  // we reserve n_elemnts*dofs for single element for each space
  std::cout<<"begin N=="<<N<<std::endl;
  std::cout<<n_elements<<" "<< NLocalDofsArray[N]<<std::endl;
  tuple_get<N>(dofmap_)=std::make_shared<std::vector<Array<Integer,NLocalDofsArray[N]>>>();
  tuple_get<N>(dofmap_)->resize(n_elements);
  space_dofs_[N]=std::make_shared<std::vector<Integer>>();
  space_dofs_[N]->resize(n_elements*NLocalDofsArray[N]);
  std::cout<<"end N=="<<N<<" dofmap tuple size="<<tuple_get<N>(dofmap_)->size()<<" spacedofs size="<< space_dofs_[N]->size()<<std::endl;
  init_dofmap_aux<N+1>(n_elements);};

 inline void  init_dofmap(const Integer n_elements){init_dofmap_aux(n_elements);};









 // template<Integer N=0>
 // inline std::enable_if_t<(N>=Nsubspaces),void>
 // init_dofmap_aux2(const Integer n_elements){};


 // template<Integer N=0>
 // inline std::enable_if_t<(N<Nsubspaces),void>
 // init_dofmap_aux2(const Integer n_elements)
 // {
 //  // we reserve n_elemnts*dofs for single element for each space
 //  std::cout<<"begin N=="<<N<<std::endl;
 //  std::cout<<n_elements<<" "<< NLocalDofsArray[N]<<std::endl;
 //  tuple_get<N>(dofmap_)=std::make_shared<std::vector<Array<Integer,NLocalDofsArray[N]>>>();
 //  tuple_get<N>(dofmap_)->resize(n_elements);
 //  space_dofs_[N]=std::make_shared<std::vector<Integer>>();
 //  space_dofs_[N]->resize(n_elements*NLocalDofsArray[N]);
 //  std::cout<<"end N=="<<N<<" dofmap tuple size="<<tuple_get<N>(dofmap_)->size()<<" spacedofs size="<< space_dofs_[N]->size()<<std::endl;
 //  init_dofmap_aux<N+1>(n_elements);};

 // inline void  init_dofmap2(const Integer n_elements){init_dofmap_aux2(n_elements);};

 private:
  DofMap dofmap_;
  SpacesDofsArray space_dofs_;
  Array<Integer,Nsubspaces> n_dofs_;
  Array<std::vector<Integer>,Nsubspaces> level_n_dofs_array_;
};


template<Integer Nsubspaces, typename...Args>
class DofsDofMap //<Nsubspaces,Args...>
{
public:
   using ElemDofMap=TupleCatType<typename Args::DofsDM::ElemDofMap...>;
   using DofMap=TupleCatType<typename Args::DofsDM::DofMap...>;
   using SpacesDofsArray=Array<std::shared_ptr<std::vector<Integer>>, Nsubspaces>;
   using NdofsArray=Array<Integer, Nsubspaces>;
   using CumulativeNdofsArray=Array<Integer, Nsubspaces+1>;
   static constexpr auto NLocalDofs=Sum(Args::DofsDM::NLocalDofs...);
   // static constexpr Integer Nsubspaces=Sum(Args::DofsDM::Nsubspaces...);
   static constexpr auto NsubspacesArray=concat(Args::DofsDM::NsubspacesArray...);
   static constexpr auto NLocalDofsArray=concat(Args::DofsDM::NLocalDofsArray...);

   template<typename...DofsDofMaps>
   DofsDofMap(const DofsDofMaps&...dms):
   dofmap_(std::tuple_cat(dms.dofmap()...)),
   space_dofs_(concat(dms.space_dofs()...)),
   n_dofs_array_(concat(dms.n_dofs()...)),
   level_n_dofs_array_(concat(dms.level_n_dofs_array()...)),
   cumulative_dofs_array_(cumulative_array(n_dofs_array_,NsubspacesArray)),
   level_cumulative_dofs_array_(cumulative_array(level_n_dofs_array_,NsubspacesArray))
   {}

   inline auto& space_dofs() const {return space_dofs_;};
   inline auto& space_dofs()       {return space_dofs_;};
   inline auto& dofmap() const {return dofmap_;};
   inline auto& dofmap()       {return dofmap_;};
   inline auto& n_dofs()       {return n_dofs_array_;}
   inline auto& n_dofs() const {return n_dofs_array_;}
   inline auto& level_n_dofs_array()       {return level_n_dofs_array_;}
   inline auto& level_n_dofs_array() const {return level_n_dofs_array_;}
   inline auto& cumulative_n_dofs()       {return cumulative_dofs_array_;}
   inline auto& cumulative_n_dofs() const {return cumulative_dofs_array_;}
   inline auto& level_cumulative_dofs_array()       {return level_cumulative_dofs_array_;}
   inline auto& level_cumulative_dofs_array() const {return level_cumulative_dofs_array_;}

   template<Integer N=0>
   inline std::enable_if_t<(N>=Nsubspaces),void>
   init_dofmap_aux(const Integer n_elements){};


   template<Integer N=0>
   inline std::enable_if_t<(N<Nsubspaces),void>
   init_dofmap_aux(const Integer n_elements)
   {
    // we reserve n_elemnts*dofs for single element for each space
    // tuple_get<N>(dofmap_)->reserve(n_elements*NLocalDofsArray[N]);
    tuple_get<N>(dofmap_)->resize(n_elements);

    init_dofmap_aux<N+1>(n_elements);};

   inline void  init_dofmap(const Integer n_elements){init_dofmap_aux(n_elements);};


   template<Integer N>
   void dofmap_get(GetType<ElemDofMap,N>& elemdm, const Integer elem_id)const
   {
    const auto& dm=tuple_get<N>(dofmap_)->operator[](elem_id);
    const auto& offset=cumulative_dofs_array_[N];
    

    for(std::size_t i=0;i<elemdm.size();i++)
        elemdm[i]=dm[i]+offset;
    }

   template<Integer N>
   auto dofmap_size() const
   {return tuple_get<N>(dofmap_)->size();}


   auto& space_dofs_get(const Integer N,const Integer id)
   {return space_dofs_[N]->operator[](id);}

   auto& space_dofs_get(const Integer N,const Integer id)const
   {return space_dofs_[N]->operator[](id);}

   auto space_dofs_size(const Integer N) const
   {return space_dofs_[N]->size();}


 private:
   DofMap dofmap_;
   SpacesDofsArray space_dofs_;
   ElemDofMap elemdofmap_;
   NdofsArray n_dofs_array_;
   Array<std::vector<Integer>, Nsubspaces> level_n_dofs_array_;
   CumulativeNdofsArray cumulative_dofs_array_;
   Array<std::vector<Integer>, Nsubspaces+1> level_cumulative_dofs_array_;
};


template<Integer Nsubspaces, typename Arg1, typename Arg2>
class DofsDofMapFullAndAux //<Nsubspaces,Args...>
{
public:
   using ElemDofMap=TupleCatType<typename Arg1::DofsDM::ElemDofMap, typename Arg2::DofsDM::ElemDofMap>;
   using DofMap=TupleCatType<typename Arg1::DofsDM::DofMap,typename Arg2::DofsDM::DofMap>;
   using SpacesDofsArray=Array<std::shared_ptr<std::vector<Integer>>, Nsubspaces>;
   using NdofsArray=Array<Integer, Nsubspaces>;
   using CumulativeNdofsArray=Array<Integer, Nsubspaces+1>;
   static constexpr auto NLocalDofs=Sum(Arg1::DofsDM::NLocalDofs,Arg2::DofsDM::NLocalDofs);
   static constexpr auto NsubspacesArray=concat(Arg1::DofsDM::NsubspacesArray,Arg2::DofsDM::NsubspacesArray);
   static constexpr auto NLocalDofsArray=concat(Arg1::DofsDM::NLocalDofsArray,Arg2::DofsDM::NLocalDofsArray);

   // using ElemDofMap=TupleCatType<typename Args::DofsDM::ElemDofMap...>;
   // using DofMap=TupleCatType<typename Args::DofsDM::DofMap...>;
   // using SpacesDofsArray=Array<std::shared_ptr<std::vector<Integer>>, Nsubspaces>;
   // using NdofsArray=Array<Integer, Nsubspaces>;
   // using CumulativeNdofsArray=Array<Integer, Nsubspaces+1>;
   // static constexpr auto NLocalDofs=Sum(Args::DofsDM::NLocalDofs...);
   // static constexpr auto NsubspacesArray=concat(Args::DofsDM::NsubspacesArray...);

   template<Integer N1, typename...Ts1, Integer N2, typename...Ts2>
   DofsDofMapFullAndAux(const DofsDofMap<N1,Ts1...>& dm1,const DofsDofMap<N2,Ts2...>& dm2,Integer M):
   dofmap_(std::tuple_cat(dm1.dofmap(),dm2.dofmap())),
   space_dofs_(concat(dm1.space_dofs(),dm2.space_dofs())),
   n_dofs_array_(concat(dm1.n_dofs(),dm2.n_dofs())),
   level_n_dofs_array_(concat(dm1.level_n_dofs_array(),dm2.level_n_dofs_array())),
   // cumulative_dofs_array_(cumulative_array_and_zero(dm1.n_dofs(),dm2.n_dofs()),Arg1::DofsDM::NsubspacesArray)
   cumulative_dofs_array_(cumulative_array_and_zero(dm1.cumulative_n_dofs(),dm2.n_dofs())),
   level_cumulative_dofs_array_(cumulative_array_and_zero(dm1.level_cumulative_dofs_array(),dm2.level_n_dofs_array()))
   {}


   inline auto& space_dofs() const {return space_dofs_;};
   inline auto& space_dofs()       {return space_dofs_;};
   inline auto& dofmap() const {return dofmap_;};
   inline auto& dofmap()       {return dofmap_;};
   inline auto& n_dofs()       {return n_dofs_array_;}
   inline auto& n_dofs() const {return n_dofs_array_;}
   inline auto& level_n_dofs_array()       {return level_n_dofs_array_;}
   inline auto& level_n_dofs_array() const {return level_n_dofs_array_;}
   inline auto& cumulative_n_dofs()       {return cumulative_dofs_array_;}
   inline auto& cumulative_n_dofs() const {return cumulative_dofs_array_;}

   template<Integer N=0>
   inline std::enable_if_t<(N>=Nsubspaces),void>
   init_dofmap_aux(const Integer n_elements){};


   template<Integer N=0>
   inline std::enable_if_t<(N<Nsubspaces),void>
   init_dofmap_aux(const Integer n_elements)
   {
    // we reserve n_elemnts*dofs for single element for each space
    std::cout<<"begin N=="<<N<<std::endl;
    // tuple_get<N>(dofmap_)->reserve(n_elements*NLocalDofsArray[N]);
    tuple_get<N>(dofmap_)->resize(n_elements);
    std::cout<<"end N=="<<N<<std::endl;
    init_dofmap_aux<N+1>(n_elements);};

   inline void  init_dofmap(const Integer n_elements){init_dofmap_aux(n_elements);};


   template<Integer N>
   void dofmap_get(GetType<ElemDofMap,N>& elemdm, const Integer elem_id)const
   {
    const auto& dm=tuple_get<N>(dofmap_)->operator[](elem_id);
    const auto& offset=cumulative_dofs_array_[N];
    

    for(std::size_t i=0;i<elemdm.size();i++)
        elemdm[i]=dm[i]+offset;
    }


   template<Integer N>
   auto dofmap_size() const
   {return tuple_get<N>(dofmap_)->size();}


   auto space_dofs_get(const Integer N,const Integer id)
   {return space_dofs_[N]->operator[](id)+cumulative_dofs_array_[N];}

   auto space_dofs_get(const Integer N,const Integer id)const
   {return space_dofs_[N]->operator[](id)+cumulative_dofs_array_[N];}

   auto space_dofs_size(const Integer N) const
   {return space_dofs_[N]->size();}


 private:
   DofMap dofmap_;
   SpacesDofsArray space_dofs_;
   ElemDofMap elemdofmap_;
   NdofsArray n_dofs_array_;
   Array<std::vector<Integer>, Nsubspaces> level_n_dofs_array_;
   CumulativeNdofsArray cumulative_dofs_array_;
   Array<std::vector<Integer>, Nsubspaces+1> level_cumulative_dofs_array_;
};



// template<Integer N,Integer Nsubspaces,typename...Args>
// auto& get(DofsDofMap<Nsubspaces,Args...>& dm,Integer elem_id)
// {
//  return dm;
// }

template<typename MeshT>
class ElemEntityCollection;

template<typename MeshT>
class MeshAndEntity;

template<typename MeshT>
class ConnectivitySimpliacialMapCollection;

template<typename MeshT_, typename BaseFunctionSpace,typename...BaseFunctionSpaces>
class FunctionSpace 
{
public:
      
      using MeshT=MeshT_;
      using Elem= typename MeshT::Elem;
      using DofsDM=SingleSpaceDofsDofMap<Elem,BaseFunctionSpace,BaseFunctionSpaces...>;
      using TupleOfSpaces=std::tuple<ElemFunctionSpace<Elem,BaseFunctionSpace>,ElemFunctionSpace<Elem,BaseFunctionSpaces>... >;
      static constexpr Integer Nsubspaces=1+sizeof...(BaseFunctionSpaces);

      static constexpr Integer Nelem_dofs=DofsPerElemNums<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;
      static constexpr Array<Integer,Nsubspaces> Nelem_dofs_array=
      concat(Array<Integer,1>{DofsPerElemNums<Elem,BaseFunctionSpace>::value},Array<Integer,1>{DofsPerElemNums<Elem,BaseFunctionSpaces>::value}...);
      using DofMapType=std::vector<std::array<Integer, Nelem_dofs>>;
      using DofMapType2=std::tuple<std::vector<Array<Integer, DofsPerElemNums<Elem,BaseFunctionSpace>::value>>,
                                   std::vector<Array<Integer, DofsPerElemNums<Elem,BaseFunctionSpaces>::value>>
                                  ...>;

      using DofMapType3=std::tuple<std::shared_ptr<std::vector<Array<Integer, DofsPerElemNums<Elem,BaseFunctionSpace>::value>>>,
                                   std::shared_ptr<std::vector<Array<Integer, DofsPerElemNums<Elem,BaseFunctionSpaces>::value>>>...>;


      static constexpr auto faces_dofs_array=std::tuple_cat(std::make_tuple(TraceDofs<ElemFunctionSpace<Elem,BaseFunctionSpace>>::dofs()),
                                                             std::make_tuple(TraceDofs<ElemFunctionSpace<Elem,BaseFunctionSpaces>>::dofs())...);
      static constexpr Array<std::size_t,Nsubspaces> Nfaces_dofs_array{TraceDofs<ElemFunctionSpace<Elem,BaseFunctionSpace>>::dofs()[0].size(),TraceDofs<ElemFunctionSpace<Elem,BaseFunctionSpaces>>::dofs()[0].size()...};
      using OffSetType=Array<std::vector<Integer>, Nsubspaces>;
      using SpacesDofsArrayType=Array<std::vector<std::vector<Integer>>, Nsubspaces>;
      using SpacesDofsArrayType3=Array<std::shared_ptr<std::vector<Integer>>, Nsubspaces>;

      using SpacesInfosArrayType=Array<Array<Integer,4>,Nsubspaces>;
      using ArrayNdofs=Array<Integer,Nsubspaces>;
      using ElemsTupleType=std::tuple<Elem>;
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



      // inline const Integer& n_dofs()const{return n_dofs_;};

      inline const Integer& n_dofs()const{return dofsdm_.n_dofs()[0];};


      inline const auto& level_n_dofs_array()const{return dofsdm_.level_n_dofs_array();}

      inline const Integer& level_n_dofs_array(const Integer level)const
      {if(level==-1)
        return dofsdm_.level_n_dofs_array()[0][dofsdm_.level_n_dofs_array()[0].size()-1];
       else
        return dofsdm_.level_n_dofs_array()[0][level];};



      inline const Integer n_dofs(const Integer& space_id)const
                                 {Integer n=0;
                                  for(std::size_t i=0;i<space_infos_[space_id][3];i++)
                                       n+=space_dofs_[space_id][i].size(); 
                                     return n;
                                   };



      inline const Integer n_dofs(const Integer& space_id,const Integer& component_id)const
                                 {return space_dofs_[space_id][component_id].size(); };

      inline const DofMapType2& dofmap2()const{return dofmap2_;};

      inline const auto& dofmap3()const{return dofsdm_.dofmap();};

      inline const auto& space_dofs3()const{return dofsdm_.space_dofs();};

      inline const auto& dofsdofmap()const{return dofsdm_;};


      // auto& array_ndofs(){return array_ndofs_;}
      // auto& array_ndofs()const{return array_ndofs_;}

      // auto& level_array_ndofs(){return level_array_ndofs_;}
      // auto& level_array_ndofs()const{return level_array_ndofs_;}

      inline const Array<std::vector<Integer>, Nsubspaces>& offset() const {return offset_;};

      inline void  offset(const std::array<std::vector<Integer>, Nsubspaces> &os)const {os=offset_;};


      inline const std::vector<Integer>& offset(const Integer& space_id)const{return offset_[space_id];};

      inline void offset(Integer space_id, const std::vector<Integer>& space_os)const {space_os=offset_[space_id];};


      inline const auto& space_dofs() const {return space_dofs_;};




      inline const std::vector<Integer>& space_dofs(const Integer& space_id,const Integer& component_id) const
                                         {return space_dofs_[space_id][component_id];};


      inline void space_dofs(const Integer& space_id, const Integer& component_id,std::vector<Integer>& spacedofs)const
                            {spacedofs.resize(n_dofs(space_id,component_id));
                             spacedofs=space_dofs_[space_id][component_id];};

      inline const SpacesInfosArrayType& space_info()const{return space_infos_;};

      inline auto mesh_ptr()const {return mesh_ptr_;};

      inline auto& bisection()const {return bisection_;};
       
      // FunctionSpace(const MeshT& mesh,const Integer dof_count_start=0):
      // mesh_ptr_(std::make_shared< MeshT >(mesh))
      // {
      // function_space_info<0,Nsubspaces,BaseFunctionSpace,BaseFunctionSpaces...>(space_infos_);
      // // dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap2_,offset_,n_dofs_,space_infos_,space_dofs_,array_ndofs_);     
      // // dofmap_fespace3<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap3_,offset_,n_dofs_,space_infos_,space_dofs3_,array_ndofs_);     
      // // dofmap_fespace3<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofsdm_.dofmap(),offset_,n_dofs_,space_infos_,dofsdm_.space_dofs(),array_ndofs_);     
      // dofmap_fespace3<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofsdm_,offset_,array_ndofs_);     

      // // for(std::size_t i=0; i<Nsubspaces ;i++)
      // // dofsdm_.n_dofs()[i]=n_dofs_;
      // }

      // FunctionSpace( MeshAndEntity<MeshT>& mesh_and_entity):
      // mesh_ptr_(mesh_and_entity.mesh_ptr())
      // {
      // function_space_info<0,Nsubspaces,BaseFunctionSpace,BaseFunctionSpaces...>(space_infos_);
      // // dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap2_,offset_,n_dofs_,space_infos_,space_dofs_,array_ndofs_);     
      // // dofmap_fespace3<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap3_,offset_,n_dofs_,space_infos_,space_dofs3_,array_ndofs_);     
      // // dofmap_fespace3<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofsdm_.dofmap(),offset_,n_dofs_,space_infos_,dofsdm_.space_dofs(),array_ndofs_);     
      // dofmap_fespace4<BaseFunctionSpace,BaseFunctionSpaces...>(mesh_and_entity,dofsdm_,offset_,array_ndofs_);     

      // // for(std::size_t i=0; i<Nsubspaces ;i++)
      // // dofsdm_.n_dofs()[i]=n_dofs_;
      // }


      FunctionSpace( ConnectivitySimpliacialMapCollection<MeshT>& entities):
      mesh_ptr_(std::make_shared<MeshT>(entities.mesh())),
      bisection_(entities.bisection())

      {
      function_space_info<0,Nsubspaces,BaseFunctionSpace,BaseFunctionSpaces...>(space_infos_);
      // dofmap_fespace<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap2_,offset_,n_dofs_,space_infos_,space_dofs_,array_ndofs_);     
      // dofmap_fespace3<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofmap3_,offset_,n_dofs_,space_infos_,space_dofs3_,array_ndofs_);     
      // dofmap_fespace3<BaseFunctionSpace,BaseFunctionSpaces...>(mesh,dofsdm_.dofmap(),offset_,n_dofs_,space_infos_,dofsdm_.space_dofs(),array_ndofs_);     
      
      dofmap_fespace5<BaseFunctionSpace,BaseFunctionSpaces...>(entities,dofsdm_,offset_,level_array_ndofs_);     

      // for(std::size_t i=0; i<Nsubspaces ;i++)
      // dofsdm_.n_dofs()[i]=n_dofs_;
      }


private:
   
      std::shared_ptr< MeshT > mesh_ptr_;
      Integer n_dofs_;
      DofMapType dofmap_;
      DofMapType2 dofmap2_;
      DofMapType3 dofmap3_;
      SpacesDofsArrayType3 space_dofs3_;
      OffSetType offset_;
      SpacesDofsArrayType space_dofs_;
      SpacesInfosArrayType space_infos_;
      ArrayNdofs array_ndofs_;
      Array<std::vector<Integer>,Nsubspaces> level_array_ndofs_;
      DofsDM dofsdm_;
      Bisection<MeshT>& bisection_;
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
  using MeshT=typename Arg::MeshT;
  using Elem=typename Arg::Elem;
  // using MeshT=Mesh<Elem::Dim,Elem::ManifoldDim>;
  using TupleOfSpaces=TupleCatType<typename Arg::TupleOfSpaces,typename Args::TupleOfSpaces...>;
  using DofMapType=std::tuple<typename Arg::DofMapType,typename Args::DofMapType...>;
  using DofMapType2=TupleCatType<typename Arg::DofMapType2,typename Args::DofMapType2...>;
  using DofMapType3=TupleCatType<typename Arg::DofMapType3,typename Args::DofMapType3...>;

  using UniqueElementFunctionSpacesTupleType=
  RemoveTupleDuplicates<TupleCatType<typename Arg::UniqueElementFunctionSpacesTupleType,typename Args::UniqueElementFunctionSpacesTupleType...>>;
  using SpacesToUniqueFEFamily=SpacesToUniqueFEFamilies2<UniqueElementFunctionSpacesTupleType>;
  using FromElementFunctionSpacesToUniqueNumbersTupleType=
        TupleAllToUniqueMap<TupleOfSpaces,UniqueElementFunctionSpacesTupleType>;
  using FromElementFunctionSpacesToFirstSpaceTupleType=
        FromElementFunctionSpacesToFirstSpaceTupleType<FromElementFunctionSpacesToUniqueNumbersTupleType>;

  static constexpr Integer Nuniquesubspaces=TupleTypeSize<UniqueElementFunctionSpacesTupleType>::value;
  


  inline const Integer& n_dofs()const{return n_dofs_;}; 

  inline const auto& level_n_dofs_array(){return dofsdm_.level_n_dofs_array();}

      inline const Integer& level_n_dofs_array(const Integer level)const
      {std::cout<<"level_n_dofs_array, level=="<<level<<std::endl;
      std::cout<<"dofsdm_.level_n_dofs_array().size()-1=="<<dofsdm_.level_n_dofs_array().size()-1<<std::endl;
        if(level==-1)
        return dofsdm_.level_n_dofs_array()[0][dofsdm_.level_n_dofs_array()[0].size()-1];
       else
        return dofsdm_.level_n_dofs_array()[0][level];};


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
        using type = TupleCatType<OtherArg,OtherArgs...>;
      };


template<typename OtherArg>
      struct tuple_type2<OtherArg>
      {
       using type = OtherArg;
     };



template<Integer N=0,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< (N==0&&0==sizeof...(OtherArgs)),typename tuple_type2<OtherArg,OtherArgs...>::type >::type  
     cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     { 
     return otherarg;//add_costant(otherarg,0);
   }

template<Integer N=0,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< (0<N&&0==sizeof...(OtherArgs)),typename tuple_type2<OtherArg,OtherArgs...>::type >::type  
     cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     { 
      // std::cout<<"addin="<<tot_n_dofs<N-1>()<<std::endl;
     return // otherarg;//
         add_costant(otherarg,tot_n_dofs<N-1>());
     }

template<Integer N=0,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< 0<N && 0<sizeof...(OtherArgs),typename tuple_type2<OtherArg,OtherArgs...>::type >::type  
     cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     {
      // std::cout<<"addin="<<tot_n_dofs<N-1>()<<std::endl;
     return std::tuple_cat(add_costant(otherarg,tot_n_dofs<N-1>()) ,
       cumulative_tuple_make<N+1,OtherArgs...>(otherargs...));
   }

template<Integer N=0,typename OtherArg, typename...OtherArgs>
     typename std::enable_if< 0==N && 0<sizeof...(OtherArgs),typename tuple_type2<OtherArg,OtherArgs...>::type >::type  
     cumulative_tuple_make(const OtherArg& otherarg, const OtherArgs&...otherargs )
     {
     return std::tuple_cat( otherarg ,
       cumulative_tuple_make<N+1,OtherArgs...>(otherargs...));
   }


      // inline auto& array_ndofs(){return array_ndofs_;}
      // inline auto& array_ndofs()const{return array_ndofs_;}
     
      // inline auto& level_array_ndofs(){return level_array_ndofs_;}
      // inline auto& level_array_ndofs()const{return level_array_ndofs_;}


      inline const auto& space_dofs() const {return space_dofs_;};

     MixedSpace(const Arg& arg,const Args&...args):
     mesh_ptr_(arg.mesh_ptr()),
     spaces_(std::make_tuple(std::make_shared<Arg>(arg),std::make_shared<Args>(args)...)),
     n_dofs_(tot_n_dofs<1+sizeof...(Args)-1>()),
    // dofmap_(cumulative_tuple_make<0,typename Arg::DofMapType,typename Args::DofMapType...>(arg.dofmap(),args.dofmap()...))
 

      // dofmap2_(cumulative_tuple_make<0,Arg,Args...>(0,arg,args...))

     // dofmap2_(cumulative_tuple_make<0,typename Arg::DofMapType2,typename Args::DofMapType2...>
     dofmap2_(cumulative_tuple_make(arg.dofmap2(),args.dofmap2()...))
     // ,
     // dofmap2_(std::tuple_cat(arg.dofmap2(),args.dofmap2()...))//cumulative_tuple_make<0,typename Arg::DofMapType2,typename Args::DofMapType2...>(arg.dofmap2(),args.dofmap2()...))
     ,
     dofmap3_(std::tuple_cat(arg.dofmap3(),args.dofmap3()...))
     // ,
     // array_ndofs_(concat(arg.array_ndofs(),args.array_ndofs()...))
     // ,
     // level_array_ndofs_(concat(arg.level_array_ndofs(),args.level_array_ndofs()...))
     ,
     space_dofs_(concat(arg.space_dofs(),args.space_dofs()...))
     ,
     space_dofs3_(concat(arg.space_dofs3(),args.space_dofs3()...))
     ,
     dofsdm_(arg.dofsdofmap(),args.dofsdofmap()...)
     ,
     bisection_(arg.bisection())
     {}


     inline auto mesh_ptr()const {return mesh_ptr_;};

     constexpr const auto& spaces(){return spaces_;}
     static constexpr Integer Nsubspaces=tot_subspaces<Arg,Args...>::value;
     static constexpr Integer Nelem_dofs=Sum(Arg::Nelem_dofs, Args::Nelem_dofs...);
     static constexpr Array<Integer,Nsubspaces> Nelem_dofs_array=concat(Arg::Nelem_dofs_array,Args::Nelem_dofs_array...);

     static constexpr auto faces_dofs_array=std::tuple_cat(Arg ::faces_dofs_array,
                                                            Args::faces_dofs_array...);
     static constexpr auto Nfaces_dofs_array=concat(Arg::Nfaces_dofs_array,Args::Nfaces_dofs_array...);
     using SpacesDofsArrayType=Array<std::vector<std::vector<Integer>>, Nsubspaces>;
     using SpacesDofsArrayType3=Array<std::shared_ptr<std::vector<Integer>>, Nsubspaces>;

     using DofsDM=DofsDofMap<Nsubspaces,Arg,Args...>;

      inline const auto& dofmap3()const{return dofmap3_;};

      inline const auto& space_dofs3()const{return space_dofs3_;};

      inline const auto& dofsdofmap()const{return dofsdm_;};

      inline auto& bisection()const {return bisection_;};
   private:
    std::shared_ptr< MeshT > mesh_ptr_;
    std::tuple<std::shared_ptr<Arg>,std::shared_ptr<Args>...> spaces_;
    Integer n_dofs_;
    DofMapType dofmap_;
    DofMapType2 dofmap2_;
    DofMapType3 dofmap3_;
    Array<Integer,Nsubspaces> array_ndofs_;
    Array<std::vector<Integer>,Nsubspaces> level_array_ndofs_;
    SpacesDofsArrayType space_dofs_;
    SpacesDofsArrayType3 space_dofs3_;
    DofsDM dofsdm_;
    Bisection<MeshT>& bisection_;
};

template<typename...Args>
MixedSpace<Args...> MixedFunctionSpace(const Args&...args){return MixedSpace<Args...>(args...);};













template<typename...Args>
class AuxMixedSpace; 

template<typename MeshT_,typename Arg,typename...Args>
class AuxMixedSpace<AuxFunctionSpace<MeshT_,Arg>,AuxFunctionSpace<MeshT_,Args>...>
{
public:
  using MeshT=MeshT_;
  using Elem=typename AuxFunctionSpace<MeshT,Arg>::Elem;
  using TupleOfSpaces=TupleCatType<typename AuxFunctionSpace<MeshT,Arg>::TupleOfSpaces,
  typename AuxFunctionSpace<MeshT,Args>::TupleOfSpaces...>;
  using DofMapType=std::tuple<typename AuxFunctionSpace<MeshT,Arg>::DofMapType,
  typename AuxFunctionSpace<MeshT,Args>::DofMapType...>;

  using DofMapType2=TupleCatType<typename AuxFunctionSpace<MeshT,Arg>:: DofMapType2,
                                 typename AuxFunctionSpace<MeshT,Args>::DofMapType2...>;
  using DofMapType3=TupleCatType<typename AuxFunctionSpace<MeshT,Arg>:: DofMapType3,
                                 typename AuxFunctionSpace<MeshT,Args>::DofMapType3...>;


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
  dofmap2_(std::tuple_cat(arg.dofmap2(),args.dofmap2()...))
  // ,
  // array_ndofs_(concat(arg.array_ndofs(),args.array_ndofs()...))
  // ,
  // level_array_ndofs_(concat(arg.level_array_ndofs(),args.level_array_ndofs()...))
  ,
  space_dofs_(concat(arg.space_dofs(),args.space_dofs()...))
  ,
  dofmap3_(std::tuple_cat(arg.dofmap3(),args.dofmap3()...))
  ,
  space_dofs3_(concat(arg.space_dofs3(),args.space_dofs3()...))
  ,
  dofsdm_(arg.dofsdofmap(),args.dofsdofmap()...)
  {}

  inline auto mesh_ptr()const {return mesh_ptr_;};

  // inline const DofMapType& dofmap()const{return dofmap_;};

  inline const DofMapType2& dofmap2()const{return dofmap2_;};
  
  // auto& array_ndofs(){return array_ndofs_;}
  // auto& array_ndofs()const{return array_ndofs_;}

  // auto& level_array_ndofs(){return level_array_ndofs_;}
  // auto& level_array_ndofs()const{return level_array_ndofs_;}
  //    template<Integer...Ns>
  // inline constexpr const auto& dofmap()const{return tuple_get<Ns...>(dofmap_);};
  constexpr const auto& spaces(){return spaces_;}

   inline const auto& space_dofs() const {return space_dofs_;};
   inline const auto& space_dofs3() const {return space_dofs3_;};


   inline const auto& dofsdofmap()const{return dofsdm_;};


  static constexpr Integer Nsubspaces=TupleTypeSize<TupleOfSpaces>::value;
  using SpacesDofsArrayType=Array<std::vector<std::vector<Integer>>, Nsubspaces>;
  using SpacesDofsArrayType3=Array<std::shared_ptr<std::vector<Integer>>, Nsubspaces>;
  using DofsDM=DofsDofMap<Nsubspaces,AuxFunctionSpace<MeshT,Arg>,AuxFunctionSpace<MeshT,Args>...>;

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
  DofMapType3 dofmap3_;
  Array<Integer,Nsubspaces> array_ndofs_;
  Array<std::vector<Integer>,Nsubspaces> level_array_ndofs_;
  SpacesDofsArrayType space_dofs_;
  SpacesDofsArrayType3 space_dofs3_;
  DofsDM dofsdm_;
};

template<typename MeshT,typename...Args>
auto AuxFunctionSpacesBuild(const AuxFunctionSpace<MeshT,Args>&...args){return AuxMixedSpace<AuxFunctionSpace<MeshT,Args>...>(args...);};





template<typename ...Spaces>
class FullSpace;

template<typename...Args>
class FullSpace<MixedSpace<Args...>>
{
public:
  using MeshT=typename MixedSpace<Args...>::MeshT;
  using FunctionSpace=MixedSpace<Args...>;
  using Elem=typename FunctionSpace::Elem;
  // using MeshT=Mesh<Elem::Dim,Elem::ManifoldDim>;
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
  using SpacesDofsArrayType=Array<std::vector<std::vector<Integer>>, Nsubspaces>;
  using DofsDM=DofsDofMap<Nsubspaces,MixedSpace<Args...>>;


  FullSpace(const FunctionSpace& W):
  spaces_ptr_(std::make_shared<FunctionSpace>(W)),
  // array_ndofs_(W.array_ndofs()),
  // level_array_ndofs_(W.level_array_ndofs()),
  space_dofs_(W.space_dofs()),
  dofsdm_(W.dofsdofmap()),
  bisection_(W.bisection())
  {}

  // auto& array_ndofs(){return array_ndofs_;}
  // auto& array_ndofs()const{return array_ndofs_;}
  // auto& level_array_ndofs(){return level_array_ndofs_;}
  // auto& level_array_ndofs()const{return level_array_ndofs_;}


  inline const auto& space_dofs() const {return space_dofs_;};

  inline auto mesh_ptr(){return spaces_ptr_->mesh_ptr();}

  inline auto mesh_ptr()const{return spaces_ptr_->mesh_ptr();}

  inline       auto  spaces_ptr()     {return spaces_ptr_;}
  
  inline const auto& spaces_ptr()const{return spaces_ptr_;}
  

  inline auto& dofmap2()const{return spaces_ptr_->dofmap2();};
  inline auto& dofmap3()const{return spaces_ptr_->dofmap3();};
  inline const auto& space_dofs3() const {return spaces_ptr_->space_dofs3();};

  inline const auto& dofsdofmap()const{return dofsdm_;};

  inline auto& bisection()const {return bisection_;};


private:
std::shared_ptr<FunctionSpace> spaces_ptr_;
Array<Integer,Nsubspaces> array_ndofs_;
Array<std::vector<Integer>,Nsubspaces> level_array_ndofs_;
SpacesDofsArrayType space_dofs_;
DofsDM dofsdm_;
Bisection<MeshT>& bisection_;
};


template<typename...Args,typename...AuxArgs>
class FullSpace<MixedSpace<Args...>,AuxMixedSpace<AuxArgs...>>
{
public:
  using MeshT=typename MixedSpace<Args...>::MeshT;
  using FunctionSpace=MixedSpace<Args...>;
  using AuxFunctionSpace=AuxMixedSpace<AuxArgs...>;
  using Elem=typename FunctionSpace::Elem;
  // using MeshT=Mesh<Elem::Dim,Elem::ManifoldDim>;
  using TupleOfSpaces=TupleCatType<typename FunctionSpace::TupleOfSpaces,
                                   typename AuxFunctionSpace::TupleOfSpaces>;
  using DofMapType=std::tuple<typename FunctionSpace::DofMapType,
                              typename AuxFunctionSpace::DofMapType>;

  using DofMapType2=std::tuple<typename FunctionSpace::DofMapType2,
                              typename AuxFunctionSpace::DofMapType2>;
  using DofMapType3=std::tuple<typename FunctionSpace::DofMapType3,
                              typename AuxFunctionSpace::DofMapType3>;

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
  using SpacesDofsArrayType=Array<std::vector<std::vector<Integer>>, Nsubspaces>;
  using SpacesDofsArrayType3=Array<std::shared_ptr<std::vector<Integer>>, Nsubspaces>;
  using DofsDM=DofsDofMapFullAndAux<Nsubspaces,MixedSpace<Args...>,AuxMixedSpace<AuxArgs...>>;

  FullSpace(const FunctionSpace& W1, const AuxFunctionSpace& W2):
  spaces_ptr_(std::make_shared<FunctionSpace>(W1)),
  aux_spaces_ptr_(std::make_shared<AuxFunctionSpace>(W2)),
  // array_ndofs_(concat(W1.array_ndofs(),W2.array_ndofs())),
  // level_array_ndofs_(concat(W1.level_array_ndofs(),W2.level_array_ndofs())),
  space_dofs_(concat(W1.space_dofs(),W2.space_dofs())),
  space_dofs3_(concat(W1.space_dofs3(),W2.space_dofs3())),
  dofsdm_(W1.dofsdofmap(),W2.dofsdofmap(),1),
  bisection_(W1.bisection())
  {}

  inline       auto  spaces_ptr()     {return spaces_ptr_;}
  inline const auto& spaces_ptr()const{return spaces_ptr_;}
  inline       auto  aux_spaces_ptr()     {return aux_spaces_ptr_;}
  inline const auto& aux_spaces_ptr()const{return aux_spaces_ptr_;}

  // inline auto& array_ndofs(){return array_ndofs_;}
  // inline auto& array_ndofs()const{return array_ndofs_;}
  // auto& level_array_ndofs(){return level_array_ndofs_;}
  // auto& level_array_ndofs()const{return level_array_ndofs_;}


  inline const auto& space_dofs() const {return space_dofs_;};

  inline auto mesh_ptr(){return spaces_ptr_->mesh_ptr();}
  inline auto mesh_ptr()const{return spaces_ptr_->mesh_ptr();}

  // inline auto& dofmap2()const{return spaces_ptr_->dofmap2();};
  // inline auto& dofmap3()const{return spaces_ptr_->dofmap3();};
  inline auto& space_dofs3() const {return space_dofs3_;};
  
  inline       auto& dofsdofmap()     {return dofsdm_;};
  inline const auto& dofsdofmap()const{return dofsdm_;};
  inline auto& bisection()const {return bisection_;};
private:
std::shared_ptr<FunctionSpace> spaces_ptr_;
std::shared_ptr<AuxFunctionSpace> aux_spaces_ptr_;
Array<Integer,Nsubspaces> array_ndofs_;
Array<std::vector<Integer>,Nsubspaces> level_array_ndofs_;
SpacesDofsArrayType space_dofs_;
SpacesDofsArrayType3 space_dofs3_;
DofsDM dofsdm_;
Bisection<MeshT>& bisection_;
};

template<typename...Args>
auto FullSpaceBuild(const MixedSpace<Args...>& W)
    {return FullSpace<MixedSpace<Args...>>(W);}

template<typename...Args,typename...AuxArgs>
auto FullSpaceBuild(const MixedSpace<Args...>& W1,const AuxMixedSpace<AuxArgs...>& W2)
    {return FullSpace<MixedSpace<Args...>,AuxMixedSpace<AuxArgs...>>(W1,W2);}





// template<Integer N, typename...Args>
// void sub_solution(std::vector<Real>& sub_sol,const std::vector<Real>& sol,const FullSpace<Args...>& W )
// {
//  std::vector<Real> found_vec;

//  const auto& array_ndofs=W.array_ndofs();
//  auto n_dofs=array_ndofs[N];

//  if(sub_sol.size()!=n_dofs)
//   sub_sol.resize(n_dofs);

//  found_vec.resize(n_dofs);
 
//  Integer offset=0;

//  for(std::size_t i=0;i<N;i++) 
//  {
//   offset+=array_ndofs[i];
//  }

//  auto&dofmap=tuple_get<N>(W.dofmap2());


//  for(std::size_t el=0;el<dofmap.size();el++)
//  {
//   for(std::size_t i=0;i<dofmap[el].size();i++)
//    {
//     const auto index=dofmap[el][i];
//     sub_sol[index] = sol[index];
//    }

//  }

// }
}


#endif