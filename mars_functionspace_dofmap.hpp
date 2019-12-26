#ifndef MARS_FUNCTIONSPACE_DOFMAP_HPP
#define MARS_FUNCTIONSPACE_DOFMAP_HPP


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// 1) EntitiesOfFunctionSpace<Integer Dim, Integer ManifoldDim, Integer FEFamily, Integer Order>:                                   ////////                                        
//////    given the function space FEFamily of order Order,                                                                             ////////
//////    it builds a tuple of the entities related to its Dim_dofs                                                                     ////////
////// 2) EntitiesOfFunctionSpaceType is the type of EntitiesOfFunctionSpace                                                            ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mars_base.hpp"
#include "mars_elementfunctionspace.hpp"
#include "mars_elem_to_sub_elem.hpp"
#include "mars_tuple_utilities.hpp"
#include "mars_array.hpp"

#include "mars_tracker.hpp"
namespace mars{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// FunctionSpaceDofsPerElem:                                                                                     //////////////
/////// For the given entity of dimension entity[N] of FunctionSpace, returns the corresponding number of dofs        //////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename FunctionSpace,Integer N=FunctionSpace::entity.size()-1>
class 
FunctionSpaceDofsPerElem
{
     public: 
     using Elem = typename FunctionSpace::Elem; 
     static_assert(N>0," FunctionSpaceDofsPerElem N >0");
     static constexpr std::size_t ManifoldDim=FunctionSpace::ManifoldDim;
     static constexpr std::size_t n_components=FunctionSpace::NComponents;
     static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
     static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
     static constexpr std::size_t dofs_per_elem=n_components * 
                                                dofs_per_entity   * 
                                                ElemEntityCombinations<Elem,entity_dim>::value;
     
     static constexpr std::size_t value=FunctionSpaceDofsPerElem<FunctionSpace,N-1>::value+dofs_per_elem;
};

template<typename FunctionSpace>
class
FunctionSpaceDofsPerElem<FunctionSpace,0>
{
     public:  
     using Elem = typename FunctionSpace::Elem; 
     static constexpr std::size_t N=0;
     static constexpr std::size_t ManifoldDim=FunctionSpace::ManifoldDim;
     static constexpr std::size_t n_components=FunctionSpace::NComponents;
     static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
     static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
     static constexpr std::size_t value=n_components * 
                                        dofs_per_entity *
                                        ElemEntityCombinations<Elem,entity_dim>::value;
};





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// FunctionSpaceDofsPerSubEntityElem:                                                                       //////////////
/////// For the given entity of dimension entity[N] of FunctionSpace, returns the corresponding number of dofs   //////////////
/////// for a given subentity of the element                                                                     //////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename FunctionSpace,Integer SubElemDim,Integer N=FunctionSpace::entity.size()-1>
class 
FunctionSpaceDofsPerSubEntityElem
{
     public: 
     using Elem = typename FunctionSpace::Elem; 
     static_assert(N>0," FunctionSpaceDofsPerElem N >0");
     static constexpr std::size_t ManifoldDim=FunctionSpace::ManifoldDim;
     static constexpr std::size_t n_components=FunctionSpace::NComponents;
     static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
     static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
     static constexpr std::size_t dofs_per_elem=n_components * 
                                                dofs_per_entity   * 
                                                ElemEntityCombinations<typename ElemToSubElemHelper<Elem,SubElemDim>::type,entity_dim>::value;
     
     static constexpr std::size_t value=FunctionSpaceDofsPerSubEntityElem<FunctionSpace,SubElemDim,N-1>::value+dofs_per_elem;
};

template<typename FunctionSpace,Integer SubElemDim>
class
FunctionSpaceDofsPerSubEntityElem<FunctionSpace,SubElemDim,0>
{
     public:  
     using Elem = typename FunctionSpace::Elem; 
     static constexpr std::size_t N=0;
     static constexpr std::size_t ManifoldDim=FunctionSpace::ManifoldDim;
     static constexpr std::size_t n_components=FunctionSpace::NComponents;
     static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
     static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
     static constexpr std::size_t value=n_components * 
                                        dofs_per_entity *
                                        ElemEntityCombinations<typename ElemToSubElemHelper<Elem,SubElemDim>::type,entity_dim>::value;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// EntitiesOfFunctionSpaceType:                                                                                  //////////////
/////// For given FunctionSpaces on a given Elem, returns the entities corresponding to its dofs                      //////////////                                        
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Elem, Integer FEFamily, Integer Order,Integer N>
struct EntitiesOfFunctionSpaceType
{
      using rest = typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N-1>::type; 
      using ens = ElemEntity<Elem,ElementFunctionSpace<Elem,FEFamily,Order>::entity[N]>;
      using tuple_ens=std::tuple<ens>;
      using type = decltype( std::tuple_cat( std::declval< rest >(), std::declval< tuple_ens >() ) );
};


template<typename Elem,Integer FEFamily, Integer Order>
struct EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,0>
{
     using ens = ElemEntity<Elem,ElementFunctionSpace<Elem,FEFamily,Order>::entity[0]>;
     using type = typename std::tuple<ens>;
};


template<typename Elem,Integer FEFamily, 
         Integer Order, 
         Integer N=(ElementFunctionSpace<Elem,FEFamily,Order>::entity.size()-1),
         typename MeshT>
typename std::enable_if<N==0, 
                        typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type>::type 
EntitiesOfFunctionSpace(const MeshT& mesh, 
                        const std::vector< std::vector<Integer> > &node_2_element)
{   
      using ens=ElemEntity<Elem,ElementFunctionSpace<Elem,FEFamily,Order>::entity[N]>;
      return std::tuple<ens>(ens(mesh,node_2_element));
}


template<typename Elem,Integer FEFamily, 
         Integer Order, 
         Integer N=(ElementFunctionSpace<Elem,FEFamily,Order>::entity.size()-1),
         typename MeshT>
typename std::enable_if< 0<N, 
                         typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type >::type  
EntitiesOfFunctionSpace(const MeshT& mesh, 
                        const std::vector< std::vector<Integer> >&node_2_element )
{

      using ens=ElemEntity<Elem,ElementFunctionSpace<Elem,FEFamily,Order>::entity[N]>;
      return std::tuple_cat(EntitiesOfFunctionSpace<Elem,FEFamily,Order,N-1>(mesh,node_2_element),
                            std::tuple<ens>(ens(mesh,node_2_element)));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// tuple2array_dofs_offset:                                                                                  //////////////
/////// Transform the dofs_offset tuple to a more suitable array for the functionspace                            //////////////                                        
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<Integer N, typename T, typename A>
typename std::enable_if< 0==N, void>::type
tuple2array_dofs_offset(const T &tuple,A& array)
{
      auto tuple_N=std::get<N>(tuple);
      auto tuple_N_size=tuple_N.size();
      array[N].resize(tuple_N_size);
      for(Integer nn=0;nn<tuple_N_size;nn++)
          array[N][nn]=tuple_N[nn];
}


template<Integer N, typename T, typename A>
typename std::enable_if< 0<N, void>::type
tuple2array_dofs_offset(const T &tuple,A& array)
{
      auto tuple_N=std::get<N>(tuple);
      auto tuple_N_size=tuple_N.size();
      array[N].resize(tuple_N_size);
      for(Integer nn=0;nn<tuple_N_size;nn++)
          array[N][nn]=tuple_N[nn];
      tuple2array_dofs_offset<N-1>(tuple,array);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// initialize_flag_tuple_entities, FlagTuple:                                                                    //////////////
/////// We create a flag tuple (for each space) of tuples (for each entity of the dofs of the FunctionSpace)          ////////////// 
/////// So we can remember, for each space, the entities already visited and we can create the dofmap                 //////////////                                       
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Elem, Integer FEFamily, Integer Order, Integer M  , typename ...Args>
typename std::enable_if< M==sizeof ...(Args),void >::type
initialize_flag_tuple_entities(std::tuple< Args...> const &tuple, 
                               std::array<std::vector< std::array<Integer,2> >, 
                               ElementFunctionSpace<Elem,FEFamily,Order>::entity.size() > & entity_found)
    {static_assert(M>=0," the tuple must have non negative length ");};

template<typename Elem, Integer FEFamily, Integer Order, Integer M = 0, typename...Args>
typename std::enable_if< M<sizeof ...(Args),void >::type
initialize_flag_tuple_entities(std::tuple< Args...> const &tuple,
                               std::array<std::vector< std::array<Integer,2> >, 
                               ElementFunctionSpace<Elem,FEFamily,Order>::entity.size() > & entity_found)
{
     static_assert(M>=0," the tuple must have non negative length ");
     Integer entity_length=std::get<M>(tuple).size();
     entity_found[M].resize(entity_length,{-1,-1});
     initialize_flag_tuple_entities<Elem,FEFamily,Order,M+1,Args...>(tuple,entity_found);       
};

template<typename Elem,typename FunctionSpace,typename...FunctionSpaces>
struct FlagTupleType
{
     using FS=ElemFunctionSpace<Elem,FunctionSpace>;
     static constexpr Integer FEFamily=FS::FEFamily;
     static constexpr Integer Order=FS::Order; 
     static constexpr Integer entities_nums=FS::entity.size();
     using rest = typename FlagTupleType<Elem,FunctionSpaces...>::type;
     using ens  = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
     using type = decltype( std::tuple_cat( std::declval< rest >(), 
                                            std::declval< std::tuple<ens> >() ) );
};


template<typename Elem,typename FunctionSpace>
struct FlagTupleType<Elem,FunctionSpace>
{
 using FS=ElemFunctionSpace<Elem,FunctionSpace>;
 static constexpr Integer FEFamily=FS::FEFamily;
 static constexpr Integer Order=FS::Order; 
 static constexpr Integer entities_nums=FS::entity.size();
 using ens =  typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
 using type = typename std::tuple<ens>;
};


template<typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename...Args>
typename std::enable_if< 0<sizeof...(FunctionSpaces), 
                         typename FlagTupleType<Elem,FunctionSpace,FunctionSpaces...>::type>::type
FlagTuple(std::tuple<Args...> tuple)
{    
      using FS=ElemFunctionSpace<Elem,FunctionSpace>;
      constexpr Integer FEFamily=FS::FEFamily;
      constexpr Integer Order=FS::Order;
      static constexpr Integer M=sizeof...(FunctionSpaces);
      static constexpr Integer entities_nums=FS::entity.size();
      using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
      type ens;

      initialize_flag_tuple_entities<Elem,FEFamily,Order>(std::get<M>(tuple),ens);
      return std::tuple_cat(FlagTuple<Elem,FunctionSpaces...>(tuple),
                            std::tuple<type>(ens)
                           );
}

template<typename Elem,typename FunctionSpace, typename...FunctionSpaces, typename...Args>
 typename std::enable_if< 0==sizeof...(FunctionSpaces), 
                          typename FlagTupleType<Elem,FunctionSpace>::type>::type
FlagTuple(std::tuple<Args...> tuple)
{   
      using FS=ElemFunctionSpace<Elem,FunctionSpace>;
      constexpr Integer FEFamily=FS::FEFamily;
      constexpr Integer Order=FS::Order;
      static constexpr Integer entities_nums=FS::entity.size();
      using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
      type ens;
      initialize_flag_tuple_entities<Elem,FEFamily,Order>(std::get<0>(tuple),ens);
      return std::tuple<type>(ens);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// EntitiesOfFunctionSpaceTuple:                                                                             //////////////
/////// Tuple of EntitiesOfFunctionSpace for more FunctionSpaces on a given Elem                                  //////////////                                        
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Elem,typename FunctionSpace,typename...FunctionSpaces>
struct EntitiesOfFunctionSpaceTupleType
{
     using FS=ElemFunctionSpace<Elem,FunctionSpace>;
     static constexpr Integer FEFamily=FS::FEFamily;
     static constexpr Integer Order=FS::Order; 
     static constexpr Integer N=FS::entity.size()-1;
     using rest = typename EntitiesOfFunctionSpaceTupleType<Elem,FunctionSpaces...>::type;
     using ens  = typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type;
     using type = decltype( std::tuple_cat( std::declval< rest >(), 
                                            std::declval< std::tuple<ens> >() ) );
};


template<typename Elem,typename FunctionSpace>
struct EntitiesOfFunctionSpaceTupleType<Elem,FunctionSpace>
{
     using FS=ElemFunctionSpace<Elem,FunctionSpace>;
     static constexpr Integer FEFamily=FS::FEFamily;
     static constexpr Integer Order=FS::Order; 
     static constexpr Integer N=FS::entity.size()-1;
     using ens =  typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type ;
     using type = typename std::tuple<ens>;
};


template<typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename MeshT>
typename std::enable_if< 0<sizeof...(FunctionSpaces), 
                         typename EntitiesOfFunctionSpaceTupleType<Elem,FunctionSpace,FunctionSpaces...>::type>::type
EntitiesOfFunctionSpaceTuple(const MeshT& mesh, 
                             const std::vector< std::vector<Integer> > &node_2_element)
{    
      using FS=ElemFunctionSpace<Elem,FunctionSpace>;
      static constexpr Integer FEFamily=FS::FEFamily;
      static constexpr Integer Order=FS::Order;
      static constexpr Integer N=FS::entity.size()-1;
      using type = typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type;
      return std::tuple_cat(EntitiesOfFunctionSpaceTuple<Elem,FunctionSpaces...>(mesh,node_2_element),
                            std::tuple<type>(EntitiesOfFunctionSpace<Elem,FEFamily,Order>(mesh,node_2_element))
                           );
}

template<typename Elem,typename FunctionSpace, typename...FunctionSpaces, typename MeshT>
 typename std::enable_if< 0==sizeof...(FunctionSpaces), 
                          typename EntitiesOfFunctionSpaceTupleType<Elem,FunctionSpace>::type>::type
EntitiesOfFunctionSpaceTuple(const MeshT& mesh, 
                             const std::vector< std::vector<Integer> > &node_2_element)
{   
      using FS=ElemFunctionSpace<Elem,FunctionSpace>;
      static constexpr Integer FEFamily=FS::FEFamily;
      static constexpr Integer Order=FS::Order;
      static constexpr Integer N=FS::entity.size()-1;
      using type = typename EntitiesOfFunctionSpaceType<Elem,FEFamily,Order,N>::type;
      return std::tuple<type>(EntitiesOfFunctionSpace<Elem,FEFamily,Order>(mesh,node_2_element));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// ElementDofMap, ElementDofMap_LoopEntities:                                                                    //////////////
/////// We create the element dofmap by looping on the functionspaces and the corresponding entities                  //////////////                                       
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




template<Integer N,typename FunctionSpace,typename FlagArray,typename...Args1, long M,typename OS, typename SpaceDofs>
// template<Integer N,typename FunctionSpace,typename FlagArray,typename...Args1, long M,typename OS>
typename std::enable_if< -1<N,void>::type
ElementDofMap_LoopEntities(const std::tuple<Args1...>& entitiestuple,
                           FlagArray& flagtuples,
                           const Integer& elem_id, 
                           std::vector<Array<Integer, M>> & dofmap_vec,
                           Integer& global_dof_count,
                           Integer& loc_dof_count,
                           const OS &dofs_offset,
                           SpaceDofs& space_dofs,
                           Integer& dofs_count)
{

using Elem = typename FunctionSpace::Elem;
constexpr auto ManifoldDim=FunctionSpace::ManifoldDim;
constexpr auto continuity=FunctionSpace::Continuity;
constexpr auto n_components=FunctionSpace::NComponents;
constexpr auto entity_dim=FunctionSpace::entity[N];
constexpr auto dofs_per_entity=FunctionSpace::dofs_per_entity[N];
constexpr auto entity_points=entity_dim+1;
constexpr auto manifold_points=ManifoldDim+1;
constexpr auto combinations_nums=ElemEntityCombinations<Elem,entity_dim>::value;
const     auto& entity=std::get<N>(entitiestuple);
const     auto& elem2entity=entity.elem_2_entity(elem_id);
          auto& flag=flagtuples[N];


// move to the entity of dimension entity[N-1]
ElementDofMap_LoopEntities<N-1,FunctionSpace>
             (entitiestuple,flagtuples,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,dofs_count);
// ElementDofMap_LoopEntities<N-1,FunctionSpace>
//              (entitiestuple,flagtuples,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,dofs_count);


// std::cout<<"ElementDofMap_LoopEntities entity_dim = "<<entity_dim<<std::endl;

// loop on all the entities of a given entity_dim
for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
   { 
    const auto& entity_id=elem2entity[entity_iter];
    // std::cout<<"entity_iter = "<<entity_iter<<std::endl;
    // std::cout<<"entity_id = "<<entity_id<<std::endl;
    // std::cout<<"entity_id already visited ="<<flag[entity_id][0]<<std::endl;

    // if the entity has not been already visited, then create n_components new dofs
    if(flag[entity_id][0]==-1 || continuity==Discontinuous)
    {

       
      flag[entity_id][0]=elem_id;
      flag[entity_id][1]=entity_iter;
      
      for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
         for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
            {            
             dofmap_vec[elem_id][loc_dof_count]=global_dof_count;
             space_dofs[fs_dim].push_back(global_dof_count);

             // std::cout<<" global_dof_count (new)= "<<global_dof_count<<std::endl;
             loc_dof_count++;
             global_dof_count++;
             dofs_count++;
             }
    }
    else
    {
     // if the entity has been already visited, find the element (elem_id_tmp) in which was found 
     // and the corresponding entity_iter (iter_tmp)
     // we do not need any offset, because dofmap_vec contains  of the entity N, we move to the iter_tmp entity and find in dofmap_vec the already numbered dof
     const auto& elem_id_tmp=flag[entity_id][0];
     const auto& iter_tmp=flag[entity_id][1];
           // auto old_dof=dofmap_vec[elem_id_tmp][dofs_per_entity*iter_tmp*FunctionSpace::NComponents];
           auto old_dof=dofmap_vec[elem_id_tmp][dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents];
     




     // for(std::size_t i=0;i<dofmap_vec[elem_id_tmp].size();i++)
      // std::cout<<"dofmap_vec[elem_id_tmp] = "<<dofmap_vec[elem_id_tmp][i]<<std::endl;
     // std::cout<<"elem_id_tmp = "<<elem_id_tmp<<std::endl;
     

     
     // std::cout<<"dofs_offset[N] = "<<dofs_offset[N]<<std::endl;
     // std::cout<<"dofs_per_entity = "<<dofs_per_entity<<std::endl;
     // std::cout<<"iter_tmp = "<<iter_tmp<<std::endl;
     // std::cout<<"NComponents = "<<FunctionSpace::NComponents<<std::endl;

     // std::cout<<dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents<<std::endl;

     // std::cout<<"old_dof = "<<old_dof<<std::endl;
     
     
     for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
        for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
           {           
            dofmap_vec[elem_id][loc_dof_count]=old_dof;
            // std::cout<<" global_dof_count (old)= "<<old_dof<<std::endl;
            loc_dof_count++;
            old_dof++;



           }
    } 
  }
}

template<Integer N,typename FunctionSpace, typename FlagArray,typename...Args1, long M,typename OS, typename SpaceDofs>
// template<Integer N,typename FunctionSpace, typename FlagArray,typename...Args1, long M,typename OS>
typename std::enable_if< -1==N,void>::type
ElementDofMap_LoopEntities(const std::tuple<Args1...>& entitiestuple,
              FlagArray& flagtuples,
              const Integer& elem_id,
              std::vector<Array<Integer, M>> & dofmap_vec,
              Integer& global_dof_count,
              Integer &loc_dof_count,
              const OS &dofs_offset,
              SpaceDofs& space_dofs,
              Integer& dofs_count  )
{
  // std::cout<<"ElementDofMap_LoopEntities N="<<N<<std::endl;
}



template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename...Args1,typename...Args2,typename...Args3,typename OS, typename SpaceDofs,typename ArrayNdofs>
// template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename...Args1,typename...Args2,typename...Args3,typename OS,typename ArrayNdofs>
typename std::enable_if< 0<sizeof...(FunctionSpaces),void>::type
ElementDofMap(const std::tuple<Args1...>& entitiestuple,
                   std::tuple<Args2...>& flagtuples,
                   const Integer& elem_id,
                   std::tuple<Args3...>& dofmap_vec,
                   Integer& global_dof_count,
                   Integer& loc_dof_count,
                   OS &dofs_offset,
                   SpaceDofs& space_dofs,
                   ArrayNdofs& array_ndofs   )
{
 static constexpr Integer M=sizeof...(FunctionSpaces);
 const auto m1=std::get<M>(entitiestuple);
 auto &m2=std::get<M>(flagtuples);
 static constexpr Integer NN=std::tuple_size<decltype(m1)>::value-1;
 // std::cout<<"ElementDofMap="<<K<<std::endl;
 loc_dof_count=0;
 ElementDofMap_LoopEntities<NN,ElemFunctionSpace<Elem,FunctionSpace>>
                            (m1,m2,elem_id,tuple_get<K>(dofmap_vec),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),space_dofs[K],array_ndofs[K]);
  ElementDofMap<K+1,Elem,FunctionSpaces...>(entitiestuple,flagtuples,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,array_ndofs);

 
 // ElementDofMap_LoopEntities<NN,ElemFunctionSpace<Elem,FunctionSpace>>
 //                            (m1,m2,elem_id,tuple_get<K>(dofmap_vec),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),array_ndofs[K]);
 
 // ElementDofMap<K+1,Elem,FunctionSpaces...>(entitiestuple,flagtuples,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,array_ndofs);

};




template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename...Args1,typename...Args2,typename...Args3,typename OS, typename SpaceDofs,typename ArrayNdofs>
// template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename...Args1,typename...Args2,typename...Args3,typename OS,typename ArrayNdofs>
typename std::enable_if< 0==sizeof...(FunctionSpaces),void>::type
ElementDofMap(const std::tuple<Args1...>& entitiestuple,
                   std::tuple<Args2...>& flagtuples,
                   const Integer& elem_id,
                   std::tuple<Args3...>& dofmap_vec,
                   Integer& global_dof_count, 
                   Integer& loc_dof_count,
                   OS &dofs_offset,
                   SpaceDofs& space_dofs,
                   ArrayNdofs& array_ndofs )
{
 const auto m1=std::get<0>(entitiestuple);
 auto& m2=std::get<0>(flagtuples);
 static constexpr Integer NN=std::tuple_size<decltype(m1)>::value-1;
 // std::cout<<"ElementDofMap="<<K<<std::endl;
 loc_dof_count=0;

 ElementDofMap_LoopEntities<NN,ElemFunctionSpace<Elem,FunctionSpace>>
                           (m1,m2,elem_id,tuple_get<K>(dofmap_vec),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),space_dofs[K],array_ndofs[K]);
 // ElementDofMap_LoopEntities<NN,ElemFunctionSpace<Elem,FunctionSpace>>
 //                           (m1,m2,elem_id,tuple_get<K>(dofmap_vec),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),array_ndofs[K]);


};


































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// ElementDofMap, ElementDofMap_LoopEntities:                                                                    //////////////
/////// We create the element dofmap by looping on the functionspaces and the corresponding entities                  //////////////                                       
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<Integer EntityDim,typename ...T>
class DofsOrdering;
template<Integer EntityDim,typename ...T>
class DofsOrdered;

template<Integer N,typename FunctionSpace,typename FlagArray,typename Element,typename...Args1, long M,typename OS, typename SpaceDofs>
// template<Integer N,typename FunctionSpace,typename FlagArray,typename...Args1, long M,typename OS>
typename std::enable_if< -1<N,void>::type
ElementDofMap_LoopEntities3(const std::tuple<Args1...>& entitiestuple,
                           FlagArray& flagtuples,
                           Element& elem,
                           const Integer& elem_id, 
                           std::shared_ptr<std::vector<Array<Integer, M>>> & dofmap_vec,
                           Integer& global_dof_count,
                           Integer& loc_dof_count,
                           const OS &dofs_offset,
                           SpaceDofs& space_dofs,
                           Integer& dofs_count)
{

using Elem = typename FunctionSpace::Elem;
constexpr auto ManifoldDim=FunctionSpace::ManifoldDim;
constexpr auto continuity=FunctionSpace::Continuity;
constexpr auto n_components=FunctionSpace::NComponents;
constexpr auto entity_dim=FunctionSpace::entity[N];
constexpr auto dofs_per_entity=FunctionSpace::dofs_per_entity[N];
constexpr auto entity_points=entity_dim+1;
constexpr auto manifold_points=ManifoldDim+1;
constexpr auto combinations_nums=ElemEntityCombinations<Elem,entity_dim>::value;
const     auto& entity=std::get<N>(entitiestuple);
const     auto& elem2entity=entity.elem_2_entity(elem_id);
          auto& flag=flagtuples[N];

Integer cont;

// move to the entity of dimension entity[N-1]
ElementDofMap_LoopEntities3<N-1,FunctionSpace>
             (entitiestuple,flagtuples,elem,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,dofs_count);

const auto& nodes=elem.nodes;

// loop on all the entities of a given entity_dim
for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
   { 
    const auto& entity_id=elem2entity[entity_iter];

    if(DofsOrdering<entity_dim,FunctionSpace>::must_reorder)
    {
      auto ordered_entity_nodes=DofsOrdering<entity_dim,FunctionSpace>::value(nodes,entity_iter);
       // for(std::size_t i=0;i<ordered_entity_nodes.size();i++)
       //      {
       //        std::cout<<ordered_entity_nodes[i]<<" ";
       //      }

    cont=0;
    if(flag[entity_id][0]==-1 || continuity==Discontinuous)
    {

      // std::cout<< " elem id = "<< elem_id <<" NEW, iter ="<<entity_iter<<std::endl;

      auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,loc_dof_count);
      // std::cout<<"ordered_dofs"<<std::endl;
      // std::cout<<ordered_dofs<<std::endl;
       
      flag[entity_id][0]=elem_id;
      flag[entity_id][1]=entity_iter;
      // std::cout<<"elem_id = "<<elem_id<<std::endl;
      // std::cout<<"entity_iter = "<<entity_iter<<std::endl;
      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
      // std::cout<<std::endl;
      for(Integer m=0;m<ordered_dofs.size();m++)
      {
        (*dofmap_vec)[elem_id][ordered_dofs[cont]]=global_dof_count;
        // std::cout<<(*dofmap_vec)[elem_id][ordered_dofs[cont]]<<std::endl;
        cont++;
        loc_dof_count++;
        // std::cout<<global_dof_count<<std::endl;
        global_dof_count++;
        dofs_count++;
      }
      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
      // std::cout<<std::endl;
     }
     else
     {
      // std::cout<< " elem id = "<< elem_id<<" NEW, iter ="<<entity_iter<<std::endl;

      
      
      const auto& elem_id_tmp=flag[entity_id][0];
      const auto& iter_tmp=flag[entity_id][1];
            auto& elem_dm=(*dofmap_vec)[elem_id_tmp];
            auto ex_loc_dof_count=(dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents);
            auto global_old_dof=elem_dm[ex_loc_dof_count];
     

      for (int i = 0; i < n_components*dofs_per_entity; ++i)
      {
        const auto m= i +ex_loc_dof_count;
        if(elem_dm[m]<global_old_dof)
        {
          global_old_dof=elem_dm[m];
        }
      }

      // std::cout<<"old offset = = "<<dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents<<std::endl;
      // std::cout<<"dofs_offset= = "<<dofs_offset[N]<<std::endl;
      // std::cout<<"dofs_per_entity= = "<<dofs_per_entity<<std::endl;
      // std::cout<<"iter_tmp= = "<<iter_tmp<<std::endl;
      // std::cout<<"loc_dof_count = = "<<loc_dof_count<<std::endl;
      // std::cout<<"global_old_dof = = "<<global_old_dof<<std::endl;
      
      auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,loc_dof_count);
      // std::cout<<"ordered_dofs"<<std::endl;
      // std::cout<<ordered_dofs<<std::endl;
      // std::cout<<std::endl;
      cont =0;
      for(Integer m=0;m<ordered_dofs.size();m++)
      {
        (*dofmap_vec)[elem_id][ordered_dofs[cont]]=global_old_dof;
        // std::cout<<ordered_dofs[cont]<<std::endl;
        // std::cout<<(*dofmap_vec)[elem_id][ordered_dofs[cont]]<<std::endl;
        // std::cout<<global_dof_count<<std::endl;
        cont++;
        loc_dof_count++;
        global_old_dof++;
        dofs_count++;
      }
      // std::cout<<std::endl;

      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
     }
    


    }
    else
    {
    // if the entity has not been already visited, then create n_components new dofs
    if(flag[entity_id][0]==-1 || continuity==Discontinuous)
    {

       
      flag[entity_id][0]=elem_id;
      flag[entity_id][1]=entity_iter;
      // std::cout<<"elem_id = "<<elem_id<<std::endl;
      // std::cout<<"entity_iter = "<<entity_iter<<std::endl;
      for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
         for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
            {            
             (*dofmap_vec)[elem_id][loc_dof_count]=global_dof_count;
             space_dofs->push_back(global_dof_count);
             // std::cout<<global_dof_count<<std::endl;
             // std::cout<<"entity_dofs_iter = "<<entity_dofs_iter<<std::endl;
             // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
             // std::cout<<"(dofsmap)="<<(*dofmap_vec)[elem_id][loc_dof_count]<<std::endl;
             // std::cout<<" global_dof_count (new)= "<<global_dof_count<<std::endl;
             loc_dof_count++;
             global_dof_count++;
             dofs_count++;
             }
      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
    }
    else
    {
     // if the entity has been already visited, find the element (elem_id_tmp) in which was found 
     // and the corresponding entity_iter (iter_tmp)
     // we do not need any offset, because dofmap_vec contains  of the entity N, we move to the iter_tmp entity and find in dofmap_vec the already numbered dof
     const auto& elem_id_tmp=flag[entity_id][0];
     const auto& iter_tmp=flag[entity_id][1];
           auto old_dof=(*dofmap_vec)[elem_id_tmp][dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents];
      // std::cout<<" elem_id = "<<elem_id<<std::endl;
      // std::cout<<" entity_iter = "<<flag[entity_id][1]<<std::endl  ;  
     
     for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
        for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
           {           
            (*dofmap_vec)[elem_id][loc_dof_count]=old_dof;
            // std::cout<<old_dof<<std::endl;
             // std::cout<<"entity_dofs_iter = "<<entity_dofs_iter<<std::endl;
             // std::cout<<"old global_dof_count = "<<global_dof_count<<std::endl;
             loc_dof_count++;
             old_dof++;
           }
      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
    }

    }
     
  }
}

template<Integer N,typename FunctionSpace, typename FlagArray,typename Element,typename...Args1, long M,typename OS, typename SpaceDofs>
typename std::enable_if< -1==N,void>::type
ElementDofMap_LoopEntities3(const std::tuple<Args1...>& entitiestuple,
              FlagArray& flagtuples,
              Element& elem,
              const Integer& elem_id,
              std::shared_ptr<std::vector<Array<Integer, M>>> & dofmap_vec,
              Integer& global_dof_count,
              Integer &loc_dof_count,
              const OS &dofs_offset,
              SpaceDofs& space_dofs,
              Integer& dofs_count  )
{
}



template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename Element, typename...Args1,typename...Args2,typename...Args3,typename OS, typename SpaceDofs,typename ArrayNdofs>
typename std::enable_if< 0<sizeof...(FunctionSpaces),void>::type
ElementDofMap3(const std::tuple<Args1...>& entitiestuple,
                   std::tuple<Args2...>& flagtuples,
                   Element& elem,
                   const Integer& elem_id,
                   std::tuple<Args3...>& dofmap_vec,
                   Integer& global_dof_count,
                   Integer& loc_dof_count,
                   OS &dofs_offset,
                   SpaceDofs& space_dofs,
                   ArrayNdofs& array_ndofs   )
{
 static constexpr Integer M=sizeof...(FunctionSpaces);
 const auto& m1=std::get<M>(entitiestuple);
 auto &m2=std::get<M>(flagtuples);
 static constexpr Integer NN=std::tuple_size<remove_all_t<decltype(m1)>>::value-1;
 loc_dof_count=0;
  // std::cout<<"ElementDofMap3 K = "<<K <<std::endl;

 ElementDofMap_LoopEntities3<NN,ElemFunctionSpace<Elem,FunctionSpace>>
                            (m1,m2,elem,elem_id,tuple_get<K>(dofmap_vec),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),space_dofs[K],array_ndofs[K]);
  ElementDofMap3<K+1,Elem,FunctionSpaces...>(entitiestuple,flagtuples,elem,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,array_ndofs);

};




template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename Element, typename...Args1,typename...Args2,typename...Args3,typename OS, typename SpaceDofs,typename ArrayNdofs>
typename std::enable_if< 0==sizeof...(FunctionSpaces),void>::type
ElementDofMap3(const std::tuple<Args1...>& entitiestuple,
                   std::tuple<Args2...>& flagtuples,
                   Element& elem,
                   const Integer& elem_id,
                   std::tuple<Args3...>& dofmap_vec,
                   Integer& global_dof_count, 
                   Integer& loc_dof_count,
                   OS &dofs_offset,
                   SpaceDofs& space_dofs,
                   ArrayNdofs& array_ndofs )
{
 const auto& m1=std::get<0>(entitiestuple);
 auto& m2=std::get<0>(flagtuples);
 static constexpr Integer NN=std::tuple_size<remove_all_t<decltype(m1)>>::value-1;
 loc_dof_count=0;
 ElementDofMap_LoopEntities3<NN,ElemFunctionSpace<Elem,FunctionSpace>>
                           (m1,m2,elem,elem_id,tuple_get<K>(dofmap_vec),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),space_dofs[K],array_ndofs[K]);

};




















////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// DofsPerElemNums:                                                                                              //////////////
/////// For given FunctionSpaces, returns the total number of dofs for the givene element                             //////////////                                        
/////// For P3=LagrangeSpace of Order=3 and P1=LagrangeSpace of Order=1 in 5D                                         //////////////
/////// arr=[0,5,25,35], because dofs_per_entity(P3)=[1,2,1] and dofs_per_entity(P1)=[1]                              //////////////
/////// We have DofsPerElemNums=5 nodes * 1, 10 edges * 2, 10 triangles * 1, 5 nodes * 1 = 40                         //////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Elem, typename FunctionSpace,typename...FunctionSpaces>
class DofsPerElemNums
{ 
 public: 
 static constexpr Integer value= (DofsPerElemNums<Elem,FunctionSpaces...>::value + FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,FunctionSpace>>::value);

};

template<typename Elem,typename FunctionSpace>
class DofsPerElemNums<Elem,FunctionSpace>
{
 public: 
 static constexpr Integer value= FunctionSpaceDofsPerElem< ElemFunctionSpace<Elem,FunctionSpace>>::value;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// FunctionSpaceOffSetDofs:                                                                                     //////////////
/////// For FunctionSpace, returns an array with offsets                                                             //////////////
/////// For LagrangeSpace of Order=3 in 5D, arr=[0,5,25]+ addConst, because dofs_per_entity=[1,2,1]                  //////////////
/////// We have 5 nodes * 1, 10 edges * 2, 10 triangles * 1                                                          //////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<Integer MaxSize,typename FunctionSpace,Integer N,Integer AddConst=0>
typename  std::enable_if< (0<N)&&(N==MaxSize-1), void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
 {   
    arr[N] = AddConst + FunctionSpaceDofsPerElem<FunctionSpace,N-1>::value;
 }

template<Integer MaxSize,typename FunctionSpace,Integer N,Integer AddConst=0>
typename  std::enable_if< (0<N)&&(N<MaxSize-1), void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
{   
    arr[N] = AddConst + FunctionSpaceDofsPerElem<FunctionSpace,N-1>::value;
    FunctionSpaceOffSetDofs<MaxSize,FunctionSpace,N+1,AddConst>(arr);    
 }



template<Integer MaxSize,typename FunctionSpace,Integer N=0,Integer AddConst=0>
typename  std::enable_if< 0==N && 1<MaxSize, void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
{   
    arr[0] = AddConst + 0;
    FunctionSpaceOffSetDofs<MaxSize,FunctionSpace,N+1,AddConst>(arr);
}
template<Integer MaxSize,typename FunctionSpace,Integer N=0,Integer AddConst=0>
typename  std::enable_if< 0==N && 1==MaxSize, void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
{   
    arr[0] = AddConst + 0;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// OffSetDofs:                                                                                                   //////////////
/////// Generalize FunctionSpaceOffSetDofs to more FunctionSpaces , returns an array with offsets                     //////////////                                        
/////// For P3=LagrangeSpace of Order=3 and P1=LagrangeSpace of Order=1 in 5D                                         //////////////
/////// arr=[0,5,25,35], because dofs_per_entity(P3)=[1,2,1] and dofs_per_entity(P1)=[1]                              //////////////
/////// We have 5 nodes * 1, 10 edges * 2, 10 triangles * 1, 5 nodes * 1                                              //////////////
/////// PS: removed the add_const (not anymore necessary, each space is considered separately)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Elem,typename FunctionSpace,typename...FunctionSpaces>
struct OffsetDofsType
{
    using FS=ElemFunctionSpace<Elem,FunctionSpace>;
    static constexpr Integer entities_nums=FS::entity.size();
    using rest = typename OffsetDofsType<Elem, FunctionSpaces...>::type;
    using ens = std::array<Integer,entities_nums+1>;
    using tuple_ens=std::tuple<ens>;
    using type = decltype( std::tuple_cat(std::declval< tuple_ens >(), std::declval< rest >() ) );
};

template<typename Elem,typename FunctionSpace>
struct OffsetDofsType<Elem,FunctionSpace>
{
 using FS=ElemFunctionSpace<Elem,FunctionSpace>;
 static constexpr Integer entities_nums=FS::entity.size();
 using ens = std::array<Integer,entities_nums+1>;
 using type = typename std::tuple<ens>;
};


template<Integer M=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces>
typename std::enable_if< 0==sizeof...(FunctionSpaces),
                         typename OffsetDofsType<Elem,FunctionSpace>::type >::type
OffsetDofs()
{   
    using FS=ElemFunctionSpace<Elem,FunctionSpace>;
    constexpr const auto entities_nums=FS::entity.size();
    std::array<Integer,entities_nums+1> arr;
    FunctionSpaceOffSetDofs<entities_nums+1,FS,0,M>(arr);
    return std::tuple<decltype(arr)>(arr);
}

template<Integer M=0,typename Elem, typename FunctionSpace,typename...FunctionSpaces>
typename std::enable_if< 0<sizeof...(FunctionSpaces),
                         typename OffsetDofsType<Elem,FunctionSpace,FunctionSpaces...>::type >::type 
OffsetDofs()
{
    using FS=ElemFunctionSpace<Elem,FunctionSpace>;
    constexpr const auto entities_nums=FS::entity.size();
    std::array<Integer,entities_nums+1> arr;

    FunctionSpaceOffSetDofs<entities_nums+1,FS,0,M>(arr);
    // return std::tuple_cat(std::tuple<decltype(arr)>(arr),
    //                       OffsetDofs<M+FunctionSpaceDofsPerElem<FS,entities_nums-1>::value,
    //                                   Elem,FunctionSpaces...>());
    return std::tuple_cat(std::tuple<decltype(arr)>(arr),
                          OffsetDofs<M,
                                      Elem,FunctionSpaces...>());
}




template<Integer N,Integer M, typename FunctionSpace, typename...FunctionSpaces>
typename std::enable_if< 0==sizeof...(FunctionSpaces) , void >::type
function_space_info(Array<Array<Integer,4>, M> &array)
{
 array[N][0]=FunctionSpace::FEFamily;
 array[N][1]=FunctionSpace::Order;
 array[N][2]=FunctionSpace::Continuity;
 array[N][3]=FunctionSpace::NComponents;
};


template<Integer N,Integer M, typename FunctionSpace, typename...FunctionSpaces>
typename std::enable_if< 0<sizeof...(FunctionSpaces) , void >::type
function_space_info(Array<Array<Integer,4>, M> &array)
{
 array[N][0]=FunctionSpace::FEFamily;
 array[N][1]=FunctionSpace::Order;
 array[N][2]=FunctionSpace::Continuity;
 array[N][3]=FunctionSpace::NComponents;
 function_space_info<N+1,M,FunctionSpaces...>(array);
};





template<Integer N=0, typename...Args>
std::enable_if_t<(N==sizeof...(Args)),void>
resize_tuple_of_vector(std::tuple<Args...>& tuple, const Integer size)
{}

template<Integer N=0, typename...Args>
std::enable_if_t<(N<sizeof...(Args)),void>
resize_tuple_of_vector(std::tuple<Args...>& tuple, const Integer size)
{
 tuple_get<N>(tuple).resize(size);
 resize_tuple_of_vector<N+1>(tuple,size);
}










// template<typename FunctionSpace, typename...FunctionSpaces, typename MeshT, typename ...Args, typename OffSetT>//, typename DofMapT, typename OffSetT>
// void dofmap_fespace(const MeshT& mesh,
//              std::tuple<Args...>& dofmap_vec,
//              OffSetT& dofs_offset_arr,
//              Integer& global_dof_count,
//              const Array<Array<Integer,4> , 1 + sizeof...(FunctionSpaces)>& space_components,
//              Array<std::vector<std::vector<Integer>>,1+sizeof...(FunctionSpaces)>& space_dofs,
//              Array<Integer,1 + sizeof...(FunctionSpaces)>& array_ndofs
//              )
// {

//     using     Elem = typename MeshT::Elem; 
//     constexpr auto Dim=MeshT::Dim;
//     constexpr auto ManifoldDim=MeshT::ManifoldDim;

//     constexpr auto dofs_per_elem=DofsPerElemNums<Elem,FunctionSpace,FunctionSpaces...>::value;
//     // compute the connection node to elem (auxiliary tool, often used)
//     NodeToElem<Elem> node_2_elem(mesh);
//     const auto& node2elem=node_2_elem.val();
//     const auto& n_elements=mesh.n_elements();
//     const auto dofs_offset=OffsetDofs<0,Elem,FunctionSpace,FunctionSpaces...>();
//     const auto n_spaces=1+sizeof...(FunctionSpaces); 
//           auto entitiestuple=EntitiesOfFunctionSpaceTuple<Elem,FunctionSpace,FunctionSpaces...>(mesh,node2elem);
//           auto flagtuples= FlagTuple<Elem,FunctionSpace,FunctionSpaces...>(entitiestuple);


//     resize_tuple_of_vector(dofmap_vec,n_elements);//.resize(n_elements);   
//     global_dof_count=0;
//     // loop on all the elements
//     for(Integer space_id=0;space_id<n_spaces;space_id++)
//         space_dofs[space_id].resize(space_components[space_id][3]);
//     for(std::size_t i=0;i<array_ndofs.size();i++)
//       array_ndofs[i]=0;


//     for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
//     {
//       if(mesh.is_active(elem_iter))
//       {
//        // change it for smarter algorithms   
//        auto &elem_id=elem_iter;
//        Integer loc_dof_count=0;
//        ElementDofMap<0,Elem,FunctionSpace,FunctionSpaces...>
//                          (entitiestuple,flagtuples,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,array_ndofs);

//        // ElementDofMap<0,Elem,FunctionSpace,FunctionSpaces...>
//        //                   (entitiestuple,flagtuples,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,array_ndofs);

//       }

//     }

// tuple2array_dofs_offset<n_spaces-1>(dofs_offset,dofs_offset_arr);
// };



































template<Integer N=0, typename...Args>
std::enable_if_t<(N==sizeof...(Args)),void>
resize_tuple_of_ptr_vector(std::tuple<std::shared_ptr<Args>...>& tuple, const Integer size)
{}

template<Integer N=0, typename...Args>
std::enable_if_t<(N<sizeof...(Args)),void>
resize_tuple_of_ptr_vector(std::tuple<std::shared_ptr<Args>...>& tuple, const Integer size)
{
  using type=GetType<std::tuple<Args...>,N>;
 tuple_get<N>(tuple)=std::make_shared<type>();
 tuple_get<N>(tuple)->resize(size);
 resize_tuple_of_ptr_vector<N+1>(tuple,size);
}


template<typename FunctionSpace, typename...FunctionSpaces, typename MeshT, typename Dofmap, typename OffSetT>//, typename DofMapT, typename OffSetT>
void dofmap_fespace3(const MeshT& mesh,
             Dofmap& dofsdm_,
             OffSetT& dofs_offset_arr,
             Array<Integer,1 + sizeof...(FunctionSpaces)>& array_ndofs
             )
{

    using     Elem = typename MeshT::Elem; 
    constexpr auto Dim=Elem::Dim;
    constexpr auto ManifoldDim=Elem::ManifoldDim;

    constexpr auto dofs_per_elem=DofsPerElemNums<Elem,FunctionSpace,FunctionSpaces...>::value;

    auto& space_dofs=dofsdm_.space_dofs();
    auto& dofmap_vec=dofsdm_.dofmap();
    std::cout<<"dofs dm reserve begin"<<std::endl;
    dofsdm_.init_dofmap(mesh.n_elements());
    std::cout<<"dofs dm reserve end"<<std::endl;
    // compute the connection node to elem (auxiliary tool, often used)
    clock_t begin = clock();
    NodeToElem<MeshT> node_2_elem(mesh);
    const auto& node2elem=node_2_elem.val();
    const auto& n_elements=mesh.n_elements();
    const auto dofs_offset=OffsetDofs<0,Elem,FunctionSpace,FunctionSpaces...>();
    const auto n_spaces=1+sizeof...(FunctionSpaces); 
          auto entitiestuple=EntitiesOfFunctionSpaceTuple<Elem,FunctionSpace,FunctionSpaces...>(mesh,node2elem);
          auto flagtuples= FlagTuple<Elem,FunctionSpace,FunctionSpaces...>(entitiestuple);
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout<<"dofmap entities tuple elapsed_secs="<<elapsed_secs<<std::endl;
    // resize_tuple_of_ptr_vector(dofmap_vec,n_elements);   
    Integer global_dof_count=0;
    // loop on all the elements
    // for(Integer space_id=0;space_id<n_spaces;space_id++)
    // {
    //     space_dofs[space_id]=std::make_shared<std::vector<Integer>>();
    //     // space_dofs[space_id]->resize(space_components[space_id][3]);
    // }
    for(std::size_t i=0;i<array_ndofs.size();i++)
      array_ndofs[i]=0;

    begin = clock();
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {
      if(mesh.is_active(elem_iter))
      {
       // change it for smarter algorithms   
       auto& elem=mesh.elem(elem_iter);
       auto &elem_id=elem_iter;

       Integer loc_dof_count=0;
       ElementDofMap3<0,Elem,FunctionSpace,FunctionSpaces...>
                         (entitiestuple,flagtuples,elem,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,array_ndofs);

 
      }

    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout<<"dofmap tuple elapsed_secs="<<elapsed_secs<<std::endl;
    for(std::size_t i=0; i<dofsdm_.n_dofs().size() ;i++)
      dofsdm_.n_dofs()[i]=global_dof_count;


// tuple2array_dofs_offset<n_spaces-1>(dofs_offset,dofs_offset_arr);
};




























































template<typename MeshT>
class ElemEntityCollection;

template<typename MeshT>
class MeshAndEntity;










template<Integer N,typename FunctionSpace,typename FlagArray,typename Element,typename...Args1, long M,typename OS, typename SpaceDofs>
// template<Integer N,typename FunctionSpace,typename FlagArray,typename...Args1, long M,typename OS>
typename std::enable_if< -1<N,void>::type
ElementDofMap_LoopEntities4(const std::tuple<Args1...>& entitiestuple,
                           FlagArray& flagtuples,
                           Element& elem,
                           const Integer& elem_id, 
                           std::vector<Array<Integer, M>>&  dofmap_vec,
                           Integer& global_dof_count,
                           Integer& loc_dof_count,
                           const OS &dofs_offset,
                           SpaceDofs& space_dofs,
                           Integer& dofs_count)
{

std::cout<< " ElementDofMap_LoopEntities4 N = "<< N <<std::endl;
using Elem = typename FunctionSpace::Elem;
constexpr auto ManifoldDim=FunctionSpace::ManifoldDim;
constexpr auto continuity=FunctionSpace::Continuity;
constexpr auto n_components=FunctionSpace::NComponents;
constexpr auto entity_dim=FunctionSpace::entity[N];
constexpr auto dofs_per_entity=FunctionSpace::dofs_per_entity[N];
constexpr auto entity_points=entity_dim+1;
constexpr auto manifold_points=ManifoldDim+1;
constexpr auto combinations_nums=ElemEntityCombinations<Elem,entity_dim>::value;
const     auto& entity=std::get<entity_dim>(entitiestuple);
const     auto& elem2entity=entity.elem_2_entity(elem_id);
          auto& flag=flagtuples[N];

Integer cont;

// move to the entity of dimension entity[N-1]
ElementDofMap_LoopEntities4<N-1,FunctionSpace>
             (entitiestuple,flagtuples,elem,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,dofs_count);

const auto& nodes=elem.nodes;




// loop on all the entities of a given entity_dim
for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
   { 
    const auto& entity_id=elem2entity[entity_iter];
    
    if(DofsOrdering<entity_dim,FunctionSpace>::must_reorder)
    {
      auto ordered_entity_nodes=DofsOrdering<entity_dim,FunctionSpace>::value(nodes,entity_iter);
       // for(std::size_t i=0;i<ordered_entity_nodes.size();i++)
       //      {
       //        std::cout<<ordered_entity_nodes[i]<<" ";
       //      }

    cont=0;
    if(flag[entity_id][0]==-1 || continuity==Discontinuous)
    {

      // std::cout<< " elem id = "<< elem_id <<" NEW, iter ="<<entity_iter<<std::endl;

      auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,loc_dof_count);
      // std::cout<<"ordered_dofs"<<std::endl;
      // std::cout<<ordered_dofs<<std::endl;
       
      flag[entity_id][0]=elem_id;
      flag[entity_id][1]=entity_iter;
      // std::cout<<"elem_id = "<<elem_id<<std::endl;
      // std::cout<<"entity_iter = "<<entity_iter<<std::endl;
      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
      // std::cout<<std::endl;
      for(Integer m=0;m<ordered_dofs.size();m++)
      {
        dofmap_vec[elem_id][ordered_dofs[cont]]=global_dof_count;
        // std::cout<<(*dofmap_vec)[elem_id][ordered_dofs[cont]]<<std::endl;
        cont++;
        loc_dof_count++;
        // std::cout<<global_dof_count<<std::endl;
        global_dof_count++;
        dofs_count++;
      }
      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
      // std::cout<<std::endl;
     }
     else
     {
      // std::cout<< " elem id = "<< elem_id<<" NEW, iter ="<<entity_iter<<std::endl;

      
      
      const auto& elem_id_tmp=flag[entity_id][0];
      const auto& iter_tmp=flag[entity_id][1];
            auto& elem_dm=dofmap_vec[elem_id_tmp];
            auto ex_loc_dof_count=(dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents);
            auto global_old_dof=elem_dm[ex_loc_dof_count];
     

      for (int i = 0; i < n_components*dofs_per_entity; ++i)
      {
        const auto m= i +ex_loc_dof_count;
        if(elem_dm[m]<global_old_dof)
        {
          global_old_dof=elem_dm[m];
        }
      }

      // std::cout<<"old offset = = "<<dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents<<std::endl;
      // std::cout<<"dofs_offset= = "<<dofs_offset[N]<<std::endl;
      // std::cout<<"dofs_per_entity= = "<<dofs_per_entity<<std::endl;
      // std::cout<<"iter_tmp= = "<<iter_tmp<<std::endl;
      // std::cout<<"loc_dof_count = = "<<loc_dof_count<<std::endl;
      // std::cout<<"global_old_dof = = "<<global_old_dof<<std::endl;
      
      auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,loc_dof_count);
      // std::cout<<"ordered_dofs"<<std::endl;
      // std::cout<<ordered_dofs<<std::endl;
      // std::cout<<std::endl;
      cont =0;
      for(Integer m=0;m<ordered_dofs.size();m++)
      {
        dofmap_vec[elem_id][ordered_dofs[cont]]=global_old_dof;
        std::cout<< "(*dofmap_vec)[elem_id][ordered_dofs[cont]]="<<dofmap_vec[elem_id][ordered_dofs[cont]]<<std::endl;
        // std::cout<<ordered_dofs[cont]<<std::endl;
        // std::cout<<(*dofmap_vec)[elem_id][ordered_dofs[cont]]<<std::endl;
        // std::cout<<global_dof_count<<std::endl;
        cont++;
        loc_dof_count++;
        global_old_dof++;
        dofs_count++;
      }
      // std::cout<<std::endl;

      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
     }
    


    }
    else
    {
      std::cout<< "ElementDofMap_LoopEntities4 = "<<entity_iter<<"  ->"<<flag[entity_id][0]<<std::endl;
    // if the entity has not been already visited, then create n_components new dofs
    if(flag[entity_id][0]==-1 || continuity==Discontinuous)
    {

       
      flag[entity_id][0]=elem_id;
      flag[entity_id][1]=entity_iter;
      // std::cout<<"elem_id = "<<elem_id<<std::endl;
      // std::cout<<"entity_iter = "<<entity_iter<<std::endl;
      for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
         for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
            {            
             dofmap_vec[elem_id][loc_dof_count]=global_dof_count;
             space_dofs->push_back(global_dof_count);
             std::cout<< "(dofmap_vec)[elem_id][loc_dof_count]="<<dofmap_vec[elem_id][loc_dof_count]<<std::endl;
             // std::cout<<global_dof_count<<std::endl;
             // std::cout<<"entity_dofs_iter = "<<entity_dofs_iter<<std::endl;
             // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
             // std::cout<<"(dofsmap)="<<(*dofmap_vec)[elem_id][loc_dof_count]<<std::endl;
             // std::cout<<" global_dof_count (new)= "<<global_dof_count<<std::endl;
             loc_dof_count++;
             global_dof_count++;
             dofs_count++;
             }
      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
    }
    else
    {
     // if the entity has been already visited, find the element (elem_id_tmp) in which was found 
     // and the corresponding entity_iter (iter_tmp)
     // we do not need any offset, because dofmap_vec contains  of the entity N, we move to the iter_tmp entity and find in dofmap_vec the already numbered dof
     const auto& elem_id_tmp=flag[entity_id][0];
     const auto& iter_tmp=flag[entity_id][1];
           auto old_dof=dofmap_vec[elem_id_tmp][dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents];
      // std::cout<<" elem_id = "<<elem_id<<std::endl;
      // std::cout<<" entity_iter = "<<flag[entity_id][1]<<std::endl  ;  
     
     for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
        for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
           {           
            dofmap_vec[elem_id][loc_dof_count]=old_dof;
            // std::cout<<old_dof<<std::endl;
             // std::cout<<"entity_dofs_iter = "<<entity_dofs_iter<<std::endl;
             // std::cout<<"old global_dof_count = "<<global_dof_count<<std::endl;
            std::cout<< "(dofmap_vec)[elem_id][loc_dof_count]="<<dofmap_vec[elem_id][loc_dof_count]<<std::endl;
             loc_dof_count++;
             old_dof++;
           }
      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
    }

    }
     
  }



   std::cout<<" DOFMAP VEC ==" <<std::endl;
 
         for(int i=0;i<dofmap_vec.size();i++)
          {
            for(int j=0;j<dofmap_vec[i].size();j++)
          std::cout<<dofmap_vec[i][j]<<" ";
          std::cout<<std::endl;
        }
}

template<Integer N,typename FunctionSpace, typename FlagArray,typename Element,typename...Args1, long M,typename OS, typename SpaceDofs>
typename std::enable_if< -1==N,void>::type
ElementDofMap_LoopEntities4(const std::tuple<Args1...>& entitiestuple,
              FlagArray& flagtuples,
              Element& elem,
              const Integer& elem_id,
              std::vector<Array<Integer, M>>&  dofmap_vec,
              Integer& global_dof_count,
              Integer &loc_dof_count,
              const OS &dofs_offset,
              SpaceDofs& space_dofs,
              Integer& dofs_count  )
{
}






template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename Element, typename...Args1,typename...Args2,typename...Args3,typename OS, typename SpaceDofs,typename ArrayNdofs>
typename std::enable_if< 0<sizeof...(FunctionSpaces),void>::type
ElementDofMap4(const std::tuple<Args1...>& entitiestuple,
                   std::tuple<Args2...>& flagtuples,
                   Element& elem,
                   const Integer& elem_id,
                   std::tuple<Args3...>& dofmap_vec,
                   Integer& global_dof_count,
                   Integer& loc_dof_count,
                   OS &dofs_offset,
                   SpaceDofs& space_dofs,
                   ArrayNdofs& array_ndofs   )
{
 static constexpr Integer M=sizeof...(FunctionSpaces);
 const auto& m1=std::get<M>(entitiestuple);
 auto &m2=std::get<M>(flagtuples);
 static constexpr Integer NN=std::tuple_size<remove_all_t<decltype(m2)>>::value-1;
 loc_dof_count=0;
  // std::cout<<"ElementDofMap3 K = "<<K <<std::endl;

 ElementDofMap_LoopEntities4<NN,ElemFunctionSpace<Elem,FunctionSpace>>
                            (entitiestuple,m2,elem,elem_id,(*tuple_get<K>(dofmap_vec)),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),space_dofs[K],array_ndofs[K]);


 auto & dm=*tuple_get<K>(dofmap_vec);
 std::cout<<" elem dofmap 4 K==" << K <<std::endl;
 std::cout<<" dm.size()==" << dm.size() <<std::endl;
 
 for(int i=0;i<dm.size();i++)
  {
    for(int j=0;j<dm[i].size();j++)
  std::cout<<dm[i][j]<<" ";
  std::cout<<std::endl;

}

  ElementDofMap4<K+1,Elem,FunctionSpaces...>(entitiestuple,flagtuples,elem,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,array_ndofs);

};




template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename Element, typename...Args1,typename...Args2,typename...Args3,typename OS, typename SpaceDofs,typename ArrayNdofs>
typename std::enable_if< 0==sizeof...(FunctionSpaces),void>::type
ElementDofMap4(const std::tuple<Args1...>& entitiestuple,
                   std::tuple<Args2...>& flagtuples,
                   Element& elem,
                   const Integer& elem_id,
                   std::tuple<Args3...>& dofmap_vec,
                   Integer& global_dof_count, 
                   Integer& loc_dof_count,
                   OS &dofs_offset,
                   SpaceDofs& space_dofs,
                   ArrayNdofs& array_ndofs )
{
 // const auto& m1=std::get<0>(entitiestuple);
 auto& m2=std::get<0>(flagtuples);
 static constexpr Integer NN=std::tuple_size<remove_all_t<decltype(m2)>>::value-1;
 loc_dof_count=0;


 ElementDofMap_LoopEntities4<NN,ElemFunctionSpace<Elem,FunctionSpace>>
                           (entitiestuple,m2,elem,elem_id,(*tuple_get<K>(dofmap_vec)),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),space_dofs[K],array_ndofs[K]);

 auto & dm=*tuple_get<K>(dofmap_vec);
 std::cout<<" elem dofmap 4 K==" << K <<std::endl;
 std::cout<<" dm.size()==" << dm.size() <<std::endl;
 
 for(int i=0;i<dm.size();i++)
  {
    for(int j=0;j<dm[i].size();j++)
  std::cout<<dm[i][j]<<" ";
  std::cout<<std::endl;

}
};




template<typename FunctionSpace,  Integer Mmax, Integer M,  typename ...Args>
typename std::enable_if< (M>Mmax),void >::type
initialize_flag_tuple_entities_aux4(std::tuple< Args...> const &tuple, 
                               std::array<std::vector< std::array<Integer,2> >, 
                               FunctionSpace::entity.size() > & entity_found)
    {static_assert(M>=0," the tuple must have non negative length ");};

template<typename FunctionSpace,  Integer Mmax,Integer M, typename...Args>
typename std::enable_if< (M<=Mmax),void >::type
initialize_flag_tuple_entities_aux4(std::tuple< Args...> const &tuple,
                               std::array<std::vector< std::array<Integer,2> >, FunctionSpace::entity.size() > & entity_found)
{
     static_assert(M>=0," the tuple must have non negative length ");
     static constexpr Integer entity_dim=FunctionSpace::entity[M];
     // maybe 
     Integer entity_length=tuple_get<entity_dim>(tuple).size();

     entity_found[M].resize(entity_length,{-1,-1});

     std::cout<<"entity_found[M] size==" <<entity_found[M].size()<<std::endl;
     initialize_flag_tuple_entities_aux4<FunctionSpace,Mmax,M+1>(tuple,entity_found);       
};


template<typename FunctionSpace, typename...Args>
 void initialize_flag_tuple_entities4(
  std::tuple< Args...> const &tuple,
  std::array<std::vector< std::array<Integer,2> >, FunctionSpace::entity.size() > & entity_found)
{
     static constexpr Integer entities_nums=FunctionSpace::entity.size();
     initialize_flag_tuple_entities_aux4<FunctionSpace,entities_nums-1,0>(tuple,entity_found);       
};




template<typename Elem,typename FunctionSpace,typename...FunctionSpaces>
struct FlagTupleType4
{
     using FS=ElemFunctionSpace<Elem,FunctionSpace>;
     static constexpr Integer FEFamily=FS::FEFamily;
     static constexpr Integer Order=FS::Order; 
     static constexpr Integer entities_nums=FS::entity.size();
     using rest = typename FlagTupleType4<Elem,FunctionSpaces...>::type;
     using ens  = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
     using type = decltype( std::tuple_cat( std::declval< rest >(), 
                                            std::declval< std::tuple<ens> >() ) );
};


template<typename Elem,typename FunctionSpace>
struct FlagTupleType4<Elem,FunctionSpace>
{
 using FS=ElemFunctionSpace<Elem,FunctionSpace>;
 static constexpr Integer FEFamily=FS::FEFamily;
 static constexpr Integer Order=FS::Order; 
 static constexpr Integer entities_nums=FS::entity.size();
 using ens =  typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
 using type = typename std::tuple<ens>;
};


template<typename Elem,typename BaseFunctionSpace,typename...BaseFunctionSpaces, typename MeshT>
typename std::enable_if< 0<sizeof...(BaseFunctionSpaces), 
                         typename FlagTupleType4<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::type>::type
FlagTuple4(ElemEntityCollection<MeshT>& entitycollection)
{    
      using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
      static constexpr Integer entities_nums=FunctionSpace::entity.size();
      using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
      type ens;

      initialize_flag_tuple_entities4<FunctionSpace>(entitycollection.tuple_entities(),ens);
      return std::tuple_cat(FlagTuple4<Elem,BaseFunctionSpaces...>(entitycollection),
                            std::tuple<type>(ens));
}

template<typename Elem,typename BaseFunctionSpace, typename...BaseFunctionSpaces, typename MeshT>
 typename std::enable_if< 0==sizeof...(BaseFunctionSpaces), 
                          typename FlagTupleType4<Elem,BaseFunctionSpace>::type>::type
FlagTuple4(ElemEntityCollection<MeshT>& entitycollection)
{   
      using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
      static constexpr Integer entities_nums=FunctionSpace::entity.size();
      using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
      type ens;
      initialize_flag_tuple_entities4<FunctionSpace>(entitycollection.tuple_entities(),ens);
      return std::tuple<type>(ens);
}



template<typename BaseFunctionSpace, typename...BaseFunctionSpaces, typename MeshT, typename Dofmap, typename OffSetT>//, typename DofMapT, typename OffSetT>
void dofmap_fespace4( MeshAndEntity<MeshT>& mesh_and_entity,
             Dofmap& dofsdm_,
             OffSetT& dofs_offset_arr,
             Array<Integer,1 + sizeof...(BaseFunctionSpaces)>& array_ndofs
             )
{

    using     Elem = typename MeshT::Elem; 
    constexpr auto Dim=MeshT::Dim;
    constexpr auto ManifoldDim=MeshT::ManifoldDim;

    constexpr auto dofs_per_elem=DofsPerElemNums<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;

    auto mesh_ptr=mesh_and_entity.mesh_ptr();
    auto& space_dofs=dofsdm_.space_dofs();
    auto& dofmap_vec=dofsdm_.dofmap();
    std::cout<<"dofs dm reserve begin"<<std::endl;
    dofsdm_.init_dofmap(mesh_ptr->n_elements());
    std::cout<<"dofs dm reserve end"<<std::endl;
    // compute the connection node to elem (auxiliary tool, often used)
    clock_t begin = clock();
    
    auto& entities_collection=mesh_and_entity.entities_collection();
    init<Elem,BaseFunctionSpace,BaseFunctionSpaces...>(mesh_and_entity);
//     NodeToElem<MeshT> node_2_elem(mesh);
//     const auto& node2elem=node_2_elem.val();
//     const auto& n_elements=mesh.n_elements();
    const auto dofs_offset=OffsetDofs<0,Elem,BaseFunctionSpace,BaseFunctionSpaces...>();
//     const auto n_spaces=1+sizeof...(BaseFunctionSpaces); 
//           auto entitiestuple=EntitiesOfFunctionSpaceTuple<Elem,BaseFunctionSpace,BaseFunctionSpaces...>(mesh,node2elem);
          auto flagtuples= FlagTuple4<Elem,BaseFunctionSpace,BaseFunctionSpaces...>(entities_collection);
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout<<"dofmap entities tuple elapsed_secs="<<elapsed_secs<<std::endl;
    // resize_tuple_of_ptr_vector(dofmap_vec,n_elements);   
    Integer global_dof_count=0;
    // loop on all the elements
    // for(Integer space_id=0;space_id<n_spaces;space_id++)
    // {
    //     space_dofs[space_id]=std::make_shared<std::vector<Integer>>();
    //     // space_dofs[space_id]->resize(space_components[space_id][3]);
    // }
    for(std::size_t i=0;i<array_ndofs.size();i++)
      array_ndofs[i]=0;

    auto bisection_ptr=mesh_and_entity.bisection_ptr();
    auto level=mesh_and_entity.level();
    begin = clock();
    for(Integer elem_iter=0;elem_iter<mesh_ptr->n_elements();elem_iter++)
    {
      if(!elem_belongs_to_level(mesh_ptr,elem_iter,level,bisection_ptr)) continue;
      // if(mesh_ptr->is_active(elem_iter))
      // {
       // change it for smarter algorithms   
      std::cout<<" dofmap elem == "<< elem_iter <<std::endl;
       auto& elem=mesh_ptr->elem(elem_iter);
       auto &elem_id=elem_iter;

       Integer loc_dof_count=0;
       ElementDofMap4<0,Elem,BaseFunctionSpace,BaseFunctionSpaces...>
                         (entities_collection.tuple_entities(),flagtuples,elem,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,array_ndofs);

 
      // }

    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout<<"dofmap tuple elapsed_secs="<<elapsed_secs<<std::endl;
    for(std::size_t i=0; i<dofsdm_.n_dofs().size() ;i++)
      dofsdm_.n_dofs()[i]=global_dof_count;


// tuple2array_dofs_offset<n_spaces-1>(dofs_offset,dofs_offset_arr);
}

















































































































template<typename MeshT,Integer N>
class ConnectivitySimpliacialMap;

template<typename MeshT>
class ConnectivitySimpliacialMapCollection;



template<Integer N,typename FunctionSpace,typename FlagArray,typename Element,typename MeshT, Integer P, long M,typename OS, typename SpaceDofs>
// template<Integer N,typename FunctionSpace,typename FlagArray,typename...Args1, long M,typename OS>
typename std::enable_if< -1<N,void>::type
ElementDofMap_LoopEntities5(const ConnectivitySimpliacialMap<MeshT,P>& entitiestuple,
                           FlagArray& flagtuples,
                           Element& elem,
                           const Integer& elem_id, 
                           std::vector<Array<Integer, M>>&  dofmap_vec,
                           Integer& global_dof_count,
                           Integer& loc_dof_count,
                           const OS &dofs_offset,
                           SpaceDofs& space_dofs,
                           Integer& dofs_count)
{

std::cout<< " ElementDofMap_LoopEntities5 N = "<< N <<std::endl;
using Elem = typename FunctionSpace::Elem;
constexpr auto ManifoldDim=FunctionSpace::ManifoldDim;
constexpr auto continuity=FunctionSpace::Continuity;
constexpr auto n_components=FunctionSpace::NComponents;
constexpr auto entity_dim=FunctionSpace::entity[N];
constexpr auto dofs_per_entity=FunctionSpace::dofs_per_entity[N];
constexpr auto entity_points=entity_dim+1;
constexpr auto manifold_points=ManifoldDim+1;
constexpr auto combinations_nums=ElemEntityCombinations<Elem,entity_dim>::value;
// const     auto& entity=std::get<entity_dim>(entitiestuple);
// const     auto& elem2entity=entity.elem_2_entity(elem_id);
          auto& flag=tuple_get<N>(flagtuples);


std::array<Integer,entity_points> entity_nodes;

Integer cont;

// move to the entity of dimension entity[N-1]
ElementDofMap_LoopEntities5<N-1,FunctionSpace>
             (entitiestuple,flagtuples,elem,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,dofs_count);

std::cout<<"pre nodes ="<<std::endl;
const auto& nodes=elem.nodes;

std::cout<<"after nodes ="<<std::endl;

Integer entity_[entity_points];

std::cout<<"after entity_ ="<<std::endl;


// loop on all the entities of a given entity_dim
for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
   { 

    Combinations<manifold_points, entity_points>::generate(entity_iter,entity_);
    for(std::size_t i=0;i<entity_points;i++)
      entity_nodes[i]=nodes[entity_[i]];
    std::sort(entity_nodes.begin(),entity_nodes.end());
    std::cout<<"entity_iter="<<entity_iter<<std::endl;
  

    auto ordered_entity_nodes=DofsOrdering<entity_dim,FunctionSpace>::value(nodes,entity_iter);

    // consider more dofs per entity (RT1)
    if(DofsOrdering<entity_dim,FunctionSpace>::must_reorder)
    {

      // std::cout<<"must_reorder=1"<<std::endl;
      auto ordered_dofs=DofsOrdering<entity_dim,FunctionSpace>::value(ordered_entity_nodes,entity_iter);
      
      cont=0;
      // if entity not found yet
      if(!flag[entity_nodes] || continuity==Discontinuous )
      {

      // std::cout<<"!flag[ordered_entity_nodes]"<<std::endl;
      flag[entity_nodes]=std::make_shared<std::array<Integer,2>>(std::array<Integer,2>{elem_id,entity_iter});
      // std::cout<<"after !flag[ordered_entity_nodes]"<<std::endl;
      // flag[ordered_entity_nodes][0]=elem_id;
      // flag[ordered_entity_nodes][1]=entity_iter;
      for(Integer m=0;m<ordered_dofs.size();m++)
      {
        // std::cout<<"m="<<m<<std::endl;
        // std::cout<<"dofmap_vec[elem_id]"<<dofmap_vec.size()<<std::endl;
        // std::cout<<"dofmap_vec[elem_id]"<<dofmap_vec[elem_id].size()<<std::endl;
        
        dofmap_vec[elem_id][ordered_dofs[cont]]=global_dof_count;
        // std::cout<<(*dofmap_vec)[elem_id][ordered_dofs[cont]]<<std::endl;
        cont++;
        loc_dof_count++;
        // std::cout<<global_dof_count<<std::endl;
        global_dof_count++;
        dofs_count++;
      }

      }
      else
      {

       // std::cout<<"flag[entity_nodes]"<<std::endl;
      const auto& elem_id_tmp=(*flag[entity_nodes])[0];
      const auto& iter_tmp=(*flag[entity_nodes])[1];
            auto& elem_dm=dofmap_vec[elem_id_tmp];
            auto ex_loc_dof_count=(dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents);
            auto global_old_dof=elem_dm[ex_loc_dof_count];
     

      for (int i = 0; i < n_components*dofs_per_entity; ++i)
      {
        const auto m= i +ex_loc_dof_count;
        if(elem_dm[m]<global_old_dof)
        {
          global_old_dof=elem_dm[m];
        }
      }

      auto ordered_dofs=DofsOrdered<entity_dim,FunctionSpace>::local_dofs(ordered_entity_nodes,ex_loc_dof_count);
      
      cont =0;
      std::cout<<"ordered_dofs="<<ordered_dofs<<std::endl;
      for(Integer m=0;m<ordered_dofs.size();m++)
      {
        std::cout<<"dofmap_vec[elem_id].size()="<<dofmap_vec[elem_id].size()<<std::endl;
        std::cout<<"ordered_dofs[cont]="<<ordered_dofs[cont]<<std::endl;
        dofmap_vec[elem_id][ordered_dofs[cont]]=global_old_dof;
        // std::cout<< "(*dofmap_vec)[elem_id][ordered_dofs[cont]]="<<dofmap_vec[elem_id][ordered_dofs[cont]]<<std::endl;
        cont++;
        loc_dof_count++;
        global_old_dof++;
        dofs_count++;
      }



      }
     


    }
    // otherwise
    else
    {
      // std::cout<<"must_reorder=0"<<std::endl;

    // if the entity has not been already visited, then create n_components new dofs
    if(!flag[entity_nodes] || continuity==Discontinuous)
    {
      // std::cout<<"!flag[ordered_entity_nodes]"<<std::endl;
      flag[entity_nodes]=std::make_shared<std::array<Integer,2>>(std::array<Integer,2>{elem_id,entity_iter});

      // flag[ordered_entity_nodes][0]=elem_id;
      // flag[ordered_entity_nodes][1]=entity_iter;


      // std::cout<<"elem_id = "<<elem_id<<std::endl;
      // std::cout<<"entity_iter = "<<entity_iter<<std::endl;
      for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
         for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
            {            
             dofmap_vec[elem_id][loc_dof_count]=global_dof_count;
             space_dofs->push_back(global_dof_count);
             // std::cout<< "(dofmap_vec)[elem_id][loc_dof_count]="<<dofmap_vec[elem_id][loc_dof_count]<<std::endl;
             // std::cout<<global_dof_count<<std::endl;
             // std::cout<<"entity_dofs_iter = "<<entity_dofs_iter<<std::endl;
             // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
             // std::cout<<"(dofsmap)="<<(*dofmap_vec)[elem_id][loc_dof_count]<<std::endl;
             // std::cout<<" global_dof_count (new)= "<<global_dof_count<<std::endl;
             loc_dof_count++;
             global_dof_count++;
             dofs_count++;
             }
      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
    }
    else
    {
     // std::cout<<"flag[ordered_entity_nodes]"<<std::endl;
     // if the entity has been already visited, find the element (elem_id_tmp) in which was found 
     // and the corresponding entity_iter (iter_tmp)
     // we do not need any offset, because dofmap_vec contains  of the entity N, we move to the iter_tmp entity and find in dofmap_vec the already numbered dof
     const auto& elem_id_tmp=(*flag[entity_nodes])[0];
     const auto& iter_tmp=(*flag[entity_nodes])[1];
           auto old_dof=dofmap_vec[elem_id_tmp][dofs_offset[N]+dofs_per_entity*iter_tmp*FunctionSpace::NComponents];
      // std::cout<<" elem_id = "<<elem_id<<std::endl;
      // std::cout<<" entity_iter = "<<flag[entity_id][1]<<std::endl  ;  
     
     for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
        for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
           {           
            dofmap_vec[elem_id][loc_dof_count]=old_dof;
            // std::cout<<old_dof<<std::endl;
             // std::cout<<"entity_dofs_iter = "<<entity_dofs_iter<<std::endl;
             // std::cout<<"old global_dof_count = "<<global_dof_count<<std::endl;
            // std::cout<< "(dofmap_vec)[elem_id][loc_dof_count]="<<dofmap_vec[elem_id][loc_dof_count]<<std::endl;
             loc_dof_count++;
             old_dof++;
           }
      // std::cout<<"global_dof_count = "<<global_dof_count<<std::endl;
    }
 
    }

  }


   // std::cout<<" DOFMAP VEC5 ==" <<std::endl;
   // std::cout<<" dofmap_vec.size() ==" <<dofmap_vec.size()<<std::endl;
 
        //  for(int i=0;i<dofmap_vec.size();i++)
        //   {
        //     for(int j=0;j<dofmap_vec[i].size();j++)
        //   std::cout<<dofmap_vec[i][j]<<" ";
        //   std::cout<<std::endl;
        // }
  // std::cout<<"after  dofmap_vec.size() ==" <<dofmap_vec.size()<<std::endl;
}

template<Integer N,typename FunctionSpace, typename FlagArray,typename Element,typename MeshT, Integer P, long M,typename OS, typename SpaceDofs>
typename std::enable_if< -1==N,void>::type
ElementDofMap_LoopEntities5(const ConnectivitySimpliacialMap<MeshT,P>& entitiestuple,
              FlagArray& flagtuples,
              Element& elem,
              const Integer& elem_id,
              std::vector<Array<Integer, M>>&  dofmap_vec,
              Integer& global_dof_count,
              Integer &loc_dof_count,
              const OS &dofs_offset,
              SpaceDofs& space_dofs,
              Integer& dofs_count  )
{
  std::cout<<"ElementDofMap_LoopEntities5 K ==-1" <<std::endl;
}






template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename Element, typename...Args1,typename...Args2,typename...Args3,typename OS, typename SpaceDofs,typename ArrayNdofs>
typename std::enable_if< 0<sizeof...(FunctionSpaces),void>::type
ElementDofMap5(const std::tuple<Args1...>& entitiestuple,
                   std::tuple<Args2...>& flagtuples,
                   Element& elem,
                   const Integer& elem_id,
                   std::tuple<Args3...>& dofmap_vec,
                   Integer& global_dof_count,
                   Integer& loc_dof_count,
                   OS &dofs_offset,
                   SpaceDofs& space_dofs,
                   ArrayNdofs& array_ndofs   )
{





    std::cout<<"__________qui0111111=2=" << K <<std::endl;

 static constexpr Integer M=sizeof...(FunctionSpaces);
 const auto& entity_tuple=std::get<M>(entitiestuple);
 auto &flagtuple=std::get<M>(flagtuples);
 static constexpr Integer NN=std::tuple_size<remove_all_t<decltype(flagtuple)>>::value-1;
 loc_dof_count=0;

 ElementDofMap_LoopEntities5<NN,ElemFunctionSpace<Elem,FunctionSpace>>
                            (entity_tuple,flagtuple,elem,elem_id,(*tuple_get<K>(dofmap_vec)),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),space_dofs[K],array_ndofs[K]);

 std::cout<<"__________qui==" << K <<std::endl;
 auto & dm=*tuple_get<K>(dofmap_vec);
 std::cout<<" elem dofmap 5 K==" << K <<std::endl;
 std::cout<<" dm.size()==" << dm.size() <<std::endl;

const auto& m3=entity_tuple.elem2dofs(elem_id);
   std::cout<<"elem2dofs=="<<std::endl;
    for(int i=0;i<m3.size();i++)
      std::cout<<m3[i]<<" "<<std::endl;
   std::cout<<std::endl;

 for(int i=0;i<dm.size();i++)
  {
    for(int j=0;j<dm[i].size();j++)
  std::cout<<dm[i][j]<<" ";
  std::cout<<std::endl;

}

  ElementDofMap5<K+1,Elem,FunctionSpaces...>(entitiestuple,flagtuples,elem,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset,space_dofs,array_ndofs);

};




template<Integer K=0,typename Elem,typename FunctionSpace,typename...FunctionSpaces, typename Element, typename...Args1,typename...Args2,typename...Args3,typename OS, typename SpaceDofs,typename ArrayNdofs>
typename std::enable_if< 0==sizeof...(FunctionSpaces),void>::type
ElementDofMap5(const std::tuple<Args1...>& entitiestuple,
                   std::tuple<Args2...>& flagtuples,
                   Element& elem,
                   const Integer& elem_id,
                   std::tuple<Args3...>& dofmap_vec,
                   Integer& global_dof_count, 
                   Integer& loc_dof_count,
                   OS &dofs_offset,
                   SpaceDofs& space_dofs,
                   ArrayNdofs& array_ndofs )
{


  std::cout<<"__________qui00000=2=" << K <<std::endl;
 const auto& entity_tuple=std::get<0>(entitiestuple);
 auto& flagtuple=std::get<0>(flagtuples);
 static constexpr Integer NN=std::tuple_size<remove_all_t<decltype(flagtuple)>>::value-1;
 loc_dof_count=0;


 ElementDofMap_LoopEntities5<NN,ElemFunctionSpace<Elem,FunctionSpace>>
                           (entity_tuple,flagtuple,elem,elem_id,(*tuple_get<K>(dofmap_vec)),global_dof_count,loc_dof_count,std::get<K>(dofs_offset),space_dofs[K],array_ndofs[K]);

 std::cout<<"__________qui2=2=" << K <<std::endl;

 auto & dm=*tuple_get<K>(dofmap_vec);
 std::cout<<" elem dofmap 4 K==" << K <<std::endl;
 std::cout<<" dm.size()==" << dm.size() <<std::endl;
 const auto& m3=entity_tuple.elem2dofs(elem_id);

   std::cout<<"elem2dofs=="<<std::endl;
    for(int i=0;i<m3.size();i++)
      std::cout<<m3[i]<<" "<<std::endl;
   std::cout<<std::endl;

 for(int i=0;i<dm.size();i++)
  {
    for(int j=0;j<dm[i].size();j++)
  std::cout<<dm[i][j]<<" ";
  std::cout<<std::endl;

}
};








template<typename FunctionSpace,  Integer Mmax, Integer M,  typename ...Args>
typename std::enable_if< (M>Mmax),void >::type
initialize_flag_tuple_entities_aux5(std::tuple< Args...> const &tuple, 
                               std::array<std::vector< std::array<Integer,2> >, 
                               FunctionSpace::entity.size() > & entity_found)
    {static_assert(M>=0," the tuple must have non negative length ");};

template<typename FunctionSpace,  Integer Mmax,Integer M, typename...Args>
typename std::enable_if< (M<=Mmax),void >::type
initialize_flag_tuple_entities_aux5(std::tuple< Args...> const &tuple,
                               std::array<std::vector< std::array<Integer,2> >, FunctionSpace::entity.size() > & entity_found)
{
     static_assert(M>=0," the tuple must have non negative length ");
     static constexpr Integer entity_dim=FunctionSpace::entity[M];
     // maybe 
     Integer entity_length=tuple_get<entity_dim>(tuple).size();

     entity_found[M].resize(entity_length,{-1,-1});

     std::cout<<"entity_found[M] size==" <<entity_found[M].size()<<std::endl;
     initialize_flag_tuple_entities_aux5<FunctionSpace,Mmax,M+1>(tuple,entity_found);       
};


template<typename FunctionSpace, typename...Args>
 void initialize_flag_tuple_entities5(
  std::tuple< Args...> const &tuple,
  std::array<std::vector< std::array<Integer,2> >, FunctionSpace::entity.size() > & entity_found)
{
     static constexpr Integer entities_nums=FunctionSpace::entity.size();
     initialize_flag_tuple_entities_aux5<FunctionSpace,entities_nums-1,0>(tuple,entity_found);       
};




template<typename Elem,typename FunctionSpace,typename...FunctionSpaces>
struct FlagTupleType5
{
     using FS=ElemFunctionSpace<Elem,FunctionSpace>;
     static constexpr Integer FEFamily=FS::FEFamily;
     static constexpr Integer Order=FS::Order; 
     static constexpr Integer entities_nums=FS::entity.size();
     using rest = typename FlagTupleType5<Elem,FunctionSpaces...>::type;
     using ens  = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
     using type = decltype( std::tuple_cat( std::declval< rest >(), 
                                            std::declval< std::tuple<ens> >() ) );
};


template<typename Elem,typename FunctionSpace>
struct FlagTupleType5<Elem,FunctionSpace>
{
 using FS=ElemFunctionSpace<Elem,FunctionSpace>;
 static constexpr Integer FEFamily=FS::FEFamily;
 static constexpr Integer Order=FS::Order; 
 static constexpr Integer entities_nums=FS::entity.size();
 using ens =  typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
 using type = typename std::tuple<ens>;
};


template<typename Elem,typename BaseFunctionSpace,typename...BaseFunctionSpaces, typename MeshT>
typename std::enable_if< 0<sizeof...(BaseFunctionSpaces), 
                         typename FlagTupleType4<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::type>::type
FlagTuple5(ConnectivitySimpliacialMapCollection<MeshT>& entities)
{    
      using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
      static constexpr Integer entities_nums=FunctionSpace::entity.size();
      using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
      type ens;

      initialize_flag_tuple_entities5<FunctionSpace>(entities.tuple_entities(),ens);
      return std::tuple_cat(FlagTuple5<Elem,BaseFunctionSpaces...>(entities),
                            std::tuple<type>(ens));
}

template<typename Elem,typename BaseFunctionSpace, typename...BaseFunctionSpaces, typename MeshT>
 typename std::enable_if< 0==sizeof...(BaseFunctionSpaces), 
                          typename FlagTupleType4<Elem,BaseFunctionSpace>::type>::type
FlagTuple5(ConnectivitySimpliacialMapCollection<MeshT>& entities)
{   
      using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
      static constexpr Integer entities_nums=FunctionSpace::entity.size();
      using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
      type ens;
      initialize_flag_tuple_entities5<FunctionSpace>(entities.tuple_entities(),ens);
      return std::tuple<type>(ens);
}










  // template <Integer N>
  // struct FlagEntitiesTuple;

  template<Integer N>
  struct FlagSimplexEntitiesTuple
   {
    using type=std::map<std::array<Integer,N+1>,std::shared_ptr<std::array<Integer,2>>>;
   };


  template<typename FunctionSpace,Integer Nmax,Integer N>
  struct FlagSimplexEntitiesFunctionSpaceTupleHelper;


  template<typename FunctionSpace,Integer Nmax>
  struct FlagSimplexEntitiesFunctionSpaceTupleHelper<FunctionSpace,Nmax,Nmax>
   {
    using single_type=typename FlagSimplexEntitiesTuple<FunctionSpace::entity[Nmax]>::type;
    using type=std::tuple<single_type>;
   };



  template<typename FunctionSpace,Integer Nmax,Integer N>
  struct FlagSimplexEntitiesFunctionSpaceTupleHelper
   {
     using single_type=typename FlagSimplexEntitiesTuple<FunctionSpace::entity[N]>::type;
     using rest=typename FlagSimplexEntitiesFunctionSpaceTupleHelper<FunctionSpace,Nmax,N+1>::type;
     using type = decltype( std::tuple_cat( std::declval< std::tuple<single_type> >() ,
                                            std::declval< rest >() 
                                            ) );

   };


  template<typename FunctionSpace>
  struct FlagSimplexEntitiesFunctionSpaceTuple
   {
    using type=typename FlagSimplexEntitiesFunctionSpaceTupleHelper<FunctionSpace,FunctionSpace::entity.size()-1,0>::type;
   };
 //  template <Integer N>
 //  struct FlagEntitiesTuple
 //  {
 //   using rest = typename FlagEntitiesTuple<N-1>::type; 
 //   using single_type=std::map<std::array<Integer,N+1>,Integer>;
 //   using tuple_single_type=std::tuple<single_type>;
 //   using type = decltype( std::tuple_cat(std::declval< rest >() , std::declval< tuple_single_type >()) );
 //  };

 //  template <Integer N>
 // using FlagEntities=typename FlagEntitiesTuple<N>::type;



template<typename Elem,typename BaseFunctionSpace,typename...BaseFunctionSpaces>
struct FlagTupleType6
{
     using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
     using rest = typename FlagTupleType6<Elem,BaseFunctionSpaces...>::type;
     using ens  = typename FlagSimplexEntitiesFunctionSpaceTuple<FunctionSpace>::type;
     using type = decltype( std::tuple_cat( std::declval< std::tuple<ens> >(),std::declval< rest >()
                                             ) );
};


template<typename Elem,typename BaseFunctionSpace>
struct FlagTupleType6<Elem,BaseFunctionSpace>
{
 using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
 using single_type = typename FlagSimplexEntitiesFunctionSpaceTuple<FunctionSpace>::type;
 using type = typename std::tuple<single_type>;
};


template<typename Elem,typename...BaseFunctionSpaces>
using FlagTuple6=typename FlagTupleType6<Elem,BaseFunctionSpaces...>::type;



template<typename BaseFunctionSpace, typename...BaseFunctionSpaces, typename MeshT, typename Dofmap, typename OffSetT>//, typename DofMapT, typename OffSetT>
void dofmap_fespace5( ConnectivitySimpliacialMapCollection<MeshT>& entities,
             Dofmap& dofsdm_,
             OffSetT& dofs_offset_arr,
             Array<Integer,1 + sizeof...(BaseFunctionSpaces)>& array_ndofs
             )
{

    using     Elem = typename MeshT::Elem; 
    constexpr auto Dim=MeshT::Dim;
    constexpr auto ManifoldDim=MeshT::ManifoldDim;

    constexpr auto dofs_per_elem=DofsPerElemNums<Elem,BaseFunctionSpace,BaseFunctionSpaces...>::value;


    // this shuld be passed from outside
    FlagTuple6<Elem,BaseFunctionSpace,BaseFunctionSpaces...> flag_entities;
    auto n_levels=entities.n_levels();
    auto mesh=entities.mesh();
    auto& space_dofs=dofsdm_.space_dofs();
    auto& dofmap_vec=dofsdm_.dofmap();
    std::cout<<"dofs dm reserve begin"<<std::endl;
    dofsdm_.init_dofmap(mesh.n_elements());
    std::cout<<"dofs dm reserve end"<<std::endl;
    // compute the connection node to elem (auxiliary tool, often used)
    clock_t begin = clock();
    
    auto& entities_collection=entities.tuple_entities();
    init<  BaseFunctionSpace,BaseFunctionSpaces...>(entities);
//     NodeToElem<MeshT> node_2_elem(mesh);
//     const auto& node2elem=node_2_elem.val();
//     const auto& n_elements=mesh.n_elements();
    const auto dofs_offset=OffsetDofs<0,Elem,BaseFunctionSpace,BaseFunctionSpaces...>();
//     const auto n_spaces=1+sizeof...(BaseFunctionSpaces); 
//           auto entitiestuple=EntitiesOfFunctionSpaceTuple<Elem,BaseFunctionSpace,BaseFunctionSpaces...>(mesh,node2elem);
          // auto flagtuples= FlagTuple5<Elem,BaseFunctionSpace,BaseFunctionSpaces...>(entities);
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout<<"dofmap entities tuple elapsed_secs="<<elapsed_secs<<std::endl;
    // resize_tuple_of_ptr_vector(dofmap_vec,n_elements);   
    std::vector<Integer> global_dof_count(n_levels,0);
    // loop on all the elements
    // for(Integer space_id=0;space_id<n_spaces;space_id++)
    // {
    //     space_dofs[space_id]=std::make_shared<std::vector<Integer>>();
    //     // space_dofs[space_id]->resize(space_components[space_id][3]);
    // }
    for(std::size_t i=0;i<array_ndofs.size();i++)
      array_ndofs[i]=0;

    auto bisection=entities.bisection();
    // auto level=entities.level();
    begin = clock();




    Integer level=0;
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {


      // if(!elem_belongs_to_level(mesh,elem_iter,level,bisection)) continue;

       // change it for smarter algorithms   
      std::cout<<" dofmap elem == "<< elem_iter <<std::endl;
       auto& elem=mesh.elem(elem_iter);
       auto &elem_id=elem_iter;

       Integer loc_dof_count=0;
       ElementDofMap5<0,Elem,BaseFunctionSpace,BaseFunctionSpaces...>
                         (entities.tuple_entities(),flag_entities,elem,elem_id,dofmap_vec,global_dof_count[level],loc_dof_count,dofs_offset,space_dofs,array_ndofs);

 
      // }

    }
    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout<<"dofmap tuple elapsed_secs="<<elapsed_secs<<std::endl;

    for(std::size_t i=0; i<dofsdm_.n_dofs().size() ;i++)
      dofsdm_.n_dofs()[i]=global_dof_count[level];


// tuple2array_dofs_offset<n_spaces-1>(dofs_offset,dofs_offset_arr);
}


}


#endif