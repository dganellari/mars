#ifndef MARS_FUNCTIONSPACE_DOFMAP_HPP
#define MARS_FUNCTIONSPACE_DOFMAP_HPP


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// 1) EntitiesOfFunctionSpace<Integer Dim, Integer ManifoldDim, Integer SpaceID, Integer Order>:                              ////////                                        ////////
//////    given the function space SpaceID of order Order,                                                                        ////////
//////    it builds a tuple of the entities related to its Dim_dofs                                                                     ////////
////// 2) EntitiesOfFunctionSpaceType is the type of EntitiesOfFunctionSpace                                                            ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mars_base.hpp"
#include "mars_functionspace.hpp"
namespace mars{
















// create a big tuple with N number of "T" types in it                           


template <Integer N, typename T, typename SeqWithArgs>      
struct Append;                                                                  

template <Integer N, typename T, template <typename...> class Seq, typename... Args >      
struct Append<N, T, Seq<Args...> >                                               
{                                                                               
    using type = typename Append<N-1, T, Seq<T,Args...> >::type;                         
};                                                                              

template <typename T, template<typename...> class Seq, typename... Args>              
struct Append<0, T, Seq<Args...> >                                               
{                                                                               
    using type = Seq<Args...>;                                                  
};                                                                              






 



























template<typename FunctionSpace,Integer N=FunctionSpace::entities_nums-1>
class 
FunctionSpaceDofsPerElem
{
 public: 
 static_assert(N>0," FunctionSpaceDofsPerElem N >0");
 static constexpr std::size_t manifold_dim=FunctionSpace::manifold_dim;
 static constexpr std::size_t functionspace_dim=FunctionSpace::functionspace_dim;
 static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
 static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
 static constexpr std::size_t dofs_per_elem=functionspace_dim * 
                                            dofs_per_entity   * 
                                            Combinations<manifold_dim+1,entity_dim+1>::value;
 
 static constexpr std::size_t value=FunctionSpaceDofsPerElem<FunctionSpace,N-1>::value+dofs_per_elem;
};

template<typename FunctionSpace>
class
FunctionSpaceDofsPerElem<FunctionSpace,0>
{
 public:  
 static constexpr std::size_t N=0;
 static constexpr std::size_t manifold_dim=FunctionSpace::manifold_dim;
 static constexpr std::size_t functionspace_dim=FunctionSpace::functionspace_dim;
 static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
 static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
 static constexpr std::size_t value=functionspace_dim * dofs_per_entity*Combinations<manifold_dim+1,entity_dim+1>::value;
};





template<Integer MaxSize,typename FunctionSpace,Integer N,Integer AddConst=0>
typename  std::enable_if< (0<N)&&(N==FunctionSpace::entities_nums-1), void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
{   
    arr[N] = AddConst + FunctionSpaceDofsPerElem<FunctionSpace,N-1>::value;
 }

template<Integer MaxSize,typename FunctionSpace,Integer N,Integer AddConst=0>
typename  std::enable_if< (0<N)&&(N<FunctionSpace::entities_nums-1), void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
{   
    arr[N] = AddConst + FunctionSpaceDofsPerElem<FunctionSpace,N-1>::value;
    FunctionSpaceOffSetDofs<MaxSize,FunctionSpace,N+1,AddConst>(arr);
    
 }



template<Integer MaxSize,typename FunctionSpace,Integer N=0,Integer AddConst=0>
typename  std::enable_if< 0==N && 1<FunctionSpace::entities_nums, void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
{   
    arr[0] = AddConst + 0;
    FunctionSpaceOffSetDofs<MaxSize,FunctionSpace,N+1,AddConst>(arr);
}
template<Integer MaxSize,typename FunctionSpace,Integer N=0,Integer AddConst=0>
typename  std::enable_if< 0==N && 1==FunctionSpace::entities_nums, void>::type
FunctionSpaceOffSetDofs(std::array<Integer , MaxSize> &arr)
{   
    arr[0] = AddConst + 0;
}





template<typename FunctionSpace,typename...FunctionSpaces>
struct TupleTemplateType
{
    using rest = typename TupleTemplateType<FunctionSpaces...>::type;
    using ens = std::array<Integer,FunctionSpace::entities_nums>;
    using tuple_ens=std::tuple<ens>;
    using type = decltype( std::tuple_cat(std::declval< tuple_ens >(), std::declval< rest >() ) );
};

//,typename...
template<typename FunctionSpace>
struct TupleTemplateType<FunctionSpace>
{
 using ens = std::array<Integer,FunctionSpace::entities_nums>;
 using type = typename std::tuple<ens>;
};


// specialization for N=0
template<Integer M=0,typename FunctionSpace,typename...FunctionSpaces>
typename 
std::enable_if<0==sizeof...(FunctionSpaces),
               typename TupleTemplateType<FunctionSpace>::type >::type
TupleTemplate()
{   
    const auto &entities_nums=FunctionSpace::entities_nums;
    std::array<Integer,entities_nums> arr;
        std::cout<<" M=="<<M<<std::endl;

    FunctionSpaceOffSetDofs<entities_nums,FunctionSpace,0,M>(arr);
    return std::tuple<decltype(arr)>(arr);
}

template<Integer M=0,typename FunctionSpace,typename...FunctionSpaces>
typename 
std::enable_if<0<sizeof...(FunctionSpaces),
               typename TupleTemplateType<FunctionSpace,FunctionSpaces...>::type >::type 
TupleTemplate()
{
    std::array<Integer,FunctionSpace::entities_nums> arr;
    FunctionSpaceOffSetDofs<FunctionSpace::entities_nums,FunctionSpace,0,M>(arr);
    std::cout<<" M=="<<M<<std::endl;
    return std::tuple_cat(std::tuple<decltype(arr)>(arr),
                          TupleTemplate<M+FunctionSpaceDofsPerElem<FunctionSpace,FunctionSpace::entities_nums-1>::value,
                                        FunctionSpaces...>());
}















template<typename FunctionSpace,typename...FunctionSpaces>
class DofsPerElemNums
{ 
 public: 
 static constexpr Integer value= (DofsPerElemNums<FunctionSpaces...>::value + FunctionSpaceDofsPerElem<FunctionSpace>::value);

};

template<typename FunctionSpace>
class DofsPerElemNums<FunctionSpace>
{
 public: 
 static constexpr Integer value= FunctionSpaceDofsPerElem<FunctionSpace>::value;
};








template<Integer Dim, Integer ManifoldDim, Integer SpaceID, Integer Order,Integer N>
struct EntitiesOfFunctionSpaceType
{
    using rest = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N-1>::type;
    using ens = Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entity[N]>;
    using tuple_ens=std::tuple<ens>;
    using type = decltype( std::tuple_cat( std::declval< rest >(), std::declval< tuple_ens >() ) );
};


template<Integer Dim, Integer ManifoldDim,Integer SpaceID, Integer Order>
struct EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,0>
{
 using ens = Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entity[0]>;
 using type = typename std::tuple<ens>;
};















// specialization for N=0
template<Integer Dim, Integer ManifoldDim,Integer SpaceID, Integer Order, Integer N=(FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entities_nums-1)>
typename std::enable_if<N==0, 
                        typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type>::type 
EntitiesOfFunctionSpace(const Mesh<Dim, ManifoldDim>&mesh, const std::vector< std::vector<Integer> > &node_2_element)
{    
    using ens=Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entity[N]>;
    return std::tuple<ens>(ens(mesh,node_2_element));
}


template<Integer Dim, Integer ManifoldDim,Integer SpaceID, Integer Order, Integer N=(FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entities_nums-1)>
typename std::enable_if< 0<N, 
                         typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type >::type  
EntitiesOfFunctionSpace(const Mesh<Dim,ManifoldDim>&mesh, const std::vector< std::vector<Integer> >&node_2_element )
{
    using ens=Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entity[N]>;
    return std::tuple_cat(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order,N-1>(mesh,node_2_element),
                          std::tuple<ens>(ens(mesh,node_2_element)));
}
















template<Integer Dim,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces>
struct EntitiesOfFunctionSpaceTupleType
{
 static constexpr Integer SpaceID=FunctionSpace::id;
 static constexpr Integer Order=FunctionSpace::order; 
 static constexpr Integer N=FunctionSpace::entities_nums-1;
 using rest = typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpaces...>::type;
 using ens  = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type;
 using type = decltype( std::tuple_cat( std::declval< rest >(), 
                                        std::declval< std::tuple<ens> >() ) );
};


template<Integer Dim,Integer ManifoldDim,typename FunctionSpace>
struct EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace>
{
 static constexpr Integer SpaceID=FunctionSpace::id;
 static constexpr Integer Order=FunctionSpace::order; 
 static constexpr Integer N=FunctionSpace::entities_nums-1;
 using ens =  typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type ;
 using type = typename std::tuple<ens>;
};




/////////////////////////////////////////////////////////////////////////////



// general case 
template<Integer Dim,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces>
typename std::enable_if< 0<sizeof...(FunctionSpaces), 
                         typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace,FunctionSpaces...>::type>::type
EntitiesOfFunctionSpaceTuple
(const Mesh<Dim, ManifoldDim>&mesh, const std::vector< std::vector<Integer> > &node_2_element)
{    
    constexpr Integer SpaceID=FunctionSpace::id;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer N=FunctionSpace::entities_nums-1;
    using type = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type;
    std::cout<<"mmm1111"<<std::endl;
    //return std::tuple<type>(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node_2_element));
    return std::tuple_cat(EntitiesOfFunctionSpaceTuple<Dim,ManifoldDim,FunctionSpaces...>(mesh,node_2_element),
                          std::tuple<type>(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node_2_element))
                         );
}

// base case
template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces>
 typename std::enable_if< 0==sizeof...(FunctionSpaces), 
                          typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace>::type>::type
//typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace>::type
EntitiesOfFunctionSpaceTuple
(const Mesh<Dim, ManifoldDim>&mesh, 
 const std::vector< std::vector<Integer> > &node_2_element)
{    
    constexpr Integer SpaceID=FunctionSpace::id;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer N=FunctionSpace::entities_nums-1;
    using type = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,SpaceID,Order,N>::type;
    std::cout<<"mmm0000"<<std::endl;
    return std::tuple<type>(EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node_2_element));
}





/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   INITIALIZE BOOL VECTOR ENTITIES   ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


// base case (M==tuple size)
template<Integer Dim, Integer ManifoldDim, Integer SpaceID, Integer Order, Integer M  , typename ...Args>
typename std::enable_if< M==sizeof ...(Args),void >::type
    initialize_vector_entities(std::tuple< Args...> const &tuple, 
                               std::array<std::vector< std::array<Integer,2> >, FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entities_nums > & entity_found)
    {
     static_assert(M>=0," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");
        
    };
// general case (M<tuple size)
template<Integer Dim,Integer ManifoldDim, Integer SpaceID, Integer Order, Integer M = 0, typename...Args>
typename std::enable_if< M<sizeof ...(Args),void >::type
    initialize_vector_entities(std::tuple< Args...> const &tuple,
                               std::array<std::vector< std::array<Integer,2> >, FunctionSpace<Dim,ManifoldDim,SpaceID,Order>::entities_nums > & entity_found)
    {
     static_assert(M>=0," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");

     Integer entity_length=std::get<M>(tuple).size();
     entity_found[M].resize(entity_length,{-1,-1});
     std::cout<<" (M,entity_length )=("<<M<<", "<<entity_length<<")"<<std::endl;
     // iterate to the next 
     initialize_vector_entities<Dim,ManifoldDim,SpaceID,Order,M+1,Args...>(tuple,entity_found);
        
    };








template<Integer Dim,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces>
struct FlagTupleType
{
 static constexpr Integer SpaceID=FunctionSpace::id;
 static constexpr Integer Order=FunctionSpace::order; 
 static constexpr Integer entities_nums=FunctionSpace::entities_nums;

 using rest = typename FlagTupleType<Dim,ManifoldDim,FunctionSpaces...>::type;
 using ens  = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
 using type = decltype( std::tuple_cat( std::declval< rest >(), 
                                        std::declval< std::tuple<ens> >() ) );
};


template<Integer Dim,Integer ManifoldDim,typename FunctionSpace>
struct FlagTupleType<Dim,ManifoldDim,FunctionSpace>
{
 static constexpr Integer SpaceID=FunctionSpace::id;
 static constexpr Integer Order=FunctionSpace::order; 
 static constexpr Integer entities_nums=FunctionSpace::entities_nums;
 using ens =  typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
 using type = typename std::tuple<ens>;
};




template<Integer Dim,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces, typename...Args>
typename std::enable_if< 0<sizeof...(FunctionSpaces), 
                         typename FlagTupleType<Dim,ManifoldDim,FunctionSpace,FunctionSpaces...>::type>::type
FlagTuple
(std::tuple<Args...> tuple)
{    
    constexpr Integer SpaceID=FunctionSpace::id;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer M=sizeof...(FunctionSpaces);
    static constexpr Integer entities_nums=FunctionSpace::entities_nums;
    using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
    type ens;

    initialize_vector_entities<Dim,ManifoldDim,SpaceID,Order>(std::get<M>(tuple),ens);
    std::cout<<"aaa1111"<<std::endl;
    return std::tuple_cat(FlagTuple<Dim,ManifoldDim,FunctionSpaces...>(tuple),
                          std::tuple<type>(ens)
                         );
}

// base case
template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces, typename...Args>
 typename std::enable_if< 0==sizeof...(FunctionSpaces), 
                          typename FlagTupleType<Dim,ManifoldDim,FunctionSpace>::type>::type
//typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace>::type
FlagTuple(std::tuple<Args...> tuple)
{    
    constexpr Integer SpaceID=FunctionSpace::id;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer entities_nums=FunctionSpace::entities_nums;

    using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
    type ens;
    initialize_vector_entities<Dim,ManifoldDim,SpaceID,Order>(std::get<0>(tuple),ens);
    std::cout<<"aaa0000"<<std::endl;
    return std::tuple<type>(ens);
}

































template<Integer ManifoldDim, Integer N,typename FunctionSpace,typename Array,typename...Args1,typename T,typename OS>
typename std::enable_if< -1<N,void>::type
Loop_Entities(const std::tuple<Args1...>& entitiestuple,
              Array& flagtuples,
              const Integer& elem_id, 
              std::vector<T> &dofmap_vec,
              Integer& global_dof_count,
              Integer& loc_dof_count,
              OS &dofs_offset
              )
{

Loop_Entities<ManifoldDim,N-1,FunctionSpace>
             (entitiestuple,flagtuples,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset);

auto&entity=std::get<N>(entitiestuple);
auto&flag=flagtuples[N];
// I have to loop on all entities of the given type


constexpr Integer functionspace_dim=FunctionSpace::functionspace_dim;
constexpr Integer entity_dim=FunctionSpace::entity[N];
constexpr Integer dofs_per_entity=FunctionSpace::dofs_per_entity[N];
constexpr Integer entity_points=entity_dim+1;
constexpr Integer manifold_points=ManifoldDim+1;
constexpr Integer combinations_nums=Combinations<manifold_points,entity_points>::value;

const auto& elem2entity=entity.elem_2_entity(elem_id);

std::cout<<std::endl;
std::cout<<"SPACE N=="<<N<<", "<<FunctionSpace::entity[N]<<std::endl;
std::cout<<std::endl;

    std::cout<<"dofs_offset =="<<dofs_offset[N]<<", "<<std::endl;
std::cout<<std::endl;

// loop on all the entities of a given entity_dim
for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
    {
     // check if the entity has been already visited
     const auto& entity_id=elem2entity[entity_iter];
     // std::cout<<"sono molto="<< elem_id<<", "<<entity_iter<<", "<<flag[entity_id][0]<<std::endl;
    // if not, then createa functionspace_dim new dofs
    if(flag[entity_id][0]==-1)
    {
       flag[entity_id][0]=elem_id;
       flag[entity_id][1]=entity_iter;
      for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
       for(Integer fs_dim=0;fs_dim<functionspace_dim;fs_dim++)
       {
        
        dofmap_vec[elem_id][loc_dof_count]=global_dof_count;
      if(N==1) 
      std::cout<<"loc_dof_count="<<loc_dof_count<<" global_dof_count="<<global_dof_count<<std::endl;
        
        loc_dof_count++;
        global_dof_count++;

       }


    }
    else
    {
    // if the entity has been already visited
    // find the element (elem_id_tmp) in which was found 
    // and the corresponding entity_iter (iter_tmp)
      const auto& elem_id_tmp=flag[entity_id][0];
      const auto& iter_tmp=flag[entity_id][1];
       auto old_dof=dofmap_vec[elem_id_tmp][dofs_offset[N]+ dofs_per_entity*iter_tmp*FunctionSpace::functionspace_dim];
      //std::cout<<"QUI "<<std::endl;
      //std::cout<<"elem_id="<<elem_id<<", elem_id_tmp="<<elem_id_tmp <<std::endl;
      //std::cout<<"entity_id="<<entity_id<<", iter_tmp="<<iter_tmp <<std::endl; 


      //dofmap_vec[elem_id_tmp][loc_dof_count]=old_dof;


    // then find in dofmap_vec the dof of the entity 
     // const auto& entity_ids=entity.elem_2_entity(elem_id_tmp);
      
     // Integer entity_iter=0;
      // for(entity_iter=0;entity_iter<entity_ids.size();entity_iter++)
      //    if(entity_ids[entity_iter]==entity_id)
      //    {
      //     break;  
      //    }
      // qui devi fare offset[N]+
    //  dofmap_vec[elem_id_tmp][entity_iter] 
      // std::cout<<"ci entro="<<elem_id_tmp<<std::endl;

      for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
       {
        for(Integer fs_dim=0;fs_dim<functionspace_dim;fs_dim++)
       {
        
        dofmap_vec[elem_id][loc_dof_count]=old_dof;
      if(N==1) 
        {std::cout<<" entity_id="<<entity_id<<"     loc_dof_count="<<loc_dof_count<<"     old dof="<<old_dof<<"   dofs_offset="<< dofs_offset[N]<<"    iter_tmp="<< iter_tmp<<"   functionspace_dim="<<FunctionSpace::functionspace_dim<<std::endl;
        }

        loc_dof_count++;
        old_dof++;
       }
   }



    }

    }
// do stuff
// for(Integer elem_iter=0;elem_iter<flag.size();elem_iter++)
//     std::cout<<"ooooooooooooo---------->>>"<<flag[elem_iter]<<std::endl;
    std::cout<<std::endl;
    //  for(Integer elem_iter=0;elem_iter<dofmap_vec.size();elem_iter++)
    // {std::cout<<"inside Loop_Entities"<<std::endl;
    //  auto &elem_id=elem_iter;
    //  for(Integer nn=0;nn<dofmap_vec[0].size();nn++)
    //  {
    //     std::cout<<dofmap_vec[elem_id][nn]<<" ";
    //  }
    // } 
}

template<Integer ManifoldDim,Integer N,typename FunctionSpace, typename Array,typename...Args1,typename T,typename OS>
typename std::enable_if< -1==N,void>::type
Loop_Entities(const std::tuple<Args1...>& entitiestuple,
              Array& flagtuples,
              const Integer& elem_id,
              std::vector<T> & dofmap_vec,
              Integer& global_dof_count,
              Integer &loc_dof_count,
              OS &dofs_offset  )
{}



template<Integer K=0,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces, typename...Args1,typename...Args2,typename T,typename OS>
typename std::enable_if< 0<sizeof...(FunctionSpaces),void>::type
Loop_functionspace(const std::tuple<Args1...>& entitiestuple,
                   std::tuple<Args2...>& flagtuples,
                   const Integer& elem_id,
                   std::vector<T>& dofmap_vec,
                   Integer& global_dof_count,
                   Integer& loc_dof_count,
                   OS &dofs_offset   )
{
 static constexpr Integer M=sizeof...(FunctionSpaces);
 const auto m1=std::get<M>(entitiestuple);
 auto &m2=std::get<M>(flagtuples);
 static constexpr Integer NN=std::tuple_size<decltype(m1)>::value-1;
 std::cout<<"Loop_functionspace M=== "<<M<<", "<<FunctionSpace::id<<", "<<FunctionSpace::order<< std::endl;
 Loop_Entities<ManifoldDim,NN,FunctionSpace>(m1,m2,elem_id,dofmap_vec,global_dof_count,loc_dof_count,std::get<K>(dofs_offset));
 Loop_functionspace<K+1,ManifoldDim,FunctionSpaces...>(entitiestuple,flagtuples,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset);
};




template<Integer K=0,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces, typename...Args1,typename...Args2,typename T,typename OS>
typename std::enable_if< 0==sizeof...(FunctionSpaces),void>::type
Loop_functionspace(const std::tuple<Args1...>& entitiestuple,
                   std::tuple<Args2...>& flagtuples,
                   const Integer& elem_id,
                   std::vector<T>& dofmap_vec,
                   Integer& global_dof_count, 
                   Integer& loc_dof_count,
                   OS &dofs_offset )
{
 const auto m1=std::get<0>(entitiestuple);
 auto& m2=std::get<0>(flagtuples);
 static constexpr Integer NN=std::tuple_size<decltype(m1)>::value-1;
 std::cout<<"Loop_functionspace M=== "<<0<<", "<<FunctionSpace::id<<", "<<FunctionSpace::order<< std::endl;
 Loop_Entities<ManifoldDim,NN,FunctionSpace>(m1,m2,elem_id,dofmap_vec,global_dof_count,loc_dof_count,std::get<K>(dofs_offset));
};




























template <typename ... Ts>
constexpr std::tuple<Ts...> tuple_transform(std::tuple<Ts...>)
{
    return std::make_tuple(Ts{}()...);
}


template<typename TupleComponent, typename...TupleComponents>
struct TupleType
{
    using rest = typename TupleType<TupleComponents...>::type;
    using tuple_component=std::tuple<TupleComponent>;
    using type = decltype( std::tuple_cat( std::declval< rest >(), std::declval< tuple_component >() ) );
};

template<typename TupleComponent>
struct TupleType<TupleComponent>
{
    using type = std::tuple< TupleComponent >;
};




// specialization for N=0
template<typename TupleComponent, typename...TupleComponents>
constexpr typename std::enable_if< 0==sizeof...(TupleComponents), 
                        typename TupleType<TupleComponent>::type>::type 
TupleCreate(TupleComponent tuplecomponent,TupleComponents...tuplecomponents)
{    
    return std::tuple<TupleComponent>(tuplecomponent);//params...));
}


template<typename TupleComponent, typename...TupleComponents>
constexpr typename std::enable_if< 0<sizeof...(TupleComponents), 
                         typename TupleType<TupleComponent,TupleComponents...>::type>::type   
TupleCreate(TupleComponent tuplecomponent,TupleComponents...tuplecomponents)
{
    return std::tuple_cat(std::tuple<TupleComponent>(tuplecomponent),
                          TupleCreate<TupleComponents...>(tuplecomponents...)//(params...),
                          );//params...)));
}










template<unsigned... args> struct ArrayHolder {
    static constexpr std::array<std::array<Integer,2 >,sizeof...(args)> data={args...};
};

template<unsigned... args> 
constexpr std::array<std::array<Integer,2 >,sizeof...(args)> ArrayHolder<args...>::data;




template<Integer N=0, Integer M,Integer P,typename T,typename...Ts>
typename  std::enable_if< 0==sizeof...(Ts), void>::type
ArrOfArr(std::array<std::array<Integer, P> , M > &arr)
{   
    arr[N] = std::array<Integer, P>{T::space_dim, T::manifold_dim, T::id, T::order, T::functionspace_dim};
 
}

template<Integer N=0,Integer M,Integer P,typename T,typename...Ts>
typename  std::enable_if< 0<sizeof...(Ts), void>::type
ArrOfArr(std::array<std::array<Integer, P> , M> &arr)
{   
    arr[N] = std::array<Integer, P>{T::space_dim, T::manifold_dim, T::id, T::order, T::functionspace_dim};
    ArrOfArr<N+1,M,P,Ts...>(arr);
    
 }

// template< typename T, typename...Ts, Integer N=sizeof...(Ts)>
// constexpr std::array<std::array<Integer, 5> , sizeof...(Ts)+1 > fillarr()
// {   
//       constexpr std::array<std::array<Integer, 5> , sizeof...(Ts)+1 > array;
//     std::array[N] = std::array<Integer, 4>{T::space_dim, T::manifold_dim, T::id, T::order, T::functionspace_dim};

//     return 
// }

template <typename ... Ts>
struct MixedFunctionSpace
 {

    const std::tuple<decltype(Ts::dofs_per_entity) & ...> tpl{ Ts::dofs_per_entity... };

    static constexpr std::array<std::size_t, sizeof...(Ts)> arr_entities_nums
    {{ sizeof(Ts::dofs_per_entity)/sizeof(Ts::dofs_per_entity[0])... }};

   template <std::size_t N>
   const decltype(std::get<N>(tpl)) dofs_per_entity ()
    { return std::get<N>(tpl); }

  // template <std::size_t N>
  //  static constexpr decltype(std::get<N>(tpl)) dofs_per_entity2 =std::get<N>(tpl); 


   static constexpr std::size_t entities_nums (std::size_t n)
    { return arr_entities_nums[n]; }
   // template<Integer n>
   // static constexpr std::size_t entities_nums=arr_entities_nums[n]; 



 };

template <typename ... Ts>
constexpr std::array<std::size_t, sizeof...(Ts)> MixedFunctionSpace<Ts...>::arr_entities_nums;

// template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces,Integer N=sizeof...(FunctionSpaces)+1>
// typename std::enable_if<0<N,void >::type
//  BuildMixedFunctionSpace(const Mesh<Dim,ManifoldDim>&mesh,std::array<FunctionSpace, 1+sizeof...(FunctionSpaces)> spaces)
//  {

//   spaces[N]=
//  }

// template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces,Integer N=sizeof...(FunctionSpaces)+1>
// typename std::enable_if<0==N,void >::type
//  BuildMixedFunctionSpace(const Mesh<Dim,ManifoldDim>&mesh,std::tuple<FunctionSpace, FunctionSpaces...> spaces)
//  {
//   spaces[N]=FunctionSpace;
//  }

// template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces>
// class MixedFunctionSpace{
// public:
// MixedFunctionSpace(const Mesh<Dim,ManifoldDim>&mesh):
// spaces(BuildMixedFunctionSpace<Dim,ManifoldDim,FunctionSpace,FunctionSpaces...>(mesh))
// {};
// private:
// std::tuple<FunctionSpace, FunctionSpaces...> spaces;
// };
















template<Integer Dim,Integer ManifoldDim, typename...Functions>  //Integer Dim, Integer ManifoldDim,Integer SpaceID,Integer Order>//, typename...FunctionSpace>
typename std::enable_if<0==sizeof...(Functions), void>::type
dofmap(const Mesh<Dim,ManifoldDim>&mesh){};


template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces>
void dofmap(const Mesh<Dim,ManifoldDim>&mesh)
{
    const auto&n_elements=mesh.n_elements();
    // compute the connection node to elem (auxiliary tool, often used)
    NodeToElem<Dim,ManifoldDim> node_2_elem(mesh);
    const std::vector< std::vector<Integer> > & node2elem=node_2_elem.val();
    // compute 
    // constexpr Integer SpaceID=FunctionSpace::id;
    // constexpr Integer Order=FunctionSpace::order;
    // //auto ehi3=EntitiesOfFunctionSpace<Dim,ManifoldDim,SpaceID,Order>(mesh,node2elem);
    // constexpr Integer NN=FunctionSpace::entities_nums-1;
  
    auto entitiestuple=EntitiesOfFunctionSpaceTuple<Dim,ManifoldDim,FunctionSpace,FunctionSpaces...>(mesh,node2elem);
//    auto flagtuples= FlagTuple<Dim,ManifoldDim,FunctionSpace>(entitiestuple);
    auto flagtuples= FlagTuple<Dim,ManifoldDim,FunctionSpace,FunctionSpaces...>(entitiestuple);

    const auto dofs_offset=TupleTemplate<0,FunctionSpace,FunctionSpaces...>();

    constexpr std::size_t dofs_per_elem=DofsPerElemNums<FunctionSpace,FunctionSpaces...>::value;
    std::vector<std::array<Integer,dofs_per_elem>> dofmap_vec(n_elements);

    Integer global_dof_count=0;
    //auto second_entity=std::get<1>(entitiestuple);
    // loop on all the elements
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {
     // change it for smarter algorithms   
     auto &elem_id=elem_iter;
     Integer loc_dof_count=0;
     // loop on all the function spaces
     Loop_functionspace<0,ManifoldDim,FunctionSpace,FunctionSpaces...>
                       (entitiestuple,flagtuples,elem_id,dofmap_vec,global_dof_count,loc_dof_count,dofs_offset);
    }

std::cout<<"ORA"<<std::endl;


     for(Integer nn2=0;nn2<std::get<0>(dofs_offset).size();nn2++)
           std::cout<<std::get<0>(dofs_offset)[nn2]<<" ";
     std::cout<<std::endl;
     for(Integer nn2=0;nn2<std::get<1>(dofs_offset).size();nn2++)
           std::cout<<std::get<1>(dofs_offset)[nn2]<<" ";
     std::cout<<std::endl;
     // for(Integer nn2=0;nn2<std::get<2>(dofs_offset).size();nn2++)
     //       std::cout<<std::get<2>(dofs_offset)[nn2]<<" ";
     // std::cout<<std::endl;

   
std::cout<<std::endl;
    for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
    {
     auto &elem_id=elem_iter;
     std::cout<<"elem_id="<<elem_id<<", "<<dofmap_vec[elem_id].size()<<std::endl;
     for(Integer nn=0;nn<dofmap_vec[0].size();nn++)
     {
        std::cout<<dofmap_vec[elem_id][nn]<<" ";
     }
     std::cout<<std::endl;
    } 



};







}


#endif