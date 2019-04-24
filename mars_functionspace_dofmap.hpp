#ifndef MARS_FUNCTIONSPACE_DOFMAP_HPP
#define MARS_FUNCTIONSPACE_DOFMAP_HPP


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// 1) EntitiesOfFunctionSpace<Integer Dim, Integer ManifoldDim, Integer FEFamily, Integer Order>:                              ////////                                        ////////
//////    given the function space FEFamily of order Order,                                                                        ////////
//////    it builds a tuple of the entities related to its Dim_dofs                                                                     ////////
////// 2) EntitiesOfFunctionSpaceType is the type of EntitiesOfFunctionSpace                                                            ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mars_base.hpp"
#include "mars_functionspace.hpp"
namespace mars{



template<typename FunctionSpace,Integer N=FunctionSpace::entities_nums-1>
class 
FunctionSpaceDofsPerElem
{
 public: 
 static_assert(N>0," FunctionSpaceDofsPerElem N >0");
 static constexpr std::size_t manifold_dim=FunctionSpace::manifold_dim;
 static constexpr std::size_t n_components=FunctionSpace::n_components;
 static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
 static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
 static constexpr std::size_t dofs_per_elem=n_components * 
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
 static constexpr std::size_t n_components=FunctionSpace::n_components;
 static constexpr std::size_t entity_dim=FunctionSpace::entity[N];
 static constexpr std::size_t dofs_per_entity=FunctionSpace::dofs_per_entity[N];
 static constexpr std::size_t value=n_components * dofs_per_entity*Combinations<manifold_dim+1,entity_dim+1>::value;
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








template<Integer Dim, Integer ManifoldDim, Integer FEFamily, Integer Order,Integer N>
struct EntitiesOfFunctionSpaceType
{
    using rest = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,FEFamily,Order,N-1>::type;
    using ens = Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,FEFamily,Order>::entity[N]>;
    using tuple_ens=std::tuple<ens>;
    using type = decltype( std::tuple_cat( std::declval< rest >(), std::declval< tuple_ens >() ) );
};


template<Integer Dim, Integer ManifoldDim,Integer FEFamily, Integer Order>
struct EntitiesOfFunctionSpaceType<Dim,ManifoldDim,FEFamily,Order,0>
{
 using ens = Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,FEFamily,Order>::entity[0]>;
 using type = typename std::tuple<ens>;
};















// specialization for N=0
template<Integer Dim, Integer ManifoldDim,Integer FEFamily, Integer Order, Integer N=(FunctionSpace<Dim,ManifoldDim,FEFamily,Order>::entities_nums-1)>
typename std::enable_if<N==0, 
                        typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,FEFamily,Order,N>::type>::type 
EntitiesOfFunctionSpace(const Mesh<Dim, ManifoldDim>&mesh, const std::vector< std::vector<Integer> > &node_2_element)
{    
    using ens=Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,FEFamily,Order>::entity[N]>;
    return std::tuple<ens>(ens(mesh,node_2_element));
}


template<Integer Dim, Integer ManifoldDim,Integer FEFamily, Integer Order, Integer N=(FunctionSpace<Dim,ManifoldDim,FEFamily,Order>::entities_nums-1)>
typename std::enable_if< 0<N, 
                         typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,FEFamily,Order,N>::type >::type  
EntitiesOfFunctionSpace(const Mesh<Dim,ManifoldDim>&mesh, const std::vector< std::vector<Integer> >&node_2_element )
{
    using ens=Entity<Dim,ManifoldDim,FunctionSpace<Dim,ManifoldDim,FEFamily,Order>::entity[N]>;
    return std::tuple_cat(EntitiesOfFunctionSpace<Dim,ManifoldDim,FEFamily,Order,N-1>(mesh,node_2_element),
                          std::tuple<ens>(ens(mesh,node_2_element)));
}
















template<Integer Dim,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces>
struct EntitiesOfFunctionSpaceTupleType
{
 static constexpr Integer FEFamily=FunctionSpace::family;
 static constexpr Integer Order=FunctionSpace::order; 
 static constexpr Integer N=FunctionSpace::entities_nums-1;
 using rest = typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpaces...>::type;
 using ens  = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,FEFamily,Order,N>::type;
 using type = decltype( std::tuple_cat( std::declval< rest >(), 
                                        std::declval< std::tuple<ens> >() ) );
};


template<Integer Dim,Integer ManifoldDim,typename FunctionSpace>
struct EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace>
{
 static constexpr Integer FEFamily=FunctionSpace::family;
 static constexpr Integer Order=FunctionSpace::order; 
 static constexpr Integer N=FunctionSpace::entities_nums-1;
 using ens =  typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,FEFamily,Order,N>::type ;
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
    constexpr Integer FEFamily=FunctionSpace::family;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer N=FunctionSpace::entities_nums-1;
    using type = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,FEFamily,Order,N>::type;
    return std::tuple_cat(EntitiesOfFunctionSpaceTuple<Dim,ManifoldDim,FunctionSpaces...>(mesh,node_2_element),
                          std::tuple<type>(EntitiesOfFunctionSpace<Dim,ManifoldDim,FEFamily,Order>(mesh,node_2_element))
                         );
}

// base case
template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces>
 typename std::enable_if< 0==sizeof...(FunctionSpaces), 
                          typename EntitiesOfFunctionSpaceTupleType<Dim,ManifoldDim,FunctionSpace>::type>::type
EntitiesOfFunctionSpaceTuple
(const Mesh<Dim, ManifoldDim>&mesh, 
 const std::vector< std::vector<Integer> > &node_2_element)
{    
    constexpr Integer FEFamily=FunctionSpace::family;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer N=FunctionSpace::entities_nums-1;
    using type = typename EntitiesOfFunctionSpaceType<Dim,ManifoldDim,FEFamily,Order,N>::type;
    return std::tuple<type>(EntitiesOfFunctionSpace<Dim,ManifoldDim,FEFamily,Order>(mesh,node_2_element));
}





/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   INITIALIZE BOOL VECTOR ENTITIES   ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


// base case (M==tuple size)
template<Integer Dim, Integer ManifoldDim, Integer FEFamily, Integer Order, Integer M  , typename ...Args>
typename std::enable_if< M==sizeof ...(Args),void >::type
    initialize_vector_entities(std::tuple< Args...> const &tuple, 
                               std::array<std::vector< std::array<Integer,2> >, FunctionSpace<Dim,ManifoldDim,FEFamily,Order>::entities_nums > & entity_found)
    {
     static_assert(M>=0," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");
        
    };
// general case (M<tuple size)
template<Integer Dim,Integer ManifoldDim, Integer FEFamily, Integer Order, Integer M = 0, typename...Args>
typename std::enable_if< M<sizeof ...(Args),void >::type
    initialize_vector_entities(std::tuple< Args...> const &tuple,
                               std::array<std::vector< std::array<Integer,2> >, FunctionSpace<Dim,ManifoldDim,FEFamily,Order>::entities_nums > & entity_found)
    {
     static_assert(M>=0," the tuple must have non negative length ");
     static_assert(ManifoldDim > 0," the ManifoldDim must be positive ");

     Integer entity_length=std::get<M>(tuple).size();
     entity_found[M].resize(entity_length,{-1,-1});
     // iterate to the next 
     initialize_vector_entities<Dim,ManifoldDim,FEFamily,Order,M+1,Args...>(tuple,entity_found);
        
    };








template<Integer Dim,Integer ManifoldDim,typename FunctionSpace,typename...FunctionSpaces>
struct FlagTupleType
{
 static constexpr Integer FEFamily=FunctionSpace::family;
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
 static constexpr Integer FEFamily=FunctionSpace::family;
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
    constexpr Integer FEFamily=FunctionSpace::family;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer M=sizeof...(FunctionSpaces);
    static constexpr Integer entities_nums=FunctionSpace::entities_nums;
    using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
    type ens;

    initialize_vector_entities<Dim,ManifoldDim,FEFamily,Order>(std::get<M>(tuple),ens);
    return std::tuple_cat(FlagTuple<Dim,ManifoldDim,FunctionSpaces...>(tuple),
                          std::tuple<type>(ens)
                         );
}

// base case
template<Integer Dim,Integer ManifoldDim,typename FunctionSpace, typename...FunctionSpaces, typename...Args>
 typename std::enable_if< 0==sizeof...(FunctionSpaces), 
                          typename FlagTupleType<Dim,ManifoldDim,FunctionSpace>::type>::type
FlagTuple(std::tuple<Args...> tuple)
{    
    constexpr Integer FEFamily=FunctionSpace::family;
    constexpr Integer Order=FunctionSpace::order;
    static constexpr Integer entities_nums=FunctionSpace::entities_nums;

    using type = typename std::array<std::vector<std::array<Integer,2>>, entities_nums>;
    type ens;
    initialize_vector_entities<Dim,ManifoldDim,FEFamily,Order>(std::get<0>(tuple),ens);
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

constexpr Integer continuity=FunctionSpace::continuity;
constexpr Integer n_components=FunctionSpace::n_components;
constexpr Integer entity_dim=FunctionSpace::entity[N];
constexpr Integer dofs_per_entity=FunctionSpace::dofs_per_entity[N];
constexpr Integer entity_points=entity_dim+1;
constexpr Integer manifold_points=ManifoldDim+1;
constexpr Integer combinations_nums=Combinations<manifold_points,entity_points>::value;

const auto& elem2entity=entity.elem_2_entity(elem_id);

// loop on all the entities of a given entity_dim
for(Integer entity_iter=0;entity_iter<combinations_nums;entity_iter++)
    {
     // check if the entity has been already visited
     const auto& entity_id=elem2entity[entity_iter];
    // if not, then createa n_components new dofs
    if(flag[entity_id][0]==-1 || continuity==false)
    {
       flag[entity_id][0]=elem_id;
       flag[entity_id][1]=entity_iter;
      for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
       for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
       {
        
        dofmap_vec[elem_id][loc_dof_count]=global_dof_count;
        
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
       auto old_dof=dofmap_vec[elem_id_tmp][dofs_offset[N]+ dofs_per_entity*iter_tmp*FunctionSpace::n_components];
   
      for(Integer entity_dofs_iter=0;entity_dofs_iter<dofs_per_entity;entity_dofs_iter++)
       {
        for(Integer fs_dim=0;fs_dim<n_components;fs_dim++)
       {
        
        dofmap_vec[elem_id][loc_dof_count]=old_dof;

        loc_dof_count++;
        old_dof++;
       }
   }



    }

    }

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
 Loop_Entities<ManifoldDim,NN,FunctionSpace>(m1,m2,elem_id,dofmap_vec,global_dof_count,loc_dof_count,std::get<K>(dofs_offset));
};


















































template<Integer Dim,Integer ManifoldDim, typename...Functions>  //Integer Dim, Integer ManifoldDim,Integer FEFamily,Integer Order>//, typename...FunctionSpace>
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
    // constexpr Integer FEFamily=FunctionSpace::id;
    // constexpr Integer Order=FunctionSpace::order;
    // //auto ehi3=EntitiesOfFunctionSpace<Dim,ManifoldDim,FEFamily,Order>(mesh,node2elem);
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