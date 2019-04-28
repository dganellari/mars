#ifndef MARS_CONNECTIVITY_IMPL_HPP
#define MARS_CONNECTIVITY_IMPL_HPP

#include "mars_connectivity.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// Connectivity<Integer Dim, Integer ManifoldDim, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>             ////////
////// Given a mesh, with elements whose dimension is ManifoldDim, defined in a Dim-dimensional space                                   ////////
////// We want to define the connectivity between different entities                                                                    ////////
////// Given an entity e1 of dimension EntityDimFrom, we want to find all the other entities e2 of dimension EntityDimTo                ////////
////// Such that e2 share at least SubEntityDimFrom+1 points with EntityDimFrom                                                         ////////
////// Example 1) Dim=ManifoldDim=3D, e1=triangle, SubEntityDimFrom=0, e2=triangle:                                                     ////////
//////            all the triangles e2 which share a node with the triangle e1                                                          ////////
////// Example 2) Dim=ManifoldDim=3D, e1=triangle, SubEntityDimFrom=1, e2=triangle:                                                     ////////
//////            all the triangles e2 which share an edge with the triangle e1                                                         ////////
////// Example 3) Dim=ManifoldDim=3D, e1=triangle, SubEntityDimFrom=2, e2=triangle:                                                     ////////
//////            all the triangles e2 which share a triangle with the triangle e1 -> e1=e2, known a priori                             ////////
////// Rules: 1) if SubEntityDimFrom > EntityDimFrom: the set is empty                                                                  ////////
//////        2) if SubEntityDimFrom = EntityDimFrom == EntityDimTo: the set is the entity e1 itself                                    ////////
//////        3) The class with EntityDimFrom=SubEntityDimFrom=0 and EntityDimTo=ManifoldDim is defined separately                      ////////
//////           Indeed it will store a vector of vectors, whose component (the node id) returns the neighborhood elements              ////////
//////        4) For all the other classes, given the index entity, the connectivity is computed on the fly                             ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace mars{


// template <Integer Dim, Integer ManifoldDim, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>
// inline std::vector<Integer> Connectivity<Dim,ManifoldDim,EntityDimFrom,SubEntityDimFrom,EntityDimTo>::compute
//                                    (const Entity<Dim,ManifoldDim,EntityDimFrom> &entity_from,
//                                     const Integer index_from,
//                                     const Entity<Dim,ManifoldDim,EntityDimTo>   &entity_to)
// {


//     std::vector<Integer> fromto;
    
//     static_assert(Dim>=0 , " the space dimension must be non negative ");
//     static_assert(ManifoldDim>=0 , " the manifold dimension must be non negative ");
//     static_assert(Dim>=ManifoldDim , " the space dimension must be greater than or equal to than the manifold dimension ");   
//     static_assert(EntityDimFrom>=0 , " the entityfrom dimension must be non negative ");
//     static_assert(ManifoldDim>=EntityDimFrom , " the manifold dimension must be greater than or equal to the entityfrom dimension ");
//     static_assert(SubEntityDimFrom>=0 , " the subentityfrom dimension must be non negative ");
//     static_assert(ManifoldDim>=SubEntityDimFrom , " the manifold dimension must be greater than or equal to the subentityfrom dimension ");   
//     static_assert(EntityDimTo>=0 , " the entityto dimension must be non negative ");
//     static_assert(ManifoldDim>=EntityDimTo , " the manifold dimension must be greater than or equal to the entityto dimension ");
    
//     // rule 1)
//     if(SubEntityDimFrom>EntityDimFrom || SubEntityDimFrom > EntityDimTo)
//          return fromto;
         
//     // rule 2)
//     else if(SubEntityDimFrom==EntityDimFrom && EntityDimFrom==EntityDimTo)
//         {
//          fromto.resize(1);
//          fromto[0]=index_from;
//          return fromto;
//          }
         
//     else
//        {
//         Integer entity_nodes_from[EntityDimFrom+1];
//         Integer entity_nodes_to[EntityDimTo+1];
//         Integer common_nodes;
//         std::vector<Integer> el_neighs;
//         const auto& elem_index_from=entity_from.entity_2_elem(index_from)[0];
//         const auto& iter_entity_from=entity_from.entity_2_elem(index_from)[1];
//         const auto& nodes_from=mesh_.elem(elem_index_from).nodes;
    
//         // entity_nodes_from contains the local indexes of the nodes of the entity
//         Combinations<ManifoldDim + 1, EntityDimFrom+1>::generate(iter_entity_from,entity_nodes_from);

//         // now entity_nodes_from contains the global indexes of the nodes of the entity
//         for(Integer nn=0;nn<EntityDimFrom+1;nn++)
//             entity_nodes_from[nn]=nodes_from[entity_nodes_from[nn]];
        
//         // loop on all the nodes of the entity with index=index_from
//         for(Integer nn=0;nn<EntityDimFrom + 1;nn++)
//             {
//              const Integer& nn_from=entity_nodes_from[nn];
//              // loop on all the elements surrounding
//              for(Integer mm=0;mm<node2elem_[nn_from].size();mm++)
//                  el_neighs.push_back(node2elem_[nn_from][mm]);
//             }
    
//         // define all the neighborhood elements in el_neighs
//         std::sort(el_neighs.begin(),el_neighs.end());
//         el_neighs.erase( std::unique( el_neighs.begin(), el_neighs.end() ), el_neighs.end() );          
    
//         // loop on all the surrounding elements
//         for(Integer elem_iter_to=0;elem_iter_to< el_neighs.size();elem_iter_to++)
//            {
//             const auto &elem_index_to=el_neighs[elem_iter_to];
//             const auto& nodes_to=mesh_.elem(elem_index_to).nodes;
//             // loop on all the entities of dimension EntityDimTo of the element with index=elem_index_to
//             for(Integer iter_entity_to=0;iter_entity_to< entity_to.entity_combinations();iter_entity_to++)
//                 {
//                 // consider the iter_entity nodes
//                 const auto& index_entity_to= entity_to.elem_2_entity(elem_index_to)[iter_entity_to];    
//                 // entity_nodes_to contains the local indexes of the nodes of the simplex   
//                 Combinations<ManifoldDim + 1, EntityDimTo+1>::generate(iter_entity_to,entity_nodes_to);
//                 // now entity_nodes_to contains the global indexes of the nodes of the entity
//                 for(Integer nn=0;nn<EntityDimTo+1;nn++)
//                     entity_nodes_to[nn]=nodes_to[entity_nodes_to[nn]];          
//                 // check if the actual index_entity_to has at least (SubEntityDimFrom + 1) nodes shared with index_entity_from
//                 common_nodes=0;
//                 for(Integer nn=0;nn<EntityDimFrom + 1;nn++)
//                    for(Integer mm=0;mm<EntityDimTo + 1;mm++)
//                       if(entity_nodes_from[nn]==entity_nodes_to[mm]) 
//                          common_nodes++;
                 
//                 if(common_nodes >= SubEntityDimFrom + 1)
//                    fromto.push_back(index_entity_to);
//                 }
//             }
        
//         std::sort(fromto.begin(),fromto.end());
//         fromto.erase( std::unique( fromto.begin(), fromto.end() ), fromto.end() );      

//         return fromto;
//         }
//     };

// template <Integer Dim, Integer ManifoldDim>
// void Connectivity<Dim, ManifoldDim,0,0,ManifoldDim>::init(const Mesh<Dim,ManifoldDim> mesh)
//         {
//             const auto & n_nodes    = mesh.n_nodes();
//             const auto & n_elements = mesh.n_elements();            
//             val_.resize(n_nodes);
//                 // loop on all the nodes of all the active elements 
//                 // add to node_2_element the corresponding active element
//                 for(Integer i = 0; i < n_elements; ++i) 
//                    {
//                         if(!mesh.is_active(i)) continue;

//                         const auto &e    = mesh.elem(i);
//                         // this can be put outside the loop, right?
//                         const auto & nn = ManifoldDim + 1;
//                         for(Integer k = 0; k < nn; ++k) 
//                             {
//                                 val_[e.nodes[k]].push_back(i);
//                             }
//                     }       
//         };



























template <Integer Dim, Integer ManifoldDim, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>
inline std::vector<Integer> 
ElemConnectivity<Simplex<Dim, ManifoldDim>,EntityDimFrom,SubEntityDimFrom,EntityDimTo>::compute
                                          (const ElemEntity<Simplex<Dim,ManifoldDim>,EntityDimFrom> &entity_from,
                                           const Integer index_from,
                                           const ElemEntity<Simplex<Dim,ManifoldDim>,EntityDimTo>   &entity_to)
{


    std::vector<Integer> fromto;
    
    static_assert(Dim>=0 , " the space dimension must be non negative ");
    static_assert(ManifoldDim>=0 , " the manifold dimension must be non negative ");
    static_assert(Dim>=ManifoldDim , " the space dimension must be greater than or equal to than the manifold dimension ");   
    static_assert(EntityDimFrom>=0 , " the entityfrom dimension must be non negative ");
    static_assert(ManifoldDim>=EntityDimFrom , " the manifold dimension must be greater than or equal to the entityfrom dimension ");
    static_assert(SubEntityDimFrom>=0 , " the subentityfrom dimension must be non negative ");
    static_assert(ManifoldDim>=SubEntityDimFrom , " the manifold dimension must be greater than or equal to the subentityfrom dimension ");   
    static_assert(EntityDimTo>=0 , " the entityto dimension must be non negative ");
    static_assert(ManifoldDim>=EntityDimTo , " the manifold dimension must be greater than or equal to the entityto dimension ");
    
    // rule 1)
    if(SubEntityDimFrom>EntityDimFrom || SubEntityDimFrom > EntityDimTo)
         return fromto;
         
    // rule 2)
    else if(SubEntityDimFrom==EntityDimFrom && EntityDimFrom==EntityDimTo)
        {
         fromto.resize(1);
         fromto[0]=index_from;
         return fromto;
         }
         
    else
       {
        Integer entity_nodes_from[EntityDimFrom+1];
        Integer entity_nodes_to[EntityDimTo+1];
        Integer common_nodes;
        std::vector<Integer> el_neighs;
        const auto& elem_index_from=entity_from.entity_2_elem(index_from)[0];
        const auto& iter_entity_from=entity_from.entity_2_elem(index_from)[1];
        const auto& nodes_from=mesh_.elem(elem_index_from).nodes;
    
        // entity_nodes_from contains the local indexes of the nodes of the entity
        Combinations<ManifoldDim + 1, EntityDimFrom+1>::generate(iter_entity_from,entity_nodes_from);

        // now entity_nodes_from contains the global indexes of the nodes of the entity
        for(Integer nn=0;nn<EntityDimFrom+1;nn++)
            entity_nodes_from[nn]=nodes_from[entity_nodes_from[nn]];
        
        // loop on all the nodes of the entity with index=index_from
        for(Integer nn=0;nn<EntityDimFrom + 1;nn++)
            {
             const Integer& nn_from=entity_nodes_from[nn];
             // loop on all the elements surrounding
             for(Integer mm=0;mm<node2elem_[nn_from].size();mm++)
                 el_neighs.push_back(node2elem_[nn_from][mm]);
            }
    
        // define all the neighborhood elements in el_neighs
        std::sort(el_neighs.begin(),el_neighs.end());
        el_neighs.erase( std::unique( el_neighs.begin(), el_neighs.end() ), el_neighs.end() );          
    
        // loop on all the surrounding elements
        for(Integer elem_iter_to=0;elem_iter_to< el_neighs.size();elem_iter_to++)
           {
            const auto &elem_index_to=el_neighs[elem_iter_to];
            const auto& nodes_to=mesh_.elem(elem_index_to).nodes;
            // loop on all the entities of dimension EntityDimTo of the element with index=elem_index_to
            for(Integer iter_entity_to=0;iter_entity_to< entity_to.entity_combinations();iter_entity_to++)
                {
                // consider the iter_entity nodes
                const auto& index_entity_to= entity_to.elem_2_entity(elem_index_to)[iter_entity_to];    
                // entity_nodes_to contains the local indexes of the nodes of the simplex   
                Combinations<ManifoldDim + 1, EntityDimTo+1>::generate(iter_entity_to,entity_nodes_to);
                // now entity_nodes_to contains the global indexes of the nodes of the entity
                for(Integer nn=0;nn<EntityDimTo+1;nn++)
                    entity_nodes_to[nn]=nodes_to[entity_nodes_to[nn]];          
                // check if the actual index_entity_to has at least (SubEntityDimFrom + 1) nodes shared with index_entity_from
                common_nodes=0;
                for(Integer nn=0;nn<EntityDimFrom + 1;nn++)
                   for(Integer mm=0;mm<EntityDimTo + 1;mm++)
                      if(entity_nodes_from[nn]==entity_nodes_to[mm]) 
                         common_nodes++;
                 
                if(common_nodes >= SubEntityDimFrom + 1)
                   fromto.push_back(index_entity_to);
                }
            }
        
        std::sort(fromto.begin(),fromto.end());
        fromto.erase( std::unique( fromto.begin(), fromto.end() ), fromto.end() );      

        return fromto;
        }
    };

template <Integer Dim, Integer ManifoldDim>
void ElemConnectivity<Simplex<Dim, ManifoldDim>,0,0,ManifoldDim>::init(const Mesh<Dim,ManifoldDim> mesh)
        {
            const auto & n_nodes    = mesh.n_nodes();
            const auto & n_elements = mesh.n_elements();            
            val_.resize(n_nodes);
                // loop on all the nodes of all the active elements 
                // add to node_2_element the corresponding active element
                for(Integer i = 0; i < n_elements; ++i) 
                   {
                        if(!mesh.is_active(i)) continue;

                        const auto &e    = mesh.elem(i);
                        // this can be put outside the loop, right?
                        const auto & nn = ManifoldDim + 1;
                        for(Integer k = 0; k < nn; ++k) 
                            {
                                val_[e.nodes[k]].push_back(i);
                            }
                    }       
        };



}

#endif