#ifndef MARS_ENTITY_HPP
#define MARS_ENTITY_HPP

#include "mars_base.hpp"
#include "mars_mesh.hpp"
#include <iostream>
#include <vector>

namespace mars {

template<Integer Dim, Integer ManifoldDim>
class Mesh;
	
template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
class Entity {
public:		
Entity(const Mesh<Dim,ManifoldDim> mesh, const std::vector< std::vector<Integer> >node_2_element); // constructor
std::vector<std::array<Integer,EntityDim> > init_entity(const Mesh<Dim,ManifoldDim> mesh,const std::vector< std::vector<Integer> > node_2_element); 
inline std::vector<std::array<Integer,EntityDim>> value(void) const {return entity_; };	
private:
std::vector<std::array<Integer,EntityDim> > entity_;
 	};



template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
Entity<Dim,ManifoldDim,EntityDim>::Entity(const Mesh<Dim,ManifoldDim> mesh,const std::vector< std::vector<Integer> > node_2_element):
     entity_(init_entity(mesh,node_2_element))
     {};


template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
std::vector<std::array<Integer,EntityDim> > Entity<Dim,ManifoldDim,EntityDim>::init_entity(const Mesh<Dim,ManifoldDim> mesh,const std::vector< std::vector<Integer> > node_2_element)
 {
    const Integer n_elements = mesh.n_elements();
	std::vector<std::array<Integer,EntityDim> > entity;	 
    std::array<Integer, EntityDim> tmp;
    Integer entity_e1[EntityDim];
    Integer entity_e2[EntityDim];
    std::array<Integer, EntityDim> entity_nodes1;
    std::array<Integer, EntityDim> entity_nodes2;
    
    constexpr Integer entity_combinations=Combinations<ManifoldDim+1,EntityDim>::value;


    // loop on all the elements
    for(Integer elem_iter = 0; elem_iter < n_elements; ++elem_iter) {
    
    const auto &el1 = mesh.elem(elem_iter);
    const auto &el_nodes1=el1.nodes;
    std::vector<Integer> el_neighs;
    // find all the elements touching the element=elem_iter
    for(Integer ii=0;ii<ManifoldDim+1;ii++)
       for(Integer jj=0;jj<node_2_element[el_nodes1[ii]].size();jj++)
       {el_neighs.push_back(node_2_element[el_nodes1[ii]][jj]);
       } 
 
    // remove all the elements larger than  elem_iter
    el_neighs.erase(std::remove_if(el_neighs.begin(), el_neighs.end(), [elem_iter](int index_tmp) { return index_tmp>=elem_iter; }), el_neighs.end());       
    // sort and make unique the neighborhood elements
    std::sort( el_neighs.begin(), el_neighs.end() );
    el_neighs.erase( std::unique( el_neighs.begin(), el_neighs.end() ), el_neighs.end() );
    
    // loop on all the entities of the actual element e1
    for(Integer iter_entity_e1=0;iter_entity_e1< entity_combinations;iter_entity_e1++)
    {     
     Combinations<ManifoldDim + 1, EntityDim>::generate(iter_entity_e1,entity_e1);
     
     std::cout<<"(elem1,entity1)==("<<elem_iter<<", "<< iter_entity_e1<<")"<<" "<<std::endl;
     // build the face nodes of iter_entity_e1
     for(int nn=0;nn<EntityDim;nn++)
         {
         entity_nodes1[nn]=el_nodes1[entity_e1[nn]];
         }
         
     std::sort(std::begin(entity_nodes1), std::end(entity_nodes1)); 
     
     bool found_entity=false; 
     
    
       // loop on all the neighborhood elements
       for(Integer elem_iter2 = 0; elem_iter2 < el_neighs.size(); ++elem_iter2) 
    	{	
    	const auto &el2 = mesh.elem(el_neighs[elem_iter2]);
        const auto &el_nodes2=el2.nodes;
        
        // loop on all the entity of the element e2
    	for(Integer iter_entity_e2=0;iter_entity_e2<entity_combinations;iter_entity_e2++)
    	    {
    	    
    	    std::cout<<"(elem2,entity2)==("<<elem_iter2<<", "<< iter_entity_e2<<")"<<" "<<std::endl;
            Combinations<ManifoldDim + 1, EntityDim>::generate(iter_entity_e2,entity_e2);
            

            for(int nn=0;nn<EntityDim;nn++)
                entity_nodes2[nn]=el_nodes2[entity_e2[nn]];   
                	    
   	        std::sort(std::begin(entity_nodes2), std::end(entity_nodes2));  
   	        
   	        if (std::equal(std::begin(entity_nodes1), std::end(entity_nodes1), std::begin(entity_nodes2)))
   	           {
   	           found_entity=true;
   	           break;
   	           }         
    	    }  
	        
	        if(found_entity==true)
	           break;
	        
			}
 	      
   	  if(found_entity==false)
	     {
	     	   std::copy(std::begin(entity_nodes1), std::end(entity_nodes1), std::begin(tmp));
   	           entity.push_back( tmp);
   	           }	  
	  }
	  
	  }
	 return entity;
	 };



     
     
}
#endif //MARS_ENTITY_HPP
