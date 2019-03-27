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

	inline const Integer entity_combinations(void) const {return entity_combinations_;};

	//std::vector<std::array<Integer,EntityDim> > 
	void init_entity(const Mesh<Dim,ManifoldDim> mesh,const std::vector< std::vector<Integer> > node_2_element, 
	                std::vector<std::array<Integer,EntityDim> > &entity_, std::vector<std::vector<Integer> > &elem_2_entity_); 

	inline std::vector<std::array<Integer,EntityDim> > value(void) const {return entity_; };
	inline std::vector<std::vector<Integer> > elem_2_entity(void) const {return elem_2_entity_; };	


private:
	Integer entity_combinations_; // number of entities with EntityDim points inside the simplex of dimension ManifoldDim
	std::vector<std::array<Integer,EntityDim> > entity_;
	std::vector<std::vector<Integer> > elem_2_entity_;

 	};



// template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
// Entity<Dim,ManifoldDim,EntityDim>::Entity(const Mesh<Dim,ManifoldDim> mesh,const std::vector< std::vector<Integer> > node_2_element):
//      entity_combinations_(Combinations<ManifoldDim+1,EntityDim>::value),
//      entity_(init_entity(mesh,node_2_element))
//      {};


template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
Entity<Dim,ManifoldDim,EntityDim>::Entity(const Mesh<Dim,ManifoldDim> mesh,const std::vector< std::vector<Integer> > node_2_element)
     {entity_combinations_=Combinations<ManifoldDim+1,EntityDim>::value;
     init_entity(mesh,node_2_element,entity_,elem_2_entity_); };
     
     
template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
//std::vector<std::array<Integer,EntityDim> > 
void Entity<Dim,ManifoldDim,EntityDim>::init_entity(const Mesh<Dim,ManifoldDim> mesh,const std::vector< std::vector<Integer> > node_2_element,
                                    std::vector<std::array<Integer,EntityDim> > &entity_, std::vector<std::vector<Integer> > &elem_2_entity_)
 {
    const Integer n_elements = mesh.n_elements();
    Integer entity_number=0;
    std::array<Integer, EntityDim> tmp;
    Integer entity_e1[EntityDim];
    Integer entity_e2[EntityDim];
    std::array<Integer, EntityDim> entity_nodes1;
    std::array<Integer, EntityDim> entity_nodes2;
    //std::vector<std::vector<Integer> > 
    //elem_2_entity_(n_elements);
    elem_2_entity_.resize(n_elements);
    // loop on all the elements
    for(Integer elem_iter1 = 0; elem_iter1 < n_elements; ++elem_iter1) {
	
		const auto &el1 = mesh.elem(elem_iter1);
		const auto &el_nodes1=el1.nodes;
		std::vector<Integer> el_neighs;
		// find all the elements touching the element=elem_iter1
		for(Integer ii=0;ii<ManifoldDim+1;ii++)
			   for(Integer jj=0;jj<node_2_element[el_nodes1[ii]].size();jj++)
			   {el_neighs.push_back(node_2_element[el_nodes1[ii]][jj]);
			   } 
 
		// remove all the elements larger than  elem_iter1
		el_neighs.erase(std::remove_if(el_neighs.begin(), el_neighs.end(), [elem_iter1](int index_tmp) { return index_tmp>=elem_iter1; }), el_neighs.end());       
		// sort and make unique the neighborhood elements
		std::sort( el_neighs.begin(), el_neighs.end() );
		el_neighs.erase( std::unique( el_neighs.begin(), el_neighs.end() ), el_neighs.end() );
	
		// loop on all the entities of the actual element e1
		for(Integer iter_entity_e1=0;iter_entity_e1< entity_combinations_;iter_entity_e1++)
		{     
		 Combinations<ManifoldDim + 1, EntityDim>::generate(iter_entity_e1,entity_e1);
	 
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
			for(Integer iter_entity_e2=0;iter_entity_e2<entity_combinations_;iter_entity_e2++){
			
				Combinations<ManifoldDim + 1, EntityDim>::generate(iter_entity_e2,entity_e2);
			

				for(int nn=0;nn<EntityDim;nn++)
					entity_nodes2[nn]=el_nodes2[entity_e2[nn]];   
						
				std::sort(std::begin(entity_nodes2), std::end(entity_nodes2));  
			
				if (std::equal(std::begin(entity_nodes1), std::end(entity_nodes1), std::begin(entity_nodes2))){
				       // all the entities already computed on elem=elem_iter2
				       auto elem2entity=elem_2_entity_[elem_iter2];
				       // loop on each of these, to find the one that equals iter_entity_e1
				       for(auto i: elem2entity)
				          {
				          
				          auto nodes_i=entity_[i];
				          if(std::equal(std::begin(nodes_i), std::end(nodes_i), std::begin(entity_nodes1)))
				             {elem_2_entity_[elem_iter1].push_back(i);
				             break;}
				          }
					   found_entity=true;
					   break;}         
				}  
			
				if(found_entity==true)
				   break;
			
				}
		  // if no neighbborhood element share the same entity, then we must add this entity   
		  if(found_entity==false){
				   std::copy(std::begin(entity_nodes1), std::end(entity_nodes1), std::begin(tmp));
				   entity_.push_back( tmp);
				   elem_2_entity_[elem_iter1].push_back(entity_number);
				   entity_number++;
				   }	  
		  }
	  
		  }
		  
	 for(int elem_iter=0;elem_iter<n_elements;elem_iter++)
	  {std::cout<<std::endl;
	  std::cout<<" elem_iter = = "<<elem_iter<<std::endl;
	  for(int ii=0;ii<elem_2_entity_[elem_iter].size();ii++)
	      std::cout<<elem_2_entity_[elem_iter][ii]<<" ";
	  }
	 };



     
     
}
#endif //MARS_ENTITY_HPP
