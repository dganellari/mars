#ifndef MARS_ENTITY_IMPL_HPP
#define MARS_ENTITY_IMPL_HPP

#include "mars_entity.hpp"

namespace mars {

template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
Entity<Dim,ManifoldDim,EntityDim>::Entity(const Mesh<Dim,ManifoldDim> mesh,const std::vector< std::vector<Integer> > node_2_element)
     {init_entity(mesh,node_2_element,entity_2_elem_,elem_2_entity_); };
     
     
template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
void Entity<Dim,ManifoldDim,EntityDim>::init_entity(const Mesh<Dim,ManifoldDim> mesh,const std::vector< std::vector<Integer> > node_2_element,
                                    std::vector<std::array<Integer,2>> &entity_2_elem_, std::vector<std::array<Integer, entity_combinations() >> &elem_2_entity_)
 {
    const Integer n_elements = mesh.n_elements();
    Integer entity_number=0;
    //std::array<Integer, EntityDim+1> tmp;
    std::array<Integer, 2> tmp;
    Integer entity_e1[EntityDim+1];
    Integer entity_e2[EntityDim+1];
    Integer entity_e3[EntityDim+1];
    std::array<Integer, EntityDim+1> entity_nodes1;
    std::array<Integer, EntityDim+1> entity_nodes2;
    std::array<Integer, EntityDim+1> entity_nodes3;
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
		for(Integer iter_entity_e1=0;iter_entity_e1< entity_combinations();iter_entity_e1++)
		{     
		 Combinations<ManifoldDim + 1, EntityDim+1>::generate(iter_entity_e1,entity_e1);
	 
		 // build the face nodes of iter_entity_e1
		 for(int nn=0;nn<EntityDim+1;nn++)
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
			for(Integer iter_entity_e2=0;iter_entity_e2<entity_combinations();iter_entity_e2++){
			
				Combinations<ManifoldDim + 1, EntityDim+1>::generate(iter_entity_e2,entity_e2);
			
				for(int nn=0;nn<EntityDim+1;nn++)
					entity_nodes2[nn]=el_nodes2[entity_e2[nn]];   
						
				std::sort(std::begin(entity_nodes2), std::end(entity_nodes2));  
			
				if (std::equal(std::begin(entity_nodes1), std::end(entity_nodes1), std::begin(entity_nodes2))){
				       // all the entities already computed on elem=elem_iter2
				       auto elem2entity=elem_2_entity_[elem_iter2];
				       // loop on each of these, to find the one that equals iter_entity_e1
				       for(auto i: elem2entity)
				          {
				          // we generate the combination of dimension=EntityDim with local numbering=entity_2_elem_[i][1]
				          Combinations<ManifoldDim + 1, EntityDim+1>::generate(entity_2_elem_[i][1],entity_e3);
				          // then we find the sorted nodes
			              for(int nn=0;nn<EntityDim+1;nn++)
		         			entity_nodes3[nn]=mesh.elem(entity_2_elem_[i][0]).nodes[entity_e3[nn]]; 
					      std::sort(std::begin(entity_nodes3), std::end(entity_nodes3));  
									          
				          //auto nodes_i=entity_2_elem_[i][0];
				          if(std::equal(std::begin(entity_nodes3), std::end(entity_nodes3), std::begin(entity_nodes1)))
				             {
				             elem_2_entity_[elem_iter1][iter_entity_e1]=i;
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
				   //std::copy(std::begin(entity_nodes1), std::end(entity_nodes1), std::begin(tmp));
				   //entity_.push_back( tmp);
				   tmp[0]=elem_iter1; 
				   tmp[1]=iter_entity_e1;
				   entity_2_elem_.push_back(tmp);
				   elem_2_entity_[elem_iter1][iter_entity_e1]=entity_number;				   
				   entity_number++;
				   }	  
		     }	  
		  }    
	 };



     
     
}
#endif //MARS_ENTITY_HPP
