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


	std::vector<std::array<Integer,EntityDim> > init_entity(Mesh<Dim,ManifoldDim> mesh)
	 {
	 std::vector<std::array<Integer,EntityDim> > entity;
	 std::array<Integer,EntityDim> tmp;
	 
	 for(Integer ii=0;ii<EntityDim;ii++)
	    tmp[ii]=ii;
	    
	 entity.push_back(tmp);
	 return entity;
	 };
	 
	 Entity(Mesh<Dim,ManifoldDim> mesh):
     entity_(init_entity(Mesh<Dim,ManifoldDim> mesh))
     {};
	 
	 
std::vector<Integer> value(void) const {return entity_; };	
private:
std::vector<std::array<Integer,EntityDim> > entity_;
 	};


//  std::vector<std::array<Integer,EntityDim> > Entity<ManifoldDim,EntityDim>::init_entity(Mesh<Dim,ManifoldDim> mesh)
//  {
//  std::vector<std::array<Integer,EntityDim>, 1> entity;
//  //entity.push_back(std::array<Integer>);
//  return entity;
//  }


// Entity::Entity(Mesh<ManifoldDim,ManifoldDim> mesh):
// entity_(init_entity(Mesh<ManifoldDim,ManifoldDim> mesh))
// {}



}
#endif //MARS_ENTITY_HPP
