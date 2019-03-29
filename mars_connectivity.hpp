#ifndef MARS_CONNECTIVITY_HPP
#define MARS_CONNECTIVITY_HPP

#include "mars_base.hpp"
#include "mars_mesh.hpp"
#include "mars_entity.hpp"
#include <iostream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)																							////////
////// We define the following class:																									////////
////// Connectivity<Integer Dim, Integer ManifoldDim, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>				////////
////// Given a mesh, with elements whose dimension is ManifoldDim, defined in a Dim-dimensional space									////////
////// We want to define the connectivity between different entities 																	////////
////// Given an entity e1 of dimension EntityDimFrom, we want to find all the other entities e2 of dimension EntityDimTo   				////////
////// Such that e2 share at least SubEntityDimFrom+1 points with EntityDimFrom															////////
////// Example 1) Dim=ManifoldDim=3D, e1=triangle, SubEntityDimFrom=0, e2=triangle:		 												////////
//////            all the triangles e2 which share a node with the triangle e1															////////
////// Example 2) Dim=ManifoldDim=3D, e1=triangle, SubEntityDimFrom=1, e2=triangle:														////////
//////            all the triangles e2 which share an edge with the triangle e1								 							////////
////// Example 3) Dim=ManifoldDim=3D, e1=triangle, SubEntityDimFrom=2, e2=triangle:														////////
//////            all the triangles e2 which share a triangle with the triangle e1 -> e1=e2, known a priori								////////
////// Rules: 1) if SubEntityDimFrom > EntityDimFrom: the set is empty																	////////
//////        2) if SubEntityDimFrom = EntityDimFrom: the set is the entity e1 itself													////////
//////        3) The class with EntityDimFrom=SubEntityDimFrom=0 and EntityDimTo=ManifoldDim is defined separately						////////
//////           Indeed it will store a vector of vectors, whose component (the node id) returns the neighborhood elements              ////////
//////        4) For all the other classes, given the index entity, the connectivity is computed on the fly   							////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Issue: I am copying mesh and node2elem_. Maybe it is better to have a make_shared pointer to these quantities
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace mars{


template <Integer Dim, Integer ManifoldDim>
class Mesh;

template <Integer Dim, Integer ManifoldDim, Integer EntityDim>
class Entity;

template <Integer Dim, Integer ManifoldDim, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>
class Connectivity
{
public:
		inline std::vector< std::vector<Integer> > node2elem() const{return node2elem_;};

		inline std::vector<Integer> compute(const Entity<Dim,ManifoldDim,EntityDimFrom> &entity_from,
											const Integer index_from,
											const Entity<Dim,ManifoldDim,EntityDimTo>   &entity_to);
					 
		Connectivity(const Mesh<Dim,ManifoldDim> mesh,std::vector<std::vector<Integer>> node2elem):
				mesh_(mesh),
				node2elem_(node2elem)
				{};
		
private:
		std::vector< std::vector<Integer> > node2elem_;
		Mesh<Dim,ManifoldDim> mesh_;
};



template <Integer Dim, Integer ManifoldDim>
class Connectivity<Dim, ManifoldDim,0,0,ManifoldDim>
{
public:
		inline std::vector< std::vector<Integer> > val() const{return val_;};

		void init(const Mesh<Dim,ManifoldDim> mesh);

		Connectivity(const Mesh<Dim,ManifoldDim> mesh) {init(mesh);};	
				
private:
		std::vector<std::vector<Integer>> val_;
};



using NodeToElem1=mars::Connectivity<1,1,0,0,1>;
using NodeToElem2=mars::Connectivity<2,2,0,0,2>;
using NodeToElem3=mars::Connectivity<3,3,0,0,3>;
using NodeToElem4=mars::Connectivity<4,4,0,0,4>;
}

#endif