#ifndef MARS_ENTITY_HPP
#define MARS_ENTITY_HPP

#include "mars_base.hpp"
#include "mars_base_entity.hpp"
#include "mars_mesh.hpp"
#include "mars_simplex.hpp"
#include <iostream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// ElemEntity<Elem, Integer EntityDim>                                                                                              ////////
////// Given a mesh, with elements of the kind Elem                                                                                     ////////
////// We want to define the entities of dimension 0<=EntityDim<=ManifoldDim                                                            ////////
////// The class has two fields:                                                                                                        ////////                                                        
//////  1) entity_2_elem_: for each entity id, returns a 2-D array A such that:                                                         ////////
//////  A[0]= id of one of the element to wich entity belongs                                                                           ////////
//////  A[1]= local number of the entity inside the element A[0]                                                                        ////////
//////  2) elem_2_entity: for each element it returns a entity_combinations-D array containing all the ids of the type of the entity    ////////                                    ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




namespace mars {

template<Integer Dim, Integer ManifoldDim>
class Mesh;


template<typename Elem, Integer EntityDim>
class ElemEntityCombinations;

template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
class ElemEntityCombinations<Simplex<Dim,ManifoldDim>,EntityDim>
{

public: 
static constexpr Integer value=Combinations<ManifoldDim+1,EntityDim+1>::value;

};


template<typename Elem, Integer EntityDim>
class ElemEntity{};


template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
class ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>
{

public: 

    constexpr static Integer num_of_points_ = EntityDim+1;
    inline constexpr static Integer entity_combinations()
    { 
        return ElemEntityCombinations<Simplex<Dim,ManifoldDim>,EntityDim>::value;
        //Combinations<ManifoldDim+1,EntityDim+1>::value;
    }

    inline constexpr static Integer num_of_points()
    { 
        return EntityDim+1;
    }
            
    ElemEntity(const Mesh<Dim,ManifoldDim>& mesh, const std::vector< std::vector<Integer> >& node_2_element); // constructor

    void init_elem_entity(const Mesh<Dim,ManifoldDim>& mesh,
                          const std::vector< std::vector<Integer> >& node_2_element, 
                          std::vector<std::array<Integer,2>> & entity_2_elem_, 
                          std::vector<std::array<Integer, entity_combinations() > > &elem_2_entity_,
                          Integer &size_); 

    inline const Integer size() const {return size_; };

    inline const std::vector<std::array<Integer,2>> entity_2_elem() const {return entity_2_elem_; };
    
    inline const std::array<Integer,2> entity_2_elem(Integer index) const {return entity_2_elem_[index]; };
    
    inline const std::vector<std::array<Integer, entity_combinations() > > elem_2_entity() const {return elem_2_entity_; };   
    
    inline const std::array<Integer, entity_combinations() > elem_2_entity(Integer index) const {return elem_2_entity_[index]; }; 


private:
    // entity_2_elem is a vector of 2D arrays: first component=1 elem, second component iter
    std::vector<std::array<Integer,2>> entity_2_elem_; 
    
    std::vector<std::array<Integer, entity_combinations() >> elem_2_entity_;
    
    Integer size_;


};








     
     
}
#endif //MARS_ENTITY_HPP
