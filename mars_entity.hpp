#ifndef MARS_ENTITY_HPP
#define MARS_ENTITY_HPP

#include "mars_base.hpp"
#include "mars_base_entity.hpp"
#include "mars_mesh.hpp"
#include <iostream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// Entity<Integer Dim, Integer ManifoldDim, Integer EntityDim>                                                                      ////////
////// Given a mesh, with elements whose dimension is ManifoldDim, defined in a Dim-dimensional space                                   ////////
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
    
template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
class Entity : public BaseEntity{

public: 

    constexpr static Integer num_of_points_ = EntityDim+1;
    inline constexpr static Integer entity_combinations()
    { 
        return Combinations<ManifoldDim+1,EntityDim+1>::value;
    }

    inline constexpr static Integer num_of_points()
    { 
        return EntityDim+1;
    }
            
    Entity(const Mesh<Dim,ManifoldDim> mesh, const std::vector< std::vector<Integer> >node_2_element); // constructor

    void init_entity(const Mesh<Dim,ManifoldDim> mesh,const std::vector< std::vector<Integer> > node_2_element, 
                    std::vector<std::array<Integer,2>> & entity_2_elem_, 
                    std::vector<std::array<Integer, entity_combinations() > > &elem_2_entity_,
                    Integer &entity_nums_); 

    inline const Integer entity_nums() const {return entity_nums_; };

    inline std::vector<std::array<Integer,2>> entity_2_elem() const {return entity_2_elem_; };
    
    inline std::array<Integer,2> entity_2_elem(Integer index) const {return entity_2_elem_[index]; };
    
    inline std::vector<std::array<Integer, entity_combinations() > > elem_2_entity() const {return elem_2_entity_; };   
    
    inline std::array<Integer, entity_combinations() > elem_2_entity(Integer index) const {return elem_2_entity_[index]; }; 


private:
    // entity_2_elem is a vector of 2D arrays: first component=1 elem, second component iter
    std::vector<std::array<Integer,2>> entity_2_elem_; 
    
    std::vector<std::array<Integer, entity_combinations() >> elem_2_entity_;
    
    Integer entity_nums_;
    };





using EdgeMap1=Entity<1,1,1>;
using EdgeMap2=Entity<2,2,1>;
using EdgeMap3=Entity<3,3,1>;
using EdgeMap4=Entity<4,4,1>;

using EdgeMap21=Entity<2,1,1>;

using EdgeMap31=Entity<3,1,1>;
using EdgeMap32=Entity<3,2,1>;

using EdgeMap41=Entity<4,1,1>;
using EdgeMap42=Entity<4,2,1>;
using EdgeMap43=Entity<4,3,1>;


using TriangleMap2=Entity<2,2,2>;
using TriangleMap3=Entity<3,3,2>;
using TriangleMap4=Entity<4,4,2>;

using TriangleMap22=Entity<2,2,2>;
using TriangleMap32=Entity<3,2,2>;
using TriangleMap42=Entity<4,2,2>;
using TriangleMap43=Entity<4,3,2>;

using TetrahedronMap3=Entity<3,3,3>;
using TetrahedronMap4=Entity<4,4,3>;


using PentatopeMap4=Entity<4,4,4>;
     
     
}
#endif //MARS_ENTITY_HPP
