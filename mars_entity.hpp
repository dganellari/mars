#ifndef MARS_ENTITY_HPP
#define MARS_ENTITY_HPP

#include "mars_base.hpp"
#include "mars_base_entity.hpp"
#include "mars_mesh.hpp"
#include "mars_simplex.hpp"
#include "mars_bisection.hpp"
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

template<Integer Dim, Integer ManifoldDim,class Implementation_>
class Mesh;

     // specialisation in case we choose non existent combinations 
    template<Integer N>
    class Combinations<N,N+1>
    {
     public:
      static const Integer value =0;  
    };

    template<Integer N>
    class Combinations<N,N+2>
    {
     public:
      static const Integer value =0;  
    };
    
template<typename Elem, Integer EntityDim>
class ElemEntityCombinations;

template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
class ElemEntityCombinations<Simplex<Dim,ManifoldDim>,EntityDim>
{

public: 
static constexpr Integer value=Combinations<ManifoldDim+1,EntityDim+1>::value;

};


template<typename Elem,Integer EntityDim>
class ElemEntityNPoints;


template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
class ElemEntityNPoints<Simplex<Dim,ManifoldDim>,EntityDim>
{
public:
    static constexpr Integer value=EntityDim+1;
};


template<typename Elem, Integer EntityDim>
class ElemEntity{};


template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
class ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>
{

public: 

    using MeshT=Mesh<Dim,ManifoldDim>;
    using BisectionT=Bisection<MeshT>;

    constexpr static Integer num_of_points_ = EntityDim+1;
    inline constexpr static Integer entity_combinations()
    { 
        return ElemEntityCombinations<Simplex<Dim,ManifoldDim>,EntityDim>::value;
        //Combinations<ManifoldDim+1,EntityDim+1>::value;
    }

    inline constexpr static Integer num_of_points()
    { 
        return ElemEntityNPoints<Simplex<Dim,ManifoldDim>,EntityDim>::value;
    }
    ElemEntity(){}
    ElemEntity(const Mesh<Dim,ManifoldDim>& mesh);
    ElemEntity(const std::shared_ptr<Mesh<Dim,ManifoldDim>>& mesh_ptr);
    ElemEntity(const Mesh<Dim,ManifoldDim>& mesh, const std::vector< std::vector<Integer> >& node_2_element, const Integer level=-1); // constructor
    ElemEntity(const std::shared_ptr<Mesh<Dim,ManifoldDim>>& mesh_ptr, const std::vector< std::vector<Integer> >& node_2_element, const Integer level=-1); 
   
    void init(const Mesh<Dim,ManifoldDim>& mesh);

    // void init_elem_entity(const Mesh<Dim,ManifoldDim>& mesh,
    //                       const std::vector< std::vector<Integer> >& node_2_element, 
    //                       std::vector<std::array<Integer,2>> & entity_2_elem_, 
    //                       std::vector<std::array<Integer, entity_combinations() > > &elem_2_entity_,
    //                       Integer &size_); 

    void init_elem_entity(const Mesh<Dim,ManifoldDim>& mesh,
                          const std::vector< std::vector<Integer> >& node_2_element,
                          const Integer level=-1); 


    inline Integer size() const {return size_; };

    inline const std::vector<std::array<Integer,2>>& entity_2_elem() const {return entity_2_elem_; };

    inline       std::vector<std::array<Integer,2>>& entity_2_elem()       {return entity_2_elem_; };
    
    inline const std::array<Integer,2>& entity_2_elem(Integer index) const {return entity_2_elem_[index]; };

    inline       std::array<Integer,2>& entity_2_elem(Integer index)       {return entity_2_elem_[index]; };
    
    inline const std::vector<std::vector<Integer > >& elem_2_entity() const {return elem_2_entity_; };

    inline       std::vector<std::vector<Integer > >& elem_2_entity()       {return elem_2_entity_; };   
    
    inline const std::vector<Integer >& elem_2_entity(Integer index) const {return elem_2_entity_[index]; };

    inline       std::vector<Integer >& elem_2_entity(Integer index)       {return elem_2_entity_[index]; }; 
    
    inline auto mesh_ptr(){return mesh_ptr_;}

    // inline void add_bisection(const BisectionT& bisection)
    // {}//bisection_ptr_=std::make_shared<BisectionT>(bisection);}

    inline void add_bisection(const std::shared_ptr<BisectionT> bisection_ptr)
    {bisection_ptr_=bisection_ptr;}

private:
    // entity_2_elem is a vector of 2D arrays: first component=1 elem, second component iter
    std::vector<std::array<Integer,2>> entity_2_elem_; 
    
    std::vector<std::vector<Integer>> elem_2_entity_;
    
    Integer size_;

    std::shared_ptr<Mesh<Dim,ManifoldDim>> mesh_ptr_;

    std::shared_ptr<BisectionT> bisection_ptr_;


};


// ENTITY == NODE 
template<Integer Dim,Integer ManifoldDim>
class ElemEntity<Simplex<Dim,ManifoldDim>,0>
{
public: 
    using MeshT=Mesh<Dim,ManifoldDim>;
    using BisectionT=Bisection<MeshT>;

    constexpr static Integer num_of_points_ = 1;
    inline constexpr static Integer entity_combinations() {return ManifoldDim+1;}
    inline constexpr static Integer num_of_points(){return 1;}
    ElemEntity(){}
    ElemEntity(const Mesh<Dim,ManifoldDim>& mesh);
    ElemEntity(const std::shared_ptr<Mesh<Dim,ManifoldDim>>& mesh_ptr);
    ElemEntity(const Mesh<Dim,ManifoldDim>& mesh, const std::vector< std::vector<Integer> >& node_2_element, const Integer level=-1);
    ElemEntity(const std::shared_ptr<Mesh<Dim,ManifoldDim>>& mesh_ptr, const std::vector< std::vector<Integer> >& node_2_element, const Integer level=-1); 
    // void init_elem_entity(const Mesh<Dim,ManifoldDim>& mesh,
    //                       const std::vector< std::vector<Integer> >& node_2_element, 
    //                       std::vector<std::array<Integer,2>> & entity_2_elem_, 
    //                       std::vector<std::array<Integer, entity_combinations() > > &elem_2_entity_,
    //                       Integer &size_); 
    
    void init(const Mesh<Dim,ManifoldDim>& mesh);

    void init_elem_entity(const Integer level);

    // void init_elem_entity(const Mesh<Dim,ManifoldDim>& mesh,
    //                       std::vector<std::array<Integer,2>> & entity_2_elem_, 
    //                       std::vector<std::array<Integer, entity_combinations() > > &elem_2_entity_,
    //                       Integer &size_); 
    void init_elem_entity(const Mesh<Dim,ManifoldDim>& mesh,
                          const std::vector< std::vector<Integer> >& node_2_element,
                          const Integer level=-1); 


    inline Integer size() const {return size_; };

    inline const std::vector<std::array<Integer,2>>& entity_2_elem() const {return entity_2_elem_; };

    inline       std::vector<std::array<Integer,2>>& entity_2_elem()       {return entity_2_elem_; };
    
    inline const std::array<Integer,2>& entity_2_elem(Integer index) const {return entity_2_elem_[index]; };

    inline       std::array<Integer,2>& entity_2_elem(Integer index)       {return entity_2_elem_[index]; };
    
    inline const std::vector<std::vector<Integer> >& elem_2_entity() const {return elem_2_entity_; };   

    inline       std::vector<std::vector<Integer> >& elem_2_entity()       {return elem_2_entity_; };   
    
    inline const std::vector<Integer>& elem_2_entity(Integer index) const {return elem_2_entity_[index]; }; 

    inline       std::vector<Integer>& elem_2_entity(Integer index)       {return elem_2_entity_[index]; }; 
    
    inline auto mesh_ptr(){return mesh_ptr_;}

    // inline void add_bisection(const BisectionT& bisection)
    // {}//bisection_ptr_=std::make_shared<BisectionT>(bisection);}

    inline void add_bisection(const std::shared_ptr<BisectionT> bisection_ptr)
    {bisection_ptr_=bisection_ptr;}

private:
    // entity_2_elem is a vector of 2D arrays: first component=1 elem, second component iter
    std::vector<std::array<Integer,2>> entity_2_elem_; 
    
    std::vector<std::vector<Integer>> elem_2_entity_;
    
    Integer size_;

    std::shared_ptr<Mesh<Dim,ManifoldDim>> mesh_ptr_;

    std::shared_ptr<BisectionT> bisection_ptr_;
};



 // ENTITY == ELEM 
template<Integer Dim,Integer ManifoldDim>
class ElemEntity<Simplex<Dim,ManifoldDim>,ManifoldDim>
{
public: 
    using MeshT=Mesh<Dim,ManifoldDim>;
    using BisectionT=Bisection<MeshT>;

    constexpr static Integer num_of_points_ = ManifoldDim+1;
    inline constexpr static Integer entity_combinations() {return 1;}
    inline constexpr static Integer num_of_points(){return ManifoldDim+1;}
    ElemEntity(){}
    ElemEntity(const Mesh<Dim,ManifoldDim>& mesh);
    ElemEntity(const std::shared_ptr<Mesh<Dim,ManifoldDim>>& mesh_ptr);
    ElemEntity(const Mesh<Dim,ManifoldDim>& mesh, const std::vector< std::vector<Integer> >& node_2_element, const Integer level=-1);
    ElemEntity(const std::shared_ptr<Mesh<Dim,ManifoldDim>>& mesh_ptr, const std::vector< std::vector<Integer> >& node_2_element, const Integer level=-1);
    
    void init(const Mesh<Dim,ManifoldDim>& mesh);

    // void init_elem_entity(const Mesh<Dim,ManifoldDim>& mesh,
    //                       const std::vector< std::vector<Integer> >& node_2_element, 
    //                       std::vector<std::array<Integer,2>> & entity_2_elem_, 
    //                       std::vector<std::array<Integer, entity_combinations() > > &elem_2_entity_,
    //                       Integer &size_); 

    void init_elem_entity(const Mesh<Dim,ManifoldDim>& mesh,
                          const std::vector< std::vector<Integer> >& node_2_element,
                          const Integer level=-1);

    inline Integer size() const {return size_; };

    inline       std::vector<std::array<Integer,2>>& entity_2_elem()       {return entity_2_elem_; };

    inline const std::vector<std::array<Integer,2>>& entity_2_elem() const {return entity_2_elem_; };
    
    inline const std::array<Integer,2>& entity_2_elem(Integer index) const {return entity_2_elem_[index]; };

    inline       std::array<Integer,2>& entity_2_elem(Integer index)       {return entity_2_elem_[index]; };
    
    inline const std::vector<std::vector<Integer> >& elem_2_entity() const {return elem_2_entity_; };  

    inline       std::vector<std::vector<Integer> >& elem_2_entity()       {return elem_2_entity_; };   
    
    inline const std::vector<Integer >& elem_2_entity(Integer index) const {return elem_2_entity_[index]; };

    inline       std::vector<Integer >& elem_2_entity(Integer index)  {return elem_2_entity_[index]; }; 
    
    inline auto mesh_ptr(){return mesh_ptr_;}

    // inline void add_bisection(const BisectionT& bisection)
    // {}//bisection_ptr_=std::make_shared<BisectionT>(bisection);}

    inline void add_bisection(const std::shared_ptr<BisectionT> bisection_ptr)
    {bisection_ptr_=bisection_ptr;}

private:
    // entity_2_elem is a vector of 2D arrays: first component=1 elem, second component iter
    std::vector<std::array<Integer,2>> entity_2_elem_; 
    
    std::vector<std::vector<Integer>> elem_2_entity_;
    
    Integer size_;

    std::shared_ptr<Mesh<Dim,ManifoldDim>> mesh_ptr_;

    std::shared_ptr<BisectionT> bisection_ptr_;
};

     
     
}
#endif //MARS_ENTITY_HPP
