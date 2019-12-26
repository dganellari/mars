#ifndef MARS_CONNECTIVITY_HPP
#define MARS_CONNECTIVITY_HPP

#include "mars_base.hpp"
#include "mars_mesh.hpp"
#include "mars_entity.hpp"
#include "mars_tracker.hpp"
#include "mars_bisection.hpp"
#include <iostream>
#include <vector>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// Written by Gabriele Rovi (April 2019)                                                                                            ////////
////// We define the following class:                                                                                                   ////////
////// Connectivity<typename Elem, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>                                ////////
////// Specialized for the simplex case:                                                                                                ////////
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
//////        2) if SubEntityDimFrom = EntityDimFrom: the set is the entity e1 itself                                                   ////////
//////        3) The class with EntityDimFrom=SubEntityDimFrom=0 and EntityDimTo=ManifoldDim is defined separately                      ////////
//////           Indeed it will store a vector of vectors, whose component (the node id) returns the neighborhood elements              ////////
//////        4) For all the other classes, given the index entity, the connectivity is computed on the fly                             ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace mars{

void connectivity_example();

template <Integer Dim, Integer ManifoldDim,class Implementation_>
class Mesh;

template<typename Elem, Integer EntityDim>
class ElemEntity;

template<Integer Dim,Integer ManifoldDim,Integer EntityDim>
class ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>;

template <typename Elem, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>
class ElemConnectivity;

template <Integer Dim, Integer ManifoldDim, Integer EntityDimFrom, Integer SubEntityDimFrom, Integer EntityDimTo>
class ElemConnectivity<Mesh<Dim,ManifoldDim>,EntityDimFrom,SubEntityDimFrom,EntityDimTo >
{
public:
    using MeshT=Mesh<Dim,ManifoldDim>;
    using BisectionT=Bisection<MeshT>;
        inline std::vector< std::vector<Integer> > node2elem() const{return node2elem_;};

        std::vector<Integer> compute(const ElemEntity<Simplex<Dim,ManifoldDim>,EntityDimFrom> &entity_from,
                                     const Integer& index_from,
                                     const ElemEntity<Simplex<Dim,ManifoldDim>,EntityDimTo>   &entity_to);
                     
        inline void add_bisection(const BisectionT& bisection){bisection_ptr_=std::make_shared<BisectionT>(bisection);}
        
        ElemConnectivity(const Mesh<Dim,ManifoldDim> &mesh,const std::vector<std::vector<Integer>> &node2elem,const Integer level=-1):
                // mesh_(mesh),
                mesh_ptr_(std::make_shared<MeshT>(mesh)),
                node2elem_(node2elem),
                bisection_ptr_(NULL),
                level_(level)
                {};
        
private:
        const std::vector< std::vector<Integer> > &node2elem_;
        // const Mesh<Dim,ManifoldDim> &mesh_;
        std::shared_ptr<MeshT> mesh_ptr_;
        std::shared_ptr<BisectionT> bisection_ptr_;
        Integer level_;
};



template <Integer Dim, Integer ManifoldDim>
class ElemConnectivity<Mesh<Dim, ManifoldDim>,0,0,ManifoldDim>
{
public:
        using MeshT=Mesh<Dim,ManifoldDim>;
        using BisectionT=Bisection<MeshT>;

        inline const std::vector< std::vector<Integer> >& val() const{return val_;};
        inline std::vector< std::vector<Integer> >&       val()      {return val_;};

        inline void add_bisection(const BisectionT& bisection){bisection_ptr_=std::make_shared<BisectionT>(bisection);}

        void init();
        
        ElemConnectivity(const Mesh<Dim,ManifoldDim>& mesh,const Integer level=-1):
        mesh_ptr_(std::make_shared<MeshT>(mesh)),
        bisection_ptr_(NULL),
        level_(level)
        {
        if(level_==-1)
        init();
        };   
                
private:
        std::shared_ptr<MeshT> mesh_ptr_;
        std::vector<std::vector<Integer>> val_;
        std::shared_ptr<BisectionT> bisection_ptr_;
        Integer level_;
};


template<typename MeshT>
using NodeToElem=mars::ElemConnectivity<MeshT,0,0,MeshT::ManifoldDim>;


}

#endif