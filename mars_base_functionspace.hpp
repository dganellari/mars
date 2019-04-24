#ifndef MARS_BASE_FUNCTIONSPACE_HPP
#define MARS_BASE_FUNCTIONSPACE_HPP

#include "mars_base.hpp"

namespace mars{

////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////---------- CONTINUITY TYPES ---------- //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
static constexpr Integer Discontinuous=0;
static constexpr Integer Continuous=1;
static constexpr Integer Mortar=2;

////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////---------- SPACE IDs ---------- //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
constexpr Integer LagrangeFE = -10;
constexpr Integer RaviartThomasFE = -11;
constexpr Integer NedelecFE = -12;
constexpr Integer GeneralSpace = -666;







template <typename Elem, Integer FEFamily, Integer Order,Integer Continuity, Integer NComponents>
class BaseElementFunctionSpace
{
    static const Integer family=FEFamily;
    static const Integer order=Order;  
    static const Integer n_components=NComponents;
    static const Integer continuity=Continuity;  
};

template <Integer Dim, Integer ManifoldDim, Integer FEFamily, Integer Order,Integer Continuity, Integer NComponents>
class BaseElementFunctionSpace<Simplex<Dim, ManifoldDim>,FEFamily,Order,Continuity,NComponents>
{public:
    static const Integer space_dim=Dim; 
    static const Integer manifold_dim=ManifoldDim;
    static const Integer family=FEFamily;
    static const Integer order=Order;  
    static const Integer n_components=NComponents;
    static const Integer continuity=Continuity;  
};



template <Integer Dim,Integer ManifoldDim, Integer FEFamily, Integer Order,Integer Continuity=Continuous, Integer NComponents=1>
class BaseFunctionSpace {
public:
   // static const Integer space_dim=Elem::Dim; 
   //  static const Integer manifold_dim=Elem::ManifoldDim;
 
    static const Integer space_dim=Dim; 
    static const Integer manifold_dim=ManifoldDim;
    static const Integer family=FEFamily;
    static const Integer order=Order;  
    static const Integer n_components=NComponents;
    static const Integer continuity=Continuity;  

};










}

#endif