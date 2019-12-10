////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
//////// Written by Gabriele Rovi (April 2019)                                              ////////
//////// We define a general base space, indipendent of the type of element                 ////////
//////// We specialize the base space to a given element with BaseElementFunctionSpace      ////////
//////// For example, elem=simplex<dim,manifoldim>                                          ////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


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

////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////----------  BASE  SPACE  ---------- //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


template <Integer FEFamily_, Integer Order_,Integer Continuity_=Continuous, Integer NComponents_=1>
class BaseFunctionSpace
{
public:
    static const Integer FEFamily=FEFamily_;
    static const Integer Order=Order_;  
    static const Integer NComponents=NComponents_;
    static const Integer Continuity=Continuity_;  
};


template<Integer Order,Integer NComponents>
using Lagrange=BaseFunctionSpace<LagrangeFE,Order,Continuous,NComponents>;

template<Integer NComponents>
using Lagrange1=BaseFunctionSpace<LagrangeFE,1,Continuous,NComponents>;
template<Integer NComponents>
using Lagrange2=BaseFunctionSpace<LagrangeFE,2,Continuous,NComponents>;
template<Integer NComponents>
using Lagrange3=BaseFunctionSpace<LagrangeFE,3,Continuous,NComponents>;


template<Integer Order,Integer NComponents>
using LagrangeDG=BaseFunctionSpace<LagrangeFE,Order,Discontinuous,NComponents>;

template<Integer NComponents>
using Lagrange1DG=BaseFunctionSpace<LagrangeFE,1,Discontinuous,NComponents>;
template<Integer NComponents>
using Lagrange2DG=BaseFunctionSpace<LagrangeFE,2,Discontinuous,NComponents>;
template<Integer NComponents>
using Lagrange3DG=BaseFunctionSpace<LagrangeFE,3,Discontinuous,NComponents>;

template<Integer Order,Integer NComponents>
using RT=BaseFunctionSpace<RaviartThomasFE,Order,Continuous,NComponents>;

template<Integer NComponents>
using RT0=BaseFunctionSpace<RaviartThomasFE,0,Continuous,NComponents>;
template<Integer NComponents>
using RT1=BaseFunctionSpace<RaviartThomasFE,1,Continuous,NComponents>;

template<Integer Order,Integer NComponents>
using RTDG=BaseFunctionSpace<RaviartThomasFE,Order,Discontinuous,NComponents>;

template<Integer NComponents>
using RT0DG=BaseFunctionSpace<RaviartThomasFE,0,Discontinuous,NComponents>;
template<Integer NComponents>
using RT1DG=BaseFunctionSpace<RaviartThomasFE,1,Discontinuous,NComponents>;

////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////------- BASE ELEMENT SPACE -------- //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename Elem, Integer FEFamily, Integer Order,Integer Continuity, Integer NComponents>
class BaseElementFunctionSpace
{
    static const Integer family=FEFamily;
    static const Integer order=Order;  
    static const Integer n_components=NComponents;
    static const Integer continuity=Continuity;  
};
////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////------- BASE SIMPLEX SPACE ---------- ////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
//template <Integer Dim, Integer ManifoldDim, Integer FEFamily, Integer Order,Integer Continuity, Integer NComponents>
template <Integer Dim_, Integer ManifoldDim_, Integer FEFamily_, Integer Order_,Integer Continuity_, Integer NComponents_>
class BaseElementFunctionSpace<Simplex<Dim_, ManifoldDim_>,FEFamily_,Order_,Continuity_,NComponents_>
{public:
    static const Integer Dim=Dim_; 
    static const Integer ManifoldDim=ManifoldDim_;
    static const Integer FEFamily=FEFamily_;
    static const Integer Order=Order_;  
    static const Integer NComponents=NComponents_;
    static const Integer Continuity=Continuity_; 
    using  Elem=mars::Simplex<Dim, ManifoldDim>; 
};


////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////---------- ELEMENT SPACE ---------- //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename Elem, Integer FEFamily, Integer Order,Integer Continuity=Continuous, Integer NComponents=1>
class ElementFunctionSpace: public BaseElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>
{};

////////////////////////////////////////////////////////////////////////////////////////////////////
//////// using ElemFunctionSpace                                                            ////////
//////// is used to go from BaseFunctionSpace to an ElementFunctionSpace                    ////////
////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename Elem, typename BaseFunctionSpace>
using ElemFunctionSpace=ElementFunctionSpace<Elem,
                                             BaseFunctionSpace::FEFamily,
                                             BaseFunctionSpace::Order,
                                             BaseFunctionSpace::Continuity,
                                             BaseFunctionSpace::NComponents >;

template<typename FunctionSpace>
using Elem2FunctionSpace=BaseFunctionSpace<FunctionSpace::FEFamily,
                                                        FunctionSpace::Order,
                                                        FunctionSpace::Continuity,
                                                        FunctionSpace::NComponents >;

}

#endif