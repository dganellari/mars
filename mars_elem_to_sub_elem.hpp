#ifndef MARS_ELEM_TO_SUB_ELEM_HPP
#define MARS_ELEM_TO_SUB_ELEM_HPP
#include "mars_base.hpp"

namespace mars{


/////////////////////////////////////////////////////////////////////////////
//// Go grom Elem<Dim,ManifoldDim> to Elem<Dim,SubElemDim>               ////
/////////////////////////////////////////////////////////////////////////////
template<typename Elem, Integer SubElemDim>
class ChangeElemDim;

template<template<Integer,Integer>class Elem_,Integer Dim, Integer ManifoldDim,Integer SubElemDim>
class ChangeElemDim<Elem_<Dim,ManifoldDim>,SubElemDim>
{
 public:
  static_assert(SubElemDim<=Dim && " In ChangeElemDim: a SubElemDim dim must be smaller or equal to the entity Dim");
  using type=Elem_<Dim,SubElemDim>; 
};

template<typename Elem>
using FromBoundaryToVolumetricElem=typename ChangeElemDim<Elem,Elem::ManifoldDim+1>::type;
template<typename Elem>
using FromVolumetricToBoundaryElem=typename ChangeElemDim<Elem,Elem::ManifoldDim-1>::type;

/////////////////////////////////////////////////////////////////////////////
//// VolumeOrSurfaceElem return the volume/surface Element if true/false ////
/////////////////////////////////////////////////////////////////////////////
template<typename Elem, bool=true>
class VolumeOrSurfaceElem;

template<template<Integer Dim,Integer ManifoldDim>class Elem_,Integer Dim, Integer ManifoldDim>
class VolumeOrSurfaceElem<Elem_<Dim,ManifoldDim>,true>
{
 public:
  using type=Elem_<Dim,ManifoldDim>; 
};

template<template<Integer Dim,Integer ManifoldDim>class Elem_,Integer Dim, Integer ManifoldDim>
class VolumeOrSurfaceElem<Elem_<Dim,ManifoldDim>,false>
{
 public:
  using type=Elem_<Dim,ManifoldDim-1>; 
};


}


#endif