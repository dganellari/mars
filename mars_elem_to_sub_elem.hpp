#ifndef MARS_ELEM_TO_SUB_ELEM_HPP
#define MARS_ELEM_TO_SUB_ELEM_HPP
#include "mars_base.hpp"

namespace mars{


/////////////////////////////////////////////////////////////////////////////
//// Go grom Elem<Dim,ManifoldDim> to Elem<Dim,SubElemDim>               ////
/////////////////////////////////////////////////////////////////////////////
template<typename Elem, Integer SubElemDim>
class ElemToSubElemHelper;

template<template<Integer,Integer>class Elem_,Integer Dim, Integer ManifoldDim,Integer SubElemDim>
class ElemToSubElemHelper<Elem_<Dim,ManifoldDim>,SubElemDim>
{
 public:
  static_assert(SubElemDim<=ManifoldDim && " In FromElemToSubElem: a SubElemDim dim must be smaller or equal to the entity ManifoldDim");
  using type=Elem_<Dim,SubElemDim>; 
};
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
