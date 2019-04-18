#ifndef MARS_BASE_FUNCTIONSPACE_HPP
#define MARS_BASE_FUNCTIONSPACE_HPP

#include "mars_base.hpp"

namespace mars{

template <Integer Dim, Integer ManifoldDim, Integer SpaceID, Integer Order, Integer FunctionSpaceDim=1>
class BaseFunctionSpace {
public:
    static const Integer space_dim=Dim; 
    static const Integer manifold_dim=ManifoldDim;
    static const Integer id=SpaceID;
    static const Integer order=Order;  
    static const Integer functionspace_dim=FunctionSpaceDim;

};


constexpr Integer Lagrange = -10;
constexpr Integer RaviartThomas = -11;
constexpr Integer Nedelec = -12;
constexpr Integer GeneralSpace = -666;
}

#endif