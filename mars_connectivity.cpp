#include "mars_connectivity_impl.hpp"

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



namespace mars{


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D -----   //////////
////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D -----   //////////
////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D ----- ////// ----- 2D -----   //////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle  /////////////
/////////////////////////////////////////////////////////////

// mesh 2D  Sub=Node, To=Node
template class Connectivity<2,2,0,0,0>;
template class Connectivity<2,2,1,0,0>;
template class Connectivity<2,2,2,0,0>;
// mesh 2D  Sub=Edge, To=Node
template class Connectivity<2,2,0,1,0>;
template class Connectivity<2,2,1,1,0>;
template class Connectivity<2,2,2,1,0>;
// mesh 2D  Sub=Triangle, To=Node
template class Connectivity<2,2,0,2,0>;
template class Connectivity<2,2,1,2,0>;
template class Connectivity<2,2,2,2,0>;

/////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle  /////////////
/////////////////////////////////////////////////////////////

// mesh 2D  Sub=Node, To=Edge
template class Connectivity<2,2,0,0,1>;
template class Connectivity<2,2,1,0,1>;
template class Connectivity<2,2,2,0,1>;
// mesh 2D  Sub=Edge, To=Edge
template class Connectivity<2,2,0,1,1>;
template class Connectivity<2,2,1,1,1>;
template class Connectivity<2,2,2,1,1>;
// mesh 2D  Sub=Triangle, To=Edge
template class Connectivity<2,2,0,2,1>;
template class Connectivity<2,2,1,2,1>;
template class Connectivity<2,2,2,2,1>;


/////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle  /////////////
/////////////////////////////////////////////////////////////

// mesh 2D  Sub=Node, To=Triangle
template class Connectivity<2,2,0,0,2>;
template class Connectivity<2,2,1,0,2>;
template class Connectivity<2,2,2,0,2>;
// mesh 2D  Sub=Edge, To=Triangle
template class Connectivity<2,2,0,1,2>;
template class Connectivity<2,2,1,1,2>;
template class Connectivity<2,2,2,1,2>;
// mesh 2D  Sub=Triangle, To=Triangle
template class Connectivity<2,2,0,2,2>;
template class Connectivity<2,2,1,2,2>;
template class Connectivity<2,2,2,2,2>;







////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D -----   //////////
////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D -----   //////////
////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D ----- ////// ----- 3D -----   //////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle, Tetrahedron  /////////////
//////////////////////////////////////////////////////////////////////////


// mesh 3D  Sub=Node, To=Node
template class Connectivity<3,3,0,0,0>;
template class Connectivity<3,3,1,0,0>;
template class Connectivity<3,3,2,0,0>;
template class Connectivity<3,3,3,0,0>;
// mesh 3D  Sub=Edge, To=Node
template class Connectivity<3,3,0,1,0>;
template class Connectivity<3,3,1,1,0>;
template class Connectivity<3,3,2,1,0>;
template class Connectivity<3,3,3,1,0>;
// mesh 3D  Sub=Triangle, To=Node
template class Connectivity<3,3,0,2,0>;
template class Connectivity<3,3,1,2,0>;
template class Connectivity<3,3,2,2,0>;
template class Connectivity<3,3,3,2,0>;
// mesh 3D  Sub=Tetrahedron, To=Node
template class Connectivity<3,3,0,3,0>;
template class Connectivity<3,3,1,3,0>;
template class Connectivity<3,3,2,3,0>;
template class Connectivity<3,3,3,3,0>;

//////////////////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle, Tetrahedron  /////////////
//////////////////////////////////////////////////////////////////////////

// mesh 3D  Sub=Node, To=Edge
template class Connectivity<3,3,0,0,1>;
template class Connectivity<3,3,1,0,1>;
template class Connectivity<3,3,2,0,1>;
template class Connectivity<3,3,3,0,1>;
// mesh 3D  Sub=Edge, To=Edge
template class Connectivity<3,3,0,1,1>;
template class Connectivity<3,3,1,1,1>;
template class Connectivity<3,3,2,1,1>;
template class Connectivity<3,3,3,1,1>;
// mesh 3D  Sub=Triangle, To=Edge
template class Connectivity<3,3,0,2,1>;
template class Connectivity<3,3,1,2,1>;
template class Connectivity<3,3,2,2,1>;
template class Connectivity<3,3,3,2,1>;
// mesh 3D  Sub=Tetrahedron, To=Edge
template class Connectivity<3,3,0,3,1>;
template class Connectivity<3,3,1,3,1>;
template class Connectivity<3,3,2,3,1>;
template class Connectivity<3,3,3,3,1>;

//////////////////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle, Tetrahedron  /////////////
//////////////////////////////////////////////////////////////////////////

// mesh 3D  Sub=Node, To=Triangle
template class Connectivity<3,3,0,0,2>;
template class Connectivity<3,3,1,0,2>;
template class Connectivity<3,3,2,0,2>;
template class Connectivity<3,3,3,0,2>;
// mesh 3D  Sub=Edge, To=Triangle
template class Connectivity<3,3,0,1,2>;
template class Connectivity<3,3,1,1,2>;
template class Connectivity<3,3,2,1,2>;
template class Connectivity<3,3,3,1,2>;
// mesh 3D  Sub=Triangle, To=Triangle
template class Connectivity<3,3,0,2,2>;
template class Connectivity<3,3,1,2,2>;
template class Connectivity<3,3,2,2,2>;
template class Connectivity<3,3,3,2,2>;
// mesh 3D  Sub=Tetrahedron, To=Triangle
template class Connectivity<3,3,0,3,2>;
template class Connectivity<3,3,1,3,2>;
template class Connectivity<3,3,2,3,2>;
template class Connectivity<3,3,3,3,2>;

//////////////////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle, Tetrahedron  /////////////
//////////////////////////////////////////////////////////////////////////

// mesh 3D  Sub=Node, To=Tetrahedron
template class Connectivity<3,3,0,0,3>;
template class Connectivity<3,3,1,0,3>;
template class Connectivity<3,3,2,0,3>;
template class Connectivity<3,3,3,0,3>;
// mesh 3D  Sub=Edge, To=Tetrahedron
template class Connectivity<3,3,0,1,3>;
template class Connectivity<3,3,1,1,3>;
template class Connectivity<3,3,2,1,3>;
template class Connectivity<3,3,3,1,3>;
// mesh 3D  Sub=Triangle, To=Tetrahedron
template class Connectivity<3,3,0,2,3>;
template class Connectivity<3,3,1,2,3>;
template class Connectivity<3,3,2,2,3>;
template class Connectivity<3,3,3,2,3>;
// mesh 3D  Sub=Tetrahedron, To=Tetrahedron
template class Connectivity<3,3,0,3,3>;
template class Connectivity<3,3,1,3,3>;
template class Connectivity<3,3,2,3,3>;
template class Connectivity<3,3,3,3,3>;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D -----   //////////
////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D -----   //////////
////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D ----- ////// ----- 4D -----   //////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle, Tetrahedron, Pentatope  /////////////
/////////////////////////////////////////////////////////////////////////////////////

// mesh 4D  Sub=Node, To=Node
template class Connectivity<4,4,0,0,0>;
template class Connectivity<4,4,1,0,0>;
template class Connectivity<4,4,2,0,0>;
template class Connectivity<4,4,3,0,0>;
template class Connectivity<4,4,4,0,0>;
// mesh 4D  Sub=Edge, To=Node
template class Connectivity<4,4,0,1,0>;
template class Connectivity<4,4,1,1,0>;
template class Connectivity<4,4,2,1,0>;
template class Connectivity<4,4,3,1,0>;
template class Connectivity<4,4,4,1,0>;
// mesh 4D  Sub=Triangle, To=Node
template class Connectivity<4,4,0,2,0>;
template class Connectivity<4,4,1,2,0>;
template class Connectivity<4,4,2,2,0>;
template class Connectivity<4,4,3,2,0>;
template class Connectivity<4,4,4,2,0>;
// mesh 4D  Sub=Tetrahedron, To=Node
template class Connectivity<4,4,0,3,0>;
template class Connectivity<4,4,1,3,0>;
template class Connectivity<4,4,2,3,0>;
template class Connectivity<4,4,3,3,0>;
template class Connectivity<4,4,4,3,0>;
// mesh 4D  Sub=Pentatope, To=Edge
template class Connectivity<4,4,0,4,0>;
template class Connectivity<4,4,1,4,0>;
template class Connectivity<4,4,2,4,0>;
template class Connectivity<4,4,3,4,0>;
template class Connectivity<4,4,4,4,0>;


/////////////////////////////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle, Tetrahedron, Pentatope  /////////////
/////////////////////////////////////////////////////////////////////////////////////

// mesh 4D  Sub=Node, To=Edge
template class Connectivity<4,4,0,0,1>;
template class Connectivity<4,4,1,0,1>;
template class Connectivity<4,4,2,0,1>;
template class Connectivity<4,4,3,0,1>;
template class Connectivity<4,4,4,0,1>;
// mesh 4D  Sub=Edge, To=Edge
template class Connectivity<4,4,0,1,1>;
template class Connectivity<4,4,1,1,1>;
template class Connectivity<4,4,2,1,1>;
template class Connectivity<4,4,3,1,1>;
template class Connectivity<4,4,4,1,1>;
// mesh 4D  Sub=Triangle, To=Edge
template class Connectivity<4,4,0,2,1>;
template class Connectivity<4,4,1,2,1>;
template class Connectivity<4,4,2,2,1>;
template class Connectivity<4,4,3,2,1>;
template class Connectivity<4,4,4,2,1>;
// mesh 4D  Sub=Tetrahedron, To=Edge
template class Connectivity<4,4,0,3,1>;
template class Connectivity<4,4,1,3,1>;
template class Connectivity<4,4,2,3,1>;
template class Connectivity<4,4,3,3,1>;
template class Connectivity<4,4,4,3,1>;
// mesh 4D  Sub=Pentatope, To=Edge
template class Connectivity<4,4,0,4,1>;
template class Connectivity<4,4,1,4,1>;
template class Connectivity<4,4,2,4,1>;
template class Connectivity<4,4,3,4,1>;
template class Connectivity<4,4,4,4,1>;

/////////////////////////////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle, Tetrahedron, Pentatope  /////////////
/////////////////////////////////////////////////////////////////////////////////////

// mesh 4D  Sub=Node, To=Triangle
template class Connectivity<4,4,0,0,2>;
template class Connectivity<4,4,1,0,2>;
template class Connectivity<4,4,2,0,2>;
template class Connectivity<4,4,3,0,2>;
template class Connectivity<4,4,4,0,2>;
// mesh 4D  Sub=Edge, To=Triangle
template class Connectivity<4,4,0,1,2>;
template class Connectivity<4,4,1,1,2>;
template class Connectivity<4,4,2,1,2>;
template class Connectivity<4,4,3,1,2>;
template class Connectivity<4,4,4,1,2>;
// mesh 4D  Sub=Triangle, To=Triangle
template class Connectivity<4,4,0,2,2>;
template class Connectivity<4,4,1,2,2>;
template class Connectivity<4,4,2,2,2>;
template class Connectivity<4,4,3,2,2>;
template class Connectivity<4,4,4,2,2>;
// mesh 4D  Sub=Tetrahedron, To=Triangle
template class Connectivity<4,4,0,3,2>;
template class Connectivity<4,4,1,3,2>;
template class Connectivity<4,4,2,3,2>;
template class Connectivity<4,4,3,3,2>;
template class Connectivity<4,4,4,3,2>;
// mesh 4D  Sub=Pentatope, To=Triangle
template class Connectivity<4,4,0,4,2>;
template class Connectivity<4,4,1,4,2>;
template class Connectivity<4,4,2,4,2>;
template class Connectivity<4,4,3,4,2>;
template class Connectivity<4,4,4,4,2>;



/////////////////////////////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle, Tetrahedron, Pentatope  /////////////
/////////////////////////////////////////////////////////////////////////////////////

// mesh 4D  Sub=Node, To=Tetrahedron
template class Connectivity<4,4,0,0,3>;
template class Connectivity<4,4,1,0,3>;
template class Connectivity<4,4,2,0,3>;
template class Connectivity<4,4,3,0,3>;
template class Connectivity<4,4,4,0,3>;
// mesh 4D  Sub=Edge, To=Tetrahedron
template class Connectivity<4,4,0,1,3>;
template class Connectivity<4,4,1,1,3>;
template class Connectivity<4,4,2,1,3>;
template class Connectivity<4,4,3,1,3>;
template class Connectivity<4,4,4,1,3>;
// mesh 4D  Sub=Triangle, To=Tetrahedron
template class Connectivity<4,4,0,2,3>;
template class Connectivity<4,4,1,2,3>;
template class Connectivity<4,4,2,2,3>;
template class Connectivity<4,4,3,2,3>;
template class Connectivity<4,4,4,2,3>;
// mesh 4D  Sub=Tetrahedron, To=Tetrahedron
template class Connectivity<4,4,0,3,3>;
template class Connectivity<4,4,1,3,3>;
template class Connectivity<4,4,2,3,3>;
template class Connectivity<4,4,3,3,3>;
template class Connectivity<4,4,4,3,3>;
// mesh 4D  Sub=Pentatope, To=Tetrahedron
template class Connectivity<4,4,0,4,3>;
template class Connectivity<4,4,1,4,3>;
template class Connectivity<4,4,2,4,3>;
template class Connectivity<4,4,3,4,3>;
template class Connectivity<4,4,4,4,3>;

/////////////////////////////////////////////////////////////////////////////////////
/////////////    Entity = Node, Edge, Triangle, Tetrahedron, Pentatope  /////////////
/////////////////////////////////////////////////////////////////////////////////////

// mesh 4D  Sub=Node, To=Pentatope
template class Connectivity<4,4,0,0,4>;
template class Connectivity<4,4,1,0,4>;
template class Connectivity<4,4,2,0,4>;
template class Connectivity<4,4,3,0,4>;
template class Connectivity<4,4,4,0,4>;
// mesh 4D  Sub=Edge, To=Pentatope
template class Connectivity<4,4,0,1,4>;
template class Connectivity<4,4,1,1,4>;
template class Connectivity<4,4,2,1,4>;
template class Connectivity<4,4,3,1,4>;
template class Connectivity<4,4,4,1,4>;
// mesh 4D  Sub=Triangle, To=Pentatope
template class Connectivity<4,4,0,2,4>;
template class Connectivity<4,4,1,2,4>;
template class Connectivity<4,4,2,2,4>;
template class Connectivity<4,4,3,2,4>;
template class Connectivity<4,4,4,2,4>;
// mesh 4D  Sub=Tetrahedron, To=Pentatope
template class Connectivity<4,4,0,3,4>;
template class Connectivity<4,4,1,3,4>;
template class Connectivity<4,4,2,3,4>;
template class Connectivity<4,4,3,3,4>;
template class Connectivity<4,4,4,3,4>;
// mesh 4D  Sub=Pentatope, To=Pentatope
template class Connectivity<4,4,0,4,4>;
template class Connectivity<4,4,1,4,4>;
template class Connectivity<4,4,2,4,4>;
template class Connectivity<4,4,3,4,4>;
template class Connectivity<4,4,4,4,4>;
}