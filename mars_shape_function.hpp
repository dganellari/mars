#ifndef MARS_SHAPE_FUNCTION_HPP
#define MARS_SHAPE_FUNCTION_HPP

#include "mars_simplex.hpp"
#include "mars_base_elementfunctionspace.hpp"
#include "mars_array.hpp"
#include "mars_vector.hpp"
#include "mars_fqpexpressions.hpp"
#include "mars_quadrature_rules.hpp"
#include "mars_operators.hpp"
#include "mars_referencemap.hpp"
#include "mars_shape_function_coefficients.hpp"
#include "mars_tuple_utilities.hpp"

#include "mars_signed_normal.hpp"


namespace mars{




 constexpr Integer simplex_face_sub_entities(const Integer& SimplexDim,const Integer& FaceNumber,const Integer& SubEntityDim, const Integer& SubEntityDimNumber)
 {
  switch(SimplexDim)
  {
   // triangles
   case 2:
   // triangle faces (edges)
   switch(FaceNumber)
   {
    // face 0
    case 0:
       switch(SubEntityDim)
       {
        // face 0, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 1;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 0, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }        
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};

       }
   
    case 1:
       switch(SubEntityDim)
       {
        // face 0, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 2;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 0, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 1;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }        
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};

       }
    case 2:
       switch(SubEntityDim)
       {
        // face 2, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 1;
          case 1: return 2;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 2, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 2;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }        
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};

       }
    default: {assert(0 &&"simplex_face_sub_entities: invalid face number");return -1;};
   }
  // tetrahedrons
  case 3:
  // tetrahedrons faces (triangles)
   switch(FaceNumber)
   {
    // face 0
    case 0:
       switch(SubEntityDim)
       {
        // face 0, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 1;
          case 2: return 2;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 0, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 1;
          case 2: return 3;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 0, triangles
        case 2:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }         
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
       }
   
    case 1:
       switch(SubEntityDim)
       {
        // face 1, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 1;
          case 2: return 3;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 1, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 2;
          case 2: return 4;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 1, triangles
        case 2:
        switch(SubEntityDimNumber)
        {
          case 0: return 1;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }         
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
       }
    case 2:
       switch(SubEntityDim)
       {
        // face 0, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 2;
          case 2: return 3;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 0, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 1;
          case 1: return 2;
          case 2: return 5;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 0, triangles
        case 2:
        switch(SubEntityDimNumber)
        {
          case 0: return 2;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }         
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
       }
    case 3:
        switch(SubEntityDim)
       {
        // face 0, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 1;
          case 1: return 2;
          case 2: return 3;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 0, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 3;
          case 1: return 4;
          case 2: return 5;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 0, triangles
        case 2:
        switch(SubEntityDimNumber)
        {
          case 0: return 3;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }         
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
       }

    default: {assert(0 &&"simplex_face_sub_entities: invalid face number");return -1;};
   }

  // pentatope
  case 4:
  // tetrahedrons faces (tetrahedrons)
   switch(FaceNumber)
   {
    // face 0
    case 0:
       switch(SubEntityDim)
       {
        // face 0, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 1;
          case 2: return 2;
          case 3: return 3;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 0, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 1;
          case 2: return 2;
          case 3: return 4;
          case 4: return 5;
          case 5: return 7;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 0, triangles
        case 2:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 1;
          case 2: return 3;
          case 3: return 6;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 0, tetrahedrons
        case 3:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }        
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
       }

    // face 1
    case 1:
       switch(SubEntityDim)
       {
        // face 1, edges
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 1;
          case 2: return 2;
          case 3: return 4;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 1, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 1;
          case 2: return 3;
          case 3: return 4;
          case 4: return 6;
          case 5: return 8;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 1, triangles
        case 2:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 2;
          case 2: return 4;
          case 3: return 7;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 1, tetrahedrons
        case 3:
        switch(SubEntityDimNumber)
        {
          case 0: return 1;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }        
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
       }
    // face 2
    case 2:
       switch(SubEntityDim)
       {
        // face 2, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 1;
          case 2: return 3;
          case 3: return 4;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 2, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 2;
          case 2: return 3;
          case 3: return 5;
          case 4: return 6;
          case 5: return 9;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 2, triangles
        case 2:
        switch(SubEntityDimNumber)
        {
          case 0: return 1;
          case 1: return 2;
          case 2: return 5;
          case 3: return 8;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 2, tetrahedrons
        case 3:
        switch(SubEntityDimNumber)
        {
          case 0: return 2;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }        
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
       }

    // face 3
    case 3:
       switch(SubEntityDim)
       {
        // face 3, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 0;
          case 1: return 2;
          case 2: return 3;
          case 3: return 4;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 3, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 1;
          case 1: return 2;
          case 2: return 3;
          case 3: return 7;
          case 4: return 8;
          case 5: return 9;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 3, triangles
        case 2:
        switch(SubEntityDimNumber)
        {
          case 0: return 3;
          case 1: return 4;
          case 2: return 5;
          case 3: return 9;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 3, tetrahedrons
        case 3:
        switch(SubEntityDimNumber)
        {
          case 0: return 3;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }        
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
       }

    // face 4
    case 4:
       switch(SubEntityDim)
       {
        // face 4, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 1;
          case 1: return 2;
          case 2: return 3;
          case 3: return 4;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 4, edges
        case 1:
        switch(SubEntityDimNumber)
        {
          case 0: return 4;
          case 1: return 5;
          case 2: return 6;
          case 3: return 7;
          case 4: return 8;
          case 5: return 9;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 4, triangles
        case 2:
        switch(SubEntityDimNumber)
        {
          case 0: return 6;
          case 1: return 7;
          case 2: return 8;
          case 3: return 9;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }  
        // face 4, tetrahedrons
        case 3:
        switch(SubEntityDimNumber)
        {
          case 0: return 4;
          default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
        }        
        default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
       }
    default: {assert(0 &&"simplex_face_sub_entities: invalid face number");return -1;};
   }
  default: {assert(0 &&"simplex_face_sub_entities: invalid simplex dimension");return -1;};
 }
}


















template<typename Elem,typename Operator,Integer FEFamily,Integer Order>
class ReferenceShapeFunctionValue;

template<typename Elem,Integer FEFamily,Integer Order>
class SingleShapeFunctionCoefficientsCollection;







// works only for simplices
template <typename Space>
constexpr Integer trace_surf_n_dofs()
  { 
    using Elem=typename Space::Elem;
    Integer n_dofs=0;
    for(Integer ii=0;ii<Space::entity.size();ii++)
    {
      if(Space::entity[ii]<=Elem::ManifoldDim)
       n_dofs+=Space::dofs_per_entity[ii]* binomial_coefficient(Elem::ManifoldDim,Space::entity[ii]+1);
   }
   // this happens for L2 functions (constant on elements, but discontinuous)
   // we still choose to have 1 dof trace
   if(n_dofs==0)
    n_dofs=1;

   return n_dofs*Space::NComponents;
 }


template<typename FunctionSpace>
constexpr Integer function_space_dofs_per_elem(const Integer Dim,const Integer N)
{   
    std::size_t dofs_per_elem=0;
    if(N>0)
    dofs_per_elem=FunctionSpace::NComponents * 
                                     FunctionSpace::dofs_per_entity[N-1]   * 
                                     binomial_coefficient(Dim+1,FunctionSpace::entity[N-1]+1);  

     switch(N)
     {
      case  -1: return 0;
      case  0: return 0;
      case  1: return dofs_per_elem;
      default: return function_space_dofs_per_elem<FunctionSpace>(Dim,N-1)+dofs_per_elem;
     }
};



template <typename Space>
constexpr auto trace_dofs(const Integer face)
{

 
 constexpr auto n=trace_surf_n_dofs<Space>();
 const auto ManifoldDim=Space::ManifoldDim;
 const auto NComponents=Space::NComponents;
 const auto entity=Space::entity;
 const auto dofs_per_entity= Space::dofs_per_entity;

 Array<Integer, n> dofs;
 Integer cont=0;
 Integer dofs_per_entity_cont=0;
 // loop on all the kind of dofs-entities
 for(Integer ii=0; ii< entity.size(); ii++)
 {
  const auto entityii=entity[ii];
  // consider the ii-th dof-entity only if it smaller than the dimension of the manifolddim
  // TODO FIXME IT SHOULD BE entity[ii]<ManifoldDim
  if(entity[ii]<ManifoldDim)
  {
    const auto dofs_per_entityii=dofs_per_entity[ii];
    const auto binomial_coeff=binomial_coefficient(ManifoldDim,entity[ii]+1);
    // loop on the entityes of entityii dim
    for(Integer jj=0; jj<binomial_coeff; jj++)
      {
        // loop on the dofs of the given entity
        dofs[cont] = NComponents*dofs_per_entityii*simplex_face_sub_entities(ManifoldDim,face,entityii,jj)+function_space_dofs_per_elem<Space>(ManifoldDim,ii);
        cont++;
        for(Integer ss=1;ss<NComponents;ss++)
        { 
          dofs[cont] = dofs[cont-1]+1;
          cont++;
        }

        for(Integer kk=1; kk<dofs_per_entityii; kk++)
          {
               for(Integer ss=0;ss<NComponents;ss++)
                { 
                  dofs[cont] = dofs[cont-1]+1;
                  cont++;
                }
          }
      }
  }
 }
 return dofs;
}


template <typename Space>
class TraceDofs;

template <typename Elem, Integer FEFamily, Integer Order,Integer Continuity, Integer NComponents>
class TraceDofs<ElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>>
{
public:
  using Space=ElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>;
  static constexpr Integer ManifoldDim=Elem::ManifoldDim;
  static constexpr auto dofs()
 {
  constexpr auto n_dofs_per_face=trace_surf_n_dofs<Space>();
  const auto n_faces=binomial_coefficient(ManifoldDim+1,ManifoldDim);
  Array< Array<Integer, n_dofs_per_face>, n_faces> vec;
  for(Integer ii=0;ii<n_faces;ii++)
    vec[ii]=trace_dofs<Space>(ii);
  return vec;
}
};


template<typename MeshT, typename...BaseFunctionSpaces>
class TraceDofs<FunctionSpace<MeshT,BaseFunctionSpaces...>>
{
 public:
 using Elem=typename MeshT::Elem;

 static constexpr auto dofs()
 {
 return std::tuple_cat(std::make_tuple(TraceDofs<ElemFunctionSpace<Elem,BaseFunctionSpaces>>::dofs())...);
 } 
};


template<typename ...Args>
class TraceDofs<FullSpace<MixedSpace<Args...>>>
{
 public:

 static constexpr auto dofs_aux=std::tuple_cat(TraceDofs<Args>::dofs()...);

 static constexpr auto dofs()
 {
 return dofs_aux;
} 
};


template<template<class...> class AuxMixedSpace_,typename ...Args,typename...AuxArgs>
class TraceDofs<FullSpace<MixedSpace<Args...>,AuxMixedSpace_<AuxArgs...>>>
{
 public:

 static constexpr auto dofs_aux=std::tuple_cat(TraceDofs<Args>::dofs()...,TraceDofs<AuxArgs>::dofs()...);

 static constexpr auto dofs()
 {
 return dofs_aux;
} 
};




// template <typename Space>
//  constexpr auto trace_dofs()
//  {
//   const auto ManifoldDim=Space::ManifoldDim;
//   constexpr auto n_dofs_per_face=trace_surf_n_dofs<Space>();
//   const auto n_faces=binomial_coefficient(ManifoldDim+1,ManifoldDim);
//   Array< Array<Integer, n_dofs_per_face>, n_faces> vec;
//   for(Integer ii=0;ii<n_faces;ii++)
//     vec[ii]=trace_dofs<Space>(ii);
//   return vec;
// }


      template<typename Space,typename Operator,typename single_type,
               Integer NQPoints,Integer Dim>
       constexpr auto reference_trace_shape_function_init(const Integer face, const Matrix<Real,NQPoints,Dim>&qp_points)
       {
        
        Array<Real,Dim> qp_point;
        const auto dofs=trace_dofs<Space>(face);
        const auto n_dofs=dofs.size();
        Array<Array<single_type,NQPoints>,n_dofs> v;
        Array<single_type,n_dofs> func;
        
            for(Integer qp=0;qp<NQPoints;qp++)
            {
             qp_point=qp_points.get_row(qp);
             // func=value<Elem,Operator,FEFamily,Order,single_type,Ndofs>(qp_point);
             ReferenceShapeFunctionValue<typename Space::Elem,Operator,Space::FEFamily,Space::Order>::apply(qp_point,func);
             // value<Space::Elem,Operator,Space::FEFamily,Space::Order>(qp_point,func);
              for(Integer n_dof = 0; n_dof < n_dofs; ++n_dof) {
                  const_cast<single_type&>
                  (static_cast<const std::array<single_type,NQPoints>& >
                   ((static_cast<const std::array<Array<single_type,NQPoints>,n_dofs>& >(v())[n_dof] )())[qp])=
                  static_cast<const std::array<single_type,n_dofs>& >(func())[n_dof];
              }
            }
       return v;
      };

// template<Integer N,Integer K>
//  constexpr void combinations_generate_aux(
//             Array<Integer, K> &data,
//             const Integer index, 
//             const Integer i,
//             Array<Array<Integer, K>, binomial_coefficient(N,K)> &combs,
//             Integer &comb_index)
//         {
//             if(index == K) {
//                 for(Integer ii=0;ii<data.size();ii++)
//                     combs[comb_index][ii]=data[ii];
//                 comb_index++;
//                 return;
//             }

//             if(i >= N) {
//                 return;
//             }

//             data[index] = i;

//             combinations_generate_aux<N,K>(data, index+1, i+1, combs, comb_index);
            
//             // current is excluded, replace it with next (Note that
//             // i+1 is passed, but index is not changed)
//             combinations_generate_aux<N,K>(data, index, i+1, combs, comb_index);
//         }

//         template<Integer N,Integer K >
//         constexpr Array<Array<Integer, K>, binomial_coefficient(N,K)> combinations_generate()
//         {
//             Array<Array<Integer, K>, binomial_coefficient(N,K)> combs;
//             Array<Integer, K> data;
//             Integer comb_index = 0;
//             combinations_generate_aux<N,K>(data, 0, 0, combs, comb_index);
//             return combs;
//         }


   template<Integer Dim, Integer ManifoldDim>
inline constexpr auto jacobian_faces()
{
    static_assert(Dim >= ManifoldDim, "Dim must be greater or equal ManifoldDim");
    // static_assert(Npoints == ManifoldDim+1, "Npoints must be equal to ManifoldDim+1");
    Array<Vector<Real, Dim>,ManifoldDim> points;
    constexpr auto n_faces=ManifoldDim+1;
    const auto combs=combinations_generate<ManifoldDim+1,ManifoldDim>(); 
    Matrix<Real, Dim, ManifoldDim-1> J;
    Vector<Matrix<Real, Dim, ManifoldDim-1>,n_faces> Jmat;
    // loop on all the faces
    for(Integer ii=0;ii<n_faces;ii++)
    {
        // take the indices of the reference simplex related to the face ii
        const auto &comb_ii=combs[ii];
        // fill points with the corresponding reference simplex face 
        for(Integer jj=0;jj<ManifoldDim;jj++)
          for(Integer kk=0;kk<ManifoldDim;kk++)
            points[jj][kk]=Simplex<Dim,ManifoldDim>::reference[comb_ii[jj]][kk];
        
        // compute the jacobian of the given face
        for(Integer ii=0;ii<ManifoldDim;ii++)
        {
            Vector<Real, Dim> v0 = points[0];
            for(Integer i = 1; i < ManifoldDim; ++i) {
                const auto &vi = points[i];
                J.col(i-1, vi - v0);
            }
        }
        Jmat[ii]=J;

    }
    return Jmat;
}


template<typename QuadratureRule>
inline constexpr auto reference_face_shape_functions()
{
  using Elem=typename QuadratureRule::Elem;
  constexpr Integer Dim=Elem::Dim;
  constexpr Integer ManifoldDim=Elem::ManifoldDim;
  const auto Jacobian_faces=jacobian_faces<Dim,ManifoldDim>();
  constexpr auto n_faces=Jacobian_faces.size();
  Vector<typename QuadratureRule::qp_points_type, n_faces> qp_points_face;

  for(Integer ii=0;ii<n_faces;ii++)
     {
        qp_points_face[ii]=Jacobian_faces[ii]*QuadratureRule::qp_points_type;
     }

 return qp_points_face;
}





/////////////////////////////////////////////////////////////////////////////

class MaxElementOrder
{
public:
    static constexpr Integer value=2;
};

template<typename Tuple>
class ElementOrder;

template<typename...Spaces>
class ElementOrder<std::tuple<Spaces...>>
{
public: 
    static constexpr auto value=Min(MaxElementOrder::value,Max(QuadratureOrder<IdentityOperator,Spaces>::value...));
};


template<typename Elem,typename Operator,Integer FEFamily,Integer Order>
class ElemPoints;


static constexpr Vector<Matrix<Real,2,1>,3> LinearTrianglePoints({0.0,0.0},
                                                                 {1.0,0.0},
                                                                 {0.0,1.0});

// static constexpr Vector<Integer,3> LinearTriangleOrdering{0,1,2};

static constexpr Vector<Matrix<Real,3,1>,4> LinearTetrahedronPoints({0.0,0.0,0.0},
                                                                    {1.0,0.0,0.0},
                                                                    {0.0,1.0,0.0},
                                                                    {0.0,0.0,1.0});
template<Integer ManifoldDim>
class LinearSimplex
{
 public:
 constexpr static auto points()
 {
  Matrix<Real,ManifoldDim+1,ManifoldDim+1> mat;
  for(Integer i=0;i<ManifoldDim+1;i++)
      {
       for(Integer j=0;j<ManifoldDim+1;j++)
          {
            if(i==j)
                mat(i,j)=1.0;
            else
                mat(i,j)=0.0;
          }       
      }
  return mat;
 }
 constexpr static auto barycenter()
 {
  Vector<Real,ManifoldDim> vec;
   for(Integer i=0;i<ManifoldDim;i++)
      {
                vec(i)=1.0/(ManifoldDim+1);      
      } 
   return vec;
 }
 
};
// static constexpr Vector<Integer,3> LinearTetrahedronOrdering{0,1,2,3};

template<Integer Dim,Integer Ndofs>
constexpr auto vecmat2mat(const Vector<Matrix<Real,Dim,1>,Ndofs>& vec)
{
 Matrix<Real,Ndofs,Dim> mat;
 for(Integer i=0;i<Ndofs;i++)
 {
  for(Integer j=0;j<Dim;j++)
    mat(i,j)=vec[i](j,0);
 }
 return mat;
}



static constexpr Vector<Matrix<Real,2,1>,6> QuadraticTrianglePoints({0.0,0.0},
                                                                    {1.0,0.0},
                                                                    {0.0,1.0},
                                                                    {0.5,0.0},
                                                                    {0.5,0.5},
                                                                    {0.0,0.5});

// static constexpr Vector<Integer,3> QuadraticTriangleOrdering{0,1,2,3,5,4};

static constexpr Vector<Matrix<Real,3,1>,10> QuadraticTetrahedronPoints( {0.0,0.0,0.0},
                                                                         {1.0,0.0,0.0},
                                                                         {0.0,1.0,0.0},
                                                                         {0.0,0.0,1.0},
                                                                         {0.5,0.0,0.0},
                                                                         {0.5,0.5,0.0},
                                                                         {0.0,0.5,0.0},
                                                                         {0.0,0.0,0.5},
                                                                         {0.5,0.0,0.5},
                                                                         {0.0,0.5,0.5});

// static constexpr Vector<Integer,3> QuadraticTetrahedronOrderingg{0,1,2,3,4, 5,4};

template<typename Elem,Integer Order>
class ElemGeometricPoints;

template<Integer Dim>
class ElemGeometricPoints<Simplex<Dim,2>,1>
{
public:
  static constexpr auto points=LinearTrianglePoints;
  using type=decltype(points);
};

template<Integer Dim>
class ElemGeometricPoints<Simplex<Dim,2>,2>
{
public:
  static constexpr auto points=QuadraticTrianglePoints;
  using type=decltype(points);
};


template<Integer Dim>
class ElemGeometricPoints<Simplex<Dim,3>,1>
{
public:
  static constexpr auto points=LinearTetrahedronPoints;
  using type=decltype(points);
};

template<Integer Dim>
class ElemGeometricPoints<Simplex<Dim,3>,2>
{
public:
  static constexpr auto points=QuadraticTetrahedronPoints;
  using type=decltype(points);
};




template<Integer Dim,typename Operator>
class ElemPoints<Simplex<Dim,2>,Operator,LagrangeFE,1>
{
public:
  static constexpr auto points=LinearTrianglePoints;
  using type=decltype(points);
};

template<Integer Dim,typename Operator>
class ElemPoints<Simplex<Dim,2>,Operator,LagrangeFE,2>
{
public:
  static constexpr auto points=QuadraticTrianglePoints;
  using type=decltype(points);
};


template<Integer Dim,typename Operator>
class ElemPoints<Simplex<Dim,3>,Operator,LagrangeFE,1>
{
public:
  static constexpr auto points=LinearTetrahedronPoints;
  using type=decltype(points);
};

template<Integer Dim,typename Operator>
class ElemPoints<Simplex<Dim,3>,Operator,LagrangeFE,2>
{
public:
  static constexpr auto points=QuadraticTetrahedronPoints;
  using type=decltype(points);
};




template<Integer Dim,typename Operator>
class ElemPoints<Simplex<Dim,2>,Operator,RaviartThomasFE,0>
{
public:
  static constexpr auto points=LinearTrianglePoints;
  using type=decltype(points);
};

template<Integer Dim,typename Operator>
class ElemPoints<Simplex<Dim,2>,Operator,RaviartThomasFE,1>
{
public:
  static constexpr auto points=QuadraticTrianglePoints;
  using type=decltype(points);
};


template<Integer Dim,typename Operator>
class ElemPoints<Simplex<Dim,3>,Operator,RaviartThomasFE,0>
{
public:
  static constexpr auto points=LinearTetrahedronPoints;
  using type=decltype(points);
};

template<Integer Dim,typename Operator>
class ElemPoints<Simplex<Dim,3>,Operator,RaviartThomasFE,1>
{
public:
  static constexpr auto points=QuadraticTetrahedronPoints;
  using type=decltype(points);
};


template<typename Elem,typename Operator,Integer FEFamily,Integer Order>
class DofsPoints;


template<typename DofsPoints_>
class DofsPointsType;


template<Integer EntityDim, typename...Ts>
class DofsOrdering;
template<Integer EntityDim, typename...Ts>
class DofsOrdered;

template<Integer EntityDim,typename Elem,Integer FEFamily,Integer Order,Integer Continuity,Integer NComponents>
class DofsOrdering<EntityDim,ElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>>
{
public:

    static constexpr bool must_reorder=false; 
    template<typename T,typename...Ts>
    static auto value(const T&t,const Ts&...ts)
    {
        // std::cout<<"DofsOrdering GENERAL"<<std::endl;


        return t;
    }


};

template<Integer EntityDim,typename Elem,Integer FEFamily,Integer Order,Integer Continuity,Integer NComponents>
class DofsOrdered<EntityDim,ElementFunctionSpace<Elem,FEFamily,Order,Continuity,NComponents>>
{
public:
static constexpr auto entity_points=EntityDim+1;

    
    template<typename Input>
    static auto local_dofs(const Input ordered, const Integer local_dof_count)
    {
     Integer cont=0;
     Array<Integer,NComponents*entity_points> result;
     
     for(Integer ii=0;ii<entity_points;ii++)
     {
        for(Integer cc=0;cc<NComponents;cc++)
           {result[ordered[ii]*NComponents+cc]=local_dof_count+cont;
            cont++;}

     }

     return result;
    }
};   
// case 2d triangolo con 6 dofs per faccia
/// Allora dovrei avere un array che mi contiene
// n0 n1 n2 n3 n4 n5 ma ordinati in base 


template <Integer Dim,Integer ManifoldDim,Integer Continuity, Integer NComponents>
class DofsOrdering<ManifoldDim-1,ElementFunctionSpace<Simplex<Dim,ManifoldDim>,RaviartThomasFE,1,Continuity,NComponents>>
{ 
public:
 using Elem=Simplex<Dim,ManifoldDim>;
 using FunctionSpace=ElementFunctionSpace<Elem,RaviartThomasFE,1,Continuity,NComponents>;
 static constexpr bool must_reorder=true; 
 static constexpr auto entity_dim=ManifoldDim-1;
 // static constexpr auto dofs_per_entity=FunctionSpace::dofs_per_entity[N];
 static constexpr auto entity_points=entity_dim+1;
 static constexpr auto manifold_points=ManifoldDim+1;
 static constexpr auto combinations_nums=ElemEntityCombinations<Elem,entity_dim>::value;
    
    template<std::size_t N>
    static auto value(const std::array<Integer,N>&nodes, Integer entity_iter)
    {
     Integer local_nodes[entity_points];
     std::array<Integer,entity_points> face_nodes;
     Combinations<manifold_points,entity_points>::generate(entity_iter,local_nodes);
     // std::cout<<"DofsOrdering"<<std::endl;

     for(std::size_t i=0;i<entity_points;i++)
          {
            face_nodes[i]=nodes[local_nodes[i]];
            // std::cout<<face_nodes[i]<<" ";
          }
     // std::cout<<std::endl;

     auto ordered_face_nodes=argsort(face_nodes);

     // for(std::size_t i=0;i<ordered_face_nodes.size();i++)
     //      {
     //        std::cout<<ordered_face_nodes[i]<<" ";
     //      }
     // std::cout<<std::endl;




     return ordered_face_nodes;
    }

};


template<typename Elem,typename Operator,Integer FEFamily,Integer Order>
class DofsPointsType<DofsPoints<Elem,Operator,FEFamily,Order>>
{
public:
  using Dofs=DofsPoints<Elem,Operator,FEFamily,Order>;
  using type=decltype(Dofs::points);
  using mapped_type=Vector<Matrix<Real,Elem::ManifoldDim,1> , type::Dim>;
};

static constexpr Vector<Matrix<Real,1,1>,1> Dofs_Simplex1_Trace_Lagrange0(0.5);

template<Integer Dim>
class DofsPoints<Simplex<Dim,2>,TraceOperator, LagrangeFE, 0>
{
public:
 static constexpr auto points=Dofs_Simplex1_Trace_Lagrange0;
 using type=decltype(points);
};

template<Integer Dim,typename Operator>
class DofsPoints<Simplex<Dim,1>,Operator, LagrangeFE, 0>
{
public:
 static constexpr auto points=Dofs_Simplex1_Trace_Lagrange0;
 using type=decltype(points);
};




static constexpr Vector<Matrix<Real,1,1>,2> Dofs_Simplex1_Trace_Lagrange1(0, 1.0);

template<Integer Dim>
class DofsPoints<Simplex<Dim,2>,TraceOperator, LagrangeFE, 1>
{
public:
 static constexpr auto points=Dofs_Simplex1_Trace_Lagrange1;
 using type=decltype(points);
};

template<Integer Dim,typename Operator>
class DofsPoints<Simplex<Dim,1>,Operator, LagrangeFE, 1>
{
public:
 static constexpr auto points=Dofs_Simplex1_Trace_Lagrange1;
 using type=decltype(points);
};







static constexpr Vector<Matrix<Real,1,1>,3> Dofs_Simplex1_Trace_Lagrange2(0, 1.0, 0.5);

template<Integer Dim>
class DofsPoints<Simplex<Dim,2>,TraceOperator, LagrangeFE, 2>
{
public:
 static constexpr auto points=Dofs_Simplex1_Trace_Lagrange2;
 using type=decltype(points);
};

template<Integer Dim,typename Operator>
class DofsPoints<Simplex<Dim,1>,Operator, LagrangeFE, 2>
{
public:
 static constexpr auto points=Dofs_Simplex1_Trace_Lagrange2;
 using type=decltype(points);
};





static constexpr Vector<Matrix<Real,1,1>,4> Dofs_Simplex1_Trace_Lagrange3(0, 1.0, 1.0/3.0, 2.0/3.0);

template<Integer Dim>
class DofsPoints<Simplex<Dim,2>,TraceOperator, LagrangeFE, 3>
{
public:
 static constexpr auto points=Dofs_Simplex1_Trace_Lagrange3;
 using type=decltype(points);
};

template<Integer Dim,typename Operator>
class DofsPoints<Simplex<Dim,1>,Operator, LagrangeFE, 3>
{
public:
 static constexpr auto points=Dofs_Simplex1_Trace_Lagrange3;
 using type=decltype(points);
};



static constexpr Vector<Matrix<Real,1,1>,1> Dofs_Simplex1_Trace_RT0(0.5);

template<Integer Dim>
class DofsPoints<Simplex<Dim,2>,TraceOperator, RaviartThomasFE, 0>
{
public:
 static constexpr auto points=Dofs_Simplex1_Trace_RT0;
 using type=decltype(points);
};





static constexpr Vector<Matrix<Real,1,1>,2> Dofs_Simplex1_Trace_RT1(0.0,1.0);

template<Integer Dim>
class DofsPoints<Simplex<Dim,2>,TraceOperator, RaviartThomasFE, 1>
{
public:
 static constexpr auto points=Dofs_Simplex1_Trace_RT1;
 using type=decltype(points);
};









static constexpr Vector<Matrix<Real,2,1>,1> Dofs_Simplex2_Trace_Lagrange0({1.0/3.0, 1.0/3.0});

template<Integer Dim>
class DofsPoints<Simplex<Dim,3>,TraceOperator, LagrangeFE, 0>
{
public:
 static constexpr auto points=Dofs_Simplex2_Trace_Lagrange0;
 using type=decltype(points);
};

template<Integer Dim,typename Operator>
class DofsPoints<Simplex<Dim,2>,Operator, LagrangeFE, 0>
{
public:
 static constexpr auto points=Dofs_Simplex2_Trace_Lagrange0;
 using type=decltype(points);
};





static constexpr Vector<Matrix<Real,2,1>,3> Dofs_Simplex2_Trace_Lagrange1({0.0,0.0},
                                                                        {1.0,0.0},
                                                                        {0.0,1.0});

template<Integer Dim>
class DofsPoints<Simplex<Dim,3>,TraceOperator, LagrangeFE, 1>
{
public:
 static constexpr auto points=Dofs_Simplex2_Trace_Lagrange1;
 using type=decltype(points);
};

template<Integer Dim,typename Operator>
class DofsPoints<Simplex<Dim,2>,Operator, LagrangeFE, 1>
{
public:
 static constexpr auto points=Dofs_Simplex2_Trace_Lagrange1;
 using type=decltype(points);
};







static constexpr Vector<Matrix<Real,2,1>,6> Dofs_Simplex2_Trace_Lagrange2({0.0,0.0},
                                                                        {1.0,0.0},
                                                                        {0.0,1.0},
                                                                        {0.5,0.0},
                                                                        {0.0,0.5},
                                                                        {0.5,0.5});

template<Integer Dim>
class DofsPoints<Simplex<Dim,3>,TraceOperator, LagrangeFE, 2>
{
public:
 static constexpr auto points=Dofs_Simplex2_Trace_Lagrange1;
 using type=decltype(points);
};


template<Integer Dim,typename Operator>
class DofsPoints<Simplex<Dim,2>,Operator, LagrangeFE, 2>
{
public:
 static constexpr auto points=Dofs_Simplex2_Trace_Lagrange1;
 using type=decltype(points);
};





static constexpr Vector<Matrix<Real,2,1>,10> Dofs_Simplex2_Trace_Lagrange3({0.0,0.0},
                                                                         {1.0,0.0},
                                                                         {0.0,1.0},
                                                                         {1.0/3.0,0.0},
                                                                         {2.0/3.0,0.0},
                                                                         {0.0,1.0/3.0},
                                                                         {0.0,2.0/3.0},
                                                                         {2.0/3.0,1.0/3.0},
                                                                         {1.0/3.0,2.0/3.0},
                                                                         {1.0/3.0,1.0/3.0});

template<Integer Dim>
class DofsPoints<Simplex<Dim,3>,TraceOperator, LagrangeFE, 3>
{
public:
 static constexpr auto points=Dofs_Simplex2_Trace_Lagrange3;
 using type=decltype(points);
};

template<Integer Dim,typename Operator>
class DofsPoints<Simplex<Dim,2>,Operator, LagrangeFE, 3>
{
public:
 static constexpr auto points=Dofs_Simplex2_Trace_Lagrange3;
 using type=decltype(points);
};









static constexpr Vector<Matrix<Real,2,1>,1> Dofs_Simplex2_Trace_RT0({1.0/3.0, 1.0/3.0});

static constexpr Vector<Matrix<Real,2,1>,3> Dofs_Simplex2_Operator_RT0({0.5, 0.0},  
                                                                       {0.0, 0.5},
                                                                       {0.5, 0.5});


template<Integer Dim>
class DofsPoints<Simplex<Dim,3>,TraceOperator, RaviartThomasFE, 0>
{
public:
 static constexpr auto points=Dofs_Simplex2_Trace_RT0;
 using type=decltype(points);
};

template<Integer Dim,typename Operator>
class DofsPoints<Simplex<Dim,2>,Operator, RaviartThomasFE, 0>
{
public:
 static constexpr auto points=Dofs_Simplex2_Operator_RT0;
 using type=decltype(points);
};




static constexpr Vector<Matrix<Real,2,1>,3> Dofs_Simplex2_Trace_RT1({0.0,0.0},
                                                                  {1.0,0.0},
                                                                  {0.0,1.0});

static constexpr Vector<Matrix<Real,2,1>,8> Dofs_Simplex2_Operator_RT1({0.0,0.0},
                                                                       {1.0,0.0},
                                                                       {0.0,0.0},
                                                                       {0.0,1.0},
                                                                       {1.0,0.0},
                                                                       {0.0,1.0},
                                                                       {1.0/3.0, 1.0/3.0},
                                                                       {1.0/3.0, 1.0/3.0});
template<Integer Dim>
class DofsPoints<Simplex<Dim,3>,TraceOperator, RaviartThomasFE, 1>
{
public:
 static constexpr auto points=Dofs_Simplex2_Trace_RT1;
 using type=decltype(points);
};


template<Integer Dim,typename Operator>
class DofsPoints<Simplex<Dim,2>,Operator, RaviartThomasFE, 1>
{
public:
 static constexpr auto points=Dofs_Simplex2_Operator_RT1;
 using type=decltype(points);
};




template<Integer ManifoldDim,Integer Order>
constexpr auto NumberOfLagrangianSimplexDofs()
{
    Real N=1;

    // dim is Integer, so divide N by dim: (N/dim) is Real
    // then (dim+Order) is again integer, but we multiply it by a real
    // do not do (dim+Order)/dim which is an integer division!
    for(Integer dim=1;dim<ManifoldDim+1;dim++)
        N= (dim+Order) * (N/dim);
    return static_cast<Integer>(N);
}

template<Integer ManifoldDim >
constexpr auto NumberOfRaviartThomas0SimplexDofs()
{
    return static_cast<Integer>(ManifoldDim+1);
}

template<Integer Order >
constexpr auto NumberOfRaviartThomasSimplex2Dofs()
{
    return static_cast<Integer>((Order+1)*(Order+3));
}

template<Integer Order >
constexpr auto NumberOfRaviartThomasSimplex3Dofs()
{
    return static_cast<Integer>(0.5*(Order+1)*(Order+2)*(Order+4));
}

template<Integer ManifoldDim,Integer Order >
constexpr auto NumberOfRaviartThomasSimplexDofs()
{
    if(Order==0)
      return NumberOfRaviartThomas0SimplexDofs<ManifoldDim>();
    else if(ManifoldDim==2)
      return NumberOfRaviartThomasSimplex2Dofs<Order>();
    else if(ManifoldDim==3)
      return NumberOfRaviartThomasSimplex3Dofs<Order>();
    else
      return static_cast<Integer>(-1);
   
}

template<Integer Dim, Integer ManifoldDim>
class SimplexRT1Coefficient
{
public:

    using OutputRT0=typename ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 0>::Output;
    static constexpr auto barycenter=LinearSimplex<ManifoldDim>::barycenter();
    
    static constexpr auto value()
    {
        OutputRT0 rt0barycenter;
        ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 0>::apply(barycenter,rt0barycenter);
        Matrix<Real,ManifoldDim,ManifoldDim> mat;

        for(Integer i=0;i<ManifoldDim;i++)
           for(Integer j=0;j<ManifoldDim;j++)
               mat(i,j)=rt0barycenter[j](i,0);
        return mat;
    }
};
//////////////////////////////////////////////////////
////////        SIMPLEX - IDENTITY - P0       ////////
//////////////////////////////////////////////////////
template<Integer Dim,Integer ManifoldDim>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, LagrangeFE, 0>
{
 public: 
  using Output=Vector<Matrix<Real, 1, 1>,1>;
constexpr inline static void 
apply(const Vector<Real,ManifoldDim>& point, Output & func, Real alpha=1.0)
{
    Output func2((alpha));       
    func=func2;
}
};


//////////////////////////////////////////////////////
////////        SIMPLEX - IDENTITY - P1       ////////
//////////////////////////////////////////////////////


template<Integer Dim,Integer ManifoldDim>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, LagrangeFE, 1>
{
 public: 
  using Output=Vector<Matrix<Real, 1, 1>,ManifoldDim+1>;
constexpr inline static void 
apply(const Vector<Real,ManifoldDim>& point, Output & func, Real alpha=1.0)
{
    func[0](0,0)=1 * alpha;
    for(Integer i=0;i<ManifoldDim;i++)
        func[0](0,0)-=point[i] * alpha;

    for(Integer i=0;i<ManifoldDim;i++)
        func[i+1](0,0)=point[i] * alpha;

}
};




// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,1>, IdentityOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,2>;
// constexpr inline static void 
// apply(const Vector<Real,1>& point, Output & func)
// {
//     Output func2((1. - point[0]), // 1 in (0)
//                   point[0]);      // 1 in (1)
//     func=func2;
// }
// };


// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,3>;
// constexpr inline static void 
// apply(const Vector<Real,2>& point, Output & func)
// {
//     Output func2((1. - point[0] - point[1]), // 1 in (0,0)
//                   point[0],                  // 1 in (1,0)
//                   point[1]);                 // 1 in (0,1)
//     func=func2;
// }
// };


// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,3>, IdentityOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,4>;
// constexpr inline static void 
// apply(const Vector<Real,3>& point, Output & func)
// {
//     Output func2((1. - point[0] - point[1]- point[2]), // 1 in (0,0,0)
//                   point[0],                            // 1 in (1,0,0)
//                   point[1],                            // 1 in (0,1,0)
//                   point[2]);                           // 1 in (1,1,1)           
//     func=func2;
// }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,4>, IdentityOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,5>;
// constexpr inline static void 
// apply(const Vector<Real,4>& point, Output & func)
// {
//     Output func2((1. - point[0] - point[1]- point[2] - point[3]), // 1 in (0,0,0,0)
//                   point[0],                                       // 1 in (1,0,0,0)
//                   point[1],                                       // 1 in (0,1,0,0)
//                   point[2],                                       // 1 in (0,0,1,0)  
//                   point[3]);                                      // 1 in (1,1,1,1)       
//     func=func2;
// }
// };

//////////////////////////////////////////////////////
////////        SIMPLEX - GRADIENT - P1       ////////
//////////////////////////////////////////////////////

template<Integer Dim,Integer ManifoldDim>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, GradientOperator, LagrangeFE, 1>
{
 public: 
  using Output=Vector<Matrix<Real, ManifoldDim, 1>,ManifoldDim+1>;
 constexpr inline static void 
 apply(const Vector<Real,ManifoldDim>& point, Output & func, Real alpha=1.0)
  {

      for(Integer i=0;i<ManifoldDim;i++)
        func[0](i,0)=-1 * alpha;

      for(Integer i=0;i<ManifoldDim;i++)
        for(Integer j=0;j<ManifoldDim;j++)
            if(i==j)
                func[i+1](i,0)= 1.0 * alpha;
            else
                func[i+1](j,0)=0.0;
  }
};


// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,1>, GradientOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,2>;
//  constexpr inline static void 
//  apply(const Vector<Real,1>& point, Output & func)
//   {
//       Output func2({-1},
//                    {+1});
//       func=func2;
//   }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,2>, GradientOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 2, 1>,3>;
//  constexpr inline static void 
//  apply(const Vector<Real,2>& point, Output & func)
//   {
//       Output func2({-1,-1},
//                    {+1, 0},
//                    { 0,+1});
//       func=func2;
//   }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,3>, GradientOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 3, 1>,4>;
//  constexpr inline static void 
//  apply(const Vector<Real,3>& point, Output & func)
//   {
//       Output func2({-1,-1,-1},
//                    {+1, 0, 0},
//                    { 0,+1, 0},
//                    { 0, 0,+1});
//       func=func2;
//   }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,4>, GradientOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 4, 1>,5>;
//  constexpr inline static void 
//  apply(const Vector<Real,4>& point, Output & func)
//   {
//       Output func2({-1,-1,-1,-1},
//                    {+1, 0, 0, 0},
//                    { 0,+1, 0, 0},
//                    { 0, 0,+1, 0},
//                    { 0, 0, 0,+1});
//       func=func2;
//   }
// };


//////////////////////////////////////////////////////
////////        SIMPLEX - IDENTITY - P2       ////////
//////////////////////////////////////////////////////

template<Integer Dim,Integer ManifoldDim>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, LagrangeFE, 2>
{
 public:
  static constexpr Integer Ndofs=NumberOfLagrangianSimplexDofs<ManifoldDim,2>();
  using Output=Vector<Matrix<Real, 1, 1>,Ndofs>;
constexpr inline static void 
apply(const Vector<Real,ManifoldDim>& point, Output & func, Real alpha=1.0)
{

    Real xend=1.0 * alpha;
    Integer cont=0;

    for(Integer i=0;i<ManifoldDim;i++)
       xend-=point[i] * alpha;

    func[0](0,0)=2*xend*(xend-0.5) * alpha;
    for(Integer i=0;i<ManifoldDim;i++)
       func[i+1](0,0)=2*point[i]*(point[i]-0.5) * alpha;

    for(Integer i=0;i<ManifoldDim;i++)
       func[ManifoldDim+1+i](0,0)=4*point[i]*xend * alpha;

    for(Integer i=0;i<ManifoldDim-1;i++)
        for(Integer j=i+1;j<ManifoldDim;j++)
           {
            func[2*ManifoldDim+1+cont](0,0)=4*point[i]*point[j] * alpha;
            cont++;
           }

}
};


// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,1>, IdentityOperator, LagrangeFE, 2>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,3>;
// constexpr inline static void 
// apply(const Vector<Real,1>& point, Output & func)
// {
//     Output func2(2*(0.5-point[0])*(1-point[0]), // 1 in (0.0)
//                  4*(point[0])*(1-point[0]),     // 1 in (1.0)
//                  2*(-point[0])*(0.5-point[0])); // 1 in (0.0)
//     func=func2;
// }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 2>
// {
//  public:
//   using Output=Vector<Matrix<Real, 1, 1>,6>;
// constexpr inline static void 
// apply(const Vector<Real,2>& point, Output & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
//     const Real zeta = 1. - xi - eta;
//     Output func2(2.*zeta*(zeta-0.5), // 1 in (0,0)
//                  2.*xi*(xi-0.5),     // 1 in (1,0)
//                  2.*eta*(eta-0.5),   // 1 in (0,1)
//                  4.*xi*zeta,         // 1 in (0.5,0)
//                  4.*eta*zeta,        // 1 in (0,0.5)
//                  4.*xi*eta);         // 1 in (0.5,0.5)
//     func=func2;
// }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,3>, IdentityOperator, LagrangeFE, 2>
// {
//  public:
//   using Output=Vector<Matrix<Real, 1, 1>,10>;
// constexpr inline static void 
// apply(const Vector<Real,3>& point, Output & func)
// {
//     const auto& x0=point[0];
//     const auto& x1=point[1];
//     const auto& x2=point[2];
//     const Real x3 = 1. - x0 - x1 - x2;
//     Output func2(2.*x3*(x3-0.5),     // 1 in (0,0,0)
//                  2.*x0*(x0-0.5),     // 1 in (1,0,0)
//                  2.*x1*(x1-0.5),     // 1 in (0,1,0)
//                  2.*x2*(x2-0.5),     // 1 in (0,0,1)
//                  4.*x0*x3,           // 1 in (0.5,0,0)
//                  4.*x1*x3,           // 1 in (0,0.5,0)
//                  4.*x2*x3,           // 1 in (0,0,0.5)
//                  4.*x0*x1,           // 1 in (0,0.5)
//                  4.*x0*x2,           // 1 in (0,0.5)
//                  4.*x1*x2);          // 1 in (0.5,0.5)
//     func=func2;
// }
// };











//////////////////////////////////////////////////////
////////        SIMPLEX - GRADIENT - P2       ////////
//////////////////////////////////////////////////////
// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,1>, GradientOperator, LagrangeFE, 2>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,3>;
// constexpr inline static void 
// apply(const Vector<Real,1>& point, Output & func)
// {

//     Output func2({4 * point[0] - 3},           
//                  {4 - 8 * point[0]}, 
//                  {4 * point[0] - 1});
//     func=func2;
// }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,2>, GradientOperator, LagrangeFE, 2>
// {
// public:
//   using Output=Vector<Matrix<Real, 2, 1>,6>;
// constexpr inline static void 
// apply(const Vector<Real,2>& point, Output & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
//     const Real zeta = 1. - xi - eta;
//     const Real dxi_dxi    = 1.;
//     const Real deta_deta  = 1.;
//     const Real dzeta_dxi  = -1.;
//     const Real dzeta_deta = -1.;
//     Output func2({2.*zeta*dzeta_dxi  + 2*dzeta_dxi *(zeta-0.5), 2.*zeta*dzeta_deta + 2*dzeta_deta*(zeta-0.5)},
//                  {2.*xi*dxi_dxi  + 2.*dxi_dxi *(xi-0.5),0},
//                  {0,2.*eta*deta_deta + 2.*deta_deta*(eta-0.5)},
//                  {4.*(dxi_dxi * zeta + xi * dzeta_dxi),4. * xi * dzeta_deta},
//                  {4. * eta * dzeta_dxi,4. * (deta_deta * zeta + eta * dzeta_deta)},
//                  {4 * eta, 4 * xi});
    
//     func=func2;
//   }
// };



// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,3>, GradientOperator, LagrangeFE, 2>
// {
// public:
//   using Output=Vector<Matrix<Real, 3, 1>,10>;
// constexpr inline static void 
// apply(const Vector<Real,3>& point, Output & func)
// {
//     const auto& x0=point[0];
//     const auto& x1=point[1];
//     const auto& x2=point[1];

//     Output func2({4*x0 + 4*x1 + 4*x2 - 3, 4*x0 + 4*x1 + 4*x2 - 3, 4*x0 + 4*x1 + 4*x2 - 3},
//                  {4*x0 - 1,0,0},
//                  {0,4*x1 - 1,0},
//                  {0,0,4*x2 - 1},
//                  {4 - 4*x1 - 4*x2 - 8*x0,-4*x0,-4*x0},
//                  {-4*x1, 4 - 8*x1 - 4*x2 - 4*x0,-4*x1},
//                  {-4*x2,-4*x2, 4 - 4*x1 - 8*x2 - 4*x0},
//                  {4*x1,4*x0,0},
//                  {4*x2,0,4*x0},
//                  {0,4*x2,4*x1});
    
//     func=func2;
//   }
// };


template<Integer Dim,Integer ManifoldDim>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, GradientOperator, LagrangeFE, 2>
{
public:
 static constexpr Integer Ndofs=NumberOfLagrangianSimplexDofs<ManifoldDim,2>();
 using Output=Vector<Matrix<Real, ManifoldDim, 1>,Ndofs>;
 constexpr inline static void 
 apply(const Vector<Real,ManifoldDim>& point, Output & func, Real alpha=1.0)
 {
   // set all elements to zero
   for(Integer i=0;i<Ndofs;i++)
    for(Integer j=0;j<ManifoldDim;j++)
        func[i](j,0)=0.0;

    Real tmp= (- 3.0) * alpha;
    for(Integer i=0;i<ManifoldDim;i++)
        tmp+= (4.0 * point[i]) * alpha;

    // set the 0-row=[4 x0 + 4 x1 + 4 x2...-3,4 x0 + 4 x1 + 4 x2...-3... ]
    for(Integer j=0;j<ManifoldDim;j++)
        func[0](j,0)= tmp * alpha;

    // set row 1...ManifoldDim
    for(Integer i=0;i<ManifoldDim;i++)
        func[i+1](i,0)=(4*point[i]-1) * alpha;    

    tmp=(4.0) * alpha;
    for(Integer i=0;i<ManifoldDim;i++)
        tmp-= (4*point[i]) * alpha;

    // set row ManifoldDim+1...2*ManifoldDim+1
    for(Integer i=0;i<ManifoldDim;i++)
        for(Integer j=0;j<ManifoldDim;j++)
            func[ManifoldDim+1+i](j,0)=(-4*point[i]) * alpha; 

    for(Integer i=0;i<ManifoldDim;i++)
        func[ManifoldDim+1+i](i,0)= (tmp-4*point[i]) * alpha;  

    Integer cont=0;
    for(Integer i=0;i<ManifoldDim-1;i++)
        {for(Integer j=i+1;j<ManifoldDim;j++)
           {
            for(Integer k=0;k<ManifoldDim;k++)
            {
                if(k==j)
                    func[2*ManifoldDim+1+cont](k,0)= (4*point[i]) * alpha;
                else if(k==i)
                    func[2*ManifoldDim+1+cont](k,0)= (4*point[j]) * alpha;
                else
                    func[2*ManifoldDim+1+cont](k,0)=0.0;
            }
            cont++;
           }
         
        }

  }
};









//////////////////////////////////////////////////////
////////        SEGMENT - TRACE - P_n         ////////
//////////////////////////////////////////////////////
template<Integer Dim,Integer Order>
class ReferenceShapeFunctionValue<Simplex<Dim,1>, TraceOperator, LagrangeFE, Order>
{
 public: 
  using Output=Vector<Matrix<Real, 1, 1>,Order>;
constexpr inline static void 

apply(const Vector<Real,0>& point, Output & func, Real alpha=1.0)
// check if Vector<Real,0> or Vector<Real,1>
// apply(const Vector<Real,1>& point, Output & func)
  {
      Output func2{1.0 * alpha};
      func=func2;
  }
};

//////////////////////////////////////////////////////
////////        SIMPLEX - TRACE - P1          ////////
//////////////////////////////////////////////////////
template<Integer Dim,Integer ManifoldDim,Integer Order>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, TraceOperator, LagrangeFE, Order>
{
 public: 
 static constexpr Integer Ndofs=NumberOfLagrangianSimplexDofs<ManifoldDim,ManifoldDim-1>();
 using Output=Vector<Matrix<Real, 1, 1>,Ndofs>;
 constexpr inline static void 
 apply(const Vector<Real,ManifoldDim-1>& point, Output & func, Real alpha=1.0)
  {
      ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim-1>, IdentityOperator, LagrangeFE, Order>::apply(point,func,alpha);
  }
};

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,2>, TraceOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,2>;
// constexpr inline static void 
// apply(const Vector<Real,1>& point, Output & func)
//   {
//       ReferenceShapeFunctionValue<Simplex<Dim,1>, IdentityOperator, LagrangeFE, 1>::apply(point,func);
//   }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,3>, TraceOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,3>;
// constexpr inline static void 
// apply(const Vector<Real,2>& point, Output & func)
//   {
//       ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 1>::apply(point,func);
//   }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,4>, TraceOperator, LagrangeFE, 1>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,4>;
// constexpr inline static void 
// apply(const Vector<Real,3>& point, Output & func)
//   {
//       ReferenceShapeFunctionValue<Simplex<Dim,3>, IdentityOperator, LagrangeFE, 1>::apply(point,func);
//   }
// };






//////////////////////////////////////////////////////
////////        SIMPLEX - TRACE - P2          ////////
//////////////////////////////////////////////////////
// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,2>, TraceOperator, LagrangeFE, 2>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,3>;
// constexpr inline static void 
// apply(const Vector<Real,1>& point, Output & func)
//   {
//       ReferenceShapeFunctionValue<Simplex<Dim,1>, IdentityOperator, LagrangeFE, 2>::apply(point,func);
//   }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,3>, TraceOperator, LagrangeFE, 2>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,6>;
// constexpr inline static void 
// apply(const Vector<Real,2>& point, Output & func)
//   {
//       ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 2>::apply(point,func);
//   }
// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,4>, TraceOperator, LagrangeFE, 2>
// {
//  public: 
//   using Output=Vector<Matrix<Real, 1, 1>,10>;
// constexpr inline static void 
// apply(const Vector<Real,3>& point, Output & func)
//   {
//       ReferenceShapeFunctionValue<Simplex<Dim,3>, IdentityOperator, LagrangeFE, 2>::apply(point,func);
//   }
// };



















//////////////////////////////////////////////////////////
////////        SIMPLEX - IDENTITY - RT0          ////////
//////////////////////////////////////////////////////////


template<Integer Dim,Integer ManifoldDim>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 0>
{
public:
 static constexpr Integer Ndofs=NumberOfRaviartThomasSimplexDofs<ManifoldDim,0>();
 using Output=Vector<Matrix<Real, ManifoldDim, 1>,Ndofs>;

  // using Output=Vector<Matrix<Real, ManifoldDim, 1>,ManifoldDim+1>;


 constexpr inline static void 
 apply(const Vector<Real,ManifoldDim>& point, Output & func)
{
    // const auto& xi=point[0];
    // const auto& eta=point[1];
    // Output func2{{xi,eta-1},{xi-1,eta},{xi,eta}};
    // func=func2;
    // loop on shape functions
    Real one_frac_vol_surface=1.0/reference_simplex_volume<ManifoldDim-1>();

    for(Integer i=0;i<ManifoldDim+1;i++)
        // loop on dimension
        for(Integer j=0;j<ManifoldDim;j++)
           func[i](j,0)=point[j] * one_frac_vol_surface;

    for(Integer i=0;i<ManifoldDim;i++)
        // loop on dimension
        for(Integer j=0;j<ManifoldDim;j++)
            if(j==ManifoldDim-1-i)
               func[i](j,0)-=(1 * one_frac_vol_surface); 

    // func=SimplexRT0Identity(point);
}

};

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, RaviartThomasFE, 0>
// {
// public:
//   using Output=Vector<Matrix<Real, 2, 1>,3>;
//  constexpr inline static void 
//  apply(const Vector<Real,2>& point, Output & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
//     Output func2{{xi,eta-1},{xi-1,eta},{xi,eta}};
//     func=func2;
// }

// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,3>, IdentityOperator, RaviartThomasFE, 0>
// {
// public:
//   using Output=Vector<Matrix<Real, 3, 1>,4>;
//  constexpr inline static void 
//  apply(const Vector<Real,3>& point, Output & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
//     const auto& zeta=point[2];
//     Output func2{{xi,eta,zeta-1},{xi,eta-1,zeta},{xi-1,eta,zeta}, {xi,eta,zeta}};
//     func=func2;
// }

// };

// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,4>, IdentityOperator, RaviartThomasFE, 0>
// {
// public:
//   using Output=Vector<Matrix<Real, 4, 1>,5>;
//  constexpr inline static void 
//  apply(const Vector<Real,3>& point, Output & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
//     const auto& zeta=point[2];
//     const auto& theta=point[3];
//     Output func2{{xi,eta,zeta,theta-1},
//                  {xi,eta,zeta-1,theta},
//                  {xi,eta-1,zeta,theta},
//                  {xi-1,eta,zeta,theta}, 
//                  {xi,eta,zeta,theta}};
//     func=func2;
// }

// };





////////////////////////////////////////////////////////////
////////        SIMPLEX - DIVERGENCE - RT0          ////////
////////////////////////////////////////////////////////////
// template<Integer Dim>
// class ReferenceShapeFunctionValue<Simplex<Dim,2>, DivergenceOperator, RaviartThomasFE, 0>
// {
// public:
//   using Output=Vector<Matrix<Real, 1, 1>,3>;
// constexpr inline static void 
// apply(const Vector<Real,2>& point, Output & func)
// {
//     Output func2{{2},{2},{2}};
//     func=func2;
// }
// };
template<Integer Dim,Integer ManifoldDim>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, DivergenceOperator, RaviartThomasFE, 0>
{
public:
 static constexpr Integer Ndofs=NumberOfRaviartThomasSimplexDofs<ManifoldDim,0>();
 using Output=Vector<Matrix<Real, 1, 1>,Ndofs>;
  // using Output=Vector<Matrix<Real, 1, 1>,ManifoldDim+1>;
constexpr inline static void 
apply(const Vector<Real,ManifoldDim>& point, Output & func)
{
    Real one_frac_vol_surface=1.0/reference_simplex_volume<ManifoldDim-1>();
    for(Integer i=0;i<ManifoldDim+1;i++)
        func[i](0,0)=ManifoldDim * one_frac_vol_surface;
}
};


template<Integer Dim,Integer ManifoldDim>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 1>
{
public:
  static constexpr Integer Order=1;
  // static constexpr Integer ManifoldDim=2;
  static constexpr Integer Ndofs=NumberOfRaviartThomasSimplexDofs<ManifoldDim,1>();
  static constexpr Integer RT0Ndofs=NumberOfRaviartThomasSimplexDofs<ManifoldDim,0>();
  static constexpr Integer PnNdofs=NumberOfLagrangianSimplexDofs<ManifoldDim,Order>();
  using Output=Vector<Matrix<Real, ManifoldDim, 1>,Ndofs>;
  // using Output=Vector<Matrix<Real, 1, 1>,ManifoldDim+1>;
  // using Output=Vector<Matrix<Real, 2, 1>,8>;
  using OutputRT0=typename ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 0>::Output;
  using OutputP1=typename ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, LagrangeFE, 1>::Output;
  // static constexpr auto barycenter=LinearSimplex<ManifoldDim>::barycenter();
  static constexpr auto mat_coeffs=SimplexRT1Coefficient<Dim,ManifoldDim>::value();

constexpr inline static void 
apply (const Vector<Real,ManifoldDim>& point, Output & func)
{


    OutputRT0 rt0;//, rt0barycenter;
    OutputP1 p1;
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 0>::apply(point,rt0);
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, LagrangeFE, 1>::apply(point,p1);
    
    // // const auto& rt0_0x=rt0[0](0,0);
    // // const auto& rt0_0y=rt0[0](1,0);
    // // const auto& rt0_1x=rt0[1](0,0);
    // // const auto& rt0_1y=rt0[1](1,0);
    // // const auto& rt0_2x=rt0[2](0,0);
    // // const auto& rt0_2y=rt0[2](1,0);

    // // const auto& p1_0=p1[0](0,0);
    // // const auto& p1_1=p1[1](0,0);
    // // const auto& p1_2=p1[2](0,0);

    Integer cont=0;
    // loop on RT0 functions
    for(Integer i=0;i<RT0Ndofs;i++)
      // loop on dimension
      for(Integer j=0;j<ManifoldDim;j++)
        {
        // loop on dimension
        for(Integer k=0;k<ManifoldDim;k++)
           func[cont](k,0)=rt0[i](k); 
         cont++;
         }
    cont=0;

    for(Integer i=0;i<RT0Ndofs;i++)
      // loop on dimension
      for(Integer j=0;j<PnNdofs;j++)
        {
        // loop on dimension
        if(j!=ManifoldDim-i)
        {for(Integer k=0;k<ManifoldDim;k++)
           func[cont](k,0)=func[cont](k,0)*p1[j](0,0); 
         cont++;
         }
         }



    // // the i-th column represents the coefficients of the i-th momenta shape function
    // // auto inv_mat=inverse(mat);
    for(Integer i=0;i<ManifoldDim;i++)
        for(Integer j=0;j<ManifoldDim;j++)
            {
                func[cont+i](j,0)=0;
                for(Integer k=0;k<ManifoldDim;k++)
                   func[cont+i](j,0)+=mat_coeffs(k,i)*rt0[ManifoldDim-k](j,0)*p1[k](0,0);
            }


    // Vector<Real,ManifoldDim> coeff;
    

    // mat(0,0)=rt0barycenter[0](0,0);
    // mat(1,0)=rt0barycenter[0](1,0);
    // mat(0,1)=rt0barycenter[1](0,0);
    // mat(1,1)=rt0barycenter[1](1,0);

        // for(Integer i=0;i<ManifoldDim;i++)
    //     // loop on dimension
    //     for(Integer j=0;j<ManifoldDim;j++)
    //         if(j==ManifoldDim-1-i)
    //            func[i](j,0)-=1; 



    // Output func2{ {rt0_0x * p1_0, rt0_0y * p1_0  },
    //               {rt0_0x * p1_1, rt0_0y * p1_1  },
    //               {rt0_1x * p1_0, rt0_1y * p1_2  },
    //               {rt0_1x * p1_0, rt0_1y * p1_2  },
    //               {rt0_2x * p1_1, rt0_2y * p1_2  },
    //               {rt0_2x * p1_1, rt0_2y * p1_2  },
    //               {rt0_0x * p1_2, rt0_0y * p1_2  },
    //               {rt0_1x * p1_1, rt0_1y * p1_1  }};






    // const auto& xi=point[0];
    // const auto& eta=point[1];
    // const Real a=0.333333333333333;
    // const Real b=-0.666666666666666;
    // func[cont](0,0)=(a*(xi)*(1. - xi - eta) - a*(xi-1)*xi);
    // func[cont](1,0)=( a*(eta)*(1. - xi - eta) - a*(eta)*xi);
    // func[cont+1](0,0)=(- (b*(xi)*(1. - xi - eta) - a*(xi-1)*xi));
    // func[cont+1](1,0)=(- (b*(eta)*(1. - xi - eta) - a*(eta)*xi));
    // Output func2{
    //     {(1. - xi - eta)*xi,(1. - xi - eta)*(eta-1)},   // 0 in (1,0), (0,1), non-zero normal on edge0
    //     {xi*xi,xi*(eta-1)},                             // 0 in (0,0), (0,1), non-zero normal on edge0
    //     {(1. - xi - eta)*(xi-1),(1. - xi - eta)*(eta)}, // 0 in (1,0), (0,1), non-zero normal on edge1
    //     {eta*(xi-1),eta*eta},                           // 0 in (0,0), (1,0), non-zero normal on edge1
    //     {xi*xi,xi*eta},                                 // 0 in (0,0), (0,1), non-zero normal on edge2
    //     {eta*xi,eta*eta},                               // 0 in (0,0), (1,0), non-zero normal on edge2
    //     // {eta*xi,eta*(eta-1)},                           // normal 0 on all edges, element-dof
    //     // {xi*(xi-1),xi*eta}                              // normal 0 on all edges, element-dof,
    //     { (a*(xi)*(1. - xi - eta) - a*(xi-1)*xi),  ( a*(eta)*(1. - xi - eta) - a*(eta)*xi)},
    //     {- (b*(xi)*(1. - xi - eta) - a*(xi-1)*xi), - (b*(eta)*(1. - xi - eta) - a*(eta)*xi)}    

    // };
    // func=func2;
}







inline static void 
apply2 (const Vector<Real,ManifoldDim> point)
{
    OutputRT0 rt0;
    OutputP1 p1;
    Output func;
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 0>::apply(point,rt0);
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, LagrangeFE, 1>::apply(point,p1);

    Integer cont=0;
    // loop on RT0 functions
    for(Integer i=0;i<RT0Ndofs;i++)
      // loop on dimension
      for(Integer j=0;j<ManifoldDim;j++)
        {
        // loop on dimension
        for(Integer k=0;k<ManifoldDim;k++)
           func[cont](k,0)=rt0[i](k); 
         // std::cout<<cont<<"   "<<std::endl;
         // std::cout<<func[cont]<<std::endl;  
         cont++;
         }
    cont=0;

    for(Integer i=0;i<RT0Ndofs;i++)
      // loop on dimension
      for(Integer j=0;j<PnNdofs;j++)
        {
        // loop on dimension
        if(j!=ManifoldDim-i)
        {for(Integer k=0;k<ManifoldDim;k++)
           func[cont](k,0)=func[cont](k,0)*p1[j](0,0); 
         // std::cout<<cont<<"   "<<std::endl;
         // std::cout<<func[cont]<<std::endl;         
         cont++;
         }
         }
    // // the i-th column represents the coefficients of the i-th momenta shape function
    // // auto inv_mat=inverse(mat);
    for(Integer i=0;i<ManifoldDim;i++)
        {
            for(Integer j=0;j<ManifoldDim;j++)
            {
                func[cont+i](j,0)=0;
                for(Integer k=0;k<ManifoldDim;k++)
                   func[cont+i](j,0)+=mat_coeffs(k,i)*rt0[ManifoldDim-k](j,0)*p1[k](0,0);
                
            }
         // std::cout<<cont+i<<"   "<<std::endl;
         // std::cout<<func[cont+i]<<std::endl;
        }

}
};







template<Integer Dim,Integer ManifoldDim>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, DivergenceOperator, RaviartThomasFE, 1>
{
public:
  // static constexpr Integer ManifoldDim=2;
  static constexpr Integer Ndofs=NumberOfRaviartThomasSimplexDofs<ManifoldDim,1>();
  static constexpr Integer RT0Ndofs=NumberOfRaviartThomasSimplexDofs<ManifoldDim,0>();
  static constexpr Integer P1Ndofs=NumberOfLagrangianSimplexDofs<ManifoldDim,1>();

  using OutputRT0=typename ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 0>::Output;
  using OutputP1=typename ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, LagrangeFE, 1>::Output;
  using OutputDivRT0=typename ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, DivergenceOperator, RaviartThomasFE, 0>::Output;
  using OutputGradP1=typename ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, GradientOperator, LagrangeFE, 1>::Output;
  // static constexpr auto barycenter=LinearSimplex<ManifoldDim>::barycenter();
  static constexpr auto mat_coeffs=SimplexRT1Coefficient<Dim,ManifoldDim>::value();



  using Output=Vector<Matrix<Real, 1, 1>,Ndofs>;
  // using Output=Vector<Matrix<Real, 1, 1>,8>;
  // using OutputRT0=typename ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, RaviartThomasFE, 0>::Output;
  // using OutputP1=typename ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 1>::Output;
  // using OutputDivRT0=typename ReferenceShapeFunctionValue<Simplex<Dim,2>, DivergenceOperator, RaviartThomasFE, 0>::Output;
  // using OutputGradP1=typename ReferenceShapeFunctionValue<Simplex<Dim,2>, GradientOperator, LagrangeFE, 1>::Output;


constexpr inline static void 
apply(const Vector<Real,ManifoldDim>& point, Output & func)
{

    // OutputRT0 rt0;
    // OutputP1 p1;
    // OutputDivRT0 div_rt0;
    // OutputGradP1 grad_p1;

    // ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, RaviartThomasFE, 0>::apply(point,rt0);
    // ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 1>::apply(point,p1);

    // ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, RaviartThomasFE, 0>::apply(point,div_rt0);
    // ReferenceShapeFunctionValue<Simplex<Dim,2>, IdentityOperator, LagrangeFE, 1>::apply(point,grad_p1);


    OutputRT0 rt0;
    OutputP1 p1;
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 0>::apply(point,rt0);
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, LagrangeFE, 1>::apply(point,p1);


    OutputDivRT0 divRT0;
    OutputGradP1 gradP1;
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, DivergenceOperator, RaviartThomasFE, 0>::apply(point,divRT0);
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, GradientOperator, LagrangeFE, 1>::apply(point,gradP1);


    Integer cont=0;
    // loop on RT0 functions
    for(Integer i=0;i<RT0Ndofs;i++)
      // loop on dimension
      for(Integer j=0;j<P1Ndofs;j++)
        {
           // loop on dimension
            if(j!=ManifoldDim-i)
            {
               // for(Integer k=0;k<ManifoldDim;k++)
               func[cont](0,0)=divRT0[i](0,0)*p1[j](0,0); 
             cont++;
             }
         }
    cont=0;
    // loop on RT0 functions
    for(Integer i=0;i<RT0Ndofs;i++)
      // loop on dimension
      for(Integer j=0;j<P1Ndofs;j++)
        {
           // loop on dimension
            if(j!=ManifoldDim-i)
            {
               for(Integer k=0;k<ManifoldDim;k++)
               func[cont](0,0)+=rt0[i](k,0)*gradP1[j](k,0); 
             cont++;
             }
         }


    for(Integer i=0;i<ManifoldDim;i++)
            {
                func[cont+i](0,0)=0;
                for(Integer k=0;k<ManifoldDim;k++)
                   func[cont+i](0,0)+=mat_coeffs(k,i)*divRT0[ManifoldDim-k](0,0)*p1[k](0,0);

                for(Integer k=0;k<ManifoldDim;k++)
                   for(Integer j=0;j<ManifoldDim;j++)
                   func[cont+i](0,0)+=mat_coeffs(k,i)*rt0[ManifoldDim-k](j,0)*gradP1[k](j,0);
            }


   
    // const Real a=0.333333333333333;
    // const Real b=-0.666666666666666;
    // const auto& xi=point[0];
    // const auto& eta=point[1];

    // func[cont](0,0)=(a * (1. - xi - eta) - a* xi -a* (2*xi-1) + a*(1. - xi - eta) - a*eta - a * eta );
    // func[cont+1](0,0)= -  ( b * (1. - xi - eta) - b * xi - a* (2*xi-1) + b*(1. - xi - eta) - b * eta  - a * xi );
 


   //  Output func2{
   //      {3*(1-xi-eta)},
   //      {3*xi},
   //      {3*(1-xi-eta)},
   //      {3*eta},
   //      {3*xi},
   //      {3*eta},
   //      {  ( a * (1. - xi - eta) - a* xi -a* (2*xi-1) + a*(1. - xi - eta) - a*eta - a * eta )},
   //      {-  ( b * (1. - xi - eta) - b * xi - a* (2*xi-1) + b*(1. - xi - eta) - b * eta  - a * xi )}
   //      // {3*eta},
   //      // {3*xi}
   //  };

   //  func=func2;
}

 inline static void 
apply2(const Vector<Real,ManifoldDim> point)
{

    OutputRT0 rt0;
    OutputP1 p1;
    Output func;
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, RaviartThomasFE, 0>::apply(point,rt0);
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, IdentityOperator, LagrangeFE, 1>::apply(point,p1);


    OutputDivRT0 divRT0;
    OutputGradP1 gradP1;
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, DivergenceOperator, RaviartThomasFE, 0>::apply(point,divRT0);
    ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, GradientOperator, LagrangeFE, 1>::apply(point,gradP1);


    Integer cont=0;
    // loop on RT0 functions
    for(Integer i=0;i<RT0Ndofs;i++)
      // loop on dimension
      for(Integer j=0;j<P1Ndofs;j++)
        {
           // loop on dimension
            if(j!=ManifoldDim-i)
            {
               // for(Integer k=0;k<ManifoldDim;k++)
               func[cont](0,0)=divRT0[i](0,0)*p1[j](0,0); 
             // std::cout<<cont<<"   "<<func[cont](0,0)<<std::endl;
             cont++;

             }
         }
    cont=0;
    // loop on RT0 functions
    for(Integer i=0;i<RT0Ndofs;i++)
      // loop on dimension
      for(Integer j=0;j<P1Ndofs;j++)
        {
           // loop on dimension
            if(j!=ManifoldDim-i)
            {
               for(Integer k=0;k<ManifoldDim;k++)
               func[cont](0,0)+=rt0[i](k,0)*gradP1[j](k,0); 
             // std::cout<<cont<<"   "<<func[cont](0,0)<<std::endl;
             cont++;
             }
         }


    for(Integer i=0;i<ManifoldDim;i++)
            {
                func[cont+i](0,0)=0;
                for(Integer k=0;k<ManifoldDim;k++)
                   func[cont+i](0,0)+=mat_coeffs(k,i)*divRT0[ManifoldDim-k](0,0)*p1[k](0,0);

                for(Integer k=0;k<ManifoldDim;k++)
                   for(Integer j=0;j<ManifoldDim;j++)
                   func[cont+i](0,0)+=mat_coeffs(k,i)*rt0[ManifoldDim-k](j,0)*gradP1[k](j,0);
               // std::cout<<cont+i<<"   "<<func[cont+i](0,0)<<std::endl;
            }

}

};


template<Integer Dim,Integer ManifoldDim, Integer Order>
class ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim>, TraceOperator, RaviartThomasFE, Order>
{
 public: 
  using Output=typename ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim-1>, IdentityOperator, LagrangeFE, Order>::Output;
constexpr inline static void 
apply(const Vector<Real,ManifoldDim-1>& point, Output & func)
  {
      Real one_frac_vol_surface=1.0/reference_simplex_volume<ManifoldDim-1>();
      ReferenceShapeFunctionValue<Simplex<Dim,ManifoldDim-1>, IdentityOperator, LagrangeFE, Order>::apply(point,func,one_frac_vol_surface);
  }
};

// template<Integer Dim>
// constexpr void value<Simplex<Dim,2>, DivergenceOperator, RaviartThomasFE, 1>
//  (const Vector<Real,2>& point, Vector<Matrix<Real, 1, 1>,8> & func)
// {
//     const auto& xi=point[0];
//     const auto& eta=point[1];
//     Vector<Matrix<Real, 1, 1>,8> func2{
//         {3*(1-xi-eta)},
//         {3*xi},
//         {3*(1-xi-eta)},
//         {3*eta},
//         {3*xi},
//         {3*eta},
//         {3*eta},
//         {3*xi}
//     };
//     func=func2;
// }








template<typename Elem,Integer FEFamily,Integer Order, typename ConstInput, typename ShapeFunctionCoefficientsCollection>
void shape_function_coefficients_init(const ConstInput& mesh_ptr,ShapeFunctionCoefficientsCollection& coeff);











template<Integer Dim,Integer ManifoldDim>
class SingleShapeFunctionCoefficientsCollection<Simplex<Dim,ManifoldDim>, RaviartThomasFE, 0>
{
 public: 
constexpr inline static  void 
apply(const Array<Real, ManifoldDim+1 >& outward,Array<Real, ManifoldDim+1 >& coeff)
{
    for(std::size_t i=0;i<ManifoldDim+1;i++)
    coeff[i]=outward[i];
    // coeff[1]=outward[1];
    // coeff[2]=outward[2];
}
};


template<Integer Dim>
class SingleShapeFunctionCoefficientsCollection<Simplex<Dim,2>, RaviartThomasFE, 1>
{
 public: 
constexpr inline static  void  
apply (const Array<Real, 3 >& outward,Array<Real, 8 >& coeff)
{
    coeff[0]=outward[0];
    coeff[1]=outward[0];
    coeff[2]=outward[1];
    coeff[3]=outward[1];
    coeff[4]=outward[2];
    coeff[5]=outward[2];
    coeff[6]=outward[0];
    coeff[7]=outward[1];
}
};

template<Integer Dim,Integer ManifoldDim>
class SingleShapeFunctionCoefficientsCollection<Simplex<Dim,ManifoldDim>, RaviartThomasFE, 1>
{
 public: 
 static constexpr Integer Ndofs=NumberOfRaviartThomasSimplexDofs<ManifoldDim,1>();
 static constexpr Integer RT0Ndofs=NumberOfRaviartThomasSimplexDofs<ManifoldDim,0>();
 static constexpr Integer P1Ndofs=NumberOfLagrangianSimplexDofs<ManifoldDim,1>();

constexpr inline static  void  
apply (const Array<Real, ManifoldDim +1 >& outward,Array<Real, Ndofs >& coeff)
{

    Integer cont=0;
    for(Integer i=0;i<RT0Ndofs;i++)
    {
           for(Integer j=0;j<ManifoldDim;j++)
            {
             coeff[cont]=outward[i];
             cont++;
            }
    
    }

    for(Integer i=0;i<ManifoldDim;i++)
       {
        coeff[cont+i]=outward[i];
        // cont++;
        }

}
};


// template<>
// void shape_function_coefficients_init<Simplex<2,2>, RaviartThomasFE, 0>
//  (const Vector<Real, 3 >& outward,Vector<Real, 3 >& coeff)
// {
//     coeff[0]=outward[0];
//     coeff[1]=outward[1];
//     coeff[2]=outward[2];
// }





// template<>
// void shape_function_coefficients_init<Simplex<2,2>, RaviartThomasFE, 1>
//  (const Vector<Real, 3 >& outward,Vector<Real, 8 >& coeff)
// {
//     coeff[0]=outward[0];
//     coeff[1]=outward[0];
//     coeff[2]=outward[1];
//     coeff[3]=outward[1];
//     coeff[4]=outward[2];
//     coeff[5]=outward[2];
//     coeff[6]=outward[0];
//     coeff[7]=outward[1];
// }

template<typename Elem,typename Operator, Integer FEFamily,Integer Order,typename single_type,Integer Ndofs,
Integer NQPoints,Integer Dim>
constexpr const Vector<Vector<single_type,NQPoints>,Ndofs> reference_shape_function_init(const Matrix<Real,NQPoints,Dim>&qp_points)
{
    Vector<Vector<single_type,NQPoints>,Ndofs> v;
    Vector<Real,Dim> qp_point;
    Vector<single_type,Ndofs> func;
    for(Integer qp=0;qp<NQPoints;qp++)
    {
        qp_point=qp_points.get_row(qp);
        // func=value<Elem,Operator,FEFamily,Order,single_type,Ndofs>(qp_point);
        ReferenceShapeFunctionValue<Elem,Operator,FEFamily,Order>::apply(qp_point,func);
        // value<Elem,Operator,FEFamily,Order>(qp_point,func);
        for(Integer n_dof = 0; n_dof < Ndofs; ++n_dof) {
            const_cast<single_type&>
            (static_cast<const std::array<single_type,NQPoints>& >
             ((static_cast<const std::array<Vector<single_type,NQPoints>,Ndofs>& >(v())[n_dof] )())[qp])=
            static_cast<const std::array<single_type,Ndofs>& >(func())[n_dof];
        }
    }
    return v;
};


template<typename Elem,typename Operator, Integer FEFamily,Integer Order,typename single_type,Integer Ndofs,
Integer NQPoints,Integer Dim>
constexpr const Vector<Vector<single_type,NQPoints>,Ndofs> reference_shape_function_init(const Vector<Matrix<Real,Dim,1>,NQPoints>&qp_points)
{
    Vector<Vector<single_type,NQPoints>,Ndofs> v;
    Vector<Real,Dim> qp_point;
    Vector<single_type,Ndofs> func;
    for(Integer qp=0;qp<NQPoints;qp++)
    {
        for(std::size_t i=0;i<Dim;i++)
          qp_point[i]=qp_points[qp](i,0);
        // func=value<Elem,Operator,FEFamily,Order,single_type,Ndofs>(qp_point);
        ReferenceShapeFunctionValue<Elem,Operator,FEFamily,Order>::apply(qp_point,func);
        // value<Elem,Operator,FEFamily,Order>(qp_point,func);
        for(Integer n_dof = 0; n_dof < Ndofs; ++n_dof) {
            const_cast<single_type&>
            (static_cast<const std::array<single_type,NQPoints>& >
             ((static_cast<const std::array<Vector<single_type,NQPoints>,Ndofs>& >(v())[n_dof] )())[qp])=
            static_cast<const std::array<single_type,Ndofs>& >(func())[n_dof];
        }
    }
    return v;
};


template<typename T,Integer NQPoints,Integer Ndofs>
constexpr FQPValues<T,NQPoints,Ndofs>
weighted_reference_shape_function_init(const FQPValues<T,NQPoints,Ndofs>&input,const Array<Real,NQPoints>& weights)
{
    FQPValues<T,NQPoints,Ndofs> output;
    
    for(Integer qp=0;qp<NQPoints;qp++)
    {
        for(Integer n_dof = 0; n_dof < Ndofs; ++n_dof) {
            
            
            const_cast<T&>
            (static_cast<const Vector<T,NQPoints>& >
             ((static_cast<const FQPValues<T,NQPoints,Ndofs>& >(output)[n_dof] ))[qp])
            // output[n_dof][qp]
            =
            static_cast<const Array<Real,NQPoints>& >(weights)[qp]*
            const_cast<T&>
            (static_cast<const Vector<T,NQPoints>& >
             ((static_cast<const FQPValues<T,NQPoints,Ndofs>& >(input)[n_dof] ))[qp]);
        }
    }
    return output;
};




template<typename FunctionSpace,typename Operator>
class SingleTypeShapeFunction;



template<typename FunctionSpace>
class SingleTypeShapeFunction<FunctionSpace,IdentityOperator>
{
public:
    static constexpr Integer NComponents=FunctionSpace::NComponents;
    static constexpr Integer ShapeFunctionDim1=FunctionSpace::ShapeFunctionDim1;
    static constexpr Integer ShapeFunctionDim2=FunctionSpace::ShapeFunctionDim2;
    using SingleType=Matrix<Real, ShapeFunctionDim1, ShapeFunctionDim2>;
    // using TotType= Matrix<Real, NComponents, ShapeFunctionDim1>;
    using TotType=typename
    std::conditional_t<(ShapeFunctionDim1==1),
    Matrix<Real, NComponents*ShapeFunctionDim1, ShapeFunctionDim2>,
    std::conditional_t<(ShapeFunctionDim1>1 && ShapeFunctionDim2==1 && NComponents==1),
    Matrix<Real, ShapeFunctionDim1, 1>,//Matrix<Real, NComponents, ShapeFunctionDim1>,
    std::conditional_t<(ShapeFunctionDim1>1 && ShapeFunctionDim2==1 && NComponents>1),
    Matrix<Real,NComponents, ShapeFunctionDim1>,
    Matrix<Real,-6,-6,-1>
    >
    >
    >;
};
template<typename FunctionSpace>
class SingleTypeShapeFunction<FunctionSpace,GradientOperator>
{
public:
    static constexpr Integer NComponents=FunctionSpace::NComponents;
    static constexpr Integer ShapeFunctionDim1=FunctionSpace::ShapeFunctionDim1;
    static constexpr Integer ShapeFunctionDim2=FunctionSpace::ShapeFunctionDim2;
    static constexpr Integer Dim=FunctionSpace::Elem::Dim;
    using SingleType= typename
    std::conditional_t<(1==ShapeFunctionDim1 && 1==ShapeFunctionDim2),
    Matrix<Real, Dim , 1 >,//Vector<Real, Dim >,
    Matrix<Real, ShapeFunctionDim1 , ShapeFunctionDim2 * Dim>
    >;
    using TotType=typename
    std::conditional_t<(1==ShapeFunctionDim1 && 1==ShapeFunctionDim2 && NComponents==1),
    Matrix<Real, Dim , 1 >,//Vector<Real, Dim >,
    // Matrix<Real, ShapeFunctionDim1 , ShapeFunctionDim2 * Dim>,
    Matrix<Real, ShapeFunctionDim1 * NComponents, ShapeFunctionDim2 * Dim >
    >;
    
    // Matrix<Real, ShapeFunctionDim1 * NComponents, ShapeFunctionDim2 * Dim >;
};


template<typename FunctionSpace>
class SingleTypeShapeFunction<FunctionSpace,DivergenceOperator>
{
public:
    static constexpr Integer NComponents=FunctionSpace::NComponents;
    static constexpr Integer ShapeFunctionDim1=1;
    static constexpr Integer ShapeFunctionDim2=1;
    using SingleType=Matrix<Real,1,1>;
    using TotType= Matrix<Real, NComponents,1 >;
};


template<typename FunctionSpace>
class SingleTypeShapeFunction<FunctionSpace,TraceOperator>
{
public:
    static constexpr Integer NComponents=FunctionSpace::NComponents;
    static constexpr Integer ShapeFunctionDim1=1;
    static constexpr Integer ShapeFunctionDim2=1;
    using SingleType=Matrix<Real,1,1>;
    using TotType= Matrix<Real, NComponents,1 >;
};


template<typename FunctionSpace, typename Operator>
class SingleTypeShapeFunction<FunctionSpace,Transposed<Expression<Operator>>>
{
public:
    static constexpr Integer NComponents=FunctionSpace::NComponents;
    static constexpr Integer ShapeFunctionDim1=FunctionSpace::ShapeFunctionDim1;
    static constexpr Integer ShapeFunctionDim2=FunctionSpace::ShapeFunctionDim2;
    using SingleTypeShape=SingleTypeShapeFunction<FunctionSpace,Operator>;
    using SingleType=Transposed<typename SingleTypeShape::SingleType>;
    using TotType=Transposed<typename SingleTypeShape::TotType>;
    // using SingleTypeTmp=typename SingleTypeShape::SingleType;
    // using TotTypeTmp=typename SingleTypeShape::TotType;
    // using SingleType=Matrix<Real,SingleTypeTmp::Cols,SingleTypeTmp::Rows>;
    // using TotType=std::conditional_t<TotTypeTmp::Cols==TotTypeTmp::Rows,
    //                                  TotTypeTmp,
    //                                  Matrix<Real,TotTypeTmp::Cols*NComponents,TotTypeTmp::Rows/NComponents>>;
    
};



template< typename Elem_,typename BaseFunctionSpace, typename Operator_, typename QuadratureRule>
class ShapeFunction;

template< typename Elem_,typename BaseFunctionSpace, typename Operator_, typename QuadratureRule>
class ShapeFunction
: public Expression<ShapeFunction<Elem_,BaseFunctionSpace,Operator_,QuadratureRule>> 
{ 
public:
    using Elem=Elem_;
    using QuadratureElem=typename QuadratureRule::Elem;

    using Operator=Operator_;
    
    using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
    static constexpr Integer Dim=Elem::Dim;
    static constexpr Integer ManifoldDim=Elem::ManifoldDim;
    static constexpr Integer NComponents=BaseFunctionSpace::NComponents;
    static constexpr Integer NQPoints=QuadratureRule::NQPoints;
    // static constexpr Integer Ntot=FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,BaseFunctionSpace>>::value;
    // if QuadratureElem==Elem, then FunctionSpaceDofsPerSubEntityElem==FunctionSpaceDofsPerElem
    // if QuadratureElem::ManifoldDim=Elem::ManifoldDim-1, 
    // then FunctionSpaceDofsPerSubEntityElem counts the total number of dofs on a boundary face
    // static constexpr Integer Ntot=FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,BaseFunctionSpace>>::value;
    static constexpr Integer Ntot=FunctionSpaceDofsPerSubEntityElem<ElemFunctionSpace<Elem,BaseFunctionSpace>,QuadratureElem::ManifoldDim>::value;
    static constexpr Integer Ndofs=Ntot/NComponents;
    static constexpr Integer Order=BaseFunctionSpace::Order;
    static constexpr Integer FEFamily=BaseFunctionSpace::FEFamily;
    
    using SingleType   = typename SingleTypeShapeFunction<FunctionSpace,Operator>::SingleType;
    using VectorSingleType   = Vector<SingleType,Ndofs>;
    using tot_type= typename SingleTypeShapeFunction<FunctionSpace,Operator>::TotType;
    using qpvalues_type= QPValues<tot_type,NQPoints>;
    using type= FQPValues<tot_type,NQPoints,Ntot>;
    using Point = Vector<Real,Dim>;
    using QP = Matrix<Real,NQPoints,Dim>;
    using qp_points_type=typename QuadratureRule::qp_points_type;
    // it can be Elem<dim,manifolddim> for volumetric quadrature rule
    //           Elem<dim,manifolddim-1> for surface quadrature rule
    using Map=MapFromReference<Operator,QuadratureElem,BaseFunctionSpace::FEFamily>;
    
    static constexpr Integer ShapeFunctionDim1=SingleTypeShapeFunction<FunctionSpace,Operator>::ShapeFunctionDim1;
    static constexpr Integer ShapeFunctionDim2=SingleTypeShapeFunction<FunctionSpace,Operator>::ShapeFunctionDim2;
    
    static constexpr FQPValues<SingleType,NQPoints,Ndofs>
    reference_values{reference_shape_function_init<Elem,Operator,FEFamily,Order,SingleType,Ndofs>(QuadratureRule::qp_points)};
    

    // static constexpr FQPValues<SingleType,NQPoints,Ndofs>
    // weighted_reference_values{  weighted_reference_shape_function_init(reference_values,QuadratureRule::qp_sqrt_abs_weights)};
    
    constexpr const type& eval()const{return func_values_;}
    
    
    void init()
    {
        const auto& map=(*map_ptr);
        const auto& mapping=map();

        for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
        {
            for(Integer n_comp=0;n_comp<NComponents;n_comp++)
            {
                n_tot_=n_dof * NComponents +  n_comp ;
                n_=n_comp*ShapeFunctionDim1;
                for(Integer qp=0;qp<NQPoints;qp++)
                {
                    func_values_[n_tot_][qp].zero();
                    // func_tmp_=  mapping * weighted_reference_values[n_dof][qp];
                    func_tmp_=  mapping * reference_values[n_dof][qp];
                    
                    // se ncompontensts >1, allora assegni alla riga
                    assign<NComponents>(func_values_[n_tot_][qp],func_tmp_,n_,0);
                }
                
            }
        }
    }
    
    void init(const Array<Real,Ndofs> &alpha)
    {

        reference_shape_function_init<Elem,Operator,FEFamily,Order,SingleType,Ndofs>(QuadratureRule::qp_points);
        const auto& map=(*map_ptr);
        const auto& mapping=map();
        for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
        {
            for(Integer n_comp=0;n_comp<NComponents;n_comp++)
            {
                n_tot_=n_dof * NComponents +  n_comp ;
                n_=n_comp;
                for(Integer qp=0;qp<NQPoints;qp++)
                {
                    func_values_[n_tot_][qp].zero();
                    // func_tmp_=alpha[n_dof] * mapping * weighted_reference_values[n_dof][qp];
                    // std::cout<< "n_dof, qp, n_tot_, n_comp, n_=("<<n_dof<<", "<<qp<<", "<< n_tot_<<", "<<n_comp<<", "<< n_<<")"<<std::endl;
                    // std::cout<< "func_values_="<<func_values_[n_tot_][qp]<<std::endl;
                    func_tmp_=alpha[n_dof] * mapping * reference_values[n_dof][qp];
                    // std::cout<< "func_tmp_="<<func_tmp_<<std::endl;
                    assign<NComponents>(func_values_[n_tot_][qp],func_tmp_,n_,0);
                    // std::cout<< "func_values_ after="<<func_values_[n_tot_][qp]<<std::endl;
                }
                
            }
        }
        
    };



    // initialization of shape functions by choosing dynamically the reference quadrature points
    // the order of the quadrature rule (and so the number of points) is chosen at compile-time
    void init(const qp_points_type& reference_qp_points)
    {
        const auto& map=(*map_ptr);
        const auto& mapping=map();
        dynamic_reference_values_=reference_shape_function_init<Elem,Operator,FEFamily,Order,SingleType,Ndofs>(reference_qp_points);
        // std::cout<<"dynamic_reference_values_"<<std::endl;
        // std::cout<<dynamic_reference_values_<<std::endl;
        // std::cout<<"mapping"<<std::endl;
        // std::cout<<mapping<<std::endl;
        for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
        {
            for(Integer n_comp=0;n_comp<NComponents;n_comp++)
            {
                n_tot_=n_dof * NComponents +  n_comp ;
                n_=n_comp*ShapeFunctionDim1;
                for(Integer qp=0;qp<NQPoints;qp++)
                {
                    func_values_[n_tot_][qp].zero();
                    // func_tmp_=  mapping * weighted_reference_values[n_dof][qp];
                    func_tmp_=  mapping * dynamic_reference_values_[n_dof][qp];
                    
                    // se ncompontensts >1, allora assegni alla riga
                    assign<NComponents>(func_values_[n_tot_][qp],func_tmp_,n_,0);
                }
            }
        }
    }


    // initialization of shape functions by choosing dynamically the reference quadrature points
    // the order of the quadrature rule (and so the number of points) is chosen at compile-time
    void init(const Array<Real,Ndofs> &alpha, const qp_points_type& reference_qp_points)
    {
        const auto& map=(*map_ptr);
        const auto& mapping=map();
        dynamic_reference_values_=reference_shape_function_init<Elem,Operator,FEFamily,Order,SingleType,Ndofs>(reference_qp_points);
        // std::cout<<"dynamic_reference_values_"<<std::endl;
        // std::cout<<dynamic_reference_values_<<std::endl;
        // std::cout<<"mapping"<<std::endl;
        // std::cout<<mapping<<std::endl;
                
        for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
        {
            for(Integer n_comp=0;n_comp<NComponents;n_comp++)
            {
                n_tot_=n_dof * NComponents +  n_comp ;
                n_=n_comp*ShapeFunctionDim1;
                for(Integer qp=0;qp<NQPoints;qp++)
                {
                    func_values_[n_tot_][qp].zero();
                    // func_tmp_=  mapping * weighted_reference_values[n_dof][qp];
                    func_tmp_=alpha[n_dof] * mapping * dynamic_reference_values_[n_dof][qp];
                    
                    // se ncompontensts >1, allora assegni alla riga
                    assign<NComponents>(func_values_[n_tot_][qp],func_tmp_,n_,0);
                }
            }
        }
    }


    constexpr void init_map(const Map& map){map_ptr=std::make_shared<Map>(map);}
    
    ShapeFunction(const Map& map):
    map_ptr(std::make_shared<Map>(map))
    {}
    
    ShapeFunction(const ShapeFunction& shape):
    map_ptr(std::make_shared<Map>(shape.map()))
    {}
    
    
    ShapeFunction(){}
    
    const auto& map()const{return (*map_ptr);}
    
private:
    SingleType func_tmp_;
    VectorSingleType func_;
    Point qp_point_;
    FQPValues<SingleType,NQPoints,Ndofs> component_func_values_;
    FQPValues<SingleType,NQPoints,Ndofs> dynamic_reference_values_;
    type func_values_;
    std::shared_ptr<Map> map_ptr;
    Integer n_tot_;
    Integer n_;
};







template<typename Elem,typename BaseFunctionSpace,typename Operator, typename QuadratureRule>
constexpr FQPValues<typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QuadratureRule>::SingleType,
ShapeFunction<Elem,BaseFunctionSpace,Operator,QuadratureRule>::NQPoints,
ShapeFunction<Elem,BaseFunctionSpace,Operator,QuadratureRule>::Ndofs>
ShapeFunction<Elem,BaseFunctionSpace,Operator,QuadratureRule>::reference_values;

// template<typename Elem,typename BaseFunctionSpace,typename Operator, typename QuadratureRule>
// constexpr FQPValues<typename ShapeFunction<Elem,BaseFunctionSpace,Operator,QuadratureRule>::SingleType,
// ShapeFunction<Elem,BaseFunctionSpace,Operator,QuadratureRule>::NQPoints,
// ShapeFunction<Elem,BaseFunctionSpace,Operator,QuadratureRule>::Ndofs>
// ShapeFunction<Elem,BaseFunctionSpace,Operator,QuadratureRule>::weighted_reference_values;













template< typename Elem_,typename BaseFunctionSpace, typename QuadratureRule>
class ShapeFunction<Elem_,BaseFunctionSpace,TraceOperator,QuadratureRule>
: public Expression<ShapeFunction<Elem_,BaseFunctionSpace,TraceOperator,QuadratureRule>> 
{ 
public:
    using Elem=Elem_;
    using QuadratureElem=typename QuadratureRule::Elem;
    
    using Operator=TraceOperator;
    
    using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
    static constexpr Integer Dim=Elem::Dim;
    static constexpr Integer ManifoldDim=Elem::ManifoldDim;
    static constexpr Integer NComponents=BaseFunctionSpace::NComponents;
    static constexpr Integer NQPoints=QuadratureRule::NQPoints;
    // static constexpr Integer Ntot=FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,BaseFunctionSpace>>::value;
    // if QuadratureElem==Elem, then FunctionSpaceDofsPerSubEntityElem==FunctionSpaceDofsPerElem
    // if QuadratureElem::ManifoldDim=Elem::ManifoldDim-1, 
    // then FunctionSpaceDofsPerSubEntityElem counts the total number of dofs on a boundary face
    // static constexpr Integer Ntot=FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,BaseFunctionSpace>>::value;
    static constexpr Integer NtotVolume=FunctionSpaceDofsPerSubEntityElem<ElemFunctionSpace<Elem,BaseFunctionSpace>,Elem::ManifoldDim>::value;
    static constexpr Integer NdofsVolume=NtotVolume/NComponents;
    static constexpr Integer Ntot=FunctionSpaceDofsPerSubEntityElem<ElemFunctionSpace<Elem,BaseFunctionSpace>,QuadratureElem::ManifoldDim>::value;
    static constexpr Integer Ndofs=Ntot/NComponents;
    static constexpr Integer FEFamily=BaseFunctionSpace::FEFamily;
    static constexpr Integer Order=BaseFunctionSpace::Order;
    static constexpr Integer Continuity=BaseFunctionSpace::Continuity;
    static constexpr auto trace=TraceDofs<FunctionSpace>::dofs();//trace_dofs<FunctionSpace>();
    using SingleType   = typename SingleTypeShapeFunction<FunctionSpace,Operator>::SingleType;
    using VectorSingleType   = Vector<SingleType,Ndofs>;
    using tot_type= typename SingleTypeShapeFunction<FunctionSpace,Operator>::TotType;
    using qpvalues_type= QPValues<tot_type,NQPoints>;
    using type= FQPValues<tot_type,NQPoints,Ntot>;
    using Point = Vector<Real,Dim>;
    using QP = Matrix<Real,NQPoints,Dim>;
    using qp_points_type=typename QuadratureRule::qp_points_type;
    // it can be Elem<dim,manifolddim> for volumetric quadrature rule
    //           Elem<dim,manifolddim-1> for surface quadrature rule
    using Map=MapFromReference<Operator,QuadratureElem,BaseFunctionSpace::FEFamily>;
    
    static constexpr Integer ShapeFunctionDim1=SingleTypeShapeFunction<FunctionSpace,Operator>::ShapeFunctionDim1;
    static constexpr Integer ShapeFunctionDim2=SingleTypeShapeFunction<FunctionSpace,Operator>::ShapeFunctionDim2;
    
    static constexpr FQPValues<SingleType,NQPoints,Ndofs>
    reference_values{reference_shape_function_init<Elem,Operator,FEFamily,Order,SingleType,Ndofs>(QuadratureRule::qp_points)};
    
    using trace_type_tmp=typename decltype(trace)::value_type;
    using trace_type=ArrayChangeType<Real,trace_type_tmp>;
    // static constexpr FQPValues<SingleType,NQPoints,Ndofs>
    // weighted_reference_values{  weighted_reference_shape_function_init(reference_values,QuadratureRule::qp_sqrt_abs_weights)};
    
    constexpr const type& eval()const{return func_values_;}
    
    
    void init(const Integer face)
    {
        const auto& map=(*map_ptr);
        const auto& mapping=map();
        // decltype(func_values_) ok1(1);
        // SingleType ok2(2);
        // std::cout<<"SHAPE FUNCTION TRACE reference_values="<<std::endl;
        
        // std::cout<<QuadratureRule::qp_points<<std::endl;
        // std::cout<<reference_values<<std::endl;
        // std::cout<<mapping<<std::endl;
        
        // std::cout<<"func_values_[0][0]"<<std::endl;
        // std::cout<<func_values_[0][0]<<std::endl;
        // std::cout<<"func_values_"<<std::endl;
        
        // std::cout<<func_values_<<std::endl;
        // std::cout<<"func_tmp_"<<std::endl;
        // std::cout<<func_tmp_<<std::endl;
        for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
        {
            for(Integer n_comp=0;n_comp<NComponents;n_comp++)
            {
                
                n_tot_=n_dof * NComponents +  n_comp ;
                n_=n_comp*ShapeFunctionDim1;
                for(Integer qp=0;qp<NQPoints;qp++)
                {
                    func_values_[n_tot_][qp].zero();
                    // func_tmp_=  mapping * weighted_reference_values[n_dof][qp];
                    func_tmp_=  mapping * reference_values[n_dof][qp];
                    
                    // se ncompontensts >1, allora assegni alla riga
                    assign<NComponents>(func_values_[n_tot_][qp],func_tmp_,n_,0);
                   // std::cout<<"func_tmp_="<<mapping<<" "<<reference_values[n_dof][qp]<<std::endl;
                }
                
            }
        }
        
     // std::cout<<"ShapeFunction  trace end"<<std::endl;

    }
    
    void init(const Array<Real,NdofsVolume> &beta, const Integer face)
    {

        //face=0
        // const Integer face=0;
        // auto alpha=subarray(beta,trace[face]);
        // std::cout<<"face"<<face<<std::endl;
        //   std::cout<<"beta"<<beta<<std::endl;
        //   std::cout<<"trace[face]"<<trace[face]<<std::endl;
        subarray(alpha_,beta,trace[face]);
        
        // std::cout<<"init RT elements (coeffs)"<<std::endl;
        // std::cout<<beta<<std::endl;
        //  std::cout<<alpha_<<std::endl;
        
        const auto& map=(*map_ptr);
        const auto& mapping=map();
        for(Integer n_dof=0;n_dof<Ndofs;n_dof++)
        {
            for(Integer n_comp=0;n_comp<NComponents;n_comp++)
            {
                
                n_tot_=n_dof * NComponents +  n_comp ;
                n_=n_comp;
                for(Integer qp=0;qp<NQPoints;qp++)
                {
                    func_values_[n_tot_][qp].zero();
                    // func_tmp_=alpha[n_dof] * mapping * weighted_reference_values[n_dof][qp];
                    // std::cout<< "n_dof, qp, n_tot_, n_comp, n_=("<<n_dof<<", "<<qp<<", "<< n_tot_<<", "<<n_comp<<", "<< n_<<")"<<std::endl;
                    // std::cout<< "func_values_="<<func_values_[n_tot_][qp]<<std::endl;
                    func_tmp_=alpha_[n_dof] * mapping * reference_values[n_dof][qp];
                    // std::cout<< "func_tmp_="<<func_tmp_<<std::endl;
                    assign<NComponents>(func_values_[n_tot_][qp],func_tmp_,n_,0);
                    // std::cout<< "func_values_ after="<<func_values_[n_tot_][qp]<<std::endl;
                }
                
            }
        }
        // std::cout<<"init end"<<std::endl;
        
    };
    
    constexpr void init_map(const Map& map){map_ptr=std::make_shared<Map>(map);}
    
    ShapeFunction(const Map& map):
    map_ptr(std::make_shared<Map>(map))
    {}
    
    ShapeFunction(const ShapeFunction& shape):
    map_ptr(std::make_shared<Map>(shape.map()))
    {}
    
    
    ShapeFunction(){}
    
    const auto& map()const{return (*map_ptr);}
    
private:
    SingleType func_tmp_;
    VectorSingleType func_;
    Point qp_point_;
    FQPValues<SingleType,NQPoints,Ndofs> component_func_values_;
    type func_values_;
    std::shared_ptr<Map> map_ptr;
    Integer n_tot_;
    Integer n_;
    trace_type alpha_;
};


















































}



#endif






