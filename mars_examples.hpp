#ifndef MARS_EXAMPLES_HPP
#define MARS_EXAMPLES_HPP

#include "mars_simplex.hpp"
#include "mars_connectivity.hpp"
#include "mars_functionspace_dofmap.hpp"
#include "mars_functionspace.hpp"
#include "mars_shape_function.hpp"
#include "mars_matrix.hpp"
#include "mars_jacobian.hpp"




#include "mars_function.hpp"

namespace mars{


  using std::cout;
  using std::endl;







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
        // face 0, nodes
        case 0:
        switch(SubEntityDimNumber)
        {
          case 0: return 1;
          case 1: return 2;
          default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
        }
        // face 0, edges
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
   return n_dofs;
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
 const auto entity=Space::entity;
 Vector<Integer, n> dofs;
 Integer cont=0;
 Integer dofs_per_entity_cont=0;
 for(Integer ii=0; ii< Space::entity.size(); ii++)
 {
  const auto entityii=Space::entity[ii];
  if(entity[ii]<=ManifoldDim)
  {
    const auto dofs_per_entity=Space::dofs_per_entity[ii];
    const auto binomial_coeff=binomial_coefficient(ManifoldDim,entity[ii]+1);
    // loop on the entityes of entityii dim
    for(Integer jj=0; jj<binomial_coeff; jj++)
      {
        // loop on the dofs of the given entity
        dofs[cont] = dofs_per_entity*simplex_face_sub_entities(ManifoldDim,face,entityii,jj)+function_space_dofs_per_elem<Space>(ManifoldDim,ii);
        cont++;
        for(Integer kk=1; kk<dofs_per_entity; kk++)
          {
            dofs[cont] = dofs[cont-1]+1;
            cont++;
          }
      }
  }
 }
 return dofs;
}


template <typename Space>
 constexpr auto trace_dofs()
 {
  const auto ManifoldDim=Space::ManifoldDim;
  constexpr auto n_dofs_per_face=trace_surf_n_dofs<Space>();
  const auto n_faces=binomial_coefficient(ManifoldDim+1,ManifoldDim);
  Vector< Vector<Integer, n_dofs_per_face>, n_faces> vec;
  for(Integer ii=0;ii<n_faces;ii++)
    vec[ii]=trace_dofs<Space>(ii);
  return vec;
}


      template<typename Space,typename Operator,typename single_type,
               Integer NQPoints,Integer Dim>
       constexpr auto reference_trace_shape_function_init(const Integer face, const Matrix<Real,NQPoints,Dim>&qp_points)
       {
        
        Vector<Real,Dim> qp_point;
        const auto dofs=trace_dofs<Space>(face);
        const auto n_dofs=dofs.size();
        Vector<Vector<single_type,NQPoints>,n_dofs> v;
        Vector<single_type,n_dofs> func;
        
            for(Integer qp=0;qp<NQPoints;qp++)
            {
             qp_point=qp_points.get_row(qp);
             // func=value<Elem,Operator,FEFamily,Order,single_type,Ndofs>(qp_point);
             value<Space::Elem,Operator,Space::FEFamily,Space::Order>(qp_point,func);
              for(Integer n_dof = 0; n_dof < n_dofs; ++n_dof) {
                  const_cast<single_type&>
                  (static_cast<const std::array<single_type,NQPoints>& >
                   ((static_cast<const std::array<Vector<single_type,NQPoints>,n_dofs>& >(v())[n_dof] )())[qp])=
                  static_cast<const std::array<single_type,n_dofs>& >(func())[n_dof];
              }
            }
       return v;
      };

template<Integer N,Integer K>
 constexpr void combinations_generate_aux(
            Vector<Integer, K> &data,
            const Integer index, 
            const Integer i,
            Vector<Vector<Integer, K>, binomial_coefficient(N,K)> &combs,
            Integer &comb_index)
        {
            if(index == K) {
                for(Integer ii=0;ii<data.size();ii++)
                    combs[comb_index][ii]=data[ii];
                comb_index++;
                return;
            }

            if(i >= N) {
                return;
            }

            data[index] = i;

            combinations_generate_aux<N,K>(data, index+1, i+1, combs, comb_index);
            
            // current is excluded, replace it with next (Note that
            // i+1 is passed, but index is not changed)
            combinations_generate_aux<N,K>(data, index, i+1, combs, comb_index);
        }

        template<Integer N,Integer K >
        constexpr Vector<Vector<Integer, K>, binomial_coefficient(N,K)> combinations_generate()
        {
            Vector<Vector<Integer, K>, binomial_coefficient(N,K)> combs;
            Vector<Integer, K> data;
            Integer comb_index = 0;
            combinations_generate_aux<N,K>(data, 0, 0, combs, comb_index);
            return combs;
        }


   template<Integer Dim, Integer ManifoldDim>
inline constexpr auto jacobian_faces()
{
    static_assert(Dim >= ManifoldDim, "Dim must be greater or equal ManifoldDim");
    // static_assert(Npoints == ManifoldDim+1, "Npoints must be equal to ManifoldDim+1");
    Vector<Vector<Real, Dim>,ManifoldDim> points;
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


void boundary_example()
{
  // triangles, face 0
  static_assert(simplex_face_sub_entities(2,0,0,0)==0,"");
  static_assert(simplex_face_sub_entities(2,0,0,1)==1,"");
  static_assert(simplex_face_sub_entities(2,0,1,0)==0,"");
  // triangles, face 1
  static_assert(simplex_face_sub_entities(2,1,0,0)==0,"");
  static_assert(simplex_face_sub_entities(2,1,0,1)==2,"");
  static_assert(simplex_face_sub_entities(2,1,1,0)==1,"");
  // triangles, face 2
  static_assert(simplex_face_sub_entities(2,2,0,0)==1,"");
  static_assert(simplex_face_sub_entities(2,2,0,1)==2,"");
  static_assert(simplex_face_sub_entities(2,2,1,0)==2,"");
  // tetrahedron, face 0
  static_assert(simplex_face_sub_entities(3,0,0,0)==0,"");
  static_assert(simplex_face_sub_entities(3,0,0,1)==1,"");
  static_assert(simplex_face_sub_entities(3,0,0,2)==2,"");
  static_assert(simplex_face_sub_entities(3,0,1,0)==0,"");
  static_assert(simplex_face_sub_entities(3,0,1,1)==1,"");
  static_assert(simplex_face_sub_entities(3,0,1,2)==3,"");
  static_assert(simplex_face_sub_entities(3,0,2,0)==0,"");
  // tetrahedron, face 1
  static_assert(simplex_face_sub_entities(3,1,0,0)==0,"");
  static_assert(simplex_face_sub_entities(3,1,0,1)==1,"");
  static_assert(simplex_face_sub_entities(3,1,0,2)==3,"");
  static_assert(simplex_face_sub_entities(3,1,1,0)==0,"");
  static_assert(simplex_face_sub_entities(3,1,1,1)==2,"");
  static_assert(simplex_face_sub_entities(3,1,1,2)==4,"");
  static_assert(simplex_face_sub_entities(3,1,2,0)==1,"");
  // tetrahedron, face 2
  static_assert(simplex_face_sub_entities(3,2,0,0)==0,"");
  static_assert(simplex_face_sub_entities(3,2,0,1)==2,"");
  static_assert(simplex_face_sub_entities(3,2,0,2)==3,"");
  static_assert(simplex_face_sub_entities(3,2,1,0)==1,"");
  static_assert(simplex_face_sub_entities(3,2,1,1)==2,"");
  static_assert(simplex_face_sub_entities(3,2,1,2)==5,"");
  static_assert(simplex_face_sub_entities(3,2,2,0)==2,"");
  // tetrahedron, face 3
  static_assert(simplex_face_sub_entities(3,3,0,0)==1,"");
  static_assert(simplex_face_sub_entities(3,3,0,1)==2,"");
  static_assert(simplex_face_sub_entities(3,3,0,2)==3,"");
  static_assert(simplex_face_sub_entities(3,3,1,0)==3,"");
  static_assert(simplex_face_sub_entities(3,3,1,1)==4,"");
  static_assert(simplex_face_sub_entities(3,3,1,2)==5,"");
  static_assert(simplex_face_sub_entities(3,3,2,0)==3,"");

  // pentatope, face 0
  static_assert(simplex_face_sub_entities(4,0,0,0)==0,"");
  static_assert(simplex_face_sub_entities(4,0,0,1)==1,"");
  static_assert(simplex_face_sub_entities(4,0,0,2)==2,"");
  static_assert(simplex_face_sub_entities(4,0,0,3)==3,"");
  static_assert(simplex_face_sub_entities(4,0,1,0)==0,"");
  static_assert(simplex_face_sub_entities(4,0,1,1)==1,"");
  static_assert(simplex_face_sub_entities(4,0,1,2)==2,"");
  static_assert(simplex_face_sub_entities(4,0,1,3)==4,"");
  static_assert(simplex_face_sub_entities(4,0,1,4)==5,"");
  static_assert(simplex_face_sub_entities(4,0,1,5)==7,"");
  static_assert(simplex_face_sub_entities(4,0,2,0)==0,"");
  static_assert(simplex_face_sub_entities(4,0,2,1)==1,"");
  static_assert(simplex_face_sub_entities(4,0,2,2)==3,"");
  static_assert(simplex_face_sub_entities(4,0,2,3)==6,"");
  static_assert(simplex_face_sub_entities(4,0,3,0)==0,"");

  // pentatope, face 1
  static_assert(simplex_face_sub_entities(4,1,0,0)==0,"");
  static_assert(simplex_face_sub_entities(4,1,0,1)==1,"");
  static_assert(simplex_face_sub_entities(4,1,0,2)==2,"");
  static_assert(simplex_face_sub_entities(4,1,0,3)==4,"");
  static_assert(simplex_face_sub_entities(4,1,1,0)==0,"");
  static_assert(simplex_face_sub_entities(4,1,1,1)==1,"");
  static_assert(simplex_face_sub_entities(4,1,1,2)==3,"");
  static_assert(simplex_face_sub_entities(4,1,1,3)==4,"");
  static_assert(simplex_face_sub_entities(4,1,1,4)==6,"");
  static_assert(simplex_face_sub_entities(4,1,1,5)==8,"");
  static_assert(simplex_face_sub_entities(4,1,2,0)==0,"");
  static_assert(simplex_face_sub_entities(4,1,2,1)==2,"");
  static_assert(simplex_face_sub_entities(4,1,2,2)==4,"");
  static_assert(simplex_face_sub_entities(4,1,2,3)==7,"");
  static_assert(simplex_face_sub_entities(4,1,3,0)==1,"");

  // pentatope, face 2
  static_assert(simplex_face_sub_entities(4,2,0,0)==0,"");
  static_assert(simplex_face_sub_entities(4,2,0,1)==1,"");
  static_assert(simplex_face_sub_entities(4,2,0,2)==3,"");
  static_assert(simplex_face_sub_entities(4,2,0,3)==4,"");
  static_assert(simplex_face_sub_entities(4,2,1,0)==0,"");
  static_assert(simplex_face_sub_entities(4,2,1,1)==2,"");
  static_assert(simplex_face_sub_entities(4,2,1,2)==3,"");
  static_assert(simplex_face_sub_entities(4,2,1,3)==5,"");
  static_assert(simplex_face_sub_entities(4,2,1,4)==6,"");
  static_assert(simplex_face_sub_entities(4,2,1,5)==9,"");
  static_assert(simplex_face_sub_entities(4,2,2,0)==1,"");
  static_assert(simplex_face_sub_entities(4,2,2,1)==2,"");
  static_assert(simplex_face_sub_entities(4,2,2,2)==5,"");
  static_assert(simplex_face_sub_entities(4,2,2,3)==8,"");
  static_assert(simplex_face_sub_entities(4,2,3,0)==2,"");

  // pentatope, face 3
  static_assert(simplex_face_sub_entities(4,3,0,0)==0,"");
  static_assert(simplex_face_sub_entities(4,3,0,1)==2,"");
  static_assert(simplex_face_sub_entities(4,3,0,2)==3,"");
  static_assert(simplex_face_sub_entities(4,3,0,3)==4,"");
  static_assert(simplex_face_sub_entities(4,3,1,0)==1,"");
  static_assert(simplex_face_sub_entities(4,3,1,1)==2,"");
  static_assert(simplex_face_sub_entities(4,3,1,2)==3,"");
  static_assert(simplex_face_sub_entities(4,3,1,3)==7,"");
  static_assert(simplex_face_sub_entities(4,3,1,4)==8,"");
  static_assert(simplex_face_sub_entities(4,3,1,5)==9,"");
  static_assert(simplex_face_sub_entities(4,3,2,0)==3,"");
  static_assert(simplex_face_sub_entities(4,3,2,1)==4,"");
  static_assert(simplex_face_sub_entities(4,3,2,2)==5,"");
  static_assert(simplex_face_sub_entities(4,3,2,3)==9,"");
  static_assert(simplex_face_sub_entities(4,3,3,0)==3,"");

  // pentatope, face 4
  static_assert(simplex_face_sub_entities(4,4,0,0)==1,"");
  static_assert(simplex_face_sub_entities(4,4,0,1)==2,"");
  static_assert(simplex_face_sub_entities(4,4,0,2)==3,"");
  static_assert(simplex_face_sub_entities(4,4,0,3)==4,"");
  static_assert(simplex_face_sub_entities(4,4,1,0)==4,"");
  static_assert(simplex_face_sub_entities(4,4,1,1)==5,"");
  static_assert(simplex_face_sub_entities(4,4,1,2)==6,"");
  static_assert(simplex_face_sub_entities(4,4,1,3)==7,"");
  static_assert(simplex_face_sub_entities(4,4,1,4)==8,"");
  static_assert(simplex_face_sub_entities(4,4,1,5)==9,"");
  static_assert(simplex_face_sub_entities(4,4,2,0)==6,"");
  static_assert(simplex_face_sub_entities(4,4,2,1)==7,"");
  static_assert(simplex_face_sub_entities(4,4,2,2)==8,"");
  static_assert(simplex_face_sub_entities(4,4,2,3)==9,"");
  static_assert(simplex_face_sub_entities(4,4,3,0)==4,"");

std::cout<<"Combinations<3>"<<std::endl;
Combinations<3,2>::print_all();
std::cout<<"Combinations<4>"<<std::endl;
Combinations<4,2>::print_all();
Combinations<4,3>::print_all();
std::cout<<"Combinations<5>"<<std::endl;
Combinations<5,2>::print_all();
Combinations<5,3>::print_all();
Combinations<5,4>::print_all();

  using Space=ElementFunctionSpace<Simplex<2,2>,LagrangeFE,3,1,1>;
  using single_type=Matrix<Real,Space::ShapeFunctionDim1,Space::ShapeFunctionDim2>;
  // constexpr const Integer m=trace_surf_n_dofs<Space>();
  constexpr const auto trace9=trace_dofs<Space>(0);
  constexpr const auto trace8=trace_dofs<Space>(1);
  constexpr const auto trace7=trace_dofs<Space>(2);

  constexpr const auto trace10=trace_dofs<Space>();


  // constexpr auto automa=reference_trace_shape_function_init<Space,IdentityOperator,single_type>(0, GaussPoints<Space::Elem,3>::qp_points_type);



std::cout<<"-------------------------------------------"<<std::endl;
     
  for(int ii=0;ii<trace10.size();ii++)
     {for(int jj=0;jj<trace10[ii].size();jj++)
          std::cout<<trace10[ii][jj]<<std::endl;
          std::cout<<std::endl;
     } 


 Simplex<3,3> simpl3;
 Simplex<3,2> simpl32;

 constexpr auto reference=Simplex<3,3>::reference;
 static_assert(reference[0][0]==0,"non ok");
 Integer vs[3];
 for(Integer ii=0;ii<4;ii++)
   {
    Combinations<3 + 1,3>::generate(ii,vs);
    for(Integer jj=0;jj<3;jj++)
      {for(Integer kk=0;kk<3;kk++)
        std::cout<<reference[vs[jj]][kk]<<" ";
          std::cout<<std::endl;}
    std::cout<<std::endl;
 
  }

  constexpr Integer dim=3;
  constexpr Integer manifolddim=2;

  Simplex<dim,manifolddim> elem32;
  Matrix<Real, dim,manifolddim> J;

  constexpr Integer Npoints=Simplex<dim,manifolddim>::Npoints;
  constexpr Vector<Vector<Real,dim>,3> points { {0,0,1},{2,3,4},{1,4,8} };

  for(Integer ii=0;ii<Npoints;ii++)
   elem32.nodes[ii]=ii;

 // points[0][0]=0;
 // points[0][1]=0;
 // points[0][2]=1;
 // points[1][0]=2;
 // points[1][1]=3;
 // points[1][2]=6;
 // points[2][0]=1;
 // points[2][1]=4;
 // points[2][2]=8;
 // constexpr auto J32=jacobian<3,3,4>();
 constexpr auto Jkml=combinations_generate<4,3>();
   std::cout<<"oo oooo  oooooo oo ooo"<<std::endl;

 for(Integer ii=0;ii<Jkml.size();ii++)
   {
   for(Integer jj=0;jj<Jkml[ii].size();jj++)
    std::cout<<Jkml[ii][jj]<<" ";
  std::cout<<std::endl;
   }

 constexpr auto Jkm=jacobian_faces<3,3>();
 std::cout<<std::endl;
  for(Integer ii=0;ii<Jkm.size();ii++)
   {
    std::cout<<Jkm[ii]<<std::endl;
   }
 // J32=describe(std::cout);
 


  constexpr Integer ManifoldDim=2;
  constexpr Integer Dim=2;
  using MeshT=Mesh<Dim, ManifoldDim>;  
  MeshT mesh;
  using Elem = typename MeshT::Elem; 
  read_mesh("../data/beam-tri.MFEM", mesh);
  mark_boundary(mesh);
  Bisection<MeshT> b(mesh);
  b.uniform_refine(1);
    // b.clear();
  Elem elem;
  Simplex<Dim,ManifoldDim-1> face;
  mesh.update_dual_graph();
  mark_boundary(mesh);
  for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
  {
   auto &elem_id=elem_iter;
   const auto & elem=mesh.elem(elem_id);
   std::cout<<std::endl<<" elem =="<<elem_iter<<std::endl;
     // elem.side();
   for(Integer  nn=0;nn<3;nn++)
     {std::cout<<elem.side_tags[nn];
      if(elem.side_tags[nn]>0)
       {  elem.side(nn,face);
        std::cout<<"---boundary: ";
        for(Integer mm=0;mm<2;mm++)
          std::cout<<face.nodes[mm]<<" ";
      }
      std::cout<<" "<<std::endl;
    }
  } 

  // Integer entity_nodes_from[3]; 
  // for(Integer ii=0;ii<Combinations<4, 3>::value;ii++)
  // {
  //   Combinations<4, 3>::generate(ii,entity_nodes_from);
  //   std::cout<<"mesh.n_boundary_sides="<<mesh.n_boundary_sides()<<std::endl;
  //   for(auto& i: entity_nodes_from)
  //     std::cout<<i<<std::endl;
  // }


  // mesh.find_side_tags();
  // const auto& sides=mesh.side_tags();
  // std::cout<<"sides="<<std::endl;
  // for(auto& i: sides)
  //   std::cout<<i<<std::endl;

}









void assembly_example()
{


  using Point=Vector<Real,3>;
  Point point{0,1,2};

  constexpr Integer ManifoldDim=2;
  constexpr Integer Dim=2;
  using MeshT=Mesh<Dim, ManifoldDim>;
  using MeshT2=Mesh<Dim+1, ManifoldDim+1>;
  MeshT2 mesh2;
  read_mesh("../data/beam-tet.MFEM", mesh2);

 

  MeshT mesh;
  using Elem = typename MeshT::Elem; 
  // read_mesh("../data/beam-tri.MFEM", mesh);
  read_mesh("../data/triangle_square.MFEM", mesh);


  constexpr Integer QPOrder=4;
  constexpr Integer NQPoints=GaussPoints<Elem,QPOrder>::NQPoints; 
  using QP=typename GaussPoints<Elem,QPOrder>::qp_points_type;

  constexpr Integer Npoints=Elem::Npoints;
  std::vector<Vector<Real,Dim>> points(Npoints);



 using FSspace1= FunctionSpace< MeshT, Lagrange1<2>>;
 FSspace1 FEspace1(mesh);
 using FSspace2= FunctionSpace< MeshT, Lagrange2<2>>;
 FSspace2 FEspace2(mesh);
 using FSspace3= FunctionSpace< MeshT, Lagrange1<1>,Lagrange2<1>>;
 FSspace3 FEspace3(mesh);
 using FSspace4= FunctionSpace< MeshT, Lagrange1<1>,Lagrange2<1>,Lagrange1<1>>;
 FSspace4 FEspace4(mesh);
 using FSspace5= FunctionSpace< MeshT, Lagrange1<2>>;
 FSspace5 FEspace5(mesh);
 using FSspace6= FunctionSpace< MeshT, Lagrange3<2>>;
 FSspace6 FEspace6(mesh);
 using FSspace7= FunctionSpace< MeshT, Lagrange1<1>,Lagrange2<2>,Lagrange3<3>,Lagrange1<2>,Lagrange2<2> >;
 FSspace7 FEspace7(mesh);
 using FSspace8= FunctionSpace< MeshT, RT0<1>,RT1<2>>;
 FSspace8 FEspace8(mesh);
 using FSspace9= FunctionSpace< MeshT, RT0<2> >;
 FSspace9 FEspace9(mesh);
 using FSspace10= FunctionSpace< MeshT, RT1<2>>;
 FSspace10 FEspace10(mesh);


 // auto W5=MixedFunctionSpace(MixedFunctionSpace(FEspace7,FEspace8),FEspace1);

 using AuxFSspace1= FunctionSpace< MeshT, Lagrange1<1> >;
 using AuxFSspace2= FunctionSpace< MeshT, Lagrange2<1> >;
 using AuxFSspace3= FunctionSpace< MeshT, Lagrange1<2> >;
 AuxFSspace1 AuxFEspace1(mesh);
 AuxFSspace2 AuxFEspace2(mesh);
 AuxFSspace3 AuxFEspace3(mesh);


 auto Wtrial=MixedFunctionSpace(FEspace1,FEspace3);
 auto Waux=AuxFunctionSpacesBuild(AuxFEspace1);//,AuxFEspace2,AuxFEspace3);
 auto W=FullSpaceBuild(Wtrial,Waux);

 auto f1 = MakeFunction<0,Function1>(W);
 // auto f2 = MakeFunction<1>(W);
 // auto f3 = MakeFunction<2>(W);

 auto u1 = MakeTrial<0>(W);
 auto u2 = MakeTrial<1>(W);
 auto u3 = MakeTrial<2>(W);

 auto v1 = MakeTest<0>(W);
 auto v2 = MakeTest<1>(W);
 auto v3 = MakeTest<2>(W);


  Jacobian<Elem> J(mesh);
  
  constexpr auto C=Constant<Scalar1>(0.5);
  constexpr auto alpha=Constant<Scalar1>(2.0);
  constexpr auto beta=Constant<Scalar1>(3.0);
  constexpr auto gamma=Constant<Scalar1>(3.0);
  constexpr auto id2=Constant<Identity2>();

  constexpr auto  alpha4=alpha.qp_eval<4>();

  constexpr auto matrix1=Constant<Mat1>(1.,2.);
  constexpr auto matrix2=Constant<Mat2>(3.,4.);


auto NewOp1=NewOperator(IdentityOperator()/alpha);
auto NewOp2=NewOperator(IdentityOperator()*alpha*f1);
// auto Epsilon=NewOperator((-f1)*(+C)*((+GradientOperator())+(+Transpose(GradientOperator()))));
// auto Epsilon=NewOperator((-f1)*Trace(C)*(+GradientOperator()+(Transpose(-GradientOperator()))));
// auto Epsilon=NewOperator(+(-GradientOperator()));
auto Epsilon=NewOperator(Trace(f1)*Trace(C)*(-GradientOperator()-GradientOperator())/Trace(Transpose(f1)-Trace(Transpose(C))));
///Trace(Transpose(Transpose(Trace(f1)))));
// auto Epsilon=NewOperator((GradientOperator()+Transpose(GradientOperator())));//+Transpose(GradientOperator())));

  auto bilinearform= L2Inner(Transpose(Transpose(u3)),-Transpose(v3))+
                    L2Inner(-Transpose(-u3),(v3))+ //+ L2Inner(Grad(u3),Grad(v1))+L2Inner(u3,v3)+
                    L2Inner(Trace(f1)*(+Transpose(u2)),(Transpose(v2))) +//+ L2Inner(Grad(u2),Grad(v1))+L2Inner(u2,v3)+
                    L2Inner(Epsilon(u1),Grad(v1));//+ L2Inner(Grad(u1),Grad(v1))+L2Inner(u1,v3);//+L2Inner(Grad(u2),Grad(v2));//+L2Inner(f3*u3,v3);
  auto linearform=
                  //L2Inner(Grad(v1),+Transpose(id2));//Transpose(Transpose((matrix1)+Transpose(matrix2))));//Transpose(-f1*(matrix1+matrix1)*Transpose(alpha*matrix1-matrix1)));//+L2Inner(f2,v2)+L2Inner(f1,v1);//+ L2Inner(f1,v1);
                  // L2Inner((+Transpose(C))*(-v2),-Transpose(f1))+
                  // L2Inner((Transpose(C)+Transpose(C))*v2,C)+
                  // L2Inner((Transpose(C)+(C))*v2,f1)+
                  // L2Inner((C+C)*v2,C)+
                  // L2Inner((-C-Transpose(C))*(-v2),Transpose(-f1))+
                  // L2Inner((-Transpose(C)-Transpose(C))*v2,C)+
                  // L2Inner((-Transpose(C)-(C))*v2,f1/C)+
                  L2Inner(Trace(id2)*id2,(Epsilon(v1)));

  // // auto bilinearform= 
  //                   L2Inner((u3),(v3))+ //+ L2Inner(Grad(u3),Grad(v1))+L2Inner(u3,v3)+
  //                   L2Inner((u2),(v2)) +//+ L2Inner(Grad(u2),Grad(v1))+L2Inner(u2,v3)+
  //                   L2Inner(alpha*beta*(u1),(v1));//+ L2Inner(Grad(u1),Grad(v1))+L2Inner(u1,v3);//+L2Inner(Grad(u2),Grad(v2));//+L2Inner(f3*u3,v3);
  // auto linearform= L2Inner(f1,(v1));//+L2Inner(f2,v2)+L2Inner(f1,v1);//+ L2Inner(f1,v1);


  auto bilinear_form=general_form(bilinearform);
  auto linear_form=general_form(linearform);



  auto shape_coefficients=shape_function_coefficients(bilinear_form,linear_form);
  auto reference_maps=reference_maps2(bilinear_form,linear_form);
  auto shape_functions=shape_functions2(shape_coefficients,reference_maps,bilinear_form,linear_form);
  shape_coefficients.init(mesh);
  
  auto eval_bilinear_form=Eval(bilinear_form,shape_functions);
  auto eval_linear_form=Eval(linear_form,shape_functions);

  J.init(0);
  reference_maps.init(J);
  shape_coefficients.init(0);
  shape_functions.init(J);///////////////<------------------------ problema qui


  std::cout<<"------_______-----_______-----_______-----_______------"<<std::endl;
  std::cout<<"------_______-----BEFORE EVALUATION-----_______--------"<<std::endl;
  eval_bilinear_form.apply(J);
  eval_linear_form.apply(J);


 // auto newform= L2Inner(Epsilon(u1),Epsilon(v1));//+ L2Inner(Grad(u1),Grad(v1))+L2Inner(u1,v3);//+L2Inner(Grad(u2),Grad(v2));//+L2Inner(f3*u3,v3);
 // decltype(shape_functions)::TupleOfTupleCompositeShapeFunctionTensor ee(4);

// auto oks=Eval(decltype(linearform)(linearform.left(),linearform.right()),shape_functions);
// local_tensor_
// OperatorAndQuadratureTupleType<decltype(bilinearform)::Left>::type e2e(5);

// OperatorAndQuadratureTupleType<decltype(linearform)::Left>::type ee(5);

// Vector<Vector<Matrix<Real,1,2>,1>,1> vecme({{6,7}});
// FQPValues<Matrix<Real,1,2>,1,1> fqp(vecme);
// FQPValues<Transposed<Matrix<Real,1,2>>,1,1> fqtrans(fqp);

// std::cout<<fqp<<std::endl;
// std::cout<<fqtrans<<std::endl;

// OperatorType<decltype(Transpose(Epsilon(v1))), GaussPoints<Simplex<2, 2>, 1> > aa;
// decltype(shape_functions)::TupleCompositeOperatorsAndQuadrature ok03(0);
// OperatorAndQuadratureTupleType<decltype(linearform)>::type ee(5);
// OperatorAndQuadratureTupleType<decltype(linearform)>::composite_type e4e(5);

// decltype(L2Inner((u3),(v3+v2))) ee1(5);
// decltype(L2Inner((u3),(v3-v2))) ee2(5);
// decltype(L2Inner((u3+u2),(v3))) ee3(5);
// decltype(L2Inner((u3-u2),(v3))) ee4(5);
// decltype(L2Inner((u3+u2),(v3+v2))) ee5(5);
// decltype(L2Inner((u3-u2),(v3-v2))) ee6(5);
// decltype(L2Inner((u3+u2),(v3-v2))) ee7(5);
// decltype(L2Inner((u3-u2),(v3+v2))) ee8(5);

// decltype(L2Inner((u3-u2-u1),(v3+v2+v1))) ee9(5);
// decltype(L2Inner((u3-u2-u1),(v3+v2-v1))) ee10(5);
// decltype(L2Inner((u3-u2-u1),(v3-v2+v1))) ee11(5);
// decltype(L2Inner((u3-u2-u1),(v3-v2-v1))) ee12(5);


// decltype(L2Inner((u3-u2-u1),(v3+v2+v1))) ee13(5);
// decltype(L2Inner((u3-u2+u1),(v3+v2-v1))) ee14(5);
// decltype(L2Inner((u3+u2-u1),(v3-v2+v1))) ee15(5);
// decltype(L2Inner((u3+u2+u1),(v3-v2-v1))) ee16(5);

// decltype(L2Inner((C*u3+C*u2+C*u1),(C*v3-C*v2-C*v1))) ee17(5);
// decltype(Epsilon) ok(1);
// decltype(reference_maps)::TupleOperatorsAndQuadrature o2k(1);
// decltype(reference_maps)::TupleCompositeOperatorsAndQuadrature o24k(1);


// decltype(shape_functions)::Form oo(5);
// OperatorAndQuadratureTupleType<decltype(reference_maps)::Form::Left>::type o2o(5);
// decltype(shape_functions)::TupleOfTupleShapeFunction rr4444(5,6,7,8,9,0,8,7,6,4);
// decltype(shape_functions)::TupleCompositeOperatorsAndQuadrature rr4(5,6,7,8,9,0,8,7,6,4);
// decltype(shape_functions)::TupleOfTupleShapeFunction rr334(5,6,7,8,9,0,8,7,6,4);
// decltype(shape_functions)::TupleOfTupleCompositeShapeFunction e2s(6,4,44,4,56,7,7);
// decltype(shape_functions)::TupleOfTupleCompositeShapeFunctionTensor rr(5,6,7,8,9,0,8,7,6,4);
  // OperatorAndQuadratureTupleType<decltype(reference_maps)::Form::Right>::L2prod::type ee(5,4);
// OperatorAndQuadratureTupleType<OperatorAndQuadratureTupleType<decltype(reference_maps)::Form::Right>::L2prod::type>::type o24o(5,4);


// OperatorType<Transposed<Expression<Test<FullSpace<MixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>,
//       BaseFunctionSpace<-10, 1, 1, 2>>, FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>, BaseFunctionSpace<-10, 2, 1, 1> > >,
//       AuxMixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>> > >, 0,
      
//       Addition<Expression<Transposed<Expression<ConstantTensor<Scalar1, double> > > >, 
//                Expression<Function<FullSpace<MixedSpace<FunctionSpace<Mesh<2, 2,DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 2>>, FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>, BaseFunctionSpace<-10, 2, 1, 1> > >,
//                                              AuxMixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>> > >, 
//                                    3, IdentityOperator, Function1> > > > > > ,
//       GaussPoints<Simplex<2, 2>, 3> > ok3(4);

// OperatorType<Transposed<Expression<Test<FullSpace<MixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>,
//       BaseFunctionSpace<-10, 1, 1, 2>>, FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>, BaseFunctionSpace<-10, 2, 1, 1> > >,
//       AuxMixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>> > >, 0,
//       Addition<Expression<Transposed<Expression<ConstantTensor<Scalar1, double> > > >, Expression<Function<FullSpace<MixedSpace<FunctionSpace<Mesh<2, 2,
//       DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 2>>, FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>, BaseFunctionSpace<-10, 2, 1, 1> > >,
//       AuxMixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>> > >, 3, IdentityOperator, Function1> > > > > > ,
//       GaussPoints<Simplex<2, 2>, 3> > ok3(4);


// OperatorType<Multiplication<Expression<Transposed<Expression<Test<FullSpace<MixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>,
//       BaseFunctionSpace<-10, 1, 1, 2>>, FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>, BaseFunctionSpace<-10, 2, 1, 1> > >,
//       AuxMixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>> > >, 0,
//       Addition<Expression<Transposed<Expression<ConstantTensor<Scalar1, double> > > >, Expression<Function<FullSpace<MixedSpace<FunctionSpace<Mesh<2, 2,
//       DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 2>>, FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>, BaseFunctionSpace<-10, 2, 1, 1> > >,
//       AuxMixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>> > >, 3, IdentityOperator, Function1> > > > > > >,
//       Expression<Test<FullSpace<MixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 2>>, FunctionSpace<Mesh<2, 2,
//       DefaultImplementation>, BaseFunctionSpace<-10, 1, 1, 1>, BaseFunctionSpace<-10, 2, 1, 1> > >, AuxMixedSpace<FunctionSpace<Mesh<2, 2, DefaultImplementation>,
//       BaseFunctionSpace<-10, 1, 1, 1>> > >, 0, IdentityOperator> > >, GaussPoints<Simplex<2, 2>, 3> > ok(1);
// //   // auto ctx = CreateContext(bilinear_form, linear_form);

// //   //  for(i--)
// //   // ctx.set_element(0);
// //   // Eval(ctx, mat, vec);
// //   // local2glonal(mat, vev, dofs, gmat, gvecv)
//   // decltype(W)::FromElementFunctionSpacesToUniqueNumbersTupleType ee(6);
//   // FromElementFunctionSpacesToFirstSpaceTupleType<decltype(W)::FromElementFunctionSpacesToUniqueNumbersTupleType> ese(6);
//   // decltype(shape_functions)::TupleOfTupleCompositeShapeFunction o444(5);
//   // decltype(bilinear_form)::FunctionSpace s5e(6);

//   // TupleOfCompositeTensor<decltype(bilinear_form)::FunctionSpace,decltype(shape_functions)::TupleOfTupleCompositeShapeFunction> edee(4,4,4,444,444,4);
//   // decltype(shape_functions)::TupleOfTupleCompositeShapeFunctionTensor ede4e(4);
  // decltype(bilinearform) eeee(6);
  // decltype(linearform) eee4e(6);
//   // OperatorAndQuadratureTupleType<decltype(NewOp2(v1))>::type ee(65,4,4,4,4);
// // decltype(reference_maps)::TupleOperatorsAndQuadrature eee(5);
// // decltype(reference_maps)::TupleCompositeOperatorsAndQuadrature ee4e(5);
// // decltype(reference_maps)::TupleOfTupleNoQuadrature ee444e(5);
// // decltype(reference_maps)::SpacesToUniqueFEFamilies ee3334e(5);
// // decltype(reference_maps)::Map ee3rr334e(5);
// // decltype(reference_maps)::UniqueMapping t34e(5);

// /////////////////////////////////////
//   ////////// PROBLEMA: SE FAI REMOVEDUPLICATES SUI TENSORI perdi informazioni rispetto che a farlo sull'eva;
//   ///////////// due eval diversi possono avere lo stesso tensore, che ha pero' un valore diverso!!!!!!!!!


//  // decltype(tup) tup0(4,55,5,5,5,5,5,5,55,555,5);
//  // decltype(tup1) tup11(4,55,5,5,5,5,5,5,55,555,5);
//  // decltype(tup2) tup22(4,55,5,5,5,5,5,5,55,555,5);
//   // decltype(tup3) tup22(4,55,5,5,5,5,5,5,55,555,5);
//   // decltype(tuple_get<0,0,0>(onlycomposite)) okoko(6);
//   // std::cout<<"qui->"<<tuple_get<0,0,0>(onlycomposite).composite_operator().left().eval()<<std::endl;
//    // decltype(shape_functions)::Form e(5);
//      // decltype(shape_functions)::TupleOfTupleShapeFunction e2s3(6,4,44,4,56,7,7);

//   // GetType<decltype(shape_functions)::TupleOfTupleShapeFunctionCombinations,1,0> e2s3(6,4,44,4,56,7,7);


//   // decltype(NewOp1) ok;
//   // decltype(reference_maps)::TupleCompositeOperatorsAndQuadrature e2s(6);
//   // decltype(reference_maps)::TupleOperatorsAndQuadrature es(6);
//   // decltype(shape_functions)::TupleOfTupleCompositeShapeFunction o3k(5);
//   // decltype(shape_functions)::TupleOfTupleCompositeShapeFunctionTensor kli(6);

// IsTestOrTrial<decltype(Transpose(alpha))>::type ok0(0);
// TupleRemoveNumber0<typename IsTestOrTrial<decltype(Transpose(+alpha))>::type> ok(1);
// TupleRemoveNumber0Helper<typename IsTestOrTrial<decltype(Transpose(+(+(+alpha))))>::type>::Tuple ok2(3);
 // const auto& tensors=shape_functions.composite_tensor();
 // std::cout<<"EVALUATION"<<std::endl;

// OperatorType<Addition<FQPValues<Transposed<Matrix<double, 1, 1, -1> >, 3L, 3L>,FQPValues<Matrix<double,1,1,-1>, 3, 3>>> ok(5);
// auto a1=Transpose(alpha); 
// std::cout<<linear_form().left().derived().qp_eval<4>()<<std::endl;

 // constexpr Matrix<Real,2,3> mat({1,2,3,4,5,6});
 //   Transposed<Matrix<Real,2,3>> transmat(mat);
 //  assert(transmat(0,0)==1 );
 //  assert(transmat(1,0)==2 );
 //  assert(transmat(2,0)==3 );
 //  assert(transmat(0,1)==4 );
 //  assert(transmat(1,1)==5 );
 //  assert(transmat(2,1)==6 );;

 //  static_assert(Transposed<Matrix<Real,2,3>>::Rows==3&&"okok");
 //  static_assert(Transposed<Matrix<Real,2,3>>::Cols==2&&"okok");

  // static_assert(transmat(2,1)==6 && "okok");

  // J.init(1);
  // reference_maps.init(J);
  // shape_coefficients.init(1);
  // shape_functions.init();

  // eval_bilinear_form.apply(J);
  // eval_linear_form.apply(J);


//  auto mm=Evaluation<Expression<decltype(alpha*u1)>, GaussPoints<Simplex<2,2>,3>>((alpha*u1));


 // decltype(shape_functions)::TupleOfTupleCompositeShapeFunctionTensor e(6);
 // decltype(shape_functions)::TupleOfTupleCompositeShapeFunction e4(6);

// FQPValues<Matrix<double, 1L, 2L>, 1L, 3L> hh;
//  auto ecc=form_of_composite_operator(NewOp1(v1));
//  auto ecc2=Evaluation<Expression<decltype(ecc)>,GaussPoints<Simplex<2,2>,1>>(ecc);
//  ecc2.apply(hh,J,shape_functions());
 
//  // OperatorType<decltype(NewOp1(v1)),GaussPoints<Simplex<2,2>,1>> e(6);

//   TupleOfCombinationFunctions<decltype(u1)::Nmax,MultipleAddition<decltype(bilinearform),decltype(linearform)>>::type onlycomposite=build_tuple_of_combination_functions<decltype(u1)::Nmax>(bilinearform,linearform);








// FormOfCompositeOperatorType<decltype(NewOp1(v1))>::type ok(1);
// decltype(shape_functions)::TupleOfTupleCompositeShapeFunction e4ee(6);
// decltype(shape_functions)::TupleOfTupleCompositeShapeFunctionTensor eee(6);



// TupleOfCombinationFunctions<decltype(u1)::Nmax,decltype(bilinearform),decltype(linearform)>::type eee(6);

 // decltype(onlycomposite) onle333ycomposite;
 // decltype(shape_functions()) e(6);
 // static_assert(IsSame<decltype(onlycomposite),typename TupleOfCombinationFunctions<decltype(u1)::Nmax,MultipleAddition<decltype(bilinearform),decltype(linearform)>>::type>::value && "they are not same");

// Number<decltype(bilinear_form)::FunctionSpace::Nuniquesubspaces > ee(5);
// Number<decltype(linear_form)::FunctionSpace::Nuniquesubspaces > e4(5);

// Number<decltype(u1)::Nmax> eee(6);



// decltype(build_tuple_of_combination_functions<decltype(u1)::Nmax>(decltype(bilinearform)(),decltype(linearform)())) ekl(6);

// decltype(shape_functions()) eee(6);

// decltype(u1)::UniqueElementFunctionSpacesTupleType oo(5);
// Number<decltype(u1)::value> o1(3);
// Number<decltype(u2)::value>o41(3);
// Number<decltype(u3)::value>o14(3);



// std::cout<<hh<<std::endl;
// OperatorType<Division<FQPValues<Matrix<double,1,2,-1>,1,3>,QPValues<Matrix<double,1,1,-1>,1>>> eee(5);
// OperatorType<Multiplication<Matrix<double, 1,
//       2, -1>, Matrix<double, 1, 1, -1> > >  eee(6,5,5,5);
// decltype(shape_functions()) oo(5);
// decltype(ecc2) l(5);

//   std::cout<<"apply2"<<std::endl;
  // auto Ku=4*Grad(u1);
  // auto Kv=4*Grad(v1);
  // auto esm=L2Inner(Ku,Kv);
  // auto generalesm=general_form(esm);
  // auto shape=shape_function_coefficients(generalesm);
  // auto referencemaps=reference_maps2(bilinear_form,generalesm);
  // auto shapefunctions=shape_functions2(bilinear_form,generalesm);


// decltype(referencemaps6)::UniqueMapping eeesa5(6);
// decltype(referencemaps7)::UniqueMapping eees6a(6);
// decltype(referencemaps8)::UniqueMapping eee7sa(6);
  // decltype(bilinear_form)::UniqueElementFunctionSpacesTupleType ok0(4);
  



    // using Left=decltype(f3);
    // using Right=decltype(v3);
    // using TestOrTrialLeft= IsTestOrTrial<Left>;
    // using TestOrTrialRight= IsTestOrTrial<Right>;
    // static constexpr Integer leftN=TestOrTrialLeft::number;
    // static constexpr Integer rightN=TestOrTrialRight::number;
    // static constexpr Integer TestOrTrialLeftValue =GetType<typename TestOrTrialLeft::type,0>::value;
    // static constexpr Integer TestOrTrialRightValue =GetType<typename TestOrTrialRight::type,0>::value;

    // using form= std::tuple<typename TypeOfForm<GetType<typename IsTestOrTrial<Left>::type,0>,
    //                                            GetType<typename IsTestOrTrial<Right>::type,0>
    //                         >::type >;
    // using TestTrialNumbers=typename FormTestTrialNumbers<GetType<form,0>::value,TestOrTrialLeftValue,TestOrTrialRightValue,leftN,rightN>::type;
    
    //  Number<GetType<form,0>::value> e3ee(6);
    // TestTrialNumbers eee(6);


  // decltype(L2Inner(f3,v3))::TestTrialNumbers eeee(6);
  // auto seee=L2Inner(f3,v3);


// Evaluation<Expression<decltype(f3*u3)>,GaussPoints<Simplex<2,2>,5>>::subtype okk(6); 


//////////////////// HERE



//  auto W4=MixedFunctionSpace(MixedFunctionSpace(FEspace1,FEspace2),FEspace5);
 
//  auto u7 = MakeTrial<0>(W4);
//  auto u8 = MakeTrial<1>(W4);
//  auto u9 = MakeTrial<2>(W4);


//  auto v7 = MakeTest<0>(W4);
//  auto v8 = MakeTest<1>(W4);
//  auto v9 = MakeTest<2>(W4);


//  Function1 func_init;
//  auto func1=Function<Function1,FSspace1>();
//  auto func2=Function<Function1,FSspace6>();
//  auto func3=Function<Function1,FSspace9>();
//  auto func4=Function<Function1,FSspace10>();


//   auto l8=  L2Inner(func1*u7,v7)+L2Inner(func2*u8,v8)+L2Inner(func3*u9,func4*v9);

//   auto generalform8=general_form(l8);
//   auto shape_coefficients8=shape_function_coefficients(generalform8);
//   // auto referencemaps8=reference_maps2(generalform8);
//   // auto shapefunctions8=shape_functions2(generalform8);

// auto dedee=generalform8();
// // decltype(dedee) oko(6);


// using unique=decltype(W5)::UniqueElementFunctionSpacesTupleType;


// using prova1=Gabriele1<decltype(l8),unique,TupleTypeSize<unique>::value-1>;
// using tipo1=typename prova1::type;

// auto ecc=Constructor<tipo1>::apply();


// // unique emskl(6);
// // UniqueElementFunctionSpacesTupleTypeUpdate<decltype(l8),unique>::UniqueElementFunctionSpacesTupleType osk(5);
// // tipo1 ok(2);
// // decltype(ecc) ok1          ke90(6);
// // auto esm=Constructor<decltype(L2Inner(Function<Function1,FSspace10,41>()*u7,func2*v7)) >::apply();
// // decltype(esm) esm2(4);

// // auto eeee=Function<Function1,FSspace10,41>();
// // unique ok00(6);
// // prova1::UniqueElementFunctionSpacesTupleType ok0(0);
// // decltype(ecc)  ok1(1);

// // std::cout<<ok.left().left().left()::eval(Vector<Real,2>{1,2})<<std::endl;
// // std::cout<<decltype(generalform8)::Form::Left::Left::Left::value<<std::endl;
// // std::cout<<decltype(generalform8)::Form::Left::Left::Right::value<<std::endl;
// // std::cout<<decltype(generalform8)::Form::Left::Right::Left::value<<std::endl;
// // std::cout<<decltype(generalform8)::Form::Left::Right::Right::value<<std::endl;
// // std::cout<<decltype(generalform8)::Form::Right::Left::Left::value<<std::endl;
// // std::cout<<decltype(generalform8)::Form::Right::Left::Right::value<<std::endl;
// // std::cout<<decltype(generalform8)::Form::Right::Right::Left::value<<std::endl;
// // std::cout<<decltype(generalform8)::Form::Right::Right::Right::value<<std::endl;


// // unique ok0(0);
// // gabtype::Left ok1(5);
// // gabtype::Right::Left ok2(5);
// // gabtype::Right::Right ok3(5);
// // Gab1::Tuple ok4(5);
// // Gab1::N ok3(5);
// //   std::cout<<"l8::Order="<<l8.Order<<std::endl;

//   // decltype(W5)::UniqueElementFunctionSpacesTupleType  ejkhd1(6);
//   // decltype(W5)::SpacesToUniqueFEFamily  ejkhd2(6);
//   // decltype(W5)::FromElementFunctionSpacesToUniqueNumbersTupleType ejkhd3(6);

//   // decltype(W5)::UniqueElementFunctionSpacesTupleType2  ejkhd4(6);
//   // decltype(W5)::SpacesToUniqueFEFamily2  ejkhd5(6);
//   // decltype(W5)::FromElementFunctionSpacesToUniqueNumbersTupleType2 ejkhd6(6);
//   auto l9= 
//            L2Inner(u8,v7)+L2Inner(u7,v9)+L2Inner(u7,v7)+
//            L2Inner(Grad(u7),Grad(v7))+L2Inner(u8,v8)+L2Inner(u9,v9)
//            ;


//   // decltype(generalform8)::TupleFunctionSpace a1(4);
//   // decltype(generalform8)::FunctionSpace a2(4);
//   // decltype(generalform8)::UniqueElementFunctionSpacesTupleType a3(4);


//   auto generalform=general_form(l9);
//   auto shape_coefficients_general=shape_function_coefficients(generalform);
//   auto referencemaps_general=reference_maps2(generalform);
//   auto shapefunctions_general=shape_functions2(generalform);
  
//   Jacobian<Elem> J(mesh);
//   shape_coefficients_general.init(mesh);


//   J.init(0);
//   shape_coefficients_general.init(0);
//   referencemaps_general.init(J);
//   shapefunctions_general.init(referencemaps_general,shape_coefficients_general);




//   std::cout<<"qui2"<<std::endl;
//   auto esc=get_spaces_ptr(L2Inner(Grad(u7),Grad(v7)));

//  std::cout<<"unweighted"<<std::endl;
//  std::cout<<ShapeFunctionDependent<Elem,Lagrange1<1>,IdentityOperator,GaussPoints<Elem,QPOrder>>::reference_values<<std::endl;
//  std::cout<<"weights"<<std::endl;
//  std::cout<<GaussPoints<Elem,QPOrder>::qp_weights<<std::endl;
//  std::cout<<"weighted"<<std::endl;
//  std::cout<<ShapeFunctionDependent<Elem,Lagrange1<1>,IdentityOperator,GaussPoints<Elem,QPOrder>>::weighted_reference_values<<std::endl;
//  // auto esm = Eval(L2Inner(mesh,u7,v7)+L2Inner(mesh,Grad(u7),Grad(v7)),shapefunctions_general);
//  // auto mmm = Eval(L2Inner(mesh,(u7),(v7)),shapefunctions_general);
//  Matrix<Real,3,3> mat_loc;
//  // mmm.apply(mat_loc,J);
//  std::cout<<"mat_loc="<<mat_loc<<std::endl;

//  auto eval_generalform=Eval(generalform,shapefunctions_general);

//  eval_generalform.apply(J);
//////////////////// HERE






















//  // using l2mass=decltype(L2Inner(mesh,u9,v9)+L2Inner(mesh,u8,v8)+L2Inner(mesh,u7,v7)+L2Inner(mesh,Grad(u7),Grad(v7)));
//  // // auto eshg=l2mass(mesh);
//  // // auto as1=L2Inner(mesh,l2mass::Left::Left(), l2mass::Left::Right());
//  // // auto as2=L2Inner(mesh,l2mass::Right::Left(),l2mass::Right::Right()) ;
//  // // auto eshg2=as1+as2;
//  // auto ecc=Constructor<l2mass>::apply(mesh);
//  // decltype(ecc) eekle(5);
//  // GetType< decltype(eval_generalform)::L2Products,0> ee1(Constructor< GetType< decltype(eval_generalform)::L2Products,0>>::apply(mesh));
//  // GetType< decltype(eval_generalform)::L2Products,1> ee2(Constructor< GetType< decltype(eval_generalform)::L2Products,1>>::apply(mesh));
//  // GetType< decltype(eval_generalform)::L2Products,2> ee3(Constructor< GetType< decltype(eval_generalform)::L2Products,2>>::apply(mesh));

//  // decltype(ee1) sslkre1(5);
//  // decltype(ee2) sslkre2(5);
//  // decltype(ee3) sslkre3(5);
//  // auto m1=Eval(ee1,shapefunctions_general);
//  // auto m2=Eval(ee2,shapefunctions_general);
//  // auto m3=Eval(ee3,shapefunctions_general);
//  // OKKAux<decltype(eval_generalform)::L2Products,0,decltype(shapefunctions_general)> eee(6);
//  using L2Products=decltype(eval_generalform)::L2Products;
//  using ShapeFunctions=decltype(shapefunctions_general);
//  // auto m4=EvalOfL2Inners<L2Products>(mesh,shapefunctions_general);


// // LocalMatrices<decltype(eval_generalform)>::type rrr(6);
// // decltype(eval_generalform)::TupleOfPairsNumbers ee5e(6);
// // std::cout<<abc<<std::endl;
// EvaluationOfL2Inners<decltype(eval_generalform),decltype(shapefunctions_general)> eee(shapefunctions_general);


// IdentityOperator token;
// eee.apply(token,J);
// EvalOfL2InnersAux<L2Products,0,ShapeFunctions> ee(5);


 // Matrix<double, 6, 6> mat0;
 // std::get<3>(m4).apply(mat0,J);

 // decltype(eval_generalform) eese(6);
 // decltype(m1) sslkre1(5);
 // decltype(m2) sslkre2(5);
 // decltype(m3) sslkre3(5);
  // auto emh=GetType< decltype(eval_generalform)::TupleOfEvals::type,0>(mesh);
// GetType<decltype(L2Inner(mesh,u7,v7))::FunctionSpace,2> rr3r(5);
// GetType<decltype(W4)::Spaces,1> rr5r(5);
// decltype(L2Inner(mesh,u7,v7))::TestTrialNumbers a1(1);
// decltype(L2Inner(mesh,u7,v7))::UniqueElementFunctionSpacesTupleType a2(1);
// decltype(L2Inner(mesh,u7,v7))::TupleFunctionSpace a3(1);
// decltype(L2Inner(mesh,u7,v7))::FunctionSpace a4(5);
// auto generaleval=Eval(generalform,shapefunctions_general);

// decltype(generalform)::FunctionSpace okeeeee44k(5);
// decltype(generaleval)::EvalType rrr(5);
// auto evalprova=Eval(L2Inner(mesh,u7,v7),shapefunctions_general);

// decltype(evalprova) sss(5);
// decltype(evalprova)::EvalLeft rrr4(5);
// decltype(evalprova)::EvalRight rrr5(5);
// evalprova.apply(mat_loc);
// generaleval.apply(mat_loc);
//  for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
//  {
//   elem=mesh.elem(elem_iter);
//   // mesh.points(elem_iter,points);

//   jacobian(elem,mesh.points(),J);

//   shape_coefficients.init(elem_iter);

//   referencemaps.init(J);

//   shapefunctions.init_map(referencemaps);
//   shapefunctions.init(shape_coefficients);
//   shapefunctions.init(referencemaps,shape_coefficients);


//   // std::cout<<"elem_iter=="<<elem_iter<<std::endl;
//   // for(Integer mm=0;mm<elem.nodes.size();mm++)
//   //   std::cout<<elem.nodes[mm]<<std::endl;
//   // std::cout<<J<<std::endl;
// }


  // std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  // std::cout<<shapefunctions.get<0,0>().reference_values<<std::endl;
  // std::cout<<shapefunctions.value<0,0>()<<std::endl;
  // std::cout<<"........."<<std::endl;
  // std::cout<<shapefunctions.get<0,1>().reference_values<<std::endl;
  // std::cout<<shapefunctions.value<0,1>()<<std::endl;
  // std::cout<<"........."<<std::endl;
  // std::cout<<shapefunctions.get<0,2>().reference_values<<std::endl;
  // std::cout<<shapefunctions.value<0,2>()<<std::endl;
  // std::cout<<"........."<<std::endl;
  // std::cout<<shapefunctions.get<0,3>().reference_values<<std::endl;
  // std::cout<<shapefunctions.value<0,3>()<<std::endl;
  // std::cout<<"........."<<std::endl;
  // std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  // std::cout<<shapefunctions.get<1,0>().reference_values<<std::endl;
  // std::cout<<shapefunctions.value<1,0>()<<std::endl;
  // std::cout<<"........."<<std::endl;
  // std::cout<<shapefunctions.get<1,1>().reference_values<<std::endl;
  // std::cout<<shapefunctions.value<1,1>()<<std::endl;
  // std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
  // std::cout<<shapefunctions.get<2,0>().reference_values<<std::endl;
  // std::cout<<shapefunctions.value<2,0>()<<std::endl;
  // std::cout<<"........."<<std::endl;
  // std::cout<<shapefunctions.get<2,1>().reference_values<<std::endl;
  // std::cout<<shapefunctions.value<2,1>()<<std::endl;

 // init_map< SpacesToUniqueFEFamiliesW1,MapTupleNumbersW1>(stuple,mapping);
 // decltype(shape_coefficients) lgy;
 // decltype(stuple) kstuple;
 // init(stuple,shape_coefficients);
// shapefunctions()=4;
 // auto evals=Eval(s);

 // evals.apply();
 // eval2.apply();


// std::cout<<"--------J"<<std::endl;
// std::cout<<J<<std::endl;
// std::cout<<"--------map 0,0:"<<std::endl;
// std::cout<<referencemaps.get<0,0>()<<std::endl;
// std::cout<<"--------map 0,1:"<<std::endl;
// std::cout<<referencemaps.get<0,1>()<<std::endl;
// std::cout<<"--------map 1,0:"<<std::endl;
// std::cout<<referencemaps.get<1,0>()<<std::endl;
// std::cout<<"--------map 1,1:"<<std::endl;
// std::cout<<referencemaps.get<1,1>()<<std::endl;

 // std::cout<<"tuple00.eval()[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[["<<std::endl;
 // auto tuple00=shapefunctions.get<0,0>();
 // std::cout<<tuple00.map()()<<std::endl;
 // tuple00.init(alpha);
 // std::cout<<tuple00.reference_values()<<std::endl;
 // std::cout<<tuple00.eval()<<std::endl;

 // SignedNormal<Elem> sn;
 // sn.init(mesh);
 // std::cout<<tuple00.eval()<<std::endl;
 // decltype(MakeTest<6>(W1))::BaseFunctionSpaces okll(4);
// TupleOfTupleShapeFunctionType<TupleSpacesW1,TupleOperatorsAndQuadratureW1> eew1;
// ShapeFunctionAssembly<TupleSpaces,TupleOperatorsAndQuadrature> ee2;

 // static_assert(ref[0][0](0,0)==2,"grad non ok");
 // static_assert(ref[1][0](0,0)==2,"grad non ok");
 // static_assert(ref[2][0](0,0)==2,"grad non ok");



 using tp5=std::tuple< Number<1>,Number<0>,Number<0>,Number<0>,Number<1>,Number<0>,Number<0>,Number<0>,Number<0>>;
 using tp6=std::tuple< Number<1>,Number<0>,Number<0>,Number<0>,Number<1>,Number<0>,Number<0>,Number<0>,Number<0>>;
 using nonzeroproduct=StaticBooleanMatrixMatrixMultiplicationFindNonZeros<3,3,tp5,3,3,tp6>;
 using booleanmatrix=TupleOfTupleToBooleanTuple<nonzeroproduct>;
 using findnonzeronumbers=FindNonZeroNumbers<booleanmatrix>;
 Matrix<Real,3,3> mat1{1,2,3,4,5,6,7,8,9};
 Matrix<Real,3,3> mat2{2,3,4,5,6,7,8,9,10};
Matrix<Real,3,3> mat3;//{0,0,0,0,0,0,0,0,0};
Matrix<Real,3,3> mat4;//{0,0,0,0,0,0,0,0,0};


constexpr Integer Rows=3;
constexpr Integer Cols=3;

using boolean1=std::tuple<One,Zero,Zero,One,One,Zero,Zero,Zero,One>;
StaticMatrix<Real,3,3,boolean1> ede;
StaticMatrix<Real,3,3> hfe;

Real nonzero;
Integer Constant=pow(10,7);
clock_t begin = clock();
for(Integer ii=0;ii<Constant;ii++)
{
    //     for(Integer i = 0; i < 3; ++i) {
    //       for(Integer j = 0; j < 3; ++j) {
    //           mat1(i, j) = mat1(i, j) /pow(1.0,1.0/ii);
    //         }
    //       }


    // StaticMatrixProduct<findnonzeronumbers,nonzeroproduct>(mat1,mat2,mat3);
    // StaticMatrixProduct<findnonzeronumbers,nonzeroproduct>(mat1,mat3,mat4);
    // mat3=mat1*mat2;


        // mat3=0;
        // for(Integer i = 0; i < 3; ++i) {
        //   for(Integer j = 0; j < 3; ++j) {
        //     for(Integer k = 0; k < 3; ++k) {
        //       mat3(i, k) += mat1(i, j) * mat2(j, k);
        //     }
        //   }
        // }

//         for(Integer i = 0; i < 3; ++i) {
//           for(Integer j = 0; j < 3; ++j) {
//             mat3(i, j) = mat1(i, 0) * mat2(0,j);
//             for(Integer k = 1; k < 1; ++k) {
//               mat3(i, j) += mat1(i, k) * mat2(k,j);
//             }
//           }
//         }


// nonzero+=mat3(0,0)+mat3(0,1)+mat3(0,2)+mat3(1,0)+mat3(1,1)+mat3(1,2)+mat3(2,0)+mat3(2,1)+mat3(2,2);

   // nonzero=StaticScalarProduct<Nonzero>(v1,v2);
    // StaticMatrixProduct<findnonzeronumbers,nonzeroproduct>(mat1,mat2,mat3,vec1,vec2);
    //     // nonzero=0;


              // mat3= mat1 * mat2;

     // nonzero=v1[3]*v2[2]+v1[4]*v2[5]+v1[5]*v2[8];
    // nonzero=0;
    // for(Integer ii=0;ii<3;ii++)
    //   nonzero+=mat1(1,ii)*mat2(ii,2);
}
clock_t end = clock();
double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

// std::cout<<"mat==="<<mat3<<std::endl;

std::cout<<" elapsed_secs------------>"<<elapsed_secs<<std::endl;

std::cout<<" nonzero------------>"<<nonzero<<std::endl;

// auto& dm =W1.dofmap();
// auto& dm0=std::get<0>(dm);
// auto& dm01=std::get<0>(dm0);
// auto& dm02=std::get<1>(dm0);
// auto& dm1=std::get<1>(dm);
// // auto& dm2=std::get<2>(dm);




// std::cout<<"dm01"<<std::endl;
// for(Integer ii=0;ii<dm01.size();ii++)
// {
//  for(Integer jj=0;jj<dm01[ii].size();jj++)
//   std::cout<<dm01[ii][jj]<<" ";
// std::cout<<std::endl;
// }
// std::cout<<"dm02"<<std::endl;

// for(Integer ii=0;ii<dm02.size();ii++)
// {
//  for(Integer jj=0;jj<dm02[ii].size();jj++)
//   std::cout<<dm02[ii][jj]<<" ";
// std::cout<<std::endl;
// }
// std::cout<<"dm1"<<std::endl;

// for(Integer ii=0;ii<dm1.size();ii++)
// {
//  for(Integer jj=0;jj<dm1[ii].size();jj++)
//   std::cout<<dm1[ii][jj]<<" ";
// std::cout<<std::endl;
// }
// std::cout<<"dm2"<<std::endl;
// for(Integer ii=0;ii<dm2.size();ii++)
// {
//    for(Integer jj=0;jj<dm2[ii].size();jj++)
//       std::cout<<dm2[ii][jj]<<" ";
//     std::cout<<std::endl;
// }





};















void functionspaces_example2D()
{

  constexpr Integer ManifoldDim=2;
  constexpr Integer Dim=2;
  using MeshT=Mesh<Dim, ManifoldDim>;
  MeshT mesh;
  using Elem = typename MeshT::Elem; 
  read_mesh("../data/beam-tri.MFEM", mesh);






  NodeToElem<Elem> node2elem3(mesh);
  auto node2elem=node2elem3.val();

  for(Integer nn=0; nn<node2elem.size();nn++)
    {std::cout<<"node=="<<nn<<"     ";
  for(Integer mm=0; mm<node2elem[nn].size();mm++)
    std::cout<<node2elem[nn][mm]<<" ";
  std::cout<<std::endl;
}
ElemEntity<Elem,0> nodes(mesh,node2elem);

std::cout<<"node 2 elem size=="<<nodes.entity_2_elem().size()<<std::endl;
std::cout<<"node 2 elem=="<<std::endl;
for(Integer nn=0; nn<nodes.entity_2_elem().size();nn++)
 std::cout<<nodes.entity_2_elem(nn)[0]<<" ";
std::cout<<std::endl;
for(Integer nn=0; nn<nodes.elem_2_entity().size();nn++)
{
 std::cout<<"elem="<<nn<< " made of nodes "<< std::endl;
 for(Integer mm=0; mm<nodes.elem_2_entity(nn).size();mm++)
  std::cout<<nodes.elem_2_entity(nn)[mm]<<" ";
std::cout<< std::endl;
} 


using FSspace= FunctionSpace< MeshT, Lagrange2<2>,RT0<1>>;
FSspace FEspace(mesh);


const auto& P2_ens0=ElemEntity<Elem,ElementFunctionSpace<Elem,LagrangeFE,2>::entity[0]>(mesh,node2elem);
const auto& P2_ens1=ElemEntity<Elem,ElementFunctionSpace<Elem,LagrangeFE,2>::entity[1]>(mesh,node2elem);

auto ens2elem20=P2_ens0.entity_2_elem();
auto ens2elem21=P2_ens1.entity_2_elem();

auto elem2ens20=P2_ens0.elem_2_entity();
auto elem2ens21=P2_ens1.elem_2_entity();

std::cout<<"ens2elem 2 0="<< std::endl;
for(Integer nn=0;nn<ens2elem20.size();nn++)
{
  for(Integer mm=0;mm<ens2elem20[nn].size();mm++)
    std::cout<<ens2elem20[nn][mm]<<" ";
  std::cout<<std::endl;
} 
std::cout<<"ens2elem 2 1="<< std::endl;
for(Integer nn=0;nn<ens2elem21.size();nn++)
{
  for(Integer mm=0;mm<ens2elem21[nn].size();mm++)
    std::cout<<ens2elem21[nn][mm]<<" ";
  std::cout<<std::endl;
} 
std::cout<<"elem2ens20 2 0="<< std::endl;
for(Integer nn=0;nn<elem2ens20.size();nn++)
{
  for(Integer mm=0;mm<elem2ens20[nn].size();mm++)
    std::cout<<elem2ens20[nn][mm]<<" ";
  std::cout<<std::endl;
} 
std::cout<<"elem2ens21 2 1="<< std::endl;
for(Integer nn=0;nn<elem2ens21.size();nn++)
{
  for(Integer mm=0;mm<elem2ens21[nn].size();mm++)
    std::cout<<elem2ens21[nn][mm]<<" ";
  std::cout<<std::endl;
} 

std::cout<<"n_dofs="<<FEspace.n_dofs()<< std::endl;
std::cout<<std::endl;
for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
{
 auto &elem_id=elem_iter;
 std::cout<<std::endl<<" elem =="<<elem_iter<<std::endl;
 for(Integer nn=0;nn<FEspace.dofmap(elem_id).size();nn++)
 {
  std::cout<<FEspace.dofmap(elem_id)[nn]<< "  ";
}
} 


   // FEspace.set_new_start(4);
   //  std::cout<<"dofmap_new_start n_dofs="<<FEspace.n_dofs()<< std::endl;
   //  std::cout<<std::endl;
   //  for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
   //  {
   //   auto &elem_id=elem_iter;
   //   std::cout<<std::endl<<" elem =="<<elem_iter<<std::endl;
   //   for(Integer nn=0;nn<FEspace.dofmap_new_start(elem_id).size();nn++)
   //   {
   //      std::cout<<FEspace.dofmap_new_start(elem_id)[nn]<< "  ";
   //   }
   //  } 


std::cout<<std::endl;
for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
 std::cout<<"components of space["<<ss<<"]=="<<FEspace.components(ss)<<std::endl;
std::cout<<std::endl;

for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
{
 std::cout<<"dofs of space ="<<ss<<std::endl;
 for(Integer cc=0;cc<FEspace.components(ss);cc++)
  {std::cout<<"component ="<<cc<<"   "<<std::endl;
auto& vec=FEspace.space_dofs(ss,cc);
for(Integer mm=0;mm<FEspace.n_dofs(ss,cc);mm++)
  std::cout<<vec[mm]<<" ";
}
std::cout<<std::endl;


}


for(Integer mm=0;mm<FEspace.offset().size();mm++)
{
  std::cout<<"OFFSET space ="<<mm<<std::endl;
  for(Integer nn=0;nn<FEspace.offset()[mm].size();nn++)
  {
    std::cout<< FEspace.offset()[mm][nn]<<" ";
  }
  std::cout<<std::endl;
}
std::cout<<std::endl;


auto spaceinf=FEspace.space_info();
for(Integer mm=0;mm<spaceinf.size();mm++)
{
  std::cout<<"Space=="<<mm<<std::endl;
  for(Integer nn=0;nn<spaceinf[mm].size();nn++)
    std::cout<<spaceinf[mm][nn]<<" ";
  std::cout<<std::endl;
}

std::cout<<"Whole dofmap=="<<std::endl;
for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
  {const auto& dm=FEspace.dofmap(0,elem_iter);
    for(auto& i:dm)
     std::cout<<i<<" ";
   std::cout<<std::endl;
 }

 std::cout<<"First Space, first component, dofmap=="<<std::endl;

 for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
  {const auto& dm=FEspace.dofmap(0,0,elem_iter);
    for(auto& i:dm)
     std::cout<<i<<" ";
   std::cout<<std::endl;
 }
 std::cout<<"First Space, second component, dofmap=="<<std::endl;
 for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
  {const auto& dm=FEspace.dofmap(0,1,elem_iter);
    for(auto& i:dm)
     std::cout<<i<<" ";
   std::cout<<std::endl;
 }

 std::cout<<"Second Space, first component, dofmap=="<<std::endl;
 for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
  {const auto& dm=FEspace.dofmap(1,0,elem_iter);
    for(auto& i:dm)
     std::cout<<i<<" ";
   std::cout<<std::endl;
 }
}

void normals_example2D()
{
  constexpr Integer ManifoldDim=2;
  constexpr Integer Dim=2;   
  using MeshT=Mesh<Dim, ManifoldDim>;
  MeshT mesh;
  using Elem = typename MeshT::Elem;   
  read_mesh("../data/beam-tri.MFEM", mesh);
  mark_boundary(mesh);
  auto sn=SignedNormal<Simplex<Dim,ManifoldDim>>(mesh);
  auto n=sn();
  auto alpha=sn.sign();
  sn.print(mesh);
}

void normals_example3D()
{
  constexpr Integer ManifoldDim=3;
  constexpr Integer Dim=3;   
  using MeshT=Mesh<Dim, ManifoldDim>;
  MeshT mesh;
  using Elem = typename MeshT::Elem;   
  read_mesh("../data/tetrahedron_2.MFEM", mesh); 
  mark_boundary(mesh);
  auto sn=SignedNormal<Simplex<Dim,ManifoldDim>>(mesh);
  // auto n=sn();
  // auto alpha=sn.sign();
  sn.print(mesh);
}

void simplex_reference_normals_3D()
{
  constexpr Integer ManifoldDim=3;
  constexpr Integer Dim=3;   
  using MeshT=Mesh<Dim, ManifoldDim>;
  MeshT mesh;
  using Elem = typename MeshT::Elem;   
  read_mesh("../data/tetrahedron_1.MFEM", mesh); 
  std::array<Vector<Real,Dim>, ManifoldDim+1> points;
  Vector<Vector<Real,Dim>, ManifoldDim+1> elempoints;
  Simplex<Dim,ManifoldDim> simplex_elem;
  for(Integer nn=0;nn<ManifoldDim+1;nn++)
    simplex_elem.nodes[nn]=nn;
  for(Integer mm=0;mm<points.size();mm++)
    points[mm]=mesh.point(mm);
    
  elempoints[0]=points[0];
  elempoints[1]=points[1];
  elempoints[2]=points[2];
  elempoints[3]=points[3];
 auto nn = normal(simplex_elem, elempoints);
 for(Integer mm=0;mm<nn.size();mm++)
  std::cout<<nn[mm]<<std::endl;
  std::cout<<std::endl;


  elempoints[0]=points[0];
  elempoints[1]=points[1];
  elempoints[2]=points[3];
  elempoints[3]=points[2];
 nn = normal(simplex_elem, elempoints);
 for(Integer mm=0;mm<nn.size();mm++)
  std::cout<<nn[mm]<<std::endl;
  std::cout<<std::endl;

  elempoints[0]=points[0];
  elempoints[1]=points[2];
  elempoints[2]=points[3];
  elempoints[3]=points[1];
 nn = normal(simplex_elem, elempoints);
 for(Integer mm=0;mm<nn.size();mm++)
  std::cout<<nn[mm]<<std::endl;
  std::cout<<std::endl;

  elempoints[0]=points[1];
  elempoints[1]=points[2];
  elempoints[2]=points[3];
  elempoints[3]=points[0];
 nn = normal(simplex_elem, elempoints);
 for(Integer mm=0;mm<nn.size();mm++)
  std::cout<<nn[mm]<<std::endl;
 std::cout<<std::endl;
}

void normals_example4D()
{
  constexpr Integer ManifoldDim=4;
  constexpr Integer Dim=4;   
  using MeshT=Mesh<Dim, ManifoldDim>;
  MeshT mesh;
  using Elem = typename MeshT::Elem;   
  read_mesh("../data/pentatope_2.MFEM", mesh);
  auto sn=SignedNormal<Simplex<Dim,ManifoldDim>>(mesh);
  auto n=sn();
  auto alpha=sn.sign();
  for(Integer e=0;e<mesh.n_elements();e++)
  {
   std::cout<<"elem=="<<e<<std::endl; 
   std::cout<<"normals:"<<std::endl; 
   for(Integer f=0;f<n[e].size();f++)
     n[e][f].describe(std::cout);
   std::cout<<"signs:"<<std::endl;
   std::cout<<alpha[e]<<std::endl; 
 }    
 sn.print(mesh);

}



void functionspaces_example4D()
{

  constexpr Integer ManifoldDim=4;
  constexpr Integer Dim=4;
  constexpr Integer FEFamily=LagrangeFE;
  constexpr Integer Order=3;
  using MeshT=Mesh<Dim, ManifoldDim>;
  MeshT mesh;
  using Elem = typename MeshT::Elem; 
  read_mesh("../data/pentatope_2.MFEM", mesh);
    //read_mesh("../data/cube_6.MFEM", mesh);
    //read_mesh("../data/square_2_def.MFEM", mesh);


  NodeToElem<Elem> node2elem3(mesh);
  auto node2elem=node2elem3.val();
    //auto vec=dofmap<ElemLagrange1<Elem>   , ElemLagrange3<Elem,1,1> >(mesh);

    // Primitive minusone(-1);
    // Primitive zero(0);
    // Composite first(1); 
    // Composite second(2); 
    // Composite third(3);
    // Composite fifth(5);  

    // first.add(zero);//std::make_shared<decltype(zero)>(zero));
    // third.add(minusone); 
    // third.add(zero); 
    // first.add(third);//std::make_shared<decltype(third)>(third));
    // second.add(first);
    // second.traverse(); 
    // fifth.add(minusone,zero,second,third);
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth.traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth[0]->traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth[1]->traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;
    // fifth[2]->traverse();
    // std::cout<<std::endl<<"-----------------------"<<std::endl;

    // auto fourth=second[0];
    // fourth.traverse(); 
    // cout<<endl<<" nnn = "<< *nnn<<endl;
     // Leaf a0(3);
     // Leaf a1(1);
     // Leaf a2(2);
     // Tree t1(-10);
     // Tree t2(-20);
     // a0.print();
     // a1.print();
     // a2.print();
     // a2[0].print();
     // t2.add(a0);
     // t2.add(a1);

     // t1.add(t2);
     // t1.add(a0);
     // t1.add(a1);
     // t1.add(a2);

     // auto ecco=t1[0][0];
     // ecco.print();
    // ecco.print();
    // tree t1;

    // t1.add(a1);
    // auto ecc=t1[0];
    // std::cout<<"-----------------------"<<a1.val()<<"   "<<ecc->val()<<std::endl;

    // a1.print();
    // auto c1=t1[0];
    // std::cout<<"MA CI ARRIVO?"<<std::endl;
    // std::cout<<a1.val()<<" "<<a2.val()<<" "<<t3.val_vec().size()<<" "<<t3.val_vec()[0]->val()<<std::endl;
    // // constexpr auto n_spaces=2;
    // constexpr auto dofs_per_elem=DofsPerElemNums1<Elem,RT0<1>,Lagrange3<1>>::value;
    // std::vector<std::array<Integer,dofs_per_elem>> dofmap_vec;
    // std::array<std::vector<Integer>,n_spaces> offset;
    // dofmap1<RT0<1>,Lagrange3<1>>(mesh,dofmap_vec,offset);


  FunctionSpace< MeshT, Lagrange1<2>, Lagrange2<1> > FEspace(mesh);

    // auto eeh=FunctionSpaceSystem(FEspace);


  std::cout<<"n_dofs="<<FEspace.n_dofs()<< std::endl;
  std::cout<<std::endl;
  for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
  {
   auto &elem_id=elem_iter;
   std::cout<<"elem_id="<<elem_id<<", number of dofs=s"<<FEspace.dofmap(elem_id).size()<<std::endl;
   for(Integer nn=0;nn<FEspace.dofmap(elem_id).size();nn++)
   {
    std::cout<<FEspace.dofmap(elem_id)[nn]<<" ";
  }
  std::cout<<std::endl;
} 
std::cout<<std::endl;
for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
 std::cout<<"components of space["<<ss<<"]=="<<FEspace.components(ss)<<std::endl;
std::cout<<std::endl;

for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
{
 std::cout<<"dofs of space ="<<ss<<std::endl;
 for(Integer cc=0;cc<FEspace.components(ss);cc++)
  {std::cout<<"component ="<<cc<<"   "<<std::endl;
auto& vec=FEspace.space_dofs(ss,cc);
for(Integer mm=0;mm<FEspace.n_dofs(ss,cc);mm++)
  std::cout<<vec[mm]<<" ";
std::cout<<std::endl;
}
std::cout<<std::endl;

}


for(Integer mm=0;mm<FEspace.offset().size();mm++)
{
  std::cout<<"offset space="<<mm<<std::endl;
  for(Integer nn=0;nn<FEspace.offset()[mm].size();nn++)
  {
    std::cout<< FEspace.offset()[mm][nn]<<" ";
  }
}
std::cout<<std::endl;


for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
  {std::cout<<std::endl;
    const auto size=FEspace.dofmap(1,elem_iter).size();
    std::cout<<"elem_iter="<<elem_iter<<std::endl;
    for(Integer nn=0;nn<size;nn++)
      std::cout<<FEspace.dofmap(1,elem_iter)[nn]<<" ";
    std::cout<<std::endl;
  }

  std::cout<<std::endl;
}





void connectivity_example5D()
{

  constexpr Integer ManifoldDim=4;
  constexpr Integer Dim=4;
  using MeshT=Mesh<Dim, ManifoldDim>;
  MeshT mesh;
  using Elem = typename MeshT::Elem; 
  read_mesh("../data/pentatope_2.MFEM", mesh);
    //read_mesh("../data/cube_6.MFEM", mesh);
    //read_mesh("../data/square_2_def.MFEM", mesh)    


  NodeToElem<Elem> node2elem3(mesh);
  const auto node2elem=node2elem3.val();
  const auto const_entities_tuple=EntitiesOfFunctionSpace<Elem,GeneralSpace,0>(mesh,node2elem);
  const auto node=std::get<0>(const_entities_tuple);
  const auto edge=std::get<1>(const_entities_tuple);
  const auto triangle=std::get<2>(const_entities_tuple);
  constexpr Integer entity_from_index=7;
  constexpr Integer entitydim_from=2;
  constexpr Integer subentitydim_from=1;
  constexpr Integer entitydim_to=1;

  Integer entity_e[3];
  Integer entity_t[3];
  ElemConnectivity<Elem,entitydim_from,subentitydim_from,entitydim_to> conn_e2t(mesh,node2elem);



  const auto &entity_2_elem_e=edge.entity_2_elem();
  const auto &elem2entity_e=edge.elem_2_entity();

  const auto &entity_2_elem_t=triangle.entity_2_elem();
  const auto &elem2entity_t=triangle.elem_2_entity();

  cout<<"ENTITY FROM #elems="<<triangle.size()<<endl;
  cout<<"ENTITY FROM dimension="<<entitydim_from<<endl;

  const auto & elemii=mesh.elem(entity_2_elem_t[entity_from_index][0]);
  Combinations<ManifoldDim + 1, triangle.num_of_points()>::generate(entity_2_elem_t[entity_from_index][1],entity_t);
  for(int jj=0;jj<triangle.num_of_points();jj++)
   cout<<elemii.nodes[entity_t[jj]]<<" ";    
 cout<<endl;


 cout<<"SUBENTITY FROM dimension="<<subentitydim_from<<endl;
 cout<<"ENTITY TO dimension="<<entitydim_to<<endl;
 cout<<"ENTITY FROM of interest family="<<entity_from_index<<endl;

 

 
 cout<<"ENTITY TO #elems="<<edge.size()<<endl;
 for(int ii=0;ii<edge.size();ii++)
 {

  cout<<endl;
  for(int jj=0;jj<1;jj++)
   cout<<entity_2_elem_e[ii][jj]<<" ";
 cout<<"--------ENTITY TO elem id="<<ii<<"   ";
 const auto & elemii=mesh.elem(entity_2_elem_e[ii][0]);
 Combinations<ManifoldDim + 1, edge.num_of_points()>::generate(entity_2_elem_e[ii][1],entity_e);
 for(int jj=0;jj<edge.num_of_points();jj++)
   cout<<elemii.nodes[entity_e[jj]]<<" ";    
}   

const auto & connection=conn_e2t.compute(triangle,entity_from_index,edge);

cout<<endl;
for(Integer ii=0;ii<connection.size();ii++)
 cout<<" CONNECTIONS: ="<<connection[ii]<<endl;



}


}



#endif




// Base Function: f(x) = x
// class Function
// {
//     protected:
//     struct Implementation
//     {
//         virtual ~Implementation() {}
//         virtual double evaluate(double x) const { return x; }
//     };

//     public:
//     Function()
//     :   self_(std::make_shared<Implementation>())
//     {}

//     double operator () (double x) const { return self_->evaluate(x); }

//     protected:
//     Function(std::shared_ptr<Implementation> self)
//     :   self_(self)
//     {}

//     private:
//     std::shared_ptr<Implementation> self_;
// };



// class Function2: public Function
// {
//      protected:
//     struct Implementation: Function::Implementation
//     {
//         Function f;
//         Implementation(Function f):   
//         f(f)
//         {};
//         virtual double evaluate(double x) const override{return 3.333*x; }
//     };  
// public: 
//     // double operator () (double x) const  { return 2*x; };
//     // Function2(Function f)
//     // :   Function(std::make_shared<Implementation>(f))
//     // {}
//     Function2()  
//     :   Function(std::make_shared<Implementation>(Function()))
//     {};
// };

// // Unary Function: u(-f(x))
// class UnaryMinus : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;


//         Implementation(Function f):   
//         f(f)
//         {};

//         virtual double evaluate(double x) const override { return -f(x); }
//     };

//     public:
//     UnaryMinus(Function f)
//     :   Function(std::make_shared<Implementation>(f))
//     {}
// };

// // Binary Function: u(f(x) + g(x))
// class BinaryAdd : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;
//         Function g;
//         Implementation(Function f, Function g)
//         :   f(f), g(g)
//         {};

//         virtual double evaluate(double x) const override { return f(x) + g(x); }
//     };

//     public:
//     BinaryAdd(Function f, Function g)
//     :   Function(std::make_shared<Implementation>(f, g))
//     {}
// };

// // Binary Function: u(f(x) * g(x))
// class BinaryMultiply : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;
//         Function g;
//         Implementation(Function f, Function g)
//         :   f(f), g(g)
//         {};

//         virtual double evaluate(double x) const override { return f(x) * g(x); }
//     };

//     public:
//     BinaryMultiply(Function f, Function g)
//     :   Function(std::make_shared<Implementation>(f, g))
//     {}
// };

// // Binary Function: u(f(x) * g(x))
// class scalarMultiply : public Function
// {
//     protected:
//     struct Implementation : Function::Implementation
//     {
//         Function f;
//         double alpha_;
//         Implementation(Function f, double alpha):   
//         f(f), 
//         alpha_(alpha)
//         {};


//         virtual double evaluate(double x) const override { return f(x) * alpha_; }
//     };

//     public:
//     scalarMultiply(Function f, double alpha)
//     :   Function(std::make_shared<Implementation>(f, alpha))
//     {}

// };

// inline scalarMultiply operator * (Function f,double alpha) { return scalarMultiply(f,alpha); }
// inline scalarMultiply operator * (double alpha,Function f) { return scalarMultiply(f,alpha); }
// inline UnaryMinus operator - (Function f) { return UnaryMinus(f); }
// inline BinaryAdd operator + (Function f, Function g) { return BinaryAdd(f, g); }
// inline BinaryMultiply operator * (Function f, Function g) { return BinaryMultiply(f, g); }
