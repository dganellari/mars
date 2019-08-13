#ifndef MARS_EXAMPLES_HPP
#define MARS_EXAMPLES_HPP

#include "mars_simplex.hpp"
#include "mars_connectivity.hpp"
#include "mars_functionspace_dofmap.hpp"
#include "mars_functionspace.hpp"
#include "mars_shape_function.hpp"
#include "mars_matrix.hpp"



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
  using MeshT=mars::Mesh<Dim, ManifoldDim>;  
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
  using MeshT=mars::Mesh<Dim, ManifoldDim>;
  using MeshT2=mars::Mesh<Dim+1, ManifoldDim+1>;
  MeshT2 mesh2;
  read_mesh("../data/beam-tet.MFEM", mesh2);
 

  MeshT mesh;
  using Elem = typename MeshT::Elem; 
  read_mesh("../data/beam-tri.MFEM", mesh);

  constexpr Integer QPOrder=4;
  constexpr Integer NQPoints=GaussPoints<Elem,QPOrder>::NQPoints; 
  using QP=typename GaussPoints<Elem,QPOrder>::qp_points_type;
  GaussPoints<Elem,QPOrder> gauss;

  Elem elem;
  Matrix<Real, Dim, ManifoldDim> J;

  constexpr Integer Npoints=Elem::Npoints;
  std::vector<Vector<Real,Dim>> points(Npoints);

  for(Integer ii=0;ii<Npoints;ii++)
   elem.nodes[ii]=ii;

 points[0][0]=0;
 points[0][1]=0;
 points[1][0]=2;
 points[1][1]=3;
 points[2][0]=1;
 points[2][1]=4;
 jacobian(elem,points,J);

 using FSspace1= FunctionSpace< MeshT, Lagrange2<1>,RT0<1>>;
 FSspace1 FEspace1(mesh);
 using FSspace2= FunctionSpace< MeshT, Lagrange2<1>,RT0<1>>;
 FSspace2 FEspace2(mesh);
 using FSspace3= FunctionSpace< MeshT, Lagrange1<1>, RT1<1>,Lagrange2<1>>;
 FSspace3 FEspace3(mesh);
 using FSspace4= FunctionSpace< MeshT, Lagrange2<1>,Lagrange1<1>>;
 FSspace4 FEspace4(mesh);
 using FSspace5= FunctionSpace< MeshT, Lagrange2<1>, Lagrange3<2>>;

 using FSspace10= FunctionSpace< MeshT, Lagrange1<1>>;
 FSspace10 FEspace10(mesh);
 using FSspace11= FunctionSpace< MeshT, RT0<1>>;
 FSspace11 FEspace11(mesh);

 using FSspaceT2= FunctionSpace< MeshT2, Lagrange1<1>>;
 FSspaceT2 FEspaceT2(mesh2);

// auto ecchime=FEspace3.set_quadrature_rule<0,GaussPoints<Elem,2>,GaussPoints<Elem,2>>();
 using fespace1= ElemFunctionSpace<Simplex<2,2>, Lagrange1<1>>;
 using fespace2= ElemFunctionSpace<Simplex<2,2>, Lagrange2<2>>;
 using fespace3= ElemFunctionSpace<Simplex<2,2>, RT0<1>>;
 
 BilinearFormIntegrator<4,FSspace1,FSspace1,GradientOperator,IdentityOperator>(FEspace1,FEspace1);



 auto W=MixedFunctionSpace(FEspace3,FEspace4);
 auto W1=MixedFunctionSpace(W,FEspace4);

decltype(W1)::ElementFunctionSpacesTupleType lkj;


 constexpr Integer NComponents=2;

 using QR=GaussPoints<Elem,QPOrder>;
 ShapeFunctionOperatorDependent<Simplex<2,2>, Lagrange1<NComponents>,IdentityOperator,QR >  sfod;
 sfod.init_reference();
 sfod.init(J);
 sfod.function();
 std::cout<<"sfod.function()=="<<sfod.function()<<std::endl;





 auto W3=MixedFunctionSpace(FEspace10,FEspace11);


 ShapeFunctionExpression<0,2,Lagrange1<1>> sfe1; 
 GradientShapeFunctionExpression<0,2,Lagrange1<1>> sfe2;
 ShapeFunctionExpression<1,2,RT0<1>> sfe3;


 auto sfen=sfe1+sfe2+sfe3;

 using TupleOperatorsAndQuadrature=TupleOfTupleChangeType<1,QR,OperatorTupleType<decltype(sfen)>::type>; 
 using TupleSpaces=typename decltype(W3)::UniqueElementFunctionSpacesTupleType;


 TupleSpaces reerree;
 TupleOperatorsAndQuadrature erre;



 auto u =     MakeTrial<0>(W1);
 auto sigma = MakeTrial<1>(W1);
 auto p =     MakeTrial<2>(W1);
 auto r =     MakeTrial<3>(W1);
 auto rr1 =     MakeTrial<4>(W1);
 auto rr2 =     MakeTrial<5>(W1);
 auto rr3 =     MakeTrial<6>(W1);

 auto v =   MakeTest<0>(W1);
 auto tau = MakeTest<1>(W1);
 auto q =   MakeTest<2>(W1);
 auto s =   MakeTest<3>(W1);
 auto ss1 =   MakeTest<4>(W1);
 auto ss2 =   MakeTest<5>(W1);
 auto ss3 =   MakeTest<6>(W1);

// bilinear form
 auto l22=L2Inner(mesh,2*Div(sigma),Div(tau))+
 L2Inner(mesh,Grad(u),Grad(v))
 -
 L2Inner(mesh,u,Div(tau))+
 L2Inner(mesh,rr1,q)+
 L2Inner(mesh,r,s)+
 L2Inner(mesh,Grad(rr1)+sigma,tau)+
 L2Inner(mesh,r,s);

 auto shape_coefficients=shape_function_coefficients(l22);
 shape_coefficients.init(mesh);
 auto referencemaps=reference_maps(l22);
 auto shapefunctions=shape_functions(l22);


 for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
 {
  elem=mesh.elem(elem_iter);
  // mesh.points(elem_iter,points);
  jacobian(elem,mesh.points(),J);
  shape_coefficients.init(elem_iter);
  referencemaps.init(J);
  shapefunctions.init_map(referencemaps);
  shapefunctions.init(shape_coefficients);
  std::cout<<"elem_iter=="<<elem_iter<<std::endl;
  for(Integer mm=0;mm<elem.nodes.size();mm++)
    std::cout<<elem.nodes[mm]<<std::endl;
  std::cout<<J<<std::endl;
}


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
 auto l2rs=L2Inner(mesh,r,s);
 auto eval2=Eval(l2rs,shapefunctions);
 // evals.apply();
 // eval2.apply(shapefunctions());


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

std::cout<<"mat==="<<mat3<<std::endl;

std::cout<<" elapsed_secs------------>"<<elapsed_secs<<std::endl;

std::cout<<" nonzero------------>"<<nonzero<<std::endl;

auto& dm =W1.dofmap();
auto& dm0=std::get<0>(dm);
auto& dm01=std::get<0>(dm0);
auto& dm02=std::get<1>(dm0);
auto& dm1=std::get<1>(dm);
// auto& dm2=std::get<2>(dm);




std::cout<<"dm01"<<std::endl;
for(Integer ii=0;ii<dm01.size();ii++)
{
 for(Integer jj=0;jj<dm01[ii].size();jj++)
  std::cout<<dm01[ii][jj]<<" ";
std::cout<<std::endl;
}
std::cout<<"dm02"<<std::endl;

for(Integer ii=0;ii<dm02.size();ii++)
{
 for(Integer jj=0;jj<dm02[ii].size();jj++)
  std::cout<<dm02[ii][jj]<<" ";
std::cout<<std::endl;
}
std::cout<<"dm1"<<std::endl;

for(Integer ii=0;ii<dm1.size();ii++)
{
 for(Integer jj=0;jj<dm1[ii].size();jj++)
  std::cout<<dm1[ii][jj]<<" ";
std::cout<<std::endl;
}
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
  using MeshT=mars::Mesh<Dim, ManifoldDim>;
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
  using MeshT=mars::Mesh<Dim, ManifoldDim>;
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
  using MeshT=mars::Mesh<Dim, ManifoldDim>;
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
  using MeshT=mars::Mesh<Dim, ManifoldDim>;
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
  using MeshT=mars::Mesh<Dim, ManifoldDim>;
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
  using MeshT=mars::Mesh<Dim, ManifoldDim>;
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
  using MeshT=mars::Mesh<Dim, ManifoldDim>;
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
