	#ifndef MARS_EXAMPLES_HPP
	#define MARS_EXAMPLES_HPP


	#include <fstream>
	#include <sstream>
	#include "mars_simplex.hpp"
	#include "mars_connectivity.hpp"
	#include "mars_functionspace_dofmap.hpp"
	#include "mars_functionspace.hpp"
	#include "mars_shape_functions_collection.hpp"

	#include "mars_shape_function.hpp"
	#include "mars_matrix.hpp"
	#include "mars_jacobian.hpp"




	#include "mars_evaluation.hpp"
	#include "mars_quadrature_order.hpp"
	#include "mars_function.hpp"
	#include "mars_general_form.hpp"
	#include "mars_l2_dot_product_integral.hpp"
	#include "mars_test_function.hpp"
	#include "mars_trial_function.hpp"
	#include "mars_shape_functions_collection.hpp"

	#include "generation/mars_mesh_generation.hpp"

	#include "mars_context.hpp"

	#include "mars_dirichlet_bc.hpp"

	#include "mars_shape_function.hpp"


	namespace mars{


		using std::cout;
		using std::endl;







	//  constexpr Integer simplex_face_sub_entities(const Integer& SimplexDim,const Integer& FaceNumber,const Integer& SubEntityDim, const Integer& SubEntityDimNumber)
	//  {
	//   switch(SimplexDim)
	//   {
	//    // triangles
	//    case 2:
	//    // triangle faces (edges)
	//    switch(FaceNumber)
	//    {
	//     // face 0
	//     case 0:
	//        switch(SubEntityDim)
	//        {
	//         // face 0, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 1;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 0, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }        
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};

	//        }

	//     case 1:
	//        switch(SubEntityDim)
	//        {
	//         // face 0, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 2;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 0, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 1;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }        
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};

	//        }
	//     case 2:
	//        switch(SubEntityDim)
	//        {
	//         // face 2, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 1;
	//           case 1: return 2;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 2, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 2;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }        
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};

	//        }
	//     default: {assert(0 &&"simplex_face_sub_entities: invalid face number");return -1;};
	//    }
	//   // tetrahedrons
	//   case 3:
	//   // tetrahedrons faces (triangles)
	//    switch(FaceNumber)
	//    {
	//     // face 0
	//     case 0:
	//        switch(SubEntityDim)
	//        {
	//         // face 0, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 1;
	//           case 2: return 2;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 0, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 1;
	//           case 2: return 3;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 0, triangles
	//         case 2:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }         
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
	//        }

	//     case 1:
	//        switch(SubEntityDim)
	//        {
	//         // face 1, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 1;
	//           case 2: return 3;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 1, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 2;
	//           case 2: return 4;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 1, triangles
	//         case 2:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 1;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }         
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
	//        }
	//     case 2:
	//        switch(SubEntityDim)
	//        {
	//         // face 0, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 2;
	//           case 2: return 3;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 0, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 1;
	//           case 1: return 2;
	//           case 2: return 5;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 0, triangles
	//         case 2:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 2;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }         
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
	//        }
	//     case 3:
	//         switch(SubEntityDim)
	//        {
	//         // face 0, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 1;
	//           case 1: return 2;
	//           case 2: return 3;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 0, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 3;
	//           case 1: return 4;
	//           case 2: return 5;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 0, triangles
	//         case 2:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 3;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }         
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
	//        }

	//     default: {assert(0 &&"simplex_face_sub_entities: invalid face number");return -1;};
	//    }

	//   // pentatope
	//   case 4:
	//   // tetrahedrons faces (tetrahedrons)
	//    switch(FaceNumber)
	//    {
	//     // face 0
	//     case 0:
	//        switch(SubEntityDim)
	//        {
	//         // face 0, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 1;
	//           case 2: return 2;
	//           case 3: return 3;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 0, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 1;
	//           case 2: return 2;
	//           case 3: return 4;
	//           case 4: return 5;
	//           case 5: return 7;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 0, triangles
	//         case 2:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 1;
	//           case 2: return 3;
	//           case 3: return 6;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 0, tetrahedrons
	//         case 3:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }        
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
	//        }

	//     // face 1
	//     case 1:
	//        switch(SubEntityDim)
	//        {
	//         // face 1, edges
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 1;
	//           case 2: return 2;
	//           case 3: return 4;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 1, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 1;
	//           case 2: return 3;
	//           case 3: return 4;
	//           case 4: return 6;
	//           case 5: return 8;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 1, triangles
	//         case 2:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 2;
	//           case 2: return 4;
	//           case 3: return 7;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 1, tetrahedrons
	//         case 3:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 1;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }        
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
	//        }
	//     // face 2
	//     case 2:
	//        switch(SubEntityDim)
	//        {
	//         // face 2, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 1;
	//           case 2: return 3;
	//           case 3: return 4;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 2, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 2;
	//           case 2: return 3;
	//           case 3: return 5;
	//           case 4: return 6;
	//           case 5: return 9;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 2, triangles
	//         case 2:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 1;
	//           case 1: return 2;
	//           case 2: return 5;
	//           case 3: return 8;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 2, tetrahedrons
	//         case 3:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 2;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }        
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
	//        }

	//     // face 3
	//     case 3:
	//        switch(SubEntityDim)
	//        {
	//         // face 3, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 0;
	//           case 1: return 2;
	//           case 2: return 3;
	//           case 3: return 4;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 3, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 1;
	//           case 1: return 2;
	//           case 2: return 3;
	//           case 3: return 7;
	//           case 4: return 8;
	//           case 5: return 9;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 3, triangles
	//         case 2:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 3;
	//           case 1: return 4;
	//           case 2: return 5;
	//           case 3: return 9;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 3, tetrahedrons
	//         case 3:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 3;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }        
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
	//        }

	//     // face 4
	//     case 4:
	//        switch(SubEntityDim)
	//        {
	//         // face 4, nodes
	//         case 0:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 1;
	//           case 1: return 2;
	//           case 2: return 3;
	//           case 3: return 4;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid sub entity number");return -1;};
	//         }
	//         // face 4, edges
	//         case 1:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 4;
	//           case 1: return 5;
	//           case 2: return 6;
	//           case 3: return 7;
	//           case 4: return 8;
	//           case 5: return 9;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 4, triangles
	//         case 2:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 6;
	//           case 1: return 7;
	//           case 2: return 8;
	//           case 3: return 9;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }  
	//         // face 4, tetrahedrons
	//         case 3:
	//         switch(SubEntityDimNumber)
	//         {
	//           case 0: return 4;
	//           default: {assert(0 &&"simplex_face_sub_entities: invalid edge number");return -1;};
	//         }        
	//         default: {assert(0 &&"simplex_face_sub_entities: invalid SubEntityDim");return -1;};
	//        }
	//     default: {assert(0 &&"simplex_face_sub_entities: invalid face number");return -1;};
	//    }
	//   default: {assert(0 &&"simplex_face_sub_entities: invalid simplex dimension");return -1;};
	//  }
	// }


	// // works only for simplices
	// template <typename Space>
	// constexpr Integer trace_surf_n_dofs()
	//   { 
	//     using Elem=typename Space::Elem;
	//     Integer n_dofs=0;
	//     for(Integer ii=0;ii<Space::entity.size();ii++)
	//     {
	//       if(Space::entity[ii]<=Elem::ManifoldDim)
	//        n_dofs+=Space::dofs_per_entity[ii]* binomial_coefficient(Elem::ManifoldDim,Space::entity[ii]+1);
	//    }
	//    return n_dofs*Space::NComponents;
	//  }


	// template<typename FunctionSpace>
	// constexpr Integer function_space_dofs_per_elem(const Integer Dim,const Integer N)
	// {   
	//     std::size_t dofs_per_elem=0;
	//     if(N>0)
	//     dofs_per_elem=FunctionSpace::NComponents * 
	//                                      FunctionSpace::dofs_per_entity[N-1]   * 
	//                                      binomial_coefficient(Dim+1,FunctionSpace::entity[N-1]+1);  

	//      switch(N)
	//      {
	//       case  -1: return 0;
	//       case  0: return 0;
	//       case  1: return dofs_per_elem;
	//       default: return function_space_dofs_per_elem<FunctionSpace>(Dim,N-1)+dofs_per_elem;
	//      }
	// };


	// template <typename Space>
	// constexpr auto trace_dofs(const Integer face)
	// {


	//  constexpr auto n=trace_surf_n_dofs<Space>();
	//  const auto ManifoldDim=Space::ManifoldDim;
	//  const auto NComponents=Space::NComponents;
	//  const auto entity=Space::entity;
	//  const auto dofs_per_entity= Space::dofs_per_entity;

	//  Array<Integer, n> dofs;
	//  Integer cont=0;
	//  Integer dofs_per_entity_cont=0;
	//  // loop on all the kind of dofs-entities
	//  for(Integer ii=0; ii< entity.size(); ii++)
	//  {
	//   const auto entityii=entity[ii];
	//   // consider the ii-th dof-entity only if it smaller than the dimension of the manifolddim
	//   // TODO FIXME IT SHOULD BE entity[ii]<ManifoldDim
	//   if(entity[ii]<ManifoldDim)
	//   {
	//     const auto dofs_per_entityii=dofs_per_entity[ii];
	//     const auto binomial_coeff=binomial_coefficient(ManifoldDim,entity[ii]+1);
	//     // loop on the entityes of entityii dim
	//     for(Integer jj=0; jj<binomial_coeff; jj++)
	//       {
	//         // loop on the dofs of the given entity
	//         dofs[cont] = NComponents*dofs_per_entityii*simplex_face_sub_entities(ManifoldDim,face,entityii,jj)+function_space_dofs_per_elem<Space>(ManifoldDim,ii);
	//         cont++;
	//         for(Integer ss=1;ss<NComponents;ss++)
	//         { 
	//           dofs[cont] = dofs[cont-1]+1;
	//           cont++;
	//         }

	//         for(Integer kk=1; kk<dofs_per_entityii; kk++)
	//           {
	//                for(Integer ss=0;ss<NComponents;ss++)
	//                 { 
	//                   dofs[cont] = dofs[cont-1]+1;
	//                   cont++;
	//                 }
	//           }
	//       }
	//   }
	//  }
	//  return dofs;
	// }


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


	//       template<typename Space,typename Operator,typename single_type,
	//                Integer NQPoints,Integer Dim>
	//        constexpr auto reference_trace_shape_function_init(const Integer face, const Matrix<Real,NQPoints,Dim>&qp_points)
	//        {

	//         Array<Real,Dim> qp_point;
	//         const auto dofs=trace_dofs<Space>(face);
	//         const auto n_dofs=dofs.size();
	//         Array<Array<single_type,NQPoints>,n_dofs> v;
	//         Array<single_type,n_dofs> func;

	//             for(Integer qp=0;qp<NQPoints;qp++)
	//             {
	//              qp_point=qp_points.get_row(qp);
	//              // func=value<Elem,Operator,FEFamily,Order,single_type,Ndofs>(qp_point);
	//              ReferenceShapeFunctionValue<typename Space::Elem,Operator,Space::FEFamily,Space::Order>::apply(qp_point,func);
	//              // value<Space::Elem,Operator,Space::FEFamily,Space::Order>(qp_point,func);
	//               for(Integer n_dof = 0; n_dof < n_dofs; ++n_dof) {
	//                   const_cast<single_type&>
	//                   (static_cast<const std::array<single_type,NQPoints>& >
	//                    ((static_cast<const std::array<Array<single_type,NQPoints>,n_dofs>& >(v())[n_dof] )())[qp])=
	//                   static_cast<const std::array<single_type,n_dofs>& >(func())[n_dof];
	//               }
	//             }
	//        return v;
	//       };

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


	//    template<Integer Dim, Integer ManifoldDim>
	// inline constexpr auto jacobian_faces()
	// {
	//     static_assert(Dim >= ManifoldDim, "Dim must be greater or equal ManifoldDim");
	//     // static_assert(Npoints == ManifoldDim+1, "Npoints must be equal to ManifoldDim+1");
	//     Array<Vector<Real, Dim>,ManifoldDim> points;
	//     constexpr auto n_faces=ManifoldDim+1;
	//     const auto combs=combinations_generate<ManifoldDim+1,ManifoldDim>(); 
	//     Matrix<Real, Dim, ManifoldDim-1> J;
	//     Vector<Matrix<Real, Dim, ManifoldDim-1>,n_faces> Jmat;
	//     // loop on all the faces
	//     for(Integer ii=0;ii<n_faces;ii++)
	//     {
	//         // take the indices of the reference simplex related to the face ii
	//         const auto &comb_ii=combs[ii];
	//         // fill points with the corresponding reference simplex face 
	//         for(Integer jj=0;jj<ManifoldDim;jj++)
	//           for(Integer kk=0;kk<ManifoldDim;kk++)
	//             points[jj][kk]=Simplex<Dim,ManifoldDim>::reference[comb_ii[jj]][kk];

	//         // compute the jacobian of the given face
	//         for(Integer ii=0;ii<ManifoldDim;ii++)
	//         {
	//             Vector<Real, Dim> v0 = points[0];
	//             for(Integer i = 1; i < ManifoldDim; ++i) {
	//                 const auto &vi = points[i];
	//                 J.col(i-1, vi - v0);
	//             }
	//         }
	//         Jmat[ii]=J;

	//     }
	//     return Jmat;
	// }


	// template<typename QuadratureRule>
	// inline constexpr auto reference_face_shape_functions()
	// {
	//   using Elem=typename QuadratureRule::Elem;
	//   constexpr Integer Dim=Elem::Dim;
	//   constexpr Integer ManifoldDim=Elem::ManifoldDim;
	//   const auto Jacobian_faces=jacobian_faces<Dim,ManifoldDim>();
	//   constexpr auto n_faces=Jacobian_faces.size();
	//   Vector<typename QuadratureRule::qp_points_type, n_faces> qp_points_face;

	//   for(Integer ii=0;ii<n_faces;ii++)
	//      {
	//         qp_points_face[ii]=Jacobian_faces[ii]*QuadratureRule::qp_points_type;
	//      }

	//  return qp_points_face;
	// }

	template<typename Mat,typename Vec>
		void ConjugateGradient(Vec& x,const Mat&A,const Vec& b,const Integer maxiter,const Real toll)
		{

			const int n=b.size();
			Vec r(n);
			Vec p(n);

			if(x.size()!=n)
			{
				x.resize(n);
				for(std::size_t i=0; i<n;i++)
					x[i]=0;
			}

			for(Integer i=0;i<n;i++)
				for(Integer j=0;j<n;j++) 
				{
					r[i]=b[i];
					r[i]-=-A[i][j]*x[j];
					p[i]=r[i];
				}

				for(Integer kk=0;kk<maxiter;kk++)
				{
	    // auto r2=r[0]*r[0];
	    // for(Integer i=1;i<n;i++)
	    //     auto r2+=r[i]*r[i];
	    // Ap=A*p;
	    // pAp=p'*Ap;
	    // alpha=r2/pAp;
	    // x=x+alpha*p;
	    // r=r-alpha*A*p;
	    // xk(kk,:)=x;
	    // rk(kk,:)=r;
	    // pk(kk,:)=p;
	    // normrk(kk)=norm(r);
	    // if(norm(r)<toll)
	    //     break
	    // end
	    // beta=r'*r/r2;
	    // p=r+beta*p;
				}


			}







	template<typename Mat,typename Vec>
			void gauss_seidel(Vec& x, const Mat& A, const Vec& b,const Integer max_iter)
			{
				const int n=b.size();

				if(x.size()!=n)
				{
					x.resize(n);
					for(std::size_t i=0; i<n;i++)
						x[i]=0;
				}
				Real tmp;

				for(std::size_t it=0;it<max_iter;it++)
				{

					A.row_gauss_seidel(x,b);

	  // if(x.size()!=n)
	  //   {
	  //    x.resize(n);
	  //    for(std::size_t i=0; i<n;i++)
	  //       x[i]=0;
	  //    }
	  // Real tmp;

	  // for(std::size_t it=0;it<max_iter;it++)
	  // {
	  // for(std::size_t i=0; i<n;i++)
	  // { 
	  //   tmp=b[i];
	  //   for(std::size_t j=0; j<i;j++)
	  //   {
	  //     // tmp-=A[i][j]*x[j];
	  //     tmp-=A(i,j)*x[j];
	  //   }
	  //   for(std::size_t j=i+1; j<n;j++)
	  //   {
	  //     // tmp-=A[i][j]*x[j];
	  //     tmp-=A(i,j)*x[j];
	  //   }
	  //   // x[i]=(tmp)/A[i][i];
	  //   x[i]=(tmp)/A(i,i);
	  // }
				}

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




				using Space2=ElementFunctionSpace<Simplex<2,2>,LagrangeFE,1,1,1>;
				constexpr const auto trace2=TraceDofs<Space2>::dofs();

				using Space3=ElementFunctionSpace<Simplex<2,2>,LagrangeFE,2,1,1>;
				constexpr const auto trace3=TraceDofs<Space3>::dofs();

				using Space4=ElementFunctionSpace<Simplex<2,2>,LagrangeFE,3,1,1>;
				constexpr const auto trace4=TraceDofs<Space4>::dofs();


				using Space5=ElementFunctionSpace<Simplex<2,2>,RaviartThomasFE,0,1,1>;
				constexpr const auto trace5=TraceDofs<Space5>::dofs();


				using Space6=ElementFunctionSpace<Simplex<2,2>,RaviartThomasFE,1,1,1>;
				constexpr const auto trace6=TraceDofs<Space6>::dofs();





				using Space7=ElementFunctionSpace<Simplex<2,2>,LagrangeFE,1,1,2>;
				constexpr const auto trace7=TraceDofs<Space7>::dofs();

				using Space8=ElementFunctionSpace<Simplex<2,2>,LagrangeFE,2,1,2>;
				constexpr const auto trace8=TraceDofs<Space8>::dofs();

				using Space9=ElementFunctionSpace<Simplex<2,2>,LagrangeFE,3,1,2>;
				constexpr const auto trace9=TraceDofs<Space9>::dofs();


				using Space10=ElementFunctionSpace<Simplex<2,2>,RaviartThomasFE,0,1,2>;
				constexpr const auto trace10=TraceDofs<Space10>::dofs();


				using Space11=ElementFunctionSpace<Simplex<2,2>,RaviartThomasFE,1,1,2>;
				constexpr const auto trace11=TraceDofs<Space11>::dofs();


	  // constexpr auto automa=reference_trace_shape_function_init<Space,IdentityOperator,single_type>(0, GaussPoints<Space::Elem,3>::qp_points_type);
				std::cout<<"Simplex<2,2>,LagrangeFE,Order=1,Continuity=1,Ncomponents=1>"<<std::endl;

				for(int ii=0;ii<trace2.size();ii++)
					{for(int jj=0;jj<trace2[ii].size();jj++)
						std::cout<<trace2[ii][jj]<<std::endl;
						std::cout<<std::endl;
					} 
					std::cout<<"Simplex<2,2>,LagrangeFE,Order=2,Continuity=1,Ncomponents=1"<<std::endl;

					for(int ii=0;ii<trace3.size();ii++)
						{for(int jj=0;jj<trace3[ii].size();jj++)
							std::cout<<trace3[ii][jj]<<std::endl;
							std::cout<<std::endl;
						} 
						std::cout<<"Simplex<2,2>,LagrangeFE,Order=3,Continuity=1,Ncomponents=1"<<std::endl;


						for(int ii=0;ii<trace4.size();ii++)
							{for(int jj=0;jj<trace4[ii].size();jj++)
								std::cout<<trace4[ii][jj]<<std::endl;
								std::cout<<std::endl;
							} 

							std::cout<<"Simplex<2,2>,RaviartThomasFE,Order=0,Continuity=1,Ncomponents=1"<<std::endl;

							for(int ii=0;ii<trace5.size();ii++)
								{for(int jj=0;jj<trace5[ii].size();jj++)
									std::cout<<trace5[ii][jj]<<std::endl;
									std::cout<<std::endl;
								} 

								std::cout<<"Simplex<2,2>,RaviartThomasFE,Order=1,Continuity=1,Ncomponents=1"<<std::endl;

								for(int ii=0;ii<trace6.size();ii++)
									{for(int jj=0;jj<trace6[ii].size();jj++)
										std::cout<<trace6[ii][jj]<<std::endl;
										std::cout<<std::endl;
									} 

									std::cout<<"Simplex<2,2>,LagrangeFE,Order=1,Continuity=1,Ncomponents=2>"<<std::endl;

									for(int ii=0;ii<trace7.size();ii++)
										{for(int jj=0;jj<trace7[ii].size();jj++)
											std::cout<<trace7[ii][jj]<<std::endl;
											std::cout<<std::endl;
										} 
										std::cout<<"Simplex<2,2>,LagrangeFE,Order=2,Continuity=1,Ncomponents=2>"<<std::endl;

										for(int ii=0;ii<trace8.size();ii++)
											{for(int jj=0;jj<trace8[ii].size();jj++)
												std::cout<<trace8[ii][jj]<<std::endl;
												std::cout<<std::endl;
											} 
											std::cout<<"Simplex<2,2>,LagrangeFE,Order=3,Continuity=1,Ncomponents=2>"<<std::endl;

											for(int ii=0;ii<trace9.size();ii++)
												{for(int jj=0;jj<trace9[ii].size();jj++)
													std::cout<<trace9[ii][jj]<<std::endl;
													std::cout<<std::endl;
												} 


												std::cout<<"Simplex<2,2>,RaviartThomasFE,Order=0,Continuity=1,Ncomponents=2>"<<std::endl;

												for(int ii=0;ii<trace10.size();ii++)
													{for(int jj=0;jj<trace10[ii].size();jj++)
														std::cout<<trace10[ii][jj]<<std::endl;
														std::cout<<std::endl;
													} 
													std::cout<<"Simplex<2,2>,RaviartThomasFE,Order=1,Continuity=1,Ncomponents=2>"<<std::endl;

													for(int ii=0;ii<trace11.size();ii++)
														{for(int jj=0;jj<trace11[ii].size();jj++)
															std::cout<<trace11[ii][jj]<<std::endl;
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





















	template<typename T>
		class SparseMatrix
		{
		private:
			std::vector<T> A_;
			std::vector<std::map<Integer,Integer>> cols_idx_; 
			Integer max_rows_;
			Integer max_cols_;
			Integer cont_; 


		public:
			SparseMatrix(){}


			SparseMatrix(const Integer max_rows,const Integer max_cols):
			max_rows_(max_rows),
			max_cols_(max_cols),
			cont_(0)
			{
				A_.resize(max_rows_*max_cols);
				cols_idx_.resize(max_rows_);
			}

			auto rows(){return max_rows_;}
			void init(const Integer max_rows,const Integer max_cols)
			{
				max_rows_=max_rows;
				max_cols_=max_cols;
				cont_=0;
				A_.resize(max_rows_*max_cols);
				cols_idx_.resize(max_rows_);
			}

			void equal(const T& value,const Integer i, const Integer j)
			{
				if(cols_idx_[i].count(j))
				{
					A_[cols_idx_[i].at(j)]=value;
// std::cout<<"already exists in "<< cols_idx_[i].count(j)<<std::endl;
// for(std::size_t i=0;i<A_.size();i++)
//   std::cout<<A_[i]<<" ";
// std::cout<<std::endl;
				}
				else
				{
					A_[cont_]=value ;
// std::cout<<"does not already exists"<<std::endl;
// for(std::size_t i=0;i<A_.size();i++)
// std::cout<<A_[i]<<" ";
// std::cout<<std::endl;
					cols_idx_[i].insert (std::pair<Integer,Integer>(j,cont_) );
					cont_++ ; 
				}

			}


			void plus_equal(const T& value,const Integer i, const Integer j)
			{
// std::cout<<"plus equal"<<std::endl;
// std::cout<<cols_idx_.size()<<std::endl;
				if(cols_idx_[i].count(j))
				{
// std::cout<<"already exists in "<< cols_idx_[i].count(j)<<std::endl;
					A_[cols_idx_[i].at(j)]+=value;

// for(std::size_t i=0;i<A_.size();i++)
//   std::cout<<A_[i]<<" ";
// std::cout<<std::endl;
				}
				else
				{
// std::cout<<"does not already exists"<<std::endl;
					A_[cont_]=value ;

// for(std::size_t i=0;i<A_.size();i++)
//   std::cout<<A_[i]<<" ";
// std::cout<<std::endl;
					cols_idx_[i].insert (std::pair<Integer,Integer>(j,cont_) );
					cont_++ ; 
				}

			}

			void set_zero_row(const Integer i)
			{
				const auto& map=cols_idx_[i]; 
				for (auto it=map.begin(); it!=map.end(); ++it)
					A_[it->second]=0;
			}



			void print()
			{
				std::cout<<"printing matrix"<<std::endl;
				for(Integer i=0;i<max_rows_;i++)
				{
					const auto& map=cols_idx_[i]; 
					for (auto it=map.begin(); it!=map.end(); ++it)
						std::cout<<"( " << i<<", "<< it->first << ", " << A_[it->second]<<" )" ;
					std::cout<<std::endl;
				}
			}

			void print_val()
			{
				Integer cont_cols=0;

				std::cout<<"printing matrix"<<std::endl;
				for(Integer i=0;i<max_rows_;i++)
				{   
					cont_cols=0;
					const auto& map=cols_idx_[i]; 
					for (auto it=map.begin(); it!=map.end(); ++it)
						{
						  if(it->first>cont_cols)
						  {
						  	for (int s=cont_cols;s<it->first;s++)
						  		std::cout<<0<<" ";
						  }
                          cont_cols=it->first+1;
					      std::cout << A_[it->second]<<" ";
					     }
				for (int s=cont_cols;s<max_rows_;s++)
					std::cout<<0<<" ";

					std::cout<<std::endl;
				}
			}

			auto multiply(const Integer i,std::vector<T>& b)
			{
				T tmp=0;
				const auto& map=cols_idx_[i]; 
// std::cout<<" multiply ";
				for (auto it=map.begin(); it!=map.end(); ++it)
				{
					tmp+=A_[it->second]*b[it->first];
// std::cout<<it->first<<"  "<<it->second<<"  A="<<A_[it->second]<<"   b="<<b[it->first]<<" "<<std::endl;
// std::cout<<tmp<<" "<<std::endl;;
				}

				return tmp;
			}

			void row_static_condensation(const Integer i,const std::vector<T>& b,const std::vector<bool>& constraint,Real & b_i)
			{
				const auto& map=cols_idx_[i]; 
// std::cout<<" multiply ";
				for (auto it=map.begin(); it!=map.end(); ++it)
				{
					if(constraint[it->first])
					{
						b_i-=A_[it->second]*b[it->first];
						A_[it->second]=0;
					}
// std::cout<<it->first<<"  "<<it->second<<"  A="<<A_[it->second]<<"   b="<<b[it->first]<<" "<<std::endl;
// std::cout<<tmp<<" "<<std::endl;;
				}
			}

			auto multiply(std::vector<T>& b)
			{

				std::vector<T> result(max_rows_);
// std::cout<<std::endl<<std::endl;
				for(Integer i=0;i<max_rows_;i++)
				{
					result[i]=multiply(i,b);
// std::cout<<std::endl<<result[i]<<std::endl;
				}
// std::cout<<std::endl<<std::endl;
				return result;
			}

			auto&  operator()(const Integer i, const Integer j)
			{
				return A_[cols_idx_[i].at(j)];
			}


			auto&  operator() (const Integer i, const Integer j)const
			{
				return A_[cols_idx_[i].at(j)];
			}


			void row_gauss_seidel(std::vector<T>& x,const std::vector<T>& b)const
			{

				for(std::size_t i=0; i<max_rows_;i++)
				{
					auto tmp=b[i];
					const auto& map=cols_idx_[i]; 
					for (auto it=map.begin(); it!=map.end(); ++it)
					{
						const auto& j=it->first;
						if(j!=i)
							tmp-=A_[cols_idx_[i].at(j)]*x[j];
// std::cout<<it->first<<"  "<<it->second<<"  A="<<A_[it->second]<<"   b="<<b[it->first]<<" "<<std::endl;
// std::cout<<tmp<<" "<<std::endl;;
					}
					x[i]=(tmp)/A_[cols_idx_[i].at(i)];
				}
			}





		};





	template<typename Elem,Integer ElementOrder, typename Operator,Integer FEFamily,Integer Order,Integer Continuity,Integer NComponents>
																	class Variable
{
public:  
// using Points=ElemPoints<Elem,Operator,FEFamily,Order>;
using Points=ElemGeometricPoints<Elem,ElementOrder>;
static constexpr Integer Npoints=Points::type::Dim;
using BaseFunctionSpace=BaseElementFunctionSpace<Elem,FEFamily,Order,Continuity, NComponents>;
using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
using SingleType   = typename SingleTypeShapeFunction<FunctionSpace,Operator>::SingleType;
static constexpr Integer Rows=SingleType::Rows;
static constexpr Integer Cols=SingleType::Cols;
static constexpr Integer single_shape_function_components_number=Rows*Cols;
static constexpr Integer solution_array_size=Npoints*Rows*Cols;
using TotType = typename SingleTypeShapeFunction<FunctionSpace,Operator>::TotType;
static constexpr Integer Ntot=FunctionSpaceDofsPerSubEntityElem<ElemFunctionSpace<Elem,BaseFunctionSpace>,Elem::ManifoldDim>::value;
static constexpr Integer Ndofs=Ntot/NComponents;
static constexpr auto 
reference_values{reference_shape_function_init<Elem,Operator,FEFamily,Order,SingleType,Ndofs>(Points::points)};
using Map=MapFromReference<Operator,Elem,FEFamily>;  
static constexpr auto NShapes=reference_values.size();

// qui devi passare i nodi dell'elemento
// l'ouput e' lungo elem:n_nods()*ncomponents

// EACH COMPONENT IS TREATED SEPARATELY, se vec is the array of dofs related to one component

void init(const FiniteElem<Elem>& FE,Array<Real,Ntot> vec,const Integer component)
{

map_.init(FE);
const auto& map=map_();
// std::cout<<"component="<<component<<std::endl;
// we compute u= sum_i u_i Map(phi_i)=Map(sum_i u_i phi_i) for each dof point on the element
Integer cont=0;
// std::cout<<"reference_values="<<std::endl;
// std::cout<<reference_values<<std::endl;
// std::cout<<"vec="<<std::endl;
// std::cout<<vec<<std::endl;
// loop on dof points and evaluate u_0 phi_0 in those points
// if i have RT0 valutato su sei punti...cosa devo fare?
// loop esterno sui punti e' ok, ma poi internamente dovrei loopare solo sule shape functions

for(Integer k=0;k<Npoints ;k++)
{

// std::cout<<"pre pre pre output_="<<output_<<std::endl;
for(Integer i=0;i<Rows;i++)
for(Integer j=0;j<Cols;j++)            
{
mat_tmp_[k](i,j)=vec[0+component]*reference_values[0][k](i,j);
// std::cout<<"k,i,j="<<k<<","<<i<<","<<j<<"="<<mat_tmp_[k](i,j)<<std::endl;
}

// for(std::size_t qp=1;qp<Npoints;qp++)
//    for(std::size_t i=0;i<Ndofs;i++)
//     for(std::size_t j=0;j<Ndofs;j++)
// std::cout<<"pre pre output_="<<output_<<std::endl;

//for(Integer qp=1;qp<NPoints;qp++) TODO FIX CHECK <---------------
for(Integer qp=1;qp<NShapes;qp++)
for(Integer i=0;i<Rows;i++)
	for(Integer j=0;j<Cols;j++)            
	{
// std::cout<<"(qp,i, j)=("<<qp<<", "<<i<<", "<<j<<")"<<std::endl;
// std::cout<<output_<<std::endl;
		mat_tmp_[k](i,j)+=vec[qp*NComponents+component]*reference_values[qp][k](i,j);
// std::cout<<output_<<std::endl;
	}


	mat_tmp_[k]=map*mat_tmp_[k];



// std::cout<<"map="<<map<<std::endl;
// for(std::size_t i=0;i<SingleType::Rows;i++)
// for(std::size_t j=0;j<SingleType::Cols;j++)
// std::cout<<"pre output_="<<output_<<std::endl;
	for(Integer i=0;i<SingleType::Rows;i++)
		for(Integer j=0;j<SingleType::Cols;j++)
		{
			output_[cont]=mat_tmp_[k](i,j);
// std::cout<<"output_["<<cont<<"]="<<output_[cont]<<std::endl;
			cont++;
		}
// std::cout<<"after output_="<<output_<<std::endl;

// }
	}

// std::cout<<"Npoints="<<Npoints<<std::endl;
// std::cout<<"NComponents="<<NComponents<<std::endl;
// std::cout<<"map="<<map<<std::endl;
// std::cout<<"component="<<component<<std::endl;
// std::cout<<"vec="<<vec<<std::endl;
// std::cout<<"mat_tmp_="<<mat_tmp_<<std::endl;
// std::cout<<" all output_="<<output_<<std::endl;
// std::cout<<"reference_values="<<reference_values<<std::endl;
}

auto& value()const{
	return output_;}

	auto& operator[](const Integer i)const{
// auto& operator[](const std::size_t i)const{
		return output_[i];}


	private:
		Map map_;
		Array<SingleType,Npoints*NComponents> mat_tmp_;
		Array<Real,solution_array_size> output_;
// std::vector<> 

	};


template<typename...Strings>
	auto variables_names(const std::string& string,const Strings&...strings)
	{

		return std::vector<std::string>{string,strings...};
	}

template<Integer Dim, Integer Order,typename Elem>
	auto vtk_cell_type(const Elem&elem)
	{
// linear triangle
		if(IsSame<Elem,Simplex<Dim,2>>::value && Order==1)
			return 5;
// linear tetrahedron
		if(IsSame<Elem,Simplex<Dim,3>>::value && Order==1)
			return 10;
// linear segment
		if(IsSame<Elem,Simplex<Dim,1>>::value && Order==2)
			return 21;
// quadratic triangle
		if(IsSame<Elem,Simplex<Dim,2>>::value && Order==2)
			return 22;
// quadratic tetrahedron
		if(IsSame<Elem,Simplex<Dim,3>>::value && Order==2)
			return 24;
// invalid element
		else
			return -1;
	}


template<Integer N, Integer Nmax, typename Space,Integer ElementOrder,typename FiniteElem,typename...Args,typename Solution,typename VariablesNamesVec>
	void print_solutions_aux_aux(std::ostream &os,FiniteElem& FE,const FullSpace<Args...>& space,const Solution&sol,const VariablesNamesVec& names,const Integer level=-1)
	{
		using Elem=typename Space::Elem;
		constexpr auto Dim=Elem::Dim;
		constexpr auto ManifoldDim=Elem::ManifoldDim;
		constexpr auto FEFamily=Space::FEFamily;
		constexpr auto Order=Space::Order;
		constexpr auto NComponents=Space::NComponents;
		constexpr auto Continuity=Space::Continuity;
		using Var=Variable<Elem,ElementOrder,IdentityOperator,FEFamily,Order,Continuity,NComponents>;
		using FunctionSpace=FullSpace<Args...>;
		using ElemDofMap=GetType<typename FunctionSpace::DofsDM::ElemDofMap,N>;
		constexpr auto Nelem_dofs=FunctionSpace::Nelem_dofs_array[N];
		constexpr auto single_shape_function_components_number=Var::single_shape_function_components_number;
		constexpr auto size=Var::solution_array_size;
		constexpr auto Npoints=Var::Npoints;
		Array<Real,Nelem_dofs> local_sol;
		auto m_max=size/Nelem_dofs-1;
		Var var;

		auto& bisection=space.bisection();
		auto& tracker=bisection.tracker();
// FE.init(0);
// var.init(FE,sol_tmp);
// std::cout<<"sol_tmp="<<std::endl;
// std::cout<<sol_tmp<<std::endl;

		std::vector<std::string> sub_scripts{"_x","_y","_z"};

		std::vector<std::string> sub_scripts_numbers{"_1","_2","_3","_4","_5","_6","_7","_8","_9","_10"};
// std::cout<<"<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>"<<std::endl;
// std::cout<<"Nelem_dofs = "<<Nelem_dofs<<std::endl;
// std::cout<<"size = "<<size<<std::endl;
// size the number of components of the finite element shape function (for 1 component)


		for(Integer comp=0;comp<NComponents;comp++)
		{
			if(single_shape_function_components_number==1)
				os << "SCALARS ";
			else if(single_shape_function_components_number>1)
				os << "VECTORS ";
			else
			{
				os << "ERROR ";  
			}

			if(NComponents==1)
				os << names[N];
			else
				os << names[N]+sub_scripts_numbers[comp];

			os << " float";
			os << "\n";
			if(single_shape_function_components_number==1)
				os << "LOOKUP_TABLE default\n";

// const auto& dofmap=tuple_get<N>(space.spaces_ptr()->dofmap2());
			const auto& dofsdofmap=space.dofsdofmap();
			ElemDofMap localdofmap;

			// auto mesh_ptr=space.mesh_ptr();
			auto& mesh=space.mesh();


// decltype(space.spaces_ptr()->dofmap2()) okl(6);
// std::cout<<"N===="<<N<<std::endl;
			Integer cont;
			// auto n_elements=mesh_ptr->n_elements();    
			auto n_elements=mesh.n_elements();    

			for(Integer i=0;i<n_elements;i++)
			{

				// if(!mesh_ptr->is_active(i)) continue;
				// if(!elem_belongs_to_level(mesh_ptr,i,level,tracker))continue;
				if(!elem_belongs_to_level(mesh,i,level,tracker))continue;


				FE.init(i,level);
				dofsdofmap.template dofmap_get<N>(localdofmap,i,level);
                // std::cout<<"elem id=="<<i<<std::endl;
				// std::cout<<"localdofmap=="<<localdofmap<<std::endl;
// const auto& localdofmap=dofmap[i];
				subarray(local_sol, sol, localdofmap);
				var.init(FE,local_sol,comp);
				const auto& val=var.value();
// std::cout<<"ElementOrder="<<ElementOrder<<std::endl;
// std::cout<<" val = "<<val<<std::endl;
// std::cout<<"localdofmap="<<localdofmap<<std::endl;
// std::cout<<"local_sol="<<local_sol<<std::endl;
// std::cout<<"var.value()="<<var.value()<<std::endl;
// std::cout<<"Nelem_dofs="<<Nelem_dofs<<std::endl;
// std::cout<<"var.value().size()="<<var.value().size()<<std::endl;
// std::cout<<"size="<<size<<std::endl;
				cont=0;

// vector value finite elements (Raviart-Thomas)
				if(single_shape_function_components_number>1)
				{
// std::cout<<"ci entro"<<std::endl;

					for(Integer j=0;j<Npoints;j++)
					{
																												for(Integer c=0;c<single_shape_function_components_number-1;c++)
																												{
	                os << val[cont];//sol[localdofmap[cont]];
	                // std::cout<<val[cont]<<" ";
	                os <<" "; 
	                cont++;
	            }
	            os << val[cont];
	            cont++;  
	            if(Dim==1)
	            {
	            	os <<" "; 
	            	os <<0.0;
	            }  
	            if(Dim==2)
	            {
	            	os <<" "; 
	            	os <<0.0;
	            } 
	            os <<"\n"; 
	             // std::cout<<std::endl;
	        } 

	    }
	         // scalar value finite elements (Lagrangian)
	    else
	    {
	    	for(Integer j=0;j<Npoints;j++)
	    	{
	          // for(Integer c=0;c<NComponents-1;c++)
	          //     {
	          //       std::cout<<" ____ "<<c<<std::endl; 
	          //       std::cout<< var[c]<<" "<<std::endl; 
	          //       os << var[c];//sol[localdofmap[cont]];
	          //       os <<" "; 
	          //       cont++;
	          //     }
	    		os<< val[cont];
	             // os << sol[localdofmap[cont]];
	    		cont++;  
	             // if(Dim==1 && NComponents==1)
	             // {
	             //  os <<" "; 
	             //  os <<0.0;
	             // }  
	             // if(Dim==2&& NComponents==2)
	             // {
	             //  os <<" "; 
	             //  os <<0.0;
	             // } 
	    		os <<"\n"; 
	    	} 
	    }

	          // os << "\n"; 
	}  
	if(N<Nmax)
		{os << "\n"; }
	}
	}


	template<Integer N, Integer Nmax,typename TupleOfSpaces,Integer ElementOrder,typename FiniteElem,typename...Args,typename Solution,typename VariablesNamesVec>
	constexpr std::enable_if_t< (N>Nmax), void>
	print_solutions_aux(std::ostream &os,FiniteElem& FE, const FullSpace<Args...>& space,const Solution&sol,const VariablesNamesVec& names,const Integer level=-1)
	{}


	template<Integer N, Integer Nmax, typename TupleOfSpaces,Integer ElementOrder,typename FiniteElem,typename...Args,typename Solution,typename VariablesNamesVec>
	constexpr std::enable_if_t< (N<=Nmax), void>
	print_solutions_aux(std::ostream &os,FiniteElem& FE, const FullSpace<Args...>& space,const Solution&sol,const VariablesNamesVec& names,const Integer level=-1)
	{
		using Space=GetType<TupleOfSpaces,N>;

		print_solutions_aux_aux<N,Nmax,Space,ElementOrder>(os,FE,space,sol,names,level);
		print_solutions_aux<N+1,Nmax,TupleOfSpaces,ElementOrder>(os,FE,space,sol,names,level);
	}

	template<typename TupleOfSpaces,Integer ElementOrder,typename...Args,typename Solution,typename VariablesNamesVec>
	constexpr void print_solutions(std::ostream &os,const FullSpace<Args...>& space,const Solution&sol,const VariablesNamesVec& names,const Integer level=-1)
	{
		using Elem=typename FullSpace<Args...>::Elem; 
		// FiniteElem<Elem> FE(space.mesh_ptr());
		FiniteElem<Elem> FE(space.mesh());
		print_solutions_aux<0,TupleTypeSize<TupleOfSpaces>::value-1,TupleOfSpaces,ElementOrder>(os,FE,space,sol,names,level);
	}



	template<Integer IsIsoparametric,typename TupleOfSpaces>
	class IsIsoparametricOrder
	{
	public:
		static constexpr Integer value=1+IsIsoparametric*(ElementOrder<TupleOfSpaces>::value-1);
	}; 

	template<Integer IsIsoparametric,typename...Args,typename Solution,typename VariablesNamesVec>
	void write_wtk_aux(std::ostream &os,const FullSpace<Args...>& space,const Solution&sol,const VariablesNamesVec& names,const Integer level=-1) {
		using FullSpace=FullSpace<Args...>;
		using TupleOfSpaces=typename FullSpace::FunctionSpace::TupleOfSpaces;
		using MeshT=typename FullSpace::MeshT;
		using Elem=typename MeshT::Elem;
		using ElemPoints=ElemGeometricPoints<Elem,IsIsoparametricOrder<IsIsoparametric,TupleOfSpaces>::value>;
		static constexpr auto Points=ElemPoints::points;
		// const auto mesh_ptr=space.mesh_ptr();
		const auto& mesh=space.mesh();
		const auto Dim=MeshT::Dim;
		const auto ManifoldDim=MeshT::ManifoldDim;
		// FiniteElem<Elem> FE(mesh_ptr);
		FiniteElem<Elem> FE(mesh);
		os << "# vtk DataFile Version 1.0\n";
		os<<Dim;
		os << "D Unstructured Grid of Linear Triangles\n";
		os << "ASCII\n";
		os << "\n";
		os << "DATASET UNSTRUCTURED_GRID\n";
		os << "POINTS ";
	// auto n_points=mesh_ptr->points().size();
		auto& bisection=space.bisection();
		auto& tracker=bisection.tracker();
		// auto n_elements=mesh_ptr->n_elements();
		auto n_elements=mesh.n_elements();

		Integer n_active_elements=0;
		for(Integer i=0;i<n_elements;i++)
		{
			// if(!mesh_ptr->is_active(i)) continue;
			// if(!elem_belongs_to_level(mesh_ptr,i,level,tracker))continue;
			if(!elem_belongs_to_level(mesh,i,level,tracker))continue;
			n_active_elements++;
		}


	auto n_points_per_elem=Points.size();//(mesh_ptr->elem(0).nodes.size());
	auto n_points=n_active_elements * n_points_per_elem;//(mesh_ptr->elem(0).nodes.size());
	// const auto& points=mesh_ptr->points();
	os << n_points;
	os << " float\n";

	// WRITE ELEMENTS
	for(Integer i=0;i<n_elements;i++)
	{
		// if(!mesh_ptr->is_active(i)) continue;
		// if(!elem_belongs_to_level(mesh_ptr,i,level,tracker))continue;
		// const auto& e=mesh_ptr->elem(i);
        if(!elem_belongs_to_level(mesh,i,level,tracker))continue;
		const auto& e=mesh.elem(i);
		FE.init(i);
		for(Integer j=0;j<Points.size();j++)
		{
			const auto point=FE.jac()*Points[j];
			const auto v0=FE.v0();
	     // std::cout<<"Points="<<Points[j]<<std::endl;
	     // std::cout<<"point="<<point<<std::endl;
	     // std::cout<<"v0="<<v0<<std::endl;

			for(Integer k=0;k<Dim-1;k++)
			{
				auto tmp=point(k,0)+v0[k];
				os << tmp ;
	      // os << point(k,0)+v0[k];
				os <<" ";
			}       
			os << point(Dim-1,0)+v0[Dim-1];

			if(Dim==1)
			{
				os <<" ";
				os << 0.0;
			}

			if(Dim==2)
			{
				os <<" ";
				os << 0.0;
			}

			os << "\n"; 
		}


	}


	os << "\n"; 
	os << "CELLS ";
	os << n_active_elements;
	os << " "; 
	os << n_active_elements * (1 + n_points_per_elem);//mesh_ptr->elem(0).nodes.size()); 
	os << "\n"; 

	// WRITE NODES AS IF WE HAVE DISCONTINUOUS FIELDS 
	// a common node between n elements result in n nodes
	Integer count=0;
	for(Integer i=0;i<n_elements;i++)
	{
		// if(!mesh_ptr->is_active(i)) continue;
		// if(!elem_belongs_to_level(mesh_ptr,i,level,tracker))continue;
		// const auto& e=mesh_ptr->elem(i);

		if(!elem_belongs_to_level(mesh,i,level,tracker))continue;
		const auto& e=mesh.elem(i);		
	    const auto size=n_points_per_elem;//e.nodes.size();

	    os << size;
	    os <<" ";   

	    for(Integer j=0;j<size-1;j++)
	    {
	    	os << count;
	    	os <<" ";
	    	count++;
	    } 
	    os << count;
	    os << "\n"; 
	    count++;
	}

	os << "\n"; 
	os << "CELL_TYPES ";
	os << n_active_elements;
	os << "\n";

	auto MaxElement=ElementOrder<TupleOfSpaces>::value;
	// std::cout<<"MaxElementOrder<TupleOfSpaces>>:value="<<MaxElement<<std::endl;
	// std::cout<<Points<<std::endl;
	for(Integer i=0;i<n_elements;i++)
	{
		// if(!mesh_ptr->is_active(i)) continue;
		// if(!elem_belongs_to_level(mesh_ptr,i,level,tracker))continue;
		// const auto& e=mesh_ptr->elem(i);

		if(!elem_belongs_to_level(mesh,i,level,tracker))continue;
		const auto& e=mesh.elem(i);

		os <<vtk_cell_type<Dim,IsIsoparametricOrder<IsIsoparametric,TupleOfSpaces>::value>(e);
		os << "\n"; 
	}
	os << "\n";

	// std::cout<<"names.size()--->"<<names.size()<<std::endl;
	os << "POINT_DATA ";
	os << n_points;
	os << "\n";

    std::cout<<"before print_solutions"<<std::endl;
	print_solutions<TupleOfSpaces,IsIsoparametricOrder<IsIsoparametric,TupleOfSpaces>::value>(os,space,sol,names,level);



	// std::ofstream file(baseFilename);
	// std::string my_string = "Hello text in file\n";
	// file << my_string;



	//   std::string gridFileName = baseFilename + ".vtk";
	//   std::ofstream outputfile;
	//   outputfile.open(gridFileName);

	// std::cin>>"# vtk DataFile Version 1.0 ">>std::endl;

	// outputfile.close();
	}



	template<typename...Args,typename Solution,typename VariablesNamesVec>
	void write_wtk(std::ostream &os,const FullSpace<Args...>& space,const Solution&sol,const VariablesNamesVec& names,const Integer level=-1) 
	{
		write_wtk_aux<0>(os,space,sol,names,level);
	}

	template<typename...Args,typename Solution,typename VariablesNamesVec>
	void write_wtk_isoparametric(std::ostream &os,const FullSpace<Args...>& space,const Solution&sol,const VariablesNamesVec& names,const Integer level=-1) 
	{
		write_wtk_aux<1>(os,space,sol,names,level);
	}

	template<typename...Args,typename Solution,typename VariablesNamesVec>
	void write_wtk_isoparametric(std::ostream &os,const std::shared_ptr<FullSpace<Args...>>& space,const Solution&sol,const VariablesNamesVec& names,const Integer level=-1) 
	{
		write_wtk_aux<1>(os,*space,sol,names,level);
	}


	template<typename Elem, Integer EntityDim>
	class ElemEntity;

	template<typename Elem, Integer EntityDimMax, Integer EntityDim>
	 class EntitiesOfMeshTupleTypeHelper;

	template<typename Elem, Integer EntityDimMax, Integer EntityDim>
	class EntitiesOfMeshTupleTypeHelper
	{
	public:
	     using rest = typename EntitiesOfMeshTupleTypeHelper<Elem,EntityDimMax,EntityDim+1>::type;
	     using ens  = ElemEntity<Elem,EntityDim>;
	     using type = decltype( std::tuple_cat( std::declval< std::tuple<ens> >(),std::declval< rest >() ) );
	};


	template<typename Elem, Integer EntityDimMax>
	class EntitiesOfMeshTupleTypeHelper<Elem, EntityDimMax,EntityDimMax>
	{
	public:
	     using ens  = ElemEntity<Elem,EntityDimMax>;
	     using type = typename std::tuple<ens>;
	};

	template<typename Elem>
	using EntitiesOfMeshTupleType= typename EntitiesOfMeshTupleTypeHelper<Elem,Elem::ManifoldDim,0>::type;








	template<typename MeshT_>
	class ElemEntityCollection
	{
	public:
		using MeshT=MeshT_;
		using BisectionT=Bisection<MeshT>;
		using Elem=typename MeshT::Elem;
		using EntitiesTuple=EntitiesOfMeshTupleType<Elem>;
		static constexpr Integer ManifoldDim=Elem::ManifoldDim;
		// ElemEntity<Elem,N>(mesh,node2elem);

		template<Integer N=0>
		inline std::enable_if_t<(N>ManifoldDim),void>
		entity_collection_init(const MeshT& mesh){}

		template<Integer N=0>
		inline std::enable_if_t<(N<=ManifoldDim),void>
		entity_collection_init(const MeshT& mesh)
		{
			tuple_get<N>(tuple_entities_).init(mesh);
			entity_collection_init<N+1>(mesh);
		}


	    ElemEntityCollection(const MeshT& mesh, const Integer level=-1):
	    node2elem_update_(false),
	    node2elem_(mesh,level),
	    mesh_ptr_(std::make_shared<MeshT>(mesh)),
	    level_(level)
	    {
	        if(level==-1)
	        	node2elem_update_=true;
	    	entity_collection_init(mesh);
	    	for(std::size_t i=0;i<is_updated_.size();i++)
	    		is_updated_[i]=false;
	    }



	    
	    template<Integer N>
	    void init()
	    {   
	    	
	    	if(!node2elem_update_)
	    	{ 
	    		std::cout<<"1 begin update N="<<N<<std::endl;
	    		node2elem_.init();
	    		node2elem_update_=true;
	    	}
	        if(!is_updated_[N]&&level_==-1)
	    	{
	    	std::cout<<"2 real update N="<<N<<std::endl;
	    	is_updated_[N]=true;
	    	auto& ens=tuple_get<N>(tuple_entities_);
	    	ens.init_elem_entity(*mesh_ptr_,node2elem_.val(),level_);
	        }
	    	else if(!is_updated_[N]&&level_!=-1)
	    	{
	    	std::cout<<"3 real update N="<<N<<std::endl;
	    	is_updated_[N]=true;
	    	auto& ens=tuple_get<N>(tuple_entities_);
	    	std::cout<<"ens="<<N<<std::endl;
	    	ens.add_bisection(bisection_ptr_);
	    	std::cout<<"added bisec="<<N<<std::endl;
	    	ens.init_elem_entity(*mesh_ptr_,node2elem_.val(),level_);
	    	std::cout<<"end init elem="<<N<<std::endl;
	        }
	        else
	        {
	        	std::cout<<" already updated "<<std::endl;
	        }
	    }



	    template<typename FunctionSpace,Integer Nmax,Integer N>
	    std::enable_if_t<(N>Nmax),void>
	    init_aux_aux()
	    {}

	    template<typename FunctionSpace,Integer Nmax,Integer N>
	    std::enable_if_t<(N<=Nmax),void>
	    init_aux_aux()
	    {   
	    	init<FunctionSpace::entity[N]>();
	    	std::cout<<"init_aux_aux "<<N<<std::endl;
	    	init_aux_aux<FunctionSpace,Nmax,N+1>();
	    }


	    template<typename TupleOfSpaces,Integer Nmax,Integer N>
	    std::enable_if_t<(N>Nmax),void>
	    init_aux(){}

	    template<typename TupleOfSpaces,Integer Nmax,Integer N>
	    std::enable_if_t<(N<=Nmax),void>
	    init_aux()
	    {   
	    	using FunctionSpace=GetType<TupleOfSpaces,N>;
	    	std::cout<<"init_aux "<<N<<std::endl;
	    	init_aux_aux<FunctionSpace,FunctionSpace::entity.size()-1,0>();
	    	init_aux<TupleOfSpaces,Nmax,N+1>();
	    }


	    template<template<class...>class TemplateClass,typename...Args>
	    void init(const TemplateClass<Args...>& W)
	    {   
	    	using TupleOfSpaces=typename TemplateClass<Args...>::TupleOfSpaces;
	    	std::cout<<"init W"<<std::endl;
	    	init_aux<TupleOfSpaces,TupleTypeSize<TupleOfSpaces>::value-1,0>();
	    }



	    template<typename Elem,typename BaseFunctionSpace, typename...BaseFunctionSpaces>
	    std::enable_if_t<(sizeof...(BaseFunctionSpaces)==0),void>
	    init_aux()
	    {
	    	using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
	    	init_aux_aux<FunctionSpace,FunctionSpace::entity.size()-1,0>();
	    }

	    template<typename Elem,typename BaseFunctionSpace, typename...BaseFunctionSpaces>
	    std::enable_if_t<(sizeof...(BaseFunctionSpaces)>0),void>
	    init_aux()
	    {   
	    	using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
	    	init_aux_aux<FunctionSpace,FunctionSpace::entity.size()-1,0>();
	    	init_aux<Elem,BaseFunctionSpaces...>();
	    }


	    template<typename Elem,typename BaseFunctionSpace, typename...BaseFunctionSpaces>
	    void init()
	    {   
	    	std::cout<<"init "<<std::endl;
	    	init_aux<Elem,BaseFunctionSpace,BaseFunctionSpaces...>();
	    }

	     
	    auto& tuple_entities(){return tuple_entities_;}

		inline void add_bisection(const BisectionT& bisection)
		{node2elem_.add_bisection(bisection);
		 bisection_ptr_=std::make_shared<BisectionT>(bisection);}
	    
	    auto& node2elem(){return node2elem_.val();}

	    inline auto bisection_ptr(){return bisection_ptr_;}
	private:
		bool node2elem_update_;
		NodeToElem<MeshT> node2elem_;
		EntitiesTuple tuple_entities_;
		std::shared_ptr<MeshT> mesh_ptr_;
		std::shared_ptr<BisectionT> bisection_ptr_;
		Array<bool, ManifoldDim+1> is_updated_;
	    Integer level_;
	};

	template<Integer N,typename MeshT>
	void init(ElemEntityCollection<MeshT>& coll){coll.template init<N>();}

	template<Integer N,typename MeshT>
	auto& get_entity(ElemEntityCollection<MeshT>& coll){return tuple_get<N>(coll.tuple_entities());}





	template<typename MeshT_>
	class MeshAndEntity
	{
	 public:
	 	using MeshT=MeshT_;
	 	using Elem=typename MeshT::Elem;
	 	MeshAndEntity(const MeshT& mesh, const Integer level=-1):
	    mesh_ptr_(std::make_shared<MeshT>(mesh)),
	    level_(level),
	    entities_(mesh,level)
	 	{}


	    inline auto mesh_ptr(){return mesh_ptr_;}

	    inline auto level(){return level_;}

	    inline auto bisection_ptr() {return entities_.bisection_ptr();}

	    inline auto& entities_collection(){return entities_;}

	    template<template<class...>class TemplateClass,typename...Args>
	    inline void init(const TemplateClass<Args...>& W){return entities_.init(W);}

	    template<typename Elem,typename FunctionSpace, typename...FunctionSpaces>
	    inline void init(){return entities_.template init<Elem,FunctionSpace,FunctionSpaces...> ();}

		inline void add_bisection(const Bisection<MeshT>& bisection) {entities_.add_bisection(bisection);}

		


	 private:
	 	std::shared_ptr<MeshT> mesh_ptr_;
	 	Integer level_;
	 	ElemEntityCollection<MeshT> entities_;
	 	
	};



	template<Integer N,typename MeshT>
	auto& get_entity(MeshAndEntity<MeshT>& mesh_and_entity){return tuple_get<N>(mesh_and_entity.entities_collection().tuple_entities());}


	template<typename Elem, typename...FunctionSpaces,typename MeshT>
	void init(MeshAndEntity<MeshT>& mesh_and_entity){mesh_and_entity.template init<Elem,FunctionSpaces...>();}


















	template<typename MeshT,Integer EntityDim>
	class EntitySimplicialMap
	{
	 public:
	 using IntegerVector=std::vector<Integer>;
	 using Key=std::array<Integer,EntityDim+1>;
	 using Value=std::shared_ptr<IntegerVector>;
	 using Map=std::map<Key, Value>;
	 using MapDof=std::map<Key, std::shared_ptr<Integer>>;
	 static constexpr Integer Dim=MeshT::Dim;
	 static constexpr Integer ManifoldDim=MeshT::ManifoldDim;
	 static constexpr auto entity_combinations= ElemEntity<Simplex<Dim,ManifoldDim>,EntityDim>::entity_combinations();
	 EntitySimplicialMap(const MeshT& mesh):
	 mesh_ptr_(std::make_shared<MeshT>(mesh))
	 {}

	 void init()
	 {
	  Integer entity[EntityDim+1];
	  Integer cont=0;
	  for(std::size_t i=0;i<mesh_ptr_->n_elements();i++)
		{
			const auto& nodes=mesh_ptr_->elem(i).nodes;

	          for(Integer iter_entity=0;iter_entity<entity_combinations;iter_entity++)
	          {

	            Combinations<ManifoldDim + 1, EntityDim+1>::generate(iter_entity,entity);

	            for(int nn=0;nn<EntityDim+1;nn++)
	              local_nodes_[nn]=nodes[entity[nn]];   

	            std::sort(std::begin(local_nodes_), std::end(local_nodes_)); 


	            auto& val=map_[local_nodes_];
					if(!val) {
						val = std::make_shared<IntegerVector>();
					}
	            auto& val2=entity_[local_nodes_];
					if(!val2)
					{

						val2 = std::make_shared<Integer>(cont);
						// entity_.insert(std::pair<Key,Integer>(local_nodes_,cont));
						// entity_[local_nodes_]=cont;
						// std::cout<<" local_nodes_="<<std::endl;
						// for(int j=0;j<local_nodes_.size();j++)
						// std::cout<<local_nodes_[j]<<" ";
					 //    std::cout<<" CONT="<< entity_[local_nodes_]<<std::endl;

						cont++;
					}
	            val->push_back(i);
	            // .insert( std::pair<Key,Value>(local_nodes_,i)); 
	        // for(std::size_t j=0;j<local_nodes_.size();j++)
	        	// std::cout<<local_nodes_[j]<<" ";
	        // std::cout<<std::endl;
	        }


		}

	 }
	 		void describe(std::ostream &os) const
			{
				os << "-----------------------------------\n";

				for(const auto &m : map_) {
					const auto& nodes=m.first;
					for(const auto n : nodes) {
						os << n << " ";
					}

					os << "-> ";

					for(const auto n : *m.second) {
						os << n << " ";
					}

					os << "\n";
				}

				os << "-----------------------------------\n";
			}
	 		void describe2(std::ostream &os) const
			{
				os << "----------------------------------- entity nums -------------------------\n";

				for(const auto &m : entity_) {
					const auto& nodes=m.first;
					os << "NODES= ";
					for(const auto n : nodes) {
						os << n << " ";
					}

					os << "-> ";

						os << *m.second << " ";
					

					os << "\n";
				}

				os << "-----------------------------------\n";
			}

	 private:
	 	std::shared_ptr<MeshT> mesh_ptr_;
	 	Map map_;
	 	std::array<Integer, EntityDim+1> local_nodes_;
	    MapDof entity_;

	  
	};




	template<typename MeshT_>
	class Node2ElemMap
	{

	public:
	 using MeshT=MeshT_;
	 using BisectionT=Bisection<MeshT>;
	 static constexpr Integer ManifoldDim=MeshT::ManifoldDim;
	 static constexpr Integer Dim=MeshT::Dim;
	 using IntegerVector=std::vector<Integer>;
	 using Key=Integer;
	 using Value=std::shared_ptr<IntegerVector>;
	 using Map=std::map<Key, Value>;


	 
	 Node2ElemMap(MeshT& mesh):
	 mesh_(mesh),
	 last_n_elements_(0)
	 {}

	 inline void add_bisection(const Bisection<MeshT>& bisection) {bisection_ptr_=std::make_shared<BisectionT>(bisection);}

	 auto get(const Integer dof,const Integer level=0 )
	 {
	 	auto& val=map_[dof];
	 	if(!val) {
	 		init();
	 	}
	 	IntegerVector out;
	 	out.reserve(val->size());
	 	for(Integer i=0;i<val->size();i++)
	 	{
	 		if(elem_belongs_to_level(mesh_,(*val)[i],level,bisection_ptr_))
	 			out.push_back((*val)[i]);
	 	}
	 	return out;
	 }

	 void init()
	 {
	 	for(std::size_t i=last_n_elements_;i<mesh_.n_elements();i++)
	 	{
	        const auto& elem=mesh_.elem(i);
	     	last_n_elements_++;
		    for(Integer j=0;j<ManifoldDim+1;j++)
		          {
		            auto& val=map_[elem.nodes[j]];
						if(!val) {
							val = std::make_shared<IntegerVector>();
						}
		            val->push_back(i);
		        } 	  
		}                 
	 }
	 
	 auto& operator()()     {return map_;}
	 auto& operator()()const{return map_;}

	private:
		MeshT& mesh_;
		std::shared_ptr<BisectionT> bisection_ptr_;
	    Map map_;
	    Integer last_n_elements_;
	};



	template<typename T>
	std::vector<T> vector_intersection(std::vector<T> &v1, std::vector<T> &v2){
	    std::vector<T> v3;

	    std::sort(v1.begin(), v1.end());
	    std::sort(v2.begin(), v2.end());

	    std::set_intersection(v1.begin(),v1.end(),
	                          v2.begin(),v2.end(),
	                          back_inserter(v3));
	    return v3;
	}

	template<typename T>
	std::vector<T> vector_intersection(std::vector<T> &v1, std::vector<T> &&v2){
	    std::vector<T> v3;

	    std::sort(v1.begin(), v1.end());
	    std::sort(v2.begin(), v2.end());

	    std::set_intersection(v1.begin(),v1.end(),
	                          v2.begin(),v2.end(),
	                          back_inserter(v3));
	    return v3;
	}

	template<typename T,typename...Ts>
	std::vector<T> vector_intersection(std::vector<T> &v1, std::vector<T> &v2, std::vector<Ts> &...vs)
	{
	   
	    return vector_intersection(vector_intersection(v1,v2),vs...);
	}



	// a simplex is containd into another simplex if all its nodes are inside it
	template<Integer Dim,Integer ManifoldDim1,Integer ManifoldDim2,Integer ManifoldDim3>
	bool is_simplex_inside_simplex(const Simplex<Dim,ManifoldDim1>& simplex1,
								   const Simplex<Dim,ManifoldDim2>& simplex2,
								   const Mesh<Dim,ManifoldDim3>& mesh)
	{
	    bool is_inside=true;
	    const auto nodes1=simplex1.nodes;
		const auto nodes2=simplex2.nodes;
		Simplex<Dim,ManifoldDim1> simplex_tmp;
		Integer nodes_tmp[ManifoldDim1];

	    // loop on the nodes of simplex2
		for(std::size_t i=0;i<nodes2.size();i++)
		{

	      const auto& node=nodes2[i];

	      std::cout<<"node == "<< node<<std::endl;

	      // for a given node, check if it inside
	      // loop on all the face of simplex1 and create a temporary simplex
	      for(std::size_t j=0;j<ManifoldDim1+1;j++)
	      {
	       // std::cout<<"nodes j == "<< j<<std::endl;
	       Combinations<ManifoldDim1 + 1,ManifoldDim1>::generate(j,nodes_tmp);
	       // std::cout<<"ManifoldDim1-j == "<< node<<std::endl;
	       simplex_tmp.nodes[ManifoldDim1-j]= node;
		   for(std::size_t k=0;k<ManifoldDim1;k++)
		      {
		      	simplex_tmp.nodes[nodes_tmp[k]]= nodes1[nodes_tmp[k]];
		      }

	       std::cout <<std::endl;
	       std::cout<<"simplex_tmp =  " <<std::endl;
		   for(std::size_t k=0;k<ManifoldDim1+1;k++)
		      {
		      	std::cout<<simplex_tmp.nodes[k]<< " ";
		      } 
		    std::cout <<std::endl;
		   const auto vol=volume(simplex_tmp,mesh.points());
		   std::cout<<"vol =  "<< vol <<std::endl;
		   if(vol<=-0.00000000000001)
		      {is_inside=false;
		       break;}
		      	
		      // }  
		      // if(is_inside==false)
		      //    break;           

	      }
		   if(is_inside==false)
		      break;  
		}

		return is_inside;
	}
















	template<typename MeshT_,Integer EntityDim>
	class ConnectivitySimpliacialMap;

	template<typename MeshT_,Integer EntityDim>
	class ConnectivitySimpliacialMap
	{

	public:
	 using MeshT=MeshT_;
	 using BisectionT=Bisection<MeshT>;
	 using Elem=typename MeshT::Elem;
	 static constexpr Integer ManifoldDim=MeshT::ManifoldDim;
	 static constexpr Integer Dim=MeshT::Dim;
	 static constexpr auto entity_combinations= ElemEntity<Elem,EntityDim>::entity_combinations();
	 static constexpr auto Npoints=EntityDim+1;
	 using IntegerVector=std::vector<Integer>;
	 using Key=std::array<Integer,Npoints>;
	 using Value=std::shared_ptr<IntegerVector>;
	 using Map=std::map<Key, Value>;
	 using MapDof=std::map<Key,Integer>;

	 
	 ConnectivitySimpliacialMap(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection):
	 mesh_(mesh),
	 node2elem_(node2elem),
	 bisection_(bisection),
	 last_n_elements_(0),
	 count_n_entity_(0),
	 init_coarse_level_(false) 
	 {
	 	count_n_entity_vec_.resize(1,0);
	 }




	 void init_elem(const Elem elem,const Integer level )
	 {

	        // std::cout<<"___________LEVEL="<<level<<std::endl;
	        // std::cout<<"count_n_entity_vec_.size()-1="<<count_n_entity_vec_.size()-1<<std::endl;
	 		if(count_n_entity_vec_.size()-1<level)
	 			count_n_entity_vec_.push_back(count_n_entity_vec_[count_n_entity_vec_.size()-1]);
	 		// std::cout<<"___________count_n_entity_vec_[END]="<<std::endl;

	 		// for(std::size_t i=0;i<count_n_entity_vec_.size();i++)
	 			// std::cout<<count_n_entity_vec_[i]<<" ";
	 		// std::cout<<std::endl;
	 		Integer cont=0;
	 		Integer cont1=0;
	 		Integer cont2=0;

	 	 	const auto& parent_id=elem.parent_id;
		 	const auto& parent_elem=mesh_.elem(parent_id);
		 	const auto& parent_nodes=parent_elem.nodes;
		 	std::array<Integer,Npoints> parent_entity_nodes;
		 	std::array<Integer,Npoints> child_entity_nodes;
		 	std::array<Integer,Npoints> entity_used;
	        bool found_parent_entity;
	 	 
		 	 // find the parent element
		 	 // loop on all its entities parent_entity of dimension EntityDim
		 	 // given entity parent_entity:
		 	 //     loop on all the child elements child_elem
		 	 //          loop on all the entities child_entity
		 	 //               if child_entity==parent_entity
	         //                  do nothing, such entity already exists
	         //      if no agreement has been found,  cont++
	         // then we start creating new dofs for the entities on this level
	         // entities which belong to both coarse and fine level are untouched
	         // the new ones are created, counting them from n_coarse_entities-cont
		 	 // we loop on all the entities of dimension EntityDim of the element
		 	 const auto& child=parent_elem.children;

	         for(std::size_t parent_entity=0;parent_entity<entity_combinations;parent_entity++)
	         	{
				 Combinations<ManifoldDim + 1, EntityDim+1>::generate(parent_entity,entity_);
	             // std::cout<<std::endl;
	         	 for(std::size_t i=0;i<Npoints;i++)
	                	{
	                		parent_entity_nodes[i]=parent_nodes[entity_[i]];
	                		// std::cout<<parent_entity_nodes[i]<<" ";
	                	}
	             // std::cout<<std::endl;
	             std::sort(parent_entity_nodes.begin(),parent_entity_nodes.end());
	            

	            if(!parent_map_[parent_entity_nodes])
	            {
	            	// std::cout<<" entro in !parent_map_"<<std::endl;
	            found_parent_entity=false;
	            // loop on all the child elements child_elem          
			 	for(std::size_t i=0;i<child.size();i++)
			 	{
			 	 const auto& child_elem=mesh_.elem(child[i]);
			 	 const auto& child_id=child_elem.id;
			 	 const auto& child_nodes=child_elem.nodes;
	             // loop on all the entities child_entity
		         for(std::size_t child_entity=0;child_entity<entity_combinations;child_entity++)
		         	{
					 Combinations<ManifoldDim + 1, EntityDim+1>::generate(child_entity,entity_);

		         	 for(std::size_t i=0;i<Npoints;i++)
		                	child_entity_nodes[i]=child_nodes[entity_[i]];
		             std::sort(child_entity_nodes.begin(),child_entity_nodes.end());
	                 // parent_map_[parent_entity_nodes]=true;
		             found_parent_entity=std::equal(child_entity_nodes.begin(),child_entity_nodes.end(),parent_entity_nodes.begin());
				     
				     if(found_parent_entity)			     	 
				     	 goto label;
				     	
				     }
				  }
				  label:
				  if(found_parent_entity)
				  {
				  	// parent_map_[parent_entity_nodes]=true;
				  	// std::cout<<"found_parent_entity="<<found_parent_entity<<std::endl;
				  	// entity_used[cont1]=map_dof_[parent_entity_nodes];
				  	// cont1++;
				  }
				  else
				  {
				     	// std::cout<<"entity_used["<<cont1<<"]="<<map_dof_[parent_entity_nodes]<<std::endl;
	                    parent_map_[parent_entity_nodes]=true;
				  	    entity_used[cont]=map_dof_[parent_entity_nodes];
				  		cont++;
				  }
				  {}
				  }
		         }


	            count_n_entity_vec_[level]=count_n_entity_vec_[level];//-cont;
	            // loop on all the child elements child_elem          
			 	for(std::size_t i=0;i<child.size();i++)
			 	{
			 	 const auto& child_elem=mesh_.elem(child[i]);
			 	 const auto& child_id=child_elem.id;
			 	 const auto& child_nodes=child_elem.nodes;			 		
			 	 if(elem2dofs_[child_id].size()==0)
			 	 {

	             // loop on all the entities child_entity
		         for(std::size_t child_entity=0;child_entity<entity_combinations;child_entity++)
		         	{
					 Combinations<ManifoldDim + 1, EntityDim+1>::generate(child_entity,entity_);
		         	 for(std::size_t i=0;i<Npoints;i++)
		                	child_entity_nodes[i]=child_nodes[entity_[i]];
		             std::sort(child_entity_nodes.begin(),child_entity_nodes.end());
	                 
	                 auto& child_map=map_[child_entity_nodes];

	                 // if the entity has not been visited yet
	                 if(!child_map)
	                 {
	                 	// create the vector
	                    child_map=std::make_shared<IntegerVector>();
	                    auto& child_dof=map_dof_[child_entity_nodes];
	                    if(cont2<cont)
	                    {
	                     child_dof=entity_used[cont2];
	                     // std::cout<<" cont2<cont "<<child_dof<<std::endl;
	                     cont2++;	
	                    }
	                    else
	                    {child_dof=count_n_entity_vec_[level];
	                     // std::cout<<" cont2>=cont "<<child_dof<<std::endl;

						 count_n_entity_vec_[level]++;}
	                 }
	                 // std::cout<<"count_n_entity_vec_[level]="<<count_n_entity_vec_[level]<<std::endl;
	                 child_map->push_back(child_id);
	                 std::cout<<" child_id = "<<child_id<<" with "<<map_dof_[child_entity_nodes]<<std::endl;
	                 elem2dofs_[child_id].push_back(map_dof_[child_entity_nodes]);
				  }
			 	 }
				}
		 	 }


	 void init_coarse()
	 {  
	 	node2elem_.init();
	 	elem2dofs_.resize(mesh_.n_elements());
	 	for(std::size_t i=last_n_elements_;i<mesh_.n_elements();i++)
	 	{
	 		const auto& elem=mesh_.elem(i);
	 		if(elem.parent_id==INVALID_INDEX)
	 		{
	 		last_n_elements_++;
			const auto& nodes=mesh_.elem(i).nodes;

	          for(Integer iter_entity=0;iter_entity<entity_combinations;iter_entity++)
	          {
	           // std::cout<<" ini coarse = "<<i<<std::endl;
	            Combinations<ManifoldDim + 1, Npoints>::generate(iter_entity,entity_);

	            for(int nn=0;nn<Npoints;nn++)
	              local_nodes_[nn]=nodes[entity_[nn]];   

	            std::sort(std::begin(local_nodes_), std::end(local_nodes_)); 


	            auto& val=map_[local_nodes_];
					if(!val) {
						val = std::make_shared<IntegerVector>();
						// parent_map_[local_nodes_]=true;
					}
	            auto& val2=map_dof_[local_nodes_];
					if(!val2)
					{

						val2 = count_n_entity_;
						count_n_entity_++;
					}
	            val->push_back(i);
	            // std::cout<<"qui elem2dofs_.size="<<elem2dofs_.size()<<std::endl;
	             // std::cout<<"elem.parent_id="<<elem.id<<std::endl;
                elem2dofs_[elem.id].push_back(val2);
                // std::cout<<" after qui elem2dofs_.size="<<elem2dofs_.size()<<std::endl;
	          }
	 		}
	      }
	     // std::cout<<"end1 ini coarse = "<<std::endl;
	     count_n_entity_vec_[0]=count_n_entity_;
	     // std::cout<<" end2 coarse = "<<std::endl;
	     // count_levels_=1;
	     // std::cout<<" end3 coarse = "<<std::endl;
		}             

	 void init()
	 {
	 	node2elem_.init();
	 	elem2dofs_.reserve(mesh_.n_elements());

        // we initialize the coarse level only the first time
	 	if(!init_coarse_level_)
	 	{
	 		init_coarse();
	 		init_coarse_level_=true;
	 	}
        
        // std::cout<<"elem2dofs_ with size="<<elem2dofs_.size()<<std::endl;
        // for(std::size_t i=0;i<elem2dofs_.size();i++)
        // {
        // 	std::cout<<std::endl;
        // 	for(std::size_t j=0;j<elem2dofs_[i].size();j++)
        // 	std::cout<<elem2dofs_[i][j]<<" ";
        // std::cout<<std::endl;
        // }



	    // we assume the finest elements are added after the coarsest ones  
	 	const auto& tracker=bisection_.tracker();
	    Integer level;
	    // std::cout<<"last_n_elements_="<<last_n_elements_<<std::endl;
	    // std::cout<<"mesh_.n_elements()="<<mesh_.n_elements()<<std::endl;    
	 	for(std::size_t i=last_n_elements_;i<mesh_.n_elements();i++)
	 	   {
	 	   	
	 	   	level=tracker.get_iterate(i);
	 	   	// std::cout<<"init elem="<<i <<" of level "<<level<<std::endl;
	 		init_elem(mesh_.elem(i),tracker.get_iterate(i));
	 		// std::cout<<count_n_entity_vec_[level]<<std::endl;
	 		last_n_elements_++;	
	        }
	    last_n_elements_=mesh_.n_elements();

	 	

	  } 


	 		void describe(std::ostream &os) const
			{
				os << "-----------------------------------\n";

				for(const auto &m : map_) {
					const auto& nodes=m.first;
					for(const auto n : nodes) {
						os << n << " ";
					}

					os << "-> ";

					for(const auto n : *m.second) {
						os << n << " ";
					}

					os << "\n";
				}

				os << "-----------------------------------\n";
			}

	 		void describe2(std::ostream &os) const
			{
				os << "-----------------------------------\n";

				for(const auto &m : map_dof_) {
					const auto& nodes=m.first;
					for(const auto n : nodes) {
						os << n << " ";
					}

					os << "-> ";

						os << m.second << " ";
					

					os << "\n";
				}

				os << "-----------------------------------\n";
			}		



	 		void describe3(std::ostream &os) const
			{
				os << "-----------------------------------\n";

				for(const auto &m : parent_map_) {
					const auto& nodes=m.first;
					for(const auto n : nodes) {
						os << n << " ";
					}

					os << "-> ";

						os << m.second << " ";
					

					os << "\n";
				}

				os << "-----------------------------------\n";
			}	



			auto& elem2dofs(const Integer i)const{return elem2dofs_[i];}


	private:
		MeshT& mesh_;
		Node2ElemMap<MeshT>& node2elem_;
		Bisection<MeshT>& bisection_;
		std::shared_ptr<BisectionT> bisection_ptr_;
	    Map map_;
	    MapDof map_dof_;
	    std::map<Key, bool> parent_map_;
	    Key local_nodes_;
	    Integer entity_[Npoints];
	    Simplex<Dim,EntityDim> simplex_;
	    Simplex<Dim,EntityDim> simplex_parent_;
	    Integer last_n_elements_;
	    Integer count_n_entity_;
	    IntegerVector count_n_entity_vec_;
	    // Integer count_levels_;
	    bool init_coarse_level_;
	    std::vector<std::vector<Integer>> elem2dofs_;
	};



















	template<typename Elem, Integer EntityDimMax, Integer EntityDim>
	 class ConnectivitySimpliacialMapOfMeshTupleTypeHelper;

	template<typename MeshT, Integer EntityDimMax, Integer EntityDim>
	class ConnectivitySimpliacialMapOfMeshTupleTypeHelper
	{
	public:
	     using rest = typename ConnectivitySimpliacialMapOfMeshTupleTypeHelper<MeshT,EntityDimMax,EntityDim+1>::type;
	     using ens  = ConnectivitySimpliacialMap<MeshT,EntityDim>;
	     using type = decltype( std::tuple_cat( std::declval< std::tuple<ens> >(),std::declval< rest >() ) );
	};


	template<typename MeshT, Integer EntityDimMax>
	class ConnectivitySimpliacialMapOfMeshTupleTypeHelper<MeshT, EntityDimMax,EntityDimMax>
	{
	public:
	     using ens  = ConnectivitySimpliacialMap<MeshT,EntityDimMax>;
	     using type = typename std::tuple<ens>;
	};

	template<typename MeshT>
	using ConnectivitySimpliacialMapOfMesh= typename ConnectivitySimpliacialMapOfMeshTupleTypeHelper<MeshT,MeshT::Elem::ManifoldDim,0>::type;



			template<typename MeshT_>
	class ConnectivitySimpliacialMapCollection
	{
	public:
		using MeshT=MeshT_;
		using BisectionT=Bisection<MeshT>;
		using Elem=typename MeshT::Elem;
		using EntitiesTuple=ConnectivitySimpliacialMapOfMesh<MeshT>;
		static constexpr Integer ManifoldDim=Elem::ManifoldDim;


			template<Integer N=0>
		struct ConstructorTupleHelper;

			template<>
		struct ConstructorTupleHelper<0>
		{
			using ens = ConnectivitySimpliacialMap<MeshT,0>;
			using type = typename std::tuple<ens>;
		};

			template <Integer N>
		struct ConstructorTupleHelper
		{
			using rest = typename ConstructorTupleHelper<N-1>::type; 
			using ens = ConnectivitySimpliacialMap<MeshT,N>;
			using tuple_ens=std::tuple<ens>;
			using type = decltype( std::tuple_cat(std::declval< rest >() , std::declval< tuple_ens >()) );
		};

			template<Integer N>
		using ConstructorTuple=typename ConstructorTupleHelper<N>::type;




			template< Integer N>
		std::enable_if_t<0==N, ConstructorTuple<N> >
		construct_tuple(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection)
		{   
			using type=ConnectivitySimpliacialMap<MeshT,N>;
			return std::tuple<type>
			(type(mesh,node2elem,bisection));
		}


			template<Integer N>
		std::enable_if_t< 0<N, ConstructorTuple<N> >
		construct_tuple(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection)
		{
			using type=ConnectivitySimpliacialMap<MeshT,N>;
			return std::tuple_cat(construct_tuple<N-1>(mesh,node2elem,bisection), 
				std::tuple<type>(type(mesh,node2elem,bisection)));
		}








		ConnectivitySimpliacialMapCollection(MeshT& mesh,Node2ElemMap<MeshT>& node2elem,Bisection<MeshT>& bisection):
		mesh_(mesh),
		node2elem_(node2elem),
		bisection_(bisection),
		tuple_entities_(construct_tuple<MeshT::ManifoldDim>(mesh,node2elem,bisection))
		{}



			    template<typename FunctionSpace,Integer Nmax,Integer N>
		std::enable_if_t<(N>Nmax),void>
		init_aux_aux()
		{}

			    template<typename FunctionSpace,Integer Nmax,Integer N>
		std::enable_if_t<(N<=Nmax),void>
		init_aux_aux()
		{   
			auto& ens=tuple_get<FunctionSpace::entity[N]>(tuple_entities_);
			ens.init();
			init_aux_aux<FunctionSpace,Nmax,N+1>();
		}


			    template<typename BaseFunctionSpace, typename...BaseFunctionSpaces>
		std::enable_if_t<(sizeof...(BaseFunctionSpaces)==0),void>
		init_aux()
		{
			using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
			init_aux_aux<FunctionSpace,FunctionSpace::entity.size()-1,0>();
		}

			    template<typename BaseFunctionSpace, typename...BaseFunctionSpaces>
		std::enable_if_t<(sizeof...(BaseFunctionSpaces)>0),void>
		init_aux()
		{   
			using FunctionSpace=ElemFunctionSpace<Elem,BaseFunctionSpace>;
			init_aux_aux<FunctionSpace,FunctionSpace::entity.size()-1,0>();
			init_aux<BaseFunctionSpaces...>();
		}


			    template<typename BaseFunctionSpace, typename...BaseFunctionSpaces>
		void init()
		{   
			init_aux<BaseFunctionSpace,BaseFunctionSpaces...>();
		}


        template<Integer Nmax,Integer N>
	 	std::enable_if_t<(N>Nmax),void> describe2_aux(std::ostream &os) const
	 	{}

        template<Integer Nmax,Integer N>
	 	std::enable_if_t<N<=Nmax,void> describe2_aux(std::ostream &os) const
		{
			os << "-----------------ENTITY N=="<<N<<"\n";
			tuple_get<N>(tuple_entities_).describe2(os);
			describe2_aux<Nmax,N+1>(os);
		}

	 	void describe2(std::ostream &os) const
		{
         describe2_aux<ManifoldDim,0>(os);
		}	


		auto& tuple_entities(){return tuple_entities_;}

		auto& node2elem(){return node2elem_;}

		auto& mesh(){return mesh_;}

		auto& bisection(){return bisection_;}

        const auto n_levels(){return bisection_.tracker().current_iterate();}
		// auto& mesh_ptr(){return std::make_shared<MeshT>(mesh_);}

		// inline auto bisection_ptr(){return bisection_ptr_;}
	private:
		MeshT& mesh_;
		Node2ElemMap<MeshT>& node2elem_;
		Bisection<MeshT>& bisection_;
		// std::shared_ptr<BisectionT> bisection_ptr_;


		bool node2elem_update_;
		EntitiesTuple tuple_entities_;
		// Array<bool, ManifoldDim+1> is_updated_;
		Integer level_;
	};






    template<typename...BaseFunctionSpaces,typename MeshT>
    void init(ConnectivitySimpliacialMapCollection<MeshT>& c) 
    {c.template init<BaseFunctionSpaces...>();}










template<typename FunctionSpace>
class FullFunctionSpaceInterpolation
{
public:

	FullFunctionSpaceInterpolation(std::shared_ptr<FunctionSpace> W_ptr):
	spaces_ptr_(W_ptr)
	{}

    inline void add_functionspace(std::shared_ptr<FunctionSpace> W_ptr){spaces_ptr_=W_ptr;}




    template<typename MeshT>
    void find_children(const MeshT& mesh, const Integer el)
    {
    	auto& fine_elem=mesh.elem(el);
      	auto& children=fine_elem.children;

      	for(std::size_t F_el=0;F_el<children.size();F_el++)
      	{
      		auto& fine_elem=mesh.elem(children[F_el]);
            find_children();

      	}
    }


	void init(const Integer C_level, const Integer F_level)
	{
      
      auto &mesh=spaces_ptr_->mesh();
      auto &bisection=spaces_ptr_->bisection();
      auto &tracker=bisection.tracker();


      for(std::size_t el=0;el<mesh.n_elements();el++)
      {
      	// if the elem does not belong to the C_level, continue
      	if(!elem_belongs_to_level(mesh,el,C_level,tracker))continue;
        
        auto& elem=mesh.elem(el);
        find_children(mesh,el);
      	// loop recursively on the fine children up to the F_level
      	// auto& children=elem.children;

      	// for(std::size_t F_el=0;F_el<children.size();F_el++)
      	// {
      	// 	auto& fine_elem=mesh.elem(children[F_el]);


      	// }


      }
	}

private:
 std::shared_ptr<FunctionSpace> spaces_ptr_;
};








	class ExactPoisson2D
	{
	public: 
	    // using Point=Matrix<Real,3,1>;
		using type=Matrix<Real,1,1>;
	    template<typename Point>
		static auto eval(const Point& p)
		{
			type func{2*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1])};
			return func; 
		}
	};

	class ExactPoisson3D
	{
	public: 
	    // using Point=Matrix<Real,3,1>;
		using type=Matrix<Real,1,1>;
	    template<typename Point>
		static auto eval(const Point& p)
		{
			type func{3*M_PI*M_PI*sin(M_PI*p[0])*sin(M_PI*p[1])*sin(M_PI*p[2])};
			return func; 
		}
	};

	template<Integer ManifoldDim>
	class ExactPoisson
	{
	public: 
		using type=Matrix<Real,1,1>;

	    template<typename Point>
		static type eval(const Point& p)
		{
			if(ManifoldDim==2)
				return ExactPoisson2D::eval(p); 
			else if(ManifoldDim==3)
				return ExactPoisson3D::eval(p); 
		}


	};



	template<Integer ManifoldDim,Integer Order>
	void Poisson_Lagrange(const Integer n)
	{
	  // constexpr Integer ManifoldDim=2;
		constexpr Integer Dim=ManifoldDim;
		using MeshT=Mesh<Dim, ManifoldDim>;
		MeshT mesh;
		if(ManifoldDim==2)
			generate_cube(mesh,n,n,0);
		else if(ManifoldDim==3)
			generate_cube(mesh,n,n,n);
		mesh.build_dual_graph();
		mark_boundary(mesh);
	   
		using AuxP_1= FunctionSpace< MeshT, Lagrange<Order,1>>;
	    clock_t begin = clock();
	    std::cout<<" init functionspaces "<<std::endl;
		AuxP_1 p1(mesh);
		auto Wtrial=MixedFunctionSpace(p1);
		auto Waux=AuxFunctionSpacesBuild(p1);
		auto W=FullSpaceBuild(Wtrial,Waux);
	    std::cout<<" end functionspaces  "<<std::endl;
	    clock_t end = clock();
	    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	    std::cout<<"elapsed_secs="<<elapsed_secs<<std::endl;

		auto u0 = MakeTrial<0>(W);
		auto v0 = MakeTest<0>(W);
		auto f1 = MakeFunction<0,ExactPoisson<ManifoldDim>>(W);


		std::cout<<W.spaces_ptr()->n_dofs()<<std::endl;
		constexpr auto C=Constant<Scalar>(1.0);


	  // 2D LSFEM POISSION
		auto bilinearform=L2Inner(Grad(u0),Grad(v0));
		auto linearform=L2Inner((v0),f1);
		auto bcs1=DirichletBC<0,FunctionZero<ManifoldDim>>(W,1);
		auto bcs2=DirichletBC<0,FunctionZero<ManifoldDim>>(W,2);
		auto bcs3=DirichletBC<0,FunctionZero<ManifoldDim>>(W,3);
		auto bcs4=DirichletBC<0,FunctionZero<ManifoldDim>>(W,4);
		auto bcs5=DirichletBC<0,FunctionZero<ManifoldDim>>(W,5);
		auto bcs6=DirichletBC<0,FunctionZero<ManifoldDim>>(W,6);
		std::cout<<" n_elements  "<<mesh.n_elements()<<std::endl;
	    std::cout<<" init context "<<std::endl;
		auto context=create_context(bilinearform,linearform,bcs1,bcs2,bcs3,bcs4,bcs5,bcs6);

		SparseMatrix<Real> A;
		std::vector<Real> b;
	    std::cout<<" init assembly A"<<std::endl;
	    begin = clock();
		context.assembly(A,b);
		std::cout<<" init assembly righthandside "<<std::endl;
		context.apply_bc(A,b);
	    end = clock();
	    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	    std::cout<<"elapsed_secs="<<elapsed_secs<<std::endl;
	 

		std::vector<Real> x;
		Integer max_iter=10000;


		std::cout<<"START SOLVING"<<std::endl;
		gauss_seidel(x,A,b,max_iter);
		std::cout<<"END SOLVING"<<std::endl;
		std::string output_file1 ="Poisson"+ std::to_string(ManifoldDim) +"D_P"+ std::to_string(Order)+"_output.vtk";
		auto var_names=variables_names("disp");

		std::ofstream os;
		os.open(output_file1.c_str());
		write_wtk(os,W,x,var_names);
		os.close();

		

	}


	template<Integer ManifoldDim,Integer Order1,Integer Order2>
	void LSFEM_Poisson(const Integer n)
	{
	  // constexpr Integer ManifoldDim=2;
		constexpr Integer Dim=ManifoldDim;
		using MeshT=Mesh<Dim, ManifoldDim>;
		MeshT mesh;
		if(ManifoldDim==2)
		{
			std::cout<<"2D lsfem poisson";
			// generate_cube(mesh,n,n,0);
			// mesh.build_dual_graph();
			// mark_boundary(mesh);
	      read_mesh("../data/triangle_square.MFEM", mesh);
	      assign_tags(mesh);


	      // Bisection<MeshT> bisection(mesh);
	      // bisection.uniform_refine(4);
	      mesh.update_dual_graph();
		}
		else if(ManifoldDim==3)
		{
			std::cout<<"3D lsfem poisson";
			generate_cube(mesh,n,n,n);
			mesh.build_dual_graph();
			mark_boundary(mesh);

	      // read_mesh("../data/tetrahedron.MFEM", mesh);
	      // assign_tags(mesh);
	      // Bisection<MeshT> bisection(mesh);
	      // bisection.uniform_refine(0);
	      // mesh.update_dual_graph();
		}
        Integer n_levels=12;
	    clock_t begin=clock();
		clock_t end = clock();
		double elapsed_secs;

		// std::cout<<"mat==="<<mat3<<std::endl;

		

        std::cout<<"n2em init="<<std::endl;
	    Node2ElemMap<MeshT> n2em(mesh);
	    n2em.init();

		std::cout<<"first child="<<mesh.elem(0).children.size()<<std::endl;
		Bisection<MeshT> bisection(mesh);

	    begin=clock();
		std::cout<<"bisec 1="<<std::endl;
		// bisection.tracking_begin();
		// bisection.uniform_refine(1);
		// bisection.tracking_end();
		end=clock();
		std::cout<<double(end - begin) / CLOCKS_PER_SEC<<std::endl;

		begin=clock();
		std::cout<<"bisec 2="<<std::endl;
		for(int i=0;i<n_levels;i++)
		{
		// bisection.tracking_begin();
		// bisection.uniform_refine(1);
		// bisection.tracking_end();			
		}

		end=clock();
		std::cout<<double(end - begin) / CLOCKS_PER_SEC<<std::endl;

		std::cout<<"dual 1="<<std::endl;
		mesh.update_dual_graph();

	    auto track=bisection.tracker();
	   
	    std::cout<<"___="<<std::endl;
	    auto& edge_node_map=bisection.edge_node_map();
	    auto& edge_element_map=bisection.edge_element_map();
	    std::cout<<"____="<<std::endl;
	    // edge_node_map.describe(std::cout);
	    // edge_element_map.describe(std::cout);

	 //    begin=clock();
	 //    std::cout<<"ens 1="<<std::endl;
	 //    EntitySimplicialMap<MeshT,1> enssim1(mesh);
		// end=clock();
		// std::cout<<double(end - begin) / CLOCKS_PER_SEC<<std::endl;


	 //    begin=clock();
	 //    std::cout<<"ens 2="<<std::endl;
	 //    EntitySimplicialMap<MeshT,2> enssim2(mesh);

		// end=clock();
		// std::cout<<double(end - begin) / CLOCKS_PER_SEC<<std::endl;


	 //    begin=clock();
	 //    std::cout<<"init 1="<<std::endl;
	 //    enssim1.init();
		// end=clock();
		// std::cout<<double(end - begin) / CLOCKS_PER_SEC<<std::endl;

	    // enssim1.describe(std::cout);
	    // enssim1.describe2(std::cout);
	 //    begin=clock();
	 //    std::cout<<"init 2="<<std::endl;
	 //    enssim2.init();
		// end=clock();
		// std::cout<<double(end - begin) / CLOCKS_PER_SEC<<std::endl;
	     
	    // enssim2.describe(std::cout);
	    // enssim2.describe2(std::cout);
	    
	  
	      
	    n2em.add_bisection(bisection);
	    n2em.init();
	    std::cout<<"--------------:mesh n elements"<<mesh.n_elements()<<std::endl;


	   //  for(int i=0;i<9;i++)
	   //  {
	   //  std::cout<<"node=="<<i<<std::endl;
	   //  std::cout<<"level=="<<0<<std::endl;
	   // auto tmp1=n2em.get(i,0);
	   //  for(int i=0;i<tmp1.size();i++)
	   //      std::cout<<tmp1[i]<<" ";
	   // std::cout<<std::endl;
	   // std::cout<<"level=="<<1<<std::endl;
	   // auto tmp2=n2em.get(i,1);
	   //  for(int i=0;i<tmp2.size();i++)
	   //      std::cout<<tmp2[i]<<" ";
	   // std::cout<<std::endl;
	   // std::cout<<"level=="<<2<<std::endl;
	   // auto tmp3=n2em.get(i,2);
	   //  for(int i=0;i<tmp3.size();i++)
	   //      std::cout<<tmp3[i]<<" "; 
	   // std::cout<<std::endl;
	   // }

	 //    const auto& simplex1=mesh.elem(0);
	 //    const auto& simplex2=mesh.elem(0);
		// const auto& simplex3=mesh.elem(2);


	    // std::cout<<is_simplex_inside_simplex(simplex1,simplex2,mesh)<<std::endl;

	    
	    // ConnectivitySimpliacialMap<MeshT,0> csm0(mesh,n2em,bisection);
	    // ConnectivitySimpliacialMap<MeshT,1> csm1(mesh,n2em,bisection);
	    // ConnectivitySimpliacialMap<MeshT,2> csm2(mesh,n2em,bisection);
	    // csm1.init();
	    // csm0.init();
	    // csm2.init();
	    // std::array<Integer,2> edge{0,1};
	    
	    // csm1.describe(std::cout);
	    // csm1.describe2(std::cout);
	    // std::cout<<csm.get(simplex3,edge,1)<<std::endl;

	    // edge[0]=0;edge[1]=2;
	    // std::cout<<csm.get(simplex3,edge,1)<<std::endl;
	    // csm1.init();
	    // csm0.init();
	    // csm2.init();
	    // edge[0]=1;edge[1]=4;
	    // std::cout<<"__________________________________________________"<<std::endl;
	    // std::cout<<csm.get(simplex3,edge,1)<<std::endl;
	    // std::cout<<"__________________________________________________"<<std::endl;
	    // edge[0]=3;edge[1]=4;
	    // const auto& simplex4=mesh.elem(4);
	    // std::cout<<csm.get(simplex4,edge,1)<<std::endl;
	    // std::cout<<"__________________________________________________"<<std::endl;
	    
	    // std::cout<<csm.get(simplex3,edge,1)<<std::endl;


	    // // std::cout<<" entity 0 describe __________________________________________________"<<std::endl;
	    // // csm0.describe(std::cout);
	    // std::cout<<" entity 0 describe 2 __________________________________________________"<<std::endl;
	    // csm0.describe2(std::cout);
	    // std::cout<<" entity 0 describe 3 __________________________________________________"<<std::endl;
	    // csm0.describe3(std::cout);

	    // // std::cout<<" entity 1 describe __________________________________________________"<<std::endl;
	    // // csm1.describe(std::cout);
	    // std::cout<<" entity 1 describe 2 __________________________________________________"<<std::endl;
	    // csm1.describe2(std::cout);
	    // std::cout<<" entity 1 describe 3 __________________________________________________"<<std::endl;
	    // csm1.describe3(std::cout);


	    // std::cout<<" entity 2 describe __________________________________________________"<<std::endl;
	    // csm2.describe(std::cout);
	    // std::cout<<" entity 2 describe 2 __________________________________________________"<<std::endl;
	    // csm2.describe2(std::cout);
	    // std::cout<<" entity 2 describe 3 __________________________________________________"<<std::endl;
	    // csm2.describe3(std::cout);

	    




	    
	    // auto kok= std::make_tuple(ConnectivitySimpliacialMap<MeshT,0>(mesh,n2em,bisection));
	    // auto kok= std::make_tuple(ConnectivitySimpliacialMap<MeshT,0>(mesh,n2em,bisection));

	    // std::cout<<"____1"<<std::endl;
	    // auto tmp1=n2em.get(0,0);
	    // for(int i=0;i<tmp1->size();i++)
	    //     std::cout<<(*tmp1)[i]<<" ";
	    // std::cout<<std::endl;
	    // std::cout<<"____2"<<std::endl;
	    // auto tmp2=n2em.get(1,0);
	    // for(int i=0;i<tmp2->size();i++)
	    //     std::cout<<(*tmp2)[i]<<" ";
	    // std::cout<<std::endl;
	    // std::cout<<"____3"<<std::endl;
	    // auto tmp3=n2em.get(2,0);
	    // for(int i=0;i<tmp3->size();i++)
	    //     std::cout<<(*tmp3)[i]<<" ";
	    // std::cout<<std::endl;
	    // std::cout<<"____4"<<std::endl;
	    // auto tmp4=n2em.get(3,0);
	    // for(int i=0;i<tmp4->size();i++)
	    //     std::cout<<(*tmp4)[i]<<" ";
	    // std::cout<<std::endl;

	    // n2em.init();
	    // std::cout<<"5"<<std::endl;
	    // auto tmp5=n2em.get(4,1);
	    // for(int i=0;i<tmp5->size();i++)
	    //     std::cout<<(*tmp5)[i]<<" ";
	    // std::cout<<std::endl;
	// faccio il count delle eentita' di interesse sulla mesh coarse
	// qui se itero sull'iteratore e ho un count=0, con count++,
	// count rappresente il numero dell'entita
	// 
	// considero un geenerico refinement di cui voglio ricalcolare il numero dell'entita'
	// devo per forza loopare dai primi agli ultimi elementi o comunque prima i padri poi i figli
	// fissato l'elemento, controllo se non ha figli. in tal caseo le eentita sono mantenute uguali.
	// se invece l'elemento ha figli, devo ricreare le enetita'. ma devo ricreare unicamente quelle nuove.
	// quindi devo capire cosa vuol dire avere dei figli.    
	// ogniqualvolta faccio local refinement, loop su tutte le entita' locali. 
	// cerco se esiste gia' la mappa corrispondente. se esiste, non faccio nulla.
	// se non esiste, creo una nuova entita' e aggiorno il contatore dell'entita' di quel livello.
	// il punto e' che se un elemento 
	  // for(int i=0;i<mesh.n_elements();i++)
	  // {
	  // 	auto el=mesh.elem(i);
	  // 	auto& child=el.children;
	  // 	auto& parent_id=el.parent_id;
	  // 	std::cout<<"parent_id="<<parent_id<<std::endl;
	  // 	std::cout<<"children="<<std::endl;
	  // 	for(int j=0;j<child.size();j++)
	  // 	{
	  //    std::cout<<child[j]<<" ";
	  // 	}
	  // 	std::cout<<std::endl;

	  // }





	 //    bisec.uniform_refine(1);
		// bisec.clear();

		// mesh.clean_up();
		// mesh.reorder_nodes();

		// Quality<MeshT> q(mesh);
		// q.compute();

		// mesh.update_dual_graph();
		// mark_boundary(mesh);

		// auto edge_select = std::make_shared<ImplicitOrderEdgeSelect<Mesh2>>();
		// bisec.set_edge_select(edge_select);
		// bisec.uniform_refine(1);


	  // mesh.build_dual_graph();
	  // mark_boundary(mesh);






	 //    MeshAndEntity<MeshT> mesh_and_entity0(mesh,0);
	 //    MeshAndEntity<MeshT> mesh_and_entity1(mesh,1);
	 //    MeshAndEntity<MeshT> mesh_and_entity2(mesh,2);
	 //    mesh_and_entity0.add_bisection(bisection);
	 //    mesh_and_entity1.add_bisection(bisection);
	 //    mesh_and_entity2.add_bisection(bisection);

		using AuxRT_n= FunctionSpace< MeshT, RT<Order1,1>>;
		using AuxP_n= FunctionSpace< MeshT, Lagrange<Order2,1>>;
		using LSFEM= FunctionSpace< MeshT, RT<Order1,1>,Lagrange<Order2,1>>;

	    ConnectivitySimpliacialMapCollection<MeshT> csmc(mesh,n2em,bisection);
	    csmc.template init<RT<Order1,1>,Lagrange<Order2,1>>();

	    csmc.describe2(std::cout);
	 //   std::vector<Integer> list_of_elems;
	 //    get_elems_of_level(list_of_elems,mesh,given_level,track);


	 //    std::cout<<"list_of_elems"<<std::endl;	
	 //    for(Integer i=0;i<list_of_elems.size();i++)
	 //    	std::cout<<list_of_elems[i]<<std::endl;	
	 //    std::cout<<"current_iterate="<<track.current_iterate()<<std::endl;



// AuxRT_n rtn(mesh);
// AuxP_n pn(mesh);
// LSFEM lsfem(mesh);
AuxRT_n rtn(mesh,bisection,n2em);//csmc);
AuxP_n pn(mesh,bisection,n2em);//csmc);
std::cout<<"FIRST PRE UPDATE="<<std::endl;
LSFEM lsfem(mesh,bisection,n2em);//csmc);
std::cout<<"FIRST POST UPDATE="<<std::endl;

	 // // AuxRT_n rtn(mesh);
		// std::cout<<" p n init  = "<<std::endl;
		// AuxP_n pn(mesh_and_entity0);
		// std::cout<<" lsfem init  = "<<std::endl;
		// LSFEM lsfem(mesh_and_entity0);

		// auto Wtrial=MixedFunctionSpace(lsfem);
	 auto Wtrial=MixedFunctionSpace(rtn,pn);

		auto Waux=AuxFunctionSpacesBuild(pn);

		auto W=FullSpaceBuild(Wtrial,Waux);
		auto W_ptr=std::make_shared<decltype(W)>(W);


		for(int i=0;i<n_levels;i++)
		{
		bisection.tracking_begin();
		bisection.uniform_refine(1);
		bisection.tracking_end();			
		}
 



	// auto &dofsdofmap=W.dofsdofmap();
	// std::cout<<"<<<< UPDATE >>>>"<<std::endl;
 //    // dofsdofmap.update();
 //    std::cout<<"PRE UPDATE="<<std::endl;
 //    W.update();

    
 //    auto& level_cumultive_n_dofs=dofsdofmap.level_cumultive_n_dofs();
 //     std::cout<<"level_n_dofs_array="<<std::endl;
 //     for(int i=0;i<level_cumultive_n_dofs.size();i++)
 //     {
 //      // for(int j=0; j<level_n_dofs_array[i].size();j++)
 //        std::cout<<level_cumultive_n_dofs[i]<<" ";
 //      std::cout<<std::endl;
 //     }   
 //    std::cout<<"POST UPDATE="<<std::endl;



	// std::cout<<"dofsdofmap.level_cumulative_dofs_array()"<<std::endl;
	// auto& h1=dofsdofmap.level_cumulative_dofs_array();
	// for(Integer i=0;i<h1.size();i++)
	// {

	// 	for(Integer j=0;j<h1[i].size();j++)
	// 		std::cout<<h1[i][j]<<" ";
	// 	std::cout<<std::endl;
	// }



	// std::cout<<"dofsdofmap.level_cumulative_dofs_array()"<<std::endl;
	// auto& h2=dofsdofmap.level_cumulative_dofs_array();
	// for(Integer i=0;i<h2.size();i++)
	// {

	// 	for(Integer j=0;j<h2[i].size();j++)
	// 		std::cout<<h2[i][j]<<" ";
	// 	std::cout<<std::endl;
	// }



	// 	begin=clock();
	// 	std::cout<<"OTHER BISECTION="<<std::endl;
		// for(int i=n_levels-1;i<n_levels;i++)
		// {
		// bisection.tracking_begin();
		// bisection.uniform_refine(1);
		// bisection.tracking_end();			
		// }
     
    // auto& m1=lsfem.dofsdofmap().level_n_dofs_array();

    //  for(std::size_t i=0; i<m1.size() ;i++)
    //   { 
    //     // loop on the levels
    //     for(std::size_t j=0;j<n_levels;j++)
    //        {std::cout<<m1[i][j]<<" ";}
    //      std::cout<<std::endl;
    //    }





	//     std::cout<<"-----dofmap of space="<<0<<std::endl;
	// GetType<typename decltype(W)::DofsDM::ElemDofMap> elemdm0;
	// for(Integer i=0;i<dofsdofmap.template dofmap_size<0>();i++)
	// {
	//          dofsdofmap.template dofmap_get<0>(elemdm0,i,level);
	//          std::cout<<elemdm0<<std::endl;
	// }



	// std::cout<<"-----dofmap of space="<<1<<std::endl;
	// GetType<typename decltype(W)::DofsDM::ElemDofMap,1> elemdm1;
	// for(Integer i=0;i<dofsdofmap.template dofmap_size<1>();i++)
	// {
	//          dofsdofmap.template dofmap_get<1>(elemdm1,i,level);
	//          std::cout<<elemdm1<<std::endl;
	// }


 //    auto& spacedofs=dofsdofmap.space_dofs();
 //    auto& n_dofs=dofsdofmap.n_dofs();
 //    auto& cumulative_n_dofs=dofsdofmap.cumulative_n_dofs();

 //    std::cout<<"spacedofs"<<std::endl;
 //   	for(Integer i=0;i<spacedofs.size();i++)
 //   		{for(Integer j=0;j<spacedofs[i]->size();j++)
	//          std::cout<<(*spacedofs[i])[j]<<" ";
	// std::cout<<std::endl;}

 //    std::cout<<"n_dofs"<<std::endl;
 //   	for(Integer i=0;i<n_dofs.size();i++)
	//          std::cout<<n_dofs[i]<<" ";
	// std::cout<<std::endl;

 //    std::cout<<"cumulative_n_dofs"<<std::endl;
 //   	for(Integer i=0;i<cumulative_n_dofs.size();i++)
 //   		{
	//          std::cout<<cumulative_n_dofs[i]<<" ";
	// std::cout<<std::endl;
 //        }






	    //  std::cout<<"elems "<<std::endl;
	    //  for(Integer i=0;i<mesh.n_elements();i++)
	    //    {std::cout<<"elem = "<<i<< std::endl;
	    //     for(Integer j=0;j<mesh.elem(i).nodes.size();j++)
	    // 	std::cout<<mesh.elem(i).nodes[j]<< " ";
	    //     std::cout<<std::endl;
	    //     // std::cout<<std::endl;
	    // }	



      // auto arr_lev=W.level_array_ndofs();
      // auto arr=W.array_ndofs();
      
      // std::cout<<"arr_lev"<<std::endl;
      // for(int i=0;i<arr_lev.size();i++)
      // {
      // 	std::cout<<std::endl;
      // 	for(int j=0;j<arr_lev[i].size();j++)
      // 		std::cout<<arr_lev[i][j]<<" ";
      // 	std::cout<<std::endl;
      // }

      // std::cout<<"arr"<<std::endl;
      // for(int i=0;i<arr.size();i++)
      // {
      // 	std::cout<<std::endl;
      // 		std::cout<<arr[i]<<" ";
      // 	std::cout<<std::endl;
      // }

	 //    init<Simplex<Dim,ManifoldDim>,Lagrange<Order2,1>>(mesh_and_entity0);
	 //    init<Simplex<Dim,ManifoldDim>,RT<Order1,1>>(mesh_and_entity0);

	    // mesh_and_entity0.init(pn);
	    // mesh_and_entity0.init(lsfem);
	    // mesh_and_entity0.init(pn);
	    // mesh_and_entity0.init(pn);
	    // mesh_and_entity0.init(lsfem);
	    // mesh_and_entity0.init(lsfem);

	    // mesh_and_entity1.init(pn);
	    // mesh_and_entity1.init(lsfem);
	    // mesh_and_entity1.init(pn);
	    // mesh_and_entity1.init(pn);
	    // mesh_and_entity1.init(lsfem);
	    // mesh_and_entity1.init(lsfem);
	    // ElemEntityCollection<MeshT> entitycollection0(mesh,0);
	    // ElemEntityCollection<MeshT> entitycollection1(mesh,1);
	    // ElemEntityCollection<MeshT> entitycollection2(mesh,2);
	    // entitycollection0.add_bisection(bisection);
	    // entitycollection1.add_bisection(bisection);
	    // entitycollection2.add_bisection(bisection);
	    // entitycollection0.init(W);
	    // entitycollection1.init(W);
	    // entitycollection2.init(W);

	    // std::cout<<"coll0 ---------- "<<std::endl;
	    // init<0>(entitycollection0);
	    // init<0>(entitycollection1);
	    // init<0>(entitycollection2);
	    // std::cout<<"coll1 ---------- "<<std::endl;
	    // init<1>(entitycollection0);
	    // init<1>(entitycollection1);
	    // init<1>(entitycollection2);
	    // std::cout<<"coll2 ---------- "<<std::endl;
	    // init<2>(entitycollection0);
	    // init<2>(entitycollection1);
	    // init<2>(entitycollection2);


	    // const auto& n2e0=entitycollection0.node2elem();
	    // std::cout<<"node 2 elems LEVEL 0"<<std::endl;	
	    // for(Integer i=0;i<n2e0.size();i++)
	    // 	{for(Integer j=0;j<n2e0[i].size();j++)
	    // 	std::cout<<n2e0[i][j]<< " ";
	    //     std::cout<<std::endl;}	

	    // const auto& n2e1=entitycollection1.node2elem();
	    // std::cout<<"node 2 elems LEVEL 1"<<std::endl;	
	    // for(Integer i=0;i<n2e1.size();i++)
	    // 	{for(Integer j=0;j<n2e1[i].size();j++)
	    // 	std::cout<<n2e1[i][j]<< " ";
	    //     std::cout<<std::endl;}	

	    // const auto& n2e2=entitycollection2.node2elem();
	    // std::cout<<"node 2 elems LEVEL 2"<<std::endl;	
	    // for(Integer i=0;i<n2e2.size();i++)
	    // 	{for(Integer j=0;j<n2e2[i].size();j++)
	    // 	std::cout<<n2e2[i][j]<< " ";
	    //     std::cout<<std::endl;}	

	    //  std::cout<<"node "<<std::endl;
	    //  for(Integer i=0;i<mesh.points().size();i++)
	    //    {std::cout<<"node = "<<i<< std::endl;
	    //     // for(Integer j=0;j<mesh.points()[0].size();j++)
	    // 	std::cout<<mesh.points()[i]<< std::endl;
	    //     // std::cout<<std::endl;
	    // }	   

	    //  std::cout<<"elems "<<std::endl;
	    //  for(Integer i=0;i<mesh.n_elements();i++)
	    //    {std::cout<<"elem = "<<i<< std::endl;
	    //     for(Integer j=0;j<mesh.elem(i).nodes.size();j++)
	    // 	std::cout<<mesh.elem(i).nodes[j]<< " ";
	    //     std::cout<<std::endl;
	    //     // std::cout<<std::endl;
	    // }		
	    // std::cout<<"define ens---------- "<<std::endl;

	    // const auto& ens00=get_entity<0>(mesh_and_entity0);
	    // const auto& ens01=get_entity<1>(mesh_and_entity0);
	    // const auto& ens02=get_entity<2>(mesh_and_entity0);
	    // const auto& ens10=get_entity<0>(mesh_and_entity1);
	    // const auto& ens11=get_entity<1>(mesh_and_entity1);
	    // const auto& ens12=get_entity<2>(mesh_and_entity1);



	    // std::cout<<"---------- "<<std::endl;
	    // const auto e2e00=ens00.elem_2_entity();
	    // const auto e2e01=ens01.elem_2_entity();
	    // const auto e2e02=ens02.elem_2_entity();
	    // std::cout<<"---------- "<<std::endl;
	    // const auto e2e10=ens10.elem_2_entity();
	    // const auto e2e11=ens11.elem_2_entity();
	    // const auto e2e12=ens12.elem_2_entity();

	    //  std::cout<<std::endl;
	    //  std::cout<<"e2e00 size="<<e2e00.size()<<std::endl;
	    //  for(Integer i=0;i<e2e00.size();i++)
	    //    {
	    //    	if(e2e00[i].size()>0)
	    //    	{std::cout<<"elem = "<<i<<std::endl;
	    //     for(Integer j=0;j<e2e00[i].size();j++)
	    // 	std::cout<<e2e00[i][j]<< " ";
	    //     std::cout<<std::endl;
	    // }
	    //     // std::cout<<std::endl;
	    // }
	    //  std::cout<<std::endl;
	    //  std::cout<<"e2e01 size="<<e2e00.size()<<std::endl;
	    //  for(Integer i=0;i<e2e01.size();i++)
	    //    {
	    //    	if(e2e01[i].size()>0)
	    //    	{std::cout<<"elem = "<<i<< std::endl;
	    //     for(Integer j=0;j<e2e01[i].size();j++)
	    // 	std::cout<<e2e01[i][j]<< " ";
	    //     std::cout<<std::endl;}
	    //     // std::cout<<std::endl;
	    // }
	    //  std::cout<<std::endl;
	    //  std::cout<<"e2e02 size="<<e2e02.size()<<std::endl;
	    //  for(Integer i=0;i<e2e02.size();i++)
	    //    {
	    //    	if(e2e02[i].size()>0){
	    //    	std::cout<<"elem = "<<i<<std::endl;
	    //     for(Integer j=0;j<e2e02[i].size();j++)
	    // 	std::cout<<e2e02[i][j]<< " ";
	    //     std::cout<<std::endl;}
	    //     // std::cout<<std::endl;
	    // }

	    //  std::cout<<std::endl;
	    //  std::cout<<"e2e10 size="<<e2e10.size()<<std::endl;
	    //  for(Integer i=0;i<e2e10.size();i++)
	    //    {
	    //    	if(e2e10[i].size()>0){
	    //    	std::cout<<"elem = "<<i<< std::endl;
	    //     for(Integer j=0;j<e2e10[i].size();j++)
	    // 	std::cout<<e2e10[i][j]<< " ";
	    //     std::cout<<std::endl;
	    // }
	    //     // std::cout<<std::endl;
	    // }
	    //  std::cout<<std::endl;
	    //  std::cout<<"e2e11 size="<<e2e11.size()<<std::endl;
	    //  for(Integer i=0;i<e2e11.size();i++)
	    //    {
	    //     if(e2e11[i].size()>0){
	    //    	std::cout<<"elem = "<<i<< std::endl;
	    //     for(Integer j=0;j<e2e11[i].size();j++)
	    // 	std::cout<<e2e11[i][j]<< " ";
	    //     std::cout<<std::endl;
	    // }
	    //     // std::cout<<std::endl;
	    // }
	    //  std::cout<<std::endl;
	    //  std::cout<<"e2e12 size="<<e2e12.size()<<std::endl;
	    //  for(Integer i=0;i<e2e02.size();i++)
	    //    {
	    //     if(e2e11[i].size()>0){
	    //    	std::cout<<"elem = "<<i<< std::endl;
	    //     for(Integer j=0;j<e2e12[i].size();j++)
	    // 	std::cout<<e2e12[i][j]<< " ";
	    //     std::cout<<std::endl;
	    // }
	    //     // std::cout<<std::endl;
	    // }


	    // std::cout<<"list_of_elems"<<std::endl;	
	    // for(Integer i=0;i<mesh.n_elements();i++)
	    // 	std::cout<<"elem="<<i<<" with n child = "<<mesh.elem(i).children.size()<<" parent id = "<< mesh.elem(i).parent_id<< std::endl;	


	    // for(Integer i=0;i<mesh.n_elements();i++)
	    // 	std::cout<<track.get_iterate(i)<<std::endl;	

	 

	// auto &dofsdofmap=W.dofsdofmap();


	//     std::cout<<"-----dofmap of space="<<0<<std::endl;
	// GetType<typename decltype(W)::DofsDM::ElemDofMap> elemdm0;
	// for(Integer i=0;i<dofsdofmap.template dofmap_size<0>();i++)
	// {
	//          dofsdofmap.template dofmap_get<0>(elemdm0,i);
	//          std::cout<<elemdm0<<std::endl;
	// }



	// std::cout<<"-----dofmap of space="<<1<<std::endl;
	// GetType<typename decltype(W)::DofsDM::ElemDofMap,1> elemdm1;
	// for(Integer i=0;i<dofsdofmap.template dofmap_size<1>();i++)
	// {
	//          dofsdofmap.template dofmap_get<1>(elemdm1,i);
	//          std::cout<<elemdm1<<std::endl;
	// }

	// std::cout<<"-----dofmap of space="<<2<<std::endl;

	// for(Integer i=0;i<dofsdofmap.template dofmap_size<2>();i++)
	// {
	//          dofsdofmap.template dofmap_get<2>(elemdm,i);
	//          std::cout<<elemdm<<std::endl;
	// }







		auto u0 = MakeTrial<0>(W_ptr);
		auto u1 = MakeTrial<1>(W_ptr);

		auto v0 = MakeTest<0>(W_ptr);
		auto v1 = MakeTest<1>(W_ptr);

		auto f1 = MakeFunction<0,ExactPoisson<ManifoldDim>>(W_ptr);

		std::cout<<W.spaces_ptr()->n_dofs()<<std::endl;

	  // 2D LSFEM POISSION
		auto bilinearform=L2Inner(Div(u0),Div(v0))+
		L2Inner((u0),(v0))+
		L2Inner(Grad(u1),Grad(v1))-
		L2Inner(Grad(u1),(v0))-
		L2Inner((u0),Grad(v1));
	  // auto bilinearform=//L2Inner(Div(u0),Div(v0))+
	  //                   L2Inner((u0),(v0));   

		auto linearform=
		L2Inner(-Div(v0),f1);




		auto bcs1=DirichletBC<1,FunctionZero2D>(W_ptr,1);
		auto bcs2=DirichletBC<1,FunctionZero2D>(W_ptr,2);
		auto bcs3=DirichletBC<1,FunctionZero2D>(W_ptr,3);
		auto bcs4=DirichletBC<1,FunctionZero2D>(W_ptr,4);
		auto bcs5=DirichletBC<1,FunctionZero2D>(W_ptr,5);
		auto bcs6=DirichletBC<1,FunctionZero2D>(W_ptr,6);








		std::cout<<"CREATE CONTEXT"<<std::endl;
		auto context=create_context(bilinearform,linearform,bcs1,bcs2,bcs3,bcs4,bcs5,bcs6);





    auto &dofsdofmap=W_ptr->dofsdofmap();
	// auto &dofsdofmap=context.full_spaces_ptr()->dofsdofmap();
	std::cout<<"<<<< UPDATE >>>>"<<std::endl;
    // dofsdofmap.update();
    std::cout<<"PRE UPDATE="<<std::endl;
    // context.full_spaces_ptr()->update();
     W_ptr->update();
    
    auto& level_cumultive_n_dofs=dofsdofmap.level_cumultive_n_dofs();
     std::cout<<"level_n_dofs_array="<<std::endl;
     for(int i=0;i<level_cumultive_n_dofs.size();i++)
     {
      // for(int j=0; j<level_n_dofs_array[i].size();j++)
        std::cout<<level_cumultive_n_dofs[i]<<" ";
      std::cout<<std::endl;
     }   
    std::cout<<"POST UPDATE="<<std::endl;






		SparseMatrix<Real> A;
		std::vector<Real> b;

		std::cout<<"ASSEMBLY"<<std::endl;
		Integer level=0;//bisection.tracker().current_iterate();
		// Integer level=4;
		std::cout<<"level---"<<level<<std::endl;
		context.assembly(A,b,level);
	   // A.print_val();
		std::cout<<"APPLY BC "<<std::endl;
		context.apply_bc(A,b);


	   // A.print_val();
		std::vector<Real> x;
		Integer max_iter=12;


		std::cout<<"START SOLVING"<<std::endl;
		gauss_seidel(x,A,b,max_iter);
		std::cout<<"END SOLVING"<<std::endl;
		std::ofstream os;
		auto var_names=variables_names("stress","disp");
		std::string output_file ="LSFEM_Poisson"+ std::to_string(ManifoldDim) +
		"D_RT" + std::to_string(Order1)+
		"_P" + std::to_string(Order2)+"_output.vtk";

		os.close();
		os.open(output_file.c_str());
		// write_wtk_isoparametric(os,*context.full_spaces_ptr(),x,var_names,level);
		write_wtk_isoparametric(os,W_ptr,x,var_names,level);

	    os.close();


	    // std::cout<<"x"<<std::endl;
	    // for(int i=0;i<x.size();i++)
	    // 	std::cout<<x[i]<<std::endl;

     //    std::cout<<std::endl;
	    // std::cout<<"----------b"<<std::endl;
	    // for(int i=0;i<b.size();i++)
	    // 	std::cout<<b[i]<<std::endl;



        


		SparseMatrix<Real> AL;
		std::vector<Real> bL;

		level=bisection.tracker().current_iterate()-1;

		context.assembly(AL,bL,level);
	   // A.print_val();
		std::cout<<"APPLY BC "<<std::endl;
		context.apply_bc(AL,bL);


	   // A.print_val();
		std::vector<Real> xL;


		std::cout<<"START SOLVING"<<std::endl;
		gauss_seidel(xL,AL,bL,max_iter);
		std::cout<<"END SOLVING"<<std::endl;
		auto var_namesL=variables_names("stress","disp");
		std::string output_fileL ="LSFEM_Poisson"+ std::to_string(ManifoldDim) +
		"D_RT" + std::to_string(Order1)+
		"_P" + std::to_string(Order2)+"_outputFINE.vtk";

		os.close();
		os.open(output_fileL.c_str());
		// write_wtk_isoparametric(os,*context.full_spaces_ptr(),xL,var_namesL,level);
		write_wtk_isoparametric(os,W_ptr,xL,var_namesL,level);

	    os.close();


	    // std::cout<<"xL"<<std::endl;
	    // for(int i=0;i<xL.size();i++)
	    // 	std::cout<<xL[i]<<std::endl;
	    // std::cout<<"bL"<<std::endl;
	    // for(int i=0;i<bL.size();i++)
	    // 	std::cout<<bL[i]<<std::endl;




	auto &dofsdofmap2=W_ptr->dofsdofmap();


    
    auto& level_cumultive_n_dofs2=dofsdofmap2.level_cumultive_n_dofs();
     std::cout<<"level_n_dofs_array="<<std::endl;
     for(int i=0;i<level_cumultive_n_dofs2.size();i++)
     {
      // for(int j=0; j<level_n_dofs_array[i].size();j++)
        std::cout<<level_cumultive_n_dofs2[i]<<" ";
      std::cout<<std::endl;
     }   




	//    // std::cout<<"W ndofs"<<W.spaces_ptr()->n_dofs()<<std::endl;
	//    // std::cout<<"W array_ndofs"<<W.array_ndofs()<<std::endl;
	//    // std::cout<<ShapeFunction<Simplex<Dim,ManifoldDim>,RT1<1> ,IdentityOperator,GaussPoints<Simplex<3,3>,3>>::reference_values()<<std::endl;
	//    // ReferenceShapeFunctionValue<Simplex<3,3>, DivergenceOperator, RaviartThomasFE, 1>::apply2(Vector<Real,ManifoldDim> {0.25,0.25,0.25});
	//    // ReferenceShapeFunctionValue<Simplex<3,3>, DivergenceOperator, RaviartThomasFE, 1>::apply2(Vector<Real,ManifoldDim> {0.0,0.0,0.0});
	//    // ReferenceShapeFunctionValue<Simplex<3,3>, DivergenceOperator, RaviartThomasFE, 1>::apply2(Vector<Real,ManifoldDim> {1.0,0.0,0.0});
	//    // ReferenceShapeFunctionValue<Simplex<3,3>, DivergenceOperator, RaviartThomasFE, 1>::apply2(Vector<Real,ManifoldDim> {0.0,1.0,0.0});
	//    // ReferenceShapeFunctionValue<Simplex<3,3>, DivergenceOperator, RaviartThomasFE, 1>::apply2(Vector<Real,ManifoldDim> {0.0,0.0,1.0});
	//    // std::cout<<"Vector<Real,ManifoldDim> {0.25,0.25,0.25}"<<std::endl;
	//    // ReferenceShapeFunctionValue<Simplex<3,3>, IdentityOperator, RaviartThomasFE, 1>::apply2(Vector<Real,ManifoldDim> {0.25,0.25,0.25});
	//    // std::cout<<"Vector<Real,ManifoldDim> {0.0,0.0,0.0}"<<std::endl;
	//    // ReferenceShapeFunctionValue<Simplex<3,3>, IdentityOperator, RaviartThomasFE, 1>::apply2(Vector<Real,ManifoldDim> {0.0,0.0,0.0});
	//    // std::cout<<"Vector<Real,ManifoldDim> {1.0,0.0,0.0}"<<std::endl;
	//    // ReferenceShapeFunctionValue<Simplex<3,3>, IdentityOperator, RaviartThomasFE, 1>::apply2(Vector<Real,ManifoldDim> {1.0,0.0,0.0});
	//    // std::cout<<"Vector<Real,ManifoldDim> {0.0,1.0,0.0}"<<std::endl;
	//    // ReferenceShapeFunctionValue<Simplex<3,3>, IdentityOperator, RaviartThomasFE, 1>::apply2(Vector<Real,ManifoldDim> {0.0,1.0,0.0});
	//    // std::cout<<"Vector<Real,ManifoldDim> {0.0,0.0,1.0}"<<std::endl;
	//    // ReferenceShapeFunctionValue<Simplex<3,3>, IdentityOperator, RaviartThomasFE, 1>::apply2(Vector<Real,ManifoldDim> {0.0,0.0,1.0});

	// auto qp=GaussPoints<Simplex<Dim,3>,3>::qp_points;
	// std::cout<<"point 0 " <<std::endl;
	// ReferenceShapeFunctionValue<Simplex<3,3>, IdentityOperator, RaviartThomasFE, 1>::apply2(Vector<Real,3>{qp(0,0),qp(0,1),qp(0,2)});
	// std::cout<<"point 1 " <<std::endl;
	// ReferenceShapeFunctionValue<Simplex<3,3>, IdentityOperator, RaviartThomasFE, 1>::apply2(Vector<Real,3>{qp(1,0),qp(1,1),qp(1,2)});
	// std::cout<<"point 2 " <<std::endl;

	// ReferenceShapeFunctionValue<Simplex<3,3>, IdentityOperator, RaviartThomasFE, 1>::apply2(Vector<Real,3>{qp(2,0),qp(2,1),qp(2,2)});
	// std::cout<<"point 3 " <<std::endl;

	// ReferenceShapeFunctionValue<Simplex<3,3>, IdentityOperator, RaviartThomasFE, 1>::apply2(Vector<Real,3>{qp(3,0),qp(3,1),qp(3,2)});
	// std::cout<<"point 4 " <<std::endl;

	// ReferenceShapeFunctionValue<Simplex<3,3>, IdentityOperator, RaviartThomasFE, 1>::apply2(Vector<Real,3>{qp(4,0),qp(4,1),qp(4,2)});

	// std::cout<<"nlocaldofs="<<decltype(W)::DofsDM::NLocalDofs<<std::endl;
	// auto &dofsdofmap=W.dofsdofmap();

	// for(Integer i=0;i<dofsdofmap.space_dofs().size();i++)
	// {
	//     std::cout<<"===== space="<<i<<std::endl;
	//     for(Integer j=0;j<dofsdofmap.space_dofs_size(i);j++)
	//     {
	//       // for(Integer k=0;k<sd3[i]->operator[](j).size();k++)
	//          std::cout<<dofsdofmap.space_dofs_get(i,j)<<std::endl;
	//     }

	// }

	//     std::cout<<"-----dofmap of space="<<0<<std::endl;
	// GetType<typename decltype(W)::DofsDM::ElemDofMap> elemdm;
	// for(Integer i=0;i<dofsdofmap.template dofmap_size<0>();i++)
	// {
	//          dofsdofmap.template dofmap_get<0>(elemdm,i);
	//          std::cout<<elemdm<<std::endl;
	// }



	// std::cout<<"-----dofmap of space="<<1<<std::endl;

	// for(Integer i=0;i<dofsdofmap.template dofmap_size<1>();i++)
	// {
	//          dofsdofmap.template dofmap_get<1>(elemdm,i);
	//          std::cout<<elemdm<<std::endl;
	// }

	// std::cout<<"-----dofmap of space="<<2<<std::endl;

	// for(Integer i=0;i<dofsdofmap.template dofmap_size<2>();i++)
	// {
	//          dofsdofmap.template dofmap_get<2>(elemdm,i);
	//          std::cout<<elemdm<<std::endl;
	// }
	// std::cout<<"-----dofmap of space="<<3<<std::endl;

	// for(Integer i=0;i<dofsdofmap.template dofmap_size<3>();i++)
	// {
	//          dofsdofmap.template dofmap_get<3>(elemdm,i);
	//          std::cout<<elemdm<<std::endl;
	// }
		// std::cout<<" cumulative_n_dofs="<<W.dofsdofmap().cumulative_n_dofs()<<std::endl;
		// std::cout<<" n_dofs="<<W.dofsdofmap().n_dofs()<<std::endl;
		// std::cout<<" n_dofs="<<W.dofsdofmap().NsubspacesArray<<std::endl;
	std::cout<<"end "<<std::endl;
	}



	// void LSFEM_Poisson2D_RT0_P1()
	// {
	// 	constexpr Integer ManifoldDim=2;
	// 	constexpr Integer Dim=2;
	// 	using MeshT=Mesh<Dim, ManifoldDim>;

	// 	MeshT mesh;
	// 	generate_square(mesh,20,20);
	// 	mesh.build_dual_graph();
	// 	mark_boundary(mesh);

	// 	using AuxRT_n= FunctionSpace< MeshT, RT0<1>>;
	// 	using AuxP_n= FunctionSpace< MeshT, Lagrange2<1>>;


	// 	AuxRT_n rtn(mesh);
	// 	AuxP_n pn(mesh);
	// 	auto Wtrial=MixedFunctionSpace(rtn,pn);
	// 	auto Waux=AuxFunctionSpacesBuild(pn);

	// 	auto W=FullSpaceBuild(Wtrial,Waux);



	// 	auto u0 = MakeTrial<0>(W);
	// 	auto u1 = MakeTrial<1>(W);

	// 	auto v0 = MakeTest<0>(W);
	// 	auto v1 = MakeTest<1>(W);

	// 	auto f1 = MakeFunction<0,ExactPoisson2D>(W);

	// 	std::cout<<W.spaces_ptr()->n_dofs()<<std::endl;

	//   // 2D LSFEM POISSION
	// 	auto bilinearform=L2Inner(Div(u0),Div(v0))+
	// 	L2Inner((u0),(v0))+
	// 	L2Inner(Grad(u1),Grad(v1))-
	// 	L2Inner(Grad(u1),(v0))-
	// 	L2Inner((u0),Grad(v1));
	// 	auto linearform=
	// 	L2Inner(-Div(v0),f1);




	// 	auto bcs1=DirichletBC<1,FunctionZero2D>(W,1);
	// 	auto bcs2=DirichletBC<1,FunctionZero2D>(W,2);
	// 	auto bcs3=DirichletBC<1,FunctionZero2D>(W,3);
	// 	auto bcs4=DirichletBC<1,FunctionZero2D>(W,4);

	// 	auto context=create_context(bilinearform,linearform,bcs1,bcs2,bcs3,bcs4);

	// 	SparseMatrix<Real> A;
	// 	std::vector<Real> b;


	// 	context.assembly(A,b);
	// 	context.apply_bc(A,b);



	// 	std::vector<Real> x;
	// 	Integer max_iter=25000;


	// 	std::cout<<"START SOLVING"<<std::endl;
	// 	gauss_seidel(x,A,b,max_iter);
	// 	std::cout<<"END SOLVING"<<std::endl;
	// 	std::string output_file1 ="LSFEM_Poisson2D_RT0_P1_output.vtk";
	// 	auto var_names=variables_names("stress","disp");
	// 	std::ofstream os;
	// 	os.open(output_file1.c_str());
	// 	write_wtk(os,W,x,var_names);
	// 	os.close();
	// }






	// void LSFEM_Poisson2D_RT1_P2()
	// {
	// 	constexpr Integer ManifoldDim=2;
	// 	constexpr Integer Dim=2;
	// 	using MeshT=Mesh<Dim, ManifoldDim>;

	// 	MeshT mesh;

	//   // generate_cube(mesh,4,4,4);
	//   // mesh.build_dual_graph();
	//   // mark_boundary(mesh);

	// 	generate_square(mesh,15,15);
	// 	mesh.build_dual_graph();
	// 	mark_boundary(mesh);

	//   // read_mesh("../data/triangle.MFEM", mesh);
	//   // assign_tags(mesh);
	//   // Bisection<MeshT> bisection(mesh);
	//   // bisection.uniform_refine(0);
	//   // mesh.update_dual_graph();


	// 	using AuxRT_n= FunctionSpace< MeshT, RT1<1>>;
	// 	using AuxP_n= FunctionSpace< MeshT, Lagrange2<1>>;
	// 	using AuxP_1= FunctionSpace< MeshT, Lagrange2<1>>;



	//  // using LSFEM_Poisson= FunctionSpace< MeshT, RT0<1>, Lagrange1<1>>;
	//  // LSFEM_Poisson LSFEM_poisson(mesh);




	// 	AuxRT_n rtn(mesh);
	// 	AuxP_n pn(mesh);
	// 	auto Wtrial=MixedFunctionSpace(rtn,pn);
	// 	auto Waux=AuxFunctionSpacesBuild(pn);
	// 	auto W=FullSpaceBuild(Wtrial,Waux);

	// 	auto u0 = MakeTrial<0>(W);
	// 	auto u1 = MakeTrial<1>(W);
	// 	auto v0 = MakeTest<0>(W);
	// 	auto v1 = MakeTest<1>(W);
	// 	auto f1 = MakeFunction<0,ExactPoisson2D>(W);

	//   // 2D LSFEM POISSION
	// 	auto bilinearform=L2Inner(Div(u0),Div(v0))+
	// 	L2Inner((u0),(v0))+
	// 	L2Inner(Grad(u1),Grad(v1))-
	// 	L2Inner(Grad(u1),(v0))-
	// 	L2Inner((u0),Grad(v1));
	// 	auto linearform=
	// 	L2Inner(-Div(v0),f1);




	// 	auto bcs1=DirichletBC<1,FunctionZero2D>(W,1);
	// 	auto bcs2=DirichletBC<1,FunctionZero2D>(W,2);
	// 	auto bcs3=DirichletBC<1,FunctionZero2D>(W,3);
	// 	auto bcs4=DirichletBC<1,FunctionZero2D>(W,4);

	// 	auto context=create_context(bilinearform,linearform,bcs1,bcs2,bcs3,bcs4);

	// 	SparseMatrix<Real> A;
	// 	std::vector<Real> b;


	// 	context.assembly(A,b);
	// 	context.apply_bc(A,b);


	// 	A.print();
	// 	std::vector<Real> x;
	// 	Integer max_iter=1000;


	// 	std::cout<<"START SOLVING"<<std::endl;
	// 	gauss_seidel(x,A,b,max_iter);
	// 	std::cout<<"END SOLVING"<<std::endl;
	// 	std::string output_file1 ="LSFEM_Poisson2D_RT1_P2_output1.vtk";
	// 	std::string output_file2 ="LSFEM_Poisson2D_RT1_P2_output1.vtk";
	// 	auto var_names=variables_names("stress","disp");

	// 	std::cout<<"---1"<<std::endl;
	// 	std::ofstream os;
	// 	os.open(output_file1.c_str());
	// 	std::cout<<"---2"<<std::endl;
	// 	write_wtk(os,W,x,var_names);
	// 	std::cout<<"---3"<<std::endl;

	// 	os.close();
	// 	std::cout<<"---4"<<std::endl;
	// 	os.open(output_file2.c_str());
	// 	std::cout<<"---5"<<std::endl;
	// 	write_wtk_isoparametric(os,W,x,var_names);
	// 	os.close();

	// // auto &dofsdofmap=W.dofsdofmap();

	// // for(Integer i=0;i<dofsdofmap.space_dofs().size();i++)
	// // {
	// //     std::cout<<"===== space="<<i<<std::endl;
	// //     for(Integer j=0;j<dofsdofmap.space_dofs_size(i);j++)
	// //     {
	// //       // for(Integer k=0;k<sd3[i]->operator[](j).size();k++)
	// //          std::cout<<dofsdofmap.space_dofs_get(i,j)<<std::endl;
	// //     }

	// // }

	// //     std::cout<<"-----dofmap of space="<<0<<std::endl;
	// // GetType<decltype(W)::DofsDM::ElemDofMap> elemdm;
	// // for(Integer i=0;i<dofsdofmap.template dofmap_size<0>();i++)
	// // {
	// //          dofsdofmap.template dofmap_get<0>(elemdm,i);
	// //          std::cout<<elemdm<<std::endl;
	// // }


	// }








	// void LSFEM_Poisson3D()
	// {
	// 	constexpr Integer ManifoldDim=3;
	// 	constexpr Integer Dim=3;
	// 	using MeshT=Mesh<Dim, ManifoldDim>;

	// 	MeshT mesh;

	// 	generate_cube(mesh,1,1,1);
	// 	mesh.build_dual_graph();
	// 	mark_boundary(mesh);

	//   // read_mesh("../data/triangle.MFEM", mesh);
	//   // read_mesh("../data/tetrahedron.MFEM", mesh);
	//   // assign_tags(mesh);
	//   // Bisection<MeshT> bisection(mesh);
	//   // bisection.uniform_refine(0);
	//   // mesh.update_dual_graph();


	// 	using AuxRT_n= FunctionSpace< MeshT, RT0<1>>;
	// 	using AuxP_n= FunctionSpace< MeshT, Lagrange2<1>>;
	// 	using AuxP_1= FunctionSpace< MeshT, Lagrange1<1>>;



	// 	using LSFEM_Poisson= FunctionSpace< MeshT, RT0<2>, Lagrange2<2>>;
	// 	LSFEM_Poisson LSFEM_poisson(mesh);




	// 	AuxRT_n rtn(mesh);
	// 	AuxP_n pn(mesh);
	// 	AuxP_1 p1(mesh);
	// 	auto Wtrial=MixedFunctionSpace(rtn,pn);
	// 	auto Waux=AuxFunctionSpacesBuild(p1,pn);

	// 	auto W=FullSpaceBuild(Wtrial,Waux);



	// 	auto u0 = MakeTrial<0>(W);
	// 	auto u1 = MakeTrial<1>(W);

	// 	auto v0 = MakeTest<0>(W);
	// 	auto v1 = MakeTest<1>(W);



	//  // auto sol0 = MakeFunction<0>(W);
	// 	std::cout<<"pre makefunct SOLVING"<<std::endl;
	// 	auto sol1 = MakeFunction<0>(W);
	// 	auto f1 = MakeFunction<1,Function1>(W);
	// 	std::cout<<"post makefunct SOLVING"<<std::endl;



	//  // auto f3 = MakeFunction<2>(W);
	// 	constexpr auto C=Constant<Scalar>(1.0);
	// 	constexpr auto Half=Constant<Scalar>(0.5);
	// 	constexpr auto alpha=Constant<Scalar>(2.0);
	// 	constexpr auto beta=Constant<Scalar>(3.0);
	// 	constexpr auto gamma=Constant<Scalar>(3.0);
	// 	constexpr auto id=Constant<Identity<Dim>>();

	// 	constexpr auto  alpha4=alpha.qp_eval<4>();

	// 	constexpr auto matrix1=Constant<Mat<1,2>>(1.,2.);
	// 	constexpr auto matrix2=Constant<Mat<2,1>>(1.,1.);


	//   // 2D LSFEM POISSION
	// 	auto bilinearform=L2Inner(Div(u0),Div(v0))+
	// 	L2Inner((u0),(v0))+
	// 	L2Inner(Grad(u1),Grad(v1))-
	// 	L2Inner(Grad(u1),(v0))-
	// 	L2Inner((u0),Grad(v1));
	// 	auto linearform=
	// 	L2Inner(-Div(v0),f1);




	// 	auto bcs1=DirichletBC<1,FunctionZero3D>(W,1);
	// 	auto bcs2=DirichletBC<1,FunctionZero3D>(W,2);
	// 	auto bcs3=DirichletBC<1,FunctionZero3D>(W,3);
	// 	auto bcs4=DirichletBC<1,FunctionZero3D>(W,4);
	// 	auto bcs5=DirichletBC<1,FunctionZero3D>(W,5);
	// 	auto bcs6=DirichletBC<1,FunctionZero3D>(W,6);

	// 	auto context=create_context(bilinearform,linearform,bcs1,bcs2,bcs3,bcs4,bcs5,bcs6);

	// 	SparseMatrix<Real> A;
	// 	std::vector<Real> b;


	// 	context.assembly(A,b);
	// 	context.apply_bc(A,b);

	// 	A.print();

	// 	std::vector<Real> x;
	// 	Integer max_iter=10000;


	// 	std::cout<<"START SOLVING"<<std::endl;
	// 	gauss_seidel(x,A,b,max_iter);
	// 	std::cout<<"END SOLVING"<<std::endl;
	// 	std::string output_file1 ="LSFEM_Poisson3D_output1.vtk";
	// 	std::string output_file2 ="LSFEM_Poisson3D_output2.vtk";
	// 	auto var_names=variables_names("stress","disp");

	// 	std::cout<<"---1"<<std::endl;
	// 	std::ofstream os;
	// 	os.open(output_file1.c_str());
	// 	std::cout<<"---2"<<std::endl;
	// 	write_wtk(os,W,x,var_names);
	// 	std::cout<<"---3"<<std::endl;

	// 	os.close();
	// 	std::cout<<"---4"<<std::endl;
	// 	os.open(output_file2.c_str());
	// 	std::cout<<"---5"<<std::endl;
	// 	write_wtk_isoparametric(os,W,x,var_names);
	// 	os.close();
	//    // std::cout<<"---6"<<std::endl;

	//    // auto ee= sol1.self3();
	//    // std::cout<<"---6.5"<<std::endl;
	//    // std::shared_ptr<std::vector<Real>> egm;
	//    // egm=std::make_shared<std::vector<Real>>(x);
	//    // std::cout<<"---6.6"<<std::endl;


	//    // auto zeroform=L2Inner((sol1),(C));
	//    // std::cout<<"---6.7"<<std::endl;




	//    // std::cout<<"---7"<<std::endl;





	//    // auto& left=zeroform.left();
	//    // auto iddd=(*left.global_dofs_ptr());
	//    // std::cout<<"---7.5"<<std::endl;
	// 	for(Integer i=0;i<sol1.aux_ptr()->vec_ptr()->size();i++)
	// 	{
	// 		std::cout<<(*(sol1.aux_ptr()->vec_ptr()))[i]<<std::endl;
	// 	}
	//    // std::cout<<"---8"<<std::endl;
	//    // auto context0=create_context(zeroform);
	//    // std::cout<<"---9"<<std::endl;
	//    // sol1.global_dofs_update(x);
	//    // Real int_x=context0.assembly();
	//    // std::cout<<"INTEGRAL WITH X="<<int_x<<std::endl;
	//    // std::cout<<"---10"<<std::endl;

	//    //  std::vector<Real> y(100,1.0);
	//    //  sol1.global_dofs_update(y);
	//    // Real int_y=context0.assembly();
	//    // std::cout<<"INTEGRAL WITH Y="<<int_y<<std::endl;
	//    // std::cout<<"---11"<<std::endl;


	//    // auto zeroform2=L2Inner(C,C);
	//    // auto context1=create_context(zeroform2);
	// //   for(int i=0;i<LSFEM_poisson.space_dofs(0,0).size();i++)
	// //   std::cout<< LSFEM_poisson.space_dofs(0,0)[i]<<std::endl;
	// // for(int i=0;i<LSFEM_poisson.space_dofs(1,0).size();i++)
	// //   std::cout<< LSFEM_poisson.space_dofs(1,0)[i]<<std::endl;

	// // std::cout<< LSFEM_poisson.n_dofs(0,0)<<std::endl;
	// // std::cout<< LSFEM_poisson.n_dofs(1,0)<<std::endl;

	// // std::cout<< LSFEM_poisson.n_dofs(0)<<std::endl;
	// // std::cout<< LSFEM_poisson.n_dofs(1)<<std::endl;
	// 	std::cout<< W.array_ndofs()<<std::endl;
	// // auto &dm=W.space_dofs();
	// // for(Integer i=0;i<dm.size();i++)
	// // {
	// //     std::cout<<"SPACE="<<i<<std::endl;
	// //     for(Integer j=0;j<dm[i].size();j++)
	// //     {
	// //       for(Integer k=0;k<dm[i][j].size();k++)
	// //          std::cout<<dm[i][j][k]<<std::endl;
	// //     }

	// // }
	// // std::cout<<"space_info="<<LSFEM_poisson.space_info()<<std::endl;



	// 	auto &dofsdofmap=W.dofsdofmap();

	//  // auto &dm3=dofsdofmap.dofmap();
	//  // auto& dm30=tuple_get<0>(dm3);
	//  // auto& dm31=tuple_get<1>(dm3);
	// // const auto& sd3=dofsdofmap.space_dofs();


	// 	for(Integer i=0;i<dofsdofmap.space_dofs().size();i++)
	// 	{
	// 		std::cout<<"===== space="<<i<<std::endl;
	// 		for(Integer j=0;j<dofsdofmap.space_dofs_size(i);j++)
	// 		{
	//       // for(Integer k=0;k<sd3[i]->operator[](j).size();k++)
	// 			std::cout<<dofsdofmap.space_dofs_get(i,j)<<std::endl;
	// 		}

	// 	}

	// 	std::cout<<"-----dofmap of space="<<0<<std::endl;

	// // for(Integer i=0;i<dofsdofmap.template dofmap_size<0>();i++)
	// // {
	// //          std::cout<<dofsdofmap.template dofmap_get<0>(i)<<" ";
	// //          std::cout<<std::endl;
	// // }
	// //     std::cout<<"-----dofmap of space="<<1<<std::endl;

	// // for(Integer i=0;i<dofsdofmap.template dofmap_size<1>();i++)
	// // {
	// //          std::cout<<dofsdofmap.template dofmap_get<1>(i)<<" ";
	// //          std::cout<<std::endl;
	// // }

	// //     std::cout<<"-----dofmap of space="<<2<<std::endl;

	// // for(Integer i=0;i<dofsdofmap.template dofmap_size<2>();i++)
	// // {
	// //          std::cout<<dofsdofmap.template dofmap_get<2>(i)<<" ";
	// //          std::cout<<std::endl;
	// // }
	// //     std::cout<<"-----dofmap of space="<<3<<std::endl;

	// // for(Integer i=0;i<dofsdofmap.template dofmap_size<3>();i++)
	// // {
	// //          std::cout<<dofsdofmap.template dofmap_get<3>(i)<<" ";
	// //          std::cout<<std::endl;
	// // }


	// 	std::cout<<"-----dofmap of space="<<dofsdofmap.n_dofs()<<std::endl;
	// 	std::cout<<"-----dofmap of space="<<dofsdofmap.cumulative_n_dofs()<<std::endl;


	// 	std::cout<<"===== solution x="<<std::endl;
	// 	for(Integer i=0;i<x.size();i++)
	// 	{

	//       // for(Integer k=0;k<sd3[i]->operator[](j).size();k++)
	// 		std::cout<<x[i]<<std::endl;

	// 	}
	//   // for(int el=0;el<mesh.n_elements();el++)
	//   // {
	//   //   std::cout<<"elem="<<el<<std::endl;
	//   //   std::cout<<"is active="<<mesh.is_active(el)<<std::endl;

	//   //   const auto& nodes=mesh.elem(el).nodes;
	//   //   auto side_tags=mesh.elem(el).side_tags;
	//   //     for(int i=0;i<side_tags.size();i++)
	//   //     std::cout<<mesh.points()[nodes[i]]<<" ";
	//   //   std::cout<<std::endl;
	//   //     for(int i=0;i<side_tags.size();i++)
	//   //     std::cout<<side_tags[i]<<" "<<std::endl;
	//   // std::cout<<std::endl;
	//   // }

	//   // Vector<Matrix<Real, 3, 1>,10> out;
	//   // Vector<Real,3> p;
	//   // const auto& points=QuadraticTetrahedronPoints;
	//   // for(Integer i=0;i<QuadraticTetrahedronPoints.size();i++)
	//   // {
	//   //   for(Integer j=0;j<3;j++)
	//   //     p[j]=points[i](j,0);
	//   // ReferenceShapeFunctionValue<Simplex<3,3>,GradientOperator, LagrangeFE, 2>::apply(p,out);
	//   // std::cout<<"p"<<std::endl;
	//   // std::cout<<p<<std::endl;
	//   // std::cout<<"out"<<std::endl;
	//   // std::cout<<out<<std::endl;
	//   // }

	// }



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



	// 	MeshT mesh;
	// 	using Elem = typename MeshT::Elem; 
	//   // read_mesh("../data/beam-tri.MFEM", mesh);
	//   // read_mesh("../data/triangle_square.MFEM", mesh);
	//   // read_mesh("../data/triangle.MFEM", mesh);

	// 	generate_square(mesh, 2, 2 );

	// 	mesh.build_dual_graph();
	//   // mesh.update_dual_graph();
	// 	mark_boundary(mesh);
	//   // assign_tags(mesh);
	//   // Bisection<MeshT> bisection(mesh);
	//   // bisection.uniform_refine(8);
	//   // mesh.update_dual_graph();

	// 	constexpr Integer QPOrder=4;
	// 	constexpr Integer NQPoints=GaussPoints<Elem,QPOrder>::NQPoints; 
	// 	using QP=typename GaussPoints<Elem,QPOrder>::qp_points_type;

	// 	constexpr Integer Npoints=Elem::Npoints;
	// 	std::vector<Vector<Real,Dim>> points(Npoints);



	//  using P1= FunctionSpace< MeshT, Lagrange1<1>>;//,Lagrange1<2>>;
	//  P1 p1(mesh);
	//  using P2= FunctionSpace< MeshT, Lagrange2<1>>;//,Lagrange1<2>>;
	//  P2 p2(mesh);

	//  using RT_0= FunctionSpace< MeshT, RT0<1>>;//,Lagrange1<2>>;
	//  RT_0 rt(mesh);

	//  using LSFEM_Poisson= FunctionSpace< MeshT, RT0<1>, Lagrange2<1>>;
	//  LSFEM_Poisson LSFEM_poisson(mesh);

	//  // using FSspace2= FunctionSpace< MeshT, Lagrange2<2>>;
	//  // FSspace2 FEspace2(mesh);
	//  using LSFEM_space= FunctionSpace< MeshT, RT0<Dim>, Lagrange2<Dim>>;
	//  LSFEM_space LS_space(mesh);
	//  using FSspace4= FunctionSpace< MeshT, RT0<1>>;
	//  FSspace4 FEspace4(mesh);
	//  using P1_2D= FunctionSpace< MeshT, Lagrange1<2>>;
	//  P1_2D p1_2D(mesh);
	//  using FSspace6= FunctionSpace< MeshT, Lagrange3<2>>;
	//  FSspace6 FEspace6(mesh);
	//  using FSspace7= FunctionSpace< MeshT, Lagrange1<1>,Lagrange2<2>,Lagrange3<3>,Lagrange1<2>,Lagrange2<2> >;
	//  FSspace7 FEspace7(mesh);
	//  using FSspace8= FunctionSpace< MeshT, RT0<1>,RT1<2>>;
	//  FSspace8 FEspace8(mesh);
	//  using FSspace9= FunctionSpace< MeshT, RT0<2> >;
	//  FSspace9 FEspace9(mesh);
	//  using FSspace10= FunctionSpace< MeshT, RT1<2>>;
	//  FSspace10 FEspace10(mesh);

	//  // auto W5=MixedFunctionSpace(MixedFunctionSpace(FEspace7,FEspace8),FEspace1);

	//  using AuxFSspace1= FunctionSpace< MeshT, Lagrange1<1> >;
	//  using AuxFSspace2= FunctionSpace< MeshT, Lagrange2<1> >;
	//  using AuxFSspace3= FunctionSpace< MeshT, Lagrange1<2> >;
	//  AuxFSspace1 AuxFEspace1(mesh);
	//  AuxFSspace2 AuxFEspace2(mesh);
	//  AuxFSspace3 AuxFEspace3(mesh);


	//  // auto Wtrial=MixedFunctionSpace(FEspace4,FEspace1);
	//  // auto Wtrial=MixedFunctionSpace(FEspace1);
	//  // auto Wtrial=MixedFunctionSpace(LS_space);
	//  auto Wtrial=MixedFunctionSpace(LSFEM_poisson);

	//  // auto Wtrial=MixedFunctionSpace(p1_2D);
	//  // auto Wtrial=MixedFunctionSpace(p1);
	//  // auto Wtrial=MixedFunctionSpace(p2);

	//  auto Waux=AuxFunctionSpacesBuild(AuxFEspace1);//,AuxFEspace2,AuxFEspace3);
	//  // auto W=FullSpaceBuild(Wtrial,Waux);
	//  auto W=FullSpaceBuild(Wtrial);



	//  // auto f1 = MakeFunction<0,Function1>(W);
	//  // auto f2 = MakeFunction<1>(W);
	//  // auto f3 = MakeFunction<2>(W);

	//  auto u0 = MakeTrial<0>(W);
	//  auto u1 = MakeTrial<1>(W);
	//  // auto u2 = MakeTrial<2>(W);

	//  auto v0 = MakeTest<0>(W);
	//  auto v1 = MakeTest<1>(W);
	//  // auto v2 = MakeTest<2>(W);

	//  FiniteElem<Elem> FE(mesh);

	//  constexpr auto C=Constant<Scalar>(1.0);
	//  constexpr auto Half=Constant<Scalar>(0.5);
	//  constexpr auto alpha=Constant<Scalar>(2.0);
	//  constexpr auto beta=Constant<Scalar>(3.0);
	//  constexpr auto gamma=Constant<Scalar>(3.0);
	//  constexpr auto id=Constant<Identity<Dim>>();

	//  constexpr auto  alpha4=alpha.qp_eval<4>();

	//  constexpr auto matrix1=Constant<Mat<1,2>>(1.,2.);
	//  constexpr auto matrix2=Constant<Mat<2,1>>(1.,1.);


	//  auto NewOp1=NewOperator(IdentityOperator()/alpha);
	// // auto NewOp2=NewOperator(IdentityOperator()*alpha*f1);

	// // auto Epsilon=NewOperator((-f1)*(+C)*((+GradientOperator())+(+Transpose(GradientOperator()))));
	//  auto Eps=NewOperator(Half*(GradientOperator()+(Transpose(GradientOperator()))));
	//  auto C_inv=NewOperator(alpha * IdentityOperator() + beta * MatTrace(IdentityOperator()) * id );
	// // auto Epsilon=NewOperator(+(-GradientOperator()));
	// // auto Epsilon=NewOperator((+C)*(-C)*(-f1)*(+f1)*(C-f1)*(C+f1)/(C+f1)*Transpose(f1)*(+MatTrace(+f1))*(-MatTrace(-C))*(-GradientOperator()-GradientOperator())/MatTrace(Transpose(f1)-MatTrace(Transpose(C))));
	//  auto NewTrace=NewOperator(TraceOperator());

	// // /MatTrace(Transpose(Transpose(MatTrace(f1)))));
	// // auto Epsilon=NewOperator((GradientOperator()+Transpose(GradientOperator())));//+Transpose(GradientOperator())));
	// // auto Epsilon=NewOperator((GradientOperator()));//+Transpose(GradientOperator())));

	//   // auto bilinearform=
	//                     // surface_integral(0,NewTrace(u2),NewTrace(v2))-
	//                     // surface_integral(1,NewTrace(u1),NewTrace(v1))-
	//                     // surface_integral(0,NewTrace(u0),NewTrace(v0))-


	//                     // surface_integral(NewTrace(u2),NewTrace(v2))-
	//                     // surface_integral(Trace(u2),Trace(v2))-
	//                     // L2Inner(-(-u2),-(v2))-
	//                     // surface_integral(Trace(u2),Trace(v2))-
	//                     // L2Inner(-(-u2),-(v2))-
	//                     // surface_integral(Trace(u2),Trace(v2))-
	//                     // surface_integral(Trace(u2),Trace(v2))-
	//                     // L2Inner(-(-u2),-(v2))-
	//                     // surface_integral(Trace(u2),Trace(v2))-
	//                     // surface_integral(Trace(u2),Trace(v2))-
	//                     // L2Inner(Transpose(u2),Transpose(v2))-
	//                     // L2Inner(-Transpose(u2),Transpose(-v2))-
	//                     // L2Inner(-(-u2),-(v2))- //+ L2Inner(Grad(u2),Grad(v0))+L2Inner(u2,v2)+
	//                     // L2Inner(MatTrace(f1)*(+Transpose(u1)),-(Transpose(v1))) -//+ L2Inner(Grad(u1),Grad(v0))+L2Inner(u1,v2)+
	//                     // surface_integral(Trace(u2),Trace(v2))-
	//                     // L2Inner(-(-u2),-(v2))-
	//                     // surface_integral(Trace(u2),Trace(v2))-
	//                     // L2Inner(-(-u2),-(v2))-
	//                     // surface_integral(Trace(u2),Trace(v2))-
	//                     // L2Inner(Epsilon(u0),-Grad(v0));//+ L2Inner(Grad(u0),Grad(v0))+L2Inner(u0,v2);//+L2Inner(Grad(u1),Grad(v1));//+L2Inner(f3*u2,v2);
	//   // auto bilinearform= 
	//                     // L2Inner(u0,v1)+L2Inner(u0,v2)+
	//                     // L2Inner(u0+u1,v0+v1)-
	//                     // L2Inner(u0,v0)+L2Inner(u0,v1)+ //L2Inner(u0,v2)+
	//                     // L2Inner(u1,v0)+L2Inner(u1,v1)+//L2Inner(u1,v2)+
	//                     // L2Inner(u2,v0)+L2Inner(u2,v1)+
	//                     // L2Inner(Grad(u0),Grad(v0))

	//                     // L2Inner(Div(u0),Div(v0))+
	//                     // L2Inner((u0),(v0))+
	//                     // L2Inner(Grad(u1),Grad(v1))
	//                     // -L2Inner(Grad(u1),v0)
	//                     // -L2Inner(u0,Grad(v1))

	//                     // +L2Inner((u2),(v2))
	//                     //     +
	//                     // L2Inner(Grad(u0),Grad(v0))+
	//                     // L2Inner(Grad(u1),Grad(v1))//+
	//                     // L2Inner(Grad(u2),Grad(v2))
	//                     // ;

	//   // 2D LSFEM POISSION
	//  auto bilinearform=L2Inner(Div(u0),Div(v0))+
	//  L2Inner((u0),(v0))+
	//  L2Inner(Grad(u1),Grad(v1))-
	//  L2Inner(Grad(u1),(v0))-
	//  L2Inner((u0),Grad(v1));
	//  auto linearform=
	//  L2Inner(-Div(v0),C);

	//   // 2D LSFEM ELASTICITY
	//   // auto bilinearform=L2Inner(Div(u0),Div(v0))+
	//   //                   L2Inner(C_inv(u0),C_inv(v0))+
	//   //                   L2Inner(Eps(u1),Eps(v1))-
	//   //                   L2Inner(Eps(u1),C_inv(v0))-
	//   //                   L2Inner(C_inv(u0),Eps(v1));
	//   //  auto linearform=L2Inner(-Div(v0),matrix2);

	//    // 2D POISSON WITH NEUMANN
	//    // auto bilinearform=L2Inner(Grad(u0),Grad(v0));
	//    // auto linearform=L2Inner((v0),matrix2)+surface_integral(Trace(v0),matrix2);

	//    // auto bilinearform=L2Inner(Grad(u0),Grad(v0));
	//    // auto linearform=
	//    //                    surface_integral(1,Trace(v0),C)+
	//    //                    L2Inner((v0),C);

	//   // auto linearform=
	//                   //L2Inner(Grad(v0),+Transpose(id2));//Transpose(Transpose((matrix1)+Transpose(matrix2))));//Transpose(-f1*(matrix1+matrix1)*Transpose(alpha*matrix1-matrix1)));//+L2Inner(f2,v1)+L2Inner(f1,v0);//+ L2Inner(f1,v0);
	//                   // L2Inner((+Transpose(C))*(-v1),-Transpose(f1))+
	//                   // L2Inner((Transpose(C)+Transpose(C))*v1,C)+
	//                   // L2Inner((Transpose(C)+(C))*v1,f1)+
	//                   // L2Inner((C+C)*v1,C)+
	//                   // L2Inner((-C-Transpose(C))*(-v1),Transpose(-f1))+
	//                   // L2Inner((-Transpose(C)-Transpose(C))*v1,C)+
	//                   // L2Inner((-Transpose(C)-(C))*v1,f1/C)+
	//                   // L2Inner(Inner(gamma,Transpose(gamma)),MatTrace(Epsilon(v0)));
	//   // -L2Inner(Inner(gamma,Transpose(gamma))*Epsilon(v0),id2)-L2Inner(-Epsilon(v0),-id2)+
	//   // L2Inner(id2,v2)
	//   // L2Inner(id2,v2)-
	//   // L2Inner(v2,C)+
	//   // L2Inner(v1,C)+
	//   // L2Inner(v0,C)
	//   // surface_integral(1,NewTrace(v0),C)
	//   // +
	//   // surface_integral(20,NewTrace(v0),C)-
	//   // surface_integral(1,NewTrace(v0),C)
	//   // +L2Inner((v0),C)

	//   // L2Inner(-Div(v0),C)

	//   // L2Inner((v0),C)
	//   // ;


	//  auto bilinear_form=general_form(bilinearform);
	//  auto linear_form=general_form(linearform);



	//  auto shape_coefficients=shape_function_coefficients(bilinear_form,linear_form);
	//  auto reference_maps=reference_maps2(bilinear_form,linear_form);
	//  auto shapefunctions=shape_functions(shape_coefficients,reference_maps,bilinear_form,linear_form);



	//   // auto spaces_ptr=bilinear_form.spaces_ptr()->spaces_ptr();
	//   // auto n_dofs=spaces_ptr->n_dofs();

	//   // std::vector<std::vector<Real>> A;
	//  SparseMatrix<Real> A;
	//  std::vector<Real> b;


	//   // A.resize(n_dofs,std::vector<Real>(n_dofs));
	//   // b.resize(n_dofs);

	//  auto bcs1=DirichletBC<1,FunctionZero2D>(W,1);
	//  auto bcs2=DirichletBC<1,FunctionZero2D>(W,2);
	//  auto bcs3=DirichletBC<1,FunctionZero2D>(W,3);
	//  auto bcs4=DirichletBC<0,Function4>(W,4);

	//    // auto bcs1=DirichletBC<0,Function2>(W,1);
	//    // auto bcs2=DirichletBC<0,Function2>(W,2);
	//    // auto bcs3=DirichletBC<0,Function2>(W,3);
	//    // auto bcs4=DirichletBC<0,Function2>(W,4);

	//    // auto bcs1=DirichletBC<0,FunctionZero1D>(W,1);
	//    // auto bcs2=DirichletBC<0,FunctionZero1D>(W,2);
	//    // auto bcs3=DirichletBC<0,FunctionZero1D>(W,3);
	//    // auto bcs4=DirichletBC<0,FunctionZero1D>(W,4);


	//    // DirichletBCMapsCollection<decltype(bcs0)> kk(5);
	//    // DirichletBCMapsCollection<decltype(bcs0),decltype(bcs1),decltype(bcs2)> k3k(5);


	//    // auto context=create_context(bilinearform,linearform,bcs2,bcs3,bcs4);
	//  auto context=create_context(bilinearform,linearform,bcs1,bcs2,bcs3,bcs4);


	//    // for(std::size_t i=0;i<mesh.boundary2elem().size();i++)
	//    //  std::cout<<mesh.boundary2elem()[i]<<std::endl;

	//    // auto bcs_maps=DirichletMapCollection(bcs0,bcs1,bcs2,bcs3);


	//    // FE.init(0);
	//    // FE.init_boundary(0);
	//    // std::cout<<"__________elem_id="<<FE.elem_id()<<" side_tag="<<FE.side_tag()<<std::endl;
	//    // bcs_maps.init(FE);
	//    // FE.init_boundary(1);
	//    // std::cout<<"__________elem_id="<<FE.elem_id()<<" side_tag="<<FE.side_tag()<<std::endl;
	//    // bcs_maps.init(FE);
	//    // FE.init_boundary(2);
	//    // std::cout<<"__________elem_id="<<FE.elem_id()<<" side_tag="<<FE.side_tag()<<std::endl;
	//    // bcs_maps.init(FE);

	//    // FE.init(1);
	//    // FE.init_boundary(0);
	//    // std::cout<<"__________elem_id="<<FE.elem_id()<<" side_tag="<<FE.side_tag()<<std::endl;
	//    // bcs_maps.init(FE);
	//    // FE.init_boundary(1);
	//    // std::cout<<"__________elem_id="<<FE.elem_id()<<" side_tag="<<FE.side_tag()<<std::endl;
	//    // bcs_maps.init(FE);
	//    // FE.init_boundary(2);
	//    // std::cout<<"__________elem_id="<<FE.elem_id()<<" side_tag="<<FE.side_tag()<<std::endl;
	//    // bcs_maps.init(FE);

	//    // Array<std::map<int,int>,2> mmm;
	//    // mmm[0].insert(std::pair<int, int>(1, 0));
	//    // mmm[0].insert(std::pair<int, int>(2, 0));
	//    // mmm[0].insert(std::pair<int, int>(3, 0));

	//    // std::map<int, int>::iterator itr; 
	//    // for (itr = mmm[0].begin(); itr != mmm[0].end(); ++itr) { 
	//    //      cout << '\t' << itr->first 
	//    //           << '\t' << itr->second << '\n'; 
	//    //  }

	//    // // bcs.assembly(FE);
	//    // auto tr=TraceDofs<decltype(W)>::dofs();
	//    // std::cout<<tuple_get<0>(tr)<<std::endl;
	//    // std::cout<<tuple_get<1>(tr)<<std::endl;
	//    // std::cout<<tuple_get<2>(tr)<<std::endl;


	//    // auto dofmap=W.dofmap2();
	//    //     auto dm1=tuple_get<0>(dofmap);
	//    //     auto dm2=tuple_get<1>(dofmap);
	//    //     auto dm3=tuple_get<2>(dofmap);
	//    //     std::cout<<std::endl;
	//    //    std::cout<<"___--dm1="<<std::endl;
	//    //  for(std::size_t i=0;i<dm1.size();i++)
	//    //      {
	//    //       std::cout<<dm1[i]<<std::endl;
	//    //       }     
	//    //    std::cout<<"___--dm2="<<std::endl;
	//    //  for(std::size_t i=0;i<dm2.size();i++)
	//    //      {
	//    //       std::cout<<dm2[i]<<std::endl;
	//    //       }   
	//    //    std::cout<<"___--dm3="<<std::endl;
	//    //  for(std::size_t i=0;i<dm3.size();i++)
	//    //      {
	//    //       std::cout<<dm3[i]<<std::endl;
	//    //       }  
	//     // std::cout<<std::endl;___--dofmap_trace=
	// //    Simplex<6,6> elem6;
	// //    Simplex<6,5> elem5;

	// //             for(Integer i = 0; i < n_nodes(elem6); ++i) {
	// //                 elem6.nodes[i] = i;
	// //             }

	// //    elem6.side(0,elem5); 
	// //    for(std::size_t i=0;i<elem5.nodes.size();i++)
	// //    std::cout<<elem5.nodes[i]<<std::endl;
	// //    std::cout<<std::endl;
	// //    elem6.side(1,elem5); 
	// //    for(std::size_t i=0;i<elem5.nodes.size();i++)
	// //    std::cout<<elem5.nodes[i]<<std::endl;
	// // std::cout<<std::endl;
	// //    elem6.side(2,elem5); 
	// //    for(std::size_t i=0;i<elem5.nodes.size();i++)
	// //    std::cout<<elem5.nodes[i]<<std::endl;
	// // std::cout<<std::endl;
	// //    elem6.side(3,elem5); 
	// //    for(std::size_t i=0;i<elem5.nodes.size();i++)
	// //    std::cout<<elem5.nodes[i]<<std::endl;
	// // std::cout<<std::endl;
	// //    elem6.side(4,elem5); 
	// //    for(std::size_t i=0;i<elem5.nodes.size();i++)
	// //    std::cout<<elem5.nodes[i]<<std::endl;
	// // std::cout<<std::endl;
	// //    elem6.side(5,elem5); 
	// //    for(std::size_t i=0;i<elem5.nodes.size();i++)
	// //    std::cout<<elem5.nodes[i]<<std::endl;
	// // std::cout<<std::endl;

	//   // auto dofmap=W.dofmap2();
	//   //      auto dm1=tuple_get<0>(dofmap);
	//   //      auto dm2=tuple_get<1>(dofmap);
	//   //      auto dm3=tuple_get<2>(dofmap);
	//   //      std::cout<<std::endl;
	//   //     std::cout<<"___--dm1="<<std::endl;
	//   //   for(std::size_t i=0;i<dm1.size();i++)
	//   //       {
	//   //        std::cout<<dm1[i]<<std::endl;
	//   //        }     
	//   //     std::cout<<"___--dm2="<<std::endl;
	//   //   for(std::size_t i=0;i<dm2.size();i++)
	//   //       {
	//   //        std::cout<<dm2[i]<<std::endl;
	//   //        }   
	//   //     std::cout<<"___--dm3="<<std::endl;
	//   //   for(std::size_t i=0;i<dm3.size();i++)
	//   //       {
	//   //        std::cout<<dm3[i]<<std::endl;
	//   //        }  
	//   //   std::cout<<std::endl;
	//    // decltype(W) kl(5);

	//    ////////////////////////////////////////////////////////////////////////////////////////////////////
	//    ////////////////////////////////////////////////////////////////////////////////////////////////////
	//    ////////////////////////////////////////////////////////////////////////////////////////////////////
	//    /// RICORDA:
	//    /// NON HA SENSO FARE SOMME DI ZERI: QUINDI CONSIDERA BENE GLI INTEGRALI DI SUPERFICIE SE HAI SOMME...
	//    /// MA SOPRATUTTO SE NON NE HAI
	//    /// IDEALE: DEVI DIRE SE UN ELEMENTO HA SIDES DI BORDO

	//    // RICORDA CHE SE HAI UN SIDE CON TAG CHE NON E' PRESENTE NEGLI L2PRODUCT (AD ES.DIRICHLET )
	//    // E' INUTILE CALCOLARE LE SHAPE FUNCTIONS
	//    // QUINDI LE FORMS DOVREBBERO RESTITUIRE I TAGS DEI LORO INTEGRALI...

	//     // RICORDA CHE GLI INTEGRALI DI BORDO VANNO MOLTIPLICATI PER DET(J) BORDO
	//     // MODIFICA EVAK GENERAL UTILS LocalTensor
	//    ////////////////////////////////////////////////////////////////////////////////////////////////////
	//    ////////////////////////////////////////////////////////////////////////////////////////////////////
	//    ////////////////////////////////////////////////////////////////////////////////////////////////////

	//    // context.assembly(A,b);
	//  auto n_dofs=W.spaces_ptr()->n_dofs();
	//  std::cout<<"PRE APPLY BC n_dofs()="<<std::endl;
	//  std::cout<<n_dofs<<std::endl;
	//   // for(Integer ii=0;ii<n_dofs;ii++)
	//   //   {
	//   //     for(Integer jj=0;jj<n_dofs;jj++)
	//   //       // std::cout<<A[ii][jj]<<" ";
	//   //       std::cout<<A(ii,jj)<<" ";
	//   //     std::cout<<std::endl;
	//   //   }
	//  A.print();

	//  for(Integer ii=0;ii<n_dofs;ii++)
	//  {
	//  	std::cout<<b[ii]<<std::endl;
	//  }

	//  context.apply_bc(A,b);



	//  std::vector<Real> x(n_dofs);
	//  Integer max_iter=10000;



	//  std::cout<<"AFTER APPLY BC A="<<std::endl;
	//   // for(Integer ii=0;ii<n_dofs;ii++)
	//   //   {
	//   //     for(Integer jj=0;jj<n_dofs;jj++)
	//   //       // std::cout<<A[ii][jj]<<" ";
	//   //       std::cout<<A(ii,jj)<<" ";
	//   //     std::cout<<std::endl;
	//   //   }
	//  A.print();
	//  std::cout<<"AFTER APPLY BC b="<<std::endl;
	//  for(Integer ii=0;ii<n_dofs;ii++)
	//  {
	//  	std::cout<<b[ii]<<std::endl;
	//  }

	//  gauss_seidel(x,A,b,max_iter);
	//  std::cout<<"solution with n dofs="<<n_dofs<<std::endl;
	//  for(Integer ii=0;ii<x.size();ii++)
	//  {
	//  	std::cout<<x[ii]<<std::endl;
	//  }

	//  std::string output_file1 ="output1.vtk";
	//  std::string output_file2 ="output2.vtk";
	//    // auto variables_names=variables_names("velocity","pressure","pressure2");
	//  auto var_names=variables_names("stress","disp");

	//    // auto variables_names=variables_names("pressure");//,"pressure","pressure2");
	//  std::ofstream os;
	//  os.open(output_file1.c_str());
	//  write_wtk(os,W,x,var_names);

	//  os.close();
	//  os.open(output_file2.c_str());
	//  write_wtk_isoparametric(os,W,x,var_names);
	//  os.close();



	//  std::cout<<"dim = 2, order=0"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<2,0>()<<std::endl;
	//  std::cout<<"dim = 2, order=1"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<2,1>()<<std::endl;
	//  std::cout<<"dim = 2, order=2"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<2,2>()<<std::endl;
	//  std::cout<<"dim = 2, order=3"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<2,3>()<<std::endl;
	//  std::cout<<"dim = 2, order=4"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<2,4>()<<std::endl;


	//  std::cout<<"dim = 3, order=0"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<3,0>()<<std::endl;
	//  std::cout<<"dim = 3, order=1"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<3,1>()<<std::endl;
	//  std::cout<<"dim = 3, order=2"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<3,2>()<<std::endl;
	//  std::cout<<"dim = 3, order=3"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<3,3>()<<std::endl;
	//  std::cout<<"dim = 3, order=4"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<3,4>()<<std::endl;

	//  std::cout<<"dim = 4, order=0"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<4,0>()<<std::endl;
	//  std::cout<<"dim = 4, order=1"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<4,1>()<<std::endl;
	//  std::cout<<"dim = 4, order=2"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<4,2>()<<std::endl;
	//  std::cout<<"dim = 4, order=3"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<4,3>()<<std::endl;
	//  std::cout<<"dim = 4, order=4"<<std::endl;
	//  std::cout<<NumberOfLagrangianSimplexDofs<4,4>()<<std::endl;
	//   // Variable<Elem,IdentityOperator,LagrangeFE,1,1,2>::Points::type eee(5);
	//   // decltype(Variable<Elem,IdentityOperator,LagrangeFE,1,1,2>::Points::points) okok(6);
	 
	//   // std::cout<<"----------------"<<std::endl;
	//   // std::cout<<Variable<Elem,IdentityOperator,RaviartThomasFE,0,1,2>::Points::points<<std::endl;
	//   //  std::cout<<"----------------"<<std::endl;
	//   // std::cout<<Variable<Elem,IdentityOperator,RaviartThomasFE,0,1,2>::reference_values<<std::endl;
	//   // std::cout<<"----------------"<<std::endl;
	//   // Array<Real,3> sol_tmp{1.0,1.0,1.0};
	//   // Variable<Elem,IdentityOperator,RaviartThomasFE,0,1,2> var;
	//   // FE.init(0);
	//   // var.init(FE,sol_tmp);
	//   // std::cout<<"sol_tmp="<<std::endl;
	//   // std::cout<<sol_tmp<<std::endl;

	//   // std::cout<<var.value()<<std::endl;
	//   // for(int el=0;el<mesh.side_nodes().size();el++)
	//   // {
	//   //   std::cout<<"side_nodes_="<<el<<std::endl;
	//   //   auto ecm=mesh.side_nodes()[el];
	//   //     for(int i=0;i<ecm.size();i++)
	//   //     std::cout<<ecm[i]<<" "<<std::endl;
	//   // std::cout<<std::endl;
	//   // }
	//   // for(int el=0;el<mesh.side_tags().size();el++)
	//   // {
	//   //   std::cout<<"side_tags="<<el<<std::endl;
	//   //   auto ecm=mesh.side_tags()[el];
	//   //     std::cout<<ecm<<" "<<std::endl;
	//   // std::cout<<std::endl;
	//   // }

	//  for(int el=0;el<mesh.n_elements();el++)
	//  {
	//  	std::cout<<"elem="<<el<<std::endl;
	//  	std::cout<<"is active="<<mesh.is_active(el)<<std::endl;

	//  	const auto& nodes=mesh.elem(el).nodes;
	//  	auto side_tags=mesh.elem(el).side_tags;
	//  	for(int i=0;i<side_tags.size();i++)
	//  		std::cout<<mesh.points()[nodes[i]]<<" ";
	//  	std::cout<<std::endl;
	//  	for(int i=0;i<side_tags.size();i++)
	//  		std::cout<<side_tags[i]<<" "<<std::endl;
	//  	std::cout<<std::endl;
	//  }

	  //  MeshT meshprova;
	  //   generate_square(meshprova,2,2);

	  // for(int el=0;el<meshprova.n_elements();el++)
	  // {
	  //   std::cout<<"elem="<<el<<std::endl;
	  //   std::cout<<"is active="<<meshprova.is_active(el)<<std::endl;

	  //   const auto& nodes=meshprova.elem(el).nodes;
	  //   auto side_tags=meshprova.elem(el).side_tags;
	  //     for(int i=0;i<side_tags.size();i++)
	  //     std::cout<<meshprova.points()[nodes[i]]<<" ";
	  //   std::cout<<std::endl;
	  //     for(int i=0;i<side_tags.size();i++)
	  //     std::cout<<side_tags[i]<<" "<<std::endl;
	  // std::cout<<std::endl;
	  // }
	  // meshprova.build_dual_graph();
	  // // mesh.update_dual_graph();
	  // mark_boundary(meshprova);
	  // for(int el=0;el<meshprova.n_elements();el++)
	  // {
	  //   std::cout<<"elem="<<el<<std::endl;
	  //   std::cout<<"is active="<<meshprova.is_active(el)<<std::endl;

	  //   const auto& nodes=meshprova.elem(el).nodes;
	  //   auto side_tags=meshprova.elem(el).side_tags;
	  //     for(int i=0;i<side_tags.size();i++)
	  //     std::cout<<meshprova.points()[nodes[i]]<<" ";
	  //   std::cout<<std::endl;
	  //     for(int i=0;i<side_tags.size();i++)
	  //     std::cout<<side_tags[i]<<" "<<std::endl;
	  // std::cout<<std::endl;
	  // }
	// auto& dm=W.spaces_ptr()->dofmap2();
	// auto& dm0=tuple_get<0>(dm);
	// auto& dm1=tuple_get<1>(dm);

	//  std::cout<<"dm0="<<std::endl;
	//  for(int el=0;el<mesh.n_elements();el++)
	//   {
	//     std::cout<<"elem="<<el<<std::endl;
	//     if(mesh.is_active(el))
	//     {
	//        std::cout<<dm0[el]<<std::endl; 
	//     }
	//   }
	//  std::cout<<"dm1="<<std::endl;
	//  for(int el=0;el<mesh.n_elements();el++)
	//   {
	//     std::cout<<"elem="<<el<<std::endl;
	//     if(mesh.is_active(el))
	//     {
	//        std::cout<<dm1[el]<<std::endl; 
	//     }
	//   }



	  // SparseMatrix<Real> M;
	  // M.init(6,7);

	  // M.plus_equal(4,0,0);
	  // M.plus_equal(3,0,0);
	  // M.plus_equal(4,0,2);
	  // M.plus_equal(3,2,0);
	  // M.plus_equal(3,4,0);
	  // M.plus_equal(3,0,4);
	  // std::vector<Real> vecme{1,2,3,4,5};
	  // M.print();
	  // std::cout<<M.multiply(0,vecme);

	 // std::cout<<std::endl;


	 //  auto res=M.multiply(vecme);
	 // for(int i=0;i<res.size();i++)
	 //    std::cout<<res[i]<<std::endl;





	  // decltype(L2Inner(u0,v0))::TestOrTrialLeftType m1(1);
	  // decltype(L2Inner(u0,v0))::TestOrTrialRightType m2(2);

	  // std::cout<<"elem_id="<<J.elem_id()<<std::endl;
	  // std::cout<<"------_______-----_______-----_______-----_______------"<<std::endl;
	  // std::cout<<"------_______-----BEGIN EVALUATION-----_______--------"<<std::endl;

	  // std::cout<<"------_______-----BILINEAR VOLUME-----_______--------"<<std::endl;
	  // eval_bilinear_form.apply(A,J);
	  // // std::cout<<"------_______-----BILINEAR SURFACE-----_______--------"<<std::endl;
	  // J.init_boundary(0);
	  // reference_maps.init_boundary(J);
	  // shapefunctions.init_boundary(J);
	  // eval_bilinear_form.apply_boundary(A,J);
	  // std::cout<<"------_______-----LINEAR VOLUME-----_______--------"<<std::endl;

	  // eval_linear_form.apply(b,J);
	  // J.init_boundary(0);
	  // reference_maps.init_boundary(J);
	  // shapefunctions.init_boundary(J);

	  // std::cout<<"------_______-----LINEAR SURFACE 0-----_______-----_______--------"<<std::endl;

	  // J.init_boundary(0);
	  // reference_maps.init_boundary(J);
	  // shapefunctions.init_boundary(J);
	  // eval_linear_form.apply_boundary(b,J);

	  // std::cout<<"------_______-----LINEAR SURFACE 1-----_______-----_______--------"<<std::endl;

	  // J.init_boundary(1);
	  // reference_maps.init_boundary(J);
	  // shapefunctions.init_boundary(J);
	  // eval_linear_form.apply_boundary(b,J);

	  // std::cout<<"------_______-----LINEAR SURFACE 2-----_______-----_______--------"<<std::endl;

	  // J.init_boundary(2);
	  // reference_maps.init_boundary(J);
	  // shapefunctions.init_boundary(J);
	  // eval_linear_form.apply_boundary(b,J);

	  // J.init_boundary(1);

	  // reference_maps.init_boundary(J);
	  // shapefunctions.init_boundary(J);
	  // eval_linear_form.apply_boundary(J);

	  // J.init_boundary(2);
	  // reference_maps.init_boundary(J);
	  // shapefunctions.init_boundary(J);
	  // eval_linear_form.apply_boundary(J);

	  // J.init_boundary(0);
	  // std::cout<<"side volume="<<std::endl;
	  // std::cout<<J.side_id()<<std::endl;
	  // std::cout<<J.side_volume()<<std::endl;
	  // J.init_boundary(1);
	  // std::cout<<"side volume="<<std::endl;
	  // std::cout<<J.side_id()<<std::endl;
	  // std::cout<<J.side_volume()<<std::endl;
	  // J.init_boundary(2);
	  // std::cout<<"side volume="<<std::endl;
	  // std::cout<<J.side_id()<<std::endl;
	  // std::cout<<J.side_volume()<<std::endl;
	  // std::cout<<shapefunctions.value_surface<1,0>();
	  // decltype(shapefunctions.get_surface<1,0>()) eek(6);

	  // std::cout<<"------_______-----_______-----_______-----_______------"<<std::endl;
	  // std::cout<<"------_______-----END EVALUATION-----_______--------"<<std::endl;



	  // std::cout<<"tags="<<std::endl;

	  // for(int el=0;el<mesh.n_elements();el++)
	  //   {
	  //     // auto side_tags=mesh.elem(el).side_tags;
	  //       for(int j=0;j<mesh.elem(el).side_tags.size();j++)
	  //         {
	  //           // mesh.elem(el).side_tags[j]=j;
	  //         }
	  //       // std::cout<<side_tags[j]<<std::endl;
	  //   }
	  // for(int el=0;el<mesh.n_elements();el++)
	  //   {
	  //     auto side_tags=mesh.elem(el).side_tags;
	  //       for(int j=0;j<mesh.elem(el).side_tags.size();j++)
	  //       std::cout<<side_tags[j]<<std::endl;
	  //   }
	  //   std::cout<<"mesh.side_nodes()="<<std::endl;
	    // std::cout<<mesh.side_nodes()<<std::endl;


	// decltype(L2Inner(v1,u1))::form mm(5);
	    // using t=std::tuple<std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>>;
	    // constexpr Integer HHH=-2;
	    // TupleChangeType<HHH,IdentityOperator,t> oo(4);
	    // SubTupleHelper<0,HHH+1,TupleTypeSize<t>::value-1,t>::type kk(5);
	    // Number<TupleTypeSize<t>::value> lk(5);
	    // using a1=TupleChangeTypeHelper<HHH,IdentityOperator,t>::typeLeft;
	    // using a2=TupleChangeTypeHelper<HHH,IdentityOperator,t>::typeRight;
	    // TupleCatType<a1,std::tuple<IdentityOperator>,a2> ok5(5);
	    // a1 oo66(1);
	    // a2 kk55(3);

	    // TupleChangeTypeHelper<HHH,IdentityOperator,t>::typeLeft k1(1);
	    // TupleChangeTypeHelper<HHH,IdentityOperator,t>::typeRight k2(1);

	    // TupleChangeType<0,IdentityOperator,std::tuple<std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>,std::tuple<>>> oo(4);

	    // TupleChangeType<1,char,std::tuple<TraceOperator,IdentityOperator,GradientOperator,CurlOperator,DivergenceOperator>> oo(4);

	// constexpr const auto face_dofs=decltype(W)::faces_dofs_array;
	// constexpr const auto Nface_dofs=decltype(W)::Nfaces_dofs_array;

	// std::cout<<"faces_dofs_array="<<tuple_get<0>(face_dofs)<<std::endl; 
	// std::cout<<"faces_dofs_array="<<tuple_get<1>(face_dofs)<<std::endl; 
	// std::cout<<"faces_dofs_array="<<tuple_get<2>(face_dofs)<<std::endl; 
	// std::cout<<"faces_dofs_array="<<tuple_get<3>(face_dofs)<<std::endl; 
	// std::cout<<"Nfaces_dofs_array="<<Nface_dofs<<std::endl; 

	// constexpr const auto _0dofs=tuple_get<0>(face_dofs)[0].size();
	// static_assert(_0dofs==4 &&"ok");
	// decltype(Wtrial)::DofMapType2 o3o(5,4,5,6,6);
	// decltype(W)::DofMapType2 oo(5,4,5,6,6);
	// auto dof0=tuple_get<0>(dof);
	// auto dof1=tuple_get<1>(dof);
	// // auto dof2=tuple_get<2>(dof);
	// std::cout<<"dofmap0="<<std::endl;
	// for(Integer i=0;i<dof0.size();i++)
	//   {for(Integer j=0;j<dof0[i].size();j++)
	//      std::cout<<dof0[i][j]<<std::endl;
	//      std::cout<<std::endl;

	//    }
	// std::cout<<"dofmap1="<<std::endl;

	// for(Integer i=0;i<dof1.size();i++)
	//   {for(Integer j=0;j<dof1[i].size();j++)
	//      std::cout<<dof1[i][j]<<std::endl;
	//      std::cout<<std::endl;
	//    }



	// auto _dof1=FEspace1.dofmap();
	// auto _dof3=FEspace3.dofmap();



	// for(Integer i=0;i<_dof1.size();i++)
	//   {for(Integer j=0;j<_dof1[i].size();j++)
	//      std::cout<<_dof1[i][j]<<std::endl;
	//      std::cout<<std::endl;
	//    }

	// for(Integer i=0;i<_dof3.size();i++)
	//   {for(Integer j=0;j<_dof3[i].size();j++)
	//      std::cout<<_dof3[i][j]<<std::endl;
	//      std::cout<<std::endl;
	//    }


	// auto dof1_=FEspace1.dofmap2();
	// auto dof3_=FEspace3.dofmap2();

	// auto dof1=tuple_get<0>(dof1_);
	// auto dof2=tuple_get<0>(dof3_);
	// auto dof3=tuple_get<1>(dof3_);



	// std::cout<<"dofmap____1="<<std::endl;

	// for(Integer i=0;i<dof1.size();i++)
	//   {for(Integer j=0;j<dof1[i].size();j++)
	//      std::cout<<dof1[i][j]<<std::endl;
	//      std::cout<<std::endl;
	//    }

	// std::cout<<"dofmap____2="<<std::endl;

	// for(Integer i=0;i<dof2.size();i++)
	//   {for(Integer j=0;j<dof2[i].size();j++)
	//      std::cout<<dof2[i][j]<<std::endl;
	//      std::cout<<std::endl;
	//    }

	// std::cout<<"dofmap____3="<<std::endl;

	// for(Integer i=0;i<dof3.size();i++)
	//   {for(Integer j=0;j<dof3[i].size();j++)
	//      std::cout<<dof3[i][j]<<std::endl;
	//      std::cout<<std::endl;
	//    }


	//    auto DF2=W.spaces_ptr()->dofmap2();

	//    auto d1=tuple_get<0>(DF2);
	//    auto d2=tuple_get<1>(DF2);
	//    auto d3=tuple_get<2>(DF2);



	// std::cout<<"------------------------d1="<<std::endl;
	//    for(Integer i=0;i<d1.size();i++)
	//   {for(Integer j=0;j<d1[i].size();j++)
	//      std::cout<<d1[i][j]<<std::endl;
	//      std::cout<<std::endl;
	//    }

	// std::cout<<"------------------------d2="<<std::endl;
	//    for(Integer i=0;i<d2.size();i++)
	//   {for(Integer j=0;j<d2[i].size();j++)
	//      std::cout<<d2[i][j]<<std::endl;
	//      std::cout<<std::endl;
	//    }

	// std::cout<<"------------------------d3="<<std::endl;
	//    for(Integer i=0;i<d3.size();i++)
	//   {for(Integer j=0;j<d3[i].size();j++)
	//      std::cout<<d3[i][j]<<std::endl;
	//      std::cout<<std::endl;
	//    }
	// for(Integer i=0;i<dof2.size();i++)
	//   for(Integer j=0;j<dof2[i].size();j++)
	//      std::cout<<dof2[i][j]<<std::endl;
	// decltype(eval_linear_form)::EvaluationOfL2InnersVolume::L2Products a4k(6);

	// decltype(bilinear_form)::FunctionSpace a444k(6);
	// auto l2prdo=L2Inner(v0,u0);
	// std::cout<<"kkkk3333="<<std::endl;
	// auto spcptr=find_spaces_ptr<decltype(bilinear_form)::FunctionSpace>(bilinearform);//C-v0+v0*C*C*(-Transpose(f1)));
	// if(spcptr==nullptr)
	//  std::cout<<"nullptr="<<std::endl;
	// else
	// {
	//   auto ecm=spcptr->dofmap2();

	//    std::cout<<"kkkk="<<std::endl;
	//    std::cout<<tuple_get<1>(ecm)[0][4]<<std::endl;
	// }

	// Number<TupleTypeSize<decltype(bilinear_form)::TupleOfPairsNumbers>::value> mb(5);
	// decltype(bilinear_form)::TupleOfPairsNumbers m3b(5);

	// decltype(eval_bilinear_form)::EvaluationOfL2InnersVolume::EvalOfL2InnersType s3k(6);
	// decltype(eval_linear_form)::EvaluationOfL2InnersSurface::EvalOfL2InnersType s32k(6);
	// Number<TupleTypeSize<decltype(eval_bilinear_form)::EvaluationOfL2InnersSurface::L2Products>::value-1> lkkljkkj(6);        
	// decltype(eval_linear_form)::EvaluationOfL2InnersSurface::EvalOfL2InnersType a44k(6);
	// decltype(eval_bilinear_form)::ShapeFunctions escho(6);
	// decltype(shapefunctions)::UniqueElementFunctionSpacesTupleType kk(5);
	// decltype(shapefunctions)::TupleOfTupleShapeFunctionSurface k3k4(5,4,54,6,7,8);
	// decltype(shapefunctions)::SpacesToUniqueFEFamilies k3444k4(5,4,54,6,7,8);
	// decltype(shapefunctions)::UniqueMapping k3eeeek4(5,4,54,6,7,8);
	// decltype(shapefunctions)::TupleOfTupleCompositeShapeFunctionSurface kk2(5);
	// decltype(shapefunctions)::SurfaceForm k3k(5);
	// Number<decltype(shapefunctions)::GeneralForm::FunctionSpace::Nuniquesubspaces> kjn(5);
	// decltype(surface_integral(NewTrace(u1),NewTrace(v1)))::QRule oki(6);
	// Number<QuadratureOrder<decltype(NewTrace(u1))>::value> ok4i5(6);
	// decltype(NewTrace(u1)) mio(5);

	// decltype(shapefunctions.surface_tuple()) jk(5);
	// decltype(shapefunctions.composite_shapes_surface()) j4k(5);
	// auto em=surface_integral(NewTrace(u2),NewTrace(v2));
	// using GeneralForm=decltype(eval_bilinear_form);
	// TupleOfL2Products2<1,typename GeneralForm::TupleOfPairsNumbers, typename GeneralForm::type::Form >::type ok(5);
	// auto escom=
	// L2Inner((u0),(v0))-
	// surface_integral(Trace(u2),Trace(v2))
	// -
	// L2Inner((u0),(v0))+
	// surface_integral(Trace(u0),Trace(v0))-
	// surface_integral(Trace(u2),Trace(v2))-
	// L2Inner((u2),(v2))
	// -
	// L2Inner((u1),(v1))
	// -
	// surface_integral(Trace(u2),Trace(v2))-
	// surface_integral(Trace(u2),Trace(v2))-
	// L2Inner((u0),(v0))-
	// surface_integral(Trace(u2),Trace(v2))-
	// surface_integral(Trace(u2),Trace(v2))-
	// L2Inner((u0),(v0))+
	// L2Inner((u0),(v0))-
	// surface_integral(Trace(u0),Trace(v0))-
	// L2Inner((u0),(v0));

	// auto ooo=ExtractForm<1>(escom);//VolumeForm<-1>(VolumeForm<-1>(bilinearform));
	// decltype(ooo) klib(5);
	  // OperatorType<ExtractFormType<1,decltype(bilinearform)>> ee(65);
	// decltype(escom) ee(6);
	 // auto bilinearform2=surface_integral(Trace(u2),Trace(v2));
	  //  // std::vector<Vector<Real,2>> points;
	  //  const auto & elem=mesh.elem(0);
	  //  Simplex<2,1> simplex_side;

	  // mesh.points(0,points);
	  // elem.side(0,simplex_side);
	  // std::cout<<"side volume="<<unsigned_volume(simplex_side,points)<<std::endl;
	  // elem.side(1,simplex_side);
	  // std::cout<<"side volume="<<unsigned_volume(simplex_side,points)<<std::endl;
	  // elem.side(2,simplex_side);
	  // std::cout<<"side volume="<<unsigned_volume(simplex_side,points)<<std::endl;
	  // for (int ii=0;ii<points.size();ii++)
	  // std::cout<<"points="<<points[ii]<<std::endl;

	// ElemToSubElemHelper<Simplex<3,3>,2>::type ok(1);
	// decltype(shapefunctions)::TupleOfTupleShapeFunction oi(1);

	// Number<FunctionSpaceDofsPerElem<ElemFunctionSpace<Elem,RT1<1>>>::value> kjh(5);

	// Number<FunctionSpaceDofsPerSubEntityElem<ElemFunctionSpace<Elem,Lagrange3<1>>,1>::value> kjh(5);
	// decltype(shapefunctions)::SpacesToUniqueFEFamilies ook(3);
	// decltype(reference_maps)::UniqueMappingVolumetric o4i(1);
	// decltype(reference_maps)::UniqueMappingSurface okkl(4);
	// std::cout<<"----"<<tuple_get<0,0>(reference_maps())()<<std::endl;
	// std::cout<<"----"<<tuple_get<1,0>(reference_maps())()<<std::endl;
	// J.init_boundary(0);
	// reference_maps.init_boundary(J);
	// std::cout<<"----"<<tuple_get<1,0>(reference_maps.surface_map())()<<std::endl;
	// J.init_boundary(1);
	// reference_maps.init_boundary(J);
	// std::cout<<"----"<<tuple_get<1,0>(reference_maps.surface_map())()<<std::endl;
	// J.init_boundary(2);
	// reference_maps.init_boundary(J);
	// std::cout<<"----"<<tuple_get<1,0>(reference_maps.surface_map())()<<std::endl;
	// decltype(shapefunctions)::MapTupleNumbers ookll9(5);
	// decltype(shapefunctions)::TupleOperatorsAndQuadrature ee(1);
	// decltype(shapefunctions)::TupleOfTupleCompositeShapeFunction klklk(6);
	// decltype(shapefunctions)::TupleOfTupleCompositeShapeFunction esde(5);


	// decltype(shapefunctions)::MapCollection::UniqueMapping klk3lk(6);
	// decltype(reference_maps)::UniqueMappingSurface ok(6);
	// decltype(shape_coefficients)::UniqueElementFunctionSpacesTupleType lh(3);
	// decltype(shape_coefficients)::SpacesToUniqueFEFamily  o43i(4);
	// decltype(shape_functions)::TupleOfTupleCompositeShapeFunction ee(1);
	// decltype(linear_form()) ee3(1);

	 // auto bilinearform2= L2Inner(u0,v0)-L2Inner(u0,v0)-L2Inner(u0,v0);//+ L2Inner(Grad(u0),Grad(v0))+L2Inner(u0,v2);//+L2Inner(Grad(u1),Grad(v1));//+L2Inner(f3*u2,v2);


	 // decltype(bilinearform) ooo(5);
	 // auto newform= L2Inner(Epsilon(u0),Epsilon(v0));//+ L2Inner(Grad(u0),Grad(v0))+L2Inner(u0,v2);//+L2Inner(Grad(u1),Grad(v1));//+L2Inner(f3*u2,v2);
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

	// OperatorType<decltype(Transpose(Epsilon(v0))), GaussPoints<Simplex<2, 2>, 1> > aa;
	// decltype(shape_functions)::TupleCompositeOperatorsAndQuadrature ok03(0);
	// OperatorAndQuadratureTupleType<decltype(linearform)>::type ee(5);
	// OperatorAndQuadratureTupleType<decltype(linearform)>::composite_type e4e(5);

	// decltype(L2Inner((u2),(v2+v1))) ee1(5);
	// decltype(L2Inner((u2),(v2-v1))) ee2(5);
	// decltype(L2Inner((u2+u1),(v2))) ee3(5);
	// decltype(L2Inner((u2-u1),(v2))) ee4(5);
	// decltype(L2Inner((u2+u1),(v2+v1))) ee5(5);
	// decltype(L2Inner((u2-u1),(v2-v1))) ee6(5);
	// decltype(L2Inner((u2+u1),(v2-v1))) ee7(5);
	// decltype(L2Inner((u2-u1),(v2+v1))) ee8(5);

	// decltype(L2Inner((u2-u1-u0),(v2+v1+v0))) ee9(5);
	// decltype(L2Inner((u2-u1-u0),(v2+v1-v0))) ee10(5);
	// decltype(L2Inner((u2-u1-u0),(v2-v1+v0))) ee11(5);
	// decltype(L2Inner((u2-u1-u0),(v2-v1-v0))) ee12(5);


	// decltype(L2Inner((u2-u1-u0),(v2+v1+v0))) ee13(5);
	// decltype(L2Inner((u2-u1+u0),(v2+v1-v0))) ee14(5);
	// decltype(L2Inner((u2+u1-u0),(v2-v1+v0))) ee15(5);
	// decltype(L2Inner((u2+u1+u0),(v2-v1-v0))) ee16(5);

	// decltype(L2Inner((C*u2+C*u1+C*u0),(C*v2-C*v1-C*v0))) ee17(5);
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
	//   // OperatorAndQuadratureTupleType<decltype(NewOp2(v0))>::type ee(65,4,4,4,4);
	// // decltype(reference_maps)::TupleOperatorsAndQuadrature eee(5);
	// // decltype(reference_maps)::TupleCompositeOperatorsAndQuadrature ee4e(5);
	// decltype(reference_maps)::TupleOfTupleNoQuadratureSurface ee444e(5);
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


	//  auto mm=Evaluation<Expression<decltype(alpha*u0)>, GaussPoints<Simplex<2,2>,3>>((alpha*u0));


	 // decltype(shape_functions)::TupleOfTupleCompositeShapeFunctionTensor e(6);
	 // decltype(shape_functions)::TupleOfTupleCompositeShapeFunction e4(6);

	// FQPValues<Matrix<double, 1L, 2L>, 1L, 3L> hh;
	//  auto ecc=form_of_composite_operator(NewOp1(v0));
	//  auto ecc2=Evaluation<Expression<decltype(ecc)>,GaussPoints<Simplex<2,2>,1>>(ecc);
	//  ecc2.apply(hh,J,shape_functions());
	 
	//  // OperatorType<decltype(NewOp1(v0)),GaussPoints<Simplex<2,2>,1>> e(6);

	//   TupleOfCombinationFunctions<decltype(u0)::Nmax,MultipleAddition<decltype(bilinearform),decltype(linearform)>>::type onlycomposite=build_tuple_of_combination_functions<decltype(u0)::Nmax>(bilinearform,linearform);








	// FormOfCompositeOperatorType<decltype(NewOp1(v0))>::type ok(1);
	// decltype(shape_functions)::TupleOfTupleCompositeShapeFunction e4ee(6);
	// decltype(shape_functions)::TupleOfTupleCompositeShapeFunctionTensor eee(6);



	// TupleOfCombinationFunctions<decltype(u0)::Nmax,decltype(bilinearform),decltype(linearform)>::type eee(6);

	 // decltype(onlycomposite) onle333ycomposite;
	 // decltype(shape_functions()) e(6);
	 // static_assert(IsSame<decltype(onlycomposite),typename TupleOfCombinationFunctions<decltype(u0)::Nmax,MultipleAddition<decltype(bilinearform),decltype(linearform)>>::type>::value && "they are not same");

	// Number<decltype(bilinear_form)::FunctionSpace::Nuniquesubspaces > ee(5);
	// Number<decltype(linear_form)::FunctionSpace::Nuniquesubspaces > e4(5);

	// Number<decltype(u0)::Nmax> eee(6);



	// decltype(build_tuple_of_combination_functions<decltype(u0)::Nmax>(decltype(bilinearform)(),decltype(linearform)())) ekl(6);

	// decltype(shape_functions()) eee(6);

	// decltype(u0)::UniqueElementFunctionSpacesTupleType oo(5);
	// Number<decltype(u0)::value> o1(3);
	// Number<decltype(u1)::value>o41(3);
	// Number<decltype(u2)::value>o14(3);



	// std::cout<<hh<<std::endl;
	// OperatorType<Division<FQPValues<Matrix<double,1,2,-1>,1,3>,QPValues<Matrix<double,1,1,-1>,1>>> eee(5);
	// OperatorType<Multiplication<Matrix<double, 1,
	//       2, -1>, Matrix<double, 1, 1, -1> > >  eee(6,5,5,5);
	// decltype(shape_functions()) oo(5);
	// decltype(ecc2) l(5);

	//   std::cout<<"apply2"<<std::endl;
	  // auto Ku=4*Grad(u0);
	  // auto Kv=4*Grad(v0);
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
	    // using Right=decltype(v2);
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


	  // decltype(L2Inner(f3,v2))::TestTrialNumbers eeee(6);
	  // auto seee=L2Inner(f3,v2);


	// Evaluation<Expression<decltype(f3*u2)>,GaussPoints<Simplex<2,2>,5>>::subtype okk(6); 


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

	   // nonzero=StaticScalarProduct<Nonzero>(v0,v1);
	    // StaticMatrixProduct<findnonzeronumbers,nonzeroproduct>(mat1,mat2,mat3,vec1,vec2);
	    //     // nonzero=0;


	              // mat3= mat1 * mat2;

	     // nonzero=v0[3]*v1[2]+v0[4]*v1[5]+v0[5]*v1[8];
	    // nonzero=0;
	    // for(Integer ii=0;ii<3;ii++)
	    //   nonzero+=mat1(1,ii)*mat2(ii,2);
	}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	// std::cout<<"mat==="<<mat3<<std::endl;

	// std::cout<<" elapsed_secs------------>"<<elapsed_secs<<std::endl;

	// std::cout<<" nonzero------------>"<<nonzero<<std::endl;

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

	// 	constexpr Integer ManifoldDim=2;
	// 	constexpr Integer Dim=2;
	// 	using MeshT=Mesh<Dim, ManifoldDim>;
	// 	MeshT mesh;
	// 	using Elem = typename MeshT::Elem; 
	// 	read_mesh("../data/beam-tri.MFEM", mesh);






	// 	NodeToElem<MeshT> node2elem3(mesh);
	// 	auto node2elem=node2elem3.val();

	// 	for(Integer nn=0; nn<node2elem.size();nn++)
	// 		{std::cout<<"node=="<<nn<<"     ";
	// 	for(Integer mm=0; mm<node2elem[nn].size();mm++)
	// 		std::cout<<node2elem[nn][mm]<<" ";
	// 	std::cout<<std::endl;
	// }
	// ElemEntity<Elem,0> nodes(mesh,node2elem);

	// std::cout<<"node 2 elem size=="<<nodes.entity_2_elem().size()<<std::endl;
	// std::cout<<"node 2 elem=="<<std::endl;
	// for(Integer nn=0; nn<nodes.entity_2_elem().size();nn++)
	// 	std::cout<<nodes.entity_2_elem(nn)[0]<<" ";
	// std::cout<<std::endl;
	// for(Integer nn=0; nn<nodes.elem_2_entity().size();nn++)
	// {
	// 	std::cout<<"elem="<<nn<< " made of nodes "<< std::endl;
	// 	for(Integer mm=0; mm<nodes.elem_2_entity(nn).size();mm++)
	// 		std::cout<<nodes.elem_2_entity(nn)[mm]<<" ";
	// 	std::cout<< std::endl;
	// } 


	// using FSspace= FunctionSpace< MeshT, Lagrange2<2>,RT0<1>>;
	// FSspace FEspace(mesh);


	// const auto& P2_ens0=ElemEntity<Elem,ElementFunctionSpace<Elem,LagrangeFE,2>::entity[0]>(mesh,node2elem);
	// const auto& P2_ens1=ElemEntity<Elem,ElementFunctionSpace<Elem,LagrangeFE,2>::entity[1]>(mesh,node2elem);

	// auto ens2elem20=P2_ens0.entity_2_elem();
	// auto ens2elem21=P2_ens1.entity_2_elem();

	// auto elem2ens20=P2_ens0.elem_2_entity();
	// auto elem2ens21=P2_ens1.elem_2_entity();

	// std::cout<<"ens2elem 2 0="<< std::endl;
	// for(Integer nn=0;nn<ens2elem20.size();nn++)
	// {
	// 	for(Integer mm=0;mm<ens2elem20[nn].size();mm++)
	// 		std::cout<<ens2elem20[nn][mm]<<" ";
	// 	std::cout<<std::endl;
	// } 
	// std::cout<<"ens2elem 2 1="<< std::endl;
	// for(Integer nn=0;nn<ens2elem21.size();nn++)
	// {
	// 	for(Integer mm=0;mm<ens2elem21[nn].size();mm++)
	// 		std::cout<<ens2elem21[nn][mm]<<" ";
	// 	std::cout<<std::endl;
	// } 
	// std::cout<<"elem2ens20 2 0="<< std::endl;
	// for(Integer nn=0;nn<elem2ens20.size();nn++)
	// {
	// 	for(Integer mm=0;mm<elem2ens20[nn].size();mm++)
	// 		std::cout<<elem2ens20[nn][mm]<<" ";
	// 	std::cout<<std::endl;
	// } 
	// std::cout<<"elem2ens21 2 1="<< std::endl;
	// for(Integer nn=0;nn<elem2ens21.size();nn++)
	// {
	// 	for(Integer mm=0;mm<elem2ens21[nn].size();mm++)
	// 		std::cout<<elem2ens21[nn][mm]<<" ";
	// 	std::cout<<std::endl;
	// } 

	// std::cout<<"n_dofs="<<FEspace.n_dofs()<< std::endl;
	// std::cout<<std::endl;
	// for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
	// {
	// 	auto &elem_id=elem_iter;
	// 	std::cout<<std::endl<<" elem =="<<elem_iter<<std::endl;
	// //  for(Integer nn=0;nn<FEspace.dofmap(elem_id).size();nn++)
	// //  {
	// //   std::cout<<FEspace.dofmap(elem_id)[nn]<< "  ";
	// // }
	// } 


	//    // FEspace.set_new_start(4);
	//    //  std::cout<<"dofmap_new_start n_dofs="<<FEspace.n_dofs()<< std::endl;
	//    //  std::cout<<std::endl;
	//    //  for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
	//    //  {
	//    //   auto &elem_id=elem_iter;
	//    //   std::cout<<std::endl<<" elem =="<<elem_iter<<std::endl;
	//    //   for(Integer nn=0;nn<FEspace.dofmap_new_start(elem_id).size();nn++)
	//    //   {
	//    //      std::cout<<FEspace.dofmap_new_start(elem_id)[nn]<< "  ";
	//    //   }
	//    //  } 


	// std::cout<<std::endl;
	// for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
	// 	std::cout<<"components of space["<<ss<<"]=="<<FEspace.components(ss)<<std::endl;
	// std::cout<<std::endl;

	// for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
	// {
	// 	std::cout<<"dofs of space ="<<ss<<std::endl;
	// 	for(Integer cc=0;cc<FEspace.components(ss);cc++)
	// 		{std::cout<<"component ="<<cc<<"   "<<std::endl;
	// 	auto& vec=FEspace.space_dofs(ss,cc);
	// 	for(Integer mm=0;mm<FEspace.n_dofs(ss,cc);mm++)
	// 		std::cout<<vec[mm]<<" ";
	// }
	// std::cout<<std::endl;


	// }


	// for(Integer mm=0;mm<FEspace.offset().size();mm++)
	// {
	// 	std::cout<<"OFFSET space ="<<mm<<std::endl;
	// 	for(Integer nn=0;nn<FEspace.offset()[mm].size();nn++)
	// 	{
	// 		std::cout<< FEspace.offset()[mm][nn]<<" ";
	// 	}
	// 	std::cout<<std::endl;
	// }
	// std::cout<<std::endl;


	// auto spaceinf=FEspace.space_info();
	// for(Integer mm=0;mm<spaceinf.size();mm++)
	// {
	// 	std::cout<<"Space=="<<mm<<std::endl;
	// 	for(Integer nn=0;nn<spaceinf[mm].size();nn++)
	// 		std::cout<<spaceinf[mm][nn]<<" ";
	// 	std::cout<<std::endl;
	// }

	// std::cout<<"Whole dofmap=="<<std::endl;
	// for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
	// {
	//    //  const auto& dm=FEspace.dofmap(0,elem_iter);
	//    //  for(auto& i:dm)
	//    //   std::cout<<i<<" ";
	//    // std::cout<<std::endl;
	// }

	// std::cout<<"First Space, first component, dofmap=="<<std::endl;

	// for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
	// {
	//    //  const auto& dm=FEspace.dofmap(0,0,elem_iter);
	//    //  for(auto& i:dm)
	//    //   std::cout<<i<<" ";
	//    // std::cout<<std::endl;
	// }
	// std::cout<<"First Space, second component, dofmap=="<<std::endl;
	// for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
	// {
	//    //  const auto& dm=FEspace.dofmap(0,1,elem_iter);
	//    //  for(auto& i:dm)
	//    //   std::cout<<i<<" ";
	//    // std::cout<<std::endl;
	// }

	// std::cout<<"Second Space, first component, dofmap=="<<std::endl;
	// for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
	// {
	//    //  const auto& dm=FEspace.dofmap(1,0,elem_iter);
	//    //  for(auto& i:dm)
	//    //   std::cout<<i<<" ";
	//    // std::cout<<std::endl;
	// }
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


		NodeToElem<MeshT> node2elem3(mesh);
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


	// 	FunctionSpace< MeshT, Lagrange1<2>, Lagrange2<1> > FEspace(mesh);

	//     // auto eeh=FunctionSpaceSystem(FEspace);


	// 	std::cout<<"n_dofs="<<FEspace.n_dofs()<< std::endl;
	// 	std::cout<<std::endl;
	// 	for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
	// 	{
	// 		auto &elem_id=elem_iter;
	//   //  std::cout<<"elem_id="<<elem_id<<", number of dofs=s"<<FEspace.dofmap(elem_id).size()<<std::endl;
	//   //  for(Integer nn=0;nn<FEspace.dofmap(elem_id).size();nn++)
	//   //  {
	//   //   std::cout<<FEspace.dofmap(elem_id)[nn]<<" ";
	//   // }
	// 		std::cout<<std::endl;
	// 	} 
	// 	std::cout<<std::endl;
	// 	for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
	// 		std::cout<<"components of space["<<ss<<"]=="<<FEspace.components(ss)<<std::endl;
	// 	std::cout<<std::endl;

	// 	for(Integer ss=0;ss<FEspace.n_subspaces();ss++)
	// 	{
	// 		std::cout<<"dofs of space ="<<ss<<std::endl;
	// 		for(Integer cc=0;cc<FEspace.components(ss);cc++)
	// 			{std::cout<<"component ="<<cc<<"   "<<std::endl;
	// 		auto& vec=FEspace.space_dofs(ss,cc);
	// 		for(Integer mm=0;mm<FEspace.n_dofs(ss,cc);mm++)
	// 			std::cout<<vec[mm]<<" ";
	// 		std::cout<<std::endl;
	// 	}
	// 	std::cout<<std::endl;

	// }


	// for(Integer mm=0;mm<FEspace.offset().size();mm++)
	// {
	// 	std::cout<<"offset space="<<mm<<std::endl;
	// 	for(Integer nn=0;nn<FEspace.offset()[mm].size();nn++)
	// 	{
	// 		std::cout<< FEspace.offset()[mm][nn]<<" ";
	// 	}
	// }
	// std::cout<<std::endl;


	// for(Integer elem_iter=0;elem_iter<mesh.n_elements();elem_iter++)
	// 	{std::cout<<std::endl;
	//     // const auto size=FEspace.dofmap(1,elem_iter).size();
	//     // std::cout<<"elem_iter="<<elem_iter<<std::endl;
	//     // for(Integer nn=0;nn<size;nn++)
	//     //   std::cout<<FEspace.dofmap(1,elem_iter)[nn]<<" ";
	//     // std::cout<<std::endl;
	// 	}

	// 	std::cout<<std::endl;
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


		NodeToElem<MeshT> node2elem3(mesh);
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
		ElemConnectivity<MeshT,entitydim_from,subentitydim_from,entitydim_to> conn_e2t(mesh,node2elem);



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
